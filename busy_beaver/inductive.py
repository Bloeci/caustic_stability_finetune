"""Inductive (bouncer) prover: sound non-halting proofs by symbolic induction.

Mechanism
---------
1. Run the accelerated simulator and record, per configuration *shape*
   (state + run symbols), the vector of run lengths at each visit.
2. If a shape recurs with a constant non-negative increment c between visits
   (v, v+c, v+2c, ...), conjecture the family C(n) with run lengths v + n*c.
3. Prove C(n) ->* C(n+1) for all n >= n0 by simulating symbolically: run
   lengths are affine polynomials a*n + b. Every operation that could depend
   on the concrete value of n (removing one cell from a run) must behave
   identically for all n >= n0; when needed, n0 is raised, and the proof is
   only accepted if the real trace already witnessed C(0..n0) (base cases).
   No halt transition may occur along the symbolic path.
4. C(n0) is reachable from the blank tape (it was observed), and
   C(n) ->* C(n+1) halt-free forever, so the machine never halts.

This is a working miniature of the "inductive rules / formula tape"
technique behind modern Busy Beaver results: simulation discovers the
pattern, induction over a symbolic tape turns it into a proof.
"""
from __future__ import annotations

from .machine import TM, HALT, R
from .simulate import simulate_accel, NONHALT_ESCAPE
from .rule_miner import ShapeRecorder

# A symbolic run length is (a, b) meaning a*n + b with integers a >= 0.


class _Fail(Exception):
    pass


class _Sym:
    """Symbolic configuration with an adjustable lower bound n >= n0."""

    def __init__(self, q, left, right, n0, n0_cap):
        self.q = q
        self.left = left    # list of [sym, (a, b)]
        self.right = right
        self.n0 = n0
        self.n0_cap = n0_cap

    def _require_pos_after_dec(self, a, b):
        """Ensure a*n + b - 1 >= 1 for all n >= n0 (raising n0 if allowed)."""
        if a == 0:
            if b - 1 < 1:
                raise _Fail()
            return
        # smallest n0 with a*n0 + b - 1 >= 1:
        need = -(-(2 - b) // a)  # ceil((2-b)/a)
        if need > self.n0:
            if need > self.n0_cap:
                raise _Fail()
            self.n0 = need

    def pop_one(self, stack):
        """Remove one cell from the top run; return its symbol (0 if sea)."""
        if not stack:
            return 0
        sym, (a, b) = stack[-1]
        if a == 0 and b == 1:
            stack.pop()
            return sym
        self._require_pos_after_dec(a, b)
        stack[-1][1] = (a, b - 1)
        return sym

    @staticmethod
    def push(stack, sym, cnt):
        if stack:
            top = stack[-1]
            if top[0] == sym:
                a, b = top[1]
                top[1] = (a + cnt[0], b + cnt[1])
                return
        elif sym == 0:
            return
        stack.append([sym, cnt])

    def key(self):
        return (
            self.q,
            tuple((s, c) for s, c in self.left),
            tuple((s, c) for s, c in self.right),
        )


def _sym_step(cfg: _Sym, tm: TM):
    """One symbolic macro step. Returns None, or 'escape' for a direct proof."""
    s = cfg.right[-1][0] if cfg.right else 0
    tr = tm.get(cfg.q, s)
    if tr is None or tr[2] == HALT:
        raise _Fail()  # a halt is (symbolically) reachable: no proof
    w, d, q2 = tr
    if d == R:
        if not cfg.right:
            if q2 == cfg.q:
                return "escape"
            cfg.push(cfg.left, w, (0, 1))
        elif q2 == cfg.q:
            sym, cnt = cfg.right.pop()
            cfg.push(cfg.left, w, cnt)
        else:
            cfg.pop_one(cfg.right)
            cfg.push(cfg.left, w, (0, 1))
    else:
        cfg.pop_one(cfg.right)
        extra = (0, 0)
        if q2 == cfg.q:
            if cfg.left and cfg.left[-1][0] == s:
                extra = cfg.left.pop()[1]
            elif not cfg.left and s == 0:
                return "escape"
        cfg.push(cfg.right, w, (extra[0], extra[1] + 1))
        y = cfg.pop_one(cfg.left)
        cfg.push(cfg.right, y, (0, 1))
    cfg.q = q2
    return None


def _try_shape(tm, shape, occs, max_sym_steps):
    """Try to prove non-halting from one shape's visit history."""
    # longest suffix with constant difference
    diffs = [
        tuple(b - a for a, b in zip(occs[i], occs[i + 1]))
        for i in range(len(occs) - 1)
    ]
    k = len(diffs) - 1
    while k > 0 and diffs[k - 1] == diffs[-1]:
        k -= 1
    c = diffs[-1]
    m = len(diffs) - k  # number of constant-diff transitions observed
    if m < 3:
        return False
    if any(x < 0 for x in c) or all(x == 0 for x in c):
        return False  # shrinking or pure cycle (cycler decider's job)
    base = occs[k]  # C(0)
    n0_cap = m - 1  # base cases C(0..n0) must all be observed
    q, lsyms, rsyms = shape
    counts = list(base)
    left = [[s, (c[i], counts[i])] for i, s in enumerate(lsyms)]
    right = [
        [s, (c[len(lsyms) + i], counts[len(lsyms) + i])]
        for i, s in enumerate(rsyms)
    ]
    cfg = _Sym(q, left, right, 0, n0_cap)
    target = (
        q,
        tuple((s, (c[i], counts[i] + c[i])) for i, s in enumerate(lsyms)),
        tuple(
            (s, (c[len(lsyms) + i], counts[len(lsyms) + i] + c[len(lsyms) + i]))
            for i, s in enumerate(rsyms)
        ),
    )
    try:
        for _ in range(max_sym_steps):
            res = _sym_step(cfg, tm)
            if res == "escape":
                return True
            if cfg.key() == target:
                return True
            if len(cfg.left) + len(cfg.right) > 64:
                return False
    except _Fail:
        return False
    return False


def _decide_inductive_one(
    machine,
    trace_macro,
    trace_base,
    max_shapes,
    min_occ,
    max_sym_steps,
):
    rec = ShapeRecorder(max_runs=10, max_shapes=2000, max_occ=200)
    status, _, _, _ = simulate_accel(
        machine, max_macro=trace_macro, max_base=trace_base, recorder=rec
    )
    if status == NONHALT_ESCAPE:
        # Sound also for a BlockMachine: an endless chain of blank blocks in a
        # fixed macro state is an endless run of the underlying machine.
        return True
    candidates = sorted(
        (
            (shape, occs)
            for shape, occs in rec.data.items()
            if len(occs) >= min_occ
        ),
        key=lambda t: -len(t[1]),
    )
    for shape, occs in candidates[:max_shapes]:
        # stride > 1 handles e.g. two-sided bouncers, whose shape visits
        # alternate sides so only every 2nd visit has a constant increment
        for stride in (1, 2, 3, 4):
            for offset in range(stride):
                sub = occs[offset::stride]
                if len(sub) >= min_occ and _try_shape(
                    machine, shape, sub, max_sym_steps
                ):
                    return True
    return False


def decide_inductive(
    tm: TM,
    trace_macro=4000,
    trace_base=10**6,
    max_shapes=6,
    min_occ=5,
    max_sym_steps=20000,
    block_sizes=(1, 2, 3, 4),
):
    """True iff a sound inductive non-halting proof is found.

    Tries the raw machine and k-block macro machines: block grouping turns
    alternating patterns such as (10)^n into single growing runs that the
    affine induction can capture.
    """
    from .blocks import BlockMachine

    for k in block_sizes:
        machine = tm if k == 1 else BlockMachine(tm, k)
        if _decide_inductive_one(
            machine, trace_macro, trace_base, max_shapes, min_occ, max_sym_steps
        ):
            return True
    return False
