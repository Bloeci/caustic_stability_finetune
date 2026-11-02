"""Counter-abstraction closure decider (a RepWL / CTL-style regular invariant).

Abstract domain: the tape is a sequence of runs (sym, cnt) where cnt is kept
exact while < T and widened to ('T', r) - meaning "any length >= T congruent
to r mod m" - once it reaches T. One abstract configuration therefore denotes
a *regular set* of concrete configurations. The abstract step relation
over-approximates the concrete one (popping a cell from a widened run
branches into "still >= T" and, when the congruence allows length exactly T,
"now exactly T-1"), so the set of abstract configurations reachable from the
blank tape, if finite and halt-free, is a regular invariant: closed under
the step relation, containing the initial configuration and excluding every
halting configuration - a sound non-halting proof.

The modulus m preserves congruence information that plain thresholding
loses (many bouncers/counters only ever produce run lengths of one parity,
and the spurious other parity is what reaches the halt transition). Block
machines add the same power along the word dimension.

This is the mechanism behind bbchallenge's RepWL and finite-automata
deciders, in its simplest form.
"""
from __future__ import annotations

from .machine import TM, HALT, R


def _push(runs, sym, T, m):
    """Push one cell of sym onto a run tuple (top = last element)."""
    if runs:
        s, c = runs[-1]
        if s == sym:
            if isinstance(c, int):
                nc = c + 1
                if nc >= T:
                    nc = ("T", nc % m)
            else:
                nc = ("T", (c[1] + 1) % m)
            return runs[:-1] + ((s, nc),)
    elif sym == 0:
        return runs  # blank sea
    return runs + ((sym, 1),)


def _pop(runs, T, m):
    """Pop one cell; yields (sym, new_runs) for every abstract possibility."""
    if not runs:
        yield 0, runs
        return
    s, c = runs[-1]
    if isinstance(c, int):
        if c == 1:
            yield s, runs[:-1]
        else:
            yield s, runs[:-1] + ((s, c - 1),)
    else:
        r = c[1]
        # remaining length may still be >= T
        yield s, runs[:-1] + ((s, ("T", (r - 1) % m)),)
        # or the length was exactly T (possible iff T = r mod m)
        if T % m == r:
            if T - 1 == 0:
                yield s, runs[:-1]
            else:
                yield s, runs[:-1] + ((s, T - 1),)


def decide_abstract_closure(tm, T=3, m=1, max_configs=50000, max_runs=48):
    """True iff the abstract reachable set is finite, within caps, and
    contains no halting configuration."""
    start = (0, (), ())  # state, left runs, right runs
    seen = {start}
    stack = [start]
    get = tm.get
    while stack:
        if len(seen) > max_configs:
            return False
        q, left, right = stack.pop()
        if len(left) + len(right) > max_runs:
            return False
        for s, right2 in _pop(right, T, m):
            tr = get(q, s)
            if tr is None or tr[2] == HALT:
                return False  # a halting configuration is abstractly reachable
            w, d, q2 = tr
            if d == R:
                nxt = (q2, _push(left, w, T, m), right2)
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
            else:
                mid = _push(right2, w, T, m)
                for y, left2 in _pop(left, T, m):
                    nxt = (q2, left2, _push(mid, y, T, m))
                    if nxt not in seen:
                        seen.add(nxt)
                        stack.append(nxt)
    return True


def decide_abstract(
    tm: TM,
    thresholds=(2, 3, 4),
    moduli=(1, 2, 3),
    max_configs=50000,
    block_sizes=(1, 2, 3),
):
    """Closure over (threshold, modulus, block size) combinations."""
    from .blocks import BlockMachine

    for k in block_sizes:
        machine = tm if k == 1 else BlockMachine(tm, k)
        for m in moduli:
            for T in thresholds:
                if decide_abstract_closure(
                    machine, T=T, m=m, max_configs=max_configs
                ):
                    return True
    return False
