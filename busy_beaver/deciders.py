"""Sound non-halting deciders: cyclers, translated cyclers, backward reasoning.

Each decider returns True only with a proof of non-halting; False means
"could not decide" (never "halts"). All three are standard techniques from
the bbchallenge BB(5) effort, reimplemented in simplified form:

- Cycler: the exact configuration (state, head position, tape) repeats, so the
  machine loops forever.
- Translated cycler: the configuration repeats up to a rightward (or, on the
  mirrored machine, leftward) translation at "record" times, when the head
  enters never-visited territory. If at two record times t1 < t2 with equal
  state the tape windows [l, p1] and [l + d, p2] agree, where d = p2 - p1 and
  l is the leftmost head position during [t1, t2], then everything ahead is
  blank in both configurations and the run from t2 replays the run from t1
  shifted by d, forever.
- Backward reasoning: explore all chains of predecessor configurations of
  every halting configuration; if every chain dies on a write/read
  contradiction within the depth budget, no halting configuration is
  reachable from anywhere, in particular not from the blank tape.
"""
from __future__ import annotations

from .machine import TM, HALT, R, L


def decide_cycler(tm: TM, max_steps=1000):
    """True iff an exact configuration repeat is found within max_steps."""
    tape = {}
    q = 0
    pos = 0
    seen = set()
    get = tm.get
    for _ in range(max_steps):
        key = (q, pos, frozenset((p, v) for p, v in tape.items() if v))
        if key in seen:
            return True
        seen.add(key)
        s = tape.get(pos, 0)
        tr = get(q, s)
        if tr is None or tr[2] == HALT:
            return False
        w, d, q2 = tr
        tape[pos] = w
        pos += d
        q = q2
    return False


def _run_record_trace(tm: TM, upto_steps):
    """Simulate, yielding (steps, q, pos, tape) at rightward record times."""
    tape = {}
    q = 0
    pos = 0
    rightmost = -1
    get = tm.get
    for steps in range(upto_steps):
        if pos > rightmost:
            rightmost = pos
            yield steps, q, pos, tape
        s = tape.get(pos, 0)
        tr = get(q, s)
        if tr is None or tr[2] == HALT:
            return
        w, d, q2 = tr
        tape[pos] = w
        pos += d
        q = q2


def _tc_verify(tm: TM, t1, p1, t2, p2, max_width=8192):
    """Re-simulate and check the translated-cycle condition for one pair."""
    d = p2 - p1
    if d <= 0:
        return False
    tape = {}
    q = 0
    pos = 0
    snap1 = None
    lmin = None
    get = tm.get
    for steps in range(t2 + 1):
        if steps == t1:
            snap1 = dict(tape)
            lmin = pos
        if snap1 is not None:
            lmin = min(lmin, pos)
        if steps == t2:
            if p1 - lmin + 1 > max_width:
                return False
            return all(
                snap1.get(i, 0) == tape.get(i + d, 0)
                for i in range(lmin, p1 + 1)
            )
        s = tape.get(pos, 0)
        tr = get(q, s)
        if tr is None or tr[2] == HALT:
            return False
        w, dd, q2 = tr
        tape[pos] = w
        pos += dd
        q = q2
    return False


def _tc_one_side(tm: TM, max_steps, window, max_verifs=12):
    # Pass 1: hash (state, last `window` cells) at record times; every
    # signature collision is a candidate, verified soundly by re-simulation.
    sigs = {}
    verifs = 0
    for steps, q, pos, tape in _run_record_trace(tm, max_steps):
        sig = (q, tuple(tape.get(pos - i, 0) for i in range(window)))
        if sig in sigs:
            t1, p1 = sigs[sig]
            verifs += 1
            if _tc_verify(tm, t1, p1, steps, pos):
                return True
            if verifs >= max_verifs:
                return False
        else:
            sigs[sig] = (steps, pos)
    return False


def decide_translated_cycler(tm: TM, max_steps=10000, window=48):
    """True iff a translated-cycle proof is found on either side."""
    return _tc_one_side(tm, max_steps, window) or _tc_one_side(
        tm.mirror(), max_steps, window
    )


def decide_backward(tm: TM, max_depth=12, max_nodes=50000):
    """True iff every halting configuration is provably unreachable.

    Works backward from each halting transition over partial tape windows
    (dict offset -> required symbol). A predecessor via transition
    (p, r) -> (w, dir, q) exists iff the cell the previous head wrote is
    compatible with w; the window then records that this cell held r before.
    If all branches contradict within max_depth, no halting configuration has
    an infinite (or start-reaching) past, so the blank-tape run never halts.
    """
    n = tm.n_states
    halt_entries = [
        (q, s)
        for q in range(n)
        for s in (0, 1)
        if tm.get(q, s) is None or tm.get(q, s)[2] == HALT
    ]
    if not halt_entries:
        return True  # no halt transition at all
    # Collect incoming transitions per target state.
    incoming = {q: [] for q in range(n)}
    for p in range(n):
        for r in (0, 1):
            tr = tm.get(p, r)
            if tr is not None and tr[2] != HALT:
                w, d, q2 = tr
                incoming[q2].append((p, r, w, d))
    frontier = [(q, 0, {0: s}) for (q, s) in halt_entries]
    nodes = 0
    for _depth in range(max_depth):
        new = []
        for q, h, win in frontier:
            # If this partial configuration is consistent with the start
            # configuration we cannot rule halting out.
            if q == 0 and all(v == 0 for v in win.values()):
                return False
            for p, r, w, d in incoming[q]:
                ph = h - d  # previous head position
                cur = win.get(ph)
                if cur is not None and cur != w:
                    continue  # the write contradicts the known tape
                nw = dict(win)
                nw[ph] = r
                new.append((p, ph, nw))
                nodes += 1
                if nodes > max_nodes:
                    return False
        if not new:
            return True  # all backward chains died: halting unreachable
        frontier = new
    return False


def classify_nonhalt(tm: TM, cfg=None):
    """Run the decider pipeline; return a label or None if undecided.

    cfg: dict with optional keys cycler_steps, tc_steps, tc_window,
    bw_depth, bw_nodes, ind_macro, ind_sym_steps.
    """
    from .inductive import decide_inductive
    from .abstract import decide_abstract

    cfg = cfg or {}
    if decide_cycler(tm, cfg.get("cycler_steps", 1000)):
        return "cycler"
    if decide_translated_cycler(
        tm, cfg.get("tc_steps", 10000), cfg.get("tc_window", 48)
    ):
        return "translated_cycler"
    if decide_backward(tm, cfg.get("bw_depth", 12), cfg.get("bw_nodes", 50000)):
        return "backward"
    if decide_inductive(
        tm,
        trace_macro=cfg.get("ind_macro", 4000),
        max_sym_steps=cfg.get("ind_sym_steps", 20000),
    ):
        return "inductive"
    if decide_abstract(
        tm,
        thresholds=cfg.get("abs_thresholds", (2, 3, 4)),
        moduli=cfg.get("abs_moduli", (1, 2, 3)),
        max_configs=cfg.get("abs_configs", 50000),
        block_sizes=cfg.get("abs_blocks", (1, 2, 3)),
    ):
        return "abstract_closure"
    return None
