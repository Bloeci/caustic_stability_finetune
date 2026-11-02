"""Naive and accelerated (run-length + chain-step) Turing machine simulation.

The accelerated simulator represents the tape as two stacks of runs
[symbol, count] (top of stack = nearest the head); the current cell is the
first cell of the top run of the right stack, and both stacks sit on an
implicit infinite sea of 0s. A transition that re-enters its own state while
moving over a run of the symbol it reads is executed over the whole run in one
"macro step" (Marxen-Buntrock chain step), so counts (and step counts) are big
integers and sweeps cost O(1) per run instead of O(run length).

A chain step into the blank sea proves non-halting ("escape"): the machine
moves in one direction forever reading 0 in a fixed state.
"""
from __future__ import annotations

from .machine import TM, HALT, R, L

HALTED = "halt"
RUNNING = "running"
NONHALT_ESCAPE = "nonhalt_escape"


def simulate_naive(tm: TM, max_steps: int):
    """Plain dict-tape simulation.

    Returns (status, steps, ones, tape, state, pos). ones is None unless
    halted.
    """
    tape = {}
    q = 0
    pos = 0
    steps = 0
    get = tm.get
    while steps < max_steps:
        s = tape.get(pos, 0)
        tr = get(q, s)
        if tr is None or tr[2] == HALT:
            w, d = (1, R) if tr is None else (tr[0], tr[1])
            tape[pos] = w
            steps += 1
            ones = sum(1 for v in tape.values() if v == 1)
            return HALTED, steps, ones, tape, q, pos + d
        w, d, q2 = tr
        tape[pos] = w
        pos += d
        q = q2
        steps += 1
    return RUNNING, steps, None, tape, q, pos


def _push(stack, sym, cnt):
    """Push cnt cells of sym; merge with top run; drop 0-runs at the sea edge."""
    if stack:
        top = stack[-1]
        if top[0] == sym:
            top[1] += cnt
            return
    elif sym == 0:
        return  # merging into the blank sea
    stack.append([sym, cnt])


def _pop_one(stack):
    """Remove one cell from the top run; return its symbol (0 if blank sea)."""
    if not stack:
        return 0
    top = stack[-1]
    sym = top[0]
    if top[1] == 1:
        stack.pop()
    else:
        top[1] -= 1
    return sym


def _ones(left, right):
    return sum(c for s, c in left if s == 1) + sum(c for s, c in right if s == 1)


def simulate_accel(tm: TM, max_macro=1_000_000, max_base=None, recorder=None):
    """Accelerated simulation.

    Returns (status, base_steps, ones, macro_steps).
    status: HALTED (ones set), RUNNING (budget exhausted), NONHALT_ESCAPE
    (proof of non-halting: endless chain into blank tape).
    recorder, if given, is called as recorder(q, left, right, base_steps)
    at the start of every macro step (used by the rule miner).
    """
    left = []   # runs away from head; top of stack = end of list
    right = []  # top run's first cell is the current cell
    q = 0
    steps = 0
    macro = 0
    get = tm.get
    while macro < max_macro and (max_base is None or steps < max_base):
        macro += 1
        if recorder is not None:
            recorder(q, left, right, steps)
        s = right[-1][0] if right else 0
        tr = get(q, s)
        if tr is None or tr[2] == HALT:
            w = 1 if tr is None else tr[0]
            _pop_one(right)
            _push(right, w, 1)
            steps += 1
            return HALTED, steps, _ones(left, right), macro
        w, d, q2 = tr
        if d == R:
            if not right:
                if q2 == q:
                    return NONHALT_ESCAPE, steps, None, macro
                _push(left, w, 1)
                steps += 1
            elif q2 == q:
                sym, k = right.pop()
                _push(left, w, k)
                steps += k
            else:
                _pop_one(right)
                _push(left, w, 1)
                steps += 1
        else:  # move Left
            _pop_one(right)  # consume current cell
            k_extra = 0
            if q2 == q:
                if left and left[-1][0] == s:
                    k_extra = left.pop()[1]
                elif not left and s == 0:
                    return NONHALT_ESCAPE, steps, None, macro
            _push(right, w, 1 + k_extra)
            steps += 1 + k_extra
            y = _pop_one(left)  # step onto next cell to the left
            _push(right, y, 1)
        q = q2
    return RUNNING, steps, None, macro
