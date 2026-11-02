"""Brady-style Tree Normal Form (TNF) enumeration of n-state 2-symbol machines.

A node of the tree is a partial transition table. The root fixes A0 -> 1RB
(write/direction/state-naming symmetry). A node is simulated from the blank
tape; the first undefined transition it reaches is (a) scored as the halt
transition (write 1, move R, halt) and (b) expanded with every possible
defined transition whose target is an already-used state or the next unused
one. A child that would leave no undefined transition can never halt and is
counted without recursion. Partial machines that exceed the node step budget
are completed with halt transitions and handed to the accelerated simulator
and the non-halting deciders; because reaching an undefined transition in the
completed machine *is* halting, a non-halting proof for the completion covers
the whole subtree.

Exactness caveat: a leaf reported "undecided" is a hole in the proof; the
driver reports these counts, and for n <= 4 they are zero with the default
budgets, which is what makes the reproduced S/Sigma values genuine.
"""
from __future__ import annotations

import random
from collections import Counter

from .machine import TM, HALT, R, L
from .simulate import simulate_accel, HALTED, NONHALT_ESCAPE
from .deciders import classify_nonhalt


def _simulate_partial(table, max_steps):
    """Run a partial table from blank tape.

    Returns ("undef", (q, s), steps, ones_if_halting) when an undefined
    transition is reached after `steps` steps (so the halting score is
    steps + 1), or ("limit", None, max_steps, None).
    """
    tape = {}
    q = 0
    pos = 0
    for steps in range(max_steps):
        s = tape.get(pos, 0)
        tr = table.get((q, s))
        if tr is None:
            ones = sum(1 for v in tape.values() if v == 1)
            if s == 0:
                ones += 1  # the halt step writes a 1 over the read 0
            return "undef", (q, s), steps, ones
        w, d, q2 = tr
        tape[pos] = w
        pos += d
        q = q2
    return "limit", None, max_steps, None


def _complete(table, n):
    trans = []
    for q in range(n):
        for s in (0, 1):
            trans.append(table.get((q, s)))
    return TM(n, trans)


def tnf_enumerate(
    n,
    node_steps=600,
    accel_macro=20000,
    accel_base=10**7,
    decider_cfg=None,
    rng: random.Random | None = None,
    max_leaves=None,
    collect_undecided=16,
):
    """Enumerate (fully, or randomly sampled when max_leaves is set).

    Returns a dict with: stats (Counter of leaf classes), best_steps,
    best_steps_tm, best_ones, best_ones_tm, leaves, undecided (sample of
    machine strings).
    """
    stats = Counter()
    best = {"steps": 0, "steps_tm": None, "ones": 0, "ones_tm": None}
    undecided_samples = []
    state = {"leaves": 0}

    def leaf(kind):
        stats[kind] += 1
        state["leaves"] += 1

    def rec(table, used, budget):
        if max_leaves is not None and state["leaves"] >= max_leaves:
            return
        kind, qs, steps, ones = _simulate_partial(table, budget)
        if kind == "limit":
            tm = _complete(table, n)
            status, bsteps, bones, _ = simulate_accel(
                tm, max_macro=accel_macro, max_base=accel_base
            )
            if status == HALTED:
                if bsteps <= 50 * budget:
                    # re-expand this node with a budget that reaches the slot
                    rec(table, used, int(bsteps) + 1)
                else:
                    leaf("halt_deep")  # sampling mode only; not expanded
                return
            if status == NONHALT_ESCAPE:
                leaf("nonhalt_escape")
                return
            label = classify_nonhalt(tm, decider_cfg)
            if label:
                leaf("nonhalt_" + label)
            else:
                leaf("undecided")
                if len(undecided_samples) < collect_undecided:
                    undecided_samples.append(tm.to_std())
            return
        # reached an undefined transition (q, s) after `steps` steps:
        # 1) score it as the halting machine
        leaf("halt")
        hsteps = steps + 1
        if hsteps > best["steps"]:
            t2 = dict(table)
            t2[qs] = None
            best["steps"], best["steps_tm"] = hsteps, _complete(t2, n).to_std()
        if ones > best["ones"]:
            t2 = dict(table)
            t2[qs] = None
            best["ones"], best["ones_tm"] = ones, _complete(t2, n).to_std()
        # 2) expand
        remaining_undef = 2 * n - len(table)
        hi = min(used, n - 1)  # new state allowed only as index `used`
        children = [
            (w, d, q2)
            for q2 in range(hi + 1)
            for d in (R, L)
            for w in (0, 1)
        ]
        if remaining_undef == 1:
            # every child machine has no undefined transition left -> no halt
            stats["fully_defined_nonhalt"] += len(children)
            return
        if rng is not None:
            rng.shuffle(children)
        for w, d, q2 in children:
            if max_leaves is not None and state["leaves"] >= max_leaves:
                return
            t2 = dict(table)
            t2[qs] = (w, d, q2)
            rec(t2, max(used, q2 + 1), node_steps)

    root = {(0, 0): (1, R, 1)}  # A0 -> 1RB, states A,B used
    rec(root, 2, node_steps)
    return {
        "stats": stats,
        "best_steps": best["steps"],
        "best_steps_tm": best["steps_tm"],
        "best_ones": best["ones"],
        "best_ones_tm": best["ones_tm"],
        "leaves": state["leaves"],
        "undecided": undecided_samples,
    }


def tnf_random_walks(
    n,
    walks,
    seed=0,
    node_steps=1000,
    accel_macro=20000,
    accel_base=10**6,
    decider_cfg=None,
    collect_undecided=24,
    min_defined=0,
):
    """Independent random root-to-leaf walks through the TNF tree.

    At every node that reaches an undefined transition, one leaf option is
    "this slot is the halt transition" and the others are its expansions;
    a walk picks uniformly among them. Leaves are classified exactly as in
    tnf_enumerate. This samples leaves independently (shallow leaves are
    over-represented relative to the whole tree, which is fine for a
    qualitative breakdown).
    """
    rng = random.Random(seed)
    stats = Counter()
    undecided_samples = []
    max_halt_steps = 0
    max_halt_tm = None
    for _ in range(walks):
        table = {(0, 0): (1, R, 1)}
        used = 2
        budget = node_steps
        while True:
            kind, qs, steps, ones = _simulate_partial(table, budget)
            if kind == "limit":
                tm = _complete(table, n)
                status, bsteps, bones, _ = simulate_accel(
                    tm, max_macro=accel_macro, max_base=accel_base
                )
                if status == HALTED:
                    if bsteps <= 10**6:
                        budget = int(bsteps) + 1
                        continue
                    stats["halt_deep"] += 1
                elif status == NONHALT_ESCAPE:
                    stats["nonhalt_escape"] += 1
                else:
                    label = classify_nonhalt(tm, decider_cfg)
                    if label:
                        stats["nonhalt_" + label] += 1
                    else:
                        stats["undecided"] += 1
                        if len(undecided_samples) < collect_undecided:
                            undecided_samples.append(tm.to_std())
                break
            remaining_undef = 2 * n - len(table)
            hi = min(used, n - 1)
            children = [
                (w, d, q2)
                for q2 in range(hi + 1)
                for d in (R, L)
                for w in (0, 1)
            ]
            allow_halt = len(table) >= min_defined
            pick = rng.randrange(len(children) + (1 if allow_halt else 0))
            if pick == len(children) or remaining_undef == 1:
                if pick == len(children):
                    stats["halt"] += 1
                    if steps + 1 > max_halt_steps:
                        t2 = dict(table)
                        t2[qs] = None
                        max_halt_steps = steps + 1
                        max_halt_tm = _complete(t2, n).to_std()
                else:
                    stats["fully_defined_nonhalt"] += 1
                break
            w, d, q2 = children[pick]
            table[qs] = (w, d, q2)
            used = max(used, q2 + 1)
            budget = node_steps
    return {
        "stats": stats,
        "undecided": undecided_samples,
        "max_halt_steps": max_halt_steps,
        "max_halt_tm": max_halt_tm,
    }
