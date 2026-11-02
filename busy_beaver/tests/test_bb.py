"""Correctness tests. Run:  python -m busy_beaver.tests.test_bb"""
from __future__ import annotations

import random
import sys

from ..machine import TM, BB5_CHAMPION, ANTIHYDRA, R, L, HALT
from ..simulate import simulate_naive, simulate_accel, HALTED, RUNNING, NONHALT_ESCAPE
from ..deciders import decide_cycler, decide_translated_cycler, decide_backward
from ..enumerate_tm import tnf_enumerate
from ..antihydra import run_abstraction


def random_tm(rng, n):
    trans = []
    for _ in range(2 * n):
        if rng.random() < 0.08:
            trans.append(None)
        else:
            trans.append(
                (rng.randint(0, 1), rng.choice((R, L)), rng.randrange(n))
            )
    return TM(n, trans)


def test_differential():
    """Accelerated simulator must agree with the naive one."""
    rng = random.Random(42)
    checked = 0
    for _ in range(400):
        tm = random_tm(rng, rng.choice((2, 3, 4)))
        st_n, steps_n, ones_n, _, _, _ = simulate_naive(tm, 3000)
        st_a, steps_a, ones_a, _ = simulate_accel(tm, max_macro=20000, max_base=3000)
        if st_n == HALTED:
            assert st_a == HALTED, (tm.to_std(), st_a)
            assert steps_a == steps_n, (tm.to_std(), steps_a, steps_n)
            assert ones_a == ones_n, (tm.to_std(), ones_a, ones_n)
            checked += 1
        elif st_a == HALTED:
            # accel may overshoot max_base within one chain; the halt must
            # then be beyond the naive budget
            assert steps_a >= 3000, tm.to_std()
    assert checked > 50, f"only {checked} halting machines sampled"
    print(f"differential test ok ({checked} halting machines compared)")


def test_escape_detection():
    # A0: write 1 move R stay A -> runs right forever over blanks.
    tm = TM(2, [(1, R, 0), None, None, None])
    st, _, _, _ = simulate_accel(tm, max_macro=100)
    assert st == NONHALT_ESCAPE
    print("escape detection ok")


def test_deciders():
    # Cycler: A0->1RB, B0->1LA writes then returns; A1->1RB, B1->1LA ->
    # ping-pong on two cells forever.
    tm = TM(2, [(1, R, 1), (1, R, 1), (1, L, 0), (1, L, 0)])
    assert decide_cycler(tm, 100)
    # Translated cycler: A0->1RB, B0->1RA : marches right writing 1s,
    # alternating states (no self-chain, so escape detection does not fire).
    tm2 = TM(2, [(1, R, 1), None, (1, R, 0), None])
    assert decide_translated_cycler(tm2, 1000, 16)
    # Backward: halt only on B reading 1. The only way into B is A0, which
    # writes 0 and moves R, so the cell left of the halt cell held 0; the only
    # ways into A (A1, B0) write 0 onto the halt cell itself and move L, which
    # contradicts the required 1 there. All backward chains die at depth 2.
    tm3 = TM(2, [(0, R, 1), (0, L, 0), (0, L, 0), None])
    assert decide_backward(tm3, 8, 10000)
    print("deciders ok")


def test_small_bb_values():
    r2 = tnf_enumerate(2, node_steps=200)
    assert r2["best_steps"] == 6 and r2["best_ones"] == 4, r2
    assert r2["stats"]["undecided"] == 0
    r3 = tnf_enumerate(3, node_steps=400)
    assert r3["best_steps"] == 21 and r3["best_ones"] == 6, (
        r3["best_steps"],
        r3["best_ones"],
    )
    # Two second-level bouncers defeat the implemented decider families; they
    # need composed parametrized rules (full "bouncers" decider). Both are
    # proven non-halting in the literature (Coq-BB5 decides all of n<=5).
    assert r3["stats"]["undecided"] <= 2, dict(r3["stats"])
    print(
        "BB(2)=6/Sigma=4 and BB(3)=21/Sigma=6 reproduced "
        f"({r3['stats']['undecided']} known-nonhalting leaves unproven)"
    )


def test_no_false_nonhalt_proofs():
    """The whole decider pipeline must stay silent on halting machines."""
    from ..deciders import classify_nonhalt
    from ..machine import BB5_CHAMPION

    assert classify_nonhalt(TM.from_std(BB5_CHAMPION)) is None
    assert classify_nonhalt(TM.from_std("1RB1LB_1LA0LC_1RZ1LD_1RD0RA")) is None
    rng = random.Random(123)
    tested = 0
    for _ in range(600):
        tm = random_tm(rng, rng.choice((2, 3, 4)))
        st, steps, _, _, _, _ = simulate_naive(tm, 5000)
        if st == HALTED:
            lab = classify_nonhalt(tm)
            assert lab is None, (tm.to_std(), lab, steps)
            tested += 1
    print(f"soundness ok: {tested} halting machines, no false non-halting proof")


def test_hydra_prefix():
    # OEIS A386792: a runs 8, 12, 18, 27, 40, ...
    res = run_abstraction(200, keep=5)
    assert res["prefix"] == [8, 12, 18, 27, 40], res["prefix"]
    assert not res["halted"]
    print("Antihydra abstraction prefix matches OEIS A386792")


def main():
    test_differential()
    test_escape_detection()
    test_deciders()
    test_small_bb_values()
    test_no_false_nonhalt_proofs()
    test_hydra_prefix()
    print("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
