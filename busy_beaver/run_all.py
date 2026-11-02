"""Experiment driver. Run:  python -m busy_beaver.run_all [--quick]

Sections:
1. Known-value validation: full TNF enumeration for n = 2, 3 (and 4 unless
   --quick) must reproduce S/Sigma, with the leaf classification breakdown.
2. BB(5) champion: exact accelerated verification of S(5) = 47,176,870.
3. Antihydra: abstraction run, automatic Collatz-rule mining from the raw
   tape trace, and the cross-check between the two.
4. BB(6) frontier: random sample of the 6-state TNF tree classified by the
   decider pipeline.
"""
from __future__ import annotations

import sys
import time

from .machine import TM, BB5_CHAMPION, ANTIHYDRA
from .simulate import simulate_accel
from .enumerate_tm import tnf_enumerate, tnf_random_walks
from .rule_miner import ShapeRecorder, mine_rules, shape_str
from .antihydra import run_abstraction


def section_validation(quick):
    print("== 1. TNF enumeration vs known values ==")
    expect = {2: (6, 4), 3: (21, 6), 4: (107, 13)}
    for n in (2, 3) if quick else (2, 3, 4):
        t0 = time.time()
        r = tnf_enumerate(n, node_steps=600)
        s, o = expect[n]
        ok = "OK" if (r["best_steps"] == s and r["best_ones"] == o) else "MISMATCH"
        print(
            f"n={n}: S={r['best_steps']} Sigma={r['best_ones']} [{ok}] "
            f"({time.time()-t0:.0f}s)"
        )
        print(f"   leaf classes: {dict(r['stats'])}")
        if r["undecided"]:
            print(f"   unproven-nonhalt leaves: {r['undecided']}")


def section_bb5():
    print("== 2. BB(5) champion ==")
    tm = TM.from_std(BB5_CHAMPION)
    t0 = time.time()
    status, steps, ones, macro = simulate_accel(tm, max_macro=10**7)
    assert status == "halt" and steps == 47176870 and ones == 4098
    print(
        f"{BB5_CHAMPION}: S=47,176,870 Sigma=4,098 verified in {macro:,} "
        f"macro steps ({steps/macro:.0f}x chain compression, "
        f"{time.time()-t0:.2f}s)"
    )


def section_antihydra():
    print("== 3. Antihydra ==")
    tm = TM.from_std(ANTIHYDRA)
    rec = ShapeRecorder(max_runs=8)
    t0 = time.time()
    status, steps, _, macro = simulate_accel(tm, max_macro=3_000_000, recorder=rec)
    print(
        f"TM trace: {steps:,} base steps in {macro:,} macro steps "
        f"({time.time()-t0:.0f}s), status={status}"
    )
    for shape, nvis, rules, occ in mine_rules(rec, min_occ=8)[:3]:
        print(f"mined shape {shape_str(shape)} ({nvis} visits):")
        for r in rules:
            print(f"    {r}")
        print(f"    visits: {occ[:8]}")
    t0 = time.time()
    res = run_abstraction(100_000)
    print(
        f"abstraction: {res['iters']:,} hydra iterations ({time.time()-t0:.1f}s), "
        f"halted={res['halted']}, b={res['b']}, b_min={res['b_min']}, "
        f"evens/odds={res['evens']}/{res['odds']}, a has {res['a_bits']:,} bits"
    )
    print(f"a_n prefix: {res['prefix'][:12]}")


def section_bb6_sample(walks=1000):
    print(f"== 4. BB(6) frontier: {walks} random TNF leaves ==")
    cfg = {
        "tc_steps": 5000,
        "abs_thresholds": (2, 3),
        "abs_moduli": (1, 2),
        "abs_blocks": (1, 2),
        "abs_configs": 20000,
        "ind_macro": 2000,
    }
    for label, min_def in (("uniform walks", 0), ("deep walks (>=10 defined)", 10)):
        t0 = time.time()
        r = tnf_random_walks(
            6, walks=walks, seed=7, node_steps=800,
            decider_cfg=cfg, min_defined=min_def,
        )
        print(f"{label} ({time.time()-t0:.0f}s): {dict(r['stats'])}")
        print(
            f"   longest halter sampled: {r['max_halt_steps']} steps  "
            f"{r['max_halt_tm']}"
        )
        if r["undecided"]:
            print("   undecided sample:")
            for m in r["undecided"][:10]:
                print("      ", m)


def main():
    quick = "--quick" in sys.argv
    section_validation(quick)
    section_bb5()
    section_antihydra()
    section_bb6_sample(400 if quick else 1000)


if __name__ == "__main__":
    main()
