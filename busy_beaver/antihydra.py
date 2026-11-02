"""Antihydra: the flagship BB(6) cryptid and its Collatz abstraction.

Machine (bbchallenge): 1RB1RA_0LC1LE_1LD1LC_1LA0LB_1LF1RE_---0RA

Known equivalence (bbchallenge wiki / OEIS A386792): with (a, b) = (8, 0),
iterate
    a even:  a <- 3a/2,      b <- b + 2
    a odd:   a <- (3a-1)/2,  b <- b - 1
Antihydra halts iff b ever reaches -1 (i.e. the orbit of a under
H(a) = floor(3a/2) eventually shows more odd terms than twice the even
terms). The parities behave empirically like fair coin flips, so b performs
a random walk with drift +1/2 per step; hitting -1 ever has small (and
rapidly decaying) probability, which is why the machine is "probviously"
non-halting - but nothing about the parity sequence of this single orbit is
provable with current number theory.

This module iterates the abstraction (exponentially faster than the tape:
one map step advances the machine through ~a tape sweeps) and provides the
cross-check data for the rule miner.
"""
from __future__ import annotations

from .machine import ANTIHYDRA, TM


def hydra_map(a):
    """One Antihydra step: returns (a', db)."""
    if a % 2 == 0:
        return 3 * a // 2, 2
    return (3 * a - 1) // 2, -1


def run_abstraction(n_iters, a0=8, keep=64):
    """Iterate the abstraction. Returns dict with trajectory stats.

    keep: how many initial a-values to retain for cross-checking.
    """
    a, b = a0, 0
    prefix = [a]
    b_min = b
    b_min_at = 0
    evens = odds = 0
    for i in range(1, n_iters + 1):
        a, db = hydra_map(a)
        if db == 2:
            evens += 1
        else:
            odds += 1
        b += db
        if b < b_min:
            b_min = b
            b_min_at = i
        if b == -1:
            return {
                "halted": True,
                "iters": i,
                "prefix": prefix[:keep],
            }
        if len(prefix) < keep:
            prefix.append(a)
    return {
        "halted": False,
        "iters": n_iters,
        "b": b,
        "b_min": b_min,
        "b_min_at": b_min_at,
        "evens": evens,
        "odds": odds,
        "a_bits": a.bit_length(),
        "prefix": prefix,
    }


def antihydra_tm() -> TM:
    return TM.from_std(ANTIHYDRA)
