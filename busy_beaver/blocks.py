"""k-block macro machines (Marxen-Buntrock).

Groups k tape cells into one symbol (an integer 0..2^k-1, MSB = leftmost
cell; 0 is the blank block). Macro states are 2*q + side, where side says
whether the head enters the block at its leftmost (0) or rightmost (1) cell.
A macro transition simulates the underlying 2-symbol machine inside the
block until the head leaves it. If the underlying machine halts inside the
block, or never leaves it, the macro transition is reported as a halt
(None): sound for every non-halting proof built on top (the proof simply
fails), and cyclers-inside-a-block are caught by the plain cycler decider.

The object exposes .get(state, symbol) and duck-types as a TM for the
accelerated simulator and the inductive prover, which turns alternating
patterns like (10)^n - hopeless for run-length encoding at block size 1 -
into single growing runs.
"""
from __future__ import annotations

from .machine import TM, HALT, R, L


class BlockMachine:
    def __init__(self, tm: TM, k: int):
        self.tm = tm
        self.k = k
        self.n_states = 2 * tm.n_states
        # upper bound on distinct in-block configurations; beyond it the
        # underlying machine provably loops inside the block
        self.max_inner = (2 ** k) * k * tm.n_states + 2
        self.cache = {}

    def get(self, mq, word):
        key = (mq, word)
        hit = self.cache.get(key, "MISS")
        if hit != "MISS":
            return hit
        k = self.k
        q, side = divmod(mq, 2)
        cells = [(word >> (k - 1 - i)) & 1 for i in range(k)]
        pos = 0 if side == 0 else k - 1
        tm = self.tm
        result = None
        for _ in range(self.max_inner):
            tr = tm.get(q, cells[pos])
            if tr is None or tr[2] == HALT:
                result = None  # halts inside the block
                break
            w, d, q2 = tr
            cells[pos] = w
            pos += d
            q = q2
            if pos < 0:
                new_word = 0
                for c in cells:
                    new_word = (new_word << 1) | c
                result = (new_word, L, 2 * q + 1)  # enters next block from its right
                break
            if pos >= k:
                new_word = 0
                for c in cells:
                    new_word = (new_word << 1) | c
                result = (new_word, R, 2 * q + 0)  # enters next block from its left
                break
        # loop exhausted without break -> loops inside block -> None (sound)
        self.cache[key] = result
        return result
