"""Turing machine representation: n states, 2 symbols, bbchallenge standard format.

Standard format example (6 states): "1RB1RA_0LC1LE_1LD1LC_1LA0LB_1LF1RE_---0RA"
Each underscore-separated block is one state (A, B, ...). Each block holds two
3-character entries: the action on reading 0, then on reading 1. An entry
"1RB" means: write 1, move Right, go to state B. "---" is an undefined
transition; reaching it halts the machine. "Z" or "H" as target state is an
explicit halt transition.

Scoring convention (Marxen-Buntrock / bbchallenge TNF): an undefined
transition behaves as "write 1, move Right, halt" and the halting step is
counted. Explicit halt transitions use their stated write/move. Under this
convention S(2)=6, S(3)=21, S(4)=107, S(5)=47176870.
"""
from __future__ import annotations

R = 1
L = -1
HALT = -1


class TM:
    __slots__ = ("n_states", "trans")

    def __init__(self, n_states, trans):
        # trans: tuple of length 2*n_states, indexed by 2*q + sym.
        # Entry: (write, dir, next_state) with next_state == HALT for explicit
        # halt, or None for undefined ('---').
        self.n_states = n_states
        self.trans = tuple(trans)

    def get(self, q, sym):
        return self.trans[2 * q + sym]

    @staticmethod
    def from_std(s: str) -> "TM":
        blocks = s.strip().split("_")
        trans = []
        for b in blocks:
            if len(b) != 6:
                raise ValueError(f"bad state block {b!r}")
            for e in (b[:3], b[3:]):
                if e == "---":
                    trans.append(None)
                    continue
                w = int(e[0])
                d = R if e[1] == "R" else L
                c = e[2]
                q2 = HALT if c in "ZH" else ord(c) - ord("A")
                trans.append((w, d, q2))
        return TM(len(blocks), trans)

    def to_std(self) -> str:
        out = []
        for q in range(self.n_states):
            b = ""
            for s in (0, 1):
                tr = self.trans[2 * q + s]
                if tr is None:
                    b += "---"
                else:
                    w, d, q2 = tr
                    b += f"{w}{'R' if d == R else 'L'}{'Z' if q2 == HALT else chr(ord('A') + q2)}"
            out.append(b)
        return "_".join(out)

    def mirror(self) -> "TM":
        """Left-right mirrored machine (same halting behaviour)."""
        trans = []
        for tr in self.trans:
            if tr is None:
                trans.append(None)
            else:
                w, d, q2 = tr
                trans.append((w, -d, q2))
        return TM(self.n_states, trans)

    def __repr__(self):
        return f"TM({self.to_std()})"


# Reference machines (strings verified against bbchallenge.org).
BB5_CHAMPION = "1RB1LC_1RC1RB_1RD0LE_1LA1LD_1RZ0LA"  # S(5)=47176870, Sigma(5)=4098
ANTIHYDRA = "1RB1RA_0LC1LE_1LD1LC_1LA0LB_1LF1RE_---0RA"  # 6-state cryptid
