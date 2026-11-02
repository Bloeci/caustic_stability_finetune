"""Automatic Collatz-abstraction mining from accelerated simulation traces.

Mechanism: the accelerated simulator sees the tape as runs. The *shape* of a
configuration is (state, left run symbols, right run symbols); the *values*
are the run lengths. Structured machines revisit a small set of shapes with
changing values, and the return map of a shape (values at one visit -> values
at the next visit) is, for a large class of machines, piecewise affine with
the branch chosen by a parity (or small modulus) test - exactly the shape of
a Collatz-like map. Mining that map from a trace turns tape-level simulation
into integer-map iteration and exposes the number-theoretic core of the
machine (for Antihydra: a -> 3a/2 resp. (3a-1)/2 on even/odd a).

Fitted rules are *verified against every recorded transition* of the shape,
so a reported rule is an empirical invariant of the whole observed trace, not
a guess from two points. It is still an induction hypothesis, not a proof;
proving it needs a (mechanisable) induction over the simulator semantics.
"""
from __future__ import annotations

from fractions import Fraction


class ShapeRecorder:
    """Recorder callback for simulate_accel; keeps value vectors per shape."""

    def __init__(self, max_runs=12, max_shapes=4000, max_occ=400):
        self.data = {}
        self.max_runs = max_runs
        self.max_shapes = max_shapes
        self.max_occ = max_occ

    def __call__(self, q, left, right, steps):
        if len(left) + len(right) > self.max_runs:
            return
        shape = (
            q,
            tuple(s for s, _ in left),
            tuple(s for s, _ in right),
        )
        lst = self.data.get(shape)
        if lst is None:
            if len(self.data) >= self.max_shapes:
                return
            lst = self.data[shape] = []
        if len(lst) < self.max_occ:
            lst.append(
                tuple(c for _, c in left) + tuple(c for _, c in right)
            )


def _fit_affine(xs, ys):
    """y = a*x + b over all samples, a,b rational; None if impossible."""
    base = None
    for x, y in zip(xs, ys):
        if base is None:
            base = (x, y)
        elif x != base[0]:
            a = Fraction(y - base[1], x - base[0])
            b = base[1] - a * base[0]
            if all(a * xi + b == yi for xi, yi in zip(xs, ys)):
                return a, b
            return None
    # all x equal; need all y equal too
    if len(set(ys)) == 1:
        return Fraction(0), Fraction(ys[0])
    return None


def _fmt_affine(a, b, var):
    if a == 0:
        return f"{b}"
    s = f"{a}*{var}" if a != 1 else var
    if b > 0:
        s += f" + {b}"
    elif b < 0:
        s += f" - {-b}"
    return s


def _fit_component(vins, outs):
    """Fit outs[k] as a function of vins[k] (vectors). Returns a rule string.

    Tries: constant; affine in one input; parity-split affine in one input
    (split on any pivot input's value mod 2).
    """
    m = len(vins[0])
    if len(set(outs)) == 1:
        return f"{outs[0]}"
    for i in range(m):
        xs = [v[i] for v in vins]
        fit = _fit_affine(xs, outs)
        if fit:
            return _fmt_affine(*fit, var=f"n{i}")
    for piv in range(m):
        even_idx = [k for k, v in enumerate(vins) if v[piv] % 2 == 0]
        odd_idx = [k for k, v in enumerate(vins) if v[piv] % 2 == 1]
        if len(even_idx) < 2 or len(odd_idx) < 2:
            continue
        for i in range(m):
            fe = _fit_affine(
                [vins[k][i] for k in even_idx], [outs[k] for k in even_idx]
            )
            fo = _fit_affine(
                [vins[k][i] for k in odd_idx], [outs[k] for k in odd_idx]
            )
            if fe and fo:
                return (
                    f"if n{piv} even: {_fmt_affine(*fe, var=f'n{i}')}; "
                    f"if n{piv} odd: {_fmt_affine(*fo, var=f'n{i}')}"
                )
    return None


def mine_rules(recorder: ShapeRecorder, min_occ=6):
    """Return list of (shape, n_occurrences, [rule per component], values)."""
    mined = []
    for shape, occ in recorder.data.items():
        if len(occ) < min_occ:
            continue
        vins = occ[:-1]
        vouts = occ[1:]
        m = len(occ[0])
        rules = []
        ok = True
        for j in range(m):
            r = _fit_component(vins, [v[j] for v in vouts])
            if r is None:
                ok = False
                break
            rules.append(f"n{j}' = {r}")
        if ok:
            mined.append((shape, len(occ), rules, occ))
    mined.sort(key=lambda t: -t[1])
    return mined


def shape_str(shape):
    q, lsyms, rsyms = shape
    state = chr(ord("A") + q)
    lt = " ".join(f"{s}^" for s in lsyms) or "(blank)"
    rt = " ".join(f"{s}^" for s in rsyms) or "(blank)"
    return f"[{lt}]  <{state}  [{rt}]"
