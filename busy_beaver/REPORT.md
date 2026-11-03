# BB(6): acceleration mechanisms, decidability mechanisms, and the honest wall

Working notes from building this package. Everything labelled *measured* was
produced by code in this directory; run `python -m busy_beaver.run_all` to
reproduce (`--quick` skips the ~40-minute 4-state enumeration;
`python -m busy_beaver.tests.test_bb` runs the correctness tests).

## 1. What was built and what it establishes

| Component | File | Result (measured) |
|---|---|---|
| Accelerated simulator (run-length tape + chain steps, big-int counts) | `simulate.py` | BB(5) champion verified **exactly**: S = 47,176,870 steps, Σ = 4,098 ones, in 67,346 macro steps (700× compression, <0.1 s). Differentially tested against the naive simulator on 400 random machines. |
| TNF enumeration (Brady tree, interleaved simulation) | `enumerate_tm.py` | Reproduces S(2)=6, Σ(2)=4, S(3)=21, Σ(3)=6, S(4)=107, Σ(4)=13 from scratch. |
| Decider pipeline: escape, cycler, translated cycler, backward reasoning | `deciders.py` | Sound non-halting proofs; every positive verdict is re-verified by deterministic re-simulation. |
| Inductive (bouncer) prover: mined affine progressions proven by symbolic induction over run-length polynomials | `inductive.py` | Closes most bouncers, incl. two-sided ones (visit strides) and alternating patterns (via k-block machines). |
| Counter-abstraction closure (RepWL/CTL-style regular invariant with congruences) | `abstract.py` | Closes counters and parity-constrained bouncers. |
| Collatz rule miner (shape return-map fitting) | `rule_miner.py` | **Automatically recovers the hydra dynamics of Antihydra from its raw tape trace** (section 5). |
| Antihydra abstraction | `antihydra.py` | 100,000 hydra iterations in 1.3 s; the TM-level equivalent is ~10^35,000 base steps (a ≈ 10^17,600 and tape work per iteration scales like a²). |

Decider residue at n = 3 (16,549 TNF leaves), as technique families were added
— *measured*:

| Pipeline | Undecided leaves |
|---|---|
| escape + cycler + translated cycler + backward | 27 |
| + inductive prover | 16 |
| + k-block machines | 5 |
| + counter-abstraction closure (+ congruences) | 2 |

The last two (`1RB1LC_0LA0RB_1LA---`, `1RB---_1LC1RA_0RA0LC`) need *composed*
parametrized rules: an inner shift rule applied an n-dependent number of
times plus a reset rule (two-level induction). That is the full "bouncers"
decider of the Coq-BB5 proof — one more rung on a ladder that provably has no
top (section 3). Both machines are non-halting (all n ≤ 5 machines are decided
in the literature); "undecided" here means *my* pipeline produces no proof.

Full 4-state TNF enumeration (2,943,669 machines classified, 99 min pure
Python; light deciders inline, heavy pipeline as a parallel second pass over
the 1,992 survivors) — *measured*:

| Class | Machines |
|---|---|
| halting (champion `1RB1LB_1LA0LC_---1LD_1RD0RA`: **S(4)=107, Σ(4)=13** — both reproduced exactly) | 183,983 |
| fully defined, no halt transition | 2,323,408 |
| translated cycler | 309,856 |
| cycler | 59,026 |
| escape (endless chain into blank tape) | 53,437 |
| backward reasoning | 8,112 |
| inductive (bouncer) prover | 3,361 |
| counter-abstraction closure | 1,945 |
| **undecided** | **541** |

The 541 unproven machines (0.12 % of those needing a decider at all) are the
same two classes as the n = 3 stragglers — two-level bouncers and
coupled/binary counters — i.e. exactly the machines whose proofs need
composed parametrized rules or learned regular invariants (the bouncers and
FAR deciders of the BB(5) effort). All 4-state machines are decided in the
literature, so the S/Σ values above are exact; "undecided" measures the reach
of *this* pipeline, and the geometric-shrinkage pattern from n = 3 repeats
one level up, with the residue two orders of magnitude larger.

Random samples of the 6-state TNF tree (independent root-to-leaf walks) —
*measured*:

| Sample | Classified | Breakdown |
|---|---|---|
| 4,000 uniform walks | 4,000 / 4,000 | 633 halt; 1,077 escape; 1,974 translated cycler; 296 cycler; 19 backward; 1 inductive |
| 1,500 deep walks (≥10 of 12 transitions defined before a halt leaf is allowed) | 1,499 / 1,500 | 1 halt; 503 escape; 850 translated cycler; 140 cycler; 4 backward; 1 fully-defined non-halter; **1 undecided** |

The one undecided machine, `1RB---_1LC0RD_0LC1RB_0RB---_------_------`, is a
binary-counter variant: its tape is an evolving binary numeral with carries,
its halting transition sits on a carry pattern that never occurs, and
neither congruence-refined closure (up to T=6, m=3, 4-blocks, 300k abstract
configs) nor backward reasoning to depth 32 proves it. It needs a learned
regular invariant (FAR-class decider). This is the shape of the whole BB(6)
problem in miniature: random sampling finds almost nothing hard — the
difficulty is concentrated in a vanishingly small set of machines, and each
survivor demands a strictly stronger proof technique than the last.

## 2. Acceleration: how far cleverness actually goes

Mechanisms, in increasing order of power; the first two are implemented here,
and each is exponential over the previous one on the machines that matter.

1. **Chain steps over run-length tape** (`simulate.py`). A transition that
   re-enters its own state crosses an entire run in O(1). Measured: 700× on
   the BB(5) champion; 17,750× and growing on Antihydra (53·10^9 base steps
   in 3·10^6 macro steps), because run lengths grow exponentially.
2. **k-block macro machines** (`blocks.py`). Group k cells into one symbol so
   alternating patterns become runs; composes with everything above.
3. **Inductive rules / formula tapes** (`inductive.py` is the one-variable
   case). Replace concrete run lengths by polynomials in n and prove
   C(n) →* C(n+1) once, symbolically. This is no longer simulation but proof:
   it covers *infinitely many* steps per rule. The 10↑↑15 and 2↑↑↑5 BB(6)
   lower bounds are established this way — nobody ran those machines.
4. **Abstraction compilation** (`rule_miner.py` automates the discovery).
   Compile the machine's return map into an integer map and iterate that:
   one hydra iteration ≈ a² tape steps for Antihydra. Each level of
   compilation converts an exponential into a linear cost.

The limit: acceleration changes *which* machines are feasible to analyze,
never *whether* the halting question is answerable. For the BB(6) champion
(≈ 2↑↑↑5 steps) even writing the step count in decimal is impossible; its
score is *proven*, not simulated. Acceleration is necessary scaffolding for
the proofs, and it is the easy half.

## 3. Decidability: every decider is a proof-system fragment

A "decider" is a sound, incomplete procedure proving non-halting. The ladder
implemented here — loop detection, translated-cycle detection, backward
reasoning, affine induction, regular closure with congruences — mirrors the
real history: each family is an invariant language (finite configurations;
shifted configurations; co-reachability; affine arithmetic progressions;
regular languages with counters mod m). The measured residues (27 → 16 → 5 →
2 at n = 3) show the characteristic pattern: **geometric shrinkage, never
zero.**

That is a theorem, not bad luck. The halting problem is undecidable, so no
recursively enumerable family of sound deciders is complete; worse, by
Conway (1972), even the *target* of my rule miner — generalized Collatz maps,
piecewise-affine with modular branching — is Turing-complete, so the
"compile to an integer map, then decide the map" pipeline cannot be a
decision procedure either. Every ladder of techniques leaves holdouts; the
only question is whether the *finite* set of holdouts at a given n happens to
consist of problems mathematics can solve.

At n = 5 the answer was yes: ~181 million TNF machines, closed completely by
automated deciders plus a finite list of hand proofs, machine-checked in Coq
(2024). At n = 6 the answer is currently no, and provably not for lack of
compute (section 4).

For precision: BB(6) is a single integer, so it is "computable" in the
trivial sense. The actual task is to *prove* each non-halting machine
non-halting in a formal system we trust. BB(745) is independent of ZFC; BB(5)
is settled; where between 5 and 745 the provability boundary sits is open,
and it is a live possibility that 6 is already past it.

## 4. The mathematician's attack plan, and where each line dies

What I would do with unlimited ambition, stated as a program:

**(a) Compile away the tape.** Run the miner over every holdout; the tape is
scaffolding, the content is an integer map. Measured proof of concept: the
miner extracts Antihydra's map automatically (section 5). Classify the
compiled maps: bounded/Presburger-flat systems are decidable (decide them);
affine maps with modular branching are Collatz-class (isolate them). This
pipeline would shrink the n = 6 frontier to a catalog of number-theory
problems — the bbchallenge "cryptids". *Dies at:* Conway's theorem
guarantees the catalog is nonempty for some n, and empirically it is
nonempty at n = 6.

**(b) Attack the cryptids with drift arguments.** Antihydra halts iff a
functional of the parity sequence of the orbit a ← ⌊3a/2⌋ hits a barrier.
If parities were i.i.d. fair coins, b would be a random walk with drift +1/2;
the probability of ever hitting −1 from 0 solves q = ½ + ½q³, i.e.
q = (√5−1)/2 ≈ 0.618 a priori — but *conditioned on the measured prefix*
(100,000 iterations, b = 50,477, never below 0) the surviving halting
probability is q^50478 ≈ 10^−10,549. So the answer is "known" with absurd
confidence. *Dies at:* the parity sequence of one specific orbit is not a
random variable. Proving anything about it (even "infinitely many even
terms" is nontrivial for such orbits) is the same wall as Mahler's Z-number
problem (open since 1968) and Collatz itself: we have measure-theoretic
statements about almost all orbits and nothing about any individual one.
A proof would be a number-theory breakthrough whose value dwarfs BB(6).

**(c) Search for invariants instead of trajectories.** Non-halting ⟺ there
exists a closed invariant set containing the start and avoiding halts. My
`abstract.py` searches regular invariants with congruence refinement;
bbchallenge's FAR/WFAR search (weighted) automata invariants by SAT. The
honest generalization is a hierarchy: regular → register/semilinear →
polynomial (Positivstellensatz / SOS certificates). Invariant search is the
practical frontier and is what actually finished BB(5). *Dies at:* an
invariant proof of Antihydra would encode a parity-equidistribution fact as
a finite certificate; no finite-state abstraction can track ⌊3a/2⌋ parity
exactly (the parity of a depends on the full binary expansion, which the map
shifts and carries), so the certificate would have to come from structure
nobody has found. Worth automating; not expected to close the list.

**(d) Fall back to independence.** If (b) and (c) fail, prove that they must:
show BB(6) (or a specific cryptid) is independent of PA or ZFC, the way
Yedidia–Aaronson/O'Rear/Riebel did for n = 745. Nothing rules out n = 6
being independent of PA; nobody has a technique for proving independence of
a *natural* small machine either. *Dies at:* current independence proofs
work by embedding a proof-searcher for the theory into the machine, which
needs hundreds of states. A 6-state independence proof needs fundamentally
new ideas about what tiny programs can encode.

Net assessment, as the person who just built the ladder: BB(6) will not fall
to a clever trick. It is lower-bounded by ≥ 1 open Collatz-class problem
(Antihydra; the wiki lists further 6-state cryptids, e.g. hydra-function
relatives), each individually resistant for structural reasons, and the
believed value (the 2↑↑↑5-level champion score) is already known — what is
missing is exactly the part cleverness cannot shortcut: proofs about
individual exponential orbits.

## 5. Antihydra case study (all measured)

Machine: `1RB1RA_0LC1LE_1LD1LC_1LA0LB_1LF1RE_---0RA`.

- Accelerated trace: 53,249,891,725 base steps in 3·10^6 macro steps (6 s).
- Miner output, with no knowledge of the literature, from the trace alone —
  shape `[1^n0 0^n1 1^n2] <E` (head at left end, 26 visits):

  ```
  n0' = if n2 even: n0 + 2;        if n2 odd: n0 - 1
  n1' = 1
  n2' = if n2 even: 3/2*n2 + 2;    if n2 odd: 3/2*n2 + 3/2
  ```

  This *is* the hydra map (an affine conjugate of a ← ⌊3a/2⌋; starting value
  8 appears as the first visit) together with the halting counter: the left
  1-run tracks b (+2 on even, −1 on odd, offset 2), and the machine reaches
  its undefined transition exactly when that run is exhausted. Automatic
  recovery of the known reduction (bbchallenge wiki / OEIS A386792).
- Abstraction run: 100,000 iterations, 1.3 s; never halted; b = 50,477,
  minimum 0; parities 50,159 even / 49,841 odd; a has 58,500 bits. The
  bbchallenge project has pushed this vastly further with no halt.
- Pseudo-randomness of the orbit (200,000 iterations, measured): even
  fraction 0.50018; parity run lengths fit the geometric(1/2) law
  (49,953/25,054/12,663/6,137 observed vs 50,034/25,017/12,509/6,254
  expected for lengths 1-4); lag-1 autocorrelation −0.0007. The orbit is
  statistically indistinguishable from the fair coin the drift argument
  assumes — and proving *any* of these statistics rigorously for this one
  deterministic orbit is the open problem.

The pair (trace → mined map → fast abstraction → drift analysis) is the
complete modern workflow; the step that does not exist is the last one:
turning the drift argument into a proof.

## 6. What would constitute real progress

1. Implement composed-rule induction (two-level bouncers) and multivariate
   closure — mechanical work; closes the known-boring residue (my 2 leftover
   n=3 machines, most of n=6's non-cryptid holdouts).
2. Invariant synthesis via SAT/SMT over weighted automata and semilinear
   sets, CEGAR-style — the highest-leverage engineering.
3. Any theorem about the parity sequence of a single ⌊3a/2⌋ orbit — a
   mathematics paper regardless of BB.
4. A formalized pipeline (as Coq-BB5 did) so that "BB(6) = value of the
   champion, conditional on cryptid list L" becomes a machine-checked
   statement. That conditional theorem is realistically the best endpoint:
   BB(6) reduced to named open problems.

## 7. Sources

- bbchallenge wiki: [Antihydra](https://wiki.bbchallenge.org/wiki/Antihydra),
  [Hydra function](https://wiki.bbchallenge.org/wiki/Hydra_function),
  [BB(6) champion listing](https://wiki.bbchallenge.org/wiki/Champion) (2↑↑↑5-level halter, mxdys 2025),
  [antihydra page on bbchallenge.org](https://bbchallenge.org/antihydra)
- [OEIS A386792](https://oeis.org/A386792) — the Antihydra sequence (a, b).
- [Coq-BB5: Determination of the fifth Busy Beaver value](https://arxiv.org/pdf/2509.12337)
  (machine-checked S(5) = 47,176,870; decider taxonomy).
- [Turing machine deciders, part I](https://arxiv.org/pdf/2504.20563).
- J. H. Conway, *Unpredictable Iterations* (1972): generalized Collatz maps
  are Turing-complete (undecidability of the compiled-map pipeline).
- S. Aaronson, *The Busy Beaver Frontier* (2020) — survey; independence of
  BB(745) from ZFC follows Yedidia–Aaronson → O'Rear → Riebel.

All machine strings above were cross-checked against bbchallenge sources
before use; the BB(5) champion string is additionally self-verifying (the
simulator reproduces the published step and ones counts exactly).
