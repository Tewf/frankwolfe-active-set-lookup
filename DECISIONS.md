# Decisions awaiting Mohamed

Nothing in this repository was posted anywhere. This file is every place a
human judgement call was made or is still needed, so any of it can be
revisited.

## The draft comment for issue #244

**Posting this is Mohamed's to do, not automated, not implied by anything
else in this repository.** Ready to paste as-is if the numbers still look
right on a re-read.

**This replaces an earlier draft that concluded no advantage.** That draft
only ever tested hashing the whole atom against a miss query; both
restrictions turned out to matter, and relaxing them reverses the
conclusion. The earlier text is in git history (`git log -p DECISIONS.md`),
not repeated here.

> I measured this rather than guessing, and my first pass got the wrong
> answer, so here is the corrected one. `find_atom`'s linear scan is called
> whenever `active_set_update!` isn't given an index, which turns out to be
> more call sites than I expected: every "add to active set" step of BPCG
> (`blended_pairwise.jl:374`), every pairwise step of plain PFW
> (`pairwise.jl:242`, since `pfw_step` never supplies one), BCG
> (`blended_cg.jl:358`), and the generic corrective/block-coordinate drivers
> (`corrective_frankwolfe.jl:240`, `block_coordinate_algorithms.jl:421`),
> but never away-step FW, which always tracks its own index.
>
> On three real BPCG runs (Birkhoff polytope at n=25 and n=60, the L∞-ball
> at d=3,000; code at github.com/Tewf/frankwolfe-active-set-lookup), active
> sets topped out at 158-389 atoms and `find_atom` took 0.02%-0.14% of total
> runtime (`measurement/results.csv`). I first checked whether hashing the
> *whole* atom would help (linear scan vs. `Dict`, sizes 1-20k, dims
> 16-8,192) and it doesn't: hashing an atom is always O(dimension), so a
> full-atom `Dict` only wins once the scan's own O(dimension) worst case is
> actually being paid, which these three runs never came close to.
>
> That was the wrong question, though: nobody has to hash the *whole* atom.
> A hash over just the first `k` coordinates costs O(k), independent of
> dimension, and is exactly as sound as a full-atom hash, since a bucket
> hit still gets confirmed against the whole atom with the same exact
> `_unsafe_equal` before being trusted; a shorter hash can only cost speed,
> never correctness. I timed that too (`microbenchmark/run_prefix.jl`),
> against `FrankWolfe.jl`'s own atom shapes this time, not generic random
> vectors: Birkhoff permutation matrices and L∞-ball box corners, generated
> by calling `BirkhoffPolytopeLMO`/`LpNormBallLMO{Inf}` directly, and
> against a hit query (the atom already present, forcing the scan through
> the whole match) as well as a miss. A `k=8` prefix hash beats the scan for
> every alphabet and query mix tested, at active-set sizes below every one
> of the three real runs' own maximum
> (`microbenchmark/results_prefix_crossover.csv`): 500.6ns vs. 949.4ns on a
> miss at active-set size 158 (Birkhoff n=25's own maximum), 350.5ns vs.
> 488.2ns on a hit at the same size.
>
> So: contrary to what I first wrote here, a prefix-hashed active set would
> help, at exactly the active-set sizes this workload reaches.
> `active_set_cleanup!`'s `deleteat!` shifts every index past the one
> removed, so a hash→index map needs upkeep on every drop, not just every
> push, the same cost a full-atom hash would have paid. Would you want a
> hash-augmented default `ActiveSet` keyed by a short prefix, or a separate
> subtype the way `ActiveSetQuadraticProductCaching`
> (`active_set_quadratic.jl`) already is? I'd lean subtype: the right `k`
> likely depends on the polytope (permutation matrices and box corners both
> needed `k=8` here; a higher-entropy vertex might need less, a more
> degenerate one more), which argues for something tunable rather than a
> fixed default.

## The draft links a private repository

The comment above cites `github.com/Tewf/frankwolfe-active-set-lookup` as where
the code lives. **This repository is private**, so that link 404s for everyone
who reads the issue, which is worse than no link at all.

Three ways out, and this one has to be settled before the comment is posted:

- **Make it public first.** It is a measurement of a public library answering a
  public issue, and there is nothing in it that is not already sayable in the
  comment. This is the option that makes the comment strongest, since the
  numbers become checkable by the maintainers rather than asserted.
- **Drop the link** and let the comment stand on the numbers alone. It still
  reads as work rather than opinion, but nobody can re-run it.
- **Post the numbers, offer the code.** "Happy to share the harness if useful":
  keeps the repository private and puts the choice on them.

**Recommendation: make it public before posting.** A measurement whose code
cannot be inspected is an assertion, and the whole point of answering a
five-year-old question this way is that the answer can be checked. The
repository contains no personal data and no unpublished work.

## Open questions, unresolved by this repository

- **Choosing `k`.** `k=8` was the smallest prefix that cleared every real
  alphabet's own observed maximum in `microbenchmark/run_prefix.jl`'s
  sweep, not a value derived from a formula; a polytope with atoms more
  degenerate than Birkhoff's permutation matrices (which needed `k=8`
  themselves) could need more, and this repository has no rule for picking
  it beyond "sweep `k` and check the collision rate."
- **Scope of "real problems".** Only Birkhoff and the L∞-ball were run.
  Both are now covered by `microbenchmark/run_prefix.jl`'s alphabet-matched
  sweep, so the collision behaviour that made the original full-hash answer
  hold (or not) is measured rather than bounded. A problem whose optimum
  genuinely needs thousands of active vertices (a spectrahedron at high
  dimension, an adversarially-designed cost) still was not tried, and a
  polytope whose vertices are more confusable than a permutation matrix's
  flattened prefix (100% bucket collision at `k=8`, and the prefix hash
  still won) has not been found or ruled out.
- **PFW and BCG were found, not measured.** `references.md` and the draft
  above name four more call sites than the task's own framing listed
  (`pairwise.jl`, `blended_cg.jl`, and two generic drivers); only BPCG was
  run end to end. Worth a second harness pass if this becomes a live PR.
- **README.fr.md** was updated for the new headline and answer, but the
  fuller "Prefix hashing" section that carries the crossover table and the
  collision-rate explanation was not translated; a native read is still
  worth doing before anything with this repository's URL in it goes out.
- **The CI workflow was written, not observed green.** No GitHub run of
  `.github/workflows/ci.yml` has happened yet; it should be watched once
  pushed, per this repository's own no-push rule.
