# Decisions awaiting Mohamed

Nothing in this repository was posted anywhere. This file is every place a
human judgement call was made or is still needed, so any of it can be
revisited.

## The draft comment for issue #244

**Posting this is Mohamed's to do — not automated, not implied by anything
else in this repository.** Ready to paste as-is if the numbers still look
right on a re-read.

> I measured this rather than guessing. `find_atom`'s linear scan is called
> whenever `active_set_update!` isn't given an index — which turns out to be
> more call sites than I expected: every "add to active set" step of BPCG
> (`blended_pairwise.jl:374`), every pairwise step of plain PFW
> (`pairwise.jl:242`, since `pfw_step` never supplies one), BCG
> (`blended_cg.jl:358`), and the generic corrective/block-coordinate drivers
> (`corrective_frankwolfe.jl:240`, `block_coordinate_algorithms.jl:421`) —
> but never away-step FW, which always tracks its own index.
>
> On three real BPCG runs (Birkhoff polytope at n=25 and n=60, the L∞-ball
> at d=3,000; code at github.com/Tewf/frankwolfe-active-set-lookup), active
> sets topped out at 158–389 atoms and `find_atom` took 0.02%–0.14% of total
> runtime (`measurement/results.csv`). An isolated microbenchmark (linear
> scan vs. `Dict`, sizes 1–20k, dims 16–8,192) shows why: hashing an atom is
> always O(dimension), so it
> only wins once the scan's own O(dimension) worst case is actually being
> paid — which needs either a large active set (crossover from ~50 atoms at
> dim 16 up to several thousand at dim 8,192, for atoms with independent
> coordinates) or atoms that are unusually easy to confuse (crossover at
> 5–10 atoms when atoms share a long common prefix). `_unsafe_equal` is
> exact, so a `Dict` keyed by the atom would be sound, not just faster where
> it helps.
>
> So: no advantage on what I measured, and I think this can be closed on
> that basis unless someone has a workload where active sets do grow into
> the thousands. If it's worth doing anyway for such a workload — would you
> want a hash-augmented default `ActiveSet`, or a separate subtype the way
> `ActiveSetQuadraticProductCaching` (`active_set_quadratic.jl`) already is?
> `active_set_cleanup!`'s `deleteat!` shifts every index past the one
> removed, so a hash→index map needs upkeep on every drop, not just every
> push — cheap enough to seem worth it only where the scan is genuinely a
> bottleneck, which argues for a subtype over complicating the default.

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
- **Post the numbers, offer the code.** "Happy to share the harness if useful"
  — keeps the repository private and puts the choice on them.

**Recommendation: make it public before posting.** A measurement whose code
cannot be inspected is an assertion, and the whole point of answering a
five-year-old question this way is that the answer can be checked. The
repository contains no personal data and no unpublished work.

## Open questions, unresolved by this repository

- **Scope of "real problems".** Only Birkhoff and the L∞-ball were run. A
  problem whose optimum genuinely needs thousands of active vertices (a
  spectrahedron at high dimension, an adversarially-designed cost) was not
  tried, and might reach the microbenchmark's generic crossover. The
  microbenchmark's `adversarial` scenario is a bound, not an observation —
  no `FrankWolfe.jl` polytope was found whose vertices are that confusable.
- **PFW and BCG were found, not measured.** `references.md` and the draft
  above name four more call sites than the task's own framing listed
  (`pairwise.jl`, `blended_cg.jl`, and two generic drivers); only BPCG was
  run end to end. Worth a second harness pass if this becomes a live PR.
- **README.fr.md** is a full draft, not a stub — but technical French
  benefits from a native read before anything with this repository's URL in
  it goes out. A pass, not a rewrite, is what it needs.
- **The CI workflow was written, not observed green.** No GitHub run of
  `.github/workflows/ci.yml` has happened yet; it should be watched once
  pushed, per this repository's own no-push rule.
