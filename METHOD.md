# How the lookup works, and why it is correct

This file states two ideas and the correctness argument for each. The first,
the absence certificate, does not search at all and is the answer wherever
the caller is a Frank-Wolfe step; the second, the folded sparse-pattern key,
is the answer for a caller that has nothing but the atom. `README.md` has
the numbers; `REJECTED.md` has what these replaced and why; `src/` is the
code. Longer than 80 lines because each idea and each correctness argument
needs a worked example to actually land, not just a definition.

## The certificate: the lookup that does not search

Every Frank-Wolfe variant that keeps an active set does the same thing at
the top of an iteration: it computes the gradient `g`, then the inner
product `<g, a>` of `g` with every active atom `a`, to find the active atom
that scores lowest (the local Frank-Wolfe vertex `s`) and highest (the away
vertex). `active_set_argminmax` in `FrankWolfe.jl` returns both, with their
values. Only then does it call the LMO for the vertex `v` of the whole
polytope that scores lowest, and it computes `<g, v>` for the dual gap.

At that moment the question `find_atom` is about to answer by scanning,
"is `v` already active?", has already been answered by two numbers in hand.
If `v` were active, `<g, v>` would be one of the values just minimised, so
`<g, v> >= <g, s>`. Therefore:

> **`<g, v> < <g, s>` proves `v` is not in the active set.**

One comparison, no index to build, nothing to insert into, nothing to
repair after `deleteat!`, and no knowledge of the atom's type at all. That
is `certified_absent` in `src/certificate.jl`; `certified_lookup` wraps it
into a drop-in for `find_atom`.

When the comparison fails, `<g, v> >= <g, s>`, and since `v` minimises `<g, .>`
over the whole polytope, exact arithmetic gives `<g, v> = <g, s>`: a tie. A
tie has two cases. Either `v` is `s` itself, which one `confirm_match` on
the best atom settles (this is the pairwise Frank-Wolfe case, where the LMO
routinely returns a vertex that is already active), or `v` ties a distinct
active atom, and only then does anything search. `certified_lookup` hands
that case to a `fallback`, the plain scan by default, an index built with
`index.jl` if the caller keeps one. With a real-valued gradient, a tie
between distinct vertices has probability zero; it is reachable with
degenerate gradients (a constant gradient makes every permutation matrix
score the same), and `test/test_certificate.jl` exercises exactly those.

## Why every measured BPCG call was a miss

`measurement/results.csv` records zero hits in every `find_atom` call the
three original BPCG runs made. That is not a property of the problems. At
the one BPCG call site (`blended_pairwise.jl`, the `active_set_update!`
with `nothing` for the index), the FW step is taken only when the pairwise
step was refused, that is when `local_gap < phi / sparsity_control` or
`local_gap < epsilon`, and the FW step itself requires `dual_gap >= epsilon`
and, lazily, `dual_gap >= phi / sparsity_control`. Here `local_gap = <g, a> -
<g, s>` over the active set and `dual_gap = <g, x> - <g, v>`. If `v` were
active, `<g, v> = <g, s>` (above), and since `x` is a convex combination of
active atoms, `<g, x> <= <g, a>`; so `dual_gap <= local_gap`, which
contradicts refusing the pairwise step while accepting the FW step, for any
`sparsity_control >= 1`. **In exact arithmetic, BPCG never reaches
`find_atom` with an atom that is present.** The scan there is redundant
work by construction, and the certificate is the floating-point-safe way to
skip it: when rounding produces a near-tie, the comparison fails and the
search runs as before. `measurement/instrumentation.jl` evaluates the
certificate at every real call and counts what it decided;
`README.md` reports the tallies for BPCG, pairwise Frank-Wolfe with and
without lazification, and blended conditional gradients.

## What the certificate rests on, in floating point

The argument compares `<g, v>` with values of the same form, never a stored
value with a recomputed one. Its only assumption is that the same `dot` on
equal inputs returns the same Float64: if `v == a` elementwise then
`dot(g, v) == dot(g, a)` bit for bit. That holds for a deterministic
kernel, which is what BLAS and SparseArrays provide; `test/test_certificate.jl`
checks it on this machine's OpenBLAS for dense vectors (fresh copies, both
argument orders, views at odd offsets, since `blended_pairwise.jl` writes
`dot(gradient, v)` while `active_set_argminmax` writes `dot(atom,
direction)`) and on permutation matrices. A NaN or Inf in `g` makes the
comparison false, which is the fall-back direction, and a search runs; the
fingerprint form (`certified_lookup` over every `dot(g, a)`) hands a NaN
fingerprint straight to the scan for the same reason. Signed zeros cannot
split an atom here: `-0.0 * g == 0.0 * g` in every sum, and the final
`confirm_match` follows `==`. The one soundness gap the dense key below has
to close does not exist for the certificate, because the certificate keys
on nothing.

## Where the certificate applies, and where the key still does

The certificate needs a caller that holds `<g, v>` and the active set's
minimum for the same `g`. Every algorithm in `FrankWolfe.jl` that reaches
`find_atom` with no index does: BPCG (`blended_pairwise.jl`, from
`active_set_argminmax`), pairwise Frank-Wolfe in both modes
(`pairwise.jl`, whose `pfw_step` discards the minimum `argminmax` computes
and could keep it), the corrective and block-coordinate drivers built on
the same step, and blended conditional gradients (`blended_cg.jl`, whose
`lp_separation_oracle` scans the active set for its best atom and value
against the same direction before it calls the LMO). BCG has a second,
separate waste the measurement made visible: when that oracle returns an
active atom rather than calling the LMO, the step passes it to
`active_set_update!` with no index and `find_atom` scans for an atom the
oracle had in hand a moment before; 98% of BCG's lookups on Birkhoff were
that, and returning the position alongside the atom removes them. Only a
caller that has nothing but the atom, `active_set_update!` invoked
directly, `weight_from_atom`, or an active set built outside a step, has
no certificate to use. For those, the folded key below is the answer, and
for a tie between distinct atoms it is the natural `fallback`.

## The folded sparse-pattern key, for a caller that has only the atom

A Birkhoff-polytope atom is a permutation matrix. Every stored entry is
1.0: the *values* carry no information at all, because they are all the
same value. What actually identifies the atom is *where* those entries
sit, and that is exactly the permutation. So the key hashes positions
(`SparseMatrixCSC.rowval`, the row each stored column maps to), not values.

An L-infinity-ball atom is a box corner, `{-1, +1}` in every coordinate. It
has no sparse structure at all: every coordinate is stored, and there is no
"where" to speak of separately from the atom itself. So there the key does
the opposite thing, and hashes the first few coordinate *values*.

One idea, `src/keys.jl`'s `atom_key`, dispatched on the atom's type, exactly
the way `FrankWolfe.jl`'s own `_unsafe_equal` (`active_set.jl`) already
dispatches equality checking: one method for a dense `Array`, one for a
sparse atom. The dispatch for the pattern key specifically is narrower than
`_unsafe_equal`'s own (`SparseMatrixCSC`, not every `AbstractSparseArray`):
`rowval` is that concrete type's field, and "first k stored indices" only
means "first k columns' row" when nonzeros are stored in column order,
which `SparseMatrixCSC` guarantees and a general sparse type need not.

## Correctness: a collision costs a comparison, never a wrong answer

Both keys are allowed to be imprecise. The pattern key's fold
(`h = hash(rowval[i], h)`, chained over the first k indices) is lossy: two
different row-index sequences can fold to the same `UInt64`. That is a real
hash collision, not merely a shortened-prefix tie.

Nothing here trusts a bucket on its own. `lookup_atom` (`src/index.jl`)
reads a key's bucket as a list of *candidates*, then checks every one
against the whole atom with exact equality (`confirm_match`, `src/confirm.jl`,
mirroring `_unsafe_equal` itself: `==` for a sparse atom, elementwise `!=`
for a dense one) before returning it. So a collision, real or forced,
costs one extra comparison. It never returns the wrong atom, because the
wrong atom always fails that final check. `microbenchmark/test_pattern_key_reps.jl`
and `test/test_public_api.jl` both force a real fold collision (by
narrowing the fold's bit width) and confirm the lookup still answers
correctly through it.

## The one soundness gap this argument does not cover on its own

Confirming a bucket hit protects against **false positives**: two
different atoms sharing a bucket. It does nothing about **false
negatives**, and the dense value key has exactly one.

`confirm_match`'s dense branch compares with `!=`, so it follows `==`
semantics. A `Dict` key follows `isequal` semantics. The two disagree
about exactly one thing in Float64: the sign of zero. `0.0 == -0.0` is
true; `isequal(0.0, -0.0)` is false. So two atoms differing only in the
sign of a zero coordinate are one atom as far as `confirm_match` is
concerned, and two atoms as far as a raw-value `Dict` key is concerned:
the lookup would miss an atom the scan would have found, and it would
miss it *before* reaching any bucket, so confirmation cannot help.

This is not hypothetical. A naive L-infinity-ball vertex formula,
`-1.0 .* sign.(gradient)`, produces `-0.0` at any zero gradient
component, since `sign(0.0) == 0.0`:

    -1.0 .* sign.([0.0, 2.0, -1.0])  ==  [-0.0, -1.0, 1.0]

(`FrankWolfe.jl`'s own bundled `LpNormBallLMO{Inf}` happens not to produce
this, so it is a gap in the contract rather than a live bug in the
bundled LMOs; a user-supplied LMO, or an atom built by scaling or
negating one elsewhere, could still hit it.)

**The fix is one addition per hashed coordinate.** `atom_key`'s dense
method builds the key from `atom[i] + 0.0`, never from `atom[i]` directly.
Adding zero turns `-0.0` into `0.0` and leaves every other Float64
bit-identical, infinities and subnormals included. The cost is k
additions, independent of the atom's dimension, so it does not change
which structure wins on speed.

NaN goes the harmless way without any special-casing: `NaN == NaN` is
`false`, so `confirm_match` already says not equal; `hash(NaN)` and
`isequal(NaN, NaN)` both treat NaN as equal to itself, so the key still
collides and reaches a bucket, and confirmation then fails, agreeing with
the scan on not-found rather than disagreeing with it.

The sparse pattern key has no equivalent gap to close. Its key is
`Vector{Int}`-typed row indices, not Float64 values, and an integer has
exactly one representation of zero, so `==` and `isequal` cannot disagree
about it. This is a structural consequence of what gets keyed, not a
canonicalisation step the pattern key remembers to take.

`microbenchmark/test_soundness.jl` demonstrates and pins both the hazard
and the fix; `test/test_public_api.jl` checks the shipped `atom_key`
closes it.

## The precondition: the positions the key reads must vary

The method is often described as working "for sparse atoms". That is the wrong
boundary and it will mislead whoever reads it next. What the key actually needs
is narrower and simpler:

> **The k positions the key reads must differ across the atoms in the set.**

Sparsity is not the property, it is a way of getting the property for free.
`nzind` lists only positions that hold something, so reading the first k of them
reads k positions that vary by construction. That is why the sparse key needs no
selection, no tuning and no knowledge of the other atoms.

A dense key has no such guarantee. It reads fixed cells and assumes they vary.
Box corners satisfy that assumption, since every coordinate is independently
plus or minus one. **Atoms with a dominant repeated value do not.** Five dense
atoms that are 7.0 everywhere except one late position all produce the same key
and land in one bucket:

    atoms = [fill(7.0, 40) with a single 3.0 at position 12, 19, 25, 31, 37]
    build_index(atoms)   ->   1 bucket for 5 atoms

A simplex vertex stored densely is the same shape, and it is the permutation
case in a different container.

**It fails safe, not wrong.** Confirmation still returns the correct answer, so
the cost is a wasted hash followed by a scan of the whole bucket: slower than
the plain scan it replaced, never incorrect. That is the right way round for a
degradation, but it is still a degradation.

**How to notice.** `bucket_health(index)` returns the mean bucket size. Near 1.0
means the key is separating atoms. Near the atom count means it is reading
constant positions, and the fix is to read different ones rather than to read
more of them: raising k does not help when the cells it adds are constant too.

Selecting positions that vary, by rarity or by spread, is the general answer and
is measured in `REJECTED.md`. It was rejected for the sparse case because the
selection has to be built from the whole atom set, maintained as that set churns,
and can go stale and silently invalidate every stored key. Those costs are real
whether or not the atoms are sparse; they are simply not worth paying when
`nzind` hands you varying positions for nothing.

## Choosing k

`DEFAULT_K = 4` (`src/keys.jl`) is measured, not guessed: the smallest k
`microbenchmark/run_pattern_key_reps.jl`'s sweep found that already
reaches 0.0% collision at both real Birkhoff sizes and gives the fold's
best total per-iteration cost at the larger one (0.812ns at n=60, this
repository's headline number). It is an ordinary `Int` keyword on every
entry point here, on purpose: `NTuple{k,Int}` keys measured faster still
in the same sweep, but only stay allocation-free when `k` is fixed as a
compile-time `Val(k)`, which would force every caller of this module to
carry a `Val` through their own API just to pick `k`. `REJECTED.md` has
the numbers for that trade and why it was declined.

## The lifecycle, briefly

`push_atom!` appends to the caller's `atoms::Vector` and records one
bucket entry; `delete_atom!` calls `deleteat!` (which shifts every later
position down by one) and then walks every bucket in the index to repair
it, since there is no way to know in advance which buckets hold a shifted
position. That repair is O(index size) regardless of bucket width, real,
and not free (see `MEASURING.md` and `README.md`'s lifecycle table for the
numbers), but it fired 2, 1, and 0 times across the three real BPCG runs
this repository measured; see `DECISIONS.md` for what a workload with a
higher deletion rate would change.

## Reading further

What was tried and refused instead, with numbers: `REJECTED.md`. What each
test in this repository actually protects: `TESTING.md`. Open judgement
calls, including scope questions this method does not resolve: `DECISIONS.md`.
