# How the folded sparse-pattern key works, and why it is correct

This file states the idea and the correctness argument. `README.md` has the
number; `REJECTED.md` has what this replaced and why; `src/` is the code.
Longer than 80 lines because both the idea and the correctness argument need
a worked example each to actually land, not just a definition.

## The key idea, in plain terms

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
