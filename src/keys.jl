# One idea, two shapes, because the two atom families put their
# information in different places. A Birkhoff-polytope atom (a permutation
# matrix) has values that are all 1.0 and carry nothing: its information is
# *where* the nonzeros sit. An L-infinity-ball atom (a box corner) has no
# sparse structure to read at all: its information is in the coordinate
# *values* themselves. `atom_key` dispatches on the atom's type to pick the
# right one, the same shape `FrankWolfe.jl`'s own `_unsafe_equal`
# (active_set.jl) already uses to dispatch on `Array` vs
# `AbstractSparseArray`; see `METHOD.md` for the full argument and
# `REJECTED.md` for what this replaced.
#
# `k` is how many leading stored indices (sparse) or coordinates (dense) get
# folded into the key. It is an ordinary `Int` keyword everywhere in this
# module, never a compile-time type parameter: `NTuple{k,Int}` measured
# faster in `microbenchmark/pattern_key_reps.jl`, but needs `k` fixed as
# `Val(k)` to stay allocation-free, which would force every caller of this
# module to thread a `Val` through their own API. `REJECTED.md` has the
# numbers for that trade.
module AtomKeys

using SparseArrays

export DEFAULT_K, atom_key

# Measured, not guessed: k=4 is the smallest prefix `run_pattern_key_reps.jl`
# swept that already reaches 0.0% collision at both real Birkhoff sizes
# (158 and 389 atoms) and gives the fold's best total per-iteration cost at
# n=60 (0.812ns, this repository's headline number); k=2 is very slightly
# faster at n=25 but carries a real 27.9%/11.3% collision rate, and k=8
# never wins the fold's own k-sweep at either size. The full comparison is
# `microbenchmark/results_pattern_key_reps_total.csv` (time) and
# `results_pattern_key_reps_collisions.csv` (collision rate).
const DEFAULT_K = 4


# Both sparse keys fold through here, and the reason it takes `nzval` as well
# as the indices is a hazard that is easy to miss.
#
# A sparse array may hold an *explicitly stored zero*: an entry that occupies
# a slot in `nzind`/`rowval` while its value is 0.0. Sparse arithmetic leaves
# these behind routinely, and nothing drops them unless someone calls
# `dropzeros!`. Two arrays can therefore be numerically equal, so `==` and
# `confirm_match` call them one atom, while their stored index arrays differ.
# Keying on the raw indices would then hand them different keys and the
# lookup would miss an atom that is present, which is a false negative that
# no amount of bucket confirmation can catch: the bucket is never reached.
#
# This is the sparse twin of the signed-zero hazard on the dense key below,
# and it is the more reachable of the two. So the fold skips stored zeros and
# keys only on positions that genuinely hold a nonzero, which is exactly the
# set `==` compares. Cost stays O(k) unless an array carries a run of stored
# zeros at its front, which is rare and which only that array pays for.
function _fold_stored_indices(indices, values, k::Int)
    h = zero(UInt64)
    taken = 0
    @inbounds for i in eachindex(indices)
        taken >= k && break
        iszero(values[i]) && continue
        h = hash(indices[i], h)
        taken += 1
    end
    return h
end

# Sparse atoms (`SparseMatrixCSC` specifically, not any `AbstractSparseArray`):
# fold the first k stored row indices (`rowval`, one Int per stored nonzero,
# in column order) into a UInt64 with an incremental hash. Scope is
# deliberately narrower than `_unsafe_equal`'s own `AbstractSparseArray`
# dispatch: `rowval` is `SparseMatrixCSC`'s own field, and "first k stored
# indices" only means "first k columns' row" when nonzeros are stored in
# column order, which `SparseMatrixCSC` guarantees and, say, a `SparseVector`
# (whose analogous field is `nzind`, not `rowval`) does not share. A future
# atom family stored as a different sparse type needs its own method here,
# not a widened one.
#
# The fold is lossy (two different row-index sequences can land on the same
# UInt64): that is a real hash collision, not a prefix tie, and it is why
# every lookup in `index.jl` confirms a bucket hit against the whole atom
# (`confirm.jl`) before trusting it. An `Int` row index has no sign of zero
# to disagree about, so this key needs no canonicalisation at all, unlike
# the dense key below.
function atom_key(atom::SparseMatrixCSC; k::Int=DEFAULT_K)
    return _fold_stored_indices(atom.rowval, atom.nzval, k)
end

# Sparse *vectors*, which is what several of FrankWolfe.jl's own LMOs return:
# `KSparseLMO` hands back a `SparseVector{Float64,Int64}`, and a
# `SparseMatrixCSC` atom flattened with `vec` becomes one too. Its stored
# indices live in `nzind` rather than `rowval`, and for a vector they are
# plain ascending positions, so "the first k stored indices" is if anything
# more directly meaningful here than in the matrix case.
#
# This is a separate method rather than a widening of the one above to
# `AbstractSparseArray`, exactly as that method's own comment prescribes: the
# two types name their index array differently, so one signature covering
# both would have to branch on the field name at run time and would break on
# the next sparse type that names it a third thing. Everything else applies
# unchanged, including that an `Int` index has no signed zero to canonicalise.
function atom_key(atom::SparseVector; k::Int=DEFAULT_K)
    return _fold_stored_indices(atom.nzind, atom.nzval, k)
end

# Dense atoms: no sparse structure to key on, so key on the first k
# coordinate *values* instead, canonicalised with `+ 0.0`. That addition
# turns `-0.0` into `0.0` and leaves every other Float64 bit-identical
# (infinities and subnormals included); without it, two atoms differing
# only in the sign of a zero would be one atom to `confirm_match` (which
# follows `==`, so `0.0 == -0.0`) and two atoms to this key's `Dict` (which
# follows `isequal`, so `0.0 !isequal -0.0`), and the lookup would miss an
# atom `confirm_match` would call equal. See METHOD.md for where a -0.0
# actually comes from and why it is not hypothetical.
function atom_key(atom::Array{<:Real}; k::Int=DEFAULT_K)
    kk = min(k, length(atom))
    return Float64[atom[i] + 0.0 for i in 1:kk]
end

end # module
