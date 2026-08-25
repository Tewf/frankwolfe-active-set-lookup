# Idea 1 of `sparse-key-and-trie` (Mohamed's brief): key a Birkhoff atom on
# *where* its nonzeros are, not on the first few flattened coordinate
# values `lookup_methods.jl`'s `sparse_prefix` reads.
#
# `results_prefix_collisions.csv` shows why the flattened prefix is weak:
# at k=8, a 389-atom Birkhoff active set (n=60) collapses into 9 buckets,
# mean bucket size 43.2, because a flattened permutation matrix's first 8
# entries are row 0, columns 0-7, which hold the single 1 only when the
# permutation sends one of columns 1-8 to row 1, about a 32% chance, and
# the other 68% of atoms all share the same all-zero bucket. That is a
# representation problem, not a hashing problem: a permutation matrix's
# actual information is which row each column's nonzero sits in, and that
# is already sitting in `SparseMatrixCSC`'s own `rowval` array, one Int
# per stored nonzero, in column order, no different in kind from what
# `sparse_prefix` already reads except *which* field of the sparse
# structure it reads. Hashing `rowval[1:k]` costs the same O(k) `getindex`
# calls `sparse_prefix` costs, no densify, and is close to perfectly
# discriminating for a uniformly random n-permutation: two independent
# permutations agree on where their first k columns map with probability
# 1/(n)_k (the falling factorial n(n-1)...(n-k+1)), already below 1 in 10
# million at n=25, k=8.
#
# The key is a `Vector{Int}`, not `Vector{Float64}`: `_unsafe_equal`'s
# false-negative hazard around signed zero (`../DECISIONS.md`) is a
# Float64 `==`-vs-`isequal` disagreement, and an integer row index has no
# sign-of-zero to disagree about, so this key sidesteps that hazard
# entirely rather than needing the `.+ 0.0` canonicalisation
# `lookup_methods.jl` and `hash_trie.jl` still carry for their Float64
# keys.
module SparsePatternLookup

using SparseArrays

export PatternIndex, build_pattern_index, pattern_lookup, pattern_key, pattern_collision_stats

# Defined only for SparseMatrixCSC (not any AbstractSparseArray): `rowval`
# is that concrete type's own field, not part of the abstract sparse
# interface, and a permutation matrix's "first k stored indices" only
# means "first k columns' row" when nonzeros are stored in column order,
# which SparseMatrixCSC guarantees and a general AbstractSparseArray does
# not.
pattern_key(atom::SparseMatrixCSC, k::Int) = atom.rowval[1:min(k, length(atom.rowval))]

struct PatternIndex
    k::Int
    buckets::Dict{Vector{Int},Vector{Int}}
end

function build_pattern_index(atoms::Vector{<:SparseMatrixCSC}, k::Int)
    buckets = Dict{Vector{Int},Vector{Int}}()
    for (idx, atom) in enumerate(atoms)
        key = pattern_key(atom, k)
        if haskey(buckets, key)
            push!(buckets[key], idx)
        else
            buckets[key] = [idx]
        end
    end
    return PatternIndex(k, buckets)
end

# Same confirm-before-trust shape as `prefix_lookup`/`sparse_prefix_lookup`:
# a bucket hit is a candidate list, checked against the whole atom with the
# real sparse `==` `_unsafe_equal` itself runs, never trusted on its own.
function pattern_lookup(index::PatternIndex, atoms::Vector{<:SparseMatrixCSC}, query::SparseMatrixCSC)
    key = pattern_key(query, index.k)
    candidates = get(index.buckets, key, nothing)
    candidates === nothing && return -1
    @inbounds for idx in candidates
        atoms[idx] == query && return idx
    end
    return -1
end

# Same shape as `lookup_methods.jl`'s `collision_stats`, kept as its own
# four-line body rather than shared: the two operate on structurally
# identical `Dict{K,Vector{Int}}` buckets but different key types
# (`Vector{Int}` here, `Vector{Float64}` there), and Julia would need a
# type parameter threaded through both modules to share one method
# without either module depending on the other, which is more machinery
# than four lines justifies.
function pattern_collision_stats(index::PatternIndex)
    sizes = length.(values(index.buckets))
    n_atoms = sum(sizes; init=0)
    n_buckets = length(sizes)
    atoms_in_collision = sum(s for s in sizes if s > 1; init=0)
    return (
        k=index.k,
        n_atoms=n_atoms,
        n_buckets=n_buckets,
        mean_bucket_size=n_buckets == 0 ? 0.0 : n_atoms / n_buckets,
        max_bucket_size=isempty(sizes) ? 0 : maximum(sizes),
        atom_collision_rate=n_atoms == 0 ? 0.0 : atoms_in_collision / n_atoms,
    )
end

end # module
