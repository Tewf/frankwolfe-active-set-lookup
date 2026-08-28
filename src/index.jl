# The index structure itself: a `Dict` from a key (`keys.jl`) to the
# positions, in a caller-owned `atoms::Vector`, of every atom that hashes to
# it. This mirrors `FrankWolfe.jl`'s own active set on purpose: an
# `ActiveSet` already owns its atoms in a plain `Vector` field, and what
# this repository proposes is an auxiliary index kept *alongside* that
# field, not a container that owns the atoms itself. So every function here
# takes `atoms` as an explicit argument rather than storing a copy.
#
# One index type per atom family, chosen by `build_index`'s own dispatch on
# `atoms`'s element type, the same way `atom_key` (keys.jl) dispatches on
# one atom at a time. `lookup_atom`/`push_atom!`/`delete_atom!` below are
# each a single generic method over `AtomIndex`: which key function runs is
# decided by `atom_key`'s own dispatch on the atom passed in, not by a
# branch here, so adding a third atom family later means adding one
# `atom_key` method, not touching this file.
module AtomIndexing

using SparseArrays
using ..AtomKeys, ..AtomConfirm

export AtomIndex, SparsePatternIndex, DenseValueIndex, bucket_health,
    build_index, lookup_atom, push_atom!, delete_atom!

abstract type AtomIndex end

struct SparsePatternIndex <: AtomIndex
    k::Int
    buckets::Dict{UInt64,Vector{Int}}
end

struct DenseValueIndex <: AtomIndex
    k::Int
    buckets::Dict{Vector{Float64},Vector{Int}}
end

# Shared by both index types and by `push_atom!` below: append `position`
# to whichever bucket `key` already has, or start a new one-element bucket.
# Buckets are never pruned back to nothing on the way down (`delete_atom!`
# leaves an empty `Vector{Int}` in place rather than removing the key): a
# stray empty bucket costs a few bytes, and `get(buckets, key, nothing)`
# still answers a lookup correctly either way, so pruning would only add
# work for a memory saving nobody asked for.
function bucket_insert!(buckets::Dict{K,Vector{Int}}, key::K, position::Int) where {K}
    if haskey(buckets, key)
        push!(buckets[key], position)
    else
        buckets[key] = [position]
    end
    return nothing
end

# Build once, over a `Vector` of atoms already in hand. `k` is fixed here,
# for the life of the index: `lookup_atom`/`push_atom!`/`delete_atom!` read
# it off `index.k` rather than taking their own `k`, the same reasoning
# `pattern_key_uint64`'s `bits` field carries in the research harness
# (`microbenchmark/pattern_key_reps.jl`): a query's key has to be folded
# the same way every stored atom's key was, and reading that off the index
# it was built with is the only way to guarantee it rather than trusting
# every caller to repeat the same keyword forever.
function build_index(atoms::AbstractVector{<:SparseMatrixCSC}; k::Int=DEFAULT_K)
    buckets = Dict{UInt64,Vector{Int}}()
    for (idx, atom) in enumerate(atoms)
        bucket_insert!(buckets, atom_key(atom; k=k), idx)
    end
    return SparsePatternIndex(k, buckets)
end

function build_index(atoms::AbstractVector{<:SparseVector}; k::Int=DEFAULT_K)
    buckets = Dict{UInt64,Vector{Int}}()
    for (idx, atom) in enumerate(atoms)
        bucket_insert!(buckets, atom_key(atom; k=k), idx)
    end
    return SparsePatternIndex(k, buckets)
end

function build_index(atoms::AbstractVector{<:Array}; k::Int=DEFAULT_K)
    buckets = Dict{Vector{Float64},Vector{Int}}()
    for (idx, atom) in enumerate(atoms)
        bucket_insert!(buckets, atom_key(atom; k=k), idx)
    end
    return DenseValueIndex(k, buckets)
end


# A one-line health check for the precondition stated in METHOD.md: the k
# positions a key reads have to differ across the atoms, and only sparse
# atoms get that guarantee for free from `nzind`. A dense key reads fixed
# cells and assumes they vary, which box corners satisfy and an atom with a
# dominant repeated value does not.
#
# Mean bucket size is the cheapest way to see it. Near 1.0 means the key is
# separating atoms and every lookup confirms about one candidate. Near the
# atom count means every atom shares a key, so each lookup hashes and then
# scans the whole set, which is slower than the plain scan this replaces
# while still being correct. Raising k will not rescue that case: if the
# cells being read are constant, reading more constant cells adds nothing.
# The fix is to read positions that vary. `REJECTED.md` has the measurements
# on selecting such positions and why that was not worth it where `nzind`
# already supplies them.
function bucket_health(index::AtomIndex)
    isempty(index.buckets) && return 0.0
    total = sum(length(b) for b in values(index.buckets))
    return total / length(index.buckets)
end

# The one lookup, for either index type: hash `query` down to a bucket
# (`atom_key` dispatches on `query`'s own type, independent of which
# concrete `AtomIndex` was passed in), then confirm every candidate against
# the whole atom before trusting it (`confirm_match`, confirm.jl). Returns
# the atom's position in `atoms`, or -1, mirroring `FrankWolfe.jl`'s own
# `find_atom`.
function lookup_atom(index::AtomIndex, atoms::AbstractVector, query)
    key = atom_key(query; k=index.k)
    candidates = get(index.buckets, key, nothing)
    candidates === nothing && return -1
    @inbounds for idx in candidates
        confirm_match(atoms[idx], query) && return idx
    end
    return -1
end

# Append `atom` to `atoms` and record it in the index, mirroring
# `active_set_update!`'s `push!` branch: every real BPCG run this
# repository measured took this branch on every single `find_atom` miss
# (`measurement/results.csv`'s `find_atom_hits` column is 0 for all three).
function push_atom!(index::AtomIndex, atoms::AbstractVector, atom)
    push!(atoms, atom)
    bucket_insert!(index.buckets, atom_key(atom; k=index.k), length(atoms))
    return nothing
end

# Remove the atom at `pos` from `atoms` and repair the index, mirroring
# `active_set_cleanup!`'s `deleteat!`. `deleteat!` shifts every later
# atom's position down by one, and there is no way to find every bucket
# entry that needs shifting short of visiting all of them (the removed
# atom's own bucket is not knowable without its key, and any *other*
# bucket may hold positions past `pos`), so this walk is O(index size)
# regardless of how narrow any one bucket is: the same total work
# `deleteat!` already does shifting the Vector itself, paid again here.
# Real, but rare in the three workloads this repository measured (`deleteat!`
# fired 2, 1, and 0 times across 8,002-20,002 iterations); see METHOD.md.
function delete_atom!(index::AtomIndex, atoms::AbstractVector, pos::Int)
    deleteat!(atoms, pos)
    for bucket in values(index.buckets)
        filter!(!=(pos), bucket)
        @inbounds for i in eachindex(bucket)
            bucket[i] > pos && (bucket[i] -= 1)
        end
    end
    return nothing
end

end # module
