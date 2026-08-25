# `pattern-key-integer-hash`'s brief: `sparse_pattern.jl`'s `pattern_key`
# already won (README.md's "Idea 1"), but its key is a `Vector{Int}`, so
# every lookup and every insert allocates a small array before the `Dict`
# is even touched. This file builds two more representations of the exact
# same pattern (`rowval[1:k]`), chosen to remove that allocation:
#
#   - `pattern_key_uint64`: folds the k row indices into one `UInt64` with
#     an incremental hash (`Base.hash` chained through each index), so
#     nothing is allocated building the key. That folding is lossy: two
#     different k-index sequences can fold to the same `UInt64`, a real
#     hash collision, unlike `pattern_key`'s `Vector{Int}`, which a `Dict`
#     compares element-by-element and never confuses with another vector.
#     Soundness survives this exactly the way `pattern_lookup` already
#     survives an ordinary `Dict` bucket holding more than one atom: the
#     bucket is a `Vector{Int}` of candidate positions, never trusted on
#     its own, and every candidate is confirmed against the whole atom
#     with `==` before being returned. That confirmation is not optional
#     here the way it is merely convenient for `Vector{Int}`: without it,
#     a folded collision would return the wrong atom, not just do
#     unnecessary work. `test_pattern_key_reps.jl` forces a real collision
#     (via the `bits` keyword below) and checks the structure still
#     answers correctly through it.
#
#   - `pattern_key_tuple`: an `NTuple{K,Int}`, Julia's fixed-length tuple.
#     Unlike `UInt64`, this has no collision hazard at all: a `Dict`
#     compares two tuples element-by-element, exactly as it compares two
#     `Vector{Int}`s, so two different k-index sequences can never share a
#     tuple key. What changes is only where the k Ints live: a `Vector{Int}`
#     is a heap object with its own header; an `NTuple{K,Int}` is stored
#     inline, in the `Dict`'s own key array, with no separate allocation,
#     *provided* `K` is known to the compiler as a type parameter rather
#     than carried as a runtime `Int`. `pattern_key_tuple` takes `Val(K)`
#     for exactly that reason: `ntuple(f, k::Int)` for a runtime `k` can
#     still avoid allocating for small `k` (Base unrolls it up to length
#     10), but dispatching on `Val(k)` freshly built from a runtime `k`
#     inside a hot loop does allocate (confirmed empirically before this
#     file was written: a fresh `Val(k)` per call cost 32-80 bytes here,
#     the same run showed 0 bytes once `Val(k)` was built once outside the
#     timed call and reused). `run_pattern_key_reps.jl` builds `Val(k)`
#     once per `k` in its sweep and threads that single object through.
#
# Both keep the same scope restriction `pattern_key` already has: defined
# only for `SparseMatrixCSC`, reading `rowval` directly, no densifying.
module PatternKeyReps

using SparseArrays

export pattern_key_uint64,
    PatternIndexU64,
    build_pattern_index_u64,
    pattern_lookup_u64,
    pattern_key_tuple,
    PatternIndexTuple,
    build_pattern_index_tuple,
    pattern_lookup_tuple,
    bucket_collision_stats

# --- UInt64: incremental hash, no allocation, real collisions possible --

# `bits` exists to support `test_pattern_key_reps.jl`'s collision test, not
# production tuning: at the default `bits=64`, the fold uses the whole
# `UInt64` range Julia's own `hash` gives, and a collision among real
# Birkhoff patterns (n up to a few hundred, k up to a few dozen) is not
# something this repository found or looked for beyond that default.
# Masking down to fewer bits forces the collision rate up on demand, which
# is what the correctness test needs: a real collision, not a hoped-for
# one.
function pattern_key_uint64(atom::SparseMatrixCSC, k::Int; bits::Int=64)
    h = zero(UInt64)
    @inbounds for i in 1:min(k, length(atom.rowval))
        h = hash(atom.rowval[i], h)
    end
    bits >= 64 && return h
    mask = (one(UInt64) << bits) - one(UInt64)
    return h & mask
end

# `bits` is carried on the index, not re-passed at lookup time: the fold a
# query's key uses has to match the fold every stored atom's key used, and
# reading it off the index it was built with is the only way to guarantee
# that rather than trusting the caller to repeat the same keyword.
struct PatternIndexU64
    k::Int
    bits::Int
    buckets::Dict{UInt64,Vector{Int}}
end

function build_pattern_index_u64(atoms::Vector{<:SparseMatrixCSC}, k::Int; bits::Int=64)
    buckets = Dict{UInt64,Vector{Int}}()
    for (idx, atom) in enumerate(atoms)
        key = pattern_key_uint64(atom, k; bits=bits)
        if haskey(buckets, key)
            push!(buckets[key], idx)
        else
            buckets[key] = [idx]
        end
    end
    return PatternIndexU64(k, bits, buckets)
end

# Same confirm-before-trust shape as `pattern_lookup`: a bucket hit is a
# candidate list, never an answer on its own. For this key specifically,
# that confirmation is what makes a folded collision harmless rather than
# a wrong answer: two atoms with different real patterns can share a
# `UInt64` bucket, and only the `==` below tells them apart.
function pattern_lookup_u64(index::PatternIndexU64, atoms::Vector{<:SparseMatrixCSC}, query::SparseMatrixCSC)
    key = pattern_key_uint64(query, index.k; bits=index.bits)
    candidates = get(index.buckets, key, nothing)
    candidates === nothing && return -1
    @inbounds for idx in candidates
        atoms[idx] == query && return idx
    end
    return -1
end

# --- NTuple{K,Int}: no allocation, no collision hazard, K fixed at the
# type level -------------------------------------------------------------

# Requires `atom` to have at least `K` stored nonzeros, unlike
# `pattern_key`'s `Vector{Int}`, which clips to `min(k, length(rowval))`
# and so degrades gracefully on a sparser atom. A fixed-length tuple can't
# do that: its length is part of its type, decided by the caller's `Val`,
# not by what `atom` actually holds. Every atom `run_pattern_key_reps.jl`
# builds this against (Birkhoff permutation matrices, n >= 25) has exactly
# n stored nonzeros, always >= the k values swept (2, 4, 8), so this is a
# real restriction on the representation, not a bug, and it does not bite
# at the sizes measured here.
pattern_key_tuple(atom::SparseMatrixCSC, ::Val{K}) where {K} = ntuple(i -> atom.rowval[i], Val(K))

struct PatternIndexTuple{K}
    buckets::Dict{NTuple{K,Int},Vector{Int}}
end

function build_pattern_index_tuple(atoms::Vector{<:SparseMatrixCSC}, ::Val{K}) where {K}
    buckets = Dict{NTuple{K,Int},Vector{Int}}()
    for (idx, atom) in enumerate(atoms)
        key = pattern_key_tuple(atom, Val(K))
        if haskey(buckets, key)
            push!(buckets[key], idx)
        else
            buckets[key] = [idx]
        end
    end
    return PatternIndexTuple{K}(buckets)
end

# Still confirms against the whole atom, for the same reason `pattern_key`
# does: two atoms can share the first K row indices and differ further out
# past position K (a real tie, not a collision). What this key can never
# do, unlike the UInt64 fold above, is disagree with itself:
# `NTuple{K,Int}` equality in a `Dict` is exact, elementwise `Int`
# comparison, so two atoms only ever land in the same bucket here because
# their first K row indices really are identical, never because two
# different index sequences hashed alike.
function pattern_lookup_tuple(index::PatternIndexTuple{K}, atoms::Vector{<:SparseMatrixCSC}, query::SparseMatrixCSC) where {K}
    key = pattern_key_tuple(query, Val(K))
    candidates = get(index.buckets, key, nothing)
    candidates === nothing && return -1
    @inbounds for idx in candidates
        atoms[idx] == query && return idx
    end
    return -1
end

# One collision-stats function for both new representations, taking the
# raw `buckets` Dict rather than either index struct: the logic (bucket
# sizes in, atom/bucket counts and rates out) does not depend on the key
# type at all, unlike `pattern_lookup_u64`/`pattern_lookup_tuple` above,
# which do. `sparse_pattern.jl` and `lookup_methods.jl` each keep their own
# near-identical copy of this rather than sharing one across files (see
# `sparse_pattern.jl`'s own comment on that choice); this file has two key
# types that both want it, so writing it once here, generic over the Dict,
# is the same "more machinery than four lines justifies" judgement landing
# the other way once there is more than one caller in the same file.
function bucket_collision_stats(buckets::Dict)
    sizes = length.(values(buckets))
    n_atoms = sum(sizes; init=0)
    n_buckets = length(sizes)
    atoms_in_collision = sum(s for s in sizes if s > 1; init=0)
    return (
        n_atoms=n_atoms,
        n_buckets=n_buckets,
        mean_bucket_size=n_buckets == 0 ? 0.0 : n_atoms / n_buckets,
        max_bucket_size=isempty(sizes) ? 0 : maximum(sizes),
        atom_collision_rate=n_atoms == 0 ? 0.0 : atoms_in_collision / n_atoms,
    )
end

end # module
