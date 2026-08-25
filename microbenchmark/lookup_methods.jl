# Lookup strategies for finding an atom in an active set of dense
# `Vector{Float64}` atoms, compared head to head: the linear scan
# `find_atom` uses today, copied from `active_set.jl` (the `Array` branch of
# `_unsafe_equal`, line 499); a `Dict` keyed by the whole atom, which hashes
# every coordinate; and a `Dict` keyed only by the atom's first `k`
# coordinates (a "prefix hash"). A full-atom `Dict` hashes first and falls
# back to exact equality on a collision, so it is sound under the same exact
# semantics `_unsafe_equal` gives an `Array`: see ../README.md's "Why
# hashing is sound here" for why that matters. A prefix hash is sound for
# the identical reason: a bucket hit is never trusted on its own, it is
# always confirmed against the whole atom with the same `!=` the scan uses,
# so shortening the hash can only change speed, never correctness.
module LookupMethods

using SparseArrays

export linear_scan,
    build_dict,
    dict_lookup,
    PrefixIndex,
    build_prefix_index,
    prefix_lookup,
    collision_stats,
    sparse_linear_scan,
    sparse_prefix,
    build_sparse_prefix_index,
    sparse_prefix_lookup

function linear_scan(atoms::Vector{Vector{Float64}}, query::Vector{Float64})
    @inbounds for idx in eachindex(atoms)
        a = atoms[idx]
        found = true
        @inbounds for j in eachindex(a)
            if a[j] != query[j]
                found = false
                break
            end
        end
        if found
            return idx
        end
    end
    return -1
end

function build_dict(atoms::Vector{Vector{Float64}})
    d = Dict{Vector{Float64},Int}()
    for (idx, atom) in enumerate(atoms)
        d[atom] = idx
    end
    return d
end

dict_lookup(d::Dict{Vector{Float64},Int}, query::Vector{Float64}) = get(d, query, -1)

# A prefix-hash index: buckets keyed by an atom's first `k` coordinates,
# each holding the indices of every atom that shares that prefix. Building
# it copies each atom's prefix (`atom[1:k]`), the same O(k) cost a real
# hash over k coordinates would pay; that copy is not timed, by the same
# convention `build_dict` follows (see microbenchmark/run.jl's header).
struct PrefixIndex
    k::Int
    buckets::Dict{Vector{Float64},Vector{Int}}
end

function build_prefix_index(atoms::Vector{Vector{Float64}}, k::Int)
    buckets = Dict{Vector{Float64},Vector{Int}}()
    for (idx, atom) in enumerate(atoms)
        key = atom[1:k]
        if haskey(buckets, key)
            push!(buckets[key], idx)
        else
            buckets[key] = [idx]
        end
    end
    return PrefixIndex(k, buckets)
end

# Hashes only `query`'s first `k` coordinates to find the candidate bucket,
# then confirms every candidate against the *whole* atom with the same
# elementwise `!=` `linear_scan` uses. This confirm step is what makes a
# prefix hash sound rather than a shortcut that trades correctness for
# speed: a bucket match is a candidate, never an answer, exactly as a
# full-atom Dict's own hash collisions are never trusted without it.
function prefix_lookup(index::PrefixIndex, atoms::Vector{Vector{Float64}}, query::Vector{Float64})
    key = query[1:index.k]
    candidates = get(index.buckets, key, nothing)
    candidates === nothing && return -1
    @inbounds for idx in candidates
        a = atoms[idx]
        found = true
        @inbounds for j in eachindex(a)
            if a[j] != query[j]
                found = false
                break
            end
        end
        found && return idx
    end
    return -1
end

# How often a prefix-hash bucket holds more than one atom, i.e. how often a
# lookup that reaches a populated bucket still has more than one full
# comparison to make. Two views on it: `bucket_collision_rate` (the
# fraction of *buckets* that hold >1 atom, what the task asked for
# literally) and `atom_collision_rate` (the fraction of *atoms* that live
# in such a bucket, what actually drives the average confirm cost per
# lookup, since a bigger bucket is hit by more queries than a singleton
# one). A short prefix over a low-entropy alphabet (box corners: 1 bit per
# coordinate) can send both rates near 1 even though the hash itself is
# cheap: that is exactly how a cheap prefix hash turns expensive.
function collision_stats(index::PrefixIndex)
    sizes = length.(values(index.buckets))
    n_atoms = sum(sizes; init=0)
    n_buckets = length(sizes)
    atoms_in_collision = sum(s for s in sizes if s > 1; init=0)
    buckets_with_collision = count(>(1), sizes)
    return (
        k=index.k,
        n_atoms=n_atoms,
        n_buckets=n_buckets,
        mean_bucket_size=n_buckets == 0 ? 0.0 : n_atoms / n_buckets,
        max_bucket_size=isempty(sizes) ? 0 : maximum(sizes),
        atom_collision_rate=n_atoms == 0 ? 0.0 : atoms_in_collision / n_atoms,
        bucket_collision_rate=n_buckets == 0 ? 0.0 : buckets_with_collision / n_buckets,
    )
end

# --- Sparse-atom variants -----------------------------------------------
#
# `FrankWolfe.jl`'s own Birkhoff-polytope atoms are `SparseMatrixCSC`
# permutation matrices, and `_unsafe_equal`'s sparse branch (active_set.jl,
# line 513) is `a == b`, SparseArrays' own comparison, not the dense
# elementwise loop above. Timed directly (see run_prefix.jl), that
# comparison costs about 20ns flat regardless of dimension on a mismatch,
# and scales with the atom's nonzero count (n for an n x n permutation
# matrix), not its dimension (n^2), on a match: running the dense
# `linear_scan` above against a *flattened* permutation matrix would
# overstate the real scan's cost, so the permutation-atom scan baseline
# uses these functions, which call the real sparse `==` FrankWolfe.jl
# itself runs, not a reimplementation of it.
#
# A prefix hash still only needs `k` coordinates, sparse or dense: reading
# them straight off the sparse structure via linear indexing is O(k)
# `getindex` calls, not an O(dimension) densify-then-slice.

function sparse_linear_scan(atoms::Vector{<:AbstractSparseArray}, query::AbstractSparseArray)
    @inbounds for idx in eachindex(atoms)
        atoms[idx] == query && return idx
    end
    return -1
end

sparse_prefix(atom::AbstractSparseArray, k::Int) = Float64[atom[i] for i in 1:k]

function build_sparse_prefix_index(atoms::Vector{<:AbstractSparseArray}, k::Int)
    buckets = Dict{Vector{Float64},Vector{Int}}()
    for (idx, atom) in enumerate(atoms)
        key = sparse_prefix(atom, k)
        if haskey(buckets, key)
            push!(buckets[key], idx)
        else
            buckets[key] = [idx]
        end
    end
    return PrefixIndex(k, buckets)
end

function sparse_prefix_lookup(
    index::PrefixIndex,
    atoms::Vector{<:AbstractSparseArray},
    query::AbstractSparseArray,
)
    key = sparse_prefix(query, index.k)
    candidates = get(index.buckets, key, nothing)
    candidates === nothing && return -1
    @inbounds for idx in candidates
        atoms[idx] == query && return idx
    end
    return -1
end

end # module
