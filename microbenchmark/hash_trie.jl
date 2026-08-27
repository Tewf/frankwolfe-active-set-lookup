# Idea 2 of the lifecycle stage: a hash trie over
# coordinate blocks. `run_prefix.jl`'s prefix hash is a single flat level:
# hash `k` coordinates, and whatever lands in one bucket pays a full
# residual scan. A trie adds levels only where atoms are genuinely
# confusable: hash `k` coordinates to pick a bucket; if that bucket still
# holds more than one atom, hash a *further* block of `k` coordinates,
# picked from whichever atoms share the first bucket, and recurse, up to
# `max_depth` levels. The common case (a bucket of one) pays one level,
# same as the flat prefix hash; only a bucket that keeps colliding pays
# more, and only as much more as it needs.
#
# What "coordinate block" a level uses is a tunable, not a fixed choice,
# and it is the whole point of building this rather than only widening
# `k`: four strategies are compared here (`first_k_order`, `strided_order`,
# `random_order`, `selectivity_order`), one of which (selectivity) looks
# at the atoms actually being indexed rather than only at `dim` and `k`.
# All four reduce to the same shape underneath: compute one coordinate
# *order*, a permutation of `1:dim`, once, at whatever point in the build
# each strategy needs (`first_k`/`strided` from `dim` and `k` alone,
# `random` from a seeded shuffle, `selectivity` from the atom pool), and
# then every level of the trie takes its next unclaimed block of `k`
# coordinates from that same order (`level_coords`). This is deliberate:
# it keeps `build_node`/`trie_insert_node`/`lookup_node` below completely
# unaware of which strategy chose the order they are walking.
#
# Confirmed at the leaf with the same exact equality every other structure
# in this repository uses (`_trie_equal`, mirroring `_unsafe_equal`'s two
# branches), so a trie can misjudge how many levels it needs, never
# whether two atoms are actually equal: the same soundness argument
# `lookup_methods.jl`'s header gives a flat prefix hash. Keys carry the
# same `.+ 0.0` signed-zero canonicalisation `run_prefix.jl` and
# `../DECISIONS.md` require for a Float64 key.
module HashTrie

using Random
using SparseArrays

export HashTrieIndex,
    build_trie,
    trie_lookup,
    trie_insert!,
    trie_delete_repair!,
    trie_stats,
    first_k_order,
    strided_order,
    random_order,
    selectivity_order

# --- Coordinate-selection strategies: each returns a permutation of
# `1:dim`, consumed in sequential k-blocks (`level_coords`) as the trie
# goes deeper. -------------------------------------------------------

first_k_order(dim::Int) = collect(1:dim)

# Spreads a level's k picks roughly evenly across `dim` rather than
# bunching them at one end: coordinate `offset` is followed by
# `offset+stride`, `offset+2*stride`, ... before the next offset starts,
# so `level_coords(order, 1, k)` returns k coordinates ~stride apart
# spanning the whole atom, and level 2 gets the next k in that same
# interleaving (a different phase, not a repeat).
function strided_order(dim::Int, k::Int)
    stride = max(1, cld(dim, k))
    order = Int[]
    sizehint!(order, dim)
    for offset in 1:stride
        c = offset
        while c <= dim
            push!(order, c)
            c += stride
        end
    end
    return order
end

random_order(dim::Int, rng::AbstractRNG) = Random.shuffle(rng, collect(1:dim))

# "What databases do when ordering a composite index", as the brief put it:
# rank coordinates by how many distinct values they take across the
# atoms actually being indexed, most first, computed once over the whole
# pool a build starts from (a real composite index picks one column
# order for the whole index, not a different one per query). Ties are
# routine for a low-entropy alphabet: a permutation matrix's flattened
# entries take exactly 2 values (0 or 1) at *every* coordinate, so raw
# cardinality alone cannot rank them at all (see ../DECISIONS.md). The
# tie-break here, the size of the smaller value-class, descending,
# refines rather than replaces "most distinct values": a coordinate split
# close to half-half separates a bucket better than one split 1-in-25 the
# same way a higher-cardinality column would, even when both nominally
# tie on distinct-value count.
function selectivity_order(atoms, dim::Int)
    n = length(atoms)
    scored = Vector{Tuple{Int,Int,Int}}(undef, dim)
    for c in 1:dim
        counts = Dict{Float64,Int}()
        for atom in atoms
            v = atom[c] + 0.0
            counts[v] = get(counts, v, 0) + 1
        end
        largest_class = isempty(counts) ? 0 : maximum(values(counts))
        balance = n - largest_class # atoms outside the majority value
        scored[c] = (c, length(counts), balance)
    end
    sort!(scored; by=x -> (-x[2], -x[3], x[1]))
    return [c for (c, _, _) in scored]
end

# Level `level` (1-based) claims `order[(level-1)*k+1 : level*k]`, clipped
# to `dim`; an empty result means the order is exhausted before `level`,
# i.e. every coordinate has already been spent on shallower levels of
# this path, so a node here has no further block to split on and must
# stay a leaf regardless of `max_depth`.
function level_coords(order::Vector{Int}, level::Int, k::Int)
    dim = length(order)
    lo = (level - 1) * k + 1
    lo > dim && return Int[]
    hi = min(level * k, dim)
    return order[lo:hi]
end

# --- The trie itself --------------------------------------------------

abstract type TrieNode end

struct TrieLeaf <: TrieNode
    members::Vector{Int}
end

struct TrieInternal <: TrieNode
    coords::Vector{Int}
    children::Dict{Vector{Float64},TrieNode}
end

# `root` is mutable because an insert can turn a `TrieLeaf` into a
# `TrieInternal` (a bucket that used to hold one atom now holds two and
# has to split), and a struct field, unlike a `Dict` entry or a `Vector`
# grown in place, cannot be replaced with a different concrete type
# without the enclosing struct itself being mutable.
mutable struct HashTrieIndex
    order::Vector{Int}
    k::Int
    max_depth::Int
    root::TrieNode
end

_trie_equal(a::AbstractSparseArray, b::AbstractSparseArray) = a == b
function _trie_equal(a::AbstractVector{Float64}, b::AbstractVector{Float64})
    @inbounds for j in eachindex(a)
        a[j] != b[j] && return false
    end
    return true
end

function build_node(atoms, members::Vector{Int}, order::Vector{Int}, k::Int, max_depth::Int, level::Int)
    if length(members) <= 1 || level > max_depth
        return TrieLeaf(members)
    end
    coords = level_coords(order, level, k)
    isempty(coords) && return TrieLeaf(members)
    groups = Dict{Vector{Float64},Vector{Int}}()
    for idx in members
        key = Float64[atoms[idx][c] + 0.0 for c in coords]
        push!(get!(() -> Int[], groups, key), idx)
    end
    children = Dict{Vector{Float64},TrieNode}()
    for (key, group_members) in groups
        children[key] = build_node(atoms, group_members, order, k, max_depth, level + 1)
    end
    return TrieInternal(coords, children)
end

function build_trie(atoms, order::Vector{Int}, k::Int, max_depth::Int)
    root = build_node(atoms, collect(eachindex(atoms)), order, k, max_depth, 1)
    return HashTrieIndex(order, k, max_depth, root)
end

function lookup_node(node::TrieLeaf, atoms, query)
    @inbounds for idx in node.members
        _trie_equal(atoms[idx], query) && return idx
    end
    return -1
end

function lookup_node(node::TrieInternal, atoms, query)
    key = Float64[query[c] + 0.0 for c in node.coords]
    child = get(node.children, key, nothing)
    child === nothing && return -1
    return lookup_node(child, atoms, query)
end

trie_lookup(index::HashTrieIndex, atoms, query) = lookup_node(index.root, atoms, query)

# `position` is `atoms`'s own new length: the caller has already
# `push!`ed the atom before calling this, the same convention
# `bucket_lifecycle.jl`'s `bucket_insert!` follows.
function trie_insert_node(node::TrieLeaf, atoms, position::Int, order, k, max_depth, level)
    push!(node.members, position)
    if length(node.members) <= 1 || level > max_depth
        return node
    end
    coords = level_coords(order, level, k)
    isempty(coords) && return node
    # This leaf just grew past one member: split it exactly as build_node
    # would have if it had seen every one of its members up front.
    return build_node(atoms, node.members, order, k, max_depth, level)
end

function trie_insert_node(node::TrieInternal, atoms, position::Int, order, k, max_depth, level)
    key = Float64[atoms[position][c] + 0.0 for c in node.coords]
    if haskey(node.children, key)
        node.children[key] =
            trie_insert_node(node.children[key], atoms, position, order, k, max_depth, level + 1)
    else
        node.children[key] = TrieLeaf([position])
    end
    return node
end

function trie_insert!(index::HashTrieIndex, atoms, position::Int)
    index.root = trie_insert_node(index.root, atoms, position, index.order, index.k, index.max_depth, 1)
    return nothing
end

# One full walk of every leaf, exactly `bucket_lifecycle.jl`'s
# `bucket_delete_repair!` reasoning: `deleteat!(atoms, pos)` shifts every
# later position down by one, there is no way to find every affected leaf
# except visiting all of them, so this costs O(size) regardless of trie
# depth, the same total work `deleteat!` itself does. Which coordinates a
# node used to get here does not matter for a repair, only that every
# stored position still points at the right atom afterwards, so this
# needs no key computation at all, unlike a lookup or an insert.
function trie_repair_node!(node::TrieLeaf, pos::Int)
    filter!(!=(pos), node.members)
    @inbounds for i in eachindex(node.members)
        node.members[i] > pos && (node.members[i] -= 1)
    end
    return nothing
end

function trie_repair_node!(node::TrieInternal, pos::Int)
    for child in values(node.children)
        trie_repair_node!(child, pos)
    end
    return nothing
end

trie_delete_repair!(index::HashTrieIndex, pos::Int) = trie_repair_node!(index.root, pos)

function collect_leaves!(node::TrieLeaf, depth::Int, out::Vector{Tuple{Int,Int}})
    push!(out, (depth, length(node.members)))
    return nothing
end

function collect_leaves!(node::TrieInternal, depth::Int, out::Vector{Tuple{Int,Int}})
    for child in values(node.children)
        collect_leaves!(child, depth + 1, out)
    end
    return nothing
end

# Depth is counted in trie levels below the root (a bucket resolved at
# level 1 has depth 1), and `mean_depth` weights by how many atoms landed
# at each depth, not by how many leaves did, since a lookup's expected
# cost tracks the former.
function trie_stats(index::HashTrieIndex)
    out = Tuple{Int,Int}[]
    collect_leaves!(index.root, 0, out)
    sizes = [s for (_, s) in out]
    n_atoms = sum(sizes; init=0)
    n_leaves = length(sizes)
    atoms_in_collision = sum(s for s in sizes if s > 1; init=0)
    mean_depth = n_atoms == 0 ? 0.0 : sum(d * s for (d, s) in out; init=0) / n_atoms
    return (
        k=index.k,
        max_depth=index.max_depth,
        n_atoms=n_atoms,
        n_leaves=n_leaves,
        mean_leaf_size=n_leaves == 0 ? 0.0 : n_atoms / n_leaves,
        max_leaf_size=isempty(sizes) ? 0 : maximum(sizes),
        atom_collision_rate=n_atoms == 0 ? 0.0 : atoms_in_collision / n_atoms,
        mean_depth=mean_depth,
        max_depth_reached=isempty(sizes) ? 0 : maximum(d for (d, _) in out),
    )
end

end # module
