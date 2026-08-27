# Every bucket-map index in this repository (`lookup_methods.jl`'s
# `PrefixIndex`, `sparse_pattern.jl`'s `PatternIndex`) is shaped the same
# way underneath: a `Dict{K,Vector{Int}}` from a key to the positions, in
# the shared `atoms` Vector, of every atom that hashes to it. Both need the
# identical answer to two questions a lookup-only benchmark never asks:
#
#   - what does `push!(atoms, atom)` cost the index (an insert)?
#   - what does `deleteat!(atoms, pos)` cost the index (a repair)?
#
# `deleteat!` on the underlying Vector does not just remove one atom: it
# shifts every later atom's position down by one, so any index that has
# cached a raw position has to move with it. There is no way to find every
# affected entry short of visiting every stored position (the removed
# atom's own bucket is not knowable in advance without recomputing its
# key, and every *other* bucket may hold positions past the one removed),
# so a repair costs O(size) regardless of how narrow any one bucket is:
# the same total work `deleteat!` already does shifting the Vector itself,
# just paid again by the index on top. Two views on the lifecycle stage's
# brief make sharing this worthwhile: keeping it in one file means
# `PrefixIndex` and `PatternIndex` cannot drift into two different repair
# costs for what is, underneath, the same bucket-map shape (see
# `ergonomic-conventions`'s "zero redundancy").
module BucketLifecycle

export bucket_insert!, bucket_delete_repair!

function bucket_insert!(buckets::Dict{K,Vector{Int}}, key::K, position::Int) where {K}
    if haskey(buckets, key)
        push!(buckets[key], position)
    else
        buckets[key] = [position]
    end
    return nothing
end

# Two passes in one: remove `pos` from whichever bucket holds it (found by
# scanning, since the key that produced it is not passed in: recomputing
# it would need the atom itself, and by the time this runs the atom is
# already gone from `atoms`), and decrement every remaining stored
# position greater than `pos`. Emptied buckets are left in the `Dict`
# rather than pruned: a stray empty `Vector{Int}` costs a few bytes and
# `get(buckets, key, nothing)` still returns it correctly on the next
# lookup, so pruning would only add work for a memory saving nobody asked
# for here.
function bucket_delete_repair!(buckets::Dict{K,Vector{Int}}, pos::Int) where {K}
    for bucket in values(buckets)
        filter!(!=(pos), bucket)
        @inbounds for i in eachindex(bucket)
            bucket[i] > pos && (bucket[i] -= 1)
        end
    end
    return nothing
end

end # module
