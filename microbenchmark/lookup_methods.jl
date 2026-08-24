# Two lookup strategies for finding an atom in an active set of dense
# `Vector{Float64}` atoms, compared head to head: the linear scan
# `find_atom` uses today, copied from `active_set.jl` (the `Array` branch of
# `_unsafe_equal`, line 499), and a `Dict` keyed by the atom itself. A `Dict`
# hashes first and falls back to exact equality on a collision, so it is
# sound under the same exact semantics `_unsafe_equal` gives an `Array` —
# see ../README.md's "Why hashing is sound here" for why that matters.
module LookupMethods

export linear_scan, build_dict, dict_lookup

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

end # module
