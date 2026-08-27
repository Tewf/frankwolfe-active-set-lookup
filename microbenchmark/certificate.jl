# The research copy of `src/certificate.jl`'s idea, kept apart from `src/`
# for the reason `what-is-where.md` gives for every other structure here:
# the harness compares candidates against each other and against a scan,
# and must not import the module its numbers are meant to justify.
#
# The certificate: `<g,v> < min_a <g,a>` proves `v` is not in the active set,
# because an active copy of `v` would have been one of the values minimised.
# `best_index`/`best_value` are what `active_set_argminmax` already returns;
# `query_value` is the `dot(gradient, v)` every FW step computes for its dual
# gap. There is no structure to build, insert into, or repair, so the only
# cost is this call.
module CertificateLookup

using SparseArrays

export certified_lookup, certified_lookup_values, scan_fallback

scan_fallback(atoms::Vector{<:AbstractSparseArray}, query) = something(findfirst(a -> a == query, atoms), -1)
function scan_fallback(atoms::Vector{Vector{Float64}}, query::Vector{Float64})
    @inbounds for i in eachindex(atoms)
        equal = true
        for j in eachindex(query)
            if atoms[i][j] != query[j]
                equal = false
                break
            end
        end
        equal && return i
    end
    return -1
end

function certified_lookup(atoms, query, query_value::Float64, best_index::Int, best_value::Float64)
    query_value < best_value && return -1
    @inbounds atoms[best_index] == query && return best_index
    return scan_fallback(atoms, query)
end

# With every `dot(g, a)` kept from the argminmax loop, the gradient is a
# fingerprint: only atoms whose value equals the query's can be the query.
function certified_lookup_values(atoms, query, query_value::Float64, values::Vector{Float64})
    @inbounds for i in eachindex(atoms)
        values[i] == query_value && atoms[i] == query && return i
    end
    return -1
end

end # module
