# The lookup that does not search. Every index in this repository asks
# "which stored atom equals this one?" and answers by reading the atom.
# This file asks a different question the algorithm has already answered:
# "can this atom be in the set at all?", using two numbers the caller holds
# before it ever thinks of looking anything up.
#
# Every Frank-Wolfe variant that keeps an active set computes, each
# iteration, the inner product of the gradient `g` with every active atom
# (`active_set_argminmax` in FrankWolfe.jl), so it knows the active atom `s`
# with the smallest value `<g,s>`. It then asks the LMO for the vertex `v`
# and computes `<g,v>` for the dual gap. If `v` were already active, `<g,v>`
# would be one of the values just minimised over, so `<g,v> >= <g,s>`.
# Therefore `<g,v> < <g,s>` proves `v` is absent, in one comparison, with no
# index to build, no insert bookkeeping and no repair after `deleteat!`, for
# any atom type at all. Only a tie needs a search, and a tie's first
# candidate is `s` itself, since an active copy of `v` scores exactly
# `<g,v>` and `s` is the atom that scores lowest.
#
# Floating point does not weaken the argument, because it never compares
# a stored value to a recomputed one: `<g,v>` and each `<g,a>` are the same
# `dot` on the same `g`, and the same `dot` on equal inputs returns the same
# Float64 (checked on OpenBLAS and SparseArrays, both orientations, in
# `test/test_certificate.jl`). A NaN or Inf in `g` makes the comparison
# false, which is the fall-back direction: a search still runs. The
# precondition is on the caller, and it is exactly what `argminmax` already
# gives: `best_value` is the minimum of `dot(g, a)` over the atoms, and
# `query_value` is `dot(g, query)` by the same `dot`. `METHOD.md` has the
# argument in full and the reason every measured BPCG call was a miss.
module AtomCertificate

include(joinpath(@__DIR__, "confirm.jl"))
using .AtomConfirm

export certified_absent, certified_lookup, scan_atoms

# The certificate on its own, for a caller who only wants the guard.
certified_absent(query_value::Real, best_value::Real) = query_value < best_value

# The fall-back every other method here was compared against: FrankWolfe.jl's
# own `find_atom`, with this module's `confirm_match` in place of
# `_unsafe_equal`. First match wins, so duplicates resolve the way `find_atom`
# resolves them.
function scan_atoms(atoms::AbstractVector, query)
    @inbounds for idx in eachindex(atoms)
        confirm_match(atoms[idx], query) && return idx
    end
    return -1
end

# Position of `query` in `atoms`, or -1, given the active set's best atom.
# `best_index` and `best_value` are `argmin` and `min` of `dot(g, a)` over
# `atoms`; `query_value` is `dot(g, query)`. `fallback(atoms, query)` runs
# only on a tie between distinct atoms, which is where `lookup_atom` over a
# built index (`index.jl`) belongs if the caller keeps one; the plain scan is
# the default because a tie is rare enough that its cost does not register
# (`microbenchmark/results_certificate_*.csv`).
#
# The answer matches `scan_atoms` whenever `atoms` holds no duplicate, which
# is the invariant `find_atom` exists to keep; with duplicates it still
# matches provided the caller's argmin picks the first of equal values, as
# `active_set_argminmax`'s strict `<` does.
function certified_lookup(
    atoms::AbstractVector,
    query,
    query_value::Real,
    best_index::Integer,
    best_value::Real;
    fallback=scan_atoms,
)
    isempty(atoms) && return -1
    certified_absent(query_value, best_value) && return -1
    @inbounds confirm_match(atoms[best_index], query) && return Int(best_index)
    return fallback(atoms, query)
end

# The same idea when the caller kept every `dot(g, a)` rather than only the
# minimum: the gradient is then a fingerprint the algorithm computed for
# free, and only atoms whose fingerprint equals `query_value` can be `query`.
# No argmin needed, no index, and the walk is one Float64 comparison per
# atom, which on sparse atoms is an order of magnitude cheaper than one `==`.
# A NaN fingerprint (a NaN or two opposite infinities in `g`) equals nothing,
# itself included, so it says nothing about any atom and the walk hands over
# to the scan rather than miss a present atom.
function certified_lookup(
    atoms::AbstractVector,
    query,
    query_value::Real,
    values::AbstractVector{<:Real},
)
    isnan(query_value) && return scan_atoms(atoms, query)
    @inbounds for idx in eachindex(atoms)
        values[idx] == query_value && confirm_match(atoms[idx], query) && return idx
    end
    return -1
end

end # module
