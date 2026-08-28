# Three answers to the one question a Frank-Wolfe step asks its active set:
# "is the vertex the oracle just returned already active, and where?" Each
# strategy answers `position_of` the same way (-1 or an index), and each
# owns the two mutations of the atom list, `append_atom!` and
# `remove_atom!`, because an index beside the atoms has bookkeeping to do on
# every insert and delete while the scan and the certificate have none.
#
# The inputs a strategy may use are what the step has in hand at that
# moment: the active set, the vertex `v`, the gradient `g`, the values
# `dot(g, a)` for every active atom (the step computed them to choose its
# local and away vertices), and the index `s` of the smallest. The scan reads
# the atoms; the index reads a few stored positions of `v`; the certificate
# reads two of the values and no atom at all.
module GuideLookups

using LinearAlgebra, SparseArrays
include(joinpath(@__DIR__, "..", "src", "ActiveSetLookup.jl"))
using .ActiveSetLookup
include(joinpath(@__DIR__, "active_set.jl"))
using .GuideActiveSet

export Lookup, ScanLookup, IndexLookup, CertificateLookup, CrossChecked,
    position_of, append_atom!, remove_atom!, purge_empty!, disagreements

abstract type Lookup end

# FrankWolfe.jl's `find_atom`: compare v with every atom until one matches.
struct ScanLookup <: Lookup end
position_of(::ScanLookup, A::ActiveSet, v, g, values, s) = scan_atoms(A.atoms, v)

# The folded structural key from src/index.jl: a Dict beside the atoms,
# built once, then maintained on every append and removal.
mutable struct IndexLookup <: Lookup
    index::Any
end
IndexLookup(A::ActiveSet) = IndexLookup(build_index(A.atoms))
position_of(l::IndexLookup, A::ActiveSet, v, g, values, s) = lookup_atom(l.index, A.atoms, v)
function append_atom!(l::IndexLookup, A::ActiveSet, v, weight::Float64)
    push_atom!(l.index, A.atoms, v)
    push!(A.weights, weight)
    return nothing
end
function remove_atom!(l::IndexLookup, A::ActiveSet, position::Int)
    delete_atom!(l.index, A.atoms, position)
    deleteat!(A.weights, position)
    return nothing
end

# The certificate from src/certificate.jl: dot(g, v) against the smallest
# dot(g, a). Strictly smaller proves v absent; equal means v ties the best
# atom, which one comparison settles.
struct CertificateLookup <: Lookup end
position_of(::CertificateLookup, A::ActiveSet, v, g, values, s) =
    certified_lookup(A.atoms, v, dot(g, v), s, values[s])

# The strategies with no bookkeeping mutate the active set directly.
function append_atom!(::Lookup, A::ActiveSet, v, weight::Float64)
    push!(A.atoms, v)
    push!(A.weights, weight)
    return nothing
end
function remove_atom!(::Lookup, A::ActiveSet, position::Int)
    deleteat!(A.atoms, position)
    deleteat!(A.weights, position)
    return nothing
end

# Drop every entry whose weight reached zero (a full Frank-Wolfe step empties
# all the others; a capped pairwise step empties the away atom), from the
# back so positions stay valid while deleting. Never empties the set.
function purge_empty!(l::Lookup, A::ActiveSet)
    for i in length(A.weights):-1:1
        A.weights[i] <= 0 && length(A.atoms) > 1 && remove_atom!(l, A, i)
    end
    return nothing
end

# All three at once, for the demonstration and the tests: every answer is
# recorded, any disagreement is kept, and the certificate's answer is the
# one returned. The index is the one strategy with state, so the
# bookkeeping calls go to it.
mutable struct CrossChecked <: Lookup
    index::IndexLookup
    answers::Vector{NTuple{3,Int}}
end
CrossChecked(A::ActiveSet) = CrossChecked(IndexLookup(A), NTuple{3,Int}[])
function position_of(l::CrossChecked, A::ActiveSet, v, g, values, s)
    answer = (
        position_of(ScanLookup(), A, v, g, values, s),
        position_of(l.index, A, v, g, values, s),
        position_of(CertificateLookup(), A, v, g, values, s),
    )
    push!(l.answers, answer)
    return answer[3]
end
append_atom!(l::CrossChecked, A::ActiveSet, v, weight::Float64) = append_atom!(l.index, A, v, weight)
remove_atom!(l::CrossChecked, A::ActiveSet, position::Int) = remove_atom!(l.index, A, position)
disagreements(l::CrossChecked) = [a for a in l.answers if !(a[1] == a[2] == a[3])]

end # module
