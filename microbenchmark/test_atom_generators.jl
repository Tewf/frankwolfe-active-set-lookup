# One real-atom generator per alphabet, and the "route by atom type"
# dispatch this repository's proposed API shape rests on. Shared by
# `test_equivalence.jl`, `test_lifecycle.jl`, and `test_dispatch.jl` so the
# three files test against literally the same generators and the same
# routing code, rather than three near-identical copies of it (see
# ergonomic-conventions' "zero redundancy"). This file has no @testset of
# its own: it is generators and glue, not an assertion.
#
# Three alphabets, matching run_pattern_key_reps.jl's own scope: Birkhoff
# permutation matrices (sparse, routes to the folded UInt64 pattern key,
# the structure the README's opening paragraph now recommends) and two
# dense alphabets, L-infinity box corners and a generic control, both of
# which route to the existing k-coordinate value-prefix hash, since the
# pattern key needs `SparseMatrixCSC`'s own `rowval` and neither dense
# alphabet has it. Dimensions are kept small (n=25 for Birkhoff, d=64 for
# the two dense alphabets, both well above this suite's k=16 ceiling) so
# property and lifecycle tests, which build many atom pools, stay fast:
# real active-set sizes (158-389) are a maximum this suite exercises, not
# a dimension it needs to reproduce, since dispatch and equivalence do not
# depend on dimension.
module TestAtomGenerators

using FrankWolfe, Random, SparseArrays

include(joinpath(@__DIR__, "lookup_methods.jl"))
include(joinpath(@__DIR__, "sparse_pattern.jl"))
include(joinpath(@__DIR__, "pattern_key_reps.jl"))
using .LookupMethods, .SparsePatternLookup, .PatternKeyReps

export AlphabetSpec, ALPHABETS, atom_pool, route_build, route_lookup, route_scan, route_key,
    PatternIndexU64, PrefixIndex, pattern_key_uint64

const BIRKHOFF_N = 25
const DENSE_DIM = 64

struct AlphabetSpec
    name::Symbol
    atom::Function       # (rng) -> one fresh atom
    sparse::Bool         # true routes to the pattern key, false to the value prefix
    eltype::Type         # concrete atom type, so an empty pool still type-checks
end

birkhoff_atom(rng) = FrankWolfe.compute_extreme_point(FrankWolfe.BirkhoffPolytopeLMO(), randn(rng, BIRKHOFF_N, BIRKHOFF_N))
linf_atom(rng) = FrankWolfe.compute_extreme_point(FrankWolfe.LpNormBallLMO{Inf}(1.0), randn(rng, DENSE_DIM))
generic_atom(rng) = rand(rng, DENSE_DIM)

const ALPHABETS = [
    AlphabetSpec(:birkhoff, birkhoff_atom, true, SparseMatrixCSC{Float64,Int}),
    AlphabetSpec(:linf, linf_atom, false, Vector{Float64}),
    AlphabetSpec(:generic, generic_atom, false, Vector{Float64}),
]

# A concretely-typed pool, built one atom at a time from an explicit rng
# (never the global one), so a whole atom set is reproducible from a
# single seed regardless of what else ran first in the same process.
# `size == 0` still produces a correctly-typed empty vector, unlike a bare
# `[spec.atom(rng) for _ in 1:0]`, which would be `Any[]`.
function atom_pool(spec::AlphabetSpec, size::Int, rng)
    pool = Vector{spec.eltype}(undef, size)
    for i in 1:size
        pool[i] = spec.atom(rng)
    end
    return pool
end

# The dispatch a real implementation would ship: one method per atom
# family, chosen by Julia's own multiple dispatch on the atoms vector's
# element type, the same shape `_unsafe_equal` already has (active_set.jl:
# one method on `Array`, one on `SparseArrays.AbstractSparseArray`).
# test_dispatch.jl checks this routing directly; test_equivalence.jl and
# test_lifecycle.jl call these same three functions rather than a
# hand-rolled `if spec.sparse` branch, so a routing bug would show up in
# every test that touches an atom, not only in test_dispatch.jl.
route_build(atoms::Vector{<:SparseMatrixCSC}, k::Int) = build_pattern_index_u64(atoms, k)
route_build(atoms::Vector{<:Array}, k::Int) = build_prefix_index(atoms, k)

route_lookup(index::PatternIndexU64, atoms::Vector{<:SparseMatrixCSC}, query) = pattern_lookup_u64(index, atoms, query)
route_lookup(index::PrefixIndex, atoms::Vector{<:Array}, query) = prefix_lookup(index, atoms, query)

route_scan(atoms::Vector{<:SparseMatrixCSC}, query) = sparse_linear_scan(atoms, query)
route_scan(atoms::Vector{<:Array}, query) = linear_scan(atoms, query)

# The key a bucket-lifecycle insert (`bucket_lifecycle.jl`'s
# `bucket_insert!`) needs to compute for one atom, matching whichever
# structure `route_build` chose for that atom's family.
route_key(::AlphabetSpec, atom::SparseMatrixCSC, k::Int) = pattern_key_uint64(atom, k)
route_key(::AlphabetSpec, atom::Array, k::Int) = atom[1:k]

end # module
