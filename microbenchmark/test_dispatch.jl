# The proposed API shape (DECISIONS.md's "The pattern key's integer
# representations" section) argues a real `ActiveSet` should route
# a sparse atom to the folded UInt64 pattern key and a dense atom to the
# existing value-prefix hash, using Julia's own multiple dispatch to
# choose, exactly the shape `_unsafe_equal` itself already has
# (active_set.jl: one method on `Array`, one on
# `SparseArrays.AbstractSparseArray`). test_atom_generators.jl's
# `route_build`/`route_lookup`/`route_scan` are that routing, written as
# two ordinary methods per generic function rather than an `if
# atom isa SparseMatrixCSC` branch anywhere, and test_equivalence.jl and
# test_lifecycle.jl already call them on every atom they touch. What this
# file checks that those do not: that the dispatch actually lands on the
# structure it is supposed to (`isa` on the built index, not just "the
# answer came out right", since a routing bug that happened to produce
# the same answer through the wrong structure would slip past a
# correctness-only check), for every alphabet, and that a mixed batch of
# lookups (sparse atoms first, then dense, interleaved) still routes each
# one correctly rather than "sticking" to whichever branch ran first.
#
# Run: julia --project=. microbenchmark/test_dispatch.jl

using Test, Random

include(joinpath(@__DIR__, "test_atom_generators.jl"))
using .TestAtomGenerators

const MASTER_SEED = 8 # distinct from run.jl(1)..run_pattern_key_reps.jl(4), test_equivalence.jl(5), test_lifecycle.jl(6)
const K = 4
const SIZE = 30

function report_and_check(cond::Bool, context::String)
    cond || println(stderr, "DISPATCH FAILURE (master_seed=$MASTER_SEED): ", context)
    @test cond
end

@testset "sparse atoms route to the pattern key" begin
    rng = Random.Xoshiro(MASTER_SEED)
    spec = only(s for s in ALPHABETS if s.sparse)
    atoms = atom_pool(spec, SIZE, rng)
    index = route_build(atoms, K)
    @test index isa PatternIndexU64

    for pos in (1, SIZE, rand(rng, 1:SIZE, 3)...)
        query = atoms[pos]
        expected = route_scan(atoms, query)
        got = route_lookup(index, atoms, query)
        report_and_check(got == expected, "birkhoff pos=$pos")
    end
end

@testset "dense atoms route to the value prefix" begin
    rng = Random.Xoshiro(MASTER_SEED + 1)
    for spec in (s for s in ALPHABETS if !s.sparse)
        atoms = atom_pool(spec, SIZE, rng)
        index = route_build(atoms, K)
        @test index isa PrefixIndex

        for pos in (1, SIZE, rand(rng, 1:SIZE, 3)...)
            query = atoms[pos]
            expected = route_scan(atoms, query)
            got = route_lookup(index, atoms, query)
            report_and_check(got == expected, "alphabet=$(spec.name) pos=$pos")
        end
    end
end

# Not just "each alphabet works in isolation": build all three indices
# together and query them interleaved (sparse, dense, dense, sparse, ...),
# so a routing bug that only shows up when Julia's method cache has both
# branches "warm" at once (unlikely for `isa`-based dispatch, but this is
# the one place in the suite that would catch it if it existed) has
# somewhere to show up.
@testset "dispatch does not stick: interleaved queries across alphabets" begin
    rng = Random.Xoshiro(MASTER_SEED + 2)
    pools = Dict(spec.name => atom_pool(spec, SIZE, rng) for spec in ALPHABETS)
    indices = Dict(spec.name => route_build(pools[spec.name], K) for spec in ALPHABETS)
    expected_types = Dict(spec.name => (spec.sparse ? PatternIndexU64 : PrefixIndex) for spec in ALPHABETS)

    order = repeat([spec.name for spec in ALPHABETS], 10)
    Random.shuffle!(rng, order)

    for name in order
        report_and_check(indices[name] isa expected_types[name], "index-type name=$name")
        pos = rand(rng, 1:SIZE)
        query = pools[name][pos]
        expected = route_scan(pools[name], query)
        got = route_lookup(indices[name], pools[name], query)
        report_and_check(got == expected, "interleaved name=$name pos=$pos")
    end
end
