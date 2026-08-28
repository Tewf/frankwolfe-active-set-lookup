# Whether the two new key representations `pattern_key_reps.jl` adds stay
# sound. Two different hazards, one per representation:
#
#   - `UInt64`: the fold is lossy by construction, so a real collision (two
#     atoms with genuinely different `rowval[1:k]` patterns landing in the
#     same bucket) has to be survivable, not just theoretically unlikely.
#     `pattern_key_uint64`'s default (`bits=64`) makes a collision among
#     real Birkhoff atoms astronomically rare, which is exactly why it
#     cannot be relied on to *demonstrate* the hazard: this file forces one
#     with `bits=1`, per `pattern_key_reps.jl`'s own comment that `bits`
#     exists for this test, not for production tuning.
#   - `NTuple{K,Int}`: no such hazard exists. A `Dict` compares two tuples
#     elementwise, exactly as it compares two `Vector{Int}`s, so this
#     testset only has to show that stays true, not defend against
#     anything.
#
# The two permutations below were found by a short random search (seed 11,
# n=6, k=2, `bits` 1 and 2), not hand-constructed: search code is not kept
# here, only the two vectors it found, since the property being tested
# (two different patterns, one folded bucket) does not depend on how they
# were found.
#
# Run: julia --project=microbenchmark microbenchmark/test_pattern_key_reps.jl

using Test, SparseArrays

include(joinpath(@__DIR__, "sparse_pattern.jl"))
include(joinpath(@__DIR__, "pattern_key_reps.jl"))
using .SparsePatternLookup, .PatternKeyReps

permutation_matrix(p::Vector{Int}) = sparse(collect(1:length(p)), p, ones(length(p)), length(p), length(p))

@testset "UInt64 fold: a real collision, survived by confirmation" begin
    p1 = [6, 1, 3, 4, 2, 5]
    p2 = [2, 3, 6, 4, 5, 1]
    a1 = permutation_matrix(p1)
    a2 = permutation_matrix(p2)
    k = 2

    # Genuinely different patterns: this is not a prefix tie, it is what
    # the test is forcing to collide.
    @test pattern_key(a1, k) != pattern_key(a2, k)

    # Narrowed to 1 bit, both fold to the same UInt64 key.
    @test pattern_key_uint64(a1, k; bits=1) == pattern_key_uint64(a2, k; bits=1)

    # So the index puts both in one bucket: a list, not a single index,
    # exactly what a folded key requires (`pattern_key_reps.jl`'s header).
    atoms = [a1, a2]
    index = build_pattern_index_u64(atoms, k; bits=1)
    @test length(index.buckets) == 1
    bucket = only(values(index.buckets))
    @test sort(bucket) == [1, 2]

    # And the lookup still answers each atom correctly, because
    # `pattern_lookup_u64` confirms every candidate against the whole atom
    # rather than trusting the bucket. This is the assertion the whole
    # testset exists for: correctness through a real, forced collision,
    # not around one that never happens to occur.
    @test pattern_lookup_u64(index, atoms, a1) == 1
    @test pattern_lookup_u64(index, atoms, a2) == 2

    # A query matching neither atom still misses cleanly, hazard or not.
    a3 = permutation_matrix([1, 2, 3, 4, 5, 6])
    @test pattern_lookup_u64(index, atoms, a3) == -1

    # At the default bits=64, this specific pair's fold does not collide:
    # the hazard is real but not the default's problem at this k, n.
    @test pattern_key_uint64(a1, k) != pattern_key_uint64(a2, k)
end

@testset "NTuple{K,Int}: no fold, so no collision to force" begin
    p1 = [6, 1, 3, 4, 2, 5]
    p2 = [2, 3, 6, 4, 5, 1]
    a1 = permutation_matrix(p1)
    a2 = permutation_matrix(p2)
    k = 2

    @test pattern_key_tuple(a1, Val(k)) != pattern_key_tuple(a2, Val(k))

    atoms = [a1, a2]
    index = build_pattern_index_tuple(atoms, Val(k))
    # Two different patterns, two different tuple keys, two buckets: the
    # same k, the same two atoms that shared one UInt64 bucket above,
    # never share a tuple bucket, at any `bits`-style narrowing, because
    # there is no narrowing step to apply.
    @test length(index.buckets) == 2
    @test pattern_lookup_tuple(index, atoms, a1) == 1
    @test pattern_lookup_tuple(index, atoms, a2) == 2
end

@testset "both new representations still agree with the Vector{Int} baseline" begin
    # Not a collision test: a sanity check that neither new representation
    # silently changed *what* gets matched, only *how* the key is stored,
    # across a slightly larger pool where bucket occupancy is closer to
    # what run_pattern_key_reps.jl actually sweeps.
    perms = [
        [1, 2, 3, 4, 5, 6],
        [2, 1, 3, 4, 5, 6],
        [6, 1, 3, 4, 2, 5],
        [2, 3, 6, 4, 5, 1],
        [3, 4, 5, 6, 1, 2],
    ]
    atoms = [permutation_matrix(p) for p in perms]
    k = 2

    index_vec = build_pattern_index(atoms, k)
    index_u64 = build_pattern_index_u64(atoms, k)
    index_tup = build_pattern_index_tuple(atoms, Val(k))

    for (i, a) in enumerate(atoms)
        @test pattern_lookup(index_vec, atoms, a) == i
        @test pattern_lookup_u64(index_u64, atoms, a) == i
        @test pattern_lookup_tuple(index_tup, atoms, a) == i
    end

    miss = permutation_matrix([4, 5, 6, 1, 2, 3])
    @test pattern_lookup(index_vec, atoms, miss) == pattern_lookup_u64(index_u64, atoms, miss)
    @test pattern_lookup_u64(index_u64, atoms, miss) == pattern_lookup_tuple(index_tup, atoms, miss)
end
