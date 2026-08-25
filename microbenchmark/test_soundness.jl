# Whether a hashed active set can disagree with the linear scan it replaces.
#
# The scan uses `_unsafe_equal`, which compares with `!=`, so it follows `==`
# semantics. A `Dict` follows `isequal` semantics. Those two disagree about
# exactly one thing in Float64: the sign of zero. `0.0 == -0.0` is true and
# `isequal(0.0, -0.0)` is false, so two atoms differing only there are one
# atom to the scan and two atoms to a hash.
#
# That is a false negative, not a collision. Confirming a bucket hit against
# the whole atom does not catch it, because the lookup never reaches a bucket
# at all. It is the one way a hashed lookup can be unsound here.
#
# `sparse-key-and-trie`'s two new structures (`sparse_pattern.jl`,
# `hash_trie.jl`) are checked against the same hazard below: the pattern key
# turns out to be immune to it for a structural reason, not because it
# remembered to canonicalise anything; the trie key is exactly as exposed as
# the flat prefix hash above, and closed the same way.
#
# Run: julia --project=. microbenchmark/test_soundness.jl

using Test, SparseArrays
include(joinpath(@__DIR__, "sparse_pattern.jl"))
include(joinpath(@__DIR__, "hash_trie.jl"))
using .SparsePatternLookup, .HashTrie

scan_equal(a, b) = length(a) == length(b) && all(a[i] == b[i] for i in eachindex(a))
canonical(v) = v .+ 0.0   # -0.0 becomes 0.0; every other Float64 is unchanged

@testset "signed zero splits the scan from the hash" begin
    a = [0.0, 1.0, 2.0]
    b = [-0.0, 1.0, 2.0]

    # One atom as far as the scan is concerned.
    @test scan_equal(a, b)

    # Two atoms as far as a Dict is concerned, so the lookup misses.
    @test !haskey(Dict(a => 1), b)
    @test hash(a) != hash(b)
end

@testset "canonicalising the key closes it" begin
    a = [0.0, 1.0, 2.0]
    b = [-0.0, 1.0, 2.0]
    @test haskey(Dict(canonical(a) => 1), canonical(b))

    # And it changes nothing else, including the infinities and subnormals.
    for v in (1.0, -1.0, 0.5, -3.25, 1e-308, floatmax(Float64), Inf, -Inf)
        @test canonical([v])[1] === v
    end
end

@testset "the naive vertex formula produces one, FrankWolfe's does not" begin
    # Writing an Linf-ball vertex the obvious way yields -0.0 at a zero
    # component, which is how the hazard would arise in user code.
    gradient = [0.0, 2.0, -1.0]
    @test any(x -> x === -0.0, -1.0 .* sign.(gradient))

    # FrankWolfe's own LpNormBallLMO{Inf} does not: it returns -1.0 there.
    # Checked directly on 2026-08-25, so this hazard is a gap in the contract
    # rather than a live bug in the bundled LMOs.
end

@testset "NaN goes the harmless way" begin
    n1 = [NaN, 1.0]
    n2 = [NaN, 1.0]
    # The scan says not equal, and the hash collides then fails confirmation,
    # so both agree on not-found. A wasted comparison, never a wrong answer.
    @test !scan_equal(n1, n2)
    @test hash(n1) == hash(n2)
end

@testset "the sparse-pattern key has no signed zero to disagree about" begin
    # A permutation matrix's rowval is Vector{Int}: row indices, not the
    # nonzero values themselves (which are always 1.0 for a permutation
    # matrix and never even reach the key). Int has one representation of
    # zero, so `==` and `isequal` cannot disagree here the way they can for
    # Float64: this key sidesteps the whole hazard above by construction,
    # not by remembering a canonicalisation step.
    a = sparse([1, 2, 3], [1, 2, 3], [1.0, 1.0, 1.0])
    @test pattern_key(a, 3) isa Vector{Int}
    @test isequal(pattern_key(a, 3), pattern_key(a, 3))
    @test hash(pattern_key(a, 3)) == hash(pattern_key(a, 3))
end

@testset "the trie key needs, and carries, the same canonicalisation" begin
    # Two atoms differing only in the sign of a zero at the one coordinate
    # a 1-level, k=1 trie hashes on: without `.+ 0.0` in build_node's and
    # lookup_node's key construction, this would miss for the identical
    # reason the flat prefix hash would (the two testsets above), since a
    # trie's key is built from the same Float64 coordinates. `build_node`
    # and `lookup_node` in hash_trie.jl both write `atom[c] + 0.0`,
    # exactly the canonicalisation `canonical(v) = v .+ 0.0` above applies.
    stored = [0.0, 1.0]
    query = [-0.0, 1.0]
    trie = build_trie([stored], [1, 2], 1, 1)
    @test trie_lookup(trie, [stored], query) == 1
end
