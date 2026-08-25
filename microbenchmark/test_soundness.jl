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
# Run: julia --project=. microbenchmark/test_soundness.jl

using Test

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

@testset "an LMO really can produce a negative zero" begin
    # The shape an Linf-ball vertex takes when a gradient component is zero.
    gradient = [0.0, 2.0, -1.0]
    vertex = -1.0 .* sign.(gradient)
    @test any(x -> x === -0.0, vertex)
end

@testset "NaN goes the harmless way" begin
    n1 = [NaN, 1.0]
    n2 = [NaN, 1.0]
    # The scan says not equal, and the hash collides then fails confirmation,
    # so both agree on not-found. A wasted comparison, never a wrong answer.
    @test !scan_equal(n1, n2)
    @test hash(n1) == hash(n2)
end
