# The guide's own algorithms, checked for the properties the walkthrough
# claims: the oracle really minimises, plain Frank-Wolfe really converges
# inside the polytope, the active set stays a valid mixture after every
# step, the three lookups agree on every question asked, and in the
# Frank-Wolfe branch of the blended algorithm the oracle's vertex is never
# already active, which is the step-rule argument METHOD.md makes. The
# worked n=3 example from guide/README.md is checked by hand at the end.
#
# Seeded with MASTER_SEED = 11 (10 is test_certificate.jl's, 9 is
# test_public_api.jl's, 5-8 the microbenchmark tests', 1 and 4 the sweeps').
#
# Run: julia --project=. test/test_guide.jl

using Test, Random, LinearAlgebra, SparseArrays

include(joinpath(@__DIR__, "..", "guide", "frank_wolfe.jl"))
using .GuideFrankWolfe, .GuideFrankWolfe.BirkhoffToy, .GuideFrankWolfe.GuideLookups,
    .GuideFrankWolfe.GuideLookups.GuideActiveSet

const MASTER_SEED = 11

# A target inside the polytope: a random mixture of `k` vertices.
function inside_target(rng, lmo, k)
    chosen = lmo.vertices[randperm(rng, length(lmo.vertices))[1:k]]
    w = rand(rng, k)
    w ./= sum(w)
    return Quadratic(Matrix(sum(wi .* v for (wi, v) in zip(w, chosen))))
end

@testset "the oracle returns the vertex with the smallest dot(g, v)" begin
    rng = Random.Xoshiro(MASTER_SEED)
    for n in (3, 4, 5), _ in 1:20
        lmo = BruteForceLMO(n)
        @test length(lmo.vertices) == factorial(n)
        @test all(in_birkhoff, lmo.vertices)
        g = randn(rng, n, n)
        v = extreme_point(lmo, g)
        @test all(dot(g, v) <= dot(g, u) for u in lmo.vertices)
    end
end

# The target is a mixture of three corners, so the optimum lies on a face of
# the polytope, not in its relative interior, and plain Frank-Wolfe zigzags
# toward it at the sublinear rate (Canon and Cullum 1968): f(x_T)·T stays
# constant, and no tolerance like 1e-8 is reached in any budget. What the
# theory does give, and what is checked here, is Jaggi (2013): with L = 1 and
# D² = 2n (two permutation matrices differ in at most 2n entries),
# f(x_T) − f* ≤ 2 L D² / (T + 2) after T steps, and some iterate among the
# first T has a dual gap at most 6.75 L D² / (T + 2). f* = 0, the target
# being a point of the polytope.
@testset "plain Frank-Wolfe meets the 1/t bound and never leaves the polytope" begin
    rng = Random.Xoshiro(MASTER_SEED + 1)
    for n in (3, 4), trial in 1:5
        lmo = BruteForceLMO(n)
        problem = inside_target(rng, lmo, 3)
        x0 = extreme_point(lmo, randn(rng, n, n))
        T = 5000
        x, gaps = plain_frank_wolfe(problem, lmo, x0; iterations=T, epsilon=1e-8)
        @test in_birkhoff(x)
        @test all(gaps[i] >= 0 for i in eachindex(gaps))
        # The dual gap bounds the suboptimality, and the exact line search
        # never increases f, so the last gap recorded still bounds the final x.
        @test objective(problem, x) <= gaps[end]
        @test objective(problem, x) <= 2 * 2n / (length(gaps) + 2)
        smallest_so_far = accumulate(min, gaps)
        @test all(smallest_so_far[t] <= 6.75 * 2n / (t + 2) for t in eachindex(gaps))
    end
end

# The active set is checked after every single step, so blended_pairwise is
# driven one iteration at a time through its own `iterations` budget.
@testset "blended pairwise: a valid mixture after every step, three lookups agree, certificate holds" begin
    rng = Random.Xoshiro(MASTER_SEED + 2)
    for n in (3, 4, 5), trial in 1:4
        lmo = BruteForceLMO(n)
        problem = inside_target(rng, lmo, 4)
        x0 = extreme_point(lmo, randn(rng, n, n))
        lookup = CrossChecked(ActiveSet(x0))
        A, records = blended_pairwise(problem, lmo, x0; iterations=3000, epsilon=1e-9, lookup=lookup)
        ctx = "n=$n trial=$trial"
        @test is_consistent(A)
        @test in_birkhoff(A.x)
        @test records[end].kind == :converged
        @test objective(problem, A.x) <= 1e-6
        @test isempty(disagreements(lookup))
        fw_steps = filter(r -> r.kind == :frank_wolfe, records)
        @test !isempty(fw_steps)
        # The step-rule argument: in the Frank-Wolfe branch the oracle's
        # vertex was never already active, and the certificate saw it.
        @test all(r -> r.position == -1, fw_steps)
        @test all(r -> r.certified, fw_steps)
        # The active set never holds a duplicate atom.
        @test length(unique(A.atoms)) == length(A.atoms)
        # The three strategies also agree when the certificate is the only one driving.
        A2, records2 = blended_pairwise(problem, lmo, x0; iterations=3000, epsilon=1e-9, lookup=CertificateLookup())
        A3, records3 = blended_pairwise(problem, lmo, x0; iterations=3000, epsilon=1e-9, lookup=ScanLookup())
        @test [r.kind for r in records2] == [r.kind for r in records3]
        @test A2.atoms == A3.atoms
    end
end

@testset "a step-by-step check of the mixture, including drops" begin
    rng = Random.Xoshiro(MASTER_SEED + 3)
    lmo = BruteForceLMO(4)
    problem = inside_target(rng, lmo, 6)
    x0 = extreme_point(lmo, randn(rng, 4, 4))
    # Rerun with growing budgets so every prefix of the trajectory is inspected.
    for budget in 1:60
        A, records = blended_pairwise(problem, lmo, x0; iterations=budget, epsilon=1e-12, lookup=IndexLookup(ActiveSet(x0)))
        @test is_consistent(A)
        @test length(A.atoms) == records[end].active_atoms
    end
    # Something was dropped along the way (a capped pairwise step or a full FW step).
    A, records = blended_pairwise(problem, lmo, x0; iterations=3000, epsilon=1e-9, lookup=IndexLookup(ActiveSet(x0)))
    sizes = [r.active_atoms for r in records]
    @test any(sizes[i+1] < sizes[i] for i in 1:length(sizes)-1)
end

@testset "the worked example from guide/README.md" begin
    g = [1.0 5.0 9.0; 4.0 2.0 8.0; 7.0 6.0 3.0]
    lmo = BruteForceLMO(3)
    identity_atom = vertex([1, 2, 3])
    cycle_atom = vertex([2, 3, 1])                      # picks g12 + g23 + g31 = 5 + 8 + 7
    @test dot(g, identity_atom) == 6.0
    @test dot(g, cycle_atom) == 20.0
    @test extreme_point(lmo, g) == identity_atom         # the cheapest assignment is the identity

    A = ActiveSet([identity_atom, cycle_atom], [0.7, 0.3], Matrix(0.7 .* identity_atom .+ 0.3 .* cycle_atom))
    values = [dot(g, a) for a in A.atoms]
    s = argmin(values)
    @test s == 1 && values[s] == 6.0 && values[argmax(values)] - values[s] == 14.0
    v = extreme_point(lmo, g)
    @test dot(g, v) == values[s]                          # a tie: v is s itself
    @test position_of(CertificateLookup(), A, v, g, values, s) == 1

    # Make the identity expensive: the cheapest assignment moves elsewhere,
    # scores below every active atom, and the certificate says "absent".
    g2 = copy(g); g2[1, 1] = 10.0
    v2 = extreme_point(lmo, g2)
    values2 = [dot(g2, a) for a in A.atoms]
    @test dot(g2, v2) == 12.0 && v2 == vertex([2, 1, 3])
    @test dot(g2, v2) < minimum(values2)
    @test position_of(CertificateLookup(), A, v2, g2, values2, argmin(values2)) == -1
    @test position_of(ScanLookup(), A, v2, g2, values2, argmin(values2)) == -1
end
