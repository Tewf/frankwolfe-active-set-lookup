"""
The problems the harness runs BPCG on. Each returns `(f, grad!, lmo, x0,
label)`. Two are the Birkhoff polytope at two sizes — permutation matrices,
so `find_atom`'s `_unsafe_equal` takes the `SparseArrays.AbstractSparseArray`
branch (`active_set.jl`, line 513) and a comparison costs O(n), the atom's
nonzero count, not O(n²) the way a dense atom's elementwise scan would. The
third is the L∞-ball, whose atoms are plain dense `Array`s (the
`_unsafe_equal(a::Array, b::Array)` branch, line 499) so a comparison really
does cost O(dimension) — the case the microbenchmark isolates.
"""
module Problems

using FrankWolfe, LinearAlgebra, Random, SparseArrays

export birkhoff_problem, linf_box_problem

"""
    birkhoff_problem(n; seed)

Quadratic distance to a random `n×n` target over the Birkhoff polytope,
following the package's own `examples/birkhoff_polytope.jl`.
"""
function birkhoff_problem(n::Int; seed::Int)
    Random.seed!(seed)
    xp = rand(n, n)
    normxp2 = dot(xp, xp)
    f(x) = (normxp2 - 2dot(x, xp) + dot(x, x)) / n^2
    function grad!(storage, x)
        @. storage = 2 * (x - xp) / n^2
        return storage
    end
    lmo = FrankWolfe.BirkhoffPolytopeLMO()
    x0 = FrankWolfe.compute_extreme_point(lmo, randn(n, n))
    return (f=f, grad! =grad!, lmo=lmo, x0=x0, label="birkhoff_n$(n)", dimension=n^2)
end

"""
    linf_box_problem(d; seed, interior_scale)

Quadratic distance to a target strictly inside the box `[-1, 1]^d`, scaled by
`interior_scale` so the optimum is a genuine combination of corners rather
than one corner alone.
"""
function linf_box_problem(d::Int; seed::Int, interior_scale::Float64=0.3)
    Random.seed!(seed)
    xp = interior_scale .* (2 .* rand(d) .- 1)
    f(x) = dot(x - xp, x - xp)
    function grad!(storage, x)
        @. storage = 2 * (x - xp)
        return storage
    end
    lmo = FrankWolfe.LpNormBallLMO{Inf}(1.0)
    x0 = FrankWolfe.compute_extreme_point(lmo, randn(d))
    return (f=f, grad! =grad!, lmo=lmo, x0=x0, label="linf_box_d$(d)", dimension=d)
end

end # module
