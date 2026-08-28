# The smallest polytope that shows everything: the Birkhoff polytope, whose
# vertices are permutation matrices. Nothing here depends on FrankWolfe.jl;
# the oracle is brute force over all n! permutations, which is exactly what
# an LMO means (the vertex minimising a linear function) with no cleverness
# in the way. Fine up to n = 7.
module BirkhoffToy

using LinearAlgebra, SparseArrays

export permutations, vertex, BruteForceLMO, extreme_point, in_birkhoff, Quadratic, objective, gradient, best_step

# Every permutation of 1:n, as vectors, built by inserting n into each slot
# of every permutation of 1:n-1.
function permutations(n::Int)
    n == 0 && return [Int[]]
    out = Vector{Int}[]
    for p in permutations(n - 1), k in 1:n
        push!(out, [p[1:k-1]; n; p[k:end]])
    end
    return out
end

# The vertex for permutation p: a 1 at (i, p[i]) for every row i, stored
# sparse, the way FrankWolfe.jl's own Birkhoff oracle returns it.
vertex(p::Vector{Int}) = sparse(collect(1:length(p)), p, ones(length(p)), length(p), length(p))

# The linear minimisation oracle: all vertices kept in a list, and a query
# is one dot product per vertex. Ties go to the first vertex in the list,
# so the answer is deterministic. The values are computed with the same
# `dot(g, vertex)` the algorithms use, so a score computed here and a score
# computed on an active atom are the same function of the same inputs.
struct BruteForceLMO
    vertices::Vector{SparseMatrixCSC{Float64,Int}}
end
BruteForceLMO(n::Int) = BruteForceLMO([vertex(p) for p in permutations(n)])

function extreme_point(lmo::BruteForceLMO, g::AbstractMatrix)
    best, best_value = 1, dot(g, lmo.vertices[1])
    for (i, v) in enumerate(lmo.vertices)
        value = dot(g, v)
        if value < best_value
            best, best_value = i, value
        end
    end
    return lmo.vertices[best]
end

# Membership in the polytope, for the tests: nonnegative, every row and
# column summing to one.
function in_birkhoff(x::AbstractMatrix; tol=1e-9)
    all(>=(-tol), x) || return false
    n = size(x, 1)
    return all(abs.(sum(x; dims=1) .- 1) .<= tol) && all(abs.(sum(x; dims=2) .- 1) .<= tol)
end

# The objective: half the squared distance to a target matrix. Its gradient
# is x - target, and along any direction d the best step is closed-form,
# which keeps the algorithms below free of a line-search routine.
struct Quadratic
    target::Matrix{Float64}
end
objective(q::Quadratic, x) = 0.5 * sum(abs2, x .- q.target)
gradient(q::Quadratic, x) = x .- q.target
function best_step(q::Quadratic, x, d, step_max)
    denominator = dot(d, d)
    denominator == 0 && return 0.0
    return clamp(dot(q.target .- x, d) / denominator, 0.0, step_max)
end

end # module
