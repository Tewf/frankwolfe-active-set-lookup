# The two algorithms, written to be read. `plain_frank_wolfe` never keeps
# an active set and never asks the question this repository is about.
# `blended_pairwise` keeps one, and takes exactly the decision
# FrankWolfe.jl's non-lazy BPCG takes each iteration: if the active set can
# make enough progress on its own (local gap at least half the dual gap,
# the library's default `sparsity_control = 2`), move weight inside it;
# otherwise step toward the oracle's vertex, which is the one moment the
# vertex has to be found or appended.
#
# Every iteration is recorded, so a reader (and `test/test_guide.jl`) can
# see which branch ran, what the lookup answered, and whether the
# certificate's comparison held.
module GuideFrankWolfe

using LinearAlgebra, SparseArrays
include(joinpath(@__DIR__, "birkhoff.jl"))
include(joinpath(@__DIR__, "lookups.jl"))
using .BirkhoffToy, .GuideLookups, .GuideLookups.GuideActiveSet

export plain_frank_wolfe, blended_pairwise, StepRecord

# x ← x + γ (v - x), nothing remembered but x. The dual gap <g, x - v> is
# the standard bound on how far f(x) is from the optimum, and the stop test.
function plain_frank_wolfe(problem::Quadratic, lmo, x0; iterations=500, epsilon=1e-8)
    x = Matrix(x0)
    gaps = Float64[]
    for _ in 1:iterations
        g = gradient(problem, x)
        v = extreme_point(lmo, g)
        gap = dot(g, x) - dot(g, v)
        push!(gaps, gap)
        gap <= epsilon && break
        d = v .- x
        x .+= best_step(problem, x, d, 1.0) .* d
    end
    return x, gaps
end

struct StepRecord
    kind::Symbol            # :pairwise, :frank_wolfe, or :converged
    dual_gap::Float64
    local_gap::Float64
    position::Int           # what the lookup answered on a Frank-Wolfe step, 0 otherwise
    certified::Bool         # dot(g, v) < dot(g, s) held on that step
    active_atoms::Int       # size of the active set after the step
end

function blended_pairwise(problem::Quadratic, lmo, x0; iterations=500, epsilon=1e-8, lookup::Lookup=CertificateLookup())
    A = ActiveSet(x0)
    records = StepRecord[]
    for _ in 1:iterations
        g = gradient(problem, A.x)
        values = [dot(g, a) for a in A.atoms]      # what active_set_argminmax computes
        s = argmin(values)                          # the local Frank-Wolfe vertex
        a = argmax(values)                          # the away vertex
        local_gap = values[a] - values[s]
        v = extreme_point(lmo, g)                   # the oracle's vertex
        dual_gap = dot(g, A.x) - dot(g, v)
        if dual_gap <= epsilon
            push!(records, StepRecord(:converged, dual_gap, local_gap, 0, false, length(A.atoms)))
            break
        end
        if local_gap >= dual_gap / 2
            # Progress is available inside the active set: shift weight from
            # the away atom to the local FW atom. Both are indices already.
            γ = best_step(problem, A.x, A.atoms[s] .- A.atoms[a], A.weights[a])
            pairwise_step!(A, s, a, γ)
            purge_empty!(lookup, A)
            push!(records, StepRecord(:pairwise, dual_gap, local_gap, 0, false, length(A.atoms)))
        else
            # Step toward v. Is v already active? This is the question.
            position = position_of(lookup, A, v, g, values, s)
            certified = dot(g, v) < values[s]
            γ = best_step(problem, A.x, v .- A.x, 1.0)
            frank_wolfe_step!(A, v, γ)
            position == -1 ? append_atom!(lookup, A, v, γ) : (A.weights[position] += γ)
            purge_empty!(lookup, A)
            push!(records, StepRecord(:frank_wolfe, dual_gap, local_gap, position, certified, length(A.atoms)))
        end
    end
    return A, records
end

end # module
