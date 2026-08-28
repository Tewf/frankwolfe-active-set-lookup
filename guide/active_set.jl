# The active set, reduced to what the question needs: the vertices that
# carry weight in the current point, their weights, and the point itself,
# kept equal to the weighted sum at all times. FrankWolfe.jl's `ActiveSet`
# has the same three fields. "Active" means nothing more than "has a
# positive weight here".
#
# The two moves an active-set algorithm makes are below. Neither adds or
# removes an atom: appending the oracle's vertex and dropping an emptied
# entry go through the lookup strategy (`lookups.jl`), because an index kept
# beside the atoms has bookkeeping to do on both, and this file should not
# know which strategy is in use.
module GuideActiveSet

using LinearAlgebra, SparseArrays

export ActiveSet, Atom, is_consistent, frank_wolfe_step!, pairwise_step!

const Atom = SparseMatrixCSC{Float64,Int}

mutable struct ActiveSet
    atoms::Vector{Atom}
    weights::Vector{Float64}
    x::Matrix{Float64}
end
ActiveSet(v::Atom) = ActiveSet([v], [1.0], Matrix(v))

# The invariant every test checks after every step: weights positive and
# summing to one, and x equal to the weighted sum of the atoms.
function is_consistent(A::ActiveSet; tol=1e-9)
    isempty(A.atoms) && return false
    length(A.atoms) == length(A.weights) || return false
    all(>(0), A.weights) || return false
    abs(sum(A.weights) - 1) <= tol || return false
    mixture = sum(w .* a for (w, a) in zip(A.weights, A.atoms))
    return maximum(abs.(A.x .- mixture)) <= tol
end

# Move toward the vertex v by γ: every existing weight shrinks by (1 - γ)
# and x moves. The γ that v itself gains is the caller's to place: on an
# existing entry if v was already active, on a new entry otherwise. That
# placement is the whole membership question, and it is decided elsewhere.
function frank_wolfe_step!(A::ActiveSet, v::Atom, γ::Float64)
    A.weights .*= (1 - γ)
    A.x .= (1 - γ) .* A.x .+ γ .* v
    return nothing
end

# Move weight γ from the away atom (index `away`) to the local FW atom
# (index `toward`), both already active, so nothing is looked up. The step
# is capped by the away atom's weight; reaching the cap empties that entry.
function pairwise_step!(A::ActiveSet, toward::Int, away::Int, γ::Float64)
    A.weights[away] -= γ
    A.weights[toward] += γ
    A.x .+= γ .* (A.atoms[toward] .- A.atoms[away])
    return nothing
end

end # module
