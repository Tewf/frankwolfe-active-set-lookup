# Run both algorithms on one small problem and narrate what happened: how
# many steps of each kind, how often the membership question was asked,
# what the three answers were, and whether the certificate's comparison
# held every time the oracle's vertex was appended. Read guide/README.md
# first; this is its numbers.
#
#   julia --project=. guide/run.jl
include("frank_wolfe.jl")
using .GuideFrankWolfe, .GuideFrankWolfe.BirkhoffToy, .GuideFrankWolfe.GuideLookups,
    .GuideFrankWolfe.GuideLookups.GuideActiveSet
using LinearAlgebra, Random

const N = 4
Random.seed!(1)
lmo = BruteForceLMO(N)
# A target strictly inside the polytope, so the optimum needs several
# vertices at once: a random mixture of a few of them.
chosen = lmo.vertices[randperm(length(lmo.vertices))[1:5]]
weights = rand(5); weights ./= sum(weights)
problem = Quadratic(Matrix(sum(w .* v for (w, v) in zip(weights, chosen))))
x0 = extreme_point(lmo, randn(N, N))

println("Birkhoff polytope, n = $N: $(length(lmo.vertices)) vertices, target inside, start at a vertex.")

x, gaps = plain_frank_wolfe(problem, lmo, x0; iterations=2000)
println("\nPlain Frank-Wolfe: $(length(gaps)) iterations, final dual gap $(round(gaps[end], sigdigits=3)), ",
    "objective $(round(objective(problem, x), sigdigits=3)): the 1/t zigzag toward an optimum on a face. ",
    "No active set, no lookup, ever.")

lookup = CrossChecked(ActiveSet(x0))
A, records = blended_pairwise(problem, lmo, x0; iterations=2000, lookup=lookup)
pairwise = count(r -> r.kind == :pairwise, records)
fw = filter(r -> r.kind == :frank_wolfe, records)
println("\nBlended pairwise: $(length(records)) iterations, final dual gap $(round(records[end].dual_gap, sigdigits=3)), ",
    "objective $(round(objective(problem, A.x), sigdigits=3)), $(length(A.atoms)) active atoms at the end.")
println("  $pairwise pairwise steps: weight moved inside the active set, nothing looked up.")
println("  $(length(fw)) Frank-Wolfe steps: the membership question asked $(length(fw)) times.")
println("    scan / index / certificate answered: ",
    isempty(disagreements(lookup)) ? "the same every time" : "DIFFERENTLY: $(disagreements(lookup))")
println("    answers: $(count(r -> r.position == -1, fw)) 'not active, append', $(count(r -> r.position > 0, fw)) 'already active'")
println("    certificate held (<g,v> < <g,s>) on $(count(r -> r.certified, fw)) of $(length(fw)) Frank-Wolfe steps",
    all(r -> r.certified, fw) ? ": every one, as the step rule predicts." : ".")
println("  active set consistent at the end: $(is_consistent(A))")
