# Sweeps active-set size and atom dimension, timing the linear scan against
# the Dict lookup for a query atom absent from the set — the linear scan's
# worst case, and the case a genuinely new vertex takes.
#
# Two atom scenarios, because the first sweep run showed they are not close:
# "generic" atoms have independent random coordinates, so `!=` almost always
# exits after the first one or two entries and the scan barely feels the
# dimension; "adversarial" atoms share one fixed random prefix and differ
# only in their last coordinate, forcing every comparison to run the prefix
# out, which is what a per-atom cost of O(dimension) actually requires.
# Hashing an atom always costs O(dimension) — it must touch every entry — so
# the two scenarios give very different crossovers. Real FrankWolfe.jl atoms
# sit somewhere between them; ../measurement/ is what locates them.
#
# Dict construction is not timed: a real hash-based active set would
# maintain the Dict incrementally as atoms are pushed and dropped, the same
# way `atoms` itself is maintained, not rebuild it on every lookup.
include("lookup_methods.jl")
include("timing.jl")
using .LookupMethods
using .Timing
using Random

Random.seed!(1)

generic_atoms(size, dim) = [rand(dim) for _ in 1:size]

function adversarial_atoms(size, dim)
    prefix = rand(dim - 1)
    return [vcat(prefix, rand()) for _ in 1:size]
end

function miss_query(scenario, atoms, dim)
    return scenario == :generic ? rand(dim) : vcat(atoms[1][1:end-1], rand())
end

sizes = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 6000, 6500, 10000, 20000]
dims = [16, 128, 1024, 8192]

rows = []
for scenario in (:generic, :adversarial), dim in dims, size in sizes
    atoms = scenario == :generic ? generic_atoms(size, dim) : adversarial_atoms(size, dim)
    dict = build_dict(atoms)
    query = miss_query(scenario, atoms, dim)

    scan_ns = time_per_call(linear_scan, atoms, query) * 1e9
    dict_ns = time_per_call(dict_lookup, dict, query) * 1e9

    push!(rows, (scenario=scenario, dim=dim, size=size, scan_ns=scan_ns, dict_ns=dict_ns))
    println(
        "$scenario dim=$dim size=$size  scan=$(round(scan_ns,digits=1))ns  dict=$(round(dict_ns,digits=1))ns",
    )
end

open(joinpath(@__DIR__, "results.csv"), "w") do io
    println(io, "scenario,dim,size,scan_miss_ns,dict_miss_ns")
    for r in rows
        println(io, "$(r.scenario),$(r.dim),$(r.size),$(round(r.scan_ns,digits=1)),$(round(r.dict_ns,digits=1))")
    end
end

println()
for scenario in (:generic, :adversarial), dim in dims
    subset = filter(r -> r.scenario == scenario && r.dim == dim, rows)
    crossing = findfirst(r -> r.dict_ns < r.scan_ns, subset)
    label = crossing === nothing ? "no crossover in [1, 5000]" : "crossover at size $(subset[crossing].size)"
    println("$scenario dim=$dim: $label")
end
