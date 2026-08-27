# Costs the absence certificate (`certificate.jl`) against the two lookups
# this repository already measured, at each alphabet's own real active-set
# size, for the four cases a Frank-Wolfe step can put it in:
#
#   miss        the LMO vertex is absent and scores below every active atom,
#               so the certificate decides alone (every BPCG call measured)
#   hit_best    the LMO vertex is present, so it ties the best atom and one
#               comparison decides (the pairwise Frank-Wolfe case)
#   tie_search  the LMO vertex is absent but ties a distinct active atom,
#               the one case that falls back to a search (never observed in
#               a real run; forced here by handing the certificate the best
#               atom's own value)
#   values      the same miss and hit answered from every `dot(g, a)`
#               instead of only the minimum, the fingerprint walk
#
# Baselines are timed in the same session on the same atoms: the folded
# pattern key at k=4 for sparse atoms (`pattern_key_reps.jl`), the k=8
# value-prefix hash for dense atoms (`lookup_methods.jl`), and the scan.
# Their insert and delete-repair costs are not re-measured here; the total
# per-iteration figures take them from the committed
# `results_pattern_key_reps_total.csv` and `results_lifecycle_total.csv`,
# which were measured the same way (`../MEASURING.md`), while the
# certificate has no insert or repair cost at all.
#
# Timing follows `timing.jl` (fastest of five batches over a floor). The
# certificate's inputs are the real ones: a random gradient, the real
# argmin over the atoms, the real `dot(gradient, query)`; a miss query is an
# LMO vertex for that gradient, drawn until it is absent, and a hit query is
# the best atom itself.
include("lookup_methods.jl")
include("pattern_key_reps.jl")
include("certificate.jl")
include("timing.jl")
using .LookupMethods, .PatternKeyReps, .CertificateLookup, .Timing
using FrankWolfe, LinearAlgebra, Random, SparseArrays

Random.seed!(4) # the same seed run_pattern_key_reps.jl uses, so the atoms match

const K_PATTERN = 4  # src/keys.jl's DEFAULT_K
const K_PREFIX = 8   # the dense prefix hash's measured choice

birkhoff_lmo = FrankWolfe.BirkhoffPolytopeLMO()
linf_lmo = FrankWolfe.LpNormBallLMO{Inf}(1.0)

scenarios = [
    (alphabet=:birkhoff_n25, real_size=158, lmo=birkhoff_lmo, direction=() -> randn(25, 25)),
    (alphabet=:birkhoff_n60, real_size=389, lmo=birkhoff_lmo, direction=() -> randn(60, 60)),
    (alphabet=:linf_box_d3000, real_size=241, lmo=linf_lmo, direction=() -> randn(3000)),
]

# What the caller holds at a lookup: every dot(g, a), the first index of the
# minimum (argminmax's strict `<`), and the minimum itself.
function caller_state(atoms, g)
    values = [dot(a, g) for a in atoms]
    best = argmin(values)
    return values, best, values[best]
end

build_baseline(atoms::Vector{<:SparseMatrixCSC}) = build_pattern_index_u64(atoms, K_PATTERN)
build_baseline(atoms::Vector{Vector{Float64}}) = build_prefix_index(atoms, K_PREFIX)
baseline_lookup(index::PatternIndexU64, atoms, q) = pattern_lookup_u64(index, atoms, q)
baseline_lookup(index::PrefixIndex, atoms, q) = prefix_lookup(index, atoms, q)
baseline_name(::Vector{<:SparseMatrixCSC}) = "pattern_k$(K_PATTERN)"
baseline_name(::Vector{Vector{Float64}}) = "prefix_k$(K_PREFIX)"
scan_lookup(atoms::Vector{<:SparseMatrixCSC}, q) = sparse_linear_scan(atoms, q)
scan_lookup(atoms::Vector{Vector{Float64}}, q) = linear_scan(atoms, q)

timing_rows = []
for sc in scenarios
    atoms = [FrankWolfe.compute_extreme_point(sc.lmo, sc.direction()) for _ in 1:sc.real_size]
    g = sc.direction()
    values, best, best_value = caller_state(atoms, g)
    miss = FrankWolfe.compute_extreme_point(sc.lmo, g)
    while scan_lookup(atoms, miss) != -1
        g = sc.direction()
        values, best, best_value = caller_state(atoms, g)
        miss = FrankWolfe.compute_extreme_point(sc.lmo, g)
    end
    miss_value = dot(g, miss)
    @assert miss_value < best_value "the LMO vertex should score below every active atom"
    hit = atoms[best]
    hit_value = dot(g, hit)
    @assert hit_value == best_value
    index = build_baseline(atoms)
    # Every structure answers both queries correctly before anything is timed.
    @assert certified_lookup(atoms, miss, miss_value, best, best_value) == -1
    @assert certified_lookup(atoms, hit, hit_value, best, best_value) == best
    @assert certified_lookup(atoms, miss, best_value, best, best_value) == -1   # the forced tie still answers right
    @assert certified_lookup_values(atoms, miss, miss_value, values) == -1
    @assert certified_lookup_values(atoms, hit, hit_value, values) == best
    @assert baseline_lookup(index, atoms, miss) == -1 && baseline_lookup(index, atoms, hit) == best
    @assert scan_lookup(atoms, miss) == -1 && scan_lookup(atoms, hit) == best
    GC.gc()

    ns(f, args...) = time_per_call(f, args...) * 1e9
    cases = [
        ("certificate", "miss", ns(certified_lookup, atoms, miss, miss_value, best, best_value)),
        ("certificate", "hit_best", ns(certified_lookup, atoms, hit, hit_value, best, best_value)),
        ("certificate", "tie_search", ns(certified_lookup, atoms, miss, best_value, best, best_value)),
        ("certificate_values", "miss", ns(certified_lookup_values, atoms, miss, miss_value, values)),
        ("certificate_values", "hit_best", ns(certified_lookup_values, atoms, hit, hit_value, values)),
        (baseline_name(atoms), "miss", ns(baseline_lookup, index, atoms, miss)),
        (baseline_name(atoms), "hit_best", ns(baseline_lookup, index, atoms, hit)),
        ("scan", "miss", ns(scan_lookup, atoms, miss)),
        ("scan", "hit_best", ns(scan_lookup, atoms, hit)),
    ]
    for (structure, case, t) in cases
        push!(timing_rows, (alphabet=sc.alphabet, size=sc.real_size, structure=structure, case=case, ns=t))
        println("  $(sc.alphabet) size=$(sc.real_size) $structure $case: $(round(t, digits=2)) ns")
    end
end

open(joinpath(@__DIR__, "results_certificate_timing.csv"), "w") do io
    println(io, "alphabet,size,structure,case,lookup_ns")
    for r in timing_rows
        println(io, "$(r.alphabet),$(r.size),$(r.structure),$(r.case),$(round(r.ns, digits=2))")
    end
end

# --- Total per-iteration cost at the real call rates, the arithmetic
# run_lifecycle.jl and run_pattern_key_reps.jl already use. The certificate's
# insert and repair cost is zero by construction; the baselines' come from
# their own committed sweeps, so those rows are a re-statement with this
# session's lookup number swapped in, not a fresh measurement. ---------------

function read_csv_dicts(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    return [Dict(header .=> split(line, ",")) for line in lines[2:end]]
end

measured = read_csv_dicts(joinpath(@__DIR__, "..", "measurement", "results.csv"))
reps_total = read_csv_dicts(joinpath(@__DIR__, "results_pattern_key_reps_total.csv"))
lifecycle_total = read_csv_dicts(joinpath(@__DIR__, "results_lifecycle_total.csv"))

function committed_costs(alphabet, structure)
    if startswith(structure, "pattern")
        row = only(r for r in reps_total if r["alphabet"] == alphabet && r["representation"] == "uint64" && r["k"] == "4")
    elseif startswith(structure, "prefix")
        row = only(r for r in reps_total if r["alphabet"] == alphabet && r["representation"] == "vector_f64" && r["k"] == "8")
    else
        row = only(r for r in lifecycle_total if r["alphabet"] == alphabet && r["structure"] == "scan")
    end
    return parse(Float64, row["insert_ns"]), parse(Float64, row["delete_repair_ns"])
end

total_rows = []
for sc in scenarios
    row = only(r for r in measured if r["problem"] == String(sc.alphabet))
    iterations = parse(Int, row["iterations_run"])
    lookup_rate = parse(Int, row["find_atom_calls"]) / iterations
    insert_rate = (parse(Int, row["find_atom_calls"]) - parse(Int, row["find_atom_hits"])) / iterations
    delete_rate = parse(Int, row["deleteat_calls"]) / iterations
    for structure in unique(r.structure for r in timing_rows if r.alphabet == sc.alphabet)
        lookup_ns = only(r.ns for r in timing_rows if r.alphabet == sc.alphabet && r.structure == structure && r.case == "miss")
        insert_ns, repair_ns = startswith(structure, "certificate") ? (0.0, 0.0) : committed_costs(String(sc.alphabet), structure)
        total = lookup_rate * lookup_ns + insert_rate * insert_ns + delete_rate * repair_ns
        push!(total_rows, (alphabet=sc.alphabet, real_size=sc.real_size, structure=structure, lookup_ns=lookup_ns, insert_ns=insert_ns, delete_repair_ns=repair_ns, total_ns_per_iter=total))
    end
end

open(joinpath(@__DIR__, "results_certificate_total.csv"), "w") do io
    println(io, "alphabet,real_size,structure,lookup_ns,insert_ns,delete_repair_ns,total_ns_per_iter")
    for r in sort(total_rows; by=r -> (r.alphabet, r.total_ns_per_iter))
        println(io, "$(r.alphabet),$(r.real_size),$(r.structure),$(round(r.lookup_ns, digits=2)),$(round(r.insert_ns, digits=2)),$(round(r.delete_repair_ns, digits=2)),$(round(r.total_ns_per_iter, digits=4))")
    end
end

println()
println("Total per-iteration cost (ns) at each alphabet's real max active-set size, miss path, real call rates:")
for r in sort(total_rows; by=r -> (r.alphabet, r.total_ns_per_iter))
    println("  $(r.alphabet) (size=$(r.real_size)) $(r.structure): $(round(r.total_ns_per_iter, digits=4)) ns/iter (lookup $(round(r.lookup_ns, digits=2)), insert $(round(r.insert_ns, digits=2)), repair $(round(r.delete_repair_ns, digits=2)))")
end
