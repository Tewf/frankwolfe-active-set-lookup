# Runs each problem in `problems.jl` under the instrumentation in
# `instrumentation.jl`, and writes one row per run: the problem, its
# dimension, the iteration budget, the maximum and mean active-set size
# seen by the callback, how many times `find_atom` was called, how many of
# those calls were hits (the atom was already present), how many individual
# atoms were removed by `deleteat!` (direct drop-step calls and
# `active_set_cleanup!`'s batch calls both counted, see instrumentation.jl),
# the share of total run time `find_atom`'s accumulated time represents,
# and the certificate tally: of the calls made for the LMO's own vertex, how
# many `<g,v> < min_a <g,a>` decided outright, how many were a tie resolved
# by the best atom alone, how many were a tie with some other atom (the only
# case that needs a search), how many had the LMO vertex scoring above the
# minimum, and how many contradicted the scan (zero, or the certificate is
# unsound).
#
# `results.csv` holds the BPCG runs, as before: `microbenchmark/run_lifecycle.jl`
# and `run_pattern_key_reps.jl` read its per-problem call rates by problem
# name, one row each, so its shape is kept and the certificate columns are
# appended. `results_algorithms.csv` holds the same measurement for pairwise
# Frank-Wolfe (with and without lazification) and blended conditional
# gradients, whose `find_atom` call sites (`pairwise.jl:242`,
# `blended_cg.jl:358`) this repository had only found by reading the source.
#
# find_atom_hits exists because microbenchmark/'s prefix-hash sweep needed a
# real "how often is a real FrankWolfe.jl lookup a hit" number rather than a
# guessed query mix (see ../DECISIONS.md and ../README.md). deleteat_calls
# exists for the same reason, for the lifecycle stage's brief: a lookup
# structure also has to survive `deleteat!` shifting every later index down
# by one, and the real rate at which that happens (not a guess) is what
# `microbenchmark/run_lifecycle.jl` amortises its measured repair cost against.
#
# Every run is done twice: once to compile (small iteration budget,
# discarded), once timed. See ../MEASURING.md for why, and for what else was
# running on the machine while this was measured.
include("instrumentation.jl")
include("problems.jl")
using .LookupInstrumentation
using .Problems
using FrankWolfe, TimerOutputs

const ALGORITHMS = [
    (name="bpcg_lazy", run=(f, g, lmo, x0; kw...) -> FrankWolfe.blended_pairwise_conditional_gradient(f, g, lmo, x0; lazy=true, kw...)),
    (name="pfw", run=(f, g, lmo, x0; kw...) -> FrankWolfe.pairwise_frank_wolfe(f, g, lmo, x0; lazy=false, kw...)),
    (name="pfw_lazy", run=(f, g, lmo, x0; kw...) -> FrankWolfe.pairwise_frank_wolfe(f, g, lmo, x0; lazy=true, kw...)),
    (name="bcg", run=(f, g, lmo, x0; kw...) -> FrankWolfe.blended_conditional_gradient(f, g, lmo, x0; kw...)),
]

function run_once(algorithm, problem; max_iteration, epsilon=1e-9)
    sizes = Int[]
    lmo = RecordingLMO(problem.lmo)
    RECORDER[] = lmo
    reset_instrumentation!()
    @timeit TIMER "total" begin
        algorithm.run(
            problem.f,
            problem.grad!,
            lmo,
            problem.x0;
            max_iteration=max_iteration,
            epsilon=epsilon,
            callback=active_set_sizes!(sizes),
            verbose=false,
        )
    end
    RECORDER[] = nothing
    return sizes
end

function measure(algorithm, problem; max_iteration, warmup_iteration=50)
    run_once(algorithm, problem; max_iteration=warmup_iteration) # compiles, discarded
    sizes = run_once(algorithm, problem; max_iteration=max_iteration)
    total_ns = run_time_of(TIMER["total"])   # the tally's own scans excluded
    cert = certificate_counts()
    return (
        algorithm=algorithm.name,
        problem=problem.label,
        dimension=problem.dimension,
        max_iteration=max_iteration,
        iterations_run=length(sizes),
        max_active_set=maximum(sizes),
        mean_active_set=sum(sizes) / length(sizes),
        find_atom_calls=lookup_calls(),
        find_atom_hits=lookup_hits(),
        deleteat_calls=deletion_calls(),
        total_seconds=total_ns / 1e9,
        lookup_share=lookup_share_of(TIMER["total"]),
        certificate_calls=cert[:calls],
        certified_absent=cert[:certified],
        tie_resolved_by_best=cert[:tie_best],
        tie_needs_search=cert[:tie_other],
        inverted=cert[:inverted],
        unsound=cert[:unsound],
    )
end

const PROBLEMS = [
    (birkhoff_problem(25; seed=1), 8000),
    (birkhoff_problem(60; seed=1), 20000),
    (linf_box_problem(3000; seed=1), 15000),
]

const COLUMNS = "problem,dimension,max_iteration,iterations_run,max_active_set,mean_active_set,find_atom_calls,find_atom_hits,deleteat_calls,total_seconds,lookup_share,certificate_calls,certified_absent,tie_resolved_by_best,tie_needs_search,inverted,unsound"

function row_string(r)
    return "$(r.problem),$(r.dimension),$(r.max_iteration),$(r.iterations_run),$(r.max_active_set),$(round(r.mean_active_set, digits=2)),$(r.find_atom_calls),$(r.find_atom_hits),$(r.deleteat_calls),$(round(r.total_seconds, digits=4)),$(round(r.lookup_share, sigdigits=4)),$(r.certificate_calls),$(r.certified_absent),$(r.tie_resolved_by_best),$(r.tie_needs_search),$(r.inverted),$(r.unsound)"
end

function report(r)
    hit_rate = r.find_atom_calls == 0 ? 0.0 : r.find_atom_hits / r.find_atom_calls
    println(
        "$(r.algorithm) $(r.problem): iters=$(r.iterations_run) max=$(r.max_active_set) mean=$(round(r.mean_active_set,digits=1)) ",
        "calls=$(r.find_atom_calls) hits=$(r.find_atom_hits) ($(round(100*hit_rate,digits=2))%) deletes=$(r.deleteat_calls) ",
        "total=$(round(r.total_seconds,digits=3))s lookup_share=$(round(100*r.lookup_share,digits=3))% ",
        "certificate: calls=$(r.certificate_calls) absent=$(r.certified_absent) tie_best=$(r.tie_resolved_by_best) ",
        "tie_search=$(r.tie_needs_search) inverted=$(r.inverted) unsound=$(r.unsound)",
    )
end

results = [measure(a, problem; max_iteration=k) for a in ALGORITHMS for (problem, k) in PROBLEMS]

open(joinpath(@__DIR__, "results.csv"), "w") do io
    println(io, COLUMNS)
    for r in results
        r.algorithm == "bpcg_lazy" && println(io, row_string(r))
    end
end

open(joinpath(@__DIR__, "results_algorithms.csv"), "w") do io
    println(io, "algorithm,", COLUMNS)
    for r in results
        r.algorithm == "bpcg_lazy" || println(io, r.algorithm, ",", row_string(r))
    end
end

foreach(report, results)
