# The one number README.md derived by arithmetic, measured end to end: the
# same pairwise Frank-Wolfe run on stock FrankWolfe.jl and on the branch that
# decides membership from the step's own minimum (`find_atom` with the
# minimum `active_set_argminmax` already computed), same machine, back to
# back. Two versions of one package cannot share a session, so this is two
# processes, one per variant:
#
#   julia --project=measurement measurement/run_end_to_end.jl
#   FRANKWOLFE_CHECKOUT=/path/to/FrankWolfe.jl FRANKWOLFE_VARIANT=master  julia --project=measurement measurement/run_end_to_end.jl
#   FRANKWOLFE_CHECKOUT=/path/to/FrankWolfe.jl FRANKWOLFE_VARIANT=patched julia --project=measurement measurement/run_end_to_end.jl
#
# The first form is the registry release in measurement/Project.toml. The
# others develop a checkout into a temporary environment and leave the
# committed Project.toml untouched; `master` is the commit the branch was
# cut from, and it is the baseline that matters, because the registry
# release can lag master (0.6.4 predates "update dual step update in BPCG
# (#647)", which changes lazy BPCG's iterates, so registry and master differ
# on that run for reasons that have nothing to do with the lookup). Each
# invocation replaces its own variant's rows in results_end_to_end.csv and
# keeps the others', so the file holds all three after three runs. Fastest
# of three timed runs after a warm-up, per MEASURING.md. The counts and the
# final primal value are what must agree between master and patched: the
# certificate changes no arithmetic, so the iterates, the active set and the
# objective are the same to the last bit, and only the time and the number
# of scans may differ.
const PATCHED = get(ENV, "FRANKWOLFE_CHECKOUT", "")
const VARIANT = isempty(PATCHED) ? "stock" : get(ENV, "FRANKWOLFE_VARIANT", "patched")
if !isempty(PATCHED)
    using Pkg
    Pkg.activate(; temp=true, io=devnull)
    Pkg.develop(path=PATCHED; io=devnull)
    Pkg.add("TimerOutputs"; io=devnull)
end

include("instrumentation.jl")
include("problems.jl")
using .LookupInstrumentation
using .Problems
using FrankWolfe, TimerOutputs, LinearAlgebra

const SOURCE = if isempty(PATCHED)
    "registry v$(pkgversion(FrankWolfe))"
else
    branch = readchomp(`git -C $PATCHED rev-parse --abbrev-ref HEAD`)
    sha = readchomp(`git -C $PATCHED rev-parse --short HEAD`)
    (branch == "HEAD" ? sha : "$branch@$sha") * " v$(pkgversion(FrankWolfe))"
end

const RUNS = [
    (algorithm="pfw", problem=birkhoff_problem(25; seed=1), max_iteration=8000,
        run=(f, g, lmo, x0; kw...) -> FrankWolfe.pairwise_frank_wolfe(f, g, lmo, x0; lazy=false, kw...)),
    (algorithm="pfw", problem=birkhoff_problem(60; seed=1), max_iteration=20000,
        run=(f, g, lmo, x0; kw...) -> FrankWolfe.pairwise_frank_wolfe(f, g, lmo, x0; lazy=false, kw...)),
    (algorithm="bpcg_lazy", problem=birkhoff_problem(60; seed=1), max_iteration=20000,
        run=(f, g, lmo, x0; kw...) -> FrankWolfe.blended_pairwise_conditional_gradient(f, g, lmo, x0; lazy=true, kw...)),
]
const TIMED_RUNS = 3

function run_once(r; max_iteration)
    sizes = Int[]
    # no RecordingLMO here: with the recorder armed, instrumentation.jl's
    # certificate tally scans the active set on every find_atom call, and
    # that is not a cost the stock package pays
    reset_instrumentation!()
    result = @timeit TIMER "total" begin
        r.run(r.problem.f, r.problem.grad!, r.problem.lmo, r.problem.x0;
            max_iteration=max_iteration, epsilon=1e-9,
            callback=active_set_sizes!(sizes), verbose=false)
    end
    total = TimerOutputs.time(TIMER["total"]) / 1e9
    return (
        seconds=total,
        lookup_seconds=lookup_share_of(TIMER["total"]) * total,
        iterations=length(sizes),
        calls=lookup_calls(),
        hits=lookup_hits(),
        active=length(result.active_set),
        primal=result.primal,
        dual_gap=result.dual_gap,
    )
end

function measure(r)
    run_once(r; max_iteration=50)                                # compiles, discarded
    timed = [run_once(r; max_iteration=r.max_iteration) for _ in 1:TIMED_RUNS]
    best = timed[argmin(t.seconds for t in timed)]
    # counts and values do not depend on the run; say so if they ever do
    for t in timed
        (t.iterations, t.calls, t.hits, t.active, t.primal) == (best.iterations, best.calls, best.hits, best.active, best.primal) ||
            error("non-deterministic run for $(r.algorithm) $(r.problem.label): $t vs $best")
    end
    return (
        variant=VARIANT, source=SOURCE, algorithm=r.algorithm, problem=r.problem.label,
        max_iteration=r.max_iteration, iterations=best.iterations, runs=TIMED_RUNS,
        seconds_min=best.seconds, seconds_all=join(round.([t.seconds for t in timed]; digits=3), ";"),
        lookup_seconds=best.lookup_seconds, calls=best.calls, hits=best.hits,
        active=best.active, primal=best.primal, dual_gap=best.dual_gap,
    )
end

const COLUMNS = "variant,source,algorithm,problem,max_iteration,iterations_run,timed_runs,seconds_min,seconds_all,find_atom_seconds,find_atom_calls,find_atom_hits,final_active_set,final_primal,final_dual_gap,julia,blas_threads"
row_string(r) = join([
    r.variant, r.source, r.algorithm, r.problem, r.max_iteration, r.iterations, r.runs,
    round(r.seconds_min; digits=3), r.seconds_all, round(r.lookup_seconds; digits=4),
    r.calls, r.hits, r.active, repr(r.primal), repr(r.dual_gap), VERSION, BLAS.get_num_threads(),
], ",")

results = map(measure, RUNS)

const OUT = joinpath(@__DIR__, "results_end_to_end.csv")
kept = isfile(OUT) ? filter(l -> !startswith(l, VARIANT * ","), readlines(OUT)[2:end]) : String[]
open(OUT, "w") do io
    println(io, COLUMNS)
    foreach(l -> println(io, l), kept)
    foreach(r -> println(io, row_string(r)), results)
end

for r in results
    println("$(r.variant) [$(r.source)] $(r.algorithm) $(r.problem): $(r.iterations) iterations, ",
        "$(round(r.seconds_min; digits=2)) s (runs $(r.seconds_all)), find_atom $(r.calls) calls / $(r.hits) hits ",
        "in $(round(r.lookup_seconds; digits=3)) s, active set $(r.active), primal $(r.primal)")
end
