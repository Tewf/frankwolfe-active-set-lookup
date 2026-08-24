# Runs BPCG on each problem in `problems.jl` under the instrumentation in
# `instrumentation.jl`, and writes one row per run to `results.csv`: the
# problem, its dimension, the iteration budget, the maximum and mean active-set
# size seen by the callback, how many times `find_atom` was called, and the
# share of total run time its accumulated time represents.
#
# Every problem is run twice: once to compile (small iteration budget,
# discarded), once timed. See ../MEASURING.md for why, and for what else was
# running on the machine while this was measured.
include("instrumentation.jl")
include("problems.jl")
using .LookupInstrumentation
using .Problems
using FrankWolfe, TimerOutputs

function run_once(problem; max_iteration, epsilon=1e-9)
    sizes = Int[]
    reset_instrumentation!()
    @timeit TIMER "total" begin
        FrankWolfe.blended_pairwise_conditional_gradient(
            problem.f,
            problem.grad!,
            problem.lmo,
            problem.x0;
            max_iteration=max_iteration,
            lazy=true,
            epsilon=epsilon,
            callback=active_set_sizes!(sizes),
            verbose=false,
        )
    end
    return sizes
end

function measure(problem; max_iteration, warmup_iteration=50)
    run_once(problem; max_iteration=warmup_iteration) # compiles, discarded
    sizes = run_once(problem; max_iteration=max_iteration)
    total_ns = TimerOutputs.time(TIMER["total"])
    return (
        problem=problem.label,
        dimension=problem.dimension,
        max_iteration=max_iteration,
        iterations_run=length(sizes),
        max_active_set=maximum(sizes),
        mean_active_set=sum(sizes) / length(sizes),
        find_atom_calls=lookup_calls(),
        total_seconds=total_ns / 1e9,
        lookup_share=lookup_share_of(TIMER["total"]),
    )
end

runs = [
    (birkhoff_problem(25; seed=1), 8000),
    (birkhoff_problem(60; seed=1), 20000),
    (linf_box_problem(3000; seed=1), 15000),
]

results = [measure(problem; max_iteration=k) for (problem, k) in runs]

open(joinpath(@__DIR__, "results.csv"), "w") do io
    println(
        io,
        "problem,dimension,max_iteration,iterations_run,max_active_set,mean_active_set,find_atom_calls,total_seconds,lookup_share",
    )
    for r in results
        println(
            io,
            "$(r.problem),$(r.dimension),$(r.max_iteration),$(r.iterations_run),$(r.max_active_set),$(round(r.mean_active_set, digits=2)),$(r.find_atom_calls),$(round(r.total_seconds, digits=4)),$(round(r.lookup_share, sigdigits=4))",
        )
    end
end

for r in results
    println(
        "$(r.problem): max=$(r.max_active_set) mean=$(round(r.mean_active_set,digits=1)) ",
        "calls=$(r.find_atom_calls) total=$(round(r.total_seconds,digits=3))s ",
        "lookup_share=$(round(100*r.lookup_share,digits=3))%",
    )
end
