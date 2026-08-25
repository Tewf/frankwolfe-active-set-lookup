# How a number here was measured

Every timing in this repository was produced the same way. This file states
it once; no measurement restates it.

## The machine, honestly

**12th Gen Intel Core i5-12450H, 12 threads, Julia 1.12.7.** Every run in
this repository was taken on this laptop **while it was not idle**: at the
time of the runs behind `measurement/results.csv` and
`microbenchmark/results.csv`, `uptime` reported a load average of **2.36**,
a Brave browser with several renderer processes was open, and two other
Claude Code agents were building unrelated repositories on the same machine.
No step was taken to quiet the machine before measuring, because none of
tensor-rank-toolkit's "one core, otherwise quiet" protocol was practical to
reproduce here without stopping work the other agents were mid-way through.
This is the disclosed exception `tensor-rank-toolkit/MEASURING.md` asks for
in place of a silent one: **the numbers here carry more noise than a single
quiet core would give**, and the one place that noise visibly shows
(whether the microbenchmark's dim=8192 crossover falls at 6,500 or 10,000)
is called out in `README.md` rather than smoothed over.

## Fastest of five, not the mean

Same reasoning as tensor-rank-toolkit's: a slow run measures what else the
machine was doing, not the code, so the minimum over repeated runs is the
closest available estimate of the work itself. Three runs is their floor;
this machine is noisier, so `microbenchmark/timing.jl` takes the **fastest
of five** batches per point instead, after one call to force compilation.

## Sub-microsecond calls are batched, not timed singly

`find_atom` on a small active set returns in tens of nanoseconds, well under
`@elapsed`'s useful resolution. `microbenchmark/timing.jl` runs the call in
a tight loop, doubling the batch size until one batch clears **1 ms**, then
takes the fastest of five such batches and divides by the batch size. This
is the same idea as `BenchmarkTools.jl`'s tuning phase, hand-rolled so the
whole sweep (96 points, two atom scenarios) finishes in under fifteen
seconds instead of the many minutes `@benchmark`'s default budget would take
per point.

## What the harness timer actually measures

`measurement/instrumentation.jl` adds a method for `find_atom` on
`FrankWolfe.jl`'s own `ActiveSet` type, more specific than the package's
method on the abstract `AbstractActiveSet`, so Julia's dispatch picks it for
every run in this repository, and every other active-set type ships
untouched. **`TimerOutputs.jl` was chosen over `Profile`** because the
question is "what share of one specific, known call's time", which a
manually-placed `@timeit` answers directly and cheaply (nanoseconds of
overhead per call); a sampling profiler would need a much longer run to
resolve a call this fast and answers a different question ("where does time
go", not "how much of it"). The instrumentation reproduces `find_atom`'s
published loop line for line (verified against
`~/.julia/packages/FrankWolfe/*/src/active_set.jl` at the version pinned in
`Manifest.toml`), so timing it does not change what it does.

## Warm-up, always

Every measured run, harness and microbenchmark alike, is preceded by an
identical, discarded call so JIT compilation lands outside the timed region.
`measurement/run.jl` runs each problem at a small iteration budget first;
`microbenchmark/timing.jl` runs the timed function once before its first
batch.

## What is reproducible, and what is not

**Counts are**: active-set sizes and `find_atom` call counts are integers
read off a real run and do not depend on the machine. **Timings are not**:
none are asserted by `.github/workflows/ci.yml`, which runs the scripts only
to confirm they still execute. Re-running `measurement/run.jl` or
`microbenchmark/run.jl` on a quieter machine will move the timing columns
and may move the noisy dim=8192 crossover; it should not move the active-set
counts or which problems keep them small.
