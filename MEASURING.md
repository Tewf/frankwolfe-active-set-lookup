# How a number here was measured

Every timing in this repository was produced the same way. This file states
it once; no measurement restates it.

## The machine, honestly

**12th Gen Intel Core i5-12450H, 12 threads, Julia 1.12.7.** Every run in
this repository was taken on this laptop **while it was not idle**: at the
time of the runs behind `measurement/results.csv` and
`microbenchmark/results.csv`, `uptime` reported a load average of **2.36**,
a Brave browser with several renderer processes was open, and two other
unrelated builds were running on the same machine.
No step was taken to quiet the machine before measuring, because none of
tensor-rank-toolkit's "one core, otherwise quiet" protocol was practical to
reproduce here without stopping the builds that were mid-way through.
This is the disclosed exception `tensor-rank-toolkit/MEASURING.md` asks for
in place of a silent one: **the numbers here carry more noise than a single
quiet core would give**, and the one place that noise visibly shows
(whether the microbenchmark's dim=8192 crossover falls at 6,500 or 10,000)
is called out in `README.md` rather than smoothed over.

The run behind `measurement/results.csv`'s `find_atom_hits` column and
`microbenchmark/results_prefix_*.csv` was taken on the same laptop, also
not idle: `uptime` reported a load average of **0.6-1.0** across the two
runs, lower than the 2.36 above but still real background load (a browser,
an editor, other processes). Both this run's noise floor and the fix for a
specific bias it uncovered are in "Index construction can bias whichever
measurement runs right after it," below.

The run behind `microbenchmark/results_pattern_key_reps_*.csv`
(the representation stage) was taken on the same laptop, also not idle:
`uptime` reported a load average of **1.00, 0.63, 0.61** immediately
before the run started. `Random.seed!(4)` (see `run_pattern_key_reps.jl`'s
own header) makes the exact atom sequence reproducible across runs on the
same machine regardless of that load, which is how the `NTuple{4,Int}`
insert-timing anomaly noted in `DECISIONS.md` was confirmed to hold at the same magnitude (177.65 ns for
Birkhoff n=25) across two independent full runs, ruling out ordinary
timer noise as its explanation even though the machine itself was not
quiet for either run.

The runs behind `measurement/results.csv`'s certificate columns,
`measurement/results_algorithms.csv` and
`microbenchmark/results_certificate_*.csv` were taken on 2026-08-27 on the
same laptop, load average **0.9-1.1** with a browser open and nothing else
computing. Two things about them belong here. First, the harness's
`RecordingLMO` copies the direction on every LMO call; that copy is inside
the timed total and outside the timed scan, so it can only shrink every
run's `lookup_share`, and equally. Second, `run_certificate.jl` re-times
the folded key, the prefix hash and the scan on its own atoms in the same
session as the certificate, and those numbers differ from the same
structures' committed sweeps (28.0 ns against 15.6 ns for the folded key's
miss at n=60): that gap is the session-to-session noise this file already
warns about, and it is why every comparison with the certificate is made
within one session rather than against a committed number. Non-lazy
pairwise Frank-Wolfe and blended conditional gradients on Birkhoff n=60
took 66 s and 130 s, pairwise Frank-Wolfe on the L-inf ball 25 s; the
other nine runs took under five seconds each.

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

## Index construction can bias whichever measurement runs right after it

`microbenchmark/run_prefix.jl` first surfaced this: building a prefix-hash
index allocates heavily (up to one bucket vector per stored atom), and a
first pass at that sweep showed the scan measured immediately afterwards
coming out slower than the same scan measured in isolation, by enough to
flip a crossover that a controlled re-test (same atoms, same query, scan
and prefix timed in both orders) showed did not exist: the garbage from
index construction was getting collected mid-measurement, and it landed on
whichever method happened to be timed next rather than on the method that
allocated it. `run_prefix.jl` now calls `GC.gc()` once after building each
group's indices and once per query type, before either method's timing
starts, so every measurement in a group starts from the same clean heap.
This does not appear in `microbenchmark/run.jl`'s original sweep because
its `Dict` is built once per point, not once per `(size, k)` pair, so the
same allocation burst is much smaller relative to the timed work.

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

**The certificate tally is instrumentation, and it is timed apart.** The
same override also evaluates the certificate on every call and compares it
with the scan, which costs one pass over the active set per call. Until
2026-08-28 that pass sat inside the `"total"` timer: pairwise Frank-Wolfe
at n=60 measured 98.8 s with it and 68.3 s without, and every share in the
README table was a ratio to the inflated figure. It is now under its own
`"certificate_tally"` timer and `run_time_of` subtracts it; the corrected
shares are in `DECISIONS.md`. **End-to-end seconds carry the session's
noise**: `measurement/run_end_to_end.jl` runs identical code twice (the
registry release and master, same iterates on pairwise n=60) and got
68.3 s and 71.4 s as fastest-of-three, so a difference between variants
below about 3 s means nothing here, and the counts and the final objective
are the part of that file to trust.

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
