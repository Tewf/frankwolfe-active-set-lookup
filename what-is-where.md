# What is where

```
Project.toml               dependencies: FrankWolfe, TimerOutputs, LinearAlgebra
Manifest.toml               the exact resolved versions; gitignored, see .gitignore
measurement/                a real BPCG run, instrumented, on problems where the
                             active set grows
  instrumentation.jl          times, counts, and counts hits for find_atom without
                               editing FrankWolfe.jl
  problems.jl                  the three problems: two Birkhoff sizes, one L∞-ball
  run.jl                        runs all three, writes results.csv
  results.csv                    max/mean active-set size, calls, hits, lookup share,
                                  per problem
microbenchmark/              the lookup itself, isolated from any solver
  lookup_methods.jl            linear scan, full-atom Dict, and prefix-hash lookup
                                (dense and sparse-atom variants), copied from find_atom
  timing.jl                     warm up, batch to clear 1 ms, fastest of five, plus a
                                 variant that cycles through a pre-generated query mix
  run.jl                        sweeps size x dimension x atom scenario, writes results.csv
  results.csv                    scan and dict time per point, generic and adversarial atoms
  run_prefix.jl                  sweeps size x k x query type x atom alphabet against
                                  FrankWolfe.jl's own atom shapes, writes the three
                                  results_prefix_*.csv files below
  results_prefix_timing.csv        scan and prefix-hash time per point
  results_prefix_collisions.csv    bucket collision stats per point
  results_prefix_crossover.csv     smallest size where each prefix-hash variant beats
                                    the scan, or none
README.md                  the question, the numbers, the answer; leads with the answer
README.fr.md                the same, in French, updated for the new answer; the
                             "Prefix hashing" section's detail was not translated,
                             see DECISIONS.md
references.md               the papers and the issue, cited never redistributed
MEASURING.md                the machine, the noise, and what a number here does not claim
DECISIONS.md                what is Mohamed's to decide, including the draft issue comment
CITATION.cff                how to cite this repository
LICENSE                     MIT
.github/workflows/ci.yml    runs all three scripts on every push; times nothing
explain/                    gitignored; HTML explainers for one reader, never shipped
```

Two directories, one house style: `measurement/` asks whether the lookup
costs anything in a real run; `microbenchmark/` asks how much the lookup
itself would cost either way, with no solver around it to hide the answer
in. None of the three scripts needs another script's own files
(`measurement/run.jl` needs `instrumentation.jl` and `problems.jl` only;
`microbenchmark/run.jl` needs `lookup_methods.jl` and `timing.jl` only;
`microbenchmark/run_prefix.jl` needs the same two, plus the `FrankWolfe`
package itself to generate atoms in the exact shape `measurement/problems.jl`
uses, not to call any of `measurement/`'s own code), so any of the three
can be read, or rerun, on its own.
