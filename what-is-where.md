# What is where

```
Project.toml               dependencies: FrankWolfe, TimerOutputs, LinearAlgebra
Manifest.toml               the exact resolved versions; gitignored, see .gitignore
measurement/                a real BPCG run, instrumented, on problems where the
                             active set grows
  instrumentation.jl          times and counts find_atom without editing FrankWolfe.jl
  problems.jl                  the three problems: two Birkhoff sizes, one L∞-ball
  run.jl                        runs all three, writes results.csv
  results.csv                    max/mean active-set size, calls, lookup share, per problem
microbenchmark/              the lookup itself, isolated from any solver
  lookup_methods.jl            the linear scan (copied from find_atom) and a Dict lookup
  timing.jl                     warm up, batch to clear 1 ms, fastest of five
  run.jl                        sweeps size x dimension x atom scenario, writes results.csv
  results.csv                    scan and dict time per point, generic and adversarial atoms
README.md                  the question, the numbers, the answer; leads with the answer
README.fr.md                the same, in French (draft — see its own first line)
references.md               the papers and the issue, cited never redistributed
MEASURING.md                the machine, the noise, and what a number here does not claim
DECISIONS.md                what is Mohamed's to decide, including the draft issue comment
CITATION.cff                how to cite this repository
LICENSE                     MIT
.github/workflows/ci.yml    runs both scripts on every push; times nothing
explain/                    gitignored; HTML explainers for one reader, never shipped
```

Two directories, one house style: `measurement/` asks whether the lookup
costs anything in a real run; `microbenchmark/` asks how much the lookup
itself would cost either way, with no solver around it to hide the answer
in. Neither file runs the other — `measurement/run.jl` needs
`instrumentation.jl` and `problems.jl` only, `microbenchmark/run.jl` needs
`lookup_methods.jl` and `timing.jl` only — so either can be read, or rerun,
on its own.
