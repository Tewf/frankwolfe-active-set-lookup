# Contributing, and how this repository is written

For two readers: a stranger who wants the lookup in their own Frank-Wolfe
code, and the FrankWolfe.jl maintainers deciding what, if anything, to take
upstream. The conventions are stated here once; other files point here
rather than restating them.

## What this repository is

A Julia package, `ActiveSetLookup` (`src/`, `test/`, `Project.toml`), and
the measurement campaign that justifies it (`measurement/`,
`microbenchmark/`, the committed `results*.csv`). The package answers one
question a Frank-Wolfe step asks its active set; the campaign shows, on
real runs of FrankWolfe.jl, that the answer is right and what it costs.
`guide/` is neither: a teaching reimplementation for a reader who has
never opened the library. `README.md` has the answer, the usage and the
reproduction commands; `what-is-where.md` names every file in one line.

## Conventions

- **One file, one role**, describable in a sentence without an "and".
  Structure (`what-is-where.md`), argument (`METHOD.md`), refusals
  (`REJECTED.md`), judgement calls (`DECISIONS.md`), test coverage
  (`TESTING.md`) and timing method (`MEASURING.md`) are separate files.
- **Every fact lives in one place.** Numbers sit in the CSVs and are quoted
  from them; a claim is stated where it is proved and linked from
  everywhere else, so correcting a fact means editing one file.
- **Markdown files stay under about 80 lines**, so each is read in one
  sitting and is split when its role grows a second clause. A file that
  exceeds the ceiling says so in its first paragraph, with the reason:
  `README.md` (a stranger's entry point needs answer, usage and
  reproduction on one page), `METHOD.md` and `REJECTED.md` (each argument
  or refusal carries its own numbers), `guide/README.md` (a walkthrough
  that links out before the picture is complete is not a walkthrough).
  Such a note marks a deliberate exception, not an apology.
- **Names are full words**; nothing that needs a glossary.
- **Nothing is claimed that is not measured or tested.** Every timing comes
  from a committed `results*.csv` whose provenance is in `MEASURING.md`.
  CI runs every script end to end but asserts no timing, since a shared
  runner cannot reproduce them; the one campaign fact it asserts is that
  the certificate never contradicted the scan.
- **Comments say why the code is as it is; commit messages say why it
  changed.** Subject under 72 characters, a body whenever the diff is real.

## Working on it

Julia 1.10 or later. `julia --project=. -e 'using Pkg; Pkg.test()'` runs
Aqua's package-quality checks and the three suites under `test/`. The
campaign runs in its own environments, `measurement/` and
`microbenchmark/`, so FrankWolfe.jl never becomes a dependency of the
package; `README.md`, "Reproducing the measurements", has every command.

A change to the method comes with a test in `test/` that would have failed
without it. A new lookup strategy goes beside the others in
`microbenchmark/lookup_methods.jl` and into the sweep that compares them.
A new measurement writes its CSV beside its script and gets a row in the
README table, with commit, machine and flags recorded in `MEASURING.md`.
`Project.toml` and `CITATION.cff` carry the same version; a release is a
tag `vX.Y.Z` on the commit that bumps both.

## For FrankWolfe.jl maintainers

What is proposed upstream is the certificate alone, a few lines at each
call site, not this package. At `blended_pairwise.jl`, `pairwise.jl`,
`blended_cg.jl` and the corrective and block-coordinate drivers: pass `-1`
to `active_set_update!` when `dot(gradient, v) < dot_forward_vertex`, both
already in scope, and fall back to `find_atom` otherwise; keep the minimum
`pfw_step` currently discards; have `lp_separation_oracle` return the
position beside the atom it removes. `METHOD.md` has the argument and the
two floating-point facts it rests on, `test/test_certificate.jl` checks
both on whatever machine runs it, and `measurement/run.jl` re-runs the
comparison against whatever FrankWolfe.jl version resolves, so every
number here can be regenerated on any release. The index (`src/index.jl`)
serves a caller that reaches `find_atom` with only an atom; `DECISIONS.md`
says why it should stay optional, offered the way
`ActiveSetQuadraticProductCaching` is.

Issues and pull requests are welcome. A measurement that contradicts a
number here is the most useful thing to send.
