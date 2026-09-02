# Changelog

Versions are `Project.toml`'s and `CITATION.cff`'s, tagged `vX.Y.Z`. Dates
are the commit that made the change; `DECISIONS.md` has the reasoning
behind each entry, `REJECTED.md` what each one replaced.

## Unreleased

- `certified_lookup` compares the best atom to the query before it consults
  the certificate. The old order certified a present atom absent whenever it
  was stored in a different representation than the query, because the two
  `dot` calls are then different methods and can differ by an ulp: 98 of 200
  seeded Float32-stored/Float64-queried lookups, and 10 of 50 on dense copies
  of sparse Birkhoff vertices. `test/test_certificate.jl` covers both.
- `confirm_match` gains the generic `isequal` method that `_unsafe_equal` has
  beside its dense and sparse ones, so a pair crossing families is compared
  instead of raising a `MethodError`.
- `measurement/run_end_to_end.jl` and `results_end_to_end.csv`: pairwise
  Frank-Wolfe and BPCG timed end to end, FrankWolfe.jl master against the
  branch that applies the certificate at its call sites. At n=60, 71.4 s
  with 7.3 s of scanning becomes 66.3 s with none, same atoms, same
  objective; the README's "about 2 ms by arithmetic" is gone.
- The harness's certificate tally is timed on its own and excluded from
  run time. It had inflated every total, so the shares in the README table
  were understated: pairwise n=60 is 10.4% of a 66 s run, not 6.4% of 99 s.
- The certificate is in FrankWolfe.jl. ZIB-IOL/FrankWolfe.jl#649, from the
  branch `certificate-244`, was merged on 2026-09-02 as `6b72bd97` and closed
  issue #244: `find_atom(active_set, atom, direction, best_index, best_value)`
  at the BPCG, pairwise, blended-CG and block-coordinate call sites,
  `lp_separation_oracle` returning the position beside the atom, and the
  library's own tests passing unchanged.
- This file, and a native pass on `README.fr.md`.

## 0.2.0, 2026-08-27 (tagged 2026-08-28)

- The absence certificate, `src/certificate.jl`: `<g,v> < <g,s>` proves the
  LMO's vertex absent from the active set in one comparison; every call of
  every measured run was decided without a search. Leads every document.
- Measured on four algorithms and three problems (`measurement/`), timed
  against the key, the prefix hash and the scan (`microbenchmark/`).
- On 2026-08-28: `guide/`, a reimplementation of both algorithms for a
  reader who has never opened the library; `CONTRIBUTING.md`; the root
  turned into the package `ActiveSetLookup` with `Pkg.test()`, Aqua and a
  docstring on every export; the harness moved to its own environments.

## 0.1.0, 2026-08-26 (never tagged)

- `src/`: the folded structural key with exact confirmation
  (`build_index`, `lookup_atom`, `push_atom!`, `delete_atom!`), dispatching
  on sparse and dense atoms, with the signed-zero canonicalisation.
- `METHOD.md`, `REJECTED.md`, `DECISIONS.md`, `TESTING.md`, `MEASURING.md`.

## Before 0.1.0, 2026-08-24 to 2026-08-25

- The instrumented harness on BPCG; full-atom hashing measured and
  refused; the sparse-pattern key, the trie and the integer
  representations swept.
