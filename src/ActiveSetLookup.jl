# The method itself, usable on its own, apart from the benchmarks that
# measured it. Three files, one idea each:
#
#   keys.jl       computing a key from an atom (sparse pattern, dense value)
#   confirm.jl    the confirmation step that makes a bucket hit trustworthy
#   index.jl      the index structure: build, lookup, insert, delete-repair
#
# `README.md` has the number, `METHOD.md` has the argument for why this is
# correct, and `REJECTED.md` has what this replaced. This file is only
# plumbing: it includes the three above and re-exports their public names
# in one place, so `using .ActiveSetLookup` is all a caller needs.
#
# Not a registered Julia package (no `uuid`/`version` in `Project.toml`):
# this repository already runs everything by `include`-ing a file and
# `using` its module, `measurement/` and `microbenchmark/` alike, and this
# module follows the same convention rather than introducing a second way
# to load code.
module ActiveSetLookup

include(joinpath(@__DIR__, "keys.jl"))
include(joinpath(@__DIR__, "confirm.jl"))
include(joinpath(@__DIR__, "index.jl"))

using .AtomKeys, .AtomConfirm, .AtomIndexing

export DEFAULT_K, atom_key, confirm_match, bucket_health,
    AtomIndex, SparsePatternIndex, DenseValueIndex,
    build_index, lookup_atom, push_atom!, delete_atom!

end # module
