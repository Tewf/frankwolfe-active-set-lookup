# The method itself, usable on its own, apart from the benchmarks that
# measured it. Four files, one idea each:
#
#   certificate.jl  the lookup that does not search: absence proved from two
#                   inner products the algorithm already holds
#   keys.jl         computing a key from an atom (sparse pattern, dense value)
#   confirm.jl      the confirmation step that makes a bucket hit trustworthy
#   index.jl        the index structure: build, lookup, insert, delete-repair
#
# `certificate.jl` is the answer wherever the caller is a Frank-Wolfe step
# that just minimised the gradient over the active set (every call site in
# FrankWolfe.jl that reaches `find_atom` with no index); `index.jl` is the
# answer for a caller that has only the atom. `README.md` has the numbers,
# `METHOD.md` the argument for why both are correct, and `REJECTED.md` what
# they replaced. This file is only plumbing: it includes the four above and
# re-exports their public names in one place, so `using .ActiveSetLookup` is
# all a caller needs.
#
# A Julia package (`Project.toml` at the repository root): `Pkg.add(url=...)`
# then `using ActiveSetLookup`. `include`-ing this file and `using
# .ActiveSetLookup` also works, for a clone with nothing installed, which is
# how `guide/` loads it. The submodules below are each defined once; they
# reach each other with `using ..Name`, never by including a file twice.
module ActiveSetLookup

include(joinpath(@__DIR__, "keys.jl"))
include(joinpath(@__DIR__, "confirm.jl"))
include(joinpath(@__DIR__, "index.jl"))
include(joinpath(@__DIR__, "certificate.jl"))

using .AtomKeys, .AtomConfirm, .AtomIndexing, .AtomCertificate

export DEFAULT_K, atom_key, confirm_match, bucket_health,
    AtomIndex, SparsePatternIndex, DenseValueIndex,
    build_index, lookup_atom, push_atom!, delete_atom!,
    certified_absent, certified_lookup, scan_atoms

end # module
