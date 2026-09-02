# The confirmation step: what makes a bucket hit trustworthy. `atom_key`
# (keys.jl) is allowed to be lossy (the UInt64 fold can collide) precisely
# because nothing here ever returns a bucket's contents as an answer on
# their own. A bucket is a list of *candidates*, and every one is checked
# against the whole atom with exact equality before `lookup_atom`
# (index.jl) returns it. A collision costs one extra comparison; it never
# costs a wrong answer. See METHOD.md for the full correctness argument.
#
# Mirrors `FrankWolfe.jl`'s own `_unsafe_equal` (active_set.jl), method for
# method: one for a dense `Array` (elementwise `!=`, so it follows `==`
# semantics, the same semantics behind `atom_key`'s `-0.0`/`0.0`
# canonicalisation), one for an `AbstractSparseArray` (`==`, SparseArrays'
# own comparison), and a generic `isequal` for every other pair, including a
# dense atom against a sparse one. Nothing here changes what "equal" means;
# it only decides who gets asked.
module AtomConfirm

using SparseArrays

export confirm_match

"""
    confirm_match(a, b) -> Bool

Exact equality of two atoms, the step that makes a bucket hit trustworthy.
One method per atom family, mirroring `FrankWolfe.jl`'s `_unsafe_equal`:
`==` for sparse arrays, elementwise `!=` for dense ones, so both follow
`==` semantics (`0.0 == -0.0`), and `isequal` for anything else, so a pair
that crosses families is compared rather than refused. Every candidate an
index returns passes through here before `lookup_atom` trusts it.
"""
confirm_match(a, b) = isequal(a, b)

confirm_match(a::AbstractSparseArray, b::AbstractSparseArray) = a == b

function confirm_match(a::Array, b::Array)
    length(a) == length(b) || return false
    @inbounds for i in eachindex(a)
        a[i] != b[i] && return false
    end
    return true
end

end # module
