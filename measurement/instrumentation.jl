"""
Times and counts every call `find_atom` (`active_set.jl`, around line 316)
makes on a plain `FrankWolfe.ActiveSet`, without editing the installed
package. Adding a method for the concrete `ActiveSet` type is more specific,
under Julia's dispatch rules, than the package's own method for the abstract
`AbstractActiveSet` it defines `find_atom` on — so this method runs instead of
that one, for exactly the runs in this repository, and every other active-set
type (`ActiveSetQuadraticProductCaching`, etc.) is untouched. The loop below
is copied from the package verbatim; nothing about what it does changes,
only that each call is now counted and timed.
"""
module LookupInstrumentation

using FrankWolfe
using TimerOutputs

export TIMER,
    reset_instrumentation!, lookup_calls, lookup_hits, lookup_share_of, active_set_sizes!

const TIMER = TimerOutput()
const CALLS = Ref(0)
# HITS counts calls where the atom was already in the active set (idx != -1):
# a re-encountered vertex, which forces `_unsafe_equal` to run to completion
# rather than exit early. CALLS - HITS is the miss count. Added to answer
# "how often is find_atom's answer actually a hit", which the original
# instrumentation never recorded and microbenchmark/'s prefix-hash sweep
# needs as its "realistic mix" ratio.
const HITS = Ref(0)

function FrankWolfe.find_atom(active_set::FrankWolfe.ActiveSet, atom)
    CALLS[] += 1
    idx = -1
    @timeit TIMER "find_atom" begin
        @inbounds for i in eachindex(active_set)
            if FrankWolfe._unsafe_equal(active_set.atoms[i], atom)
                idx = i
                break
            end
        end
    end
    if idx != -1
        HITS[] += 1
    end
    return idx
end

function reset_instrumentation!()
    CALLS[] = 0
    HITS[] = 0
    reset_timer!(TIMER)
    return nothing
end

lookup_calls() = CALLS[]
lookup_hits() = HITS[]

"""
    lookup_share_of(total_section)

Fraction of `total_section`'s accumulated time spent inside the nested
`"find_atom"` timer, or `0.0` if `find_atom` was never called under it.
"""
function lookup_share_of(total_section)
    if !haskey(total_section.inner_timers, "find_atom")
        return 0.0
    end
    return TimerOutputs.time(total_section["find_atom"]) / TimerOutputs.time(total_section)
end

"""
    active_set_sizes!(sizes::Vector{Int})

Returns a BPCG-compatible callback (`(state, active_set, args...) -> true`)
that appends `length(active_set)` on every iteration. Both BPCG call shapes —
`callback(state, active_set)` on an FW step and `callback(state, active_set,
a)` on a pairwise step — match `args...`.
"""
function active_set_sizes!(sizes::Vector{Int})
    return (state, active_set, args...) -> begin
        push!(sizes, length(active_set))
        return true
    end
end

end # module
