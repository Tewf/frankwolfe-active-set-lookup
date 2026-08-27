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

Also counts `deleteat!(active_set, idx::Int)` calls the same way, for the
same reason: the lifecycle stage's brief points out that every
`find_atom` miss is followed by a `push!`, and the active set also has to
survive `active_set_cleanup!`'s `deleteat!`, which shifts every later index
down by one. `active_set.jl`'s own scalar `deleteat!` (line 66) is what
both the direct drop-step call (`active_set.jl:177`) and the batch cleanup
path (`active_set.jl:58`'s loop over a sorted index vector) end up calling,
so overriding just the scalar method the same way `find_atom` is overridden
catches both call sites without editing the installed package.

Also evaluates, at every `find_atom` call made for the LMO's own vertex, the
absence certificate `src/certificate.jl` proposes in place of the search:
`<g,v> < min_a <g,a>` proves `v` is not active. `find_atom` never sees the
gradient, so `RecordingLMO` below wraps the problem's LMO and keeps a copy of
the last direction it was asked to minimise, which is the gradient of the
iteration the lookup happens in. The certificate is then computed, outside
the timed scan, with the package's own orientations (`dot(gradient, v)` as in
`blended_pairwise.jl`, `dot(atom, direction)` as in `active_set_argminmax`),
and compared with what the scan answered. `certificate_counts()` reports how
often it decided, how often a tie was resolved by the best atom alone, how
often a tie needed a search, and whether it ever contradicted the scan.
"""
module LookupInstrumentation

using FrankWolfe
using LinearAlgebra
using TimerOutputs

export TIMER,
    RecordingLMO,
    RECORDER,
    reset_instrumentation!,
    lookup_calls,
    lookup_hits,
    lookup_share_of,
    active_set_sizes!,
    deletion_calls,
    certificate_counts

const TIMER = TimerOutput()
const CALLS = Ref(0)
# HITS counts calls where the atom was already in the active set (idx != -1):
# a re-encountered vertex, which forces `_unsafe_equal` to run to completion
# rather than exit early. CALLS - HITS is the miss count. Added to answer
# "how often is find_atom's answer actually a hit", which the original
# instrumentation never recorded and microbenchmark/'s prefix-hash sweep
# needs as its "realistic mix" ratio.
const HITS = Ref(0)
# One count per scalar deleteat!(active_set, idx::Int) call, which is one
# atom actually removed and every later atom's stored position shifted
# down by one: active_set_cleanup!'s batch delete (a sorted index vector)
# calls this once per index it removes, so DELETIONS counts individual
# removals, not cleanup invocations.
const DELETIONS = Ref(0)

# The certificate's own tally, one counter per outcome. `certified` calls
# are the ones where `<g,v> < min` held and no search was needed;
# `tie_best` are ties where the best atom itself was the query (a hit found
# in one comparison); `tie_other` are ties with a different atom, the only
# case that needs a search; `inverted` means `<g,v>` exceeded the minimum,
# which an exact LMO cannot produce against the direction it was given and
# is counted so a stale or inexact direction shows up rather than hides;
# `other_caller` counts `find_atom` calls that were not for the recorded
# vertex; `unsound` counts any certified call the scan contradicted, and
# must stay zero.
const CERT = Dict{Symbol,Int}(
    :calls => 0, :certified => 0, :tie_best => 0, :tie_other => 0,
    :inverted => 0, :other_caller => 0, :unsound => 0,
)

"""
    RecordingLMO(inner)

Delegates every `compute_extreme_point` to `inner` and keeps a copy of the
direction it was called with and the vertex it returned, so the certificate
can be evaluated later with the gradient the step actually used.
"""
mutable struct RecordingLMO{L} <: FrankWolfe.LinearMinimizationOracle
    inner::L
    direction::Any
    vertex::Any
end
RecordingLMO(inner) = RecordingLMO(inner, nothing, nothing)

function FrankWolfe.compute_extreme_point(lmo::RecordingLMO, direction; kwargs...)
    v = FrankWolfe.compute_extreme_point(lmo.inner, direction; kwargs...)
    lmo.direction = copy(direction)
    lmo.vertex = v
    return v
end

# The run script points this at the `RecordingLMO` it built, since
# `find_atom` receives the active set and the atom and nothing else.
const RECORDER = Ref{Any}(nothing)

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
    tally_certificate!(active_set, atom, idx)
    return idx
end

function tally_certificate!(active_set, atom, scan_idx)
    recorder = RECORDER[]
    (recorder === nothing || recorder.direction === nothing || isempty(active_set)) && return nothing
    if atom !== recorder.vertex
        CERT[:other_caller] += 1
        return nothing
    end
    g = recorder.direction
    query_value = dot(g, atom)
    best, best_value = -1, typemax(Float64)
    @inbounds for i in eachindex(active_set)
        val = dot(active_set.atoms[i], g)
        if val < best_value
            best_value, best = val, i
        end
    end
    CERT[:calls] += 1
    if query_value < best_value
        CERT[:certified] += 1
        scan_idx == -1 || (CERT[:unsound] += 1)
    elseif query_value == best_value
        if FrankWolfe._unsafe_equal(active_set.atoms[best], atom)
            CERT[:tie_best] += 1
            scan_idx == best || (CERT[:unsound] += 1)
        else
            CERT[:tie_other] += 1
        end
    else
        CERT[:inverted] += 1
    end
    return nothing
end

function Base.deleteat!(active_set::FrankWolfe.ActiveSet, idx::Int)
    DELETIONS[] += 1
    deleteat!(active_set.atoms, idx)
    deleteat!(active_set.weights, idx)
    return active_set
end

function reset_instrumentation!()
    CALLS[] = 0
    HITS[] = 0
    DELETIONS[] = 0
    for key in keys(CERT)
        CERT[key] = 0
    end
    reset_timer!(TIMER)
    return nothing
end

lookup_calls() = CALLS[]
lookup_hits() = HITS[]
deletion_calls() = DELETIONS[]
certificate_counts() = copy(CERT)

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

Returns a callback (`(state, active_set, args...) -> true`) that appends
`length(active_set)` on every iteration. Every call shape the three
algorithms use — `callback(state, active_set)` on an FW step,
`callback(state, active_set, a)` on a BPCG pairwise step,
`callback(state, active_set, non_simplex_iter)` in BCG — matches `args...`.
"""
function active_set_sizes!(sizes::Vector{Int})
    return (state, active_set, args...) -> begin
        push!(sizes, length(active_set))
        return true
    end
end

end # module
