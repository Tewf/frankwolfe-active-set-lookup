# A minimal warm-up-then-batch timer, in the spirit of tensor-rank-toolkit's
# "fastest of three runs" (../MEASURING.md): one call to compile, then
# repeated batches sized to clear a floor so timer resolution and function
# call overhead do not dominate a single-call measurement, and the minimum
# per-call time over the batches is kept — the same reasoning as fastest-of-N
# for wall clock: nothing makes a run finish sooner than the work it does.
module Timing

export time_per_call, time_per_call_seq

function time_per_call(f, args...; floor_seconds=1e-3, repeats=5)
    f(args...) # compile, discard
    batch = 1
    elapsed = 0.0
    while true
        elapsed = @elapsed for _ in 1:batch
            f(args...)
        end
        elapsed >= floor_seconds && break
        batch *= 4
    end
    best = elapsed / batch
    for _ in 2:repeats
        elapsed = @elapsed for _ in 1:batch
            f(args...)
        end
        best = min(best, elapsed / batch)
    end
    return best
end

# Like `time_per_call`, but the last argument varies across a pre-generated
# sequence instead of being fixed: needed for a "realistic mix" of hit and
# miss queries, where what matters is the *average* cost over a query
# distribution, not one query's cost repeated. `queries` is consumed
# round-robin (`mod1`) so a batch can run longer than the sequence without
# regenerating it mid-timing, which would make the fastest-of-five
# comparison noisier, not less so. The `mod1` indexing adds a few
# nanoseconds of constant overhead to every call; that overhead is the same
# for every method timed this way, so it cannot favour one over another,
# only shift all their absolute numbers up slightly.
function time_per_call_seq(f, fixed_args, queries; floor_seconds=1e-3, repeats=5)
    n = length(queries)
    f(fixed_args..., queries[1]) # compile, discard
    batch = 1
    elapsed = 0.0
    while true
        elapsed = @elapsed for i in 1:batch
            f(fixed_args..., queries[mod1(i, n)])
        end
        elapsed >= floor_seconds && break
        batch *= 4
    end
    best = elapsed / batch
    for _ in 2:repeats
        elapsed = @elapsed for i in 1:batch
            f(fixed_args..., queries[mod1(i, n)])
        end
        best = min(best, elapsed / batch)
    end
    return best
end

end # module
