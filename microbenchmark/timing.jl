# A minimal warm-up-then-batch timer, in the spirit of tensor-rank-toolkit's
# "fastest of three runs" (../MEASURING.md): one call to compile, then
# repeated batches sized to clear a floor so timer resolution and function
# call overhead do not dominate a single-call measurement, and the minimum
# per-call time over the batches is kept — the same reasoning as fastest-of-N
# for wall clock: nothing makes a run finish sooner than the work it does.
module Timing

export time_per_call

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

end # module
