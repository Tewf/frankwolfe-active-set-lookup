# Nobody tested the lifecycle before this branch. `active_set_cleanup!`
# removes atoms with `deleteat!`, which shifts every later index down by
# one, and `bucket_delete_repair!` (bucket_lifecycle.jl) is the naive
# whole-index walk that repairs a bucket map for that shift. This file is
# the test that would catch someone "optimising" that repair later (an
# indirection layer instead of a full walk, say) and getting the shift
# wrong: it runs randomised sequences of `push!`/`deleteat!` against a
# real index (built via test_atom_generators.jl's `route_build`, the same
# routing test_equivalence.jl and test_dispatch.jl use) and, after *every*
# single operation, checks the structure still agrees with a fresh linear
# scan over whatever the atoms vector currently holds.
#
# Deletions are forced to about half of all operations, deliberately far
# above the real BPCG rate `measurement/results.csv` recorded (2, 1, and 0
# deletions across 8,002-20,002 iterations): the point here is to exercise
# `bucket_delete_repair!` hard, not to reproduce how rarely a real run
# calls it (README.md's "Lookup is not the whole cost" already covers
# that question).
#
# On failure, the exact operation sequence is not replayable from the
# printed message alone, but the master seed is: rerunning this file
# reproduces every subseed, and the printed step number/last op pins down
# where in that one sequence it broke.
#
# Run: julia --project=. microbenchmark/test_lifecycle.jl

using Test, Random

include(joinpath(@__DIR__, "test_atom_generators.jl"))
include(joinpath(@__DIR__, "bucket_lifecycle.jl"))
using .TestAtomGenerators, .BucketLifecycle

const MASTER_SEED = 6 # distinct from run.jl(1)..run_pattern_key_reps.jl(4), test_equivalence.jl(5)

const NUM_OPS = 200
const K_VALUES = [4, 8] # the recommended k and the earlier default; see README.md's k-sweep
const N_SEQUENCES = 2
const START_SIZE = 5

function report_and_check(cond::Bool, context::String)
    cond || println(stderr, "LIFECYCLE FAILURE (master_seed=$MASTER_SEED): ", context)
    @test cond
end

# One randomised insert/delete sequence for one (alphabet, k), checked
# after every single operation. Returns nothing; all assertions go through
# `report_and_check` above.
function run_sequence(spec, k::Int, seed::UInt64)
    rng = Random.Xoshiro(seed)
    atoms = atom_pool(spec, START_SIZE, rng)
    index = route_build(atoms, k)
    last_op = "initial build, size=$START_SIZE"

    for step in 1:NUM_OPS
        do_insert = isempty(atoms) || rand(rng) < 0.5
        local check_positions
        if do_insert
            new_atom = spec.atom(rng)
            push!(atoms, new_atom)
            bucket_insert!(index.buckets, route_key(spec, new_atom, k), length(atoms))
            last_op = "step=$step insert -> size=$(length(atoms))"
            n = length(atoms)
            check_positions = unique(vcat(1, n, rand(rng, 1:n, min(3, n))))
        else
            pos = rand(rng, 1:length(atoms))
            deleteat!(atoms, pos)
            bucket_delete_repair!(index.buckets, pos)
            last_op = "step=$step delete pos=$pos -> size=$(length(atoms))"
            n = length(atoms)
            check_positions = n == 0 ? Int[] : unique(filter(p -> 1 <= p <= n, vcat(1, n, pos - 1, pos, rand(rng, 1:n, min(3, n)))))
        end

        ctx(label) = "alphabet=$(spec.name) k=$k seed=$seed last_op=[$last_op] ($label)"

        for pos in check_positions
            query = atoms[pos]
            expected = route_scan(atoms, query)
            report_and_check(expected != -1, ctx("post-op present pos=$pos scan-sanity"))
            got = route_lookup(index, atoms, query)
            report_and_check(got == expected, ctx("post-op present pos=$pos"))
        end

        if !isempty(atoms)
            absent = spec.atom(rng)
            expected_absent = route_scan(atoms, absent)
            got_absent = route_lookup(index, atoms, absent)
            report_and_check(got_absent == expected_absent, ctx("post-op absent"))
        end
    end
end

@testset "lifecycle: insert/delete sequences never desync from the scan" begin
    rng = Random.Xoshiro(MASTER_SEED)
    for spec in ALPHABETS, k in K_VALUES, rep in 1:N_SEQUENCES
        seed = rand(rng, UInt64)
        run_sequence(spec, k, seed)
    end
end
