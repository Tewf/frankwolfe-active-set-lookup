# The single invariant that matters for a hashed active-set lookup: for
# any atom set and any query, the structure returns exactly what
# `find_atom`'s linear scan would return. Everything else in this
# repository (the timing, the allocation counts) is only worth trusting
# once this holds, and holds beyond the handful of hand-built examples
# test_soundness.jl and test_pattern_key_reps.jl already check. This file
# checks it as a property: random seeds, active-set sizes from 0 to 500,
# every k this repository has swept (2, 4, 8, 16), and all three atom
# alphabets (test_atom_generators.jl's Birkhoff/L-infinity/generic),
# routed to whichever structure that alphabet actually gets
# (test_dispatch.jl checks the routing itself).
#
# Four query shapes per atom set, all required by the brief: present
# (first atom, last atom, and a few atoms from the middle), and absent (a
# fresh atom confirmed absent by the scan itself, not assumed). A second
# testset forces duplicates into the atom pool and checks the structure
# agrees with the scan on *which* index wins: the scan returns the first
# occurrence, and nothing about a bucketed structure guarantees that for
# free unless its bucket order (or its own confirm loop) preserves it.
#
# Run: julia --project=. microbenchmark/test_equivalence.jl

using Test, Random

include(joinpath(@__DIR__, "test_atom_generators.jl"))
using .TestAtomGenerators

const MASTER_SEED = 5 # distinct from run.jl(1)..run_pattern_key_reps.jl(4)

const SIZES = [0, 1, 2, 3, 5, 10, 25, 50, 100, 250, 500]
const KS = [2, 4, 8, 16]
const N_REPS = 3 # independent atom draws per (alphabet, size, k)

# Prints the failing case before the @test itself fails, so a CI log shows
# exactly which alphabet/size/k/seed to rerun rather than just "1 test
# failed". Kept as a plain function, not a macro: `@test` does not need
# `cond`'s own text to be meaningful here, since the context string below
# already says what mattered.
function report_and_check(cond::Bool, context::String)
    cond || println(stderr, "EQUIVALENCE FAILURE (master_seed=$MASTER_SEED): ", context)
    @test cond
end

@testset "structure agrees with the scan: seeds x sizes x k x alphabet" begin
    rng = Random.Xoshiro(MASTER_SEED)
    for spec in ALPHABETS, size in SIZES, k in KS, rep in 1:N_REPS
        subseed = rand(rng, UInt64)
        srng = Random.Xoshiro(subseed)
        atoms = atom_pool(spec, size, srng)
        index = route_build(atoms, k)
        ctx(label) = "alphabet=$(spec.name) size=$size k=$k subseed=$subseed ($label)"

        if size > 0
            # First, last, and up to three atoms from the middle: the
            # brief's explicit "equal to the first atom, equal to the last
            # atom" cases, plus a few more for coverage at larger sizes.
            positions = unique(vcat(1, size, rand(srng, 1:size, min(3, size))))
            for pos in positions
                query = atoms[pos]
                expected = route_scan(atoms, query)
                report_and_check(expected != -1, ctx("present pos=$pos scan-sanity"))
                got = route_lookup(index, atoms, query)
                report_and_check(got == expected, ctx("present pos=$pos"))
            end
        end

        # Absent: a fresh atom, confirmed absent by the scan itself first
        # (these alphabets are high-entropy enough that a false "absent"
        # is essentially impossible, but the check is against ground
        # truth, not an assumption about the alphabet).
        absent = spec.atom(srng)
        tries = 0
        while size > 0 && route_scan(atoms, absent) != -1 && tries < 5
            absent = spec.atom(srng)
            tries += 1
        end
        expected_absent = route_scan(atoms, absent)
        got_absent = route_lookup(index, atoms, absent)
        report_and_check(got_absent == expected_absent, ctx("absent"))
    end
end

@testset "duplicates: the structure agrees with the scan on which index wins" begin
    rng = Random.Xoshiro(MASTER_SEED + 1_000_000)
    size = 40
    for spec in ALPHABETS, k in KS
        subseed = rand(rng, UInt64)
        srng = Random.Xoshiro(subseed)
        atoms = atom_pool(spec, size, srng)
        # atoms[2] is duplicated into position 30 (the middle); atoms[5]
        # is duplicated into position `size` (the very end). The scan must
        # keep returning 2 and 5, never the later duplicate, and the
        # structure must agree.
        atoms[30] = atoms[2]
        atoms[size] = atoms[5]
        index = route_build(atoms, k)
        ctx(label) = "alphabet=$(spec.name) k=$k subseed=$subseed ($label)"

        for (pos, query) in ((2, atoms[2]), (5, atoms[5]))
            expected = route_scan(atoms, query)
            report_and_check(expected == pos, ctx("scan-sanity pos=$pos"))
            got = route_lookup(index, atoms, query)
            report_and_check(got == expected, ctx("pos=$pos"))
        end
    end
end
