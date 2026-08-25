# Whether `pattern_key_uint64`'s hand-written incremental fold behaves
# like a good hash or clusters. At 389 real atoms and the fold's real
# 64-bit width, no collision will ever be observed (the birthday count is
# astronomically small), so this file narrows the fold artificially, via
# `pattern_key_uint64`'s own `bits` keyword (a test-only knob the shipped
# function already exposes for exactly this, per `pattern_key_reps.jl`'s
# header), and counts real collisions at each width against two
# independent baselines: Julia's own `hash` over the identical key, masked
# the same way, and the birthday approximation for n items into 2^b
# buckets (`n^2 / 2^(b+1)`). This file calls the shipped
# `pattern_key_uint64` directly, never a reimplementation of the fold, so
# a change to the real recipe is what gets measured here.
#
# k=4 is used throughout, not swept: it is the headline recommendation
# (README.md's opening paragraph, 0.812ns at Birkhoff n=60), and this file
# asks whether *that* fold clusters, not how fold quality varies with k.
# n=60 (the larger real Birkhoff size) is used for the same reason.
#
# A subtlety that would otherwise contaminate the count: at k=4, n=60,
# two distinct real atoms occasionally share the same *real* pattern
# (`rowval[1:4]`) by genuine coincidence, a tie in the underlying data,
# not a hash collision. Atoms are deduplicated on their real pattern
# before hashing, so every masked collision counted below is a real fold
# (or hash) collision between two genuinely different patterns, never a
# restatement of a tie already present in the data.
#
# "Many random atom sets": five independent trials, each drawing 8,000
# fresh Birkhoff atoms, pooled by summing collisions and the birthday
# prediction across trials rather than using one giant pool, so the
# result does not depend on one particular draw. Generation is cheap
# (measured at roughly 1.8s for 8,000 atoms at n=60 on this machine); the
# whole file runs in well under a minute.
#
# A finding, not just a check: Julia's `hash(::Vector{Int})` and this
# fold turn out to be *provably* related for a fixed key length, not just
# empirically similar. `hash(pattern) - fold(pattern)` was found, while
# writing this file, to be an exact constant (mod 2^64) across thousands
# of distinct real k=4 patterns. Masking to the low b bits preserves that
# relationship (AND distributes over addition mod a power of two), so the
# two are collision-*isomorphic* at every width tested here: comparing
# this fold against `hash` is a genuine architectural sanity check
# (confirming Julia's own hash is not being reimplemented worse), but it
# is not an independent randomness test the way the birthday comparison
# is, and the printed numbers bear this out (fold and hash pair counts
# come out identical, or nearly so, every run, not just approximately).
# TESTING.md records this as a finding, not a weakness papered over: the
# birthday comparison is the one doing the real discriminating here.
#
# Run: julia --project=. microbenchmark/test_fold_quality.jl

using Test, Random, FrankWolfe

include(joinpath(@__DIR__, "sparse_pattern.jl"))
include(joinpath(@__DIR__, "pattern_key_reps.jl"))
using .SparsePatternLookup, .PatternKeyReps

const MASTER_SEED = 7 # distinct from run.jl(1)..run_pattern_key_reps.jl(4), test_equivalence.jl(5), test_lifecycle.jl(6)
const N_TRIALS = 5
const ATOMS_PER_TRIAL = 8_000
const BIRKHOFF_N = 60
const K = 4
const BIT_WIDTHS = [8, 12, 16, 20, 24]

# Generous on purpose: this is a "does it cluster" gate, not a precise
# statistical test, and the low bit widths' expected counts are Poisson
# enough at these trial sizes that a tight bound would be flaky. A
# genuinely bad (clustering) fold would blow past this by a large margin,
# not graze it; that is the failure mode this is written to catch.
const FACTOR = 4.0
const SLACK = 20.0

within_small_factor(observed, expected) = observed <= FACTOR * expected + SLACK

pair_collisions(sizes) = sum(s * (s - 1) ÷ 2 for s in sizes; init=0)

function masked_bucket_sizes(keys, mask::UInt64)
    counts = Dict{UInt64,Int}()
    for key in keys
        m = key & mask
        counts[m] = get(counts, m, 0) + 1
    end
    return values(counts)
end

# One trial's worth of distinct real atoms: `ATOMS_PER_TRIAL` fresh
# Birkhoff atoms at dimension `n`, deduplicated on their real (unmasked,
# un-narrowed) `rowval[1:k]` pattern so a masked collision below is never
# a restatement of a real tie already present in the data (see this
# file's header).
function distinct_atoms(rng, n::Int, k::Int, count::Int)
    lmo = FrankWolfe.BirkhoffPolytopeLMO()
    atoms = [FrankWolfe.compute_extreme_point(lmo, randn(rng, n, n)) for _ in 1:count]
    seen = Set{Vector{Int}}()
    kept = similar(atoms, 0)
    for a in atoms
        p = pattern_key(a, k)
        if p in seen
            continue
        end
        push!(seen, p)
        push!(kept, a)
    end
    return kept
end

@testset "the UInt64 fold does not cluster, at every narrowed width" begin
    rng = Random.Xoshiro(MASTER_SEED)
    fold_totals = Dict(b => 0 for b in BIT_WIDTHS)
    hash_totals = Dict(b => 0 for b in BIT_WIDTHS)
    birthday_totals = Dict(b => 0.0 for b in BIT_WIDTHS)
    n_atoms_total = 0

    for trial in 1:N_TRIALS
        subseed = rand(rng, UInt64)
        srng = Random.Xoshiro(subseed)
        atoms = distinct_atoms(srng, BIRKHOFF_N, K, ATOMS_PER_TRIAL)
        n = length(atoms)
        n_atoms_total += n

        # Julia's own hash, over the identical key `pattern_key_uint64`
        # folds, computed once at full 64 bits and masked per width below,
        # exactly the same masking the fold's own `bits` keyword applies
        # internally.
        hash_keys_full = [hash(pattern_key(a, K)) for a in atoms]

        for b in BIT_WIDTHS
            mask = (UInt64(1) << b) - UInt64(1)
            fold_keys = [pattern_key_uint64(a, K; bits=b) for a in atoms]
            fold_totals[b] += pair_collisions(masked_bucket_sizes(fold_keys, mask))
            hash_totals[b] += pair_collisions(masked_bucket_sizes(hash_keys_full, mask))
            birthday_totals[b] += n^2 / 2.0^(b + 1)
        end
    end

    println("Fold-quality comparison (master_seed=$MASTER_SEED, $N_TRIALS trials, $n_atoms_total distinct patterns total, k=$K, n=$BIRKHOFF_N):")
    println(rpad("bits", 6), rpad("fold_pairs", 12), rpad("hash_pairs", 12), rpad("birthday", 12), "fold/birthday")
    for b in BIT_WIDTHS
        ratio = birthday_totals[b] > 0 ? round(fold_totals[b] / birthday_totals[b], digits=3) : NaN
        println(rpad(b, 6), rpad(fold_totals[b], 12), rpad(hash_totals[b], 12), rpad(round(birthday_totals[b], digits=2), 12), ratio)

        cond_vs_hash = within_small_factor(fold_totals[b], hash_totals[b])
        cond_vs_hash || println(stderr, "FOLD QUALITY FAILURE (master_seed=$MASTER_SEED): bits=$b fold_pairs=$(fold_totals[b]) is not within a small factor of hash_pairs=$(hash_totals[b])")
        @test cond_vs_hash

        cond_vs_birthday = within_small_factor(fold_totals[b], birthday_totals[b])
        cond_vs_birthday || println(stderr, "FOLD QUALITY FAILURE (master_seed=$MASTER_SEED): bits=$b fold_pairs=$(fold_totals[b]) is not within a small factor of the birthday prediction=$(round(birthday_totals[b], digits=2))")
        @test cond_vs_birthday
    end
end
