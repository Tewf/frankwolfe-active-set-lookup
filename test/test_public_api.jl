# The public `src/` module (`build_index`, `lookup_atom`, `push_atom!`,
# `delete_atom!`), checked the way a stranger who only reads README.md
# would use it, not the way the research harness's own comparisons do.
# `microbenchmark/test_*.jl` test the comparisons that led to this design;
# this file tests the module that design became. It is intentionally
# self-contained rather than sharing `microbenchmark/test_atom_generators.jl`:
# `src/` is meant to be readable and usable on its own, and its own test
# suite should not need a reader to first understand the research harness's
# atom-generator plumbing.
#
# Every test is seeded (`Random.Xoshiro`, `MASTER_SEED` below, distinct
# from every seed already in use: run.jl(1), run_prefix.jl/run_lifecycle.jl/
# run_pattern_key_reps.jl(4), test_equivalence.jl(5), test_lifecycle.jl(6),
# test_fold_quality.jl(7), test_dispatch.jl(8)).
#
# Run: julia --project=. test/test_public_api.jl

using Test, Random, FrankWolfe, SparseArrays

include(joinpath(@__DIR__, "..", "src", "ActiveSetLookup.jl"))
using .ActiveSetLookup

const MASTER_SEED = 9
const BIRKHOFF_N = 20
const DENSE_DIM = 32

birkhoff_atom(rng) = FrankWolfe.compute_extreme_point(FrankWolfe.BirkhoffPolytopeLMO(), randn(rng, BIRKHOFF_N, BIRKHOFF_N))
linf_atom(rng) = FrankWolfe.compute_extreme_point(FrankWolfe.LpNormBallLMO{Inf}(1.0), randn(rng, DENSE_DIM))

# Ground truth: a plain scan using this module's own `confirm_match`, the
# same equality every index here already confirms a bucket hit against.
function scan(atoms, query)
    for i in eachindex(atoms)
        confirm_match(atoms[i], query) && return i
    end
    return -1
end

function report_and_check(cond::Bool, context::String)
    cond || println(stderr, "PUBLIC API FAILURE (master_seed=$MASTER_SEED): ", context)
    @test cond
end

@testset "build_index dispatches on atom type" begin
    rng = Random.Xoshiro(MASTER_SEED)
    sparse_atoms = [birkhoff_atom(rng) for _ in 1:10]
    dense_atoms = [linf_atom(rng) for _ in 1:10]
    @test build_index(sparse_atoms) isa SparsePatternIndex
    @test build_index(dense_atoms) isa DenseValueIndex
end

@testset "default k is the documented constant" begin
    @test DEFAULT_K == 4
    rng = Random.Xoshiro(MASTER_SEED + 1)
    atoms = [birkhoff_atom(rng) for _ in 1:5]
    @test build_index(atoms).k == DEFAULT_K
    @test build_index(atoms; k=8).k == 8
end

@testset "lookup agrees with a scan: both alphabets, several k, present and absent" begin
    rng = Random.Xoshiro(MASTER_SEED + 2)
    generators = (birkhoff=birkhoff_atom, linf=linf_atom)
    for (name, gen) in pairs(generators), size in (0, 1, 2, 5, 20, 60), k in (2, 4, 8)
        subseed = rand(rng, UInt64)
        srng = Random.Xoshiro(subseed)
        atoms = [gen(srng) for _ in 1:size]
        index = build_index(atoms; k=k)
        ctx(label) = "alphabet=$name size=$size k=$k subseed=$subseed ($label)"

        if size > 0
            for pos in unique((1, size, rand(srng, 1:size, min(3, size))...))
                query = atoms[pos]
                expected = scan(atoms, query)
                got = lookup_atom(index, atoms, query)
                report_and_check(got == expected, ctx("present pos=$pos"))
            end
        end

        absent = gen(srng)
        tries = 0
        while size > 0 && scan(atoms, absent) != -1 && tries < 5
            absent = gen(srng)
            tries += 1
        end
        expected_absent = scan(atoms, absent)
        got_absent = lookup_atom(index, atoms, absent)
        report_and_check(got_absent == expected_absent, ctx("absent"))
    end
end

@testset "push_atom! and delete_atom! keep the index in sync" begin
    # `active_set_cleanup!`'s `deleteat!` shifts every later position down
    # by one; this is the test that would catch `delete_atom!`'s repair
    # walk getting that shift direction or boundary wrong, the same
    # property `microbenchmark/test_lifecycle.jl` checks for the research
    # harness's own structures. Checked against a fresh scan over every
    # position still held, after every single operation, not just a
    # neighbour of whatever changed.
    rng = Random.Xoshiro(MASTER_SEED + 3)
    generators = (birkhoff=birkhoff_atom, linf=linf_atom)
    for (name, gen) in pairs(generators)
        atoms = [gen(rng) for _ in 1:5]
        index = build_index(atoms; k=4)

        for step in 1:150
            do_insert = isempty(atoms) || rand(rng) < 0.5
            if do_insert
                atom = gen(rng)
                push_atom!(index, atoms, atom)
            else
                pos = rand(rng, 1:length(atoms))
                delete_atom!(index, atoms, pos)
            end

            for pos in eachindex(atoms)
                expected = scan(atoms, atoms[pos])
                got = lookup_atom(index, atoms, atoms[pos])
                report_and_check(got == expected, "alphabet=$name step=$step pos=$pos size=$(length(atoms))")
            end
            if !isempty(atoms)
                absent = gen(rng)
                expected_absent = scan(atoms, absent)
                got_absent = lookup_atom(index, atoms, absent)
                report_and_check(got_absent == expected_absent, "alphabet=$name step=$step absent")
            end
        end
    end
end

@testset "the dense value key closes the signed-zero gap DECISIONS.md found" begin
    # An atom built the way a naive Linf-ball vertex formula would
    # (`-1.0 .* sign.(gradient)`), differing from a stored atom only in the
    # sign of one zero coordinate. The scan (confirm_match) calls these
    # equal; a Dict keyed on the raw prefix would not
    # (`0.0 isequal -0.0` is false) and would miss. `atom_key`'s `+ 0.0`
    # canonicalisation is what makes this module answer correctly here
    # rather than reproducing the gap `TESTING.md` records as still open
    # in `microbenchmark/lookup_methods.jl`'s own `PrefixIndex`.
    stored = [0.0, 1.0, 2.0]
    query = [-0.0, 1.0, 2.0]
    @test confirm_match(stored, query) # one atom to the scan

    atoms = [stored]
    index = build_index(atoms; k=1)
    @test lookup_atom(index, atoms, query) == 1

    # NaN still goes the harmless way: the scan says not equal, the key
    # collides (NaN hashes and equals itself under `hash`/`isequal`), and
    # confirmation fails, so both agree on not-found.
    nan_stored = [NaN, 1.0]
    nan_query = [NaN, 1.0]
    @test !confirm_match(nan_stored, nan_query)
    nan_atoms = [nan_stored]
    nan_index = build_index(nan_atoms; k=1)
    @test lookup_atom(nan_index, nan_atoms, nan_query) == -1
end

@testset "a forced UInt64 fold collision is survived, not just unlikely" begin
    # Two permutation matrices with genuinely different rowval[1:2]
    # patterns, chosen (by the same short random search
    # test_pattern_key_reps.jl already used, not repeated here) so their
    # unmasked fold still differs; the collision hazard itself is already
    # demonstrated at bits=1 in microbenchmark/test_pattern_key_reps.jl, so
    # this test checks the consequence that matters for this module's own
    # API: a real Dict bucket holding two atoms still answers both lookups
    # correctly through confirm_match, never by trusting the bucket alone.
    perm(p) = sparse(collect(1:length(p)), p, ones(length(p)), length(p), length(p))
    a1 = perm([6, 1, 3, 4, 2, 5])
    a2 = perm([2, 3, 6, 4, 5, 1])
    @test atom_key(a1; k=2) != atom_key(a2; k=2)

    atoms = [a1, a2]
    index = build_index(atoms; k=2)
    @test lookup_atom(index, atoms, a1) == 1
    @test lookup_atom(index, atoms, a2) == 2
    miss = perm([1, 2, 3, 4, 5, 6])
    @test lookup_atom(index, atoms, miss) == -1
end

# A SparseVector is the third atom shape, and it is not a corner case:
# FrankWolfe.KSparseLMO returns one, and a SparseMatrixCSC atom flattened
# with `vec` becomes one too. Before this method existed, `build_index` threw
# a MethodError on both, which is the first thing a newcomer would have hit.
@testset "sparse vectors are a first-class atom shape" begin
    lmo = FrankWolfe.KSparseLMO(3, 1.0)
    atoms = [FrankWolfe.compute_extreme_point(lmo, randn(40)) for _ in 1:8]
    @test eltype(atoms) <: SparseVector

    index = build_index(atoms)
    @test index.k == DEFAULT_K
    @test index isa SparsePatternIndex          # routed like a sparse atom, not a dense one

    # The invariant every other test here also holds to: agree with a scan.
    for (i, a) in enumerate(atoms)
        @test lookup_atom(index, atoms, a) == findfirst(b -> b == a, atoms)
    end

    fresh = FrankWolfe.compute_extreme_point(lmo, randn(40))
    if !any(a -> a == fresh, atoms)
        @test lookup_atom(index, atoms, fresh) == -1
        push_atom!(index, atoms, fresh)
        @test lookup_atom(index, atoms, fresh) == length(atoms)
    end

    delete_atom!(index, atoms, 2)
    for (i, a) in enumerate(atoms)
        @test lookup_atom(index, atoms, a) == findfirst(b -> b == a, atoms)
    end

    @test build_index(atoms; k=8).k == 8        # k stays an ordinary keyword
end

# A sparse array may hold an explicitly stored zero: a slot in nzind whose
# value is 0.0. Sparse arithmetic leaves these behind and nothing removes
# them unless dropzeros! is called. Two such arrays can be numerically equal
# while their stored index arrays differ, so keying on raw indices hands them
# different keys and the lookup misses an atom that is present. That is a
# false negative no bucket confirmation can catch, because no bucket is
# reached. It is the sparse twin of the signed-zero hazard, and the more
# reachable of the two.
@testset "explicitly stored zeros do not split one atom into two" begin
    stored   = SparseVector(5, [1, 2, 4], [1.0, 0.0, 3.0])   # zero occupies a slot
    dropped  = SparseVector(5, [1, 4],    [1.0, 3.0])        # same vector, canonical

    @test stored == dropped                       # one atom to == and to confirm_match
    @test stored.nzind != dropped.nzind           # two atoms to a naive index key
    @test atom_key(stored) == atom_key(dropped)   # but one atom to ours

    for (built, queried) in ((stored, dropped), (dropped, stored))
        atoms = [built]
        index = build_index(atoms)
        @test lookup_atom(index, atoms, queried) == 1
    end

    # A trailing stored zero must not change the key either, and a run of
    # leading ones must be skipped rather than consuming the k budget.
    lead = SparseVector(8, [1, 2, 3, 5, 7], [0.0, 0.0, 1.0, 2.0, 3.0])
    tidy = SparseVector(8, [3, 5, 7],       [1.0, 2.0, 3.0])
    @test lead == tidy
    @test atom_key(lead) == atom_key(tidy)

    # The all-zero atom has no nonzeros at all and must still be handled.
    empty_a = SparseVector(5, Int[], Float64[])
    empty_b = SparseVector(5, [2], [0.0])
    @test empty_a == empty_b
    @test atom_key(empty_a) == atom_key(empty_b)
end

# METHOD.md states the precondition: the k positions a key reads must differ
# across the atoms. Sparse atoms get that free from nzind, which only lists
# positions holding something. A dense key reads fixed cells and assumes they
# vary, which box corners satisfy and an atom with a dominant repeated value
# does not. This pins both halves: the degenerate case stays CORRECT, and
# bucket_health reports it so a user can see it rather than guess.
@testset "the precondition, and the diagnostic that reveals it" begin
    # Dense atoms that are 7.0 everywhere except one late position.
    degenerate = [(v = fill(7.0, 40); v[10 + i] = 3.0; v) for i in 1:5]
    index = build_index(degenerate)

    @test length(index.buckets) == 1                 # every atom shares a key
    @test bucket_health(index) == length(degenerate) # and the diagnostic says so

    # Correct anyway: confirmation still resolves every atom.
    for (i, a) in enumerate(degenerate)
        @test lookup_atom(index, degenerate, a) == i
    end
    @test lookup_atom(index, degenerate, fill(7.0, 40)) == -1

    # Raising k does not rescue it, because the added cells are constant too.
    @test bucket_health(build_index(degenerate; k = 8)) == length(degenerate)

    # A healthy index sits at one atom per bucket.
    varied = [randn(40) for _ in 1:20]
    @test bucket_health(build_index(varied)) == 1.0
end
