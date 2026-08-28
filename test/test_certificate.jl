# `src/certificate.jl`, checked against a plain scan the way
# `test_public_api.jl` checks the index: same alphabets, same ground truth,
# plus the two things the certificate rests on that an index never needed.
# First, that `dot` is a function of its inputs' values: the same gradient
# against two equal atoms, in either argument order, in fresh or misaligned
# memory, returns the same Float64. Second, that a tie between distinct
# atoms, which is the only case the certificate cannot decide, is handed to
# the fall-back and never guessed.
#
# Seeded with `MASTER_SEED = 10`, distinct from every seed already in use
# (run.jl 1; run_prefix/run_lifecycle/run_pattern_key_reps 4;
# test_equivalence 5; test_lifecycle 6; test_fold_quality 7; test_dispatch 8;
# test_public_api 9).
#
# Run, with the other suites: julia --project=. -e 'using Pkg; Pkg.test()'

using Test, Random, LinearAlgebra, FrankWolfe, SparseArrays

using ActiveSetLookup

const MASTER_SEED = 10
const BIRKHOFF_N = 20
const DENSE_DIM = 32

birkhoff_atom(rng) = FrankWolfe.compute_extreme_point(FrankWolfe.BirkhoffPolytopeLMO(), randn(rng, BIRKHOFF_N, BIRKHOFF_N))
linf_atom(rng) = FrankWolfe.compute_extreme_point(FrankWolfe.LpNormBallLMO{Inf}(1.0), randn(rng, DENSE_DIM))
generic_atom(rng) = rand(rng, DENSE_DIM)
gradient_like(rng, atom) = randn(rng, size(atom)...)

# Ground truth, independent of src/: first position whose atom == query.
scan(atoms, query) = something(findfirst(a -> a == query, atoms), -1)

# What a caller holds: every dot(g, a), and the first index of the minimum,
# the way `active_set_argminmax`'s strict `<` picks it.
function caller_state(atoms, g)
    values = [dot(g, a) for a in atoms]
    best = isempty(values) ? -1 : argmin(values)
    return values, best, (best == -1 ? NaN : values[best])
end

function report_and_check(cond::Bool, context::String)
    cond || println(stderr, "CERTIFICATE FAILURE (master_seed=$MASTER_SEED): ", context)
    @test cond
end

# Runs every query through both certified_lookup methods and the scan, with a
# fall-back that records whether it was reached.
function agree_with_scan(atoms, g, query, ctx::String)
    values, best, best_value = caller_state(atoms, g)
    query_value = dot(g, query)
    expected = scan(atoms, query)
    reached = Ref(false)
    fb = (a, q) -> (reached[] = true; scan(a, q))
    got = certified_lookup(atoms, query, query_value, best, best_value; fallback=fb)
    report_and_check(got == expected, "$ctx argmin form: got $got expected $expected")
    got_values = certified_lookup(atoms, query, query_value, values)
    report_and_check(got_values == expected, "$ctx values form: got $got_values expected $expected")
    # Two things the fall-back must never do: run after the certificate has
    # decided, or run when the best atom was the query. (An arbitrary query
    # that scores above the minimum does reach it, correctly: only an LMO
    # vertex is promised to score at or below every active atom, and the
    # testset below checks that promise.)
    if !isempty(atoms)
        certified = certified_absent(query_value, best_value)
        report_and_check(!(certified && reached[]), "$ctx fallback ran after the certificate decided")
        report_and_check(!(atoms[best] == query && reached[]), "$ctx fallback ran although the best atom was the query")
    end
    return nothing
end

@testset "dot is a function of values: equal inputs, either order, any alignment" begin
    rng = Random.Xoshiro(MASTER_SEED)
    for d in (17, 3000), _ in 1:200
        g = randn(rng, d)
        a = sign.(randn(rng, d))
        buffer = zeros(d + 7)
        buffer[4:d+3] .= a
        @test dot(g, a) == dot(g, copy(a))
        @test dot(g, a) == dot(a, g)
        @test dot(g, a) == dot(g, view(buffer, 4:d+3))
        @test dot(g, a) == dot(g, buffer[4:d+3])
    end
    for n in (20, 60), _ in 1:200
        g = randn(rng, n, n)
        a = FrankWolfe.compute_extreme_point(FrankWolfe.BirkhoffPolytopeLMO(), g)
        @test dot(g, a) == dot(g, copy(a))
        @test dot(g, a) == dot(a, g)
    end
end

@testset "agrees with a scan: three alphabets, present, absent, random gradients" begin
    rng = Random.Xoshiro(MASTER_SEED + 1)
    generators = (birkhoff=birkhoff_atom, linf=linf_atom, generic=generic_atom)
    for (name, gen) in pairs(generators), size in (0, 1, 2, 5, 20, 60)
        subseed = rand(rng, UInt64)
        srng = Random.Xoshiro(subseed)
        atoms = eltype([gen(srng)])[]
        while length(atoms) < size                       # no duplicates: find_atom's own invariant
            candidate = gen(srng)
            scan(atoms, candidate) == -1 && push!(atoms, candidate)
        end
        g = gradient_like(srng, gen(srng))
        ctx(label) = "alphabet=$name size=$size subseed=$subseed ($label)"
        for pos in unique((1, size, rand(srng, 1:max(size, 1), 3)...))
            1 <= pos <= size || continue
            agree_with_scan(atoms, g, atoms[pos], ctx("present pos=$pos"))
        end
        absent = gen(srng)
        while scan(atoms, absent) != -1
            absent = gen(srng)
        end
        agree_with_scan(atoms, g, absent, ctx("absent"))
    end
end

@testset "an LMO vertex never reaches the fall-back, except on a tie" begin
    # The situation the certificate is built for: the query is the LMO's
    # answer for the very gradient the active set was just minimised over,
    # so it scores at or below every active atom. Absent, it scores strictly
    # below (a continuous random gradient makes ties measure-zero) and the
    # certificate alone decides; present, it is the best atom and one
    # comparison decides. Either way the search is never run.
    rng = Random.Xoshiro(MASTER_SEED + 5)
    lmos = (
        birkhoff=(FrankWolfe.BirkhoffPolytopeLMO(), () -> randn(rng, BIRKHOFF_N, BIRKHOFF_N)),
        linf=(FrankWolfe.LpNormBallLMO{Inf}(1.0), () -> randn(rng, DENSE_DIM)),
    )
    for (name, (lmo, direction)) in pairs(lmos), size in (1, 5, 20, 60)
        atoms = [FrankWolfe.compute_extreme_point(lmo, direction()) for _ in 1:size]
        for trial in 1:40
            g = direction()
            query = FrankWolfe.compute_extreme_point(lmo, g)
            values, best, best_value = caller_state(atoms, g)
            query_value = dot(g, query)
            expected = scan(atoms, query)
            reached = Ref(false)
            fb = (a, q) -> (reached[] = true; scan(a, q))
            got = certified_lookup(atoms, query, query_value, best, best_value; fallback=fb)
            ctx = "alphabet=$name size=$size trial=$trial expected=$expected"
            report_and_check(got == expected, "$ctx got=$got")
            report_and_check(!reached[], "$ctx fallback reached for an LMO vertex")
            report_and_check(certified_lookup(atoms, query, query_value, values) == expected, "$ctx values form")
            # And the shape of the decision: absent means certified, present means the best atom.
            if expected == -1
                report_and_check(certified_absent(query_value, best_value), "$ctx absent but not certified")
            else
                report_and_check(best == expected, "$ctx present but best=$best")
            end
        end
    end
end

@testset "ties between distinct atoms go to the fall-back, and the answer stays right" begin
    # A gradient with few distinct values makes many distinct atoms score the
    # same, which is the one case the certificate cannot decide. A constant
    # gradient is the extreme: every permutation matrix scores identically.
    rng = Random.Xoshiro(MASTER_SEED + 2)
    for (name, gen) in pairs((birkhoff=birkhoff_atom, linf=linf_atom))
        atoms = eltype([gen(rng)])[]
        while length(atoms) < 30
            candidate = gen(rng)
            scan(atoms, candidate) == -1 && push!(atoms, candidate)
        end
        for (glabel, g) in (
            ("constant", fill(1.0, size(atoms[1])...)),
            ("few values", Float64.(rand(rng, -1:1, size(atoms[1])...))),
            ("zero", zeros(size(atoms[1])...)),
        )
            for pos in (1, 7, 30)
                agree_with_scan(atoms, g, atoms[pos], "alphabet=$name gradient=$glabel present pos=$pos")
            end
            absent = gen(rng)
            while scan(atoms, absent) != -1
                absent = gen(rng)
            end
            agree_with_scan(atoms, g, absent, "alphabet=$name gradient=$glabel absent")
        end
    end
end

@testset "duplicates resolve to the first copy, as find_atom does" begin
    rng = Random.Xoshiro(MASTER_SEED + 3)
    base = [linf_atom(rng) for _ in 1:6]
    atoms = [base[1], base[2], base[3], base[2], base[4], base[3], base[3]]
    g = randn(rng, DENSE_DIM)
    for pos in eachindex(atoms)
        agree_with_scan(atoms, g, atoms[pos], "duplicates pos=$pos")
    end
end

@testset "the certificate never fires on an empty set or a broken gradient" begin
    rng = Random.Xoshiro(MASTER_SEED + 4)
    atoms = [linf_atom(rng) for _ in 1:5]
    @test certified_lookup(eltype(atoms)[], atoms[1], -1.0, -1, NaN) == -1
    @test certified_lookup(eltype(atoms)[], atoms[1], -1.0, Float64[]) == -1
    for g in (fill(NaN, DENSE_DIM), fill(Inf, DENSE_DIM), [Inf; zeros(DENSE_DIM - 1)])
        for pos in (1, 5)
            agree_with_scan(atoms, g, atoms[pos], "gradient=$(g[1]) present pos=$pos")
        end
        absent = linf_atom(rng)
        while scan(atoms, absent) != -1
            absent = linf_atom(rng)
        end
        agree_with_scan(atoms, g, absent, "gradient=$(g[1]) absent")
    end
end

@testset "no signed-zero gap: the certificate keys on nothing" begin
    # The dense value key needed `+ 0.0` to stop -0.0 and 0.0 from splitting
    # one atom into two (METHOD.md). The certificate compares dot products
    # and then atoms with `==`, and `-0.0 == 0.0`, so there is nothing to
    # canonicalise: the two atoms are one atom at every step.
    stored = [0.0, 1.0, 2.0]
    query = [-0.0, 1.0, 2.0]
    g = [3.0, -1.0, 0.5]
    atoms = [stored]
    values, best, best_value = caller_state(atoms, g)
    @test certified_lookup(atoms, query, dot(g, query), best, best_value) == 1
    @test certified_lookup(atoms, query, dot(g, query), values) == 1
end

@testset "certified_absent is the guard, and only the guard" begin
    @test certified_absent(-2.0, -1.0)
    @test !certified_absent(-1.0, -1.0)
    @test !certified_absent(-0.5, -1.0)
    @test !certified_absent(NaN, -1.0)
    @test scan_atoms([[1.0, 2.0], [3.0, 4.0]], [3.0, 4.0]) == 2
    @test scan_atoms([[1.0, 2.0], [3.0, 4.0]], [9.0, 9.0]) == -1
end
