# Closes the two gaps ../DECISIONS.md and ../README.md's original answer
# left open (see ../REJECTED.md's section 2 for the write-up):
#
# Gap 1: run.jl's Dict hashes the *whole* atom, always O(dimension). A hash
# over only the atom's first `k` coordinates costs O(k), flat in dimension,
# and is exactly as sound: a bucket hit still has to be confirmed with the
# same exact equality the scan uses, so a shorter hash can only cost speed,
# never correctness (see lookup_methods.jl's `prefix_lookup`).
#
# Gap 2: run.jl only ever times a miss (a fresh, never-stored atom), which
# is the scan's best case: `!=` exits at the first coordinate almost
# surely. A hit forces the scan through the whole matching atom. This file
# times miss, hit, and a mixed workload for every k.
#
# The decisive axis is the atom alphabet, not the mechanism: run.jl's
# "generic" atoms are `rand(dim)` Float64, where one coordinate almost
# always discriminates. `FrankWolfe.jl`'s real atoms are not like that.
# Birkhoff-polytope vertices are permutation matrices (1 bit of
# information per entry, 0 or 1) and L-infinity-ball vertices are box
# corners (also 1 bit per entry, -1 or +1); see ../measurement/problems.jl.
# Both are generated here by calling FrankWolfe.jl's own LMOs
# (`BirkhoffPolytopeLMO`, `LpNormBallLMO{Inf}`), the same functions
# ../measurement/problems.jl calls, not a hand-rolled approximation.
#
# Permutation atoms are kept sparse (`SparseMatrixCSC`), not flattened to
# dense, because `_unsafe_equal`'s sparse branch (active_set.jl:513) is
# `==`, which timing showed costs about 20ns flat on a mismatch and scales
# with the atom's nonzero count on a match, not its dimension (see
# lookup_methods.jl's "Sparse-atom variants" section). Scanning a
# *flattened* permutation matrix with the dense `linear_scan` would
# overstate what `find_atom` really costs on this alphabet, so the scan
# baseline for permutation atoms is the real sparse `==`, not that loop.
#
# Dimensions are the exact dimensions ../measurement/'s three real runs
# used (625 = 25x25 Birkhoff, 3600 = 60x60 Birkhoff, 3000 = the
# L-infinity ball), so a crossover found here is directly comparable to
# that run's own active-set sizes (158, 389, 241).
include("lookup_methods.jl")
include("timing.jl")
using .LookupMethods
using .Timing
using FrankWolfe, Random, SparseArrays

Random.seed!(2) # distinct from run.jl's seed!(1): this is a separate sweep

const SIZES = [1, 2, 5, 10, 20, 50, 100, 158, 200, 300, 389, 500, 750, 1000, 1500, 2000]
const MAX_SIZE = maximum(SIZES)
const K_BASE = [1, 2, 4, 8, 16, 32, 64]
const MIX_HIT_PROBABILITY = 0.5 # a stress mix; ../measurement/results.csv's
# three real BPCG runs saw 0 hits out of 159/389/240 calls (0%), so a
# 50/50 mix is deliberately harder on the scan than anything measured, to
# show how the answer would move if hits became common. See ../README.md.
const QUERY_N = 200 # pre-generated queries every timing cycles through

k_values(dim) = sort(unique(vcat(filter(<=(dim), K_BASE), [dim])))

# Every query type is timed the same way: QUERY_N pre-generated queries,
# cycled through by `time_per_call_seq` (see timing.jl). This matters most
# for :miss and :hit on a low-entropy alphabet (box corners, permutation
# matrices): a *single* fixed query, repeated, only samples one point in a
# heavily skewed bucket-size distribution, so its timing reflects that one
# query's luck, not the method's real average cost. QUERY_N distinct draws
# average that out, the same way ../microbenchmark/run.jl's own "fastest of
# five" averages out timer noise, not sampling noise.
function make_queries(sc, atoms_slice, size, querytype)
    if querytype == :miss
        return [sc.fresh() for _ in 1:QUERY_N]
    elseif querytype == :hit
        return [copy(atoms_slice[rand(1:size)]) for _ in 1:QUERY_N]
    else # :mix
        return [
            rand() < MIX_HIT_PROBABILITY ? copy(atoms_slice[rand(1:size)]) : sc.fresh() for
            _ in 1:QUERY_N
        ]
    end
end

generic_pool(size, dim) = [rand(dim) for _ in 1:size]
generic_fresh(dim) = () -> rand(dim)

# Permutation matrices are kept sparse here, not flattened to a dense
# vector, so the scan baseline below is FrankWolfe.jl's real comparison
# cost, not an inflated one: see run_prefix.jl's header and
# lookup_methods.jl's "Sparse-atom variants" section for why. The atoms
# are still exactly what FrankWolfe.BirkhoffPolytopeLMO.compute_extreme_point
# produces (see ../measurement/problems.jl), and sparse_prefix reads a
# prefix straight off that structure without densifying it first.
function permutation_pool(size, n)
    lmo = FrankWolfe.BirkhoffPolytopeLMO()
    return [FrankWolfe.compute_extreme_point(lmo, randn(n, n)) for _ in 1:size]
end
permutation_fresh(n) = () -> FrankWolfe.compute_extreme_point(FrankWolfe.BirkhoffPolytopeLMO(), randn(n, n))

# Box corners exactly as FrankWolfe.LpNormBallLMO{Inf}.compute_extreme_point
# produces them (see ../measurement/problems.jl's linf_box_problem):
# independent +-1.0 entries, one bit of information per coordinate.
function box_pool(size, d)
    lmo = FrankWolfe.LpNormBallLMO{Inf}(1.0)
    return [FrankWolfe.compute_extreme_point(lmo, randn(d)) for _ in 1:size]
end
box_fresh(d) = () -> FrankWolfe.compute_extreme_point(FrankWolfe.LpNormBallLMO{Inf}(1.0), randn(d))

scenarios = [
    (alphabet=:generic, dim=625, sparse=false, atoms=generic_pool(MAX_SIZE, 625), fresh=generic_fresh(625)),
    (alphabet=:generic, dim=3600, sparse=false, atoms=generic_pool(MAX_SIZE, 3600), fresh=generic_fresh(3600)),
    (alphabet=:generic, dim=3000, sparse=false, atoms=generic_pool(MAX_SIZE, 3000), fresh=generic_fresh(3000)),
    (
        alphabet=:permutation,
        dim=625,
        sparse=true,
        atoms=permutation_pool(MAX_SIZE, 25),
        fresh=permutation_fresh(25),
    ),
    (
        alphabet=:permutation,
        dim=3600,
        sparse=true,
        atoms=permutation_pool(MAX_SIZE, 60),
        fresh=permutation_fresh(60),
    ),
    (alphabet=:box, dim=625, sparse=false, atoms=box_pool(MAX_SIZE, 625), fresh=box_fresh(625)),
    (alphabet=:box, dim=3600, sparse=false, atoms=box_pool(MAX_SIZE, 3600), fresh=box_fresh(3600)),
    (alphabet=:box, dim=3000, sparse=false, atoms=box_pool(MAX_SIZE, 3000), fresh=box_fresh(3000)),
]

timing_rows = []
collision_rows = []

for sc in scenarios
    ks = k_values(sc.dim)
    build_index = sc.sparse ? build_sparse_prefix_index : build_prefix_index
    scan = sc.sparse ? sparse_linear_scan : linear_scan
    lookup = sc.sparse ? sparse_prefix_lookup : prefix_lookup

    for size in SIZES
        atoms_slice = sc.atoms[1:size]

        # Built once per (scenario, size, k), reused below for collision
        # stats and every query type: index construction is not timed, the
        # same convention run.jl's Dict construction follows.
        indices = Dict(k => build_index(atoms_slice, k) for k in ks)

        for k in ks
            stats = collision_stats(indices[k])
            push!(collision_rows, (alphabet=sc.alphabet, dim=sc.dim, size=size, stats...))
        end

        # Index construction just above allocates heavily (up to `size`
        # bucket vectors per k). Left uncollected, that garbage would sit
        # between the group's measurements and could get swept by the
        # allocator mid-timing, landing on whichever measurement happens to
        # run right after the burst rather than on the method actually
        # responsible for it. A `GC.gc()` here, and again before each
        # timing call below, gives every measurement in this group the same
        # clean starting heap instead of favouring whichever ran later.
        GC.gc()

        for querytype in (:miss, :hit, :mix)
            queries = make_queries(sc, atoms_slice, size, querytype)
            GC.gc()
            scan_ns = time_per_call_seq(scan, (atoms_slice,), queries) * 1e9
            push!(
                timing_rows,
                (
                    alphabet=sc.alphabet,
                    dim=sc.dim,
                    size=size,
                    querytype=querytype,
                    method="scan",
                    k=0,
                    ns=scan_ns,
                ),
            )

            for k in ks
                prefix_ns = time_per_call_seq(lookup, (indices[k], atoms_slice), queries) * 1e9
                push!(
                    timing_rows,
                    (
                        alphabet=sc.alphabet,
                        dim=sc.dim,
                        size=size,
                        querytype=querytype,
                        method="prefix",
                        k=k,
                        ns=prefix_ns,
                    ),
                )
            end
        end
        println("$(sc.alphabet) dim=$(sc.dim) size=$size done")
    end
end

open(joinpath(@__DIR__, "results_prefix_timing.csv"), "w") do io
    println(io, "alphabet,dim,size,querytype,method,k,ns")
    for r in timing_rows
        println(io, "$(r.alphabet),$(r.dim),$(r.size),$(r.querytype),$(r.method),$(r.k),$(round(r.ns,digits=1))")
    end
end

open(joinpath(@__DIR__, "results_prefix_collisions.csv"), "w") do io
    println(
        io,
        "alphabet,dim,size,k,n_atoms,n_buckets,mean_bucket_size,max_bucket_size,atom_collision_rate,bucket_collision_rate",
    )
    for r in collision_rows
        println(
            io,
            "$(r.alphabet),$(r.dim),$(r.size),$(r.k),$(r.n_atoms),$(r.n_buckets),$(round(r.mean_bucket_size,digits=2)),$(r.max_bucket_size),$(round(r.atom_collision_rate,sigdigits=4)),$(round(r.bucket_collision_rate,sigdigits=4))",
        )
    end
end

# Crossover: for each (alphabet, dim, querytype, k), the smallest tested
# size at which the prefix-hash timing beats that same group's scan
# timing, or "none" if it never does within SIZES.
crossover_rows = []
for sc in scenarios, querytype in (:miss, :hit, :mix)
    scan_by_size = Dict(
        r.size => r.ns for
        r in timing_rows if r.alphabet == sc.alphabet && r.dim == sc.dim && r.querytype == querytype && r.method == "scan"
    )
    for k in k_values(sc.dim)
        prefix_by_size = sort(
            [
                r for r in timing_rows if
                r.alphabet == sc.alphabet && r.dim == sc.dim && r.querytype == querytype && r.method == "prefix" && r.k == k
            ];
            by=r -> r.size,
        )
        crossing = findfirst(r -> r.ns < scan_by_size[r.size], prefix_by_size)
        crossover_size = crossing === nothing ? nothing : prefix_by_size[crossing].size
        push!(
            crossover_rows,
            (alphabet=sc.alphabet, dim=sc.dim, querytype=querytype, k=k, crossover_size=crossover_size),
        )
    end
end

open(joinpath(@__DIR__, "results_prefix_crossover.csv"), "w") do io
    println(io, "alphabet,dim,querytype,k,crossover_size")
    for r in crossover_rows
        label = r.crossover_size === nothing ? "none" : string(r.crossover_size)
        println(io, "$(r.alphabet),$(r.dim),$(r.querytype),$(r.k),$label")
    end
end

println()
for r in crossover_rows
    label = r.crossover_size === nothing ? "no crossover in [$(minimum(SIZES)), $(maximum(SIZES))]" : "crossover at size $(r.crossover_size)"
    println("$(r.alphabet) dim=$(r.dim) query=$(r.querytype) k=$(r.k): $label")
end
