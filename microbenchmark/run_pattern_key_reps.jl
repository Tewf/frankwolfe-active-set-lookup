# The representation stage's brief: `sparse_pattern.jl`'s `pattern_key`
# already won on total per-iteration cost (DECISIONS.md's "The sparse-pattern key and the trie"), but its
# key is a `Vector{Int}`, and every lookup and every insert allocates one
# before the `Dict` is even touched. `pattern_key_reps.jl` builds two
# allocation-free representations of the same key (`UInt64`, `NTuple{K,Int}`);
# this file measures whether removing that allocation actually moves the
# total, the way `run_lifecycle.jl` measured whether the sparse-pattern key
# itself moved the total over the flat value-prefix hash.
#
# Allocation is the point of this file (`@allocated` bytes, `@allocations`
# count), measured directly on one warmed-up call, not batched the way
# sub-microsecond *time* has to be (`../MEASURING.md`): `@allocated` has no
# timer-resolution floor to clear, so a single clean call after warm-up is
# the exact number, not an estimate.
#
# Two allocation-measurement pitfalls, found while writing this file, that
# a careless version of it would have reported as the *representations'*
# cost when they were really the harness's own:
#
#   1. A closure stored in a struct field typed `::Function` (an abstract
#      type) loses its concrete type there, so calling it through that
#      field forces a dynamic dispatch, which allocates, on every call:
#      `pattern_lookup_u64` measured 0 bytes called directly, 80 bytes
#      called through a `struct Rep; lookup::Function; end` field, with
#      nothing about the lookup itself different. Every representation
#      below is passed as a closure argument straight into a `where`
#      function, never stored in an abstractly-typed field, so it keeps
#      its own concrete type end to end.
#   2. `@allocated f(make_args()...)` counts `make_args()`'s own
#      allocation together with `f`'s: building a fresh copy of a bucket
#      Dict to hand a mutating function (so it doesn't corrupt state a
#      later measurement needs) has to happen *before* the `@allocated`
#      expression starts, or the copy's cost is charged to the call it
#      was only there to protect.
#
# Time is measured the same way `run_lifecycle.jl` already does (lookup:
# `time_per_call_seq` over a fixed query sequence; insert: marginal cost of
# building at `size` vs. `size+INSERT_BATCH`, differenced, never a
# structure repeatedly mutated in place; delete-repair: `time_per_call` on
# one built index's buckets, the same accepted quirk `bucket_lifecycle.jl`'s
# header documents, i.e. repeated calls keep walking the same amount of
# work even though the stored positions drift after the first), so a
# number here is directly comparable to `results_lifecycle_timing.csv`'s.
#
# Scope: the three pattern-key representations only make sense for
# Birkhoff's sparse atoms (`pattern_key_reps.jl` reads `SparseMatrixCSC`'s
# own `rowval`), so they are swept there, at k in {2,4,8} (the brief: "the
# pattern key already gave perfect discrimination at k=4, so include k=2
# and k=4 as well as k=8"). The L∞-ball and the generic control have no
# sparse structure to key on, so they are swept once each with the
# existing `Vector{Float64}` value-prefix hash (`lookup_methods.jl`'s
# `PrefixIndex`), included per the brief ("for... the Vector{Float64}
# value-prefix baseline on the dense alphabet") as the allocation
# comparison point: this structure is not new, but its allocation cost was
# never measured directly before this file, only its time.
include("lookup_methods.jl")
include("sparse_pattern.jl")
include("pattern_key_reps.jl")
include("bucket_lifecycle.jl")
include("timing.jl")
using .LookupMethods, .SparsePatternLookup, .PatternKeyReps, .BucketLifecycle, .Timing
using FrankWolfe, Random, SparseArrays

Random.seed!(4) # distinct from run.jl's seed!(1), run_prefix.jl's seed!(2), run_lifecycle.jl's seed!(3)

const K_GRID = [2, 4, 8]
const QUERY_N = 200
const INSERT_BATCH = 100 # atoms added between the two builds insert time is differenced from

permutation_pool(size, n) = [FrankWolfe.compute_extreme_point(FrankWolfe.BirkhoffPolytopeLMO(), randn(n, n)) for _ in 1:size]
permutation_fresh(n) = () -> FrankWolfe.compute_extreme_point(FrankWolfe.BirkhoffPolytopeLMO(), randn(n, n))
box_pool(size, d) = [FrankWolfe.compute_extreme_point(FrankWolfe.LpNormBallLMO{Inf}(1.0), randn(d)) for _ in 1:size]
box_fresh(d) = () -> FrankWolfe.compute_extreme_point(FrankWolfe.LpNormBallLMO{Inf}(1.0), randn(d))
generic_pool(size, dim) = [rand(dim) for _ in 1:size]
generic_fresh(dim) = () -> rand(dim)

# real_size is each alphabet's own max active-set size from
# measurement/results.csv, the same three sizes run_lifecycle.jl swept;
# generic_d3000 has no real BPCG run of its own, so it reuses the L∞-ball's
# real_size as a sizing stand-in, the same convention run_lifecycle.jl's
# own `scenarios` list already uses.
sparse_scenarios = [
    (alphabet=:birkhoff_n25, dim=625, real_size=158, pool=permutation_pool(158 + INSERT_BATCH, 25), fresh=permutation_fresh(25)),
    (alphabet=:birkhoff_n60, dim=3600, real_size=389, pool=permutation_pool(389 + INSERT_BATCH, 60), fresh=permutation_fresh(60)),
]
dense_scenarios = [
    (alphabet=:linf_box_d3000, dim=3000, real_size=241, pool=box_pool(241 + INSERT_BATCH, 3000), fresh=box_fresh(3000)),
    (alphabet=:generic_d3000, dim=3000, real_size=241, pool=generic_pool(241 + INSERT_BATCH, 3000), fresh=generic_fresh(3000)),
]

# Marginal per-atom insert TIME: identical shape to run_lifecycle.jl's
# `insert_cost_ns`, restated here rather than shared, since that file is a
# script, not a module, and has nothing to `include` from. `build_fn` is a
# plain positional argument, never routed through an abstractly-typed
# field, so it keeps the concrete closure type Julia specializes on.
function insert_cost_ns(build_fn, sc, size)
    lo = sc.pool[1:size]
    hi = sc.pool[1:min(size + INSERT_BATCH, length(sc.pool))]
    batch = length(hi) - length(lo)
    batch <= 0 && return NaN
    t_lo = time_per_call(build_fn, lo)
    t_hi = time_per_call(build_fn, hi)
    return max(t_hi - t_lo, 0.0) / batch * 1e9
end

# A bucket-map index's own Dict, copied deeply enough that inserting or
# repairing the copy cannot mutate the original: a shallow `copy(buckets)`
# would still share every stored `Vector{Int}` bucket, so `push!`ing into a
# copy's bucket would corrupt the original's. Generic over the key type, so
# the same function serves all four representations' `Dict{K,Vector{Int}}`.
deepcopy_buckets(buckets::Dict{K,Vector{Int}}) where {K} = Dict{K,Vector{Int}}(key => copy(v) for (key, v) in buckets)

# Single-call allocation of a *mutating* function, not batched: `@allocated`
# has no timer floor to clear, so one clean call after warm-up is the exact
# number. `make_args()` is called three times, once per call below, each
# time producing a *fresh* copy the coming call is free to mutate; crucially,
# each call to `make_args()` happens on its own line, before the timed
# expression, so its own allocation (building that fresh copy) is never
# counted as part of `f`'s (pitfall 2 in this file's header comment).
function mutating_alloc_bytes_and_count(f::F, make_args::A) where {F,A}
    args = make_args()
    f(args...)
    args = make_args()
    bytes = @allocated f(args...)
    args = make_args()
    allocs = @allocations f(args...)
    return bytes, allocs
end

# Same idea for a read-only function: no copy needed, since nothing here
# mutates `index` or `atoms`, so the three calls can share one fixed query.
function readonly_alloc_bytes_and_count(f::F, args...) where {F}
    f(args...)
    bytes = @allocated f(args...)
    allocs = @allocations f(args...)
    return bytes, allocs
end

timing_rows = []
alloc_rows = []
collision_rows = []

# `build`, `lookup`, `key_fn` are plain closure arguments, constrained by
# `where {B,L,KF}` so each keeps its own concrete type through this whole
# function and everything it calls (`insert_cost_ns`, `time_per_call_seq`):
# pitfall 1 in this file's header comment is a struct field losing that
# concreteness, not a `where`-bound argument, which never does.
function measure!(name::String, sc, k::Int, build::B, lookup::L, key_fn::KF) where {B,L,KF}
    size = sc.real_size
    atoms_slice = sc.pool[1:size]
    GC.gc()

    index = build(atoms_slice)
    buckets = index.buckets

    # --- collisions: how many atoms still share a bucket, for context
    # against the timing/allocation numbers below (see
    # results_sparse_pattern_collisions.csv for the k=4/8/16 sweep the
    # vector representation already has; this repeats it at k=2 too and
    # alongside the two new representations). ---------------------------
    stats = bucket_collision_stats(buckets)
    push!(collision_rows, (alphabet=sc.alphabet, representation=name, k=k, size=size, stats...))

    # --- lookup: time (query sequence) + allocation (one warmed-up call) -
    queries = [sc.fresh() for _ in 1:QUERY_N]
    lookup_ns = time_per_call_seq(lookup, (index, atoms_slice), queries) * 1e9
    lookup_bytes, lookup_allocs = readonly_alloc_bytes_and_count(lookup, index, atoms_slice, queries[1])
    GC.gc()

    # --- insert: time (marginal build diff, unchanged atoms) + allocation
    # (one combined "compute this atom's key, then bucket_insert! it" call
    # on a fresh copy, keying the pool's own next atom, size+1: at k=8 the
    # vector representation's own collision rate is 0.0 for both real
    # Birkhoff sizes (results_sparse_pattern_collisions.csv), so "insert a
    # genuinely new key" is the representative case here, not an edge case
    # picked to flatter one representation. Key computation has to be
    # *inside* the measured call, not done once beforehand and reused: an
    # earlier version of this script computed `key_fn(fresh_atom)` outside
    # the timed closure, which hid the one place the three representations
    # actually differ on insert (`bucket_insert!` itself only ever
    # allocates a length-1 bucket Vector{Int}, identical for all three;
    # building the *key* is where `Vector{Int}` pays and `UInt64`/`NTuple`
    # do not), the same way `lookup` above computes its key inside the
    # measured call rather than being handed one pre-built. ---------------
    insert_ns = insert_cost_ns(build, sc, size)
    fresh_atom = sc.pool[size+1]
    insert_bytes, insert_allocs = mutating_alloc_bytes_and_count(
        (buckets_copy, atom, position) -> bucket_insert!(buckets_copy, key_fn(atom), position),
        () -> (deepcopy_buckets(buckets), fresh_atom, size + 1),
    )
    GC.gc()

    # --- delete-repair: allocation first (needs `buckets` still pristine
    # to copy from), then time (which mutates the shared `buckets` in
    # place across its batched repeats, the accepted quirk noted above;
    # nothing downstream still needs `buckets` clean after this). ---------
    pos = max(1, size ÷ 2)
    delete_repair_bytes, delete_repair_allocs = mutating_alloc_bytes_and_count(
        bucket_delete_repair!,
        () -> (deepcopy_buckets(buckets), pos),
    )
    delete_repair_ns = time_per_call(bucket_delete_repair!, buckets, pos) * 1e9
    GC.gc()

    push!(timing_rows, (alphabet=sc.alphabet, representation=name, k=k, size=size, lookup_ns=lookup_ns, insert_ns=insert_ns, delete_repair_ns=delete_repair_ns))
    push!(
        alloc_rows,
        (
            alphabet=sc.alphabet, representation=name, k=k, size=size,
            lookup_bytes=lookup_bytes, lookup_allocs=lookup_allocs,
            insert_bytes=insert_bytes, insert_allocs=insert_allocs,
            delete_repair_bytes=delete_repair_bytes, delete_repair_allocs=delete_repair_allocs,
        ),
    )
    println("$(sc.alphabet) $(name) k=$k: lookup=$(round(lookup_ns,digits=2))ns/$(lookup_bytes)B/$(lookup_allocs)a insert=$(round(insert_ns,digits=2))ns/$(insert_bytes)B/$(insert_allocs)a repair=$(round(delete_repair_ns,digits=2))ns/$(delete_repair_bytes)B/$(delete_repair_allocs)a")
end

for sc in sparse_scenarios, k in K_GRID
    kv = Val(k) # built once per k, not per call: see pattern_key_reps.jl's
    # header on why a freshly-built Val(k) inside a hot loop allocates and
    # a reused one does not.
    measure!("vector", sc, k, atoms -> build_pattern_index(atoms, k), pattern_lookup, atom -> pattern_key(atom, k))
    measure!("uint64", sc, k, atoms -> build_pattern_index_u64(atoms, k), pattern_lookup_u64, atom -> pattern_key_uint64(atom, k))
    measure!("tuple", sc, k, atoms -> build_pattern_index_tuple(atoms, kv), pattern_lookup_tuple, atom -> pattern_key_tuple(atom, kv))
end

for sc in dense_scenarios, k in K_GRID
    measure!("vector_f64", sc, k, atoms -> build_prefix_index(atoms, k), prefix_lookup, atom -> atom[1:k])
end

open(joinpath(@__DIR__, "results_pattern_key_reps_timing.csv"), "w") do io
    println(io, "alphabet,representation,k,size,lookup_ns,insert_ns,delete_repair_ns")
    for r in timing_rows
        println(io, "$(r.alphabet),$(r.representation),$(r.k),$(r.size),$(round(r.lookup_ns,digits=2)),$(round(r.insert_ns,digits=2)),$(round(r.delete_repair_ns,digits=2))")
    end
end

open(joinpath(@__DIR__, "results_pattern_key_reps_allocations.csv"), "w") do io
    println(io, "alphabet,representation,k,size,lookup_bytes,lookup_allocs,insert_bytes,insert_allocs,delete_repair_bytes,delete_repair_allocs")
    for r in alloc_rows
        println(io, "$(r.alphabet),$(r.representation),$(r.k),$(r.size),$(r.lookup_bytes),$(r.lookup_allocs),$(r.insert_bytes),$(r.insert_allocs),$(r.delete_repair_bytes),$(r.delete_repair_allocs)")
    end
end

open(joinpath(@__DIR__, "results_pattern_key_reps_collisions.csv"), "w") do io
    println(io, "alphabet,representation,k,size,n_atoms,n_buckets,mean_bucket_size,max_bucket_size,atom_collision_rate")
    for r in collision_rows
        println(io, "$(r.alphabet),$(r.representation),$(r.k),$(r.size),$(r.n_atoms),$(r.n_buckets),$(round(r.mean_bucket_size,digits=3)),$(r.max_bucket_size),$(round(r.atom_collision_rate,sigdigits=4))")
    end
end

# --- Total per-iteration cost, weighted by the REAL per-iteration call
# rates measurement/results.csv recorded, the same arithmetic
# run_lifecycle.jl's own step 3 does, extended here to also total
# allocated bytes per iteration: the point of this file. Only the three
# alphabets with a real BPCG run behind them get a rate; the generic
# control has none, same as run_lifecycle.jl's own `rate_for`. -----------

function read_csv_dicts(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    return [Dict(header .=> split(line, ",")) for line in lines[2:end]]
end

measured = read_csv_dicts(joinpath(@__DIR__, "..", "measurement", "results.csv"))
rate_for = Dict(:birkhoff_n25 => "birkhoff_n25", :birkhoff_n60 => "birkhoff_n60", :linf_box_d3000 => "linf_box_d3000")

total_rows = []
for sc in vcat(sparse_scenarios, dense_scenarios)
    haskey(rate_for, sc.alphabet) || continue
    row = only(r for r in measured if r["problem"] == rate_for[sc.alphabet])
    iterations_run = parse(Int, row["iterations_run"])
    find_atom_calls = parse(Int, row["find_atom_calls"])
    find_atom_hits = parse(Int, row["find_atom_hits"])
    deleteat_calls = parse(Int, row["deleteat_calls"])
    lookup_rate = find_atom_calls / iterations_run
    insert_rate = (find_atom_calls - find_atom_hits) / iterations_run
    delete_rate = deleteat_calls / iterations_run

    by_rep = Dict((r.representation, r.k) => r for r in timing_rows if r.alphabet == sc.alphabet)
    alloc_by_rep = Dict((r.representation, r.k) => r for r in alloc_rows if r.alphabet == sc.alphabet)

    for ((repname, k), t) in by_rep
        a = alloc_by_rep[(repname, k)]
        total_ns = lookup_rate * t.lookup_ns + insert_rate * t.insert_ns + delete_rate * t.delete_repair_ns
        total_bytes = lookup_rate * a.lookup_bytes + insert_rate * a.insert_bytes + delete_rate * a.delete_repair_bytes
        push!(
            total_rows,
            (
                alphabet=sc.alphabet, representation=repname, k=k, real_size=sc.real_size,
                lookup_rate=lookup_rate, insert_rate=insert_rate, delete_rate=delete_rate,
                lookup_ns=t.lookup_ns, insert_ns=t.insert_ns, delete_repair_ns=t.delete_repair_ns,
                total_ns_per_iter=total_ns,
                lookup_bytes=a.lookup_bytes, insert_bytes=a.insert_bytes, delete_repair_bytes=a.delete_repair_bytes,
                total_bytes_per_iter=total_bytes,
            ),
        )
    end
end

open(joinpath(@__DIR__, "results_pattern_key_reps_total.csv"), "w") do io
    println(
        io,
        "alphabet,representation,k,real_size,lookup_rate_per_iter,insert_rate_per_iter,delete_rate_per_iter,lookup_ns,insert_ns,delete_repair_ns,total_ns_per_iter,lookup_bytes,insert_bytes,delete_repair_bytes,total_bytes_per_iter",
    )
    for r in sort(total_rows; by=r -> (r.alphabet, r.k, r.total_ns_per_iter))
        println(
            io,
            "$(r.alphabet),$(r.representation),$(r.k),$(r.real_size),$(round(r.lookup_rate,sigdigits=4)),$(round(r.insert_rate,sigdigits=4)),$(round(r.delete_rate,sigdigits=6)),$(round(r.lookup_ns,digits=2)),$(round(r.insert_ns,digits=2)),$(round(r.delete_repair_ns,digits=2)),$(round(r.total_ns_per_iter,digits=3)),$(r.lookup_bytes),$(r.insert_bytes),$(r.delete_repair_bytes),$(round(r.total_bytes_per_iter,digits=2))",
        )
    end
end

println()
println("Total per-iteration cost (ns, and allocated bytes), at each alphabet's own real max active-set size:")
for r in sort(total_rows; by=r -> (r.alphabet, r.k, r.total_ns_per_iter))
    println("  $(r.alphabet) (size=$(r.real_size)) $(r.representation) k=$(r.k): $(round(r.total_ns_per_iter,digits=3)) ns/iter, $(round(r.total_bytes_per_iter,digits=2)) bytes/iter")
end
