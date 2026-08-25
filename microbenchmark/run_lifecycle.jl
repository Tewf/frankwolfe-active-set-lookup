# `sparse-key-and-trie`'s brief: every measurement up to this point timed
# lookups only, and `measurement/results.csv` shows the three real BPCG
# runs produced zero `find_atom` hits, i.e. every single call was
# followed by a `push!`. The real per-iteration cost is lookup plus
# insert, and a structure also has to survive `active_set_cleanup!`
# calling `deleteat!`, which shifts every later index down by one, so any
# auxiliary index built on top has to repair itself too. This file times
# all three, separately, for four candidate structures:
#
#   - scan: today's `find_atom`, no auxiliary index at all. Its "insert"
#     is `push!` on the atoms Vector and nothing else; it has no repair
#     cost because it has nothing to repair.
#   - prefix (k=8): `run_prefix.jl`'s flat value-prefix hash, the
#     existing recommendation, included so the new structures below are
#     measured against it on the same footing, not just against the scan.
#   - pattern (k=8, Birkhoff only): idea 1, keyed on
#     `SparseMatrixCSC.rowval` (`sparse_pattern.jl`) instead of flattened
#     coordinate values.
#   - trie: idea 2, `hash_trie.jl`'s recursive coordinate-block index, at
#     whichever (k, max_depth, coordinate order) this file's own sweep
#     picks best per alphabet.
#
# Insert cost is measured by differencing, not by repeatedly mutating one
# structure: `build_*` is a pure function of the atoms slice it is given,
# so timing it at size N and at size N+INSERT_BATCH and dividing the
# difference by INSERT_BATCH gives the marginal per-atom insert cost
# without ever repeating a call on state a previous call already changed
# (see ../DECISIONS.md for why a naive repeated-insert loop would not be
# stable to time). Delete-repair cost has no such problem: every repair
# function here walks every stored position doing a compare-and-maybe-
# decrement, which does the same amount of work whether or not a
# previous call in the same timed batch already ran, so it is timed by
# direct repeated calls, the same way a lookup is.
#
# Atom pools are generated the same way `run_prefix.jl`'s are (real LMOs,
# Birkhoff kept sparse), so a number here is directly comparable to that
# file's own crossover and collision tables.
include("lookup_methods.jl")
include("sparse_pattern.jl")
include("hash_trie.jl")
include("bucket_lifecycle.jl")
include("timing.jl")
using .LookupMethods, .SparsePatternLookup, .HashTrie, .BucketLifecycle, .Timing
using FrankWolfe, Random, SparseArrays

Random.seed!(3) # distinct from run_prefix.jl's seed!(2) and run.jl's seed!(1)

const SIZES = [50, 100, 158, 200, 241, 300, 389]
const MAX_SIZE = maximum(SIZES)
const K_GRID = [4, 8, 16]
const DEPTH_GRID = [1, 2, 4]
const QUERY_N = 200
const INSERT_BATCH = 100 # atoms added between the two builds insert cost is differenced from
# Pools are generated INSERT_BATCH atoms larger than MAX_SIZE, so
# insert_cost_ns always has room to slice size+INSERT_BATCH even when
# size == MAX_SIZE (birkhoff_n60's own real_size), without which the
# largest SIZES point would silently have no atoms left to grow into.
const POOL_SIZE = MAX_SIZE + INSERT_BATCH

permutation_pool(size, n) = [FrankWolfe.compute_extreme_point(FrankWolfe.BirkhoffPolytopeLMO(), randn(n, n)) for _ in 1:size]
permutation_fresh(n) = () -> FrankWolfe.compute_extreme_point(FrankWolfe.BirkhoffPolytopeLMO(), randn(n, n))
box_pool(size, d) = [FrankWolfe.compute_extreme_point(FrankWolfe.LpNormBallLMO{Inf}(1.0), randn(d)) for _ in 1:size]
box_fresh(d) = () -> FrankWolfe.compute_extreme_point(FrankWolfe.LpNormBallLMO{Inf}(1.0), randn(d))
generic_pool(size, dim) = [rand(dim) for _ in 1:size]
generic_fresh(dim) = () -> rand(dim)

# alphabet => (dim, sparse?, real max active-set size the matching real
# BPCG run reached, per measurement/results.csv)
scenarios = [
    (alphabet=:birkhoff_n25, dim=625, sparse=true, real_size=158, pool=permutation_pool(POOL_SIZE, 25), fresh=permutation_fresh(25)),
    (alphabet=:birkhoff_n60, dim=3600, sparse=true, real_size=389, pool=permutation_pool(POOL_SIZE, 60), fresh=permutation_fresh(60)),
    (alphabet=:linf_box_d3000, dim=3000, sparse=false, real_size=241, pool=box_pool(POOL_SIZE, 3000), fresh=box_fresh(3000)),
    (alphabet=:generic_d3000, dim=3000, sparse=false, real_size=241, pool=generic_pool(POOL_SIZE, 3000), fresh=generic_fresh(3000)),
]

# --- Step 1: sweep (k, max_depth, coordinate order) by collision/depth
# stats alone (no timing) at every SIZES point, then pick the config with
# the fewest atoms still sharing a leaf at that alphabet's own real_size,
# tie-broken by shallower max_depth then smaller k (both cheaper to
# build and to walk). ------------------------------------------------

# first_k, random and selectivity do not depend on k, so each is built
# once per (alphabet, size); strided's own stride is derived from k
# (`cld(dim, k)`), so it is rebuilt inside the k loop below instead of
# being shared across K_GRID the way the other three are.
fixed_orders(sc, atoms_slice) = (
    first_k=first_k_order(sc.dim),
    random=random_order(sc.dim, MersenneTwister(7)),
    selectivity=selectivity_order(atoms_slice, sc.dim),
)

trie_sweep_rows = []
best_trie_config = Dict{Symbol,NamedTuple}()

for sc in scenarios
    best = nothing
    for size in SIZES
        atoms_slice = sc.pool[1:size]
        fixed = fixed_orders(sc, atoms_slice)
        for k in K_GRID
            k > sc.dim && continue
            orders = merge(fixed, (strided=strided_order(sc.dim, k),))
            for (order_name, order) in pairs(orders), max_depth in DEPTH_GRID
                trie = build_trie(atoms_slice, order, k, max_depth)
                stats = trie_stats(trie)
                push!(
                    trie_sweep_rows,
                    (alphabet=sc.alphabet, dim=sc.dim, size=size, order=order_name, k=k, max_depth=max_depth, stats...),
                )
                if size == sc.real_size
                    score = (stats.atom_collision_rate, max_depth, k) # fewest stragglers first, then cheapest
                    if best === nothing || score < best.score
                        best = (score=score, order_name=order_name, order=order, k=k, max_depth=max_depth, stats=stats)
                    end
                end
            end
        end
    end
    best_trie_config[sc.alphabet] = best
    println(
        "$(sc.alphabet): best trie = order=$(best.order_name) k=$(best.k) max_depth=$(best.max_depth) ",
        "atom_collision_rate=$(round(best.stats.atom_collision_rate,sigdigits=4)) mean_leaf_size=$(round(best.stats.mean_leaf_size,digits=2)) ",
        "mean_depth=$(round(best.stats.mean_depth,digits=2))",
    )
end

open(joinpath(@__DIR__, "results_lifecycle_collisions.csv"), "w") do io
    println(
        io,
        "alphabet,dim,size,order,k,max_depth,n_atoms,n_leaves,mean_leaf_size,max_leaf_size,atom_collision_rate,mean_depth,max_depth_reached",
    )
    for r in trie_sweep_rows
        println(
            io,
            "$(r.alphabet),$(r.dim),$(r.size),$(r.order),$(r.k),$(r.max_depth),$(r.n_atoms),$(r.n_leaves),$(round(r.mean_leaf_size,digits=3)),$(r.max_leaf_size),$(round(r.atom_collision_rate,sigdigits=4)),$(round(r.mean_depth,digits=3)),$(r.max_depth_reached)",
        )
    end
end

# --- Step 2: lifecycle timing (lookup, insert, delete-repair) for scan,
# prefix (k=8), pattern (k=8, Birkhoff only), and each alphabet's best
# trie config, at every SIZES point. ----------------------------------

lifecycle_rows = []

function time_lookup_scan(sc, atoms_slice, size)
    queries = [sc.fresh() for _ in 1:QUERY_N]
    scan = sc.sparse ? sparse_linear_scan : linear_scan
    return time_per_call_seq(scan, (atoms_slice,), queries) * 1e9
end

function time_lookup_prefix(sc, atoms_slice, size, k)
    build = sc.sparse ? build_sparse_prefix_index : build_prefix_index
    lookup = sc.sparse ? sparse_prefix_lookup : prefix_lookup
    index = build(atoms_slice, k)
    queries = [sc.fresh() for _ in 1:QUERY_N]
    return time_per_call_seq(lookup, (index, atoms_slice), queries) * 1e9
end

function time_lookup_pattern(sc, atoms_slice, size, k)
    index = build_pattern_index(atoms_slice, k)
    queries = [sc.fresh() for _ in 1:QUERY_N]
    return time_per_call_seq(pattern_lookup, (index, atoms_slice), queries) * 1e9
end

function time_lookup_trie(sc, atoms_slice, size, order, k, max_depth)
    trie = build_trie(atoms_slice, order, k, max_depth)
    queries = [sc.fresh() for _ in 1:QUERY_N]
    return time_per_call_seq(trie_lookup, (trie, atoms_slice), queries) * 1e9
end

# Marginal per-atom insert cost: build at size and at size+INSERT_BATCH
# off the SAME pool (so the extra atoms are the pool's own next atoms,
# not resampled), difference, divide by INSERT_BATCH. `build_fn` takes
# only the atoms slice; extra positional args are passed through.
function insert_cost_ns(build_fn, sc, size, args...)
    lo = sc.pool[1:size]
    hi = sc.pool[1:min(size + INSERT_BATCH, POOL_SIZE)]
    batch = length(hi) - length(lo)
    batch <= 0 && return NaN
    t_lo = time_per_call(build_fn, lo, args...)
    t_hi = time_per_call(build_fn, hi, args...)
    return max(t_hi - t_lo, 0.0) / batch * 1e9
end

function delete_repair_cost_ns(build_fn, repair_fn, sc, atoms_slice, size, args...)
    index = build_fn(atoms_slice, args...)
    pos = max(1, size ÷ 2)
    return time_per_call(repair_fn, index, pos) * 1e9
end

# `deleteat!` on a plain Vector shrinks it in place, unlike the bucket and
# trie repairs above (which only touch stored Int position values, never
# resize anything): a fixed Vector cannot be timed by `time_per_call`'s
# repeated-call batching, which would eventually run it out of elements.
# Each of `n` independent copies is used exactly once instead, and the
# copying itself happens outside the timed region, so it does not inflate
# the number the way `deleteat!(copy(a), pos)` inside a timed closure
# would (an earlier version of this script did exactly that).
function time_raw_deleteat_ns(atoms_slice, pos; n=2000, repeats=5)
    best = Inf
    for _ in 1:repeats
        copies = [copy(atoms_slice) for _ in 1:n]
        elapsed = @elapsed for c in copies
            deleteat!(c, pos)
        end
        best = min(best, elapsed / n)
    end
    return best * 1e9
end

for sc in scenarios
    build_scan_vec(atoms) = begin
        v = similar(atoms, 0)
        for a in atoms
            push!(v, a)
        end
        v
    end
    build = sc.sparse ? build_sparse_prefix_index : build_prefix_index
    lookup = sc.sparse ? sparse_prefix_lookup : prefix_lookup
    best = best_trie_config[sc.alphabet]

    for size in SIZES
        atoms_slice = sc.pool[1:size]
        GC.gc()

        # scan
        scan_lookup_ns = time_lookup_scan(sc, atoms_slice, size)
        scan_insert_ns = insert_cost_ns(build_scan_vec, sc, size)
        push!(lifecycle_rows, (alphabet=sc.alphabet, size=size, structure="scan", k=0, lookup_ns=scan_lookup_ns, insert_ns=scan_insert_ns, delete_repair_ns=0.0))
        GC.gc()

        # prefix hash, k=8 (the existing recommendation)
        prefix_lookup_ns = time_lookup_prefix(sc, atoms_slice, size, 8)
        prefix_insert_ns = insert_cost_ns(build, sc, size, 8)
        prefix_repair_ns = delete_repair_cost_ns(build, (idx, pos) -> bucket_delete_repair!(idx.buckets, pos), sc, atoms_slice, size, 8)
        push!(lifecycle_rows, (alphabet=sc.alphabet, size=size, structure="prefix_k8", k=8, lookup_ns=prefix_lookup_ns, insert_ns=prefix_insert_ns, delete_repair_ns=prefix_repair_ns))
        GC.gc()

        # sparse-pattern key, k=8 (Birkhoff only: needs SparseMatrixCSC.rowval)
        if sc.sparse
            pattern_lookup_ns = time_lookup_pattern(sc, atoms_slice, size, 8)
            pattern_insert_ns = insert_cost_ns(build_pattern_index, sc, size, 8)
            pattern_repair_ns = delete_repair_cost_ns(build_pattern_index, (idx, pos) -> bucket_delete_repair!(idx.buckets, pos), sc, atoms_slice, size, 8)
            push!(lifecycle_rows, (alphabet=sc.alphabet, size=size, structure="pattern_k8", k=8, lookup_ns=pattern_lookup_ns, insert_ns=pattern_insert_ns, delete_repair_ns=pattern_repair_ns))
            GC.gc()
        end

        # this alphabet's best trie config
        trie_lookup_ns = time_lookup_trie(sc, atoms_slice, size, best.order, best.k, best.max_depth)
        trie_insert_ns = insert_cost_ns((atoms) -> build_trie(atoms, best.order, best.k, best.max_depth), sc, size)
        trie_repair_ns = delete_repair_cost_ns(
            (atoms) -> build_trie(atoms, best.order, best.k, best.max_depth),
            trie_delete_repair!,
            sc,
            atoms_slice,
            size,
        )
        push!(
            lifecycle_rows,
            (alphabet=sc.alphabet, size=size, structure="trie_$(best.order_name)_k$(best.k)_d$(best.max_depth)", k=best.k, lookup_ns=trie_lookup_ns, insert_ns=trie_insert_ns, delete_repair_ns=trie_repair_ns),
        )
        GC.gc()

        # raw deleteat! on the underlying Vector alone: shared infrastructure
        # every structure pays regardless of its own auxiliary index, timed
        # here for context, not folded into any structure's own total.
        pos = max(1, size ÷ 2)
        raw_deleteat_ns = time_raw_deleteat_ns(atoms_slice, pos)
        push!(lifecycle_rows, (alphabet=sc.alphabet, size=size, structure="raw_deleteat", k=0, lookup_ns=NaN, insert_ns=NaN, delete_repair_ns=raw_deleteat_ns))

        println("$(sc.alphabet) size=$size lifecycle timed")
    end
end

open(joinpath(@__DIR__, "results_lifecycle_timing.csv"), "w") do io
    println(io, "alphabet,size,structure,k,lookup_ns,insert_ns,delete_repair_ns")
    for r in lifecycle_rows
        fmt(x) = isnan(x) ? "" : string(round(x, digits=2))
        println(io, "$(r.alphabet),$(r.size),$(r.structure),$(r.k),$(fmt(r.lookup_ns)),$(fmt(r.insert_ns)),$(fmt(r.delete_repair_ns))")
    end
end

# --- Step 3: combine with the REAL per-iteration call rates
# `measurement/run.jl` observed (find_atom_calls and deleteat_calls over
# iterations_run; see instrumentation.jl), read from its own results.csv
# rather than reimplementing it here. This is the one place this
# repository reads a sibling directory's *data*: what-is-where.md notes
# it as the deliberate exception to "no script needs another script's own
# files", since combining a microbenchmark cost with a measured real rate
# is arithmetic on two already-committed CSVs, not a code dependency. ---

function read_csv_dicts(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    return [Dict(header .=> split(line, ",")) for line in lines[2:end]]
end

measured = read_csv_dicts(joinpath(@__DIR__, "..", "measurement", "results.csv"))
rate_for = Dict(
    :birkhoff_n25 => "birkhoff_n25",
    :birkhoff_n60 => "birkhoff_n60",
    :linf_box_d3000 => "linf_box_d3000",
)

total_rows = []
for sc in scenarios
    haskey(rate_for, sc.alphabet) || continue
    row = only(r for r in measured if r["problem"] == rate_for[sc.alphabet])
    iterations_run = parse(Int, row["iterations_run"])
    find_atom_calls = parse(Int, row["find_atom_calls"])
    find_atom_hits = parse(Int, row["find_atom_hits"])
    deleteat_calls = parse(Int, row["deleteat_calls"])
    lookup_rate = find_atom_calls / iterations_run
    insert_rate = (find_atom_calls - find_atom_hits) / iterations_run # every miss is a push!
    delete_rate = deleteat_calls / iterations_run

    for structure_rows in Dict(
        r.structure => r for r in lifecycle_rows if r.alphabet == sc.alphabet && r.size == sc.real_size && r.structure != "raw_deleteat"
    ) |> values
        total_ns = lookup_rate * structure_rows.lookup_ns + insert_rate * structure_rows.insert_ns + delete_rate * structure_rows.delete_repair_ns
        push!(
            total_rows,
            (
                alphabet=sc.alphabet,
                real_size=sc.real_size,
                structure=structure_rows.structure,
                lookup_rate=lookup_rate,
                insert_rate=insert_rate,
                delete_rate=delete_rate,
                lookup_ns=structure_rows.lookup_ns,
                insert_ns=structure_rows.insert_ns,
                delete_repair_ns=structure_rows.delete_repair_ns,
                total_ns_per_iter=total_ns,
            ),
        )
    end
end

open(joinpath(@__DIR__, "results_lifecycle_total.csv"), "w") do io
    println(
        io,
        "alphabet,real_size,structure,lookup_rate_per_iter,insert_rate_per_iter,delete_rate_per_iter,lookup_ns,insert_ns,delete_repair_ns,total_ns_per_iter",
    )
    for r in sort(total_rows; by=r -> (r.alphabet, r.total_ns_per_iter))
        println(
            io,
            "$(r.alphabet),$(r.real_size),$(r.structure),$(round(r.lookup_rate,sigdigits=4)),$(round(r.insert_rate,sigdigits=4)),$(round(r.delete_rate,sigdigits=6)),$(round(r.lookup_ns,digits=2)),$(round(r.insert_ns,digits=2)),$(round(r.delete_repair_ns,digits=2)),$(round(r.total_ns_per_iter,digits=3))",
        )
    end
end

println()
println("Total per-iteration cost (lookup_rate*lookup + insert_rate*insert + delete_rate*repair), at each alphabet's own real max active-set size:")
for r in sort(total_rows; by=r -> (r.alphabet, r.total_ns_per_iter))
    println("  $(r.alphabet) (size=$(r.real_size)) $(r.structure): $(round(r.total_ns_per_iter,digits=3)) ns/iter")
end
