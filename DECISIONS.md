# Design decisions, and the calls still open

Every place this repository made a human judgement rather than followed a
measurement, recorded so any of it can be revisited or argued with. It runs
past eighty lines deliberately: it is a running log kept across every stage
of the work, not a document written once at the end.

If you want the short version of what was tried and refused, read
[REJECTED.md](REJECTED.md). If you want how the method works, read
[METHOD.md](METHOD.md). This file is the reasoning underneath both.

## The soundness argument was incomplete, and here is the gap

The prefix-hash write-up says a bucket hit is confirmed against the whole atom
with `_unsafe_equal`, so a shorter hash "can only cost speed, never
correctness". That covers **false positives**, two different atoms landing in
one bucket. It does not cover **false negatives**, and there is one.

`_unsafe_equal` compares with `!=`, so it follows `==` semantics. A `Dict`
follows `isequal` semantics. Those disagree about exactly one thing in
Float64: **the sign of zero**. `0.0 == -0.0` is true, `isequal(0.0, -0.0)` is
false. So two atoms differing only in the sign of a zero are **one atom to the
scan and two atoms to a hash**, and the lookup misses an atom the scan would
have found. Confirming the bucket cannot help, because no bucket is reached.

In FrankWolfe that means `find_atom` returns -1, `active_set_update!` takes
its `push!` branch, and the active set silently gains a duplicate the scan
would have merged, splitting one atom's weight across two entries.

**This is not hypothetical for the very case that was measured.** An
L-infinity ball vertex has the shape `-1.0 .* sign.(gradient)`, and
`sign(0.0)` is `0.0`, so a zero gradient component yields `-0.0`:

    -1.0 .* sign.([0.0, 2.0, -1.0])  ==  [-0.0, -1.0, 1.0]

**The fix is one add per hashed coordinate.** Key the bucket on `prefix .+ 0.0`
rather than on `prefix`: adding zero turns `-0.0` into `0.0` and leaves every
other Float64 bit-identical, infinities and subnormals included. Cost is k
additions, not dimension additions, so it does not disturb the timing result.

NaN goes the harmless way: the scan says not equal, the hash collides and then
fails confirmation, so both agree on not-found.

Demonstrated and pinned in `microbenchmark/test_soundness.jl`, fifteen
assertions, all passing.

**Recommendation:** say this in the issue comment. A maintainer would find it
within a day of trying the idea, and proposing a hashed active set without
mentioning it would be proposing a silent bug.

## sparse-key-and-trie: what changed, and what is still a judgement call

The brief for this stage pointed at the prefix hash's weak spot directly: on Birkhoff
atoms, `k=8` gives 9 buckets for 389 atoms, and it wins despite the key being
nearly useless, not because of it. Two structures were built to fix that, and
a gap in every measurement so far (lookup only, never insert or deletion) was
closed at the same time. This section is the reasoning behind the numbers
README.md's "Lookup is not the whole cost" reports; treat it as the audit
trail for that section, not a restatement of it.

**Idea 1, the sparse-pattern key, is a clean, low-risk win for Birkhoff-like
atoms and should be the headline recommendation for them.** It measures
exactly as the representation argument predicts (9 buckets become 158 or
389, one atom per bucket, `microbenchmark/results_sparse_pattern_collisions.csv`),
it costs the same O(k) as the flattened prefix it replaces, and it has no
signed-zero hazard at all (the key is `Vector{Int}`, not `Vector{Float64}`,
see `microbenchmark/test_soundness.jl`'s new testset). The only judgement
call is scope: it needs `SparseMatrixCSC` and a one-nonzero-per-row-and-column
structure to be this discriminating, so it is a recommendation for
permutation-matrix-like atoms specifically, not a general replacement for the
value-prefix hash.

**Idea 2, the hash trie, does not earn adoption, and this repository is not
recommending it.** It was swept honestly (four coordinate-selection
strategies, `k in {4,8,16}`, `max_depth in {1,2,4}`, at every alphabet's own
real active-set size, `microbenchmark/results_lifecycle_collisions.csv`), not
tuned down to a config chosen to make it look bad. The best config it found
never won on total per-iteration cost at any of the three real sizes, and at
Birkhoff n=25's own size it lost even to the plain scan, because its own
insert cost (building a 16-coordinate key and walking up to 4 dict levels)
outweighs whatever it saves on lookup. The reason is structural, not a
tuning miss: a permutation matrix's flattened coordinates each carry the same
weak, skewed signal (about a 1-in-n chance of being 1) regardless of which
ones a level picks, so no coordinate-*selection* strategy over that
representation can substitute for the representation change idea 1 makes. If
this repository is asked "should the trie's depth just be pushed further for
Birkhoff n=60", the honest answer is: probably, eventually it would fully
resolve (any two distinct permutation matrices must disagree somewhere in
their flattened form), but idea 1 already gets there in one level at lower
cost, so this was not chased further. That is a scope boundary, not a claim
that no depth would work.

**The deletion-rate finding is the most useful thing this branch found
before comparing any structures, and it is not the finding the brief
anticipated.** The brief's framing ("if maintenance under `deleteat!` sinks
all of them, that is the most useful finding available") expected repair
cost to possibly dominate. It does not, for these three runs: `deleteat!` on
an individual atom happened 2, 1, and 0 times across 8,002/20,002/15,002
iterations (`measurement/instrumentation.jl`'s new `DELETIONS` counter,
`measurement/results.csv`'s new `deleteat_calls` column). Every structure's
repair function is real, measured, and in every case costs more per
operation than the raw `deleteat!` shift it sits alongside
(`microbenchmark/results_lifecycle_timing.csv`'s `raw_deleteat` rows), but at
a rate this low it cannot move the total ranking. **This is scoped to BPCG
with `lazy=true` run to `epsilon=1e-9`, the only algorithm actually run
end to end in this repository** (see "Open questions" below, PFW and BCG
were only found by grep, never run). A workload with more frequent drop
steps, a looser tolerance, or a different algorithm could see repair cost
matter far more; nothing here rules that out, and nothing here measured it.

**The repair implementation itself is naive, and a smarter one is an open
question, not something this branch built.** Every repair function here
(`microbenchmark/bucket_lifecycle.jl`'s `bucket_delete_repair!`,
`microbenchmark/hash_trie.jl`'s `trie_repair_node!`) walks every stored
position in the whole index, because there is no way to find which entries
moved without an added layer of indirection this branch did not build (a
stable atom-identity to current-position map, updated once per deletion,
rather than every affected bucket walked separately). Whether that
indirection would be worth the extra bookkeeping on every insert, for a
workload where deletion is not this rare, was not measured.

**A file-length judgement call:** `microbenchmark/hash_trie.jl` is 285
lines, past the ~80-logical-line trigger the ergonomic conventions set.
Its role, "the hash trie: build it, query it, and keep it valid under
`push!`/`deleteat!`", was judged to need the "and" the conventions ask a
long file to avoid, because the four pieces share one recursive
`TrieNode`/`TrieLeaf`/`TrieInternal` structure that splitting across files
would either duplicate or expose across a file boundary for no real
separation of concerns, the way `lookup_methods.jl` already keeps its scan,
full-`Dict`, and prefix-hash lookups (plus their sparse-atom variants) in
one file. Flagging it here rather than silently letting it stand.

## pattern-key-integer-hash: what changed, and what is still a judgement call

The brief this time pointed at the sparse-pattern key's own weak
spot: it is a `Vector{Int}`, so every lookup and every insert allocates
one before the `Dict` is even touched. `microbenchmark/pattern_key_reps.jl`
built two allocation-free representations of the identical key (`UInt64`,
`NTuple{K,Int}`); `microbenchmark/run_pattern_key_reps.jl` measured all
three the same way `run_lifecycle.jl` measured Idea 1 against the flat
prefix hash. This section is the reasoning and the caveats behind the
numbers README.md's "Idea 1, tightened" reports; treat it as the audit
trail, not a restatement.

**The fold pays off, cleanly, and this is not a marginal call.** Zero
bytes and zero allocations per lookup for both new representations, at
every k swept, at both real Birkhoff sizes, against the `Vector{Int}`
key's 80-128 bytes and 2 allocations (scaling with k, since the key
itself is what is being allocated). Lookup and insert time both improved
too, not just allocation count: at k=8, `UInt64` and `NTuple{8,Int}` beat
`Vector{Int}` by roughly 1.6-2x on total per-iteration cost at both real
sizes (1.57-1.62x for `UInt64`/`NTuple` at n=60, 1.81-2.00x at n=25;
README.md's table), and the win gets larger once k is also
lowered, since neither new representation's lookup allocation depends on
k at all. **Recommend `UInt64` over `NTuple{k,Int}` for the actual
proposal**, not because it measured consistently faster (at k=8, `NTuple`
wins insert by 21.3%/6.7% at n=25/n=60 and lookup by 1.6% at n=25, but
loses lookup by 1.5% at n=60, so it is not a clean sweep) but because
`UInt64` takes `k` as an ordinary `Int`, the same as `Vector{Int}` does
today, while
`NTuple{K,Int}` needs `k` fixed as a compile-time type parameter
(`Val(k)`) to actually be allocation-free; a real `ActiveSet` where `k`
might reasonably be a per-call or per-LMO tunable, not a value baked into
a type across a whole codebase, is a better fit for the representation
that does not ask for that. This is a judgement call about API shape, not
a correctness or performance one: if a maintainer is comfortable fixing
`k` as a `Val`, `NTuple{k,Int}` is a legitimate, slightly faster
alternative, and this repository would not object to it being chosen
instead.

**The `UInt64` fold's collision hazard is real, not hypothetical, and
survivable, not by luck.** `test_pattern_key_reps.jl`'s first testset
forces one, on two real Birkhoff atoms whose actual `k=2` patterns are
verifiably different, by narrowing the fold to 1 bit (`bits=1`, a
test-only keyword on the shipped `pattern_key_uint64`, not a separate
reimplementation), and checks the resulting index puts both atoms in one
bucket and still answers every query correctly, because every candidate
is confirmed against the whole atom before being trusted. At the fold's
real, unnarrowed default, this same pair does not collide; no attempt was
made to find a real 64-bit collision among actual Birkhoff atoms at the
sizes and k values this repository measured (158-389 atoms, k in
{2,4,8}), because the birthday bound on a 64-bit space makes that
astronomically unlikely to find by search, and the test does not need one
to demonstrate the property that matters: that the structure survives a
collision when one occurs, not that one never will. **The `NTuple{K,Int}`
representation has no equivalent hazard, and the code says so** (both
`pattern_key_reps.jl`'s header and inline comments, and this file's
"Recommend `UInt64`..." paragraph above), which is itself part of why
`UInt64` was chosen as the primary recommendation despite carrying a
hazard `NTuple` does not: the hazard is closed by an argument
(confirm-before-trust) this repository already relies on for every other
hashed structure here, not by hoping it doesn't come up.

**Two measurement pitfalls were found and fixed while building this,
worth recording because a careless version of this same file would have
reported them as properties of the representations being measured, not
of the harness measuring them:**

1. A closure stored in a struct field typed `::Function` loses its
   concrete type there, so calling it through that field forces a dynamic
   dispatch on every call, which allocates. An early version of
   `run_pattern_key_reps.jl` used a `struct Representation; lookup::Function;
   ...; end` to avoid repeating four near-identical measurement blocks;
   `pattern_lookup_u64` measured 0 bytes called directly and 80 bytes
   called through that field, with nothing about the lookup itself
   different between the two calls. The fix was to pass each
   representation's closures as plain arguments into a `where {L,...}`
   -bound function instead of storing them in a field, which keeps
   Julia's specialization intact end to end.
2. `@allocated f(make_args()...)` counts `make_args()`'s own allocation
   together with `f`'s. Measuring insert/delete-repair allocation needs a
   *fresh* copy of the bucket `Dict` per call (so a mutating function
   never corrupts state a later measurement needs), and building that
   copy inside the same expression `@allocated` times charges the copy's
   cost to the call it was only there to protect: an early version of
   this file reported ~5,700-13,600 bytes for a single `bucket_insert!`
   call, when the real number (confirmed once the copy was moved outside
   the timed expression) is 64 bytes. `run_pattern_key_reps.jl`'s
   `mutating_alloc_bytes_and_count` now builds the fresh copy as its own
   statement, before the `@allocated`/`@allocations` expression starts.

Both are recorded in `run_pattern_key_reps.jl`'s own header comment too,
next to the code that fixes them; repeated here because a mistake this
easy to make silently, in a file whose entire point is reporting
allocation numbers accurately, seemed worth a second, more visible
record.

**A reproducible timing anomaly was found, not smoothed over.**
`NTuple{4,Int}`'s marginal insert time (the same build-at-N-vs-N+100
differencing `run_lifecycle.jl` already uses) measured 135.57 ns at n=60
and 177.65 ns at n=25, against 19-38 ns at k=2 and k=8 for the same two
alphabets. Both figures held at the same magnitude across two independent
full runs of `run_pattern_key_reps.jl` (177.65 ns for n=25 both times,
identically, since `Random.seed!(4)` regenerates the exact same atom
sequence each run), which rules out ordinary timer noise; the leading
hypothesis, not confirmed, is a Julia `Dict`
resize boundary landing inside the differenced atom-count range
specifically for `NTuple{4,Int}`'s 32-byte key width and not for
`NTuple{2,Int}`'s 16 or `NTuple{8,Int}`'s 64. This was not chased further
(would need inspecting `Dict`'s internal table size across the exact
insert sequence, not just timing it), so README.md's k-sweep table marks
that one figure rather than either hiding it or repeating a wrong
conclusion the way trusting it uncritically would.

**A file-length judgement call, same shape as `hash_trie.jl`'s above:**
`microbenchmark/run_pattern_key_reps.jl` is around 320 lines total (about
144 non-comment, non-blank lines), the longest file this branch adds,
past the ~80-logical-line trigger. Its role, "sweep all three pattern-key
representations plus the dense value-prefix baseline for time,
allocation, and collision rate, at every k, then combine with the real
call rates into a total," was judged to need the "and" the same way
`run_lifecycle.jl`'s own, even longer equivalent (213 non-comment lines,
not flagged when that branch added it) already does: one script running
one coherent sweep end to end, where splitting the collision sweep, the
timing sweep, and the final rate-weighted combination into separate files
would mean either re-building every index three times (once per file) or
threading a lot of intermediate state across a file boundary for no real
separation of concerns. Flagging it here for the same reason
`hash_trie.jl`'s entry above does, and noting the inconsistency: this
branch does not know why `run_lifecycle.jl` itself was not flagged when
it was added, only that the same reasoning that exempted it applies here
too.

## public-ready: what changed, and what is still a judgement call

This branch's brief was to turn a private research harness into something
a stranger could read and use, without changing any measured number. It
added `src/` (a small, decomposable module: `keys.jl`, `confirm.jl`,
`index.jl`, the design README.md's opening paragraph and METHOD.md now
describe), `test/test_public_api.jl`, `METHOD.md`, and `REJECTED.md`, and
rewrote `README.md` and `README.fr.md` around a stranger's entry point
rather than the research narrative. This section is the audit trail for
the judgement calls that restructuring made, the same role the two
sections above serve for their own branches.

**`src/` duplicates rather than shares code with `microbenchmark/`, on
purpose.** `sparse_pattern.jl`'s `pattern_key`, `pattern_key_reps.jl`'s
`pattern_key_uint64`, `lookup_methods.jl`'s prefix-hash reads, and
`bucket_lifecycle.jl`'s insert/repair are, underneath, the same ideas
`src/keys.jl` and `src/index.jl` now ship. They were not unified into one
shared file, because the two have different jobs: `microbenchmark/`'s
copies exist to be compared against a scan, a full-Dict hash, and a
trie, and stay free to keep changing shape as new ideas get swept;
`src/`'s copies exist to be the one thing a caller actually uses, and
importing `src/` from a benchmark script (or vice versa) would either
make the module depend on files whose whole purpose is churn, or make the
research harness's own numbers depend on the module they are supposed to
be justifying from outside. `what-is-where.md` states this explicitly
rather than leaving it to be noticed. The cost is real: a future fix to
the fold or the repair walk has to be made in two places if both are
meant to keep agreeing, and nothing here enforces that they do beyond
`test/test_public_api.jl` and `microbenchmark/test_*.jl` each checking
their own copy against a scan independently.

**Nothing under `measurement/` or `microbenchmark/` was edited.** Every
file this branch added is new; `git diff` against `main` touches no
existing script, so every committed `results*.csv` is exactly the file
that was already there, and CI's existing steps needed no changes to keep
covering them. Only a new CI step, running `test/test_public_api.jl`, was
added for the new module.

**Naming: `push_atom!`/`delete_atom!`, not `Base.insert!`/`Base.delete!`.**
The public API's "insert" always appends (mirroring
`active_set_update!`'s `push!` branch, the one every real BPCG run this
repository measured actually took), never inserts at an arbitrary
position the way `Base.insert!(vector, i, x)` does; naming it `push_atom!`
says that directly rather than relying on a reader to check the
signature. `k` stays out of `lookup_atom`/`push_atom!`/`delete_atom!`'s
own arguments entirely, read off the index instead: letting a caller pass
a different `k` at lookup time than the index was built with would silently
break every query, the same reasoning `pattern_key_uint64`'s `bits` field
already carries for the research harness's own `PatternIndexU64`.

**This does the restructuring that going public needs, but not the
decision itself.** Making the repository public was the recommended way
out; this branch is
what "public" would need to look like (a stranger has something to read
and a small module to use, not only benchmark scripts and CSVs), but it
does not flip the repository's visibility, and nothing in it was posted
anywhere. Both of those remain deliberate human decisions rather than
anything this work performs on its own.

**A `REJECTED.md`/`METHOD.md` file-length call, same shape as
`hash_trie.jl`'s and `run_pattern_key_reps.jl`'s above.** Both run past
the ~80-logical-line trigger (126 and 124 lines respectively, most of it
prose rather than code); both state their own reason at the top, per this
branch's understanding of the convention, rather than being flagged only
here. `REJECTED.md`'s role, "every rejected method, with the number that
killed it," was judged to need four worked examples in one place rather
than a folder of four short files, since a reader deciding whether to
re-propose one of these needs to compare it against the others on the
same page, not follow four links to do it.

## Open questions, unresolved by this repository

- **Choosing `k`.** `k=8` was the smallest prefix that cleared every real
  alphabet's own observed maximum in `microbenchmark/run_prefix.jl`'s
  sweep, not a value derived from a formula; a polytope with atoms more
  degenerate than Birkhoff's permutation matrices (which needed `k=8`
  themselves) could need more, and this repository has no rule for picking
  it beyond "sweep `k` and check the collision rate."
- **Scope of "real problems".** Only Birkhoff and the L∞-ball were run.
  Both are now covered by `microbenchmark/run_prefix.jl`'s alphabet-matched
  sweep, so the collision behaviour that made the original full-hash answer
  hold (or not) is measured rather than bounded. A problem whose optimum
  genuinely needs thousands of active vertices (a spectrahedron at high
  dimension, an adversarially-designed cost) still was not tried, and a
  polytope whose vertices are more confusable than a permutation matrix's
  flattened prefix (100% bucket collision at `k=8`, and the prefix hash
  still won) has not been found or ruled out.
- **PFW and BCG were found, not measured.** `references.md` and the draft
  above name four more call sites than the task's own framing listed
  (`pairwise.jl`, `blended_cg.jl`, and two generic drivers); only BPCG was
  run end to end. Worth a second harness pass if this becomes a live PR.
  This now also scopes the deletion-rate finding below: the 2/1/0
  deletions counted are BPCG's, at `lazy=true` and `epsilon=1e-9`, not a
  general property of these problems under every algorithm or tolerance.
- **README.fr.md was rewritten by `public-ready`, not by a native
  speaker.** The two-branches-behind gap this item used to describe
  (the plain `k=8` prefix-hash headline, untranslated "Prefix hashing"
  and "Lookup is not the whole cost" sections) no longer exists:
  `public-ready` replaced `README.fr.md` wholesale, following the new
  `README.md`'s own shortened structure rather than translating the old
  one section by section. What is still open is quality, not staleness:
  the translation is this branch's own French, flagged at the top of the
  file itself as wanting a native pass, the same way this item asked for
  one. Read it against `README.md` before either goes out to anyone who
  would notice a wrong idiom.
- **The CI workflow was written, not observed green.** No GitHub run of
  `.github/workflows/ci.yml` has happened yet, now including the two new
  steps this branch added (the lifecycle sweep, `test_soundness.jl`'s
  assertions); it should be watched once pushed, per this repository's own
  no-push rule.
- **The trie's depth was capped at 4 for Birkhoff n=60, and it was not
  pushed further.** 28.8% of atoms still share a leaf at the best swept
  config (`microbenchmark/results_lifecycle_collisions.csv`); depth 8 or
  higher might fully resolve it eventually, since two distinct permutation
  matrices must disagree somewhere in their flattened form, but idea 1
  already resolves the same case in one level at far lower cost, so this
  was never tried. If a future polytope needs the trie specifically
  (idea 1's sparse-structure trick does not generalize to it), how deep is
  actually needed there is still open.
- **The trie's tie-break (prefer shallower `max_depth` over smaller `k`
  when collision rate ties) was a design choice, not something checked
  against its reverse.** A `k=8, max_depth=2` config was never compared
  head to head against the `k=16, max_depth=1` configs the sweep actually
  picked, at equal collision rate, to confirm shallower really is cheaper
  in this codebase's Dict-of-Dict implementation.
- **Repair is O(size) by construction, and no cheaper alternative was
  built or measured.** Every repair function here re-walks the whole
  index on every `deleteat!`; an indirection layer (a stable atom identity
  mapped to a mutable current-position cell, shared between the index and
  a parallel positions array) could in principle make repair proportional
  to how many entries actually moved rather than to the whole index's
  size, at the cost of an extra pointer per stored atom on every insert.
  Whether that trade is worth it depends on a deletion rate this
  repository's own measurement found to be near zero for BPCG, so it was
  not built.
- **`NTuple{4,Int}`'s marginal insert-time anomaly was observed, not
  explained.** 135-178 ns at `k=4`, reproducible across two independent
  runs, against 19-38 ns at `k=2` and `k=8` for the same alphabet and
  size; the leading hypothesis (a `Dict` resize boundary landing inside
  the differenced size range specifically at `NTuple{4,Int}`'s 32-byte key
  width) was not checked against `Dict`'s actual internal table size
  across the insert sequence, only guessed from the shape of the numbers.
  If `NTuple{k,Int}` is ever proposed at `k=4` specifically, this should
  be resolved first, since it currently makes that one config look worse
  than it may actually be.
- **Delete-repair showed no consistent winner among the three pattern-key
  representations**, and no attempt was made to explain the specific
  numbers beyond the structural argument that `bucket_delete_repair!`
  never touches a key at all, only stored position values
  (`bucket_lifecycle.jl`), so key type has no principled reason to move
  this cost. `Vector{Int}` won at n=60 against both new representations;
  `UInt64` won at n=25 but lost at n=60; `NTuple{8,Int}` lost at both.
  Whether this is genuinely representation-independent noise (Dict
  iteration order, which can depend on a key's hash) or a real, smaller
  effect this repository's noise floor (`MEASURING.md`) is too coarse to
  resolve was not determined either way.
- **Whether `UInt64` fold collisions are actually reachable at the sizes
  and `k` values this repository cares about was not searched for.** The
  collision test forces one with a narrowed `bits` keyword rather than
  finding a naturally-occurring one, because the birthday bound on a
  64-bit space makes an exhaustive or random search over 158-389 atoms at
  `k` in {2,4,8} astronomically unlikely to find one. That is an argument
  from probability, not a search that came back empty; if a future
  polytope's atom alphabet is large enough (many thousands of active
  atoms, or a much larger `k`) that the birthday bound stops being
  reassuring, this repository has not measured where that boundary is.
- **Whether `Val(k)`-based dispatch could be made ergonomic enough for
  `FrankWolfe.jl` to actually adopt the `NTuple{k,Int}` representation was
  not explored.** The measurement code threads a precomputed `Val(k)`
  through by hand at each sweep point; a real integration would need `k`
  fixed once per `ActiveSet` construction (plausible, since an active
  set's atom shape does not change mid-run) or accept the representation
  is only viable where that is true. Neither this repository's own
  `PatternIndexTuple{K}` nor its measurement harness tried threading `k`
  through anything resembling `FrankWolfe.jl`'s actual `ActiveSet` type.
