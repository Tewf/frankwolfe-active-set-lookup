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

## The sparse-pattern key and the trie: what changed, and what is still a judgement call

The brief for this stage pointed at the prefix hash's weak spot directly: on Birkhoff
atoms, `k=8` gives 9 buckets for 389 atoms, and it wins despite the key being
nearly useless, not because of it. Two structures were built to fix that, and
a gap in every measurement so far (lookup only, never insert or deletion) was
closed at the same time. This section is the reasoning behind the numbers
`microbenchmark/results_lifecycle_total.csv` reports.

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
end to end in this repository** (PFW and BCG were run later; see the certificate section below). A workload with more frequent drop
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

## The pattern key's integer representations: what changed, and what is still a judgement call

The brief this time pointed at the sparse-pattern key's own weak
spot: it is a `Vector{Int}`, so every lookup and every insert allocates
one before the `Dict` is even touched. `microbenchmark/pattern_key_reps.jl`
built two allocation-free representations of the identical key (`UInt64`,
`NTuple{K,Int}`); `microbenchmark/run_pattern_key_reps.jl` measured all
three the same way `run_lifecycle.jl` measured Idea 1 against the flat
prefix hash. This section is the reasoning and the caveats behind the
numbers `microbenchmark/results_pattern_key_reps_total.csv` reports.

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
insert sequence, not just timing it), so `microbenchmark/results_pattern_key_reps_total.csv` carries
that one figure as measured rather than either hiding it or repeating a wrong
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

## The public src/ module: what changed, and what is still a judgement call

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

## certificate: what changed, and what is still a judgement call

The brief this time was open: the folded key was measured and shipped, and
the question was whether a better method existed, hash or not. This
section is the audit trail for the answer `README.md` now leads with.

**The finding came from reading the call site, not from a data structure.**
Every earlier stage took `find_atom`'s question as given and made the
search faster. Reading `blended_pairwise.jl` around the one BPCG call
showed the step already held `<g, v>` and the active set's minimum
`<g, s>`, and that the branch is entered only when the first is smaller,
which is a proof of absence. The zero hits `measurement/results.csv` had
recorded from the start were that theorem showing up as data; nobody had
asked why. The literature sweep (`references.md`) found the observation
stated nowhere, which is a statement about what was searched.

**Ship the certificate as the headline and keep the index, rather than
replace one by the other.** The certificate needs a caller that holds the
gradient's products, and every algorithm in `FrankWolfe.jl` that reaches
`find_atom` without an index does. A caller that has only the atom does
not, and for it the index is still the measured best. Both are in `src/`,
the index is the certificate's natural `fallback` on a tie, and neither
file imports the other.

**API shape: mirror `find_atom`.** `certified_lookup` returns a position
or -1, takes the active set's own `Vector`, and takes the two products as
plain numbers rather than a gradient and an argmin function, so the
upstream change at each call site is one line and the module never
recomputes what the step already has. The fingerprint form is a second
method on the same name because it answers the same question from a
different caller state, and `REJECTED.md` section 6 says why it is not the
headline.

**The floating-point assumption is stated and tested, not assumed.** The
certificate compares `dot(g, v)` with `dot(g, a)` values and needs equal
atoms to give equal Float64s. That was checked empirically on this
machine's OpenBLAS and SparseArrays (`test/test_certificate.jl`, 2,400
cases, both argument orders since the package writes them both ways), and
CI runs the same test on its own runner. It is a property of deterministic
kernels, not of IEEE arithmetic in general; a BLAS that reorders a
reduction by thread count would break it, and the fallback direction (a
failed comparison searches) means the failure mode would be a duplicate
atom, never a wrong index.

**The BPCG theorem's scope.** It needs `sparsity_control >= 1` (the default
is 2, and the lazy and non-lazy branches were both checked); it is exact
arithmetic, and the certificate is the guard that makes it safe in floating
point. `inverted` was zero in every run, meaning the LMO vertex never scored
above the active minimum even in floating point on these problems.

**Pairwise Frank-Wolfe was measured for the lookup, not end to end with
the certificate in place.** The 2 ms figure in `README.md` is 20,001 calls
priced at this session's per-call costs; an end-to-end number needs
`pairwise.jl` changed, which is the upstream pull request's job, not this
repository's.

**`results.csv` kept its shape; the other algorithms got their own file.**
`run_lifecycle.jl` and `run_pattern_key_reps.jl` read one row per problem
from `results.csv` by name, so the certificate columns were appended and
PFW, lazy PFW and BCG went to `results_algorithms.csv`, rather than
breaking two readers to keep one file.

**The harness's `RecordingLMO` copies the direction on every LMO call.**
That copy sits inside the timed total and outside the timed scan, so it
can only make every run's `lookup_share` slightly smaller, equally. The
BPCG counts came back identical to the earlier session's (159/389/240
calls, 158/389/241 atoms); the shares moved from 0.02-0.14% to 0.02-0.07%,
inside the noise `MEASURING.md` describes.

**BCG's index discard is reported, not fixed.** `lp_separation_oracle`
returns an active atom without its position and the step scans for it;
that is 98% of BCG's lookups on Birkhoff. It is a one-line upstream change
and belongs in the same pull request, but nothing here builds it.

**A file-length call, same shape as the earlier ones.**
`microbenchmark/run_certificate.jl` is one sweep script like the others,
past the ~80-logical-line trigger for the same reason `run_lifecycle.jl`
was: one coherent sweep whose timing and its weighting by real rates would
otherwise share state across a file boundary for no separation of
concerns.

## A package, and the campaign in its own environments: what changed, and what is still a judgement call

Until 2026-08-28 `src/` was loaded by `include`, the way every script here
loads its files, and `Project.toml` was the environment the scripts ran in.
Making the repository usable from outside a clone meant turning the root
into a package, which is a different object with different rules.

**What changed.** `Project.toml` names the package `ActiveSetLookup`, its
version (the same number as `CITATION.cff`), its one dependency
(`SparseArrays`), compat bounds for everything including the test target,
and that target: `Test`, `Random`, `LinearAlgebra`, `FrankWolfe` for real
atoms, and `Aqua` for the package-quality checks a registry reviewer
would run. `test/runtests.jl` makes `Pkg.test()` real. Every exported name
carries the comment that stood above it as a docstring, so `?certified_lookup`
answers in a session. `index.jl` and `certificate.jl` no longer `include`
`keys.jl` and `confirm.jl` a second time: each submodule is defined once
and reached with `using ..Name`, where before the package held three
copies of `AtomConfirm`. CI runs `Pkg.test()` on Julia 1.10 and the newest
release before the harness job.

**FrankWolfe.jl is not a dependency of the package, by decision.** The
package proposes code to FrankWolfe.jl; a package that depended on it
could never be depended on by it, and `src/` uses nothing from it (`dot`
and the atoms come from the caller). So the campaign, which needs the
library to instrument, runs in its own environments,
`measurement/Project.toml` and `microbenchmark/Project.toml`. Two files
with nearly the same three lines, against the rule that a fact lives in
one place: the alternative, one shared environment somewhere neither
folder owns, would have each folder reaching sideways for its
dependencies. Two small declarations beside the scripts they serve won.

**Still a judgement call.** The version number lives in two files
(`Project.toml`, `CITATION.cff`) because each format demands its own; a
release bumps both in one commit. The Julia 1.10 CI job is what the compat
bound promises and was not run locally, where only 1.12 is installed; the
first CI run on a push decides it. `microbenchmark/` keeps its own copies
of the lookup methods rather than importing the package (the reason is in
`what-is-where.md`), so a change to `src/` is caught by `test/`, never by
the sweeps.

## The end-to-end run: what changed, and what is still a judgement call

The branch `certificate-244` of FrankWolfe.jl (`CONTRIBUTING.md`) applies
the certificate at every call site that reached `find_atom` with no index,
and `measurement/run_end_to_end.jl` times the stock package against it.
Running it caught two things about this repository's own numbers before it
produced the one it was written for.

**The harness had been timing itself.** `instrumentation.jl`'s certificate
tally scans the active set once per `find_atom` call to compare the
certificate with the scan, and it sat inside the `"total"` timer. On
pairwise Frank-Wolfe at n=60 that is 20,001 extra scans of up to 9,368
atoms: the run took 98.8 s with the tally and 68.3 s without. Every share
in the README table was a ratio to an inflated total; the tally is now
timed on its own and `run_time_of` excludes it, and the shares moved from
6.4% to 10.4% at n=60, from 6.05% to 11.8% at n=25, from 0.81% to 1.56% on
the L-inf ball, and by a hundredth of a percent elsewhere. The counts, the
`find_atom` timer itself and every microbenchmark were never affected. The
end-to-end script never arms the tally.

**The registry release is not master.** Registry 0.6.4 predates "update
dual step update in BPCG (#647)"; lazy BPCG's iterates diverge from master
at iteration 2 for that reason alone (389 atoms and one objective against
426 and another). Unpatched master and the branch are byte-identical over
the whole run. So the baseline is master, developed from a checkout into a
temporary environment, and the CSV keeps all three variants so the
difference is on record.

**What the run says, and how far.** Master: 71.4 s, 7.3 s of it in 20,001
scans; the branch: 66.3 s and no scan, same atoms, same objective. The
5.1 s saved is the scan less the certificate's own inner products. Two
runs of identical code in the same session (registry and master on this
problem, same iterates) came out 3 s apart across fastest-of-three, so a
5 s difference is real but its second decimal is not. Lazy BPCG, whose
scans cost 2 ms, came out 0.1 s slower on the branch, inside that noise.

**Still a judgement call.** Upstream, the certificate is enabled for
`ActiveSet` only: the caching active-set types compute their minimum
differently from `dot(a, direction)`, and the argument needs the same
function on the same inputs. `corrective_frankwolfe.jl` keeps the scan, no
minimum being in scope at its call site. The end-to-end script is not in
CI: it needs a checkout of the branch, and its numbers are timings, which
CI asserts nothing about.

## Open questions, unresolved by this repository

- **Pairwise Frank-Wolfe end to end: measured on 2026-08-28**, see the
  section above; it stays listed because the seconds carry the session's
  3 s of noise and only the counts and the iterates are exact.
- **`sparsity_control` below 1.** The BPCG argument needs it at or above
  1; the default is 2 and nothing here checks what the package does with a
  smaller value, or whether any caller uses one.
- **Determinism of `dot` on other BLAS builds.** Checked on this machine's
  OpenBLAS and on CI's runner, both single-configuration. A build that
  splits a reduction across threads by size could give equal atoms
  different sums; the consequence would be a duplicate atom, and
  `test/test_certificate.jl` is the check to run on any new platform.
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
- **The deletion-rate finding is BPCG's.** The 2/1/0 deletions counted
  are BPCG's at `lazy=true` and `epsilon=1e-9`; `results_algorithms.csv`
  now shows BCG on the L-infinity ball removing 299 atoms in 14,757
  iterations, so a workload where an index's repair cost matters does
  exist. The certificate has no repair cost, which is one more reason it
  leads.
- **README.fr.md was rewritten by the public-module stage, not by a native
  speaker.** The two-branches-behind gap this item used to describe
  (the plain `k=8` prefix-hash headline, untranslated "Prefix hashing"
  and "Lookup is not the whole cost" sections) no longer exists:
  The public-module stage replaced `README.fr.md` wholesale, following the new
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
