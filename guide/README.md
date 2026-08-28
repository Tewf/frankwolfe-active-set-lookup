# The guide: Frank-Wolfe from zero to the one question this repository answers

For a reader who has never opened FrankWolfe.jl. Four small files, no
dependency on the library, each readable in a few minutes; then a run and
a test. Longer than the 80-line ceiling `CONTRIBUTING.md` sets, because a walkthrough that stops to link out
before the reader has the picture is not a walkthrough.

```
julia --project=. guide/run.jl        # both algorithms on one problem, narrated
julia --project=. test/test_guide.jl  # every claim below, checked
```

## 1. The problem

Minimise a smooth convex `f(x)` over a **polytope**, a convex set with
finitely many corners. Here: the **Birkhoff polytope** of `n×n` matrices
with nonnegative entries whose rows and columns each sum to one. Its
corners are the **permutation matrices** (one 1 per row and column), and
`f` is half the squared distance to a target matrix (`birkhoff.jl`,
`Quadratic`).

Every point of a polytope is a weighted average of corners:
`x = Σ λᵢ aᵢ`, `λᵢ > 0`, `Σ λᵢ = 1`. Frank-Wolfe only ever produces points
of that form, which is its appeal: the answer comes with the few corners
that make it up.

## 2. The vocabulary

| Word | Means | In `guide/` |
|---|---|---|
| **atom**, vertex | a corner of the polytope | `vertex(p)`, a sparse permutation matrix |
| **gradient `g`** | `∇f(x)` at the current point | `gradient(problem, x)` |
| **LMO** | given `g`, the corner `v` minimising the *linear* function `⟨g, ·⟩` | `extreme_point(lmo, g)`, brute force over all `n!` corners |
| **active set** | the corners currently carrying weight in `x = Σ λᵢ aᵢ` | `ActiveSet` in `active_set.jl`: `atoms`, `weights`, `x` |
| **active** | "has a positive weight in the current mixture", nothing more | `a ∈ A.atoms` |

`⟨g, a⟩` is `dot(g, a)`: for a permutation matrix, the sum of the `n`
entries of `g` the permutation picks. It says how far downhill that corner
lies, seen from `x`, to first order. Lower is better.

## 3. Algorithm one: plain Frank-Wolfe (`frank_wolfe.jl`, `plain_frank_wolfe`)

Each iteration: `g = ∇f(x)`; `v = LMO(g)`; `x ← x + γ (v − x)`. The
**dual gap** `⟨g, x⟩ − ⟨g, v⟩` bounds `f(x) − f*` and is the stop test.
Nothing is remembered but `x`, so no question about membership ever
arises. The price is slow convergence and no explicit list of corners:
`f(x_t) − f* ≤ 2LD²/(t+2)` (Jaggi 2013), and when the optimum lies on a
face rather than at a corner, as here (the target mixes a few corners),
the iterates zigzag between the face's corners and that `1/t` rate is what
they get (Canon and Cullum 1968); only an optimum in the polytope's
relative interior lets plain Frank-Wolfe converge linearly (Guélat and
Marcotte 1986). All three are in `references.md`; `test/test_guide.jl`
checks the bound at every iteration.

## 4. Algorithm two: blended pairwise conditional gradients (`blended_pairwise`)

Keep the active set. Each iteration, compute `⟨g, a⟩` for **every active
atom** (`values`), and keep two of them:

- `s`, the active atom with the **smallest** value: the best corner the
  mixture already holds (the local Frank-Wolfe vertex);
- `a`, the active atom with the **largest** value: the worst one (the away
  vertex).

Two gaps follow. `local_gap = ⟨g, a⟩ − ⟨g, s⟩` is what shifting weight
from `a` to `s` can gain, with no oracle call. `dual_gap = ⟨g, x⟩ − ⟨g, v⟩`
is what the oracle's corner offers. The rule, FrankWolfe.jl's own with its
default `sparsity_control = 2`:

- if `local_gap ≥ dual_gap / 2`: a **pairwise step**, weight moves `a → s`.
  Both are indices already. Nothing to look up.
- otherwise: a **Frank-Wolfe step** toward `v`. Every weight shrinks by
  `(1 − γ)` and `v` gains `γ`. **If `v` is already active, its entry grows;
  if not, it is appended.** Appending an active corner a second time would
  split its weight across two entries and grow the list for nothing, and
  every later loop over the active set would pay for the duplicate. That is
  the membership question, and it is the only place it arises.

## 5. Three ways to answer it (`lookups.jl`)

At that moment the step holds: the atoms, `v`, `g`, the `values`, and `s`.

- **Scan** (`ScanLookup`, FrankWolfe.jl's `find_atom`): compare `v` with
  every atom until one matches. Reads every atom.
- **Index** (`IndexLookup`, `src/index.jl`): hash a few stored positions of
  `v` to a bucket, confirm the candidates exactly. Reads a little of `v`,
  and must be told about every append and removal to stay in sync.
- **Certificate** (`CertificateLookup`, `src/certificate.jl`): compare two
  numbers already computed, `⟨g, v⟩` and `⟨g, s⟩`. Reads no atom at all.

Why the certificate is correct: if `v` were active, `⟨g, v⟩` would be one
of the `values` just minimised, so `⟨g, v⟩ ≥ ⟨g, s⟩`. Therefore
`⟨g, v⟩ < ⟨g, s⟩` proves `v` absent. And `v` minimises `⟨g, ·⟩` over the
whole polytope, which contains the active set, so `⟨g, v⟩ ≤ ⟨g, s⟩` always:
the only other case is a tie, in which `v` is almost surely `s` itself,
settled by one exact comparison at a known index.

Why, in this algorithm, the answer is always "absent": the Frank-Wolfe
branch runs only when `local_gap < dual_gap / 2`. If `v` were active then
`⟨g, v⟩ = ⟨g, s⟩`, and `⟨g, x⟩ ≤ ⟨g, a⟩` because an average never exceeds
its maximum, so `dual_gap ≤ local_gap`, contradicting the branch condition.
`test/test_guide.jl` asserts exactly that on every Frank-Wolfe step of
every run, and `CrossChecked` asserts the three strategies never disagree.

## 6. The worked example, n = 3

`g = [1 5 9; 4 2 8; 7 6 3]`. Active: the identity (`⟨g, I⟩ = 1+2+3 = 6`,
weight 0.7) and the cycle `1→2→3→1` (`5+8+7 = 20`, weight 0.3). So `s` is
the identity with value 6 and `local_gap = 14`. The oracle's cheapest
assignment is the identity, value 6: a tie, `v == s`, found at its index
with one comparison. Now set `g₁₁ = 10`. The cheapest assignment becomes
the swap of rows 1 and 2 (`5+4+3 = 12`), and `12 < 15 = ⟨g, I⟩`, below
every active atom: certified absent, appended, no atom compared. The last
testset in `test/test_guide.jl` is this paragraph, checked.

## 7. What the real library adds, and where to read next

FrankWolfe.jl's BPCG has the same skeleton with a lazy variant (call the
oracle only when the active set cannot deliver), renormalisation, vertex
storage and generic line searches; `active_set_argminmax` is the `values`
loop above, and `active_set_update!` with `nothing` for the index is the
`position_of` call, which today is the scan. The measured runs
(`README.md`), the argument in full (`METHOD.md`) and the two shipped
implementations (`src/`) are the rest of this repository.
