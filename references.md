# References

Papers are cited, never redistributed.

**`besancon2025`**: Mathieu Besançon, Alejandro Carderera, Sebastian
Pokutta, Elias Wirth. *FrankWolfe.jl: A high-performance and flexible
toolbox for Frank-Wolfe algorithms and Conditional Gradients.* ACM
Transactions on Mathematical Software (2025). arXiv:
[2501.14613](https://arxiv.org/abs/2501.14613). The library this repository
measures from outside: the paper the package citation asks for, and the
description of the active set as the package's own maintainers give it. The
`ActiveSetQuadraticProductCaching` type used in the package's own
`examples/birkhoff_polytope.jl` (not used in this repository's harness,
which measures the plain `ActiveSet` the issue is about) is the paper's
specialised-subtype precedent this repository's `DECISIONS.md` points back
to.

**`canon1968`**: Michael D. Canon, Clifton D. Cullum. *A tight upper bound
on the rate of convergence of Frank-Wolfe algorithm.* SIAM Journal on
Control **6** (1968), no. 4, 509-516. The `1/t` rate of plain Frank-Wolfe
is tight when the optimum lies on a face of the polytope rather than at a
corner: the iterates zigzag between that face's corners. Why `guide/`'s
plain run, whose target mixes a few corners, is still `1e-4` from the
optimum after thousands of iterations, and why `test/test_guide.jl`
checks the rate rather than a tolerance.

**`frankwolfe1956`**: Marguerite Frank, Philip Wolfe. *An algorithm for
quadratic programming.* Naval Research Logistics Quarterly **3** (1956),
no. 1-2, 95-110. The original method: minimise a convex function over a
polytope by repeatedly moving toward the vertex the current gradient likes
best. Every algorithm this repository runs (and the active set every one of
them keeps a convex combination of vertices in) descends from this paper.

**`guelat1986`**: Jacques Guélat, Patrice Marcotte. *Some comments on
Wolfe's 'away step'.* Mathematical Programming **35** (1986), 110-119.
DOI: [10.1007/BF01589445](https://doi.org/10.1007/BF01589445). Plain
Frank-Wolfe converges linearly when the optimum lies in the polytope's
relative interior, and the away step recovers a linear rate otherwise:
the contrast `guide/README.md` draws, and the start of the line of
active-set variants whose membership question this repository is about.

**`jaggi2013`**: Martin Jaggi. *Revisiting Frank-Wolfe: Projection-Free
Sparse Convex Optimization.* ICML 2013, PMLR **28** (1), 427-435.
<https://proceedings.mlr.press/v28/jaggi13.html>. The bounds
`test/test_guide.jl` asserts on plain Frank-Wolfe: `f(x_t) - f* <=
2 C_f / (t + 2)` (Theorem 1) and a dual gap at most `6.75 C_f / (t + 2)`
somewhere among the first `t` iterates (Theorem 2, with `beta = 27/8`),
where the curvature constant `C_f` is at most `L` times the squared
diameter, `2n` for permutation matrices. The dual gap as a stopping
certificate, which every algorithm in `guide/` and in FrankWolfe.jl uses,
is this paper's framing.

**`lacoste2015`**: Simon Lacoste-Julien, Martin Jaggi. *On the Global Linear
Convergence of Frank-Wolfe Optimization Variants.* NeurIPS 2015. arXiv:
[1511.05932](https://arxiv.org/abs/1511.05932). Introduces the away-step and
pairwise variants whose active set is what `find_atom` searches: the
paper `afw.jl` and `pairwise.jl` implement, and the reason a Frank-Wolfe run
keeps an active set of vertices at all rather than just the current iterate.

**`tsuji2022`**: Kazuma K. Tsuji, Ken'ichiro Tanaka, Sebastian Pokutta.
*Pairwise Conditional Gradients without Swap Steps and Sparser Kernel
Herding.* ICML 2022. arXiv:
[2110.12650](https://arxiv.org/abs/2110.12650). Blended pairwise conditional
gradient (BPCG): the algorithm `measurement/run.jl` actually runs, chosen
because it is the one call site (`blended_pairwise.jl`, line 374) that
passes `find_atom` an explicit `nothing` rather than an index it already
has, so its lookup is the one this repository can measure honestly rather
than assume.

## How other implementations decide membership

Read for the certificate, to place it against what exists rather than
against nothing. Each is a code base, not a paper, because no paper found
says how its active set is stored.

**`copt`**: Fabian Pedregosa et al., *copt: composite optimization in
Python*, <https://github.com/openopt/copt>. `copt/frank_wolfe.py` keeps the
pairwise active set as a `dict` keyed on a *vertex representation* the LMO
itself returns ("a hashable representation of s, for active set
management", `copt/constraint.py`); for the L1 ball that key is the pair
`(sign, index)`. No vector is ever compared: identity comes from the
oracle, and only oracles that supply one support the pairwise variant.

**`linearFW`**: Simon Lacoste-Julien, the code behind `lacoste2015`,
<https://github.com/Simon-Lacoste-Julien/linearFW>. `AFW.m` and `PFW.m`
map a string built from the atom (`hashing`, "a for 0, b for 1", 0-1
vectors only) to its position with a `containers.Map`, "to see whether
already seen vertex"; the README asks anyone with other atoms to "modify
the hashing internal function ... so that you can properly encode the
atoms in your domain in a unique string". Dropped atoms keep their slot at
weight zero.

**`TorchFW`**: <https://github.com/marcojira/TorchFW>. Appends the FW
vertex to its atom list on every step with no membership test at all, so
duplicates simply accumulate: the failure mode a lookup exists to prevent,
tolerated outright.

**FrankWolfe.jl's own history**: `find_atom`'s scan dates from the
`ActiveSet` type's first commit (2021-01-20). Pull request
[#89](https://github.com/ZIB-IOL/FrankWolfe.jl/pull/89) (2021-02-17),
"replace equal operator by isequal", made the comparison itself fast, and
its author noted "I thought about further improving it using the hash of
atoms, but it seems overkill now that's it's not that limiting". Issue #244
below followed eight months later.

**Hash consing**: Eiichi Goto, *Monocopy and associative algorithms in
extended Lisp* (1974); Jean-Christophe Filliâtre, Sylvain Conchon,
*Type-Safe Modular Hash-Consing* (ML Workshop 2006). The general name for
what `copt` and `linearFW` do: give each structurally distinct value one
canonical identity so that equality becomes an identity test. Named here
so the design is called by the field's word rather than described from
scratch.

## The issue this repository answers

**ZIB-IOL/FrankWolfe.jl#244**, "Consider ordered sets for active set
management", opened 2021-10-07:
<https://github.com/ZIB-IOL/FrankWolfe.jl/issues/244>. Labelled `enhancement`
and `good first issue`; no comments in the five years before this repository
was written. Its entire text asks whether
[`OrderedCollections.OrderedSet`](https://github.com/JuliaCollections/OrderedCollections.jl)
or, more generally, a hash kept for the active set's atoms, would help. Still
open with no comments as of 2026-08-27, and no commit or pull request since
mentions `find_atom`, `OrderedSet` or a lookup. This repository's answer,
with the numbers behind it, is `README.md`'s.

## What was searched for and not found

The certificate rests on an observation about BPCG's step rule (METHOD.md):
in the branch that adds the LMO vertex, that vertex cannot already be
active. It was not found stated anywhere: not in `tsuji2022` (whose
Algorithm 1 writes the update as a set union, `S_{t+1} <- S_t ∪ {w_t}`,
which hides the question), not in the "Conditional Gradient Methods"
survey (Braun, Carderera, Combettes, Hassani, Karbasi, Mokhtari, Pokutta,
arXiv:[2211.14103](https://arxiv.org/abs/2211.14103), whose away-step
listing does split the coefficient update on `v ∉ S_t` versus `v ∈ S_t`),
not on the authors' blog post introducing BPCG, and not in
`blended_pairwise.jl`, whose branch is commented only "add to active set".
Queries run: the phrases "already in the active set", "already belongs to
the active set", "already contained in the active set", "cannot be in the
active set", with and without "pairwise", "blended", and the authors'
names, plus greps for the same and for `∉` over the three papers' text.
That is a statement about what was searched, not a claim that nobody has
noticed it.
