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

**`frankwolfe1956`**: Marguerite Frank, Philip Wolfe. *An algorithm for
quadratic programming.* Naval Research Logistics Quarterly **3** (1956),
no. 1-2, 95-110. The original method: minimise a convex function over a
polytope by repeatedly moving toward the vertex the current gradient likes
best. Every algorithm this repository runs (and the active set every one of
them keeps a convex combination of vertices in) descends from this paper.

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

## The issue this repository answers

**ZIB-IOL/FrankWolfe.jl#244**, "Consider ordered sets for active set
management", opened 2021-10-07:
<https://github.com/ZIB-IOL/FrankWolfe.jl/issues/244>. Labelled `enhancement`
and `good first issue`; no comments in the five years before this repository
was written. Its entire text asks whether
[`OrderedCollections.OrderedSet`](https://github.com/JuliaCollections/OrderedCollections.jl)
or, more generally, a hash kept for the active set's atoms, would help. This
repository's answer, with the numbers behind it, is `README.md`'s.
