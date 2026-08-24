# frankwolfe-active-set-lookup

[![CI](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml/badge.svg)](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml)
[![Licence](https://img.shields.io/badge/licence-MIT-lightgrey)](LICENSE)

> [Read in English](README.md)

**Non — pour les problèmes et tailles d'ensemble actif mesurés ici, hacher
l'ensemble actif de `FrankWolfe.jl` n'apporterait pas d'avantage mesurable.**
Le balayage linéaire de `find_atom` a pris moins de **0,15 % du temps
d'exécution** sur trois exécutions réelles de blended pairwise conditional
gradient (BPCG), ensembles actifs jusqu'à 389 atomes. Un micro-benchmark
isolé montre pourquoi : hacher un atome coûte toujours O(dimension), tandis
que le `!=` du balayage s'arrête dès la première coordonnée différente —
en général la première ou la deuxième pour des atomes indépendants. Un
`Dict` ne l'emporte qu'une fois l'ensemble actif, ou la ressemblance des
atomes, assez grand pour combler cet écart ; ni l'un ni l'autre ici.

## La question

[`ZIB-IOL/FrankWolfe.jl#244`](https://github.com/ZIB-IOL/FrankWolfe.jl/issues/244),
ouverte depuis 2021 sans aucun commentaire : un
[`OrderedCollections.OrderedSet`](https://github.com/JuliaCollections/OrderedCollections.jl),
ou plus généralement un hachage des atomes de l'ensemble actif, aiderait-il ?
`find_atom` (`active_set.jl:316`) balaie linéairement, appelée par
`active_set_update!` chaque fois qu'aucun indice n'est fourni : à chaque
étape « ajout » de BPCG (`blended_pairwise.jl:374`, `nothing` explicite) et
à chaque étape pairwise du PFW simple (`pairwise.jl:242`), mais jamais dans
l'away-step FW (`afw.jl`), qui suit toujours son propre indice.
`_unsafe_equal` (ligne 499/513) est **exacte** — `!=` par coordonnée pour un
`Array`, `==` pour un tableau creux — ce qui rend un hachage valide : aucune
tolérance n'est sacrifiée.

## Ce qui a été mesuré

`measurement/` exécute BPCG (`lazy=true`) sur trois problèmes, instrumentés
sans toucher au paquet : une méthode chronométrée ajoutée pour `find_atom`
sur le type `ActiveSet` (`measurement/instrumentation.jl`) :

| Problème | Dimension | Itérations | Ensemble actif max | Moyenne | Appels `find_atom` | Part du balayage |
|---|---|---|---|---|---|---|
| Birkhoff, n=25 | 625 | 8 000 | 158 | 138,0 | 159 | 0,135 % |
| Birkhoff, n=60 | 3 600 | 20 000 | 389 | 300,3 | 389 | 0,057 % |
| Boule L∞, d=3 000 | 3 000 | 15 000 | 241 | 170,2 | 240 | 0,022 % |

`microbenchmark/` isole la recherche : balayage linéaire contre `Dict`, sur
des tailles de 1 à 20 000 et des dimensions de 16 à 8 192, selon deux
scénarios — **générique** (coordonnées indépendantes, le `!=` s'arrête
presque aussitôt) et **adversarial** (atomes partageant un même préfixe,
ne différant qu'en dernière coordonnée, forçant chaque comparaison à
parcourir tout le préfixe). Le point de croisement (`Dict` gagnant) :

| Dimension | Atomes génériques | Atomes adversariaux |
|---|---|---|
| 16 | 50 | 5 |
| 128 | 500 | 5 |
| 1 024 | 2 000 | 5 |
| 8 192 | entre 6 500 et 10 000 (bruité — voir `MEASURING.md`) | 10 |

Les deux scripts écrivent leur tableau dans un `results.csv` versionné, à
côté d'eux (`what-is-where.md`).

## La réponse

Les trois exécutions réelles restent loin d'un croisement générique — le
plus grand ensemble actif observé (389) reste sous le seuil de la dimension
128 (500) — et les sommets de `FrankWolfe.jl` (permutations creuses, coins
de boîte denses) se comportent plus comme le cas générique que
l'adversarial, d'où une part toujours sous 0,15 %. **Un `Dict` ne se
rentabiliserait que si l'ensemble actif dépassait quelques centaines à
plusieurs milliers d'atomes, ou si la géométrie rendait les sommets
anormalement faciles à confondre.** Rien de tel ici, sur les polytopes
mesurés — pas une affirmation sur tout polytope possible. `DECISIONS.md`
contient le commentaire prêt à poster et la question ouverte.

## Pour aller plus loin

Exactitude, sites d'appel, l'issue : `references.md`. La machine et son
bruit : `MEASURING.md`. Chaque fichier, en une ligne : `what-is-where.md`.
Ce qui revient à Mohamed de décider : `DECISIONS.md`.
