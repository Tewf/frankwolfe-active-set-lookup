# frankwolfe-active-set-lookup

[![CI](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml/badge.svg)](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml)
[![Licence](https://img.shields.io/badge/licence-MIT-lightgrey)](LICENSE)

> [Read in English](README.md)

**Oui : pour les atomes que `FrankWolfe.jl` génère réellement, un hachage
de préfixe court l'emporte sur le balayage linéaire de `find_atom`, aux
tailles d'ensemble actif que les trois exécutions réelles ont atteintes.**
La première version de ce dépôt n'avait testé que le hachage de l'atome
*entier* (toujours en O(dimension)) contre une requête d'échec (le
meilleur cas pour le balayage : le `!=` s'arrête aussitôt) ; aucune de ces
deux restrictions n'est nécessaire, et lever les deux inverse la réponse.
Un hachage sur seulement les `k=8` premières coordonnées, toujours
confirmé par la même égalité exacte que le balayage utilise déjà à chaque
collision de case, l'emporte sur le balayage dès une taille d'ensemble
actif de 100, pour chaque alphabet d'atomes testé (matrices de permutation
de Birkhoff, coins de boîte de la boule L∞, vecteurs aléatoires génériques)
et chaque mélange de requêtes (échec, succès, un mélange de stress), bien
en deçà de la plage 158-389 que les trois exécutions réelles de blended
pairwise conditional gradient (BPCG) ont atteinte. **La réponse à la
question #244 change : un ensemble actif haché par préfixe aiderait.**

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
`_unsafe_equal` (ligne 499/513) est **exacte** (`!=` par coordonnée pour un
`Array`, `==` pour un tableau creux), ce qui rend un hachage valide : aucune
tolérance n'est sacrifiée.

## Ce qui a été mesuré

`measurement/` exécute BPCG (`lazy=true`) sur trois problèmes, instrumentés
sans toucher au paquet : une méthode chronométrée ajoutée pour `find_atom`
sur le type `ActiveSet` (`measurement/instrumentation.jl`) :

| Problème | Dimension | Itérations | Ensemble actif max | Moyenne | Appels `find_atom` | Part du balayage |
|---|---|---|---|---|---|---|
| Birkhoff, n=25 | 625 | 8 000 | 158 | 138,0 | 159 | 0,144 % |
| Birkhoff, n=60 | 3 600 | 20 000 | 389 | 300,3 | 389 | 0,058 % |
| Boule L∞, d=3 000 | 3 000 | 15 000 | 241 | 170,2 | 240 | 0,021 % |

`microbenchmark/` isole la recherche : balayage linéaire contre `Dict`, sur
des tailles de 1 à 20 000 et des dimensions de 16 à 8 192, selon deux
scénarios : **générique** (coordonnées indépendantes, le `!=` s'arrête
presque aussitôt) et **adversarial** (atomes partageant un même préfixe,
ne différant qu'en dernière coordonnée, forçant chaque comparaison à
parcourir tout le préfixe). Le point de croisement (`Dict` gagnant) :

| Dimension | Atomes génériques | Atomes adversariaux |
|---|---|---|
| 16 | 50 | 5 |
| 128 | 500 | 5 |
| 1 024 | 2 000 | 5 |
| 8 192 | entre 6 500 et 10 000 (bruité, voir `MEASURING.md`) | 10 |

Les deux scripts écrivent leur tableau dans un `results.csv` versionné, à
côté d'eux (`what-is-where.md`).

**Ce tableau ne hachait que l'atome entier, et n'interrogeait qu'un échec.**
`microbenchmark/run_prefix.jl` referme les deux : un hachage de préfixe sur
les `k` premières coordonnées seulement, confirmé par la même égalité
exacte à chaque collision de case ; trois types de requête (échec, succès,
un mélange 50/50) ; et trois alphabets d'atomes, dont les vraies formes de
`FrankWolfe.jl` (permutations de Birkhoff, coins de boîte de la boule L∞,
générés en appelant directement ses propres LMO). Le tableau des points de
croisement par `k` et par alphabet, et l'explication du taux de collision
qui les sous-tend, ne sont pas traduits ici : voir `README.md`
("Prefix hashing") et `microbenchmark/results_prefix_crossover.csv`.

## La réponse

Cela dépend entièrement de ce qui est haché. **Hacher l'atome entier reste
une mauvaise idée** : cela coûte toujours O(dimension), et le `Dict`
complet ne l'emporte qu'au-delà de quelques centaines à plusieurs milliers
d'atomes pour des coordonnées indépendantes (de 50 à la dimension 16,
jusqu'à 6 500-10 000 à la dimension 8 192) ; les propres atomes de
`FrankWolfe.jl` donnent en fait un avantage plus tôt à un hachage complet
(500 à 1 500 pour les permutations de Birkhoff, 300 à 1 500 pour les coins
de boîte), sans que cela suffise : aucune des trois exécutions réelles (389
atomes au maximum) n'en a approché, d'où une part toujours sous 0,15 %.

**Hacher un court préfixe est une question différente, avec une réponse
différente : oui.** Un hachage de préfixe `k=8`, confirmé par la même
égalité exacte que le balayage, l'emporte sur `find_atom` pour chaque
alphabet d'atomes que `FrankWolfe.jl` génère réellement et chaque mélange
de requêtes testé, à des tailles d'ensemble actif inférieures au maximum
observé par chacune des trois exécutions réelles. **Cela change la réponse
à la question #244 : un ensemble actif haché par un court préfixe
aiderait.** `DECISIONS.md` contient le commentaire réécrit et la question
ouverte.

## Pour aller plus loin

Exactitude, sites d'appel, l'issue : `references.md`. La machine et son
bruit : `MEASURING.md`. Chaque fichier, en une ligne : `what-is-where.md`.
Ce qui revient à Mohamed de décider : `DECISIONS.md`.
