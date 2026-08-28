# frankwolfe-active-set-lookup

[![CI](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml/badge.svg)](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml)
[![Licence](https://img.shields.io/badge/licence-MIT-lightgrey)](LICENSE)

> [Read in English](README.md)

Suit la structure de `README.md`, donc dépasse aussi les 80 lignes que ce
dépôt impose à ses fichiers Markdown (`CONTRIBUTING.md`, en anglais,
énonce la règle et ses exceptions), pour la même raison : la réponse, un
exemple utilisable et les commandes de reproduction tiennent ensemble sur
une seule page. En cas de doute, la version anglaise fait foi.

**Le sommet que l'oracle linéaire (LMO) vient de renvoyer n'est pas dans
l'ensemble actif dès que `<g,v> < <g,s>`, où `s` est l'atome actif que le
pas a déjà trouvé meilleur pour le gradient `g`. Cette seule comparaison
remplace le balayage linéaire de `find_atom` à chaque site d'appel de
`FrankWolfe.jl`, en 10 ns environ contre 2 042 ns pour le balayage sur
Birkhoff n=60, sans index à construire, à alimenter ni à réparer ; une
égalité se règle par une comparaison avec `s`. Sur des exécutions réelles
elle a tranché chaque appel : 788 sur 788 en gradient conditionnel par
paires mélangé (BPCG), 20 001 sur 20 001 en Frank-Wolfe par paires (PFW)
sur Birkhoff n=60, où l'ensemble actif atteint 9 368 atomes et où le
balayage coûte 6,4 % de l'exécution.** Pour un appelant qui ne dispose que
de l'atome, une clé structurelle repliée avec confirmation exacte bat le
balayage à toutes les tailles mesurées et reste ici comme solution de
repli. Ceci répond à
[`ZIB-IOL/FrankWolfe.jl#244`](https://github.com/ZIB-IOL/FrankWolfe.jl/issues/244),
ouverte depuis 2021 sans réponse. Plusieurs retournements ont été
nécessaires, dont un premier brouillon concluant que hacher ne sert à rien
et un deuxième qui s'arrêtait au hachage ; `REJECTED.md` (en anglais) dit
ce qui a été essayé puis refusé, chiffres à l'appui.

**Jamais ouvert FrankWolfe.jl ?** Commencez par [`guide/`](guide/README.md)
(en anglais) : les deux algorithmes réimplémentés en quatre courts
fichiers sans dépendre de la bibliothèque, la question d'appartenance
répondue de trois façons côte à côte, une exécution commentée et un test
de chaque affirmation.

## La méthode, en trois phrases

Tout pas de Frank-Wolfe qui maintient un ensemble actif commence par
minimiser `<g, a>` sur les atomes actifs `a` (c'est ainsi qu'il trouve ses
sommets local et d'éloignement), puis demande au LMO le sommet `v` et
calcule `<g, v>` pour l'écart dual ; si `v` était déjà actif, `<g, v>`
serait l'une des valeurs que l'on vient de minimiser, donc `<g, v> < <g, s>`
prouve qu'il ne l'est pas. Quand la comparaison échoue, `v` fait égalité
avec le meilleur atome `s`, et soit `v == s`, ce qu'une comparaison exacte
règle (c'est le cas de Frank-Wolfe par paires, où le LMO renvoie
couramment un sommet déjà actif), soit `v` fait égalité avec un *autre*
atome, seul cas où l'on cherche vraiment, de probabilité nulle pour un
gradient à valeurs réelles. En BPCG le pas de Frank-Wolfe n'est pris que si
`v` bat tout l'ensemble actif, donc le balayage n'y a jamais été
nécessaire en arithmétique exacte, et le certificat est la façon sûre en
virgule flottante de s'en passer ; `METHOD.md` (en anglais) donne les deux
arguments et ce sur quoi ils reposent.

## L'utiliser

```julia
using Pkg; Pkg.add(url="https://github.com/Tewf/frankwolfe-active-set-lookup")
using ActiveSetLookup
# Depuis un clone sans rien installer, le même module se charge par fichier :
# include("src/ActiveSetLookup.jl"); using .ActiveSetLookup

# Dans un pas. `atoms` est le Vector de l'ensemble actif, `g` le gradient ;
# le pas vient de minimiser dot(g, a) sur atoms (active_set_argminmax), il
# tient donc `best` et `best_value` ; le LMO a renvoyé `v`, et dot(g, v) a
# été calculé pour l'écart dual.
pos = certified_lookup(atoms, v, dot(g, v), best, best_value)  # -1 si absent, comme find_atom

# Un appelant qui n'a que l'atome garde un index à la place.
index = build_index(atoms)              # k vaut DEFAULT_K = 4 par défaut
pos = lookup_atom(index, atoms, query)  # -1 si absent
push_atom!(index, atoms, new_atom)      # garde atoms et index synchronisés
delete_atom!(index, atoms, pos)         # répare l'index après deleteat!

# L'index est aussi le repli naturel pour la rare égalité du certificat.
pos = certified_lookup(atoms, v, dot(g, v), best, best_value;
                       fallback=(a, q) -> lookup_atom(index, a, q))
```

`certified_lookup` ne sait rien du type de l'atome : le même appel sert
aux matrices de permutation, aux coins de boîte et à tout ce que `dot`
accepte. L'index choisit sa clé selon le type des atomes (`SparseMatrixCSC`
et `SparseVector` vers une clé sur les *positions* stockées, un `Array`
dense vers une clé sur les premières *valeurs*), et chaque candidat d'une
case est confirmé contre l'atome entier avant d'être accepté, si bien
qu'une collision coûte une comparaison et jamais une mauvaise réponse.
Derrière `src/ActiveSetLookup.jl`, quatre petits fichiers, `src/certificate.jl`,
`src/keys.jl`, `src/confirm.jl`, `src/index.jl`, chacun lisible seul.

## Ce qui a été mesuré

`measurement/` exécute quatre algorithmes sur trois problèmes,
instrumentés sans modifier `FrankWolfe.jl` : à chaque appel de
`find_atom` il note ce que le balayage a répondu et, pour le sommet du LMO,
ce que le certificat aurait décidé, puis vérifie que les deux concordent.

| Algorithme | Problème | Ensemble actif max | Appels `find_atom` | Trouvés | Part du balayage | Absence certifiée | Égalité, `v == s` | Égalité, recherche |
|---|---|---|---|---|---|---|---|---|
| BPCG paresseux | Birkhoff n=25 | 158 | 159 | 0 | 0,07 % | 159 | 0 | 0 |
| BPCG paresseux | Birkhoff n=60 | 389 | 389 | 0 | 0,06 % | 389 | 0 | 0 |
| BPCG paresseux | Boule L-inf d=3 000 | 241 | 240 | 0 | 0,02 % | 240 | 0 | 0 |
| PFW | Birkhoff n=25 | 2 463 | 8 001 | 5 535 | 6,05 % | 2 466 | 5 535 | 0 |
| PFW | Birkhoff n=60 | 9 368 | 20 001 | 10 628 | 6,41 % | 9 373 | 10 628 | 0 |
| PFW | Boule L-inf d=3 000 | 2 613 | 15 001 | 12 387 | 0,81 % | 2 614 | 12 387 | 0 |
| PFW paresseux | Birkhoff n=25 | 181 | 204 | 0 | 0,17 % | 204 | 0 | 0 |
| PFW paresseux | Birkhoff n=60 | 576 | 598 | 0 | 0,12 % | 598 | 0 | 0 |
| PFW paresseux | Boule L-inf d=3 000 | 295 | 323 | 0 | 0,03 % | 323 | 0 | 0 |
| BCG | Birkhoff n=25 | 148 | 7 948 | 7 800 | 0,14 % | 148 | 38 | 0 |
| BCG | Birkhoff n=60 | 623 | 19 988 | 19 365 | 0,14 % | 623 | 126 | 0 |
| BCG | Boule L-inf d=3 000 | 348 | 1 827 | 1 360 | 0,06 % | 467 | 17 | 0 |

(`measurement/results.csv`, `measurement/results_algorithms.csv` ; 8 000,
20 000 et 15 000 itérations, `epsilon=1e-9`.) Le certificat n'a jamais
contredit le balayage, et aucun appel n'a eu besoin d'une recherche. Trois
choses que la mesure BPCG seule ne pouvait pas dire : le balayage est un
vrai coût en Frank-Wolfe par paires non paresseux, dont l'ensemble actif
grossit jusqu'à des milliers d'atomes (6,4 % d'une exécution de 99 s à
n=60) ; chacun de ses 53 à 83 % d'appels trouvés était le meilleur atome
lui-même, qu'une comparaison suffit à reconnaître ; et les appels de BCG
concernent surtout un atome qu'il vient de prendre dans son propre
ensemble actif (`lp_separation_oracle` renvoie l'atome sans sa position),
ce qui est un gaspillage distinct, propre à l'index.

Par appel, sur les mêmes atomes dans la même session
(`microbenchmark/results_certificate_timing.csv`, ns) :

| Taille | Cas | Certificat | Clé repliée / hachage de préfixe | Balayage |
|---|---|---|---|---|
| Birkhoff n=25, 158 atomes | absent | 10,5 | 28,0 | 777,6 |
| | présent (meilleur atome) | 87,4 | 106,9 | 200,3 |
| Birkhoff n=60, 389 atomes | absent | 10,3 | 27,4 | 2 041,8 |
| | présent (meilleur atome) | 177,8 | 196,3 | 436,9 |
| Boule L-inf d=3 000, 241 atomes | absent | 10,2 | 41,8 | 322,1 |
| | présent (meilleur atome) | 1 392,7 | 1 448,1 | 1 619,9 |

Un absent coûte au certificat une comparaison de Float64 (les 10 ns sont
le plancher du chronomètre pour un appel, le même pour chaque colonne) ; un
présent coûte à chaque méthode l'unique comparaison exacte, soit 1,4 µs
sur un atome dense de 3 000 coordonnées. Une égalité forcée avec un atome
distinct retombe sur le balayage et en coûte le prix (798, 2 091 et
357 ns), et aucune exécution réelle n'en a produit. Pondéré par les taux
d'appel réels de BPCG, le coût total par itération est de 0,16 à 0,21 ns
pour le certificat contre 1,0 à 3,1 ns pour l'index (ses coûts
d'insertion et de réparation venant de sa propre campagne) et 5 à 40 ns
pour le balayage (`results_certificate_total.csv`). Pour Frank-Wolfe par
paires sur Birkhoff n=60, les 6,35 s que le balayage a passées sur 20 001
appels deviennent, aux coûts par appel ci-dessus, environ 2 ms : ce chiffre
est une arithmétique sur des parties mesurées, non une exécution de bout
en bout, car l'obtenir suppose de modifier `pairwise.jl` lui-même.

## Reproduire les mesures

Avec Julia 1.10 ou plus récent, depuis la racine du dépôt :

```
julia --project=. -e 'using Pkg; Pkg.test()'                       # Aqua, puis le certificat, l'index et le guide
julia --project=measurement -e 'using Pkg; Pkg.instantiate()'
julia --project=measurement measurement/run.jl                     # results.csv, results_algorithms.csv
julia --project=microbenchmark -e 'using Pkg; Pkg.instantiate()'
julia --project=microbenchmark microbenchmark/run_certificate.jl   # le certificat face à la clé, au hachage et au balayage
julia --project=microbenchmark microbenchmark/run.jl               # la campagne hachage-complet-contre-balayage
julia --project=microbenchmark microbenchmark/run_prefix.jl        # la campagne du hachage de préfixe
julia --project=microbenchmark microbenchmark/run_lifecycle.jl     # clé de motif creux + trie, cycle de vie complet
julia --project=microbenchmark microbenchmark/run_pattern_key_reps.jl # représentations UInt64/NTuple/Vector{Int}
```

L'environnement du paquet ne contient que `SparseArrays` ; les deux
dossiers de scripts portent chacun un `Project.toml` avec FrankWolfe.jl,
si bien que la bibliothèque mesurée ne devient jamais une dépendance du
code qu'on lui propose. Chaque script écrit un `results*.csv` committé à
côté de lui ;
`MEASURING.md` dit, une fois, comment chaque temps a été pris et ce qu'il
n'affirme pas (aucun n'est vérifié par la CI, qui ne fait que les
exécuter ; la seule assertion est que le certificat ne contredit jamais le
balayage). `TESTING.md` couvre chaque test de correction et ce qu'il
protège.

## Pour aller plus loin

Comment les deux méthodes fonctionnent et pourquoi elles sont correctes, y
compris pourquoi BPCG n'a jamais eu besoin du balayage et la subtilité du
zéro signé de la clé dense : `METHOD.md`. Ce qui a été essayé puis refusé,
chiffres à l'appui, et comment les deux implémentations qui résolvent ce
problème ailleurs (`copt`, `linearFW`) s'y prennent : `REJECTED.md`,
`references.md`. Ce que chaque test protège : `TESTING.md`. La machine et
son bruit : `MEASURING.md`. Chaque choix et chaque question ouverte :
`DECISIONS.md`. Chaque fichier, en une ligne : `what-is-where.md`. Comment
ce dépôt est écrit, comment y travailler et ce qui est proposé en amont,
à l'adresse des mainteneurs : `CONTRIBUTING.md`. Tous en anglais.
