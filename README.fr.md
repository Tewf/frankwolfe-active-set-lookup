# frankwolfe-active-set-lookup

[![CI](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml/badge.svg)](https://github.com/Tewf/frankwolfe-active-set-lookup/actions/workflows/ci.yml)
[![Licence](https://img.shields.io/badge/licence-MIT-lightgrey)](LICENSE)

> [Read in English](README.md)

Suit la structure de `README.md`, donc dépasse aussi les 80 lignes pour
la même raison : la réponse, un exemple utilisable et les commandes de
reproduction tiennent ensemble sur une seule page.

**Cette traduction attend une relecture par un locuteur natif** : elle
suit la nouvelle structure de `README.md` mais n'a pas eu la même
relecture attentive que la version anglaise. En cas de doute, la version
anglaise fait foi.

**Oui : une clé de motif creux repliée bat le balayage linéaire
`find_atom` de `FrankWolfe.jl`, aux tailles d'ensemble actif que des
exécutions réelles atteignent effectivement : 0,812 ns par itération
(recherche, insertion, et réparation occasionnelle après suppression,
pondérées par leur fréquence réelle) contre 57,49 ns pour le balayage, à
Birkhoff n=60 (k=4, la valeur par défaut).** Ceci répond à
[`ZIB-IOL/FrankWolfe.jl#244`](https://github.com/ZIB-IOL/FrankWolfe.jl/issues/244),
ouverte depuis 2021 sans réponse. Plusieurs retournements ont été
nécessaires pour y arriver, dont un premier brouillon concluant
l'inverse ; `REJECTED.md` (en anglais) dit ce qui a été essayé puis
refusé, chiffres à l'appui.

## La méthode, en trois phrases

Pour un atome creux (une matrice de permutation du polytope de Birkhoff,
par exemple), les valeurs valent toutes 1,0 et ne portent aucune
information : la clé hache donc les `k` premières *positions* stockées,
repliées en un seul `UInt64`. Pour un atome dense (un coin de boîte de la
boule L∞), il n'y a pas de structure creuse à lire : la clé hache alors
les `k` premières *valeurs* de coordonnées. Dans les deux cas, la clé ne
fait que désigner une case : chaque candidat qu'elle contient est vérifié
contre l'atome entier par égalité exacte avant d'être accepté, si bien
qu'une collision coûte une comparaison et jamais une mauvaise réponse.
`METHOD.md` (en anglais) donne l'argument complet, y compris une vraie
faille de correction autour du signe de zéro et comment elle est fermée.

## Utilisation

```julia
include("src/ActiveSetLookup.jl")
using .ActiveSetLookup

# atoms est un Vector déjà possédé par l'appelant.
index = build_index(atoms)              # k vaut DEFAULT_K = 4 par défaut
pos = lookup_atom(index, atoms, query)  # -1 si absent, comme find_atom
push_atom!(index, atoms, new_atom)      # garde atoms et index synchronisés
delete_atom!(index, atoms, pos)         # répare l'index après deleteat!
```

`k` est un simple mot-clé partout (`build_index(atoms; k=8)`), jamais un
paramètre de type figé à la compilation. Le module tient en trois
petits fichiers, `src/keys.jl`, `src/confirm.jl`, `src/index.jl`, chacun
lisible seul.

## Ce qui a été mesuré

`measurement/` exécute BPCG (`lazy=true`) sur trois problèmes réels ;
`microbenchmark/` chiffre chaque structure candidate séparément
(recherche, insertion, réparation après suppression), puis les pondère
par le taux d'appel réel de ces trois exécutions. La clé de motif creux
l'emporte sur les deux tailles de Birkhoff (2,72 ns à 158 atomes, 1,88 ns
à 389, avant le repliement `UInt64` sans allocation qui accélère encore) ;
le hachage de préfixe de valeurs existant reste la bonne réponse pour la
boule L∞ (3,31 ns à 241 atomes). Le tableau complet et les trois
exécutions réelles sont dans `README.md`.

## Reproduire les mesures

```
source ~/miniforge3/etc/profile.d/conda.sh && conda activate frankwolfe
julia --project=. measurement/run.jl
julia --project=. microbenchmark/run.jl
julia --project=. microbenchmark/run_prefix.jl
julia --project=. microbenchmark/run_lifecycle.jl
julia --project=. microbenchmark/run_pattern_key_reps.jl
julia --project=. test/test_public_api.jl
```

`MEASURING.md` explique une fois pour toutes comment chaque temps a été
mesuré ; `TESTING.md` couvre les six fichiers de test de correction sous
`microbenchmark/`.

## Pour aller plus loin

Comment la méthode fonctionne et pourquoi elle est correcte : `METHOD.md`.
Ce qui a été essayé puis refusé, chiffres à l'appui : `REJECTED.md`. Ce
que chaque test protège : `TESTING.md`. La machine et son bruit :
`MEASURING.md`. Chaque décision, les questions ouvertes, le commentaire
d'issue en brouillon : `DECISIONS.md`. Articles, sites d'appel,
l'issue elle-même : `references.md`. Chaque fichier, en une ligne :
`what-is-where.md`.
