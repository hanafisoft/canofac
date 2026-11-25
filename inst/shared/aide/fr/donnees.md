# Données

## Format d'échange (R et Python)

Un fichier CSV : première colonne, noms des lignes ; puis les variables nommées `<bloc>.x.<variable>` (format
des données potato du package R multiblock). Les colonnes sans `.x.` sont des variables supplémentaires
(par exemple `subtype`).

```
R       cf_write_blocks(blocks, "fichier.csv") ; cf_read_blocks("fichier.csv", blocks = c("A", "B"))
Python  cf_write_blocks(blocks, "fichier.csv") ; cf_read_blocks("fichier.csv", blocks=["A", "B"])
```

## Jeux des articles

| Jeu | Article | R | Python | Prétraitement |
|---|---|---|---|---|
| Friday (16 étangs, 10 blocs) | A | `cf_friday()` (ade4) | `cf_friday(fichier)` | centrage, inertie 1 |
| breast.TCGA (150 tumeurs, 3 blocs) | A | `cf_tcga()` (mixOmics) | `cf_tcga(fichier)` | centrage, réduction, inertie 1 |
| potato (26 échantillons ; NIRraw, FIDraw, CPMGraw, mesurés sur les tubercules crus) | A | `cf_potato(blocks = "A")` (multiblock) ou `cf_potato(fichier, "A")` | `cf_potato(fichier, blocks="A")` | SNV des spectres NIR (déjà dans multiblock), centrage, inertie 1 |
| potato (26 échantillons ; Chemical, Compression, Sensory) | BB | `cf_potato()` (multiblock) ou `cf_potato(fichier)` | `cf_potato(fichier)` | centrage, inertie 1 |

Les données ne sont pas copiées dans canofac : elles sont lues dans les packages R qui les distribuent, ou
dans des fichiers. Pour Python, les exporter depuis R, par exemple
`cf_write_blocks(cf_friday(preprocess = FALSE), "friday.csv")`.
