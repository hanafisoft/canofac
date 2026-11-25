# cf_lang() — langue des sorties

```
R       cf_lang("fr") ; cf_lang("en") ; cf_lang()        # ou options(canofac.lang = "fr")
Python  cf_lang("fr") ; cf_lang("en") ; cf_lang()        # ou variable d'environnement CANOFAC_LANG=fr
```

Fixe ou renvoie la langue des impressions, des titres, des axes, des légendes et des issues du test
(« tous les blocs », « plusieurs blocs », « un seul bloc », « aucun bloc »). Les valeurs stockées dans les
objets (par exemple `link = "all"`) ne sont pas traduites : le code qui les utilise ne dépend pas de la langue.
R et Python lisent le même dictionnaire (`shared/labels.csv`).

L'aide existe aussi dans les deux langues : `cf_help(sujet, "fr")` ou `cf_help(sujet, "en")`.
