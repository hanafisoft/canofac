# cf_help() — aide en français ou en anglais

```
R       cf_help("cf_test", "fr") ; cf_help("cf_test", "en") ; cf_help("topics")
Python  cf_help("cf_test", "fr") ; cf_help("cf_test", "en") ; cf_help("topics")
```

En anglais : la page d'aide standard (`?cf_test` en R, `help(cf_test)` en Python). En français : les pages
du dossier `shared/aide/fr`, communes à R et à Python. Sans langue, `cf_help()` suit `cf_lang()`.
La vignette « Prise en main de canofac » (R) existe en anglais et en français.
