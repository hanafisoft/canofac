"""Help in English or French.

English: the docstrings (``help(cf_test)``). French: the pages of
``shared/aide/fr``, the same files as the French help of the R package
(``cf_help(topic, "fr")`` in both languages).
"""
import os
import inspect
from .lang import _tr, cf_lang

_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "shared", "aide")


def _topics():
    d = os.path.join(_DIR, "fr")
    return sorted(f[:-3] for f in os.listdir(d) if f.endswith(".md")) if os.path.isdir(d) else []


def cf_help(topic="canofac", lang=None):
    """Print the help page of a topic in English (docstring) or French.

    Parameters
    ----------
    topic : str
        Function name (``"cf_test"``...) or ``"canofac"`` for the overview;
        ``cf_help("topics")`` lists the topics.
    lang : {"en", "fr"}, optional
        Language; default: the current language (:func:`cf_lang`).
    """
    lang = lang or cf_lang()
    if topic == "topics":
        print(", ".join(_topics()))
        return
    if lang == "fr":
        f = os.path.join(_DIR, "fr", f"{topic}.md")
        if not os.path.exists(f):
            print(_tr("help_missing") % (topic, ", ".join(_topics())))
            return
        with open(f, encoding="utf-8") as h:
            print(h.read())
        return
    import canofac
    obj = canofac if topic in ("canofac-package", "overview") else getattr(canofac, topic, None)
    if obj is None:
        print(_tr("help_missing") % (topic, ", ".join(_topics())))
        return
    print(inspect.getdoc(obj))
