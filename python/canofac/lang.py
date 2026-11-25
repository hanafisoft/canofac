"""Language of the printed and plotted results (English or French).

The labels come from the file ``shared/labels.csv``, the same file as in the R
package (``inst/shared/labels.csv``), so that both versions print and plot the
same words.
"""
import csv
import os

_LANG = {"value": os.environ.get("CANOFAC_LANG", "en") if os.environ.get("CANOFAC_LANG", "en") in ("en", "fr") else "en"}
_HERE = os.path.dirname(os.path.abspath(__file__))
_DICT = {}


def _load():
    if _DICT:
        return _DICT
    with open(os.path.join(_HERE, "shared", "labels.csv"), encoding="utf-8") as f:
        for row in csv.DictReader(f):
            _DICT[row["key"]] = (row["en"].replace("\\n", "\n"), row["fr"].replace("\\n", "\n"))
    return _DICT


def cf_lang(lang=None):
    """Set or return the language of the printed and plotted results.

    Parameters
    ----------
    lang : {"en", "fr"}, optional
        ``"en"`` (English, default) or ``"fr"`` (French). If omitted, the
        current language is returned.

    Returns
    -------
    str
        The previous language when ``lang`` is given, else the current one.

    Notes
    -----
    Printed summaries, plot titles, axis labels, legends and the outcomes of
    the test are translated. The values stored in the result objects (for
    example ``link = "all"``) are not, so that code that uses them does not
    depend on the language. The language can also be set with the
    environment variable ``CANOFAC_LANG``. French help: ``cf_help("cf_lang", "fr")``.
    """
    cur = _LANG["value"]
    if lang is None:
        return cur
    if lang not in ("en", "fr"):
        raise ValueError("lang must be 'en' or 'fr'")
    _LANG["value"] = lang
    return cur


def _tr(key):
    d = _load()
    if key not in d:
        return key
    return d[key][1 if _LANG["value"] == "fr" else 0]


def _link_label(link):
    return [_tr(str(x)) for x in link]
