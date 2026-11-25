"""Minimal roxygen -> Rd converter for canofac (used where R/roxygen2 is unavailable).

Handles the tags used in canofac: title (first line), description (first paragraph), details (other
paragraphs), @param, @inheritParams, @return, @examples, @references, @seealso, @rdname, @export, @name.
The authoritative Rd files are those produced by roxygen2 (devtools::document()); this script only keeps
man/ consistent with the R sources between two runs of roxygen2.
Usage: python3 tools/roxy2rd.py <function> [<function> ...]   (run from the package root)
"""
import os, re, sys, glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def blocks_of(path):
    """Yield (roxygen lines, object name, usage) for each documented object of an R file."""
    lines = open(path, encoding="utf-8").read().split("\n")
    i = 0
    while i < len(lines):
        if lines[i].startswith("#'"):
            rox = []
            while i < len(lines) and lines[i].startswith("#'"):
                rox.append(lines[i][2:][1:] if lines[i][2:3] == " " else lines[i][2:])
                i += 1
            j = i
            while j < len(lines) and not lines[j].strip():
                j += 1
            m = re.match(r"^([\w.]+)\s*<-\s*function\s*\(", lines[j]) if j < len(lines) else None
            if m:
                name = m.group(1)
                txt, depth, k, started = "", 0, j, False
                while k < len(lines):
                    seg = lines[k]
                    for c in seg:
                        if c == "(":
                            depth += 1; started = True
                        elif c == ")":
                            depth -= 1
                    txt += seg.strip() + " "
                    if started and depth == 0:
                        break
                    k += 1
                args = txt[txt.index("(") + 1: txt.rindex(")")]
                usage = f"{name}({re.sub(r'\s+', ' ', args).strip()})"
                usage = usage.replace("( ", "(").replace(" )", ")")
                yield rox, name, usage
            else:
                yield rox, None, None
        else:
            i += 1


def parse(rox):
    tags, cur = {"_text": []}, "_text"
    params = []
    for l in rox:
        m = re.match(r"^@(\w+)\s*(.*)$", l)
        if m:
            t, rest = m.group(1), m.group(2)
            if t == "param":
                nm, _, d = rest.partition(" ")
                params.append([nm, d])
                cur = ("param", len(params) - 1)
            else:
                tags.setdefault(t, []).append(rest)
                cur = t
        else:
            if isinstance(cur, tuple):
                params[cur[1]][1] += "\n" + l
            else:
                tags.setdefault(cur, []).append(l)
    tags["_params"] = params
    return tags


def all_docs():
    docs = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "R", "*.R"))):
        for rox, name, usage in blocks_of(f):
            t = parse(rox)
            key = name or (t.get("name", [None])[0])
            if key:
                docs[key] = dict(tags=t, usage=usage, file=os.path.basename(f))
    return docs


def wrap_usage(u, width=78):
    out, line = [], ""
    for tok in re.split(r"(?<=,) ", u):
        if len(line) + len(tok) + 1 > width and line:
            out.append(line.rstrip())
            line = "  " + tok + " "
        else:
            line += tok + " "
    out.append(line.rstrip())
    return "\n".join(out)


def para(lines):
    txt = "\n".join(lines).strip("\n")
    return [p.strip() for p in re.split(r"\n\s*\n", txt) if p.strip()]


def build(target, docs):
    d = docs[target]
    t = d["tags"]
    group = [n for n, v in docs.items() if v["tags"].get("rdname", [None])[0] == target]
    names = [target] + group
    params, seen = [], set()

    def add_params(tg):
        for nm, ds in tg["_params"]:
            for n in nm.split(","):
                pass
            if nm not in seen:
                seen.add(nm); params.append((nm, ds))
    for n in names:
        add_params(docs[n]["tags"])
    used = set()
    for n in names:
        if docs[n]["usage"]:
            a = docs[n]["usage"][docs[n]["usage"].index("(") + 1:-1]
            for x in re.findall(r"(?:^|,)\s*([\w.]+)\s*(?==|,|$)", a):
                used.add(x)
            if "..." in a:
                used.add("...")
    def inherited(name, depth=0):
        tg = docs[name]["tags"]
        out = list(tg["_params"])
        if depth < 5:
            for s2 in tg.get("inheritParams", []):
                out += inherited(s2.strip(), depth + 1)
        return out
    for src in t.get("inheritParams", []):
        for nm, ds in inherited(src.strip()):
            parts = nm.split(",")
            if nm not in seen and any(p in used for p in parts) and not any(p in {q for s in seen for q in s.split(",")} for p in parts):
                seen.add(nm); params.append((nm, ds))
    text = t["_text"]
    title = text[0].strip()
    ps = para(text[1:])
    desc, details = (ps[0], ps[1:]) if ps else ("", [])
    out = ["% Generated by roxygen2: do not edit by hand",
           f"% Please edit documentation in R/{d['file']}", f"\\name{{{target}}}"]
    for n in names:
        out.append(f"\\alias{{{n}}}")
    out.append(f"\\title{{{title}}}")
    us = [wrap_usage(docs[n]["usage"]) for n in names if docs[n]["usage"]]
    if us:
        out.append("\\usage{\n" + "\n\n".join(us) + "\n}")
    if params:
        out.append("\\arguments{")
        for nm, ds in params:
            out.append(f"\\item{{{nm}}}{{{' '.join(ds.split())}}}")
        out.append("}")
    if "return" in t:
        out.append("\\value{\n" + "\n".join(t["return"]).strip() + "\n}")
    out.append("\\description{\n" + desc + "\n}")
    if details:
        out.append("\\details{\n" + "\n\n".join(details) + "\n}")
    if "references" in t:
        out.append("\\references{\n" + "\n".join(t["references"]).strip() + "\n}")
    if "seealso" in t:
        out.append("\\seealso{\n" + "\n".join(t["seealso"]).strip() + "\n}")
    ex = []
    for n in names:
        ex += docs[n]["tags"].get("examples", [])
    if ex:
        out.append("\\examples{\n" + "\n".join(ex).strip("\n") + "\n}")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    docs = all_docs()
    for target in sys.argv[1:]:
        rd = build(target, docs)
        open(os.path.join(ROOT, "man", f"{target}.Rd"), "w", encoding="utf-8").write(rd)
        print("wrote man/%s.Rd" % target)
