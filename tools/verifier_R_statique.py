"""Static checks of the R part of canofac, for an environment without R.

It does not replace R CMD check; it catches, without running R, the errors that R CMD check would report most
often: unbalanced brackets, non-ASCII characters in R code, exports without definition or documentation,
usage sections that do not match the function arguments, internal functions or labels that do not exist,
packages used but not declared in DESCRIPTION.

Usage (from the package root): python3 tools/verifier_R_statique.py
Exit code 0 if no error; warnings are printed but do not fail.
"""
import csv
import glob
import os
import re
import sys

errors, warnings = [], []


def strip_code(src):
    """Remove comments and string contents from R code, keeping line structure."""
    out, i, n = [], 0, len(src)
    while i < n:
        c = src[i]
        if c == "#":
            while i < n and src[i] != "\n":
                i += 1
            continue
        if c in "\"'`":
            q = c
            out.append(q)
            i += 1
            while i < n and src[i] != q:
                if src[i] == "\\":
                    i += 2
                    continue
                out.append("\n" if src[i] == "\n" else " ")
                i += 1
            out.append(q)
            i += 1
            continue
        if c == "r" and i + 2 < n and src[i + 1] in "\"'" and src[i + 2] in "([{":
            # raw string r"(...)"
            q, o = src[i + 1], src[i + 2]
            cl = {"(": ")", "[": "]", "{": "}"}[o]
            end = src.find(cl + q, i + 3)
            out.append('""')
            i = end + 2
            continue
        out.append(c)
        i += 1
    return "".join(out)


def check_brackets(path, code):
    stack, pairs = [], {")": "(", "]": "[", "}": "{"}
    line = 1
    for ch in code:
        if ch == "\n":
            line += 1
        elif ch in "([{":
            stack.append((ch, line))
        elif ch in ")]}":
            if not stack or stack[-1][0] != pairs[ch]:
                errors.append(f"{path}:{line}: unexpected '{ch}'")
                return
            stack.pop()
    for ch, ln in stack:
        errors.append(f"{path}:{ln}: '{ch}' never closed")


def split_args(s):
    """Split a formal argument list at top-level commas."""
    parts, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        parts.append(cur)
    return [p.split("=")[0].strip() for p in parts]


def formals_of(src_raw):
    """Map function name -> list of argument names (definitions name <- function(...))."""
    res = {}
    for m in re.finditer(r"^([A-Za-z.][\w.]*)\s*(?:<-|=)\s*function\s*\(", src_raw, re.M):
        i, depth = m.end(), 1
        while depth:
            depth += {"(": 1, ")": -1}.get(src_raw[i], 0)
            i += 1
        res[m.group(1)] = split_args(src_raw[m.end():i - 1])
    return res


# ---------------------------------------------------------------- R files
rfiles = sorted(glob.glob("R/*.R"))
defs, calls_internal, used_pkgs, tr_keys = {}, set(), set(), set()
for f in rfiles:
    raw = open(f, encoding="utf-8").read()
    for ln, l in enumerate(raw.splitlines(), 1):
        if any(ord(c) > 127 for c in l):
            where = "roxygen comment" if l.lstrip().startswith("#'") else "code"
            (warnings if where == "roxygen comment" else errors).append(
                f"{f}:{ln}: non-ASCII character in {where}" + ("" if where == "roxygen comment" else
                " (R CMD check: use \\uxxxx escapes)"))
    code = strip_code(raw)
    check_brackets(f, code)
    defs.update({k: (f, v) for k, v in formals_of(raw).items()})
    calls_internal |= set(re.findall(r"(?<![\w.])(\.[A-Za-z][\w.]*)\s*\(", code))
    used_pkgs |= set(re.findall(r"\b([A-Za-z][\w.]*)::", code))
    used_pkgs |= set(re.findall(r"requireNamespace\(\s*[\"']([\w.]+)", raw))
    tr_keys |= set(re.findall(r"\.tr\(\s*\"([\w.]+)\"", raw))

base_internal = {".Machine", ".libPaths", ".Platform", ".Call", ".External", ".Internal", ".Primitive",
                 ".subset2", ".subset", ".Deprecated", ".Defunct", ".packageName", ".GlobalEnv"}
for fn in sorted(calls_internal - set(defs) - base_internal):
    errors.append(f"internal function {fn}() called but not defined in R/")

# ---------------------------------------------------------------- NAMESPACE
ns = open("NAMESPACE", encoding="utf-8").read()
exports = re.findall(r"^export\(([^)]+)\)", ns, re.M)
s3 = re.findall(r"^S3method\(([^,]+),\s*([^)]+)\)", ns, re.M)
imports = set(re.findall(r"^import(?:From)?\(([\w.]+)", ns, re.M))
for e in exports:
    if e not in defs:
        errors.append(f"NAMESPACE exports {e} but no definition in R/")
for g, c in s3:
    if f"{g}.{c}" not in defs:
        errors.append(f"NAMESPACE registers S3method({g}, {c}) but {g}.{c} is not defined")
for name, (f, _) in defs.items():
    m = re.match(r"(print|plot|summary|format|as\.data\.frame|coef|predict)\.(\w+)$", name)
    if m and (m.group(1), m.group(2)) not in s3:
        warnings.append(f"{name} ({f}) looks like an S3 method but is not registered in NAMESPACE")

# ---------------------------------------------------------------- DESCRIPTION
desc = open("DESCRIPTION", encoding="utf-8").read()
declared = set(re.findall(r"([A-Za-z][\w.]+)\s*(?:\(|,|\n)", " ".join(
    re.findall(r"^(?:Imports|Suggests|Depends):(.*?)(?=^\S)", desc, re.M | re.S))))
basepk = {"base", "R"}
for p in sorted(used_pkgs - declared - basepk - {"canofac"}):
    errors.append(f"package {p} used in R/ but not in Imports/Suggests of DESCRIPTION")
for line in desc.splitlines():
    if any(ord(c) > 127 for c in line):
        warnings.append("DESCRIPTION contains non-ASCII characters (Encoding: UTF-8 is declared)")
        break

# ---------------------------------------------------------------- Rd files
rd_alias, rd_usage = {}, {}
for rd in sorted(glob.glob("man/*.Rd")):
    t = open(rd, encoding="utf-8").read()
    for a in re.findall(r"\\alias\{([^}]+)\}", t):
        rd_alias[a] = rd
    u = re.search(r"\\usage\{(.*?)\n\}", t, re.S)
    if u:
        for m in re.finditer(r"^(?:\\method\{(\w+)\}\{(\w+)\}|([A-Za-z.][\w.]*))\(", u.group(1), re.M):
            name = f"{m.group(1)}.{m.group(2)}" if m.group(1) else m.group(3)
            i, depth, s = m.end(), 1, u.group(1)
            while depth and i < len(s):
                depth += {"(": 1, ")": -1}.get(s[i], 0)
                i += 1
            rd_usage[name] = (rd, split_args(s[m.end():i - 1]))
    args_doc = set()
    for a in re.findall(r"\\item\{([^}]+)\}\{", t.split("\\value")[0]):
        args_doc |= {x.strip() for x in a.split(",")}
    if u:
        for name, (rdf, args) in list(rd_usage.items()):
            if rdf == rd:
                for a in args:
                    if a and a not in args_doc and a != "...":
                        warnings.append(f"{rd}: argument '{a}' of {name} in \\usage but not in \\arguments")
for e in exports:
    if e not in rd_alias:
        errors.append(f"exported {e} has no Rd page (no \\alias)")
for name, (rd, args) in rd_usage.items():
    if name in defs:
        real = defs[name][1]
        if args != real:
            errors.append(f"{rd}: \\usage of {name} {args} differs from the code {real}")

# ---------------------------------------------------------------- labels
lab_file = "inst/shared/labels.csv"
if os.path.exists(lab_file):
    keys = {r[0] for r in csv.reader(open(lab_file, encoding="utf-8"))}
    for k in sorted(tr_keys - keys):
        errors.append(f"label key '{k}' used by .tr() but absent from {lab_file}")

# ---------------------------------------------------------------- tests, scripts, vignettes
for f in glob.glob("tests/testthat/*.R") + glob.glob("inst/scripts/*.R") + glob.glob("vignettes/*.Rmd"):
    raw = open(f, encoding="utf-8").read()
    check_brackets(f, strip_code(raw if f.endswith(".R") else "\n".join(
        re.findall(r"```\{r[^}]*\}\n(.*?)```", raw, re.S))))
    for fn in set(re.findall(r"\b(cf_\w+|canofac)\s*\(", raw)):
        if fn not in defs:
            errors.append(f"{f}: calls {fn}() which is not defined in R/")

print(f"R files: {len(rfiles)}, functions: {len(defs)}, exports: {len(exports)}, S3 methods: {len(s3)}, "
      f"Rd pages: {len(glob.glob('man/*.Rd'))}, label keys used: {len(tr_keys)}")
for w in warnings:
    print("WARNING", w)
for e in errors:
    print("ERROR  ", e)
print(f"{len(errors)} error(s), {len(warnings)} warning(s)")
sys.exit(1 if errors else 0)
