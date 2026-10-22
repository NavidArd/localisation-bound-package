"""Compare regenerated tables against the pristine deposited copies.

After re-running any analysis script, the regenerated JSON lands in
`multi_geometry_campaign/data/<phase>/`.  The pristine copy of the same file is in
`results_deposited/<phase>/`, and is never overwritten.  This script diffs the two,
key by key, and reports every difference.

Usage:
    python compare_to_deposited.py                 # every deposited JSON that was re-run
    python compare_to_deposited.py verify_Z2xZ4_1_1.json

Exit status is 0 when everything that was re-run matches, 1 otherwise.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.environ.get("PAPERA_PKG") or os.path.dirname(HERE)
DEP = os.path.join(PKG, "results_deposited")
WORK = os.path.join(PKG, "multi_geometry_campaign", "data")

# deposited subdirectory -> workspace subdirectory, relative to their two roots.
# `two_classes` is the one entry whose workspace copy is NOT under
# multi_geometry_campaign/data/, because that code imports nothing from the campaign and
# lives beside the package root; WORK_OVERRIDE gives its location.
MAP = {
    "su5_ten_bound": "su5_ten_bound",
    "gate_su5_bound": "gate_su5_bound",
    "su5_construction": "su5_construction",
    "gate_su5_construction": "gate_su5_construction",
    "root_lift": "root_lift",
    "two_classes": "two_classes",
}

WORK_OVERRIDE = {
    "two_classes": os.path.join(PKG, "two_classes"),
}

# keys that legitimately differ between runs and are not compared
IGNORE_KEYS = {"seconds", "elapsed", "timestamp", "date", "_meta"}

# Two classes of difference are expected and are reported separately rather than
# silently normalised:
#
#   PATH  the deposited tables record the absolute path of the author's working tree
#         ("F:\physics\...") in their provenance fields; a regenerated copy records the
#         package root instead.  The file NAME, the line number and the quoted text all
#         still have to match -- only the directory prefix may differ.
#
#   NOTE  a few keys are hand-written commentary added to a table after the fact by its
#         author, not emitted by any script.  No re-run can reproduce them.
#
# Every commentary key in this package is written into `d8_classwide.json` and begins with
# `_gate_`.  There are two such blocks: `_gate_notes` (45 keys) and
# `_gate_verdict_vs_phase_7_1` (102 keys), the second of which records the gate's own
# agreement check against a superseded table.  Both are prose, not data, and no claim rests
# on either -- `claims.json` reads the budgets out of the per-geometry keys, which the
# script does emit and which reproduce exactly.  Earlier versions of this list named only
# the first block, so the second was reported as a substantive difference; it is not one.
NOTE_PREFIXES = ("_gate_",)


def path_only(va, vb):
    """True if two values are the same file path under different roots."""
    if not (isinstance(va, str) and isinstance(vb, str)):
        return False
    if "\\" not in va and "/" not in va:
        return False
    return os.path.basename(va.replace("\\", "/")) == os.path.basename(vb.replace("\\", "/"))


def classify(key, va, vb):
    if key.startswith(NOTE_PREFIXES) and vb == "<missing>":
        return "NOTE"
    if path_only(va, vb):
        return "PATH"
    return "REAL"


def flat(node, pre=""):
    """Flatten a JSON tree to {dotted key: scalar}."""
    out = {}
    if isinstance(node, dict):
        for k, v in node.items():
            if k in IGNORE_KEYS:
                continue
            out.update(flat(v, pre + str(k) + "."))
    elif isinstance(node, list):
        out[pre.rstrip(".") + ".__len__"] = len(node)
        for i, v in enumerate(node[:2000]):
            out.update(flat(v, pre + str(i) + "."))
    else:
        out[pre.rstrip(".")] = node
    return out


def work_path(rel):
    """Where the regenerated copy of a deposited table lives."""
    sub = rel.split("/", 1)[0]
    if sub in WORK_OVERRIDE:
        return os.path.join(WORK_OVERRIDE[sub], rel.split("/", 1)[1])
    return os.path.join(WORK, rel)


def compare(rel):
    a = os.path.join(DEP, rel)
    b = work_path(rel)
    if not os.path.exists(b):
        return "NOT-RERUN", []
    try:
        fa, fb = flat(json.load(open(a, encoding="utf-8"))), \
                 flat(json.load(open(b, encoding="utf-8")))
    except Exception as e:
        return "READ-ERROR: %s" % e, []
    diffs = []
    for k in sorted(set(fa) | set(fb)):
        va, vb = fa.get(k, "<missing>"), fb.get(k, "<missing>")
        if va != vb:
            diffs.append((k, va, vb, classify(k, va, vb)))
    if not diffs:
        return "MATCH", diffs
    if any(d[3] == "REAL" for d in diffs):
        return "DIFFERS", diffs
    return "MATCH*", diffs


def main():
    want = sys.argv[1:] or None
    rels = []
    for sub in MAP:
        d = os.path.join(DEP, sub)
        if not os.path.isdir(d):
            continue
        for name in sorted(os.listdir(d)):
            if name.endswith(".json") and (not want or name in want):
                rels.append(os.path.join(sub, name).replace("\\", "/"))

    n_match = n_star = n_diff = n_skip = 0
    bad = []
    for rel in rels:
        status, diffs = compare(rel)
        if status == "MATCH":
            n_match += 1
            print("MATCH      %s" % rel)
        elif status == "MATCH*":
            n_star += 1
            kinds = {}
            for d in diffs:
                kinds[d[3]] = kinds.get(d[3], 0) + 1
            print("MATCH*     %s  (%s -- no substantive difference)"
                  % (rel, ", ".join("%d %s" % (v, k) for k, v in sorted(kinds.items()))))
            for k, va, vb, kind in diffs[:4]:
                print("             [%s] %-46s deposited=%r  regenerated=%r"
                      % (kind, k, va, vb))
            if len(diffs) > 4:
                print("             ... and %d more of the same kinds" % (len(diffs) - 4))
        elif status == "NOT-RERUN":
            n_skip += 1
        else:
            n_diff += 1
            bad.append((rel, status, diffs))
            real = [d for d in diffs if d[3] == "REAL"]
            print("DIFFERS    %s  (%d substantive, %d path/commentary)"
                  % (rel, len(real), len(diffs) - len(real)))
            for k, va, vb, kind in real[:25]:
                print("             %-56s deposited=%r  regenerated=%r" % (k, va, vb))
            if len(real) > 25:
                print("             ... and %d more" % (len(real) - 25))

    print("\n%d matched exactly, %d matched up to absolute paths / non-regenerated "
          "commentary, %d differ substantively, %d not re-run (of %d deposited tables)"
          % (n_match, n_star, n_diff, n_skip, len(rels)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
