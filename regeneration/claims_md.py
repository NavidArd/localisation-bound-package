"""Render claims.md from claims.json, so the two cannot drift apart.

`claims.json` is the machine-readable claim map; `claims.md` is the same table for a human
reader. Writing the second by hand invites exactly the failure the map exists to prevent --
a number corrected in one file and left standing in the other -- so it is generated.

Run `claims_check.py --write` first: this script prints the `actual` and `status` fields
that resolver puts in, and refuses to run if any row is still PENDING.

    python regeneration/claims_md.py

Sections are ordered by the `paper_section` field, in the order the sections first appear
in claims.json, so adding a row in a new section needs no edit here.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

HEADINGS = {
    "verification": "Scope of the verification",
    "projection norms": "The projection norms",
    "theorem": "The theorem, tested on the full catalogues",
    "mechanism": "The mechanism, tested step by step",
    "Gram ceiling table": "The Gram ceiling against measurement",
    "occupancy": "Occupancy patterns",
    "budgets table": "The budgets table (seven geometries)",
    "construction": "The construction",
    "two-class theorem": "The two-class boundary, re-enumerated independently",
    "denominators": "Denominators and units",
    "signature (scoped)": "The hidden-SU(2) signature, at its stated scope",
    "control (a)": "Control (a) - below threshold",
    "control (b)": "Control (b) - above threshold, doubled (CB1534)",
    "control (c)": "Control (c) - at the ceiling",
    "limitations": "Limitations",
    "exclusions": "Exclusions",
}


def fmt(v):
    if isinstance(v, bool):
        return "**%s**" % v
    if isinstance(v, int):
        return "**%s**" % format(v, ",")
    if isinstance(v, float):
        return "**%s**" % v
    if isinstance(v, dict):
        return "`%s`" % json.dumps(v)
    return "**%s**" % v


def refs(row):
    rs = row.get("derivation") or ([row["source"]] if row.get("source") else [])
    if not rs:
        return "-"
    body = " + ".join("`%s`" % r for r in rs)
    if row.get("reduce"):
        body += " *(%s%s)*" % (row["reduce"],
                              (", key contains `%s`" % row["key_contains"])
                              if row.get("key_contains") else "")
    return body


def main() -> int:
    with open(os.path.join(ROOT, "claims.json"), encoding="utf-8") as fh:
        doc = json.load(fh)
    rows = doc["claims"]

    pending = [r["id"] for r in rows if r.get("status") in (None, "PENDING")]
    if pending:
        print("refusing to render: %d rows unresolved (%s ...). Run claims_check.py --write"
              % (len(pending), ", ".join(pending[:5])))
        return 1

    n = len(rows)
    ok = sum(1 for r in rows if str(r["status"]).startswith("VERIFIED"))
    sup = sum(1 for r in rows if r.get("superseded_by"))
    bad = n - ok

    out = []
    w = out.append
    w("# Claim-to-script map\n")
    w("Every number Paper A prints, the script that regenerates it, the deposited file and")
    w("key that store it, and whether the stored value matches.\n")
    w("**This file is generated from [`claims.json`](claims.json) by")
    w("`regeneration/claims_md.py`, and the `value` column is read out of the deposited")
    w("files by `regeneration/claims_check.py` rather than transcribed.** `status` compares")
    w("that value with the one the paper will print. Rows that are a sum or another")
    w("reduction across the two catalogues show their component keys and the reduction, and")
    w("it is recomputed when this file is built.\n")
    w("| | |")
    w("|---|---:|")
    w("| claims mapped | **%d** |" % n)
    w("| resolved and matching (`VERIFIED`) | **%d** |" % ok)
    w("| mismatched or unresolved | **%d** |" % bad)
    w("| of the matching rows, superseded and **not** printed by the paper | **%d** |" % sup)
    w("")
    if bad == 0:
        w("**All %d rows resolved against the deposited files and matched.**\n" % n)
    else:
        w("**%d rows did NOT resolve or did not match. See `claims_check.py` output.**\n"
          % bad)
    w("A row marked `SUPERSEDED` is still true of the file it came from -- superseded")
    w("numbers are kept and labelled rather than deleted, so that a correction can be")
    w("audited -- but it is not the number the paper prints. The row that replaces it is")
    w("named in the status column.\n")
    w("File paths below are relative to `results_deposited/`; the same tables live in the")
    w("workspace at `multi_geometry_campaign/data/`, where re-runs overwrite them.")
    w("Scripts are relative to `multi_geometry_campaign/code/`, except the two-class")
    w("re-enumeration, which lives in `two_classes/`.\n")
    w("---\n")

    seen = []
    for r in rows:
        s = r.get("paper_section") or "other"
        if s not in seen:
            seen.append(s)

    for sec in seen:
        w("## %s\n" % HEADINGS.get(sec, sec))
        w("| id | claim | value | unit | script | file :: key | status |")
        w("|---|---|---:|---|---|---|---|")
        for r in rows:
            if (r.get("paper_section") or "other") != sec:
                continue
            claim = r["claim"].replace("|", "\\|")
            if r.get("note"):
                claim += " <br/>*%s*" % r["note"].replace("|", "\\|")
            w("| `%s` | %s | %s | %s | `%s` | %s | %s |"
              % (r["id"], claim, fmt(r["actual"]), r.get("unit", ""),
                 r.get("produced_by", ""), refs(r), r["status"]))
        w("")

    w("---\n")
    w("## Regenerating any row\n")
    w("Run the script in the `script` column from its own directory; it overwrites the")
    w("table in `multi_geometry_campaign/data/`. Then\n")
    w("```")
    w("python regeneration/compare_to_deposited.py    # workspace vs pristine reference")
    w("python regeneration/claims_check.py            # every row re-resolved")
    w("python regeneration/claims_md.py               # rebuild this file")
    w("```")

    path = os.path.join(ROOT, "claims.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("wrote %s: %d rows, %d sections, %d superseded" % (path, n, len(seen), sup))
    return 0


if __name__ == "__main__":
    sys.exit(main())
