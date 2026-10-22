"""Resolve every row of claims.json against the deposited tables and re-derive its status.

`claims.json` states, for each number the paper prints, where that number lives on disk.
This script reads it back out and compares. It is the script the claim map's `note` field
names, and it is what makes the `status` column a measurement rather than a transcription.

How a row is resolved
---------------------
A row carries either

  "source":     "<file>::<dotted.key.path>"                        one location, or
  "derivation": ["<file>::<key>", "<file>::<key>", ...]            several, summed,

where `<file>` is relative to `results_deposited/` and `<dotted.key.path>` walks the JSON
tree. A key segment may itself contain a slash (fractions such as `6/5` are dictionary
keys throughout the deposited tables), so the walk tries progressively longer joins of the
remaining segments before giving up -- `w2_by_fam.10.6/5` resolves correctly.

A row may also carry a `reduce` field naming the one operation applied to what was read,
so that no reduction is implicit:

  "len"            the number printed is the length of the list stored there
  "max_int_key"    ... the largest integer key of the multiplicity histogram stored there
  "sum_values"     ... the sum of one histogram's values
  "sum_hist_values"  ... the sum of the values of each of several histograms
  "difference"     ... the first value minus the second
  "ratio3"         ... the first divided by the second, to three decimals
  "sum_matching"   ... the sum of the histogram values whose key contains the
                       row's `key_contains` string

A row may carry `superseded_by: "<id>"`. Such a row is still resolved and still compared --
a superseded number must remain true of the file it came from -- but it is reported
separately, and it is not the number the paper prints.

Rows with neither `source` nor `derivation` are reported as UNRESOLVED rather than
silently skipped.

Usage
-----
    python regeneration/claims_check.py               # check, print a summary
    python regeneration/claims_check.py --write       # also write the resolved
                                                      # actual/status back into claims.json

Exit code is 0 when every resolvable row matches, 1 otherwise.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DEPOSITED = os.path.join(ROOT, "results_deposited")
CLAIMS = os.path.join(ROOT, "claims.json")

_cache: dict[str, object] = {}


def load(rel: str):
    """Load one deposited JSON table, by path relative to results_deposited/."""
    if rel not in _cache:
        with open(os.path.join(DEPOSITED, rel.replace("/", os.sep)), encoding="utf-8") as fh:
            _cache[rel] = json.load(fh)
    return _cache[rel]


class Unresolved(Exception):
    pass


def walk(node, segments):
    """Walk a JSON tree by dotted segments, allowing a key to contain '.' or '/'.

    At each step the longest run of remaining segments that names an existing key wins;
    this is what lets `w2_by_fam.10.6/5` and `B_by_sector.(1, 3).7/16` resolve."""
    if not segments:
        return node
    if isinstance(node, list):
        try:
            return walk(node[int(segments[0])], segments[1:])
        except (ValueError, IndexError):
            raise Unresolved("list index %r" % (segments[0],))
    if not isinstance(node, dict):
        raise Unresolved("not a container at %r" % (segments[0],))
    for n in range(len(segments), 0, -1):
        key = ".".join(segments[:n])
        if key in node:
            return walk(node[key], segments[n:])
    raise Unresolved("no key for %r" % (segments[0],))


def resolve_one(ref: str):
    """Resolve one '<file>::<dotted.key>' reference to its stored value."""
    if "::" not in ref:
        raise Unresolved("malformed reference %r" % ref)
    rel, keypath = ref.split("::", 1)
    return walk(load(rel), keypath.split("."))


def apply_reduce(values, how: str, key_contains: str | None = None):
    """Apply the row's declared reduction to the list of values read."""
    if how == "sum_matching":
        if key_contains is None:
            raise Unresolved("reduce sum_matching needs a key_contains")
        return sum(v for hist in values
                   for k, v in hist.items() if key_contains in k)
    if how == "len":
        (v,) = values
        return len(v)
    if how == "max_int_key":
        (v,) = values
        return max(int(k) for k in v)
    if how == "sum_values":
        (v,) = values
        return sum(v.values())
    if how == "sum_hist_values":
        return sum(sum(v.values()) for v in values)
    if how == "difference":
        a, b = values
        return a - b
    if how == "ratio3":
        a, b = values
        return round(a / b, 3)
    raise Unresolved("unknown reduce %r" % how)


def resolve(row) -> tuple[object, str]:
    """Return (actual, how) for one claim row."""
    refs = row.get("derivation") or ([row["source"]] if row.get("source") else [])
    if not refs:
        raise Unresolved("row carries neither source nor derivation")
    values = [resolve_one(r) for r in refs]
    red = row.get("reduce")
    if red:
        return (apply_reduce(values, red, row.get("key_contains")),
                "%s over %d" % (red, len(values)))
    if len(values) == 1:
        return values[0], "read"
    if all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in values):
        return sum(values), "sum of %d" % len(values)
    raise Unresolved("cannot combine %d non-numeric values" % len(values))


def main(write: bool = False) -> int:
    with open(CLAIMS, encoding="utf-8") as fh:
        doc = json.load(fh)

    counts = {"VERIFIED": 0, "MISMATCH": 0, "UNRESOLVED": 0}
    n_superseded = 0
    problems = []

    for row in doc["claims"]:
        try:
            actual, how = resolve(row)
        except (Unresolved, FileNotFoundError, KeyError) as exc:
            row["actual"] = None
            row["status"] = "UNRESOLVED"
            counts["UNRESOLVED"] += 1
            problems.append("%-10s UNRESOLVED  %s" % (row["id"], exc))
            continue
        row["actual"] = actual
        if actual == row["expected"]:
            if row.get("superseded_by"):
                row["status"] = "VERIFIED (SUPERSEDED by %s)" % row["superseded_by"]
                n_superseded += 1
            else:
                row["status"] = "VERIFIED"
            counts["VERIFIED"] += 1
        else:
            row["status"] = "MISMATCH"
            counts["MISMATCH"] += 1
            problems.append("%-10s MISMATCH    expected %r, deposited %r  (%s)"
                            % (row["id"], row["expected"], actual, how))

    print("claims resolved against %s" % os.path.relpath(DEPOSITED, ROOT))
    for k in ("VERIFIED", "MISMATCH", "UNRESOLVED"):
        print("  %-11s %4d" % (k, counts[k]))
    print("  %-11s %4d" % ("total", len(doc["claims"])))
    print("  (of the verified rows, %d are marked superseded and are not printed by the "
          "paper)" % n_superseded)
    if problems:
        print("\n".join([""] + problems))

    if write:
        with open(CLAIMS, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, indent=1, ensure_ascii=False)
            fh.write("\n")
        print("\nwrote resolved actual/status back into claims.json")

    return 0 if counts["MISMATCH"] == 0 and counts["UNRESOLVED"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main(write="--write" in sys.argv[1:]))
