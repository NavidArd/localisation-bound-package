#!/bin/bash
# Are the duplicated `W = ` lines within one model IDENTICAL?
# `create coupling(F F F)` forks a child that inherits the redirected stdin, so in a
# MULTI-model run the child re-executes the tail of the command stream and the
# superpotential is printed twice.  Harmless iff the two prints agree exactly.
python3 - "$1" <<'PY'
import re, sys
t = open(sys.argv[1], errors="replace").read()
# segment by model using the couplings prompt
parts = re.split(r"/(MSSM\d+)/couplings>", t)
cur = {}
for i in range(1, len(parts), 2):
    lab, body = parts[i], parts[i+1]
    m = re.match(r"\s*\n\s+W = (.*?)(?:\n\s*\n|\n/)", body, re.S)
    if m:
        cur.setdefault(lab, []).append(" ".join(m.group(1).split()))
dups = {k: v for k, v in cur.items() if len(v) > 1}
diff = {k: v for k, v in dups.items() if len(set(v)) > 1}
print("models seen: %d | models printing W more than once: %d | DIFFERING duplicates: %d"
      % (len(cur), len(dups), len(diff)))
if diff:
    for k in list(diff)[:3]:
        print("  DIFFER", k, [x[:80] for x in diff[k]])
PY
