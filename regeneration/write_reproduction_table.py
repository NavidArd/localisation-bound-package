"""Emit the run table for REPRODUCTION.md from the artifacts of a run_all.sh pass.

The prose of `REPRODUCTION.md` is written by hand; the table of what ran, for how long,
and with what exit status is generated here from `run_all_status.txt`, so that the record
cannot drift from the run it records.

    python regeneration/write_reproduction_table.py            # print the markdown table
    python regeneration/write_reproduction_table.py --which    # also list which deposited
                                                               # tables this pass rewrote

`--which` compares file modification times in `multi_geometry_campaign/data/` against the
start of the run, and is the honest answer to "which tables did this pass actually
regenerate, as opposed to merely compare?" -- a table can compare equal because it was
never rewritten.
"""
from __future__ import annotations

import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
STATUS = os.path.join(HERE, "run_all_status.txt")
LOGS = os.path.join(HERE, "logs2")
DATA = os.path.join(ROOT, "multi_geometry_campaign", "data")

ROW = re.compile(r"^(\S+)\s*(\S*)\s+exit=(\d+)\s+(\d+)s\s*$")


def human(seconds: int) -> str:
    if seconds < 90:
        return "%d s" % seconds
    if seconds < 5400:
        return "%d min" % round(seconds / 60)
    return "%.1f h" % (seconds / 3600)


def rows():
    out = []
    with open(STATUS, encoding="utf-8") as fh:
        for line in fh:
            m = ROW.match(line.strip())
            if m:
                script, arg, code, secs = m.groups()
                out.append((script, arg, int(code), int(secs)))
    return out


def main() -> int:
    rs = rows()
    if not rs:
        print("no completed rows in %s" % STATUS)
        return 1

    total = sum(r[3] for r in rs)
    nonzero = [r for r in rs if r[2] != 0]

    print("| script | argument | exit | wall time | transcript |")
    print("|---|---|---:|---:|---|")
    for script, arg, code, secs in rs:
        name = os.path.basename(script)
        stem = name[:-3] if name.endswith(".py") else name
        log = "%s%s.log" % (stem, ("_" + arg) if arg else "")
        print("| `%s` | %s | %d | %s | `regeneration/logs2/%s` |"
              % (script, arg or "-", code, human(secs), log))
    print()
    print("%d scripts, %d exited cleanly, %d exited non-zero, %s of wall time in total."
          % (len(rs), len(rs) - len(nonzero), len(nonzero), human(total)))
    if nonzero:
        print("NON-ZERO EXITS: " + ", ".join("%s %s" % (r[0], r[1]) for r in nonzero))

    # The machine-readable twin of the table above, for anyone who would rather parse than
    # read. Regenerated from run_all_status.txt on every invocation, so the two cannot drift.
    import json
    summary = {
        "scripts_run": len(rs),
        "exited_cleanly": len(rs) - len(nonzero),
        "exited_non_zero": len(nonzero),
        "total_seconds": total,
        "transcripts": "regeneration/logs2/",
        "runs": [{"script": a, "argument": b or None, "exit": c, "seconds": d}
                 for a, b, c, d in rs],
    }
    with open(os.path.join(HERE, "run_summary.json"), "w", encoding="utf-8",
              newline="\n") as fh:
        json.dump(summary, fh, indent=1)
        fh.write("\n")
    print("wrote regeneration/run_summary.json")

    if "--which" in sys.argv[1:]:
        start = os.path.getmtime(STATUS) - total - 300
        print("\nDeposited tables this pass actually rewrote (mtime after the run began):")
        n = 0
        for dirpath, _dirs, files in os.walk(DATA):
            for f in sorted(files):
                if not f.endswith(".json"):
                    continue
                p = os.path.join(dirpath, f)
                if os.path.getmtime(p) >= start:
                    n += 1
                    print("   %s" % os.path.relpath(p, ROOT).replace(os.sep, "/"))
        two = os.path.join(ROOT, "two_classes", "paperA_two_classes_enum.json")
        if os.path.exists(two) and os.path.getmtime(two) >= start:
            n += 1
            print("   two_classes/paperA_two_classes_enum.json")
        print("   -- %d tables rewritten" % n)
        print("\nRun began about %s."
              % time.strftime("%Y-%m-%d %H:%M", time.localtime(start)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
