"""Rebuild manifest.json and CHECKSUMS.sha256 over every file in the package.

This is the last step of assembling or amending the deposit: it walks the tree, hashes
every file, assigns each a role, and writes both integrity files. Deposited so that a
referee can regenerate them and get the same bytes, and so that amending the package is a
scripted step rather than a hand edit.

    python regeneration/build_manifest.py             # write both files
    python regeneration/build_manifest.py --check     # verify, write nothing (exit 1 on drift)

What is skipped
---------------
`__pycache__/` and `*.pyc` only. Byte-compiled Python depends on the interpreter build and
on absolute paths, so including it would make CHECKSUMS.sha256 unverifiable on any other
machine. `.gitignore` lists the same exclusion. Everything else in the tree is hashed --
including manifest.json's own companion CHECKSUMS.sha256, which is written first so that
the manifest can cover it.

manifest.json cannot cover itself, and does not claim to.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
MANIFEST = os.path.join(ROOT, "manifest.json")
CHECKSUMS = os.path.join(ROOT, "CHECKSUMS.sha256")

SKIP_DIRS = {"__pycache__"}
SKIP_SUFFIX = (".pyc", ".pyo")
SELF = {"manifest.json", "CHECKSUMS.sha256"}

# Longest matching prefix wins; the final fallbacks are by extension.
ROLES = [
    ("orbifolder-1.2.1/", "upstream input, GPL v3 (redistributed unmodified)"),
    ("environment/orbifolder_patches/", "geometry patches (NOT on the critical path)"),
    ("environment/orbifolder_drivers/", "orbifolder run driver"),
    ("environment/", "environment specification"),
    ("multi_geometry_campaign/data/literature/",
     "upstream catalogue, arXiv:2003.01732 ancillary"),
    ("multi_geometry_campaign/data/root_lift/", "parsed table (workspace copy; re-runs overwrite)"),
    ("multi_geometry_campaign/code/", "analysis code"),
    ("multi_geometry_campaign/data/", "parsed table (workspace copy; re-runs overwrite)"),
    ("alphabet_survey/code/", "analysis code (certified engine, imported by the parsers)"),
    ("alphabet_survey/data/", "engine input table"),
    ("results_deposited/", "pristine reference copy of a parsed table"),
    ("regeneration/logs", "transcript of a reproduction run made for this package"),
    ("regeneration/", "reproduction record and tooling"),
    ("two_classes/", "independent re-enumeration of the two-class boundary"),
    ("construction/", "the constructed model: deck, vectors, raw output"),
    ("controls/", "positive controls"),
    ("exclusions/", "exclusion and recovery log"),
]

LICENCES = {
    "orbifolder-1.2.1/": "GPL v3 (see orbifolder-1.2.1/NOTICE.md)",
    "environment/orbifolder_patches/patched_geometry_files/":
        "GPL v3, modified (see environment/orbifolder_patches/NOTICE.md)",
    "multi_geometry_campaign/data/literature/parr2020_full/":
        "arXiv ancillary files of arXiv:2003.01732 (see NOTICE.md there)",
    "code (*.py)": "MIT (LICENSE-CODE.txt)",
    "everything else": "CC BY 4.0 (LICENSE-DATA.txt)",
}


SUPERSEDED = ("multi_geometry_campaign/data/representatives/brothers/superseded_v2.2/",
              "superseded v2.2 file, kept and labelled (see NOTICE.md there); enters no number")


def role_for(rel: str) -> str:
    if rel.startswith(SUPERSEDED[0]):
        return SUPERSEDED[1]
    if rel.endswith(".txt.gz"):
        return "raw orbifolder output (gzipped)"
    for prefix, role in ROLES:
        if rel.startswith(prefix):
            return role
    return "package documentation"


def walk():
    """Every deposited file, as a path relative to the package root, POSIX-separated."""
    out = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS)
        for name in sorted(filenames):
            if name.endswith(SKIP_SUFFIX):
                continue
            full = os.path.join(dirpath, name)
            rel = os.path.relpath(full, ROOT).replace(os.sep, "/")
            out.append(rel)
    return sorted(out)


def sha256(path: str) -> tuple[str, int]:
    h = hashlib.sha256()
    n = 0
    with open(path, "rb") as fh:
        while True:
            block = fh.read(1 << 20)
            if not block:
                break
            n += len(block)
            h.update(block)
    return h.hexdigest(), n


def build(check: bool) -> int:
    # CHECKSUMS first, over everything except the two integrity files themselves,
    # then the manifest, which covers CHECKSUMS too.
    rels = [r for r in walk() if r not in SELF]

    entries = []
    lines = []
    for rel in rels:
        digest, nbytes = sha256(os.path.join(ROOT, rel))
        entries.append({"path": rel, "bytes": nbytes, "sha256": digest,
                        "role": role_for(rel)})
        lines.append("%s *%s\n" % (digest, rel))
    checksum_text = "".join(lines)

    if check:
        try:
            with open(CHECKSUMS, encoding="utf-8") as fh:
                on_disk = fh.read()
        except FileNotFoundError:
            print("CHECKSUMS.sha256 missing")
            return 1
        if on_disk != checksum_text:
            have = dict(l.split(" *", 1) for l in on_disk.splitlines() if " *" in l)
            want = dict(l.split(" *", 1) for l in checksum_text.splitlines() if " *" in l)
            for p in sorted(set(want) | set(have)):
                pass
            missing = sorted(set(want.values()) - set(have.values()))
            extra = sorted(set(have.values()) - set(want.values()))
            print("CHECKSUMS.sha256 is out of date: %d new/changed, %d stale"
                  % (len(missing), len(extra)))
            for p in (missing + extra)[:20]:
                print("   ", p)
            return 1
        print("CHECKSUMS.sha256 matches all %d files" % len(entries))
        return 0

    with open(CHECKSUMS, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(checksum_text)

    digest, nbytes = sha256(CHECKSUMS)
    entries.append({"path": "CHECKSUMS.sha256", "bytes": nbytes, "sha256": digest,
                    "role": "package documentation"})
    entries.sort(key=lambda e: e["path"])

    total = sum(e["bytes"] for e in entries)
    summary: dict[str, dict] = {}
    for e in entries:
        s = summary.setdefault(e["role"], {"files": 0, "mb": 0.0})
        s["files"] += 1
        s["mb"] += e["bytes"] / 1e6
    for s in summary.values():
        s["mb"] = round(s["mb"], 2)

    doc = {
        "package": "paperA_repro_package",
        "version": "2.2",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "doi": "NOT YET ASSIGNED (placeholder -- see CITATION.cff)",
        "file_count": len(entries),
        "total_bytes": total,
        "total_mb": round(total / 1e6, 1),
        "note": ("Every file in the deposit except manifest.json itself, which cannot "
                 "cover its own hash. __pycache__ and *.pyc are excluded and are listed "
                 "in .gitignore. CHECKSUMS.sha256 carries the same hashes in sha256sum "
                 "format; .gitattributes disables end-of-line conversion so both verify "
                 "on Windows clones."),
        "summary_by_role": dict(sorted(summary.items())),
        "licences": LICENCES,
        "files": entries,
    }
    with open(MANIFEST, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(doc, fh, indent=1)
        fh.write("\n")

    print("manifest.json  %d files, %s bytes (%.1f MB)"
          % (len(entries), format(total, ","), total / 1e6))
    print("CHECKSUMS.sha256  %d lines" % len(lines))
    for role, s in sorted(summary.items(), key=lambda kv: -kv[1]["mb"]):
        print("   %-58s %4d  %8.2f MB" % (role, s["files"], s["mb"]))
    return 0


if __name__ == "__main__":
    sys.exit(build(check="--check" in sys.argv[1:]))
