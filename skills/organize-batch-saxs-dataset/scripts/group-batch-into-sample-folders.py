#!/usr/bin/env python3
"""Group a flat beamline batch directory into one folder per sample.

A series is named ``<name>_<run:04d>_<frame:05d>.<ext>`` (BL19U2 batch mode).
Sample identity is read from the ``Description:`` header of ``<series>_00001.txt``;
series whose description exactly matches a name in ``--background`` are buffer
measurements and are COPIED into the folder of the nearest sample(s) --- both
neighbours when the buffer is bracketed between two samples.

Defaults target a dry run; pass --apply to actually copy, then every copied file
is verified against its source by sha256 and the folder's file SET is compared
with the plan.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import sys
from collections import defaultdict

DEFAULT_BACKGROUND = "pb7,ddh2o"
DEFAULT_EXT = "tif"
SERIES_RE_TEMPLATE = r"^(?P<name>.+)_(?P<run>\d{{4}})_(?P<frame>\d{{5}})\.{ext}$"


class Collision(Exception):
    """A planned folder name collides with an existing entry (case-insensitively)."""


def sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_description(path: str) -> str | None:
    """Return the ``Description:`` header of a frame .txt file, if present."""
    if not os.path.exists(path):
        return None
    with open(path, errors="replace") as handle:
        for line in handle:
            if line.startswith("Description:"):
                return line.split(":", 1)[1].strip()
    return None


def collect(src: str, ext: str):
    """Return (run_index, description, files_by_series, orphans)."""
    names = sorted(n for n in os.listdir(src) if not n.startswith("."))
    pattern = re.compile(SERIES_RE_TEMPLATE.format(ext=re.escape(ext)))

    run: dict[str, int] = {}
    for name in names:
        match = pattern.match(name)
        if match:
            run[f"{match['name']}_{match['run']}"] = int(match["run"])
    if not run:
        raise SystemExit(f"no series matching '*_0000_00000.{ext}' found in {src}")

    ids = sorted(run, key=len, reverse=True)
    files: dict[str, list[str]] = defaultdict(list)
    orphans: list[str] = []
    for name in names:
        for sid in ids:  # pass 1: series id at the start of the name
            if name == sid or name.startswith(sid + "_") or name.startswith(sid + "."):
                files[sid].append(name)
                break
        else:  # pass 2: derived artefacts such as S_A_<series>_00001.out / _01.log
            for sid in ids:
                if sid + "_" in name:
                    files[sid].append(name)
                    break
            else:
                orphans.append(name)

    desc = {sid: read_description(os.path.join(src, f"{sid}_00001.txt")) for sid in run}
    return run, desc, files, orphans


def build_groups(run, desc, backgrounds):
    """Split samples from backgrounds; assign each background to its nearest sample(s)."""
    samples = [s for s in run if (desc.get(s) or "").strip().casefold() not in backgrounds]
    buffers = [s for s in run if (desc.get(s) or "").strip().casefold() in backgrounds]
    groups: dict[str, list[str]] = {s: [s] for s in samples}
    for buf in buffers:
        distance = {s: abs(run[s] - run[buf]) for s in samples}
        nearest = min(distance.values()) if distance else None
        for sample, dist in distance.items():
            if dist == nearest:  # ties go to every equally-near sample
                groups[sample].append(buf)
    for sample in groups:
        groups[sample] = sorted(set(groups[sample]), key=lambda s: run[s])
    return samples, buffers, groups


def folder_name(sid: str, desc) -> str:
    return (desc.get(sid) or sid.rsplit("_", 1)[0]).strip()


def check_collisions(dst: str, planned: dict[str, str]) -> None:
    """Refuse case-variant folder names: on APFS/NTFS they silently MERGE."""
    existing = {n.casefold(): n for n in os.listdir(dst)} if os.path.isdir(dst) else {}
    problems = []
    for folder, sample in planned.items():
        probe = folder.casefold()
        if probe in existing and existing[probe] != folder:
            problems.append(
                f"  {folder!r} (sample {sample}) vs existing {existing[probe]!r}"
            )
    if problems:
        raise Collision(
            "destination folder names collide case-insensitively with existing "
            "entries; copying would merge two experiments into one folder:\n"
            + "\n".join(problems)
            + "\nRename the new folders (or confirm the merge with the user) and retry."
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--src", required=True, help="flat batch directory to read (files)")
    parser.add_argument("--dst", required=True, help="directory to create one folder per sample in")
    parser.add_argument(
        "--background",
        default=DEFAULT_BACKGROUND,
        help="comma-separated Description values meaning 'background, not a sample' (no spaces)",
    )
    parser.add_argument("--ext", default=DEFAULT_EXT, help="frame extension that defines a series")
    parser.add_argument(
        "--apply", action="store_true", help="actually copy; default is a dry run (plan only)"
    )
    parser.add_argument(
        "--yes", action="store_true", help="skip the interactive confirmation before copying"
    )
    args = parser.parse_args()

    backgrounds = {b.strip().casefold() for b in args.background.split(",") if b.strip()}
    run, desc, files, orphans = collect(args.src, args.ext)
    samples, buffers, groups = build_groups(run, desc, backgrounds)

    planned = {folder_name(s, desc): s for s in samples}
    if len(planned) != len(samples):
        raise SystemExit(
            "two samples share a Description; folder names would collide -- resolve manually"
        )

    print(
        f"series={len(run)}  samples={len(samples)}  background={len(buffers)}  "
        f"files={sum(len(v) for v in files.values())}"
    )
    if orphans:
        print(f"UNASSIGNED FILES (fix the matcher before copying): {orphans}")

    placements = defaultdict(int)
    for sample in groups:
        for sid in groups[sample]:
            placements[sid] += 1
    shared = sorted(s for s in buffers if placements[s] == 2)

    total = 0
    print(f"\n{'sample folder':<20} {'run':>4}  {'files':>5}  contains")
    for sample in sorted(groups, key=lambda s: run[s]):
        members = groups[sample]
        count = sum(len(files[s]) for s in members)
        total += count
        inside = ", ".join(
            f"{s}{' [BG]' if s in buffers else ''}" for s in members
        )
        print(f"{folder_name(sample, desc):<20} {run[sample]:>4}  {count:>5}  {inside}")
    print(f"\ntotal files after copying (duplicates counted): {total}")
    print(f"backgrounds landing in 2 folders: {len(shared)}")
    single = sorted(
        (f"{s} [BG] -> {folder_name(next(x for x in groups if s in groups[x]), desc)}"
         for s in buffers if placements[s] == 1)
    )
    print(f"backgrounds landing in 1 folder ({len(single)}): {single}")
    if not args.apply:
        print("\ndry run only -- re-run with --apply to copy")
        return 0
    if orphans:
        return 2

    os.makedirs(args.dst, exist_ok=True)
    try:
        check_collisions(args.dst, planned)
    except Collision as exc:
        print(f"\nCOLLISION: {exc}", file=sys.stderr)
        return 3

    if not args.yes:
        reply = input("\nproceed with the copy? [y/N] ").strip().casefold()
        if reply not in {"y", "yes"}:
            print("aborted")
            return 1

    source_hash = {
        name: sha256(os.path.join(args.src, name)) for name in os.listdir(args.src)
        if not name.startswith(".")
    }
    copied = 0
    for sample in sorted(groups, key=lambda s: run[s]):
        folder = os.path.join(args.dst, folder_name(sample, desc))
        os.makedirs(folder, exist_ok=True)
        for sid in groups[sample]:
            for name in files[sid]:
                target = os.path.join(folder, name)
                if not os.path.exists(target):
                    shutil.copy2(os.path.join(args.src, name), target)
                    copied += 1

    bad = []
    for sample in sorted(groups, key=lambda s: run[s]):
        folder = os.path.join(args.dst, folder_name(sample, desc))
        want = sorted(f for sid in groups[sample] for f in files[sid])
        have = sorted(f for f in os.listdir(folder) if not f.startswith("."))
        extra = sorted(set(have) - set(want))
        missing = sorted(set(want) - set(have))
        corrupt = [f for f in have if source_hash.get(f) and sha256(os.path.join(folder, f)) != source_hash[f]]
        if extra or missing or corrupt:
            bad.append((folder, extra[:3], missing[:3], corrupt[:3]))
    print(f"\ncopied {copied} file(s); verification: {'OK' if not bad else 'FAILED'}")
    for entry in bad:
        print(f"  {entry}")
    return 0 if not bad else 4


if __name__ == "__main__":
    raise SystemExit(main())
