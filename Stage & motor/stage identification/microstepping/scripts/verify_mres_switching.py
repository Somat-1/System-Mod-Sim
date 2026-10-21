#!/usr/bin/env python3
"""Did the microstep size really change when MRES was written?

    python3 scripts/verify_mres_switching.py data/terminal_*.log

A CHOPCONF readback only proves the register holds the code. It does not prove
the driver is stepping on the grid that code asks for - if GCONF.mstep_reg_select
were clear the pins would govern and the readback would still look right, and if
MicroPlyer (intpol) were on every size would interpolate to 1/256.

MSCNT is the driver's own index into the sine table: 1024 entries over one
electrical period, 256 per full step. It can only ever land on the grid the
driver is actually stepping on, so it reports the real microstep size:

    expected spacing = 256 / mres    1/16 -> 16, 1/32 -> 8, 1/64 -> 4, 1/256 -> 1

The test is on RESIDUES, not on the GCD of differences. Within a block every
MSCNT value should share one residue mod (256/mres) - one coset of the correct
grid. A handful of samples taken across an MRES change sit on a different
coset, and those few are enough to drag a GCD down to half the spacing and fake
a failure. Counting residues separates "the grid is wrong" from "a few samples
straddle a phase change", which a GCD cannot do.

Attribution is exact rather than inferred: every `# S,` sample carries the plan
row it belongs to, and every `# MOVE,` line gives that row's MRES. Samples with
row = -1 are interactive (homing, warm-up) and are left out - the warm-up menu
could change MRES without a plan row recording it.

Both telemetry formats are read:

    # S,ms,pos_u256,enc,mscnt,row      mres_test  (5 fields)
    # S,ms,pos_u256,mscnt,row          settling_test (4 fields)
"""

from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

SAMPLE = re.compile(r"# S,(-?\d+(?:,-?\d+){3,4})")
MOVE = re.compile(r"# MOVE,(\d+),.*?mres=(\d+)")
MRESLINE = re.compile(r"# MRES,(\d+)")

# The share of samples that must sit on one coset for the grid to be sound.
# The rest are the samples straddling an MRES change.
DOMINANT_MIN = 0.95


def parse(path: Path):
    """Return {mres: [mscnt, ...]} attributed by plan row."""
    row_mres: dict[int, int] = {}
    by_mres: dict[int, list[int]] = defaultdict(list)
    moved: dict[int, list[int]] = defaultdict(list)
    for line in path.open(encoding="utf-8", errors="replace"):
        m = MOVE.search(line)
        if m:
            row_mres[int(m.group(1))] = int(m.group(2))
            continue
        s = SAMPLE.search(line)
        if not s:
            continue
        f = [int(x) for x in s.group(1).split(",")]
        if len(f) == 5:            # ms, pos, enc, mscnt, row
            pos, mscnt, row = f[1], f[3], f[4]
        elif len(f) == 4:          # ms, pos, mscnt, row
            pos, mscnt, row = f[1], f[2], f[3]
        else:
            continue
        if row in row_mres:
            by_mres[row_mres[row]].append(mscnt)
            moved[row_mres[row]].append(pos)
    return by_mres, moved


def report(path: Path) -> bool:
    by_mres, moved = parse(path)
    if not by_mres:
        print(f"\n=== {path.name} ===\n  no plan rows in this log "
              "(homing / warm-up only)")
        return True
    print(f"\n=== {path.name} ===")
    print(f"{'MRES':>7} {'samples':>8} {'moved_mm':>9} {'spacing':>8} "
          f"{'cosets':>7} {'dominant':>9}  verdict")
    ok_all = True
    for mres in sorted(by_mres):
        mscnt = by_mres[mres]
        pos = moved[mres]
        expect = 256 // mres
        span = (max(pos) - min(pos)) / 25600.0
        if expect == 1:
            # 1/256 can land anywhere; there is no grid to check.
            print(f"  1/{mres:<4} {len(mscnt):8d} {span:9.3f} {expect:8d} "
                  f"{'-':>7} {'-':>9}  n/a (finest)")
            continue
        res = Counter(v % expect for v in mscnt)
        top, count = res.most_common(1)[0]
        share = count / len(mscnt)
        ok = share >= DOMINANT_MIN
        ok_all &= ok
        print(f"  1/{mres:<4} {len(mscnt):8d} {span:9.3f} {expect:8d} "
              f"{len(res):7d} {share * 100:8.1f}%  "
              f"{'OK' if ok else 'MISMATCH'}")
        if not ok:
            print(f"          residues mod {expect}: {res.most_common(6)}")
    print(f"  -> {'microstep switching verified' if ok_all else 'SWITCHING SUSPECT'}")
    return ok_all


def main(argv: list[str]) -> int:
    paths = [Path(a) for a in argv]
    if not paths:
        print(__doc__)
        return 2
    ok = all(report(p) for p in paths if p.is_file())
    print()
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
