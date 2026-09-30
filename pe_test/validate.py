#!/usr/bin/env python3
"""Validate STRAND output against simulation truth.

The expected set of downgraded positions is derived from molecule geometry
(truth.tsv), NOT from SAM flags, so it is independent of the code under test:

  expected(read) = { ref positions covered by the read's aligned bases
                     that are a transition SNP of the molecule's damage class
                     (C/T for + molecules, G/A for - molecules)
                     and lie within W5 bases of the molecule's 5' end
                     or within W3 bases of its 3' end }

Usage: validate.py <truth.tsv> <snps.txt> <output.sam> <W5> <W3> <label>
"""
import re
import sys
from collections import defaultdict


def load_truth(path):
    mols = {}
    with open(path) as f:
        next(f)
        for line in f:
            mid, start, length, strand, mu, dmg = line.rstrip("\n").split("\t")
            mols[int(mid)] = dict(start=int(start), length=int(length), rev=strand == "-",
                                  mate_unmapped=mu == "1",
                                  damaged=set() if dmg == "." else set(map(int, dmg.split(","))))
    return mols


def load_snps(path):
    ct, ga = set(), set()
    for line in open(path):
        _, pos, ref, alt = line.split()[:4]
        pair = {ref, alt}
        if pair == {"C", "T"}:
            ct.add(int(pos))
        elif pair == {"G", "A"}:
            ga.add(int(pos))
    return ct, ga


def aligned_positions(pos, cigar):
    """Yield (read_index, ref_pos) for aligned bases."""
    ri, rp = 0, pos
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            for _ in range(n):
                yield ri, rp
                ri += 1
                rp += 1
        elif op in "IS":
            ri += n
        elif op in "DN":
            rp += n


def main(truth_p, snps_p, sam_p, w5, w3, label):
    w5, w3 = int(w5), int(w3)
    mols = load_truth(truth_p)
    ct, ga = load_snps(snps_p)

    stats = defaultdict(lambda: defaultdict(int))
    down = defaultdict(lambda: defaultdict(set))   # mid -> "M"/"P" -> downgraded ref positions
    covered = defaultdict(lambda: defaultdict(set))  # mid -> "M"/"P" -> aligned ref positions
    examples = []

    for line in open(sam_p):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        name, flag, pos, cigar, qual = f[0], int(f[1]), int(f[3]), f[5], f[10]
        if flag & 4:
            continue
        mid, kind = name[1:].split("_")
        m = mols[int(mid)]
        mol_lo, mol_hi = m["start"], m["start"] + m["length"] - 1
        target = ga if m["rev"] else ct

        if kind == "M":
            cat = "merged"
        else:
            mate = "R1" if flag & 64 else "R2"
            cat = f"unmerged {mate}{' (mate unmapped)' if flag & 8 else ''}"
        cat += " -" if m["rev"] else " +"
        fallback = kind == "P" and bool(flag & 8)   # TLEN unknown: extra downgrades allowed

        expected, observed, damaged_at_snp = set(), set(), set()
        for ri, rp in aligned_positions(pos, cigar):
            covered[int(mid)][kind].add(rp)
            d5 = (mol_hi - rp) if m["rev"] else (rp - mol_lo)
            d3 = (rp - mol_lo) if m["rev"] else (mol_hi - rp)
            in_window = d5 < w5 or d3 < w3
            if rp in target and in_window:
                expected.add(rp)
                if rp in m["damaged"]:
                    damaged_at_snp.add(rp)
            if qual[ri] == "!":
                observed.add(rp)
        down[int(mid)][kind] |= observed

        s = stats[cat]
        s["reads"] += 1
        s["expected"] += len(expected)
        s["downgraded"] += len(observed)
        s["missed"] += len(expected - observed)
        s["extra"] += len(observed - expected)
        s["damage_at_snp"] += len(damaged_at_snp)
        s["damage_leaked"] += len(damaged_at_snp - observed)
        bad = expected - observed or (observed - expected and not fallback)
        if bad:
            s["reads_wrong"] += 1
            if len(examples) < 5:
                examples.append(f"{name} flag={flag} {cigar}: missed={sorted(expected - observed)} "
                                f"extra={sorted(observed - expected)}")

    # Merged vs unmerged consistency (molecules present in both forms with both mates mapped),
    # over the positions covered by both representations (mates may leave an unsequenced gap)
    consistent = inconsistent = 0
    for mid, kinds in down.items():
        if mols[mid]["mate_unmapped"] or not covered[mid]["M"] or not covered[mid]["P"]:
            continue
        shared = covered[mid]["M"] & covered[mid]["P"]
        if kinds["M"] & shared == kinds["P"] & shared:
            consistent += 1
        else:
            inconsistent += 1

    print(f"\n=== {label} ===")
    hdr = f"{'category':<34}{'reads':>6}{'expect':>8}{'downgr':>8}{'missed':>8}{'extra':>7}{'dmg@SNP':>9}{'LEAKED':>8}{'wrong':>7}"
    print(hdr)
    print("-" * len(hdr))
    tot = defaultdict(int)
    for cat in sorted(stats):
        s = stats[cat]
        for k, v in s.items():
            tot[k] += v
        print(f"{cat:<34}{s['reads']:>6}{s['expected']:>8}{s['downgraded']:>8}{s['missed']:>8}"
              f"{s['extra']:>7}{s['damage_at_snp']:>9}{s['damage_leaked']:>8}{s['reads_wrong']:>7}")
    print("-" * len(hdr))
    print(f"{'TOTAL':<34}{tot['reads']:>6}{tot['expected']:>8}{tot['downgraded']:>8}{tot['missed']:>8}"
          f"{tot['extra']:>7}{tot['damage_at_snp']:>9}{tot['damage_leaked']:>8}{tot['reads_wrong']:>7}")
    if consistent + inconsistent:
        print(f"Merged vs unmerged downgraded positions identical for {consistent}/{consistent + inconsistent} molecules")
    for e in examples:
        print("  e.g.", e)
    ok = tot["reads_wrong"] == 0 and tot["damage_leaked"] == 0 and inconsistent == 0
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(*sys.argv[1:7]) else 1)
