#!/usr/bin/env python3
"""Simulate single-stranded aDNA molecules and emit each one twice into one SAM:
as a merged read (m<id>_M) and as an unmerged read pair (m<id>_P, R1 + R2).

Damage model (single-stranded library): C->T at the first D5 and last D3 bases of
the original molecule. On a reverse-strand molecule this shows up as G->A in
reference orientation.

Outputs (in the output directory):
  ref.fa      reference (chr1)
  snps.txt    SNP panel: chr pos ref alt (1-based), dense, transitions + transversions
  test.sam    merged + unmerged records (unsorted)
  truth.tsv   id, start (1-based), length, strand, mate_unmapped, damaged positions
"""
import random
import sys
import os

SEED = 42
REF_LEN = 20000
N_MOLECULES = 3000
READ_LEN = 75            # 2x75 paired-end sequencing
D5, D3 = 5, 3            # damage placed only within these molecule-end windows
DAMAGE_P = 0.5
SNP_DENSITY = 0.35
ADAPTER = "AGATCGGAAGAGC"
P_MATE_UNMAPPED = 0.05

COMP = str.maketrans("ACGT", "TGCA")
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def revcomp(s):
    return s.translate(COMP)[::-1]


def main(outdir):
    rng = random.Random(SEED)
    os.makedirs(outdir, exist_ok=True)
    ref = "".join(rng.choice("ACGT") for _ in range(REF_LEN))

    with open(os.path.join(outdir, "ref.fa"), "w") as f:
        f.write(">chr1\n")
        for i in range(0, REF_LEN, 60):
            f.write(ref[i:i + 60] + "\n")

    with open(os.path.join(outdir, "snps.txt"), "w") as f:
        for i, b in enumerate(ref):
            if rng.random() < SNP_DENSITY:
                if rng.random() < 0.7:
                    alt = TRANSITION[b]
                else:
                    alt = rng.choice([x for x in "ACGT" if x not in (b, TRANSITION[b])])
                f.write(f"chr1\t{i + 1}\t{b}\t{alt}\n")

    sam = [f"@HD\tVN:1.6\tSO:unsorted", f"@SQ\tSN:chr1\tLN:{REF_LEN}"]
    truth = ["id\tstart\tlength\tstrand\tmate_unmapped\tdamaged_positions"]

    def rec(name, flag, pos, cigar, rnext, pnext, tlen, seq):
        sam.append("\t".join(map(str, [name, flag, "chr1", pos, 60, cigar, rnext, pnext,
                                       tlen, seq, "I" * len(seq)])))

    for mid in range(N_MOLECULES):
        # ~40% short inserts (< read length, so reads run into the adapter)
        length = rng.randint(25, READ_LEN - 1) if rng.random() < 0.4 else rng.randint(READ_LEN, 200)
        start = rng.randint(0, REF_LEN - length)          # 0-based
        rev = rng.random() < 0.5
        seg = ref[start:start + length]

        # Original molecule 5'->3', apply C->T damage near both ends
        mol = list(revcomp(seg) if rev else seg)
        damaged_idx = []
        for k in range(length):
            if (k < D5 or k >= length - D3) and mol[k] == "C" and rng.random() < DAMAGE_P:
                mol[k] = "T"
                damaged_idx.append(k)
        mol = "".join(mol)
        obs = revcomp(mol) if rev else mol               # reference orientation
        damaged_ref = sorted((start + (length - 1 - k if rev else k)) + 1 for k in damaged_idx)

        # ---- merged read: the whole molecule ----
        rec(f"m{mid}_M", 16 if rev else 0, start + 1, f"{length}M", "*", 0, 0, obs)

        # ---- unmerged pair ----
        # R1 reads the original strand from its 5' end; R2 reads the complement from the 3' end.
        l1 = l2 = min(READ_LEN, length)
        if not rev:
            r1_lo, r1_rev = start, False
            r2_lo, r2_rev = start + length - l2, True
        else:
            r1_lo, r1_rev = start + length - l1, True
            r2_lo, r2_rev = start, False

        def read_record(lo, n, is_rev):
            """SEQ/CIGAR for an aligned mate; short inserts may keep a soft-clipped adapter tail."""
            seq, cigar = obs[lo - start: lo - start + n], f"{n}M"
            if length < READ_LEN and rng.random() < 0.5:
                tail = ADAPTER[:min(READ_LEN - length, len(ADAPTER))]
                # adapter follows the sequenced 3' end: right side if forward, left if reverse
                if is_rev:
                    seq, cigar = revcomp(tail) + seq, f"{len(tail)}S{n}M"
                else:
                    seq, cigar = seq + tail, f"{n}M{len(tail)}S"
            return seq, cigar

        s1, c1 = read_record(r1_lo, l1, r1_rev)
        s2, c2 = read_record(r2_lo, l2, r2_rev)
        mate_unmapped = rng.random() < P_MATE_UNMAPPED

        if mate_unmapped:
            f1 = 1 | 8 | 64 | (16 if r1_rev else 0)
            f2 = 1 | 4 | 128 | (32 if r1_rev else 0)
            rec(f"m{mid}_P", f1, r1_lo + 1, c1, "=", r1_lo + 1, 0, s1)
            rec(f"m{mid}_P", f2, r1_lo + 1, "*", "=", r1_lo + 1, 0, s2)
        else:
            f1 = 1 | 2 | 64 | (16 if r1_rev else 0) | (32 if r2_rev else 0)
            f2 = 1 | 2 | 128 | (16 if r2_rev else 0) | (32 if r1_rev else 0)
            # TLEN: + for the leftmost mate (R1 on ties), - for the other
            t1 = length if (r1_lo < r2_lo or (r1_lo == r2_lo)) else -length
            rec(f"m{mid}_P", f1, r1_lo + 1, c1, "=", r2_lo + 1, t1, s1)
            rec(f"m{mid}_P", f2, r2_lo + 1, c2, "=", r1_lo + 1, -t1, s2)

        truth.append(f"{mid}\t{start + 1}\t{length}\t{'-' if rev else '+'}\t{int(mate_unmapped)}\t"
                     + (",".join(map(str, damaged_ref)) or "."))

    with open(os.path.join(outdir, "test.sam"), "w") as f:
        f.write("\n".join(sam) + "\n")
    with open(os.path.join(outdir, "truth.tsv"), "w") as f:
        f.write("\n".join(truth) + "\n")
    print(f"Simulated {N_MOLECULES} molecules -> {len(sam) - 2} SAM records in {outdir}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sim")
