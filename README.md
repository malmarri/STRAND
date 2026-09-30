# STRAND
**S**ingle-strand **T**ransition **R**ecalibration for **A**ncient **N**ucleic **D**NA

**Version 1.2**

## Purpose
In single-stranded ancient DNA libraries, post-mortem deamination (C &rarr; T) is strand-specific. 
`STRAND` targets these specific transition types at a user-defined set of SNP coordinates. It identifies bases that match the expected strand-specific damage profile, and downgrades their base quality score to **0** (`!`). This keeps a substantial amount of data that is lost in comparison to trimming bases from each end of a fragment to remove damage.

## Context & Attribution
This tool allows users to define how many bases from the 5′ and 3′ ends of the original molecule should be targeted for recalibration, instead of removing them across the whole fragment which is overly conservative. Position arithmetic is fully CIGAR-aware: soft-clipped, inserted and deleted bases are handled so that downstream reference coordinates are never shifted.

*The original concept was inspired by bassdrop [https://github.com/pontussk/baSSdrop]

## Installation
First, clone the repository to your local machine:
```bash
git clone https://github.com/malmarri/STRAND.git
cd STRAND
```

## Compilation
Compile the tool using any standard C++ compiler:
```bash
g++ -std=c++11 -O3 strand.cpp -o strand
```

## Usage
The tool reads SAM text from `stdin` and writes SAM text to `stdout`. It is typically used by piping `samtools view` into it:

```bash
# End windows: recalibrate SNP bases near the molecule ends only
samtools view -h input.bam | ./strand [--unmerged] <path_to_SNP_file> <bases_from_start> <bases_from_end> | samtools view -bS - > output_recalibrated.bam

# Whole fragment: recalibrate strand-appropriate SNP bases along the entire fragment
samtools view -h input.bam | ./strand [--unmerged] --all <path_to_SNP_file> | samtools view -bS - > output_recalibrated.bam
```

> [!IMPORTANT]
> **Input reads must be merged (or single-end, adapter-trimmed) by default.** STRAND infers the strand of the original molecule from the read's orientation, which is only valid when each alignment represents a whole molecule. For unmerged paired-end reads, R2 is sequenced from the complementary strand, so its damage would be missed; the same limitation applies to [pileupCaller's `--singleStrandMode`](https://github.com/stschiff/sequenceTools). STRAND therefore stops with an error if it encounters a paired-end alignment (FLAG `0x1`) unless `--unmerged` is given.

Run `./strand --version` to print the installed version.

### Parameters
1. **`<path_to_SNP_file>`**: A tab-separated file (no header) containing known SNPs: `[Chr] [Pos] [Ref] [Alt]`.
2. **`<bases_from_start>`**: Number of bases from the **5' end of the original molecule** to recalibrate. For forward reads this corresponds to the start of the BAM sequence string; for reverse reads it corresponds to the end (since BAM stores the reverse complement). The window is anchored to the first/last aligned (non-soft-clipped) base; inserted bases inside the window count toward its length (distance is measured along the read, i.e. the molecule), while deleted reference bases do not. Inserted bases themselves are never recalibrated, since they have no reference coordinate.
3. **`<bases_from_end>`**: Number of bases from the **3' end of the original molecule** to recalibrate.

> [!NOTE]
> Damage patterns are often not symmetrical across an ancient fragment, the parameters above allow you to control the number of bases from each end to recalibrate. It is recommended to assess the damage patterns for your sample empirically using a method like mapdamage/damageprofiler and then choose an appropriate number of bases for recalibration.

### Whole-fragment mode: `--all`
With `--all`, the end windows are replaced by the whole fragment: every aligned base at a C/T SNP on a forward-strand molecule, and at a G/A SNP on a reverse-strand molecule, is set to quality 0. `<bases_from_start>` and `<bases_from_end>` are not given in this mode. This is the same strand-aware principle as pileupCaller's `--singleStrandMode`, applied to the BAM instead of at genotype calling.

| Data | Suggested mode |
|---|---|
| UDG-treated / partial-UDG libraries (damage restricted to the terminal bases) | End windows, which retain more data |
| Non-UDG libraries, especially low-coverage pseudohaploid calling (residual C→T damage also occurs in the fragment interior) | `--all` |

Only molecules of the damage-prone strand are affected at each SNP, so transition SNPs keep roughly half their coverage (the opposite-strand molecules) and transversion SNPs keep all of it. `--all` can be combined with `--unmerged`.

> [!IMPORTANT]
> In both modes STRAND only lowers base qualities; it does not remove bases. The downstream genotype caller must filter on base quality (e.g. `samtools mpileup -Q 20`, or pileupCaller's minimum base quality) for the recalibration to take effect.

### Experimental: `--unmerged`
> [!WARNING]
> This option is **experimental**. It has been validated on simulated data only, not yet on real sequencing data.

With `--unmerged`, STRAND also handles unmerged paired-end reads. Merged and single-end reads in the same BAM are processed exactly as without the option.
- **Strand:** R1 is assumed to sequence the original molecule; R2 sequences its complement, so its molecule strand is the opposite of its mapping orientation (a reverse-mapped R2 is checked at C/T SNPs, a forward-mapped R2 at G/A SNPs).
- **Windows:** the molecule's ends are taken from the fragment extent (`min(POS, PNEXT)` to `+|TLEN|`), so each base's distance to *both* molecule ends is exact: a mate that stops short of a molecule end gets no window there, while short or overlapping inserts are handled correctly.
- **Unknown fragment length** (mate unmapped or `TLEN = 0`): the alignment's own edges are treated as molecule ends, which is conservative (may downgrade some undamaged bases, never misses damage).

> [!CAUTION]
> The R1 = original strand assumption holds for standard single-stranded protocols (Gansauge & Meyer), but other library preparations may differ. Before using this option, confirm it on your data, e.g. with DamageProfiler run separately on R1 and R2: C→T should be elevated at the 5′ end of R1.

### Checks and diagnostics
- Chromosome names are matched with or without a `chr` prefix (e.g. BAM `chr1` ↔ SNP file `1`); other naming schemes (e.g. `23` for X) must match exactly.
- The run aborts if no C/T or G/A transitions could be loaded from the SNP file (e.g. an unconverted 6-column EIGENSTRAT `.snp`), if the window sizes are not non-negative integers, or if window sizes are given together with `--all`.
- Alignments with a missing (`*`) QUAL string are passed through unchanged.
- A summary is printed to `stderr` at the end, with a warning if no alignments fell on a contig present in the SNP file.

## Input SNP File Format
The file should be tab-separated and contain at least 4 columns: `[Chr] [Pos] [Ref] [Alt]`. All known SNPs (e.g. 1240K, 1000G, HGDP data) can be included; the program will only focus on the appropriate transitions to downgrade base quality.
```text
chr1    101    C    T
chr1    108    G    A
chr2    505    T    C
```

> [!NOTE]
> A ready-to-use, pre-formatted SNP file based on the AADR v66 2 million SNP dataset is provided directly in this repository as [v66.ADDR.2M.snp](v66.ADDR.2M.snp).
> 
> In addition, a custom Y-chromosome specific SNP list is provided as [Y_chrom_yleaf.snp](Y_chrom_yleaf.snp). This contains **909,439** unique Y-chromosome SNPs integrated from the base set and all Yleaf databases (FTDNA, ISOGG, YFull v14 and v10). It is optimized to protect diagnostic SNPs required for Y-haplogroup classification during base quality recalibration.

## Pipeline Example
To recalibrate damage-prone transitions only within the terminal **5 bp** and **5 bp** of each read:

```bash
samtools view -h input.bam | \
  ./strand known_snps.txt 5 5 | \
  samtools view -bS - > output_recalibrated.bam
```

To recalibrate strand-appropriate transitions along the whole fragment instead:

```bash
samtools view -h input.bam | \
  ./strand --all known_snps.txt | \
  samtools view -bS - > output_recalibrated.bam
```

## Testing

Two test files are provided in `test_data/`. Both use `known_snps.txt` (C→T transitions at refs 101, 104, 105, 115, 118) with trim values `5 2`.

**Basic test** (`test.sam` — simple `20M` reads):
```bash
./strand test_data/known_snps.txt 5 2 < test_data/test.sam
```
- `read1_forward`: quality `1!11!1111111111111!1` (recalibrates refs 101, 104, 118)
- `read2_reverse`: quality `1!1111111111111!11!1` (recalibrates refs 101, 115, 118)

**CIGAR test** (`test_indels.sam` — soft-clips, insertions, deletions):
```bash
./strand test_data/known_snps.txt 5 2 < test_data/test_indels.sam
```
- `read3_fwd_softclip` (`3S17M`, pos=103): quality `!!!1!!111111111111!1` — 5′ window anchored to first mapped base (ref=103); recalibrates refs 104, 105, 118
- `read4_fwd_insertion` (`2M1I17M`, pos=100): quality `1!11111111111111111!` — the inserted base occupies one window slot, so the 5′ window covers refs 100–103 and ref 104 is not recalibrated; recalibrates refs 101, 118
- `read5_fwd_deletion` (`2M1D17M`, pos=100): quality `1!1!!111111111111!1` — deletion advances reference correctly; recalibrates refs 101, 104, 105, 118

**Whole-fragment test** (`--all`):
```bash
./strand --all test_data/known_snps.txt < test_data/test.sam
```
- `read1_forward`: quality `1!11!!111111111!11!1` (all C/T SNPs: refs 101, 104, 105, 115, 118)
- `read2_reverse`: quality `1!!11111111111!!11!1` (all G/A SNPs: refs 101, 102, 114, 115, 118)

**Simulation test** (experimental `--unmerged` and `--all`; requires `samtools` and `python3`):
```bash
pe_test/run_test.sh
```
Simulates 3,000 damaged single-stranded molecules, writes each as both a merged read and an unmerged pair into one BAM, and checks STRAND's output against the simulation truth: every damaged SNP base in a molecule-end window must be downgraded, nothing else may be (except where fragment length is unknown), and merged and unmerged representations of the same molecule must give identical results. It also checks that `--all` downgrades every strand-appropriate SNP base on every read, and gives the same output as windows longer than any read.

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.
