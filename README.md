# STRAND
**S**ingle-strand **T**ransition **R**ecalibration for **A**ncient **N**ucleic **D**NA

Tool to recalibrate (downgrade) base quality scores in a strand-aware manner for single-stranded ancient DNA library-derived data.

## Purpose
In single-stranded ancient DNA libraries, post-mortem deamination (C &rarr; T) is strand-specific. 
`STRAND` targets these specific transition types at a user-defined set of SNP coordinates. It identifies bases that match the expected strand-specific damage profile, and downgrades their base quality score to **0** (`!`). This keeps a substantial amount of data that is lost in comparison to trimming bases from each end of a fragment to remove damage.

## Context & Attribution
This tool allows users to define how many bases from the 5′ and 3′ ends of the original molecule should be targeted for recalibration, instead of removing them across the whole fragment which is overly conservative. Position arithmetic is fully CIGAR-aware: soft-clipped, inserted and deleted bases are handled so that downstream reference coordinates are never shifted.

*The original concept was inspired by [pontussk]([https://github.com/pontussk](https://github.com/pontussk/baSSdrop)).*

## Installation
First, clone the repository to your local machine:
```bash
git clone https://github.com/malmarri/STRAND.git
cd STRAND
```

## Compilation
Compile the tool using any standard C++ compiler:
```bash
g++ -O3 strand.cpp -o strand
```

## Usage
The tool reads SAM text from `stdin` and writes SAM text to `stdout`. It is typically used by piping `samtools view` into it:

```bash
samtools view -h input.bam | ./strand <path_to_SNP_file> <bases_from_start> <bases_from_end> | samtools view -bS - > output_recalibrated.bam
```

### Parameters
1. **`<path_to_SNP_file>`**: A tab-separated file (no header) containing known SNPs: `[Chr] [Pos] [Ref] [Alt]`.
2. **`<bases_from_start>`**: Number of bases from the **5' end of the original molecule** to recalibrate. For forward reads this corresponds to the start of the BAM sequence string; for reverse reads it corresponds to the end (since BAM stores the reverse complement). Soft-clipped and inserted bases are excluded from the window count.
3. **`<bases_from_end>`**: Number of bases from the **3' end of the original molecule** to recalibrate.

> [!NOTE]
> Damage patterns are often not symmetrical across an ancient fragment, the parameters above allow you to control the number of bases from each end to recalibrate. It is recommended to assess the damage patterns for your sample empirically using a method like mapdamage/damageprofiler and then choose an appropriate number of bases for recalibration.

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
- `read4_fwd_insertion` (`2M1I17M`, pos=100): quality `1!11111111111111111!` — inserted base skipped; recalibrates refs 101, 118
- `read5_fwd_deletion` (`2M1D17M`, pos=100): quality `1!1!!111111111111!1` — deletion advances reference correctly; recalibrates refs 101, 104, 105, 118

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.
