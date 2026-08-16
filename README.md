# STRAND
**S**ingle-strand **T**ransition **R**ecalibration for **A**ncient **N**ucleic **D**NA

A high-performance CIGAR- and strand-aware tool to recalibrate base quality scores for single-stranded ancient DNA (ssDNA / aDNA) libraries at known polymorphic sites without data loss from hard trimming.

---

## Overview

In single-stranded ancient DNA libraries (prepared without UDG treatment), post-mortem cytidine deamination ($\text{C} \rightarrow \text{U} \rightarrow \text{T}$) creates elevated transition rates concentrated near the terminal ends of the sequenced fragments.

Traditional approaches either:
1. **Hard-trim read ends**: Discards all bases (including informative transversions and non-damaged bases), shortening reads and losing substantial sequence coverage and phasing power.
2. **Whole-read masking**: Overly conservative, discarding variation across the entire length of the fragment.

**STRAND** solves this by:
* Targeting only candidate damage transitions ($\text{C} \rightarrow \text{T}$ on forward strand, $\text{G} \rightarrow \text{A}$ on reverse strand) at user-defined SNP coordinates.
* Confining recalibration to user-specified terminal windows from the 5′ and 3′ ends of the original molecule.
* Downgrading targeted base quality scores to **Phred 0** (`!`). Genotype likelihood calculators and imputation engines (e.g., **GLIMPSE**, **QUILT**, **bcftools**) naturally discount these bases while retaining the rest of the read and all authentic coverage.
* Full **CIGAR awareness**: Insertion, deletion, soft-clipping, and hard-clipping operations are correctly accounted for so reference coordinates are never shifted.

---

## Installation & Compilation

### Requirements
* C++ compiler (`g++` or `clang++` with C++11 support)
* `samtools` (for SAM/BAM streaming)

### Clone & Build
```bash
git clone https://github.com/malmarri/STRAND.git
cd STRAND
g++ -O3 strand.cpp -o strand
```

---

## Usage

STRAND streams SAM records from `stdin` and writes modified SAM records to `stdout`:

```bash
samtools view -h input.bam \
  | ./strand <path_to_SNP_file> <bases_from_5prime> <bases_from_3prime> \
  | samtools view -bS - > output_recalibrated.bam
```

### Parameters
1. **`<path_to_SNP_file>`**: A headerless, tab-separated file containing known SNP coordinates:
   ```text
   [Chr]    [Pos]    [Ref]    [Alt]
   ```
2. **`<bases_from_5prime>`**: Number of mapped bases from the **5′ end of the original molecule** to target for recalibration.
3. **`<bases_from_3prime>`**: Number of mapped bases from the **3′ end of the original molecule** to target for recalibration.

---

## Example Pipeline

```bash
# Recalibrate 3 bp from the 5' end and 6 bp from the 3' end
samtools view -h -@ 4 deduplicated.bam \
  | ./strand hgdp1kgp_snps.tsv 3 6 \
  | samtools view -bS -@ 4 - > strand_recalibrated.bam

# Index output BAM
samtools index strand_recalibrated.bam
```

---

## Acknowledgments & Attribution

STRAND is developed and maintained by **Mohamed Almarri** ([@malmarri](https://github.com/malmarri)). The foundational concept was inspired by an initial prototype by [pontussk](https://github.com/pontussk), which has been significantly expanded and generalized with asymmetric windowing, complete CIGAR/indel alignment handling, and stream optimization.

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.
