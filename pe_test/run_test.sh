#!/usr/bin/env bash
# Tests for STRAND's experimental --unmerged mode.
# Simulates single-stranded aDNA molecules, emits each as a merged read AND an unmerged
# pair in one BAM, runs STRAND, and validates against the simulation truth.
set -euo pipefail
cd "$(dirname "$0")"
W5=5; W3=3
OUT=sim
status=0

g++ -std=c++11 -O3 ../strand.cpp -o strand

python3 simulate.py "$OUT"
samtools sort -o "$OUT/test.bam" "$OUT/test.sam" 2>/dev/null
samtools index "$OUT/test.bam"
samtools view -h -e '!(flag & 1)' -o "$OUT/merged_only.sam" "$OUT/test.bam"
echo "Test BAM: $(samtools view -c "$OUT/test.bam") records," \
     "$(samtools view -c -f 1 "$OUT/test.bam") unmerged, $(samtools view -c -F 1 "$OUT/test.bam") merged"

echo -e "\n[1] Without --unmerged, a BAM containing paired reads must be rejected"
if samtools view -h "$OUT/test.bam" | ./strand "$OUT/snps.txt" $W5 $W3 > /dev/null 2> "$OUT/reject.log"; then
    echo "FAIL: STRAND accepted unmerged reads"; status=1
else
    grep "^Error" "$OUT/reject.log"; echo "PASS"
fi

echo -e "\n[2] --unmerged on mixed merged + unmerged BAM"
samtools view -h "$OUT/test.bam" | ./strand --unmerged "$OUT/snps.txt" $W5 $W3 2> "$OUT/unmerged.log" > "$OUT/unmerged.sam"
grep -vE "^Note" "$OUT/unmerged.log"
python3 validate.py "$OUT/truth.tsv" "$OUT/snps.txt" "$OUT/unmerged.sam" $W5 $W3 "--unmerged, mixed BAM" || status=1

echo -e "\n[3] Merged reads are unaffected by --unmerged"
./strand "$OUT/snps.txt" $W5 $W3 < "$OUT/merged_only.sam" 2>/dev/null | grep -v '^@' | sort > "$OUT/a.txt"
grep -v '^@' "$OUT/unmerged.sam" | grep "_M	" | sort > "$OUT/b.txt"
if cmp -s "$OUT/a.txt" "$OUT/b.txt"; then echo "PASS: identical"; else echo "FAIL: merged output differs"; status=1; fi

echo -e "\n[4] Window-size sweep with --unmerged"
for w in "10 2" "0 7" "3 0" "1 1" "20 20"; do
    set -- $w
    samtools view -h "$OUT/test.bam" | ./strand --unmerged "$OUT/snps.txt" "$1" "$2" 2>/dev/null > "$OUT/w.sam"
    printf "W5=%-3s W3=%-3s " "$1" "$2"
    python3 validate.py "$OUT/truth.tsv" "$OUT/snps.txt" "$OUT/w.sam" "$1" "$2" x | grep RESULT || status=1
done

echo -e "\n[5] --all: every strand-appropriate SNP base on every read is downgraded"
BIG=1000000   # validator window larger than any molecule == whole fragment
./strand --all "$OUT/snps.txt" < "$OUT/merged_only.sam" 2>/dev/null > "$OUT/all_merged.sam"
python3 validate.py "$OUT/truth.tsv" "$OUT/snps.txt" "$OUT/all_merged.sam" $BIG $BIG "--all, merged only" | grep -E "TOTAL|RESULT" || status=1
samtools view -h "$OUT/test.bam" | ./strand --unmerged --all "$OUT/snps.txt" 2>/dev/null > "$OUT/all_mixed.sam"
python3 validate.py "$OUT/truth.tsv" "$OUT/snps.txt" "$OUT/all_mixed.sam" $BIG $BIG "--all --unmerged, mixed BAM" | grep -E "TOTAL|identical|RESULT" || status=1

echo -e "\n[6] --all matches windows larger than any read"
samtools view -h "$OUT/test.bam" | ./strand --unmerged "$OUT/snps.txt" $BIG $BIG 2>/dev/null > "$OUT/big.sam"
if cmp -s "$OUT/all_mixed.sam" "$OUT/big.sam"; then echo "PASS: identical"; else echo "FAIL: outputs differ"; status=1; fi

echo -e "\n[7] --all together with window sizes is rejected"
if ./strand --all "$OUT/snps.txt" 5 3 < /dev/null 2>/dev/null; then echo "FAIL: accepted"; status=1; else echo "PASS"; fi

echo; [ $status -eq 0 ] && echo "ALL TESTS PASSED" || echo "SOME TESTS FAILED"
exit $status
