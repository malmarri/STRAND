#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <cctype>
#include <cstdlib>
#include <algorithm>

using namespace std;

static const char* STRAND_VERSION = "1.2";

// CIGAR-aware mapping: maps read index -> 1-based reference position.
// Returns -1 for soft-clips and insertions.
vector<int> parseCigar(const string& cigar, int startPos) {
    vector<int> refMap;
    refMap.reserve(150);
    int refPos = startPos;
    size_t i = 0;
    while (i < cigar.size()) {
        int len = 0;
        while (i < cigar.size() && isdigit((unsigned char)cigar[i])) {
            len = len * 10 + (cigar[i] - '0');
            i++;
        }
        if (i >= cigar.size()) break;
        char op = cigar[i++];
        switch (op) {
            case 'M': case 'X': case '=':
                for (int j = 0; j < len; j++) refMap.push_back(refPos++);
                break;
            case 'I': case 'S':
                for (int j = 0; j < len; j++) refMap.push_back(-1);
                break;
            case 'D': case 'N':
                refPos += len;
                break;
            case 'H': case 'P':
                break;
            default:
                break;
        }
    }
    return refMap;
}

// Map chromosome names to numeric IDs to avoid string hashing in inner loops
unordered_map<string, int> chr2id;
vector<string> id2chr;

int getChromId(const string& chrom) {
    auto it = chr2id.find(chrom);
    if (it != chr2id.end()) return it->second;
    int new_id = (int)id2chr.size();
    chr2id[chrom] = new_id;
    id2chr.push_back(chrom);
    return new_id;
}

// Resolve a SAM contig name to a SNP-panel chromosome ID, tolerating a
// "chr" prefix mismatch (e.g. BAM "chr1" vs panel "1"). Returns -1 if absent.
unordered_map<string, int> contigCache;

int resolveContig(const string& chrom) {
    auto cached = contigCache.find(chrom);
    if (cached != contigCache.end()) return cached->second;

    int cid = -1;
    auto it = chr2id.find(chrom);
    if (it != chr2id.end()) {
        cid = it->second;
    } else {
        string alias = (chrom.compare(0, 3, "chr") == 0) ? chrom.substr(3) : "chr" + chrom;
        auto ait = chr2id.find(alias);
        if (ait != chr2id.end()) {
            cid = ait->second;
            cerr << "Note: BAM contig '" << chrom << "' matched to SNP panel contig '" << alias << "'.\n";
        }
    }
    contigCache[chrom] = cid;
    return cid;
}

// Parse a non-negative integer argument, rejecting trailing garbage.
bool parseWindow(const char* s, int& out) {
    char* endp;
    long v = strtol(s, &endp, 10);
    if (endp == s || *endp != '\0' || v < 0 || v > 1000000) return false;
    out = (int)v;
    return true;
}

int main(int argc, char* argv[]) {
    bool unmerged_mode = false;
    bool all_mode = false;
    bool bad_option = false;
    vector<char*> args;
    for (int a = 1; a < argc; a++) {
        string arg = argv[a];
        if (arg == "--version" || arg == "-v") { cout << "STRAND " << STRAND_VERSION << "\n"; return 0; }
        if (arg == "--unmerged") unmerged_mode = true;
        else if (arg == "--all") all_mode = true;
        else if (arg.compare(0, 2, "--") == 0) { cerr << "Error: unknown option " << arg << "\n"; bad_option = true; }
        else args.push_back(argv[a]);
    }
    if (all_mode && args.size() == 3) {
        cerr << "Error: --all recalibrates whole fragments; do not also give <bases_from_start> <bases_from_end>.\n";
        bad_option = true;
    }
    if (bad_option || args.size() != (all_mode ? 1u : 3u)) {
        cerr << "STRAND " << STRAND_VERSION << " - Single-strand Transition Recalibration for Ancient Nucleic DNA\n";
        cerr << "Usage: " << argv[0] << " [--unmerged] <path_to_SNP_file> <bases_from_start> <bases_from_end>\n";
        cerr << "       " << argv[0] << " [--unmerged] --all <path_to_SNP_file>\n";
        cerr << "Example: " << argv[0] << " known_snps.txt 5 2\n";
        cerr << "  --all       recalibrate strand-appropriate SNP bases along the whole fragment\n";
        cerr << "              instead of only within the end windows\n";
        cerr << "  --unmerged  EXPERIMENTAL: also handle unmerged paired-end reads\n";
        cerr << "              (assumes R1 sequences the original molecule strand)\n";
        cerr << "  --version   print version and exit\n";
        return 1;
    }
    const char* snpPath = args[0];

    int bases_from_start = 0, bases_from_end = 0;
    if (!all_mode && (!parseWindow(args[1], bases_from_start) || !parseWindow(args[2], bases_from_end))) {
        cerr << "Error: <bases_from_start> and <bases_from_end> must be non-negative integers.\n";
        return 1;
    }

    cerr << "STRAND " << STRAND_VERSION << "\n";

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    ifstream snpFile(snpPath);
    if (!snpFile.is_open()) {
        cerr << "Error opening SNP file: " << snpPath << "\n";
        return 1;
    }

    // Separate integer sets for CT and GA transition positions indexed by chromosome ID
    // snp_ct[chr_id] and snp_ga[chr_id]
    vector<unordered_set<int>> snp_ct;
    vector<unordered_set<int>> snp_ga;

    string line;
    long n_transitions = 0;
    while (getline(snpFile, line)) {
        if (line.empty() || line[0] == '#') continue;

        // Fast delimiter searching
        size_t p1 = line.find_first_of(" \t");
        if (p1 == string::npos) continue;
        size_t p2 = line.find_first_not_of(" \t", p1);
        if (p2 == string::npos) continue;
        size_t p3 = line.find_first_of(" \t", p2);
        if (p3 == string::npos) continue;
        size_t p4 = line.find_first_not_of(" \t", p3);
        if (p4 == string::npos) continue;
        size_t p5 = line.find_first_of(" \t", p4);
        if (p5 == string::npos) continue;
        size_t p6 = line.find_first_not_of(" \t", p5);
        if (p6 == string::npos) continue;
        size_t p7 = line.find_first_of(" \t\r\n", p6);

        string chrom = line.substr(0, p1);
        int pos = atoi(line.substr(p2, p3 - p2).c_str());
        string ref = line.substr(p4, p5 - p4);
        string alt = (p7 == string::npos) ? line.substr(p6) : line.substr(p6, p7 - p6);
        for (char& c : ref) c = (char)toupper((unsigned char)c);
        for (char& c : alt) c = (char)toupper((unsigned char)c);

        bool is_ct = (ref == "C" && alt == "T") || (ref == "T" && alt == "C");
        bool is_ga = (ref == "G" && alt == "A") || (ref == "A" && alt == "G");

        if (!is_ct && !is_ga) continue;

        int cid = getChromId(chrom);
        if (cid >= (int)snp_ct.size()) {
            snp_ct.resize(cid + 1);
            snp_ga.resize(cid + 1);
        }

        if (is_ct) snp_ct[cid].insert(pos);
        if (is_ga) snp_ga[cid].insert(pos);
        n_transitions++;
    }
    snpFile.close();

    if (n_transitions == 0) {
        cerr << "Error: no C/T or G/A transitions loaded from " << snpPath
             << ". Expected tab-separated columns [Chr] [Pos] [Ref] [Alt]"
             << " (EIGENSTRAT 6-column .snp files must be converted first).\n";
        return 1;
    }

    cerr << "SNP list loaded: " << n_transitions << " transitions across " << id2chr.size() << " contigs.\n";
    if (all_mode)
        cerr << "Recalibrating strand-appropriate SNP bases along the whole fragment (--all).\n";
    else
        cerr << "Restricting recalibration to first " << bases_from_start << " and last " << bases_from_end << " mapped bases.\n";
    if (unmerged_mode)
        cerr << "EXPERIMENTAL --unmerged mode: assuming R1 = original molecule strand, R2 = its complement.\n";

    vector<string> fields;
    fields.reserve(30);

    long n_records = 0, n_on_panel = 0, n_no_qual = 0, n_bases_downgraded = 0;
    long n_paired = 0, n_paired_partial = 0, n_paired_unknown_tlen = 0;

    while (getline(cin, line)) {
        if (line.empty()) continue;
        if (line[0] == '@') {
            cout << line << '\n';
            continue;
        }

        // Fast tab splitting for SAM records
        fields.clear();
        size_t start = 0, end = 0;
        while ((end = line.find('\t', start)) != string::npos) {
            fields.push_back(line.substr(start, end - start));
            start = end + 1;
        }
        fields.push_back(line.substr(start));

        if (fields.size() < 11) {
            cout << line << '\n';
            continue;
        }

        n_records++;
        int flag = atoi(fields[1].c_str());
        if ((flag & 1) && !unmerged_mode) {
            cerr << "Error: alignment '" << fields[0] << "' is an unmerged paired-end read. Without"
                 << " --unmerged, STRAND assumes merged or single-end reads and would miss damage on R2.\n"
                 << "Merge read pairs before mapping, or use the experimental --unmerged option.\n";
            return 1;
        }
        // Skip unmapped reads directly
        if (flag & 4) {
            cout << line << '\n';
            continue;
        }

        int position = atoi(fields[3].c_str());
        // Strand of the original molecule. Merged/single-end reads (and R1) share it;
        // an unmerged R2 is read from the complementary copy, so its strand is flipped.
        bool is_paired = flag & 1;
        bool is_r2 = unmerged_mode && is_paired && (flag & 128);
        bool reverse = is_r2 ? !(flag & 16) : (flag & 16);
        if (is_paired) n_paired++;
        const string& chrom = fields[2];
        const string& cigar = fields[5];
        const string& read  = fields[9];
        string& qual        = fields[10];

        int cid = resolveContig(chrom);
        if (cid == -1) {
            // Contig not present in SNP panel, pass through
            cout << line << '\n';
            continue;
        }
        n_on_panel++;
        const auto& target_snps = reverse ? snp_ga[cid] : snp_ct[cid];

        vector<int> refPos = parseCigar(cigar, position);
        if (refPos.size() != read.size()) {
            cout << line << '\n';
            continue;
        }
        // QUAL may be '*' (absent); nothing to downgrade, and indexing it would overflow
        if (qual.size() != read.size()) {
            n_no_qual++;
            cout << line << '\n';
            continue;
        }

        int firstMapped = -1, lastMapped = -1;
        for (int i = 0; i < (int)refPos.size(); i++) {
            if (refPos[i] != -1) {
                if (firstMapped == -1) firstMapped = i;
                lastMapped = i;
            }
        }

        if (firstMapped == -1) {
            cout << line << '\n';
            continue;
        }

        // A merged read spans the whole molecule. An unmerged mate may stop short of one
        // molecule end; the fragment extent is [min(POS, PNEXT), +|TLEN|-1], so the distance
        // from this alignment's edges to the fragment's edges is added to each base's
        // distance. If TLEN is unknown, treat the alignment edges as molecule ends (conservative).
        int gapLeft = 0, gapRight = 0;
        if (unmerged_mode && is_paired) {
            long tlen = labs(atol(fields[8].c_str()));
            long pnext = atol(fields[7].c_str());
            if (tlen == 0 || (flag & 8) || pnext <= 0) {
                n_paired_unknown_tlen++;
            } else {
                long fragLo = min((long)position, pnext);
                long fragHi = fragLo + tlen - 1;
                gapLeft  = (int)max(0L, refPos[firstMapped] - fragLo);
                gapRight = (int)max(0L, fragHi - refPos[lastMapped]);
                if (gapLeft > 0 || gapRight > 0) n_paired_partial++;
            }
        }

        for (int i = 0; i < (int)read.size(); i++) {
            if (refPos[i] == -1) continue;

            int distFrom5prime = i - firstMapped + gapLeft;
            int distFrom3prime = lastMapped - i + gapRight;

            // distFrom5prime/distFrom3prime are in alignment orientation; map to molecule ends
            int molDist5 = reverse ? distFrom3prime : distFrom5prime;
            int molDist3 = reverse ? distFrom5prime : distFrom3prime;
            bool in_trim_window = all_mode || molDist5 < bases_from_start || molDist3 < bases_from_end;

            if (in_trim_window) {
                if (target_snps.find(refPos[i]) != target_snps.end()) {
                    qual[i] = '!'; // Downgrade to Phred 0
                    n_bases_downgraded++;
                }
            }
        }

        for (size_t i = 0; i < fields.size(); i++) {
            cout << fields[i] << (i == fields.size() - 1 ? "" : "\t");
        }
        cout << '\n';
    }

    cerr << "Processed " << n_records << " alignments: " << n_on_panel << " on SNP panel contigs, "
         << n_bases_downgraded << " bases downgraded.\n";
    if (unmerged_mode)
        cerr << "Unmerged: " << n_paired << " paired-end alignments, " << n_paired_partial
             << " not reaching both molecule ends, " << n_paired_unknown_tlen << " with unknown TLEN (both windows kept).\n";
    if (n_no_qual > 0)
        cerr << "Warning: " << n_no_qual << " alignments had no/mismatched QUAL and were passed through unchanged.\n";
    if (n_records > 0 && n_on_panel == 0)
        cerr << "WARNING: no alignments fell on a contig in the SNP panel. Check that chromosome names"
             << " (e.g. '1' vs 'chr1', 'X' vs '23') match between the BAM and the SNP file.\n";

    return 0;
}