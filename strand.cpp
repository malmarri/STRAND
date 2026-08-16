#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <sstream>
#include <iterator>
#include <cctype>
#include <stdexcept>

using namespace std;

// Parse a CIGAR string and return a vector mapping each read-string index
// to its reference position (1-based). Returns -1 for unmapped bases
// (soft-clips 'S' and insertions 'I'), which should be skipped.
// Hard-clips ('H') and padding ('P') do not appear in the read string and
// are handled without producing entries in the map.
vector<int> parseCigar(const string& cigar, int startPos) {
    vector<int> refMap;
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
                // Aligned bases: consume both read and reference
                for (int j = 0; j < len; j++) refMap.push_back(refPos++);
                break;
            case 'I':
                // Insertion into reference: consume read but not reference
                for (int j = 0; j < len; j++) refMap.push_back(-1);
                break;
            case 'S':
                // Soft-clip: present in read string but not aligned to reference
                for (int j = 0; j < len; j++) refMap.push_back(-1);
                break;
            case 'D': case 'N':
                // Deletion from reference or intron: consume reference but not read
                refPos += len;
                break;
            case 'H': case 'P':
                // Hard-clip / padding: not in read string at all; no entries added
                break;
            default:
                break;
        }
    }
    return refMap;
}

int main(int argc, char* argv[]) {
    // Check if the SNP file and the two end-distances are provided
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <path_to_SNP_file> <bases_from_start> <bases_from_end>\n";
        cerr << "Example: " << argv[0] << " known_snps.txt 5 2\n";
        cerr << "  (Recalibrates the first 5 and last 2 MAPPED bases of the original molecule)\n";
        return 1;
    }

    int bases_from_start = stoi(argv[2]);
    int bases_from_end   = stoi(argv[3]);

    // Use fast IO
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    unordered_set<string> snp_positions;

    ifstream snpFile(argv[1]);
    if (!snpFile.is_open()) {
        cerr << "Error opening SNP file.\n";
        return 1;
    }

    string line;
    while (getline(snpFile, line)) {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        vector<string> cols((istream_iterator<string>(iss)), istream_iterator<string>());
        if (cols.size() < 4) continue;

        string chrom = cols[0];
        int pos;
        try { pos = stoi(cols[1]); } catch (...) { continue; }

        // Target C/T transitions for forward strands, G/A for reverse strands
        string snp = (cols[2] == "C" && cols[3] == "T") || (cols[2] == "T" && cols[3] == "C") ? "CT" :
                     (cols[2] == "G" && cols[3] == "A") || (cols[2] == "A" && cols[3] == "G") ? "GA" : "";

        if (!snp.empty()) {
            snp_positions.insert(chrom + "_" + to_string(pos) + "_" + snp);
        }
    }
    snpFile.close();

    cerr << "SNP list loaded. Restricting recalibration to the first " << bases_from_start
         << " and last " << bases_from_end << " mapped bases of the original molecule.\n";

    vector<string> cols;
    string qual, newQual;

    while (getline(cin, line)) {
        // Pass header lines through unchanged
        if (!line.empty() && line[0] == '@') {
            cout << line << '\n';
            continue;
        }

        istringstream iss(line);
        cols.assign((istream_iterator<string>(iss)), istream_iterator<string>());

        // Require at least 11 mandatory SAM fields; pass malformed lines through
        if (cols.size() < 11) {
            cout << line << '\n';
            continue;
        }

        int flag, position;
        try {
            flag     = stoi(cols[1]);
            position = stoi(cols[3]);
        } catch (...) {
            cout << line << '\n';
            continue;
        }

        bool   reverse  = flag & 16;
        string chrom    = cols[2];
        string cigar    = cols[5];
        string read     = cols[9];
        qual            = cols[10];
        newQual         = qual;

        // Build CIGAR-aware reference position map for every read-string index.
        // refPos[i] == -1 for soft-clipped or inserted (unmapped) bases.
        vector<int> refPos = parseCigar(cigar, position);

        // Guard: refPos must cover the entire read string
        if (refPos.size() != read.size()) {
            // CIGAR/read-length mismatch (malformed record); pass through unchanged
            for (int i = 0; i < (int)cols.size(); i++)
                cout << cols[i] << (i == (int)cols.size() - 1 ? "" : "\t");
            cout << '\n';
            continue;
        }

        // Find the indices of the first and last MAPPED (aligned) base in the
        // read string. Soft-clips and insertions are excluded.
        // The trim window is anchored to the mapped region so that
        // bases_from_start always counts from the 5' end of the original molecule.
        int firstMapped = -1, lastMapped = -1;
        for (int i = 0; i < (int)refPos.size(); i++) {
            if (refPos[i] != -1) {
                if (firstMapped == -1) firstMapped = i;
                lastMapped = i;
            }
        }

        // If the read has no mapped bases at all, output as-is
        if (firstMapped == -1) {
            for (int i = 0; i < (int)cols.size(); i++)
                cout << cols[i] << (i == (int)cols.size() - 1 ? "" : "\t");
            cout << '\n';
            continue;
        }

        string snpType = reverse ? "GA" : "CT";

        for (int i = 0; i < (int)read.size(); i++) {
            // Skip unmapped bases (soft-clips, insertions): they have no reference
            // coordinate and must never be recalibrated
            if (refPos[i] == -1) continue;

            // Offsets within the mapped region from each end
            int distFrom5prime = i - firstMapped;   // distance from first mapped base
            int distFrom3prime = lastMapped - i;     // distance from last mapped base

            // For forward reads: 5' end of molecule = firstMapped (low index in BAM)
            // For reverse reads: 5' end of molecule = lastMapped  (high index in BAM,
            //   because BAM stores the reverse complement)
            // bases_from_start always refers to the 5' end of the original molecule;
            // bases_from_end always refers to the 3' end.
            bool in_trim_window;
            if (!reverse) {
                in_trim_window = (distFrom5prime < bases_from_start ||
                                  distFrom3prime < bases_from_end);
            } else {
                in_trim_window = (distFrom3prime < bases_from_start ||
                                  distFrom5prime < bases_from_end);
            }

            if (in_trim_window) {
                string key = chrom + "_" + to_string(refPos[i]) + "_" + snpType;
                if (snp_positions.count(key)) {
                    newQual[i] = '!'; // Downgrade to Phred 0
                }
            }
        }

        cols[10] = newQual;

        for (int i = 0; i < (int)cols.size(); i++) {
            cout << cols[i] << (i == (int)cols.size() - 1 ? "" : "\t");
        }
        cout << '\n';
    }

    return 0;
}