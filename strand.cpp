#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <cctype>
#include <cstdlib>

using namespace std;

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

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <path_to_SNP_file> <bases_from_start> <bases_from_end>\n";
        cerr << "Example: " << argv[0] << " known_snps.txt 5 2\n";
        return 1;
    }

    int bases_from_start = atoi(argv[2]);
    int bases_from_end   = atoi(argv[3]);

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    ifstream snpFile(argv[1]);
    if (!snpFile.is_open()) {
        cerr << "Error opening SNP file: " << argv[1] << "\n";
        return 1;
    }

    // Separate integer sets for CT and GA transition positions indexed by chromosome ID
    // snp_ct[chr_id] and snp_ga[chr_id]
    vector<unordered_set<int>> snp_ct;
    vector<unordered_set<int>> snp_ga;

    string line;
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
    }
    snpFile.close();

    cerr << "SNP list loaded across " << id2chr.size() << " contigs.\n";
    cerr << "Restricting recalibration to first " << bases_from_start << " and last " << bases_from_end << " mapped bases.\n";

    vector<string> fields;
    fields.reserve(30);

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

        int flag = atoi(fields[1].c_str());
        // Skip unmapped reads directly
        if (flag & 4) {
            cout << line << '\n';
            continue;
        }

        int position = atoi(fields[3].c_str());
        bool reverse = flag & 16;
        const string& chrom = fields[2];
        const string& cigar = fields[5];
        const string& read  = fields[9];
        string& qual        = fields[10];

        auto cit = chr2id.find(chrom);
        if (cit == chr2id.end()) {
            // Contig not present in SNP panel, pass through
            cout << line << '\n';
            continue;
        }
        int cid = cit->second;
        const auto& target_snps = reverse ? snp_ga[cid] : snp_ct[cid];

        vector<int> refPos = parseCigar(cigar, position);
        if (refPos.size() != read.size()) {
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

        for (int i = 0; i < (int)read.size(); i++) {
            if (refPos[i] == -1) continue;

            int distFrom5prime = i - firstMapped;
            int distFrom3prime = lastMapped - i;

            bool in_trim_window;
            if (!reverse) {
                in_trim_window = (distFrom5prime < bases_from_start || distFrom3prime < bases_from_end);
            } else {
                in_trim_window = (distFrom3prime < bases_from_start || distFrom5prime < bases_from_end);
            }

            if (in_trim_window) {
                if (target_snps.find(refPos[i]) != target_snps.end()) {
                    qual[i] = '!'; // Downgrade to Phred 0
                }
            }
        }

        for (size_t i = 0; i < fields.size(); i++) {
            cout << fields[i] << (i == fields.size() - 1 ? "" : "\t");
        }
        cout << '\n';
    }

    return 0;
}