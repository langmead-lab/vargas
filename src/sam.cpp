/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date Nov 4, 2016
 *
 * @brief
 * Provides a C++ wrapper for htslib handling SAM/BAM files.
 * @details
 * Both file types are handled transparently by htslib. Records
 * are loaded into Reads.
 *
 * Implementation.
 *
 * @file
 */


#include "sam.h"
#include "doctest.h"
#include <assert.h>
#include <numeric>

const std::string vargas::SAM::Record::REQUIRED_POS = "POS";
const std::string vargas::SAM::Record::REQUIRED_QNAME = "QNAME";
const std::string vargas::SAM::Record::REQUIRED_RNEXT = "RNEXT";
const std::string vargas::SAM::Record::REQUIRED_RNAME = "RNAME";
const std::string vargas::SAM::Record::REQUIRED_SEQ = "SEQ";
const std::string vargas::SAM::Record::REQUIRED_CIGAR = "CIGAR";
const std::string vargas::SAM::Record::REQUIRED_FLAG = "FLAG";
const std::string vargas::SAM::Record::REQUIRED_PNEXT = "PNEXT";
const std::string vargas::SAM::Record::REQUIRED_MAPQ = "MAPQ";
const std::string vargas::SAM::Record::REQUIRED_TLEN = "TLEN";
const std::string vargas::SAM::Record::REQUIRED_QUAL = "QUAL";

void vargas::SAM::Optional::add(std::string a) {
    const std::vector<std::string> s = rg::split(a, ':');
    assert(s[0].length() == 2);
    if (s.size() == 2) {
        aux[s[0]] = s[1];
        aux_fmt[s[0]] = 'Z';
    } else if (s.size() >= 3) {
        assert(s[1].length() == 1);
        aux_fmt[s[0]] = s[1].at(0);
        aux[s[0]] = a.substr(5);
    } else throw std::invalid_argument("Invalid format: " + a);
}

std::string vargas::SAM::Optional::to_string() const {
    std::ostringstream ss;
    for (auto &pair : aux) {
        ss << '\t' << pair.first
           << ':' << aux_fmt.at(pair.first)
           << ':' << pair.second;
    }
    return ss.str();
}

std::string vargas::SAM::Header::Sequence::to_string() const {
    std::ostringstream ss;
    ss << "@SQ"
       << "\tSN:" << name
       << "\tLN:" << std::to_string(len)
       << (genome_assembly.length() > 0 ? "\tAS:" + genome_assembly : "")
       << (md5.length() > 0 ? "\tM5:" + md5 : "")
       << (species.length() > 0 ? "\tSP:" + species : "")
       << (URI.length() > 0 ? "\tUR:" + URI : "")
       << aux.to_string();
    return ss.str();
}

void vargas::SAM::Header::Sequence::parse(std::string line) {
    genome_assembly = "";
    md5 = "";
    species = "";
    URI = "";
    aux.clear();

    std::vector<std::string> tags = rg::split(line, '\t');
    for (auto &p : tags) {
        if (p.at(0) == '@') continue;
        std::string tag = p.substr(0, 2);
        std::string val = p.substr(3);
        if (tag == "SN") {
            name = val;
        } else if (tag == "LN") {
            len = std::stoi(val);
        } else if (tag == "AS") {
            genome_assembly = val;
        } else if (tag == "M5") {
            md5 = val;
        } else if (tag == "SP") {
            species = val;
        } else if (tag == "UR") {
            URI = val;
        } else {
            aux.add(p);
        }
    }
}

std::string vargas::SAM::Header::ReadGroup::to_string() const {
    std::ostringstream ss;
    ss << "@RG" <<
       "\tID:" << id <<
       (seq_center.length() > 0 ? "\tCN:" + seq_center : "") <<
       (desc.length() > 0 ? "\tDS:" + desc : "") <<
       (date.length() > 0 ? "\tDT:" + date : "") <<
       (flow_order.length() > 0 ? "\tFO:" + flow_order : "") <<
       (key_seq.length() > 0 ? "\tKS:" + key_seq : "") <<
       (library.length() > 0 ? "\tLB:" + library : "") <<
       (programs.length() > 0 ? "\tPG:" + programs : "") <<
       (insert_size.length() > 0 ? "\tPI:" + insert_size : "") <<
       (platform.length() > 0 ? "\tPL:" + platform : "") <<
       (platform_model.length() > 0 ? "\tPM:" + platform_model : "") <<
       (platform_unit.length() > 0 ? "\tPU:" + platform_unit : "") <<
       (sample.length() > 0 ? "\tSM:" + sample : "") << aux.to_string();
    return ss.str();
}

void vargas::SAM::Header::ReadGroup::parse(std::string line) {
    seq_center = "";
    desc = "";
    date = "";
    flow_order = "";
    key_seq = "";
    library = "";
    programs = "";
    insert_size = "";
    platform = "";
    platform_model = "";
    platform_unit = "";
    sample = "";
    aux.clear();

    std::vector<std::string> tags = rg::split(line, '\t');
    for (auto &p : tags) {
        if (p.at(0) == '@') continue;
        std::string tag = p.substr(0, 2);
        std::string val = p.substr(3);
        if (tag == "ID") {
            id = val;
        } else if (tag == "CN") {
            seq_center = val;
        } else if (tag == "DS") {
            desc = val;
        } else if (tag == "DT") {
            date = val;
        } else if (tag == "FO") {
            flow_order = val;
        } else if (tag == "KS") {
            key_seq = val;
        } else if (tag == "LB") {
            library = val;
        } else if (tag == "PG") {
            programs = val;
        } else if (tag == "PI") {
            insert_size = val;
        } else if (tag == "PL") {
            platform = val;
        } else if (tag == "PM") {
            platform_model = val;
        } else if (tag == "PU") {
            platform_unit = val;
        } else if (tag == "SM") {
            sample = val;
        } else {
            aux.add(p);
        }
    }
}

std::string vargas::SAM::Header::Program::to_string() const {
    std::ostringstream ss;
    ss << "@PG" <<
       "\tID:" << id <<
       (name.length() > 0 ? "\tPN:" + name : "") <<
       (command_line.length() > 0 ? "\tCL:" + command_line : "") <<
       (prev_pg.length() > 0 ? "\tPP:" + prev_pg : "") <<
       (desc.length() > 0 ? "\tDS:" + desc : "") <<
       (version.length() > 0 ? "\tVN:" + version : "") << aux.to_string();
    return ss.str();
}

void vargas::SAM::Header::Program::parse(std::string line) {
    name = "";
    command_line = "";
    prev_pg = "";
    desc = "";
    version = "";
    aux.clear();

    std::vector<std::string> tags = rg::split(line, '\t');
    for (auto &p : tags) {
        if (p.at(0) == '@') continue;
        std::string tag = p.substr(0, 2);
        std::string val = p.substr(3);
        if (tag == "ID") {
            id = val;
        } else if (tag == "PN") {
            name = val;
        } else if (tag == "CL") {
            command_line = val;
        } else if (tag == "PP") {
            prev_pg = val;
        } else if (tag == "DS") {
            desc = val;
        } else if (tag == "VN") {
            version = val;
        } else {
            aux.add(p);
        }
    }
}

std::string vargas::SAM::Header::to_string() const {
    std::ostringstream ret;
    ret << "@HD" <<
        "\tVN:" << version <<
        ((sorting_order.length() > 0) ? std::string("\tSO:") + sorting_order : "") <<
        ((grouping.length() > 0) ? std::string("\tGO:") + grouping : "") <<
        "\n";
    for (auto &seq : sequences) {
        ret << seq.second.to_string() << "\n";
    }
    for (auto &rg : read_groups) {
        ret << rg.second.to_string() << "\n";
    }
    for (auto &pg : programs) {
        ret << pg.second.to_string() << "\n";
    }
    return ret.str();
}

void vargas::SAM::Header::parse(std::string hdr) {
    sorting_order = "";
    grouping = "";
    sequences.clear();
    read_groups.clear();
    programs.clear();

    std::vector<std::string> lines = rg::split(hdr, '\n');

    // @HD line
    std::vector<std::string> tags = rg::split(lines[0], '\t');
    if (tags[0] != "@HD") throw std::invalid_argument("First line must start with \"@HD\"");
    std::vector<std::string> pair;
    for (auto &p : tags) {
        if (p.at(0) == '@') continue;
        rg::split(p, ':', pair);
        if (pair.size() != 2) continue;
        if (pair[0] == "VN") {
            version = pair[1];
        } else if (pair[0] == "SO") {
            sorting_order = pair[1];
        } else if (pair[0] == "GO") {
            grouping = pair[1];
        }
    }

    // Other lines
    for (auto &l : lines) {
        std::string tag = l.substr(0, 3);
        if (tag == "@SQ") add(Sequence(l));
        else if (tag == "@RG") add(ReadGroup(l));
        else if (tag == "@PG") add(Program(l));
    }
}

int vargas::SAM::Record::Flag::encode() const {
    return (multiple ? 0x001 : 0) +
    (aligned ? 0x002 : 0) +
    (unmapped ? 0x004 : 0) +
    (next_unmapped ? 0x008 : 0) +
    (rev_complement ? 0x010 : 0) +
    (next_rev_complement ? 0x020 : 0) +
    (first ? 0x040 : 0) +
    (last ? 0x080 : 0) +
    (secondary ? 0x100 : 0) +
    (pass_fail ? 0x200 : 0) +
    (duplicate ? 0x400 : 0) +
    (supplementary ? 0x800 : 0);
}

void vargas::SAM::Record::Flag::decode(unsigned int f) {
    multiple = ((f & 0x001) != 0u);
    aligned = ((f & 0x002) != 0u);
    unmapped = ((f & 0x004) != 0u);
    next_unmapped = ((f & 0x008) != 0u);
    rev_complement = ((f & 0x010) != 0u);
    next_rev_complement = ((f & 0x020) != 0u);
    first = ((f & 0x040) != 0u);
    last = ((f & 0x080) != 0u);
    secondary = ((f & 0x100) != 0u);
    pass_fail = ((f & 0x200) != 0u);
    duplicate = ((f & 0x400) != 0u);
    supplementary = ((f & 0x800) != 0u);
}

std::string vargas::SAM::Record::to_string() const {
    std::ostringstream ss;
    ss << (query_name.size() ? query_name : "*") << '\t'
       << flag.encode() << '\t'
       << (ref_name.size() ? ref_name : "*") << '\t'
       << pos << '\t'
       << mapq << '\t'
       << cigar.to_string() << '\t'
       << (ref_next.size() ? ref_next : "*") << '\t'
       << pos_next << '\t'
       << tlen << '\t'
       << seq << '\t'
       << (qual.size() ? qual : "*")
       << aux.to_string();
    return ss.str();
}

void vargas::SAM::Record::parse(std::string line) {
    std::vector<std::string> cols = rg::split(line, '\t');
    if (cols.size() < 11) throw std::invalid_argument("Record should have at least 11 columns");
    try {
        query_name = cols[0];
        flag = std::stoi(cols[1]);
        ref_name = cols[2];
        pos = std::stoi(cols[3]);
        mapq = std::stoi(cols[4]);
        cigar = cols[5];
        ref_next = cols[6];
        pos_next = std::stoi(cols[7]);
        tlen = std::stoi(cols[8]);
        seq = cols[9];
        qual = cols[10];
        aux.clear();
    } catch (std::exception &e) {
        throw std::invalid_argument("Error parsing SAM record:\n" + line);
    }

    // Aux fields
    for (size_t i = 11; i < cols.size(); ++i) {
        aux.add(cols[i]);
    }

}
bool vargas::SAM::Record::get_required(const std::string &tag, std::string &val) const {
    if (tag == REQUIRED_POS) {
        val = rg::to_string(pos);
    } else if (tag == REQUIRED_QNAME) {
        val = query_name;
    } else if (tag == REQUIRED_RNEXT) {
        val = ref_next;
    } else if (tag == REQUIRED_RNAME) {
        val = ref_name;
    } else if (tag == REQUIRED_SEQ) {
        val = seq;
    } else if (tag == REQUIRED_CIGAR) {
        val = cigar.to_string();
    } else if (tag == REQUIRED_FLAG) {
        val = rg::to_string(flag.encode());
    } else if (tag == REQUIRED_PNEXT) {
        val = rg::to_string(pos_next);
    } else if (tag == REQUIRED_MAPQ) {
        val = rg::to_string(mapq);
    } else if (tag == REQUIRED_TLEN) {
        val = rg::to_string(tlen);
    } else if (tag == REQUIRED_QUAL) {
        val = rg::to_string(qual);
    } else return false;
    return true;
}

void vargas::isam::open(std::istream &is) {
    std::ostringstream hdr;
    while (std::getline(is, _curr_line) && _curr_line.at(0) == '@') {
        hdr << _curr_line << '\n';
    }
    if (hdr.str().length() > 0) _hdr << hdr.str();
    _pprec << _curr_line;
}

void vargas::isam::open(std::string file_name) {
    close();
    if (file_name.length() == 0) {
        _use_stdio = true;
        open(std::cin);
    }
    else {
        _use_stdio = false;
        in.open(file_name);
        if (!in.good()) throw std::invalid_argument("Error opening file \"" + file_name + "\"");
        open(in);
    }
}

vargas::isam vargas::isam::subset(size_t n) {
    if (!good()) throw std::invalid_argument("No records available.");
    std::vector<Record> pending;
    do { pending.push_back(_pprec); } while (next());
    if (pending.size() == 0) throw std::invalid_argument("No records available.");
    std::vector<size_t> idx(pending.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::random_shuffle(idx.begin(), idx.begin());
    isam ss;
    if (n >= pending.size() || n == 0) {
        ss._buff = std::move(pending);
    } else {
        ss._buff.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            ss._buff.push_back(pending[i]);
        }
    }
    ss.next();
    ss._hdr = _hdr;
    return ss;
}

bool vargas::isam::next() {
    if (_buff.size() > 0) {
        _pprec = _buff.back();
        _buff.pop_back();
        return true;
    }

    if (!std::getline((_use_stdio ? std::cin : in), _curr_line)) return false;
    _pprec.parse(_curr_line);
    return true;
}

void vargas::osam::open(std::string file_name) {
    close();
    if (file_name.length() == 0) _use_stdio = true;
    else {
        _use_stdio = false;
        out.open(file_name);
        if (!out.good()) throw std::invalid_argument("Error opening output file \"" + file_name + "\"");
    }
    (_use_stdio ? std::cout : out) << _hdr.to_string() << std::flush;
}

vargas::Cigar vargas::Cigar::operator=(const std::string &s) {
    parse(s);
    return *this;
}

void vargas::Cigar::parse(const std::string &s) {
    _cigar.clear();
    if (s.length() < 2) return;
    size_t prev = 0, curr;
    while (prev != s.length()) {
        curr = s.find_first_of(CIGAR_OPERATORS, prev);
        _cigar.push_back(std::pair<size_t, char>(prev == curr ? 1 : std::stoi(s.substr(prev, curr)), s[curr]));
        prev = ++curr;
    }
}


TEST_SUITE("SAM Parser");

TEST_CASE ("SAM Optional") {
    vargas::SAM::Optional o;
    o.set("a", "b");
    o.set("b", 1);
    o.set("c", 1.0f);
    o.set("d", "DD");

    {
        std::string val;
        CHECK(o.get("a", val));
        CHECK(val == "b");
    }

    {
        int val;
        CHECK(o.get("b", val));
        CHECK(val == 1);
    }

    {
        float val;
        CHECK(o.get("c", val));
        CHECK(val == 1.0);
    }

    {
        std::string val;
        CHECK(o.get("d", val));
        CHECK(val == "DD");
    }
}

TEST_CASE ("SAM File") {
    {
        std::ofstream ss("tmp_s.sam");
        ss << "@HD\tVN:1.0\tSO:coordinate\n"
           << "@SQ\tSN:1\tLN:249250621\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta"
           << "\tM5:1b22b98cdeb4a9304cb5d48026a85128\n"
           << "@SQ\tSN:2S\tLN:243199373\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta"
           << "\tM5:a0d9851da00400dec1098a9255ac712e\n"
           << "@SQ\tSN:3Q\tLN:198022430\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta"
           << "\tM5:fdfd811849cc2fadebc929bb925902e5\n"
           << "@RG\tID:UM0098:1\tPL:ILLUMINA\tPU:HWUSI-EAS1707-615LHAAXX-L001\tLB:80\tDT:2010-05-05T20:00:00-0400"
           << "\tSM:SD37743\tCN:UMCORE\n"
           << "@RG\tID:UM0098:2\tPL:ILLUMINA\tPU:HWUSI-EAS1707-615LHAAXX-L002\tLB:80\tDT:2010-05-05T20:00:00-0400"
           << "\tSM:SD37743\tCN:UMCORE\n"
           << "@PG\tID:bwa\tVN:0.5.4\n"
           << "1:497:R:-272+13M17D24M\t113\t1\t497\t37\t37M\t15\t100338662\t0"
           << "\tCGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAG\t0;==-==9;>>>>>=>>>>>>>>>>>=>>>>>>>>>>"
           << "\tXT:A:U\tNM:i:0\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:37\n"
           << "19:20389:F:275+18M2D19M\t99\t1\t17644\t0\t37M\t=\t17919\t314"
           << "\tTATGACTGCTAATAATACCTACACATGTTAGAACCAT\t>>>>>>>>>>>>>>>>>>>><<>>><<>>4::>>:<9"
           << "\tRG:Z:UM0098:1\tXT:A:R\tNM:i:0\tSM:i:0\tAM:i:0\tX0:i:4\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:37\n"
           << "19:20389:F:275+18M2D19M\t147\t1\t17919\t0\t18M2D19M\t=\t17644\t-314"
           << "\tGTAGTACCAACTGTAAGTCCTTATCTTCATACTTTGT\t;44999;499<8<8<<<8<<><<<<><7<;<<<>><<"
           << "\tXT:A:R\tNM:i:2\tSM:i:0\tAM:i:0\tX0:i:4\tX1:i:0\tXM:i:0\tXO:i:1\tXG:i:2\tMD:Z:18^CA19\n"
           << "9:21597+10M2I25M:R:-209\t83\t1\t21678\t0\t8M2I27M\t=\t21469\t-244"
           << "\tCACCACATCACATATACCAAGCCTGGCTGTGTCTTCT\t<;9<<5><<<<><<<>><<><>><9>><>>>9>>><>"
           << "\tXT:A:R\tNM:i:2\tSM:i:0\tAM:i:0\tX0:i:5\tX1:i:0\tXM:i:0\tXO:i:1\tXG:i:2\tMD:Z:35\n";
    }

    SUBCASE("File Wrapper") {
        SUBCASE("SAM IO") {
            try {

                {
                    vargas::isam sf("tmp_s.sam");
                    vargas::osam os("osam.sam", sf.header());
                    do {
                        os.add_record(sf.record());
                    } while (sf.next());
                }

                vargas::isam a("tmp_s.sam");
                vargas::isam b("osam.sam");
                const auto &ah = a.header();
                const auto &bh = b.header();

                std::string v1, v2;

                REQUIRE(a.record().get(ah, "SEQ", v1));
                CHECK(v1 == "CGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAG");
                REQUIRE(b.record().get(bh, vargas::SAM::Record::REQUIRED_SEQ, v1));
                CHECK(v1 == "CGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAG");

                {
                    int v;
                    REQUIRE(a.record().get(ah, "NM", v));
                    CHECK(v == 0);
                    REQUIRE(a.record().get(ah, "SM", v));
                    CHECK(v == 37);
                    REQUIRE(b.record().get(bh, "NM", v));
                    CHECK(v == 0);
                    REQUIRE(b.record().get(bh, "SM", v));
                    CHECK(v == 37);
                }
                {
                    std::string v;
                    REQUIRE(a.record().get(ah, "XT", v));
                    CHECK(v == "U");
                    REQUIRE(a.record().get(ah, "XT", v));
                    CHECK(v == "U");
                }
                {
                    REQUIRE(a.record().get(ah, "MD", v1));
                    CHECK(v1 == "37");
                    REQUIRE(a.record().get(ah, "MD", v1));
                    CHECK(v1 == "37");
                }

                do {
                    const auto &ar = a.record();
                    const auto &br = b.record();
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_POS, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_POS, v2));
                    CHECK(v1 == v2);
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_QUAL, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_QUAL, v2));
                    CHECK(v1 == v2);
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_TLEN, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_TLEN, v2));
                    CHECK(v1 == v2);
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_MAPQ, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_MAPQ, v2));
                    CHECK(v1 == v2);
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_CIGAR, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_CIGAR, v2));
                    CHECK(v1 == v2);
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_PNEXT, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_PNEXT, v2));
                    CHECK(v1 == v2);
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_FLAG, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_FLAG, v2));
                    CHECK(v1 == v2);
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_RNAME, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_RNAME, v2));
                    CHECK(v1 == v2);
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_RNEXT, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_RNEXT, v2));
                    CHECK(v1 == v2);
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_SEQ, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_SEQ, v2));
                    CHECK(v1 == v2);
                    REQUIRE(ar.get(ah, vargas::SAM::Record::REQUIRED_QUAL, v1));
                    REQUIRE(br.get(bh, vargas::SAM::Record::REQUIRED_QUAL, v2));
                    CHECK(v1 == v2);
                } while (a.next() && b.next());
                CHECK(!b.next());
                CHECK(!a.next());
            } catch (std::exception &e) {
                std::cerr << e.what() << std::endl;
                throw;
            }
        }

        SUBCASE("Subset") {
            {
                vargas::isam orig("tmp_s.sam");
                auto ss = orig.subset(2);
                CHECK(ss.record().query_name.length()); // First record is loaded
                CHECK(ss.next());
                CHECK(ss.record().query_name.length());
                CHECK_FALSE(ss.next());
                try {
                    // ss is consumed so this should throw
                    auto sss = ss.subset(1);
                    CHECK(0);
                } catch (std::exception &e) { CHECK(1); }
            }

            {
                vargas::isam orig("tmp_s.sam");
                auto s1 = orig.subset(2);
                auto s2 = s1.subset(1);
                CHECK(s2.record().query_name.length());
                CHECK_FALSE(s2.next());
            }
        }

    }

    remove("tmp_s.sam");
    remove("osam.sam");
}

TEST_CASE ("Cigar") {
    std::string s = "MI10M1D10M";
    vargas::Cigar c = s;

    REQUIRE(c.size() == 5);
    CHECK(c[0].first == 1);
    CHECK(c[0].second == 'M');
    CHECK(c[1].first == 1);
    CHECK(c[1].second == 'I');
    CHECK(c[2].first == 10);
    CHECK(c[2].second == 'M');
    CHECK(c[3].first == 1);
    CHECK(c[3].second == 'D');
    CHECK(c[4].first == 10);
    CHECK(c[4].second == 'M');
}

TEST_SUITE_END();
