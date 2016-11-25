/**
 * @file
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date Nov 4, 2-16
 *
 * @brief
 * FASTAFile provides an interface to a FASTA formatted file through htslib. Implementation.
 * @details
 * An index is created for the opened file if it does not exist.
 */

#include "fasta.h"
#include "doctest.h"

void vargas::ofasta::open(std::string file_name) {
    close();
    if (file_name.length() == 0) {
        _use_stdio = true;
    } else {
        _use_stdio = false;
        _o.open(file_name);
        if (!_o.good()) throw std::invalid_argument("Error opening file \"" + file_name + "\"");
    }
}

void vargas::ofasta::write(const std::string &name, const std::string &sequence) {
    (_use_stdio ? std::cout : _o) << '>' << name << '\n';
    size_t pos = 0;
    while (pos < sequence.length()) {
        (_use_stdio ? std::cout : _o) << sequence.substr(pos, _char_per_line) << '\n';
        pos += _char_per_line;
    }
}

int vargas::ifasta::open(const std::string &file_name) {
    // Check if a Fasta index exists. If it doesn't build it.
    if (!file_exists(file_name + ".fai")) {
        if (fai_build(file_name.c_str()) != 0) {
            return -1;
        }
    }
    _index = fai_load(file_name.c_str());
    if (!_index) return -2;
    _file_name = file_name;

    _seq_names.clear();
    for (size_t i = 0; i < num_seq(); ++i) {
        _seq_names.push_back(seq_name(i));
    }
    return 0;
}

std::vector<std::pair<std::string, std::string>> vargas::ifasta::sequences() const {
    if (!_index) throw std::invalid_argument("No file loaded.");
    std::vector<std::pair<std::string, std::string>> ret;
    for (size_t i = 0; i < num_seq(); ++i) {
        std::string name = std::string(faidx_iseq(_index, i));
        ret.push_back(std::pair<std::string, std::string>(name, seq(name)));
    }
    return ret;
}

std::string vargas::ifasta::subseq(const std::string &name, int beg, int end) const {
    int len;
    char *ss = faidx_fetch_seq(_index, name.c_str(), beg, end, &len);
    if (len < 0) {
        if (len == -2) throw std::invalid_argument("Sequence \"" + name + "\" does not exist.");
        throw std::invalid_argument("htslib general error");
    }
    std::string ret(ss, len);
    free(ss);
    return ret;
}

std::string vargas::ifasta::seq_name(const size_t i) const {
    if (i > num_seq()) throw std::range_error("Out of sequence index range.");
    return std::string(faidx_iseq(_index, i));
}

TEST_SUITE("FASTA Parser");

TEST_CASE ("FASTA Reading") {
    using std::endl;
    std::string tmpfa = "tmp_tc.fa";

    SUBCASE("Basic read") {
        {
            std::ofstream fao(tmpfa);
            fao
            << ">x" << endl
            << "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGTTCCTGGTGCTATGTGTAACTAGTAATGG" << endl
            << "TAATGGATATGTTGGGCTTTTTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAA" << endl
            << "TGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTTCAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGA" << endl
            << "ACAAGTTAGTTAATCTCTCTGAACTTCAGTTTAATTATCTCTAATATGGAGATGATACTACTGACAGCAGAGGTTTGCTG" << endl
            << "TGAAGATTAAATTAGGTGATGCTTGTAAAGCTCAGGGAATAGTGCCTGGCATAGAGGAAAGCCTCTGACAACTGGTAGTT" << endl
            << "ACTGTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCCTCCTTGAGGCTCT" << endl
            << "TCTGGCTTTTCATTGTCAACACAGTCAACGCTCAATACAAGGGACATTAGGATTGGCAGTAGCTCAGAGATCTCTCTGCT" << endl
            << ">y" << endl
            << "GGAGCCAGACAAATCTGGGTTCAAATCCTGGAGCCAGACAAATCTGGGTTCAAATCCTGGAGCCAGACAAATCTGGGTTC" << endl;
        }
        vargas::ifasta fa(tmpfa);

        CHECK(fa.num_seq() == 2);
        REQUIRE(fa.sequence_names().size() == 2);
        CHECK(fa.seq_name(0) == "x");
        CHECK(fa.seq_name(1) == "y");
        CHECK(fa.subseq("x", 0, 3) == "CAAA");
        CHECK(fa.subseq("y", 0, 2) == "GGA");
        CHECK(fa.sequence_names()[0] == "x");
        CHECK(fa.sequence_names()[1] == "y");

    }

    SUBCASE("iterator") {
        {
            std::ofstream o(tmpfa);
            o << ">a\nAAA\nAA\n>b\nCCC\nCC\n>c c\nTTT\nTT\n";
        }

        SUBCASE("Normal iterator") {
            vargas::ifasta fin(tmpfa);
            auto i = fin.begin();

            CHECK(i->first == "a");
            CHECK(i->second == "AAAAA");
            ++i;

            CHECK(i->first == "b");
            CHECK(i->second == "CCCCC");
            ++i;

            CHECK(i->first == "c");
            CHECK(i->second == "TTTTT");
            ++i;

            CHECK(i == fin.end());
            ++i;
            CHECK(i == fin.end());
        }

        SUBCASE("Resuming iterator") {
            vargas::ifasta fin(tmpfa);
            {
                auto i = fin.begin("B");
                CHECK(i == fin.end());
            }

            auto i = fin.begin("b");
            CHECK(i->first == "b");
            CHECK(i->second == "CCCCC");
            ++i;

            CHECK(i->first == "c");
            CHECK(i->second == "TTTTT");
            ++i;

            CHECK(i == fin.end());
            ++i;
            CHECK(i == fin.end());

        }
    }

    remove(tmpfa.c_str());
    remove((tmpfa + ".fai").c_str());
}
TEST_CASE ("FASTA Writing") {
    SUBCASE("open constructor") {
        {
            vargas::ofasta fa("tmp_tc_wr.fa");
            fa.char_per_line(5);
            fa.write("a", "AAAAA");
            fa.write("b", "TT");
            fa.write("c", "CCCCCCCCCCCC");
        }
        std::ifstream in("tmp_tc_wr.fa");
        std::string line;

        std::getline(in, line);
        CHECK(line == ">a");
        std::getline(in, line);
        CHECK(line == "AAAAA");
        std::getline(in, line);
        CHECK(line == ">b");
        std::getline(in, line);
        CHECK(line == "TT");
        std::getline(in, line);
        CHECK(line == ">c");
        std::getline(in, line);
        CHECK(line == "CCCCC");
        std::getline(in, line);
        CHECK(line == "CCCCC");
        std::getline(in, line);
        CHECK(line == "CC");

        remove("tmp_tc_wr.fa");
    }
}

TEST_SUITE_END();