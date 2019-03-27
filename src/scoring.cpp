/**
 * Ravi Gaddipati
 * Jan 10, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Program CL parsing and scoring structs.
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#include "scoring.h"

std::string vargas::ScoreProfile::to_string() const {
    std::ostringstream ss;
    ss << "ma=" << (int) match
       << ";mp=" << (int) mismatch_max << ',' << (int) mismatch_min
       << ";rdg=" << (int) read_gopen << ',' << (int) read_gext
       << ";rfg=" << (int) ref_gopen << ',' << (int) ref_gext
       << ";np=" << (int) ambig
       << ";ete=" << end_to_end;
    return ss.str();
}

void vargas::Results::resize(size_t size) {
    max_pos.resize(size);
    sub_pos.resize(size);
    max_count.resize(size);
    sub_count.resize(size);
    max_score.resize(size);
    sub_score.resize(size);
    max_strand.resize(size);
    sub_strand.resize(size);
}

std::vector<std::string> vargas::tokenize_cl(std::string cl) {
    std::replace_if(cl.begin(), cl.end(), isspace, ' ');
    cl.erase(std::unique(cl.begin(), cl.end(), [](char a, char b) { return a == b && a == '-'; }), cl.end());
    return rg::split(cl, " =");
}

vargas::ScoreProfile vargas::bwt2(const std::string &cl) {
    auto sp = tokenize_cl(cl);

    if (std::find(sp.begin(), sp.end(), "-U") == sp.end()) {
        throw std::invalid_argument("Bowtie2/HISAT: Unpaired read alignment expected.");
    }

    ScoreProfile ret;

    auto f = std::find(sp.begin(), sp.end(), "-local");
    if (f != sp.end()) ret.end_to_end = false;
    else ret.end_to_end = true;

    f = std::find(sp.begin(), sp.end(), "-np");
    if (f != sp.end()) ret.ambig = std::stoi(*++f);
    else ret.ambig = 1;

    ret.match = ret.end_to_end ? 0 : 2;
    if (!ret.end_to_end) {
        // Match always 0 in end to end
        f = std::find(sp.begin(), sp.end(), "-ma");
        if (f != sp.end()) ret.match = std::stoi(*++f);
    }

    f = std::find(sp.begin(), sp.end(), "-mp");
    if (f != sp.end()) {
        auto a = rg::split(*++f, ",");
        if (sp.size() == 1) ret.mismatch_max = ret.mismatch_min = std::stoi(a[0]);
        else {
            ret.mismatch_max = std::stoi(a[0]);
            ret.mismatch_min = std::stoi(a[1]);
        }
    }
    else {
        ret.mismatch_min = 2;
        ret.mismatch_max = 6;
    }

    f = std::find(sp.begin(), sp.end(), "-rfg");
    if (f != sp.end()) {
        auto a = rg::split(*++f, ",");
        if (a.size() != 2) throw std::invalid_argument("Expected gap argument format <open,extend>");
        ret.ref_gopen = std::stoi(a[0]);
        ret.ref_gext = std::stoi(a[1]);
    } else {
        ret.ref_gopen = 5;
        ret.ref_gext = 3;
    }

    f = std::find(sp.begin(), sp.end(), "-rdg");
    if (f != sp.end()) {
        auto a = rg::split(*++f, ",");
        if (a.size() != 2) throw std::invalid_argument("Expected gap argument format <open,extend>");
        ret.read_gopen = std::stoi(a[0]);
        ret.read_gext = std::stoi(a[1]);
    } else {
        ret.read_gopen = 5;
        ret.read_gext = 3;
    }

    return ret;
}

vargas::ScoreProfile vargas::bwa_mem(const std::string &cl) {
    auto sp = tokenize_cl(cl);

    ScoreProfile ret;
    ret.end_to_end = false;
    ret.ambig = 1;

    auto f = std::find(sp.begin(), sp.end(), "-A");
    if (f != sp.end()) ret.match = std::stoi(*++f);
    else ret.match = 1;

    f = std::find(sp.begin(), sp.end(), "-B");
    if (f != sp.end()) ret.mismatch_max = ret.mismatch_min = std::stoi(*++f);
    else ret.mismatch_max = ret.mismatch_min = 4;

    f = std::find(sp.begin(), sp.end(), "-O");
    if (f != sp.end()) ret.read_gopen = std::stoi(*++f);
    else ret.read_gopen = 6;

    f = std::find(sp.begin(), sp.end(), "-E");
    if (f != sp.end()) ret.read_gext = std::stoi(*++f);
    ret.read_gext = 1;

    ret.ref_gopen = ret.read_gopen;
    ret.ref_gext = ret.read_gext;

    return ret;

}

vargas::ScoreProfile vargas::program_profile(const std::string &cl) {
    if (cl.find("bowtie2") != std::string::npos) return bwt2(cl);
    else if (cl.find("hisat2") != std::string::npos) return bwt2(cl);
    else if (cl.find("bwa mem") != std::string::npos) return bwa_mem(cl);
    else throw std::invalid_argument("Unsupported program ID: " + cl);
}
