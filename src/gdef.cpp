/**
 * @author Ravi Gaddipati
 * @date July 31, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Defines a set of subgraphs from a given reference
 * and vcf file.
 *
 * @file
 */

#include "gdef.h"
#include "utils.h"

const std::string vargas::GraphManager::GDEF_FILE_MARKER = "@gdef";
const std::string vargas::GraphManager::GDEF_REF = "ref";
const std::string vargas::GraphManager::GDEF_VAR = "var";
const std::string vargas::GraphManager::GDEF_REGION = "reg";
const std::string vargas::GraphManager::GDEF_NODELEN = "nlen";
const std::string vargas::GraphManager::GDEF_BASEGRAPH = "BASE";
const std::string vargas::GraphManager::GDEF_REFGRAPH = "REF";
const std::string vargas::GraphManager::GDEF_MAXAFGRAPH = "MAXAF";
const std::string vargas::GraphManager::GDEF_SAMPLE_FILTER = "FILTER";
const std::string vargas::GraphManager::GDEF_NEGATE_FILTER = "INVERT";
const char vargas::GraphManager::GDEF_NEGATE = '~';
const char vargas::GraphManager::GDEF_SCOPE = ':';
const char vargas::GraphManager::GDEF_ASSIGN = '=';
const char vargas::GraphManager::GDEF_DELIM = ';';

vargas::GraphManager::GraphManager(std::string gdef_file) {
    if (!open(gdef_file)) throw std::invalid_argument("Invalid GDEF file \"" + gdef_file + "\"");
}


void vargas::GraphManager::close() {
    _subgraph_filters.clear();
    _subgraphs.clear();
    _ref_file = _variant_file = _region = "";
}


bool vargas::GraphManager::open(std::string file_name, bool build_base) {
    if (file_name.length() == 0) open(std::cin);
    std::ifstream in(file_name);
    if (!in.good()) return false;
    return open(in, build_base);
}


bool vargas::GraphManager::open(std::istream &in, bool build_base) {
    close();
    std::string line;

    // Check file type, get next line
    if (!std::getline(in, line) || line != GDEF_FILE_MARKER || !std::getline(in, line)) return false;

    // Pull meta info
    {
        std::vector<std::string> meta_split = rg::split(line, GDEF_DELIM);
        std::vector<std::string> tv_pair;
        for (const std::string &tv : meta_split) {
            rg::split(tv, GDEF_ASSIGN, tv_pair);
            if (tv_pair.size() != 2) throw std::invalid_argument("Invalid token: \"" + tv + "\"");
            const std::string &tag = tv_pair[0];
            const std::string &val = tv_pair[1];
            if (tag == GDEF_REF) _ref_file = val;
            else if (tag == GDEF_VAR) _variant_file = val;
            else if (tag == GDEF_REGION) _region = val;
            else if (tag == GDEF_NODELEN) _node_len = std::stoi(val);
            else if (tag == GDEF_SAMPLE_FILTER) _sample_filter = val;
            else if (tag == GDEF_NEGATE_FILTER) _invert_filter = val == "1";
        }
    }

    // Build base graph
    size_t nsamps;
    {
        GraphFactory gb(_ref_file);
        nsamps = gb.open_vcf(_variant_file);
        if (_sample_filter != "-") {
            nsamps = gb.add_sample_filter(_sample_filter, _invert_filter);
        }


        gb.set_region(_region);
        gb.node_len(_node_len);

        if (build_base) {
            _subgraphs[GDEF_BASEGRAPH] = std::make_shared<Graph>(gb.build());
            if (!_subgraphs.at(GDEF_BASEGRAPH)->validate()) {
                throw std::domain_error("Invalid graph- invalid node ordering.");
            }
        }
    }

    // subgraphs
    {
        std::vector<std::string> p_pair;
        Graph::Population pop(nsamps);
        while (std::getline(in, line)) {
            rg::split(line, GDEF_ASSIGN, p_pair);

            if (p_pair.size() != 2) {
                throw std::invalid_argument("Invalid token: \"" + line + "\"");
            }

            if (p_pair[1] == "-") p_pair[1] = "";

            if (p_pair[1].length() != nsamps)
                throw std::range_error("Population length does not match variant file: \"" + p_pair[0] + "\","
                " expected " + std::to_string(nsamps) + " got " + std::to_string(p_pair[1].length()));

            pop.reset();
            for (size_t i = 0; i < p_pair[1].length(); ++i) {
                if (p_pair[1][i] == '1') pop.set(i);
            }
            if (_subgraph_filters.count(p_pair[0]))
                throw std::invalid_argument("Duplicate definition: \"" + p_pair[0] + "\"");

            _subgraph_filters[p_pair[0]] = pop;


        }
    }

    return true;
}


std::shared_ptr<const vargas::Graph> vargas::GraphManager::make_subgraph(std::string label) {
    if (!_subgraphs.count(GDEF_BASEGRAPH)) throw std::invalid_argument("No base graph built.");
    if (label == GDEF_BASEGRAPH) return base();
    label = GDEF_BASEGRAPH + GDEF_SCOPE + label;

    if (rg::ends_with(label, GDEF_REFGRAPH)) return make_ref(label);
    if (rg::ends_with(label, GDEF_MAXAFGRAPH)) return make_maxaf(label);

    if (_subgraphs.count(label)) return _subgraphs.at(label);

    if (!_subgraph_filters.count(label)) throw std::invalid_argument("Label \"" + label + "\" does not exist.");
    auto sub = std::make_shared<const Graph>(*(_subgraphs[GDEF_BASEGRAPH]), _subgraph_filters.at(label));
    #pragma omp critical(_gdef_make_subgraph)
    {
        _subgraphs[label] = sub;
    }
    return _subgraphs.at(label);
}


std::shared_ptr<const vargas::Graph> vargas::GraphManager::subgraph(std::string label) const {
    if (label == GDEF_BASEGRAPH) return base();
    label = GDEF_BASEGRAPH + GDEF_SCOPE + label;
    if (_subgraphs.count(label) == 0) return nullptr;
    return _subgraphs.at(label);
}


std::shared_ptr<const vargas::Graph> vargas::GraphManager::make_ref(std::string const &label) {
    if (!_subgraphs.count(GDEF_BASEGRAPH)) throw std::invalid_argument("No base graph built.");
    if (_subgraphs.count(label)) return _subgraphs.at(label);
    std::string root = label.substr(0, label.length() - GDEF_REFGRAPH.length() - 1);
    #pragma omp critical(_gdef_make_subgraph)
    {
        _subgraphs[label] = std::make_shared<const Graph>(*make_subgraph(root), Graph::Type::REF);
    }
    return _subgraphs[label];
}


std::shared_ptr<const vargas::Graph> vargas::GraphManager::make_maxaf(std::string const &label) {
    if (!_subgraphs.count(GDEF_BASEGRAPH)) throw std::invalid_argument("No base graph built.");
    if (_subgraphs.count(label)) return _subgraphs.at(label);
    std::string root = label.substr(0, label.length() - GDEF_REFGRAPH.length() - 1);
    #pragma omp critical(_gdef_make_subgraph)
    {
        _subgraphs[label] = std::make_shared<const Graph>(*make_subgraph(root), Graph::Type::MAXAF);
    }
    return _subgraphs[label];
}


std::shared_ptr<const vargas::Graph> vargas::GraphManager::base() const {
    if (!_subgraphs.count(GDEF_BASEGRAPH)) return nullptr;
    return _subgraphs.at(GDEF_BASEGRAPH);
}


vargas::Graph::Population vargas::GraphManager::filter(std::string label) const {
    if (label == GDEF_BASEGRAPH) return base()->filter();
    if (rg::ends_with(label, GDEF_REFGRAPH)) return vargas::Graph::Population(0);
    if (rg::ends_with(label, GDEF_MAXAFGRAPH)) return vargas::Graph::Population(0);
    label = GDEF_BASEGRAPH + GDEF_SCOPE + label;
    if (!_subgraph_filters.count(label)) throw std::invalid_argument("Label \"" + label + "\" does not exist.");
    return _subgraph_filters.at(label);
}


bool vargas::GraphManager::write(std::string ref_file, std::string variant_file, std::string region,
                                 const std::string &defs, int node_len, std::string out_file,
                                 bool build_base) {
    if (out_file.length() == 0)
        return write(ref_file, variant_file, region, defs, node_len, std::cout, build_base);

    std::ofstream out(out_file);
    if (!out.good()) throw std::invalid_argument("Invalid output file: \"" + out_file + "\".");
    return write(ref_file, variant_file, region, defs, node_len, out, build_base);
}


bool vargas::GraphManager::write(std::string ref_file, std::string variant_file, std::string region,
                                 std::string defs_str, int node_len, std::ostream &out,
                                 bool build_base, int nsamps) {

    if (variant_file.length() == 0) variant_file = "-";
    vargas::ifasta f(ref_file);
    if (region.length() == 0) {
        region = f.seq_name(0);
        std::cerr << "No region defined, using \"" << region << "\".\n";
    }
    if (region.find(':') == std::string::npos) region += ":0-0";

    const auto &snames = f.sequence_names();
    const auto rname = rg::split(region, ':')[0];
    if (std::find(snames.begin(), snames.end(), rname) == snames.end()) {
        std::cerr << "Available sequences:\n";
        rg::print_vec(snames);
        std::cerr << '\n';
        throw std::invalid_argument("Sequence \"" + rname + "\" does not exist.");
    }

    std::string out_str = GDEF_FILE_MARKER + "\n"
    + GDEF_REF + GDEF_ASSIGN + ref_file + GDEF_DELIM
    + GDEF_VAR + GDEF_ASSIGN + variant_file + GDEF_DELIM
    + GDEF_REGION + GDEF_ASSIGN + region + GDEF_DELIM
    + GDEF_NODELEN + GDEF_ASSIGN + std::to_string(node_len) + GDEF_DELIM
    + GDEF_NEGATE_FILTER + GDEF_ASSIGN + (_invert_filter ? "1" : "0") + GDEF_DELIM
    + GDEF_SAMPLE_FILTER + GDEF_ASSIGN + _sample_filter + '\n';

    // Replace new lines with the delim, remove any spaces
    std::replace(defs_str.begin(), defs_str.end(), '\n', GDEF_DELIM);
    defs_str.erase(std::remove_if(defs_str.begin(), defs_str.end(), isspace), defs_str.end());
    std::vector<std::string> defs = rg::split(defs_str, GDEF_DELIM);

    // Get number of samples from VCF file
    if (nsamps == 0) {
        GraphFactory gb(ref_file);
        if (gb.open_vcf(variant_file)) nsamps = gb.add_sample_filter(_sample_filter, _invert_filter);
    }

    std::unordered_map<std::string, Graph::Population> populations;

    {
        std::vector<std::string> pair;
        Graph::Population pop(nsamps);
        std::vector<int> avail_set;
        std::set<int> added;
        size_t count, r;
        std::string parent;
        size_t parent_end;

        // Base graph "BASE" uses the full filter
        {
            Graph::Population base(nsamps);
            base.set();
            populations.insert(std::pair<std::string, Graph::Population>(GDEF_BASEGRAPH, base));
            base.reset();
        }

        if (nsamps > 0) {
            for (const auto &def : defs) {
                rg::split(def, GDEF_ASSIGN, pair);
                if (pair.size() != 2) throw std::invalid_argument("Invalid assignment: \"" + def + "\".");

                pair[0] = GDEF_BASEGRAPH + GDEF_SCOPE + pair[0];
                parent_end = pair[0].find_last_of(GDEF_SCOPE);
                parent = pair[0].substr(0, parent_end);

                if (populations.count(parent) == 0)
                    throw std::invalid_argument("Parent \"" + parent + "\" not yet defined.");

                if (pair[0].at(parent_end + 1) == '~')
                    throw std::invalid_argument("Complement graphs cannot be defined explicitly: \"" + def + "\".");

                bool top_n = false;
                if (pair[1].at(pair[1].length() - 1) == '%') {
                    count = (size_t) (((double) populations.at(parent).count() / 100) *
                    std::stoi(pair[1].substr(0, pair[1].length() - 1)));
                } else if (pair[1].at(pair[1].length() - 1) == 't') {
                    top_n = true;
                    count = std::stoul(pair[1].substr(0, pair[1].length() - 1));
                } else count = std::stoul(pair[1]);

                if (count > populations.at(parent).count())
                    throw std::invalid_argument("Not enough samples available to pick " +
                    std::to_string(count) + " in definition \"" + def + "\", "
                                                + std::to_string(populations.at(parent).count()) + " available.");

                pop.reset();

                avail_set.clear();
                for (int j = 0; j < nsamps; ++j) {
                    if (populations.at(parent).at(j)) avail_set.push_back(j);
                }

                if (top_n) {
                    size_t k = 0;
                    for (int i : avail_set) {
                        pop.set(i);
                        if (++k == count) break;
                    }
                } else {
                    added.clear();
                    for (size_t k = 0; k < count;) {
                        r = rand() % avail_set.size();
                        if (added.count(avail_set[r]) == 0) {
                            ++k;
                            pop.set(avail_set[r]);
                            added.insert(avail_set[r]);
                        }
                    }
                }
                populations[parent + GDEF_SCOPE + pair[0].substr(parent_end + 1)] = pop;
                populations[parent + GDEF_SCOPE + GDEF_NEGATE + pair[0].substr(parent_end + 1)] =
                ~pop & populations.at(parent);
            }
        }
    }

    for (auto &p : populations) {
        out_str += p.first + GDEF_ASSIGN + p.second.to_string() + "\n";
    }

    out << out_str << std::flush;
    std::istringstream instream{out_str};
    instream.seekg(0);
    open(instream, build_base);

    return true;
}


std::string vargas::GraphManager::to_DOT(std::string name) const {
    std::ostringstream dot;
    dot << "digraph " << name << " {\n";
    std::string node_label;
    std::unordered_map<std::string, int> id_map;
    int ids = 0;
    for (const auto &l : _subgraph_filters) {
        node_label = l.first.substr(l.first.find_last_of(GDEF_SCOPE) + 1);
        dot << ++ids << "[ label=\"" << node_label << " : " << l.second.count() << "\" ";
        if (node_label.at(0) == GDEF_NEGATE) dot << "style=dotted ";
        dot << "];\n";
        id_map[node_label] = ids;
    }

    for (const auto &l : _subgraph_filters) {
        size_t last_scope = l.first.find_last_of(GDEF_SCOPE);
        if (last_scope == std::string::npos) continue;
        std::string root = l.first.substr(0, last_scope);
        dot << id_map.at(root.substr(root.find_last_of(GDEF_SCOPE) + 1))
            << " -> " << id_map.at(l.first.substr(last_scope + 1)) << ";\n";
    }
    dot << "labelloc=\"t\";\nlabel=\"Subgraph Name : Population Size\";\n}\n";
    return dot.str();
}

void vargas::GraphManager::set_filter(std::string filter, bool invert) {
    if (filter.length() == 0) filter = "-";
    else {
        filter.erase(std::unique(filter.begin(), filter.end(),
                                 [](const char &l, const char &r) { return std::isspace(l) && std::isspace(r); }),
                     filter.end());
        std::replace_if(filter.begin(), filter.end(), isspace, ',');
    }
    _sample_filter = filter;
    _invert_filter = invert;
}

TEST_SUITE("Graph Def");

TEST_CASE ("Graph Manager") {
    using std::endl;
    std::string tmpfa = "tmp_tc.fa";
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

    std::srand(12345);

    CHECK(vargas::GraphManager::is_definition("a=b") == true);
    CHECK(vargas::GraphManager::is_definition("a.txt") == false);

    SUBCASE("File Write Wrapper") {

        SUBCASE("VCF") {

            std::string tmpvcf = "tmp_tc.vcf";

            // Write temp VCF file
            {
                std::ofstream vcfo(tmpvcf);
                vcfo
                << "##fileformat=VCFv4.1" << endl
                << "##phasing=true" << endl
                << "##contig=<ID=x>" << endl
                << "##contig=<ID=y>" << endl
                << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
                << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Freq\">" << endl
                << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele count\">" << endl
                << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Num samples at site\">" << endl
                << "##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Num alt alleles\">" << endl
                << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"Length of each alt\">" << endl
                << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"type of variant\">" << endl
                << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2" << endl
                << "x\t9\t.\tG\tA,C,T\t99\t.\tAF=0.01,0.6,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t0|1\t2|3" << endl
                << "x\t10\t.\tC\t<CN7>,<CN0>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1"
                << endl
                << "x\t14\t.\tG\t<DUP>,<BLAH>\t99\t.\tAF=0.01,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t1|1"
                << endl
                << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1"
                << endl
                << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
            }

            std::stringstream ss;
            vargas::GraphManager gm;
            gm.write_from_vcf(tmpfa,
                              tmpvcf,
                              "x:0-10",
                              "ingroup = 2;~ingroup:1_1=1;ingroup:1_2=1;top=2t",
                              100000,
                              ss);
            gm.open(ss);

            CHECK_THROWS(gm.filter("sdf"));

            CHECK(gm.filter("ingroup").count() == 2);
            CHECK((gm.filter("ingroup") && gm.filter("~ingroup")) == 0);
            CHECK(gm.filter("ingroup:1_2").count() == 1);
            CHECK((gm.filter("ingroup:1_2") && gm.filter("ingroup:~1_2")) == 0);
            CHECK((gm.filter("ingroup") && gm.filter("~ingroup:1_1")) == 0);
            CHECK(gm.filter("~ingroup:1_1").count() == 1);

            CHECK((gm.filter("ingroup") | gm.filter("~ingroup")).count() == 4);
            CHECK((gm.filter("ingroup:1_2") | gm.filter("ingroup:~1_2")) == gm.filter("ingroup"));

            CHECK(gm.filter("top").at(0) == 1);
            CHECK(gm.filter("top").at(1) == 1);
            for (size_t i = 2; i < gm.filter("top").size(); ++i) {
                CHECK(gm.filter("top").at(i) == 0);
            }

            {
                auto in_graph = gm.make_subgraph("ingroup");
                CHECK(in_graph.use_count() == 2);
                gm.destroy("ingroup");
                CHECK(in_graph.use_count() == 1);
            }

            {
                auto in_graph = gm.make_subgraph("ingroup");
                CHECK(in_graph.use_count() == 2);
                gm.clear();
                CHECK(in_graph.use_count() == 1);
            }

            remove(tmpvcf.c_str());
        }

    }

    remove(tmpfa.c_str());
    remove((tmpfa + ".fai").c_str());

}

TEST_SUITE_END();