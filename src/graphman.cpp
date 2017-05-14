/**
 * Ravi Gaddipati
 * Feb 27, 2017
 * rgaddip1@jhu.edu
 *
 * @brief
 * Tools to construct a graph and a set of subgraphs.
 * Graphs are stored in a JSON format (binarized by default)
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#include <iomanip>
#include <iterator>
#include "graphman.h"


std::shared_ptr<vargas::Graph>
vargas::GraphMan::create_base(const std::string fasta, const std::string vcf, std::vector<vargas::Region> region,
                              std::string sample_filter, size_t limvar) {

    if (_nodes == nullptr) _nodes = std::make_shared<Graph::nodemap_t>();
    else _nodes->clear();

    _graphs.clear();

    // Default regions
    if (region.size() == 0) {
        vargas::ifasta ref(fasta);
        if (!ref.good()) throw std::invalid_argument("Invalid reference: " + fasta);
        for (const auto &r : ref.sequence_names()) {
            region.emplace_back(r, 0, 0);
        }
    }

    // Meta block
    _aux["vargas-build"] = __DATE__;
    _aux["date"] = rg::current_date();
    _aux["fasta"] = fasta;

    size_t nhaplo = 0;
    if (vcf.size()) {
        _aux["vcf"] = vcf;
        vargas::VCF v(vcf);
        if (!v.good()) throw std::invalid_argument("Invalid VCF: " + vcf);
        sample_filter.erase(std::remove_if(sample_filter.begin(), sample_filter.end(), isspace), sample_filter.end());
        auto vec = rg::split(sample_filter, ',');
        if (vec.size()) v.create_ingroup(vec);
        if (v.samples().size()) _aux["samples"] = rg::vec_to_str(v.samples(), ",");
        nhaplo = v.num_haplotypes();
    }

    _graphs["base"] = std::make_shared<Graph>(_nodes);

    unsigned offset = 0;

    for (auto reg : region) {
        if (_print) std::cerr << "Building \"" << reg.seq_name << "\" (offset: " << offset << ")..." << std::endl;
        GraphFactory gf(fasta, vcf);
        gf.add_sample_filter(sample_filter);
        gf.limit_variants(limvar);
        gf.set_region(reg);
        if (_assume_contig) gf.assume_contig_chr();
        auto g = gf.build(offset);
        if (_print) std::cerr << g.statistics().to_string() << "\n";
        _resolver._contig_offsets[offset] = reg.seq_name;
        offset = g.rbegin()->end_pos() + 1;
        _graphs["base"]->assimilate(g);
    }

    _graphs["base"]->set_filter(Graph::Population(nhaplo, true));
    _graphs["base"]->set_popsize(nhaplo);

    if (nhaplo) {
        _graphs["maxaf"] = std::make_shared<Graph>(*_graphs["base"], Graph::Type::MAXAF);
        _graphs["ref"] = std::make_shared<Graph>(*_graphs["base"], Graph::Type::REF);
    }
    return _graphs["base"];
}

void vargas::GraphMan::write(const std::string &filename) {
    std::ios::sync_with_stdio(false);
    std::ofstream of(filename);
    if (!of.good()) throw std::invalid_argument("Error opening file: " + filename);

    // Meta
    of << "@vgraph\n";
    for (const auto &pair : _aux) {
        of << pair.first << '\t' << pair.second << '\n';
    }

    // Contigs
    of << "\n@contigs\n";
    for (auto &o : _resolver._contig_offsets) {
        of << o.first << '\t' << o.second << '\n';
    }

    // graphs
    // Label    [node_id_list]  [edge-list a:b,c;d:b,c;]
    of << "\n@graphs\n";
    if (_print) std::cerr << "Flushing " << _graphs.size() << " graphs...\n";
    for (auto &g : _graphs) {
        of << g.first << '\t' << rg::vec_to_str(g.second->order(), ",") << '\t';
        for (auto &p : g.second->next_map()) {
            of << p.first << ':' << rg::vec_to_str(p.second, ",") << ';';
        }
        of << '\n';
    }

    // Nodes
    of << "\n@nodes\n";
    if (_print) std::cerr << "Flushing " << _nodes->size() << " nodes...\n";
    for (auto &p : *_nodes) {
        of << p.first << '\t' << p.second.end_pos() << '\t' << p.second.freq()
           << '\t' << p.second.is_pinched() << '\t' << p.second.is_ref() << '\t' << p.second.seq().size() << '\n';
        std::for_each(p.second.seq().begin(), p.second.seq().end(), [&of](rg::Base &b){of << rg::num_to_base(b);});
        of << '\n';
    }
    std::ios::sync_with_stdio(true);
}

void vargas::GraphMan::open(const std::string &filename) {
    std::ifstream in(filename);
    if (!in.good()) throw std::invalid_argument("Error opening file: " + filename);

    std::string line;
    while(!line.size() || line[0] == '#') std::getline(in, line);
    if (line != "@vgraph") throw std::invalid_argument(filename + " is not a graph file.");

    std::vector<std::string> tokens;
    _aux.clear();
    _graphs.clear();
    _resolver._contig_offsets.clear();
    _nodes = std::make_shared<Graph::nodemap_t>();

    while (std::getline(in, line) && line[0] != '@') {
        if (!line.size()) continue;
        rg::split(line, '\t', tokens);
        if (tokens.size() == 1) tokens.push_back("");
        if (tokens.size() != 2) throw std::domain_error("Invalid meta-line: " + line);
        _aux[tokens[0]] = tokens[1];
    }

    assert(line == "@contigs");

    while(std::getline(in, line) && line[0] != '@') {
        if (!line.size()) continue;
        rg::split(line, '\t', tokens);
        if (tokens.size() != 2) throw std::domain_error("Invalid contig def: " + line);
        _resolver._contig_offsets[std::stoul(tokens[0])] = tokens[1];
    }

    assert(line == "@graphs");
    if (_print) std::cerr << "Loading graphs...\n";
    std::vector<std::string> unparsed;
    while(std::getline(in, line) && line[0] != '@') {
        if (!line.size()) continue;
        rg::split(line, '\t', tokens);
        if (tokens.size() < 2) throw std::domain_error("Invalid graph definition.");
        _graphs[tokens[0]] = std::make_shared<Graph>(_nodes);
        rg::split(tokens[1], ',', unparsed);
        std::vector<unsigned> order;
        order.reserve(unparsed.size());
        std::transform(unparsed.begin(), unparsed.end(), std::back_inserter(order), [](const std::string &c){return std::stoul(c);});
        _graphs[tokens[0]]->set_order(order);

        if (tokens.size() > 2) {
            unparsed = rg::split(tokens[2], ';');
            std::vector<std::string> edge;
            for (const auto &epair : unparsed) {
                rg::split(epair, ':', edge);
                assert(edge.size() == 2);
                const unsigned from = std::stoul(edge[0]);
                for (auto &to : rg::split(edge[1], ',')) {
                    _graphs[tokens[0]]->add_edge_unchecked(from, std::stoul(to));
                }
            }
        }
    }

    assert(line =="@nodes");
    if (_print) std::cerr << "Loading nodes...\n";
    while(std::getline(in, line)) {
        if (!line.size()) continue;
        // node meta
        rg::split(line, '\t', tokens);
        if (tokens.size() != 6) throw std::invalid_argument("Invalid node definition: " + line);
        const unsigned id = std::stoul(tokens[0]);
        _nodes->emplace(id, Graph::Node{});
        auto &n = _nodes->at(id);
        n.set_id(id);
        n.set_endpos(std::stoul(tokens[1]));
        n.set_af(std::stof(tokens[2]));
        if (tokens[3] == "1") n.pinch();
        if (tokens[4] == "1") n.set_as_ref();
        const size_t seqsize = std::stoul(tokens[5]);
        std::vector<rg::Base> &seq = n.seq();
        seq.reserve(seqsize);
        // Load sequence char by char
        char c;
        while(in.get(c)) {
            if (c == '\n') break;
            seq.push_back(rg::base_to_num(c));
        }
    }
}

std::string vargas::GraphMan::derive(std::string def) {
    std::transform(def.begin(), def.end(), def.begin(), tolower);

    std::string ancestor, label, assignment;
    {
        auto pair = rg::split(def, '=');
        if (pair.size() != 2) throw std::invalid_argument("Malformed graph definition: " + def);
        label = pair[0];
        assignment = pair[1];
        auto d = pair[0].find_last_of(":");
        if (d == std::string::npos) {
            ancestor = "base";
        } else {
            ancestor = pair[0].substr(0, d);
        }

        if (_graphs.count(ancestor) == 0) {
            throw std::logic_error("Encountered ancestor \"" + ancestor + "\" before it was defined.");
        }

        if(_graphs.count(label)) {
            throw std::domain_error("Label \"" + label + "\" is already defined.");
        }

    }

    auto &&parent_population = _graphs.at(ancestor)->filter();
    if (parent_population.size() < 2) throw std::domain_error("Cannot derive from \"" + ancestor + "\". Less than 2 samples available.");
    size_t avail = parent_population.count(), amount;
    if (assignment.back() == '%') {
        assignment.pop_back();
        amount = size_t((double(avail) / 100.) * std::stod(assignment));
    } else {
        amount = std::stoul(assignment);
        if (amount > avail) throw std::invalid_argument(def + " : requests more samples than available (" + std::to_string(avail) + ").");
    }

    std::vector<size_t> idx;
    for (size_t i = 0; i < parent_population.size(); ++i) {
        if (parent_population.at(i)) idx.push_back(i);
    }

    std::shuffle(idx.begin(), idx.end(),
                 std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));

    Graph::Population newpop(parent_population.size());
    for (size_t i = 0; i < amount; ++i) newpop.set(idx[i], true);

    _graphs[label] = std::make_shared<Graph>(*_graphs.at(ancestor), newpop);
    return label;
}

TEST_CASE("Load graph") {
    const std::string jfile = "tmp.vgraph";
    const std::string jstr = R"(
# Test file
@vgraph
aux	null

@contigs
0	chr1
13	chr2

@graphs
base	0,1,2,3,4,5	0:1;1:2,3;2:4;3:4;4:5;

@nodes
0	5	1.0	1	5	1
AAAAA
1	8	1	1	3	1
GGG
2	9	0.5	0	1	1
C
3	9	0.5	0	1	1
T
4	13	1.0	1	4	1
GCGC
5	22	1	1	9	1
ACGTACGAC
)";

    std::ofstream o(jfile);
    o << jstr;
    o.close();

    vargas::GraphMan gg;
    gg.open(jfile);
    {
        REQUIRE(gg.count("base"));


        auto &g = *gg["base"];
        auto it = g.begin();

        CHECK(it->end_pos() == 5);
        CHECK(it->seq_str() == "AAAAA");
        CHECK(it->freq() == 1);
        CHECK(it->is_pinched() == true);

        ++it;
        CHECK(it->end_pos() == 8);
        CHECK(it->seq_str() == "GGG");
        CHECK(it->freq() == 1);
        CHECK(it->is_pinched() == true);

        ++it;
        CHECK(it->end_pos() == 9);
        CHECK(it->seq_str() == "C");
        CHECK(it->is_pinched() == false);

        ++it;
        CHECK(it->end_pos() == 9);
        CHECK(it->seq_str() == "T");
        CHECK(it->is_pinched() == false);

        ++it;
        CHECK(it->end_pos() == 13);
        CHECK(it->seq_str() == "GCGC");
        CHECK(it->is_pinched() == true);

        ++it;
        CHECK(it->end_pos() == 22);
        CHECK(it->seq_str() == "ACGTACGAC");
        CHECK(it->is_pinched() == true);

        ++it;
        CHECK(it == g.end());

        auto p = gg.absolute_position(13);
        CHECK(p.first == "chr1");
        CHECK(p.second == 13);
        p = gg.absolute_position(14);
        CHECK(p.first == "chr2");
        CHECK(p.second == 1);
        p = gg.absolute_position(20);
        CHECK(p.first == "chr2");
        CHECK(p.second == 7);
    }
    remove(jfile.c_str());
    gg.write(jfile);
    gg.open(jfile);
    {
        REQUIRE(gg.count("base"));


        auto &g = *gg["base"];
        auto it = g.begin();

        CHECK(it->end_pos() == 5);
        CHECK(it->seq_str() == "AAAAA");
        CHECK(it->freq() == 1);
        CHECK(it->is_pinched() == true);

        ++it;
        CHECK(it->end_pos() == 8);
        CHECK(it->seq_str() == "GGG");
        CHECK(it->freq() == 1);
        CHECK(it->is_pinched() == true);

        ++it;
        CHECK(it->end_pos() == 9);
        CHECK(it->seq_str() == "C");
        CHECK(it->is_pinched() == false);

        ++it;
        CHECK(it->end_pos() == 9);
        CHECK(it->seq_str() == "T");
        CHECK(it->is_pinched() == false);

        ++it;
        CHECK(it->end_pos() == 13);
        CHECK(it->seq_str() == "GCGC");
        CHECK(it->is_pinched() == true);

        ++it;
        CHECK(it->end_pos() == 22);
        CHECK(it->seq_str() == "ACGTACGAC");
        CHECK(it->is_pinched() == true);

        ++it;
        CHECK(it == g.end());

        auto p = gg.absolute_position(13);
        CHECK(p.first == "chr1");
        CHECK(p.second == 13);
        p = gg.absolute_position(20);
        CHECK(p.first == "chr2");
        CHECK(p.second == 7);
        p = gg.absolute_position(14);
        CHECK(p.first == "chr2");
        CHECK(p.second == 1);
    }
    remove(jfile.c_str());
}

TEST_CASE("Write graph") {
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
        << "x\t10\t.\tC\t<CN7>,<CN0>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t5\t.\tC\tT,G\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t34\t.\tC\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t39\t.\tC\tT,G\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
    }

    SUBCASE("Derive") {
        vargas::GraphMan gg;
        gg.create_base(tmpfa, tmpvcf);
        auto label = gg.derive("a=1");
        auto &&g = *gg.at(label);

        CHECK(label == "a");
        CHECK(g.pop_size() == 4);

        auto base_s = gg.at("base")->statistics();
        auto dev_s = gg.at(label)->statistics();
        CHECK(dev_s.total_length < base_s.total_length);

    }

    SUBCASE("All regions") {
        vargas::GraphMan gg;
        const std::vector<vargas::Region> reg = {vargas::Region("x", 0, 15), vargas::Region("y", 0, 15)};
        auto base = gg.create_base(tmpfa, tmpvcf, reg);
        auto giter = base->begin();

        CHECK(giter->seq_str() == "CAAATAAG");
        CHECK(giter->is_ref());

        ++giter;
        CHECK(giter->seq_str() == "G");

        ++giter;
        CHECK(giter->seq_str() == "A");

        ++giter;
        CHECK(giter->seq_str() == "C");

        ++giter;
        CHECK(giter->seq_str() == "T");
        CHECK(!giter->is_ref());

        ++giter;
        CHECK(giter->seq_str() == "C");
        CHECK(giter->is_ref());

        ++giter;
        CHECK(giter->seq_str() == "CCCCCCC");
        CHECK(!giter->is_ref());

        ++giter;
        CHECK(giter->seq_str() == "");

        ++giter;
        CHECK(giter->seq_str() == "TTGGA");

        ++giter;
        CHECK(giter->seq_str() == "GGAG");
        CHECK(giter->begin_pos() == 15);

        ++giter;
        CHECK(giter->seq_str() == "C");
        auto p = gg.absolute_position(giter->end_pos() + 1);
        CHECK(p.first == "y");
        CHECK(p.second == 5);

        ++giter;
        CHECK(giter->seq_str() == "T");
        p = gg.absolute_position(giter->end_pos() + 1);
        CHECK(p.first == "y");
        CHECK(p.second == 5);

        ++giter;
        CHECK(giter->seq_str() == "G");
        p = gg.absolute_position(giter->end_pos() + 1);
        CHECK(p.first == "y");
        CHECK(p.second == 5);

        ++giter;
        CHECK(giter->seq_str() == "CAGACAAATC");

        ++giter;
        CHECK(giter == base->end());

        p = gg.absolute_position(16);
        CHECK(p.first == "y");
        CHECK(p.second == 1);

        p = gg.absolute_position(1);
        CHECK(p.first == "x");
        CHECK(p.second == 1);

    }

    remove(tmpfa.c_str());
    remove(tmpvcf.c_str());
}
