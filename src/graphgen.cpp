/**
 * Ravi Gaddipati
 * Feb 27, 2017
 * rgaddip1@jhu.edu
 *
 * @brief
 * Tools to construct a graph and a set of subgraphs.
 * Graphs are stored in a JSON format (binarized by default)
 *
 * @file
 */

#include "graphgen.h"

void vargas::to_json(json &j, const vargas::Graph::Node &n) {
    j = json{{"pos", n.end_pos()},
             {"af", n.freq()},
             {"z", n.is_pinched()},
             {"seq", n.seq_str()}};
}

void vargas::from_json(const json &j, vargas::Graph::Node &n) {
    n.set_endpos(j["pos"].get<unsigned>());
    n.set_af(j["af"].get<float>());
    n.set_pinch(j["z"].get<bool>());
    n.set_seq(j["seq"].get<std::string>());
}
void vargas::to_json(json &j, const vargas::Region &reg) {
    j = json{{"chr", reg.seq_name}, {"min", reg.min}, {"max", reg.max}};
}

void vargas::to_json(json &j, const vargas::Graph::Type &t) {
    switch (t) {
        case Graph::Type::MAXAF:
            j = "maxaf";
            break;
        case Graph::Type::REF:
            j = "ref";
            break;
    }
}
void vargas::from_json(const json &j, vargas::Region &reg) {
    reg.seq_name = j["chr"];
    reg.min = j["min"];
    reg.max = j["max"];
}

void vargas::from_json(const json &j, vargas::Graph::Type &t) {
    if (j == "ref") t = Graph::Type::REF;
    else if (j == "maxaf") t = Graph::Type::MAXAF;
    else throw std::domain_error("Invalid graph type.");
}

void vargas::to_json(json &j, const vargas::GraphDef &gd) {
    j = json {{"invert",gd.invert},
              {"linear",gd.linear},
              {"population", gd.population},
              {"region", gd.region}};
    if (gd.parent.size() == 0) j["parent"] = nullptr;
    else j["parent"] = gd.parent;
    if (gd.linear) {
        j["filter"] = gd.type;
    }
    else {
        j["snp"] = gd.snp_only;
        j["min-af"] = gd.min_af;
        j["max-af"] = gd.max_af;
    }
}

void vargas::from_json(const json &j, vargas::GraphDef &gd) {
    gd.invert = j["invert"];
    gd.linear = j["linear"];
    if (j["parent"].is_null()) gd.parent = "";
    else gd.parent = j["parent"];
    gd.region = j["region"].get<std::vector<Region>>();
    gd.population = j["population"].get<VCF::Population>();

    if (gd.linear) {
        gd.type = j["filter"];
    }
    else {
        gd.snp_only = j["snp"];
        gd.min_af = j["min-af"];
        gd.max_af = j["max-af"];
    }
}

std::shared_ptr<vargas::Graph>
vargas::GraphGen::create_base(const std::string fasta, const std::string vcf, std::vector<vargas::Region> region,
                              std::string sample_filter, bool print) {

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
    _j["meta"]["vargas-build"] = __DATE__;
    _j["meta"]["date"] = rg::current_date();
    _j["meta"]["fasta"] = fasta;
    VCF::Population pop;
    if (vcf.size()) {
        _j["meta"]["vcf"] = vcf;
        vargas::VCF v(vcf);
        if (!v.good()) throw std::invalid_argument("Invalid VCF: " + vcf);
        sample_filter.erase(std::remove_if(sample_filter.begin(), sample_filter.end(), isspace),
                            sample_filter.end());
        auto vec = rg::split(sample_filter, ',');
        v.create_ingroup(vec);
        pop = VCF::Population(v.samples().size(), true);
        _j["meta"]["samples"] = v.samples();
    }

    // Build base graph
    GraphDef base;
    base.parent = "";
    base.invert = false;
    base.population = pop;
    base.region = region;
    base.snp_only = false;
    base. min_af = 0;
    base.max_af = 0;
    base.linear = vcf.size() == 0;
    base.type = Graph::Type::REF;
    _j["graphs"]["base"]["def"] = base;

    _graphs["base"] = std::make_shared<Graph>(_nodes);
    unsigned offset = 0;
    for (auto reg : region) {
        if (print) std::cerr << "Building \"" << reg.seq_name << "\"...\n";
        GraphFactory gf(fasta, vcf);
        gf.add_sample_filter(sample_filter);
        gf.set_region(reg);
        auto g = gf.build(offset);
        if (print) std::cerr << g.statistics().to_string() << "\n";
        _contig_offsets[offset] = reg.seq_name;
        offset = g.rbegin()->end_pos() + 1;
        _graphs["base"]->assimilate(g);
    }
    return _graphs["base"];
}

std::shared_ptr<vargas::Graph>
vargas::GraphGen::generate_subgraph(std::string label, const vargas::GraphDef &def) {
    if (_graphs.count(def.parent) == 0) {
        throw std::domain_error("Parent graph \"" + def.parent + "\" does not exist.");
    }

    if (def.linear) {
        _graphs[label] = std::make_shared<Graph>(*_graphs[def.parent], def.type);
    }

    return _graphs[label];
}

void vargas::GraphGen::write(const std::string &filename, vargas::json_mode mode) {

    // Don't modify nodes if the nodemap is in use
    bool purge = true;
    if (_nodes.use_count() > 1) purge = false;

    // Nodes
    const std::string empty;
    for (auto &p : *_nodes) {
        _j["nodes"][std::to_string(p.first)] = p.second;
        if (purge) p.second.set_seq(empty); // Free up some memory
    }
    _nodes.reset(); // Release node map

    // Contigs
    for (auto &o : _contig_offsets) {
        _j["contigs"][o.second] = o.first;
    }

    // graphs
    for (auto &g : _graphs) {
        _j["graphs"][g.first]["nodes"] = g.second->order();
        for (auto &p : g.second->next_map()) {
            _j["graphs"][g.first]["fwd"][std::to_string(p.first)] = p.second;
        }
        g.second.reset();
    }
    _graphs.clear();

    // Write
    if (filename.size() == 0) std::cout << _j;
    else {
        std::ofstream out;
        if (mode == json_mode::plaintext) out.open(filename);
        else throw std::domain_error("Unimplemented");
        if (!out.good()) throw std::invalid_argument("Error opening file: \"" + filename + "\"");
        out << _j.dump(-1);
    }
}

void vargas::GraphGen::open(const std::string &filename, vargas::json_mode mode) {
    _j.clear();
    if (filename.size() == 0) {
        std::cin >> _j;
    }
    else {
        std::ifstream in;
        if (mode == json_mode::plaintext) in.open(filename);
        else throw std::domain_error("Unimplemented");
        if (!in.good()) throw std::invalid_argument("Error opening file: \"" + filename + "\"");
        in >> _j;
    }

    // Validate fields
    if (_j.find("contigs") == _j.end() || _j.find("graphs") == _j.end() || _j.find("nodes") == _j.end()) {
        throw std::invalid_argument("Invalid graph file.");
    }

    // Load nodes
    _nodes = std::make_shared<Graph::nodemap_t>();
    auto &nodemap = *_nodes;
    unsigned node_count = 0;
    for (json::iterator it = _j["nodes"].begin(); it != _j["nodes"].end(); ++it) {
        unsigned id = std::stoul(it.key());
        Graph::Node n = it.value();
        n.set_id(id);
        nodemap[id] = n;
        it.value() = nullptr;
        ++node_count;
    }
    _j["nodes"] = nullptr;

    // contigs
    for (json::iterator it = _j["contigs"].begin(); it != _j["contigs"].end(); ++it) {
        _contig_offsets[it.value().get<unsigned>()] = it.key();
    }
    _j["contigs"] = nullptr;

    // Load graphs
    _graphs.clear();
    for (json::iterator it = _j["graphs"].begin(); it != _j["graphs"].end(); ++it) {
        _graphs[it.key()] = std::make_shared<Graph>(_nodes);
        auto &graph = *_graphs[it.key()];
        json &jgraph = it.value();
        graph.set_order(jgraph["nodes"].get<std::vector<unsigned>>());
        jgraph["nodes"] = nullptr;
        for (json::iterator e = jgraph["fwd"].begin(); e != jgraph["fwd"].end(); ++e) {
            unsigned from = std::stoul(e.key());
            std::vector<unsigned> tovec = e.value();
            for (auto to : tovec) {
                graph.add_edge(from, to);
            }
        }
        jgraph["fwd"] = nullptr;
    }
}

std::pair<std::string, unsigned> vargas::GraphGen::absolute_position(unsigned pos) const {
    std::map<unsigned, std::string>::const_iterator lb = _contig_offsets.lower_bound(pos);
    if (lb != _contig_offsets.begin()) --lb; // For rare case that pos = 0
    return {lb->second, pos - lb->first};
}

TEST_CASE("Load graph") {
    const std::string jfile = "tmpjson.vargas";
    const std::string jstr = "{\n"
    "    \"meta\" : null,\n"
    "    \"nodes\" : {\n"
    "        \"0\" : {\"pos\": 5, \"af\" : 1.0, \"z\" : true, \"seq\" : \"AAAAA\"},\n"
    "        \"1\" : {\"pos\": 8, \"af\" : 1.0, \"z\" : true, \"seq\" : \"GGG\"},\n"
    "        \"2\" : {\"pos\": 9, \"af\" : 0.5, \"z\" : false, \"seq\" : \"C\"},\n"
    "        \"3\" : {\"pos\": 9, \"af\" : 0.5, \"z\" : false, \"seq\" : \"T\"},\n"
    "        \"4\" : {\"pos\": 13, \"af\" : 1.0, \"z\" : true, \"seq\" : \"GCGC\"},\n"
    "        \"5\" : {\"pos\" : 22, \"af\" : 1.0, \"z\" : true, \"seq\" : \"ACGTACGAC\"}\n"
    "    },\n"
    "\n"
    "    \"contigs\": {\n"
    "        \"chr1\" : 0,\n"
    "        \"chr2\" : 13\n"
    "    },\n"
    "    \n"
    "    \"graphs\": {\n"
    "        \"base\" : {\n"
    "            \"def\" : {\n"
    "                \"parent\" : null, \n"
    "                \"invert\" : false, \n"
    "                \"population\" : [1,1], \n"
    "                \"region\" : [{\"chr\": \"chr1\", \"min\": 0, \"max\": 0}, {\"chr\": \"chr2\", \"min\": 0, \"max\": 0}], \n"
    "                \"snp\" : false, \n"
    "                \"min-af\" : 0, \n"
    "                \"max-af\" : 0, \n"
    "                \"linear\" : false\n"
    "            },\n"
    "            \n"
    "            \"nodes\" : [0,1,2,3,4,5],\n"
    "            \n"
    "            \"fwd\" : {\n"
    "                \"0\" : [1], \n"
    "                \"1\" : [2,3],\n"
    "                \"2\" : [4],\n"
    "                \"3\" : [4],\n"
    "                \"4\" : [5]\n"
    "            }\n"
    "        }\n"
    "    }\n"
    "}";

    std::ofstream o(jfile);
    o << jstr;
    o.close();

    vargas::GraphGen gg;
    gg.open(jfile);
    {
        REQUIRE(gg.count("base"));

        auto def = gg.definition("base");
        CHECK(def.parent == "");
        CHECK(def.invert == false);
        REQUIRE(def.population.size() == 2);
        CHECK(def.population[0] == 1);
        CHECK(def.population[1] == 1);
        REQUIRE(def.region.size() == 2);
        auto r1 = def.region[0];
        CHECK(r1.seq_name == "chr1");
        CHECK(r1.min == 0);
        CHECK(r1.max == 0);
        r1 = def.region[1];
        CHECK(r1.seq_name == "chr2");
        CHECK(r1.min == 0);
        CHECK(r1.max == 0);
        CHECK(def.snp_only == false);
        CHECK(def.linear == false);
        CHECK(def.min_af == 0);
        CHECK(def.max_af == 0);

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

        auto def = gg.definition("base");
        CHECK(def.parent == "");
        CHECK(def.invert == false);
        REQUIRE(def.population.size() == 2);
        CHECK(def.population[0] == 1);
        CHECK(def.population[1] == 1);
        REQUIRE(def.region.size() == 2);
        auto r1 = def.region[0];
        CHECK(r1.seq_name == "chr1");
        CHECK(r1.min == 0);
        CHECK(r1.max == 0);
        r1 = def.region[1];
        CHECK(r1.seq_name == "chr2");
        CHECK(r1.min == 0);
        CHECK(r1.max == 0);
        CHECK(def.snp_only == false);
        CHECK(def.linear == false);
        CHECK(def.min_af == 0);
        CHECK(def.max_af == 0);

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

    SUBCASE("All regions") {
        vargas::GraphGen gg;
        const std::vector<vargas::Region> reg = {vargas::Region("x", 0, 15), vargas::Region("y", 0, 15)};
        auto base = gg.create_base(tmpfa, tmpvcf, reg);
        auto giter = base->begin();

        CHECK((*giter).seq_str() == "CAAATAAG");
        CHECK((*giter).is_ref());

        ++giter;
        CHECK((*giter).seq_str() == "G");

        ++giter;
        CHECK((*giter).seq_str() == "A");

        ++giter;
        CHECK((*giter).seq_str() == "C");

        ++giter;
        CHECK((*giter).seq_str() == "T");
        CHECK(!(*giter).is_ref());

        ++giter;
        CHECK((*giter).seq_str() == "C");
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