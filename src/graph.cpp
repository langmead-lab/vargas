/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date November 20, 2015
 *
 * @brief
 * Vargas::Graph is a DAG representation of a reference and its variants.
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include "graph.h"


unsigned vargas::Graph::Node::_newID = 0;


vargas::Graph::Graph(const std::string &ref_file, const std::string &vcf_file, const std::string &region) {
    _IDMap = std::make_shared<std::unordered_map<unsigned, Node>>();
    GraphFactory gb(ref_file);
    gb.open_vcf(vcf_file);
    gb.set_region(region);
    gb.build(*this);
}


vargas::Graph::Graph(const vargas::Graph &g, const Population &filter) {
    _IDMap = g._IDMap;
    _pop_size = g.pop_size();
    _filter = filter;

    // Add all nodes
    std::unordered_set<unsigned> includedNodes;
    for (auto &nid : g._add_order) {
        auto &n = (*_IDMap)[nid];
        if (n.belongs(filter)) {
            includedNodes.insert(nid);
            _add_order.push_back(nid);
        }
    }

    _build_derived_edges(g, includedNodes);
}


vargas::Graph::Graph(const Graph &g, Type type) {
    _IDMap = g._IDMap;
    _pop_size = g.pop_size();
    _filter = Population(_pop_size, true);

    std::unordered_set<unsigned> includedNodes;

    if (type == Type::REF) {
        for (auto &nid : g._add_order) {
            auto &n = (*_IDMap)[nid];
            if (n.is_ref()) {
                includedNodes.insert(nid);
                _add_order.push_back(nid);
            }
        }
    } else if (type == Type::MAXAF) {

        std::vector<unsigned> graphstarts;
        for (auto id : g._add_order) {
            if (g._prev_map.count(id) == 0) graphstarts.push_back(id);
        }

        for (auto start : graphstarts) {
            unsigned curr = start;
            unsigned maxid, size;
            while (true) {
                includedNodes.insert(curr);
                _add_order.push_back(curr);
                if (g._next_map.count(curr) == 0) break; // end of graph
                maxid = g._next_map.at(curr).at(0);
                size = g._next_map.at(curr).size();
                for (unsigned i = 1; i < size; ++i) {
                    const unsigned &id = g._next_map.at(curr).at(i);
                    if (g._IDMap->at(id).freq() > g._IDMap->at(maxid).freq())
                        maxid = id;
                }
                curr = maxid;
            }
        }
    }

    _build_derived_edges(g, includedNodes);

}


void vargas::Graph::_build_derived_edges(const vargas::Graph &g, const std::unordered_set<unsigned> &includedNodes) {
    // Add all edges for included nodes
    for (auto &n : includedNodes) {
        if (g._next_map.count(n) == 0) continue;
        for (auto &e : g._next_map.at(n)) {
            if (includedNodes.count(e)) {
                add_edge(n, e);
            }
        }
    }

}


unsigned vargas::Graph::add_node(const Node &n) {
    if (_IDMap->find(n.id()) != _IDMap->end()) {
        throw std::invalid_argument("Duplicate node insertion.");
    }

    _IDMap->emplace(n.id(), n);
    _add_order.push_back(n.id());
    return n.id();
}


bool vargas::Graph::add_edge(const unsigned n1, const unsigned n2) {
    // Check if the nodes exist
    if (_IDMap->count(n1) == 0 || _IDMap->count(n2) == 0) return false;

    // init if first edge to be added
    if (_next_map.count(n1) == 0) {
        _next_map[n1] = std::vector<unsigned>();
    }
    if (_prev_map.count(n2) == 0) {
        _prev_map[n2] = std::vector<unsigned>();
    }
    _next_map[n1].push_back(n2);
    _prev_map[n2].push_back(n1);
    return true;
}

void vargas::Graph::add_edge_unchecked(const unsigned n1, const unsigned n2) {
    if (_next_map.count(n1) == 0) {
        _next_map[n1] = std::vector<unsigned>();
    }
    if (_prev_map.count(n2) == 0) {
        _prev_map[n2] = std::vector<unsigned>();
    }
    _next_map[n1].push_back(n2);
    _prev_map[n2].push_back(n1);
}


std::string vargas::Graph::to_DOT(std::string name) const {
    std::ostringstream dot;
    dot << "// Each node has the sequence, followed by end_pos,allele_freq,REF\n";
    dot << "digraph \"" << name << "\" {\n";

    for (auto n = begin(); n != end(); ++n) {
        auto seq = n->seq_str();
        if (seq.size() >    9) {
            seq = seq.substr(0, 4) + ".." + seq.substr(seq.size() - 5, 4);
        }
        dot << n->id() << "[label=\"" << seq
            << "\nP:" << n->end_pos() << ", F:" << std::setprecision(3) <<   n->freq() << ", R:" << n->is_ref()
            << "\"];\n";

        for (auto x : n.outgoing()) dot << n->id() << " -> " << x << ";\n";
    }
    dot << "}\n";
    return dot.str();
}


void vargas::GraphFactory::build(vargas::Graph &g, pos_t pos_offset) {
    if (_vf == nullptr) throw std::invalid_argument("No VCF file opened.");
    g = vargas::Graph(g.node_map());
    _fa.open(_fa_file);

    auto &vf = *_vf;

    if (!_fa.good()) throw std::invalid_argument("Invalid FASTA file: " + _fa_file);

    // If no region is specified, the default is the first sequence in the FASTA file
    if (vf.region().seq_name.length() == 0) {
        vf.set_region(_fa.sequence_names()[0] + ":0-0");
    }
    else {
        const auto &snames = _fa.sequence_names();
        if (std::find(snames.begin(), snames.end(), _vf->region().seq_name) == snames.end()) {
            throw std::invalid_argument("Sequence \"" + _vf->region().seq_name + "\" not in " + _fa_file);
        }
    }

    rg::pos_t curr = vf.region().min; // The Graph has been built up to this position, exclusive
    std::unordered_set<unsigned> prev_unconnected; // ID's of nodes at the end of the Graph left unconnected
    std::unordered_set<unsigned> curr_unconnected; // ID's of nodes added that are unconnected

    size_t num_samples = vf.num_haplotypes() > 0 ? vf.num_haplotypes() : 1;
    g.set_popsize(num_samples);
    const Graph::Population all_pop(g.pop_size(), true);
    g.set_filter(all_pop);

    rg::pos_t prevpos = curr; // Used to validate that VCF is sorted

    while (vf.next()) {
        if (vf.pos() < prevpos) throw std::invalid_argument("VCF file should be sorted by position.");
        else prevpos = vf.pos();

        auto &af = vf.frequencies();
        curr = _build_linear_ref(g, prev_unconnected, curr_unconnected, curr, vf.pos(), pos_offset);
        assert(_fa.subseq(_vf->region().seq_name, curr, curr + vf.ref().length() - 1) == vf.ref());

        curr += vf.ref().length();

        // Positions of variant nodes are referenced to ref node
        {
            Graph::Node n;
            n.set_endpos(curr - 1 + pos_offset);
            n.set_seq(vf.ref());
            n.set_as_ref();
            n.set_population(vf.allele_pop(vf.ref()));
            n.set_af(af[0]);
            curr_unconnected.insert(g.add_node(n));
        }

        //alt nodes
        for (unsigned i = 1; i < vf.alleles().size(); ++i) {
            const std::string &allele = vf.alleles()[i];
            if (allele == vf.ref()) continue; // Remove duplicate nodes, REF is substituted in for unknown tags
            Graph::Population pop = vf.allele_pop(allele);
            if (g.pop_size() == 1 || (pop && all_pop)) { // Only add if someone has the allele. == 1 for KSNP
                Graph::Node n;
                n.set_endpos(curr - 1 + pos_offset);
                n.set_population(pop);
                n.set_seq(allele);
                if (af.size() > i) n.set_af(af[i]);
                n.set_not_ref();
                curr_unconnected.insert(g.add_node(n));

            }
        }
        _build_edges(g, prev_unconnected, curr_unconnected);

    }
    // Nodes after last variant
    _build_linear_ref(g, prev_unconnected, curr_unconnected, curr, vf.region().max, pos_offset);

    _fa.close();
    _vf.reset();
}


void vargas::GraphFactory::_build_edges(vargas::Graph &g, std::unordered_set<unsigned> &prev,
                                        std::unordered_set<unsigned> &curr) {
    for (unsigned pID : prev) {
        for (unsigned cID : curr) {
            g.add_edge(pID, cID);
        }
    }
    prev = curr;
    curr.clear();
}


rg::pos_t vargas::GraphFactory::_build_linear_ref(Graph &g, std::unordered_set<unsigned> &prev, std::unordered_set<unsigned> &curr,
                                            pos_t pos, pos_t target, pos_t pos_offset) {

    if (target == 0) target = _fa.seq_len(_vf->region().seq_name);
    if (pos == target) return target; // For adjacent var positions

    Graph::Node n;
    n.pinch();
    n.set_population(g.pop_size(), true);
    n.set_as_ref();
    n.set_seq(_fa.subseq(_vf->region().seq_name, pos, target - 1));
    n.set_endpos(target - 1 + pos_offset);
    curr.insert(g.add_node(n));
    _build_edges(g, prev, curr);
    return target;
}


unsigned vargas::GraphFactory::add_sample_filter(std::string filter, bool invert) {
    if (!_vf) throw std::invalid_argument("No VCF file opened, cannot add filter.");
    if (filter.length() == 0 || filter == "-") return _vf->num_haplotypes();
    std::vector<std::string> filt;
    filter.erase(std::remove_if(filter.begin(), filter.end(), isspace), filter.end());
    if (invert) {
        const auto s = _vf->samples();
        auto vcf_samples = std::unordered_set<std::string>(s.begin(), s.end());
        const auto filter_samples = rg::split(filter, ',');
        for (const auto &f : filter_samples) vcf_samples.erase(f);
        filt = std::vector<std::string>(vcf_samples.begin(), vcf_samples.end());
    } else {
        rg::split(filter, ',', filt);
    }
    _vf->create_ingroup(filt);
    return _vf->num_haplotypes();
}

unsigned vargas::GraphFactory::open_vcf(std::string const &file_name) {
    _vf.reset();
    _vf = std::unique_ptr<VCF>(new VCF(file_name));
    if (file_name.length() == 0 || file_name == "-") return 0;
    if (!_vf->good()) throw std::invalid_argument("Invalid VCF/BCF file: \"" + file_name + "\"");
    return _vf->num_haplotypes();
}


vargas::Graph::Population vargas::Graph::subset(int ingroup) const {
    vargas::Graph::Population p(_pop_size);
    for (unsigned i = 0; i < _pop_size; ++i) {
        if (rand() % 100 < ingroup) p.set(i);
    }
    return p;
}

vargas::Graph vargas::Graph::subgraph(const pos_t min, const pos_t max) const {
    Graph ret;
    std::unordered_map<unsigned, unsigned> new_to_old, old_to_new;
    unsigned new_id;
    for (const Node &n : *this) {
        if (n.end_pos() < min) continue;
        if (n.is_pinched() && n.begin_pos() > max) break;

        // begin in range
        if (n.begin_pos() >= min) {
            if (n.end_pos() <= max) {
                new_id = ret.add_node(n);
                new_to_old[new_id] = n.id();
                old_to_new[n.id()] = new_id;
            } else {
                Node cpy = n;
                std::vector<rg::Base> cropped = n.seq();
                cropped.resize(n.end_pos() - max + 1);
                cpy.set_seq(cropped);
                cpy.set_endpos(max);
                new_id = ret.add_node(cpy);
                new_to_old[new_id] = n.id();
                old_to_new[n.id()] = new_id;
            }
        }

            // Begin out of range
        else if (n.end_pos() <= max) {
            Node cpy = n;
            std::vector<rg::Base> seq = n.seq();
            std::vector<rg::Base> cropped(seq.begin() + min - n.begin_pos(), seq.end());
            cpy.set_seq(cropped);
            new_id = ret.add_node(cpy);
            new_to_old[new_id] = n.id();
            old_to_new[n.id()] = new_id;
        }
    }

    // Build Edges
    for (const auto &ids : new_to_old) {
        const auto &nid = ids.first;
        const auto &oid = ids.second;
        if (_next_map.count(oid)) {
            for (const auto &next : _next_map.at(oid)) {
                if (old_to_new.count(next)) ret.add_edge(nid, old_to_new[next]);
            }
        }
    }

    return ret;
}

bool vargas::Graph::validate() const {
    std::unordered_set<unsigned> filled;
    for (auto gi = begin(); gi != end(); ++gi) {
        filled.insert(gi->id());
        for (auto i : gi.incoming()) {
            if (filled.count(i) == 0) {
                std::cerr << "Node (ID:" << gi->id() << ", POS:" << gi->end_pos() << ")"
                          << " hit before previous node " << i << std::endl;
                return false;
            }
        }
    }
    return true;
}



TEST_SUITE("Graphs");

TEST_CASE ("Node class") {
    vargas::Graph::Node::_newID = 0;
    vargas::Graph::Node n1;
    vargas::Graph::Node n2;
    CHECK(n1.id() == 0);
    CHECK(n2.id() == 1);

    SUBCASE("Set Node params") {
        using rg::Base;
        n1.set_seq("ACGTN");
        std::vector<bool> a = {0, 0, 1};
        n1.set_population(a);
        n1.set_endpos(100);

        REQUIRE(n1.seq().size() == 5);

        CHECK(n1.seq()[0] == Base::A);
        CHECK(n1.seq()[1] == Base::C);
        CHECK(n1.seq()[2] == Base::G);
        CHECK(n1.seq()[3] == Base::T);
        CHECK(n1.seq()[4] == Base::N);
        CHECK(n1.end_pos() == 100);
        CHECK(!n1.is_ref());
        CHECK(!n1.belongs(0));
        CHECK(!n1.belongs(1));
        CHECK(n1.belongs(2));

        n1.set_as_ref();
        n1.set_population(3, true);
        CHECK(n1.is_ref());
        CHECK(n1.belongs(0) == true);
        CHECK(n1.belongs(1) == true);
        CHECK(n1.belongs(2) == true);
    }

}
TEST_CASE ("Graph class") {
    vargas::Graph::Node::_newID = 0;
    vargas::Graph g;

    /**
     *     GGG
     *    /   \
     * AAA     TTT
     *    \   /
     *     CCC(ref)
     */

    {
        vargas::Graph::Node n;
        n.set_endpos(3);
        n.set_as_ref();
        std::vector<bool> a = {0, 1, 1};
        n.set_population(a);
        n.set_seq("AAA");
        g.add_node(n);
    }

    {
        vargas::Graph::Node n;
        n.set_endpos(6);
        n.set_as_ref();
        std::vector<bool> a = {0, 0, 1};
        n.set_population(a);
        n.set_af(0.4);
        n.set_seq("CCC");
        g.add_node(n);
    }

    {
        vargas::Graph::Node n;
        n.set_endpos(6);
        n.set_not_ref();
        std::vector<bool> a = {0, 1, 0};
        n.set_population(a);
        n.set_af(0.6);
        n.set_seq("GGG");
        g.add_node(n);
    }

    {
        vargas::Graph::Node n;
        n.set_endpos(9);
        n.set_as_ref();
        std::vector<bool> a = {0, 1, 1};
        n.set_population(a);
        n.set_seq("TTT");
        n.set_af(0.3);
        g.add_node(n);
    }

    g.add_edge(0, 1);
    g.add_edge(0, 2);
    g.add_edge(1, 3);
    g.add_edge(2, 3);

    REQUIRE(g.node_map()->size() == 4);
    REQUIRE(g.prev_map().size() == 3);
    REQUIRE(g.next_map().size() == 3);

    // Check forward edges
    REQUIRE(g.next_map().at(0).size() == 2);
    REQUIRE(g.next_map().at(1).size() == 1);
    REQUIRE(g.next_map().at(2).size() == 1);
    REQUIRE(g.next_map().count(3) == 0);

    // Check prev edges
    REQUIRE(g.prev_map().count(0) == 0);
    REQUIRE(g.prev_map().at(1).size() == 1);
    REQUIRE(g.prev_map().at(2).size() == 1);
    REQUIRE(g.prev_map().at(3).size() == 2);

    SUBCASE("Proper Graph setup") {
        CHECK(num_to_seq(g.node(0).seq()) == "AAA");
        CHECK(num_to_seq(g.node(1).seq()) == "CCC");
        CHECK(num_to_seq(g.node(2).seq()) == "GGG");
        CHECK(num_to_seq(g.node(3).seq()) == "TTT");
    }

    SUBCASE("Graph ops") {
        /**
         * Adding three nodes:
         *[2]   GGG ----- GTGT  [4]
         *    /    \         \
         *  AAA    TTT      ATAT -- ACAC
         *    \    /        [5]     [6]
         *[1]  CCC(ref)
         */

        vargas::Graph b;
        {
            vargas::Graph::Node n;
            n.set_seq("GTGT");
            b.add_node(n);
        }
        {
            vargas::Graph::Node n;
            n.set_seq("ATAT");
            b.add_node(n);
        }
        {
            vargas::Graph::Node n;
            n.set_seq("ACAC");
            b.add_node(n);
        }
        b.add_node(g.node(2));
        b.add_edge(2,4);
        b.add_edge(4,5);
        b.add_edge(5,6);

        auto c = g + b;

        auto i = c.begin();

        CHECK(num_to_seq(i->seq()) == "AAA");
        ++i;

        // Order of these two don't matter
        bool mid = (i->seq_str() == "CCC") || (i->seq_str() == "GGG");
        CHECK(mid);
        ++i;
        mid = (i->seq_str() == "CCC") || (i->seq_str() == "GGG");
        CHECK(mid);
        ++i;

        CHECK(i->seq_str() == "TTT");
        ++i;
        CHECK(i->seq_str() == "GTGT");
        ++i;
        CHECK(i->seq_str() == "ATAT");
        ++i;
        CHECK(i->seq_str() == "ACAC");
        ++i;
        CHECK(i == c.end());

    }

    SUBCASE("Iterators") {
        /**         (ref)
         *     GGG   TTT
         *    /   \ /   \
         * AAA     \     CCA
         *    \   / \   /
         *     CCC   ACA
         *    (ref)
         */
        {
            vargas::Graph::Node n;
            std::vector<bool> pop = {1, 0, 0};
            n.set_population(pop);
            n.set_af(0.7);
            n.set_seq("ACA");
            n.set_endpos(9);
            n.set_not_ref();
            g.add_node(n);
            g.add_edge(1, 4);
            g.add_edge(2, 4);
        }
        {
            vargas::Graph::Node n;
            std::vector<bool> pop = {1, 1, 1};
            n.set_population(pop);
            n.set_af(1);
            n.set_seq("CCA");
            n.set_endpos(12);
            n.set_as_ref();
            g.add_node(n);
            g.add_edge(3, 5);
            g.add_edge(4, 5);
        }

        SUBCASE("Filtering Iterator") {
            vargas::Graph::Population filter(3, false);
            filter.set(2);
            vargas::Graph g2(g, filter);
            auto i = g2.begin();
            CHECK(num_to_seq((*i).seq()) == "AAA");
            ++i;
            CHECK(num_to_seq((*i).seq()) == "CCC");
            ++i;
            CHECK(num_to_seq((*i).seq()) == "TTT");
            ++i;
            CHECK(num_to_seq((*i).seq()) == "CCA");
            ++i;
            CHECK(i == g2.end());
        }

        SUBCASE("Filtering Ierator #2") {
            vargas::Graph::Population filter(3, false);
            filter.set(2);
            filter.set(1);
            vargas::Graph g2(g, filter);
            auto i = g2.begin();
            CHECK(num_to_seq((*i).seq()) == "AAA");
            ++i;
            // Order of these two don't matter
            bool mid = (num_to_seq((*i).seq()) == "CCC") || (num_to_seq((*i).seq()) == "GGG");
            CHECK(mid);
            ++i;
            mid = (num_to_seq((*i).seq()) == "CCC") || (num_to_seq((*i).seq()) == "GGG");
            CHECK(mid);

            ++i;
            CHECK(num_to_seq((*i).seq()) == "TTT");
            ++i;
            CHECK(num_to_seq((*i).seq()) == "CCA");
            ++i;
            CHECK(i == g2.end());
        }

        SUBCASE("Filtering Ierator: REF") {
            vargas::Graph g2(g, vargas::Graph::Type::REF);
            auto i = g2.begin();
            CHECK(num_to_seq((*i).seq()) == "AAA");
            ++i;
            CHECK(num_to_seq((*i).seq()) == "CCC");
            ++i;
            CHECK(num_to_seq((*i).seq()) == "TTT");
            ++i;
            CHECK(num_to_seq((*i).seq()) == "CCA");
            ++i;
            CHECK(i == g2.end());
        }

        SUBCASE("Filtering Ierator: MAXAF") {
            vargas::Graph g2(g, vargas::Graph::Type::MAXAF);
            auto i = g2.begin();
            CHECK(num_to_seq((*i).seq()) == "AAA");
            ++i;
            CHECK(num_to_seq((*i).seq()) == "GGG");
            ++i;
            CHECK(num_to_seq((*i).seq()) == "ACA");
            ++i;
            CHECK(num_to_seq((*i).seq()) == "CCA");
            ++i;
            CHECK(i == g2.end());
        }
    }


    SUBCASE("Graph iterator") {
        // Node visit order should be topological
        auto i = g.begin();

        CHECK(num_to_seq((*i).seq()) == "AAA");
        ++i;

        // Order of these two don't matter
        bool mid = (num_to_seq((*i).seq()) == "CCC") || (num_to_seq((*i).seq()) == "GGG");
        CHECK(mid);
        ++i;
        mid = (num_to_seq((*i).seq()) == "CCC") || (num_to_seq((*i).seq()) == "GGG");
        CHECK(mid);
        ++i;

        CHECK(num_to_seq((*i).seq()) == "TTT");
        ++i;
        CHECK(i == g.end());
        ++i;
        CHECK(i == g.end());
    }

    SUBCASE("Derived Graph") {
        std::vector<bool> filter = {0, 0, 1};
        vargas::Graph g2(g, filter);

        CHECK(g2.node_map()->size() == 4);
        CHECK(&(*g.node_map()) == &(*g2.node_map())); // Underlying node map unchanged
        CHECK(g2.next_map().size() == 2);
        CHECK(g2.prev_map().size() == 2);

        CHECK(g2.next_map().at(0).size() == 1);
        CHECK(g2.next_map().at(1).size() == 1);
        CHECK(g2.next_map().count(2) == 0); // This node shouldn't be included
        CHECK(g2.next_map().count(3) == 0);
        CHECK(g2.prev_map().count(0) == 0);
        CHECK(g2.prev_map().at(1).size() == 1);
        CHECK(g2.prev_map().at(3).size() == 1);
    }

    SUBCASE("REF graph") {
        vargas::Graph g2(g, vargas::Graph::Type::REF);
        auto iterator = g2.begin();

        CHECK((*iterator).seq_str() == "AAA");
        ++iterator;
        CHECK((*iterator).seq_str() == "CCC");
        ++iterator;
        CHECK((*iterator).seq_str() == "TTT");
        ++iterator;
        CHECK(iterator == g2.end());
    }

    SUBCASE("MAXAF graph") {
        vargas::Graph g2(g, vargas::Graph::Type::MAXAF);
        vargas::Graph::const_iterator iterator(g2);

        CHECK((*iterator).seq_str() == "AAA");
        ++iterator;
        CHECK((*iterator).seq_str() == "GGG");
        ++iterator;
        CHECK((*iterator).seq_str() == "TTT");
        ++iterator;
        CHECK(iterator == g2.end());

    }

    SUBCASE("Subgraph") {
        auto g2 = g.subgraph(2, 8);
        auto iterator = g2.begin();
        CHECK(iterator->seq_str() == "AA");
        ++iterator;
        CHECK(iterator->seq_str() == "CCC");
        ++iterator;
        CHECK(iterator->seq_str() == "GGG");
        ++iterator;
        CHECK(iterator->seq_str() == "TT");
        ++iterator;
        CHECK(iterator == g2.end());

    }

}
TEST_CASE ("Graph Factory") {
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
        << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
    }

    SUBCASE("Empty VCF") {
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
            << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2" << endl;
        }
        vargas::GraphFactory gb(tmpfa);
        gb.open_vcf(tmpvcf);
        gb.set_region("x:0-5");
        auto g = gb.build();
        auto giter = std::begin(g);
        CHECK(giter->seq_str() == "CAAAT");
        ++giter;
        CHECK(giter == std::end(g));

    }

    SUBCASE("File write wrapper") {

        SUBCASE("Basic Graph") {
            vargas::GraphFactory gb(tmpfa);
            gb.open_vcf(tmpvcf);
            gb.set_region("x:0-15");

            vargas::Graph g;
            gb.build(g);

            SUBCASE("Forward iterator") {
                auto giter = g.begin();

                CHECK((*giter).seq_str() == "CAAATAAG");
                CHECK((*giter).belongs(0) == true); // its a ref
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
                CHECK(!(*giter).belongs(0));
                CHECK(!(*giter).belongs(1));
                CHECK(!(*giter).belongs(2));
                CHECK((*giter).belongs(3));

                ++giter;
                CHECK((*giter).seq_str() == "C");
                CHECK(giter->is_ref());

                ++giter;
                CHECK((*giter).seq_str() == "CCCCCCC");
                CHECK(!giter->is_ref());

                ++giter;
                CHECK(giter->seq_str() == "");

                ++giter;
                CHECK(giter->seq_str() == "TTGGA");

                ++giter;
                CHECK(giter == g.end());


            }

            SUBCASE("Reverse iterator") {
                auto giter = g.rbegin();
                CHECK(giter->seq_str() == "TTGGA");
                ++giter;

                CHECK(giter->seq_str() == "");
                ++giter;

                CHECK(!giter->is_ref());
                CHECK(giter->seq_str() == "CCCCCCC");
                ++giter;

                CHECK(giter->is_ref());
                CHECK(giter->seq_str() == "C");
                ++giter;

                CHECK(giter->belongs(3));
                CHECK(!giter->belongs(2));
                CHECK(!giter->belongs(1));
                CHECK(!giter->belongs(0));
                CHECK(!giter->is_ref());
                CHECK(giter->seq_str() == "T");
                ++giter;

                CHECK(giter->seq_str() == "C");
                ++giter;

                CHECK(giter->seq_str() == "A");
                ++giter;

                CHECK(giter->seq_str() == "G");
                ++giter;

                CHECK(giter->is_ref());
                CHECK(giter->belongs(0) == true);
                CHECK(giter->seq_str() == "CAAATAAG");

                ++giter;
                CHECK(giter == g.rend());
                ++giter;
                CHECK(giter == g.rend());
            }

        }

        SUBCASE("Deriving a Graph") {
            vargas::GraphFactory gb(tmpfa);
            gb.open_vcf(tmpvcf);
            gb.set_region("x:0-15");

            vargas::Graph g;
            gb.build(g);

            std::vector<bool> filter = {0, 0, 0, 1};
            vargas::Graph g2(g, filter);
            auto iterator = g2.begin();

            CHECK((*iterator).seq_str() == "CAAATAAG");
            ++iterator;
            CHECK((*iterator).seq_str() == "T");
            ++iterator;
            CHECK((*iterator).seq_str() == "CCCCCCC");

        }

        SUBCASE("Sample filtering") {
            vargas::GraphFactory gb(tmpfa);
            CHECK_THROWS(gb.add_sample_filter("-"));
            gb.open_vcf(tmpvcf);
            gb.set_region("x:0-15");

            SUBCASE("No inversion") {
                gb.add_sample_filter("s1");
                auto g = gb.build();

                auto giter = g.begin();
                CHECK((*giter).seq_str() == "CAAATAAG");
                CHECK((*giter).belongs(0) == true); // its a ref
                CHECK((*giter).is_ref());

                ++giter;
                CHECK((*giter).seq_str() == "G");

                ++giter;
                CHECK((*giter).seq_str() == "A");

            }

            SUBCASE("No inversion 2") {
                gb.add_sample_filter("s2");
                auto g = gb.build();

                auto giter = g.begin();
                CHECK((*giter).seq_str() == "CAAATAAG");
                CHECK((*giter).belongs(0) == true); // its a ref
                CHECK((*giter).is_ref());

                ++giter;
                CHECK((*giter).seq_str() == "G");

                ++giter;
                CHECK((*giter).seq_str() == "C");

                ++giter;
                CHECK((*giter).seq_str() == "T");
            }

            SUBCASE("Inversion") {
                gb.add_sample_filter("s2", true);
                auto g = gb.build();

                auto giter = g.begin();
                CHECK((*giter).seq_str() == "CAAATAAG");
                CHECK((*giter).belongs(0) == true); // its a ref
                CHECK((*giter).is_ref());

                ++giter;
                CHECK((*giter).seq_str() == "G");

                ++giter;
                CHECK((*giter).seq_str() == "A");

            }

            SUBCASE("Inversion 2") {
                gb.add_sample_filter("s1", true);
                auto g = gb.build();

                auto giter = g.begin();
                CHECK((*giter).seq_str() == "CAAATAAG");
                CHECK((*giter).belongs(0) == true); // its a ref
                CHECK((*giter).is_ref());

                ++giter;
                CHECK((*giter).seq_str() == "G");

                ++giter;
                CHECK((*giter).seq_str() == "C");

                ++giter;
                CHECK((*giter).seq_str() == "T");
            }
        }

    }
    remove(tmpfa.c_str());
    remove(tmpvcf.c_str());
    remove((tmpfa + ".fai").c_str());
}

TEST_SUITE_END();
