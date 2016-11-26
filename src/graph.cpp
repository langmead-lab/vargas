/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date November 20, 2015
 *
 * @brief
 * Vargas::Graph is a DAG representation of a reference and its variants.
 *
 * @file
 */

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "graph.h"


size_t vargas::Graph::Node::_newID = 0;

std::string vargas::Graph::GID::to_string() const {
    std::ostringstream ss;
    ss << (outgroup ? 'o' : 'i') << ',' << num << ',' << id << ',' << pct;
    return ss.str();
}

vargas::Graph::GID::GID(std::string s) {
    try {
        std::vector<std::string> s_split = split(s, ',');
        outgroup = s[0] == 'o';
        num = std::stoi(s_split[1]);
        id = std::stoi(s_split[2]);
        pct = s_split[3][0] == '1';
    }
    catch (std::exception &e) {
        std::cerr << "Invalid GID string format: " << s << std::endl;
        GID();
    }

}

vargas::Graph::Graph(const std::string &ref_file,
                     const std::string &vcf_file,
                     const std::string &region,
                     const int max_node_len) {
    _IDMap = std::make_shared<std::unordered_map<size_t, nodeptr>>();
    GraphFactory gb(ref_file);
    gb.open_vcf(vcf_file);
    gb.set_region(region);
    gb.node_len(max_node_len);
    gb.build(*this);
}


vargas::Graph::Graph(const vargas::Graph &g,
                     const Population &filter) {
    _IDMap = g._IDMap;
    _pop_size = g.pop_size();
    _filter = filter;

    // Add all nodes
    std::unordered_map<size_t, nodeptr> includedNodes;
    for (auto &nid : g._add_order) {
        auto &n = (*_IDMap)[nid];
        if (n->belongs(filter)) {
            includedNodes[nid] = n;
            _add_order.push_back(nid);
        }
    }

    _build_derived_edges(g, includedNodes);

    // Use the same description but add the filter we used
    _desc = g.desc() + "\n#Sample filter: ";
    _desc += filter.to_string();
    _max_node_len = g._max_node_len;
}


vargas::Graph::Graph(const Graph &g, Type type) {
    _IDMap = g._IDMap;
    _pop_size = g.pop_size();
    _filter = Population(_pop_size, true);
    std::unordered_map<size_t, nodeptr> includedNodes;

    if (type == Type::REF) {
        for (auto &nid : g._add_order) {
            auto &n = (*_IDMap)[nid];
            if (n->is_ref()) {
                includedNodes[nid] = n;
                _add_order.push_back(nid);
            }
        }
        _desc = g.desc() + "\n#filter: REF";
    } else if (type == Type::MAXAF) {
        size_t curr = g.root();
        size_t maxid, size;
        while (true) {
            includedNodes[curr] = (*g._IDMap).at(curr);
            _add_order.push_back(curr);
            if (g._next_map.count(curr) == 0) break; // end of graph
            maxid = g._next_map.at(curr).at(0);
            size = g._next_map.at(curr).size();
            for (size_t i = 1; i < size; ++i) {
                const size_t &id = g._next_map.at(curr).at(i);
                if ((*g._IDMap).at(id)->freq() > (*g._IDMap).at(maxid)->freq())
                    maxid = id;
            }
            curr = maxid;
        }
        _desc = g.desc() + "\n#filter: MAXAF";
    }

    _build_derived_edges(g, includedNodes);
    _max_node_len = g._max_node_len;

}


void vargas::Graph::_build_derived_edges(const vargas::Graph &g,
                                         const std::unordered_map<size_t, nodeptr> &includedNodes) {
    // Add all edges for included nodes
    for (auto &n : includedNodes) {
        if (g._next_map.count(n.second->id()) == 0) continue;
        for (auto &e : g._next_map.at(n.second->id())) {
            if (includedNodes.count(e)) {
                add_edge(n.second->id(), e);
            }
        }
    }

    // Set the new root
    if (includedNodes.find(g.root()) == includedNodes.end()) {
        throw std::invalid_argument("Currently the root must be common to all graphs.");
    }
    _root = g.root();
}


size_t vargas::Graph::add_node(const Node &n) {
    if (_IDMap->find(n.id()) != _IDMap->end()) return 0; // make sure node isn't duplicate
    if (_IDMap->size() == 0) _root = n.id(); // first node added is default root

    _IDMap->emplace(n.id(), std::make_shared<Node>(n));
    _add_order.push_back(n.id());
    return n.id();
}


bool vargas::Graph::add_edge(const size_t n1,
                             const size_t n2) {
    // Check if the nodes exist
    if (_IDMap->count(n1) == 0 || _IDMap->count(n2) == 0) return false;

    // init if first edge to be added
    if (_next_map.count(n1) == 0) {
        _next_map[n1] = std::vector<size_t>();
    }
    if (_prev_map.count(n2) == 0) {
        _prev_map[n2] = std::vector<size_t>();
    }
    _next_map[n1].push_back(n2);
    _prev_map[n2].push_back(n1);
    return true;
}


std::string vargas::Graph::to_DOT(std::string name) const {
    std::ostringstream dot;
    dot << "// Each node has the sequence, followed by end_pos,allele_freq\n";
    dot << "digraph " << name << " {\n";
    for (const auto n : *_IDMap) {
        dot << n.second->id() << "[label=\"" << n.second->seq_str()
            << "\nP:" << n.second->end_pos() << ", F:" << n.second->freq() << ", R:" << n.second->is_ref()
            << "\n[" << n.second->individuals().to_string() << "]"
            << "\"];\n";
    }
    for (const auto &n : _next_map) {
        for (const auto e : n.second) {
            dot << n.first << " -> " << e << ";\n";
        }
    }
    dot << "}\n";
    return dot.str();
}


void vargas::GraphFactory::build(vargas::Graph &g) {
    if (!_vf) throw std::invalid_argument("No variant file opened.");
    g = vargas::Graph();
    _fa.open(_fa_file);

    auto &vf = *_vf;

    if (!_fa.good()) throw std::invalid_argument("Invalid FASTA file: " + _fa_file);

    // If no region is specified, the default is the first sequence in the FASTA file
    if (vf.region_chr().length() == 0) {
        vf.set_region(_fa.sequence_names()[0] + ":0-0");
    }

    int curr = vf.region_lower(); // The Graph has been built up to this position, exclusive
    std::unordered_set<size_t> prev_unconnected; // ID's of nodes at the end of the Graph left unconnected
    std::unordered_set<size_t> curr_unconnected; // ID's of nodes added that are unconnected
    std::unordered_map<size_t, size_t> chain;

    g.set_popsize(vf.num_samples());
    const Graph::Population all_pop(g.pop_size(), true);
    g.set_filter(all_pop);

    while (vf.next()) {
        auto &af = vf.frequencies();

        curr = _build_linear_ref(g, prev_unconnected, curr_unconnected, curr, vf.pos());
        assert(_fa.subseq(_vf->region_chr(), curr, curr + vf.ref().length() - 1) == vf.ref() &&
        ("Variant and FASTA Reference does not match at position " + std::to_string(curr)) != "");

        curr += vf.ref().length();

        // Add variant nodes
        // Positions of variant nodes are referenced to ref node
        // ref node
        {
            Graph::Node n;
            n.set_endpos(curr - 1);
            n.set_seq(vf.ref());
            n.set_as_ref();
            n.set_population(vf.allele_pop(vf.ref()));
            n.set_af(af[0]);
            curr_unconnected.insert(g.add_node(n));
        }

        //alt nodes
        size_t prev_split, curr_split, chain_origin;
        chain.clear();
        for (size_t i = 1; i < vf.alleles().size(); ++i) {
            const std::string &allele = vf.alleles()[i];
            if (allele == vf.ref()) continue; // Remove duplicate nodes, REF is substituted in for unknown tags
            auto allele_split = _split_seq(allele);
            Graph::Population pop(vf.allele_pop(allele));
            if (g.pop_size() == 1 || (pop && all_pop)) { // Only add if someone has the allele. == 1 for KSNP
                {
                    Graph::Node n;
                    n.set_endpos(curr - 1);
                    n.set_population(pop);
                    n.set_seq(allele_split[0]);
                    n.set_af(af[i]);
                    n.set_not_ref();
                    prev_split = chain_origin = g.add_node(n);
                    curr_unconnected.insert(prev_split);
                }
                for (size_t i = 1; i < allele_split.size(); ++i) {
                    Graph::Node n;
                    n.set_endpos(curr - 1);
                    n.set_population(pop);
                    n.set_seq(allele_split[i]);
                    n.set_af(af[i]);
                    n.set_not_ref();
                    curr_split = g.add_node(n);
                    g.add_edge(prev_split, curr_split);
                    prev_split = curr_split;
                    chain[chain_origin] = prev_split;
                }
            }
        }
        _build_edges(g, prev_unconnected, curr_unconnected, &chain);

    }
    // Nodes after last variant
    _build_linear_ref(g, prev_unconnected, curr_unconnected, curr, vf.region_upper());

    std::string desc = "#Reference: " + _fa.file();
    desc += "\n#Region: " + vf.region_chr() + ":" + std::to_string(vf.region_lower()) + "-"
    + std::to_string(vf.region_upper());
    g.set_desc(desc);
    g.node_len(_max_node_len);

    _fa.close();
    _vf.reset();
}


void vargas::GraphFactory::_build_edges(vargas::Graph &g,
                                        std::unordered_set<size_t> &prev,
                                        std::unordered_set<size_t> &curr,
                                        std::unordered_map<size_t, size_t> *chain) {
    for (size_t pID : prev) {
        for (size_t cID : curr) {
            g.add_edge(pID, cID);
        }
    }
    prev = curr;
    curr.clear();

    if (chain) {
        for (auto &r : *chain) {
            prev.erase(r.first);
            prev.insert(r.second);
        }
    }
}


int vargas::GraphFactory::_build_linear_ref(Graph &g,
                                            std::unordered_set<size_t> &prev,
                                            std::unordered_set<size_t> &curr,
                                            size_t pos,
                                            size_t target) {

    if (target == 0) target = _fa.seq_len(_vf->region_chr());
    if (pos == target) return target; // For adjacent var positions
    auto split_seq = _split_seq(_fa.subseq(_vf->region_chr(), pos, target - 1));
    for (const auto s : split_seq) {
        Graph::Node n;
        n.pinch();
        n.set_population(g.pop_size(), true);
        n.set_as_ref();
        n.set_seq(s);
        pos += s.length();
        n.set_endpos(pos - 1);
        curr.insert(g.add_node(n));
        _build_edges(g, prev, curr);
    }
    return target;
}


std::vector<std::string> vargas::GraphFactory::_split_seq(std::string seq) {
    std::vector<std::string> split;
    if (seq.length() <= _max_node_len) {
        split.push_back(seq);
        return split;
    }
    size_t num_nodes = seq.length() / _max_node_len;
    size_t rem = seq.length() % _max_node_len;
    for (size_t i = 0; i < num_nodes; ++i) {
        split.push_back(seq.substr(i * _max_node_len, _max_node_len));
    }
    if (rem > 0) {
        split.push_back(seq.substr(num_nodes * _max_node_len, rem));
    }
    return split;

}
size_t vargas::GraphFactory::add_sample_filter(std::string filter, bool invert) {
    if (!_vf) throw std::invalid_argument("No VCF file opened, cannot add filter.");
    if (filter.length() == 0 || filter == "-") return _vf->num_samples();
    std::vector<std::string> filt;
    filter.erase(std::remove_if(filter.begin(), filter.end(), isspace), filter.end());
    if (invert) {
        const auto s = _vf->samples();
        auto vcf_samples = std::unordered_set<std::string>(s.begin(), s.end());
        const auto filter_samples = split(filter, ',');
        for (const auto &f : filter_samples) vcf_samples.erase(f);
        filt = std::vector<std::string>(vcf_samples.begin(), vcf_samples.end());
    } else {
        split(filter, ',', filt);
    }
    _vf->create_ingroup(filt);
    return _vf->num_samples();
}

size_t vargas::GraphFactory::open_vcf(std::string const &file_name) {
    _vf.reset();
    _vf = std::unique_ptr<VariantFile>(new VCF(file_name));
    if (!_vf->good()) throw std::invalid_argument("Invalid VCF/BCF file: \"" + file_name + "\"");
    return _vf->num_samples();
}


vargas::Graph::Population vargas::Graph::subset(int ingroup) const {
    vargas::Graph::Population p(_pop_size);
    for (size_t i = 0; i < _pop_size; ++i) {
        if (rand() % 100 < ingroup) p.set(i);
    }
    return p;
}

vargas::Graph vargas::Graph::subgraph(const size_t min, const size_t max) const {
    Graph ret;
    std::unordered_map<size_t, size_t> new_to_old, old_to_new;
    size_t new_id;
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
                std::vector<Base> cropped = n.seq();
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
            std::vector<Base> seq = n.seq();
            std::vector<Base> cropped(seq.begin() + min - n.begin_pos(), seq.end());
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

    ret.set_desc(_desc + "\n subgraph min:" + std::to_string(min) + " max:" + std::to_string(max));
    return ret;
}

bool vargas::Graph::validate() const {
    std::unordered_set<size_t> filled;
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

bool ::vargas::operator==(const vargas::Graph::GID &a, const vargas::Graph::GID &b) {
    if (a.outgroup != b.outgroup) return false;
    if (a.pct != b.pct) return false;
    if (a.num != b.num) return false;
    return a.id == b.id;
}

std::ostream &::vargas::operator<<(std::ostream &os, const vargas::Graph::GID &gid) {
    os << (gid.outgroup ? 'o' : 'i') << ',' << gid.num << ',' << gid.id << ',' << gid.pct;
    return os;
}

bool ::vargas::operator<(const vargas::Graph::GID &a, const vargas::Graph::GID &b) {
    if (a.outgroup != b.outgroup) return a.outgroup < b.outgroup;
    if (a.pct != b.pct) return a.pct < b.pct;
    if (a.num != b.num) return a.num < b.num;
    return a.id < b.id;
}

TEST_SUITE("Graphs");

TEST_CASE ("Node class") {
    vargas::Graph::Node::_newID = 0;
    vargas::Graph::Node n1;
    vargas::Graph::Node n2;
    CHECK(n1.id() == 0);
    CHECK(n2.id() == 1);

    SUBCASE("Node ID change") {
        n1.setID(1);
        CHECK(n1.id() == 0);
        n1.setID(2);
        CHECK(n1.id() == 2);
    }

    SUBCASE("Set Node params") {
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

    SUBCASE("Filtering iterators") {
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
        vargas::Graph::iterator i = g.begin();

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
        vargas::Graph::iterator iter(g2);

        CHECK((*iter).seq_str() == "AAA");
        ++iter;
        CHECK((*iter).seq_str() == "CCC");
        ++iter;
        CHECK((*iter).seq_str() == "TTT");
        ++iter;
        CHECK(iter == g2.end());
    }

    SUBCASE("MAXAF graph") {
        vargas::Graph g2(g, vargas::Graph::Type::MAXAF);
        vargas::Graph::iterator iter(g2);

        CHECK((*iter).seq_str() == "AAA");
        ++iter;
        CHECK((*iter).seq_str() == "GGG");
        ++iter;
        CHECK((*iter).seq_str() == "TTT");
        ++iter;
        CHECK(iter == g2.end());

    }

    SUBCASE("Subgraph") {
        auto g2 = g.subgraph(2, 8);
        auto iter = g2.begin();
        CHECK(iter->seq_str() == "AA");
        ++iter;
        CHECK(iter->seq_str() == "CCC");
        ++iter;
        CHECK(iter->seq_str() == "GGG");
        ++iter;
        CHECK(iter->seq_str() == "TT");
        ++iter;
        CHECK(iter == g2.end());

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
        << "x\t14\t.\tG\t<DUP>,<BLAH>\t99\t.\tAF=0.01,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t1|1" << endl
        << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
    }

    SUBCASE("File write wrapper") {

        SUBCASE("Basic Graph") {
            vargas::GraphFactory gb(tmpfa);
            gb.open_vcf(tmpvcf);
            gb.node_len(5);
            gb.set_region("x:0-15");

            vargas::Graph g;
            gb.build(g);

            auto giter = g.begin();

            CHECK((*giter).seq_str() == "CAAAT");
            CHECK((*giter).belongs(0) == true); // its a ref
            CHECK((*giter).is_ref());

            ++giter;
            CHECK((*giter).seq_str() == "AAG");
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
            CHECK((*giter).seq_str() == "CCCCC");
            CHECK(!giter->is_ref());
            ++giter;
            CHECK((*giter).seq_str() == "CC");
            CHECK(!giter->is_ref());

            ++giter;
            CHECK(giter->seq_str() == "");

            ++giter;
            CHECK(giter->seq_str() == "TTG");

        }

        SUBCASE("Deriving a Graph") {
            vargas::GraphFactory gb(tmpfa);
            gb.open_vcf(tmpvcf);
            gb.node_len(5);
            gb.set_region("x:0-15");

            vargas::Graph g;
            gb.build(g);

            std::vector<bool> filter = {0, 0, 0, 1};
            vargas::Graph g2(g, filter);
            auto iter = g2.begin();

            CHECK((*iter).seq_str() == "CAAAT");
            ++iter;
            CHECK((*iter).seq_str() == "AAG");
            ++iter;
            CHECK((*iter).seq_str() == "T");
            ++iter;
            CHECK((*iter).seq_str() == "CCCCC");
            ++iter;
            CHECK((*iter).seq_str() == "CC");

        }

        SUBCASE("Sample filtering") {
            vargas::GraphFactory gb(tmpfa);
            CHECK_THROWS(gb.add_sample_filter("-"));
            gb.open_vcf(tmpvcf);
            gb.node_len(5);
            gb.set_region("x:0-15");

            SUBCASE("No inversion") {
                gb.add_sample_filter("s1");
                auto g = gb.build();

                auto giter = g.begin();
                CHECK((*giter).seq_str() == "CAAAT");
                CHECK((*giter).belongs(0) == true); // its a ref
                CHECK((*giter).is_ref());

                ++giter;
                CHECK((*giter).seq_str() == "AAG");
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
                CHECK((*giter).seq_str() == "CAAAT");
                CHECK((*giter).belongs(0) == true); // its a ref
                CHECK((*giter).is_ref());

                ++giter;
                CHECK((*giter).seq_str() == "AAG");
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
                CHECK((*giter).seq_str() == "CAAAT");
                CHECK((*giter).belongs(0) == true); // its a ref
                CHECK((*giter).is_ref());

                ++giter;
                CHECK((*giter).seq_str() == "AAG");
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
                CHECK((*giter).seq_str() == "CAAAT");
                CHECK((*giter).belongs(0) == true); // its a ref
                CHECK((*giter).is_ref());

                ++giter;
                CHECK((*giter).seq_str() == "AAG");
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
