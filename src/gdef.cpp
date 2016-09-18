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

Vargas::GraphManager::GraphManager(std::string gdef_file) {
    if (!open(gdef_file)) throw std::invalid_argument("Invalid GDEF file \"" + gdef_file + "\"");
}


void Vargas::GraphManager::close() {
    _subgraph_filters.clear();
    _subgraphs.clear();
    _ref_file = _variant_file = _region = "";
}


bool Vargas::GraphManager::open(std::string file_name, bool build_base) {
    if (file_name.length() == 0) open(std::cin);
    std::ifstream in(file_name);
    if (!in.good()) return false;
    return open(in, build_base);
}


bool Vargas::GraphManager::open(std::istream &in, bool build_base) {
    close();
    std::string line;

    // Check file type, get next line
    if (!std::getline(in, line) || line != GDEF_FILE_MARKER || !std::getline(in, line)) return false;

    // Pull meta info
    {
        std::vector<std::string> meta_split = split(line, GDEF_DELIM);
        std::vector<std::string> tv_pair;
        for (const std::string &tv : meta_split) {
            split(tv, GDEF_ASSIGN, tv_pair);
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
    int nsamps;
    {
        GraphBuilder gb(_ref_file);
        nsamps = gb.open_vcf(_variant_file);
        if (_sample_filter != "-") {
            nsamps = gb.add_sample_filter(_sample_filter, _invert_filter);
        }


        gb.region(_region);
        gb.node_len(_node_len);

        if (build_base) _subgraphs[GDEF_BASEGRAPH] = std::make_shared<Graph>(gb.build());
    }

    // subgraphs
    {
        std::vector<std::string> p_pair;
        Graph::Population pop(nsamps);
        while (std::getline(in, line)) {
            split(line, GDEF_ASSIGN, p_pair);

            if (p_pair.size() != 2) {
                throw std::invalid_argument("Invalid token: \"" + line + "\"");
            }

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


std::shared_ptr<const Vargas::Graph> Vargas::GraphManager::make_subgraph(std::string label) {
    if (!_subgraphs.count(GDEF_BASEGRAPH)) throw std::invalid_argument("No base graph built.");
    if (label == GDEF_BASEGRAPH) return base();
    label = GDEF_BASEGRAPH + GDEF_SCOPE + label;

    if (_ends_with(label, GDEF_REFGRAPH)) return make_ref(label);
    if (_ends_with(label, GDEF_MAXAFGRAPH)) return make_maxaf(label);

    if (_subgraphs.count(label)) return _subgraphs.at(label);

    if (!_subgraph_filters.count(label)) throw std::invalid_argument("Label \"" + label + "\" does not exist.");
    auto sub = std::make_shared<const Graph>(*(_subgraphs[GDEF_BASEGRAPH]), _subgraph_filters.at(label));
    #pragma omp critical(_gdef_make_subgraph)
    {
        _subgraphs[label] = sub;
    }
    return _subgraphs.at(label);
}


std::shared_ptr<const Vargas::Graph> Vargas::GraphManager::subgraph(std::string label) const {
    if (label == GDEF_BASEGRAPH) return base();
    label = GDEF_BASEGRAPH + GDEF_SCOPE + label;
    if (_subgraphs.count(label) == 0) return nullptr;
    return _subgraphs.at(label);
}


std::shared_ptr<const Vargas::Graph> Vargas::GraphManager::make_ref(std::string const &label) {
    if (!_subgraphs.count(GDEF_BASEGRAPH)) throw std::invalid_argument("No base graph built.");
    if (_subgraphs.count(label)) return _subgraphs.at(label);
    std::string root = label.substr(0, label.length() - GDEF_REFGRAPH.length() - 1);
    #pragma omp critical(_gdef_make_subgraph)
    {
        _subgraphs[label] = std::make_shared<const Graph>(*make_subgraph(root), Graph::REF);
    }
    return _subgraphs[label];
}


std::shared_ptr<const Vargas::Graph> Vargas::GraphManager::make_maxaf(std::string const &label) {
    if (!_subgraphs.count(GDEF_BASEGRAPH)) throw std::invalid_argument("No base graph built.");
    if (_subgraphs.count(label)) return _subgraphs.at(label);
    std::string root = label.substr(0, label.length() - GDEF_REFGRAPH.length() - 1);

    #pragma omp critical(_gdef_make_subgraph)
    {
        _subgraphs[label] = std::make_shared<const Graph>(*make_subgraph(root), Graph::MAXAF);
    }
    return _subgraphs[label];
}


std::shared_ptr<const Vargas::Graph> Vargas::GraphManager::base() const {
    if (!_subgraphs.count(GDEF_BASEGRAPH)) return nullptr;
    return _subgraphs.at(GDEF_BASEGRAPH);
}


Vargas::Graph::Population Vargas::GraphManager::filter(std::string label) const {
    label = GDEF_BASEGRAPH + GDEF_SCOPE + label;
    if (!_subgraph_filters.count(label)) throw std::invalid_argument("Label \"" + label + "\" does not exist.");
    return _subgraph_filters.at(label);
}


bool Vargas::GraphManager::write(std::string ref_file,
                                 std::string variant_file,
                                 std::string region,
                                 const std::string &defs,
                                 int node_len,
                                 std::string out_file,
                                 bool build_base) {
    if (out_file.length() == 0)
        return write(ref_file, variant_file, region, defs, node_len, std::cout, build_base);

    std::ofstream out(out_file);
    if (!out.good()) throw std::invalid_argument("Invalid output file: \"" + out_file + "\".");
    return write(ref_file, variant_file, region, defs, node_len, out, build_base);
}


bool Vargas::GraphManager::write(std::string ref_file,
                                 std::string variant_file,
                                 std::string region,
                                 std::string defs_str,
                                 int node_len,
                                 std::ostream &out,
                                 bool build_base,
                                 int nsamps) {

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
    std::vector<std::string> defs = split(defs_str, GDEF_DELIM);

    // Get number of samples from VCF file
    if (nsamps == 0) {
        GraphBuilder gb(ref_file);
        gb.open_vcf(variant_file);
        nsamps = gb.add_sample_filter(_sample_filter, _invert_filter);
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

        for (auto def : defs) {
            split(def, GDEF_ASSIGN, pair);
            if (pair.size() != 2) throw std::invalid_argument("Invalid assignment: \"" + def + "\".");

            pair[0] = GDEF_BASEGRAPH + GDEF_SCOPE + pair[0];
            parent_end = pair[0].find_last_of(GDEF_SCOPE);
            parent = pair[0].substr(0, parent_end);

            if (populations.count(parent) == 0)
                throw std::invalid_argument("Parent \"" + parent + "\" not yet defined.");

            if (pair[0].at(parent_end + 1) == '~')
                throw std::invalid_argument("Negative graphs cannot be defined explicitly: \"" + def + "\".");

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
                int k = 0;
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

    for (auto &p : populations) {
        out_str += p.first + GDEF_ASSIGN + p.second.to_string() + "\n";
    }

    out << out_str << std::flush;
    std::istringstream instream{out_str};
    instream.seekg(0);
    open(instream, build_base);

    return true;
}


std::string Vargas::GraphManager::to_DOT(std::string name) const {
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