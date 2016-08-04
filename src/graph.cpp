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


uint32_t Vargas::Graph::Node::_newID = 0;

Vargas::Graph::Graph(std::string ref_file, std::string vcf_file, std::string region, int max_node_len) {
    _IDMap = std::make_shared<std::unordered_map<uint32_t, nodeptr>>();
    GraphBuilder gb(ref_file, vcf_file);
    gb.region(region);
    gb.node_len(max_node_len);
    gb.build(*this);
}

Vargas::Graph::Graph(const Vargas::Graph &g,
                     const Population &filter) {
  _IDMap = g._IDMap;
  _pop_size = g.pop_size();
    _filter = filter;

  // Add all nodes
  std::unordered_map<uint32_t, nodeptr> includedNodes;
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


Vargas::Graph::Graph(const Graph &g,
                     Type type) {
  _IDMap = g._IDMap;
  _pop_size = g.pop_size();
    _filter = Population(_pop_size, true);
  std::unordered_map<uint32_t, nodeptr> includedNodes;

  if (type == REF) {
    for (auto &nid : g._add_order) {
      auto &n = (*_IDMap)[nid];
      if (n->is_ref()) {
        includedNodes[nid] = n;
        _add_order.push_back(nid);
      }
    }
    _desc = g.desc() + "\n#filter: REF";
  }

  else if (type == MAXAF) {
    uint32_t curr = g.root();
    uint32_t maxid, size;
    while (true) {
      includedNodes[curr] = (*g._IDMap).at(curr);
      _add_order.push_back(curr);
      if (g._next_map.count(curr) == 0) break; // end of graph
      maxid = g._next_map.at(curr).at(0);
      size = g._next_map.at(curr).size();
      for (uint32_t i = 1; i < size; ++i) {
        const uint32_t &id = g._next_map.at(curr).at(i);
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


void Vargas::Graph::_build_derived_edges(const Vargas::Graph &g,
                                         const std::unordered_map<uint32_t, nodeptr> &includedNodes) {
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

uint32_t Vargas::Graph::add_node(Node &n) {
  if (_IDMap->find(n.id()) != _IDMap->end()) return 0; // make sure node isn't duplicate
  if (_IDMap->size() == 0) _root = n.id(); // first node added is default root

  _IDMap->emplace(n.id(), std::make_shared<Node>(n));
  _add_order.push_back(n.id());
  return n.id();
}


bool Vargas::Graph::add_edge(uint32_t n1,
                             uint32_t n2) {
  // Check if the nodes exist
  if (_IDMap->count(n1) == 0 || _IDMap->count(n2) == 0) return false;

  // init if first edge to be added
  if (_next_map.count(n1) == 0) {
    _next_map[n1] = std::vector<uint32_t>();
  }
  if (_prev_map.count(n2) == 0) {
    _prev_map[n2] = std::vector<uint32_t>();
  }
  _next_map[n1].push_back(n2);
  _prev_map[n2].push_back(n1);
  return true;
}

std::string Vargas::Graph::to_DOT(std::string name) const {
  std::ostringstream dot;
  dot << "// Each node has the sequence, followed by end_pos,allele_freq\n";
  dot << "digraph " << name << " {\n";
  for (auto n : *_IDMap) {
    dot << n.second->id() << "[label=\"" << n.second->seq_str()
        << "\nP:" << n.second->end() << ", F:" << n.second->freq() << ", R:" << n.second->is_ref()
        << "\n[" << n.second->individuals().to_string() << "]"
        << "\"];\n";
  }
  for (auto &n : _next_map) {
    for (auto e : n.second) {
      dot << n.first << " -> " << e << ";\n";
    }
  }
  dot << "}\n";
  return dot.str();
}


void Vargas::GraphBuilder::build(Vargas::Graph &g) {
  g = Vargas::Graph();
  _fa.open(_fa_file);
  _vf.open(_vf_file);
  if (!_fa.good()) throw std::invalid_argument("Invalid FASTA file: " + _fa_file);
  if (!_vf.good()) throw std::invalid_argument("Invalid B/VCF file: " + _vf_file);

  _vf.create_ingroup(100);

  // If no region is specified, the default is the first sequence in the FASTA file
  if (_vf.region_chr().length() == 0) {
    _vf.set_region(_fa.sequence_names()[0] + ":0-0");
  }

  int curr = _vf.region_lower(); // The Graph has been built up to this position, exclusive
  std::unordered_set<uint32_t> prev_unconnected; // ID's of nodes at the end of the Graph left unconnected
  std::unordered_set<uint32_t> curr_unconnected; // ID's of nodes added that are unconnected
  std::unordered_map<uint32_t, uint32_t> chain;

  g.set_popsize(_vf.samples().size() * 2);
  const Graph::Population all_pop(g.pop_size(), true);
    g.set_filter(all_pop);

  while (_vf.next()) {
    _vf.genotypes();
    auto &af = _vf.frequencies();

    curr = _build_linear_ref(g, prev_unconnected, curr_unconnected, curr, _vf.pos());

    // Add variant nodes
    // Positions of variant nodes are referenced to ref node
    curr += _vf.ref().length();

    // ref node
    {
      Graph::Node n;
      n.set_endpos(curr - 1);
      n.set_seq(_vf.ref());
      n.set_as_ref();
      n.set_population(_vf.allele_pop(_vf.ref()));
      n.set_af(af[0]);
      curr_unconnected.insert(g.add_node(n));
    }

    //alt nodes
    uint32_t prev_split, curr_split, chain_origin;
    chain.clear();
    for (size_t i = 1; i < _vf.alleles().size(); ++i) {
      const std::string &allele = _vf.alleles()[i];
      if (allele == _vf.ref()) continue; // Remove duplicate nodes, REF is substituted in for unknown tags
      auto allele_split = _split_seq(allele);
      Graph::Population pop(_vf.allele_pop(allele));
      if (pop && all_pop) { // Only add if someone has the allele
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
  _build_linear_ref(g, prev_unconnected, curr_unconnected, curr, _vf.region_upper());

  std::string desc = "#Reference: " + _fa.file();
  desc += "\n#B/VCF: " + _vf.file();
  desc += "\n#Region: " + _vf.region_chr() + ":" + std::to_string(_vf.region_lower()) + "-"
      + std::to_string(_vf.region_upper());
  desc += "\n#Ingroup pct: " + _vf.ingroup_str();
  desc += "\n#Number of samples: " + std::to_string(_vf.num_samples());
  g.set_desc(desc);
  g.node_len(_max_node_len);

  _fa.close();
  _vf.close();
}

void Vargas::GraphBuilder::_build_edges(Vargas::Graph &g,
                                        std::unordered_set<uint32_t> &prev,
                                        std::unordered_set<uint32_t> &curr,
                                        std::unordered_map<uint32_t, uint32_t> *chain) {
  for (uint32_t pID : prev) {
    for (uint32_t cID : curr) {
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


int Vargas::GraphBuilder::_build_linear_ref(Graph &g,
                                            std::unordered_set<uint32_t> &prev,
                                            std::unordered_set<uint32_t> &curr,
                                            uint32_t pos,
                                            uint32_t target) {

  if (target == 0) target = _fa.seq_len(_vf.region_chr());
  if (pos == target) return target; // For adjacent var positions
  auto split_seq = _split_seq(_fa.subseq(_vf.region_chr(), pos, target - 1));
  for (auto s : split_seq) {
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

std::vector<std::string> Vargas::GraphBuilder::_split_seq(std::string seq) {
  std::vector<std::string> split;
  if (seq.length() <= _max_node_len && seq.size() > 0) {
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


Vargas::Graph::Population Vargas::Graph::subset(int ingroup) const {
  Vargas::Graph::Population p(_pop_size);
  for (size_t i = 0; i < _pop_size; ++i) {
    if (rand() % 100 < ingroup) p.set(i);
  }
  return p;
}

bool Vargas::Graph::FilteringIter::operator==(const Vargas::Graph::FilteringIter &other) const {
  if (_type == END && other._type == END) return true; // All ends are equal
  if (_type != other._type) return false; // Same type of iterator
  if (&_graph != &other._graph) return false; // same base graph
  return _currID == other._currID;
}

bool Vargas::Graph::FilteringIter::operator!=(const Vargas::Graph::FilteringIter &other) const {
  if (_type != other._type) return true;
  if (_type == other._type) return false;
  if (&_graph != &other._graph) return true;
  return _currID != other._currID;
}


Vargas::Graph::FilteringIter &Vargas::Graph::FilteringIter::operator++() {
  // If end of graph has been reached
  if (_type == END) return *this;

  if (_type == TOPO) {
    ++_currID;
    if (_currID == _add_order_size) _type = END;
    return *this;
  }

  if (_graph._next_map.count(_currID) == 0) {
    _type = END;
    return *this;
  }

  const auto &next_vec = _graph._next_map.at(_currID);
  const auto &graph_map = *(_graph._IDMap);

  switch (_type) {
    case REF:
      for (const uint32_t &nextID : next_vec) {
        if (graph_map.at(nextID)->is_ref()) {
          _insert_queue(nextID);
          break; // Assuming there is only one REF node per branch
        }
      }
          break;

    case FILTER:
      for (const uint32_t &nextID : next_vec) {
        // Add all nodes that intersect with filter
        if (graph_map.at(nextID)->belongs(_filter)) _insert_queue(nextID);
      }
          break;

    case MAXAF: {
      uint32_t max_id = graph_map.at(next_vec.at(0))->id();
      float max_af = graph_map.at(next_vec.at(0))->freq();
      float freq;
      for (size_t i = 1; i < next_vec.size(); ++i) {
        freq = graph_map.at(next_vec.at(i))->freq();
        if (freq > max_af) {
          max_af = freq;
          max_id = graph_map.at(next_vec.at(i))->id();
        }
      }
      _insert_queue(max_id);
    }
          break;

    default:
      throw std::logic_error("Invalid type.");
  }


  if (_queue.empty()) {
    _type = END;
    return *this;
  }

  _traversed.insert(_currID); // Insert previous iterator position (node id)
  _currID = _queue.front();
  _queue.pop();
  _queue_unique.erase(_currID);
  return *this;
}

const Vargas::Graph::Node &Vargas::Graph::FilteringIter::operator*() const {
  if (_type == TOPO) return *(_graph._IDMap->at(_graph._add_order.at(_currID)));
  return *(_graph._IDMap->at(_currID));
}

const std::vector<uint32_t> &Vargas::Graph::FilteringIter::incoming() {
  if (_type == TOPO) {
    const uint32_t nid = _graph._add_order.at(_currID);
    if (_graph._prev_map.count(nid) == 0) return _incoming;
    return _graph._prev_map.at(nid);
  }
  if (_graph._prev_map.count(_currID) == 0) return _incoming;
  _incoming.clear();
  for (auto &id : _graph._prev_map.at(_currID)) {
    if (_traversed.count(id)) _incoming.push_back(id);
  }
  return _incoming;
}

const std::vector<uint32_t> &Vargas::Graph::FilteringIter::outgoing() {
  if (_graph._prev_map.count(_currID) == 0) return _outgoing;
  return _graph._next_map.at(_graph._add_order.at(_currID));
}
