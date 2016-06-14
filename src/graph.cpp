/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date November 20, 2015
 *
 * vargas::Graph is a DAG representation of a reference and its variants.
 * The class wraps a gssw_graph from gssw and provides a way to construct
 * graphs from a FASTA and a VCF file with various options.
 *
 * GSSW was originally written by Erik Garrison and was moderately modified.
 *
 * graph.cpp
 */

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "../include/graph.h"


uint32_t vargas::Graph::Node::_newID = 0;


vargas::Graph::Graph(const vargas::Graph &g,
                     const Population &filter) {

  _IDMap = g._IDMap;
  _pop_size = g.pop_size();

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


vargas::Graph::Graph(const Graph &g,
                     Type type) {
  _IDMap = g._IDMap;
  _pop_size = g.pop_size();
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


void vargas::Graph::_build_derived_edges(const vargas::Graph &g,
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

uint32_t vargas::Graph::add_node(Node &n) {
  if (_IDMap->find(n.id()) != _IDMap->end()) return 0; // make sure node isn't duplicate
  if (_IDMap->size() == 0) _root = n.id(); // first node added is default root

  _IDMap->emplace(n.id(), std::make_shared<Node>(n));
  _add_order.push_back(n.id());
  return n.id();
}


bool vargas::Graph::add_edge(uint32_t n1,
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

std::string vargas::Graph::to_DOT(std::string name) const {
  std::stringstream dot;
  dot << "// Each node has the sequence, followed by end_pos,allele_freq\n";
  dot << "digraph " << name << " {\n";
  for (auto n : *_IDMap) {
    dot << n.second->id() << "[label=\"" << n.second->seq_str()
        << "\nPOS:" << n.second->end() << ", AF:" << n.second->freq() << ", REF:" << n.second->is_ref()
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


void vargas::GraphBuilder::build(vargas::Graph &g) {
  g = vargas::Graph();
  _fa.open(_fa_file);
  _vf.open(_vf_file);
  if (!_fa.good()) throw std::invalid_argument("Invalid FASTA file: " + _fa.file());
  if (!_vf.good()) throw std::invalid_argument("Invalid B/VCF file: " + _vf.file());

  _vf.create_ingroup(_ingroup);

  // If no region is specified, the default is the first sequence in the FASTA file
  if (_vf.region_chr().length() == 0) {
    _vf.set_region(_fa.sequences()[0] + ":0-0");
  }

  int curr = _vf.region_lower(); // The Graph has been built up to this position, exclusive
  std::vector<uint32_t> prev_unconnected; // ID's of nodes at the end of the Graph left unconnected
  std::vector<uint32_t> curr_unconnected; // ID's of nodes added that are unconnected

  g.set_popsize(_vf.samples().size() * 2);
  const Graph::Population all_pop(g.pop_size(), true);

  while (_vf.next()) {
    _vf.genotypes();
    auto &af = _vf.frequencies();

    curr = _build_linear_ref(g, prev_unconnected, curr_unconnected, curr, _vf.pos());

    // Add variant nodes
    // Positions of variant nodes are referenced to ref node
    curr += _vf.ref().length();

    // ref pos
    {
      Graph::Node n;
      n.set_endpos(curr - 1);
      n.set_seq(_vf.ref());
      n.set_as_ref();
      n.set_af(af[0]);
      curr_unconnected.push_back(g.add_node(n));
    }

    //alt nodes
    for (size_t i = 1; i < _vf.alleles().size(); ++i) {

      const std::string &allele = _vf.alleles()[i];
      auto allele_split = _split_seq(allele);
      Graph::Population pop(_vf.allele_pop(allele));
      if (pop && all_pop) { // Only add if someone has the allele
        for (auto &part : allele_split) {
          Graph::Node n;
          n.set_population(pop);
          n.set_seq(part);
          n.set_af(af[i]);
          n.set_not_ref();
          curr_unconnected.push_back(g.add_node(n));
          _build_edges(g, prev_unconnected, curr_unconnected);
        }
      }
    }


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

void vargas::GraphBuilder::_build_edges(vargas::Graph &g,
                                        std::vector<uint32_t> &prev,
                                        std::vector<uint32_t> &curr) {
  for (uint32_t pID : prev) {
    for (uint32_t cID : curr) {
      g.add_edge(pID, cID);
    }
  }
  prev = curr;
  curr.clear();
}


int vargas::GraphBuilder::_build_linear_ref(Graph &g,
                                            std::vector<uint32_t> &prev,
                                            std::vector<uint32_t> &curr,
                                            uint32_t pos,
                                            uint32_t target) {

  if (target == 0) target = _fa.seq_len(_vf.region_chr());
  if (pos == target) return target; // For adjacent var positions
  auto split_seq = _split_seq(_fa.subseq(_vf.region_chr(), pos, target - 1));
  for (auto s : split_seq) {
    Graph::Node n;
    n.set_as_ref();
    n.set_seq(s);
    pos += s.length();
    n.set_endpos(pos - 1);
    curr.push_back(g.add_node(n));
    _build_edges(g, prev, curr);
  }
  return target;
}

std::vector<std::string> vargas::GraphBuilder::_split_seq(std::string seq) {
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


void vargas::GraphBuilder::ingroup(int percent) {
  if (percent < 0 || percent > 100) return;
  _ingroup = percent;
}
