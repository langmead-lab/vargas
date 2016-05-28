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
#include "../include/graph.h"


long vargas::graph::Node::_newID = 0;


vargas::graph::graph(const vargas::graph &g,
                     const std::vector<bool> &filter) {
  _popSize = g._popSize;
  if (filter.size() != _popSize) {
    throw std::invalid_argument("Filter size should match graph population size.");
  }
  _IDMap = g._IDMap;
  std::vector<long> indexes; // indexes of the individuals that are included in the filter
  for (long i = 0; i < filter.size(); ++i) {
    if (filter[i]) {
      indexes.push_back(i);
    }
  }

  // Add all nodes
  std::unordered_map<long, nodeptr> includedNodes;
  for (auto &n : *(g._IDMap)) {
    for (long i : indexes) {
      if (n.second->belongs(i)) {
        includedNodes[n.first] = n.second;
        break;
      }
    }
  }

  // Add all edges for included nodes
  for (auto &n : includedNodes) {
    if (g._next_map.find(n.second->id()) == g._next_map.end()) continue;
    for (auto e : g._next_map.at(n.second->id())) {
      if (includedNodes.find(e) != includedNodes.end()) {
        add_edge(n.second->id(), e);
      }
    }
  }

  // Set the new root
  if (includedNodes.find(g.root()) == includedNodes.end()) {
    throw std::invalid_argument("Currently the root must be common to all graphs.");
  }
  _root = g.root();
  finalize();

  // Use the same description but add the filter we used
  _desc = g.desc() + "\nfilter: ";
  for (auto b : filter) {
    _desc += b == true ? "1" : "0";
    _desc += ",";
  }

}

void vargas::graph::finalize() {
  _toposort.clear();
  std::set<long> unmarked, tempmarked, permmarked;
  for (auto &n : *_IDMap) {
    unmarked.insert(n.first);
  }
  while (!unmarked.empty()) {
    _visit(*unmarked.begin(), unmarked, tempmarked, permmarked);
  }
  std::reverse(_toposort.begin(), _toposort.end());
}

long vargas::graph::add_node(Node &n) {
  if (_popSize < 0) _popSize = n.popSize(); // first node dictates graph population size
  if (_IDMap->find(n.id()) != _IDMap->end()) return 0; // make sure node isn't duplicate
  if (_root < 0) _root = n.id(); // first node added is default root

  _IDMap->emplace(n.id(), std::make_shared<Node>(n));
  return n.id();
}

bool vargas::graph::add_edge(long n1, long n2) {
  // Check if the nodes exist
  if (_IDMap->find(n1) == _IDMap->end() || _IDMap->find(n2) == _IDMap->end()) return false;

  // init if first edge to be added
  if (_next_map.find(n1) == _next_map.end()) {
    _next_map[n1] = std::vector<long>();
  }
  if (_prev_map.find(n2) == _prev_map.end()) {
    _prev_map[n2] = std::vector<long>();
  }
  _next_map[n1].push_back(n2);
  _prev_map[n2].push_back(n1);
  _toposort.clear(); // any ordering is invalidated
  return true;
}

void vargas::GraphBuilder::build(vargas::graph &g) {
  g = vargas::graph();
  _fa.open(_fa_file);
  _vf.open(_vf_file);
  if (!_fa.good()) throw std::invalid_argument("Invalid FASTA file: " + _fa.file());
  if (!_vf.good()) throw std::invalid_argument("Invalid B/VCF file: " + _vf.file());

  _vf.create_ingroup(_ingroup);

  // If no region is specified, the default is the first sequence in the FASTA file
  if (_chr.length() == 0) {
    _chr = _fa.sequences()[0];
    _vf.set_region(_chr, 0, 0);
  }

  int curr = _vf.region_lower(); // The graph has been built up to this position, exclusive
  std::vector<bool> pop(_vf.num_samples() * 2); // Population subset that has the node. *2 for genotypes
  std::vector<int> prev_unconnected; // ID's of nodes at the end of the graph left unconnected
  std::vector<int> curr_unconnected; // ID's of nodes added that are unconnected

  while (_vf.next()) {

    curr = _build_linear(g, prev_unconnected, curr_unconnected, curr, _vf.pos());

    // Add variant nodes
    // Positions of variant nodes are referenced to ref node
    curr += _vf.ref().length();

    // ref pos
    {
      graph::Node n;
      n.set_endpos(curr - 1);
      n.set_seq(_vf.ref());
      n.set_as_ref();
      curr_unconnected.push_back(g.add_node(n));
    }

    //alt nodes
    for (int i = 1; i < _vf.alleles().size(); ++i) {
      graph::Node n;
      n.set_not_ref();
      const std::string &allele = _vf.alleles()[i];
      n.set_seq(allele);
      _vf.genotypes();
      n.set_population(_vf.allele_pop(allele));
      curr_unconnected.push_back(g.add_node(n));
    }


    _build_edges(g, prev_unconnected, curr_unconnected);

  }
  // Nodes after last variant
  _build_linear(g, prev_unconnected, curr_unconnected, curr, _max_pos);


  _fa.close();
  _vf.close();
  g.finalize();
}


void vargas::GraphBuilder::_build_edges(vargas::graph &g,
                                        std::vector<int> &prev,
                                        std::vector<int> &curr) {
  for (int pID : prev) {
    for (int cID : curr) {
      g.add_edge(pID, cID);
    }
  }
  prev = curr;
  curr.clear();
}

int vargas::GraphBuilder::_build_linear(graph &g,
                                        std::vector<int> &prev,
                                        std::vector<int> &curr,
                                        int pos,
                                        int target) {

  if (target <= 0) target = _fa.seq_len(_chr);
  while (pos < target) {
    graph::Node n;
    n.set_as_ref();

    if (pos + _max_node_len >= target) {
      n.set_seq(_fa.subseq(_chr, pos, target - 1));
      pos = target;
    }

    else {
      n.set_seq(_fa.subseq(_chr, pos, pos + _max_node_len - 1));
      pos += _max_node_len;
    }

    n.set_endpos(pos - 1);
    n.set_as_ref();
    curr.push_back(g.add_node(n));
    _build_edges(g, prev, curr);
  }
  return pos;
}

void vargas::GraphBuilder::ingroup(int percent) {
  if (percent < 0 || percent > 100) return;
  _ingroup = percent;
}


