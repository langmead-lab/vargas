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

std::shared_ptr<vargas::graph> vargas::GraphBuilder::build() {
  if (g != nullptr) return g;
  _fa.open(_fa_file);
  _vf.open(_vf_file);
  if (!_fa.good()) throw std::invalid_argument("Invalid FASTA file: " + _fa.file());
  if (!_vf.good()) throw std::invalid_argument("Invalid B/VCF file: " + _vf.file());
  g = std::make_shared<graph>();

  _vf.create_ingroup(_ingroup);

  // If no region is specified, the default is the first sequence in the FASTA file
  if (_chr.length() == 0) {
    _chr = _fa.sequences()[0];
    _vf.set_region(_chr, 0, 0);
  }

  int curr = _vf.region_lower(); // The graph has been built up to this position, exclusive
  graph::Node n; // Node to be added to graph
  std::vector<bool> pop(_vf.num_samples() * 2); // Population subset that has the node. *2 for genotypes
  std::vector<int> prev_unconnected; // ID's of nodes at the end of the graph left unconnected
  std::vector<int> curr_unconnected; // ID's of nodes added that are unconnected

  while (_vf.next()) {

    curr = _build_linear(*g, prev_unconnected, curr_unconnected, curr, _vf.pos());

    // Add variant nodes
    // Positions of variant nodes are referenced to ref node
    curr += _vf.ref().length();
    n.setEndPos(curr - 1);

    // ref node
    n.set_seq(_vf.ref());
    n.set_as_ref();
    curr_unconnected.push_back(g->add_node(n));

    //alt nodes
    n.set_not_ref();
    for (int i = 1; i < _vf.alleles().size(); ++i) {
      const std::string &allele = _vf.alleles()[i];
      n.set_seq(allele);
      n.setPopulation(_vf.allele_pop(allele));
      curr_unconnected.push_back(g->add_node(n));
    }
    n.setPopulation({});

    _build_edges(*g, prev_unconnected, curr_unconnected);

  }
  // Nodes after last variant
  _build_linear(*g, prev_unconnected, curr_unconnected, curr, _max_pos);


  _fa.close();
  _vf.close();
  return g;
}

void vargas::GraphBuilder::_build_edges(vargas::graph &g, std::vector<int> &prev, std::vector<int> &curr) {
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
  graph::Node n;
  n.set_as_ref();
  if (target <= 0) target = _fa.seq_len(_chr);
  while (pos < target) {

    if (pos + _max_node_len >= target) {
      n.set_seq(_fa.subseq(_chr, pos, target - 1));
      pos = target;
    }

    else {
      n.set_seq(_fa.subseq(_chr, pos, pos + _max_node_len - 1));
      pos += _max_node_len;
    }

    n.setEndPos(pos - 1);
    n.set_as_ref();
    curr.push_back(g.add_node(n));
    _build_edges(g, prev, curr);
  }
  return pos;
}

std::shared_ptr<vargas::graph> vargas::GraphBuilder::rebuild() {
  if (g == nullptr) g.reset();
  return build();
}
void vargas::GraphBuilder::ingroup(int percent) {
  if (percent < 0 || percent > 100) return;
  _ingroup = percent;
}












