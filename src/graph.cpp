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

  _fa.close();
  _vf.close();
  return g;
}

std::shared_ptr<vargas::graph> vargas::GraphBuilder::rebuild() {
  if (g == nullptr) g.reset();
  return build();
}
void vargas::GraphBuilder::ingroup(int percent) {
  if (percent < 0 || percent > 100) return;
  _ingroup = percent;
}





