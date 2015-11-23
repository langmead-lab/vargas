/**
 * Ravi Gaddipati
 * November 20, 2015
 * rgaddip1@jhu.edu
 *
 * vmatch::Graph is a DAG representation of a reference and its variants.
 * The class wraps a gssw_graph from gssw and provides a way to construct
 * graphs from a FASTA and a VCF file with various options.
 *
 * GSSW was originally written by Erik Garrison and was moderately modified.
 */

#ifndef VMATCH_GRAPH_H
#define VMATCH_GRAPH_H

#include "../gssw/src/gssw.h"
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include "../include/utils.h"
#include <stdexcept>

namespace vmatch {

//TODO Treat each haploid as an individual

class Graph {

 public:

  /** Struct defining the main paramemters and options used when building a graph **/
  struct GraphParams {
    uint32_t maxNodeLen = 50000; // Maximum length of a single node graph
    int32_t ingroup = 100; // percent of individuals to include
    std::string region = ""; // Specify a specific region to build, a:b
    std::string buildfile = ""; // existing buildfile to build the complement
    bool genComplement = false; // Generate a graph with individuals not in buildfile
    bool maxAF = false; // Linear graph with only the variants/ref with maximum allele frequency
    int8_t *nt_table = NULL; // Table of nt mappings
    int8_t *mat = NULL; // table of scores
    int32_t match = 2, mismatch = 2; // default scores
    uint8_t gap_open = 3, gap_extension = 1; // default gap scores
  };

  /** Empty graph uses default parameters **/
  Graph() {
    params.nt_table = gssw_create_nt_table();
    params.mat = gssw_create_score_matrix(params.match, params.mismatch);
  }

  /** Create a graph with given parameters **/
  Graph(GraphParams p) {
    if (p.nt_table == NULL) p.nt_table = gssw_create_nt_table();
    if (p.mat == NULL) p.mat = gssw_create_score_matrix(p.match, p.mismatch);
    params = p;
  }

  /**  Build a graph from a FASTA and a VCF **/
  Graph(std::istream &reference, std::istream &vcf) {
    Graph();
    buildGraph(buildGraph(reference, vcf));
  }

  /** Build a graph from a buildfile, as made by exportBuildfile **/
  Graph(std::ifstream &buildfile) {
    Graph();
    buildGraph(buildfile);
  }

  /** Delete gssw graph on destruction **/
  ~Graph() {
    if (graph != NULL) gssw_graph_destroy(graph);
  }

  /** Export the graph as a DOT representation **/
  void exportDOT(std::ostream &out) const;

  /** Export a quickbuild file to avoid processing the FASTA and VCF again **/
  void exportBuildfile(std::ostream &out) const;
  std::iostream &buildGraph(std::istream &reference, std::istream &variants);
  void buildGraph(std::istream &buildfile);

  gssw_graph *getGSSWGraph() const {
    return graph;
  }
  GraphParams getParamsCopy() const {
    return params;
  }
  void setParams(GraphParams p) {
    if (p.nt_table == NULL) p.nt_table = gssw_create_nt_table();
    if (p.mat == NULL) p.mat = gssw_create_score_matrix(p.match, p.mismatch);
    params = p;
  }
  void setScores(int32_t m = 2, int32_t mm = 2, uint8_t go = 3, uint8_t ge = 1) {
    params.match = m;
    params.mismatch = mm;
    params.gap_extension = ge;
    params.gap_open = go;
  }

 protected:
  GraphParams params;
  gssw_graph *graph = NULL; // The graph

  void parseRegion(std::string region, uint32_t *min, uint32_t *max);

};

}
#endif //VMATCH_GRAPH_H
