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
#include "../include/vcfstream.h"
#include <stdexcept>

namespace vmatch {

//TODO Treat each haploid as an individual

class Graph {

 public:

  /** Struct defining the main parameters and options used when building a graph **/
  struct GraphParams {
    uint32_t maxNodeLen = 50000; // Maximum length of a single node graph
    int32_t ingroup = 100; // percent of individuals to include
    std::string region = ""; // Specify a specific region to build, a:b
    std::string buildfile = ""; // existing buildfile to build the complement
    bool genComplement = false; // Generate a graph with individuals not in buildfile
    bool maxAF = false; // Linear graph with only the variants/ref with maximum allele frequency
    int32_t match = 2, mismatch = 2; // default scores
    uint8_t gap_open = 3, gap_extension = 1; // default gap scores
    int8_t *nt_table = gssw_create_nt_table(); // Table of nt mappings
    int8_t *mat = gssw_create_score_matrix(match, mismatch); // table of scores
  };

  /** Empty graph uses default parameters **/
  Graph() { }

  /** Create a graph with given parameters **/
  Graph(GraphParams p) {
    if (p.nt_table == NULL) p.nt_table = gssw_create_nt_table();
    if (p.mat == NULL) p.mat = gssw_create_score_matrix(p.match, p.mismatch);
    params = p;
  }

  Graph(std::string refFile, std::string vcfFile, std::string buildFile) {
    std::ifstream ref(refFile);
    std::ofstream buildOut(buildFile);
    if (!ref.good() || !buildOut.good()) throw std::invalid_argument("Error opening files.");
    vcfstream vcf(vcfFile);
    buildGraph(ref, vcf, buildOut);
    buildOut.close();
    std::ifstream buildIn(buildFile);
    buildGraph(buildIn);
  }

  /** Build a graph from a buildfile **/
  Graph(std::string buildfile) {
    std::ifstream build(buildfile);
    if (!build.good()) throw std::invalid_argument("Error opening buildfile.");
    buildGraph(build);
  }

  /** Delete gssw graph on destruction **/
  ~Graph() {
    if (graph != NULL) gssw_graph_destroy(graph);
  }

  /** Export the graph as a DOT representation **/
  void exportDot(std::string file) {
    std::ofstream out(file);
    exportDOT(out);
  }
  void exportDOT(std::ostream &out) const;

  /** Build a graph from a vcf and variant for a buildfile **/
  void buildGraph(std::istream &reference, vcfstream &variants, std::ostream &buildout);
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
  void generateIngroup(vcfstream &variants);
};

}
#endif //VMATCH_GRAPH_H
