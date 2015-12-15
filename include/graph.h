/**
 * Ravi Gaddipati
 * November 20, 2015
 * rgaddip1@jhu.edu
 *
 * vargas::Graph is a DAG representation of a reference and its variants.
 * The class wraps a gssw_graph from gssw and provides a way to construct
 * graphs from a FASTA and a VCF file with various options.
 *
 * GSSW was originally written by Erik Garrison and was moderately modified.
 *
 * graph.h
 */

#ifndef VARGAS_GRAPH_H
#define VARGAS_GRAPH_H

#include "../gssw/src/gssw.h"
#include "readsource.h"
#include "vcfstream.h"

namespace vargas {

// Object to hold alignment results
struct Alignment {
  Read read;

  // Optimal alignment
  uint16_t optScore;
  int32_t optAlignEnd;
  int32_t optCount;

  // Suboptimal alignment
  uint16_t subOptScore;
  int32_t subOptAlignEnd;
  int32_t subOptCount;

  // alignment flag
  int8_t corflag;
};


inline std::ostream &operator<<(std::ostream &os, const Alignment &a) {
  os << a.read << ',' << a.optScore << ',' << a.optAlignEnd << ',' << a.optCount
      << ',' << a.subOptScore << ',' << a.subOptAlignEnd << ',' << a.subOptCount
      << ',' << int32_t(a.corflag);
  return os;
}


class Graph {


 public:

  /** Struct defining the main parameters and options used when building a graph **/
  struct GraphParams {
    uint32_t maxNodeLen = 50000; // Maximum length of a single node graph
    int32_t ingroup = 100; // percent of individuals to include
    std::string region = ""; // Specify a specific region to build, a:b
    std::string complementSource = ""; // existing buildfile to build the complement
    bool genComplement = false; // Generate a graph with individuals not in buildfile
    bool maxAF = false; // Linear graph with only the variants/ref with maximum allele frequency
    bool includeRefIndivs = true; // Include individuals for reference nodes
    int32_t match = 2, mismatch = 2; // default scores
    uint8_t gap_open = 3, gap_extension = 1; // default gap scores
    int8_t *nt_table = gssw_create_nt_table(); // Table of nt mappings
    int8_t *mat = gssw_create_score_matrix(match, mismatch); // table of scores

    bool includeIndividuals = false;
  };

  /** Empty graph uses default parameters **/
  Graph() { }

  /** Create a graph with given parameters **/
  Graph(GraphParams p) : params(p) {
    if (!params.nt_table) params.nt_table = gssw_create_nt_table();
    if (!params.mat) params.mat = gssw_create_score_matrix(p.match, p.mismatch);
  }

  /** Exports a buildfile and builds a graph with given fasta and vcf **/
  Graph(std::string refFile, std::string vcfFile, std::string buildFile) {
    std::ifstream ref(refFile);
    std::ofstream buildOut(buildFile);
    if (!ref.good() || !buildOut.good()) throw std::invalid_argument("Error opening files.");
    vcfstream vcf(vcfFile);
    exportBuildfile(ref, vcf, buildOut);
    buildOut.close();
    std::ifstream buildIn(buildFile);
    buildGraph(buildIn);
  }

  /** Build a graph from a buildfile **/
  Graph(std::string buildfile) {
    buildGraph(buildfile);
  }

  /** Delete gssw graph on destruction **/
  ~Graph() {
    if (graph != NULL) gssw_graph_destroy(graph);
    free(params.nt_table);
    free(params.mat);
  }

  /** Export the graph as a DOT representation **/
  void exportDOT(std::string file) {
    std::ofstream out(file);
    exportDOT(out);
    out.close();
  }
  void exportDOT(std::ostream &out) const;

  /** Export a buildfile **/
  void exportBuildfile(std::string ref, std::string vcf, std::string build = "");
  void exportBuildfile(std::istream &reference, vcfstream &variants) {
    exportBuildfile(reference, variants, std::cout);
  }
  void exportBuildfile(std::istream &reference, vcfstream &variants, std::ostream &buildout);

  /** Build the graph in memory **/
  void buildGraph(std::string build) {
    std::ifstream b(build);
    if (!b.good()) throw std::invalid_argument("Error opening file.");
    buildGraph(b);
    b.close();
  }
  void buildGraph(std::istream &buildfile);

  /** Align to the graph **/
  Alignment *align(const Read &read);
  void align(const Read &r, Alignment &a);

  /** Setters and getters **/
  gssw_graph *getGSSWGraph() const { return graph; }
  GraphParams getParamsCopy() const { return params; }
  void setParams(GraphParams p) {
    if (p.nt_table == NULL) p.nt_table = gssw_create_nt_table();
    if (p.mat == NULL) p.mat = gssw_create_score_matrix(p.match, p.mismatch);
    params = p;
  }
  void setIngroup(int32_t i) { params.ingroup = i; }
  void setComplement(bool b) { params.genComplement = b; }
  void setComplementSource(std::string s) { params.complementSource = s; }
  void setScores(int32_t m = 2, int32_t mm = 2, uint8_t go = 3, uint8_t ge = 1) {
    if (graph) std::cerr << "Graph must be rebuilt for new scores to take effect." << std::endl;
    params.match = m;
    params.mismatch = mm;
    params.gap_extension = ge;
    params.gap_open = go;
  }
  void setMaxAF(bool b) { params.maxAF = b; }
  void useIndividuals(bool b) { params.includeIndividuals = b; }

 protected:
  GraphParams params;
  gssw_graph *graph = NULL; // The graph

  /** splits a region argument **/
  void parseRegion(std::string region, uint32_t *min, uint32_t *max);
  /** generate an ingroup based on params **/
  void generateIngroup(vcfstream &variants);
};

}
#endif //VARGAS_GRAPH_H
