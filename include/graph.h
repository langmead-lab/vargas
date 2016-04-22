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
 * @file graph.h
 */

#ifndef VARGAS_GRAPH_H
#define VARGAS_GRAPH_H

#include "../gssw/src/gssw.h"
#include "readsource.h"
#include "vcfstream.h"
#include "loadfile.h"
#include "alignment.h"


namespace vargas {


class Graph {


 public:

  /**
   * A structure to store various parameters that control graph creation.
   *
   * @param maxNodeLen Limits each maximum node length. Default 50,000
   * @param ingroup Percent of individuals to include in the graph, default 100
   * @param region A string of the region to use, of the form "A:B" where
   *        A is the min position, and B is the max position.
   * @param complementSource Buildfile to generate a complement graph of
   * @param genComplement Creates a complement of complementSource if true
   * @param maxAF If true, only keeps alleles with the maximum frequency
   * @param includeRefIndivs Default true, mark reference nodes with individuals that have that allele
   * @param match Default match score 2
   * @param mismatch Default mismatch score 2
   * @param gap_open Default gap opening penalty of 3
   * @param gap_extension Default gap extension penalty of 1
   * @param *nt_table Table that maps nucleotides to integers
   * @param *mat Scoring matrix, created with match and mismatch.
   * @param inMemoryRef Load the entire reference into memory instead of streaming.
   * @param includeIndividuals Default false, loads compressed indivudual lists into memory
   */
  struct GraphParams {
    uint32_t maxNodeLen = 50000;
    int32_t ingroup = 100;
    std::string region = "";
    std::string complementSource = "";
    bool genComplement = false;
    bool maxAF = false;
    bool includeRefIndivs = true;
    int32_t match = 2, mismatch = 2;
    uint8_t gap_open = 3, gap_extension = 1;
    int8_t *nt_table = gssw_create_nt_table();
    int8_t *mat = gssw_create_score_matrix(match, mismatch);
    bool inMemoryRef = false;
    bool includeIndividuals = false;
  };


  /*********************************** CONSTRUCTORS ***********************************/

  Graph() { }

  /**
   * Creates a Graph object with a given gssw_graph.
   * @param g gssw_graph pointer
   */
  Graph(gssw_graph *g) : graph(g) { }

  /**
   * Generate a left context graph of a given alignment and graph.
   * The locality of > read length base pairs is loaded.
   * @param g An existing Graph to pull a subgraph from
   * @param a An alignment to copy a graph around.
   */
  Graph(const Graph &g, const Alignment &a);

  /**
   * Create a graph with the given parameters.
   * @param p GraphParams to use
   */
  Graph(GraphParams p) : params(p) {
    if (!params.nt_table) params.nt_table = gssw_create_nt_table();
    if (!params.mat) params.mat = gssw_create_score_matrix(p.match, p.mismatch);
  }

  /**
   * Create a graph from a given reference and VCF file.
   * A buildfile is also exported.
   * @param refFile Reference FASTA filename.
   * @param vcfFile VCF filename.
   * @param buildFile buildfile output fileName.
   */
  Graph(std::string refFile, std::string vcfFile, std::string buildFile) {
    std::istream *ref;
    if (params.inMemoryRef) {
      ref = loadFile(refFile);
    } else {
      ref = new std::ifstream(refFile);
    }
    std::ofstream buildOut(buildFile);
    if (!ref || !buildOut.good()) throw std::invalid_argument("Error opening files.");
    vcfstream vcf(vcfFile);
    exportBuildfile(ref, vcf, buildOut);
    buildOut.close();
    std::ifstream buildIn(buildFile);
    buildGraph(buildIn);
    delete ref;
  }

  /**
   * Build a graph from a given buildfile.
   * @param buildfile Buildfile filename.
   */
  Graph(std::string buildfile) {
    buildGraph(buildfile);
  }

  ~Graph() {
    // Delete gssw_graph on destruction
    if (graph != NULL) gssw_graph_destroy(graph);
    free(params.nt_table);
    free(params.mat);
  }

  /*********************************** FUNCTIONS ***********************************/

  /**
   * Exports the graph in the DOT format.
   * @param file Output filename.
   */
  void exportDOT(std::string file) {
    std::ofstream out(file);
    exportDOT(out);
    out.close();
  }

  /**
   * Exports the graph in the DOT format.
   * @param out Output stream to export DOT graph on.
   * @param name Default "vargraph", the name of the digraph.
   */
  void exportDOT(std::ostream &out, std::string name = "vargraph") const;

  /**
   * Exports a buildfile using file names.
   * @param ref Reference FASTA filename
   * @param vcf VCF filename
   * @param build filename of output file, default std::cout
   */
  void exportBuildfile(std::string ref, std::string vcf, std::string build = "");

  /**
   * Export a buildfile to std::cout.
   * @param reference stream to intepret as the referense sequence
   * @param variants a vcfstream to obtain variant records from
   */
  void exportBuildfile(std::istream *reference, vcfstream &variants) {
    exportBuildfile(reference, variants, std::cout);
  }

  /**
   * Export a buildfile to a output stream.
   * @param _reference input stream to use as input sequence.
   * @param variants input vcfstream to obtain variant records
   * @param buildout output stream for buildfile.
   */
  void exportBuildfile(std::istream *_reference, vcfstream &variants, std::ostream &buildout);

  /**
   * Build a graph using a buildfile.
   * @build filename of a graph buildfile.
   */
  void buildGraph(std::string build) {
    std::ifstream b(build);
    if (!b.good()) throw std::invalid_argument("Error opening file.");
    buildGraph(b);
    b.close();
  }

  /**
   * Build a graph using a buildfile input stream.
   * @param buildfile buildfile input stream.
   */
  void buildGraph(std::istream &buildfile);

  /**
   * Align a read to the graph.
   * @param read Read object to align.
   * @return an Alignment object
   */
  Alignment *align(const Read &read);

  /**
   * Align a read to an existing alignment object.
   * @param r Read to align
   * @param a Alignment to store the alignment result in.
   */
  void align(const Read &r, Alignment &a);

  /*********************************** SETTERS & GETTERS ***********************************/

  /**
   * Get a pointer to the internal gssw_graph object.
   * @returns gssw_graph pointer
   */
  gssw_graph *getGSSWGraph() const { return graph; }

  /**
   * Get a copy of the graph parameters.
   * @returns GraphParams copy
   */
  GraphParams getParamsCopy() const { return params; }

  /**
   * Set the parameters of the graph.
   * @param p GraphParams to set.
   */
  void setParams(GraphParams p) {
    if (p.nt_table == NULL) {
      p.nt_table = gssw_create_nt_table();
    }
    if (p.mat != NULL) {
      free(p.mat);
    }
    p.mat = gssw_create_score_matrix(p.match, p.mismatch);
    params = p;
  }

  /**
   * Set the ingroup percentage.
   * @param i percent of individuals to include
   */
  void setIngroup(int32_t i) { params.ingroup = i; }

  /**
   * Select if a complement graph should be built.
   * @param b if true generate a complement graph.
   */
  void setComplement(bool b) { params.genComplement = b; }

  /**
   * Set the source of a complenment graph
   * @param s filename of the complement origin buildfile.
   */
  void setComplementSource(std::string s) { params.complementSource = s; }

  /**
   * Set the scores for alignment.
   * @param m match score
   * @param mm mismatch score
   * @param go gap_open score
   * @param ge gap_extend score
   */
  void setScores(int32_t m = 2, int32_t mm = 2, uint8_t go = 3, uint8_t ge = 1) {
    if (graph) std::cerr << "Graph must be rebuilt for new scores to take effect." << std::endl;
    params.match = m;
    params.mismatch = mm;
    params.gap_extension = ge;
    params.gap_open = go;

    if (params.mat != NULL) {
      free(params.mat);
    }
    params.mat = gssw_create_score_matrix(params.match, params.mismatch);
  }

  /**
   * Keep only the allele with the max allele frequency.
   * @param b True to only keep maximum allele frequency.
   */
  void setMaxAF(bool b) { params.maxAF = b; }

  /**
   * If true, load the individuals into memory.
   * @param b
   */
  void useIndividuals(bool b) { params.includeIndividuals = b; }

 protected:
  GraphParams params; // Graph construction parameters
  gssw_graph *graph = NULL; // The raw graph

  /**
   * Parse a string into an upper and lower bound.
   * @param region String of format "a:b"
   * @param min parsed "a"
   * @param max parsed "b"
   */
  void parseRegion(std::string region, uint32_t *min, uint32_t *max);

  /**
   * Create an ingroup as specified by local params.
   * @param variants vcfstream to set the ingroup for.
   */
  void generateIngroup(vcfstream &variants);
};

inline std::ostream &operator<<(std::ostream &os, const Graph::GraphParams &gp) {
  std::stringstream ss;
  ss << "maxNodeLen:" << gp.maxNodeLen;
  if (gp.region.length() > 0) ss << ", Region:" << gp.region;
  ss << ",maxAF:" << gp.maxAF << ",refIndivs:" << gp.includeRefIndivs;
  if (gp.genComplement) ss << ",ComplementSource:" << gp.complementSource;
  else ss << ",ingroup:" << gp.ingroup;
  os << ss.str();
  return os;
}

}

#endif //VARGAS_GRAPH_H
