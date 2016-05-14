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
#include "alignment.h"
#include <memory>
#include <set>

namespace vargas {

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned long ulong;

/**
 * Converts a character to a numeral representation.
 * @param c character
 * @return numeral repersentation
 */
inline uchar baseToNum(char c) {
  switch (c) {
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
    default:
      return 4;
  }
}


/**
 * Convert a numeric form to a char, upper case.
 * All ambigious bases are represented as 'N'
 * @param num numeric form
 * @return char in [A,C,G,T,N]
 */
inline char numToBase(uchar num) {
  switch (num) {
    case 0:
      return 'A';
    case 1:
      return 'C';
    case 2:
      return 'G';
    case 3:
      return 'T';
    default:
      return 'N';
  }
}


/**
 * Convert a character sequence in a numeric sequence
 * @param seq Sequence string
 * @return vector of numerals
 */
inline std::vector<uchar> seqToNum(const std::string &seq) {
  std::vector<uchar> num(seq.length());
  std::transform(seq.begin(), seq.end(), num.begin(), baseToNum);
  return num;
}


/**
 * Convert a numeric vector to a sequence of bases.
 * @param num Numeric vector
 * @return sequence string, Sigma={A,G,T,C,N}
 */
inline std::string numToSeq(const std::vector<uchar> &num) {
  std::stringstream builder;
  for (auto &n : num) {
    builder << numToBase(n);
  }
  return builder.str();
}


/**
 * Represents a graph of the genome. The graph is backed by a map of Graph::Nodes, and edges
 * are backed by a map of node ID's.
 */
class graph {

 public:

  /**
   * Represents a node in the directed graphs. Edges stored externally.
   */
  class Node {
   public:
    // Assign a unique ID to each node
    Node() : _id(newID++) { }

    // Access functions
    ulong length() const { return _seq.size(); } // Length of sequence
    ulong end() const { return _endPos; } // Sequence end position in genome
    bool belongs(uint ind) const { return _individuals[ind]; } // Check if a certain individual has this node
    const std::vector<uchar> &seq() const { return _seq; } // Sequence in numeric form
    ulong popSize() const { return _individuals.size(); } // How many individuals are represented in the node
    long id() const { return _id; } // Node ID
    bool isRef() const { return _ref; } // True if part of the reference seq


    // Setup functions
    void setID(long id) { this->_id = id; }
    void setEndPos(ulong pos) { this->_endPos = pos; }
    void setPopulation(std::vector<bool> &pop) { _individuals = pop; }
    void setSeq(std::string seq) { _seq = seqToNum(seq); }
    void setSeq(std::vector<uchar> &seq) { this->_seq = seq; }
    void setAsRef() { _ref = true; }
    void setAsNotRef() { _ref = false; }

   protected:
    static long newID; // ID of the next instance to be created

   private:
    long _id;
    ulong _endPos; // End position of the sequence
    std::vector<bool> _individuals; // Each bit marks an individual, 1 if they have this node
    std::vector<uchar> _seq;
    bool _ref = false; // Part of the reference sequence
  };

  typedef std::shared_ptr<Node> nodeptr;

  /**
   * Default constructor inits a new graph, including a new node map.
   */
  graph() : _IDMap(std::make_shared<std::map<long, nodeptr>>(std::map<long, nodeptr>())) { }

  /**
   * Create a graph with another graph and a population filter. The new graph will only
   * contain nodes if any of the individuals in filter possess the node. The actual nodes
   * are shared_ptr's to the parent graph, as to prevent duplication of Nodes.
   * @param g Graph to derive the new graph from
   * @param filter population filter, only include nodes representative of this population
   */
  graph(const graph &g, const std::vector<bool> &filter) {
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
    std::map<long, nodeptr> includedNodes;
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
      if (g._nextMap.find(n.second->id()) == g._nextMap.end()) continue;
      for (auto e : g._nextMap.at(n.second->id())) {
        if (includedNodes.find(e) != includedNodes.end()) {
          addEdge(n.second->id(), e);
        }
      }
    }

    // Set the new root
    if (includedNodes.find(g.root()) == includedNodes.end()) {
      throw std::invalid_argument("Currently the root must be common to all graphs.");
    }
    _root = g.root();
    finalize();
  }

  void finalize() {
    _toposort.clear();
    std::set<long> unmarked, tempmarked, permmarked;
    for (auto &n : *_IDMap) {
      unmarked.insert(n.first);
    }
    while (!unmarked.empty()) {
      visit(*unmarked.begin(), unmarked, tempmarked, permmarked);
    }
    std::reverse(_toposort.begin(), _toposort.end());

  }

  void visit(long n, std::set<long> &unmarked, std::set<long> &temp, std::set<long> &perm) {
    if (temp.count(n) != 0) throw std::out_of_range("Graph is not acyclic.");
    if (unmarked.count(n)) {
      unmarked.erase(n);
      temp.insert(n);
      for (auto m : _nextMap[n]) {
        visit(m, unmarked, temp, perm);
      }
      temp.erase(n);
      perm.insert(n);
      _toposort.push_back(n);
    }
  }

  /**
   * Add a new node to the graph. A new node is created so the original can be destroyed.
   * The first node added is set as the graph root.
   */
  long addNode(const Node &n) {
    if (_popSize < 0) _popSize = n.popSize(); // first node dictates graph population size
    else if (_popSize != n.popSize()) return 0; // graph should be for same population size
    if (_IDMap->find(n.id()) != _IDMap->end()) return 0; // make sure node isn't duplicate
    if (_root < 0) _root = n.id(); // first node added is default root

    _IDMap->emplace(n.id(), std::make_shared<Node>(n));
    return n.id();
  }

  /**
   * Create an edge linking two nodes. Previous and Next edges are added.
   * @param n1 Node one ID
   * @param n2 Node two ID
   */
  bool addEdge(long n1, long n2) {
    // Check if the nodes exist
    if (_IDMap->find(n1) == _IDMap->end() || _IDMap->find(n2) == _IDMap->end()) return false;

    // init if first edge to be added
    if (_nextMap.find(n1) == _nextMap.end()) {
      _nextMap[n1] = std::vector<long>();
    }
    if (_prevMap.find(n2) == _prevMap.end()) {
      _prevMap[n2] = std::vector<long>();
    }
    _nextMap[n1].push_back(n2);
    _prevMap[n2].push_back(n1);
    _toposort.clear(); // any ordering is invalidated
    return true;
  }

  /**
   * Sets the root of the graph.
   * @param id ID of root node
   */
  void setRoot(long id) {
    _root = id;
    _toposort.clear(); // any ordering is invalidated
  }

  // Return root node ID
  long root() const { return _root; }

  // const reference to node map
  const std::shared_ptr<std::map<long, nodeptr>> &getNodeMap() const { return _IDMap; }

  const Node &getNode(long id) const { return *(*_IDMap)[id]; }

  // Export the graph in DOT format.
  std::string toDOT(std::string name = "graph") const {
    std::stringstream dot;
    dot << "digraph " << name << " {\n";
    for (auto &n : _nextMap) {
      for (auto e : n.second) {
        dot << n.first << " -> " << e << ";\n";
      }
    }
    dot << "}\n";
    return dot.str();
  }


  class graphIter {

   public:
    graphIter(const graph &g) : _graph(g), _idx(0) { }
    graphIter(const graph &g, long index) : _graph(g), _idx(index) { }
    ~graphIter() { }

    graphIter &operator=(const graphIter &other) {
      _idx = other._idx;
      return *this;
    }

    bool operator==(const graphIter &other) const {
      if (other._graph._toposort != _graph._toposort) return false;
      return _idx == other._idx;
    }

    bool operator!=(const graphIter &other) const {
      return _idx != other._idx;
    }

    graphIter &operator++() {
      if (_idx < _graph._toposort.size()) {
        _idx++;
      }
      return *this;
    }

    const graph::Node &operator*() const { return _graph.getNode(_graph._toposort[_idx]); }

   private:
    const graph &_graph;
    long _idx;
  };

  graphIter begin() const {
    if (_toposort.size() == 0 && _IDMap->size() > 0) {
      throw std::logic_error("Graph must be finalized before iteration.");
    }
    return graphIter(*this);
  }

  graphIter end() const {
    return graphIter(*this, _toposort.size());
  }


 private:
  long _root = -1; // Root of the graph
  // maps a node ID to a nodeptr. Any derived graphs use the same base node ID map.
  std::shared_ptr<std::map<long, nodeptr>> _IDMap;
  // maps a node ID to the vector of nodes it points to
  std::map<long, std::vector<long>> _nextMap;
  // maps a node ID to a vector of node ID's that point to it
  std::map<long, std::vector<long>> _prevMap;
  long _popSize = -1; // Used to make sure we have the same population size for all nodes in the graph
  std::vector<long> _toposort;

};





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
    int8_t *scoreMat = gssw_create_score_matrix(match, mismatch);
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
    if (!params.scoreMat) params.scoreMat = gssw_create_score_matrix(p.match, p.mismatch);
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
    ref = new std::ifstream(refFile);

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
    if (graph != nullptr) gssw_graph_destroy(graph);
    free(params.nt_table);
    free(params.scoreMat);
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
  GraphParams getParams() const { return params; }

  /**
   * Set the parameters of the graph.
   * @param p GraphParams to set.
   */
  void setParams(GraphParams p) {
    if (p.nt_table == nullptr) {
      p.nt_table = gssw_create_nt_table();
    }
    if (p.scoreMat != nullptr) {
      free(p.scoreMat);
    }
    p.scoreMat = gssw_create_score_matrix(p.match, p.mismatch);
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
   * Set the scores for alignment. If set after the graph is made, it must be rebuilt.
   * @param m match score
   * @param mm mismatch score
   * @param go gap_open score
   * @param ge gap_extend score
   */
  void setScores(int32_t m = 2, int32_t mm = 2, uint8_t go = 3, uint8_t ge = 1) {
    params.match = m;
    params.mismatch = mm;
    params.gap_extension = ge;
    params.gap_open = go;

    if (params.scoreMat != nullptr) {
      free(params.scoreMat);
    }
    params.scoreMat = gssw_create_score_matrix(params.match, params.mismatch);
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
  gssw_graph *graph = nullptr; // The raw graph

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
