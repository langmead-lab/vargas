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


#include "../doctest/doctest/doctest.h"
#include "../gssw/src/gssw.h"
#include "readsource.h"
#include "vcfstream.h"
#include "alignment.h"
#include "utils.h"
#include "../htslib/htslib/faidx.h"
#include "../htslib/htslib/vcfutils.h"
#include <cstdio>
#include <memory>
#include <set>
#include <htslib/vcf.h>

namespace vargas {

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned long ulong;


/**
 * Converts a character to a numeral representation.
 * @param c character
 * @return numeral repersentation
 */
inline uchar base_to_num(char c) {
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
inline char num_to_base(uchar num) {
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
inline std::vector<uchar> seq_to_num(const std::string &seq) {
  std::vector<uchar> num(seq.length());
  std::transform(seq.begin(), seq.end(), num.begin(), base_to_num);
  return num;
}
TEST_CASE ("Sequence to Numeric") {
  std::vector<uchar> a = seq_to_num("ACGTN");
      REQUIRE(a.size() == 5);
      CHECK(a[0] == 0);
      CHECK(a[1] == 1);
      CHECK(a[2] == 2);
      CHECK(a[3] == 3);
      CHECK(a[4] == 4);
}

/**
 * Convert a numeric vector to a sequence of bases.
 * @param num Numeric vector
 * @return sequence string, Sigma={A,G,T,C,N}
 */
inline std::string num_to_seq(const std::vector<uchar> &num) {
  std::stringstream builder;
  for (auto &n : num) {
    builder << num_to_base(n);
  }
  return builder.str();
}
TEST_CASE ("Numeric to Sequence") {
  std::string a = num_to_seq({0, 1, 2, 3, 4});
      REQUIRE(a.length() == 5);
      CHECK(a[0] == 'A');
      CHECK(a[1] == 'C');
      CHECK(a[2] == 'G');
      CHECK(a[3] == 'T');
      CHECK(a[4] == 'N');
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
    Node() : _id(_newID++) { }

    // Access functions
    ulong length() const { return _seq.size(); } // Length of sequence
    ulong end() const { return _endPos; } // Sequence end position in genome
    bool belongs(uint ind) const { return _individuals[ind]; } // Check if a certain individual has this node
    const std::vector<uchar> &seq() const { return _seq; } // Sequence in numeric form
    ulong popSize() const { return _individuals.size(); } // How many individuals are represented in the node
    long id() const { return _id; } // Node ID
    bool isRef() const { return _ref; } // True if part of the reference seq


    static long _newID; // ID of the next instance to be created

    // Setup functions
    void setID(long id) {
      if (id >= _newID) {
        this->_id = id;
        _newID = ++id;
      }
    }
    void setEndPos(ulong pos) { this->_endPos = pos; }
    void setPopulation(std::vector<bool> &pop) { _individuals = pop; }
    void setSeq(std::string seq) { _seq = seq_to_num(seq); }
    void setSeq(std::vector<uchar> &seq) { this->_seq = seq; }
    void setAsRef() { _ref = true; }
    void setAsNotRef() { _ref = false; }

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

  /**
   * Builds the topographical sort of the graph, used for graph iteration.
   */
  void finalize() {
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

  /**
   * Add a new node to the graph. A new node is created so the original can be destroyed.
   * The first node added is set as the graph root.
   */
  long add_node(const Node &n) {
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
  bool add_edge(long n1, long n2) {
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

  /**
   * Sets the root of the graph.
   * @param id ID of root node
   */
  void set_root(long id) {
    _root = id;
  }

  void set_desc(std::string description) { _desc = description; }

  // Return root node ID
  long root() const { return _root; }

  // const reference to maps
  const std::shared_ptr<std::map<long, nodeptr>> &node_map() const { return _IDMap; }
  const std::map<long, std::vector<long>> &next_map() const { return _next_map; }
  const std::map<long, std::vector<long>> &prev_map() const { return _prev_map; }
  const Node &node(long id) const { return *(*_IDMap)[id]; }

  ulong pop_size() const {
    if (_root >= 0) return _IDMap->at(0)->popSize();
    return 0;
  }

  std::string desc() const { return _desc; }

  // Export the graph in DOT format.
  std::string to_DOT(std::string name = "graph") const {
    std::stringstream dot;
    dot << "digraph " << name << " {\n";
    for (auto &n : _next_map) {
      for (auto e : n.second) {
        dot << n.first << " -> " << e << ";\n";
      }
    }
    dot << "}\n";
    return dot.str();
  }

  /**
   * const forward iterator to traverse the graph topologically.
   */
  class GraphIter {

   public:
    GraphIter(const graph &g) : _graph(g), _idx(0) { }
    GraphIter(const graph &g, long index) : _graph(g), _idx(index) { }
    ~GraphIter() { }

    GraphIter &operator=(const GraphIter &other) {
      _idx = other._idx;
      return *this;
    }

    bool operator==(const GraphIter &other) const {
      // Check if comparing like-graphs (weak check)
      if (other._graph._toposort != _graph._toposort) return false;
      return _idx == other._idx;
    }

    bool operator!=(const GraphIter &other) const {
      return _idx != other._idx;
    }

    GraphIter &operator++() {
      if (_idx < _graph._toposort.size()) {
        _idx++;
      }
      return *this;
    }

    GraphIter &operator--() {
      if (_idx > 0) {
        _idx--;
      }
      return *this;
    }

    const graph::Node &operator*() const { return _graph.node(_graph._toposort[_idx]); }

   private:
    const graph &_graph;
    long _idx;
  };

  /**
   * Provides an iterator to a topological sorting of the graph.
   */
  GraphIter begin() const {
    if (_toposort.size() == 0 && _IDMap->size() > 0) {
      throw std::logic_error("Graph must be finalized before iteration.");
    }
    return GraphIter(*this);
  }

  /**
   * Iterator to end of topological sorting.
   */
  GraphIter end() const {
    return GraphIter(*this, _toposort.size());
  }

 private:
  long _root = -1; // Root of the graph
  // maps a node ID to a nodeptr. Any derived graphs use the same base node ID map.
  std::shared_ptr<std::map<long, nodeptr>> _IDMap;
  // maps a node ID to the vector of nodes it points to
  std::map<long, std::vector<long>> _next_map;
  // maps a node ID to a vector of node ID's that point to it
  std::map<long, std::vector<long>> _prev_map;
  long _popSize = -1; // Used to make sure we have the same population size for all nodes in the graph
  std::vector<long> _toposort; // Sorted graph
  // Description, used by the builder to store construction params
  std::string _desc;

  /**
   * Recursive depth first search to find dependencies. Used to topological sort.
   * @param n current node ID
   * @param unmarked set of unvisited nodes
   * @param temp set of visited but unadded nodes
   * @param perm set of completed nodes
   */
  void _visit(long n, std::set<long> &unmarked, std::set<long> &temp, std::set<long> &perm) {
    if (temp.count(n) != 0) throw std::domain_error("Graph contains a cycle.");
    if (unmarked.count(n)) {
      unmarked.erase(n);
      temp.insert(n);
      if (_next_map.find(n) != _next_map.end()) {
        for (auto m : _next_map[n]) {
          _visit(m, unmarked, temp, perm);
        }
      }
      temp.erase(n);
      perm.insert(n);
      _toposort.push_back(n);
    }
  }

};

TEST_CASE ("Node tests") {
  vargas::graph::Node::_newID = 0;
  vargas::graph::Node n1;
  vargas::graph::Node n2;
      CHECK(n1.id() == 0);
      CHECK(n2.id() == 1);

      SUBCASE("Node ID change") {
    n1.setID(1);
        CHECK(n1.id() == 0);
    n1.setID(2);
        CHECK(n1.id() == 2);
  }

      SUBCASE("Set Node params") {
    n1.setSeq("ACGTN");
    std::vector<bool> a = {0, 0, 1};
    n1.setPopulation(a);
    n1.setAsRef();
    n1.setEndPos(100);

        REQUIRE(n1.seq().size() == 5);

        CHECK(n1.seq()[0] == 0);
        CHECK(n1.seq()[1] == 1);
        CHECK(n1.seq()[2] == 2);
        CHECK(n1.seq()[3] == 3);
        CHECK(n1.seq()[4] == 4);
        CHECK(n1.isRef());
        CHECK(n1.end() == 100);
        CHECK(!n1.belongs(0));
        CHECK(!n1.belongs(1));
        CHECK(n1.belongs(2));
  }

}

TEST_CASE ("graph class") {
  vargas::graph::Node::_newID = 0;
  vargas::graph g;

  /**   GGG
  *    /   \
  * AAA     TTT
  *    \   /
  *     CCC(ref)
  */

  {
    vargas::graph::Node n;
    n.setEndPos(3);
    n.setAsRef();
    std::vector<bool> a = {0, 1, 1};
    n.setPopulation(a);
    n.setSeq("AAA");
    g.add_node(n);
  }

  {
    vargas::graph::Node n;
    n.setEndPos(6);
    n.setAsRef();
    std::vector<bool> a = {0, 0, 1};
    n.setPopulation(a);
    n.setSeq("CCC");
    g.add_node(n);
  }

  {
    vargas::graph::Node n;
    n.setEndPos(6);
    n.setAsNotRef();
    std::vector<bool> a = {0, 1, 0};
    n.setPopulation(a);
    n.setSeq("GGG");
    g.add_node(n);
  }

  {
    vargas::graph::Node n;
    n.setEndPos(9);
    n.setAsRef();
    std::vector<bool> a = {0, 1, 1};
    n.setPopulation(a);
    n.setSeq("TTT");
    g.add_node(n);
  }

  g.add_edge(0, 1);
  g.add_edge(0, 2);
  g.add_edge(1, 3);
  g.add_edge(2, 3);

      CHECK_THROWS(g.begin());

  g.finalize();

      REQUIRE(g.node_map()->size() == 4);
      REQUIRE(g.prev_map().size() == 3);
      REQUIRE(g.next_map().size() == 3);

  // Check forward edges
      REQUIRE(g.next_map().at(0).size() == 2);
      REQUIRE(g.next_map().at(1).size() == 1);
      REQUIRE(g.next_map().at(2).size() == 1);
      REQUIRE(g.next_map().count(3) == 0);

  // Check prev edges
      REQUIRE(g.prev_map().count(0) == 0);
      REQUIRE(g.prev_map().at(1).size() == 1);
      REQUIRE(g.prev_map().at(2).size() == 1);
      REQUIRE(g.prev_map().at(3).size() == 2);

      SUBCASE("Proper graph setup") {
        CHECK(num_to_seq(g.node(0).seq()) == "AAA");
        CHECK(num_to_seq(g.node(1).seq()) == "CCC");
        CHECK(num_to_seq(g.node(2).seq()) == "GGG");
        CHECK(num_to_seq(g.node(3).seq()) == "TTT");
  }

      SUBCASE("Topographical invalidation") {
    g.add_edge(1, 2);
        CHECK_THROWS(g.begin());
    g.finalize();
        CHECK_NOTHROW(g.begin());
  }

      SUBCASE("Cyclic graph") {
    g.add_edge(3, 1);
        CHECK_THROWS(g.finalize());
  }

      SUBCASE("graph iterator") {
    // Node visit order should be topological
    vargas::graph::GraphIter i = g.begin();

        CHECK(num_to_seq((*i).seq()) == "AAA");
    ++i;

    // Order of these two don't matter
    bool mid = (num_to_seq((*i).seq()) == "CCC") || (num_to_seq((*i).seq()) == "GGG");
        CHECK(mid);
    ++i;
    mid = (num_to_seq((*i).seq()) == "CCC") || (num_to_seq((*i).seq()) == "GGG");
        CHECK(mid);
    ++i;

        CHECK(num_to_seq((*i).seq()) == "TTT");
    ++i;
        CHECK(i == g.end());
    ++i;
        CHECK(i == g.end());
  }

      SUBCASE("Derived graph") {
    std::vector<bool> filter = {0, 0, 1};
    vargas::graph g2(g, filter);

        CHECK(g2.node_map()->size() == 4); // Underlying node map unchanged
        CHECK(g2.next_map().size() == 2);
        CHECK(g2.prev_map().size() == 2);

        CHECK(g2.next_map().at(0).size() == 1);
        CHECK(g2.next_map().at(1).size() == 1);
        CHECK(g2.next_map().count(2) == 0); // This node shouldn't be included
        CHECK(g2.next_map().count(3) == 0);
        CHECK(g2.prev_map().count(0) == 0);
        CHECK(g2.prev_map().at(1).size() == 1);
        CHECK(g2.prev_map().at(3).size() == 1);
  }

}


/**
 * Provides an interface for a FASTA File. An index is built if one does not
 * already exist.
 */
class FASTAFile {

 public:
  /**
   * Load a given FASTA and index, create if none exists.
   * @param file filename
   */
  FASTAFile(std::string file) : _file_name(file) { _init(); }

  FASTAFile() { }

  ~FASTAFile() {
    if (_index) fai_destroy(_index);
  }

  /**
   * Open a specified FASTA file and make an index.
   * @param file filename
   * @return -1 on index build error, -2 on open error, 0 otherwise
   */
  int open(std::string file) {
    _file_name = file;
    return _init();
  }

  /**
   * Check the number of sequences in the index.
   * @return number of sequences in the FASTA file
   */
  int num_seq() const { return faidx_nseq(_index); }

  /**
   * Return a subsequence of the FASTA file, absolute 0 based indexing.
   * The position is not relative to the min/max params.
   * @param beg beginning index
   * @param end ending index, inclusive
   * @return subsequence string
   */
  std::string subseq(std::string chr, int beg, int end) const {
    int len;
    char *ss = faidx_fetch_seq(_index, chr.c_str(), beg, end, &len);
    std::string ret(ss);
    free(ss);
    return ret;
  }

  /**
   * Sequence name given a sequence ID.
   * @param ID of sequence
   * @return sequence name
   */
  std::string seq_name(int i) const {
    if (i > num_seq()) return "";
    return std::string(faidx_iseq(_index, i));
  }

  /**
   * Get all sequences in the File
   * @return vector of sequence names.
   */
  std::vector<std::string> sequences() const {
    std::vector<std::string> ret;
    for (int i = 0; i < num_seq(); ++i) {
      ret.push_back(seq_name(i));
    }
    return ret;
  }

  bool good() const { return _index != nullptr; }

  std::string file() const { return _file_name; }


 protected:

  /**
   * Loads a FASTA index. If the index does not exist, one is created.
   */
  int _init() {
    // Check if a Fasta index exists. If it doesn't build it.
    bool exists;
    {
      std::ifstream fidx(_file_name + ".fai");
      exists = fidx.good();
    }

    if (!exists) {
      if (fai_build(_file_name.c_str()) < 0) {
        _file_name = "";
        return -1;
      }
    }

    _index = fai_load(_file_name.c_str());
    if (!_index) return -2;
    return 0;
  }

 private:
  std::string _file_name;
  faidx_t *_index = nullptr;
};
TEST_CASE ("FASTA Handler") {
  using std::endl;
  std::string tmpfa = "tmp_tc.fa";
  {
    std::ofstream fao(tmpfa);
    fao
        << ">x" << endl
        << "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGTTCCTGGTGCTATGTGTAACTAGTAATGG" << endl
        << "TAATGGATATGTTGGGCTTTTTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAA" << endl
        << "TGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTTCAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGA" << endl
        << "ACAAGTTAGTTAATCTCTCTGAACTTCAGTTTAATTATCTCTAATATGGAGATGATACTACTGACAGCAGAGGTTTGCTG" << endl
        << "TGAAGATTAAATTAGGTGATGCTTGTAAAGCTCAGGGAATAGTGCCTGGCATAGAGGAAAGCCTCTGACAACTGGTAGTT" << endl
        << "ACTGTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCCTCCTTGAGGCTCT" << endl
        << "TCTGGCTTTTCATTGTCAACACAGTCAACGCTCAATACAAGGGACATTAGGATTGGCAGTAGCTCAGAGATCTCTCTGCT" << endl
        << ">y" << endl
        << "GGAGCCAGACAAATCTGGGTTCAAATCCTGGAGCCAGACAAATCTGGGTTCAAATCCTGGAGCCAGACAAATCTGGGTTC" << endl;
  }
  FASTAFile fa(tmpfa);

      CHECK(fa.num_seq() == 2);
      REQUIRE(fa.sequences().size() == 2);
      CHECK(fa.seq_name(0) == "x");
      CHECK(fa.seq_name(1) == "y");
      CHECK(fa.subseq("x", 0, 3) == "CAAA");
      CHECK(fa.subseq("y", 0, 2) == "GGA");
      CHECK(fa.sequences()[0] == "x");
      CHECK(fa.sequences()[1] == "y");

  remove(tmpfa.c_str());
  remove((tmpfa + ".fai").c_str());
}

class VarFile {
 public:

  /**
  * Get the specified format field from the record.
  * @param T Valid types are int, char, float.
  */
  template<typename T>
  class FormatField {
   public:
    /**
     * Get the specified field.
     * @param hdr VCF Header
     * @param rec current record
     * @param tag Field to get, e.g. "GT"
     */
    FormatField(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag) : tag(tag) {
      if (!hdr || !rec || tag.length() == 0) throw std::invalid_argument("Invalid header, rec, or tag.");

      T *dst = NULL;
      int n_arr = 0;

      int n = _get_vals(hdr, rec, tag, &dst, n_arr);
      if (n == -1) throw std::invalid_argument("No such tag in header: " + tag);
      else if (n == -2) throw std::invalid_argument("Header and tag type clash: " + tag);
      else if (n == -3) throw std::invalid_argument(tag + " does not exist in record.");

      for (int i = 0; i < n; ++i) {
        values.push_back(dst[i]);
      }

      free(dst); // get_format_values allocates
    }

    std::vector<T> values;
    std::string tag; // Type of FORMAT or INFO field.

   private:
    // Change the parse type based on what kind of type we have
    int _get_vals(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, int32_t **dst, int &ndst) {
      return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_INT);
    }
    int _get_vals(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, float **dst, int &ndst) {
      return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_REAL);
    }
    int _get_vals(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, char **dst, int &ndst) {
      return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_STR);
    }
  };

  /**
   * Get the specified info field from the record.
   * @param T Valid types are int, char, float.
   */
  template<typename T>
  class InfoField {
   public:
    /**
     * Get the specified field.
     * @param hdr VCF Header
     * @param rec current record
     * @param tag Field to get, e.g. "GT"
     */
    InfoField(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag) : tag(tag) {
      if (!hdr || !rec || tag.length() == 0) throw std::invalid_argument("Invalid header, rec, or tag.");

      T *dst = NULL;
      int n_arr = 0;

      int n = _bcf_get_info_values(hdr, rec, tag, &dst, n_arr);
      if (n == -1) throw std::invalid_argument("No such tag in header: " + tag);
      else if (n == -2) throw std::invalid_argument("Header and tag type clash: " + tag);
      else if (n == -3) throw std::invalid_argument(tag + " does not exist in record.");

      for (int i = 0; i < n; ++i) {
        values.push_back(dst[i]);
      }

      free(dst);
    }

    std::vector<T> values;
    std::string tag; // Type of FORMAT or INFO field.

   private:
    // Change the parse type based on what kind of type we have
    int _bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, int32_t **dst, int &ndst) {
      return bcf_get_info_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_INT);
    }
    int _bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, float **dst, int &ndst) {
      return bcf_get_info_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_REAL);
    }
    int _bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, char **dst, int &ndst) {
      return bcf_get_info_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_STR);
    }
  };

  VarFile() { }
  VarFile(std::string file) : _file_name(file) { _init(); }
  VarFile(std::string file, std::string chr, int min, int max) :
      _file_name(file), _chr(chr), _min_pos(min), _max_pos(max) { _init(); }
  ~VarFile() {
    if (_bcf) hts_close(_bcf);
    if (_header) bcf_hdr_destroy(_header);
    bcf_destroy(_curr_rec);
  }

  /**
   * Open the specified VCF or BCF file. The first record is loaded.
   * Do not use this to load a new file after another one, rather create a new VarFile object.
   * @param file filename
   * @return -1 on file open error, -2 on header load error, 0 otherwise
   */
  int open(std::string file) {
    _file_name = file;
    return _init();
  }
  /**
 * Set the minimum and maximum position
 * @param min Minimum position, 0 indexed
 * @param max Max position, inclusive, 0 indexed
 */
  void set_region(std::string chr, int min, int max) {
    _min_pos = min;
    _max_pos = max;
    _chr = chr;
  }

  /**
   * Parse a region string in the format
   * CHR:XX,XXX-YY,YYY
   * commas are stripped, range is inclusive. 0 indexed.
   * @param region region string
   */
  void set_region(std::string region) {
    std::vector<std::string> regionSplit = split(region, ':');

    // Name
    if (regionSplit.size() != 2) throw std::invalid_argument("Invalid region format, should be CHR:XX,XXX-YY,YYY");
    _chr = regionSplit[0];

    // Strip commas
    regionSplit.erase(std::remove(regionSplit.begin(), regionSplit.end(), ","), regionSplit.end());

    // Range
    regionSplit = split(regionSplit[1], '-');
    if (regionSplit.size() != 2) throw std::invalid_argument("Invalid region format, should be CHR:XX,XXX-YY,YYY");
    _min_pos = std::stoi(regionSplit[0]);
    _max_pos = std::stoi(regionSplit[1]);

  }

  /**
   * Get a list of sequences in the VCF file.
   * @return vector of sequence names
   */
  std::vector<std::string> sequences() const {
    std::vector<std::string> ret;
    int num;
    const char **seqnames = bcf_hdr_seqnames(_header, &num);
    if (!seqnames) {
      throw std::invalid_argument("Error reading VCF header.");
    }
    for (int i = 0; i < num; ++i) {
      ret.push_back(std::string(seqnames[i]));
    }
    free(seqnames);
    return ret;
  }

  /**
   * @return Number of samples the VCF has. Each sample represents two genotypes.
   */
  int num_samples() const { return bcf_hdr_nsamples(_header); }

  /**
   * Load the next VCF record. Only shared info is unpacked.
   * @return false on read error or if outside restriction range.
   */
  bool next() {
    int status = bcf_read(_bcf, _header, _curr_rec);
    if (_max_pos > 0 && _curr_rec->pos > _max_pos) return false;
    unpack_shr();
    // Skip it if its the wrong CHROM.
    if (_chr.length() != 0 && _chr != std::string(bcf_hdr_id2name(_header, _curr_rec->rid))) {
      next();
    }
    return status >= 0;
  }

  /**
   * Unpack only the shared information, and loads ref and allele info.
   */
  void unpack_shr() {
    bcf_unpack(_curr_rec, BCF_UN_SHR);
    _load_shared();
  }

  /**
   * Unpacks shared information as well as all sample information.
   * Subject to sample set restrictions.
   */
  void unpack_all() {
    bcf_unpack(_curr_rec, BCF_UN_ALL);
  }

  std::string ref() const { return _alleles[0]; }
  const std::vector<std::string> alleles() const { return _alleles; }
  /**
   * 0 based position, i.e. the VCF pos - 1.
   * @return position.
   */
  int pos() const { return _curr_rec->pos; }

  // TODO Better handling policy of unknown elements?
  /**
   * Get a list of alleles for all samples (subject to sample set restriction).
   * Consecutive alleles represent phasing, e.g. all odd indexes are one phase,
   * all even indexes are the other. Call will unpack the full record.
   * Explicit copy number variations are replaced, other ambigous types are replaced
   * ambiguous.
   * @return Vector of alleles, ordered by sample.
   */
  const std::vector<std::string> &genotypes() {
    FormatField<int> gt(_header, _curr_rec, "GT");
    _genotypes.clear();
    for (int o : gt.values) {
      _genotypes.push_back(_alleles[bcf_gt_allele(o)]);
    }
    return _genotypes;
  }


  bool good() const { return !_header && !_bcf; }

  std::string file() const { return _file_name; }

 protected:
  int _init() {
    _bcf = bcf_open(_file_name.c_str(), "r");
    if (_bcf == NULL) return -1;
    _header = bcf_hdr_read(_bcf);
    if (!_header) return -2;

    next();
    // Skip records of wrong CHROM and POS as imposed by filter
    while (_curr_rec->pos < _min_pos) {
      next();
    }
    unpack_shr();
    return 0;
  }

  void _load_shared() {
    _alleles.clear();
    for (int i = 0; i < _curr_rec->n_allele; ++i) {
      std::string allele(_curr_rec->d.allele[i]);
      // Some replacement tag
      if (allele.at(0) == '<') {
        std::string ref = _curr_rec->d.allele[0];
        // Copy number
        if (allele.substr(1, 2) == "CN" && allele.at(3) != 'V') {
          int copy = std::stoi(allele.substr(3, allele.length() - 4));
          allele = "";
          for (int i = 0; i < copy; ++i) allele += ref;
        } else {
          // Other types are just subbed with the ref.
          allele = ref;
        }
      }
      _alleles.push_back(allele);
    }
  }


 private:
  std::string _file_name, // VCF file name
      _chr; // Sequence restriction
  int _min_pos = -1, _max_pos = -1; // Region of _chr restriction

  htsFile *_bcf = NULL;
  bcf_hdr_t *_header = NULL;
  bcf1_t *_curr_rec = bcf_init();

  std::vector<std::string> _genotypes;
  std::vector<std::string> _alleles;

};
TEST_CASE ("VCF File handler") {
  using std::endl;
  std::string tmpvcf = "tmp_tc.vcf";

  // Write temp VCF file
  {
    std::ofstream vcfo(tmpvcf);
    vcfo
        << "##fileformat=VCFv4.1" << endl
        << "##phasing=true" << endl
        << "##contig=<ID=x>" << endl
        << "##contig=<ID=y>" << endl
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Freq\">" << endl
        << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternete Allele count\">" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Num samples at site\">" << endl
        << "##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Num alt alleles\">" << endl
        << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"Length of each alt\">" << endl
        << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"type of variant\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2" << endl
        << "x\t9\t.\tGG\tA,C,T\t99\t.\tAF=0.01,0.6,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t0|1\t2|3" << endl
        << "x\t10\t.\tC\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "x\t14\t.\tG\t<DUP>,<BLAH>\t99\t.\tAF=0.01,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t1|1" << endl
        << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
  }

      SUBCASE("Unfiltered") {
    VarFile vcf(tmpvcf);
        CHECK(vcf.num_samples() == 2);
        CHECK(vcf.sequences().size() == 2);
        CHECK(vcf.sequences()[0] == "x");
        CHECK(vcf.sequences()[1] == "y");

    // On load, first record is already loaded
        REQUIRE(vcf.genotypes().size() == 4);
        CHECK(vcf.genotypes()[0] == "GG");
        CHECK(vcf.genotypes()[1] == "A");
        CHECK(vcf.genotypes()[2] == "C");
        CHECK(vcf.genotypes()[3] == "T");
        REQUIRE(vcf.alleles().size() == 4);
        CHECK(vcf.alleles()[0] == "GG");
        CHECK(vcf.alleles()[1] == "A");
        CHECK(vcf.alleles()[2] == "C");
        CHECK(vcf.alleles()[3] == "T");
        CHECK(vcf.ref() == "GG");
        CHECK(vcf.pos() == 8);

    // Copy number alleles
    vcf.next();
        REQUIRE(vcf.genotypes().size() == 4);
        CHECK(vcf.genotypes()[0] == "CC");
        CHECK(vcf.genotypes()[1] == "CC");
        CHECK(vcf.genotypes()[2] == "");
        CHECK(vcf.genotypes()[3] == "CC");
        REQUIRE(vcf.alleles().size() == 3);
        CHECK(vcf.alleles()[0] == "C");
        CHECK(vcf.alleles()[1] == "CC");
        CHECK(vcf.alleles()[2] == "");
        CHECK(vcf.ref() == "C");
        CHECK(vcf.pos() == 9);

    // Invalid tags
    vcf.next();
        REQUIRE(vcf.alleles().size() == 3);
        CHECK(vcf.alleles()[0] == "G");
        CHECK(vcf.alleles()[1] == "G");
        CHECK(vcf.alleles()[2] == "G");
        CHECK(vcf.ref() == "G");
        CHECK(vcf.pos() == 13);

    // Next y contig should still load
    vcf.next();
        CHECK(vcf.alleles()[0] == "TATA");
  }

      SUBCASE("CHROM Filtering") {
    VarFile vcf;
    vcf.set_region("y:0-0");
    vcf.open(tmpvcf);

        CHECK(vcf.ref() == "TATA");
    vcf.next();
        CHECK(vcf.ref() == "T");
        CHECK(vcf.next() == 0); // File end
  }

      SUBCASE("Region filtering") {
    VarFile vcf;
    vcf.set_region("x:0-14");
    vcf.open(tmpvcf);

        CHECK(vcf.ref() == "GG");
    vcf.next();
        CHECK(vcf.ref() == "C");
    vcf.next();
        CHECK(vcf.ref() == "G");
        CHECK(vcf.next() == 0); // Region end
  }

  remove(tmpvcf.c_str());
}

class GraphBuilder {

 public:
  GraphBuilder(std::string reffile, std::string vcffile) : _fa_file(reffile), _vf_file(vcffile) { }

  void open(std::string ref, std::string vcf) {
    _fa_file = ref;
    _vf_file = vcf;
  }

  void region(std::string region) {
    _vf.set_region(region);
  }

  void region(std::string chr, int min, int max) {
    _vf.set_region(chr, min, max);
  }

  /**
   * Use a certain percentage of individuals. Reference nodes are always included.
   * @param percent, 0 - 100
   */
  void ingroup(int percent);

  /**
   * Only include nodes above this frequency. If >=1, only use
   * the node with the highest AF.
   * @param af minimum allele frequency
   */
  void min_af(float af) {
    if (af < 0) return;
    _min_af = af;
  }

  /**
   * Set maximum node length. If <= 0, length is unbounded.
   * @param max maximum node length.
   */
  void node_len(int max) { _max_node_len = max; }

  /**
   * Apply the various parameters and build the graph.
   * @return pointer to graph.
   */
  std::shared_ptr<graph> build();

  /**
   * After initial build, if parameters are changed, use this
   * to recreate a graph. Otherwise old graph will be returned.
   * @return graph reflecting any construction changes
   */
  std::shared_ptr<graph> rebuild();

 private:
  std::string _fa_file, _vf_file;
  VarFile _vf;
  FASTAFile _fa;
  std::shared_ptr<graph> g = nullptr;

  // Graph construction parameters
  int _ingroup = 100; // percent of individuals to use. Ref nodes always included
  float _min_af = 0; // Minimum AF
  int _max_node_len = 100000;
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
