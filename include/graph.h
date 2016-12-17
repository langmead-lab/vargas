/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date June 26, 2016
 *
 * @brief
 * Implementation of a directed graph.
 * @details
 * Each node stores a sequence and relevant
 * information. Graphs can be derived from other graphs with a filter, allowing
 * the extraction of population subsets.
 *
 * @file
 */

#ifndef VARGAS_GRAPH_H
#define VARGAS_GRAPH_H

#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>
#include "fasta.h"
#include "varfile.h"
#include "utils.h"
#include "dyn_bitset.h"

namespace vargas {

/**
 * @brief
 * Represents a Graph of the genome.
 * @details
 * The Graph is backed by a map of Graph::Node, and edges
 * are backed by a map of node ID's. If a Graph is built using a previous graph, the underlying
 * nodes of the origin graph are used by the derived graph. shared_ptr's are used to preserve
 * lifetimes. \n
 * Usage: \n
 * @code{.cpp}
 * #include "graph.h"
 *
 * //     GGG
 * //    /   \
 * // AAA     TTT
 * //    \   /
 * //     CCC(ref)
 *
 * // Node ID's are incremented each instantiation
 * Vargas::Graph::Node::_newID = 0;
 * Vargas::Graph g;

 * {
 * 	Vargas::Graph::Node n;
 * 	n.set_endpos(3);
 * 	n.set_as_ref();
 * 	std::vector<bool> a = {1, 1, 1};
 * 	n.set_population(a);
 * 	n.set_seq("AAA");
 * 	g.add_node(n);
 * }

 * {
 * 	Vargas::Graph::Node n;
 * 	n.set_endpos(6);
 * 	n.set_as_ref();
 * 	std::vector<bool> a = {1, 1, 1};
 * 	n.set_population(a);
 * 	n.set_af(0.4);
 * 	n.set_seq("CCC");
 * 	g.add_node(n);
 * }

 * {
 * 	Vargas::Graph::Node n;
 * 	n.set_endpos(6);
 * 	n.set_not_ref();
 * 	std::vector<bool> a = {0, 1, 0}; // Individual 1 has this allele
 * 	n.set_population(a);
 * 	n.set_af(0.6);
 * 	n.set_seq("GGG");
 * 	g.add_node(n);
 * }

 * {
 * 	Vargas::Graph::Node n;
 * 	n.set_endpos(9);
 * 	n.set_as_ref();
 * 	std::vector<bool> a = {1, 1, 1};
 * 	n.set_population(a);
 * 	n.set_seq("TTT");
 * 	n.set_af(0.3);
 * 	g.add_node(n);
 * }

 * g.add_edge(0, 1);
 * g.add_edge(0, 2);
 * g.add_edge(1, 3);
 * g.add_edge(2, 3);
 *
 * // Traverse nodes
 * for(Vargas::Graph::FilteringIter i = g.begin(); i != g.end(); ++i)
 *  std::cout << i->id() << " "; // 0 1 2 3 4
 *
 * // Derive a graph
 * std::vector<bool> filter = {1, 0, 0};
 * Vargas::Graph g2(g, filter); // AAA-CCC-TTT
 * Vargas::Graph g2(g, Vargas::Graph::REF); // Also AAA-CCC-TTT
 *
 * @endcode
 */
  class Graph {

    public:

      /**
       * @brief
       * Uniquely identifies a subgraph when mapped to a Population.
       * Default constructor indicates whole graph
       */
      struct GID {

          /**
           * @brief
           * Default GID assumes 100% ingroup.
           */
          GID() : num(100), id(0), pct(true), outgroup(false) {}

          /**
           * @param num if percent, set pct=true. If number of individuals, set pct=false.
           * @param id Unique ID for a given num
           * @param pct True if num is a percentage
           */
          GID(int num, int id, bool pct = false) : num(num), id(id), pct(pct), outgroup(false) {}

          /**
           * @brief
           * Build a GID from a string output from to_string
           * @param s string form
           */
          GID(std::string s);

          int num; /**< Percent or number of individuals included in the graph. */
          int id; /**< unique id if multiple graphs of num exist. */
          bool pct; /**< if true, num is a percentage. Otherwise number of individuals.*/
          bool outgroup; /**< True if the origin was an outgroup graph.*/

          std::string to_string() const;
      };

      using Population = VCF::Population;

      /**
        * @enum Type
        * Indicate a special type of graph, one that only incudes REF nodes
        * or MAXAF nodes.
        */
      enum class Type: char {
          REF, /**< Only include reference nodes. */
          MAXAF, /**< Keep the node with the highest AF at each branch. */
      };

      /**
       * @brief
       * Represents a node in the directed graphs.
       * @details
       * Sequences are stored numerically.
       * populations are stored as bitsets, where 1 indicates that individual posesses
       * the given allele.
       */
      class Node {
        public:
          /**
           * @brief
           * Make new node, and assign a unique ID.
           */
          Node() : _id(_newID++) {}

          Node(size_t pos, const std::string &seq, const Population &pop, bool ref, float af) :
          _endPos(pos), _seq(seq_to_num(seq)), _individuals(pop), _ref(ref), _af(af), _id(_newID++) {}

          /**
           * @return length of the sequence
           */
          size_t length() const { return _seq.size(); }

          /**
           * @return position of last base in seq, 0 indexed
           */
          size_t end_pos() const { return _endPos; }

          size_t begin_pos() const { return _endPos - _seq.size() + 1; }

          /**
           * @brief
           * Return true if the node is a ref node, or if the bit is set.
           * @param idx bit index
           * @return belongs
           */
          bool belongs(uint idx) const { return _individuals.at(idx); }

          /**
           * @brief
           * Returns true if the node is a ref node or if the population intersects the node population.
           * @param pop Population filter
           * @return belongs
           */
          bool belongs(const Population &pop) const { return pop && _individuals; }

          /**
           * @brief
           * Sequence as a vector of unsigned chars.
           * @return seq
           */
          const std::vector<Base> &seq() const { return _seq; }

          /**
           * @brief
           * Sequence is stored numerically. Return as a string.
           * @return seq
           */
          std::string seq_str() const { return num_to_seq(_seq); }

          /**
           * @brief
           * Size of the Population of the node. This should be consistant throughout the graph.
           * @return pop_size
           */
          size_t pop_size() const { return _individuals.size(); }

          /**
           * @brief
           * Node id.
           * @return unique node ID
           */
          size_t id() const { return _id; }

          /**
           * @brief
           * True if node common to all individuals, or REF.
           * @return is_ref
           */
          bool is_ref() const { return _ref; }

          /**
           * @brief
           * Allele frequency.
           * @return af
           */
          float freq() const { return _af; }

          /**
           * @brief
           * Reference to the raw Population data member.
           * @return individuals
           */
          const Population &individuals() const { return _individuals; }

          /**
           * @brief
           * Set the id of the node, should rarely be used as unique ID's are generated.
           * @param id
           */
          void setID(const size_t id) {
              if (id >= _newID) {
                  this->_id = id;
                  _newID = id + 1;
              }
          }

          /**
           * @brief
           * Set the position of the last base in the sequence.
           * @param pos 0-indexed
           */
          void set_endpos(const size_t pos) { this->_endPos = pos; }

          /**
           * @brief
           * Set the population from an existing Population
           * @param pop
           */
          void set_population(const Population &pop) { _individuals = pop; }

          /**
           * @brief
           * Set the population from a vector. Each individual is set if pop[i] evaluates true.
           * @param pop
           */
          template<typename T>
          void set_population(const std::vector<T> &pop) { _individuals = pop; }

          /**
           * @brief
           * Set the population.
           * @param len number of genotypes
           * @param val true/false for each individual
           */
          void set_population(size_t len, bool val) { _individuals = Population(len, val); }

          /**
           * @brief
           * Set the stored node sequence. Sequence is converted to numeric form.
           * @param seq
           */
          void set_seq(const std::string &seq) { _seq = seq_to_num(seq); }

          /**
           * @brief
           * Set the stored node sequence
           * @param seq
           */
          void set_seq(std::vector<Base> &seq) { this->_seq = seq; }

          /**
           * @brief
           * Sets the node as a reference node, and sets all bits in the population.
           */
          void set_as_ref() {
              _ref = true;
              _individuals.set();
          }

          /**
           * @brief
           * Deselects node as a referene node, does not modify population.
           */
          void set_not_ref() { _ref = false; }

          /**
           * @brief
           * Set the allele frequency of the node. Used for MAXAF filtering.
           * @param af float frequency, between 0 and 1
           */
          void set_af(float af) { _af = af; }

          /**
           * @brief
           * Set as pinch node. If this node is removed, then each node before
           * it does not have an edge to a node after this node. All paths traverse
           * through this node.
           */
          void pinch() { _pinch = true; }

          void unpinch() { _pinch = false; }
          /**
           * @brief
           * If true, then previous alignment seeds can be cleared as no nodes after the current
           * one depends on nodes before current one,in a topographical ordering.
           * @return true if pinched
           */
          bool is_pinched() const { return _pinch; }

          static size_t _newID; /**< ID of the next instance to be created */

        private:
          size_t _endPos; // End position of the sequence
          std::vector<Base> _seq; // sequence in numeric form
          Population _individuals; // Each bit marks an individual, 1 if they have this node
          bool _ref = false; // Part of the reference sequence if true
          bool _pinch = false; // If this node is removed, the graph will split into two distinct subgraphs
          float _af = 1;
          size_t _id;

      };

      typedef std::shared_ptr<Node> nodeptr;
      typedef std::shared_ptr<const Node> const_nodeptr;

      /**
       * @brief
       * Forward iterator to traverse the Graph in insertion order. Checks assume underlying graphs are equal.
       */
      template<typename T, typename Unqualified_T = typename std::remove_cv<T>::type>
      class GraphIterator: public std::iterator<std::forward_iterator_tag, Unqualified_T, std::ptrdiff_t, T *, T &> {
        public:

          /**
           * @brief
           * @param g Graph
           * @param idx Node in the insertion order to begin iterator at.
           */
          GraphIterator(const Graph &g, const size_t idx = 0) : _graph(g), _currID(idx) {}

          /**
           * @brief
           * Reference to the underlying Graph.
           * @return const ref to graph
           */
          const Graph &graph() const { return _graph; }

          /**
           * @return True if iterators point to same node index.
           */
          template<typename O>
          bool operator==(const GraphIterator<O> &other) const {
              return _currID == other._currID;
          }

          /**
           * @return true when current Node ID's are not the same or if the
           * underlying graph is not the same.
           */
          template<typename O>
          bool operator!=(const GraphIterator<O> &other) const {
              return _currID != other._currID;
          }

          /**
           * @brief
           * Goes to next Node.
           * @return iterator to the next Node.
           */
          GraphIterator &operator++() {
              if (_currID < _graph._add_order.size()) ++_currID;
              return *this;
          }

          /**
           * @brief
           * Goes to next Node, returns previous node.
           * @return iterator to the next Node.
           */
          GraphIterator operator++(int) {
              auto ret = *this;
              if (_currID < _graph._add_order.size()) ++_currID;
              return ret;
          }

          /**
           * @brief
           * Const reference to the current node. Undefined for end iterator.
           * @return Node
           */
          T &operator*() const {
              return *(_graph._IDMap->at(_graph._add_order[_currID]));
          }

          /**
           * @return pointer to underlying node
           */
          T *operator->() const {
              return &operator*();
          }

          /**
           * @brief
           * All nodes that we've traversed that have incoming edges to the current node.
           * @return vector of previous nodes
           */
          const std::vector<size_t> &incoming() const {
              try {
                  // Assume the previous edge existing is the common case (DAG w 1 start node)
                  return _graph._prev_map.at(_graph._add_order[_currID]);
              } catch (std::exception &e) { return _empty_vec; }
          }

          /**
           * @return vector of all outgoing edges
           */
          const std::vector<size_t> &outgoing() const {
              try {
                  // Assume the previous edge existing is the common case (DAG w 1 start node)
                  return _graph._next_map.at(_graph._add_order[_currID]);
              } catch (std::exception &e) { return _empty_vec; }
          }

          //TODO icc has a problem with this
          /**
           * @brief
           * Allow conversion from iterator to const_iterator
           */
          operator GraphIterator<const T>() const {
              return GraphIterator<const T>(_graph, _currID);
          }


        private:

          const Graph &_graph;
          size_t _currID;
          const std::vector<size_t> _empty_vec;

      };

      using iterator = GraphIterator<Graph::Node>;
      using const_iterator = GraphIterator<const Graph::Node>;

      /**
       * @brief
       * Default constructor inits a new Graph, including a new node map.
       */
      Graph() : _IDMap(std::make_shared<std::unordered_map<size_t, nodeptr>>()) {}

      /**
       * @brief
       * Build a graph given a reference FASTA file and a variant VCF or BCF file.
       * @param ref_file FASTA file name
       * @param vcf_file VCF or BCF file name
       * @param region region in the format chromosome:start-end
       * @param max_node_len Maximum graph node length
       */
      Graph(const std::string &ref_file, const std::string &vcf_file,
            const std::string &region, const int max_node_len = 1000000);

      /**
       * @brief
       * Create a Graph with another Graph and a population filter.
       * @details
       * The new Graph will only
       * contain nodes if any of the individuals in filter possess the node. The actual nodes
       * are shared_ptr's to the parent Graph, as to prevent duplication of Nodes.
       * @param g Graph to derive the new Graph from
       * @param filter population filter, only include nodes representative of this population
       */
      Graph(const Graph &g, const Population &filter);

      /**
        * @brief
        * Construct a graph using a base graph and a filter.
        * @details
        * Constructs a graph using filter tags, one of:
        * Graph::Type::REF, or Graph::Type::MAXAF
        * The former keeps reference nodes, the later picks the node with the highest allele
        * frequency. Both result in linear graphs.
        * @param g Graph to derive from
        * @param t one of Graph::Type::REF, Graph::Type::MAXAF
        */
      Graph(const Graph &g, Type t);

      /**
       * @brief
       * Add a new node to the Graph.
       * @details
       * A new node is created so the original can be destroyed.
       * The first node added is set as the Graph root. Nodes must be added in topographical order.
       * @param n node to add, ID of original node is preserved.
       * @return ID of the inserted node
       */
      size_t add_node(const Node &n);

      /**
       * @brief
       * Create an edge linking two nodes. Previous and Next edges are added.\n
       * n1->n2
       * @param n1 Node one ID
       * @param n2 Node two ID
       */
      bool add_edge(const size_t n1, const size_t n2);

      /**
       * @brief
       * Sets the root of the Graph.
       * @param id ID of root node
       */
      void set_root(const size_t id) {
          _root = id;
      }

      /**
       * @brief
       * Set the graph description.
       * @param description
       */
      void set_desc(const std::string &description) { _desc = description; }

      /**
       * @brief
       * Return root node ID
       * @return root
       */
      size_t root() const { return _root; }

      /**
       * @brief
       * Maps a ndoe ID to a shared node object
       * @return map of ID, shared_ptr<Node> pairs
       */
      std::shared_ptr<const std::unordered_map<size_t, nodeptr>> node_map() const { return _IDMap; }

      /**
       * @brief
       * Maps a node ID to a vector of all next nodes (outgoing edges)
       * @return map of ID, outgoing edge vectors
       */
      const std::unordered_map<size_t, std::vector<size_t>> &next_map() const { return _next_map; }

      /**
       * @brief
       *  Maps a node ID to a vector of all incoming edge nodes
       *  @return map of ID, incoming edges
       */
      const std::unordered_map<size_t, std::vector<size_t>> &prev_map() const { return _prev_map; }

      /**
       * @brief
       * Const reference to a node
       * @param id ID of the node
       * @return shared node object
       * @throws std::invalid_argument if node ID does not exist
       */
      const Node &node(size_t id) const {
          if (_IDMap->count(id) == 0) throw std::domain_error("Invalid Node ID.");
          return *(*_IDMap).at(id);
      }

      /**
       * @return description of graph
       */
      std::string desc() const { return _desc; }

      /**
       * @brief
       * Exports the graph in DOT format.
       * @param name graph name
       */
      std::string to_DOT(const std::string name = "g") const;

      /**
       * @brief
       * Export the graph in DOT format
       * @param filename to export to
       * @param name of the graph
       * @throws std::invalid_argument if output file cannot be opened
       */
      void to_DOT(const std::string filename, const std::string name) const {
          std::ofstream out(filename);
          if (!out.good()) throw std::invalid_argument("Error opening file: \"" + filename + "\"");
          out << to_DOT(name);
      }

      /**
       * @return DOT repersentation of graph.
       */
      std::string to_string() const {
          return to_DOT();
      }


      /**
       * @brief
       * Define the graph population size.
       * @param popsize number of genotypes
       */
      void set_popsize(const size_t popsize) { _pop_size = popsize; }

      /**
       * @brief
       * Associate a filter with the graph.
       * @param filter
       */
      void set_filter(const Population &filter) { _filter = filter; }

      /**
       * @return assosciated filter
       */
      const Population &filter() const { return _filter; }

      /**
       * @return population size
       */
      size_t pop_size() const { return _pop_size; }

      /**
       * @brief
       * Return a Population of a subset of the graph.
       * @return Population with ingroup % indivduals set.
       */
      Population subset(const int ingroup) const;

      /**
       * @return begin iterator.
       */
      iterator begin() const {
          return iterator(*this, 0);
      }

      /**
       * @return end iterator.
       */
      iterator end() const {
          return iterator(*this, _add_order.size());
      }

      /**
     * @return const begin iterator.
     */
      const_iterator cbegin() const {
          return const_iterator(*this, 0);
      }

      /**
       * @return const end iterator.
       */
      const_iterator cend() const {
          return const_iterator(*this, _add_order.size());
      }

      /**
       * @brief
       * Create a subgraph including bases from min to max.
       * @param min
       * @param max
       * @return Graph
       */
      Graph subgraph(const size_t min, const size_t max) const;

      /**
       * @brief
       * Set the maximum node length of the graph.
       * @param len max node length
       */
      void node_len(size_t len) { _max_node_len = len; }

      /**
       * @return max node len
       */
      size_t max_node_len() const { return _max_node_len; }

      /**
       * @brief
       * Ensures that the graph is topologically sorted.
       * @param begin Graph begin iterator
       * @param end Graph end iterator
       */
      bool validate() const;


    private:
      size_t _root = 0; // Root of the Graph
      // maps a node ID to a nodeptr. Any derived graphs use the same base node ID map.
      std::shared_ptr<std::unordered_map<size_t, nodeptr>> _IDMap;
      // maps a node ID to the vector of nodes it points to
      std::unordered_map<size_t, std::vector<size_t>> _next_map;
      // maps a node ID to a vector of node ID's that point to it
      std::unordered_map<size_t, std::vector<size_t>> _prev_map;
      std::vector<size_t> _add_order; // Order nodes were added
      // Description, used by the builder to store construction params
      std::string _desc;
      size_t _pop_size = 0;
      size_t _max_node_len;
      Population _filter;

      /**
       * Given a subset of nodes from Graph g, rebuild all applicable edges in the new graph.
       * @param g underlying parent graph
       * @param includedNodes subset of g's nodes to include
       */
      void _build_derived_edges(const Graph &g, const std::unordered_map<size_t, nodeptr> &includedNodes);

  };

  /**
   * @brief
   * Takes a reference sequence and a variant file and builds a graph.
   * @details
   * The base graph can
   * include a subset of samples, or a full graph can be built and subsequent graphs derived
   * from the base graph. \n
   * Usage: \n
   * @code{.cpp}
   * #include "graph.h"
   *
   * Vargas::GraphFactory gb("reference.fa", "var.bcf");
   * gb.node_len(5);
   * gb.ingroup(100);
   * gb.region("x:0-15");
   *
   * Vargas::Graph g = gb.build();
   * @endcode
   */
  class GraphFactory {

    public:

      /**
       * @brief
       * Create a graph builder for the reference file.
       * @param reffile
       */
      GraphFactory(std::string const &reffile) : _fa_file(reffile) {}

      /**
       * @param fafile FASTA file
       * @param varfile VCF or BCF file.
       */
      GraphFactory(std::string const &fafile, std::string const &varfile) : _fa_file(fafile) {
          open_bcf(varfile);
      }

      /**
       * @brief
       * Set the region of the graph to build. Format should be
       * CHR:XX,XXX-YY,YYY
       * Where CHR is the sequence name, XX,XXX is the min pos and YY,YYY is the max pos.
       * Both are inclusive.
       * @param region
       */
      void set_region(std::string region) {
          if (!_vf) throw std::invalid_argument("No variant file opened.");
          _vf->set_region(region);
      }

      /**
       * @brief
       * Set the region of the sequence to build.
       * @param chr Chromosome/sequence of build
       * @param min min pos, inclusive
       * @param max max pos, inclusive
       */
      void set_region(std::string chr, int min, int max) {
          if (!_vf) throw std::invalid_argument("No variant file opened.");
          _vf->set_region(chr, min, max);
      }

      /**
       * @brief
       * Set maximum node length. If <= 0, length is unbounded.
       * @param max maximum node length.
       */
      void node_len(size_t max) { _max_node_len = max; }

      /**
       * Open the given file
       * @param file_name
       * @return Number of samples
       */
      size_t open_vcf(std::string const &file_name);

      /**
       * @brief
       * Limit the samples used from the VCF file.
       * @param filter CSV list of sample names to keep.
       * @param invert Use the samples not specified in filter
       * @return Number of samples in the filter.
       */
      size_t add_sample_filter(std::string filter, bool invert = false);

      /**
       * Open the given file
       * @param file_name
       * @return Number of samples
       */
      size_t open_bcf(std::string const &file_name) {
          return open_vcf(file_name);
      }

      /**
       * @brief
       * Apply the various parameters and build the Graph from a VCF and a FASTA.
       * @param g Graph to build into
       */
      void build(Graph &g);

      /**
       * @brief
       * Build the graph using the specified params from a VCF and a FASTA.
       * @return Graph Built Graph
       */
      Graph build() {
          Graph g;
          build(g);
          return g;
      }

    protected:
      /**
       * @brief
       * Builds edges between pending previous nodes and the current level. Each
       * prev is connected to each curr, and curr becomes the new prev.
       * @param g build edges for Graph g
       * @param prev previous unconnected nodes (linked to main graph already)
       * @param curr current unconnected nodes
       * @param chain if a single linear node is split, map the beginning of the sequence to the end.
       */
      __RG_STRONG_INLINE__
      void _build_edges(Graph &g,
                        std::unordered_set<size_t> &prev,
                        std::unordered_set<size_t> &curr,
                        std::unordered_map<size_t, size_t> *chain = nullptr);

      /**
       * @brief
       * Builds a linear sequence of nodes set as reference nodes.
       * @param g Graph to build linear ref in
       * @param prev previous unconnected nodes (linked to main graph already)
       * @param curr current unconnected nodes
       * @param pos current position, inclusive
       * @param target build linear sequence up to this position, exclusive
       * @return ending position
       */
      __RG_STRONG_INLINE__
      int _build_linear_ref(Graph &g,
                            std::unordered_set<size_t> &prev,
                            std::unordered_set<size_t> &curr,
                            size_t pos,
                            size_t target);

      /**
       * @brief
       * Splits a sequence into multiple sequences if neccessary to conform to the
       * max_node_len spec.
       * @return vector of split sequences
       */
      __RG_STRONG_INLINE__
      std::vector<std::string> _split_seq(std::string seq);

    private:
      std::string _fa_file;
      std::unique_ptr<VariantFile> _vf;
      ifasta _fa;
      Graph g;

      // Graph construction parameters
      size_t _max_node_len = 10000000; // If a node is longer, split into multiple nodes
  };

  /**
 * @brief
 * Operator used to map Graph::GID's.
 * @param a Graph::GID a
 * @param b Graph::GID b
 */
  bool operator<(const Graph::GID &a, const Graph::GID &b);

  /**
   * @brief
   * Outputs the Graph::GID in a CSV format:\n
   * [o,i],num,id,pct\n
   * @param os Output stream
   * @param gid gid to print
   */
  std::ostream &operator<<(std::ostream &os, const Graph::GID &gid);

  /**
   * @brief
   * Check if two GID's are equal.
   */
  bool operator==(const Graph::GID &a, const Graph::GID &b);

}

#endif //VARGAS_GRAPH_H
