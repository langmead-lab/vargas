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

#include <cstdio>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <bitset>
#include <thread>
#include "fasta.h"
#include "varfile.h"
#include "doctest.h"
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
          GID() : num(100), id(0), pct(true), outgroup(false) {}

          /**
           * @brief
           * @param num if percent, set pct=true. If number of individuals, set pct=false.
           * @param id Unique ID for a given num
           * @param pct True if num is a percentage
           */
          GID(int num,
              int id,
              bool pct = false) :
              num(num), id(id), pct(pct), outgroup(false) {}

          /**
           * @brief
           * Build a GID from a string output from to_string
           * @param s string form
           */
          GID(std::string s) {
              try {
                  std::vector<std::string> s_split = split(s, ',');
                  outgroup = s[0] == 'o';
                  num = std::stoi(s_split[1]);
                  id = std::stoi(s_split[2]);
                  pct = s_split[3][0] == '1';
              }
              catch (std::exception &e) {
                  std::cerr << "Invalid GID string format: " << s << std::endl;
                  GID();
              }

          }

          int num; /**< Percent or number of individuals included in the graph. */
          int id; /**< unique id if multiple graphs of num exist. */
          bool pct; /**< if true, num is a percentage. Otherwise number of individuals.*/
          bool outgroup; /**< True if the origin was an outgroup graph.*/

          std::string to_string() const {
              std::ostringstream ss;
              ss << (outgroup ? 'o' : 'i') << ',' << num << ',' << id << ',' << pct;
              return ss.str();
          }
      };

      /** Alias for a dyn_bitset */
      typedef VCF::Population Population;

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

          Node(int pos,
               const std::string &seq,
               const Population &pop,
               bool ref,
               float af) :
              _endPos(pos), _seq(seq_to_num(seq)), _individuals(pop), _ref(ref), _af(af), _id(_newID++) {}

          /**
           * @return length of the sequence
           */
          size_t length() const {
              return _seq.size();
          }

          /**
           * @return position of last base in seq, 0 indexed
           */
          int end() const {
              return _endPos;
          }

          /**
           * @brief
           * Return true if the node is a ref node, or if the bit is set.
           * @param idx bit index
           * @return belongs
           */
          bool belongs(uint idx) const {
              return _individuals.at(idx);
          }

          /**
           * @brief
           * Returns true if the node is a ref node or if the population intersects the node population.
           * @param pop Population filter
           * @return belongs
           */
          bool belongs(const Population &pop) const {
              return pop && _individuals;
          }

          /**
           * @brief
           * Sequence as a vector of unsigned chars.
           * @return seq
           */
          const std::vector<Base> &seq() const {
              return _seq;
          }

          /**
           * @brief
           * Sequence is stored numerically. Return as a string.
           * @return seq
           */
          std::string seq_str() const {
              return num_to_seq(_seq);
          }

          /**
           * @brief
           * Size of the Population of the node. This should be consistant throughout the graph.
           * @return pop_size
           */
          size_t pop_size() const {
              return _individuals.size();
          }

          /**
           * @brief
           * Node id.
           * @return unique node ID
           */
          uint32_t id() const {
              return _id;
          }

          /**
           * @brief
           * True if node common to all individuals, or REF.
           * @return is_ref
           */
          bool is_ref() const {
              return _ref;
          }

          /**
           * @brief
           * Allele frequency.
           * @return af
           */
          float freq() const {
              return _af;
          }

          /**
           * @brief
           * Reference to the raw Population data member.
           * @return individuals
           */
          const Population &individuals() const {
              return _individuals;
          }

          /**
           * @brief
           * Set the id of the node, should rarely be used as unique ID's are generated.
           * @param id
           */
          void setID(uint32_t id) {
              if (id >= _newID) {
                  this->_id = id;
                  _newID = ++id;
              }
          }

          /**
           * @brief
           * Set the position of the last base in the sequence.
           * @param pos 0-indexed
           */
          void set_endpos(int pos) {
              this->_endPos = pos;
          }

          /**
           * @brief
           * Set the population from an existing Population
           * @param pop
           */
          void set_population(const Population &pop) {
              _individuals = pop;
          }

          /**
           * @brief
           * Set the population from a vector. Each individual is set if pop[i] evaluates true.
           * @param pop
           */
          template<typename T>
          void set_population(const std::vector<T> &pop) {
              _individuals = pop;
          }

          /**
           * @brief
           * Set the population.
           * @param len number of genotypes
           * @param val true/false for each individual
           */
          void set_population(size_t len,
                              bool val) {
              _individuals = Population(len, val);
          }

          /**
           * @brief
           * Set the stored node sequence. Sequence is converted to numeric form.
           * @param seq
           */
          void set_seq(const std::string &seq) {
              _seq = seq_to_num(seq);
          }

          /**
           * @brief
           * Set the stored node sequence
           * @param seq
           */
          void set_seq(std::vector<Base> &seq) {
              this->_seq = seq;
          }

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
          void set_not_ref() {
              _ref = false;
          }

          /**
           * @brief
           * Set the allele frequency of the node. Used for MAXAF filtering.
           * @param af float frequency, between 0 and 1
           */
          void set_af(float af) {
              _af = af;
          }

          /**
           * @brief
           * Set as pinch node. If this node is removed, then each node before
           * it does not have an edge to a node after this node. All paths traverse
           * through this node.
           */
          void pinch() {
              _pinch = true;
          }

          void unpinch() {
              _pinch = false;
          }
          /**
           * @brief
           * If true, then previous alignment seeds can be cleared as no nodes after the current
           * one depends on nodes before current one,in a topographical ordering.
           * @return true if pinched
           */
          bool is_pinched() const {
              return _pinch;
          }

          static uint32_t _newID; /**< ID of the next instance to be created */

        private:
          int _endPos; // End position of the sequence
          std::vector<Base> _seq; // sequence in numeric form
          Population _individuals; // Each bit marks an individual, 1 if they have this node
          bool _ref = false; // Part of the reference sequence if true
          bool _pinch = false; // If this node is removed, the graph will split into two distinct subgraphs
          float _af = 1;
          uint32_t _id;

      };

      typedef std::shared_ptr<Node> nodeptr;

      /**
       * @brief
       * Default constructor inits a new Graph, including a new node map.
       */
      Graph() : _IDMap(std::make_shared<std::unordered_map<uint32_t, nodeptr>>(std::unordered_map<uint32_t,
                                                                                                  nodeptr>())) {}

      /**
       * @brief
       * Build a graph given a reference FASTA file and a variant VCF or BCF file.
       * @param ref_file FASTA file name
       * @param vcf_file VCF or BCF file name
       * @param region region in the format chromosome:start-end
       * @param max_node_len Maximum graph node length
       */
      Graph(std::string ref_file,
            std::string vcf_file,
            std::string region,
            int max_node_len = 1000000);

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
      Graph(const Graph &g,
            const Population &filter);

      /**
       * @brief
       * Add a new node to the Graph.
       * @details
       * A new node is created so the original can be destroyed.
       * The first node added is set as the Graph root. Nodes must be added in topographical order.
       * @param n node to add, ID of original node is preserved.
       */
      uint32_t add_node(Node &n);

      /**
       * @brief
       * Create an edge linking two nodes. Previous and Next edges are added.\n
       * n1->n2
       * @param n1 Node one ID
       * @param n2 Node two ID
       */
      bool add_edge(uint32_t n1,
                    uint32_t n2);

      /**
       * @brief
       * Sets the root of the Graph.
       * @param id ID of root node
       */
      void set_root(uint32_t id) {
          _root = id;
      }

      /**
       * @brief
       * Set the graph description.
       * @param description
       */
      void set_desc(std::string description) { _desc = description; }

      /**
       * @brief
       * Return root node ID
       * @return root
       */
      uint32_t root() const { return _root; }

      /**
       * @brief
       * Maps a ndoe ID to a shared node object
       * @return map of ID, shared_ptr<Node> pairs
       */
      const std::shared_ptr<std::unordered_map<uint32_t, nodeptr>> &node_map() const { return _IDMap; }

      /**
       * @brief
       * Maps a node ID to a vector of all next nodes (outgoing edges)
       * @return map of ID, outgoing edge vectors
       */
      const std::unordered_map<uint32_t, std::vector<uint32_t>> &next_map() const { return _next_map; }

      /**
       * @brief
       *  Maps a node ID to a vector of all incoming edge nodes
       *  @return map of ID, incoming edges
       */
      const std::unordered_map<uint32_t, std::vector<uint32_t>> &prev_map() const { return _prev_map; }

      /**
       * @brief
       * Const reference to a node
       * @param id ID of the node
       * @return shared node object
       * @throws std::invalid_argument if node ID does not exist
       */
      const Node &node(uint32_t id) const {
          if (_IDMap->count(id) == 0) throw std::invalid_argument("Invalid Node ID.");
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
      std::string to_DOT(std::string name = "g") const;

      /**
       * @brief
       * Export the graph in DOT format
       * @param filename to export to
       * @param name of the graph
       * @throws std::invalid_argument if output file cannot be opened
       */
      void to_DOT(std::string filename, std::string name) const {
          std::ofstream out(filename);
          if (!out.good()) throw std::invalid_argument("Error opening file: \"" + filename + "\"");
          out << to_DOT(name);
      }

      /**
       * @brief
       * Define the graph population size.
       * @param popsize number of genotypes
       */
      void set_popsize(int popsize) { _pop_size = popsize; }

      void set_filter(const Population &filter) { _filter = filter; }
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
      Population subset(int ingroup) const;

      /**
       * @brief
       * const forward iterator to traverse the Graph while applying a filter.
       * @details
       * When incrementing the iterator, only nodes that match a condition will
       * return.
       * Options include: \n
       * Filtering: Provided a population, a node is returned if there is an intersection \n
       * REF: Only return reference nodes \n
       * MAXAF: Return the node with the highest allele frequency. \n
       */
      class GraphIterator {
        public:

          /**
         * @enum Type
         * When a normal population filter is not used, a flag can be used. REF includes
         * only reference alleles, MAXAF picks the allele with the highest frequency.
         * Both result in linear graphs.
         * Filter used as a placeholder, it should never be passed as a param.
         */
          enum class Type : char {
              TOPO, /**< Use a topographical ordering.*/
              REF, /**< Only include reference nodes. */
              MAXAF, /**< Keep the node with the highest AF at each branch. */
              FILTER, /**< Use a given Population filter. */
              END  /**< End iterator. */
          };

          /**
           * @brief
           * Accept all nodes.
           * @ param g Graph
           */
          explicit GraphIterator(const Graph &g) :
              _graph(g), _type(Type::TOPO), _currID(0) {}

          /**
           * @brief
           * Traverse nodes when filter has at least one individual in common with the
           * Node Population.
           * @param g graph
           * @param filter Population the node is checked against. If there is an intersection,
           * the Node is returned.
           */
          explicit GraphIterator(const Graph &g, const Population &filter) :
              _graph(g), _filter(filter), _type(Type::FILTER), _currID(g._root) {}

          /**
           * @brief
           * Traverse a linear subgraph.
           * @param g graph
           * @param type one of Graph::MAXAF, Graph::REF
           */
          explicit GraphIterator(const Graph &g, Type type)
              : _graph(g), _type(type), _currID(g._root) {
              if (_type == Type::TOPO) _currID = 0;
              if (_type == Type::END) _currID = g._add_order.size() - 1;
          }

          /**
           * @brief
           * Reference to the underlying Graph.
           * @return const ref to graph
           */
          const Graph &graph() const { return _graph; }

          /**
           * @return true when underlying Graph address is the same and current Node ID's
           * are the same. Two end iterators always compare equal.
           */
          bool operator==(const GraphIterator &other) const;

          /**
           * @return true when current Node ID's are not the same or if the
           * underlying graph is not the same. Two end iterators always compare
           * false.
           */
          bool operator!=(const GraphIterator &other) const;

          /**
           * @brief
           * Goes to next Node. Nodes are only included if it satisfies the filter.
           * Once the end of the graph is reached, _end is set.
           * @return iterator to the next Node.
           */
          GraphIterator &operator++();

          /**
           * @return iterator filtering Type
           */
          Type type() const { return _type; }

          /**
           * @brief
           * Const reference to the current node. Undefined for end iterator.
           * @return Node
           */
          const Graph::Node &operator*() const;

          /**
           * @brief
           * @return pointer to underlying node
           */
          const Graph::Node *operator->() const {
              return &operator*();
          }

          /**
           * @brief
           * All nodes that we've traversed that have incoming edges to the current node.
           * @return vector of previous nodes
           */
          const std::vector<uint32_t> &incoming();

          /**
           * @return vector of all outgoing edges
           */
          const std::vector<uint32_t> &outgoing();


        private:
          /**
           * @brief
           * Inserts the ID into the queue if it's unique.
           * @param id node id to insert
           */
          __RG_STRONG_INLINE__
          void _insert_queue(uint32_t id) {
              if (_queue_unique.count(id) == 0) {
                  _queue_unique.insert(id);
                  _queue.push(id);
              }
          }

          const Graph &_graph; // Underlying graph
          Population _filter; // Nodes that intersect with _filter are included if _type == FILTER
          Graph::GraphIterator::Type _type; // Set to MAXAF or REF for linear subgraph traversals
          uint32_t _currID;
          std::queue<uint32_t> _queue;
          std::unordered_set<uint32_t> _queue_unique; // Used to make sure we don't repeat nodes
          std::unordered_set<uint32_t> _traversed; // Set of all nodes we've passed
          std::vector<uint32_t> _incoming; // Set of incoming edges, for use by incoming()
          const std::vector<uint32_t> _outgoing; // empty vec
          size_t _add_order_size = _graph._add_order.size();

      };

      /**
     * @brief
     * Construct a graph using a base graph and a filter.
     * @details
     * Constructs a graph using filter tags, one of:
     * Graph::REF, or Graph::MAXAF
     * The former keeps reference nodes, the later picks the node with the highest allele
     * frequency. Both result in linear graphs.
     * @param g Graph to derive from
     * @param t one of Graph::REF, Graph::MAXAF
     */
      Graph(const Graph &g, GraphIterator::Type t);

      /**
       * @brief
       * Provides an iterator to the whole graph.
       * @return reference to the root node.
       */
      GraphIterator begin() const {
          return GraphIterator(*this);
      }

      /**
       * @brief
       * Iterator when conditions are applied.
       * @param filter population to compare nodes to
       */
      GraphIterator begin(const Population &filter) const {
          return GraphIterator(*this, filter);
      }

      /**
       * @brief
       * Linear subgraph interator.
       * @param type one of Graph::REF, Graph::MAXAF
       */
      GraphIterator begin(GraphIterator::Type type) const {
          return GraphIterator(*this, type);
      }

      /**
       * @return end iterator.
       */
      GraphIterator end() const {
          return GraphIterator(*this, Graph::GraphIterator::Type::END);
      }

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


    private:
      uint32_t _root = -1; // Root of the Graph
      // maps a node ID to a nodeptr. Any derived graphs use the same base node ID map.
      std::shared_ptr<std::unordered_map<uint32_t, nodeptr>> _IDMap;
      // maps a node ID to the vector of nodes it points to
      std::unordered_map<uint32_t, std::vector<uint32_t>> _next_map;
      // maps a node ID to a vector of node ID's that point to it
      std::unordered_map<uint32_t, std::vector<uint32_t>> _prev_map;
      std::vector<uint32_t> _add_order; // Order nodes were added
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
      void _build_derived_edges(const Graph &g,
                                const std::unordered_map<uint32_t, nodeptr> &includedNodes);

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
   * Vargas::GraphBuilder gb("reference.fa", "var.bcf");
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

      ~GraphFactory() { _vf.reset(); }

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
      void set_region(std::string chr,
                      int min,
                      int max) {
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
      int open_vcf(std::string const &file_name) {
          _vf.reset();
          _vf = std::unique_ptr<VariantFile>(new VCF(file_name));
          if (!_vf->good()) throw std::invalid_argument("Invalid VCF/BCF file: \"" + file_name + "\"");
          return _vf->num_samples();
      }

      /**
       * @brief
       * Limit the samples used from the VCF file.
       * @param filter CSV list of sample names to keep.
       * @param invert Use the samples not specified in filter
       * @return Number of samples in the filter.
       */
      int add_sample_filter(std::string filter, bool invert = false) {
          if (!_vf) throw std::invalid_argument("No VCF file opened, cannot add filter.");
          if (filter.length() == 0 || filter == "-") return _vf->num_samples();
          std::vector<std::string> filt;
          filter.erase(std::remove_if(filter.begin(), filter.end(), isspace), filter.end());
          if (invert) {
              const auto s = _vf->samples();
              auto vcf_samples = std::unordered_set<std::string>(s.begin(), s.end());
              const auto filter_samples = split(filter, ',');
              for (const auto &f : filter_samples) vcf_samples.erase(f);
              filt = std::vector<std::string>(vcf_samples.begin(), vcf_samples.end());
          } else {
              split(filter, ',', filt);
          }
          _vf->create_ingroup(filt);
          return _vf->num_samples();
      }

      /**
       * Open the given file
       * @param file_name
       * @return Number of samples
       */
      int open_bcf(std::string const &file_name) {
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
                        std::unordered_set<uint32_t> &prev,
                        std::unordered_set<uint32_t> &curr,
                        std::unordered_map<uint32_t, uint32_t> *chain = NULL);

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
                            std::unordered_set<uint32_t> &prev,
                            std::unordered_set<uint32_t> &curr,
                            uint32_t pos,
                            uint32_t target);

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

      int _abort_after_warns = 10; // Abort after N warnings, -1 means never.

      // Graph construction parameters
      size_t _max_node_len = 10000000; // If a node is longer, split into multiple nodes
  };

  /**
 * @brief
 * Operator used to map Graph::GID's.
 * @param a Graph::GID a
 * @param b Graph::GID b
 */
  inline bool operator<(const Graph::GID &a,
                        const Graph::GID &b) {
      if (a.outgroup != b.outgroup) return a.outgroup < b.outgroup;
      if (a.pct != b.pct) return a.pct < b.pct;
      if (a.num != b.num) return a.num < b.num;
      return a.id < b.id;
  }

  /**
   * @brief
   * Outputs the Graph::GID in a CSV format:\n
   * [o,i],num,id,pct\n
   * @param os Output stream
   * @param gid gid to print
   */
  inline std::ostream &operator<<(std::ostream &os,
                                  const Graph::GID &gid) {
      os << (gid.outgroup ? 'o' : 'i') << ',' << gid.num << ',' << gid.id << ',' << gid.pct;
      return os;
  }

  /**
   * @brief
   * Check if two GID's are equal.
   */
  inline bool operator==(const Graph::GID &a,
                         const Graph::GID &b) {
      if (a.outgroup != b.outgroup) return false;
      if (a.pct != b.pct) return false;
      if (a.num != b.num) return false;
      return a.id == b.id;
  }

}

TEST_CASE ("Node class") {
    vargas::Graph::Node::_newID = 0;
    vargas::Graph::Node n1;
    vargas::Graph::Node n2;
        CHECK(n1.id() == 0);
        CHECK(n2.id() == 1);

        SUBCASE("Node ID change") {
        n1.setID(1);
            CHECK(n1.id() == 0);
        n1.setID(2);
            CHECK(n1.id() == 2);
    }

        SUBCASE("Set Node params") {
        n1.set_seq("ACGTN");
        std::vector<bool> a = {0, 0, 1};
        n1.set_population(a);
        n1.set_endpos(100);

            REQUIRE(n1.seq().size() == 5);

            CHECK(n1.seq()[0] == Base::A);
            CHECK(n1.seq()[1] == Base::C);
            CHECK(n1.seq()[2] == Base::G);
            CHECK(n1.seq()[3] == Base::T);
            CHECK(n1.seq()[4] == Base::N);
            CHECK(n1.end() == 100);
            CHECK(!n1.is_ref());
            CHECK(!n1.belongs(0));
            CHECK(!n1.belongs(1));
            CHECK(n1.belongs(2));

        n1.set_as_ref();
        n1.set_population(3, true);
            CHECK(n1.is_ref());
            CHECK(n1.belongs(0) == true);
            CHECK(n1.belongs(1) == true);
            CHECK(n1.belongs(2) == true);
    }

}

TEST_CASE ("Graph class") {
    vargas::Graph::Node::_newID = 0;
    vargas::Graph g;

    /**
     *     GGG
     *    /   \
     * AAA     TTT
     *    \   /
     *     CCC(ref)
     */

    {
        vargas::Graph::Node n;
        n.set_endpos(3);
        n.set_as_ref();
        std::vector<bool> a = {0, 1, 1};
        n.set_population(a);
        n.set_seq("AAA");
        g.add_node(n);
    }

    {
        vargas::Graph::Node n;
        n.set_endpos(6);
        n.set_as_ref();
        std::vector<bool> a = {0, 0, 1};
        n.set_population(a);
        n.set_af(0.4);
        n.set_seq("CCC");
        g.add_node(n);
    }

    {
        vargas::Graph::Node n;
        n.set_endpos(6);
        n.set_not_ref();
        std::vector<bool> a = {0, 1, 0};
        n.set_population(a);
        n.set_af(0.6);
        n.set_seq("GGG");
        g.add_node(n);
    }

    {
        vargas::Graph::Node n;
        n.set_endpos(9);
        n.set_as_ref();
        std::vector<bool> a = {0, 1, 1};
        n.set_population(a);
        n.set_seq("TTT");
        n.set_af(0.3);
        g.add_node(n);
    }

    g.add_edge(0, 1);
    g.add_edge(0, 2);
    g.add_edge(1, 3);
    g.add_edge(2, 3);

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

        SUBCASE("Proper Graph setup") {
            CHECK(num_to_seq(g.node(0).seq()) == "AAA");
            CHECK(num_to_seq(g.node(1).seq()) == "CCC");
            CHECK(num_to_seq(g.node(2).seq()) == "GGG");
            CHECK(num_to_seq(g.node(3).seq()) == "TTT");
    }

        SUBCASE("Filtering iterators") {
        /**         (ref)
         *     GGG   TTT
         *    /   \ /   \
         * AAA     \     CCA
         *    \   / \   /
         *     CCC   ACA
         *    (ref)
         */
        {
            vargas::Graph::Node n;
            std::vector<bool> pop = {1, 0, 0};
            n.set_population(pop);
            n.set_af(0.7);
            n.set_seq("ACA");
            n.set_endpos(9);
            n.set_not_ref();
            g.add_node(n);
            g.add_edge(1, 4);
            g.add_edge(2, 4);
        }
        {
            vargas::Graph::Node n;
            std::vector<bool> pop = {1, 1, 1};
            n.set_population(pop);
            n.set_af(1);
            n.set_seq("CCA");
            n.set_endpos(12);
            n.set_as_ref();
            g.add_node(n);
            g.add_edge(3, 5);
            g.add_edge(4, 5);
        }

            SUBCASE("Filtering Iterator") {
            vargas::Graph::Population filter(3, false);
            filter.set(2);
            vargas::Graph::GraphIterator i = g.begin(filter);
                CHECK(num_to_seq((*i).seq()) == "AAA");
            ++i;
                CHECK(num_to_seq((*i).seq()) == "CCC");
            ++i;
                CHECK(num_to_seq((*i).seq()) == "TTT");
            ++i;
                CHECK(num_to_seq((*i).seq()) == "CCA");
            ++i;
                CHECK(i == g.end());
        }

            SUBCASE("Filtering Ierator #2") {
            vargas::Graph::Population filter(3, false);
            filter.set(2);
            filter.set(1);
            vargas::Graph::GraphIterator i = g.begin(filter);
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
                CHECK(num_to_seq((*i).seq()) == "CCA");
            ++i;
                CHECK(i == g.end());
        }

            SUBCASE("Filtering Ierator: REF") {
            vargas::Graph::GraphIterator i = g.begin(vargas::Graph::GraphIterator::Type::REF);
                CHECK(num_to_seq((*i).seq()) == "AAA");
            ++i;
                CHECK(num_to_seq((*i).seq()) == "CCC");
            ++i;
                CHECK(num_to_seq((*i).seq()) == "TTT");
            ++i;
                CHECK(num_to_seq((*i).seq()) == "CCA");
            ++i;
                CHECK(i == g.end());
        }

            SUBCASE("Filtering Ierator: MAXAF") {
            vargas::Graph::GraphIterator i = g.begin(vargas::Graph::GraphIterator::Type::MAXAF);
                CHECK(num_to_seq((*i).seq()) == "AAA");
            ++i;
                CHECK(num_to_seq((*i).seq()) == "GGG");
            ++i;
                CHECK(num_to_seq((*i).seq()) == "ACA");
            ++i;
                CHECK(num_to_seq((*i).seq()) == "CCA");
            ++i;
                CHECK(i == g.end());
        }
    }


        SUBCASE("Graph iterator") {
        // Node visit order should be topological
        vargas::Graph::GraphIterator i = g.begin();

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

        SUBCASE("Derived Graph") {
        std::vector<bool> filter = {0, 0, 1};
        vargas::Graph g2(g, filter);

            CHECK(g2.node_map()->size() == 4);
            CHECK(&(*g.node_map()) == &(*g2.node_map())); // Underlying node map unchanged
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

        SUBCASE("REF graph") {
        vargas::Graph g2(g, vargas::Graph::GraphIterator::Type::REF);
        vargas::Graph::GraphIterator iter(g2);

            CHECK((*iter).seq_str() == "AAA");
        ++iter;
            CHECK((*iter).seq_str() == "CCC");
        ++iter;
            CHECK((*iter).seq_str() == "TTT");
        ++iter;
            CHECK(iter == g2.end());
    }

        SUBCASE("MAXAF graph") {
        vargas::Graph g2(g, vargas::Graph::GraphIterator::Type::MAXAF);
        vargas::Graph::GraphIterator iter(g2);

            CHECK((*iter).seq_str() == "AAA");
        ++iter;
            CHECK((*iter).seq_str() == "GGG");
        ++iter;
            CHECK((*iter).seq_str() == "TTT");
        ++iter;
            CHECK(iter == g2.end());

    }

}

TEST_CASE ("Graph Factory") {
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
            << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele count\">" << endl
            << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Num samples at site\">" << endl
            << "##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Num alt alleles\">" << endl
            << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"Length of each alt\">" << endl
            << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"type of variant\">" << endl
            << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2" << endl
            << "x\t9\t.\tG\tA,C,T\t99\t.\tAF=0.01,0.6,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t0|1\t2|3" << endl
            << "x\t10\t.\tC\t<CN7>,<CN0>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
            << "x\t14\t.\tG\t<DUP>,<BLAH>\t99\t.\tAF=0.01,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t1|1" << endl
            << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
            << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
    }

        SUBCASE("File write wrapper") {

            SUBCASE("Basic Graph") {
            vargas::GraphFactory gb(tmpfa);
            gb.open_vcf(tmpvcf);
            gb.node_len(5);
            gb.set_region("x:0-15");

            vargas::Graph g;
            gb.build(g);

            auto giter = g.begin();

                CHECK((*giter).seq_str() == "CAAAT");
                CHECK((*giter).belongs(0) == true); // its a ref
                CHECK((*giter).is_ref());

            ++giter;
                CHECK((*giter).seq_str() == "AAG");
                CHECK((*giter).belongs(0) == true); // its a ref
                CHECK((*giter).is_ref());

            ++giter;
                CHECK((*giter).seq_str() == "G");

            ++giter;
                CHECK((*giter).seq_str() == "A");

            ++giter;
                CHECK((*giter).seq_str() == "C");

            ++giter;
                CHECK((*giter).seq_str() == "T");
                CHECK(!(*giter).is_ref());
                CHECK(!(*giter).belongs(0));
                CHECK(!(*giter).belongs(1));
                CHECK(!(*giter).belongs(2));
                CHECK((*giter).belongs(3));

            ++giter;
                CHECK((*giter).seq_str() == "C");
                CHECK(giter->is_ref());

            ++giter;
                CHECK((*giter).seq_str() == "CCCCC");
                CHECK(!giter->is_ref());
            ++giter;
                CHECK((*giter).seq_str() == "CC");
                CHECK(!giter->is_ref());

            ++giter;
                CHECK(giter->seq_str() == "");

            ++giter;
                CHECK(giter->seq_str() == "TTG");

        }

            SUBCASE("Deriving a Graph") {
            vargas::GraphFactory gb(tmpfa);
            gb.open_vcf(tmpvcf);
            gb.node_len(5);
            gb.set_region("x:0-15");

            vargas::Graph g;
            gb.build(g);

            std::vector<bool> filter = {0, 0, 0, 1};
            vargas::Graph g2(g, filter);
            auto iter = g2.begin();

                CHECK((*iter).seq_str() == "CAAAT");
            ++iter;
                CHECK((*iter).seq_str() == "AAG");
            ++iter;
                CHECK((*iter).seq_str() == "T");
            ++iter;
                CHECK((*iter).seq_str() == "CCCCC");
            ++iter;
                CHECK((*iter).seq_str() == "CC");

        }

            SUBCASE("Sample filtering") {
            vargas::GraphFactory gb(tmpfa);
                CHECK_THROWS(gb.add_sample_filter("-"));
            gb.open_vcf(tmpvcf);
            gb.node_len(5);
            gb.set_region("x:0-15");

                SUBCASE("No inversion") {
                gb.add_sample_filter("s1");
                auto g = gb.build();

                auto giter = g.begin();
                    CHECK((*giter).seq_str() == "CAAAT");
                    CHECK((*giter).belongs(0) == true); // its a ref
                    CHECK((*giter).is_ref());

                ++giter;
                    CHECK((*giter).seq_str() == "AAG");
                    CHECK((*giter).belongs(0) == true); // its a ref
                    CHECK((*giter).is_ref());

                ++giter;
                    CHECK((*giter).seq_str() == "G");

                ++giter;
                    CHECK((*giter).seq_str() == "A");

            }

                SUBCASE("No inversion 2") {
                gb.add_sample_filter("s2");
                auto g = gb.build();

                auto giter = g.begin();
                    CHECK((*giter).seq_str() == "CAAAT");
                    CHECK((*giter).belongs(0) == true); // its a ref
                    CHECK((*giter).is_ref());

                ++giter;
                    CHECK((*giter).seq_str() == "AAG");
                    CHECK((*giter).belongs(0) == true); // its a ref
                    CHECK((*giter).is_ref());

                ++giter;
                    CHECK((*giter).seq_str() == "G");

                ++giter;
                    CHECK((*giter).seq_str() == "C");

                ++giter;
                    CHECK((*giter).seq_str() == "T");
            }

                SUBCASE("Inversion") {
                gb.add_sample_filter("s2", true);
                auto g = gb.build();

                auto giter = g.begin();
                    CHECK((*giter).seq_str() == "CAAAT");
                    CHECK((*giter).belongs(0) == true); // its a ref
                    CHECK((*giter).is_ref());

                ++giter;
                    CHECK((*giter).seq_str() == "AAG");
                    CHECK((*giter).belongs(0) == true); // its a ref
                    CHECK((*giter).is_ref());

                ++giter;
                    CHECK((*giter).seq_str() == "G");

                ++giter;
                    CHECK((*giter).seq_str() == "A");

            }

                SUBCASE("Inversion 2") {
                gb.add_sample_filter("s1", true);
                auto g = gb.build();

                auto giter = g.begin();
                    CHECK((*giter).seq_str() == "CAAAT");
                    CHECK((*giter).belongs(0) == true); // its a ref
                    CHECK((*giter).is_ref());

                ++giter;
                    CHECK((*giter).seq_str() == "AAG");
                    CHECK((*giter).belongs(0) == true); // its a ref
                    CHECK((*giter).is_ref());

                ++giter;
                    CHECK((*giter).seq_str() == "G");

                ++giter;
                    CHECK((*giter).seq_str() == "C");

                ++giter;
                    CHECK((*giter).seq_str() == "T");
            }
        }

    }
    remove(tmpfa.c_str());
    remove(tmpvcf.c_str());
    remove((tmpfa + ".fai").c_str());
}

#endif //VARGAS_GRAPH_H
