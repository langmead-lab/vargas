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
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#ifndef VARGAS_GRAPH_H
#define VARGAS_GRAPH_H

#include "fasta.h"
#include "varfile.h"
#include "utils.h"
#include "dyn_bitset.h"

#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>
#include <stdexcept>

namespace vargas {

  using rg::pos_t;

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

      using Population = VCF::Population;

      /**
        * @enum Type
        * Indicate a special type of subgraph, one that only includes REF nodes
        * or nodes with the maximum allele frequency at a given variant set.
        */
      enum class Type {
          REF, /**< Only include reference nodes. */
          MAXAF /**< Keep the node with the highest AF at each branch. */
      };

      /**
       * @brief
       * Represents a node in the directed graphs.
       * @details
       * Sequences are stored numerically.
       * populations are stored as bitsets, where 1 indicates that individual has
       * the given allele.
       */
      class Node {
        public:
          /**
           * @brief
           * Make new node, and assign a unique ID.
           */
          Node() : _id(_newID++) {}

          Node(const Node &n) : _end_pos(n._end_pos), _seq(n._seq), _individuals(n._individuals),
                                _ref(n._ref), _pinch(n._pinch), _af(n._af), _id(n._id) {}

          Node(unsigned pos, const std::string &seq, const Population &pop, bool ref, float af) :
          _end_pos(pos), _seq(rg::seq_to_num(seq)), _individuals(pop), _ref(ref), _af(af), _id(_newID++) {}

          Node &operator=(const Node &n) {
              _end_pos = n._end_pos;
              _seq = n._seq;
              _individuals = n._individuals;
              _ref = n._ref;
              _pinch = n._pinch;
              _af = n._af;
              _id = n._id;
              return *this;
          }

          /**
           * @return length of the sequence
           */
          unsigned length() const { return _seq.size(); }

          /**
           * @return position of last base in seq, 0 indexed
           */
          pos_t end_pos() const { return _end_pos; }

          pos_t begin_pos() const { return _end_pos - _seq.size() + 1; }

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
          const std::vector<rg::Base> &seq() const { return _seq; }

          std::vector<rg::Base> &seq() { return _seq; }

          /**
           * @brief
           * Sequence is stored numerically. Return as a string.
           * @return seq
           */
          std::string seq_str() const { return rg::num_to_seq(_seq); }

          /**
           * @brief
           * Node id.
           * @return unique node ID
           */
          unsigned id() const { return _id; }

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
           * Set the position of the last base in the sequence.
           * @param pos 0-indexed
           */
          void set_endpos(pos_t pos) { this->_end_pos = pos; }

          /**
           * @brief
           * Set the population from an existing Population
           * @param pop
           */
          void set_population(const Population &pop) { _individuals = pop; }

          /**
           * @brief
           * Set the population.
           * @param len number of genotypes
           * @param val true/false for each individual
           */
          void set_population(unsigned len, bool val) { _individuals = Population(len, val); }

          /**
           * @brief
           * Set the stored node sequence. Sequence is converted to numeric form.
           * @param seq
           */
          void set_seq(const std::string &seq) { _seq = rg::seq_to_num(seq); }

          /**
           * @brief
           * Set the stored node sequence
           * @param seq
           */
          void set_seq(const std::vector<rg::Base> &seq) { this->_seq = seq; }

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

          void set_id(unsigned id) {_id = id;}

          void set_pinch(bool p) { _pinch = p;}

          /**
           * @brief
           * Set as pinch node. If this node is removed, then each node before
           * it does not have an edge to a node after this node. All paths traverse
           * through this node.
           */
          void pinch() { _pinch = true; }

          /**
           * @brief
           * If true, then previous alignment seeds can be cleared as no nodes after the current
           * one depends on nodes before current one in a topological-sort ordering.
           * @return true if pinched
           */
          bool is_pinched() const { return _pinch; }

          std::vector<rg::Base>::const_iterator begin() const {
              return _seq.cbegin();
          }

          std::vector<rg::Base>::const_iterator end() const {
              return _seq.cend();
          }

          std::vector<rg::Base>::const_reverse_iterator rbegin() const {
              return _seq.crbegin();
          }

          std::vector<rg::Base>::const_reverse_iterator rend() const {
              return _seq.crend();
          }

          static unsigned _newID; /**< ID of the next instance to be created */

        private:
          pos_t _end_pos; // End position of the sequence
          std::vector<rg::Base> _seq; // sequence in numeric form
          Population _individuals; // Each bit marks an individual, 1 if they have this node
          bool _ref = false; // Part of the reference sequence if true
          bool _pinch = false; // If this node is removed, the graph will split into two distinct subgraphs
          float _af = 1;
          unsigned _id;

      };

      using nodemap_t = std::unordered_map<unsigned, Node>; // Map an ID to a node
      using edgemap_t = std::unordered_map<unsigned, std::vector<unsigned>>; // Map an ID to a vec of next ID's

      /**
       * @brief
       * Forward iterator to traverse the Graph in insertion order. Checks assume underlying graphs are equal.
       */
      template<typename T, bool FWD, typename Unqualified_T = typename std::remove_cv<T>::type>
      class GraphIterator: public std::iterator<std::forward_iterator_tag, Unqualified_T, std::ptrdiff_t, T*, T&> {
        public:

          GraphIterator(const GraphIterator &gi) : _graph(gi._graph), _currID(gi._currID), _empty(0) {}

          /**
           * @param g Graph
           * @param idx Node in the insertion order to begin iterator at.
           */
          GraphIterator(const Graph &g, const unsigned idx = 0) : _graph(g), _currID(idx), _empty(0) {}

          GraphIterator operator=(const GraphIterator &gi) {
              _graph = gi._graph;
              _currID = gi._currID;
          }
          /**
           * @brief
           * Reference to the underlying Graph.
           * @return const ref to graph
           */
          const Graph &graph() const { return _graph; }

          /**
           * @return True if iterators point to same node index.
           */
          bool operator==(const GraphIterator &other) const {
              return _currID == other._currID;
          }

          /**
           * @return true when current Node ID's are not the same or if the
           * underlying graph is not the same.
           */
          bool operator!=(const GraphIterator &other) const {
              return _currID != other._currID;
          }

          /**
           * @brief
           * Goes to next Node.
           * @return iterator to the next Node.
           */
          GraphIterator &operator++() {
              if (FWD) {
                  if (_currID < _graph.get()._add_order.size()) ++_currID;
              }
              else {
                  // reverse iterator
                  const auto s = _graph.get()._add_order.size();
                  if (_currID == 0) _currID = s;
                  else if (_currID != s) --_currID;
              }
              return *this;
          }

          /**
           * @brief
           * Goes to next Node, returns previous node.
           * @return iterator to the next Node.
           */
          GraphIterator operator++(int) {
              auto ret = *this;
              if (FWD) {
                  if (_currID < _graph.get()._add_order.size()) ++_currID;
              }
              else {
                  const auto s = _graph.get()._add_order.size();
                  if (_currID == 0) _currID = s;
                  else if (_currID != s) --_currID;
              }
              return ret;
          }

          /**
           * @brief
           * Const reference to the current node. Undefined for end iterator.
           * @return Node
           */
          T &operator*() const {
              return _graph.get()._IDMap->at(_graph.get()._add_order[_currID]);
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
          const std::vector<unsigned> &incoming() const {
              try {
                  // Assume the previous edge existing is the common case (DAG w 1 start node)
                  return _graph.get()._prev_map.at(_graph.get()._add_order[_currID]);
              } catch (std::exception &e) { return _empty; }
          }

          /**
           * @return vector of all outgoing edges
           */
          const std::vector<unsigned> &outgoing() const {
              try {
                  // Assume the previous edge existing is the common case (DAG w 1 start node)
                  return _graph.get()._next_map.at(_graph.get()._add_order[_currID]);
              } catch (std::exception &e) { return _empty; }
          }

          //TODO icc has a problem with this
          /**
           * @brief
           * Allow conversion from iterator to const_iterator
           */
          operator GraphIterator<const T, FWD>() const {
              return GraphIterator<const T, FWD>(_graph, _currID);
          }

        private:
          std::reference_wrapper<const Graph> _graph;
          unsigned _currID;
          const std::vector<unsigned> _empty;
      };

      using const_iterator = GraphIterator<const Graph::Node, true>;
      using const_reverse_iterator = GraphIterator<const Graph::Node, false>;

      /**
       * @return begin iterator.
       */
      const_iterator begin() const {
          return const_iterator(*this, 0);
      }

      /**
       * @return end iterator.
       */
      const_iterator end() const {
          return const_iterator(*this, _add_order.size());
      }

      const_reverse_iterator rbegin() const {
          return const_reverse_iterator(*this, _add_order.size() - 1);
      }

      const_reverse_iterator rend() const {
          return const_reverse_iterator(*this, _add_order.size());
      }

      /**
       * Seek a sequence position.
       * @param pos 1 based
       * @return pair of iterator and sequence offset
       */
      std::pair<const_iterator, pos_t> seek(pos_t pos) const {
          auto it = begin();
          --pos; // graph pos are stored as 0 indexed
          while (it->end_pos() < pos) ++it;
          return {it, pos - it->begin_pos()};
      };

      /**
       * @brief
       * Default constructor inits a new Graph, including a new node map.
       */
      Graph() : _IDMap(std::make_shared<nodemap_t>()) {}

      Graph(std::shared_ptr<nodemap_t> nodes) : _IDMap(nodes) {}

      /**
       * @brief
       * Create a graph and bind it to an existing node map with given edges.
       */
       Graph(std::shared_ptr<nodemap_t> nodes, const edgemap_t &fwd, const edgemap_t &rev,
             const std::vector<unsigned> &node_order)
       : _IDMap(nodes), _next_map(fwd), _prev_map(rev), _add_order(node_order) {}

      /**
       * @brief
       * Build a graph given a reference FASTA file and a variant VCF or BCF file.
       * @param ref_file FASTA file name
       * @param vcf_file VCF or BCF file name
       * @param region region in the format chromosome:start-end
       */
      Graph(const std::string &ref_file, const std::string &vcf_file, const std::string &region);

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
       * The first node added is set as the Graph root. Nodes must be added in topographical order.
       * @param n node to add, ID of original node is preserved.
       * @return ID of the inserted node
       */
      unsigned add_node(const Node &n);

      /**
       * @brief
       * Create an edge linking two nodes. Previous and Next edges are added.\n
       * n1->n2
       * @param n1 Node one ID
       * @param n2 Node two ID
       */
      bool add_edge(const unsigned n1, const unsigned n2);

      void add_edge_unchecked(const unsigned n1, const unsigned n2);

      /**
       * @brief
       * Maps a ndoe ID to a shared node object
       * @return map of ID, shared_ptr<Node> pairs
       */
      std::shared_ptr<nodemap_t> node_map() const { return _IDMap; }

      /**
       * @brief
       * Maps a node ID to a vector of all next nodes (outgoing edges)
       * @return map of ID, outgoing edge vectors
       */
      const edgemap_t &next_map() const { return _next_map; }

      /**
       * @brief
       *  Maps a node ID to a vector of all incoming edge nodes
       *  @return map of ID, incoming edges
       */
      const edgemap_t &prev_map() const { return _prev_map; }

      /**
       * @brief
       * Const reference to a node
       * @param id ID of the node
       * @return shared node object
       * @throws std::invalid_argument if node ID does not exist
       */
      const Node &node(unsigned id) const {
          if (_IDMap->count(id) == 0) throw std::domain_error("Invalid Node ID.");
          return _IDMap->at(id);
      }

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
      void set_popsize(const unsigned popsize) { _pop_size = popsize; }

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
      unsigned pop_size() const { return _pop_size; }

      /**
       * @brief
       * Return a Population of a subset of the graph.
       * @return Population with ingroup % of samples set.
       */
      Population subset(const int ingroup) const;

      /**
       * @brief
       * Create a subgraph including bases from min to max.
       * @param min
       * @param max
       * @return Graph
       */
      Graph subgraph(const pos_t min, const pos_t max) const;

      /**
       * @brief
       * Ensures that the graph is topologically sorted.
       */
      bool validate() const;

      const std::vector<unsigned> &order() const {
          return _add_order;
      }

      /**
       * Set node order.
       * @param ids Node ID's, topographically ordered
       */
      void set_order(const std::vector<unsigned> &ids) {
          _add_order = ids;
      }

      /**
       * @brief
       * Set forward edges and generate reverse edges.
       * @param edges forward edges
       */
      void set_edges(const edgemap_t &edges) {
          for (auto &p : edges) {
              for (auto to : p.second) {
                  add_edge(p.first, to);
              }
          }
      }

      /**
       * @brief
       * Merges g into self. There should be no node ID conflicts, and g positions should be
       * greater than the maximum position of this.
       * @details
       * if a node ID in g is already present, it will keep the original node.
       * @param g
       */
      void assimilate(const Graph &g) {
          // Insert new nodes
          std::set<unsigned> shared;
          for (auto &p : *g._IDMap) {
              if (!_IDMap->count(p.first)) _IDMap->insert(p);
              else shared.insert(p.first);
          }

          _add_order.reserve(_add_order.size() + g._add_order.size());
          for (auto i : g._add_order) {
              if (!shared.count(i)) _add_order.push_back(i);
          }
          _add_order.shrink_to_fit();

          _merge_edges(_next_map, g._next_map);
          _merge_edges(_prev_map, g._prev_map);
      }

      /**
       * @brief
       * Statistics about the current graph size.
       */
      struct Stats {
          Stats() = default;
          unsigned num_nodes = 0;
          unsigned num_edges = 0;
          unsigned total_length = 0;
          unsigned num_snps = 0;
          unsigned num_dels = 0;

          std::string to_string() const {
              std::stringstream ss;
              ss << "Length: " << total_length << ", Nodes: " << num_nodes << ", Edges: " << num_edges
                 << ", SNPs: " << num_snps << ", Dels: " << num_dels;
              return ss.str();
          }
      };

      /**
       * @return Counted statistics about the current graph.
       */
      Stats statistics() const {
          Stats ret;
          for (auto i : *this) {
              ++ret.num_nodes;
              ret.total_length += i.length();
              ret.num_snps += (i.length() == 1 && !i.is_ref());
              ret.num_dels += (i.length() == 0);
              if (_next_map.count(i.id())) {
                  ret.num_edges += _next_map.at(i.id()).size();
              }
          }
          return ret;
      }


    private:
      // maps a node ID to a nodeptr. Any derived graphs use the same base node ID map.
      std::shared_ptr<nodemap_t> _IDMap;
      // maps a node ID to the vector of nodes it points to
      edgemap_t _next_map;
      // maps a node ID to a vector of node ID's that point to it
      edgemap_t _prev_map;
      std::vector<unsigned> _add_order; // Order nodes were added
      unsigned _pop_size = 0;
      Population _filter;

      /**
       * Given a subset of nodes from Graph g, rebuild all applicable edges in the new graph.
       * @param g underlying parent graph
       * @param includedNodes subset of g's nodes to include
       */
      void _build_derived_edges(const Graph &g, const std::unordered_set<unsigned> &includedNodes);

      /**
       * @brief
       * Inserts new edges from b into a. If a and b have edges from some node id, merge them into a.
       * @param a
       * @param b
       */
      void _merge_edges(edgemap_t &a, const edgemap_t &b) {
          // For each edge pair in b
          for (const auto &p : b) {
              // If we already have that in the map, merge the mapped vectors
              if (a.count(p.first)) {
                  // Insert new edges
                  auto &host = a[p.first];
                  for (auto k : p.second) {
                      if (std::find(host.begin(), host.end(), k) == host.end()) host.push_back(k);
                  }
              } else {
                  a[p.first] = p.second;
              }
          }
      }

  };

  inline Graph &operator+=(Graph &a, const Graph &b) {
      a.assimilate(b);
      return a;
  }

  inline Graph operator+(const Graph &a, const Graph &b) {
      Graph ret(a.node_map(), a.next_map(), a.prev_map(), a.order());
      ret.assimilate(b);
      return ret;
  }

  inline std::ostream &operator<<(std::ostream &os, const vargas::Graph::Stats s) {
      os << s.to_string();
      return os;
  }

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
       * @param ref file
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
       * Set the region of the graph to build. Format should be
       * CHR:XX,XXX-YY,YYY
       * Where CHR is the sequence name, XX,XXX is the min pos and YY,YYY is the max pos.
       * Both are inclusive.
       * @param region
       */
      void set_region(Region region) {
          if (!_vf) throw std::invalid_argument("No variant file opened.");
          _vf->set_region(region);
      }

      /**
       * Only include upto n varaint records. The first n are processed.
       * @param n
       */
      void limit_variants(size_t n) {
          if (!_vf) throw std::invalid_argument("No variant file opened.");
          _vf->limit_num_variants(n);
      }

      /**
       * @brief
       * Assume VCF records for a given contig occur sequentially
       */
      void assume_contig_chr() {
          if (!_vf) throw std::invalid_argument("No variant file opened.");
          _vf->assume_contig_chr();
      }

      /**
       * Open the given file
       * @param file_name
       * @return Number of samples
       */
      unsigned open_vcf(std::string const &file_name);

      /**
       * @brief
       * Limit the samples used from the VCF file.
       * @param filter CSV list of sample names to keep.
       * @param invert Use the samples not specified in filter
       * @return Number of samples in the filter.
       */
      unsigned add_sample_filter(std::string filter, bool invert = false);

      /**
       * Open the given file
       * @param file_name
       * @return Number of samples
       */
      unsigned open_bcf(std::string const &file_name) {
          return open_vcf(file_name);
      }

      /**
       * @brief
       * Apply the various parameters and build the Graph from a VCF and a FASTA.
       * @param g Graph to build into
       */
      void build(Graph &g, pos_t pos_offset=0);

      /**
       * @brief
       * Build the graph using the specified params from a VCF and a FASTA.
       * @return Graph Built Graph
       */
      Graph build(pos_t pos_offset=0) {
          Graph g;
          build(g, pos_offset);
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
       */
      __RG_STRONG_INLINE__
      void _build_edges(Graph &g, std::unordered_set<unsigned> &prev, std::unordered_set<unsigned> &curr);

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
      rg::pos_t _build_linear_ref(Graph &g, std::unordered_set<unsigned> &prev, std::unordered_set<unsigned> &curr,
                                  pos_t pos, pos_t target, pos_t pos_offset);


    private:
      std::string _fa_file;
      std::unique_ptr<VCF> _vf;
      ifasta _fa;

  };

}

#endif //VARGAS_GRAPH_H
