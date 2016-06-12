/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date May 28, 2016
 *
 * Implementation of a directed graph. Each node stores a sequence and relevant
 * information. Graphs can be derived from other graphs with a filter, allowing
 * the extraction of population subsets.
 *
 * @file graph.h
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
#include "doctest/doctest.h"
#include "utils.h"
#include "dyn_bitset.h"

namespace vargas {

/**
 * Represents a Graph of the genome. The Graph is backed by a map of Graph::Nodes, and edges
 * are backed by a map of node ID's.
 */
  class Graph {

    public:

      typedef VarFile::Population Population;

      /**
       * When a normal population filter is not used, a flag can be used. REF includes
       * only reference alleles, MAXAF picks the allele with the highest frequency.
       * Both result in linear graphs.
       * Filter used as a placeholder, it should never be passed as a param.
       */
      enum Type { TOPO, REF, MAXAF, FILTER, END };

      /**
       * Represents a node in the directed graphs. Sequences are stored numerically.
       * populations are stored as bitsets, where 1 indicates that indivudal posseses
       * the allele.
       */
      class Node {
        public:
          // Assign a unique ID to each node
          Node() : _id(_newID++) { }
          Node(int pos,
               const std::string &seq,
               const std::vector<bool> &pop,
               bool ref,
               float af) :
              _endPos(pos), _seq(seq_to_num(seq)), _individuals(pop), _ref(ref), _af(af), _id(_newID++) { }
          Node(int pos,
               const std::string &seq,
               const Population &pop,
               bool ref,
               float af) :
              _endPos(pos), _seq(seq_to_num(seq)), _individuals(pop), _ref(ref), _af(af), _id(_newID++) { }

          /**
           * @return length of the sequence
           */
          size_t length() const { return _seq.size(); } // Length of sequence

          /**
           * @return position of last base in seq, 0 indexed
           */
          int end() const { return _endPos; } // Sequence end position in genome

          /**
           * Return true if the node is a ref node, or if the bit is set.
           * @param idx bit index
           * @return belongs
           */
          bool belongs(uint idx) const {
              if (_ref) return true;
              return _individuals.at(idx);
          }

          /**
           * Returns true if the node is a ref node or if the population intersects the node population.
           * @param pop Population filter
           * @return belongs
           */
          bool belongs(const Population &pop) {
              if (_ref) return true;
              return pop && _individuals;
          }

          /**
           * Sequence as a vector of unsigned chars.
           * @return seq
           */
          const std::vector<Base> &seq() const { return _seq; }

          /**
           * Sequence is stored numerically. Return as a string.
           * @return seq
           */
          std::string seq_str() const { return num_to_seq(_seq); }

          /**
           * Size of the Population of the node. This should be consistant throughout the graph.
           * @return pop_size
           */
          ulong pop_size() const { return _individuals.size(); }

          /**
           * Node id.
           * @return unique node ID
           */
          uint32_t id() const { return _id; }

          /**
           * True if node common to all individuals, or REF.
           * @return is_ref
           */
          bool is_ref() const { return _ref; }

          /**
           * Allele frequency.
           * @return af
           */
          float freq() const { return _af; }

          /**
           * Reference to the raw Population data member.
           * @return individuals
           */
          const Population &individuals() const { return _individuals; }

          static uint32_t _newID; // ID of the next instance to be created

          /**
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
           * Set the position of the last base in the sequence.
           * @param pos 0-indexed
           */
          void set_endpos(int pos) { this->_endPos = pos; }

          /**
           * Set the population from an existing Population
           * @param pop
           */
          void set_population(const Population &pop) { _individuals = pop; }

          /**
           * Set the population from a vector. Each individual is set if pop[i] evaluates true.
           * @param pop
           */
          template<typename T>
          void set_population(const std::vector<T> &pop) { _individuals = pop; }

          /**
           * Set the population.
           * @param len number of genotypes
           * @param val true/false for each individual
           */
          void set_population(size_t len, bool val) { _individuals = Population(len, val); }

          /**
           * Set the stored node sequence. Sequence is converted to numeric form.
           * @param seq
           */
          void set_seq(std::string seq) { _seq = seq_to_num(seq); }

          /**
           * Set the stored node sequence
           * @param seq
           */
          void set_seq(std::vector<Base> &seq) { this->_seq = seq; }

          /**
           * Sets the node as a reference node, and sets all bits in the population.
           */
          void set_as_ref() {
              _ref = true;
              _individuals.set();
          }

          /**
           * Deselects node as a referene node, does not modify population.
           */
          void set_not_ref() { _ref = false; }

          /**
           * Set the allele frequency of the node. Used for MAXAF filtering.
           * @param af float frequency, between 0 and 1
           */
          void set_af(float af) { _af = af; }

        private:
          int _endPos; // End position of the sequence
          std::vector<Base> _seq; // sequence in numeric form
          Population _individuals; // Each bit marks an individual, 1 if they have this node
          bool _ref = false; // Part of the reference sequence if true
          float _af = 1;
          uint32_t _id;

      };

      typedef std::shared_ptr<Node> nodeptr;

      /**
       * Default constructor inits a new Graph, including a new node map.
       */
      Graph() : _IDMap(std::make_shared<std::unordered_map<uint32_t, nodeptr>>(std::unordered_map<uint32_t,
                                                                                                  nodeptr>())) { }

      /**
       * Create a Graph with another Graph and a population filter. The new Graph will only
       * contain nodes if any of the individuals in filter possess the node. The actual nodes
       * are shared_ptr's to the parent Graph, as to prevent duplication of Nodes.
       * This constructor is slow for large graphs, when possible a filtering iterator should be
       * used instead.
       * @param g Graph to derive the new Graph from
       * @param filter population filter, only include nodes representative of this population
       */
      Graph(const Graph &g,
            const Population &filter);

      /**
       * Constructs a graph using filter tags, one of:
       * Graph::REF, or Graph::MAXAF
       * The former keeps reference nodes, the later picks the node with the highest allele
       * frequency. Both result in linear graphs.
       * This constructor is slow for large graphs, when possible a filtering iterator should be
       * used instead.
       * @param g Graph to derive from
       * @aram t one of Graph::REF, Graph::MAXAF
       */
      Graph(const Graph &g,
            Type t);

      /**
       * Add a new node to the Graph. A new node is created so the original can be destroyed.
       * The first node added is set as the Graph root. Nodes must be added in topographical order.
       * @param n node to add, ID is preserved.
       */
      uint32_t add_node(Node &n);

      /**
       * Create an edge linking two nodes. Previous and Next edges are added.
       * n1->n2
       * @param n1 Node one ID
       * @param n2 Node two ID
       */
      bool add_edge(uint32_t n1,
                    uint32_t n2);

      /**
       * Sets the root of the Graph.
       * @param id ID of root node
       */
      void set_root(uint32_t id) {
          _root = id;
      }

      /**
       * Set the graph description.
       * @param description
       */
      void set_desc(std::string description) { _desc = description; }

      /**
       * Return root node ID
       * @return root
       */
      uint32_t root() const { return _root; }

      /**
       * Maps a ndoe ID to a shared node object
       * @return map of ID, shared_ptr<Node> pairs
       */
      const std::shared_ptr<std::unordered_map<uint32_t, nodeptr>> &node_map() const { return _IDMap; }

      /**
       * Maps a node ID to a vector of all next nodes (outgoing edges)
       * @return map of ID, outgoing edge vectors
       */
      const std::unordered_map<uint32_t, std::vector<uint32_t>> &next_map() const { return _next_map; }

      /**
       *  Maps a node ID to a vector of all incoming edge nodes
       *  @return map of ID, incoming edges
       */
      const std::unordered_map<uint32_t, std::vector<uint32_t>> &prev_map() const { return _prev_map; }

      /**
       * Const reference to a node
       * @param ID of the node
       * @return shared node object
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
       * Exports the graph in DOT format.
       * @param graph name
       */
      std::string to_DOT(std::string name = "g") const;

      /**
       * Export the graph in DOT format
       * @param filename to export to
       * @param name of the graph
       */
      void to_DOT(std::string filename, std::string name) {
          std::ofstream out(filename);
          if (!out.good()) throw std::invalid_argument("Error opening file: \"" + filename + "\"");
          out << to_DOT(name);
      }

      /**
       * Define the graph population size.
       * @param popsize, number of genotypes
       */
      void set_popsize(int popsize) { _pop_size = popsize; }

      /**
       * @return population size
       */
      size_t pop_size() const { return _pop_size; }

      /**
       * Return a Population of a subset of the graph.
       * @return Population with ingroup % indivduals set.
       */
      Population subset(int ingroup) const {
          Population p(_pop_size);
          for (size_t i = 0; i < _pop_size; ++i) {
              if (rand() % 100 < ingroup) p.set(i);
          }
          return p;
      }

      /**
       * const forward iterator to traverse the graph while applying a filter.
       * When incrementing the iterator, only nodes that match a condition will
       * return.
       * Options include:
       * Filtering: Provided a population, a node is returned if there is an intersection
       * REF: Only return reference nodes
       * MAXAF: Return the node with the highest allele frequency.
       */
      class FilteringIter {
        public:

          /**
           * Accept all nodes.
           * @ param g Graph
           */
          explicit FilteringIter(const Graph &g) :
              _graph(g), _type(TOPO), _currID(0) { }

          /**
           * Traverse nodes when filter has at least one individual in common with the
           * node population.
           * @param g graph
           * @param filter Population the node is checked against. If there is an intersection,
           * the node is returned.
           */
          explicit FilteringIter(const Graph &g, const Population &filter) :
              _graph(g), _filter(filter), _type(FILTER), _currID(g._root) { }

          /**
           * Traverse a linear subgraph.
           * @param g graph
           * @param type one of Graph::MAXAF, Graph::REF
           */
          explicit FilteringIter(const Graph &g, Graph::Type type)
              : _graph(g), _type(type), _currID(g._root) {
              if (_type == TOPO) _currID = 0;
          }

          /**
           * Reference to the underlying graph.
           * @return const ref to graph
           */
          const Graph &graph() const { return _graph; }

          /**
           * @return true when underlying graph address is the same and current node ID's
           * are the same. Two end iterators always compare equal.
           */
          bool operator==(const FilteringIter &other) const {
              if (_type == END && other._type == END) return true; // All ends are equal
              if (_type != other._type) return false; // Same type of iterator
              if (&_graph != &other._graph) return false; // same base graph
              return _currID == other._currID;
          }

          /**
           * @return true when current node ID's are not the same or if the
           * underlying graph is not the same. Two end iterators always compare
           * false.
           */
          bool operator!=(const FilteringIter &other) const {
              if (_type == END && other._type == END) return false;
              if (_type != other._type) return true;
              if (&_graph != &other._graph) return true;
              return _currID != other._currID;
          }

          /**
           * Goes to next node. Nodes are only included if it satisfies the filter.
           * Once the end of the graph is reached, _end is set.
           * @return iterator to the next node.
           */
          FilteringIter &operator++() {
              // If end of graph has been reached
              if (_type == END) return *this;

              if (_type == TOPO) {
                  ++_currID;
                  if (_currID == _add_order_size) _type = END;
                  return *this;
              }

              if (_graph._next_map.count(_currID) == 0) {
                  _type = END;
                  return *this;
              }

              const auto &next_vec = _graph._next_map.at(_currID);
              const auto &graph_map = *(_graph._IDMap);

              switch (_type) {
                  case REF:
                      for (uint32_t nextID : next_vec) {
                          if (graph_map.at(nextID)->is_ref()) {
                              _insert_queue(nextID);
                              break; // Assuming there is only one REF node per branch
                          }
                      }
                      break;

                  case FILTER:
                      for (uint32_t nextID : next_vec) {
                          // Add all nodes that intersect with filter
                          if (graph_map.at(nextID)->belongs(_filter)) _insert_queue(nextID);
                      }
                      break;

                  case MAXAF: {
                      uint32_t max_id = graph_map.at(next_vec.at(0))->id();
                      float max_af = graph_map.at(next_vec.at(0))->freq();
                      float freq;
                      for (size_t i = 1; i < next_vec.size(); ++i) {
                          freq = graph_map.at(next_vec.at(i))->freq();
                          if (freq > max_af) {
                              max_af = freq;
                              max_id = graph_map.at(next_vec.at(i))->id();
                          }
                      }
                      _insert_queue(max_id);
                  }
                      break;

                  default:
                      throw std::logic_error("Invalid type.");
                      break;
              }


              if (_queue.empty()) {
                  _type = END;
                  return *this;
              }

              _traversed.insert(_currID); // Insert previous iterator position (node id)
              _currID = _queue.front();
              _queue.pop();
              _queue_unique.erase(_currID);
              return *this;
          }

          Type type() const { return _type; }

          /**
           * Inserts the ID into the queue if it's unique.
           * @param id node id to insert
           */
          __attribute__((always_inline))
          inline void _insert_queue(uint32_t id) {
              if (_queue_unique.count(id) == 0) {
                  _queue_unique.insert(id);
                  _queue.push(id);
              }
          }

          /**
           * Const reference to the current node.
           * @return Node
           */
          const Graph::Node &operator*() const {
              if (_type == TOPO) return *(_graph._IDMap->at(_graph._add_order.at(_currID)));
              if (_type == END) return *(_graph._IDMap->at(*(_graph._add_order.end())));
              return *(_graph._IDMap->at(_currID));
          }

          /**
           * All nodes that we've traversed that have incoming edges to the current node.
           * @return vector of previous nodes
           */
          const std::vector<uint32_t> &incoming() {
              _incoming.clear();
              if (_graph._prev_map.count(_currID) == 0) return _incoming;
              if (_type == TOPO) return _graph._prev_map.at(_graph._add_order.at(_currID));
              for (auto &id : _graph._prev_map.at(_currID)) {
                  if (_traversed.count(id)) _incoming.push_back(id);
              }
              return _incoming;
          }

          const std::vector<uint32_t> &outgoing() {
              if (_graph._prev_map.count(_currID) == 0) return _outgoing;
              return _graph._next_map.at(_graph._add_order.at(_currID));
          }


        private:
          const Graph &_graph; // Underlying graph
          Population _filter; // Nodes that intersect with _filter are included if _type == FILTER
          Graph::Type _type; // Set to MAXAF or REF for linear subgraph traversals
          uint32_t _currID;
          std::queue<uint32_t> _queue;
          std::unordered_set<uint32_t> _queue_unique; // Used to make sure we don't repeat nodes
          std::unordered_set<uint32_t> _traversed; // Set of all nodes we've passed
          std::vector<uint32_t> _incoming; // Set of incoming edges, for use by incoming()
          const std::vector<uint32_t> _outgoing; // empty vec
          size_t _add_order_size = _graph._add_order.size();

      };

      /**
       * Provides an iterator to the whole graph.
       * @return reference to the root node.
       */
      FilteringIter begin() const {
          return FilteringIter(*this);
      }


      /**
       * Iterator when conditions are applied.
       * @param filter population to compare nodes to
       */
      FilteringIter begin(const Population &filter) const {
          return FilteringIter(*this, filter);
      }

      /**
       * Linear subgraph interator.
       * @param type one of Graph::REF, Graph::MAXAF
       */
      FilteringIter begin(Graph::Type type) const {
          return FilteringIter(*this, type);
      }

      /**
       * @return end iterator.
       */
      FilteringIter end() const {
          return FilteringIter(*this, END);
      }

      void node_len(size_t len) { _max_node_len = len; }
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

      /**
       * Given a subset of nodes from Graph g, rebuild all applicable edges in the new graph.
       * @param g underlying parent graph
       * @param includedNodes subset of g's nodes to include
       */
      void _build_derived_edges(const Graph &g,
                                const std::unordered_map<uint32_t, nodeptr> &includedNodes);

  };

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

      /**   GGG
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
          g.to_DOT("tmp.dot", "g");

              SUBCASE("Filtering Iterator") {
              Graph::Population filter(3, false);
              filter.set(2);
              vargas::Graph::FilteringIter i = g.begin(filter);
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
              Graph::Population filter(3, false);
              filter.set(2);
              filter.set(1);
              vargas::Graph::FilteringIter i = g.begin(filter);
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
              vargas::Graph::FilteringIter i = g.begin(Graph::REF);
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
              vargas::Graph::FilteringIter i = g.begin(Graph::MAXAF);
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
          vargas::Graph::FilteringIter i = g.begin();

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
          vargas::Graph g2(g, vargas::Graph::REF);
          vargas::Graph::FilteringIter iter(g2);

              CHECK((*iter).seq_str() == "AAA");
          ++iter;
              CHECK((*iter).seq_str() == "CCC");
          ++iter;
              CHECK((*iter).seq_str() == "TTT");
          ++iter;
              CHECK(iter == g2.end());
      }

          SUBCASE("MAXAF graph") {
          vargas::Graph g2(g, vargas::Graph::MAXAF);
          vargas::Graph::FilteringIter iter(g2);

              CHECK((*iter).seq_str() == "AAA");
          ++iter;
              CHECK((*iter).seq_str() == "GGG");
          ++iter;
              CHECK((*iter).seq_str() == "TTT");
          ++iter;
              CHECK(iter == g2.end());

      }

  }

  class GraphBuilder {

    public:
      GraphBuilder(std::string reffile,
                   std::string vcffile) :
          _fa_file(reffile), _vf_file(vcffile) { }

      void open(std::string ref,
                std::string vcf) {
          _fa_file = ref;
          _vf_file = vcf;
      }

      bool good() const {
          return !_fa.good() || !_vf.good();
      }

      void region(std::string region) {
          _vf.set_region(region);
      }

      void region(std::string chr,
                  int min,
                  int max) {
          _vf.set_region(chr, min, max);
      }

      /**
       * Use a certain percentage of individuals. Reference nodes are always included.
       * @param percent, 0 - 100
       */
      void ingroup(int percent);

      /**
       * Set maximum node length. If <= 0, length is unbounded.
       * @param max maximum node length.
       */
      void node_len(int max) { _max_node_len = max; }

      /**
       * Apply the various parameters and build the Graph.
       * @return pointer to Graph.
       */
      void build(Graph &g);

    protected:
      __attribute__((always_inline))
      inline void _build_edges(Graph &g, std::vector<uint32_t> &prev,
                               std::vector<uint32_t> &curr);

      __attribute__((always_inline))
      inline int _build_linear_ref(Graph &g, std::vector<uint32_t> &prev,
                                   std::vector<uint32_t> &curr,
                                   uint32_t pos,
                                   uint32_t target);

      __attribute__((always_inline))
      inline std::vector<std::string> _split_seq(std::string seq);

    private:
      std::string _fa_file, _vf_file;
      VarFile _vf;
      FASTAFile _fa;
      Graph g;

      // Graph construction parameters
      int _ingroup = 100; // percent of individuals to use. Ref nodes always included
      int _max_node_len = 10000000;
  };

}

TEST_CASE ("Graph Builder") {
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
            << "x\t10\t.\tC\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
            << "x\t14\t.\tG\t<DUP>,<BLAH>\t99\t.\tAF=0.01,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t1|1" << endl
            << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
            << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
    }

        SUBCASE("Basic Graph") {
        vargas::GraphBuilder gb(tmpfa, tmpvcf);
        gb.node_len(5);
        gb.ingroup(100);
        gb.region("x:0-15");

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

    }

        SUBCASE("Deriving a Graph") {
        vargas::GraphBuilder gb(tmpfa, tmpvcf);
        gb.node_len(5);
        gb.ingroup(100);
        gb.region("x:0-15");

        vargas::Graph g;
        gb.build(g);

        std::vector<bool> filter = {0, 0, 0, 1};
        vargas::Graph g2(g, filter);

        g2.to_DOT("tmp.dot", "gr");

        vargas::Graph::FilteringIter iter(g2);

            CHECK((*iter).seq_str() == "CAAAT");
        ++iter;
            CHECK((*iter).seq_str() == "AAG");
        ++iter;
            CHECK((*iter).seq_str() == "G");
        ++iter;
            CHECK((*iter).seq_str() == "T");
        ++iter;
            CHECK((*iter).seq_str() == "C");
        ++iter;
            CHECK((*iter).seq_str() == "CC");

    }


    remove(tmpfa.c_str());
    remove(tmpvcf.c_str());
    remove((tmpfa + ".fai").c_str());
}

#endif //VARGAS_GRAPH_H
