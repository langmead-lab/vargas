/**
 * Ravi Gaddipati
 * Feb 27, 2017
 * rgaddip1@jhu.edu
 *
 * @brief
 * Tools to construct a graph and a set of subgraphs.
 * Graphs are stored in a JSON format (binarized by default)
 *
 * @file
 */

#ifndef VARGAS_GRAPHGEN_H
#define VARGAS_GRAPHGEN_H

#include "json.hpp"
#include "graph.h"
#include "varfile.h"
#include "fasta.h"

namespace vargas {
  using json = nlohmann::json;
  enum class json_mode {cbor, messagepack, plaintext};
  constexpr json_mode default_json_mode = json_mode::cbor;

  // ADL conversions
  inline void to_json(json &j, const Graph::Node &n) {
      j = json{{"id", n.id()},
               {"pos", n.end_pos()},
               {"af", n.freq()},
               {"z", n.is_pinched()},
               {"seq", n.seq_str()}};
  }
  inline void from_json(const json &j, Graph::Node &n) {
      n.set_id(j["id"].get<unsigned>());
      n.set_endpos(j["pos"].get<unsigned>());
      n.set_af(j["af"].get<float>());
      n.set_pinch(j["z"].get<bool>());
      n.set_seq(j["seq"].get<std::string>());
  }

  /**
   * @brief
   * Generate and manage a family of graphs from a FASTA and VCF file.
   * @details
   * GraphGen makes its own structure of nodes. Any provided graphs reference its node map.
   * By default store in a binary JSON format.
   * @code{.json}
   *
    // Though there are comments in this example, comments are not supported.
    {
        "meta": {
            "version" : "vargas VER",
            "fasta" : "Fasta filename",
            "vcf" : "vcf filename",
            "command-line" : "Graph generation CL",
            "samples" : ["sample1", "sample2"]
        },

        "contigs": {
            "chr1" : {
                "min" : 0, // inclusive
                "max" : null, // null is maximum length, inclusive
                "offset" : 0, // Positions in nodes are offset by this amount.
                "nodes": [
                    {"id": 0, "pos": 5, "af" : 1.0, "seq" : [0,1,2,3,2,1]} // numeric rep of seq
                    // ... for all nodes
                ]
            }
            // ... for all chromosomes
        },

        "graphs": {
            // At minimum, base is defined. Keys should correspond to those found in
            // ["meta"]["subgraph-def"]
            "base" : {
                "def" : {
                    "parent" : null, // Derived from parent graph
                    "invert" : false, // Derived from complement of parent
                    "population" : [1,1,1,1], // Population, ength is len ["samples"]
                    "snp" : true, // Only include snps
                    "min-af" : 0, // min allele frequency
                    "max-af" : null,
                    "filter" : "ref OR maxaf"
                },
                "nodes" : [
                    // Nodes from each contig
                    {"contig" : "chr1", "nid" : []} //.. List of node ID nums
                    // ...
                ],
                "edges" : [[0,1],[1,2]] // List of edge pairs
            }
            // ... For all subgraphs
        }
    }
   * @endcode
   */
  class GraphGen {
    public:
      GraphGen() = default;
      GraphGen(const std::string &filename, json_mode mode = default_json_mode) {
          open(filename, mode);
      }

      void write(const std::string &filename, json_mode mode=default_json_mode);
      void open(const std::string &filename, json_mode mode=default_json_mode) {
          json j;
          if (filename.size() == 0) std::cin >> j;
          else {
              std::ifstream in;
              if (mode == json_mode::plaintext) in.open(filename);
              else in.open(filename, std::ios::binary);
              if (!in.good()) throw std::invalid_argument("Error opening file: \"" + filename + "\"");
              in >> j;
          }

          // Validate fields
          if (j.find("contigs") == j.end() || j.find("graphs") == j.end() ||
          !j["contigs"].size() || j["graphs"].find("base") == j["graphs"].end()) {
              throw std::invalid_argument("Invalid graph file.");
          }

          // Load nodes
          _nodes.reset();
          _nodes = std::make_shared<Graph::nodemap_t>();
          for (json::iterator it = j["contigs"].begin(); it != j["contigs"].end(); ++it) {
              assert(it.value().find("nodes") != it.value().end());
              assert(it.value().find("offset") != it.value().end());
              _contig_offsets[it.value()["offset"]] = it.key();
              std::unordered_map<std::string, Graph::Node> i = it.value()["nodes"];
              it.value() = nullptr; // Clear as we load to conserve mem
              std::for_each(i.begin(), i.end(), [this](const std::pair<std::string, Graph::Node> &p) {
                  this->_nodes->insert({std::stoul(p.first), p.second});
              });
          }

          // Load graphs
          _graphs.clear();
          for (json::iterator it = j["graphs"].begin(); it != j["graphs"].end(); ++it) {
              auto &graph = it.value();
              assert(graph.find("nodes") != graph.end());
              assert(graph.find("fwd") != graph.end());
              assert(graph.find("rev") != graph.end());

              std::vector<unsigned> nids;
              for (json::iterator n = graph["nodes"].begin(); n != graph["nodes"].end(); ++n) {
                  std::vector<unsigned> contig_ids = n.value();
                  nids.insert(nids.end(), contig_ids.begin(), contig_ids.end());
              }

              Graph::edgemap_t fwd, rev;
              for (json::iterator i = graph["fwd"].begin(); i != graph["fwd"].end(); ++i) {
                  fwd[std::stoul(i.key())] = i.value().get<Graph::edgemap_t::mapped_type>();
              }
              for (json::iterator i = graph["rev"].begin(); i != graph["rev"].end(); ++i) {
                  rev[std::stoul(i.key())] = i.value().get<Graph::edgemap_t::mapped_type>();
              }

              _graphs[it.key()] = std::make_shared<Graph>(_nodes, fwd, rev);
              graph = nullptr;
          }

      }

      /**
       * @brief
       * Return the contig and position relative to the contig beginning.
       * @param pos offset position
       * @return pair <contig name, position>
       */
      std::pair<std::string, unsigned> absolute_position(unsigned pos) {
          const auto lb = _contig_offsets.lower_bound(pos);
          return {lb->second, pos - lb->first};
      };

    protected:


    private:
      std::shared_ptr<Graph::nodemap_t> _nodes;
      std::unordered_map<std::string, std::shared_ptr<vargas::Graph>> _graphs;
      std::map<unsigned, std::string> _contig_offsets; // Maps an offset to contig
  };
}

TEST_CASE("Load graph") {
    const std::string jfile = "tmpjson.vargas";
    {
        /**
         * AAAAA->GGGGG
         */
        const std::string jstr =
        "{\"contigs\":{\"chr1\":{\"min\":0,\"max\":null,\"offset\":0,\"nodes\":{\"0\":{\"id\":0,\"pos\":5,\"af\":1.0,\"z\":true,\"seq\":\"AAAAA\"},\"1\":{\"id\":1,\"pos\":10,\"af\":1.0,\"z\":true,\"seq\":\"GGGGG\"}}}},\"graphs\":{\"base\":{\"nodes\":{\"chr1\":[0,1]},\"fwd\":{\"0\":[1]},\"rev\":{\"1\":[0]}}}}";
        std::ofstream o(jfile);
        o << jstr;
    }
    SUBCASE("File write wrapper") {
        vargas::GraphGen gg;
        gg.open(jfile);
        gg.absolute_position(13);
    }

    remove(jfile.c_str());
}

#endif //VARGAS_GRAPHGEN_H
