/**
 * Ravi Gaddipati
 * Feb 27, 2017
 * rgaddip1@jhu.edu
 *
 * @brief
 * Tools to construct a graph and a set of subgraphs.
 * Graphs are stored in a JSON format (binarized by default)
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#ifndef VARGAS_GRAPHGEN_H
#define VARGAS_GRAPHGEN_H

#include "graph.h"
#include "varfile.h"
#include "fasta.h"
#include "dyn_bitset.h"

#include <stdexcept>
#include <random>
#include <chrono>


namespace vargas {

  /**
   * @brief
   * Resolve graph coords to a contig and position.
   */
  struct coordinate_resolver {
      /**
       * @brief
       * Return the contig and position relative to the contig beginning.
       * @param pos offset position, 1 indexed
       * @return pair <contig name, position>
       */
      std::pair<std::string, unsigned> resolve(unsigned pos) const {
          if (_contig_offsets.size() == 0) return {"", pos};
          std::map<unsigned, std::string>::const_iterator lb = _contig_offsets.lower_bound(pos);
          if (lb != _contig_offsets.begin()) --lb; // For rare case that pos = 0
          return {lb->second, pos - lb->first};
      }

      std::map<unsigned, std::string> _contig_offsets; // Maps an offset to contig
  };

  /*
   * @brief
   * Handle graph generation and file IO
   *
   * @details
   * Graph definition file format (tab delimited):\n
   *
   * @code{.txt}
   * @vgraph
   * date <graph build date>
   * fasta <reference fasta name>
   * vargas-build <vargas build date>
   * vcf <vcf file name>
   * [other meta info]
   *
   * @contigs
   * <offset> <contig name>
   * ...
   *
   * @graphs
   * <name> <node id list> <edges>
   * ...
   *
   * @nodes
   * <ID> <endpos> <frequency> <pinched> <ref> <seqsize>
   * <node sequence>
   * ...
   *
   * @endcode
   */
  class GraphMan {
    public:
      GraphMan() = default;
      GraphMan(const std::string &filename) {
          open(filename);
      }

      /**
       * @brief
       * Create a base graph. This replaces any current graphs. Also creates REF and MAXAF graphs.
       * @param fasta filename
       * @param vcf filename
       * @param region list of regions to include. Default all.
       * @param sample_filter Use only these samples. Default all ("")
       * @return pointer to new base graph
       */
      std::shared_ptr<Graph>
      create_base(const std::string fasta, const std::string vcf="", std::vector<Region> region={},
                  std::string sample_filter="", bool print=false);


      /**
       * @brief
       * Write graphs to a file. Graphs are consumed while written.
       * @param filename Output file
       */
      void write(const std::string &filename);

      /**
       * @brief
       * Open a graph definition file.
       * @param filename
       */
      void open(const std::string &filename);

      /**
       * @brief
       * Return the contig and position relative to the contig beginning.
       * @param pos offset position, 1 indexed
       * @return pair <contig name, position>
       */
      std::pair<std::string, unsigned> absolute_position(unsigned pos) const {
          return _resolver.resolve(pos);
      };

      /**
       * @return Object to resolve coordinates with.
       */
      coordinate_resolver resolver() const {
          return _resolver;
      }

      size_t count(std::string label) const {
          return _graphs.count(label);
      }

      std::shared_ptr<Graph> at(std::string label) const {
          std::transform(label.begin(), label.end(), label.begin(), tolower);
          if (!count(label)) throw std::domain_error("No graph named \"" + label + "\"");
          return _graphs.at(label);
      }

      std::shared_ptr<Graph> operator[](std::string label) {
          std::transform(label.begin(), label.end(), label.begin(), tolower);
          return _graphs[label];
      }

      /**
       * @return Available graph labels.
       */
      std::vector<std::string> labels() const {
          std::vector<std::string> ret;
          for (const auto &i : _graphs) ret.push_back(i.first);
          return ret;
      }

      /**
       * @brief
       * Parse a subgraph definition and create the child graph. This should be called after building the base.
       * Note that samples are only loaded from original VCF files, and are not persisted. You cannot derive
       * a graph from a previous GDEF. Labels are not case sensitive.
       * @details
       * Format : [<parent>:]*<label>=[0-9]{1,2}[%]
       * Implied root parent is the base graph, or all the samples.
       * @return label of the graph
       */
      std::string derive(std::string def);


    private:
      std::shared_ptr<Graph::nodemap_t> _nodes;
      std::map<std::string, std::shared_ptr<vargas::Graph>> _graphs; // Map label to a graph
      coordinate_resolver _resolver;
      std::map<std::string, std::string> _aux;
  };
}

#endif //VARGAS_GRAPHGEN_H
