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

#include "graph.h"
#include "varfile.h"
#include "fasta.h"
#include "dyn_bitset.h"


namespace vargas {

  /**
   * @brief
   * Defintion of a subgraph
   * @param parent label of parent graph
   * @param invert use the complement of the parents population
   * @param snp_only Only consider SNP variants
   * @param Use a subset of the parent
   * @param linear Graph is linear
   * @param type if linear, is it REF or MAXAF
   * @param population if not linear, population subset
   * @param min_af if not linear, minimum allele frequency
   * @param max_af if not linear, max allele frequency
   */
  struct GraphDef {
      std::string parent;
      bool invert;
      bool snp_only;
      std::vector<Region> region;
      bool linear;
      VCF::Population population;
      Graph::Type type;
      float min_af, max_af;
  };

  class GraphGen {
    public:
      GraphGen() = default;
      GraphGen(const std::string &filename, bool build=true) {
          open(filename, build);
      }

      /**
       * @brief
       * Create a base graph. This replaces any current graphs.
       * @param fasta filename
       * @param vcf filename
       * @param region list of regions to include. Default all.
       * @param sample_filter Use only these samples. Default all ("")
       * @return pointer to new base graph
       */
      std::shared_ptr<Graph> create_base(const std::string fasta, const std::string vcf="",
                                         std::vector<Region> region={}, std::string sample_filter="",
                                         bool print=false);

      std::shared_ptr<Graph> generate_subgraph(std::string label, const GraphDef &def);

      /**
       * @brief
       * Write graphs to a file. Graphs are consumed while written.
       * @details
       * Meta information is maintained in JSON object. To write
       * nodes, edges, and contigs are converted.
       * @param filename Output file
       * @param mode Output mode
       */
      void write(const std::string &filename);

      void open(const std::string &filename, bool build=true);



      /**
       * @brief
       * Return the contig and position relative to the contig beginning.
       * @param pos offset position, 1 indexed
       * @return pair <contig name, position>
       */
      std::pair<std::string, unsigned> absolute_position(unsigned pos) const;

      size_t count(std::string label) const {
          return _graphs.count(label);
      }

      std::shared_ptr<Graph> at(std::string label) {
          if (!count(label)) throw std::domain_error("No graph named \"" + label + "\"");
          return _graphs[label];
      }

      std::shared_ptr<Graph> operator[](std::string label) {
          return _graphs[label];
      }

      std::vector<std::string> labels() const {
          std::vector<std::string> ret;
          for (const auto &i : _graphs) ret.push_back(i.first);
          return ret;
      }


    private:
      std::shared_ptr<Graph::nodemap_t> _nodes;
      std::map<std::string, std::shared_ptr<vargas::Graph>> _graphs; // Map label to a graph
      std::map<std::string, GraphDef> _graph_def; // Map label to graph def
      std::map<unsigned, std::string> _contig_offsets; // Maps an offset to contig
      std::map<std::string, std::string> _aux;
  };
}

#endif //VARGAS_GRAPHGEN_H
