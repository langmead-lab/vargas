/**
 * @author Ravi Gaddipati
 * @date June 24, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Defines a set of subgraphs from a given reference
 * and vcf file.
 *
 * @file
 */

#ifndef VARGAS_GDEF_H
#define VARGAS_GDEF_H

#include <string>
#include <map>
#include <fstream>
#include "graph.h"

namespace Vargas {
  /**
 * @brief
 * Provides an interface for working with GDEF files.
 */
  class GDEF {
    public:
      /**
       * @brief
       * Create a new gdef using the given parameters.
       * @param fasta_file file name of reference seq
       * @param var_file file name of VCF or BCF file
       * @param region region of graph, CHR:XX,XXX-YY,YYY
       * @param node_len maximum node length
       */
      GDEF(std::string fasta_file, std::string var_file, std::string region, int node_len) :
          _fasta_file(fasta_file), _var_file(var_file), _region(region), _node_len(node_len) {
          Vargas::VCF vf(var_file);
          _num_samples = vf.num_samples();
      }

      /**
       * @brief
       * Load an existing GDEF file
       * @param file_name GDEF to load
       */
      GDEF(std::string file_name) { load(file_name); }

      /**
       * @brief
       * Load a given GDEF file.
       * @param file_name GDEF file to load
       */
      void load(std::string file_name) {
          std::ifstream in(file_name);
          if (!in.good()) throw std::invalid_argument("Error opening file \"" + file_name + "\"");

          std::string line;
          std::getline(in, line);
          if (line != "@gdef") throw std::invalid_argument("Expected a GDEF file (" + file_name + ")");

          std::getline(in, line);
          std::vector<std::string> tv_pairs = split(line, GDEF_META_DELIM);
          for (auto &tv : tv_pairs) {
              std::vector<std::string> tv_pair = split(tv, GDEF_TAG_VAL_DELIM);
              if (tv_pair.size() != 2) throw std::invalid_argument("Malformed tag-value pair \"" + tv + "\"");
              auto &tag = tv_pair[0];
              auto &val = tv_pair[1];

              if (tag == REF_FILE_TAG) {
                  _fasta_file = val;
              }
              else if (tag == VAR_FILE_TAG) {
                  _var_file = val;
              }
              else if (tag == REGION_TAG) {
                  _region = val;
              }
              else if (tag == NODE_LEN_TAG) {
                  _node_len = std::stoi(val);
              }
          }

          _pops.clear();
          while (std::getline(in, line)) {
              std::vector<std::string> p_pair = split(line, GDEF_POP_DELIM);
              if (p_pair.size() != 2) throw std::invalid_argument("Invalid Graph::GID:Population formation");
              std::vector<std::string> gid = split(p_pair[0], ',');
              if (gid.size() != 4) throw std::invalid_argument("Invalid Graph::GID format");
              Graph::Population pop(p_pair[1].length());
              for (size_t i = 0; i < p_pair[1].length(); ++i) {
                  if (p_pair[1][i] == '1') pop.set(i);
              }
              Graph::GID g(std::stoi(gid[1]), std::stoi(gid[2]), gid[3] == "1");
              g.outgroup = gid[0] == "o";
              add_population(g, pop);
          }
      }

      /**
       * @brief
       * After params have been set, write the gdef to a file.
       * @param file_name output file name.
       */
      void write(std::string file_name) {
          //Validate args
          {
              GraphBuilder gb(_fasta_file, _var_file);
              gb.region(_region);
              gb.node_len(_node_len);
          }

          std::ofstream out(file_name);
          if (!out.good()) throw std::invalid_argument("Error opening output file\"" + file_name + "\"");
          write(out);
      }

      /**
       * @brief
       * Output to a stream.
       * @param os Output stream.
       */
      void write(std::ostream &os) {
          os << "@gdef\n"
              << REF_FILE_TAG << GDEF_TAG_VAL_DELIM << _fasta_file << GDEF_META_DELIM
              << VAR_FILE_TAG << GDEF_TAG_VAL_DELIM << _var_file << GDEF_META_DELIM
              << REGION_TAG << GDEF_TAG_VAL_DELIM << _region << GDEF_META_DELIM
              << NODE_LEN_TAG << GDEF_TAG_VAL_DELIM << _node_len << std::endl;

          for (auto &tv : _pops) {
              os << tv.first << GDEF_POP_DELIM << tv.second.to_string() << std::endl;
          }
      }

      /**
       * @brief
       * Add a Population to the gdef. All added graphs should
       * be ingroup graphs since outgroups are derived.
       * @param key Graph::GID
       * @param pop Population included in the graph
       */
      bool add_population(const Graph::GID &key, const Graph::Population &pop) {
          if (key.outgroup) throw std::invalid_argument("Only ingroup graphs can be explicitly added.");
          if (_pops.find(key) == _pops.end()) _pops[key] = pop;
          else return false;
          return true;
      }

      /**
       * @brief
       * Creates outgroup versions of all graphs that exist.
       */
      void include_outgroups() {
          auto pop_cpy = _pops;
          for (auto &p : pop_cpy) {
              Graph::GID g = p.first;
              g.outgroup = true;
              _pops[g] = ~p.second;
          }
      }

      /**
       * @brief
       * Set the reference file.
       * @param fasta Relative path to fasta file.
       */
      void set_fasta(std::string fasta) { _fasta_file = fasta; }

      /**
       * @brief
       * Set the variant file.
       * @param varfile Relative path to variant file.
       */
      void set_variant(std::string varfile) {
          _var_file = varfile;
          Vargas::VCF vf(_var_file);
          _num_samples = vf.num_samples();
      }

      /**
       * @brief
       * Set the region of the graph to build.
       * @param region format: CHR:XX,XXX-YY,YYY
       */
      void set_region(std::string region) { _region = region; }

      /**
       * @param node_len maximum graph node length
       */
      void set_node_length(int node_len) { _node_len = node_len; }

      /**
       * @brief
       * Get the number of samples present in the given VCF.
       * @return number of samples, 0 if no var file.
       */
      size_t num_samples() const { return _num_samples; }

      /**
       * @return the specified fasta file.
       */
      std::string fasta() const { return _fasta_file; }

      /**
       * @return the specified variant file.
       */
      std::string var() const { return _var_file; }

      /**
       * @return the specified region.
       */
      std::string region() const { return _region; }

      /**
       * @brief
       * Return a map of all subgraphs defined in the file/added.
       * Outgroup graphs are also included if include_outgroups() was called.
       */
      const std::map<Graph::GID, Graph::Population> &populations() const { return _pops; }


    protected:
      // Used for printing and parsing of files
      const std::string REF_FILE_TAG = "ref";
      const std::string VAR_FILE_TAG = "var";
      const std::string REGION_TAG = "reg";
      const std::string NODE_LEN_TAG = "nlen";
      const char GDEF_META_DELIM = ';';
      const char GDEF_TAG_VAL_DELIM = '=';
      const char GDEF_POP_DELIM = ':';

    private:
      std::string _fasta_file, _var_file;
      std::string _region;
      size_t _node_len = 1000000;
      std::map<Graph::GID, Graph::Population> _pops;
      size_t _num_samples = 0;
  };

}
#endif //VARGAS_GDEF_H
