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
  enum class json_mode {plaintext}; // binary not well implemented yet
  constexpr json_mode default_json_mode = json_mode::plaintext;

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
      GraphGen(const std::string &filename, json_mode mode = default_json_mode) {
          open(filename, mode);
      }

      std::shared_ptr<Graph> &graph(std::string label) {
          if (_graphs.count(label) == 0) throw std::domain_error("Invalid graph label: " + label);
          return _graphs[label];
      }

      void create_base(const std::string fasta, const std::string vcf="",
                       std::vector<Region> region={}, std::string sample_filter="") {

          if (_nodes == nullptr) _nodes = std::make_shared<Graph::nodemap_t>();
          else _nodes->clear();

          if (region.size() == 0) {
              vargas::ifasta ref(fasta);
              for (const auto &r : ref.sequence_names()) {
                  std::cerr << r << ',';
                  region.emplace_back(r, 0, 0);
              }
          }


          _graphs["base"] = std::make_shared<Graph>(_nodes);
          unsigned offset = 0;
          for (auto reg : region) {
              std::cerr << "Building contig " << reg.seq_name << "\n";
              GraphFactory gf(fasta, vcf);
              gf.add_sample_filter(sample_filter);
              gf.set_region(reg);
              auto g = gf.build(offset);
              std::cerr << g.statistics().to_string() << "\n\n";
              offset = g.rbegin()->end_pos();
              _contig_offsets[offset] = reg.seq_name;
              _graphs["base"]->assimilate(g);
          }

          _j["meta"]["fasta"] = fasta;
          if (vcf.size()) _j["meta"]["vcf"] = vcf;
          sample_filter.erase(std::remove_if(sample_filter.begin(), sample_filter.end(), isspace),
                              sample_filter.end());
          _j["meta"]["samples"] = rg::split(sample_filter, ',');

      }

      void generate_graph(std::string label, const GraphDef &def) {
          if (_graphs.count(def.parent) == 0) throw std::domain_error("Parent graph \"" + def.parent + "\" does not exist.");
      }

      /**
       * @brief
       * Write graphs to a file. Loaded graphs are consumed.
       * @param filename Output file
       * @param mode Output mode
       */
      void write(const std::string &filename, json_mode mode=default_json_mode) {

          // Don't modify nodes if the nodemap is in use
          bool purge = true;
          if (_nodes.use_count() > 1) purge = false;

          // Nodes
          const std::string empty;
          for (auto &p : *_nodes) {
              _j["contigs"][absolute_position(p.second.begin_pos()).first][std::to_string(p.first)] = p.second;
              if (purge) p.second.set_seq(empty); // Free up some memory
          }
          _nodes.reset();

          for (auto &o : _contig_offsets) {
              _j["contigs"][o.second]["offset"] = o.first;
              _j["contigs"][o.second][o.second] = json::array();
          }

          for (auto &g : _graphs) {
              for (auto &p : g.second->order()) {
                  _j["contigs"][absolute_position(g.second->node(p).begin_pos()).first].push_back(p);
              }
              for (auto &f : g.second->next_map()) {
                  _j["graphs"][g.first]["fwd"][std::to_string(f.first)] = f.second;
              }
              for (auto &p : g.second->prev_map()) {
                  _j["graphs"][g.first]["rev"][std::to_string(p.first)] = p.second;
              }
              g.second.reset();
          }
          _graphs.clear();

          if (filename.size() == 0) std::cout << _j;
          else {
              std::ofstream out;
              if (mode == json_mode::plaintext) out.open(filename);
              else std::domain_error("Unimplemented");
              if (!out.good()) throw std::invalid_argument("Error opening file: \"" + filename + "\"");
              out << _j.dump(-1);
          }
      }

      void open(const std::string &filename, json_mode mode=default_json_mode) {
          _j.clear();
          if (filename.size() == 0) {
              std::cin >> _j;
          }
          else {
              std::ifstream in;
              if (mode == json_mode::plaintext) in.open(filename);
              else throw std::domain_error("Unimplemented");
              if (!in.good()) throw std::invalid_argument("Error opening file: \"" + filename + "\"");
              in >> _j;
          }

          // Validate fields
          if (_j.find("contigs") == _j.end() || _j.find("graphs") == _j.end() || !_j["contigs"].size()) {
              throw std::invalid_argument("Invalid graph file.");
          }

          // Load nodes
          _nodes.reset();
          _nodes = std::make_shared<Graph::nodemap_t>();
          for (json::iterator it = _j["contigs"].begin(); it != _j["contigs"].end(); ++it) {
              auto &contig = it.value();
              _contig_offsets[contig["offset"].get<unsigned>()] = it.key();
              for (json::iterator n = contig["nodes"].begin(); n != contig["nodes"].end(); ++n) {
                  Graph::Node node = n.value().get<Graph::Node>();
                  _nodes->insert({node.id(), node});
                  n.value() = nullptr;
              }
              contig["nodes"] = nullptr;
          }

          // Load graphs
          _graphs.clear();
          for (json::iterator it = _j["graphs"].begin(); it != _j["graphs"].end(); ++it) {
              auto &graph = it.value();

              std::vector<unsigned> nids;
              for (json::iterator n = graph["nodes"].begin(); n != graph["nodes"].end(); ++n) {
                  std::vector<unsigned> contig_ids = n.value();
                  nids.insert(nids.end(), contig_ids.begin(), contig_ids.end());
              }

              Graph::edgemap_t fwd, rev;
              for (json::iterator i = graph["fwd"].begin(); i != graph["fwd"].end(); ++i) {
                  fwd[std::stoul(i.key())] = i.value().get<Graph::edgemap_t::mapped_type>();
              }
              graph["fwd"] = nullptr;
              for (json::iterator i = graph["rev"].begin(); i != graph["rev"].end(); ++i) {
                  rev[std::stoul(i.key())] = i.value().get<Graph::edgemap_t::mapped_type>();
              }
              graph["rev"] = nullptr;

              _graphs[it.key()] = std::make_shared<Graph>(_nodes, std::move(fwd), std::move(rev), std::move(nids));
          }

      }

      /**
       * Clear JSON data.
       */
      void dump() {
          _j.clear();
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

      json &get_json() {
          return _j;
      }

    private:
      std::shared_ptr<Graph::nodemap_t> _nodes;
      std::unordered_map<std::string, std::shared_ptr<vargas::Graph>> _graphs; // Map label to a graph
      std::map<unsigned, std::string> _contig_offsets; // Maps an offset to contig
      json _j; // Maintains all info except for nodes, graph edges for space
  };
}

TEST_CASE("Load graph") {
    const std::string jfile = "tmpjson.vargas";
    /**
     * AAAAA->GGGGG
     */
    const std::string jstr =
    "{\"contigs\":{\"chr1\":{\"min\":0,\"max\":null,\"offset\":0,\"nodes\":{\"0\":{\"id\":0,\"pos\":5,\"af\":1.0,"
    "\"z\":true,\"seq\":\"AAAAA\"},\"1\":{\"id\":1,\"pos\":10,\"af\":1.0,\"z\":true,\"seq\":\"GGGGG\"}}}},\"graphs\":"
    "{\"base\":{\"nodes\":{\"chr1\":[0,1]},\"fwd\":{\"0\":[1]},\"rev\":{\"1\":[0]}}}}";
    std::ofstream o(jfile);
    o << jstr;
    o.close();

    SUBCASE("File write wrapper") {
        vargas::GraphGen gg;
        gg.open(jfile);

    }

    remove(jfile.c_str());
}

TEST_CASE("Write graph") {
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
        << "y\t34\t.\tC\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t39\t.\tC\tT,G\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
    }

    SUBCASE("Concat graphs") {
        std::vector<vargas::Region> regions;
        regions.reserve(24);
        for (unsigned i = 1; i < 23; ++i) {
            regions.emplace_back(std::to_string(i), 0, 0);
        }
        regions.emplace_back("X", 0, 0);
        regions.emplace_back("Y", 0, 0);

        {
            vargas::GraphGen gg;
            gg.create_base("GRCh38_renamed.fa", "", regions);
            gg.graph("base")->to_DOT("g.dot", "g");
            gg.write("graph.json");
        }
        {
            vargas::GraphGen gg;
            gg.open("graph.json");
            gg.graph("base")->to_DOT("g2.dot", "g");
            gg.write("graph2.json");
        }
    }

    remove(tmpfa.c_str());
    remove(tmpvcf.c_str());
}

#endif //VARGAS_GRAPHGEN_H
