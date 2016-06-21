//
// Created by gaddra on 6/17/16.
//

#ifndef VARGAS_SIM_H
#define VARGAS_SIM_H

#define READ_SIM_RATE -1
#define READ_SIM_FIXED -2

#include "readsource.h"
#include "graph.h"

namespace vargas {

  /**
   * Parameter list controlling the types of reads created.
   * @param len Nominal length of the read
   * @param enforce_len reject reads that do not match len
   * @param num_mut Number of mutation errors
   * @param num_indel number of insertions/deletions
   * @param rand Introduce mutations and indels at a random rate
   * @param mut_rate mutation rate, used with rand=true
   * @param indel_rate indel generation rate, used with rand=true
   */
  struct ReadProfile {
      int len = 50;
      bool rand = false;
      float mut = 0;
      float indel = 0;
  };

  inline std::ostream &operator<<(std::ostream &os, const ReadProfile &rp) {
      os << "len=" << rp.len
          << " mut=" << rp.mut
          << " indel=" << rp.indel
          << " rand=" << rp.rand;
      return os;
  }

  /**
   * Generate reads from a graph using a given profile. srand() must be called.
   */
  class ReadSim: public ReadSource {

    public:
      ReadSim(const Graph &g) : _graph(g) { _init(); }

      ReadSim(const Graph &_graph, const ReadProfile &prof) : _graph(_graph), _prof(prof) { _init(); }

      /**
       * generate and store an updated read.
       * @return true if successful
       */
      virtual bool update_read() override {
          auto &nodes = *_graph.node_map();
          auto &next = _graph.next_map();

          size_t curr_node;
          int curr_indiv = -1;
          std::string read_str = "";

          curr_node = next_keys[rand() % next_keys.size()]; // Random start node ID
          size_t curr_pos = rand() % nodes[curr_node]->length(); // current pos relative to node origin

          read.var_bases = 0;
          read.var_nodes = 0;

          while (true) {
              // The first time a branch is hit, pick an individual to track
              if (curr_indiv < 0 && !nodes[curr_node]->is_ref()) {
                  do {
                      curr_indiv = rand() % _graph.pop_size();
                  } while (nodes[curr_node]->individuals()[curr_indiv] == 0);
              }

              // Extract len subseq
              size_t len = _prof.len - read_str.length();
              if (len > nodes.at(curr_node)->length() - curr_pos) len = nodes.at(curr_node)->length() - curr_pos;
              read_str += nodes.at(curr_node)->seq_str().substr(curr_pos, len);
              curr_pos += len;

              if (!nodes[curr_node]->is_ref()) {
                  ++read.var_nodes;
                  read.var_bases += len;
              }
              if (read_str.length() >= _prof.len) break; // Done

              // Pick random next node.
              if (next.find(curr_node) == next.end()) return update_read(); // End of graph
              std::vector<uint32_t> valid_next;
              for (auto n : next.at(curr_node)) {
                  if (curr_indiv < 0 || nodes[n]->belongs(curr_indiv)) valid_next.push_back(n);
              }
              if (valid_next.size() == 0) return update_read();
              curr_node = valid_next[rand() % valid_next.size()];
              curr_pos = 0;
          }

          if (std::count(read_str.begin(), read_str.end(), 'N') >= _prof.len / 2) return update_read();

          // Introduce errors
          read.sub_err = 0;
          read.indel_err = 0;
          std::string read_mut = "";


          // Rate based errors
          if (_prof.rand) {
              for (size_t i = 0; i < read_str.length(); ++i) {
                  char m = read_str[i];
                  // Mutation error
                  if (rand() % 10000 < 10000 * _prof.mut) {
                      do {
                          m = rand_base();
                      } while (m == read_str[i]);
                      ++read.sub_err;
                  }

                  // Insertion
                  if (rand() % 10000 < 5000 * _prof.indel) {
                      read_mut += rand_base();
                      ++read.indel_err;
                  }

                  // Deletion (if we don't enter)
                  if (rand() % 10000 > 5000 * _prof.indel) {
                      read_mut += m;
                      ++read.indel_err;
                  }
              }
          }
          else {
              // Fixed number of errors
              read.sub_err = std::round(_prof.mut);
              read.indel_err = std::round(_prof.indel);
              std::vector<int> indel_pos;
              for (int i = 0; i < read.indel_err;) {
                  int r = rand() % read_str.length();
                  if (std::find(indel_pos.begin(), indel_pos.end(), r) == indel_pos.end()) {
                      indel_pos.push_back(r);
                      ++i;
                  }
              }
              std::sort(indel_pos.begin(), indel_pos.end());
              int prev = 0;
              for (int p : indel_pos) {
                  read_mut += read_str.substr(prev, p);
                  if (rand() % 2) {
                      read_mut += rand_base();
                      read_mut += read_str[p];
                  }
                  prev = p + 1;
              }
              read_mut += read_str.substr(prev, std::string::npos);

              for (int i = 0; i < read.sub_err;) {
                  int r = rand() % read_str.length();
                  if (read_str[r] == read_mut[r]) { // Make sure we don't double mutate same base
                      do {
                          read_mut[r] = rand_base();
                      } while (read_str[r] == read_mut[r]);
                      ++i;
                  }
              }
          }


          read.read_orig = read_str;
          read.read = read_mut;
          read.read_num = seq_to_num(read_mut);
          read.end_pos = nodes[curr_node]->end() - nodes[curr_node]->length() + curr_pos;

          return true;
      }

      /**
       * Get the profile being used to generate reads (if use_prof)
       * @return ReadProfile
       */
      const ReadProfile &prof() const {
          return _prof;
      }

      /**
       * Crete reads following prof as a template
       * @param prof Read Profile
       */
      void set_prof(const ReadProfile &prof) {
          _prof = prof;
      }

      virtual std::string get_header() const override {
          std::stringstream ss;
          ss << _prof;
          return ss.str();
      }

      ReadProfile get_profile() const { return _prof; }


    private:
      const vargas::Graph &_graph;
      std::vector<uint32_t> next_keys;
      ReadProfile _prof;

      void _init() {
          for (auto &n : _graph.next_map()) {
              next_keys.push_back(n.first);
          }
      }

  };


}

TEST_CASE ("Read sim") {
    srand(1);
    vargas::Graph::Node::_newID = 0;
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

    vargas::GraphBuilder gb(tmpfa, tmpvcf);
    gb.node_len(5);
    gb.ingroup(100);
    gb.region("x:0-20");
    vargas::Graph g = gb.build();

        SUBCASE("Mutations") {
        for (int i = 0; i < 3; ++i) {
            vargas::ReadProfile prof;
            prof.len = 5;
            prof.mut = i;
            vargas::ReadSim rs(g, prof);
            rs.update_read();
            auto read = rs.get_read();
            // CHECK(levenshtein_distance(read.read_orig, read.read) <= i + 1);
        }
    }

        SUBCASE("Indels") {
        for (int i = 0; i < 3; ++i) {
            vargas::ReadProfile prof;
            prof.len = 5;
            prof.indel = i;
            vargas::ReadSim rs(g, prof);
            rs.update_read();
            auto read = rs.get_read();
            // CHECK(levenshtein_distance(read.read_orig, read.read) <= i + 1);
        }
    }
}
#endif //VARGAS_SIM_H
