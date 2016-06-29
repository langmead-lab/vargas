/**
 * @author Ravi Gaddipati
 * @date June 26, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Simulates random reads from a graph, returning reads that follow a specified Sim::Profile.
 *
 * @file
 */

#ifndef VARGAS_SIM_H
#define VARGAS_SIM_H

#include "readsource.h"
#include "graph.h"

namespace Vargas {

  /**
  * @brief
  * Generate reads from a graph using a given profile. srand() should be called externally.
  */
  class Sim: public ReadSource {

    public:

      /**
       * @brief
       * Parameter list controlling the types of reads created.
       */
      struct Profile {
          unsigned int len = 50;
          /**< Nominal length of the read */
          bool rand = false;
          /**< Number of mutation errors, or rate */
          float mut = 0;
          /**< number of insertions/deletions, or rate */
          float indel = 0;
          /**< Introduce mutations and indels at a random rate */
          int var_nodes = -1;
          /**< Number of variant nodes */
          int var_bases = -1; /**< number of total variant bases */
      };

      /**
       * @param g Graph to simulate from
       */
      Sim(const Graph &g) : _graph(g) { _init(); }

      /**
       * @param _graph Graph to simulate from
       * @param prof accept reads following this profile
       */
      Sim(const Graph &_graph, const Profile &prof) : _graph(_graph), _prof(prof) { _init(); }

      /**
       * @brief
       * Generate and store an updated read.
       * @return true if successful
       */
      virtual bool update_read() override;

      /**
       * @brief
       * Get the profile being used to generate reads (if use_prof)
       * @return Sim::Profile
       */
      const Profile &prof() const {
          return _prof;
      }

      /**
       * @brief
       * Crete reads following prof as a template
       * @param prof Read Profile
       */
      void set_prof(const Profile &prof) {
          _prof = prof;
      }

      /**
       * @return the profile used to generate the reads
       */
      virtual std::string get_header() const override {
          std::ostringstream ss;
          return ss.str();
      }

      /**
       * @return current profile being used to filter reads.
       */
      Profile get_profile() const { return _prof; }


    private:
      const Vargas::Graph &_graph;
      std::vector<uint32_t> next_keys;
      Profile _prof;

      /**
       * Creates a vector of keys of all outgoing edges. Allows for random node selection
       * This precludes the possibility of having reads begin in the last node of the graph.
       */
      void _init() {
          for (auto &n : _graph.next_map()) {
              next_keys.push_back(n.first);
          }
      }

  };

  /**
   * @brief
   * Output the profile.
   * @details
   * Form: \n
   * len=READ_LEN mut=NUM_MUT indel=NUM_INDEL vnode=NUM_VAR_NODES vbase=NUM_VAR_BASE rand=USE_RATES \n
   * @param os output stream
   * @param rp Read Profile
   * @return output stream
   */
  inline std::ostream &operator<<(std::ostream &os, const Sim::Profile &rp) {
      os << "len=" << rp.len
          << " mut=" << rp.mut
          << " indel=" << rp.indel
          << " vnode=" << rp.var_nodes
          << " vbase=" << rp.var_bases
          << " rand=" << rp.rand;
      return os;
  }

}

TEST_CASE ("Read sim") {
    srand(1);
    Vargas::Graph::Node::_newID = 0;
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

    Vargas::GraphBuilder gb(tmpfa, tmpvcf);
    gb.node_len(5);
    gb.ingroup(100);
    gb.region("x:0-20");
    Vargas::Graph g = gb.build();

        SUBCASE("Mutations") {
        for (int i = 0; i < 3; ++i) {
            Vargas::Sim::Profile prof;
            prof.len = 5;
            prof.mut = i;
            Vargas::Sim rs(g, prof);
            rs.update_read();
            auto read = rs.get_read();
            // CHECK(levenshtein_distance(read.read_orig, read.read) <= i + 1);
        }
    }

        SUBCASE("Indels") {
        for (int i = 0; i < 3; ++i) {
            Vargas::Sim::Profile prof;
            prof.len = 5;
            prof.indel = i;
            Vargas::Sim rs(g, prof);
            rs.update_read();
            auto read = rs.get_read();
            // CHECK(levenshtein_distance(read.read_orig, read.read) <= i + 1);
        }
    }
}

#endif //VARGAS_SIM_H
