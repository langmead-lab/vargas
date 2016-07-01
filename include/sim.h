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

#define SIM_SAM_READ_ORIG_TAG "ro"
#define SIM_SAM_INDIV_TAG "nd"
#define SIM_SAM_SUB_ERR_TAG  "se"
#define SIM_SAM_VAR_NODES_TAG "vd"
#define SIM_SAM_VAR_BASE_TAG "vb"
#define SIM_SAM_INDEL_ERR_TAG "ne"
#define SIM_SAM_END_POS_TAG "ep"
#define SIM_SAM_GID_TAG "gd"
#define SIM_SAM_USE_RATE_TAG "rt"
#define SIM_SAM_POPULATION "po"
#define SIM_SAM_REF_TAG "fa"
#define SIM_SAM_VCF_TAG "vf"

#include "graph.h"

namespace Vargas {

  /**
  * @brief
  * Generate reads from a graph using a given profile. srand() should be called externally.
  */
  class Sim {

    public:

      // Tags defining meta information in FASTA read names
      const std::string READ_META_END = "pos";
      const std::string READ_META_MUT = "sub";
      const std::string READ_META_INDEL = "ind";
      const std::string READ_META_VARNODE = "vnd";
      const std::string READ_META_VARBASE = "vbs";
      const std::string READ_META_SRC = "src";
      const char READ_FASTA_META_DELIM = ';';
/**
 * @brief
 * Struct to represent a Read.
 */
      struct Read {
          Read() : read_orig(""), read(""),
                   end_pos(-1), indiv(-1), sub_err(-1), var_nodes(-1), var_bases(-1), indel_err(-1) { }
          Read(std::string r) : read_orig(""), read(r), read_num(seq_to_num(r)),
                                end_pos(-1), indiv(-1), sub_err(-1), var_nodes(-1), var_bases(-1), indel_err(-1) { }

          std::string read_orig;
          /**< unmutated read sequence */
          std::string read;
          /**< base sequence. */
          std::vector<Base> read_num;
          /**< Numeric read representation */
          int32_t end_pos;
          /**< position of last base in seq. */
          int32_t indiv;
          /**< Individual the read was taken from. */
          int32_t sub_err;
          /**< Number of substitiution errors introduced. */
          int32_t var_nodes;
          /**< Number of variant nodes the read traverses. */
          int32_t var_bases;
          /**< Number of bases that are in variant nodes. */
          int32_t indel_err;
          /**< Number of insertions and deletions introduced. */
          Graph::GID src; /**< Read origin graph, as defined in GDEF file. */

      };

      /**
       * @brief
       * Output two lines in FASTA format.
       * @details
       * Output two lines given the form: \n
       * > Meta information \n
       * read_sequence \n
       * @param r Read to print
       * @return two-line string
       */
      inline std::string to_fasta(const Read &r) {
          std::ostringstream ss;
          ss << ">"
              << READ_META_END << ":" << r.end_pos << READ_FASTA_META_DELIM
              << READ_META_MUT << ":" << r.sub_err << READ_FASTA_META_DELIM
              << READ_META_INDEL << ":" << r.indel_err << READ_FASTA_META_DELIM
              << READ_META_VARNODE << ":" << r.var_nodes << READ_FASTA_META_DELIM
              << READ_META_VARBASE << ":" << r.var_bases << READ_FASTA_META_DELIM
              << READ_META_SRC << ":" << r.src
              << std::endl
              << r.read;
          return ss.str();
      }

      /**
       * @brief
       * Convert the read to a single line CSV.
       * @details
       * Output form: \n
       * src,read_seq,end_pos,sub_err,indel_err,var_nodes,var_bases \n
       * @param r Read to print
       * @return single line string
       */
      inline std::string to_csv(const Read &r) {
          std::ostringstream ss;
          ss << r.src << ','
              << r.read << ','
              << r.end_pos << ','
              << r.sub_err << ','
              << r.indel_err << ','
              << r.var_nodes << ','
              << r.var_bases;
          return ss.str();
      }

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

          std::string to_string() const {
              std::ostringstream os;
              os << "len=" << len
                  << ";mut=" << mut
                  << ";indel=" << indel
                  << ";vnode=" << var_nodes
                  << ";vbase=" << var_bases
                  << ";rand=" << rand;
              return os.str();
          }
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
      bool update_read();

      /**
 * @brief
 * Get size reads. If more reads are not available, a undersized
 * batch is returned.
 * @param size nominal number of reads to get.
 */
      const std::vector<Read> &get_batch(int size) {
          if (size <= 0) size = 1;
          _batch.clear();
          for (int i = 0; i < size; ++i) {
              if (!update_read()) break;
              _batch.push_back(_read);
          }
          return _batch;
      }

      /**
       * @brief
       * Get the stored batch of reads.
       * @return vector of Reads
       */
      const std::vector<Read> &batch() const {
          return _batch;
      }

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
      std::string get_header() const {
          return _prof.to_string();
      }

      /**
       * @return current profile being used to filter reads.
       */
      Profile get_profile() const { return _prof; }

      Read &get_read() { return _read; };

    private:
      const Vargas::Graph &_graph;
      std::vector<uint32_t> next_keys;
      std::vector<Read> _batch;
      Profile _prof;
      Read _read;

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
      os << rp.to_string();
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
