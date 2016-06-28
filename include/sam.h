/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date May 26, 2016
 *
 * @brief
 * Provides a C++ wrapper for htslib handling SAM/BAM files.
 * @details
 * Both file types are handled transparently by htslib. Records
 * are loaded into Reads.
 *
 * @file
 */

#ifndef VARGAS_SAM_H
#define VARGAS_SAM_H

#include "htslib/sam.h"
#include "doctest.h"
#include "readsource.h"

namespace vargas {

  /**
   * @brief
   * Provides an interface for a SAM/BAM File.
   */
  class SAMFile: public ReadSource {
    public:
      SAMFile() { }
      SAMFile(std::string file_name) : _file_name(file_name) { _init(); }
      ~SAMFile() { _deinit(); }

      void open(std::string file_name) {
          _deinit();
          _file_name = file_name;
          _init();
      }

      void close() { _deinit(); }

      bool good() const {
          return _bam && _header;
      }

      std::string get_header() const {
          if (!_header) return "";
          return std::string(_header->text, _header->l_text);
      }

      bool update_read() {
          if (!_bam || !_header) return false;
          if (sam_read1(_bam, _header, _curr_rec) != 0) return false;
          read.end_pos = bam_endpos(_curr_rec) - 1;

          uint8_t *seq = bam_get_seq(_curr_rec);
          read.read = "";
          for (int i = 0; i < _curr_rec->core.l_qseq; ++i) {
              int b = bam_seqi(seq, i);
              switch (b) {
                  case 1:
                      read.read += "A";
                      break;
                  case 2:
                      read.read += "C";
                      break;
                  case 4:
                      read.read += "G";
                      break;
                  case 8:
                      read.read += "T";
                      break;
                  default:
                      read.read += "N";
              }
          }
          read.read_num = seq_to_num(read.read);
          return true;
      }

      std::vector<std::pair<std::string, uint32_t>> ref_names() const {
          std::vector<std::pair<std::string, uint32_t>> ret(0);
          if (!_header) return ret;

          for (int32_t i = 0; i < _header->n_targets; ++i) {
              ret.push_back(std::pair<std::string, uint32_t>(std::string(_header->target_name[i],
                                                                         _header->target_len[i]),
                                                             _header->target_len[i]));
          }
          return ret;
      }

    protected:
      void _init() {
          _bam = sam_open(_file_name.c_str(), "r");
          if (!_bam) throw std::invalid_argument("Error opening file \"" + _file_name + "\"");
          _header = sam_hdr_read(_bam);
          if (!_header) throw std::invalid_argument("Error reading SAM header");
      }

      void _deinit() {
          if (_bam) hts_close(_bam);
          if (_header) bam_hdr_destroy(_header);
          if (_curr_rec) bam_destroy1(_curr_rec);
          _bam = nullptr;
          _header = nullptr;
          _curr_rec = bam_init1();
          _file_name = "";
      }

    private:
      std::string _file_name = "";
      samFile *_bam = nullptr;
      bam_hdr_t *_header = nullptr;
      bam1_t *_curr_rec = bam_init1();
  };

}

TEST_CASE ("SAM File") {
    {
        std::ofstream ss("tmp_s.sam");
        ss << "@HD\tVN:1.0\tSO:coordinate\n";
        ss
            << "@SQ\tSN:1\tLN:249250621\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta\tM5:1b22b98cdeb4a9304cb5d48026a85128\n";
        ss
            << "@SQ\tSN:2\tLN:243199373\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta\tM5:a0d9851da00400dec1098a9255ac712e\n";
        ss
            << "@SQ\tSN:3\tLN:198022430\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta\tM5:fdfd811849cc2fadebc929bb925902e5\n";
        ss
            << "@RG\tID:UM0098:1\tPL:ILLUMINA\tPU:HWUSI-EAS1707-615LHAAXX-L001\tLB:80\tDT:2010-05-05T20:00:00-0400\tSM:SD37743\tCN:UMCORE\n";
        ss
            << "@RG\tID:UM0098:2\tPL:ILLUMINA\tPU:HWUSI-EAS1707-615LHAAXX-L002\tLB:80\tDT:2010-05-05T20:00:00-0400\tSM:SD37743\tCN:UMCORE\n";
        ss << "@PG\tID:bwa\tVN:0.5.4\n";
        ss
            << "1:497:R:-272+13M17D24M\t113\t1\t497\t37\t37M\t15\t100338662\t0\tCGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAG\t0;==-==9;>>>>>=>>>>>>>>>>>=>>>>>>>>>>\tXT:A:U\tNM:i:0\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:37\n";
        ss
            << "19:20389:F:275+18M2D19M\t99\t1\t17644\t0\t37M\t=\t17919\t314\tTATGACTGCTAATAATACCTACACATGTTAGAACCAT\t>>>>>>>>>>>>>>>>>>>><<>>><<>>4::>>:<9\tRG:Z:UM0098:1\tXT:A:R\tNM:i:0\tSM:i:0\tAM:i:0\tX0:i:4\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:37\n";
        ss
            << "19:20389:F:275+18M2D19M\t147\t1\t17919\t0\t18M2D19M\t=\t17644\t-314\tGTAGTACCAACTGTAAGTCCTTATCTTCATACTTTGT\t;44999;499<8<8<<<8<<><<<<><7<;<<<>><<\tXT:A:R\tNM:i:2\tSM:i:0\tAM:i:0\tX0:i:4\tX1:i:0\tXM:i:0\tXO:i:1\tXG:i:2\tMD:Z:18^CA19\n";
        ss
            << "9:21597+10M2I25M:R:-209\t83\t1\t21678\t0\t8M2I27M\t=\t21469\t-244\tCACCACATCACATATACCAAGCCTGGCTGTGTCTTCT\t<;9<<5><<<<><<<>><<><>><9>><>>>9>>><>\tXT:A:R\tNM:i:2\tSM:i:0\tAM:i:0\tX0:i:5\tX1:i:0\tXM:i:0\tXO:i:1\tXG:i:2\tMD:Z:35\n";
    }

    vargas::SAMFile sf("tmp_s.sam");
    sf.update_read();
    int x;
}
#endif //VARGAS_SAM_H
