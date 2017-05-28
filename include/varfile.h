/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date May 26, 2016
 *
 * @brief
 * Provides a C++ wrapper for htslib handling of VCF and BCF files.
 * @details
 * Both file types are handled transparently by htslib. The records
 * are parsed to substitute in copy number variations, and skip
 * records outside of a defined range. A subset of individuals can be
 * defined using create_ingroup.
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#ifndef VARGAS_VARFILE_H
#define VARGAS_VARFILE_H

#include "dyn_bitset.h"
#include "utils.h"
#include "htslib/vcfutils.h"
#include "htslib/hts.h"

#include <string>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <set>
#include <unordered_map>
#include <map>

namespace vargas {

  using rg::pos_t;

  struct Region;

  /**
     * @brief
     * Parse a region in the format:
     * SEQ_NAME:MIN_POS-MAX_POS
     * Commas are stripped.
     * @param region_str
     * @return Region packet.
     */
  Region parse_region(const std::string &region_str);

  /**
   * @brief
   * Packet to represent a parsed region string.
   */
  struct Region {
      Region() : min(0), max(0) {}
      Region(const std::string &seq, unsigned min, unsigned max) : seq_name(seq), min(min), max(max) {}
      Region(std::string str) {
          auto r = parse_region(str);
          seq_name = r.seq_name;
          min = r.min;
          max = r.max;
      }
      std::string seq_name; /**< Contig name */
      pos_t min, /**< Min position, inclusive. */
      max; /**< Max position, inclusive. */
  };

/**
 * @brief
 * Provides an interface to a VCF/BCF file. Core processing
 * provided by htslib. \n
 * Usage: \n
 * @code{.cpp}
 * #include "varfile.h"
 *
 * Vargas::VCF vcf("variants.bcf");
 * std::vector<std::string> seqs = vcf.sequences(); // Names of CHROM's present
 *
 * vcf.set_region("22:0-10,000,000");
 * vcf.create_ingroup(50); // Randomly select 50% of individuals to include
 *
 * vcf.samples().size(); // 4 individuals
 *
 * // Allele with the maximum allele frequency for each VCF record
 * while(vcf.next()) {
 *  float max_af = 0;
 *  std::string max_allele;
 *  std::vector<float> &freqs = vcf.frequencies();
 *  for (size_t i = 0; i < freqs.size(); ++i) {
 *    if (f > max_af) {
 *     max_af = f;
 *     max_allele = vcf.alleles[i];
 *     }
 *   }
 *   std::cout << "Variant position:" << vcf.pos() << ", Ref:" << vcf.ref()
 *             << ", Max AF allele:" << max_allele << ", Freq:" << max_af
 *             << ", Population:" << vcf.allele_pop(max_allele).to_string();
 *
 *   // 4 samples/individuals = 8 genotypes
 *   // Variant position: 100 Ref: A, Max AF allele: G, Freq: 0.6, Population: 01100101
 * }
 *
 * @endcode
 */
  class VCF {
    public:

      typedef dyn_bitset<64> Population;

      VCF() {}

      /**
       * @param file VCF/BCF File name
       */
      VCF(std::string file) : _file_name(file) {
          _init();
      }

      /**
       * @param file VCF/BCF File name
       * @param chr Chromosome
       * @param min Min position, 0 indexed
       * @param max Max position, 0 indexed, inclusive
       */
      VCF(std::string file, std::string chr, pos_t min, pos_t max) : _file_name(file), _region(chr, min, max) {
          _init();
      }

      ~VCF() {
          close();
      }

      /**
       * @brief
      * Get the specified format field from the record.
      * @tparam T Valid types are int, char, float.
      */
      template<typename T>
      class FormatField {
        public:
          /**
           * Get the specified field tag.
           * @param hdr VCF Header
           * @param rec current record
           * @param tag Field to get, e.g. "GT"
           */
          FormatField(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag) : tag(tag) {
              if (!hdr || !rec || tag.length() == 0) throw std::invalid_argument("Invalid header, rec, or tag.");

              T *dst = nullptr;
              int n_arr = 0;

              int n = _get_vals(hdr, rec, tag, &dst, n_arr);

              /*
              if (n == -1) throw std::invalid_argument("No such tag in header: " + tag);
              else if (n == -2) throw std::invalid_argument("Header and tag type clash: " + tag);
              else if (n == -3) throw std::invalid_argument(tag + " does not exist in record.");
               */

              for (int i = 0; i < n; ++i) {
                  values.push_back(dst[i]);
              }

              if (dst) free(dst); // get_format_values allocates
          }

          std::vector<T> values; /**< Retrieved values. */
          std::string tag; /**< Type of FORMAT or INFO field. */

        private:
          // Change the parse type based on what kind of type we have
          inline int _get_vals(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, int32_t **dst, int &ndst) {
              return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_INT);
          }

          inline int _get_vals(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, float **dst, int &ndst) {
              return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_REAL);
          }

          inline int _get_vals(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, char **dst, int &ndst) {
              return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_STR);
          }
      };

      /**
       * @brief
       * Get the specified info field from the record.
       * @tparam T Valid types are int, char, float.
       */
      template<typename T>
      class InfoField {
        public:
          /**
           * @brief
           * Get the specified field.
           * @param hdr VCF Header
           * @param rec current record
           * @param tag Field to get, e.g. "GT"
           */
          InfoField(bcf_hdr_t *hdr,
                    bcf1_t *rec,
                    std::string tag) : tag(tag) {
              if (!hdr || !rec || tag.length() == 0) throw std::invalid_argument("Invalid header, rec, or tag.");

              T *dst = nullptr;
              int n_arr = 0;

              int n = _bcf_get_info_values(hdr, rec, tag, &dst, n_arr);
              /*
              if (n == -1) throw std::invalid_argument("No such tag in header: " + tag);
              else if (n == -2) throw std::invalid_argument("Header and tag type clash: " + tag);
              else if (n == -3) throw std::invalid_argument(tag + " does not exist in record.");
              */
              for (int i = 0; i < n; ++i) {
                  values.push_back(dst[i]);
              }

              if (dst) free(dst);
          }

          std::vector<T> values;
          /**< Retrieved values */
          std::string tag; /**< Type of FORMAT or INFO field. */

        private:
          // Change the parse type based on what kind of type we have
          inline int _bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, int32_t **dst, int &ndst) {
              return bcf_get_info_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_INT);
          }

          inline int _bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, float **dst, int &ndst) {
              return bcf_get_info_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_REAL);
          }

          inline int _bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, char **dst, int &ndst) {
              return bcf_get_info_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_STR);
          }
      };

      /**
       * @brief
       * Open the specified VCF or BCF file and load the header. If empty, load an "empty" vcf.
       * @param file filename
       * @return -1 on file open error, -2 on header load error, 0 otherwise
       */
      int open(std::string file) {
          _file_name = file;
          return _init();
      }

      void close();

      bool good() {
          return _header && _bcf;
      }

      /**
       * @brief
       * Ingroup parameter on BCF reading. Empty string indicates none, "-" indicates all.
       * @return string of ingroup samples
       */
      std::string ingroup_str() const {
          if (!_ingroup_cstr) return "-";
          return std::string(_ingroup_cstr);
      }

      /**
       * @brief
       * Get a list of sequences in the VCF file.
       * @return vector of sequence names
       */
      std::vector<std::string> seq_names() const;

      /**
       * @brief
       * counts each haplotype as distinct. num_haplotypes() = samples().size() * 2
       * @return Number of samples the VCF has. Each sample represents two genotypes.
       */
      size_t num_haplotypes() const;

      /**
       * @brief
       * Get a vector of sample names.
       * @return vector of samples
       */
      const std::vector<std::string> &samples() const {
          return _samples;
      }

      /**
       * @brief
       * Load the next VCF record. All information is unpacked,
       * subject to sample set restrictions.
       * @return false on read error or if outside restriction range.
       */
      bool next();

      /**
       * @brief
       * Unpacks shared information as well as all sample information.
       * Subject to sample set restrictions.
       */
      void unpack_all() {
          bcf_unpack(_curr_rec, BCF_UN_ALL);
          _load_shared();
      }

      /**
       * @brief
       * Reference allele of the current record.
       * @return reference allele
       */
      std::string ref() const {
          return _alleles[0];
      }

      /**
       * @brief
       * List of all the alleles in the current record.
       * @details
       * The first is the reference. Allele copy number variant tags are converted,
       * whereas other tags are substituted for the reference.
       * @return vector of alleles
       */
      const std::vector<std::string> &alleles() const {
          return _alleles;
      }

      /**
       * @brief
       * 0 based position, i.e. the VCF pos - 1.
       * @return position.
       */
      pos_t pos() const {
          return _curr_rec->pos;
      }

      /**
       * @brief
       * Get a list of alleles for all samples (subject to sample set restriction).
       * @details
       * Consecutive alleles represent phasing, e.g. all odd indexes are one phase,
       * all even indexes are the other. Call will unpack the full record.
       * Explicit copy number variations are replaced, other ambiguous types are replaced.
       * A map is also built, mapping each allele to the subpopulation that has it.
       * @return Vector of alleles, ordered by sample.
       */
      const std::vector<std::string> &gen_genotypes();

      /**
       * @brief
       * Get the allele frequencies of the ref and alt alleles.
       * The ref freq is computed with 1-sum(alt_frequencies).
       * @return const ref to vector of frequencies
       */
      const std::vector<float> &frequencies() const;

      /**
       * @brief
       * Get values of an arbitrary INFO tag.
       * @tparam T Field format
       * @param tag tag to retrieve values for
       * @return vector of values.
       */
      template<typename T>
      std::vector<T> info_tag(std::string tag) {
          return InfoField<T>(_header, _curr_rec, tag).values;
      }

      /**
       * @brief
      * Get values of an arbitrary FORMAT tag.
      * @tparam T Field format
       * @param tag tag to extract
      * @return vector of values.
      */
      template<typename T>
      std::vector<T> format_tag(std::string tag) {
          return FormatField<T>(_header, _curr_rec, tag).values;
      }

      /**
       * @brief
       * Return the population set that has the allele.
       * @details
       * The returned vector has the same size as number of genotypes (samples * 2).
       * When true, that individual/phase has that allele.
       * @param allele allele to get the population of
       * @return Population of individuals that have the allele
       */
      const Population allele_pop(const std::string allele) const {
          if (_genotype_indivs.count(allele)) return _genotype_indivs.at(allele);
          else return Population(0);
      }

      /**
       * @brief
       * Check if the file is properly loaded.
       * @return true if file is open and has a valid header.
       */
      bool good() const {
          return _header && _bcf;
      }

      /**
       * @brief
       * File name of VCF/BCF file.
       * @return file name
       */
      std::string file() const {
          return _file_name;
      }


      /**
       * @brief
       * Create a random subset of the samples.
       * @param percent of samples to keep.
       */
      void create_ingroup(int percent);

      /**
       * @brief
       * Include only the provided sample names in the Graph.
       * @param samples vector of sample names
       */
      void create_ingroup(const std::vector<std::string> &samples) {
          _ingroup = samples;
          _apply_ingroup_filter();
      }

      /**
       * @return string passed to htslib for ingroup filtering.
       */
      const std::vector<std::string> &ingroup() const {
          return _ingroup;
      }

      void set_region(const Region &region);
      const Region &region() const { return _region; }

      /**
       * Limit the number of variants to first num records.
       * @param num
       */
      void limit_num_variants(size_t num) {
          _limit = num;
      }

      /**
       * @brief
       * Assume that all variants for a given contig will appear consecutively in the file.
       */
      void assume_contig_chr() {
          _assume_contig = true;
      }

    protected:

      /**
       * @brief
       * Open the provided file and load the header.
       * @return -1 on file open error, -2 on header read error.
       */
      int _init();

      /**
       * @brief
       * Loads the list of alleles for the current record.
       */
      void _load_shared();

      /**
       * @brief
       * Applies the contents of _ingroup to the header. The filter
       * impacts all following unpacks.
       */
      void _apply_ingroup_filter();


    private:
      std::string _file_name; // VCF/BCF file name
      Region _region;

      htsFile *_bcf = nullptr;
      bcf_hdr_t *_header = nullptr;
      bcf1_t *_curr_rec = bcf_init();

      std::vector<std::string> _genotypes; // restricted to _ingroup
      std::unordered_map<std::string, Population> _genotype_indivs;
      std::vector<std::string> _alleles;
      std::vector<std::string> _samples;
      std::vector<std::string> _ingroup; // subset of _samples
      char *_ingroup_cstr = nullptr;

      size_t _limit, _counter;

      bool _assume_contig, _entered_contig;

  };

  inline std::ostream &operator<<(std::ostream &os, const VCF &vcf) {
      os << "POS: " << vcf.pos() << " REF: " << vcf.ref() << " ALTS: ";
      for (const auto &al : vcf.alleles()) os << al << ' ';
      return os;
  }

}

#endif //VARGAS_VARFILE_H
