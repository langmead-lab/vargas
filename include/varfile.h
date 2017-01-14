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
 * @file
 */

#ifndef VARGAS_VARFILE_H
#define VARGAS_VARFILE_H

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
#include "dyn_bitset.h"
#include "utils.h"
#include "htslib/vcfutils.h"
#include "htslib/hts.h"

namespace vargas {

  /**
   * @brief
   * Packet to represent a parsed region string.
   */
  struct Region {
      Region() : min(0), max(0) {}
      Region(const std::string &seq, unsigned min, unsigned max) : seq_name(seq), min(min), max(max) {}
      std::string seq_name; /**< Contig name */
      unsigned min, /**< Min position, inclusive. */
      max; /**< Max position, inclusive. */
  };

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
   * Base class for files representing variants in a reference. Provides an interface
   * to iterate through variants, and encapsulates a reference region.
   */
  class VariantFile {
    public:
      VariantFile() {}

      /**
       * @brief
       * Specify a region of the reference sequence.
       * @param chr Chromosome
       * @param min Minimum position, 0 indexed
       * @param max Max position, inclusive, 0 indexed
       */
      VariantFile(std::string const &chr, unsigned min, unsigned max) : _region(chr, min, max) {}

      virtual ~VariantFile() = default;

      typedef dyn_bitset<64> Population;

      /**
       * @brief
        * Set the minimum and maximum position, inclusive.
        * If max is == 0, go until end.
        * @param chr contig ID
        * @param min Minimum position, 0 indexed
        * @param max Max position, inclusive, 0 indexed
     */
      void set_region(std::string chr, unsigned min, unsigned max);

      /**
     * @brief
     * Parse a region string in the format: \n
     * CHR:XX,XXX-YY,YYY \n
     * commas are stripped, range is inclusive. 0 indexed.
     * If max is <= 0, go until end.
     * @param region region string
     */
      void set_region(const std::string &region);
      void set_region(const Region &region);

      const Region &region() const { return _region; }

      virtual void create_ingroup(const std::vector<std::string> &) {
          throw std::invalid_argument("No default impl.");
      };

      /**
       * @brief
       * Load the next VCF record. All information is unpacked,
       * subject to sample set restrictions.
       * Upon initilization, the first record is loaded.
       * @return false on read error or if outside restriction range.
       */
      virtual bool next() = 0;

      /**
       * @return true if the opened file is good.
       */
      virtual bool good() = 0;

      /**
       * @return Current reference sequence
       */
      virtual std::string ref() const = 0;

      /**
       * @return Vector of alleles at the current position
       */
      virtual const std::vector<std::string> &alleles() const = 0;

      /**
       * @return Current position of variant
       */
      virtual int pos() const = 0;

      /**
       * @brief
       * Allele frequencies corresponding to alleles()
       * @return vector of freqs, 0-1
       */
      virtual const std::vector<float> &frequencies() const = 0;

      /**
       * @brief
       * Samples in the file. If samples are unknown, an empty vector is returned.
       * @return vector of sample names
       */
      virtual const std::vector<std::string> &samples() const = 0;

      /**
       * @brief
       * Number of samples in the file.
       * @return
       */
      virtual size_t num_samples() const = 0;

      /**
       * @brief
       * Population that posseses the the given allele. If the variant source does not
       * have this information, {1, 1} is returned.
       * @return Population of the given allele at the current position.
       */
      virtual const Population &allele_pop(const std::string allele = "") const = 0;


    protected:
      Region _region;
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
  class VCF: public VariantFile {
    public:

      VCF() {}

      /**
       * @param file VCF/BCF File name
       */
      explicit VCF(std::string file) : _file_name(file) {
          _init();
      }

      /**
       * @param file VCF/BCF File name
       * @param chr Chromosome
       * @param min Min position, 0 indexed
       * @param max Max position, 0 indexed, inclusive
       */
      VCF(std::string file, std::string chr, int min, int max) : VariantFile(chr, min, max), _file_name(file) {
          _init();
      }

      ~VCF() override {
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
          FormatField(bcf_hdr_t *hdr,
                      bcf1_t *rec,
                      std::string tag) : tag(tag) {
              if (!hdr || !rec || tag.length() == 0) throw std::invalid_argument("Invalid header, rec, or tag.");

              T *dst = nullptr;
              int n_arr = 0;

              int n = _get_vals(hdr, rec, tag, &dst, n_arr);
              if (n == -1) throw std::invalid_argument("No such tag in header: " + tag);
              else if (n == -2) throw std::invalid_argument("Header and tag type clash: " + tag);
              else if (n == -3) throw std::invalid_argument(tag + " does not exist in record.");

              for (int i = 0; i < n; ++i) {
                  values.push_back(dst[i]);
              }

              free(dst); // get_format_values allocates
          }

          std::vector<T> values;
          /**< Retrieved values. */
          std::string tag; /**< Type of FORMAT or INFO field. */

        private:
          // Change the parse type based on what kind of type we have
          inline int _get_vals(bcf_hdr_t *hdr,
                               bcf1_t *rec,
                               std::string tag,
                               int32_t **dst,
                               int &ndst) {
              return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_INT);
          }

          inline int _get_vals(bcf_hdr_t *hdr,
                               bcf1_t *rec,
                               std::string tag,
                               float **dst,
                               int &ndst) {
              return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_REAL);
          }

          inline int _get_vals(bcf_hdr_t *hdr,
                               bcf1_t *rec,
                               std::string tag,
                               char **dst,
                               int &ndst) {
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
              if (n == -1) throw std::invalid_argument("No such tag in header: " + tag);
              else if (n == -2) throw std::invalid_argument("Header and tag type clash: " + tag);
              else if (n == -3) throw std::invalid_argument(tag + " does not exist in record.");

              for (int i = 0; i < n; ++i) {
                  values.push_back(dst[i]);
              }

              free(dst);
          }

          std::vector<T> values;
          /**< Retrieved values */
          std::string tag; /**< Type of FORMAT or INFO field. */

        private:
          // Change the parse type based on what kind of type we have
          inline int _bcf_get_info_values(bcf_hdr_t *hdr,
                                          bcf1_t *rec,
                                          std::string tag,
                                          int32_t **dst,
                                          int &ndst) {
              return bcf_get_info_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_INT);
          }

          inline int _bcf_get_info_values(bcf_hdr_t *hdr,
                                          bcf1_t *rec,
                                          std::string tag,
                                          float **dst,
                                          int &ndst) {
              return bcf_get_info_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_REAL);
          }

          inline int _bcf_get_info_values(bcf_hdr_t *hdr,
                                          bcf1_t *rec,
                                          std::string tag,
                                          char **dst,
                                          int &ndst) {
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

      bool good() override {
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
      std::vector<std::string> sequences() const;

      /**
       * @brief
       * num_samples() counts each haplotype as distinct. num_samples() = samples().size() * 2
       * @return Number of samples the VCF has. Each sample represents two genotypes.
       */
      size_t num_samples() const override;

      /**
       * @brief
       * Get a vector of sample names.
       * @return vector of samples
       */
      const std::vector<std::string> &samples() const override {
          return _samples;
      }

      /**
       * @brief
       * Load the next VCF record. All information is unpacked,
       * subject to sample set restrictions.
       * @return false on read error or if outside restriction range.
       */
      bool next() override;

      /**
       * @brief
       * Unpack only the shared information, and loads ref and allele info.
       */
      void unpack_shr() {
          bcf_unpack(_curr_rec, BCF_UN_SHR);
          _load_shared();
      }

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
      int pos() const {
          return _curr_rec->pos;
      }

      /**
       * @brief
       * Get a list of alleles for all samples (subject to sample set restriction).
       * @details
       * Consecutive alleles represent phasing, e.g. all odd indexes are one phase,
       * all even indexes are the other. Call will unpack the full record.
       * Explicit copy number variations are replaced, other ambigous types are replaced
       * ambiguous.
       * A map is also built, mapping each allele to the subpopulation that has it.
       * @return Vector of alleles, ordered by sample.
       */
      const std::vector<std::string> &genotypes();

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
      std::vector<T> fmt_tag(std::string tag) {
          return FormatField<T>(_header, _curr_rec, tag).values;
      }

      /**
       * @brief
       * Return the population set that has the allele.
       * @details
       * The returned vector has the same size as number of genotypes (samples * 2).
       * When true, that individual/phase possed that allele.
       * @param allele allele to get the population of
       * @return Population of indviduals that have the allele
       */
      const Population &allele_pop(const std::string allele) const override {
          return _genotype_indivs.at(allele);
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
      void create_ingroup(const std::vector<std::string> &samples) override {
          _ingroup = samples;
          _apply_ingroup_filter();
      }

      /**
       * @return string passed to htslib for ingroup filtering.
       */
      const std::vector<std::string> &ingroup() const {
          return _ingroup;
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

      htsFile *_bcf = nullptr;
      bcf_hdr_t *_header = nullptr;
      bcf1_t *_curr_rec = bcf_init();

      std::vector<std::string> _genotypes; // restricted to _ingroup
      std::unordered_map<std::string, Population> _genotype_indivs;
      std::vector<std::string> _alleles;
      std::vector<std::string> _samples;
      std::vector<std::string> _ingroup; // subset of _samples
      char *_ingroup_cstr = nullptr;

  };

  inline std::ostream &operator<<(std::ostream &os, const VCF &vcf) {
      os << "POS: " << vcf.pos() << " REF: " << vcf.ref() << " ALTS: ";
      for (const auto &al : vcf.alleles()) os << al << ' ';
      return os;
  }

}

#endif //VARGAS_VARFILE_H
