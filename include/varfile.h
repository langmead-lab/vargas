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
#include <time.h>
#include <unordered_map>
#include "dyn_bitset.h"
#include "utils.h"
#include "htslib/vcfutils.h"
#include "htslib/hts.h"
#include "doctest.h"

namespace vargas {

/**
 * @brief
 * Provides an interface to a VCF/BCF file. Core processing
 * provided by htslib.
 */
class VarFile {
 public:
    typedef dyn_bitset<64> Population;
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
    InfoField(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag) : tag(tag) {
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

  VarFile() { }
  VarFile(std::string file) : _file_name(file) { _init(); }
  VarFile(std::string file, std::string chr, int min, int max) :
      _file_name(file), _chr(chr), _min_pos(min), _max_pos(max) { _init(); }
  ~VarFile() {
    close();
  }

    /**
     * @brief
     * Open the specified VCF or BCF file and load the header.
     * @param file filename
     * @return -1 on file open error, -2 on header load error, 0 otherwise
     */
  int open(std::string file) {
    _file_name = file;
    return _init();
  }

  void close() {
    if (_bcf) bcf_close(_bcf);
    if (_header) bcf_hdr_destroy(_header);
    if (_curr_rec) bcf_destroy(_curr_rec);
    if (_ingroup_cstr) free(_ingroup_cstr);
    _bcf = nullptr;
    _header = nullptr;
    _curr_rec = nullptr;
    _ingroup_cstr = nullptr;
  }

    /**
     * @brief
      * Set the minimum and maximum position, inclusive.
      * If max is <= 0, go until end.
      * @param chr contig ID
      * @param min Minimum position, 0 indexed
      * @param max Max position, inclusive, 0 indexed
   */
  void set_region(std::string chr, int min, int max) {
    _min_pos = min;
    _max_pos = max;
    _chr = chr;
  }

    /**
     * @brief
     * Parse a region string in the format: \n
     * CHR:XX,XXX-YY,YYY \n
     * commas are stripped, range is inclusive. 0 indexed.
     * If max is <= 0, go until end.
     * @param region region string
     */
  void set_region(std::string region);

    /**
     * @brief
     * Get the minimum position.
     * @return 0 indexed minimum pos
     */
  int region_lower() const { return _min_pos; }

    /**
     * @brief
     * Get the maximum position.
     * @return 0 indexed maximum position
     */
  int region_upper() const { return _max_pos; }

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
     * Current contig filter.
     * @return ID of current CHROM filter
     */
  std::string region_chr() const { return _chr; }

    /**
     * @brief
     * Get a list of sequences in the VCF file.
     * @return vector of sequence names
     */
  std::vector<std::string> sequences() const;

  /**
   * @return Number of samples the VCF has. Each sample represents two genotypes.
   */
  int num_samples() const {
    if (!_header) return -1;
    return bcf_hdr_nsamples(_header);
  }

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
  std::string ref() const { return _alleles[0]; }

    /**
     * @brief
     * List of all the alleles in the current record.
     * @details
     * The first is the reference. Allele copy number variant tags are converted,
     * whereas other tags are substituted for the reference.
     * @return vector of alleles
     */
  const std::vector<std::string> &alleles() const { return _alleles; }

    /**
     * @brief
     * 0 based position, i.e. the VCF pos - 1.
     * @return position.
     */
  int pos() const { return _curr_rec->pos; }

  // TODO Better handling policy of unknown elements?
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
  const std::vector<float> &frequencies();

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
    const Population &allele_pop(std::string allele) const { return _genotype_indivs.at(allele); }

    /**
     * @brief
     * Check if the file is properly loaded.
     * @return true if file is open and has a valid header.
     */
  bool good() const { return _header && _bcf; }

    /**
     * @brief
     * File name of VCF/BCF file.
     * @return file name
     */
  std::string file() const { return _file_name; }


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
  const std::vector<std::string> &ingroup() const { return _ingroup; }

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
  std::string _file_name, // VCF/BCF file name
      _chr; // Sequence restriction (CHROM)
  int _min_pos = -1, _max_pos = -1; // Region of _chr restriction

    htsFile *_bcf = nullptr;
    bcf_hdr_t *_header = nullptr;
  bcf1_t *_curr_rec = bcf_init();

  std::vector<std::string> _genotypes; // restricted to _ingroup
  std::vector<float> _allele_freqs;
    std::unordered_map<std::string, Population> _genotype_indivs;
  std::vector<std::string> _alleles;
  std::vector<std::string> _samples;
  std::vector<std::string> _ingroup; // subset of _samples
    char *_ingroup_cstr = nullptr;

};
TEST_CASE ("VCF File handler") {
  using std::endl;
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
        << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternete Allele count\">" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Num samples at site\">" << endl
        << "##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Num alt alleles\">" << endl
        << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"Length of each alt\">" << endl
        << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"type of variant\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2" << endl
        << "x\t9\t.\tG\tA,C,T\t99\t.\tAF=0.01,0.6,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t0|1\t2|3" << endl
        << "x\t10\t.\tC\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "x\t14\t.\tG\t<DUP>,<BLAH>\t99\t.\tAF=0.01,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t1|1" << endl
        << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
  }

      SUBCASE("Unfiltered") {
    VarFile vcf(tmpvcf);
    vcf.next();
        CHECK(vcf.num_samples() == 2);
        CHECK(vcf.sequences().size() == 2);
        CHECK(vcf.sequences()[0] == "x");
        CHECK(vcf.sequences()[1] == "y");
        REQUIRE(vcf.samples().size() == 2);
        CHECK(vcf.samples()[0] == "s1");
        CHECK(vcf.samples()[1] == "s2");

    // On load, first record is already loaded
        REQUIRE(vcf.genotypes().size() == 4);
        CHECK(vcf.genotypes()[0] == "G");
        CHECK(vcf.genotypes()[1] == "A");
        CHECK(vcf.genotypes()[2] == "C");
        CHECK(vcf.genotypes()[3] == "T");
        REQUIRE(vcf.alleles().size() == 4);
        CHECK(vcf.alleles()[0] == "G");
        CHECK(vcf.alleles()[1] == "A");
        CHECK(vcf.alleles()[2] == "C");
        CHECK(vcf.alleles()[3] == "T");
        CHECK(vcf.ref() == "G");
        CHECK(vcf.pos() == 8);

    // Copy number alleles
    vcf.next();
        REQUIRE(vcf.genotypes().size() == 4);
        CHECK(vcf.genotypes()[0] == "CC");
        CHECK(vcf.genotypes()[1] == "CC");
          CHECK(vcf.genotypes()[2] == "C");
        CHECK(vcf.genotypes()[3] == "CC");
        REQUIRE(vcf.alleles().size() == 3);
        CHECK(vcf.alleles()[0] == "C");
        CHECK(vcf.alleles()[1] == "CC");
          CHECK(vcf.alleles()[2] == "C");
        CHECK(vcf.ref() == "C");
        CHECK(vcf.pos() == 9);

    // Invalid tags
    vcf.next();
        REQUIRE(vcf.alleles().size() == 3);
        CHECK(vcf.alleles()[0] == "G");
        CHECK(vcf.alleles()[1] == "G");
        CHECK(vcf.alleles()[2] == "G");
        CHECK(vcf.ref() == "G");
        CHECK(vcf.pos() == 13);

    // Next y contig should still load
    vcf.next();
        CHECK(vcf.alleles()[0] == "TATA");
  }

      SUBCASE("CHROM Filtering") {
    VarFile vcf;
    vcf.set_region("y:0-0");
    vcf.open(tmpvcf);

    vcf.next();
        CHECK(vcf.ref() == "TATA");
    vcf.next();
        CHECK(vcf.ref() == "T");
        CHECK(vcf.next() == 0); // File end
  }

      SUBCASE("Region filtering") {
    VarFile vcf;
    vcf.set_region("x:0-14");
    vcf.open(tmpvcf);

    vcf.next();
        CHECK(vcf.ref() == "G");
    vcf.next();
        CHECK(vcf.ref() == "C");
    vcf.next();
        CHECK(vcf.ref() == "G");
        CHECK(vcf.next() == 0); // Region end
  }

      SUBCASE("Ingroup generation") {
    VarFile vcf;
      srand(12345);
    vcf.open(tmpvcf);
    vcf.create_ingroup(50);

        CHECK(vcf.ingroup().size() == 1);
        CHECK(vcf.ingroup()[0] == "s2");

    vcf.next();
        REQUIRE(vcf.genotypes().size() == 2);
        CHECK(vcf.genotypes()[0] == "C");
        CHECK(vcf.genotypes()[1] == "T");

    vcf.next();
        REQUIRE(vcf.genotypes().size() == 2);
          CHECK(vcf.genotypes()[0] == "C");
        CHECK(vcf.genotypes()[1] == "CC");

    // Allele set should be complete, ingroup should reflect minimized set
        CHECK(vcf.alleles().size() == 3);
        CHECK(vcf.ingroup().size() == 1);
  }

      SUBCASE("Allele populations") {
    VarFile vcf;
    vcf.open(tmpvcf);
    vcf.next();
    vcf.genotypes();

        REQUIRE(vcf.allele_pop("G").size() == 4);
        CHECK(vcf.allele_pop("G")[0]);
        CHECK(!vcf.allele_pop("G")[1]);
        CHECK(!vcf.allele_pop("G")[2]);
        CHECK(!vcf.allele_pop("G")[3]);

        REQUIRE(vcf.allele_pop("A").size() == 4);
        CHECK(!vcf.allele_pop("A")[0]);
        CHECK(vcf.allele_pop("A")[1]);
        CHECK(!vcf.allele_pop("A")[2]);
        CHECK(!vcf.allele_pop("A")[3]);

        REQUIRE(vcf.allele_pop("C").size() == 4);
        CHECK(!vcf.allele_pop("C")[0]);
        CHECK(!vcf.allele_pop("C")[1]);
        CHECK(vcf.allele_pop("C")[2]);
        CHECK(!vcf.allele_pop("C")[3]);

        REQUIRE(vcf.allele_pop("T").size() == 4);
        CHECK(!vcf.allele_pop("T")[0]);
        CHECK(!vcf.allele_pop("T")[1]);
        CHECK(!vcf.allele_pop("T")[2]);
        CHECK(vcf.allele_pop("T")[3]);

  }

      SUBCASE("Filtered allele populations") {
    VarFile vcf;
    vcf.open(tmpvcf);
    vcf.create_ingroup({"s1"});
    vcf.next();
    vcf.genotypes();

        REQUIRE(vcf.allele_pop("G").size() == 2);
        CHECK(vcf.allele_pop("G")[0]);
        CHECK(!vcf.allele_pop("G")[1]);

        REQUIRE(vcf.allele_pop("A").size() == 2);
        CHECK(!vcf.allele_pop("A")[0]);
        CHECK(vcf.allele_pop("A")[1]);

        REQUIRE(vcf.allele_pop("C").size() == 2);
        CHECK(!vcf.allele_pop("C")[0]);
        CHECK(!vcf.allele_pop("C")[1]);

        REQUIRE(vcf.allele_pop("T").size() == 2);
        CHECK(!vcf.allele_pop("T")[0]);
        CHECK(!vcf.allele_pop("T")[1]);
  }

      SUBCASE("Allele frequencies") {
    VarFile vcf;
    vcf.open(tmpvcf);
    vcf.next();

    auto af = vcf.frequencies();
        REQUIRE(af.size() == 4);
        CHECK(af[0] > 0.289f); // af[0] should be 0.29
        CHECK(af[0] < 0.291f);
        CHECK(af[1] == 0.01f);
        CHECK(af[2] == 0.6f);
        CHECK(af[3] == 0.1f);
  }

  remove(tmpvcf.c_str());
}
}

#endif //VARGAS_VARFILE_H
