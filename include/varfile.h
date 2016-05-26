//
// Created by gaddra on 5/26/16.
//

#ifndef VARGAS_VARFILE_H
#define VARGAS_VARFILE_H

#include <string>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <htslib/vcf.h>
#include <stdlib.h>
#include <time.h>

#include "utils.h"
#include "../htslib/htslib/vcfutils.h"
#include "../doctest/doctest/doctest.h"

namespace vargas {
class VarFile {
 public:

  /**
  * Get the specified format field from the record.
  * @param T Valid types are int, char, float.
  */
  template<typename T>
  class FormatField {
   public:
    /**
     * Get the specified field.
     * @param hdr VCF Header
     * @param rec current record
     * @param tag Field to get, e.g. "GT"
     */
    FormatField(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag) : tag(tag) {
      if (!hdr || !rec || tag.length() == 0) throw std::invalid_argument("Invalid header, rec, or tag.");

      T *dst = NULL;
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
    std::string tag; // Type of FORMAT or INFO field.

   private:
    // Change the parse type based on what kind of type we have
    int _get_vals(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, int32_t **dst, int &ndst) {
      return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_INT);
    }
    int _get_vals(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, float **dst, int &ndst) {
      return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_REAL);
    }
    int _get_vals(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, char **dst, int &ndst) {
      return bcf_get_format_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_STR);
    }
  };

  /**
   * Get the specified info field from the record.
   * @param T Valid types are int, char, float.
   */
  template<typename T>
  class InfoField {
   public:
    /**
     * Get the specified field.
     * @param hdr VCF Header
     * @param rec current record
     * @param tag Field to get, e.g. "GT"
     */
    InfoField(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag) : tag(tag) {
      if (!hdr || !rec || tag.length() == 0) throw std::invalid_argument("Invalid header, rec, or tag.");

      T *dst = NULL;
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
    std::string tag; // Type of FORMAT or INFO field.

   private:
    // Change the parse type based on what kind of type we have
    int _bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, int32_t **dst, int &ndst) {
      return bcf_get_info_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_INT);
    }
    int _bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, float **dst, int &ndst) {
      return bcf_get_info_values(hdr, rec, tag.c_str(), (void **) dst, &ndst, BCF_HT_REAL);
    }
    int _bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *rec, std::string tag, char **dst, int &ndst) {
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
   * Open the specified VCF or BCF file.
   * Do not use this to load a new file after another one, rather create a new VarFile object.
   * @param file filename
   * @return -1 on file open error, -2 on header load error, 0 otherwise
   */
  int open(std::string file) {
    _file_name = file;
    return _init();
  }

  void close() {
    if (_bcf) hts_close(_bcf);
    if (_header) bcf_hdr_destroy(_header);
    bcf_destroy(_curr_rec);
    _curr_rec = bcf_init();
    _alleles.clear();
    _file_name = "";
    _chr = "";
    _min_pos = _max_pos = -1;
    _genotypes.clear();
    if (_ingroup_cstr) free(_ingroup_cstr);
  }

  /**
 * Set the minimum and maximum position
 * @param min Minimum position, 0 indexed
 * @param max Max position, inclusive, 0 indexed
 */
  void set_region(std::string chr, int min, int max) {
    _min_pos = min;
    _max_pos = max;
    _chr = chr;
  }

  /**
   * Parse a region string in the format
   * CHR:XX,XXX-YY,YYY
   * commas are stripped, range is inclusive. 0 indexed.
   * @param region region string
   */
  void set_region(std::string region) {
    std::vector<std::string> regionSplit = split(region, ':');

    // Name
    if (regionSplit.size() != 2) throw std::invalid_argument("Invalid region format, should be CHR:XX,XXX-YY,YYY");
    _chr = regionSplit[0];

    // Strip commas
    regionSplit.erase(std::remove(regionSplit.begin(), regionSplit.end(), ","), regionSplit.end());

    // Range
    regionSplit = split(regionSplit[1], '-');
    if (regionSplit.size() != 2) throw std::invalid_argument("Invalid region format, should be CHR:XX,XXX-YY,YYY");
    _min_pos = std::stoi(regionSplit[0]);
    _max_pos = std::stoi(regionSplit[1]);

  }

  /**
   * Get a list of sequences in the VCF file.
   * @return vector of sequence names
   */
  std::vector<std::string> sequences() const {
    if (!_header) return std::vector<std::string>();
    std::vector<std::string> ret;
    int num;
    const char **seqnames = bcf_hdr_seqnames(_header, &num);
    if (!seqnames) {
      throw std::invalid_argument("Error reading VCF header.");
    }
    for (int i = 0; i < num; ++i) {
      ret.push_back(std::string(seqnames[i]));
    }
    free(seqnames);
    return ret;
  }

  /**
   * @return Number of samples the VCF has. Each sample represents two genotypes.
   */
  int num_samples() const {
    if (!_header) return -1;
    return bcf_hdr_nsamples(_header);
  }

  const std::vector<std::string> &samples() const {
    return _samples;
  }

  /**
   * Load the next VCF record. Only shared info is unpacked.
   * @return false on read error or if outside restriction range.
   */
  bool next() {
    if (!_header || !_bcf) return false;
    int status = bcf_read(_bcf, _header, _curr_rec);
    unpack_all();

    // Check if its within the filter range
    if (_max_pos > 0 && _curr_rec->pos > _max_pos) return false;
    if (_curr_rec->pos < _min_pos ||
        (_chr.length() != 0 && strcmp(_chr.c_str(), bcf_hdr_id2name(_header, _curr_rec->rid)))) {
      next();
    }

    return status >= 0;
  }

  /**
   * Unpack only the shared information, and loads ref and allele info.
   */
  void unpack_shr() {
    bcf_unpack(_curr_rec, BCF_UN_SHR);
    _load_shared();
  }

  /**
   * Unpacks shared information as well as all sample information.
   * Subject to sample set restrictions.
   */
  void unpack_all() {
    bcf_unpack(_curr_rec, BCF_UN_ALL);
    _load_shared();
  }

  /**
   * Reference allele of the current record.
   * @param reference allele
   */
  std::string ref() const { return _alleles[0]; }

  /**
   * List of all the alleles in the current record. The first
   * is the reference.
   * @return vector of alleles
   */
  const std::vector<std::string> &alleles() const { return _alleles; }

  /**
   * 0 based position, i.e. the VCF pos - 1.
   * @return position.
   */
  int pos() const { return _curr_rec->pos; }

  // TODO Better handling policy of unknown elements?
  /**
   * Get a list of alleles for all samples (subject to sample set restriction).
   * Consecutive alleles represent phasing, e.g. all odd indexes are one phase,
   * all even indexes are the other. Call will unpack the full record.
   * Explicit copy number variations are replaced, other ambigous types are replaced
   * ambiguous.
   * @return Vector of alleles, ordered by sample.
   */
  const std::vector<std::string> &genotypes() {
    FormatField<int> gt(_header, _curr_rec, "GT");
    _genotypes.clear();
    for (int o : gt.values) {
      _genotypes.push_back(_alleles[bcf_gt_allele(o)]);
    }
    return _genotypes;
  }

  /**
   * Check if the file is properly loaded.
   * @return true if file is open and has a valid header.
   */
  bool good() const { return !_header && !_bcf; }

  /**
   * File name of VCF/BCF file.
   * @return file name
   */
  std::string file() const { return _file_name; }

  /**
   * Random seed used for ingroup generation.
   * @return seed
   */
  time_t seed() const { return _seed; }

  /**
   * Set the random seed.
   * @param seed random generator seed.
   */
  void seed(time_t seed) {
    _seed = seed;
    srand(seed);
  }

  /**
   * Create a random subset of the samples.
   * @param percent of samples to keep.
   */
  void create_ingroup(int percent) {
    _ingroup.clear();

    if (percent == 100) {
      _ingroup = _samples;
    }
    else if (percent != 0) {
      for (auto s : _samples) {
        if (rand() % 100 < percent) _ingroup.push_back(s);
      }
    }

    _apply_ingroup_filter();
  }

  /**
   * Include only the provided sample names in the graph.
   * @param vector of sample names
   */
  void create_ingroup(const std::vector<std::string> &samples) {
    _ingroup = samples;
    _apply_ingroup_filter();
  }

  const std::vector<std::string> &ingroup() const { return _ingroup; }

 protected:
  int _init() {
    srand(_seed);

    _bcf = bcf_open(_file_name.c_str(), "r");
    if (_bcf == NULL) return -1;
    _header = bcf_hdr_read(_bcf);
    if (!_header) return -2;

    // Load samples
    for (int i = 0; i < num_samples(); ++i) {
      _samples.push_back(_header->samples[i]);
    }

    return 0;
  }

  void _load_shared() {
    _alleles.clear();
    for (int i = 0; i < _curr_rec->n_allele; ++i) {
      std::string allele(_curr_rec->d.allele[i]);
      // Some replacement tag
      if (allele.at(0) == '<') {
        std::string ref = _curr_rec->d.allele[0];
        // Copy number
        if (allele.substr(1, 2) == "CN" && allele.at(3) != 'V') {
          int copy = std::stoi(allele.substr(3, allele.length() - 4));
          allele = "";
          for (int i = 0; i < copy; ++i) allele += ref;
        } else {
          // Other types are just subbed with the ref.
          allele = ref;
        }
      }
      _alleles.push_back(allele);
    }
  }

  void _apply_ingroup_filter() {
    if (!_header) {
      throw std::logic_error("Ingroup filter should only be applied after loading header!");
    }

    if (_ingroup.size() == 0) {
      if (_ingroup_cstr) free(_ingroup_cstr);
      _ingroup_cstr = NULL;
    }

    else if (_ingroup.size() == _samples.size()) {
      if (_ingroup_cstr) free(_ingroup_cstr);
      _ingroup_cstr = (char *) malloc(1);
      strcpy(_ingroup_cstr, "-");
    }

    else {
      std::stringstream ss;
      for (auto s : _ingroup) ss << s << ',';
      std::string smps = ss.str().substr(0, ss.str().length() - 1);
      if (_ingroup_cstr) free(_ingroup_cstr);
      _ingroup_cstr = (char *) malloc(smps.length() + 1);
      strcpy(_ingroup_cstr, smps.c_str());
    }

    bcf_hdr_set_samples(_header, _ingroup_cstr, 0);
  };


 private:
  time_t _seed = time(NULL);
  std::string _file_name, // VCF file name
      _chr; // Sequence restriction
  int _min_pos = -1, _max_pos = -1; // Region of _chr restriction

  htsFile *_bcf = NULL;
  bcf_hdr_t *_header = NULL;
  bcf1_t *_curr_rec = bcf_init();

  std::vector<std::string> _genotypes;
  std::vector<std::string> _alleles;
  std::vector<std::string> _samples;
  std::vector<std::string> _ingroup;
  char *_ingroup_cstr = NULL;

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
        << "x\t9\t.\tGG\tA,C,T\t99\t.\tAF=0.01,0.6,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t0|1\t2|3" << endl
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
        CHECK(vcf.genotypes()[0] == "GG");
        CHECK(vcf.genotypes()[1] == "A");
        CHECK(vcf.genotypes()[2] == "C");
        CHECK(vcf.genotypes()[3] == "T");
        REQUIRE(vcf.alleles().size() == 4);
        CHECK(vcf.alleles()[0] == "GG");
        CHECK(vcf.alleles()[1] == "A");
        CHECK(vcf.alleles()[2] == "C");
        CHECK(vcf.alleles()[3] == "T");
        CHECK(vcf.ref() == "GG");
        CHECK(vcf.pos() == 8);

    // Copy number alleles
    vcf.next();
        REQUIRE(vcf.genotypes().size() == 4);
        CHECK(vcf.genotypes()[0] == "CC");
        CHECK(vcf.genotypes()[1] == "CC");
        CHECK(vcf.genotypes()[2] == "");
        CHECK(vcf.genotypes()[3] == "CC");
        REQUIRE(vcf.alleles().size() == 3);
        CHECK(vcf.alleles()[0] == "C");
        CHECK(vcf.alleles()[1] == "CC");
        CHECK(vcf.alleles()[2] == "");
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
        CHECK(vcf.ref() == "GG");
    vcf.next();
        CHECK(vcf.ref() == "C");
    vcf.next();
        CHECK(vcf.ref() == "G");
        CHECK(vcf.next() == 0); // Region end
  }

      SUBCASE("Ingroup generation") {
    VarFile vcf;
    vcf.seed(12345);
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
        CHECK(vcf.genotypes()[0] == "");
        CHECK(vcf.genotypes()[1] == "CC");

    // Allele set should be complete, ingroup should reflect minimized set
        CHECK(vcf.alleles().size() == 3);
        CHECK(vcf.ingroup().size() == 1);
  }

  remove(tmpvcf.c_str());
}
}

#endif //VARGAS_VARFILE_H
