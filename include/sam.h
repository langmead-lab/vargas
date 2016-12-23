/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date May 26, 2016
 *
 * @brief
 * Provides an interface to SAM files..
 *
 * @file
 */

#ifndef VARGAS_SAM_H
#define VARGAS_SAM_H

#include "utils.h"
#include <unordered_map>
#include <vector>
#include <sstream>
#include <iostream>

namespace vargas {

  /**
   * @brief
   * Provides a structured representation of SAM data.
   * @details
   * SAM files can be either read with SAM::isam, or written with SAM::osam. Both use
   * the same core structures to store data. \n
   * Each SAM file consists of:
   * - 1 Header with:
   *    - Any number of ReadGroups
   *    - Any number of Sequences
   *    - Any number of Programs
   * - Any number of Records
   *
   * Tags defined in the spec can be directly accessed. Other tags are stored in an
   * aux Optional struct. \n
   * Usage: \n
   * @code{.cpp}
   *    // Keep records only if they belong to a certain read group.
   *    #include "sam.h"
   *
   *    std::string keep_rg = "ID_to_keep";
   *
   *    SAM::isam input("alignments.sam");
   *    SAM::Header header = input.header();
   *    auto &rg = header.read_groups;
   *
   *    rg.erase(std::remove_if(rg.begin(), rg.end(),
   *                            [](SAM::Header::ReadGroup _r){_r.id != keep_rg;}),
   *             rg.end());
   *
   *    SAM::osam output("alignments_mod.sam", header);
   *
   *    std::string rec_rg;
   *    do {
   *        auto &rec = input.record();
   *        if (rec.aux.get("RG", rec_rg) && rec_rg == keep_rg)
   *            output << red;
   *    } while(input.next());
   * @endcode
   *
   */
  class SAM {
    public:

      /**
       * @brief
       * Represents optional data fields for any Header or Record type row.
       */
      struct Optional {
          std::unordered_map<std::string, std::string> aux;
          /**< Map the two char tag to value */
          std::unordered_map<std::string, char> aux_fmt; /**< Map the two char tag to format */

          /**
           * @brief
           * Clear all tags.
           */
          void clear() {
              aux.clear();
              aux_fmt.clear();
          }

          /**
           * @brief
           * Add a aux field of form X:Y:Z, where X is the tag
           * Y is the format, Z is the value. If Y is not provided, assume string.
           * @param a aux field
           */
          void add(std::string a);

          template<typename T>
          void set(const std::string &tag, const T &val) {
              char fmt_tag = 'Z'; // Default string type
              if (std::is_floating_point<T>::value) fmt_tag = 'f';
              else if (std::is_integral<T>::value) fmt_tag = 'i';
              aux[tag] = rg::to_string(val);
              aux_fmt[tag] = fmt_tag;
          }

          template<typename T>
          bool get(const std::string &tag, T &val) const {
              if (aux.count(tag) == 0) return false;
              rg::from_string(aux.at(tag), val);
              return true;
          }


          /**
           * @brief
           * returns all the optional fields formatted into a line. Tab leading.
           * @return formatted fields
           */
          std::string to_string() const;

      };

      /**
       * @brief
       * Represents all stored information in a SAM/BAM header.
      */
      struct Header {

          /**
           * @brief
           * All tags in the header lines came from SAMv1 spec.
           * @details
           * https://samtools.github.io/hts-specs/SAMv1.pdf
           */
          Header() : version("1.0") {};

          /************************* TAGS *************************/

          /**
           * @brief
           * Defines a reference sequence
           */
          struct Sequence {

              Sequence() { }

              /**
               * @brief
               * Construct a sequence from a sequence header line.
               * @param line seq header line to parse
               */
              Sequence(std::string line) {
                  parse(line);
              }

              int len; /**< Length of the sequence */
              std::string genome_assembly = "", /**< Genome assembly identifier */
              name, /**< name of of the sequence */
              md5 = "", /**< MD5 checksum of sequence */
              species = "", /**< Species of sequence */
              URI = ""; /**< URL of sequence, can start with http or ftp */

              Optional aux;

              /**
               * @brief
               * Prints a header-formatted line containing all populated tags
               * @return formatted line
               */
              std::string to_string() const;

              /**
               * @brief
               * Parse a line and populate tags.
               * @param line space delimited tag-pairs, colon delimited pairs
               */
              void parse(std::string line);

              /**
               * @brief
               * Parse the given line.
               * @param line to parse.
               */
              void operator<<(std::string line) {
                  parse(line);
              }

              /**
               * @brief print the header line into the string.
               * @param line string to populate
               */
              void operator>>(std::string &line) const {
                  line = to_string();
              }
          };

          /**
           * @brief
           * Defines a grouping of reads.
           */
          struct ReadGroup {
              ReadGroup() { }

              /**
               * @brief
               * Construct a read group from a RG header line.
               * @param line RG line
               */
              ReadGroup(std::string line) {
                  parse(line);
              }

              std::string
              seq_center = "", /**< Name of sequencing center producing the read */
              desc = "", /**< Description */
              date = "", /**< Date run was produced */
              flow_order = "", /**< Array of bases that corresponds to flow of each read */
              key_seq = "", /**< bases that correspond to the key sequence of the read */
              library = "", /**< Library */
              programs = "", /**< Programs used to process read group */
              insert_size = "", /**< Predicted median insert size */
              platform = "", /**< Platform used to produce reads */
              platform_model = "", /**< Platform desc */
              platform_unit = "", /**< Unique unit ID */
              sample = "", /**< Sample, or pool name */
              id = ""; /**< Unique ID */

              Optional aux;

              /**
     * @brief
     * Prints a header-formatted line containing all populated tags
     * @return formatted line
     */
              std::string to_string() const;

              /**
               * @brief
               * Parse a line and populate tags.
               * @param line space delimited tag-pairs, colon delimited pairs
               */
              void parse(std::string line);

              /**
               * @brief
               * Parse the given line.
               * @param line to parse.
               */
              void operator<<(std::string line) {
                  parse(line);
              }

              /**
               * @brief print the header line into the string.
               * @param line string to populate
               */
              void operator>>(std::string &line) const {
                  line = to_string();
              }
          };

          /**
           * @brief
           * Describes a program used on the data.
           */
          struct Program {
              Program() { }

              Program(std::string line) {
                  parse(line);
              }

              std::string name = "",
              id = "",
              command_line = "",
              prev_pg = "",
              desc = "",
              version = "";

              Optional aux;

              /**
     * @brief
     * Prints a header-formatted line containing all populated tags
     * @return formatted line
     */
              std::string to_string() const;

              /**
               * @brief
               * Parse a line and populate tags.
               * @param line space delimited tag-pairs, colon delimited pairs
               */
              void parse(std::string line);

              /**
               * @brief
               * Parse the given line.
               * @param line to parse.
               */
              void operator<<(std::string line) {
                  parse(line);
              }

              /**
               * @brief print the header line into the string.
               * @param line string to populate
               */
              void operator>>(std::string &line) const {
                  line = to_string();
              }
          };

          /************************END TAGS************************/


          /**
           * @brief
           * Add a new sequence line.
           * @param s Sequence
           * @throws std::out_of_range if sequence name already exists
           */
          void add(const Sequence &s) {
              if (sequences.count(s.name) != 0) throw std::out_of_range("Sequence name already exists.");
              sequences[s.name] = s;
          }

          /**
           * @brief
           * Add a new ReadGroup line.
           * @param rg ReadGroup
           * @throws std::out_of_range if read group ID already exists
           */
          void add(const ReadGroup &rg) {
              if (read_groups.count(rg.id) != 0) throw std::out_of_range("Read group ID already exists.");
              read_groups[rg.id] = rg;
          }

          /**
           * @brief
           * Add a new Program line.
           * @throws std::out_of_range if program ID already exists
           * @param pg Program
           */
          void add(const Program &pg) {
              if (programs.count(pg.id) != 0) throw std::out_of_range("Program ID already exists.");
              programs[pg.id] = pg;
          }

          /**
            * @brief
            * Output a multi-line header containing all tags. Terminated with newline.
            * @return formatted line
            */
          std::string to_string() const;

          /**
           * @brief
           * Parse a header and populate tags.
           * @param hdr space delimited tag-pairs, colon delimited pairs
           * @throws std::invalid_argument if invalid SAM file
           */
          void parse(std::string hdr);

          /**
           * @brief
           * Parse the header.
           * @param hdr header to parse.
           */
          void operator<<(std::string hdr) {
              parse(hdr);
          }

          /**
           * @brief print the header into the string.
           * @param hdr string to populate
           */
          void operator>>(std::string &hdr) const {
              hdr = to_string();
          }

          /**
           * @brief
           * Load a string.
           * @param line header string
           */
          void operator=(std::string line) {
              parse(line);
          }

          std::string version; /**< SAM spec version */
          std::string sorting_order = "", /**< Type of alignment sorting, default unknown */
          grouping = ""; /**< Grouping of alignments. Default None */
          std::unordered_map<std::string, Sequence> sequences; /**< All sequence lines */
          std::unordered_map<std::string, ReadGroup> read_groups; /**< All Read Group lines */
          std::unordered_map<std::string, Program> programs; /**< ALl program lines */

      };

      /**
       * @brief
       * Represents an alignment record.
       */
      struct Record {

          Record() { }

          /**
           * @brief
           * Parse the given alignment line.
           * @param line
           */
          Record(std::string line) {
              parse(line);
          }

          /**
           * @param
           * Represents a bitwise flag of an alignment.
           */
          struct Flag {

              Flag() { }

              /**
               * @brief
               * Decode bitwise flag.
               * @param f bitwise flag
               */
              Flag(unsigned int f) : multiple((f & 0x001) != 0u),
                                     aligned((f & 0x002) != 0u),
                                     unmapped((f & 0x004) != 0u),
                                     next_unmapped((f & 0x008) != 0u),
                                     rev_complement((f & 0x010) != 0u),
                                     next_rev_complement((f & 0x020) != 0u),
                                     first((f & 0x040) != 0u),
                                     last((f & 0x080) != 0u),
                                     secondary((f & 0x100) != 0u),
                                     pass_fail((f & 0x200) != 0u),
                                     duplicate((f & 0x400) != 0u),
                                     supplementary((f & 0x800) != 0u) {}

              bool multiple = false, /**< template having multiple segments in sequencing */
              aligned = false, /**< each segment properly aligned according to the aligner */
              unmapped = false, /**< segment unmapped */
              next_unmapped = false, /**< next segment in the template unmapped */
              rev_complement = false, /**< SEQ being reverse complemented */
              next_rev_complement =
              false, /**< SEQ of the next segment in the template being reverse complemented */
              first = false, /**< the first segment in the template */
              last = false, /**< he last segment in the template */
              secondary = false, /**< secondary alignment */
              pass_fail = false, /**< not passing lters, such as platform/vendor quality controls */
              duplicate = false, /**< PCR or optical duplicate */
              supplementary = false; /**< supplementary alignment */

              /**
               * @return integer representation of flags.
               */
              int encode() const;

              /**
               * @brief
               * Decode bit flags
               * @param f bit flags
               */
              void decode(unsigned int f);

              /**
               * @brief
               * Encode flags into an int
               * @param i int to store in.
               */
              void operator>>(unsigned int &i) const {
                  i = encode();
              }

              /**
               * Decode bit flag into the struct.
               * @param i bits to decode.
               */
              void operator<<(unsigned int i) {
                  decode(i);
              }

              /**
               * @brief
               * Decode the integer into the struct.
               * @param i int to decode.
               */
              void operator=(unsigned int i) {
                  decode(i);
              }
          };

          // Mandatory fields
          std::string query_name = "*", /**< Query template name */
          ref_name = "*", /**< Reference sequence name */
          cigar = "*", /**< Alignment CIGAR */
          ref_next = "*", /**< Refrence name of next mate/read */
          seq = "*", /**< segment sequence */
          qual = "*"; /**< Phred Qual+33 */

          int pos = 0, /**< 1 baset leftmost mapping position */
          mapq = 255, /**< mapping quality */
          pos_next = 0, /**< position of next mate/read */
          tlen = 0; /**< template length */

          Flag flag; /**< Bitwise flag */

          Optional aux;


          /**
           * @brief
           * Output the record in single line format
           * @return single line string
           */
          std::string to_string() const;

          /**
           * @brief
           * Parse the line and populate fields.
           * @throws std::invalid_argument if record has incorrect number of fields
           */
          void parse(std::string line);

          /**
           * @brief
           * Output the record in single line format
           * @param s single line string
           */
          void operator>>(std::string &s) const {
              s = to_string();
          }

          /**
           * @brief
           * Populate fields from a line
           * @param line SAM record line
           */
          void operator<<(std::string line) {
              parse(line);
          }

          /**
           * @brief
           * Populate fields from a line
           * @param line SAM record line
           */
          void operator=(std::string line) {
              parse(line);
          }

          /**
           * @brief
           * Get a tag from the associated header read group.
           * @param hdr SAM header
           * @param tag
           * @param val result
           * @return true if tag exists
           */
          template<typename T>
          bool read_group(const Header &hdr,
                          const std::string &tag,
                          T &val) const {
              std::string rg;
              if (!aux.get("RG", rg))return false;
              if (hdr.read_groups.count(rg) == 0) return false;
              return hdr.read_groups.at(rg).aux.get(tag, val);
          }

          bool get_required(const std::string &tag, std::string &val) const;

          /**
           * @brief
           * Get the value of a tag. If prefixed with RG:, get the tag from the
           * associated read group. Otherwise check required fields before optional fields.
           * @param hdr SAM Header
           * @param tag tag
           * @param val output
           * @return true if tag was found and val updated
           */
          template<typename T>
          bool get(const Header &hdr, const std::string &tag, T &val) const {
              if (tag.length() > 3 && tag.substr(0, 3) == "RG:") {
                  return read_group(hdr, tag.substr(3), val);
              }
              std::string ret;
              if (get_required(tag, ret) || aux.get(tag, ret)) {
                  rg::from_string(ret, val);
                  return true;
              } else return false;
          }

          static const std::string REQUIRED_POS;
          static const std::string REQUIRED_QNAME;
          static const std::string REQUIRED_RNEXT;
          static const std::string REQUIRED_RNAME;
          static const std::string REQUIRED_SEQ;
          static const std::string REQUIRED_CIGAR;
          static const std::string REQUIRED_FLAG;
          static const std::string REQUIRED_PNEXT;
          static const std::string REQUIRED_MAPQ;
          static const std::string REQUIRED_TLEN;
          static const std::string REQUIRED_QUAL;

      };

    protected:
      bool _use_stdio = false;
      SAM::Header _hdr;
  };

  /**
 * @brief
 * Provides an interface to read a SAM file.
 * @details
 * Minimal error checking is done, no error will be raised if there is a mismatch (e.g. a ref name
 * that is not defined in the header). \n
 * See Vargas::SAM for more usage information.
 */
  class isam: public SAM {
    public:

      isam() {}

      isam(std::string file_name) {
          open(file_name);
      }

      isam(isam &&o) {
          _buff = std::move(o._buff);
          _pprec = std::move(o._pprec);
          _hdr = std::move(o._hdr);
          _use_stdio = o._use_stdio;
      }

      ~isam() {
          close();
      }

      /**
       * @brief
       * Close any open file and open the given stream.
       * @param is input stream
       */
      void open(std::istream &is);

      /**
       * @brief
       * Close any open file and open the given file.
       * @param file_name SAM file to open
       * @throws std::invalid_argument if file cannot be opened
       */
      void open(std::string file_name);

      /**
       * @brief
       * Clear data and close any open handles.
       */
      void close() {
          in.close();
          _hdr = SAM::Header();
          _pprec = SAM::Record();
      }

      /**
       * @return true if file is open.
       */
      bool good() const {
          return in.good() || _use_stdio;
      }

      /**
       * @brief
       * Load the next record.
       * @details
       * Since the first record is loaded at init, to loop through all the records
       * a do{}while() loop should be used.
       * @return true is a new record was obtained.
       */
      bool next();

      /**
       * @brief
       * Load the rest of the records in the file and return a subset of them.
       * @param n Number of records to keep.
       * @return isam with subset of reads.
       */
      isam subset(size_t n);

      /**
       * @brief
       * First record is available after opening the file.
       * @return current Record
       */
      const SAM::Record &record() const {
          return _pprec;
      }

      /**
       * Get the SAM Header.
       * @return SAM::Header
       */
      const SAM::Header &header() const {
          return _hdr;
      }

    private:
      std::string _curr_line;
      std::ifstream in;
      std::vector<Record> _buff;

      SAM::Record _pprec;
  };

  /**
 * @brief
 * Provides an interface to write a SAM file.
 * @details
 * Outputs are not checked to see if they conform with the standard format. \n
 * See Vargas::SAM for more information about usage.
 */
  class osam: public SAM {
    public:

      /**
       * @brief
       * Use stdout
       */
      osam(const SAM::Header &hdr) {
          _hdr = hdr;
          open("");
      }

      /**
       * @brief
       * Create a SAM file with the given header and file name.
       * @details
       * To allow streaming operation, the SAM::Header must be specified
       * at file open. This also allows for added alignments to be validated
       * against header information.
       * @param file_name file to write
       * @param hdr SAM::Header of the file
       */
      osam(std::string file_name,
           const SAM::Header &hdr) {
          _hdr = hdr;
          open(file_name);
      }

      ~osam() {
          close();
      }

      /**
       * @brief
       * Open a new file.
       * @details
       * Any added alignments are flushed to the previous file (if any). The header
       * is written to the new file.
       * @param file_name file to open
       * @throws std::invalid_argument if file cannot be opened
       */
      void open(std::string file_name);

      /**
       * @brief
       * Flush any data, and close the output file.
       */
      void close() {
          if (out.is_open()) {
              out.close();
          }
      }

      /**
       * @return true of output open.
       */
      bool good() const {
          return out.good() || _use_stdio;
      }

      /**
       * @brief
       * Writes a record.
       * @param r record to add
       * @throws std::invalid_argument if no output file open
       */
      void add_record(const SAM::Record &r) {
          if (!good()) throw std::invalid_argument("No valid file open.");
          (_use_stdio ? std::cout : out) << r.to_string() << '\n';
      }

      /**
       * @param rec write record to output
       * @return SAM::osam
       */
      osam &operator<<(const SAM::Record &rec) {
          add_record(rec);
          return *this;
      }

    private:
      std::ofstream out;
  };

}

#endif //VARGAS_SAM_H
