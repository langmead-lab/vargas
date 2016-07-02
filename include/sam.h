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

#include "doctest.h"

namespace Vargas {

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
           * Y is the format, Z is the value. If the format is not correct, it is ignored.
           * @param a aux field
           */
          void add(std::string a) {
              if (a.at(2) != ':') return;
              if (a.at(4) != ':') return;
              std::string tag = a.substr(0, 2);
              aux[tag] = a.substr(5);
              aux_fmt[tag] = a.at(3);
          }

          /**
           * @brief
           * Add a char type aux tag.
           * @param tag key
           * @param val value
           */
          void set(std::string tag, char val) {
              aux[tag] = std::string(1, val);
              aux_fmt[tag] = 'A';
          }

          /**
           * @brief
           * Add a int type aux tag.
           * @param tag key
           * @param val value
           */
          void set(std::string tag, int val) {
              aux[tag] = std::to_string(val);
              aux_fmt[tag] = 'i';
          }

          /**
           * @brief
           * Add a float type aux tag.
           * @param tag key
           * @param val value
           */
          void set(std::string tag, float val) {
              aux[tag] = std::to_string(val);
              aux_fmt[tag] = 'f';
          }

          /**
           * @brief
           * Add a string type aux tag.
           * @param tag key
           * @param val value
           */
          void set(std::string tag, std::string val) {
              aux[tag] = val;
              aux_fmt[tag] = 'Z';
          }

          /**
           * Get a char type tag.
           * @param tag key
           * @param val to store result in
           * @return false if key does not exist, or is the wrong type.
           */
          bool get(std::string tag, char &val) const {
              if (aux.count(tag) == 0 || aux_fmt.at(tag) != 'A') return false;
              val = aux.at(tag).at(0);
              return true;
          }

          /**
           * Get a int type tag.
           * @param tag key
           * @param val to store result in
           * @return false if key does not exist, or is the wrong type.
           */
          bool get(std::string tag, int &val) const {
              if (aux.count(tag) == 0 || aux_fmt.at(tag) != 'i') return false;
              val = std::stoi(aux.at(tag));
              return true;
          }

          /**
           * Get a float type tag.
           * @param tag key
           * @param val to store result in
           * @return false if key does not exist, or is the wrong type.
           */
          bool get(std::string tag, float &val) const {
              if (aux.count(tag) == 0 || aux_fmt.at(tag) != 'f') return false;
              val = std::stof(aux.at(tag));
              return true;
          }

          /**
           * Get a string type tag.
           * @param tag key
           * @param val to store result in
           * @return false if key does not exist, or is the wrong type.
           */
          bool get(std::string tag, std::string &val) const {
              if (aux.count(tag) == 0 || aux_fmt.at(tag) != 'Z') return false;
              val = aux.at(tag);
              return true;
          }

          /**
           * @brief
           * returns all the optional fields formatted into a line. Tab leading.
           * @return formatted fields
           */
          std::string to_string() const {
              std::ostringstream ss;
              for (auto &pair : aux) {
                  ss << '\t' << pair.first
                      << ':' << aux_fmt.at(pair.first)
                      << ':' << pair.second;
              }
              return ss.str();
          }

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
          Header() : version("1.0") { }

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
              Sequence(std::string line) { parse(line); }

              int len;
              /**< Length of the sequence */
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
              std::string to_string() const {
                  std::ostringstream ss;
                  ss << "@SQ" <<
                      "\tSN:" << name <<
                      "\tLN:" << std::to_string(len) <<
                      (genome_assembly.length() > 0 ? "\tAS:" + genome_assembly : "") <<
                      (md5.length() > 0 ? "\tM5:" + md5 : "") <<
                      (species.length() > 0 ? "\tSP:" + species : "") <<
                      (URI.length() > 0 ? "\tUR:" + URI : "") << aux.to_string();
                  return ss.str();
              }

              /**
               * @brief
               * Parse a line and populate tags.
               * @param line space delimited tag-pairs, colon delimited pairs
               */
              void parse(std::string line) {
                  genome_assembly = "";
                  md5 = "";
                  species = "";
                  URI = "";
                  aux.clear();

                  std::vector<std::string> tags = split(line, '\t');
                  for (auto &p : tags) {
                      std::string tag = p.substr(0, 2);
                      std::string val = p.substr(3);
                      if (tag == "SN") {
                          name = val;
                      }
                      else if (tag == "LN") {
                          len = std::stoi(val);
                      }
                      else if (tag == "AS") {
                          genome_assembly = val;
                      }
                      else if (tag == "M5") {
                          md5 = val;
                      }
                      else if (tag == "SP") {
                          species = val;
                      }
                      else if (tag == "UR") {
                          URI = val;
                      }
                      else {
                          aux.add(p);
                      }
                  }
              }

              /**
               * @brief
               * Parse the given line.
               * @param line to parse.
               */
              void operator<<(std::string line) { parse(line); }

              /**
               * @brief print the header line into the string.
               * @param line string to populate
               */
              void operator>>(std::string &line) const { line = to_string(); }
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
              ReadGroup(std::string line) { parse(line); }

              std::string seq_center = "", /**< Name of sequencing center producing the read */
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
              std::string to_string() const {
                  std::ostringstream ss;
                  ss << "@RG" <<
                      "\tID:" << id <<
                      (seq_center.length() > 0 ? "\tCN:" + seq_center : "") <<
                      (desc.length() > 0 ? "\tDS:" + desc : "") <<
                      (date.length() > 0 ? "\tDT:" + date : "") <<
                      (flow_order.length() > 0 ? "\tFO:" + flow_order : "") <<
                      (key_seq.length() > 0 ? "\tKS:" + key_seq : "") <<
                      (library.length() > 0 ? "\tLB:" + library : "") <<
                      (programs.length() > 0 ? "\tPG:" + programs : "") <<
                      (insert_size.length() > 0 ? "\tPI:" + insert_size : "") <<
                      (platform.length() > 0 ? "\tPL:" + platform : "") <<
                      (platform_model.length() > 0 ? "\tPM:" + platform_model : "") <<
                      (platform_unit.length() > 0 ? "\tPU:" + platform_unit : "") <<
                      (sample.length() > 0 ? "\tSM:" + sample : "") << aux.to_string();
                  return ss.str();
              }

              /**
               * @brief
               * Parse a line and populate tags.
               * @param line space delimited tag-pairs, colon delimited pairs
               */
              void parse(std::string line) {
                  seq_center = "";
                  desc = "";
                  date = "";
                  flow_order = "";
                  key_seq = "";
                  library = "";
                  programs = "";
                  insert_size = "";
                  platform = "";
                  platform_model = "";
                  platform_unit = "";
                  sample = "";
                  aux.clear();

                  std::vector<std::string> tags = split(line, '\t');
                  for (auto &p : tags) {
                      std::string tag = p.substr(0, 2);
                      std::string val = p.substr(3);
                      if (tag == "ID") {
                          id = val;
                      }
                      else if (tag == "CN") {
                          seq_center = val;
                      }
                      else if (tag == "DS") {
                          desc = val;
                      }
                      else if (tag == "DT") {
                          date = val;
                      }
                      else if (tag == "FO") {
                          flow_order = val;
                      }
                      else if (tag == "KS") {
                          key_seq = val;
                      }
                      else if (tag == "LB") {
                          library = val;
                      }
                      else if (tag == "PG") {
                          programs = val;
                      }
                      else if (tag == "PI") {
                          insert_size = val;
                      }
                      else if (tag == "PL") {
                          platform = val;
                      }
                      else if (tag == "PM") {
                          platform_model = val;
                      }
                      else if (tag == "PU") {
                          platform_unit = val;
                      }
                      else if (tag == "SM") {
                          sample = val;
                      }
                      else {
                          aux.add(p);
                      }
                  }
              }

              /**
               * @brief
               * Parse the given line.
               * @param line to parse.
               */
              void operator<<(std::string line) { parse(line); }

              /**
               * @brief print the header line into the string.
               * @param line string to populate
               */
              void operator>>(std::string &line) const { line = to_string(); }
          };

          /**
           * @brief
           * Describes a program used on the data.
           */
          struct Program {
              Program() { }
              Program(std::string line) { parse(line); }

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
              std::string to_string() const {
                  std::ostringstream ss;
                  ss << "@PG" <<
                      "\tID:" << id <<
                      (name.length() > 0 ? "\tPN:" + name : "") <<
                      (command_line.length() > 0 ? "\tCL:" + command_line : "") <<
                      (prev_pg.length() > 0 ? "\tPP:" + prev_pg : "") <<
                      (desc.length() > 0 ? "\tDS:" + desc : "") <<
                      (version.length() > 0 ? "\tVN:" + version : "") << aux.to_string();
                  return ss.str();
              }

              /**
               * @brief
               * Parse a line and populate tags.
               * @param line space delimited tag-pairs, colon delimited pairs
               */
              void parse(std::string line) {
                  name = "";
                  command_line = "";
                  prev_pg = "";
                  desc = "";
                  version = "";
                  aux.clear();

                  std::vector<std::string> tags = split(line, '\t');
                  for (auto &p : tags) {
                      std::string tag = p.substr(0, 2);
                      std::string val = p.substr(3);
                      if (tag == "ID") {
                          id = val;
                      }
                      else if (tag == "PN") {
                          name = std::stoi(val);
                      }
                      else if (tag == "CL") {
                          command_line = val;
                      }
                      else if (tag == "PP") {
                          prev_pg = val;
                      }
                      else if (tag == "DS") {
                          desc = val;
                      }
                      else if (tag == "VN") {
                          version = val;
                      }
                      else {
                          aux.add(p);
                      }
                  }
              }

              /**
               * @brief
               * Parse the given line.
               * @param line to parse.
               */
              void operator<<(std::string line) { parse(line); }

              /**
               * @brief print the header line into the string.
               * @param line string to populate
               */
              void operator>>(std::string &line) const { line = to_string(); }
          };

          /************************END TAGS************************/


          /**
           * @brief
           * Add a new sequence line.
           * @param s Sequence
           */
          void add(const Sequence &s) { sequences.push_back(s); }

          /**
           * @brief
           * Add a new ReadGroup line.
           * @param rg ReadGroup
           */
          void add(const ReadGroup &rg) { read_groups.push_back(rg); }

          /**
           * @brief
           * Add a new Program line.
           * @param pg Program
           */
          void add(const Program &pg) { programs.push_back(pg); }

          /**
            * @brief
            * Output a multi-line header containing all tags. Terminated with newline.
            * @return formatted line
            */
          std::string to_string() const {
              std::ostringstream ret;
              ret << "@HD" <<
                  "\tVN:" << version <<
                  ((sorting_order.length() > 0) ? std::string("\tSO:") + sorting_order : "") <<
                  ((grouping.length() > 0) ? std::string("\tGO:") + grouping : "") <<
                  "\n";
              for (auto &seq : sequences) {
                  ret << seq.to_string() << "\n";
              }
              for (auto &rg : read_groups) {
                  ret << rg.to_string() << "\n";
              }
              for (auto &pg : programs) {
                  ret << pg.to_string() << "\n";
              }
              return ret.str();
          }

          /**
           * @brief
           * Parse a header and populate tags.
           * @param hdr space delimited tag-pairs, colon delimited pairs
           */
          void parse(std::string hdr) {
              sorting_order = "";
              grouping = "";
              sequences.clear();
              read_groups.clear();
              programs.clear();

              std::vector<std::string> lines = split(hdr, '\n');

              // @HD line
              std::vector<std::string> tags = split(lines[0], '\t');
              if (tags[0] != "@HD") throw std::invalid_argument("First line must start with \"@HD\"");
              for (auto &p : tags) {
                  std::string tag = p.substr(0, 2);
                  std::string val = p.substr(3);
                  if (tag == "VN") {
                      version = val;
                  }
                  else if (tag == "SO") {
                      sorting_order = val;
                  }
                  else if (tag == "GO") {
                      grouping = val;
                  }
              }

              // Other lines
              for (auto &l : lines) {
                  std::string tag = l.substr(0, 3);
                  if (tag == "@SQ") add(Sequence(l));
                  else if (tag == "@RG") add(ReadGroup(l));
                  else if (tag == "@PG") add(Program(l));
              }
          }

          /**
           * @brief
           * Parse the header.
           * @param hdr header to parse.
           */
          void operator<<(std::string hdr) { parse(hdr); }

          /**
           * @brief print the header into the string.
           * @param hdr string to populate
           */
          void operator>>(std::string &hdr) const { hdr = to_string(); }

          /**
           * @brief
           * Load a string.
           * @param line header string
           */
          void operator=(std::string line) { parse(line); }

          std::string version;
          /**< SAM spec version */
          std::string sorting_order = "", /**< Type of alignment sorting, default unknown */
              grouping = "";
          /**< Grouping of alignments. Default None */
          std::vector<Sequence> sequences;
          /**< All sequence lines */
          std::vector<ReadGroup> read_groups;
          /**< All Read Group lines */
          std::vector<Program> programs; /**< ALl program lines */

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
          Record(std::string line) { parse(line); }

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
              Flag(unsigned int f) : multiple(f & 0x001),
                                     aligned(f & 0x002),
                                     unmapped(f & 0x004),
                                     next_unmapped(f & 0x008),
                                     rev_complement(f & 0x010),
                                     next_rev_complement(f & 0x020),
                                     first(f & 0x040),
                                     last(f & 0x080),
                                     secondary(f & 0x100),
                                     pass_fail(f & 0x200),
                                     duplicate(f & 0x400),
                                     supplementary(f & 0x800) { }

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
              int encode() const {
                  return (multiple ? 0x001 : 0) +
                      (aligned ? 0x002 : 0) +
                      (unmapped ? 0x004 : 0) +
                      (next_unmapped ? 0x008 : 0) +
                      (rev_complement ? 0x010 : 0) +
                      (next_rev_complement ? 0x020 : 0) +
                      (first ? 0x040 : 0) +
                      (last ? 0x080 : 0) +
                      (secondary ? 0x100 : 0) +
                      (pass_fail ? 0x200 : 0) +
                      (duplicate ? 0x400 : 0) +
                      (supplementary ? 0x800 : 0);
              }

              /**
               * @brief
               * Decode bit flags
               * @param f bit flags
               */
              void decode(unsigned int f) {
                  multiple = f & 0x001;
                  aligned = f & 0x002;
                  unmapped = f & 0x004;
                  next_unmapped = f & 0x008;
                  rev_complement = f & 0x010;
                  next_rev_complement = f & 0x020;
                  first = f & 0x040;
                  last = f & 0x080;
                  secondary = f & 0x100;
                  pass_fail = f & 0x200;
                  duplicate = f & 0x400;
                  supplementary = f & 0x800;
              }

              /**
               * @brief
               * Encode flags into an int
               * @param i int to store in.
               */
              void operator>>(unsigned int &i) const { i = encode(); }

              /**
               * Decode bit flag into the struct.
               * @param i bits to decode.
               */
              void operator<<(unsigned int i) { decode(i); }

              /**
               * @brief
               * Decode the integer into the struct.
               * @param i int to decode.
               */
              void operator=(unsigned int i) { decode(i); }
          };

          // Mandatory fields
          std::string query_name = "*", /**< Query template name */
              ref_name = "*", /**< Reference sequence name */
              cigar = "*", /**< Alignment CIGAR */
              ref_next = "*", /**< Refrence name of next mate/read */
              seq = "*", /**< segment sequence */
              qual = "*";
          /**< Phred Qual+33 */

          int pos = 0, /**< 1 baset leftmost mapping position */
              mapq = 255, /**< mapping quality */
              pos_next = 0, /**< position of next mate/read */
              tlen = 0;
          /**< template length */

          Flag flag; /**< Bitwise flag */

          Optional aux;



          /**
           * @brief
           * Output the record in single line format
           * @return single line string
           */
          std::string to_string() const {
              std::ostringstream ss;
              ss << query_name << '\t'
                  << flag.encode() << '\t'
                  << ref_name << '\t'
                  << pos << '\t'
                  << mapq << '\t'
                  << cigar << '\t'
                  << ref_next << '\t'
                  << pos_next << '\t'
                  << tlen << '\t'
                  << seq << '\t'
                  << qual << aux.to_string();
              return ss.str();
          }

          /**
           * @brief
           * Parse the line and populate fields.
           */
          void parse(std::string line) {
              std::vector<std::string> cols = split(line, '\t');
              if (cols.size() < 11) throw std::invalid_argument("Record should have at least 11 columns");
              query_name = cols[0];
              flag = std::stoi(cols[1]);
              ref_name = cols[2];
              pos = std::stoi(cols[3]);
              mapq = std::stoi(cols[4]);
              cigar = cols[5];
              ref_next = cols[6];
              pos_next = std::stoi(cols[7]);
              tlen = std::stoi(cols[8]);
              seq = cols[9];
              qual = cols[10];
              aux.clear();

              // Aux fields
              for (size_t i = 11; i < cols.size(); ++i) {
                  aux.add(cols[i]);
              }

          }

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
          void operator<<(std::string line) { parse(line); }

          /**
           * @brief
           * Populate fields from a line
           * @param line SAM record line
           */
          void operator=(std::string line) { parse(line); }

      };

    protected:
      bool _use_stdio;
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

      /**
       * Read
       */
      isam() { }
      isam(std::string file_name) { open(file_name); }
      ~isam() { close(); }

      /**
       * @brief
       * Close any open file and open the given file.
       * @param file_name SAM file to open
       */
      void open(std::string file_name) {
          close();
          if (file_name.length() == 0) _use_stdio = true;
          else {
              _use_stdio = false;
              in.open(file_name);
              if (!in.good()) throw std::invalid_argument("Error opening file \"" + file_name + "\"");
          }
          std::ostringstream hdr;
          while (std::getline((_use_stdio ? std::cin : in), _curr_line) && _curr_line.at(0) == '@') {
              hdr << _curr_line << '\n';
          }
          _hdr << hdr.str();
          _pprec << _curr_line;
      }

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
      bool next() {
          if (!std::getline((_use_stdio ? std::cin : in), _curr_line)) return false;
          _pprec.parse(_curr_line);
          return true;
      }

      /**
       * @brief
       * First record is available after opening the file.
       * @return current Record
       */
      const SAM::Record &record() const { return _pprec; }

      /**
       * Get the SAM Header.
       * @return SAM::Header
       */
      const SAM::Header &header() const { return _hdr; }

    private:
      std::string _curr_line;
      std::ifstream in;

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
      osam(std::string file_name, const SAM::Header &hdr) {
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
       */
      void open(std::string file_name) {
          close();
          if (file_name.length() == 0) _use_stdio = true;
          else {
              _use_stdio = false;
              out.open(file_name);
              if (!out.good()) throw std::invalid_argument("Error opening output file \"" + file_name + "\"");
          }
          (_use_stdio ? std::cout : out) << _hdr.to_string() << std::flush;
      }

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
      bool good() const { return out.good() || _use_stdio; }

      /**
       * @brief
       * Writes a record.
       * @param r
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

TEST_CASE ("SAM File") {
    {
        std::ofstream ss("tmp_s.sam");
        ss << "@HD\tVN:1.0\tSO:coordinate\n"
            << "@SQ\tSN:1\tLN:249250621\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta"
            << "\tM5:1b22b98cdeb4a9304cb5d48026a85128\n"
            << "@SQ\tSN:2S\tLN:243199373\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta"
            << "\tM5:a0d9851da00400dec1098a9255ac712e\n"
            << "@SQ\tSN:3Q\tLN:198022430\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta"
            << "\tM5:fdfd811849cc2fadebc929bb925902e5\n"
            << "@RG\tID:UM0098:1\tPL:ILLUMINA\tPU:HWUSI-EAS1707-615LHAAXX-L001\tLB:80\tDT:2010-05-05T20:00:00-0400"
            << "\tSM:SD37743\tCN:UMCORE\n"
            << "@RG\tID:UM0098:2\tPL:ILLUMINA\tPU:HWUSI-EAS1707-615LHAAXX-L002\tLB:80\tDT:2010-05-05T20:00:00-0400"
            << "\tSM:SD37743\tCN:UMCORE\n"
            << "@PG\tID:bwa\tVN:0.5.4\n"
            << "1:497:R:-272+13M17D24M\t113\t1\t497\t37\t37M\t15\t100338662\t0"
            << "\tCGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAG\t0;==-==9;>>>>>=>>>>>>>>>>>=>>>>>>>>>>"
            << "\tXT:A:U\tNM:i:0\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:37\n"
            << "19:20389:F:275+18M2D19M\t99\t1\t17644\t0\t37M\t=\t17919\t314"
            << "\tTATGACTGCTAATAATACCTACACATGTTAGAACCAT\t>>>>>>>>>>>>>>>>>>>><<>>><<>>4::>>:<9"
            << "\tRG:Z:UM0098:1\tXT:A:R\tNM:i:0\tSM:i:0\tAM:i:0\tX0:i:4\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:37\n"
            << "19:20389:F:275+18M2D19M\t147\t1\t17919\t0\t18M2D19M\t=\t17644\t-314"
            << "\tGTAGTACCAACTGTAAGTCCTTATCTTCATACTTTGT\t;44999;499<8<8<<<8<<><<<<><7<;<<<>><<"
            << "\tXT:A:R\tNM:i:2\tSM:i:0\tAM:i:0\tX0:i:4\tX1:i:0\tXM:i:0\tXO:i:1\tXG:i:2\tMD:Z:18^CA19\n"
            << "9:21597+10M2I25M:R:-209\t83\t1\t21678\t0\t8M2I27M\t=\t21469\t-244"
            << "\tCACCACATCACATATACCAAGCCTGGCTGTGTCTTCT\t<;9<<5><<<<><<<>><<><>><9>><>>>9>>><>"
            << "\tXT:A:R\tNM:i:2\tSM:i:0\tAM:i:0\tX0:i:5\tX1:i:0\tXM:i:0\tXO:i:1\tXG:i:2\tMD:Z:35\n";
    }

    {
        Vargas::isam sf("tmp_s.sam");
        Vargas::osam os("osam.sam", sf.header());
        do {
            os.add_record(sf.record());
        } while (sf.next());
        int x;
    }

    std::ifstream a("tmp_s.sam");
    std::ifstream b("osam.sam");
    std::string A, B;

    // Orderings of tags might be different.
    while (std::getline(a, A) && std::getline(b, B)) CHECK (A.length() == B.length());

    //   remove("tmp_s.sam");
    //   remove("osam.sam");
}
#endif //VARGAS_SAM_H
