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

namespace vargas {

  /**
   * @brief
   * Provides an interface for a SAM/BAM File.
   */
  class SAM {
    public:

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

              Sequence(std::string line) { parse(line); }

              int len;
              /**< Length of the sequence */
              std::string genome_assembly = "", /**< Genome assembly identifier */
                  name, /**< name of of the sequence */
                  md5 = "", /**< MD5 checksum of sequence */
                  species = "", /**< Species of sequence */
                  URI = ""; /**< URL of sequence, can start with http or ftp */

              /**
               * @brief
               * Prints a header-formatted line containing all populated tags
               * @return formatted line
               */
              std::string to_string() const {
                  std::stringstream ss;
                  ss << "@SQ" <<
                      " SN:" << name <<
                      " LN:" << std::to_string(len) <<
                      ((genome_assembly.length() > 0) ? std::string(" AS:") + genome_assembly : "") <<
                      ((md5.length() > 0) ? std::string(" M5:" + md5) : "") <<
                      ((species.length() > 0) ? std::string(" SP:" + species) : "") <<
                      ((URI.length() > 0) ? std::string(" UR:" + URI) : "");
                  return ss.str();
              }

              /**
               * @brief
               * Parse a line and populate tags.
               * @param line space delimited tag-pairs, colon delimited pairs
               */
              void parse(std::string line) {
                  std::vector<std::string> tags = split(line, '\t');
                  for (auto &p : tags) {
                      std::vector<std::string> pair(2);
                      pair[0] = p.substr(0, 2);
                      pair[1] = p.substr(3);
                      if (pair[0] == "SN") {
                          name = pair[1];
                      }
                      else if (pair[0] == "LN") {
                          len = std::stoi(pair[1]);
                      }
                      else if (pair[0] == "AS") {
                          genome_assembly = pair[1];
                      }
                      else if (pair[0] == "M5") {
                          md5 = pair[1];
                      }
                      else if (pair[0] == "SP") {
                          species = pair[1];
                      }
                      else if (pair[0] == "UR") {
                          URI = pair[1];
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
               * @param string string to populate
               */
              void operator>>(std::string &line) const { line = to_string(); }
          };

          /**
           * @brief
           * Defines a grouping of reads.
           */
          struct ReadGroup {
              ReadGroup(std::string line) { parse(line); }

              std::string seq_center = "",
                  desc = "",
                  date = "",
                  flow_order = "",
                  key_seq = "",
                  library = "",
                  programs = "",
                  insert_size = "",
                  platform = "",
                  platform_model = "",
                  platform_unit = "",
                  sample = "",
                  id = "";

              /**
     * @brief
     * Prints a header-formatted line containing all populated tags
     * @return formatted line
     */
              std::string to_string() const {
                  std::stringstream ss;
                  ss << "@RG" <<
                      " ID:" << id <<
                      ((seq_center.length() > 0) ? std::string(" CN:") + seq_center : "") <<
                      ((desc.length() > 0) ? std::string(" DS:") + desc : "") <<
                      ((date.length() > 0) ? std::string(" DT:") + date : "") <<
                      ((flow_order.length() > 0) ? std::string(" FO:") + flow_order : "") <<
                      ((key_seq.length() > 0) ? std::string(" KS:") + key_seq : "") <<
                      ((library.length() > 0) ? std::string(" LB:") + library : "") <<
                      ((programs.length() > 0) ? std::string(" PG:") + programs : "") <<
                      ((insert_size.length() > 0) ? std::string(" PI:") + insert_size : "") <<
                      ((platform.length() > 0) ? std::string(" PL:") + platform : "") <<
                      ((platform_model.length() > 0) ? std::string(" PM:") + platform_model : "") <<
                      ((platform_unit.length() > 0) ? std::string(" PU:") + platform_unit : "") <<
                      ((sample.length() > 0) ? std::string(" SM:") + sample : "");
                  return ss.str();
              }

              /**
               * @brief
               * Parse a line and populate tags.
               * @param line space delimited tag-pairs, colon delimited pairs
               */
              void parse(std::string line) {
                  std::vector<std::string> tags = split(line, '\t');
                  for (auto &p : tags) {
                      std::vector<std::string> pair(2);
                      pair[0] = p.substr(0, 2);
                      pair[1] = p.substr(3);
                      if (pair[0] == "ID") {
                          id = pair[1];
                      }
                      else if (pair[0] == "CN") {
                          seq_center = pair[1];
                      }
                      else if (pair[0] == "DS") {
                          desc = pair[1];
                      }
                      else if (pair[0] == "DT") {
                          date = pair[1];
                      }
                      else if (pair[0] == "FO") {
                          flow_order = pair[1];
                      }
                      else if (pair[0] == "KS") {
                          key_seq = pair[1];
                      }
                      else if (pair[0] == "LB") {
                          library = pair[1];
                      }
                      else if (pair[0] == "PG") {
                          programs = pair[1];
                      }
                      else if (pair[0] == "PI") {
                          insert_size = pair[1];
                      }
                      else if (pair[0] == "PL") {
                          platform = pair[1];
                      }
                      else if (pair[0] == "PM") {
                          platform_model = pair[1];
                      }
                      else if (pair[0] == "PU") {
                          platform_unit = pair[1];
                      }
                      else if (pair[0] == "SM") {
                          sample = pair[1];
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
               * @param string string to populate
               */
              void operator>>(std::string &line) const { line = to_string(); }
          };

          /**
           * @brief
           * Describes a program used on the data.
           */
          struct Program {
              Program(std::string line) { parse(line); }

              std::string name = "",
                  id = "",
                  command_line = "",
                  prev_pg = "",
                  desc = "",
                  version = "";

              /**
     * @brief
     * Prints a header-formatted line containing all populated tags
     * @return formatted line
     */
              std::string to_string() const {
                  std::stringstream ss;
                  ss << "@PG" <<
                      " ID:" << id <<
                      ((name.length() > 0) ? std::string(" PN:") + name : "") <<
                      ((command_line.length() > 0) ? std::string(" CL:") + command_line : "") <<
                      ((prev_pg.length() > 0) ? std::string(" PP:") + prev_pg : "") <<
                      ((desc.length() > 0) ? std::string(" DS:") + desc : "") <<
                      ((version.length() > 0) ? std::string(" VN:") + version : "");
                  return ss.str();
              }

              /**
               * @brief
               * Parse a line and populate tags.
               * @param line space delimited tag-pairs, colon delimited pairs
               */
              void parse(std::string line) {
                  std::vector<std::string> tags = split(line, '\t');
                  for (auto &p : tags) {
                      std::vector<std::string> pair(2);
                      pair[0] = p.substr(0, 2);
                      pair[1] = p.substr(3);
                      if (pair[0] == "ID") {
                          id = pair[1];
                      }
                      else if (pair[0] == "PN") {
                          name = std::stoi(pair[1]);
                      }
                      else if (pair[0] == "CL") {
                          command_line = pair[1];
                      }
                      else if (pair[0] == "PP") {
                          prev_pg = pair[1];
                      }
                      else if (pair[0] == "DS") {
                          desc = pair[1];
                      }
                      else if (pair[0] == "VN") {
                          version = pair[1];
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
               * @param string string to populate
               */
              void operator>>(std::string &line) const { line = to_string(); }
          };

          /************************END TAGS************************/


          /**
           * @brief
           * Add a new sequence line.
           * @param Sequence
           */
          void add(const Sequence &s) { sequences.push_back(s); }

          /**
           * @brief
           * Add a new ReadGroup line.
           * @param ReadGroup
           */
          void add(const ReadGroup &rg) { read_groups.push_back(rg); }

          /**
           * @brief
           * Add a new Program line.
           * @param Program
           */
          void add(const Program &pg) { programs.push_back(pg); }

          /**
            * @brief
            * Output a multi-line header containing all tags.
            * @return formatted line
            */
          std::string to_string() const {
              std::stringstream ret;
              ret << "@HD" <<
                  " VN:" << version <<
                  ((sorting_order.length() > 0) ? std::string(" SO:") + sorting_order : "") <<
                  ((grouping.length() > 0) ? std::string(" GO:") + grouping : "") <<
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
           * @param line space delimited tag-pairs, colon delimited pairs
           */
          void parse(std::string hdr) {
              std::vector<std::string> lines = split(hdr, '\n');

              // @HD line
              std::vector<std::string> tags = split(lines[0], '\t');
              if (tags[0] != "@HD") throw std::invalid_argument("First line must start with \"@HD\"");
              for (auto &p : tags) {
                  std::vector<std::string> pair(2);
                  pair[0] = p.substr(0, 2);
                  pair[1] = p.substr(3);
                  if (pair[0] == "VN") {
                      version = pair[1];
                  }
                  else if (pair[0] == "SO") {
                      sorting_order = pair[1];
                  }
                  else if (pair[0] == "GO") {
                      grouping = pair[1];
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
           * @param header to parse.
           */
          void operator<<(std::string hdr) { parse(hdr); }

          /**
           * @brief print the header into the string.
           * @param string string to populate
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
              Flag(int f) : multiple(f & 0x001),
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
                  return
                      (multiple & 0x1) +
                          (aligned & 0x2) +
                          (unmapped & 0x4) +
                          (next_unmapped & 0x8) +
                          (rev_complement & 0x10) +
                          (next_rev_complement & 0x20) +
                          (first & 0x40) +
                          (last & 0x80) +
                          (secondary & 0x100) +
                          (pass_fail & 0x200) +
                          (duplicate & 0x400) +
                          (supplementary & 0x800);
              }

              /**
               * @brief
               * Decode bit flags
               * @param f bit flags
               */
              void decode(int f) {
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
              void operator>>(int &i) const { i = encode(); }

              /**
               * Decode bit flag into the struct.
               * @param i bits to decode.
               */
              void operator<<(int i) { decode(i); }

              /**
               * @brief
               * Decode the integer into the struct.
               * @param i int to decode.
               */
              void operator=(int i) { decode(i); }
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

          // Optional fields
          std::unordered_map<std::string, std::string> aux;
          /**< Map the two char tag to value */
          std::unordered_map<std::string, char> aux_fmt; /**< Map the two char tag to format */

          /**
           * @brief
           * Add a aux field of form X:Y:Z, where X is the tag
           * Y is the format, Z is the value.
           * @param a aux field
           */
          void add(std::string a) {
              std::vector<std::string> s = split(a, ':');
              if (s.size() != 3) throw std::invalid_argument("Invalid aux field \"" + a + "\"");
              aux[s[0]] = s[2];
              aux_fmt[s[0]] = s[1].at(0);
          }

          /**
           * @brief
           * Add a char type aux tag.
           * @param tag key
           * @param val value
           */
          void add(std::string tag, char val) {
              aux[tag] = aux[std::string(1, val)];
              aux_fmt[tag] = 'A';
          }

          /**
           * @brief
           * Add a int type aux tag.
           * @param tag key
           * @param val value
           */
          void add(std::string tag, int val) {
              aux[tag] = aux[std::to_string(val)];
              aux_fmt[tag] = 'i';
          }

          /**
           * @brief
           * Add a float type aux tag.
           * @param tag key
           * @param val value
           */
          void add(std::string tag, float val) {
              aux[tag] = aux[std::to_string(val)];
              aux_fmt[tag] = 'f';
          }

          /**
           * @brief
           * Add a string type aux tag.
           * @param tag key
           * @param val value
           */
          void add(std::string tag, std::string val) {
              aux[tag] = aux[val];
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
           * Output the record in single line format
           * @return single line string
           */
          std::string to_string() const {
              std::stringstream ss;
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
                  << qual;
              for (auto &pair : aux) {
                  ss << '\t' << pair.first
                      << ':' << aux_fmt.at(pair.first)
                      << ':' << pair.second;
              }
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

              // Aux fields
              for (size_t i = 11; i < cols.size(); ++i) {
                  add(cols[i]);
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
  };

  class isam: public SAM {
    public:
      isam() { }
      isam(std::string file_name) : _file_name(file_name) { _init(); }
      ~isam() { _deinit(); }

      void open(std::string file_name) {
          _deinit();
          _file_name = file_name;
          _init();
      }

      void close() { _deinit(); }

      /**
       * @return true if file is open.
       */
      bool good() const {
          return in.good();
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
          if (!std::getline(in, _curr_line)) return false;
          _pprec << _curr_line;
          return true;
      }

      /**
       * @brief
       * First record is available after opening the file.
       * @return current Record
       */
      const Record &record() const { return _pprec; }

    protected:
      void _init() {
          in.open(_file_name);
          if (!in.good()) throw std::invalid_argument("Error opening file \"" + _file_name + "\"");
          std::stringstream hdr;
          while (std::getline(in, _curr_line) && _curr_line.at(0) == '@') {
              hdr << _curr_line << '\n';
          }
          _pphdr << hdr.str();
          _pprec << _curr_line;
      }

      void _deinit() {
          in.close();
          _pphdr = Header();
          _pprec = Record();
          _file_name = "";
      }

    private:
      std::string _file_name = "", _curr_line;
      std::ifstream in;

      Header _pphdr;
      Record _pprec;
  };

}

TEST_CASE ("SAM File") {
    {
        std::ofstream ss("tmp_s.sam");
        ss << "@HD\tVN:1.0\tSO:coordinate\n";
        ss
            << "@SQ\tSN:1\tLN:249250621\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta\tM5:1b22b98cdeb4a9304cb5d48026a85128\n";
        ss
            << "@SQ\tSN:2S\tLN:243199373\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta\tM5:a0d9851da00400dec1098a9255ac712e\n";
        ss
            << "@SQ\tSN:3Q\tLN:198022430\tAS:NCBI37\tUR:file:/data/local/ref/GATK/human_g1k_v37.fasta\tM5:fdfd811849cc2fadebc929bb925902e5\n";
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

    vargas::isam sf("tmp_s.sam");
    int x;
}
#endif //VARGAS_SAM_H
