/**
 * Ravi Gaddipati
 * April 22, 2016
 * rgaddip1@jhu.edu
 *
 * Abstract class for objects that can be used as a
 * source of reads (i.e. ReadSim and ReadFile).
 *
 * readsource.h
 */

#ifndef VARGAS_READS_H
#define VARGAS_READS_H

#include <string>
#include <sstream>
#include "doctest/doctest.h"
#include "utils.h"
#include "simdpp/simd.h"

namespace vargas {

/**
 * Struct to represent a Read.
 * @param read base sequence.
 * @param end_pos position of last base in seq.
 * @param indiv Individual the read was taken from.
 * @param sub_err Number of substitiution errors introduced.
 * @param var_nodes Number of variant nodes the read traverses.
 * @param var_bases Number of bases that are in variant nodes.
 * @param indel_err Number of insertions and deletions introduced.
 */
  struct Read {
      Read() { }
      Read(std::string r) : read(r), read_num(seq_to_num(r)),
                            end_pos(-1), indiv(-1), sub_err(-1), var_nodes(-1), var_bases(-1), indel_err(-1) { }

      std::string read;
      std::vector<Base> read_num;
      int32_t end_pos;
      int32_t indiv;
      int32_t sub_err;
      int32_t var_nodes;
      int32_t var_bases;
      int32_t indel_err;

  };

  inline std::ostream &operator<<(std::ostream &os, const Read &r) {
      std::stringstream ss;
      ss << r.read
          << '#' << r.end_pos
          << ',' << r.indiv
          << ',' << r.sub_err
          << ',' << r.indel_err
          << ',' << r.var_nodes
          << ',' << r.var_bases;
      os << ss.str();
      return os;
  }

/**
 * Abstract class defining functions for read sources. A read source encapsulates
 * one read at a time. The stored read is updated with update_read(), and obtained
 * with get_read().
 */
  class ReadSource {

    public:

      ReadSource() { }
      virtual ~ReadSource() { }

      /**
       * Updates the stored and and returns the read.
       */
      virtual std::string update_and_get() {
          if (!update_read()) {
              read.read = "";
          }
          return to_string();
      }

      // Returns a string representation
      virtual std::string to_string() {
          std::stringstream ss;
          Read r = get_read();
          ss << r.read
              << '#' << r.end_pos
              << ',' << r.indiv
              << ',' << r.sub_err
              << ',' << r.indel_err
              << ',' << r.var_nodes
              << ',' << r.var_bases;
          return ss.str();
      };

      // Get the current read object
      virtual Read &get_read() = 0;

      // Get read file header
      virtual std::string get_header() const = 0;

      // Update the current read, return false if none are available
      virtual bool update_read() = 0;

      /**
       * Get size reads. If more reads are not available, a undersized
       * batch is returned.
       * @param size nominal number of reads to get.
       */
      const std::vector<Read> &get_batch(int size) {
          if (size <= 0) size = 1;
          _batch.clear();
          for (int i = 0; i < size; ++i) {
              if (!update_read()) break;
              _batch.push_back(read);
          }
          return _batch;
      }

      inline std::ostream &operator<<(std::ostream &os) {
          os << update_and_get();
          return os;
      }


    protected:
      Read read;
      std::string header = "";
      std::vector<Read> _batch;

  };

  /**
   * Container for a packaged batch of reads. Reads are interleaved so each SIMD vector
   * contains bases from all reads, respective to the base number. For example ReadBatch[0]
   * would contain the first bases of every read. Short reads or missing reads are padded
   * with Base::N.
   * @param read_len maximum length of the read
   * @param num_reads max number of reads. If a non-default T is used, this should be set to
   *    SIMDPP_FAST_T_SIZE where T corresponds to the width of T. For ex. Default T=simdpp::uint8 uses
   *    SIMDPP_FAST_INT8_SIZE
   * @param T element type
   */

  template<unsigned int read_len,
      unsigned int num_reads = SIMDPP_FAST_INT8_SIZE,
      template<unsigned int, typename=void> class T=simdpp::uint8>
  class ReadBatch {
    public:

      // Optimal number of elements in the vector for the current architecture
      //const static int num_reads = SIMDPP_FAST_INT8_SIZE;

      ReadBatch() { }

      /**
       * @param batch package the given vector of reads
       */
      ReadBatch(const std::vector<Read> &batch) : _reads(batch) { _package_reads(); }

      /**
       * @param obtain a batch of reads from the Read source and package them.
       */
      ReadBatch(ReadSource &rs) : _reads(rs.get_batch(num_reads)) {
          _package_reads();
      }

      /**
       * Load reads from a read source.
       * @param Read source to load a batch from.
       */
      bool load_reads(ReadSource &rs) {
          return load_reads(rs.get_batch(num_reads));
      }

      /**
       * @param bath load the given vector of reads.
       */
      bool load_reads(const std::vector<Read> &batch) {
          if (batch.size() == 0) return false;
          _reads = batch;
          _package_reads();
          return true;
      }

      /**
       * Return the i'th base of every read in a simdpp vector.
       * @param base index.
       */
      const T<num_reads> &at(int i) const {
          // let vector handle out of range errors
          return _packaged_reads.at(i);
      }

      /**
       * Non const version of at(i).
       * @param base index
       */
      T<num_reads> &operator[](int i) {
          return _packaged_reads.at(i);
      }

      /**
       * @return max read length. Echos template arg
       */
      size_t max_len() const { return read_len; }

      /**
       * Returns optimal number of reads in a batch based on SIMD architecture.
       * @return batch size.
       */
      size_t batch_size() const { return num_reads; }

      /**
       * @return Reads used to build the batch.
       */
      const std::vector<Read> &reads() const { return _reads; }

      /**
       * Get the utilization of the batch capacity. In effect how much
       * padding was used.
       * @return fill, between 0 and 1.
       */
      float fill() const {
          float f = 0;
          for (auto &r : _reads) f += r.read_num.size();
          return f / (num_reads * read_len);
      }

      typename std::vector<T<num_reads>>::const_iterator begin() const { return _packaged_reads.begin(); }
      typename std::vector<T<num_reads>>::const_iterator end() const { return _packaged_reads.end(); }

    private:

      /**
       * _packaged_reads[i] contains all i'th bases.
       * The length of _packaged_reads is the length of the read,
       * where as the length of _packaged_reads[i] is the number
       * of reads.
       */
      std::vector<T<num_reads>> _packaged_reads;

      // Unpackaged reads
      std::vector<Read> _reads;

      /**
       * Interleaves reads so all same-index base positions are in one
       * vector. Empty spaces are padded with Base::N.
       */
      inline void _package_reads() {
          _packaged_reads.resize(read_len);
          if (_reads.size() > num_reads) throw std::range_error("Too many reads for batch size.");

// allocate memory
          uchar **pckg = (uchar **) malloc(read_len * sizeof(uchar *));
          for (int i = 0; i < read_len; ++i) {
              pckg[i] = (uchar *) malloc(num_reads * sizeof(uchar));
          }

// Interleave reads
// For each read (read[i] is in _packaged_reads[0..n][i]
          for (size_t r = 0; r < _reads.size(); ++r) {

              if (_reads.at(r).read_num.size() > read_len) throw std::range_error("Read too long for batch size.");

// Put each base in the appropriate vector element
              for (size_t p = 0; p < _reads[r].read_num.size(); ++p) {
                  pckg[p][r] = _reads[r].read_num[p];
              }

// Pad the shorter reads
              for (size_t p = _reads[r].read_num.size(); p < read_len; ++p) {
                  pckg[p][r] = Base::N;
              }
          }

// Pad underful batches
          for (size_t r = _reads.size(); r < num_reads; ++r) {
              for (size_t p = 0; p < read_len; ++p) {
                  pckg[p][r] = Base::N;
              }
          }

// Load into vectors
          for (int i = 0; i < read_len; ++i) {
              _packaged_reads[i] = simdpp::load(pckg[i]);
          }

// Free memory
          for (int i = 0; i < read_len; ++i) {
              free(pckg[i]);
          }
          free(pckg);

      }

  };
}

TEST_CASE ("Read Batch") {
    std::vector<vargas::Read> reads;
    for (int i = 0; i < 15; ++i) {
        reads.push_back(vargas::Read("ACGTACGTCAGCCNNNCTAGTANCGTACTNGGCTAGAACGTACGTCAGCC"));
        }

        SUBCASE ("packaging") {
        vargas::ReadBatch<64> rb(reads);

            CHECK(rb.batch_size() == SIMDPP_FAST_INT8_SIZE);
            CHECK(rb.max_len() == 64);

        auto N = base_to_num('N');
        for (size_t i = 0; i < reads[0].read_num.size(); ++i) {
            auto b = rb[i];
            auto n = reads[0].read_num[i];
                CHECK(simdpp::extract<0>(b) == n);
                CHECK(simdpp::extract<1>(b) == n);
                CHECK(simdpp::extract<2>(b) == n);
                CHECK(simdpp::extract<3>(b) == n);
                CHECK(simdpp::extract<4>(b) == n);
                CHECK(simdpp::extract<5>(b) == n);
                CHECK(simdpp::extract<6>(b) == n);
                CHECK(simdpp::extract<7>(b) == n);
                CHECK(simdpp::extract<8>(b) == n);
                CHECK(simdpp::extract<9>(b) == n);
                CHECK(simdpp::extract<10>(b) == n);
                CHECK(simdpp::extract<11>(b) == n);
                CHECK(simdpp::extract<12>(b) == n);
                CHECK(simdpp::extract<13>(b) == n);
                CHECK(simdpp::extract<14>(b) == n);
                CHECK(simdpp::extract<15>(b) == N); // Since only 15 reads in the batch
        }

        // since read_len = 64 but len(read) is 50.
        for (size_t i = reads[0].read_num.size(); i < 64; ++i) {
            auto b = rb[i];
                CHECK(simdpp::extract<0>(b) == N);
                CHECK(simdpp::extract<1>(b) == N);
                CHECK(simdpp::extract<2>(b) == N);
                CHECK(simdpp::extract<3>(b) == N);
                CHECK(simdpp::extract<4>(b) == N);
                CHECK(simdpp::extract<5>(b) == N);
                CHECK(simdpp::extract<6>(b) == N);
                CHECK(simdpp::extract<7>(b) == N);
                CHECK(simdpp::extract<8>(b) == N);
                CHECK(simdpp::extract<9>(b) == N);
                CHECK(simdpp::extract<10>(b) == N);
                CHECK(simdpp::extract<11>(b) == N);
                CHECK(simdpp::extract<12>(b) == N);
                CHECK(simdpp::extract<13>(b) == N);
                CHECK(simdpp::extract<14>(b) == N);
                CHECK(simdpp::extract<15>(b) == N);
        }
    }

    std::remove("rds_tc.reads");

}

#endif //VARGAS_READS_H
