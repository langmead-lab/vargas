//
// Created by gaddra on 6/17/16.
//

#ifndef VARGAS_SIM_H
#define VARGAS_SIM_H

#include "readsource.h"
#include "graph.h"

namespace vargas {

  /**
   * Paramter list controlling the types of reads created.
   * @param len Nominal length of the read
   * @param enforce_len reject reads that do not match len
   * @param num_mut Number of mutation errors
   * @param num_indel number of insertions/deletions
   * @param rand Introduce mutations and indels at a random rate
   * @param mut_rate mutation rate, used with rand=true
   * @param indel_rate indel generation rate, used with rand=true
   * @param real Follow a single individual throughout reads
   * @param ambig allow reads with a significant (>25%) ambigous bases
   */
  struct ReadProfile {
      int len = 50;
      int num_mut = 0;
      int num_indel = 0;

      bool rand = false;
      float mut_rate = 0.02f;
      float indel_rate = 0.02f;

      bool real = true;
      bool ambig = false;
      bool enforce_len = false;
  };

  inline std::ostream &operator<<(std::ostream &os, const ReadProfile &rp) {
      os << "len=" << rp.len
          << " num_mut=" << rp.num_mut
          << " num_indel=" << rp.num_indel
          << " rand=" << rp.rand
          << " mut_rate=" << rp.mut_rate
          << " indel_rate=" << rp.indel_rate
          << " real=" << rp.real;
      return os;
  }

  class ReadSim: public ReadSource {

    public:
      ReadSim(const Graph &g) : _graph(g) { }

      ReadSim(const Graph &_graph, const ReadProfile &prof) : _graph(_graph), _prof(prof) { }

      /**
       * generate and store an updated read.
       * @return true if successful
       */
      virtual bool update_read() override {
          auto &nodes = *_graph.node_map();
          auto &next = _graph.next_map();

      }

      /**
       * Get the profile being used to generate reads (if use_prof)
       * @return ReadProfile
       */
      const ReadProfile &prof() const {
          return _prof;
      }

      /**
       * Crete reads following prof as a template
       * @param prof Read Profile
       */
      void set_prof(const ReadProfile &prof) {
          _prof = prof;
      }

      virtual std::string get_header() const override {
          std::stringstream ss;
          ss << _prof;
          return ss.str();
      }

      std::ostream &operator<<(std::ostream &os) const {
          os << '>' << _prof << std::endl << read;
          return os;
      }


    private:
      const vargas::Graph &_graph;
      ReadProfile _prof;

  };

}

#endif //VARGAS_SIM_H
