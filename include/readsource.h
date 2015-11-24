//
// Created by gaddra on 11/23/15.
//

#ifndef VMATCH_READS_H
#define VMATCH_READS_H

#include <string>

namespace vmatch {

class readsource {
 public:
  virtual std::string get() {
    updateRead();
    return getFull();
  }
  virtual std::string getFull() {
    return getRead() + "#" + getMeta();
  };
  virtual std::string getRead() = 0;
  virtual std::string getMeta() = 0;
  virtual bool updateRead() = 0;

  std::ostream &operator<<(std::ostream &os) {
    os << get();
    return os;
  }

 protected:
  std::string read;
  std::string meta;

};

}

#endif //VMATCH_READS_H
