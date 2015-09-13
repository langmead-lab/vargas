//
// Aligns reads using gssw, traceback is not possible.
// Created by Ravi Gaddipati (rgaddip1@jhu.edu).
//

#include "main.h"

using std::cout;
using std::endl;
using std::cerr;


int main(int argc, char *argv[]) {

  if (argc > 1) {
    if (!strcmp(argv[1], "build")) {
      build_main(argc, argv);
      exit(0);
    }
    else if (!strcmp(argv[1], "sim")) {
      sim_main(argc, argv);
      exit(0);
    }
    else if (!strcmp(argv[1], "align")) {
      align_main(argc, argv);
      exit(0);
    }
    else if (!strcmp(argv[1], "export")) {
      exportDOT(argc, argv);
      exit(0);
    }
  }
  GetOpt::GetOpt_pp args(argc, argv);
  if (args >> GetOpt::OptionPresent('h', "help")) {
    printMainHelp();
    exit(0);
  }
  cerr << "Error: Please define a valid mode of operation. See options with -h." << endl;
  exit(1);

}


void printMainHelp() {

  cout << endl << "---------------------- VMatch, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------" << endl;
  cout << "Operating modes \'vmatch MODE\':" << endl;
  cout << "build     Generate graph build file from reference and VCF files." << endl;
  cout << "sim       Simulate reads from a graph." << endl;
  cout << "align     Align reads to a graph." << endl;
  cout << "export    Export graph in DOT format." << endl << endl;

}