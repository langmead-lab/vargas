//
// Created by gaddra on 10/2/15.
//

#include <iostream>
#include "jobMain.h"

using std::cout;
using std::endl;

int job_main(int argc, char *argv[]) {

}

void printJobHelp() {

  cout << endl << "---------------------- VMatch, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------" << endl;
  cout << "Operating modes \'vmatch MODE\':" << endl;
  cout << "build     Generate graph build file from reference and VCF files." << endl;
  cout << "sim       Simulate reads from a graph." << endl;
  cout << "align     Align reads to a graph." << endl;
  cout << "export    Export graph in DOT format." << endl << endl;

}