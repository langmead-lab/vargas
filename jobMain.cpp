//
// Created by gaddra on 10/2/15.
//

#include <iostream>
#include <sstream>
#include "jobMain.h"
#include "getopt_pp.h"
#include "utils.h"

using std::cout;
using std::endl;

int job_main(int argc, char *argv[]) {
  using std::string;
  using GetOpt::Option;

  string refFile = "", vcfFile = "", ingroups = ".*", numvarnodes = ".*", numvarbase = ".*", numsuberr = ".*";
  std::vector<std::string> ingroupsSplit(0), numvarnodesSplit(0), numvarbaseSplit(0), numsuberrSplit(0);
  int32_t numReads = 0;
  std::ofstream out;
  std::streambuf *coutbuf;
  std::stringstream tmp;


  GetOpt::GetOpt_pp args(argc, argv);

  if ((args >> GetOpt::OptionPresent('h', "help"))) {
    printJobHelp();
    return 0;
  }

  args >> Option('r', "ref", refFile)
      >> Option('v', "vcf", vcfFile)
      >> Option('i', "ingroups", ingroups)
      >> Option('n', "numreads", numReads)
      >> Option('a', "numvarnodes", numvarnodes)
      >> Option('s', "numsubserr", numsuberr)
      >> Option('g', "numvarbase", numvarbase);

  split(ingroups, ',', ingroupsSplit);
  split(numvarnodes, ',', numvarnodesSplit);
  split(numvarbase, ',', numvarbaseSplit);
  split(numsuberr, ',', numsuberrSplit);

  for (auto ig : ingroupsSplit) {
    tmp << ig << "inGroup.build";
//    out.open(tmp);
    tmp.clear();
    std::cout.rdbuf(out.rdbuf());
    for (auto nvn : numvarnodesSplit) {
      for (auto nvb : numvarbaseSplit) {
        for (auto nvs : numsuberrSplit) {
        }
      }
    }
  }


}

template<typename t>
std::string getRegex(t subErr, t numVarNode, t numVarBase) {
  std::stringstream ss;
  ss << "^.*#.*,.*," << subErr << "," << numVarNode << "," << numVarBase;
  return ss.str();
}

void printJobHelp() {

  cout << endl << "---------------------- VMatch job, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------"
      << endl;
  cout << "Export a PBS job for the specified run profile. Default is any." << endl << endl;
  cout << "-r\t--ref           Reference FASTA" << endl;
  cout << "-v\t--vcf           Reference VCF" << endl;
  cout << "-n\t--numreads      Number of reads per combination of conditions" << endl;
  cout << "-i\t--ingroups      [0-100],...,[0-100] Set of ingroup percentage graphs to generate" << endl;
  cout << "-a\t--numvarnodes   [0-10],...,[0-10] Set of variant nodes per read." << endl;
  cout << "-s\t--numsuberr     [0-10],...,[0-10] Set of substitution errors per read." << endl;
  cout << "-g\t--numvarbase    [0-100],...,[0-100] Set of variant bases per read" << endl;
}