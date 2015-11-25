//
// Aligns reads using gssw, traceback is not possible.
// Created by Ravi Gaddipati (rgaddip1@jhu.edu).
//

#include <string.h>
#include "../include/main.h"


int main(int argc, char *argv[]) {

  if (argc > 1) {
    if (!strcmp(argv[1], "build")) {
      exit(build_main(argc, argv));
    }
    else if (!strcmp(argv[1], "sim")) { ;//      exit(sim_main(argc, argv));
    }
    else if (!strcmp(argv[1], "align")) { ;//     exit(align_main(argc, argv));
    }
    else if (!strcmp(argv[1], "export")) { ;//      exit(exportDOT(argc, argv));
    }
  }
  GetOpt::GetOpt_pp args(argc, argv);
  if (args >> GetOpt::OptionPresent('h', "help")) {
    printMainHelp();
    exit(0);
  }
  std::cerr << "Error: Please define a valid mode of operation. See options with -h." << std::endl;
  exit(1);

}

int build_main(int argc, char *argv[]) {
  std::string VCF = "", REF = "";
  std::string setString;

  vmatch::Graph g;

  GetOpt::GetOpt_pp args(argc, argv);

  if (args >> GetOpt::OptionPresent('h', "help")) {
    printBuildHelp();
    return 0;
  }

  if (!(args >> GetOpt::Option('v', "vcf", VCF))
      || !(args >> GetOpt::Option('r', "ref", REF))) {
    std::cerr << "No input files specified!" << std::endl;
    return 1;
  }

  vmatch::Graph::GraphParams p = g.getParamsCopy();

  args >> GetOpt::Option('l', "maxlen", p.maxNodeLen)
      >> GetOpt::Option('R', "region", p.region)
      >> GetOpt::Option('g', "ingroup", p.ingroup)
      >> GetOpt::OptionPresent('c', "complement", p.genComplement)
      >> GetOpt::Option('c', "complement", p.complementSource)
      >> GetOpt::OptionPresent('m', "maxref", p.maxAF);

  g.setParams(p);

  if (args >> GetOpt::Option('s', "set", setString)) {
    std::vector<std::string> setSplit(0);
    split(setString, ',', setSplit);
    // Gen a set of ingroup build files
    std::stringstream fileName;

    for (auto &ingrp : setSplit) {
      // Generates a buildfile for all % ingroup specified, as well as the outgroup buildfiles
      g.setIngroup(stoi(ingrp));

      // Ingroup
      fileName.str(std::string());
      fileName << ingrp << "In.build";
      g.setComplement(false);
      g.exportBuildfile(REF, VCF, fileName.str());

      // Outgroup
      g.setComplementSource(fileName.str());
      fileName.str(std::string());
      fileName << ingrp << "Out.build";
      g.setComplement(true);
      g.exportBuildfile(REF, VCF, fileName.str());
    }

  } else g.exportBuildfile(REF, VCF);

  return 0;
}

void printMainHelp() {
  using std::cout;
  using std::endl;
  cout << endl << "---------------------- VMatch, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------" << endl;
  cout << "Operating modes \'vmatch MODE\':" << endl;
  cout << "build     Generate graph build file from reference and VCF files." << endl;
  cout << "sim       Simulate reads from a graph." << endl;
  cout << "align     Align reads to a graph." << endl;
  cout << "export    Export graph in DOT format." << endl << endl;
}
void printBuildHelp() {
  using std::cout;
  using std::endl;
  cout << endl << "------------------- VMatch build, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" <<
      endl;
  cout << "-v\t--vcf           (required) VCF file, uncompressed." << endl;
  cout << "-r\t--ref           (required) reference single record FASTA" << endl;
  cout << "-l\t--maxlen        Maximum node length" << endl;
  cout << "-R\t--region        [min:max] Ref region, inclusive. Default is entire graph." << endl;
  cout << "-g\t--ingroup       Percent of individuals to build graph from, default all." << endl;
  cout << "-c\t--complement    Generate a complement of the specified graph" << endl;
  cout << "-s\t--set           <#,#,..,#> Generate a buildfile for a list of ingroup %'s and their complements."
      << endl;
  cout << "-m\t--maxref        Generate a graph using allele's w/ the highest frequency. Overrides other opts." << endl;
  cout << "\t                  -s outputs to files." << endl;

  cout << endl << "Buildfile is printed on stdout." << endl;
}