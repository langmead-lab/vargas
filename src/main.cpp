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
    else if (!strcmp(argv[1], "sim")) {
      exit(sim_main(argc, argv));
    }
    else if (!strcmp(argv[1], "align")) {
      exit(align_main(argc, argv));
    }
    else if (!strcmp(argv[1], "export")) {
      exit(export_main(argc, argv));
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
  bool makeComplements = false;

  vargas::Graph g;

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

  vargas::Graph::GraphParams p = g.getParamsCopy();

  args >> GetOpt::Option('l', "maxlen", p.maxNodeLen)
      >> GetOpt::Option('R', "region", p.region)
      >> GetOpt::OptionPresent('c', "complement", makeComplements)
      >> GetOpt::Option('c', "complement", p.complementSource)
      >> GetOpt::OptionPresent('m', "maxref", p.maxAF);

  g.setParams(p);

  std::vector<std::string> setSplit(0);
  if (args >> GetOpt::Option('s', "set", setString)) {
    split(setString, ',', setSplit);
  } else setSplit.push_back("100"); // Default includes everything

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
      if (makeComplements) {
        g.setComplementSource(fileName.str());
        fileName.str(std::string());
        fileName << ingrp << "Out.build";
        g.setComplement(true);
        g.exportBuildfile(REF, VCF, fileName.str());
      }

    }

  return 0;
}

int export_main(int argc, char *argv[]) {

  GetOpt::GetOpt_pp args(argc, argv);

  if (args >> GetOpt::OptionPresent('h', "help")) {
    printExportHelp();
    return 0;
  }

  std::string buildfile;

  if (args >> GetOpt::Option('b', "build", buildfile)) {
    vargas::Graph g(buildfile);
    g.exportDOT(std::cout);
  }
  else {
    std::cerr << "Error: Buildfile required." << std::endl;
    return 1;
  }
  return 0;
}

int align_main(int argc, char *argv[]) { }

int sim_main(int argc, char *argv[]) { }

void printMainHelp() {
  using std::cout;
  using std::endl;
  cout << endl << "---------------------- VMatch, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------" << endl;
  cout << "Operating modes \'vargas MODE\':" << endl;
  cout << "\tbuild     Generate graph build file from reference and VCF files." << endl;
  cout << "\tsim       Simulate reads from a graph." << endl;
  cout << "\talign     Align reads to a graph." << endl;
  cout << "\texport    Export graph in DOT format." << endl << endl;
}
void printBuildHelp() {
  using std::cout;
  using std::endl;
  cout << endl << "------------------- VMatch build, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" <<
      endl;
  cout << "-v\t--vcf           (required) VCF file, uncompressed." << endl;
  cout << "-r\t--ref           (required) reference, single record FASTA" << endl;
  cout << "-l\t--maxlen        Maximum node length" << endl;
  cout << "-R\t--region        <min:max> Ref region, inclusive. Default is the entire graph." << endl;
  cout << "-m\t--maxref        Generate a graph using alleles in the ingroup w/ the highest frequency." << endl;
  cout << "-s\t--set           <#,#,..,#> Generate a buildfile for a list of ingroup %'s and their complements."
      << endl;
  cout << "-c\t--complement    <graph.build> Generate a complement of all graphs in -s" << endl;

  cout << endl << "Buildfile is output to [s][In/Out].build" << endl << endl;
}
void printExportHelp() {
  std::cout << std::endl << "------------------- VMatch export, " << __DATE__
      << ". rgaddip1@jhu.edu -------------------" << std::endl;
  std::cout << "-b\t--buildfile    (required) Graph to export to DOT." << std::endl;
  std::cout << std::endl << "DOT file printed to stdout." << std::endl << std::endl;
}
void printSimHelp() {
  using std::cout;
  using std::endl;
  cout << endl << "------------------- VMatch sim, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     quick rebuild file, generate with vargas build" << endl;
  cout << "-n\t--numreads      Number of reads to simulate" << endl;
  cout << "-m\t--muterr        Simulated read mutation error rate" << endl;
  cout << "-i\t--indelerr      Simulated read Indel error rate" << endl;
  cout << "-r\t--rand          Use a random mutation error rate, up to the value specified by -m." << endl;
  cout << "-l\t--readlen       Nominal read length" << endl;
  cout << "-e\t--regex         Match regex expressions. Produces -n of each, discard others." << endl;
  cout << "-f\t--filesplit     Each regex put in seperate files, with prefix specified. Otherwise all to stdout"
      << endl;
  cout << "  \t                List of expressions is space delimited -e \"exp1 exp2\"." << endl << endl;

  cout << "NOTE: End of line anchor may not work in regex depending on C++ version. " << endl;
  cout << "Reads are printed on stdout." << endl;
  cout << "Read Format:" << endl;
  cout << "READ#READ_END_POSITION,INDIVIDUAL,NUM_SUB_ERR,NUM_VAR_NODE,NUM_VAR_BASES" << endl << endl;
}
void printAlignHelp() {
  using std::cout;
  using std::endl;
  cout << endl << "------------------- VMatch align, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     quick rebuild file, required if -v, -r are not defined." << endl;
  cout << "-m\t--match         Match score, default  " << 2 << endl;
  cout << "-n\t--mismatch      Mismatch score, default " << 2 << endl;
  cout << "-o\t--gap_open      Gap opening score, default " << 3 << endl;
  cout << "-e\t--gap_extend    Gap extend score, default " << 1 << endl;
  cout << "-r\t--reads         Reads to align. If not specified, read from stdin" << endl;

  cout << endl << "Alignments output to stdout. Reads read from stdin or -r, 1 per line." << endl;
  cout << "Lines beginning with \'#\' are ignored." << endl;
  cout << "Output format:" << endl;
  cout << "READ,OPTIMAL_SCORE,OPTIMAL_ALIGNMENT_END,NUM_OPTIMAL_ALIGNMENTS,SUBOPTIMAL_SCORE," << endl;
  cout << "SUBOPTIMAL_ALIGNMENT_END,NUM_SUBOPTIMAL_ALIGNMENTS,ALIGNMENT_MATCH" << endl << endl;
  cout << "ALIGNMENT_MATCH: 0- optimal match, 1- suboptimal match, 2- no match" << endl << endl;

}