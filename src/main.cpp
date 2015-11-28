/**
 * Ravi Gaddipati
 * November 25, 2015
 * rgaddip1@jhu.edu
 *
 * Interface for simulating and aligning reads from/to a DAG.
 * Uses a modified gssw from Erik Garrison.
 *
 * main.cpp
 */

#include "../include/main.h"


int main(const int argc, const char *argv[]) {

  try {
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
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    exit(1);
  }

  GetOpt::GetOpt_pp args(argc, argv);
  if (args >> GetOpt::OptionPresent('h', "help")) {
    printMainHelp();
    exit(0);
  }
  std::cerr << "Define a valid mode of operation." << std::endl;
  printMainHelp();
  exit(1);

}


int build_main(const int argc, const char *argv[]) {
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

  vargas::Graph::GraphParams p;

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
  for (auto ingrp : setSplit) {
      // Generates a buildfile for all % ingroup specified, as well as the outgroup buildfiles
      g.setIngroup(stoi(ingrp));

      // Ingroup
      fileName.str(std::string());
      fileName << ingrp << "In.build";
      g.exportBuildfile(REF, VCF, fileName.str());

      // Outgroup
      if (makeComplements) {
        g.setComplementSource(fileName.str());
        fileName.str(std::string());
        fileName << ingrp << "Out.build";
        g.setComplement(true);
        g.exportBuildfile(REF, VCF, fileName.str());
        g.setComplement(false);
      }

    }

  return 0;
}


int export_main(const int argc, const char *argv[]) {

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


int align_main(const int argc, const char *argv[]) {
  vargas::Graph g;
  GetOpt::GetOpt_pp args(argc, argv);

  if (args >> GetOpt::OptionPresent('h', "help")) {
    printAlignHelp();
    return 0;
  }

  // Load parameters
  vargas::Graph::GraphParams p;
  args >> GetOpt::Option('m', "match", p.match)
      >> GetOpt::Option('n', "mismatch", p.mismatch)
      >> GetOpt::Option('o', "gap_open", p.gap_open)
      >> GetOpt::Option('e', "gap_extend", p.gap_extension);
  g.setParams(p);

  // Read file handler
  std::string readsfile;
  if (!(args >> GetOpt::Option('r', "reads", readsfile))) {
    std::cerr << "No buildfile defined!" << std::endl;
    return 1;
  }
  vargas::ReadFile reads(readsfile);


  std::string buildfile;
  if (!(args >> GetOpt::Option('b', "build", buildfile))) {
    std::cerr << "No buildfile defined!" << std::endl;
    return 1;
  }
  g.buildGraph(buildfile);

  // Align until there are not more reads
  vargas::Alignment align;
  while (reads.updateRead()) {
    g.align(reads.getRead(), align);
    std::cout << align << std::endl;
  }
  std::cout << reads.getHeader() << std::endl;

  return 0;
}


int sim_main(const int argc, const char *argv[]) {
  GetOpt::GetOpt_pp args(argc, argv);

  if (args >> GetOpt::OptionPresent('h', "help")) {
    printSimHelp();
    return 0;
  }

  vargas::SimParams p;
  args >> GetOpt::Option('n', "numreads", p.maxreads)
      >> GetOpt::Option('r', "randwalk", p.randWalk)
      >> GetOpt::Option('m', "muterr", p.muterr)
      >> GetOpt::Option('i', "indelerr", p.indelerr)
      >> GetOpt::Option('l', "readlen", p.readLen);

  std::string buildfile;
  if (!(args >> GetOpt::Option('b', "buildfile", buildfile))) throw std::invalid_argument("Buildfile required.");
  vargas::Graph g;
  g.useIndividuals(!p.randWalk); // Don't need individuals if we're doing a random walk
  g.buildGraph(buildfile);

  vargas::ReadSim sim(p);
  sim.setGraph(g);

  std::string profiles;
  if (args >> GetOpt::Option('e', "profile", profiles)) {
    std::string prefix = "sim";
    args >> GetOpt::Option('p', "prefix", prefix);

    std::vector<std::string> splitProfiles = split(profiles, ' ');
    std::vector<std::string> splitProf;
    vargas::ReadProfile prof;
    for (int i = 0; i < splitProfiles.size(); ++i) {
      split(splitProfiles[i], ',', splitProf);
      if (splitProf.size() != 4) throw std::invalid_argument("Profile must have 4 fields (" + splitProfiles[i] + ").");
      prof.indiv = (splitProf[0] == "*") ? -1 : std::stoi(splitProf[0]);
      prof.numSubErr = (splitProf[1] == "*") ? -1 : std::stoi(splitProf[1]);
      prof.numVarNodes = (splitProf[2] == "*") ? -1 : std::stoi(splitProf[2]);
      prof.numVarBases = (splitProf[3] == "*") ? -1 : std::stoi(splitProf[3]);
      sim.addProfile(prof, prefix + std::to_string(i) + ".reads");
      while(sim.updateRead());
    }

  } else {
    for (int i = 0; i < p.maxreads; ++i) {
      std::cout << sim.updateAndGet() << std::endl;
    }
  }

  return 0;
}


void printMainHelp() {
  using std::cout;
  using std::endl;
  cout << endl
      << "---------------------- vargas, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------" << endl;
  cout << "Operating modes \'vargas MODE\':" << endl;
  cout << "\tbuild     Generate graph build file from reference and VCF files." << endl;
  cout << "\tsim       Simulate reads from a graph." << endl;
  cout << "\talign     Align reads to a graph." << endl;
  cout << "\texport    Export graph in DOT format." << endl << endl;
}


void printBuildHelp() {
  using std::cout;
  using std::endl;
  cout << endl
      << "------------------- vargas build, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-v\t--vcf           (required) VCF file, uncompressed." << endl;
  cout << "-r\t--ref           (required) reference, single record FASTA" << endl;
  cout << "-l\t--maxlen        Maximum node length" << endl;
  cout << "-R\t--region        <min:max> Ref region, inclusive. Default is the entire graph." << endl;
  cout << "-m\t--maxref        Generate a graph using alleles in the ingroup w/ the highest frequency." << endl;
  cout << "-s\t--set           <#,#,..,#> Generate a buildfile for a list of ingroup %'s and their complements."
      << endl;
  cout << "-c\t--complement    <graph.build> Generate a complement of all graphs in -s" << endl;


  cout << endl << "--maxref is applied after ingroup filter" << endl;
  cout << "Buildfile is output to [s][In/Out].build" << endl << endl;
}


void printExportHelp() {
  using std::cout;
  using std::endl;
  cout << endl
      << "------------------ vargas export, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile    (required) Graph to export to DOT." << endl;
  cout << endl << "DOT file printed to stdout." << endl << endl;
}


void printSimHelp() {
  using std::cout;
  using std::endl;
  cout << endl
      << "-------------------- vargas sim, " << __DATE__ << ". rgaddip1@jhu.edu --------------------" << endl;
  cout << "-b\t--buildfile     quick rebuild file, generate with vargas build" << endl;
  cout << "-n\t--numreads      Number of reads to simulate" << endl;
  cout << "-m\t--muterr        Simulated read mutation error rate" << endl;
  cout << "-i\t--indelerr      Simulated read Indel error rate" << endl;
  cout << "-l\t--readlen       Nominal read length" << endl;
  cout << "-e\t--profile       <p1 p2 .. p3> Match read profiles, space delimited. Produces -n of each." << endl;
  cout << "-p\t--prefix        Prefix to use for read files generated with -e" << endl;
  cout << "-r\t--randwalk      Random walk, read may change individuals at branches." << endl << endl;

  cout << "Read Profile format (use \'*\' for any): " << endl;
  cout << "\tindiv,numSubErr,numVarNodes,numVarBases" << endl;
  cout << "Read Format:" << endl;
  cout << "\tREAD#READ_END_POSITION,INDIVIDUAL,NUM_SUB_ERR,NUM_VAR_NODE,NUM_VAR_BASES" << endl << endl;
}


void printAlignHelp() {
  using std::cout;
  using std::endl;
  cout << endl
      << "------------------- vargas align, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     Quick rebuild file." << endl;
  cout << "-m\t--match         Match score, default  " << 2 << endl;
  cout << "-n\t--mismatch      Mismatch score, default " << 2 << endl;
  cout << "-o\t--gap_open      Gap opening score, default " << 3 << endl;
  cout << "-e\t--gap_extend    Gap extend score, default " << 1 << endl;
  cout << "-r\t--reads         Reads to align." << endl;

  cout << endl << "Alignments output to stdout. Reads read from stdin or -r, 1 per line." << endl;
  cout << "Lines beginning with \'#\' are ignored." << endl;
  cout << "Output format:" << endl;
  cout << "\tREAD,OPTIMAL_SCORE,OPTIMAL_ALIGNMENT_END,NUM_OPTIMAL_ALIGNMENTS,SUBOPTIMAL_SCORE," << endl;
  cout << "\tSUBOPTIMAL_ALIGNMENT_END,NUM_SUBOPTIMAL_ALIGNMENTS,ALIGNMENT_MATCH" << endl << endl;
  cout << "ALIGNMENT_MATCH:\n\t0- optimal match, 1- suboptimal match, 2- no match" << endl << endl;

}