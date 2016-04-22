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

#include "../include/getopt_pp.h"
#include "../include/graph.h"
#include "../include/readsim.h"
#include "../include/readfile.h"
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
      else if (!strcmp(argv[1], "stat")) {
        exit(stat_main(argc, argv));
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


int stat_main(const int argc, const char *argv[]) {
  std::string buildfile;

  GetOpt::GetOpt_pp args(argc, argv);

  if (args >> GetOpt::OptionPresent('h', "help")) {
    printStatHelp();
    return 0;
  }

  if (!(args >> GetOpt::Option('b', "buildfile", buildfile))) {
    std::cerr << "No buildfile specified" << std::endl;
    return 1;
  }

  std::ifstream bf(buildfile);
  if (!bf.good()) {
    std::cerr << "Error opening buildfile" << std::endl;
    return 1;
  }

  uint32_t numTotalNodes = 0, numVarNodes = 0, numEdges = 0;
  std::string line;
  std::vector<std::string> splitLine;

  while (getline(bf, line)) {
    if (line.at(0) == '#') continue;
    if (line.at(0) == ':') {
      numVarNodes++;
      continue;
    }
    split(line, ',', splitLine);

    switch (splitLine.size()) {
      case 2:
        numEdges++;
        break;
      case 3:
        numTotalNodes++;
        break;
      default:
        std::cerr << "Line split length of " << splitLine.size() << " unexpected." << std::endl << line << std::endl;
        break;
    }

  }

  std::cout << std::endl;
  std::cout << buildfile << " counts:" << std::endl;
  std::cout << "\tTotal number of nodes: " << numTotalNodes << std::endl;
  std::cout << "\tNumber of variant nodes: " << numVarNodes << std::endl;
  std::cout << "\tTotal number of edges: " << numEdges << std::endl;
  std::cout << std::endl;

  return 0;
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
      >> GetOpt::OptionPresent('m', "maxref", p.maxAF)
      >> GetOpt::Option('e', "exref", p.includeRefIndivs)
      >> GetOpt::Option('m', "loadref", p.inMemoryRef);

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

  vargas::Graph g;
  if (args >> GetOpt::Option('b', "build", buildfile)) {
    g.buildGraph(buildfile);
  }
  else {
    std::cerr << "Error: Buildfile required." << std::endl;
    return 1;
  }

  std::string inputAligns = "";
  if (args >> GetOpt::Option('c', "context", inputAligns)) {
    std::string line;
    std::ifstream input(inputAligns);
    if (!input.good()) {
      throw std::invalid_argument("Invalid file: " + inputAligns);
    }

    while (std::getline(input, line)) {
      if (line.length() == 0) continue;
      vargas::Alignment a(line);
      vargas::Graph contextGraph(g, a);
      contextGraph.exportDOT(std::cout, line);
    }
  }

    // Export whole graph
  else {
    g.exportDOT(std::cout);
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
  std::string alignOutFile = "";
  args >> GetOpt::Option('m', "match", p.match)
      >> GetOpt::Option('n', "mismatch", p.mismatch)
      >> GetOpt::Option('o', "gap_open", p.gap_open)
      >> GetOpt::Option('e', "gap_extend", p.gap_extension)
      >> GetOpt::Option('f', "outfile", alignOutFile);

  // Set output
  bool useFile = false;
  std::ofstream aOutStream;
  if (alignOutFile.length() != 0) {
    aOutStream.open(alignOutFile);
    useFile = true;
  }

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

  // Align until there are no more reads
  vargas::Alignment align;
  while (reads.updateRead()) {
    g.align(reads.getRead(), align);
    if (useFile) {
      aOutStream << align << std::endl;
    } else {
      std::cout << align << std::endl;
    }
  }

  if (useFile) {
    aOutStream << reads.getHeader() << std::endl;
  } else {
    std::cout << reads.getHeader() << std::endl;
  }

  return 0;
}


int sim_main(const int argc, const char *argv[]) {

  //TODO this can be made more efficient by creating all lists of ingroups first, rather than parsing the input every single time

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
      >> GetOpt::Option('l', "readlen", p.readLen)
      >> GetOpt::Option('a', "ambiguity", p.ambiguity);

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
      prof.numSubErr = (splitProf[0] == "*") ? -1 : std::stoi(splitProf[0]);
      prof.numVarNodes = (splitProf[1] == "*") ? -1 : std::stoi(splitProf[1]);
      prof.numVarBases = (splitProf[2] == "*") ? -1 : std::stoi(splitProf[2]);
      prof.numIndelErr = (splitProf[3] == "*") ? -1 : std::stoi(splitProf[3]);
      sim.addProfile(prof, prefix + std::to_string(i) + ".reads");
    }
    sim.populateProfiles();

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
  cout << "\tstat      Count nodes and edges of a given graph." << endl;
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
  cout << "-e\t--exref         Exclude the list of individuals from the reference alleles." << endl;
  cout << "-s\t--set           <#,#,..,#> Generate a buildfile for a list of ingroup percents." << endl;
  cout << "-c\t--complement    <graph.build> Generate a complement of all graphs in -s" << endl;
  cout << "-m\t--loadref       Load entire reference into memory" << endl;


  cout << endl << "--maxref is applied after ingroup filter" << endl;
  cout << "Buildfile is output to [s][In/Out].build" << endl << endl;
}


void printExportHelp() {
  using std::cout;
  using std::endl;
  cout << endl
      << "------------------ vargas export, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile    (required) Graph to export to DOT." << endl;
  cout << "-c\t--context      [InputFile] Export the context of alignments." << endl;

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
  cout << "-l\t--readlen       Read length, strictly enforced." << endl;
  cout << "-e\t--profile       <p1 p2 .. p3> Match read profiles, space delimited. Produces -n of each." << endl;
  cout << "-p\t--prefix        Prefix to use for read files generated with -e" << endl;
  cout << "-r\t--randwalk      Random walk, read may change individuals at branches." << endl;
  cout << "-a\t--ambiguity     Max number of ambiguous bases to allow in reads" << endl << endl;

  cout << "Read Profile format (use \'*\' for any): " << endl;
  cout << "\tnumSubErr,numVarNodes,numVarBases,numIndelErr" << endl;
  cout << "\tExample: Any read with 1 sub error and 1 variant node." << endl;
  cout << "\t\t1,1,*,*" << endl;
  cout << "Read Format:" << endl;
  cout << "\tREAD#READ_END_POSITION,INDIVIDUAL,NUM_SUB_ERR,NUM_INDEL_ERR,NUM_VAR_NODE,NUM_VAR_BASES" << endl << endl;
}


void printAlignHelp() {
  using std::cout;
  using std::endl;

  vargas::Graph::GraphParams p;

  cout << endl
      << "------------------- vargas align, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     Quick rebuild file." << endl;
  cout << "-m\t--match         Match score, default " << int(p.match) << endl;
  cout << "-n\t--mismatch      Mismatch score, default " << int(p.mismatch) << endl;
  cout << "-o\t--gap_open      Gap opening penalty, default " << int(p.gap_open) << endl;
  cout << "-e\t--gap_extend    Gap extend penalty, default " << int(p.gap_extension) << endl;
  cout << "-r\t--reads         Reads to align. Use stdin if not defined." << endl;
  cout << "-f\t--outfile       Alignment output file. If not defined, use stdout." << endl;

  cout << "Lines beginning with \'#\' are ignored." << endl;
  cout << "Output format:" << endl;
  cout << "\tREAD,OPTIMAL_SCORE,OPTIMAL_ALIGNMENT_END,NUM_OPTIMAL_ALIGNMENTS,SUBOPTIMAL_SCORE," << endl;
  cout << "\tSUBOPTIMAL_ALIGNMENT_END,NUM_SUBOPTIMAL_ALIGNMENTS,ALIGNMENT_MATCH" << endl << endl;
  cout << "ALIGNMENT_MATCH:\n\t0- optimal match, 1- suboptimal match, 2- no match" << endl << endl;

}


void printStatHelp() {
  using std::cout;
  using std::endl;

  cout << endl
      << "------------------- vargas stat, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     Quick rebuild file." << endl << endl;

}