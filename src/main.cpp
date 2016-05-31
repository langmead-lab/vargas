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
#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/doctest.h"
#include "getopt_pp.h"
#include "graph.h"
#include "alignment.h"
//#include "../include/readsim.h"


int main(const int argc, const char *argv[]) {

  std::clock_t start;
  std::cout << "Initial Build: ";
  start = std::clock();

  vargas::GraphBuilder gb("hs37d5_22.fa", "chr22.bcf");
  gb.region("22:22,000,000-30,000,000");
  gb.ingroup(100);
  vargas::Graph g;
  gb.build(g);
  std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
  std::cout << g.next_map().size() << ":" << (*(g.begin())).end() << ", " << (*(--g.end())).end() << std::endl;

  srand(time(NULL));
  std::vector<bool> filter;
  for (int i = 0; i < g.pop_size(); ++i) filter.push_back(rand() % 100 > 50);
  std::cout << "FILTER: ";
  start = std::clock();
  vargas::Graph g1(g, filter);
  std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
  std::cout << g1.next_map().size() << ":" << (*(g.begin())).end() << ", " << (*(--g.end())).end() << std::endl;

  std::cout << "REF: ";
  start = std::clock();
  vargas::Graph g2(g, vargas::Graph::REF);
  std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
  std::cout << g2.next_map().size() << ":" << (*(g.begin())).end() << ", " << (*(--g.end())).end() << std::endl;

  std::cout << "MAXAF: ";
  start = std::clock();
  vargas::Graph g3(g, vargas::Graph::MAXAF);
  std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
  std::cout << g3.next_map().size() << ":" << (*(g.begin())).end() << ", " << (*(--g.end())).end() << std::endl;

  return 0;

  GetOpt::GetOpt_pp args(argc, argv);

  // Run test cases
  if (args >> GetOpt::OptionPresent("runtests")) {
    doctest::Context doc(argc, argv);
    doc.setOption("no-breaks", true);
    doc.setOption("abort-after", 5);
    doc.setOption("sort", "name");
    int res = doc.run();
    if (doc.shouldExit()) return res;
    std::cerr << std::endl;
  }

  try {
    if (argc > 1) {
      if (!strcmp(argv[1], "build")) {
//        exit(build_main(argc, argv));
      }
      else if (!strcmp(argv[1], "sim")) {
//        exit(sim_main(argc, argv));
      }
      else if (!strcmp(argv[1], "align")) {
//        exit(align_main(argc, argv));
      }
      else if (!strcmp(argv[1], "export")) {
        //       exit(export_main(argc, argv));
      }
      else if (!strcmp(argv[1], "stat")) {
        //       exit(stat_main(argc, argv));
      }
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    exit(1);
  }

  if (args >> GetOpt::OptionPresent('h', "help")) {
//    printMainHelp();
    exit(0);
  }
  std::cerr << "Define a valid mode of operation." << std::endl;
//  printMainHelp();
  exit(1);

}

/**
int stat_main(const int argc, const char *argv[]) {
  std::string buildfile;

  GetOpt::GetOpt_pp args(argc, argv);

  if (args >> GetOpt::OptionPresent('h', "help")) {
    printStatHelp();
    return 0;
  }

  if (!(args >> GetOpt::Option('b', "buildfile", buildfile))) {
    printStatHelp();
    throw std::invalid_argument("No buildfile specified!");
  }

  std::ifstream bf(buildfile);
  if (!bf.good()) {
    throw std::invalid_argument("Error opening buildfile: " + buildfile);
  }

  uint32_t numTotalNodes = 0, var_nodes = 0, numEdges = 0;
  std::string line;
  std::vector<std::string> splitLine;

  while (getline(bf, line)) {
    if (line.at(0) == '#') std::cout << line << std::endl;
    if (line.at(0) == ':') {
      var_nodes++;
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
  std::cout << "\tNumber of variant nodes: " << var_nodes << std::endl;
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
    printBuildHelp();
    throw std::invalid_argument("No input files specified!");
  }

  vargas::Graph::GraphParams p;

  args >> GetOpt::Option('l', "maxlen", p.maxNodeLen)
      >> GetOpt::Option('R', "region", p.region)
      >> GetOpt::OptionPresent('c', "complement", makeComplements)
      >> GetOpt::Option('c', "complement", p.complementSource)
      >> GetOpt::OptionPresent('m', "maxref", p.maxAF)
      >> GetOpt::Option('e', "exref", p.includeRefIndivs);

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
    printExportHelp();
    throw std::invalid_argument("No buildfile defined.");
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

    // Export whole Graph
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
  std::string rawRead = "";
  if (alignOutFile.length() != 0) {
  	// Find the last alignment
    std::ifstream resume(alignOutFile);
    if (resume.good()) {
      std::string lastLine = getLastLine(resume);
      if (lastLine.length() > 0) rawRead = split(lastLine, '#')[0];
    }
    resume.close();
    aOutStream.open(alignOutFile, std::ios_base::app);
    useFile = true;
  }

  g.setParams(p);

  // Read file handler
  std::string readsfile;
  if (!(args >> GetOpt::Option('r', "reads", readsfile))) {
    printAlignHelp();
    throw std::invalid_argument("No reads file defined.");
  }

  // Check if we need to resume this alignment job.
  vargas::ReadFile reads(readsfile);
  if (rawRead.length() != 0) {
  	reads.resume_from(rawRead);
  }


  std::string buildfile;
  if (!(args >> GetOpt::Option('b', "build", buildfile))) {
    printAlignHelp();
    throw std::invalid_argument("No buildfile defined.");
  }
  g.buildGraph(buildfile);

  // Align until there are no more reads
  vargas::Alignment align;
  while (reads.update_read()) {
    g.align(reads.get_read(), align);
    if (useFile) {
      aOutStream << align << std::endl;
    } else {
      std::cout << align << std::endl;
    }
  }

  if (useFile) {
    aOutStream << reads.get_header() << std::endl;
  } else {
    std::cout << reads.get_header() << std::endl;
  }

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
      >> GetOpt::Option('l', "readlen", p.readLen)
      >> GetOpt::Option('a', "ambiguity", p.ambiguity);

  std::string buildfile;
  if (!(args >> GetOpt::Option('b', "buildfile", buildfile))) {
    printBuildHelp();
    throw std::invalid_argument("Buildfile required.");
  }
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
      prof.sub_err = (splitProf[0] == "*") ? -1 : std::stoi(splitProf[0]);
      prof.indel_err = (splitProf[1] == "*") ? -1 : std::stoi(splitProf[1]);
      prof.var_nodes = (splitProf[2] == "*") ? -1 : std::stoi(splitProf[2]);
      prof.var_bases = (splitProf[3] == "*") ? -1 : std::stoi(splitProf[3]);
      sim.addProfile(prof, prefix + std::to_string(i) + ".reads");
    }
    sim.populateProfiles();

  } else {
    for (int i = 0; i < p.maxreads; ++i) {
      std::cout << sim.update_and_get() << std::endl;
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
  cout << "\tbuild     Generate Graph build file from reference FASTA and VCF files." << endl;
  cout << "\tsim       Simulate reads from a Graph." << endl;
  cout << "\talign     Align reads to a Graph." << endl;
  cout << "\tstat      Count nodes and edges of a given Graph." << endl;
  cout << "\texport    Export Graph in DOT format." << endl << endl;
}


void printBuildHelp() {
  using std::cout;
  using std::endl;
  vargas::Graph g;
  vargas::Graph::GraphParams p = g.getParams();

  cout << endl
      << "------------------- vargas build, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-v\t--vcf           <string> VCF file, uncompressed." << endl;
  cout << "-r\t--ref           <string> reference, single record FASTA" << endl;
  cout << "-l\t--maxlen        <int> Maximum node length, default " << p.maxNodeLen << endl;
  cout << "-R\t--region        <<int>:<int>> Ref region, inclusive. Default is the entire Graph." << endl;
  cout << "-m\t--maxref        Generate a linear Graph using maximum allele frequency nodes." << endl;
  cout << "-e\t--exref         Exclude the list of individuals from the reference alleles." << endl;
  cout << "-s\t--set           <<int>,<int>,..,<int>> Generate a buildfile for a list of ingroup percents." << endl;
  cout << "-c\t--complement    <string> Generate a complement of all graphs in -s, or of provided Graph." << endl;

  cout << endl << "--maxref is applied after ingroup filter." << endl;
  cout << "Buildfile is output to [s][In Out].build" << endl << endl;
}


void printExportHelp() {
  using std::cout;
  using std::endl;
  cout << endl
      << "------------------ vargas export, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile    <string> Graph to export to DOT." << endl;
  cout << "-c\t--context      <string> Export the local context Graph of these alignments." << endl;

  cout << endl << "DOT file printed to stdout." << endl << endl;
}


void printSimHelp() {
  using std::cout;
  using std::endl;
  vargas::ReadSim s;
  vargas::SimParams p = s.getParams();


  cout << endl
      << "-------------------- vargas sim, " << __DATE__ << ". rgaddip1@jhu.edu --------------------" << endl;
  cout << "-b\t--buildfile     <string> Graph build file, generate with \'vargas build\'" << endl;
  cout << "-n\t--numreads      <int> Number of reads to simulate, default " << p.maxreads << endl;
  cout << "-m\t--muterr        <float> Read mutation error rate, default " << p.muterr << endl;
  cout << "-i\t--indelerr      <float> Read Indel error rate, default " << p.indelerr << endl;
  cout << "-l\t--readlen       <int> Read length, default " << p.readLen << endl;
  cout << "-e\t--profile       <p1 p2 .. p3> Space delimited read profiles. Produces -n of each" << endl;
  cout << "-p\t--prefix        Prefix to use for read files, default \'sim\'" << endl;
  cout << "-r\t--randwalk      Random walk, read may change individuals at branches" << endl;
  cout << "-a\t--ambiguity     Max number of ambiguous bases to allow in reads, default " << p.ambiguity << endl;
  cout << endl;

  cout << "Outputs to \'[prefix][n].reads\' where [n] is the profile number." << endl;
  cout << "Read Profile format (use \'*\' for any): " << endl;
  cout << "\tsub_err,indel_err,var_nodes,var_bases" << endl;
  cout << "\tExample: Any read with 1 substitution error and 1 variant node." << endl;
  cout << "\t\tvargas sim -b BUILD -e \"1,*,1,*\"" << endl;
  cout << "Read Format:" << endl;
  cout << "\tREAD#READ_END_POSITION,INDIVIDUAL,NUM_SUB_ERR,NUM_INDEL_ERR,NUM_VAR_NODE,NUM_VAR_BASES" << endl << endl;
}


void printAlignHelp() {
  using std::cout;
  using std::endl;

  vargas::Graph::GraphParams p;

  cout << endl
      << "------------------- vargas align, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     <string> Graph build file" << endl;
  cout << "-m\t--match         <int> Match score, default " << int(p.match) << endl;
  cout << "-n\t--mismatch      <int> Mismatch score, default " << int(p.mismatch) << endl;
  cout << "-o\t--gap_open      <int> Gap opening penalty, default " << int(p.gap_open) << endl;
  cout << "-e\t--gap_extend    <int> Gap extend penalty, default " << int(p.gap_extension) << endl;
  cout << "-r\t--reads         <string> Reads to align. Use stdin if not defined." << endl;
  cout << "-f\t--outfile       <string> Alignment output file. If not defined, use stdout." << endl;

  cout << "Lines beginning with \'#\' are ignored." << endl;
  cout << "With -f, alignments will begin after the last read in -f in -r, i.e. \'resume\'." << endl;
  cout << "Output format:" << endl;
  cout << "\tREAD,OPTIMAL_SCORE,OPTIMAL_ALIGNMENT_END,NUM_OPTIMAL_ALIGNMENTS,SUBOPTIMAL_SCORE,";
  cout << "SUBOPTIMAL_ALIGNMENT_END,NUM_SUBOPTIMAL_ALIGNMENTS,ALIGNMENT_MATCH" << endl << endl;
  cout << "ALIGNMENT_MATCH:\n\t0- optimal match, 1- suboptimal match, 2- no match" << endl << endl;

}


void printStatHelp() {
  using std::cout;
  using std::endl;

  cout << endl
      << "------------------- vargas stat, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     <string> Graph build file." << endl << endl;

}

*/