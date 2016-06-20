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

#include <iostream>
#include <thread>
#include <algorithm>
#include "doctest/doctest.h"
#include "main.h"
#include "getopt_pp.h"


int main(const int argc, const char *argv[]) {


    //   try {
        if (argc > 1) {
            if (!strcmp(argv[1], "test")) {
                doctest::Context doc(argc, argv);
                doc.setOption("no-breaks", true);
                doc.setOption("abort-after", 5);
                doc.setOption("sort", "name");
                exit(doc.run());
            }
            else if (!strcmp(argv[1], "profile")) {
                exit(profile(argc, argv));
            }
            else if (!strcmp(argv[1], "build")) {
//        exit(build_main(argc, argv));
            }
            else if (!strcmp(argv[1], "sim")) {
                exit(sim_main(argc, argv));
            }
            else if (!strcmp(argv[1], "align")) {
                exit(align_main(argc, argv));
            }
            else if (!strcmp(argv[1], "export")) {
//       exit(export_main(argc, argv));
            }
            else if (!strcmp(argv[1], "stat")) {
//       exit(stat_main(argc, argv));
            }
        }
//    } catch (std::exception &e) {
//        std::cerr << "EXCEPTION: " << e.what() << std::endl;
//        exit(1);
//    }

    GetOpt::GetOpt_pp args(argc, argv);
    if (args >> GetOpt::OptionPresent('h', "help")) {
        main_help();
        exit(0);
    }

    std::cerr << "Define a valid mode of operation." << std::endl;
    main_help();
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

void printStatHelp() {
  using std::cout;
  using std::endl;

  cout << endl
      << "------------------- vargas stat, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     <string> Graph build file." << endl << endl;

}

*/

int profile(const int argc, const char *argv[]) {
    std::string bcf = "chr22.bcf", fasta = "hs37d5_22.fa";
    std::string region = "22:25,000,000-25,500,000";
    std::string read;
    int ingroup = 100;

    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        profile_help();
        return 0;
    }

    args >> GetOpt::Option('f', "fasta", fasta)
        >> GetOpt::Option('v', "var", bcf)
        >> GetOpt::Option('g', "region", region)
        >> GetOpt::Option('i', "ingroup", ingroup)
        >> GetOpt::Option('s', "string", read);


    if (!file_exists(fasta) || !file_exists(bcf)) {
        throw std::invalid_argument("File does not exist.");
    }

    vargas::GraphBuilder gb(fasta, bcf);
    gb.region(region);
    if (!gb.good()) throw std::invalid_argument("Error opening files\n" + fasta + "\nand/or:\n" + bcf);

    srand(time(NULL));
    std::clock_t start = std::clock();

    std::cout << "Initial Build:\n\t";
    vargas::Graph g;
    gb.build(g);
    std::vector<bool> filter;
    for (size_t i = 0; i < g.pop_size(); ++i) filter.push_back(rand() % 100 > 95);
    std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s, " << "Nodes: " << g.node_map()->size()
        << std::endl;

    size_t num = 0;

    {
        std::cout << "Insertion order traversal:\n\t";
        start = std::clock();
        for (auto i = g.begin(); i != g.end(); ++i) {
            ++num;
        }
        std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        num = 0;
        std::cout << "Filtering traversal, 100% in:\n\t";
        start = std::clock();

        for (auto i = g.begin(g.subset(100)); i != g.end(); ++i) {
            ++num;
        }
        std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        num = 0;
        std::cout << "Filtering traversal, 5% in:\n\t";
        vargas::Graph::Population filt(filter);
        start = std::clock();

        for (auto i = g.begin(filt); i != g.end(); ++i) {
            ++num;
        }
        std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s, " << "Nodes: " << num << std::endl;
    }

    {
        num = 0;
        std::cout << "REF traversal:\n\t";
        start = std::clock();
        for (auto i = g.begin(vargas::Graph::REF); i != g.end(); ++i) {
            ++num;
        }
        std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        num = 0;
        std::cout << "MAXAF traversal:\n\t";
        start = std::clock();
        for (auto i = g.begin(vargas::Graph::MAXAF); i != g.end(); ++i) { ;
        }
        std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        std::cout << "Filter constructor:\n\t";
        auto pop_filt = g.subset(ingroup);
        start = std::clock();
        vargas::Graph g2(g, pop_filt);
        std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        std::cout << "REF constructor:\n\t";
        start = std::clock();
        vargas::Graph g2(g, vargas::Graph::REF);
        std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        std::cout << "MAXAF constructor:\n\t";
        start = std::clock();
        vargas::Graph g2(g, vargas::Graph::MAXAF);
        std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }


    node_fill_profile();


    {
        std::vector<vargas::Read> reads;

        for (int i = 0; i < SIMDPP_FAST_INT8_SIZE; ++i) {
            std::stringstream rd;
            for (int r = 0; r < 50; ++r) rd << rand_base();
            reads.push_back(rd.str());
        }

        std::vector<std::string> split_str;
        if (read.length() > 0) {
            split(read, ',', split_str);
            if (split_str.size() > SIMDPP_FAST_INT8_SIZE) split_str.resize(SIMDPP_FAST_INT8_SIZE);
            for (size_t i = 0; i < split_str.size(); ++i) reads[i] = vargas::Read(split_str[i]);
        }

        vargas::ReadBatch<> rb(reads, 50);
        vargas::Aligner<> a(g.max_node_len(), 50);

        std::cout << SIMDPP_FAST_INT8_SIZE << " read alignment:\n\t";
        std::cout << "Filtering iterator:\n\t";
        auto pop_filt = g.subset(ingroup);

        {
            start = std::clock();
            std::vector<vargas::Alignment> aligns = a.align(rb, g.begin(pop_filt), g.end());
            std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
            for (size_t i = 0; i < split_str.size(); ++i) std::cout << aligns[i] << std::endl;
        }

        {
            start = std::clock();
            vargas::Graph g2(g, g.subset(ingroup));
            std::cout << "\tDerived Graph (" << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s)\n\t";
            start = std::clock();
            std::vector<vargas::Alignment> aligns = a.align(rb, g2);
            std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
            for (size_t i = 0; i < split_str.size(); ++i) std::cout << aligns[i] << std::endl;
        }

    }
    return 0;
}

int align_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        align_help();
        return 0;
    }

    // Load parameters
    uint8_t match = 2, mismatch = 2, gopen = 3, gext = 1;
    unsigned int threads = 1, read_len = 50;
    std::string outfile, readsfile, reffile, varfile, region, ingroups = "100";
    bool outgroups = false;

    args >> GetOpt::Option('m', "match", match)
        >> GetOpt::Option('n', "mismatch", mismatch)
        >> GetOpt::Option('o', "gap_open", gopen)
        >> GetOpt::Option('e', "gap_extend", gext)
        >> GetOpt::Option('t', "outfile", outfile)
        >> GetOpt::Option('r', "reads", readsfile)
        >> GetOpt::Option('f', "fasta", reffile)
        >> GetOpt::Option('v', "var", varfile)
        >> GetOpt::Option('g', "region", region)
        >> GetOpt::Option('i', "ingroup", ingroups)
        >> GetOpt::Option('j', "threads", threads)
        >> GetOpt::Option('l', "rlen", read_len)
        >> GetOpt::OptionPresent('x', "outgroup", outgroups);

    if (threads == 0) threads = std::thread::hardware_concurrency();

    std::ofstream out(outfile);
    if (!out.good()) throw std::invalid_argument("Error opening output file " + outfile);

    // Build base graph, all graphs are derived from this
    vargas::GraphBuilder gb(reffile, varfile);
    gb.region(region);
    vargas::Graph base_graph;
    gb.build(base_graph);

    out << base_graph.desc() << std::endl;

    std::vector<vargas::Read> read_batch;
    std::vector<std::shared_ptr<vargas::Aligner<>>> aligners;
    for (uint8_t i = 0; i < threads; ++i)
        aligners.push_back(std::make_shared<vargas::Aligner<>>(base_graph.max_node_len(),
                                                               read_len,
                                                               match,
                                                               mismatch,
                                                               gopen,
                                                               gext));
    auto readfile_split = split(readsfile, ',');
    auto ingroups_split = split(ingroups, ',');

    // Create populations
    std::unordered_map<std::string, vargas::Graph::Population> pops;
    for (auto igrp : ingroups_split) {
        int pct = std::stoi(igrp);
        if (pct < 0) continue;
        auto pop = base_graph.subset(pct);
        pops.insert(std::make_pair(igrp + 'i', pop));
        if (pct >= 0 && outgroups) {
            pops.insert(std::make_pair(igrp + 'o', ~pop));
        }
    }
    for (auto &p : pops) out << "#" << p.first << ":" << p.second.to_string() << std::endl;

    // For each read file
    for (auto rfile : readfile_split) {
        // For each subgraph type
        for (auto &pop : pops) {
            vargas::ReadFile reads(rfile);

            // Create subgraph
            vargas::Graph subgraph;
            subgraph = vargas::Graph(base_graph, pop.second);

            while (true) {
                std::vector<std::thread> jobs;
                std::vector<std::vector<vargas::Alignment>> aligns(threads);
                for (unsigned int i = 0; i < threads; ++i) {
                    read_batch = reads.get_batch(SIMDPP_FAST_INT8_SIZE);
                    if (read_batch.size() == 0) break;
                    jobs.emplace(jobs.end(),
                                 &vargas::Aligner<>::align_into,
                                 aligners[i],
                                 read_batch,
                                 subgraph.begin(),
                                 subgraph.end(),
                                 std::ref(aligns[i]));
                }
                if (jobs.size() == 0) break; // no more reads

                std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });

                for (auto &aset : aligns) {
                    for (auto &a : aset) {
                        out << pop.first << ',' << a << std::endl;
                    }
                }
            }
        }
    }
    return 0;
}

int sim_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        sim_help();
        return 0;
    }

    // Load parameters
    unsigned int threads = 1, read_len = 50, num_sub = 0, num_indel = 0, num_reads = 1000;
    float sub_err = 0.02, indel_err = 0;
    std::string outfile, reffile, varfile, region, ingroups = "100";
    bool outgroups = false, use_rate = false;

    args >> GetOpt::Option('t', "outfile", outfile)
        >> GetOpt::Option('f', "fasta", reffile)
        >> GetOpt::Option('v', "var", varfile)
        >> GetOpt::Option('g', "region", region)
        >> GetOpt::Option('i', "ingroup", ingroups)
        >> GetOpt::Option('j', "threads", threads)
        >> GetOpt::Option('l', "rlen", read_len)
        >> GetOpt::OptionPresent('x', "outgroup", outgroups)
        >> GetOpt::OptionPresent('a', "rate", use_rate)
        >> GetOpt::Option('n', "numreads", num_reads);

    if (use_rate) {
        args >> GetOpt::Option('m', "mut", sub_err)
            >> GetOpt::Option('d', "indel", indel_err);
    } else {
        args >> GetOpt::Option('m', "mut", num_sub)
            >> GetOpt::Option('d', "indel", num_indel);
    }

    if (threads == 0) threads = std::thread::hardware_concurrency();

    std::ofstream out(outfile);
    if (!out.good()) throw std::invalid_argument("Error opening output file " + outfile);

    // Build base graph, all graphs are derived from this
    vargas::GraphBuilder gb(reffile, varfile);
    gb.region(region);
    vargas::Graph base_graph = gb.build();

    out << base_graph.desc() << std::endl;

    std::vector<vargas::Read> read_batch;
    std::vector<std::shared_ptr<vargas::ReadSim>> sims;
    for (uint8_t i = 0; i < threads; ++i)
        sims.push_back(std::make_shared<vargas::ReadSim>(base_graph));
    auto ingroups_split = split(ingroups, ',');

    // Create populations
    std::unordered_map<std::string, vargas::Graph::Population> pops;
    for (auto igrp : ingroups_split) {
        int pct = std::stoi(igrp);
        if (pct < 0) continue;
        auto pop = base_graph.subset(pct);
        pops.insert(std::make_pair(igrp + 'i', pop));
        if (pct >= 0 && outgroups) {
            pops.insert(std::make_pair(igrp + 'o', ~pop));
        }
    }
    for (auto &p : pops) out << "#" << p.first << ":" << p.second.to_string() << std::endl;

    vargas::ReadProfile prof;
    prof.len = read_len;
    prof.num_mut = num_sub;
    prof.num_indel = num_indel;
    prof.mut_rate = sub_err;
    prof.indel_rate = indel_err;
    prof.rand = use_rate;

    out << "#prof:" << prof << std::endl;

    // For each subgraph type
    for (auto &pop : pops) {
        // Create subgraph
        vargas::Graph subgraph;
        subgraph = vargas::Graph(base_graph, pop.second);
        std::vector<std::thread> jobs;
        for (unsigned int i = 0; i < threads; ++i) {
            sims[i]->set_prof(prof);
            jobs.emplace(jobs.end(),
                         &vargas::ReadSim::get_batch,
                         sims[i],
                         num_reads / threads);
        }

        std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });

        for (auto &s : sims) {
            for (auto &r : s->batch()) {
                out << '>' << pop.first
                    << " end=" << r.end_pos
                    << " mut=" << r.sub_err
                    << " indel=" << r.indel_err
                    << " vnode=" << r.var_nodes
                    << " vbase=" << r.var_bases
                    << " desc=" << r.desc
                    << std::endl
                    << r.read
                    << std::endl;
            }
        }
    }

    return 0;
}

void main_help() {
    using std::cout;
    using std::endl;
    cout << endl
        << "---------------------- vargas, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------" << endl;
    cout << "Operating modes \'vargas MODE\':" << endl;
    cout << "\ttest        Run doctests." << endl;
    cout << "\tprofile     Run profiles." << endl;
    cout << "\tsim         Simulate reads from a graph." << endl;
    cout << "\talign       Align reads to a graph." << endl;
    cout << "\texport      Export Graph in DOT format." << endl << endl;
}

void profile_help() {
    using std::cout;
    using std::endl;
    cout << endl
        << "---------------------- vargas profile, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------" << endl;
    cout << "-f\t--fasta         <string> Reference filename." << endl;
    cout << "-v\t--var           <string> VCF/BCF filename." << endl;
    cout << "-g\t--region        <string> Region of graph, format CHR:MIN-MAX." << endl;
    cout << "-i\t--ingroup       <int> Percent of genotypes to include in alignment" << endl;
    cout << "-s\t--string        <string,string..> Include reads in alignment. Rest will be random." << endl << endl;
}

void align_help() {
    using std::cout;
    using std::endl;

    cout << endl
        << "------------------- vargas align, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
    cout << "-f\t--fasta         <string> Reference filename." << endl;
    cout << "-v\t--var           <string> VCF/BCF filename." << endl;
    cout << "-g\t--region        <string> Region of graph, format CHR:MIN-MAX." << endl;
    cout << "-i\t--ingroup       <int, int...> Align to percent ingroup subgraphs. -1 for max AF. Default 100." << endl;
    cout << "-x\t--outgroup      Align to outgroups for all -i" << endl;
    cout << "-m\t--match         <int> Match score, default 2" << endl;
    cout << "-n\t--mismatch      <int> Mismatch penalty, default 2" << endl;
    cout << "-o\t--gap_open      <int> Gap opening penalty, default 3" << endl;
    cout << "-e\t--gap_extend    <int> Gap extend penalty, default 1" << endl;
    cout << "-r\t--reads         <string, string...> Read file to align" << endl;
    cout << "-l\t--rlen          <int> Max read length. Default 50." << endl;
    cout << "-t\t--outfile       <string> Alignment output file." << endl;
    cout << "-j\t--threads       Number of threads. 0 for maximum hardware concurrency." << endl << endl;

    cout << "Lines beginning with \'#\' are ignored." << endl;
    cout << "Output format:" << endl;
    cout << "\tINGROUP,READ,OPTIMAL_SCORE,OPTIMAL_ALIGNMENT_END,NUM_OPTIMAL_ALIGNMENTS,SUBOPTIMAL_SCORE,";
    cout << "SUBOPTIMAL_ALIGNMENT_END,NUM_SUBOPTIMAL_ALIGNMENTS,ALIGNMENT_MATCH" << endl << endl;
    cout << "ALIGNMENT_MATCH:\n\t0- optimal match, 1- suboptimal match, 2- no match" << endl << endl;

}

void sim_help() {
    using std::cout;
    using std::endl;

    cout << endl
        << "-------------------- vargas sim, " << __DATE__ << ". rgaddip1@jhu.edu --------------------" << endl;
    cout << "-n\t--numreads      <int> Number of reads to simulate, default " << endl;
    cout << "-m\t--muterr        <float> Read mutation error rate, default " << endl;
    cout << "-i\t--indelerr      <float> Read Indel error rate, default " << endl;
    cout << "-l\t--readlen       <int> Read length, default " << endl;
    cout << "-e\t--profile       <p1 p2 .. p3> Space delimited read profiles. Produces -n of each" << endl;
    cout << "-p\t--prefix        Prefix to use for read files, default \'sim\'" << endl;
    cout << "-r\t--randwalk      Random walk, read may change individuals at branches" << endl;
    cout << "-a\t--ambiguity     Max number of ambiguous bases to allow in reads, default " << endl;
    cout << endl;

    cout << "Outputs to \'[prefix][n].reads\' where [n] is the profile number." << endl;
    cout << "Read Profile format (use \'*\' for any): " << endl;
    cout << "\tsub_err,indel_err,var_nodes,var_bases" << endl;
    cout << "\tExample: Any read with 1 substitution error and 1 variant node." << endl;
    cout << "\t\tvargas sim -b BUILD -e \"1,*,1,*\"" << endl;
    cout << "Read Format:" << endl;
    cout << "\tREAD#READ_END_POSITION,INDIVIDUAL,NUM_SUB_ERR,NUM_INDEL_ERR,NUM_VAR_NODE,NUM_VAR_BASES" << endl << endl;
}