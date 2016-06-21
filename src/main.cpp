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
        else if (!strcmp(argv[1], "define")) {
            exit(define_main(argc, argv));
        }
        else if (!strcmp(argv[1], "sim")) {
            exit(sim_main(argc, argv));
        }
        else if (!strcmp(argv[1], "align")) {
            exit(align_main(argc, argv));
        }
    }

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
    std::string outfile, readsfile, gdeffile;
    bool R = false, X = false, O = false, I = false;

    args >> GetOpt::Option('m', "match", match)
        >> GetOpt::Option('n', "mismatch", mismatch)
        >> GetOpt::Option('o', "gap_open", gopen)
        >> GetOpt::Option('e', "gap_extend", gext)
        >> GetOpt::Option('t', "outfile", outfile)
        >> GetOpt::Option('r', "reads", readsfile)
        >> GetOpt::Option('g', "gdef", gdeffile)
        >> GetOpt::Option('j', "threads", threads)
        >> GetOpt::Option('l', "rlen", read_len)
        >> GetOpt::OptionPresent('R', R)
        >> GetOpt::OptionPresent('X', X)
        >> GetOpt::OptionPresent('O', O)
        >> GetOpt::OptionPresent('I', I);

    if (threads == 0) threads = std::thread::hardware_concurrency();

    std::ofstream out(outfile);
    if (!out.good()) throw std::invalid_argument("Error opening output file " + outfile);


    std::unordered_map<std::string, vargas::Graph::Population> pop_defs;
    auto gb = load_gdef(gdeffile, pop_defs);

    vargas::Graph base_graph = gb.build();

    out << base_graph.desc() << std::endl;

    std::vector<std::shared_ptr<vargas::Aligner<>>> aligners;
    for (uint8_t i = 0; i < threads; ++i)
        aligners.push_back(std::make_shared<vargas::Aligner<>>(base_graph.max_node_len(),
                                                               read_len,
                                                               match,
                                                               mismatch,
                                                               gopen,
                                                               gext));

    // Load reads
    std::unordered_map<std::string, std::vector<vargas::Read>> read_origins;
    auto readfile_split = split(readsfile, ',');
    for (auto &rfile : readfile_split) {
        vargas::ReadFile reads(rfile);
        while (reads.update_read()) {
            vargas::Read r = reads.get_read();
            read_origins[r.desc].push_back(r);
        }
    }

    if (R) {
        vargas::Graph subgraph = vargas::Graph(base_graph, vargas::Graph::REF);
        for (auto &reads : read_origins) align_to_graph("R", subgraph, reads.second, aligners, out, threads);
    }
    if (X) {
        vargas::Graph subgraph = vargas::Graph(base_graph, vargas::Graph::MAXAF);
        for (auto &reads : read_origins) align_to_graph("X", subgraph, reads.second, aligners, out, threads);
    }
    if (I) {
        for (auto &reads : read_origins) {
            vargas::Graph subgraph = vargas::Graph(base_graph, pop_defs.at(reads.first));
            align_to_graph("I", subgraph, reads.second, aligners, out, threads);
        }
    }
    if (O) {
        for (auto &reads : read_origins) {
            vargas::Graph subgraph = vargas::Graph(base_graph, ~(pop_defs.at(reads.first)));
            align_to_graph("O", subgraph, reads.second, aligners, out, threads);
        }
    }
    return 0;
}

void align_to_graph(std::string label,
                    const vargas::Graph &subgraph,
                    const std::vector<vargas::Read> &reads,
                    const std::vector<std::shared_ptr<vargas::Aligner<>>> &aligners,
                    std::ostream &out,
                    int threads) {
    size_t batches = reads.size() / SIMDPP_FAST_INT8_SIZE;

    std::vector<std::thread> jobs;
    std::vector<std::vector<vargas::Alignment>> aligns(threads);

    for (size_t i = 0; i <= batches; ++i) {
        auto end_iter = ((i + 1) * SIMDPP_FAST_INT8_SIZE) > reads.size() ?
                        reads.end() : reads.begin() + ((i + 1) * SIMDPP_FAST_INT8_SIZE);

        jobs.emplace(jobs.end(),
                     &vargas::Aligner<>::align_into,
                     aligners[jobs.size()],
                     std::vector<vargas::Read>(reads.begin() + (i * SIMDPP_FAST_INT8_SIZE), end_iter),
                     subgraph.begin(),
                     subgraph.end(),
                     std::ref(aligns[jobs.size()]));

        if (jobs.size() == threads || i == batches) {
            std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });
            jobs.clear();

            for (auto &aset : aligns) {
                for (auto &a : aset) {
                    out << label << ',' << a << std::endl;
                }
            }
        }
    }
}

int sim_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        sim_help();
        return 0;
    }

    // Load parameters
    unsigned int threads = 1, read_len = 50, num_reads = 1000;
    std::string mut = "0", indel = "0";
    std::string outfile, gdeffile, sources = "";
    bool use_rate = false, exclude_comment = false;

    args >> GetOpt::Option('t', "outfile", outfile)
        >> GetOpt::Option('g', "gdef", gdeffile)
        >> GetOpt::Option('s', "source", sources)
        >> GetOpt::Option('j', "threads", threads)
        >> GetOpt::Option('l', "rlen", read_len)
        >> GetOpt::OptionPresent('a', "rate", use_rate)
        >> GetOpt::Option('n', "numreads", num_reads)
        >> GetOpt::OptionPresent('e', "exclude", exclude_comment)
        >> GetOpt::Option('m', "mut", mut)
        >> GetOpt::Option('d', "indel", indel);

    if (threads == 0) threads = std::thread::hardware_concurrency();

    std::ofstream out(outfile);
    if (!out.good()) throw std::invalid_argument("Error opening output file " + outfile);

    std::unordered_map<std::string, vargas::Graph::Population> pop_defs, pops;
    auto gb = load_gdef(gdeffile, pop_defs);

    if (sources.length() == 0) pops = pop_defs;
    else {
        std::vector<std::string> source_split = split(sources, ',');
        for (auto &i : source_split) {
            if (pop_defs.find(i) == pop_defs.end())
                throw std::invalid_argument("\"" + i + "\" not found in GDEF file.");
            else pops[i] = pop_defs[i];
        }
    }



    vargas::Graph base_graph = gb.build();

    if (!exclude_comment) out << base_graph.desc() << std::endl;
    if (!exclude_comment) for (auto &p : pops) out << "#" << p.first << ":" << p.second.to_string() << std::endl;

    auto mut_split = split(mut, ',');
    auto indel_split = split(indel, ',');

    // For each subgraph type
    for (auto &pop : pops) {
        // Create subgraph
        vargas::Graph subgraph;
        subgraph = vargas::Graph(base_graph, pop.second);
        std::vector<std::thread> jobs;

        std::vector<std::shared_ptr<vargas::ReadSim>> sims;
        vargas::ReadSim rs_template(subgraph);
        for (int i = 0; i < threads; ++i) sims.push_back(std::make_shared<vargas::ReadSim>(rs_template));

        for (auto &ind : indel_split) {
            for (auto &m : mut_split) {
                vargas::ReadProfile prof;
                prof.len = read_len;
                prof.mut = std::stof(m);
                prof.indel = std::stof(ind);
                prof.rand = use_rate;

                if (!exclude_comment) out << "#prof:" << prof << std::endl;

                for (unsigned int i = 0; i < threads; ++i) {
                    sims[i]->set_prof(prof);
                    jobs.emplace(jobs.end(),
                                 &vargas::ReadSim::get_batch,
                                 sims[i],
                                 num_reads / threads);
                }

                std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });

                for (auto &s : sims) {
                    for (auto r : s->batch()) {
                        r.desc = pop.first;
                        out << r << std::endl;

                    }
                }
            }
        }
    }

    return 0;
}

int define_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        define_help();
        return 0;
    }

    std::string fastafile = "", varfile = "", region = "", ingroups = "", outfile = "";
    int base_ingroup = 100, num = 1, node_len = 1000000;
    bool outgroup = true;

    args >> GetOpt::Option('f', "fasta", fastafile)
        >> GetOpt::Option('v', "var", varfile)
        >> GetOpt::Option('g', "region", region)
        >> GetOpt::Option('b', "base", base_ingroup)
        >> GetOpt::Option('i', "ingroup", ingroups)
        >> GetOpt::Option('n', "num", num)
        >> GetOpt::Option('l', "nodelen", node_len)
        >> GetOpt::Option('t', "out", outfile);

    std::ofstream out(outfile);
    if (!out.good()) throw std::invalid_argument("Unable to open output file \"" + outfile + "\"");

    // Check if files are valid, create fai if necessary
    vargas::VarFile vf(varfile);
    vargas::FASTAFile fa(fastafile);

    out << "@gdef" << std::endl
        << "ref=" << fastafile << ';'
        << "var=" << varfile << ';'
        << "base=" << base_ingroup << ';'
        << "reg=" << region << ';'
        << "nlen=" << node_len << ';'
        << "outgroup=" << outgroup << std::endl;

    std::vector<std::string> ingroup_split = split(ingroups, ',');

    for (auto &ingrp_str : ingroup_split) {
        int ingrp = std::stoi(ingrp_str);
        for (int n = 0; n < num; ++n) {
            dyn_bitset<64> filter(vf.num_samples() * 2);
            for (int i = 0; i < filter.size(); ++i) {
                if (rand() % 100 < ingrp) filter.set(i);
            }
            out << std::to_string(ingrp) + "i" + std::to_string(n) << ":" << filter.to_string() << std::endl;
        }
    }

}

vargas::GraphBuilder load_gdef(std::string file_name,
                               std::unordered_map<std::string, vargas::Graph::Population> &pset) {
    std::ifstream in(file_name);
    if (!in.good()) throw std::invalid_argument("Unable to open file \"" + file_name + "\"");

    std::string line;
    std::getline(in, line);
    if (line != "@gdef") throw std::invalid_argument("Expected a GDEF file.");

    std::getline(in, line);
    auto line_split = split(line, ';');

    std::string ref, var, reg, base, nlen;
    bool outgroup;

    for (auto &tpair : line_split) {
        auto tpair_split = split(tpair, '=');
        auto tag = tpair_split[0];
        auto val = tpair_split[1];

        if (tag == "ref") {
            ref = val;
        }
        else if (tag == "var") {
            var = val;
        }
        else if (tag == "base") {
            base = val;
        }
        else if (tag == "reg") {
            reg = val;
        }
        else if (tag == "nlen") {
            nlen = val;
        }
        else if (tag == "outgroup") {
            outgroup = std::stoi(val);
        }
    }

    vargas::GraphBuilder gb(ref, var);
    gb.node_len(std::stoi(nlen));
    gb.region(reg);
    gb.ingroup(std::stoi(base));

    pset.clear();
    while (std::getline(in, line)) {
        split(line, ':', line_split);
        vargas::Graph::Population pop(line_split[1].size());
        for (int i = 0; i < line_split[1].length(); ++i) {
            if (line_split[1][i] == '1') pop.set(i);
        }
        pset[line_split[0]] = pop;
        if (outgroup) {
            line_split[0][line_split[0].find('i')] = 'o';
            pset[line_split[0]] = ~pop;
        }
    }

    return gb;
}

void define_help() {
    using std::cout;
    using std::endl;

    cout << endl
        << "-------------------- vargas define, " << __DATE__ << ". rgaddip1@jhu.edu --------------------" << endl;
    cout << "-f\t--fasta         <string> Reference filename." << endl;
    cout << "-v\t--var           <string> VCF/BCF filename." << endl;
    cout << "-t\t--out           <string> Output filename." << endl;
    cout << "-g\t--region        <string> Region of graph, format CHR:MIN-MAX." << endl;
    cout << "-l\t--nodelen       <int> Max node length, default 1,000,000" << endl;
    cout << "-b\t--base          <int> Base graph ingroup. -i graphs derived from this. Default 100%." << endl;
    cout << "-i\t--ingroup       <int, int...> Percent ingroup subgraphs." << endl;
    cout << "-n\t--num           <int> Number of unique graphs to create for all -i, -x. Default 1." << endl;
}

void main_help() {
    using std::cout;
    using std::endl;
    cout << endl
        << "---------------------- vargas, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------" << endl;
    cout << "Operating modes \'vargas MODE\':" << endl;
    cout << "\tdefine      Define a set of graphs for use with sim/align." << endl;
    cout << "\tsim         Simulate reads from a graph." << endl;
    cout << "\talign       Align reads to a graph." << endl;
    cout << "\ttest        Run doctests." << endl;
    cout << "\tprofile     Run profiles." << endl << endl;
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
    cout << "-g\t--gdef          <string> Graph definition file." << endl;
    cout << "-r\t--reads         <string, string...> Read files to align" << endl;
    cout << "-R\t                Align to reference graph." << endl;
    cout << "-X\t                Align to maximum allele frequency graph." << endl;
    cout << "-I\t                Align to ingroup graph." << endl;
    cout << "-O\t                Align to outgroup graph." << endl;
    cout << "-m\t--match         <int> Match score, default 2" << endl;
    cout << "-n\t--mismatch      <int> Mismatch penalty, default 2" << endl;
    cout << "-o\t--gap_open      <int> Gap opening penalty, default 3" << endl;
    cout << "-e\t--gap_extend    <int> Gap extend penalty, default 1" << endl;
    cout << "-l\t--rlen          <int> Max read length. Default 50." << endl;
    cout << "-t\t--outfile       <string> Alignment output file." << endl;
    cout << "-j\t--threads       <int> Number of threads. 0 for maximum hardware concurrency." << endl << endl;

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
    cout << "-g\t--gdef          <string> Graph definition file." << endl;
    cout << "-s\t--source       <string, string...> Simulate from specified subgraphs, default all." << endl;
    cout << "-n\t--numreads      <int> Number of reads to simulate from each subgraph." << endl;
    cout << "-m\t--muterr        <float, float...> Read mutation error. Default 0." << endl;
    cout << "-i\t--indelerr      <float, float...> Read indel error. Default 0." << endl;
    cout << "-v\t--vnodes        <int, int...> Number of variant nodes, default 0." << endl;
    cout << "-l\t--rlen          <int> Read length, default " << endl;
    cout << "-t\t--outfile       <string> Alignment output file." << endl;
    cout << "-a\t--rate          Interpret -m, -i as rates, instead of exact number of errors." << endl;
    cout << "-e\t--exclude       Exclude comments, include only FASTA lines" << endl;
    cout << "-j\t--threads       <int> Number of threads. 0 for maximum hardware concurrency." << endl << endl;

    cout << "Comments preceded by \'#\'." << endl;
    cout << "-n reads are produced for each -s, -m, -i, -v combination." << endl << endl;
}