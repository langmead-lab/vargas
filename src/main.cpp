/**
 * Ravi Gaddipati
 * November 25, 2015
 * rgaddip1@jhu.edu
 *
 * @brief
 * Interface for simulating and aligning reads from/to a DAG.
 *
 * @file
 */


#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <iostream>
#include <thread>
#include <algorithm>
#include "main.h"
#include "getopt_pp.h"
#include "sam.h"


int main(const int argc, const char *argv[]) {

    srand(time(NULL)); // Rand used in profiles and sim.h

    if (argc > 1) {
        if (!strcmp(argv[1], "test")) {
            doctest::Context doc(argc, argv);
            doc.setOption("no-breaks", true);
            doc.setOption("abort-after", 5);
            doc.setOption("sort", "name");
            return doc.run();
        }
        else if (!strcmp(argv[1], "profile")) {
            return profile(argc, argv);
        }
        else if (!strcmp(argv[1], "define")) {
            return define_main(argc, argv);
        }
        else if (!strcmp(argv[1], "sim")) {
            return sim_main(argc, argv);
        }
        else if (!strcmp(argv[1], "align")) {
            return align_main(argc, argv);
        }
    }

    GetOpt::GetOpt_pp args(argc, argv);
    if (args >> GetOpt::OptionPresent('h', "help")) {
        main_help();
        return 0;
    }

    std::cerr << "Define a valid mode of operation." << std::endl;
    main_help();
    return 1;

}

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
    std::string out_file, read_file, gdf_file;
    bool R = false, // Align to reference
        X = false, // Align to max-AF graph
        O = false, // Align to outgroup graph
        I = false; // Align to ingroup graph

    args >> GetOpt::Option('m', "match", match)
        >> GetOpt::Option('n', "mismatch", mismatch)
        >> GetOpt::Option('o', "gap_open", gopen)
        >> GetOpt::Option('e', "gap_extend", gext)
        >> GetOpt::Option('t', "outfile", out_file)
        >> GetOpt::Option('r', "reads", read_file)
        >> GetOpt::Option('g', "gdef", gdf_file)
        >> GetOpt::Option('j', "threads", threads)
        >> GetOpt::Option('l', "rlen", read_len)
        >> GetOpt::OptionPresent('R', R)
        >> GetOpt::OptionPresent('X', X)
        >> GetOpt::OptionPresent('O', O)
        >> GetOpt::OptionPresent('I', I);

    if (!R && !X && !O && !I) throw std::invalid_argument("No alignment mode (R,X,I,O) selected.");
    if (threads == 0) threads = std::thread::hardware_concurrency();

    std::ofstream out(out_file);
    if (!out.good()) throw std::invalid_argument("Error opening output file " + out_file);

    vargas::Gdef gdf(gdf_file);
    vargas::GraphBuilder gb(gdf.fasta(), gdf.var());
    gb.region(gdf.region());
    std::cout << "Loading graph..." << std::endl;
    vargas::Graph base_graph = gb.build();

    std::vector<std::shared_ptr<vargas::Aligner<>>> aligners;
    for (uint8_t i = 0; i < threads; ++i)
        aligners.push_back(std::make_shared<vargas::Aligner<>>(base_graph.max_node_len(),
                                                               read_len,
                                                               match,
                                                               mismatch,
                                                               gopen,
                                                               gext));

    // Load reads. Maps graph origin label to a vector of reads
    std::cout << "Loading reads..." << std::endl;
    std::map<vargas::GID, std::vector<vargas::Read>> read_origins;
    auto readfile_split = split(read_file, ',');
    for (auto &rfile : readfile_split) {
        vargas::ReadFile reads(rfile);
        while (reads.update_read()) {
            vargas::Read r = reads.get_read();
            read_origins[r.src].push_back(r);
        }
    }

    auto &pops = gdf.populations();

    if (R) {
        std::cout << "Aligning to reference..." << std::endl;
        vargas::Graph subgraph = vargas::Graph(base_graph, vargas::Graph::REF);
        for (auto &reads : read_origins) {
            std::cout << "\tGID: " << reads.first << std::endl;
            align_to_graph("r", subgraph, reads.second, aligners, out, threads);
        }
        std::cout << "Done." << std::endl;
    }
    if (X) {
        std::cout << "Aligning to MAX AF..." << std::endl;
        vargas::Graph subgraph = vargas::Graph(base_graph, vargas::Graph::MAXAF);
        for (auto &reads : read_origins) {
            std::cout << "\tGID: " << reads.first << std::endl;
            align_to_graph("x", subgraph, reads.second, aligners, out, threads);
        }
        std::cout << "Done." << std::endl;
    }
    if (I) {
        std::cout << "Aligning to Ingroup..." << std::endl;
        for (auto &reads : read_origins) {
            /**
             * Align to the same graph that the read came from, but enforce
             * that it be the ingroup version.
             */
            std::cout << "\tGID: " << reads.first << std::endl;
            vargas::GID g = reads.first;
            g.outgroup = false;
            vargas::Graph subgraph = vargas::Graph(base_graph, pops.at(g));
            align_to_graph("i", subgraph, reads.second, aligners, out, threads);
        }
        std::cout << "Done." << std::endl;
    }
    if (O) {
        std::cout << "Aligning to Outgroup..." << std::endl;
        gdf.include_outgroups();
        for (auto &reads : read_origins) {
            /**
             * Align to the same graph that the read came from, but enforce
             * that it be the outgroup version.
             */
            std::cout << "\tGID: " << reads.first << std::endl;
            vargas::GID g = reads.first;
            g.outgroup = true;
            vargas::Graph subgraph = vargas::Graph(base_graph, pops.at(g));
            align_to_graph("o", subgraph, reads.second, aligners, out, threads);
        }
        std::cout << "Done." << std::endl;
    }
    return 0;
}

void align_to_graph(std::string label,
                    const vargas::Graph &subgraph,
                    const std::vector<vargas::Read> &reads,
                    const std::vector<std::shared_ptr<vargas::Aligner<>>> &aligners,
                    std::ostream &out,
                    unsigned int threads) {

    // Number of full parallel batches
    size_t batches = reads.size() / SIMDPP_FAST_INT8_SIZE;

    std::vector<std::thread> jobs;
    std::vector<std::vector<vargas::Alignment>> aligns(threads);

    // <= to catch remainder reads
    for (size_t i = 0; i <= batches; ++i) {
        auto end_iter = ((i + 1) * SIMDPP_FAST_INT8_SIZE) > reads.size() ?
                        reads.end() : reads.begin() + ((i + 1) * SIMDPP_FAST_INT8_SIZE);

        if (threads > 1) {
            jobs.emplace(jobs.end(),
                         &vargas::Aligner<>::align_into,
                         aligners[jobs.size()],
                         std::vector<vargas::Read>(reads.begin() + (i * SIMDPP_FAST_INT8_SIZE), end_iter),
                         subgraph.begin(),
                         subgraph.end(),
                         std::ref(aligns[jobs.size()]));

            // Once we use the allocated threads let them all finish.
            if (jobs.size() == threads || i == batches) {
                std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });

                // Print alignments
                for (size_t i = 0; i < jobs.size(); ++i) {
                    auto &aset = aligns[i];
                    for (auto &a : aset) {
                        out << label << ',' << a << std::endl;
                    }
                }

                jobs.clear();
            }
        } else {
            // Single threaded
            aligners[0]->align_into(std::vector<vargas::Read>(reads.begin() + (i * SIMDPP_FAST_INT8_SIZE), end_iter),
                                    subgraph.begin(),
                                    subgraph.end(),
                                    aligns[0]);
            for (auto &a : aligns[0]) {
                out << label << ',' << a << std::endl;
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
    std::string mut = "0", indel = "0", vnodes = "*", vbases = "*";
    std::string out_file, gdf_file;
    bool use_rate = false, outgroup = false;

    args >> GetOpt::Option('t', "outfile", out_file)
        >> GetOpt::Option('g', "gdef", gdf_file)
        >> GetOpt::OptionPresent('o', "outgroup", outgroup)
        >> GetOpt::Option('v', "vnodes", vnodes)
        >> GetOpt::Option('b', "vbases", vbases)
        >> GetOpt::Option('j', "threads", threads)
        >> GetOpt::Option('l', "rlen", read_len)
        >> GetOpt::Option('n', "numreads", num_reads)
        >> GetOpt::Option('m', "mut", mut)
        >> GetOpt::Option('i', "indel", indel)
        >> GetOpt::OptionPresent('a', "rate", use_rate);

    if (threads == 0) threads = std::thread::hardware_concurrency();

    std::ofstream out(out_file);
    if (!out.good()) throw std::invalid_argument("Error opening output file " + out_file);

    vargas::Gdef gdf(gdf_file);
    if (outgroup) gdf.include_outgroups();
    vargas::GraphBuilder gb(gdf.fasta(), gdf.var());
    gb.region(gdf.region());
    auto &pops = gdf.populations();

    std::cout << "Loading graph..." << std::endl;
    vargas::Graph base_graph = gb.build();

    auto mut_split = split(mut, ',');
    auto indel_split = split(indel, ',');
    auto vnode_split = split(vnodes, ',');
    auto vbase_split = split(vbases, ',');

    std::cout << "Building profiles..." << std::endl;
    std::vector<vargas::ReadProfile> pending_sims;
    for (auto &vbase : vbase_split) {
        for (auto &vnode : vnode_split) {
            for (auto &ind : indel_split) {
                for (auto &m : mut_split) {
                    vargas::ReadProfile prof;
                    prof.len = read_len;
                    prof.mut = m == "*" ? -1 : std::stof(m);
                    prof.indel = ind == "*" ? -1 : std::stof(ind);
                    prof.rand = use_rate;
                    prof.var_bases = vbase == "*" ? -1 : std::stoi(vbase);
                    prof.var_nodes = vnode == "*" ? -1 : std::stoi(vnode);
                    pending_sims.push_back(prof);
                }
            }
        }
    }
    std::cout << pending_sims.size() << " profiles, " << pops.size() << " subgraphs." << std::endl;

    // For each subgraph type
    for (auto &pop : pops) {
        std::cout << "\nSubgraph: " << pop.first;
        // Create subgraph
        vargas::Graph subgraph;
        subgraph = vargas::Graph(base_graph, pop.second);
        std::vector<std::thread> jobs;

        std::vector<std::shared_ptr<vargas::ReadSim>> sims;
        vargas::ReadSim rs_template(subgraph);
        for (unsigned int i = 0; i < threads; ++i) {
            sims.push_back(std::make_shared<vargas::ReadSim>(rs_template));
        }

        for (size_t p = 0; p < pending_sims.size(); ++p) {
            auto &prof = pending_sims[p];
            std::cout << "; (" << prof << ")" << std::flush;
            if (threads > 1) {
                sims[jobs.size()]->set_prof(prof);
                jobs.emplace(jobs.end(), &vargas::ReadSim::get_batch, sims[jobs.size()], num_reads);

                if (jobs.size() == threads || p == pending_sims.size() - 1) {
                    std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });

                    for (size_t i = 0; i < jobs.size(); ++i) {
                        auto &s = sims[i];
                        for (auto r : s->batch()) {
                            r.src = pop.first;
                            out << r << std::endl;
                        }
                    }
                    jobs.clear();
                }
            } else {
                // Single thread
                sims[0]->set_prof(prof);
                sims[0]->get_batch(num_reads);
                for (auto r : sims[0]->batch()) {
                    r.src = pop.first;
                    out << r << std::endl;
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

    std::string fasta_file = "", varfile = "", region = "", ingroups = "", outfile = "";
    int num = 1, node_len = 1000000;

    args >> GetOpt::Option('f', "fasta", fasta_file)
        >> GetOpt::Option('v', "var", varfile)
        >> GetOpt::Option('g', "region", region)
        >> GetOpt::Option('i', "ingroup", ingroups)
        >> GetOpt::Option('n', "num", num)
        >> GetOpt::Option('l', "nodelen", node_len)
        >> GetOpt::Option('t', "out", outfile);

    vargas::Gdef gdef(fasta_file, varfile, region, node_len);

    std::vector<std::string> ingroup_split = split(ingroups, ',');
    for (auto &ingrp_str : ingroup_split) {
        if (ingrp_str.at(0) == 'c') {
            // Explicit number of individuals
            size_t ingrp = std::stoi(ingrp_str.substr(1));
            for (int n = 0; n < num; ++n) {
                std::unordered_set<size_t> picked;
                vargas::Graph::Population filter(gdef.num_samples() * 2);
                while (picked.size() < ingrp) {
                    size_t ind = rand() % (gdef.num_samples() * 2);
                    if (picked.insert(ind).second) {
                        filter.set(ind);
                    }
                }
                gdef.add_population(vargas::GID(ingrp, n, false), filter);
            }
        } else {
            // Percentage
            int ingrp = std::stoi(ingrp_str);
            for (int n = 0; n < num; ++n) {
                vargas::Graph::Population filter(gdef.num_samples() * 2);
                for (unsigned int i = 0; i < filter.size(); ++i) {
                    if (rand() % 100 < ingrp) filter.set(i);
                }
                gdef.add_population(vargas::GID(ingrp, n, true), filter);
            }
        }
    }

    gdef.write(outfile);
    std::cout << gdef.populations().size() << " subgraph definitions generated." << std::endl;
    return 0;
}

void main_help() {
    using std::cout;
    using std::endl;
    cout << endl
        << "---------------------- vargas, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------\n";
    cout << "Operating modes \'vargas MODE\':" << endl;
    cout << "\tdefine      Define a set of graphs for use with sim/align.\n";
    cout << "\tsim         Simulate reads from a set of graphs.\n";
    cout << "\talign       Align reads to a set of graphs.\n";
    cout << "\ttest        Run doctests.\n";
    cout << "\tprofile     Run profiles.\n" << endl;
}

void define_help() {
    using std::cout;
    using std::endl;

    cout << endl
        << "-------------------- vargas define, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cout << "-f\t--fasta         *<string> Reference filename.\n";
    cout << "-v\t--var           *<string> VCF/BCF filename.\n";
    cout << "-t\t--out           *<string> Output filename.\n";
    cout << "-g\t--region        *<string> Region of graph, format CHR:MIN-MAX.\n";
    cout << "-l\t--nodelen       <int> Max node length, default 1,000,000\n";
    cout << "-i\t--ingroup       <int, int...> Percent ingroup subgraphs.\n";
    cout << "\t                  Prefix with \'c\' to interpret as number of individuals.\n";
    cout << "-n\t--num           <int> Number of unique graphs to create for all -i. Default 1.\n" << endl;
}

void profile_help() {
    using std::cout;
    using std::endl;
    cout << endl
        << "---------------------- vargas profile, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------" << endl;
    cout << "-f\t--fasta         *<string> Reference filename." << endl;
    cout << "-v\t--var           *<string> VCF/BCF filename." << endl;
    cout << "-g\t--region        *<string> Region of graph, format CHR:MIN-MAX." << endl;
    cout << "-i\t--ingroup       *<int> Percent of genotypes to include in alignment." << endl;
    cout << "-s\t--string        <string,string..> Include reads in alignment. Rest will be random." << endl << endl;
}

void align_help() {
    using std::cout;
    using std::endl;

    cout << endl
        << "------------------- vargas align, " << __DATE__ << ". rgaddip1@jhu.edu -------------------\n";;
    cout << "-g\t--gdef          *<string> Graph definition file.\n";
    cout << "-r\t--reads         *<string, string...> Read files to align.\n";
    cout << "-t\t--outfile       *<string> Alignment output file.\n";
    cout << "-l\t--rlen          <int> Max read length. Default 50.\n";
    cout << "-R\t                Align to reference graph.\n";
    cout << "-X\t                Align to maximum allele frequency graph.\n";
    cout << "-I\t                Align to ingroup graph.\n";
    cout << "-O\t                Align to outgroup graph.\n";
    cout << "-m\t--match         <int> Match score, default 2.\n";
    cout << "-n\t--mismatch      <int> Mismatch penalty, default 2.\n";
    cout << "-o\t--gap_open      <int> Gap opening penalty, default 3.\n";
    cout << "-e\t--gap_extend    <int> Gap extend penalty, default 1.\n";
    cout << "-j\t--threads       <int> Number of threads. 0 for maximum hardware concurrency.\n" << endl;

    cout << "Lines beginning with \'#\' are ignored.\n";
    cout << "Output format:\n\t[RXIO], read origin ingroup, read origin graph number, read sequence,\n"
        << " read origin pos, read sub errors, read indel errors, read variant nodes, read variant bases,\n"
        << " best score, best pos, best count, sub score, sub pos, sub count, corflag\n" << endl;
}

void sim_help() {
    using std::cout;
    using std::endl;

    cout << endl
        << "-------------------- vargas sim, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cout << "-g\t--gdef          *<string> Graph definition file. Reads are simulated from Ingroups.\n";
    cout << "-t\t--outfile       *<string> Alignment output file.\n";
    cout << "-o\t--outgroup      Simulate from outgroup graphs.\n";
    cout << "-n\t--numreads      <int> Number of reads to simulate from each subgraph, default 1000.\n";
    cout << "-m\t--muterr        <float, float...> Read mutation error. Default 0.\n";
    cout << "-i\t--indelerr      <float, float...> Read indel error. Default 0.\n";
    cout << "-v\t--vnodes        <int, int...> Number of variant nodes, default any (*).\n";
    cout << "-b\t--vbases        <int, int...> Number of variant bases, default any (*).\n";
    cout << "-l\t--rlen          <int> Read length, default 50.\n";
    cout << "-a\t--rate          Interpret -m, -i as rates, instead of exact number of errors.\n";
    cout << "-j\t--threads       <int> Number of threads. 0 for maximum hardware concurrency.\n" << endl;

    cout << "Comments preceded by \'#\'.\n";
    cout << "-n reads are produced for each -m, -i, -v, -b combination. If set to \'*\', any value is accepted."
        << endl << endl;
}