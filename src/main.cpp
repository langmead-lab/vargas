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
#include "gdef.h"
#include "main.h"
#include "getopt_pp.h"
#include "sam.h"


int main(const int argc, const char *argv[]) {

    srand(time(NULL)); // Rand used in profiles and sim.h

    if (argc > 1) {
        if (!strcmp(argv[1], "test")) {
            doctest::Context doc(argc, argv);
            doc.setOption("no-breaks", true);
            doc.setOption("abort-after", 10);
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
//            return align_main(argc, argv);
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

    Vargas::GraphBuilder gb(fasta, bcf);
    gb.region(region);
    if (!gb.good()) throw std::invalid_argument("Error opening files\n" + fasta + "\nand/or:\n" + bcf);


    std::clock_t start = std::clock();

    std::cerr << "Initial Build:\n\t";
    Vargas::Graph g;
    gb.build(g);
    std::vector<bool> filter;
    for (size_t i = 0; i < g.pop_size(); ++i) filter.push_back(rand() % 100 > 95);
    std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s, " << "Nodes: " << g.node_map()->size()
        << std::endl;

    size_t num = 0;

    {
        std::cerr << "Insertion order traversal:\n\t";
        start = std::clock();
        for (auto i = g.begin(); i != g.end(); ++i) {
            ++num;
        }
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        num = 0;
        std::cerr << "Filtering traversal, 100% in:\n\t";
        start = std::clock();

        for (auto i = g.begin(g.subset(100)); i != g.end(); ++i) {
            ++num;
        }
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        num = 0;
        std::cerr << "Filtering traversal, 5% in:\n\t";
        Vargas::Graph::Population filt(filter);
        start = std::clock();

        for (auto i = g.begin(filt); i != g.end(); ++i) {
            ++num;
        }
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s, " << "Nodes: " << num << std::endl;
    }

    {
        num = 0;
        std::cerr << "REF traversal:\n\t";
        start = std::clock();
        for (auto i = g.begin(Vargas::Graph::REF); i != g.end(); ++i) {
            ++num;
        }
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        num = 0;
        std::cerr << "MAXAF traversal:\n\t";
        start = std::clock();
        for (auto i = g.begin(Vargas::Graph::MAXAF); i != g.end(); ++i) { ;
        }
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        std::cerr << "Filter constructor:\n\t";
        auto pop_filt = g.subset(ingroup);
        start = std::clock();
        Vargas::Graph g2(g, pop_filt);
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        std::cerr << "REF constructor:\n\t";
        start = std::clock();
        Vargas::Graph g2(g, Vargas::Graph::REF);
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        std::cerr << "MAXAF constructor:\n\t";
        start = std::clock();
        Vargas::Graph g2(g, Vargas::Graph::MAXAF);
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        std::vector<std::string> reads;

        for (int i = 0; i < SIMDPP_FAST_INT8_SIZE; ++i) {
            std::ostringstream rd;
            for (int r = 0; r < 50; ++r) rd << rand_base();
            reads.push_back(rd.str());
        }

        std::vector<std::string> split_str;
        if (read.length() > 0) {
            split(read, ',', split_str);
            if (split_str.size() > SIMDPP_FAST_INT8_SIZE) split_str.resize(SIMDPP_FAST_INT8_SIZE);
            for (size_t i = 0; i < split_str.size(); ++i) reads[i] = split_str[i];
        }

        Vargas::Aligner<> a(g.max_node_len(), 50);

        std::cerr << SIMDPP_FAST_INT8_SIZE << " read alignment:\n";


        {
            start = std::clock();
            Vargas::Graph g2(g, g.subset(ingroup));
            std::cerr << "\tDerived Graph (" << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s)\n\t";
            start = std::clock();
            Vargas::Aligner<>::Results aligns = a.align(reads, g2.begin(), g2.end());
            std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
        }

    }
    return 0;
}

//int align_main(const int argc, const char *argv[]) {
//    GetOpt::GetOpt_pp args(argc, argv);
//
//    if (args >> GetOpt::OptionPresent('h', "help")) {
//        align_help();
//        return 0;
//    }
//
//    // Load parameters
//    uint8_t match = 2, mismatch = 2, gopen = 3, gext = 1;
//    unsigned int threads = 1, read_len = 50;
//    std::string read_file, gdf_file;
//    bool R = false, // Align to reference
//        X = false, // Align to max-AF graph
//        O = false, // Align to outgroup graph
//        I = false; // Align to ingroup graph
//
//    args >> GetOpt::Option('m', "match", match)
//        >> GetOpt::Option('n', "mismatch", mismatch)
//        >> GetOpt::Option('o', "gap_open", gopen)
//        >> GetOpt::Option('e', "gap_extend", gext)
//        >> GetOpt::Option('r', "reads", read_file)
//        >> GetOpt::Option('g', "gdef", gdf_file)
//        >> GetOpt::Option('j', "threads", threads)
//        >> GetOpt::Option('l', "rlen", read_len)
//        >> GetOpt::OptionPresent('R', R)
//        >> GetOpt::OptionPresent('X', X)
//        >> GetOpt::OptionPresent('O', O)
//        >> GetOpt::OptionPresent('I', I);
//
//    if (!R && !X && !O && !I) throw std::invalid_argument("No alignment mode (R,X,I,O) selected.");
//    if (threads == 0) threads = std::thread::hardware_concurrency();
//
//    std::ostream &out = std::cout;
//
//    Vargas::GDEF gdf(gdf_file);
//    Vargas::GraphBuilder gb(gdf.fasta(), gdf.var());
//    gb.region(gdf.region());
//    std::cerr << "Loading graph..." << std::endl;
//    Vargas::Graph base_graph = gb.build();
//
//    std::vector<std::shared_ptr<Vargas::Aligner<>>> aligners;
//    for (uint8_t i = 0; i < threads; ++i)
//        aligners.push_back(std::make_shared<Vargas::Aligner<>>(base_graph.max_node_len(),
//                                                               read_len,
//                                                               match,
//                                                               mismatch,
//                                                               gopen,
//                                                               gext));
//
//    // Load reads. Maps graph origin label to a vector of reads
//    std::cerr << "Loading reads..." << std::endl;
//    std::map<Vargas::Graph::GID, std::vector<Vargas::Read>> read_origins;
//    if (read_file.length() == 0) {
//        // Use stdin
//        Vargas::ReadFile reads;
//        while (reads.update_read()) {
//            Vargas::Read r = reads.get_read();
//            read_origins[r.src].push_back(r);
//        }
//    }
//    else {
//        auto readfile_split = split(read_file, ',');
//        for (auto &rfile : readfile_split) {
//            Vargas::ReadFile reads(rfile);
//            while (reads.update_read()) {
//                Vargas::Read r = reads.get_read();
//                read_origins[r.src].push_back(r);
//            }
//        }
//    }
//
//    auto &pops = gdf.populations();
//
//    if (R) {
//        std::cerr << "Aligning to reference..." << std::endl;
//        Vargas::Graph subgraph = Vargas::Graph(base_graph, Vargas::Graph::REF);
//        for (auto &reads : read_origins) {
//            std::cerr << "\tGraph::GID: " << reads.first << std::endl;
//            align_to_graph("r", subgraph, reads.second, aligners, out, threads);
//        }
//        std::cerr << "Done." << std::endl;
//    }
//    if (X) {
//        std::cerr << "Aligning to MAX AF..." << std::endl;
//        Vargas::Graph subgraph = Vargas::Graph(base_graph, Vargas::Graph::MAXAF);
//        for (auto &reads : read_origins) {
//            std::cerr << "\tGraph::GID: " << reads.first << std::endl;
//            align_to_graph("x", subgraph, reads.second, aligners, out, threads);
//        }
//        std::cerr << "Done." << std::endl;
//    }
//    if (I) {
//        std::cerr << "Aligning to Ingroup..." << std::endl;
//        for (auto &reads : read_origins) {
//            /**
//             * Align to the same graph that the read came from, but enforce
//             * that it be the ingroup version.
//             */
//            std::cerr << "\tGraph::GID: " << reads.first << std::endl;
//            Vargas::Graph::GID g = reads.first;
//            g.outgroup = false;
//            Vargas::Graph subgraph = Vargas::Graph(base_graph, pops.at(g));
//            align_to_graph("i", subgraph, reads.second, aligners, out, threads);
//        }
//        std::cerr << "Done." << std::endl;
//    }
//    if (O) {
//        std::cerr << "Aligning to Outgroup..." << std::endl;
//        gdf.include_outgroups();
//        for (auto &reads : read_origins) {
//            /**
//             * Align to the same graph that the read came from, but enforce
//             * that it be the outgroup version.
//             */
//            std::cerr << "\tGraph::GID: " << reads.first << std::endl;
//            Vargas::Graph::GID g = reads.first;
//            g.outgroup = true;
//            Vargas::Graph subgraph = Vargas::Graph(base_graph, pops.at(g));
//            align_to_graph("o", subgraph, reads.second, aligners, out, threads);
//        }
//        std::cerr << "Done." << std::endl;
//    }
//    return 0;
//}

void align_to_graph(std::string label,
                    const Vargas::Graph &subgraph,
                    const std::vector<std::string> &reads,
                    const std::vector<std::shared_ptr<Vargas::Aligner<>>> &aligners,
                    std::ostream &out,
                    unsigned int threads) {



}

int sim_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        sim_help();
        return 0;
    }

    // Load parameters
    unsigned int read_len = 50, num_reads = 1000;
    std::string mut = "0", indel = "0", vnodes = "*", vbases = "*";
    std::string gdf_file, out_file;
    bool use_rate = false, outgroup = false;

    Vargas::SAM::Header sam_hdr;

    {
        Vargas::SAM::Header::Program pg;
        std::ostringstream ss;
        ss << "vargas sim ";
        for (int i = 0; i < argc; ++i) ss << std::string(argv[i]) << " ";
        pg.command_line = ss.str();
        pg.name = "vargas_sim";
        pg.id = "VS";
        pg.version = __DATE__;
        std::replace_if(pg.version.begin(), pg.version.end(), isspace, ' '); // rm tabs
        sam_hdr.add(pg);
    }

    args >> GetOpt::Option('g', "gdef", gdf_file)
        >> GetOpt::OptionPresent('o', "outgroup", outgroup)
        >> GetOpt::Option('v', "vnodes", vnodes)
        >> GetOpt::Option('b', "vbases", vbases)
        >> GetOpt::Option('l', "rlen", read_len)
        >> GetOpt::Option('n', "numreads", num_reads)
        >> GetOpt::Option('m', "mut", mut)
        >> GetOpt::Option('i', "indel", indel)
        >> GetOpt::Option('t', "out", out_file)
        >> GetOpt::OptionPresent('a', "rate", use_rate);


    time_t t = time(0);
    struct tm *now = localtime(&t);
    std::string curr_date;
    {
        std::ostringstream ss;
        ss << (now->tm_year + 1900) << '-'
            << (now->tm_mon + 1) << '-'
            << now->tm_mday;
        curr_date = ss.str();
    }

    Vargas::GDEF gdf(gdf_file);
    if (outgroup) gdf.include_outgroups();
    Vargas::GraphBuilder gb(gdf.fasta(), gdf.var());
    gb.region(gdf.region());
    auto &pops = gdf.populations();

    auto mut_split = split(mut, ',');
    auto indel_split = split(indel, ',');
    auto vnode_split = split(vnodes, ',');
    auto vbase_split = split(vbases, ',');

    std::cerr << "Building profiles..." << std::endl;
    std::vector<std::pair<int, Vargas::Sim::Profile>> pending_sims;
    int rg_id = 0;
    Vargas::SAM::Header::ReadGroup rg;
    rg.seq_center = "vargas_sim";
    rg.date = curr_date;
    rg.aux.set(SIM_SAM_REF_TAG, gdf.fasta());
    rg.aux.set(SIM_SAM_VCF_TAG, gdf.var());
    for (auto &vbase : vbase_split) {
        for (auto &vnode : vnode_split) {
            for (auto &ind : indel_split) {
                for (auto &m : mut_split) {
                    Vargas::Sim::Profile prof;
                    prof.len = read_len;
                    prof.mut = m == "*" ? -1 : std::stof(m);
                    prof.indel = ind == "*" ? -1 : std::stof(ind);
                    prof.rand = use_rate;
                    prof.var_bases = vbase == "*" ? -1 : std::stoi(vbase);
                    prof.var_nodes = vnode == "*" ? -1 : std::stoi(vnode);

                    pending_sims.push_back(std::pair<int, Vargas::Sim::Profile>(++rg_id, prof));

                    rg.aux.set(SIM_SAM_INDEL_ERR_TAG, prof.indel);
                    rg.aux.set(SIM_SAM_VAR_NODES_TAG, prof.var_nodes);
                    rg.aux.set(SIM_SAM_VAR_BASE_TAG, prof.var_bases);
                    rg.aux.set(SIM_SAM_SUB_ERR_TAG, prof.mut);
                    rg.aux.set(SIM_SAM_USE_RATE_TAG, (int) prof.rand);


                    // Each profile and subgraph combination is a unique set of reads
                    for (auto p : pops) {
                        rg.aux.set(SIM_SAM_GID_TAG, p.first.to_string());
                        rg.aux.set(SIM_SAM_POPULATION, p.second.to_string());
                        rg.id = std::to_string(rg_id) + "_" + p.first.to_string();
                        sam_hdr.add(rg);
                    }
                }
            }
        }
    }
    std::cerr << pending_sims.size() << " profile(s), " << pops.size() << " subgraphs." << std::endl;

    Vargas::osam out(out_file, sam_hdr);
    if (!out.good()) throw std::invalid_argument("Error opening output file \"" + out_file + "\"");

    std::cerr << "Loading graph..." << std::endl;
    Vargas::Graph base_graph = gb.build();

    // For each subgraph type
    for (auto &pop : pops) {
        std::cerr << "\nSubgraph: " << pop.first;

        // Create subgraph
        Vargas::Graph subgraph;
        subgraph = Vargas::Graph(base_graph, pop.second);

        Vargas::Sim sim(subgraph);

        Vargas::SAM::Record rec;
        rec.flag.unmapped = true;

        for (size_t p = 0; p < pending_sims.size(); ++p) {

            sim.set_prof(pending_sims[p].second);
            sim.get_batch(num_reads);

            for (auto r : sim.batch()) {
                rec.seq = r.read;
                rec.aux.set(SIM_SAM_INDEL_ERR_TAG, r.indel_err);
                rec.aux.set(SIM_SAM_VAR_NODES_TAG, r.var_nodes);
                rec.aux.set(SIM_SAM_VAR_BASE_TAG, r.var_bases);
                rec.aux.set(SIM_SAM_SUB_ERR_TAG, r.sub_err);
                rec.aux.set(SIM_SAM_END_POS_TAG, r.end_pos);
                rec.aux.set(SIM_SAM_INDIV_TAG, r.indiv);
                rec.aux.set(SIM_SAM_READ_ORIG_TAG, r.read_orig);
                rec.aux.set(SIM_SAM_VAR_BASE_TAG, r.var_bases);
                rec.aux.set("RG", std::to_string(pending_sims[p].first) + "_" + pop.first.to_string());
                out.add_record(rec);
            }

        }
    }
    std::cerr << std::endl;
    return 0;
}

int define_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        define_help();
        return 0;
    }

    std::string fasta_file = "", varfile = "", region = "", ingroups = "";
    int num = 1, node_len = 1000000;

    args >> GetOpt::Option('f', "fasta", fasta_file)
        >> GetOpt::Option('v', "var", varfile)
        >> GetOpt::Option('g', "region", region)
        >> GetOpt::Option('i', "ingroup", ingroups)
        >> GetOpt::Option('n', "num", num)
        >> GetOpt::Option('l', "nodelen", node_len);

    Vargas::GDEF gdef(fasta_file, varfile, region, node_len);

    std::vector<std::string> ingroup_split = split(ingroups, ',');
    for (auto &ingrp_str : ingroup_split) {
        if (ingrp_str.at(0) == 'c') {
            // Explicit number of individuals
            size_t ingrp = std::stoi(ingrp_str.substr(1));
            for (int n = 0; n < num; ++n) {
                std::unordered_set<size_t> picked;
                Vargas::Graph::Population filter(gdef.num_samples() * 2);
                while (picked.size() < ingrp) {
                    size_t ind = rand() % (gdef.num_samples() * 2);
                    if (picked.insert(ind).second) {
                        filter.set(ind);
                    }
                }
                gdef.add_population(Vargas::Graph::GID(ingrp, n, false), filter);
            }
        } else {
            // Percentage
            int ingrp = std::stoi(ingrp_str);
            for (int n = 0; n < num; ++n) {
                Vargas::Graph::Population filter(gdef.num_samples() * 2);
                for (unsigned int i = 0; i < filter.size(); ++i) {
                    if (rand() % 100 < ingrp) filter.set(i);
                }
                gdef.add_population(Vargas::Graph::GID(ingrp, n, true), filter);
            }
        }
    }

    gdef.write(std::cout);
    std::cerr << gdef.populations().size() << " subgraph definitions generated." << std::endl;
    return 0;
}

void main_help() {
    using std::cerr;
    using std::endl;
    cerr << endl
        << "---------------------- Vargas, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------\n";
    cerr << "Operating modes \'Vargas MODE\':" << endl;
    cerr << "\tdefine      Define a set of graphs for use with sim/align.\n";
    cerr << "\tsim         Simulate reads from a set of graphs.\n";
    cerr << "\talign       Align reads to a set of graphs.\n";
    cerr << "\ttest        Run doctests.\n";
    cerr << "\tprofile     Run profiles.\n" << endl;
}

void define_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
        << "-------------------- Vargas define, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-f\t--fasta         *<string> Reference filename.\n";
    cerr << "-v\t--var           *<string> VCF/BCF filename.\n";
    cerr << "-g\t--region        *<string> Region of graph, format CHR:MIN-MAX.\n";
    cerr << "-l\t--nodelen       <int> Max node length, default 1,000,000\n";
    cerr << "-i\t--ingroup       <int, int...> Percent ingroup subgraphs.\n";
    cerr << "\t                  Prefix with \'c\' to interpret as number of individuals.\n";
    cerr << "-n\t--num           <int> Number of unique graphs to create for all -i. Default 1.\n" << endl;

    cerr << "Definition output to stdout.\n" << endl;
}

void profile_help() {
    using std::cerr;
    using std::endl;
    cerr << endl
        << "---------------------- Vargas profile, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------" << endl;
    cerr << "-f\t--fasta         *<string> Reference filename." << endl;
    cerr << "-v\t--var           *<string> VCF/BCF filename." << endl;
    cerr << "-g\t--region        *<string> Region of graph, format CHR:MIN-MAX." << endl;
    cerr << "-i\t--ingroup       *<int> Percent of genotypes to include in alignment." << endl;
    cerr << "-s\t--string        <string,string..> Include reads in alignment. Rest will be random." << endl << endl;
}

void align_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
        << "------------------- Vargas align, " << __DATE__ << ". rgaddip1@jhu.edu -------------------\n";;
    cerr << "-g\t--gdef          *<string> Graph definition file.\n";
    cerr << "-r\t--reads         <string, string...> Read files to align. Default stdin.\n";
    cerr << "-l\t--rlen          <int> Max read length. Default 50.\n";
    cerr << "-R\t                Align to reference graph.\n";
    cerr << "-X\t                Align to maximum allele frequency graph.\n";
    cerr << "-I\t                Align to ingroup graph.\n";
    cerr << "-O\t                Align to outgroup graph.\n";
    cerr << "-m\t--match         <int> Match score, default 2.\n";
    cerr << "-n\t--mismatch      <int> Mismatch penalty, default 2.\n";
    cerr << "-o\t--gap_open      <int> Gap opening penalty, default 3.\n";
    cerr << "-e\t--gap_extend    <int> Gap extend penalty, default 1.\n";
    cerr << "-j\t--threads       <int> Number of threads. 0 for maximum hardware concurrency.\n" << endl;

    cerr << "Lines beginning with \'#\' are ignored. Alignments output to stdout, reads input from stdin\n";
    cerr << "if no -r specified.\n";
    cerr << "Output format:\n\t[RXIO], read origin ingroup, read origin graph number, read sequence,\n"
        << " read origin pos, read sub errors, read indel errors, read variant nodes, read variant bases,\n"
        << " best score, best pos, best count, sub score, sub pos, sub count, corflag\n" << endl;
}

void sim_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
        << "-------------------- Vargas sim, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-g\t--gdef          *<string> Graph definition file. Reads are simulated from Ingroups.\n";
    cerr << "-t\t--out           <string> Output file. Default stdout.\n";
    cerr << "-o\t--outgroup      Simulate from outgroup graphs.\n";
    cerr << "-n\t--numreads      <int> Number of reads to simulate from each subgraph, default 1000.\n";
    cerr << "-m\t--muterr        <float, float...> Read mutation error. Default 0.\n";
    cerr << "-i\t--indelerr      <float, float...> Read indel error. Default 0.\n";
    cerr << "-v\t--vnodes        <int, int...> Number of variant nodes, default any (*).\n";
    cerr << "-b\t--vbases        <int, int...> Number of variant bases, default any (*).\n";
    cerr << "-l\t--rlen          <int> Read length, default 50.\n";
    cerr << "-a\t--rate          Interpret -m, -i as rates, instead of exact number of errors.\n" << endl;

    cerr << "-n reads are produced for each -m, -i, -v, -b combination. If set to \'*\', any value is accepted."
        << endl << endl;
}