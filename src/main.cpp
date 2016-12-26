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

#define DOCTEST_CONFIG_IMPLEMENT // User controlled test execution

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gdef.h"
#include "main.h"
#include "align_main.h"
#include <iostream>
#include <thread>
#include <algorithm>


int main(int argc, char *argv[]) {

    srand(time(0)); // Rand used in profiles and sim

    try {
        if (argc > 1) {
            if (!strcmp(argv[1], "test")) {
                doctest::Context doc(argc, argv);
                doc.setOption("abort-after", 10);
                return doc.run();
            } else if (!strcmp(argv[1], "profile")) {
                return profile(argc, argv);
            } else if (!strcmp(argv[1], "define")) {
                return define_main(argc, argv);
            } else if (!strcmp(argv[1], "sim")) {
                return sim_main(argc, argv);
            } else if (!strcmp(argv[1], "align")) {
                return align_main(argc, argv);
            } else if (!strcmp(argv[1], "export")) {
                return export_main(argc, argv);
            } else if (!strcmp(argv[1], "convert")) {
                return convert_main(argc, argv);
            } else if (!strcmp(argv[1], "query")) {
                return query_main(argc, argv);
            }
        }
    } catch (std::exception &e) {
        std::cerr << "\033[1;31m"
                  << "\nFatal Error: " << e.what()
                  << "\033[0m\n" << std::endl;
        return 1;
    }

    std::cerr << "Define a valid mode of operation." << std::endl;
    main_help();
    return 1;

}

int define_main(int argc, char *argv[]) {
    std::string fasta_file, varfile, region, subgraph_def, out_file, dot_file, sample_filter;
    bool invert_filter = false;
    int node_len;

    cxxopts::Options opts("vargas define", "Define subgraphs deriving from a reference and VCF file.");
    try {
        opts.add_options()
        ("f,fasta", "<str> *Reference FASTA.", cxxopts::value(fasta_file))
        ("v,vcf", "<str> VCF/BCF File.", cxxopts::value<std::string>(varfile))
        ("g,region", "<str> Region of format \"CHR:MIN-MAX\". \"CHR:0-0\" for all.", cxxopts::value(region))
        ("s,subgraph", "<str> File or definitions of format \"[~]NAME=N[%t],...\".",
         cxxopts::value(subgraph_def))
        ("l,nodelen", "<N> Maximum node length.", cxxopts::value(node_len)->default_value("100000"))
        ("p,filter", "<str> Filter by sample names in file.", cxxopts::value(sample_filter))
        ("x,invert", "Invert sample filter, exclude -p samples.", cxxopts::value(invert_filter))
        ("t,out", "<str> Output filename. (default: stdout)", cxxopts::value(out_file))
        ("d,dot", "<str> Export hierarchy to a DOT file.", cxxopts::value(dot_file))
        ("h,help", "Display this message.");
        opts.parse(argc, argv);
    } catch (std::exception &e) { throw std::invalid_argument("Error parsing options."); }
    if (opts.count("h")) {
        define_help(opts);
        return 0;
    }
    if (!opts.count("f")) throw std::invalid_argument("FASTA file required.");

    // Load subgraph defs
    std::string subgraph_str;
    if (subgraph_def.length() == 0) {
        subgraph_str = "";
    } else if (vargas::GraphManager::is_definition(subgraph_def)) {
        subgraph_str = subgraph_def;
    } else {
        std::ifstream in(subgraph_def);
        if (!in.good()) throw std::invalid_argument("Error opening file \"" + subgraph_def + "\".");
        std::stringstream buff;
        buff << in.rdbuf();
        subgraph_str = buff.str();
    }

    vargas::GraphManager gm;
    if (sample_filter.length()) {
        std::ifstream in(sample_filter);
        if (!in.good()) throw std::invalid_argument("Error opening file: \"" + sample_filter + "\"");
        std::stringstream ss;
        ss << in.rdbuf();
        sample_filter = ss.str();
    }
    gm.set_filter(sample_filter, invert_filter);
    gm.write(fasta_file, varfile, region, subgraph_str, node_len, out_file, false);
    if (dot_file.length() > 0) gm.to_DOT(dot_file, "subgraphs");

    std::cerr << gm.size() << " subgraph definitions generated." << std::endl;

    return 0;
}

int sim_main(int argc, char *argv[]) {
    std::string cl;
    {
        std::ostringstream ss;
        for (int i = 0; i < argc; ++i) ss << std::string(argv[i]) << " ";
        cl = ss.str();
    }

    int read_len, num_reads, threads;
    std::string mut, indel, vnodes, vbases, gdf_file, out_file, sim_src;
    bool use_rate = false, sim_src_isfile = false;

    cxxopts::Options opts("vargas sim", "Simulate reads from genome graphs.");
    try {
        opts.add_options()
        ("g,gdef", "<str> *Graph definition file.", cxxopts::value(gdf_file))
        ("s,sub", "<S1,S2..> Subgraphs to simulate from. (default: all)", cxxopts::value(sim_src))
        ("l,rlen", "<N> Read length.", cxxopts::value(read_len)->default_value("50"))
        ("n,numreads", "<N> Number of reads to generate.", cxxopts::value(num_reads)->default_value("1000"))
        ("f,file", "-s specifies a filename.", cxxopts::value(sim_src_isfile))
        ("d,vnodes", "<N1,N2...> Number of variant nodes. \'*\' for any.",
         cxxopts::value(vnodes)->default_value("*"))
        ("b,vbases", "<N1,N2...> Number of variant bases. \'*\' for any.",
         cxxopts::value(vbases)->default_value("*"))
        ("m,mut", "<N1,N2...> Number of mutations. \'*\' for any.", cxxopts::value(mut)->default_value("0"))
        ("i,indel", "<N1,N2...> Number of insertions/deletions. \'*\' for any.",
         cxxopts::value(indel)->default_value("0"))
        ("a,rate", "<N1,N2...> Interpret -m, -i as error rates.", cxxopts::value(use_rate))
        ("t,out", "<str> Output file. (default: stdout)", cxxopts::value(out_file))
        ("j,threads", "<N> Number of threads.", cxxopts::value(threads)->default_value("1"))
        ("h,help", "Display this message.");
        opts.parse(argc, argv);
    } catch (std::exception &e) { throw std::invalid_argument("Error parsing options."); }
    if (opts.count("h")) {
        sim_help(opts);
        return 0;
    }
    if (!opts.count("g")) throw std::invalid_argument("Graph definition file required.");

    #ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
    #endif

    vargas::SAM::Header sam_hdr;

    {
        vargas::SAM::Header::Program pg;
        pg.command_line = cl;
        pg.name = "vargas_sim";
        pg.id = "VS";
        pg.version = __DATE__;
        std::replace_if(pg.version.begin(), pg.version.end(), isspace, ' '); // rm tabs
        sam_hdr.add(pg);
    }

    vargas::GraphManager gm;

    const std::vector<std::string>
    mut_split = rg::split(mut, ','),
    indel_split = rg::split(indel, ','),
    vnode_split = rg::split(vnodes, ','),
    vbase_split = rg::split(vbases, ',');

    if (sim_src_isfile) {
        std::ifstream in(sim_src);
        if (!in.good()) throw std::invalid_argument("Error opening file: \"" + sim_src + "\".");
        std::stringstream ss;
        ss << in.rdbuf();
        sim_src = ss.str();
    }

    std::cerr << "Loading base graph... " << std::flush;
    auto start_time = std::chrono::steady_clock::now();
    gm.open(gdf_file);
    std::cerr << rg::chrono_duration(start_time) << " seconds." << std::endl;
    std::cerr << gdf_file
              << "- FASTA: " << gm.reference()
              << " VCF: " << gm.variants()
              << " REGION: " << gm.region() << std::endl;

    std::vector<std::string> subdef_split;
    if (sim_src.length() == 0) {
        // Use all subgraphs. If there are none, use the base graph
        subdef_split = gm.labels();
        if (subdef_split.size() == 0) subdef_split.push_back(vargas::GraphManager::GDEF_BASEGRAPH);
    } else {
        std::replace(sim_src.begin(), sim_src.end(), '\n', gm.GDEF_DELIM);
        sim_src.erase(std::remove_if(sim_src.begin(), sim_src.end(), isspace), sim_src.end());
        subdef_split = rg::split(sim_src, gm.GDEF_DELIM);

        // validate graph labels
        for (const std::string &l : subdef_split) {
            gm.filter(l); // Will throw if l does not exist
        }

    }

    std::cerr << "Building profiles... " << std::flush;
    start_time = std::chrono::steady_clock::now();

    // Map a graph label to a vector of read group ID's/profiles
    std::unordered_map<std::string, // Graph label
                       std::vector<std::pair<std::string, // Read Group ID
                                             vargas::Sim::Profile>>> // Sim profile
    queue;

    int rg_id = 0;
    vargas::SAM::Header::ReadGroup rg;
    rg.seq_center = "vargas_sim";
    rg.date = rg::current_date();
    rg.aux.set(SIM_SAM_REF_TAG, gm.reference());
    rg.aux.set(SIM_SAM_VCF_TAG, gm.variants());

    vargas::Sim::Profile prof;
    prof.len = read_len;
    prof.rand = use_rate;
    for (const std::string &vbase : vbase_split) {
        for (const std::string &vnode : vnode_split) {
            for (const std::string &ind : indel_split) {
                for (const std::string &m : mut_split) {
                    try {
                        prof.mut = m == "*" ? -1 : std::stof(m);
                        prof.indel = ind == "*" ? -1 : std::stof(ind);
                        prof.var_bases = vbase == "*" ? -1 : std::stoi(vbase);
                        prof.var_nodes = vnode == "*" ? -1 : std::stoi(vnode);
                    } catch (std::exception &e) {
                        std::cerr << "Invalid profile argument: " << e.what() << std::endl;
                        return 1;
                    }

                    rg.aux.set(SIM_SAM_INDEL_ERR_TAG, prof.indel);
                    rg.aux.set(SIM_SAM_VAR_NODES_TAG, prof.var_nodes);
                    rg.aux.set(SIM_SAM_VAR_BASE_TAG, prof.var_bases);
                    rg.aux.set(SIM_SAM_SUB_ERR_TAG, prof.mut);
                    rg.aux.set(SIM_SAM_USE_RATE_TAG, (int) prof.rand);

                    // Each profile and subgraph combination is a unique set of reads
                    for (const std::string &p : subdef_split) {
                        rg.aux.set(SIM_SAM_SRC_TAG, p);
                        rg.aux.set(SIM_SAM_POPULATION, gm.filter(p).to_string());
                        rg.id = std::to_string(++rg_id);
                        sam_hdr.add(rg);
                        queue[p].emplace(queue[p].end(), rg.id, prof);
                    }
                }
            }
        }
    }
    std::cerr << rg::chrono_duration(start_time)
              << " seconds. "
              << sam_hdr.read_groups.size() << " read groups over "
              << subdef_split.size() << " subgraphs. " << std::endl;

    vargas::osam out(out_file, sam_hdr);
    if (!out.good()) throw std::invalid_argument("Error opening output file \"" + out_file + "\"");

    std::cerr << "Simulating... " << std::flush;

    start_time = std::chrono::steady_clock::now();

    std::vector<std::pair<std::string, // Graph label
                          std::pair<std::string, // RG ID
                                    vargas::Sim::Profile>>> // sim prof
    task_list;

    for (size_t k = 0; k < subdef_split.size(); ++k) {
        for (size_t i = 0; i < queue.at(subdef_split[k]).size(); ++i) {
            auto &p = queue.at(subdef_split[k]).at(i);
            task_list.push_back(std::pair<std::string, std::pair<std::string, vargas::Sim::Profile>>
                                (subdef_split[k], std::pair<std::string, vargas::Sim::Profile>(p.first, p.second)));
        }
    }

    const size_t num_tasks = task_list.size();
    #pragma omp parallel for
    for (size_t n = 0; n < num_tasks; ++n) {
        const std::string label = task_list.at(n).first;
        const auto subgraph_ptr = gm.make_subgraph(label);
        vargas::Sim sim(*subgraph_ptr, task_list.at(n).second.second);
        auto results = sim.get_batch(num_reads);
        gm.destroy(label);
        for (auto &r : results) r.aux.set("RG", task_list.at(n).second.first);
        #pragma omp critical(sam_out)
        {
            for (const auto &r : results) out.add_record(r);
        }
    }

    std::cerr << rg::chrono_duration(start_time) << " seconds." << std::endl;

    return 0;
}

int convert_main(int argc, char **argv) {
    std::string sam_file, format;

    cxxopts::Options opts("vargas convert", "Export a SAM file as a CSV file.");
    try {
        opts.add_options()
        ("f,format", "<str> *Output format.", cxxopts::value<std::string>(format))
        ("s,sam", "<str> SAM file. (default: stdin)", cxxopts::value<std::string>(sam_file))
        ("h,help", "Display this message.");
        opts.parse(argc, argv);
    } catch (std::exception &e) { throw std::invalid_argument("Error parsing options: " + std::string(e.what())); }
    if (opts.count("h")) {
        convert_help(opts);
        return 0;
    }
    if (!opts.count("f")) throw std::invalid_argument("Format specifier required.");

    auto start_time = std::chrono::steady_clock::now();

    format.erase(std::remove(format.begin(), format.end(), ' '), format.end());
    std::vector<std::string> fmt_split = rg::split(format, ',');
    std::unordered_set<std::string> warned;
    vargas::isam input(sam_file);
    std::string buff, val;
    do {
        buff = "";
        for (auto &tag : fmt_split) {
            val = "*";
            if (!input.record().get(input.header(), tag, val) && warned.count(tag) == 0) {
                std::cerr << "WARN: Tag \"" << tag << "\" not present." << std::endl;
            }
            buff += val;
            buff += ",";
        }
        std::cout << buff.substr(0, buff.length() - 1) << '\n'; // Crop leading comma
    } while (input.next());

    std::cerr << std::chrono::duration_cast<std::chrono::duration<double>>(
    std::chrono::steady_clock::now() - start_time).count()
              << " seconds." << std::endl;

    return 0;
}

int profile(int argc, char *argv[]) {
    std::string bcf, fasta, region;
    size_t nreads, read_len, ingroup;

    cxxopts::Options opts("vargas convert", "Export a SAM file as a CSV file.");
    try {
        opts.add_options()
        ("f,fasta", "<str> *Reference FASTA.", cxxopts::value(fasta))
        ("v,vcf", "<str> *Variant File.", cxxopts::value(bcf))
        ("g,region", "<str> *Region of format \"CHR:MIN-MAX\". \"CHR:0-0\" for all.", cxxopts::value(region))
        ("i,ingroup", "<N> Ingroup percentage.", cxxopts::value(ingroup)->default_value("100"))
        ("n,nreads", "<N> Number of reads.", cxxopts::value(nreads)->default_value("32"))
        ("l,len", "<N> Number of reads.", cxxopts::value(read_len)->default_value("50"))
        ("h,help", "Display this message.");
        opts.parse(argc, argv);
    } catch (std::exception &e) { throw std::invalid_argument("Error parsing options: " + std::string(e.what())); }
    if (opts.count("h")) {
        profile_help(opts);
        return 0;
    }
    if (!opts.count("f")) throw std::invalid_argument("FASTA file required.");


    vargas::GraphFactory gb(fasta);
    gb.open_vcf(bcf);
    gb.set_region(region);

    auto start = std::clock();

    std::cerr << "Initial Graph Build:\n\t";
    vargas::Graph g;
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
        std::cerr << "Filter constructor (" << ingroup << "):\n\t";
        auto pop_filt = g.subset(ingroup);
        start = std::clock();
        vargas::Graph g2(g, pop_filt);
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        std::cerr << "REF constructor:\n\t";
        start = std::clock();
        vargas::Graph g2(g, vargas::Graph::Type::REF);
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }

    {
        std::cerr << "MAXAF constructor:\n\t";
        start = std::clock();
        vargas::Graph g2(g, vargas::Graph::Type::MAXAF);
        std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    }
    return 0;
}

int export_main(int argc, char *argv[]) {

    // Load parameters
    std::string subgraph, file, out;

    cxxopts::Options opts("vargas export", "Export a graph in DOT format.");
    try {
        opts.add_options()
        ("g,gdef", "<str> *Graph definition file. (default: stdin)", cxxopts::value(file))
        ("s,subgraph", "<str> Subgraph to export.",
         cxxopts::value(subgraph)->default_value(vargas::GraphManager::GDEF_BASEGRAPH))
        ("t,out", "<str> Output file. (default: stdout)", cxxopts::value(out))
        ("h,help", "Display this message.");
        opts.parse(argc, argv);
    } catch (std::exception &e) { throw std::invalid_argument("Error parsing options: " + std::string(e.what())); }
    if (opts.count("h")) {
        export_help(opts);
        return 0;
    }

    vargas::GraphManager gm;
    if (file.size()) gm.open(file);
    else gm.open(std::cin);
    auto g = gm.make_subgraph(subgraph);
    if (out.size()) g->to_DOT(out, "g");
    std::cout << g->to_DOT();
    return 0;
}

int query_main(int argc, char *argv[]) {
    std::string region, gdef, fasta, vcf, subgraph, out;

    cxxopts::Options opts("vargas query", "Query VCF, Graph, or FASTA files.");
    try {
        opts.add_options()
        ("d,gdef", "<str> Export graph region as a DOT file.", cxxopts::value(gdef))
        ("a,stat", "Print statistics about subgraphs.")
        ("s,subgraph", "<str> Subgraph to export.",
         cxxopts::value(subgraph)->default_value(vargas::GraphManager::GDEF_BASEGRAPH))
        ("g,region", "<str> *Region of format \"CHR:MIN-MAX\". \"CHR:0-0\" for all.", cxxopts::value(region))
        ("t,out", "<str> -d output file. (default: stdout)", cxxopts::value(out))
        ("f,fasta", "<str> Reference FASTA.", cxxopts::value(fasta))
        ("v,vcf", "<str> Variant File.", cxxopts::value(vcf))
        ("h,help", "Display this message.");
        opts.parse(argc, argv);
    } catch (std::exception &e) { throw std::invalid_argument("Error parsing options: " + std::string(e.what())); }
    if (opts.count("h")) {
        query_help(opts);
        return 0;
    }

    if (opts.count("a")) {
        vargas::GraphManager gm(gdef);

        auto q = gm.labels();

        #pragma omp parallel for
        for (size_t k = 0; k < q.size(); ++k) {
            const std::string s = q[k];
            auto res = gm.make_subgraph(s)->statistics();

            #pragma omp critical
            {
                std::cerr << s << " Graph"
                          << "\n\tGraph Nodes: " << res.num_nodes
                          << "\n\tTotal Length: " << res.total_length
                          << "\n\tNumber of Edges: " << res.num_edges
                          << "\n\tNumber of SNPs: " << res.num_snps
                          << "\n\tNumber of Deletions: " << res.num_dels
                          << std::endl;
                gm.destroy(s);
            }
        }
        return 0;
    }


    if (!opts.count("g")) throw std::invalid_argument("Region specifier required.");
    auto reg = vargas::parse_region(region);

    if (gdef.size()) {
        vargas::GraphManager gm(gdef);

        auto sg = gm.make_subgraph(subgraph)->subgraph(reg.min, reg.max);
        if (out.length() == 0) std::cout << sg.to_DOT();
        else {
            std::ofstream o(out);
            o << sg.to_DOT();
        }

        std::cout << std::endl;

        vargas::ifasta in(gm.reference());
        std::cout << gm.reference() << ", " << region
                  << "\n----------------------------------------\n"
                  << in.subseq(reg.seq_name, reg.min, reg.max) << "\n\n";

        vargas::VCF vcf(gm.variants());
        vcf.set_region(region);
        std::cout << gm.variants() << ", "
                  << region
                  << "\n----------------------------------------\n";

        while (vcf.next()) {
            std::cout << vcf << '\n';
        }

        std::cout << std::endl;

    }

    if (fasta.size()) {
        vargas::ifasta in(fasta);
        std::cout << fasta << ", " << region
                  << "\n----------------------------------------\n"
                  << in.subseq(reg.seq_name, reg.min, reg.max)
                  << '\n' << std::endl;
    }

    if (vcf.size()) {
        vargas::VCF v(vcf);
        v.set_region(region);
        std::cout << vcf << ", " << region
                  << "\n----------------------------------------\n";
        while (v.next()) {
            std::cout << v << '\n';
        }
        std::cout << '\n' << std::endl;
    }

    return 0;

}

void main_help() {
    using std::cerr;
    using std::endl;
    cerr << endl
         << "---------------------- vargas, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------\n";
    cerr << "define          Define a set of graphs for use with sim and align.\n";
    cerr << "sim             Simulate reads from a set of graphs.\n";
    cerr << "align           Align reads to a set of graphs.\n";
    cerr << "export          Export graph to DOT format.\n";
    cerr << "convert         Convert a SAM file to a CSV file.\n";
    cerr << "query           Pull a region from a GDEF/VCF/FASTA file.\n";
    cerr << "test            Run unit tests.\n";
    cerr << "profile         Run profiles (debug).\n" << endl;
    cerr << "Compiled for architecture: ";

    std::string arch = "None: no SIMD instruction set specified!";

    #if SIMDPP_USE_SSE2
    arch = "SSE2";
    #endif
    #if SIMDPP_USE_SSE3
    arch = "SSE3";
    #endif
    #if SIMDPP_USE_SSSE3
    arch = "SSSE3";
    #endif
    #if SIMDPP_USE_SSE4_1
    arch = "SSE4.1";
    #endif
    #if SIMDPP_USE_AVX
    arch = "AVX";
    #endif
    #if SIMDPP_USE_AVX2
    arch = "AVX2";
    #endif
    #if SIMDPP_USE_FMA3
    arch = "FMA3";
    #endif
    #if SIMDPP_USE_FMA4
    arch = "FMA4";
    #endif
    #if SIMDPP_USE_XOP
    arch = "XOP";
    #endif
    #if SIMDPP_USE_AVX512F
    arch = "AVX512F";
    #endif
    #if SIMDPP_USE_NEON
    arch = "NEON";
    #endif

    cerr << arch << "\n" << std::endl;
}

void query_help(const cxxopts::Options &opts) {
    using std::cerr;
    using std::endl;
    cerr << opts.help() << "\n" << endl;
}

void export_help(const cxxopts::Options &opts) {
    using std::cerr;
    using std::endl;

    cerr << opts.help() << "\n" << endl;
}

void define_help(const cxxopts::Options &opts) {
    using std::cerr;
    using std::endl;

    cerr << opts.help() << "\n\n"
         << "Subgraphs are defined using the format \"label=N[%t]\",\n"
         << "where \'N\' is the number of samples / percentage of samples selected.\n"
         << "The samples are selected from the parent graph, scoped with \':\'.\n"
         << "The BASE graph is implied as the root for all labels. Example:\n"
         << "\ta=50;a:b=10%;~a:c=5\n"
         << "\'~\' indicates the complement graph. \'"
         << vargas::GraphManager::GDEF_BASEGRAPH
         << "\' refers the the whole graph.\n" << endl;
}

void profile_help(const cxxopts::Options &opts) {
    using std::cerr;
    using std::endl;
    cerr << opts.help() << "\n" << endl;
}

void sim_help(const cxxopts::Options &opts) {
    using std::cerr;
    using std::endl;

    cerr << opts.help() << "\n\n";
    cerr << "-n reads are produced for each -m, -i, -v, -b combination.\nIf set to \'*\', any value is accepted.\n"
         << endl;
}

void convert_help(const cxxopts::Options &opts) {
    using std::cerr;
    using std::endl;

    cerr << opts.help() << "\n\n";
    cerr << "Required column names:\n\tQNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL\n";
    cerr << "Prefix with \"RG:\" to obtain a value from the associated read group.\n" << endl;

}