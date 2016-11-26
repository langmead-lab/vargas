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

#define SIMDPP_ARCH_X86_SSE4_1 // SSE support
#define DOCTEST_CONFIG_IMPLEMENT // User controlled test execution

#include <iostream>
#include <thread>
#include <algorithm>
#include <omp.h>
#include "gdef.h"
#include "main.h"


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
        ("v,vcf", "<str> *VCF/BCF File.", cxxopts::value<std::string>(varfile))
        ("g,region", "<str> *Region of format \"CHR:MIN-MAX\". \"CHR:0-0\" for all.", cxxopts::value(region))
        ("s,subgraph", "<str> File or definitions of format \"[~]NAME=N[%t],...\".",
         cxxopts::value(subgraph_def))
        ("l,nodelen", "<N> Maximum node length.", cxxopts::value(node_len)->default_value("1000000"))
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
    if (!opts.count("v")) throw std::invalid_argument("VCF file required.");
    if (!opts.count("g")) throw std::invalid_argument("Region required.");

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

    if (threads > 0) omp_set_num_threads(threads);

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
    mut_split = split(mut, ','),
    indel_split = split(indel, ','),
    vnode_split = split(vnodes, ','),
    vbase_split = split(vbases, ',');

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
    std::cerr << chrono_duration(start_time) << " seconds." << std::endl;
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
        subdef_split = split(sim_src, gm.GDEF_DELIM);

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
    rg.date = current_date();
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
    std::cerr << chrono_duration(start_time)
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

    std::cerr << chrono_duration(start_time) << " seconds." << std::endl;

    return 0;
}

int align_main(int argc, char *argv[]) {
    std::string cl;
    {
        std::ostringstream ss;
        for (int i = 0; i < argc; ++i) ss << std::string(argv[i]) << " ";
        cl = ss.str();
    }

    // Load parameters
    // hisat similar params: match = 2, mismatch = 6, open = 5, extend = 3
    size_t match, mismatch, gopen, gext, threads, read_len, tolerance, chunk_size;
    std::string read_file, gdf_file, align_targets, out_file;
    bool align_targets_isfile = false;

    cxxopts::Options opts("vargas align", "Align reads to a graph.");
    try {
        opts.add_options()
        ("g,gdef", "<str> *Graph definition file.", cxxopts::value(gdf_file))
        ("r,reads", "<str> SAM reads file. (default: stdin)", cxxopts::value(read_file))
        ("a,align", "<str> Alignment targets/file of form \"RG:[ID][gd],target\"", cxxopts::value(align_targets))
        ("f,file", " -a specifies a file name.", cxxopts::value(align_targets_isfile))
        ("l,rlen", "<N> Maximum read length.", cxxopts::value(read_len)->default_value("50"))
        ("m,match", "<N> Match score.", cxxopts::value(match)->default_value("2"))
        ("n,mismatch", "<N> Mismatch penalty.", cxxopts::value(mismatch)->default_value("2"))
        ("o,gap_open", "<N> Gap opening penalty.", cxxopts::value(gopen)->default_value("3"))
        ("e,gap_extend", "<N> Gap extension penalty.", cxxopts::value(gext)->default_value("1"))
        ("c,tolerance", "<N> Correct if within readlen/N.",
         cxxopts::value(tolerance)->default_value(std::to_string(vargas::Aligner::default_tolerance())))
        ("h,chunk", "<N> Partition tasks into chunks with max size N.",
         cxxopts::value(chunk_size)->default_value("1024"))
        ("t,out", "<str> Output file. (default: stdout)", cxxopts::value(out_file))
        ("j,threads", "<N> Number of threads.", cxxopts::value(threads)->default_value("1"))
        ("h,help", "Display this message.");
        opts.parse(argc, argv);
    } catch (std::exception &e) { throw std::invalid_argument("Error parsing options: " + std::string(e.what())); }
    if (opts.count("h")) {
        align_help(opts);
        return 0;
    }
    if (!opts.count("g")) throw std::invalid_argument("Graph definition file required.");

    if (read_len * match > 255) {
        throw std::invalid_argument("Score matrix overflow with read length " + std::to_string(read_len) +
        " and match score " + std::to_string((int) match) + ".");
    }

    if (threads) omp_set_num_threads(threads);

    if (align_targets_isfile) {
        std::ifstream in(align_targets);
        if (!in.good()) throw std::invalid_argument("Invalid alignment targets file \"" + align_targets + "\".");
        std::stringstream ss;
        ss << in.rdbuf();
        align_targets = ss.str();
    }

    std::vector<std::string> alignment_pairs;
    if (align_targets.length() != 0) {
        std::replace(align_targets.begin(), align_targets.end(), '\n', ';');
        alignment_pairs = split(align_targets, ';');
    }

    std::cerr << "Match=" << match
              << " Mismatch=" << mismatch
              << " GapOpen=" << gopen
              << " GapExtend=" << gext
              << " MaxReadLen=" << read_len
              << " CorrectnessTol=" << tolerance
              << "\nLoading reads... " << std::flush;

    auto start_time = std::chrono::steady_clock::now();

    std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>> task_list;
    vargas::SAM::Header reads_hdr;
    size_t total = 0;
    {
        // Maps a read group ID to a vector of reads
        std::unordered_map<std::string, std::vector<vargas::SAM::Record>> read_groups;
        {
            vargas::isam reads(read_file);
            reads_hdr = reads.header();
            std::string read_group;
            vargas::SAM::Record rec;
            do {
                rec = reads.record();
                if (rec.seq.length() != read_len) {
                    throw std::invalid_argument("Expected read of length " +
                    std::to_string(read_len) + ", got " + std::to_string(rec.seq.length()));
                }
                if (!rec.aux.get("RG", read_group)) {
                    read_group = "NULL";
                    rec.aux.set("RG", "VA-UNGROUPED");
                    if (!reads_hdr.read_groups.count("VA-UNGROUPED")) {
                        reads_hdr.add(vargas::SAM::Header::ReadGroup("@RG\tID:VA-UNGROUPED"));
                    }
                }
                read_groups[read_group].push_back(rec);
            } while (reads.next());
        }

        if (alignment_pairs.size() == 0) {
            for (const auto &p : read_groups) {
                alignment_pairs.push_back("RG:ID:" + p.first + "\t" + vargas::GraphManager::GDEF_BASEGRAPH);
            }
        }

        // Maps target graph to read group ID's
        std::unordered_map<std::string, std::vector<std::string>> alignment_rg_map;
        {
            std::vector<std::string> pair;
            std::string tag, val, target_val;
            for (const std::string &p : alignment_pairs) {
                split(p, pair);
                if (pair.size() != 2)
                    throw std::invalid_argument("Malformed alignment pair \"" + p + "\".");
                if (pair[0].at(2) != ':')
                    throw std::invalid_argument("Expected source format Read_group_tag:value in \"" + pair[0] + "\".");
                if (pair[0].substr(0, 2) != "RG")
                    throw std::invalid_argument("Expected a read group tag \'RG:xx:\', got \"" + pair[0] + "\"");

                tag = pair[0].substr(3, 2);
                target_val = pair[0].substr(6);

                for (const auto &rg_pair : reads_hdr.read_groups) {
                    if (tag == "ID") val = rg_pair.second.id;
                    else if (rg_pair.second.aux.get(tag, val));
                    else continue;
                    if (val == target_val) alignment_rg_map[pair[1]].push_back(rg_pair.first);
                }

            }
        }

        std::cerr << chrono_duration(start_time) << " seconds." << std::endl;

        // graph label to vector of reads

        for (const auto &sub_rg_pair : alignment_rg_map) {
            for (const std::string &rgid : sub_rg_pair.second) {
                // If there is a header line that there are no reads associated with, skip
                if (read_groups.count(rgid) == 0) continue;

                const auto beg = std::begin(read_groups.at(rgid));
                const auto end = std::end(read_groups.at(rgid));
                const size_t nrecords = read_groups.at(rgid).size();
                const size_t n_chunks = (nrecords / chunk_size) + 1;
                total += read_groups.at(rgid).size();

                for (size_t i = 0; i < n_chunks; ++i) {
                    const auto safe_beg = beg + (i * chunk_size);
                    const auto safe_end = (i + 1) * chunk_size > nrecords ? end : safe_beg + chunk_size;
                    task_list.emplace_back(sub_rg_pair.first, std::vector<vargas::SAM::Record>(safe_beg, safe_end));
                }
            }
        }

        std::cerr << '\t' << read_groups.size() << " Read groups.\n"
                  << '\t' << alignment_rg_map.size() << " Subgraphs.\n"
                  << '\t' << task_list.size() << " Tasks.\n"
                  << '\t' << total << " Total alignments.\n";
    }

    std::cerr << "Loading graphs... " << std::flush;
    start_time = std::chrono::steady_clock::now();
    vargas::GraphManager gm(gdf_file);
    std::cerr << "(" << gm.base()->node_map()->size() << " nodes), ";
    std::cerr << chrono_duration(start_time) << " seconds." << std::endl;


    {
        vargas::SAM::Header::Program pg;
        pg.command_line = cl;
        pg.name = "vargas_align";
        pg.id = "VA";
        pg.version = __DATE__;
        std::replace_if(pg.version.begin(), pg.version.end(), isspace, ' '); // rm tabs
        reads_hdr.add(pg);
    }

    const size_t num_tasks = task_list.size();
    if (num_tasks < threads) {
        std::cerr << "Warning: Number of threads is greater than number of tasks. Try decreasing --chunk.\n";
    }

    std::cerr << "Aligning with " << threads << " thread(s)..." << std::endl;
    start_time = std::chrono::steady_clock::now();
    auto start_cpu = std::clock();


    vargas::osam aligns_out(out_file, reads_hdr);

    #pragma omp parallel for
    for (size_t l = 0; l < num_tasks; ++l) {
        const size_t num_reads = task_list.at(l).second.size();
        std::vector<std::string> read_seqs(num_reads);
        std::vector<size_t> targets(num_reads);
        for (size_t i = 0; i < num_reads; ++i) {
            const auto &r = task_list.at(l).second.at(i);
            read_seqs[i] = r.seq;
            targets[i] = r.pos + r.seq.length() - 1;
        }
        vargas::Aligner aligner(gm.node_len(), read_len, match, mismatch, gopen, gext);
        aligner.set_correctness_tolerance(tolerance);
        task_list.at(l).first;
        auto subgraph = gm.make_subgraph(task_list.at(l).first);
        auto aligns = aligner.align(read_seqs, targets, subgraph->begin(), subgraph->end());
        for (size_t j = 0; j < task_list.at(l).second.size(); ++j) {
            vargas::SAM::Record &rec = task_list.at(l).second.at(j);
            rec.ref_name = task_list.at(l).first;
            rec.aux.set(ALIGN_SAM_MAX_POS_TAG, (int) aligns.max_pos[j]);
            rec.aux.set(ALIGN_SAM_MAX_SCORE_TAG, aligns.max_score[j]);
            rec.aux.set(ALIGN_SAM_MAX_COUNT_TAG, aligns.max_count[j]);
            rec.aux.set(ALIGN_SAM_SUB_POS_TAG, (int) aligns.sub_pos[j]);
            rec.aux.set(ALIGN_SAM_SUB_SCORE_TAG, aligns.sub_score[j]);
            rec.aux.set(ALIGN_SAM_SUB_COUNT_TAG, aligns.sub_count[j]);
            rec.aux.set(ALIGN_SAM_COR_FLAG_TAG, aligns.correctness_flag[j]);
        }
        gm.destroy(task_list.at(l).first);

        #pragma omp critical(aligns_out)
        {
            for (size_t j = 0; j < task_list.at(l).second.size(); ++j) {
                aligns_out.add_record(task_list.at(l).second.at(j));
            }
        }
    }

    auto end_time = std::chrono::steady_clock::now();
    auto cput = (std::clock() - start_cpu) / (double) CLOCKS_PER_SEC;
    std::cerr << chrono_duration(start_time, end_time) << " seconds, "
              << cput << " CPU seconds, "
              << cput / total << " CPU s/alignment.\n" << std::endl;

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
    std::vector<std::string> fmt_split = split(format, ',');

    std::unordered_set<std::string> warned;

    vargas::isam input(sam_file);

    std::string buff, val;
    do {
        buff = "";
        for (auto &tag : fmt_split) {
            val = "*";
            if (!input.record().get(input.header(), tag, val) && warned.count(tag) == 0) {
                std::cerr << "Warning: Tag \"" << tag << "\" not present." << std::endl;
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
    if (!opts.count("v")) throw std::invalid_argument("VCF file required.");
    if (!opts.count("g")) throw std::invalid_argument("Region specifier required.");


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

    {
        std::vector<std::string> reads;

        for (size_t i = 0; i < nreads; ++i) {
            std::ostringstream rd;
            for (size_t r = 0; r < read_len; ++r) rd << rand_base();
            reads.push_back(rd.str());
        }

        std::vector<std::string> split_str;

        vargas::Aligner a(g.max_node_len(), 50);

        std::cerr << nreads << " read alignment:\n";


        {
            start = std::clock();
            vargas::Graph g2(g, g.subset(ingroup));
            std::cerr << "\tDerived Graph (" << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s)\n\t";
            start = std::clock();
            vargas::Aligner::Results aligns = a.align(reads, g2.begin(), g2.end());
            std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
        }

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

void align_help(const cxxopts::Options &opts) {
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

TEST_SUITE("System");

TEST_CASE ("Coordinate System Matches") {
    srand(1);
    vargas::Graph::Node::_newID = 0;
    using std::endl;
    std::string tmpfa = "tmp_tc.fa";
    {
        std::ofstream fao(tmpfa);
        fao
        << ">x" << endl
        << "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGTTCCTGGTGCTATGTGTAACTAGTAATGG" << endl
        << "TAATGGATATGTTGGGCTTTTTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAA" << endl
        << "TGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTTCAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGA" << endl
        << "ACAAGTTAGTTAATCTCTCTGAACTTCAGTTTAATTATCTCTAATATGGAGATGATACTACTGACAGCAGAGGTTTGCTG" << endl
        << "TGAAGATTAAATTAGGTGATGCTTGTAAAGCTCAGGGAATAGTGCCTGGCATAGAGGAAAGCCTCTGACAACTGGTAGTT" << endl
        << "ACTGTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCCTCCTTGAGGCTCT" << endl
        << "TCTGGCTTTTCATTGTCAACACAGTCAACGCTCAATACAAGGGACATTAGGATTGGCAGTAGCTCAGAGATCTCTCTGCT" << endl
        << ">y" << endl
        << "GGAGCCAGACAAATCTGGGTTCAAATCCTGGAGCCAGACAAATCTGGGTTCAAATCCTGGAGCCAGACAAATCTGGGTTC" << endl;
    }
    std::string tmpvcf = "tmp_tc.vcf";

    {
        std::ofstream vcfo(tmpvcf);
        vcfo
        << "##fileformat=VCFv4.1" << endl
        << "##phasing=true" << endl
        << "##contig=<ID=x>" << endl
        << "##contig=<ID=y>" << endl
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Freq\">" << endl
        << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele count\">" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Num samples at site\">" << endl
        << "##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Num alt alleles\">" << endl
        << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"Length of each alt\">" << endl
        << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"type of variant\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2" << endl
        << "x\t9\t.\tG\tA,CC,T\t99\t.\tAF=0.01,0.6,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t0|1\t2|3" << endl
        << "x\t10\t.\tC\t<CN7>,<CN0>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "x\t14\t.\tG\t<DUP>,<BLAH>\t99\t.\tAF=0.01,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t1|1" << endl
        << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
    }

    vargas::GraphFactory gb(tmpfa);
    gb.open_vcf(tmpvcf);
    gb.node_len(5);
    gb.set_region("x:0-50");
    vargas::Graph g = gb.build();

    vargas::Sim::Profile prof;
    prof.len = 5;
    vargas::Sim sim(g, prof);

    vargas::Aligner aligner(g.max_node_len(), 5);
    auto reads = sim.get_batch(aligner.read_capacity());

    std::vector<std::string> seqs;
    std::vector<size_t> targets;
    for (auto &r : reads) {
        seqs.push_back(r.seq);
        targets.push_back(r.pos + r.seq.length() - 1);
    }

    auto results = aligner.align(seqs, targets, g.begin(), g.end());

    for (auto i : results.correctness_flag) CHECK ((int) i == 1);

    remove(tmpfa.c_str());
    remove(tmpvcf.c_str());
    remove((tmpfa + ".fai").c_str());
}
TEST_CASE ("Correctness flag") {
    srand(1);
    vargas::Graph::Node::_newID = 0;
    using std::endl;
    std::string tmpfa = "tmp_tc.fa";
    {
        std::ofstream fao(tmpfa);
        fao
        << ">x" << endl
        << "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGTTCCTGGTGCTATGTGTAACTAGTAATGG" << endl
        << "TAATGGATATGTTGGGCTTTTTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAA" << endl
        << "TGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTTCAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGA" << endl
        << "ACAAGTTAGTTAATCTCTCTGAACTTCAGTTTAATTATCTCTAATATGGAGATGATACTACTGACAGCAGAGGTTTGCTG" << endl
        << "TGAAGATTAAATTAGGTGATGCTTGTAAAGCTCAGGGAATAGTGCCTGGCATAGAGGAAAGCCTCTGACAACTGGTAGTT" << endl
        << "ACTGTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCCTCCTTGAGGCTCT" << endl
        << "TCTGGCTTTTCATTGTCAACACAGTCAACGCTCAATACAAGGGACATTAGGATTGGCAGTAGCTCAGAGATCTCTCTGCT" << endl
        << ">y" << endl
        << "GGAGCCAGACAAATCTGGGTTCAAATCCTGGAGCCAGACAAATCTGGGTTCAAATCCTGGAGCCAGACAAATCTGGGTTC" << endl;
    }
    std::string tmpvcf = "tmp_tc.vcf";

    {
        std::ofstream vcfo(tmpvcf);
        vcfo
        << "##fileformat=VCFv4.1" << endl
        << "##phasing=true" << endl
        << "##contig=<ID=x>" << endl
        << "##contig=<ID=y>" << endl
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Freq\">" << endl
        << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele count\">" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Num samples at site\">" << endl
        << "##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Num alt alleles\">" << endl
        << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"Length of each alt\">" << endl
        << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"type of variant\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2" << endl
        << "x\t9\t.\tG\tA,CC,T\t99\t.\tAF=0.01,0.6,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t0|1\t2|3" << endl
        << "x\t10\t.\tC\t<CN7>,<CN0>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "x\t14\t.\tG\t<DUP>,<BLAH>\t99\t.\tAF=0.01,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t1|1" << endl
        << "x\t20\t.\tTTC\t<CN3>,<CN2>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
    }

    std::string reads_file("tmp_rd.sam");
    {
        std::ofstream ro(reads_file);
        ro << "@HD\tVN:1.0\n*\t4\t*\t14\t255\t*\t*\t0\t0\tGAAATT\t*\n*\t4\t*\t17\t255\t*\t*\t0\t0\tATTTTC\t*";
    }

    vargas::GraphFactory gb(tmpfa);
    gb.open_vcf(tmpvcf);
    gb.set_region("x:0-100");
    vargas::Graph g = gb.build();

    vargas::Aligner aligner(g.max_node_len(), 6);
    vargas::isam reads(reads_file);

    std::vector<vargas::SAM::Record> records;
    std::vector<std::string> read_seq;
    std::vector<size_t> targets;
    do {
        records.push_back(reads.record());
        read_seq.push_back(reads.record().seq);
        targets.push_back(reads.record().pos + read_seq.back().length() - 1);
    } while (reads.next());

    auto res = aligner.align(read_seq, targets, g.begin(), g.end());

    for (size_t i = 0; i < records.size(); ++i) {
        records[i].aux.set(ALIGN_SAM_MAX_POS_TAG, (int) res.max_pos[i]);
        records[i].aux.set(ALIGN_SAM_MAX_SCORE_TAG, (int) res.max_score[i]);
        records[i].aux.set(ALIGN_SAM_MAX_COUNT_TAG, (int) res.max_count[i]);

        records[i].aux.set(ALIGN_SAM_SUB_POS_TAG, (int) res.sub_pos[i]);
        records[i].aux.set(ALIGN_SAM_SUB_SCORE_TAG, (int) res.sub_score[i]);
        records[i].aux.set(ALIGN_SAM_SUB_COUNT_TAG, (int) res.sub_count[i]);

        records[i].aux.set(ALIGN_SAM_COR_FLAG_TAG, (int) res.correctness_flag[i]);
    }

    vargas::osam align_out("tmp_aout.sam", reads.header());
    for (auto &r : records) align_out.add_record(r);

    remove(tmpfa.c_str());
    remove((tmpfa + ".fai").c_str());
    remove(tmpvcf.c_str());
    remove(reads_file.c_str());
    remove("tmp_aout.sam");
}

TEST_SUITE_END();