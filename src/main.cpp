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
#include "getopt_pp.h"


int main(const int argc, const char *argv[]) {

    srand(time(NULL)); // Rand used in profiles and sim

    try {
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
                return align_main(argc, argv);
            }
            else if (!strcmp(argv[1], "export")) {
                return export_main(argc, argv);
            } else if (!strcmp(argv[1], "convert")) {
                return sam2csv(argc, argv);
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

    GetOpt::GetOpt_pp args(argc, argv);
    if (args >> GetOpt::OptionPresent('h', "help")) {
        main_help();
        return 0;
    }

    std::cerr << "Define a valid mode of operation." << std::endl;
    main_help();
    return 1;

}

int define_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        define_help();
        return 0;
    }

    std::string
        fasta_file = "",
        varfile = "",
        region = "",
        subgraph_def = "",
        out_file = "",
        dot_file = "",
        sample_filter = "";

    bool invert_filter = false, build_base = false;

    int node_len = 1000000;

    args >> GetOpt::Option('g', "region", region)
         >> GetOpt::Option('l', "nodelen", node_len)
         >> GetOpt::Option('s', "subgraph", subgraph_def)
         >> GetOpt::Option('t', "out", out_file)
         >> GetOpt::Option('d', "dot", dot_file)
         >> GetOpt::Option('p', "filter", sample_filter)
         >> GetOpt::OptionPresent('x', "invert", invert_filter)
         >> GetOpt::OptionPresent('b', "base", build_base);

    if (!(args >> GetOpt::Option('f', "fasta", fasta_file) >> GetOpt::Option('v', "vcf", varfile))) {
        define_help();
        throw std::invalid_argument("No FASTA or VCF file provided!");
    }

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
    gm.write(fasta_file, varfile, region, subgraph_str, node_len, out_file, build_base);
    if (dot_file.length() > 0) gm.to_DOT(dot_file, "subgraphs");

    std::cerr << gm.size() << " subgraph definitions generated." << std::endl;

    return 0;
}

int sim_main(const int argc, const char *argv[]) {

    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        sim_help();
        return 0;
    }

    vargas::SAM::Header sam_hdr;

    {
        vargas::SAM::Header::Program pg;
        std::ostringstream ss;
        for (int i = 0; i < argc; ++i) ss << std::string(argv[i]) << " ";
        pg.command_line = ss.str();
        pg.name = "vargas_sim";
        pg.id = "VS";
        pg.version = __DATE__;
        std::replace_if(pg.version.begin(), pg.version.end(), isspace, ' '); // rm tabs
        sam_hdr.add(pg);
    }

    // Load parameters
    int
        read_len = 50,
        num_reads = 1000,
        threads = 1;

    std::string
        mut = "0",
        indel = "0",
        vnodes = "*",
        vbases = "*",
        gdf_file = "",
        out_file = "",
        sim_src = "";

    bool use_rate = false, sim_src_isfile = false;

    args >> GetOpt::Option('g', "gdef", gdf_file)
         >> GetOpt::Option('s', "sub", sim_src)
         >> GetOpt::OptionPresent('f', "file", sim_src_isfile)
         >> GetOpt::Option('d', "vnodes", vnodes)
         >> GetOpt::Option('b', "vbases", vbases)
         >> GetOpt::Option('l', "rlen", read_len)
         >> GetOpt::Option('n', "numreads", num_reads)
         >> GetOpt::Option('m', "mut", mut)
         >> GetOpt::Option('i', "indel", indel)
         >> GetOpt::Option('t', "out", out_file)
         >> GetOpt::Option('j', "threads", threads)
         >> GetOpt::OptionPresent('a', "rate", use_rate);

    omp_set_num_threads(threads);


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
        subdef_split = gm.labels();
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

    std::cerr << std::endl << chrono_duration(start_time) << " seconds." << std::endl;

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
    std::string read_file, gdf_file, align_targets, out_file;
    bool align_targets_isfile = false;

    args >> GetOpt::Option('m', "match", match)
         >> GetOpt::Option('n', "mismatch", mismatch)
         >> GetOpt::Option('o', "gap_open", gopen)
         >> GetOpt::Option('e', "gap_extend", gext)
         >> GetOpt::Option('j', "threads", threads)
         >> GetOpt::Option('l', "rlen", read_len)
         >> GetOpt::Option('t', "out", out_file)
         >> GetOpt::Option('a', "align", align_targets)
         >> GetOpt::OptionPresent('f', "file", align_targets_isfile);

    if (!(args >> GetOpt::Option('g', "gdef", gdf_file))) {
        align_help();
        throw std::invalid_argument("No GDEF file provided.");
    }

    size_t tolerance = read_len / 2;
    args >> GetOpt::Option('c', "tolerance", tolerance);

    if (threads) omp_set_num_threads(threads);

    if (align_targets_isfile) {
        std::ifstream in(align_targets);
        if (!in.good()) throw std::invalid_argument("Invalid alignment targets file \"" + align_targets + "\".");
        std::stringstream ss;
        ss << in.rdbuf();
        align_targets = ss.str();
    }

    std::replace(align_targets.begin(), align_targets.end(), '\n', ';');
    const std::vector<std::string> alignment_pairs = split(align_targets, ';');

    std::cerr << "\nLoading reads... " << std::flush;
    auto start_time = std::chrono::steady_clock::now();

    std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>> task_list;
    vargas::SAM::Header reads_hdr;
    size_t total = 0;
    {
        // Maps a read group ID to a vector of reads
        std::unordered_map<std::string, std::vector<vargas::SAM::Record>> alignment_reads;
        {
            vargas::isam reads(read_file);
            reads_hdr = reads.header();
            std::string read_group;
            do {
                if (!reads.record().aux.get("RG", read_group)) read_group = "NULL";
                alignment_reads[read_group].push_back(reads.record());
            } while (reads.next());
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
                    if (rg_pair.second.aux.get(tag, val)) {
                        if (val == target_val) alignment_rg_map[pair[1]].push_back(rg_pair.first);
                    }
                }

            }
        }

        std::cerr << chrono_duration(start_time) << " seconds." << std::endl;
        std::cerr << "\tSubgraph\t# Reads\tRG ID" << std::endl;

        // graph label to vector of reads
        for (const auto &sub_rg_pair : alignment_rg_map) {
            for (const std::string &rgid : sub_rg_pair.second) {
                if (alignment_reads.count(rgid)) {
                    // If there is a header line that there are no reads associated with, skip
                    task_list.push_back(std::pair<std::string,
                                                  std::vector<vargas::SAM::Record>>(sub_rg_pair.first,
                                                                                                 alignment_reads.at(rgid)));
                    std::cerr << '\t' << sub_rg_pair.first << '\t' << alignment_reads.at(rgid).size()
                              << '\t' << rgid << '\n';
                    total += alignment_reads.at(rgid).size();
                }
            }

        }

        std::cerr << '\t' << alignment_reads.size() << " Read groups.\n"
                  << '\t' << alignment_rg_map.size() << " Subgraphs.\n"
                  << '\t' << task_list.size() << " Tasks.\n"
                  << '\t' << total << " Total alignments.\n";
    }

    std::cerr << "Loading graphs... " << std::flush;
    start_time = std::chrono::steady_clock::now();
    vargas::GraphManager gm(gdf_file);
    std::cerr << chrono_duration(start_time) << " seconds." << std::endl;


    {
        vargas::SAM::Header::Program pg;
        std::ostringstream ss;
        for (int i = 0; i < argc; ++i) ss << std::string(argv[i]) << " ";
        pg.command_line = ss.str();
        pg.name = "vargas_align";
        pg.id = "VA";
        pg.version = __DATE__;
        std::replace_if(pg.version.begin(), pg.version.end(), isspace, ' '); // rm tabs
        reads_hdr.add(pg);
    }


    std::cerr << "Aligning with " << threads << " thread(s)..." << std::endl;
    start_time = std::chrono::steady_clock::now();
    auto start_cpu = std::clock();

    const size_t num_tasks = task_list.size();

    vargas::osam aligns_out(out_file, reads_hdr);

    #pragma omp parallel for
    for (size_t l = 0; l < num_tasks; ++l) {
        const size_t num_reads = task_list.at(l).second.size();
        std::vector<std::string> read_seqs(num_reads);
        std::vector<uint32_t> targets(num_reads);
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

int sam2csv(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        sam2csv_help();
        return 0;
    }

    std::string sam_file = "", format = "";

    args >> GetOpt::Option('s', "sam", sam_file);
    if (!(args >> GetOpt::Option('f', "format", format))) {
        sam2csv_help();
        throw std::invalid_argument("Not format provided!");
    }

    auto start_time = std::chrono::steady_clock::now();

    format.erase(std::remove(format.begin(), format.end(), ' '), format.end());
    std::vector<std::string> fmt_split = split(format, ',');

    std::unordered_set<std::string> warned;

    vargas::isam input(sam_file);

    std::string buff, val;
    do {
        buff = "";
        for (auto &tag : fmt_split) {
            buff += ",";
            if (tag == "POS") {
                buff += std::to_string(input.record().pos);
            }
            else if (tag == "QNAME") {
                buff += input.record().query_name;
            }
            else if (tag == "RNEXT") {
                buff += input.record().ref_next;
            }
            else if (tag == "RNAME") {
                buff += input.record().ref_name;
            }
            else if (tag == "SEQ") {
                buff += input.record().seq;
            }
            else if (tag == "CIGAR") {
                buff += input.record().cigar;
            }
            else if (tag == "FLAG") {
                buff += std::to_string(input.record().flag.encode());
            }
            else if (tag == "PNEXT") {
                buff += std::to_string(input.record().pos_next);
            }
            else if (tag == "MAPQ") {
                buff += std::to_string(input.record().mapq);
            }
            else if (tag == "TLEN") {
                buff += std::to_string(input.record().tlen);
            }
            else if (tag == "QUAL") {
                buff += input.record().qual;
            }
            else if (tag.substr(0, 3) == "RG:") {
                val = "*";
                if (!input.record().read_group(input.header(), tag.substr(3), val) && warned.count(tag) == 0) {
                    std::cerr << "Warning: \"" << tag << "\" tag not present in read group." << std::endl;
                    warned.insert(tag);
                }
                buff += val;
            }
            else {
                val = "*";
                if (!input.record().aux.get(tag, val) && warned.count(tag) == 0) {
                    std::cerr << "Warning: \"" << tag << "\" tag not present." << std::endl;
                    warned.insert(tag);
                }
                buff += val;
            }
        }
        std::cout << buff.substr(1) << '\n'; // Crop leading comma
    } while (input.next());

    std::cerr << std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::steady_clock::now() - start_time).count()
              << " seconds." << std::endl;

    return 0;
}

int profile(const int argc, const char *argv[]) {
    std::string bcf, fasta;
    std::string region = "22:25,000,000-25,500,000";
    std::string read;
    int ingroup = 100;

    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        profile_help();
        return 0;
    }

    if (!(args >> GetOpt::Option('f', "fasta", fasta)
               >> GetOpt::Option('v', "var", bcf)
               >> GetOpt::Option('g', "region", region))) {
        throw std::invalid_argument("FASTA & Variant File with region required.");
    }

    args >> GetOpt::Option('i', "ingroup", ingroup)
         >> GetOpt::Option('s', "string", read);

    if (!file_exists(fasta) || !file_exists(bcf)) {
        throw std::invalid_argument("File does not exist.");
    }

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
        std::cerr << "Filter constructor:\n\t";
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

        vargas::Aligner a(g.max_node_len(), 50);

        std::cerr << SIMDPP_FAST_INT8_SIZE << " read alignment:\n";


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

int export_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        export_help();
        return 0;
    }

    // Load parameters
    std::string subgraph = "BASE";
    std::string file = "";

    args >> GetOpt::Option('g', "graph", subgraph)
         >> GetOpt::Option('t', "out", file);

    if (file.length() == 0) {
        export_help();
        throw std::invalid_argument("No output file specified.");
    }

    vargas::GraphManager gm(std::cin);
    auto g = gm.make_subgraph(subgraph);
    g->to_DOT(file, "g");
    return 0;
}

int query_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        query_help();
        return 0;
    }

    std::string region, file;

    if (!(args >> GetOpt::Option('g', "region", region))) {
        query_help();
        throw std::invalid_argument("Region argument required!");
    }

    auto reg = vargas::parse_region(region);

    if (args >> GetOpt::Option('d', "gdef", file)) {
        vargas::GraphManager gm(file);
        std::string subgraph = gm.GDEF_BASEGRAPH, out = "";

        args >> GetOpt::Option('s', "subgraph", subgraph)
             >> GetOpt::Option('t', "out", out);

        auto sg = gm.make_subgraph(subgraph)->subgraph(reg.min, reg.max);
        if (out.length() == 0) std::cout << sg.to_DOT();
        else {
            std::ofstream o(out);
            o << sg.to_DOT();
        }

        std::cout << std::endl;

        vargas::ifasta in(gm.reference());
        std::cout << gm.reference() << ", " << region
                  << " \n----------------------------------------\n"
                  << in.subseq(reg.seq_name, reg.min, reg.max) << "\n\n";

        vargas::VCF vcf(gm.variants());
        vcf.set_region(region);
        std::cout << gm.variants() << ", "
                  << region << "\n----------------------------------------\n";

        while (vcf.next()) {
            std::cout << vcf << '\n';
        }

        std::cout << std::endl;

    }

    if (args >> GetOpt::Option('f', "fasta", file)) {
        vargas::ifasta in(file);
        std::cout << file << ", " << region
                  << "\n----------------------------------------\n"
                  << in.subseq(reg.seq_name, reg.min, reg.max)
                  << '\n' << std::endl;
    }

    if (args >> GetOpt::Option('v', "vcf", file)) {
        vargas::VCF vcf(file);
        vcf.set_region(region);
        std::cout << file << ", " << region
                  << "\n----------------------------------------\n";
        while (vcf.next()) {
            std::cout << vcf << '\n';
        }
        std::cout << '\n' << std::endl;
    }

    return 0;

}

void query_help() {
    using std::cerr;
    using std::endl;
    cerr << endl
         << "-------------------- vargas query, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-g\t--region        *<string> Region to export, format CHR:MIN-MAX\n";
    cerr << "-f\t--fasta         <string> Get a subsequence.\n";
    cerr << "-v\t--vcf           <string> VCF or BCF file.\n";
    cerr << "-d\t--gdef          <string> Query a graph, export DOT format.\n";
    cerr << "-t\t--out           <string> DOT output file for -g.\n";
    cerr << "-s\t--subgraph      <string> Subgraph of GDEF to query, default is the whole graph.\n" << endl;
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

void export_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
         << "-------------------- vargas export, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-t\t--out           Output filename.\n";
    cerr << "-g\t--graph         Subgraph to export, default BASE\n\n";
    cerr << "vargas export -g \"IN\" -t out.dot < graph_definition.gdef\n" << std::endl;

}

void define_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
         << "-------------------- vargas define, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-f\t--fasta         *<string> Reference filename.\n";
    cerr << "-v\t--vcf           *<string> VCF or BCF file.\n";
    cerr << "-g\t--region        *<string> Region of graph, format CHR:MIN-MAX.\n";
    cerr << "-s\t--subgraph      <string> Subgraph definition or filename.\n";
    cerr << "-p\t--filter        <string> Filename of sample filter.\n";
    cerr << "-x\t--invert        Invert sample filter.\n";
    cerr << "-l\t--nodelen       <int> Max node length, default 1,000,000\n";
    cerr << "-t\t--out           <string> Output file, default stdout.\n";
    cerr << "-d\t--dot           <string> Output DOT graph of subgraph hierarchy to file.\n\n";
    //cerr << "-b\t--base          Build base graph, used for debugging.\n" << endl;

    cerr << "Subgraphs are defined using the format \"label=N[%t]\".\n"
         << "\t where \'N\' is the number of samples / percentage of samples selected.\n"
         << "The samples are selected from the parent graph, scoped with \':\'.\n"
         << "The BASE graph is implied as the root for all labels. Example:\n"
         << "\ta=50;a:b=10%;~a:c=5\n"
         << "\'~\' indicates the complement graph. \'BASE\' refers the the whole graph.\n" << endl;
}

void profile_help() {
    using std::cerr;
    using std::endl;
    cerr << endl
         << "---------------------- vargas profile, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------"
         << endl;
    cerr << "-f\t--fasta         *<string> Reference filename." << endl;
    cerr << "-v\t--var           *<string> VCF/BCF filename." << endl;
    cerr << "-g\t--region        *<string> Region of graph, format CHR:MIN-MAX." << endl;
    cerr << "-i\t--ingroup       <int> Percent of genotypes to include in alignment." << endl;
    cerr << "-s\t--string        <string,string..> Include reads in alignment. Rest will be random." << endl << endl;
}

void align_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
         << "------------------- vargas align, " << __DATE__ << ". rgaddip1@jhu.edu -------------------\n";;
    cerr << "-g\t--gdef          *<string> Graph definition file.\n";
    cerr << "-r\t--reads         *<string> SAM file to align. Default stdin.\n";
    cerr << "-a\t--align         *<string:string> Alignment targets, origin graph : target graph.\n";
    cerr << "-t\t--out           *<string> Alignment output file, default stdout.\n";
    cerr << "-f\t--file          -a specifies a file name.\n";
    cerr << "-l\t--rlen          <int> Max read length. Default 50.\n";
    cerr << "-m\t--match         <int> Match score, default 2.\n";
    cerr << "-n\t--mismatch      <int> Mismatch penalty, default 2.\n";
    cerr << "-o\t--gap_open      <int> Gap opening penalty, default 3.\n";
    cerr << "-e\t--gap_extend    <int> Gap extend penalty, default 1.\n";
    cerr << "-c\t--tolerance     <int> Count an alignment as correct if within -c, default read_len/2\n";
    cerr << "-j\t--threads       <int> Number of threads. 0 for maximum hardware concurrency.\n" << endl;
}

void sim_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
         << "-------------------- vargas sim, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-g\t--gdef          *<string> Graph definition file. Default stdin.\n";
    cerr << "-s\t--sub           <string(;string)*> list of graphs to simulate from. Default all.\n";
    cerr << "-f\t--file          -s specifies a file name.\n";
    cerr << "-t\t--out           <string> Output file. Default stdout.\n";
    cerr << "-n\t--numreads      <int> Number of reads to simulate from each profile, default 1000.\n";
    cerr << "-m\t--muterr        <int/float, int/float...> Read mutation error. Default 0.\n";
    cerr << "-i\t--indelerr      <int/float, int/float...> Read indel error. Default 0.\n";
    cerr << "-d\t--vnodes        <int, int...> Number of variant nodes, default any (*).\n";
    cerr << "-b\t--vbases        <int, int...> Number of variant bases, default any (*).\n";
    cerr << "-l\t--rlen          <int> Read length, default 50.\n";
    cerr << "-a\t--rate          Interpret -m, -i as rates, instead of exact number of errors.\n";
    cerr << "-j\t--threads       <int> Number of threads. 0 for maximum hardware concurrency.\n" << endl;

    cerr << "-n reads are produced for each -m, -i, -v, -b combination. If set to \'*\', any value is accepted."
         << endl << endl;
}

void sam2csv_help() {
    using std::cerr;
    using std::endl;

    cerr << endl << "-------------------- vargas convert, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-s\t--sam          <string> SAM input file. Default stdin.\n";
    cerr << "-f\t--format       *<string,string...> Specify tags per column. Case sensitive.\n";
    cerr << "\nOutput printed to stdout.\n";
    cerr << "Required column names:\n\tQNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL\n";
    cerr << "Prefix with \"RG:\" to obtain a value from the associated read group.\n" << endl;

}