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

#define TIME_ALIGNMENT 0 // Disable output and report time to align

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
            else if (!strcmp(argv[1], "split")) {
                return split_main(argc, argv);
            }
            else if (!strcmp(argv[1], "convert")) {
                return sam2csv(argc, argv);
            }
            else if (!strcmp(argv[1], "merge")) {
                return merge_main(argc, argv);
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
        dot_file = "";

    int node_len = 1000000;

    args >> GetOpt::Option('f', "fasta", fasta_file)
         >> GetOpt::Option('g', "region", region)
         >> GetOpt::Option('l', "nodelen", node_len)
         >> GetOpt::Option('s', "subgraph", subgraph_def)
         >> GetOpt::Option('t', "out", out_file)
         >> GetOpt::Option('d', "dot", dot_file)
         >> GetOpt::Option('v', "vcf", varfile);


    std::string subgraph_str = "";

    if (subgraph_def.length() == 0) {
        std::string line;
        while (std::getline(std::cin, line)) subgraph_str += line + "\n";
    } else {
        std::ifstream in(subgraph_def);
        if (!in.good()) throw std::invalid_argument("Error opening file \"" + subgraph_def + "\".");
        std::stringstream buff;
        buff << in.rdbuf();
        subgraph_str = buff.str();
    }

    Vargas::GraphManager gm;
    gm.write(fasta_file, varfile, region, subgraph_str, node_len, out_file);
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

    Vargas::SAM::Header sam_hdr;

    {
        Vargas::SAM::Header::Program pg;
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

    if (threads) omp_set_num_threads(threads);


    Vargas::GraphManager gm;

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
    std::replace(sim_src.begin(), sim_src.end(), '\n', gm.GDEF_DELIM);
    sim_src.erase(std::remove_if(sim_src.begin(), sim_src.end(), isspace), sim_src.end());
    const std::vector<std::string> subdef_split = split(sim_src, gm.GDEF_DELIM);

    std::cerr << "Loading base graph... " << std::flush;
    auto start_time = std::chrono::steady_clock::now();
    gm.open(gdf_file);
    std::cerr << chrono_duration(start_time) << " seconds." << std::endl;

    // validate graph labels
    for (const std::string &l : subdef_split) {
        gm.filter(l); // Will throw if l does not exist
    }

    std::cerr << "Building profiles... " << std::flush;
    start_time = std::chrono::steady_clock::now();

    // Map a graph label to a vector of read group ID's/profiles
    std::unordered_map<std::string, // Graph label
                       std::vector<std::pair<std::string, // Read Group ID
                                             Vargas::Sim::Profile>>> // Sim profile
        queue;

    int rg_id = 0;
    Vargas::SAM::Header::ReadGroup rg;
    rg.seq_center = "vargas_sim";
    rg.date = current_date();
    rg.aux.set(SIM_SAM_REF_TAG, gm.reference());
    rg.aux.set(SIM_SAM_VCF_TAG, gm.variants());

    Vargas::Sim::Profile prof;
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
    std::cerr << sam_hdr.read_groups.size() << " read groups over " << subdef_split.size() << " subgraphs. "
              << chrono_duration(start_time)
              << " seconds." << std::endl;

    Vargas::osam out(out_file, sam_hdr);
    if (!out.good()) throw std::invalid_argument("Error opening output file \"" + out_file + "\"");

    std::cerr << "Simulating... " << std::flush;

    start_time = std::chrono::steady_clock::now();

    std::vector<std::pair<std::string, // Graph label
                          std::pair<std::string, // RG ID
                                    Vargas::Sim::Profile>>> // sim prof
        task_list;

    for (size_t k = 0; k < subdef_split.size(); ++k) {
        for (size_t i = 0; i < queue.at(subdef_split[k]).size(); ++i) {
            auto &p = queue.at(subdef_split[k]).at(i);
            task_list.push_back(std::pair<std::string, std::pair<std::string, Vargas::Sim::Profile>>
                                    (subdef_split[k], std::pair<std::string, Vargas::Sim::Profile>(p.first, p.second)));
        }
    }

    const size_t num_tasks = task_list.size();
    #pragma omp parallel for
    for (size_t n = 0; n < num_tasks; ++n) {
        const std::string label = task_list.at(n).first;
        const auto subgraph_ptr = gm.make_subgraph(label);
        Vargas::Sim sim(*subgraph_ptr, task_list.at(n).second.second);
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
         >> GetOpt::Option('r', "reads", read_file)
         >> GetOpt::Option('g', "gdef", gdf_file)
         >> GetOpt::Option('j', "threads", threads)
         >> GetOpt::Option('l', "rlen", read_len)
         >> GetOpt::Option('t', "out", out_file)
         >> GetOpt::Option('a', "align", align_targets)
         >> GetOpt::OptionPresent('f', "file", align_targets_isfile);

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

    std::cerr << "Loading reads... " << std::endl;
    auto start_time = std::chrono::steady_clock::now();

    std::vector<std::pair<std::string, std::vector<Vargas::SAM::Record>>> task_list;
    Vargas::SAM::Header reads_hdr;
    size_t total = 0;
    {
        // Maps a read group ID to a vector of reads
        std::unordered_map<std::string, std::vector<Vargas::SAM::Record>> alignment_reads;
        {
            Vargas::isam reads(read_file);
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
                split(p, '\t', pair);
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

        std::cerr << '\t' << "Graph\t# Reads" << std::endl;

        // graph label to vector of reads
        for (const auto &sub_rg_pair : alignment_rg_map) {
            for (const std::string &rgid : sub_rg_pair.second) {
                if (alignment_reads.count(rgid)) {
                    // If there is a header line that there are no reads associated with, skip
                    task_list.push_back(std::pair<std::string, std::vector<Vargas::SAM::Record>>(sub_rg_pair.first,
                                                                                                 alignment_reads.at(rgid)));
                    std::cerr << '\t' << sub_rg_pair.first << '\t' << alignment_reads.at(rgid).size() << '\n';
                    total += alignment_reads.at(rgid).size();
                }
            }

        }

        std::cerr << '\t' << alignment_reads.size() << " Read groups.\n"
                  << '\t' << alignment_rg_map.size() << " Subgraphs.\n"
                  << '\t' << task_list.size() << " Tasks.\n"
                  << '\t' << total << " Total alignments.\n"
                  << chrono_duration(start_time) << " seconds.\n" << std::endl;
    }

    std::cerr << "Loading graphs... \n" << std::endl;
    start_time = std::chrono::steady_clock::now();
    Vargas::GraphManager gm(gdf_file);
    std::cerr << chrono_duration(start_time) << " seconds." << std::endl;


    {
        Vargas::SAM::Header::Program pg;
        std::ostringstream ss;
        for (int i = 0; i < argc; ++i) ss << std::string(argv[i]) << " ";
        pg.command_line = ss.str();
        pg.name = "vargas_align";
        pg.id = "VA";
        pg.version = __DATE__;
        std::replace_if(pg.version.begin(), pg.version.end(), isspace, ' '); // rm tabs
        reads_hdr.add(pg);
    }


    std::cerr << "Aligning... " << std::endl;
    start_time = std::chrono::steady_clock::now();
    auto start_cpu = std::clock();

    const size_t num_tasks = task_list.size();

    Vargas::osam aligns_out(out_file, reads_hdr);

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
        Vargas::ByteAligner aligner(gm.node_len(), read_len, match, mismatch, gopen, gext);
        task_list.at(l).first;
        auto subgraph = gm.make_subgraph(task_list.at(l).first);
        auto aligns = aligner.align(read_seqs, targets, subgraph->begin(), subgraph->end());
        for (size_t j = 0; j < task_list.at(l).second.size(); ++j) {
            Vargas::SAM::Record &rec = task_list.at(l).second.at(j);
            rec.ref_name = task_list.at(l).first;
            rec.aux.set(ALIGN_SAM_MAX_POS_TAG, (int) aligns.max_pos[j]);
            rec.aux.set(ALIGN_SAM_MAX_SCORE_TAG, aligns.max_score[j]);
            rec.aux.set(ALIGN_SAM_MAX_COUNT_TAG, aligns.max_count[j]);
            rec.aux.set(ALIGN_SAM_SUB_POS_TAG, (int) aligns.sub_pos[j]);
            rec.aux.set(ALIGN_SAM_SUB_SCORE_TAG, aligns.sub_score[j]);
            rec.aux.set(ALIGN_SAM_SUB_COUNT_TAG, aligns.sub_count[j]);
            rec.aux.set(ALIGN_SAM_COR_FLAG_TAG, aligns.cor_flag[j]);
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
    std::cerr << chrono_duration(start_time, end_time) << " s\n"
              << cput << " CPU s\n"
              << cput / total << " CPU s / alignment.\n" << std::endl;

    return 0;
}

int split_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        split_help();
        return 0;
    }

    std::string sam_file = "", prefix = "";
    size_t lines = 0, num_files = 0;
    args >> GetOpt::Option('s', "sam", sam_file)
         >> GetOpt::Option('p', "prefix", prefix)
         >> GetOpt::Option('l', "lines", lines)
         >> GetOpt::Option('n', "num", num_files);

    if ((lines && num_files) || (!lines && !num_files))
        throw std::invalid_argument("One of -n and -l should be defined.");

    if (prefix.length() == 0 && sam_file.length() == 0) prefix = "sam_split.";
    else if (prefix.length() == 0) prefix = sam_file + ".";

    size_t suffix = 0;
    Vargas::isam input(sam_file);

    auto start_time = std::chrono::steady_clock::now();

    if (lines) {
        std::ofstream out;
        size_t l = lines;
        do {
            if (l == lines) {
                out.close();
                out.open(prefix + std::to_string(suffix));
                std::cerr << "\"" << prefix << std::to_string(suffix) << "\" ..." << std::endl;
                out << input.header().to_string();
                ++suffix;
                l = 0;
            }
            out << input.record().to_string() << "\n";
            ++l;
        } while (input.next());

    }

    if (num_files) {
        std::vector<std::ofstream *> outs;
        for (size_t i = 0; i < num_files; ++i) {
            outs.push_back(new std::ofstream(prefix + std::to_string(suffix++)));
            if (!outs.back()->good())
                throw std::invalid_argument("Error opening output file \"" + prefix + std::to_string(suffix) + "\"");
            *(outs.back()) << input.header().to_string();
        }
        suffix = 0;
        do {
            *(outs[suffix++]) << input.record().to_string();
            if (suffix == outs.size()) suffix = 0;
        } while (input.next());

        for (auto p : outs) {
            p->close();
            delete p;
        }
    }

    std::cerr << std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::steady_clock::now() - start_time).count()
              << " seconds." << std::endl;

    return 0;
}

int merge_main(const int argc, const char *argv[]) {

    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        merge_help();
        return 0;
    }

    if (argc < 3) throw std::invalid_argument("No SAM files specified.");

    std::string out_file = "";
    args >> GetOpt::Option('t', "out", out_file);

    auto start_time = std::chrono::steady_clock::now();

    Vargas::isam in(argv[2]);
    Vargas::osam out(out_file, in.header());

    for (int j = 2; j < argc; ++j) {
        Vargas::isam in(argv[j]);
        if (!in.good()) throw std::invalid_argument("Error opening SAM file \"" + std::string(argv[j]) + "\"");
        do {
            out.add_record(in.record());
        } while (in.next());
    }

    std::cerr << std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::steady_clock::now() - start_time).count()
              << " seconds." << std::endl;

    return 0;
}

int sam2csv(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        sam2csv_help();
        return 0;
    }

    std::string sam_file = "", format = "";

    args >> GetOpt::Option('s', "sam", sam_file)
         >> GetOpt::Option('f', "format", format);

    auto start_time = std::chrono::steady_clock::now();

    format.erase(std::remove(format.begin(), format.end(), ' '), format.end());
    std::vector<std::string> fmt_split = split(format, ',');

    std::unordered_set<std::string> warned;

    Vargas::isam input(sam_file);

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

    Vargas::GraphBuilder gb(fasta);
    gb.open_vcf(bcf);
    gb.region(region);

    auto start = std::clock();

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

        Vargas::ByteAligner a(g.max_node_len(), 50);

        std::cerr << SIMDPP_FAST_INT8_SIZE << " read alignment:\n";


        {
            start = std::clock();
            Vargas::Graph g2(g, g.subset(ingroup));
            std::cerr << "\tDerived Graph (" << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s)\n\t";
            start = std::clock();
            Vargas::ByteAligner::Results aligns = a.align(reads, g2.begin(), g2.end());
            std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
        }

    }
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
    cerr << "\tsplit       Split a SAM file into multiple files.\n";
    cerr << "\tmerge       Merge SAM files.\n";
    cerr << "\tconvert     Convert a SAM file to a CSV file.\n";
    cerr << "\ttest        Run doctests.\n";
    cerr << "\tprofile     Run profiles.\n" << endl;
}

void define_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
         << "-------------------- Vargas define, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-f\t--fasta         *<string> Reference filename.\n";
    cerr << "-v\t--vcf           <string> VCF or BCF file.\n";
    cerr << "-k\t--ksnp          <string> KSNP File. Only one of -v or -k should be defined.\n";
    cerr << "-g\t--region        *<string> Region of graph, format CHR:MIN-MAX.\n";
    cerr << "-l\t--nodelen       <int> Max node length, default 1,000,000\n";
    cerr << "-s\t--subgraph      *<string> Subgraph definition file. Default stdin.\n";
    cerr << "-t\t--out           <string> Output file, default stdout.\n";
    cerr << "-d\t--dot           <string> Output DOT graph to file.\n" << endl;

    cerr << "Subgraphs are defined using the format \"label=N[%]\".\n"
         << "\'N\' is the number of samples / percentage of samples selected.\n"
         << "The samples are selected from the parent graph, scoped with \':\'.\n"
         << "The base graph is implied as the root for all labels. Example:\n"
         << "a=50;a:b=10;~a:c=5\n"
         << "\'~\' indicates the complement graph. \'BASE\' refers the the whole graph.\n" << endl;
}

void profile_help() {
    using std::cerr;
    using std::endl;
    cerr << endl
         << "---------------------- Vargas profile, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------"
         << endl;
    cerr << "-f\t--fasta         <string> Reference filename." << endl;
    cerr << "-v\t--var           <string> VCF/BCF filename." << endl;
    cerr << "-g\t--region        <string> Region of graph, format CHR:MIN-MAX." << endl;
    cerr << "-i\t--ingroup       <int> Percent of genotypes to include in alignment." << endl;
    cerr << "-s\t--string        <string,string..> Include reads in alignment. Rest will be random." << endl << endl;
}

void align_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
         << "------------------- Vargas align, " << __DATE__ << ". rgaddip1@jhu.edu -------------------\n";;
    cerr << "-g\t--gdef          *<string> Graph definition file.\n";
    cerr << "-r\t--reads         *<string, string...> Read files to align. Default stdin.\n";
    cerr << "-a\t--align         *<string> Alignment targets.\n";
    cerr << "-f\t--file          -a specifies a file name.\n";
    cerr << "-t\t--out           *<string> Alignment output file, default stdout.\n";
    cerr << "-l\t--rlen          <int> Max read length. Default 50.\n";
    cerr << "-m\t--match         <int> Match score, default 2.\n";
    cerr << "-n\t--mismatch      <int> Mismatch penalty, default 2.\n";
    cerr << "-o\t--gap_open      <int> Gap opening penalty, default 3.\n";
    cerr << "-e\t--gap_extend    <int> Gap extend penalty, default 1.\n";
    cerr << "-j\t--threads       <int> Number of threads. 0 for maximum hardware concurrency.\n";
    cerr << "            \t            Optimal reads per subgraph: n * j * " << SIMDPP_FAST_INT8_SIZE << endl << endl;
}

void sim_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
         << "-------------------- Vargas sim, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-g\t--gdef          *<string> Graph definition file. Default stdin.\n";
    cerr << "-s\t--sub           *<string;string; ...> list of graphs to simulate from.\n";
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

void split_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
         << "-------------------- Vargas split, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-s\t--sam           <string> SAM input file. Default stdin.\n";
    cerr << "-l\t--lines         <int> Maximum number of lines per file.\n";
    cerr << "-n\t--num           <int> Number of files to distribute records to.\n";
    cerr << "-p\t--prefix        <string> Prefix of output files. Numeric suffix will be appended.\n\n";
    cerr << "-n distributes records round-robin style.\n" << endl;

}

void merge_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
         << "-------------------- Vargas merge, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "Usage:\n";
    cerr << "\tvargas merge [-t out.sam] [f1.sam f2.sam ... fn.sam]\n";
    cerr << "-t\t--out           <string> Output file. Default stdout.\n";
    cerr << "\nNote: It is assumed that all SAM files share the same header, i.e. files were produced with\n"
         << "        vargas split.\n" << endl;
}

void sam2csv_help() {
    using std::cerr;
    using std::endl;

    cerr << endl
         << "-------------------- Vargas convert, " << __DATE__ << ". rgaddip1@jhu.edu --------------------\n";
    cerr << "-s\t--sam          <string> SAM input file. Default stdin.\n";
    cerr << "-f\t--format       *<string,string...> Specify tags per column. Case sensitive.\n";
    cerr << "\nOutput printed to stdout.\n";
    cerr << "Requred column names: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL\n";
    cerr << "Prefix with \"RG:\" to obtain a value from the associated read group.\n" << endl;

}