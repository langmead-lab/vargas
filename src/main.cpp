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

#define TIME_ALIGNMENT 1

#include "doctest.h"

#include <iostream>
#include <thread>
#include <algorithm>
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

    if (lines) {
        std::ofstream out;
        size_t l = lines;
        do {
            if (l == lines) {
                out.close();
                out.open(prefix + std::to_string(suffix));
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

    Vargas::isam in(argv[2]);
    Vargas::osam out(out_file, in.header());

    for (int j = 2; j < argc; ++j) {
        Vargas::isam in(argv[j]);
        if (!in.good()) throw std::invalid_argument("Error opening SAM file \"" + std::string(argv[j]) + "\"");
        do {
            out.add_record(in.record());
        } while (in.next());
    }

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
    std::string read_file, gdf_file;
    bool R = false, // Align to reference
        X = false, // Align to max-AF graph
        O = false, // Align to outgroup graph
        I = false; // Align to ingroup graph

    args >> GetOpt::Option('m', "match", match)
         >> GetOpt::Option('n', "mismatch", mismatch)
         >> GetOpt::Option('o', "gap_open", gopen)
         >> GetOpt::Option('e', "gap_extend", gext)
         >> GetOpt::Option('r', "reads", read_file)
         >> GetOpt::Option('g', "gdef", gdf_file)
         >> GetOpt::Option('j', "threads", threads)
         >> GetOpt::Option('l', "rlen", read_len)
         >> GetOpt::OptionPresent('R', R)
         >> GetOpt::OptionPresent('X', X)
         >> GetOpt::OptionPresent('O', O)
         >> GetOpt::OptionPresent('I', I);

    if (!R && !X && !O && !I) throw std::invalid_argument("No alignment mode [RXIO] selected.");
    if (threads == 0) threads = std::thread::hardware_concurrency();

    std::cerr << "Loading reads..." << std::endl;
    std::map<Vargas::Graph::GID, std::vector<Vargas::SAM::Record>> read_queue;
    Vargas::SAM::Header hdr;
    {
        Vargas::isam reads_in(read_file);
        hdr = reads_in.header();
        std::string gid_str;
        Vargas::Graph::GID gid;
        do {
            if (!reads_in.record().read_group(hdr, SIM_SAM_GID_TAG, gid_str)) gid = Vargas::Graph::GID();
            else gid = Vargas::Graph::GID(gid_str);
            read_queue[gid].push_back(reads_in.record());
            assert(reads_in.record().seq.length() == read_len);
        } while (reads_in.next());
    }
    size_t num_target_graphs = read_queue.size() * ((R ? 1 : 0) + (X ? 1 : 0) + (I ? 1 : 0) + (O ? 1 : 0));
    std::cerr << hdr.read_groups.size() << " read groups, " << num_target_graphs << " target graphs." << std::endl;


    std::map<Vargas::Graph::GID, Vargas::Graph::Population> pops;
    Vargas::Graph base_graph;

    {
        Vargas::GDEF gdf(gdf_file);
        if (O) gdf.include_outgroups();
        Vargas::GraphBuilder gb(gdf.fasta(), gdf.var());
        gb.region(gdf.region());

        std::cerr << "Loading base graph..." << std::endl;
        base_graph = gb.build();
        pops = gdf.populations();
    }

    {
        Vargas::SAM::Header::Program pg;
        std::ostringstream ss;
        ss << "vargas align ";
        for (int i = 0; i < argc; ++i) ss << std::string(argv[i]) << " ";
        pg.command_line = ss.str();
        pg.name = "vargas_align";
        pg.id = "VA";
        pg.version = __DATE__;
        std::replace_if(pg.version.begin(), pg.version.end(), isspace, ' '); // rm tabs
        hdr.add(pg);
    }


    std::vector<std::shared_ptr<Vargas::ByteAligner>> aligners;
    for (uint8_t i = 0; i < threads; ++i) {
        aligners.push_back(std::make_shared<Vargas::ByteAligner>(base_graph.max_node_len(),
                                                                 read_len,
                                                                 match,
                                                                 mismatch,
                                                                 gopen,
                                                                 gext));
    }


    std::cerr << "Aligning..." << std::endl;

    Vargas::osam aligns_out(hdr);
    std::vector<std::thread> jobs;
    std::vector<Vargas::ByteAligner::Results> aligns(threads);

    int prog_bar_scale = (num_target_graphs + 49) / 50;
    int scale_counter = 0;
    for (size_t i = 0; i < num_target_graphs / prog_bar_scale; ++i) std::cerr << '_';
    std::cerr << std::endl;

    Vargas::Graph ref_graph;
    Vargas::Graph af_graph;
    if (R) ref_graph = Vargas::Graph(base_graph, Vargas::Graph::REF);
    if (X) af_graph = Vargas::Graph(base_graph, Vargas::Graph::MAXAF);
    int alignment_type;

    #if TIME_ALIGNMENT
    time_t start = std::clock();
    #endif

    for (const auto &pair : read_queue) {
        auto gid = pair.first;

        std::vector<std::vector<Vargas::SAM::Record>> aligner_source(threads);
        std::vector<std::vector<std::string>> aligner_seqs(threads);
        std::vector<std::vector<uint32_t>> aligner_targets(threads);
        size_t aligner_idx = 0;

        for (const auto &s : pair.second) {
            // Distribute reads round-robin
            if (aligner_idx == threads) aligner_idx = 0;
            aligner_seqs[aligner_idx].push_back(s.seq);
            // SAM is left based coord. -1 gives last base coord
            aligner_targets[aligner_idx].push_back(s.pos + s.seq.length() - 1);
            aligner_source[aligner_idx].push_back(s);
            ++aligner_idx;
        }

        for (int i = 0; i < 4; ++i) {
            if (R && i == 0) {
                for (size_t j = 0; j < threads; ++j) {
                    jobs.emplace(jobs.end(),
                                 &Vargas::ByteAligner::align_into,
                                 aligners[jobs.size()],
                                 std::ref(aligner_seqs[jobs.size()]),
                                 std::ref(aligner_targets[jobs.size()]),
                                 ref_graph.begin(),
                                 ref_graph.end(),
                                 std::ref(aligns[jobs.size()]));
                }
                alignment_type = 0;
                std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });
            }
            else if (X && i == 1) {
                for (size_t j = 0; j < threads; ++j) {
                    jobs.emplace(jobs.end(),
                                 &Vargas::ByteAligner::align_into,
                                 aligners[jobs.size()],
                                 std::ref(aligner_seqs[jobs.size()]),
                                 std::ref(aligner_targets[jobs.size()]),
                                 af_graph.begin(),
                                 af_graph.end(),
                                 std::ref(aligns[jobs.size()]));
                }
                alignment_type = 1;
                std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });
            }
            else if (I && i == 2) {
                gid.outgroup = false;
                Vargas::Graph sub;
                try {
                    sub = Vargas::Graph(base_graph, pops.at(gid));
                } catch (std::out_of_range &e) {
                    throw std::out_of_range("Graph definition does not exist: \"" + gid.to_string() + "\"");
                }
                for (size_t j = 0; j < threads; ++j) {
                    jobs.emplace(jobs.end(),
                                 &Vargas::ByteAligner::align_into,
                                 aligners[jobs.size()],
                                 std::ref(aligner_seqs[jobs.size()]),
                                 std::ref(aligner_targets[jobs.size()]),
                                 sub.begin(),
                                 sub.end(),
                                 std::ref(aligns[jobs.size()]));
                }
                alignment_type = 2;
                std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });
            }
            else if (O && i == 3) {
                gid.outgroup = true;
                Vargas::Graph sub;
                try {
                    sub = Vargas::Graph(base_graph, pops.at(gid));
                } catch (std::out_of_range &e) {
                    throw std::out_of_range("Graph definition does not exist: \"" + gid.to_string() + "\"");
                }
                for (size_t j = 0; j < threads; ++j) {
                    jobs.emplace(jobs.end(),
                                 &Vargas::ByteAligner::align_into,
                                 aligners[jobs.size()],
                                 std::ref(aligner_seqs[jobs.size()]),
                                 std::ref(aligner_targets[jobs.size()]),
                                 sub.begin(),
                                 sub.end(),
                                 std::ref(aligns[jobs.size()]));
                }
                alignment_type = 3;
                std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });
            }

#if !TIME_ALIGNMENT
            for (size_t k = 0; k < threads; ++k) {
                for (size_t j = 0; j < aligner_source[k].size(); ++j) {
                    auto &rec = aligner_source[k][j];
                    switch (alignment_type) {
                        case 0:
                            rec.aux.set(ALIGN_SAM_TYPE_TAG, ALIGN_SAM_TYPE_REF);
                            break;
                        case 1:
                            rec.aux.set(ALIGN_SAM_TYPE_TAG, ALIGN_SAM_TYPE_MAXAF);
                            break;
                        case 2:
                            rec.aux.set(ALIGN_SAM_TYPE_TAG, ALIGN_SAM_TYPE_IN);
                            break;
                        case 3:
                            rec.aux.set(ALIGN_SAM_TYPE_TAG, ALIGN_SAM_TYPE_OUT);
                            break;
                        default:
                            throw std::invalid_argument("Invalid type: " + std::to_string(alignment_type));
                    }

                    rec.aux.set(ALIGN_SAM_MAX_POS_TAG, (int) aligns[k].max_pos[j]);
                    rec.aux.set(ALIGN_SAM_MAX_SCORE_TAG, aligns[k].max_score[j]);
                    rec.aux.set(ALIGN_SAM_MAX_COUNT_TAG, aligns[k].max_count[j]);
                    rec.aux.set(ALIGN_SAM_SUB_POS_TAG, (int) aligns[k].sub_pos[j]);
                    rec.aux.set(ALIGN_SAM_SUB_SCORE_TAG, aligns[k].sub_score[j]);
                    rec.aux.set(ALIGN_SAM_SUB_COUNT_TAG, aligns[k].sub_count[j]);
                    rec.aux.set(ALIGN_SAM_COR_FLAG_TAG, aligns[k].cor_flag[j]);
                    aligns_out.add_record(rec);
                }
            }

           if (++scale_counter == prog_bar_scale) {
                std::cerr << "\u2588" << std::flush;
                scale_counter = 0;
            }
#endif
            jobs.clear();

        }
    }

#if TIME_ALIGNMENT
    std::cerr << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
#endif

    std::cerr << std::endl;

    return 0;
}

int sim_main(const int argc, const char *argv[]) {
    GetOpt::GetOpt_pp args(argc, argv);

    if (args >> GetOpt::OptionPresent('h', "help")) {
        sim_help();
        return 0;
    }

    // Load parameters
    unsigned int read_len = 50, num_reads = 1000, threads = 1;
    std::string mut = "0", indel = "0", vnodes = "*", vbases = "*";
    std::string gdf_file, out_file = "";
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
         >> GetOpt::Option('j', "threads", threads)
         >> GetOpt::OptionPresent('a', "rate", use_rate);


    Vargas::GDEF gdf(gdf_file);
    if (outgroup) gdf.include_outgroups();
    auto &pops = gdf.populations();

    Vargas::GraphBuilder gb(gdf.fasta(), gdf.var());
    gb.region(gdf.region());

    auto mut_split = split(mut, ',');
    auto indel_split = split(indel, ',');
    auto vnode_split = split(vnodes, ',');
    auto vbase_split = split(vbases, ',');

    std::cerr << "Building profiles..." << std::endl;

    // Map a read group ID to a unique read group
    std::unordered_map<std::string, std::pair<Vargas::Graph::Population, Vargas::Sim::Profile>> pending_sims;

    int rg_id = 0;
    Vargas::SAM::Header::ReadGroup rg;
    rg.seq_center = "vargas_sim";
    rg.date = current_date();
    rg.aux.set(SIM_SAM_REF_TAG, gdf.fasta());
    rg.aux.set(SIM_SAM_VCF_TAG, gdf.var());

    for (auto &vbase : vbase_split) {
        for (auto &vnode : vnode_split) {
            for (auto &ind : indel_split) {
                for (auto &m : mut_split) {
                    Vargas::Sim::Profile prof;

                    try {
                        prof.len = read_len;
                        prof.mut = m == "*" ? -1 : std::stof(m);
                        prof.indel = ind == "*" ? -1 : std::stof(ind);
                        prof.rand = use_rate;
                        prof.var_bases = vbase == "*" ? -1 : std::stoi(vbase);
                        prof.var_nodes = vnode == "*" ? -1 : std::stoi(vnode);
                    } catch (std::out_of_range &e) {
                        std::cerr << "Invalid profile argument." << std::endl;
                        return 1;
                    }

                    rg.aux.set(SIM_SAM_INDEL_ERR_TAG, prof.indel);
                    rg.aux.set(SIM_SAM_VAR_NODES_TAG, prof.var_nodes);
                    rg.aux.set(SIM_SAM_VAR_BASE_TAG, prof.var_bases);
                    rg.aux.set(SIM_SAM_SUB_ERR_TAG, prof.mut);
                    rg.aux.set(SIM_SAM_USE_RATE_TAG, (int) prof.rand);

                    // Each profile and subgraph combination is a unique set of reads
                    for (auto &p : pops) {
                        rg.aux.set(SIM_SAM_GID_TAG, p.first.to_string());
                        rg.aux.set(SIM_SAM_POPULATION, p.second.to_string());
                        rg.id = std::to_string(++rg_id);
                        sam_hdr.add(rg);

                        pending_sims[rg.id] = std::pair<Vargas::Graph::Population,
                                                        Vargas::Sim::Profile>(p.second, prof);
                    }
                }
            }
        }
    }
    std::cerr << pending_sims.size() << " read groups over " << pops.size() << " subgraphs." << std::endl;

    Vargas::osam out(out_file, sam_hdr);
    if (!out.good()) throw std::invalid_argument("Error opening output file \"" + out_file + "\"");

    std::cerr << "Loading base graph..." << std::endl;
    Vargas::Graph base_graph = gb.build();

    int prog_bar_scale = ((pending_sims.size() / threads) + 49) / 50;
    int scale_counter = 0;

    for (size_t i = 0; i < (pending_sims.size() / threads) / prog_bar_scale; ++i) std::cerr << "_";
    std::cerr << std::endl;

    // Threading
    std::vector<std::vector<Vargas::SAM::Record>> reads(threads);
    std::vector<std::thread> jobs;

    // For each readgroup
    //TODO Building the same graph multiple times (one for each prof)
    unsigned int num = 0;
    for (auto &ps : pending_sims) {
        ++num;
        jobs.emplace(jobs.end(),
                     &derive_and_sim,
                     ps.first, std::ref(base_graph),
                     std::ref(ps.second.first),
                     std::ref(ps.second.second),
                     num_reads,
                     std::ref(reads[jobs.size()]));

        if (jobs.size() == threads || num == pending_sims.size()) {
            std::for_each(jobs.begin(), jobs.end(), [](std::thread &t) { t.join(); });
            if (++scale_counter == prog_bar_scale) {
                std::cerr << "\u2588" << std::flush;
                scale_counter = 0;
            }
            for (auto &t : reads) {
                for (auto &r : t) {
                    out.add_record(r);
                }
                jobs.clear();
            }
        }
    }
    std::cerr << std::endl;
    return 0;
}

void derive_and_sim(std::string rg,
                    const Vargas::Graph &base,
                    const Vargas::Graph::Population &pop,
                    const Vargas::Sim::Profile &prof,
                    unsigned int num_reads,
                    std::vector<Vargas::SAM::Record> &results) {
    Vargas::Sim sim(base, pop, prof);
    results = sim.get_batch(num_reads);
    for (auto &r : results) r.aux.set("RG", rg);
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
         << "---------------------- Vargas profile, " << __DATE__ << ". rgaddip1@jhu.edu ----------------------"
         << endl;
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
    cerr << "-g\t--gdef          *<string> Graph definition file. Reads are simulated from the ingroups.\n";
    cerr << "-t\t--out           <string> Output file. Default stdout.\n";
    cerr << "-o\t--outgroup      Simulate from outgroup graphs.\n";
    cerr << "-n\t--numreads      <int> Number of reads to simulate from each subgraph, default 1000.\n";
    cerr << "-m\t--muterr        <int/float, int/float...> Read mutation error. Default 0.\n";
    cerr << "-i\t--indelerr      <int/float, int/float...> Read indel error. Default 0.\n";
    cerr << "-v\t--vnodes        <int, int...> Number of variant nodes, default any (*).\n";
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