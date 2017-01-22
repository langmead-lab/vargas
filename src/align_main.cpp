/**
 * Ravi Gaddipati
 * Dec 23, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Main aligner interface
 *
 * @file
 */

#include "align_main.h"
#include "alignment.h"
#include "sam.h"
#include "gdef.h"
#include "doctest.h"
#include "sim.h"
#include "cxxopts.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


int align_main(int argc, char *argv[]) {
    std::string cl;
    {
        std::ostringstream ss;
        for (int i = 0; i < argc; ++i) ss << std::string(argv[i]) << " ";
        cl = ss.str();
    }

    // Load parameters
    unsigned match, mismatch, gopen, gext, threads, tolerance, chunk_size, subsample = 0;
    std::string read_file, gdf, align_targets, out_file, pgid;
    bool align_targets_isfile = false, end_to_end = false;

    cxxopts::Options opts("vargas align", "Align reads to a graph.");
    try {
        opts.add_options()
        ("g,gdef", "<str> *Graph definition file.", cxxopts::value(gdf))
        ("r,reads", "<str> SAM reads file. (default: stdin)", cxxopts::value(read_file))
        ("p,subsample", "<N> Subsample N random reads.", cxxopts::value(subsample))
        ("a,align", "<str> Alignment targets/file of form \"RG:[ID][gd],target\"", cxxopts::value(align_targets))
        ("f,file", " -a specifies a file name.", cxxopts::value(align_targets_isfile))
        ("m,match", "<N> Match score.", cxxopts::value(match)->default_value("2"))
        ("n,mismatch", "<N> Mismatch penalty.", cxxopts::value(mismatch)->default_value("2"))
        ("o,gap_open", "<N> Gap opening penalty.", cxxopts::value(gopen)->default_value("3"))
        ("e,gap_extend", "<N> Gap extension penalty.", cxxopts::value(gext)->default_value("1"))
        ("x,endtoend", "Perform end to end alignment", cxxopts::value(end_to_end))
        ("c,tolerance", "<N> Correct if within readlen/N.",
         cxxopts::value(tolerance)->default_value(std::to_string(vargas::Aligner::default_tolerance())))
        ("u,chunk", "<N> Partition tasks into chunks with max size N.",
         cxxopts::value(chunk_size)->default_value("2048"))
        ("t,out", "<str> Output file. Can safely be same as -r. (default: stdout)", cxxopts::value(out_file))
        ("s,assess", "[ID] Get score profile from PG with ID.", cxxopts::value(pgid)->implicit_value("-"))
        ("j,threads", "<N> Number of threads.", cxxopts::value(threads)->default_value("1"))
        ("h,help", "Display this message.");
        opts.parse(argc, argv);
    } catch (std::exception &e) { throw std::invalid_argument("Error parsing options: " + std::string(e.what())); }
    if (opts.count("h")) {
        align_help(opts);
        return 0;
    }
    if (!opts.count("gdef")) throw std::invalid_argument("Graph definition file required.");

    if (chunk_size < vargas::Aligner::read_capacity() ||
    chunk_size % vargas::Aligner::read_capacity() != 0) {
        std::cerr << "WARN: Chunk size is not a multiple of SIMD vector length: "
                  << vargas::Aligner::read_capacity() << std::endl;
    }

    #ifndef _OPENMP
    // Disable threads if no openMP.
    if (threads != 1) {
        std::cerr << "WARN: Threads specified without OpenMP Compilation." << std::endl;
    }
    threads = 1;
    #endif

    if (align_targets_isfile) {
        std::ifstream in(align_targets);
        if (!in.good()) throw std::invalid_argument("Invalid alignment targets file \"" + align_targets + "\".");
        std::stringstream ss;
        ss << in.rdbuf();
        align_targets = ss.str();
    }

    auto reads = vargas::isam(read_file).subset(subsample);
    auto &reads_hdr = reads.header();

    vargas::ScoreProfile prof(match, mismatch, gopen, gext);
    if (pgid == "-") {
        bool check = false;
        for (const auto &i : reads.header().programs) {
            if (std::find(vargas::supported_pgid.begin(), vargas::supported_pgid.end(), i.first)
            != vargas::supported_pgid.end()) {
                pgid = i.first;
                prof = vargas::program_profile(i.second.command_line);
                check = true;
                break;
            }
        }
        if (!check) throw std::invalid_argument("No suitable scoring profile found in SAM program header.");
        std::cerr << "Using profile for: " << pgid << "\n";
    } else if (pgid.length()) {
        try {
            prof = vargas::program_profile(reads_hdr.programs.at(pgid).command_line);
        } catch (std::exception &e) { throw std::invalid_argument("Unrecognized PG ID: " + pgid); }
    } else {
        prof.end_to_end = end_to_end;
        prof.ambig = 0;
    }
    prof.tol = tolerance;
    std::cerr << "Scoring profile: " << prof.to_string() << "\n";

    vargas::SAM::Header::Program pg;
    pg.command_line = cl;
    pg.name = "vargas_align";
    pg.id = "VA";
    pg.version = __DATE__;
    std::replace_if(pg.version.begin(), pg.version.end(), isspace, ' '); // rm tabs
    const auto assigned_pgid = reads_hdr.add(pg);

    size_t read_len;
    bool padded;
    auto task_list = create_tasks(reads, align_targets, chunk_size, read_len, padded);
    if (prof.end_to_end && prof.ambig && padded) {
        std::cerr << "WARN: ETE alignment with N penalty will impact padded read scores.\n";
    }

    const size_t num_tasks = task_list.size();
    if (num_tasks < threads) {
        std::cerr << "WARN: Number of threads is greater than number of tasks. Try decreasing --chunk.\n";
    }

    #ifdef _OPENMP
    if (threads) threads = threads > task_list.size() ? task_list.size() : threads;
    omp_set_num_threads(threads);
    #endif

    const bool use_wide = read_len * match > 255;
    if (use_wide) {
        std::cerr << "Maximum possible score: " << read_len * match << ". Using 16-bit aligner ("
                  << vargas::WordAligner::read_capacity() << " reads/vector).\n";
    }


    std::vector<std::unique_ptr<vargas::AlignerBase>> aligners(threads);
    for (size_t k = 0; k < threads; ++k) {
        aligners[k] = make_aligner(prof, read_len, use_wide);
    }

    auto files = rg::split(gdf, ',');
    for (const auto &gdef : files) {
        std::cerr << "\nLoading \"" << gdef << "\"... ";
        auto start_time = std::chrono::steady_clock::now();
        vargas::GraphManager gm(gdef);
        std::cerr << "(" << gm.base()->node_map()->size() << " nodes), "
                  << rg::chrono_duration(start_time) << "s.\n";

        align(gm, task_list, aligners);

        std::string fname = out_file;
        if (fname.length() > 0 && files.size() > 1) {
            auto ld = out_file.find_last_of('.');
            fname = fname.substr(0, ld) + "_" + gdef.substr(0, gdef.find_last_of('.'));
            if (ld != std::string::npos) fname += out_file.substr(ld);
        }
        std::cerr << "Writing to \"" << fname << "\".\n";
        reads_hdr.programs[assigned_pgid].aux.set(ALIGN_SAM_PG_GDF, gdef);
        vargas::osam aligns_out(fname, reads_hdr);
        for (size_t l = 0; l < num_tasks; ++l) {
            for (size_t j = 0; j < task_list.at(l).second.size(); ++j) {
                aligns_out.add_record(task_list.at(l).second.at(j));
            }
        }
    }

    return 0;
}

void align(vargas::GraphManager &gm, std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>> &task_list,
           const std::vector<std::unique_ptr<vargas::AlignerBase>> &aligners) {

    std::cerr << "Aligning... " << std::flush;
    auto start_cpu = std::clock();
    auto start_time = std::chrono::steady_clock::now();

    const auto num_tasks = task_list.size();

    #pragma omp parallel for
    for (size_t l = 0; l < num_tasks; ++l) {
        #ifdef _OPENMP
        const int tid = omp_get_thread_num();
        #else
        const int tid = 0;
        #endif
        const size_t num_reads = task_list.at(l).second.size();
        std::vector<std::string> read_seqs(num_reads);
        std::vector<unsigned> targets(num_reads);
        for (size_t i = 0; i < num_reads; ++i) {
            const auto &r = task_list.at(l).second.at(i);
            read_seqs[i] = r.seq;
            targets[i] = r.pos;
            if (r.pos > 0) {
                targets[i] = r.pos - 1;
                if (r.cigar.size()) {
                    for (const auto &p : r.cigar) {
                        if (p.second == 'M' || p.second == 'D' || p.second == '=' || p.second == 'X')
                            targets[i] += p.first;
                    }
                    // If no cigar is present, use whole read length offset
                } else targets[i] = r.pos + r.seq.length() - 1;
            } else targets[i] = 0;
        }
        auto subgraph = gm.make_subgraph(task_list.at(l).first);
        const auto aligns = aligners[tid]->align(read_seqs, targets, subgraph->begin(), subgraph->end());
        for (size_t j = 0; j < task_list.at(l).second.size(); ++j) {
            vargas::SAM::Record &rec = task_list.at(l).second.at(j);
            rec.aux.set(ALIGN_SAM_MAX_POS_TAG, aligns.max_pos[j]);
            rec.aux.set(ALIGN_SAM_MAX_SCORE_TAG, aligns.max_score[j]);
            rec.aux.set(ALIGN_SAM_SCORE_PROFILE, aligns.profile.to_string());
            rec.aux.set(ALIGN_SAM_SEQ, subgraph->region().seq_name);
            rec.aux.set(ALIGN_SAM_MAX_COUNT_TAG, aligns.max_count[j]);
            rec.aux.set(ALIGN_SAM_SUB_POS_TAG, aligns.sub_pos[j]);
            rec.aux.set(ALIGN_SAM_SUB_SCORE_TAG, aligns.sub_score[j]);
            rec.aux.set(ALIGN_SAM_SUB_COUNT_TAG, aligns.sub_count[j]);
            rec.aux.set(ALIGN_SAM_COR_FLAG_TAG, aligns.correct[j]);
            rec.aux.set(ALIGN_SAM_TARGET_SCORE, aligns.target_score[j]);
        }
    }

    auto end_time = std::chrono::steady_clock::now();
    auto cput = (std::clock() - start_cpu) / (double) CLOCKS_PER_SEC;
    std::cerr << rg::chrono_duration(start_time, end_time) << "s, "
              << cput << " CPUs.\n";

}


std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>>
create_tasks(vargas::isam &reads, std::string &align_targets,
             const size_t chunk_size, size_t &read_len, bool &resized) {
    std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>> task_list;
    std::unordered_map<std::string, std::vector<vargas::SAM::Record>> read_groups;

    std::vector<std::string> alignment_pairs;
    if (align_targets.length() != 0) {
        std::replace(align_targets.begin(), align_targets.end(), '\n', ';');
        alignment_pairs = rg::split(align_targets, ';');
    }

    std::cerr << "Loading reads... " << std::flush;
    auto start_time = std::chrono::steady_clock::now();

    size_t total = 0;
    auto &reads_hdr = reads.header();
    std::string read_group;
    vargas::SAM::Record rec;
    read_len = reads.record().seq.length();
    do {
        rec = reads.record();
        if (rec.seq.length() > read_len) read_len = rec.seq.length();
        if (!rec.aux.get("RG", read_group)) {
            read_group = UNGROUPED_READGROUP;
            rec.aux.set("RG", UNGROUPED_READGROUP);
            if (!reads_hdr.read_groups.count(UNGROUPED_READGROUP)) {
                reads_hdr.add(vargas::SAM::Header::ReadGroup("@RG\tID:" + std::string(UNGROUPED_READGROUP)));
            }
        }
        read_groups[read_group].push_back(rec);
    } while (reads.next());

    // Pad short reads
    resized = false;
    for (auto &rg : read_groups) {
        for (auto &rd : rg.second) {
            if (rd.seq.size() != read_len) {
                rd.seq.resize(read_len, 'N');
                resized = true;
            }
        }
    }

    if (alignment_pairs.size() == 0) {
        for (const auto &p : read_groups) {
            alignment_pairs.push_back("RG:ID:" + p.first + "\t" + vargas::GraphManager::GDEF_BASEGRAPH);
        }
    }

    // Maps target graph to read group ID's
    std::unordered_map<std::string, std::vector<std::string>> alignment_rg_map;

    std::vector<std::string> pair;
    std::string tag, val, target_val;
    for (const std::string &p : alignment_pairs) {
        rg::split(p, pair);
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


    std::cerr << rg::chrono_duration(start_time) << "s." << std::endl;

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
                if (safe_beg != safe_end)
                    task_list.emplace_back(sub_rg_pair.first, std::vector<vargas::SAM::Record>(safe_beg, safe_end));
            }
        }
    }

    std::cerr << read_groups.size() << "\tRead group(s).\n"
              << alignment_rg_map.size() << "\tSubgraph(s).\n"
              << task_list.size() << "\tTask(s).\n"
              << total << "\tTotal alignments.\n"
              << read_len << "\tMax read length.\n";

    return task_list;
}

std::unique_ptr<vargas::AlignerBase> make_aligner(const vargas::ScoreProfile &prof, size_t read_len,
                                                  bool use_wide) {
    if (prof.end_to_end) {
        if (use_wide) return rg::make_unique<vargas::WordAlignerETE>(read_len, prof);
        else return rg::make_unique<vargas::AlignerETE>(read_len, prof);
    } else {
        if (use_wide) return rg::make_unique<vargas::WordAligner>(read_len, prof);
        else return rg::make_unique<vargas::Aligner>(read_len, prof);
    }
}

void align_help(const cxxopts::Options &opts) {
    using std::cerr;
    using std::endl;

    cerr << opts.help() << "\n" << endl;
    cerr << "Elements per SIMD vector: " << vargas::Aligner::read_capacity() << endl;
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

    vargas::Aligner aligner(5);
    auto reads = sim.get_batch(aligner.read_capacity());

    std::vector<std::string> seqs;
    std::vector<unsigned> targets;
    for (auto &r : reads) {
        seqs.push_back(r.seq);
        targets.push_back(r.pos + r.seq.length() - 1);
    }

    auto results = aligner.align(seqs, targets, g.begin(), g.end());

    for (auto i : results.correct) CHECK ((int) i == 1);

    remove(tmpfa.c_str());
    remove(tmpvcf.c_str());
    remove((tmpfa + ".fai").c_str());
}
TEST_CASE ("Correctness flag") {
    try {
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

        vargas::Aligner aligner(6);
        vargas::isam reads(reads_file);

        std::vector<vargas::SAM::Record> records;
        std::vector<std::string> read_seq;
        std::vector<unsigned> targets;
        do {
            records.push_back(reads.record());
            read_seq.push_back(reads.record().seq);
            targets.push_back(reads.record().pos + read_seq.back().length() - 1);
        } while (reads.next());

        auto res = aligner.align(read_seq, targets, g.begin(), g.end());

        vargas::osam align_out("tmp_aout.sam", reads.header());
        for (auto &r : records) align_out.add_record(r);

        remove(tmpfa.c_str());
        remove((tmpfa + ".fai").c_str());
        remove(tmpvcf.c_str());
        remove(reads_file.c_str());
        remove("tmp_aout.sam");
    } catch (std::exception &e) {
        std::cerr << e.what();
        throw;
    }
}

TEST_SUITE_END();