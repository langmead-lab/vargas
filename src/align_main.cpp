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
    unsigned match, mismatch, gopen, gext, threads, tolerance, chunk_size, subsample;
    std::string read_file, gdf, align_targets, out_file, pgid;
    bool end_to_end = false;

    cxxopts::Options opts("vargas align", "Align reads to a graph.");
    try {
        opts.add_options()
        ("g,gdef", "<str> *Graph definition file.", cxxopts::value(gdf))
        ("r,reads", "<str> *Unpaired reads in SAM, FASTQ, or FASTA format.", cxxopts::value(read_file))
        ("p,subsample", "<N> Subsample N random reads, 0 for all.", cxxopts::value(subsample)->default_value("0"))
        ("a,align", "<str> Alignment targets/file of form \"RG:[ID][gd],target;...\"", cxxopts::value(align_targets))
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
    if (!opts.count("reads")) throw std::invalid_argument("No read file provided.");
    ReadFmt format = read_fmt(read_file);

    if (chunk_size < vargas::Aligner::read_capacity() ||
    chunk_size % vargas::Aligner::read_capacity() != 0) {
        std::cerr << "WARN: Chunk size is not a multiple of SIMD vector length: "
                  << vargas::Aligner::read_capacity() << std::endl;
    }

    if (align_targets.size() > 0 && format != ReadFmt::SAM) {
        throw std::invalid_argument("Alignment targets only available for SAM inputs.");
    }
    // If targets doesn't have a comma, assume its a filename
    if (align_targets.size() > 0 && std::find(align_targets.begin(), align_targets.end(), ',') == align_targets.end()) {
        std::ifstream in(align_targets);
        if (!in.good()) throw std::invalid_argument("Invalid alignment targets file \"" + align_targets + "\".");
        std::stringstream ss;
        ss << in.rdbuf();
        align_targets = ss.str();
    }

    vargas::isam reads;
    if (format == ReadFmt::FASTQ) {
        load_fast(read_file, true, reads);
    } else if (format == ReadFmt::FASTA) {
        load_fast(read_file, false, reads);
    } else {
        reads.open(read_file);
    }
    reads.subset(subsample);
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
    #else
    threads = 1;
    #endif

    const bool use_wide = read_len * match > 255;
    if (use_wide) {
        std::cerr << "Maximum possible score: " << read_len * match << ". Using 16-bit aligner ("
                  << vargas::WordAligner::read_capacity() << " reads/vector).\n";
    }

    prof.tol = read_len / tolerance;
    std::cerr << "Scoring profile: " << prof.to_string() << "\n";

    std::vector<std::unique_ptr<vargas::AlignerBase>> aligners(threads);
    for (size_t k = 0; k < threads; ++k) {
        aligners[k] = make_aligner(prof, read_len, use_wide);
    }

    auto files = rg::split(gdf, ',');
    for (const auto &gdef : files) {
        std::cerr << "\nLoading \"" << gdef << "\"... ";
        auto start_time = std::chrono::steady_clock::now();
        vargas::GraphGen gm(gdef);
        gm.dump();
        std::cerr << rg::chrono_duration(start_time) << "s.\n";

        align(gm, task_list, aligners);

        std::string fname = out_file;
        if (fname.length() > 0 && files.size() > 1) {
            auto ld = out_file.find_last_of('.');
            fname = fname.substr(0, ld) + "_" + gdef.substr(0, gdef.find_last_of('.'));
            if (ld != std::string::npos) fname += out_file.substr(ld);
        }
        if (fname.length()) std::cerr << "Writing to \"" << (fname == "" ? "stdout" : fname) << "\".\n";
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

void align(vargas::GraphGen &gm, std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>> &task_list,
           const std::vector<std::unique_ptr<vargas::AlignerBase>> &aligners) {

    std::cerr << "Aligning... " << std::flush;
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
        auto subgraph = gm.at(task_list.at(l).first);

        const auto aligns = aligners[tid]->align(read_seqs, targets, subgraph->begin(), subgraph->end());
        for (size_t j = 0; j < task_list.at(l).second.size(); ++j) {
            vargas::SAM::Record &rec = task_list.at(l).second.at(j);
            auto abs = gm.absolute_position(aligns.max_pos[j]);
            rec.aux.set(ALIGN_SAM_MAX_POS_TAG, abs.second);
            rec.aux.set(ALIGN_SAM_MAX_SEQ, abs.first);

            abs = gm.absolute_position(aligns.sub_pos[j]);
            rec.aux.set(ALIGN_SAM_SUB_SEQ, abs.first);
            rec.aux.set(ALIGN_SAM_SUB_POS_TAG, abs.second);

            rec.aux.set(ALIGN_SAM_MAX_COUNT_TAG, aligns.max_count[j]);
            rec.aux.set(ALIGN_SAM_SUB_COUNT_TAG, aligns.sub_count[j]);

            rec.aux.set(ALIGN_SAM_MAX_SCORE_TAG, aligns.max_score[j]);
            rec.aux.set(ALIGN_SAM_SCORE_PROFILE, aligns.profile.to_string());
            rec.aux.set(ALIGN_SAM_SUB_SCORE_TAG, aligns.sub_score[j]);
            if (targets.at(j) != 0) {
                rec.aux.set(ALIGN_SAM_COR_FLAG_TAG, aligns.correct[j]);
                rec.aux.set(ALIGN_SAM_TARGET_SCORE, aligns.target_score[j]);
            }
        }
    }

    std::cerr << rg::chrono_duration(start_time) << "s.\n";

}

std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>>
create_tasks(vargas::isam &reads, std::string &align_targets, const int chunk_size, size_t &read_len, bool &resized) {
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
            alignment_pairs.push_back("RG:ID:" + p.first + ",base");
        }
    }

    // Maps target graph to read group ID's
    std::unordered_map<std::string, std::vector<std::string>> alignment_rg_map;

    std::string tag, val, target_val;
    for (const std::string &p : alignment_pairs) {
        auto pair = rg::split(p, ',');
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

std::unique_ptr<vargas::AlignerBase>
make_aligner(const vargas::ScoreProfile &prof, size_t read_len, bool use_wide) {
    if (prof.end_to_end) {
        if (use_wide) return rg::make_unique<vargas::WordAlignerETE>(read_len, prof);
        else return rg::make_unique<vargas::AlignerETE>(read_len, prof);
    } else {
        if (use_wide) return rg::make_unique<vargas::WordAligner>(read_len, prof);
        else return rg::make_unique<vargas::Aligner>(read_len, prof);
    }
}

void load_fast(std::string &file, const bool fastq, vargas::isam &ret) {
    std::string input;
    if (file.size() == 0) {
        std::stringstream ss;
        ss << std::cin.rdbuf();
        input = ss.str();
    } else {
        std::ifstream in(file);
        if (!in.good()) throw std::invalid_argument("Unable to open file \"" + file + "\"");
        std::stringstream ss;
        ss << in.rdbuf();
        input = ss.str();
    }
    const auto lines = rg::split(input, '\n');
    try {
        for (unsigned i = 0; i < lines.size(); i += (fastq ? 4 : 2)) {
            vargas::SAM::Record rec;
            rec.query_name = std::string(lines.at(i).begin() + 1,
                                         std::find_if(lines.at(i).begin() + 1, lines.at(i).end(), isspace));
            rec.seq = lines.at(i + 1);
            if (fastq) rec.qual = lines.at(i + 3);
            ret.push(rec);
        }
    } catch (std::exception &e) {
        throw std::runtime_error("Invalid FASTA/Q file.");
    }
    ret.next();
}

void align_help(const cxxopts::Options &opts) {
    using std::cerr;
    using std::endl;

    cerr << opts.help() << "\n" << endl;
    cerr << "Elements per SIMD vector: " << vargas::Aligner::read_capacity() << endl;
}

ReadFmt read_fmt(const std::string filename) {
    std::ifstream in(filename);
    if (!in.good()) throw std::invalid_argument("Invalid read file: " + filename);
    std::string line;
    std::getline(in, line);
    if (line.size() && line.at(0) == '>') return ReadFmt::FASTA;
    std::getline(in, line); // Read or some SAM line
    std::getline(in, line); // either + or some SAM line
    if (line.at(0) == '+') return ReadFmt::FASTQ;
    else return ReadFmt::SAM;
}

TEST_SUITE("System");
TEST_CASE ("Load FASTQ") {
    std::string tmpfq = "tmp_fastq.va";
    {
        std::ofstream o(tmpfq);
        o << "@name desc\nAAAAACCCCC\n+\n!!!!!!!!!!";
    }
    CHECK(read_fmt(tmpfq) == ReadFmt::FASTQ);
    vargas::isam ss;
    load_fast(tmpfq, true, ss);

    CHECK(ss.record().query_name == "name");
    CHECK(ss.record().seq == "AAAAACCCCC");
    CHECK(ss.record().qual == "!!!!!!!!!!");
    CHECK_FALSE(ss.next());
    remove(tmpfq.c_str());
}
TEST_CASE ("Load FASTA") {
    std::string tmpfq = "tmp_fastq.va";
    {
        std::ofstream o(tmpfq);
        o << ">name desc\nAAAAACCCCC\n>x\nGGGGGTTTTT";
    }
    CHECK(read_fmt(tmpfq) == ReadFmt::FASTA);
    vargas::isam ss;
    load_fast(tmpfq, false, ss);

    CHECK(ss.record().query_name == "x");
    CHECK(ss.record().seq == "GGGGGTTTTT");
    ss.next();

    CHECK(ss.record().query_name == "name");
    CHECK(ss.record().seq == "AAAAACCCCC");
    CHECK_FALSE(ss.next());
    remove(tmpfq.c_str());
}
TEST_CASE ("Coordinate System Matches") {
    srand(12345);
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
    gb.set_region("x:0-50");
    vargas::Graph g = gb.build();

    vargas::Sim::Profile prof;
    prof.len = 5;
    vargas::Sim sim(g, prof);

    vargas::Aligner aligner(5);
    auto reads = sim.get_batch(aligner.read_capacity()*2);

    std::vector<std::string> seqs;
    std::vector<unsigned> targets;
    for (auto &r : reads) {
        seqs.push_back(r.seq);
        targets.push_back(r.pos + r.seq.length() - 1);
    }

    auto results = aligner.align(seqs, targets, g.begin(), g.end());

    int sum = 0;
    for (auto i : results.correct) if (i == 1) ++sum;
    // Sometimes an ambig read is made, so we just want most to be right
    CHECK(sum > 0.8*reads.size());

    remove(tmpfa.c_str());
    remove(tmpvcf.c_str());
    remove((tmpfa + ".fai").c_str());
}
TEST_SUITE_END();