/**
 * Ravi Gaddipati
 * Dec 23, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Main aligner interface
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#include "align_main.h"
#include "alignment.h"
#include "sim.h"
#include "threadpool.h"
#include <mutex>

int align_main(int argc, char *argv[]) {
    std::string cl = "vargas ";
    {
        std::ostringstream ss;
        for (int i = 0; i < argc; ++i) ss << std::string(argv[i]) << " ";
        cl = ss.str();
    }

    // Load parameters
    unsigned match, npenalty, threads, chunk_size, subsample;
    std::string read_file, gdf, align_targets, out_file, pgid, mismatch, rdg, rfg;
    bool end_to_end = false, fwdonly = false, p64=false, msonly=false, maxonly=false, notraceback=false;

    cxxopts::Options opts("vargas align", "Align reads to a graph.");
    try {
        opts.add_options("Input")
        ("g,gdef", "<str> *Graph definition file.", cxxopts::value(gdf))
        ("U,reads", "<str> *Unpaired reads in SAM, FASTQ, or FASTA format.", cxxopts::value(read_file));

        opts.add_options("Optional")
        ("S,sam", "<str> Output file.", cxxopts::value(out_file))
        ("msonly", "Only report max score. Improves speed.", cxxopts::value(msonly)->implicit_value("1"))
        ("maxonly", "Only report max score, location, and count. Improves speed.", cxxopts::value(maxonly)->implicit_value("1"))
        ("phred64", "Qualities are Phred+64, not Phred+33.", cxxopts::value(p64)->implicit_value("1"))
        ("p,subsample", "<N> Sample N random reads, 0 for all.", cxxopts::value(subsample)->default_value("0"))
        ("a,alignto", "<str> Target graph, or SAM Read Group -> graph mapping.\"(RG:ID:<group>,<target_graph>;)+|<graph>\"", cxxopts::value(align_targets))
        ("s,assess", "[ID] Use score profile from a previous alignment.", cxxopts::value(pgid)->implicit_value("."))
        ("f,forward", "Only align to forward strand.", cxxopts::value(fwdonly))
        ("notraceback", "If graph contains no variants, do not compute traceback", cxxopts::value(notraceback)->implicit_value("1"));

        opts.add_options("Scoring")
        ("ete", "End to end alignment.", cxxopts::value(end_to_end))
        ("ma", "<N> Match bonus.", cxxopts::value(match)->default_value("2"))
        ("mp", "<MX,MN> Mismatch penalty. Lower qual=lower penalty.", cxxopts::value(mismatch)->default_value("6,2"))
        ("np", "<N> Penalty for non-A/C/G/T.", cxxopts::value(npenalty)->default_value("1"))
        ("rdg", "<N1,N2> Read gap open/extension penalty.", cxxopts::value(rdg)->default_value("1,3"))
        ("rfg", "<N1,N2> Ref gap open/extension penalty.", cxxopts::value(rfg)->default_value("1,3"));

        opts.add_options("Threading")
        ("j,threads", "<N> Number of threads.", cxxopts::value(threads)->default_value("1"))
        ("u,chunk", "<N> Partition into tasks of max size N.", cxxopts::value(chunk_size)->default_value("64"));

        opts.add_options()("h,help", "Display this message.");

        opts.parse(argc, argv);
    } catch (std::exception &e) {
        throw std::invalid_argument("Error parsing options: " + std::string(e.what()));
    }

    if (opts.count("h")) {
        align_help(opts);
        return 0;
    }

    if (!opts.count("gdef")) {
        align_help(opts);
        throw std::invalid_argument("Graph definition file required.");
    }

    if (!opts.count("reads")) {
        align_help(opts);
        throw std::invalid_argument("No read file provided.");
    }
    ReadFmt format = read_fmt(read_file);

    if (chunk_size < vargas::Aligner::read_capacity() || chunk_size % vargas::Aligner::read_capacity() != 0) {
        std::cerr << "[warn] Chunk size is not a multiple of SIMD vector length: "
                  << vargas::Aligner::read_capacity() << std::endl;
    }

    if (opts.count("assess") && format != ReadFmt::SAM) {
        throw std::invalid_argument("Assess is only available for SAM inputs.");
    }
    if (!align_targets.empty() && format != ReadFmt::SAM) {
        throw std::invalid_argument("Alignment targets only available for SAM inputs.");
    }

    if(opts.count("msonly") && opts.count("maxonly")) {
        throw std::invalid_argument("At most one of msonly and maxonly can be specified.");
    }

    vargas::isam reads;
    if (format == ReadFmt::FASTQ) {
        load_fast(read_file, true, reads, p64);
    } else if (format == ReadFmt::FASTA) {
        load_fast(read_file, false, reads, p64);
    } else {
        reads.open(read_file);
    }
    reads.subset(subsample);
    auto &reads_hdr = reads.header();

    vargas::ScoreProfile prof;
    {
        prof.match = match;
        prof.ambig = npenalty;

        auto sp = rg::split(mismatch, ',');
        if (sp.size() == 2) {
            prof.mismatch_max = std::stoi(sp[0]);
            prof.mismatch_min = std::stoi(sp[1]);
        }
        else if (sp.size() == 1) {
            prof.mismatch_max = prof.mismatch_min = std::stoi(sp[0]);
        }
        else {
            throw std::invalid_argument("Invalid --mp argument.");
        }

        sp = rg::split(rdg, ',');
        if (sp.size() == 2) {
            prof.read_gopen = std::stoi(sp[0]);
            prof.read_gext = std::stoi(sp[1]);
        }
        else {
            throw std::invalid_argument("Invalid --rdg argument.");
        }

        sp = rg::split(rfg, ',');
        if (sp.size() == 2) {
            prof.ref_gopen = std::stoi(sp[0]);
            prof.ref_gext = std::stoi(sp[1]);
        }
        else {
            throw std::invalid_argument("Invalid --rfg argument.");
        }
    }

    if (pgid == ".") {
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
        } catch (std::exception &e) {
            throw std::invalid_argument("Unrecognized PG ID: " + pgid);
        }
    } else {
        prof.end_to_end = end_to_end;
    }

    vargas::SAM::Header::Program pg;
    pg.command_line = cl;
    pg.name = "vargas_align";
    pg.id = "VA";
    pg.version = __DATE__;
    std::replace_if(pg.version.begin(), pg.version.end(), isspace, ' '); // rm tabs
    const auto assigned_pgid = reads_hdr.add(pg);

    size_t read_len;
    auto task_list = create_tasks(reads, align_targets, chunk_size, read_len);

    const size_t num_tasks = task_list.size();
    if (num_tasks < threads) {
        std::cerr << "[warn] Number of threads is greater than number of tasks. Try decreasing -u.\n";
    }

    threads = threads ? threads > task_list.size() ? task_list.size() : threads
                      : 1;

    int bias = 255 - (read_len * match);
    const bool use_wide = (bias < 0) or (end_to_end and (prof.ref_gopen + (prof.ref_gext * (read_len - 1)) > bias || read_len * prof.mismatch_max > bias));
    if (use_wide) {
        std::cerr << "Score range: " << read_len * match << " to -" << std::min(prof.ref_gopen + (prof.ref_gext * (read_len - 1)), read_len * prof.mismatch_max) <<
        ". Using 16-bit aligner (" << vargas::WordAligner::read_capacity() << " reads/vector).\n";
    }
    std::cerr << "Scoring profile: " << prof.to_string() << "\n";

    std::vector<std::unique_ptr<vargas::AlignerBase>> aligners(threads);
    for (size_t k = 0; k < threads; ++k) {
        aligners[k] = make_aligner(prof, read_len, use_wide, msonly, maxonly);
    }


    std::cerr << "\nLoading \"" << gdf << "\"...\n";
    auto start_time = std::chrono::steady_clock::now();
    vargas::GraphMan gm(gdf);
    if (gm.labels().size() != 1 && maxonly) {
        std::cerr << "[warn] With --maxonly, max score position and count may be incorrect because the genome is a graph." << std::endl;
    }
    if (gm.labels().size() != 1 && !maxonly && !msonly) {
        throw std::invalid_argument("Cannot calculate 2nd-max score when the genome is a graph. Use --msonly or --maxonly.");
    }
    std::cerr << rg::chrono_duration(start_time) << "s.\n";

    if (out_file.length()) std::cerr << "Writing to \"" << (out_file.empty() ? "stdout" : out_file) << "\".\n";
    reads_hdr.programs[assigned_pgid].aux.set(ALIGN_SAM_PG_GDF, gdf);
    vargas::osam aligns_out(out_file, reads_hdr);
    align(gm, task_list, aligns_out, aligners, fwdonly, msonly, maxonly, notraceback);

    return 0;
}

struct align_helper {
    vargas::GraphMan &gm;
    std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>> &task_list;
    vargas::osam &out;
    const std::vector<std::unique_ptr<vargas::AlignerBase>> &aligners;
    bool fwdonly, msonly, maxonly, notraceback;
    std::mutex &mut;
};

void align_helper_func(void *data, long index, int tid) {
    align_helper &help(*(align_helper *)data);
    auto &aligners = help.aligners;
    auto &out = help.out;
    auto &task_list = help.task_list;
    auto &gm = help.gm;
    auto fwdonly = help.fwdonly;
    auto msonly = help.msonly;
    auto maxonly = help.maxonly;
    auto notraceback = help.notraceback;

    const size_t num_reads = task_list.at(index).second.size();
    std::vector<std::string> read_seqs(num_reads);
    std::vector<std::vector<char>> quals(num_reads);

    for (size_t i = 0; i < num_reads; ++i) {
        const auto &r = task_list.at(index).second.at(i);
        read_seqs[i] = r.seq;
        if (r.qual.size() == r.seq.size()) {
            std::transform(r.qual.begin(),
                           r.qual.end(),
                           std::back_inserter(quals[i]),
                           [](char c){ return c - 33; });
        }
    }
    auto subgraph = gm.at(task_list.at(index).first);
    vargas::Results aligns;
    aligners[tid]->align_into(read_seqs, quals, subgraph->begin(), subgraph->end(), aligns, fwdonly);

    //If no variants (# nodes == # contigs) compute the alignment traceback
    bool not_graph = subgraph->node_map()->size() == gm.resolver()._contig_hdr_order.size();

    for (size_t j = 0; j < task_list.at(index).second.size(); ++j) {
        vargas::SAM::Record &rec = task_list.at(index).second.at(j);
        auto abs = gm.absolute_position(aligns.max_pos[j]);
        rec.aux.set("AS", aligns.max_score[j]);
        if (!msonly) {
            // Can only guess start position for end to end
            // if (aligns.profile.end_to_end) rec.pos = abs.second - rec.seq.size() + 1;

            rec.ref_name = abs.first;
            rec.flag.rev_complement = aligns.max_strand[j] == vargas::Strand::REV;
            if (rec.flag.rev_complement) {
                rg::reverse_complement_inplace(rec.seq);
                std::reverse(rec.qual.begin(), rec.qual.end());
            }
            rec.aux.set(ALIGN_SAM_MAX_POS_TAG, abs.second);
            rec.aux.set(ALIGN_SAM_MAX_COUNT_TAG, aligns.max_count[j]);

            if (not_graph & !notraceback) {
                int nodeID = gm.nodeID_from_contig(rec.ref_name);
                vargas::Graph::nodemap_t _node_map = *(subgraph->node_map());
                //TODO upper-bound the length of reference slice needed based on the score or scoring function
                int ref_len = 2*rec.seq.length() < abs.second ? 2*rec.seq.length() : abs.second ;
                int ref_start = abs.second-ref_len;
                auto ref_iter = std::next(_node_map[nodeID].begin(),ref_start); //TODO this is linear time I think?

                // Allocate the three DP score matrixes: M (match) D (deletion) I (insertion), initialize with zero
                std::vector<std::vector<int>> M;
                M.resize(rec.seq.length()+1);
                for (auto &&r : M) r.resize(ref_len+1,0);
                std::vector<std::vector<int>> D;
                D.resize(rec.seq.length()+1);
                for (auto &&r : D) r.resize(ref_len+1,0);
                std::vector<std::vector<int>> I;
                I.resize(rec.seq.length()+1);
                for (auto &&r : I) r.resize(ref_len+1,0);

                // Allocate the three DP traceback matrixes: tM (match) tD (deletion) tI (insertion), initialize with zero
                std::vector<std::vector<int>> tM;
                tM.resize(rec.seq.length()+1);
                for (auto &&r : tM) r.resize(ref_len+1,0);
                std::vector<std::vector<int>> tD;
                tD.resize(rec.seq.length()+1);
                for (auto &&r : tD) r.resize(ref_len+1,0);
                std::vector<std::vector<int>> tI;
                tI.resize(rec.seq.length()+1);
                for (auto &&r : tI) r.resize(ref_len+1,0);

                // TODO check initializations
                // initialize first row and first column of each score matrix as appropriate
                // gaps in beginning of reference (first row) ending in match doesn't make sense
                std::fill(M[0].begin(),M[0].end(),-aligns.profile.ref_gext*ref_len);
                M[0][0] = 0; //except zero characters of each is free
                // gaps in beginning of reference (first row) ending in gap in query doesn't make sense
                std::fill(I[0].begin(),I[0].end(),-aligns.profile.ref_gext*ref_len);
                // gaps in beginning of query (first col) ending in gap in reference doesn't make sense
                // gaps in beginning of query (first col) ending in match doesn't make sense
                for (unsigned row = 1; row <= rec.seq.length(); ++row) {
                    D[row][0] = -aligns.profile.ref_gext*ref_len;
                    M[row][0] = -aligns.profile.ref_gext*ref_len;
                }
                if (aligns.profile.end_to_end) { //semiglobal
                    // gaps in beginning of query (first col) ending in gap in query accumulate
                    for (unsigned row = 1; row <= rec.seq.length(); ++row) {
                        I[row][0] = 0 - row * aligns.profile.read_gext - aligns.profile.read_gopen;
                    }
                }

                // Fill the matrixes
                bool has_quality = rec.seq.size() == quals[j].size();
                for (unsigned col = 1; col <= ref_len; ++col) {
                    char ref_char = rg::num_to_base(*ref_iter);
                    for (unsigned row = 1; row <= rec.seq.length(); ++row) {
                        char query_char = rec.seq[row-1];
                        char query_qual = has_quality ? rec.qual[row-1] : 40;
                        // Compute the M matrix entry. Force a match or a mismatch between the last characters
                        // Traceback always goes "diagonally" but to which other matrix?
                        int possibleM = M[row-1][col-1];
                        int possibleD = D[row-1][col-1];
                        int possibleI = I[row-1][col-1];
                        if (ref_char == 'N' or query_char == 'N') { //ambiguous query and/or reference
                            possibleM -= aligns.profile.ambig;
                            possibleD -= aligns.profile.ambig;
                            possibleI -= aligns.profile.ambig;
                        }
                        else if (ref_char != query_char) { //mismatch
                            possibleM -= aligns.profile.penalty(query_qual);
                            possibleD -= aligns.profile.penalty(query_qual);
                            possibleI -= aligns.profile.penalty(query_qual);
                        } else { //match
                            possibleM += aligns.profile.match;
                            possibleD += aligns.profile.match;
                            possibleI += aligns.profile.match;
                        }
                        int best = std::max({possibleM, possibleD, possibleI});
                        if (aligns.profile.end_to_end or best > 0) { //local mode nothing happens if it's <= 0
                            if (possibleM == best) {
                                M[row][col] = possibleM;
                                tM[row][col] = 0;
                            } else if (possibleD == best) {
                                M[row][col] = possibleD;
                                tM[row][col] = 1;
                            } else {
                                M[row][col] = possibleI;
                                tM[row][col] = 2;
                            }
                        }

                        // Compute the D matrix entry. Force a gap in end of reference.
                        // Traceback always goes "left" but to which other matrix?
                        //TODO check read ref gap
                        possibleM = M[row][col-1] - aligns.profile.read_gopen - aligns.profile.read_gext;
                        possibleD = D[row][col-1] - aligns.profile.read_gext;
                        possibleI = I[row][col-1] - aligns.profile.read_gopen - aligns.profile.read_gext;
                        best = std::max({possibleM, possibleD, possibleI});
                        if (aligns.profile.end_to_end or best > 0) { //local mode nothing happens if it's <= 0
                            if (possibleM == best) {
                                D[row][col] = possibleM;
                                tD[row][col] = 0;
                            } else if (possibleD == best) {
                                D[row][col] = possibleD;
                                tD[row][col] = 1;
                            } else {
                                D[row][col] = possibleI;
                                tD[row][col] = 2;
                            }
                        }

                        // Compute the I entry. Force a gap in end of reference.
                        // Traceback always goes "up" but to which other matrix?
                        //TODO check read ref gap
                        possibleM = M[row-1][col] - aligns.profile.ref_gopen - aligns.profile.ref_gext;
                        possibleD = D[row-1][col] - aligns.profile.ref_gopen - aligns.profile.ref_gext;
                        possibleI = I[row-1][col] - aligns.profile.ref_gext;
                        best = std::max({possibleM, possibleD, possibleI});
                        if (aligns.profile.end_to_end or best > 0) { //local mode nothing happens if it's <= 0
                            if (possibleM == best) {
                                I[row][col] = possibleM;
                                tI[row][col] = 0;
                            } else if (possibleD == best) {
                                I[row][col] = possibleD;
                                tI[row][col] = 1;
                            } else {
                                I[row][col] = possibleI;
                                tI[row][col] = 2;
                            }
                        }
                    }
                    ref_iter++;
                }
                // Compute traceback: CIGAR string and start position
                std::vector<char> aln; //reverse order of operations
                if (aligns.profile.end_to_end) { //semiglobal
                    //best score is in last row and last column because we ended the reference at the max-scoring position
                    int currCol = ref_len;
                    int currRow = rec.seq.length();
                    int bestMatrix = 0;
                    int best = M[currRow][currCol];
                    if (D[currRow][currCol] > best) {
                        bestMatrix = 1;
                        best = D[currRow][currCol];
                    }
                    if (I[currRow][currCol] > best) {
                        bestMatrix = 2;
                        best = I[currRow][currCol];
                    }
                    if (best != aligns.max_score[j]) {
                        std::cerr << "[WARNING] " << rec.query_name << " DP optimal score " << best << " and SIMD optimal score " << aligns.max_score[j] << " not equal\n";
                    }
                    while(currRow > 0 and currCol > 0) {
                        if (bestMatrix == 0) {
                            aln.push_back('M');
                            bestMatrix = tM[currRow][currCol];
                            --currCol;
                            --currRow;
                        } else if (bestMatrix == 1) {
                            aln.push_back('D');
                            bestMatrix = tD[currRow][currCol];
                            --currCol;
                        } else {
                            aln.push_back('I');
                            bestMatrix = tI[currRow][currCol];
                            --currRow;
                        }
                    }
                    rec.pos = abs.second - ref_len + currCol;
                    for (int row = 0; row < currRow; ++row) {
                        aln.push_back('I'); //unaligned bases in beginning of query
                    }
                } else { //local
                    //best score is somewhere in last column because we ended the reference at the max-scoring position
                    int bestRow = -1;
                    int best = -1;
                    int bestMatrix = -1;
                    for (int row = 0; row <= rec.seq.length(); ++row) {
                        if (M[row][ref_len] > best) {
                            best = M[row][ref_len];
                            bestRow = row;
                            bestMatrix = 0;
                        }
                        if (D[row][ref_len] > best) {
                            best = D[row][ref_len];
                            bestRow = row;
                            bestMatrix = 1;
                        }
                        if (I[row][ref_len] > best) {
                            best = I[row][ref_len];
                            bestRow = row;
                            bestMatrix = 2;
                        }
                    }
                    if (best != aligns.max_score[j]) {
                        std::cerr << "[WARNING] " << rec.query_name << " DP optimal score " << best << " and SIMD optimal score " << aligns.max_score[j] << " not equal\n";
                    }
                    int currCol = ref_len;
                    int currRow = bestRow;
                    for (int row = bestRow; row < rec.seq.length(); ++row) {
                        aln.push_back('S'); //unaligned bases in end of query
                    }
                    while(currRow > 0 and currCol > 0) {
                        if (bestMatrix == 0 and M[currRow][currCol] > 0) {
                            aln.push_back('M');
                            bestMatrix = tM[currRow][currCol];
                            --currCol;
                            --currRow;
                        } else if (bestMatrix == 1 and D[currRow][currCol] > 0) {
                            aln.push_back('D');
                            bestMatrix = tD[currRow][currCol];
                            --currCol;
                        } else if (I[currRow][currCol] > 0){
                            aln.push_back('I');
                            bestMatrix = tI[currRow][currCol];
                            --currRow;
                        } else { break; } //if score goes to or below zero
                    }
                    rec.pos = abs.second - ref_len + currCol;
                    for (int row = 0; row < currRow; ++row) {
                        aln.push_back('S'); //unaligned bases in beginning of query
                    }

                }
                // Reverse and run-length-collapse the sequence of operations
                char last_seen = aln.back();
                unsigned count = 1;
                std::string cigar;
                auto rit = std::next(aln.rbegin(),1);
                for (; rit != aln.rend(); ++rit) {
                    if (*rit != last_seen) {
                        cigar.append(std::to_string(count));
                        cigar.push_back(last_seen);
                        count = 1;
                        last_seen = *rit;
                    } else {count += 1;}
                }
                cigar.append(std::to_string(count));
                cigar.push_back(last_seen);
                rec.cigar = cigar;
            }

            // Flags for 2nd max
            if (!maxonly) {
                abs = gm.absolute_position(aligns.sub_pos[j]);
                rec.aux.set(ALIGN_SAM_SUB_SEQ, abs.first);
                rec.aux.set(ALIGN_SAM_SUB_POS_TAG, abs.second);
                rec.aux.set(ALIGN_SAM_SUB_COUNT_TAG, aligns.sub_count[j]);
                rec.aux.set(ALIGN_SAM_SUB_SCORE_TAG, aligns.sub_score[j]);
                rec.aux.set(ALIGN_SAM_SUB_STRAND_TAG, aligns.sub_strand[j] == vargas::Strand::FWD ? "fwd" : "rev");
            }
        }
    }

    {
        std::lock_guard<std::mutex> lock(help.mut);
        for (const auto & j : task_list.at(index).second) out.add_record(j);
    }
}

#if !NDEBUG
#else
#define at operator[]
#endif

void align(vargas::GraphMan &gm,
           std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>> &task_list,
           vargas::osam &out,
           const std::vector<std::unique_ptr<vargas::AlignerBase>> &aligners,
           bool fwdonly, bool msonly, bool maxonly, bool notraceback) {
    std::cerr << "Aligning... " << std::flush;
    rg::ForPool fp(aligners.size());
    auto start_time = std::chrono::steady_clock::now();

    const auto num_tasks = task_list.size();
    std::mutex mut;
    align_helper help{gm, task_list, out, aligners, fwdonly, msonly, maxonly, notraceback, mut};
    fp.forpool(&align_helper_func, (void *)&help, num_tasks);

    std::cerr << rg::chrono_duration(start_time) << "s.\n";

}

std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>>
create_tasks(vargas::isam &reads, std::string &align_targets, const int chunk_size, size_t &read_len) {
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

    if (alignment_pairs.empty()) {
        for (const auto &p : read_groups) {
            alignment_pairs.push_back("RG:ID:" + p.first + ",base");
        }
    }

    // Specify a single target graph
    if (alignment_pairs.size() == 1) {
        auto target = rg::split(alignment_pairs[0], ',');
        if (target.size() == 1) {
            alignment_pairs.clear();
            for (const auto &p : read_groups) {
                alignment_pairs.push_back("RG:ID:" + p.first + "," + target[0]);
            }
        }
    }

    // Maps target graph to read group ID's
    std::unordered_map<std::string, std::vector<std::string>> alignment_rg_map;

    std::string tag, val, target_val;
    for (const std::string &p : alignment_pairs) {
        auto pair = rg::split(p, ',');
        if (pair.size() != 2)
            throw std::invalid_argument("Malformed alignment pair \"" + p + "\".");
        if (pair[0].substr(0, 2) != "RG")
            throw std::invalid_argument(R"(Expected a read group tag 'RG:xx:', got ")" + pair[0] + "\"");
        if (pair[0].at(2) != ':')
            throw std::invalid_argument("Expected source format Read_group_tag:value in \"" + pair[0] + "\".");


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
                if (safe_beg != safe_end) {
                    task_list.emplace_back(sub_rg_pair.first, std::vector<vargas::SAM::Record>(safe_beg, safe_end));
                }
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
make_aligner(const vargas::ScoreProfile &prof, size_t read_len, bool use_wide, bool msonly, bool maxonly) {
    if (msonly) {
        if (prof.end_to_end) {
            if (use_wide) return rg::make_unique<vargas::MSWordAlignerETE>(read_len, prof);
            else return rg::make_unique<vargas::MSAlignerETE>(read_len, prof);
        } else {
            if (use_wide) return rg::make_unique<vargas::MSWordAligner>(read_len, prof);
            else return rg::make_unique<vargas::MSAligner>(read_len, prof);
        }
    }
    else {
        if (prof.end_to_end) {
            if (use_wide) return rg::make_unique<vargas::WordAlignerETE>(read_len, prof);
            else return rg::make_unique<vargas::AlignerETE>(read_len, prof);
        } else {
            if (use_wide) return rg::make_unique<vargas::WordAligner>(read_len, prof);
            else return rg::make_unique<vargas::Aligner>(read_len, prof);
        }
    }
}

void load_fast(std::string &file, const bool fastq, vargas::isam &ret, bool p64) {
    std::string input;
    if (file.empty()) {
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
            if (p64) std::transform(rec.qual.begin(), rec.qual.end(), rec.qual.begin(), [](char c){return c-31;});

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

    cerr << opts.help(opts.groups()) << "\n" << endl;
    cerr << "Elements per SIMD vector: " << vargas::Aligner::read_capacity() << endl;
}

ReadFmt read_fmt(const std::string& filename) {
    std::ifstream in(filename);
    if (!in.good()) throw std::invalid_argument("Invalid read file: " + filename);

    std::string line;
    if (!std::getline(in, line)) throw std::invalid_argument("Empty Read File."); // @SAM or fasta/q name
    if (line.substr(0,3) == "@HD") return ReadFmt::SAM;
    if (!std::getline(in, line)) throw std::invalid_argument("Invalid Read File."); // SAM comment/header, or read
    if (!std::getline(in, line)) return ReadFmt::FASTA; // Single record fasta, or SAM line, or +, or name
    if (!line.empty() && line[0] == '+') return ReadFmt::FASTQ;
    if (!line.empty() && (line[0] == '>' || line[0] == '@')) return ReadFmt::FASTA;
    return ReadFmt::SAM;
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
