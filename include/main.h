/**
 * @author Ravi Gaddipati
 * @date June 26, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Interface for simulating and aligning reads from/to a DAG.
 *
 * @file
 */

#ifndef VARGAS_MAIN_H
#define VARGAS_MAIN_H


#include "alignment.h"
#include "graph.h"
#include "sim.h"

/**
 * Simulate reads from given graph definitions.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int sim_main(const int argc, const char *argv[]);

/**
 * Align given reads to specified target graphs.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int align_main(const int argc, const char *argv[]);

/**
 * Define graphs from a FASTA and a VCF/BCF file. Allows graphs
 * to remain consistent between simulating and aligning steps.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int define_main(const int argc, const char *argv[]);

/**
 * Profile aligner and graph construction.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int profile(const int argc, const char *argv[]);

/**
 * Split a SAM file into multiple files. All files have the same header.
 * @param argc CL arg count
 * @param argv CL args
 */
int split_main(const int argc, const char *argv[]);

/**
 * Merge multiple SAM files that have the same header.
 * @param argc CL arg count
 * @param argv CL args
 */
int merge_main(const int argc, const char *argv[]);

/**
 * Extract fields from a SAM file and export them to a CSV file.
 * @param argc CL arg count
 * @param argv CL args
 */
int sam2csv(const int argc, const char *argv[]);


/**
 * Aligns the given vector of reads to the given graph,
 * using the provided aligners.
 * @param label Label to prepend to output
 * @param subgraph Graph to align to
 * @param reads Reads to align
 * @param aligners Use the given aligners. Size should be equal to number of threads
 * @param out Stream to output result to
 * @param threads number of execution threads.
 */
void align_to_graph(std::string label,
                    const Vargas::Graph &subgraph,
                    const std::vector<std::string> &reads,
                    const std::vector<std::shared_ptr<Vargas::ByteAligner>> &aligners,
                    std::ostream &out,
                    unsigned int threads);

// Menus
void main_help();
void profile_help();
void align_help();
void sim_help();
void define_help();
void split_help();
void merge_help();
void sam2csv_help();

template<typename T>
inline double chrono_duration(const std::chrono::time_point<T> &start_time) {
    return std::chrono::duration_cast<std::chrono::duration<double>>
        (std::chrono::steady_clock::now() - start_time).count();
}

template<typename T>
inline double chrono_duration(const std::chrono::time_point<T> &start_time, const std::chrono::time_point<T> &end) {
    return std::chrono::duration_cast<std::chrono::duration<double>>
        (end - start_time).count();
}

#endif //VARGAS_MAIN_H

// Checks to see if coordinate systems between the simulator and aligner line up.
TEST_CASE ("Coordinate matches") {
    srand(1);
    Vargas::Graph::Node::_newID = 0;
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

    Vargas::GraphBuilder gb(tmpfa);
    gb.open_vcf(tmpvcf);
    gb.node_len(5);
    gb.region("x:0-50");
    Vargas::Graph g = gb.build();

    Vargas::Sim::Profile prof;
    prof.len = 5;
    Vargas::Sim sim(g, prof);

    Vargas::ByteAligner aligner(g.max_node_len(), 5);
    auto reads = sim.get_batch(aligner.read_capacity());

    std::vector<std::string> seqs;
    std::vector<uint32_t> targets;
    for (auto &r : reads) {
        seqs.push_back(r.seq);
        targets.push_back(r.pos + r.seq.length() - 1);
    }

    auto results = aligner.align(seqs, targets, g.begin(), g.end());

    for (auto i : results.cor_flag) CHECK ((int) i == 1);

    remove(tmpfa.c_str());
    remove(tmpvcf.c_str());
}

TEST_CASE ("Cor flag") {
    srand(1);
    Vargas::Graph::Node::_newID = 0;
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

    Vargas::GraphBuilder gb(tmpfa);
    gb.open_vcf(tmpvcf);
    gb.region("x:0-100");
    Vargas::Graph g = gb.build();

    Vargas::ByteAligner aligner(g.max_node_len(), 6);
    Vargas::isam reads(reads_file);

    std::vector<Vargas::SAM::Record> records;
    std::vector<std::string> read_seq;
    std::vector<uint32_t> targets;
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

        records[i].aux.set(ALIGN_SAM_COR_FLAG_TAG, (int) res.cor_flag[i]);
    }

    Vargas::osam align_out("tmp_aout.sam", reads.header());
    for (auto &r : records) align_out.add_record(r);

    remove(tmpfa.c_str());
    remove(tmpvcf.c_str());
    remove(reads_file.c_str());
    remove("tmp_aout.sam");
}