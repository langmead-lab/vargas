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

int split_main(const int argc, const char *argv[]);

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

void derive_and_sim(std::string rg,
                    const Vargas::Graph &base,
                    const Vargas::Graph::Population &pop,
                    const Vargas::Sim::Profile &prof,
                    unsigned int num_reads,
                    std::vector<Vargas::SAM::Record> &results);


// Menus
void main_help();
void profile_help();
void align_help();
void sim_help();
void define_help();
void split_help();
void sam2csv_help();
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

    Vargas::GraphBuilder gb(tmpfa, tmpvcf);
    gb.node_len(5);
    gb.ingroup(100);
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
