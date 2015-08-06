//
// Created by gaddra on 8/6/15.
//

#include "main.h"

using std::string;
using std::cout;
using std::endl;
using std::cerr;

int main(int argc, char* argv[])
{
    /** Scores **/
    int32_t match = 2, mismatch = 2, gap_open = 3, gap_extension = 1;
    string VCF, REF;

    GetOpt::GetOpt_pp args(argc, argv);
    if(args >> GetOpt::OptionPresent('h', "help"))
    {
        cout << "VMatch options:" << endl;
        cout << "-v\t--vcf           VCF file, uncompressed" << endl;
        cout << "-r\t--ref           reference FASTA" << endl;
        cout << "-m\t--match         Match score, default  " << match << endl;
        cout << "-n\t--mismatch      Mismatch score, default " << mismatch << endl;
        cout << "-o\t--gap_open      Gap opening score, default " << gap_open << endl;
        cout << "-e\t--gap_extend    Gap extend score, default " << gap_extension << endl;
        exit(0);
    }
    if (!(args >> GetOpt::Option('v', "vcf", VCF))
        || !(args >> GetOpt::Option('r', "ref", REF)))
    {
        cerr << "VCF and/or ref file not specified, use -v and -r." << endl;
        exit(1);
    }
    args >> GetOpt::Option('m', "match", match)
         >> GetOpt::Option('n', "mismatch", mismatch)
         >> GetOpt::Option('o', "gap_open", gap_open)
         >> GetOpt::Option('e', "gap_extend", gap_extension);

    char *ref0 = "GGGCG", *ref1 = "TAAAAA", *ref2 = "ACGTGCA";
    char *query = "ACGT";

    int8_t *nt_table = gssw_create_nt_table(); // Nucleotide -> Num
    int8_t *mat = gssw_create_score_matrix(match, mismatch);

    gssw_node* nodes[3];
    nodes[0] = gssw_node_create(ref0, 0, ref0, nt_table, mat);
    nodes[1] = gssw_node_create(ref1, 1, ref1, nt_table, mat);
    nodes[2] = gssw_node_create(ref2, 2, ref1, nt_table, mat);

    gssw_nodes_add_edge(nodes[0], nodes[1]);
    gssw_nodes_add_edge(nodes[1], nodes[2]);

    gssw_graph* graph = gssw_graph_create(3);
    gssw_graph_add_node(graph, nodes[0]);
    gssw_graph_add_node(graph, nodes[1]);
    gssw_graph_add_node(graph, nodes[2]);

    gssw_graph_fill(graph, query, nt_table, mat, gap_open, gap_extension, 15, 2);

    printNode(graph->nodes[0]);
    printNode(graph->nodes[1]);
    printNode(graph->nodes[2]);
    std::cout << graph->max_node->id;


    gssw_graph_destroy(graph);
    delete[] nt_table;
    delete[] mat;

    return(0);
}
void printNode(gssw_node *node) {
    using namespace std;
    cout << "optimal score: " << node->alignment->score1 << " end: " << node->alignment->ref_end1 << endl;
    cout << "suboptimal score: " << node->alignment->score2 << " end: " << node->alignment->ref_end2 << endl << endl;
}