//
// Created by gaddra on 8/6/15.
//

#include <iostream>
#include "gssw/src/gssw.h"
#include "main.h"

using std::string;

int main(int argc, char** argv)
{
    /** Scores **/
    int32_t match = 2, mismatch = 2, gap_open = 3, gap_extension = 1;

    char *ref = "GCCCGCGTAGCAGCATA";
    char *query = "CATAGCCAACG";

    int8_t *nt_table = gssw_create_nt_table(); // Nucleotide -> Num
    int8_t *mat = gssw_create_score_matrix(match, mismatch);

    gssw_node* nodes[2];
    nodes[0] = gssw_node_create(ref, 1, ref, nt_table, mat);
    nodes[1] = gssw_node_create(ref, 2, ref, nt_table, mat);

    gssw_nodes_add_edge(nodes[0], nodes[1]);

    gssw_graph* graph = gssw_graph_create(2);
    gssw_graph_add_node(graph, nodes[0]);
    gssw_graph_add_node(graph, nodes[1]);

    gssw_graph_fill(graph, query, nt_table, mat, gap_open, gap_extension, 15, 2);

    gssw_graph_mapping* gm = gssw_graph_trace_back (graph,
                                                    query,
                                                    strlen(query),
                                                    match,
                                                    mismatch,
                                                    gap_open,
                                                    gap_extension);

    gssw_print_graph_mapping(gm);
    gssw_graph_mapping_destroy(gm);
    // note that nodes which are referred to in this graph are destroyed as well
    gssw_graph_destroy(graph);

    free(nt_table);
    free(mat);

    return(0);
}