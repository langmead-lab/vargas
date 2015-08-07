//
// Created by gaddra on 8/6/15.
//

#include "main.h"
#include "gssw/src/gssw.h"

using std::string;
using std::cout;
using std::endl;
using std::cerr;

int main(int argc, char *argv[]) {
    /** Scores **/
    int32_t match = 2, mismatch = 2;
    uint8_t gap_open = 3, gap_extension = 1;
    string VCF, REF;

    GetOpt::GetOpt_pp args(argc, argv);
    if (args >> GetOpt::OptionPresent('h', "help")) {
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
        || !(args >> GetOpt::Option('r', "ref", REF))) {
        cerr << "VCF and/or ref file not specified, use -v and -r." << endl;
        exit(1);
    }
    args >> GetOpt::Option('m', "match", match)
    >> GetOpt::Option('n', "mismatch", mismatch)
    >> GetOpt::Option('o', "gap_open", gap_open)
    >> GetOpt::Option('e', "gap_extend", gap_extension);


    char *query = "CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTT";

    int8_t *nt_table = gssw_create_nt_table(); // Nucleotide -> Num
    int8_t *mat = gssw_create_score_matrix(match, mismatch);

    gssw_graph *graph = generateGraph(REF, VCF, nt_table, mat);
    gssw_graph_fill(graph, query, nt_table, mat, gap_open, gap_extension, 30, 2);
    printNode(graph->max_node);

    gssw_graph_destroy(graph);

/**    gssw_node* nodes[3];
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
 **/
    delete[] nt_table;
    delete[] mat;

    return (0);
}

/// <summary>
/// Generates a graph from the given reference and variant file.
/// </summary>
/// <param name="REF">Reference FASTA</param>
/// <param name="VCF">Variant File, uncompressed</param>
/// <param name="nt_table">base table, construct using gssw_create_nt_table()</param>
/// <param name="mat">Score matrix, use gssw_create_score_matrix()</param>
/// <returns>Constructed gssw_graph</returns>
gssw_graph* generateGraph(std::string REF, std::string VCF, int8_t *nt_table, int8_t *mat, std::string outputFile) {
    using namespace std;

    /** reference line, variant line **/
    string fline, vline;
    /** Parsed lines, VCF header row **/
    vector<string> vline_split(0), valt_split(0), fline_split(0), header(0);
    /** Vector of all the nodes in the graph **/
    vector<gssw_node *> nodes(0);
    /** columns in VCF file, current ref pos **/
    int pos = -1, ref = -1, alt = -1, rpos = 0;
    uint32_t nodenum = 0;
    /** Pos from variant file **/
    int vpos;
    /** To track edges that need to be built **/
    int numalts = 0, numprev;
    /** strings that represent node contents **/
    string vref, valt, nodestring;
    char base;
    bool write = false;

    ifstream variants(VCF.c_str(), ios_base::in | ios_base::binary);
    ifstream reference(REF.c_str());
    ofstream out;

    if(outputFile.size() > 0) {
        write = true;
        out.open(outputFile.c_str());
    }

    if (!variants.good() || !reference.good()) {
        boolalpha(cout);
        cerr << "Error in opening files." << endl;
        cerr << VCF << ": " << variants.good() << endl;
        cerr << REF << ": " << reference.good() << endl;
        exit(1);
    }

    getline(reference, fline);
    if (fline.at(0) != '>') cerr << "Error in ref file, first char should be >" << endl;

    /** Go to first VCF record **/
    do { getline(variants, vline); } while (vline.substr(0, 2) == "##");
    transform(vline.begin(), vline.end(), vline.begin(), ::tolower);
    split(vline, '\t', header);
    pos = int(find(header.begin(), header.end(), "pos") - header.begin());
    ref = int(find(header.begin(), header.end(), "ref") - header.begin());
    alt = int(find(header.begin(), header.end(), "alt") - header.begin());

    /** Find the POS, REF, ALT cols **/
    if (pos < 0 || ref < 0 || alt < 0) {
        cerr << "POS, REF, and/or ALT not found in VCF header." << endl;
        exit(1);
    }

    /** Generate Nodes **/
    cout << "Generating nodes..." << endl;
    while (getline(variants, vline)) {
        nodestring = "";
        split(vline, '\t', vline_split);
        vpos = atoi(vline_split[pos].c_str());
        vref = vline_split[ref];
        valt = vline_split[alt];

        /** build node string up to var pos **/
        while (rpos < vpos - 1) {
            if (!reference.get(base)) {
                cerr << "End of ref found while looking for variant pos " << vpos << endl;
                exit(1);
            }
            if (!isspace(base)) {
                nodestring += base;
                rpos++;
            }
        }

        /** If there is space between the variants, add a new node **/
        if (nodestring.length() > 0) {
            nodes.push_back(gssw_node_create(rpos, nodenum, nodestring.c_str(), nt_table, mat));
            if (write) out << rpos << "," << nodenum << "," << nodestring.c_str() << endl;
#if debug > 4
            cout << "Node: " << rpos << ", ID: " << nodenum << ", " << nodestring << endl;
#endif
            /** Only connect with edge if it's not the first node **/
            if (nodenum != 0) {
                /** Connect to all of the previous alt/ref nodes **/
                for (int i = 0; i < numalts; i++) {
                    gssw_nodes_add_edge(nodes.end()[-2 - i], nodes.end()[-1]);
                    if (write) out << nodes.end()[-2 - i] << "," << nodes.end()[-1] << endl;
#if debug > 4
                    cout << "Edge: " << nodes.end()[-2 - i]->id << ", " << nodes.end()[-1]->id << endl;
#endif
                }
            }
            nodenum++;
            numprev = 1;
        }
        else numprev = numalts;

        /** Ref node **/
        for (int i = 0; i < vref.length(); i++) {
            reference.get(base);
            rpos++;
        }
        nodes.push_back(gssw_node_create(rpos, nodenum, vref.c_str(), nt_table, mat));
        if (write) out << rpos << "," << nodenum << "," << vref.c_str() << endl;
#if debug > 4
        cout << "Node: " << rpos << ", ID: " << nodenum << ", " << vref << endl;
#endif
        nodenum++;
        numalts = 1;

        /** Variants **/
        split(valt, ',', valt_split);
        for (int i = 0; i < valt_split.size(); i++) {
            nodes.push_back(gssw_node_create(rpos, nodenum, valt_split[i].c_str(), nt_table, mat));
            if (write) out << rpos << "," << nodenum << "," << valt_split[i].c_str() << endl;
#if debug > 4
                cout << "Node: " << rpos << ", ID: " << nodenum << ", " << valt_split[i].c_str() << endl;
#endif
            nodenum++;
            numalts++;
        }

        /** Build edges **/
        for (int p = 0; p < numprev; p++) {
            for (int a = 0; a < numalts; a++) {
                gssw_nodes_add_edge(nodes.end()[-1 - numalts - p], nodes.end()[-1 - a]);
                if (write) out << nodes.end()[-1 - numalts - p] << "," << nodes.end()[-1 - a] << endl;
#if debug > 4
                cout << "Edge: " << nodes.end()[-1 - numalts - p]->id << ", " << nodes.end()[-1 - a]->id << endl;
#endif
            }
        }
    }

    cout << nodes.size() << " nodes generated. Building graph..." << endl;
    /** Buffer node at the end, alignment doesn't seem to look at the last node. **/
    nodes.push_back(gssw_node_create(NULL, nodenum, "", nt_table, mat));
    gssw_nodes_add_edge(nodes.end()[-2], nodes.end()[-1]);
    if (nodes.size() > 4294967294){
        cerr << "Too many nodes to generate graph." << endl;
        exit(1);
    }

    /** Add nodes to graph **/
    gssw_graph* graph = gssw_graph_create(nodes.size());
    for (int n = 0; n < nodes.size(); n++) {
        gssw_graph_add_node(graph, nodes[n]);
    }
    return graph;
}

/// <summary>
/// Prints node information.
/// </summary>
/// <param name="node">gssw_node to print</param>
void printNode(gssw_node *node) {
    using namespace std;
    cout << "Node sequence: " << node->seq << endl;
    cout << "Ending position: " << node->data << endl;
    cout << "optimal score: " << node->alignment->score1 << " end: " << node->alignment->ref_end1 << endl;
    cout << "suboptimal score: " << node->alignment->score2 << " end: " << node->alignment->ref_end2 << endl;
}


/// <summary>
/// Splits the specified string, resets elems and returns with split string.
/// </summary>
/// <param name="s">The string</param>
/// <param name="delim">The delimiter</param>
/// <param name="elems">Vector to store results in. Vector is replaced!</param>
/// <returns>Vector of split string.</returns>
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    elems = *new std::vector<std::string>(0);

    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}