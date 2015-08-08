//
// Aligns reads using gssw, traceback is not possible.
// Created by Ravi Gaddipati (rgaddip1@jhu.edu).
//

#include "main.h"

using std::string;
using std::cout;
using std::endl;
using std::cerr;

int main(int argc, char *argv[]) {
    /** Scores **/
    int32_t match = 2, mismatch = 2, readnum = 0, minpos = 0, maxpos = 2147483640;
    uint8_t gap_open = 3, gap_extension = 1;
    string VCF = "", REF = "", outfile = "", buildfile = "",
           readfile = "", alignfile = "", query, read, region;
    std::ifstream reads;
    std::ofstream aligns;
    bool print = true;

    GetOpt::GetOpt_pp args(argc, argv);
    if (args >> GetOpt::OptionPresent('h', "help")) {
        cout << "VMatch options:" << endl;
        cout << "-v\t--vcf           VCF file, uncompressed" << endl;
        cout << "-r\t--ref           reference FASTA" << endl;
        cout << "-b\t--buildfile     quick rebuild file, required if -v, -r are not defined. Takes priority." << endl;
        cout << "-m\t--match         Match score, default  " << match << endl;
        cout << "-n\t--mismatch      Mismatch score, default " << mismatch << endl;
        cout << "-o\t--gap_open      Gap opening score, default " << gap_open << endl;
        cout << "-e\t--gap_extend    Gap extend score, default " << gap_extension << endl;
        cout << "-t\t--outfile       Output file for quick rebuild of graph" << endl;
        cout << "-d\t--reads         Reads to align, one per line" << endl;
        cout << "-a\t--aligns        Outputfile for alignments" << endl;
        cout << "-s\t--string        Align a single string to stdout, overrides reads file" << endl;
        cout << "-p\t--noprint       0 to disable stdout printing" << endl;
        cout << "-R\t--region        Ref region, inclusive. min:max" << endl;
        exit(0);
    }

    if(!(args >> GetOpt::Option('b', "buildfile", buildfile))) {
        if (!(args >> GetOpt::Option('v', "vcf", VCF))
            || !(args >> GetOpt::Option('r', "ref", REF))) {
            cerr << "VCF and/or ref file not specified, use -v and -r or quick build with -b." << endl;
            exit(1);
        }
    }

    args >> GetOpt::Option('m', "match", match)
    >> GetOpt::Option('n', "mismatch", mismatch)
    >> GetOpt::Option('o', "gap_open", gap_open)
    >> GetOpt::Option('e', "gap_extend", gap_extension)
    >> GetOpt::Option('t', "outfile", outfile)
    >> GetOpt::Option('d', "reads", readfile)
    >> GetOpt::Option('a', "aligns", alignfile)
    >> GetOpt::Option('s', "string", query)
    >> GetOpt::Option('p', "noprint", print)
    >> GetOpt::Option('R', "region", region);


    int8_t *nt_table = gssw_create_nt_table(); // Nucleotide -> Num
    int8_t *mat = gssw_create_score_matrix(match, mismatch);
    gssw_graph *graph;
    std::vector<string> region_split(0);

    /** Parse region **/
    if(region.length() > 0) {
        split(region, ':', region_split);
        if (region_split.size() < 1) {
            std::cerr << "Malformed region, must be in the form a:b" << endl;
            exit(1);
        }
        minpos = std::atoi(region_split[0].c_str());
        maxpos = std::atoi(region_split[1].c_str());
    }

    if (buildfile.length() > 0) graph = buildGraph(buildfile, nt_table, mat);
    else graph = generateGraph(REF, VCF, nt_table, mat, minpos, maxpos, outfile);

    if(query.length() > 0){
        gssw_graph_fill(graph, query.c_str(), nt_table, mat, gap_open, gap_extension, query.length()/2, 2);
        printNode(graph->max_node);
    }
    else {
        reads.open(readfile.c_str());
        aligns.open(alignfile.c_str());
        if (!reads.good() || !aligns.good()) {
            cerr << "Error in opening files, specify -d and -a, or -s to align a string." << endl;
            exit(1);
        }
        aligns << "#node ID, node max position, score1, alignment end, score 2, alignment 2 end" << endl;
        while(std::getline(reads, read)){
            readnum++;
            cout << std::setw(12) << readnum << '\r' << std::flush;
            gssw_graph_fill(graph, read.c_str(), nt_table, mat, gap_open, gap_extension, read.length()/2, 2);
            aligns << ">" << read << ","
                   << graph->max_node->id << ","
                   << graph->max_node->data << ","
                   << graph->max_node->alignment->score1 << ","
                   << graph->max_node->alignment->ref_end1 << ","
                   << graph->max_node->alignment->score2 << ","
                   << graph->max_node->alignment->ref_end2 << endl;
        }
        reads.close();
        aligns.close();
    }

    gssw_graph_destroy(graph);
    delete[] nt_table;
    delete[] mat;

    return (0);
}


gssw_graph* buildGraph(std::string buildfile, int8_t *nt_table, int8_t *mat) {
    using namespace std;

    string line;
    ifstream graphDat(buildfile.c_str());
    vector<string> lineSplit(0);
    vector<gssw_node*> nodes(0);
    uint32_t curr = 0;

    while(getline(graphDat, line)) {
        split(line, ',', lineSplit);
        switch (lineSplit.size()) {
            case 3: // New node
                curr = uint32_t(strtol(lineSplit[1].c_str(), NULL, 10));
                nodes.push_back(gssw_node_create(int(strtol(lineSplit[0].c_str(), NULL, 10)),
                                                 curr,
                                                 lineSplit[2].c_str(), nt_table, mat));

                break;
            case 2: // New edge
                gssw_nodes_add_edge(nodes.end()[strtol(lineSplit[0].c_str(), NULL, 10)],
                                    nodes.end()[strtol(lineSplit[1].c_str(), NULL, 10)]);
                break;
            default:
                cerr << "Unexpected line in buildfile: " << endl << line << endl;
                break;
        }
    }

    nodes.push_back(gssw_node_create(NULL, ++curr, "", nt_table, mat));
    gssw_nodes_add_edge(nodes.end()[-2], nodes.end()[-1]);

    /** Add nodes to graph **/
    gssw_graph* graph = gssw_graph_create(nodes.size());
    for (int n = 0; n < nodes.size(); n++) {
        gssw_graph_add_node(graph, nodes[n]);
    }

    graphDat.close();
    return graph;
}


gssw_graph* generateGraph(std::string REF, std::string VCF, int8_t *nt_table, int8_t *mat, int32_t minpos, int32_t maxpos, std::string outputFile) {
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

    // Go to minimum position
    if(minpos > 0) {
        while (rpos < minpos - 1) {
            reference.get(base);
            if (!isspace(base)) rpos++;
        }
    }

    while (getline(variants, vline)) {
        nodestring = "";
        split(vline, '\t', vline_split);
        vpos = atoi(vline_split[pos].c_str());
        if (vpos <= rpos) goto endvar;
        if (vpos > maxpos) break;
        vref = vline_split[ref];
        valt = vline_split[alt];
        cout << setw(12) << nodenum << '\r' << flush;

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
                    if (write) out << -2 - i << "," << -1 << endl;
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
                if (write) out << -1 - numalts - p << "," << -1 - a << endl;
#if debug > 4
                cout << "Edge: " << nodes.end()[-1 - numalts - p]->id << ", " << nodes.end()[-1 - a]->id << endl;
#endif
            }
        }
        endvar:;
    }

    /** The remaining bases after the last variant **/
    nodestring = "";
    while((rpos < maxpos || maxpos < 0) && reference.get(base)) {
        if(!isspace(base)) {
            nodestring += base;
            rpos++;
        }
    }
    nodes.push_back(gssw_node_create(rpos, nodenum, nodestring.c_str(), nt_table, mat));
    if (write) out << rpos << "," << nodenum << "," << nodestring.c_str() << endl;
#if debug > 4
    cout << "Node: " << rpos << ", ID: " << nodenum << ", " << nodestring.c_str() << endl;
#endif
    nodenum++;
    for (int p = 0; p < numalts; p++) {
        gssw_nodes_add_edge(nodes.end()[-2 - p], nodes.end()[-1]);
        if (write) out << -2 - p << "," << -1 << endl;
#if debug > 4
        cout << "Edge: " << nodes.end()[-2 - p]->id << ", " << nodes.end()[-1]->id << endl;
#endif
    }
    cout << endl << nodes.size() << " nodes generated. Building graph..." << endl;

    /** Buffer node at the end, alignment doesn't seem to look at the last node. **/
    nodes.push_back(gssw_node_create(rpos, nodenum, "", nt_table, mat));
    gssw_nodes_add_edge(nodes.end()[-2], nodes.end()[-1]);
#if debug > 4
    cout << "Node: " << rpos << ", ID: " << nodenum << ", " << "" << endl;
    cout << "Edge: " << nodes.end()[-2]->id << ", " << nodes.end()[-1]->id << endl;
#endif
    if (nodes.size() > 4294967294){
        cerr << "Too many nodes to generate graph." << endl;
        exit(1);
    }

    /** Add nodes to graph **/
    gssw_graph* graph = gssw_graph_create(nodes.size());
    for (int n = 0; n < nodes.size(); n++) {
        gssw_graph_add_node(graph, nodes[n]);
    }

    variants.close();
    reference.close();
    out.close();
    return graph;
}


void printNode(gssw_node *node) {
    using namespace std;
    cout << "Node sequence: " << node->seq << endl;
    cout << "Node end position: " << node->data << endl;
    cout << "optimal score: " << node->alignment->score1 << " end: " << node->alignment->ref_end1 << endl;
    cout << "suboptimal score: " << node->alignment->score2 << " end: " << node->alignment->ref_end2 << endl;
}


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    elems = *new std::vector<std::string>(0);

    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
