//
// Aligns reads using gssw, traceback is not possible.
// Created by Ravi Gaddipati (rgaddip1@jhu.edu).
//

#include "main.h"

using std::string;
using std::cout;
using std::endl;
using std::cerr;

bool print = true;
bool novar = false;

int main(int argc, char *argv[]) {
    /** Scores **/
    int32_t match = 2, mismatch = 2;
    uint8_t gap_open = 3, gap_extension = 1;

    /** read number, default region **/
    int32_t readnum = 0, minpos = 0, maxpos = 2147483640;
    std::vector<string> region_split(0);

    /** Default sim read error **/
    float muterr = 0.01f, indelerr = 0.01f;

    /** Read generation defaults **/
    int32_t numreads = 1000, readlen = 100, readEnd;
    bool simreads = false;

    /** File names and arguments **/
    string VCF = "", REF = "", outfile = "", buildfile = "", simfile = "",
            readfile = "", alignfile = "", query, read, region,
            NVbuildfile;

    /** File streams **/
    std::ifstream reads;
    std::ofstream aligns, NValigns, simout;

    /** Graph score and conversion table **/
    int8_t *nt_table = gssw_create_nt_table();
    int8_t *mat = gssw_create_score_matrix(match, mismatch);

    /** Alignment graph and max node length **/
    gssw_graph *graph, *NVgraph;
    int32_t maxNodelen = 50000;
    bool dual = false;

    srand(time(NULL));

    GetOpt::GetOpt_pp args(argc, argv);
    if (args >> GetOpt::OptionPresent('h', "help")) {
        cout << "---------------------------- VMatch, August 2015. rgaddip1@jhu.edu ----------------------------" << endl;
        cout << "-v\t--vcf           VCF file, uncompressed." << endl;
        cout << "-r\t--ref           reference single record FASTA" << endl;
        cout << "-b\t--buildfile     quick rebuild file, required if -v, -r are not defined. Takes priority." << endl;
        cout << "-B\t--NVbuildfile   quick rebuild file for no variant graph, use with -D." << endl;
        cout << "-g\t--maxlen        Maximum node length, default " << maxNodelen << endl;
        cout << "-m\t--match         Match score, default  " << match << endl;
        cout << "-n\t--mismatch      Mismatch score, default " << mismatch << endl;
        cout << "-o\t--gap_open      Gap opening score, default " << int32_t(gap_open) << endl;
        cout << "-e\t--gap_extend    Gap extend score, default " << int32_t(gap_extension) << endl;
        cout << "-t\t--outfile       Graph output file for quick rebuild" << endl;
        cout << "-d\t--reads         Reads to align, one per line. Reads can be appended with # and a comment." << endl;
        cout << "-s\t--string        Align a single string to stdout, overrides read file arguments" << endl;
        cout << "-a\t--aligns        Output file for alignments, one per line. With -D, '.nv' will be appended." << endl;
        cout << "-R\t--region        Ref region, inclusive: min:max. Default is entire graph." << endl;
        cout << "-p\t--noprint       Disable stdout printing" << endl;
        cout << "-x\t--novar         Generate and align to a no-variant graph. VCF still required." << endl;
        cout << "-D\t--dual          Align to both variant and non-variant graphs." << endl;
        cout << "-i\t--simreads      Simulate reads and write to the file specified by -T" << endl;
        cout << "-T\t--readout       Simulated reads output file, 1 per line followed by position information." << endl;
        cout << "-N\t--numreads      Number of reads to simulate, default " << numreads << endl;
        cout << "-M\t--muterr        Simulated read mutation error rate, default " << muterr << endl;
        cout << "-I\t--indelerr      Simulated read Indel error rate, default " << indelerr << endl;
        cout << "-L\t--readlen       Nominal read length, default " << readlen << endl << endl;

        cout << "Sim read format:    READ#NODE_ID, NODE_LEN, NODE_MAX_POSITION, READ_END_POSITION" << endl;
        cout << "Alignment format:   #READ; 1_NODE_ID, NODE_LEN, 1_NODE_MAX, 1_SCORE, 1_END_POS; " << endl;
        cout << "                    2_NODE_ID, 2_NODE_LEN, 2_NODE_MAX, 2_SCORE, 2_END_POS" << endl;
        cout << "Note: Suboptimal score only applies to nodes different then the best alignment " <<
                endl << "node. Control granularity with --maxlen." << endl << endl;
        exit(0);
    }

    /** make sure there's a valid input **/
    if (!(args >> GetOpt::Option('b', "buildfile", buildfile))) {
        if (!(args >> GetOpt::Option('v', "vcf", VCF))
            || !(args >> GetOpt::Option('r', "ref", REF))) {
            cerr << "No inputs specified, see options with -h" << endl;
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
    >> GetOpt::OptionPresent('p', "noprint", print)
    >> GetOpt::Option('R', "region", region)
    >> GetOpt::OptionPresent('x', "novar", novar)
    >> GetOpt::Option('N', "numreads", numreads)
    >> GetOpt::Option('M', "muterr", muterr)
    >> GetOpt::Option('I', "indelerr", indelerr)
    >> GetOpt::Option('L', "readlen", readlen)
    >> GetOpt::OptionPresent('i', "simreads", simreads)
    >> GetOpt::Option('T', "readout", simfile)
    >> GetOpt::Option('g', "maxlen", maxNodelen)
    >> GetOpt::OptionPresent('D', "dual", dual)
    >> GetOpt::Option('B', "NVbuildfile", NVbuildfile);

    print = !print;

    /** Parse region **/
    if (region.length() > 0) {
        split(region, ':', region_split);
        if (region_split.size() < 1) {
            std::cerr << "Malformed region, must be in the form a:b" << endl;
            exit(1);
        }
        minpos = std::atoi(region_split[0].c_str());
        maxpos = std::atoi(region_split[1].c_str());
    }

    /** Build graph **/
    if(dual) {
        if(buildfile.length() < 1 || NVbuildfile.length() < 1) {
            cerr << "-b and -B need to be defined with -D." << endl;
            exit(1);
        }
        graph = buildGraph(buildfile, nt_table, mat);
        NVgraph = buildGraph(NVbuildfile,nt_table, mat);
    } else {
        if (buildfile.length() > 0) graph = buildGraph(buildfile, nt_table, mat);
        else graph = generateGraph(REF, VCF, nt_table, mat, minpos, maxpos, maxNodelen, outfile);
    }

    /** Simulate reads **/
    if (simreads) {
        if(simfile.length() < 1) {
            cerr << "Specify an output file for simulated reads with -T" << endl;
            exit(1);
        }
        simout.open(simfile.c_str());
        if (print) cout << "Generating reads..." << endl;
        for(int i = 0; i < numreads; i++) {
            if (print) cout << std::setw(12) << i + 1 << '\r' << std::flush;
            simout << generateRead(*graph, readlen, muterr, indelerr) << endl;
        }
        if(print) cout << endl;
        simout.close();
    }

    /** Align to graph **/
    if (query.length() > 0) {
        /** If a single query is specified **/
        gssw_graph_fill(graph, query.c_str(), nt_table, mat, gap_open, gap_extension, query.length() / 2, 2);
        printNode(graph->max_node);
    }
    else {
        /** Read file **/
        reads.open(readfile.c_str());
        aligns.open(alignfile.c_str());
        if(dual) NValigns.open((alignfile + ".nv").c_str());
        if (!reads.good() || !aligns.good() || (!NValigns.good() && dual)) {
            cerr << "Error opening reads file or alignment files. No alignment will be done." << endl;
            exit(0);
        }
        aligns << "#READ; 1_NODE_ID, NODE_LEN, 1_NODE_MAX, 1_SCORE, 1_END_POS; " <<
                "2_NODE_ID, 2_NODE_LEN, 2_NODE_MAX, 2_SCORE, 2_END_POS" << endl;
        if (dual) NValigns << "#READ; 1_NODE_ID, NODE_LEN, 1_NODE_MAX, 1_SCORE, 1_END_POS; " <<
                    "2_NODE_ID, 2_NODE_LEN, 2_NODE_MAX, 2_SCORE, 2_END_POS" << endl;
        if (print) cout << "Aligning reads..." << endl;
        while (std::getline(reads, read)) {
            readnum++;
            if (print) cout << std::setw(12) << readnum << '\r' << std::flush;

            readEnd = read.find('#');
            if(readEnd == string::npos)  readEnd = read.length();

            gssw_graph_fill(graph, read.substr(0, readEnd).c_str(), nt_table, mat, gap_open, gap_extension, read.length() / 2, 2);
            aligns << read << ";"
            << graph->max_node->id << ","
            << graph->max_node->len << ","
            << graph->max_node->data << ","
            << graph->max_node->alignment->score1 << ","
            << graph->max_node->alignment->ref_end1 << ";"
            << graph->nodes[graph->max_node->alignment->prevmax]->id << ","
            << graph->nodes[graph->max_node->alignment->prevmax]->len << ","
            << graph->nodes[graph->max_node->alignment->prevmax]->data << ","
            << graph->nodes[graph->max_node->alignment->prevmax]->alignment->score1 << ","
            << graph->nodes[graph->max_node->alignment->prevmax]->alignment->ref_end1 << endl;

            /** Align to second graph if -D is specified **/
            if(dual) {
                gssw_graph_fill(NVgraph, read.substr(0, readEnd).c_str(), nt_table, mat, gap_open, gap_extension, read.length() / 2, 2);
                NValigns << read << ";"
                << NVgraph->max_node->id << ","
                << NVgraph->max_node->len << ","
                << NVgraph->max_node->data << ","
                << NVgraph->max_node->alignment->score1 << ","
                << NVgraph->max_node->alignment->ref_end1 << ";"
                << NVgraph->nodes[NVgraph->max_node->alignment->prevmax]->id << ","
                << NVgraph->nodes[NVgraph->max_node->alignment->prevmax]->len << ","
                << NVgraph->nodes[NVgraph->max_node->alignment->prevmax]->data << ","
                << NVgraph->nodes[NVgraph->max_node->alignment->prevmax]->alignment->score1 << ","
                << NVgraph->nodes[NVgraph->max_node->alignment->prevmax]->alignment->ref_end1 << endl;
            }
        }
        reads.close();
        aligns.close();
        if(dual) NValigns.close();
    }

    gssw_graph_destroy(graph);
    delete[] nt_table;
    delete[] mat;

    return (0);
}


std::string generateRead(gssw_graph &graph, int32_t readLen, float muterr, float indelerr) {
    gssw_node *node;
    int base, RAND;
    char mut;
    std::stringstream readmut;
    std::string read = "";

    /** initial random node and base **/
    node = graph.nodes[rand() % (graph.size - 1)];
    base = rand() % (node->len);

    for(int i = 0; i < readLen; i++) {
        read += node->seq[base];
        base++;
        /** Go to next random node **/
        if(base == node->len) {
            node = node->next[rand() % node->count_next];
            if(node->count_next == 0) break; // End of graph reached
            base = 0;
        }
    }

    /** Mutate string **/
    for(int i = 0; i < read.length(); i++) {
        RAND = rand() % 100000;
        mut = read.at(i);
        if (RAND < (100000 - (100000 * indelerr / 2))) { // represents del
            if (RAND > (100000 - (100000 * indelerr))) RAND = rand() % int32_t(100000 * muterr); // insert rand base
            /** Mutation **/
            if (RAND < (100000 * muterr) / 4) mut = 'A';
            else if (RAND < 2 * (100000 * muterr) / 4) mut = 'G';
            else if (RAND < 3 * (100000 * muterr) / 4) mut = 'C';
            else if (RAND < (100000 * muterr)) mut = 'T';
            readmut << mut;
        }
    }
    /** Append suffix recording read position **/
    readmut << "#" << node->id << "," << node->len << "," << node->data << "," << base;
    return readmut.str();
}


gssw_graph *buildGraph(std::string buildfile, int8_t *nt_table, int8_t *mat) {
    using namespace std;

    string line;
    ifstream graphDat(buildfile.c_str());
    vector<string> lineSplit(0);
    vector<gssw_node *> nodes(0);
    uint32_t curr = 0;

    /** Build nodes and edges from buildfile **/
    if(print) cout << "Generating Nodes (" << buildfile << ")" << "..." << endl;
    while (getline(graphDat, line)) {
        split(line, ',', lineSplit);
        switch (lineSplit.size()) {
            case 3: // New node
                curr = uint32_t(strtol(lineSplit[1].c_str(), NULL, 10));
                nodes.push_back(gssw_node_create(int(strtol(lineSplit[0].c_str(), NULL, 10)),
                                                 curr,
                                                 lineSplit[2].c_str(), nt_table, mat));
                if (print) cout << setw(12) << curr << '\r' << flush;
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

    if(print) cout << endl << "Building Graph..." << endl;

    /** Buffer node **/
    nodes.push_back(gssw_node_create(NULL, ++curr, "", nt_table, mat));
    gssw_nodes_add_edge(nodes.end()[-2], nodes.end()[-1]);

    /** Add nodes to graph **/
    gssw_graph *graph = gssw_graph_create(nodes.size());
    for (int n = 0; n < nodes.size(); n++) {
        gssw_graph_add_node(graph, nodes[n]);
    }

    graphDat.close();
    return graph;
}


gssw_graph *generateGraph(
        std::string REF, std::string VCF,
        int8_t *nt_table, int8_t *mat,
        int32_t minpos, int32_t maxpos, // Region to parse
        int32_t maxNodeLen,
        std::string outputFile) {

    using namespace std;

    /** reference line, variant line **/
    string fline, vline;
    /** Parsed lines, VCF header row **/
    vector<string> vline_split(0), valt_split(0), fline_split(0), header(0);
    /** Vector of all the nodes in the graph **/
    vector<gssw_node *> nodes(0);
    /** columns in VCF file, current ref pos **/
    int32_t pos = -1, ref = -1, alt = -1, rpos = 0;
    int32_t nodenum = 0;
    /** Pos from variant file **/
    int vpos;
    /** To track edges that need to be built **/
    int numalts = 0, numprev;
    /** strings that represent node contents **/
    string vref, valt, nodestring;
    int32_t nodelen = 0;
    char base;
    bool write = false;

    /** File stream **/
    ifstream variants(VCF.c_str(), ios_base::in | ios_base::binary);
    ifstream reference(REF.c_str());
    ofstream out;

    if (outputFile.size() > 0) {
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
    pos = int32_t(find(header.begin(), header.end(), "pos") - header.begin());
    ref = int32_t(find(header.begin(), header.end(), "ref") - header.begin());
    alt = int32_t(find(header.begin(), header.end(), "alt") - header.begin());

    /** Find the POS, REF, ALT cols **/
    if (pos < 0 || ref < 0 || alt < 0) {
        cerr << "POS, REF, and/or ALT not found in VCF header." << endl;
        exit(1);
    }

    /** Generate Nodes **/
    if (print) cout << "Generating nodes..." << endl;

    /** Go to minimum position **/
    if (minpos > 0) {
        while (rpos < minpos - 1) {
            reference.get(base);
            if (!isspace(base)) rpos++;
        }
    }

    /** Process variants **/
    while (getline(variants, vline)) {
        nodestring = "";
        nodelen = 0;
        split(vline, '\t', vline_split);
        vpos = atoi(vline_split[pos].c_str());
        if (vpos <= rpos) goto endvar;
        if (vpos > maxpos) break;
        vref = vline_split[ref];
        valt = vline_split[alt];
        if (print) cout << setw(12) << nodenum << '\r' << flush;

        /** build node string up to var pos **/
        while (rpos < vpos - 1) {
            if (!reference.get(base)) {
                cerr << "End of ref found while looking for variant pos " << vpos << endl;
                exit(1);
            }
            if (!isspace(base)) {
                nodestring += base;
                rpos++;
                nodelen++;
            }

            /** Max node length reached, split node **/
            if (nodelen == maxNodeLen) {
                nodes.push_back(gssw_node_create(rpos, nodenum, nodestring.c_str(), nt_table, mat));
                if (write) out << rpos << "," << nodenum << "," << nodestring.c_str() << endl;
#if debug > 4
            cout << "Node: " << rpos << ", ID: " << nodenum << ", " << nodestring << endl;
#endif
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
                numalts = 1;
                nodestring = "";
                nodelen = 0;
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
            if(isspace(base)) {
                reference.get(base);
            }
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
        if (!novar) {
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
    while ((rpos < maxpos || maxpos < 0) && reference.get(base)) {
        if (!isspace(base)) {
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
    if (print) cout << endl << nodes.size() << " nodes generated. Building graph..." << endl;

    /** Buffer node at the end, alignment doesn't seem to look at the last node. **/
    nodes.push_back(gssw_node_create(rpos, nodenum, "", nt_table, mat));
    gssw_nodes_add_edge(nodes.end()[-2], nodes.end()[-1]);
#if debug > 4
    cout << "Node: " << rpos << ", ID: " << nodenum << ", " << "" << endl;
    cout << "Edge: " << nodes.end()[-2]->id << ", " << nodes.end()[-1]->id << endl;
#endif
    if (nodes.size() > 4294967294) {
        cerr << "Too many nodes to generate graph." << endl;
        exit(1);
    }

    /** Add nodes to graph **/
    gssw_graph *graph = gssw_graph_create(nodes.size());
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
    /** Split string with delim, return a vector **/
    std::stringstream ss(s);
    std::string item;
    elems = *new std::vector<std::string>(0);

    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
