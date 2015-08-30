//
// Aligns reads using gssw, traceback is not possible.
// Created by Ravi Gaddipati (rgaddip1@jhu.edu).
//

#include "main.h"
#include "gssw/src/gssw.h"

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
  int32_t regionMin = 0, regionMax = 2147483640;
  std::vector<string> region_split(0);

  /** Default sim read error **/
  float muterr = 0.01f, indelerr = 0.01f;

  /** Read generation defaults **/
  int32_t numreads = 1000, readlen = 100, tol = 5;

  /** File names and arguments **/
  string VCF = "", REF = "", outfile = "", buildfile = "", simfile = "",
      readfile = "", alignfile = "", query, region,
      NVbuildfile, line;

  /** File streams **/
  std::ifstream alignsin, alignsinNV;
  std::ofstream simout;

  /** Graph score and conversion table **/
  int8_t *nt_table = gssw_create_nt_table();
  int8_t *mat = gssw_create_score_matrix(match, mismatch);

  /** Alignment graph and max node length **/
  gssw_graph *graph, *NVgraph;
  int32_t maxNodelen = 50000;
  int32_t ingroup = -1;

  /** Modes of operation **/
  bool dual = false;
  bool statmode = false;
  bool simreads = false;
  bool dot = false;

  srand(uint32_t(time(NULL)));

  GetOpt::GetOpt_pp args(argc, argv);
  if (args >> GetOpt::OptionPresent('h', "help")) {
    cout << "---------------------------- VMatch, August 2015. rgaddip1@jhu.edu ----------------------------" <<
        endl;
    cout << "-v\t--vcf           VCF file, uncompressed." << endl;
    cout << "-r\t--ref           reference single record FASTA" << endl;
    cout << "-b\t--buildfile     quick rebuild file, required if -v, -r are not defined." << endl;
    cout << "-B\t--NVbuildfile   quick rebuild file for no variant graph, use with -D." << endl;
    cout << "-G\t--ingroup       Percent of individuals to build graph from, default all." << endl;
    cout << "-g\t--maxlen        Maximum node length, default " << maxNodelen << endl;
    cout << "-m\t--match         Match score, default  " << match << endl;
    cout << "-n\t--mismatch      Mismatch score, default " << mismatch << endl;
    cout << "-o\t--gap_open      Gap opening score, default " << int32_t(gap_open) << endl;
    cout << "-e\t--gap_extend    Gap extend score, default " << int32_t(gap_extension) << endl;
    cout << "-t\t--outfile       Graph output file for quick rebuild" << endl;
    cout << "-d\t--reads         Reads to align, one per line. Symbols after '#' are ignored." << endl;
    cout << "-s\t--string        Align a single string to stdout, overrides read file arguments" << endl;
    cout << "-a\t--aligns        Output alignments, one per line. With -D, '.nv' will be appended." << endl;
    cout << "-R\t--region        Ref region, inclusive: min:max. Default is entire graph." << endl;
    cout << "-p\t--noprint       Disable stdout printing" << endl;
    cout << "-x\t--novar         Generate and align to a no-variant graph. VCF still required." << endl;
    cout << "-D\t--dual          Align to both variant and non-variant graphs." << endl;
    cout << "-i\t--simreads      Simulate reads and write to the file specified by -T" << endl;
    cout << "-T\t--readout       Simulated reads output, 1 per line with position information." << endl;
    cout << "-N\t--numreads      Number of reads to simulate, default " << numreads << endl;
    cout << "-M\t--muterr        Simulated read mutation error rate, default " << muterr << endl;
    cout << "-I\t--indelerr      Simulated read Indel error rate, default " << indelerr << endl;
    cout << "-L\t--readlen       Nominal read length, default " << readlen << endl << endl;

    cout << "-S\t--stat          Alignment stats of file(s), aligns1[,aligns2,...]" << endl;
    cout << "-X\t--tol           Alignment within x bases to be as correct, default " << tol << endl << endl;

    cout << "-P\t--dot           Output graph in dot format." << endl << endl;

    cout << "Sim read format:    READ#NODE_ID, NODE_LEN, NODE_MAX_POSITION, READ_END_POSITION" << endl;
    cout << "Alignment format:   #READ; OPTIMAL_SCORE, OPTIMAL_END_POS, NUM_OPTIMAL_POSITIONS;" << endl;
    cout << "                    SUBOPTIMAL_SCORE, SUBOPTIMAL_END_POS, NUM_SUBOPTIMAL_POSITIONS" << endl;
    cout << "Note: Suboptimal score only applies to nodes different then the best alignment " <<
        endl << "node. Control granularity with --maxlen." << endl << endl;

    exit(0);
  }

  /** make sure there's a valid input **/
  if (!(args >> GetOpt::Option('b', "buildfile", buildfile)) && !(args >> GetOpt::Option('S', "stat", alignfile))) {
    if (!(args >> GetOpt::Option('v', "vcf", VCF))
        || !(args >> GetOpt::Option('r', "ref", REF))) {
      cerr << "No inputs specified, see options with -h" << endl;
      exit(1);
    }
  }

  if (alignfile.length() > 0) statmode = true;

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
      >> GetOpt::Option('B', "NVbuildfile", NVbuildfile)
      >> GetOpt::Option('X', "tol", tol)
      >> GetOpt::OptionPresent('P', "dot", dot)
      >> GetOpt::Option('G', "ingroup", ingroup);

  print = !print;

  /** Stat mode **/
  if (statmode) {
    stat_main(alignfile, tol);
    exit(0);
  }

  /** Parse region **/
  if (region.length() > 0) {
    split(region, ':', region_split);
    if (region_split.size() < 1) {
      std::cerr << "Malformed region, must be in the form a:b" << endl;
      exit(1);
    }
    regionMin = std::atoi(region_split[0].c_str());
    regionMax = std::atoi(region_split[1].c_str());
  }

  /** Build graph **/
  if (dual) {
    if (buildfile.length() < 1 || NVbuildfile.length() < 1) {
      cerr << "-b and -B need to be defined with -D." << endl;
      exit(1);
    }
    graph = buildGraph(buildfile, nt_table, mat);
    NVgraph = buildGraph(NVbuildfile, nt_table, mat);
  } else {
    if (buildfile.length() > 0) graph = buildGraph(buildfile, nt_table, mat);
    else graph = generateGraph(REF, VCF, nt_table, mat, regionMin, regionMax, maxNodelen, outfile, ingroup);
  }

  /** Output to dot format **/
  if (dot) graphToDot(graph, "graph.dot");
  if (dot && dual) graphToDot(NVgraph, "NVgraph.dot");

  /** Simulate reads **/
  if (simreads) {
    if (simfile.length() < 1) {
      cerr << "Specify an output file for simulated reads with -T" << endl;
      exit(1);
    }
    simout.open(simfile.c_str());
    if (print) cout << "Generating reads..." << endl;
    for (int i = 0; i < numreads; i++) {
      if (print) cout << std::setw(12) << i + 1 << '\r' << std::flush;
      simout << generateRead(*graph, readlen, muterr, indelerr) << endl;
    }
    if (print) cout << endl;
    simout.close();
  }


  /** Align to graph **/
  if (query.length() > 0) {
    /** If a single query is specified **/
    gssw_graph_fill(graph, query.c_str(), nt_table, mat, gap_open, gap_extension, uint32_t(query.length()), 2);
    printNode(graph->max_node);
    printNode(graph->submax_node);
  } else {
    align_main(graph, NVgraph, mat, nt_table, gap_open, gap_extension, readfile, alignfile, dual);
  }

  gssw_graph_destroy(graph);
  if (dual)gssw_graph_destroy(NVgraph);
  delete[] nt_table;
  delete[] mat;

  return (0);
}

void graphToDot(gssw_graph *graph, std::string output) {
  using namespace std;
  ofstream out;
  out.open(output.c_str());

  out << "digraph gssw {\n";
  out << "rankdir=\"LR\";\n";

  for (int i = 0; i < graph->size; i++) {
    out << graph->nodes[i]->id << " [label=\"" << graph->nodes[i]->data << ":" << graph->nodes[i]->seq << "\"];\n";
  }
  for (int i = 0; i < graph->size; i++) {
    for (int n = 0; n < graph->nodes[i]->count_next; n++)
      out << graph->nodes[i]->id << " -> " << graph->nodes[i]->next[n]->id << ";\n";
  }
  out << "}";
}

void align_main(gssw_graph *graph, gssw_graph *NVgraph,
                int8_t *mat, int8_t *nt_table,
                uint8_t gap_open, uint8_t gap_extension,
                std::string readfile, std::string alignfile, bool dual) {

  std::ifstream reads;
  std::ofstream aligns, NValigns;
  string read;
  int32_t readnum = 0, readEndPos;

  reads.open(readfile.c_str());
  aligns.open(alignfile.c_str());
  if (dual) NValigns.open((alignfile + ".nv").c_str());
  if (!reads.good() || !aligns.good() || (!NValigns.good() && dual)) {
    cerr << "Error opening reads file or alignment files. No alignment will be done." << endl;
    exit(0);
  }

  if (print) cout << "Aligning reads..." << endl;
  while (std::getline(reads, read)) {
    readnum++;
    if (print) cout << std::setw(12) << readnum << '\r' << std::flush;

    readEndPos = int32_t(read.find('#'));
    if (readEndPos == string::npos) readEndPos = int32_t(read.length());

    gssw_graph_fill(graph, read.substr(0, readEndPos).c_str(), nt_table, mat,
                    gap_open, gap_extension, uint32_t(read.length()), 2);

    aligns << read << ";"
        << graph->max_node->alignment->score << ","
        << (graph->max_node->data + 1 - graph->max_node->len + graph->max_node->alignment->ref_end) << ","
        << graph->maxCount << ";";
    if (graph->submax_node) {
      aligns << graph->submax_node->alignment->score << ","
          << (graph->submax_node->data + 1 -
              graph->submax_node->len +
              graph->submax_node->alignment->ref_end) << ","
          << graph->submaxCount;
    } else {
      aligns << "-1,-1,-1";
    }
    aligns << endl;

    /** Align to second graph if -D is specified **/
    if (dual) {
      gssw_graph_fill(NVgraph, read.substr(0, readEndPos).c_str(), nt_table, mat,
                      gap_open, gap_extension, int32_t(read.length()), 2);

      NValigns << read << ";"
          << NVgraph->max_node->alignment->score << ","
          << (NVgraph->max_node->data + 1 - NVgraph->max_node->len + NVgraph->max_node->alignment->ref_end) << ","
          << NVgraph->maxCount << ";";
      if (NVgraph->submax_node) {
        NValigns << NVgraph->submax_node->alignment->score << ","
            << (NVgraph->submax_node->data + 1 -
                NVgraph->submax_node->len +
                NVgraph->submax_node->alignment->ref_end) << ","
            << NVgraph->submaxCount;
      } else {
        NValigns << "-1,-1,-1";
      }
      NValigns << endl;

    }
  }
  reads.close();
  aligns.close();
  if (dual) NValigns.close();
}


void stat_main(std::string alignfile, int32_t tol) {

  string line;
  std::ifstream alignsin;
  std::vector<string> line_split(0), read_split(0), opt_split(0), subopt_split(0), files(0);
  int32_t readPos, optPos, correct, total,
      tCorrect = 0, tReads = 0, tCorrectNV = 0, tReadsNV = 0;

  split(alignfile, ',', files);

  for (int i = 0; i < files.size(); i++) {
    correct = 0;
    total = 0;
    alignsin.open(files[i].c_str());
    if (!alignsin.good()) {
      cerr << "Error opening input file " << alignfile << endl;
      exit(1);
    }

    while (std::getline(alignsin, line)) {
      if (line.at(0) != '#') {
        split(line, ';', line_split);
        split(line_split[0], ',', read_split);
        split(line_split[1], ',', opt_split);
        split(line_split[2], ',', subopt_split);
        readPos = atoi(read_split[2].c_str()) - atoi(read_split[1].c_str()) + atoi(read_split[3].c_str());
        optPos = atoi(opt_split[2].c_str()) - atoi(opt_split[1].c_str()) + atoi(opt_split[4].c_str());
        if (optPos > readPos - tol && optPos < readPos + tol) correct++;
        total++;
      }
    }
    cout << files[i] << ": " << correct << "/" << total << " correct." << endl;
    alignsin.close();
    tCorrect += correct;
    tReads += total;
    correct = 0;
    total = 0;

    alignsin.open((files[i] + ".nv").c_str());
    if (alignsin.good()) {
      while (std::getline(alignsin, line)) {
        if (line.at(0) != '#') {
          split(line, ';', line_split);
          split(line_split[0], ',', read_split);
          split(line_split[1], ',', opt_split);
          split(line_split[2], ',', subopt_split);
          readPos = atoi(read_split[2].c_str()) - atoi(read_split[1].c_str()) + atoi(read_split[3].c_str());
          optPos = atoi(opt_split[2].c_str()) - atoi(opt_split[1].c_str()) + atoi(opt_split[4].c_str());
          if (optPos > readPos - tol && optPos < readPos + tol) correct++;
          total++;
        }
      }
      cout << files[i] << ".nv: " << correct << "/" << total << " correct." << endl;
    } else {
      cout << "No-variant aligns " << files[i] << ".nv found." << endl;
    }
    alignsin.close();
    tCorrectNV += correct;
    tReadsNV += total;
    cout << endl;
  }
  cout << "Total: " << tCorrect << "/" << tReads << "(" << float(tCorrect) / tReads << "%) correct." << endl;
  cout << "Total NV: " << tCorrectNV << "/" << tReadsNV << "(" << float(tCorrectNV) / tReadsNV << "%) correct." <<
      endl;
  cout << "Tolerance: " << tol << endl;
}


std::string generateRead(gssw_graph &graph, int32_t readLen, float muterr, float indelerr) {
  gssw_node *node, *nodeCandidate;
  int32_t base, RAND, ambig = 0, currIndiv = -1;
  char mut;
  std::stringstream readmut;
  std::string read = "";
  bool valid;

  /** initial random node and base **/
  do {
    node = graph.nodes[rand() % (graph.size - 1)];
  } while(node->len < 1);
  base = rand() % (node->len);
  if (node->indivSize > 0) currIndiv = node->indiv[rand() % node->indivSize];

  for (int i = 0; i < readLen; i++) {
    read += node->seq[base];
    if (node->seq[base] == 'N') ambig++;
    base++;

    /** Go to next random node **/
    if (base == node->len) {
      do {
        nodeCandidate = node->next[rand() % node->count_next];
        valid = false;
        if (nodeCandidate->indivSize == 0) break;
        if(currIndiv < 0){
          RAND = rand() % nodeCandidate->indivSize;
          currIndiv = nodeCandidate->indiv[RAND];
          break;
        }
        for (int i = 0; i < nodeCandidate->indivSize; i++) {
          if (currIndiv == nodeCandidate->indiv[i]) {
            valid = true;
            break;
          }
        }
      } while (node->indivSize == 0 || !valid);
      node = nodeCandidate;
      if (node->count_next == 0) break; // End of graph reached
      base = 0;
    }
  }

  if (ambig > readLen / 2 || read.length() < readLen / 2) return generateRead(graph, readLen, muterr, indelerr);

  /** Mutate string **/
  for (int i = 0; i < read.length(); i++) {
    RAND = rand() % 100000;
    mut = read.at(i);

    if (RAND < (100000 - (100000 * indelerr / 2))) { // represents del
      /** Mutation **/
      if (RAND < (100000 * muterr) / 4) mut = 'A';
      else if (RAND < 2 * (100000 * muterr) / 4) mut = 'G';
      else if (RAND < 3 * (100000 * muterr) / 4) mut = 'C';
      else if (RAND < (100000 * muterr)) mut = 'T';
      readmut << mut;

      /* Insertion **/
      if (RAND > (100000 - (100000 * indelerr))) {
        RAND = rand() % int32_t(100000 * muterr);
        if (RAND < (100000 * muterr) / 4) mut = 'A';
        else if (RAND < 2 * (100000 * muterr) / 4) mut = 'G';
        else if (RAND < 3 * (100000 * muterr) / 4) mut = 'C';
        else if (RAND < (100000 * muterr)) mut = 'T';
        readmut << mut;
      }
    }
  }

  /** Append suffix recording read position **/
  readmut << '#' << node->data - node->len + base << ',' << currIndiv;
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
  if (print) cout << "Generating Nodes (" << buildfile << ")" << "..." << endl;
  while (getline(graphDat, line)) {
    if (line.at(0) != '#') {
      if (line.at(0) == '[') {
        line = line.substr(1, line.length() - 2);
        split(line, ',', lineSplit);
        for (int i = 0; i < lineSplit.size(); i++) {
          gssw_node_add_indiv(nodes.back(), strtol(lineSplit[i].c_str(), NULL, 10));
        }
      } else {
        split(line, ',', lineSplit);
        switch (lineSplit.size()) {
          case 3: // New node
            curr = uint32_t(strtol(lineSplit[1].c_str(), NULL, 10));
            nodes.push_back(gssw_node_create(int32_t(strtol(lineSplit[0].c_str(), NULL, 10)),
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
    }
  }

  if (print) cout << endl << "Building Graph..." << endl;

  /** Buffer node **/
  nodes.push_back(gssw_node_create(-1, ++curr, "", nt_table, mat));
  gssw_nodes_add_edge(nodes.end()[-2], nodes.end()[-1]);

  /** Add nodes to graph **/
  gssw_graph *graph = gssw_graph_create(uint32_t(nodes.size()));
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
    std::string outputFile,
    int32_t inGroup) {

  using namespace std;

  /** Used to track which indivs have a variant **/
  vector<int16_t> inVar(0);

  /** reference line, variant line **/
  string ref_line, vcf_line;
  /** Parsed lines, VCF header row **/
  vector<string> vline_split(0), altList_split(0), header(0);
  /** Vector of all the nodes in the graph **/
  vector<gssw_node *> nodes(0);
  /** Columns used to build graph **/
  vector<int32_t> inGroupCols(0);
  /** columns in VCF file, current ref pos **/
  int32_t posColumn = -1, refColumn = -1, altColumn = -1, formatColumn = -1, ref_position = 0;
  int32_t nodenum = 0, numIndivs;
  /** Pos from variant file **/
  int vpos;
  int32_t cn;
  /** To track edges that need to be built **/
  int numalts = 0, numprev;
  /** strings that represent node contents **/
  string variantRef, variantAlt, nodestring;
  int32_t nodelen = 0;
  char base;
  bool write = false;
  int32_t randtemp;

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

  getline(reference, ref_line);
  if (ref_line.at(0) != '>') cerr << "Error in ref file, first char should be >" << endl;

  /** Go to first VCF record **/
  do { getline(variants, vcf_line); } while (vcf_line.substr(0, 2) == "##");
  transform(vcf_line.begin(), vcf_line.end(), vcf_line.begin(), ::tolower);
  split(vcf_line, '\t', header);
  posColumn = int32_t(find(header.begin(), header.end(), "pos") - header.begin());
  refColumn = int32_t(find(header.begin(), header.end(), "ref") - header.begin());
  altColumn = int32_t(find(header.begin(), header.end(), "alt") - header.begin());
  formatColumn = int32_t(find(header.begin(), header.end(), "format") - header.begin());
  numIndivs = int32_t(header.size() - formatColumn - 1);

  /** Construct the in group, the graph will be built with these individuals **/
  if (write) out << '#';
  if (inGroup >= 0) {
    for (int i = 0; i < int32_t((numIndivs / 100.0f) * inGroup); i++) {
      randtemp = rand() % numIndivs + formatColumn + 1;
      if(find(inGroupCols.begin(), inGroupCols.end(), randtemp) == inGroupCols.end()) {
        inGroupCols.push_back(randtemp);
        out << inGroupCols[i] << ',';
      } else i--;
    }
    sort(inGroupCols.begin(), inGroupCols.end());
  } else {
    for (int i = 0; i < numIndivs; i++) {
      inGroupCols.push_back(i + formatColumn + 1);
      out << inGroupCols[i] << ',';
    }
  }
  if (write) out << endl;

  /** Find the POS, REF, ALT cols **/
  if (posColumn < 0 || refColumn < 0 || altColumn < 0 || formatColumn < 0) {
    cerr << "POS, REF, ALT, and/or INFO not found in VCF header." << endl;
    exit(1);
  }

  /** Generate Nodes **/
  if (print) cout << "Generating nodes..." << endl;

  /** Go to minimum position **/
  if (minpos > 0) {
    while (ref_position < minpos - 1) {
      reference.get(base);
      if (!isspace(base)) ref_position++;
    }
  }

  /** Process variants **/
  while (getline(variants, vcf_line)) {
    nodestring = "";
    nodelen = 0;
    split(vcf_line, '\t', vline_split);
    vpos = atoi(vline_split[posColumn].c_str());
    if (vpos <= ref_position) goto endvar;
    if (vpos > maxpos) break;
    variantRef = vline_split[refColumn];
    variantAlt = vline_split[altColumn];
    if (print) cout << setw(12) << nodenum << '\r' << flush;

    /** build node string up to var pos **/
    while (ref_position < vpos - 1) {
      if (!reference.get(base)) {
        cerr << "End of ref found while looking for variant pos " << vpos << endl;
        exit(1);
      }
      if (!isspace(base)) {
        nodestring += base;
        ref_position++;
        nodelen++;
      }

      /** Max node length reached, split node **/
      if (nodelen == maxNodeLen) {
        nodes.push_back(gssw_node_create(ref_position, nodenum, nodestring.c_str(), nt_table, mat));
        if (write) out << ref_position << "," << nodenum << "," << nodestring.c_str() << endl;
#if debug > 4
        cout << "Node: " << ref_position << ", ID: " << nodenum << ", " << nodestring << endl;
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
      nodes.push_back(gssw_node_create(ref_position, nodenum, nodestring.c_str(), nt_table, mat));
      if (write) out << ref_position << "," << nodenum << "," << nodestring.c_str() << endl;
#if debug > 4
      cout << "Node: " << ref_position << ", ID: " << nodenum << ", " << nodestring << endl;
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
    for (int i = 0; i < variantRef.length(); i++) {
      reference.get(base);
      if (isspace(base)) {
        reference.get(base);
      }
      ref_position++;
    }
    nodes.push_back(gssw_node_create(ref_position, nodenum, variantRef.c_str(), nt_table, mat));
    if (write) out << ref_position << "," << nodenum << "," << variantRef.c_str() << endl;
#if debug > 4
    cout << "Node: " << ref_position << ", ID: " << nodenum << ", " << variantRef << endl;
#endif
    nodenum++;
    numalts = 1;

    /** Variants **/
    if (!novar) {
      split(variantAlt, ',', altList_split);
      for (int i = 0; i < altList_split.size(); i++) {
        inVar.clear();
        for (int c = 0; c < inGroupCols.size(); c++) {
          /** Check if it is in the ingroup **/
          //TODO parse rather than at, check diploid
          if (strtol(&vline_split[inGroupCols[c]].at(0), NULL, 10) == i + 1) {
            inVar.push_back(inGroupCols[c]);
          }
        }
        if (inVar.size() > 0) {
          if(altList_split[i].substr(0, 3) == "<CN") {
            cn = int32_t(strtol(altList_split[i].substr(3, altList_split[i].length() - 4).c_str(), NULL, 10));
            altList_split[i] = "";
            for (int v = 0; v < cn; v++){
              altList_split[i] += variantRef;
            }
          }
          nodes.push_back(gssw_node_create(ref_position, nodenum, altList_split[i].c_str(), nt_table, mat));
          if (write) out << ref_position << "," << nodenum << "," << altList_split[i].c_str() << endl;
#if debug > 4
          cout << "Node: " << ref_position << ", ID: " << nodenum << ", " << altList_split[i].c_str() << endl;
#endif
          nodenum++;
          numalts++;

          /** Add the individuals that have this particular variant **/
          if(write) out << '[';
          for (int n = 0; n < inVar.size(); n++) {
            gssw_node_add_indiv(nodes.back(), inVar[n]);
#if debug > 5
            cout << "Add Indiv(" << int32_t(nodes.back()->indivSize) << "): " << inVar[n] << endl;
#endif
            if (write) out << inVar[n];
            if (n != inVar.size() - 1 && write) out << ',';
          }
          if (write) out << ']' << endl;

        }
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
  while ((ref_position < maxpos || maxpos < 0) && reference.get(base)) {
    if (!isspace(base)) {
      nodestring += base;
      ref_position++;
    }
  }
  nodes.push_back(gssw_node_create(ref_position, nodenum, nodestring.c_str(), nt_table, mat));
  if (write) out << ref_position << "," << nodenum << "," << nodestring.c_str() << endl;
#if debug > 4
  cout << "Node: " << ref_position << ", ID: " << nodenum << ", " << nodestring.c_str() << endl;
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
  nodes.push_back(gssw_node_create(ref_position, nodenum, "", nt_table, mat));
  gssw_nodes_add_edge(nodes.end()[-2], nodes.end()[-1]);
#if debug > 4
  cout << "Node: " << ref_position << ", ID: " << nodenum << ", " << "" << endl;
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
  using std::cout;
  using std::endl;
  cout << "Node sequence: " << node->seq << endl;
  cout << "Score: " << node->alignment->score << " end: " << node->data + 1 - node->len + node->alignment->ref_end <<
      endl;
}


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  /** Split string with delim, return a vector **/
  std::stringstream ss(s);
  std::string item;
  elems = *new std::vector<std::string>(0);

  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  if (s.at(s.length() - 1) == ',') elems.push_back("");
  return elems;
}
