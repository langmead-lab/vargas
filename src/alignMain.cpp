//
// Created by gaddra on 9/12/15.
//

#include "../include/alignMain.h"

int align_main(int argc, char *argv[]) {
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::string;

  /** Scores **/
  int32_t match = 2, mismatch = 2;
  uint8_t gap_open = 3, gap_extension = 1;

  /** Graph score and conversion table **/
  int8_t *nt_table = gssw_create_nt_table();
  int8_t *mat = gssw_create_score_matrix(match, mismatch);

  string buildfile = "", readsfile = "";
  gssw_graph *graph;
  bool readsin = false;

  GetOpt::GetOpt_pp args(argc, argv);

  if (args >> GetOpt::OptionPresent('h', "help")) {
    printAlignHelp();
    return 0;
  }
  /** make sure there's a valid input **/
  if (!(args >> GetOpt::Option('b', "buildfile", buildfile))) {
    cerr << "Error: no build file defined." << endl;
    return 1;
  }

  args >> GetOpt::Option('m', "match", match)
      >> GetOpt::Option('n', "mismatch", mismatch)
      >> GetOpt::Option('o', "gap_open", gap_open)
      >> GetOpt::Option('e', "gap_extend", gap_extension);

  if ((args >> GetOpt::Option('r', "reads", readsfile))) readsin = true;

  ;//  graph = buildGraph(buildfile, nt_table, mat);
  align(graph, mat, nt_table, gap_open, gap_extension, readsfile, readsin);

  gssw_graph_destroy(graph);
  free(nt_table);
  free(mat);

  return 0;
}

void align(gssw_graph *graph,
           int8_t *mat, int8_t *nt_table,
           uint8_t gap_open, uint8_t gap_extension,
           std::string readsfile, bool readsin) {
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::string;

  std::ifstream reads;
  if (readsin) reads.open(readsfile.c_str());

  string read;
  int32_t readEndPos, optAlignEnd, suboptAlignEnd;
  std::vector<string> readMeta(0);
  int32_t readSampleLoc; //TODO May overflow with whole genome
  int32_t tol;

  cerr << "Aligning reads..." << endl;
  while ((readsin ? std::getline(reads, read) : std::getline(std::cin, read))) {
    if (read.length() > 0 && read.at(0) != '#') {
      readSampleLoc = -1;
      readEndPos = int32_t(read.find('#'));
      if (readEndPos == string::npos) readEndPos = int32_t(read.length());
      else {
        // Get the original read location, else -1
        readMeta = split(read.substr(readEndPos + 1, read.length() - readEndPos - 1), ',');
        readSampleLoc = readMeta[0].length() > 0 ? atoi(readMeta[0].c_str()) : -1;
      }
      tol = (readEndPos / 10) + 1;

      gssw_graph_fill(graph, read.substr(0, readEndPos).c_str(), nt_table, mat,
                      gap_open, gap_extension, uint32_t(read.length()), 2, readSampleLoc);

      optAlignEnd = graph->max_node->data + 1 - graph->max_node->len + graph->max_node->alignment->ref_end;
      suboptAlignEnd = graph->submax_node->data + 1 - graph->submax_node->len + graph->submax_node->alignment->ref_end;

      if (readSampleLoc > 0) {
        if (readSampleLoc > optAlignEnd - tol && readSampleLoc < optAlignEnd + tol) readSampleLoc = 0;
        else if (readSampleLoc > suboptAlignEnd - tol && readSampleLoc < suboptAlignEnd + tol) readSampleLoc = 1;
        else readSampleLoc = 2;
      }
      cout << read << ","
          << graph->max_node->alignment->score << ","
          << optAlignEnd << ","
          << graph->maxCount << ",";
      if (graph->submax_node) {
        cout << graph->submax_node->alignment->score << ","
            << suboptAlignEnd << ","
            << graph->submaxCount << ",";
      } else {
        cout << "-1,-1,-1,";
      }
      cout << readSampleLoc << endl;
    }
  }
}

void printAlignHelp() {
  using std::cout;
  using std::endl;
  cout << endl << "------------------- VMatch align, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     quick rebuild file, required if -v, -r are not defined." << endl;
  cout << "-m\t--match         Match score, default  " << 2 << endl;
  cout << "-n\t--mismatch      Mismatch score, default " << 2 << endl;
  cout << "-o\t--gap_open      Gap opening score, default " << 3 << endl;
  cout << "-e\t--gap_extend    Gap extend score, default " << 1 << endl;
  cout << "-r\t--reads         Reads to align. If not specified, read from stdin" << endl;

  cout << endl << "Alignments output to stdout. Reads read from stdin or -r, 1 per line." << endl;
  cout << "Lines beginning with \'#\' are ignored." << endl;
  cout << "Output format:" << endl;
  cout << "READ,OPTIMAL_SCORE,OPTIMAL_ALIGNMENT_END,NUM_OPTIMAL_ALIGNMENTS,SUBOPTIMAL_SCORE," << endl;
  cout << "SUBOPTIMAL_ALIGNMENT_END,NUM_SUBOPTIMAL_ALIGNMENTS,ALIGNMENT_MATCH" << endl << endl;
  cout << "ALIGNMENT_MATCH: 0- optimal match, 1- suboptimal match, 2- no match" << endl << endl;

}