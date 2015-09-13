//
// Created by gaddra on 9/12/15.
//

#include "alignMain.h"

void align_main(int argc, char *argv[]) {
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

  GetOpt::GetOpt_pp args(argc, argv);

  if (args >> GetOpt::OptionPresent('h', "help")) {
    printAlignHelp();
    exit(0);
  }
  /** make sure there's a valid input **/
  if (!(args >> GetOpt::Option('b', "buildfile", buildfile))) {
    cerr << "Error: no build file defined." << endl;
    exit(1);
  }

  args >> GetOpt::Option('m', "match", match)
      >> GetOpt::Option('n', "mismatch", mismatch)
      >> GetOpt::Option('o', "gap_open", gap_open)
      >> GetOpt::Option('e', "gap_extend", gap_extension)
      >> GetOpt::Option('r', "reads", readsfile);

  graph = buildGraph(buildfile, nt_table, mat);
  align(graph, mat, nt_table, gap_open, gap_extension, readsfile);

  gssw_graph_destroy(graph);
  delete[] nt_table;
  delete[] mat;

}

void align(gssw_graph *graph,
           int8_t *mat, int8_t *nt_table,
           uint8_t gap_open, uint8_t gap_extension,
           std::string readfile) {
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::string;

  std::ifstream reads;
  string read;
  int32_t readnum = 0, readEndPos;

  reads.open(readfile.c_str());
  if (!reads.good()) {
    cerr << "Error opening reads file. No alignment will be done." << endl;
    exit(0);
  }

  cerr << "Aligning reads..." << endl;
  while (std::getline(reads, read)) {
    readnum++;

    readEndPos = int32_t(read.find('#'));
    if (readEndPos == string::npos) readEndPos = int32_t(read.length());

    gssw_graph_fill(graph, read.substr(0, readEndPos).c_str(), nt_table, mat,
                    gap_open, gap_extension, uint32_t(read.length()), 2);

    cout << read << ";"
        << graph->max_node->alignment->score << ","
        << (graph->max_node->data + 1 - graph->max_node->len + graph->max_node->alignment->ref_end) << ","
        << graph->maxCount << ";";
    if (graph->submax_node) {
      cout << graph->submax_node->alignment->score << ","
          << (graph->submax_node->data + 1 -
              graph->submax_node->len +
              graph->submax_node->alignment->ref_end) << ","
          << graph->submaxCount;
    } else {
      cout << "-1,-1,-1";
    }
    cout << endl;
  }
  reads.close();
}

void printAlignHelp() {
  using std::cout;
  using std::endl;
  cout << endl << "------------------- VMatch align, September 2015. rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     quick rebuild file, required if -v, -r are not defined." << endl;
  cout << "-m\t--match         Match score, default  " << 2 << endl;
  cout << "-n\t--mismatch      Mismatch score, default " << 2 << endl;
  cout << "-o\t--gap_open      Gap opening score, default " << 3 << endl;
  cout << "-e\t--gap_extend    Gap extend score, default " << 1 << endl;
  cout << "-r\t--reads         Reads to align, one per line. Symbols after '#' are ignored." << endl;
  cout << endl << "Alignments output to stdout." << endl;

}