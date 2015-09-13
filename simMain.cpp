//
// Created by gaddra on 9/12/15.
//


#include "simMain.h"

void sim_main(int argc, char *argv[]) {
  /** Default sim read error **/
  float muterr = 0.01f, indelerr = 0.01f;

  /** Graph to build reads from **/
  std::string buildfile = "";

  /** Read generation defaults **/
  int32_t numreads = 1000, readlen = 100, tol = 5;

  /** Scores **/
  int32_t match = 2, mismatch = 2;

  /** Graph score and conversion table **/
  int8_t *nt_table = gssw_create_nt_table();
  int8_t *mat = gssw_create_score_matrix(match, mismatch);

  srand(uint32_t(time(NULL)));
  GetOpt::GetOpt_pp args(argc, argv);

  if ((args >> GetOpt::OptionPresent('h', "help"))) {
    printSimHelp();
    exit(0);
  }

  if (!(args >> GetOpt::Option('b', "buildfile", buildfile))) {
    std::cerr << "Error: no build file defined." << std::endl;
    exit(1);
  }

  args >> GetOpt::Option('n', "numreads", numreads)
      >> GetOpt::Option('m', "muterr", muterr)
      >> GetOpt::Option('i', "indelerr", indelerr)
      >> GetOpt::Option('l', "readlen", readlen);

  gssw_graph *graph = buildGraph(buildfile, nt_table, mat);

  /** Simulate reads **/

  std::cerr << "Generating reads..." << std::endl;
  for (int i = 0; i < numreads; i++) {
    std::cout << generateRead(*graph, readlen, muterr, indelerr) << std::endl;
  }

  gssw_graph_destroy(graph);
  delete[] nt_table;
  delete[] mat;
}


std::string generateRead(gssw_graph &graph, int32_t readLen, float muterr, float indelerr) {
  gssw_node *node, *nodeCandidate;
  int32_t base, RAND, ambig = 0, currIndiv = -1;
  char mut;
  std::stringstream readmut, path;
  std::string read = "";
  bool valid;

  /** initial random node and base **/
  do {
    node = graph.nodes[rand() % (graph.size - 1)];
  } while (node->len < 1);
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
        if (currIndiv < 0) {
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


void printSimHelp() {
  using std::cout;
  using std::endl;
  cout << endl << "------------------- VMatch sim, September 2015. rgaddip1@jhu.edu -------------------" << endl;
  cout << "-b\t--buildfile     quick rebuild file, required if -v, -r are not defined." << endl;
  cout << "-n\t--numreads      Number of reads to simulate" << endl;
  cout << "-m\t--muterr        Simulated read mutation error rate" << endl;
  cout << "-i\t--indelerr      Simulated read Indel error rate" << endl;
  cout << "-l\t--readlen       Nominal read length" << endl << endl;
}