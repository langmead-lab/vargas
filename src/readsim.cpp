/**
 * Ravi Gaddipati
 * November 24, 2015
 * rgaddip1@jhu.edu
 *
 * ReadSim simulates reads given a gssw_graph or a Graph object.
 * Reads can either be from a single individual or a random walk
 * with adjustable error rates.
 *
 * readsim.cpp
 */

#include "../include/readsim.h"

bool vargas::ReadSim::updateRead() {
  if (!graph) throw std::invalid_argument("No graph assigned.");
  if (readProfiles.size() > 0) {
    generateRead();
    for (auto &prof : readProfiles) {
      if (counters.at(&prof) < p.maxreads) {
        if (prof == read) {
          counters[&prof]++;
          totalreads++;
          *(logs[&prof]) << str() << std::endl;
          break;
        }
      }
      if (totalreads >= readProfiles.size() * p.maxreads) return false;
    }
  }
  else generateRead();
  return true;
}

void vargas::ReadSim::generateRead() {
  //TODO segfault when read length is too short
  gssw_node *node, *nodeCandidate;
  int32_t base, RAND, ambig = 0, currIndiv = -1, numSubErr = 0, numVarNodes = 0, numVarBases = 0;
  char mut;
  std::stringstream readmut;
  std::string read = "";
  bool valid;

  /** initial random node and base **/
  do {
    node = (*graph).nodes[rand() % ((*graph).size - 1)];
  } while (node->len < 1);
  base = rand() % (node->len);
  if (node->indivSize > 0) {
    currIndiv = node->indiv[rand() % node->indivSize];
    numVarNodes++;
  }

  for (int i = 0; i < p.readLen; i++) {
    if (node->len != 0) read += node->seq[base];
    if (node->indivSize > 0) numVarBases++;
    if (node->seq[base] == 'N') ambig++;
    base++;

    /** Go to next random node **/
    if (base == node->len) {
      if (node->count_next == 0) break;
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
          if (p.randWalk || currIndiv == nodeCandidate->indiv[i]) {
            valid = true;
            break;
          }
        }
      } while (node->indivSize == 0 || !valid);
      node = nodeCandidate;
      if (node->indivSize > 0) numVarNodes++;
      base = 0;
    }
  }

  //Replace read if it has less than non-ambiguous/ half bases
  if (ambig > p.readLen / 2 || read.length() < p.readLen / 2) {
    generateRead();
    return;
  }

  /** Mutate string **/
  for (uint32_t i = 0; i < read.length(); i++) {
    RAND = rand() % 1000;
    mut = read.at(i);

    // If its below (1000*error), there's a deletion
    if (RAND > 1000.0 * p.indelerr) {
      // Substitution errors
      if (RAND < 250.0 * p.muterr && mut != 'A') {
        mut = 'A';
        numSubErr++;
      }
      else if (RAND < 500.0 * p.muterr && mut != 'G') {
        mut = 'G';
        numSubErr++;
      }
      else if (RAND < 750.0 * p.muterr && mut != 'C') {
        mut = 'C';
        numSubErr++;
      }
      else if (RAND < 1000.0 * p.muterr && mut != 'T') {
        mut = 'T';
        numSubErr++;
      }
      readmut << mut;
    }

    // Insertion
    RAND = rand() % 1000;
    if (RAND < 1000.0 * p.indelerr) {
      RAND = rand() % 1000;
      if (RAND < 250) mut = 'A';
      else if (RAND < 500) mut = 'G';
      else if (RAND < 750) mut = 'C';
      else mut = 'T';
      readmut << mut;
    }

  }

  /** populate read info **/
  this->read.read = readmut.str();
  readmut.str(std::string());
  this->read.readEnd = uint32_t(node->data - node->len + base);
  this->read.indiv = currIndiv;
  this->read.numSubErr = numSubErr;
  this->read.numVarNodes = numVarNodes;
  this->read.numVarBases = numVarBases;
}