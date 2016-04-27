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

#include <map>
#include <random>
#include "../include/graph.h"
#include "../include/readsim.h"


void vargas::ReadSim::populateProfiles() {
  if (!graph) throw std::invalid_argument("No graph assigned.");
  while (true) {
    updateRead();
    if (readProfiles.size() > 0) {
      for (auto profile : readProfiles) {
        if (counters[profile] < p.maxreads) {
          if (*profile == read) {
            counters[profile]++;
            totalreads++;
            *(logs[profile]) << this->toString() << std::endl;
            break;
          }
        }
      }
      if (totalreads >= readProfiles.size() * p.maxreads) break;
    }
  }
}


bool vargas::ReadSim::updateRead() {
  //TODO segfault when read length is too short
  if (!graph) return false;
  gssw_node *node, *nodeCandidate;
  int32_t base, RAND, ambig = 0, currIndiv = -1, numSubErr = 0, numVarNodes = 0, numVarBases = 0, numIndelErr = 0;
  char mutatedBase;
  std::stringstream mutatedRead;
  std::string read = "";
  std::vector<uint32_t> individuals;
  vargas::Xcoder coder;
  bool valid;

  /** initial random node and base **/
  do {
    node = graph->nodes[rand() % graph->size];
  } while (node->len < 1); // The graph allows empty nodes

  base = rand() % node->len;
  // Pick an individual to follow
  if (node->indivCompressedSize > 0) {
    coder.inflate(node->indivCompressed, node->indivCompressedSize, individuals);
    currIndiv = individuals[rand() % individuals.size()];
    numVarNodes++;
  }

  for (int i = 0; i < p.readLen; i++) {
    if (node->len != 0) {
      read += node->seq[base];
    }
    if (node->indivCompressedSize > 0) {
      numVarBases++;
    }
    if (node->seq[base] != 'A' && node->seq[base] != 'T' && node->seq[base] != 'G' && node->seq[base] != 'C') {
      ambig++;
    }

    base++;

    // Reached end of node, go to next random node
    if (base == node->len) {
      if (node->count_next == 0) break;
      do {
        // Find a node that the individual possess
        nodeCandidate = node->next[rand() % node->count_next];
        // Get the list of individuals for the node candidate
        coder.inflate(nodeCandidate->indivCompressed, nodeCandidate->indivCompressedSize, individuals);
        if (nodeCandidate->indivCompressedSize == 0) break; // A node common to all individuals is valid

        valid = false;

        // If we haven't already picked an individual, pick one
        if (currIndiv < 0) {
          currIndiv = individuals[rand() % individuals.size()];
          valid = true;
          break;
        }

//        Can probs replace lower block with this, need to confirm indiv. is always
//        sorted. Can also get rid of valid flag.
//
//        if(p.randWalk || std::binary_search(individuals.begin(), individuals.end(), currIndiv)){
//          valid = true;
//        }

        for (int i = 0; i < individuals.size(); i++) {
          if (p.randWalk || currIndiv == individuals[i]) {
            valid = true;
            break;
          }
        }

      } while (!valid);
      node = nodeCandidate;
      if (node->indivCompressedSize > 0) numVarNodes++;
      base = 0;
    }
  }

  // Replace read if it has less than half non-ambiguous bases
  if (ambig > p.ambiguity || read.length() != p.readLen) {
    return updateRead();
  }

  // Mutate string
  for (uint32_t i = 0; i < read.length(); i++) {
    RAND = rand() % 1000;
    mutatedBase = read.at(i);


    // If its below (1000*error), there's a deletion
    if (RAND >= 1000.0 * p.indelerr) {
      // Substitution errors
      if (RAND < 250.0 * p.muterr && mutatedBase != 'A') {
        mutatedBase = 'A';
        numSubErr++;
      }
      else if (RAND < 500.0 * p.muterr && mutatedBase != 'G') {
        mutatedBase = 'G';
        numSubErr++;
      }
      else if (RAND < 750.0 * p.muterr && mutatedBase != 'C') {
        mutatedBase = 'C';
        numSubErr++;
      }
      else if (RAND < 1000.0 * p.muterr && mutatedBase != 'T') {
        mutatedBase = 'T';
        numSubErr++;
      }
      mutatedRead << mutatedBase;

      /* Insertion
       * The insertion block is in here to prevent a deletion + insertion (aka substitution)
       */
      RAND = rand() % 1000;
      if (RAND < 1000.0 * p.indelerr) {
        RAND = rand() % 1000;
        if (RAND < 250) mutatedBase = 'A';
        else if (RAND < 500) mutatedBase = 'G';
        else if (RAND < 750) mutatedBase = 'C';
        else mutatedBase = 'T';
        mutatedRead << mutatedBase;
        numIndelErr++;
      }
    } else {
      // Skipped the append, so deletion
      numIndelErr++;
    }

  }

  // populate read info
  this->read.read = mutatedRead.str();
  mutatedRead.str(std::string());
  this->read.readEndPos = uint32_t(node->data - node->len + base);
  this->read.indiv = currIndiv;
  this->read.numSubErr = numSubErr;
  this->read.numVarNodes = numVarNodes;
  this->read.numVarBases = numVarBases;
  this->read.numIndelErr = numIndelErr;

  return true;
}
