/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date November 20, 2015
 *
 * vargas::Graph is a DAG representation of a reference and its variants.
 * The class wraps a gssw_graph from gssw and provides a way to construct
 * graphs from a FASTA and a VCF file with various options.
 *
 * GSSW was originally written by Erik Garrison and was moderately modified.
 *
 * graph.cpp
 */

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "../include/utils.h"
#include "../include/graph.h"

/*********************************** CONSTRUCTOR ***********************************/

vargas::Graph::Graph(const Graph &g, const Alignment &a) {
  std::vector<gssw_node *> N;
  std::vector<std::pair<gssw_node *, gssw_node *>> oldLinks;
  std::map<gssw_node *, gssw_node *> oldToNewMap;

  /*
   * Get all the nodes that are in the locality of the alignment
   * Since Graph should be topographically ordered, inserting as we crawl the graph
   *
   * EndCounter sets a limit to how far we can go after a first 'cap' node, to prevent crawling the entire graph.
   */
  //TODO can look for beginning node via binary search variant
  gssw_node *n;
  int id = 0, endCounter = 0, len = a.read.read.length();
  len += len / 10;
  for (uint32_t i = 0; i < g.graph->size; ++i) {
    n = g.graph->nodes[i];
    if (endCounter) ++endCounter;
    if (endCounter > 10) break; // Don't look 10 nodes past the end

    // If the node contains the start of the local region
    if (n->data > a.optAlignEnd - len
        && n->data - n->len <= a.optAlignEnd - len) {
      uint32_t nSeqLen = n->data - (a.optAlignEnd - len); // The number of bases in this node we need
      char *seq = (char *) malloc(nSeqLen + 1);
      strncpy(seq, n->seq + n->len - nSeqLen, nSeqLen);
      seq[nSeqLen] = 0;
      N.push_back(gssw_node_create(n->data, id++, seq, this->params.nt_table, this->params.mat));
      free(seq);

      oldToNewMap[n] = N.back();
      for (int i = 0; i < n->count_next; ++i) {
        oldLinks.push_back(std::pair<gssw_node *, gssw_node *>(n, n->next[i]));
      }
    }

      // If we need an entire node
    else if (a.optAlignEnd > n->data
        && n->data - n->len > a.optAlignEnd - len) {
      N.push_back(gssw_node_create(n->data, id++, n->seq, this->params.nt_table, this->params.mat));

      oldToNewMap[n] = N.back();
      for (int i = 0; i < n->count_next; ++i) {
        oldLinks.push_back(std::pair<gssw_node *, gssw_node *>(n, n->next[i]));
      }
    }

      // The end node
    else if (n->data >= a.optAlignEnd
        && n->data - n->len < a.optAlignEnd) {
      uint32_t nSeqLen = a.optAlignEnd - (n->data - n->len);
      char *seq = (char *) malloc(nSeqLen + 1);
      strncpy(seq, n->seq, nSeqLen);
      seq[nSeqLen] = 0;
      N.push_back(gssw_node_create(n->data - n->len + nSeqLen, id++, seq, this->params.nt_table, this->params.mat));
      free(seq);
      ++endCounter;
      oldToNewMap[n] = N.back();
    }

  }

  // Add links
  for (auto &lnk : oldLinks) {
    if (oldToNewMap.count(lnk.first) && oldToNewMap.count(lnk.second)) {
      gssw_nodes_add_edge(oldToNewMap[lnk.first], oldToNewMap[lnk.second]);
    }
  }

  this->graph = gssw_graph_create(N.size());
  for (auto node : N) {
    gssw_graph_add_node(this->graph, node);
  }

}

/*********************************** BUILD GRAPH ***********************************/

void vargas::Graph::buildGraph(std::istream &graphDat) {

  using std::string;
  using std::endl;
  using std::vector;
  using std::cerr;

  if (graph != NULL) {
    gssw_graph_destroy(graph);
    graph = NULL;
  }

  string line, seq;
  vector<string> lineSplit(0);
  vector<gssw_node *> nodes(0);
  uint32_t curr = 0;

  vargas::Xcoder coder; // Used to store individuals

  uint8_t *data = NULL;
  size_t dataLen;

  /** Build nodes and edges from buildfile **/
  while (getline(graphDat, line)) {
    if (line.at(0) == '#') continue; // skip comment line
    if (line.at(0) == ':') {
      // Compressed individual data
        if (params.includeIndividuals) {
          // Add individuals to the last node
          line = line.substr(1); // Remove indiv line marker
          dataLen = coder.decode(line, &data); // Allocates data
          gssw_node_set_indivs(nodes.back(), data, dataLen); // stores data pointer
        }
      } else {
        split(line, ',', lineSplit);
        switch (lineSplit.size()) {
          case 3: // New node
            curr = uint32_t(strtoul(lineSplit[1].c_str(), NULL, 10));
            seq = (lineSplit[2] == "-") ? "" : lineSplit[2].c_str();
            nodes.push_back(gssw_node_create(int32_t(strtol(lineSplit[0].c_str(), NULL, 10)),
                                             curr,
                                             seq.c_str(), params.nt_table, params.mat));
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

  /** Add nodes to graph **/
  graph = gssw_graph_create(nodes.size());
  for (auto &n : nodes) {
    gssw_graph_add_node(graph, n);
  }

}

/*********************************** EXPORTERS *************************************/

void vargas::Graph::exportDOT(std::ostream &out, std::string name) const {

  if (graph == NULL) {
    std::cerr << "Error: No graph has been built. Aborting export." << std::endl;
    return;
  }

  out << "digraph \"" << name << "\" {\n";
  out << "rankdir=\"LR\";\n";

  for (uint32_t i = 0; i < graph->size; i++) {
    out << graph->nodes[i]->id << " [label=\"" << graph->nodes[i]->data << ":" << graph->nodes[i]->seq
        << "\" shape=plaintext];\n";
  }
  for (uint32_t i = 0; i < graph->size; i++) {
    for (int32_t n = 0; n < graph->nodes[i]->count_next; n++)
      out << graph->nodes[i]->id << " -> " << graph->nodes[i]->next[n]->id << ";\n";
  }
  out << "labelloc=\"t\";" << std::endl;
  out << "label=\"" << name << "\";" << std::endl;
  out << "}\n";
}


void vargas::Graph::exportBuildfile(std::istream *_reference, vcfstream &variants, std::ostream &buildout) {
  using std::vector;
  using std::endl;
  using std::string;

  std::istream &reference = *_reference;

  vargas::Xcoder encoder; // Used for compressing individual list

  // Get region
  uint32_t minpos, maxpos;
  parseRegion(params.region, &minpos, &maxpos);

  // Create the ingroup
  generateIngroup(variants);
  variants.printIngroup(buildout);

  // Check FASTA header
  std::string refLine;
  getline(reference, refLine);
  if (refLine.at(0) != '>') {
    throw std::invalid_argument("Error in ref file, first char should be >");
  }

  // Buildfile comments
  buildout << "##Generated by vargas, built on " << __DATE__ << endl;
  buildout << "##Reference sequence: " << refLine.substr(1) << endl;
  buildout << "##Graph params: "
  << "Max Node Len: " << params.maxNodeLen
  << ", Use max AF? " << params.maxAF
  << ", Include ref indivs? " << params.includeRefIndivs;
  if (params.region.length() > 0) buildout << ", Region: " << params.region;
  if (params.genComplement) buildout << ", Complement source: " << params.complementSource;
  else buildout << ", In-group Percent: " << params.ingroup;
  buildout << endl;


  /** Go to minimum position **/
  //TODO this can be made faster
  uint32_t currentRefPosition = 0;
  char base;
  if (minpos > 0) {
    while (currentRefPosition < minpos - 1) {
      reference.get(base);
      if (!isspace(base)) currentRefPosition++;
    }
  }

  /** Process variants **/
  vcfrecord variantRecord;
  string nodestring = "";
  uint32_t nodenum = 0;
  int32_t numUnconnectedPrev = 0, numUnconnectedCurr = 0;
  // For max AF graphs
  string maxAFAllele;
  double_t maxAF;

  while (variants.getRecord(variantRecord)) {
    if (variantRecord.pos <= currentRefPosition) continue;
    if (variantRecord.pos > maxpos) break;

    /** build node string up to variant position **/
    while (currentRefPosition < variantRecord.pos - 1) {
      do {
        if (!reference.get(base)) {
          throw std::range_error("End of ref found while looking for pos " + std::to_string(variantRecord.pos));
        }
      } while (isspace(base));

      nodestring += base;
      currentRefPosition++;

      /** Max node length reached, split node **/
      if (nodestring.length() == params.maxNodeLen) {
        buildout << currentRefPosition << "," << nodenum << "," << nodestring << endl;
          /** Connect to all of the previous alt/ref nodes **/
        for (int i = 0; i < numUnconnectedPrev; ++i) {
          buildout << -2 - i << "," << -1 << endl;
          }
        numUnconnectedPrev = 1; // We have one node at the tail of the graph
        nodenum++;
        nodestring = "";
      }
    }

    /** Make a node for the rest of the common string **/
    if (nodestring.length() > 0) {
      buildout << currentRefPosition << "," << nodenum << "," << nodestring << endl;
        /** Connect to all of the previous alt/ref nodes **/
      for (int i = 0; i < numUnconnectedPrev; ++i) {
        buildout << -2 - i << "," << -1 << endl;
        }
      numUnconnectedPrev = 1;
      nodenum++;
      nodestring = "";
    }

    std::string faref = "";
    /** progress reference file position past the ref that's stored in the variant file **/
    for (unsigned int i = 0; i < variantRecord.ref.length(); ++i) {
      do {
        reference.get(base);
      } while (isspace(base));

      faref += base;
      currentRefPosition++;
    }
    if (variantRecord.ref != faref)
      throw std::invalid_argument(
          "VCF reference does not match at " + std::to_string(variantRecord.pos) + ". REF: " + faref + ", VCF: "
              + variantRecord.ref);

    /** Variants and ref **/
    numUnconnectedCurr = 0;
    if (params.maxAF) {
      maxAFAllele = variantRecord.ref;
      maxAF = variantRecord.freqs[variantRecord.ref];
      for (auto &alt : variantRecord.freqs) {
        if (alt.second > maxAF) {
          maxAF = alt.second;
          maxAFAllele = alt.first;
        }
      }

      buildout << currentRefPosition << "," << nodenum << "," << maxAFAllele << endl;
      nodenum++;
      numUnconnectedCurr++;
      // Print individuals that have the maxAF allele
      if (params.includeRefIndivs) printEncodedVec(buildout, variantRecord.indivs[maxAFAllele], encoder);

    } else {
      for (auto &var : variantRecord.indivs) {
        buildout << currentRefPosition << "," << nodenum << "," << var.first << endl;
        nodenum++;
        numUnconnectedCurr++;
        // Print individuals that have this variant
        if (params.includeRefIndivs || var.first != variantRecord.ref) printEncodedVec(buildout, var.second, encoder);
      }
    }

    /** Build edges **/
    for (int p = 0; p < numUnconnectedPrev; p++) {
      for (int a = 0; a < numUnconnectedCurr; a++) {
        buildout << -1 - numUnconnectedCurr - p << "," << -1 - a << endl;
      }
    }
    numUnconnectedPrev = numUnconnectedCurr;
  }

  /** The remaining bases after the last variant **/
  nodestring = "";
  while ((currentRefPosition < maxpos) && reference.get(base)) {
    if (!isspace(base)) {
      nodestring += base;
      currentRefPosition++;
    }
    if (nodestring.length() == params.maxNodeLen) {
      buildout << currentRefPosition << "," << nodenum << "," << nodestring.c_str() << endl;
      for (int p = 0; p < numUnconnectedPrev; p++) {
        buildout << -2 - p << "," << -1 << endl;
      }
      nodenum++;
      numUnconnectedPrev = 1;
      nodestring = "";
    }
  }
  if (nodestring != "") {
    buildout << currentRefPosition << "," << nodenum << "," << nodestring.c_str() << endl;
    for (int p = 0; p < numUnconnectedPrev; p++) {
      buildout << -2 - p << "," << -1 << endl;
    }
  }
}


void vargas::Graph::exportBuildfile(std::string ref, std::string vcf, std::string build) {
  std::istream *r;
  if (params.inMemoryRef) {
    r = loadFile(ref);
  }
  else {
    r = new std::ifstream(ref);
  }

  vcfstream v(vcf);
  if (!r) throw std::invalid_argument("Error opening file: " + ref);

  if (build.length() > 0) {
    std::ofstream b(build);
    if (!b.good()) throw std::invalid_argument("Error opening file: " + build);
    exportBuildfile(r, v, b);
    b.close();
  } else {
    exportBuildfile(r, v, std::cout);
  }

  delete r;
}

/*********************************** ALIGN *****************************************/

vargas::Alignment *vargas::Graph::align(const vargas::Read &r) {
  Alignment *align = new Alignment;
  Alignment &a = *align;
  vargas::Graph::align(r, a);
  return align;
}


void vargas::Graph::align(const vargas::Read &r, vargas::Alignment &a) {

  int32_t tol = r.read.length();
  gssw_graph_fill(graph, r.read.c_str(), params.nt_table, params.mat,
                  params.gap_open, params.gap_extension, tol, 2, r.readEnd);

  // Absolute alignment positions
  a.optAlignEnd = graph->max_node->data + 1 - graph->max_node->len + graph->max_node->alignment->ref_end;
  a.subOptAlignEnd = graph->submax_node->data + 1 - graph->submax_node->len + graph->submax_node->alignment->ref_end;

  if (r.readEnd > 0) {
    if (r.readEnd > a.optAlignEnd - tol && r.readEnd < a.optAlignEnd + tol) a.corflag = 0;
    else if (r.readEnd > a.subOptAlignEnd - tol && r.readEnd < a.subOptAlignEnd + tol) a.corflag = 1;
    else a.corflag = 2;
  }
  else a.corflag = 2;

  a.optCount = graph->maxCount;
  a.subOptCount = graph->submaxCount;
  a.optScore = graph->max_node->alignment->score;
  a.subOptScore = graph->submax_node->alignment->score;
  a.read = r;
}

/*********************************** TOOLS *****************************************/

void vargas::Graph::generateIngroup(vcfstream &variants) {
  using std::string;
  using std::vector;

  if (params.genComplement) {
    /** Ingroup is the indivs not included in the specified file **/
    if (params.complementSource.length() == 0) {
      throw std::invalid_argument("No buildfile specified, complement cannot be built.");
    }
    std::ifstream complementSource(params.complementSource.c_str());
    string inputGroupLine;
    getline(complementSource, inputGroupLine);
    inputGroupLine = inputGroupLine.substr(1); // Remove #
    vector<string> inputGroup = split(inputGroupLine, ',');
    vector<uint32_t> inputGroupIndivs(0);
    for (uint32_t i = 0; i < inputGroup.size(); i++) {
      inputGroupIndivs.push_back(std::stoul(inputGroup.at(i)));
    }
    variants.createComplementIngroup(inputGroupIndivs);
    complementSource.close();
  }
  else if (params.ingroup >= 0) {
    // Use a percentage of individuals
    variants.createIngroup(params.ingroup);
  }
  else {
    // use all individuals
    variants.createIngroup(100);
  }
}

void vargas::Graph::parseRegion(std::string region, uint32_t *min, uint32_t *max) {
  std::vector<std::string> region_split(0);
  if (region.length() > 0) {
    region_split = split(region, ':');
    if (region_split.size() < 2) {
      std::cerr << "Malformed region, must be in the form a:b" << std::endl;
      *min = 0;
      *max = UINT32_MAX;
    } else {
      *min = (uint32_t) std::stoi(region_split[0]);
      *max = (uint32_t) std::stoi(region_split[1]);
    }
  } else {
    *min = 0;
    *max = UINT32_MAX;
  }
}