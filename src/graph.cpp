/**
 * Ravi Gaddipati
 * November 20, 2015
 * rgaddip1@jhu.edu
 *
 * vargas::Graph is a DAG representation of a reference and its variants.
 * The class wraps a gssw_graph from gssw and provides a way to construct
 * graphs from a FASTA and a VCF file with various options.
 *
 * GSSW was originally written by Erik Garrison and was moderately modified.
 *
 * graph.cpp
 */

#include "../include/graph.h"

void vargas::Graph::exportDOT(std::ostream &out) const {

  if (graph == NULL) {
    std::cerr << "Error: No graph has been built. Aborting export." << std::endl;
    return;
  }

  out << "digraph gssw_graph {\n";
  out << "rankdir=\"LR\";\n";

  for (uint32_t i = 0; i < graph->size; i++) {
    out << graph->nodes[i]->id << " [label=\"" << graph->nodes[i]->data << ":" << graph->nodes[i]->seq << "\"];\n";
  }
  for (uint32_t i = 0; i < graph->size; i++) {
    for (int32_t n = 0; n < graph->nodes[i]->count_next; n++)
      out << graph->nodes[i]->id << " -> " << graph->nodes[i]->next[n]->id << ";\n";
  }
  out << "}";
}

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

  /** Build nodes and edges from buildfile **/
  while (getline(graphDat, line)) {
    if (line.at(0) == '#') continue; // Commment line
      if (line.at(0) == '[') {
        if (params.includeIndividuals) {
          // Add individuals to the last node
          line = line.substr(1, line.length() - 2);
          split(line, ',', lineSplit);
          for (uint32_t i = 0; i < lineSplit.size(); i++) {
            gssw_node_add_indiv(nodes.back(), std::stoi(lineSplit[i].c_str(), NULL, 10));
          }
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
      inputGroupIndivs.push_back(uint32_t(atoi(inputGroup.at(i).c_str())));
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

void vargas::Graph::exportBuildfile(std::istream &reference, vcfstream &variants, std::ostream &buildout) {
  using std::vector;
  using std::endl;
  using std::string;

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
  buildout << "##Reference sequence: " << refLine.substr(1) << endl;
  buildout << "##Graph params: "
      << "Max Node Len: " << params.maxNodeLen
      << ", Use max AF? " << params.maxAF;
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
      if (!reference.get(base)) {
        throw std::range_error("End of ref found while looking for pos " + std::to_string(variantRecord.pos));
      }
      if (!isspace(base)) {
        nodestring += base;
        currentRefPosition++;
      }

      /** Max node length reached, split node **/
      if (nodestring.length() == params.maxNodeLen) {
        buildout << currentRefPosition << "," << nodenum << "," << nodestring.c_str() << endl;
          /** Connect to all of the previous alt/ref nodes **/
        for (int i = 0; i < numUnconnectedPrev; ++i) {
          buildout << -2 - i << "," << -1 << endl;
          }
        numUnconnectedPrev = 1; // We have one node at the tail of the graph
        nodenum++;
        nodestring = "";
      }
    }

    /** If there is space between the variants, add a new node **/
    if (nodestring.length() > 0) {
      buildout << currentRefPosition << "," << nodenum << "," << nodestring.c_str() << endl;
        /** Connect to all of the previous alt/ref nodes **/
      for (int i = 0; i < numUnconnectedPrev; ++i) {
        buildout << -2 - i << "," << -1 << endl;
        }
      numUnconnectedPrev = 1;
      nodenum++;
      nodestring = "";
    }

    std::string faref = "";
    /** progress reference position past the ref that's stored in the variant file **/
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
      printVec<uint32_t>(buildout, variantRecord.indivs[maxAFAllele]);

    } else {
      for (auto &var : variantRecord.indivs) {
        buildout << currentRefPosition << "," << nodenum << "," << var.first << endl;
        nodenum++;
        numUnconnectedCurr++;
        // Print individuals that have this variant
        printVec<uint32_t>(buildout, var.second);
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

void vargas::Graph::exportBuildfile(std::string ref, std::string vcf, std::string build) {
  std::ifstream r(ref);
  vcfstream v(vcf);
  if (!r.good()) throw std::invalid_argument("Error opening file: " + ref);

  if (build.length() > 0) {
    std::ofstream b(build);
    if (!b.good()) throw std::invalid_argument("Error opening file: " + build);
    exportBuildfile(r, v, b);
    b.close();
  } else {
    exportBuildfile(r, v, std::cout);
  }
  r.close();
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
      *min = (uint32_t) std::stoi(region_split[0].c_str());
      *max = (uint32_t) std::atoi(region_split[1].c_str());
    }
  } else {
    *min = 0;
    *max = UINT32_MAX;
  }
}