//
// Created by gaddra on 11/12/15.
//

#include "../include/graph.h"


void vmatch::Graph::exportDOT() {
  using std::cout;
  using std::cin;

  if (graph == NULL) {
    std::cerr << "Error: No graph has been build. Aborting export." << std::endl;
    return;
  }

  cout << "digraph gssw_graph {\n";
  cout << "rankdir=\"LR\";\n";

  for (uint32_t i = 0; i < graph->size; i++) {
    cout << graph->nodes[i]->id << " [label=\"" << graph->nodes[i]->data << ":" << graph->nodes[i]->seq << "\"];\n";
  }
  for (uint32_t i = 0; i < graph->size; i++) {
    for (int32_t n = 0; n < graph->nodes[i]->count_next; n++)
      cout << graph->nodes[i]->id << " -> " << graph->nodes[i]->next[n]->id << ";\n";
  }
  cout << "}";
}

void vmatch::Graph::exportBuildfile(std::ofstream &out) {

}

void vmatch::Graph::buildGraph(std::istream &graphDat) {

  using std::string;
  using std::endl;
  using std::vector;
  using std::cerr;

  if (graph != NULL) {
    delete graph;
  }

  string line;
  vector<string> lineSplit(0);
  vector<gssw_node *> nodes(0);
  uint32_t curr = 0;

  /** Build nodes and edges from buildfile **/
  cerr << "Building graph..." << endl;
  while (getline(graphDat, line)) {
    if (line.at(0) != '#') {
      if (line.at(0) == '[') {
        line = line.substr(1, line.length() - 2);
        lineSplit = split(line, ',');
        for (uint32_t i = 0; i < lineSplit.size(); i++) {
          gssw_node_add_indiv(nodes.back(), strtol(lineSplit[i].c_str(), NULL, 10));
        }
      } else {
        lineSplit = split(line, ',');
        switch (lineSplit.size()) {
          case 3: // New node
            curr = uint32_t(strtol(lineSplit[1].c_str(), NULL, 10));
            nodes.push_back(gssw_node_create(int32_t(strtol(lineSplit[0].c_str(), NULL, 10)),
                                             curr,
                                             lineSplit[2].c_str(), params.nt_table, params.mat));
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

  /** Add nodes to graph **/
  graph = gssw_graph_create(uint32_t(nodes.size()));
  for (uint32_t n = 0; n < nodes.size(); n++) {
    gssw_graph_add_node(graph, nodes[n]);
  }

}

std::iostream &vmatch::Graph::buildGraph(std::istream &reference, std::istream &variants) {
  using std::vector;
  using std::endl;
  using std::cerr;
  using std::string;

  std::stringstream *cout_p = new std::stringstream(); // Output stream
  std::stringstream &cout = *cout_p;
  std::ifstream build; // Used for complement graph

  /** Used to track which indivs have a variant **/
  vector<int16_t> inVar(0);

  /** To process complement graph **/
  vector<string> inputGroup(0);
  vector<int32_t> inputGroupInt(0);
  string inputGroupLine;

  /** reference line, variant line **/
  string ref_line, vcf_line;
  /** Parsed lines, VCF header row **/
  vector<string> vline_split(0), altList_split(0), header(0), info_split(0);
  /** Columns used to build graph **/
  vector<int32_t> inGroupCols(0);
  /** columns in VCF file, current ref pos **/
  int32_t posColumn, refColumn, altColumn, formatColumn, infoColumn, ref_position = 0;
  int32_t nodenum = 0, numIndivs;
  /** Pos from variant file **/
  int vpos;
  int32_t cn;
  /** To track edges that need to be built **/
  int numalts = 0, numprev;
  /** strings that represent node contents **/
  string variantRef, variantAlt, nodestring, info;
  int32_t nodelen = 0;
  char base;
  int32_t randtemp;

  double altAFsum;
  double maxAltAF;
  int32_t maxAltAFIdx;
  double af;

  // Get region
  uint32_t minpos, maxpos;
  parseRegion(params.region, &minpos, &maxpos);

  getline(reference, ref_line);
  if (ref_line.at(0) != '>') cerr << "Error in ref file, first char should be >" << endl;

  /** Go to first VCF record **/
  do { getline(variants, vcf_line); } while (vcf_line.substr(0, 2) == "##");
  transform(vcf_line.begin(), vcf_line.end(), vcf_line.begin(), ::tolower);
  header = split(vcf_line, '\t');
  posColumn = int32_t(find(header.begin(), header.end(), "pos") - header.begin());
  refColumn = int32_t(find(header.begin(), header.end(), "ref") - header.begin());
  altColumn = int32_t(find(header.begin(), header.end(), "alt") - header.begin());
  infoColumn = int32_t(find(header.begin(), header.end(), "info") - header.begin());
  formatColumn = int32_t(find(header.begin(), header.end(), "format") - header.begin());
  numIndivs = int32_t(header.size() - formatColumn - 1);

  /** Construct the in group, the graph will be built with these individuals **/
  cout << '#';
  if (params.genComplement) {
    /** Ingroup is the indivs not included in the specified file **/
    if (params.buildfile.length() == 0) {
      cerr << "Error: No buildfile specified, complement cannot be built. Aborting." << endl;
      exit(1);
    }
    build.open(params.buildfile.c_str());
    getline(build, inputGroupLine);
    inputGroupLine = inputGroupLine.substr(1, inputGroupLine.length() - 2);
    inputGroup = split(inputGroupLine, ',');
    for (uint32_t i = 0; i < inputGroup.size(); i++) {
      inputGroupInt.push_back(atoi(inputGroup.at(i).c_str()));
    }
    for (int32_t i = 0; i < numIndivs; i++) {
      if (find(inputGroupInt.begin(), inputGroupInt.end(), i + formatColumn + 1) == inputGroupInt.end()) {
        inGroupCols.push_back(i + formatColumn + 1);
        cout << inGroupCols[i] << ',';
      }
    }

  } else {
    if (params.ingroup >= 0) {
      for (int32_t i = 0; i < int32_t((numIndivs / 100.0f) * params.ingroup); i++) {
        randtemp = rand() % numIndivs + formatColumn + 1;
        if (find(inGroupCols.begin(), inGroupCols.end(), randtemp) == inGroupCols.end()) {
          inGroupCols.push_back(randtemp);
          cout << inGroupCols[i] << ',';
        } else i--;
      }
      sort(inGroupCols.begin(), inGroupCols.end());
    } else {
      for (int32_t i = 0; i < numIndivs; i++) {
        inGroupCols.push_back(i + formatColumn + 1);
        cout << inGroupCols[i] << ',';
      }
    }
  }
  cout << endl;

  /** Find the POS, REF, ALT cols **/
  if (posColumn < 0 || refColumn < 0 || altColumn < 0 || formatColumn < 0) {
    cerr << "POS, REF, ALT, FORMAT, and/or INFO not found in VCF header." << endl;
    exit(1);
  }

  /** Generate Nodes **/

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
    vline_split = split(vcf_line, '\t');
    vpos = atoi(vline_split[posColumn].c_str());
    if (vpos <= ref_position) continue;
    if (vpos > maxpos) break;
    variantRef = vline_split[refColumn];
    variantAlt = vline_split[altColumn];

    // Get list of allele frequencies
    if (params.maxAF) {
      info = vline_split[infoColumn];
      info_split = split(info, ';');
      maxAltAF = 0;
      altAFsum = 0;
      maxAltAFIdx = -1;
      params.maxAF = false; // Using as a flag to make sure 'AF' is found in INFO
      for (auto infItem : info_split) {
        std::transform(infItem.begin(), infItem.end(), infItem.begin(), ::tolower);
        if (infItem.substr(0, 2) == "af") {
          // Find the list of AF's (allele frequencies)
          info = infItem.substr(3, string::npos);
          params.maxAF = true;
          break;
        }
      }
      if (!params.maxAF) {
        std::cerr << "Error: AF not found in INFO field." << std::endl;
        std::cerr << "At variant position " << vline_split[1] << std::endl;
        exit(1);
      }
      // Get the raw list of AF's
      info_split = split(info, ',');
      // Find the max AF and its index, as well as the sum of all alternate AF's
      for (uint16_t i = 0; i < info_split.size(); i++) {
        af = stod(info_split[i]);
        if (af > maxAltAF) {
          maxAltAF = af;
          maxAltAFIdx = i;
        }
        altAFsum += af;
      }
      if (maxAltAFIdx < 0) {
        std::cerr << "Error: max AF not found." << std::endl;
        std::cerr << "At VCF pos " << vline_split[1] << std::endl;
        std::cerr << "AF list: " << info << std::endl;
        std::cerr << "Continuing with reference node." << std::endl;
      }
    }


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
      if (nodelen == params.maxNodeLen) {
        cout << ref_position << "," << nodenum << "," << nodestring.c_str() << endl;
#if debug > 4
        cerr << "Node: " << ref_position << ", ID: " << nodenum << ", " << nodestring << endl;
#endif
        if (nodenum != 0) {
          /** Connect to all of the previous alt/ref nodes **/
          for (int i = 0; i < numalts; i++) {
            cout << -2 - i << "," << -1 << endl;
#if debug > 4
            cerr << "Edge: " << nodes.end()[-2 - i]->id << ", " << nodes.end()[-1]->id << endl;
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
      cout << ref_position << "," << nodenum << "," << nodestring.c_str() << endl;
#if debug > 4
      cerr << "Node: " << ref_position << ", ID: " << nodenum << ", " << nodestring << endl;
#endif
      /** Only connect with edge if it's not the first node **/
      if (nodenum != 0) {
        /** Connect to all of the previous alt/ref nodes **/
        for (int i = 0; i < numalts; i++) {
          cout << -2 - i << "," << -1 << endl;
#if debug > 4
          cerr << "Edge: " << nodes.end()[-2 - i]->id << ", " << nodes.end()[-1]->id << endl;
#endif
        }
      }
      nodenum++;
      numprev = 1;
    }
    else numprev = numalts;

    /** Ref node **/
    for (uint32_t i = 0; i < variantRef.length(); i++) {
      reference.get(base);
      if (isspace(base)) {
        reference.get(base);
      }
      ref_position++;
    }
    numalts = 0;

    /** Variants **/
    if (params.maxAF) {
      // Only use the max AF variant
      altList_split.clear();
      // The reference has the highest AF, so add a ref node (moved below)
      if (1 - altAFsum >= maxAltAF || maxAltAFIdx < 0);
        // an Alt is the highest AF, so add that to the cleared var list
      else altList_split.push_back(split(variantAlt, ',')[maxAltAFIdx]);
    } else {
      altList_split = split(variantAlt, ',');
    }
    // Add a single node if not doing max AF or the max AF is the ref
    if (!params.maxAF || (params.maxAF && (1 - altAFsum >= maxAltAF || maxAltAFIdx < 0))) {
      cout << ref_position << "," << nodenum << "," << variantRef.c_str() << endl;
#if debug > 4
      cerr << "Node: " << ref_position << ", ID: " << nodenum << ", " << variantRef << endl;
#endif
      nodenum++;
      numalts = 1;
    }

    for (uint32_t i = 0; i < altList_split.size(); i++) {
      inVar.clear();
      for (uint32_t c = 0; c < inGroupCols.size(); c++) {
        /** Check if it is in the ingroup **/
        //TODO parse rather than at, check diploid
        if (strtol(&vline_split[inGroupCols[c]].at(0), NULL, 10) == i + 1) {
          inVar.push_back(inGroupCols[c]);
        }
      }

      // if there are variants, add node otherwise continue (i.e. ref is merged with previous node)
      if (inVar.size() > 0) {

        if (altList_split[i].substr(0, 3) == "<CN") {
          cn = int32_t(strtol(altList_split[i].substr(3, altList_split[i].length() - 4).c_str(), NULL, 10));
          altList_split[i] = "";
          for (int v = 0; v < cn; v++) {
            altList_split[i] += variantRef;
          }
        }
        cout << ref_position << "," << nodenum << "," << altList_split[i].c_str() << endl;
#if debug > 4
        cerr << "Node: " << ref_position << ", ID: " << nodenum << ", " << altList_split[i].c_str() << endl;
#endif
        nodenum++;
        numalts++;

        /** Add the individuals that have this particular variant **/
        cout << '[';
        for (uint32_t n = 0; n < inVar.size(); n++) {
#if debug > 5
          cerr << "Add Indiv(" << int32_t(nodes.back()->indivSize) << "): " << inVar[n] << endl;
#endif
          cout << inVar[n];
          if (n != inVar.size() - 1) cout << ',';
        }
        cout << ']' << endl;

      }
    }

    /** Build edges **/
    for (int p = 0; p < numprev; p++) {
      for (int a = 0; a < numalts; a++) {
        cout << -1 - numalts - p << "," << -1 - a << endl;
#if debug > 4
        cerr << "Edge: " << nodes.end()[-1 - numalts - p]->id << ", " << nodes.end()[-1 - a]->id << endl;
#endif
      }
    }
  }

  /** The remaining bases after the last variant **/
  nodestring = "";
  while ((ref_position < maxpos) && reference.get(base)) {
    if (!isspace(base)) {
      nodestring += base;
      ref_position++;
    }
  }
  cout << ref_position << "," << nodenum << "," << nodestring.c_str() << endl;
#if debug > 4
  cerr << "Node: " << ref_position << ", ID: " << nodenum << ", " << nodestring.c_str() << endl;
#endif
  for (int p = 0; p < numalts; p++) {
    cout << -2 - p << "," << -1 << endl;
#if debug > 4
    cerr << "Edge: " << nodes.end()[-2 - p]->id << ", " << nodes.end()[-1]->id << endl;
#endif
  }

  return cout;
}

void vmatch::Graph::parseRegion(std::string region, uint32_t *min, uint32_t *max) {
  std::vector<std::string> region_split(0);
  if (region.length() > 0) {
    region_split = split(region, ':');
    if (region_split.size() < 1) {
      std::cerr << "Malformed region, must be in the form a:b" << std::endl;
      *min = 0;
      *max = UINT32_MAX;
    }
    *min = (uint32_t) std::atoi(region_split[0].c_str());
    *max = (uint32_t) std::atoi(region_split[1].c_str());
  } else {
    *min = 0;
    *max = UINT32_MAX;
  }
}