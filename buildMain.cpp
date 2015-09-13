//
// Created by gaddra on 9/12/15.
//

#include "buildMain.h"


void build_main(int argc, char *argv[]) {
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::string;

  string VCF = "", REF = "", region, buildfile = "";
  int32_t maxNodelen = 50000, ingroup = -1;
  bool genComplement = false;

  /**default region **/
  int32_t regionMin = 0, regionMax = 2147483640;
  std::vector<string> region_split(0);

  GetOpt::GetOpt_pp args(argc, argv);

  if (args >> GetOpt::OptionPresent('h', "help")) {
    printBuildHelp();
    exit(0);
  }

  if (!(args >> GetOpt::Option('v', "vcf", VCF))
      || !(args >> GetOpt::Option('r', "ref", REF))) {
    cerr << "No inputs specified!" << endl;
    exit(1);
  }

  args >> GetOpt::Option('l', "maxlen", maxNodelen)
      >> GetOpt::Option('R', "region", region)
      >> GetOpt::Option('g', "ingroup", ingroup)
      >> GetOpt::OptionPresent('c', "complement", genComplement)
      >> GetOpt::Option('b', "buildfile", buildfile);

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
  generateGraph(REF, VCF, regionMin, regionMax, maxNodelen, ingroup, genComplement, buildfile);
}

void generateGraph(
    std::string REF, std::string VCF,
    int32_t minpos, int32_t maxpos, // Region to parse
    int32_t maxNodeLen,
    int32_t inGroup,
    bool genComplement,
    std::string buildfile) {

  using namespace std;

  /** Used to track which indivs have a variant **/
  vector<int16_t> inVar(0);

  /** To process complement graph **/
  vector<string> inputGroup(0);
  vector<int32_t> inputGroupInt(0);
  string inputGroupLine;
  int32_t ingroupSize;

  /** reference line, variant line **/
  string ref_line, vcf_line;
  /** Parsed lines, VCF header row **/
  vector<string> vline_split(0), altList_split(0), header(0);
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
  int32_t randtemp;

  /** File stream **/
  ifstream variants(VCF.c_str(), ios_base::in | ios_base::binary);
  ifstream reference(REF.c_str());
  ifstream build;

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
  cout << '#';
  if (genComplement) {
    /** Ingroup is the indivs not included in the specified file **/
    if (buildfile.length() == 0) {
      cerr << "Error: No buildfile specified, complement cannot be built. Aborting." << endl;
      exit(1);
    }
    build.open(buildfile.c_str());
    getline(build, inputGroupLine);
    inputGroupLine = inputGroupLine.substr(1, inputGroupLine.length() - 2);
    split(inputGroupLine, ',', inputGroup);
    for (int32_t i = 0; i < inputGroup.size(); i++) {
      inputGroupInt.push_back(atoi(inputGroup.at(i).c_str()));
    }
    for (int32_t i = 0; i < numIndivs; i++) {
      if (find(inputGroupInt.begin(), inputGroupInt.end(), i + formatColumn + 1) == inputGroupInt.end()) {
        inGroupCols.push_back(i + formatColumn + 1);
        cout << inGroupCols[i] << ',';
      }
    }

  } else {
    if (inGroup >= 0) {
      for (int32_t i = 0; i < int32_t((numIndivs / 100.0f) * inGroup); i++) {
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
    cerr << "POS, REF, ALT, and/or INFO not found in VCF header." << endl;
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
    split(vcf_line, '\t', vline_split);
    vpos = atoi(vline_split[posColumn].c_str());
    if (vpos <= ref_position) goto endvar;
    if (vpos > maxpos) break;
    variantRef = vline_split[refColumn];
    variantAlt = vline_split[altColumn];

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
    for (int i = 0; i < variantRef.length(); i++) {
      reference.get(base);
      if (isspace(base)) {
        reference.get(base);
      }
      ref_position++;
    }
    cout << ref_position << "," << nodenum << "," << variantRef.c_str() << endl;
#if debug > 4
    cerr << "Node: " << ref_position << ", ID: " << nodenum << ", " << variantRef << endl;
#endif
    nodenum++;
    numalts = 1;

    /** Variants **/
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
        for (int n = 0; n < inVar.size(); n++) {
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

  variants.close();
  reference.close();

}

void printBuildHelp() {
  using std::cout;
  using std::endl;
  cout << endl << "------------------- VMatch build, " << __DATE__ << ". rgaddip1@jhu.edu -------------------" <<
      endl;
  cout << "-v\t--vcf           (required) VCF file, uncompressed." << endl;
  cout << "-r\t--ref           (required) reference single record FASTA" << endl;
  cout << "-l\t--maxlen        Maximum node length" << endl;
  cout << "-R\t--region        [min:max] Ref region, inclusive. Default is entire graph." << endl;
  cout << "-g\t--ingroup       Percent of individuals to build graph from, default all." << endl;
  cout << "-c\t--complement    Generate a complement of the graph specified with -b." << endl;
  cout << "-b\t--buildfile     Buildfile of the graph to generate a complement of." << endl;

  cout << endl << "Buildfile is printed on stdout." << endl;
}