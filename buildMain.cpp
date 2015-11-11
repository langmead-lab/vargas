//
// Created by gaddra on 9/12/15.
//

#include "buildMain.h"


int build_main(int argc, char *argv[]) {
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::string;

  string VCF = "", REF = "", region, buildfile = "";
  int32_t maxNodelen = 50000, ingroup = -1;
  bool genComplement = false;
  bool set = false;
  bool maxAF = false; // Max allele frequency graph
  std::string setString;
  std::vector<std::string> setSplit(0);

  /**default region **/
  int32_t regionMin = 0, regionMax = 2147483640;
  std::vector<string> region_split(0);

  GetOpt::GetOpt_pp args(argc, argv);

  if (args >> GetOpt::OptionPresent('h', "help")) {
    printBuildHelp();
    return 0;
  }

  if (!(args >> GetOpt::Option('v', "vcf", VCF))
      || !(args >> GetOpt::Option('r', "ref", REF))) {
    cerr << "No inputs specified!" << endl;
    return 1;
  }

  if ((args >> GetOpt::Option('s', "set", setString))) {
    set = true;
    setSplit = split(setString, ',');
  }

  args >> GetOpt::Option('l', "maxlen", maxNodelen)
      >> GetOpt::Option('R', "region", region)
      >> GetOpt::Option('g', "ingroup", ingroup)
      >> GetOpt::OptionPresent('c', "complement", genComplement)
      >> GetOpt::Option('c', "complement", buildfile)
      >> GetOpt::OptionPresent('m', "maxref", maxAF);

  /** Parse region **/
  if (region.length() > 0) {
    region_split = split(region, ':');
    if (region_split.size() < 1) {
      std::cerr << "Malformed region, must be in the form a:b" << endl;
      return 1;
    }
    regionMin = std::atoi(region_split[0].c_str());
    regionMax = std::atoi(region_split[1].c_str());
  }
  // max AF ref overrides other opts
  if (maxAF) {
    set = false;
    genComplement = false;
    ingroup = -1;
  }
  if (set) {
    // Gen a set of ingroup build files
    std::ofstream out(0);
    auto oldBuf = std::cout.rdbuf(out.rdbuf());
    std::stringstream fileName, compFilename;
    for (auto ingrp : setSplit) {
      // Generates a buildfile for all % ingroup specified, as well as the outgroup buildfiles
      fileName << ingrp << "In.build";
      compFilename << ingrp << "Out.build";
      out.open(fileName.str());
      generateGraph(REF, VCF, regionMin, regionMax, maxNodelen, stoi(ingrp), false, "", false);
      out.close();
      out.open(compFilename.str());
      generateGraph(REF, VCF, regionMin, regionMax, maxNodelen, stoi(ingrp), true, fileName.str(), false);
      out.close();
      fileName.str("");
      compFilename.str("");
    }
    std::cout.rdbuf(oldBuf);

  } else generateGraph(REF, VCF, regionMin, regionMax, maxNodelen, ingroup, genComplement, buildfile, maxAF);

  return 0;
}

void generateGraph(
    std::string REF, std::string VCF,
    int32_t minpos, int32_t maxpos, // Region to parse
    int32_t maxNodeLen,
    int32_t inGroup,
    bool genComplement,
    std::string buildfile,
    bool maxAF) {

  using namespace std;

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
  header = split(vcf_line, '\t');
  posColumn = int32_t(find(header.begin(), header.end(), "pos") - header.begin());
  refColumn = int32_t(find(header.begin(), header.end(), "ref") - header.begin());
  altColumn = int32_t(find(header.begin(), header.end(), "alt") - header.begin());
  infoColumn = int32_t(find(header.begin(), header.end(), "info") - header.begin());
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
    inputGroup = split(inputGroupLine, ',');
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
    if (maxAF) {
      info = vline_split[infoColumn];
      info_split = split(info, ';');
      maxAltAF = 0;
      altAFsum = 0;
      maxAltAFIdx = -1;
      maxAF = false; // Using as a flag to make sure 'AF' is found in INFO
      for (auto infItem : info_split) {
        std::transform(infItem.begin(), infItem.end(), infItem.begin(), ::tolower);
        if (infItem.substr(0, 2) == "af") {
          // Find the list of AF's (allele frequencies)
          info = infItem.substr(3, string::npos);
          maxAF = true;
          break;
        }
      }
      if (!maxAF) {
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
      // The reference has the highest AF, so go to the next VCF line. As a result, the reference is used.
      if (1 - altAFsum >= maxAltAF || maxAltAFIdx < 0) continue;
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
    numalts = 0;
    if (!maxAF) {
      cout << ref_position << "," << nodenum << "," << variantRef.c_str() << endl;
#if debug > 4
      cerr << "Node: " << ref_position << ", ID: " << nodenum << ", " << variantRef << endl;
#endif
      nodenum++;
      numalts = 1;
    }
    /** Variants **/
    if (maxAF) {
      // Only use the max AF variant
      altList_split.clear();
      altList_split.push_back(split(variantAlt, ',')[maxAltAFIdx]);
    } else {
      altList_split = split(variantAlt, ',');
    }
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
  cout << "-c\t--complement    Generate a complement of the specified graph" << endl;
  cout << "-s\t--set           <#,#,..,#> Generate a buildfile for a list of ingroup %'s and their complements."
      << endl;
  cout << "-m\t--maxref        Generate a graph using allele's w/ the highest frequency. Overrides other opts." << endl;
  cout << "\t                  -s outputs to files." << endl;

  cout << endl << "Buildfile is printed on stdout." << endl;
}