//
// Created by gaddra on 8/6/15.
//

#ifndef PROJECT_MAIN_H
#define PROJECT_MAIN_H

#define debug 1

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include "gssw/src/gssw.h"
#include "getopt_pp.h"


void printNode(gssw_node *node);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
gssw_graph* buildGraph(std::string buildfile, int8_t *nt_table, int8_t *mat);
gssw_graph* generateGraph(std::string REF, std::string VCF, int8_t *nt_table, int8_t *mat, std::string outputFile="");

#endif //PROJECT_MAIN_H
