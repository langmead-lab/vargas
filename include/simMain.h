//
// Created by gaddra on 9/12/15.
//

#ifndef VMATCH_SIMMAIN_H
#define VMATCH_SIMMAIN_H

#include "getopt_pp.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <regex>
#include "utils.h"

int sim_main(int argc, char *argv[]);


/// <summary>
/// Generates a random read from the given graph. Edges are taken at random.
/// </summary>
/// <param name="graph">Graph to simulate a read from.</param>
/// <param name="readLen">Length of the read to simulate</param>
/// <param name="muterr">Mutation error rate</param>
/// <param name="indelerr">Indel error rate</param>
/// <returns>Constructed gssw_graph</returns>
std::string generateRead(gssw_graph &graph, int32_t readLen, float muterr, float indelerr);


void printSimHelp();


#endif //VMATCH_SIMMAIN_H
