/**
 * Ravi Gaddipati
 * November 25, 2015
 * rgaddip1@jhu.edu
 *
 * Interface for simulating and aligning reads from/to a DAG.
 * Uses a modified gssw from Erik Garrison.
 *
 * main.h
 */

#ifndef VARGAS_H
#define VARGAS_H

#include <string>
#include <stdlib.h>
#include <iostream>
#include "getopt_pp.h"
#include "graph.h"
#include "readsim.h"
#include "readfile.h"


int build_main(const int argc, const char *argv[]);
int export_main(const int argc, const char *argv[]);
int sim_main(const int argc, const char *argv[]);
int align_main(const int argc, const char *argv[]);

void printMainHelp();
void printBuildHelp();
void printSimHelp();
void printAlignHelp();
void printExportHelp();

#endif //VARGAS_H
