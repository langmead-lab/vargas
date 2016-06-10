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

#ifndef VARGAS_MAIN_H
#define VARGAS_MAIN_H

// Operational modes
int build_main(const int argc, const char *argv[]);
int export_main(const int argc, const char *argv[]);
//int sim_main(const int argc, const char *argv[]);
int align_main(const int argc, const char *argv[]);
int stat_main(const int argc, const char *argv[]);

// Program menus

void printBuildHelp();
void printSimHelp();
void printAlignHelp();
void printExportHelp();
void printStatHelp();

void main_help();
void profile_help();

int profile(const int argc, const char *argv[]);

#endif //VARGAS_MAIN_H
