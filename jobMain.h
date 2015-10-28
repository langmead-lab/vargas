//
// Created by gaddra on 10/2/15.
//

#ifndef VMATCH_JOB_MAIN_H
#define VMATCH_JOB_MAIN_H

int job_main(int argc, char *argv[]);

template<typename t>
std::string getRegex(t subErr = ".*", t numVarNode = ".*", t numVarBase = ".*");

void printJobHelp();

#endif //VMATCH_JOB_MAIN_H
