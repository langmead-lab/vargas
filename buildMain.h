//
// Created by gaddra on 9/12/15.
//

#ifndef VMATCH_BUILDMAIN_H
#define VMATCH_BUILDMAIN_H

void build_main(int argc, char *argv[]);
void generateGraph(
    std::string REF, std::string VCF,
    int32_t minpos, int32_t maxpos, // Region to parse
    int32_t maxNodeLen,
    int32_t inGroup);
void printBuildHelp();

#endif //VMATCH_BUILDMAIN_H
