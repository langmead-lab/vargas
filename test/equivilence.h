//
// Created by gaddra on 11/25/15.
//

#ifndef VARGAS_EQUIVILANCE_H
#define VARGAS_EQUIVILANCE_H

#include <string>
#include <fstream>
#include <algorithm>
#include <set>
#include <iostream>

bool fileEquiv(std::string afile, std::string bfile, bool orderImp = false);
bool fileEquiv(std::istream &a, std::istream &b, bool orderImp = false);
bool fileEquiv(std::string afile, std::istream &b, bool orderImp = false);
bool fileEquiv(std::istream &a, std::string bfile, bool orderImp = false);

#endif //VARGAS_EQUIVILANCE_H
