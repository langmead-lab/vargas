/**
 * Ravi Gaddipati
 * Dec 24, 2015
 * rgaddip1@jhu.edu
 *
 * Loads a file into a string.
 *
 * loadfile.cpp
 */


#include <stdexcept>
#include <sstream>
#include <fstream>
#include <vector>

std::stringstream *loadFile(const std::string filename) {

    std::ifstream file(filename, std::ios::in | std::ios::binary);

    if (!file) {
        throw std::invalid_argument("Invalid file: " + std::string(filename));
    }


    file.seekg(0, std::ios::end);
    std::streampos length = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(length);
    file.read(&buffer[0], length);

    std::stringstream *localStream = new std::stringstream;
    (*localStream).rdbuf()->pubsetbuf(&buffer[0], length);

    return localStream;
}
