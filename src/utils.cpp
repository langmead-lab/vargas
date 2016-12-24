/**
 * @author Ravi Gaddipati
 * @date November 23, 2015
 * rgaddip1@jhu.edu
 *
 * @brief
 * Contains common functions.
 *
 * @file
 */

#include "utils.h"

std::string rg::current_date() {
    time_t t = time(0);
    struct tm *now = localtime(&t);
    std::ostringstream ss;
    ss << (now->tm_year + 1900) << '-'
       << (now->tm_mon + 1) << '-'
       << now->tm_mday;
    return ss.str();
}

std::vector<std::string> rg::split(const std::string &str, const std::string &delims, bool skip_empty) {
    std::vector<std::string> output;
    for_each_token(str.cbegin(), str.cend(), delims.cbegin(), delims.cend(),
                   [&](std::string::const_iterator first, std::string::const_iterator second) {
                       if (first != second || !skip_empty) {
                           output.emplace_back(first, second);
                       }
                   });
    return output;
}


char rg::guess_delim(const std::string &line) {
    static const std::string options = "\n\t:;,=|/";
    for (const char d : options) {
        if (split(line, d).size() > 1) return d;
    }
    throw std::logic_error("Unable to determine delimiter in line: " + line);
}

#if RG_UTIL_INCLUDE_DOCTESET
#include "doctest.h"


TEST_CASE ("Split strings") {
    std::string s("a b\tccc  d "), delim(" \t");
    {
        auto o = rg::split(s, delim);
        REQUIRE(o.size() == 4);
        CHECK(o[0] == "a");
        CHECK(o[1] == "b");
        CHECK(o[2] == "ccc");
        CHECK(o[3] == "d");
    }

    {
        delim = " ";
        auto o = rg::split(s, delim);
        REQUIRE(o.size() == 3);
        CHECK(o[0] == "a");
        CHECK(o[1] == "b\tccc");
        CHECK(o[2] == "d");

    }

    {
        delim = " \t";
        auto o = rg::split(s, delim, false);
        REQUIRE(o.size() == 5);
        CHECK(o[0] == "a");
        CHECK(o[1] == "b");
        CHECK(o[2] == "ccc");
        CHECK(o[3] == "");
        CHECK(o[4] == "d");
    }

    {
        delim = " ";
        auto o = rg::split(s, delim, false);
        REQUIRE(o.size() == 4);
        CHECK(o[0] == "a");
        CHECK(o[1] == "b\tccc");
        CHECK(o[2] == "");
        CHECK(o[3] == "d");

    }

}

#endif