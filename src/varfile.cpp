/**
 * Ravi Gaddipati
 * November 23, 2015
 * rgaddip1@jhu.edu
 *
 * @brief
 * Provides tools to interact with VCF/BCF files.
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#include "varfile.h"

vargas::Region vargas::parse_region(const std::string &region_str) {
    vargas::Region ret;

    std::vector <std::string> regionSplit = rg::split(region_str, ':');

    if(regionSplit.size() == 1) {
        regionSplit.push_back("0-0");
    }
    else if (regionSplit.size() != 2) throw std::invalid_argument("Invalid region: " + region_str);

    ret.seq_name = regionSplit[0];

    // Strip commas
    regionSplit[1].erase(std::remove(regionSplit[1].begin(), regionSplit[1].end(), ','), regionSplit[1].end());

    // Range
    regionSplit = rg::split(regionSplit[1], '-');
    if (regionSplit.size() != 2)
        throw std::invalid_argument("Invalid region format, should be CHR:XX,XXX-YY,YYY\n\t" + region_str);

    ret.min = std::stoul(regionSplit[0]);
    ret.max = std::stoul(regionSplit[1]);

    if (ret.min > ret.max) {
        throw std::invalid_argument("Invalid region, min > max.");
    }

    return ret;
}

void vargas::VCF::set_region(const Region &region) {
    _region = region;
}


std::vector<std::string> vargas::VCF::sequences() const {
    if (!_header) return std::vector<std::string>(0);
    std::vector<std::string> ret;
    int num;
    // As per htslib, only top level is freed.
    const char **seqnames = bcf_hdr_seqnames(_header, &num);
    if (!seqnames) {
        throw std::invalid_argument("Error reading VCF header.");
    }
    for (int i = 0; i < num; ++i) {
        ret.push_back(std::string(seqnames[i]));
    }
    free(seqnames);
    return ret;
}


bool vargas::VCF::next() {
    if (_limit > 0 && _counter >= _limit) return false;
    if (!_header || !_bcf) return false;
    do {
        if (bcf_read(_bcf, _header, _curr_rec) != 0) return false;
    } while ((_region.seq_name.size() && strcmp(_region.seq_name.c_str(), bcf_seqname(_header, _curr_rec))) ||
        unsigned(_curr_rec->pos) < _region.min ||
        (_region.max > 0 && unsigned(_curr_rec->pos) > _region.max));

    unpack_all();
    gen_genotypes();
    ++_counter;
    return true;
}


const std::vector<std::string> &vargas::VCF::gen_genotypes() {
    FormatField<int> gt(_header, _curr_rec, "GT");
    _genotypes.resize(gt.values.size());
    for (size_t i = 0; i < _genotypes.size(); ++i) {
        _genotypes[i] = (_alleles[bcf_gt_allele(gt.values.at(i))]);
    }

    // Map of which indivs have each allele
    _genotype_indivs.clear();
    for (auto &allele : alleles()) {
        _genotype_indivs[allele] = Population(_genotypes.size(), false);
    }
    for (size_t s = 0; s < _genotypes.size(); ++s) {
        _genotype_indivs[_genotypes[s]].set(s);
    }

    return _genotypes;
}


const std::vector<float> &vargas::VCF::frequencies() const {
    InfoField<float> af(_header, _curr_rec, "AF");
    static std::vector<float> _allele_freqs;
    const auto &val = af.values;
    _allele_freqs.resize(val.size() + 1); // make room for the ref
    float sum = 0;

    // Get the ref frequency and add to result vector +1, so we can put ref at index 0
    for (size_t i = 0; i < val.size(); ++i) {
        sum += val[i];
        _allele_freqs[i + 1] = val[i];
    }
    _allele_freqs[0] = 1 - sum;

    return _allele_freqs;
}


void vargas::VCF::create_ingroup(int percent) {
    _ingroup.clear();

    if (percent == 100) {
        _ingroup = _samples;
    } else if (percent != 0) {
        for (auto s : _samples) {
            if (rand() % 100 < percent) _ingroup.push_back(s);
        }
    }

    _apply_ingroup_filter();
}


int vargas::VCF::_init() {
    _counter = 0;
    _limit = 0;
    if (_file_name.length() && _file_name != "-") {
        _bcf = bcf_open(_file_name.c_str(), "r");
        if (!_bcf) return -1;

        _header = bcf_hdr_read(_bcf);
        if (_header == nullptr) {
            return -2;
        }

        // Load samples
        for (size_t i = 0; i < num_haplotypes() / 2; ++i) {
            _samples.push_back(_header->samples[i]);
        }
        create_ingroup(100);
    }
    return 0;
}


void vargas::VCF::_load_shared() {
    _alleles.clear();
    for (int i = 0; i < _curr_rec->n_allele; ++i) {
        std::string allele(_curr_rec->d.allele[i]);
        // Some replacement tag
        if (allele.at(0) == '<') {
            std::string ref = _curr_rec->d.allele[0];
            // Copy number
            if (allele.substr(1, 2) == "CN" && allele.at(3) != 'V') {
                int copy = std::stoi(allele.substr(3, allele.length() - 4));
                allele = "";
                for (int i = 0; i < copy; ++i) allele += ref;
            } else {
                // Other types are just subbed with the ref.
                allele = ref;
            }
        }
        _alleles.push_back(allele);
    }
}


void vargas::VCF::_apply_ingroup_filter() {
    if (_header == nullptr) {
        throw std::logic_error("Ingroup filter should only be applied after loading header!");
    }

    if (_ingroup.size() == 0) {
        if (_ingroup_cstr != nullptr) {
            free(_ingroup_cstr);
        }
        _ingroup_cstr = nullptr;
    } else if (_ingroup.size() == _samples.size()) {
        if (_ingroup_cstr != nullptr) {
            free(_ingroup_cstr);
        }
        _ingroup_cstr = (char *) malloc(2);
        strcpy(_ingroup_cstr, "-");
    } else {
        std::ostringstream ss;
        for (auto s : _ingroup) ss << s << ',';
        std::string smps = ss.str().substr(0, ss.str().length() - 1);
        if (_ingroup_cstr != nullptr) {
            free(_ingroup_cstr);
        }
        _ingroup_cstr = (char *) malloc(smps.length() + 1);
        strcpy(_ingroup_cstr, smps.c_str());
    }

    bcf_hdr_set_samples(_header, _ingroup_cstr, 0);
}

size_t vargas::VCF::num_haplotypes() const {
    if (_header == nullptr) {
        return 0;
    }
    return (size_t) bcf_hdr_nsamples(_header) * 2;
}

void vargas::VCF::close() {
    if (_bcf) bcf_close(_bcf);
    if (_header != nullptr) {
        bcf_hdr_destroy(_header);
    }
    if (_curr_rec != nullptr) {
        bcf_destroy(_curr_rec);
    }
    if (_ingroup_cstr != nullptr) {
        free(_ingroup_cstr);
    }
    _bcf = nullptr;
    _header = nullptr;
    _curr_rec = nullptr;
    _ingroup_cstr = nullptr;
}

TEST_SUITE("VCF Parser");

TEST_CASE ("VCF File handler") {
    using std::endl;
    std::string tmpvcf = "tmp_tc.vcf";

    // Write temp VCF file
    {
        std::ofstream vcfo(tmpvcf);
        vcfo
        << "##fileformat=VCFv4.1" << endl
        << "##phasing=true" << endl
        << "##contig=<ID=x>" << endl
        << "##contig=<ID=y>" << endl
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Freq\">" << endl
        << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternete Allele count\">" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Num samples at site\">" << endl
        << "##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Num alt alleles\">" << endl
        << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"Length of each alt\">" << endl
        << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"type of variant\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2" << endl
        << "x\t9\t.\tG\tA,C,T\t99\t.\tAF=0.01,0.6,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t0|1\t2|3" << endl
        << "x\t10\t.\tC\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.01;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "x\t14\t.\tG\t<DUP>,<BLAH>\t99\t.\tAF=0.01,0.1;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t1|1" << endl
        << "y\t34\t.\tTATA\t<CN2>,<CN0>\t99\t.\tAF=0.01,0.1;AC=2;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|1\t2|1" << endl
        << "y\t39\t.\tT\t<CN0>\t99\t.\tAF=0.01;AC=1;LEN=1;NA=1;NS=1;TYPE=snp\tGT\t1|0\t0|1" << endl;
    }

    SUBCASE("File write wrapper") {

        SUBCASE("Unfiltered") {
            vargas::VCF vcf(tmpvcf);
            vcf.next();
            CHECK(vcf.num_haplotypes() == 4);
            CHECK(vcf.sequences().size() == 2);
            CHECK(vcf.sequences()[0] == "x");
            CHECK(vcf.sequences()[1] == "y");
            REQUIRE(vcf.samples().size() == 2);
            CHECK(vcf.samples()[0] == "s1");
            CHECK(vcf.samples()[1] == "s2");

            REQUIRE(vcf.gen_genotypes().size() == 4);
            CHECK(vcf.gen_genotypes()[0] == "G");
            CHECK(vcf.gen_genotypes()[1] == "A");
            CHECK(vcf.gen_genotypes()[2] == "C");
            CHECK(vcf.gen_genotypes()[3] == "T");
            REQUIRE(vcf.alleles().size() == 4);
            CHECK(vcf.alleles()[0] == "G");
            CHECK(vcf.alleles()[1] == "A");
            CHECK(vcf.alleles()[2] == "C");
            CHECK(vcf.alleles()[3] == "T");
            CHECK(vcf.ref() == "G");
            CHECK(vcf.pos() == 8);

            // Copy number alleles
            vcf.next();
            REQUIRE(vcf.gen_genotypes().size() == 4);
            CHECK(vcf.gen_genotypes()[0] == "CC");
            CHECK(vcf.gen_genotypes()[1] == "CC");
            CHECK(vcf.gen_genotypes()[2] == "");
            CHECK(vcf.gen_genotypes()[3] == "CC");
            REQUIRE(vcf.alleles().size() == 3);
            CHECK(vcf.alleles()[0] == "C");
            CHECK(vcf.alleles()[1] == "CC");
            CHECK(vcf.alleles()[2] == "");
            CHECK(vcf.ref() == "C");
            CHECK(vcf.pos() == 9);

            // Invalid tags
            vcf.next();
            REQUIRE(vcf.alleles().size() == 3);
            CHECK(vcf.alleles()[0] == "G");
            CHECK(vcf.alleles()[1] == "G");
            CHECK(vcf.alleles()[2] == "G");
            CHECK(vcf.ref() == "G");
            CHECK(vcf.pos() == 13);

            // Next y contig should still load
            vcf.next();
            CHECK(vcf.alleles()[0] == "TATA");
        }

        SUBCASE("CHROM Filtering") {
            vargas::VCF vcf;
            vcf.set_region(std::string("y:0-0"));
            vcf.open(tmpvcf);

            vcf.next();
            CHECK(vcf.ref() == "TATA");
            vcf.next();
            CHECK(vcf.ref() == "T");
            CHECK(vcf.next() == 0); // File end
        }

        SUBCASE("Region filtering") {
            vargas::VCF vcf;
            vcf.set_region(std::string("x:0-14"));
            vcf.open(tmpvcf);

            vcf.next();
            CHECK(vcf.ref() == "G");
            vcf.next();
            CHECK(vcf.ref() == "C");
            vcf.next();
            CHECK(vcf.ref() == "G");
            CHECK(vcf.next() == 0); // Region end
        }

        SUBCASE("Ingroup generation") {
            vargas::VCF vcf;
            srand(12345);
            vcf.open(tmpvcf);
            vcf.create_ingroup(50);

            CHECK(vcf.ingroup().size() == 1);
            CHECK(vcf.ingroup()[0] == "s2");

            vcf.next();
            REQUIRE(vcf.gen_genotypes().size() == 2);
            CHECK(vcf.gen_genotypes()[0] == "C");
            CHECK(vcf.gen_genotypes()[1] == "T");

            vcf.next();
            REQUIRE(vcf.gen_genotypes().size() == 2);
            CHECK(vcf.gen_genotypes()[0] == "");
            CHECK(vcf.gen_genotypes()[1] == "CC");

            // Allele set should be complete, ingroup should reflect minimized set
            CHECK(vcf.alleles().size() == 3);
            CHECK(vcf.ingroup().size() == 1);
        }

        SUBCASE("Allele populations") {
            vargas::VCF vcf;
            vcf.open(tmpvcf);
            vcf.next();
            vcf.gen_genotypes();

            REQUIRE(vcf.allele_pop("G").size() == 4);
            CHECK(vcf.allele_pop("G")[0]);
            CHECK(!vcf.allele_pop("G")[1]);
            CHECK(!vcf.allele_pop("G")[2]);
            CHECK(!vcf.allele_pop("G")[3]);

            REQUIRE(vcf.allele_pop("A").size() == 4);
            CHECK(!vcf.allele_pop("A")[0]);
            CHECK(vcf.allele_pop("A")[1]);
            CHECK(!vcf.allele_pop("A")[2]);
            CHECK(!vcf.allele_pop("A")[3]);

            REQUIRE(vcf.allele_pop("C").size() == 4);
            CHECK(!vcf.allele_pop("C")[0]);
            CHECK(!vcf.allele_pop("C")[1]);
            CHECK(vcf.allele_pop("C")[2]);
            CHECK(!vcf.allele_pop("C")[3]);

            REQUIRE(vcf.allele_pop("T").size() == 4);
            CHECK(!vcf.allele_pop("T")[0]);
            CHECK(!vcf.allele_pop("T")[1]);
            CHECK(!vcf.allele_pop("T")[2]);
            CHECK(vcf.allele_pop("T")[3]);

        }

        SUBCASE("Filtered allele populations") {
            vargas::VCF vcf;
            vcf.open(tmpvcf);
            vcf.create_ingroup({"s1"});
            vcf.next();
            vcf.gen_genotypes();

            REQUIRE(vcf.allele_pop("G").size() == 2);
            CHECK(vcf.allele_pop("G")[0]);
            CHECK(!vcf.allele_pop("G")[1]);

            REQUIRE(vcf.allele_pop("A").size() == 2);
            CHECK(!vcf.allele_pop("A")[0]);
            CHECK(vcf.allele_pop("A")[1]);

            REQUIRE(vcf.allele_pop("C").size() == 2);
            CHECK(!vcf.allele_pop("C")[0]);
            CHECK(!vcf.allele_pop("C")[1]);

            REQUIRE(vcf.allele_pop("T").size() == 2);
            CHECK(!vcf.allele_pop("T")[0]);
            CHECK(!vcf.allele_pop("T")[1]);
        }

        SUBCASE("Allele frequencies") {
            vargas::VCF vcf;
            vcf.open(tmpvcf);
            vcf.next();

            auto af = vcf.frequencies();
            REQUIRE(af.size() == 4);
            CHECK(af[0] > 0.289f); // af[0] should be 0.29
            CHECK(af[0] < 0.291f);
            CHECK(af[1] == 0.01f);
            CHECK(af[2] == 0.6f);
            CHECK(af[3] == 0.1f);
        }

    }

    remove(tmpvcf.c_str());
}

TEST_SUITE_END();