/**
 * Ravi Gaddipati
 * November 23, 2015
 * rgaddip1@jhu.edu
 *
 * @brief
 * Provides tools to interact with VCF/BCF files.
 *
 * @file
 */

#include "varfile.h"

vargas::Region vargas::parse_region(const std::string &region_str) {
    vargas::Region ret;

    std::vector<std::string> regionSplit = split(region_str, ':');

    // Name
    if (regionSplit.size() != 2)
        throw std::invalid_argument("Invalid region format, should be CHR:XX,XXX-YY,YYY");

    ret.seq_name = regionSplit[0];

    // Strip commas
    regionSplit[1].erase(std::remove(regionSplit[1].begin(), regionSplit[1].end(), ','), regionSplit[1].end());

    // Range
    regionSplit = split(regionSplit[1], '-');
    if (regionSplit.size() != 2)
        throw std::invalid_argument("Invalid region format, should be CHR:XX,XXX-YY,YYY");

    ret.min = std::stoul(regionSplit[0]);
    ret.max = std::stoul(regionSplit[1]);

    if (ret.min > ret.max) {
        throw std::invalid_argument("Invalid region, min > max.");
    }
    return ret;
}

void vargas::VariantFile::set_region(std::string region) {
    auto parsed = parse_region(region);
    _chr = parsed.seq_name;
    _min_pos = parsed.min;
    _max_pos = parsed.max;
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
    if (!_header || !_bcf) return false;
    if (bcf_read(_bcf, _header, _curr_rec) != 0) return false;
    unpack_all();
    // Check if its within the filter range
    if (_max_pos > 0 && _curr_rec->pos > _max_pos) {
        return false;
    }
    if (_curr_rec->pos < _min_pos ||
        (_chr.length() != 0 && strcmp(_chr.c_str(), bcf_hdr_id2name(_header, _curr_rec->rid)))) {
        return next();
    }
    genotypes();
    return true;
}


const std::vector<std::string> &vargas::VCF::genotypes() {
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
    _bcf = bcf_open(_file_name.c_str(), "r");
    if (!_bcf) return -1;

    _header = bcf_hdr_read(_bcf);
    if (_header == nullptr) {
        return -2;
    }

    // Load samples
    for (size_t i = 0; i < num_samples() / 2; ++i) {
        _samples.push_back(_header->samples[i]);
    }
    create_ingroup(100);
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