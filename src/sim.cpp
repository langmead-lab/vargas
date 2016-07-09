/**
 * @author Ravi Gaddipati
 * @date June 26, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Simulates random reads from a graph, returning reads that follow a specified Sim::Profile.
 *
 * @file
 */

#include "sim.h"

bool Vargas::Sim::update_read() {

    if (_prof.var_nodes == 0 && _prof.var_bases > 0)
        throw std::invalid_argument("Invalid profile option var_nodes = 0, var_bases > 0.");

    auto &nodes = *_graph.node_map();
    auto &next = _graph.next_map();

    unsigned int curr_node;
    int curr_indiv = -1;
    std::string read_str = "";

    curr_node = next_keys[rand() % next_keys.size()]; // Random start node ID
    size_t curr_pos = rand() % nodes[curr_node]->length(); // current pos relative to node origin

    int var_bases = 0;
    int var_nodes = 0;

    while (true) {
        // The first time a branch is hit, pick an individual to track
        if (curr_indiv < 0 && !nodes[curr_node]->is_ref()) {
            do {
                curr_indiv = rand() % _graph.pop_size();
            } while (nodes[curr_node]->individuals()[curr_indiv] == 0);
        }

        // Extract len subseq
        size_t len = _prof.len - read_str.length();
        if (len > nodes.at(curr_node)->length() - curr_pos) len = nodes.at(curr_node)->length() - curr_pos;
        read_str += nodes.at(curr_node)->seq_str().substr(curr_pos, len);
        curr_pos += len;

        if (!nodes[curr_node]->is_ref()) {
            ++var_nodes;
            var_bases += len;
        }

        if (read_str.length() >= _prof.len) break; // Done

        // Pick random next node.
        if (next.find(curr_node) == next.end()) return update_read(); // End of graph
        std::vector<uint32_t> valid_next;
        for (auto n : next.at(curr_node)) {
            if (curr_indiv < 0 || nodes[n]->belongs(curr_indiv)) valid_next.push_back(n);
        }
        if (valid_next.size() == 0) return update_read();
        curr_node = valid_next[rand() % valid_next.size()];
        curr_pos = 0;
    }

    if (std::count(read_str.begin(), read_str.end(), 'N') >= _prof.len / 2) return update_read();
    if (_prof.var_nodes >= 0 && var_nodes != _prof.var_nodes) return update_read();
    if (_prof.var_bases >= 0 && var_bases != _prof.var_bases) return update_read();

    // Introduce errors
    int sub_err = 0;
    int indel_err = 0;
    std::string read_mut = "";


    // Rate based errors
    if (_prof.rand) {
        for (size_t i = 0; i < read_str.length(); ++i) {
            char m = read_str[i];
            // Mutation error
            if (rand() % 10000 < 10000 * _prof.mut) {
                do {
                    m = rand_base();
                } while (m == read_str[i]);
                ++sub_err;
            }

            // Insertion
            if (rand() % 10000 < 5000 * _prof.indel) {
                read_mut += rand_base();
                ++indel_err;
            }

            // Deletion (if we don't enter)
            if (rand() % 10000 > 5000 * _prof.indel) {
                read_mut += m;
                ++indel_err;
            }
        }
    }
    else {
        // Fixed number of errors
        sub_err = (int) std::round(_prof.mut);
        indel_err = (int) std::round(_prof.indel);
        std::vector<int> indel_pos;
        for (int i = 0; i < indel_err;) {
            int r = rand() % read_str.length();
            if (std::find(indel_pos.begin(), indel_pos.end(), r) == indel_pos.end()) {
                indel_pos.push_back(r);
                ++i;
            }
        }
        std::sort(indel_pos.begin(), indel_pos.end());
        int prev = 0;
        for (int p : indel_pos) {
            read_mut += read_str.substr(prev, p);
            if (rand() % 2) {
                // Insertion
                read_mut += rand_base();
                read_mut += read_str[p];
            }
            prev = p + 1;
        }
        read_mut += read_str.substr(prev, std::string::npos);

        for (int i = 0; i < sub_err;) {
            int r = rand() % read_str.length();
            if (read_str[r] == read_mut[r]) { // Make sure we don't double mutate same base
                do {
                    read_mut[r] = rand_base();
                } while (read_str[r] == read_mut[r]);
                ++i;
            }
        }
    }

    _read = SAM::Record();

    _read.flag.unmapped = false;
    _read.flag.aligned = true;

    _read.seq = read_mut;
    _read.aux.set(SIM_SAM_INDIV_TAG, curr_indiv);
    _read.aux.set(SIM_SAM_INDEL_ERR_TAG, indel_err);
    _read.aux.set(SIM_SAM_VAR_BASE_TAG, var_bases);
    _read.aux.set(SIM_SAM_VAR_NODES_TAG, var_nodes);
    _read.aux.set(SIM_SAM_SUB_ERR_TAG, sub_err);

    _read.pos = nodes[curr_node]->end() - nodes[curr_node]->length() + curr_pos - _prof.len;

    //TODO set cigar

    _read.aux.set(SIM_SAM_READ_ORIG_TAG, read_str);

    return true;
}