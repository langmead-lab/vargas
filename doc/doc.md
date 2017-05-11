# Vargas
_Updated: May 11, 2017_

[![Build Status](https://travis-ci.org/RaviGaddipati/vargas.svg?branch=master)](https://travis-ci.org/RaviGaddipati/vargas)

Vargas aligns short reads to a directed acyclic graph (DAG). Reads are aligned using a SIMD vectorized version of Smith-Waterman, aligning multiple reads at once. The the aligner reports the best positions of reads without traceback to keep memory usage manageable. Additionally, vargas can be used to simulate reads from a DAG, and reads can be stratified by various parameters.

Source documentation via Doxygen. Usage documantation is on [Github](https://github.com/RaviGaddipati/vargas).