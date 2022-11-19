#ifndef CODE_H
#define CODE_H

#include "base.h"

const map<char, set<char>> COMPARISON = {
    {'A', {'A', 'N'}},
    {'C', {'C', 'N'}},
    {'G', {'G', 'N'}},
    {'T', {'T', 'N'}},
    {'R', {'A', 'G', 'N'}},
    {'Y', {'C', 'T', 'N'}},
    {'S', {'C', 'G', 'N'}},
    {'W', {'A', 'T', 'N'}},
    {'K', {'G', 'T', 'N'}},
    {'M', {'A', 'C', 'N'}},
    {'B', {'C', 'G', 'T', 'N'}},
    {'D', {'A', 'G', 'T', 'N'}},
    {'H', {'A', 'C', 'T', 'N'}},
    {'V', {'A', 'C', 'G', 'N'}},
    {'N', {'A', 'C', 'G', 'T', 'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'}}};

const set<char> AMBIGUOUS_LETTERS = {'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'};

#endif // CODE_H