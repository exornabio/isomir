#ifndef MIRNA_H
#define MIRNA_H

#include "base.h"

struct Mirna
{
    string id_;
    string motif_;
    string consensus_;

    Mirna(const string &id, const string& motif, const string& consensus)
        : id_(id), motif_(motif), consensus_(consensus) {}
};

#endif // MIRNA_H