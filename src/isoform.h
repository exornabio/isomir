#ifndef ISOFORM_H
#define ISOFORM_H

#include "base.h"

struct Isoform
{
    string mirna_id_;
    string read_id_;
    string read_seq_;
    int read_num_;
    int dist_;

    Isoform(const string &mirna_id, const string &read_id,
            const string &read_seq, int read_num, int dist)
    {
        mirna_id_ = mirna_id;
        read_id_ = read_id;
        read_seq_ = read_seq;
        read_num_ = read_num;
        dist_ = dist;
    }
};

#endif // IOSFORM_H