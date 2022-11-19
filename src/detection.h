#ifndef ETECTION_H
#define DETECTION_H

#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include "base.h"
#include "config.h"
#include "code.h"
#include "mirna.h"
#include "isoform.h"

struct Detection
{
    Config *config_;
    set<string> motifs_;
    map<string, Mirna *> motif2mirna_;
    map<string, Mirna *> motif2mirna_amb_;
    map<string, Mirna *> motif2mirna_no_amb_;
    map<string, vector<Isoform *>> id2isoforms_;
    vector<Isoform *> isoforms_;

    Detection(Config *config) : config_(config) {}

    void detect();

    void clear();

private:
    /**
     * @brief fill motif2mirna_ from mirna_file_
     */
    void load_mirna_();

    /**
     * @brief check if consensus contain ambiguous nt
     * @retval true contained
     * @retval false NOT contained
     */
    bool consensus_contain_amb_(const string &consensus);

    /**
     * @brief compare the mirna nt and read nt using COMPARISON
     * @param mirna_nt mirna sequence nt
     * @param seq_nt read sequence nt
     * @retval true same
     * @retval false diffrent
     */
    bool compare_mirna_and_seq_nt_(char motif_nt, char seq_nt);

    /**
     * @brief compare the mirna and read sequence with same length
     * @param mirna mirna sequence
     * @param seq read sequence
     * @retval true same
     * @retval false diffrent
     */
    bool compare_mirna_and_seq_(const string &mirna, const string &seq);

    /**
     * @brief calculate edit distance
     * @param seq_5p 5' part of read sequence
     * @param seq_3p 3' part of read sequence
     * @param consensus_5p 5' part of consensus
     * @param consensus_3p 3' part of consensus
     *
     * @return the distance
     *        -5 higher than cutoff on 5' end
     *        -3 lower than cutoff on 3' end
     */
    int calc_edit_dist_(const string &seq_5p, const string &seq_3p,
                        const string &consensus_5p, const string &consensus_3p);

    /**
     * @brief check if read sequence contain N
     * @retval true
     * @retval false
     */
    bool seq_contain_n_(const string &seq);

    /**
     * @brief find isoform on one seq, put into isoforms_
     * @param read_id read ID
     * @param seq read sequence
     * @param read_num duplicated number of the seq
     * @param motif2mirna
     * @param is_amb
     */
    void detect_one_seq_(const string& read_id, const string &seq, int read_num,
                         const map<string, Mirna *> motif2mirna, bool is_amb);

    /**
     * @brief find the positon of motif on ambiguous seq
     * @param motif
     * @param seq
     * @return position of motif on seq
     *
     *    -1 not found
     */
    int find_motif_in_seq_amb_(const string &motif, const string &seq);

    /**
     * @brief find the positon of motif on seq
     * @param motif
     * @param seq
     * @return position of motif on seq 
     *      --1 not found
     */
    int find_motif_in_seq_(const string &motif, const string &seq);

    /**
     * @brief check if read sequence contain the moitf
     * @param seq
     * @param motif
     * @retval true contained
     * @retval false not contained
     */
    bool seq_contains_motif_(const string &motif, const string &seq);

    /**
     * @brief check if reads contain the moitf with ambiguous nt
     * @param motif
     * @param seq
     *
     * @return
     *  -true: contained
     *  -false: not contained
     */
    bool seq_contains_motif_amb_(const string &motif, const string &seq);

    void output_isoform_();
};

#endif // DETECTION_H
