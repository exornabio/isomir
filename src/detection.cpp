#include "detection.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include "levenshtein.h"
#include "code.h"

using boost::bad_lexical_cast;
using boost::lexical_cast;

void Detection::detect()
{
    load_mirna_();

    for (auto it = motif2mirna_.cbegin(); it != motif2mirna_.cend(); it++)
    {
        string motif = it->first;
        Mirna *mirna = it->second;
        motifs_.insert(it->first);
        if (consensus_contain_amb_(mirna->consensus_))
        {
            motif2mirna_amb_.insert({motif, mirna});
        }
        else
        {
            motif2mirna_no_amb_.insert({motif, mirna});
        }
    }

    ifstream fin(config_->reads_file_);
    if (!fin.is_open())
    {
        cerr << "Could not open " << config_->reads_file_ << endl;
        exit(1);
    }

    string line;
    while (getline(fin, line))
    {
        boost::trim(line);
        vector<string> ss;
        boost::split(ss, line, boost::is_any_of("\t"));
        if (ss.size() != 3)
            continue;

        string read_id = ss[0];
        int reads_num;
        
        try
        {
            reads_num = lexical_cast<int>(ss[1]);
        }
        catch (bad_lexical_cast &)
        {
        }
        string seq = ss[2];
        int seq_len = seq.size();
        detect_one_seq_(read_id, seq, reads_num, motif2mirna_amb_, true);
        detect_one_seq_(read_id, seq, reads_num, motif2mirna_no_amb_, false);
    }

    output_isoform_();

    cout << "Done" << endl;
}

void Detection::load_mirna_()
{
    ifstream fin;
    fin.open(config_->mirna_file_.c_str());
    if (!fin.is_open())
    {
        cerr << "Could not open " << config_->mirna_file_ << endl;
        exit(1);
    }

    string line;
    // skip header
    getline(fin, line);
    while (getline(fin, line))
    {
        boost::trim(line);
        vector<string> ss;
        boost::split(ss, line, boost::is_any_of("\t"));
        // hsa-let-7a-5p	TGAGGTAGTAGGTTGTATAGTT	GTAGTAGGTTGTA	2
        if (ss.size() != 4)
            continue;
        string mirna_id = ss[0];
        string consensus = ss[1];
        string motif = ss[2];
        motif2mirna_.insert({motif, new Mirna(mirna_id, motif, consensus)});
    }
}

bool Detection::consensus_contain_amb_(const string &consensus)
{
    for (auto it = AMBIGUOUS_LETTERS.cbegin(); it != AMBIGUOUS_LETTERS.cend(); it++)
    {
        if (consensus.find(*it) != string::npos)
            return true;
    }
    return false;
}

bool Detection::compare_mirna_and_seq_nt_(char mirna_nt, char seq_nt)
{
    set<char> codes = COMPARISON.at(mirna_nt);
    return codes.find(seq_nt) != codes.cend();
}

bool Detection::compare_mirna_and_seq_(const string &mirna, const string &seq)
{
    int seq_len = seq.size();
    if (seq_len != mirna.size())
        return false;

    for (int i = 0; i < seq_len; i++)
    {
        if (!compare_mirna_and_seq_nt_(mirna[i], seq[i]))
        {
            return false;
        }
    }
    return true;
}

int Detection::calc_edit_dist_(const string &seq_5p, const string &seq_3p,
                               const string &consensus_5p, const string &consensus_3p)
{

    int lev_5p = edit_dist(seq_5p, consensus_5p);
    int lev_3p = edit_dist(seq_3p, consensus_3p);
    int dist = lev_3p + lev_5p;
    if (config_->max_edit_dist_5p_ != -1 && config_->max_edit_dist_3p_ == -1)
    {
        if (lev_5p > config_->max_edit_dist_5p_)
            dist = -5;
    }
    else if (config_->max_edit_dist_5p_ == -1 && config_->max_edit_dist_3p_ != -1)
    {
        if (lev_3p > config_->max_edit_dist_3p_)
            dist = -3;
    }
    else if (config_->max_edit_dist_5p_ != -1 && config_->max_edit_dist_3p_ != -1)
    {
        if (lev_5p > config_->max_edit_dist_5p_)
            dist = -5;

        if (lev_3p > config_->max_edit_dist_3p_)
            dist = -3;
    }
    return dist;
}

int Detection::find_motif_in_seq_amb_(const string &motif, const string &seq)
{
    int motif_len = motif.size();
    for (int i = 0; i < seq.size() - motif_len + 1; i++)
    {
        if (compare_mirna_and_seq_(motif, seq.substr(i, motif_len)))
        {
            return i;
        }
    }
    return -1;
}

int Detection::find_motif_in_seq_(const string &motif, const string &seq)
{
    int pos = seq.find(motif);
    if (pos == string::npos)
    {
        pos = -1;
    }
    return pos;
}

bool Detection::seq_contains_motif_(const string &motif, const string &seq)
{
    if (find_motif_in_seq_(motif, seq) == -1)
    {
        return false;
    }
    return true;
}

bool Detection::seq_contains_motif_amb_(const string &motif, const string &seq)
{
    if (find_motif_in_seq_amb_(motif, seq) == -1)
    {
        return false;
    }
    return true;
}

bool Detection::seq_contain_n_(const string &seq)
{
    if (seq.find('N') == string::npos)
        return false;

    return true;
}

void Detection::detect_one_seq_(const string& read_id, const string &seq, int read_num,
                                const map<string, Mirna *> motif2mirna, bool is_amb)
{

    vector<Isoform *> hits;
    if (seq_contain_n_(seq))
        is_amb = true;

    for (auto it = motif2mirna.cbegin(); it != motif2mirna.cend(); it++)
    {
        string motif = it->first;
        string mirna_id = it->second->id_;
        string consensus = it->second->consensus_;

        if (is_amb)
        {
            if (!seq_contains_motif_amb_(motif, seq))
                continue;
        }
        else
        {
            if (!seq_contains_motif_(motif, seq))
                continue;
        }

        int consensus_index_5p = consensus.find(motif);
        if (consensus_index_5p == string::npos)
        {
            cerr << "Motif: " << motif << " is not found in consensus" << endl;
            exit(1);
        }
        int consensus_index_3p = consensus_index_5p + motif.size();

        int seq_index_5p = -1;
        if (is_amb)
        {
            seq_index_5p = find_motif_in_seq_amb_(motif, seq);
        }
        else
        {
            seq_index_5p = find_motif_in_seq_(motif, seq);
        }

        if (seq_index_5p == -1) 
            continue;

        int seq_index_3p = seq_index_5p + motif.size();
        
        string consensus_5p = consensus.substr(0, consensus_index_5p);
        string consensus_3p = consensus.substr(consensus_index_3p);
        string seq_5p = seq.substr(0, seq_index_5p);
        string seq_3p = seq.substr(seq_index_3p);

        int dist = calc_edit_dist_(seq_5p, seq_3p, consensus_5p, consensus_3p);

        if (dist == -3)
        {
            continue;
        }
        if (dist == -5)
        {
            continue;
        }

        Isoform *isoform = new Isoform(mirna_id, read_id, seq, read_num, dist);
        hits.push_back(isoform);
    }

    if (!hits.empty())
    {
        int min_dist = hits[0]->dist_;

        for (int i = 1; i < hits.size(); i++)
        {
            Isoform *hit = hits[i];
            if (hit->dist_ < min_dist)
            {
                min_dist = hit->dist_;
            }
        }

        for (int i = hits.size() - 1; i >= 0; i--)
        {
            if (hits[i]->dist_ == min_dist)
            {
                isoforms_.push_back(hits[i]);
            }
            else
            {
                delete hits[i];
            }
        }
    }
}

void Detection::output_isoform_()
{
    ofstream fout(config_->isoform_file_.c_str());
    if (!fout.is_open())
    {
        cerr << "Cannot open " << config_->isoform_file_ << "file for output" << endl;
    }
    fout << "mirna_id\tread_id\tread_seq\tread_num\tdist" << endl;
    for (auto it = isoforms_.cbegin(); it != isoforms_.cend(); it++)
    {
        Isoform *isoform = *it;
        fout.precision(2);
        fout << isoform->mirna_id_ << "\t"  << isoform->read_id_ << "\t" << isoform->read_seq_ << "\t" << isoform->read_num_ << "\t" << isoform->dist_ << endl;
    }
    fout.close();
}

void Detection::clear()
{
    for (auto it = isoforms_.begin(); it != isoforms_.end(); it++)
    {
        delete (*it);
    }

    for (auto it = motif2mirna_.begin(); it != motif2mirna_.end(); it++)
    {
        delete (it->second);
    }
}
