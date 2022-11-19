#ifndef CONFIG_H
#define CONFIG_H

#include "base.h"

struct Config {
    int max_edit_dist_5p_ = 2;
    int max_edit_dist_3p_ = 3;
    string reads_file_;
    string mirna_file_;
    string isoform_file_;

    void parse_arg(int argc, char *argv[]);

private:
    void print_usag_();
    void check_arg_();
    void check_file_();
};

#endif //CONFIG_H
