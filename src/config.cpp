#include <cstdlib>
#include "config.h"
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include "base.h"

using boost::bad_lexical_cast;
using boost::lexical_cast;

void Config::parse_arg(int argc, char *argv[])
{
    if (argc < 2)
        print_usag_();

    int c;
    while ((c = getopt(argc, argv, "l:r:s:m:o:h")) != -1)
    {
        switch (c)
        {
        case 'l':
            try
            {
                max_edit_dist_5p_ = lexical_cast<int>(optarg);
            }
            catch (bad_lexical_cast &)
            {
                print_usag_();
            }
            break;
        case 'r':
            try
            {
                max_edit_dist_3p_ = lexical_cast<int>(optarg);
            }
            catch (bad_lexical_cast &)
            {
                print_usag_();
            }
            break;
        case 's':
            reads_file_ = optarg;
            break;
        case 'm':
            mirna_file_ = optarg;
            break;
        case 'o':
            isoform_file_ = optarg;
            break;
        case 'h':
            print_usag_();
            break;
        default:
            print_usag_();
        }
    }
    check_arg_();
    check_file_();
}

void Config::print_usag_()
{
    const string USAGE = "isomir [option] args\n"
        "-l max edit distance of 5' end (optional, default = 2)\n"
        "-r max edit distance of 3' end (optional, default = 3)\n"
        "-s reads file (required)\n"
        "-m miRNA file (required)\n"
        "-o output isoform file (required)\n";

    cout << USAGE << endl;
    exit(1);
}

void Config::check_arg_()
{
    if (reads_file_.empty() || mirna_file_.empty() || isoform_file_.empty())
        print_usag_();
}

void Config::check_file_()
{
    boost::filesystem::path reads_path(reads_file_);
    if (!boost::filesystem::exists(reads_path))
    {
        cerr << "Reads file doesn't exist" << endl;
        exit(1);
    }

    boost::filesystem::path mirna_path(mirna_file_);
    if (!boost::filesystem::exists(mirna_path))
    {
        cerr << "miRNA file doesn't exist" << endl;
        exit(1);
    }
}
