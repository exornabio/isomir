#!/usr/bin/python3

import argparse
import os
from typing import List

from util import Isoform, Mirna, edit_dist, find_motif_in_seq, seq_contain_n

parser = argparse.ArgumentParser()
parser.add_argument("reads_file", help="reads file")
parser.add_argument("mirna_file", help="miRNA file")
parser.add_argument("isoform_file", help="output isoform file")
parser.add_argument("-l", "--max_edit_dist_5p", help="max edit distance of 5' end",
                    type=int, default=2)
parser.add_argument("-r", "--max_edit_dist_3p", help="max edit distance of 3' end",
                    type=int, default=3)
parser.add_argument("-n", "--min_read_num",
                    help="minimum read numbers", type=int, default=100)
args = parser.parse_args()

reads_file = args.reads_file
mirna_file = args.mirna_file
isoform_file = args.isoform_file
max_edit_dist_5p = args.max_edit_dist_5p
max_edit_dist_3p = args.max_edit_dist_3p
min_read_num = args.min_read_num

if not os.path.exists(reads_file):
    print("Reads file doesn't exist")
    exit(1)

if not os.path.exists(mirna_file):
    print("miRNA file doesn't exist")
    exit(1)

mirnas: List[Mirna] = []
isoforms: List[Isoform] = []

def load_mirna():
    """ Fill mirnas from mirna_file
    """
    with open(mirna_file, "rt") as fin:
        fin.readline()  # skip header
        for line in fin:
            # mmu   hsa-let-7a-5p	TGAGGTAGTAGGTTGTATAGTT	GTAGTAGGTTGTA
            ss = line.strip().split("\t")
            if len(ss) == 3:
                mirna_id = ss[0]
                consensus = ss[1]
                motif = ss[2]
                mirnas.append(Mirna(mirna_id, motif, consensus))


def calc_edit_dist(seq_5p: str, seq_3p: str, consensus_5p: str, consensus_3p: str) -> int:
    lev_5p = edit_dist(seq_5p, consensus_5p)
    lev_3p = edit_dist(seq_3p, consensus_3p)
    dist = lev_3p + lev_5p
    if max_edit_dist_5p != -1 and max_edit_dist_3p == -1:
        if lev_5p > max_edit_dist_5p:
            dist = -5
    elif max_edit_dist_5p == -1 and max_edit_dist_3p != -1:
        if lev_3p > max_edit_dist_3p:
            dist = -3
    elif max_edit_dist_5p != -1 and max_edit_dist_3p != -1:
        if lev_5p > max_edit_dist_5p:
            dist = -5
        if lev_3p > max_edit_dist_3p:
            dist = -3
    return dist


def detect_one_seq(read_id: str, seq: str, read_num: int, mirnas: List[Mirna]):
    """Find isoform on one seq, put into isoforms
    
    Parameters
    ----------
    read_id: str, read ID

    seq read: str, sequence

    read_num: int, duplicated number of the seq

    mirnas: List[Mirna]
    """

    is_amb = False
    hits: List[Isoform] = []
    if seq_contain_n(seq):
        is_amb = True

    for mirna in mirnas:
        motif = mirna.motif
        mirna_id = mirna.id
        consensus = mirna.consensus
        seq_index_5p = find_motif_in_seq(motif, seq, is_amb)
        if seq_index_5p == -1:
            continue

        consensus_index_5p = consensus.find(motif)
        if consensus_index_5p == -1:
            print("Motif: " + motif + " is not found in consensus")
            exit(1)

        consensus_index_3p = consensus_index_5p + len(motif)
        seq_index_3p = seq_index_5p + len(motif)
        consensus_5p = consensus[0:consensus_index_5p]
        consensus_3p = consensus[consensus_index_3p:]
        seq_5p = seq[0:seq_index_5p]
        seq_3p = seq[seq_index_3p:]
        dist = calc_edit_dist(seq_5p, seq_3p, consensus_5p, consensus_3p)
        if dist == -3:
            continue
        if dist == -5:
            continue
        hits.append(Isoform(mirna_id, read_id, seq, read_num, dist))

    if len(hits) != 0:
        min_dist = hits[0].dist

        for i in range(1, len(hits)):
            hit = hits[i]
            if hit.dist < min_dist:
                min_dist = hit.dist

        for i in range(len(hits)):
            if hits[i].dist == min_dist:
                isoforms.append(hits[i])


if __name__ == '__main__':
    load_mirna()

    with open(reads_file, "rt") as fin:
        for line in fin:
            ss = line.strip().split("\t")
            if len(ss) != 3:
                continue

            read_id = ss[0]
            read_num = int(ss[1])
            seq = ss[2]
            seq_len = len(seq)
            if read_num >= min_read_num:
                detect_one_seq(read_id, seq, read_num, mirnas)

    with open(isoform_file, "wt") as fout:
        fout.write("mirna_id\tread_id\tread_seq\tread_num\tdist\n")
        for isoform in isoforms:
            fout.write(
                f"{isoform.mirna_id}\t{isoform.read_id}\t{isoform.read_seq}\t{isoform.read_num}\t{isoform.dist}\n")

    print("Done")
