# IsoMir: A Pipeline for Detection of isomiRs

IsomiRs are isoforms of the same canonical mature miRNA with alternative length and/or sequence variants.
Growing body of evidence suggested that some isomiRs appear biologically relevant.

This pipeline is designed to organize functions for detecting isomiRs from miRseq data with the Snakemake workflow management system.

The algorithm was borrowed from QuagmiR <https://github.com/Gu-Lab-RBL-NCI/QuagmiR> with some changes:

* We implement the algorithm with **C++**, which make it run faster.
* The workflow was controlled by **Snakemake**.
* Add some helper functions, including **visualization**.
* Beside the algorithm from QuagmiR, we align the short reads to pre-miRNA from  miRBase directly with Smith–Waterman algorithm.

Particularly, this pipeline could be adapted to analyze the kinetics of chimeric miRNA-siRNA, for example the miHTT design in pre-miR-451 backbone (Cells. 2022; 11(17): 2748), and the siHTT design in pre-miR-155 backbone (Brain. 2021; 144(11): 3421–3435). 

## Installation
### Prerequisites
* [R 4.2.2](https://mirrors.tuna.tsinghua.edu.cn/CRAN/)
  + tidyverse
  + yaml
* [Bioconductor 3.16](https://www.bioconductor.org/install/)
  + Biostrings
* samtools
```sh
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar xjvf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure
make
make install
```

* Python 3.10
  + pandas
  + snakemake

### Download and compile
```sh
git clone https://github.com/exornabio/isomir.git
# Move to c++ source code directory and compile
cd isomir/src
make
# The compiled binary is called `isomir`
```

## Workflow

### Mark duplicates in the FASTQ file
To speed up the computation, we first extract sequences of each reads, and then mark and count duplicates. 
The sequence with its amount are stored in a tab-separated values (TSV) file.

### Split reads file
Split the reads file into multiple files, which can be scattered to parallelize the task of identifing isomiRs.

### Identify reads as certein miRNA based on the motif
Reads matching a certain motif were considered as potential isomiRs for the corresponding miRNA.

### Caculate edit distance of 5′ or 3′ ends
For the potential isomiR reads, the 5' part and 3' part (preced and follow the motif) sequences are further compared with reference miRNA respetively.
The pairwise sequence similarity between 5' part and 3' part with the reference miRNA is calculated using the Levenshtein distances (edit distances).
The maximum allowed edit distances for the 5' and 3' regions can be set 
independently to capture the asymmetrical sequence heterogeneity.

### Visualize 
The output of IsomiRs are stored as SAM/BAM format, which can be interactively visuallized by genomics viewer tools, such as [The Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/)
![](img/igv.png)

## Usage
Download the miRNA file (`miRNA.xls`) from miRBase <https://mirbase.org/ftp/> and transform it to TSV file (e.g. `mibase22.1.tsv`).
Set paramaters in `config/config.yaml`:

```yaml
samples: config/sample.txt
spe: mmu
mibase: mibase22.1.tsv
kmer_len: 13
max_edit_dist_5p: 2
max_edit_dist_3p: 3
min_identity: 0.9
```

* `max_edit_dist_5p` is the maximum distance between 5' part of reads with reference miRNA.
* `max_edit_dist_3p` is the maximum distance between 3' part of reads with reference miRNA.

After that, run
```sh
snakemake --cores 4
```

## Result files
The result files were stored in `data/result/`

### sample_isoform.sam
Sequence Alignment Map (SAM) format file of reads located on pre-miRNA.

```
read4868	0	hsa-mir-151a	11	30	22M	*	0	0	TCGAGGAGCTCACAGTCTAGTA	*	RN:i:75	DT:i:1
read4867	0	hsa-mir-151a	11	30	21M	*	0	0	TCGAGGAGCTCACAGTCTAGT	*	RN:i:64	DT:i:0
read4869	0	hsa-mir-151a	11	30	23M	*	0	0	TCGAGGAGCTCACAGTCTAGTAA	*	RN:i:14	DT:i:2
```

In the SAM format, each alignment line represents the alignment of a short reads. 
Each line consists of 11 or more TAB-separated columns. 
For more details see <https://samtools.github.io/hts-specs/SAMv1.pdf>. 
Here we just list some import columns in this case.

	* Col. 1: read ids
	* Col. 3: pre-miRNA
	* Col. 4: 1-based leftmost mapping positon of the reads
	* Col. 6: CIGAR string, reprenting match, insertiong, deletion, et al. in an alignment.
	* Col. 10: sequence of the read
	* Optonal columns: `RN` means amonut of the read in raw sequencing FASTQ file; `ID` repesents identity/similarity of the local alignment between read and reference.

### sample_isoform.bam
The compressed binary veriosn of SAM (BAM) of `sample_isoform.sam`.

### sample_hit.sam
SAM format file of reads located on pre-miRNA.
The mapping was implemented with Smith–Waterman algorithm.

### sample_hit_bam
The compressed binary veriosn of SAM (BAM) file of `sample_hit.sam`.

### mirna.sam and mirna.bam
SAM/BAM files of canonical miRNA.
