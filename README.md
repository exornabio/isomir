# IsoMir: A Pipeline for Detection of isomiRs

IsomiRs are isoforms of the same canonical mature miRNA with alternative length and sequence variants.
A growing body of evidence suggested that some isomiRs appear biologically relevant.

This pipeline is designed to organize functions for detecting isomiRs from miRseq data with the Snakemake workflow management system.

The algorithm was borrowed from QuagmiR <https://github.com/Gu-Lab-RBL-NCI/QuagmiR> with some changes:

* The workflow was controlled by **Snakemake**.
* Add some helper functions, including **visualization**.
* This pipeline can be used to analyze the kinetics of chimeric miRNA-siRNA.

## Installation
### Prerequisites
* Ubuntu 22.04 LTS
* [R 4.3.0](https://mirrors.tuna.tsinghua.edu.cn/CRAN/)
    + tidyverse
* [Bioconductor 3.18](https://www.bioconductor.org/install/)
    + Biostrings
* Python 3.10+
    + snakemake
    + biopython
    + pandas
* [samtools](https://github.com/samtools/samtools)
* fastp

### Download and compile

## Workflow

### Prepare miRNA dataset
Download the miRNA file (`miRNA.dat`) from miRBase <https://mirbase.org/download/>.  
Extract mature miRNA with the attribute, including precursor sequences and the coordinates of the reference miRNAs on the precursors for the species the user selected.

### Trim adapter and control quality 
Trim adapter and control quality of NGS reads using fastp.

### Mark duplicates in the FASTQ file
To speed up the computation, we first extract sequences of each read and then mark and count duplicates. 
The sequence and amount are stored in a tab-separated values (TSV) file.

### Identity reads as candidate miRNA based on the motif
Reads matching a specific motif were considered as potential isomiRs for the corresponding miRNA.

### Caculate edit distance of 5′ or 3′ ends
For the potential isomiR reads, the 5' part and 3' part (precede and follow the motif) sequences are further compared with reference miRNA, respectively.
The pairwise sequence similarity between the 5' part and 3' part with the reference miRNA is calculated using the Levenshtein distances (edit distances).
The maximum allowed edit distances for the 5' and 3' regions can be set 
independently to capture the asymmetrical sequence heterogeneity.

### Visualize 
The output of IsomiRs is stored in SAM/BAM format, which can be interactively visualized by genomics viewer tools, such as [The Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/)
![](img/igv.png)


## Usage
Set parameters in `data/config/config.yaml`:

```yaml
samples: sample.ids
spe: mmu
kmer_len: 13
max_edit_dist_5p: 2
max_edit_dist_3p: 3
min_read_num: 100
fastp:
  thread: 4
```

* `samples`: the file that contains sample names per line
* `spe`: the species (e.g. hsa, mmu)
* `kmer_len`: length of motif in mature miRNA
* `max_edit_dist_5p`: the maximum distance between 5 part of reads with reference miRNA
* `max_edit_dist_3p`: the maximum distance between 3' part of reads with reference miRNA
* `min_read_num`: is the minimum read count for an isoform

After that, run
```sh
snakemake --cores 4 --forceall
```

## Result files
The result files were stored in `data/result/`

### sample_isoform.sam
Sequence Alignment Map (SAM) format file of reads located on pre-miRNA.

```
read4868    0   hsa-mir-151a    11  30  22M *   0   0   TCGAGGAGCTCACAGTCTAGTA  *   RN:i:75 DT:i:1
read4867    0   hsa-mir-151a    11  30  21M *   0   0   TCGAGGAGCTCACAGTCTAGT   *   RN:i:64 DT:i:0
read4869    0   hsa-mir-151a    11  30  23M *   0   0   TCGAGGAGCTCACAGTCTAGTAA *   RN:i:14 DT:i:2
```

In the SAM format, each alignment line represents the alignment of short reads. 
Each line consists of 11 or more TAB-separated columns; for more details, see <https://samtools.github.io/hts-specs/SAMv1.pdf>. 
Here, we list some important columns in this case.

* Col. 1: read id
* Col. 3: pre-miRNA
* Col. 4: 1-based leftmost mapping position of the reads
* Col. 6: CIGAR string, representing match, insertion, deletion, et al. in an alignment.
* Col. 10: sequence of the read
* Optional columns: `RN` means the amount of the read in raw sequencing FASTQ file; `ID` represents identity/similarity of the local alignment between read and reference.

### sample_isoform.bam
The compressed binary version of SAM (BAM) of `sample_isoform.sam`.

### mirna.sam and mirna.bam
SAM/BAM files of canonical miRNA.

## Docker

```sh
docker build -t caibinperl/isomir .
docker push caibinperl/isomir
```

```sh
# Mac / Linux
docker run --mount type=bind,src="$(pwd)/data",target=/isomir/data caibinperl/isomir
# PowerShell
docker run --mount "type=bind,src=$pwd/data,target=/isomir/data" caibinperl/isomir
```
