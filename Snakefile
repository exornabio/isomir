import os
import pandas as pd

configfile: "config/config.yaml"

os.system("Rscript script/process_mibase.R")


def start_on():
    data_dir = "data"
    resource_dir = os.path.join(data_dir, "resource")
    read_dir = os.path.join(data_dir, "read")
    result_dir = os.path.join(data_dir, "result")
    if not os.path.exists(resource_dir):
        os.mkdir(resource_dir)
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
  

def get_samples():
    samples = set()
    with open(config["samples"], "rt") as fin:
        for line in fin:
            line = line.strip()
            if (not line.startswith("#")) and len(line) != 0:
                samples.add(line)
    return samples


start_on()
samples = get_samples()


rule all:
    input:
        expand("data/result/{sample}_isoform.sam", sample=samples),
        expand("data/result/{sample}_hit.sam", sample=samples),
        "data/result/mirna.bam",
        "data/result/pre.fasta"


rule index_mirna:
    input:
        f"data/resource/{config['spe']}_mirna.sam",
        f"data/resource/{config['spe']}_pre.fasta"
    output:
        "data/result/mirna.bam",
        "data/result/pre.fasta"
    shell:
        "samtools view {input[0]} -b | samtools sort - -o {output[0]} && "
        "samtools index {output[0]} && "
        "cp {input[1]} {output[1]} && "
        "samtools faidx {output[1]}"


rule trim_reads:
    input:
        "data/fastq/{sample}.fastq"
    output:
        "data/trim/{sample}.fastq",
        "data/trim/{sample}_fastp.html",
        "data/trim/{sample}_fastp.json"
    shell:
        "fastp -i {input[0]} -o {output[0]} -h {output[1]} -j {output[2]}"


rule mark_duplicates:
    input: 
        "data/trim/{sample}.fastq"
    output:
        "data/read/{sample}_reads.tsv"
    shell:
        "awk \"NR%4==2\" {input} |sort |uniq -c |nl -s' ' |sed 's/^ \+//g' "
        "|sed 's/ \+/\t/g' >{output}"


rule find_isoforms:
    input:
        "data/read/{sample}_reads.tsv",
        f"data/resource/{config['spe']}_mirna.tsv"
    output:
        "data/result/{sample}_isoform.tsv"
    shell:
        "script/isomir -l {config[max_edit_dist_5p]} -r {config[max_edit_dist_3p]} "
        "-s {input[0]} -m {input[1]} -o {output[0]}"


rule align_pre:
    input:
        "data/read/{sample}_reads.tsv",
        f"data/resource/{config['spe']}_pre.fasta"
    output:
        "data/result/{sample}_hit.tsv"
    script:
        "script/align_pre.R"


rule isoform_to_sam:
    input:
        "data/result/{sample}_isoform.tsv",
        "data/read/{sample}_reads.tsv",
        f"data/resource/{config['spe']}_mibase.tsv",
        f"data/resource/{config['spe']}_pre.fasta"
    output:
        "data/result/{sample}_isoform.sam"
    script:
        "script/isoform_to_sam.R"


rule hit_to_sam:
    input:
        "data/result/{sample}_hit.tsv",
        "data/read/{sample}_reads.tsv",
        f"data/resource/{config['spe']}_pre.fasta"
    output:
        "data/result/{sample}_hit.sam"
    script:
        "script/hit_to_sam.R"

   
