import os

configfile: "data/config/config.yaml"

def start_on():
    data_dir = "data"
    temp_dir = os.path.join(data_dir,"temp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    result_dir = os.path.join(data_dir,"result")
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)


def get_samples():
    samples = set()
    with open(os.path.join("data/config",config["samples"]),"rt") as fin:
        for line in fin:
            line = line.strip()
            if not line.startswith("#") and len(line) != 0:
                samples.add(line)
    return samples


start_on()
samples = get_samples()


onsuccess:
    shell("chmod -R 777 data")


rule all:
    input:
        expand("data/result/{sample}_isoform.bam",sample=samples),
        "data/result/mirna.bam",
        "data/result/pre.fasta"


rule format_mibase:
    input:
        "data/resource/miRNA.dat"
    output:
        "data/resource/mibase.tsv"
    script:
        "script/format_mibase.py"


rule ex_mibase:
    input:
        "data/resource/mibase.tsv",
        "data/resource/custom_mibase.tsv",
    output:
        f"data/resource/{config['spe']}_mibase.tsv",
        f"data/resource/{config['spe']}_mirna.tsv",
        f"data/resource/{config['spe']}_pre.fasta",
        f"data/resource/{config['spe']}_mirna.sam"
    script:
        "script/ex_mibase.R"


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
        "data/fastq/{sample}.fastq.gz"
    output:
        "data/temp/{sample}.fastq",
        "data/temp/{sample}_fastp.html",
        "data/temp/{sample}_fastp.json"
    threads:
        config["fastp"]["thread"]
    shell:
        "fastp -i {input[0]} -o {output[0]} -h {output[1]} -j {output[2]} -w {threads}"


rule mark_duplicates:
    input:
        'data/temp/{sample}.fastq'
    output:
        "data/temp/{sample}_reads.tsv"
    shell:
        "awk \"NR%4==2\" {input} |sort |uniq -c |nl -s' ' |sed 's/^ \+//g' "
        "|sed 's/ \+/\t/g' >{output}"


rule find_isoforms:
    input:
        "data/temp/{sample}_reads.tsv",
        f"data/resource/{config['spe']}_mirna.tsv"
    output:
        "data/result/{sample}_isoform.tsv"
    shell:
        "core/detect.py  "
        "-l {config[max_edit_dist_5p]} -r {config[max_edit_dist_3p]} -n {config[min_read_num]} "
        "{input[0]} {input[1]} {output[0]}"


rule isoform_to_sam:
    input:
        "data/result/{sample}_isoform.tsv",
        "data/temp/{sample}_reads.tsv",
        f"data/resource/{config['spe']}_mibase.tsv",
        f"data/resource/{config['spe']}_pre.fasta"
    output:
        "data/result/{sample}_isoform.sam"
    script:
        "script/isoform_to_sam.R"


rule index_isoform:
    input:
        "data/result/{sample}_isoform.sam",
    output:
        "data/result/{sample}_isoform.bam",
    shell:
        "samtools view {input[0]} -b | samtools sort - -o {output[0]} && "
        "samtools index {output[0]}"
