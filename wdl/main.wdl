version 1.0


workflow isomir {
    input {
        String script_dir
        String result_dir
        String mibase_file
        String fq_file
        String mirna_file
        String pre_file
        String mirna_sam_file
        Int split_num = 10
        Int max_edit_dist_5p = 2
        Int max_edit_dist_3p = 3
        Float min_identity = 0.9
        Boolean is_pre = false
    }

    call get_sample_name {
        input:
            fq_file = fq_file
    }

    call mk_dir {
        input:
            result_dir = result_dir,
            sample_name = get_sample_name.sample_name_
    }

    call index_mirna {
        input:
            mirna_sam_file = mirna_sam_file,
            pre_file = pre_file,
            sample_dir = mk_dir.sample_dir_
    }

    call mark_duplicates {
        input: 
            fq_file = fq_file,
    }

    call split_reads {
        input:
            script_dir = script_dir,
            read_file = mark_duplicates.read_file_,
            split_num = split_num
    }
    
    scatter (split_read_file in split_reads.split_read_files_) {
        call find_isoforms {
            input:
                script_dir = script_dir,
                split_read_file = split_read_file,
                mirna_file = mirna_file,              
                max_edit_dist_5p = max_edit_dist_5p,
                max_edit_dist_3p = max_edit_dist_3p
        }
    }

    call merge_isoform {
        input:
            script_dir = script_dir,
            mibase_file = mibase_file,
            pre_file = pre_file,
            split_isoform_files = find_isoforms.split_isoform_file_,
            read_file = mark_duplicates.read_file_,
            sample_dir = mk_dir.sample_dir_
    }

    if (is_pre) {
        scatter (split_read_file in split_reads.split_read_files_) {
            call align_pre {
                input:
                    script_dir = script_dir,
                    split_read_file = split_read_file,
                    pre_file = pre_file,
                    min_identity = min_identity
            }
        }

        call merge_hit {
            input:
                script_dir = script_dir,
                pre_file = pre_file,
                split_hit_files = align_pre.split_hit_file_,
                read_file = mark_duplicates.read_file_,
                sample_dir = mk_dir.sample_dir_
        }
    }    

    output {
        File isoform_bam_file_ = merge_isoform.isoform_bam_file_
        File? hit_bam_file_ = merge_hit.hit_bam_file_
        File mirna_bam_file_ = index_mirna.mirna_bam_file_
        File out_pre_file_ = index_mirna.out_pre_file_
    }
}

# ========================================================
# tasks
# ========================================================
task get_sample_name {
    input {
        String fq_file
    }

    command <<<
        file_name=$(basename ~{fq_file})
        echo ${file_name%%.*}
    >>>

    output {
        String sample_name_ = read_string(stdout())
    }
}

task mk_dir {
    input {
        String result_dir
        String sample_name
    }

    String sample_dir = result_dir + "/" + sample_name

    command <<<
        if [ ! -d "~{result_dir}" ]
        then
            mkdir ~{result_dir}
        fi

        rm -rf ~{sample_dir} && mkdir ~{sample_dir}
    >>>

    output {
        File sample_dir_ = sample_dir
    }
}

task index_mirna {
    input {
        String mirna_sam_file
        String pre_file
        String sample_dir
    }

    String mirna_bam_file = sample_dir + "/mirna.bam"
    String out_pre_file = sample_dir + "/pre.fa"

    command <<<
        samtools view ~{mirna_sam_file} -b | samtools sort - -o ~{mirna_bam_file}
        samtools index ~{mirna_bam_file}
        cp ~{pre_file} ~{out_pre_file} && samtools faidx ~{out_pre_file}
    >>>

    output {
        File mirna_bam_file_ = mirna_bam_file
        File out_pre_file_ = out_pre_file
    }
}

# mark  duplicates in FASTQ file
task mark_duplicates {
    input {
        String fq_file
    }

    String read_file = "read.tsv"

    command <<<
        awk "NR%4==2" ~{fq_file} |sort |uniq -c |nl -s' ' |sed 's/^ \+//g' |sed 's/ \+/\t/g' >~{read_file}
    >>>

    output {
        File read_file_ = read_file
    }
}

task split_reads {
    input {
        String script_dir
        String read_file
        Int split_num
    }

    String read_dir = "read"

    command <<<
        mkdir ~{read_dir}
        Rscript ~{script_dir}/split_reads.R ~{read_file} ~{split_num} ~{read_dir} 
    >>>

    output {
        Array[File] split_read_files_ = glob(read_dir + "/*.tsv")
    }
}

task find_isoforms {
    input {
        String script_dir
        String split_read_file
        String mirna_file
        Int max_edit_dist_5p
        Int max_edit_dist_3p
    }

    String split_isoform_file = "isoform.tsv"

    command <<<
        ~{script_dir}/isomir \
        -l ~{max_edit_dist_5p} \
        -r ~{max_edit_dist_3p} \
        -s ~{split_read_file} \
        -m ~{mirna_file} \
        -o ~{split_isoform_file}
    >>>

    output {
        File split_isoform_file_ = split_isoform_file
    }
}

task align_pre {
    input {
        String script_dir
        String split_read_file
        String pre_file
        Float min_identity
    }

    String split_hit_file = "hit.tsv"

    command <<<
        Rscript ~{script_dir}/align_pre.R ~{script_dir} ~{split_read_file} \
        ~{pre_file} ~{min_identity} ~{split_hit_file} 
    >>>

    output {
        File split_hit_file_ = split_hit_file
    }
}


task merge_isoform {
    input {
        String script_dir
        String mibase_file
        String pre_file
        Array[String] split_isoform_files
        String read_file
        String sample_dir
    }

    String isoform_file = sample_dir + "/isoform.sam"
    String isoform_bam_file = sample_dir + "/isoform.bam"

    command <<<
        Rscript ~{script_dir}/merge_isoform.R ~{script_dir} ~{mibase_file} \
                ~{pre_file} ~{sep="," split_isoform_files} ~{read_file} \
                ~{isoform_file}

        samtools view ~{isoform_file} -b | samtools sort - -o ~{isoform_bam_file}
        samtools index ~{isoform_bam_file}
    >>>

    output {
        File isoform_bam_file_ = isoform_bam_file
    }
}


task merge_hit {
    input {
        String script_dir
        String pre_file
        Array[String] split_hit_files
        String read_file
        String sample_dir
    }

    String hit_file = sample_dir + "/hit.sam"
    String hit_bam_file = sample_dir + "/hit.bam"

    command <<<
        Rscript ~{script_dir}/merge_hit.R ~{pre_file} ~{sep="," split_hit_files} \
        ~{read_file} ~{hit_file}
        samtools view ~{hit_file} -b | samtools sort - -o ~{hit_bam_file}
        samtools index ~{hit_bam_file}
    >>>

    output {
        File hit_bam_file_ = hit_bam_file
    }
}