version 1.0


workflow isomir {
    input {
        String script_dir
        String analysis_dir
        String mibase_file
        String fq_file
        String mirna_file
        String pre_file
        String mirna_sam_file
        Int split_line_num = 1000
        Int max_edit_dist_5p = 2
        Int max_edit_dist_3p = 3
        Float min_identity = 0.9
    }

    # call check {
    #     input: 
    #         script_dir = script_dir,
    #         analysis_dir = analysis_dir,
    #         mibase_file = mibase_file,
    #         fq_file = fq_file,
    #         mirna_file = mirna_file,
    #         pre_file = pre_file,
    #         mirna_sam_file = mirna_sam_file
    # }

    call get_sample_name {
        input:
            fq_file = fq_file
    }

    call mk_dir {
        input:
            analysis_dir = analysis_dir,
            sample_name = get_sample_name.sample_name_
    }

    call index_mirna {
        input:
            mirna_sam_file = mirna_sam_file,
            pre_file = pre_file,
            out_dir = mk_dir.out_dir_
    }

    call mark_duplicates {
        input: 
            fq_file = fq_file,
            out_dir = mk_dir.out_dir_
    }

    call split_reads {
        input:
            script_dir = script_dir,
            read_file = mark_duplicates.read_file_,
            split_line_num = split_line_num,
            read_dir = mk_dir.read_dir_
    }
    
    scatter (split_read_file in split_reads.split_read_files_) {
        call find_isoforms {
            input:
                script_dir = script_dir,
                split_read_file = split_read_file,
                mirna_file = mirna_file,
                split_isoform_file = mk_dir.isoform_dir_ + "/" + basename(split_read_file),
                max_edit_dist_5p = max_edit_dist_5p,
                max_edit_dist_3p = max_edit_dist_3p
        }
    }


    scatter (split_read_file in split_reads.split_read_files_) {
        call align_pre {
            input:
                script_dir = script_dir,
                split_read_file = split_read_file,
                pre_file = pre_file,
                split_hit_file = mk_dir.hit_dir_ + "/" + basename(split_read_file),
                min_identity = min_identity
        }
    }

    call merge_isoform {
        input:
            script_dir = script_dir,
            mibase_file = mibase_file,
            pre_file = pre_file,
            split_isoform_files = find_isoforms.split_isoform_file_,
            read_file = mark_duplicates.read_file_,
            out_dir = mk_dir.out_dir_
    }

    call merge_hit {
        input:
            script_dir = script_dir,
            pre_file = pre_file,
            split_hit_files = align_pre.split_hit_file_,
            read_file = mark_duplicates.read_file_,
            out_dir = mk_dir.out_dir_
    }


    output {
        String isoform_bam_file_ = merge_isoform.isoform_bam_file_
        String hit_bam_file_ = merge_hit.hit_bam_file_
        String mirna_bam_file_ = index_mirna.mirna_bam_file_
        String out_pre_file_ = index_mirna.out_pre_file_
    }

}

# ========================================================
# tasks
# ========================================================

task check {
    input {
        String script_dir
        String analysis_dir
        String mibase_file
        String fq_file
        String mirna_file
        String pre_file
        String mirna_sam_file
    }
    
    command <<<
        if [ ! -d "~{script_dir}" ]; then
            echo "script dir does not exist!"
            exit 1
        fi
    
        if [ ! -d "~{analysis_dir}" ]; then
            mkdir analysis_dir
        fi

        if [ ! -f "~{mibase_file}" ]; then
            echo "mibase file does not exist!"
            exit 1
        fi

        if [ ! -f "~{fq_file}" ]; then
            echo "FASTQ file does not exist!"
            exit 1
        fi

        if [ ! -f "~{mirna_file}" ]; then
            echo "miRNA file does not exist!"
            exit 1
        fi

        if [ ! -f "~{pre_file}" ]; then
            echo "pre file does not exist!"
            exit 1
        fi

        if [ ! -f "~{mirna_sam_file}" ]; then
            echo "mirna sam file file does not exist!"
            exit 1
        fi
    >>>
}

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
        String analysis_dir
        String sample_name
    }

    String out_dir = analysis_dir + "/" + sample_name
    String read_dir = out_dir + "/" + "read"
    String isoform_dir = out_dir + "/" + "isoform"
    String hit_dir = out_dir + "/" + "hit"

    command <<<
        if [ ! -d "~{analysis_dir}" ]
        then
            mkdir ~{analysis_dir}
        fi

        rm -rf ~{out_dir} && mkdir ~{out_dir}
        rm -rf ~{read_dir} && mkdir ~{read_dir}
        rm -rf ~{isoform_dir} && mkdir ~{isoform_dir}
        rm -rf ~{hit_dir} && mkdir ~{hit_dir}
    >>>

    output {
        String out_dir_ = out_dir
        String read_dir_ = read_dir
        String isoform_dir_ = isoform_dir
        String hit_dir_ = hit_dir
    }
}

task index_mirna {
    input {
        String mirna_sam_file
        String pre_file
        String out_dir
    }

    String mirna_bam_file = out_dir + "/mirna.bam"
    String out_pre_file = out_dir + "/pre.fa"

    command <<<
        samtools view ~{mirna_sam_file} -b | samtools sort - -o ~{mirna_bam_file}
        samtools index ~{mirna_bam_file}
        cp ~{pre_file} ~{out_pre_file} && samtools faidx ~{out_pre_file}
    >>>

    output {
        String mirna_bam_file_ = mirna_bam_file
        String out_pre_file_ = out_pre_file
    }
}

# mark  duplicates in FASTQ file
task mark_duplicates {
    input {
        String fq_file
        String out_dir
    }

    String read_file = out_dir + "/read.tsv"

    command <<<
        awk "NR%4==2" ~{fq_file} |sort |uniq -c |nl -s' ' |sed 's/^ \+//g' |sed 's/ \+/\t/g' >~{read_file}
    >>>

    output {
        String read_file_ = read_file
    }
}

task split_reads {
    input {
        String script_dir
        String read_file
        Int split_line_num
        String read_dir
    }

    command <<<
        Rscript ~{script_dir}/split_reads.R ~{read_file} ~{split_line_num} ~{read_dir} 
    >>>

    output {
        Array[String] split_read_files_ = glob(read_dir + "/*.tsv")
    }
}

task find_isoforms {
    input {
        String script_dir
        String split_read_file
        String mirna_file
        Int max_edit_dist_5p
        Int max_edit_dist_3p
        String split_isoform_file
    }

    command <<<
        ~{script_dir}/isomir \
        -l ~{max_edit_dist_5p} \
        -r ~{max_edit_dist_3p} \
        -s ~{split_read_file} \
        -m ~{mirna_file} \
        -o ~{split_isoform_file}
    >>>

    output {
        String split_isoform_file_ = split_isoform_file
    }
}

task align_pre {
    input {
        String script_dir
        String split_read_file
        String pre_file
        Float min_identity
        String split_hit_file
    }

    command <<<
        Rscript ~{script_dir}/align_pre.R ~{script_dir} ~{split_read_file} \
        ~{pre_file} ~{min_identity} ~{split_hit_file} 
    >>>

    output {
        String split_hit_file_ = split_hit_file
    }
}

task merge_hit {
    input {
        String script_dir
        String pre_file
        Array[String] split_hit_files
        String read_file
        String out_dir
    }

    String hit_file = out_dir + "/hit.sam"
    String hit_bam_file = out_dir + "/hit.bam"

    command <<<
        Rscript ~{script_dir}/merge_hit.R ~{pre_file} ~{sep="," split_hit_files} \
        ~{read_file} ~{hit_file}
        samtools view ~{hit_file} -b | samtools sort - -o ~{hit_bam_file}
        samtools index ~{hit_bam_file}
    >>>

    output {
        String hit_bam_file_ = hit_bam_file
    }
}

task merge_isoform {
    input {
        String script_dir
        String mibase_file
        String pre_file
        Array[String] split_isoform_files
        String read_file
        String out_dir
    }

    String isoform_file = out_dir + "/isoform.sam"
    String isoform_bam_file = out_dir + "/isoform.bam"

    command <<<
        Rscript ~{script_dir}/merge_isoform.R ~{script_dir} ~{mibase_file} \
                ~{pre_file} ~{sep="," split_isoform_files} ~{read_file} \
                ~{isoform_file}

        samtools view ~{isoform_file} -b | samtools sort - -o ~{isoform_bam_file}
        samtools index ~{isoform_bam_file}
    >>>

    output {
        String isoform_bam_file_ = isoform_bam_file
    }
}
