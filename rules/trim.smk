rule trimmomatic_pe:
    input:
        r1 = "results/Clean_Human_removed/{sample}_paired_1.fastq.gz",
        r2 = "results/Clean_Human_removed/{sample}_paired_2.fastq.gz",
    output:
        r1 = "results/trimmed/{sample}_R1_001.fastq.gz",
        r2 = "results/trimmed/{sample}_R2_001.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired = "results/trimmed/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired = "results/trimmed/{sample}_R2_unpaired.fastq.gz",
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        adapter_filepath=config["adapter_file"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"],

    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}"