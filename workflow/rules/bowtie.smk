rule bowtie2Build:
    input:
        config["reference_path"]
    params:
        basename=config["reference_path"]
    output:
        output1=config["reference_path"] + ".1.bt2",
        output2=config["reference_path"] + ".2.bt2",
        output3=config["reference_path"] + ".3.bt2",
        output4=config["reference_path"] + ".4.bt2",
        outputrev1=config["reference_path"] + ".rev.1.bt2",
        outputrev2=config["reference_path"] + ".rev.2.bt2",
    shell: "bowtie2-build {input} {params.basename}"


rule bowtie2Align:
    input:
        r1 = "trimmed/{sample}_R1_001.fastq.gz",
        r2 = "trimmed/{sample}_R2_001.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired = "trimmed/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired = "trimmed/{sample}_R2_unpaired.fastq.gz",
        indexoutput1=config["reference_path"] + ".1.bt2",

    output:
        samout=temp('mapped/{sample}_aln.sam'),
        bamout=temp('mapped/{sample}_aln.bam'),
        bamsortout="mapped/{sample}_aln_sort.bam",

    log:
        "logs/bowtie2/{sample}.log"
    params:
        index=config["reference_path"],  # prefix of reference genome index (built with bowtie2-build)
        # rgid=expand("{sample}", sample=SAMPLES),
        # rgSM=expand("{sample}", sample=SAMPLES),
        # rgid={sample},
        # rgSM={sample},
    threads: 8  # Use at least two threads
    shell:
        # "bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} -U {input.r1_unpaired} -U {input.r2_unpaired} -S {output.samout} -t -p 8 --non-deterministic --end-to-end --rg-id {params.rgid} --rg SM:{params.rgSM} --rg LB:1 --rg PL:Illumina && samtools view -Sb {output.samout} > {output.bamout} && samtools sort -O BAM -o {output.bamsortout} {output.bamout}"
        "bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} -U {input.r1_unpaired} -U {input.r2_unpaired} -S {output.samout} -t -p 8 --non-deterministic --end-to-end && samtools view -Sb {output.samout} > {output.bamout} && samtools sort -O BAM -o {output.bamsortout} {output.bamout}"
