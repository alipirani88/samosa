rule bwaindex:
    input:
        config["reference_path"]
    params:
        basename=config["reference_path"]
    output:
        amb,ann,bwt,pac,sa
        output1=config["reference_path"] + ".bwt",
        output2=config["reference_path"] + ".amb",
        output3=config["reference_path"] + ".ann",
        output4=config["reference_path"] + ".pac",
        output5=config["reference_path"] + ".sa",
        
    conda:
        "envs/bwa.yml"
    shell: 
        "bwa index {input}"


rule bwamem:
    input:
        r1 = "results/trimmed/{sample}_R1_001.fastq.gz",
        r2 = "results/trimmed/{sample}_R2_001.fastq.gz",
        r1_unpaired = "results/trimmed/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired = "results/trimmed/{sample}_R2_unpaired.fastq.gz",
        indexoutput1=config["reference_path"] + ".bwt",

    output:
        samout=temp('results/mapped/{sample}_aln.sam'),
        samoutunpaired1=temp('results/mapped/{sample}_unpaired1_aln.sam'),
        samoutunpaired2=temp('results/mapped/{sample}_unpaired2_aln.sam'),
        bamout=temp('results/mapped/{sample}_aln.bam'),
        bamsortout="results/mapped/{sample}_aln_sort.bam",

    log:
        "logs/bwa/{sample}.log"
    params:
        index=config["reference_path"],  # prefix of reference genome index (built with bowtie2-build)

    threads: 8  # Use at least two threads
    conda:
        "envs/bwa.yml"
    shell:
        "rgid=`echo {input.r1} | cut -d'/' -f3 | sed 's/_R1_001.fastq.gz//g'` && rgsm=`echo {input.r1} | cut -d'/' -f3 | sed 's/_R1_001.fastq.gz//g'` && bwa mem -M -R '@RG\tID:$rgid\tSM:$rgsm' -t {threads} {params.index} {input.r1} {input.r2} > {output.samout} && bwa mem -M -R '@RG\tID:$rgid\tSM:$rgsm' -t {threads} {params.index} {input.r1_unpaired} >> {output.samout} && bwa mem -M -R '@RG\tID:$rgid\tSM:$rgsm' -t {threads} {params.index} {input.r2_unpaired} >> {output.samout} && samtools view -Sb {output.samout} > {output.bamout} && samtools sort -O BAM -o {output.bamsortout} {output.bamout}"

