rule samtools_index:
    input:
        bamdedup = "results/DedupReads/{sample}_aln_marked.bam",
    output:
        bamdedupindex = "results/DedupReads/{sample}_aln_marked.bam.bai",
    log:
        "logs/samtools/{sample}_samtools.log"
    conda:
        "envs/samtools.yml"
    shell:
        "samtools index {input.bamdedup}"

rule samtools_faidx:
    input:
        config["reference_path"],
    output:
        samtoolsreferenceindex = config["reference_path"] + ".fai",
    conda:
        "envs/samtools.yml"
    shell:
        "samtools faidx {input}"
