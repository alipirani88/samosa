rule freebayes:
    input:
        bamdedup = "results/DedupReads/{sample}_aln_marked.bam",
        bamdedupbai = "results/DedupReads/{sample}_aln_marked.bam.bai",
        samtoolsreferenceindex = config["reference_path"] + ".fai",
        reference=config["reference_path"],
    output:
        rawvcf="results/freebayes/{sample}_raw.vcf",
    log:
        "logs/freebayes/{sample}_freebayes.log"
    shell:
        "freebayes-parallel <(fasta_generate_regions.py {input.samtoolsreferenceindex} 100000) 4 -f {input.reference} --min-alternate-count 2 --min-alternate-fraction 0 --ploidy 1 --pooled-continuous --report-monomorphic {input.bamdedup} > {output.rawvcf}"
