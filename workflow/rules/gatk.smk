rule GATK_MarkDuplicates:
    input:
        bamsortout="mapped/{sample}_aln_sort.bam",
    output:
        bamdedup="DedupReads/{sample}_aln_marked.bam",
        bamdedupmetrics="DedupReads/{sample}_marked_dup_metrics.txt",
    log:
        "logs/gatk/{sample}_GATK_MarkDuplicates.log"
    shell:
        "picard MarkDuplicates REMOVE_DUPLICATES=true I={input.bamsortout} O={output.bamdedup} M={output.bamdedupmetrics} 2>{log}"
