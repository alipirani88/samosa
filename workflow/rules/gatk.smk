rule GATK_MarkDuplicates:
    input:
        bamsortout="results/mapped/{sample}_aln_sort.bam",
    output:
        bamdedup="results/DedupReads/{sample}_aln_marked.bam",
        bamdedupmetrics="results/DedupReads/{sample}_marked_dup_metrics.txt",
    log:
        "logs/gatk/{sample}_GATK_MarkDuplicates.log"
    shell:
        "picard MarkDuplicates REMOVE_DUPLICATES=true I={input.bamsortout} O={output.bamdedup} M={output.bamdedupmetrics} 2>{log}"

rule picard_SeqDict:
    input:
        reference=config["reference_path"],
    output:
        ref_seqdict=reference_seqdict,
    shell:
        "picard CreateSequenceDictionary R={input.reference} O={output.ref_seqdict}"


rule GATK_DepthofCoverage:
    input:
        bamdedup = "results/DedupReads/{sample}_aln_marked.bam",
        reference=config["reference_path"],
        ref_seqdict=reference_seqdict,
        bamdedupbai = "results/DedupReads/{sample}_aln_marked.bam.bai",

    output:
        gatk_depth_summary="results/DedupReads/{sample}_depth_of_coverage",
    log:
        "logs/gatk/{sample}_GATK_DepthofCoverage.log"
    shell:
        "java -Xmx5G -jar ./bin/GenomeAnalysisTK.jar -T DepthOfCoverage -R {input.reference} -o {output.gatk_depth_summary} -I {input.bamdedup} --summaryCoverageThreshold 1 --summaryCoverageThreshold 2 --summaryCoverageThreshold 3 --summaryCoverageThreshold 4 --summaryCoverageThreshold 5 --summaryCoverageThreshold 6 --summaryCoverageThreshold 7 --summaryCoverageThreshold 8 --summaryCoverageThreshold 9 --summaryCoverageThreshold 10 --summaryCoverageThreshold 15 --summaryCoverageThreshold 20 --summaryCoverageThreshold 25 --ignoreDeletionSites 2>{log}"
        # "java -Xmx5G -jar /nfs/esnitkin/bin_group/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar DepthOfCoverage -R {input.reference} -O {output.gatk_depth_summary} -I {input.bamdedup} --summary-coverage-threshold 1 --summary-coverage-threshold 2 --summary-coverage-threshold 3 --summary-coverage-threshold 4 --summary-coverage-threshold 5 --summary-coverage-threshold 6 --summary-coverage-threshold 7 --summary-coverage-threshold 8 --summary-coverage-threshold 9 --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 25 --ignore-deletion-sites true --intervals true 2>{log}"


rule GATK_HaplotypeCaller:
    input:
        bamdedup = "results/DedupReads/{sample}_aln_marked.bam",
        reference=config["reference_path"],
        ref_seqdict=reference_seqdict,
        bamdedupbai = "results/DedupReads/{sample}_aln_marked.bam.bai",

    output:
        gatk_vcf="results/gatk/{sample}_raw.vcf",
    log:
        "logs/gatk/{sample}_GATK_HaplotypeCaller.log"
    shell:
        "java -Xmx5G -jar /nfs/esnitkin/bin_group/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar HaplotypeCaller -R {input.reference} -O {output.gatk_vcf} -I {input.bamdedup} 2>{log}"