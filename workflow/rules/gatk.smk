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


rule GATK_HaplotypeCaller:
    input:
        bamdedup = "results/DedupReads/{sample}_aln_marked.bam",
        reference=config["reference_path"],
        ref_seqdict=reference_seqdict,
        bamdedupbai = "results/DedupReads/{sample}_aln_marked.bam.bai",

    output:
        gatk_vcf="results/gatk/{sample}_HaplotypeCaller_raw.vcf",
        gatk_haplotypecaller_bam="results/gatk/{sample}_HaplotypeCaller.bam",
    log:
        "logs/gatk/{sample}_GATK_HaplotypeCaller.log"
    shell:
        "java -Xmx5G -jar /nfs/esnitkin/bin_group/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar HaplotypeCaller -R {input.reference} -O {output.gatk_vcf} -I {input.bamdedup} -ploidy 1 --annotate-with-num-discovered-alleles true --annotation-group AlleleSpecificAnnotation --annotation AlleleFraction -ERC BP_RESOLUTION --bam-output {output.gatk_haplotypecaller_bam} 2>{log}"

rule GATK_Mutect2:
    input:
        bamdedup = "results/DedupReads/{sample}_aln_marked.bam",
        reference=config["reference_path"],
        ref_seqdict=reference_seqdict,
        bamdedupbai = "results/DedupReads/{sample}_aln_marked.bam.bai",

    output:
        gatk_mutect_vcf="results/gatk_mutect/{sample}_Mutect_raw.vcf",
        gatk_mutect_bam="results/gatk_mutect/{sample}_HaplotypeCaller.bam",
        #gatk_mutect_alleles_vcf="results/gatk_mutect/{sample}_Mutect_alleles.vcf",
    log:
        "logs/gatk/{sample}_GATK_Mutect.log"
    shell:
        "java -Xmx5G -jar /nfs/esnitkin/bin_group/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar Mutect2 -R {input.reference} -O {output.gatk_mutect_vcf} -I {input.bamdedup} -af-of-alleles-not-in-resource 0.33 --annotation-group AlleleSpecificAnnotation --annotation AlleleFraction -ERC BP_RESOLUTION --bam-output {output.gatk_mutect_bam} 2>{log}"


rule GATK_MergeVCF:
    input:
        rawvcf="results/freebayes/{sample}_raw.vcf",
        gatk_vcf="results/gatk/{sample}_HaplotypeCaller_raw.vcf",
        gatk_mutect_vcf="results/gatk_mutect/{sample}_Mutect_raw.vcf",

    output:
        mergedvcf="results/MergedVCF/{sample}_all_merged.vcf"

    log:
        "logs/gatk/{sample}_GATK_MergeVCF.log"
    shell:
        "java -Xmx5G -jar /nfs/esnitkin/bin_group/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar MergeVcfs --INPUT {input.rawvcf} --INPUT {input.gatk_vcf} --INPUT {input.gatk_mutect_vcf} --OUTPUT {output.mergedvcf} -R /scratch/esnitkin_root/esnitkin/apirani/Project_VRE_metagenomics_analysis/2021_04_09_VREfm_variant_calling/Reference_genome/Aus0004//Aus0004.fasta --CREATE_INDEX 2>{log}"
