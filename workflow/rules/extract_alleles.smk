rule extract_alleles:
    input:
        freebayescf="results/freebayes/{sample}_raw.vcf",
        gatk_mutect_vcf="results/gatk_mutect/{sample}_Mutect_raw.vcf",
        #instrain_vcf="results/instrain/{sample}_instrain/output/{sample}_instrain_SNVs.tsv",
	instrain_out=directory("results/instrain/{sample}_instrain"),
	instrain_outfile="results/instrain/{sample}_instrain/output/{sample}_instrain_SNVs.tsv",	
    output:
        LDV_abund_frequency="results/extract_allele/{sample}_LDV_abund_frequency.csv",
    conda:
        "envs/extract_allele.yml"
    log:
        "logs/extract_allele/{sample}_instrain.log"
    shell:
        #--skip_plot_generation
        "python bin/extract_alleles.py -Freebayes_vcf {input.freebayescf} -Mutect_vcf {input.gatk_mutect_vcf} -instrain_vcf {input.instrain_outfile} -LDV_abund_frequency {output.LDV_abund_frequency}"
