rule extract_alleles:
    input:
        rawvcf="results/freebayes/{sample}_raw.vcf",
        gatk_mutect_vcf="results/gatk_mutect/{sample}_Mutect_raw.vcf",
        instrain_vcf="results/instrain/{sample}_instrain/output/{sample}_instrain_SNVs.tsv",
    output:
        LDV_abund_frequency="results/instrain/{sample}_instrain/output/{sample}_LDV_abund_frequency.csv",
    conda:
        "envs/instrain.yml"
    log:
        "logs/instrain/{sample}_instrain.log"
    shell:
        #--skip_plot_generation
        "inStrain profile {input.bamdedup} {input.reference} -o {output.instrain_out} -s {input.reference_stb}"