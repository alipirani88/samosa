rule instrain:
    input:
        bamdedup = "results/DedupReads/{sample}_aln_marked.bam",
        bamdedupbai = "results/DedupReads/{sample}_aln_marked.bam.bai",
        samtoolsreferenceindex = config["reference_path"] + ".fai",
        reference=config["reference_path"],
        reference_stb=config["reference_stb"],
    output:
        instrain_out=directory("results/instrain/{sample}_instrain"),
	instrain_outfile="results/instrain/{sample}_instrain/output/{sample}_instrain_SNVs.tsv",
    conda:
        "envs/instrain.yml"
    log:
        "logs/instrain/{sample}_instrain.log"
    shell:
        #--skip_plot_generation
        "inStrain profile {input.bamdedup} {input.reference} -o {output.instrain_out} -s {input.reference_stb} --pairing_filter non_discordant --detailed_mapping_info &>{log}"
