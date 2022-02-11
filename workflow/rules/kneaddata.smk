rule Remove_human:
	input:
		r1 = config["fastq_dir"] + "/{sample}_R1_001.fastq.gz",
		r2 = config["fastq_dir"] + "/{sample}_R2_001.fastq.gz",
	output:
		fq1 = "results/Clean_Human_removed/{sample}_paired_1.fastq.gz",
		fq2 = "results/Clean_Human_removed/{sample}_paired_2.fastq.gz",
	params:
		p = 7,
		db = "resources/kneaddata/",
		outputdir = "results/Clean_Human_removed",
		prefix = "{sample}",
		unzip_fq1 = "results/Clean_Human_removed/{sample}_paired_1.fastq",
		unzip_fq2 = "results/Clean_Human_removed/{sample}_paired_2.fastq",
	threads: 4
 	conda:
		"envs/kneaddata.yml"
	shell:
		"kneaddata -i {input.r1} -i {input.r2} -db {params.db} -t {threads} -p {params.p} -o {params.outputdir} --output-prefix {params.prefix} --trimmomatic /nfs/esnitkin/bin_group/variant_calling_bin/Trimmomatic/ && gzip {params.unzip_fq1} {params.unzip_fq2}"
