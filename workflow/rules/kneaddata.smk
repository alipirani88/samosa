rule Remove_human:
	input:
		r1 = config["fastq_dir"] + "/{sample}_R1.fastq.gz",
		r2 = config["fastq_dir"] + "/{sample}_R2.fastq.gz",
	output:
		fq1 = "results/Clean_Human_removed/{sample}_paired_1.fastq.gz",
		fq2 = "results/Clean_Human_removed/{sample}_paired_2.fastq.gz",
	params:
		db = config["hgbtdb"],
		outputdir = "results/Clean_Human_removed",
	threads: 4
 	conda:
		"envs/kneaddata.yml"
	shell:
		#"kneaddata -i {input.r1} -i {input.r2} -db {params.db} -t {threads} -p {params.p} -o {params.outputdir} --output-prefix {params.prefix} && gzip {params.unzip_fq1} {params.unzip_fq2}"
        #"trimmomatic_bin=`which trimmomatic` && trimmomatic_dir=`dirname $trimmomatic_bin` && trimmomatic_conda=`dirname $trimmomatic_dir` && kneaddata -i {input.r1} -i {input.r2} -db {params.db} -t {threads} -p {params.p} -o {params.outputdir} --output-prefix {params.prefix} --trimmomatic $trimmomatic_conda/share/trimmomatic/ && gzip {params.unzip_fq1} {params.unzip_fq2}"
		bowtie2 -p {threads} -x {params.db} -1 /nfs/esnitkin/Project_VRE_metagenomics_analysis/Sequence_data/fastq/3381-CB/fastqs_3381-CB/MBR04877_R1_001.fastq.gz -2 /nfs/esnitkin/Project_VRE_metagenomics_analysis/Sequence_data/fastq/3381-CB/fastqs_3381-CB/MBR04877_R2_001.fastq.gz --un-conc-gz results/Clean_Human_removed_bowtie/MBR04877_pair --un-gz results/Clean_Human_removed_bowtie/MBR04877_unpair -S /tmp/MBR04877_aln.sam --non-deterministic --end-to-end