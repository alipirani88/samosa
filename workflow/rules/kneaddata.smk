rule Remove_human:
	input:
		r1 = config["fastq_dir"] + "/{sample}_R1.fastq.gz",
		r2 = config["fastq_dir"] + "/{sample}_R2.fastq.gz",
	output:
		fq1 = "results/Clean_Human_removed/{sample}.1.paired",
		fq2 = "results/Clean_Human_removed/{sample}.2.paired",
		unp = "results/Clean_Human_removed/{sample}.unpaired",
		samout=temp('results/Clean_Human_removed/{sample}_aln.sam'),
	params:
		db = config["hgbtdb"],
		outputdir = "results/Clean_Human_removed",
		prefix = "results/Clean_Human_removed/{sample}.paired",
		threads = config["ncores"],
	conda:
		"envs/kneaddata.yml"
	log:
		"logs/kneaddata/{sample}.log"
	shell:
		#"kneaddata -i {input.r1} -i {input.r2} -db {params.db} -t {threads} -p {params.p} -o {params.outputdir} --output-prefix {params.prefix} && gzip {params.unzip_fq1} {params.unzip_fq2}"
        #"trimmomatic_bin=`which trimmomatic` && trimmomatic_dir=`dirname $trimmomatic_bin` && trimmomatic_conda=`dirname $trimmomatic_dir` && kneaddata -i {input.r1} -i {input.r2} -db {params.db} -t {threads} -p {params.p} -o {params.outputdir} --output-prefix {params.prefix} --trimmomatic $trimmomatic_conda/share/trimmomatic/ && gzip {params.unzip_fq1} {params.unzip_fq2}"
		"bowtie2 -p {params.threads} -x {params.db} -1 {input.r1} -2 {input.r2} --un-conc {params.prefix} --un {output.unp} -S {output.samout} --non-deterministic --end-to-end &>{log}"
