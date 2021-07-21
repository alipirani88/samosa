configfile: "config/config.yaml"

rule bowtie2Build:
    input:
        config["reference_path"]
    params:
        basename=config["reference_path"]
    output:
        output1="output/reference.1.bt2",
        output2="output/reference.2.bt2",
        output3="output/reference.3.bt2",
        output4="output/reference.4.bt2",
        outputrev1="output/reference.rev1.bt2",
        outputrev2="output/reference.rev2.bt2",
    shell: "bowtie2-build {input} {params.basename}"




# rule bowtie2_build:
#     input:
#         config["reference_path"]
#     output:
#         multiext(".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
#     log:
#         "logs/bowtie2_build/build.log"
#     params:
#         extra=""  # optional parameters
#     threads: 8
#     wrapper:
#         "0.77.0/bio/bowtie2/build"