#files = glob_wildcards("data/{file}.fastq.gz")
configfile: "config.yaml"

#sample = config["sample"]

from glob import glob

rule all:
    input:
        "index_result.txt"

rule merge_fastq:
    input:
        files = glob("data/*")

    output:
        "all_sequences.fastq"
    shell:
        "zcat {input.files} > {output}"

rule filter_all_sequence:
    input:
        "all_sequences.fastq"
    output:
        "filter_all_sequence.fastq"
    shell:
        "seqkit seq -m 700 -M 1200 {input} > {output}"

rule barcode_index:
    input:
        #index = expand("{sample}", sample = config["samples"]),
        index = config["index"],
        seq =  "filter_all_sequence.fastq"
    output:
        "index_result.txt"
    shell:
        "seqkit locate --pattern-file {input.index} -m 0 {input.seq} > {output}"
