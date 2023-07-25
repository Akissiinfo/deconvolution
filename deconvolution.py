
from Bio import SeqIO

def getFile_name():
        records =[]
        for seqrecord in SeqIO.parse(config["index"], "fasta"):
            records.append(seqrecord.id)
            File_list = records

def getInput():
    my_file = getFile_name()
    return my_file
#files = glob_wildcards("data/{file}.fastq.gz")
configfile: "config.yaml"

#sample = config["sample"]

from glob import glob

rule all:
    input: expand("{file}.fastq", file=getInput())

#rule all:
    #input:
        #"barcode_found.output.result"

rule merge_fastq:
    input:
        #files = glob("data_1/*")
        files = config["data"]
    output:
        "all_sequences.fastq"
    shell:
        "cat {input.files} > {output}"

rule filter_all_sequence:
    input:
        rules.merge_fastq.output
    output:
        "filter_all_sequence.fastq"
    shell:
        "seqkit seq -m 1000 {input} > {output}"
         #define amplicon length

rule barcode_index:
    input:
        #index = expand("{sample}", sample = config["samples"]),
        index = config["index"],
        seq =  rules.filter_all_sequence.output
    output:
        "index_result.txt"
    log:
        "/logs/index.log"
    shell:
        "seqkit locate --pattern-file {input.index} -m 0 {input.seq} > {output} 2> {log}"

rule barcode_forward:
    input:
        barcodeforward = config["forward"],
        sequence = rules.filter_all_sequence.output
    output:
        "forward_result.txt"
    log:
        "/logs/forward.log"
    shell:
        "seqkit locate  --degenerate --pattern-file {input.barcodeforward} {input.sequence} > {output} 2> {log}"

rule barcode_reverse:
    input:
       barcodereverse = config["reverse"],
        sequence1 = rules.filter_all_sequence.output
    output:
        "reverse_result.txt"
    log:
        "/logs/reverse.log"
    shell:
         "seqkit locate  --degenerate --pattern-file {input.barcodereverse}  {input.sequence1} > {output} 2> {log}"

rule barcode_found:
    input:
        ind = rules.barcode_index.output,#"index_result.txt",
        forw = rules.barcode_forward.output,#"forward_result.txt",
        rev = rules.barcode_reverse.output,#"reverse_result.txt",#rules.barcode_reverse.output,
        seqs = rules.filter_all_sequence.output,#"filter_all_sequence.fastq",
        script = "pipeline_of_deconvolution.py"
    output:
      touch(expand("{file}.csv", file=getInput()))
    shell:
      "python3 {input.script} -i {input.ind} -f {input.forw} -r {input.rev} -t {input.seqs} > {output}"

rule split_reads:
  input:
    rules.barcode_found.output   
  output:
     touch(expand("{file}.txt", file=getInput())) 
  shell:
   "./script_split.sh {input} > {output}"


rule catch_reads:
  input:
    filename = rules.split_reads.output,  
    target =  rules.filter_all_sequence.output
  output: 
    touch(expand("{file}.fastq", file=getInput())) #"RPI1.fastq"

  shell:
    #"grep -f {input.filename} -A 3 {input.target} > {output}"
    #"for file in {input.filename}; do if [ !-s $file ]; then grep -f $file -A 3 {input.target}; fi; done > {output}"
   #"for i in {input.filename}; do grep -f $i -A 3 {input.target}; done" #> {output}"
     "./script_catch.sh {input.filename}  {input.target} > {output}"
