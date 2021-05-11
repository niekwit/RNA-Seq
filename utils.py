#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os
import subprocess
import yaml

def fastqc(work_dir,threads):
    threads=threads
    work_dir=work_dir
    if not os.path.isdir(work_dir + "/fastqc") or len(os.listdir(work_dir + "/fastqc")) == 0:
        os.makedirs(work_dir+"/fastqc",exist_ok=True)
        fastqc_command="fastqc --threads "+str(threads)+" --quiet -o fastqc/ raw-data/*.fastq.gz"
        multiqc_command=["multiqc","-o","fastqc/","fastqc/"]
        #log commands
        with open(work_dir+"/commands.log","w") as file:
            file.write("FastQC: ")
            print(fastqc_command, file=file)
            file.write("MultiQC: ")
            print(*multiqc_command, sep=" ", file=file)
        print("Running FastQC on raw data")
        subprocess.run(fastqc_command, shell=True)
        print("Running MultiQC")
        subprocess.run(multiqc_command)
    else:
        print("Skipping FastqC/MultiQC (already performed)")

def trim(threads, work_dir):
    threads=threads
    work_dir=work_dir

    #cap threads at 4 for trim_galore
    if int(threads) > 4:
        threads="4"

    if not os.path.isdir(work_dir + "/trim") or len(os.listdir(work_dir + "/trim")) == 0:
        fastq_list=glob.glob(work_dir+"/raw-data/*R1_001.fastq.gz")
        for read1 in fastq_list:
            read2=read1.replace("R1","R2")
            trim_galore_command=["trim_galore","-j",threads,"-o","./trim", "--paired",read1,read2]
            #log commands
            with open(work_dir+"/commands.log", "a") as file:
                file.write("Trim Galore: ")
                print(*trim_galore_command, sep=" ",file=file)
            subprocess.run(trim_galore_command)
    elif os.path.isdir(work_dir + "/trim") == True or len(os.listdir(work_dir + "/trim")) > 0:
        print("Skipping adapter trimming (already performed)")

def salmon(salmon_index,threads,work_dir,gtf,fasta,script_dir,settings):
    salmon_index=salmon_index
    threads=str(threads)
    work_dir=work_dir
    gtf=gtf
    fasta=fasta
    script_dir=script_dir
    settings=settings

    if salmon_index == "": #Salmon index not found, make on the fly
        print("No Salmon index found: generating Salmon index")
        if os.path.isfile(fasta):
            index_dir=os.path.join(script_dir,"salmon-index")
            os.mkdir(index_dir)
            salmon_index_command=["salmon","index","-t",fasta,"-i",index_dir, "--gencode"]
            #log commands
            with open(work_dir+"/commands.log", "a") as file:
                file.write("Salmon index: ")
                print(*salmon_index_command, sep=" ",file=file)

            subprocess.run(salmon_index_command) #run Salmon index

            #Write salmon index file location to settings.yaml
            with open(script_dir+"/settings.yaml") as f:
                doc=yaml.safe_load(f)
            doc["salmon_index"]["gencode-v35"]=index_dir
            with open(script_dir+"/settings.yaml", "w") as f:
                yaml.dump(doc,f)
        else:
            print("ERROR: no FASTA file specified in settings.yaml")

    if not os.path.isdir(work_dir + "/salmon") or len(os.listdir(work_dir + "/salmon")) == 0:
        print("Mapping reads with Salmon:")
        trim_list=glob.glob(work_dir+"/trim/*R1_001_val_1.fq.gz")
        for read1 in trim_list:
            print("Mapping sample " + read1.replace("_R1_001_val_1.fq.gz", ""))
            os.makedirs(work_dir+"/salmon", exist_ok=True)
            read2=read1.replace("R1_001_val_1.fq.gz", "R2_001_val_2.fq.gz")
            out_file=os.path.basename(read1.replace("_R1_001_val_1.fq.gz",""))
            salmon_output_file=work_dir+"/salmon/"+out_file+"-quant"
            salmon_index=settings["salmon_index"]["gencode-v35"] #reload index
            salmon_command=["salmon","quant","--index",salmon_index,"-l","A",
            "-g", gtf,"-p",threads,"-1", read1,"-2",read2,"--validateMappings",
            "--gcBias","-o", salmon_output_file]
            with open(work_dir+"/commands.log", "a") as file:
                file.write("Salmon quant: ")
                print(*salmon_command, sep=" ",file=file)
            subprocess.run(salmon_command) #Run Salmon quant
    else:
        print("Skipping Salmon (already performed)")

def hisat2():
    print("Mapping reads with HISAT2 (UNDER CONSTRUCTION)")

def diff_expr(work_dir,gtf,script_dir):
    work_dir=work_dir
    gtf=gtf
    script_dir=script_dir

    deseq2_command=["Rscript",script_dir+"/deseq2.R",work_dir,gtf,script_dir]
    with open(work_dir+"/commands.log", "a") as file:
        file.write("DESeq2: ")
        print(*deseq2_command, sep=" ",file=file)
    print("Running differential expression analysis with DESeq2")
    os.makedirs(work_dir+"/DESeq2",exist_ok=True)
    subprocess.run(deseq2_command)

def go():
    print("Running GO analysis with DAVID")
