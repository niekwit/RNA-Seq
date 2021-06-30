#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os
import subprocess
import yaml
import sys
import pkg_resources

def install_packages(): #check for required python packages; installs if absent
    required = {"pyyaml"}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        python = sys.executable
        print("Installing missing required Python3 packages")
        subprocess.check_call([python, '-m', 'pip3', 'install', *missing], stdout=subprocess.DEVNULL)

def file_exists(file): #check if file exists/is not size zero
    if os.path.exists(file):
            if os.path.getsize(file) > 0:
                print("Skipping "+file+" (already exists/analysed)")
                return(True)
    else:
        return(False)

def write2log(work_dir,command,name):
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep=" ",file=file)

def fastqc(work_dir,threads):
    threads=threads
    work_dir=work_dir
    if not os.path.isdir(os.path.join(work_dir,"fastqc")) or len(os.listdir(os.path.join(work_dir,"fastqc"))) == 0:
        os.makedirs(work_dir+"/fastqc",exist_ok=True)
        fastqc_command="fastqc --threads "+str(threads)+" --quiet -o fastqc/ raw-data/*.fastq.gz"
        multiqc_command=["multiqc","-o","fastqc/","fastqc/"]
        #log commands
        with open(os.path.join(work_dir,"commands.log"),"w") as file:
            file.write("FastQC: ")
            print(fastqc_command, file=file)
            file.write("MultiQC: ")
            print(*multiqc_command, sep=" ", file=file)
        print("Running FastQC on raw data")
        subprocess.run(fastqc_command, shell=True)
        print("Running MultiQC")
        subprocess.run(multiqc_command)
    else:
        print("Skipping FastQC/MultiQC (already performed)")

def trim(threads, work_dir):
    #cap threads at 4 for trim_galore
    if int(threads) > 4:
        threads="4"

    print("Trimming fastq.gz files")
    fastq_list=glob.glob(work_dir+"/raw-data/*R1_001.fastq.gz")
    for read1 in fastq_list:
        out_dir=os.path.dirname(read1)
        out_dir=out_dir.replace("raw-data","trim")
        out_file1=read1.split(".",1)[0]+"_val_1.fq.gz"
        out_file1=os.path.basename(out_file1)
        out_file1=os.path.join(out_dir,out_file1)
        if not file_exists(out_file1):
            read2=read1.replace("R1","R2")
            trim_galore_command=["trim_galore","-j",threads,"-o","./trim", "--paired",read1,read2]
            #log commands
            with open(work_dir+"/commands.log", "a") as file:
                file.write("Trim Galore: ")
                print(*trim_galore_command, sep=" ",file=file)
            subprocess.run(trim_galore_command)

def salmon(salmon_index,threads,work_dir,gtf,fasta,script_dir,settings):
    if salmon_index == "": #Salmon index not found, make on the fly
        print("No Salmon index found: generating Salmon index")
        if os.path.isfile(fasta):
            index_dir=os.path.join(script_dir,"salmon-index")
            os.mkdir(index_dir)
            salmon_index_command=["salmon","index","-t",fasta,"-i",index_dir, "--gencode"]
            #log commands
            with open(os.path.join(work_dir,"commands.log"), "a") as file:
                file.write("Salmon index: ")
                print(*salmon_index_command, sep=" ",file=file)

            subprocess.run(salmon_index_command) #run Salmon index

            #Write salmon index file location to settings.yaml
            with open(os.path.join(script_dir,"settings.yaml")) as f:
                doc=yaml.safe_load(f)
            doc["salmon_index"]["gencode-v35"]=index_dir
            with open(os.path.join(script_dir,"settings.yaml"), "w") as f:
                yaml.dump(doc,f)
        else:
            print("ERROR: no FASTA file specified in settings.yaml")
            sys.exit()

    print("Mapping reads with Salmon:")
    salmon_output_dir=os.path.join(work_dir,"salmon")
    os.makedirs(salmon_output_dir, exist_ok=True)

    trim_list=glob.glob(os.path.join(work_dir,"trim/*_R1_001_val_1.fq.gz"))
    for read1 in trim_list:
        base_read1=os.path.basename(read1).replace("_R1_001_val_1.fq.gz","")+"-quant"
        salmon_folder_test=os.path.join(salmon_output_dir,base_read1)
        if not file_exists(salmon_folder_test):
            print("Mapping sample " + read1.replace("_R1_001_val_1.fq.gz", ""))
            read2=read1.replace("R1_001_val_1.fq.gz", "R2_001_val_2.fq.gz")
            out_file=os.path.basename(read1.replace("_R1_001_val_1.fq.gz",""))
            salmon_output_file=os.join.path(work_dir,"salmon",out_file)+"-quant"
            salmon_index=settings["salmon_index"]["gencode-v35"] #reload index
            salmon_command=["salmon","quant","--index",salmon_index,"-l","A",
            "-g", gtf,"-p",threads,"-1", read1,"-2",read2,"--validateMappings",
            "--gcBias","-o", salmon_output_file]
            with open(os.join.path(work_dir,"commands.log"), "a") as file:
                file.write("Salmon quant: ")
                print(*salmon_command, sep=" ",file=file)
            subprocess.run(salmon_command) #Run Salmon quant

def hisat2():
    print("Mapping reads with HISAT2 (UNDER CONSTRUCTION)")
    sys.exit()

def diff_expr(work_dir,gtf,script_dir,species):
    deseq2_command=["Rscript",script_dir+"/deseq2.R",work_dir,gtf,script_dir,species]
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write("DESeq2: ")
        print(*deseq2_command, sep=" ",file=file)
    print("Running differential expression analysis with DESeq2")
    os.makedirs(work_dir+"/DESeq2",exist_ok=True)
    subprocess.run(deseq2_command)

def go():
    print("Running GO analysis with DAVID")
