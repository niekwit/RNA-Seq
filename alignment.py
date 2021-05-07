#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
#import urllib.request
#from pathlib import Path
#import stat

def trim(threads, work_dir):
    threads=threads
    work_dir=work_dir

    if not os.path.isdir(work_dir + "/trim") or len(os.listdir(work_dir + "/trim")) == 0:
        fastq_list=glob.glob(work_dir+"/raw-data/*R1_001.fastq.gz")
        for read1 in fastq_list:
            read2=read1.replace("R1","R2")
            trim_galore_command=["trim_galore","-j","4","-o","./trim", "--paired",read1,read2]
            print(*trim_galore_command, sep=" ")
            subprocess.run(trim_galore_command)
    elif os.path.isdir(work_dir + "/trim") == True or len(os.listdir(work_dir + "/trim")) > 0:
        print("Skipping adapter trimming (already performed)")

def salmon(salmon_index,threads,work_dir,gtf,fasta,script_dir):
    salmon_index=salmon_index
    threads=str(threads)
    work_dir=work_dir
    gtf=gtf
    fasta=fasta
    script_dir=script_dir

    if salmon_index == "None": #Salmon index not found, make on the fly
        print("Generating Salmon index")
        if os.path.isfile(fasta):
            index_name=fasta.replace(".fa", "")
            salmon_index_command=["salmon","index","-t",fasta,"-i","salmon_index/",index_name]
            print(*salmon_index_command, sep=" ")
            os.chdir(script_dir)
            os.mkdir("salmon_index")
            subprocess.run(salmon_index_command)
            os.chdir(work_dir)
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
            salmon_command=["salmon","quant","--index",salmon_index,"-l","A",
            "-g", gtf,"-p",threads,"-1", read1,"-2",read2,"--validateMappings",
            "--gcBias","-o", salmon_output_file]
            print(*salmon_command, sep=" ")
            subprocess.run(salmon_command)


def hisat2():
    print("Mapping reads with HISAT2")
