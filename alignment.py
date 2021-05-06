#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess

def trim(threads, work_dir):
    threads=threads
    work_dir=work_dir

    if os.path.isdir(work_dir + "/trim") == False or len(os.listdir(work_dir + "/trim")) == 0:
        fastq_list=glob.glob(work_dir+"/raw-data/*R1_001.fastq.gz")
        for read1 in fastq_list:
            read2=read1.replace("R1","R2")
            trim_galore_command=["trim_galore","-j","4","-o","./trim", "--paired",read1,read2]
            subprocess.run(trim_galore_command)
    elif os.path.isdir(work_dir + "/trim") == True or len(os.listdir(work_dir + "/trim")) > 0:
        print("Skipping adapter trimming (already performed)")

def salmon(salmon_index,threads,work_dir,gtf):
    salmon_index=salmon_index
    threads=threads
    work_dir=work_dir
    gtf=gtf

    if os.path.isdir(work_dir + "/salmon") == False or len(os.listdir(work_dir + "/salmon")) == 0:
        print("Mapping reads with Salmon:")
        trim_list=glob.glob(work_dir+"/trim_galore/*R1_001_val_1.fq.gz")
        for read1 in trim_list:
            print("Mapping" + read1.replace("_R1_001_val_1.fq.gz", ""))
            os.mkdir(work_dir+"/salmon")
            read2=read1.replace("R1_001_val_1.fq.gz", "R2_001_val_2.fq.gz")
            out_file=read1.replace("_R1_001_val_1.fq.gz","")
            salmon_output_file=work_dir+"/salmon/"+out_file
            salmon_command=["salmon","quant","--index",salmon_index,"-l","A",
            "-g", gtf,"-p",threads,"-1", read1,"-2",read2,"--validateMappings",
            "--gcBias","-o", salmon_output_file]

def hisat2():
    print("Mapping reads with HISAT2")
