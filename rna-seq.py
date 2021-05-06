#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import multiprocessing
import yaml
import timeit
import time
import glob

def main():
    start = timeit.default_timer()#initiate timing of run

    script_dir=os.path.abspath(os.path.dirname(__file__))
    work_dir=os.getcwd()

    sys.path.append(script_dir)#adds script directory to runtime (for importing modules)

    ###loads available RNA-Seq settings
    if os.path.exists(script_dir+"settings.yaml") == True:
        with open(script_dir+"/settings.yaml") as file: settings=yaml.full_load(file)
    else:
        print("ERROR: settings.yaml not found in analysis folder. Please provide this file for further analysis.")
        sys.exit()

    genome_list=["gencode-v35"]#######

    ###command line argument parser
    ap = argparse.ArgumentParser()
    ap.add_argument("-t", "--threads",
                    required=False,
                    default=1,
                    metavar="<int>",
                    help="Number of CPU threads to use (default is 1). Use max to apply all available CPU threads. For Salmon 8-12 threads are optimal")
    ap.add_argument("-r", "--genome",
                    required=False,
                    choices=genome_list,
                    help="Reference genome")
    ap.add_argument("-a", "--align",
                    required=False,
                    choices=["salmon","hisat2"],
                    default="salmon",
                    help="Choose aligner. Default is Salmon.")
    ap.add_argument("-g", "--go",
                    required=False,
                    action='store_true',
                    help="GO analysis with DAVID")
    args = vars(ap.parse_args())

    ####set thread count for processing
    max_threads=str(multiprocessing.cpu_count())
    threads=args["threads"]
    if threads == "max":
        threads=max_threads

    ###Run FastQC/MultiQC

    def fastqc():
        fastqc_command="fastqc --threads "+threads+" --quiet -o fastqc/ raw-data/*.fastq.gz"
        multiqc_command=["multiqc","-o","fastqc/","fastqc/"]
        print("Running FastQC on raw data")
        subprocess.run(fastqc_command, shell=True)
        print("Running MultiQC")
        subprocess.run(multiqc_command)

    if os.path.isdir(work_dir + "/fastqc") == False:
        os.mkdir(work_dir+"/fastqc")
        fastqc()
    elif len(os.listdir(work_dir + "/fastqc")) == 0:
        fastqc()

    ###trim and align
    align=args["align"]
    if align.lower() == "salmon":
        from alignment import trim,salmon
        trim(threads,work_dir)
        salmon_index=settings.get("salmon_index", {}).get('gencode-V35')
        gtf=settings.get("salmon_gtf", {}).get('gencode-V35')
        fasta=settings.get("FASTA", {}).get('hg19')
        salmon(salmon_index,threads,work_dir,gtf,fasta,script_dir)
    elif align.lower() == "hisat2":
        from alignment import trim,hisat2
        trim(threads,work_dir)
        hisat2()

    ###print total run time
    stop = timeit.default_timer()
    total_time = stop - start
    ty_res = time.gmtime(total_time)
    res = time.strftime("%H:%M:%S",ty_res)
    print('Total run time: ', res)


if __name__ == "__main__":
   main()
