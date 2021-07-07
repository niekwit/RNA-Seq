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
#import glob

def main(script_dir):
    
    work_dir=os.getcwd()

    ###loads available RNA-Seq settings
    if os.path.exists(script_dir+"/settings.yaml") == True:
        with open(script_dir+"/settings.yaml") as file:
            settings=yaml.full_load(file)
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
    ap.add_argument("-r", "--reference",
                    required=False,
                    choices=genome_list,
                    help="Reference genome")
    ap.add_argument("-s", "--species",
                    required=True,
                    choices=["mouse","human"],
                    help="Set species.")
    ap.add_argument("-a", "--align",
                    required=False,
                    choices=["salmon","hisat2"],
                    default="salmon",
                    help="Choose aligner. Default is Salmon.")
    ap.add_argument("--go",
                    required=False,
                    action='store_true',
                    help="Gene set enrichment analysis with Enrichr")
    ap.add_argument("--skip-fastqc",
                    required=False,
                    action='store_true',
                    default=False,
                    help="Skip FastQC/MultiQC")
    args = vars(ap.parse_args())

    ####set thread count for processing
    max_threads=str(multiprocessing.cpu_count())
    threads=args["threads"]
    if threads == "max":
        threads=max_threads

    ###Run FastQC/MultiQC
    skip_fastqc=args["skip_fastqc"]
    if not skip_fastqc:
        utils.fastqc(work_dir,threads,file_extension,exe_dict)
    else:
        print("Skipping FastQC/MultiQC analysis")

    ###Set species variable
    species=args["species"]

    ###trim and align
    align=args["align"]
    if align.lower() == "salmon":
        utils.trim(threads,work_dir)
        salmon_index=settings["salmon_index"]["gencode-v35"]
        gtf=settings["salmon_gtf"]["gencode-v35"]
        fasta=settings["FASTA"]["gencode-v35"]
        utils.salmon(salmon_index,str(threads),work_dir,gtf,fasta,script_dir,settings)
        utils.diff_expr(work_dir,gtf,script_dir,species)
    elif align.lower() == "hisat2":
        from alignment import trim,hisat2
        utils.trim(threads,work_dir)
        #hisat2()


if __name__ == "__main__":
    #start run timer
    start = timeit.default_timer()
    
    script_dir=os.path.abspath(os.path.dirname(__file__))
    sys.path.append(script_dir)#adds script directory to runtime (for importing modules)
    import utils as utils

    #utils.install_packages()
    main(script_dir)

    #print total run time
    stop = timeit.default_timer()
    total_time = stop - start
    ty_res = time.gmtime(total_time)
    res = time.strftime("%H:%M:%S",ty_res)
    print('Total run time: ', res)
