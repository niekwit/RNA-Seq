# RNA-Seq pipeline

## Software dependencies:
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)
- [Python](https://www.python.org/)
	- [PyYAML](https://python.land/data-processing/python-yaml)
- [R](https://www.r-project.org/)
	- [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)
	- [readr](https://www.rdocumentation.org/packages/readr/versions/1.3.1)
	- [GenomicFeatures](https://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
	- [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
	- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)



## Installation:
In the command line, navigate to a directory where you want to install the scripts and run:
> `git clone https://github.com/niekwit/RNA-Seq.git`

`rna-seq.py` can be permamently added to $PATH by adding the following line to `~/.bashrc`:
> `export PATH=/path/to/RNA-Seq:$PATH`


## Configuration:

## Usage:

```
usage: rna-seq.py [-h] [-t <int>] [-r {gencode-v35}] [-a {salmon,hisat2}] [-g]

optional arguments:
  -h, --help            show this help message and exit
  -t <int>, --threads <int>
                        Number of CPU threads to use (default is 1). Use max
                        to apply all available CPU threads. For Salmon 8-12
                        threads are optimal
  -r {gencode-v35}, --reference {gencode-v35}
                        Reference genome
  -a {salmon,hisat2}, --align {salmon,hisat2}
                        Choose aligner. Default is Salmon.
  -g, --go              GO analysis with DAVID
```

## Output:
