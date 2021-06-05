---
title: Training
layout: default
filename: train
---
## Training

We provide support to train a DeepAccess model with either bed files or fasta files (and labels) as input. For training from bed files, you will also need to download a reference genome and chrom.sizes file, which are available on UCSC:

- mm10 
    - [genome](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz) 
    - [chrom.sizes](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes)
- hg38 
    - [genome](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) 
    - [chrom.sizes](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes)

### Usage
```markdown
usage: deepaccess train [-h] 
            -l LABELS [LABELS ...] -out OUT [-ref REFFASTA]
           [-g GENOME] [-beds BEDFILES [BEDFILES ...]]
           [-fa FASTA] [-fasta_labels FASTA_LABELS]
           [-f FRAC_RANDOM] [-bg BG]
           [-nepochs NEPOCHS] [-model MODEL] [-ho HOLDOUT]
           [-verbose]

optional arguments:
  -h, --help            show this help message and exit
  -l LABELS [LABELS ...], --labels LABELS [LABELS ...]
  -out OUT, --out OUT
  -ref REFFASTA, --refFasta REFFASTA
  -g GENOME, --genome GENOME
                        chrom.sizes file
  -beds BEDFILES [BEDFILES ...], --bedfiles BEDFILES [BEDFILES ...]
  -fa FASTA, --fasta FASTA
  -fasta_labels FASTA_LABELS, --fasta_labels FASTA_LABELS
  -f FRAC_RANDOM, --frac_random FRAC_RANDOM
  -bg BG, --bg BG
  -nepochs NEPOCHS, --nepochs NEPOCHS
  -model MODEL, --model MODEL
  -ho HOLDOUT, --holdout HOLDOUT
                        chromosome to holdout
  -verbose, --verbose   Print training progress
  ```
### Arguments

| Argument   | Description | Example |
| ---------  | ----------- | -------- |
| -h | show this help message and exit | NA |
| -l | list of labels for each bed file | C1 C2 C3 |
| -out | output folder name | myoutput |
| -ref  | reference fasta; required with bed input | mm10.fa |
| -g | genome chromosome sizes; required with bed input | default/mm10.chrom.sizes |
| -beds | list of bed files; one of beds or fa input required | C1.bed C2.bed C3.bed |
| -fa  | fasta file;  one of beds or fa input required | C1C2C3.fa |
| -fasta_labels  | text file containing tab delimited labels (0 or 1) for each fasta line with one column for each class | C1C2C3.txt |
| -f   | for bed file input fraction of random outgroup regions to add to training | 0.1 |
| -bg | fasta file containning background sequences | default/backgrounds.fa |
| -nepochs  | number of training iterations | 1 |
| -model | folder containing base model to begin training | default/DeepAccessMultiMouse |
| -ho | chromosome name to hold out (only with bed input) | chr19 |
| -verbose | print training and evaluation progress | NA |

