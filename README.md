# deepaccess-package
[![PyPI version](https://badge.fury.io/py/deepaccess.svg)](https://badge.fury.io/py/deepaccess)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/deepaccess/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

This is the code for training and interpretation of an ensemble of convolutional neural networks for multi-task classification. Instructions for downloading and getting started with the current release are available at [https://cgs.csail.mit.edu/deepaccess-package/](https://cgs.csail.mit.edu/deepaccess-package/). deepaccess is available via [pip](https://pypi.org/project/pip/) and [bioconda](https://bioconda.github.io/). The DeepAccess model trained on ATAC-seq data from 10 mouse cell types is available as a [zenodo record](https://zenodo.org/record/4908895#.YL6YpR0pDfY).

## Dependencies
* [bedtools](https://bedtools.readthedocs.io/en/latest/) (v2.29.2)

To run DeepAccess with regions (bedfile format) you must install bedtools and add it to your path. Bedtools binaries are available [here](https://github.com/arq5x/bedtools2/releases).

After installation, you can add bedtools to your path via the terminal or modifying your ~/.bashrc
```
export PATH="/path/to/bedtools:$PATH"
```

## Installation
deepaccess is available on the Python Package Index (PyPI) and can be installed with pip:
```
pip install deepaccess
```
and via bioconda:
```
conda install -c bioconda deepaccess
```

## Training
To train a DeepAccess model for a new task
```
usage: deepaccess train [-h] -l LABELS [LABELS ...]
       		  -out OUT [-ref REFFASTA]
		  [-g GENOME] [-beds BEDFILES [BEDFILES ...]]
		  [-fa FASTA] [-fasta_labels FASTA_LABELS]
                  [-f FRAC_RANDOM] [-nepochs NEPOCHS]
		  [-ho HOLDOUT] [-seed SEED] [-verbose]

optional arguments:
  -h, --help            show this help message and exit
  -l LABELS [LABELS ...], --labels LABELS [LABELS ...]
  -out OUT, --out OUT
  -ref REFFASTA, --refFasta REFFASTA
  -g GENOME, --genome GENOME
                        genome chrom.sizes file
  -beds BEDFILES [BEDFILES ...], --bedfiles BEDFILES [BEDFILES ...]
  -fa FASTA, --fasta FASTA
  -fasta_labels FASTA_LABELS, --fasta_labels FASTA_LABELS
  -f FRAC_RANDOM, --frac_random FRAC_RANDOM
  -nepochs NEPOCHS, --nepochs NEPOCHS
  -ho HOLDOUT, --holdout HOLDOUT
                        chromosome to holdout
  -seed SEED, --seed SEED
  -verbose, --verbose   Print training progress
```
### Arguments
| Argument   | Description | Example |
| ---------  | ----------- | -------- |
| -h, --help | show this help message and exit | NA |
| -l --labels | list of labels for each bed file | C1 C2 C3 |
| -out --out  | output folder name | myoutput |
| -ref --ref  | reference fasta; required with bed input | mm10.fa |
| -g --genome | genome chromosome sizes; required with bed input | default/mm10.chrom.sizes |
| -beds --bedfiles | list of bed files; one of beds or fa input required | C1.bed C2.bed C3.bed |
| -fa --fasta | fasta file;  one of beds or fa input required | C1C2C3.fa |
| -fasta_labels --fasta_labels | text file containing tab delimited labels (0 or 1) for each fasta line with one column for each class | C1C2C3.txt |
| -f  --frac_random | for bed file input fraction of random outgroup regions to add to training | 0.1 |
| -nepochs --nepochs | number of training iterations | 1 |
| -ho --holdout | chromosome name to hold out (only with bed input) | chr19 |
| -verbose --verbose | print training and evaluation progress | NA |
| -seed --seed | set tensorflow seed | 2021 |

## Interpretation
To run interpretation of a DeepAccess model
```
usage: deepaccess interpret [-h] -trainDir TRAINDIR
       		  [-fastas FASTAS [FASTAS ...]]
		  [-l LABELS [LABELS ...]] [
		  -c COMPARISONS [COMPARISONS ...]]
		  [-evalMotifs EVALMOTIFS]
                  [-evalPatterns EVALPATTERNS]
		  [-p POSITION] [-saliency]
		  [-subtract] [-bg BACKGROUND] [-vis]

optional arguments:
  -h, --help            show this help message and exit
  -trainDir TRAINDIR, --trainDir TRAINDIR
  -fastas FASTAS [FASTAS ...], --fastas FASTAS [FASTAS ...]
  -l LABELS [LABELS ...], --labels LABELS [LABELS ...]
  -c COMPARISONS [COMPARISONS ...], --comparisons COMPARISONS [COMPARISONS ...]
  -evalMotifs EVALMOTIFS, --evalMotifs EVALMOTIFS
  -evalPatterns EVALPATTERNS, --evalPatterns EVALPATTERNS
  -p POSITION, --position POSITION
  -saliency, --saliency
  -subtract, --subtract
  -bg BACKGROUND, --background BACKGROUND
  -vis, --makeVis
```
### Arguments 
| Argument   | Description | Example |
| ---------  | ----------- | -------- |
| -h, --help | show this help message and exit | NA |
| -trainDir --trainDir | directory containing trained DeepAccess model | test/ASCL1vsCTCF |
| -fastas --fastas | list of fasta files to evaulate | test/ASCL1vsCTCF/test.fa |
| -l --labels | list of labels for each bed file | C1 C2 C3 |
| -c --comparisons | list of comparisons between different labels | ASCL1vsCTCF ASCL1vsNone runs differential EPE between ASCL1 and CTCF and EPE on ASCL1; C1,C2vsC3 runs differential EPE for (C1 and C2) vs C3 |
| -evalMotifs --evalMotifs | PWM or PCM data base of DNA sequence motifs | default/HMv11_MOUSE.txt |
| -evalPatterns --evalPatterns | fasta file containing DNA sequence patterns | data/ASCL1_space.fa |
| -bg --bg | fasta file containning background sequences | default/backgrounds.fa |
| -saliency --saliency | calculate per base nucleotide importance | NA |
| -subtract --subtract | use subtraction instead of ratio for EPE / DEPE | False |
| -vis --makeVis | to be used with saliency to make plot visualizing results | NA |
