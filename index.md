## Quick Start

These are instructions for how to run the packaged version of DeepAccess training and interpretation. If you would like to be able to make modifications to the code or run DeepAccess on GPU, you will want to use instructions for how to run DeepAccess from [source](https://github.com/gifford-lab/DeepAccessTransfer) code.  

To run DeepAccess with regions (bedfile format) you must install bedtools and add it to your path. Bedtools binaries are available [here](https://github.com/arq5x/bedtools2/releases).
After installation, you can add bedtools to your path via the terminal or modifying your ~/.bashrc
```markdown
export PATH="/path/to/bedtools:$PATH"
```

After installing bedtools, you can install DeepAccessTransfer by downloading the tarball containing the executable with supporting data and files.

We provide an example of how to run training and interpretation of DeepAccess in the shell script 
```markdown
sh run_ASCL1vsCTCF_DeepAccessExec.sh 
```
## Citation
[Discovering differential genome sequence activity with interpretable and efficient deep learning](https://www.biorxiv.org/content/10.1101/2021.02.26.433073v1.abstract?%3Fcollection=) <br />  by Jennifer Hammelman and David K. Gifford

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
usage: dist/DeepAccessTrainTransfer [-h] 
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

## Interpretation
We provide two methods of interpretation of trained DeepAccess models: 
1. ExpectedPatternEffect and DifferentialExpectedPatternEffect
2. Per-nucleotide importance

#### ExpectedPatternEffect and DifferentialExpectedPatternEffect
ExpectedPatternEffect and DifferentialExpectedPatternEffect can be run using either motifs in PWM or PCM representation (from a database like HOCOMOCO or HOMER) or patterns in a fasta representation which can be used to test spacing or combinations of motifs.

#### Per-nucleotide importance 
Per-nucleotide importance is run using an input of one or more fastas and returns the model-derived importance of each nucleotide within each fasta sequence. 

### Usage
```markdown
usage: dist/DeepAccessInterpret [-h] -trainDir TRAINDIR
                [-fastas FASTAS [FASTAS ...]]
                [-l LABELS [LABELS ...]]
                [-comparisons COMPARISONS [COMPARISONS ...]]
                [-evalMotifs EVALMOTIFS]
                [-evalPatterns EVALPATTERNS] [-saliency]
                [-bg BACKGROUND] [-vis]

optional arguments:
  -h, --help            show this help message and exit
  -trainDir TRAINDIR, --trainDir TRAINDIR
  -fastas FASTAS [FASTAS ...], --fastas FASTAS [FASTAS ...]
  -l LABELS [LABELS ...], --labels LABELS [LABELS ...]
  -c COMPARISONS [COMPARISONS ...], --comparisons COMPARISONS [COMPARISONS ...]
  -evalMotifs EVALMOTIFS, --evalMotifs EVALMOTIFS
  -evalPatterns EVALPATTERNS, --evalPatterns EVALPATTERNS
  -saliency, --saliency
  -bg BACKGROUND, --background BACKGROUND
  -vis, --makeVis
```
### Arguments 

| Argument   | Description | Example |
| ---------  | ----------- | -------- |
| -h  | show this help message and exit | NA |
| -trainDir  | directory containing trained DeepAccess model | test/ASCL1vsCTCF |
| -fastas | list of fasta files to evaulate | test/ASCL1vsCTCF/test.fa |
| -l | list of labels for each bed file | C1 C2 C3 |
| -c  | list of comparisons between different labels | ASCL1vsCTCF ASCL1vsNone runs differential EPE between ASCL1 and CTCF and EPE on ASCL1; C1,C2vsC3 runs differential EPE for (C1 and C2) vs C3 |
| -evalMotifs  | PWM or PCM data base of DNA sequence motifs | default/HMv11_MOUSE.txt |
| -evalPatterns  | fasta file containing DNA sequence patterns | data/ASCL1_space.fa |
| -bg | fasta file containning background sequences | default/backgrounds.fa |
| -saliency | calculate per base nucleotide importance | NA |
| -vis | to be used with saliency to make plot visualizing results | NA |
