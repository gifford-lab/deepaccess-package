# DeepAccessTransfer
## Data
The pretrained DeepAccess model used in the paper and executable versions of DeepAccess training and interpretation can be downloaded here.
## Environment setup
Install and activate a conda environment with all necessary Python dependencies by:
```
conda env create -f deepaccessaccess.yml
source activate deepaccessaccess
```
## Training
To train a DeepAccess model for a new task
```
usage: python DeepAccessTrainTransfer.py [-h] -comparisons COMPARISONS [COMPARISONS ...]
                               -l LABELS [LABELS ...] -out OUT [-ref REFFASTA]
                               [-g GENOME] [-beds BEDFILES [BEDFILES ...]]
                               [-fa FASTA] [-fasta_labels FASTA_LABELS]
                               [-f FRAC_RANDOM] [-motifDB MOTIFDB]
                               [-patternFA PATTERNFA] [-bg BG]
                               [-nepochs NEPOCHS] [-model MODEL] [-ho HOLDOUT]
                               [-verbose]

optional arguments:
  -h, --help            show this help message and exit
  -comparisons COMPARISONS [COMPARISONS ...], --comparisons COMPARISONS [COMPARISONS ...]
  -l LABELS [LABELS ...], --labels LABELS [LABELS ...]
  -out OUT, --out OUT
  -ref REFFASTA, --refFasta REFFASTA
  -g GENOME, --genome GENOME
                        genome string or chrom.sizes file
  -beds BEDFILES [BEDFILES ...], --bedfiles BEDFILES [BEDFILES ...]
  -fa FASTA, --fasta FASTA
  -fasta_labels FASTA_LABELS, --fasta_labels FASTA_LABELS
  -f FRAC_RANDOM, --frac_random FRAC_RANDOM
  -motifDB MOTIFDB, --motifDB MOTIFDB
  -patternFA PATTERNFA, --patternFA PATTERNFA
  -bg BG, --bg BG
  -nepochs NEPOCHS, --nepochs NEPOCHS
  -model MODEL, --model MODEL
  -ho HOLDOUT, --holdout HOLDOUT
                        chromosome to holdout
  -verbose, --verbose   Print training progress
```

## Interpretation
To run interpretation of a DeepAccess model
```
usage: python DeepAccessInterpret.py [-h] -trainDir TRAINDIR
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
  -comparisons COMPARISONS [COMPARISONS ...], --comparisons COMPARISONS [COMPARISONS ...]
  -evalMotifs EVALMOTIFS, --evalMotifs EVALMOTIFS
  -evalPatterns EVALPATTERNS, --evalPatterns EVALPATTERNS
  -saliency, --saliency
  -bg BACKGROUND, --background BACKGROUND
  -vis, --makeVis
```