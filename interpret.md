---
title: Interpretation
layout: default
filename: interpret
---
## Interpretation
We provide two methods of interpretation of trained DeepAccess models: 
1. ExpectedPatternEffect and DifferentialExpectedPatternEffect
2. Per-nucleotide importance
3. Prediction

#### ExpectedPatternEffect and DifferentialExpectedPatternEffect
ExpectedPatternEffect and DifferentialExpectedPatternEffect can be run using either motifs in PWM or PCM representation (from a database like HOCOMOCO or HOMER) or patterns in a fasta representation which can be used to test spacing or combinations of motifs.

#### Per-nucleotide importance 
Per-nucleotide importance is run using an input of one or more fastas and returns the model-derived importance of each nucleotide within each fasta sequence. 

#### Prediction
For predicting from a trained model, input the directory of the model and one or more fasta sequences:
```markdown
deepaccess interpret -trainDir trained_deepaccessmodel -fastas seqs_1.fa seqs_2.fa
```

The output will be files trained_deepaccessmodel/seqs_1.prediction, trained_deepaccessmodel/seqs_2.prediction, containing the predicted accessibility for each class.


### Usage
```markdown
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