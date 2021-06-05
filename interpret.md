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
usage: deepaccess interpret [-h] -trainDir TRAINDIR
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
