mkdir test

python DeepAccessTrainTransfer.py -beds  data/mm10_trimmed_ASCL1.bed data/mm10_trimmed_CTCF.bed \
       -l ASCL1 CTCF -c ASCL1vsCTCF ASCL1vsNone CTCFvsNone  \
       -ref /archive/gl/shared/genomes/mm10/mm10.fa \
       -g /archive/gl/shared/user/jhammelm/old_cluster/jhammelm/genomes/atac_mm10/mm10/mm10.chrom.sizes \
       -out test/ASCL1vsCTCF -nepochs 1 -motifDB HMv11_MOUSE.txt

python DeepAccessInterpret.py -trainDir test/ASCL1vsCTCF \
       -l ASCL1 CTCF -c ASCL1vsCTCF \
       -fastas test/ASCL1vsCTCF_example.fa -saliency -vis

