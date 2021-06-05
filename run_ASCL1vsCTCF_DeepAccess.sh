mkdir test

#pip install deepaccess

wget https://raw.githubusercontent.com/jhammelman/deepaccess-package/main/data/mm10_trimmed_CTCF.bed
wget https://raw.githubusercontent.com/jhammelman/deepaccess-package/main/data/mm10_trimmed_ASCL1.bed

wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz

wget https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
pip install bedtools

#example train from bed file
deepaccess train \
	   -beds mm10_trimmed_ASCL1.bed mm10_trimmed_CTCF.bed \
	   -l ASCL1 CTCF \
	   -ref mm10.fa -g mm10.chrom.sizes \
	   -out test/ASCL1vsCTCF -nepochs 1

wget https://raw.githubusercontent.com/jhammelman/deepaccess-package/main/default/HMv11_MOUSE.txt

#example test motif activity
deepaccess interpret -trainDir test/ASCL1vsCTCF  \
	   -l ASCL1 CTCF \
	   -c ASCL1vsCTCF ASCL1vsNone \
	   -evalMotifs HMv11_MOUSE.txt

wget https://raw.githubusercontent.com/jhammelman/deepaccess-package/main/data/ASCL1_space.fa

#example run test motif spacing on differential and ascl1 activity
deepaccess interpret -trainDir test/ASCL1vsCTCF \
	   -l ASCL1 CTCF \
	   -c ASCL1vsCTCF ASCL1vsNone \
	   -evalPatterns ASCL1_space.fa

wget https://raw.githubusercontent.com/jhammelman/deepaccess-package/main/data/ASCL1vsCTCF_example.fa

#example run saliency visualization
deepaccess interpret -trainDir test/ASCL1vsCTCF \
       -l ASCL1 CTCF -c ASCL1vsCTCF ASCL1vsNone CTCFvsNone \
       -fastas ASCL1vsCTCF_example.fa -saliency -vis
       
wget https://raw.githubusercontent.com/jhammelman/deepaccess-package/main/data/ASCL1vsCTCF_train.fa
wget https://raw.githubusercontent.com/jhammelman/deepaccess-package/main/data/ASCL1vsCTCF_train.txt

#example train from fasta
deepaccess train --fasta data/ASCL1vsCTCF_train.fa  \
       -fasta_labels data/ASCL1vsCTCF_train.txt \
       -l ASCL1 CTCF \
       -ref mm10.fa \
       -g mm10.chrom.sizes  \
       -out test/ASCL1vsCTCF_from_fasta -nepochs 1 

