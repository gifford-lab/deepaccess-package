---
title: Quick Start
layout: default
filename: quick-start
---
## Quick Start

These are instructions for how to run the packaged version of DeepAccess training and interpretation. We provide a tutorial for running DeepAccess training and interpretation as a [google colab notebook](https://colab.research.google.com/drive/14q8-qO93-S4SkIwKJaC5WOSJSEQ8OZYI?usp=sharing). 

To run DeepAccess with regions (bedfile format) you must install bedtools and add it to your path. Bedtools binaries are available [here](https://github.com/arq5x/bedtools2/releases).
After installation, you can add bedtools to your path via the terminal or modifying your ~/.bashrc
```markdown
export PATH="/path/to/bedtools:$PATH"
```

DeepAccess can be installed as a command line tool with pip:
```markdown
pip install deepaccess==0.1.0
```

We provide an example of how to run training and interpretation of DeepAccess as a shell script 
```markdown
sh run_ASCL1vsCTCF_DeepAccess.sh 
```
