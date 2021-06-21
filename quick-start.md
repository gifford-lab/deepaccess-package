---
title: Quick Start
layout: default
filename: quick-start
---
## Quick Start

These are instructions for how to run the packaged version of DeepAccess training and interpretation. We provide a tutorial for running DeepAccess training and interpretation as a [google colab notebook](https://colab.research.google.com/drive/14q8-qO93-S4SkIwKJaC5WOSJSEQ8OZYI?usp=sharing). You may also download deepaccess which has been tested for Ubuntu 18.04.3 (will not work on Mac or Windows) as a complete binary using zip or tarball on our [github releases page](https://github.com/gifford-lab/deepaccess-package/releases). 

To run DeepAccess with regions (bedfile format) you must install bedtools and add it to your path. Bedtools binaries are available [here](https://github.com/arq5x/bedtools2/releases).
After installation, you can add bedtools to your path via the terminal or modifying your ~/.bashrc
```markdown
export PATH="/path/to/bedtools:$PATH"
```
If you choose to install deepaccess through PyPI, the default will install a version of tensorflow for GPU usage. If you do not intend to use GPUs, please first install tensorflow for cpu:
```markdown
pip install tensorflow-cpu
```
DeepAccess can be installed as a command line tool with pip:
```markdown
pip install deepaccess
```
or with bioconda:
```markdown
conda install -c bioconda deepaccess 
```
If deepaccess is properly installed, you should be able to run the following without errors:
```markdown
deepaccess -h
```

