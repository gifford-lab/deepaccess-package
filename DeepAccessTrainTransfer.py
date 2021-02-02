#!/bin/env python
import os
#set local keras environment
os.environ["KERAS_HOME"] = ".keras"

import keras
import numpy as np
import subprocess
import argparse
from CNN import *
from DeepAccessTransferModel import *
from ensemble_utils import *
from ExpectedPatternEffect import *
import pickle
from sklearn.metrics import average_precision_score

#decided by trained network size
WINDOWSIZE=100

def create_train_test(bedfile,classes,out,holdout='chr19'):
    #merge categories code from merge2cat
    #plus hold out chr19 for testing
    train_classifications = []
    test_classifications = []
    train_beds = []
    test_beds = []

    for line in open(bedfile):
        data = line.strip().split()
        region_len = int(data[2]) - int(data[1])
        if region_len < WINDOWSIZE:
            data[2] = str(int(data[1]) + (WINDOWSIZE-region_len))
        if 'random' in data[0] or '_' in data[0]:
            continue
        if data[0] == holdout:
            test_classifications.append(np.zeros((len(classes,))))
            if data[3] != ".":
                for c in data[3].split(","):
                    test_classifications[-1][classes.index(c)] = 1
            with open(out+'/test.bed','a') as f:
                f.write('\t'.join(data[:3])+'\n')
        else:
            train_classifications.append(np.zeros((len(classes,))))
            if data[3] != ".":
                for c in data[3].split(","):
                    train_classifications[-1][classes.index(c)] = 1
            with open(out+'/train.bed','a') as f:
                f.write('\t'.join(data[:3])+'\n')
    np.savetxt(out+'/test.txt',np.array(test_classifications))
    np.savetxt(out+'/train.txt',np.array(train_classifications))
    

parser=argparse.ArgumentParser()
parser.add_argument('-comparisons','--comparisons',nargs="+",required=True)
parser.add_argument('-l','--labels',nargs="+",required=True)
parser.add_argument('-out','--out',required=True)
parser.add_argument('-ref','--refFasta',required=False,default=None)
parser.add_argument('-g','--genome',required=False,default=None,
                    help='genome string or chrom.sizes file')
parser.add_argument('-beds','--bedfiles',nargs='+',required=False)
parser.add_argument('-fa','--fasta',default=None,required=False)
parser.add_argument('-fasta_labels','--fasta_labels',default=None,
                    required=False)
parser.add_argument('-f','--frac_random',required=False,default=0.1,
                    type=float)
parser.add_argument('-motifDB','--motifDB',default="HMv11_MOUSE.txt",
                    required=False)
parser.add_argument('-patternFA','--patternFA',default=None,required=False)
parser.add_argument('-bg','--bg',default="backgrounds.fa",required=False)
parser.add_argument('-nepochs','--nepochs',default=5,type=int,required=False)
parser.add_argument('-model','--model',default='DeepAccessMultiMouse',
                    required=False)
parser.add_argument('-ho','--holdout',default='chr19',required=False,
                    help='chromosome to holdout')
parser.add_argument('-verbose','--verbose',action='store_true',
                    default=False,required=False,
                    help='Print training progress')
opts=parser.parse_args()

ensure_dir(opts.out,exit_if_exists=True)
        
if opts.fasta == None:
    assert(opts.refFasta != None)
    assert(opts.genome != None)
    assert(len(opts.labels) == len(opts.bedfiles))
    total_peaks = 0
    with open(opts.out+'/all_peaks.bed','a') as f:
        for i,bedfile in enumerate(opts.bedfiles):
            for line in open(bedfile):
                f.write("\t".join(line.strip().split('\t')[:3])+"\t"+opts.labels[i]+'\n')
                total_peaks += 1
        
    print('Combining peak files')
    subprocess.call(["sort -k1,1 -k2,2n  "+opts.out+"/all_peaks.bed | bedtools merge > "+opts.out+"/sloppeaks.bed"],shell=True)
    subprocess.call(["bedtools makewindows -w "+str(WINDOWSIZE)+" -b "+opts.out+'/sloppeaks.bed | sort -k1,1 -k2,2n | uniq > '+opts.out+"/windows.bed"],shell=True)
    subprocess.call(["bedtools random -l "+str(WINDOWSIZE)+" -n "+str(int(total_peaks*opts.frac_random))+" -g "+opts.genome+" > "+opts.out+"/randompeaks.bed"],shell=True)
    subprocess.call(["cat "+opts.out+"/randompeaks.bed "+opts.out+"/windows.bed | sort -k1,1 -k2,2n | cut -f1,2,3 > "+opts.out+"/windows_and_random.bed"],shell=True)
    subprocess.call(["bedtools intersect -loj -a "+opts.out+'/windows_and_random.bed -b '+opts.out+'/all_peaks.bed > '+opts.out+'/windows_overlap_peaks.bed'],shell=True)
    subprocess.call(["bedtools groupby -grp 1,2,3 -o distinct -c 7 -i "+opts.out+'/windows_overlap_peaks.bed > '+opts.out+'/windows_grouped.bed'],shell=True)
    create_train_test(opts.out+'/windows_grouped.bed',opts.labels,opts.out,opts.holdout)
    subprocess.call(['bedtools getfasta -bed '+opts.out+'/train.bed -fi '+opts.refFasta+' -fo '+opts.out+'/train.fa'],shell=True)
    subprocess.call(['bedtools getfasta -bed '+opts.out+'/test.bed -fi '+opts.refFasta+' -fo '+opts.out+'/test.fa'],shell=True)
    trainX = fa_to_onehot(opts.out+'/train.fa')

    
    trainY = np.loadtxt(opts.out+'/train.txt')
    testY = np.loadtxt(opts.out+'/test.txt')
    
else:
    fasta_all = [l for l in open(opts.fasta).read().split('>')[1:]]
    label_all = open(opts.fasta_labels).readlines()
    assert(len(fasta_all) == len(label_all))
    assert(len(label_all[0].strip().split('\t')) == len(opts.labels))
    random_index = np.random.permutation(len(fasta_all))
    test_index = random_index[:int(len(fasta_all)*0.1)]
    train_index = random_index[int(len(fasta_all)*0.1):]
    with open(opts.out+'/train.fa','w') as f:
        for i in train_index:
            f.write('>'+fasta_all[i])
    with open(opts.out+'/train.txt','w') as f:
        for i in train_index:
            f.write(label_all[i])
    with open(opts.out+'/test.fa','w') as f:
        for i in test_index:
            f.write('>'+fasta_all[i])
    with open(opts.out+'/test.txt','w') as f:
        for i in test_index:
            f.write(label_all[i])
    trainX = fa_to_onehot(opts.out+'/train.fa')
    trainY = np.loadtxt(opts.out+'/train.txt')
    testY = np.loadtxt(opts.out+'/test.txt')
    
sample_weights = np.ones((trainX.shape[0],))
testX = fa_to_onehot(opts.out+'/test.fa')


#load keras model
print('-------------------------------------')
print('             Training Model          ')
print('-------------------------------------')
DAModel = DeepAccessTransferModel(opts.out)
DAModel.load(opts.model)
if opts.verbose:
    with_verbose=1
else:
    with_verbose=0

DAModel.retrain(trainX,trainY,
              n_epochs=opts.nepochs,
              verbose=with_verbose)

print('-------------------------------------')
print('           Evaluating Model          ')
print('-------------------------------------')
train_pred = DAModel.predict(trainX)
test_pred = DAModel.predict(testX)
np.savetxt(opts.out+'/train_predictions.txt',train_pred)
np.savetxt(opts.out+'/test_predictions.txt',test_pred)

metric = average_precision_score

with open(opts.out+'/performance.txt','w') as f:
    f.write('Train performance: '+str(metric(trainY,train_pred))+'\n')
with open(opts.out+'/performance.txt','a') as f:
    f.write('Test performance: '+str(metric(testY,test_pred))+'\n')
    
print('----------------------------------------')
print('Performing Differential Motif Evaluation')
print('----------------------------------------')
X,X_bg,seqsamples = motif2test(opts.motifDB,opts.bg)
comps = [(comp.strip().split('vs')[0].split(','),
          comp.strip().split('vs')[1].split(','))
          for comp in opts.comparisons]
for comp in comps:
    c1 = np.array([li for li,l in enumerate(opts.labels) if l in comp[0]])
    c2 = np.array([li for li,l in enumerate(opts.labels) if l in comp[1]])
    if c1.shape[0] == 0:
        _,_,EPEdata = ExpectedPatternEffect(DAModel.predict,
                                            c2,X,X_bg,seqsamples)
        valuecol = 'ExpectedPatternEffect'
    elif c2.shape[0] == 0:
        _,_,EPEdata = ExpectedPatternEffect(DAModel.predict,
                                            c1,X,X_bg,seqsamples)
        valuecol = 'ExpectedPatternEffect'
    else:
        _,_,EPEdata = DifferentialExpectedPatternEffect(DAModel.predict,
                                                        c1,c2,X,X_bg,seqsamples)
        valuecol = 'DifferentialExpectedPatternEffect'
    with open(opts.out+'_EPE_'+'-'.join(comp[0])+'vs'+'-'.join(comp[1])+'.txt','w') as f:
        f.write('\t'.join(EPEdata.keys())+'\n')
        for index in np.argsort(EPEdata[valuecol])[::-1]:
            f.write('\t'.join([str(EPEdata[k][index]) for k in EPEdata.keys()])+'\n')
