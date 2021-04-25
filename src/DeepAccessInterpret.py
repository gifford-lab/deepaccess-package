#!/bin/env python
import os
#set local keras environment
os.environ["KERAS_HOME"] = ".keras"

import keras
import numpy as np
import  subprocess
import argparse
import os
from CNN import *
from ensemble_utils import *
from DeepAccessTransferModel import *
from ExpectedPatternEffect import *
import pickle
import matplotlib.pyplot as plt


parser=argparse.ArgumentParser()
parser.add_argument('-trainDir','--trainDir',required=True)
parser.add_argument('-fastas','--fastas',nargs="+",required=False,default=[])
parser.add_argument('-l','--labels',nargs="+",required=False,default=[])
parser.add_argument('-c','--comparisons',nargs="+",required=False,default=[])
parser.add_argument('-evalMotifs','--evalMotifs',default=None,required=False)
parser.add_argument('-evalPatterns','--evalPatterns',default=None,required=False)
parser.add_argument('-saliency','--saliency',default=False,action='store_true')
parser.add_argument('-bg','--background',default="default/backgrounds.fa",required=False)
parser.add_argument('-vis','--makeVis',action='store_true',default=False,required=False)
opts=parser.parse_args()

print('-------------------------------------')
print('         Making Predictions          ')
print('-------------------------------------')
DAModel = DeepAccessTransferModel(opts.trainDir)
DAModel.load(opts.trainDir)

for fasta in opts.fastas:
    X=fa_to_onehot(fasta)
    fname = fasta.split('/')[-1].split('.fa')[0]
    DA_pred = DAModel.predict(X)
    if opts.saliency:
        saliencyX = []
        comps = [[comp.strip().split('vs')[0].split('-'),
                  comp.strip().split('vs')[1].split('-')] for comp in opts.comparisons]
        for ci,comp in enumerate(comps):
            c1 = [li for li,l in enumerate(opts.labels) if l in comp[0]]
            c2 = [li for li,l in enumerate(opts.labels) if l in comp[1]]
            saliencyX.append(DAModel.saliency_input(X,c1,c2))
            if opts.makeVis:
                seqs = []
                for seq in open(fasta).read().split('>')[1:]:
                    seqs.append(seq.strip().split('\n')[1].upper())

                for si in range(saliencyX[ci].shape[0]):
                    f=plt.figure(figsize=(13,4))
                    plt.rcParams.update({'font.size': 16})
                    ax=f.add_subplot(1,1,1)
                    ax.bar(range(saliencyX[ci].shape[1]),saliencyX[ci][si,:].reshape((-1,)))
                    for i in range(100):
                        if len(seqs[si]) <= i:
                            break
                        if seqs[si][i] == 'T':
                            ax.get_children()[i].set_color('r')
                        elif seqs[si][i] == 'A':
                            ax.get_children()[i].set_color('g')
                        elif seqs[si][i] == 'C':
                            ax.get_children()[i].set_color('b')
                        elif seqs[si][i] == 'G':
                            ax.get_children()[i].set_color('gold')
                        else:
                            ax.get_children()[i].set_color('k')
                    plt.xticks(range(len(seqs[si])),list(seqs[si]),fontsize=10)
                    plt.ylabel('Saliency')
                    plt.savefig(opts.trainDir+'_'+'-'.join(comp[0])+'vs'+'-'.join(comp[1])+'-saliency'+fname+'_'+str(si)+'.svg')
                    
            np.savetxt(opts.trainDir+'_'+'-'.join(comp[0])+'vs'+'-'.join(comp[1])+'-saliency'+fname+'.txt',saliencyX[ci])
    
    np.savetxt(opts.trainDir+'/'+fname+'.prediction',DA_pred)
    
if opts.evalMotifs != None:
    motifDB = opts.evalMotifs.split('/')[-1].split('.txt')[0]
    print('----------------------------------------')
    print('Performing Differential Motif Evaluation')
    print('----------------------------------------')
    X,X_bg,seqsamples = motif2test(opts.evalMotifs,opts.background)
    
    comps = [(comp.strip().split('vs')[0].split(','),
              comp.strip().split('vs')[1].split(','))
             for comp in opts.comparisons]
    for comp in comps:
        c1 = np.array([li for li,l in enumerate(opts.labels) if l in comp[0]])
        c2 = np.array([li for li,l in enumerate(opts.labels) if l in comp[1]])
        if c1.shape[0] == 0 and c2.shape[0] == 0:
            sys.exit('Error: invalid comparison '+str(comp))
        elif c1.shape[0] == 0:
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
        with open(opts.trainDir+'_EPE_'+'-'.join(comp[0])+'vs'+'-'.join(comp[1])+'-'+motifDB+'.txt','w') as f:
            f.write('\t'.join(EPEdata.keys())+'\n')
            for index in np.argsort(EPEdata[valuecol])[::-1]:
                f.write('\t'.join([str(EPEdata[k][index]) for k in EPEdata.keys()])+'\n')


if opts.evalPatterns != None:
    patternDB = opts.evalPatterns.split('/')[-1].split('.txt')[0]
    print('----------------------------------------')
    print('Performing Differential Motif Evaluation')
    print('----------------------------------------')
    X,X_bg,seqsamples = fasta2test(opts.evalPatterns,opts.background)
    
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
        with open(opts.trainDir+'_EPE_'+'-'.join(comp[0])+'vs'+'-'.join(comp[1])+'-'+patternDB+'.txt','w') as f:
            f.write('\t'.join(EPEdata.keys())+'\n')
            for index in np.argsort(EPEdata[valuecol])[::-1]:
                f.write('\t'.join([str(EPEdata[k][index]) for k in EPEdata.keys()])+'\n')

