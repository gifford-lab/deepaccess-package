#!/bin/env python
import keras
import numpy as np
import  subprocess
import argparse
import os
from CNN import *
from ensemble_utils import *
from importance_utils import saliency
import pickle
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser()
parser.add_argument('-trainDir','--trainDir',required=True)
parser.add_argument('-fastas','--fastas',nargs="+",required=False,default=[])
parser.add_argument('-l','--labels',nargs="+",required=False,default=[])
parser.add_argument('-comparisons','--comparisons',nargs="+",required=False,default=[])
parser.add_argument('-evalMotifs','--evalMotifs',default=None,required=False)
parser.add_argument('-saliency','--saliency',default=False,action='store_true')
parser.add_argument('-bg','--background',default="backgrounds.fa",required=False)
parser.add_argument('-vis','--makeVis',action='store_true',default=False,required=False)
opts=parser.parse_args()

print('-------------------------------------')
print('         Making Predictions          ')
print('-------------------------------------')
#TODO make work for full models
model_folders = [opts.trainDir+"/"+d for d in os.listdir(opts.trainDir) if os.path.isdir(opts.trainDir+"/"+d)]

with open(opts.trainDir+"/model_acc.pkl","rb") as f:
    initial_accuracies = pickle.load(f)

prevkeys  = list(initial_accuracies.keys())
accuracies = {}

for key in prevkeys:
    accuracies[key.split('/')[-1]] = initial_accuracies[key]
Y = np.loadtxt(opts.trainDir+'/test.txt')
for fasta in opts.fastas:
    fname = fasta.split('.fa')[0].split('/')[-1]
    X = fa_to_onehot(fasta)
    pred_mat = np.zeros((X.shape[0],Y.shape[1]))
    saliencyX = [np.zeros((X.shape[0],X.shape[1])) for _ in opts.comparisons]
    for mi,model in enumerate(model_folders):
        #load keras model
        cnn = keras.models.load_model(model+'/model.h5')

        predictions = cnn.predict(X)
        pred_mat += predictions*accuracies[model.split('/')[-1]]
        if opts.saliency:
            comps = [[comp.strip().split('vs')[0].split('-'),
                      comp.strip().split('vs')[1].split('-')] for comp in opts.comparisons]
            for ci,comp in enumerate(comps):
                c1 = [li for li,l in enumerate(opts.labels) if l in comp[0]]
                c2 = [li for li,l in enumerate(opts.labels) if l in comp[1]]
                saliencyX[ci] += np.sum(saliency(0,model+'/model.h5',X,c1,c2,n=5,batch_size=528)*X,axis=2)*accuracies[model.split('/')[-1]]
        del cnn
    if opts.saliency:
        for ci,comp in enumerate(comps):
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
    pred_mat = pred_mat/sum(accuracies.values())
    np.savetxt(opts.trainDir+'/'+fname+'.prediction',pred_mat)
    
if opts.evalMotifs != None:
    print('----------------------------------------')
    print('Performing Differential Motif Evaluation')
    print('----------------------------------------')
    X,X_bg,seqsamples = motif2test(opts.evalMotifs,opts.background)
    
    pred_mat = np.zeros((X.shape[0],len(opts.labels)))
    pred_null_mat = np.zeros((X_bg.shape[0],len(opts.labels)))
    for mi,model in enumerate(model_folders):
        model_folder = opts.trainDir + "/" + model.split('/')[-1]
        trained_cnn = keras.models.load_model(model_folder+"/model.h5")
        pred_mat += accuracies[model.split('/')[-1]]*trained_cnn.predict(X)
        pred_null_mat += accuracies[model.split('/')[-1]]*trained_cnn.predict(X_bg)
        del trained_cnn
    pred_mat = pred_mat/sum(accuracies.values())
    pred_null_mat = pred_null_mat/sum(accuracies.values())

    motifs_seqs = np.array([s.split('_bg')[0] for s in seqsamples])
    scores = np.zeros((len(motifs_seqs),pred_mat.shape[1]))
    motifs = list(set([s.split('_bg')[0] for s in seqsamples]))
    comps = [[comp.strip().split('vs')[0].split('-'),
          comp.strip().split('vs')[1].split('-')] for comp in opts.comparisons]
    for comp in comps:
        c1 = np.array([li for li,l in enumerate(opts.labels) if l in comp[0]])
        c2 = np.array([li for li,l in enumerate(opts.labels) if l in comp[1]])
        with open(opts.trainDir+'_'+'-'.join(comp[0])+'vs'+'-'.join(comp[1])+'-'+opts.evalMotifs.split('/')[-1].split('.txt')[0]+'-subtraction.txt','w') as f:
            f.write('TF\tL2FC\tWilcoxonES\t-log10(p-value)\t-log10(adj p-value)\n')
            for i,motif in enumerate(motifs):
                ind_motif_seqs = np.where(motif == motifs_seqs)[0]
                #for a given motif:
                # fold change relative to background
                # compared across two cell contexts
                if c1.shape[0] == 1:
                    num = pred_mat[ind_motif_seqs,c1]
                    denom = pred_null_mat[:,c1].reshape((-1,))
                    fc1 = num-denom
                elif c1.shape[0] > 1:
                    num = np.mean([pred_mat[ind_motif_seqs,ind] for ind in c1],axis=0).reshape((-1,))
                    denom = np.mean([pred_null_mat[:,ind] for ind in c1],axis=0).reshape((-1,))
                    fc1 = num-denom
                if c2.shape[0] == 0:
                    stat,pval = wilcoxon(num,denom)
                    adjpval = pval*len(motifs)
                    f.write(motif+'\t'+str(np.mean(fc1))+'\t'+str(stat)+'\t'+str(-np.log10(pval))+'\t'+str(-np.log10(adjpval))+'\n')     
                if c2.shape[0] == 1:
                    num = pred_mat[ind_motif_seqs,c2]
                    denom = pred_null_mat[:,c2].reshape((-1,))
                    fc2 = num-denom
                elif c2.shape[0] > 1:
                    num = np.mean([pred_mat[ind_motif_seqs,ind] for ind in c2],axis=0).reshape((-1,))
                    denom = np.mean([pred_null_mat[:,ind] for ind in c2],axis=0).reshape((-1,))            
                    fc2 = num-denom
                if c1.shape[0] == 0:
                    stat,pval = wilcoxon(num,denom)
                    adjpval = pval*len(motifs)
                    f.write(motif+'\t'+str(np.mean(fc2))+'\t'+str(stat)+'\t'+str(-np.log10(pval))+'\t'+str(-np.log10(adjpval))+'\n')            
                if c1.shape[0] != 0 and c2.shape[0] != 0:
                    stat,pval = wilcoxon(fc1,fc2)
                    adjpval = pval*len(motifs)
                    f.write(motif+'\t'+str(np.mean(fc1-fc2))+'\t'+str(stat)+'\t'+str(-np.log10(pval))+'\t'+str(-np.log10(adjpval))+'\n')
