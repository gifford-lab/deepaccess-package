#!/bin/env python
import os
os.environ["KERAS_HOME"] = ".keras"
import keras
import numpy as np
import  subprocess
import argparse
from CNN import *
from ensemble_utils import *
from scipy.stats import wilcoxon
import pickle
from sklearn.metrics import average_precision_score
from sklearn.metrics.pairwise import cosine_similarity

def retrain(model,X,y,sample_weights,n_epochs,verbose=0):
    callbacks = [keras.callbacks.EarlyStopping(monitor='val_loss', 
                                                   patience=3),
                 keras.callbacks.History()]
    history = model.fit(x=X,
                             y=y,epochs=n_epochs,
                             shuffle=True,
                             validation_split=0.2,
                             batch_size=250,verbose=verbose,
                             callbacks=callbacks)
    return history,model



def ensure_dir(fdir):
    if not os.path.exists(fdir):
        os.makedirs(fdir)
    
def create_train_test(bedfile,classes,out,holdout='chr19'):
    #merge categories code from merge2cat
    #plus hold out chr19 for testing
    train_classifications = []
    test_classifications = []
    train_beds = []
    test_beds = []

    for line in open(bedfile):
        data = line.strip().split()
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
    
#decided by trained network size
WINDOWSIZE=100
#TODO make so if not mm10 or hg38, user can make
# file for additional genomes
MOUSEBEDFILE="mm10.100.sorted.bed"
HUMANBEDFILE="hg38.100.sorted.bed"

parser=argparse.ArgumentParser()
parser.add_argument('-comparisons','--comparisons',nargs="+",required=True)
parser.add_argument('-l','--labels',nargs="+",required=True)
parser.add_argument('-out','--out',required=True)
parser.add_argument('-ref','--refFasta',required=False,default=None)
parser.add_argument('-g','--genome',required=False,default=None,help='genome string or chrom.sizes file')
parser.add_argument('-beds','--bedfiles',nargs='+',required=False)
parser.add_argument('-fa','--fasta',default=None,required=False)
parser.add_argument('-fasta_labels','--fasta_labels',default=None,required=False)
parser.add_argument('-f','--frac_random',required=False,default=0.1,type=float)
parser.add_argument('-bams','--bamfiles',nargs='+',default=[],required=False)
parser.add_argument('-motifDB','--motifDB',default="HMv11_consensus.txt",required=False)
parser.add_argument('-bg','--bg',default="backgrounds.fa",required=False)
parser.add_argument('-nepochs','--nepochs',default=5,type=int,required=False)
parser.add_argument('-model','--model',default='DeepAccessMultiMouse',required=False)
parser.add_argument('-multiclass','--multiclass',default=False,action='store_true',required=False,help='set to multiclass when only 1 positive label exists for region')
parser.add_argument('-ho','--holdout',default='chr19',required=False)
parser.add_argument('-verbose','--verbose',action='store_true',default=False,required=False)
opts=parser.parse_args()

model_folders = sorted([opts.model+"/"+d for d in os.listdir(opts.model) if os.path.isdir(opts.model+"/"+d)],key = lambda x: int(x.split('_')[-1]))
ensure_dir(opts.out)
if opts.fasta == None:
    assert(opts.refFasta != None)
    assert(opts.genome != None)
    assert(len(opts.labels) == len(opts.bedfiles))
    if len(opts.bamfiles) != 0:
        assert(len(opts.bedfiles) == len(opts.bamfiles))


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
    if opts.multiclass:
        with open(opts.out+'/all_peaks.bed','a') as f:
            for i,bedfile in enumerate(opts.bedfiles):
                for line in open(opts.out+'/randompeaks.bed'):
                    f.write("\t".join(line.strip().split('\t')[:3])+"\trandom\n")
                    total_peaks += 1
    subprocess.call(["bedtools intersect -loj -a "+opts.out+'/windows_and_random.bed -b '+opts.out+'/all_peaks.bed > '+opts.out+'/windows_overlap_peaks.bed'],shell=True)
    subprocess.call(["bedtools groupby -grp 1,2,3 -o distinct -c 7 -i "+opts.out+'/windows_overlap_peaks.bed > '+opts.out+'/windows_grouped.bed'],shell=True)
    if  opts.multiclass:
        create_train_test(opts.out+'/windows_grouped.bed',opts.labels+['random'],opts.out,opts.holdout)
    else:
        create_train_test(opts.out+'/windows_grouped.bed',opts.labels,opts.out,opts.holdout)



    if len(opts.bamfiles) > 0:
        #get read counts over bed train
        subprocess.call(['bedtools multicov -bams '+' '.join(opts.bamfiles)+' -bed '+opts.out+'/train.bed > '+opts.out+'/train.regress'],shell=True)
        subprocess.call(['bedtools multicov -bams '+' '.join(opts.bamfiles)+' -bed '+opts.out+'/test.bed > '+opts.out+'/test.regress'],shell=True)


    subprocess.call(['bedtools getfasta -bed '+opts.out+'/train.bed -fi '+opts.refFasta+' -fo '+opts.out+'/train.fa'],shell=True)
    subprocess.call(['bedtools getfasta -bed '+opts.out+'/test.bed -fi '+opts.refFasta+' -fo '+opts.out+'/test.fa'],shell=True)
    trainX = fa_to_onehot(opts.out+'/train.fa')

    if len(opts.bamfiles) == 0:
        trainY = np.loadtxt(opts.out+'/train.txt')
        testY = np.loadtxt(opts.out+'/test.txt')
    else:
        trainY = []
        for line in open(opts.out+'/train.regress'):
            trainY.append(np.array([float(v) for v in line.strip().split('\t')[2:]]))
        trainY = np.array(trainY)
        testY = []
        for line in open(opts.out+'/test.regress'):
            testY.append(np.array([float(v) for v in line.strip().split('\t')[2:]]))
        testY = np.array(testY)
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
accuracies = {}
#load keras model

print('-------------------------------------')
print('             Training Model          ')
print('-------------------------------------')
for mi,model in enumerate(model_folders):
    cnn = keras.models.load_model(model+'/model.h5')
    #add up to last layer
    new_cnn = Sequential()
    for layer in cnn.layers[:-1]:
        new_cnn.add(layer)
    #add new layer with new output
    if len(opts.bamfiles) > 0 or np.any(np.logical_and(testY != 0,testY != 1)):
        new_cnn.add(Dense(trainY.shape[1],activation='softplus'))
        adam = optimizers.Adam(lr=1e-5,clipnorm=0.1,decay=(1e-5/100.0))
        new_cnn.compile(optimizer=adam,loss='poisson',metrics=['cosine_similarity'])
    elif opts.multiclass:
        new_cnn.add(Dense(trainY.shape[1],activation='softmax'))
        adam = optimizers.Adam(lr=1e-4,clipnorm=0.5,decay=(1e-4/100.0))
        new_cnn.compile(optimizer=adam,loss='categorical_crossentropy',metrics=['accuracy'])
    else:
        new_cnn.add(Dense(trainY.shape[1],activation='sigmoid'))
        adam = optimizers.Adam(lr=1e-4,clipnorm=0.5,decay=(1e-4/100.0))
        new_cnn.compile(optimizer=adam,loss='binary_crossentropy',metrics=['accuracy'])

    model_folder = opts.out + "/" + model.split('/')[-1]
    ensure_dir(model_folder)
    if opts.verbose:
        with_verbose=1
    else:
        with_verbose=0
    history,trained_cnn = retrain(new_cnn,trainX,trainY,sample_weights,n_epochs=opts.nepochs,verbose=with_verbose)
    loss,train_acc = trained_cnn.evaluate(trainX,trainY,batch_size=250)
    accuracies[model_folder] = train_acc
    trained_cnn.save(model_folder+"/model.h5")
    np.save(model_folder+'/history.npy',history.history)
    np.save(model_folder+'/sample_weights.npy',sample_weights)
    with open(model_folder+'/model_summary.txt','w') as f:
        f.write(trained_cnn.to_yaml())
        f.write("\nTraining Accuracy: "+str(train_acc)+"\n")        

    sample_weights += np.linalg.norm(trainY-trained_cnn.predict(trainX,batch_size=250))
    del cnn
    del new_cnn
    del trained_cnn

with open(opts.out+'/model_acc.pkl','wb') as f:
    pickle.dump(accuracies,f)


print('-------------------------------------')
print('           Evaluating Model          ')
print('-------------------------------------')
train_pred = np.zeros(trainY.shape)
test_pred = np.zeros(testY.shape)
for mi,model in enumerate(model_folders):
    model_folder = opts.out + "/" + model.split('/')[-1]
    trained_cnn = keras.models.load_model(model_folder+"/model.h5")
    test_pred += accuracies[model_folder]*trained_cnn.predict(testX,batch_size=250)
    train_pred += accuracies[model_folder]*trained_cnn.predict(trainX,batch_size=250)
    del trained_cnn
train_pred = train_pred/sum(accuracies.values())
test_pred = test_pred/sum(accuracies.values())

if len(opts.bamfiles) or np.any(np.logical_and(testY != 0,testY != 1)):
    metric = r2_score
else:
    metric = average_precision_score

np.savetxt(opts.out+'/train_predictions.txt',train_pred)
np.savetxt(opts.out+'/test_predictions.txt',test_pred)
with open(opts.out+'/performance.txt','w') as f:
    f.write('Train performance: '+str(metric(trainY,train_pred))+'\n')
with open(opts.out+'/performance.txt','a') as f:
    f.write('Test performance: '+str(metric(testY,test_pred))+'\n')
    
print('----------------------------------------')
print('Performing Differential Motif Evaluation')
print('----------------------------------------')
X,X_bg,seqsamples = motif2test(opts.motifDB,opts.bg)
pred_mat = np.zeros((X.shape[0],trainY.shape[1]))
pred_null_mat = np.zeros((X_bg.shape[0],trainY.shape[1]))
for mi,model in enumerate(model_folders):
    model_folder = opts.out + "/" + model.split('/')[-1]
    trained_cnn = keras.models.load_model(model_folder+"/model.h5")
    pred_mat += accuracies[model_folder]*trained_cnn.predict(X)
    pred_null_mat += accuracies[model_folder]*trained_cnn.predict(X_bg)
    del trained_cnn
pred_mat = pred_mat/sum(accuracies.values())
pred_null_mat = pred_null_mat/sum(accuracies.values())


motifs_seqs = np.array([s.split('_bg')[0] for s in seqsamples])
scores = np.zeros((len(motifs_seqs),pred_mat.shape[1]))
motifs = list(set([s.split('_bg')[0] for s in seqsamples]))
comps = [(comp.strip().split('vs')[0].split(','),
          comp.strip().split('vs')[1].split(','))
          for comp in opts.comparisons]
for comp in comps:
    c1 = np.array([li for li,l in enumerate(opts.labels) if l in comp[0]])
    c2 = np.array([li for li,l in enumerate(opts.labels) if l in comp[1]])
    with open(opts.out+'_'+'-'.join(comp[0])+'vs'+'-'.join(comp[1])+'.txt','w') as f:
        f.write('TF\tL2FC\tWilcoxonES\t-log10(p-value)\t-log10(adj p-value)\n')
        for i,motif in enumerate(motifs):
            ind_motif_seqs = np.where(motif == motifs_seqs)[0]
            #for a given motif:
            # fold change relative to background
            # and across cell contexts
            if c1.shape[0] == 1:
                num = pred_mat[ind_motif_seqs,c1]
                denom = pred_null_mat[:,c1].reshape((-1,))
                fc1 = num/denom
            elif c1.shape[0] > 1:
                num = np.mean([pred_mat[ind_motif_seqs,ind] for ind in c1],axis=0).reshape((-1,))
                denom = np.mean([pred_null_mat[:,ind] for ind in c1],axis=0).reshape((-1,))
                fc1 = num/denom
            if c2.shape[0] == 0:
                stat,pval = wilcoxon(num,denom)
                adjpval = pval*len(motifs)
                f.write(motif+'\t'+str(np.log2(np.mean(fc1)))+'\t'+str(stat)+'\t'+str(-np.log10(pval))+'\t'+str(-np.log10(adjpval))+'\n')
            if c2.shape[0] == 1:
                num = pred_mat[ind_motif_seqs,c2]
                denom = pred_null_mat[:,c2].reshape((-1,))
                fc2 = num/denom
            elif c2.shape[0] > 1:
                num = np.mean([pred_mat[ind_motif_seqs,ind] for ind in c2],axis=0).reshape((-1,))
                denom = np.mean([pred_null_mat[:,ind] for ind in c2],axis=0).reshape((-1,))            
                fc2 = num/denom
            if c1.shape[0] == 0:
                stat,pval = wilcoxon(num,denom)
                adjpval = pval*len(motifs)
                f.write(motif+'\t'+str(np.log2(np.mean(fc2)))+'\t'+str(stat)+'\t'+str(-np.log10(pval))+'\t'+str(-np.log10(adjpval))+'\n')
            if c1.shape[0] != 0 and c2.shape[0] != 0:
                stat,pval = wilcoxon(fc1,fc2)
                adjpval = pval*len(motifs)
                f.write(motif+'\t'+str(np.log2(np.mean(fc1/fc2)))+'\t'+str(stat)+'\t'+str(-np.log10(pval))+'\t'+str(-np.log10(adjpval))+'\n')
