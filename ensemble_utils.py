import os
import numpy as np 


def ensure_dir(file_path):
    #directory = os.path.dirname(file_path)
    if not os.path.exists(file_path):
        os.makedirs(file_path)

def act_to_class(act):
    y = []
    header = True
    for line in open(act):
        if header:
            header = False
            continue
        data  = line.strip().split()
        y.append([int(d) for d in data[1:]])
    return np.array(y)

def fa_to_onehot(fa):
    alpha = ['A','C','G','T']
    sequences = open(fa).read().split(">")[1:]
    seqdict = [seq.strip().split("\n")[1] for seq in sequences]
    seq_mat = []
    slen = max([len(seq) for seq in seqdict])
    for i,seqc in enumerate(seqdict):
        seq = np.zeros((slen,4))
        for j,c in enumerate(seqc.upper()):
            if c not in alpha:
                seq[j,:] = 0.25
            else:
                aind = alpha.index(c)
                seq[j,aind] = 1
        seq_mat.append(seq)
    return np.array(seq_mat)

def motif2test(motiffile,backgroundfile):
    motifs = open(motiffile).read().split('>')[1:]
    motifmats = {}
    for motif in motifs:
        motifname = motif.strip().split('\n')[0]
        motiflines = motif.strip().split('\n')[1:]
        
        motif_mat = np.zeros((len(motiflines),4))
        for i,line in enumerate(motiflines):
            counts = np.array([float(c) for c in line.split('\t')])
            #counts = counts + np.min(counts)
            #probs = counts/np.sum(counts)
            motif_mat[i,:] = counts
        motifmats[motifname] = motif_mat
    names = []
    all_backgrounds = []
    for motifname in motifmats.keys():
        motif = motifmats[motifname]
        backgrounds = fa_to_onehot(backgroundfile)
        start = int(backgrounds.shape[1]/2 - motif.shape[0]/2)
        for bi in range(backgrounds.shape[0]):
            for pos in range(motif.shape[0]):
                consensus_char = np.argmax(motif[pos,:])
                backgrounds[bi,pos+start,:] = 0
                backgrounds[bi,pos+start,consensus_char] = 1.0
            names.append(motifname+'_bg'+str(bi))
        all_backgrounds.append(backgrounds)
        
    null_backgrounds = fa_to_onehot(backgroundfile)
    return np.concatenate(all_backgrounds,axis=0),np.array(null_backgrounds),names

def motif2testCombinatorial(motiffile,backgroundfile,primary_motif,motifsubset=None):
    motifs = open(motiffile).read().split('>')[1:]
    motifmats = {}
    for motif in motifs:
        motifname = motif.strip().split('\n')[0]
        motiflines = motif.strip().split('\n')[1:]
        
        motif_mat = np.zeros((len(motiflines),4))
        for i,line in enumerate(motiflines):
            counts = np.array([float(c) for c in line.split('\t')])
            probs = counts/np.sum(counts)
            motif_mat[i,:] = probs
        motifmats[motifname] = motif_mat
    if motifsubset != None:
        motifset = [line.strip() for line in open(motifsubset)]
        motifmats = {motif:motifmats[motif] for motif in motifset}
        
    names = []
    all_backgrounds = []
    for motifname1 in [primary_motif]:
        for motifname2 in motifmats.keys():
            motif1 = motifmats[motifname1]            
            motif2 = motifmats[motifname2]

            backgrounds = fa_to_onehot(backgroundfile)
            start = int(backgrounds.shape[1]/2 - motif1.shape[0]/2)
            end = start + motif1.shape[0]
            
            for motif2_start in range(0,start-motif2.shape[0]):
                backgrounds = fa_to_onehot(backgroundfile)
                for bi in range(backgrounds.shape[0]):
                    for pos in range(motif1.shape[0]):
                        consensus_char = np.argmax(motif1[pos,:])
                        backgrounds[bi,pos+start,:] = 0
                        backgrounds[bi,pos+start,consensus_char] = 1.0
                    for pos in range(motif2.shape[0]):
                        consensus_char = np.argmax(motif2[pos,:])
                        backgrounds[bi,pos+motif2_start,:] = 0
                        backgrounds[bi,pos+motif2_start,consensus_char] = 1.0
                    names.append(motifname1+'_bg'+str(bi)+'_'+motifname2+'_'+str(motif2_start))
                all_backgrounds.append(backgrounds)
            for motif2_start in range(end,backgrounds.shape[1]-motif2.shape[0]):
                backgrounds = fa_to_onehot(backgroundfile)
                for bi in range(backgrounds.shape[0]):
                    for pos in range(motif1.shape[0]):
                        consensus_char = np.argmax(motif1[pos,:])
                        backgrounds[bi,pos+start,:] = 0
                        backgrounds[bi,pos+start,consensus_char] = 1.0
                    for pos in range(motif2.shape[0]):
                        consensus_char = np.argmax(motif2[pos,:])
                        backgrounds[bi,pos+motif2_start,:] = 0
                        backgrounds[bi,pos+motif2_start,consensus_char] = 1.0
                    names.append(motifname1+'_bg'+str(bi)+'_'+motifname2+'_'+str(motif2_start))
                all_backgrounds.append(backgrounds)


            for motif2_start in range(0,start-motif2.shape[0]):
                backgrounds = fa_to_onehot(backgroundfile)
                for bi in range(backgrounds.shape[0]):
                    for pos in range(motif2.shape[0]):
                        consensus_char = np.argmax(motif2[pos,:])
                        backgrounds[bi,pos+motif2_start,:] = 0
                        backgrounds[bi,pos+motif2_start,consensus_char] = 1.0
                    names.append('None_bg'+str(bi)+'_'+motifname2+'_'+str(motif2_start))
                all_backgrounds.append(backgrounds)
            for motif2_start in range(end,backgrounds.shape[1]-motif2.shape[0]):
                backgrounds = fa_to_onehot(backgroundfile)
                for bi in range(backgrounds.shape[0]):
                    for pos in range(motif2.shape[0]):
                        consensus_char = np.argmax(motif2[pos,:])
                        backgrounds[bi,pos+motif2_start,:] = 0
                        backgrounds[bi,pos+motif2_start,consensus_char] = 1.0
                    names.append('None_bg'+str(bi)+'_'+motifname2+'_'+str(motif2_start))
                all_backgrounds.append(backgrounds)

    return np.concatenate(all_backgrounds,axis=0),names
