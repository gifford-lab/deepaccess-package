import os
import numpy as np
import sys


def ensure_dir(file_path, exit_if_exists=False):
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    elif exit_if_exists:
        sys.exit("Error: Folder " + file_path + " already exists")
    else:
        pass


def act_to_class(act):
    y = []
    header = True
    for line in open(act):
        if header:
            header = False
            continue
        data = line.strip().split()
        y.append([int(d) for d in data[1:]])
    return np.array(y)


def fa_to_onehot(fa, make_uniform_length=True):
    alpha = ["A", "C", "G", "T"]
    sequences = open(fa).read().split(">")[1:]
    seqdict = [seq.strip().split("\n")[1] for seq in sequences]
    seq_mat = []
    if make_uniform_length:
        slen = max([len(seq) for seq in seqdict])
    for i, seqc in enumerate(seqdict):
        if not make_uniform_length:
            slen = len(seqc)
        seq = np.zeros((slen, 4))
        for j, c in enumerate(seqc.upper()):
            if c not in alpha:
                seq[j, :] = 0.25
            else:
                aind = alpha.index(c)
                seq[j, aind] = 1
        seq_mat.append(seq)
    if make_uniform_length:
        return np.array(seq_mat)
    else:
        return seq_mat
