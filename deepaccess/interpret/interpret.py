#!/bin/env python
import os

# set local keras environment
os.environ["KERAS_HOME"] = ".keras"

import numpy as np
import subprocess
import argparse
import sys
import matplotlib.pyplot as plt
from deepaccess.train.DeepAccessModel import *
from deepaccess.interpret.ExpectedPatternEffect import *
from deepaccess.interpret.importance_utils import *

def main(args):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-trainDir", "--trainDir", required=True
    )
    parser.add_argument(
        "-fastas", "--fastas", nargs="+", required=False, default=[]
    )
    parser.add_argument(
        "-l", "--labels", nargs="+", required=False, default=[]
    )
    parser.add_argument(
        "-c", "--comparisons", nargs="+", required=False, default=[]
    )
    parser.add_argument(
        "-evalMotifs", "--evalMotifs", default=None, required=False
    )
    parser.add_argument(
        "-evalPatterns", "--evalPatterns", default=None, required=False
    )
    parser.add_argument(
        '-p','--position',default=None, required=False, type=int
    )
    parser.add_argument(
        "-saliency", "--saliency", default=False, action="store_true"
    )
    parser.add_argument(
        "-subtract", "--subtract", default=False, action="store_true"
    )
    parser.add_argument(
        "-bg", "--background", default="default", required=False
    )
    parser.add_argument(
        "-vis", "--makeVis", action="store_true", default=False, required=False
    )
    opts = parser.parse_args(args)

    if opts.background == "default":
        opts.background = pkg_resources.resource_filename(__name__,
                                                          "data/backgrounds.fa")
    elif not os.path.exists(opts.background):
        opts.background = pkg_resources.resource_filename(__name__,
                                                          "data/"+opts.background)
        
    print("-------------------------------------")
    print("         Making Predictions          ")
    print("-------------------------------------")
    DAModel = DeepAccessModel(opts.trainDir)
    DAModel.load(opts.trainDir)

    fasta_hots = []
    for fasta in opts.fastas:
        fasta_hots.append(fa_to_onehot(fasta))
        
    N = 0
    for fasta_ind, fasta in enumerate(opts.fastas):
        fname = fasta.split("/")[-1].split(".fa")[0]
        X = fasta_hots[fasta_ind]
        DA_fa_pred = DAModel.predict(X)
        if opts.saliency:
            saliencyX = []
            comps = [
                [
                    comp.strip().split("vs")[0].split("-"),
                    comp.strip().split("vs")[1].split("-"),
                ]
                for comp in opts.comparisons
            ]
            for ci, comp in enumerate(comps):
                c1 = [li for li, l in enumerate(opts.labels) if l in comp[0]]
                c2 = [li for li, l in enumerate(opts.labels) if l in comp[1]]
                saliencyX.append(DAModel.saliency_input(X, c1, c2))
                if opts.makeVis:
                    seqs = []
                    for seq in open(fasta).read().split(">")[1:]:
                        seqs.append(seq.strip().split("\n")[1].upper())

                    for si in range(saliencyX[ci].shape[0]):
                        f = plt.figure(figsize=(13, 4))
                        plt.rcParams.update({"font.size": 16})
                        ax = f.add_subplot(1, 1, 1)
                        ax.bar(
                            range(saliencyX[ci].shape[1]),
                            saliencyX[ci][si, :].reshape((-1,)),
                        )
                        for i in range(100):
                            if len(seqs[si]) <= i:
                                break
                            if seqs[si][i] == "T":
                                ax.get_children()[i].set_color("r")
                            elif seqs[si][i] == "A":
                                ax.get_children()[i].set_color("g")
                            elif seqs[si][i] == "C":
                                ax.get_children()[i].set_color("b")
                            elif seqs[si][i] == "G":
                                ax.get_children()[i].set_color("gold")
                            else:
                                ax.get_children()[i].set_color("k")
                        plt.xticks(
                                    range(len(seqs[si])), list(seqs[si]), fontsize=10
                        )
                        plt.ylabel("Saliency")
                        plt.savefig(
                            opts.trainDir
                            + "_"
                            + "-".join(comp[0])
                            + "vs"
                            + "-".join(comp[1])
                            + "-saliency"
                            + fname
                            + "_"
                            + str(si)
                            + ".svg"
                        )

                np.savetxt(
                    opts.trainDir
                    + "_"
                    + "-".join(comp[0])
                    + "vs"
                    + "-".join(comp[1])
                    + "-saliency"
                    + fname
                    + ".txt",
                    saliencyX[ci],
                )

                np.savetxt(opts.trainDir + "/" + fname + ".prediction", DA_fa_pred)

    if opts.evalMotifs is not None:
        if opts.evalMotifs == "HMv11_MOUSE":
            opts.evalMotifs = 'data/HMv11_MOUSE.txt'
        motifDB = opts.evalMotifs.split("/")[-1].split(".txt")[0]
        print("----------------------------------------")
        print("Performing Differential Motif Evaluation")
        print("----------------------------------------")
        X, X_bg, seqsamples = motif2test(
            opts.evalMotifs, opts.background, opts.position
        )

        comps = [
            (comp.strip().split("vs")[0].split(","),
             comp.strip().split("vs")[1].split(","))
            for comp in opts.comparisons
        ]
        if opts.subtract:
            method = SubtractExpectedPatternEffect
            diffmethod = SubtractDifferentialExpectedPatternEffect
        else:
            method = ExpectedPatternEffect
            diffmethod = DifferentialExpectedPatternEffect
        if opts.subtract:
            motifDB = "_subtraction_" + motifDB
            
        for comp in comps:
            c1 = np.array([li for li, l in enumerate(opts.labels) if l in comp[0]])
            c2 = np.array([li for li, l in enumerate(opts.labels) if l in comp[1]])
            if c1.shape[0] == 0 and c2.shape[0] == 0:
                sys.exit("Error: invalid comparison " + str(comp))
            elif c1.shape[0] == 0:
                _, _, EPEdata = method(
                    DAModel.predict, c2, X, X_bg, seqsamples
                )
                valuecol = "ExpectedPatternEffect"
            elif c2.shape[0] == 0:
                _, _, EPEdata = method(
                    DAModel.predict, c1, X, X_bg, seqsamples
                )
                valuecol = "ExpectedPatternEffect"
            else:
                _, _, EPEdata = diffmethod(
                    DAModel.predict, c1, c2, X, X_bg, seqsamples
                )
                valuecol = "DifferentialExpectedPatternEffect"
            with open(
                    opts.trainDir
                    + "_EPE_"
                    + "-".join(comp[0])
                    + "vs"
                    + "-".join(comp[1])
                    + "-"
                    + motifDB
                    + ".txt",
                    "w",
            ) as f:
                f.write("\t".join(EPEdata.keys()) + "\n")
                for index in np.argsort(EPEdata[valuecol])[::-1]:
                    f.write(
                        "\t".join([str(EPEdata[k][index])
                                   for k in EPEdata.keys()]) + "\n"
                    )


    if opts.evalPatterns != None:
        patternDB = opts.evalPatterns.split("/")[-1].split(".txt")[0]
        print("----------------------------------------")
        print("Performing Differential Motif Evaluation")
        print("----------------------------------------")
        X, X_bg, seqsamples = fasta2test(
            opts.evalPatterns, opts.background, opts.position
        )

        comps = [
            (comp.strip().split("vs")[0].split(","),
             comp.strip().split("vs")[1].split(","))
            for comp in opts.comparisons
        ]
        if opts.subtract:
            method = SubtractExpectedPatternEffect
            diffmethod = SubtractDifferentialExpectedPatternEffect
        else:
            method = ExpectedPatternEffect
            diffmethod = DifferentialExpectedPatternEffect

        for comp in comps:
            c1 = np.array([li for li, l in enumerate(opts.labels) if l in comp[0]])
            c2 = np.array([li for li, l in enumerate(opts.labels) if l in comp[1]])
            if c1.shape[0] == 0:
                _, _, EPEdata = method(
                    DAModel.predict, c2, X, X_bg, seqsamples
                )
                valuecol = "ExpectedPatternEffect"
            elif c2.shape[0] == 0:
                _, _, EPEdata = method(
                    DAModel.predict, c1, X, X_bg, seqsamples
                )
                valuecol = "ExpectedPatternEffect"
            else:
                _, _, EPEdata = diffmethod(
                    DAModel.predict, c1, c2, X, X_bg, seqsamples
                )
                valuecol = "DifferentialExpectedPatternEffect"
            if opts.subtract:
                patternDB = "_subtraction_" + patternDB
                
            with open(
                    opts.trainDir
                    + "_EPE_"
                    + "-".join(comp[0])
                    + "vs"
                    + "-".join(comp[1])
                    + "-"
                    + patternDB
                    + ".txt",
                    "w",
            ) as f:
                f.write("\t".join(EPEdata.keys()) + "\n")
                for index in np.argsort(EPEdata[valuecol])[::-1]:
                    f.write(
                        "\t".join([str(EPEdata[k][index])
                                   for k in EPEdata.keys()]) + "\n"
                    )

if __name__ == "__main__":
    import sys
    main(sys.argv)
