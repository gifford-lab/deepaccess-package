from scipy.stats import wilcoxon, norm
import numpy as np
from deepaccess.ensemble_utils import *

def motif2test(motiffile, backgroundfile, p=None):
    motifs = open(motiffile).read().split(">")[1:]
    motifmats = {}
    for motif in motifs:
        motifname = motif.strip().split("\n")[0]
        motiflines = motif.strip().split("\n")[1:]

        motif_mat = np.zeros((len(motiflines), 4))
        for i, line in enumerate(motiflines):
            counts = np.array([float(c) for c in line.split("\t")])
            motif_mat[i, :] = counts
        motifmats[motifname] = motif_mat
    names = []
    all_backgrounds = []
    for motifname in motifmats.keys():
        motif = motifmats[motifname]
        backgrounds = fa_to_onehot(backgroundfile)
        if p is None:
            start = int(backgrounds.shape[1] / 2 - motif.shape[0] / 2)
        else:
            start = p
        for bi in range(backgrounds.shape[0]):
            for pos in range(motif.shape[0]):
                consensus_char = np.argmax(motif[pos, :])
                backgrounds[bi, pos + start, :] = 0
                backgrounds[bi, pos + start, consensus_char] = 1.0
            names.append(motifname + "_EPEDistribution_" + str(bi))
        all_backgrounds.append(backgrounds)

    null_backgrounds = fa_to_onehot(backgroundfile)
    return(
        np.concatenate(all_backgrounds, axis=0),
        np.array(null_backgrounds),
        names
        )


def fasta2test(fastafile, backgroundfile, p=None):
    fastanames = [f.split("\n")[0]
                  for f in open(fastafile).read().split(">")[1:]]
    fastas = fa_to_onehot(fastafile, make_uniform_length=False)
    names = []
    all_backgrounds = []
    for fi, fastaname in enumerate(fastanames):
        fasta = fastas[fi]
        backgrounds = fa_to_onehot(backgroundfile)
        if p is None:
            start = int(backgrounds.shape[1] / 2 - fasta.shape[0] / 2)
        else:
            start = p
        for bi in range(backgrounds.shape[0]):
            for pos in range(fasta.shape[0]):
                # if fasta is not N replace background,
                # otherwise keep background the same
                if np.max(fasta[pos, :]) != 0.25:
                    consensus_char = np.argmax(fasta[pos, :])
                    backgrounds[bi, pos + start, :] = 0
                    backgrounds[bi, pos + start, consensus_char] = 1.0
            names.append(fastaname + "_EPEDistribution_" + str(bi))
        all_backgrounds.append(backgrounds)
    null_backgrounds = fa_to_onehot(backgroundfile)
    return(
        np.concatenate(all_backgrounds, axis=0),
        np.array(null_backgrounds),
        names
        )


def rank_ratio_statistic(num, denom):
    sorted_by_ranks = np.argsort(np.abs(np.log2(num / denom)))
    Wp = 0
    Wn = 0
    n = num.shape[0]
    for rank, ri in enumerate(sorted_by_ranks):
        if np.sign(num[ri] - denom[ri]) == 1:
            Wp += rank + 1
        elif np.sign(num[ri] - denom[ri]) == -1:
            Wn += rank + 1
    W = min(Wp, Wn)
    num = (W - n * (n + 1) / 4)
    denom = np.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    return (
        norm.cdf(x=(num / denom)) * 2
    )


def ExpectedPatternEffect(predict_function, class_ind, X_p, X, seqsets):
    fx_p = predict_function(X_p)
    fx = predict_function(X)
    patternseqs = np.array([s.split("_EPEDistribution")[0] for s in seqsets])
    patterns = list(
        sorted(set([s.split("_EPEDistribution")[0] for s in seqsets]))
    )
    stats = []
    pvals = []
    ratio_pvals = []
    adj_pvals = []
    EPEs = []
    for pattern in patterns:
        ind_pattern_seqs = np.where(pattern == patternseqs)[0]
        if len(class_ind) == 1:
            num = fx_p[ind_pattern_seqs, class_ind[0]].reshape((-1,))
            denom = fx[:, class_ind[0]].reshape((-1,))
        else:
            num = np.mean(
                [fx_p[ind_pattern_seqs, ind] for ind in class_ind], axis=0
            ).reshape((-1,))
            denom = np.mean(
                [fx[:, ind] for ind in class_ind], axis=0
            ).reshape((-1,))
        fc = num / denom
        ratio_pvals.append(rank_ratio_statistic(num, denom))
        stat, pval = wilcoxon(num, denom)

        EPEs.append(np.mean(np.log2(fc)))
        pvals.append(pval)
        adj_pvals.append(pval * len(patterns))
        stats.append(stat)
    return (
        fx_p,
        fx,
        {
            "Pattern": patterns,
            "ExpectedPatternEffect": EPEs,
            "Significance": ratio_pvals,
            "AdjustedSignificance": [r * len(patterns) for r in ratio_pvals],
            "WilcoxonTestStatistic": stats,
            "WilcoxonSignificance": pvals,
            "WilcoxonSignificanceAdj": adj_pvals,
        },
    )


def DifferentialExpectedPatternEffect(
    predict_function, class1_ind, class2_ind, X_p, X, seqsets
):
    fx_p, fx, c1_EPEdata = ExpectedPatternEffect(
        predict_function, class1_ind, X_p, X, seqsets
    )
    patternseqs = np.array([s.split("_EPEDistribution")[0] for s in seqsets])
    patterns = list(
        set(sorted([s.split("_EPEDistribution")[0] for s in seqsets]))
    )
    stats = []
    pvals = []
    adj_pvals = []
    ratio_pvals = []
    DiffEPEs = []
    for pattern in patterns:
        ind_pattern_seqs = np.where(pattern == patternseqs)[0]
        if len(class1_ind) == 1:
            num = fx_p[ind_pattern_seqs, class1_ind].reshape((-1,))
            denom = fx[:, class1_ind].reshape((-1,))
        else:
            num = np.mean(
                [fx_p[ind_pattern_seqs, ind] for ind in class1_ind], axis=0
            ).reshape((-1,))
            denom = np.mean(
                [fx[:, ind] for ind in class1_ind], axis=0
            ).reshape((-1,))
        fc1 = num / denom
        if len(class2_ind) == 1:
            num = fx_p[ind_pattern_seqs, class2_ind].reshape((-1,))
            denom = fx[:, class2_ind].reshape((-1,))
        else:
            num = np.mean(
                [fx_p[ind_pattern_seqs, ind] for ind in class2_ind], axis=0
            ).reshape((-1,))
            denom = np.mean(
                [fx[:, ind] for ind in class2_ind], axis=0
            ).reshape((-1,))
        fc2 = num / denom
        stat, pval = wilcoxon(fc1, fc2)

        DiffEPEs.append(np.mean(np.log2(fc1 / fc2)))
        ratio_pvals.append(rank_ratio_statistic(fc1, fc2))

        pvals.append(pval)
        adj_pvals.append(pval * len(patterns))
        stats.append(stat)
    return (
        fx_p,
        fx,
        {
            "Pattern": patterns,
            "DifferentialExpectedPatternEffect": DiffEPEs,
            "Significance": ratio_pvals,
            "AdjustedSignificance": [r * len(patterns) for r in ratio_pvals],
            "WilcoxonTestStatistic": stats,
            "WilcoxonSignificance": pvals,
            "WilcoxonSignificanceAdj": adj_pvals,
        },
    )


def SubtractDifferentialExpectedPatternEffect(
    predict_function, class1_ind, class2_ind, X_p, X, seqsets
):
    fx_p, fx, c1_EPEdata = ExpectedPatternEffect(
        predict_function, class1_ind, X_p, X, seqsets
    )
    patternseqs = np.array([s.split("_EPEDistribution")[0] for s in seqsets])
    patterns = list(
        set(sorted([s.split("_EPEDistribution")[0] for s in seqsets]))
    )
    stats = []
    pvals = []
    adj_pvals = []
    ratio_pvals = []
    DiffEPEs = []
    for pattern in patterns:
        ind_pattern_seqs = np.where(pattern == patternseqs)[0]
        if len(class1_ind) == 1:
            num = fx_p[ind_pattern_seqs, class1_ind].reshape((-1,))
            denom = fx[:, class1_ind].reshape((-1,))
        else:
            num = np.mean(
                [fx_p[ind_pattern_seqs, ind] for ind in class1_ind], axis=0
            ).reshape((-1,))
            denom = np.mean(
                [fx[:, ind] for ind in class1_ind], axis=0
            ).reshape((-1,))
        fc1 = num / denom
        if len(class2_ind) == 1:
            num = fx_p[ind_pattern_seqs, class2_ind].reshape((-1,))
            denom = fx[:, class2_ind].reshape((-1,))
        else:
            num = np.mean(
                [fx_p[ind_pattern_seqs, ind] for ind in class2_ind], axis=0
            ).reshape((-1,))
            denom = np.mean(
                [fx[:, ind] for ind in class2_ind], axis=0
            ).reshape((-1,))
        fc2 = num / denom
        stat, pval = wilcoxon(fc1, fc2)

        DiffEPEs.append(np.mean(fc1 - fc2))
        ratio_pvals.append(rank_ratio_statistic(fc1, fc2))

        pvals.append(pval)
        adj_pvals.append(pval * len(patterns))
        stats.append(stat)
    return (
        fx_p,
        fx,
        {
            "Pattern": patterns,
            "DifferentialExpectedPatternEffect": DiffEPEs,
            "Significance": ratio_pvals,
            "AdjustedSignificance": [r * len(patterns) for r in ratio_pvals],
            "WilcoxonTestStatistic": stats,
            "WilcoxonSignificance": pvals,
            "WilcoxonSignificanceAdj": adj_pvals,
        },
    )


def SubtractExpectedPatternEffect(
        predict_function, class_ind, X_p, X, seqsets
):
    fx_p = predict_function(X_p)
    fx = predict_function(X)
    patternseqs = np.array([s.split("_EPEDistribution")[0] for s in seqsets])
    patterns = list(
        sorted(set([s.split("_EPEDistribution")[0] for s in seqsets]))
    )
    stats = []
    pvals = []
    ratio_pvals = []
    adj_pvals = []
    EPEs = []
    for pattern in patterns:
        ind_pattern_seqs = np.where(pattern == patternseqs)[0]
        if len(class_ind) == 1:
            num = fx_p[ind_pattern_seqs, class_ind[0]].reshape((-1,))
            denom = fx[:, class_ind[0]].reshape((-1,))
        else:
            num = np.mean(
                [fx_p[ind_pattern_seqs, ind] for ind in class_ind], axis=0
            ).reshape((-1,))
            denom = np.mean(
                [fx[:, ind] for ind in class_ind], axis=0
            ).reshape((-1,))
        fc = num - denom
        ratio_pvals.append(rank_ratio_statistic(num, denom))
        stat, pval = wilcoxon(num, denom)

        EPEs.append(np.mean(fc))
        pvals.append(pval)
        adj_pvals.append(pval * len(patterns))
        stats.append(stat)
    return (
        fx_p,
        fx,
        {
            "Pattern": patterns,
            "ExpectedPatternEffect": EPEs,
            "Significance": ratio_pvals,
            "AdjustedSignificance": [r * len(patterns) for r in ratio_pvals],
            "WilcoxonTestStatistic": stats,
            "WilcoxonSignificance": pvals,
            "WilcoxonSignificanceAdj": adj_pvals,
        },
    )
