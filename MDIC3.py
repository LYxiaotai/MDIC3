import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import scipy.stats as stat
import math
import os
import scipy
import scipy.stats as stat
import time
import sys
import numpy.linalg as la
import scipy.sparse as sp
import warnings
warnings.filterwarnings("ignore")


begin=time.asctime()


param = {}
for i in range(1, len(sys.argv)):
    t = sys.argv[i].split("=")
    param[t[0].lower()] = t[1]

help_msg = """
usage1: 
    python MDIC3.py -exp=scRNA_expression_file -label=cell_label_file -grnchoose=’other’ -grn=grn_file -out=results_output_fold
usage2:
    python MDIC3.py -exp=scRNA_expression_file -label=cell_label_file -grnchoose=’GNIPLR’ -process= process_value -step=GRN_calculation_step -out=results_output_fold

Options and arguments:
    -exp: the input single-cell gene expression data.
    -label: the input scRNAseq metadata.
    -grnchoose: availability of gene regulatory networks.
    -process: if -grnchoose =='GNIPLR', the user must select the number of work processes used.
    -step: if -grnchoose =='GNIPLR', the user must select the step size for GNIPLR block calculation.
    -grn: if -grnchoose =='other', the user must provide the gene regulatory network adjacency matrix *.txt file.
    -out: the directory to store the MDIC3 results. 
"""

if "-help" in param.keys() or "-h" in param.keys():
    print(help_msg)

if "-exp" not in param.keys() or "-label" not in param.keys() or "-grnchoose" not in param.keys():
    print("Parameter missing!")
    print(param.keys())
    print(help_msg)
    exit()

if param["-grnchoose"] == 'other':
    GRNtxt = param["-grn"]
    if "-grn" not in param.keys():
        print("Parameter missing!")
        print('The user must must provide the gene regulatory network adjacency matrix txt file.')
        print(help_msg)
        exit()

if param["-grnchoose"] == 'GNIPLR':
    import multiprocessing
    from multiprocessing import Lock, Queue, Process, Manager
    process = int(param["-process"])
    step = int(param["-step"])
    if "-process" not in param.keys() or "-step" not in param.keys():
        print("Parameter missing!")
        print('The user must select the number of work processes used and the step size for GNIPLR block calculation.')
        print(help_msg)
        exit()


fold = param["-out"]
if not os.path.exists(fold):
    os.mkdir(fold)

expression_txt = param["-exp"]
labelname = param["-label"]
choose = param["-grnchoose"]




def LassoRegression(degree, alpha):
    return Pipeline([
        ("poly", PolynomialFeatures(degree=degree)),
        ("std_scaler", StandardScaler()),
        ("lasso_reg", Lasso(alpha=alpha))])


def lasso_regress(gene_1, gene_2, degree, alpha):

    le_g1 = len(gene_1)
    t = np.linspace(min(gene_1), max(gene_1), num=le_g1)
    a1, b1 = (np.array(x).reshape(-1, 1) for x in zip(*sorted(zip(gene_1, gene_2), key=lambda x: x[0])))
    lasso1_reg = LassoRegression(degree, alpha)
    lasso1_reg.fit(a1, t)
    gene1 = lasso1_reg.predict(a1)
    gene2 = lasso1_reg.predict(b1)

    return gene1, gene2


def rss_calculate(model, right_v, left1_v, left2_v):

    model.fit(left1_v, right_v)
    rssu_1 = np.sum((model.predict(left1_v) - right_v) ** 2)
    model.fit(left2_v, right_v)
    rssr_1 = np.sum((model.predict(left2_v) - right_v) ** 2)

    return rssu_1, rssr_1


def granger(gene1, gene2, cellnum):

    sort_1 = list(sorted(zip(gene1, gene2), key=lambda x: x[0]))

    a1, b1 = (list(x) for x in zip(*sort_1))
    a1_cha = ((300 * pd.Series(a1)) / 90001).tolist()
    b1_cha = ((300 * pd.Series(b1)) / 90001).tolist()

    model = LinearRegression()
    b1_t_1 = np.array(b1_cha[0:cellnum - 1])
    a1_t1_t2_b1_1 = sm.add_constant(np.array(pd.concat([pd.Series(a1_cha),
                                                        pd.Series(b1_cha)], axis=1)[1:]))
    b1_t1_t2_1 = sm.add_constant(np.array(b1_cha[1:cellnum]))
    b1_rssu_1, b1_rssr_1 = rss_calculate(model, b1_t_1, a1_t1_t2_b1_1, b1_t1_t2_1)
    f1_1 = ((b1_rssr_1 - b1_rssu_1) / 1) / (b1_rssu_1 / (cellnum - 3))
    p1_1 = scipy.stats.f.sf(f1_1, 1, cellnum - 3)

    b1_t = np.array(b1_cha[0:cellnum - 2])
    a1_t1_t2_b1 = sm.add_constant(np.array(pd.concat([pd.Series(a1_cha)[1:cellnum - 1].reset_index(drop=True),
                                                      pd.Series(a1_cha)[2:cellnum].reset_index(drop=True),
                                                      pd.Series(b1_cha)[1:cellnum - 1].reset_index(drop=True),
                                                      pd.Series(b1_cha)[2:cellnum].reset_index(drop=True)], axis=1)))
    b1_t1_t2 = sm.add_constant(np.array(pd.concat([pd.Series(b1_cha)[1:cellnum - 1].reset_index(drop=True),
                                                   pd.Series(b1_cha)[2:cellnum].reset_index(drop=True)], axis=1)))
    b1_rssu_2, b1_rssr_2 = rss_calculate(model, b1_t, a1_t1_t2_b1, b1_t1_t2)
    f1_2 = ((b1_rssr_2 - b1_rssu_2) / 2) / (b1_rssu_2 / (cellnum - 2 - 2 - 1))
    p1_2 = scipy.stats.f.sf(f1_2, 2, cellnum - 2 - 2 - 1)

    b1_t_3 = np.array(b1_cha[0:cellnum - 3])
    a1_t1_t2_b1_3 = sm.add_constant(np.array(pd.concat([pd.Series(a1_cha)[1:cellnum - 2].reset_index(drop=True),
                                                        pd.Series(a1_cha)[3:cellnum].reset_index(drop=True),
                                                        pd.Series(a1_cha)[2:cellnum - 1].reset_index(drop=True),
                                                        pd.Series(b1_cha)[1:cellnum - 2].reset_index(drop=True),
                                                        pd.Series(b1_cha)[3:cellnum].reset_index(drop=True),
                                                        pd.Series(b1_cha)[2:cellnum - 1].reset_index(drop=True)],
                                                       axis=1)))
    b1_t1_t2_3 = sm.add_constant(np.array(pd.concat([pd.Series(b1_cha)[1:cellnum - 2].reset_index(drop=True),
                                                     pd.Series(b1_cha)[3:cellnum].reset_index(drop=True),
                                                     pd.Series(b1_cha)[2:cellnum - 1].reset_index(drop=True)], axis=1)))
    b1_rssu_3, b1_rssr_3 = rss_calculate(model, b1_t_3, a1_t1_t2_b1_3, b1_t1_t2_3)
    f1_3 = ((b1_rssr_3 - b1_rssu_3) / 3) / (b1_rssu_3 / (cellnum - 7))
    p1_3 = scipy.stats.f.sf(f1_3, 3, cellnum - 7)

    p_1 = min([p1_1, p1_2, p1_3])

    return p_1


def calculate_block(spaceblocki, gene, AA, gene_exp, ax_nl, cellnum, gene_num,pear):
    i1, i2, j1, j2 = [int(spaceblocki[i]) for i in range(len(spaceblocki))]
    # print(i1,i2,j1,j2)
    PP = sp.dok_matrix((gene_num,gene_num))
    for i in range(i1, i2):
        for j in range(j1, j2):
            if i == j:
                continue
            entry_11 = gene[i]
            entry_21 = gene[j]
            mp = list(pear[i])
            if 1 in pear[i]:
                mp.remove(1)
            mp = abs(np.array(mp))
            mpmax = max(mp)
            if abs(stat.pearsonr(AA[i],AA[j])[0]) < mpmax * 0.3:
                continue

            Gene1, Gene2 = lasso_regress(gene_exp[entry_11], gene_exp[entry_21], 30, ax_nl)

            if (Gene1 == Gene2).all() == True:
                continue

            p = granger(Gene1, Gene2, cellnum)

            Pi = gene.index(entry_11)
            Pj = gene.index(entry_21)
            PP[Pi,Pj] += p

    return PP


def MDIC3(AA, GRN):

    U1, S1, V1 = np.linalg.svd(AA)
    shapeA = np.shape(AA)
    d = min(shapeA) - len(S1)
    dd = []
    for i in range(d):
        dd.append(0)
    S11 = np.append(S1, dd)

    S11 = np.diag(S11)
    if shapeA[0] < shapeA[1]:
        b = np.zeros((shapeA[0], shapeA[1] - shapeA[0]))
        S11 = np.column_stack((S11, b))
    if shapeA[0] > shapeA[1]:
        b = np.zeros((shapeA[0] - shapeA[1], shapeA[1]))
        S11 = np.row_stack((S11, b))

    T = np.dot(GRN, S11)
    Tp = np.linalg.pinv(T)
    M = np.dot(Tp, AA)

    return M


def cellular_score(M):

    llen = len(M)
    N = np.zeros((llen, llen))
    
    for i in range(llen):
        for j in range(llen):
            if M[i][j] >= 0:
                N[i][j] += M[i][j]
            if M[i][j] < 0:
                N[j][i] += abs(M[i][j])

    return N


def celltype_score(labelsl, indexa, M):

    llen = len(labelsl)
    N = np.zeros((llen, llen))
    for l1 in labelsl:
        I = indexa[l1]
        for l2 in labelsl:
            J = indexa[l2]
            if l1 == l2:
                for i in I:
                    for j in J:
                        N[labelsl.index(l1)][labelsl.index(l2)] += abs(M[i][j])
            if l1 != l2:
                for i in I:
                    for j in J:
                        if M[i][j] >= 0:
                            N[labelsl.index(l1)][labelsl.index(l2)] += M[i][j]
                        if M[i][j] < 0:
                            N[labelsl.index(l2)][labelsl.index(l1)] += abs(M[i][j])
    Nn = np.zeros((llen, llen))
    for i in range(len(N)):
        for j in range(len(N)):
            if N[i][j] == 0:
                continue
            if math.log(N[i][j], 10)<=0:
                continue
            Nn[i][j] += math.log(N[i][j], 10)

    return Nn


if __name__ == '__main__':

    print("Begin time: " + begin)

    if choose == 'GNIPLR':

        print('Start calculating the GRN adjacency matrix...')

        f = open(expression_txt)
        gene = []
        cellname = []
        gene_exp = {}
        AA = []
        flag = 0
        for p in f:
            t = p.split()
            flag += 1
            if flag == 1:
                cellnum = len(t)
                cellname = t
                continue
            gene.append(t[0])
            gene_exp[t[0]] = [float(t[i]) for i in range(1, len(t))]
            tt = list(map(float, t[1:]))
            AA.append(tt)
        f.close()
        gene_num = len(gene_exp)

        Pear = np.zeros((len(AA),len(AA)))
        for i in range(len(AA)):
            for j in range(i,len(AA)):
                Pear[i][j] += stat.pearsonr(AA[i], AA[j])[0]
        Pear = np.where(Pear,Pear,Pear.T)


        space = list(np.linspace(start=0, stop=step * int(gene_num / step), num=int(gene_num / step) + 1))
        if (gene_num / step == int(gene_num / step)) == False:
            space.append(gene_num)
        space1 = []
        for i in range(len(space) - 1):
            space1.append((space[i], space[i + 1]))
        space_block = []
        for i in space1:
            for j in space1:
                space_block.append(i + j)

        pool = multiprocessing.Pool(process)

        GRN = sp.dok_matrix((gene_num, gene_num))

        rest = []

        for si in space_block:
            rest.append(pool.apply_async(func=calculate_block, args=(si, gene, AA, gene_exp, 0.1, cellnum, gene_num,Pear,)))

        for r in rest:
            GRN+=r.get()

        print('Waiting for all subprocesses done...')
        pool.close()
        pool.join()

        print("Complete the GRN calculation.")
        print("Complete the GRN calculation time: " + time.asctime())

        GRN = GRN.toarray()
        np.savetxt(fold + os.sep + 'GRN_GNIPLR.txt', GRN)
      



    if choose == 'other':

        f = open(expression_txt)
        AA = []
        cellname = []
        flag = 0
        for p in f:
            t = p.split()
            flag += 1
            if flag == 1:
                cellname = t
                continue
            a = t[1:]
            tt = list(map(float, a))
            AA.append(tt)
        f.close()

        print('The GRN txt file is being read...')
        f = open(GRNtxt)
        GRN = []
        for p in f:
            t = p.split()
            tt = [float(i) for i in t]
            GRN.append(tt)
        f.close()
        GRN = np.array(GRN)
        
    print('Start calculating cell-cell communication...')

    f = open(labelname)
    label = []
    labels = set()
    for p in f:
        t = p.split()
        if len(t)>2:
            print('May cause error calculation in further steps, please check the cell type name in metadata file.')
            sys.exit()
        label.append(t[1])
        labels.add(str(t[1]))
    labelsl = list(labels)
    indexa = {}
    for i in labelsl:
        indexa[i] = []
    flag = 0
    for i in label:
        flag += 1
        indexa[i].append(flag - 1)
    f.close()

    if len(cellname)!=len(label):
        print('Please check the exp file and the metadata file, the two files contain different numbers of cells.')
        sys.exit()

    cellular_adjacency = cellular_score(MDIC3(AA, GRN))
    type_adjacency = celltype_score(labelsl, indexa, cellular_adjacency)
    
    print("Complete the cell-cell communication calculation.")

    print("Complete the cell-cell communication time: " + time.asctime())

    print('Saving MDIC3 results...')

    fw = open(fold + os.sep + 'celltype_communication.txt', 'w')
    fw.write('\t')
    for tt in range(len(labelsl)):
        if tt < len(labelsl) - 1:
            fw.write(labelsl[tt] + '\t')
        if tt == len(labelsl) - 1:
            fw.write(labelsl[tt] + '\n')
    for i in range(len(type_adjacency)):
        Ni = list(type_adjacency[i])
        fw.write(str(labelsl[i]) + '\t')
        for ni in range(len(Ni)):
            if ni < len(Ni) - 1:
                fw.write(str(Ni[ni]) + '\t')
            if ni == len(Ni) - 1:
                fw.write(str(Ni[ni]) + '\n')
    fw.close()
    
    
    fw = open(fold + os.sep + 'cellular_communication.txt', 'w')
    fw.write('\t')
    for tt in range(len(cellname)):
        if tt < len(cellname) - 1:
            fw.write(cellname[tt] + '\t')
        if tt == len(cellname) - 1:
            fw.write(cellname[tt] + '\n')
    for i in range(len(cellular_adjacency)):
        Ni = list(cellular_adjacency[i])
        fw.write(str(cellname[i]) + '\t')
        for ni in range(len(Ni)):
            if ni < len(Ni) - 1:
                fw.write(str(Ni[ni]) + '\t')
            if ni == len(Ni) - 1:
                fw.write(str(Ni[ni]) + '\n')
    fw.close()


    print("End time: " + time.asctime())

