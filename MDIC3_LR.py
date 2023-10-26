import numpy as np
import os
from scipy.stats import pearsonr
from operator import itemgetter
from statsmodels.stats.multitest import multipletests
import sys


param = {}
for i in range(1, len(sys.argv)):
    t = sys.argv[i].split("=")
    param[t[0].lower()] = t[1]

help_msg = """
usage: 
    python MDIC3_LR.py -exp=scRNA_expression_file -label=cell_label_file -lrdb=LR_information_file -ltype='ligandtype_name' -rtype='receptortype_name' -out=results_output_fold

Options and arguments:
    -exp: the input single-cell gene expression data.
    -label: the input scRNAseq metadata.
    -lrdb: the input ligand-receptor information file.
    -ltype: the cell type name that sends cell-cell communication signals during the target communication progress you want to analyze.
    -rtype: the cell type name that receives cell-cell communication signals during the target communication progress you want to analyze.
    -out: the directory to store the L-R identification results. 

"""


if "-help" in param.keys() or "-h" in param.keys():
    print(help_msg)

if "-exp" not in param.keys() or "-label" not in param.keys() or "-lrdb" not in param.keys() or "-ltype" not in param.keys() or "-rtype" not in param.keys():
    print("Parameter missing!")
    print(param.keys())
    print(help_msg)
    exit()


if "-out" not in param.keys():
    fold = "."
else:
    fold = param["-out"]

if not os.path.exists(fold):
    os.mkdir(fold)

expression_txt = param["-exp"]
labelname = param["-label"]
LRDB_txt = param["-lrdb"]
Ltypename = param["-ltype"]
Rtypename = param["-rtype"]


def check_zero_majority(nums):
    count_zeros = sum(1 for num in nums if num == 0)

    return count_zeros > len(nums) * 0.7


def pearson(x, y):
    if np.std(x) == 0 or np.std(y) == 0:
        return 0
    else:
        corr, p = pearsonr(x, y)
        return corr, p


def FDRP(lrcorr):
    p_values = [v[1] for v in lrcorr.values()]
    corrected_p_values = multipletests(p_values, method='bonferroni')[1]
    p_LR_corr = {}
    index = 0
    for key, value in lrcorr.items():
        p_LR_corr[key] = [value[0], corrected_p_values[index]]
        index += 1
    p2_LR_corr = {key: value for key, value in p_LR_corr.items() if value[1] <= 0.01}

    return p2_LR_corr


def SortLR(target, type_cell, cellname, LRexp):
    cell0 = type_cell[target[0]]
    cell0index = []
    for i in cell0:
        cell0index.append(cellname.index(i))
    cell1 = type_cell[target[1]]
    cell1index = []
    for i in cell1:
        cell1index.append(cellname.index(i))

    targetexp = {}
    for key, value in LRexp.items():
        new_lists = []

        list0 = [value[0][i] for i in cell0index]
        list1 = [value[1][i] for i in cell1index]

        list0.sort(reverse=True)
        list1.sort(reverse=True)

        result0 = check_zero_majority(list0)
        result1 = check_zero_majority(list1)

        if result0 == True:
            continue
        if result1 == True:
            continue

        if len(list0) < len(list1):
            list1 = list1[:len(list0)]
        else:
            list0 = list0[:len(list1)]

        new_lists.append(list0)
        new_lists.append(list1)

        targetexp[key] = new_lists

    LRcorr = {}
    LRPcorr = {}
    for key, value in targetexp.items():
        x = value[0]
        y = value[1]
        if np.std(x) == 0 or np.std(y) == 0:
            continue
        correlation, p = pearsonr(x, y)
        LRPcorr[key] = [correlation, p]
        if p > 0.05:
            continue
        LRcorr[key] = [correlation, p]

    p2_LR_corr = FDRP(LRcorr)

    sorted_LRcorr = dict(sorted(p2_LR_corr.items(), key=itemgetter(1), reverse=True))

    return sorted_LRcorr


if __name__ == '__main__':

    f = open(expression_txt)
    cellname = []
    geneexp = {}
    flag = 0
    for p in f:
        flag += 1
        t = p.split()
        if flag == 1:
            cellname = t
            continue
        geneexp[t[0]] = [float(i) for i in t[1:]]
    f.close()

    L = []
    R = []
    LR = []
    f = open(LRDB_txt)
    flag = 0
    for p in f:
        flag += 1
        if flag == 1:
            continue
        t1 = p.split('\n')[0]
        LR.append(t1)
        tlr = t1.split(' - ')
        tlr[0] = tlr[0].split(' ')[0]
        l = []
        r = []
        l.append(tlr[0])
        L.append(l)
        if tlr[1][0] == '(':
            tr = tlr[1].split('+')
            r.append(tr[0][1:])
            r.append(tr[1][:-1])
            R.append(r)
        else:
            r.append(tlr[1])
            R.append(r)

    f.close()

    type_cell = {}
    f = open(labelname)
    for p in f:
        t = p.split('\t')
        t[1] = t[1].split('\n')[0]
        if t[1] not in type_cell.keys():
            type_cell[t[1]] = []
        type_cell[t[1]].append(t[0])

    f.close()

    celltype = list(type_cell.keys())
    typeindex = {}
    for i in celltype:
        cell0 = type_cell[i]
        cell0index = []
        for j in cell0:
            cell0index.append(cellname.index(j))
        typeindex[i] = cell0index

    gene2exp = {}
    for i in range(len(LR)):

        if len(L[i]) == 1 and L[i][0] not in geneexp.keys():
            continue

        if len(R[i]) == 1 and R[i][0] not in geneexp.keys():
            continue

        if len(R[i]) > 1:
            for r in R[i]:
                r0 = R[i][0]
                r1 = R[i][1]
            if r0 not in geneexp.keys():
                continue
            if r1 not in geneexp.keys():
                continue

        if len(L[i]) == 1:
            lexp = geneexp[L[i][0]]
            gene2exp[L[i][0]] = geneexp[L[i][0]]

        if len(R[i]) == 1:
            rexp = geneexp[R[i][0]]
            gene2exp[R[i][0]] = geneexp[R[i][0]]

        if len(R[i]) > 1:
            for r in R[i]:
                r0 = R[i][0]
                r1 = R[i][1]
            gene2exp[r0] = geneexp[r0]
            gene2exp[r1] = geneexp[r1]

    geneexp2 = {}
    for k in gene2exp.keys():
        n = 0
        exp = gene2exp[k]
        kexp = []
        for tt in celltype:
            for t in typeindex[tt]:
                kexp.append(exp[t])
            count = sum(1 for num in kexp if num > 0)
            # if k=='Col4a2':
            #    print(n)
            # count2 = sum(1 for num in kexp if num == 0)

            if count > len(kexp) * 1 / 2:
                n += 1
                # if k=='Col4a2':
            #   print(count,len(kexp),len(kexp)*1 / 2,n)
        if n >= len(type_cell) / 2:
            continue
        geneexp2[k] = gene2exp[k]

    LRexp = {}
    for i in range(len(LR)):
        if len(L[i]) == 1 and L[i][0] not in geneexp2.keys():
            continue
        if len(R[i]) == 1 and R[i][0] not in geneexp2.keys():
            continue
        if len(R[i]) > 1:
            for r in R[i]:
                r0 = R[i][0]
                r1 = R[i][1]
            if r0 not in geneexp2.keys():
                continue
            if r1 not in geneexp2.keys():
                continue
        LRexp[LR[i]] = []
        if len(L[i]) == 1:
            lexp = geneexp2[L[i][0]]
        if len(R[i]) == 1:
            rexp = geneexp2[R[i][0]]
        if len(R[i]) > 1:
            for r in R[i]:
                r0 = R[i][0]
                r1 = R[i][1]
            r01 = geneexp2[r0] + geneexp2[r1]
            rexp = [x / 2 for x in r01]
        LRexp[LR[i]].append(lexp)
        LRexp[LR[i]].append(rexp)

    target = []
    target.append(Ltypename)
    target.append(Rtypename)
    sorted_LRcorr = SortLR(target, type_cell, cellname, LRexp)

    print('Saving...')

    fw = open(fold + os.sep + 'target_LR.txt', 'w')
    fw.write('L_R_pair' + '\t' + 'pearson_corr' + '\t' + 'P_value' + '\n')
    for k in sorted_LRcorr.keys():
        fw.write(str(k) + '\t')
        kvalue = sorted_LRcorr[k]
        fw.write(str(kvalue[0]) + '\t' + str(kvalue[1]) + '\n')
    fw.close()

    print('Done')


