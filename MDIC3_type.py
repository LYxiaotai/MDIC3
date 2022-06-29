import numpy as np


f2 = open('cell_label.txt','r')
label = []
labels = set()
flag = 0
for p in f2:
    flag += 1
    if flag == 1:
        continue
    t = p.split('\t')
    tt = t[1].split('\n')
    label.append(tt[0])
    labels.add(str(tt[0]))
labelsl = list(labels)
llen = len(labelsl)
indexa = {}
for i in labelsl:
    indexa[i] = []
flag =0 
for i in label:
    flag +=1
    indexa[i].append(flag-1)
f2.close()


f1 = open('CCC_results.txt','r')
M = []
for p in f1:
    t = p.split()
    tt = list(map(float,t))
    M.append(tt)
M = np.array(M)
f1.close()


N = np.zeros((llen,llen))
for l1 in labelsl:
    I = indexa[l1]
    for l2 in labelsl:
        J = indexa[l2]
        if l1 == l2:
            for i in I:
                for j in J:
                    N[labelsl.index(l1)][labelsl.index(l2)]+=abs(M[i][j])
        if l1 != l2:
            for i in I:
                for j in J:
                    if M[i][j]>=0:
                        N[labelsl.index(l1)][labelsl.index(l2)]+=M[i][j]
                    if M[i][j]<0:
                        N[labelsl.index(l2)][labelsl.index(l1)]+=abs(M[i][j])


fw = open('type_results.txt','a')
for tt in range(len(labelsl)):
    if tt <len(labelsl)-1:
        fw.write(labelsl[tt]+'\t')
    if tt == len(labelsl)-1:
        fw.write(labelsl[tt]+'\n')
for i in range(len(N)):
    Ni = list(N[i])
    fw.write(str(labelsl[i])+'\t')
    for ni in range(len(Ni)):
        if ni <len(Ni)-1:
            fw.write(str(Ni[ni])+'\t')
        if ni == len(Ni)-1:
            fw.write(str(Ni[ni])+'\n')
fw.close()



