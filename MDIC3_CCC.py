import os
import scipy.stats as stat
import math
import sys
import numpy as np
import numpy.linalg as la



f = open('gene_exp.txt')
AA = []
flag = 0
for p in f:
    t =p.split()
    flag += 1
    if flag == 1:
        continue
    a = t[1:]
    tt = list(map(float, a))
    AA.append(tt)
f.close()


AA = np.array(AA)
shapeA = np.shape(AA)

U1,S1,V1 = np.linalg.svd(AA)

d = min(shapeA)-len(S1)
dd = []
for i in range(d):
    dd.append(0)
S11 = np.append(S1,dd)

S11 = np.diag(S11)
if shapeA[0]<shapeA[1]:
    b = np.zeros((shapeA[0],shapeA[1]-shapeA[0]))
    S11 = np.column_stack((S11,b)) 
if shapeA[0]>shapeA[1]:
    b = np.zeros((shapeA[0]-shapeA[1],shapeA[1]))
    S11 = np.row_stack((S11,b))  



GRN = np.loadtxt('GRN.txt')
T = np.dot(GRN,S11)
Tp = np.linalg.pinv(T) 
V = np.dot(Tp,AA)   


print('Saving the calculation results......')
savename = 'CCC_result.txt'
np.savetxt(savename,V)

