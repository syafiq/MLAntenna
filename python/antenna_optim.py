import scipy.io
import numpy as np
import scipy.linalg as la
import math

f = scipy.io.loadmat('data.mat')
Vin = f.get('Vin')
Z0 = np.matrix(f.get('Z0'))
Zid = f.get('Zid')
bas = f.get('bas')
etaZ = f.get('etaZ')
evh = f.get('evh')
k = f.get('k')
kvh = f.get('kvh')
meshp = f.get('meshp')
X0 = (Z0 - Z0.getH())/2j
Rr = ((Z0 + Z0.getH())/2)*10000 # idk why it has to be multiplied with 10000
X0 = X0.real
Rr = Rr.real
Rs = 0.001
D, V = la.eig(Rr) # https://stackoverflow.com/questions/51247998/numpy-equivalents-of-eig-and-eigs-functions-of-matlab
D[D<0] = 0;
Rr = V*D/V
Rd = Rr + Zid*Rs
Z = Rd + 1j*X0
index, RE = la.eig(Rr,Rd)
index = np.max(index) # python standard max would not work for numpy array
RE = np.max(RE)
N = len(Z)
V = np.zeros(N)
V[math.floor(N/4)-17] = 1
