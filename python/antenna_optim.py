import scipy.io
import numpy as np

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
Rs = 0.001
V, D = np.linalg.eig(Rr)
print(V)
