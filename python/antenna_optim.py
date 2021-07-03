import scipy.io
import numpy as np

f = scipy.io.loadmat('data.mat')
Vin = f.get('Vin')
Z0 = f.get('Z0')
Zid = f.get('Zid')
bas = f.get('bas')
etaZ = f.get('etaZ')
evh = f.get('evh')
k = f.get('k')
kvh = f.get('kvh')
meshp = f.get('meshp')

print(Z0)
print("===========")
print(np.transpose(Z0))
