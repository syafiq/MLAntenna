import scipy.io

f = scipy.io.loadmat('data.mat')
data = f.get('Vin')
