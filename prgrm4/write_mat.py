import numpy as np

n = 10

with open('mp_mat', 'wb') as f:
    A = np.random.rand(n,n)

    # Write file as binary
    A.astype('float32').tofile(f)

    print A

'''
with open('mat', 'rb') as f:
    A = np.fromfile(f, dtype=np.float64)

    print A
'''
