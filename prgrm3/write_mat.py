import numpy as np

n = 10

with open('mat', 'wb') as f:
    # Generate random nxn matrix
    A = np.random.rand(n,n)

    # Write file as binary
    A.tofile(f)

    print A

'''
with open('mat', 'rb') as f:
    A = np.fromfile(f, dtype=np.float64)

    print A
'''
