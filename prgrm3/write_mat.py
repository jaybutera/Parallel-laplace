import numpy as np

n = 10

with open('mat', 'wb') as f:
    # Generate random nxn matrix
    A = np.random.rand(n,n)

    # Make diagonal 0
    for i in range(n):
        A[i][i] = 0.

    # Write file as binary
    A.tofile(f)

    print A

'''
with open('mat', 'rb') as f:
    A = np.fromfile(f, dtype=np.float64)

    print A
'''
