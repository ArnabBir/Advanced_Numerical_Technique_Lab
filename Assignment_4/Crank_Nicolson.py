import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

def thomas_algorithm(a, b, c, d):

    a = list(a)
    b = list(b)
    c = list(c)
    d = list(d)
    assert len(a) == len(b) == len(c) == len(d)
    N = len(c)
    c_ = [None for i in range(N)]
    d_ = [None for i in range(N)]
    f = [None for i in range(N)]
    c_[0] = c[0]/b[0]
    d_[0] = d[0]/b[0]

    for i in range(1, N):
        if i is not N-1:
            c_[i] = c[i]/(b[i] - a[i]*c_[i-1])
        d_[i] = (d[i] - a[i]*d_[i-1])/(b[i] - a[i]*c_[i-1])

    f[N-1] = d_[N-1]
    for i in range(N-2, -1, -1):
        f[i] = d_[i] - c_[i]*f[i+1]

    return f

def crank_nikolson(v):
    dx = 0.2
    dt = 0.04
    rX = 1
    rT = 0.1
    n = int(rX/dx)
    m = int(rT/dt)
    X = np.linspace(0, rX, n)
    Y = np.linspace(0, rT, m)

    U = np.zeros((n, m), dtype=np.float64)
    U[0, :] = 1.0
    U[n-1, :] = 1.0
    U[:, 0] = 1.0
    for k in range(1, m/2):
        A = [None] + [v/dx**2 for i in range(1, n)]
        B = [-v/dx**2 - 2/dt] + [-v/dx**2 - 2/dt for i in range(1, n-1)] + [-v/dx**2 - 2/dt]
	#B = [-v/dx**2 - 2/dt for i in range(1, n+1)]        
	C = [v/dx**2 for i in range(1, n)] + [None]
        a = U[0:n-2, k-1]*(v/dx**2)
        b = U[1:n-1, k-1]*(v/dx**2)
        c = U[2:n, k-1]*(v/dx**2)
	B[0] = B[0] + 4.0 / 3.4 * A[0]
	C[0] = C[0] - 1.0 / 3.4 * A[0]
	B[n-1] = B[n-1] + 4.0 / 3.4 * C[n-1]  
	A[n-1] = A[n-1] - 1.0 / 3.4 * C[n-1]
        mid_vals =  a + b + c
	print A
        D = [U[0, k-1]*(-v/dx**2 - 2/dt + c/(2*dx))] + \
            list(mid_vals) \
                        + [U[n-1, k-1]*(-v/dx**2 - 2/dt - c/(2*dx))]
        U[:, k] = np.array(thomas_algorithm(A, B, C, D))[:, 0]

    return U, X, m

U, X, n = crank_nikolson(1)

plt.plot(X, U[:,n-1])
plt.show()
