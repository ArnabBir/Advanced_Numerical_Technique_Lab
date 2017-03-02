import sys
import math
import numpy as np
import matplotlib.pyplot as plt

Yk=[]

def ai(x,h,i):
    return 1
def bi(x,h,i):
    return (-2-2*Yk[i]*h*h)
def ci(x,h,i):
    return 1
def di(x,h,i):
    return (2*h*h) + Yk[i]*Yk[i]*h*h + 2*Yk[i] - Yk[i+1] -Yk[i-1]

def initialFunc(x):
    return x * (1 - x)

def block_diagonal(a, b, c, d):
    n = len(d)
    b_ = np.zeros(a.shape)
    c_ = np.zeros(a.shape)
    d_ = np.zeros((n,2,1))

    c_[0] = np.linalg.inv(b[0]).dot(c[0]) 
    d_[0] = np.linalg.inv(b[0]).dot(d[0]) 
    
    for i in range(1,n):
        b_[i] = b[i] - a[i].dot(c_[i-1])
        c_[i] = np.linalg.inv(b_[i]).dot(c[i])
        d_[i] = np.linalg.inv(b_[i]).dot((d[i] - a[i].dot(d_[i-1])))
    
    for i in range(n-2, -1, -1):
        d_[i] = d_[i] - c_[i].dot(d_[i+1])

    return d_

def ThomasAlgorithm(a, b, c, d, n):
     c_dash = np.zeros(n-1)
     d_dash = np.zeros(n-1)
     c_dash[0] = c[0] / b[0]
     d_dash[0] = d[0] / b[0]
     for itr in xrange(1, n-1):
        c_dash[itr] = c[itr] / (b[itr] - a[itr] * c_dash[itr-1])
        d_dash[itr] = (d[itr] - a[itr]*d_dash[itr-1]) / (b[itr] - a[itr] * c_dash[itr-1])
     y = np.zeros(n-1)
     y[n-2] = d_dash[n-2]
     for itr in reversed(xrange(n-2)):
        y[itr] = d_dash[itr] - c_dash[itr] * y[itr+1]
     return y

def TridiagonalMatrix(x0, y0, xn, yn, h, n):

    a = [ai(x0,h,i) for i in xrange(1, n)]
    b = [bi(x0,h,i) for i in xrange(1, n)]
    c = [ci(x0,h,i) for i in xrange(1, n)]
    d = [di(x0,h,i) for i in xrange(1, n)]
    return ThomasAlgorithm(a, b, c, d, n)   

def bvp(l, r, h):
    n = int ((r-l) / (1.0*h))
    a = np.zeros((n-1,2,2))
    b = np.zeros((n-1,2,2))
    c = np.zeros((n-1,2,2))
    d = np.zeros((n-1,2,1))

    for i in range(n-1):
        # As i starts from 0, we define x = l + (i+1)*h
        x = l+(i+1)*h             
        a[i] = np.array(([[1.0 / h**2 - 2/h, 0.0],[-h/2, -1.0]]))
        #b[i] = np.array(([[- 2.0 / h**2 + 1.0, 1/h**2 + 2.0/h],[-h/2, 1.0]]))
        b[i] = np.array(([[- 2.0 / h**2 + 1.0, -6],[-h/2, 1.0]]))
        c[i] = np.array(([[1.0 / h**2 + 2/h, 0],[0, 0]]))
        d[i] = np.array(([[1],[0]]))

    d[0] = d[0] - a[0].dot(np.array(([[0.0], [0.0]])))
    #d[n-2] = d[n-2] -  c[n-2].dot(np.array(([[a2], [a4]])))
    d[n-2] = d[n-2] - np.array([[1/h**2 + 2/h], [0.0]])
    w = block_diagonal(a, b, c, d)
    a1 = 0.0
    a3 = 0.0
    a4 = d[n-2][1] + h/2 + h/2*d[n-2][0]
    a2 = 1.0

    w = np.vstack(([np.array(([[a1], [a3]]))], w))
    return np.vstack((w, [np.array(([[a2], a4]))]))

def main():
    h_values = [0.1]
    x0 = 0.0
    xn = 1.0
    y0 = 0.0
    yn = 0.0

    for h in h_values :

        n = int(math.ceil((xn - x0) / h))
        for i in range(n+1):
            Yk.append(0)

        Yk[0] = y0
        for i in range(1,n,1):
            Yk[i] = initialFunc(x0+i*h)
        Yk[n] = yn

        for i in range(1,50,1):
            dyk = TridiagonalMatrix(x0,y0,xn,yn,h,n)
            for j in range(1, n):
                Yk[j] = Yk[j] + dyk[j-1]
        

        x = [(x0 + h * i) for i in xrange(n+1)]

        file = open("result(h="+str(h)+").txt", 'w')
        for i in xrange(n+1):
            file.write(str((x0+i*h)) + "\t" + str(Yk[i]) + "\n")
        file.close()

        plt.plot(x, Yk,  label="h={0}".format(h))
        plt.xlabel('x')
        plt.ylabel('y(x)')
        plt.savefig('Plot(h=' + str(h) + ').png')
        plt.clf()

main()