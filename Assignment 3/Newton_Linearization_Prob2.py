import sys
import math
import numpy as np
import matplotlib.pyplot as plt

Yk=[]

def initialFunc(x):
    return (0.5 - float(x)/np.pi)

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

    a = np.array([4 + 2*(Yk[i+1] - Yk[i-1]) for i in range(1, n)])
    b = np.array([4 + 2*(Yk[i-1] - Yk[i+1]) for i in range(1, n)])
    c = np.array([(-8*h**2*Yk[i] + 4*h**2 - 8) for i in range(1, n)])
    d = np.array([4*h**2*(-Yk[i]**2 + Yk[i] + 1) - (Yk[i+1] - Yk[i-1])**2 for i in range(1, n)])
    return ThomasAlgorithm(a, b, c, d, n)   

def main():
    h_values = [0.01]
    x0 = 0.0
    xn = np.pi
    y0 = 0.5
    yn = -0.5

    for h in h_values :
        n = n = int(np.pi/h)
        print n
        n = int(math.ceil((xn - x0) / h))
        print n
        for i in range(n+1):
            Yk.append(0)
        Yk[0] = y0
        for i in range(1,n,1):
            Yk[i] = initialFunc(x0+i*h)
        Yk[n] = yn
        #for i in range(1,50,1):
        for i in range(5):
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