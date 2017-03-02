import sys
import math
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

def ThomasAlgorithm(a, b, c, d, n):
     c_prime = [0]*(n-1)
     d_prime = [0]*(n-1)
     y = [0]*(n-1)

     c_prime[0] = c[0] / b[0]
     d_prime[0] = d[0] / b[0]

     for i in xrange(1, n-1):
        c_prime[i] = c[i] / (b[i] - a[i] * c_prime[i-1])
        d_prime[i] = (d[i] - a[i]*d_prime[i-1]) / (b[i] - a[i] * c_prime[i-1])
     
     y[n-2] = d_prime[n-2]

     for i in reversed(xrange(n-2)):
        y[i] = d_prime[i] - c_prime[i] * y[i+1]
     return y

def TridiagonalMatrix(x0, y0, xn, yn, h, n):

    a = [ai(x0,h,i) for i in xrange(1, n)]
    b = [bi(x0,h,i) for i in xrange(1, n)]
    c = [ci(x0,h,i) for i in xrange(1, n)]
    d = [di(x0,h,i) for i in xrange(1, n)]

    return ThomasAlgorithm(a, b, c, d, n)   

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
