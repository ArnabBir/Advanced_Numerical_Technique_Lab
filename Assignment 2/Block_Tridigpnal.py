import sys
import numpy as np
import matplotlib.pyplot as plt

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

    a = 0
    b = 1
    
    h1 = 0.1
    h2 = 0.05
    h3 = 0.02
    '''
    x1 = np.linspace(a,b,(b - a)/h1+1)
    n = (b - a)/h1 + 1
    #x1 = x1[0:(b - a)/h1]
    w1 = bvp(a, b, h1)
    #print w1
    z1 = w1[[range(w1.shape[0])],[0],[0]]
    y1 = w1[[range(w1.shape[0])],[1],[0]]
    y1 = y1[0]
    z1 = z1[0]
    #print y1, n
    y1[n-1] = h1*(z1[n-2]+1) / 2 + y1[n-2]
    print y1
    
    x2 = np.linspace(a,b,(b - a)/h2 +1)
    #x2 = x2[0:(b - a)/h2]
    w2 = bvp(a, b, h2)
    y2 = w2[[range(w2.shape[0])],[1],[0]]
    y2 = y2[0]
    
    
    x3 = np.linspace(a,b,(b - a)/h3 +1)
    #x3 = x3[0:(b - a)/h3]
    w3 = bvp(a, b, h3)
    y3 = w3[[range(w3.shape[0])],[1],[0]]
    y3 = y3[0]
    #print len(x3), len(y3)
    
    plt.ylabel('y')
    plt.xlabel('x')
    
    p1, p2, p3 = plt.plot(x3, np.interp(x3, x1, y1), 
                          x3, np.interp(x3, x2, y2), x3, y3)
    
    plt.legend([p1, (p1, p2), (p1,p2,p3)], ["h = 0.1", "h =0.05", "h = 0.005"], loc =4)
   
    plt.show()
    '''

    h = [0.5, 0.1, 0.05, 0.01]
    for step in h:
    	n = int((b - a)/step) + 1
    	file = open("Resut_h_" + str(step) + ".txt", 'w')
    	x = np.linspace(a,b,(b - a)/step+1)
    	w = bvp(a, b, step)
    	y= w[[range(w.shape[0])],[1],[0]]
    	z = w[[range(w.shape[0])],[0],[0]]
    	y = y[0]
    	z = z[0]
    	y[n-1] = step*(z[n-2]+1) / 2 + y[n-2]
    	file.write("Value of y(x) with respect to x:\n\n\tx\ty(x)\n\n")
    	for i in xrange(n):
    		file.write("\t" + str(x[i]) + "\t" + str(y[i]) + "\n")
		plt.plot(x, y,  label="h={0}".format(step))
		plt.xlabel('x')
		plt.ylabel('y(x)')
		plt.legend(bbox_to_anchor=(2,2), bbox_transform=plt.gcf().transFigure)
		plt.legend(loc="upper left", bbox_to_anchor=[0, 1],
           ncol=2, shadow=True, title="Legend", fancybox=True)
		plt.savefig('Plot_h_' + str(step) + '.png')
		plt.clf()
		#file.close()


	plt.legend(bbox_to_anchor=(2,2), bbox_transform=plt.gcf().transFigure)
	plt.legend(loc="upper left", bbox_to_anchor=[0, 1],
           ncol=2, shadow=True, title="Legend", fancybox=True)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, bbox_transform=plt.gcf().transFigure)
	plt.savefig('Plot_h.png')
	plt.show()

main()