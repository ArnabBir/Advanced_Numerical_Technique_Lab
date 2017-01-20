import sys
import numpy as np 
import math
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches

# y'' - 2xy' - 2y = -4x
# Hence A(x) =  -2x ; B(x) = - 2 ; C(x) = -4x

def A(x):
	return -2 * x
def B(x):
	return -2
def C(x):
	return -4 * x

def get_a(x, h):
	return 1/(h**2) - A(x)/(2.0 * h)
def get_b(x, h):
	return - 2/(h**2) + B(x)
def get_c(x, h):
	return 1/(h**2) + A(x)/(2.0 * h)

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

def TridiagonalBVP(x0, xn, h, n):
	x = [(x0 + itr * h) for itr in xrange(1, n)]
	#print x[n-1]
	a = [get_a(itr, h) for itr in x]
	b = [get_b(itr, h) for itr in x]
	c = [get_c(itr, h) for itr in x]
	d = [C(itr) for itr in x]
	
	b[0] += 4 / (2*h + 3) * a[0]
	c[0] += (-1) / (2*h + 3) * a[0]

	b[n-2] += 4 / (3 - 4*h) * c[n-2]
	a[n-2] += (-1) / (3 - 4*h) * c[n-2]
	d[n-2] += (2*h) / (3 - 4*h) * c[n-2]

	return ThomasAlgorithm(a, b, c, d, n)

def main():
	# h = 0.1, 0.05, 0.01
	#h = [0.1, 0.05, 0.01]
	stepsizes = [0.1, 0.05, 0.025, 0.001]
	# Boundary conditions y(1) = 0, y(1.4) = 0.0566
	x0 = 0.0
	xn = 1.0

	for step in stepsizes:
		file = open("Resut_h_" + str(step) + ".txt", 'w')
		n = int(math.ceil((xn - x0) / step))
		x = [(x0 + step * itr) for itr in xrange(n+1)]
		y = [0 for itr in xrange(n+1)]
		y[1:n] = TridiagonalBVP(x0, xn, step, n)
		y[0] = 4 / (2*step + 3) * y[1] - 1 / (2*step + 3) * y[2]
		y[n] = (4 * y[n-1] - y[n-2] - 2*step) / (3 - 4*step)
		#print y
		file.write("Value of y(x) with respect to x:\n\n\tx\ty(x)\n\n")
		for i in xrange(n+1):
			file.write("\t" + str(x[i]) + "\t" + str(y[i]) + "\n")
		file.close()
		plt.plot(x, y,  label="h={0}".format(step))
		plt.xlabel('x')
		plt.ylabel('y(x)')
		plt.savefig('Plot_h_' + str(step) + '.png')
	#red_patch = mpatches.Patch(color='red', label='h = 0.01')
	#green_patch = mpatches.Patch(color='green', label='h = 0.05')
	#blue_patch = mpatches.Patch(color='blue', label='h = 0.1')
	#plt.legend(handles=[red_patch, blue_patch, green_patch])
	plt.legend(bbox_to_anchor=(2,2), bbox_transform=plt.gcf().transFigure)
	plt.legend(loc="upper left", bbox_to_anchor=[0, 1],
           ncol=2, shadow=True, title="Legend", fancybox=True)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, bbox_transform=plt.gcf().transFigure)
	plt.savefig('Plot_h.png')
	plt.show()
main()


