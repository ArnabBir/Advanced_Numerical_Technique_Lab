import  numpy as np 
import math
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches

# x^2y'' + xy' = 1
#The equation can be rewritten as y'' + 1/x y' = 1/x^2
# Hence A(x) = 1/x ; B(x) = 0 ; C(x) = 1/x^2

def A(x):
	return 1 / x
def B(x):
	return 0
def C(x):
	return 1/x**2

def get_a(x, h):
	return 1/h**2 - A(x)/(2 * h)
def get_b(x, h):
	return - 2/h**2 + B(x)
def get_c(x, h):
	return 1/h**2 + A(x)/(2 * h)

def ThomasAlgo(a, b, c, d, n):
	 c_dash = np.zeros(n-1)
	 d_dash = np.zeros(n-1)
	 c_dash[0] = c[0] / b[0]
	 d_dash[0] = d[0] / b[0]
	 for itr in xrange(1, n-1):
	 	c_dash[itr] = c[itr] / (b[itr] - a[itr] * c_dash[itr-1])
	 	d_dash[itr] = d[itr] / (b[itr] - a[itr] * c_dash[itr-1])
	 y = np.zeros(n-1)
	 y[n-2] = d_dash[n-2]
	 for itr in reversed(xrange(n-2)):
	 	y[itr] = d_dash[itr] - c_dash[itr] * y[itr+1]
	 return y

def TridiagonalBVP(x0, y0, xn, yn, h, n):
	x = [(x0 + itr * h) for itr in xrange(1, n)]
	#print x
	a = [get_a(itr, h) for itr in x]
	b = [get_b(itr, h) for itr in x]
	c = [get_c(itr, h) for itr in x]
	d = np.zeros(n-1)
	d[0] = C(x[0]) - a[0]* y0
	for itr in xrange(1,n-2):
		d[itr] = C(x[itr])
	d[n-2] = C(x[n-2]) - c[n-2] * yn
	#print a, b, c, d
	return ThomasAlgo(a, b, c, d, n)

def main():
	# h = 0.1, 0.05, 0.01
	#h = [0.1, 0.05, 0.01]
	stepsizes = [0.1, 0.05, 0.01]
	# Boundary conditions y(1) = 0, y(1.4) = 0.0566
	x0 = 1.0
	xn = 1.4
	y0 = 0.0
	yn = 0.0566
	for step in stepsizes:
		n = int(math.ceil((xn - x0) / step))
		x = [(x0 + step * itr) for itr in xrange(n+1)]
		print x
		y = np.zeros(n+1)
		y[0] = y0
		y[1:n] = TridiagonalBVP(x0, y0, xn, yn, step, n)
		y[n] = yn
		print y
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