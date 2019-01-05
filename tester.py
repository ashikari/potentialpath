import numpy.linalg as la
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt


def Gauss(x,y,x0,y0,sig):
	sigd = la.det(sig)
	sigi = la.inv(sig)
	mu = np.array([[x-x0],[y-y0]])
	prefactor = 1/math.sqrt(sigd*(2*math.pi)**2)
	return math.exp(-0.5*np.dot(np.dot(mu.transpose(),sigi),mu))


sig =[[1,0],[0,1]]


X = np.linspace(0,36.58,150)
Y = np.linspace(0,24.38,150)

X, Y = np.meshgrid(X, Y)
Z = np.zeros((len(X),len(Y)))


P =[0,1]

#hieght of potential field
for k in range(len(X)):
	for m in range(len(Y)):
		Z[k,m] = Gauss(X[k,m],Y[k,m],20,20,sig)


plt.contour(X,Y,Z,50)
plt.show()