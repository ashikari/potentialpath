from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import math
import numpy as np
import numpy.linalg as la


#height normalized gaussian function
def Gauss(x,y,x0,y0,sig):
	sigd = la.det(sig)
	sigi = la.inv(sig)
	mu = np.array([[x-x0],[y-y0]])
	prefactor = 1/math.sqrt(sigd*(2*math.pi)**2)
	return prefactor*math.exp(-0.5*np.dot(np.dot(mu.transpose(),sigi),mu))

def Gradient(x,y,x0,y0,sig):
	mu = np.array([[x-x0],[y-y0]])
	return -Gauss(x,y,x0,y0,sig)*np.dot(la.inv(sig),mu)

#generate potential field
def Potfunc(x,y,P,Dim,typ):
	i = 0
	out = 0
	for i in range(len(P)):
		sig = [[Dim[i][0], 0],[0, Dim[i][1]]]
		out = typ[i]*Gauss(x,y,P[i][0],P[i][1],sig) + out
		i+=1
	return out

def GField(x,y,P,Dim,typ):
	i = 0
	out = np.zeros((2,1))
	for i in range(len(P)):
		sig = [[Dim[i][0], 0],[0, Dim[i][1]]]
		sigi = la.det(sig)
		mu = np.array([[x-P[i][0]],[y-P[i][1]]])
		k = typ[i]*Gradient(x,y,P[i][0],P[i][1],sig)
		out = k + out
		i+=1
	return out

#where are the obstacles?
P = []
Dim = []
typ = []

#borders of stage!
L = 130
Th = 5
xmax = 36.58
ymax = 24.38
mag = 1

P.append([xmax/2,0])
Dim.append([L,Th])
typ.append(mag)
P.append([xmax/2,ymax])
Dim.append([L,Th])
typ.append(mag)
P.append([0,ymax/2])
Dim.append([Th,L])
typ.append(mag)
P.append([xmax,ymax/2])
Dim.append([Th,L])
typ.append(mag)




#goal!
goalx = 3
goaly = 17
P.append([goalx,goaly])
Dim.append([50,50])
typ.append(-3)


#plotting

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.linspace(0,xmax,100)
Y = np.linspace(0,ymax,100)

X, Y = np.meshgrid(X, Y)
Z = np.zeros((len(X),len(Y)))

#hieght of potential field
for k in range(len(X)):
	for m in range(len(Y)):
		Z[k,m] = Potfunc(X[k,m],Y[k,m],P,Dim,typ)

x = []
y = []
x.append(25)
y.append(5)
m = 0

J = []
while math.sqrt((x[m]-goalx)**2+(y[m]-goaly)**2)>=0.005:

	G = GField(x[m],y[m],P,Dim,typ)
	J.append(G/math.sqrt(G[0]**2+G[1]**2))

	x.append(x[m] - 0.005*np.asscalar(J[m][0]))
	y.append(y[m] - 0.005*np.asscalar(J[m][1]))
	
	# if Potfunc(x[m],y[m],P,Dim,typ)>Potfunc(x[m-1],y[m-1],P,Dim,typ):
	# 	print 'error'
	# 	print m


	m+=1

	if m >9999:
		break

plt.plot(x,y)
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.RdBu,linewidth=0, antialiased=False)
plt.show()

plt.contour(X,Y,Z,100)
plt.plot(x,y)
plt.show()

print(P)
print(Dim)
print(typ)


