import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
point=np.loadtxt('tmp')
npo=2000
x=point[0:npo,0]
y=point[0:npo,1]
z=point[0:npo,2]
tri = point[npo:]

#normal = np.array([0.504135, 0.211107, 2.611719])
## a plane is a*x+b*y+c*z+d=0
## [a,b,c] is the normal. Thus, we have to calculate
## d and we're set
#d = -(normal[0]*point[7][0]+normal[1]*point[7][1]+normal[2]*point[7][2])
## create x,y
#xx, yy = np.meshgrid(np.linspace(-1,1,10), np.linspace(-1,1,10))
## calculate corresponding z
#z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]
## plot the surface
plt3d = plt.figure().gca(projection='3d')
#plt3d.plot_surface(xx, yy, z, alpha=0.2)

ax = plt.gca()
ax.scatter(x,y,z, color='green')
#ax.scatter(x[2],y[2],z[2], color='red')
#ax.scatter([-0.287578],[-0.129665],[0.099514], color='red',alpha=1)
#ax.scatter([-0.923377],[0.647341],[1.271675] , color='red',alpha=1)
#ax.scatter([-0.914716],[-0.057341],[0.839994], color='red',alpha=1)
#ax.scatter([-0.476860],[-0.037487],[0.228801], color='red',alpha=1)
#ax.scatter([-0.923975],[0.360424],[0.983636] ,color='red',alpha=1)
#ax.scatter([-0.387279],[-0.051990],[0.152688], color='red',alpha=1)
#ax.scatter([-0.495078],[0.325809],[0.351254] ,color='red',alpha=1)
#ax.scatter([-0.387279],[-0.051990],[0.152688] , color='red',alpha=1)
#ax.scatter([-0.495078],[0.325809],[0.351254]  , color='red',alpha=1)
#ax.scatter([-0.914716],[ -0.057341],[ 0.839994], color='red',alpha=1)


for i in range(len(tri)):
#for i in [5,22,23]:
    xp=[x[int(tri[i][0])],x[int(tri[i][1])],x[int(tri[i][2])]]
    yp=[y[int(tri[i][0])],y[int(tri[i][1])],y[int(tri[i][2])]]
    zp=[z[int(tri[i][0])],z[int(tri[i][1])],z[int(tri[i][2])]]
    verts=[zip(xp,yp,zp)]
    collection = Poly3DCollection(verts, linewidths=1, alpha=a)
    collection.set_facecolor([0.5,0.5,1])
    collection.set_edgecolor([0.,0.,0.])
    ax.add_collection3d(collection)

#for i in [-2,-1]:
#    xp=[x[int(tri[i][0])],x[int(tri[i][1])],x[int(tri[i][2])]]
#    yp=[y[int(tri[i][0])],y[int(tri[i][1])],y[int(tri[i][2])]]
#    zp=[z[int(tri[i][0])],z[int(tri[i][1])],z[int(tri[i][2])]]
#    verts=[zip(xp,yp,zp)]
#    collection = Poly3DCollection(verts, linewidths=1, alpha=a)
#    collection.set_facecolor([0.4,0.4,0.4])
#    collection.set_edgecolor([0.,0,0.])
#    ax.add_collection3d(collection)

#for i in [21]:
#    xp=[x[int(tri[i][0])],x[int(tri[i][1])],x[int(tri[i][2])]]
#    yp=[y[int(tri[i][0])],y[int(tri[i][1])],y[int(tri[i][2])]]
#    zp=[z[int(tri[i][0])],z[int(tri[i][1])],z[int(tri[i][2])]]
#    verts=[zip(xp,yp,zp)]
#    collection = Poly3DCollection(verts, linewidths=1, alpha=a)
#    collection.set_facecolor([1.,0.,0.2])
#    collection.set_edgecolor([0.,0.,0.])
#    ax.add_collection3d(collection)


i=11
#vec=np.loadtxt('tmp2')
#for i in range(len(tri)):
#    l = 10*sqrt(vec[i][0]*vec[i][0]+vec[i][1]*vec[i][1]+vec[i][2]*vec[i][2])
#    vec[i][0]=-vec[i][0]/l
#    vec[i][1]=-vec[i][1]/l
#    vec[i][2]=-vec[i][2]/l

#vec=[0.055159, 0.555867, 0.175401]

#for i in range(len(tri)):
#    soa = np.array([[x[int(tri[i][0])],y[int(tri[i][0])],z[int(tri[i][0])],vec[0],vec[1],vec[2]],\
#                    [x[int(tri[i][1])],y[int(tri[i][1])],z[int(tri[i][1])],vec[0],vec[1],vec[2]],\
#                    [x[int(tri[i][2])],y[int(tri[i][2])],z[int(tri[i][2])],vec[0],vec[1],vec[2]]])
#
#    X, Y, Z, U, V, W = zip(*soa)
#    ax.quiver(X, Y, Z, U, V, W)



plt.show()
