import numpy as np
data=np.loadtxt('tmp')
npoints=2000
x=data[0:npoints,0]
y=data[0:npoints,1]
z=data[0:npoints,2]
tri=data[npoints:]

fig, ax = plt.subplots()
ax.set_xlim((-1.5, 1.5))
ax.set_ylim((-1.5, 1.5))
ax.plot(x,y,'.')
#for i in [0]:
for i in range(len(tri)):
#for i in range(10):
#    print i
    cl='black'
 #   if i==21:
 #       cl='green'
    ax.plot([x[int(tri[i][0])],x[int(tri[i][1])]],[y[int(tri[i][0])],y[int(tri[i][1])]],color=cl)
    ax.plot([x[int(tri[i][1])],x[int(tri[i][2])]],[y[int(tri[i][1])],y[int(tri[i][2])]],color=cl)
    ax.plot([x[int(tri[i][2])],x[int(tri[i][0])]],[y[int(tri[i][2])],y[int(tri[i][0])]],color=cl)

#circ=plt.Circle((-0.706149, -0.489123), 0.479517, color='red')
#ax.add_artist(circ)
#circ=plt.Circle((0.169623,-0.023845), 0.523879, color='orange')
#ax.add_artist(circ)
