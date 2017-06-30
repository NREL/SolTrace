import pylab as plt
import matplotlib.patches as patches
from math import *

#rectss = [line.replace("\n","").split(",") for line in open("meshpolar.txt",'r').readlines()[1:]]
rectss = [line.replace("\n","").split(",") for line in open("meshxy.txt",'r').readlines()[1:]]

rects = [[float(v) for v in L[1:]] for L in rectss]



ax = plt.gca()
#ax = plt.subplot(111, projection='polar')
#ax.set_theta_zero_location("N")

minx = 9e9
miny = 9e9
maxx = -9e9
maxy = -9e9

for i,R in enumerate(rects):
#    if len(rectss[i][0]) < 7:
#        continue
#    print (R[0],R[2]), R[1]-R[0], R[3]-R[2]
    xn,xx,yn,yx = R
    ax.add_patch(patches.Rectangle((R[0],R[2]), R[1]-R[0], R[3]-R[2], alpha=0.1 ) )
#    plt.text(R[0],R[2],rectss[i][0],fontsize='x-small')
    
    minx = min([minx,xn])
    maxx = max([maxx,xx])
    miny = min([miny,yn])
    maxy = max([maxy,yx])
    
#    break
    
plt.xlim(minx,maxx)
plt.ylim(miny,maxy)

#plt.xlim(-pi,pi)
#plt.ylim(-pi/2.,0.)