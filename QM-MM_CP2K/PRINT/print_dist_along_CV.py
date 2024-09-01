import matplotlib.pyplot as plt
import math as math
import numpy as np

def media_movil(a, n):
        a = np.array(a)
        out=[]
        for j in range(len(a)):
                out.append(np.mean(a[j : j + n]))
        return np.array(out)

time=[]; dist1=[]; dist2=[]; dist3=[]; dist4=[]
with open("reaction_coordinate_dist.txt") as f:
        for linea in f:
                fields = linea.split()
                if fields[0][0] != '#':
                        time.append(float(fields[0])); dist1.append(float(fields[1])); dist2.append(float(fields[2])); dist3.append(float(fields[3])); dist4.append(float(fields[4]))
f.close()

dist1_print=media_movil(dist1, 75)
dist2_print=media_movil(dist2, 75)
dist3_print=media_movil(dist3, 75)
dist4_print=media_movil(dist4, 75)

xlabel="CV ($\AA$)"
ylabel="Distance ($\AA$)"
output = "DistancewithCV.svg"

fig,ax= plt.subplots()

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.vlines(x=0.92, ymin=0.5, ymax=3.6, linewidth=1,color='black', linestyle='dashed')
plt.plot(time, dist1_print, linewidth=2,color='orange')
plt.plot(time, dist2_print, linewidth=2,color='tab:pink')
plt.plot(time, dist3_print, linewidth=2,color='tab:brown')
plt.plot(time, dist4_print, linewidth=2,color='tab:blue')
plt.ylabel(ylabel, fontsize=14)
plt.xlabel(xlabel, fontsize=14)
plt.ylim(0.5,4)
plt.xlim(-4,4)
plt.savefig(output, transparent=True)
plt.show()
plt.close()
