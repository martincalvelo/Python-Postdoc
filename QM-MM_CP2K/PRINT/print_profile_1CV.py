import matplotlib.pyplot as plt
import math as math
import numpy as np

def media_movil(a, n):
        a = np.array(a)
        out=[]
        for j in range(len(a)):
                out.append(np.mean(a[j : j + n]))
        return np.array(out)

time=[]; dist1=[]; dist2=[]
with open("fes.dat") as f:
        for linea in f:
                fields = linea.split()
                if fields[0][0] != '#':
                        time.append(float(fields[0])); dist1.append(float(fields[1]))
f.close()



xlabel="CV ($\AA$)"
ylabel="$\Delta$G (kcal/mol)"
output = "step1_profile.svg"

fig= plt.figure(figsize=(8,5))
ax= fig.add_subplot(111)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

dist1_average=media_movil(dist1, 10)

plt.plot(time, dist1_average, linewidth=2.5,color='black', alpha=0.75)
plt.ylabel(ylabel, fontsize=14)
plt.xlabel(xlabel, fontsize=14)
plt.yticks(np.arange(0, 25.1, 2.5))
plt.ylim(0,20)
plt.xlim(-7.5,7.5)
plt.savefig(output, transparent=True)
plt.show()
plt.close()
