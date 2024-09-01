import matplotlib.pyplot as plt
import math as math
import numpy as np

time=[]; dist1=[]
i=0
with open("path.spath") as f:
        for linea in f:
                fields = linea.split()
                if fields[0][0] != '#':
                        time.append(i); dist1.append(float(fields[3])); i=i+1
f.close()

xlabel="Time (ps)"
ylabel="CV ($\AA$)"
output = "path.svg"

fig,ax= plt.subplots()

ax.tick_params(bottom=False,labelbottom=False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.plot(time, dist1, linewidth=1.5,color='black')
ax.hlines(y=4.46149, xmin=-10, xmax=0, linewidth=1.0,color='black',linestyle='dashed')
ax.hlines(y=0.00, xmin=105, xmax=115, linewidth=1.0,color='black',linestyle='dashed')
plt.ylim(-0.1,17)
plt.xlim(-10,115)
plt.savefig(output, transparent=True)
plt.show()
plt.close()
