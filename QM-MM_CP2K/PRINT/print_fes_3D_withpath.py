import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from scipy.interpolate import griddata
import numpy as np

data = np.loadtxt('fes_24.7ps.dat')
path_x=[]; path_y=[]
with open("path.spath") as f:
        for linea in f:
                fields = linea.split()
                if len(fields) > 0 :
                        path_x.append(float(fields[1])); path_y.append(float(fields[2]))
f.close()


x=np.unique(data[:, 0])
y=np.unique(data[:, 1])
xy= np.stack([np.tile(x, y.shape[0]), np.repeat(y, x.shape[0])]).T

Z = griddata(data[:, :2], data[:, 2], xy, method='cubic').reshape((x.shape[0], y.shape[0]))

cmap="RdYlBu_r"
fig,ax= plt.subplots()
plt.contourf(x,y,Z, levels=np.arange(0,1.1,1), colors="navy", alpha=0.15)
plt.imshow(Z, cmap=cmap, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower', aspect='auto', vmin=1, vmax=round(max(data[:, 2])))
cbar=plt.colorbar(shrink=0.85, aspect=16, ticks=np.arange(0,max(data[:, 2])+2,5))
plt.contour(x,y,Z, levels=np.arange(0,round(max(data[:, 2]))+2,1), extend='both', colors="black", linewidths=0.50, alpha=0.50, linestyles='dashed')
cbar.ax.set_ylabel("$\Delta$G (kcal/mol)", rotation=270, fontsize=12, fontweight='bold', labelpad=28)

plt.plot(path_x, path_y, linewidth=2.0,color='black', linestyle=':')

plt.plot(-1.385, 4.689, marker="o", markersize=5, markeredgecolor="black", markerfacecolor="white")
plt.plot(-4.387,-1.088, marker="o", markersize=5, markeredgecolor="black", markerfacecolor="white")
plt.plot(-4.007, 0.390, marker="o", markersize=5, markeredgecolor="black", markerfacecolor="white")
plt.plot(-4.156, 2.942, marker="o", markersize=5, markeredgecolor="black", markerfacecolor="white")
plt.plot(-4.156, 2.271, marker="o", markersize=5, markeredgecolor="black", markerfacecolor="white")
plt.plot(-3.500, 3.503, marker="o", markersize=5, markeredgecolor="black", markerfacecolor="white")

ax.tick_params(bottom=False,labelbottom=False)
ax.tick_params(left=False,labelleft=False)
output="surface.svg"
plt.savefig(output, dpi=300, transparent=True)
plt.show()
plt.close()
