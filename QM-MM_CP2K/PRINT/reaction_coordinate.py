import numpy as np

diff1=[];dists=[]
with open("DISTS") as f:
        for linea in f:
                fields = linea.split()
                if fields[0] != '#!':
                        dists.append([float(fields[1]),float(fields[2]), float(fields[3]), float(fields[4])])
f.close()

with open("DIFFS") as f:
        for linea in f:
                fields = linea.split()
                if fields[0] != '#!':
                        diff1.append(float(fields[1]))
f.close()

todo=[]
for i in range(len(diff1)):
        todo.append([diff1[i], dists[i][0], dists[i][1], dists[i][2], dists[i][3]])

todo.sort()

out=open("reaction_coordinate_dist.txt", 'w')
for i in range(len(todo)):
        out.write("%f %f %f %f %f \n" % (todo[i][0], todo[i][1], todo[i][2], todo[i][3], todo[i][4]))
out.close()
