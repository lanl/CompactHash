from matplotlib import pyplot
import string
from array import *

fp = open("memopt.dat",'r')

X1  = array ('d', [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ])
Y1  = array ('d', [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ])
Y2  = array ('d', [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ])
Y8  = array ('d', [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ])
Y9  = array ('d', [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ])
Y10 = array ('d', [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ])
Y11 = array ('d', [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ])

line = fp.readline()

for i in range(0,7):
   line = fp.readline()
   value = line.split(",")
   X1[i]  = float(value[0])
   Y1[i]  = float(value[1])
   Y2[i]  = float(value[2])
   Y8[i]  = float(value[8])
   Y9[i]  = float(value[9])
   Y10[i] = float(value[10])
   Y11[i] = float(value[11])
   # Convert from secs to milliseconds
   Y8[i]  *= 1000
   Y9[i]  *= 1000
   Y10[i]  *= 1000
   Y11[i]  *= 1000
   # Normalize by number of cells
   Y8[i]  /= Y1[i]
   Y9[i]  /= Y1[i]
   Y10[i] /= Y1[i]
   Y11[i] /= Y1[i]
   #print (X1[i]*pow(2,i+1) * X1[i]*pow(2,i+1))/Y1[i]

fig = pyplot.figure()

fig.patch.set_facecolor('white')
pyplot.xlabel("Compressibility",fontdict={'fontsize':16})
pyplot.ylabel("Runtime/numcells (millisecs/cell)",fontdict={'fontsize':16})

pyplot.semilogy()

pyplot.plot(Y2, Y8,  color = 'green',   linestyle='-', linewidth=1.0, marker='D', label='Full Perfect on GPU')
pyplot.plot(Y2, Y9,  color = 'orange',  linestyle='-', linewidth=1.0, marker='D', label='7 write, 1 read on GPU')
pyplot.plot(Y2, Y10, color = 'blue',    linestyle='-', linewidth=1.0, marker='D', label='1 write, 3 read on GPU')
pyplot.plot(Y2, Y11, color = 'red',     linestyle='-', linewidth=1.0, marker='D', label='Compact GPU')

pyplot.grid()
pyplot.legend(loc = 2)
pyplot.suptitle("Comparison of Neighbor Finding Algorithms on the GPU",fontdict={'fontsize':16})

pyplot.savefig('neighbors_gpu.pdf',format='pdf')
pyplot.show()
