from matplotlib import pyplot
from pylab import genfromtxt  


mat0 = genfromtxt("iterate.txt")
pyplot.plot(mat0[:,0], mat0[:,1], label = "Gauss-Seidel")
pyplot.yscale('log')
pyplot.grid(True) #grid(color='r', linestyle='-', linewidth=1)
pyplot.legend()
pyplot.title('Multigrid: F cycle')
pyplot.show()
