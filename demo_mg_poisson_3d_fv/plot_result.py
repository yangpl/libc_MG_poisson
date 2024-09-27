import matplotlib.pyplot as plt
import numpy as np


y1, y2= np.loadtxt('result.txt', skiprows=1, unpack=True) #skip 0 lines
x = range(y1.size)

plt.plot(x, y1, label='x_true')
plt.plot(x, y2, label='x_rec')

plt.legend()
plt.show()
