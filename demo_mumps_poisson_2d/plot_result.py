import matplotlib.pyplot as plt
import numpy as np


y1, y2= np.loadtxt('result.txt', skiprows=1, unpack=True) #skip 0 lines
x = range(y1.size)
x_true = y1.reshape([129, 129])
x_rec= y2.reshape([129, 129])

plt.figure(figsize=(10, 4)) 
plt.subplot(121)
plt.imshow(x_true)
plt.colorbar()

plt.subplot(122)
plt.imshow(x_rec-x_true)
plt.colorbar()

plt.tight_layout()
plt.legend()
plt.show()
