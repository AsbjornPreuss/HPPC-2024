import numpy as np
import matplotlib.pyplot as plt
x, y, z = np.loadtxt('sir_out.txt', skiprows=1, unpack=True)
plt.plot(x)
