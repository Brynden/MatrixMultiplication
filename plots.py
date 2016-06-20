import numpy as np
import matplotlib.pyplot as plt


plt.bar(*np.loadtxt("test2.txt", unpack=True), linewidth=2.0)
plt.show()

