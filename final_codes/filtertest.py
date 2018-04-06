import numpy as np 
import matplotlib.pyplot as plt

el = np.arange(100000)
Fel = (1. - np.exp(-el**2/(2*1000.**2)))*np.exp(-el**2/(2*8000.**2))
plt.plot(el, Fel)
plt.show()