import numpy as np
import matplotlib.pyplot as plt

grid = np.genfromtxt('output/fourier_transform.csv', delimiter=',')

fig, ax = plt.subplots()
print(np.min(grid))
pos = ax.imshow(np.log(grid - np.min(grid)))
ax.title.set_text("grid")
fig.colorbar(pos, ax=ax)

plt.show()
