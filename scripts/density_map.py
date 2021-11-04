import numpy as np
import matplotlib.pyplot as plt
mass_grid_modulo = np.genfromtxt('mass_grid_modulo_2d.csv', delimiter=',')
mass_grid_if_else = np.genfromtxt('mass_grid_if_else_2d.csv', delimiter=',')
mass_grid_margins = np.genfromtxt('mass_grid_with_margins_2d.csv', delimiter=',')

datas = [np.log(mass_grid_modulo + 1),
         np.log(mass_grid_if_else + 1),
         np.log(mass_grid_margins + 1)
         ]

for data in datas:
    print(f'{data.shape = }')

fig, axes = plt.subplots(1, 3, figsize=(15, 15))
titles = ['modulo', 'if else', 'with margins']

for ax, title, data in zip(axes, titles, datas):
    pos = ax.imshow(data)
    ax.title.set_text(title)
    fig.colorbar(pos, ax=ax, shrink=0.5)


# Diffs

diff_modulo_if_else = mass_grid_modulo - mass_grid_if_else
diff_modulo_margins = mass_grid_modulo - mass_grid_margins
diffs = [
    diff_modulo_if_else,
    diff_modulo_margins
]

print(f'{np.max(np.abs((diff_modulo_if_else))) = }')
print(f'{np.max(np.abs((diff_modulo_margins))) = }')

fig, axes = plt.subplots(1, 2, figsize=(15, 15))
titles = ['modulo - if-else', 'modulo - margins']

for ax, title, data in zip(axes, titles, diffs):
    pos = ax.imshow(data, cmap='PiYG')
    ax.title.set_text(title)
    fig.colorbar(pos, ax=ax, shrink=0.5)