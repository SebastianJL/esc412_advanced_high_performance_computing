import numpy as np
import matplotlib.pyplot as plt
my_data_1 = np.genfromtxt('mass_grid_ngp_2d.csv', delimiter=',')
my_data_2 = np.genfromtxt('mass_grid_2nd_2d.csv', delimiter=',')
my_data_3 = np.genfromtxt('mass_grid_3rd_2d.csv', delimiter=',')
my_data_4 = np.genfromtxt('mass_grid_4th_2d.csv', delimiter=',')

my_data_1 = my_data_1[:, :-1]
my_data_2 = my_data_2[:, :-1]
my_data_3 = my_data_3[:, :-1]
my_data_4 = my_data_4[:, :-1]

plt.imsave('mass_grid_ngp_2d.png', np.log(my_data_1 + 1))  
plt.imsave('mass_grid_2nd_2d.png', np.log(my_data_2 + 1))  
plt.imsave('mass_grid_3rd_2d.png', np.log(my_data_3 + 1))  
plt.imsave('mass_grid_4th_2d.png', np.log(my_data_4 + 1))  
