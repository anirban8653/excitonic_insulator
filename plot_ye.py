# creating color map and manipulating color bar

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

cdw = -0.8
eta = 0.002
path = 'D:/phd_projects/NTU/exciton_NTU/6-fold-DOS/separate_cosine_terms/data_y_edge_eta_0.002/'
nx = 80
ny = 6
interp = 'bilinear'
mu = -0.23
data1 = np.loadtxt(f'unitcell_avg_ldos_data_sem_yedge_2pi6_{mu}_eta_{eta}_cdw{cdw}.dat')
layer = int(nx/4)
plotting_layer = int(nx/4)

# en1 = -50
# en2 = -20

# for i in range(len(data1)):
#     if data1[i, 0] == en1:
#         e1_num = i

# for j in range(len(data1)):
#     if data1[j, 0] == en2:
#         e2_num = j

en_line = 0
for j in range(len(data1)):
    if data1[j, 0] == en_line:
        eline_num = j


def forceaspect(ax, aspect=1):
    """define function for custom aspect ratio"""
    im = ax.get_images()
    extent = im[0].get_extent()
    ax.set_aspect(abs((extent[1] - extent[0]) / (extent[3] - extent[2])) / aspect)


"""creating customize color bar"""
colors = ['navy', 'c', 'r', '#FFFF14']
cmap_name = 'my_list'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=1000)
norm = mpl.colors.Normalize(vmin=0, vmax=0.004)
# %%

plt.rcParams["figure.figsize"] = [8, 4]
plt.rcParams["figure.autolayout"] = True

# data2 = sum([data1[i, 1:6 * layer + 1] for i in range(e1_num, e2_num)]) / (e2_num - e1_num)

"""for a perticular energy (line cut)"""
data2 = data1[eline_num, 1:ny * layer + 1]

data = data2.reshape(layer, ny)

new_data = []
for i in range(plotting_layer, 0, -1):
    for j in range(0, len(data[0])):
        new_data.append([j + 1, i, data[i-1, j]])

# print(new_data)
# output data file
info = 'unit_cell     layer   ldos'
np.savetxt(f"{path}Topographic_DOS_{en_line}meV.dat", new_data, header=info, fmt="%10.6f")

colors = ['navy', 'blue', 'red', '#FFFF14']
cmap_name = 'my_list'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=1000)
norm = mpl.colors.Normalize(vmin=np.min(data), vmax=np.max(data))

# plt.rcParams["figure.figsize"] = [6, 5]
# plt.rcParams["figure.autolayout"] = True
fig = plt.figure()
ax = fig.add_subplot(121)


x, y, z = np.loadtxt(f'{path}Topographic_DOS_{en_line}meV.dat', skiprows=1, unpack=True)
z = z.reshape(plotting_layer, ny).T
fliped_z = np.flip(z)
im=ax.pcolormesh(fliped_z,shading='gouraud')
ax.set_xlabel('X(a)', fontsize=14)
ax.set_ylabel('Y(b)', fontsize=14)
ax.set_title(f"Energy = {en_line} meV")
# ax.xticks(fontsize=18)
# ax.yticks(fontsize=18)
# plt.ylim(1, 6)
# plt.imshow(z, interpolation=interp, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
#             cmap=cm, norm=norm)
# forceaspect(ax, aspect=1)

ax.tick_params(direction='out', length=5, width=1, colors='k', bottom=True,
                top=False, left=True, right=False)
fig.colorbar(im)
# plt.savefig(f"topo_ldos_sum_all_density_xe_OP_pi_sem_yedge_2pi6_{mu}.png", dpi=600)
# plt.show()

ax2 = fig.add_subplot(122)

x_list = np.linspace(1,plotting_layer,plotting_layer)
dos_line_data = [x_list]
for i in range(ny):
    ax2.plot(x_list,fliped_z[i])
    dos_line_data.append(fliped_z[i])
dos_line_data=np.array(dos_line_data).T
ax2.set_xlabel("X (a)",fontsize=14)
ax2.set_ylabel("DOS",fontsize=14)
ax2.set_title("DOS along X direction")
plt.savefig(f"{path}topographic_DOS_{en_line}meV.png", dpi=600)
# np.savetxt(f"{path}topographic_DOS_{en_line}meV.dat", new_data, header=info, fmt="%10.6f")
np.savetxt(f"{path}DOS_line_data_{en_line}meV.dat", dos_line_data, fmt="%10.6f")
# plt.savefig(f"topo_ldos_plot_{en_line}.png", dpi=600)
plt.show()







