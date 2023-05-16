# creating color map and manipulating color bar

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

interp = 'bilinear'
cdw = -0.65
eta = 0.003
mu = -0.23
save_path = f'D:/phd_projects/NTU/exciton_NTU/6-fold-DOS/separate_cosine_terms/data_xedge_eta_{eta}/'

data1 = np.loadtxt(f'unitcell_avg_ldos_data_sm_2pi6_mu{mu}_eta_{eta}_cdw{cdw}.dat')

nx = 6
ny = 40
layer = int(ny/2)
plotting_layer = int(ny/2)

# en1 = -50
# en2 = -20

# for i in range(len(data1)):
#     if data1[i, 0] == en1:
#         e1_num = i

# for j in range(len(data1)):
#     if data1[j, 0] == en2:
#         e2_num = j

en_line = -30
for j in range(len(data1)):
    if data1[j, 0] == en_line:
        eline_num = j


# print(e1_num,e2_num)


# define function for custom aspect ratio
def forceaspect(ax, aspect=1):
    im = ax.get_images()
    extent = im[0].get_extent()
    ax.set_aspect(abs((extent[1] - extent[0]) / (extent[3] - extent[2])) / aspect)


"""for a perticular energy (line cut)"""
data2 = data1[eline_num, 1:nx * layer + 1]


data = data2.reshape(layer, nx)
# print(data2)
new_data = []
for i in range(plotting_layer, 0, -1):
    for j in range(0, len(data[0])):
        new_data.append([j + 1, i, data[i - 1, j]])
info = 'X(a)     Y(b)     DOS'
np.savetxt(f"{save_path}Topographic_DOS_data_E_{en_line}.dat", new_data, header=info, fmt="%10.6f")
# colors = ['navy', 'blue', 'red', '#FFFF14']
# colors = ['navy', 'dodgerblue', 'red','white']
# cmap_name = 'my_list'
# cm = LinearSegmentedColormap.from_list(
#     cmap_name, colors, N=1000)
# norm = mpl.colors.Normalize(vmin=np.min(data), vmax=np.max(data))

plt.rcParams["figure.figsize"] = [12, 4]
plt.rcParams["figure.autolayout"] = True
fig = plt.figure()


ax1 = fig.add_subplot(131)

xtic = np.linspace(1,nx,nx)
ytic = np.round(np.array([i for i in range(1,plotting_layer,int(plotting_layer/5))]))

x, y, z = np.loadtxt(f'{save_path}Topographic_DOS_data_E_{en_line}.dat', skiprows=1, unpack=True)
N = int(nx)
z = z.reshape(plotting_layer, N)
ax1.set_xlabel('X(a)', fontsize=18)
ax1.set_ylabel('Y(b)', fontsize=18)
ax1.set_xticks(xtic,fontsize=18)
ax1.set_yticks(ytic,fontsize=18)
ax1.set_title(f"{en_line} meV", fontsize=15)
# plt.ylim(2, plotting_layer)
# plt.subplot(1,3,1)
im=ax1.imshow(z, interpolation=interp, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)))#,cmap=cm, norm=norm)
forceaspect(ax1, aspect=0.8)

ax1.tick_params(direction='out', length=5, width=1, colors='k', bottom=True, top=False, left=True, right=False)
fig.colorbar(im)
# plt.savefig(f"topo_ldos_sum_all_density_xe_OP_pi_sm_2pi6_mu{mu}.png", dpi=600)
# plt.savefig(f"topograph_DOS_en{en_line}.png", dpi=600)

# plt.show()

#%%
#x-edge-line-cut-slong-edge
plt.rcParams["figure.figsize"] = [5, 7]
plt.rcParams["figure.autolayout"] = True

ax2 = fig.add_subplot(132)
x_list = np.linspace(1,nx,nx)
ytic2 = np.linspace(0.03,0.1,int(5))
# plt.subplot(1,3,2)
ax2.plot(x_list,z[plotting_layer-1],'k')
ax2.set_xlabel('X(a)', fontsize=18)
ax2.set_ylabel('DOS', fontsize=18)
ax2.set_title("line cut of DOS along X edge\n", fontsize=15)
ax2.set_xlim(0,7)
# ax2.set_xticks(x_list,fontsize=18)
# ax2.set_yticks(ytic2,fontsize=18)
# forceaspect(ax2, aspect=0.4)

ax2.tick_params(direction='in', length=5, width=1, colors='k', bottom=True,
                top=False, left=True, right=False)
# plt.savefig(f"x_edge_DOS_line_cut_{en_line}.png",dpi = 300)
# plt.show()



#%%
#x-edge-line-cut-across-edge
# c_list = ['red','blue','green','yellow','pink','black']
# label_list = ['X = 1','X = 2','X = 3','X = 4','X = 5','X = 6']
# plt.rcParams["figure.figsize"] = [5, 5]
# plt.rcParams["figure.autolayout"] = True
y_list = np.linspace(1,plotting_layer,plotting_layer)
ytic3 = np.linspace(0.03,0.1,5)
xtic3 = np.round(np.linspace(1,plotting_layer,int(plotting_layer/5)))


ax3 = fig.add_subplot(133)
# plt.subplot(1,3,3)
for i in range(int(nx)):
    ax3.plot(y_list,np.flip(z.T[i]))#,c_list[i])#,label=label_list[i])
    ax3.legend(loc=(0.2,0.7),fontsize=12,ncol=2,frameon=False)
    
ax3.set_xlabel('Y(b)', fontsize=18)
ax3.set_ylabel('DOS', fontsize=18)
ax3.set_title("line cut of DOS across X edge\n", fontsize=15)
# ax3.set_xticks(xtic3,fontsize=18)
# ax3.set_yticks(ytic3,fontsize=18)


ax3.tick_params(direction='in', length=5, width=1, colors='k', bottom=True,
                top=False, left=True, right=False)
# plt.savefig(f"x_edge_DOS_across_line_cut_{en_line}.png",dpi = 300)
plt.savefig(f"{save_path}x_edge_{en_line}.png",dpi = 300)
plt.show()







line_cut = np.array([x_list,z[plotting_layer-1]]).T
np.savetxt(f"{save_path}x_edge_DOS_along_edge{en_line}.dat",line_cut, fmt="%10.6f")

line_cut = np.array([y_list,np.flip(z.T[0]),np.flip(z.T[1]),np.flip(z.T[2]),np.flip(z.T[3]),np.flip(z.T[4]),np.flip(z.T[5])]).T
np.savetxt(f"{save_path}x_edge_DOS_across_edge{en_line}.dat",line_cut, fmt="%10.6f")
#%%
# """CDW order parameter"""

# # define function for custom aspect ratio
# def forceaspect(ax, aspect=1):
#     im = ax.get_images()
#     extent = im[0].get_extent()
#     ax.set_aspect(abs((extent[1] - extent[0]) / (extent[3] - extent[2])) / aspect)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# # colors = ['navy', 'blue', 'red', '#FFFF14']
# colors = ['navy', 'dodgerblue', 'red','white']
# # colors = ['r','y']
# cmap_name = 'my_list'
# cm = LinearSegmentedColormap.from_list(
#     cmap_name, colors, N=1000)
# norm = mpl.colors.Normalize(vmin=np.min(data), vmax=np.max(data))

# # plt.rcParams["figure.figsize"] = [6, 5]
# # plt.rcParams["figure.autolayout"] = True

# # den_data = np.loadtxt("density_data_semi_metallic.dat")
# den_data = np.loadtxt("WTe2_cdw_gap_data.dat")
# dd = den_data.reshape(int(len(den_data)/4),4)
# x_list = np.linspace(1,nx,nx) 
# y_list = np.linspace(1,ny,ny) 
# density_grid = (sum(dd.T)/4).reshape(ny,nx)

# plt.imshow(density_grid, interpolation=interp,extent=(np.amin(x_list), np.amax(x_list), np.amin(y_list), np.amax(y_list)),
#            cmap=cm)
# forceaspect(ax, aspect=1)

# plt.xlabel('X', fontsize=18)
# plt.ylabel('Y', fontsize=18)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# plt.tick_params(direction='out', length=5, width=1, colors='k', bottom=True,
#                 top=False, left=True, right=False)
# plt.colorbar()
# plt.show()

# """CDW order parameter"""

# # define function for custom aspect ratio
# # def forceaspect(ax, aspect=1):
# #     im = ax.get_images()
# #     extent = im[0].get_extent()
# #     ax.set_aspect(abs((extent[1] - extent[0]) / (extent[3] - extent[2])) / aspect)
# # fig = plt.figure()
# # ax = fig.add_subplot(111)
# # colors = ['navy', 'blue', 'red', '#FFFF14']
# colors = ['navy', 'dodgerblue', 'red','white']
# # colors = ['r','y']
# cmap_name = 'my_list'
# cm = LinearSegmentedColormap.from_list(
#     cmap_name, colors, N=1000)
# norm = mpl.colors.Normalize(vmin=np.min(data), vmax=np.max(data))

# plt.rcParams["figure.figsize"] = [6, 5]
# plt.rcParams["figure.autolayout"] = True

# # den_data = np.loadtxt("density_data_semi_metallic.dat")
# den_data = np.loadtxt("WTe2_cdw_gap_data.dat")
# dd = den_data.reshape(int(len(den_data)/4),4)
# x_list = np.linspace(1,nx,nx) 
# y_list = np.linspace(1,ny,ny) 
# density_grid = (sum(dd.T)/4).reshape(ny,nx)

# plt.imshow(density_grid, interpolation=interp,extent=(np.amin(x_list), np.amax(x_list), np.amin(y_list), np.amax(y_list)),
#             cmap=cm)
# # forceaspect(ax, aspect=1)

# plt.xlabel('X', fontsize=18)
# plt.ylabel('Y', fontsize=18)
# plt.xlim(1,6)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# plt.tick_params(direction='out', length=5, width=1, colors='k', bottom=True,
#                 top=False, left=True, right=False)
# plt.colorbar()
# plt.show()


