import numpy as np
import matplotlib.pyplot as plt

en, ldos = np.loadtxt('yedge_mu-0.25_cdw-0.08_wmax0.25Ny6Nx80eta0.002Ns4000Ldos_wte2_cdw.dat'
                      , unpack=True)

plt.rcParams["figure.figsize"] = [6, 5]
plt.rcParams["figure.autolayout"] = True

# plt.subplot(1,2,2)
# for i in range(0,48,2):
ny = 6  # y length
nx = 80 # x length
mu = 0.1

lyr = nx/2  # number of layers for the topograph
ival = int(4*ny*nx/2)
nw = 501
x = en[0::ival] * 1000

y_all_cell_spin_up = []
for k in range(0, ival, 4):
    # y_cell_spin_up = sum([ldos[i::ival] for i in range(k + 0, k + 4)]) / 4
    y_cell_spin_up = ldos[k+3::ival]
    y_all_cell_spin_up.append(y_cell_spin_up)
y_all_cell_spin_up = np.array(y_all_cell_spin_up)

for i in range(0, len(y_all_cell_spin_up)):
    x = np.append(x, y_all_cell_spin_up[i])
# print(len(y))
topo_data = x.reshape(int(nx * ny / 2) + 1, nw).T
# print(topo_data)
np.savetxt(f"unitcell_avg_ldos_data_sem_yedge_2pi6_{mu}.dat", topo_data, fmt="%12.6f")
# %%
data = np.loadtxt(f'unitcell_avg_ldos_data_sem_yedge_2pi6_{mu}.dat')
for i in range(1, int(nx * ny / 2) + 1, 6*24):
    plt.plot(data[:, 0], data[:, i], linewidth=2)
plt.axvline(x=0, color='black', linewidth=0.8)
# plt.axvline(x=-70, color='red', linewidth=0.8)
plt.ylim(0, 1.9)
# plt.xlim(-3,3)
#    plt.legend(loc=(0.62,0.72),frameon=False,fontsize=15)
plt.xlabel('$\omega (meV)$', fontsize=15)
plt.ylabel('DOS', fontsize=15)
plt.tick_params(direction='in', length=7, width=1, colors='k', bottom=True,
                top=True, left=True, right=True)
plt.yticks()
# plt.xticks([i for i in range(-240,280,40)])
# plt.grid()
plt.savefig(f"xe_dos_OP_sem_yedge{mu}.png",dpi = 300)
plt.show()

dos_data = np.array([[data[:, 0][i], data[:, 1][i], data[:, int(nx * ny / 2)][i]] for i in range(len(data[:, 0]))])
np.savetxt('energy_edge_bulk_DOS.dat',dos_data,fmt="%12.6f")
