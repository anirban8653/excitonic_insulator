import numpy as np
import matplotlib.pyplot as plt

cdw = -0.65
eta = 0.003
ny = 80  # y length
nx = 6
mu = -0.23
lyr = ny  # number of layers for the topograph
ival = len(np.array([i for i in range(int(nx * ny * 4/2 ))]))
nw = 501

save_path = f'D:/phd_projects/NTU/exciton_NTU/6-fold-DOS/separate_cosine_terms/data_xedge_eta_{eta}/'

en, ldos = np.loadtxt('xedge_sem_cdw-0.65_mu-0.23_wmax0.25Ny80Nx6eta0.003Ns4000Ldos_wte2.dat', unpack=True)

plt.rcParams["figure.figsize"] = [5, 6]
plt.rcParams["figure.autolayout"] = True
# plt.subplot(1,2,2)
# for i in range(0,48,2):

x = en[0::ival] * 1000

y_all_cell_spin_up = []
for k in range(0, ival, 4):
    y_cell_spin_up = sum([ldos[i::ival] for i in range(k + 0, k + 4)]) / 4
    # y_cell_spin_up = ldos[k+3::ival]
    y_all_cell_spin_up.append(y_cell_spin_up)
y_all_cell_spin_up = np.array(y_all_cell_spin_up)

for i in range(0, len(y_all_cell_spin_up)):
    x = np.append(x, y_all_cell_spin_up[i])
# print(len(y))
topo_data = x.reshape(int(nx * ny/2 ) + 1, nw).T
# print(topo_data)
np.savetxt(f"unitcell_avg_ldos_data_sm_2pi6_mu{mu}_eta_{eta}_cdw{cdw}.dat", topo_data, fmt="%12.6f")
# %%
data = np.loadtxt(f'unitcell_avg_ldos_data_sm_2pi6_mu{mu}_eta_{eta}_cdw{cdw}.dat')

for s in range(1,7):
    for i in range(s,int(nx*ny*2/4+1),6*24):
    # for i in range(1,7):
        plt.plot(data[:, 0], data[:, i],lw=2)
    # plt.axvline(x=0, color='black', lw=1, ls ='dotted')
    plt.hlines(0, -80,50 ,color='black', lw=1, ls ='dotted')
    # plt.ylim(0, 0.8)
    plt.xlim(-200,200)
    #    plt.legend(loc=(0.62,0.72),frameon=False,fontsize=15)
    plt.xlabel('$E-E_F$ (meV)', fontsize=15)
    plt.ylabel('DOS', fontsize=15)
    plt.tick_params(direction='in', length=7, width=1, colors='k', bottom=True,
                    top=True, left=True, right=True)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    # plt.xticks(np.array([i for i in range(-200,200,100)]),fontsize=15)     
    # plt.grid(True)
    plt.savefig(f"{save_path}xe_dos_OP_sm_mu{mu}_{s}.png", dpi=300)
    plt.show()
    dos_data = np.array([[data[:, 0][j], data[:, s][j], data[:, int(nx*ny*2/4)][j]] for j in range(len(data[:, 0])) ] )
    np.savetxt(f"{save_path}DOS_profile_data_x{s}.dat",dos_data, fmt="%12.6f")
