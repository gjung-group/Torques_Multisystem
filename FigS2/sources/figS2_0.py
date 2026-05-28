
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

params = {
        'font.family' : 'Times'
        }
matplotlib.rcParams.update(params)
plt.rcParams['text.usetex'] = True

datadir2write = "/Users/yjpark/YJ/Projects/Multisystem/Figures/figS2/sources/"

# fully rigid & fully relaxed ; corners are fixed ; 
data1   = np.genfromtxt(f"{datadir2write}energycurve_GBN.dat", 
                        skip_header=1) # (ang12, ang32, Natom, E_initail(eV), E_final(eV))

# bottom-hBN rigid   ; with corners are fixed
data2   = np.genfromtxt(f"{datadir2write}energycurve_GBNrigid.dat", 
                        skip_header=1) # (ang12, ang32, Natom, E_initail(eV), E_final(eV))

# To plot |a_1^{(3)}| ; top graphene layer's lattice constant
data4alat3  = np.genfromtxt(f"{datadir2write}Hermann_Indices.txt", 
                        skip_header=2, usecols=(1, 8))  # (ang32, alat^{(3)})

# To estimate the corrections in the number of top-layer atoms and their energies
data2adjust   = np.genfromtxt(f"{datadir2write}Hermann_Indices.txt", 
                        skip_header=2, usecols=(1, 6,7, 10))  # (ang32, a", b", lambda_multiple)

# Set reference as the double-moire commensurate condition
angref12, angref32 =  1.202325 , 0.579824  # Reference for double-moire commensuration
print(f"Reference angle : {data2adjust[data2adjust[:,0]==angref32, 0]} degree.")

# --- Number-of-atom correction
Natom_L3          = 2*(data2adjust[:,1]**2 + data2adjust[:,1]*data2adjust[:,2] + data2adjust[:,2]**2)
Natom_L3_ref      = Natom_L3[data2adjust[:,0]==angref32] 
lambda_ratio      = data2adjust[:,-1]  # ratio between the simulation cell length with respect to their reference (double-moire commensurate)
Natom_corr        = (lambda_ratio**2) * Natom_L3_ref - Natom_L3
print(f"Natom_corr  =  (lambda/lambda^comm.)^2 * Natom_L3_ref - Natom_L3 ", f"ref angle = {data2adjust[data2adjust[:,0]==angref32, 0]}")

# --- Energy correction
Emin_graphene      = -7.3949972515997 ; print(f"Emin_graphene = {Emin_graphene:16.12f} eV/atom at 2.46019 \\AA")
# Emin_hBN           = -6.68997268334775; print(f"Emin_hBN      = {Emin_hBN:16.12f} eV/atom at 2.505759 \\AA")
E_correct         = Natom_corr  * Emin_graphene 

print(f"Natom_L3_correct", Natom_corr)
print(f"E_correct = E_ref^3 * Natom_L3_correct ", E_correct)

ang32s_rigid   , etots_rigid     = data1[:,1] , (data1[:,-2]+E_correct)/(data1[:,-3]+Natom_corr)
ang32s_rigidhBN, etots_rigidhBN  = data2[:,1] , (data2[:,-1]+E_correct)/(data2[:,-3]+Natom_corr)
ang32s_relax   , etots_relax     = data1[:,1] , (data1[:,-1]+E_correct)/(data1[:,-3]+Natom_corr)

eref1 = etots_relax[0]

fig      = plt.figure(figsize=(6,6))
msize, fsize = 4, 16
ax1      = fig.add_axes([0.2, 0.58, 0.62, 0.35]) # for figsize=(6,6) two figures up and down
ax2      = fig.add_axes([0.2, 0.1 , 0.62, 0.35])

for datas, mcolor in zip([
    [ang32s_rigid   , etots_rigid    - eref1],
    [ang32s_rigidhBN, etots_rigidhBN - eref1],
    [ang32s_relax   , etots_relax    - eref1]]
    ["C0", "C1", "C2"]):

    ax1.plot(datas[    0], datas[    1]*1000, 'o-' , lw=0.8, color=mcolor,ms=msize)#, alpha=0.3)

    print('Emin (meV/atom)', np.min(datas[    1]*1000), 'at ', datas[0][np.argmin(datas[    1]*1000)])
 
ax2.plot(data4alat3[:,0], data4alat3[:,1], 'ko-', ms=msize)

for ang_for_a3 in [-0.010739, 0.579824, 0.973401, 1.288158]:
    ax1.axvline(x=ang_for_a3, ls='-', color=f'gray', alpha=0.5)
    ax2.axvline(x=ang_for_a3, ls='-', color=f'gray', alpha=0.5)

for a3 in np.linspace(2.4594, 2.4610, 5, endpoint=True):
    ax2.axhline(y=a3, ls='-', color=f'gray', alpha=0.5)

    
ax1.axhline(y=0, ls='-', color='gray', alpha=0.5)

ax1.set_xlabel("$\\theta_{32}\\ (^{\\circ}$)", fontsize=fsize, fontname='times')
ax1.set_xlim(-0.05,1.55); ax1.set_xticks(np.linspace(0.0, 1.5,4, endpoint=True))
ax1.set_ylim(-0.25,2); ax1.set_yticks(np.linspace(0,2, 5, endpoint=True))

ax2.set_xlabel("$\\theta_{32}\\ (^{\\circ}$)", fontsize=fsize, fontname='times')
ax2.set_xlim(-0.05,1.55); ax1.set_xticks(np.linspace(0.0, 1.5,4, endpoint=True))
ax2.set_ylim(2.4592, 2.4612); ax2.set_yticks(np.linspace(2.4594, 2.4610, 5, endpoint=True))

ax1.tick_params(labelsize=fsize) 
ax1.set_ylabel("$E_{\\rm tot}$ (meV/atom)", fontsize=fsize, fontname='times')

ax2.tick_params(labelsize=fsize) 
ax2.set_ylabel("$|{\\bf a}_1^{(3)}|\\ \\AA$", fontsize=fsize, fontname='times')

fig.savefig(f"{datadir2write}figS2_0.pdf")

plt.rcParams['text.usetex'] = False