
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

params = {
        'font.family' : 'Times'
        }
matplotlib.rcParams.update(params)
plt.rcParams['text.usetex'] = True

datadir2write = "/Users/yjpark/YJ/Projects/Multisystem/Figures/figS2/sources/"

ensref  = -65354.6796285206/9244  # G/BN (1.202325, 0.579824) # fully relaxed

# fully relax ; corners are fixed ; 
data1   = np.genfromtxt(f"{datadir2write}energycurve_GBN.dat", 
                        skip_header=1) # (ang12, ang32, Natom, E_initail(eV), E_final(eV))

# hBN rigid   ; with corners are fixed
data2   = np.genfromtxt(f"{datadir2write}energycurve_GBNrigid.dat", 
                        skip_header=1)

# To plot |a_1^{(3)}| ; top graphene layer's lattice constant
data4   = np.genfromtxt(f"{datadir2write}Hermann_Indices.txt", 
                        skip_header=1, usecols=(1, 8))  # (ang32, alat^{(3)})

angref12, angref32 =  1.202325 , -0.010739


cond1 = np.logical_and( data1[:,0] == angref12 , data1[:,1] == angref32 ); print(data1[cond1,0:2])
cond2 = np.logical_and( data2[:,0] == angref12 , data2[:,1] == angref32 ); print(data2[cond2,0:2])

eref1 = (data1[cond1,-1]/data1[cond1,-3])[0]*1000; print(f'G/BN             case : {eref1/1000:18.10f} eV/atom; E_ref         = {eref1-ensref*1000:12.6f} meV/atom; {eref1-ensref*1000-(eref1-ensref*1000):12.6f} meV/atom')
eref2 = eref1                                    ; print(f'G/BN (rigid hBN) case : {eref2/1000:18.10f} eV/atom; E_ref^rigidBN = {eref2-ensref*1000:12.6f} meV/atom; {eref2-ensref*1000-(eref1-ensref*1000):12.6f} meV/atom')

fig      = plt.figure(figsize=(6,6))
msize, fsize = 4, 16
ax1      = fig.add_axes([0.2, 0.58, 0.62, 0.35]) # for figsize=(6,6) two figures up and down
ax2      = fig.add_axes([0.2, 0.1 , 0.62, 0.35])

for iangg, angg in enumerate([1.202325]):
    cond1 = data1[:,0] == angg 
    cond2 = data2[:,0] == angg 

    if iangg == 0:
        
        ax1.plot(data1[ cond1 ,1], (data1[ cond1 ,-2]/data1[ cond1 ,-3])*1000-eref1, 'o-' , label="G/BN    (rigid    )",  lw=0.8, color=f'C0',ms=msize)#, alpha=0.3)
        ax1.plot(data1[ cond1 ,1], (data1[ cond1 ,-1]/data1[ cond1 ,-3])*1000-eref1, 'o-' , label="G/BN    (suspended)",  lw=0.8, color=f'C2',ms=msize)#, alpha=0.3)
        ax1.plot(data2[ cond2 ,1], (data2[ cond2 ,-1]/data2[ cond2 ,-3])*1000-eref2, 'o-' , label="G/BN    (rigid hBN)",  lw=0.8, color=f'C1',ms=msize)#, alpha=0.3)
        
        ax2.plot(data4[:,0], data4[:,1], 'ko-', ms=msize)
    else:
        ax1.plot(data1[ cond1 ,1], (data1[ cond1 ,-2]/data1[ cond1 ,-3])*1000-eref1, 'o-' ,lw=0.8, color=f'C0',ms=msize)#, alpha=0.3)
        ax1.plot(data1[ cond1 ,1], (data1[ cond1 ,-1]/data1[ cond1 ,-3])*1000-eref1, 'o-' ,lw=0.8, color=f'C2',ms=msize)#, alpha=0.3)
        ax1.plot(data2[ cond2 ,1], (data2[ cond2 ,-1]/data2[ cond2 ,-3])*1000-eref2, 'o-' ,lw=0.8, color=f'C1',ms=msize)#, alpha=0.3)
        
        ax2.plot(data4[:,0], data4[:,1], 'ko-', ms=msize)

for ang_for_a3 in [-0.010739, 0.579824, 0.973401, 1.288158]:
    ax1.axvline(x=ang_for_a3, ls='-', color=f'gray', alpha=0.5)
    ax2.axvline(x=ang_for_a3, ls='-', color=f'gray', alpha=0.5)

for a3 in np.linspace(2.4594, 2.4610, 5, endpoint=True):
    ax2.axhline(y=a3, ls='-', color=f'gray', alpha=0.5)

    
ax1.axhline(y=0, ls='-', color='gray', alpha=0.5)

ax1.set_xlabel("$\\theta_{32}\\ (^{\\circ}$)", fontsize=fsize, fontname='times')
ax1.set_xlim(-0.05,1.55); ax1.set_xticks(np.linspace(0.0, 1.5,4, endpoint=True))

ax2.set_xlabel("$\\theta_{32}\\ (^{\\circ}$)", fontsize=fsize, fontname='times')
ax2.set_xlim(-0.05,1.55); ax1.set_xticks(np.linspace(0.0, 1.5,4, endpoint=True))
ax2.set_ylim(2.4592, 2.4612); ax2.set_yticks(np.linspace(2.4594, 2.4610, 5, endpoint=True))

ax1.set_ylim(-0.25,2); ax1.set_yticks(np.linspace(0,2, 5, endpoint=True))

# ax1.legend(fontsize=fsize-6, loc='lower center' , ncol=1, framealpha=1)#, title="No Coulomb")

ax1.tick_params(labelsize=fsize) 
ax1.set_ylabel("$E_{\\rm tot}$ (meV/atom)", fontsize=fsize, fontname='times')

ax2.tick_params(labelsize=fsize) 
ax2.set_ylabel("$|{\\bf a}_1^{(3)}|\\ \\AA$", fontsize=fsize, fontname='times')

fig.savefig(f"{datadir2write}figS2_0.pdf")

plt.rcParams['text.usetex'] = False