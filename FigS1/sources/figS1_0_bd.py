
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import splev, splrep #, griddata # interp1d


params = {
        'font.family' : 'Times'
        }
matplotlib.rcParams.update(params)
plt.rcParams['text.usetex'] = True

datadir2write = "/Users/yjpark/YJ/Projects/Multisystem/Figures/figS1/sources/"

ensref = -6.7150814  # BNBNBN (1.08, 1.08)


data1   = np.genfromtxt(f"{datadir2write}data_efield0.00/data_ref/energy_wo_coul.dat", 
                        skip_header=1)
data1_ref= np.genfromtxt(f"{datadir2write}data_efield0.00/data_ref/energy_hBN_ref.dat")
print(data1_ref[:,0], data1_ref[:,-1]/data1_ref[:,-3])
cnd1 = data1_ref[:,0]==2.505759
print(data1_ref[cnd1,-1]/data1_ref[cnd1,-3])
Emin_hBN_ref    = - 6.68997268
Emin_hBN_shield = - 6.68997268
Emin_hBN_long   = -11.8854476


rcut, lambda0, TapOX = 32, 3.2, 'True'
data2   = np.genfromtxt(f"{datadir2write}data_efield0.00/data_coul_shield/energy_rcut{rcut}_lambda{lambda0}_Tap{TapOX}.dat", 
                        skip_header=1)
data2_ref= np.genfromtxt(f"{datadir2write}data_efield0.00/data_coul_shield/energy_hBN_rcut{rcut}_lambda{lambda0}_Tap{TapOX}.dat")
print(data2_ref[:,0], data2_ref[:,-1]/data2_ref[:,-3])
cnd2 = data2_ref[:,0]==2.505759
print(data2_ref[cnd2,-1]/data2_ref[cnd2,-3])

rcut, fcut = 16, -7
data3   = np.genfromtxt(f"{datadir2write}data_efield0.00/data_coul_long/energy_rcut{rcut}_fcut{fcut}.dat", 
                        skip_header=1)

data3_ref= np.genfromtxt(f"{datadir2write}data_efield0.00/data_coul_long/energy_hBN_rcut{rcut}_fcut{fcut}.dat")
print(data3_ref[:,0], data3_ref[:,-1]/data3_ref[:,-3])
cnd3 = data3_ref[:,0]==2.505759
print(data3_ref[cnd3,-1]/data3_ref[cnd3,-3])


angref = 1.084549


data2adjust = np.genfromtxt(f"{datadir2write}data_efield0.00/Hermann_Indices_t3BN.txt", 
                        skip_header=2, usecols=(0,1, 6,7, 10))  # (ang12, ang32, a", b", lambda_multiple)




# Emin_hBN_ref  = -6.68997268334775; print(f"Emin_hBN      = {Emin_hBN_ref:16.12f} eV/atom at 2.505759 \\AA")


E_L3_ref      = Emin_hBN_ref   # take the reference energy for the top-layer 



cond1 = np.logical_and( data1[:,0] == angref , data1[:,1] == angref ); print(data1[cond1,0:2])
cond2 = np.logical_and( data2[:,0] == angref , data2[:,1] == angref ); print(data2[cond2,0:2])
cond3 = np.logical_and( data3[:,0] == angref , data3[:,1] == angref ); print(data3[cond3,0:2])

eref1 = (data1[cond1,-1]/data1[cond1,-3])[0]*1000; print(f'w/o coulomb      case : {eref1/1000:18.10f} eV/atom; E_ref       = {eref1-ensref*1000:12.6f} meV/atom; {eref1-ensref*1000-(eref1-ensref*1000):12.6f} meV/atom')
eref2 = (data2[cond2,-1]/data2[cond2,-3])[0]*1000; print(f'coul/shield      case : {eref2/1000:18.10f} eV/atom; E_ref^shield= {eref2-ensref*1000:12.6f} meV/atom; {eref2-ensref*1000-(eref1-ensref*1000):12.6f} meV/atom')
eref3 = (data3[cond3,-1]/data3[cond3,-3])[0]*1000; print(f'coul/long        case : {eref3/1000:18.10f} eV/atom; E_ref^long  = {eref3-ensref*1000:12.6f} meV/atom; {eref3-ensref*1000-(eref1-ensref*1000):12.6f} meV/atom')

fig      = plt.figure(figsize=(6,6))
msize, fsize = 7, 16
msize, fsize = 4, 16

ax1 = fig.add_axes([0.16, 0.63, 0.62, 0.35])
ax2 = fig.add_axes([0.16, 0.1 , 0.62, 0.35])

for iangg, angg in enumerate( [1.084549, 1.538500, 2.004628] ):
    cond1 = data1[:,0] == angg 
    cond2 = data2[:,0] == angg 
    cond3 = data3[:,0] == angg

    ang12s_1, ang32s_1 =  data1[ cond1 ,0], data1[ cond1 ,1]
    ang12s_2, ang32s_2 =  data2[ cond2 ,0], data2[ cond2 ,1]
    ang12s_3, ang32s_3 =  data3[ cond3 ,0], data3[ cond3 ,1]

    # --- Just to ensure to use the corresponding Hermann indices
    data4corr0 = data2adjust[data2adjust[:,0]==angg,:] 
    data4corr1, data4corr2, data4corr3 = [],[],[]
    for ang12_1, ang32_1, ang12_2, ang32_2, ang12_3, ang32_3 in zip(
        ang12s_1, ang32s_1,
        ang12s_2, ang32s_2,
        ang12s_3, ang32s_3 ):

        cond2check1 = np.logical_and( np.isclose( data4corr0[:,0], ang12_1, atol=1e-4),
                                      np.isclose( data4corr0[:,1], ang32_1, atol=1e-4) )
        cond2check2 = np.logical_and( np.isclose( data4corr0[:,0], ang12_2, atol=1e-4),
                                      np.isclose( data4corr0[:,1], ang32_2, atol=1e-4) )
        cond2check3 = np.logical_and( np.isclose( data4corr0[:,0], ang12_3, atol=1e-4),
                                      np.isclose( data4corr0[:,1], ang32_3, atol=1e-4) )

        if np.sum(cond2check1) == 1:
            data4corr1.append(data4corr0[cond2check1,:].squeeze(axis=0))
        if np.sum(cond2check2) == 1:
            data4corr2.append(data4corr0[cond2check2,:].squeeze(axis=0))
        if np.sum(cond2check3) == 1:
            data4corr3.append(data4corr0[cond2check3,:].squeeze(axis=0))
    data4corr1 = np.array(data4corr1); print(data4corr1)
    data4corr2 = np.array(data4corr2); print(data4corr2)
    data4corr3 = np.array(data4corr3); print(data4corr3)

    # --- Number-of-atom & Energy corrections
    app, bpp, lambda_ratio = data4corr1[:,2], data4corr1[:,3],  data4corr1[:,-1]
    Natom_L3     = 2*(app**2 + app*bpp + bpp**2)
    Natom_L3_ref = Natom_L3[data4corr1[:,1]==angg] 
    Natom_corr_1 = (lambda_ratio**2) * Natom_L3_ref - Natom_L3
    E_corr_1     = Natom_corr_1  * Emin_hBN_ref ; print(data1[cond1,-3], )
    etots1_wocorr= (data1[cond1,-1]         )/(data1[cond1,-3]             )
    etots1       = (data1[cond1,-1]+E_corr_1)/(data1[cond1,-3]+Natom_corr_1)

    app, bpp, lambda_ratio = data4corr2[:,2], data4corr2[:,3],  data4corr2[:,-1]
    Natom_L3     = 2*(app**2 + app*bpp + bpp**2)
    Natom_L3_ref = Natom_L3[data4corr2[:,1]==angg] 
    Natom_corr_2 = (lambda_ratio**2) * Natom_L3_ref - Natom_L3
    E_corr_2     = Natom_corr_2  * Emin_hBN_shield 

    etots2_wocorr= (data2[cond2,-1]         )/(data2[cond2,-3])
    etots2       = (data2[cond2,-1]+E_corr_2)/(data2[cond2,-3]+Natom_corr_2)

    app, bpp, lambda_ratio = data4corr3[:,2], data4corr3[:,3],  data4corr3[:,-1]
    Natom_L3     = 2*(app**2 + app*bpp + bpp**2)
    Natom_L3_ref = Natom_L3[data4corr3[:,1]==angg] 
    Natom_corr_3 = (lambda_ratio**2) * Natom_L3_ref - Natom_L3
    E_corr_3     = Natom_corr_3  * Emin_hBN_long 

    etots3_wocorr= (data3[cond3,-1]         )/(data3[cond3,-3]             )
    etots3       = (data3[cond3,-1]+E_corr_3)/(data3[cond3,-3]+Natom_corr_3)
        


    if iangg == 0:
        
        ax1.plot(ang32s_1, (etots1       )*1000-eref1, 'o:' , label="Without \\texttt{coul/shield}",  lw=0.8, color=f'C{iangg}',ms=msize, alpha=0.3)
        ax1.plot(ang32s_2, (etots2       )*1000-eref2, '^-' , label="With    \\texttt{coul/shield}",  lw=0.8, color=f'C{iangg}',ms=msize-1)#, alpha=0.3)
        ax2.plot(ang32s_1, (etots1       )*1000-eref1, 'o:' , label="Without \\texttt{coul/long}"  ,  lw=0.8, color=f'C{iangg}',ms=msize, alpha=0.3)
        ax2.plot(ang32s_3, (etots3       )*1000-eref3, 's-' , label="With    \\texttt{coul/long}"  ,  lw=0.8, color=f'C{iangg}',ms=msize-1)#, alpha=0.3)

        th,ens = ang32s_1, (etots1       )*1000-eref1
        ind    = list(th).index(angg)
        spl0   = splrep(th[[0,-1]], ens[[0,-1]], k=1)
        ax1.fill_between(th, ens, splev(th,spl0) , color=f'C{iangg}', alpha=0.2)  # ax.fill_betweenx(y, x1, x2, color='k', alpha=0.3)

        th,ens = ang32s_2, (etots2       )*1000-eref2
        ind    = list(th).index(angg)
        spl0   = splrep(th[[0,-1]], ens[[0,-1]], k=1)
        ax1.fill_between(th, ens, splev(th,spl0) , color=f'C{iangg}', alpha=0.4)  # ax.fill_betweenx(y, x1, x2, color='k', alpha=0.3)

        th,ens = ang32s_1, (etots1       )*1000-eref1
        ind    = list(th).index(angg)
        spl0   = splrep(th[[0,-1]], ens[[0,-1]], k=1)
        ax2.fill_between(th, ens, splev(th,spl0) , color=f'C{iangg}', alpha=0.2)  # ax.fill_betweenx(y, x1, x2, color='k', alpha=0.3)

        th,ens = ang32s_3, (etots3       )*1000-eref3
        ind    = list(th).index(angg)
        spl0   = splrep(th[[0,-1]], ens[[0,-1]], k=1)
        ax2.fill_between(th, ens, splev(th,spl0) , color=f'C{iangg}', alpha=0.4)  # ax.fill_betweenx(y, x1, x2, color='k', alpha=0.3)

    else:
        
        ax1.plot(ang32s_1, (etots1       )*1000-eref1, 'o:', lw=0.8, color=f'C{iangg}',ms=msize, alpha=0.3)
        ax1.plot(ang32s_2, (etots2       )*1000-eref2, '^-', lw=0.8, color=f'C{iangg}',ms=msize-1)#, alpha=0.3)
        ax2.plot(ang32s_1, (etots1       )*1000-eref1, 'o:', lw=0.8, color=f'C{iangg}',ms=msize, alpha=0.3)
        ax2.plot(ang32s_3, (etots3       )*1000-eref3, 's-', lw=0.8, color=f'C{iangg}',ms=msize-1)#, alpha=0.3)

        th,ens = ang32s_1, (etots1       )*1000-eref1
        ind    = list(th).index(angg)
        spl0   = splrep(th[[0,-1]], ens[[0,-1]], k=1)
        
        ax1.fill_between(th, ens, splev(th,spl0) , color=f'C{iangg}', alpha=0.2)  # ax.fill_betweenx(y, x1, x2, color='k', alpha=0.3)

        th,ens = ang32s_2, (etots2       )*1000-eref2
        ind    = list(th).index(angg)
        spl0   = splrep(th[[0,-1]], ens[[0,-1]], k=1)
        
        ax1.fill_between(th, ens, splev(th,spl0) , color=f'C{iangg}', alpha=0.4)  # ax.fill_betweenx(y, x1, x2, color='k', alpha=0.3)

        th,ens = ang32s_1, (etots1       )*1000-eref1
        ind    = list(th).index(angg)
        spl0   = splrep(th[[0,-1]], ens[[0,-1]], k=1)
        
        ax2.fill_between(th, ens, splev(th,spl0) , color=f'C{iangg}', alpha=0.2)  # ax.fill_betweenx(y, x1, x2, color='k', alpha=0.3)

        th,ens = ang32s_3, (etots3       )*1000-eref3
        ind    = list(th).index(angg)
        spl0   = splrep(th[[0,-1]], ens[[0,-1]], k=1)
        
        ax2.fill_between(th, ens, splev(th,spl0) , color=f'C{iangg}', alpha=0.4)  # ax.fill_betweenx(y, x1, x2, color='k', alpha=0.3)

    ax1.axvline(x=angg, ls='-', color=f'C{iangg}', alpha=0.5)
    ax2.axvline(x=angg, ls='-', color=f'C{iangg}', alpha=0.5)

    
ax1.axhline(y=0, ls='-', color='gray', alpha=0.5)
ax2.axhline(y=0, ls='-', color='gray', alpha=0.5)

ax1.set_xlabel("$\\theta_{32}\\ (^{\\circ}$)", fontsize=fsize, fontname='times')
ax2.set_xlabel("$\\theta_{32}\\ (^{\\circ}$)", fontsize=fsize, fontname='times')


ax1.set_xlim(0.5,2.5); ax1.set_xticks(np.linspace(0.5, 2.5,5, endpoint=True))
ax2.set_xlim(0.5,2.5); ax2.set_xticks(np.linspace(0.5, 2.5,5, endpoint=True))


ax1.set_ylim(-0.5,2.5); ax1.set_yticks(np.linspace(-0.5, 2.5, 7, endpoint=True))
ax2.set_ylim(-0.5,2.5); ax2.set_yticks(np.linspace(-0.5, 2.5, 7, endpoint=True))

ax1.legend(  fontsize=fsize-3, loc='upper left' , ncol=1, framealpha=1)#, title="No Coulomb")
ax2.legend(  fontsize=fsize-3, loc='upper left' , ncol=1, framealpha=1)#, title="No Coulomb")

ax1.tick_params(labelsize=fsize) 
ax2.tick_params(labelsize=fsize) 
ax1.set_ylabel("$E_{\\rm tot}$ (meV/atom)", fontsize=fsize, fontname='times')
ax2.set_ylabel("$E_{\\rm tot}$ (meV/atom)", fontsize=fsize, fontname='times')

fig.savefig(f"{datadir2write}figS1_0_bd.pdf")

plt.rcParams['text.usetex'] = False