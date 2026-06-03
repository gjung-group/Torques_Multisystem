import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import splev, splrep, griddata 

params = {
            'font.family' : 'Times'
        }
matplotlib.rcParams.update(params)


datadir  = '/Users/yjpark/YJ/Projects/Multisystem/Figures/fig1/sources/'
datadir1 = '/Users/yjpark/YJ/Projects/Multisystem/Figures/fig1/sources/data_energymap/'
datadir2 = '/Users/yjpark/YJ/Projects/Multisystem/Figures/fig1/sources/data_energycurve/'


eref = -6.7150814  # BNBNBN (1.08, 1.08)


# ----- Figure Generation -----
fig_width   = 6.9 # inches
golden_mean = (np.sqrt(5)-1.0)/2.0    # Aesthetic ratio
fig_height = fig_width*golden_mean*1.3 # height in inches

fig = plt.figure(figsize=(fig_width, fig_height))
msize, fsize = 2, 8

gs0 = gridspec.GridSpec(3, 2, figure=fig)
gs0.update(wspace=0.3, hspace=0.5) # set the spacing between axes.

# ----- Panel (a) : ENERGY MAP -----
ax1s = [fig.add_subplot(gs00)  for gs00 in gs0[0, 0].subgridspec(1, 3)]
ax2s = [fig.add_subplot(gs00)  for gs00 in gs0[1, 0].subgridspec(1, 3)]
ax3s = [fig.add_subplot(gs00)  for gs00 in gs0[2, 0].subgridspec(1, 3)]

print("ENERGY MAP:")
print("Type, theta12(eV), Emin(eV/atom), Emax(eV/atom)")
for stype, axs in zip(['BNBNBN', 'NBBNBN', 'BNNBBN'],
                       [ax1s, ax2s, ax3s]):
    for th12_deg, natom, ax in zip([1.084549, 1.538500, 2.004628],
                                 [16746, 8322, 4902],
                            axs):
        val = np.genfromtxt(f"{datadir1}energymap_{stype}_{th12_deg:.2f}.dat")#, skip_header=1)
        
        XX, YY, data0 = val[:,0], val[:,1]*np.sqrt(3), (val[:,-1]/natom-eref)*1000
        vvmin, vvmax  = np.min(data0), np.max(data0) 
        print(stype, f"{th12_deg:.6f}", f"{vvmin:.2f}", f"{vvmax:.2f}")
        # ax.scatter(XX,YY, s=msize, c=data0, vmin=vvmin, vmax=vvmax, cmap='jet')

        num       = 1601
        margin    = 0.1
        L         = np.sqrt(3)
        x         = np.linspace(0, 1         ,num, endpoint=True)
        y         = np.linspace(0, np.sqrt(3),num, endpoint=True)

        val_sp    = griddata((XX, YY),data0,
                                (x[None,:], y[:,None]), method='linear') 

        im = ax.imshow(val_sp, vmin=vvmin, vmax = vvmax,  
                        extent=[0, 1,0, np.sqrt(3)],origin="lower",cmap='jet')

        ax.set_xlim(0,1)
        ax.set_ylim(0, np.sqrt(3))
        ax.set_xticks([0, 1])
        ax.set_yticks(np.linspace(0, np.sqrt(3), 4, endpoint=True))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if stype=='BNBNBN': ax.set_title("$\\theta_{12}="+f"{th12_deg:.2f}"+"^{\\circ}$", fontsize=fsize, fontname='times')


# ----- Panel (b) : ENERGY CURVE -----
axRs = [fig.add_subplot(gs0[i, 1]) for i in range(3)] 

inds_BNBNBN   = [[0,1,2,3,-3,-2,-1], # 1.08
                 [0,1,2,3,-2,-1], # 1.54
                 [0,1,2,-3,-2,-1]]   # 2.00
inds_NBBNBN   = [[0,1,-3,-2,-1],  # 1.08
                 [0,1,2,-2,-1],   # 1.54
                 [0,1,-3,-1]]     # 2.00
inds_BNNBBN   = [[0,1,2,-3,-2,-1],   # 1.08
                 [0,1,2,3,4,-2,-1], # 1.54
                 [0,1,2,-3,-2,-1]]   # 2.00


print("ENERGY CURVE:")
# ----- Added To estimate the corrections in the number of top-layer atoms and their energies
data2adjust   = np.genfromtxt(f"{datadir2}Hermann_Indices_t3BN.txt", 
                        skip_header=2, usecols=(0,1, 6,7, 10))  # (ang12, ang32, a", b", lambda_multiple)

Emin_graphene = -7.3949972515997 ; print(f"Emin_graphene = {Emin_graphene:16.12f} eV/atom at 2.46019 \\AA")
Emin_hBN      = -6.68997268334775; print(f"Emin_hBN      = {Emin_hBN:16.12f} eV/atom at 2.505759 \\AA")
E_L3_ref      = Emin_hBN   # take the reference energy for the top-layer 

for stype, ax, mtype, stack_name, spl_inds in zip(['BNBNBN', 'NBBNBN', 'BNNBBN'],
                            axRs, 
                            ['o','^','s'],
                            ["$\\overline{\\rm AAA}$", "$\\overline{\\rm AAC^\\prime}$", "$\\overline{\\rm AA^\\prime A}$"],
                            [inds_BNBNBN, inds_NBBNBN, inds_BNNBBN]):
    
    vals = np.genfromtxt(f"{datadir2}energycurve_{stype}.dat", skip_header=1)

    for th12_deg, th32_deg_comm, mcolor, spl_ind in zip(
        [1.084549, 1.538500, 2.004628],
        [1.084549, 1.538500, 2.004628],
        ['C0', 'C1', 'C2'],
        spl_inds):
        
        # --- Selects the specific \theta_{12}  
        cond = vals[:,0] == th12_deg #np.isclose( vals[:,0] , th12_deg)
        val  = vals[cond,:].copy()
        if np.sum(cond)==0: continue  # Just in case that no one is selected.

        ang12_degs, ang32_degs = val[:,0], val[:,1]

        # ---------------------------------------------------
        # --- Just to ensure to use the corresponding Hermann indices
        data4corr0 = data2adjust[data2adjust[:,0]==th12_deg,:] 
        data4corr  = []
        for ang12, ang32 in zip(ang12_degs, ang32_degs):
            cond2check = np.logical_and( np.isclose( data4corr0[:,0], ang12, atol=1e-4),
                                         np.isclose( data4corr0[:,1], ang32, atol=1e-4) )
            if np.sum(cond2check) == 1:
                data4corr.append(data4corr0[cond2check,:].squeeze(axis=0))
        data4corr = np.array(data4corr)

        # --- Number-of-atom correction
        app, bpp          = data4corr[:,2], data4corr[:,3]
        Natom_L3          = 2*(app**2 + app*bpp + bpp**2)
        Natom_L3_ref      = Natom_L3[data4corr[:,1]==th32_deg_comm] 
        lambda_ratio      = data4corr[:,-1]  # ratio between the simulation cell length with respect to their reference (double-moire commensurate)
        Natom_corr        = (lambda_ratio**2) * Natom_L3_ref - Natom_L3
        
        # --- Energy correction
        E_corr  = Natom_corr  * E_L3_ref 

        # --- Total energy per atom with Correction
        etots_wocorr = (val[:,-1])/(val[:,-3])
        etots        = (val[:,-1]+E_corr)/(val[:,-3]+Natom_corr)
        # ---------------------------------------------------

        # ax.plot(ang32_degs, (etots_wocorr-eref)*1000,'*', ms= msize, color=mcolor, alpha=0.5)# label="$\\theta_{12} = "+f"{th12_deg:.2f}"+"^{\\circ}$ - "+f"{stack_name}")
        ax.plot(ang32_degs, (etots-eref)*1000,'-', ms= msize, marker=mtype, color=mcolor, lw=0.5, label="$\\theta_{12} = "+f"{th12_deg:.2f}"+"^{\\circ}$ - "+f"{stack_name}")

        ind   = list(ang32_degs).index(th12_deg)

        # Fitting curves
        spl    = splrep(ang32_degs[spl_ind],etots[spl_ind])
        theta  = np.linspace(np.min(ang32_degs),np.max(ang32_degs))
        fitted = splev(theta,spl)
        
        ax.plot(theta,(fitted-eref)*1000,"--",color=mcolor,linewidth=0.9)
        ax.fill_between(ang32_degs, (etots-eref)*1000, (splev(ang32_degs,spl)-eref)*1000 , color=mcolor, alpha=0.3)  # ax.fill_betweenx(y, x1, x2, color='k', alpha=0.3)


        print(f"{stype}_{th12_deg:.2f}")
        ind=list(ang32_degs).index(th12_deg)
        ang = ang32_degs
        val = (etots-eref)*1000
        slope_L = (val[ind]-val[ind-1])/np.deg2rad(ang[ind]-ang[ind-1])
        slope_R = (val[ind+1]-val[ind])/np.deg2rad(ang[ind+1]-ang[ind])
        print(f"Slope: {ang[ind-1]},{ang[ind]},{ang[ind+1]}: {slope_L:.2f}, {slope_R:.2f}")

        ind  = np.where(ang32_degs == th12_deg)[0][0]
        val1 = (splev(np.array([th12_deg]),spl)-eref)*1000
        print(f"Binding E:  {ang[ind]:.6f}, {val[ind]:.3f}, {val1[0]:.3f}, {val1[0]-val[ind]:.3f} \n")

        ax.set_ylabel("$E^{tot}$ (meV/atom)", fontsize=fsize, fontname='times')
        if stype=='BNNBBN': ax.set_xlabel("$\\theta_{32}$", fontsize=fsize, fontname='times')

        ax.set_xlim(0,3.1)
        ax.set_xticks(np.linspace(0, 3, 7, endpoint=True))

        if stype=='BNBNBN': 
            ax.set_ylim(-1,2)
            ax.set_yticks(np.linspace(-1, 2, 7, endpoint=True))
        else:
            ax.set_ylim(-0.5,2.5)
            ax.set_yticks(np.linspace(-0.5, 2.5, 7, endpoint=True))
             
        ax.set_xticklabels([f"{i:.1f}"+"$^{\\circ}$" for i in np.linspace(0, 3, 7, endpoint=True)])
        ax.set_title(f'{stype[:2]}/{stype[2:4]}/{stype[4:]}', fontsize=fsize, fontname='times')
        ax.legend(loc='lower right', fontsize=fsize-3, framealpha=1)
        ax.tick_params(labelsize=fsize)
        ax.vlines(th12_deg,
                   ymin=-1, 
                   ymax=(etots[ind]-eref)*1000,  
                   alpha=0.5, lw=2,
                   color=mcolor)
        ax.grid('on', alpha=0.5)


fig.savefig(f"{datadir}fig1_0.pdf")
