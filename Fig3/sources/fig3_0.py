import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import splev, splrep, griddata 

params = {
              'font.family' : 'Times'
}
matplotlib.rcParams.update(params)

datadir  = '/Users/yjpark/YJ/Projects/Multisystem/Figures/fig3/sources/'
datadir1 = '/Users/yjpark/YJ/Projects/Multisystem/Figures/fig3/sources/data_energymap/'
datadir2 = '/Users/yjpark/YJ/Projects/Multisystem/Figures/fig3/sources/data_energycurve/'

erefs = [
    -7.193096248220224, # GBNG
    -6.958824046577896, # BNGNB 0.53deg 기준
    -6.958824046577896  # BNGBN
    ]


# ----- Figure Generation -----
fig_width   = 6.9 # inches
golden_mean = (np.sqrt(5)-1.0)/2.0    # Aesthetic ratio
fig_height = fig_width*golden_mean*1.3 # height in inches

fig = plt.figure(figsize=(fig_width, fig_height))
msize, fsize = 3, 8
msize, fsize = 2, 8

gs0 = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1,2.4])
gs0.update(wspace=0.3, hspace=0.4) # set the spacing between axes.

# ----- Panels (a) & (c) : ENERGY MAP -----
ax1s = [fig.add_subplot(gs00)  for gs00 in gs0[0, 0].subgridspec(1,3, wspace=0.3)]
ax2s = [fig.add_subplot(gs00)  for gs00 in gs0[1, 0].subgridspec(2,1, hspace=0.4)[0].subgridspec(1,3, wspace=0.3)]
ax3s = [fig.add_subplot(gs00)  for gs00 in gs0[1, 0].subgridspec(2,1, hspace=0.4)[1].subgridspec(1,3, wspace=0.3)]

for stype, th12_degs, axs, natoms, eref in zip(['GBNG', 'BNGNB', 'BNGBN'],
                                 [[0.601428, 1.276839, 1.782679],
                                  [0.537754, 1.080027, 1.622768],
                                  [0.537754, 1.080027, 1.622768]],
                                 [ax1s, ax2s, ax3s],
                                 [[13450,7274,13950],
                                  [14310,8770,15870],
                                  [14310,8770,15870]], 
                                  erefs):
    for th12_deg, ax, natom in zip(th12_degs,
                            axs, natoms):
        val = np.genfromtxt(f"{datadir1}energies_{stype}{th12_deg:.2f}.dat")#, skip_header=1)
        XX, YY, data0 = val[:,0], val[:,1]*np.sqrt(3), (val[:,-1]/natom-eref)*1000


        min1 = np.min(data0)
        max1 = np.max(data0)
        max2 = np.max(data0[np.logical_and(abs( val[:,0]-0.5)<0.2, abs(val[:,1]-0.5)<0.2)])
        # print(min1,max2, (max2-min1))

        vvmin, vvmax  = np.min(data0), np.max(data0) 
        if stype == "BNGNB":
            vvmin, vvmax  = min1, max1
        else:
            vvmin, vvmax  = min1, max2
        # ax.scatter(XX,YY, s=msize, c=data0, vmin=vvmin, vmax=vvmax, cmap='jet')
        print(stype, th12_deg, f"{vvmin:.2f}", f"{vvmax:.2f}")

        num       = 1601
        margin    = 0.1
        L         = np.sqrt(3)
        x         = np.linspace(0, 1         ,num, endpoint=True)
        y         = np.linspace(0, np.sqrt(3),num, endpoint=True)

        val_sp    = griddata((XX, YY),data0,
                                (x[None,:], y[:,None]), method='linear') # , method='cubic')

        im = ax.imshow(val_sp, vmin=vvmin, vmax = vvmax,  #, cmap='gist_ncar')
                        extent=[0, 1,0, np.sqrt(3)],origin="lower",cmap='jet')

        ax.set_xlim(0,1)
        ax.set_ylim(0, np.sqrt(3))
        ax.set_xticks([0, 1])
        ax.set_yticks(np.linspace(0, np.sqrt(3), 4, endpoint=True))
        ax.set_xticklabels([])
        ax.set_yticklabels([])


# ----- Panels (b) & (d) : ENERGY CURVE -----
axRs = [fig.add_subplot(gs00)  for gs00 in [gs0[0, 1], gs0[1, 1].subgridspec(2, 1, hspace=0.4)[0], gs0[1, 1].subgridspec(2, 1, hspace=0.4)[1] ]]

inds_GBNG    = [[0,1,2,3,4,5,6,7,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1], # 0.60
                [0,1,2,3,4,5,6,7,8,9,-3,-2,-1], # 1.28
                [0,1,2,3,4,5,-5,-4,-3,-2,-1]]   # 1.78
inds_BNGNB   = [[0,1,2,3,4,5,-4,-3,-2,-1],  # 0.54
                [0,1,2,3,-6,-5,-4,-3,-2,-1],  # 1.08
                [0,1,2,3,-5,-4,-3,-2,-1]]  # 1.62
inds_BNGBN   = [[0,1,2,3,4,5,-4,-3,-2,-1],  # 0.54
                [0,1,2,3,-6,-5,-4,-3,-2,-1],  # 1.08
                [0,1,2,3,-5,-4,-3,-2,-1]]  # 1.62


print("ENERGY CURVE:")
# ----- Added To estimate the corrections in the number of top-layer atoms and their energies
data2adjust_GBNG   = np.genfromtxt(f"{datadir2}Hermann_Indices_GBNG.txt", 
                        skip_header=2, usecols=(0,1, 6,7, 10))  # (ang12, ang32, a", b", lambda_multiple)
data2adjust_BNGBN  = np.genfromtxt(f"{datadir2}Hermann_Indices_BNGBN.txt", 
                        skip_header=2, usecols=(0,1, 6,7, 10))  # (ang12, ang32, a", b", lambda_multiple)
data2adjusts       = [data2adjust_GBNG, data2adjust_BNGBN, data2adjust_BNGBN]


Emin_graphene = -7.3949972515997 ; print(f"Emin_graphene = {Emin_graphene:16.12f} eV/atom at 2.46019 \\AA")
Emin_hBN      = -6.68997268334775; print(f"Emin_hBN      = {Emin_hBN:16.12f} eV/atom at 2.505759 \\AA")
E_L3_refs     = [Emin_graphene, Emin_hBN, Emin_hBN]   # take the reference energy for the top-layer 

for stype, ax, mtype, stack_name, spl_inds, eref, th12_degs, th32_deg_comms, E_L3_ref, data2adjust  in zip(['GBNG', 'BNGNB', 'BNGBN'],
                            axRs, 
                            ['o','^','s'],
                            ["$\\overline{\\text{AAC}}$", "$\\overline{\\rm \\text{A}^{\\prime}\\text{AA}}$", "$\\overline{\\text{AAC}}$"],
                            [inds_GBNG, inds_BNGNB, inds_BNGBN],
                            erefs,
                            [[0.601428, 1.276839, 1.782679],
                             [0.537754, 1.080027, 1.622768],
                             [0.537754, 1.080027, 1.622768]],
                            [[0.601428, 1.276839, 1.782679],
                             [0.537754, 1.080027, 1.622768],
                             [0.537754, 1.080027, 1.622768]],
                             E_L3_refs,
                             data2adjusts
                            ):
    
    vals = np.genfromtxt(f"{datadir2}energycurve_{stype}.dat", skip_header=1)

    for th12_deg, th32_deg_comm, mcolor, spl_ind in zip(th12_degs,th32_deg_comms,
                                ['C0', 'C1', 'C2'],
                                spl_inds):

        # --- Selects the specific \theta_{12}  
        cond = vals[:,0] == th12_deg 
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
        data4corr = np.array(data4corr) #;  print("data4corr.shape",data4corr.shape)

        # --- Number-of-atom correction
        app, bpp          = data4corr[:,2], data4corr[:,3]
        Natom_L3          = 2*(app**2 + app*bpp + bpp**2)
        Natom_L3_ref      = Natom_L3[data4corr[:,1]==th32_deg_comm] 
        lambda_ratio      = data4corr[:,-1]  # ratio between the simulation cell length with respect to their reference (double-moire commensurate)
        Natom_corr        = (lambda_ratio**2) * Natom_L3_ref - Natom_L3
        # print(f"Natom_corr  =  (lambda/lambda^comm.)^2 * Natom_L3_ref - Natom_L3 ", f"ref angle = {data4corr[data4corr[:,1]==th32_deg_comm, 0]}")
        
        # --- Energy correction
        E_corr  = Natom_corr  * E_L3_ref 

        # --- Total energy per atom with Correction
        etots_wocorr = (val[:,-1])/(val[:,-3])
        etots        = (val[:,-1]+E_corr)/(val[:,-3]+Natom_corr)
        # ---------------------------------------------------

        # ax.plot(ang32_degs, (etots_wocorr-eref)*1000,'*', ms= msize, color=mcolor, alpha=0.5)# label="$\\theta_{12} = "+f"{th12_deg:.2f}"+"^{\\circ}$ - "+f"{stack_name}")
        ax.plot(ang32_degs, (etots-eref)*1000,'-', ms= msize, marker=mtype, color=mcolor, lw=0.5, label="$\\theta_{12} = "+f"{th12_deg:.2f}"+"^{\\circ}$ - "+f"{stack_name}")



        ind   = list(ang32_degs).index(th12_deg)
        if th12_deg < 0.8:
            thL   = np.linspace(ang32_degs[ind]-1,ang32_degs[ind])
            thR   = np.linspace(ang32_degs[ind],ang32_degs[ind]+1)
        else:
            thL   = np.linspace(ang32_degs[ind]-0.6,ang32_degs[ind])
            thR   = np.linspace(ang32_degs[ind],ang32_degs[ind]+0.6)

        zL  = np.polyfit(np.array([ang32_degs[0],ang32_degs[int((ind-0)/2)],ang32_degs[ind]]),
                        np.array([ etots[0], etots[int((ind-0)/2)], etots[ind]]), 2)
        pL  = np.poly1d(zL)
        fL  = pL(thL)
        ax.plot(thL,(fL-eref)*1000,"-",color=mcolor,linewidth=0.4)


        zR  = np.polyfit(np.array([ang32_degs[ind],ang32_degs[int((len(ang32_degs)-1+ind)/2)],ang32_degs[-1]]),
                        np.array([ etots[ind], etots[int((len(ang32_degs)-1+ind)/2)], etots[-1]]), 2)
        pR  = np.poly1d(zR)        
        fR  = pR(thR)
        ax.plot(thR,(fR-eref)*1000,"-",color=mcolor,linewidth=0.4)

        thLR  = np.hstack((thL[:5],thR[-5:]))
        ensLR = np.hstack((fL[:5],fR[-5:]))

        zLR   = np.polyfit(thLR,ensLR, 2)

        pLR   = np.poly1d(zLR)
        fitted = pLR(ang32_degs)

        theta1 = np.linspace(thLR[0], thLR[-1])
        theta2 = np.linspace(ang32_degs[0], ang32_degs[-1])
        ax.plot(theta1,(pLR(theta1)-eref)*1000,":",color=mcolor,linewidth=0.9)
        ax.plot(theta2,(pLR(theta2)-eref)*1000,"--",color=mcolor,linewidth=0.9)

        pLR0_val = np.hstack((pL(theta1[theta1<th12_deg]), pR(theta1[theta1>=th12_deg])))#; print(pLR0_val.shape, theta1.shape)
        ax.fill_between(theta1, (pLR0_val-eref)*1000, (pLR(theta1)-eref)*1000 , color=mcolor, alpha=0.3)  # ax.fill_betweenx(y, x1, x2, color='k', alpha=0.3)


        spl    = splrep(theta2,pLR(theta2))

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
        # if stype!='BNGNB': ax.set_xlabel("$\\theta_{32}$", fontsize=fsize, fontname='times')

        ax.set_xlim(0,2.3)
        ax.set_ylim(-0.05,0.85)
        ax.set_xticks(np.linspace(0, 2, 5, endpoint=True))
        ax.set_yticks(np.linspace(0, 0.8, 5, endpoint=True))
        ax.set_xticklabels([f"{i:.1f}"+"$^{\\circ}$" for i in np.linspace(0, 2, 5, endpoint=True)])
        if stype=='GBNG': 
            ax.set_title(f'{stype[0]}/{stype[1:3]}/{stype[3]}', fontsize=fsize+1, fontname='times')
        else:
            ax.set_title(f'{stype[:2]}/{stype[2:3]}/{stype[3:]}', fontsize=fsize+1, fontname='times')
        ax.legend(loc=0, fontsize=fsize-3, framealpha=1)
        ax.tick_params(labelsize=fsize)
        ax.vlines(th12_deg,
                   ymin=0, 
                   ymax=0.85,  
                   alpha=0.5, lw=2,
                   color=mcolor)
        ax.grid('on', alpha=0.5)


fig.savefig(f"{datadir}fig3_0.pdf")
