### Functions written by research students for analysis of MADNT simulation
### Things such as plotting or simplifying loops are here instead of __init_student__ to avoid bloat

###loops
def save_to_np():
    '''Written by Shu Xin Wu:
    this is our for loop to gather all of our data
    adjust filename'''
    # 347 or 158
    outputs = np.zeros((347, 14))
    fnumber = 10333
    for x in range(347):
        outputs[x] = Rstressdecompvtime(fnumber + x)
    np.save("Rstressvtime1.npy", outputs)

def make_diskstress(fstart, fend, rng):
    ''' Written by Shu Xin Wu
    cycles through diskstress '''
    for f in range(fstart, fend):
        diskstress(f, rng)
    print("****** process completed ******")
    
###plotting
def ups_plot(fname):
    ''' Written by Shu Xin Wu
        Creates a plot of upsilon, takes the filename of output plot.
        Adjust file number selection and name of npy file. Has options
        to plot amag as well. '''

    mdot = np.load("mdotuph_vs_t.npz")
    t = mdot["ts"][10333:10838]
    ups = mdot["phibh"][10333:10838]
    print(len(t))
    print(len(ups))
    outputs = np.load("stressvtime10333.npy")
    amag = outputs[:,1]
    u_max = max(ups)
    u_min = min(ups)
    u_norm = (ups - u_min)/(u_max - u_min)

    plt.clf()
    plt.plot(t, ups)
    #plt.plot(t,u_norm)
    #plt.plot(t,amag)
    plt.xlabel("Time [$r_g$/c]")
    plt.ylabel(r"$\Upsilon$")
    #plt.title("Stress and Magnetic Flux")
    plt.title("Upsilon")
    #plt.legend(["Upsilon", "Total alpha"])
    plt.savefig(fname + ".png")
    plt.close()

def plot_a_plot(plotees, foutname, plottitle):
    ''' Written by Shu Xin Wu
    Takes an array of indices corresponding to outputs from stressdecompvtime ->
    these will be plotted against time, with title plottitle. Outputs to a file foutname.
    key: 1 - alpha tot, 2 - alpha mean, 3 - alpha cross 1, 4 - alpha cross 2, 5 - alpha turb
    6 - bubble total, 7 - bubble cross 1, 8 - bubble cross 2, 9 - bubble turb, 14 - bubble alpha mean
    10 - disk total, 11 - disk cross 1, 12 - disk cross 2, 13 - disk turb, 15 - disk alpha mean '''

    labels = ["Total alpha", "Alpha mean", "Total cross 1", "Total cross 2", "Total turbulent", "Bubble total", "Bubble cross 1", "Bubble cross 2", "Bubble turbulent", "Disk total", "Disk cross 1", "Disk cross 2", "Disk turbulent", "Bubble alpha mean", "Disk alpha mean"]
    symbols = [".","","s","+","x","*"] # different markers for visibiility, alpha mean must be ""
    outputs = np.load("stressvtime1.npy")
    legend = []
    plt.clf()
    time = outputs[:,0]

    # to help w visibility, plot every few (4 or 8) points
    plot_time = []
    for i in range(0,len(time), 8): # edit step size bw points
        plot_time.append(time[i]*4)
    for item in plotees:
        if item == 14: #bubble mean
            plt.plot(time, outputs[:,6] - outputs[:,7]- outputs[:,8] - outputs[:,9])
        elif item == 15: #disk mean
            plt.plot(time, outputs[:,10] - outputs[:,11]- outputs[:,12] - outputs[:,13])
        else:
            full_item = outputs[:,item]
            to_plot = []
            for i in range(0,len(full_item), 8): # edit step size bw points
                to_plot.append(full_item[i])
            plt.plot(plot_time, to_plot, marker=symbols[item-1])
        legend.append(labels[item - 1])

    plt.xlabel("Time [$r_g$/c]")
    plt.ylabel(r"$\alpha$")
    plt.ylim([-0.3,1]) # defining range for consistent y axes - amag
    # plt.ylim([-2,4]) #-- r
    plt.title(plottitle)
    plt.legend(legend)
    plt.tight_layout()
    plt.savefig(foutname + ".png")
    plt.close()

def double_plot(fln,quad=False):
    '''Written by Calen Carron
    Function to create side by side (double) or 2x2 (quad) plots
    Manually input the file locations that you would like to use to create the
    double or quad plots in the mpimg.imread() calls labeled below. '''
    import matplotlib.image as mpimg
    plt.clf()
    if quad:
        fig=plt.figure(figsize=(8, 6)) # quad plots
        rows = 2
    else:
        fig=plt.figure(figsize=(7, 2.5)) # double plots
        rows = 1
    columns = 2
    plt.subplots_adjust(wspace=0, hspace=0)
    for i in range(rows*columns):
        if i == 0:
            # first image location
            img = mpimg.imread(os.getcwd()+'/plots/totalComp.png')
        if i == 1:
            # second image location
            img = mpimg.imread(os.getcwd()+'/plots/turbComp.png')
        if i == 2:
            # third image location (for quad plot)
            img = mpimg.imread(os.getcwd()+'/plots/crossComp.png')
        if i == 3:
            # fourth image location (for quad plot)
            img = mpimg.imread(os.getcwd()+'/plots/plotComp.png')

        fig.add_subplot(rows, columns, i+1)
        plt.imshow(img)

        # if i == 0:
        #     plt.title('Plasma Density   ')
        # if i == 1:
        #     plt.title('Plasma Velocity   ')
        plt.axis('off')

    # save the new image - adjust filename/location to your organizational liking
    plt.savefig(os.getcwd()+'/plots/test1'+str(fln)+'.png', dpi=800, bbox_inches='tight')
    plt.close()

def diskstress(fnumber, rng, ncell=400):
    ''' Written by Shu Xin Wu
        Looking at reynolds stress and diff components of it throughout the disk
        Like zoomvideo but for rstress '''

    # first calculate the same way as decomp
    # first load grid file
    plt.clf()
    grid3d("gdump.bin", use2d=True)

    #load the fieldline file for a given time and compute standard quantities
    rfd("fieldline"+str(fnumber).zfill(4)+".bin")
    ###############################
    (rhoclean,ugclean,uublob,maxbsqorhonear,maxbsqorhofar,condmaxbsqorho,condmaxbsqorhorhs,rinterp)=getrhouclean(rho,ug,uu)
    cvel()
    rhor=1+(1-a**2)**0.5
    ihor=int(iofr(rhor))
    pg=(gam-1.0)*ug
    #
    diskcondition=condmaxbsqorho
    #only around equator, not far away from equator
    diskcondition=diskcondition*(bsq/rho<1.0)*(np.fabs(h-np.pi*0.5)<0.1)
    diskeqcondition=diskcondition
    ##############################

    ###choose the radial extent of the plot
    nxin=int(iofr(2.5))
    nxout=int(iofr(25))
    ###choose extent in r,theta:
    hoverr=0.1
    hmin=np.pi/2 - hoverr
    hmax=np.pi/2 + hoverr
    mhin=int(jofh(hmin,nxout))
    mhout=int(jofh(hmax,nxout))
    iny = int(ny/2)

    loadavg()
    ###integrated stress
    # finding all terms, avg and perturbations
    ibeta=0.5*bsq/pg #for masking
    ##velocity decomposition into mean field and turbulent terms
    vr=uu[1]*np.sqrt(gv3[1,1]) # b upper
    vrmean=avg_uu[1]*np.sqrt(gv3[1,1])
    vrpert=(uu[1]-avg_uu[1])*np.sqrt(gv3[1,1])

    vz=-uu[2]*np.sqrt(gv3[2,2]) #minus sign because theta hat points in -z hat direction
    vzmean=-avg_uu[2]*np.sqrt(gv3[2,2])
    vzpert=-(uu[2]-avg_uu[2])*np.sqrt(gv3[2,2])

    vphi=ud[3]*np.sqrt(gn3[3,3])
    vphimean=avg_ud[3]*np.sqrt(gn3[3,3])
    vphipert=(ud[3]-avg_ud[3])*np.sqrt(gn3[3,3])

    ## needs editing
    # from zoomvideo
    #make the spherical polar grid fully 3D
    myr3d=mk2d3d(r)
    myh3d=mk2d3d(h)
    myph3d=mk2d3d(ph)
    for k in range(0,nz):
        myph3d[:,:,k]=(0.5+k)*2*np.pi/nz
    #compute cartesian coordinates for each grid point of the array
    myx=myr3d*np.sin(myh3d)*np.cos(myph3d)
    myy=myr3d*np.sin(myh3d)*np.sin(myph3d)

    ncell = ncell
    extent = (-rng, rng, -rng, rng)

    #rstress = reynolds stress
    rstress = vrpert*vphipert # choose radial for now

    #ibeta=0.5*bsq/pg

    rstress=reinterpxy(vr,extent,ncell,domask=1,interporder='linear')
    #i_ibeta_h=reinterpxyhor(ibeta,extent,ncell,domask=1,interporder='linear')
    #i_ibeta[i_ibeta.mask==True]=i_ibeta_h[i_ibeta.mask==True]

    grid = np.linspace(-rng, rng, ncell)
    plt.clf()
    plt.figure(1)

    plt.pcolor(grid,grid, rstress, vmin=-1, vmax=1, cmap='inferno')
    c = plt.colorbar()
    #c.ax.set_yticklabels(['$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^1$','$10^2$'])

    plt.xlim(-rng,rng)
    plt.ylim(-rng,rng)
    plt.xlabel(r"$x [r_g]$",ha='center',labelpad=0,fontsize=14)
    plt.ylabel(r"$y [r_g]$",ha='left',labelpad=0,fontsize=14)
    plt.title(str(fnumber*4).zfill(4))
    plt.savefig("rstress_pert"+str(rng)+"rg_"+str(fnumber).zfill(4)+".png")
    
def CCzoomvideo(fnumber,rng,ncell=400,deg=0,choice=1):
    ''' Use this function to generate color plots.
        Based of zoomvideo (written by Megan) Modified by Calen Carron
        Inputs:
            fnumber: field line file number
            rng: max radius desired for plot in rg
            ncell: level of detail of the plot
            deg: radial cell number for xz plots (0-208). degree value = 360*(deg/208)
            choice: option value for function sub sections with more than one set of
                    values to plot (0 or 1). see subsection header comments for details
    '''
    plt.clf()
    grid3d("gdump.bin",use2d=True) #load the gdump.bin file - use2d=True to save memory
    #make the spherical polar grid fully 3D
    myr3d=mk2d3d(r)
    myh3d=mk2d3d(h)
    myph3d=mk2d3d(ph)
    for k in range(0,nz):
        myph3d[:,:,k]=(0.5+k)*2*np.pi/nz
    #compute cartesian coordinates for each grid point of the array
    myx=myr3d*np.sin(myh3d)*np.cos(myph3d)
    myy=myr3d*np.sin(myh3d)*np.sin(myph3d)
    #load the fieldline file for a given time and compute standard quantities
    rfd("fieldline"+str(fnumber).zfill(4)+".bin")
    cvel()
    rhor=1+(1-a**2)**0.5
    ihor=np.floor(iofr(rhor)+0.5)
    pg=(gam-1.0)*ug
    #set parameters for the interpolation routine
    #rng=40.0
    ncell=ncell
    extent=(-rng,rng,-rng,rng)
    #compute quantities used for masking (currently magnetic flux and inverse beta 2/24/16)
    Br = dxdxp[1,1]*B[1]+dxdxp[1,2]*B[2]
    Bh = dxdxp[2,1]*B[1]+dxdxp[2,2]*B[2]
    Bp = B[3]*dxdxp[3,3]
    # Spherical:
    Brnorm=Br
    Bhnorm=Bh*np.abs(r)
    Bpnorm=Bp*np.abs(r*np.sin(h))
    # Cylindrical
    BRnorm=Brnorm*np.sin(h)+Bhnorm*np.cos(h)
    Bznorm=Brnorm*np.cos(h)-Bhnorm*np.sin(h)

################################################################################
    # Density and Beta Plots - Section 1
    """ This section creates color plots of density and beta values in the xy plane
        Use choice=0 for i_ibeta
        Use choice=1 for ilrho
    """
################################################################################

    irho=reinterpxy(rho,extent,ncell,domask=1,interporder='linear')
    irho_h=reinterpxyhor(rho,extent,ncell,domask=1,interporder='linear')
    irho[irho.mask==True]=irho_h[irho.mask==True]

    lrho=np.log10(rho)
    ilrho=reinterpxy(lrho,extent,ncell,domask=1,interporder='linear')
    ilrho_hor=reinterpxyhor(lrho,extent,ncell,domask=1, interporder='linear')
    ilrho[ilrho.mask==True]=ilrho_hor[ilrho.mask==True]

    ibeta=0.5*bsq/pg
    i_ibeta=reinterpxy(ibeta,extent,ncell,domask=1,interporder='linear')
    i_ibeta_h=reinterpxyhor(ibeta,extent,ncell,domask=1,interporder='linear')
    i_ibeta[i_ibeta.mask==True]=i_ibeta_h[i_ibeta.mask==True]


    grid=np.linspace(-rng, rng, ncell)
    #Make plots to use as movie frames
    plt.clf()
    plt.figure(1)

    if choice == 0:
        plt.pcolor(grid,grid,i_ibeta,vmax=80,cmap='jet') #beta - change vmax for proper mask
    if choice == 1:
        plt.pcolor(grid,grid,ilrho,vmin=-4,vmax=2,cmap='inferno') #density

    cbar = plt.colorbar()
    if choice == 1:
        cbar.ax.set_yticklabels(['$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^1$','$10^2$'])
    horcir=plt.Circle((0,0),2,color='black',fill=True)
    plt.gcf().gca().add_artist(horcir)
    plt.xlim(-rng,rng)
    plt.ylim(-rng,rng)
    plt.xlabel(r"$x [r_g]$",ha='center',labelpad=0,fontsize=14)
    plt.ylabel(r"$y [r_g]$",ha='left',labelpad=0,fontsize=14)
    plt.title(str(fnumber*4).zfill(4))
    #
    # if choice == 0:
    #     plt.savefig('ibeta_'+str(rng)+'rg_'+str(fnumber).zfill(4)+'.png')
    # if choice == 1:
    #     plt.savefig('ilrho_'+str(rng)+'rg_'+str(fnumber).zfill(4)+'.png')


################################################################################
    #section 2

    """
    This section creates color plots of an xz slice of the ilrho value.
    Specify deg in zoomvideo function call if value other than 0 is wanted
    """
################################################################################
    """
    #deg = 30
    lrho=np.log10(rho)
    ilrho_z = reinterp(lrho,extent,ncell,domask=1,interporder='linear')

    grid=np.linspace(-rng, rng, ncell)
    plt.figure(1)
    plt.pcolor(grid,grid,ilrho_z,vmin=-4,vmax=2,cmap='inferno') #density
    cbar = plt.colorbar()
    cbar.ax.set_yticklabels(['$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^1$','$10^2$'])

    horcir=plt.Circle((0,0),2,color='black',fill=True)
    plt.gcf().gca().add_artist(horcir)
    plt.xlim(-rng,rng)
    plt.ylim(-rng,rng)
    plt.xlabel(r"$x [r_g]$",ha='center',labelpad=0,fontsize=14)
    plt.ylabel(r"$z [r_g]$",ha='left',labelpad=0,fontsize=14)

    plt.savefig('ilrho_xz'+str(deg)+'_'+str(rng)+'rg_'+str(fnumber).zfill(4)+'.png')
    """

################################################################################
    #reinterpolate for magnetic field strength plots in xy or xz slices - Section 3
    """ This section creates color plots for magnetic field strength
        Comment/uncomment sections marked with xy vs xz slices as needed
    """
################################################################################
    """
    # xy slices
    iBR=reinterpxy(BRnorm,extent,ncell,domask=1,interporder='linear')
    iBR_h=reinterpxyhor(BRnorm,extent,ncell,domask=1,interporder='linear')
    iBR[iBR.mask==True]=iBR_h[iBR.mask==True]

    iBz=reinterpxy(Bznorm,extent,ncell,domask=1,interporder='linear')
    iBz_h=reinterpxyhor(Bznorm,extent,ncell,domask=1,interporder='linear')
    iBz[iBz.mask==True]=iBz_h[iBz.mask==True]

    '''
    # xz slices
    iBR_xz = reinterp(BRnorm,extent,ncell,deg,domask=1,interporder='linear')
    iBz_xz = reinterp(Bznorm,extent,ncell,deg,domask=1,interporder='linear')
    '''

    grid=np.linspace(-rng, rng, ncell)
    #Make plots to use as movie frames
    plt.figure(1)

    # xy slices
    if choice == 0:
        plt.pcolor(grid,grid,iBR,vmin=-0.1,vmax=0.1,cmap='nipy_spectral')
    if choice == 1:
        plt.pcolor(grid,grid,iBz,vmin=-0.1,vmax=0.1,cmap='nipy_spectral')
    '''
    # xz slices
    if choice == 0:
        plt.pcolor(grid,grid,iBR_xz,vmin=-0.1,vmax=0.1,cmap='nipy_spectral')
    if choice == 1:
        plt.pcolor(grid,grid,iBz_xz,vmin=-0.1,vmax=0.1,cmap='nipy_spectral')
    '''

    cbar = plt.colorbar()
    horcir=plt.Circle((0,0),2,color='black',fill=True)
    plt.gcf().gca().add_artist(horcir)
    plt.xlim(-rng,rng)
    plt.ylim(-rng,rng)
    plt.xlabel(r"$x [r_g]$",ha='center',labelpad=0,fontsize=14)

    # xy slices
    plt.ylabel(r"$y [r_g]$",ha='left',labelpad=0,fontsize=14)
    if choice == 0:
        plt.savefig('iBR_'+str(rng)+'rg_'+str(fnumber).zfill(4)+'.png', bbox_inches='tight')
    if choice == 1:
        plt.savefig('iBz_'+str(rng)+'rg_'+str(fnumber).zfill(4)+'.png', bbox_inches='tight')
    '''
    # xz slices
    plt.ylabel(r"$z [r_g]$",ha='left',labelpad=0,fontsize=14)
    if choice == 0:
        plt.savefig('iBR_xz'+str(deg)+'_'+str(rng)+'rg_'+str(fnumber).zfill(4)+'.png', bbox_inches='tight')
    if choice == 1:
        plt.savefig('iBz_xz'+str(deg)+'_'+str(rng)+'rg_'+str(fnumber).zfill(4)+'.png', bbox_inches='tight')
    '''
    """
################################################################################
    #reinterpolate for magnetic field ratio in disk compared to BH - Section 4
    """ This section creates color plots for Bp_slice and Bzfake1 values
        Bp_slice has a built-in floor.
        Bzfake1 is more detailed and you can use an adjustable ceiling when calling plt.color()
    """
################################################################################
    """
    imyx=reinterpxy(myx,extent,ncell,domask=0,interporder='linear')
    imyy=reinterpxy(myy,extent,ncell,domask=0,interporder='linear')

    iBz=reinterpxy(Bznorm,extent,ncell,domask=1,interporder='linear')
    iBr=reinterpxyhor(Bznorm,extent,ncell,domask=1, interporder='linear') #Brnorm is spherical, BRnorm is cylindrical; switching to Bznorm to keep transition to the disk continuous, has sign info
    ahor=2.0*np.pi*(a**2+3*rhor**2)/3 #surface area of half of the horizon
    BzH=5.0/ahor #Magnetic field per unit area on BH horizon, calculated from Upsilon value in Avara 2015
    Bzfake=iBz*np.sqrt(imyx**2+imyy**2)/(rhor*np.sqrt(5.75)) #vertical magnetic field in the disk per unit area with radial dependence correction, Mdot=5.75 as in Avara 2015
    Br=iBr/np.sqrt(5.75) #radial magnetic field on the horizon per unit area, Mdot=5.75 as in Avara 2015
    Bzfake[iBz.mask==True]=Br[iBz.mask==True]#replace values that are inside the black hole in the equatorial plane with those in the upper half of the horizon
    Bzfake1=Bzfake/BzH
    Bzfake1=np.abs(Bzfake1)
    #create an unnormalized probability function for B/sqrt(Mdot) to better choose relevant seedpoints; Bprob is 1D array from get_coordsw() while Bp_slice is 2D for python frames
    Bp_slice=np.copy(Bzfake1)
    bmin=0.1 #value of Bz that I want to set to 0 in probability distribution
    bmax=np.max(Bp_slice) #value where the probability distribution goes to 1; value before 7/16 was 1.0
    sl=1/(bmax-bmin)
    yint=-bmin*sl
    Bp_slice=sl*Bp_slice+yint #linear probability distribution between Bz=0.1 and Bz=1.0
    Bp_slice[Bp_slice>1.0]=1.0
    Bp_slice[Bp_slice<0.0]=0.0

    grid=np.linspace(-rng, rng, ncell)
    #Make plots to use as movie frames
    plt.figure(1)

    plt.pcolor(grid,grid,Bzfake1,vmax=1.5,cmap='jet') #Change vmax for ceiling modification

    cbar = plt.colorbar()
    # if choice == 1:
    #     cbar.ax.set_yticklabels(['$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^1$','$10^2$'])
    # plt.scatter(cs[:,0],cs[:,1],color='fuchsia')
    # plt.scatter(hc[:,0],hc[:,1],color='yellow')
    horcir=plt.Circle((0,0),2,color='black',fill=True)
    plt.gcf().gca().add_artist(horcir)
    plt.xlim(-rng,rng)
    plt.ylim(-rng,rng)
    plt.xlabel(r"$x [r_g]$",ha='center',labelpad=0,fontsize=14)
    plt.ylabel(r"$y [r_g]$",ha='left',labelpad=0,fontsize=14)
    #plt.title('Inner Region '+str(fnumber).zfill(4))

    plt.savefig('Bzfake_'+str(rng)+'rg_'+str(fnumber).zfill(4)+'.png')
    """

def make_zoomvideo():
    ''' Written by Calen Carron
    A function to call zoomvideo recursively '''
    fnstart=10660
    fnend=10740
    rng = 30
    ncell = 800
    deg = 0
    choice = 1
    for fnumber in range(fnstart, fnend + 1):
        CCzoomvideo(fnumber,rng,ncell,deg,choice)
    print('********** ********** THE PLOTS ARE FINISHED ********** **********')

################################################################################
# Plasma fluid velocity vectors
# started by Karina, modified by Calen
""" The first section of the fluid_vector function can be used to create plots of
    specified rings of vectors (specify in 'vector_rings'). Run zoomvideo first
    and comment out plt.clf() to plot the rings on top of a color plot.

    The second section of the fluid_vector function plots a field of vectors,
    scaled so that general trends in the motion of the plasma can be observed.
    This is best plotted without a color plot underneath.
    Make sure to use plt.clf() so that vector fields don't stack in
    consecutive plots.
"""
################################################################################
def scale_vector(vector_array, scale):
    """ Helper function for fluid_vector """
    new_array = []
    for v in vector_array:
        new_array = new_array + [v*scale]
    return new_array

def fluid_vector(fnumber,rng_limit):
    #plt.clf() ### comment this out when using section one with color underlay!!! ###
    grid3d("gdump.bin",use2d=True)
    #load the fieldline file for a given time and compute standard quantities
    rfd("fieldline"+str(fnumber)+".bin")
    cvel()
    rhor=1+(1-a**2)**0.5
    ihor=np.floor(iofr(rhor)+0.5)
    #compute the 3-velocities in the equatorial slice
    vr = dxdxp[1,1]*uu[1]/uu[0]+dxdxp[1,2]*uu[2]/uu[0]
    vh = dxdxp[2,1]*uu[1]/uu[0]+dxdxp[2,2]*uu[2]/uu[0]
    vp = uu[3]/uu[0]*dxdxp[3,3]
    #
    vrnorm=vr
    vhnorm=vh*np.abs(r)
    vpnorm=vp*np.abs(r*np.sin(h))
    #
    vznorm=vrnorm*np.cos(h)-vhnorm*np.sin(h)
    vRnorm=vrnorm*np.sin(h)+vhnorm*np.cos(h)
    vxnorm=vRnorm*np.cos(ph)-vpnorm*np.sin(ph)
    vynorm=vRnorm*np.sin(ph)+vpnorm*np.cos(ph)

    nxin=int(iofr(2.5))   #2.5 or 5
    nxout=int(iofr(rng_limit))
    nx = int(header[1])
    ny = int(header[2])
    my_x = vxnorm[nxin:nxout, int(ny/2),:]
    my_y = vynorm[nxin:nxout, int(ny/2),:]

    # Get grid ready
    #make the spherical polar grid fully 3D
    myr3d=mk2d3d(r)
    myh3d=mk2d3d(h)
    myph3d=mk2d3d(ph)
    for k in range(0,nz):
        myph3d[:,:,k]=(0.5+k)*2*np.pi/nz
    #compute cartesian coordinates for each grid point of the array
    myx=myr3d*np.sin(myh3d)*np.cos(myph3d)
    myy=myr3d*np.sin(myh3d)*np.sin(myph3d)
    x_coors = myx[nxin:nxout, int(ny/2),:]
    y_coors = myy[nxin:nxout, int(ny/2),:]
    num_rows = x_coors.shape[0]
    num_cols = x_coors.shape[1]


    ############################################################################
    # First Section
    # Specified Vector Rings: manually select which rings you want displayed
    ############################################################################

    # Radial distances for rings of velocity vectors
    vector_rings = [6,9,12,15,18,21,24,27,30] ### MANUALLY ENTER THESE ###
    if min(vector_rings) < 2.6:
        print("\n\n\nERROR: vector ring radius too small, enter number greater than 2.6\n\n")
        return

    # Create each desired vector ring
    for v_rad in vector_rings:
        # Create empty arrays for vector parameters
        ring_x_coors = []
        ring_y_coors = []
        ring_my_x = []
        ring_my_y = []
        leniency = v_rad * 0.016 #might need adjustments based on rad values -- uncomment below to print rad values

        # Iterate through coordinates and select those with the leniency of v_rad
        for i in range(num_rows):
            for j in range(num_cols):
                rad = sqrt((x_coors[i][j])*(x_coors[i][j]) + (y_coors[i][j])*(y_coors[i][j]))
                if (v_rad-leniency < rad < v_rad+leniency):
                    ring_x_coors = ring_x_coors + [x_coors[i][j]]
                    ring_y_coors = ring_y_coors + [y_coors[i][j]]
                    ring_my_x = ring_my_x + [my_x[i][j]]
                    ring_my_y = ring_my_y + [my_y[i][j]]

        # Add the vector ring to the plot
        size = len(ring_x_coors)
        # Calen's secret formula to make pretty vector rings
        if rng_limit <= 10:
            step = int(20 - sqrt(5*v_rad)) #for rng <= 10
        else:
            step = int(((20 - sqrt(10*v_rad))**2)/10)+6 #for rng > 10
        sub_x_coors = ring_x_coors[0:size:step]
        sub_y_coors = ring_y_coors[0:size:step]
        sub_my_x = ring_my_x[0:size:step]
        sub_my_y = ring_my_y[0:size:step]
        plt.quiver(sub_x_coors, sub_y_coors, sub_my_x, sub_my_y, width=0.005, color='white')

    plt.savefig('fluid_vector_'+str(rng_limit)+'rg_'+str(fnumber)+'.png')
    # plt.savefig(os.getcwd()+'/saved_files/fluid_vector/fluid_vector_'+str(fnumber)+'.png')
    return


    ############################################################################
    # Second Section
    # Display all vectors within rng_limit - use without color plot
    ############################################################################
"""
    # makes plot square
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')

    # Iterate through coordinates and select those within the bounds
    all_radii = [0]
    # 3-30 edit to reduce num of vectors
    for i in range(0,num_rows,7):
        for j in range(0,num_cols,7):
            rad = sqrt((x_coors[i][j])*(x_coors[i][j]) + (y_coors[i][j])*(y_coors[i][j]))
            last_rad = all_radii[-1]
            approx = 0.0001*last_rad
            # print(sqrt((x_coors[i][j])*(x_coors[i][j]) + (y_coors[i][j])*(y_coors[i][j])))
            if (last_rad-approx) < rad < (last_rad+approx):
                r_x_coors = r_x_coors + [x_coors[i][j]]
                r_y_coors = r_y_coors + [y_coors[i][j]]
                r_my_x = r_my_x + [my_x[i][j]]
                r_my_y = r_my_y + [my_y[i][j]]
            else:
                if not (all_radii == [0]):
                    ratio = 0.0000001*rng_limit # for 40rg => 0.000004, for 20rg => 0.000002, might need adjusting
                    # 3-30 edit: changing the width from /10 to /5
                    curr_width = ratio*rad*rad*(41-rng_limit)/5 # if the biggest rng_limit is beyond 40, the 41 value needs to be rng_max + 1
                    curr_scale = 1
                    r_my_xn = scale_vector(r_my_x, curr_scale)
                    r_my_yn = scale_vector(r_my_y, curr_scale)
                    plt.quiver(r_x_coors, r_y_coors, r_my_xn, r_my_yn, width=curr_width) # 3-30 edit: took out headwidth=....

                r_x_coors = [x_coors[i][j]]
                r_y_coors = [y_coors[i][j]]
                r_my_x = [my_x[i][j]]
                r_my_y = [my_y[i][j]]
                all_radii = all_radii + [rad]

    horcir=plt.Circle((0,0),2,color='white',fill=True)
    plt.gcf().gca().add_artist(horcir)
    plt.xlim(-rng_limit,rng_limit)
    plt.ylim(-rng_limit,rng_limit)
    plt.xlabel(r"$x [r_g]$",ha='center',labelpad=0,fontsize=14)
    plt.ylabel(r"$y [r_g]$",ha='left',labelpad=0,fontsize=14)
    plt.title(str(fnumber).zfill(4))
    plt.savefig('full_vector_'+str(rng_limit)+'rg_'+str(fnumber)+'.png')
    # plt.savefig(os.getcwd()+'/saved_files/full_vector/full_vector_'+str(fnumber)+'.png')
    return
"""

################################################################################
# Runs fluid vector for a given range of field line files
################################################################################
def fluid_vector_rings(fn_start, fn_end, rng_limit):
    """ Use this function to run fluid_vector section 1 over a certain range of
        field line files and at a specified max radius (rng_limit).
        Note: be sure to comment out 'plt.clf()' at the beginning of fluid_vector """
    for fn in range(fn_start, fn_end+1):
        print("Making file " + str(fn))
        CCzoomvideo(fn, rng_limit)
        fluid_vector(fn, rng_limit)

def full_fluid_vector(fn_start, fn_end, rng_limit):
    """ Use this function to run fluid_vector section 2 over a certain range of
        field line files and at a specified max radius (rng_limit).
        Note: be sure to call plt.clf() in fluid_vector between each plot """
    for fn in range(fn_start, fn_end+1):
        print("Making file " + str(fn))
        fluid_vector(fn, rng_limit)

    
