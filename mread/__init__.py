from matplotlib import rc
#rc('verbose', level='debug')
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('mathtext',fontset='cm')
#rc('mathtext',rm='stix')
#rc('text', usetex=True)

#from pylab import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show

import numpy as np
import array
import scipy as sc
from scipy.interpolate import griddata
#from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy import ma
import matplotlib.colors as colors
import os,glob


#global rho, ug, vu, uu, B, CS
#global nx,ny,nz,_dx1,_dx2_,_dx3,ti,tj,tk,x1,x2,x3,r,h,ph,gdet,conn,gn3,gv3,ck,dxdxp

def plc(myvar): #plc
    global x1, x2
    res = plt.contour(x1[:,:,0],x2[:,:,0],myvar[:,:,0],15)
    plt.colorbar()

def mainfunc(imgname):
    global xi,yi,zi,CS
    #grid3d("gdump")
    #jrdp3d("fieldline0250.bin")
    xraw=r*np.sin(h)
    yraw=r*np.cos(h)
    lrhoraw=np.log10(rho)
    x=xraw[:,:,0].view().reshape(-1)
    y=yraw[:,:,0].view().reshape(-1)
    lrho=lrhoraw[:,:,0].view().reshape(-1)
    #mirror
    x=np.concatenate((-x,x))
    y=np.concatenate((y,y))
    lrho=np.concatenate((lrho,lrho))
    extent=(-41,41,-41,41)
    # define grid.
    xi = np.linspace(-41.0, 41.0, 800)
    yi = np.linspace(-41.0, 41.0, 800)
    # grid the data.
    zi = griddata((x, y), lrho, (xi[None,:], yi[:,None]), method='cubic')
    interior = np.sqrt((xi[None,:]**2) + (yi[:,None]**2)) < 1+np.sqrt(1-a**2)
    #zi[interior] = np.ma.masked
    zim = ma.masked_where(interior, zi)
    palette=cm.jet
    palette.set_bad('k', 1.0)
    palette.set_over('r', 1.0)
    palette.set_under('g', 1.0)
    # contour the gridded data, plotting dots at the randomly spaced data points.
    cset2 = plt.contour(zi,15,linewidths=0.5,colors='k', extent=extent,hold='on',origin='lower')
    #for c in cset2.collections:
    #    c.set_linestyle('solid')
    #CS = plt.contourf(xi,yi,zi,15,cmap=palette)
    #CS = plt.imshow(zim, extent=[0.01,80,-40,40], cmap = palette, norm = colors.Normalize(vmin=-1,vmax=-0.2,clip = False))
    CS = plt.imshow(zim, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower')
    #CS.cmap=cm.jet
    #CS.set_axis_bgcolor("#bdb76b")
    plt.colorbar(CS) # draw colorbar
    plt.xlim(-40,40)
    plt.ylim(-40,40)
    plt.title('t = %f' % t)
    #plt.show()
    # rbf = Rbf(x[0:288:8,:,0].view().reshape(-1),y[0:288:8,:,0].view().reshape(-1),rho[0:288:8,:,0].view().reshape(-1),epsilon=2)
    # ZI = rbf( XI, YI )
    # # plot the result
    # n = plt.normalize(0.0, 40.0)
    # plt.subplot(1, 1, 1)
    # plt.imshow(XI, YI, ZI, cmap=cm.jet)
    # #plt.scatter(x, y, 100, z, cmap=cm.jet)
    # plt.title('RBF interpolation - multiquadrics')
    # plt.xlim(0, 40.0)
    # plt.ylim(0, 40.0)
    # plt.colorbar()
    # plt.figure()
    # plt.plot(r[:,0,0],np.log10(absflux),'b')
    # plt.legend(['Grid'])
    # plt.axis([r[0,0,0],100,-5,5])
    # plt.title('Grid plot')
    # plt.show()

def jrdp3d(fieldlinefilename):
    global t,nx,ny,nz,_dx1,_dx2,_dx3,gam,a,rho,ug,uu,uut,uu,B,uux
    #read image
    fin = open( "dumps/" + fieldlinefilename, "rb" )
    header = fin.readline().split()
    t = np.float64(header[0])
    nx = int(header[1])
    ny = int(header[2])
    nz = int(header[3])
    _dx1=float(header[7])
    _dx2=float(header[8])
    _dx3=float(header[9])
    gam=float(header[11])
    a=float(header[12])
    body = np.fromfile(fin,dtype=np.single,count=nx*ny*nz*11)
    fin.close()
    d=body.view().reshape((-1,nx,ny,nz),order='F')
    #rho, u, -hu_t, -T^t_t/U0, u^t, v1,v2,v3,B1,B2,B3
    rho=d[0,:,:,:]
    ug=d[1,:,:,:]
    uu=d[4:8,:,:,:]  #note uu[i] are 3-velocities (as read from the fieldline file)
    #uut=np.copy(d[4,:,:,:])
    #multiply by u^t to get 4-velocities: u^i = u^t v^i
    #uux=np.copy(uu)
    #for i in range(1,4):
    #    uux[i,:,:,:] = uux[i,:,:,:] * uux[0,:,:,:]  
    uu[1:4]=uu[1:4] * uu[0]
    B = sc.zeros_like(uu)
    B[1:4,:,:,:]=d[8:11,:,:,:]


def cvel():
    global ud,etad, etau, gamma, vu, vd, bu, bd, bsq
    ud = mdot(gv3,uu)                  #g_mn u^n
    etad = sc.zeros_like(uu)
    etad[0] = -1/(-gn3[0,0])**0.5      #ZAMO frame velocity (definition)
    etau = mdot(gn3,etad)
    gamma=mdot(uu,etad)                #Lorentz factor as measured by ZAMO
    vu = uu - gamma*etau               #u^m = v^m + gamma eta^m
    vd = mdot(gv3,vu)
    bu=sc.empty_like(uu)              #allocate memory for bu
    #set component per component
    bu[0]=mdot(B[1:4], ud[1:4])             #B^i u_i
    bu[1:4]=(B[1:4] + bu[0]*uu[1:4])/uu[0]  #b^i = (B^i + b^t u^i)/u^t
    bd=mdot(gv3,bu)
    bsq=mdot(bu,bd)



def grid3d(dumpname): #read gdump: header and body
    global nx,ny,nz,_dx1,_dx2_,_dx3,ti,tj,tk,x1,x2,x3,r,h,ph,conn,gn3,gv3,ck,dxdxp,gdet
    print( "Reading grid from " + "dumps/" + dumpname + " ..." )
    gin = open( "dumps/" + dumpname, "rb" )
    header = gin.readline().split()
    nx = int(header[1])
    ny = int(header[2])
    nz = int(header[3])
    _dx1=float(header[7])
    _dx2=float(header[8])
    _dx3=float(header[9])
    gin.close()
    #read gdump
    
    gd = np.loadtxt( "dumps/gdump", 
                      dtype=float, 
                      skiprows=1, 
                      unpack = True ).view().reshape((-1,nx,ny,nz), order='F')
    #gd = np.genfromtxt( "dumps/gdump", 
    #                 dtype=float, 
    #                 skip_header=1, 
    #                 skip_footer=nx*ny*(nz-1),
    #                 unpack = True ).view().reshape((137,nx,ny,nz), order='F')
    ti,tj,tk,x1,x2,x3,r,h,ph = gd[0:9,:,:,:].view() 
    #get the right order of indices by reversing the order of indices i,j(,k)
    #conn=gd[9:73].view().reshape((4,4,4,nx,ny,nz), order='F').transpose(2,1,0,3,4,5)
    gn3 = gd[73:89].view().reshape((4,4,nx,ny,nz), order='F').transpose(1,0,2,3,4)
    gv3 = gd[89:105].view().reshape((4,4,nx,ny,nz), order='F').transpose(1,0,2,3,4)
    gdet = gd[105]
    ck = gd[106:110].view().reshape((4,nx,ny,nz), order='F')
    dxdxp = gd[110:136].view().reshape((4,4,nx,ny,nz), order='F').transpose(1,0,2,3,4)
    print( "Done!" )

def mdot(a,b):
    """
    Computes a contraction of two tensors/vectors.  Assumes
    the following structure: tensor[m,n,i,j,k] OR vector[m,i,j,k], 
    where i,j,k are spatial indices and m,n are variable indices. 
    """
    if a.ndim == 4 and b.ndim == 4:
          c = (a*b).sum(0)
    elif a.ndim == 5 and b.ndim == 4:
          c = sc.empty_like(a[:,0,:,:,:])      
          for i in range(a.shape[0]):
                c[i,:,:,:] = (a[i,:,:,:,:]*b).sum(0)
    elif a.ndim == 4 and b.ndim == 5:
          c = sc.empty_like(b[0,:,:,:,:])      
          for i in range(b.shape[1]):
                c[i,:,:,:] = (a*b[:,i,:,:,:]).sum(0)
    elif a.ndim == 5 and b.ndim == 5:
          c = sc.empty((a.shape[0],b.shape[1],a.shape[2],a.shape[3],a.shape[4]))
          for i in range(c.shape[0]):
                for j in range(c.shape[1]):
                      c[i,j,:,:,:] = (a[i,:,:,:,:]*b[:,j,:,:,:]).sum(0)
    else:
           raise Exception('mdot', 'wrong dimensions')
    return c

def fieldcalc():
    """
    Computes the field vector potential
    """
    daphi = (gdet*B[1]).sum(2)*_dx2*_dx3
    aphi=daphi.cumsum(axis=1)
    return(aphi[:,:,None])

def fieldcalc2():
    """
    Computes the field vector potential
    """
    #daphi = (gdet*B[1]).sum(2)*_dx2*_dx3
    #aphi=daphi.cumsum(axis=1)
    daphi = (-gdet*B[2]).sum(2)*_dx2*_dx3
    #sc.zeros_like(B[1,:,:,0:1])
    daphi1 = (gdet[0]*B[1,0]).sum(1).cumsum(axis=0)*_dx2*_dx3
    daphi[0,:] += daphi1
    aphi=daphi.cumsum(axis=0)
    return(aphi[:,:,None])

def horfluxcalc(ihor):
    """
    Computes the absolute flux through the sphere i = ihor
    """
    #1D function of theta only:
    dfabs = (gdet[ihor]*np.abs(B[1,ihor])).sum(1)*_dx2*_dx3
    fabs = dfabs.sum(axis=0)
    return(fabs)

def mdotcalc(ihor):
    """
    Computes the absolute flux through the sphere i = ihor
    """
    #1D function of theta only:
    mdot = (-gdet[ihor]*rho[ihor]*uu[1,ihor]).sum()*_dx2*_dx3
    return(mdot)


def diskfluxcalc(jmid):
    """
    Computes the absolute flux through the disk midplane at j = jmid
    """
    #1D function of theta only:
    dfabs = (gdet[:,jmid,:]*np.abs(B[2,:,jmid,:])).sum(1)*_dx1*_dx3
    fabs = dfabs.sum(axis=0)
    return(fabs)

def fhorvstime(ihor):
    """
    Returns a tuple (ts,fs,mdot): lists of times, horizon fluxes, and Mdot
    """
    flist = glob.glob( os.path.join("dumps/", "fieldline*.bin") )
    ts=np.empty(len(flist),dtype=float)
    fs=np.empty(len(flist),dtype=float)
    mdot=np.empty(len(flist),dtype=float)
    for findex, fname in enumerate(flist):
        print( "Reading " + fname + " ..." )
        jrdp3d("../"+fname)
        fs[findex]=horfluxcalc(ihor)
        mdot[findex]=mdotcalc(ihor)
        ts[findex]=t
    print( "Done!" )
    return((ts,fs,mdot))

def plotit(ts,fs,md):
    #rc('font', family='serif')
    #plt.figure( figsize=(12,9) )
    fig,plotlist=plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(12,9))
    plotlist[0].plot(ts,fs,label=r'$\Phi_{\rm h}/0.5\Phi_{\rm i}$: Normalized Horizon Magnetic Flux')
    plotlist[0].plot(ts,fs,'r+', label=r'$\Phi_{\rm h}/0.5\Phi_{\rm i}$: Data Points')
    plotlist[0].legend(loc='lower right')
    #plt.xlabel(r'$t\;(GM/c^3)$')
    plotlist[0].set_ylabel(r'$\Phi_{\rm h}$',fontsize=16)
    plt.setp( plotlist[0].get_xticklabels(), visible=False)
    plotlist[0].grid(True)
    #
    #plotlist[1].subplot(212,sharex=True)
    plotlist[1].plot(ts,md,label=r'$\dot M_{\rm h}$: Horizon Accretion Rate')
    plotlist[1].plot(ts,md,'r+', label=r'$\dot M_{\rm h}$: Data Points')
    plotlist[1].legend(loc='upper right')
    plotlist[1].set_xlabel(r'$t\;(GM/c^3)$')
    plotlist[1].set_ylabel(r'$\dot M_{\rm h}$',fontsize=16)

    #title("\TeX\ is Number $\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!", 
    #      fontsize=16, color='r')
    plotlist[1].grid(True)
    fig.savefig('test.pdf')


def test():
    t=np.arange(10)
    f=np.arange(10)**2
    plt.plot(t,f,label='$\Phi$')
    plt.title(r"This is a title")
    plt.legend(loc='upper right')
    plt.xlabel(r'$t (GM/c^3)$')
    plt.ylabel(r'$\Phi_{\rm h}$',fontsize=16)
    #plt.legend()

if __name__ == "__main__":
    import sys
    #mainfunc()
    grid3d("gdump")
    jrdp3d("fieldline0250.bin")
    #cvel()
    #plc(rho)
    #jrdp3d("fieldline0000.bin")
    #diskflux=diskfluxcalc(ny/2)
    #ts,fs,md=fhorvstime(10)
    #plotit(ts,fs/(0.5*diskflux),md)
    #jrdp3d("fieldline0320.bin")
    plt.figure(1)
    aphi=fieldcalc()
    plc(aphi)
    plt.figure(2)
    aphi2=fieldcalc2()
    plc(aphi2)
    #test()


