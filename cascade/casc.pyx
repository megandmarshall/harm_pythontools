from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport log
from libc.math cimport exp


DTYPE = np.float
ctypedef np.float_t DTYPE_t


def fg_p( Eg not None, Ee not None, SeedPhoton seed not None):
    return fgvec( Eg, Ee, seed )

@cython.boundscheck(False) # turn off bounds-checking for entire function
cdef np.ndarray[double, ndim=1] fgvec( np.ndarray[double, ndim=1] Eg, np.ndarray[double, ndim=1] Ee, SeedPhoton seed):
    cdef int i
    cdef int dim = Ee.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=1] Eg1 = np.zeros_like(Ee)
    for i from 0 <= i < dim:
        Eg1[i] = fg( Eg[i], Ee[i], seed )
    return( Eg1 )

cdef double fg( double Eg, double Ee, SeedPhoton seed):
    cdef double Ep = Ee-Eg
    cdef double fgval = ( (seed.f(Eg/(2*Ee*Ep))/(2*Ep**2)) if (Ep>0 and Ee>0 and Eg>0) else (0) )
    return( fgval )

cdef double K1( double Enew, double Eold, SeedPhoton seed ):
    cdef double K = (4*fg(2*Enew,Eold,seed)) if (2*Enew>=seed.Egmin) else (0)
    return( K )

cdef double K2( double Enew, double Eold, SeedPhoton seed ):
    cdef double K = fg(Eold-Enew,Eold,seed)
    return( K )

@cython.boundscheck(False) # turn off bounds-checking for entire function
cdef public double* get_data( np.ndarray[double, ndim=1] nparray ):
    return <double *>nparray.data

def flnew( Evec not None, flold not None, seed not None ):
    return flnew_c( Evec, flold, seed )

@cython.boundscheck(False) # turn off bounds-checking for entire function
cdef public np.ndarray[double, ndim=1] flnew_c( Grid grid, np.ndarray[double, ndim=1] flold, SeedPhoton seed ):
    """Expect E and flold defined on a regular log grid, Evec"""
    cdef int i
    cdef int j
    cdef np.ndarray[double, ndim=1] flnew = np.zeros_like(flold)
    cdef double *flnew_data = get_data(flnew)
    cdef double *Evec_data = get_data(grid.Egrid)
    cdef double *flold_data = get_data(flold)
    cdef int dim = flnew.shape[0]
    cdef Grid newgrid = Grid.empty(dim)
    cdef Func newfunc = Func.fromGrid(newgrid)
    
    for i from 0 <= i < dim:
        #flnew_data[i] = 0
        for j from 0 <= j < dim:
            flnew_data[i] += K1(Evec_data[i],Evec_data[j],seed)*(flold_data[j]*Evec_data[j])
            flnew_data[i] += K2(Evec_data[i],Evec_data[j],seed)*(flold_data[j]*Evec_data[j])
        flnew_data[i] *= grid.dx
    return( flnew )


###############################
#
#  CLASSES
#
###############################        
    

###############################
#
#  SEED PHOTON
#
###############################        

cdef public class SeedPhoton [object CSeedPhoton, type TSeedPhoton ]:
    """our seed photon class"""
    cdef public double Emin
    cdef public double Emax
    cdef public double s
    cdef public double Egmin
    cdef public double Nprefactor

    def __init__(self, double Emin, double Emax, double s):
        self.Emin = Emin
        self.Emax = Emax
        self.s = s
        #minimum energy gamma-ray to be able to pair produce
        self.Egmin = 2./Emax
        self.Nprefactor = (1.-s)/(Emax**(1.-s)-Emin**(1.-s))

    cpdef int canPairProduce(self, double E):
        return( E > self.Egmin )

    cpdef double f(self, double E):
        return( self.Nprefactor*E**(-self.s) if (E >= self.Emin and E <= self.Emax) else 0 )


###############################
#
#  GRID
#
###############################        

cdef public class Grid [object CGrid, type TGrid ]:
    """grid class"""
    cdef public double Emin
    cdef public double Emax
    cdef public double E0
    cdef public double xmin
    cdef public double xmax
    cdef public int Ngrid
    cdef public double dx
    cdef public Egrid
    cdef public xgrid
    cdef public dEdxgrid
    cdef double *xgrid_data
    cdef double *Egrid_data
    cdef double *dEdxgrid_data

    def __init__(self, double Emin, double Emax, double E0, int Ngrid):
        """ Full constructor: allocates memory and generates the grid """
        self.Ngrid = Ngrid
        self.xgrid = np.zeros((self.Ngrid),dtype=DTYPE)
        self.Egrid = np.zeros((self.Ngrid),dtype=DTYPE)
        self.dEdxgrid = np.zeros((self.Ngrid),dtype=DTYPE)
        self.set_grid( Emin, Emax, E0 )

    @classmethod
    def fromGrid(cls, Grid grid):
        return cls( grid.Emin, grid.Emax, grid.E0, grid.Ngrid )

    @classmethod
    def empty(cls, int Ngrid):
        return cls( 0, 1, 0.5, Ngrid )

    @cython.boundscheck(False) # turn off bounds-checking for entire function
    cpdef set_grid(self, double Emin, double Emax, double E0 ):
        """ Same as Grid() but without reallocation of memory """
        cdef int i
        cdef int dim = self.Ngrid
        self.Emin = Emin
        self.Emax = Emax
        self.E0 = E0
        self.xmax = log(self.Emax-self.E0)
        self.xmin = log(self.Emin-self.E0)
        self.dx = (self.xmax - self.xmin) / (dim-1)
        #get direct C pointers to numpy arrays' data fields
        self.xgrid_data = get_data(self.xgrid)
        self.Egrid_data = get_data(self.Egrid)
        self.dEdxgrid_data = get_data(self.dEdxgrid)
        for i from 0 <= i < dim:
            self.xgrid_data[i] = self.xmin + self.dx*i
            self.Egrid_data[i] = self.E0 + exp(self.xgrid_data[i])
            self.dEdxgrid_data[i] = self.Egrid_data[i] - self.E0

    @cython.boundscheck(False) # turn off bounds-checking for entire function
    cpdef int iofx(self, double xval):
        """ Returns the index of the cell containing xval """
        cdef int ival
        ival = int( (self.xval-self.xmin)/self.dx )
        return ival

    @cython.boundscheck(False) # turn off bounds-checking for entire function
    cpdef double xofE(self, double Eval):
        """ Returns the value of x corresponding to Eval """
        cdef double xval
        xval = log(Eval - self.E0)
        return xval

    @cython.boundscheck(False) # turn off bounds-checking for entire function
    cpdef int iofE(self, double Eval):
        """ Returns the index of the cell containing Eval """
        return self.iofx( self.xofE(Eval) )


###############################
#
#  FUNCTION
#
###############################        

cdef public class Func(Grid)  [object CFunc, type TFunc ]:
    """ Function class derived from Grid class """
    
    cdef public func_vec
    cdef double *func_vec_data

    def __init__(self, double Emin, double Emax, double E0, int Ngrid, func_vec = None):
        Grid.__init__(self, Emin, Emax, E0, Ngrid)
        if func_vec is None:
            self.func_vec = np.zeros((self.Ngrid),dtype=DTYPE)
            self.func_vec_data = get_data(self.func_vec)
        else:
            self.func_vec = np.copy(func_vec)
            self.func_vec_data = get_data(self.func_vec)

    cpdef set_grid(self, double Emin, double Emax, double E0):
        """ Same as Grid() but without reallocation of memory """
        Grid.set_grid( self, Emin, Emax, E0 )

    @classmethod
    def fromGrid(cls, Grid grid):
        return cls( grid.Emin, grid.Emax, grid.E0, grid.Ngrid )

    @classmethod
    def empty(cls, int Ngrid):
        return cls( 0, 1, 0.5, Ngrid )

    @cython.boundscheck(False) # turn off bounds-checking for entire function
    cpdef double fofE(self, double Eval):
        """ Linearly interpolates f(E) in log-log """
        cdef int i = Grid.iofE( self, Eval )
        cdef double logfl, logfr, logxl, logxr, logf, f
        if i < 0:
            return self.func_vec_data[0]
        if i >= Grid.Ngrid-1:
            return self.func_vec_data[Grid.Ngrid-1]
        logx  = log(Eval)
        logxl = log(self.Evec[i])
        logxr = log(self.Evec[i+1])
        logfl = log(self.func_vec_data[i])
        logfr = log(self.func_vec_data[i+1])
        logf  = (logfr * (logx - logxl) + logfl * (logx - logxr)) / (logxr - logxl)
        f = exp(logf)
        return( f )
        
    def set_func(self, func_vec):
        self.set_func_c( get_data(func_vec) )

    @cython.boundscheck(False) # turn off bounds-checking for entire function
    cdef set_func_c(self, double *func_vec_data):
        cdef int i
        for i from 0 <= i < Grid.Ngrid:
            self.func_vec_data[i] = func_vec_data[i]


