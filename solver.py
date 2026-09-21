#solver 
import numpy as np 
#FCDS = finite central difference scheme 

def FCDS(T_0,pe_cell,N):
    """
    SUMMARY.
    solves the 1D convection diffusion equation
    rho u dt/dx = gamma d2T/dx2
    using : 
        -finite central difference scheme for convection term 
        -finite central difference scheme for diffusive term 
        
    PARAMETERS: 
        ----
        T_0 : array    -initial condition 
        pe_cell : array  -cell peclet number
        N: integer   -number of nodes
        
    Returns
    T :array   -Numerical solution for 1D temperature distribution 
    -------
    None
    """
    #parameters for convergence criterion
    num_err=1. #numerical error
    epsilon=1.E-8
    #initialise temperature 
    T=T_0.copy()
    #initialise new temperature
    T_new=T_0.copy()
    while num_err >epsilon:
        #perform iterations over all interior grid points 
        for i in range(1,N-1):
            T_new[i] = (0.5 - 0.25*pe_cell)*T[i+1] + (0.5 + 0.25*pe_cell)*T[i-1]
        
        #num_err = 0

        #for j in range(1, N-1):
             #num_err = num_err + np.abs(T_new[j]-T[j])     
        
        
        #calculate numerical error 
        num_err=np.max(
            np.abs(T_new[1:N-1]-T[1:N-1])
            )
        T=T_new.copy()
        
    return T 

def FBDS1(T_0,pe_cell,N):
    """
    SUMMARY.
    solves the 1D convection diffusion equation
    rho u dt/dx = gamma d2T/dx2
    using : 
        -finite Bacward difference 1st order scheme for the convective term  
        -finite central difference scheme for diffusive term 
        
    PARAMETERS: 
        ----
        T_0 : array    -initial condition 
        pe_cell : array  -cell peclet number
        N: integer   -number of nodes
        
    Returns
    T :array   -Numerical solution for 1D temperature distribution 
    -------
    None
    """
    num_err=1.
    epsilon=1.E-8
    #initialise the temperature
    T=T_0.copy()
    #initialise new temperature
    T_new=T_0.copy()
    while num_err > epsilon:
        #perform iteration over all grid points 
        for i in range(1,N-1):
            T_new[i]=(T[i+1]+T[i-1]*(1+pe_cell))/(2+pe_cell)
            
        #calculate numerical error 
        num_err=np.max(
            np.abs(T_new[1:N-1]-T[1:N-1])
            )
        T=T_new.copy()    
    return T


def FBDS2(T_0,pe_cell,N):
    """
    SUMMARY.
    solves the 1D convection diffusion equation
    rho u dt/dx = gamma d2T/dx2
    using : 
        -finite Bacward difference 2ND order scheme for the convective term  
        -finite central difference scheme for diffusive term 
        
    PARAMETERS: 
        ----
        T_0 : array    -initial condition 
        pe_cell : array  -cell peclet number
        N: integer   -number of nodes
        
    Returns
    T :array   -Numerical solution for 1D temperature distribution 
    -------
    None
    """
    num_err=1.
    epsilon=1.E-8
    #initialise the temperature
    T=T_0.copy()
    #initialise new temperature
    T_new=T_0.copy()
    while num_err > epsilon:
        #perform iteration over all grid points 
        T_new[1]=((1+pe_cell)*T[0]+T[2])/(2+pe_cell)
        for i in range(2,N-1):
            T_new[i]=(2*T[i+1]+(2+4*pe_cell)*T[i-1]-pe_cell*T[i-2])/(3*pe_cell + 4 )
        #calculate numerical error 
        num_err=np.max(
            np.abs(T_new[1:N-1]-T[1:N-1])
            )
        T=T_new.copy()    
    return T
