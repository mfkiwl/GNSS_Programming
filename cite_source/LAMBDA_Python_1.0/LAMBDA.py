"""
###############################################################################
#                                                                             #
#          Least-squares AMBiguity Decorrelation Adjustment (LAMBDA)          #
#                                                                             #
###############################################################################
#                                                                             #
# The LAMBDA method was introduced by Teunissen (1993) and developed at the   # 
# Delft University of Technology in the 1990s. With the LAMBDA method, the    #
# integer least squares (ILS) ambiguity estimates are computed in two steps.  #
# First the ambiguities are decorrelated, by means of the Z-transformation.   #
# The decorrelation is essential to allow for an efficient search for the     #
# integer candidate which solves the integer minimization problem. This       #
# search is the second step in the LAMBDA method.                             #
#                                                                             #                 
# LAMBDA was first implemented in Fortran-77, which is Version 1.0 of the     #
# LAMBDA software. The implementation aspects are extensively described in    #
# De Jonge and Tiberius (1996). The code was first translated into Matlab by  #
# Borre, see (Strang and Borre 1997), and further improved to Version 2.1 by  #
# Joosten (2001). In 2012 Li and Verhagen implemented the LAMBDA Version 3.0, #
# including the search-and-shrink technique, partial ambiguity resolution and # 
# a discrimination test to decide on acceptance of the fixed solution. This   #
# Python implementation is a replica of the Matlab Version 3.0 aiming to      #
# support users developing their GNSS processing software in Python.          #
#                                                                             #
#-----------------------------------------------------------------------------#
#----------------------------- Python Version 1.0 ----------------------------#
# Release Date  : 01-NOV-2019                                                 #
# Author        : Dimitrios Psychas                                           #
# Affiliation   : Mathematical Geodesy and Positioning,                       #
#                 Delft University of Technology                              #
# Contact       : D.V.Psychas-1@tudelft.nl ; d.psychas@fugro.com              #
#                                                                             #
#----------------------------- Matlab Version 3.0 ----------------------------#
# Release Date  : 01-SEP-2012                                                 #
# Authors       : Bofeng Li and Sandra Verhagen                               #
# Affiliation   : GNSS Research Centre, Curtin University of Technology       #
#                 Mathematical Geodesy and Positioning,                       #
#                 Delft University of Technology                              #
# Contact       : bofeng.li@curtin.edu.au ; A.A.Verhagen@tudelft.nl           #
#                                                                             #
#-----------------------------------------------------------------------------#
#                                                                             #
# REFERENCES:                                                                 #
#   1. LAMBDA Software Package: Matlab implementation, Version 3.0.           #
#      Documentation provided with this software package.                     #
#   2. Teunissen PJG (1993) Least-squares estimation of the integer GPS       #
#      ambiguities. In: Invited lecture, section IV theory and methodology,   #
#      IAG General Meeting, Beijing, China                                    #
#   3. Teunissen PJG (1995) The least-squares ambiguity decorrelation         #
#      adjustment: a method for fast GPS ambiguity estimation. J Geod         #
#      70:651-7                                                               #
#   4. De Jonge P, Tiberius C (1996) The LAMBDA method of integer ambiguity   #
#      estimation: implementation aspects.                                    #
#   5. Chang X, Yang X, Zhou T (2005) MLAMBDA: a modified LAMBDA method for   #
#      integer least-squares estimation                                       #
#   6. Strang G, Borre K (1997) Linear Algebra, Geodesy, and GPS.             #
#                                                                             #
###############################################################################
"""

# Import libraries
import numpy as np
from numpy.linalg import inv
from scipy.stats import norm
from numpy.linalg import eig
import numpy.matlib

def ldldecom(Qahat):
    """ 
    # ----------------------------------------------------------------------- #
    #   ldldecom: Find LtDL-decomposition of Qahat-matrix                     #
    #                                                                         #
    #             L,D = ldldecom(Qahat)                                       #
    #                                                                         #             
    #   This function finds the LtDL decomposition of a given ambiguity       #
    #   variance-covariance matrix.                                           #
    #                                                                         #
    #   Input arguments:                                                      #
    #       Qahat: Variance-covariance matrix of ambiguities to be factored   #
    #              (symmetric n-by-n matrix; 2D numpy array)                  #
    #                                                                         #
    #   Output arguments:                                                     #
    #       L    : n-by-n factor matrix (2D numpy array)                      #
    #              (unit lower triangular)                                    #
    #       D    : Diagonal n-vector (1D numpy array)                         #
    #                                                                         #
    # ----------------------------------------------------------------------- #
    #   Function : ldldecom                                                   #
    #   Date     : 19-MAY-1999                                                #
    #   Author   : Peter Joosten                                              #
    #              Mathematical Geodesy and Positioning,                      #
    #              Delft University of Technology                             #
    # ----------------------------------------------------------------------- #
    """ 
    
    # Find the total number of ambiguities
    n = len(Qahat)
    
    # Initialize the n-by-n lower triangular matrix and the diagonal n-vector
    L    = np.zeros((n,n))
    D    = np.empty((1,n))
    D[:] = np.nan

    # Decorrelate the input vc-matrix starting from its last row
    for i in range(n-1,-1,-1):
    
        D[0,i] = Qahat[i,i].copy()
        
        L[i,0:i+1] = Qahat[i,0:i+1]/np.sqrt(Qahat[i,i])
    
        for j in range(0,i,1):
            Qahat[j,0:j+1] -= L[i,0:j+1]*L[i,j]
        
        L[i,0:i+1] /= L[i,i]
        
    # Terminate in case the input vc-matrix is not positive definite      
    if np.any((D<1e-10)) == True:
        print("Error: the input vc-matrix (in function ldldecom) is not " +  
               "positive definite!")
        raise SystemExit
        
    return L, D

def decorrel(Qahat,ahat):
    """ 
    # ----------------------------------------------------------------------- #
    #   decorrel: Decorrelate a (co)variance matrix of ambiguities            #
    #                                                                         #
    #             Qzhat,Z,L,D,zhat,iZt = decorrel(Qahat,ahat)                 #
    #                                                                         #             
    #   This function creates a decorrelated Q-matrix, by finding the         #
    #   Z-matrix and performing the corresponding transformation.             #
    #                                                                         #
    #   The routine is based on Fortran routines written by Paul de Jonge     #
    #   and on Matlab-routines written by Kai Borre.                          #
    #   The resulting Z-matrix can be used as follows:                        # 
    #       zhat = Z.T.dot(ahat) | \hat{z} = Z^T * \hat{a}                    # 
    #       Qzhat = Z.T.dot(Qahat).dot(Z) | Q_{\hat{z}} = Z^T*Q_{\hat{a}}*Z   #
    #                                                                         #
    #   Input arguments:                                                      #
    #       Qahat: Variance-covariance matrix of ambiguities (original)       #
    #              (symmetric n-by-n matrix; 2D numpy array)                  #
    #       ahat : Original ambiguities (1D numpy array)                      #
    #                                                                         #
    #   Output arguments:                                                     #
    #       Qzhat: Variance-covariance matrix of decorrelated ambiguities     #
    #       Z    : Z-transformation matrix                                    #
    #       L    : L matrix (from LtDL-decomposition of Qzhat)                #
    #       D    : D matrix (from LtDL-decomposition of Qzhat)                #
    #       zhat : Transformed ambiguities                                    #
    #       iZt  : inv(Z')-transformation matrix                              #
    #                                                                         #
    # ----------------------------------------------------------------------- #
    #   Function : decorrel                                                   #
    #   Date     : 19-MAY-1999 / modified 12-APR-2012                         #
    #   Author   : Peter Joosten / Sandra Verhagen                            #
    #              Mathematical Geodesy and Positioning,                      #
    #              Delft University of Technology                             #
    # ----------------------------------------------------------------------- #
    """    
                                                                            
    # Tests on inputs ahat and Qahat
    
    # Is the Qahat symmetric ?
    if ((Qahat-Qahat.T)<1e-8).all()==False:
        print("Variance-covariance matrix of float estimated ambiguities " + \
              "is not symmetric!")
        raise SystemExit
        
    # Is the Qahat positive-definite ?    
    
    # Get the eigenvalues (D) and eigenvectors (V) of Qahat
    D,V = eig(Qahat)
    
    # Check the positive-definiteness    
    if sum(D>0) != len(Qahat):
        print("Variance-covariance matrix of float estimated ambiguities " + \
              "is not positive definite!")
        raise SystemExit   

    # Initialization
    n   = len(Qahat)
    iZt = np.eye(n)
    i1  = n-2
    sw  = 1

    # LtDL decomposition
    L,D = ldldecom(Qahat.copy())

    # Decorrelation procedure
    while sw:
    
        i  = n-1 # loop for column from n-1 to 0
        sw = 0    
    
        while (not sw) and (i > 0):
        
            i-=1 # the i-th column
        
            if i <= i1:
            
                for j in range(i+1,n):
                    
                    mu = np.round(L[j,i])
                
                    if mu != 0.0: # if mu not equal to 0
                        L[j:n,i] -= mu*L[j:n,j]
                        iZt[:,j] += mu*iZt[:,i] # iZt is inv(Zt) matrix

            delta = D[0,i] + L[i+1,i]**2 * D[0,i+1]
        
            if delta < D[0,i+1]:
            
                lamda    = D[0,i+1] * L[i+1,i] / delta
                eta      = D[0,i] / delta
                D[0,i]   = eta*D[0,i+1]
                D[0,i+1] = delta            

                L[i:i+2,0:i] = \
                     np.array([[-L[i+1,i],1],[eta, lamda]]).dot(L[i:i+2,0:i])
                L[i+1,i] = lamda
            
                # Swap rows i and i+1            
                if i==0:
                
                    L[i+2:n,i:i+2:1] = L[i+2:n,i+1::-1].copy()
                    iZt[:,i:i+2:1]   = iZt[:,i+1::-1].copy()
                
                else:
                
                    L[i+2:n,i:i+2:1] = L[i+2:n,i+1:i-1:-1].copy()
                    iZt[:,i:i+2:1]   = iZt[:,i+1:i-1:-1].copy()

                i1 = i
                sw = 1
    
    # Return the transformed Q-matrix and the transformation matrix     
    Z     = np.round(inv(iZt.T))
    Qzhat = Z.T.dot(Qahat).dot(Z)

    # Return the decorrelated ambiguities
    zhat  = Z.T.dot(ahat)                                                                          
    
    return Qzhat,Z,L,D,zhat,iZt
    
def bootstrap(ahat,L):
    """ 
    # ----------------------------------------------------------------------- #
    #   bootstrap: Bootstrapping for integer ambiguity resolution             #
    #                                                                         #
    #              afixed = bootstrap(ahat,L)                                 #
    #                                                                         #
    #   This function performs the boostrapping technique for integer         #
    #   ambiguity resolution.                                                 #
    #                                                                         #
    #   Input arguments:                                                      #
    #       ahat : float ambiguities (1D numpy array)                         #
    #       L    : unit lower triangular matrix from LtDL-decomposition of    #
    #              the vc-matrix                                              #
    #                                                                         #        
    #   Output arguments:                                                     #
    #      afixed: fixed ambiguities (1D numpy array)                         #
    #                                                                         #        
    #   Note:                                                                 #
    #       The success rate with boostrapping depends on the parametrization #
    #       of the float ambiguities ahat. It should be applied to the        #
    #       decorrelated ambiguities. The complete procedure is then:         #
    #           Qzhat,Z,L,D,zhat,iZt = decorrel(Qahat,ahat)                   # 
    #           zfixed = bootstrap(zhat,L)                                    #
    #           afixed = iZt.dot(zfixed)                                      #
    #                                                                         #        
    # ----------------------------------------------------------------------- #
    #   Function : boostrap                                                   #
    #   Date     : 06-SEP-2010                                                #
    #   Author   : Bofeng Li                                                  #
    #              GNSS Research Center, Department of Spatial Sciences,      #
    #              Curtin University of Technology                            #
    # ----------------------------------------------------------------------- #
    """   
                                                                            
    # Number of ambiguities
    n = len(ahat)

    # Initialize vectors to store conditional and fixed ambiguities
    afcond = np.zeros(n,)
    afixed = np.zeros(n,)

    # Initialize vector to compute conditional a
    S = np.zeros(n,)

    # Start from last ambiguity (should preferably be the most precise one)
    afcond[n-1] = ahat[n-1]

    # Rounding of last ambiguity
    afixed[n-1] = round(afcond[n-1])

    for i in range(n-2,-1,-1):
    
        # Compute the i-th cond. ambiguity on the ambiguities from i+1 to n
        S[0:i+1] = S[0:i+1] + (afixed[i+1]-afcond[i+1])*L[i+1,0:i+1]
    
        afcond[i] = ahat[i] + S[i]
        afixed[i] = round(afcond[i])                                                                            
        
    return afixed
    
def SR_IB_1(Qahat):
    """
    # ----------------------------------------------------------------------- #
    #   SR_IB_1: Compute the bootstrapped success rate (function #1)          #
    #                                                                         #
    #            Ps = SR_IB_1(Qahat)                                          # 
    #                                                                         #
    #   This function computes the exact bootstrapped success rate given the  #
    #   variance-covariance matrix of float estimated ambiguities.            #
    #                                                                         #    
    #   Input arguments:                                                      #
    #       Qahat: Variance-covariance matrix of ambiguities                  #
    #              (symmetric n-by-n matrix; 2D numpy array)                  #
    #                                                                         #        
    #   Output arguments:                                                     #
    #       Ps   : Exact bootstrapped success rate                            #
    #                                                                         #        
    # ----------------------------------------------------------------------- #
    #   Function : SR_IB_1                                                    #
    #   Date     : 08-FEB-2019                                                #
    #   Author   : Dimitrios Psychas                                          #
    #              Mathematical Geodesy and Positioning,                      #
    #              Delft University of Technology                             # 
    # ----------------------------------------------------------------------- #
    """     
    
    # Perform the LtDL decomposition                                                                  
    L,D = ldldecom(Qahat)           

    # Compute the exact bootstrapped success rate using the D matrix
    Ps = np.prod(2*norm.cdf(0.5/np.sqrt(D))-1)     

    return Ps

def SR_IB_2(D):
    """
    # ----------------------------------------------------------------------- #
    #   SR_IB_2: Compute the bootstrapped success rate (function #2)          #
    #                                                                         #
    #            Ps = SR_IB_2(D)                                              # 
    #                                                                         #
    #   This function computes the exact bootstrapped success rate given the  #
    #   diagonal matrix D from LtDL-decomposition of the variance-covariance  #
    #   matrix of float ambiguities                                           #
    #                                                                         #    
    #   Input arguments:                                                      #
    #       D    : Diagonal matrix D from LtDL-decomposition                  #
    #                                                                         #        
    #   Output arguments:                                                     #
    #       Ps   : Exact bootstrapped success rate                            #
    #                                                                         #        
    # ----------------------------------------------------------------------- #
    #   Function : SR_IB_2                                                    #
    #   Date     : 08-FEB-2019                                                #
    #   Author   : Dimitrios Psychas                                          #
    #              Mathematical Geodesy and Positioning,                      #
    #              Delft University of Technology                             # 
    # ----------------------------------------------------------------------- #
    """            

    # Compute the exact bootstrapped success rate using the D matrix
    Ps = np.prod(2*norm.cdf(0.5/np.sqrt(D))-1)     

    return Ps      
    
def ssearch(ahat,L,D,ncands):
    """
    # ----------------------------------------------------------------------- #
    #   ssearch: Search-and-shrink technique for integer-least squares (ILS)  #
    #            ambiguity estimation                                         #
    #                                                                         #
    #            afixed,sqnorm = ssearch(ahat,L,D,ncands)                     # 
    #                                                                         #
    #   This function performs the ILS method based on search-and-shrink      #
    #                                                                         #    
    #   Input arguments:                                                      #
    #       ahat  : Float ambiguities (should be decorrelated for             #
    #               computational efficiency; 1D numpy array)                 #
    #       L,D   : LtDL-decomposition of the variance-covariance matrix of   #
    #               the float ambiguities                                     #
    #       ncands: Number of requested candidates                            #
    #                                                                         #        
    #   Output arguments:                                                     #
    #       afixed: Estimated integers (n-by-ncands numpy array)              #
    #       sqnorm: Corresponding squared norms (ncands-vector, ascendantly   #
    #               ordered, 1D numpy array)                                  #
    #                                                                         #        
    # ----------------------------------------------------------------------- #
    #   Function  : ssearch                                                   #
    #   Date      : 02-SEP-2010                                               #
    #   Author    : Bofeng Li                                                 #
    #               GNSS Research Centre, Department of Spatial Sciences,     #
    #               Curtin University of Technology                           #
    # ----------------------------------------------------------------------- #
    """    
    
    # Is the ambiguity vector an array with entries ?
    if ahat.size>0:
        if not isinstance(ahat,np.ndarray):
            print("The ambiguity vector is not an array!")
            raise SystemExit
    else:
        print("The ambiguity vector is empty!")
        raise SystemExit
        
    # Initializing outputs
    n      = len(ahat)
    afixed = np.zeros((n,ncands))
    sqnorm = np.zeros((1,ncands))

    Chi2      = 1.0e+18 # start with an infinite chi^2
    dist      = np.zeros(n,) # dist(k)=sum_{j=k+1}^{n}(a_j-acond_j)^2/d_j 
    endsearch = 0
    count     = -1 # the number of candidates
    
    # Initialize arrays
    acond = np.zeros(n,)
    zcond = np.zeros(n,)
    step  = np.zeros(n,)

    acond[-1] = ahat[-1]
    zcond[-1] = round(acond[-1])

    left = acond[-1]-zcond[-1]
    step[-1] = np.sign(left)

    if step[-1] == 0:
        step[-1] = 1
    
    imax=ncands-1 # initially, the maximum F(z) is at ncands

    S=np.zeros((n,n)) # used to compute conditional ambiguities

    k=n-1

    # Start the main search-loop
    while not endsearch:
        newdist = dist[k]+left**2/D[0,k]
    
        if newdist < Chi2:
        
            if k != 0:
            
                k-=1
                dist[k]=newdist
                S[k,0:k+1] = S[k+1,0:k+1] + (zcond[k+1]-acond[k+1])*L[k+1,0:k+1]
            
                acond[k] = ahat[k] + S[k,k]
                zcond[k] = round(acond[k])
                left     = acond[k] - zcond[k]
                step[k]  = np.sign(left)
            
                if step[k] == 0:
                    step[k] = 1
                
            else:
            
                # Case 2: store the found candidate and try next valid integer
                if count < ncands-2:
                
                    # Store the first ncands-1 initial points as candidates
                    count += 1
                    afixed[:,count] = zcond[0:n]
                    sqnorm[0,count] = newdist # store F(zcond)     
                
                else:
                
                    afixed[:,imax] = zcond[0:n]
                    sqnorm[0,imax] = newdist               
                
                    imax = np.argmax(sqnorm[0,:])
                    Chi2 = sqnorm[0,imax]
                
                zcond[0] += step[0] # next valid integer
                left      = acond[0] - zcond[0]
                step[0]   = -step[0] - np.sign(step[0])
            
        else:
            
            # Case 3: exit or move up
            if k == n-1:
                endsearch = 1
            else:
                k += 1 # move up
                zcond[k] += step[k] # next valid integer
                left      = acond[k] - zcond[k]
                step[k]   = -step[k] - np.sign(step[k])
    
    # Order the arrays
    order  = np.argsort(sqnorm)     
    sqnorm = sqnorm[0,order]
    afixed = afixed[:,order]
    
    # Return the best integers
    afixed = afixed[:,0]
    sqnorm = sqnorm[0]
    
    return afixed, sqnorm    
    
def parsearch(zhat,Qzhat,Z,L,D,P0=0.995,ncands=2):
    """
    # ----------------------------------------------------------------------- #
    #   parsearch: Partial ambiguity resolution (PAR)                         #
    #                                                                         #
    #              zpar,sqnorm,Qzpar,Zpar,Ps,nfixed,zfixed = \                #
    #                                   parsearch(zhat,Qzhat,Z,L,D,P0,ncands) #
    #                                                                         #
    #   This function performs an integer bootstrapping procedure for partial #
    #   ambiguity resolution with user-defined success-rate P0                #
    #                                                                         #    
    #   Input arguments:                                                      #
    #       zhat  : Decorrelated float ambiguities (1D numpy array)           #
    #       Qzhat : Variance-covariance matrix of decorrelated float          #
    #               ambiguities (2D numpy array)                              #
    #       Z     : Z-transformation matrix from decorrel                     #
    #       L,D   : Lower-triangular and diagonal matrix from LtDL            #
    #               decomposition of Qzhat                                    #
    #       P0    : Minimum required success rate [default=0.995]             #
    #       ncands: Number of requested integer candidate vectors             #
    #               [default=2]                                               #
    #                                                                         #        
    #   Output arguments:                                                     #
    #       zpar  : Subset of fixed ambiguities (nfixed-by-ncands 2D numpy    #
    #               array)                                                    #
    #       sqnorm: Squared norms corresponding to fixed subsets              #
    #       Qzpar : Variance-covariance matrix of float ambiguities for the   #
    #               subset that is fixed                                      #
    #       Zpar  : Z-matrix corresponding to the fixed subset                #
    #       Ps    : Bootstrapped success rate of PAR                          #
    #       nfixed: The number of fixed ambiguities                           #
    #       zfixed: Complete 'fixed' ambiguity vector where the remaining     #
    #               (non-fixed) ambiguities are adjusted according to their   #
    #               correlation with the fixed subset                         #
    #                                                                         #
    #   Note:                                                                 #
    #       This PAR algorithm should be applied to decorrelated ambiguities, #
    #       which can be obtained from the original ahat and its variance-    #
    #       covariance matrix Qahat as:                                       #
    #                                                                         #
    #           Qzhat,Z,L,D,zhat,iZt = decorrel(Qahat,ahat)                   #
    #                                                                         #
    #       The fixed baseline solution is obtained with (see documentation): #
    #                                                                         #
    #           s       = len(zhat)-nfixed                                    #
    #           Qbz     = Qba.dot(Zpar)                                       #
    #           bfixed  = bhat-Qbz.dot(inv(Qzpar)).dot(zhat[s:]-zpar[:,0])    #
    #           Qbfixed = Qbhat-Qbz.dot(inv(Qzpar)).dot(Qbz.T)                #
    #                                                                         #
    #       Hence, zfixed is not required. LAMBDA, however, does give the     #
    #       solution in terms of the full ambiguity vector.                   #
    #                                                                         #        
    # ----------------------------------------------------------------------- #
    #   Function  : parsearch                                                 #
    #   Date      : 04-APR-2012                                               #
    #   Author    : Sandra Verhagen                                           #
    #               GNSS Research Centre, Department of Spatial Sciences,     #
    #               Curtin University of Technology                           #
    #               Mathematical Geodesy and Positioning,                     #
    #               Delft University of Technology                            #
    # ----------------------------------------------------------------------- #
    """   
    
    # Get the number of ambiguities from the dimension of the vc-matrix
    n = len(Qzhat)

    # Compute the bootstrapped success rate if all ambiguities would be fixed
    Ps = SR_IB_2(D)

    k=0
    while Ps<P0 and k<(n-1):
        k+=1
        
        # Compute the bootstrapped success rate if the last n-k+1 ambiguities 
        # would be fixed
        Ps = SR_IB_2(D[0][k:])

    if Ps > P0:
    
        zpar, sqnorm = ssearch(zhat[k:],L[k:,k:],D[[0],k:],2)
    
        Qzpar = Qzhat[k:,k:]
        Zpar  = Z[:,k:]
    
        # First k-1 ambiguities are adjusted based on correlation with the 
        # fixed ambiguities
        QP = Qzhat[:k,k:].dot(inv(Qzhat[k:,k:])) 
    
        zfixed=np.zeros((k,ncands))
        for i in range(0,ncands):
            zfixed[:,i] = zhat[:k] - QP.dot(zhat[k:]-zpar[:,i])
    
        zfixed = np.concatenate((zfixed,zpar),axis=0)
    
        nfixed = n-k
    
    else:
    
        zpar   = []
        Qzpar  = []
        Zpar   = []
        sqnorm = []
        Ps     = np.nan
        zfixed = zhat
        nfixed = 0
        
    return zpar,sqnorm,Qzpar,Zpar,Ps,nfixed,zfixed    
    
def ratioinv(Pf_FIX,Pf_ILS,n):
    """
    # ----------------------------------------------------------------------- #
    #   ratioinv: Get threshold value for Fixed Failure-rate Ratio Test       #
    #                                                                         #
    #             mu = ratioinv(Pf_FIX,Pf_ILS,n)                              # 
    #                                                                         #
    #   This function determines the appropriate threshold value 'mu' for     #
    #   Ratio Test with Fixed Failure rate Pf. Use tabulated values of 'mu'   #
    #   depending on the ILS failure rate and the number of float ambiguities.#
    #                                                                         #    
    #   Input arguments:                                                      #
    #       Pf_FIX: Fixed failure rate (maximum tolerable failure rate)       #
    #               Possible values are 0.010 and 0.001                       #
    #       Pf_ILS: ILS failure rate                                          #
    #       n     : Number of float ambiguities                               #
    #                                                                         #        
    #   Output arguments:                                                     #
    #       mu    : Threshold value for Fixed Failure-rate Ratio Test         #
    #                                                                         #        
    #   Note:                                                                 #
    #       The function loads tables for different values of Pf from text    #
    #       files (table1.txt, table10.txt), which must be in the same        #
    #       directory with this file.                                         #
    # ----------------------------------------------------------------------- #
    #   Function  : ratioinv                                                  #
    #   Date      : 20-APR-2007                                               #
    #   Author    : Sandra Verhagen                                           #
    #               Mathematical Geodesy and Positioning,                     #
    #               Delft University of Technology                            #
    #   Modified  : 20-JAN-2010 by Hans van der Marel                         #
    # ----------------------------------------------------------------------- #
    """      

    # Select the right table for the given fixed failure rate
    kPf = round(Pf_FIX*1000)

    if kPf == 1: # fixed failure rate of 0.1%
        f_in         = open('./table1.txt','r')
        ratio_tab    = np.zeros((31,64))
        ratio_tab[:] = np.nan        
    elif kPf == 10: # fixed failure rate of 1%
        f_in         = open('./table10.txt','r')
        ratio_tab    = np.zeros((31,41))
        ratio_tab[:] = np.nan        
    else:
        print("Incorrect value for Pf_FIX.")
        raise SystemExit
    
    # Read the table
    for line_idx,line in enumerate(f_in,0):
    
        if line:
        
            # Split every read line
            try:
                line = line.split()
            except:
                continue
        
            for i in range(0,len(line)):
                ratio_tab[line_idx,i]=float(line[i])  

    # Make sure n is within the table's range
    if n < 1:
        print('n must be larger than 0')
        raise SystemExit
    if n > len(ratio_tab[0])-1:
        n = len(ratio_tab[0])-2

    # Use linear interpolation to find the treshhold value for the given Pf_ILS
    mu = np.interp(Pf_ILS,ratio_tab[:,0],ratio_tab[:,n])
    
    return mu    
    
def main(ahat,Qahat,method=1,ncands=2,P0=0.995,mu=1):
    """
    # ----------------------------------------------------------------------- #
    #   main: Main LAMBDA function                                            #
    #                                                                         #
    #         afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu = \                          #
    #                                 main(ahat,Qahat,method,ncands,P0,mu)    # 
    #                                                                         #
    #   This is the main function of the LAMBDA software package. By default, #
    #   the ILS method will be used for integer estimation based on the       #
    #   provided float ambiguity vector ahat and associated variance-         #
    #   covariance matrix Qahat. However, the user may also select other      #
    #   methods: integer rounding, bootstrapping or Partial Ambiguity         #
    #   Resolution (PAR). Furthermore, there is the option to apply the Ratio #
    #   Test to decide on acceptance of the fixed solution.                   #
    #                                                                         #
    #   Note 1: LAMBDA always first applies a decorrelation before the        #
    #           integer estimation (for ILS this is required to guarantee an  #
    #           efficient search, for rounding and bootstrapping it is        #
    #           required in order to get higher success rates).               #
    #                                                                         #
    #   Note 2: In this Python version, only one ILS search strategy is       #
    #           implemented (search-and-shrink), instead of the two provided  #
    #           in the Matlab Version 3.0 software. The solution will be the  #
    #           same with both methods, but the search-and-shrink method is   #
    #           known to be faster.                                           #
    #                                                                         #               
    #   Input arguments:                                                      #
    #       ahat  : Float ambiguities (should be decorrelated for             #
    #               computational efficiency; 1D numpy array)                 #
    #       Qahat : Variance-covariance matrix of float ambiguities           #
    #               (2D numpy array)                                          #
    #       method: 1 -> ILS method based on search-and-shrink [default]      #
    #               2 -> Integer rounding method                              #
    #               3 -> Integer bootstrapping method                         #
    #               4 -> PAR with the input P0 of user-defined success rate   #
    #               5 -> ILS method with Ratio Test (uses search-and-shrink)  #
    #       ncands: Number of requested integer candidate vectors             #
    #               (only used with ILS/PAR) [default=2]                      #
    #       P0    : - with method 4 (PAR): minimum required success rate      #
    #                 [default=0.995]                                         #
    #               - with method 5 (ILS + Ratio Test): fixed failure rate    #
    #                 (available options 0.010 or 0.001) [default=0.001]      #
    #       mu    : Fixed threshold value for Ratio Test                      #
    #               (value must be between 0 and 1)                           #
    #                                                                         #        
    #   Output arguments:                                                     #
    #       afixed: Estimated integers (n-by-ncands 2D numpy array), sorted   #
    #               according to the corresponding squared norms,             #
    #               best candidate first.                                     #
    #       sqnorm: Distance between integer candidate and float ambiguity    #
    #               vectors in the metric of the variance-covariance matrix   #
    #               Qahat (ncands vector, 1D numpy array). Only available     #
    #               for ILS.                                                  #
    #       Ps    : Bootstrapped success rate.                                #
    #               If ILS/PAR is used, Ps is its lower bound;                #
    #               If rounding is used, Ps is its upper bound;               #
    #               If bootstrapping is used, Ps is the exact success rate.   #
    #       Qzhat : Variance-covariance matrix of decorrelated ambiguities    #
    #               (corresponding to fixed subset in case of PAR).           #
    #       Z     : Transformation matrix with                                #
    #                   - dimension n-by-n for methods 1-3, 5                 #
    #                   - dimension n-by-nfixed for method 4 (PAR)            #
    #       nfixed: Number of fixed ambiguities                               #
    #                   - with methods 1 to 3: will always be equal to n      #
    #                   - with method 4 (PAR): will be equal to the number of #
    #                     fixed decorrelated ambiguities                      #
    #                   - with method 5 (ILS + Ratio Test): will be equal to  #
    #                     n if fixed solution is accepted, and 0 otherwise    #
    #       mu    : Threshold value used for Ratio Test                       #
    #                                                                         #
    #   Note 3: The Ratio Test will only be applied if method 5 is selected.  #
    #           By default, the Fixed Failure-rate Ratio Test will be applied.#
    #           The thresold value 'mu' is then determined such that the      #
    #           failure rate will not exceed the specified P0 (default=0.001).#
    #           If 'mu' is specified, the value for P0 is ignored.            #
    #                                                                         #
    #   Note 4: The Ratio Test used here is:                                  #
    #                                                                         #
    #           Accept afixed iff: sqnorm[0]/sqnorm[1] <= mu                  #
    #                                                                         #
    #           Hence, the squared norm of the best (ILS) integer solution is #
    #           in the numerator. In literature often the reciprocal is used; #
    #           the corresponding critical value is then c=1/mu.              #
    #                                                                         #        
    # ----------------------------------------------------------------------- #
    #   Function  : main                                                      #
    #   Date      : 01-SEP-2012                                               #
    #   Author    : Bofeng Li and Sandra Verhagen                             #
    #               GNSS Research Centre, Department of Spatial Sciences,     #
    #               Curtin University of Technology                           #
    #               Mathematical Geodesy and Positioning,                     #
    #               Delft University of Technology                            #
    # ----------------------------------------------------------------------- #
    """       
    
    # If ILS/PAR method is not used, only a unique solution is available
    if method in [2,3]:
        ncands = 1
        
    # Default values for Ratio Test
    if method == 5:
        FFRT = 1
        if P0 == 0.995:
            P0 = 0.001        
    else:
        mu = 1
        
    # Get the number of ambiguities from the dimension of the vc-matrix
    n      = len(Qahat)
    nfixed = n
    sqnorm = []
    
    # Tests on inputs ahat and Qahat
    
    # Is the Qahat symmetric ?
    if ((Qahat-Qahat.T)<1e-8).all()==False:
        print("Variance-covariance matrix of float estimated ambiguities " + \
              "is not symmetric!")
        raise SystemExit
        
    # Is the Qahat positive-definite ?    
    
    # Get the eigenvalues (D) and eigenvectors (V) of Qahat
    D,V = eig(Qahat)
    
    # Check the positive-definiteness    
    if sum(D>0) != n:
        print("Variance-covariance matrix of float estimated ambiguities " + \
              "is not positive definite!")
        raise SystemExit
        
    # Do the Qahat and ahat have identical dimensions ?
    if len(ahat) != n:
        print("Variance-covariance matrix and vector of float estimated " + \
              "ambiguities do not have identical dimensions")
        raise SystemExit
        
    # Is the ambiguity vector an array ?
    if ahat.size > 0:
        if not isinstance(ahat,np.ndarray):
            print("The ambiguity vector is not an array!")
            raise SystemExit
    else:
        print("The ambiguity vector is empty!")
        raise SystemExit
        
    # Remove integer numbers from float solution, so that all values are 
    # between -1 and 1 (for computational convenience only)
    ahat,incr = np.modf(ahat)
    
    # Compute the Z matrix based on the decomposition Q=L^T*D*L. The transformed
    # float solution: \hat{a}=Z^T*ahat, Qzhat=Z^T*Qahat*Z
    Qzhat,Z,L,D,zhat,iZt = decorrel(Qahat,ahat)
    
    # Compute the bootstrapped success rate
    Ps = SR_IB_2(D)
    
    # Integer ambiguity estimation (and validation)    
    if method == 1: # ILS with shrinking search
        zfixed,sqnorm = ssearch(zhat,L,D,ncands)
    elif method == 2: # Integer rounding    
        zfixed = np.round(zhat)
    elif method == 3: # Integer bootstrapping
        zfixed = bootstrap(zhat,L)
    elif method == 4: # PAR
        zpar,sqnorm,Qzhat,Z,Ps,nfixed,zfixed = \
                                          parsearch(zhat,Qzhat,Z,L,D,P0,ncands)
        if nfixed == 0:
            ncands = 1
    elif method == 5: # ILS with Ratio Test
        zfixed,sqnorm = ssearch(zhat,L,D,ncands)
        if FFRT == 1:
            if 1-Ps>P0:
                mu = ratioinv(P0,1-Ps,n)
            else:
                mu = 1
                
        # Perform Ratio Test
        if sqnorm[0]/sqnorm[1] > mu: # rejection: keep float solution
            zfixed = zhat
            nfixed = 0
            ncands = 1
    
    # Perform the back-transformation and add the increments
    afixed = iZt.dot(zfixed)    
    if ncands > 1:
        afixed += numpy.matlib.repmat(incr.reshape(n,1),1,ncands)
    else:
        afixed += incr

    return afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu                                                                 