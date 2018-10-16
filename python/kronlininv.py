#!/usr/bin/env python


import numpy as NP
import scipy.linalg as LA 
from kronlininv_f2py import wrapforpy as KLwr
import sys


##########################################################
##########################################################

def calcfactors( G1,G2,G3,Cm1,Cm2,Cm3,Cd1,Cd2,Cd3 ) :
    """
       Compute the factors involved to solve the inverse problem
       These factors are the ones to be stored to compute 
          the posterior covariance or mean
    """

    print 'Preliminary checks'
    lsm = [Cm1,Cm2,Cm3,Cd1,Cd2,Cd3]
    lsn = ["Cm1","Cm2","Cm3","Cd1","Cd2","Cd3"]
    for i,mat in enumerate(lsm) :
        checkmat( mat,lsn[i] )

    print 'Calculating preliminary stuff'
    iCdG1 = LA.solve(Cd1,G1, sym_pos=True)
    iCdG2 = LA.solve(Cd2,G2, sym_pos=True)
    iCdG3 = LA.solve(Cd3,G3, sym_pos=True)
    GtiCd1  = NP.dot( G1.T, iCdG1 )
    GtiCd2  = NP.dot( G2.T, iCdG2 ) 
    GtiCd3  = NP.dot( G3.T, iCdG3 )

    ###################################
    ###  Eigendecomposition Python    #
    ###################################
    print 'Calculating  eigendecomposition'
    lambda1,U1 = LA.eigh(GtiCd1,b=Cm1, lower=False, eigvals_only=False,
                         overwrite_a=False, overwrite_b=False, turbo=True,
                         eigvals=None, type=3, check_finite=True)
    lambda2,U2 = LA.eigh(GtiCd2,b=Cm2, lower=False, eigvals_only=False,
                         overwrite_a=False, overwrite_b=False, turbo=True,
                         eigvals=None, type=3, check_finite=True)
    lambda3,U3 = LA.eigh(GtiCd3,b=Cm3, lower=False, eigvals_only=False,
                         overwrite_a=False, overwrite_b=False, turbo=True,
                         eigvals=None, type=3, check_finite=True)

    ###################################
    ### calculating the 3 factors

    print 'Calculating fa'
    n1=lambda1.size
    n2=lambda2.size
    n3=lambda3.size

    print 'Calculating fb'
    argfb = 1.00 + krondiag(lambda1,krondiag(lambda2,lambda3))
    vecdfac = 1.00/argfb
    
    print 'Calculating fc'
    iUCm1 = LA.solve(U1,Cm1) 
    iUCm2 = LA.solve(U2,Cm2) 
    iUCm3 = LA.solve(U3,Cm3) 

    print 'Calculating fd'
    fd11 = LA.solve( Cd1, NP.identity(Cd1.shape[0])) 
    fd22 = LA.solve( Cd2, NP.identity(Cd2.shape[0])) 
    fd33 = LA.solve( Cd3, NP.identity(Cd3.shape[0])) 
    iUCmGtiCd1 = NP.dot( iUCm1, NP.dot(G1.T, fd11 )) 
    iUCmGtiCd2 = NP.dot( iUCm2, NP.dot(G2.T, fd22 )) 
    iUCmGtiCd3 = NP.dot( iUCm3, NP.dot(G3.T, fd33 ))
          
    return U1,U2,U3, vecdfac, iUCm1,iUCm2,iUCm3, iUCmGtiCd1,iUCmGtiCd2,iUCmGtiCd3

##########################################################

def postmean(U1,U2,U3, diaginvlambda, Z1,Z2,Z3, G1,G2,G3, mprior, dobs ) :
    """
       Posterior mean

    """      
    # sizes
    Na = mprior.size
    Nb = dobs.size
    Ni = Z1.shape[0]
    Nl = Z1.shape[1]
    Nj = Z2.shape[0]
    Nm = Z2.shape[1]
    Nk = Z3.shape[0]
    Nn = Z3.shape[1]

    ## allocate stuff
    ddiff = NP.zeros(Nb)
    Zh = NP.zeros(Na)
    elUDZh = NP.zeros(Na)

    av = NP.arange(0,Na) 
    bv = NP.arange(0,Nb) 

    iv = av/(Nk*Nj) 
    jv = (av-iv*Nk*Nj)/Nk 
    kv = av-jv*Nk-iv*Nk*Nj 
    
    lv = bv/(Nn*Nm) 
    mv = (bv-lv*Nn*Nm)/Nn 
    nv = bv-mv*Nn-lv*Nn*Nm
    
    ##---------------------
    if Nb>20 :
        everynit = Nb/20
    else :
        everynit = 1

    postm = NP.zeros(Na)
    
    for b in range(Nb) :
        ddiff[b] = dobs[b] - sum( mprior * G1[lv[b],iv] * G2[mv[b],jv] * G3[nv[b],kv] )

    for a in range(Na) :
        if a%everynit==0 :
            line = 'postmean(): 1st loop {} of {}   \r'.format(a,Na)
            sys.stdout.write(line)
            sys.stdout.flush()

        Zh[a] = sum( Z1[iv[a],lv] * Z2[jv[a],mv] * Z3[kv[a],nv] * ddiff )


    for a in range(Na) :
        if a%everynit==0 :
            line = 'postmean(): 2nd loop {} of {}   \r'.format(a,Na)
            sys.stdout.write(line)
            sys.stdout.flush()
       
        elUDZh[a] = NP.sum( U1[iv[a],iv] * U2[jv[a],jv] * U3[kv[a],kv]  * diaginvlambda * Zh )
      
        # element of the posterior mean
        postm[a] = mprior[a] + elUDZh[a]       

    return postm

##########################################################

def blockpostcov(U1,U2,U3, diaginvlambda, iUCm1,iUCm2,iUCm3, astart,aend,bstart,bend ) :
    
    s1 = aend+1-astart
    s2 = bend+1-bstart
    
    postC = NP.zeros((s1,s2))

    # sizes
    Ni = U1.shape[0]
    Nl = U1.shape[1]
    Nj = U2.shape[0]
    Nm = U2.shape[1]
    Nk = U3.shape[0]
    Nn = U3.shape[1]
    Na = Ni*Nj*Nk
    Nb = Nl*Nm*Nn

    ## -----
    av = NP.arange(0,Na) 
    bv = NP.arange(0,Nb) 

    iv = av/(Nk*Nj) 
    jv = (av-iv*Nk*Nj)/Nk 
    kv = av-jv*Nk-iv*Nk*Nj 
    
    lv = bv/(Nn*Nm) 
    mv = (bv-lv*Nn*Nm)/Nn 
    nv = bv-mv*Nn-lv*Nn*Nm

    ## -----
    if s1>20 :
        everynit = s1/20
    else :
        everynit = 1


    for a in range(astart,aend+1)  :
        
        if a%everynit==0 :
            line = 'blockpostcov(): loop {} of {}   \r'.format(a,aend)
            sys.stdout.write(line)
            sys.stdout.flush()

        # row first two factors
        # a row x diag matrix 
        row2 =  U1[iv[a],iv] * U2[jv[a],jv] * U3[kv[a],kv] * diaginvlambda
        
        for b in range(bstart,bend+1) :
            
            ## calculate one row of first TWO factors
            col1 = iUCm1[iv,iv[b]] * iUCm2[jv,jv[b]] * iUCm3[kv,kv[b]]
            
            ## calculate one element of the posterior covariance
            postC[a,b] = NP.sum(row2*col1)

    return postC

##########################################################
  
def krondiag(a,b) :
    """
      Kronecker product of two diagonal matrices
      Returns only the diagonal as a vector
    """
    assert a.ndim==1
    assert b.ndim==1
    ni=a.size
    nj=b.size
    c = NP.zeros(ni*nj,dtype=a.dtype) ###,dtype=NP.complex128)
    k=0
    for i in range(ni) :
        c[k:k+nj] = a[i]*b
        k += nj
    return c

##########################################################

def fortran_postmean(Us,diaginvlambda,Zs,Gs,mprior,dobs) :

    ## Posterior mean
    print 'Calculating posterior mean'
    postm = NP.zeros(mprior.size,dtype=NP.float64)
    postm = KLwr.c_posteriormean(U1,U2,U3,invlambda,iUCmGtiCd1,iUCmGtiCd2,
                                 iUCmGtiCd3,G1,G2,G3,mprior,dobs)
    return postm
                                           
##########################################################

def fortran_blockpostcov( U1,U2,U3,invlambda,iUCm1,iUCm2,iUCm3,astart,aend,bstart,bend ) :
    
    ## Fortran indices start from 1 !!
    postcov = NP.zeros((aend-astart+1,bend-bstart+1),dtype=NP.float64)
    print 'tens_covm_post',tens_covm_post.shape
    postcov = KLwr.c_blockpostcov(U1,U2,U3,invlambda,iUCm1,iUCm2,iUCm3,astart,aend,bstart,bend)

    return postcov

##########################################################
        
def fortran_calcfactors( G1,G2,G3,Cm1,Cm2,Cm3,Cd1,Cd2,Cd3) :

    print 'Preliminary checks'
    lsm = [Cm1,Cm2,Cm3,Cd1,Cd2,Cd3]
    lsn = ["Cm1","Cm2","Cm3","Cd1","Cd2","Cd3"]
    for i,mat in enumerate(lsm) :
        checkmat( mat,lsn[i] )
        
    ## calculate factors
    U1,U2,U3,invlambda,iUCm1,iUCm2,iUCm3,Z1,Z2,Z3 = KLwr.c_calcfactors(G1,G2,G3,Cm1,Cm2,Cm3,Cd1,Cd2,Cd3)
    return U1,U2,U3,invlambda,iUCm1,iUCm2,iUCm3,Z1,Z2,Z3

##########################################################
        
def checkmat( mat,matname ):
    """
      Check some properties of input matrix
    """
    print 'Checking {}'.format(matname)
    issym = NP.allclose(mat,mat.T)
    if not issym :
        print 'Matrix {} is not symmetric'.format(matname)
        exit() 
    try:
        NP.linalg.cholesky(mat)
    except NP.linalg.linalg.LinAlgError:
        print 'Matrix {} is not positive definite'.format(matname)
        exit()
    return True
    
##########################################################


