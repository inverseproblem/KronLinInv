
#------------------------------------------------------------------------
#
#    Copyright 2018, Andrea Zunino 
#
#    This file is part of KronLinInv.
#
#    KronLinInv is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    KronLinInv is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with KronLinInv.  If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------


##==========================================================
module kronlininv

using Distributed
using Libdl
using LinearAlgebra

export calcfactors,posteriormean,blockpostcov
export fortran_calcfactors,fortran_posteriormean,fortran_blockpostcov

##fortranlibfl = "../src/julia/kron_utils.so"
## fortranlibfl = string(@__DIR__,"./kron_utils.so")

##==========================================================

#-------------------------------------------------------------------------
function fortran_calcfactors(G1::Array{Float64,2},G2::Array{Float64,2},G3::Array{Float64,2},
                             Cm1::Array{Float64,2},Cm2::Array{Float64,2},Cm3::Array{Float64,2},
                             Cd1::Array{Float64,2},Cd2::Array{Float64,2},Cd3::Array{Float64,2},
                             fortranlibfl::String)
    
    forlib = Libdl.dlopen(fortranlibfl)
    wcfac = Libdl.dlsym(forlib,:c_calcfactors)

    #---
    n1,m1 = size(G1)
    n2,m2 = size(G2)
    n3,m3 = size(G3)
    ### Ref{T}  Behaves like a Ptr{T} that owns its memory.
    U1 = Array{Float64,2}(undef,m1,m1)
    U2 = Array{Float64,2}(undef,m2,m2)
    U3 = Array{Float64,2}(undef,m3,m3)
    diaginvlambda = Array{Float64,1}(undef,m1*m2*m3)
    iUCm1 = Array{Float64,2}(undef,m1,m1)
    iUCm2 = Array{Float64,2}(undef,m2,m2)
    iUCm3 = Array{Float64,2}(undef,m3,m3)
    iUCmGtiCd1 = Array{Float64,2}(undef,m1,n1)
    iUCmGtiCd2 = Array{Float64,2}(undef,m2,n2)
    iUCmGtiCd3 = Array{Float64,2}(undef,m3,n3)

    # In Julia code wrapping calls to external Fortran routines, all
    # input arguments should be declared as of type Ref{T}, as Fortran
    # passes all variables by reference. The return type should either be
    # Void for Fortran subroutines, or a T for Fortran functions returning the type T.
    ccall( wcfac, Cvoid,
          (Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint},
           Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},
           Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},
           Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},
           Ref{Cdouble},Ref{Cdouble},
           Ref{Cdouble},Ref{Cdouble},Ref{Cdouble}),
           n1,n2,n3,m1,m2,m3,
           G1,G2,G3,Cm1,Cm2,Cm3,Cd1,Cd2,Cd3,
           U1,U2,U3,diaginvlambda,iUCm1,iUCm2,iUCm3,
           iUCmGtiCd1,iUCmGtiCd2,iUCmGtiCd3)
    
    return U1,U2,U3,diaginvlambda,iUCm1,iUCm2,iUCm3,iUCmGtiCd1,iUCmGtiCd2,iUCmGtiCd3
end
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
function fortran_posteriormean(U1::Array{Float64,2},U2::Array{Float64,2},U3::Array{Float64,2},
                               diaginvlambda::Array{Float64,1},iUCmGtiCd1::Array{Float64,2},
                               iUCmGtiCd2::Array{Float64,2},iUCmGtiCd3::Array{Float64,2},
                               G1::Array{Float64,2},G2::Array{Float64,2},G3::Array{Float64,2},
                               mprior::Array{Float64,1},dobs::Array{Float64,1}, fortranlibfl::String)
    
    forlib = Libdl.dlopen(fortranlibfl)
    wpmean = Libdl.dlsym(forlib,:c_posteriormean)

    #------------------
    n1,m1 = size(G1)
    n2,m2 = size(G2)
    n3,m3 = size(G3)
    postm = Array{Float64,1}(undef,m1*m2*m3)
    # In Julia code wrapping calls to external Fortran routines, all
    # input arguments should be declared as of type Ref{T}, as Fortran
    # passes all variables by reference. The return type should either be
    # Void for Fortran subroutines, or a T for Fortran functions returning the type T.
    ccall( wpmean, Cvoid,
          (Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint},
           Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},
           Ref{Cdouble},Ref{Cdouble},
           Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},
           Ref{Cdouble}),
           n1,n2,n3,m1,m2,m3,
           U1,U2,U3,diaginvlambda,iUCmGtiCd1,iUCmGtiCd2,iUCmGtiCd3,
           G1,G2,G3, mprior,dobs, postm)
    
    return postm
end
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
function fortran_blockpostcov(U1::Array{Float64,2},U2::Array{Float64,2},U3::Array{Float64,2},
                              diaginvlambda::Array{Float64,1},iUCm1::Array{Float64,2},
                              iUCm2::Array{Float64,2},iUCm3::Array{Float64,2},
                              astart::Int64,aend::Int64,bstart::Int64,bend::Int64,
                              fortranlibfl::String)
    
    forlib = Libdl.dlopen(fortranlibfl)
    wblcov = Libdl.dlsym(forlib,:c_blockpostcov)

    #------------------
    m1 = size(U1,1)
    m2 = size(U2,1)
    m3 = size(U3,1)
    postC = Array{Float64,2}(undef,aend-astart+1,bend-bstart+1)
    
    # In Julia code wrapping calls to external Fortran routines, all
    # input arguments should be declared as of type Ref{T}, as Fortran
    # passes all variables by reference. The return type should either be
    # Void for Fortran subroutines, or a T for Fortran functions returning the type T.
    ccall( wblcov, Cvoid,
          (Ref{Cint},Ref{Cint},Ref{Cint},
           Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},
           Ref{Cdouble},Ref{Cdouble},
           Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cdouble}),
           m1,m2,m3,
           U1,U2,U3,diaginvlambda,iUCm1,iUCm2,iUCm3,
           astart,aend,bstart,bend, postC)
    
    return postC
end
#-------------------------------------------------------------------------


##==========================================================

function krondiag(a::Array{Float64,1},b::Array{Float64,1}) 
    """
        Kronecker product of two diagonal matrices
        Returns only the diagonal as a vector
    """
    ni=size(a,1)
    nj=size(b,1)
    c=zeros(Float64,ni*nj)
    k=1
    for i=1:ni 
        c[k:k+nj-1] = a[i]*b
        k += nj
    end
    return c
end

##==========================================================

function calcfactors(G1::Array{Float64,2},G2::Array{Float64,2},G3::Array{Float64,2},
                     Cm1::Array{Float64,2},Cm2::Array{Float64,2},Cm3::Array{Float64,2},
                     Cd1::Array{Float64,2},Cd2::Array{Float64,2},Cd3::Array{Float64,2} )
    
    ##----------------
    lsn = ["Cm1","Cm2","Cm3","Cd1","Cd2","Cd3"]
    for (i,C) in enumerate( [Cm1,Cm2,Cm3,Cd1,Cd2,Cd3] )
        #assert( isposdef( C ) )
        if isposdef( C )==false
            println("\n $(lsn[i]) is not positive definite. Aborting \n")
            exit()
        end
    end
    
    ##----------------
    iCdG1 = Cd1 \ G1
    iCdG2 = Cd2 \ G2
    iCdG3 = Cd3 \ G3
    GtiCd1 = transpose(G1) * iCdG1
    GtiCd2 = transpose(G2) * iCdG2
    GtiCd3 = transpose(G3) * iCdG3
    GtiCdG1 = transpose(G1) * iCdG1 
    GtiCdG2 = transpose(G2) * iCdG2 
    GtiCdG3 = transpose(G3) * iCdG3 

    
    ## If itype = 3, the problem to solve is B * A * x = lambda * x.
    itype = 3 ## the problem to solve is B * A * x = lambda * x.
    uplo = 'L'
    jobz = 'V'
    
    ## Finds the generalized eigenvalues (jobz = N) or eigenvalues and eigenvectors (jobz = V)
    ## of a symmetric matrix A and symmetric positive-definite matrix B. If uplo = U, the upper
    ## triangles of A and B are used. If uplo = L, the lower triangles of A and B are used. If
    ## itype = 1, the problem to solve is A * x = lambda * B * x. If itype = 2, the problem to
    ## solve is A * B * x = lambda * x. If itype = 3, the problem to solve is
    ## B * A * x = lambda * x.
    ## sygvd!(itype::Integer,jobz::Char,uplo::Char,
    ##        A::StridedMatrix{$elty},B::StridedMatrix{$elty})
    ## print 'Calculating fa' etc

    lambda1,U1,LcholB1 = LAPACK.sygvd!(itype,jobz,uplo,GtiCdG1,copy(Cm1))
    lambda2,U2,LcholB2 = LAPACK.sygvd!(itype,jobz,uplo,GtiCdG2,copy(Cm2))
    lambda3,U3,LcholB3 = LAPACK.sygvd!(itype,jobz,uplo,GtiCdG3,copy(Cm3))

    
    ## print 'Calculating fb'
    argfb = 1.0 .+ krondiag(lambda1,krondiag(lambda2,lambda3))
    vecdfac = 1.0./argfb
    
    ## print 'Calculating fc'
    iUCm1 = U1 \ Cm1 
    iUCm2 = U2 \ Cm2
    iUCm3 = U3 \ Cm3

    ## print 'Calculating fd'
    fd11 = Cd1 \ Matrix{Float64}(I,size(Cd1,1),size(Cd1,2)) #eye(Cd1)
    fd22 = Cd2 \ Matrix{Float64}(I,size(Cd2,1),size(Cd2,2)) ##eye(Cd2) 
    fd33 = Cd3 \ Matrix{Float64}(I,size(Cd3,1),size(Cd3,2)) ##eye(Cd3)

    iUCmGtiCd1 = iUCm1 * transpose(G1) * fd11 
    iUCmGtiCd2 = iUCm2 * transpose(G2) * fd22
    iUCmGtiCd3 = iUCm3 * transpose(G3) * fd33

    return  U1,U2,U3, vecdfac, iUCm1,iUCm2,iUCm3, iUCmGtiCd1,iUCmGtiCd2,iUCmGtiCd3
end

##======================================================================

function posteriormean(U1::Array{Float64,2},U2::Array{Float64,2},U3::Array{Float64,2},
                       diaginvlambda::Array{Float64,1},
                       Z1::Array{Float64,2},Z2::Array{Float64,2},Z3::Array{Float64,2},
                       G1::Array{Float64,2},G2::Array{Float64,2},G3::Array{Float64,2},
                       mprior::Array{Float64,1}, dobs::Array{Float64,1};
                       runparallel=false)

    ## sizes
    Ni = size(Z1,1)
    Nl = size(Z1,2)
    Nj = size(Z2,1)
    Nm = size(Z2,2)
    Nk = size(Z3,1)
    Nn = size(Z3,2)

    ## sizes
    # Nr12 = size(U1,2)*size(U2,2)*size(U3,2)
    # Nc1  = size(Z1,1)*size(Z2,1)*size(Z3,1)
    Na   = size(mprior,1)
    Nb   = size(dobs,1)

    ##-------------
    av = 1:Na
    bv = 1:Nb
    ## vectors containing all possible indices for 
    ##    row calculations of Kron prod AxBxC
    iv = div.( (av.-1), (Nk.*Nj) ) .+1 
    jv = div.( (av.-1 .-(iv.-1).*Nk.*Nj), Nk ) .+1 
    kv = av.-(jv.-1).*Nk.-(iv.-1).*Nk.*Nj 
    ## vectors containing all possible indices for
    ##    column calculations of Kron prod AxBxC
    lv =  div.( (bv .-1), (Nn.*Nm) ) .+ 1 
    mv =  div.( (bv.-1 .-(lv.-1).*Nn.*Nm), Nn ) .+ 1 
    nv =  bv.-(mv.-1).*Nn-(lv.-1).*Nn.*Nm 
    ##  Gs have different shape than Us !!
    
    ####===================================

    if runparallel==true 
        ####================================================
        ##     Parallel version
        ####================================================
        startt = time()

        ## get the ids of cores
        idcpus = workers()
        firstwork = idcpus[1]
        numwork = nworkers()
        println("posteriormean(): parallel run using $numwork cpus")
        
        ##################################################
        #             loop 1                             #
        ##################################################
        ## Nb
        scheduling,looping = spreadwork(Nb,numwork,1) ## Nb !!
        everynit = looping[1,2]>100 ? div(looping[1,2],100) : 2

        ddiff = Array{Float64}(undef,Nb)
        @sync begin
            for ip=1:numwork 
                bstart,bend = looping[ip,1],looping[ip,2]

                @async ddiff[bstart:bend] = remotecall_fetch(comp_ddiff,idcpus[ip],
                                                             everynit,firstwork,
                                                             iv,lv,jv,mv,kv,nv,
                                                             G1,G2,G3,mprior,dobs,
                                                             bstart,bend)
            end
        end 

        ##################################################
        #             loop 2                             #
        ##################################################
        ### need to re-loop because full Zh is needed
        startt = time()
        ## Na
        scheduling,looping = spreadwork(Na,numwork,1) ## Na!!
        
        Zh = Array{Float64}(undef,Na)
        @sync begin
            for ip=1:numwork 
                astart,aend = looping[ip,1],looping[ip,2]
                ## distribute work to specific cores (spawnAT)
                @async Zh[astart:aend] = remotecall_fetch(comp_Zh,idcpus[ip],
                                                          everynit,firstwork,
                                                          iv,lv,jv,mv,kv,nv,
                                                          Z1,Z2,Z3,ddiff,
                                                          astart,aend)
            end
        end
        
        ##################################################
        #             loop 3                             #
        ##################################################
        startt = time()
        ## Na
        scheduling,looping = spreadwork(Na,numwork,1) ## Na!!

        postm = Array{Float64}(undef,Na)
        @sync begin
            for ip=1:numwork 
                astart,aend = looping[ip,1],looping[ip,2]
                ## distribute work to specific cores 
                @async postm[astart:aend] = remotecall_fetch(comp_postm,idcpus[ip],
                                                            everynit,firstwork,
                                                            iv,lv,jv,mv,kv,nv,
                                                            U1,U2,U3,
                                                            diaginvlambda,Zh,
                                                            mprior,astart,aend)
            end
        end
        
        ##############################################################


    else 
        ####================================================
        ##     Serial version
        ####================================================
        postm  = Array{Float64}(undef,Na)
        ddiff  = Array{Float64}(undef,Nb)    
        Zh     = Array{Float64}(undef,Na)
        elUDZh = Array{Float64}(undef,Na)

        everynit = Na>20 ? div(Na,20) : 1
        
        startt = time()
        ## dobs - dcalc(mprior)
        for b=1:Nb
            if (b%everynit==0) | (b==2)
                eta =  ( (time()-startt)/float(b-1) * (Na-b+1) ) /60.0
                reta = round(eta,digits=3)
                print("loop 1/3: $b of $Nb; ETA: $reta min \r")
                flush(stdout)
            end

            # !!---------------------------------------------------------------------------
            # ddiff(b) = dobs(b) - sum(mprior * G1(lv(b),iv) * G2(mv(b),jv) * G3(nv(b),kv))
            # !!---------------------------------------------------------------------------

            datp = 0.0
            for j=1:Na
                elG = G1[lv[b],iv[j]] * G2[mv[b],jv[j]] * G3[nv[b],kv[j]]
                datp = datp +  mprior[j] * elG
            end        
            ddiff[b] = dobs[b] - datp
            #println("$b $(ddiff[b]) $(dobs[b])")
        end

        startt = time()
        
        for i=1:Na

            if (i%everynit==0) | (i==2)
                eta =  ( (time()-startt)/float(i-1) * (Na-i+1) ) /60.0
                reta = round(eta,digits=3)
                print("loop 2/3: $i of $Na; ETA: $reta min \r")
                flush(stdout)
            end                     
            
            ## compute Zh
            Zh[i]=0.0
            for j=1:Nb
                tZZ = Z1[iv[i],lv[j]] * Z2[jv[i],mv[j]] * Z3[kv[i],nv[j]]
                Zh[i] = Zh[i] + tZZ * ddiff[j]
            end
        end
        
        startt = time()
        ### need to re-loop because full Zh is needed
        for i=1:Na

            if (i%everynit==0) | (i==2)
                eta =  ( (time()-startt)/float(i+1) * (Na-i+1) ) /60.0
                reta = round(eta,digits=3)
                print("loop 3/3: $i of $Na; ETA: $reta min \r")
                flush(stdout)
            end
            
            ## UD times Zh
            elUDZh[i] = 0.0
            for j=1:Na
                # element of row of UD
                elrowUD = U1[iv[i],iv[j]] * U2[jv[i],jv[j]] *
                    U3[kv[i],kv[j]] * diaginvlambda[j]
                # element of final vector
                elUDZh[i] = elUDZh[i] + elrowUD * Zh[j]
            end

            ## element of the posterior mean
            postm[i] = mprior[i] + elUDZh[i] # sum(bigmatrow.*ddiff)
        end
        print("\n")
    end
    
    return postm
end

##============================================================================

function blockpostcov(U1::Array{Float64,2},U2::Array{Float64,2},U3::Array{Float64,2},
                      diaginvlambda::Array{Float64,1},
                      iUCm1::Array{Float64,2},iUCm2::Array{Float64,2},iUCm3::Array{Float64,2},
                      astart::Int64,aend::Int64,bstart::Int64,bend::Int64 )

    Ni = size(U1,1)
    Nl = size(U1,2)
    Nj = size(U2,1)
    Nm = size(U2,2)
    Nk = size(U3,1)
    Nn = size(U3,2)
    Na = Ni*Nj*Nk
    Nb = Nl*Nm*Nn
    
    ## sizes----
    # Nr12 = size(U1,2)*size(U2,2)*size(U3,2)
    # Nc1  = size(Z1,1)*size(Z2,1)*size(Z3,1)

    ##-----------------
    av = 1:Na
    ## vectors containing all possible indices for 
    ##    row calculations of Kron prod AxBxC
    iv =  div( (av-1), (Nk*Nj) ) +1 
    jv =  div( (av-1-(iv-1)*Nk*Nj), Nk ) + 1 
    kv =  av-(jv-1)*Nk-(iv-1)*Nk*Nj 

    nci = aend-astart+1
    ncj = bend-bstart+1
    postC = Array{Float64}(undef,nci,ncj)
    row2  = Array{Float64}(undef,Na)
    col1  = Array{Float64}(undef,Nb)

    #tic()
    for a=astart:aend
        if a%100==0
            println("a: $a of $(astart):$(aend)")
        end
        ## calculate one row of first two factors
        ## row of  Kron prod AxBxC times a diag matrix (fb)
        for q=1:Na
            row2[q] = U1[iv[a],iv[q]] * U2[jv[a],jv[q]] * U3[kv[a],kv[q]] * diaginvlambda[q]
        end
        for b=bstart:bend
            ## : vectorize
            ## : turns off subscript checking  
            postC[a,b] = 0.0
            for p=1:Na
                ## calculate one column of fc
                col1[p] = iUCm1[iv[p],iv[b]] * iUCm2[jv[p],jv[b]] * iUCm3[kv[p],kv[b]]
                ## calculate one element 
                postC[a,b] = postC[a,b] + row2[p] * col1[p]  
            end
        end
    end
    #toc()
    return postC
end

##==========================================================
##==========================================================

##==========================================================

function comp_ddiff(everynit::Int64,firstwork::Int64,
                    iv::Array{Int64,1},lv::Array{Int64,1},jv::Array{Int64,1},
                    mv::Array{Int64,1},kv::Array{Int64,1},nv::Array{Int64,1},
                    G1::Array{Float64,2},G2::Array{Float64,2},G3::Array{Float64,2},
                    mprior::Array{Float64,1},dobs::Array{Float64,1},
                    bstart::Int64,bend::Int64)
    ## dobs - dcalc(mprior)
    @assert bend>=bstart
    startt = time()
    myNb = bend-bstart+1
    ddiff = zeros(Float64,myNb)
    Na = length(mprior)
    # loop on chunk
    for b=bstart:bend # b=1:Nb
        myb = b-bstart+1

        # print info
        if (myb%everynit==0) | (myb==5)
            if myid()==firstwork
                eta = ( (time()-startt)/float(myb-1) * (myNb-myb+1) ) /60.0
                reta = round(eta,digits=3)
                println("[parallel] loop 1/3 b: $myb of $myNb; ETA: $reta min ") #  \r")
                #flush(stdout)
            end
        end

        datp = 0.0
        for j=1:Na
            elG = G1[lv[b],iv[j]] * G2[mv[b],jv[j]] * G3[nv[b],kv[j]]
            datp = datp +  mprior[j] * elG
        end        
        ddiff[myb] = dobs[b] - datp
    end
    return ddiff
end

##==========================================================

function comp_Zh(everynit::Int64,firstwork::Int64,
                 iv::Array{Int64,1},lv::Array{Int64,1},jv::Array{Int64,1},
                 mv::Array{Int64,1},kv::Array{Int64,1},nv::Array{Int64,1},
                 Z1::Array{Float64,2},Z2::Array{Float64,2},Z3::Array{Float64,2},
                 ddiff::Array{Float64,1},
                 astart::Int64,aend::Int64)
    ## compute Zh
    @assert aend>=astart
    startt = time()
    myNa = aend-astart+1
    Nb = length(ddiff)
    Zh = zeros(Float64,myNa)
    # loop on chunk
    for i=astart:aend #i=1:Na
        mya = i-astart+1

        # print info
        if (mya%everynit==0) | (mya==5)
            if myid()==firstwork            
                eta = ( (time()-startt)/float(mya-1) * (myNa-mya+1) ) /60.0
                reta = round(eta,digits=3)
                println("[parallel] loop 2/3 a: $mya of $myNa; ETA: $reta min ") #  \r")
                #flush(stdout)
            end
        end

        Zh[mya]=0.0
        for j=1:Nb
            tZZ = Z1[iv[i],lv[j]] * Z2[jv[i],mv[j]] * Z3[kv[i],nv[j]]
            Zh[mya] = Zh[mya] + tZZ * ddiff[j]
        end
    end
    return Zh
end        

##==========================================================


function comp_postm(everynit::Int64,firstwork::Int64,
                    iv::Array{Int64,1},lv::Array{Int64,1},jv::Array{Int64,1},
                    mv::Array{Int64,1},kv::Array{Int64,1},nv::Array{Int64,1},
                    U1::Array{Float64,2},U2::Array{Float64,2},U3::Array{Float64,2},
                    diaginvlambda::Array{Float64,1},Zh::Array{Float64,1},
                    mprior::Array{Float64,1},
                    astart::Int64,aend::Int64)
    ## compute Zh
    @assert aend>=astart
    startt = time()
    Na = length(mprior)
    myNa = aend-astart+1
    postm = zeros(Float64,myNa)
    elUDZh = zeros(Float64,myNa)
    # loop on chunk
    ### need to re-loop because full Zh is needed
    for i=astart:aend  #i=1:Na
        mya = i-astart+1

        # print info
        if (mya%everynit==0) | (mya==5)
            if myid()==firstwork
                eta = ( (time()-startt)/float(mya-1) * (myNa-mya+1) ) /60.0
                reta = round(eta,digits=3)
                println("[parallel] loop 3/3 a: $mya of $myNa; ETA: $reta min ") #  \r")
                #flush(stdout)
            end
        end


        ## UD times Zh
        elUDZh[mya] = 0.0
        for j=1:Na
            # element of row of UD
            elrowUD = U1[iv[i],iv[j]] * U2[jv[i],jv[j]] *
                U3[kv[i],kv[j]] * diaginvlambda[j]
            # element of final vector
            elUDZh[mya] = elUDZh[mya] + elrowUD * Zh[j]
        end

        ## element of the posterior mean
        postm[mya] = mprior[i] + elUDZh[mya] # sum(bigmatrow.*ddiff)
    end
    return postm
end

##===============================================================

function spreadwork(nit::Int64,nunits::Int64,startpoint::Int64)

    if (nit<=nunits)
        println(myid(),": spredwork(): nit<=nunits")
        println(nit,nunits)
        error()
    end 
    # allocate
    scheduling = Array{Int64,1}(undef,nunits)
    looping = Array{Int64,2}(undef,nunits,2)

    ## compute distribution of workload for Nb
    nitcpu = div(nit,nunits) # nit/nunits
    scheduling[1:end] .= nitcpu
    resto = mod(nit,nunits)
    
    ## spread the remaining workload
    scheduling[1:resto] = scheduling[1:resto] .+ 1
    looping[1,:] = [startpoint,startpoint+scheduling[1]-1] 
    for i=2:nunits
        looping[i,1] = sum(scheduling[1:i-1]) + startpoint
        looping[i,2] = sum(scheduling[1:i]) + startpoint - 1
    end 
    
    #!!print*,myrank,": loop 1",looping(:,1),&
    ## "| loop 2",looping(:,2)," | sched ",scheduling
    return scheduling,looping
end 

##==========================================================

##==========================================================
end # module
#########################################
