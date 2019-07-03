
!------------------------------------------------------------------------
!
!    Copyright 2018, Andrea Zunino 
!
!    This file is part of KronLinInv.
!
!    KronLinInv is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    KronLinInv is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with KronLinInv.  If not, see <http://www.gnu.org/licenses/>.
!
!------------------------------------------------------------------------


!> @section Serial/shared memory (OpenMP) version of KronLinInv
!>
!> @brief This file contains the serial/parallel on shared memory (OpenMP) version of KronLinInv.
!>  LAPACK libraries are required to be installed in the system. Optionally, OpenMP
!>  features from the compiler are needed, depending on compilation options.



!-------------------------------------------------------------
!> @author
!> Andrea Zunino \n
!> Niels Bohr Institute, University of Copenhagen \n
!>
!> @brief Sets precision of real variables
!
!-------------------------------------------------------------
module realprec

  ! http://fortranwiki.org/fortran/show/Real+precision
  integer,parameter  :: pdigits=15  !15 !! used also by MPI
  integer,parameter  :: rexprange=307  !307 !! used also by MPI
    
  ! RESULT = SELECTED_REAL_KIND([P, R, RADIX]) 
  integer,parameter  :: dp = selected_real_kind(pdigits, rexprange)

  ! integer, parameter :: sp = selected_real_kind(6, 37)
  ! integer, parameter :: dp = selected_real_kind(15, 307)
  ! integer, parameter :: qp = selected_real_kind(33, 4931)

  private :: pdigits,rexprange
  public  :: dp

end module realprec


!################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!################################################################

!-------------------------------------------------------------------------
!> @author
!> Andrea Zunino \n
!> Niels Bohr Institute, University of Copenhagen \n
!
!
!> @brief Procedures to perform linear inversion under gaussian assumptions using
!>  the Kronecker-product approach.
!
!-------------------------------------------------------------------------
module kronlininv
  
  use realprec
  use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT

  implicit none
!#ifdef omppara  
  include 'omp_lib.h'
!#endif

  
  private :: solvels,symgeneigvv,symsolvels
  public :: calcfactors,blockpostcov,posteriormean,bandpostcov
  
contains

  !!!=====================================================
  !--------------------------------------------------
  ! DESCRIPTION: 
  !> @brief
  !>  <b> Computes the factors necessary to solve the inverse problem. </b>
  !> 
  !> @detail
  !> The factors are the ones to be stored to subsequently calculate posterior
  !> mean and covariance. First an eigen decomposition is performed, to get
  !>
  !> \f[
  !>  \mathbf{U}_1 \mathbf{\Lambda}_1  \mathbf{U}_1^{-1}  
  !>  = \mathbf{C}_{\rm{M}}^{\rm{x}} (\mathbf{G}^{\rm{x}})^{\sf{T}}
  !> (\mathbf{C}_{\rm{D}}^{\rm{x}})^{-1} \mathbf{G}^{\rm{x}}
  !> \f]
  !>
  !> \f[
  !>\mathbf{U}_2 \mathbf{\Lambda}_2  \mathbf{U}_2^{-1}
  !> =  \mathbf{C}_{\rm{M}}^{\rm{y}} (\mathbf{G}^{\rm{y}})^{\sf{T}}
  !> (\mathbf{C}_{\rm{D}}^{\rm{y}})^{-1} \mathbf{G}^{\rm{y}}
  !> \f]
  !>
  !> \f[
  !> \mathbf{U}_3 \mathbf{\Lambda}_3  \mathbf{U}_3^{-1}
  !> = \mathbf{C}_{\rm{M}}^{\rm{z}} (\mathbf{G}^{\rm{z}})^{\sf{T}}
  !> (\mathbf{C}_{\rm{D}}^{\rm{z}})^{-1} \mathbf{G}^{\rm{z}} 
  !> \f]
  !>
  !> The principal factors involved in the computation of the posterior covariance and
  !> mean are
  !>
  !> \f[
  !> F_{\sf{A}} =  \mathbf{U}_1 \otimes \mathbf{U}_2 \otimes \mathbf{U}_3 
  !> \f]
  !>
  !> \f[
  !> F_{\sf{B}} = \big( 
  !> \mathbf{I} + \mathbf{\Lambda}_1 \! \otimes \! \mathbf{\Lambda}_2 \!
  !> \otimes \! \mathbf{\Lambda}_3 
  !> \big)^{-1} 
  !> \f]
  !>
  !> \f[
  !> F_{\sf{C}} =
  !> \mathbf{U}_1^{-1}  \mathbf{C}_{\rm{M}}^{\rm{x}} \otimes 
  !> \mathbf{U}_2^{-1} \mathbf{C}_{\rm{M}}^{\rm{y}} \otimes 
  !> \mathbf{U}_3^{-1} \mathbf{C}_{\rm{M}}^{\rm{z}} 
  !> \f]
  !>
  !> \f[
  !> F_{\sf{D}} =   
  !> \left( \mathbf{U}_1^{-1}  \mathbf{C}_{\rm{M}}^{\rm{x}} (\mathbf{G}^{\rm{x}})^{\sf{T}}
  !> (\mathbf{C}_{\rm{D}}^{\rm{x}})^{-1} \right) \!    \otimes 
  !> \left( \mathbf{U}_2^{-1} \mathbf{C}_{\rm{M}}^{\rm{y}}  (\mathbf{G}^{\rm{y}})^{\sf{T}}
  !> (\mathbf{C}_{\rm{D}}^{\rm{y}})^{-1}  \right)   \!   
  !> \otimes  \left( \mathbf{U}_3^{-1} \mathbf{C}_{\rm{M}}^{\rm{z}}
  !> (\mathbf{G}^{\rm{z}})^{\sf{T}} (\mathbf{C}_{\rm{D}}^{\rm{z}})^{-1} \right)
  !>  \f]
  !>
  !> @detail http://www.netlib.org/lapack/lug/node54.html
  !> ?sygst
  !> Reduces a real symmetric-definite generalized eigenvalue problem to
  !> the standard form. \n
  !> \f$B A z = \lambda z\f$ 	B = LLT 	C = LT A L 	z = L y
  !> - A is symmetric
  !> - B is symmetric, positive definite
  !> 
  !>
  !>
  !>
  !>
  !> This routines calculates
  !> @param[in] G1,G2,G3  The 3 forward modeling matrices \f$ \mathbf{G} =
  !>        \mathbf{G_1} \otimes \mathbf{G_2} \otimes \mathbf{G_3} \f$
  !> @param[in] Cm1,Cm2,Cm3,Cd1,Cd2,Cd3 The 6 covariance matrices
  !> @param[out] U1,U2,U3  \f$ \mathbf{U}_1 \f$, \f$ \mathbf{U}_2 \f$,
  !>     \f$ \mathbf{U}_3  \f$  of  \f$ F_{\sf{A}} \f$ 
  !> @param[out] diaginvlambda  \f$ F_{\sf{B}} \f$
  !> @param[out] iUCm1,iUCm2,iUCm3  \f$\mathbf{U}_1^{-1} \mathbf{C}_{\rm{M}}^{\rm{x}}\f$,
  !>      \f$\mathbf{U}_2^{-1}  \mathbf{C}_{\rm{M}}^{\rm{y}}\f$,
  !>       \f$\mathbf{U}_2^{-1}  \mathbf{C}_{\rm{M}}^{\rm{z}}\f$ of  \f$ F_{\sf{C}} \f$ 
  !> @param[out] iUCmGtiCd1,iUCmGtiCd1,iUCmGtiCd1  \f$ \mathbf{U}_1^{-1}
  !>       \mathbf{C}_{\rm{M}}^{\rm{x}}
  !>       (\mathbf{G}^{\rm{x}})^{\sf{T}}(\mathbf{C}_{\rm{D}}^{\rm{x}})^{-1}  \f$,
  !>       \f$ \mathbf{U}_2^{-1} \mathbf{C}_{\rm{M}}^{\rm{y}}
  !>       (\mathbf{G}^{\rm{y}})^{\sf{T}} (\mathbf{C}_{\rm{D}}^{\rm{y}})^{-1}\f$,
  !>       \f$ \mathbf{U}_3^{-1} \mathbf{C}_{\rm{M}}^{\rm{z}}
  !>       (\mathbf{G}^{\rm{z}})^{\sf{T}} (\mathbf{C}_{\rm{D}}^{\rm{z}})^{-1} \f$
  !>       of  \f$ F_{\sf{D}} \f$ 
  !>
  !>
  !  REVISION HISTORY:
  !> @date 10/8/2016 - Initial Version
  !
  !--------------------------------------------------
  subroutine calcfactors(G1,G2,G3,Cm1,Cm2,Cm3,Cd1,Cd2,Cd3,&
       U1,U2,U3,diaginvlambda,iUCm1,iUCm2,iUCm3,&
       iUCmGtiCd1,iUCmGtiCd2,iUCmGtiCd3 ) 
 
    real(dp),intent(in) :: G1(:,:),G2(:,:),G3(:,:),&
         Cm1(:,:),Cm2(:,:),Cm3(:,:),&
         Cd1(:,:),Cd2(:,:),Cd3(:,:)
    
    real(dp),intent(out) :: U1(:,:),U2(:,:),U3(:,:),&
         iUCm1(:,:),iUCm2(:,:),iUCm3(:,:),&
         iUCmGtiCd1(:,:),iUCmGtiCd2(:,:),iUCmGtiCd3(:,:),&
         diaginvlambda(:)
    
    real(dp),allocatable :: lambda1(:),lambda2(:),lambda3(:)
    real(dp),allocatable :: iCdG1(:,:),iCdG2(:,:),iCdG3(:,:),&
         GtiCdG1(:,:),GtiCdG2(:,:),GtiCdG3(:,:),&
         ide1(:,:),ide2(:,:),ide3(:,:),&
         iCd1(:,:),iCd2(:,:),iCd3(:,:)
    integer :: i,j,k,p,nm,nm1,nm2,nm3,nd1,nd2,nd3
    character :: uplo
   
    !! Lapack routines
    external DSYGV

    nm1 = size(Cm1,1)
    nm2 = size(Cm2,1)
    nm3 = size(Cm3,1)
    nd1 = size(Cd1,1)
    nd2 = size(Cd2,1)
    nd3 = size(Cd3,1)
    nm = nm1*nm2*nm3

    !! compute preliminary stuff
    print*,'calcfactors(): compute preliminary stuff'
    allocate(iCdG1(nd1,nm1),iCdG2(nd2,nm2),iCdG3(nd3,nm3))
    call symsolvels( Cd1, G1, iCdG1 )
    call symsolvels( Cd2, G2, iCdG2 )
    call symsolvels( Cd3, G3, iCdG3 )

    !! avoid automatic reallocation from matmul...
    !! so the code works with -fno-realloc-lhs in gfrotran
    allocate(GtiCdG1(nm1,nm1),GtiCdG2(nm2,nm2),GtiCdG3(nm3,nm3))
    GtiCdG1(:,:)  = matmul( transpose(G1), iCdG1 )
    GtiCdG2(:,:)  = matmul( transpose(G2), iCdG2 ) 
    GtiCdG3(:,:)  = matmul( transpose(G3), iCdG3 )

    !!--------------------------
    !!          fa
    !!--------------------------
    !! compute eigendecomposition
    print*,'calcfactors(): compute fa'
    allocate(lambda1(nm1),lambda2(nm2),lambda3(nm3))
    uplo = 'L'
    call symgeneigvv(GtiCdG1,uplo,Cm1,lambda1,U1)
    call symgeneigvv(GtiCdG2,uplo,Cm2,lambda2,U2)
    call symgeneigvv(GtiCdG3,uplo,Cm3,lambda3,U3)

    !!--------------------------
    !!          fb
    !!--------------------------
    !! calculates diagonal central factor
    print*,'calcfactors():compute fb'
    p=1
    do i=1,nm1
       do j=1,nm2
          do k=1,nm3
             diaginvlambda(p) = 1.0_dp/(1.0_dp+lambda1(i)*lambda2(j)*lambda3(k))
             p = p+1
          end do 
       end do
    end do

    !!--------------------------
    !!          fc
    !!--------------------------
    !! computes U1^-1 Cmx, U2^-1 Cmy, U3^-1 Cmz
    print*,'calcfactors(): compute fc'
    call solvels( U1, Cm1, iUCm1 )
    call solvels( U2, Cm2, iUCm2 )
    call solvels( U3, Cm3, iUCm3 )
    
    !!--------------------------
    !!          fd
    !!--------------------------
    print*,'calcfactors(): compute fd'
    allocate(ide1(nd1,nd1),ide2(nd2,nd2),ide3(nd3,nd3))
    ide1 = 0.0_dp  ! Initialize the array.
    ide2 = 0.0_dp
    ide3 = 0.0_dp
    forall(i = 1:nd1) ide1(i,i) = 1.0_dp  ! Set the diagonal.
    forall(i = 1:nd2) ide2(i,i) = 1.0_dp 
    forall(i = 1:nd3) ide3(i,i) = 1.0_dp
    allocate(iCd1(nd1,nd1),iCd2(nd2,nd2),iCd3(nd3,nd3))
    call symsolvels( Cd1, ide1, iCd1 )
    call symsolvels( Cd2, ide2, iCd2 )
    call symsolvels( Cd3, ide3, iCd3 )

    iUCmGtiCd1 = matmul( iUCm1, matmul( transpose(G1), iCd1 ))     
    iUCmGtiCd2 = matmul( iUCm2, matmul( transpose(G2), iCd2 )) 
    iUCmGtiCd3 = matmul( iUCm3, matmul( transpose(G3), iCd3 ))

    ! print*,
    ! print*,shape(U1),shape(U2),shape(U3)
    ! print*,shape(iUCm1),shape(iUCm2),shape(iUCm3)
    ! print*,shape(iUCmGtiCd1),shape(iUCmGtiCd2),shape(iUCmGtiCd3)
    ! print*,
    ! print*,'calcfactors(): end '
  end subroutine calcfactors

  !!!==================================================
  !-------------------------------------------------
  !>  @brief <b> Computes a block of the posterior covariance. </b>
  !
  !> @param[in] U1, U2, U3 \f$ \mathbf{U}_1 \f$, \f$ \mathbf{U}_2 \f$,
  !>     \f$ \mathbf{U}_3  \f$  of  \f$ F_{\sf{A}} \f$ 
  !> @param[in] diaginvlambda  \f$ F_{\sf{B}} \f$
  !> @param[in] iUCm1, iUCm2, iUCm3  \f$\mathbf{U}_1^{-1} \mathbf{C}_{\rm{M}}^{\rm{x}}\f$,
  !>      \f$\mathbf{U}_2^{-1}  \mathbf{C}_{\rm{M}}^{\rm{y}}\f$,
  !>       \f$\mathbf{U}_2^{-1}  \mathbf{C}_{\rm{M}}^{\rm{z}}\f$ of  \f$ F_{\sf{C}} \f$ 
  !> @param[in] astart,aend index of first/last row of the requested block
  !> @param[in] bstart,bend index of first/last column of the requested block
  !> @param[out] postC block of the posterior covariance
  !
  !
  ! REVISION HISTORY:
  !>  @date 5/8/2016 - Initial Version
  !>  @date 4/11/2016 - Removed indices lv,mv,kv
  !>  @date 27/13/2016 - New ETA calculations
  !
  !-------------------------------------------------  
  subroutine blockpostcov(U1,U2,U3, diaginvlambda, &
       iUCm1,iUCm2,iUCm3, astart,aend,bstart,bend, postC) 
    !!
    !! Calculates a block of the posterior covariance
    !! 
    real(dp),intent(in) :: U1(:,:),U2(:,:),U3(:,:),iUCm1(:,:), &
         iUCm2(:,:),iUCm3(:,:)
    !! diaginvlambda = (I + lam1 x lam2 x lam3 )^-1
    real(dp),intent(in) :: diaginvlambda(:) !! diagonal/vector central factor
    real(dp),intent(out) :: postC(:,:)
    integer,intent(in) :: astart,aend,bstart,bend
    
    integer :: a,b
    real(dp),allocatable :: row1(:),row2(:),col1(:)
    integer :: Nr12,Nc1
    integer :: Ni,Nj,Nk,Nl,Nm,Nn,Na,Nb    

    integer,allocatable :: av(:)!,bv(:)
    integer,allocatable :: iv(:),jv(:),kv(:) !,lv(:),mv(:),nv(:) 
    integer :: p,csha1,csha2,everynit,nthr,privcount,totcount 
    real(dp) :: eta,frac,startt,endt
    
    Ni = size(U1,1)
    Nl = size(U1,2)
    Nj = size(U2,1)
    Nm = size(U2,2)
    Nk = size(U3,1)
    Nn = size(U3,2)
    Na = Ni*Nj*Nk
    Nb = Nl*Nm*Nn
    
    !! check shape of output array
    if (Na /= Nb) then
       write(*,*) '(Na /= Nb)', Na,Nb
       stop
    end if
    csha1 = aend-astart+1
    csha2 = bend-bstart+1
    if ( (size(postC,1)/=csha1) .or. (size(postC,2)/=csha2)) then
       write(*,*) "Wrong size of the posterior covariance block array."
       stop
    end if
    !! check limits of requested block
    if ( (astart<1) .or. (aend>Na) .or. (astart>aend) .or. &
         (bstart<1) .or. (bend>Na) .or. (bstart>bend) ) then
       write(*,*) "Wrong size of the requested block array."
       stop
    end if
        
    !! vectorize row and col calculations for Kron prod AxBxC
    allocate(av(Na),iv(Na),jv(Na),kv(Na))
    !allocate(bv(Nb),lv(Nb),mv(Nb),nv(Nb))
    forall(p = 1:Na) av(p) = p
    !forall(p = 1:Nb) bv(p) = p

    !! vectors containing all possible indices for 
    !!    row calculations of Kron prod AxBxC
    iv(:) =  (av-1)/(Nk*Nj)+1 
    jv(:) =  (av-1-(iv-1)*Nk*Nj)/Nk + 1 
    kv(:) =  av-(jv-1)*Nk-(iv-1)*Nk*Nj 
     
    !! allocate stuff
    Nr12 = size(U1,2)* size(U2,2)* size(U3,2)
    Nc1  = size(iUCm1,1)* size(iUCm2,1)* size(iUCm3,1)
    allocate(row1(Nr12),row2(Nr12),col1(Nc1))

    !!------------------------
    if (aend-astart+1<20) then
       everynit = 1
    else if (aend-astart+1<100000) then
       everynit = (aend-astart+1)/100
    else
       everynit = (aend-astart+1)/1000
    end if

    
    !!------------------------
    call cpu_time(startt) 

    totcount = aend-astart+1
    privcount = -1
    do a=astart,aend
       privcount = privcount + 1 
       if  (mod(privcount,everynit)==0) then
          call cpu_time(endt)
          frac = real(privcount,dp)/(real(totcount,dp))
          eta = ( real(endt-startt,dp) / real(privcount,dp) ) * &
               (real(totcount-privcount,dp))
          write(OUTPUT_UNIT,fmt='(a19,f7.3,6x,a4,f12.2,1x,a3,a5)') 'blockpostcov():  %',&
               frac*100_dp,"ETA:",eta/60.0,"min",char(27)//'[1A'//achar(13)
          flush(OUTPUT_UNIT) !! to make sure it prints immediately          
       end if

       !! row first two factors
       !! a row x diag matrix 
       !!row2 =  U1(iv(a),lv) * U2(jv(a),mv) * U3(kv(a),nv) * diaginvlambda
       row2(:) = U1(iv(a),iv) * U2(jv(a),jv) * U3(kv(a),kv) * diaginvlambda
       

       !$OMP PARALLEL 
       !$OMP PARALLEL DO PRIVATE(b,col1)
       do b=bstart,bend
          
          !! calculate one row of first TWO factors
          !!call columnAxBxC(Ni,Nj,Nk,Nl,Nm,Nn, iUCm1,iUCm2,iUCm3,b,col1)
          
          col1(:) = iUCm1(iv,iv(b)) * iUCm2(jv,jv(b)) * iUCm3(kv,kv(b))
          
          !! calculate one element of the posterior covariance
          postC(a,b) = sum(row2*col1)
          
       end do       
       !! !$OMP END DO
       !$OMP END PARALLEL

    end do
    call cpu_time(endt)
    write(OUTPUT_UNIT,*)
    print*,"blockpostcov():",endt-startt,"s elapsed"

  end subroutine blockpostcov

  !!!==================================================
  !-------------------------------------------------
  !>  @brief <b> Computes a band of the posterior covariance. </b>
  !> See http://www.netlib.org/lapack/lug/node124.html
  !> @param[in] U1,U2,U3  \f$ \mathbf{U}_1 \f$, \f$ \mathbf{U}_2 \f$,
  !>     \f$ \mathbf{U}_3  \f$  of  \f$ F_{\sf{A}} \f$ 
  !> @param[in] diaginvlambda  \f$ F_{\sf{B}} \f$
  !> @param[in] iUCm1,iUCm2,iUCm3  \f$\mathbf{U}_1^{-1} \mathbf{C}_{\rm{M}}^{\rm{x}}\f$,
  !>      \f$\mathbf{U}_2^{-1}  \mathbf{C}_{\rm{M}}^{\rm{y}}\f$,
  !>       \f$\mathbf{U}_2^{-1}  \mathbf{C}_{\rm{M}}^{\rm{z}}\f$ of  \f$ F_{\sf{C}} \f$ 
  !> @param[in] lowdiag,updiag  Lower and upper diagonal number of requested band
  !> @param[out] postC  band of the posterior covariance stored following Lapack convention
  !
  !
  ! REVISION HISTORY:
  !>  @date 11/8/2016 - Initial Version
  !>  @date 4/11/2016 - Removed indices lv,mv,kv
  !
  !-------------------------------------------------  
  subroutine bandpostcov(U1,U2,U3, diaginvlambda, &
       iUCm1,iUCm2,iUCm3, lowdiag, updiag, bandpostC) 
    !!
    !! Calculate a band of the posterior covariance
    !! 
    real(dp),intent(in) :: U1(:,:),U2(:,:),U3(:,:),iUCm1(:,:), &
         iUCm2(:,:),iUCm3(:,:)
    !! diaginvlambda = (I + lam1 x lam2 x lam3 )^-1
    real(dp),intent(in) :: diaginvlambda(:) !! diagonal/vector central factor
    real(dp),intent(inout) :: bandpostC(:,:)
    integer,intent(in) :: lowdiag, updiag
    
    integer :: a,b
    real(dp),allocatable :: row1(:),row2(:),col1(:)
    integer :: Nr12,Nc1
    integer :: Ni,Nj,Nk,Nl,Nm,Nn,Na,Nb    

    integer,allocatable :: av(:)!,bv(:)
    integer,allocatable :: iv(:),jv(:),kv(:) !,lv(:),mv(:),nv(:) 
    integer :: p,aband,aend,bband,astart,d
    integer :: nthr,privcount,everynit,totcount
    real(dp) :: eta,frac,startt,firststartt,endt    
    
    Ni = size(U1,1)
    Nl = size(U1,2)
    Nj = size(U2,1)
    Nm = size(U2,2)
    Nk = size(U3,1)
    Nn = size(U3,2)
    Na = Ni*Nj*Nk
    Nb = Nl*Nm*Nn    

    if (Na /= Nb) then
       write(*,*) '(Na /= Nb)', Na,Nb
       stop
    end if
    if ( (updiag>=Na) .or. (lowdiag>=Na) .or. (lowdiag<0) .or. (updiag<0) ) then
       write(*,*) "(updiag<Na) .or. (lowdiag<Na)"
       write(*,*) "updiag",updiag,"Na",Na,"lowdiag",lowdiag,"Na",Na
       stop
    end if

    !! vectorize row and col calculations for Kron prod AxBxC
    allocate(av(Na),iv(Na),jv(Na),kv(Na))
    !allocate(bv(Nb),lv(Nb),mv(Nb),nv(Nb))
    forall(p = 1:Na) av(p) = p
    !forall(p = 1:Nb) bv(p) = p

    !! vectors containing all possible indices for 
    !!    row calculations of Kron prod AxBxC
    iv(:) =  (av-1)/(Nk*Nj)+1 
    jv(:) =  (av-1-(iv-1)*Nk*Nj)/Nk+1 
    kv(:) =  av-(jv-1)*Nk-(iv-1)*Nk*Nj 
    !! vectors containing all possible indices for
    !!    column calculations of Kron prod AxBxC
    ! lv =  (bv-1)/(Nn*Nm) + 1 
    ! mv =  (bv-1-(lv-1)*Nn*Nm)/Nn + 1 
    ! nv =  bv-(mv-1)*Nn-(lv-1)*Nn*Nm
    
    !! allocate stuff
    Nr12 = size(U1,2)* size(U2,2)* size(U3,2)
    Nc1  = size(iUCm1,1)* size(iUCm2,1)* size(iUCm3,1)
    allocate(row1(Nr12),row2(Nr12),col1(Nc1))
   
    !!everynit = 250
    
    ! Lapack: http://www.netlib.org/lapack/lug/node124.html
    ! aij is stored in AB(ku+1+i-j,j) for max(1,j-ku) <= i <= \min(m,j+kl).
    !------------
    ! Diagonals of a matrix
    ! i + d = j
    ! main diag d = 0
    ! upper d > 0
    ! lower d < 0
    ! diagonals of a matrix and indices of related band matrix
    ! ONLY for square matrix
   
    ! initialize postC
    bandpostC = 0.0_dp

    !!=====================================
    firststartt = omp_get_wtime()
    
    totcount = lowdiag+updiag+1
    
    do d=-lowdiag,updiag
       if (d<0) then
          astart = abs(d)+1
       else
          astart = 1
       end if
       if (d>0) then
          aend = Na-abs(d)
       else
          aend = Na
       end if
       !print*,'diagonal',d
       !! indices of normal matrix
       
       !!print*,"Diagonal ",d," from range [",-lowdiag,",",updiag,"]"
       !-------------------------
       !$OMP PARALLEL
       !privcount = -1
       startt = omp_get_wtime()
       nthr = omp_get_num_threads()
       
       if ( (omp_get_thread_num()==0 ) ) then 
          frac = real(d+lowdiag,dp)/real(totcount,dp)
          eta = (omp_get_wtime()-firststartt)/frac * (1.0_dp - frac)
          write(OUTPUT_UNIT,fmt='(a19,f6.3,6x,a4,f12.2,1x,a3,a5)') ' bandpostcov():  % ',&
               frac*100_dp,"ETA:",eta/60.0,"min",char(27)//'[1A'//achar(13)
          flush(OUTPUT_UNIT) !! to make sure it prints immediately
       end if

       
       !$OMP DO PRIVATE(row2,col1,a,b,aband,bband) !!,privcount)
       do a=astart,aend
          !privcount = privcount + 1 
          
          b = a+d
          !! indices of the band matrix
          aband = updiag+1+a-b
          bband = b

          !! row first two factors
          row2(:) = diaginvlambda * U1(iv(a),iv) * U2(jv(a),jv) * U3(kv(a),kv)

          !! calculate one row of first TWO factors
          !!call columnAxBxC(Ni,Nj,Nk,Nl,Nm,Nn, iUCm1,iUCm2,iUCm3,b,col1)
          col1(:) = iUCm1(iv,iv(b)) * iUCm2(jv,jv(b)) * iUCm3(kv,kv(b))

          !! calculate one element of the posterior covariance
          !! store it in the band storage format
          bandpostC(aband,bband) = sum(row2*col1)
          !print*, a,b,aband,bband,bandpostC(aband,bband)

       end do ! a=astart,aend
       !$OMP END DO
       endt = omp_get_wtime()
       !$OMP END PARALLEL
    
    end do
    write(OUTPUT_UNIT,*)
    print*,"bandpostcov():",endt-firststartt," OMP wall clock time"
    
    
  end subroutine bandpostcov

  !!!===============================================================
  !-------------------------------------------------
  !>  @brief <b> Computes the posterior mean  </b>
  !
  !> @param[in] U1,U2,U3  \f$ \mathbf{U}_1 \f$, \f$ \mathbf{U}_2 \f$,
  !>     \f$ \mathbf{U}_3  \f$  of  \f$ F_{\sf{A}} \f$ 
  !> @param[in] diaginvlambda  \f$ F_{\sf{B}} \f$
  !> @param[in] Z1,Z2,Z3  \f$ \mathbf{U}_1^{-1}
  !>       \mathbf{C}_{\rm{M}}^{\rm{x}}
  !>       (\mathbf{G}^{\rm{x}})^{\sf{T}}(\mathbf{C}_{\rm{D}}^{\rm{x}})^{-1}  \f$,
  !>       \f$ \mathbf{U}_2^{-1} \mathbf{C}_{\rm{M}}^{\rm{y}}
  !>       (\mathbf{G}^{\rm{y}})^{\sf{T}} (\mathbf{C}_{\rm{D}}^{\rm{y}})^{-1}\f$,
  !>       \f$ \mathbf{U}_3^{-1} \mathbf{C}_{\rm{M}}^{\rm{z}}
  !>       (\mathbf{G}^{\rm{z}})^{\sf{T}} (\mathbf{C}_{\rm{D}}^{\rm{z}})^{-1} \f$
  !>       of  \f$ F_{\sf{D}} \f$
  !> @param[in] G1,G2,G3 The 3 forward modeling matrices \f$ \mathbf{G} =
  !>        \mathbf{G_1} \otimes \mathbf{G_2} \otimes \mathbf{G_3} \f$
  !> @param[in] mprior Prior model (vector)
  !> @param[in] dobs  Observed data (vector)
  !> @param[out] postm Calculated posterior mean model (vector)
  !
  !
  ! REVISION HISTORY:
  !>  @date 5/8/2016 - Initial Version
  !>  @date 27/12/2016 - New ETA calculations
  ! 
  !-------------------------------------------------
  subroutine posteriormean(U1,U2,U3, diaginvlambda, Z1,Z2,Z3,&
       G1,G2,G3, mprior, dobs, postm) 
    !!
    !! Calculate the posterior mean model
    !! 
    real(dp),intent(in)  :: U1(:,:),U2(:,:),U3(:,:), G1(:,:),G2(:,:),G3(:,:)
    real(dp),intent(in)  :: Z1(:,:),Z2(:,:),Z3(:,:)
    real(dp),intent(in)  :: mprior(:),dobs(:),diaginvlambda(:)
    real(dp),intent(out) :: postm(:)

    real(dp),allocatable :: row2(:),col1(:),bigmatrow(:),datrow(:)
    real(dp),allocatable :: ddiff(:)
    integer :: Nr12,Nc1,a,b,Na,Nb
    integer :: Ni,Nj,Nk,Nl,Nm,Nn

    integer,allocatable :: av(:),bv(:)
    integer,allocatable :: iv(:),jv(:),kv(:),lv(:),mv(:),nv(:)
        !! ivo(:),jvo(:),kvo(:)
    integer :: p,j,i,everynit

    real(dp),allocatable :: Zh(:),elUDZh(:)
    real(dp) :: datp,elg,elrowud,tZZ

    integer :: nsteps,hr,min,nthread,privcount,nthr
    real(dp) :: eta,et1,et2,sec,startt,firststartt,endt,finish,frac
    character(8)  :: curdate
    character(10) :: curtime

    write(OUTPUT_UNIT,*) "posteriormean(): calculating posterior mean... "
    
    !! sizes
    Ni = size(Z1,1)
    Nl = size(Z1,2)
    Nj = size(Z2,1)
    Nm = size(Z2,2)
    Nk = size(Z3,1)
    Nn = size(Z3,2)

    !! sizes
    ! Nr12 = size(U1,2)*size(U2,2)*size(U3,2)
    ! Nc1  = size(Z1,1)*size(Z2,1)*size(Z3,1)
    Na   = size(mprior,1)
    Nb   = size(dobs,1)
    
    !! allocate stuff
    allocate(ddiff(Nb))
    allocate(Zh(Na),elUDZh(Na))
 
    !! vectorize row and col calculations for Kron prod AxBxC
    !!allocate(ivo(Nb),jvo(Nb),kvo(Nb))
    allocate(av(Na),iv(Na),jv(Na),kv(Na))
    allocate(bv(Nb),lv(Nb),mv(Nb),nv(Nb))
    forall(p = 1:Na) av(p) = p
    forall(p = 1:Nb) bv(p) = p

    !! vectors containing all possible indices for 
    !!    row calculations of Kron prod AxBxC
    iv(:) =  (av-1)/(Nk*Nj)+1 
    jv(:) =  (av-1-(iv-1)*Nk*Nj)/Nk+1 
    kv(:) =  av-(jv-1)*Nk-(iv-1)*Nk*Nj 
    !! vectors containing all possible indices for
    !!    column calculations of Kron prod AxBxC
    lv(:) =  (bv-1)/(Nn*Nm) + 1 
    mv(:) =  (bv-1-(lv-1)*Nn*Nm)/Nn + 1 
    nv(:) =  bv-(mv-1)*Nn-(lv-1)*Nn*Nm 
    !!  Gs have different shape than Us !!

    !!#######################
    if (Na<1000)then
       everynit = Na/50
    else if (Na<100000) then
       everynit = Na/100
    else
       everynit = Na/1000
    end if

    
    !!---------------------------------
    !$OMP PARALLEL PRIVATE(privcount)
    firststartt = omp_get_wtime()
    nthr = omp_get_num_threads()

    !!#######################
    !!#    dobs - dcalc     #
    !!#######################
    startt= omp_get_wtime()
    privcount = -1
    !$OMP DO PRIVATE(b)
    !! difference obs-calc data
    do b=1,Nb
       privcount = privcount + 1 
       if ( (omp_get_thread_num()==0 ) .and. (mod(privcount,everynit/nthr)==0) ) then          
          frac = real(privcount,dp)/(real(Nb,dp)/real(nthr,dp))
          eta = ( (omp_get_wtime()-startt) / real(privcount,dp) ) * &
               (real(Nb-nthr*privcount,dp)/real(nthr,dp))
          write(OUTPUT_UNIT,fmt='(a28,f7.3,6x,a4,f12.2,1x,a3,a5)') ' posteriormean() loop 1/3: %',&
               frac*100_dp,"ETA:",eta/60.0,"min",char(27)//'[1A'//achar(13)
          flush(OUTPUT_UNIT) !! to make sure it prints immediately          
       end if
       
       !!---------------------------------------------------------------------------
       ddiff(b) = dobs(b) - sum(mprior * G1(lv(b),iv) * G2(mv(b),jv) * G3(nv(b),kv))
       !!---------------------------------------------------------------------------

    end do
    !$OMP END DO
       
    !!#######################
    !!#   === U d Z h ===   #
    !!#######################
    privcount = -1
    startt= omp_get_wtime()
    !$OMP DO PRIVATE(a)
    do a=1,Na
       privcount = privcount + 1 
       if ( (omp_get_thread_num()==0 ) .and. (mod(privcount,everynit/nthr)==0) ) then          
          frac = real(privcount,dp)/(real(Na,dp)/real(nthr,dp))
          eta = ( (omp_get_wtime()-startt) / real(privcount,dp) ) * &
               (real(Na-nthr*privcount,dp)/real(nthr,dp))
          write(OUTPUT_UNIT,fmt='(a28,f7.3,6x,a4,f12.2,1x,a3,a5)') ' posteriormean() loop 2/3: %',&
               frac*100_dp,"ETA:",eta/60.0,"min",char(27)//'[1A'//achar(13)
          flush(OUTPUT_UNIT) !! to make sure it prints immediately          
       end if

       !!---------------------------------------------------------------------------
       Zh(a) = sum( Z1(iv(a),lv) * Z2(jv(a),mv) * Z3(kv(a),nv) * ddiff )
       !!---------------------------------------------------------------------------
       
    end do
    !$OMP END DO
        
    !!-----------------------------------------------
    !! Second loop
    !!### need to re-loop because full Zh is needed
    !!#######################
    !!#       post(a)       #
    !!#######################
    privcount = -1
    startt= omp_get_wtime()
    !$OMP DO PRIVATE(a)
    do a=1,Na
       privcount = privcount + 1
       if ( (omp_get_thread_num()==0 ) .and. (mod(privcount,everynit/nthr)==0) ) then          
          frac = real(privcount,dp)/(real(Na,dp)/real(nthr,dp))
          eta = ( (omp_get_wtime()-startt) / real(privcount,dp) ) * &
               (real(Na-nthr*privcount,dp)/real(nthr,dp))
          write(OUTPUT_UNIT,fmt='(a28,f7.3,6x,a4,f12.2,1x,a3,a5)') ' posteriormean() loop 3/3: %',&
               frac*100_dp,"ETA:",eta/60.0,"min",char(27)//'[1A'//achar(13)
          flush(OUTPUT_UNIT) !! to make sure it prints immediately          
       end if

       !!--------------------------------------------------------------------------------
       elUDZh(a) = sum( U1(iv(a),iv) * U2(jv(a),jv) * U3(kv(a),kv) * diaginvlambda * Zh )
       !! element of the posterior mean
       postm(a) = mprior(a) + elUDZh(a)
       !!--------------------------------------------------------------------------------
       
    end do
    !$OMP END DO
    endt = omp_get_wtime()
    !$OMP END PARALLEL
    write(OUTPUT_UNIT,*)
    print*,"posteriormean():",endt-firststartt," OMP wall clock time"

  end subroutine posteriormean
  
  !!!==================================================
  !-------------------------------------------------
  !>  @brief <b> Computes eigenvalues and eigenvectors of
  !>    the generalized symmetric definite eigenproblem. </b>
  !>  See http://www.netlib.org/lapack/lug/node54.html
  !>
  !> No test to check symmetry or positive definiteness is performed.
  !>
  !> @param[in] A symmetric matrix
  !> @param[in] uplo upper or lower triangle stored
  !> @param[in] Bpd symmetric positive-definite matrix
  !> @param[out] lambda eigenvalues
  !> @param[out] U eigenvectors
  !
  !
  ! REVISION HISTORY:
  !>  @date 10/8/2016 - Initial Version
  ! 
  !-------------------------------------------------  
  subroutine symgeneigvv(A,uplo,Bpd,lambda,U)
    real(dp),intent(in) :: A(:,:),Bpd(:,:)
    character,intent(in) :: uplo
    real(dp),intent(out) :: lambda(:),U(:,:)
    integer :: n,lwork,itype,info
    character :: jobz
    real(dp),allocatable :: work(:),tmpB(:,:)
    ! DSYGV computes all the eigenvalues, and optionally, the eigenvectors
    ! of a real generalized symmetric-definite eigenproblem, of the form
    ! A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
    ! Here A and B are assumed to be symmetric and B is also
    ! positive definite.
    itype = 3
    jobz = 'V'
    n=size(A,1)
    if ( (size(A,2)/=n) .or. (size(Bpd,1)/=n) .or. (size(Bpd,2)/=n) ) then
       write(*,*) "symgeneigvv(): Matrices have wrong sizes"
       stop 
    end if
    allocate(tmpB(n,n))
    U = A
    tmpB(:,:) = Bpd
    lwork=3*n-1
    ! lwork=-1 
    ! allocate(work(1))
    ! call DSYGV(itype,jobz,uplo,nm1,U1,nm1,choCm1,nm1,lambda1,work,lwork,info)
    ! deallocate(work)
    allocate(work(lwork))
    call DSYGV(itype,jobz,uplo,n,U,n,tmpB,n,lambda,work,lwork,info)
    if (info /= 0) stop 'symgeneigvv(): Matrix eigendecomposition failed!'
  end subroutine symgeneigvv

  !!!==================================================
  !-------------------------------------------------
  !>  @brief <b> Solves a linear system AX = B,
  !>      real numbers </b>
  !
  !> @param[in] A coefficients matrix (linear operator)
  !> @param[in] B matrix
  !> @param[out] sol solution to the linear system (matrix)
  !
  !
  ! REVISION HISTORY:
  !>  @date 5/8/2016 - Initial Version
  ! 
  !-------------------------------------------------
  subroutine solvels( A, B, sol )
    real(dp),intent(in) :: A(:,:),B(:,:)
    real(dp),intent(out) :: sol(:,:)
    integer,allocatable :: ipiv(:)
    integer :: n,nrhs,lda,ldb
    integer :: info
    real(dp),allocatable :: A2(:,:)
    ! External procedures defined in LAPACK
    external DGESV
    n = size(A,1)
    nrhs = size(B,2)
    allocate(A2(n,n))
    A2(:,:) = A
    lda = size(A,1)
    allocate(ipiv(n))
    ldb = size(b,1)
    sol = B
    call DGESV(n,nrhs,A2,lda,ipiv,sol,ldb,info)
    if (info/=0) then
       write(*,*) "linear system solver failed...'"
       print*,'info: ',info
       stop
    end if
  end subroutine solvels

  !!!============================================================
  !-------------------------------------------------
  !>  @brief <b> Solves a linear system AX = B
  !>      for symmetric A, real numbers </b>
  !
  !> @param[in] A coefficients matrix (linear operator)
  !> @param[in] B matrix
  !> @param[out] sol solution to the linear system (matrix)
  !
  ! REVISION HISTORY:
  !>  @date 5/8/2016 - Initial Version
  ! 
  !-------------------------------------------------
  subroutine symsolvels( A, B, sol )
    real(dp),intent(in) :: A(:,:),B(:,:)
    real(dp),intent(out) :: sol(:,:)
    integer,parameter :: wmax=3000 ! max size allowed for work
    integer,allocatable :: ipiv(:)
    integer :: n,nrhs,lda,ldb,lwork
    character :: uplo
    real(dp),allocatable :: work(:),A2(:,:)
    integer :: info
    ! External procedures defined in LAPACK
    external DSYSV
    uplo = 'L'
    n = size(a,1)
    nrhs = size(b,2)
    allocate(A2(n,n))
    A2(:,:) = A
    lda = size(A,1)
    allocate(ipiv(n))
    ldb = size(b,1)
    ! copy B to sol to avoid it to be overwritten
    sol = B
    ! apparently work must be allocated for the query to dsysv to work...
    allocate(work(1)) 
    lwork = -1
    call DSYSV(uplo,n,nrhs,A2,lda,ipiv,sol,ldb,work,lwork,info)
    lwork = min(wmax,int(work(1)))
    !! now lwork is assigned, so actually calculate eigvec/val
    deallocate(work)
    allocate(work(lwork))
    call DSYSV(uplo,n,nrhs,A2,lda,ipiv,sol,ldb,work,lwork,info)
    if (info/=0) then
       write(*,*) "linear system solver failed...'"
       print*,'info: ',info
       stop
    end if
  end subroutine symsolvels

  !!====================================================
 
end module kronlininv

!!====================================================
!!====================================================
