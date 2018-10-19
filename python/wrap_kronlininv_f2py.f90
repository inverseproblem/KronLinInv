

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


module wrapforpy

  use realprec
  use kronlininv
  use iso_c_binding, only : c_float,c_int,c_double, c_ptr
  implicit none

contains

  !========================================================

  subroutine c_calcfactors(nd1,nd2,nd3,nm1,nm2,nm3,&
       G1,G2,G3,Cm1,Cm2,Cm3,Cd1,Cd2,Cd3,&
       U1,U2,U3,diaginvlambda,iUCm1,iUCm2,iUCm3,&
       iUCmGtiCd1,iUCmGtiCd2,iUCmGtiCd3 )  bind(c) 
    integer(c_int),intent(in) :: nd1,nd2,nd3,nm1,nm2,nm3
    real(c_double) :: Cm1(nm1,nm1)
    real(c_double),intent(in) :: G1(nd1,nm1),G2(nd2,nm2),G3(nd3,nm3),&
         Cm2(nm2,nm2),Cm3(nm3,nm3),&
         Cd1(nd1,nd1),Cd2(nd2,nd2),Cd3(nd3,nd3)    
    real(c_double),intent(out) :: U1(nm1,nm1),U2(nm2,nm2),U3(nm3,nm3),&
         iUCm1(nm1,nm1),iUCm2(nm2,nm2),iUCm3(nm3,nm3),&
         iUCmGtiCd1(nm1,nm1),iUCmGtiCd2(nm2,nm2),iUCmGtiCd3(nm3,nm3),&
         diaginvlambda(nm1*nm2*nm3)
    call calcfactors(G1,G2,G3,Cm1,Cm2,Cm3,Cd1,Cd2,Cd3,&
         U1,U2,U3,diaginvlambda,iUCm1,iUCm2,iUCm3,&
         iUCmGtiCd1,iUCmGtiCd2,iUCmGtiCd3 )
  end subroutine c_calcfactors

  !========================================================

  subroutine c_posteriormean(nd1,nd2,nd3,nm1,nm2,nm3,U1,U2,U3, &
       diaginvlambda, Z1,Z2,Z3,&
       G1,G2,G3, mprior, dobs, postm) bind(c)
    integer(c_int),intent(in) :: nd1,nd2,nd3,nm1,nm2,nm3
    real(c_double),intent(in) :: U1(nm1,nm1),U2(nm2,nm2),U3(nm3,nm3),&
         G1(nd1,nm1),G2(nd2,nm2),G3(nd3,nm3),&
         diaginvlambda(nm1*nm2*nm3),Z1(nm1,nm1),Z2(nm2,nm2),Z3(nm3,nm3),&
         mprior(nm1*nm2*nm3),dobs(nd1*nd2*nd3)
    real(c_double),intent(out) :: postm(nm1*nm2*nm3)
    call posteriormean(U1,U2,U3, diaginvlambda, Z1,Z2,Z3,&
       G1,G2,G3, mprior, dobs, postm)
  end subroutine c_posteriormean

  !========================================================

  subroutine c_blockpostcov(nm1,nm2,nm3,&
       U1,U2,U3, diaginvlambda, &
       iUCm1,iUCm2,iUCm3, astart,aend,bstart,bend, postC) bind(c)
    integer(c_int),intent(in) :: nm1,nm2,nm3
    real(c_double),intent(in) :: U1(nm1,nm1),U2(nm2,nm2),U3(nm3,nm3),&
         iUCm1(nm1,nm1),iUCm2(nm2,nm2),iUCm3(nm3,nm3)
    real(c_double),intent(in) :: diaginvlambda(nm1*nm2*nm3) 
    integer(c_int),intent(in) :: astart,aend,bstart,bend
    real(c_double),intent(out) :: postC(aend-astart+1,bend-bstart+1)
    integer :: astart2,aend2,bstart2,bend2
    astart2=astart+1
    aend2=aend+1
    bstart2=bstart+1
    bend2=bend+1
    print*,"Implicit conversion of indices from Python to Fortran"
    print*,'astart:',astart2,'aend:',aend2
    print*,'bstart:',bstart2,'bend:',bend2
    call blockpostcov(U1,U2,U3, diaginvlambda, &
         iUCm1,iUCm2,iUCm3, astart2,aend2,bstart2,bend2, postC)
  end subroutine c_blockpostcov
  
  !========================================================
  
end module wrapforpy
