!-------------------------------------------
!
!> @brief Test program for KronLinInv
!
!> @author
!> Andrea Zunino \n
!> Niels Bohr Institute, University of Copenhagen \n
!
!-------------------------------------------
program test

  use realprec
  use kronlininv
  use readwriteh5

  implicit none

  real(dp),allocatable :: G1(:,:),G2(:,:),G3(:,:),&
         Cm1(:,:),Cm2(:,:),Cm3(:,:),&
         Cd1(:,:),Cd2(:,:),Cd3(:,:),&
         mprior(:),dobs(:),mtrue(:)
  real(dp),allocatable :: U1(:,:),U2(:,:),U3(:,:),&
         iUCm1(:,:),iUCm2(:,:),iUCm3(:,:),&
         iUCmGtiCd1(:,:),iUCmGtiCd2(:,:),iUCmGtiCd3(:,:),&
         invlambda(:)  
  character(len=1024) :: inpfile,outfile
  integer :: nm,nd,nm1,nm2,nm3,nd1,nd2,nd3
  real(dp),allocatable :: postm(:),postcov(:,:),postcov_diag(:,:)
  real(dp),allocatable :: bandpostC(:,:)
  integer :: i1,i2,j1,j2, lowdiag, updiag
  real(dp) :: start,finish

  integer :: iar
  character(len=256) :: simname
  !----------------------------------
  
  !! simname -> name of the input file (without the .hdf5)

  iar = 1 
  call get_command_argument(iar, simname)

  inpfile = trim(simname)//'.h5'
  outfile = trim(simname)//'_output.h5'

  print*,"inpfile:",trim(inpfile)
  print*,"outfile:",trim(outfile)

  !----------------------------------
  
  !! read all input arrays
  print*,'reading input arrays'
  call readreal2Darrh5(inpfile,'G1',G1)
  call readreal2Darrh5(inpfile,'G2',G2)
  call readreal2Darrh5(inpfile,'G3',G3)
  call readreal2Darrh5(inpfile,'Cm1',Cm1)
  call readreal2Darrh5(inpfile,'Cm2',Cm2)
  call readreal2Darrh5(inpfile,'Cm3',Cm3)
  call readreal2Darrh5(inpfile,'Cd1',Cd1)
  call readreal2Darrh5(inpfile,'Cd2',Cd2)
  call readreal2Darrh5(inpfile,'Cd3',Cd3)
  call readreal1Darrh5(inpfile,'mtrue',mtrue)
  call readreal1Darrh5(inpfile,'dobs',dobs)
  call readreal1Darrh5(inpfile,'mprior',mprior)
  

  ! write arrays to output file
  call writereal2Darrh5(outfile,'G1',G1)
  call writereal2Darrh5(outfile,'G2',G2)
  call writereal2Darrh5(outfile,'G3',G3)
  call writereal2Darrh5(outfile,'Cm1',Cm1)
  call writereal2Darrh5(outfile,'Cm2',Cm2)
  call writereal2Darrh5(outfile,'Cm3',Cm3)
  call writereal2Darrh5(outfile,'Cd1',Cd1)
  call writereal2Darrh5(outfile,'Cd2',Cd2)
  call writereal2Darrh5(outfile,'Cd3',Cd3)
  call writereal1Darrh5(outfile,'mtrue',mtrue)
  call writereal1Darrh5(outfile,'dobs',dobs)
  call writereal1Darrh5(outfile,'mprior',mprior)
  
  !! sizes
  nm1 = size(Cm1,1)
  nm2 = size(Cm2,1)
  nm3 = size(Cm3,1)
  nd1 = size(Cd1,1)
  nd2 = size(Cd2,1)
  nd3 = size(Cd3,1)
  nm = nm1*nm2*nm3
  nd = nd1*nd2*nd3

  print*,"Sizes:"
  print*,"M:",nm1,nm2,nm3,nm
  print*,"D:",nd1,nd2,nd3,nd
  print*
    
  !! allocate working arrays
  print*,'Allocating working arrays'
  allocate(U1(nm1,nm1),U2(nm2,nm2),U3(nm3,nm3))
  allocate(invlambda(nm))
  allocate(iUCm1(nm1,nm1),iUCm2(nm2,nm2),iUCm3(nm3,nm3))
  allocate(iUCmGtiCd1(nm1,nd1),iUCmGtiCd2(nm2,nd2),iUCmGtiCd3(nm3,nd3))
                    
  !! compute the factors for post mean and covariance
  print*,'Computing the factors for post mean and covariance'
  call calcfactors(G1,G2,G3,Cm1,Cm2,Cm3,Cd1,Cd2,Cd3,&
       U1,U2,U3,invlambda,iUCm1,iUCm2,iUCm3,&
       iUCmGtiCd1,iUCmGtiCd2,iUCmGtiCd3 )


  call writereal2Darrh5(outfile,'U1',U1)
  call writereal2Darrh5(outfile,'U2',U2)
  call writereal2Darrh5(outfile,'U3',U3)
  call writereal1Darrh5(outfile,'invlambda',invlambda)
  call writereal2Darrh5(outfile,'iUCm1',iUCm1)
  call writereal2Darrh5(outfile,'iUCm2',iUCm2)
  call writereal2Darrh5(outfile,'iUCm3',iUCm3)
  call writereal2Darrh5(outfile,'iUCmGtiCd1',iUCmGtiCd1)
  call writereal2Darrh5(outfile,'iUCmGtiCd2',iUCmGtiCd2)
  call writereal2Darrh5(outfile,'iUCmGtiCd3',iUCmGtiCd3)
  ! call writereal1Darrh5(outfile,'nm',[nm1,nm2,nm3])
  ! call writereal1Darrh5(outfile,'nd',[nd1,nd2,nd3])

  

  allocate(postm(nm))
  postm=0.0_dp

  call cpu_time(start)
  print*,'Computing post mean'
  call posteriormean(U1,U2,U3, invlambda, iUCmGtiCd1,iUCmGtiCd2,iUCmGtiCd3,&
       G1,G2,G3,mprior,dobs, postm ) 
  call cpu_time(finish)
  print*,'postm time:',finish-start
  !! write posterior mean
  call writereal1Darrh5(outfile,'postm',postm)


  
  ! print*,'Computing postcov'
  ! i1 = 1
  ! i2 = 1000
  ! j1 = 1
  ! j2 = 1000
  ! allocate(postcov(i2-i1+1,j2-j1+1))
  ! postcov = 0.0_dp
  ! call cpu_time(start)
  ! call blockpostcov(U1,U2,U3, invlambda,iUCm1,iUCm2,iUCm3, &
  !                i1,i2,j1,j2, postcov(i1:i2,j1:j2))     
  ! call cpu_time(finish)
  ! print*,'postCov block time:',finish-start
  ! !! write posterior covariance
  ! call writereal2Darrh5(outfile,'postcov',postcov)

  
  ! print*,
  ! print*,'computing postcov_diag'
  ! !!http://www.netlib.org/lapack/lug/node124.html
  ! lowdiag = 0
  ! updiag  = 0
  ! allocate(postcov_diag(1,nm))
  ! postcov_diag = 0.0_dp
  ! call cpu_time(start)
  ! call bandpostcov(U1,U2,U3, invlambda,iUCm1,iUCm2,iUCm3, &
  !                lowdiag,updiag, postcov_diag)     
  ! call cpu_time(finish)
  ! print*,'postCov band time:',finish-start
  ! !! write posterior covariance
  ! call writereal2Darrh5(outfile,'postcov_diag',postcov_diag)


  ! print*,'Computing band of post cov'
  ! lowdiag = 10
  ! updiag  = 10
  ! !! allocate the array for band storage
  ! allocate(bandpostC(lowdiag+updiag+1,nm))
  ! call bandpostcov(U1,U2,U3, invlambda,iUCm1,iUCm2,iUCm3, &
  !                 lowdiag, updiag, bandpostC) 
  ! !! write posterior covariance
  ! call writereal2Darrh5(outfile,'bandpostcov',bandpostC)

  

end program test
