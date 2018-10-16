!-------------------------------------------
!
!> @brief Module to read and write HDF5 files.
!
!> @author
!> Andrea Zunino \n
!> Niels Bohr Institute, University of Copenhagen \n
!> zunino@nbi.dk
!
!-------------------------------------------
module readwriteh5

  ! HDF5:
  ! When a C application reads data stored from a Fortran program,
  ! the data will appear to be transposed due to the difference in
  ! the C and Fortran storage orders. For example, if Fortran writes
  ! a 4x6 two-dimensional dataset to the file, a C program will read
  ! it as a 6x4 two-dimensional dataset into memory. The HDF5 C
  ! utilities h5dump and h5ls will also display transposed data, if
  ! data is written from a Fortran program.
  
  use realprec
  use hdf5
  implicit none
  
  
  logical,private :: h5firsttimeread = .true.
  
contains
  
  !!=============================================================
  
  subroutine writereal1Darrh5(outfile,dsetname,arr)
     
    character(len=1024),intent(in) :: outfile
    real(dp),intent(in) :: arr(:)
    character(len=*),intent(in) :: dsetname
    
    ! filename length must be the same than in dummy variable 
    integer(hid_t)     :: file_id      ! file identifier
    integer(hid_t)     :: dset_id      ! dataset identifier
    integer            :: error        ! error flag
    integer(hid_t)     :: dspace_id
    integer(hsize_t), dimension(1) :: data_dims1
    integer :: rank
       
    !! copmpression stuff
    logical :: avail
    integer ::  filter_info
    integer :: filter_info_both
    integer(hid_t)  ::  dcpl
    integer  :: chunkingfactor
    integer(hsize_t), dimension(1) :: chunk1

    
    chunkingfactor = 10

    !##############################
    !#    initialize hdf5         #
    !##############################
    ! initialize fortran interface.
    call h5open_f(error)
      !-------------------------------------------------------------
    !  check if gzip compression is available and can be used for both
    !  compression and decompression.  normally we do not perform error
    !  checking in these examples for the sake of clarity, but in this
    !  case we will make an exception because this filter is an
    !  optional part of the hdf5 library.
    call h5zfilter_avail_f(h5z_filter_deflate_f, avail, error)
    if (.not.avail) then
       write(*,'("gzip filter not available.",/)')
       stop
    endif
    call h5zget_filter_info_f(h5z_filter_deflate_f, filter_info, error)
    filter_info_both=ior(h5z_filter_encode_enabled_f,h5z_filter_decode_enabled_f)
    if (filter_info .ne. filter_info_both) then
       write(*,'("gzip filter not available for encoding and decoding.",/)')
       stop
    endif

    !-------------------------------------------------------------
    if (h5firsttimeread.eqv..true.) then
       call h5fcreate_f(trim(adjustl(outfile)), h5f_acc_trunc_f,file_id,error)
       h5firsttimeread=.false.
    else
       ! open an existing file using the default properties.
       call h5fopen_f(trim(adjustl(outfile)), h5f_acc_rdwr_f, file_id, error)
    end if

    
    !##############################
    !#  write arr                 #
    !##############################
    rank=1
    data_dims1=shape(arr)
    chunk1 = data_dims1/chunkingfactor
    call checkchunking(chunk1)
    ! create dataspace.  setting maximum size to null sets the maximum
    ! size to be the current size.
    call h5screate_simple_f (rank, data_dims1, dspace_id, error)
    !-------------------------------------------------------------
    ! create the dataset creation property list, add the gzip
    ! compression filter and set the chunk size.
    call h5pcreate_f(h5p_dataset_create_f, dcpl, error)
    call h5pset_deflate_f(dcpl, 9, error)
    call h5pset_chunk_f(dcpl, rank, chunk1, error)
    !-------------------------------------------------------------
    ! create the dataset with default properties.
    call h5dcreate_f(file_id, dsetname, H5T_IEEE_F64LE, dspace_id, &
         dset_id, error, dcpl)
    !-------------------------------------------------------------
    ! write the data to the dataset.
    call h5dwrite_f(dset_id,H5T_IEEE_F64LE, arr, data_dims1, error)
    !-------------------------------------------------------------
    ! close and release resources.
    call h5pclose_f(dcpl,error)
    !-------------------------------------------------------------
    call h5sclose_f(dspace_id, error)
    ! close the dataset.
    call h5dclose_f(dset_id, error)


    ! ##############################
    ! #    terminate hdf5 stuff    #
    ! ##############################
    ! close the file.
    call h5fclose_f(file_id, error)
    ! close fortran interface.
    call h5close_f(error)
    
    return
  end subroutine writereal1Darrh5

 !!=============================================================
  
  subroutine writereal2Darrh5(outfile,dsetname,arr)
     
    character(len=1024),intent(in) :: outfile
    real(dp),intent(in) :: arr(:,:)
    character(len=*),intent(in) :: dsetname
    
    ! filename length must be the same than in dummy variable 
    integer(hid_t)     :: file_id      ! file identifier
    integer(hid_t)     :: dset_id      ! dataset identifier
    integer            :: error        ! error flag
    integer(hid_t)     :: dspace_id
    integer(hsize_t), dimension(2) :: data_dims2
    integer :: rank
       
    !! copmpression stuff
    logical :: avail
    integer ::  filter_info
    integer :: filter_info_both
    integer(hid_t)  ::  dcpl
    integer  :: chunkingfactor
    integer(hsize_t), dimension(2) :: chunk2

    
    chunkingfactor = 10

    !##############################
    !#    initialize hdf5         #
    !##############################
    ! initialize fortran interface.
    call h5open_f(error)
      !-------------------------------------------------------------
    !  check if gzip compression is available and can be used for both
    !  compression and decompression.  normally we do not perform error
    !  checking in these examples for the sake of clarity, but in this
    !  case we will make an exception because this filter is an
    !  optional part of the hdf5 library.
    call h5zfilter_avail_f(h5z_filter_deflate_f, avail, error)
    if (.not.avail) then
       write(*,'("gzip filter not available.",/)')
       stop
    endif
    call h5zget_filter_info_f(h5z_filter_deflate_f, filter_info, error)
    filter_info_both=ior(h5z_filter_encode_enabled_f,h5z_filter_decode_enabled_f)
    if (filter_info .ne. filter_info_both) then
       write(*,'("gzip filter not available for encoding and decoding.",/)')
       stop
    endif

    !-------------------------------------------------------------
    if (h5firsttimeread.eqv..true.) then
       call h5fcreate_f(trim(adjustl(outfile)), h5f_acc_trunc_f,file_id,error)
       h5firsttimeread=.false.
    else
       ! open an existing file using the default properties.
       call h5fopen_f(trim(adjustl(outfile)), h5f_acc_rdwr_f, file_id, error)
    end if

    
    !##############################
    !#  write arr                 #
    !##############################
    rank=2
    data_dims2=shape(arr)
    chunk2 = data_dims2/chunkingfactor
    call checkchunking(chunk2)
    ! create dataspace.  setting maximum size to null sets the maximum
    ! size to be the current size.
    call h5screate_simple_f(rank, data_dims2, dspace_id, error)
    !-------------------------------------------------------------
    ! create the dataset creation property list, add the gzip
    ! compression filter and set the chunk size.
    call h5pcreate_f(h5p_dataset_create_f, dcpl, error)
    call h5pset_deflate_f(dcpl, 9, error)
    call h5pset_chunk_f(dcpl, rank, chunk2, error)
    !-------------------------------------------------------------
    ! create the dataset with default properties.
    call h5dcreate_f(file_id, dsetname, H5T_IEEE_F64LE, dspace_id, &
         dset_id, error, dcpl)
    !-------------------------------------------------------------
    ! write the data to the dataset.
    call h5dwrite_f(dset_id,H5T_IEEE_F64LE, arr, data_dims2, error)
    !-------------------------------------------------------------
    ! close and release resources.
    call h5pclose_f(dcpl,error)
    !-------------------------------------------------------------
    call h5sclose_f(dspace_id, error)
    ! close the dataset.
    call h5dclose_f(dset_id, error)


    ! ##############################
    ! #    terminate hdf5 stuff    #
    ! ##############################
    ! close the file.
    call h5fclose_f(file_id, error)
    ! close fortran interface.
    call h5close_f(error)
    
    return
  end subroutine writereal2Darrh5

  !!============================================================

  subroutine readreal1Darrh5(inpfile,dsetname,arr)

    implicit none
    character(len=1024),intent(in) :: inpfile
    real(dp),allocatable,intent(inout) :: arr(:)
    character(len=*),intent(in) :: dsetname
    
    ! filename length must be the same than in dummy variable 
    integer(hid_t)     :: file_id,dset_id   
    integer            :: rank,error       
    integer(hid_t)     :: dataspace_id
    integer(hsize_t),allocatable :: data_dims(:),maxdims(:)

    ! Initialize FORTRAN interface.
    CALL h5open_f(error)
    ! Open an existing file.
    !print*, 'Reading ',trim(inpfile)
    CALL h5fopen_f(trim(inpfile), H5F_ACC_RDONLY_F, file_id, error)
    if (error/=0) then
       write(*,*)
       write(*,*) "Error opening file from hdf5. "
       write(*,*) "filename: ",trim(inpfile)
       stop
    endif
  
    ! Open an existing dataset.
    call h5dopen_f(file_id, trim(dsetname), dset_id, error)
    if (error/=0) then
       write(*,*)
       write(*,*) "Error opening file from hdf5. "
       write(*,*) "err dsetname: ",trim(dsetname)
       stop
    endif
    ! get the dataspace id
    call h5dget_space_f(dset_id, dataspace_id, error) 
    ! get the rank of the dataset
    call h5sget_simple_extent_ndims_f(dataspace_id, rank, error)
    allocate(data_dims(rank),maxdims(rank))
    ! get the dimensions of the dataset
    call h5sget_simple_extent_dims_f(dataspace_id, data_dims, maxdims, error)
    !print*,'rank:',rank,'maxdims:',maxdims
    ! get the datatype
    !call h5dget_type_f(dset_id, datatype_id, error)
    allocate(arr(data_dims(1)))
    print*,"Reading ",trim(dsetname), " with shape ", data_dims
    ! read data
    CALL h5dread_f(dset_id, H5T_IEEE_F64LE, arr, data_dims, error)
    ! Close the dataset.
    CALL h5dclose_f(dset_id, error)
    ! close the datatype
  
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    ! Close FORTRAN interface.
    CALL h5close_f(error)

    return
  end subroutine readreal1Darrh5

  !!===============================================================
  
  subroutine readreal2Darrh5(inpfile,dsetname,arr)

    implicit none
    character(len=1024),intent(in) :: inpfile
    real(dp),allocatable,intent(inout) :: arr(:,:)
    character(len=*),intent(in) :: dsetname
    
    ! filename length must be the same than in dummy variable 
    integer(hid_t)     :: file_id,dset_id   
    integer            :: rank,error       
    integer(hid_t)     :: dataspace_id
    integer(hsize_t),allocatable :: data_dims(:),maxdims(:)

    ! Initialize FORTRAN interface.
    CALL h5open_f(error)
    ! Open an existing file.
    !print*, 'Reading ',trim(inpfile)
    CALL h5fopen_f(trim(inpfile), H5F_ACC_RDONLY_F, file_id, error)
    if (error/=0) then
       write(*,*)
       write(*,*) "Error opening file from hdf5. "
       write(*,*) "filename: ",trim(inpfile)
       stop
    endif
  
    ! Open an existing dataset.
    call h5dopen_f(file_id, trim(dsetname), dset_id, error)
    if (error/=0) then
       write(*,*)
       write(*,*) "Error opening file from hdf5. "
       write(*,*) "err dsetname: ",trim(dsetname)
       stop
    endif
    ! get the dataspace id
    call h5dget_space_f(dset_id, dataspace_id, error) 
    ! get the rank of the dataset
    call h5sget_simple_extent_ndims_f(dataspace_id, rank, error)
    allocate(data_dims(rank),maxdims(rank))
    ! get the dimensions of the dataset
    call h5sget_simple_extent_dims_f(dataspace_id, data_dims, maxdims, error)
    !print*,'rank:',rank,'maxdims:',maxdims
    ! get the datatype
    !call h5dget_type_f(dset_id, datatype_id, error)
    allocate(arr(data_dims(1),data_dims(2)))
    print*,"Reading ",trim(dsetname), " with shape ", data_dims
    ! read data
    CALL h5dread_f(dset_id, H5T_IEEE_F64LE, arr, data_dims, error)
    ! Close the dataset.
    CALL h5dclose_f(dset_id, error)
    ! close the datatype
  
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    ! Close FORTRAN interface.
    CALL h5close_f(error)

    return
  end subroutine readreal2Darrh5

  !!!===========================================================

  subroutine checkchunking(chunk)

    use hdf5 ! this module contains all necessary modules
    use realprec
    implicit none
   
    integer(hsize_t), dimension(:), intent(inout) :: chunk

    where ( chunk < 1 ) chunk=1
           
    return
  end subroutine checkchunking
  
  !!=============================================================
  
  
end module readwriteh5
