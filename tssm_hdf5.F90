#ifdef _QUADPRECISION_
module tssmq_hdf5
#else
module tssm_hdf5
#endif

#ifndef _NO_HDF5_

#ifdef _QUADPRECISION_
use tssmq
use tssmq_grid
#else
use tssm
use tssm_grid
#endif

use hdf5

implicit none

contains

function hdf5_open_gridfun(filename) result(file_id)
    character(len=*), intent(in) :: filename
    integer(HID_T) :: file_id       ! File identifier 

    integer :: error ! Error flag

#ifdef _MPI_
    if (this_proc/=0) return

    CALL h5open_f(error) 
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
    !immediately close and reopen, because
    !MPI version needs to read from just created file
    CALL h5fclose_f(file_id, error)
    CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
#else
    CALL h5open_f(error) 
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
#endif
end function hdf5_open_gridfun



function hdf5_open_gridfun_read_only(filename) result(file_id)
    character(len=*), intent(in) :: filename
    integer(HID_T) :: file_id       ! File identifier 

    integer :: error ! Error flag

#ifdef _MPI_
    if (this_proc/=0) return
#endif    
    CALL h5open_f(error) 
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
end function hdf5_open_gridfun_read_only



subroutine hdf5_close_gridfun(file_id)
    integer(HID_T), intent(in) :: file_id       ! File identifier 

    integer :: error ! Error flag
#ifdef _MPI_
    if (this_proc/=0) return
#endif    
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
end subroutine hdf5_close_gridfun


subroutine hdf5_write_real_gridfun(file_id, g, u, dset_name)
    integer(HID_T), intent(in) :: file_id       ! File identifier 
    class(grid) :: g
    real(prec), intent(inout), target :: u(*)
    character(len=*), intent(in) :: dset_name
#ifdef _MPI_
    call  hdf5_write_gridfun_MPI_0(file_id, g, u, .true., dset_name)
#else
    call  hdf5_write_gridfun_0(file_id, g, u, .true., dset_name)
#endif    
end subroutine hdf5_write_real_gridfun


subroutine hdf5_write_complex_gridfun(file_id, g, u, dset_name_real, dset_name_imag)
    use, intrinsic :: iso_c_binding
    integer(HID_T), intent(in) :: file_id       ! File identifier 
    class(grid) :: g
    complex(prec), intent(inout), target :: u(*)
    character(len=*), intent(in) :: dset_name_real
    character(len=*), intent(in) :: dset_name_imag
    real(prec), pointer :: ur(:)
    call c_f_pointer(c_loc(u), ur, [1])
#ifdef _MPI_
    call  hdf5_write_gridfun_MPI_0(file_id, g, ur, .false., dset_name_real, dset_name_imag)
#else
    call  hdf5_write_gridfun_0(file_id, g, ur, .false., dset_name_real, dset_name_imag)
#endif    
end subroutine hdf5_write_complex_gridfun



subroutine hdf5_save_real_gridfun(g, u, filename, dset_name)
    class(grid) :: g
    real(prec), intent(inout), target :: u(*)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dset_name

    integer(HID_T) :: file_id       ! File identifier 

    file_id = hdf5_open_gridfun(filename)
    call hdf5_write_real_gridfun(file_id, g, u, dset_name)
    call hdf5_close_gridfun(file_id)
end subroutine hdf5_save_real_gridfun

subroutine hdf5_save_complex_gridfun(g, u, filename, dset_name_real, dset_name_imag)
    use, intrinsic :: iso_c_binding
    class(grid) :: g
    complex(prec), intent(inout), target :: u(*)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dset_name_real
    character(len=*), intent(in) :: dset_name_imag

    integer(HID_T) :: file_id       ! File identifier 

    file_id = hdf5_open_gridfun(filename)
    call  hdf5_write_complex_gridfun(file_id, g, u, dset_name_real, dset_name_imag)
    call hdf5_close_gridfun(file_id)
end subroutine hdf5_save_complex_gridfun


subroutine hdf5_read_real_gridfun(file_id, g, u, dset_name, offset)
    integer(HID_T), intent(in) :: file_id       ! File identifier 
    class(grid) :: g
    real(prec), intent(inout), target :: u(*)
    character(len=*), intent(in) :: dset_name
    integer, intent(in), optional :: offset(:)
#ifdef _MPI_
    call  hdf5_read_gridfun_MPI_0(file_id, g, u, .true., dset_name, offset=offset)
#else
    call  hdf5_read_gridfun_0(file_id, g, u, .true., dset_name, offset=offset)
#endif    
end subroutine hdf5_read_real_gridfun


subroutine hdf5_read_complex_gridfun(file_id, g, u, dset_name_real, dset_name_imag, offset)
    use, intrinsic :: iso_c_binding
    integer(HID_T), intent(in) :: file_id       ! File identifier 
    class(grid) :: g
    complex(prec), intent(inout), target :: u(*)
    character(len=*), intent(in) :: dset_name_real
    character(len=*), intent(in) :: dset_name_imag
    integer, intent(in), optional :: offset(:)
    real(prec), pointer :: ur(:)
    call c_f_pointer(c_loc(u), ur, [1])
#ifdef _MPI_
    call  hdf5_read_gridfun_MPI_0(file_id, g, ur, .false., dset_name_real, dset_name_imag, offset)
#else
    call  hdf5_read_gridfun_0(file_id, g, ur, .false., dset_name_real, dset_name_imag, offset)
#endif    
end subroutine hdf5_read_complex_gridfun




subroutine hdf5_load_real_gridfun(g, u, filename, dset_name, offset)
    class(grid) :: g
    real(prec), intent(inout), target :: u(*)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dset_name
    integer, intent(in), optional :: offset(:)

    integer(HID_T) :: file_id       ! File identifier 

    file_id = hdf5_open_gridfun_read_only(filename)
    call hdf5_read_real_gridfun(file_id, g, u, dset_name, offset=offset)
    call hdf5_close_gridfun(file_id)
end subroutine hdf5_load_real_gridfun



subroutine hdf5_load_complex_gridfun(g, u, filename, dset_name_real, dset_name_imag, offset)
    use, intrinsic :: iso_c_binding
    class(grid) :: g
    complex(prec), intent(inout), target :: u(*)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dset_name_real
    character(len=*), intent(in), optional :: dset_name_imag
    integer, intent(in), optional :: offset(:)

    integer(HID_T) :: file_id       ! File identifier 

    file_id = hdf5_open_gridfun_read_only(filename)
    call  hdf5_read_complex_gridfun(file_id, g, u, dset_name_real, dset_name_imag, offset)
    call hdf5_close_gridfun(file_id)
end subroutine hdf5_load_complex_gridfun



subroutine hdf5_write_gridfun_0(file_id, g, u, real_wavefunction, &
                               dset_name_real, dset_name_imag)
    use, intrinsic :: iso_c_binding
    integer(HID_T), intent(in) :: file_id       ! File identifier 
    class(grid) :: g
    real(kind=prec), intent(inout), target :: u(*)
    logical, intent(in) :: real_wavefunction 
    character(len=*), intent(in) :: dset_name_real
    character(len=*), intent(in), optional :: dset_name_imag
    
    integer :: dimensions

    logical :: file_exists
    integer :: error ! Error flag
    integer :: dset_exists ! Error flag
    integer :: k

    integer(HID_T) :: dset_id       ! Dataset identifier 
    integer(HID_T) :: filespace     ! Dataspace identifier in file 
    integer(HID_T) :: memspace      ! Dataspace identifier in memory

    integer(HSIZE_T) :: dims(3) 
    integer(HSIZE_T) :: dims_file(3) 
    integer(HSIZE_T) :: count(3)
    integer(HSIZE_T) :: offset_file(3)
    integer(HSIZE_T) :: offset_real(3) = (/ 0, 0, 0 /)
    integer(HSIZE_T) :: offset_imag(3) = (/ 1, 0, 0 /)
    integer(HSIZE_T) :: stride(3) = (/ 1, 1, 1 /)
#ifdef _QUADPRECISION_    
    integer(hid_t) :: QUAD_TYPE   
    real(8), pointer :: ur(:)

    call c_f_pointer(c_loc(u), ur, [1])

    !call H5Tcopy_f(H5T_NATIVE_DOUBLE, QUAD_TYPE, error)
    call H5Tcopy_f(H5T_IEEE_F64LE, QUAD_TYPE, error)
    call h5Tset_size_f(QUAD_TYPE, 16_8, error)
    call h5Tset_precision_f(QUAD_TYPE, 128_8, error)
    call H5Tset_fields_f(QUAD_TYPE, 127_8, 112_8, 15_8, 0_8, 112_8, error)
    call H5Tset_ebias_f(QUAD_TYPE, 16383_8, error)
#endif

    if (real_wavefunction) then
        stride(1) = 1 
    else
        stride(1) = 2
    end if

    select type (g)
    class is (grid_1D)
        dimensions = 1
        dims_file(1) = g%nn1max-g%nn1min+1
        offset_file(1) = g%n1min-g%nn1min
        dims(1) = g%m1max-g%m1min+1
        count(1) = g%n1max-g%n1min+1
    class is (grid_2D)
        dimensions = 2
        dims_file(1:2) = (/ g%nn1max-g%nn1min+1,  g%nn2max-g%nn2min+1 /)
        offset_file(1:2) = (/ g%n1min-g%nn1min, g%n2min-g%nn2min /)
        dims(1:2) = (/ g%m1max-g%m1min+1, g%m2max-g%m2min+1 /)
        count(1:2) = (/ g%n1max-g%n1min+1, g%n2max-g%n2min+1 /)
    class is (grid_3D)
        dimensions = 3
        dims_file(1:3) = (/ g%nn1max-g%nn1min+1, g%nn2max-g%nn2min+1, & 
                            g%nn3max-g%nn3min+1 /)
        offset_file(1:3) = (/ g%n1min-g%nn1min, g%n2min-g%nn2min, g%n3min-g%nn3min /)
        dims(1:3) = (/ g%m1max-g%m1min+1, g%m2max-g%m2min+1, &
                            g%m3max-g%m3min+1 /)
        count(1:3) = (/ g%n1max-g%n1min+1, g%n2max-g%n2min+1, &
                            g%n3max-g%n3min+1 /)
    class default    
        stop "XXXX default"
    end select 

    if (.not.real_wavefunction) then
        dims(1) = 2*dims(1)
    end if


    !!!!!!!!!!! real/imag part !!!!!!!!!!!!!!!!!!!
    do k=1,2
          if ((real_wavefunction).and.k>=2) exit
          if ((.not.present(dset_name_imag)).and.k>=2) exit
          !
          ! Create the data space for the  dataset. 
          !
          CALL h5screate_simple_f(dimensions, dims_file, filespace, error)
     
          !
          ! Create the dataset with default properties.
          !
          CALL h5eset_auto_f(0, error)

          if (k==1) then
#ifdef _QUADPRECISION_    
               CALL h5dcreate_f(file_id, dset_name_real,  QUAD_TYPE, filespace, &               
                      dset_id, dset_exists)
#else
               CALL h5dcreate_f(file_id, dset_name_real,  H5T_IEEE_F64LE, filespace, &               
                      dset_id, dset_exists)
#endif                      
               CALL h5eset_auto_f(1, error)
               if (dset_exists<0) then
                      CALL h5dopen_f(file_id, dset_name_real, dset_id, error)
               end if      
          else
#ifdef _QUADPRECISION_    
               CALL h5dcreate_f(file_id, dset_name_imag,  QUAD_TYPE, filespace, &
                           dset_id, dset_exists)
#else
               CALL h5dcreate_f(file_id, dset_name_imag,  H5T_IEEE_F64LE, filespace, &
                           dset_id, dset_exists)
#endif                           
               CALL h5eset_auto_f(1, error)
               if (dset_exists<0) then
                      CALL h5eset_auto_f(1, error)
                      CALL h5dopen_f(file_id, dset_name_imag, dset_id, error)
               end if      
          end if
          
          !
          ! Each process defines dataset in memory and writes it to the hyperslab
          ! in the file. 
          !
          CALL h5screate_simple_f(dimensions, dims, memspace, error) 
          ! 
          ! Select hyperslab in memory.
          !
          if (k==1) then
               CALL h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, offset_real, count, error, stride)
          else
               CALL h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, offset_imag, count, error, stride)
          end if
   
          ! 
          ! Select hyperslab in the file.
          !
          CALL h5dget_space_f(dset_id, filespace, error)
          CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset_file, count, error)
          !
          ! Write the dataset.
!print *, "DIMS", dims
!print *, "DIMS_FILE", dims_file
!print *, "COUND", count
!print *, "STRIDE", stride
!print *, "OFFSET_file", offset_file
!print *, "OFFSET_real", offset_real
!print *, "OFFSET_imag", offset_imag
#ifdef _QUADPRECISION_    
          CALL h5dwrite_f(dset_id, QUAD_TYPE, ur, dims, error, &
                      file_space_id = filespace, mem_space_id = memspace)
#else
          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, u, dims, error, &
                      file_space_id = filespace, mem_space_id = memspace)
#endif                          
          !
          ! Close dataspaces.
          !
          CALL h5sclose_f(memspace, error)
  
          !
          ! Close dataset.
          !
          CALL h5dclose_f(dset_id, error)
          CALL h5sclose_f(filespace, error)
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _QUADPRECISION_    
    CALL h5tclose_f(QUAD_TYPE, error)
#endif     
end subroutine hdf5_write_gridfun_0





#ifdef _MPI_
subroutine hdf5_write_gridfun_MPI_0(file_id, g, u, real_wavefunction, &
                               dset_name_real, dset_name_imag)
    integer(HID_T), intent(in) :: file_id       ! File identifier 
    class(grid) :: g
    real(kind=prec), intent(inout), target :: u(*)
    logical, intent(in) :: real_wavefunction 
    character(len=*), intent(in) :: dset_name_real
    character(len=*), intent(in), optional :: dset_name_imag

    integer :: proc, dims, error
    integer :: mpi_stat(MPI_STATUS_SIZE)
    integer :: buf(0:12)
 
    type(grid_1D) :: g1
    type(grid_2D) :: g2
    type(grid_3D) :: g3   

    select type(g)
    class is (grid_1D)
        dims = 1
    class is (grid_2D)
        dims = 2
    class is (grid_3D)
        dims = 3
    end select

     do  proc=0,n_proc-1
         if (proc==this_proc) then

             buf(0) = g%alloc_size
             if (.not.real_wavefunction) then
                 buf(0) = 2*buf(0)
             end if    

             select type(g)
             class is (grid_1D)
                 dims = 1
                 buf(1) = g%m1min
                 buf(2) = g%m1max
                 buf(3) = g%n1min
                 buf(4) = g%n1max
             class is (grid_2D)
                 dims = 2
                 buf(1) = g%m1min
                 buf(2) = g%m1max
                 buf(3) = g%n1min
                 buf(4) = g%n1max
                 buf(5) = g%m2min
                 buf(6) = g%m2max
                 buf(7) = g%n2min
                 buf(8) = g%n2max
             class is (grid_3D)
                 dims = 3
                 buf(1) = g%m1min
                 buf(2) = g%m1max
                 buf(3) = g%n1min
                 buf(4) = g%n1max
                 buf(5) = g%m2min
                 buf(6) = g%m2max
                 buf(7) = g%n2min
                 buf(8) = g%n2max
                 buf(9) = g%m3min
                 buf(10) = g%m3max
                 buf(11) = g%n3min
                 buf(12) = g%n3max
             end select
        end if     
        
        if (this_proc==0) then
             if (proc/=0) then
                 call MPI_Recv(buf, dims*4+1, MPI_INTEGER, proc, 40, MPI_COMM_WORLD, mpi_stat, error)   
                 call MPI_Recv(u, buf(0), MPI_DOUBLE_PRECISION, proc, 41, MPI_COMM_WORLD, mpi_stat, error)
             end if
             select type(g)
             class is (grid_1D)
                 g1%nn1min = g%nn1min
                 g1%nn1max = g%nn1max
                 g1%m1min = buf(1)
                 g1%m1max = buf(2)
                 g1%n1min = buf(3)
                 g1%n1max = buf(4)
                 call hdf5_write_gridfun_0(file_id, g1, u, real_wavefunction, &
                                           dset_name_real, dset_name_imag)
             class is (grid_2D)
                 g2%nn1min = g%nn1min
                 g2%nn1max = g%nn1max
                 g2%nn2min = g%nn2min
                 g2%nn2max = g%nn2max
                 g2%m1min = buf(1)
                 g2%m1max = buf(2)
                 g2%n1min = buf(3)
                 g2%n1max = buf(4)
                 g2%m2min = buf(5)
                 g2%m2max = buf(6)
                 g2%n2min = buf(7)
                 g2%n2max = buf(8)
                 call hdf5_write_gridfun_0(file_id, g2, u, real_wavefunction, &
                                           dset_name_real, dset_name_imag)
             class is (grid_3D)
                 g3%nn1min = g%nn1min
                 g3%nn1max = g%nn1max
                 g3%nn2min = g%nn2min
                 g3%nn2max = g%nn2max
                 g3%nn3min = g%nn3min
                 g3%nn3max = g%nn3max
                 g3%m1min = buf(1)
                 g3%m1max = buf(2)
                 g3%n1min = buf(3)
                 g3%n1max = buf(4)
                 g3%m2min = buf(5)
                 g3%m2max = buf(6)
                 g3%n2min = buf(7)
                 g3%n2max = buf(8)
                 g3%m3min = buf(9)
                 g3%m3max = buf(10)
                 g3%n3min = buf(11)
                 g3%n3max = buf(12)
                 call hdf5_write_gridfun_0(file_id, g3, u, real_wavefunction, &
                                           dset_name_real, dset_name_imag)
             end select
         else if (proc==this_proc) then
             call MPI_Send(buf, dims*4+1, MPI_INTEGER, 0, 40, MPI_COMM_WORLD, error)
             call MPI_Send(u, buf(0), MPI_DOUBLE_PRECISION, 0, 41, MPI_COMM_WORLD, error)
         end if
     end do

     if (this_proc==0) then  !reread data of processor 0, because these were overwritten...
         call hdf5_read_gridfun_0(file_id, g, u, real_wavefunction, &
                               dset_name_real, dset_name_imag)
     end if                               
end subroutine hdf5_write_gridfun_MPI_0
#endif



subroutine hdf5_write_attributes(filename, snames, svalues, inames, ivalues, dnames, dvalues)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: snames(:)
    character(len=*), intent(in), optional :: svalues(:)
    character(len=*), intent(in), optional :: inames(:)
    integer, intent(in), optional :: ivalues(:)
    character(len=*), intent(in), optional :: dnames(:)
    real(prec), intent(in), optional :: dvalues(:)
    integer :: j

    integer(HID_T) :: file_id       ! File identifier 
    integer(HID_T) :: ir_attr       
    integer(HID_T) :: attr
    integer(HID_T) :: stringtype_mem, stringtype_file 
    integer(HSIZE_T) :: dummy(7)
    integer(SIZE_T) :: str_len 
    integer :: error ! Error flag
    integer :: not_exists ! Error flag
#ifdef _QUADPRECISION_
    real(kind=8) :: v
#endif     
#ifdef _MPI_
    if (this_proc/=0) return
#endif    

    CALL h5open_f(error) 
      CALL h5eset_auto_f(0, error)
        CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, not_exists)
      CALL h5eset_auto_f(1, error)
      if (not_exists/=0) then
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      end if
      CALL h5screate_f(H5S_SCALAR_F, ir_attr, error)  

        if (present(snames).and.present(svalues)) then
            CALL H5Tcopy_f(H5T_C_S1, stringtype_file, error)
            CALL H5Tcopy_f( H5T_FORTRAN_S1, stringtype_mem, error)
            do j=1, min(size(snames),size(svalues))
                str_len = len(svalues(j))
                CALL H5Tset_size_f(stringtype_file, str_len+1_size_t, error)
                CALL H5Tset_size_f(stringtype_mem, str_len, error)
                call h5acreate_f (file_id, snames(j), stringtype_file, ir_attr, attr, error)
                call h5awrite_f (attr, stringtype_mem, svalues(j), dummy, error)
                call h5aclose_f(attr, error)
            end do
            CALL h5tclose_f(stringtype_file, error)
            CALL h5tclose_f(stringtype_mem, error)
        end if

        if (present(inames).and.present(ivalues)) then
            do j=1, min(size(inames),size(ivalues))
                call h5acreate_f (file_id, inames(j), H5T_STD_I32LE, ir_attr, attr, error)
                call h5awrite_f (attr, H5T_NATIVE_INTEGER, ivalues(j), dummy, error)
                call h5aclose_f(attr, error)
            end do
        end if

        if (present(dnames).and.present(dvalues)) then
            do j=1, min(size(dnames),size(dvalues))
                call h5acreate_f (file_id, dnames(j), H5T_IEEE_F64LE, ir_attr, attr, error)
#ifdef _QUADPRECISION_
                v = real(dvalues(j), kind=8)
                call h5awrite_f (attr, H5T_NATIVE_DOUBLE, v, dummy, error)
#else
                call h5awrite_f (attr, H5T_NATIVE_DOUBLE, dvalues(j), dummy, error)
#endif                 
                call h5aclose_f(attr, error)
            end do
        end if

      CALL h5sclose_f(ir_attr, error)
      CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
end subroutine hdf5_write_attributes



subroutine hdf5_read_string_attribute(filename, name, value, len)
    USE ISO_C_BINDING
    character(len=*), intent(out), target :: value 
    integer, intent(out) :: len
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: name

    integer(HID_T) :: file_id       ! File identifier 
    integer(HID_T) :: filetype, memtype
    integer(HID_T) :: attr
    integer(HSIZE_T) :: dummy(7)
    INTEGER(SIZE_T) :: string_size
    TYPE(C_PTR) :: f_ptr
    integer :: error ! Error flag

#ifdef _MPI_
    if (this_proc/=0) return
#endif    

    CALL h5open_f(error) 
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
    call h5aopen_f (file_id, name, attr, error)
    !
    ! Get the datatype and its size.
    !
    CALL H5Aget_type_f(attr, filetype, error)
    CALL H5Tget_size_f(filetype, string_size, error)
    !
    ! Create the memory datatype.
    !
    CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, error)
    CALL H5Tset_size_f(memtype, string_size, error)
    !
    ! Read the data.
    !
    !f_ptr = C_LOC(value)
    !CALL H5Aread_f(attr, memtype, f_ptr, dummy, error)
    CALL H5Aread_f(attr, memtype, value, dummy, error)
    
    CALL H5tclose_f(filetype, error)
    CALL H5tclose_f(memtype, error)
    call h5aclose_f(attr, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
    len = string_size
end subroutine hdf5_read_string_attribute



function hdf5_read_integer_attribute(filename, name) result(value)
    integer :: value
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: name

    integer(HID_T) :: file_id       ! File identifier 
    integer(HID_T) :: attr
    integer(HSIZE_T) :: dummy(7)
    integer :: error ! Error flag

    CALL h5open_f(error) 
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
    call h5aopen_f (file_id, name, attr, error)
    call h5aread_f (attr, H5T_NATIVE_INTEGER, value, dummy, error)
    call h5aclose_f(attr, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
end function hdf5_read_integer_attribute



function hdf5_read_double_attribute(filename, name) result(value)
    real(prec) :: value
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: name

    integer(HID_T) :: file_id       ! File identifier 
    integer(HID_T) :: attr
    integer(HSIZE_T) :: dummy(7)
    integer :: error ! Error flag
    
#ifdef _QUADPRECISION_
    real(kind=8) :: value1
#endif     
#ifdef _MPI_
    if (this_proc/=0) return
#endif    

    CALL h5open_f(error) 
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
    call h5aopen_f (file_id, name, attr, error)
#ifdef _QUADPRECISION_
    call h5aread_f (attr, H5T_NATIVE_DOUBLE, value1, dummy, error)
    value = real(value1, kind=prec)
#else
    call h5aread_f (attr, H5T_NATIVE_DOUBLE, value, dummy, error)
#endif     
    call h5aclose_f(attr, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
end function hdf5_read_double_attribute



subroutine hdf5_get_dimensions_0(file_id,  dset_name, rank, dims)
    integer(HID_T), intent(in) :: file_id       ! File identifier 
    character(len=*), intent(in) :: dset_name
    integer, intent(out) :: rank
    integer(HSIZE_T), intent(out) :: dims(*)

    integer(HID_T) :: space_id     ! Dataspace identifier in file 
    integer(HID_T) :: dset_id       ! Dataset identifier 
    
    integer(HID_T) :: attr
    integer(HSIZE_T) :: n(3), max_n(3)
    integer :: k

    integer :: error, not_exists ! Error flag
#ifdef _MPI_
    if (this_proc/=0) return
#endif    

    !CALL h5open_f(error) 
    !CALL  h5eset_auto_f(0, error)
    !call  h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, not_exists)
    !CALL  h5eset_auto_f(1, error)
    !      if (not_exists/=0) then
    !         rank = -1 
    !         CALL h5close_f(error)
    !         return
    !      end if    
    CALL   h5eset_auto_f(0, error)
    call   h5dopen_f(file_id, dset_name, dset_id, not_exists)
           if (not_exists/=0) then
              rank = -2 
              CALL h5close_f(error)
              return
           end if    
    CALL   h5eset_auto_f(1, error)
    call    h5screate_f(H5S_SIMPLE_F, space_id, error)
    call      h5dget_space_f(dset_id, space_id, error)
    call      h5sget_simple_extent_ndims_f(space_id, rank, error) 
    call      h5sget_simple_extent_dims_f(space_id, n, max_n, error)
              do k=1,rank
                dims(k) = n(k)
              end do  
              dims(1) = dims(1)
    call    h5sclose_f(space_id, error)
    CALL   h5dclose_f(dset_id, error)
    !CALL  h5fclose_f(file_id, error)
    !CALL h5close_f(error)
end subroutine hdf5_get_dimensions_0


subroutine hdf5_get_dimensions(filename,  dset_name, rank, dims)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dset_name
    integer, intent(out) :: rank
    integer(HSIZE_T), intent(out) :: dims(*)

    integer(HID_T) :: file_id       ! File identifier 

#ifdef _MPI_
    if (this_proc/=0) return
#endif    
    file_id =  hdf5_open_gridfun_read_only(filename)
    call hdf5_get_dimensions_0(file_id, dset_name, rank, dims)
    call hdf5_close_gridfun(file_id)
end subroutine hdf5_get_dimensions


subroutine hdf5_write_grid_attributes(g, filename)
    class(grid) :: g
    character(len=*), intent(in) :: filename

#ifdef _MPI_
    if (this_proc/=0) return
#endif    
    select type (g)
    class is (grid_equidistant_1D)
        call hdf5_write_attributes(filename, (/ "grid" /), &
          (/ "equidistant_1D" /), &
          (/ "nx ", "ix0" /), (/ g%nx, g%nn1min /), &
          (/ "xmin", "xmax" /), (/ g%xmin, g%xmax /) ) 
    class is (grid_equidistant_2D)
        call hdf5_write_attributes(filename, (/ "grid" /), &
          (/ "equidistant_2D" /), &
          (/ "nx ", "ix0", "ny ", "iy0" /), & 
          (/ g%nx, g%nn1min, g%ny, g%nn2min /), &
          (/ "xmin", "xmax", "ymin", "ymax" /), &
          (/ g%xmin, g%xmax, g%ymin, g%ymax /) ) 
    class is (grid_equidistant_3D)
        call hdf5_write_attributes(filename, (/ "grid" /), &
          (/ "equidistant_3D" /), &
          (/ "nx ", "ix0", "ny ", "iy0", "nz ", "iz0" /), & 
          (/ g%nx, g%nn1min, g%ny, g%nn2min, g%nz, g%nn3min /), &
          (/ "xmin", "xmax", "ymin", "ymax", "zmin", "zmax" /), &
          (/ g%xmin, g%xmax, g%ymin, g%ymax, g%zmin, g%zmax /) ) 
    class default    
        print *, "W: hdf5_write_grid_attributes: grid type not implemented."
    end select 
end subroutine hdf5_write_grid_attributes

subroutine hdf5_read_grid_attributes(g, filename, dset_name)
    class(grid), intent(out) :: g
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dset_name
    integer :: rank
    integer(HSIZE_T) :: dims(3) 
    character(len=20) :: gridname
    integer :: len
#ifdef _MPI_
    integer :: ierr
    real(kind=8) :: rbuf(9)
    integer :: ibuf(9)
#endif
    
    ! TODO error handling !!!

#ifdef _MPI_
    if (this_proc==0) then 
#endif
    call hdf5_get_dimensions(filename, dset_name, rank, dims)
    call hdf5_read_string_attribute(filename, "grid", gridname, len)
#ifdef _MPI_
    end if
#endif     

    select type (g)
    class is (grid_equidistant_1D)
#ifdef _MPI_
        if (this_proc==0) then 
#endif        
            if (gridname(1:len)/="equidistant_1D") then
                stop "E: wrong grid type"
            end if
            g%xmin = hdf5_read_double_attribute(filename, "xmin")
            g%xmax = hdf5_read_double_attribute(filename, "xmax")
            g%nx = hdf5_read_integer_attribute(filename, "nx")
            g%nn1min = hdf5_read_integer_attribute(filename, "ix0")
            g%nn1max = g%nn1min + dims(1) - 1
            g%dx = (g%xmax - g%xmin)/g%nx
#ifdef _MPI_
            rbuf(1) = g%xmin
            rbuf(2) = g%xmax
            rbuf(3) = g%dx
            ibuf(1) = g%nx
            ibuf(2) = g%nn1min
            ibuf(3) = g%nn1max
        end if    
        call MPI_Bcast(rbuf, 3, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(ibuf, 3, MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
        if (this_proc/=0) then 
            g%xmin = rbuf(1)
            g%xmax = rbuf(2)
            g%dx = rbuf(3)
            g%nx = ibuf(1)
            g%nn1min = ibuf(2)
            g%nn1max = ibuf(3)
        end if
#endif 
    class is (grid_equidistant_2D)
#ifdef _MPI_
        if (this_proc==0) then 
#endif        
            if (gridname(1:len)/="equidistant_2D") then
                stop "E: wrong grid type"
            end if
            g%xmin = hdf5_read_double_attribute(filename, "xmin")
            g%xmax = hdf5_read_double_attribute(filename, "xmax")
            g%nx = hdf5_read_integer_attribute(filename, "nx")
            g%nn1min = hdf5_read_integer_attribute(filename, "ix0")
            g%nn1max = g%nn1min + dims(1) - 1
            g%ymin = hdf5_read_double_attribute(filename, "ymin")
            g%ymax = hdf5_read_double_attribute(filename, "ymax")
            g%ny = hdf5_read_integer_attribute(filename, "ny")
            g%nn2min = hdf5_read_integer_attribute(filename, "iy0")
            g%nn2max = g%nn2min + dims(2) - 1
            g%dx = (g%xmax - g%xmin)/g%nx
            g%dy = (g%ymax - g%ymin)/g%ny
#ifdef _MPI_
            rbuf(1) = g%xmin
            rbuf(2) = g%xmax
            rbuf(3) = g%dx
            rbuf(4) = g%ymin
            rbuf(5) = g%ymax
            rbuf(6) = g%dy
            ibuf(1) = g%nx
            ibuf(2) = g%nn1min
            ibuf(3) = g%nn1max
            ibuf(4) = g%ny
            ibuf(5) = g%nn2min
            ibuf(6) = g%nn2max
        end if    
        call MPI_Bcast(rbuf, 6, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(ibuf, 6, MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
        if (this_proc/=0) then 
            g%xmin = rbuf(1)
            g%xmax = rbuf(2)
            g%dx = rbuf(3)
            g%ymin = rbuf(4)
            g%ymax = rbuf(5)
            g%dy = rbuf(6)
            g%nx = ibuf(1)
            g%nn1min = ibuf(2)
            g%nn1max = ibuf(3)
            g%ny = ibuf(4)
            g%nn2min = ibuf(5)
            g%nn2max = ibuf(6)
        end if
#endif 
    class is (grid_equidistant_3D)
#ifdef _MPI_
        if (this_proc==0) then 
#endif        
            if (gridname(1:len)/="equidistant_3D") then
                stop "E: wrong grid type"
            end if
            g%xmin = hdf5_read_double_attribute(filename, "xmin")
            g%xmax = hdf5_read_double_attribute(filename, "xmax")
            g%nx = hdf5_read_integer_attribute(filename, "nx")
            g%nn1min = hdf5_read_integer_attribute(filename, "ix0")
            g%nn1max = g%nn1min + dims(1) - 1
            g%ymin = hdf5_read_double_attribute(filename, "ymin")
            g%ymax = hdf5_read_double_attribute(filename, "ymax")
            g%ny = hdf5_read_integer_attribute(filename, "ny")
            g%nn2min = hdf5_read_integer_attribute(filename, "iy0")
            g%nn2max = g%nn2min + dims(2) - 1
            g%zmin = hdf5_read_double_attribute(filename, "zmin")
            g%zmax = hdf5_read_double_attribute(filename, "zmax")
            g%nz = hdf5_read_integer_attribute(filename, "nz")
            g%nn3min = hdf5_read_integer_attribute(filename, "iz0")
            g%nn3max = g%nn3min + dims(3) - 1
            g%dx = (g%xmax - g%xmin)/g%nx
            g%dy = (g%ymax - g%ymin)/g%ny
            g%dz = (g%zmax - g%zmin)/g%nz
#ifdef _MPI_
            rbuf(1) = g%xmin
            rbuf(2) = g%xmax
            rbuf(3) = g%dx
            rbuf(4) = g%ymin
            rbuf(5) = g%ymax
            rbuf(6) = g%dy
            rbuf(7) = g%zmin
            rbuf(8) = g%zmax
            rbuf(9) = g%dz
            ibuf(1) = g%nx
            ibuf(2) = g%nn1min
            ibuf(3) = g%nn1max
            ibuf(4) = g%ny
            ibuf(5) = g%nn2min
            ibuf(6) = g%nn2max
            ibuf(7) = g%nz
            ibuf(8) = g%nn3min
            ibuf(9) = g%nn3max
        end if    
        call MPI_Bcast(rbuf, 9, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(ibuf, 9, MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
        if (this_proc/=0) then 
            g%xmin = rbuf(1)
            g%xmax = rbuf(2)
            g%dx = rbuf(3)
            g%ymin = rbuf(4)
            g%ymax = rbuf(5)
            g%dy = rbuf(6)
            g%zmin = rbuf(7)
            g%zmax = rbuf(8)
            g%dz = rbuf(9)
            g%nx = ibuf(1)
            g%nn1min = ibuf(2)
            g%nn1max = ibuf(3)
            g%ny = ibuf(4)
            g%nn2min = ibuf(5)
            g%nn2max = ibuf(6)
            g%nz = ibuf(7)
            g%nn3min = ibuf(8)
            g%nn3max = ibuf(9)
        end if
#endif 
        
    class default    
        stop "E: hdf5_read_grid_attributes: grid type not implemented."
    end select 
end subroutine hdf5_read_grid_attributes


!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

subroutine hdf5_read_gridfun_0(file_id, g, u, real_wavefunction, &
                               dset_name_real, dset_name_imag, offset)
    use, intrinsic :: iso_c_binding
    integer(HID_T), intent(in) :: file_id       ! File identifier 
    class(grid) :: g
    real(prec), target, intent(inout) :: u(*)
    logical, intent(in) :: real_wavefunction 
    character(len=*), intent(in) :: dset_name_real
    character(len=*), intent(in), optional :: dset_name_imag
    integer, intent(in), optional :: offset(:)
    
    integer :: dimensions

    logical :: file_exists
    integer :: error ! Error flag
    integer :: dset_exists ! Error flag
    integer :: k

    integer(HID_T) :: dset_id       ! Dataset identifier 
    integer(HID_T) :: filespace     ! Dataspace identifier in file 
    integer(HID_T) :: memspace      ! Dataspace identifier in memory

    integer(HSIZE_T) :: dims_mem(3) 
    integer(HSIZE_T) :: dims_file(3) 
    integer(HSIZE_T) :: count_mem(3)
    integer(HSIZE_T) :: offset0(3)
    integer(HSIZE_T) :: offset_file(3)
    integer(HSIZE_T) :: offset_real(3) 
    integer(HSIZE_T) :: offset_imag(3) 
    integer(HSIZE_T) :: stride(3)
    integer(HSIZE_T) :: mem_size, i
#ifdef _QUADPRECISION_    
    integer(hid_t) :: QUAD_TYPE   
    real(8), pointer :: ur(:)

    call c_f_pointer(c_loc(u), ur, [1])
#endif    


    !call hdf5_get_dimensions(filename, dset_name_real, dimensions, dims_file)
    call hdf5_get_dimensions_0(file_id, dset_name_real, dimensions, dims_file)

    !if (dimensions==-1) then
    !     print *,  "E: file '", filename, "' not found"
    !     stop
    !else 
    if (dimensions==-2) then
         print *,  "E: dataset '", dset_name_real, "' not found"
         stop
    end if

    stride(1:3) = (/ 1, 1, 1 /)
    if (real_wavefunction) then
        stride(1) = 1 
    else
        stride(1) = 2
    end if


    select type (g)
    class is (grid_1D)
        if (dimensions/=1) then
           !print *,  "E: rank of '", dset_name_real, "' in '", filename, "' should be 1'"
           print *,  "E: rank of '", dset_name_real, "' should be 1'"
           stop
        end if   
        dims_mem(1) = g%m1max-g%m1min+1
        offset0(1) = 0
        if (present(offset)) then
            offset0(1) = offset(1) 
        end if
        offset_real(1) = max(0, offset0(1) - (g%m1min-g%nn1min))
        offset_file(1) = max(0, g%n1min-g%nn1min - offset0(1))
        !count_mem(1) = min(dims_file(1), max(0,g%n1max-g%n1min-offset0(1)+1))
        count_mem(1) = min(max(0,dims_file(1)-offset_file(1)), max(0,g%n1max-g%n1min+1-offset_real(1)))
        mem_size = dims_mem(1)
    class is (grid_2D)
        if (dimensions/=2) then
           !print *, "E: rank of '", dset_name_real, "' in '", filename, "' should be 2'"
           print *,  "E: rank of '", dset_name_real, "' should be 2'"
           stop
        end if   
        dims_mem(1:2) = (/ g%m1max-g%m1min+1, g%m2max-g%m2min+1 /)
        offset0(1:2) = (/ 0, 0 /)
        if (present(offset)) then
            offset0(1:2) =offset(1:2)
        end if
        offset_real(1) = max(0, offset0(1) - (g%m1min-g%nn1min))
        offset_real(2) = max(0, offset0(2) - (g%m2min-g%nn2min))
        offset_file(1) = max(0, g%n1min-g%nn1min - offset0(1))
        offset_file(2) = max(0, g%n2min-g%nn2min - offset0(2))
        !count_mem(1) = min(dims_file(1), max(0,g%n1max-g%n1min-offset0(1)+1))
        !count_mem(2) = min(dims_file(2), max(0,g%n2max-g%n2min-offset0(2)+1))
        count_mem(1) = min(max(0,dims_file(1)-offset_file(1)), max(0,g%n1max-g%n1min+1-offset_real(1)))
        count_mem(2) = min(max(0,dims_file(2)-offset_file(2)), max(0,g%n2max-g%n2min+1-offset_real(2)))
        mem_size = dims_mem(1)*dims_mem(2)
    class is (grid_3D)
        if (dimensions/=3) then
           !print *,"E: rank of '", dset_name_real, "' in '", filename, "' should be 3'"
           print *,  "E: rank of '", dset_name_real, "' should be 3'"
           stop
        end if   
        dims_mem(1:3) = (/ g%m1max-g%m1min+1, g%m2max-g%m2min+1, &
                            g%m3max-g%m3min+1 /)
        offset0(1:3) = (/ 0, 0, 0 /)
        if (present(offset)) then
            offset0(1:3) =offset(1:3)
        end if
        offset_real(1) = max(0, offset0(1) - (g%m1min-g%nn1min))
        offset_real(2) = max(0, offset0(2) - (g%m2min-g%nn2min))
        offset_real(3) = max(0, offset0(3) - (g%m3min-g%nn3min))
        offset_file(1) = max(0, g%n1min-g%nn1min - offset0(1))
        offset_file(2) = max(0, g%n2min-g%nn2min - offset0(2))
        offset_file(3) = max(0, g%n3min-g%nn3min - offset0(3))
        !count_mem(1) = min(dims_file(1), max(0,g%n1max-g%n1min-offset0(1)+1))
        !count_mem(2) = min(dims_file(2), max(0,g%n2max-g%n2min-offset0(2)+1))
        !count_mem(3) = min(dims_file(3), max(0,g%n3max-g%n3min-offset0(3)+1))
        count_mem(1) = min(max(0,dims_file(1)-offset_file(1)), max(0,g%n1max-g%n1min+1-offset_real(1)))
        count_mem(2) = min(max(0,dims_file(2)-offset_file(2)), max(0,g%n2max-g%n2min+1-offset_real(2)))
        count_mem(3) = min(max(0,dims_file(3)-offset_file(3)), max(0,g%n3max-g%n3min+1-offset_real(3)))
        mem_size = dims_mem(1)*dims_mem(2)*dims_mem(3)
    class default    
        stop "E: wrong grid"
    end select 

    if (.not.real_wavefunction) then
        mem_size = 2*mem_size
    end if
    do i=1,mem_size
        u(i) = 0.0_prec
    end do

    if (.not.real_wavefunction) then
        dims_mem(1) = 2*dims_mem(1)
        offset_real(1) = 2*offset_real(1)
        offset_imag = offset_real
        offset_imag(1) = offset_imag(1) + 1
    end if

    !
    ! Initialize HDF5
    !
!    CALL h5open_f(error) 
    !
    ! Open file
    !

#ifdef _QUADPRECISION_    
    !call H5Tcopy_f(H5T_NATIVE_DOUBLE, QUAD_TYPE, error)
    call H5Tcopy_f(H5T_IEEE_F64LE, QUAD_TYPE, error)
    call h5Tset_size_f(QUAD_TYPE, 16_8, error)
    call h5Tset_precision_f(QUAD_TYPE, 128_8, error)
    call H5Tset_fields_f(QUAD_TYPE, 127_8, 112_8, 15_8, 0_8, 112_8, error)
    call H5Tset_ebias_f(QUAD_TYPE, 16383_8, error)
#endif

!    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

    !!!!!!!!!!! real/imag part !!!!!!!!!!!!!!!!!!!
    do k=1,2
          if ((real_wavefunction).and.k>=2) exit
          if ((.not.present(dset_name_imag)).and.k>=2) exit
          !
          ! Create the dataset with default properties.
          !
          CALL h5eset_auto_f(0, error)

          if (k==1) then
               CALL h5dopen_f(file_id, dset_name_real, dset_id, dset_exists)
          else
               CALL h5eset_auto_f(0, error)
               CALL h5dopen_f(file_id, dset_name_imag, dset_id, dset_exists)
               CALL h5eset_auto_f(1, error)
               if (dset_exists<0) exit 
          end if
          !
          ! Create the data space for the  dataset. 
          !
          CALL h5screate_simple_f(dimensions, dims_file, filespace, error)
          !
          ! Each process defines dataset in memory and writes it to the hyperslab
          ! in the file. 
          !
          CALL h5screate_simple_f(dimensions, dims_mem, memspace, error) 
          ! 
          ! Select hyperslab in memory.
          !
          if (k==1) then
               CALL h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, offset_real, count_mem, error, stride)
          else
               CALL h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, offset_imag, count_mem, error, stride)
          end if
          ! 
          ! Select hyperslab in the file.
          !
          CALL h5dget_space_f(dset_id, filespace, error)
          CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset_file, count_mem, error)
          !
          ! Read the dataset.
          !
#ifdef _QUADPRECISION_    
          CALL h5dread_f(dset_id, QUAD_TYPE, ur, dims_mem, error, &
                      file_space_id = filespace, mem_space_id = memspace)
#else
          CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, u, dims_mem, error, &
                      file_space_id = filespace, mem_space_id = memspace)
#endif                      

          !
          ! Close dataspaces.
          !
          CALL h5sclose_f(memspace, error)
          CALL h5sclose_f(filespace, error)
   
          !
          ! Close dataset.
          !
          CALL h5dclose_f(dset_id, error)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _QUADPRECISION_    
     CALL h5tclose_f(QUAD_TYPE, error)
#endif     
!    !
!    ! Close the file.
!    !
!    CALL h5fclose_f(file_id, error)
!    !
!    ! Close HDF5
!    !
!    CALL h5close_f(error)

end subroutine hdf5_read_gridfun_0


#ifdef _MPI_
subroutine hdf5_read_gridfun_MPI_0(file_id, g, u, real_wavefunction, &
                               dset_name_real, dset_name_imag, offset)
    integer(HID_T), intent(in) :: file_id       ! File identifier 
    class(grid) :: g
    real(kind=prec), intent(inout), target :: u(*)
    logical, intent(in) :: real_wavefunction 
    character(len=*), intent(in) :: dset_name_real
    character(len=*), intent(in), optional :: dset_name_imag
    integer, intent(in), optional :: offset(:)

    integer :: proc, dims, error
    integer :: mpi_stat(MPI_STATUS_SIZE)
    integer :: buf(0:12)
 
    type(grid_1D) :: g1
    type(grid_2D) :: g2
    type(grid_3D) :: g3

    select type(g)
    class is (grid_1D)
        dims = 1
    class is (grid_2D)
        dims = 2
    class is (grid_3D)
        dims = 3
    end select

     do  proc=n_proc-1,0,-1 !processor #0 last...
         if (proc==this_proc) then

             buf(0) = g%alloc_size
             if (.not.real_wavefunction) then
                 buf(0) = 2*buf(0)
             end if    

             select type(g)
             class is (grid_1D)
                 buf(1) = g%m1min
                 buf(2) = g%m1max
                 buf(3) = g%n1min
                 buf(4) = g%n1max
             class is (grid_2D)
                 buf(1) = g%m1min
                 buf(2) = g%m1max
                 buf(3) = g%n1min
                 buf(4) = g%n1max
                 buf(5) = g%m2min
                 buf(6) = g%m2max
                 buf(7) = g%n2min
                 buf(8) = g%n2max
             class is (grid_3D)
                 buf(1) = g%m1min
                 buf(2) = g%m1max
                 buf(3) = g%n1min
                 buf(4) = g%n1max
                 buf(5) = g%m2min
                 buf(6) = g%m2max
                 buf(7) = g%n2min
                 buf(8) = g%n2max
                 buf(9) = g%m3min
                 buf(10) = g%m3max
                 buf(11) = g%n3min
                 buf(12) = g%n3max
             end select
        end if     
        
        if (this_proc==0) then
             if (proc/=0) then
                 call MPI_Recv(buf, dims*4+1, MPI_INTEGER, proc, 40, MPI_COMM_WORLD, mpi_stat, error)   
                ! call MPI_Recv(u, buf(0), MPI_DOUBLE_PRECISION, proc, 41, MPI_COMM_WORLD, mpi_stat, error)
             end if
             select type(g)
             class is (grid_1D)
                 g1%nn1min = g%nn1min
                 g1%nn1max = g%nn1max
                 g1%m1min = buf(1)
                 g1%m1max = buf(2)
                 g1%n1min = buf(3)
                 g1%n1max = buf(4)
                 call hdf5_read_gridfun_0(file_id, g1, u, real_wavefunction, &
                                           dset_name_real, dset_name_imag, offset)
             class is (grid_2D)
                 g2%nn1min = g%nn1min
                 g2%nn1max = g%nn1max
                 g2%nn2min = g%nn2min
                 g2%nn2max = g%nn2max
                 g2%m1min = buf(1)
                 g2%m1max = buf(2)
                 g2%n1min = buf(3)
                 g2%n1max = buf(4)
                 g2%m2min = buf(5)
                 g2%m2max = buf(6)
                 g2%n2min = buf(7)
                 g2%n2max = buf(8)
                 call hdf5_read_gridfun_0(file_id, g2, u, real_wavefunction, &
                                           dset_name_real, dset_name_imag, offset)
             class is (grid_3D)
                 g3%nn1min = g%nn1min
                 g3%nn1max = g%nn1max
                 g3%nn2min = g%nn2min
                 g3%nn2max = g%nn2max
                 g3%nn3min = g%nn3min
                 g3%nn3max = g%nn3max
                 g3%m1min = buf(1)
                 g3%m1max = buf(2)
                 g3%n1min = buf(3)
                 g3%n1max = buf(4)
                 g3%m2min = buf(5)
                 g3%m2max = buf(6)
                 g3%n2min = buf(7)
                 g3%n2max = buf(8)
                 g3%m3min = buf(9)
                 g3%m3max = buf(10)
                 g3%n3min = buf(11)
                 g3%n3max = buf(12)
                 call hdf5_read_gridfun_0(file_id, g3, u, real_wavefunction, &
                                           dset_name_real, dset_name_imag, offset)
             end select
             if (proc/=0) then
                 call MPI_Send(u, buf(0), MPI_DOUBLE_PRECISION, proc, 41, MPI_COMM_WORLD, error)
             end if   
         else if (proc==this_proc) then
             call MPI_Send(buf, dims*4+1, MPI_INTEGER, 0, 40, MPI_COMM_WORLD, error)
             call MPI_Recv(u, buf(0), MPI_DOUBLE_PRECISION, 0, 41, MPI_COMM_WORLD, mpi_stat, error)
         end if
     end do
end subroutine hdf5_read_gridfun_MPI_0
#endif


#endif

#ifdef _QUADPRECISION_
end module tssmq_hdf5
#else
end module tssm_hdf5
#endif
