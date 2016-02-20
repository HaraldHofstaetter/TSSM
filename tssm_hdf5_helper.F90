#ifdef _QUADPRECISION_
module tssmq_hdf5_helper
use tssmq
#else
module tssm_hdf5_helper
use tssm
#endif

use hdf5 
!use h5lt

implicit none

contains

#ifdef _QUADPRECISION_    
function get_quad_type() result(QUAD_TYPE)
    integer(HID_T) :: QUAD_TYPE
    integer :: error 
    call H5Tcopy_f(H5T_IEEE_F64LE, QUAD_TYPE, error)
    call h5Tset_size_f(QUAD_TYPE, 16_8, error)
    call h5Tset_precision_f(QUAD_TYPE, 128_8, error)
    call H5Tset_fields_f(QUAD_TYPE, 127_8, 112_8, 15_8, 0_8, 112_8, error)
    call H5Tset_ebias_f(QUAD_TYPE, 16383_8, error)
end function get_quad_type()
#endif    

subroutine create_file(filename)
    character(len=*), intent(in) :: filename
    integer(HID_T) :: file_id  
    integer :: error
    call h5open_f(error) 
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
end subroutine create_file


subroutine write_array_0(file_id, arrayname, HDF5_FLOAT_TYPE, u, rank, dimensions)
    use, intrinsic :: iso_c_binding
    integer(HID_T), intent(in) :: file_id  
    character(len=*), intent(in) :: arrayname
    integer(hid_t), intent(in) :: HDF5_FLOAT_TYPE
    real(kind=prec), intent(inout), target :: u(*)
    integer, intent(in) :: rank
    integer, intent(in)  :: dimensions(1:rank)
    integer(HSIZE_T)  :: dims(1:rank)
    integer(HID_T) :: space_id  
    integer(HID_T) :: dset_id  
    real(8), pointer :: ur(:)
    integer :: error
    CALL h5dcreate_f(file_id, arrayname, HDF5_FLOAT_TYPE, space_id, dset_id, error)
    call c_f_pointer(c_loc(u), ur, [1])
    dims = dimensions ! convert integer -> integer(HSIZE_T)
    CALL h5dwrite_f(dset_id, HDF5_FLOAT_TYPE, ur, dims, error)
    call hsdclose_f(space_id, error)
    call h5dclose_f(dset_id, error)
end subroutine write_array_0


subroutine write_array(filename, arrayname, u, rank, dimensions)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: arrayname
    real(kind=prec), intent(inout) :: u(*)
    integer, intent(in) :: rank
    integer, intent(in)  :: dimensions(1:rank)
    integer(HID_T) :: file_id  
    integer :: error
    integer(HID_T) :: HDF5_FLOAT_TYPE 
#ifdef _QUADPRECISION_    
    HDF5_FLOAT_TYPE = get_quad_type()  
#else    
    HDF5_FLOAT_TYPE = H5T_IEEE_F64LE 
#endif    
    call h5open_f(error) 
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    call write_array_0(file_id, arrayname, HDF5_FLOAT_TYPE, u, rank, dimensions)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
end subroutine write_array


subroutine write_integer_attr_0(file_id, attrname, u)
    use, intrinsic :: iso_c_binding
    integer(HID_T), intent(in) :: file_id  
    character(len=*), intent(in) :: attrname
    integer, intent(in), target :: u
    integer(HID_T) :: space_id       
    integer(HID_T) :: attr_id 
    integer :: error
    call h5acreate_f (file_id, attrname, H5T_STD_I32LE, space_id, attr_id, error)
    call h5awrite_f (attr_id, H5T_NATIVE_INTEGER, c_loc(u), error)
    call h5sclose_f(space_id, error)
    call h5aclose_f(attr_id, error)
end subroutine write_integer_attr_0


subroutine write_integer_attr(filename, attrname, u)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: attrname
    integer, intent(in) :: u
    integer :: error
    integer(HID_T) :: file_id  
    call h5open_f(error) 
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    call write_integer_attr_0(file_id, attrname, u)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
end subroutine write_integer_attr


subroutine write_real_attr_0(file_id, attrname, u, HDF5_FLOAT_TYPE)
    use, intrinsic :: iso_c_binding
    integer(HID_T), intent(in) :: file_id  
    character(len=*), intent(in) :: attrname
    integer(hid_t), intent(in) :: HDF5_FLOAT_TYPE
    real(kind=prec), intent(in), target :: u
    integer(HID_T) :: space_id       
    integer(HID_T) :: attr_id 
    integer :: error
    call h5acreate_f (file_id, attrname, HDF5_FLOAT_TYPE, space_id, attr_id, error)
    call h5awrite_f (attr_id, HDF5_FLOAT_TYPE, c_loc(u), error)
    call h5sclose_f(space_id, error)
    call h5aclose_f(attr_id, error)
end subroutine write_real_attr_0


subroutine write_real_attr(filename, attrname, u)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: attrname
    real(kind=prec), intent(in) :: u
    integer :: error
    integer(HID_T) :: file_id  
    integer(HID_T) :: HDF5_FLOAT_TYPE
#ifdef _QUADPRECISION_    
    HDF5_FLOAT_TYPE = get_quad_type()  
#else    
    HDF5_FLOAT_TYPE = H5T_IEEE_F64LE 
#endif    
    call h5open_f(error) 
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    call write_real_attr_0(file_id, attrname, u, HDF5_FLOAT_TYPE)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
end subroutine write_real_attr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_array_0(file_id, arrayname, HDF5_FLOAT_TYPE, u, rank, dimensions)
    use, intrinsic :: iso_c_binding
    integer(HID_T), intent(in) :: file_id  
    character(len=*), intent(in) :: arrayname
    integer(hid_t), intent(in) :: HDF5_FLOAT_TYPE
    real(kind=prec), intent(inout), target :: u(*)
    integer, intent(in) :: rank
    integer, intent(in)  :: dimensions(1:rank)
    integer(HSIZE_T) :: dims(1:rank)
    integer(HID_T) :: dset_id  
    real(8), pointer :: ur(:)
    integer :: error
    call c_f_pointer(c_loc(u), ur, [1])
    CALL h5dopen_f(file_id, arrayname, dset_id, error)
    dims = dimensions ! convert integer -> integer(HSIZE_T)
    CALL h5dread_f(dset_id, HDF5_FLOAT_TYPE, ur, dims, error)
    call h5dclose_f(dset_id, error)
end subroutine read_array_0


subroutine read_array(filename, arrayname, u, rank, dimensions)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: arrayname
    real(kind=prec), intent(inout) :: u(*)
    integer, intent(in) :: rank
    integer, intent(in)  :: dimensions(1:rank)
    integer :: error
    integer(HID_T) :: file_id  
    integer(HID_T) :: HDF5_FLOAT_TYPE
#ifdef _QUADPRECISION_    
    HDF5_FLOAT_TYPE = get_quad_type()  
#else    
    HDF5_FLOAT_TYPE = H5T_IEEE_F64LE 
#endif    
    call h5open_f(error) 
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
    call read_array_0(file_id, arrayname, HDF5_FLOAT_TYPE, u, rank, dimensions)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
end subroutine read_array

function read_integer_attr_0(file_id, attrname) result(u)
    integer(HID_T), intent(in) :: file_id  
    character(len=*), intent(in) :: attrname
    integer :: u
    integer(HID_T) :: space_id       
    integer(HID_T) :: attr_id 
    integer(HSIZE_T) :: dummy(1)
    integer :: error
    call h5aopen_f (file_id, attrname, attr_id, error)
    call h5aread_f (attr_id, H5T_NATIVE_INTEGER, u, dummy, error)
    call h5aclose_f(attr_id, error)
end function read_integer_attr_0


function read_integer_attr(filename, attrname) result(u)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: attrname
    integer :: u
    integer(HID_T) :: file_id  
    integer :: error
    call h5open_f(error) 
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
    u =  read_integer_attr_0(file_id, attrname)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
end function read_integer_attr


#ifdef _QUADPRECISION_
end module tssmq_hdf5_helper
#else
end module tssm_hdf5_helper
#endif
