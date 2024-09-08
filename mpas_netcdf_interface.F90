module mpas_netcdf_interface

use netcdf
use mpi
use esmf
use utils_mod, only  : error_handler

implicit none

private

! public subroutines
public :: netcdf_err
public :: open_netcdf, close_netcdf
public :: get_netcdf_dims
public :: get_netcdf_var

interface get_netcdf_var
   module procedure get_netcdf_var_1d_integer
!  module procedure get_netcdf_var_1d_integer_esmf_kind
   module procedure get_netcdf_var_1d_real
   module procedure get_netcdf_var_1d_real_esmf_kind
   module procedure get_netcdf_var_2d_integer
   module procedure get_netcdf_var_2d_real
   module procedure get_netcdf_var_2d_real_esmf_kind
end interface

! variables visible to this module only
logical :: fexist
integer :: error

contains

subroutine open_netcdf(fname,ncfileid)
   character(len=*), intent(in) :: fname
   integer, intent(out) :: ncfileid

   inquire(file=trim(adjustl(fname)),exist=fexist)
   if ( .not. fexist ) then
      call error_handler(trim(adjustl(fname))//' does not exist',-1)
   endif

!  write(*,fmt='(a)') 'Opening '//trim(adjustl(fname))

   error = nf90_open(path=trim(adjustl(fname)),mode=nf90_nowrite,ncid=ncfileid)
   call netcdf_err(error,'error reading '//trim(adjustl(fname)) )
      
end subroutine open_netcdf

subroutine close_netcdf(fname,ncfileid)
   character(len=*), intent(in) :: fname
   integer, intent(in) :: ncfileid
   error = nf90_close(ncfileid) ! close file
   call netcdf_err(error,'error closing '//trim(adjustl(fname)) )
end subroutine close_netcdf

subroutine get_netcdf_dims(fileid,variable,output)
   integer, intent(in) :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(out) :: output

   integer :: ncdimid

   error = nf90_inq_dimid(fileid,trim(adjustl(variable)),ncdimid)
   call netcdf_err(error,'Error nf90_inq_dimid '//trim(adjustl(variable)) )
   error = nf90_inquire_dimension(fileid,ncdimid, len=output)
   call netcdf_err(error,'Error nf90_inquire_dimension '//trim(adjustl(variable)) )
end subroutine get_netcdf_dims

subroutine get_netcdf_var_1d_integer(fileid,variable,start_index,num_points,output)

   integer, intent(in)          :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(in)          :: start_index(1),num_points(1) ! these are arrays (e.g., (/cell_start/), (/nCells/))
   integer, intent(inout), dimension(num_points(1))  :: output

   integer :: ncvarid

   error = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   call netcdf_err(error,'Error nf90_inq_varid '//trim(adjustl(variable)) )
   error = nf90_get_var(fileid,ncvarid,start=start_index, count=num_points, values=output)
   call netcdf_err(error,'Error nf90_get_var '//trim(adjustl(variable)) )

end subroutine get_netcdf_var_1d_integer

subroutine get_netcdf_var_1d_integer_esmf_kind(fileid,variable,start_index,num_points,output)

   integer, intent(in)          :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(in)          :: start_index(1),num_points(1) ! these are arrays (e.g., (/cell_start/), (/nCells/))
   integer(ESMF_KIND_I4), intent(inout), dimension(num_points(1))  :: output

   integer :: ncvarid

   error = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   call netcdf_err(error,'Error nf90_inq_varid '//trim(adjustl(variable)) )
   error = nf90_get_var(fileid,ncvarid,start=start_index, count=num_points, values=output)
   call netcdf_err(error,'Error nf90_get_var '//trim(adjustl(variable)) )

end subroutine get_netcdf_var_1d_integer_esmf_kind

subroutine get_netcdf_var_1d_real(fileid,variable,start_index,num_points,output)

   integer, intent(in)          :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(in)          :: start_index(1),num_points(1) ! these are arrays (e.g., (/cell_start/), (/nCells/))
   real, intent(inout), dimension(num_points(1))  :: output

   integer :: ncvarid

   error = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   call netcdf_err(error,'Error nf90_inq_varid '//trim(adjustl(variable)) )
   error = nf90_get_var(fileid,ncvarid,start=start_index, count=num_points, values=output)
   call netcdf_err(error,'Error nf90_get_var '//trim(adjustl(variable)) )

end subroutine get_netcdf_var_1d_real

subroutine get_netcdf_var_1d_real_esmf_kind(fileid,variable,start_index,num_points,output)

   integer, intent(in)          :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(in)          :: start_index(1),num_points(1) ! these are arrays (e.g., (/cell_start/), (/nCells/))
   real(esmf_kind_r8), intent(inout), dimension(num_points(1))  :: output

   integer :: ncvarid

   error = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   call netcdf_err(error,'Error nf90_inq_varid '//trim(adjustl(variable)) )
   error = nf90_get_var(fileid,ncvarid,start=start_index, count=num_points, values=output)
   call netcdf_err(error,'Error nf90_get_var '//trim(adjustl(variable)) )

end subroutine get_netcdf_var_1d_real_esmf_kind

subroutine get_netcdf_var_2d_integer(fileid,variable,start_index,num_points,output)

   integer, intent(in)          :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(in)          :: start_index(2),num_points(2) ! these are arrays (e.g., (/cell_start/), (/nCells/))
   integer, intent(inout), dimension(num_points(1),num_points(2))  :: output

   integer :: ncvarid

   error = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   call netcdf_err(error,'Error nf90_inq_varid '//trim(adjustl(variable)) )
   error = nf90_get_var(fileid,ncvarid,start=start_index, count=num_points, values=output)
   call netcdf_err(error,'Error nf90_get_var '//trim(adjustl(variable)) )

end subroutine get_netcdf_var_2d_integer

subroutine get_netcdf_var_2d_real(fileid,variable,start_index,num_points,output)

   integer, intent(in)          :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(in)          :: start_index(2),num_points(2) ! these are arrays (e.g., (/cell_start/), (/nCells/))
   real, intent(inout), dimension(num_points(1),num_points(2))  :: output

   integer :: ncvarid

   error = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   call netcdf_err(error,'Error nf90_inq_varid '//trim(adjustl(variable)) )
   error = nf90_get_var(fileid,ncvarid,start=start_index, count=num_points, values=output)
   call netcdf_err(error,'Error nf90_get_var '//trim(adjustl(variable)) )

end subroutine get_netcdf_var_2d_real

subroutine get_netcdf_var_2d_real_esmf_kind(fileid,variable,start_index,num_points,output)

   integer, intent(in)          :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(in)          :: start_index(2),num_points(2) ! these are arrays (e.g., (/cell_start/), (/nCells/))
   real(esmf_kind_r8), intent(inout), dimension(num_points(1),num_points(2))  :: output

   integer :: ncvarid

   error = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   call netcdf_err(error,'Error nf90_inq_varid '//trim(adjustl(variable)) )
   error = nf90_get_var(fileid,ncvarid,start=start_index, count=num_points, values=output)
   call netcdf_err(error,'Error nf90_get_var '//trim(adjustl(variable)) )

end subroutine get_netcdf_var_2d_real_esmf_kind
!! Error handler for netcdf
!! @param[in] err     error status code
!! @param[in] string  error message
subroutine netcdf_err( err, string )

   integer, intent(in) :: err
   character(len=*), intent(in) :: string

   character(len=256) :: errmsg
   integer :: iret

   if( err.EQ.NF90_NOERR )return
   errmsg = NF90_STRERROR(err)
   print*,''
   print*,'FATAL ERROR: ', trim(string), ': ', trim(errmsg)
   print*,'STOP.'
   call mpi_abort(mpi_comm_world, 999, iret)

end subroutine netcdf_err
 
end module mpas_netcdf_interface
