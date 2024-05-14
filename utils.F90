module utils_mod

use netcdf
use mpi

implicit none

private

public :: error_handler

!integer, parameter,public :: DEBUG=0, INFORM=1, WARN=2, ERROR_CODE=3

contains

subroutine error_handler(string, rc)

   character(len=*), intent(in) :: string
   integer,          intent(in) :: rc

   integer :: ierr

   print*,"- FATAL ERROR: ", string
   print*,"- IOSTAT IS: ", rc
   call mpi_abort(mpi_comm_world, 999, ierr)

end subroutine error_handler

end module utils_mod
