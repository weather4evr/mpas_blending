module utils_mod

use netcdf
use mpi
use esmf

implicit none

private

public :: error_handler, sphere_distance, nearest_cell, mark_neighbors_from_source
public :: gather_on_all_procs

!integer, parameter,public :: DEBUG=0, INFORM=1, WARN=2, ERROR_CODE=3

interface gather_on_all_procs
   module procedure gather_on_all_procs_real_1d
   module procedure gather_on_all_procs_real_2d
!  module procedure gather_on_all_procs_int_1d
end interface

contains

subroutine error_handler(string, rc)

   character(len=*), intent(in) :: string
   integer,          intent(in) :: rc

   integer :: ierr

   print*,"- FATAL ERROR: ", string
   print*,"- IOSTAT IS: ", rc
   call mpi_abort(mpi_comm_world, 999, ierr)

end subroutine error_handler

! slightly adapted from a fortran version of limited_area written by Michael Duda (?)
recursive subroutine mark_neighbors_from_source(inside_cell, tval, UNMARKED, mask, cellsOnCell, nEdgesOnCell)
   ! Starting at some cell (inside_cell), expand radially, marking all unmarked cells as tval
   ! until the boundary, which is "marked", is reached.
   ! This assumes that initially, only the boundary is "marked" to value tval, and everything else is "UNMARKED"
   integer, intent(in) :: inside_cell, tval, UNMARKED
   integer, dimension(:), intent(inout) :: mask
   integer, dimension(:), intent(in) :: nEdgesOnCell
   integer, dimension(:,:), intent(in) :: cellsOnCell

   integer :: i, iCell

   do i=1, nEdgesOnCell(inside_cell)
      iCell = cellsOnCell(i, inside_cell)
      if (mask(iCell) == UNMARKED) then
         mask(iCell) = tval
         call mark_neighbors_from_source(iCell, tval, UNMARKED, mask, cellsOnCell, nEdgesOnCell)
      end if
   end do
end subroutine mark_neighbors_from_source

integer function nearest_cell(target_lat, target_lon, start_cell, nCells, maxEdges, &
                              nEdgesOnCell, cellsOnCell, latCell, lonCell)

   real, intent(in) :: target_lat, target_lon
   integer, intent(in) :: start_cell
   integer, intent(in) :: nCells, maxEdges
   integer, dimension(nCells), intent(in) :: nEdgesOnCell
   integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
   real, dimension(nCells), intent(in) :: latCell, lonCell

   integer :: i, iCell
   integer :: current_cell
   real :: current_distance, d, nearest_distance

   nearest_cell = start_cell
   current_cell = -1

   do while (nearest_cell /= current_cell)
      current_cell = nearest_cell
      current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, &
                                           target_lon, 1.0)
      nearest_cell = current_cell
      nearest_distance = current_distance
      do i = 1, nEdgesOnCell(current_cell)
         iCell = cellsOnCell(i,current_cell)
         if (iCell <= nCells) then
         if (iCell > 0 ) then !CSS added this condition
            d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0)
            if (d < nearest_distance) then
               nearest_cell = iCell
               nearest_distance = d
            end if
         end if
         end if ! end CSS condition 
      end do
   end do

end function nearest_cell

real function sphere_distance(lat1, lon1, lat2, lon2, radius)

   real, intent(in) :: lat1, lon1, lat2, lon2, radius
   real :: arg1

   arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
              cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
   sphere_distance = 2.*radius*asin(arg1)

end function sphere_distance

subroutine gather_on_all_procs_real_1d(localpet,input_field,dims,output_field)
   integer, intent(in) :: localpet
   type(ESMF_Field), intent(in) :: input_field
   integer, dimension(1), intent(in) :: dims ! ncells
   real*8, dimension(:), intent(inout) :: output_field

   real(esmf_kind_r8), allocatable :: dum(:)
   integer :: rc, nCells

   nCells = dims(1)

   if (localpet==0) then
      allocate(dum(nCells))
   else
      allocate(dum(0))
   endif

   ! gather onto process 0
   call ESMF_FieldGather(input_field, dum, rootPet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

   if (localpet==0) write(*,*)'min/max0 = ',minval(dum),maxval(dum)

   ! dum is real(esmf_kind_r8). We need to get dum onto all processors
   ! (currently it's just on process 0)
   ! So, using mpi_bcast is the easiest way to do it, but mpi_bcast doesn't
   ! have a real "type" of esmf_kind_r8.  So we need to use a different
   ! variable to do the broadcasting that is mpi_real, hence,
   ! output_field is type real*8
   if (localpet==0) output_field(:) = dum(:)

   call mpi_bcast(output_field,nCells,mpi_real8,0,mpi_comm_world,rc)
   !write(*,*)'pe/min/max = ',localpet,minval(output_field),maxval(output_field)
   
   deallocate(dum)

end subroutine gather_on_all_procs_real_1d

subroutine gather_on_all_procs_real_2d(localpet,input_field,dims,output_field)
   integer, intent(in) :: localpet
   type(ESMF_Field), intent(in) :: input_field
   integer, dimension(2), intent(in) :: dims ! (ncells, nz)
   real*8, dimension(:,:), intent(inout) :: output_field

   real(esmf_kind_r8), allocatable :: dum(:,:)
   integer :: rc, nCells, nz

   nCells = dims(1)
   nz = dims(2)

   if (localpet==0) then
      allocate(dum(nCells,nz))
   else
      allocate(dum(0,0))
   endif

   ! gather onto process 0
   call ESMF_FieldGather(input_field, dum, rootPet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

   if (localpet==0) write(*,*)'min/max0 = ',minval(dum),maxval(dum)

   ! dum is real(esmf_kind_r8). We need to get dum onto all processors
   ! (currently it's just on process 0)
   ! So, using mpi_bcast is the easiest way to do it, but mpi_bcast doesn't
   ! have a real "type" of esmf_kind_r8.  So we need to use a different
   ! variable to do the broadcasting that is mpi_real, hence,
   ! output_field is type real*8
   if (localpet==0) output_field(:,:) = dum(:,:)

   call mpi_bcast(output_field,nCells*nz,mpi_real8,0,mpi_comm_world,rc)
   !write(*,*)'pe/min/max = ',localpet,minval(output_field),maxval(output_field)
   
   deallocate(dum)

end subroutine gather_on_all_procs_real_2d

!subroutine gather_on_all_procs_int_1d(localpet,dims,input_field,output_field)
!   integer, intent(in) :: localpet
!   integer, dimension(2), intent(in) :: dims ! (ncells, nz) ; nz is optional
!   type(ESMF_Field), intent(in) :: input_field
!   integer, intent(inout), dimension(dims(1)) :: output_field
!
!   integer, allocatable :: dum1d(:)
!   integer :: rc, nCells
!
!   nCells = dims(1)
!
!   if (localpet==0) then
!      allocate(dum1d(nCells))
!   else
!      allocate(dum1d(0))
!   endif
!
!   ! gather onto process 0
!   call ESMF_FieldGather(input_field, dum1d, rootPet=0, rc=rc)
!   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
!      call error_handler("IN FieldGather", rc)
!
!   if (localpet==0) write(*,*)'min/max0 = ',minval(dum1d),maxval(dum1d)
!
!   if (localpet==0) output_field(:) = dum1d(:)
!
!   call mpi_bcast(output_field,nCells,mpi_integer,0,mpi_comm_world,rc)
!   !write(*,*)'pe/min/max = ',localpet,minval(output_field),maxval(output_field)
!   
!   deallocate(dum1d)
!
!end subroutine gather_on_all_procs_int_1d

end module utils_mod
