module utils_mod

use netcdf
use mpi

implicit none

private

public :: error_handler, sphere_distance, nearest_cell, mark_neighbors_from_source

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
            d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0)
            if (d < nearest_distance) then
               nearest_cell = iCell
               nearest_distance = d
            end if
         end if
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

end module utils_mod
