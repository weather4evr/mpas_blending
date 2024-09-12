module unused_stuff

implicit none

private

public :: find_interior_cell, extract_boundary_cells, find_mesh_boundary

type mpas_mesh_type
   integer :: nCells        ! number of cells on the mesh
   integer :: nCellsPerPET  ! number of cells on this processor
   integer :: cell_start    ! starting cell on this processor
   integer :: cell_end      ! ending cell on this processor
   integer :: nEdges        ! number of edges on the mesh
   integer :: maxEdges      ! maximum number of edges on a cell in the mesh
   logical :: regional_mesh ! is it a regional mesh (true) or a global mesh (false)?
   integer, allocatable  :: elemIDs(:) ! IDs of the elements (cell centers) on this processor
   real    :: mpas_mesh_radius ! meters; radius of earth used in MPAS mesh
   real, allocatable :: latCell(:), lonCell(:), areaCell(:) ! each processor has full field
   integer, allocatable :: nEdgesOnCell(:)  ! each processor has full field
   integer, allocatable :: cellsOnCell(:,:) ! each processor has full field
   integer, allocatable :: bdyMaskCell(:)  ! each processor has full field
   integer, allocatable :: cellsOnEdge(:,:) ! each processor has full field
   real, allocatable :: dcEdge(:), dvEdge(:)  ! each processor has full field
end type

public :: mpas_mesh_type

contains

subroutine find_interior_cell(localpet,mpas_mesh,interior_cell)
   integer, intent(in)              :: localpet
   type(mpas_mesh_type), intent(in) :: mpas_mesh
   integer, intent(inout)           :: interior_cell

   integer :: iCell, i, j, k, jCell
   logical :: all_inside

   ! Loop over all Cells.  As soon as we find a point where all its neighbors are also in the
   !  mesh interior, we know that point is in the interior
   do i = 1,mpas_mesh%nCells
      if ( mpas_mesh%bdyMaskCell(i) /= 0 ) cycle ! check the current point; if not interior, go to next one
      all_inside = .true. ! assume all neighbors of the point are interior
      do j = 1, mpas_mesh%nEdgesOnCell(i)
         iCell = mpas_mesh%cellsOnCell(j,i)
         if ( mpas_mesh%bdyMaskCell(iCell) /= 0 ) then
            all_inside = .false.
            exit ! no need to check other cells
         endif
         do k = 1, mpas_mesh%nEdgesOnCell(iCell) ! check another layer
            jCell = mpas_mesh%cellsOnCell(k,iCell)
            if ( mpas_mesh%bdyMaskCell(jCell) /= 0 ) then
               all_inside = .false.
               exit
            endif
         enddo
      enddo
      if ( all_inside ) then
         interior_cell = i
         exit ! we found an interior cell. exit loop.
      endif
   enddo

   if ( localpet == 0 ) then
      write(*,*)'interior cell = ', interior_cell
      write(*,*)'interior cell lat/lon = ', mpas_mesh%latCell(interior_cell)*180.0/3.14, &
                                            mpas_mesh%lonCell(interior_cell)*180.0/3.14
   endif

   return

end subroutine find_interior_cell

subroutine extract_boundary_cells(mpas_mesh,boundary_cell,nbdyCells,bdyCells)
   type(mpas_mesh_type), intent(in) :: mpas_mesh
   integer, intent(in)  :: boundary_cell, nbdyCells
   integer, intent(inout), dimension(nbdyCells) :: bdyCells ! cell IDs where bdyMaskCell == boundary_cell

   integer :: i, n

   n = 0
   do i = 1,mpas_mesh%nCells
      if (mpas_mesh%bdyMaskCell(i) == boundary_cell ) then
         n = n + 1
         bdyCells(n) = i
      endif
   enddo

  ! more succinct way of going it to avoid the loop
  !bdyCells(:) = pack(mpas_mesh%bdyMaskCell, mpas_mesh%bdyMaskCell == boundary_cell)

  ! quick sanity check
  if ( n /= nbdyCells ) then
     call error_handler("Wrong number of boundary cells", -250)
  endif

end subroutine extract_boundary_cells

! given the cell indices of where the boundary is on mesh_source, find the corresponding indices on mesh_target
!  and "mark" them in mesh_target
! the boundary on mesh_source can be find by using subroutine extract_boundary_cells (above)
subroutine find_mesh_boundary(localpet,interior_cell, nbdyCells, bdyCells, mark_neighbors, mesh_source, mesh_target)

   integer, intent(in)                       :: localpet,interior_cell, nbdyCells
   logical, intent(in)                       :: mark_neighbors !fill_holes ! mark_neighbors
   integer, dimension(nbdyCells), intent(in) :: bdyCells
   type(mpas_mesh_type), intent(in)          :: mesh_source
   type(mpas_mesh_type), intent(inout)       :: mesh_target

   integer :: start_cell, i, iCell, nearest_cell_to_source_point, nc, j, jCell, x
   real    :: d
   integer, dimension(nbdyCells) :: nearest_cells, unique_cells
   integer, dimension(nbdyCells*mesh_target%maxEdges) :: indices

   mesh_target%bdyMaskCell(:) = UNMARKED ! these have been filled in define_grid_mpas...but re-initialize them to UNMARKED

   start_cell = 1
   do i = 1,nbdyCells
      iCell = bdyCells(i)
      nearest_cell_to_source_point = nearest_cell( mesh_source%latCell(iCell), mesh_source%lonCell(iCell), &
                                                   start_cell, mesh_target%nCells, mesh_target%maxEdges, &
                                                   mesh_target%nEdgesOnCell, mesh_target%cellsOnCell , &
                                                   mesh_target%latCell, mesh_target%lonCell)
      nc = nearest_cell_to_source_point ! shorter variable name
      nearest_cells(i) = nc

      ! make sure distances bewteen the current and newly-found points are close enough?
     !d = sphere_distance(mesh_source%latCell(iCell), mesh_source%lonCell(iCell), mesh_target%latCell(nc), &
      !                                    mesh_target%lonCell(nc), mesh_source%mpas_mesh_radius )
      ! check the distance between the source lat/lon and target lat/lon
     !if ( d .gt. latlon_tolerance_threshold ) then
     !   call error_handler('oops', 175)
     !endif

      mesh_target%bdyMaskCell(nc) = MARKED ! this is the location of a boundary cell on the target mesh. mark it.

      ! it's possible that there could be a hole in the boundary on one of the meshes
      ! if so, no cells will be masked
      ! fill in the hole by marking the neighbors of points supposedly on the boundary
      if ( mark_neighbors ) then
         do j=1, mesh_target%nEdgesOnCell(nc)
            jCell = mesh_target%cellsOnCell(j, nc)
            mesh_target%bdyMaskCell(jCell) = MARKED
         enddo
      endif

      ! go to the next cell, and make a guess about a good starting cell
      start_cell = nc

   enddo

   ! we now know where the boundary is on the target mesh (probably the next coarsest mesh in this context)
   ! starting from the known interior lat/lon, fill in the mesh
   ! to start, only the boundary is "marked"; everything else is unmarked
   mesh_target%bdyMaskCell(interior_cell) = MARKED
   call mark_neighbors_from_source(interior_cell, MARKED, UNMARKED, mesh_target%bdyMaskCell, mesh_target%cellsOnCell, mesh_target%nEdgesOnCell)

   where ( mesh_target%bdyMaskCell == UNMARKED ) mesh_target%bdyMaskCell = mask_value

end subroutine find_mesh_boundary

end module unused_stuff
