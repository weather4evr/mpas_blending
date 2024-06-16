module model_grid

use mpi
use netcdf
use esmf
use ESMF_LogPublicMod
use utils_mod, only      : error_handler, sphere_distance, nearest_cell, mark_neighbors_from_source
use mpas_netcdf_interface, only : open_netcdf, close_netcdf, get_netcdf_dims, &
                                  get_netcdf_var, netcdf_err
use program_setup, only  : nvars_to_blend, &
                           nlat, nlon, lat_ll, lon_ll, dx_in_degrees, &
                           is_regional, output_latlon_grid, max_boundary_layer_to_keep
implicit none

private

! Public variables accessible to other modules
character(len=500), allocatable, public :: grid_files_heirarchy(:)
real, allocatable, public              :: nominal_horizontal_cell_spacing(:) ! km; nominal dx of MPAS mesh
real, allocatable, public              :: weights(:) ! weight to give to model providing large-scale data
type(esmf_mesh),  allocatable, public  :: meshes(:) ! esmf_mesh that will hold the heirarchy of MPAS meshes
type(esmf_grid),  public               :: latlon_grid ! used if interpolating to lat-lon grid
integer, public               :: nmeshes
integer, allocatable, public  :: elemIDs(:) !< IDs of the elements on present PET
real, allocatable, public     :: lons_output_grid(:,:), lats_output_grid(:,:)
integer, parameter, public    :: mask_value = -9999
integer, parameter, public :: INSIDE = 0
integer, parameter, public :: UNMARKED = -1
integer, parameter, public :: MARKED = 1


! Public subroutines
public :: read_grid_info_file
public :: define_grid_mpas, define_grid_latlon
public :: find_interior_cell, extract_boundary_cells, find_mesh_boundary
public :: cleanup

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
type(mpas_mesh_type), public, allocatable :: mpas_meshes(:)

contains

subroutine define_grid_mpas(localpet,npets,the_file,grid_mesh,mpas_mesh)

   integer, intent(in)            :: localpet, npets
   character(len=500), intent(in) :: the_file
   type(esmf_mesh), intent(inout) :: grid_mesh ! esmf_mesh that is output
   type(mpas_mesh_type), intent(inout) :: mpas_mesh

   integer                      :: error, i, j, k, n
   integer                      :: ncid,id_var, id_dim, nVertThis
   integer                      :: nCells, nVertices, maxEdges, nEdges
   integer                      :: cell_start, cell_end, temp(1)
   integer                      :: nCellsPerPET ! Number of cells on this PET
   integer, allocatable         :: elemTypes2(:), vertOnCell(:,:), &
                                            nodesPET(:), nodeIDs(:), nodeIDs_temp(:), &
                                            elementConn_temp(:), elementConn(:)
   integer, allocatable             :: bdyMaskCell(:)
   real(esmf_kind_r8), allocatable  :: latCell(:), lonCell(:), &
                                       latVert(:), lonVert(:), &
                                       nodeCoords(:), &
                                       nodeCoords_temp(:), &
                                       elemCoords(:)
   real(esmf_kind_r8), parameter    :: PI=4.D0*DATAN(1.D0)

   if (localpet==0) print*,'- OPEN MPAS INPUT FILE: ',trim(the_file)
   call open_netcdf(trim(the_file),ncid)

   !Get some dimensions
   call get_netcdf_dims(ncid,'nCells',nCells)
   call get_netcdf_dims(ncid,'nVertices',nVertices)
   call get_netcdf_dims(ncid,'maxEdges',maxEdges)
   call get_netcdf_dims(ncid,'nEdges',nEdges)

   ! original algorithm from Larissa, which can lead to
   !  nCellsPerPET < 0 for some meshes,
   !  probably because of using "ceiling()"
!  nCellsPerPET = ceiling(real(nCells)/real(npets))
!  cell_start = localpet*nCellsPerPET+1
!  cell_end = min(localpet*nCellsPerPET+nCellsPerPET,nCells)
!  nCellsPerPET = cell_end - cell_start + 1

   ! Number of MPAS cells per processor, and starting/ending cells.
   call para_range(1, nCells, npets, localpet, cell_start, cell_end)
   nCellsPerPET = cell_end - cell_start + 1

   ! Fill up the structure...this information is useful to be able to carry around
   mpas_mesh%nCells        = nCells
   mpas_mesh%nCellsPerPET  = nCellsPerPET
   mpas_mesh%cell_start    = cell_start
   mpas_mesh%cell_end      = cell_end
   mpas_mesh%maxEdges      = maxEdges
   mpas_mesh%nEdges        = nEdges
!  error = nf90_get_att(ncid, NF90_GLOBAL, 'sphere_radius', mpas_mesh%mpas_mesh_radius) ! meters
!  call netcdf_err(error, 'reading sphere_radius global attribute')
   mpas_mesh%mpas_mesh_radius = 6371229.0 ! meters
   if ( localpet == 0 ) write(*,*)' sphere_radius = ',mpas_mesh%mpas_mesh_radius

   allocate(mpas_mesh%latCell(1:nCells))
   allocate(mpas_mesh%lonCell(1:nCells))
   allocate(latCell(cell_start:cell_end))
   allocate(lonCell(cell_start:cell_end))
   allocate(latVert(nVertices))
   allocate(lonVert(nVertices))
   allocate(vertOnCell(maxEdges,cell_start:cell_end))

   ! GET CELL CENTER LAT/LON
   call get_netcdf_var(ncid,'latCell',(/cell_start/),(/nCellsPerPET/),latCell)
   call get_netcdf_var(ncid,'lonCell',(/cell_start/),(/nCellsPerPET/),lonCell)

   ! -- Put full latitude/longitude into the mpas_mesh structure
   call get_netcdf_var(ncid,'latCell',(/1/),(/nCells/),mpas_mesh%latCell)
   call get_netcdf_var(ncid,'lonCell',(/1/),(/nCells/),mpas_mesh%lonCell)

   ! GET VERTEX LAT/LON
   call get_netcdf_var(ncid,'lonVertex',(/1/),(/nVertices/),lonVert)
   call get_netcdf_var(ncid,'latVertex',(/1/),(/nVertices/),latVert)

   if (localpet==0) print*,"- NUMBER OF CELLS ON INPUT GRID ", nCells
   if (localpet==0) print*,"- NUMBER OF NODES ON INPUT GRID ", nVertices

  !error=nf90_get_var(ncid, id_var, start=(/1,cell_start/),count=(/maxEdges,nCellsPerPET/),values=vertOnCell)
   call get_netcdf_var(ncid,'verticesOnCell',(/1,cell_start/),(/maxEdges,nCellsPerPET/),vertOnCell)

   allocate(mpas_mesh%cellsOnCell(maxEdges,nCells))
   allocate(mpas_mesh%nEdgesOnCell(nCells))
   allocate(mpas_mesh%areaCell(nCells))
   allocate(mpas_mesh%cellsOnEdge(2,nEdges))
   allocate(mpas_mesh%dcEdge(nEdges))
   allocate(mpas_mesh%dvEdge(nEdges))
   call get_netcdf_var(ncid,'cellsOnCell',(/1,1/),(/maxEdges,nCells/),mpas_mesh%cellsOnCell)
   call get_netcdf_var(ncid,'nEdgesOnCell',(/1/),(/nCells/),mpas_mesh%nEdgesOnCell)
   call get_netcdf_var(ncid,'areaCell',(/1/),(/nCells/),mpas_mesh%areaCell)
   call get_netcdf_var(ncid,'cellsOnEdge',(/1,1/),(/2,nEdges/),mpas_mesh%cellsOnEdge)
   call get_netcdf_var(ncid,'dcEdge',(/1/),(/nEdges/),mpas_mesh%dcEdge)
   call get_netcdf_var(ncid,'dvEdge',(/1/),(/nEdges/),mpas_mesh%dvEdge)

   ! Determine if the mesh is regional. If so, 
   !  fill it and mask out the boundary points.
   !  Use get variable ID--use direct nf90 call here because we just want the id
   ! Want to get the field both for this processor only, and the full field to all processors
   allocate(bdyMaskCell(nCellsPerPET)) ! integer
   allocate(mpas_mesh%bdyMaskCell(nCells)) ! integer
   error=nf90_inq_varid(ncid, 'bdyMaskCell',id_var)
   if ( error == 0 ) then
      if (localpet==0) print*,"- Found bdyMaskCell"
      call get_netcdf_var(ncid,'bdyMaskCell',(/cell_start/),(/nCellsPerPET/),bdyMaskCell)
      where (bdyMaskCell > max_boundary_layer_to_keep ) bdyMaskCell = mask_value

      call get_netcdf_var(ncid,'bdyMaskCell',(/1/),(/nCells/),mpas_mesh%bdyMaskCell)
      where (mpas_mesh%bdyMaskCell > max_boundary_layer_to_keep ) mpas_mesh%bdyMaskCell = mask_value
   else
      bdyMaskCell = 0
      mpas_mesh%bdyMaskCell = 0
   endif

   ! Global meshes could have bdyMaskCell, but it will be 0 everywhere.  Figure out if 
   ! it's a regional mesh
   if ( all(mpas_mesh%bdyMaskCell == 0 )) then
      mpas_mesh%regional_mesh = .false.
   else
      mpas_mesh%regional_mesh = .true.
   endif
   if (localpet == 0) print*,"- Regional mesh is ",mpas_mesh%regional_mesh

   ! If true, we just let ESMF do its thing in the interpolation,
   !   and it will not interpolate to points that it can't map to
  !if ( use_esmf_instinctive_masking_for_regional ) then
  !   bdyMaskCell = 0
  !   mpas_mesh%bdyMaskCell = 0
  !endif

   ! Close netCDF file--no longer needed
   call close_netcdf(trim(the_file),ncid)
     
   ! Allocate and fill element corner coordinate array.
   allocate(nodeCoords_temp(2*maxEdges*nCellsPerPET))
   allocate(nodeIDs_temp(maxEdges*nCellsPerPET))
   allocate(elementConn_temp(maxEdges*nCellsPerPET))
   allocate(elemCoords(2*nCellsPerPET))
   allocate(elemTypes2(nCellsPerPET))
   allocate(elemIDs(nCellsPerPET))

   nVertThis = 0
   k = 0

   !------------------------------------
   ! Define the MPAS mesh in ESMF terms
   !------------------------------------
   if (localpet==0) print*, "- Create PET-local element connectivity "
   do i = cell_start,cell_end
       j = i - cell_start + 1
       elemIDs(j) = i
       elemTypes2(j) = 0
       elemCoords(2*j-1) = lonCell(i)*180.0_esmf_kind_r8/PI
       if (elemCoords(2*j-1) > 180.0_esmf_kind_r8) then
           elemCoords(2*j-1) = elemCoords(2*j-1) - 360.0_esmf_kind_r8
       endif
       elemCoords(2*j) = latCell(i)*180.0_esmf_kind_r8/PI
       do n = 1,maxEdges
           if (vertOnCell(n,i)>0) then
               nVertThis = nVertThis + 1

               ! Make sure we don't duplicate nodeIDs or nodeCoords on any PET
               if (.not. any(nodeIDs_temp(1:k)==vertOnCell(n,i))) then
                   k = k+1
                   nodeCoords_temp(2*k-1) = &
                       lonVert(vertOnCell(n,i)) *180.0_esmf_kind_r8/PI
                   if (nodeCoords_temp(2*k-1) > 180.0_esmf_kind_r8) then
                       nodeCoords_temp(2*k-1)  =  &
                           nodeCoords_temp(2*k-1) - 360.0_esmf_kind_r8
                   endif
                   nodeCoords_temp(2*k) = &
                       latVert(vertOnCell(n,i))*180.0_esmf_kind_r8/PI

                   nodeIDs_temp(k) = vertOnCell(n,i)
               endif

               elemTypes2(j) = elemTypes2(j) + 1

               temp = FINDLOC(nodeIDS_temp, vertOnCell(n,i))
               elementConn_temp(nVertThis) = temp(1)

           endif

       enddo
   enddo

   allocate(nodeCoords(2*k), nodeIDs(k), elementConn(nVertThis))
   nodeCoords = nodeCoords_temp(1:k*2)
   nodeIDs = nodeIDs_temp(1:k)
   elementConn = elementConn_temp(1:nVertThis)

   !-----------------------------------------------
   ! Create ESMF mesh object for the model grid.
   !-----------------------------------------------
   if (localpet==0) print*, "- CREATE MESH -"
   grid_mesh  = ESMF_MeshCreate(parametricDim=2, &
                       spatialDim=2, &
                       nodeIDs= nodeIDs, &
                       nodeCoords = nodeCoords, &
                       elementIDs = elemIDs, &
                       elementTypes=elemTypes2, &
                       elementConn = elementConn, &
                       elementCoords=elemCoords, &
                       coordSys=ESMF_COORDSYS_SPH_DEG, &
                       elementMask=bdyMaskCell, &
                       rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN MeshCreate", error)

   ! also save elemIDs, in the mpas_mesh_type structure
   allocate(mpas_mesh%elemIDs(nCellsPerPET))
   mpas_mesh%elemIDs = elemIDs

   ! After mesh creation we are through with the arrays, so they may be deallocated.
   deallocate(elemCoords,elemTypes2)
   deallocate(nodeCoords_temp, nodeCoords)
   deallocate(elementConn_temp, elementConn)
   deallocate(nodeIDs_temp)
   deallocate(lonCell, latCell, lonVert, latVert)
   deallocate(vertOnCell)
   deallocate(elemIDs) ! we can deallocate this b/c we are saving it in mpas_mesh
   deallocate(nodeIDs) ! this isn't used later so we can deallocate
   deallocate(bdyMaskCell)

   if (localpet==0) then 
      write(*,*)'Done with define_grid_mpas for '//trim(adjustl(the_file))
      write(*,*)''
   endif

end subroutine define_grid_mpas

subroutine define_grid_latlon(localpet,npets)

   integer, intent(in)          :: localpet, npets

   character(len=500)           :: the_file

   integer                      :: id_dim, id_var, i, j
   integer                      :: i_input, j_input
   integer                      :: ip1_input, jp1_input
   integer                      :: error
   integer                      :: clb(2), cub(2), starts(2), counts(2)
   real(esmf_kind_r8)           :: half_dx_in_degrees
   real(esmf_kind_r8), allocatable :: templat(:), templon(:)
   real(esmf_kind_r8), pointer     :: lat_src_ptr(:,:), lon_src_ptr(:,:), &
                                      clat_src_ptr(:,:), clon_src_ptr(:,:)
   type(esmf_polekind_flag)        :: polekindflag(2)

   i_input = nlon ! number of x points
   j_input = nlat ! number of y points
   allocate(templon(1:i_input))
   allocate(templat(1:j_input))
   do i = 1,i_input
      templon(i) = lon_ll + dx_in_degrees*(i-1)
   enddo
   do j = 1,j_input
      templat(j) = lat_ll + dx_in_degrees*(j-1)
   enddo
   if(localpet==0) print*,"- I/J DIMENSIONS OF OUTPUT LAT/LON FILE ", i_input, j_input
   if(localpet==0) write(*,*)'min/max lat',minval(templat),maxval(templat)
   if(localpet==0) write(*,*)'min/max lon',minval(templon),maxval(templon)

   ip1_input = i_input + 1
   jp1_input = j_input + 1

   ! basically follow subroutine define_input_grid_gaussian in chgres_cube/model_grid.F90
   if (is_regional) then
      if (localpet==0) print*,"- CALL GridCreateNoPeriDim FOR INPUT GRID"
      latlon_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/i_input,j_input/), &
                                         indexflag=ESMF_INDEX_GLOBAL, &
                                         regDecomp=(/1,npets/), &
                                         rc=error)
   else
      polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE
      if (localpet==0)print*,"- CALL GridCreate1PeriDim FOR INPUT GRID."
      latlon_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
             maxIndex=(/i_input,j_input/), &
             polekindflag=polekindflag, &
             periodicDim=1, &
             poleDim=2,  &
             coordSys=ESMF_COORDSYS_SPH_DEG, &
             regDecomp=(/1,npets/),  &
             indexflag=ESMF_INDEX_GLOBAL, rc=error)

   endif
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN GridCreateNoPeriDim", error)

   !-----------------------------------------------------------------------
   ! Create needed grid coordinates
   !-----------------------------------------------------------------------
      
   if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID CENTER."
   call ESMF_GridAddCoord(latlon_grid, &
                          staggerloc=ESMF_STAGGERLOC_CENTER, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN GridAddCoord", error)
      
   if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID CORNER."
   call ESMF_GridAddCoord(latlon_grid, &
                          staggerloc=ESMF_STAGGERLOC_CORNER, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN GridAddCoord", error)
      
   !----------------------------------------------------------
   !  Read in coordinate values and set on grid
   !----------------------------------------------------------
      
   !------- Grid center coordinates
      
   if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."

   nullify(lon_src_ptr)
   nullify(lat_src_ptr)

   call ESMF_GridGetCoord(latlon_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=1, &
                        farrayPtr=lon_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
     call error_handler("IN GridGetCoord", error)

   if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
   call ESMF_GridGetCoord(latlon_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
     call error_handler("IN GridGetCoord", error)

   !allocate(templon(clb(1):cub(1)))
   !allocate(templat(clb(2):cub(2)))
   !starts = (/clb(1),clb(2)/)
   !counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)

   do j = clb(2),cub(2)
     do i = clb(1), cub(1)
       lon_src_ptr(i,j)=real(templon(i),esmf_kind_r8)
       if (lon_src_ptr(i,j) > 360.0_esmf_kind_r8) lon_src_ptr(i,j) = lon_src_ptr(i,j) - 360.0_esmf_kind_r8
       lat_src_ptr(i,j)=real(templat(j),esmf_kind_r8)
     enddo
   enddo

   ! define global variables; only needed on localpet == 0
   !  in case we output a latlon grid
   if (localpet==0) then
      allocate(lons_output_grid(nlon,nlat))
      allocate(lats_output_grid(nlon,nlat))
      do j = 1,nlat
         do i = 1,nlon
            lons_output_grid(i,j) = templon(i)
            if ( lons_output_grid(i,j) > 360.0 ) lons_output_grid(i,j) = lons_output_grid(i,j) - 360.0
            lats_output_grid(i,j) = templat(j)
         enddo
      enddo
   endif

   nullify(lon_src_ptr)
   nullify(lat_src_ptr)
    
   !---------- Grid corners coordinate creation
   nullify(clat_src_ptr)
   nullify(clon_src_ptr)

   if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT CORNERS GRID X-COORD."
     call ESMF_GridGetCoord(latlon_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CORNER, &
                          coordDim=1, &
                          farrayPtr=clon_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
     call error_handler("IN GridGetCoord", error)

   if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT CORNERS GRID Y-COORD."
   call ESMF_GridGetCoord(latlon_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CORNER, &
                          coordDim=2, &
                          computationalLBound=clb, &
                          computationalUBound=cub, &
                          farrayPtr=clat_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN GridGetCoord", error)
        
   ! For a lat/lon grid, it's easy to get the corner lat/lons from the cell center lat/lons
   ! clb, cub are bounds are for the corner points
   !print *,localpet, clb(1), cub(1), clb(2), cub(2)

   half_dx_in_degrees = 0.5*dx_in_degrees

   do j = clb(2),cub(2)
     do i = clb(1), cub(1)
       if ( i .eq. ip1_input .and. j .eq. jp1_input ) then ! top-right corner
          clon_src_ptr(i,j)=real(templon(i_input)+half_dx_in_degrees,esmf_kind_r8)
          clat_src_ptr(i,j)=real(templat(j_input)+half_dx_in_degrees,esmf_kind_r8)
       else if ( i .eq. ip1_input ) then ! points on right edge
          clon_src_ptr(i,j)=real(templon(i_input)+half_dx_in_degrees,esmf_kind_r8)
          clat_src_ptr(i,j)=real(templat(j)-half_dx_in_degrees,esmf_kind_r8)
       else if ( j .eq. jp1_input ) then ! points on top edge
          clon_src_ptr(i,j)=real(templon(i)-half_dx_in_degrees,esmf_kind_r8)
          clat_src_ptr(i,j)=real(templat(j_input)+half_dx_in_degrees,esmf_kind_r8)
       else ! all other points
          clon_src_ptr(i,j)=real(templon(i)-half_dx_in_degrees,esmf_kind_r8)
          clat_src_ptr(i,j)=real(templat(j)-half_dx_in_degrees,esmf_kind_r8)
       endif
       if (clon_src_ptr(i,j) > 360.0_esmf_kind_r8) clon_src_ptr(i,j) = clon_src_ptr(i,j) - 360.0_esmf_kind_r8
       if (clon_src_ptr(i,j) <   0.0_esmf_kind_r8) clon_src_ptr(i,j) = clon_src_ptr(i,j) + 360.0_esmf_kind_r8
       if (clat_src_ptr(i,j) >  90.0_esmf_kind_r8) clat_src_ptr(i,j) = 90.0_esmf_kind_r8
       if (clat_src_ptr(i,j) < -90.0_esmf_kind_r8) clat_src_ptr(i,j) = -90.0_esmf_kind_r8
     enddo
   enddo

   nullify(clon_src_ptr)
   nullify(clat_src_ptr)
   deallocate(templat,templon)

end subroutine define_grid_latlon

subroutine cleanup(localpet)

   integer, intent(in)              :: localpet
   integer                          :: rc, i

   if (localpet==0) print*,"- DESTROY MODEL DATA."

   do i = 1, nmeshes
      call ESMF_MeshDestroy(meshes(i), rc=rc)
   enddo
   deallocate(meshes)

   if ( output_latlon_grid ) call ESMF_GridDestroy(latlon_grid, rc=rc)

   deallocate(grid_files_heirarchy,nominal_horizontal_cell_spacing,weights)
   if (localpet==0 .and. output_latlon_grid ) then
      deallocate(lons_output_grid, lats_output_grid)
   endif

!  deallocate(elemIDs, nodeIDs) ! deallocated in define_mpas_grid
!  call ESMF_FieldDestroy(latitude_target_grid, rc=rc)

end subroutine cleanup

! this fills public variables defined in this module
subroutine read_grid_info_file(localpet,the_file)

   integer, intent(in)           :: localpet
   character(len=*), intent(in)  :: the_file

   integer :: i, k, istat, nlines
   character(len=200) :: line
   logical :: fexist

   ! First make sure the file is there
   inquire(file=trim(adjustl(the_file)),exist=fexist)
   if ( .not. fexist ) then
      call error_handler("File "//trim(adjustl(the_file))//" is missing",-110)
   endif

   ! Now open the file and read it
   open(14, file=trim(adjustl(the_file)), form='formatted', iostat=istat)
   if (istat /= 0) call error_handler("OPENING GRID INFO FILE", istat)

   nmeshes = 0
   nlines = 0

   ! Loop over lines of file to count the number of lines and meshes
   do
     read(14, '(A)', iostat=istat) line
     if (istat/=0) exit
     nlines = nlines + 1 
     if ( trim(line) .eq. '' ) cycle
     if ( trim(line(1:1)) .eq. '!' .or. trim(line(1:1)) .eq. '#' ) cycle
     nmeshes = nmeshes+1
   enddo

   if (localpet==0) print*, "READING ", nlines, " LINES ACCORDING TO ", trim(adjustl(the_file))
   if (localpet==0) print*, "READING ", nmeshes, " MESHES ACCORDING TO ", trim(adjustl(the_file))
   if ( nmeshes == 0) call error_handler("GRID_INFO FILE IS EMPTY.", -1)

   allocate(grid_files_heirarchy(nmeshes))
   allocate(nominal_horizontal_cell_spacing(nmeshes))
   allocate(weights(nmeshes))
   weights = -1.0 ! initialize it to something that's nonsense

   ! Now read the file again and pull out the grid information
   rewind(14)
   i = 0
   do k = 1,nlines
      read(14, '(A)', iostat=istat) line
      if (istat/=0) exit
      if ( trim(line) .eq. '' ) cycle
      if ( trim(line(1:1)) .eq. '!' .or. trim(line(1:1)) .eq. '#' ) cycle
      backspace(14) ! this line is valid and has a grid. go back so we can read it again.
      i = i + 1
      read(14, *, iostat=istat) grid_files_heirarchy(i), nominal_horizontal_cell_spacing(i), weights(i)
      if (istat /= 0) call error_handler("READING GRID_INFO FILE", istat)
      grid_files_heirarchy(i) = trim(adjustl(grid_files_heirarchy(i)))
      if (localpet == 0 ) write(*,*) trim(adjustl(grid_files_heirarchy(i))), nominal_horizontal_cell_spacing(i), weights(i)
   enddo
   close(14)

   ! Quick sanity check
   if ( i .ne. nmeshes ) then
      call error_handler("Weird error reading GRID_INFO file", -111)
   endif

   ! Make sure all the input files exist
   do i = 1,nmeshes
      inquire(file=grid_files_heirarchy(i),exist=fexist)
      if ( .not. fexist ) then
         deallocate(grid_files_heirarchy,nominal_horizontal_cell_spacing,weights) ! clean-up
         call error_handler("File "//grid_files_heirarchy(i)//" is missing",-112) ! exit
      endif
   enddo

   if ( localpet == 0 ) write(*,*)' Done with subroutine read_grid_info_file'

end subroutine read_grid_info_file

subroutine para_range(n1, n2, nprocs, irank, ista, iend)

  integer, intent(in) :: n1, n2, nprocs, irank
  integer, intent(out) :: ista, iend

  integer :: iwork1, iwork2

  iwork1 = (n2 - n1 + 1) / nprocs
  iwork2 = mod(n2 - n1 + 1, nprocs)
  ista = irank * iwork1 + n1 + min(irank, iwork2)
  iend = ista + iwork1 - 1
  if (iwork2 > irank) iend = iend + 1
  return
end subroutine para_range

subroutine find_interior_cell(localpet,mpas_mesh,interior_cell)
   integer, intent(in)              :: localpet
   type(mpas_mesh_type), intent(in) :: mpas_mesh
   integer, intent(inout)           :: interior_cell

   integer :: iCell, i, j, k, jCell
   logical :: all_inside, all_inside2

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
      if ( all_inside ) then ! now check one more ring, just to make extra sure the point is in the interior
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
   integer, intent(inout), dimension(nbdyCells) :: bdyCells

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

subroutine find_mesh_boundary(localpet,interior_cell, nbdyCells, bdyCells, mark_neighbors, mesh_source, mesh_target)

   integer, intent(in)                       :: localpet,interior_cell, nbdyCells
   logical, intent(in)                       :: mark_neighbors
   integer, dimension(nbdyCells), intent(in) :: bdyCells
   type(mpas_mesh_type), intent(in)          :: mesh_source
   type(mpas_mesh_type), intent(inout)       :: mesh_target

   integer :: start_cell, i, iCell, nearest_cell_to_source_point, nc, j, jCell
   real    :: d

   mesh_target%bdyMaskCell(:) = UNMARKED ! these have been filled in define_grid_mpas...but re-initialize them to UNMARKED

   start_cell = 1
   do i = 1,nbdyCells
      iCell = bdyCells(i)
      nearest_cell_to_source_point = nearest_cell( mesh_source%latCell(iCell), mesh_source%lonCell(iCell), &
                                                   start_cell, mesh_target%nCells, mesh_target%maxEdges, &
                                                   mesh_target%nEdgesOnCell, mesh_target%cellsOnCell , &
                                                   mesh_target%latCell, mesh_target%lonCell)
      nc = nearest_cell_to_source_point ! shorter variable name

      ! make sure distances bewteen the current and newly-found points are close enough?
     !d = sphere_distance(mesh_source%latCell(iCell), mesh_source%lonCell(iCell), mesh_target%latCell(nc), &
      !                                    mesh_target%lonCell(nc), mesh_source%mpas_mesh_radius )
      ! check the distance between the source lat/lon and target lat/lon
     !if ( d .gt. latlon_tolerance_threshold ) then
     !   call error_handler('oops', 175)
     !endif

      mesh_target%bdyMaskCell(nc) = MARKED ! this is the location of a boundary cell on the target mesh. mark it.

      ! also mark the neighbors 
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

end module model_grid
