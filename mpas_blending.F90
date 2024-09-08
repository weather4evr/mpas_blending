! Blend two MPAS fields together on the unstructured mesh.
! Code by Craig Schwartz, based in large part on the "MPASSIT" code by Larissa Reames.
program mpas_blending

use mpi
use esmf
use utils_mod, only      : error_handler, nearest_cell

use program_setup, only  : LogType, read_setup_namelist,  &
                           large_scale_file, small_scale_file, grid_info_file, &
                           average_upscale_before_interp, &
                           output_intermediate_files_up, output_intermediate_files_down, &
                           interp_method, output_latlon_grid, output_blended_filename, &
                           dx_in_degrees, nvars_to_blend, smooth_going_downscale, smoother_dimensionless_coefficient, &
                           max_boundary_layer_to_keep, regional_masking_method
use model_grid, only     : define_grid_mpas, define_grid_latlon, mpas_mesh_type, mpas_meshes, &
                           nmeshes, meshes, grid_files_heirarchy, read_grid_info_file, weights, &
                           latlon_grid, cleanup, cleanup_mpas_mesh_type, nominal_horizontal_cell_spacing, &
                           find_interior_cell, extract_boundary_cells, find_mesh_boundary, mask_value, ignore_value
                          !mask_info_type , update_mesh_mask
use input_data, only     : read_input_data, set_nVertLevelsPerVariable, &
                           input_data_error_checks, nVertLevelsPerVariable
use interp, only         : make_rh, destroy_rh, interp_data, radially_average, apply_smoothing_filter, force_masked_data_to_missing, &
                           extrapolate, &
                           force_masked_data_to_value, force_certain_data_to_value, update_mask, extrapolate_nearest_neigh_on_impacted_points ! should in be bundles_mod
use bundles_mod, only    : add_subtract_bundles, cleanup_bundle, define_bundle, &
                           large_scale_data_going_up, small_scale_data_going_up, &
                           large_scale_data_going_down, small_scale_data_going_down, &
                           large_scale_data_perts, small_scale_data_perts, &
                           large_scale_data_going_up_avg, small_scale_data_going_up_avg, &
                           latlon_bundle, blending_bundle, &
                           large_scale_data_going_up_nn, small_scale_data_going_up_nn
use write_data, only     : write_to_file, write_to_file_latlon

implicit none

integer                   :: ierr, localpet, npets, i, j, interior_cell, my_interior_cell
integer                   :: nbdyCells !, cell_start, cell_end
integer, allocatable      :: bdyCells(:), dstStatus(:)
real                      :: w1, w2
character(len=500)        :: my_output_name
character(len=5)          :: cell_dx, cell_degrees, iii
type(esmf_vm)             :: vm
type(esmf_fieldbundle), allocatable :: tmp_bundle(:)
type(ESMF_RouteHandle)              :: rh, rh_nn
type(ESMF_RouteHandle), allocatable :: rh_latlon(:)
integer(ESMF_KIND_I4), pointer :: unmappedDstList(:)

!------------------
! Initialize mpi
!------------------
call mpi_init(ierr)

!----------------------
! Initialize ESMF stuff
!----------------------
!  print*,"- INITIALIZE ESMF"
call ESMF_Initialize(rc=ierr, logkindflag=LogType)
if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
    call error_handler("INITIALIZING ESMF", ierr)

!print*,"- CALL VMGetGlobal"
call ESMF_VMGetGlobal(vm, rc=ierr)
if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGetGlobal", ierr)

!print*,"- CALL VMGet"
! localpet is the process number; npets is the total number of processors
call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)
if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGet", ierr)

if (localpet==0) print*,'- NPETS IS  ',npets
!print*,'- LOCAL PET ',localpet

!--------------------------------
! All processors read the namelist
!--------------------------------
call read_setup_namelist('input.nml',localpet)

!--------------------------------------------------------------------------------------------
! Get the mesh information for the hierarchy.  meshes(:) has the complete hierarchy, where
!  the ith element of meshes() is the mesh for the ith entry in grid_info_file
! The mesh is decomposed across the various processors, so each processor only has a
!  part of the mesh.
!--------------------------------------------
call read_grid_info_file(localpet,grid_info_file) ! sets grid_files_heirarchy, nmeshes
allocate(meshes(nmeshes))      ! esmf_mesh type
allocate(mpas_meshes(nmeshes)) ! mpas_mesh_type type defined in model_grid
do i = 1,nmeshes
   call define_grid_mpas(localpet, npets, grid_files_heirarchy(i), meshes(i), mpas_meshes(i)) ! output is meshes(i), mpas_meshes(i)
enddo

!-----------------------------------------------------------------------------------------
! Error checking to make sure large_scale_file and small_scale_file have the same meshes.
! Also ensure that input data has the same dimensions as the first file in grid_info_file.
! If there's an error, program will abort.
!-----------------------------------------------------------------------------------------
call input_data_error_checks(large_scale_file,small_scale_file)
call input_data_error_checks(large_scale_file,grid_files_heirarchy(1))

! Get the number of vertical levels for each variable that we are blending
allocate(nVertLevelsPerVariable(nvars_to_blend)) ! from input_data module
call set_nVertLevelsPerVariable(localpet,small_scale_file) ! sets nVertLevelsPerVariable and other vertical dimensions

!-----------------------------------------------------------------------------------------
! Define the bundles for i = 1,nmeshes, but don't associate any data with them yet.
! The ith bundle will be associated with the ith mesh.
! Data for the 1st mesh gets populated by read_input_data (below), and data for the
!  2nd to ... Nth meshes will get filled in calls to interp_data.
! More broadly, define_bundle can be used to define any bundle and 
!  associate that bundle with any "esmf_mesh type" mesh.
!-----------------------------------------------------------------------------------------
allocate(large_scale_data_going_up(nmeshes)) ! esmf_fieldbundle types
allocate(small_scale_data_going_up(nmeshes))
allocate(large_scale_data_going_up_nn(nmeshes))
allocate(small_scale_data_going_up_nn(nmeshes))
do i = 1,nmeshes
   call define_bundle(localpet,meshes(i),large_scale_data_going_up(i))
   call define_bundle(localpet,meshes(i),small_scale_data_going_up(i))
   call define_bundle(localpet,meshes(i),large_scale_data_going_up_nn(i))
   call define_bundle(localpet,meshes(i),small_scale_data_going_up_nn(i))
enddo

!----------------------------------------------
! Read data from input files
!  output is e.g., large_scale_data_going_up(1)
!----------------------------------------------
call read_input_data(localpet,meshes(1),mpas_meshes(1),large_scale_file,large_scale_data_going_up(1)) ! data from file providing large scales
call read_input_data(localpet,meshes(1),mpas_meshes(1),small_scale_file,small_scale_data_going_up(1)) ! data from file providing small scales

! If our input data is regional, make an initial route-handle where we interpolate from the finest mesh (i=1)
! to all other meshes in the heirarchy, and find those cells where the interpolation won't work
! we then mask out all of those cells in the coarser resolution meshses and they are not 
! active in interpolation
if ( mpas_meshes(1)%regional_mesh .and. trim(adjustl(regional_masking_method)) .ne. "none" ) then
   if ( trim(adjustl(regional_masking_method)) == "mpas_bdy" ) then
      call find_interior_cell(localpet, mpas_meshes(1), interior_cell) ! return a global cellID (across whole mesh) we are pretty sure is in the interior
      nbdyCells = count(mpas_meshes(1)%bdyMaskCell == max_boundary_layer_to_keep)
      allocate(bdyCells(nbdyCells))
      if (localpet==0) write(*,fmt='(a26,i5,a3,i5)')' mesh 1 num bdyMaskCell = ', max_boundary_layer_to_keep, ' : ',nbdyCells
      call extract_boundary_cells(mpas_meshes(1), max_boundary_layer_to_keep, nbdyCells, bdyCells) ! get the cells on the mesh that make up the boundary
   endif
   if (localpet==0) write(*,fmt='(a31,i10)')' num masked cells for mesh    1 = ', count(mpas_meshes(1)%bdyMaskCell == mask_value)
   do i = 1,nmeshes-1
      if ( trim(adjustl(regional_masking_method)) == "esmf" ) then
         allocate(dstStatus(mpas_meshes(i+1)%nCells))
         ! make the ESMF routehandle for interpolation between the ith and (i+1)th mesh.
         ! Main input is an esmf_fieldbundle on the ith (source) and (i+1)th (target) mesh. rh is the output (that we don't use...could not compute).
         ! Also, optionally return a field that has the status of each destination point on the full grid (dstStatus).
         call make_rh(localpet,large_scale_data_going_up(1), large_scale_data_going_up(i+1),trim(adjustl(interp_method)),rh,unmappedDstList,dstStatus)
         !See 54.50 ESMF_REGRIDSTATUS for meaning of ESMF_REGRIDSTATUS_MAPPED
         where (dstStatus < ESMF_REGRIDSTATUS_MAPPED ) mpas_meshes(i+1)%bdyMaskCell = mask_value
         call destroy_rh(rh)
         if (localpet == 0 ) then
            do j = 0,8
               write(*,*)'num status = ',j,' = ',count(dstStatus == j)
            enddo
         endif
         deallocate(dstStatus)
      else if ( trim(adjustl(regional_masking_method)) == "mpas_bdy" ) then
         ! find the boundary on each mesh, based on the finest mesh (mpas_meshes(1))
         ! then fill the region within the boundary, and set area outside of region to mask_value

         ! find the point on the (i+1)th mesh corresponding to an interior_cell on the finest resolution mesh
         my_interior_cell = nearest_cell( mpas_meshes(1)%latCell(interior_cell), mpas_meshes(1)%lonCell(interior_cell), &
                                                      1, mpas_meshes(i+1)%nCells, mpas_meshes(i+1)%maxEdges, &
                                                      mpas_meshes(i+1)%nEdgesOnCell, mpas_meshes(i+1)%cellsOnCell , &
                                                      mpas_meshes(i+1)%latCell, mpas_meshes(i+1)%lonCell)
         if (localpet==0) then
            write(*,fmt='(a32,i4,a3,2f12.4)')'interior point lat/lon for mesh ',(i+1),' = ',mpas_meshes(i+1)%latCell(my_interior_cell)*180.0/3.14, &
                                                                                            mpas_meshes(i+1)%lonCell(my_interior_cell)*180.0/3.14
         endif

         ! 5th entry .false. means don't also "mark" neighbors of the boundary
         call find_mesh_boundary(localpet,my_interior_cell, nbdyCells, bdyCells, .false., mpas_meshes(1), mpas_meshes(i+1)) ! updates mpas_meshes(i+1)%bdyMaskCell (last argument)

         ! it's possible that there could be a hole in the boundary on one of the meshes
         ! if so, no cells will be masked
         ! fill in the hole by marking the neighbors of points supposedly on the boundary
         if ( count(mpas_meshes(i+1)%bdyMaskCell == mask_value) == 0 ) then
            ! 5th entry .true. means DO "mark" neighbors of the boundary
            if (localpet==0) write(*,fmt='(a,i4)')' filling a hole for mesh ',(i+1)
            call find_mesh_boundary(localpet,my_interior_cell, nbdyCells, bdyCells, .true., mpas_meshes(1), mpas_meshes(i+1)) ! updates mpas_meshes(i+1)%bdyMaskCell (last argument)
         endif
      endif
      ! mpas_meshes(i+1)%bdyMaskCell = mask_value at areas outside of the region determined by max_boundary_layer_to_keep
      ! now update the mask on the esmf_mesh using mpas_meshes(i+1)%bdyMaskCell
     !call ESMF_MeshTurnOnCellMask(meshes(i+1), mask_value, rc=ierr)
     !call update_mesh_mask(meshes(i+1),mpas_meshes(i+1)) ! probably doesn't do anything b/c elementMasking may be broken
      if (localpet==0) write(*,fmt='(a27,i4,a3,i10)')' num masked cells for mesh ',(i+1),' = ', count(mpas_meshes(i+1)%bdyMaskCell == mask_value)
   enddo ! loop over the ith mesh
   if ( trim(adjustl(regional_masking_method)) == "mpas_bdy" ) deallocate(bdyCells)
endif ! is mpas_meshes(1)%regional_mesh == .true.

!------------------------------------------------------------------------------------------------------
! Upscale the fields through the heirarchy for both large_scale_file and small_scale_file.
! Upscale from the ith grid to the (i+1)th grid.  Both input and output are esmf_fieldbundle type.
! large_scale_data_going_up(1) and small_scale_data_going_up(1) were filled in call to read_input_data
! Interpolation is a 3 step process:
! 1) Make the routehandle.  Helpful to make it outside of the interpolation routine in case we want to
!    use it repeatedly.
! 2) Use the routehandle to do the interpolation.
! 3) Destroy the routehandle once we're done with it.
!------------------------------------------------------------------------------------------------------
if ( average_upscale_before_interp ) then
   allocate(large_scale_data_going_up_avg(nmeshes-1)) ! esmf_fieldbundle types
   allocate(small_scale_data_going_up_avg(nmeshes-1))
endif
do i = 1,nmeshes-1
   ! make the ESMF routehandle for interpolation between the ith and (i+1)th mesh.
   ! Main input is an esmf_fieldbundle on the ith (source) and (i+1)th (target) mesh. rh is the output.
   allocate(dstStatus(mpas_meshes(i+1)%nCells))
   call make_rh(localpet,large_scale_data_going_up(i), large_scale_data_going_up(i+1),trim(adjustl(interp_method)),rh, unmappedDstList, dstStatus)
  !call make_rh(localpet,large_scale_data_going_up(i), large_scale_data_going_up(i+1),'nearest',rh_nn, unmappedDstList)

  ! If requested, average data on the current mesh to the scale of the next mesh, which is given by
  !   nominal_horizontal_cell_spacing(i+1)
  ! output is large_scale_data_going_up_avg(i), small_scale_data_going_up_avg(i), which is the spatially-averaged
  ! field on the ith mesh
   if ( average_upscale_before_interp ) then
      call define_bundle(localpet,meshes(i),large_scale_data_going_up_avg(i))
      call define_bundle(localpet,meshes(i),small_scale_data_going_up_avg(i))
      call radially_average(localpet,mpas_meshes(i),nominal_horizontal_cell_spacing(i+1),large_scale_data_going_up(i),large_scale_data_going_up_avg(i))
      call radially_average(localpet,mpas_meshes(i),nominal_horizontal_cell_spacing(i+1),small_scale_data_going_up(i),small_scale_data_going_up_avg(i))
      call interp_data(localpet, rh, large_scale_data_going_up_avg(i), large_scale_data_going_up(i+1))
      call interp_data(localpet, rh, small_scale_data_going_up_avg(i), small_scale_data_going_up(i+1))
     !call interp_data(localpet, rh_nn, large_scale_data_going_up_avg(i), large_scale_data_going_up_nn(i+1)) ! nearest-neighbor interp
     !call interp_data(localpet, rh_nn, small_scale_data_going_up_avg(i), small_scale_data_going_up_nn(i+1))
   else 
      call interp_data(localpet, rh, large_scale_data_going_up(i), large_scale_data_going_up(i+1))
      call interp_data(localpet, rh, small_scale_data_going_up(i), small_scale_data_going_up(i+1))
     !call interp_data(localpet, rh_nn, large_scale_data_going_up(i), large_scale_data_going_up_nn(i+1)) ! nearest-neighbor interp
     !call interp_data(localpet, rh_nn, small_scale_data_going_up(i), small_scale_data_going_up_nn(i+1))
   endif
   call destroy_rh(rh)
  !call destroy_rh(rh_nn)

   if ( mpas_meshes(1)%regional_mesh ) then
      ! The challenge with interpolating through a heirarchy is dealing with ESMF's missing values after an interpolation.
      ! First, if a point cannot be interpolated to, ESMF sets the value to 0 on the target grid, which is annoying, b/c 0 is a physical
      !  value.  So, force those points that can't be interpolated to to mask_value.
      ! Second, there is the issue where certain points will be influenced by a masked-out point, but the points themselves aren't 
      !  masked out. These points are identifiable, as they will be somehow related to a mask-out point, as defined by mask_value, 
      !  and will be a large negative number.

      ! output of this next subroutine is e.g., updated large_scale_data_going_up(i+1) [third argument]
      !    and mpas_mesh%bdyMaskCell (second argument) if fourth argument is true
      ! set last argument to .false. in first call so as not to mess-up input to second call
     !call force_masked_data_to_missing(localpet,mpas_meshes(i+1), large_scale_data_going_up(i+1), .false. , mpas_meshes(i), large_scale_data_going_up(i) )
     !if (localpet==0) write(*,fmt='(i4,a24,i4,a3,i10)')localpet,' num bdy=missing before ',(i+1),' = ', count(mpas_meshes(i+1)%bdyMaskCell == mask_value)
     !call force_masked_data_to_missing(localpet,mpas_meshes(i+1), small_scale_data_going_up(i+1), .true. , mpas_meshes(i),  small_scale_data_going_up(i) )
     !if (localpet==0) write(*,fmt='(i4,a23,i4,a3,i10)')localpet,' num bdy=missing after ',(i+1),' = ', count(mpas_meshes(i+1)%bdyMaskCell == mask_value)

      ! Where the field is masked, set to mask_value. If ESMF can't interpolate to a point, it will give values of 0,
      !  and most of these will be where the field is masked
      ! When we create the masks, we are looking at how the fields on i = 1 (the finest mesh) map to all other meshes.
      ! But, the interpolation is from mesh i to (i+1), and we haven't dealt with masking in that case.
      ! When we made the routehandle, we returned dstStatus, which is the status of each interpolated point.  If dstStatus
      !  is certain values, that indicates the interpolation was not successful.  So, set those points to a masked value.
      !See 54.50 ESMF_REGRIDSTATUS for meaning of ESMF_REGRIDSTATUS_MAPPED
      where (dstStatus < ESMF_REGRIDSTATUS_MAPPED ) mpas_meshes(i+1)%bdyMaskCell = mask_value ! update the mask...useful for averaging

      ! mpas_meshes(i+1)%bdyMaskCell now reflects the initial masks (from i = 1 to all other meshes), as well as from
      ! interpolating from the i to i+1)th mesh.
      ! This next call will ensure that where the field is masked, the fields in the bundle will be set to mask_val. 
     !call force_masked_data_to_value(localpet,mpas_meshes(i+1),mpas_meshes(i+1)%bdyMaskCell,mask_value,mask_value,large_scale_data_going_up(i+1))
     !call force_masked_data_to_value(localpet,mpas_meshes(i+1),mpas_meshes(i+1)%bdyMaskCell,mask_value,mask_value,small_scale_data_going_up(i+1))

      if ( 1 == 1 ) then ! extrapolate
         ! where points in the bundle are < ignore_value, but NOT equal to mask_value, use a nearest neighbor approach to fill it
        !call extrapolate_nearest_neigh_on_impacted_points(localpet,mpas_meshes(i+1),large_scale_data_going_up(i+1),mpas_meshes(i), large_scale_data_going_up(i))
       !!call extrapolate_nearest_neigh_on_impacted_points(localpet,mpas_meshes(i+1),small_scale_data_going_up(i+1),mpas_meshes(i), small_scale_data_going_up(i))
       !!call extrapolate_nearest_neigh_on_impacted_points(localpet,mpas_meshes(i+1),large_scale_data_going_up(i+1),mpas_meshes(i), large_scale_data_going_up_nn(i+1))
        !call extrapolate_nearest_neigh_on_impacted_points(localpet,mpas_meshes(i+1),small_scale_data_going_up(i+1),mpas_meshes(i), small_scale_data_going_up_nn(i+1))
         call extrapolate(localpet,mpas_meshes(i+1),large_scale_data_going_up(i+1),mpas_meshes(i),large_scale_data_going_up(i),unmappedDstList)
         call extrapolate(localpet,mpas_meshes(i+1),small_scale_data_going_up(i+1),mpas_meshes(i),small_scale_data_going_up(i),unmappedDstList)
      else
         ! where the bundle is < ignore_value, set the bundle to mask_value
        !call force_certain_data_to_value(localpet,large_scale_data_going_up(i+1),'lt',ignore_value,mask_value)
        !call force_certain_data_to_value(localpet,small_scale_data_going_up(i+1),'lt',ignore_value,mask_value)
      endif

      ! where a field in the bundle large_scale_data_going_up(i+1) < ignore_value, set mpas_meshes(i+1)%bdyMaskCell = mask_value
      ! masks are the same for all meshes, so just need to do for one field/bundle
     !call update_mask(localpet,mpas_meshes(i+1),large_scale_data_going_up(i+1),ignore_value, mpas_meshes(i+1)%bdyMaskCell, mask_value)

      ! some non-masked points may have been influenced by masked points.
      ! force those influenced points to missing, and then reset the elementMask in the mesh
    !!call ESMF_MeshTurnOnCellMask(meshes(i+1), mask_value, rc=ierr)
    !!call update_mesh_mask(meshes(i+1),mpas_meshes(i+1)) ! probably doesn't do anything b/c elementMasking may be broken
   endif

   deallocate(dstStatus)

   call mpi_barrier(mpi_comm_world, ierr) ! sync-up before going to next mesh

enddo

!-----------------------------------------
! Output the upscaled fields if requested
!-----------------------------------------
if ( output_intermediate_files_up ) then
   ! Output one grid per file, because all meshes are different.
   ! In write_to_file, .false. means define the output mesh from input data and the next entry is ignored.
   !  If it's .true., then the next entry provides a template for the output file
   do i = 1, nmeshes
      write(cell_dx,fmt='(f5.1)') nominal_horizontal_cell_spacing(i)
      my_output_name = 'mpas_mesh_largeScaleData_goingUp_'//trim(adjustl(cell_dx))//'km.nc'
      call write_to_file(localpet, mpas_meshes(i), large_scale_data_going_up(i), .false., large_scale_file, my_output_name)
      my_output_name = 'mpas_mesh_smallScaleData_goingUp_'//trim(adjustl(cell_dx))//'km.nc'
      call write_to_file(localpet, mpas_meshes(i), small_scale_data_going_up(i), .false., large_scale_file, my_output_name)
      if ( average_upscale_before_interp .and. ( i .lt. nmeshes) ) then ! spatially smoothed/averaged fields only 1:nmeshes-1
         my_output_name = 'mpas_mesh_largeScaleData_goingUp_avg_'//trim(adjustl(cell_dx))//'km.nc'
         call write_to_file(localpet, mpas_meshes(i), large_scale_data_going_up_avg(i), .false., large_scale_file, my_output_name)
        !my_output_name = 'mpas_mesh_smallScaleData_goingUp_avg'//trim(adjustl(cell_dx))//'km.nc'
        !call write_to_file(localpet, mpas_meshes(i), small_scale_data_going_up_avg(i), .false., large_scale_file, my_output_name)
      endif
   enddo

   ! Output 1 file containing all the meshes interpolated to the same lat-lon grid
   if ( output_latlon_grid ) then
      write(cell_degrees,fmt='(f5.3)') dx_in_degrees

      ! define the lat-lon grid and latlon_bundle
      call define_grid_latlon(localpet,npets) ! defines latlon_grid in model_grid
      allocate(latlon_bundle(nmeshes))  ! esmf_fieldbundle type
      do i = 1,nmeshes
         call define_bundle(localpet,latlon_grid,latlon_bundle(i)) ! output is latlon_bundle(i)
      enddo

      ! First output data from file providing large scales. Need to interpolate the data onto lat-lon grid first.
      ! The same routehandle can be used for the large- and small-scale files.
      allocate(rh_latlon(nmeshes))
      do i = 1,nmeshes
         call make_rh(localpet,large_scale_data_going_up(i), latlon_bundle(i), 'bilinear', rh_latlon(i), unmappedDstList )
         call interp_data(localpet, rh_latlon(i), large_scale_data_going_up(i), latlon_bundle(i))
      enddo
      my_output_name = 'latlon_mesh_largeScaleData_goingUp_'//trim(adjustl(cell_degrees))//'degrees.nc'
      call write_to_file_latlon(localpet,nmeshes,latlon_bundle, my_output_name)

      ! Then output data from file providing small scales. Need to interpolate the data onto lat-lon grid first.
      do i = 1,nmeshes
         call interp_data(localpet, rh_latlon(i), small_scale_data_going_up(i), latlon_bundle(i))
         call destroy_rh(rh_latlon(i))
      enddo
      deallocate(rh_latlon)
      my_output_name = 'latlon_mesh_smallScaleData_goingUp_'//trim(adjustl(cell_degrees))//'degrees.nc'
      call write_to_file_latlon(localpet,nmeshes,latlon_bundle, my_output_name)
   endif
endif

! Given that processes might have been writing, sync up before printing
call mpi_barrier(mpi_comm_world, ierr)
if ( localpet == 0 ) write(*,*)'Done with upscaling and outputting optional files' 

!--------------------------
! Now downscale the fields
!--------------------------
! Downscale the fields from the (i+1)th mesh to the ith mesh. We need this to compute perturbations.
!  The ith entry is for the data on the ith mesh that has been downscaled from the (i+1)th mesh
allocate(large_scale_data_going_down(nmeshes)) ! esmf_fieldbundle types
allocate(small_scale_data_going_down(nmeshes))
do i = 1,nmeshes-1
   call define_bundle(localpet,meshes(i),large_scale_data_going_down(i))
   call define_bundle(localpet,meshes(i),small_scale_data_going_down(i))

   call make_rh(localpet, large_scale_data_going_up(i+1), large_scale_data_going_down(i), trim(adjustl(interp_method)), rh, unmappedDstList)
   call interp_data(localpet, rh, large_scale_data_going_up(i+1), large_scale_data_going_down(i))
   call interp_data(localpet, rh, small_scale_data_going_up(i+1), small_scale_data_going_down(i))
   call destroy_rh(rh)

   if ( mpas_meshes(1)%regional_mesh ) then
      ! output is e.g., updated large_scale_data_going_down(i) [third argument]
     !call force_masked_data_to_missing(localpet,mpas_meshes(i), large_scale_data_going_down(i), .false. )
     !call force_masked_data_to_missing(localpet,mpas_meshes(i), small_scale_data_going_down(i), .false. )

      ! force 0s from the interpolation to mask_value
     !call force_masked_data_to_value(localpet,mpas_meshes(i),mpas_meshes(i)%bdyMaskCell,mask_value,mask_value,large_scale_data_going_down(i))
     !call force_masked_data_to_value(localpet,mpas_meshes(i),mpas_meshes(i)%bdyMaskCell,mask_value,mask_value,small_scale_data_going_down(i))

      ! where the bundle is < ignore_value, set to mask_value
     !call force_certain_data_to_value(localpet,large_scale_data_going_down(i),'lt',ignore_value,mask_value)
     !call force_certain_data_to_value(localpet,small_scale_data_going_down(i),'lt',ignore_value,mask_value)
   endif

enddo
! fill in the last (coarsest) element from the upscaled data
large_scale_data_going_down(nmeshes) = large_scale_data_going_up(nmeshes) ! LHS is basically a pointer to RHS
small_scale_data_going_down(nmeshes) = small_scale_data_going_up(nmeshes)

!--------------------------
! Compute perturbations
!--------------------------
!function add_subtract_bundles(localpet,operation,bundle1,bundle2,w1,w2) result(bundle_out)
!  add_subtract_bundles either adds or subracts based on value of operation (2nd entry):
!  if 'subtract', it computes bundle1 minus bundle2 (third argument minus fourth argument)
!  if 'add',      it does w1*bundle1 + w2*bundle2
!  if 'subtract', w1 and w2 are ignored
allocate(large_scale_data_perts(nmeshes-1)) ! esmf_fieldbundle types
allocate(small_scale_data_perts(nmeshes-1))
do i = 1, nmeshes-1
   call define_bundle(localpet,meshes(i),large_scale_data_perts(i))
   call define_bundle(localpet,meshes(i),small_scale_data_perts(i))
   large_scale_data_perts(i) = add_subtract_bundles(localpet, 'subtract', &
         large_scale_data_going_up(i), large_scale_data_going_down(i), 1.0, 1.0)
   small_scale_data_perts(i) = add_subtract_bundles(localpet, 'subtract', &
         small_scale_data_going_up(i), small_scale_data_going_down(i), 1.0, 1.0)
   if ( mpas_meshes(1)%regional_mesh ) then
     !call force_masked_data_to_missing(localpet,mpas_meshes(i), large_scale_data_perts(i), .false. )
     !call force_masked_data_to_missing(localpet,mpas_meshes(i), small_scale_data_perts(i), .false. )

      ! force perturbations in masked area to 0
     !call force_masked_data_to_value(localpet,mpas_meshes(i),mpas_meshes(i)%bdyMaskCell,mask_value,mask_value,large_scale_data_perts(i))
     !call force_masked_data_to_value(localpet,mpas_meshes(i),mpas_meshes(i)%bdyMaskCell,mask_value,mask_value,small_scale_data_perts(i))

      ! where the bundle is < ignore_value, set to 0
     !call force_certain_data_to_value(localpet,large_scale_data_perts(i),'lt',ignore_value,mask_value)
     !call force_certain_data_to_value(localpet,small_scale_data_perts(i),'lt',ignore_value,mask_value)

      ! where the bundle is > -1*mask_value, set to 0
     !call force_certain_data_to_value(localpet,large_scale_data_perts(i),'gt',abs(mask_value+1),mask_value)
     !call force_certain_data_to_value(localpet,small_scale_data_perts(i),'gt',abs(mask_value+1),mask_value)
   endif
enddo

!------------------
! Do the blending
!------------------
! We need to reconstruct a field on the finest-resolution mesh (i=1) from
!  coarser meshes, starting from the coarsest mesh, and working our way down.
! To do this, we successively add perturbations, which are weighted between
!  the two datasets per blending specifications.
call mpi_barrier(mpi_comm_world, ierr)
if ( localpet == 0 ) write(*,*)'About to blend'

allocate(blending_bundle(nmeshes))
allocate(tmp_bundle(nmeshes)) ! for some reason, this can't just be tmp_bundle(1)....
do i = 1,nmeshes
   call define_bundle(localpet,meshes(i),blending_bundle(i))
   call define_bundle(localpet,meshes(i),tmp_bundle(i))
enddo

! First, blend on the coarsest mesh. This is just a weighted sum.
i = nmeshes
w1 = real(weights(i)) ! weights to data providing large-scales
w2 = real(1.0 - w1)

!add_subtract_bundles(localpet,operation,bundle1,bundle2,w1,w2)
!  if operation == 'add',      it does w1*bundle1 + w2*bundle2
blending_bundle(i) = add_subtract_bundles(localpet,'add', large_scale_data_going_down(i), &
                                          small_scale_data_going_down(i), w1, w2)

! Now start adding perturbations to the coarsest field to go downscale:
!  1) Interpolate the blended field from the ith mesh to the (i-1)th
!  2) Add the perturbations of the two datasets together with
!  the specified weights.
!  3) Add the blended perturbations to the field on the (i-1)th mesh
!       (step 1)
!  4) Iterate; blending_bundle(1) is the reconstructed/blended field
!     on the finest-resolution mesh we want to output
do i = nmeshes, 2, -1 ! loop from high to low
   if ( localpet == 0 ) write(*,*)'Blending for mesh ',i
   call make_rh(localpet, blending_bundle(i),blending_bundle(i-1),trim(adjustl(interp_method)), rh, unmappedDstList)
   call interp_data(localpet,rh,blending_bundle(i),blending_bundle(i-1))
   call destroy_rh(rh)
   if ( mpas_meshes(1)%regional_mesh ) then
     !call force_masked_data_to_missing(localpet,mpas_meshes(i-1), blending_bundle(i-1), .false. )

      ! force 0s from the interpolation to mask_value
     !call force_masked_data_to_value(localpet,mpas_meshes(i-1),mpas_meshes(i-1)%bdyMaskCell,mask_value,mask_value,blending_bundle(i-1))

      ! where the bundle is < ignore_value, set to mask_value
     !call force_certain_data_to_value(localpet,blending_bundle(i-1),'lt',ignore_value,mask_value)
   endif
   if ( smooth_going_downscale ) then
      call apply_smoothing_filter(localpet,mpas_meshes(i-1),smoother_dimensionless_coefficient,blending_bundle(i-1))
   endif

   w1 = real(weights(i-1)) ! weights to data providing large-scales
   w2 = real(1.0 - w1)
   tmp_bundle(i-1) = add_subtract_bundles(localpet, 'add', &
                       large_scale_data_perts(i-1),small_scale_data_perts(i-1), w1, w2)

   blending_bundle(i-1) = add_subtract_bundles(localpet, 'add', &
                    blending_bundle(i-1),tmp_bundle(i-1), 1.0, 1.0)

enddo

!-----------------------------------------
! Output the downscaled fields if requested
!-----------------------------------------
if ( output_intermediate_files_down ) then
   do i = 1, nmeshes-1
      write(cell_dx,fmt='(f5.1)') nominal_horizontal_cell_spacing(i)
      my_output_name = 'mpas_mesh_largeScalePerts_'//trim(adjustl(cell_dx))//'km.nc'
      call write_to_file(localpet, mpas_meshes(i), large_scale_data_perts(i), .false., large_scale_file, my_output_name)
   enddo
   do i = 2, nmeshes
      write(cell_dx,fmt='(f5.1)') nominal_horizontal_cell_spacing(i)
      my_output_name = 'mpas_mesh_blendingBundleStep_'//trim(adjustl(cell_dx))//'km.nc'
      call write_to_file(localpet, mpas_meshes(i), blending_bundle(i), .false., large_scale_file, my_output_name)
   enddo
   if ( output_latlon_grid ) then
      write(cell_degrees,fmt='(f5.3)') dx_in_degrees
      do i = 1,nmeshes
         call make_rh(localpet, small_scale_data_going_down(i), latlon_bundle(i), 'bilinear', rh, unmappedDstList)
         call interp_data(localpet, rh, small_scale_data_going_down(i), latlon_bundle(i))
         call destroy_rh(rh)
      enddo
      my_output_name = 'latlon_mesh_smallScaleData_goingDown_'//trim(adjustl(cell_degrees))//'degrees.nc'
      call write_to_file_latlon(localpet,nmeshes-1,latlon_bundle, my_output_name)

      do i = 1,nmeshes-1
         call make_rh(localpet, large_scale_data_perts(i), latlon_bundle(i), 'bilinear', rh, unmappedDstList)
         call interp_data(localpet, rh, large_scale_data_perts(i), latlon_bundle(i))
         call destroy_rh(rh)
      enddo
      my_output_name = 'latlon_mesh_largeScalePerts_'//trim(adjustl(cell_degrees))//'degrees.nc'
      call write_to_file_latlon(localpet,nmeshes-1,latlon_bundle, my_output_name)
   endif
endif
 
!--------------------------
! Write blended data to file
!--------------------------
! use small_scale_file as template for output file (because 4th argument is .true.)
call write_to_file(localpet, mpas_meshes(1), blending_bundle(1), .true., small_scale_file, output_blended_filename)

call mpi_barrier(mpi_comm_world, ierr)

if ( output_latlon_grid ) then
   call make_rh(localpet, blending_bundle(1), latlon_bundle(1), 'bilinear', rh, unmappedDstList)
   call interp_data(localpet, rh, blending_bundle(1), latlon_bundle(1))
   call destroy_rh(rh)
   write(cell_degrees,fmt='(f5.3)') dx_in_degrees
   my_output_name = 'latlon_mesh_BlendedFields_'//trim(adjustl(cell_degrees))//'degrees.nc'
   call write_to_file_latlon(localpet,1,latlon_bundle(1), my_output_name)
endif
 
!-------------
! Clean up 
!-------------
! degbugging
do i = 1,nmeshes
   !write(iii,fmt='(i5)') i
   call cleanup_bundle(large_scale_data_going_up(i))
   call cleanup_bundle(small_scale_data_going_up(i))
   call cleanup_bundle(blending_bundle(i))
   call cleanup_bundle(tmp_bundle(i))
   call cleanup_bundle(large_scale_data_going_up_nn(i))
   call cleanup_bundle(small_scale_data_going_up_nn(i))
   if ( output_latlon_grid ) call cleanup_bundle(latlon_bundle(i))
   call cleanup_mpas_mesh_type(mpas_meshes(i))
enddo

! degbugging
!call ESMF_LogWrite("Log Write", ESMF_LOGMSG_INFO, rc=ierr)

do i = 1,nmeshes-1
   call cleanup_bundle(large_scale_data_perts(i))
   call cleanup_bundle(small_scale_data_perts(i))
   call cleanup_bundle(large_scale_data_going_down(i)) ! this is length nmeshes, but the entry for "nmeshes" is a pointer to large_scale_data_going_up ...
   call cleanup_bundle(small_scale_data_going_down(i)) !  which is deallocated above, so only need to deallocate these from 1:nmeshes-1
enddo

if ( average_upscale_before_interp ) then
   do i = 1,nmeshes-1
      call cleanup_bundle(large_scale_data_going_up_avg(i))
      call cleanup_bundle(small_scale_data_going_up_avg(i))
   enddo
   deallocate(large_scale_data_going_up_avg, small_scale_data_going_up_avg)
endif

deallocate(large_scale_data_going_up, small_scale_data_going_up)
deallocate(large_scale_data_going_down, small_scale_data_going_down)
deallocate(blending_bundle, tmp_bundle)
if ( output_latlon_grid ) deallocate(latlon_bundle)
deallocate(mpas_meshes)
deallocate(large_scale_data_perts, small_scale_data_perts)
deallocate(nVertLevelsPerVariable)

call cleanup(localpet) ! deallocates variables allocated or defined in model_grid module

if (localpet==0) print*,"- CALL ESMF_finalize"
call ESMF_finalize(endflag=ESMF_END_KEEPMPI, rc=ierr)

call mpi_finalize(ierr)

print*,"- DONE."

end program mpas_blending
