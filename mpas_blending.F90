! Blend two MPAS fields together on the unstructured mesh.
! Code by Craig Schwartz, based in large part on the "MPASSIT" code by Larissa Reames.
program mpas_blending

use mpi
use esmf
use utils_mod, only      : error_handler, print_dstStatus
use program_setup, only  : LogType, read_setup_namelist,  &
                           large_scale_file, small_scale_file, grid_info_file, &
                           average_upscale_before_interp, &
                           output_intermediate_files_up, output_intermediate_files_down, &
                           interp_method, extrap_method, extrap_method_latlon, output_latlon_grid, output_blended_filename, &
                           dx_in_degrees, nvars_to_blend, smooth_going_downscale, smoother_dimensionless_coefficient
use model_grid, only     : define_grid_mpas, define_grid_latlon, mpas_mesh_type, mpas_meshes, &
                           nmeshes, meshes, grid_files_heirarchy, read_grid_info_file, weights, &
                           latlon_grid, cleanup, cleanup_mpas_mesh_type, nominal_horizontal_cell_spacing, &
                           mask_value, ignore_value, regional_mesh
use input_data, only     : read_input_data, set_nVertLevelsPerVariable, &
                           input_data_error_checks, nVertLevelsPerVariable
use interp, only         : make_rh, destroy_rh, interp_data, radially_average, apply_smoothing_filter, &
                           set_global_extrapolation, extrapolate, &
                           force_masked_data_to_value, force_certain_data_to_value, update_mask  ! should in be bundles_mod
use bundles_mod, only    : add_subtract_bundles, cleanup_bundle, define_bundle
use write_data, only     : write_to_file, write_to_file_latlon

implicit none

integer                   :: ierr, localpet, npets, i, j, ii, cell_start, cell_end, nCellsPerPET
integer                   :: num_unmapped, num_unmapped_tot, my_num_unmapped
integer, allocatable      :: dstStatus(:), my_dstStatus(:)
real                      :: w1, w2
character(len=500)        :: my_output_name
character(len=5)          :: cell_dx, cell_degrees, iii
type(esmf_vm)             :: vm
type(esmf_fieldbundle), allocatable :: tmp_bundle(:)
type(ESMF_RouteHandle)              :: rh
type(ESMF_RouteHandle), allocatable :: rh_latlon(:)
integer(ESMF_KIND_I4), pointer :: unmappedDstList(:)

type(esmf_fieldbundle), allocatable :: large_scale_data_going_up(:)
type(esmf_fieldbundle), allocatable :: small_scale_data_going_up(:)
type(esmf_fieldbundle), allocatable :: large_scale_data_going_down(:)
type(esmf_fieldbundle), allocatable :: small_scale_data_going_down(:)
type(esmf_fieldbundle), allocatable :: large_scale_data_perts(:)
type(esmf_fieldbundle), allocatable :: small_scale_data_perts(:)
type(esmf_fieldbundle), allocatable :: latlon_bundle(:)
type(esmf_fieldbundle), allocatable :: blending_bundle(:)
type(esmf_fieldbundle), allocatable :: large_scale_data_going_up_avg(:)
type(esmf_fieldbundle), allocatable :: small_scale_data_going_up_avg(:)

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

if (localpet==0) write(*,*)'- NPETS IS  ',npets

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
regional_mesh = mpas_meshes(1)%regional_mesh ! set the public variable defined in model_grid

!--------------------------------------------
! Check the extrapolation method set in the namelist.
! The global variable 'extrap_method' could be modified
!--------------------------------------------
call set_global_extrapolation(localpet)

!-----------------------------------------------------------------------------------------
! Error checking to make sure large_scale_file and small_scale_file have the same meshes.
! Also ensure that input data has the same dimensions as the first file in grid_info_file.
! If there's an error, program will abort.
!-----------------------------------------------------------------------------------------
call input_data_error_checks(large_scale_file,small_scale_file)
call input_data_error_checks(large_scale_file,grid_files_heirarchy(1))

! Get the number of vertical levels for each variable that we are blending
allocate(nVertLevelsPerVariable(nvars_to_blend)) ! from input_data module
call set_nVertLevelsPerVariable(localpet,small_scale_file) ! sets nVertLevelsPerVariable and vertical dimensions

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
do i = 1,nmeshes
   call define_bundle(localpet,meshes(i),large_scale_data_going_up(i))
   call define_bundle(localpet,meshes(i),small_scale_data_going_up(i))
enddo

!----------------------------------------------
! Read data from input files
!  output is e.g., large_scale_data_going_up(1)
!----------------------------------------------
call read_input_data(localpet,meshes(1),mpas_meshes(1),large_scale_file,large_scale_data_going_up(1)) ! data from file providing large scales
call read_input_data(localpet,meshes(1),mpas_meshes(1),small_scale_file,small_scale_data_going_up(1)) ! data from file providing small scales

!------------------------------------------------------------------------------------------------------
! Upscale the fields through the heirarchy for both large_scale_file and small_scale_file.
! Upscale from the ith grid to the (i+1)th grid.  Both input and output are esmf_fieldbundle type.
! large_scale_data_going_up(1) and small_scale_data_going_up(1) were filled in call to read_input_data
! Interpolation is a 3 step process:
! 1) Make the routehandle.  Helpful to make it outside of the interpolation routine so we can use it repeatedly.
! 2) Use the routehandle to do the interpolation.
! 3) Destroy the routehandle once we're done with it.
!------------------------------------------------------------------------------------------------------
if ( average_upscale_before_interp ) then
   allocate(large_scale_data_going_up_avg(nmeshes-1)) ! esmf_fieldbundle types
   allocate(small_scale_data_going_up_avg(nmeshes-1))
endif
do i = 1,nmeshes-1
   ! make the ESMF routehandle for interpolation between the ith and (i+1)th mesh.
   ! Main input is an esmf_fieldbundle on the ith (source) and (i+1)th (target) mesh. rh is the output. Also return some information
   !  about status of the interpolation at all the points on the target mesh. unmappedDstList is a list of the global MPAS cell IDs
   !  where the target cell couldn't be mapped to the source grid, and is unique for each processor, while dstStatus is a field
   !  on the full mesh
   allocate(dstStatus(mpas_meshes(i+1)%nCells))
   nullify(unmappedDstList)
   call make_rh(localpet,large_scale_data_going_up(i), large_scale_data_going_up(i+1),trim(adjustl(interp_method)),trim(adjustl(extrap_method)), rh, unmappedDstList, dstStatus)

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
   else 
      call interp_data(localpet, rh, large_scale_data_going_up(i), large_scale_data_going_up(i+1))
      call interp_data(localpet, rh, small_scale_data_going_up(i), small_scale_data_going_up(i+1))
   endif
   call destroy_rh(rh)

   if ( regional_mesh ) then

      ! Where the field is masked, set to mask_value. If ESMF can't interpolate to a point via interpolation or extrapolation.
      ! When we made the routehandle, we returned dstStatus, which is the status of each interpolated point.  If dstStatus
      !  is certain values, that indicates the interpolation was not successful.  So, set those points to a masked value.
      !See 54.50 ESMF_REGRIDSTATUS for meaning of ESMF_REGRIDSTATUS_MAPPED
      ! dstStatus = 8 for values that were extrapolated, and anything < 4 (ESMF_REGRIDSTATUS_MAPPED) means the point was not mapped
      where (dstStatus < ESMF_REGRIDSTATUS_MAPPED ) mpas_meshes(i+1)%bdyMaskCell = mask_value ! update the mask...useful for averaging
      call print_dstStatus(localpet,dstStatus)

      ! for conservative interpolation, there's probably an ESMF bug and unmappedDstList will be incorrect
      !  (it will have too many entries--possibly more than the size of the mesh!)
      ! we can find the correct points by finding those with values < ESMF_REGRIDSTATUS_MAPPED.
      num_unmapped = size(unmappedDstList)
      call mpi_allreduce(num_unmapped, num_unmapped_tot, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
      if (localpet==0) write(*,*)'total num unmapped = ',num_unmapped_tot

         ! Use our own MPAS-based extrapolation method to do nearest neighbor extrapolation
         ! If a point is unmapped on the (i+1)th mesh, force the value at that point to the
         !  value at the nearest point on the ith mesh.
         ! Output is the 3rd argument, which is intent(inout) and is the updated bundle
      if ( trim(adjustl(extrap_method)).eq."mpas_nearest" ) then

         ! because of a likely ESMF bug about unmappedDstList when using conservative interpolation,
         !  we need to redefine unmappedDstList, and find the unique points per processor where
         ! the target points didn't map to the source grid
         if ( trim(adjustl(interp_method)) .eq. "conserve1" .or. &
              trim(adjustl(interp_method)) .eq. "conserve2" ) then 
            cell_start = mpas_meshes(i+1)%cell_start
            cell_end   = mpas_meshes(i+1)%cell_end
            nCellsPerPET = mpas_meshes(i+1)%nCellsPerPET
            allocate(my_dstStatus(1:nCellsPerPET))
            my_dstStatus(1:nCellsPerPET) = dstStatus(cell_start:cell_end)
            deallocate(unmappedDstList)
            my_num_unmapped = count(my_dstStatus < ESMF_REGRIDSTATUS_MAPPED)
            allocate(unmappedDstList(my_num_unmapped))
            call mpi_allreduce(my_num_unmapped, num_unmapped_tot, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
            if (localpet==0) write(*,*)'NEW total num unmapped = ',num_unmapped_tot
            ii = 0
            do j = 1,nCellsPerPET
               if ( my_dstStatus(j) < ESMF_REGRIDSTATUS_MAPPED ) then
                  ii = ii + 1
                  unmappedDstList(ii) = mpas_meshes(i+1)%elemIDs(j)
               endif
            enddo
         endif
         call extrapolate(localpet,mpas_meshes(i+1),large_scale_data_going_up(i+1),mpas_meshes(i),large_scale_data_going_up(i),unmappedDstList)
         call extrapolate(localpet,mpas_meshes(i+1),small_scale_data_going_up(i+1),mpas_meshes(i),small_scale_data_going_up(i),unmappedDstList)
         if ( trim(adjustl(interp_method)) .eq. "conserve1" .or.  trim(adjustl(interp_method)) .eq. "conserve2" ) deallocate(unmappedDstList,my_dstStatus)

         ! All interpolation options except esmf_creep fill the entire destination grid.
         ! If using esmf_creep, that means some points will be unmapped.
         ! Unfortunately, if a point cannot be interpolated to, ESMF sets the value to 0 on the target grid, which is annoying, b/c 0 is a physical
         !  value.  So, force those points that can't be interpolated to --> mask_value.
         ! Also, there is the issue where certain points will be influenced by a masked-out point, but the points themselves aren't 
         !  masked out. These points are identifiable, as they will be somehow related to a masked-out point and will be a large 
         ! negative number.  Force those points to mask_value, as well.
         ! Note that for esmf_creep interpolation, the bigger extrap_num_levels_creep, the more areas outside of the region will be filled
         !   and as extrap_num_levels_creep goes up, more values have dstStatus = 8.
      else if ( trim(adjustl(extrap_method)) .eq. "esmf_creep" ) then
         ! mpas_meshes(i+1)%bdyMaskCell reflects interpolating from the ith to (i+1)th mesh.
         ! This call will ensure that where the field is masked, the fields in the bundle will be set to mask_val. 
         ! Output is the final argument, which is intent(inout) and is the updated bundle
         call force_masked_data_to_value(localpet,mpas_meshes(i+1),mpas_meshes(i+1)%bdyMaskCell,mask_value,mask_value,large_scale_data_going_up(i+1))
         call force_masked_data_to_value(localpet,mpas_meshes(i+1),mpas_meshes(i+1)%bdyMaskCell,mask_value,mask_value,small_scale_data_going_up(i+1))

         ! Where the bundle is < ignore_value, set the bundle to mask_value
         ! Output is the 2nd argument, which is intent(inout) and is the updated bundle
         call force_certain_data_to_value(localpet,large_scale_data_going_up(i+1),'lt',ignore_value,mask_value)
         call force_certain_data_to_value(localpet,small_scale_data_going_up(i+1),'lt',ignore_value,mask_value)

         ! where a field in the bundle large_scale_data_going_up(i+1) < ignore_value, set mpas_meshes(i+1)%bdyMaskCell = mask_value
         ! masks are the same for all meshes, so just need to do for one field/bundle
         call update_mask(localpet,mpas_meshes(i+1),large_scale_data_going_up(i+1),ignore_value, mpas_meshes(i+1)%bdyMaskCell, mask_value)
      endif
   endif ! if ( regional_mesh )

   deallocate(dstStatus)

   call mpi_barrier(mpi_comm_world, ierr) ! sync-up before going to next mesh

enddo

!----------------------------------------------------
! If outputting on a lat/lon grid, set that up here
!----------------------------------------------------
if ( output_latlon_grid ) then
   ! define the lat-lon grid and latlon_bundle
   call define_grid_latlon(localpet,npets) ! defines latlon_grid in model_grid
   allocate(latlon_bundle(nmeshes))  ! esmf_fieldbundle type
   allocate(rh_latlon(nmeshes)) ! route-handle for the ith mesh to a common lat-lon grid
   ! The same routehandle can be used for all fields/bundles when going from the ith mesh to 
   !   the common lat-lon grid, which is why me make rh_latlon(i) up here and deallocate at the very end
   do i = 1,nmeshes
      call define_bundle(localpet,latlon_grid,latlon_bundle(i)) ! output is latlon_bundle(i)
      call make_rh(localpet,large_scale_data_going_up(i), latlon_bundle(i), 'bilinear', trim(adjustl(extrap_method_latlon)), rh_latlon(i), unmappedDstList )
   enddo
   write(cell_degrees,fmt='(f5.3)') dx_in_degrees ! lat-lon grid cell size in degrees
endif

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
      endif
   enddo

   ! Output 1 file containing all the meshes interpolated to the same lat-lon grid
   if ( output_latlon_grid ) then
      ! First output data from file providing large scales. Need to interpolate the data onto lat-lon grid first.
      do i = 1,nmeshes
         call interp_data(localpet, rh_latlon(i), large_scale_data_going_up(i), latlon_bundle(i))
      enddo
      my_output_name = 'latlon_mesh_largeScaleData_goingUp_'//trim(adjustl(cell_degrees))//'degrees.nc'
      call write_to_file_latlon(localpet,nmeshes,latlon_bundle, my_output_name)

      ! Then output data from file providing small scales. Need to interpolate the data onto lat-lon grid first.
      do i = 1,nmeshes
         call interp_data(localpet, rh_latlon(i), small_scale_data_going_up(i), latlon_bundle(i))
      enddo
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

   allocate(dstStatus(mpas_meshes(i)%nCells))
   call make_rh(localpet, large_scale_data_going_up(i+1), large_scale_data_going_down(i), trim(adjustl(interp_method)), trim(adjustl(extrap_method)), rh, unmappedDstList,dstStatus)
   call interp_data(localpet, rh, large_scale_data_going_up(i+1), large_scale_data_going_down(i))
   call interp_data(localpet, rh, small_scale_data_going_up(i+1), small_scale_data_going_down(i))
   call destroy_rh(rh)

   if ( regional_mesh ) then

      call print_dstStatus(localpet,dstStatus)
      ! force 0s from the interpolation to mask_value
     !call force_masked_data_to_value(localpet,mpas_meshes(i),mpas_meshes(i)%bdyMaskCell,mask_value,mask_value,large_scale_data_going_down(i))
     !call force_masked_data_to_value(localpet,mpas_meshes(i),mpas_meshes(i)%bdyMaskCell,mask_value,mask_value,small_scale_data_going_down(i))

      ! where the bundle is < ignore_value, set to mask_value
     !call force_certain_data_to_value(localpet,large_scale_data_going_down(i),'lt',ignore_value,mask_value)
     !call force_certain_data_to_value(localpet,small_scale_data_going_down(i),'lt',ignore_value,mask_value)
   endif

   deallocate(dstStatus)

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

   if ( regional_mesh ) then

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
   call make_rh(localpet, blending_bundle(i),blending_bundle(i-1),trim(adjustl(interp_method)), trim(adjustl(extrap_method)), rh, unmappedDstList)
   call interp_data(localpet,rh,blending_bundle(i),blending_bundle(i-1))
   call destroy_rh(rh)
   if ( regional_mesh ) then
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

!---------------------------
! Write blended data to file
!-----------------------------
! use small_scale_file as template for output file (because 4th argument is .true.)
call write_to_file(localpet, mpas_meshes(1), blending_bundle(1), .true., small_scale_file, output_blended_filename)
if ( localpet == 0 )write(*,*)'Done outputting blended file!' ! once you see this, the blended file has been successfully written

call mpi_barrier(mpi_comm_world, ierr)

! write blended data on lat-lon grid to file
if ( output_latlon_grid ) then
   call interp_data(localpet, rh_latlon(1), blending_bundle(1), latlon_bundle(1))
   my_output_name = 'latlon_mesh_BlendedFields_'//trim(adjustl(cell_degrees))//'degrees.nc'
   call write_to_file_latlon(localpet,1,latlon_bundle(1), my_output_name)
endif

!-----------------------------------------
! Output the downscaled fields if requested
!-----------------------------------------
if ( output_intermediate_files_down ) then
   do i = 1, nmeshes-1
      write(cell_dx,fmt='(f5.1)') nominal_horizontal_cell_spacing(i)
      my_output_name = 'mpas_mesh_largeScalePerts_'//trim(adjustl(cell_dx))//'km.nc'
      call write_to_file(localpet, mpas_meshes(i), large_scale_data_perts(i), .false., large_scale_file, my_output_name)
   enddo
  !do i = 2, nmeshes
  !   write(cell_dx,fmt='(f5.1)') nominal_horizontal_cell_spacing(i)
  !   my_output_name = 'mpas_mesh_blendingBundleStep_'//trim(adjustl(cell_dx))//'km.nc'
  !   call write_to_file(localpet, mpas_meshes(i), blending_bundle(i), .false., large_scale_file, my_output_name)
  !enddo
   if ( output_latlon_grid ) then
      do i = 1,nmeshes
         call interp_data(localpet, rh_latlon(i), small_scale_data_going_down(i), latlon_bundle(i))
      enddo
      my_output_name = 'latlon_mesh_smallScaleData_goingDown_'//trim(adjustl(cell_degrees))//'degrees.nc'
      call write_to_file_latlon(localpet,nmeshes-1,latlon_bundle, my_output_name)

      do i = 1,nmeshes-1
         call interp_data(localpet, rh_latlon(i), large_scale_data_perts(i), latlon_bundle(i))
      enddo
      my_output_name = 'latlon_mesh_largeScalePerts_'//trim(adjustl(cell_degrees))//'degrees.nc'
      call write_to_file_latlon(localpet,nmeshes-1,latlon_bundle, my_output_name)
   endif
endif
 
!-------------
! Clean up 
!-------------
do i = 1,nmeshes
   !write(iii,fmt='(i5)') i
   call cleanup_bundle(large_scale_data_going_up(i))
   call cleanup_bundle(small_scale_data_going_up(i))
   call cleanup_bundle(blending_bundle(i))
   call cleanup_bundle(tmp_bundle(i))
   call cleanup_mpas_mesh_type(mpas_meshes(i))
   if ( output_latlon_grid ) then
      call cleanup_bundle(latlon_bundle(i))
      call destroy_rh(rh_latlon(i))
   endif
enddo

! degbugging
!call ESMF_LogWrite("Log Write", ESMF_LOGMSG_INFO, rc=ierr)

do i = 1,nmeshes-1
   call cleanup_bundle(large_scale_data_perts(i))
   call cleanup_bundle(small_scale_data_perts(i))
   call cleanup_bundle(large_scale_data_going_down(i)) ! this is length nmeshes, but the entry for "nmeshes" is a pointer to large_scale_data_going_up ...
   call cleanup_bundle(small_scale_data_going_down(i)) !  which is deallocated above, so only need to deallocate these from 1:nmeshes-1
   if ( average_upscale_before_interp ) then
      call cleanup_bundle(large_scale_data_going_up_avg(i))
      call cleanup_bundle(small_scale_data_going_up_avg(i))
   endif
enddo

if ( average_upscale_before_interp ) deallocate(large_scale_data_going_up_avg, small_scale_data_going_up_avg)
deallocate(large_scale_data_going_up, small_scale_data_going_up)
deallocate(large_scale_data_going_down, small_scale_data_going_down)
deallocate(blending_bundle, tmp_bundle)
deallocate(mpas_meshes)
deallocate(large_scale_data_perts, small_scale_data_perts)
deallocate(nVertLevelsPerVariable)
if ( output_latlon_grid ) then 
   deallocate(latlon_bundle)
   deallocate(rh_latlon)
endif

call cleanup(localpet) ! deallocates variables allocated or defined in model_grid module

if (localpet==0) write(*,*)"- CALL ESMF_finalize"
call ESMF_finalize(endflag=ESMF_END_KEEPMPI, rc=ierr)

call mpi_finalize(ierr)

write(*,*)"- DONE."

end program mpas_blending
