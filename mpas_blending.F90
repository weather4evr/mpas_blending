! Blend two MPAS fields together on the unstructured mesh.
! Code by Craig Schwartz, based on code by Larissa Reames.
program mpas_blending

use mpi
use esmf
use utils_mod, only      : error_handler

use program_setup, only  : LogType, read_setup_namelist,  &
                           large_scale_file, small_scale_file, grid_info_file, &
                           average_upscale_before_interp, &
                           output_intermediate_files_up, output_intermediate_files_down, &
                           interp_method, output_latlon_grid, output_blended_filename, &
                           dx_in_degrees, nvars_to_blend
use model_grid, only     : define_grid_mpas, define_grid_latlon, mpas_mesh_type, mpas_meshes, &
                           nmeshes, meshes, grid_files_heirarchy, read_grid_info_file, weights, &
                           latlon_grid, cleanup, nominal_horizontal_cell_spacing
use input_data, only     : read_input_data, &
                           input_data_error_checks, nVertLevelsPerVariable
use interp, only         : interp_data, radially_average
use bundles_mod, only    : add_subtract_bundles, cleanup_bundles, define_bundle, &
                           large_scale_data_going_up, small_scale_data_going_up, &
                           large_scale_data_going_down, small_scale_data_going_down, &
                           large_scale_data_perts, small_scale_data_perts, &
                           latlon_bundle, blending_bundle
use write_data, only     : write_to_file, write_to_file_latlon

implicit none

integer                   :: ierr, localpet, npets, i
real                      :: w1, w2
character(len=500)        :: my_output_name
character(len=5)          :: cell_dx, cell_degrees
type(esmf_vm)             :: vm
type(esmf_fieldbundle), allocatable :: tmp_bundle(:)

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
print*,'- LOCAL PET ',localpet

!--------------------------------
! All processors read the namelist
!--------------------------------
call read_setup_namelist('input.nml')
allocate(nVertLevelsPerVariable(nvars_to_blend)) ! try to get rid of this variable somehow....

!--------------------------------------------
! Get the mesh information for the hierarchy
!  meshes(:) has the complete hierarchy, where
!  the ith element of meshes() is the mesh for the 
!  ith entry in grid_info_file
!--------------------------------------------
call read_grid_info_file(localpet,grid_info_file) ! sets grid_files_heirarchy, nmeshes
allocate(meshes(nmeshes))      ! esmf_mesh type
allocate(mpas_meshes(nmeshes)) ! mpas_mesh_type type defined in model_grid
do i = 1,nmeshes
   call define_grid_mpas(localpet, npets, grid_files_heirarchy(i), meshes(i),mpas_meshes(i)) ! output is meshes(i), mpas_meshes(i)
enddo

!---------------------------------------------
! Error checking to make sure large_scale_file 
!   and small_scale_file have the same meshes--they should. 
! Also ensure that input data has the same dimensions as the first file in grid_info_file
!  If there's an error, program will abort.
!---------------------------------------------
call input_data_error_checks(large_scale_file,small_scale_file)
call input_data_error_checks(large_scale_file,grid_files_heirarchy(1))

!----------------------------------------------
! Read data from input files
!  output is e.g., large_scale_data_going_up(1)
!----------------------------------------------
allocate(large_scale_data_going_up(nmeshes)) ! esmf_fieldbundle types
allocate(small_scale_data_going_up(nmeshes))
call read_input_data(localpet,meshes(1),mpas_meshes(1),large_scale_file,large_scale_data_going_up(1)) ! data from file providing large scales
call read_input_data(localpet,meshes(1),mpas_meshes(1),small_scale_file,small_scale_data_going_up(1)) ! data from file providing small scales

! read_input_data defines the bundle for i = 1, which holds the data for the highest resolution mesh
! now we need to define the bundle for i = 2,nmeshes, but we don't associate any data with them yet;
! the ith bundle will be associated with the ith mesh
!  ... the data will get filled in calls to interp_data
! More broadly, define_bundle can be used to define any bundle and 
!  associate that bundle with any mesh
do i = 2,nmeshes
   call define_bundle(localpet,meshes(i),large_scale_data_going_up(i))
   call define_bundle(localpet,meshes(i),small_scale_data_going_up(i))
enddo

!--------------------------------------------------
! Upscale the fields through the heirarchy 
!  for both large_scale_file and small_scale_file
! Upscale from the ith grid to the (i+1)th grid.
!  both input and output are esmf_fieldbundle type
! large_scale_data_going_up(1) and small_scale_data_going_up(1)
!  were filled in call to read_input_data
do i = 1,nmeshes-1
  ! If requested, average data on the current mesh to the scale of the next mesh, which is given by
  !   nominal_horizontal_cell_spacing(i+1)
  ! output is large_scale_data_going_up(i), small_scale_data_going_up(i), which is the spatially-averaged
  ! field on the ith mesh
   if ( average_upscale_before_interp ) then
      call radially_average(localpet,mpas_meshes(i),nominal_horizontal_cell_spacing(i+1),large_scale_data_going_up(i))
      call radially_average(localpet,mpas_meshes(i),nominal_horizontal_cell_spacing(i+1),small_scale_data_going_up(i))
   endif
   call interp_data(localpet, large_scale_data_going_up(i), large_scale_data_going_up(i+1), interp_method)
   call interp_data(localpet, small_scale_data_going_up(i), small_scale_data_going_up(i+1), interp_method)
enddo

!-----------------------------------------
! Output the upscaled fields if requested
!-----------------------------------------
if ( output_intermediate_files_up ) then
   ! output one grid per file, because all meshes are different.
   ! In write_to_file, .false. means define the output mesh from input data and the next entry is ignored
   !  If it's .true., then the next entry provides a template for the output file
   do i = 1, nmeshes
      write(cell_dx,fmt='(f5.1)') nominal_horizontal_cell_spacing(i)
      my_output_name = 'mpas_mesh_largeScaleData_goingUp_'//trim(adjustl(cell_dx))//'km.nc'
      call write_to_file(localpet, mpas_meshes(i), large_scale_data_going_up(i), .false., large_scale_file, my_output_name)
      my_output_name = 'mpas_mesh_smallScaleData_goingUp_'//trim(adjustl(cell_dx))//'km.nc'
      call write_to_file(localpet, mpas_meshes(i), small_scale_data_going_up(i), .false., large_scale_file, my_output_name)
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
      do i = 1,nmeshes
         call interp_data(localpet, large_scale_data_going_up(i), latlon_bundle(i), 'bilinear')
      enddo
      my_output_name = 'latlon_mesh_largeScaleData_goingUp_'//trim(adjustl(cell_degrees))//'degrees.nc'
      call write_to_file_latlon(localpet,nmeshes,latlon_bundle, my_output_name)

      ! Then output data from file providing small scales. Need to interpolate the data onto lat-lon grid first.
      do i = 1,nmeshes
         call interp_data(localpet, small_scale_data_going_up(i), latlon_bundle(i), 'bilinear')
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
   call interp_data(localpet, large_scale_data_going_up(i+1), large_scale_data_going_down(i), interp_method)
   call interp_data(localpet, small_scale_data_going_up(i+1), small_scale_data_going_down(i), interp_method)
enddo
! fill in the last (coarsest) element from the upscaled data
small_scale_data_going_down(nmeshes) = small_scale_data_going_up(nmeshes) ! LHS is basically a pointer to RHS
large_scale_data_going_down(nmeshes) = large_scale_data_going_up(nmeshes)

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
w1 = real(weights(i)) ! weights to data providing large-scale data
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
do i = nmeshes, 2, -1
   if ( localpet == 0 ) write(*,*)'Blending for mesh ',i
   w1 = real(weights(i)) ! weights to data providing large-scale data
   w2 = real(1.0 - w1)
   call interp_data(localpet,blending_bundle(i),blending_bundle(i-1),interp_method) !'bilinear')

   tmp_bundle(i-1) = add_subtract_bundles(localpet, 'add', &
                       large_scale_data_perts(i-1),small_scale_data_perts(i-1), w1, w2)

   blending_bundle(i-1) = add_subtract_bundles(localpet, 'add', &
                    blending_bundle(i-1),tmp_bundle(i-1), 1.0, 1.0)
enddo

!-----------------------------------------
! Output the downscaled fields if requested
!-----------------------------------------
if ( output_intermediate_files_down ) then
   if ( output_latlon_grid ) then
      write(cell_degrees,fmt='(f5.3)') dx_in_degrees
      do i = 1,nmeshes
         call interp_data(localpet, small_scale_data_going_down(i), latlon_bundle(i), 'bilinear')
      enddo
      my_output_name = 'latlon_mesh_smallScaleData_goingDown_'//trim(adjustl(cell_degrees))//'degrees.nc'
      call write_to_file_latlon(localpet,nmeshes-1,latlon_bundle, my_output_name)
   endif
endif
 
!--------------------------
! Write blended data to file
!--------------------------
! use small_scale_file as template for output file (because 4th argument is .true.)
call write_to_file(localpet, mpas_meshes(1), blending_bundle(1), .true., small_scale_file, output_blended_filename)

if ( output_latlon_grid ) then
    call interp_data(localpet, blending_bundle(1), latlon_bundle(1), 'bilinear')
    write(cell_degrees,fmt='(f5.3)') dx_in_degrees
    my_output_name = 'latlon_mesh_BlendedFields_'//trim(adjustl(cell_degrees))//'degrees.nc'
    call write_to_file_latlon(localpet,1,latlon_bundle(1), my_output_name)
endif
 
!-------------
! Clean up 
!-------------
do i = 1,nmeshes
   call cleanup_bundles(large_scale_data_going_up(i))
   call cleanup_bundles(small_scale_data_going_up(i))
   call cleanup_bundles(large_scale_data_going_down(i))
   call cleanup_bundles(small_scale_data_going_down(i))
   call cleanup_bundles(blending_bundle(i))
   call cleanup_bundles(tmp_bundle(i))
   if ( output_latlon_grid ) call cleanup_bundles(latlon_bundle(i))
   deallocate(mpas_meshes(i)%ElemIDs)
   deallocate(mpas_meshes(i)%latCell)
   deallocate(mpas_meshes(i)%lonCell)
   deallocate(mpas_meshes(i)%nEdgesOnCell)
   deallocate(mpas_meshes(i)%cellsOnCell)
   deallocate(mpas_meshes(i)%areaCell)
   deallocate(mpas_meshes(i)%bdyMaskCell)
enddo

do i = 1,nmeshes-1
   call cleanup_bundles(large_scale_data_perts(i))
   call cleanup_bundles(small_scale_data_perts(i))
enddo

deallocate(large_scale_data_going_up, small_scale_data_going_up)
deallocate(large_scale_data_going_down, small_scale_data_going_down)
deallocate(large_scale_data_perts, small_scale_data_perts)
deallocate(blending_bundle, tmp_bundle)
deallocate(mpas_meshes)
deallocate(nVertLevelsPerVariable)
if ( output_latlon_grid ) deallocate(latlon_bundle)

call cleanup(localpet) ! deallocates variables allocated in model_grid module

if (localpet==0) print*,"- CALL ESMF_finalize"
call ESMF_finalize(endflag=ESMF_END_KEEPMPI, rc=ierr)

call mpi_finalize(ierr)

print*,"- DONE."

end program mpas_blending
