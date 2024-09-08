! Code to read the namelist file
module program_setup

use esmf
use ESMF_LogPublicMod
use utils_mod, only : error_handler

implicit none

private

! Private variables only visible to this module
integer, parameter :: nvars_max  = 200

! Public variables accessible ty other modules
type(ESMF_LogKind_Flag), public :: LogType
integer, public :: nvars_to_blend
 
! Public namelist &share variables
character(len=500), public  :: large_scale_file = "NULL" ! Full path of MPAS file with large-scale information
character(len=500), public  :: small_scale_file = "NULL" ! Full path of MPAS file with small-scale information
character(len=500), public  :: output_blended_filename = "./blended.nc" ! Full path of output file to be created with blended fields
character(len=500), public  :: grid_info_file      = "NULL"   ! Full path of text file with information about meshes
character(len=500), public  :: interp_method = "bilinear" ! Interpolation method
character(len=500), public  :: variables_to_blend(nvars_max) = "NULL" ! netCDF variable names to blend
logical, public             :: average_upscale_before_interp = .false. ! If true, average to the next coarsest mesh before interpolating to that mesh
logical, public             :: output_intermediate_files_up = .false. ! If true, output a file for each intermediate mesh when progressively upscaling
logical, public             :: output_intermediate_files_down = .false. ! If true, output a file for each intermediate mesh when progressively downscaling
logical, public             :: output_blended_edge_normal_wind = .false. ! If true, output blended edge normal wind (u)
logical, public             :: esmf_log = .false. ! logging information
logical, public             :: smooth_going_downscale = .false. ! If true, run a smoother after interpolating downscale when reconstructing fields
real,    public             :: smoother_dimensionless_coefficient = 0.0 ! Not used if smooth_going_downscale = .false.; parameter for smoothing
integer, public             :: max_boundary_layer_to_keep = 7 ! The largest boundary cell (bdyMaskCell) for which we'll use in the interpolation
character(len=50), public   :: regional_masking_method = "none" ! If "none" don't do anything special. If "esmf", let ESMF figure out the points it can't interpolate to for regional meshes. If "mpas_bdy", use MPAS mesh geometry to try and do regional masking.

! Public namelist &latlon_output variables
logical, public            :: output_latlon_grid = .false. ! Switch to turn on lat-lon output (for diagnostics purposes)
logical, public            :: is_regional = .false. ! If lat-lon mesh is regional, set to .true.  If global, set to .false.
integer, public            :: nlat = -1 ! Number of latitude points on output lat-lon grid
integer, public            :: nlon = -1 ! Number of longitude points on output lat-lon grid
real, public               :: lat_ll = -1.0 ! Latitude (degrees) of lower left corner of output lat-lon grid
real, public               :: lon_ll = -1.0 ! Longitude (degrees) of lower left corner of output lat-lon grid
real, public               :: dx_in_degrees = -9999. ! Horizontal grid spacing (degrees) of output lat-lon grid

! Public subroutines
public :: read_setup_namelist

! Namelists
namelist /share/ large_scale_file, small_scale_file, output_blended_filename, &
         variables_to_blend, grid_info_file, interp_method, &
         average_upscale_before_interp, output_intermediate_files_up, &
         output_intermediate_files_down, output_blended_edge_normal_wind, esmf_log, &
         smooth_going_downscale, smoother_dimensionless_coefficient, &
         max_boundary_layer_to_keep, regional_masking_method

namelist /latlon_output/ output_latlon_grid, is_regional, &
                         nlat, nlon, lat_ll, lon_ll, dx_in_degrees

contains

subroutine read_setup_namelist(filename,localpet)

   character(len=*), intent(in) :: filename
   integer, intent(in) :: localpet

   integer :: i, ierr, unum
 
   ! Read the namelist
   unum = 41
   open(unum, file=filename, iostat=ierr)
   if (ierr /= 0) call error_handler("OPENING SETUP NAMELIST.", ierr)
   read(unum, nml=share, iostat=ierr)
   if (ierr /= 0) call error_handler("READING SETUP NAMELIST SHARE.", ierr)
   read(unum, nml=latlon_output, iostat=ierr)
   if (ierr /= 0) call error_handler("READING SETUP NAMELIST LATLON_OUTPUT.", ierr)
   close (unum)

   ! Figure out how many variables there are to blend
   nvars_to_blend = 0
   do i = 1,nvars_max
     if ( variables_to_blend(i) == "NULL" ) exit
     nvars_to_blend = nvars_to_blend + 1
   enddo
   if ( nvars_to_blend == 0 ) call error_handler("NO VARIABLES TO BLEND.", 1)

   ! Define the ESMF log type
   if (esmf_log) then
     LogType = ESMF_LOGKIND_MULTI_ON_ERROR
   else
     LogType = ESMF_LOGKIND_NONE
   endif 

   ! Check bounds on max_boundary_layer_to_keep
  !if ( max_boundary_layer_to_keep > 7 ) then
  !   if (localpet==0) write(*,*) 'Resetting max_boundary_layer_to_keep = 7'
  !   max_boundary_layer_to_keep = 7
  !endif

  !if ( max_boundary_layer_to_keep < 1 ) then
  !   if (localpet==0) write(*,*) 'Resetting max_boundary_layer_to_keep = 1'
  !   max_boundary_layer_to_keep = 1
  !endif

   if ( trim(adjustl(regional_masking_method)) == "esmf" ) then
      if (localpet==0) write(*,*)'using ESMF routines to try regional masking'
   else if ( trim(adjustl(regional_masking_method)) == "mpas_bdy" ) then
      if (localpet==0) write(*,*)'using MPAS routines to try regional masking'
   else if ( trim(adjustl(regional_masking_method)) == "none" ) then
      if (localpet==0) write(*,*)'no attempt at regional masking'
   else
      if (localpet==0) write(*,*)'regional_masking_method = '//trim(adjustl(regional_masking_method))//' is invalid. exit'
      call error_handler("INVALID REGIONAL_MASKING_METHOD", ierr)
   endif
 
end subroutine read_setup_namelist

end module program_setup
