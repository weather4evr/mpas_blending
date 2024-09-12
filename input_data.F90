module input_data

use esmf
use netcdf
use mpas_netcdf_interface, only : open_netcdf, close_netcdf, get_netcdf_dims, netcdf_err, get_netcdf_var
use utils_mod, only      : error_handler
use program_setup, only  : nvars_to_blend, variables_to_blend
use model_grid, only     : mpas_mesh_type

implicit none

private

! Public variables accessible to other modules
integer, allocatable, public :: nVertLevelsPerVariable(:)
integer, public              :: nz_input = -1    ! number of input grid atm layers
integer, public              :: nzp1_input = -1  ! number of input grid atm layer interfaces
integer, public              :: nsoil_input = -1 ! number of input soil levels

! Public subroutines
public :: read_input_data, input_data_error_checks, set_nVertLevelsPerVariable

contains

subroutine set_nVertLevelsPerVariable(localpet,the_file)

   integer, intent(in)                   :: localpet
   character(len=500), intent(in)        :: the_file

   character(len=50)            :: vname
   integer                      :: error, ncid, rc
   integer                      :: id_var, i, j, k, nz
   integer                      :: xtype, ndims_var, natts, dimids(NF90_MAX_VAR_DIMS), dims(4)
   character(len=500)           :: varname
   character(len=NF90_MAX_NAME) :: dim_names(4)

   ! open the netCDF file with the data we want to read
   call open_netcdf(trim(the_file),ncid)

   if (localpet==0) write(*,*) "GETTING INPUT DATA FROM"//trim(adjustl(the_file))

   do i = 1,nvars_to_blend
       
      vname = trim(adjustl(variables_to_blend(i)))

      ! get variable ID--use direct nf90 call here because we just want the id
      error=nf90_inq_varid(ncid, trim(vname), id_var)
      call netcdf_err(error, 'reading field id' )
      
      ! get dimension information about the variable
      dims(:) = 1
      dim_names(:) = ''
      error = nf90_Inquire_Variable(ncid, id_var, varname, xtype, ndims_var, dimids, natts) ! Output is varname, xtype,ndims_var, dimids
      call netcdf_err(error, 'inquiring variable')
      do j = 1,ndims_var
         error = nf90_inquire_dimension( ncid, dimids(j), name=dim_names(j), len=dims(j) )
         dim_names(j) = trim(adjustl(dim_names(j))) ! can't do trim on an array, so do it here
      enddo
      if ( localpet == 0 ) write(*,fmt='(a28, a30,i2,4i8)')'    variable, ndims, dims = ',trim(adjustl(varname)),ndims_var,dims(1:ndims_var)

      if ( any(dim_names == 'nEdges') ) then
         write(*,*)'Variable '//trim(adjustl(vname))//' is defined on cell edges'
         write(*,*)'This code can only deal with variables on cell centers'
         call error_handler("Not a cell-centered variable", -23)
      endif

      if ( any(dim_names.eq."nVertLevels")) then
         call get_netcdf_dims(ncid,'nVertLevels',nz_input) ! define public variable
         nz = nz_input
      else if ( any(dim_names.eq."nSoilLevels")) then
         call get_netcdf_dims(ncid,'nSoilLevels',nsoil_input) ! define public variable
         nz = nsoil_input
      else if ( any(dim_names.eq."nVertLevelsP1")) then
         call get_netcdf_dims(ncid,'nVertLevelsP1',nzp1_input) ! define public variable
         nz = nzp1_input
      else ! 2D variable
         nz = 1
      endif

      nVertLevelsPerVariable(i) = nz ! public from this module, allocated in mpas_blending.f90

   end do

   call close_netcdf(trim(the_file),ncid)

end subroutine set_nVertLevelsPerVariable

subroutine read_input_data(localpet,input_mesh,mpas_mesh,the_file,bundle)
   integer, intent(in)                   :: localpet
   type(esmf_mesh), intent(in)           :: input_mesh
   type(mpas_mesh_type), intent(in)      :: mpas_mesh
   character(len=500), intent(in)        :: the_file
   type(esmf_fieldbundle), intent(inout) :: bundle

   character(len=50)               :: vname
   integer                         :: nCells_input, nCellsPerPET
   integer                         :: lowBound, nz, i, j
   integer                         :: error, rc, ncid
   type(esmf_field),allocatable    :: fields(:)
   real(esmf_kind_r8), allocatable :: dummy2(:,:), dummy3(:,:,:)
   real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)

   nCells_input = mpas_mesh%nCells
   nCellsPerPET = mpas_mesh%nCellsPerPET

   allocate(fields(nvars_to_blend))
   call ESMF_FieldBundleGet(bundle, fieldList=fields, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet", rc)

   !-------------------------------------
   ! Get the data from the input file
   !-------------------------------------

   ! open the netCDF file with the data we want to read
   call open_netcdf(trim(the_file),ncid)

   do i = 1,nvars_to_blend

      vname = trim(adjustl(variables_to_blend(i)))

      if (localpet==0) write(*,*) "- INIT FIELD ", vname

      lowBound = 1
      nz = nVertLevelsPerVariable(i)

      if ( nz.eq.1 ) then

         allocate(dummy2(nCells_input,1))
         call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         call get_netcdf_var(ncid,trim(vname),(/1/),(/nCells_input/),dummy2(:,1))
         do j = 1, nCellsPerPET
            varptr(j) = dummy2(mpas_mesh%elemIDs(j),1)
         enddo
         nullify(varptr)
         deallocate(dummy2)
      else

         allocate(dummy3(nz,nCells_input,1))
         call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         call get_netcdf_var(ncid,trim(vname),(/1,1/),(/nz,nCells_input/),dummy3(:,:,1))
         do j = 1, nCellsPerPET
            varptr2(j,:) = dummy3(:,mpas_mesh%elemIDs(j),1)
         enddo
         nullify(varptr2)
         deallocate(dummy3)
      endif
           
      if (localpet==0) write(*,*) "- SET ON MESH ", trim(vname)

   enddo ! end loop over number of variables

   ! update the fields in target_bundle
   call ESMF_FieldBundleAddReplace(bundle, fields, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN ESMF_FieldBundleAddReplace", rc)

   deallocate(fields)

   call close_netcdf(trim(the_file),ncid)

   if (localpet==0) write(*,*) "DONE GETTING DATA FROM ",trim(adjustl(the_file))
   if (localpet==0) write(*,*)''

end subroutine read_input_data

subroutine input_data_error_checks(file1, file2)
   character(len=*), intent(in) :: file1, file2

   integer :: ncid1, ncid2, nCells1, nCells2

   ! Open file1 and get info
   call open_netcdf(trim(file1),ncid1)
   call get_netcdf_dims(ncid1,'nCells',nCells1)

   ! Open file2 and get info
   call open_netcdf(trim(file2),ncid2)
   call get_netcdf_dims(ncid2,'nCells',nCells2)

   ! Close files--no longer needed
   call close_netcdf(trim(file1),ncid1)
   call close_netcdf(trim(file2),ncid2)

   if ( nCells1 .ne. nCells2 ) then
      call error_handler("Mismatched horizontal dimensions in input data", -111)
   endif

end subroutine input_data_error_checks
 
end module input_data
