module input_data

use esmf
use netcdf
use mpas_netcdf_interface, only : open_netcdf, close_netcdf, get_netcdf_dims, netcdf_err
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
public :: read_input_data, input_data_error_checks

contains

subroutine read_input_data(localpet,input_mesh,mpas_mesh,the_file,bundle)

   include 'mpif.h'

   integer, intent(in)                   :: localpet
   type(esmf_mesh), intent(in)           :: input_mesh
   type(mpas_mesh_type), intent(in)      :: mpas_mesh
   character(len=500), intent(in)        :: the_file
   type(esmf_fieldbundle), intent(inout) :: bundle

   character(len=50)               :: vname
   integer                         :: error, ncid, rc
   integer                         :: id_var, i, j, k
   integer                         :: lowBound, upperBound, nz
   integer                         :: xtype, ndims_var, natts, dimids(NF90_MAX_VAR_DIMS), dims(4)
   character(len=500)              :: varname
   character(len=NF90_MAX_NAME)    :: dim_names(4)
   integer, allocatable            :: elemIDs(:) !, nodeIDs
   integer                         :: nCells_input, nCellsPerPET

   type(esmf_field),allocatable    :: fields(:)
    
   real(esmf_kind_r8), allocatable :: dummy2(:,:), dummy3(:,:,:)
   real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)

   ! open the netCDF file with the data we want to read
   call open_netcdf(trim(the_file),ncid)

   if (localpet==0) print*, "GETTING INPUT DATA FROM", trim(adjustl(the_file))

   nCells_input = mpas_mesh%nCells
   nCellsPerPET = mpas_mesh%nCellsPerPET
   allocate(elemIDs(nCellsPerPET))
   elemIDs = mpas_mesh%elemIDs

   !-------------------------------------
   ! Initialize esmf atmospheric fields 
   !-------------------------------------

   allocate(fields(nvars_to_blend))

   do i = 1,nvars_to_blend
       
      vname = trim(adjustl(variables_to_blend(i)))

      ! get variable ID--use direct nf90 call here because we just want the id
      error=nf90_inq_varid(ncid, trim(vname), id_var)
      call netcdf_err(error, 'reading field id' )
      
      ! get dimension information about the variable
      dims(:) = 1
      dim_names(:) = ''
      error = nf90_Inquire_Variable(ncid, id_var, varname, xtype, ndims_var, dimids, natts) ! Output is varname, xtype,ndims_var, dimids
      call netcdf_err(error, 'inquiring variable' )
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
         upperBound = nz_input
      else if ( any(dim_names.eq."nSoilLevels")) then
         call get_netcdf_dims(ncid,'nSoilLevels',nsoil_input) ! define public variable
         upperBound = nsoil_input
      else if ( any(dim_names.eq."nVertLevelsP1")) then
         call get_netcdf_dims(ncid,'nVertLevelsP1',nzp1_input) ! define public variable
         upperBound = nzp1_input
      else ! 2D variable
         upperBound = 1
      endif

      lowBound = 1 
      nz = upperBound
      nVertLevelsPerVariable(i) = nz ! public from this module, allocated in mpas_blending.f90

      if (localpet==0) print*, "- INIT FIELD ", vname

      if ( nVertLevelsPerVariable(i).eq.1 ) then
         fields(i) = ESMF_FieldCreate(input_mesh, & 
                         typekind=ESMF_TYPEKIND_R8, &
                         meshloc=ESMF_MESHLOC_ELEMENT, &
                         name=variables_to_blend(i), rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
             call error_handler("IN FieldCreate", rc)

         allocate(dummy2(nCells_input,1))
         call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         error=nf90_get_var(ncid, id_var, dummy2) ! use direct nf90 routine b/c we already have id_var
         call netcdf_err(error, 'reading field' )
         do j = 1, nCellsPerPET
            varptr(j) = dummy2(elemIDs(j),1)
         enddo
         nullify(varptr)
         deallocate(dummy2)
      else
         fields(i) = ESMF_FieldCreate(input_mesh, & 
                         typekind=ESMF_TYPEKIND_R8, &
                         meshloc=ESMF_MESHLOC_ELEMENT, &
                         name=variables_to_blend(i), &
                         ungriddedLBound=(/lowBound/), &
                         ungriddedUBound=(/upperBound/), rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
             call error_handler("IN FieldCreate", rc)

         allocate(dummy3(nz,nCells_input,1))
         call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         error=nf90_get_var(ncid, id_var, dummy3) ! use direct nf90 routine b/c we already have id_var
         call netcdf_err(error, 'reading field' )
         do j = 1, nCellsPerPET
            varptr2(j,:) = dummy3(:,elemIDs(j),1)
         enddo
         nullify(varptr2)
         deallocate(dummy3)
      endif
           
      if (localpet==0) print*,"- SET ON MESH ", trim(vname)

   enddo ! end loop over number of variables

   ! create bundle now that we've defined all the fields
   bundle = ESMF_FieldBundleCreate(fieldList=fields, name="input bundle", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)

   deallocate(fields, elemIDs)

   call close_netcdf(trim(the_file),ncid)

   if (localpet==0) print*, "DONE GETTING DATA FROM ",trim(adjustl(the_file))
   if (localpet==0) print*, ""

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
