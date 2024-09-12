module bundles_mod

use esmf
use netcdf
use utils_mod, only  : error_handler
use input_data, only : nVertLevelsPerVariable
use program_setup, only  : nvars_to_blend, variables_to_blend

implicit none

private
 
! Public subroutines
public :: add_subtract_bundles, cleanup_bundle, define_bundle
public :: bundle_get_num_fields, bundle_get_fields, bundle_get_field_names

interface define_bundle
   module procedure define_bundle_mesh
   module procedure define_bundle_grid
end interface

contains

subroutine define_bundle_mesh(localpet,the_mesh,bundle)

   integer, intent(in)          :: localpet
   type(esmf_mesh), intent(in)  :: the_mesh
   type(esmf_fieldbundle), intent(inout)  :: bundle

   integer                      :: i, rc
   integer                      :: lowBound, upperBound
   type(esmf_field),allocatable :: fields(:)

   allocate(fields(nvars_to_blend))

   do i = 1, nvars_to_blend
      if ( nVertLevelsPerVariable(i) .eq. 1) then ! filled in call to read_input_data
         fields(i) = ESMF_FieldCreate(the_mesh, &
                     typekind=ESMF_TYPEKIND_R8, &
                     meshloc=ESMF_MESHLOC_ELEMENT, &
                     name=trim(adjustl(variables_to_blend(i))), rc=rc)
      else
         lowBound = 1
         upperBound = nVertLevelsPerVariable(i)
         fields(i) = ESMF_FieldCreate(the_mesh, &
                     typekind=ESMF_TYPEKIND_R8, &
                     meshloc=ESMF_MESHLOC_ELEMENT, &
                     name=trim(adjustl(variables_to_blend(i))), &
                     ungriddedLBound=(/lowBound/), &
                     ungriddedUBound=(/upperBound/), rc=rc)
      endif
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldCreate", rc)
   enddo

   ! create bundle
   bundle = ESMF_FieldBundleCreate(fieldList=fields, name="target bundle", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleCreate", rc)

   deallocate(fields)

end subroutine define_bundle_mesh

subroutine define_bundle_grid(localpet,the_mesh,bundle)

   integer, intent(in)          :: localpet
   type(esmf_grid), intent(in)  :: the_mesh
   type(esmf_fieldbundle), intent(inout)  :: bundle

   integer                      :: i, rc
   integer                      :: lowBound, upperBound
   type(esmf_field),allocatable :: fields(:)

   allocate(fields(nvars_to_blend))

   do i = 1, nvars_to_blend

      if ( nVertLevelsPerVariable(i) .eq. 1) then ! filled in call to read_input_data
         fields(i) = ESMF_FieldCreate(the_mesh, &
                     typekind=ESMF_TYPEKIND_R8, &
                     staggerloc=ESMF_STAGGERLOC_CENTER, &
                     name=trim(adjustl(variables_to_blend(i))), rc=rc)
      else
         lowBound = 1
         upperBound = nVertLevelsPerVariable(i)
         fields(i) = ESMF_FieldCreate(the_mesh, &
                     typekind=ESMF_TYPEKIND_R8, &
                     staggerloc=ESMF_STAGGERLOC_CENTER, &
                     name=trim(adjustl(variables_to_blend(i))), &
                     ungriddedLBound=(/lowBound/), &
                     ungriddedUBound=(/upperBound/), rc=rc)
      endif
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldCreate", rc)
   enddo

   ! create bundle
   bundle = ESMF_FieldBundleCreate(fieldList=fields, name="target bundle", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleCreate", rc)

   deallocate(fields)

end subroutine define_bundle_grid

! If operation == 'subtract', then return
!   bundle1 minus bundle2
! If operation == 'add', then return
!   w1*bundle1 + w2*bundle2
function add_subtract_bundles(localpet,operation,bundle1,bundle2,w1,w2) result(bundle_out)

   integer, intent(in)                   :: localpet
   character(len=*), intent(in)          :: operation
   type(esmf_fieldbundle), intent(in)    :: bundle1, bundle2 
   real, intent(in)                      :: w1, w2 ! weights on bundle1 and bundle2, respectively for addition

   type(esmf_fieldbundle) :: bundle_out ! output

   integer                       :: rc, i, nfields
   type(ESMF_Field), allocatable :: fields1(:), fields2(:), fields_out(:)
   type(ESMF_Mesh)               :: the_mesh
   real(esmf_kind_r8), pointer   :: varptr_1d1(:), varptr_1d2(:), varptr_1d_out(:)
   real(esmf_kind_r8), pointer   :: varptr_2d1(:,:), varptr_2d2(:,:), varptr_2d_out(:,:)

   !if(localpet==0)write(*,*)'In add_subtract_bundles'

    ! get the number of fields from the bundle
    call bundle_get_num_fields(bundle1, nfields)

    allocate(fields1(nfields))
    allocate(fields2(nfields))
    allocate(fields_out(nfields))

    call bundle_get_fields(bundle1, fields1)
    call bundle_get_fields(bundle2, fields2)

    ! define bundle_out
    call ESMF_FieldBundleGet(bundle1, mesh=the_mesh, rc=rc)
    call define_bundle(localpet, the_mesh, bundle_out)
    call bundle_get_fields(bundle_out, fields_out)

   ! compute perturbations
   do i = 1,nfields
      
      if ( nVertLevelsPerVariable(i).eq.1 ) then
         nullify(varptr_1d1)
         nullify(varptr_1d2)
         nullify(varptr_1d_out)

         call ESMF_FieldGet(fields1(i), farrayPtr=varptr_1d1, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN EMSF_FieldGet for bundle1", rc)

         call ESMF_FieldGet(fields2(i), farrayPtr=varptr_1d2, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN EMSF_FieldGet for bundle2", rc)

         call ESMF_FieldGet(fields_out(i), farrayPtr=varptr_1d_out, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN EMSF_FieldGet for bundle_out", rc)
    
         if (trim(adjustl(operation)) == 'subtract' ) then
            varptr_1d_out = varptr_1d1 - varptr_1d2
         else if (trim(adjustl(operation)) == 'add' ) then
            varptr_1d_out = w1*varptr_1d1 + w2*varptr_1d2
         else
            call error_handler("Invalid operation "//trim(operation)//" in add_subtract_bundles", rc)
         endif
      else

         nullify(varptr_2d1)
         nullify(varptr_2d2)
         nullify(varptr_2d_out)

         call ESMF_FieldGet(fields1(i), farrayPtr=varptr_2d1, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN EMSF_FieldGet for bundle1", rc)

         call ESMF_FieldGet(fields2(i), farrayPtr=varptr_2d2, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN EMSF_FieldGet for bundle2", rc)

         call ESMF_FieldGet(fields_out(i), farrayPtr=varptr_2d_out, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN EMSF_FieldGet for bundle_out", rc)

         if (trim(adjustl(operation)) == 'subtract' ) then
            varptr_2d_out = varptr_2d1 - varptr_2d2
         else if (trim(adjustl(operation)) == 'add' ) then
            varptr_2d_out = w1*varptr_2d1 + w2*varptr_2d2
         else
            call error_handler("Invalid operation "//trim(operation)//" in add_subtract_bundles", rc)
         endif

      endif

    enddo

    ! update the fields in target_bundle
    call ESMF_FieldBundleAddReplace(bundle_out, fields_out, rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
       call error_handler("IN ESMF_FieldBundleAddReplace", rc)

    ! clean-up
    deallocate(fields1, fields2, fields_out)

    return

end function add_subtract_bundles

subroutine cleanup_bundle(input_bundle)
   type(esmf_fieldbundle), intent(inout) :: input_bundle

   type(esmf_field), allocatable    :: fields(:)
   integer :: i, rc, nvars

   call bundle_get_num_fields(input_bundle,nvars)

   allocate(fields(nvars))
   call bundle_get_fields(input_bundle,fields)
  
   do i = 1, nvars
      call ESMF_FieldDestroy(fields(i), rc=rc)
   enddo
   call ESMF_FieldBundleDestroy(input_bundle)
   deallocate(fields)

end subroutine cleanup_bundle

! get the number of fields from the bundle
subroutine bundle_get_num_fields(my_bundle,nfields)
   type(esmf_fieldbundle), intent(in) :: my_bundle
   integer, intent(inout) :: nfields
   integer :: rc

   call ESMF_FieldBundleGet(my_bundle, fieldCount=nfields, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldCount", rc)
end subroutine bundle_get_num_fields

! get the esmf fields (esmf_field type) from the bundle
subroutine bundle_get_fields(my_bundle,fields)
   type(esmf_fieldbundle), intent(in) :: my_bundle
   type(esmf_field), dimension(:), intent(inout) :: fields
   integer :: rc

   call ESMF_FieldBundleGet(my_bundle, fieldList=fields, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldList", rc)
end subroutine bundle_get_fields

! get the names of the fields from the bundle
subroutine bundle_get_field_names(my_bundle,field_names)
   type(esmf_fieldbundle), intent(in) :: my_bundle
   character(len=*), dimension(:), intent(inout) :: field_names
   integer :: rc

   call ESMF_FieldBundleGet(my_bundle, fieldNameList=field_names, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldNameList", rc)
end subroutine bundle_get_field_names

end module bundles_mod
