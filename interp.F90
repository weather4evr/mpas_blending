module interp

use mpi
use esmf
use netcdf
use utils_mod, only : error_handler, sphere_distance, nearest_cell, gather_on_all_procs
use model_grid, only : mpas_mesh_type, mask_value, ignore_value, INSIDE, UNMARKED, regional_mesh
use input_data, only : nVertLevelsPerVariable
use bundles_mod, only : bundle_get_num_fields, bundle_get_fields, bundle_get_field_names
use program_setup, only  : variables_to_blend, extrap_num_levels_creep, interp_method, extrap_method

implicit none

private
 
! Public subroutines
public :: interp_data, radially_average, apply_smoothing_filter
public :: make_rh, destroy_rh
public :: force_masked_data_to_value, force_certain_data_to_value, update_mask
public :: set_global_extrapolation, extrapolate

! Variables private to this module
integer ::  isrctermprocessing = 1
integer, parameter :: num_neigh_max = 1000 ! max number of points in a neighborhood about a point
real,    parameter :: PI = 3.14159265359
real :: start, finish
real :: mpas_mesh_radius ! put this up here so all subroutines can see it
 
contains
 
subroutine interp_data(localpet,my_rh,input_bundle,target_bundle)
 
   integer, intent(in)                   :: localpet
   type(ESMF_RouteHandle), intent(inout) :: my_rh
   type(esmf_fieldbundle), intent(in)    :: input_bundle
   type(esmf_fieldbundle), intent(inout) :: target_bundle

   integer                          :: rc, i, nfields, ungriddedDimCount, k
   integer                          :: ungriddedUBound(10), ungriddedLBound(10)
   character(len=500), allocatable    :: field_names(:)
   type(ESMF_RegridMethod_Flag)     :: method
   type(ESMF_Field), allocatable    :: fields_input_grid(:), fields_target_grid(:)

   real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)
   type(ESMF_Mesh) :: mesh
   integer, allocatable :: elementMask(:)
   type(ESMF_GeomType_Flag) :: geomtype

   !if ( ESMF_RouteHandleIsCreated(routehandle, rc=rc)) then   ... 
   ! this will make  copy of a route handle (Not used)
   !my_rh = ESMF_RouteHandleCreate(rh_bilin, rc=rc)

   ! get the number of fields from the bundle
   call bundle_get_num_fields(input_bundle,nfields)

   ! get the names of the fields from the bundle (just for printing purposes)
   allocate(field_names(nfields))
   call bundle_get_field_names(input_bundle,field_names)

   ! get the input and target fields
   allocate(fields_input_grid(nfields))
   allocate(fields_target_grid(nfields))
   call bundle_get_fields(input_bundle,fields_input_grid)
   call bundle_get_fields(target_bundle,fields_target_grid)
   
   do i = 1,nfields
   
      ! do the interpolation
      call ESMF_FieldRegrid(fields_input_grid(i), fields_target_grid(i), my_rh, rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN FieldRegrid", rc)

      ! testing about getting number of ungridded dimensions...can probably remove someday
      !call ESMF_FieldGet(fields_input_grid(i),ungriddedDimCount=ungriddedDimCount, &
      !          ungriddedLBound=ungriddedLBound, ungriddedUBound=ungriddedUBound ,rc=rc)
   enddo

  !if(localpet==0)write(*,*)'Done interpolating fields'

   ! update the fields in target_bundle
   call ESMF_FieldBundleAddReplace(target_bundle, fields_target_grid, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN ESMF_FieldBundleAddReplace", rc)

   ! clean-up
   deallocate(fields_input_grid,fields_target_grid)
   deallocate(field_names)

end subroutine interp_data

! at locations where mask_field == this_mask_value, set the fields in 'bundle' to target_value
subroutine force_masked_data_to_value(localpet,mpas_mesh,mask_field,this_mask_value,target_value,bundle)

   integer, intent(in)                   :: localpet
   type(mpas_mesh_type), intent(in)      :: mpas_mesh
   integer, dimension(mpas_mesh%nCells), intent(in) :: mask_field ! on full field
   integer, intent(in)                   :: this_mask_value, target_value
   type(esmf_fieldbundle), intent(inout) :: bundle

   integer                          :: rc, i, nfields, j, k, elementCount
   integer                          :: cell_start, cell_end
   integer, allocatable             :: my_mask(:)
   type(ESMF_Field), allocatable    :: fields(:)
   real(esmf_kind_r8), pointer      :: varptr(:), varptr2(:,:)
   type(ESMF_Mesh)          :: mesh
   type(ESMF_GeomType_Flag) :: geomtype

   call ESMF_FieldBundleGet(bundle, geomtype=geomtype)
   if ( geomtype == ESMF_GEOMTYPE_MESH ) then ! only do for meshes, not ESMF grids
      if(localpet==0)write(*,*)'In force_masked_data_to_value'

      ! get the number of fields from the bundle
      call bundle_get_num_fields(bundle,nfields)

      ! get the fields from the bunlde
      allocate(fields(nfields))
      call bundle_get_fields(bundle,fields)

      ! sanity check
      call ESMF_FieldBundleGet(bundle, mesh=mesh)
      call ESMF_MeshGet(mesh, elementCount=elementCount, rc=rc)
      if ( elementCount /= mpas_mesh%nCellsPerPET ) then
         call error_handler("Wrong number of cells in force_masked_data_to_value", -43)
      endif

      ! Input data "mask_field" is over the whole mesh, so we need to know the per-PE indices
      cell_start = mpas_mesh%cell_start
      cell_end = mpas_mesh%cell_end
      allocate(my_mask(1:mpas_mesh%nCellsPerPET))
      my_mask(1:mpas_mesh%nCellsPerPET) = mask_field(cell_start:cell_end)

      do i = 1,nfields
         if ( nVertLevelsPerVariable(i) == 1 ) then
            call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc) ! varptr is decomposed across all processors
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
               call error_handler("IN FieldGet", rc)
            ! usually, this will set 0 --> mask_value
            where( my_mask == this_mask_value)  varptr = real(target_value*1.0,esmf_kind_r8) ! set to target value where my_mask == this_mask_value
            nullify(varptr)
         else
            call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc) ! varptr2 is decomposed across all processors
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
               call error_handler("IN FieldGet", rc)
            ! usually, this will set 0 --> mask_value
            do k = 1,nVertLevelsPerVariable(i)
               where( my_mask == this_mask_value) varptr2(:,k) = real(target_value*1.0,esmf_kind_r8)
            enddo
            nullify(varptr2)
         endif
      enddo

      ! update the fields in the bundle
      call ESMF_FieldBundleAddReplace(bundle, fields, rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN ESMF_FieldBundleAddReplace", rc)

      ! clean-up
      deallocate(fields, my_mask)

   endif ! endif geomtype == ESMF_GEOMTYPE_MESH

end subroutine force_masked_data_to_value

! at locations where fields in 'bundle' are < or > threshold_value (depending on 'operation'),
!set to target_value
subroutine force_certain_data_to_value(localpet,bundle,operation,threshold_value,target_value)

   integer, intent(in)                   :: localpet
   type(esmf_fieldbundle), intent(inout) :: bundle
   character(len=*), intent(in)          :: operation
   integer, intent(in)                   :: threshold_value, target_value

   integer                          :: rc, i, nfields, j, k, elementCount
   type(ESMF_Field), allocatable    :: fields(:)
   real(esmf_kind_r8), pointer      :: varptr(:), varptr2(:,:)
   type(ESMF_GeomType_Flag) :: geomtype

   call ESMF_FieldBundleGet(bundle, geomtype=geomtype)
   if ( geomtype == ESMF_GEOMTYPE_MESH ) then ! only do for meshes, not ESMF grids
      if(localpet==0)write(*,*)'In force_certain_data_to_value'

      ! get the number of fields from the bundle
      call bundle_get_num_fields(bundle,nfields)

      ! get the fields from the bunlde
      allocate(fields(nfields))
      call bundle_get_fields(bundle,fields)

      do i = 1,nfields
         if ( nVertLevelsPerVariable(i) == 1 ) then
            call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc) ! varptr is decomposed across all processors
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
               call error_handler("IN FieldGet", rc)
           ! usually used to find points influenced by masked out values
            if ( trim(adjustl(operation)) == 'lt' ) then
               where( varptr < real(threshold_value*1.0,esmf_kind_r8)) varptr = real(target_value*1.0,esmf_kind_r8) ! points influenced by masked out values
            else if ( trim(adjustl(operation)) == 'gt' ) then
               where( varptr > real(threshold_value*1.0,esmf_kind_r8)) varptr = real(target_value*1.0,esmf_kind_r8) ! points influenced by masked out values
            endif
            nullify(varptr)
         else
            call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc) ! varptr2 is decomposed across all processors
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
               call error_handler("IN FieldGet", rc)
            do k = 1,nVertLevelsPerVariable(i)
               if ( trim(adjustl(operation)) == 'lt' ) then
                  where( varptr2(:,k) < real(threshold_value*1.0,esmf_kind_r8) ) varptr2(:,k) = real(target_value*1.0,esmf_kind_r8)
               else if ( trim(adjustl(operation)) == 'gt' ) then
                  where( varptr2(:,k) > real(threshold_value*1.0,esmf_kind_r8) ) varptr2(:,k) = real(target_value*1.0,esmf_kind_r8)
               endif
            enddo
            nullify(varptr2)
         endif
      enddo

      ! update the fields in the bundle
      call ESMF_FieldBundleAddReplace(bundle, fields, rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN ESMF_FieldBundleAddReplace", rc)

      ! clean-up
      deallocate(fields)

   endif ! endif geomtype == ESMF_GEOMTYPE_MESH

end subroutine force_certain_data_to_value

! if the data in the bundle is < threshold_value, set mask_field to target_value
subroutine update_mask(localpet,mpas_mesh,bundle,threshold_value,mask_field,target_value)

   integer, intent(in)                   :: localpet
   type(mpas_mesh_type), intent(in)      :: mpas_mesh
   type(esmf_fieldbundle), intent(in)    :: bundle
   integer, intent(in)                   :: threshold_value
   integer, dimension(mpas_mesh%nCells), intent(inout) :: mask_field ! on full field
   integer, intent(in)                   :: target_value

   integer                          :: rc, i, nfields, j, k, elementCount
   integer                          :: cell_start, cell_end
   integer, allocatable             :: my_mask(:)
   real(esmf_kind_r8), pointer      :: varptr(:), varptr2(:,:)
   integer(ESMF_KIND_I4), pointer   :: varptr_mask(:)
   type(ESMF_Field)                 :: input_field, bdyMaskField
   type(ESMF_Mesh)          :: mesh
   type(ESMF_GeomType_Flag) :: geomtype

   call ESMF_FieldBundleGet(bundle, geomtype=geomtype)
   if ( geomtype == ESMF_GEOMTYPE_MESH ) then ! only do for meshes, not ESMF grids

      ! mpas_mesh%bdyMaskCell is over the whole mesh, so we need to know the per-PE indices
      cell_start = mpas_mesh%cell_start
      cell_end = mpas_mesh%cell_end
      allocate(my_mask(1:mpas_mesh%nCellsPerPET))
      my_mask(1:mpas_mesh%nCellsPerPET) = mask_field(cell_start:cell_end)

      ! create an ESMF field for bdyMaskCell; this has the ESMF parallelization, so the field
      !  is defined on each processor
      ! define the field, then initialize it to the mask; it may be modified later
      call ESMF_FieldBundleGet(bundle, mesh=mesh)
      bdyMaskField = ESMF_FieldCreate(mesh, &
                     typekind=ESMF_TYPEKIND_I4, &
                     meshloc=ESMF_MESHLOC_ELEMENT, &
                     name='bdyMaskField', rc=rc)
      call ESMF_FieldGet(bdyMaskField, farrayPtr=varptr_mask, rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
         call error_handler("IN FieldGet bdyMaskField", rc)
      do j = 1,mpas_mesh%nCellsPerPET
         varptr_mask(j) = my_mask(j) ! initialize
      enddo

      ! Get a field from the bundle. it doesn't matter which one, since all fields have identical masks
      call ESMF_FieldBundleGet(bundle, trim(adjustl(variables_to_blend(1))), field=input_field, rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN ESMF_FieldBundleGet", rc)

      if ( nVertLevelsPerVariable(1) == 1 ) then
         call ESMF_FieldGet(input_field, farrayPtr=varptr, rc=rc) ! varptr is decomposed across all processors
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         do j = 1,mpas_mesh%nCellsPerPET
            if ( varptr(j) < real(threshold_value*1.0,esmf_kind_r8)) varptr_mask(j) = target_value
         enddo
         nullify(varptr)
      else
         call ESMF_FieldGet(input_field, farrayPtr=varptr2, rc=rc) ! varptr2 is decomposed across all processors
         do j = 1,mpas_mesh%nCellsPerPET
            if ( varptr2(j,1) < real(threshold_value*1.0,esmf_kind_r8) ) varptr_mask(j) = target_value
         enddo
         nullify(varptr2)
      endif

      ! we want mask_field on the full grid on all processors
      call ESMF_FieldGather(bdyMaskField, mask_field, rootPet=0, rc=rc) ! gather on process 0
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGather", rc)
      call mpi_bcast(mask_field,mpas_mesh%nCells,mpi_integer,0,mpi_comm_world,rc) ! broadcast to all processors

      ! clean-up
      nullify(varptr_mask)
      deallocate(my_mask)

      ! done with bdyMaskField, so get rid of it.
      call ESMF_FieldDestroy(bdyMaskField, rc=rc)

   endif

end subroutine update_mask

subroutine extrapolate(localpet,mpas_mesh,bundle,mpas_mesh2,bundle2,unmappedPoints)
 
   integer, intent(in)                   :: localpet
   type(mpas_mesh_type), intent(in)      :: mpas_mesh
   type(esmf_fieldbundle), intent(inout) :: bundle
   type(mpas_mesh_type), intent(in)      :: mpas_mesh2
   type(esmf_fieldbundle), intent(in) :: bundle2
   integer(ESMF_KIND_I4), pointer, intent(in)  :: unmappedPoints(:)

   integer                          :: rc, i, nfields, j, k,nz
   integer                          :: my_point, closest_cell
   integer                          :: num_unmapped
   type(ESMF_Field), allocatable    :: fields(:), fields2(:)
   real(esmf_kind_r8), pointer      :: varptr(:), varptr2(:,:)
   real*8, allocatable :: dummy1(:), dummy2(:,:)
   type(ESMF_GeomType_Flag) :: geomtype

   call ESMF_FieldBundleGet(bundle, geomtype=geomtype)
   if ( geomtype == ESMF_GEOMTYPE_MESH ) then ! only do for meshes, not ESMF grids
      if(localpet==0)write(*,*)'In extrapolate'

      ! get the number of fields from the bundle
      call bundle_get_num_fields(bundle,nfields)

      ! get the fields from the bundle
      allocate(fields(nfields))
      call bundle_get_fields(bundle,fields)

      ! get the fields from the bundle on the other mesh
      allocate(fields2(nfields))
      call bundle_get_fields(bundle2,fields2)

      num_unmapped = size(unmappedPoints) ! unique to this processor

      do i = 1,nfields
         if ( nVertLevelsPerVariable(i) == 1 ) then
            call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc) ! varptr is decomposed across all processors
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
               call error_handler("IN FieldGet", rc)

            ! get the full-field data from the second input mesh--probably a higher-resolution mesh
            allocate(dummy1(mpas_mesh2%nCells)) ! full field on all processors
            call gather_on_all_procs(localpet, fields2(i), (/mpas_mesh2%nCells/), dummy1)

            ! loop over all the unmapped points, which should be the same for all fields, but unique for each processor
            do j = 1,num_unmapped
               my_point = unmappedPoints(j) ! on this processor only, but global cell IDs
               ! get the closest cell on mpas_mesh2 (across the whole mesh) to my_point
               closest_cell  = nearest_cell( mpas_mesh%latCell(my_point), mpas_mesh%lonCell(my_point), &
                       1, mpas_mesh2%nCells, mpas_mesh2%maxEdges, &
                       mpas_mesh2%nEdgesOnCell, mpas_mesh2%cellsOnCell, &
                       mpas_mesh2%latCell, mpas_mesh2%lonCell)

               ! find the index on this processor corresponding to my_point, and set the value
               !  there to the value at the nearest cell on mpas_mesh2
               do k = 1,mpas_mesh%nCellsPerPET
                  if ( mpas_mesh%elemIDs(k) == my_point ) then
                    !write(*,*)'closest_cell/lat/lon = ',closest_cell,mpas_mesh%latCell(my_point)*180/3.14, mpas_mesh%lonCell(my_point)*180/3.14
                     varptr(k) = dummy1(closest_cell)
                     exit
                  endif
               enddo
            enddo
           !if(localpet == 0 ) write(*,*)'done extrapolating'
            deallocate(dummy1)
            nullify(varptr)
         else
            call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc) ! varptr2 is decomposed across all processors
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
               call error_handler("IN FieldGet", rc)

            nz = nVertLevelsPerVariable(i)

            ! get the full-field data from the second input mesh--probably a higher-resolution mesh
            allocate(dummy2(mpas_mesh2%nCells,nz)) ! on all processors
            call gather_on_all_procs(localpet, fields2(i), (/mpas_mesh2%nCells,nz/), dummy2)

            ! loop over all the unmapped points, which should be the same for all fields
            do j = 1,num_unmapped
               my_point = unmappedPoints(j) ! on this processor only, but global cell IDs, which is nice
               ! get the closest cell on mpas_mesh2 (across the whole mesh)
               closest_cell  = nearest_cell( mpas_mesh%latCell(my_point), mpas_mesh%lonCell(my_point), &
                       1, mpas_mesh2%nCells, mpas_mesh2%maxEdges, &
                       mpas_mesh2%nEdgesOnCell, mpas_mesh2%cellsOnCell , &
                       mpas_mesh2%latCell, mpas_mesh2%lonCell)

               ! find the index on this processor corresponding to my_point, and set the value
               !  there to the value at the nearest cell on mpas_mesh2
               do k = 1,mpas_mesh%nCellsPerPET
                  if ( mpas_mesh%elemIDs(k) == my_point ) then
                     varptr2(k,:) = dummy2(closest_cell,:)
                     exit
                  endif
               enddo
            enddo
           !if(localpet == 0 ) write(*,*)'done extrapolating'
            deallocate(dummy2)
            nullify(varptr2)
         endif
      enddo ! loop over fields

      ! update the fields in bundle
      call ESMF_FieldBundleAddReplace(bundle, fields, rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN ESMF_FieldBundleAddReplace", rc)

      ! clean-up
      deallocate(fields,fields2)

   endif ! endif geomtype == ESMF_GEOMTYPE_MESH

end subroutine extrapolate

subroutine make_rh(localpet,input_bundle,target_bundle,my_interp_method,my_extrap_method,rh,unmappedDstList,dstStatus)

   integer, intent(in)                   :: localpet
   type(esmf_fieldbundle), intent(in)    :: input_bundle
   type(esmf_fieldbundle), intent(in)    :: target_bundle
   character(len=*), intent(in)          :: my_interp_method
   character(len=*), intent(in)          :: my_extrap_method
   type(ESMF_RouteHandle), intent(inout) :: rh
   integer(ESMF_KIND_I4), pointer, intent(inout)  :: unmappedDstList(:)
   integer, dimension(:), intent(inout), optional :: dstStatus

   integer :: rc, nCells_total
   type(ESMF_RegridMethod_Flag) :: method
   type(ESMF_Field)             :: input_field, target_field, dstStatusField
   type(ESMF_Mesh)              :: mesh
   type(ESMF_Grid)              :: grid
   type(ESMF_GeomType_Flag)     :: geomtype
   type(ESMF_NormType_Flag)     :: normType
   type(ESMF_ExtrapMethod_Flag) :: extrapMethod

   ! get a field from the input and target bundles
   call ESMF_FieldBundleGet(input_bundle, trim(adjustl(variables_to_blend(1))), field=input_field, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldCount", rc)

   call ESMF_FieldBundleGet(target_bundle, trim(adjustl(variables_to_blend(1))), field=target_field, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldCount", rc)

   ! make the interpolation route handle
   if ( my_interp_method.eq."patch") then
      method = ESMF_REGRIDMETHOD_PATCH
   else if ( my_interp_method.eq."nearest") then
      method = ESMF_REGRIDMETHOD_NEAREST_STOD
   else if (my_interp_method.eq."conserve1") then
      method = ESMF_REGRIDMETHOD_CONSERVE
   else if (my_interp_method.eq."conserve2") then
      method = ESMF_REGRIDMETHOD_CONSERVE_2ND
   else if (my_interp_method.eq."bilinear") then
      method = ESMF_REGRIDMETHOD_BILINEAR
   else
      call error_handler("Invalid interpolation method = "//my_interp_method, -223)
   endif
   if(localpet==0)write(*,*)'In make_rh; using interp method '//trim(my_interp_method)//' and extrap method '//trim(my_extrap_method)

   ! figure out the extrapolation option (set extrapMethod)
   call set_local_extrapolation_option(my_interp_method,my_extrap_method,extrapMethod)

   ! make an ESMF field that will hold the status for each destination point
   !  in the interpolation. We can use it to find points that could not be
   !  interpolated to.
   call ESMF_FieldBundleGet(target_bundle, geomtype=geomtype)
   if ( geomtype == ESMF_GEOMTYPE_MESH ) then
      call ESMF_FieldBundleGet(target_bundle, mesh=mesh)
      dstStatusField = ESMF_FieldCreate(mesh, &
                       typekind=ESMF_TYPEKIND_I4, &
                       meshloc=ESMF_MESHLOC_ELEMENT, &
                       name='dstStatusField', rc=rc)
   else
      call ESMF_FieldBundleGet(target_bundle, grid=grid)
      dstStatusField = ESMF_FieldCreate(grid, &
                       typekind=ESMF_TYPEKIND_I4, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, &
                       name='dstStatusField', rc=rc)
   endif

   ! now make the routehandle
   call ESMF_FieldRegridStore(input_field, target_field, &
                               regridmethod=method, &
                               routehandle=rh, &
                               srcTermProcessing=isrctermprocessing, &
                               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                               dstStatusField=dstStatusField, &
                               unmappedDstList=unmappedDstList, &
                               extrapMethod=extrapMethod, &
                               extrapNumLevels=extrap_num_levels_creep, &
                               normType=ESMF_NORMTYPE_FRACAREA, &
                               rc=rc)
                             ! srcMaskValues=(/mask_value/), &
                             ! dstMaskValues=(/mask_value/), &
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldRegridStore", rc)

   ! values in unmappedDstList are the GLOBAL cell IDs across the whole MPAS mesh, which is nice
  !write(*,*)'mype, size(list), min/min = ',localpet,size(unmappedDstList),minval(unmappedDstList),maxval(unmappedDstList)

   ! we want dstStatusField on the full grid on all processors
   if ( present(dstStatus)) then
      nCells_total = size(dstStatus)
      call ESMF_FieldGather(dstStatusField, dstStatus, rootPet=0, rc=rc) ! gather on process 0
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGather", rc)
      call mpi_bcast(dstStatus,nCells_total,mpi_integer,0,mpi_comm_world,rc) ! broadcast to all processors
   endif

   ! done with dstStatusField, so get rid of it.
   call ESMF_FieldDestroy(dstStatusField, rc=rc)

end subroutine make_rh

! get rid of the route handle
subroutine destroy_rh(rh)
   type(ESMF_RouteHandle), intent(inout) :: rh
   integer :: rc

   call ESMF_FieldRegridRelease(routehandle=rh, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridRelease", rc)
end subroutine destroy_rh

subroutine radially_average(localpet,mpas_mesh,averaging_radius,input_bundle,output_bundle)
 
   integer, intent(in)                   :: localpet
   type(mpas_mesh_type),   intent(in)    :: mpas_mesh
   real, intent(in)                      :: averaging_radius ! input is in km--radially average within this distance of each point
   type(esmf_fieldbundle), intent(in)    :: input_bundle
   type(esmf_fieldbundle), intent(inout) :: output_bundle

   integer                         :: rc, i, j, k, iCell, nfields, nz
   real                            :: averaging_radius_meters, area_sum
   character(len=500), allocatable :: field_names(:)
   type(ESMF_Field), allocatable   :: input_fields(:), target_fields(:)
   real*8, allocatable :: dummy1(:), dummy2(:,:)
   real(esmf_kind_r8), allocatable :: mean_1d(:), mean_2d(:,:)
   real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)
   integer, allocatable :: mask(:), num_points(:)
   integer, allocatable :: indices(:,:)

   ! First, figure out the needed geometry to radially average
   ! For each point, determine the number of points within the radius
   !  num_points(:) and then indices of those cells indices(:,:)
   ! This searching can be a bit slow for all meshes, but
   !  each processor can do its own number of cells so long as
   !  each processor has access to full geometry/mesh fields, which makes it faster

   averaging_radius_meters = averaging_radius * 1000.0 ! need meters
   mpas_mesh_radius = mpas_mesh%mpas_mesh_radius ! meters

   allocate(mask(mpas_mesh%nCells)) ! needs to be full field
   allocate(num_points(mpas_mesh%cell_start:mpas_mesh%cell_end))
   allocate(indices(mpas_mesh%cell_start:mpas_mesh%cell_end,num_neigh_max))

   if(localpet==0)write(*,*)'Radially averaging within ',averaging_radius_meters, 'meters'

   call cpu_time(start)
   do i = mpas_mesh%cell_start, mpas_mesh%cell_end
      mask(:) = UNMARKED ! reset for each cell on the full grid to get the unique neighborhood about each point
      num_points(i) = 0
      ! output is num_points(i) and indices(i,:)
      call find_neighborhoods(i, mpas_mesh%latCell(i), mpas_mesh%lonCell(i), averaging_radius_meters, .true., &
                                 INSIDE, mask, mpas_mesh%cellsOnCell, mpas_mesh%nEdgesOnCell, mpas_mesh%bdyMaskCell, &
                                 mpas_mesh%latCell, mpas_mesh%lonCell, num_points(i),indices(i,:))
   enddo
   call cpu_time(finish)
   if(localpet==0) print '("Time to get indices = ",f6.3," seconds.")',finish-start

   ! Now, loop through all fields and radially average,
   !  using the indices we have figured out,
   !  which are the same for all fields

   ! get the number of fields from the input bundle
   call bundle_get_num_fields(input_bundle,nfields)

   ! get the names of the fields from the bundle (just for printing purposes)
   allocate(field_names(nfields))
   call bundle_get_field_names(input_bundle,field_names)

   ! get the fields from the input bundle
   allocate(input_fields(nfields))
   call bundle_get_fields(input_bundle,input_fields)

   ! get the fields from the target bundle that we will fill
   allocate(target_fields(nfields))
   call bundle_get_fields(output_bundle,target_fields)

   do i = 1,nfields
   
      if (localpet==0) write(*,*)"- AVERAGING "//trim(adjustl(field_names(i)))

      if ( nVertLevelsPerVariable(i) == 1 ) then

         ! get the full field [input_fields(i)] on all processors
         allocate(dummy1(mpas_mesh%nCells))
         call gather_on_all_procs(localpet, input_fields(i), (/mpas_mesh%nCells/), dummy1)

         ! average. each processor does this and has the full field (dummy1)
         ! this could be memory intensive...we'll see
         allocate(mean_1d(mpas_mesh%cell_start:mpas_mesh%cell_end))
         mean_1d = 0.0
         do j = mpas_mesh%cell_start, mpas_mesh%cell_end
            area_sum = 0.0
            do k = 1, num_points(j)
               iCell = indices(j,k)
               mean_1d(j) = mean_1d(j) + dummy1(iCell)*mpas_mesh%areaCell(iCell) ! area-weighted average
               area_sum = area_sum + mpas_mesh%areaCell(iCell)
            end do
           !mean_1d(j) = mean_1d(j)/num_points(j) ! average over all points
           !if ( num_points(j) == 1 ) then ! testing
           !   mean_1d(j) = dummy1(j)
           !else
               mean_1d(j) = mean_1d(j)/area_sum      ! area-weighted average
           !endif
         end do

         ! now update the field
         call ESMF_FieldGet(target_fields(i), farrayPtr=varptr, rc=rc) ! varptr is decomposed across all processors
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         do j = 1, mpas_mesh%nCellsPerPET
            varptr(j) = mean_1d(mpas_mesh%elemIDs(j))
         enddo
         nullify(varptr)

         deallocate(mean_1d, dummy1)

      else ! 3d variable

         nz = nVertLevelsPerVariable(i)

         ! get the full field [input_fields(i)] on all processors
         allocate(dummy2(mpas_mesh%nCells,nz)) ! this is the dimension order things are assigned to the ESMF field
         call gather_on_all_procs(localpet, input_fields(i), (/mpas_mesh%nCells,nz/), dummy2)

         ! average. each processor does this and has the full field (dummy2)
         ! this could be memory intensive...we'll see
         allocate(mean_2d(mpas_mesh%cell_start:mpas_mesh%cell_end,nz))
         mean_2d = 0.0
         do j = mpas_mesh%cell_start, mpas_mesh%cell_end
            area_sum = 0.0
            do k = 1, num_points(j)
               iCell = indices(j,k)
               mean_2d(j,:) = mean_2d(j,:) + dummy2(iCell,:)*mpas_mesh%areaCell(iCell) ! area-weighted average
               area_sum = area_sum + mpas_mesh%areaCell(iCell)
            end do
           !mean_2d(j,:) = mean_2d(j,:)/num_points(j) ! average over all points
            mean_2d(j,:) = mean_2d(j,:)/area_sum      ! area-weighted average
         end do

         ! now update the field
         call ESMF_FieldGet(target_fields(i), farrayPtr=varptr2, rc=rc) ! varptr2 is decomposed across all processors
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         do j = 1, mpas_mesh%nCellsPerPET
            varptr2(j,:) = mean_2d(mpas_mesh%elemIDs(j),:)
         enddo
         nullify(varptr2)

         deallocate(mean_2d, dummy2)
      endif

   enddo

   ! update the fields in the output bundle
   call ESMF_FieldBundleAddReplace(output_bundle, target_fields, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN ESMF_FieldBundleAddReplace", rc)

   deallocate(mask,num_points,indices)
   deallocate(input_fields, target_fields, field_names)

end subroutine radially_average

recursive subroutine find_neighborhoods(my_cell, lat_center, lon_center, averaging_radius, central_cell, tval, mask, cellsOnCell, nEdgesOnCell, bdyMaskCell, latCell, lonCell, n, indices)
   ! Starting at some cell (my_cell), expand radially, marking all unmarked cells within a radius of my_cell as tval
   ! Save the number of points in the radius (n) and the cell center indices (indices)
   integer, intent(in) :: my_cell, tval
   logical, intent(in) :: central_cell ! true for the first time this function is called for a given cell
   real,    intent(in) :: lat_center, lon_center ! lat/lon of the central cell; always compute distances with respect to this
   real,    intent(in) :: averaging_radius ! meters
   integer, dimension(:,:), intent(in) :: cellsOnCell
   integer, dimension(:), intent(in) :: nEdgesOnCell
   integer, dimension(:), intent(in) :: bdyMaskCell
   real,    dimension(:), intent(in) :: latCell, lonCell
   integer, dimension(:), intent(inout) :: mask
   integer, intent(inout) :: n ! number of points within averaging_radius of my_cell (including my_cell itself)
   integer, intent(inout) :: indices(num_neigh_max) ! indices of points within averaging_radius of my_cell (including index of my_cell)
                                                    !could also do dimension(:) instead of declaring as num_neigh_max)
   integer :: i, j, iCell

   if ( central_cell ) then ! only true for the central cell
      mask(my_cell) = tval ! mark the central cell
      n = n + 1
      indices(n) = my_cell

      ! if the central cell is masked, return now and don't average for that cell; just use the point value
      if ( bdyMaskCell(my_cell) == mask_value ) return  ! get out of the subroutine
   endif

   do i=1, nEdgesOnCell(my_cell)
      iCell = cellsOnCell(i, my_cell)
      if (mask(iCell) == tval) cycle ! already marked, so don't mark again
      if (bdyMaskCell(iCell) == mask_value ) cycle ! don't include points we are masking out in the averaging
      if (sphere_distance(lat_center, lon_center, latCell(iCell), lonCell(iCell), mpas_mesh_radius) <= averaging_radius ) then
         mask(iCell) = tval
         n = n + 1
         indices(n) = iCell
         if ( n .gt. num_neigh_max ) then
            write(*,*)' there are > ',num_neigh_max,' points in the neighborhood'
            write(*,*)' n/num_max = ',n,num_neigh_max
            write(*,*)' change num_neigh_max and try again.'
            call error_handler("In find_neighborhoods", -3)
         endif
         call find_neighborhoods(iCell, lat_center, lon_center, averaging_radius, .false., tval, mask, cellsOnCell, nEdgesOnCell, bdyMaskCell, latCell, lonCell, n, indices)
      end if
   end do

end subroutine find_neighborhoods

subroutine apply_smoothing_filter(localpet,mpas_mesh,smoother_dimensionless_coefficient,bundle)
   integer, intent(in)                   :: localpet
   type(mpas_mesh_type),   intent(in)    :: mpas_mesh
   real,                   intent(in)    :: smoother_dimensionless_coefficient
   type(esmf_fieldbundle), intent(inout) :: bundle

   integer                         :: rc, i, j, k, iCell, nfields, nz
   integer                         :: iEdge, cell1, cell2
   real                            :: averaging_radius_meters, area_sum
   real                            :: smoother_coefficient, smoother_flux
   character(len=500), allocatable :: field_names(:)
   type(ESMF_Field), allocatable   :: fields(:)
   real(esmf_kind_r8), allocatable :: field1d(:), field2d(:,:)
   real*8, allocatable :: field_filtered1d(:), field_filtered2d(:,:)
   real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)

   ! get the number of fields from the bundle
   call bundle_get_num_fields(bundle,nfields)

   ! get the names of the fields from the bundle (just for printing purposes)
   allocate(field_names(nfields))
   call bundle_get_field_names(bundle,field_names)

   ! get the input fields
   allocate(fields(nfields))
   call bundle_get_fields(bundle,fields)

   do i = 1,nfields
   
      if (localpet==0) write(*,*)"- SMOOTHING "//trim(adjustl(field_names(i)))

      if ( nVertLevelsPerVariable(i) == 1 ) then
        !if (localpet==0) then
            allocate(field1d(mpas_mesh%nCells))
        !else
        !   allocate(field1d(0))
        !endif

         ! gather onto process 0
         call ESMF_FieldGather(fields(i), field1d, rootPet=0, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather", rc)

         ! do the smoothing on process 0
         allocate(field_filtered1d(mpas_mesh%nCells))
         if ( localpet == 0 ) then
            field_filtered1d(:) =  field1d(:)
            do iEdge=1,mpas_mesh%nEdges
               cell1 = mpas_mesh%cellsOnEdge(1,iEdge)
               cell2 = mpas_mesh%cellsOnEdge(2,iEdge)
               if (cell1 <= mpas_mesh%nCells .or. cell2 <= mpas_mesh%nCells) then
                  smoother_coefficient = 0.5*smoother_dimensionless_coefficient*(mpas_mesh%areaCell(cell1)+mpas_mesh%areaCell(cell2))
                  smoother_flux = smoother_coefficient*mpas_mesh%dvEdge(iEdge)*(field1d(cell2) - field1d(cell1))/mpas_mesh%dcEdge(iEdge)
                  field_filtered1d(cell1) = field_filtered1d(cell1) + smoother_flux/mpas_mesh%areaCell(cell1)
                  field_filtered1d(cell2) = field_filtered1d(cell2) - smoother_flux/mpas_mesh%areaCell(cell2)
               end if
            end do
         endif ! localpet == 0
         call mpi_bcast(field_filtered1d,mpas_mesh%nCells,mpi_real8,0,mpi_comm_world,rc)
         field1d(:) =  field_filtered1d(:) ! rename, field1d is esmf_kind_r8, which is needed

         ! now update the field
         call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc) ! varptr is decomposed across all processors
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         do j = 1, mpas_mesh%nCellsPerPET
            varptr(j) = field1d(mpas_mesh%elemIDs(j))
         enddo
         nullify(varptr)

         deallocate(field1d,field_filtered1d)

      else ! 3d variable

         nz = nVertLevelsPerVariable(i)

        !if (localpet==0) then
            allocate(field2d(mpas_mesh%nCells,nz)) ! this is the dimension order things are assigned to the ESMF field
        !else
        !   allocate(field2d(0,0))
        !endif

         ! gather onto process 0
         call ESMF_FieldGather(fields(i), field2d, rootPet=0, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather", rc)

         ! do the smoothing on process 0
         allocate(field_filtered2d(mpas_mesh%nCells,nz))
         if ( localpet == 0 ) then
            field_filtered2d(:,:) =  field2d(:,:)
            do iEdge=1,mpas_mesh%nEdges
               cell1 = mpas_mesh%cellsOnEdge(1,iEdge)
               cell2 = mpas_mesh%cellsOnEdge(2,iEdge)
               if (cell1 <= mpas_mesh%nCells .or. cell2 <= mpas_mesh%nCells) then
                  smoother_coefficient = 0.5*smoother_dimensionless_coefficient*(mpas_mesh%areaCell(cell1)+mpas_mesh%areaCell(cell2))
                  do k = 1,nz
                     smoother_flux = smoother_coefficient*mpas_mesh%dvEdge(iEdge)*(field2d(cell2,k) - field2d(cell1,k))/mpas_mesh%dcEdge(iEdge)
                     field_filtered2d(cell1,k) = field_filtered2d(cell1,k) + smoother_flux/mpas_mesh%areaCell(cell1)
                     field_filtered2d(cell2,k) = field_filtered2d(cell2,k) - smoother_flux/mpas_mesh%areaCell(cell2)
                  enddo
               end if
            end do
         endif ! localpet == 0
         call mpi_bcast(field_filtered2d,mpas_mesh%nCells*nz,mpi_real8,0,mpi_comm_world,rc)
         field2d(:,:) = field_filtered2d(:,:) ! rename, field2d is esmf_kind_r8, which is needed

         ! now update the field
         call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc) ! varptr2 is decomposed across all processors
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         do j = 1, mpas_mesh%nCellsPerPET
            varptr2(j,:) = field2d(mpas_mesh%elemIDs(j),:)
         enddo
         nullify(varptr2)

         deallocate(field2d,field_filtered2d)

      endif ! 2d or 3d variable

   enddo ! loop over variables

   ! update the fields in the bundle
   call ESMF_FieldBundleAddReplace(bundle, fields, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN ESMF_FieldBundleAddReplace", rc)

   deallocate(fields, field_names)

end subroutine apply_smoothing_filter

subroutine set_local_extrapolation_option(my_interp_method,my_extrap_method,extrapMethod)
   character(len=*), intent(in) :: my_interp_method
   character(len=*), intent(in) :: my_extrap_method
   type(ESMF_ExtrapMethod_Flag), intent(inout) :: extrapMethod

   if (trim(adjustl(my_extrap_method)).eq."esmf_creep") then
      extrapMethod = ESMF_EXTRAPMETHOD_CREEP
   else if (trim(adjustl(my_extrap_method)).eq."esmf_creep_nrst_d") then
      extrapMethod = ESMF_EXTRAPMETHOD_CREEP_NRST_D
   else if (trim(adjustl(my_extrap_method)).eq."esmf_nearest_d") then
      extrapMethod = ESMF_EXTRAPMETHOD_NEAREST_D
   else if (trim(adjustl(my_extrap_method)).eq."esmf_nearest_stod") then
      extrapMethod = ESMF_EXTRAPMETHOD_NEAREST_STOD
   else if (trim(adjustl(my_extrap_method)).eq."mpas_nearest") then ! we'll use our own, non-esmf extrapolation
      extrapMethod = ESMF_EXTRAPMETHOD_NONE
   else if (trim(adjustl(my_extrap_method)).eq."none") then
      extrapMethod = ESMF_EXTRAPMETHOD_NONE
   else
      call error_handler("Invalid extrapolation method = "//trim(adjustl(my_extrap_method)), -223)
   endif

end subroutine set_local_extrapolation_option

! make sure the namelist-provided 'extrap_method' is appropriate
! potentially updates namelist variable extrap_method
subroutine set_global_extrapolation(localpet)
   integer, intent(in) :: localpet
   if ( regional_mesh ) then
      if ( trim(adjustl(extrap_method)) .eq. "NULL" .or. &
           trim(adjustl(extrap_method)) .eq. "NONE" .or.  &
           trim(adjustl(extrap_method)) .eq. "none" ) then
         if (localpet==0) write(*,*)'You must choose an extrapolation option for a regional mesh.'
         call error_handler("Invalid extrapolation method = "//trim(adjustl(extrap_method)), -223)
      endif
      ! Can't use ESMF extrapolation with conservative interpolation
      if ( trim(adjustl(interp_method)) .eq. "conserve1" .or. &
           trim(adjustl(interp_method)) .eq. "conserve2" ) then
         extrap_method = 'mpas_nearest'
         if ( localpet == 0 ) write(*,*)'Resetting extrap_method = mpas_nearest b/c interp_method = '//trim(adjustl(interp_method))
      endif
      if ( localpet == 0 ) write(*,*)'Regional mesh: Using extrapolation method = '//trim(adjustl(extrap_method))
      if ( trim(adjustl(extrap_method)) .eq. "esmf_creep" ) then
         if (localpet == 0)write(*,*)'For esmf_creep, extrap_num_levels_creep = ',extrap_num_levels_creep
      endif
   else
      extrap_method = 'none'
      if ( localpet == 0 ) write(*,*)'Global mesh: No extrapolation needed.'
   endif
end subroutine set_global_extrapolation

end module interp
