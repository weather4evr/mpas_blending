module interp

use mpi
use esmf
use netcdf
use utils_mod, only : error_handler, sphere_distance
use model_grid, only : mpas_mesh_type, mask_value, INSIDE, UNMARKED
use input_data, only : nVertLevelsPerVariable

implicit none

private
 
! Public subroutines
public :: interp_data, radially_average, apply_smoothing_filter

! Variables private to this module
integer ::  isrctermprocessing = 1
integer, parameter :: num_neigh_max = 1000 ! max number of points in a neighborhood about a point
real,    parameter :: PI = 3.14159265359
real :: start, finish
real :: mpas_mesh_radius ! put this up here so all subroutines can see it
 
contains
 
subroutine interp_data(localpet,input_bundle,target_bundle,my_interp_method)
 
   integer, intent(in)                   :: localpet
   type(esmf_fieldbundle), intent(in)    :: input_bundle
   type(esmf_fieldbundle), intent(inout) :: target_bundle
   character(len=*), intent(in)          :: my_interp_method
   integer                          :: rc, i, nfields, ungriddedDimCount, k
   integer                          :: ungriddedUBound(10), ungriddedLBound(10)
   character(len=500), allocatable    :: field_names(:)
   type(ESMF_RegridMethod_Flag)     :: method
   type(ESMF_RouteHandle)           :: my_rh
   type(ESMF_Field), allocatable    :: fields_input_grid(:), fields_target_grid(:)

   real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)
   type(ESMF_Mesh) :: mesh
   integer :: elementCount
   integer, allocatable :: elementMask(:)
   type(ESMF_GeomType_Flag) :: geomtype
    
   !if ( ESMF_RouteHandleIsCreated(routehandle, rc=rc)) then   ... 

   if(localpet==0)write(*,*)'In interp_data; using method '//trim(my_interp_method)

   ! get the number of fields from the bundle
   call ESMF_FieldBundleGet(input_bundle, fieldCount=nfields, &
                             itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldCount", rc)

   ! get the names of the fields from the bundle (just for printing purposes)
   allocate(field_names(nfields))
   call ESMF_FieldBundleGet(input_bundle, fieldNameList=field_names, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldNameList", rc)

   ! get the input and target fields
   allocate(fields_input_grid(nfields))
   allocate(fields_target_grid(nfields))
   call ESMF_FieldBundleGet(input_bundle, fieldList=fields_input_grid, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet", rc)

   call ESMF_FieldBundleGet(target_bundle, fieldList=fields_target_grid, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet", rc)

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

   !if(localpet==0)write(*,*)'making routehandle for '//trim(my_interp_method)//' regridding'

   ! make the route-handle for interpolation
   call make_rh(fields_input_grid(1), fields_target_grid(1),method, my_rh)

   ! this will make  copy of a route handle (Not used)
   !my_rh = ESMF_RouteHandleCreate(rh_bilin, rc=rc)
   
   do i = 1,nfields
   
     !if (localpet==0) write(*,*)"- INTERPOLATING "//trim(adjustl(field_names(i)))//" with method "//trim(adjustl(my_interp_method))

      ! testing about getting number of ungridded dimensions...can probably remove someday
     !call ESMF_FieldGet(fields_input_grid(i),ungriddedDimCount=ungriddedDimCount, &
     !           ungriddedLBound=ungriddedLBound, ungriddedUBound=ungriddedUBound ,rc=rc)
     !if (localpet==0) write(*,*)'ungriddedDimCount = ',ungriddedDimCount
     !if (localpet==0) write(*,*)'ungriddedLBound = ',ungriddedLBound
     !if (localpet==0) write(*,*)'ungriddedUBound = ',ungriddedUBound
      ! end test

      ! do the interpolation
      call ESMF_FieldRegrid(fields_input_grid(i), fields_target_grid(i), my_rh, rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN FieldRegrid", rc)

      ! kind of a test...where masked, force the interpolated field to missing
      !  .... doesn't seem to be working ... 
      call ESMF_FieldBundleGet(target_bundle, geomtype=geomtype)
      if ( geomtype == ESMF_GEOMTYPE_MESH ) then ! only do if it's a mesh-to-mesh interpolation.
     !if ( 1 == 2 ) then
         if(localpet==0) write(*,*)'resetting interpolated values outside mask to mask_value'
         call ESMF_FieldBundleGet(target_bundle, mesh=mesh)
         call ESMF_MeshGet(mesh, elementCount=elementCount)
         allocate(elementMask(elementCount))
         call ESMF_MeshGet(mesh, elementMask=elementMask)
         if ( nVertLevelsPerVariable(i) == 1 ) then
            call ESMF_FieldGet(fields_target_grid(i), farrayPtr=varptr, rc=rc) ! varptr is decomposed across all processors
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
               call error_handler("IN FieldGet", rc)
            where( elementMask == mask_value ) varptr = real(mask_value,esmf_kind_r8)
            nullify(varptr)
         else
            call ESMF_FieldGet(fields_target_grid(i), farrayPtr=varptr2, rc=rc) ! varptr2 is decomposed across all processors
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
               call error_handler("IN FieldGet", rc)
            do k = 1,nVertLevelsPerVariable(i)
               where( elementMask == mask_value ) varptr2(:,k) = real(mask_value,esmf_kind_r8)
            enddo
            nullify(varptr2)
         endif
         deallocate(elementMask)
      endif
      ! end test

   enddo

   ! update the fields in target_bundle
   call ESMF_FieldBundleAddReplace(target_bundle, fields_target_grid, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN ESMF_FieldBundleAddReplace", rc)

   ! clean-up
   deallocate(fields_input_grid,fields_target_grid)
   deallocate(field_names)

   ! get rid of the route handle
   call ESMF_FieldRegridRelease(routehandle=my_rh, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridRelease", rc)

end subroutine interp_data
                
subroutine make_rh(input_field,target_field,method,rh)

   type(ESMF_Field), intent(in)              :: input_field
   type(ESMF_Field), intent(inout)           :: target_field
   type(ESMF_RegridMethod_Flag), intent(in)  :: method
   type(ESMF_RouteHandle), intent(inout)        :: rh

   integer :: rc

   call ESMF_FieldRegridStore(input_field, target_field, &
                               regridmethod=method, &
                               routehandle=rh, &
                               srcTermProcessing=isrctermprocessing, &
                               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                               srcMaskValues=(/mask_value/), &
                               dstMaskValues=(/mask_value/), &
                               rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldRegridStore", rc)
   return 

end subroutine make_rh

subroutine radially_average(localpet,mpas_mesh,averaging_radius,bundle)
 
   integer, intent(in)                   :: localpet
   type(mpas_mesh_type),   intent(in)    :: mpas_mesh
   real, intent(in)                      :: averaging_radius ! input is in km--radially average within this distance of each point
   type(esmf_fieldbundle), intent(inout) :: bundle

   integer                         :: rc, i, j, k, iCell, nfields, nz
   real                            :: averaging_radius_meters, area_sum
   character(len=500), allocatable :: field_names(:)
   type(ESMF_Field), allocatable   :: fields(:)
   real*8, allocatable :: dummy1(:), dummy2(:,:)
   real(esmf_kind_r8), allocatable :: dum1d(:), dum2d(:,:)
   real(esmf_kind_r8), allocatable :: mean_1d(:), mean_2d(:,:)
   real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)
   integer, allocatable :: mask(:), num_points(:)
   integer, allocatable :: indices(:,:)

   ! First, figure out the needed geometry to radially average
   ! For each point, determine the number of points within the radius
   !  num_points(:) and then indices of those cells indices(:,:)
   ! This searching can be a bit slow for all meshes, but
   !  each processor can do its own number of cells so long as
   !  each processor has access to full geometry/mesh fields, which
   !  makes it faster

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

   ! get the number of fields from the bundle
   call ESMF_FieldBundleGet(bundle, fieldCount=nfields, &
                             itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldCount", rc)

   ! get the names of the fields from the bundle (just for printing purposes)
   allocate(field_names(nfields))
   call ESMF_FieldBundleGet(bundle, fieldNameList=field_names, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldNameList", rc)

   ! get the input fields
   allocate(fields(nfields))
   call ESMF_FieldBundleGet(bundle, fieldList=fields, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet", rc)

   do i = 1,nfields
   
      if (localpet==0) write(*,*)"- AVERAGING "//trim(adjustl(field_names(i)))

      if ( nVertLevelsPerVariable(i) == 1 ) then
         if (localpet==0) then
            allocate(dum1d(mpas_mesh%nCells))
         else
            allocate(dum1d(0))
         endif

         ! gather onto process 0
         call ESMF_FieldGather(fields(i), dum1d, rootPet=0, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather", rc)

         if (localpet==0) write(*,*)'min/max0 = ',minval(dum1d),maxval(dum1d)

         ! dum1d is real(esmf_kind_r8). We need to get dum1d onto all processors
         ! (currently it's just on process 0)
         ! So, using mpi_bcast is the easiest way to do it, but mpi_bcast doesn't
         ! have a real "type" of esmf_kind_r8.  So we need to use a different
         ! variable to do the broadcasting that is mpi_real, hence, introducing dummy2
         allocate(dummy1(mpas_mesh%nCells))
         if (localpet==0) dummy1 = dum1d

         call mpi_bcast(dummy1,mpas_mesh%nCells,mpi_real8,0,mpi_comm_world,rc)
        !write(*,*)'pe/min/max = ',localpet,minval(dummy1),maxval(dummy1)

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
            mean_1d(j) = mean_1d(j)/area_sum      ! area-weighted average
         end do

         ! now update the field
         call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc) ! varptr is decomposed across all processors
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         do j = 1, mpas_mesh%nCellsPerPET
            varptr(j) = mean_1d(mpas_mesh%elemIDs(j))
         enddo
         nullify(varptr)

         deallocate(dum1d, mean_1d)
         deallocate(dummy1)

      else ! 3d variable

         nz = nVertLevelsPerVariable(i)

         if (localpet==0) then
            allocate(dum2d(mpas_mesh%nCells,nz)) ! this is the dimension order things are assigned to the ESMF field
         else
            allocate(dum2d(0,0))
         endif

         ! gather onto process 0
         call ESMF_FieldGather(fields(i), dum2d, rootPet=0, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather", rc)

         if (localpet==0) write(*,*)'min/max0 = ',minval(dum2d),maxval(dum2d)

         ! dum2d is real(esmf_kind_r8). We need to get dum2d onto all processors
         ! (currently it's just on process 0)
         ! So, using mpi_bcast is the easiest way to do it, but mpi_bcast doesn't
         ! have a real "type" of esmf_kind_r8.  So we need to use a different
         ! variable to do the broadcasting that is mpi_real, hence, introducing dummy2
         allocate(dummy2(mpas_mesh%nCells,nz)) ! this is the dimension order things are assigned to the ESMF field
         if (localpet==0) dummy2 = dum2d

         call mpi_bcast(dummy2,mpas_mesh%nCells*nz,mpi_real8,0,mpi_comm_world,rc)
        !write(*,*)'pe/min/max = ',localpet,minval(dummy2),maxval(dummy2)

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
         call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc) ! varptr2 is decomposed across all processors
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
            call error_handler("IN FieldGet", rc)
         do j = 1, mpas_mesh%nCellsPerPET
            varptr2(j,:) = mean_2d(mpas_mesh%elemIDs(j),:)
         enddo
         nullify(varptr2)

         deallocate(dum2d, mean_2d)
         deallocate(dummy2)
      endif

   enddo

   ! update the fields in the bundle
   call ESMF_FieldBundleAddReplace(bundle, fields, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN ESMF_FieldBundleAddReplace", rc)

   deallocate(mask,num_points,indices)
   deallocate(fields, field_names)

end subroutine radially_average

recursive subroutine find_neighborhoods(my_cell, lat_center, lon_center, averaging_radius, central_cell, tval, mask, cellsOnCell, nEdgesOnCell, bdyMaskCell, latCell, lonCell,n,indices)
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
      if (bdyMaskCell(iCell) == mask_value ) cycle ! don't include points we are masking out in the averaging
      if (mask(iCell) == tval) cycle
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
         call find_neighborhoods(iCell, lat_center, lon_center, averaging_radius, .false., tval, mask, cellsOnCell, nEdgesOnCell, bdyMaskCell, latCell, lonCell,n,indices)
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
   call ESMF_FieldBundleGet(bundle, fieldCount=nfields, &
                             itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldCount", rc)

   ! get the names of the fields from the bundle (just for printing purposes)
   allocate(field_names(nfields))
   call ESMF_FieldBundleGet(bundle, fieldNameList=field_names, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet fieldNameList", rc)

   ! get the input fields
   allocate(fields(nfields))
   call ESMF_FieldBundleGet(bundle, fieldList=fields, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet", rc)

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

end module interp
