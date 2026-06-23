!> Building the state vector
!!
!! This module provides variables & routines for
!! defining the state vector.
!!
!! The module contains three routines
!! - **init_id** - initialize the array `id`
!! - **init_sfields** - initialize the array `sfields`
!! - **setup_statevector** - generic routine controlling the initialization
!!
!! The declarations of **id** and **sfields** as well as the
!! routines **init_id** and **init_sfields** might need to be
!! adapted to a particular modeling case. However, for most
!! parts also the configruation using the namelist is possible.
!!
!! __Revision history__
!! ~2022 Frauke - initial functionality distributed over different routines
!! 2025-12 - Lars Nerger - restructuring code introducing module statevector_pdaf
!! 2026-04 - Anna Broschke - restructing code intoducing varible state vector setup
module statevector_pdaf


  implicit none
  save

  !---- `field_ids` and `state_field` can be adapted for a DA case -----

  ! Declare Fortran type holding the indices of model fields in the state vector
  ! This can be extended to any number of fields - it serves to give each field a name
  ! This is grouped but the groups can be changed
  type field_ids
  ! Physics
     integer :: ssh = 0         ! Sea Surface Hight
     integer :: u = 0           ! Zonal velocity (interpolated on nodes)
     integer :: v = 0           ! Meridional velocity (interpolated on nodes)
     integer :: w = 0           ! Vertical velocity
     integer :: temp = 0        ! Themperature
     integer :: salt = 0        ! Salinity
     integer :: MLD1 = 0        ! Mixed Layer Depth !diagnostic
     integer :: MLD2 = 0        ! Mixed Layer Depth !diagnostic
  ! ice
     integer :: a_ice = 0       ! Sea-ice concentration
  ! Carbon
     integer :: DIC = 0         ! Dissolved inorganic carbon
     integer :: DOC = 0         ! Dissolved organic carbon
     integer :: Alk = 0         ! Alkalinity
     integer :: pCO2s = 0       ! Partial pressure CO2 surface ocean ! diagnostic
     integer :: CO2f = 0        ! CO2 flux from atmosphere into ocean ! diagnostic
     integer :: alphaCO2 = 0    ! solubility of surface CO2 ! diagnostic
     integer :: PistonVel = 0   ! air-sea piston velocity ! diagnostic
  ! Nutrints
     integer :: DIN = 0         ! Dissolved inorganic nitrogen
     integer :: DSi = 0         ! Dissolved inorganic Silicate
     integer :: Fe = 0          ! Iron
  ! Phytoplankton
     ! small Phytoplankton
     integer :: PhyChl = 0      ! intracell chlorophyll small phytoplankton
     integer :: PhyN = 0        ! intracell nitrogen small phytoplankton
     integer :: PhyC = 0        ! intracell carbon small phytoplankton
     integer :: PhyCalc = 0     ! intracell carbonate small phytoplankton
     ! Diatoms
     integer :: DiaChl = 0      ! intracell chlorophyll diatoms
     integer :: DiaN = 0        ! intracell nitrogen diatoms
     integer :: DiaC = 0        !intracell carbon diatoms
     integer :: DiaSi = 0       ! intracell silicate diatoms
     ! ToDo insert 4p Phytoplankton
     !integer :: CoccoN    ! Coccos
     !integer :: CoccoC
     !integer :: CoccoChl
     !integer :: PhaeoN     !Phaeocystis
     !integer :: PhaeoC
     !integer :: PhaeoChl
  ! Zooplankton
     integer :: Zo1C = 0        ! carbon in small zooplankton
     integer :: Zo1N = 0        ! nitrogen in small zooplankton
     integer :: Zo2C = 0        ! carbon in macrozooplankton
     integer :: Zo2N = 0        ! Nitrogen in macrozooplankton
     integer :: Zo3C = 0        ! Microzooplankton carbon
     integer :: Zo3N = 0        ! Microzooplankton nitrogen
  ! Detritus
     integer :: DetC = 0        ! carbon in small detritus
     integer :: DetCalc = 0     ! calcite in small detritus
     integer :: DetSi = 0       ! silicate in small detritus
     integer :: DetN = 0        ! nitrogen in small detritus
     integer :: Det2C = 0       ! carbon in large detritus
     integer :: Det2Calc = 0    ! calcite in large detritus
     integer :: Det2Si = 0      ! Silicate in large detritus
     integer :: Det2N = 0       ! nitrogen in large detritus
  ! Other
     integer :: DON = 0         ! Disolves organic carbon
     integer :: O2 = 0          ! Oxygen
     integer :: PAR = 0         ! photosynthetically active radiation
     integer :: sigma = 0       ! potential density
  ! Diagnostics
     integer :: NPPn = 0        ! mean net primary production small phytoplankton
     integer :: NPPd = 0        ! mean net primary production diatoms
     integer :: export = 0      ! export through particle sinking at 190m
     !INTEGER :: TChl = 0       ! Total chlorophyll = PhyChl + DiaChl
     !INTEGER :: TDN = 0        ! Total dissolved N = DIN + DON
     !INTEGER :: TOC = 0        ! Total organic carbon: PhyC + DiaC + DetC + DOC + HetC
  ! Reflectance
     integer, allocatable  :: Edz3D (:)        ! Downwelling direct stream of light 
     integer, allocatable  :: Esz3D (:)        ! Downwelling diffuse stream of light
     integer, allocatable  :: Euz3D (:)        ! Upwelling stream of light
     integer, allocatable  :: Eutop3D (:)      ! Direct stream of light on top of layer
     !integer, allocatable  :: Estop3D (:,:,:)      ! Diffuse stream of light on top of layer

!     INTEGER :: TChl   ! Total chlorophyll = PhyChl + DiaChl
!     INTEGER :: TDN    ! Total dissolved N = DIN + DON
!     INTEGER :: TOC    ! Total organic carbon: PhyC + DiaC + DetC + DOC + HetC
  end type field_ids

  ! Declare Fortran type holding the definitions for model fields
  type state_field
     integer :: ndims = 0                  !< Number of field dimensions (2 or 3)
     integer :: dim = 0                    !< Dimension of the field
     integer :: off = 0                    !< Offset of field in state vector
     logical :: nz1 = .true.               !< Vertical coordinates (on levels / on layers)
     character(len=10) :: variable = ''    !< Name of field
     character(len=50) :: long_name = ''   !< Long name of field
     character(len=20) :: units = ''       !< Unit of variable
     integer :: varid(9)                   !< To write to netCDF file
     logical :: updated = .true.           !< Whether variable is updated through assimilation
     logical :: output(8,3) = .false.      !< How frequently output is written
     logical :: bgc = .false.              !< Whether variable is biogeochemistry (or physics)
     integer :: trnumfesom = -1            !< Tracer index in FESOM-REcoM
     integer :: tridfesom = -1             !< Tracer ID in FESOM-REcoM
     integer :: id_tr = -1                 !< Field index in list of 3D model tracer fields
  end type state_field


  ! Declare Fortran type holding the definitions for local model fields
  ! This is separate from state_field to support OpenMP
  type state_field_l
     integer :: dim = 0                    !< Dimension of the field
     integer :: off = 0                    !< Offset of field in state vector
  end type state_field_l

  integer :: phymin, phymax   ! First and last physics field in state vector
  integer :: bgcmin, bgcmax   ! First and last biogeochemistry field in state vector

  ! Variables to activate a field from the namelist
  
  logical :: sv_physics = .false.
  logical :: sv_ice = .false.
  logical :: sv_carbon = .false.
  logical :: sv_nutrients = .false.
  logical :: sv_phytoplankton = .false.
  logical :: sv_zooplankton = .false.
  logical :: sv_detritus = .false.
  logical :: sv_other = .false.
  logical :: sv_diagnostics = .false.
  logical :: sv_reflectance = .false.

  !---- The next variables usually do not need editing -----

  ! Type variable holding field IDs in state vector
  type(field_ids) :: id

  ! Type variable holding the definitions of model fields
  type(state_field), allocatable :: sfields(:)

  ! Type variable holding the definitions of local model fields
  ! This is separate from sfields to support OpenMP
  type(state_field_l), allocatable :: sfields_l(:)

!$OMP THREADPRIVATE(sfields_l)

  ! Variables to handle multiple fields in the state vector
  integer :: nfields           !< number of fields in state vector

contains


! ===================================================================================

!> Calculate the dimension of the process-local statevector.
!!
!! This routine is generic. case-specific adaptions should only
!! by done in the routines init_id and init_sfields.
!!
  subroutine setup_statevector(dim_state, dim_state_p, screen)

    use parallel_pdaf_mod, &
         only: mype=>mype_ens, npes=>npes_ens, task_id, comm_ensemble, &
         comm_model, MPI_SUM, MPI_INTEGER, MPIerr

    implicit none

! *** Arguments ***
    integer, intent(out) :: dim_state    !< Global dimension of state vector
    integer, intent(out) :: dim_state_p  !< Local dimension of state vector
    integer, intent(in)  :: screen       !< Verbosity flag

! *** Local variables ***
    integer :: i                 ! Counters


! ***********************************
! *** Initialize the state vector ***
! ***********************************

! *** Initialize array `id` ***

    call init_id(nfields)

! *** Initialize array `sfields` ***

    call init_sfields()
    call set_field_types(screen)

! *** Set state vector dimension ***

    dim_state_p = sum(sfields(:)%dim)

! *** Write information about the state vector ***

    if (mype==0) then
       write (*,'(/a,2x,a)') 'FESOM-PDAF', '*** Setup of state vector ***'
       write (*,'(a,5x,a,i5)') 'FESOM-PDAF', '--- Number of fields in state vector:', nfields
       write (*,'(a,a4,3x,a2,2x,a8,4x,a5,6x,a3,7x,a6,4x,a6,2x,a3,1x,a6)') &
            'FESOM-PDAF','pe','ID', 'variable', 'ndims', 'dim', 'offset', 'update', 'BGC', 'tracer'
    end if

    if (mype==0 .or. (task_id==1 .and. screen>2)) then
       do i = 1, nfields
          write (*,'(a, i4, i5,3x,a10,2x,i3,2x,i10,3x,i10,4x,l,4x,l,2x,i4)') 'FESOM-PDAF', &
               mype, i, sfields(i)%variable, sfields(i)%ndims, sfields(i)%dim, sfields(i)%off, sfields(i)%updated, &
               sfields(i)%bgc, sfields(i)%trnumfesom
       end do
    end if

    if (npes==1) then
       write (*,'(a,2x,a,1x,i10)') 'FESOM-PDAF', 'Full state dimension: ',dim_state_p
    else
       if (task_id==1) then
          if (screen>2 .or. mype==0) &
               write (*,'(a,2x,a,1x,i4,2x,a,1x,i10)') &
               'FESOM-PDAF', 'PE', mype, 'PE-local full state dimension: ',dim_state_p

          call MPI_Reduce(dim_state_p, dim_state, 1, MPI_INTEGER, MPI_SUM, 0, COMM_model, MPIerr)
          if (mype==0) then
             write (*,'(a,2x,a,1x,i10)') 'FESOM-PDAF', 'Global state dimension: ',dim_state
          end if
       end if
    end if
    call MPI_Barrier(comm_ensemble, MPIerr)

  end subroutine setup_statevector


!> Calculate the dimension of the process-local statevector.
!> This routine initializes the array `id`
!!
  subroutine init_id(nfields)

    use fesom_pdaf, only: tlam
    use assim_pdaf_mod, &
         only: nmlfile
    implicit none

! *** Arguments ***
    integer, intent(out) :: nfields
    integer :: cnt, i, l
    namelist /state_vector/ sv_physics, sv_ice, sv_carbon, sv_nutrients, &
            sv_phytoplankton, sv_zooplankton, sv_detritus, sv_other, &
            sv_diagnostics, sv_reflectance

#ifdef RECOM_WAVEBANDS
  !  allocate( Edz3D (nl -1,node_size,tlam))  
  !  allocate( Esz3D (nl -1,node_size,tlam))
  !  allocate( Euz3D (nl -1,node_size,tlam))
  !  allocate( Eutop3D (nl -1,node_size,tlam))
#endif

    open  (20,file=nmlfile)
    read  (20,NML=state_vector)
    close (20)

! Set field IDs
    cnt = 0
    if (sv_physics) then
            cnt = cnt +1
            id%ssh    =  cnt ! sea surface height
            cnt = cnt +1
            id%u      =  cnt ! zonal velocity
            cnt = cnt +1
            id%v      =  cnt ! meridional velocity
            cnt = cnt +1
            id%w      =  cnt ! vertical velocity
            cnt = cnt +1
            id%temp   =  cnt ! temperature
            cnt = cnt +1
            id%salt   =  cnt ! salinity
     end if

     if (sv_ice) then
             cnt = cnt +1
             id%a_ice  =  cnt ! sea-ice concentration
     end if

     if (sv_carbon) then
             cnt = cnt +1
             id%DIC    = cnt ! dissolved tracers
             cnt = cnt +1
             id%DOC    = cnt
             cnt = cnt +1
             id%Alk    = cnt
             cnt = cnt +1
             id%pCO2s = cnt ! surface carbon diags
             cnt = cnt +1
             id%CO2f   = cnt
             cnt = cnt +1
             id%alphaCO2  = cnt
             cnt = cnt +1
             id%PistonVel = cnt
      end if

      if (sv_nutrients) then
              cnt = cnt +1
              id%DIN    = cnt
              cnt = cnt +1
              id%DSi    =cnt
              cnt = cnt +1
              id%Fe     =cnt
      end if

      if (sv_phytoplankton) then
              ! small phytoplankton
              cnt = cnt +1
              id%PhyN   = cnt
              cnt = cnt +1
              id%PhyC   = cnt
              cnt = cnt +1
              id%PhyCalc= cnt
              cnt = cnt +1
              id%PhyChl = cnt ! chlorophyll-a small phytoplankton
              cnt = cnt +1
              ! Diatoms
              id%DiaN   = cnt ! diatoms
              cnt = cnt +1
              id%DiaC   = cnt
              cnt = cnt +1
              id%DiaSi  = cnt
              cnt = cnt +1
              id%DiaChl = cnt ! chlorophyll-a diatomsi

              !id%CoccoN = 45 !Coccos
              !id%CoccoC = 46
              !id%CoccoChl = 47
 
              !id%PhaeoN = 48 !Phaeocystis
              !id%PhaeoC = 49
              !id%PhaeoChl = 50
       end if

       if (sv_zooplankton) then
               cnt = cnt +1
               id%Zo1N   = cnt ! zooplankton
               cnt = cnt +1
               id%Zo1C   = cnt
               cnt = cnt +1
               id%Zo2N   = cnt
               cnt = cnt +1
               id%Zo2C   = cnt
               cnt = cnt +1
               id%Zo3N   = cnt
               cnt = cnt +1
               id%Zo3C   = cnt 
       end if

       if (sv_detritus) then
               cnt = cnt +1
               id%DetC      = cnt ! detritus
               cnt = cnt +1
               id%DetCalc   = cnt
               cnt = cnt +1
               id%DetSi     = cnt
               cnt = cnt +1
               id%DetN      = cnt
               cnt = cnt +1
               id%Det2C     = cnt
               cnt = cnt +1
               id%Det2Calc  = cnt
               cnt = cnt +1
               id%Det2Si    = cnt
               cnt = cnt +1
               id%Det2N     = cnt
       end if

       if (sv_other) then 
               cnt = cnt +1
               id%DON    = cnt
               cnt = cnt +1
               id%O2     = cnt
               cnt = cnt +1
               id%PAR    = cnt
               cnt = cnt +1
               id%sigma  = cnt
       end if

       if (sv_diagnostics) then
               cnt = cnt +1
               id%NPPn   = cnt
               cnt = cnt +1
               id%NPPd   = cnt
               cnt = cnt +1
               id%export = cnt
       end if

       if (sv_reflectance) then
               !allocate spectral varibles
              allocate(id%Edz3D(tlam))
              allocate(id%Esz3D(tlam))
              allocate(id%Euz3D(tlam))
              allocate(id%Eutop3D(tlam))

               do l = 1, tlam
                   cnt = cnt +1
                   id%Edz3D(l) = cnt
               end do                   !a lot of loops so the order is Edz3D for all wavelegnth and then next varible for all waavelegth
               do l = 1, tlam
                   cnt = cnt +1
                   id%Esz3D(l) = cnt
               end do
               do l = 1, tlam
                   cnt = cnt +1
                   id%Euz3D(l) = cnt
               end do
               do l = 1, tlam
                   cnt = cnt +1
                   id%Eutop3D(l) = cnt
               end do
       end if



! Total number of fields
    nfields = cnt

! physics part of state vector, specify start and end:
    !phymin = 1
    !phymax = 10
  
! BGC part of state vector, specify start and end:
    !bgcmin = 11
    !bgcmax = nfields


  end subroutine init_id
! ===================================================================================

!> This initializes the array sfields
!!
!! This routine initializes the sfields array with specifications
!! of the fields in the state vector.
!! 
!! ndims        - 1 for 2dim field (only surface of ocean) 
!!                2 for 3dim fiels (surface and depth)
!! nz1          - logical if z values ore one shorter like w
!! varibele     - short name for varible like in the code
!! long_name    - long name of varible
!! units        - units in latex code
!! updated      - will be updated during assimilation read in from namelist
!! bgc          - varible is a biogeochemestry varible
!! trnumfesom   - tracer index from FESOM-REcoM
!! tridfesom    - tracer ID from FESOM-REcoM
!! 
!!
  subroutine init_sfields()

    use fesom_pdaf, &
         only: myDim_nod2D, nlmax, tlam
    use assim_pdaf_mod, &
         only: nmlfile
    use parallel_pdaf_mod, &
         only: mype_world

    implicit none

! *** Local variables ***
    integer :: i, cnt           ! Counter
    integer :: id_var           ! varible for id number of varible which fills sfields 
    character(len=10), dimension(tlam) :: lams  !hard coded for wavelength names in spectral varibles

    ! logical varibles if varible should be updated defined in namelist
    logical :: upd_ssh = .false. 
    logical :: upd_u  = .false.
    logical :: upd_v = .false.
    logical :: upd_w = .false.
    logical :: upd_temp = .false.
    logical :: upd_salt = .false.

    logical :: upd_ice = .false.

    logical :: upd_DIC = .false.
    logical :: upd_DOC = .false.
    logical :: upd_Alk = .false.
    logical :: upd_pCO2s = .false.
    logical :: upd_CO2f = .false.
    logical :: upd_alphaCO2 = .false.
    logical :: upd_PistonVel = .false.

    logical :: upd_DIN = .false.
    logical :: upd_DSi =.false.
    logical :: upd_Fe = .false.

    logical :: upd_PhyCalc = .false.
    logical :: upd_PhyC = .false.
    logical :: upd_PhyN = .false.
    logical :: upd_PhyChl = .false.
    logical :: upd_DiaN = .false.
    logical :: upd_DiaC = .false.
    logical :: upd_DiaSi = .false.
    logical :: upd_DiaChl = .false.

    logical :: upd_Zo1C = .false.
    logical :: upd_Zo1N = .false.
    logical :: upd_Zo2C = .false.
    logical :: upd_Zo2N = .false.
    logical :: upd_Zo3C = .false.
    logical :: upd_Zo3N = .false.

    logical :: upd_DetC = .false.
    logical :: upd_DetCalc = .false.
    logical :: upd_DetSi = .false.
    logical :: upd_DetN = .false.
    logical :: upd_Det2C = .false.
    logical :: upd_Det2N = .false.
    logical :: upd_Det2Si = .false.
    logical :: upd_Det2Calc = .false.

    logical :: upd_DON = .false.
    logical :: upd_O2 = .false.
    logical :: upd_PAR = .false.
    logical :: upd_sigma = .false.

    logical :: upd_NPPn = .false.
    logical :: upd_NPPd = .false.
    logical :: upd_export = .false.

    logical :: upd_Edz3D = .false.
    logical :: upd_Esz3D = .false.
    logical :: upd_Euz3D = .false.
    logical :: upd_Eutop3D = .false.

! *** Allocate ***

    allocate(sfields(nfields))

! *** Read namelist file ***
    if (mype_world==0) write(*,*) 'Read namelist file for updated variables: ',nmlfile

    namelist /updated/ &
         upd_ssh, upd_u, upd_v, upd_w, upd_temp, upd_salt,&     ! Physic 
         upd_ice, &                                             ! ice
         upd_DIC, upd_DOC, upd_Alk,upd_pCO2s, upd_CO2f, upd_alphaCO2, upd_PistonVel, &  ! carbon
         upd_DIN, upd_DSi, upd_Fe, &                            ! nutrients
         upd_PhyCalc, upd_PhyC, upd_PhyN, upd_PhyChl, &         ! small phyto
         upd_DiaN, upd_DiaC, upd_DiaSi, upd_DiaChl, &           ! diatoms
         upd_Zo1C, upd_Zo1N, upd_Zo2C, upd_Zo2N, upd_Zo3C, upd_Zo3N, &              ! zooplankton
         upd_DetC, upd_DetCalc, upd_DetSi, upd_DetN     , &     ! small det
         upd_Det2C, upd_Det2N, upd_Det2Si, upd_Det2Calc , &     ! large det
         upd_DON, upd_O2, upd_PAR, upd_sigma, &                 ! other
         upd_NPPn, upd_NPPd,upd_export, &                       ! diagnostics
         upd_Edz3D, upd_Esz3D, upd_Euz3D, upd_Eutop3D           ! reflectance

    open  (20,file=nmlfile)
    read  (20,NML=updated)
    close (20)

! *** Define spectral bands for naming sfields correctly ***
     lams = [character(len=10) :: "400","425","450","475","500","525","550","575","600","625","650","675","700"]


! ****************
! *** Physics ****
! ****************

! SSH
        id_var = id%ssh
        if (id_var > 0) then
                sfields(id_var)%ndims = 1
                sfields(id_var)%variable = 'SSH'
                sfields(id_var)%long_name = 'Sea surface height'
                sfields(id_var)%units = 'm'
                sfields(id_var)%updated = upd_ssh
                sfields(id_var)%bgc = .false.
        endif

! u
        id_var = id%u
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'u'
                sfields(id_var)%long_name = 'Zonal velocity (interpolated on nodes)'
                sfields(id_var)%units = 'm/s'
                sfields(id_var)%updated = upd_u
                sfields(id_var)%bgc = .false.
        endif

! v
        id_var = id%v
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'v'
                sfields(id_var)%long_name = 'Meridional velocity (interpolated on nodes)'
                sfields(id_var)%units = 'm/s'
                sfields(id_var)%updated = upd_v
                sfields(id_var)%bgc = .false.
        endif

! w
        id_var = id%w
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .false.
                sfields(id_var)%variable = 'w'
                sfields(id_var)%long_name = 'Vertical velocity'
                sfields(id_var)%units = 'm/s'
                sfields(id_var)%updated = upd_w
                sfields(id_var)%bgc = .false.
        endif

! temp
        id_var = id%temp
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'T'
                sfields(id_var)%long_name = 'Temperature'
                sfields(id_var)%units = 'degC'
                sfields(id_var)%updated = upd_temp
                sfields(id_var)%bgc = .false.
                sfields(id_var)%trnumfesom = 1  
                sfields(id_var)%tridfesom = 0 
        endif

! salt
        id_var = id%salt
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'S'
                sfields(id_var)%long_name = 'Salinity'
                sfields(id_var)%units = 'psu'
                sfields(id_var)%updated = upd_salt
                sfields(id_var)%bgc = .false.
                sfields(id_var)%trnumfesom = 2
                sfields(id_var)%tridfesom = 1 
        endif

! **********************
! ***      ICE       ***
! **********************

!a_ice
        id_var = id%a_ice
        if (id_var > 0) then
                sfields(id_var)%ndims = 1
                sfields(id_var)%variable = 'ice'
                sfields(id_var)%long_name = 'Sea-ice concentration'
                sfields(id_var)%units = '1'
                sfields(id_var)%updated = upd_ice
                sfields(id_var)%bgc = .false.
        endif

    if (mype_world==0) write(*,*) 'a_ice set up'

! **************************
! ***     carbon        ****
! **************************

! DIC
        id_var = id%DIC
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'DIC'
                sfields(id_var)%long_name = 'Dissolved inorganic carbon'
                sfields(id_var)%units = 'mmol(C)* m^{-3}'
                sfields(id_var)%updated = upd_DIC
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 4
                sfields(id_var)%tridfesom = 1002 
        endif

! DOC
        id_var = id%DOC
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'DOC'
                sfields(id_var)%long_name = 'Dissolved organic carbon'
                sfields(id_var)%units = 'mmol(C)* m^{-3}'
                sfields(id_var)%updated = upd_DOC
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 14
                sfields(id_var)%tridfesom = 1012
        endif

! Alkalinity
        id_var = id%Alk
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'Alk'
                sfields(id_var)%long_name = 'Alkalinity'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Alk
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 5
                sfields(id_var)%tridfesom = 1003  
        endif

! pCO2
        id_var = id%pCO2s
        if (id_var > 0) then
                sfields(id_var)%ndims = 1
                sfields(id_var)%variable = 'pCO2s'
                sfields(id_var)%long_name = 'Partial pressure CO2 surface ocean'
                sfields(id_var)%units = 'micro atm'
                sfields(id_var)%updated = upd_pCO2s
                sfields(id_var)%bgc = .true.
        endif

! CO2f
        id_var = id%CO2f
        if (id_var > 0) then
                sfields(id_var)%ndims = 1
                sfields(id_var)%variable = 'CO2f'
                sfields(id_var)%long_name = 'CO2 flux from atmosphere into ocean'
                sfields(id_var)%units = 'mmol(C)* m^{-2}* d^{-1}'
                sfields(id_var)%updated = upd_CO2f
                sfields(id_var)%bgc = .true.
        endif

! Solubility of CO2
        id_var = id%alphaCO2
        if (id_var > 0) then
                sfields(id_var)%ndims = 1
                sfields(id_var)%variable = 'alphaCO2'
                sfields(id_var)%long_name = 'solubility of surface CO2'
                sfields(id_var)%units = 'mol * kg^{-1}* atm^{-1}'
                sfields(id_var)%updated = upd_alphaCO2
                sfields(id_var)%bgc = .true.
        endif

! Piston velocity
        id_var = id%PistonVel
        if (id_var > 0) then
                sfields(id_var)%ndims = 1
                sfields(id_var)%variable = 'Kw660'
                sfields(id_var)%long_name = 'air-sea piston velocity'
                sfields(id_var)%units = 'm/s'
                sfields(id_var)%updated = upd_PistonVel
                sfields(id_var)%bgc = .true.
        endif

! **********************
! ***  Nutrients    ****
! **********************    
    
! DIN
        id_var = id%DIN
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'DIN'
                sfields(id_var)%long_name = 'Dissolved inorganic nitrogen'
                sfields(id_var)%units = 'mmol(N)* m^{-3}'
                sfields(id_var)%updated = upd_DIN
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 3
                sfields(id_var)%tridfesom = 1001  
        endif

! DSi
        id_var = id%DSi
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'DSi'
                sfields(id_var)%long_name = 'Dissolved inorganic Silicate'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_DSi
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 20
                sfields(id_var)%tridfesom = 1018
        endif

! Fe
        id_var = id%Fe
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'Fe'
                sfields(id_var)%long_name = 'Iron'
                sfields(id_var)%units = 'umol* m^{-3}'
                sfields(id_var)%updated = upd_Fe
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 21
                sfields(id_var)%tridfesom = 1019
        endif

! **********************
! *** Chlorophyll   ****
! **********************

! *** Small Phyto 
! PhyN
        id_var = id%PhyN
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'PhyN'
                sfields(id_var)%long_name = 'intracell nitrogen small phytoplankton'
                sfields(id_var)%units = 'mmol(N)* m^{-3}'
                sfields(id_var)%updated = upd_PhyN
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 6
                sfields(id_var)%tridfesom = 1004  
        endif

! PhyC
        id_var = id%PhyC
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'PhyC'
                sfields(id_var)%long_name = 'intracell carbon small phytoplankton'
                sfields(id_var)%units = 'mmol(C)* m^{-3}'
                sfields(id_var)%updated = upd_PhyC
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 7
                sfields(id_var)%tridfesom = 1005   
        endif

! PhyCalc
        id_var = id%PhyCalc
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'PhyCalc'
                sfields(id_var)%long_name = 'calcium carbonate small phytoplankton'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_PhyCalc
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 22
                sfields(id_var)%tridfesom = 1020
        endif

! chlorophyll 
        id_var = id%PhyChl
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'PhyChl'
                sfields(id_var)%long_name = 'Chlorophyll-a small phytoplankton'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_PhyChl
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 8
                sfields(id_var)%tridfesom = 1006 
        endif

! *** diatoms            
! DiaN
        id_var = id%DiaN
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'DiaN'
                sfields(id_var)%long_name = 'intracell nitrogen diatoms'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_DiaN
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 15
                sfields(id_var)%tridfesom = 1013
        endif

! DiaC
        id_var = id%DiaC
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'DiaC'
                sfields(id_var)%long_name = 'intracell carbon diatoms'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_DiaC
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 16
                sfields(id_var)%tridfesom = 1014

        endif

! DiaSi
        id_var = id%DiaSi
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'DiaSi'
                sfields(id_var)%long_name = 'intracell Si diatoms'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_DiaSi
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 18
                sfields(id_var)%tridfesom = 1016
        endif

! chlorophyll
        id_var = id%DiaChl
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'DiaChl'
                sfields(id_var)%long_name = 'Chlorophyll-a diatoms'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_DiaChl
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 17
                sfields(id_var)%tridfesom = 1015
        endif

! *** coccos 
! CoccoN
    !sfields(id%CoccoN)%ndims = 2
    !sfields(id%CoccoN)%variable = 'CoccoN'
    !sfields(id%CoccoN)%long_name = 'intracell nitrogen Coccolithophore'
    !sfields(id%CoccoN)%units = 'mmol m-3'
    !sfields(id%CoccoN)%updated = .false.
    !sfields(id%CoccoN)%bgc = .true.
    !sfields(id%CoccoN  )%tridfesom = 1029 ! CoccoN

! CoccoC
    !sfields(id%CoccoC)%ndims = 2
    !sfields(id%CoccoC)%variable = 'CoccoC'
    !sfields(id%CoccoC)%long_name = 'intracell carbon Coccolithophore'
    !sfields(id%CoccoC)%units = 'mmol C m-3'
    !sfields(id%CoccoC)%updated = .false.
    !sfields(id%CoccoC)%bgc = .true.
    !sfields(id%CoccoC  )%tridfesom = 1030 ! CoccoC

! chlorophyll coccos
    !sfields(id%CoccoChl)%ndims = 2
    !sfields(id%CoccoChl)%nz1 = .true.
    !sfields(id%CoccoChl)%variable = 'CoccoChl'
    !sfields(id%CoccoChl)%long_name = 'Chlorophyll-a Coccolithophore'
    !sfields(id%CoccoChl)%units = 'mg chl m-3'
    !sfields(id%CoccoChl)%updated = .false.
    !sfields(id%CoccoChl)%bgc = .true.
    !sfields(id%CoccoChl)%tridfesom = 1031 ! CoccoChl

! *** Phaeocystis
! PhaeoN
    !sfields(id%PhaeoN)%ndims = 2
    !sfields(id%PhaeoN)%variable = 'PhaeoN'
    !sfields(id%PhaeoN)%long_name = 'intracell nitrogen Phaeocystis'
    !sfields(id%PhaeoN)%units = 'mmol m-3'
    !sfields(id%PhaeoN)%updated = .false.
    !sfields(id%PhaeoN)%bgc = .true.
    !sfields(id%PhaeoN  )%tridfesom = 1032 ! PheaoN

! PhaeoC
    !sfields(id%PhaeoC)%ndims = 2
    !sfields(id%PhaeoC)%variable = 'PhaeoC'
    !sfields(id%PhaeoC)%long_name = 'intracell carbon Phaeocystis'
    !sfields(id%PhaeoC)%units = 'mmol C m-3'
    !sfields(id%PhaeoC)%updated = .false.
    !sfields(id%PhaeoC)%bgc = .true.
    !sfields(id%PhaeoC  )%tridfesom = 1033 ! PheaoC

! chlorophyll Phaeocystis
    !sfields(id%PhaeoChl)%ndims = 2
    !sfields(id%PhaeoChl)%nz1 = .true.
    !sfields(id%PhaeoChl)%variable = 'PhaeoChl'
    !sfields(id%PhaeoChl)%long_name = 'Chlorophyll-a Phaeocystis'
    !sfields(id%PhaeoChl)%units = 'mg chl m-3'
    !sfields(id%PhaeoChl)%updated = .false.
    !sfields(id%PhaeoChl)%bgc = .true.
    !sfields(id%PhaeoChl)%tridfesom = 1034 ! PheaoChl

! *****************************
! ***     Zooplankton      ****
! *****************************

! Zoo1
! Zo1C
        id_var = id%Zo1C
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Zo1C'
                sfields(id_var)%long_name = 'carbon in small zooplankton'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Zo1C
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 12
                sfields(id_var)%tridfesom = 1010
        endif

! Zo1N
        id_var = id%Zo1N
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Zo1N'
                sfields(id_var)%long_name = 'nitrogen in small zooplankton'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Zo1N
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 11
                sfields(id_var)%tridfesom = 1009
        endif

! Zoo2
! Zo2C
        id_var = id%Zo2C
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Zo2C'
                sfields(id_var)%long_name = 'carbon in macrozooplankton'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Zo2C
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 26
                sfields(id_var)%tridfesom = 1024
        endif

! Zo2N
        id_var = id%Zo2N
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Zo2N'
                sfields(id_var)%long_name = 'nitrogen in macrozooplankton'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Zo2N
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 25
                sfields(id_var)%tridfesom = 1023
        endif

! Zoo3
! Zo3C
        id_var = id%Zo3C
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Zo3C'
                sfields(id_var)%long_name = 'carbon in Microzooplankton'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Zo3C
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 30
                sfields(id_var)%tridfesom = 1036
        endif

! Zo3N
        id_var = id%Zo3N
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Zo3N'
                sfields(id_var)%long_name = 'nitrogen in Microzooplankton'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Zo2N
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 29
                sfields(id_var)%tridfesom = 1035
        endif

! *****************************
! ***       Detritus       ****
! *****************************

! small detritus
! DetC
        id_var = id%DetC
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Det2'
                sfields(id_var)%long_name = 'carbon in small detritus'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_DetC
                sfields(id_var)%bgc = .true.
        endif

! DetCalc
        id_var = id%DetCalc
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'DetCalc'
                sfields(id_var)%long_name = 'calcite in small detritus'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_DetCalc
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 23
                sfields(id_var)%tridfesom = 1021 
        endif

! DetN
        id_var = id%DetN
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'DetN'
                sfields(id_var)%long_name = 'nitrogen in small detritus'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_DetN
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 9
                sfields(id_var)%tridfesom = 1007
        endif

! DetSi
        id_var = id%DetSi
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'DetSi'
                sfields(id_var)%long_name = 'silicate in small detritus'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_DetSi
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 19
                sfields(id_var)%tridfesom = 1017 
        endif

! large detritus
! Det2 C
        id_var = id%Det2C
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Det2C'
                sfields(id_var)%long_name = 'carbon in large detritus'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Det2C
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 28
                sfields(id_var)%tridfesom = 1026
        endif

! Det2 Calc
        id_var = id%Det2Calc
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Det2Calc'
                sfields(id_var)%long_name = 'calcite in large detritus'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Det2Calc
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 30
                sfields(id_var)%tridfesom = 1028
        endif

! Det2 N
        id_var = id%Det2N
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Det2N'
                sfields(id_var)%long_name = 'nitrogen in large detritus'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Det2N
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 27
                sfields(id_var)%tridfesom = 1025
        endif

! Det2 Si
        id_var = id%Det2Si
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'Det2Si'
                sfields(id_var)%long_name = 'silicate in large detritus'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_Det2Si
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 29
                sfields(id_var)%tridfesom = 1027
        endif

! *****************************
! ***       other          ****
! *****************************

! DON
        id_var = id%DON
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'DON'
                sfields(id_var)%long_name = 'Dissolved organic nitrogen'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_DON
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 13
                sfields(id_var)%tridfesom = 1011
        endif

! Oxygen
        id_var = id%O2
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%nz1 = .true.
                sfields(id_var)%variable = 'O2'
                sfields(id_var)%long_name = 'Oxygen'
                sfields(id_var)%units = 'mmol* m^{-3}'
                sfields(id_var)%updated = upd_O2
                sfields(id_var)%bgc = .true.
                sfields(id_var)%trnumfesom = 24
                sfields(id_var)%tridfesom = 1022
        endif

! PAR
        id_var = id%PAR
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'PAR'
                sfields(id_var)%long_name = 'photosynthetically active radiation'
                sfields(id_var)%units = 'W *m^{-2}'
                sfields(id_var)%updated = upd_PAR
                sfields(id_var)%bgc = .true.
        endif

! Potential density
        id_var = id%sigma
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'sigma'
                sfields(id_var)%long_name = 'potential density'
                sfields(id_var)%units = 'kg *l^{-1}'
                sfields(id_var)%updated = upd_sigma
                sfields(id_var)%bgc = .false.
        endif

! *****************************
! ***    diagnostics       ****
! *****************************

! NPPn
        id_var = id%NPPn
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'NPPn'
                sfields(id_var)%long_name = 'mean net primary production small phytoplankton'
                sfields(id_var)%units = 'mmol(C)* m^{-2}* d^{-1}'
                sfields(id_var)%updated = upd_NPPn
                sfields(id_var)%bgc = .true.
        endif

! NPPd
        id_var = id%NPPd
        if (id_var > 0) then
                sfields(id_var)%ndims = 2
                sfields(id_var)%variable = 'NPPd'
                sfields(id_var)%long_name = 'mean net primary production diatoms'
                sfields(id_var)%units = 'mmol(C)* m^{-2}* d^{-1}'
                sfields(id_var)%updated = upd_NPPd
                sfields(id_var)%bgc = .true.
        endif

! Export production
        id_var = id%export
        if (id_var > 0) then
                sfields(id_var)%ndims = 1 
                sfields(id_var)%variable = 'export'
                sfields(id_var)%long_name = 'export through particle sinking at 190m'
                sfields(id_var)%units = 'mmol* m^{-2} *day^{-1}'
                sfields(id_var)%updated = upd_export
                sfields(id_var)%bgc = .true.
        endif

    !~ ! TChl
    !~ sfields(id%TChl)%ndims = 2
    !~ sfields(id%TChl)%variable = 'TChl'
    !~ sfields(id%TChl)%long_name = 'Total chlorophyll (PhyChl+DiaChl)'
    !~ sfields(id%TChl)%units = 'mg chl m-3'
    !~ sfields(id%TChl)%updated = .false.
    !~ sfields(id%TChl)%bgc = .true.

    !~ ! TDN
    !~ sfields(id%TDN)%ndims = 2
    !~ sfields(id%TDN)%variable = 'TDN'
    !~ sfields(id%TDN)%long_name = 'Total dissolved nitrogen (DIN+DON)'
    !~ sfields(id%TDN)%units = 'mmol m-3'
    !~ sfields(id%TDN)%updated = .false.
    !~ sfields(id%TDN)%bgc = .true.

    !~ ! TOC
    !~ sfields(id%TOC)%ndims = 2
    !~ sfields(id%TOC)%variable = 'TOC'
    !~ sfields(id%TOC)%long_name = 'Total Organic Carbon (PhyC+DiaC+DetC+DOC+HetC)'
    !~ sfields(id%TOC)%units = 'mmol C m-3'
    !~ sfields(id%TOC)%updated = .false.
    !~ sfields(id%TOC)%bgc = .true.

! *****************************
! *** Reflectance          ****
! *****************************

! Downwelling direct
        do cnt = 1, tlam
                id_var = id%Edz3D(cnt)
                if (id_var > 0) then
                        sfields(id_var)%ndims = 2
                        sfields(id_var)%variable = 'Edz3D' //' '// lams(cnt)
                        sfields(id_var)%long_name = 'Downwelling direct stream of light'
                        sfields(id_var)%units = 'W * m^{-2}'
                        sfields(id_var)%updated = upd_Edz3D
                        sfields(id_var)%bgc = .true.
               endif
        end do
        
! Downwelling diffuse
        do cnt = 1, tlam
                id_var = id%Esz3D(cnt)
                if (id_var > 0) then
                        sfields(id_var)%ndims = 2
                        sfields(id_var)%variable = 'Esz3D' //' '// lams(cnt)
                        sfields(id_var)%long_name = 'Downwelling diffuse stream of light'
                        sfields(id_var)%units = 'W * m^{-2}'
                        sfields(id_var)%updated = upd_Esz3D
                        sfields(id_var)%bgc = .true.
               endif
        end do

! Upwelling        
        do cnt = 1, tlam
                id_var = id%Euz3D(cnt)
                if (id_var > 0) then
                        sfields(id_var)%ndims = 2
                        sfields(id_var)%variable = 'Euz3D' //' '// lams(cnt)
                        sfields(id_var)%long_name = 'Upwelling stream of light'
                        sfields(id_var)%units = 'W * m^{-2}'
                        sfields(id_var)%updated = upd_Euz3D
                        sfields(id_var)%bgc = .true.
               endif
        end do

! Upwelling top of layer        
        do cnt = 1, tlam
                id_var = id%Eutop3D(cnt)
                if (id_var > 0) then
                        sfields(id_var)%ndims = 2
                        sfields(id_var)%variable = 'Eutop3D' //' '// lams(cnt)
                        sfields(id_var)%long_name = 'Upwelling stream of light on the surface of each layer'
                        sfields(id_var)%units = 'W * m^{-2}'
                        sfields(id_var)%updated = upd_Eutop3D
                        sfields(id_var)%bgc = .true.
               endif
        end do
! **************************************
! ***   Set dimensions and offsets   ***
! **************************************

    ! *** Dimensions ***
    do i = 1, nfields
       if (sfields(i)%ndims == 1) then
          sfields(i)%dim = myDim_nod2D
       else if (sfields(i)%ndims == 2) then
          sfields(i)%dim = myDim_nod2D*nlmax
       else
          write (*, '(a,i2,a)') 'FESOM-PDAF: cannot handle', sfields(i)%ndims, ' number of dimensions.'
          write (*,*) sfields(i)%variable
          WRITE(*,*) sfields(i)%long_name
          WRITE(*,*) i
       end if
    end do

! *** Specify offset of fields in pe-local state vector ***
!
!    . . . A . . . . . B . . . . . C . . . . . D
!         . .         / .         / .         . .
!        .   .   2   /   .   3   /   .   5   .   .
! .     .     .     /     .     /     .     .     .
!  .   .   1   .   /   3   .   /   4   .   .       
!   . .         . /         . /         . .        
!    A . . . . . B . . . . . C . . . . . D . . . . 
!
!  A:    Internal nodes of left PE
!  B:    Internal nodes of left PE, simultanesously external nodes of right PE
!  C:    External nodes of left PE, simultanesously internal nodes of right PE
!  D:    Internal nodes of right PE
!  1:    Internal element of left PE, simultanesously wide-halo element of right PE (shares node B with right PE)
!  2:    Internal element of left PE, simultanesously small-halo element of right PE (shares edge BB with right PE)
!  3:    Internal (shared) elements of both PEs
!  4:    Small-halo element of left PE (shares edge CC with left PE), simultanesously internal element of right PE
!  5:    Wide-halo element of left PE (shares node C with left PE), simultanesously internal element of right PE
!
!  myDim_nod2D:      Number of internal nodes (A+B)
!  eDim_nod2D:       Number of external nodes (C)
!
!  myDim_elem2D:     Number of internal elements (1+2+3)
!  eDim_elem2D:      Number of small-halo elements (4)
!  xDim_elem2D:      Number of wide-halo elements (5)
!
!  mesh_fesom%nl:    Maximum number of fesom levels (1 is air-sea interface)
!  mesh_fesom%nl-1:  Maximum number of fesom layers (1 is surface layer, 0-5m)

!  CORE2 mesh: deepest wet cells at mesh_fesom%nl-2 or nlmax=46

    ! *** Define offsets in state vector ***
    sfields(1)%off = 0
    do i = 2, nfields
       sfields(i)%off = sfields(i-1)%off + sfields(i-1)%dim
    end do

! **************************************
! ***      updated during DA         ***
! **************************************
    
! The logical "updated" describes whether a variables is updated in at least one sweep

  ! In case of diagnostic variables, "updated" is False. Setting diagnostics variables to False
  ! and the others to True, is set in namelist.
  ! In case of weak coupling and only physics or BGC assimilation, "updated" is False for the
  ! other type of fields. This is done below.
  ! "updated" is used in init_dim_l_pdaf: only updated fields are included in local state
  ! "updated" is used in the output routine: option to write out only updated fields


! *** General settings ***
!from previus set up without groups

    !! Physics not assimilated and coupling weak: No update to physics
    !if ((.not. assimilatePHY) .and. (cda_bio=='weak')) then
    !   sfields(phymin: phymax)%updated = .false.
    !endif

    !! BGC not assimilated and coupling weak: No update to BGC
    !if ((.not. assimilateBGC) .and. (cda_phy=='weak')) then
    !   sfields(bgcmin: bgcmax)%updated = .false.
    !endif


  end subroutine init_sfields
! ===================================================================================

!> Count 2D and 3D fields and initialize index arrays
!!
!! This functioanlity of optional and not used elsewhere
!! in the code.
!!
  subroutine set_field_types(verbose)

    use parallel_pdaf_mod, &
         only: mype_world

    implicit none

! *** Arguments ***
    integer, intent(in) :: verbose     ! Verbosity level

! *** Local variables ***
    integer :: i, cnt      ! Counters
    integer, allocatable :: ids_3D(:)         ! List of 3D-field IDs
    integer, allocatable :: ids_2D(:)         ! """       2D fields """
    integer, allocatable :: ids_phy(:)        ! """       physics fields """
    integer, allocatable :: ids_bgc(:)        ! """       biogeochem. fields """
    integer, allocatable :: ids_tr3D(:)       ! List of 3D model tracer IDs
    integer :: nfields_3D                     ! Number of 3D fields in state vector
    integer :: nfields_2D                     ! """       2D fields """
    integer :: nfields_phy                    ! """       physics fields """
    integer :: nfields_bgc                    ! """       biogeochem. fields """
    integer :: nfields_tr3D                   ! Number of 3D model tracer fields in state vector


! **************************************
! ***  Indices of by type of field   ***
! **************************************
  
! ***2D/3D***

    ! count number of 3D and 2D fields (tracers and diagnostics)
    nfields_3D = 0
    nfields_2D = 0
    do i=1,nfields
       if (sfields(i)%ndims == 2) nfields_3D = nfields_3D + 1
       if (sfields(i)%ndims == 1) nfields_2D = nfields_2D + 1
    enddo
    
    ! Store indices of 3D and 2D fields
    allocate(ids_3D(nfields_3D))
    cnt = 1
    do i=1,nfields
       if (sfields(i)%ndims == 2) then
          ids_3D(cnt) = i
          cnt = cnt+1
       endif
    enddo
    if (mype_world==0 .and. verbose>2) then
       do cnt=1,nfields_3D
          write (*,'(a, 10x,3a,1x,7a)') &
               'FESOM-PDAF', '3D fields in state vector: ', sfields(ids_3D(cnt))%variable
       enddo
    endif
    
    allocate(ids_2D(nfields_2D))
    cnt = 1
    do i=1,nfields
       if (sfields(i)%ndims == 1) then
          ids_2D(cnt) = i
          cnt = cnt+1
       endif
    enddo
    if (mype_world==0 .and. verbose>2) then
       do cnt=1,nfields_2D
          write (*,'(a, 10x,3a,1x,7a)') &
               'FESOM-PDAF', '2D fields in state vector: ', sfields(ids_2D(cnt))%variable
       enddo
    endif
    
! *** model 3D tracers ***

    ! count number of 3D model tracer fields in state vector (only tracers)
    nfields_tr3D = 0
    do i=1,nfields
       if (sfields(i)%trnumfesom > 0) nfields_tr3D = nfields_tr3D + 1
    enddo
    
    ! Store indices of 3D model tracer fields in state vector
    allocate(ids_tr3D(nfields_tr3D))
    cnt = 1
    do i=1,nfields
       if (sfields(i)%trnumfesom > 0) then
          ids_tr3D(cnt)=i
          sfields(i)%id_tr = cnt
          cnt = cnt+1
       endif
    enddo
    if (mype_world==0 .and. verbose>2) then
       do cnt=1,nfields_tr3D
          write (*,'(a, 10x,3a,1x,7a)') &
               'FESOM-PDAF', '3D model tracer fields in state vector: ', sfields(ids_tr3D(cnt))%variable
       enddo
    endif
  
! *** phy/bgc fields  ***_

    ! count number of phy/bgc fields
    nfields_phy = 0
    nfields_bgc = 0
    do i=1,nfields
       if (      sfields(i)%bgc) nfields_bgc = nfields_bgc + 1
       if (.not. sfields(i)%bgc) nfields_phy = nfields_phy + 1
    enddo
    
    ! Store indices of BGC and physics fields in state vector
    allocate(ids_bgc(nfields_bgc))
    cnt = 1
    do i=1,nfields
       if (sfields(i)%bgc) then
          ids_bgc(cnt) = i
          cnt = cnt+1
       endif
    enddo
    if (mype_world==0 .and. verbose>2) then
       do cnt=1,nfields_bgc
          write (*,'(a, 10x,3a,1x,7a)') &
               'FESOM-PDAF', 'bgc fields in state vector: ', sfields(ids_bgc(cnt))%variable
       enddo
    endif
    
    allocate(ids_phy(nfields_phy))
    cnt = 1
    do i=1,nfields
       if (.not. sfields(i)%bgc) then
          ids_phy(cnt) = i
          cnt = cnt+1
       endif
    enddo
    if (mype_world==0 .and. verbose>2) then
       do cnt=1,nfields_phy
          write (*,'(a, 10x,3a,1x,7a)') &
               'FESOM-PDAF', 'phy fields in state vector: ', sfields(ids_phy(cnt))%variable
       enddo
    endif


  end subroutine set_field_types

  end module statevector_pdaf
