MODULE mod_atmos_ens_stochasticity

! Description: Adds stochastic synoptic variability to the atmospheric
! forcing fields for an ensemble of atmospheric forcings.

! Routines for this modules are:
! --- init_atmos_ens_stochasticity()
! --- add_atmos_ens_stochasticity(istep)
! --- init_atmos_stochasticity_output()
! --- write_atmos_stochasticity_output(istep)
! --- write_atmos_stoch_restart
! --- read_atmos_stochasticity_restart()
! --- compute_ipsr()
!these are in the file atmos_ens_stochasticity.F90

! Covariance file contains statistical information on all 9 atmospheric
! forcing fields. Which of these fields shall be perturbed is set
! in namelist.fesom.pdaf.

  USE mpi
  USE parallel_pdaf_mod, &
       ONLY: mype_filter, mype_model, mype_world, MPIerr, &
             COMM_filter, filterpe, task_id, COMM_model
  USE assim_pdaf_mod, &
       ONLY: dim_ens, dim_state_p, forget
  USE fesom_pdaf, &
       ONLY: partit, cyearnew, cyearold, step_per_day, exchange_nod, &
       myDim_nod2D, eDim_nod2D, &
       atmdata, i_xwind, i_ywind, i_humi, &
       i_qsr, i_qlw, i_tair, i_prec, i_mslp, i_snow
       
  IMPLICIT NONE
  
 ! INCLUDE 'netcdf.inc'

  INTEGER                 :: cnt                     ! counters
  REAL(8),ALLOCATABLE, save  :: eof_p_cvrfile(:,:)      ! Matrix of eigenvectors of covariance matrix (all fields read from file)
  REAL(8),ALLOCATABLE, save  :: eof_p(:,:)              ! Matrix of eigenvectors of covariance matrix (only activated fields)
  REAL(8),ALLOCATABLE, save  :: svals(:)                ! Singular values
  REAL(8),ALLOCATABLE        :: omega(:,:)           ! Transformation matrix Omega
  REAL(8),ALLOCATABLE        :: omega_v(:)           ! Transformation vector for local ensemble member
  INTEGER                 :: rank                    ! Rank stored in cov.-file
  CHARACTER(len=4)        :: mype_string             ! String for process rank
  CHARACTER(len=110)      :: filename                ! Name of covariance netCDF file
  INTEGER,          save  :: nfields_cvrfile         ! Number of atmospheric forcing fields in covariance netCDF file
  INTEGER,          save  :: nfields                 ! Number of activated atmospheric forcing fields
  REAL(8),ALLOCATABLE        :: perturbation(:)         ! Vector containing perturbation field for local ensemble member
  character(len=150) :: path_atm_cov
  
  REAL,ALLOCATABLE, save  :: perturbation_humi (:)   ! Final perturbations for each variable
  REAL,ALLOCATABLE, save  :: perturbation_prec (:)
  REAL,ALLOCATABLE, save  :: perturbation_snow (:)
  REAL,ALLOCATABLE, save  :: perturbation_mslp (:)
  REAL,ALLOCATABLE, save  :: perturbation_qlw  (:)
  REAL,ALLOCATABLE, save  :: perturbation_qsr  (:)
  REAL,ALLOCATABLE, save  :: perturbation_tair (:)
  REAL,ALLOCATABLE, save  :: perturbation_xwind(:)
  REAL,ALLOCATABLE, save  :: perturbation_ywind(:)
  
  REAL(8),ALLOCATABLE, save  :: ipsr(:)                 ! instantaneous potential solar radiation
  
  REAL,ALLOCATABLE, save  :: atmdata_debug(:,:)
  
  TYPE field_ids
     INTEGER :: humi
     INTEGER :: prec
     INTEGER :: snow
     INTEGER :: mslp
     INTEGER :: qlw 
     INTEGER :: qsr
     INTEGER :: tair
     INTEGER :: xwind
     INTEGER :: ywind
  END TYPE field_ids
  
  ! Type variable holding field IDs in atmospheric state vector
  TYPE(field_ids)    , save :: id_atm               ! field IDs of perturbed fields
  TYPE(field_ids)    , save :: id_cvrf              ! field IDs of all fields in covariance file
  
  INTEGER,ALLOCATABLE, save :: atm_offset(:)        ! offset of perturbed fields in atmospheric state vector
  INTEGER,ALLOCATABLE, save :: atm_offset_cvrf(:)   ! offset of fields hold in covariance fields
  
  CHARACTER(len=200) :: fname_atm      ! filename to write atmospheric stochasticity at time step
  CHARACTER(len=200) :: fname_restart  ! filename to write restart information

    ! which atmospheric fields to be perturbed (set in namelist)
    logical :: disturb_xwind=.false.
    logical :: disturb_ywind=.false.
    logical :: disturb_humi=.false.
    logical :: disturb_qlw=.false.
    logical :: disturb_qsr=.false.
    logical :: disturb_tair=.false.
    logical :: disturb_prec=.false.
    logical :: disturb_snow=.false.
    logical :: disturb_mslp=.false.

LOGICAL :: atmos_stochasticity_ON   ! if any atmospheric fields to be perturbed

REAL :: varscale_wind = 0.2         ! scaling factors
REAL :: varscale_humi = 1.0
REAL :: varscale_qlw  = 1.0
REAL :: varscale_qsr  = 1.0
REAL :: varscale_tair = 1.0
REAL :: varscale_prec = 1.0
REAL :: varscale_snow = 1.0
REAL :: varscale_mslp = 0.2

LOGICAL :: write_atmos_st = .false. ! wether to protocol the perturbed atmospheric fields,
                                    ! i.e. writing at every time step

REAL :: stable_rmse = 0 ! (ocean temperature) ensemble spread after 16 months of assimilation 

INTERFACE
   SUBROUTINE DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
     CHARACTER(len=1), INTENT(IN) :: TRANS
     INTEGER,          INTENT(IN) :: M, N, LDA, INCX, INCY
     REAL(8),          INTENT(IN) :: ALPHA, BETA
     REAL(8),          INTENT(IN) :: A(LDA,*), X(*)
     REAL(8),          INTENT(INOUT) :: Y(*)
   END SUBROUTINE DGEMV
END INTERFACE
END MODULE mod_atmos_ens_stochasticity
