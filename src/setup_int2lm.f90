!+ External procedure for the setup of the interpolation program
!==============================================================================

SUBROUTINE setup_int2lm (ierror, yerror)

!==============================================================================
!
! Description:
!   This external subroutine organizes the setup of the unified interpolation
!   program.
!
! Method:
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2005/04/11 Ulrich Schaettler
!  Initial release for INT2LM
! V1_5         2007/07/09 Guy de Morsier
!  Added some comments
! V1_6         2007/09/07 Burkhardt Rockel, Ulrich Schaettler
!  Added option lcm2lm to interpolate data from a climate model
!  Added calls to SR decompose_cm and coor_cm_lm
! V1_8         2008/05/29 Ulrich Schaettler
!  Eliminated ldebug and replaced it with idbg_level
! V1_9         2009/09/03 Guy DeMorsier
!  Adapted definitions of soil properties
! V1_10        2009/12/17 Ulrich Schaettler
!  Added treatment of Unified Model data (decompose, grib-table)
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler
!  Unified variable table for all input models and dwdlib/grib_api
!  Added treatment of JMA and NCEP data (allocation, decomposition)
! V1_19        2012/06/06 Ulrich Schaettler, Burkhardt Rockel
!  Call setup_vartab_grid_in also for HIRLAM input
!  Adaptations to environment.f90 (here: init_environment, init_procgrid),
!    because of unification with COSMO-Model 4.23
!  Call new subroutine setup_clm (CLM)
! V1_22        2013/07/11 Ulrich Schaettler
!  Enlarged interface to init_environment to initialize MPI type for special
!   grib_api integer
!  Setup some I/O organizational variables
!  Allocate sendbuf in any case
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters , ONLY :   &
    ireals,    & ! KIND-type parameters for real variables
    iintegers, & ! kind-type parameter for "normal" integer variables
    irealgrib, & ! KIND parameter for the Reals in the GRIB library
    int_ga       ! integer precision for grib_api: length of message in bytes

!------------------------------------------------------------------------------

USE data_int2lm_constants,   ONLY :   &
    T0,      & ! 0 degree Celsius                              [Kelvin]
    R_d,     & ! gas constant for dry air                      [J/K*kg]
    R_v,     & ! gas constant for water vapor                  [J/K*kg]
    Rdv,     & ! = R_d/R_v,
    Rvd_m_o, & ! = R_v/R_d - 1.0,
    O_m_rdv, & ! = 1. - Rdv
    Cp_d,    & ! specific heat of dry air at constant pressure [J/K*kg]
    lh_v,    & ! Latent heat of vaporization                   [J/kg]
    lh_f,    & ! Latent heat of fusion                         [J/kg]
    lh_s,    & ! Latent heat of sublimation                    [J/kg]
    G,       & ! gravity at sea level                          [ms-2]
    r_earth, & ! mean radius of the earth                      [m]
    Day_len, & ! sidereal day (Sterntag)                       [s]
    B1,      & !  a
    B2_w,    & !  b
    B2_i,    & !
    B3,      & !  c/b (0 degree Celsius [Kelvin])
    B4_i,    & !
    B4_w,    & !  d
    Pi,      & !
    raddeg,  & !
    degrad,  & !
    Pid5,    & !
    Omcor,   & !
    wimax,   & ! used in himbla to minimize wi
    porb,       & ! pore volume of COSMO soil types
    porb_ec_1s, & ! pore volume of IFS 1 soil type
    adpb,       & ! air dryness point of COSMO soil types
    pwpb,       & ! permanent wilting point of COSMO soil types
    pwpb_ec,    & ! permanent wilting point of IFS soil types with CY32r3
    pwpb_ec_1s, & ! permanent wilting point of IFS 1 soil type
    fcb,        & ! field capacity of COSMO soil types
    fcb_ec,     & ! field capacity of IFS soil types after CY32r3
    fcb_ec_1s     ! field capacity of IFS 1 soil type

!------------------------------------------------------------------------------

USE data_int2lm_control,     ONLY :   &
    ndebug,       & ! unit number for file with debug information
    noutput,      & ! unit number for output file
    lgme2lm,      & ! if .TRUE., gme->lm
    lgfs2lm,      & ! if .TRUE., gfs->lm
    lgsm2lm,      & ! if .TRUE., gsm->lm
    lec2lm,       & ! if .TRUE., ec ->lm
    llm2lm,       & ! if .TRUE., lm ->lm
    lum2lm,       & ! if .TRUE., um ->lm
    lhir2lm,      & ! if .TRUE., hirlam ->lm
    lcm2lm,       & ! if .TRUE., cm ->lm   !_br
    idbg_level,   & ! to control verbosity of output
    lprintdeb_all,& ! whether all PEs print debug output
    lclock          ! if .TRUE., system clock is present

!------------------------------------------------------------------------------

USE data_int2lm_io,          ONLY :   &
    ylevltypes_in,  & ! to convert GRIB1 level types to grib_api string typeOfLevel
                      ! for incoming data
    ylevltypes_out, & ! to convert GRIB1 level types to grib_api string typeOfLevel
                      ! for COSMO-Model data that is written
    ysteptypes,     & ! to convert GRIB1 time range indicator to grib_api stinf stepType
    rscalefac         ! Array to convert GRIB2 scale factors to real numbers

!------------------------------------------------------------------------------

USE data_grid_lm,     ONLY :   &
    ie2lm_max,   & ! Max. of ie2lm on all processors
    je2lm_max,   & ! Max. of je2lm on all processors
    kelm_tot,    & ! number of vertical levels for LM
    kelm,        & ! ke for LM
    ie2lm_tot,   & ! = ielm_tot + 2
    je2lm_tot,   & ! = jelm_tot + 2
    ie2lm,       & !
    je2lm,       & !
    istartpar,   & ! start index for computations in the parallel program
    iendpar,     & ! end index for computations in the parallel program
    jstartpar,   & ! start index for computations in the parallel program
    jendpar        ! end index for computations in the parallel program

!------------------------------------------------------------------------------

USE data_grid_in,     ONLY :   &
    ids,         & ! start index of diamonds (ids = 1)
    ide,         & ! end index of diamonds (ide = 10)
    ni_gme,      & ! resolution of GME
    ilim1,       & ! decomposition limits in x-direction (formerly 1)
    ilim2,       & ! decomposition limits in y-direction (formerly 2)
    igg1s,       & ! start index of global array-dimension in x-direction
    igg1e,       & ! end index of global array-dimension in x-direction
    igg2s,       & ! start index of global array-dimension in y-direction
    igg2e,       & ! end index of global array-dimension in y-direction
    ig1s,        & ! start index of array-dimension in x-direction
    ig1e,        & ! end index of array-dimension in x-direction
    ig2s,        & ! start index of array-dimension in y-direction
    ig2e           ! end index of array-dimension in y-direction

!------------------------------------------------------------------------------

USE data_fields_lm,   ONLY :   &
    latlm_m,     & ! latitudes of the LM grid points
    lonlm_m,     & ! longitudes of the LM grid points
    latlm_u,     & ! latitudes of the LM u grid points
    lonlm_u,     & ! longitudes of the LM u grid points
    latlm_v,     & ! latitudes of the LM v grid points
    lonlm_v,     & ! longitudes of the LM v grid points
    i_index,     & ! i-index of coarse mesh grid point which is lower left
                   ! to a given LM (fine mesh) grid point
    j_index,     & ! j-index of coarse mesh grid point which is lower left
                   ! to a given LM (fine mesh) grid point
    x_wght,      & ! relative distance between x- (i-) coarse mesh and
                   ! fine mesh grid points
    y_wght         ! relative distance between y- (j-) coarse mesh and
                   ! fine mesh grid points

!------------------------------------------------------------------------------

USE data_fields_in,    ONLY :   &
    lat_coarse_m,& ! latitudes of the LM grid points
    lon_coarse_m,& ! longitudes of the LM grid points
    lat_coarse_u,& ! latitudes of the LM u grid points
    lon_coarse_u,& ! longitudes of the LM u grid points
    lat_coarse_v,& ! latitudes of the LM v grid points
    lon_coarse_v   ! longitudes of the LM v grid points

!------------------------------------------------------------------------------

USE data_int2lm_parallel,    ONLY :   &
    lasync_io,       & ! if .TRUE.: the model runs with extra PEs for
                       ! asynchronous IO
    nprocx,          & ! number of processors in x-direction
    nprocy,          & ! number of processors in y-direction
    nprocio,         & ! number of extra processors for asynchronous IO
    num_compute,     & ! number of compute PEs
    num_io,          & ! number of IO PEs
    nproc1,          & ! number of processors in first direction for GME
                       ! = nprocy (!)
    nproc2,          & ! number of processors in second direction for GME
                       ! = nprocx (!)
    nproc,           & ! total number of processors: nprocx * nprocy + nprocio

    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains

    my_world_id,     & ! rank of this subdomain in the global communicator
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_pos,     & ! position of this subdomain in the cartesian grid
                       ! in x- and y-direction
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    isubpos,         & ! positions of the subdomains in the total domain.
    isubpos_coarse,  & ! positions of the subdomains in the total domain.
    icomm_compute,   & ! communicator for the group of compute PEs
    igroup_cart,     & ! group of the compute PEs
    icomm_cart,      & ! communicator for the virtual cartesian topology
    icomm_row,       & ! communicator for a east-west row of processors
    igroup_world,    & ! group that belongs to MPI_COMM_WORLD, i.e. all
                       ! processors
    icomm_world,     & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_grib,        & ! determines the REAL type for the GRIB library
    imp_integers,    & ! determines the correct INTEGER type used in the
                       ! model for MPI
    imp_integ_ga,    & ! determines the correct INTEGER type used for grib_api
    imp_byte,        & ! determines the correct BYTE type used in the model
                       ! for MPI
    imp_character,   & ! determines the correct CHARACTER type used in the
                       ! model for MPI
    imp_logical,     & ! determines the correct LOGICAL   type used in the
                       ! model for MPI
    lcompute_pe,     & ! indicates whether this is a compute PE or not
    lreorder,        & ! during the creation of the virtual topology the
                       ! ranking of the processors may be reordered
    sendbuf,         & ! sending buffer for boundary exchange
    isendbuflen        ! length of one column of sendbuf.

!------------------------------------------------------------------------------

USE src_namelists,      ONLY:  read_namelists
USE src_decomposition,  ONLY:  decompose_lm, decompose_gme, decompose_cm, &
                               decompose_coarse_grid, check_decomposition
USE src_grids,          ONLY:  compute_geo_coord, coor_coarse_grid_lm,    &
                               coor_gme_lm, coor_cm_lm
USE src_memory,         ONLY:  alloc_lm, alloc_gme, alloc_coarse_grid
USE src_gribtabs,       ONLY:  setup_vartab_lm, setup_vartab_grid_in
USE utilities,          ONLY:  elapsed_time
USE environment,        ONLY:  init_environment, init_procgrid,           &
                               get_free_unit, release_unit
USE parallel_utilities, ONLY:  init_par_utilities
USE gme_utilities,      ONLY:  setup_xd
USE mpe_io,             ONLY:  mpe_io_init, mpe_io_reconfig

#ifdef NETCDF
USE clm_utilities,      ONLY:  setup_clm
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Parameterlist
CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

! Local variables
INTEGER (KIND=iintegers)              ::  &
  nzstat, izerror, izbuflen, izdebug,     &
  iztest_compute     ! test communicator for MPE_IO

REAL    (KIND=ireals)                 ::  &
  zdreal

CHARACTER (LEN=30) yzroutine
CHARACTER (LEN=80) yzerror

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!- Begin Subroutine setup_int2lm
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierror    = 0_iintegers
  yerror    = '         '
  izerror   = 0_iintegers
  yzerror   = '         '
  nzstat    = 0_iintegers

  yzroutine = 'setup_int2lm'

  !----------------------------------------------------------------------------
  ! Section 1.1: Initialization of the desired environment
  !----------------------------------------------------------------------------

  CALL init_environment (nproc, my_world_id, icomm_world, igroup_world,       &
                         imp_integers, imp_reals, imp_grib, imp_byte,         &
                         imp_character, imp_logical, imp_integ_ga,            &
                         iexch_req, irealgrib, yzerror, izerror )

  IF (izerror /= 0) THEN
    yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
    ierror = 1_iintegers
    RETURN
  ENDIF

  IF (my_world_id == 0) THEN
    PRINT *,'  SETUP OF INT2LM'
    PRINT *,'    INITIALIZATIONS '
    PRINT *,'       Info about KIND-parameters:   iintegers / MPI_INT = ', iintegers, imp_integers
    PRINT *,'                                     int_ga    / MPI_INT = ', int_ga,    imp_integ_ga
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 1.2: Initialization of the timing
  !----------------------------------------------------------------------------

  CALL elapsed_time (zdreal, izerror)
  IF (izerror == 0) THEN
    lclock = .TRUE.
  ELSE
    IF (my_world_id == 0) THEN
      PRINT *, ' WARNING:     !!! NO SYSTEM CLOCK PRESENT !!! '
    ENDIF
    lclock = .FALSE.
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 1.3: Initializations of physical constants
  !----------------------------------------------------------------------------

  ! mathematical constants
  ! ----------------------
  pi       = 4.0_ireals * ATAN (1.0_ireals)
  
  ! physical constants and related variables
  ! ----------------------------------------
  T0      =   273.15_ireals
  R_d     =     2.8705E2_ireals
  R_v     =     4.6151E2_ireals
  Rdv     = R_d/R_v
  Rvd_m_o = R_v/R_d - 1.0_ireals
  O_m_rdv = 1.0_ireals - Rdv
  
  Cp_d    =     1.005E3_ireals
  lh_v    =     2.501E6_ireals
  lh_f    =     0.334E6_ireals
  lh_s    = lh_v + lh_f
  G       =     9.80665_ireals
  r_earth =     6.371229E6_ireals
  Day_len = 86164.09054_ireals
  B1      =   610.78_ireals
  B2_w    =    17.2693882_ireals
  B2_i    =    21.8745584_ireals
  B3      =   273.16_ireals
  B4_i    =     7.66_ireals
  B4_w    =    35.86_ireals
  raddeg  =    57.29577951_ireals
  degrad  =     0.0174532925_ireals
  Pid5    =  Pi/5._ireals
  Omcor   = 2.0_ireals*Pi/Day_len
  
  ! field capacity FC of soil types
  fcb = (/0.0_ireals,   0.0_ireals,   0.196_ireals, 0.260_ireals,   &
          0.340_ireals, 0.370_ireals, 0.463_ireals, 0.763_ireals, 1.0_ireals/)

  ! poren volume POR of soil types
  porb= (/0.0_ireals,   0.0_ireals,   0.364_ireals, 0.445_ireals,   &
          0.455_ireals, 0.475_ireals, 0.507_ireals, 0.863_ireals, 1.0_ireals/)

  ! permanent wilting point PWP of soil types
  pwpb = (/0.0_ireals,  0.0_ireals,   0.042_ireals, 0.100_ireals,   &
          0.110_ireals, 0.185_ireals, 0.257_ireals, 0.265_ireals, 1.0_ireals/)

  ! air dryness point of gr. types
  adpb= (/0.0_ireals,   0.0_ireals,   0.012_ireals, 0.030_ireals,   &
          0.035_ireals, 0.060_ireals, 0.065_ireals, 0.098_ireals, 0.0_ireals/)

  ! FC, POR and PWP of IFS 1 soil type
  ! From IFS documentation Cy28r1, Physical Processes,
  ! Land Surface parametrization: between Equs. 7.64 and 7.65
  fcb_ec_1s  = 0.323_ireals
  porb_ec_1s = 0.472_ireals
  pwpb_ec_1s = 0.171_ireals

  ! From IFS documentation Cy32r3, Physical Processes,
  ! http://www.ecmwf.int/products/changes/soil_hydrology_cy32r3/new_soil_hydrology.pdf
  fcb_ec  = (/ 0.242_ireals, 0.346_ireals, 0.382_ireals, 0.448_ireals, 0.541_ireals, 0.662_ireals/)
  pwpb_ec = (/ 0.059_ireals, 0.151_ireals, 0.133_ireals, 0.279_ireals, 0.335_ireals, 0.267_ireals/)

!------------------------------------------------------------------------------
! Section 2: Read all NAMELIST-variables and initialize I/O org variables
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 2.1: Read all NAMELIST-variables
  !----------------------------------------------------------------------------

  CALL read_namelists (izerror, yzerror)
  IF (izerror /= 0) THEN
    yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
    ierror = izerror
    RETURN
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 2.2: Initialize I/O organizational variables
  !----------------------------------------------------------------------------

  ! initialize ylevltypes_in to convert GRIB1 level types to grib_api string typeOfLevel
  ! for input data: GRIB1
  ylevltypes_in(  :,:) = 'dummy                         '
  ylevltypes_in(  1,1) = 'surface'
  ylevltypes_in(  2,1) = 'cloudBase'
  ylevltypes_in(  3,1) = 'cloudTop'
  ylevltypes_in(  4,1) = 'isothermZero'
  ylevltypes_in(  8,1) = 'nominalTop'
  ylevltypes_in(100,1) = 'isobaricInhPa'
  ylevltypes_in(102,1) = 'meanSea'
  ylevltypes_in(105,1) = 'heightAboveGround'
  ylevltypes_in(109,1) = 'hybrid'
  ylevltypes_in(110,1) = 'hybridLayer'
  ylevltypes_in(111,1) = 'depthBelowLand'
  ylevltypes_in(112,1) = 'depthBelowLandLayer'
  ylevltypes_in(200,1) = 'entireAtmosphere'
  ylevltypes_in(211,1) = 'snowLayer'

  ! initialize ylevltypes_in to convert GRIB1 level types to grib_api string typeOfLevel
  ! for input data: GRIB2
  ylevltypes_in(  1,2) = 'surface'
  ylevltypes_in(  2,2) = 'cloudBase'
  ylevltypes_in(  3,2) = 'cloudTop'
  ylevltypes_in(  4,2) = 'isothermZero'
  ylevltypes_in(  8,2) = 'nominalTop'
  ylevltypes_in(100,2) = 'isobaricInhPa'
  ylevltypes_in(102,2) = 'meanSea'
  ylevltypes_in(105,2) = 'heightAboveGround'
  IF (llm2lm) THEN
    ylevltypes_in(109,2) = 'generalVertical'
    ylevltypes_in(110,2) = 'generalVerticalLayer'
  ELSE
    ylevltypes_in(109,2) = 'hybrid'
    ylevltypes_in(110,2) = 'hybridLayer'
  ENDIF
  ylevltypes_in(111,2) = 'depthBelowLand'
  ylevltypes_in(112,2) = 'depthBelowLandLayer'
  ylevltypes_in(200,2) = 'entireAtmosphere'
  ylevltypes_in(211,2) = 'snowLayer'

  ! initialize ylevltypes_out to convert GRIB1 level types to grib_api string typeOfLevel
  ! for outgoing COSMO data: GRIB1
  ylevltypes_out(  :,:) = 'dummy                         '
  ylevltypes_out(  1,1) = 'surface'
  ylevltypes_out(  2,1) = 'cloudBase'
  ylevltypes_out(  3,1) = 'cloudTop'
  ylevltypes_out(  4,1) = 'isothermZero'
  ylevltypes_out(  8,1) = 'nominalTop'
  ylevltypes_out(100,1) = 'isobaricInhPa'
  ylevltypes_out(102,1) = 'meanSea'
  ylevltypes_out(105,1) = 'heightAboveGround'
  ylevltypes_out(109,1) = 'hybrid'
  ylevltypes_out(110,1) = 'hybridLayer'
  ylevltypes_out(111,1) = 'depthBelowLand'
  ylevltypes_out(112,1) = 'depthBelowLandLayer'
  ylevltypes_out(200,1) = 'entireAtmosphere'
  ylevltypes_out(211,1) = 'snowLayer'

  ! initialize ylevltypes_out to convert GRIB1 level types to grib_api string typeOfLevel
  ! for outgoing COSMO data: GRIB2
  ylevltypes_out(  1,2) = 'surface'
  ylevltypes_out(  2,2) = 'cloudBase'
  ylevltypes_out(  3,2) = 'cloudTop'
  ylevltypes_out(  4,2) = 'isothermZero'
  ylevltypes_out(  8,2) = 'nominalTop'
  ylevltypes_out(100,2) = 'isobaricInhPa'
  ylevltypes_out(102,2) = 'meanSea'
  ylevltypes_out(105,2) = 'heightAboveGround'
  ylevltypes_out(109,2) = 'generalVertical'
  ylevltypes_out(110,2) = 'generalVerticalLayer'
  ylevltypes_out(111,2) = 'depthBelowLand'
  ylevltypes_out(112,2) = 'depthBelowLandLayer'
  ylevltypes_out(200,2) = 'entireAtmosphere'
  ylevltypes_out(211,2) = 'snowLayer'

  ! initialize steptype to convert GRIB1 time range indicator to grib_api stinf stepType
  ysteptypes( : ) = 'dummy'
  ysteptypes(  0) = 'instant'
  ysteptypes(  2) = 'diff'
  ysteptypes(  3) = 'avg'
  ysteptypes(  4) = 'accum'
  ysteptypes( 10) = 'instant'

  ! initialize scalefactors
  rscalefac(:)     = -1.0_ireals
  rscalefac(0)     = 1.0_ireals
  rscalefac(1)     = 1.0E-1_ireals
  rscalefac(2)     = 1.0E-2_ireals
  rscalefac(3)     = 1.0E-3_ireals
  rscalefac(4)     = 1.0E-4_ireals
  rscalefac(5)     = 1.0E-5_ireals
  rscalefac(6)     = 1.0E-6_ireals
  rscalefac(7)     = 1.0E-7_ireals
  rscalefac(8)     = 1.0E-8_ireals
  rscalefac(9)     = 1.0E-9_ireals

#ifdef NETCDF   
  ! This works presently only with netCDF I/O
  IF (lcm2lm) CALL setup_clm
#endif

!------------------------------------------------------------------------------
! Section 3: Domain Decomposition and Grid Intersections
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 3.1: Allocation of memory
  !----------------------------------------------------------------------------

  ! Allocate space for organization variables and the geographical coordinates
  ALLOCATE ( isubpos(0:num_compute-1,4)   , STAT=nzstat )
  isubpos(:,:) = 0

  IF (lgme2lm) THEN
    ! Allocate and Set decomposition limits
    ALLOCATE (ilim1(0:nprocy), ilim2(0:nprocx), STAT=nzstat)
  ELSE  ! all other input models
    ALLOCATE ( isubpos_coarse(0:num_compute-1,4)   , STAT=nzstat )
    isubpos_coarse(:,:) = 0
  ENDIF

  IF (nzstat /= 0) THEN
    ierror    = 1101
    yerror    = ' ERROR    *** Allocation of space failed ***'
    RETURN
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.2: Initialize the cartesian processor grid
  !----------------------------------------------------------------------------

  CALL init_procgrid (                                                       &
        nproc, nprocx, nprocy, nprocio, .FALSE., .FALSE., .FALSE., lreorder, &
        icomm_world, igroup_world, my_world_id, icomm_compute,               &
        icomm_cart, igroup_cart, my_cart_id, my_cart_pos, my_cart_neigh,     &
        icomm_row, lcompute_pe, yzerror, izerror )

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  IF (izerror /= 0) THEN
    ierror = 1105
    yerror = ' ERROR   *** Processor grid could not be initialized ***'
    RETURN
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.3: Decomposition and Grid correspondences
  !----------------------------------------------------------------------------

  IF (lcompute_pe) THEN
    IF (my_cart_id == 0) THEN
      PRINT *, '    GRID ORGANIZATION'
    ENDIF

    IF (izdebug > 10) THEN
      PRINT *, '    Decompose COSMO-Model'
    ENDIF

    CALL decompose_lm

    ! Allocate fields for the geographical coordinates of the fine LM grid
    ALLOCATE (latlm_m(ie2lm,je2lm), lonlm_m(ie2lm,je2lm),                &
              latlm_u(ie2lm,je2lm), lonlm_u(ie2lm,je2lm),                &
              latlm_v(ie2lm,je2lm), lonlm_v(ie2lm,je2lm),                &
              STAT=nzstat)

    IF (.NOT. lgme2lm) THEN
    !!US  this should be done for all other models but GME
      ! fields for latitudes/longitudes of the fine LM grid 
      ! in the coarse grid coordinates
      ALLOCATE (lat_coarse_m(ie2lm,je2lm), lon_coarse_m(ie2lm,je2lm),    &
                lat_coarse_u(ie2lm,je2lm), lon_coarse_u(ie2lm,je2lm),    &
                lat_coarse_v(ie2lm,je2lm), lon_coarse_v(ie2lm,je2lm),    &
                STAT=nzstat)

      ! fields for the indices and weights that give the correspondence
      ! between the fine and the regular coarse grid
      ALLOCATE (i_index(ie2lm,je2lm,5), j_index(ie2lm,je2lm,5),          &
                x_wght (ie2lm,je2lm,5), y_wght (ie2lm,je2lm,5),          &
                STAT=nzstat)
    ENDIF

    IF (nzstat /= 0) THEN
      ierror    = 1101
      yerror    = ' ERROR    *** Allocation of space failed ***'
      RETURN
    ENDIF

    CALL compute_geo_coord

    IF (lgme2lm) THEN
      IF (izdebug > 10) THEN
        PRINT *, '    Decompose GME'
      ENDIF

      CALL decompose_gme
      CALL setup_xd (ig1s , ig1e , ig2s , ig2e , ids, ide, ni_gme,         &
                     igg1s, igg1e, igg2s, igg2e,                           &
                     nproc1, nproc2, my_cart_id, icomm_cart, ilim1, ilim2, &
                     yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1110
        yerror = ' ERROR   *** GME Communication could not be initialized ***'
        RETURN
      ENDIF

      CALL coor_gme_lm (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1111
        yerror = ' ERROR   *** GME grid correspondence setup failed ***'
        RETURN
      ENDIF
    ENDIF

    IF (lec2lm .OR. llm2lm .OR. lum2lm .OR. lgfs2lm .OR. lgsm2lm .OR. lhir2lm) THEN
      IF (izdebug > 10) THEN
        PRINT *, '    Decompose coarse grid'
      ENDIF

      CALL decompose_coarse_grid (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1111 * izerror
        yerror = ' ERROR   *** Coarse grid decomposition failed ***'
        RETURN
      ENDIF

      CALL coor_coarse_grid_lm (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1112
        yerror = ' ERROR   *** Coarse grid correspondence setup failed ***'
        RETURN
      ENDIF
    ENDIF

    IF (lcm2lm) THEN
      IF (izdebug > 10) THEN
        PRINT *, '    Decompose climate grid'
      ENDIF

      CALL decompose_cm (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1111 * izerror
        yerror = ' ERROR   *** Climate grid decomposition failed ***'
        RETURN
      ENDIF

      CALL coor_cm_lm (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1112
        yerror = ' ERROR   *** Climate grid correspondence setup failed ***'
        RETURN
      ENDIF
    ENDIF
  ENDIF

  ! Initialize the utility module parallel_utilities
  IF (izdebug > 10) THEN
    PRINT *, '    Initialize parallel utilities'
  ENDIF

  CALL init_par_utilities                                                    &
   (ie2lm, je2lm, kelm, ie2lm_tot, je2lm_tot, kelm_tot, ie2lm_max, je2lm_max,&
    istartpar, iendpar, jstartpar, jendpar, nproc, nprocx, nprocy, nprocio,  &
    isubpos, nboundlines, icomm_cart, my_cart_id, imp_reals, imp_integers)

  ! Open file for ASCII Debug output in my_cart_id == 0
  IF (my_cart_id == 0) THEN
    CALL get_free_unit (ndebug)
    OPEN(ndebug, FILE='YUDEBUG', FORM=  'FORMATTED', STATUS='NEW',      &
         IOSTAT=izerror)
    IF(izerror /= 0) THEN
      ierror  = 1113
      yerror  = ' ERROR    *** Error while opening file YUDEBUG *** '
      RETURN
    ENDIF
  ENDIF

  IF ( (idbg_level > 4) .AND. (num_compute > 1) .AND. (lcompute_pe)) THEN
    izbuflen = 26+nprocx+nprocy+200
    CALL check_decomposition (izbuflen, ndebug, yzerror, izerror)
    IF (izerror /= 0) THEN
      ierror = 1115
      yerror = ' ERROR   *** Decomposition check failed ***'
      RETURN
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Memory allocation and organization of Grib tables
!------------------------------------------------------------------------------

  IF (lcompute_pe) THEN
    IF (my_cart_id == 0) THEN
      PRINT *, '    MEMORY ALLOCATION'
    ENDIF
  
    ! for the output LM
    CALL alloc_lm (yzerror, izerror)
    IF (izerror /= 0) THEN
      ierror = 1301
      yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
      RETURN
    ENDIF
  
    IF (izdebug > 10) THEN
      PRINT *, '    Setup variable table for COSMO-Model'
    ENDIF

    CALL setup_vartab_lm (yzerror, izerror)
    IF (izerror /= 0) THEN
      ierror = 1302
      yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
      RETURN
    ENDIF
  
    ! for the input models
    IF (lgme2lm) THEN
      CALL alloc_gme (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1311
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
  
      IF (izdebug > 10) THEN
        PRINT *, '    Setup variable table for GME'
      ENDIF

      CALL setup_vartab_grid_in (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1312
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
    ENDIF
  
    IF (lec2lm) THEN
      CALL alloc_coarse_grid (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1311
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
  
      IF (izdebug > 10) THEN
        PRINT *, '    Setup variable table for IFS'
      ENDIF

      CALL setup_vartab_grid_in (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1312
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
    ENDIF
  
    IF (llm2lm) THEN
      CALL alloc_coarse_grid (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1311
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
  
      IF (izdebug > 10) THEN
        PRINT *, '    Setup variable table for input COSMO'
      ENDIF

      CALL setup_vartab_grid_in (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1312
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
    ENDIF

    IF (lum2lm) THEN
      CALL alloc_coarse_grid (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1311
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF

      IF (izdebug > 10) THEN
        PRINT *, '    Setup variable table for Unified Model'
      ENDIF

      CALL setup_vartab_grid_in (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1312
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
    ENDIF

    IF (lgfs2lm .OR. lgsm2lm .OR. lhir2lm) THEN
      CALL alloc_coarse_grid (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1311
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF

      IF (izdebug > 10) THEN
        PRINT *, '    Setup variable table for grib_api'
      ENDIF

      ! these models definitely have Grib2 input
      CALL setup_vartab_grid_in  (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1312
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
    ENDIF

    IF (lcm2lm) THEN
      CALL alloc_coarse_grid (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1311
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF

      IF (izdebug > 10) THEN
        PRINT *, '    Setup variable table for climate model'
      ENDIF

      CALL setup_vartab_grid_in (yzerror, izerror)
      IF (izerror /= 0) THEN
        ierror = 1312
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 5: Further initializations
!------------------------------------------------------------------------------

  IF (lasync_io .OR. (num_compute > 1) ) THEN
    CALL mpe_io_init (icomm_world, num_compute, num_io, iztest_compute,   &
                      izerror)
    IF (izerror /= 0) THEN
      ierror  = 1205
      yerror  = ' ERROR   *** MPE_IO initialization failed ***'
      RETURN
    ENDIF
  ENDIF

  IF (lcompute_pe .AND. (lasync_io .OR. (num_compute > 1) )) THEN
    CALL mpe_io_reconfig (icomm_compute, izerror)
    IF (izerror /= 0) THEN
      ierror  = 1205
      yerror  = ' ERROR   *** MPE_IO communicators need reconfiguration ***'
      RETURN
    ENDIF
  ENDIF

  ! Allocate the sendbuffer with the maximal size.
  ! also for non-parallel runs, because it is passed as argument
!US  IF (lcompute_pe .AND. num_compute > 1) THEN
    isendbuflen =                                                            &
     (MAX(ie2lm_tot/nprocx+1+2*nboundlines,je2lm_tot/nprocy+1+2*nboundlines) &
           *nboundlines*(kelm_tot+1)) * 24
    ALLOCATE (sendbuf(isendbuflen,8) , STAT=nzstat )
    sendbuf(:,:) = 0.0_ireals
!US  ENDIF

  IF (my_world_id == 0) THEN
    ! close OUTPUT and open it later on the compute PEs
    CLOSE (noutput, STATUS='KEEP')
    CALL release_unit (noutput)
  ENDIF

!------------------------------------------------------------------------------
! End of external subroutine setup_int2lm
!------------------------------------------------------------------------------

END SUBROUTINE setup_int2lm
