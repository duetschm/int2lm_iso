!+ Source module for cleaning up at the end of the program
!==============================================================================

MODULE src_cleanup

!==============================================================================
!
! Description:
!   This module cleans up at the end of the program. Special tasks are
!   - deallocation of memory
!   - cleaning up the environment (parallel)
!   - collect and print the timings, if required
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
! V1_5         2007/07/09 Ulrich Schaettler
!  Added fields qr_lm, qs_lm, qg_lm, qr_in, qs_in, qg_in
!  Eliminated akhlm, bkhlm
!  Added additional chemistry fields
! V1_7         2007/11/26 Ulrich Schaettler
!  Deallocate monthly and actual ndvi rates, if allocated
! V1_8         2008/05/29 Ulrich Schaettler
!  Renamed all gz0-variables to z0 and only work with z0-values
!  Deallocate additional field rootdp_mx
! V1_9         2009/09/03 Ulrich Schaettler
!  Eliminated GME fields for ndvi-values
!  Adapted deallocation of memory at the end of the program
! V1_10        2009/12/17 Oliver Fuhrer, Ulrich Schaettler
!  Added field p_in
!  Adaptations for lum2lm
!  Deallocation of fields for lake_coldstart
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler
!  Deallocate qv_2m_in, t_2m_in, qv_2m_gl, t_2m_gl and lake-values
! V1_19        2012/06/06 Ulrich Schaettler, Juergen Helmert
!  Renamed ak_in_uv, bk_in_uv to ak_in_rho, bk_in_rho according to UM conventions
!  Deallocation albedo fields for option itype_albedo=3
! V1_20        2012/09/03 Anne Roches
!  Allocate ztotal_times on all PEs, because it is used as subroutine argument
!  in all PEs
! V1_22        2013/07/11 Ulrich Schaettler, KIT
!  Removed variables vcoord_in, vcoord, aklm, bklm, sigmalm_sl
!  Removed ART parts, which are now in an own component (KIT)
!
!  Added fields for water isotope simulation (Stephan Pfahl)
!  Added field for height correction of soil isotopes. Hui Tang 2013-11-20
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters , ONLY :   &
    ireals,     & ! kind-type parameter for "normal" integer variables
    iintegers     ! KIND-type parameters for real variables

!------------------------------------------------------------------------------

USE data_fields_lm,   ONLY :   &
  fis_lm    ,      & ! orography * G                                    (m2/s2)
  fis_gl    ,      & ! orography * G interpolated from GME              (m2/s2)
  hsurf_lm  ,      & ! orography                                        (  m  )
  fr_land_lm,      & ! land fraction of grid element                    (  1  )
  fr_lake_lm,      & ! lake fraction of grid element                    (  1  )
  depth_lk_lm,     & ! lake depth                                       (  m  )
  salt_lk_lm,      & ! lake salinity                                    ( g/kg)
  z0_lm     ,      & ! roughness length                                 (  m  )
  soiltyp_lm,      & ! type of the soil (keys 0-9)                      (  1  )
  plcov_lm  ,      & ! fraction covered by plants                       (  1  )
  lai_lm    ,      & ! leaf area index                                  (  1  )
  rootdp_lm ,      & ! depth of the roots                               (  m  )
  rootdp_mx ,      & ! depth of the roots from external parameters      (  m  )
  t_s_lm    ,      & ! temperature of the ground surface                (  K  )
  t_snow_lm ,      & ! temperature of the snow surface                  (  K  )
  t_m_lm    ,      & ! temperature between upper and medium soil layer  (  K  )
  t_cl_lm   ,      & ! temperature between medium and lower soil layer  (  K  )
  qv_s_lm   ,      & ! specific water vapor content on the surface      (kg/kg)
  w_snow_lm ,      & ! water content of the snow                        (m H2O)
  w_i_lm    ,      & ! water content of the interception storage        (m H2O)
  w_g1_lm   ,      & ! water content of the upper soil layer            (m H2O)
  w_g2_lm   ,      & ! water content of the medium soil layer           (m H2O)
  w_g3_lm   ,      & ! water content of the lower soil layer            (m H2O)
                     ! (if nl_soil_lm = 3, unused otherwise)
  w_cl_lm   ,      & ! climatological deep soil water content           (m H2O)
  hmo3_lm   ,      & ! height of maximum ozone concentration            ( Pa  )
  vio3_lm   ,      & ! total vertically integrated ozone content        (Pa O3)
  lmask_lm  ,      & ! mask of points on the frame
  latlm_m   ,      & ! latitudes of the LM grid points
  lonlm_m   ,      & ! longitudes of the LM grid points
  latlm_u   ,      & ! latitudes of the LM u grid points
  lonlm_u   ,      & ! longitudes of the LM u grid points
  latlm_v   ,      & ! latitudes of the LM v grid points
  lonlm_v   ,      & ! longitudes of the LM v grid points
  ps_gl     ,      & ! surface pressure on the interpol. GME orogr.     ( Pa  )
  p0_gl     ,      & ! ref. pres. on full levels + interpol. COARSE LM oro.(Pa)
  dp0_gl    ,      & ! reference pressure thickness of layers           ( Pa  )
  rho0_gl   ,      & ! reference density at the full model levels       (kg/m3)
  hhl_gl    ,      & ! height of half-levels on the interpol. COARSE LM oro.(m)
  hsurf_gl  ,      & ! height of orography interpolated from coarse grid(  m  )
  hsurfs_gl          ! interpolated splitted parts of coarse topo       (  m  )

USE data_fields_lm,   ONLY :   &
  for_e_lm     ,   & ! ground fraction covered by evergreen forest      (  -  )
  for_d_lm     ,   & ! ground fraction covered by deciduous forest      (  -  )
  sso_stdh_lm  ,   & ! standard deviation of subgrid scale orography    (  m  )
  sso_gamma_lm ,   & ! anisotropy of the orography                      (  -  )
  sso_theta_lm ,   & ! angle betw. principal axis of orography and E    ( rad )
  sso_sigma_lm ,   & ! mean slope of subgrid scale orography            (  -  )
  skyview_lm   ,   & ! sky view                                         (  1  )
  slo_asp_lm   ,   & ! slope aspect                                     ( rad )
  slo_ang_lm   ,   & ! slope angle                                      ( rad )
  horizon_lm   ,   & ! horizon                                          ( rad )
  emis_rad_lm  ,   & ! thermal radiative surface emissivity             (  1  )
  prs_min_lm   ,   & ! minimum stomata resistance of plants             ( s/m )
  alb_dif12_lm ,   & ! solar surface albedo - diffuse                   (  1  )
  alb_dif_lm   ,   & ! solar surface albedo - diffuse                   (  1  )
  ndvi_mrat_lm ,   & !ratio of monthly mean normalized differential     (  1  )
                     ! vegetation index to annual maximum for 12 months
  ndviratio_lm ,   & ! actual value of ndvi for a special day           (  1  )
  aer_su12_lm  ,   & ! Tegen (1997) aerosol type sulfate drops          (  -  )
  aer_du12_lm  ,   & ! Tegen (1997) aerosol type mineral dust coarse    (  -  )
  aer_or12_lm  ,   & ! Tegen (1997) aerosol type organic(water solub)   (  -  )
  aer_bc12_lm  ,   & ! Tegen (1997) aerosol type black carbon           (  -  )
  aer_ss12_lm  ,   & ! Tegen (1997) aerosol type sea salt               (  -  )
  aer_su_lm    ,   & ! Tegen (1997) aerosol type sulfate drops          (  -  )
  aer_du_lm    ,   & ! Tegen (1997) aerosol type mineral dust coarse    (  -  )
  aer_or_lm    ,   & ! Tegen (1997) aerosol type organic(water solub)   (  -  )
  aer_bc_lm    ,   & ! Tegen (1997) aerosol type black carbon           (  -  )
  aer_ss_lm    ,   & ! Tegen (1997) aerosol type sea salt               (  -  )
  plcov12      ,   & ! monthly climatology of plant cover               (  1  )
  z012         ,   & ! monthly climatology of roughness length          (  m  )
  lai12        ,   & ! monthly climatology of leaf area index           (  1  )
  fic_gl       ,   & ! check level of geopotential                      (m2/s2)
  t_2m_gl      ,   & ! 2m temperature                                   (  K  )
  qv_2m_gl     ,   & ! 2m humidity                                      (kg/kg)
  dpsdt_gl     ,   & ! surface pressure tendency                        (Pa/s )
  t_s_gl       ,   & ! temperature of the ground surface                (  K  )
  t_skin_gl    ,   & ! skin temperature of the ground surface           (  K  )
  rh_s_gl      ,   & ! relative humidity at the surface                 (kg/kg)
  dtms_gl      ,   & ! t_m_lm    - t_s_lm                               (  K  )
  dtkes_gl     ,   & ! t(ke)_lm  - t_s_lm                               (  K  )
  dtssnow_gl         ! t_s_lm    - t_snow_lm                            (  K  )

USE data_fields_lm,   ONLY :   &
  u_lm      ,      & ! zonal wind speed                                 ( m/s )
  v_lm      ,      & ! meridional wind speed                            ( m/s )
  w_lm      ,      & ! vertical wind speed (defined on half levels)     ( m/s )
  t_lm      ,      & ! temperature                                      (  K  )
  qv_lm     ,      & ! specific water vapor content                     (kg/kg)
  qc_lm     ,      & ! specific cloud water content                     (kg/kg)
  qi_lm     ,      & ! cloud ice content                                (kg/kg)
  qr_lm     ,      & ! rain      content                                (kg/kg)
  qs_lm     ,      & ! snow      content                                (kg/kg)
  qg_lm     ,      & ! graupel   content                                (kg/kg)
  p0_lm     ,      & ! reference pressure                               ( Pa  )
  dp0_lm    ,      & ! reference pressure thickness of layers           ( Pa  )
  rho0_lm   ,      & ! reference density at the full model levels       (kg/m3)
  t0_lm     ,      & ! reference temperature                            (  K  )
  p_lm      ,      & ! full pressure (needed for lum2lm)                ( Pa  )
  pp_lm     ,      & ! deviation from the reference pressure            ( Pa  )
  hhl_lm    ,      & ! height of half-levels of LM                      (  m  )
  lolp_lm   ,      & ! Land Sea Mask of LM for 'M'atch Interpolation
  w_intpol  ,      & ! interpolation weights for gme2lm                 (  -  )
  n_intpol     ,   & ! nearest GME gridpoint for nearest neighbor interpolation
  m_intpol     ,   & ! nearest GME gridpoint with same lsm for match interpolation
  l_intpol     ,   & ! to use a far away GME gridpoint with same lsm
  cgas_lm   ,      & !
  caero_lm  ,      & !
! iso code
  riso_lm   ,      & ! isotope ratios in water vapor
  risosoil_lm,     & ! isotope ratios in soil water
! Hui Tang 2013-11-20
  drisoke_gl          ! riso_lm(:,:,ke,1-2) - risosoil_lm (:,:,1-2) (1: 18O; 2: 2H);
! end iso code
! end iso code

USE data_fields_lm, ONLY : &
  t_mnw_lk_lm  ,   & ! mean temperature of the water column          (  K  )
  t_wml_lk_lm  ,   & ! mixed-layer temperature                       (  K  )
  t_bot_lk_lm  ,   & ! temperature at the water-bottom sediment
                     ! interface                                     (  K  )
  t_b1_lk_lm   ,   & ! temperature at the bottom of the upper layer
                     ! of the sediments                              (  K  )
  c_t_lk_lm    ,   & ! shape factor with respect to the
                     ! temperature profile in lake thermocline       (  -  )
  h_ml_lk_lm   ,   & ! thickness of the mixed-layer                  (  m  )
  h_b1_lk_lm   ,   & ! thickness of the upper layer
                     ! of bottom sediments                           (  m  )
  t_ice_lm     ,   & ! temperature of ice/water surface              (  K  )
  h_ice_lm           ! lake/sea ice thickness                        (  m  )

!------------------------------------------------------------------------------

USE data_grid_in,      ONLY :   &
  ak_in ,      & ! vertical coordinate parameters for half levels
  bk_in ,      & !                  - " -
  akh_in,      & ! vertical coordinate parameters for main levels
  bkh_in,      & !                  - " -
  akh_in_rho,  & ! coefficients for main levels for u, v (lum2lm)
  bkh_in_rho,  & ! coefficients for main levels for u, v (lum2lm)
  dak_in,      & ! difference of coordinate parameters
  dbk_in,      & !                  - " -
  ilim1,       & ! decomposition limits in x-direction (formerly 1)
  ilim2          ! decomposition limits in y-direction (formerly 2)

!------------------------------------------------------------------------------

USE data_fields_in,    ONLY :   &
 fis_gme            ,   & ! orography * g                               (m2/s2)
 soiltyp_gme        ,   & ! type of the soil (keys 0-9)                   --
 ps_gme             ,   & ! surface pressure                            ( Pa  )
 t_s_gme            ,   & ! temperature of the ground surface           (  K  )
 t_snow_gme         ,   & ! temperature of the snow-surface             (  K  )
 t_m_gme            ,   & ! temp. between upper and medium soil layer   (  K  )
 qv_s_gme           ,   & ! specific water vapor content on the surface (kg/kg)
 w_g1_gme           ,   & ! water content of the upper soil layer       (m H2O)
 w_g2_gme           ,   & ! water content of the medium soil layer      (m H2O)
 w_g3_gme           ,   & ! water content of the deepest soil layer     (m H2O)
 u_gme              ,   & ! zonal wind speed                            ( m/s )
 v_gme              ,   & ! meridional wind speed                       ( m/s )
 t_gme              ,   & ! temperature                                 (  K  )
 qv_gme             ,   & ! specific water vapor content                (kg/kg)
 qc_gme             ,   & ! specific cloud water content                (kg/kg)
 fr_land_gme        ,   & ! land fraction of grid element               (  1  )
 root_gme           ,   & ! depth of the roots                          (  m  )
 vio3_gme           ,   & ! total vertically integrated ozone content   (Pa O3)
 hmo3_gme           ,   & ! height of maximum ozone concentration       ( Pa  )
 t_g_gme            ,   & ! temperature                                 (  K  )
 t_cl_gme           ,   & ! temp.  between medium and lower soil layer  (  K  )
 w_snow_gme         ,   & ! water content of the snow                   (m H2O)
 w_i_gme            ,   & ! water content of the interception storage   (m H2O)
 w_cl_gme           ,   & ! climatological deep soil water content      (m H2O)
 dpsdt_gme          ,   & ! surface pressure tendency                   (Pa/s )
 lolp_gme                 ! Land Sea Mask of GME for 'M'atch Interpolation

USE data_fields_in,    ONLY :   &
 fis_in      ,          & ! orography * g                               (m2/s2)
 hsurf_in    ,          & ! orography                                   (  m  )
 soiltyp_in  ,          & ! type of the soil (keys 0-9)                   --
 ps_in       ,          & ! surface pressure                            ( Pa  )
 t_s_in      ,          & ! temperature of the ground surface           (  K  )
 t_2m_in     ,          & ! 2m temperature                              (  K  )
 t_snow_in   ,          & ! temperature of the snow-surface             (  K  )
 t_g1_in     ,          & ! temperature of first soil layer             (  K  )
 t_g2_in     ,          & ! temperature of second soil layer            (  K  )
 t_g3_in     ,          & ! temperature of third soil layer             (  K  )
 qv_s_in     ,          & ! specific water vapor content on the surface (kg/kg)
 qv_2m_in    ,          & ! specific water vapor content in 2m          (kg/kg)
 w_g1_in     ,          & ! water content of the upper soil layer       (m H2O)
 w_g2_in     ,          & ! water content of the medium soil layer      (m H2O)
 w_g3_in     ,          & ! water content of the deepest soil layer     (m H2O)
 t_ke_in     ,          & ! temperature lowest layer                    (  K  )
 grh_in      ,          & ! generalized relative humidity at one level  (  %  )
 lolp_in     ,          & ! Land Sea Mask of input for 'M'atch Interpolation
 u_in        ,          & ! zonal wind speed                            ( m/s )
 v_in        ,          & ! meridional wind speed                       ( m/s )
 w_in        ,          & ! vertical wind speed                         ( m/s )
 t_in        ,          & ! temperature                                 (  K  )
 p_in        ,          & ! standard pressure                           ( Pa  )
 pp_in       ,          & ! deviation from standard pressure            ( Pa  )
 qv_in       ,          & ! specific water vapor content                (kg/kg)
 qc_in       ,          & ! specific cloud water content                (kg/kg)
 qi_in       ,          & ! specific cloud ice content                  (kg/kg)
 qr_in       ,          & ! specific rain      content                  (kg/kg)
 qs_in       ,          & ! specific snow      content                  (kg/kg)
 qg_in       ,          & ! specific graupel   content                  (kg/kg)
 p0_in       ,          & ! reference pressure for coarse LM grid       ( Pa  )
 hhl_in      ,          & ! height of half levels for coarse LM grid    ( Pa  )
 fr_land_in  ,          & ! land fraction of grid element               (  1  )
 root_in     ,          & ! depth of the roots                          (  m  )
 vio3_in     ,          & ! total vertically integrated ozone content   (Pa O3)
 hmo3_in                  ! height of maximum ozone concentration       ( Pa  )

USE data_fields_in,    ONLY :   &
 t_g_in      ,          & ! temperature                                 (  K  )
 t_m_in      ,          & ! temperature                                 (  K  )
 t_cl_in     ,          & ! temp.  between medium and lower soil layer  (  K  )
 w_snow_in   ,          & ! water content of the snow                   (m H2O)
 w_i_in      ,          & ! water content of the interception storage   (m H2O)
 w_cl_in     ,          & ! climatological deep soil water content      (m H2O)
 dpsdt_in    ,          & ! surface pressure tendency                   (Pa/s )
 fic_in      ,          & ! control level for geopotential              (m2/s2)
 cgas_in     ,          & !
 caero_in    ,          & !
! iso code
 riso_in     ,          & ! isotope ratios in water vapor
 risosoil_in              ! isotope ratios in soil water
! end iso code

!------------------------------------------------------------------------------

USE data_int2lm_io,          ONLY :   &
    nuchkdat,          & ! checking the I/O data
    var_lm,            & ! variable table for LM
    var_in               ! variable table for input model

!------------------------------------------------------------------------------

USE data_int2lm_control,     ONLY :   &
    lgme2lm,      & ! if .TRUE., gme->lm, if .FALSE. gme->hm
    lec2lm,       & ! if .TRUE., ec ->lm
    llm2lm,       & ! if .TRUE., lm ->lm
    llake,        & ! if .TRUE., run with lake  
    llake_coldstart,& ! if .TRUE., initialize prognostic lake variables for cold start
    ndebug,       & ! unit number for file with debug information
    noutput,      & ! unit number for output file
    itype_ndvi,   & ! to choose treatment of surface parameters (plcov, lai)
    itype_aerosol,& ! to choose treatment of surface parameters (plcov, lai)
    itype_albedo ,& ! to choose treatment of surface albedo
    linitial,     & ! if .TRUE., initial data for LM
    lboundaries,  & ! if .TRUE., lateral boundaries for LM
    lclock,       & ! if .TRUE., system clock is present
    ltime,        & ! detailled timings of the program are given
    ltime_mean,   & ! if .TRUE., mean values of the timings are printed
    ltime_proc,   & ! if .TRUE., timings for each processor are printed
    timings,      & ! for storing the times for different parts of the program
! iso code
    liso            ! if .TRUE., include variables for water isotope simulation
! end iso code

!------------------------------------------------------------------------------

USE data_int2lm_parallel,    ONLY :   &
    num_compute,     & ! number of compute PEs
    my_cart_id,      & ! rank of this subdomain in the global communicator
    icomm_cart,      & ! communicator for the virtual cartesian topology
    isubpos,         & ! positions of the subdomains in the total domain
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    sendbuf            ! sending buffer for boundary exchange:

!------------------------------------------------------------------------------

USE environment,        ONLY :   model_abort

!------------------------------------------------------------------------------

USE parallel_utilities, ONLY :   gather_values, remark

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Public and Private Subroutines

PUBLIC   org_cleanup

PRIVATE  collect_timings, free_memory

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in src_cleanup for cleaning up
!------------------------------------------------------------------------------

SUBROUTINE org_cleanup

!------------------------------------------------------------------------------
!
! Description:
!  org_cleanup is the driver routine for cleaning up at the end of the
!  interpolation program. The main tasks are:
!   - Collect timings
!   - Free memory
!   - finalize the environment
!
! Output files:
!  File YUTIMING containing the timings for the program.
!
!------------------------------------------------------------------------------
!
! Local scalars:

LOGICAL                    ::       &
  lzopen             ! to check whether files are still open

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Collect the timings
!------------------------------------------------------------------------------

  IF (ltime) THEN
    CALL collect_timings
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Close files that are still open
!------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
    ! Close files which are still open
    INQUIRE (FILE='YUCHKDAT', OPENED=lzopen)
    IF (lzopen) CLOSE (nuchkdat, STATUS='KEEP')

    INQUIRE (FILE='YUDEBUG' , OPENED=lzopen)
    IF (lzopen) CLOSE (ndebug, STATUS='KEEP')

    INQUIRE (FILE='OUTPUT' , OPENED=lzopen)
    IF (lzopen) CLOSE (noutput, STATUS='KEEP')
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Deallocate the memory
!------------------------------------------------------------------------------

  CALL free_memory

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE org_cleanup

!==============================================================================
!+ Subroutine that collects the timings from all nodes and prints them
!------------------------------------------------------------------------------

SUBROUTINE collect_timings

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine collects the timings from all nodes (if running in 
!   parallel mode) and prints them to an ASCII file YUTIMING. The measured
!   real times for different program parts are given for every file processed.
!
! Method:
!   In parallel mode collect values from all nodes with gather_values.
!
!------------------------------------------------------------------------------

! Declarations:
! Local variables

  REAL (KIND=ireals), ALLOCATABLE       ::       &
    ztotal_times(:,:)      ! for gathering the times from all nodes

  REAL (KIND=ireals)                    ::       &
    zsetup  (4), zgridgen(4), zreadext(4), zreadgme(4), zunpack (4),      &
    zinterpo(4), zcomm1  (4), zlmfield(4), zwritelm(4), zrest   (4),      &
    zsum1   (4), zsum2   (4)

  INTEGER (KIND=iintegers)   ::       &
    izerror, niostat, istat,                                  &
    nutiming, npr, nfpr, nlpr, nzsendcount, nzrecvcount,      &
    n

  CHARACTER (LEN= 8) yutiming
  CHARACTER (LEN=25) yzroutine
  CHARACTER (LEN=75) yzerrmsg

!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE collect_timings
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations and Computing sums
!------------------------------------------------------------------------------

  ! Initializations
  yzroutine   = 'collect_timings'
  yzerrmsg    = '  '
  izerror     = 0
  nzsendcount = 40
  nzrecvcount = nzsendcount

  ! Allocate the buffers (has to be done in all PEs)
  ALLOCATE ( ztotal_times (40, 0:num_compute-1), STAT=istat )
  ztotal_times(:,:) = 0.0_ireals

  IF (istat /= 0) THEN
    yzerrmsg = 'Buffer allocation failed'
    CALL model_abort (my_cart_id, 1011, yzerrmsg, yzroutine)
  ENDIF

  ! Compute sums:
  ! timings(1) = sum over all other timings
  DO n = 2, 40
    timings(1) = timings(1) + timings(n)
  ENDDO

  ! timings(5) = sum (6,7,8,9)
  timings(5) = timings(6) + timings(7) + timings(8) + timings(9)

!------------------------------------------------------------------------------
!- Section 2:  Gather the values from all nodes
!------------------------------------------------------------------------------

  IF (num_compute > 1) THEN
    CALL gather_values (timings, ztotal_times, nzsendcount, num_compute,     &
                        imp_reals, 0, icomm_cart, yzerrmsg, izerror)

    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 10531, 'gather_values', yzroutine)
    ENDIF
  ELSE
    ztotal_times(:,0) = timings(:)
  ENDIF

!------------------------------------------------------------------------------
!- Section 3: Print the data to the file YUTIMING
!------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
    nutiming = 9
    yutiming = 'YUTIMING'

    OPEN(nutiming, FILE=yutiming, FORM=  'FORMATTED', STATUS='UNKNOWN',  &
         IOSTAT=niostat)
    IF(niostat /= 0) THEN
      yzerrmsg = ' ERROR    *** opening file YUTIMING failed *** '
      CALL model_abort (my_cart_id, 9702, yzerrmsg, yzroutine)
    ENDIF

    IF (ltime_proc) THEN    ! type of timing
      ! Print a headline in file YUDEBUG
      WRITE (nutiming, '(A)')  '  '
      WRITE (nutiming, '(A)')                                                 &
                    '      Detailed timings per processor:'
      WRITE (nutiming, '(A)')                                                 &
                    '      ==============================='
      WRITE (nutiming, '(A)')  '  '

      ! Print the information from all processes
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Total time for INT2LM:        ', ztotal_times(1,0)
      WRITE (nutiming, '(A,F9.2)')                                            &
        '  Time for the setup of INT2LM:          ', ztotal_times(2,0)
      WRITE (nutiming, '(A,F9.2)')                                            &
        '  Time for reading external parameters:  ', ztotal_times(3,0)
      WRITE (nutiming, '(A)')  '  '

      DO npr = 0 , num_compute-1 , 8
        nfpr = npr
        nlpr = MIN (nfpr+7 , num_compute-1)
 
        WRITE (nutiming, '(A,8F7.2)')                                       &
          'Interpolation:          ',(ztotal_times( 5,n),n=nfpr,nlpr)
        WRITE (nutiming, '(A,8F7.2)')                                       &
          '  Reading coarse files: ',(ztotal_times( 6,n),n=nfpr,nlpr)
        WRITE (nutiming, '(A,8F7.2)')                                       &
          '  Unpack coarse records:',(ztotal_times( 7,n),n=nfpr,nlpr)
        WRITE (nutiming, '(A,8F7.2)')                                       &
          '  Interpolation:        ',(ztotal_times( 8,n),n=nfpr,nlpr)
        WRITE (nutiming, '(A,8F7.2)')                                       &
          '  Add. Communication:   ',(ztotal_times( 9,n),n=nfpr,nlpr)
        WRITE (nutiming, '(A,8F7.2)')                                       &
          'Computing LM fields     ',(ztotal_times(10,n),n=nfpr,nlpr)
        WRITE (nutiming, '(A,8F7.2)')                                       &
          'Writing LM fields       ',(ztotal_times(11,n),n=nfpr,nlpr)
        WRITE (nutiming, '(A,8F7.2)')                                       &
          'All the rest            ',(ztotal_times(12,n),n=nfpr,nlpr)
        WRITE (nutiming, '(A)')  '  '
        WRITE (nutiming, '(A)')  '  '
      ENDDO

    ELSEIF (ltime_mean) THEN    

      ! Print a headline in file YUDEBUG
      WRITE (nutiming, '(A)')  '  '
      WRITE (nutiming, '(A)')                                                 &
                    '      Mean values of the timings for all processors:'  
      WRITE (nutiming, '(A)')                                                 &
                    '      =============================================='
      WRITE (nutiming, '(A)')  '  '

      ! Print the information from all processes
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Total time for INT2LM:        ', ztotal_times(1,0)
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the setup of INT2LM:          ', ztotal_times(2,0)
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for reading external parameters:  ', ztotal_times(3,0)
      WRITE (nutiming, '(A)')  '  '

      ! Compute the meanvalues
      zsum1   (:) = 0.0_ireals
      zsum2   (:) = 0.0_ireals
      zsetup  (:) = 0.0_ireals
      zgridgen(:) = 0.0_ireals
      zreadext(:) = 0.0_ireals
      zreadgme(:) = 0.0_ireals
      zunpack (:) = 0.0_ireals
      zinterpo(:) = 0.0_ireals
      zcomm1  (:) = 0.0_ireals
      zlmfield(:) = 0.0_ireals
      zwritelm(:) = 0.0_ireals
      zrest   (:) = 0.0_ireals

      zsum1   (1) = MINVAL (ztotal_times( 1,:))
      zsetup  (1) = MINVAL (ztotal_times( 2,:))
      zgridgen(1) = MINVAL (ztotal_times( 3,:))
      zreadext(1) = MINVAL (ztotal_times( 4,:))
      zsum2   (1) = MINVAL (ztotal_times( 5,:))
      zreadgme(1) = MINVAL (ztotal_times( 6,:))
      zunpack (1) = MINVAL (ztotal_times( 7,:))
      zinterpo(1) = MINVAL (ztotal_times( 8,:))
      zcomm1  (1) = MINVAL (ztotal_times( 9,:))
      zlmfield(1) = MINVAL (ztotal_times(10,:))
      zwritelm(1) = MINVAL (ztotal_times(11,:))
      zrest   (1) = MINVAL (ztotal_times(12,:))

      zsum1   (3) = MAXVAL (ztotal_times( 1,:))
      zsetup  (3) = MAXVAL (ztotal_times( 2,:))
      zgridgen(3) = MAXVAL (ztotal_times( 3,:))
      zreadext(3) = MAXVAL (ztotal_times( 4,:))
      zsum2   (3) = MAXVAL (ztotal_times( 5,:))
      zreadgme(3) = MAXVAL (ztotal_times( 6,:))
      zunpack (3) = MAXVAL (ztotal_times( 7,:))
      zinterpo(3) = MAXVAL (ztotal_times( 8,:))
      zcomm1  (3) = MAXVAL (ztotal_times( 9,:))
      zlmfield(3) = MAXVAL (ztotal_times(10,:))
      zwritelm(3) = MAXVAL (ztotal_times(11,:))
      zrest   (3) = MAXVAL (ztotal_times(12,:))

      DO n=0,num_compute-1
        zsum1   (4) = zsum1   (4) + ztotal_times( 1,n)
        zsetup  (4) = zsetup  (4) + ztotal_times( 2,n)
        zgridgen(4) = zgridgen(4) + ztotal_times( 3,n)
        zreadext(4) = zreadext(4) + ztotal_times( 4,n)
        zsum2   (4) = zsum2   (4) + ztotal_times( 5,n)
        zreadgme(4) = zreadgme(4) + ztotal_times( 6,n)
        zunpack (4) = zunpack (4) + ztotal_times( 7,n)
        zinterpo(4) = zinterpo(4) + ztotal_times( 8,n)
        zcomm1  (4) = zcomm1  (4) + ztotal_times( 9,n)
        zlmfield(4) = zlmfield(4) + ztotal_times(10,n)
        zwritelm(4) = zwritelm(4) + ztotal_times(11,n)
        zrest   (4) = zrest   (4) + ztotal_times(12,n)
      ENDDO

      zsum1   (2) = zsum1   (4) / num_compute              
      zsetup  (2) = zsetup  (4) / num_compute              
      zgridgen(2) = zgridgen(4) / num_compute              
      zreadext(2) = zreadext(4) / num_compute              
      zsum2   (2) = zsum2   (4) / num_compute              
      zreadgme(2) = zreadgme(4) / num_compute              
      zunpack (2) = zunpack (4) / num_compute              
      zinterpo(2) = zinterpo(4) / num_compute              
      zcomm1  (2) = zcomm1  (4) / num_compute              
      zlmfield(2) = zlmfield(4) / num_compute
      zwritelm(2) = zwritelm(4) / num_compute
      zrest   (2) = zrest   (4) / num_compute

      WRITE (nutiming, '(A,A)') '                                 ',        &
                             '     Min.     Avg.     Max.    Total'
      WRITE (nutiming, '(A,4F9.2)')                                       &
        'Interpolation:                   ',zsum2    
      WRITE (nutiming, '(A,4F9.2)')                                       &
        '  Reading coarse files:          ',zreadgme
      WRITE (nutiming, '(A,4F9.2)')                                       &
        '  Unpack coarse records:         ',zunpack
      WRITE (nutiming, '(A,4F9.2)')                                       &
        '  Interpolation:                 ',zinterpo
      WRITE (nutiming, '(A,4F9.2)')                                       &
        '  Add. Communication:            ',zcomm1
      WRITE (nutiming, '(A,4F9.2)')                                       &
        'Computing LM fields              ',zlmfield
      WRITE (nutiming, '(A,4F9.2)')                                       &
        'Writing LM fields                ',zwritelm
      WRITE (nutiming, '(A,4F9.2)')                                       &
        'All the rest                     ',zrest

    ENDIF    ! type of timing

  ENDIF  ! printing in my_cart_id = 0

!------------------------------------------------------------------------------
!- Section 4: Cleanup
!------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
    CLOSE (nutiming, STATUS='KEEP')
  ENDIF

  DEALLOCATE ( ztotal_times, STAT=istat)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  
END SUBROUTINE collect_timings

!==============================================================================
!+ Deallocates the allocated fields at the end of the program.
!------------------------------------------------------------------------------

SUBROUTINE free_memory

!------------------------------------------------------------------------------
!
! Description:
!   This routine deallocates the allocated fields at the end of the program.
!
! Method:
!   DEALLOCATE statement
!
!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers) :: istat              ! for local error-code

!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE free_memory
!------------------------------------------------------------------------------

istat = 0

! Vertical coordinate parameters
!-------------------------------
  DEALLOCATE (ak_in,bk_in, akh_in,bkh_in, dak_in,dbk_in, STAT=istat)
  IF (ALLOCATED(akh_in_rho)) THEN
    DEALLOCATE (akh_in_rho, bkh_in_rho, STAT=istat)
  ENDIF

! Organization variables for parallel program
!--------------------------------------------
  DEALLOCATE (ilim1, ilim2, timings, isubpos, sendbuf, STAT=istat)

! Variables for GME and LM variable table
!----------------------------------------
  DEALLOCATE (var_lm, var_in, STAT=istat)

! LM fields from reading external parameters
!-------------------------------------------
  DEALLOCATE (hhl_lm    , p0_lm     , dp0_lm    , rho0_lm   , t0_lm     , &
                                                              STAT=istat)
  DEALLOCATE (fis_lm    , hsurf_lm  , fr_land_lm, z0_lm     , soiltyp_lm, &
              plcov_lm  , lai_lm    , rootdp_lm , rootdp_mx , lolp_lm   , &
              vio3_lm   , hmo3_lm   ,                         STAT=istat)

  IF (llake) THEN
    DEALLOCATE (fr_lake_lm, depth_lk_lm, salt_lk_lm,          STAT=istat)

    IF (llake_coldstart) THEN
      DEALLOCATE (t_mnw_lk_lm, t_wml_lk_lm, t_bot_lk_lm, t_b1_lk_lm,      &
                  c_t_lk_lm,   h_ml_lk_lm,  h_b1_lk_lm,                   &
                  t_ice_lm,    h_ice_lm,                 STAT=istat)
    ENDIF
  ENDIF

  IF (itype_aerosol == 2) THEN
    DEALLOCATE (aer_su12_lm, aer_du12_lm, aer_or12_lm, aer_bc12_lm, aer_ss12_lm)
    DEALLOCATE (aer_su_lm, aer_du_lm, aer_or_lm, aer_bc_lm, aer_ss_lm)
  ENDIF

  IF (itype_albedo == 3) THEN
    DEALLOCATE (alb_dif12_lm, alb_dif_lm)
    ! the fields from option itype_albedo=2 are only used for initial fields
    ! they are deallocated in int2lm_org after the initial fields-loop
  ENDIF

  IF     (itype_ndvi == 1) THEN
    DEALLOCATE (ndvi_mrat_lm, ndviratio_lm)
  ELSEIF (itype_ndvi == 2) THEN
    DEALLOCATE (plcov12, lai12, z012)
  ENDIF

! LM fields from the interpolation
!---------------------------------
  DEALLOCATE (u_lm      , v_lm      , w_lm      , t_lm      , qv_lm     , &
              qc_lm     , pp_lm     ,                         STAT=istat)

  IF (ALLOCATED(qi_lm)) THEN
    DEALLOCATE (qi_lm   ,                                     STAT=istat)
  ENDIF
  IF (ALLOCATED(qr_lm)) THEN
    DEALLOCATE (qr_lm   ,                                     STAT=istat)
  ENDIF
  IF (ALLOCATED(qs_lm)) THEN
    DEALLOCATE (qs_lm   ,                                     STAT=istat)
  ENDIF
  IF (ALLOCATED(qg_lm)) THEN
    DEALLOCATE (qg_lm   ,                                     STAT=istat)
  ENDIF

  IF (ALLOCATED( p_lm)) THEN
    DEALLOCATE ( p_lm   ,                                     STAT=istat)
  ENDIF

  DEALLOCATE (fis_gl    , ps_gl     , fic_gl    , dpsdt_gl  , t_s_gl    , &
              dtkes_gl  , rh_s_gl   , dtssnow_gl, dtms_gl   , t_2m_gl   , &
              qv_2m_gl  ,                                     STAT=istat)

  IF (ALLOCATED(dp0_gl)) THEN
    ! lm2lm case
    DEALLOCATE (hhl_gl  , p0_gl     , dp0_gl    , rho0_gl   , STAT=istat)
    DEALLOCATE (hsurfs_gl                                   , STAT=istat)
  ENDIF

  IF (ALLOCATED(hsurf_gl)) THEN
    ! um2lm case
    DEALLOCATE (hsurf_gl, p0_gl                             , STAT=istat)
  ENDIF

  DEALLOCATE (t_s_lm    , t_m_lm    , t_snow_lm , t_cl_lm   ,             &
              w_cl_lm   , w_i_lm    , w_snow_lm , w_g1_lm   , w_g2_lm   , &
              qv_s_lm   ,                                     STAT=istat)

  IF (ALLOCATED(w_g3_lm)) THEN
    DEALLOCATE (w_g3_lm ,                                     STAT=istat)
  ENDIF

  IF (ALLOCATED(lmask_lm)) THEN
    DEALLOCATE (lmask_lm,                                     STAT=istat)
  ENDIF

  IF (ALLOCATED(t_skin_gl)) THEN
    DEALLOCATE (t_skin_gl                                   , STAT=istat)
  ENDIF

  IF (lgme2lm) THEN

    ! GME fields from reading external parameters
    !--------------------------------------------
    DEALLOCATE (lolp_gme, fis_gme, soiltyp_gme,    STAT=istat)
  ENDIF

  IF (llm2lm .OR. lec2lm) THEN
    DEALLOCATE (soiltyp_in, vio3_in, hmo3_in, lolp_in,         STAT=istat)
    IF (ALLOCATED(hsurf_in)) THEN
      DEALLOCATE (hsurf_in, t_m_in,                            STAT=istat)
    ENDIF

    DEALLOCATE (ps_in, t_snow_in, t_s_in, t_2m_in, t_cl_in, w_snow_in,       &
                w_g1_in, w_g2_in, w_g3_in, w_cl_in, qv_s_in, qv_2m_in,       &
                dpsdt_in, fic_in, t_ke_in, grh_in, t_g_in, w_i_in,  STAT=istat)

    IF (ALLOCATED(t_g1_in)) THEN
      DEALLOCATE (t_g1_in, t_g2_in, t_g3_in,                   STAT=istat)
    ENDIF

    DEALLOCATE (u_in, v_in, t_in, qv_in, qc_in,                STAT=istat)

    IF (ALLOCATED(w_in)) THEN
      DEALLOCATE (w_in, pp_in, p_in,                           STAT=istat)
    ENDIF

    IF (ALLOCATED(qi_in)) THEN
      DEALLOCATE (qi_in,                                       STAT=istat)
    ENDIF
    IF (ALLOCATED(qr_in)) THEN
      DEALLOCATE (qr_in,                                       STAT=istat)
    ENDIF
    IF (ALLOCATED(qs_in)) THEN
      DEALLOCATE (qs_in,                                       STAT=istat)
    ENDIF
    IF (ALLOCATED(qg_in)) THEN
      DEALLOCATE (qg_in,                                       STAT=istat)
    ENDIF

    IF (lgme2lm) THEN
      DEALLOCATE (w_intpol, n_intpol, m_intpol, l_intpol)
    ENDIF
  ENDIF

! iso code
  IF (liso) THEN
    IF (ALLOCATED(riso_lm))  DEALLOCATE (riso_lm, STAT=istat)
    IF (ALLOCATED(risosoil_lm))  DEALLOCATE (risosoil_lm, STAT=istat)
! Hui Tang 2013-11-20
    IF (ALLOCATED(drisoke_gl))  DEALLOCATE (drisoke_gl, STAT=istat)
    IF (ALLOCATED(riso_in))  DEALLOCATE (riso_in, STAT=istat)
    IF (ALLOCATED(risosoil_in))  DEALLOCATE (risosoil_in, STAT=istat)
  ENDIF
! end iso code

!------------------------------------------------------------------------------
!- End SUBROUTINE free_memory
!------------------------------------------------------------------------------

END SUBROUTINE free_memory

!==============================================================================

END MODULE src_cleanup
