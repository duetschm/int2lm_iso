!+ Module for the allocation of memory
!==============================================================================

MODULE src_memory

!==============================================================================
!
! Description:
!   This module provides routines for the allocation of memory for the 
!   different models. At the end of the program, memory is deallocated again.
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
! 1.2        2005/07/22 Ulrich Schaettler
!  Added fields for_e, for_d
! 1.3        2005/12/12 Ulrich Schaettler
!  Added fields for prognostic treatment of snow density in the LM
! 1.4        2006/01/03 Ulrich Schaettler
!  Bug correction in the dimension of rho_snow_in
! V1_5         2007/07/09 Ulrich Schaettler
!  Eliminated akhlm, bkhlm
!  Replaced ke_soil to ke_soil_lm, ke_soil_in
!  Added fields qr_lm, qs_lm, qg_lm, qr_in, qs_in, qg_in
!  Added options and fields for llake, lcm2lm
!  Transition from volumetric to relative input soil moisture
!  Added additional chemistry fields
! V1_6         2007/09/07 Ulrich Schaettler
!  Bug correction for dimensions of cgas_in, caero_in
! V1_7         2007/11/26 Ulrich Schaettler
!  Added additional fields for vegetation and rest of plcov, lai for input model
! V1_8         2008/05/29 Ulrich Schaettler
!  Renamed all gz0-variables to z0 and only work with z0-values
!  Allocate additional field rootdp_mx
!  Allocate fields for subgrid scale orography and topographical corrections
!  Introduce debug output
! V1_9         2009/09/03 Ulrich Schaettler, et al.
!  Allocate memory for new external parameters
!  Removed option to interpolate ndvi-values from GME 
!   (now they are available in the COSMO external parameters)
! V1_10        2009/12/17 Oliver Fuhrer, Ulrich Schaettler
!  Added field p_in for reading standard pressure
!  Allocation of fields for lum2lm
!  Allocation of fields for lake_coldstart
!  Allocation of fields for seaice (JPS)
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler, Burkhardt Rockel
!  Allocate qv_2m_in, t_2m_in, qv_2m_gl, t_2m_gl (for JMA data ) and rlai_in
!  Added field for salt_lk_lm
! V1_17        2011/03/11 Ulrich Schaettler
!  Added field fr_urban_lm for urban fraction data (K. Trusilova)
! V1_19        2012/06/06 Ulrich Schaettler, Susanne Brienen, Daniel Luethi
!                         Burkhardt Rockel
!  Renamed ak_in_uv, bk_in_uv to ak_in_rho, bk_in_rho according to UM conventions
!  Allocate t_so_gl and w_so_gl with ke_soil_in (SB)
!  Allocate new fields for surface albedo
!  Allocate fields for hybrid height coordinate also for lcm_hgt_coor
! V1_20        2012/09/03 Anne Roches
!  Allocate lmask_lm in all cases, because it is used as argument to
!   subroutine calls (Anne Roches)
! V1_22        2013/07/11 Ulrich Schaettler, Davide Cesari, KIT
!  Renamed lprog_qrqs to lprog_qr_qs to be consistent with other names
!  Allocate pv_in big enough
!  Removed variables vcoord_in, vcoord, aklm, bklm, sigmalm_sl
!  Extend the bounds of dt_so_gl (Davide)
!  Removed the ART parts, which are now in a separate component
!
! Added fields for water isotope simulation (Stephan Pfahl)
! Added variable for height correction of soil isotopes. Hui Tang 2013-11-20
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters , ONLY :   &
  ireals,    & ! KIND-type parameters for real variables
  iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_int2lm_control,       ONLY:                              &
    nhori,        & ! number of sectors for the horizon array by the
    lgme2lm,      & ! if .TRUE., gme->lm,
    llm2lm,       & ! if .TRUE., lm ->lm
    lum2lm,       & ! if .TRUE., um ->lm
    lec2lm,       & ! if .TRUE., ec ->lm
    lcm2lm,       & ! if .TRUE., cm ->lm
    lprog_qi,     & ! if .TRUE., interpolate qi to LM grid
    lprog_qr_qs,  & ! if .TRUE., interpolate qr,qs to LM grid
    lprog_qg,     & ! if .TRUE., interpolate qg to LM grid
    lprog_rho_snow,&! if .TRUE., interpolate rho_snow from GME to LM grid
    itype_ndvi,   & ! to choose treatment of surface parameters (plcov, lai)
    lforest,      & ! if .TRUE., run with forest (evergreen and deciduous)
    lurban,       & ! if .TRUE., run the urban module
    lsso,         & ! process parameters for sso scheme
    lradtopo,     & ! process parameters for topographic correction of radiation
    lseaice,      & ! if .TRUE., run with sea ice model
    llake,        & ! if .TRUE., run with lake  !_br
    llake_coldstart,& ! if .TRUE., initialize prognostic lake variables for cold start
    itype_albedo, & ! choose treatment of solar surface albedo
    lemiss,       & ! if .TRUE., run with external parameter for surface emissivity
    lstomata,     & ! if .TRUE., run with external parameter for stomata resistance
    lbd_frame,    & ! if .TRUE., boundary fields include only frames
    npstrframe,   & ! thickness of output frames
    lmulti_layer_lm,&! if .TRUE., compute soil fields for multi-layer soil
                     ! model in the outgoing data
    lmulti_layer_in,&! if .TRUE., incoming data have soil fields from the
                     ! multi-layer soil model
    nl_soil_in,   & ! number of soil layers in input model
    nl_soil_lm,   & ! number of soil layers in LM, resp. HM
    itype_aerosol,& ! to choose treatment of surface parameters (plcov, lai)
! SP, 201405
!    itype_t_cl,   & ! to choose origin and treatment of deep soil temperature
    idbg_level,   & ! to control verbosity of output
    lprintdeb_all,& ! whether all PEs print debug output
! iso code
    liso            ! if .TRUE., include variables for water isotope simulation
! end iso code

!------------------------------------------------------------------------------

USE data_int2lm_io,     ONLY:                              &
    inrvert_in,   & ! number of vertical coordinate parameters of input data
    pv_in           ! array for vertical coordinate parameters of input data

!------------------------------------------------------------------------------

USE data_grid_lm,       ONLY:                              &
    ielm_tot,     & ! ie for LM, whole area
    jelm_tot,     & ! je for LM, whole area
    ie2lm,        & ! ie for LM, local domain
    je2lm,        & ! je for LM, local domain
    kelm,         & !
    kedim,        & !
    ke1lm,        & !
    ke_soil_lm      ! number of levels in multi-layer soil model in output

!------------------------------------------------------------------------------

USE data_grid_in,       ONLY:                              &
    ie_in,       & ! ie for input grid, local domain
    je_in,       & ! je for input grid, local domain
    ke_in,       & ! ke for input grid
    ke_in_tot,   & ! total number of levels for input grid (including nlevskip levels)
    ke1in,       & ! ke1 for input grid
    kedim_in,    & ! MAX (ke_in, ke_hybrid)
    ke_soil_in,  & ! number of levels in multi-layer soil model in input
    jd_min,      & ! smallest index for diamonds for a LM subdomain
    jd_max,      & ! biggest index for diamonds for a LM subdomain
    igg1sm2,     & ! = igg1s - 2
    igg1ep2,     & ! = igg1e + 2
    igg2sm2,     & ! = igg2s - 2
    igg2ep2,     & ! = igg2e + 2
    ak_in ,      & ! vertical coordinate parameters for half levels
    bk_in ,      & !                  - " -
    akh_in,      & ! vertical coordinate parameters for main levels
    bkh_in,      & !                  - " -
    akh_in_rho,  & ! coefficients for main levels for u, v (lum2lm)
    bkh_in_rho,  & ! coefficients for main levels for u, v (lum2lm)
    dak_in,      & ! difference of coordinate parameters
    dbk_in,      & !                  - " -
    lcm_hgt_coor   ! Input data has hybrid height coordinates 

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY:                              &
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    isubpos            ! positions of the subdomains in the total domain.

!------------------------------------------------------------------------------

USE data_fields_lm,     ONLY:                              &
  p0_lm        ,   & ! reference pressure                               ( Pa  )
  dp0_lm       ,   & ! reference pressure thickness of layers           ( Pa  )
  rho0_lm      ,   & ! reference density at the full model levels       (kg/m3)
  t0_lm        ,   & ! reference temperature                            (  K  )
  hhl_lm       ,   & ! height of half-levels of LM                      (  m  )
  fis_lm       ,   & ! orography * G                                    (m2/s2)
  ps_lm        ,   & ! surface pressure                                 ( Pa  )
  hsurf_lm     ,   & ! orography                                        (  m  )
  fr_land_lm   ,   & ! land fraction of grid element                    (  1  )
  fr_lake_lm   ,   & ! lake fraction of grid element                    (  1  )
  depth_lk_lm  ,   & ! lake depth                                       (  m  )
  salt_lk_lm   ,   & ! lake salinity                                    ( 1E03)
  z0_lm        ,   & ! roughness length                                 (  m  )
  soiltyp_lm   ,   & ! type of the soil (keys 0-9)                      (  1  )
  plcov_lm     ,   & ! fraction covered by plants                       (  1  )
  plcov_mx_lm  ,   & ! plant cover during vegetation time               (  1  )
  plcov_mn_lm  ,   & ! plant cover during time of rest                  (  1  )
  lai_mx_lm    ,   & ! leaf area index during vegetation time           (  1  )
  lai_mn_lm    ,   & ! leaf area index during time of rest              (  1  )
  lai_lm       ,   & ! leaf area index                                  (  1  )
  rootdp_lm    ,   & ! depth of the roots                               (  m  )
  rootdp_mx    ,   & ! depth of the roots from external parameters      (  m  )
  for_e_lm     ,   & ! ground fraction covered by evergreen forest      (  -  )
  for_d_lm     ,   & ! ground fraction covered by deciduous forest      (  -  )
  fr_urban_lm  ,   & ! urban fraction of grid element                   (  -  )
  sso_stdh_lm  ,   & ! standard deviation of subgrid scale orography    (  m  )
  sso_gamma_lm ,   & ! anisotropy of the orography                      (  -  )
  sso_theta_lm ,   & ! angle betw. principal axis of orography and E    ( rad )
  sso_sigma_lm ,   & ! mean slope of subgrid scale orography            (  -  )
  skyview_lm   ,   & ! sky view                                         (  1  )
  slo_asp_lm   ,   & ! slope aspect                                     ( rad )
  slo_ang_lm   ,   & ! slope angle                                      ( rad )
  horizon_lm         ! horizon                                          ( rad )

USE data_fields_lm, ONLY : &
  alb_dry_lm   ,   & ! surface albedo field for dry soil                (  1  )
  alb_sat_lm   ,   & ! surface albedo field for saturated soil          (  1  )
  alb_dif_lm   ,   & ! solar surface albedo - diffuse                   (  1  )
  alb_dif12_lm ,   & ! solar surface albedo - diffuse                   (  1  )
  emis_rad_lm  ,   & ! thermal radiative surface emissivity             (  1  )
  prs_min_lm   ,   & ! minimum stomata resistance of plants             ( s/m )
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
  lai12              ! monthly climatology of leaf area index           (  1  )

USE data_fields_lm, ONLY : &
  t_cl_lm      ,   & ! temperature between medium and lower soil layer  (  K  )
  t_s_lm       ,   & ! temperature of the ground surface                (  K  )
  t_snow_lm    ,   & ! temperature of the snow surface                  (  K  )
  t_m_lm       ,   & ! temperature between upper and medium soil layer  (  K  )
  qv_s_lm      ,   & ! specific water vapor content on the surface      (kg/kg)
  t_2m_gl      ,   & ! 2m temperature                                   (  K  )
  qv_2m_gl     ,   & ! 2m humidity                                      (kg/kg)
  w_snow_lm    ,   & ! water content of the snow                        (m H2O)
  w_i_lm       ,   & ! water content of the interception storage        (m H2O)
  w_g1_lm      ,   & ! water content of the upper soil layer            (m H2O)
  w_g2_lm      ,   & ! water content of the medium soil layer           (m H2O)
  w_g3_lm      ,   & ! water content of the lower soil layer            (m H2O)
                     ! (if nl_soil_lm = 3, unused otherwise)
  w_cl_lm      ,   & ! climatological deep soil water content           (m H2O)
  t_so_lm      ,   & ! multi-layer soil temperature                     (  K  )
  dt_so_gl     ,   & ! multi-layer soil temperature (for interpolation)
  w_so_lm      ,   & ! multi-layer soil moisture                        (m H2O)
  w_so_gl      ,   & ! multi-layer soil moisture (for interpolation)    (m H2O)
  freshsnw_lm  ,   & ! weighting function indicating 'freshness' of snow
  rho_snow_lm  ,   & ! for prognostic treatment of snow density         (kg/m3)
  hmo3_lm      ,   & ! height of maximum ozone concentration            ( Pa  )
  vio3_lm      ,   & ! total vertically integrated ozone content        (Pa O3)
  lolp_lm      ,   & ! Land Sea Mask of LM for 'M'atch Interpolation
  lmask_lm     ,   & ! mask of points on the frame
  latlm_m      ,   & ! latitudes of the LM grid points
  lonlm_m      ,   & ! longitudes of the LM grid points
  fis_gl       ,   & ! GME interpolated orography * G                   (m2/s2)
  ps_gl        ,   & ! surface pressure on the interpol. GME orogr.     ( Pa  )
  p0_gl        ,   & ! ref. pres. on full levels + interpol. COARSE LM oro.(Pa)
  dp0_gl       ,   & ! reference pressure thickness of layers           ( Pa  )
  rho0_gl      ,   & ! reference density at the full model levels       (kg/m3)
  hhl_gl       ,   & ! height of half-levels on the interpol. COARSE LM oro.(m)
  hsurf_gl     ,   & ! height of orography interpolated from coarse grid(  m  )
  hsurfs_gl    ,   & ! interpolated splitted parts of coarse topo       (  m  )
  fic_gl             ! check level of geopotential                      (m2/s2)

USE data_fields_lm, ONLY : &
  dpsdt_gl     ,   & ! surface pressure tendency                        (Pa/s )
  t_s_gl       ,   & ! temperature of the ground surface                (  K  )
  t_skin_gl    ,   & ! skin temperature of the ground surface           (  K  )  !_br
  rh_s_gl      ,   & ! relative humidity at the surface                 (kg/kg)
  dtms_gl      ,   & ! t_m_lm    - t_s_lm                               (  K  )
  dtkes_gl     ,   & ! t(ke)_lm  - t_s_lm                               (  K  )
  dtssnow_gl   ,   & ! t_s_lm    - t_snow_lm                            (  K  )
  u_lm         ,   & ! zonal wind speed                                 ( m/s )
  v_lm         ,   & ! meridional wind speed                            ( m/s )
  w_lm         ,   & ! vertical   wind speed                            ( m/s )
  t_lm         ,   & ! temperature                                      (  K  )
  p_lm         ,   & ! full pressure (needed for lum2lm)                ( Pa  )
  pp_lm        ,   & ! deviation from the reference pressure            ( Pa  )
  qv_lm        ,   & ! specific water vapor content                     (kg/kg)
  qc_lm        ,   & ! specific cloud water content                     (kg/kg)
  qi_lm        ,   & ! cloud ice content                                (kg/kg)
  qr_lm        ,   & ! rain      content                                (kg/kg)
  qs_lm        ,   & ! snow      content                                (kg/kg)
  qg_lm        ,   & ! graupel   content                                (kg/kg)
  w_intpol     ,   & ! interpolation weights for gme2lm                 (  -  )
  n_intpol     ,   & ! nearest GME gridpoint for nearest neighbor interpolation
  m_intpol     ,   & ! nearest GME gridpoint with same lsm for match interpolation
  l_intpol     ,   & ! to use a far away GME gridpoint with same lsm
  cgas_lm      ,   & !
  caero_lm     ,   & !
! iso code
  riso_lm      ,   & ! isotope ratios in water vapor
  risosoil_lm  ,   & ! isotope ratios in soil water
! Hui Tang 2013-11-20
  drisoke_gl         ! riso_lm(:,:,ke,1-2) - risosoil_lm (:,:,1-2) (1: 18O; 2: 2H);
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

USE data_fields_in,  ONLY: &
 fis_gme     ,          & ! orography * g                               (m2/s2)
 soiltyp_gme ,          & ! type of the soil (keys 0-9)                   --
 fr_land_gme ,          & ! land fraction of grid element               (  1  )
 z0_gme      ,          & ! surface roughness                           (  m  )
 plcmx_gme   ,          & ! vegetation: fraction covered by plants      (  1  )
 plcmn_gme   ,          & ! rest:       fraction covered by plants      (  1  )
 rlaimx_gme  ,          & ! vegetation: leaf area index                 (  1  )
 rlaimn_gme  ,          & ! rest:       leaf area index                 (  1  )
 root_gme    ,          & ! depth of the roots                          (  m  )
 lolp_gme    ,          & ! Land Sea Mask of GME for 'M'atch Interpolation
 fis_in      ,          & ! orography * g                               (m2/s2)
 hsurf_in    ,          & ! orography                                   (  m  )
 soiltyp_in  ,          & ! type of the soil (keys 0-9)                   --
 ps_in       ,          & ! surface pressure                            ( Pa  )
 sst_in      ,          & ! sea surface temperature                     (  K  )  !_br
 t_s_in      ,          & ! temperature of the ground surface           (  K  )
 t_2m_in     ,          & ! 2m temperature                              (  K  )
 t_skin_in   ,          & ! skin temperature of the ground surface      (  K  )  !_br
 t_snow_in   ,          & ! temperature of the snow-surface             (  K  )
 t_g1_in     ,          & ! temperature of first soil layer             (  K  )
 t_g2_in     ,          & ! temperature of second soil layer            (  K  )
 t_g3_in     ,          & ! temperature of third soil layer             (  K  )
 qv_s_in     ,          & ! specific water vapor content on the surface (kg/kg)
 qv_2m_in    ,          & ! specific water vapor content in 2m          (kg/kg)
 w_g1_in     ,          & ! water content of the upper soil layer       (m H2O)
 w_g2_in     ,          & ! water content of the medium soil layer      (m H2O)
 w_g3_in     ,          & ! water content of the deepest soil layer     (m H2O)
 t_so_in     ,          & ! temperature for new multi layer soil model  (  K  )
 w_so_in     ,          & ! soil moisture for multi layer soil model    (m H2O)
 w_so_rel_in ,          & ! multi-layer relative   soil moisture        (  1  )
 freshsnw_in ,          & ! weighting function indicating 'freshness' of snow
 rho_snow_in ,          & ! for prognostic treatment of snow density    (kg/m3)
 t_ice_in    ,          & ! temperature of sea ice surface              (  K  )
 h_ice_in    ,          & ! sea ice thickness                           (  m  )
 t_ke_in     ,          & ! temperature lowest layer                    (  K  )
 grh_in      ,          & ! generalized relative humidity at one level  (  %  )
 lolp_in                  ! Land Sea Mask of input for 'M'atch Interpolation

USE data_fields_in,  ONLY: &
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
 z0_in       ,          & ! surface roughness                           (  m  )
 fr_land_in  ,          & ! land fraction of grid element               (  1  )
 root_in     ,          & ! depth of the roots                          (  m  )
 plcov_in    ,          & ! actual fraction covered by plants           (  1  )
 plcmx_in    ,          & ! vegetation: fraction covered by plants      (  1  )
 plcmn_in    ,          & ! rest:       fraction covered by plants      (  1  )
 rlai_in     ,          & ! vegetation: leaf area index                 (  1  )
 rlaimx_in   ,          & ! vegetation: leaf area index                 (  1  )
 rlaimn_in   ,          & ! rest:       leaf area index                 (  1  )
 vio3_in     ,          & ! total vertically integrated ozone content   (Pa O3)
 hmo3_in     ,          & ! height of maximum ozone concentration       ( Pa  )
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

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Public and Private Subroutines

PUBLIC   alloc_lm, alloc_gme, alloc_coarse_grid

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure to allocate memory for the fine LM grid
!------------------------------------------------------------------------------

SUBROUTINE alloc_lm (yerror, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the allocation of memory for the LM output
!   fields on the fine grid.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist
CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

! Local variables
INTEGER (KIND=iintegers) :: &
  izerrstat, istat, l1, l2, izdebug, k

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine alloc_lm
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierror    = 0_iintegers
  yerror    = '         '
  istat     = 0_iintegers
  izerrstat = 0_iintegers

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '    Starting Memory Allocation'
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Fields for the reference atmosphere and I/O
!------------------------------------------------------------------------------

  inrvert_in = 2*ke_in_tot + 20   ! is set correctly later by input data
  ALLOCATE (pv_in(inrvert_in)            , STAT=istat);  pv_in  = 0.0_ireals
                                               izerrstat = izerrstat + istat

  ALLOCATE (hhl_lm    (ie2lm,je2lm,ke1lm), STAT=istat);  hhl_lm = 0.0_ireals
                                               izerrstat = izerrstat + istat
  ALLOCATE (p0_lm     (ie2lm,je2lm,kelm),  STAT=istat);  p0_lm  = 0.0_ireals
                                               izerrstat = izerrstat + istat
  ALLOCATE (dp0_lm    (ie2lm,je2lm,kelm),  STAT=istat);  dp0_lm = 0.0_ireals
                                               izerrstat = izerrstat + istat
  ALLOCATE (rho0_lm   (ie2lm,je2lm,kelm),  STAT=istat); rho0_lm = 0.0_ireals
                                               izerrstat = izerrstat + istat
  ALLOCATE (t0_lm     (ie2lm,je2lm,kelm),  STAT=istat);  t0_lm  = 0.0_ireals
                                               izerrstat = izerrstat + istat

!------------------------------------------------------------------------------
! Section 3: Fields for external parameters
!------------------------------------------------------------------------------

  ! Fields for external parameters
  ALLOCATE (fis_lm     (ie2lm,je2lm), STAT=istat);  fis_lm     = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (hsurf_lm   (ie2lm,je2lm), STAT=istat);  hsurf_lm   = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (fr_land_lm (ie2lm,je2lm), STAT=istat);  fr_land_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
  IF (llake) THEN
    ALLOCATE (fr_lake_lm (ie2lm,je2lm), STAT=istat);  fr_lake_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (depth_lk_lm(ie2lm,je2lm), STAT=istat);  depth_lk_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (salt_lk_lm(ie2lm,je2lm), STAT=istat);  salt_lk_lm = 0.0_ireals  
                                              izerrstat = izerrstat + istat
    IF (llake_coldstart) THEN
      ALLOCATE (t_mnw_lk_lm(ie2lm,je2lm), STAT=istat);  t_mnw_lk_lm= 0.0_ireals
                                              izerrstat = izerrstat + istat
      ALLOCATE (t_wml_lk_lm(ie2lm,je2lm), STAT=istat);  t_wml_lk_lm= 0.0_ireals
                                              izerrstat = izerrstat + istat
      ALLOCATE (t_bot_lk_lm(ie2lm,je2lm), STAT=istat);  t_bot_lk_lm= 0.0_ireals
                                              izerrstat = izerrstat + istat
      ALLOCATE (t_b1_lk_lm(ie2lm,je2lm),  STAT=istat);  t_b1_lk_lm= 0.0_ireals
                                              izerrstat = izerrstat + istat
      ALLOCATE (c_t_lk_lm  (ie2lm,je2lm), STAT=istat);  c_t_lk_lm  = 0.0_ireals
                                              izerrstat = izerrstat + istat
      ALLOCATE (h_ml_lk_lm (ie2lm,je2lm), STAT=istat);  h_ml_lk_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
      ALLOCATE (h_b1_lk_lm (ie2lm,je2lm), STAT=istat);  h_b1_lk_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ENDIF
  ENDIF
  ALLOCATE (z0_lm      (ie2lm,je2lm), STAT=istat);  z0_lm      = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (soiltyp_lm (ie2lm,je2lm), STAT=istat);  soiltyp_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (plcov_lm   (ie2lm,je2lm), STAT=istat);  plcov_lm   = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (plcov_mx_lm(ie2lm,je2lm), STAT=istat);  plcov_mx_lm= 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (plcov_mn_lm(ie2lm,je2lm), STAT=istat);  plcov_mn_lm= 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (lai_mx_lm  (ie2lm,je2lm), STAT=istat);  lai_mx_lm  = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (lai_mn_lm  (ie2lm,je2lm), STAT=istat);  lai_mn_lm  = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (lai_lm     (ie2lm,je2lm), STAT=istat);  lai_lm     = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (rootdp_lm  (ie2lm,je2lm), STAT=istat);  rootdp_lm  = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (rootdp_mx  (ie2lm,je2lm), STAT=istat);  rootdp_mx  = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (lolp_lm    (ie2lm,je2lm), STAT=istat);  lolp_lm    = .FALSE.
                                              izerrstat = izerrstat + istat
  ALLOCATE (vio3_lm    (ie2lm,je2lm), STAT=istat);  vio3_lm    = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ALLOCATE (hmo3_lm    (ie2lm,je2lm), STAT=istat);  hmo3_lm    = 0.0_ireals
                                              izerrstat = izerrstat + istat
  IF (lforest) THEN
    ALLOCATE (for_e_lm (ie2lm,je2lm), STAT=istat);  for_e_lm   = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (for_d_lm (ie2lm,je2lm), STAT=istat);  for_d_lm   = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ENDIF

  IF (lurban) THEN
    ALLOCATE (fr_urban_lm (ie2lm,je2lm), STAT=istat);  fr_urban_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ENDIF

  IF (lsso) THEN
    ALLOCATE (sso_stdh_lm (ie2lm,je2lm), STAT=istat); sso_stdh_lm  = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (sso_gamma_lm(ie2lm,je2lm), STAT=istat); sso_gamma_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (sso_theta_lm(ie2lm,je2lm), STAT=istat); sso_theta_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (sso_sigma_lm(ie2lm,je2lm), STAT=istat); sso_sigma_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ENDIF
  IF (lradtopo) THEN
    ALLOCATE (skyview_lm(ie2lm,je2lm), STAT=istat); skyview_lm  = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (slo_ang_lm(ie2lm,je2lm), STAT=istat); slo_ang_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (slo_asp_lm(ie2lm,je2lm), STAT=istat); slo_asp_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (horizon_lm(ie2lm,je2lm,nhori),STAT=istat); horizon_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ENDIF
  IF (lemiss) THEN
    ALLOCATE (emis_rad_lm(ie2lm,je2lm),STAT=istat); emis_rad_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ENDIF

  IF     (itype_albedo == 2) THEN
    ALLOCATE (alb_dry_lm(ie2lm,je2lm),STAT=istat); alb_dry_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (alb_sat_lm(ie2lm,je2lm),STAT=istat); alb_sat_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ELSEIF (itype_albedo == 3) THEN
    ALLOCATE (alb_dif_lm(ie2lm,je2lm),STAT=istat); alb_dif_lm = 0.0_ireals
    ALLOCATE (alb_dif12_lm(ie2lm,je2lm,12), STAT=istat);  alb_dif12_lm = 0.0_ireals
                                       izerrstat = izerrstat + istat
  ENDIF

  IF (lstomata) THEN
    ALLOCATE (prs_min_lm(ie2lm,je2lm), STAT=istat); prs_min_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ENDIF

  IF (itype_aerosol == 2) THEN
    ALLOCATE (aer_su12_lm(ie2lm,je2lm,12), STAT=istat);  aer_su12_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (aer_du12_lm(ie2lm,je2lm,12), STAT=istat);  aer_du12_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (aer_or12_lm(ie2lm,je2lm,12), STAT=istat);  aer_or12_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (aer_bc12_lm(ie2lm,je2lm,12), STAT=istat);  aer_bc12_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (aer_ss12_lm(ie2lm,je2lm,12), STAT=istat);  aer_ss12_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat

    ALLOCATE (aer_su_lm(ie2lm,je2lm), STAT=istat);  aer_su_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (aer_du_lm(ie2lm,je2lm), STAT=istat);  aer_du_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (aer_or_lm(ie2lm,je2lm), STAT=istat);  aer_or_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (aer_bc_lm(ie2lm,je2lm), STAT=istat);  aer_bc_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (aer_ss_lm(ie2lm,je2lm), STAT=istat);  aer_ss_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ENDIF
  IF     (itype_ndvi == 1) THEN
    ALLOCATE (ndvi_mrat_lm(ie2lm,je2lm,12), STAT=istat);  ndvi_mrat_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (ndviratio_lm(ie2lm,je2lm),    STAT=istat);  ndviratio_lm = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ELSEIF (itype_ndvi == 2) THEN
    ALLOCATE (plcov12(ie2lm,je2lm,12), STAT=istat); plcov12 = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (lai12  (ie2lm,je2lm,12), STAT=istat); lai12   = 0.0_ireals
                                              izerrstat = izerrstat + istat
    ALLOCATE (z012   (ie2lm,je2lm,12), STAT=istat); z012    = 0.0_ireals
                                              izerrstat = izerrstat + istat
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '    Fields for reference atmosphere and external parameters allocated'
  ENDIF

!------------------------------------------------------------------------------
! Section 4: 3D atmospheric fields
!------------------------------------------------------------------------------

  ALLOCATE (u_lm (ie2lm,je2lm,kedim),  STAT=istat);   u_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (v_lm (ie2lm,je2lm,kedim),  STAT=istat);   v_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (w_lm (ie2lm,je2lm,kedim+1),STAT=istat);   w_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (t_lm (ie2lm,je2lm,kedim), STAT=istat);    t_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (pp_lm(ie2lm,je2lm,kedim),  STAT=istat);  pp_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (qv_lm(ie2lm,je2lm,kedim), STAT=istat);   qv_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (qc_lm(ie2lm,je2lm,kedim), STAT=istat);   qc_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  IF (lprog_qi) THEN
    ALLOCATE (qi_lm(ie2lm,je2lm,kedim),STAT=istat);  qi_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ENDIF
  IF (lprog_qr_qs) THEN
    ALLOCATE (qr_lm(ie2lm,je2lm,kedim),STAT=istat);  qr_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
    ALLOCATE (qs_lm(ie2lm,je2lm,kedim),STAT=istat);  qs_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ENDIF
  IF (lprog_qg) THEN
    ALLOCATE (qg_lm(ie2lm,je2lm,kedim),STAT=istat);  qg_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ENDIF

  IF (lum2lm .OR. lcm_hgt_coor) THEN
    ALLOCATE (p_lm (ie2lm,je2lm,kedim),STAT=istat);   p_lm(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '    3D atmospheric fields allocated with dims:  ', ie2lm, je2lm, kedim
  ENDIF

!------------------------------------------------------------------------------
! Section 5: intermediate fields for interpolation results
!------------------------------------------------------------------------------

  ALLOCATE (fis_gl    (ie2lm,je2lm), STAT=istat); fis_gl    (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (ps_gl     (ie2lm,je2lm), STAT=istat); ps_gl     (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (fic_gl    (ie2lm,je2lm), STAT=istat); fic_gl    (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (dpsdt_gl  (ie2lm,je2lm), STAT=istat); dpsdt_gl  (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (t_s_gl    (ie2lm,je2lm), STAT=istat); t_s_gl    (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (t_2m_gl   (ie2lm,je2lm), STAT=istat); t_2m_gl   (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (qv_2m_gl  (ie2lm,je2lm), STAT=istat); qv_2m_gl  (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (rh_s_gl   (ie2lm,je2lm), STAT=istat); rh_s_gl   (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (dtssnow_gl(ie2lm,je2lm), STAT=istat); dtssnow_gl(:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (dtkes_gl  (ie2lm,je2lm), STAT=istat); dtkes_gl  (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (dtms_gl   (ie2lm,je2lm), STAT=istat); dtms_gl   (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  IF (llm2lm) THEN
    ALLOCATE (p0_gl(ie2lm,je2lm,ke_in), STAT=istat); p0_gl(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
    ALLOCATE(dp0_gl(ie2lm,je2lm,ke_in), STAT=istat);dp0_gl(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
    ALLOCATE(rho0_gl(ie2lm,je2lm,ke_in),STAT=istat);rho0_gl(:,:,:)= 0.0_ireals
                                                 izerrstat = izerrstat + istat
    ALLOCATE(hhl_gl(ie2lm,je2lm,ke1in), STAT=istat);hhl_gl(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
    ALLOCATE(hsurfs_gl(ie2lm,je2lm,2), STAT=istat); hsurfs_gl  (:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ENDIF

  IF (lum2lm .OR. lcm_hgt_coor) THEN
    ALLOCATE (p0_gl(ie2lm,je2lm,ke_in), STAT=istat); p0_gl(:,:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
    ALLOCATE(hsurf_gl(ie2lm,je2lm),  STAT=istat); hsurf_gl  (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ENDIF

  IF (lcm2lm) THEN
    ALLOCATE(t_skin_gl (ie2lm,je2lm), STAT=istat);t_skin_gl (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '    Fields for interpolation results allocated'
  ENDIF

!------------------------------------------------------------------------------
! Section 6: 2D LM fields
!------------------------------------------------------------------------------

  ALLOCATE (t_s_lm    (ie2lm,je2lm), STAT=istat); t_s_lm    (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (t_snow_lm (ie2lm,je2lm), STAT=istat); t_snow_lm (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (w_i_lm    (ie2lm,je2lm), STAT=istat); w_i_lm    (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (w_snow_lm (ie2lm,je2lm), STAT=istat); w_snow_lm (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ALLOCATE (qv_s_lm   (ie2lm,je2lm), STAT=istat); qv_s_lm   (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  IF (lmulti_layer_lm) THEN
    ALLOCATE (t_so_lm (ie2lm,je2lm,0:ke_soil_lm+1), STAT=istat)
              t_so_lm  (:,:,:) = 0.0_ireals;     izerrstat = izerrstat + istat
    ALLOCATE (w_so_lm (ie2lm,je2lm,1:ke_soil_lm+1), STAT=istat)
              w_so_lm  (:,:,:) = 0.0_ireals;     izerrstat = izerrstat + istat
    ALLOCATE (freshsnw_lm(ie2lm,je2lm), STAT=istat)
              freshsnw_lm(:,:) = 0.0_ireals;     izerrstat = izerrstat + istat
    IF (lprog_rho_snow) THEN
      ALLOCATE (rho_snow_lm(ie2lm,je2lm), STAT=istat)
              rho_snow_lm(:,:) = 0.0_ireals;     izerrstat = izerrstat + istat
    ENDIF
  ENDIF

  IF (lseaice .OR. llake_coldstart) THEN
    ALLOCATE (t_ice_lm(ie2lm,je2lm), STAT=istat)
              t_ice_lm(:,:) = 0.0_ireals;     izerrstat = izerrstat + istat
    ALLOCATE (h_ice_lm(ie2lm,je2lm), STAT=istat)
              h_ice_lm(:,:) = 0.0_ireals;     izerrstat = izerrstat + istat
  ENDIF

  IF (lmulti_layer_in) THEN
    IF (lmulti_layer_lm) THEN ! useful mainly for ec2lm, it is safe in other cases
      k = MAX(ke_soil_in, ke_soil_lm)
    ELSE
      k = ke_soil_in
    ENDIF

    ALLOCATE (dt_so_gl(ie2lm,je2lm,1:k+1), STAT=istat)
              dt_so_gl (:,:,:) = 0.0_ireals;     izerrstat = izerrstat + istat
    ALLOCATE ( w_so_gl(ie2lm,je2lm,1:ke_soil_in+1), STAT=istat)
               w_so_gl (:,:,:) = 0.0_ireals;     izerrstat = izerrstat + istat
    IF (lcm2lm .OR. lgme2lm) THEN
      ! in case of lgme2lm, these fields are used for optional vertical interpolation
      ALLOCATE (t_cl_lm (ie2lm,je2lm), STAT=istat); t_cl_lm   (:,:) = 0.0_ireals
                                                   izerrstat = izerrstat + istat
      ALLOCATE (w_cl_lm (ie2lm,je2lm), STAT=istat); w_cl_lm   (:,:) = 0.0_ireals
                                                   izerrstat = izerrstat + istat
    ENDIF
  ELSE
    ALLOCATE (t_m_lm  (ie2lm,je2lm), STAT=istat); t_m_lm    (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
    ALLOCATE (t_cl_lm (ie2lm,je2lm), STAT=istat); t_cl_lm   (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
    ALLOCATE (w_g1_lm (ie2lm,je2lm), STAT=istat); w_g1_lm   (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
    ALLOCATE (w_g2_lm (ie2lm,je2lm), STAT=istat); w_g2_lm   (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
    IF (nl_soil_lm == 3) THEN
      ALLOCATE (w_g3_lm(ie2lm,je2lm),STAT=istat); w_g3_lm   (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
    ENDIF
    ALLOCATE (w_cl_lm (ie2lm,je2lm), STAT=istat); w_cl_lm   (:,:) = 0.0_ireals
                                                 izerrstat = izerrstat + istat
  ENDIF

  ALLOCATE (lmask_lm  (ie2lm,je2lm), STAT=istat); lmask_lm(:,:) = .TRUE.

  IF (lbd_frame) THEN
    ! Set to .FALSE. inside the frame
    DO l2 = MAX(1,1+npstrframe-(isubpos(my_cart_id,2)-2*nboundlines-1)), &
     MIN(je2lm,jelm_tot-npstrframe-(isubpos(my_cart_id,2)-2*nboundlines-1))
      DO l1 = MAX(1,1+npstrframe-(isubpos(my_cart_id,1)-2*nboundlines-1)), &
       MIN(ie2lm,ielm_tot-npstrframe-(isubpos(my_cart_id,1)-2*nboundlines-1))
        lmask_lm(l1,l2) = .FALSE.
      ENDDO
    ENDDO
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '    2D COSMO Model fields allocated'
  ENDIF

!------------------------------------------------------------------------------
! Section 7: And all the rest
!------------------------------------------------------------------------------

  IF (lgme2lm) THEN
    ALLOCATE (w_intpol  (ie2lm,je2lm,19), STAT=istat); 
              w_intpol(:,:,:) = 0.0_ireals;    izerrstat = izerrstat + istat
    ALLOCATE (n_intpol  (ie2lm,je2lm, 3), STAT=istat); 
              n_intpol(:,:,:) = 0_iintegers;   izerrstat = izerrstat + istat
    ALLOCATE (m_intpol  (ie2lm,je2lm, 3), STAT=istat); 
              m_intpol(:,:,:) = 0_iintegers;   izerrstat = izerrstat + istat
    ALLOCATE (l_intpol  (ie2lm,je2lm),    STAT=istat); 
              l_intpol(:,:) = .FALSE.;         izerrstat = izerrstat + istat

    IF (izdebug > 10) THEN
      PRINT *, '    Additional fields allocated'
    ENDIF
  ENDIF

! iso code
  IF (liso) THEN
!   Change the last dimension of riso_lm to 6. Hui Tang 2013-11-20
    ALLOCATE (riso_lm (ie2lm,je2lm,kedim,6), STAT=istat)
    riso_lm = 0.0_ireals
    izerrstat = izerrstat + istat
    ALLOCATE (risosoil_lm (ie2lm,je2lm,2), STAT=istat)
    risosoil_lm = 0.0_ireals
    izerrstat = izerrstat + istat
!   Hui Tang 2013-11-20
    ALLOCATE (drisoke_gl (ie2lm,je2lm,2), STAT=istat)
    drisoke_gl = 0.0_ireals
    izerrstat = izerrstat + istat
  ENDIF
! end iso code

!------------------------------------------------------------------------------
! Section 8: Error code 
!------------------------------------------------------------------------------

  IF (izerrstat /= 0) THEN
    ierror  = 1
    yerror  = 'Memory allocation failed for LM fields'
    RETURN
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '    Memory Allocation finished'
  ENDIF

!------------------------------------------------------------------------------

END SUBROUTINE alloc_lm

!==============================================================================
!+ Module procedure to allocate memory for the coarse input grid (LM or EC)
!------------------------------------------------------------------------------

SUBROUTINE alloc_coarse_grid (yerror, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the allocation of memory for the input
!   fields on the coarse grid.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist
CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

! Local variables
INTEGER (KIND=iintegers) :: &
  izerrstat, istat, izdebug

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine alloc_coarse_grid
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierror    = 0_iintegers
  yerror    = '         '
  istat     = 0_iintegers
  izerrstat = 0_iintegers

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '    Starting Memory Allocation for coarse grid'
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Fields for the reference atmosphere
!------------------------------------------------------------------------------

  ! Fields for the reference atmosphere
  IF (llm2lm) THEN
    ALLOCATE (hhl_in (ie_in,je_in,ke1in), STAT=istat); hhl_in   = 0.0
                                        izerrstat = izerrstat + istat
    ALLOCATE (p0_in  (ie_in,je_in,ke_in), STAT=istat); p0_in    = 0.0
                                        izerrstat = izerrstat + istat
    IF (izdebug > 10) THEN
      PRINT *, '    Fields for reference atmosphere in case of llm2lm allocated'
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Fields for external parameters
!------------------------------------------------------------------------------

  ! Fields for external parameters
  ALLOCATE (hsurf_in    (ie_in,je_in), STAT=istat); hsurf_in   = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (fis_in      (ie_in,je_in), STAT=istat); fis_in     = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (z0_in       (ie_in,je_in), STAT=istat); z0_in      = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (fr_land_in  (ie_in,je_in), STAT=istat); fr_land_in = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (soiltyp_in  (ie_in,je_in), STAT=istat); soiltyp_in = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (plcov_in    (ie_in,je_in), STAT=istat); plcov_in   = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (plcmx_in    (ie_in,je_in), STAT=istat); plcmx_in   = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (plcmn_in    (ie_in,je_in), STAT=istat); plcmn_in   = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (rlai_in     (ie_in,je_in), STAT=istat); rlai_in    = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (rlaimx_in   (ie_in,je_in), STAT=istat); rlaimx_in  = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (rlaimn_in   (ie_in,je_in), STAT=istat); rlaimn_in  = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (root_in     (ie_in,je_in), STAT=istat); root_in    = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (vio3_in     (ie_in,je_in), STAT=istat); vio3_in    = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (hmo3_in     (ie_in,je_in), STAT=istat); hmo3_in    = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (lolp_in     (ie_in,je_in), STAT=istat); lolp_in    = .FALSE.
                                       izerrstat = izerrstat + istat

  IF (izdebug > 10) THEN
    PRINT *, '    Fields for external parameters from the coarse grid allocated'
  ENDIF

!------------------------------------------------------------------------------
! Section 4: 2D LM fields
!------------------------------------------------------------------------------

  ALLOCATE (ps_in       (ie_in,je_in), STAT=istat); ps_in      = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (t_snow_in   (ie_in,je_in), STAT=istat); t_snow_in  = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (t_s_in      (ie_in,je_in), STAT=istat); t_s_in     = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (t_2m_in     (ie_in,je_in), STAT=istat); t_2m_in    = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (w_snow_in   (ie_in,je_in), STAT=istat); w_snow_in  = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (qv_s_in     (ie_in,je_in), STAT=istat); qv_s_in    = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (qv_2m_in    (ie_in,je_in), STAT=istat); qv_2m_in   = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (dpsdt_in    (ie_in,je_in), STAT=istat); dpsdt_in   = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (fic_in      (ie_in,je_in), STAT=istat); fic_in     = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (t_ke_in     (ie_in,je_in), STAT=istat); t_ke_in    = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (grh_in      (ie_in,je_in), STAT=istat); grh_in     = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (t_g_in      (ie_in,je_in), STAT=istat); t_g_in     = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (w_i_in      (ie_in,je_in), STAT=istat); w_i_in     = 0.0
                                       izerrstat = izerrstat + istat

  IF (llm2lm) THEN
    IF (lmulti_layer_in) THEN
      ALLOCATE (t_so_in   (ie_in,je_in,0:ke_soil_in+1), STAT=istat)
                t_so_in   = 0.0;       izerrstat = izerrstat + istat
      ALLOCATE (w_so_in   (ie_in,je_in,1:ke_soil_in+1), STAT=istat)
                w_so_in   = 0.0;       izerrstat = izerrstat + istat
      ALLOCATE (freshsnw_in(ie_in,je_in)           , STAT=istat)
                freshsnw_in = 0.0;     izerrstat = izerrstat + istat
      IF (lprog_rho_snow) THEN
        ALLOCATE (rho_snow_in(ie_in,je_in), STAT=istat)
                rho_snow_in = 0.0;     izerrstat = izerrstat + istat
      ENDIF
    ELSE
      ALLOCATE (t_m_in  (ie_in,je_in), STAT=istat); t_m_in     = 0.0
                                       izerrstat = izerrstat + istat
      ALLOCATE (t_cl_in (ie_in,je_in), STAT=istat); t_cl_in    = 0.0
                                       izerrstat = izerrstat + istat
      ALLOCATE (w_g1_in (ie_in,je_in), STAT=istat); w_g1_in    = 0.0
                                       izerrstat = izerrstat + istat
      ALLOCATE (w_g2_in (ie_in,je_in), STAT=istat); w_g2_in    = 0.0
                                       izerrstat = izerrstat + istat
      IF (nl_soil_in == 3) THEN
        ALLOCATE (w_g3_in(ie_in,je_in),STAT=istat); w_g3_in    = 0.0
                                       izerrstat = izerrstat + istat
      ENDIF
      ALLOCATE (w_cl_in (ie_in,je_in), STAT=istat); w_cl_in    = 0.0
                                       izerrstat = izerrstat + istat
    ENDIF

    IF (lseaice) THEN
      ALLOCATE (t_ice_in(ie_in,je_in), STAT=istat)
                t_ice_in = 0.0;        izerrstat = izerrstat + istat
      ALLOCATE (h_ice_in(ie_in,je_in), STAT=istat)
                h_ice_in = 0.0;        izerrstat = izerrstat + istat
    ENDIF
  ENDIF

  IF (lec2lm) THEN
    ALLOCATE (t_g1_in   (ie_in,je_in), STAT=istat); t_g1_in    = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (t_g2_in   (ie_in,je_in), STAT=istat); t_g2_in    = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (t_g3_in   (ie_in,je_in), STAT=istat); t_g3_in    = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (t_cl_in   (ie_in,je_in), STAT=istat); t_cl_in    = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (w_g1_in   (ie_in,je_in), STAT=istat); w_g1_in    = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (w_g2_in   (ie_in,je_in), STAT=istat); w_g2_in    = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (w_g3_in   (ie_in,je_in), STAT=istat); w_g3_in    = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (w_cl_in   (ie_in,je_in), STAT=istat); w_cl_in    = 0.0
                                       izerrstat = izerrstat + istat
  ENDIF

  IF (lcm2lm) THEN
    ALLOCATE (t_skin_in (ie_in,je_in), STAT=istat); t_skin_in = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (sst_in    (ie_in,je_in), STAT=istat); sst_in  = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (w_so_rel_in(ie_in,je_in,1:ke_soil_in+1), STAT=istat)
              w_so_rel_in   = 0.0;   izerrstat = izerrstat + istat
! SP, 201405
!    IF (itype_t_cl == 0) THEN
!      ALLOCATE (t_cl_in (ie_in,je_in), STAT=istat); t_cl_in = 0.0
!                                       izerrstat = izerrstat + istat
!    ENDIF
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '    2D Fields for the coarse grid allocated'
  ENDIF

!------------------------------------------------------------------------------
! Section 5: 3D atmospheric fields
!------------------------------------------------------------------------------

  ALLOCATE (u_in   (ie_in,je_in,kedim_in), STAT=istat);     u_in  = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (v_in   (ie_in,je_in,kedim_in), STAT=istat);     v_in  = 0.0
                                       izerrstat = izerrstat + istat
  IF (llm2lm .OR. lum2lm .OR. lcm_hgt_coor) THEN
    ALLOCATE (w_in (ie_in,je_in,kedim_in+1),STAT=istat);    w_in  = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (p_in (ie_in,je_in,kedim_in), STAT=istat);      p_in = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE (pp_in(ie_in,je_in,kedim_in), STAT=istat);     pp_in = 0.0
                                       izerrstat = izerrstat + istat
  ENDIF
  ALLOCATE (t_in   (ie_in,je_in,kedim_in), STAT=istat);     t_in  = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (qv_in  (ie_in,je_in,kedim_in), STAT=istat);     qv_in = 0.0
                                       izerrstat = izerrstat + istat
  ALLOCATE (qc_in  (ie_in,je_in,kedim_in), STAT=istat);     qc_in = 0.0
                                       izerrstat = izerrstat + istat
  IF (lprog_qi) THEN
    ALLOCATE(qi_in (ie_in,je_in,kedim_in), STAT=istat);     qi_in = 0.0
                                       izerrstat = izerrstat + istat
  ENDIF
  IF (lprog_qr_qs) THEN
    ALLOCATE(qr_in (ie_in,je_in,kedim_in), STAT=istat);     qr_in = 0.0
                                       izerrstat = izerrstat + istat
    ALLOCATE(qs_in (ie_in,je_in,kedim_in), STAT=istat);     qs_in = 0.0
                                       izerrstat = izerrstat + istat
  ENDIF
  IF (lprog_qg) THEN
    ALLOCATE(qg_in (ie_in,je_in,kedim_in), STAT=istat);     qg_in = 0.0
                                       izerrstat = izerrstat + istat
  ENDIF

! iso code
  IF (liso) THEN
! increase the last dimension of riso_in to 6. Hui Tang 2013-11-20
    ALLOCATE (riso_in (ie_in,je_in,ke_in,6), STAT=istat)
    riso_in = 0.0_ireals
    izerrstat = izerrstat + istat
    ALLOCATE (risosoil_in (ie_in,je_in,2), STAT=istat)
    risosoil_in = 0.0_ireals
    izerrstat = izerrstat + istat
  ENDIF
! end iso code

  IF (izdebug > 10) THEN
    PRINT *, '    3D Fields for the coarse grid allocated'
  ENDIF

!------------------------------------------------------------------------------
! Section 6: fields for vertical coordinate parameters
!------------------------------------------------------------------------------

  ALLOCATE (ak_in (kedim_in+1),   bk_in (kedim_in+1),                     &
            akh_in(kedim_in  ),   bkh_in(kedim_in  ),                     &
            dak_in(kedim_in  ),   dbk_in(kedim_in  ),                     &
            STAT = istat)

  IF (lum2lm .OR. lcm_hgt_coor) THEN
    ALLOCATE (akh_in_rho (ke_in), bkh_in_rho (ke_in), STAT = istat)
  ENDIF

  izerrstat = izerrstat + istat

  ak_in       = 0.0
  bk_in       = 0.0
  akh_in      = 0.0
  bkh_in      = 0.0
  dak_in      = 0.0
  dbk_in      = 0.0

  IF (izdebug > 10) THEN
    PRINT *, '    Vertical coordinate parameters allocated'
  ENDIF

!------------------------------------------------------------------------------
! Section 7: Error code 
!------------------------------------------------------------------------------

  IF (izerrstat /= 0) THEN
    ierror  = 2
    yerror  = 'Memory allocation failed for coarse grid input fields'
    RETURN
  ELSE
    IF (izdebug > 10) THEN
      PRINT *, '    Memory allocation for coarse grid finished successfully'
    ENDIF
  ENDIF

!------------------------------------------------------------------------------

END SUBROUTINE alloc_coarse_grid

!==============================================================================
!==============================================================================
!+ Module procedure to allocate memory for the coarse input grid (LM or EC)
!------------------------------------------------------------------------------

SUBROUTINE alloc_gme (yerror, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the allocation of memory for the GME input
!   fields. Because of the memory amount necessary for the GME fields,
!   only few fields are allocated in long term memory. Most fields are only
!   used in the horizontal interpolation.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist
CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

! Local variables
INTEGER (KIND=iintegers) :: &
  izerrstat, istat

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine alloc_gme
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierror    = 0_iintegers
  yerror    = '         '
  istat     = 0_iintegers
  izerrstat = 0_iintegers

!------------------------------------------------------------------------------
! Section 2: Some external fields
!------------------------------------------------------------------------------

  ALLOCATE (lolp_gme    (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),  &
                                                             STAT=istat)
                            lolp_gme = .FALSE.;  izerrstat = izerrstat + istat
  ALLOCATE (fr_land_gme (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),  &
                                                            STAT=istat)
                       fr_land_gme = 0.0_ireals; izerrstat = izerrstat + istat
  ALLOCATE (z0_gme      (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),  &
                                                              STAT=istat)
                           z0_gme  = 0.0_ireals; izerrstat = izerrstat + istat
  ALLOCATE (fis_gme (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),      &
                                                        STAT=istat)
                       fis_gme = 0.0_ireals;     izerrstat = izerrstat + istat
  ALLOCATE (soiltyp_gme (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),  &
                                                            STAT=istat)
                       soiltyp_gme = 0.0_ireals; izerrstat = izerrstat + istat
  ALLOCATE (root_gme    (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),  &
                                                              STAT=istat)
                          root_gme = 0.0_ireals; izerrstat = izerrstat + istat
  ALLOCATE (plcmx_gme   (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),  &
                                                              STAT=istat)
                         plcmx_gme = 0.0_ireals; izerrstat = izerrstat + istat
  ALLOCATE (plcmn_gme   (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),  &
                                                              STAT=istat)
                         plcmn_gme = 0.0_ireals; izerrstat = izerrstat + istat
  ALLOCATE (rlaimx_gme   (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),  &
                                                              STAT=istat)
                         rlaimx_gme = 0.0_ireals; izerrstat = izerrstat + istat
  ALLOCATE (rlaimn_gme   (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),  &
                                                              STAT=istat)
                         rlaimn_gme = 0.0_ireals; izerrstat = izerrstat + istat

!------------------------------------------------------------------------------
! Section 3: Error code 
!------------------------------------------------------------------------------

  IF (izerrstat /= 0) THEN
    ierror  = 3
    yerror  = 'Memory allocation failed for GME fields'
    RETURN
  ENDIF

!------------------------------------------------------------------------------

END SUBROUTINE alloc_gme

!==============================================================================

END MODULE src_memory
