!+ Module for the definition of the Grib tables
!==============================================================================

MODULE src_gribtabs

!==============================================================================
!
! Description:
!   This module provides the definitions of the Grib tables for all different
!   models used.
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
!  Added fields for_e and for_d to the LM table;
!  Read T_SO instead of T_S, if lmulti_layer_in=.TRUE.
! 1.3        2005/12/12 Ulrich Schaettler
!  Added field rho_snow for prognostic treatment of snow density in the LM
! V1_5         2007/07/09 Ulrich Schaettler
!  Reading W_I only for initial data
!  Replaced ke_soil to ke_soil_in, ke_soil_lm
!  Added new table entries for qr, qs, qg
!  Added new table entries for chemistry variables
!  Changed interpolation type for w_snow to positive definite
!  Added SR setup_vartab_cm
!  Transition from volumetric to relative input soil moisture
!     (vw_so_in to w_so_rel_in)
!  Added new external field for IFS, to accept ECMWF orography on model levels
!     (Davide Cesari)
! V1_6         2007/09/07 Ulrich Schaettler
!  Set 'O' for optional variable for qi at the end, depending on lprog_qi
! V1_7         2007/11/26 Ulrich Schaettler
!  Introduced MX- and MN-fields for LAI, PLCOV for input models
!  Introduced field for monthly NDVI ratio in GME grib table
! V1_8         2008/05/29 Ulrich Schaettler
!  Interpolation of qr,qs also for GME variables
!  Added new external parameters for COSMO-Model grid for lsso and lradtopo
!  Renamed all gz0-variables to z0 and only work with z0-values
!  Removed bias- and factor-values for Z0
!  Bug correction in setup_vartab_cm
! V1_9         2009/09/03 Ulrich Schaettler, et al.
!  Added fields for additional external parameters and output fields: 
!  ndvi ratios, aerosol values (12 monthly and actual fields)
!  surface emissivity and stomata resistance of plants
!  Added field for depth_lk_lm  (Burkhardt Rockel)
!  Added fields for 12-monthly values for vegetation (Daniel Luethi)
!  Changed interpolation code for QV,QC,QI in IFS2LM (from QFT to LFF)
!    (Anne Roches, MCH)
!  correct the unit for PP for NetCDF output; it must be Pa (HJP)
! V1_10        2009/12/17 Oliver Fuhrer, Ulrich Schaettler, Jan-Peter Schulz
!  Link P to the pointer for standard pressure
!  Implement grib table for Unified Model data
!  Introduced additional fields for seaice and prognostic FLake variables
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler, Burkhardt Rockel, Anne Roches
!  Unified variable table for all input models and dwdlib / grib_api
!  Added fields for salt_lk_lm (BR)
!  Added parameter 'i' for inland water lakes (BR)
!  Grib leveltype for HORIZON put to 110 (multidimensional) and removed 'l'
!    from 4 external parameters HORIZON, SLO_ANG, SLO_ASP, SKYVIEW  (AR)
! V1_17        2011/03/11 Ulrich Schaettler, Davide Cesari
!  Included vartab definition of QR, QS also for IFS data, if lec2lm
!  Changed level type of LNPS, FIS_SH to hybrid for IFS products
! V1_18        2011/03/11 Ulrich Schaettler
!  Included vartab definition of fr_urban
! V1_19        2012/06/06 Ulrich Schaettler, Burkhardt Rockel, Daniel Luethi
!  Some extensions for UM data
!  Extensions for HIRLAM data: atmospheric data are decoded as hybrid instead of
!    hybridLayer (as is done in IFS)
!  Added grib definitions for alb_sat, alb_dry, alb_dif (as in COSMO),
!                             alb_dif12 (with old iee=84, itabtyp=2)
!  Corrections/additions in table definitions for lcm2lm
!    - Removed Z0 from external parameters
!    - deleted settings for dattyp for t_so, w_so, lnps, t_g1, t_g2, t_g3, w_g3
!    - changed dattyp for t_skin from IBO to IB
!    - changed dattyp for w_i    from I   to I O
!    - added dattyp for freshsnw: I O
!    - changed interpolation code for P in case of lcm_hgt_coor to LFF
! V1_20        2012/09/03 Burkhardt Rockel
!  Correct standard name for W
!  Changes to read T_S for both climate and forecast version in case of netCDF input
! V1_21        2013/03/25 Ulrich Schaettler
!  Added string for level type in variable table var_lm for COSMO-Model
!  Modified GRIB settings for albedo variables (ALB_DIF12)
! V1_22        2013/07/11 Ulrich Schaettler, Davide Cesari, KIT
!  Renamed lprog_qrqs to lprog_qr_qs to be consistent with other names
!  Renamed PRS_MIN to RSMIN, which is the official shortName
!  Removed unused variables from the use-lists
!  For IFS data: do not add W_SO and T_SO to input lists, even if 
!    lmulti_layer_in is set to .TRUE. (Davide)
!  Removed definitions for ART variables, which are now in a separate component
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
    itype_w_so_rel,&! type of relative soil moisture input
    lbdclim,      & ! if .TRUE., special boundary data for climate mode
    lante_0006,   & ! if .TRUE., force to use ECMWF dataset before 27 June 2000
    lpost_0006,   & ! if .TRUE., force to use ECMWF dataset after 27 June 2000
    lgme2lm,      & ! if .TRUE., gme->lm
    lgfs2lm,      & ! if .TRUE., gfs->lm
    lgsm2lm,      & ! if .TRUE., gsm->lm
    lec2lm,       & ! if .TRUE., ec ->lm
    llm2lm,       & ! if .TRUE., lm ->lm
    lum2lm,       & ! if .TRUE., um ->lm
    lhir2lm,      & ! if .TRUE., hirlam ->lm
    lcm2lm,       & ! if .TRUE., cm ->lm   !_br
    lvertwind_ini,& ! if .TRUE., compute vertical wind for LM for initial data
    lvertwind_bd, & ! if .TRUE., compute vertical wind for LM for boundary data
    lprog_qi,     & ! if .TRUE., interpolate qi to LM grid
    lprog_qr_qs,  & ! if .TRUE., interpolate qr,qs to LM grid
    lprog_qg,     & ! if .TRUE., interpolate qg to LM grid
    lprog_rho_snow,&! if .TRUE., interpolate rho_snow from GME to LM grid
    lseaice,      & ! if .TRUE., run with sea ice model
    lmulti_layer_in, & ! if .TRUE., incoming data have soil fields from the
                    ! multi-layer soil model
    nl_soil_lm,   & ! number of soil layers in LM, resp. HM
    nl_soil_in,   & ! number of soil layers in GME
    luse_t_skin,  & ! if .TRUE., use ECMWF skin temperature for surface
    pcontrol_fi,  & ! pressure of control level for geopotential
    l_smi,        & ! if .TRUE., interpolate soil moisture with SMI
! SP, 201405
    itype_t_cl,   & ! to choose origin and treatment of deep soil temperature
! iso code
    liso,         & ! if .TRUE., include variables for water isotope simulation
! end iso code
    l_art,        & ! if .TRUE., interpolate additional ART fields
    l_art_nested

!------------------------------------------------------------------------------

USE data_grid_lm,       ONLY:                              &
    ke_soil_lm      ! number of levels in multi-layer soil model in output  !_br

!------------------------------------------------------------------------------

USE data_grid_in,       ONLY:                              &
    ke_in,        & ! ke for input grid
    ke1in,        & ! ke1 for input grid
    ke_soil_in,   & ! number of input levels in multi-layer soil model !_br
    lcm_hgt_coor      ! Input data has hybrid height coordinates

!------------------------------------------------------------------------------

USE data_int2lm_io,            ONLY:                              &
    nvar_lm,      & ! maximum number of variables in fine grid LM variable table
    nvar_in,      & ! maximum number of variables in coarse grid variable table
    nvar_lm_chem, & ! maximum number of variables with chemistry
    nvar_in_chem, & ! maximum number of variables with chemistry
    nvar_lm_norm, & ! maximum number of variables without chemistry
    nvar_in_norm, & ! maximum number of variables without chemistry
    ydate_ini,    & ! start of the forecast yyyymmddhh (year,month,day,hour)
    yin_form_read,& ! input format of boundary data
    ar_des_lm,    & ! structure for LM variable table
    ar_des_input, & ! structure for GME variable table
    var_lm,       & ! variable for fine grid LM variable table
    var_in          ! variable for input model variable table

!------------------------------------------------------------------------------

USE data_fields_lm,     ONLY:                              &
  hhl_lm       ,   & ! height of half-levels of LM                      (  m  )
  fis_lm       ,   & ! orography * G                                    (m2/s2)
  hsurf_lm     ,   & ! orography                                        (  m  )
  fr_land_lm   ,   & ! land fraction of grid element                    (  1  )
  fr_lake_lm   ,   & ! lake fraction of grid element                    (  1  )
  depth_lk_lm  ,   & ! lake depth                                       (  m  )
  salt_lk_lm   ,   & ! lake salinity                                    ( g/kg)
  z0_lm        ,   & ! roughness length                                 (  m  )
  soiltyp_lm   ,   & ! type of the soil (keys 0-9)                      (  1  )
  plcov_lm     ,   & ! fraction covered by plants                       (  1  )
  plcov_mx_lm  ,   & ! plant cover during vegetation time               (  1  )
  plcov_mn_lm  ,   & ! plant cover during time of rest                  (  1  )
  lai_mx_lm    ,   & ! leaf area index during vegetation time           (m2 m-2)
  lai_mn_lm    ,   & ! leaf area index during time of rest              (m2 m-2)
  lai_lm       ,   & ! leaf area index                                  (m2 m-2)
  rootdp_lm    ,   & ! depth of the roots                               (  m  )
  for_e_lm     ,   & ! ground fraction covered by evergreen forest     --
  for_d_lm     ,   & ! ground fraction covered by deciduous forest     --
  fr_urban_lm  ,   & ! ground fraction covered by urban land            (  -  )
  sso_stdh_lm  ,   & ! standard deviation of subgrid scale orography    (  m  )
  sso_gamma_lm ,   & ! anisotropy of the orography                      (  -  )
  sso_theta_lm ,   & ! angle betw. principal axis of orography and E    ( rad )
  sso_sigma_lm ,   & ! mean slope of subgrid scale orography            (  -  )
  skyview_lm   ,   & ! sky view                                         (  1  )
  slo_asp_lm   ,   & ! slope aspect                                     ( rad )
  slo_ang_lm   ,   & ! slope angle                                      ( rad )
  horizon_lm   ,   & ! horizon                                          ( rad )
  t_cl_lm      ,   & ! temperature between medium and lower soil layer  (  K  )
  t_s_lm       ,   & ! temperature of the ground surface                (  K  )
  t_snow_lm    ,   & ! temperature of the snow surface                  (  K  )
  t_m_lm       ,   & ! temperature between upper and medium soil layer  (  K  )
  qv_s_lm            ! specific water vapor content on the surface      (kg/kg)

USE data_fields_lm,     ONLY:                              &
  alb_dry_lm   ,   & ! surface albedo field for dry soil                (  1  )
  alb_sat_lm   ,   & ! surface albedo field for saturated soil          (  1  )
  alb_dif_lm   ,   & ! solar surface albedo - diffuse                   (  1  )
  alb_dif12_lm ,   & ! solar surface albedo - diffuse                   (  1  )
  emis_rad_lm  ,   & ! thermal radiative surface emissivity             (  1  )
  prs_min_lm   ,   & ! minimum stomata resistance of plants             ( s/m )
  ndvi_mrat_lm ,   & ! ratio of monthly mean normalized differential    (  1  )
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
  w_snow_lm    ,   & ! water content of the snow                        (m H2O)
  w_i_lm       ,   & ! water content of the interception storage        (m H2O)
  w_g1_lm      ,   & ! water content of the upper soil layer            (m H2O)
  w_g2_lm      ,   & ! water content of the medium soil layer           (m H2O)
  w_g3_lm      ,   & ! water content of the lower soil layer            (m H2O)
                     ! (if nl_soil_lm = 3, unused otherwise)
  w_cl_lm      ,   & ! climatological deep soil water content           (m H2O)
  t_so_lm      ,   & ! multi-layer soil temperature                     (  K  )
  w_so_lm      ,   & ! multi-layer soil moisture                        (m H2O)
  freshsnw_lm  ,   & ! weighting function indicating 'freshness' of snow
  rho_snow_lm        ! for prognostic treatment of snow density         (kg/m3)

USE data_fields_lm, ONLY : &
  hmo3_lm      ,   & ! height of maximum ozone concentration            ( Pa  )
  vio3_lm      ,   & ! total vertically integrated ozone content        (Pa O3)
  latlm_m      ,   & ! latitudes of the LM grid points
  lonlm_m      ,   & ! longitudes of the LM grid points
  ps_gl        ,   & ! surface pressure on the interpol. GME orogr.     ( Pa  )
  u_lm         ,   & ! zonal wind speed                                 ( m/s )
  v_lm         ,   & ! meridional wind speed                            ( m/s )
  w_lm         ,   & ! vertical   wind speed                            ( m/s )
  t_lm         ,   & ! temperature                                      (  K  )
  p_lm         ,   & ! full pressure (needed for lum2lm)                ( Pa  )
  pp_lm        ,   & ! deviation from the reference pressure            ( Pa  )
  qv_lm        ,   & ! specific water vapor content                     (kg/kg)
  qc_lm        ,   & ! specific cloud water content                     (kg/kg)
  qi_lm        ,   & ! cloud ice content                                (kg/kg)
  qr_lm        ,   & ! rain content                                     (kg/kg)
  qs_lm        ,   & ! snow content                                     (kg/kg)
  qg_lm        ,   & ! graupel content                                  (kg/kg)
  cgas_lm      ,   & !
  caero_lm     ,   & !
! iso code
  riso_lm      ,   & ! isotope ratios in water vapor
  risosoil_lm        ! isotope ratios in soil water
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
 fis_gme      ,         & ! orography * g                               (m2/s2)
 soiltyp_gme  ,         & ! type of the soil (keys 0-9)                   --
 fr_land_gme  ,         & ! land fraction of grid element               (  1  )
 z0_gme       ,         & ! surface roughness                           (  m  )
 plcmx_gme    ,         & ! vegetation: fraction covered by plants      (  1  )
 plcmn_gme    ,         & ! rest:       fraction covered by plants      (  1  )
 rlaimx_gme   ,         & ! vegetation: leaf area index                 (  1  )
 rlaimn_gme   ,         & ! rest:       leaf area index                 (  1  )
 root_gme     ,         & ! depth of the roots                          (  m  )
 vio3_gme     ,         & ! total vertically integrated ozone content   (Pa O3)
 hmo3_gme     ,         & ! height of maximum ozone concentration       ( Pa  )
 u_gme        ,         & ! zonal wind speed                            ( m/s )
 v_gme        ,         & ! meridional wind speed                       ( m/s )
 t_gme        ,         & ! temperature                                 (  K  )
 qv_gme       ,         & ! specific water vapor content                (kg/kg)
 qc_gme       ,         & ! specific cloud water content                (kg/kg)
 qi_gme       ,         & ! cloud ice content                           (kg/kg)
 ps_gme       ,         & ! surface pressure                            ( Pa  )
 t_s_gme      ,         & ! temperature of the ground surface           (  K  )
 t_snow_gme   ,         & ! temperature of the snow-surface             (  K  )
 t_m_gme      ,         & ! temperature between upper and medium
                          ! soil layer                                  (  K  )
 t_cl_gme     ,         & ! temp.  between medium and lower soil layer  (  K  )
 w_snow_gme   ,         & ! water content of the snow                   (m H2O)
 w_i_gme      ,         & ! water content of the interception storage   (m H2O)
 w_cl_gme     ,         & ! climatological deep soil water content      (m H2O)
 w_g1_gme     ,         & ! water content of the upper soil layer       (m H2O)
 w_g2_gme     ,         & ! water content of the medium soil layer      (m H2O)
 w_g3_gme     ,         & ! water content of the deepest soil layer     (m H2O)
 qv_s_gme     ,         & ! specific water vapor content on the surface (kg/kg)
 dpsdt_gme    ,         & ! surface pressure tendency                   (Pa/s )
 t_so_gme     ,         & ! temperature for new multi layer soil model  (  K  )
 w_so_gme     ,         & ! soil moisture for multi layer soil model    (m H2O)
 freshsnw_gme ,         & ! weighting function indicating 'freshness' of snow
 rho_snow_gme ,         & ! for prognostic treatment of snow density    (kg/m3)
 t_ice_gme    ,         & ! temperature of sea ice surface              (  K  )
 h_ice_gme                ! sea ice thickness                           (  m  )

!------------------------------------------------------------------------------

USE data_fields_in,  ONLY: &
 fis_in      ,          & ! orography * g                               (m2/s2)
 hsurf_in    ,          & ! orography                                   (  m  )
 soiltyp_in  ,          & ! type of the soil (keys 0-9)                   --
 ps_in       ,          & ! surface pressure                            ( Pa  )
 sst_in      ,          & ! sea surface temperature                     (  K  )
 t_s_in      ,          & ! temperature of the ground surface           (  K  )
 t_2m_in     ,          & ! 2m temperature                              (  K  )
 t_skin_in   ,          & ! skin temperature of the ground surface      (  K  )
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
 t_so_in     ,          & ! temperature for new multi layer soil model  (  K  )
 w_so_in     ,          & ! soil moisture for multi layer soil model    (m H2O)
 w_so_rel_in ,          & ! multi-layer relative   soil moisture        (  1  )
 freshsnw_in ,          & ! weighting function indicating 'freshness' of snow
 rho_snow_in ,          & ! for prognostic treatment of snow density    (kg/m3)
 t_ice_in    ,          & ! temperature of sea ice surface              (  K  )
 h_ice_in    ,          & ! sea ice thickness                           (  m  )
 u_in        ,          & ! zonal wind speed                            ( m/s )
 v_in        ,          & ! meridional wind speed                       ( m/s )
 w_in        ,          & ! vertical wind speed                         ( m/s )
 t_in        ,          & ! temperature                                 (  K  )
 p_in        ,          & ! standard pressure                           ( Pa  )
 pp_in       ,          & ! deviation from standard pressure            ( Pa  )
 qv_in       ,          & ! specific water vapor content                (kg/kg)
 qc_in       ,          & ! specific cloud water content                (kg/kg)
 qi_in       ,          & ! specific cloud ice water content            (kg/kg)
 qr_in       ,          & ! specific rain      content                  (kg/kg)
 qs_in       ,          & ! specific snow      content                  (kg/kg)
 qg_in                    ! specific graupel   content                  (kg/kg)

USE data_fields_in,  ONLY: &
 z0_in       ,          & ! surface roughness                           (  m  )
 fr_land_in  ,          & ! land fraction of grid element               (  1  )
 root_in     ,          & ! depth of the roots                          (  m  )
 plcov_in    ,          & ! fraction covered by plants                  (  1  )
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

PUBLIC   setup_vartab_grid_in, setup_vartab_lm

!==============================================================================

CONTAINS

!==============================================================================
!+ Creates the Grib 1 variable table for the LM output fields
!------------------------------------------------------------------------------

SUBROUTINE setup_vartab_lm (yerror, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  The table that is created in this routine is necessary for dealing with
!  Grib 1 I/O for the LM output. For every variable it defines the name, 
!  grib table type, level !  type, grib number, top and bottom level, 
!  grib factor, bias, data type, rang and corresponding targets (i.e. the 
!  arrays for the variables).
!
! Method:
!  Simple specifications.
!------------------------------------------------------------------------------

! Parameterlist
CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

INTEGER (KIND=iintegers)             :: istat

REAL(KIND=ireals), POINTER           :: dum3(:,:,:)
REAL(KIND=ireals), POINTER           :: dum2(:,:)

!- End of header
!------------------------------------------------------------------------------

  ierror    = 0_iintegers
  yerror    = '         '

!------------------------------------------------------------------------------
! Begin Subroutine setup_vartab_lm
!------------------------------------------------------------------------------

! Allocate the array var_lm
istat = 0

IF (l_art) THEN
  nvar_lm = nvar_lm_chem
ELSE
  nvar_lm = nvar_lm_norm
ENDIF

ALLOCATE(var_lm(nvar_lm), STAT=istat)
IF(istat /= 0) THEN
  ierror = 5
  yerror = ' Cannot allocate variable table for fine grid LM'
ENDIF

! the entries of the following structure are:
!   name, nr. of grib table, leveltype, element-numer, top-level, bottom-level,
!   factor, bias, rang of variable, pointer connection to actual field
!   units, standard_name, long_name, land sea dependency
var_lm( :)= ar_des_lm('          ','                              ', 0, 0, 0, 0, 0, 0.0, 0.0, 0, dum3, dum2, &
                      ' ',' ',' ',' ')

! external parameters
!------------------------------------------------------------------------------
var_lm( 1)= ar_des_lm('HSURF     ','surface'          ,   2,   1,   8,   0,   0,    1.0, 0.0, 2,dum3,   hsurf_lm, &
                     'm                  ','surface_altitude',                                &
                     'surface height',' ')
var_lm( 2)= ar_des_lm('FIS       ','surface'          ,   2,   1,   6,   0,   0,    1.0, 0.0, 2,dum3,     fis_lm, &
                     'm2 s-2              ','surface_geopotential_height',                    &
                     'surface geopotential height',' ')
var_lm( 3)= ar_des_lm('Z0        ','surface'          ,   2,   1,  83,   0,   0,    1.0, 0.0, 2,dum3,      z0_lm, &
                     'm                 ','surface_roughness_length',                         &
                     'surface roughness length','l')
var_lm( 4)= ar_des_lm('FR_LAND   ','surface'          ,   2,   1,  81,   0,   0,    1.0, 0.0, 2,dum3, fr_land_lm, &
                     '1                 ','land_area_fraction',                               &
                     'land-sea fraction',' ')
var_lm( 5)= ar_des_lm('SOILTYP   ','surface'          , 202,   1,  57,   0,   0,    1.0, 0.0, 2,dum3, soiltyp_lm, &
                     '1               ','soil_type',                                          &
                     'soil type',' ')
var_lm( 6)= ar_des_lm('PLCOV     ','surface'          ,   2,   1,  87,   0,   0,  100.0, 0.0, 2,dum3,   plcov_lm, &
                     '1                 ','vegetation_area_fraction',                         &
                     'vegetation area fraction','l')
var_lm( 7)= ar_des_lm('PLCOV_MX  ','surface'          , 202,   1,  67,   0,   0,    1.0, 0.0, 2,dum3,plcov_mx_lm, &
                     '1                 ','vegetation_area_fraction',                         &
                     'vegetation area fraction vegetation period','l')
var_lm( 8)= ar_des_lm('PLCOV_MN  ','surface'          , 202,   1,  68,   0,   0,    1.0, 0.0, 2,dum3,plcov_mn_lm, &
                     '1                 ','vegetation_area_fraction',                         &
                     'vegetation area fraction resting period','l')
var_lm( 9)= ar_des_lm('LAI       ','surface'          , 202,   1,  61,   0,   0,    1.0, 0.0, 2,dum3,     lai_lm, &
                     'm2 m-2            ','leaf_area_index',                                  &
                     'leaf area index','l')
var_lm(10)= ar_des_lm('LAI_MX    ','surface'          , 202,   1,  69,   0,   0,    1.0, 0.0, 2,dum3,  lai_mx_lm, &
                     'm2 m-2            ','leaf_area_index',                                  &
                     'leaf area index','l')
var_lm(11)= ar_des_lm('LAI_MN    ','surface'          , 202,   1,  70,   0,   0,    1.0, 0.0, 2,dum3,  lai_mn_lm, &
                     'm2 m-2            ','leaf_area_index_resting_period',                   &
                     'leaf area index resting period','l')
var_lm(12)= ar_des_lm('ROOTDP    ','surface'          , 202,   1,  62,   0,   0,    1.0, 0.0, 2,dum3,  rootdp_lm, &
                     'm               ','root_depth',                                         &
                     'root depth','l')
var_lm(13)= ar_des_lm('LATLM     ','surface'          , 202,   1, 114,   0,   0,    1.0, 0.0, 2,dum3,    latlm_m, &
                     'degrees_north   ','latitude',                                           &
                     'latitude',' ')
var_lm(14)= ar_des_lm('LONLM     ','surface'          , 202,   1, 115,   0,   0,    1.0, 0.0, 2,dum3,    lonlm_m, &
                     'degrees_east   ','longitude',                                           &
                     'longitude',' ')
var_lm(15)= ar_des_lm('FOR_E     ','surface'          , 202,   1,  75,   0,   0,    1.0, 0.0, 2,dum3,   for_e_lm, &
                     '1              ','-',                                                   &
                     'ground fraction covered by evergreen forest','l')
var_lm(16)= ar_des_lm('FOR_D     ','surface'          , 202,   1,  76,   0,   0,    1.0, 0.0, 2,dum3,   for_d_lm, &
                     '1              ','-',                                                   &
                     'ground fraction covered by deciduous forest','l')
var_lm(17)= ar_des_lm('SSO_STDH  ','surface'          , 202,   1,  46,   0,   0,    1.0, 0.0, 2,dum3,sso_stdh_lm, &
                     'm              ','-',                                                   &
                     'standard deviation of subgrid scale height','l')
var_lm(18)= ar_des_lm('SSO_GAMMA ','surface'          , 202,   1,  47,   0,   0,    1.0, 0.0, 2,dum3,sso_gamma_lm,&
                     '-              ','-',                                                   &
                     'anisotropy of topography', 'l')
var_lm(19)= ar_des_lm('SSO_THETA ','surface'          , 202,   1,  48,   0,   0,    1.0, 0.0, 2,dum3,sso_theta_lm,&
                     '-              ','-',                                                   &
                     'angle between principal axis of orography and rotated east','l')
var_lm(20)= ar_des_lm('SSO_SIGMA ','surface'          , 202,   1,  49,   0,   0,    1.0, 0.0, 2,dum3,sso_sigma_lm,&
                     '-              ','-',                                                   &
                     'mean slope of subgrid scale orography','l')
var_lm(21)= ar_des_lm('HORIZON   ','hybridLayer'      , 202, 110,  96,   0,   0,    1.0, 0.0, 3, horizon_lm, dum2,&
                     '-              ','-',                                                   &
                     'horizon angles for sectors around grid cell','l')
var_lm(22)= ar_des_lm('SLO_ANG   ','surface'          , 202,   1,  98,   0,   0,    1.0, 0.0, 2,dum3,  slo_ang_lm,&
                     '-              ','-',                                                   &
                     'slope angle - topography','l')
var_lm(23)= ar_des_lm('SLO_ASP   ','surface'          , 202,   1,  99,   0,   0,    1.0, 0.0, 2,dum3,  slo_asp_lm,&
                     '-              ','-',                                                   &
                     'slope aspect - topography','l')
var_lm(24)= ar_des_lm('SKYVIEW   ','surface'          , 202,   1, 100,   0,   0,    1.0, 0.0, 2,dum3,  skyview_lm,&
                     '-              ','-',                                                   &
                     'sky-view factor','l')
var_lm(25)= ar_des_lm('EMIS_RAD  ','surface'          , 202,   1,  56,   0,   0,    1.0, 0.0, 2,dum3, emis_rad_lm,&
                     '1              ','-',                                                   &
                     'thermal radiative surface emissivity','l')
var_lm(26)= ar_des_lm('RSMIN     ','surface'          , 201,   1, 212,   0,   0,    1.0, 0.0, 2,dum3,  prs_min_lm,&
                     's m-1          ','-',                                                   &
                     'minimum stomata resistance of plants','l')
var_lm(27)= ar_des_lm('NDVI_MRAT ','surface'          , 202,   1,  79,   0,   0,    1.0, 0.0, 3,ndvi_mrat_lm,dum2,&
                     '1              ','-',                                                   &
                     'monthly means of ratio of ndvi to annual maximum ndvi','l')
var_lm(28)= ar_des_lm('NDVIRATIO ','surface'          , 202,   1,  79,   0,   0,    1.0, 0.0, 2,dum3,ndviratio_lm,&
                     '1              ','-',                                                   &
                     'actual value of ratio of ndvi to annual maximum ndvi','l')
var_lm(29)= ar_des_lm('PLCOV12   ','surface'          ,   2,   1,  87,   0,   0,  100.0, 0.0, 3,plcov12   , dum2, &
                     '1                 ','vegetation_area_fraction',                         &
                     'vegetation area fraction','l')
var_lm(30)= ar_des_lm('LAI12     ','surface'          , 202,   1,  61,   0,   0,    1.0, 0.0, 3,lai12   ,   dum2, &
                     'm2 m-2            ','leaf_area_index',                                  &
                     'leaf area index','l')
var_lm(31)= ar_des_lm('Z012      ','surface'          ,   2,   1,  83,   0,   0,    1.0, 0.0, 3,z012   ,    dum2, &
                     'm                 ','surface_roughness_length',                         &
                     'surface roughness length','l')
var_lm(32)= ar_des_lm('AER_SO412 ','surface'          , 202,   1,  84,   0,   0,    1.0, 0.0, 3,aer_su12_lm,dum2, &
                     '1                 ','-',                                                &
                     'monthly mean for aerosol type sulfate drops',' ')
var_lm(33)= ar_des_lm('AER_DUST12','surface'          , 202,   1,  86,   0,   0,    1.0, 0.0, 3,aer_du12_lm,dum2, &
                     '1                 ','-',                                                &
                     'monthly mean for aerosol type mineral dust',' ')
var_lm(34)= ar_des_lm('AER_ORG12 ','surface'          , 202,   1,  91,   0,   0,    1.0, 0.0, 3,aer_or12_lm,dum2, &
                     '1                 ','-',                                                &
                     'monthly mean for aerosol type organic',' ')
var_lm(35)= ar_des_lm('AER_BC12  ','surface'          , 202,   1,  92,   0,   0,    1.0, 0.0, 3,aer_bc12_lm,dum2, &
                     '1                 ','-',                                                &
                     'monthly mean for aerosol type black carbon',' ')
var_lm(36)= ar_des_lm('AER_SS12  ','surface'          , 202,   1,  93,   0,   0,    1.0, 0.0, 3,aer_ss12_lm,dum2, &
                     '1                 ','-',                                                &
                     'monthly mean for aerosol type sea salt',' ')
var_lm(37)= ar_des_lm('AER_SO4   ','surface'          , 202,   1,  84,   0,   0,    1.0, 0.0, 2,dum3,   aer_su_lm,&
                     '1                 ','-',                                                &
                     'aerosol type sulfate drops',' ')
var_lm(38)= ar_des_lm('AER_DUST  ','surface'          , 202,   1,  86,   0,   0,    1.0, 0.0, 2,dum3,   aer_du_lm,&
                     '1                 ','-',                                                &
                     'aerosol type mineral dust',' ')
var_lm(39)= ar_des_lm('AER_ORG   ','surface'          , 202,   1,  91,   0,   0,    1.0, 0.0, 2,dum3,   aer_or_lm,&
                     '1                 ','-',                                                &
                     'aerosol type organic',' ')
var_lm(40)= ar_des_lm('AER_BC    ','surface'          , 202,   1,  92,   0,   0,    1.0, 0.0, 2,dum3,   aer_bc_lm,&
                     '1                 ','-',                                                &
                     'aerosol type black carbon',' ')
var_lm(41)= ar_des_lm('AER_SS    ','surface'          , 202,   1,  93,   0,   0,    1.0, 0.0, 2,dum3,   aer_ss_lm,&
                     '1                 ','-',                                                &
                     'aerosol type sea salt',' ')
var_lm(42)= ar_des_lm('URBAN     ','surface'          , 250,   1,  85,   0,   0,    1.0, 0.0, 2,dum3, fr_urban_lm,&
                     '1              ','-',                                                   &
                     'ground fraction covered by urban land','l')
var_lm(43)= ar_des_lm('ALB_DRY   ','surface'          , 202,   1, 127,   0,   0,  100.0, 0.0, 2,dum3, alb_dry_lm, &
                     '1           ','soil_albedo',                                            &
                      'dry soil albedo',  'l')
var_lm(44)= ar_des_lm('ALB_SAT   ','surface'          , 202,   1, 128,   0,   0,  100.0, 0.0, 2,dum3, alb_sat_lm, &
                     '1           ','soil_albedo',                                            &
                     'saturated soil albedo',  'l')
var_lm(45)= ar_des_lm('ALB_DIF12 ','surface'          , 202,   1, 129,   0,   0,  100.0, 0.0, 3,alb_dif12_lm,dum2,&
                     '1                 ','-',                                                &
                     'monthly values of solar diffuse surface albedo',' ')
var_lm(46)= ar_des_lm('ALB_DIF   ','surface'          , 202,   1, 129,   0,   0,  100.0, 0.0, 2,dum3, alb_dif_lm, &
                     '1           ','solar_albedo',                                           &
                     'solar surface albedo',  'l')


! multi level fields
var_lm(51)= ar_des_lm('U         ','hybridLayer'      ,   2, 110,  33,   0,   0,    1.0, 0.0, 3,   u_lm,    dum2, &
                     'm s-1              ','grid_eastward_wind',                              &
                     'U-component of wind',' ')
var_lm(52)= ar_des_lm('V         ','hybridLayer'      ,   2, 110,  34,   0,   0,    1.0, 0.0, 3,   v_lm,    dum2, &
                     'm s-1              ','grid_northward_wind',                             &
                     'V-component of wind',' ')
var_lm(53)= ar_des_lm('W         ','hybrid'           ,   2, 109,  40,   0,   0,    1.0, 0.0, 3,   w_lm,    dum2, &
                     'm s-1             ', 'upward_air_velocity',                             &
                     'vertical wind velocity',' ')
var_lm(54)= ar_des_lm('T         ','hybridLayer'      ,   2, 110,  11,   0,   0,    1.0, 0.0, 3,   t_lm,    dum2, &
                     'K                  ','air_temperature',                                 &
                     'temperature',' ')
var_lm(55)= ar_des_lm('P         ','hybridLayer'      ,   2, 110,   1,   0,   0,    1.0, 0.0, 3,   p_lm,    dum2, &
                     'Pa                  ','air_pressure',                                   &
                     'pressure',' ')
var_lm(56)= ar_des_lm('PP        ','hybridLayer'      , 201, 110, 139,   0,   0,   0.01, 0.0, 3,  pp_lm,    dum2, &
                     'Pa              ','difference_of_air_pressure_from_model_reference',    &
                     'deviation from reference pressure',' ')
var_lm(57)= ar_des_lm('QV        ','hybridLayer'      ,   2, 110,  51,   0,   0,    1.0, 0.0, 3,  qv_lm,    dum2, &
                     'kg kg-1           ','specific_humidity',                                &
                     'specific humidity',' ')
var_lm(58)= ar_des_lm('QC        ','hybridLayer'      , 201, 110,  31,   0,   0,    1.0, 0.0, 3,  qc_lm,    dum2, &
                     'kg kg-1          ','mass_fraction_of_cloud_liquid_water_in_air',        &
                     'cloud liquid water content',' ')
var_lm(59)= ar_des_lm('QI        ','hybridLayer'      , 201, 110,  33,   0,   0,    1.0, 0.0, 3,  qi_lm,    dum2, &
                     'kg kg-1          ','mass_fraction_of_cloud_ice_in_air',                 &
                     'cloud ice content',' ')
var_lm(60)= ar_des_lm('QR        ','hybridLayer'      , 201, 110,  35,   0,   0,    1.0, 0.0, 3,  qr_lm,    dum2, &
                     'kg kg-1            ','mass_fraction_of_rain_in_air',                    &
                     'specific rain content',' ')
var_lm(61)= ar_des_lm('QS        ','hybridLayer'      , 201, 110,  36,   0,   0,    1.0, 0.0, 3,  qs_lm,    dum2, &
                     'kg kg-1            ','mass_fraction_of_snow_in_air',                    &
                     'specific snow content',' ')
var_lm(62)= ar_des_lm('QG        ','hybridLayer'      , 201, 110,  39,   0,   0,    1.0, 0.0, 3,  qg_lm,    dum2, &
                     'kg kg-1            ','mass_fraction_of_graupel_in_air',                 &
                     'specific graupel content',' ')
var_lm(63)= ar_des_lm('HHL       ','hybrid'           ,   2, 109,   8,   0,   0,    1.0, 0.0, 3, hhl_lm,    dum2, &
                     'm                  ','altitude',                                        &
                     'height',' ')
var_lm(64)= ar_des_lm('T_SO      ','depthBelowLand'   , 201, 111, 197,   0,   0,    1.0, 0.0, 3,t_so_lm,    dum2, &
                     'K               ','soil_temperature',                                   &
                     'soil temperature','l')
var_lm(65)= ar_des_lm('W_SO      ','depthBelowLand'   , 201, 111, 198,   0,   0, 1000.0, 0.0, 3,w_so_lm,    dum2, &
                     'm                 ','lwe_thickness_of_moisture_content_of_soil_layer',  &
                     'soil water content','l')

! single level fields
var_lm(71)= ar_des_lm('PS        ','surface'          ,   2,   1,   1,   0,   0,    1.0, 0.0, 2,dum3,      ps_gl, &
                     'Pa                  ','surface_air_pressure',                           &
                     'surface pressure',' ')
var_lm(72)= ar_des_lm('T_SNOW    ','surface'          , 201,   1, 203,   0,   0,    1.0, 0.0, 2,dum3,  t_snow_lm, &
                     'K               ','surface_temperature_where_snow',                     &
                     'snow temperature','l')
var_lm(73)= ar_des_lm('T_S       ','surface'          ,   2, 111,  85,   0,   0,    1.0, 0.0, 2,dum3,     t_s_lm, &
                     'K                 ','-',                                                &
                     'soil surface temperature',' ')
var_lm(74)= ar_des_lm('T_M       ','depthBelowLand'   ,   2, 111,  85,   0,   9,    1.0, 0.0, 2,dum3,     t_m_lm, &
                     'K                 ','soil_temperature',                                 &
                     'temperature of 1. soil layer','l')
var_lm(75)= ar_des_lm('T_CL      ','depthBelowLand'   ,   2, 111,  85,   0,  41,    1.0, 0.0, 2,dum3,    t_cl_lm, &
                     'K                 ','soil_temperature',                                 &
                     'deep soil temperature','l')
var_lm(76)= ar_des_lm('W_SNOW    ','surface'          ,   2,   1,  65,   0,   0, 1000.0, 0.0, 2,dum3,  w_snow_lm, &
                     'm                 ','lwe_thickness_of_surface_snow_amount',             &
                     'surface snow amount','l')
var_lm(77)= ar_des_lm('W_I       ','surface'          , 201,   1, 200,   0,   0, 1000.0, 0.0, 2,dum3,     w_i_lm, &
                     'm                 ','lwe_thickness_of_canopy_water_amount',             &
                     'canopy water amount','l')

IF (nl_soil_lm == 2) THEN
var_lm(78)= ar_des_lm('W_G1      ','depthBelowLandLayer', 2, 112,  86,   0,  10, 1000.0, 0.0, 2,dum3,    w_g1_lm, &
                     'm                 ','lwe_thickness_of_moisture_content_of_soil_layer',  &
                     'water content of 1. soil layer','l')
var_lm(79)= ar_des_lm('W_G2      ','depthBelowLandLayer', 2, 112,  86,  10, 100, 1000.0, 0.0, 2,dum3,    w_g2_lm, &
                     'm                 ','lwe_thickness_of_moisture_content_of_soil_layer',  &
                     'water content of 2. soil layer','l')
ELSEIF (nl_soil_lm == 3) THEN
var_lm(78)= ar_des_lm('W_G1      ','depthBelowLandLayer', 2, 112,  86,   0,   2, 1000.0, 0.0, 2,dum3,    w_g1_lm, &
                     'm                 ','lwe_thickness_of_moisture_content_of_soil_layer',  &
                     'water content of 1. soil layer','l')
var_lm(79)= ar_des_lm('W_G2      ','depthBelowLandLayer', 2, 112,  86,   2,  10, 1000.0, 0.0, 2,dum3,    w_g2_lm, &
                     'm                 ','lwe_thickness_of_moisture_content_of_soil_layer',  &
                     'water content of 2. soil layer','l')
ENDIF

var_lm(80)= ar_des_lm('W_G3      ','depthBelowLandLayer', 2, 112,  86,  10, 100, 1000.0, 0.0, 2,dum3,    w_g3_lm, &
                     'm                 ','lwe_thickness_of_moisture_content_of_soil_layer',  &
                     'water content of 3. soil layer','l')
var_lm(81)= ar_des_lm('W_CL      ','depthBelowLandLayer', 2, 112,  86, 100, 190, 1000.0, 0.0, 2,dum3,    w_cl_lm, &
                     'm                 ','lwe_thickness_of_moisture_content_of_soil_layer',  &
                     'deep soil water','l')
var_lm(82)= ar_des_lm('QV_S      ','surface'          ,   2,   1,  51,   0,   0,    1.0, 0.0, 2,dum3,    qv_s_lm, &
                     'kg kg-1           ','surface_specific_humidity',                        &
                     'surface specific humidity',' ')
var_lm(83)= ar_des_lm('VIO3      ','surface'          , 202,   1,  65,   0,   0,    1.0, 0.0, 2,dum3,    vio3_lm, &
                     'Pa              ','equivalent_pressure_of_atmosphere_ozone_content',    &
                     'vertical integrated ozone amount',' ')
var_lm(84)= ar_des_lm('HMO3      ','surface'          , 202,   1,  64,   0,   0,    1.0, 0.0, 2,dum3,    hmo3_lm, &
                     'Pa              ','air_pressure',                                       &
                     'air pressure at ozone maximum',' ')
var_lm(85)= ar_des_lm('FRESHSNW  ','surface'          , 201,   1, 129,   0,   0,    1.0, 0.0, 2,dum3,freshsnw_lm, &
                     '1              ','-',                                                   &
                     'weighting function indicating freshness of snow','l')
var_lm(86)= ar_des_lm('RHO_SNOW  ','surface'          , 201,   1, 133,   0,   0,    1.0, 0.0, 2,dum3,rho_snow_lm, &
                     'kg m-3         ','surface_snow_density',                                &
                     'density of snow',' ')
var_lm(87)= ar_des_lm('FR_LAKE   ','surface'          , 202,   1,  55,   0,   0,    1.0, 0.0, 2,dum3,fr_lake_lm,  &
                     '1                 ','-',                                                &
                     'lake area fraction',' ')
var_lm(88)= ar_des_lm('DEPTH_LK  ','surface'          , 201,   1,  96,   0,   0,    1.0, 0.0, 2,dum3,depth_lk_lm, &
                     'm                 ','sea_floor_depth_below_sea_level',                  &
                     'lake depth','i')
var_lm(89)= ar_des_lm('T_ICE     ','surface'          , 201,   1, 215,   0,   0,    1.0, 0.0, 2,dum3,   t_ice_lm, &
                     'K                 ','sea_ice_temperature',                              &
                     'temperature of sea ice surface',' ')
var_lm(90)= ar_des_lm('H_ICE     ','surface'          ,   2,   1,  92,   0,   0,    1.0, 0.0, 2,dum3,   h_ice_lm, &
                     'm                 ','sea_ice_thickness',                                &
                     'sea ice thickness',' ')
var_lm(91)= ar_des_lm('T_B1_LK   ','surface'          , 201,   1, 192,   0,   0,    1.0, 0.0, 2,dum3, t_b1_lk_lm, &
                     'K                 ','-',                                                &
                     'temperature at bottom of upper layer of sediments','i')
var_lm(92)= ar_des_lm('H_B1_LK   ','surface'          , 201,   1,  94,   0,   0,    1.0, 0.0, 2,dum3, h_b1_lk_lm, &
                     'm                 ','-',                                                &
                     'thickness of the upper layer of bottom sediments','i')
var_lm(93)= ar_des_lm('T_WML_LK  ','surface'          , 201,   1, 193,   0,   0,    1.0, 0.0, 2,dum3,t_wml_lk_lm, &
                     'K                 ','-',                                                &
                     'mixed layer temperature','i')
var_lm(94)= ar_des_lm('T_MNW_LK  ','surface'          , 201,   1, 194,   0,   0,    1.0, 0.0, 2,dum3,t_mnw_lk_lm, &
                     'K                 ','-',                                                &
                     'mean temperature of water column','i')
var_lm(95)= ar_des_lm('T_BOT_LK  ','surface'          , 201,   1, 191,   0,   0,    1.0, 0.0, 2,dum3,t_bot_lk_lm, &
                     'K                 ','-',                                                &
                     'temperature at water bottom sediment interface','i')
var_lm(96)= ar_des_lm('C_T_LK    ','surface'          , 201,   1,  91,   0,   0,    1.0, 0.0, 2,dum3,  c_t_lk_lm, &
                     '1                 ','-',                                                &
                     'shape factor of temperature profile in lake thermocline','i')
var_lm(97)= ar_des_lm('H_ML_LK   ','surface'          , 201,   1,  95,   0,   0,    1.0, 0.0, 2,dum3, h_ml_lk_lm, &
                     'm                 ','-',                                                &
                     'thickness of mixed layer','i')
var_lm(98)= ar_des_lm('SALT_LK   ','surface'          , 2,   1,  88,   0,   0,    1.0, 0.0, 2,dum3,salt_lk_lm,  &
                     'g kg-1            ','sea_water_salinity',                  &
                     'lake salinity','i')
! iso code
IF (liso) THEN
var_lm(99) = ar_des_lm('QV18O  ', 'hybridLayer'     ,206, 110,   2,   0,   0,    1.0, 0.0, 3, riso_lm(:,:,:,1),&
                       dum2, '1/SMOW          ','-', 'vapor 18O content', ' ')
var_lm(100) = ar_des_lm('QV2H  ', 'hybridLayer'     ,206, 110,   3,   0,   0,    1.0, 0.0, 3, riso_lm(:,:,:,2),&
                       dum2, '1/SMOW          ','-', 'vapor 2H content',  ' ')
var_lm(101)= ar_des_lm('R18OSOIL ', 'surface'         ,206,   1, 162,   0,   0,    1.0, 0.0, 2, dum3,            &
                       risosoil_lm(:,:,1), '1/SMOW          ','-', 'soil moist. 18O content', &
                       ' ')
var_lm(102)= ar_des_lm('R2HSOIL  ', 'surface'         ,206,   1, 163,   0,   0,    1.0, 0.0, 2, dum3,            &
                       risosoil_lm(:,:,2), '1/SMOW          ','-', 'soil moist. 2H content',  &
                       ' ')
ENDIF
! end iso code

END SUBROUTINE setup_vartab_lm

!==============================================================================
!+ Creates the variable table for the input models
!------------------------------------------------------------------------------

SUBROUTINE setup_vartab_grid_in (yerror, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  The table that is created in this routine is necessary for dealing with
!  Grib 1 I/O. For every variable it defines the name, grib table type, level
!  type, grib number, top and bottom level, grib factor, bias, data type,
!  time range indicator, rang and corresponding targets (i.e. the arrays for
!  the variables).
!
! Method:
!  Simple specifications.
!------------------------------------------------------------------------------

! Parameterlist
CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

INTEGER (KIND=iintegers)             :: istat, izdate, i
INTEGER (KIND=iintegers)             :: iec=128  ! EC Grib table

REAL(KIND=ireals), POINTER           :: dum4(:,:,:,:)
REAL(KIND=ireals), POINTER           :: dum3(:,:,:)
REAL(KIND=ireals), POINTER           :: dum2(:,:)

!- End of header
!------------------------------------------------------------------------------

  ierror    = 0_iintegers
  yerror    = '         '

  ! Get initial date from ydate_ini='yyyymmddhh'
  READ(ydate_ini(1:8),'(I8)') izdate

  ! Allocate the array var_in
  istat = 0
  IF (l_art .AND. l_art_nested) THEN
    nvar_in = nvar_in_chem
  ELSE
    nvar_in = nvar_in_norm
  ENDIF

  ALLOCATE(var_in(nvar_in), STAT=istat)
  IF(istat /= 0) THEN
    ierror = 1
    yerror = 'Error allocating var_in for input grid with dimension:  '
    WRITE (yerror(60:64),'(I5)') nvar_in
    WRITE (yerror(66:70),'(I5)') istat
    RETURN
  ENDIF

! the entries of the following structure are:
!   name, nr. of grib table, leveltype, element-numer, top-level, bottom-level,
!   factor, bias, rang of variable, pointer connection to actual field,
!   data-type, interpolation code, logical indicator whether variable has been read,
!   nr. of levels to be read, nr. of levels that have been read

var_in(:) = ar_des_input('        ', '                              ', 0, 0,  &
             0, 0, 0, 0.0, 0.0, 0, dum4, dum3, dum2,'   ','   ',.FALSE., 0 , 0)

!------------------------------------------------------------------------------
! Section 1: Basic settings
!------------------------------------------------------------------------------

! constant fields
var_in( 1) = ar_des_input('FIS       ','surface'          ,   2,   1,   6,   0,   0,    1.0, 0.0, 2, &
                          dum4,      fis_gme,     fis_in,'  E','LFF',.FALSE.,     1,      0)
var_in( 2) = ar_des_input('FIS_SH    ','hybrid'           , iec, 109, 129,   0,   0,    1.0, 0.0, 2, &
                          dum4,         dum3,     fis_in,'   ','LFF',.FALSE.,     1,      0)
var_in( 3) = ar_des_input('HSURF     ','surface'          ,   2,   1,   8,   0,   0,    1.0, 0.0, 2, &
                          dum4,         dum3,   hsurf_in,'  E','LFF',.FALSE.,     1,      0)
var_in( 4) = ar_des_input('Z0        ','surface'          ,   2,   1,  83,   0,   0,    1.0, 0.0, 2, &
                          dum4,       z0_gme,      z0_in,'   ','MFF',.FALSE.,     1,      0)
var_in( 5) = ar_des_input('FR_LAND   ','surface'          ,   2,   1,  81,   0,   0,    1.0, 0.0, 2, &
                          dum4,  fr_land_gme, fr_land_in,'  E','MFF',.FALSE.,     1,      0)
var_in( 6) = ar_des_input('SOILTYP   ','surface'          , 202,   1,  57,   0,   0,    1.0, 0.0, 2, &
                          dum4,  soiltyp_gme, soiltyp_in,'   ','NFF',.FALSE.,     1,      0)
var_in( 7) = ar_des_input('PLCOV     ','surface'          ,   2,   1,  87,   0,   0,  100.0, 0.0, 2, &
                          dum4,         dum3,   plcov_in,'   ','MFF',.FALSE.,     1,      0)
var_in( 8) = ar_des_input('PLCOV_MX  ','surface'          , 202,   1,  67,   0,   0,  100.0, 0.0, 2, &
                          dum4,    plcmx_gme,   plcmx_in,'   ','MFF',.FALSE.,     1,      0)
var_in( 9) = ar_des_input('PLCOV_MN  ','surface'          , 202,   1,  68,   0,   0,  100.0, 0.0, 2, &
                          dum4,    plcmn_gme,   plcmn_in,'   ','MFF',.FALSE.,     1,      0)
var_in(10) = ar_des_input('LAI       ','surface'          , 202,   1,  61,   0,   0,    1.0, 0.0, 2, &
                          dum4,         dum3,    rlai_in,'   ','MFF',.FALSE.,     1,      0)
var_in(11) = ar_des_input('LAI_MX    ','surface'          , 202,   1,  69,   0,   0,  100.0, 0.0, 2, &
                          dum4,   rlaimx_gme,  rlaimx_in,'   ','MFF',.FALSE.,     1,      0)
var_in(12) = ar_des_input('LAI_MN    ','surface'          , 202,   1,  70,   0,   0,  100.0, 0.0, 2, &
                          dum4,   rlaimn_gme,  rlaimn_in,'   ','MFF',.FALSE.,     1,      0)
var_in(13) = ar_des_input('ROOTDP    ','surface'          , 202,   1,  62,   0,   0,    1.0, 0.0, 2, &
                          dum4,     root_gme,    root_in,'   ','MFF',.FALSE.,     1,      0)
var_in(14) = ar_des_input('VIO3      ','surface'          , 202,   1,  65,   0,   0,    1.0, 0.0, 2, &
                          dum4,     vio3_gme,    vio3_in,'   ','LFF',.FALSE.,     1,      0)
var_in(15) = ar_des_input('HMO3      ','surface'          , 202,   1,  64,   0,   0,    1.0, 0.0, 2, &
                          dum4,     hmo3_gme,    hmo3_in,'   ','LFF',.FALSE.,     1,      0)

! multi level fields
var_in(21) = ar_des_input('U         ','hybridLayer'      ,   2, 110,  33,   0,   0,    1.0, 0.0, 3, &
                          u_gme,        u_in,       dum2,'IB ','QFF',.FALSE., ke_in,      0)
var_in(22) = ar_des_input('V         ','hybridLayer'      ,   2, 110,  34,   0,   0,    1.0, 0.0, 3, &
                          v_gme,        v_in,       dum2,'IB ','QFF',.FALSE., ke_in,      0)
var_in(23) = ar_des_input('W         ','hybrid'           ,   2, 109,  40,   0,   0,    1.0, 0.0, 3, &
                          dum4,         w_in,       dum2,'   ','QFF',.FALSE., ke1in,      0)
var_in(24) = ar_des_input('T         ','hybridLayer'      ,   2, 110,  11,   0,   0,    1.0, 0.0, 3, &
                          t_gme,        t_in,       dum2,'IB ','QFF',.FALSE., ke_in,      0)
var_in(25) = ar_des_input('P         ','hybridLayer'      ,   2, 110,   1,   0,   0,    1.0, 0.0, 3, &
                          dum4,         p_in,       dum2,'   ','QFF',.FALSE., ke_in,      0)
var_in(26) = ar_des_input('PP        ','hybridLayer'      , 201, 110, 139,   0,   0,   0.01, 0.0, 3, &
                          dum4,        pp_in,       dum2,'   ','QFF',.FALSE., ke_in,      0)
var_in(27) = ar_des_input('QV        ','hybridLayer'      ,   2, 110,  51,   0,   0,    1.0, 0.0, 3, &
                          qv_gme,      qv_in,       dum2,'IB ','QFT',.FALSE., ke_in,      0)
var_in(28) = ar_des_input('QC        ','hybridLayer'      , 201, 110,  31,   0,   0,    1.0, 0.0, 3, &
                          qc_gme,      qc_in,       dum2,'IB ','QFT',.FALSE., ke_in,      0)
var_in(29) = ar_des_input('QI        ','hybridLayer'      , 201, 110,  33,   0,   0,    1.0, 0.0, 3, &
                          qi_gme,      qi_in,       dum2,'   ','QTT',.FALSE., ke_in,      0)
var_in(30) = ar_des_input('QR        ','hybridLayer'      , 201, 110,  35,   0,   0,    1.0, 0.0, 3, &
                          dum4,        qr_in,       dum2,'   ','QTT',.FALSE., ke_in,      0)
var_in(31) = ar_des_input('QS        ','hybridLayer'      , 201, 110,  36,   0,   0,    1.0, 0.0, 3, &
                          dum4,        qs_in,       dum2,'   ','QTT',.FALSE., ke_in,      0)
var_in(32) = ar_des_input('QG        ','hybridLayer'      , 201, 110,  39,   0,   0,    1.0, 0.0, 3, &
                          dum4,        qg_in,       dum2,'   ','QTT',.FALSE., ke_in,      0)
var_in(33) = ar_des_input('T_SO      ','depthBelowLand'   , 201, 111, 197,   0,   0,    1.0, 0.0, 3, &
                          t_so_gme,  t_so_in,       dum2,'   ','MFF',.FALSE., ke_soil_in, 0)
var_in(34) = ar_des_input('W_SO      ','depthBelowLand'   , 201, 111, 198,   0,   0, 1000.0, 0.0, 3, &
                          w_so_gme,  w_so_in,       dum2,'   ','MFF',.FALSE., ke_soil_in, 0)
var_in(35) = ar_des_input('W_SO_REL  ','depthBelowLand'   , iec, 111,  39,   0,   0,    1.0, 0.0, 3, &
                          dum4,      w_so_rel_in,   dum2,'   ','MFF',.FALSE., ke_soil_in, 0)

! single level fields
var_in(41) = ar_des_input('PS        ','surface'          ,   2,   1,   1,   0,   0,    1.0, 0.0, 2, &
                          dum4,       ps_gme,      ps_in,'   ','LFF',.FALSE.,     1,      0)
var_in(42) = ar_des_input('LNPS      ','hybrid'           , iec, 109,   1,   0,   0,    1.0, 0.0, 2, &
                          dum4,         dum3,      ps_in,'   ','LFF',.FALSE.,     1,      0)
var_in(43) = ar_des_input('T_SNOW    ','surface'          , 201,   1, 203,   0,   0,    1.0, 0.0, 2, &
                          dum4,   t_snow_gme,  t_snow_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(44) = ar_des_input('T_S       ','depthBelowLand'   ,   2, 111,  85,   0,   0,    1.0, 0.0, 2, &
                          dum4,      t_s_gme,     t_s_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(45) = ar_des_input('SST       ','surface'          , iec,   1,  34,   0,   0,    1.0, 0.0, 2, &
                          dum4,         dum3,     t_s_in,'   ','MFF',.FALSE.,     1,      0)
var_in(46) = ar_des_input('T_SKIN    ','surface'          , iec,   1, 235,   0,   0,    1.0, 0.0, 2, &
                          dum4,         dum3,     t_s_in,'   ','MFF',.FALSE.,     1,      0)
var_in(47) = ar_des_input('T_2M      ','heightAboveGround',   2, 105,  11,   0,   0,    1.0, 0.0, 2, &
                          dum4,         dum3,    t_2m_in,'   ','MFF',.FALSE.,     1,      0)
var_in(48) = ar_des_input('T_M       ','depthBelowLand'   ,   2, 111,  85,   0,   9,    1.0, 0.0, 2, &
                          dum4,      t_m_gme,     t_m_in,'   ','MFF',.FALSE.,     1,      0)
var_in(49) = ar_des_input('T_G1      ','depthBelowLandLayer', iec, 112, 139, 0,   7,    1.0, 0.0, 2, &
                          dum4,         dum3,    t_g1_in,'   ','MFF',.FALSE.,     1,      0)
var_in(50) = ar_des_input('T_G2      ','depthBelowLandLayer', iec, 112, 170, 7,  28,    1.0, 0.0, 2, &
                          dum4,         dum3,    t_g2_in,'   ','MFF',.FALSE.,     1,      0)
var_in(51) = ar_des_input('T_G3      ','depthBelowLandLayer', iec, 112, 183,28, 100,    1.0, 0.0, 2, &
                          dum4,         dum3,    t_g3_in,'   ','MFF',.FALSE.,     1,      0)
var_in(52) = ar_des_input('T_CL      ','depthBelowLand',      2, 111,  85,   0,  41,    1.0, 0.0, 2, &
                          dum4,     t_cl_gme,    t_cl_in,'   ','MFF',.FALSE.,     1,      0)
var_in(53) = ar_des_input('W_SNOW    ','surface'          ,   2,   1,  65,   0,   0, 1000.0, 0.0, 2, &
                          dum4,   w_snow_gme,  w_snow_in,'IB ','MFT',.FALSE.,     1,      0)
var_in(54) = ar_des_input('W_I       ','surface'          , 201,   1, 200,   0,   0, 1000.0, 0.0, 2, &
                          dum4,      w_i_gme,     w_i_in,'I  ','MFF',.FALSE.,     1,      0)
var_in(55) = ar_des_input('W_G1      ','depthBelowLandLayer', 2, 112,  86,   0,  10, 1000.0, 0.0, 2, &
                          dum4,     w_g1_gme,    w_g1_in,'   ','MFF',.FALSE.,     1,      0)
var_in(56) = ar_des_input('W_G2      ','depthBelowLandLayer', 2, 112,  86,  10, 100, 1000.0, 0.0, 2, &
                          dum4,     w_g2_gme,    w_g2_in,'   ','MFF',.FALSE.,     1,      0)
var_in(57) = ar_des_input('W_G3      ','depthBelowLandLayer', 2, 112,  86,  10, 100, 1000.0, 0.0, 2, &
                          dum4,     w_g3_gme,    w_g3_in,'   ','MFF',.FALSE.,     1,      0)
var_in(58) = ar_des_input('W_CL      ','depthBelowLandLayer', 2, 112,  86, 100, 190, 1000.0, 0.0, 2, &
                          dum4,     w_cl_gme,    w_cl_in,'   ','MFF',.FALSE.,     1,      0)
var_in(59) = ar_des_input('QV_S      ','surface'          ,   2,   1,  51,   0,   0,    1.0, 0.0, 2, &
                          dum4,     qv_s_gme,    qv_s_in,'IB ','MFT',.FALSE.,     1,      0)
var_in(60) = ar_des_input('QV_2M     ','heightAboveGround',   2, 105,  51,   0,   0,    1.0, 0.0, 2, &
                          dum4,         dum3,   qv_2m_in,'   ','MFF',.FALSE.,     1,      0)
var_in(61) = ar_des_input('FRESHSNW  ','surface'          , 201,   1, 129,   0,   0,    1.0, 0.0, 2, &
                          dum4, freshsnw_gme,freshsnw_in,'   ','MFF',.FALSE.,     1,      0)
var_in(62) = ar_des_input('FI        ','hybrid'           ,   2, 109,   6,   0,   0,    1.0, 0.0, 2, &
                          dum4,         dum3,     fic_in,'IB ','LFF',.FALSE.,     1,      0)
var_in(63) = ar_des_input('RHO_SNOW  ','surface'          , 201,   1, 133,   0,   0,    1.0, 0.0, 2, &
                          dum4, rho_snow_gme,rho_snow_in,'   ','MFF',.FALSE.,     1,      0)
var_in(64) = ar_des_input('T_ICE     ','surface'          , 201,   1, 215,   0,   0,    1.0, 0.0, 2, &
                          dum4,    t_ice_gme,   t_ice_in,'   ','MFF',.FALSE.,     1,      0)
var_in(65) = ar_des_input('H_ICE     ','surface'          ,   2,   1,  92,   0,   0,    1.0, 0.0, 2, &
                          dum4,    h_ice_gme,   h_ice_in,'   ','MFF',.FALSE.,     1,      0)

!US Question:
!   The ECMWF products were coded with ilevtyp=109 (hybrid layers) in the former Grib1 tables.
!   But I do not see that these are 3D fields, I think they are surface fields?!?!

!------------------------------------------------------------------------------
! Section 2: Change tabtyp, levtyp, element number for lec2lm, lcm2lm
!------------------------------------------------------------------------------

! for lec2lm and lcm2lm, many things have to be changed:
IF (lec2lm .OR. lcm2lm) THEN
var_in( 1) = ar_des_input('FIS       ','surface'          , iec,   1, 129,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,     fis_in,'  E','LFF',.FALSE.,     1,      0)
var_in( 2) = ar_des_input('FIS_SH    ','hybrid'           , iec, 109, 129,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,     fis_in,'  E','LFF',.FALSE.,     1,      0)
var_in( 4) = ar_des_input('Z0        ','surface'          , iec,   1, 173,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,      z0_in,'   ','MFF',.FALSE.,     1,      0)
var_in( 5) = ar_des_input('FR_LAND   ','surface'          , iec,   1, 172,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3, fr_land_in,'  E','MFF',.FALSE.,     1,      0)
var_in( 6) = ar_des_input('SOILTYP   ','surface'          , iec,   1,  43,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3, soiltyp_in,'   ','NFF',.FALSE.,     1,      0)
var_in( 7) = ar_des_input('PLCOV     ','surface'          , iec,   1, 199,   0,   0,  100.0, 0.0, 2, &
                          dum4,    dum3,   plcov_in,'   ','MFF',.FALSE.,     1,      0)
var_in(21) = ar_des_input('U         ','hybrid'           , iec, 109, 131,   0,   0,    1.0, 0.0, 3, &
                          dum4,    u_in,       dum2,'IB ','QFF',.FALSE., ke_in,      0)
var_in(22) = ar_des_input('V         ','hybrid'           , iec, 109, 132,   0,   0,    1.0, 0.0, 3, &
                          dum4,    v_in,       dum2,'IB ','QFF',.FALSE., ke_in,      0)
var_in(24) = ar_des_input('T         ','hybrid'           , iec, 109, 130,   0,   0,    1.0, 0.0, 3, &
                          dum4,    t_in,       dum2,'IB ','QFF',.FALSE., ke_in,      0)
var_in(27) = ar_des_input('QV        ','hybrid'           , iec, 109, 133,   0,   0,    1.0, 0.0, 3, &
                          dum4,   qv_in,       dum2,'IB ','LFF',.FALSE., ke_in,      0)
var_in(28) = ar_des_input('QC        ','hybrid'           , iec, 109, 246,   0,   0,    1.0, 0.0, 3, &
                          dum4,   qc_in,       dum2,'IB ','LFF',.FALSE., ke_in,      0)
var_in(29) = ar_des_input('QI        ','hybrid'           , iec, 109, 247,   0,   0,    1.0, 0.0, 3, &
                          dum4,   qi_in,       dum2,'   ','LFF',.FALSE., ke_in,      0)
var_in(30) = ar_des_input('QR        ','hybrid'           , iec, 109,  75,   0,   0,    1.0, 0.0, 3, &
                          dum4,   qr_in,       dum2,'   ','LFF',.FALSE., ke_in,      0)
var_in(31) = ar_des_input('QS        ','hybrid'           , iec, 109,  76,   0,   0,    1.0, 0.0, 3, &
                          dum4,   qs_in,       dum2,'   ','LFF',.FALSE., ke_in,      0)
var_in(35) = ar_des_input('W_SO_REL  ','depthBelowLand'   , iec, 111,  39,   0,   0,    1.0, 0.0, 3, &
                          dum4, w_so_rel_in,   dum2,'I O','MFF',.FALSE., ke_soil_in+1,0)
var_in(41) = ar_des_input('PS        ','surface'          , iec,   1, 134,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,      ps_in,'IB ','LFF',.FALSE.,     1,      0)
var_in(42) = ar_des_input('LNPS      ','hybrid'           , iec, 109, 152,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,      ps_in,'IB ','LFF',.FALSE.,     1,      0)
var_in(43) = ar_des_input('T_SNOW    ','surface'          , iec,   1, 238,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,  t_snow_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(44) = ar_des_input('T_S       ','depthBelowLand'   ,   2, 111,  85,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,     t_s_in,'   ','MFF',.FALSE.,     1,      0)
var_in(45) = ar_des_input('SST       ','surface'          , iec,   1,  34,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,     t_s_in,'   ','MFF',.FALSE.,     1,      0)
var_in(46) = ar_des_input('T_SKIN    ','surface'          , iec,   1, 235,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,     t_s_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(49) = ar_des_input('T_G1      ','depthBelowLandLayer', iec, 112, 139, 0,   7,    1.0, 0.0, 2, &
                          dum4,    dum3,    t_g1_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(50) = ar_des_input('T_G2      ','depthBelowLandLayer', iec, 112, 170, 7,  28,    1.0, 0.0, 2, &
                          dum4,    dum3,    t_g2_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(51) = ar_des_input('T_G3      ','depthBelowLandLayer', iec, 112, 183,28, 100,    1.0, 0.0, 2, &
                          dum4,    dum3,    t_g3_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(52) = ar_des_input('T_CL      ','depthBelowLandLayer', iec, 112, 236,100,255,    1.0, 0.0, 2, &
                          dum4,    dum3,    t_cl_in,'I  ','MFF',.FALSE.,     1,      0)
var_in(53) = ar_des_input('W_SNOW    ','surface'          , iec,   1, 141,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,  w_snow_in,'IB ','MFT',.FALSE.,     1,      0)
var_in(54) = ar_des_input('W_I       ','surface'          , iec,   1, 198,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,     w_i_in,'I  ','MFF',.FALSE.,     1,      0)

IF ((izdate < 20000627 .AND. .NOT. lpost_0006) .OR. lante_0006) THEN
var_in(55) = ar_des_input('W_G1      ','depthBelowLandLayer'   , iec, 112, 140,   0,   7,    1.0, 0.0, 2, &
                          dum4,    dum3,    w_g1_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(56) = ar_des_input('W_G2      ','depthBelowLandLayer'   , iec, 112, 171,   7,  28,0.07/0.21,0.0,2, &
                          dum4,    dum3,    w_g2_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(57) = ar_des_input('W_G3      ','depthBelowLandLayer'   , iec, 112, 184,  28, 100,0.07/0.72,0.0,2, &
                          dum4,    dum3,    w_g3_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(58) = ar_des_input('W_CL      ','depthBelowLandLayer'   , iec, 112, 237, 100, 255,0.07/1.55,0.0,2, &
                          dum4,    dum3,    w_cl_in,'I  ','MFF',.FALSE.,     1,      0)
ELSE
var_in(55) = ar_des_input('W_G1      ','depthBelowLandLayer'   , iec, 112,  39,   0,   7,1.0/0.07,0.0, 2, &
                          dum4,    dum3,    w_g1_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(56) = ar_des_input('W_G2      ','depthBelowLandLayer'   , iec, 112,  40,   7,  28,1.0/0.21,0.0, 2, &
                          dum4,    dum3,    w_g2_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(57) = ar_des_input('W_G3      ','depthBelowLandLayer'   , iec, 112,  41,  28, 100,1.0/0.72,0.0, 2, &
                          dum4,    dum3,    w_g3_in,'IB ','MFF',.FALSE.,     1,      0)
var_in(58) = ar_des_input('W_CL      ','depthBelowLandLayer'   , iec, 112,  42, 100, 255,1.0/1.55,0.0, 2, &
                          dum4,    dum3,    w_cl_in,'I  ','MFF',.FALSE.,     1,      0)
ENDIF
var_in(59) = ar_des_input('QV_S      ','surface'          , iec,   1, 233,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,    qv_s_in,'IB ','MFT',.FALSE.,     1,      0)
var_in(62) = ar_des_input('FI        ','isobaricInhPa'    , iec, 100, 129,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,     fic_in,'IB ','LFF',.FALSE.,     1,      0)

IF (lcm2lm) THEN
var_in(46) = ar_des_input('T_SKIN    ','surface'          , iec,   1, 235,   0,   0,    1.0, 0.0, 2, &
                          dum4,    dum3,  t_skin_in,'IBO','MFF',.FALSE.,     1,      0)
ENDIF

! iso code
IF (liso) THEN
var_in(66) = ar_des_input('QV18O  ', 'hybrid',  206, 109, 160,   0,   0,    1.0, 0.0, 3,  &
                          dum4, riso_in(:,:,:,1), dum2,   'IB ', 'LFT', .FALSE., ke_in, 0)
var_in(67) = ar_des_input('QV2H   ', 'hybrid',  206, 109, 161,   0,   0,    1.0, 0.0, 3,  &
                          dum4, riso_in(:,:,:,2), dum2,   'IB ', 'LFT', .FALSE., ke_in, 0)
var_in(68) = ar_des_input('R18OSOIL ', 'surface', 206,   1, 162,   0,   0,    1.0, 0.0, 2,  &
                          dum4, dum3, risosoil_in(:,:,1), 'IB ', 'LFT', .FALSE.,     1, 0)
var_in(69) = ar_des_input('R2HSOIL  ', 'surface', 206,   1, 163,   0,   0,    1.0, 0.0, 2,  &
                          dum4, dum3, risosoil_in(:,:,2), 'IB ', 'LFT', .FALSE.,     1, 0)
ENDIF
! end iso code

!US Question: in case lec2lm, no grib-factors are set for w_snow, w_i? is this correct?
!             in case lcm2lm, the factors were set, but then I/O will be mainly
!                             by netcdf, which does not take factors
ENDIF

!------------------------------------------------------------------------------
! Section 3: Change settings depending on Namelist Input
!------------------------------------------------------------------------------

! set additional necessary external parameters
IF (lgme2lm .OR. llm2lm) THEN
  ! soiltype needed:
  var_in( 6)%dattyp = '  E'   ! SOILTYP
ENDIF

! get surface pressure from hydrostatic fields
IF (lgme2lm .OR. lec2lm .OR. lgfs2lm .OR. lgsm2lm .OR. lhir2lm) THEN
  var_in(41)%dattyp = 'IB '   ! PS
  var_in(62)%dattyp = 'IB '   ! FI control level
ELSE ! llm2lm, lum2lm
  var_in(41)%dattyp = '   '   ! PS
  var_in(62)%dattyp = '   '   ! FI control level
ENDIF

! pressure or pressure deviation for COSMO and UM
IF (llm2lm) THEN
  var_in(25)%dattyp = 'IB '   ! P
  var_in(26)%dattyp = 'IB '   ! PP
ENDIF

IF (lum2lm) THEN
  var_in(24)%ylevtyp = 'hybrid'   ! T
  var_in(25)%dattyp  = 'IB '   ! P
  var_in(25)%ipc     = 'LFF'   ! P  for linear interpolation
  var_in(27)%ylevtyp = 'hybrid'   ! QV
  var_in(29)%ylevtyp = 'hybrid'   ! QI

  ! UM has level type 120 in grib1
  var_in(21:32)%levtyp = 120

  ! for vertical interpolations we need also PS, T_S
  var_in(41)%dattyp = 'IB '   ! PS
  var_in(44)%dattyp = 'IB '   ! T_S
  var_in(46)%dattyp = 'IB '   ! T_SKIN
ENDIF

! adjust the dattyp of W, depending on lvertwind_ini/lvertwind_bd
IF (llm2lm) THEN
  IF (lvertwind_ini) THEN
    var_in(23)%dattyp(1:1) = 'I'   ! W
  ELSE
    var_in(23)%dattyp(1:1) = ' '
  ENDIF

  IF (lvertwind_bd) THEN
    var_in(23)%dattyp(2:2) = 'B'   ! W
  ELSE
    var_in(23)%dattyp(2:2) = ' '
  ENDIF
ENDIF

! adjust the dattyp of Qx, depending on lprog_qx
IF (lprog_qi) THEN
  var_in(29)%dattyp = 'IB '   ! QI
ELSE
  var_in(29)%dattyp = '   '
ENDIF

IF (lprog_qr_qs) THEN
  var_in(30)%dattyp = 'IB '   ! QR
  var_in(31)%dattyp = 'IB '   ! QS
ELSE
  var_in(30)%dattyp = '   '
  var_in(31)%dattyp = '   '
ENDIF

IF (lprog_qg) THEN
  var_in(32)%dattyp = 'IB '   ! QG
ELSE
  var_in(32)%dattyp = '   '
ENDIF

! adjust the dattyp of T_SO and W_SO and the fields for the old soil model,
! depending on lmulti_layer_in
! T_S is treated separate, because there are more dependencies
IF (lmulti_layer_in .AND. .NOT. lec2lm) THEN
  var_in(33)%dattyp = 'IB '       ! t_so
  var_in(34)%dattyp = 'I  '       ! w_so
  var_in(61)%dattyp = 'I  '       ! freshsnw
  IF (lprog_rho_snow) THEN
    var_in(63)%dattyp = 'I  '     ! rho_snow
  ENDIF
  var_in(48)%dattyp = '   '       ! t_m
  var_in(52)%dattyp = '   '       ! t_cl
  var_in(55)%dattyp = '   '       ! w_g1
  var_in(56)%dattyp = '   '       ! w_g2
  var_in(58)%dattyp = '   '       ! w_cl
! SP, 201111
  IF (lcm2lm .AND. (itype_t_cl == 0)) THEN
    var_in(33)%dattyp = '   '     ! t_so
    var_in(52)%dattyp = 'I  '     ! t_cl
  ENDIF
ELSE
  IF (lgme2lm .OR. llm2lm) THEN
    var_in(33)%dattyp = '   '     ! t_so
    var_in(34)%dattyp = '   '     ! w_so
    var_in(61)%dattyp = '   '     ! freshsnw
    var_in(63)%dattyp = '   '     ! rho_snow
    var_in(48)%dattyp = 'IB '     ! t_m
    var_in(52)%dattyp = 'I  '     ! t_cl
    var_in(55)%dattyp = 'IB '     ! w_g1
    var_in(56)%dattyp = 'IB '     ! w_g2
    IF (nl_soil_in == 3) THEN
      var_in(57)%dattyp = 'IB '   ! w_g3
      var_in(55)%levtop =   0.0   ! w_g1
      var_in(55)%levbot =   2.0   ! w_g1
      var_in(56)%levtop =   2.0   ! w_g2
      var_in(56)%levbot =  10.0   ! w_g2
    ENDIF
    var_in(58)%dattyp = 'I  '     ! w_cl
  ENDIF
ENDIF

! Treatment of T_S
IF (lgme2lm .OR. llm2lm) THEN
  IF (lmulti_layer_in) THEN
    IF (yin_form_read == 'ncdf') THEN
      var_in(44)%dattyp = 'IB '
    ELSE
      var_in(44)%dattyp = '   '
    ENDIF
  ELSE
    var_in(44)%dattyp = 'IB '
  ENDIF
ELSEIF (lum2lm) THEN
  var_in(44)%dattyp = 'IB '
ELSEIF (lgfs2lm) THEN
  var_in(44)%levtyp  = 1 ! surface
  var_in(44)%ylevtyp = 'surface'
  var_in(44)%dattyp  = 'IB '
ELSEIF (lcm2lm) THEN
  var_in(44)%dattyp = 'IBO'
ELSEIF (lhir2lm) THEN
  var_in(44)%dattyp  = 'IBO'
  var_in(44)%levtyp  = 1 ! surface
  var_in(44)%ylevtyp = 'surface'
ENDIF

! Treatment of QV_2M
IF (lhir2lm) THEN
  var_in(62)%dattyp  = 'IBO'
ENDIF

! adjust the dattyp of T_ICE and H_ICE, depending on lseaice
IF (lseaice) THEN
  var_in(64)%dattyp = 'I  '     ! t_ice
  var_in(65)%dattyp = 'I  '     ! h_ice
ENDIF

IF (lec2lm) THEN
  ! adjust the level for the control geopotential
  var_in(62)%levbot = pcontrol_fi / 100.0_ireals      ! fic_in
  IF (l_smi) THEN
    var_in( 6)%dattyp = '  E'     ! soiltyp
  ENDIF
ENDIF

! adjust the dattyp for lcm2lm
IF (lcm2lm) THEN
  var_in( 3)%dattyp = '  E'     ! hsurf
  var_in( 6)%dattyp = '  E'     ! soiltyp
  var_in(28)%dattyp = 'IBO'     ! qc
  var_in(29)%dattyp = 'IBO'     ! qi
  var_in(33)%dattyp = '   '     ! t_so
  var_in(34)%dattyp = '   '     ! w_so
  var_in(41)%dattyp = 'IB '     ! ps
  var_in(42)%dattyp = '   '     ! lnps
  var_in(43)%dattyp = 'IBO'     ! t_snow
  IF (luse_t_skin) var_in(46)%dattyp = 'IB '     ! t_skin
  var_in(53)%dattyp = 'IB '     ! w_snow
  var_in(59)%dattyp = 'IBO'     ! qv_s
  var_in(54)%dattyp = 'I O'     ! w_i
  var_in(45)%dattyp = 'IBO'     ! sst
  var_in(62)%dattyp = 'IBO'     ! fic_in
  var_in(44)%dattyp = 'IBO'     ! t_s
  var_in(49)%dattyp = '   '     ! t_g1   !_br 10.02.12
  var_in(50)%dattyp = '   '     ! t_g2   !_br 10.02.12
  var_in(51)%dattyp = '   '     ! t_g3   !_br 10.02.12
  var_in(57)%dattyp = '   '     ! w_g3   !_br 10.02.12
  var_in(61)%dattyp = 'I O'     ! freshsnw   !_br 10.02.12

  IF (itype_w_so_rel > 0) THEN
    var_in(35)%dattyp = 'I  '   ! w_so_rel
  ENDIF

  IF (lcm_hgt_coor) THEN
    var_in(25)%dattyp = 'IB '   ! P from UKMO, HadGEM2
    var_in(25)%ipc    = 'LFF'   ! P from UKMO, HadGEM2
  ENDIF
ENDIF

IF (lgfs2lm) THEN
  DO i = 21, 32
    var_in(i)%levtyp  = 100
    var_in(i)%ylevtyp = 'isobaricInhPa                 '
  ENDDO
  var_in(62)%dattyp  = 'IB '   ! FI
  var_in(62)%levtyp  = 100
  var_in(62)%ylevtyp = 'isobaricInhPa                 '
ENDIF

IF (lgsm2lm .OR. lhir2lm) THEN
  DO i = 21, 32
    var_in(i)%levtyp  = 109
    var_in(i)%ylevtyp = 'hybrid                        '
  ENDDO
  var_in(23)%levtyp  = 110
  var_in(23)%ylevtyp = 'hybridLayer                   '
ENDIF

END SUBROUTINE setup_vartab_grid_in

!==============================================================================
!==============================================================================

END MODULE src_gribtabs
