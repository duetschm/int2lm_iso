!+ Source Module for computine 2d fields
!==============================================================================

MODULE src_2d_fields

!==============================================================================
!
! Description:
!   This module contains routines for computing additional two-dimensional
!   fields for HM/LM. The computations include
!     - ozone distribution in HM/LM-domain
!     - ground temperatures
!     - specific and relative humidities at the surface
!     - water content of soil layers
!
! Current Code Owner: DWD, Ulrich Schaettler
!    phone:  +49  69  8062 2739
!    fax:    +49  69  8062 3721
!    email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2005/04/11 Ulrich Schaettler
!  Initial release for INT2LM
! 1.2        2005/07/22 Ulrich Schaettler
!  Use t_snow_lm instead of t_s_lm for consistency check of snow
! V1_5         2007/07/09 Ulrich Schaettler
!  Renamed czmls to czmls_lm;  ke_soil to ke_soil_lm
!  Changed computation of t_s_lm in case of lm2lm
! V1_6         2007/09/07 Ulrich Schaettler, Burkhardt Rockel, Uwe Boehm
!  Added lcm2lm and Subroutine init_multi_layer_cm
!  Transition from volumetric to relative input soil moisture vw_so_in to w_so_rel_in
!  Multiplication of w_so_lm by pore volume is skipped in SR ground_fields for lcm2lm
!  Modified calculation of initial soil moisture in case of lcm2lm in SR init_multi_layer_cm
! V1_7         2007/11/26 Ulrich Schaettler, Uwe Boehm
!  Renamed iw_so_rel_type to itype_w_so_rel to be consistent with other names
!  Adaptation in SR init_multi_layer_cm for itype_w_so_rel
!  Added call to SR init_multi_layer_cm if artificial soil moisture profile is used
!  Added SR plant_characteristics to compute actual values for slowly varying
!   external parameters for initial and boundary (climate mode) data
!  Added option to use actual ndvi ratio (interpolated from GME monthly values)
! V1_8         2008/05/29 Ulrich Schaettler, Hans-Juergen Panitz, Burkhardt Rockel
!  Corrected computation of rootdp_lm for options itype_rootdp = 1/2
!  Added another option itype_rootdp = 3
!  Introduced debug output
! V1_9         2009/09/03 Ulrich Schaettler, et al.
!  Eliminated GME fields for ndvi-values
!  Implemented SR month2hour to compute the months and the weights to 
!  interpolate monthly values to a certain day of the year
!  Implemented further options for itype_ndvi and itype_aerosol
!   (computation from montly values or ratio of monthly mean values)
!  Implemented 2 more options for itype_w_so_rel (=3,4; by Uwe Boehm)
! V1_10        2009/12/17 Heike Vogel
!  Set freshsnw to 0 only for lmulti_layer_lm
!  Initialize prognostic FLake variables for cold start
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_12        2010/06/14 Burkhardt Rockel, Uwe Boehm (CLM)
!  Modifications for itype_w_so_rel
! V1_14        2010/11/19 Ulrich Schaettler, Burkhardt Rockel
!  Modifications to allow processing of JMA and NCEP data
!  Extension of an IF-clause for the FLake Model (BR)
! V1_19        2012/06/06 Ulrich Schaettler, Dmitrii Mironov, Susanne Brienen
!                         Burkhardt Rockel, Juergen Helmert
!  Added lhir2lm as internal logical flag
!  Changes to initialization of FLake variables in case of cold-/warm-start (DM)
!  New SR inti_multi_layer_gme_ml for vertical interpolation of GME multi-layer
!   soil variables to a different number of levels (SB)
!  Correction in init_multi_layer_cm:
!   - include W_SO for deepest layer (no longer == 0)
!   - correct multiplication of layer thickness
!   - correct calculation of water content in soil layers deeper than the
!     global model soil layers in case of itype_w_so_rel=2
!  Added treatment of 12-monthly surface albedo fields in case itype_albedo==3
! V1_20        2012/09/03 Ulrich Schaettler, Burkhardt Rockel
!  Enlarged strings for date variables to 14 characters
!  For the correction of snow variables a dependency on the latitude was added
!    Now it is only done for lat > -60.0 (Burkhardt Rockel)
! V1_21        2013/03/25 Burkhardt Rockel
!  t_so(layer 0) has to be set to t_s_lm in SR init_multi_layer_cm
!   (needed for GRIB format)
! V1_22        2013/07/11 Ulrich Schaettler, Davide Cesari, Ulrich Blahak
!  Eliminated use of aklm, bklm in SR ground_fields (which was unnecessary)
!  Allow ke_soil_lm /= ke_soil_in for IFS input (Davide)
!  Corrected interpolation weights when interpolating soil moisture from input
!    soil layers to COSMO soil layers in subroutine init_multi_layer_gme_ml (UB)
!
!
! Added height correction of isotopes of soil. Hui Tang 2013-11-20
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters , ONLY :  &
  ireals,    & ! KIND-type parameters for real variables
  iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_fields_lm, ONLY : &
  fis_lm    ,      & ! orography * G                                    (m2/s2)
  hsurf_lm  ,      & ! orography                                        (  m  )
  soiltyp_lm,      & ! type of the soil (keys 0-9)                      (  1  )
  fr_land_lm,      & ! land fraction of grid element                    (  1  )
  plcov_lm     ,   & ! fraction covered by plants                       (  1  )
  plcov_mx_lm  ,   & ! plant cover during vegetation time               (  1  )
  plcov_mn_lm  ,   & ! plant cover during time of rest                  (  1  )
  lai_mx_lm    ,   & ! leaf area index during vegetation time           (  1  )
  lai_mn_lm    ,   & ! leaf area index during time of rest              (  1  )
  lai_lm       ,   & ! leaf area index                                  (  1  )
  rootdp_lm    ,   & ! depth of the roots                               (  m  )
  rootdp_mx    ,   & ! depth of the roots from external parameters      (  m  )
  ps_lm     ,      & ! surface pressure                                 ( Pa  )
  t_s_lm    ,      & ! temperature of the ground surface                (  K  )
  t_snow_lm ,      & ! temperature of the snow surface                  (  K  )
  t_m_lm    ,      & ! temperature between upper and medium soil layer  (  K  )
  t_cl_lm   ,      & ! temperature between medium and lower soil layer  (  K  )
  qv_s_lm   ,      & ! specific water vapor content on the surface      (kg/kg)
  w_snow_lm ,      & ! water content of the snow                        (m H2O)
  w_g1_lm   ,      & ! water content of the upper soil layer            (m H2O)
  w_g2_lm   ,      & ! water content of the medium soil layer           (m H2O)
  w_g3_lm   ,      & ! water content of the lower soil layer            (m H2O)
                     ! (if nl_soil_lm = 3, unused otherwise)
  w_cl_lm   ,      & ! climatological deep soil water content           (m H2O)
  t_so_lm   ,      & ! multi-layer soil temperature                     (  K  )
  dt_so_gl  ,      & ! multi-layer soil temperature diff. (down-up for interp.)
  w_so_lm   ,      & ! multi-layer soil moisture                        (m H2O)
  w_so_gl   ,      & ! multi-layer soil moisture  (for interpolation)   (m H2O)
  freshsnw_lm,     & ! weighting function indicating 'freshness' of snow
  hmo3_lm   ,      & ! height of maximum ozone concentration            ( Pa  )
  vio3_lm            ! total vertically integrated ozone content        (Pa O3)

USE data_fields_lm, ONLY : &
  latlm_m   ,      & ! latitudes of the LM grid points
  lonlm_m   ,      & ! longitudes of the LM grid points
  baryll_m  ,      & !
  w_intpol  ,      & ! interpolation weights for gme2lm                 (  -  )
  n_intpol  ,      & ! nearest GME gridpoint for nearest neighbor interpolation
  m_intpol  ,      & ! nearest GME gridpoint with same lsm for match interpolation
  l_intpol  ,      & ! to use a far away GME gridpoint with same lsm
  index_m   ,      & !
  ndviratio_lm,    & ! actual ndvi ratio                                (  1  )
  ndvi_mrat_lm,    & ! actual ndvi ratio                                (  1  )
  alb_dif_lm  ,    & ! solar surface albedo - diffuse                   (  1  )
  alb_dif12_lm,    & ! solar surface albedo - diffuse                   (  1  )
  aer_su12_lm,     & ! Tegen (1997) aerosol type sulfate drops          (  -  )
  aer_du12_lm,     & ! Tegen (1997) aerosol type mineral dust coarse    (  -  )
  aer_or12_lm,     & ! Tegen (1997) aerosol type organic(water solub)   (  -  )
  aer_bc12_lm,     & ! Tegen (1997) aerosol type black carbon           (  -  )
  aer_ss12_lm,     & ! Tegen (1997) aerosol type sea salt               (  -  )
  aer_su_lm ,      & ! Tegen (1997) aerosol type sulfate drops          (  -  )
  aer_du_lm ,      & ! Tegen (1997) aerosol type mineral dust coarse    (  -  )
  aer_or_lm ,      & ! Tegen (1997) aerosol type organic(water solub)   (  -  )
  aer_bc_lm ,      & ! Tegen (1997) aerosol type black carbon           (  -  )
  aer_ss_lm ,      & ! Tegen (1997) aerosol type sea salt               (  -  )
  ps_gl     ,      & ! surface pressure on the interpol. GME orogr.     ( Pa  )
  fis_gl    ,      & ! GME interpolated orography * G                   (m2/s2)
  t_s_gl    ,      & ! temperature of the ground surface                (  K  )
  rh_s_gl   ,      & ! relative humidity at the surface                 (kg/kg)
  dtms_gl   ,      & ! t_m_lm    - t_s_lm                               (  K  )
  dtkes_gl  ,      & ! t(ke)_lm  - t_s_lm                               (  K  )
  dtssnow_gl,      & ! t_s_lm    - t_snow_lm                            (  K  )
  t_lm      ,      & ! temperature                                      (  K  )
  hhl_gl    ,      & ! height of half-levels on the interpol. COARSE LM oro.(m)
  hhl_lm    ,      & ! height of half-levels of LM                      (  m  )
  i_index,         & ! i-index of coarse mesh grid point which is lower left
                     ! to a given LM (fine mesh) grid point
  j_index,         & ! j-index of coarse mesh grid point which is lower left
                     ! to a given LM (fine mesh) grid point
  x_wght,          & ! relative distance between x- (i-) coarse mesh and
                     ! fine mesh grid points
  y_wght,          & ! relative distance between y- (j-) coarse mesh and
                     ! fine mesh grid points
  lmask_lm  ,      & ! mask of points on the frame
  lolp_lm,         & ! Land Sea Mask of LM for 'M'atch Interpolation
! iso code  Hui Tang 2013-11-20
  riso_lm,         & ! isotope ratios in water vapor;
  risosoil_lm,     & ! isotope ratios in soil moisture;
  drisoke_gl         ! riso_lm(:,:,ke,1-2) - risosoil_lm (:,:,1-2) (1: 18O; 2: 2H);
! end iso code


USE data_fields_lm, ONLY : &
  z0_lm        ,   & ! effectively used roughness length in CCLM
  plcov12      ,   & ! monthly climatology of fractional plant cover
  z012         ,   & ! monthly climatology of rooting depth
  lai12        ,   & ! monthly climatology of leaf area index
  fr_lake_lm   ,   & ! lake fraction of grid element                 (  1  ) !_br 13.08.10
  depth_lk_lm  ,   & ! lake depth                                    (  m  )
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

USE data_fields_in, ONLY : &
  lat_coarse_m,   & ! latitudes of the LM grid points
  lon_coarse_m,   & ! longitudes of the LM grid points
  w_so_rel_in,    & ! multi-layer relative   soil moisture      (  1  )
  lolp_gme,       & ! Land Sea Mask of GME for 'M'atch Interpolation
  lolp_in           ! Land Sea Mask of input fields for 'M'atch Interpolation

!------------------------------------------------------------------------------

USE data_grid_lm,   ONLY : &
  ie2lm,       & !
  je2lm,       & !
  kelm,        & !
  kedim,       & !
  ke_soil_lm,  & ! number of levels in multi-layer soil model in output
  cw_so_rel_lm,& ! artificial volumetric soil water content profile      !_br
  czmls_lm,    & ! depth of the main soil layers in meters in output
  czhls_lm       ! depth of the half soil layers in meters in output

!------------------------------------------------------------------------------

USE data_grid_in,   ONLY : &
  ak_in ,      & ! vertical coordinate parameters for half levels
  bk_in ,      & !                  - " -
  ie_in,       & ! ie for input grid, local domain
  je_in,       & ! je for input grid, local domain
  ie_in_tot,   & ! ie for input grid, total domain
  je_in_tot,   & ! je for input grid, total domain
  startlat_in, & ! transformed latitude of the lower left grid point
                 ! of the local domain (in degrees, N>0)
  startlon_in, & ! transformed longitude of the lower left grid point
                 ! of the local domain (in degrees, E>0)
  ke_in ,      & ! number of levels in the vertical
  ke_soil_in,  & ! number of input levels in multi-layer soil model
  czmls_in,    & ! depth of the input soil layers in meters
  czhls_in,    & ! depth of the input half soil layers in meters
  latitudes_in,& ! latitudes of the input data
  longitudes_in,&! longitudes of the input data
  grdpt_rel_in,& ! relation between input and lm grid for cressman scheme
  jd_min,      & ! smallest index for diamonds for a LM subdomain
  jd_max,      & ! biggest index for diamonds for a LM subdomain
  igg1s,       & ! start index of global array-dimension in x-direction
  igg1sm2,     & ! = igg1s - 2
  igg1e,       & ! end index of global array-dimension in x-direction
  igg1ep2,     & ! = igg1e + 2
  igg2s,       & ! start index of global array-dimension in y-direction
  igg2sm2,     & ! = igg2s - 2
  igg2e,       & ! end index of global array-dimension in y-direction
  igg2ep2,     & ! = igg2e + 2
  ispoke         ! offsets of the 6 (5) neighbouring gridpoints relative to
                 ! i1-direction use ispoke(m), m=1,6 (5), in i2-direction
                 ! use ispoke(m+6), m=1,6 (5); phys. dim. ( - )

!------------------------------------------------------------------------------

USE data_int2lm_io,        ONLY : &
  nvar_lm, var_lm, undef         !

!------------------------------------------------------------------------------

USE data_int2lm_constants, ONLY : &
    G,       & ! gravity at sea level                          [ms-2]
    T0,      & ! 0 degree Celsius                              [Kelvin]
    Rdv,     & ! = R_d/R_v,
    O_m_rdv, & ! = 1. - Rdv
    Pi,      & ! circle constant
    degrad,  & !
    B1,      & !  a
    B2_w,    & !  b
    B3,      & !  c/b (0 degree Celsius [Kelvin])
    B4_w,    & !  d
    porb,    & ! poren volume of ground types
    pwpb,    & ! permanent wilting point of COSMO soil types
    fcb        ! field capacity of COSMO soil types

!------------------------------------------------------------------------------

USE data_int2lm_control,       ONLY :  &
    lgme2lm,         & ! if .TRUE., gme->lm
    lgfs2lm,         & ! if .TRUE., gfs->lm
    lgsm2lm,         & ! if .TRUE., gsm->lm
    lec2lm,          & ! if .TRUE., ec ->lm
    lhir2lm,         & ! if .TRUE., hirlam ->lm
    llm2lm,          & ! if .TRUE., lm ->lm
    lcm2lm,          & ! if .TRUE., cm ->lm    !_br
    lbdclim,         & ! if .TRUE., special boundary data for climate mode
    l_cressman,      & ! logical switch for controling the use of a Cressman scheme
    l_smi,           & ! if .TRUE., interpolate soil moisture with FC-PWP index
    lseaice,         & ! if .TRUE., run with sea ice mode
    llake_coldstart, & ! if .TRUE., initialize prognostic lake variables for cold start
    itype_w_so_rel,  & ! type of relative soil moisture input (0,1,2)
    itype_rootdp,    & ! to choose treatment of root depth
    itype_ndvi,      & ! to choose treatment of surface parameters (plcov, lai)
    itype_t_cl,      & ! to choose origin and treatment of deep soil temperature
    itype_aerosol,   & ! to choose treatment of surface parameters (plcov, lai)
    itype_albedo,    & ! to choose treatment of surface parameters (albedo)
    lt_cl_corr,      & ! if .TRUE., a height-correction of T_CL is performed
    lbd_frame_cur,   & ! if .TRUE., current boundary fields include only frames
    linitial,        & ! if .TRUE., initial data for LM
    lcomp_bound,     & ! compute fields for boundaries
    lmulti_layer_lm, & ! if .TRUE., compute soil fields for multi-layer soil model
    lmulti_layer_in, & ! if .TRUE., incoming data from new multi-layer soil model
    nl_soil_lm,      & ! number of soil layers in LM, resp. HM
    idbg_level,      & ! to control the verbosity of debug output
    lprintdeb_all      ! .TRUE.:  all tasks print debug output
                       ! .FALSE.: only task 0 prints debug output

!------------------------------------------------------------------------------

USE data_int2lm_parallel,    ONLY:   my_cart_id

!------------------------------------------------------------------------------

USE meteo_utilities,  ONLY:   psat_w, tgcom
USE interp_utilities, ONLY:   interp_l
USE gme_utilities,    ONLY:   pp_interp2ls, get_ndvi

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================
!+ organizes the computation of additional 2D fields
!------------------------------------------------------------------------------

SUBROUTINE org_2d_fields (ydate, nactday, acthour)

!------------------------------------------------------------------------------
!
! Description:
!   org_2d_fields organizes the computation of additional two-dimensional 
!   fields for LM/HM. These include
!     - ozone distribution in LM/HM-domain
!     - ground temperatures
!     - specific and relative humidities at the surface
!     - water content of soil layers
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
CHARACTER (LEN=14), INTENT(IN)    ::  &
  ydate           ! actual date in the form   yyyymmddhhmmss

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  nactday         ! actual day in the year

REAL (KIND=ireals), INTENT(IN)    ::  &
  acthour         ! hour of the day

!------------------------------------------------------------------------------

! Local arrays
REAL (KIND=ireals)         ::  &
  zwei1, zwei2,                        &  ! for monthly interpolation
  tpl_T_r, tpl_T_f, C_T_min, rflk_depth_bs_ref  ! for FLake cold start

INTEGER  (KIND=iintegers)  ::  &
  imo1, imo2, i,j,      & ! for monthly interpolation
  izdebug,              & ! for debug print outs
  izerror,              & ! status and error status variable
  kezvert                 ! to choose the vertical dimension

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling

!
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  izerror   = 0
  yzerrmsg  = '  '
  yzroutine = 'org_2d_fields'

  ! Choose, whether HM- or LM-fields shall be computed
  IF (lgme2lm .OR. lgfs2lm .OR. lgsm2lm .OR. lec2lm .OR. lcm2lm .OR. lhir2lm) THEN
    kezvert = ke_in
  ELSEIF (llm2lm) THEN
    ! atmospheric fields are already interpolated to LM output levels
    kezvert = kelm
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Compute ozone distribution
!------------------------------------------------------------------------------

  ! in NWP mode, the ozone distributions are only computed for initial fields
  ! in climate mode they are computed in every step

  IF (.NOT. lcomp_bound .OR. lbdclim) THEN
    IF (izdebug > 5) THEN
      PRINT *, '   Compute ozone distributions: vio3, hmo3'
    ENDIF

    CALL ozone_compute (ydate, nactday, acthour)
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Compute ground fields
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, '   Compute ground fields'
  ENDIF

  CALL ground_fields (kezvert)

!------------------------------------------------------------------------------
! Section 4: Compute values for multi-layer soil model (if wanted)
!------------------------------------------------------------------------------

  IF (.NOT. lcomp_bound) THEN

    IF (lcm2lm) THEN

      ! When processing data from climate models, T_SO and W_SO are always
      ! initialized here. W_SO depending on itype_w_so_rel and T_SO is done
      ! with a "weighted interpolation" between t_s_lm and t_cl_lm onto the
      ! used vertical levels.

      IF (izdebug > 5) THEN
        PRINT *, '   Compute values for (multi-layer) soil model'
      ENDIF
      CALL init_multi_layer_cm

    ELSEIF (lgme2lm .AND. lmulti_layer_in .AND. (ke_soil_in /= ke_soil_lm) ) THEN

      ! new case for interpolating to a different number of output layers
      IF (izdebug > 5) THEN
        PRINT *, '   Compute values for multi-layer soil model from GME multi-layer model'
        PRINT *, '   (vertical interpolation to a different number of layers)'
      ENDIF
      CALL init_multi_layer_gme_ml

    ELSE   ! all other int2lms

      IF (lmulti_layer_lm .AND. lmulti_layer_in .AND. itype_w_so_rel == 0) THEN
          CALL init_multi_layer_cm
      END IF

      IF (lmulti_layer_lm .AND. (.NOT. lmulti_layer_in) ) THEN
        IF (itype_w_so_rel == 0) THEN
          ! artificial profile for w_so as when using a climate model 
          ! also t_so is set in this routine: for that we need the values 
          ! of t_cl_lm!
          IF (izdebug > 5) THEN
            PRINT *, '   Compute artificial values for (multi-layer) soil model'
          ENDIF
          CALL init_multi_layer_cm
        ELSE
          ! the multi layer soil variables for the LM have to be splitted from
          ! the 2 (or 3) layers from the incoming data
          ! (only works for GME2LM and LM2LM !!!)
          IF (izdebug > 5) THEN
            PRINT *, '   Compute values for multi-layer soil model from old 2-layer model'
          ENDIF
          CALL init_multi_layer
        ENDIF
      ENDIF

    ENDIF

    ! Initialize FLake variables for cold start.
    ! Now also values for t_ice and h_ice are set:
    !   - over lake grid boxes they are set to some reference values
    !     (otherwise errors have been reported during COSMO-Model runs)
    !   - over ocean boxes they are only set to reference values, 
    !     if seaice-module is switched off (lseaice=.FALSE.)
    !     otherwise values of t_ice and h_ice are interpolated from 
    !     the coarse grid (or are taken from the analysis).

    IF (llake_coldstart) THEN
      tpl_T_r           = 277.13_ireals     ! Temperature of maximum density [K]
      tpl_T_f           = 273.15_ireals     ! Fresh water freezing point     [K]
      C_T_min           = 0.5_ireals        ! Minimum value of C_T_LK        [-]
      rflk_depth_bs_ref = 10.0_ireals       ! Reference value                [m]

      DO j = 1, je2lm
        DO i = 1, ie2lm

          ! For coldstarts, the bottom-sediment module is switched off everywhere
          t_b1_lk_lm (i,j) = tpl_T_r
          h_b1_lk_lm (i,j) = rflk_depth_bs_ref

!US       IF ( (depth_lk_lm(i,j) < 0.0_ireals) .OR. (fr_lake_lm(i,j) <= 0.5_ireals) ) THEN
!US           this can go wrong for tile approach. Therefore we must not ask for fr_lake_lm

          IF   (depth_lk_lm(i,j) < 0.0_ireals) THEN

            ! Ocean or land point, set FLake variables to their reference values
            ! (only checking depth_lk_lm < 0 is not sufficient, because there 
            !  are still some external parameter data sets around where lake 
            !  depths are given, if fr_lake_lm < 0.5 and the grid box is 
            !  considered as land)
            ! At DWD there now is the convention, that depth_lk = -1, 
            !  if and only if fr_lake <= 0.5
            ! In external_data, depth_lk is set according to this convention,
            !  if it is not available as external data set.

            t_wml_lk_lm(i,j) = tpl_T_r
            t_mnw_lk_lm(i,j) = tpl_T_r
            t_bot_lk_lm(i,j) = tpl_T_r
            c_t_lk_lm  (i,j) = C_T_min
            h_ml_lk_lm (i,j) = 0._ireals

            IF (.NOT. lseaice) THEN
              ! If sea ice parameterisation scheme is not used,
              ! set t_ice and h_ice to their reference values;
              ! otherwise do not touch t_ice and h_ice (for sea and land grid boxes),
              ! assuming that their values are interpolated from the coarse grid
              ! or are taken from the analysis.
              t_ice_lm   (i,j) = tpl_T_f
              h_ice_lm   (i,j) = 0._ireals
            ENDIF

          ELSE                    

            ! Lake point, initialize FLake variables 

            ! Reference value (if less than the lake depth)
            h_ml_lk_lm (i,j) = MIN(depth_lk_lm(i,j), 10._ireals)

            ! Water surface temperature from coarse grid (but no ice at cold start,
            ! hence t_s should be non-negative)
            t_s_lm     (i,j) = MAX(t_s_lm(i,j), tpl_T_f)       

            IF(lmulti_layer_lm) THEN
              ! Save t_s for multi-layer soil model
              ! Setting t_so_lm(:,:,0)=t_s_lm(:,:) is absolutely necessary here
              ! as it is t_so_lm, but not t_s_lm, that is written to output
              t_so_lm(i,j,0) = t_s_lm(i,j)
            ENDIF

            t_wml_lk_lm(i,j) = t_s_lm(i,j)         ! Set mixed-layer temperature equal to t_s
            c_t_lk_lm  (i,j) = C_T_min             ! Minimum value
            t_bot_lk_lm(i,j) = tpl_T_r             ! Temperature of maximum density,
            IF (h_ml_lk_lm (i,j) >= depth_lk_lm(i,j)) THEN
              ! or mixed-layer temperature
              t_bot_lk_lm(i,j) = t_wml_lk_lm(i,j)
            ENDIF
            t_mnw_lk_lm(i,j) = t_wml_lk_lm(i,j) - c_t_lk_lm(i,j)                &
                 * MAX(0._ireals, (1._ireals-h_ml_lk_lm(i,j)/depth_lk_lm(i,j))) &
                 * (t_wml_lk_lm(i,j)-t_bot_lk_lm(i,j))

            ! At cold start, t_ice and h_ice must be set to their reference values
            ! for lake grid boxes (otherwise errors may occur during COSMO-model run)
            t_ice_lm   (i,j) = tpl_T_f             ! Assume no ice
            h_ice_lm   (i,j) = 0._ireals           ! Assume no ice

          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 5: Compute external parameters
!------------------------------------------------------------------------------

  IF (linitial .OR. lbdclim) THEN
    IF (izdebug > 5) THEN
      PRINT *, '   Compute plant characteristics'
    ENDIF

    CALL plant_characteristics (nactday, ydate, izdebug)

    IF (itype_aerosol == 2) THEN
      IF (izdebug > 5) THEN
        PRINT *, '   Compute aerosol values'
      ENDIF

      CALL month2hour (ydate, imo1, imo2, zwei1, zwei2, izdebug)

      aer_su_lm(:,:) = aer_su12_lm(:,:,imo1) +                    &
                      (aer_su12_lm(:,:,imo2) - aer_su12_lm(:,:,imo1)) * zwei2
      aer_du_lm(:,:) = aer_du12_lm(:,:,imo1) +                    &
                      (aer_du12_lm(:,:,imo2) - aer_du12_lm(:,:,imo1)) * zwei2
      aer_or_lm(:,:) = aer_or12_lm(:,:,imo1) +                    &
                      (aer_or12_lm(:,:,imo2) - aer_or12_lm(:,:,imo1)) * zwei2
      aer_bc_lm(:,:) = aer_bc12_lm(:,:,imo1) +                    &
                      (aer_bc12_lm(:,:,imo2) - aer_bc12_lm(:,:,imo1)) * zwei2
      aer_ss_lm(:,:) = aer_ss12_lm(:,:,imo1) +                    &
                      (aer_ss12_lm(:,:,imo2) - aer_ss12_lm(:,:,imo1)) * zwei2
    ENDIF

    IF (itype_albedo == 3) THEN
      IF (izdebug > 5) THEN
        PRINT *, '   Compute albedo values'
      ENDIF

      CALL month2hour (ydate, imo1, imo2, zwei1, zwei2, izdebug)

      alb_dif_lm(:,:) = alb_dif12_lm(:,:,imo1) +                    &
                       (alb_dif12_lm(:,:,imo2) - alb_dif12_lm(:,:,imo1)) * zwei2
    ENDIF

  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE org_2d_fields

!==============================================================================
!+ Computes the distribution of the ozone.
!------------------------------------------------------------------------------

SUBROUTINE ozone_compute (ydate, nactday, acthour)

!------------------------------------------------------------------------------
!
! Description:
!   Computes the ozone distribution in the LM/HM domain with annual cycle.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
CHARACTER (LEN=14), INTENT(IN)    ::  &
  ydate           ! actual date in the form   yyyymmddhhmmss

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  nactday         ! actual day in the year

REAL (KIND=ireals), INTENT(IN)    ::  &
  acthour         ! hour of the day

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers) ::   &
  i, j, jj, imm, imnc, imns, jmm, jnn

REAL (KIND=ireals)       ::   &
  ztwo,ztho,zcos1,zcos2,zcos3,zcos4,zcos5,zsin1,zsin2,zsin3,zsin4,zsin5

! Local arrays:
REAL (KIND=ireals)       ::   &
  zcqc (21), zchc(21), zcqs(15), zchs(15), zalp(66), zfq(11), zfh(11)

!- End of header
!------------------------------------------------------------------------------

  ! changed from 2 digits (jj-49) to 4 digits (jj-1949)   !?!?!?!?
  READ(ydate(1:4), '(I4)') jj
  ztwo = 0.681 + 0.2422 * (jj-1949) - (jj-1949) * 0.25
  ztho = 2.0*Pi*(nactday + acthour/24.0 - 1 + ztwo)/365.2422

  ! Call *ozone* to obtain the Fourier coefficients
  ! of the ozone content and the height of the maximum
  CALL ozone(ztho, zcqc, zcqs, zchc, zchs)

  DO j = 1, je2lm
    DO i = 1, ie2lm
      ! Call *legtri* to compute the Legendre coefficients
      ! separately for each gridpoint of LM/HM
      CALL legtri (SIN(latlm_m(i,j)*degrad),6,zalp)

      zfq(:) = 0.0
      zfh(:) = 0.0
      imm  = 0
      imnc = 0
      imns = 0
      DO jmm = 1, 6
        imm = imm+1
        DO jnn = jmm, 6
          imnc     = imnc+1
          zfq(imm) = zfq(imm) + zalp(imnc) * zcqc(imnc)
          zfh(imm) = zfh(imm) + zalp(imnc) * zchc(imnc)
        ENDDO
        IF (jmm /= 1) THEN
          imm = imm+1
          DO jnn = jmm, 6
            imns     = imns+1
            zfq(imm) = zfq(imm) + zalp(imns+6) * zcqs(imns)
            zfh(imm) = zfh(imm) + zalp(imns+6) * zchs(imns)
          ENDDO
        ENDIF
      ENDDO ! jmm

      ! inverse Fourier-Transformation.
      zcos1 = COS(lonlm_m(i,j)*degrad)
      zsin1 = SIN(lonlm_m(i,j)*degrad)
      zcos2 = zcos1*zcos1-zsin1*zsin1
      zsin2 = zsin1*zcos1+zcos1*zsin1
      zcos3 = zcos2*zcos1-zsin2*zsin1
      zsin3 = zsin2*zcos1+zcos2*zsin1
      zcos4 = zcos3*zcos1-zsin3*zsin1
      zsin4 = zsin3*zcos1+zcos3*zsin1
      zcos5 = zcos4*zcos1-zsin4*zsin1
      zsin5 = zsin4*zcos1+zcos4*zsin1
      vio3_lm(i,j) = zfq(1)+  2.*(zfq(2)*zcos1+zfq(3)*zsin1+zfq(4)*zcos2 + &
                     zfq(5)*zsin2+zfq(6)*zcos3+zfq(7)*zsin3+zfq(8)*zcos4 + &
                     zfq(9)*zsin4+zfq(10)*zcos5+zfq(11)*zsin5)
      hmo3_lm(i,j) = zfh(1)+  2.*(zfh(2)*zcos1+zfh(3)*zsin1+zfh(4)*zcos2 + &
                     zfh(5)*zsin2+zfh(6)*zcos3+zfh(7)*zsin3+zfh(8)*zcos4 + &
                     zfh(9)*zsin4+zfh(10)*zcos5+zfh(11)*zsin5)
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE ozone_compute

!==============================================================================
!+ For the yearly cycle of the ozone distibution.
!------------------------------------------------------------------------------

SUBROUTINE ozone (pytime, pqc, pqs, phc, phs)
 
!------------------------------------------------------------------------------
!
! Description:
!   This routine computes instantaneous values of a T5 spectral
!   distribution for two ozone parameters (total quantity and pressure
!   at the maximum of concentration, both in pascal) from the time
!   of the year.
!
! Method:
!  Straightforward, a second order Fourier development for the time of
!  the year.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
REAL (KIND=ireals), INTENT(IN)  ::   &
  pytime

! Array arguments with intent(out):
REAL (KIND=ireals), INTENT(OUT) ::   &
  pqc(21), pqs(15), phc(21), phs(15)

!------------------------------------------------------------------------------
!
! Local scalars:
REAL    (KIND=ireals)           ::   &
  zytime, zc1yt, zs1yt, zc2yt, zs2yt

INTEGER (KIND=iintegers)        ::   &
  jmn

! Local arrays:
REAL    (KIND=ireals)           ::                                      &
          zqc0(21),zqc1(21),zqc2(21),zqc3(21),zqc4(21),                 &
          zqs0(15),zqs1(15),zqs2(15),zqs3(15),zqs4(15),                 &
          zhc0(21),zhc1(21),zhc2(21),zhc3(21),zhc4(21),                 &
          zhs0(15),zhs1(15),zhs2(15),zhs3(15),zhs4(15)
!
!- End of header
!------------------------------------------------------------------------------

!
!*    DATA STATEMENTS.
!     ---- -----------
!
!          *ZOZ Q/H C/S N* (N=0,4) CORRESPONDS TO THE *POZ Q/H C/S*
!     (SEE ABOVE) AND TO THE FIVE TERMS OF THE *FOURIER DEVELOPMENT.
!
DATA zqc0/                                                              &
     &+.6012e-01,+.1887e-02,+.7410e-02,+.9950e-03,-.1426e-02,-.2072e-03,&
     &           -.4954e-03,+.7955e-05,-.3701e-03,+.4116e-04,-.4163e-04,&
     &                      -.2933e-03,+.2154e-04,-.2849e-03,-.1604e-03,&
     &                                 -.1054e-03,+.4974e-03,+.1047e-03,&
     &                                            +.8323e-04,+.2874e-03,&
     &                                                       +.1333e-03/
DATA zqs0/                                                              &
     &           +.4210e-03,-.9591e-03,+.2811e-03,-.2257e-03,-.1713e-03,&
     &                      -.3538e-03,+.1095e-03,-.4390e-03,-.5605e-05,&
     &                                 +.1478e-03,+.2849e-03,+.3430e-03,&
     &                                            +.8248e-04,+.1442e-03,&
     &                                                       -.1375e-04/
DATA zhc0/                                                              &
     &+.3166e+04,+.8663e+02,+.9401e+03,+.1999e+02,-.3530e+03,-.3311e+02,&
     &           -.4903e+02,-.4015e+00,-.1333e+02,+.5675e+01,+.7221e+01,&
     &                      -.3001e+02,+.7570e+01,-.1142e+02,-.1365e+02,&
     &                                 -.1502e+02,+.4911e+02,+.1425e+02,&
     &                                            +.8983e+01,+.3064e+02,&
     &                                                       +.1693e+02/
DATA zhs0/                                                              &
     &           +.4231e+02,-.7391e+02,+.1273e+02,+.2086e+02,-.1597e+02,&
     &                      -.3591e+02,+.1059e+02,-.2779e+02,-.6923e+01,&
     &                                 +.1397e+02,+.2387e+02,+.2883e+02,&
     &                                            +.8626e+01,+.1607e+02,&
     &                                                       -.2676e+01/
DATA zqc1/                                                              &
     &+.7090e-04,+.4930e-05,+.6829e-03,+.1897e-03,+.7226e-04,-.2807e-03,&
     &           +.4970e-04,-.1753e-03,-.7843e-04,-.1649e-03,-.1037e-03,&
     &                      -.4830e-04,-.6304e-04,-.1100e-03,-.7952e-04,&
     &                                 +.1326e-04,+.2599e-04,+.9926e-05,&
     &                                            -.9247e-05,-.3521e-05,&
     &                                                       -.1780e-04/
DATA zqs1/                                                              &
     &           +.6333e-04,+.1145e-03,+.1192e-03,+.4934e-04,+.2699e-04,&
     &                      +.3684e-04,-.2395e-05,+.2045e-04,-.8684e-04,&
     &                                 +.5301e-04,-.4176e-05,+.4103e-04,&
     &                                            +.2783e-04,+.1754e-04,&
     &                                                       +.1116e-04/
DATA zhc1/                                                              &
     &-.3450e+02,+.2148e+03,+.3376e+02,+.6535e+02,-.1564e+02,-.4273e+02,&
     &           +.9553e+01,-.4647e+01,-.6129e+01,-.6727e+01,-.6761e+01,&
     &                      -.2467e+01,-.2181e+01,-.5361e+01,-.2395e+01,&
     &                                 +.5952e+00,+.2106e+01,-.1367e+01,&
     &                                            -.2349e+01,+.3532e+00,&
     &                                                       -.3169e+01/
DATA zhs1/                                                              &
     &           +.3977e+01,+.5032e+01,+.6226e+01,-.3625e+00,-.1373e+01,&
     &                      +.4600e+01,+.4312e+01,+.2882e+01,-.6351e+01,&
     &                                 +.5731e+01,-.2574e+01,+.3235e+00,&
     &                                            +.2806e+01,+.8133e+00,&
     &                                                       +.2032e+01/
DATA zqc2/                                                              &
     &+.8571e-03,+.3086e-02,+.9287e-03,+.2787e-03,+.1826e-03,-.1006e-03,&
     &           +.1092e-03,-.1266e-03,+.5372e-04,-.1188e-03,-.3285e-04,&
     &                      -.1783e-04,-.3018e-05,-.8709e-04,-.8707e-04,&
     &                                 +.8633e-04,+.3530e-04,+.4863e-04,&
     &                                            +.3917e-05,-.3252e-04,&
     &                                                       -.1936e-06/
DATA zqs2/                                                              &
     &           -.8822e-04,+.1341e-03,+.3095e-04,+.8230e-04,+.2735e-04,&
     &                      +.1714e-04,-.9406e-04,+.1912e-04,-.5402e-04,&
     &                                 +.3571e-04,+.3897e-04,+.4487e-04,&
     &                                            +.3079e-04,+.3196e-04,&
     &                                                       -.2391e-05/
DATA zhc2/                                                              &
     &+.5216e+02,+.1613e+03,+.3284e+02,-.7670e+02,-.9548e+01,+.1608e+02,&
     &           +.1023e+02,-.1090e+02,+.2748e+01,-.3846e+01,-.4135e+01,&
     &                      +.1255e+01,-.3301e-01,-.5273e+01,-.7247e+01,&
     &                                 +.1387e+02,+.4184e+01,+.6495e+01,&
     &                                            +.2944e+01,-.1947e+01,&
     &                                                       +.1132e+01/
DATA zhs2/                                                              &
     &           -.1968e+02,+.1192e+02,-.1194e+01,+.1084e+01,+.2946e+01,&
     &                      +.2630e+01,-.1256e+02,+.1395e+01,-.2222e+01,&
     &                                 +.4864e+01,+.6450e+01,+.5568e+01,&
     &                                            +.5292e+01,+.4876e+01,&
     &                                                       -.7579e+00/
DATA zqc3/                                                              &
     &-.2759e-03,-.2781e-03,-.1087e-03,-.1633e-03,-.3627e-04,-.4242e-04,&
     &           +.6045e-05,-.1703e-04,+.4562e-04,-.1009e-04,+.2663e-04,&
     &                      -.1786e-04,+.1550e-04,-.9135e-06,+.2372e-04,&
     &                                 +.1100e-05,+.2299e-04,+.4659e-05,&
     &                                            +.2423e-05,+.7321e-05,&
     &                                                       +.8852e-05/
DATA zqs3/                                                              &
     &           -.3678e-04,-.2219e-04,-.3911e-04,-.4398e-04,-.1142e-04,&
     &                      -.9121e-05,-.2011e-04,+.4711e-06,-.3775e-05,&
     &                                 +.3866e-05,+.2400e-04,+.2043e-04,&
     &                                            -.1824e-05,-.5550e-05,&
     &                                                       +.2506e-05/
DATA zhc3/                                                              &
     &-.1534e+03,-.2095e+02,-.1006e+03,-.7385e+01,+.5203e+01,+.9434e+00,&
     &           -.3814e+00,-.3175e+01,+.3366e+01,+.3378e+00,+.2740e+00,&
     &                      -.2669e+01,+.8452e+00,+.3498e+00,+.2192e+01,&
     &                                 -.4024e+00,+.1544e+01,-.4588e+00,&
     &                                            +.6998e+00,+.6263e+00,&
     &                                                       +.1228e+01/
DATA zhs3/                                                              &
     &           -.3588e+01,+.2076e+00,-.2088e+01,-.4159e+01,+.2244e+00,&
     &                      -.7751e+00,-.2749e+01,+.7234e+00,+.4390e+00,&
     &                                 -.1646e+00,+.1700e+01,+.1046e+01,&
     &                                            -.7856e+00,-.1644e+01,&
     &                                                       +.2648e+00/
DATA zqc4/                                                              &
     &-.1460e-03,+.3422e-03,-.3529e-04,+.1791e-03,-.1917e-03,-.2558e-04,&
     &           +.6547e-04,+.6401e-04,+.4823e-04,+.7084e-05,+.2895e-04,&
     &                      -.1561e-04,+.8179e-06,+.1028e-04,-.7667e-05,&
     &                                 -.4347e-05,+.7293e-05,-.5735e-05,&
     &                                            +.7838e-05,-.2933e-05,&
     &                                                       +.3686e-05/
DATA zqs4/                                                              &
     &           -.4560e-05,-.5292e-04,-.1252e-04,+.1850e-04,-.2273e-04,&
     &                      +.6552e-05,+.1422e-04,-.6545e-05,+.7998e-06,&
     &                                 +.2845e-04,+.2497e-04,+.2844e-04,&
     &                                            +.3855e-06,-.1487e-04,&
     &                                                       +.1954e-05/
DATA zhc4/                                                              &
     &+.9260e+01,-.9055e+01,+.5460e+01,-.7603e+01,-.3329e+02,-.1048e+02,&
     &           +.9328e+01,+.4597e+01,+.3827e+01,-.3201e+01,+.1708e+01,&
     &                      -.1548e+01,-.5323e+00,+.3039e+01,+.5740e+00,&
     &                                 +.1353e+00,-.2354e+01,+.2818e+00,&
     &                                            +.1113e+01,-.1891e+01,&
     &                                                       -.3074e+00/
DATA zhs4/                                                              &
     &           -.2446e+01,+.4199e+01,-.2571e+01,+.8194e+01,+.4206e+00,&
     &                      +.3856e+01,+.1159e+01,+.2547e+01,-.1314e+01,&
     &                                 +.2331e+01,+.1144e+01,-.4408e+00,&
     &                                            -.6797e+00,-.2598e+01,&
     &                                                       +.8953e+00/
!
!     ------------------------------------------------------------------
!
!*         1.     PRELIMINARY SETTING.
!                 ----------- --------
zytime=pytime
!
!     ------------------------------------------------------------------
!
!*         2.     COMPUTATIONS.
!                 -------------
!
zc1yt=COS(zytime)
zs1yt=SIN(zytime)
zc2yt=zc1yt**2-zs1yt**2
zs2yt=2.*zs1yt*zc1yt

DO jmn=1,21
  pqc(jmn)=zqc0(jmn)+2.*(zqc1(jmn)*zc1yt+zqc2(jmn)*zs1yt &
                        +zqc3(jmn)*zc2yt+zqc4(jmn)*zs2yt)
  phc(jmn)=zhc0(jmn)+2.*(zhc1(jmn)*zc1yt+zhc2(jmn)*zs1yt &
                        +zhc3(jmn)*zc2yt+zhc4(jmn)*zs2yt)
END DO

DO jmn=1,15
  pqs(jmn)=zqs0(jmn)+2.*(zqs1(jmn)*zc1yt+zqs2(jmn)*zs1yt &
                        +zqs3(jmn)*zc2yt+zqs4(jmn)*zs2yt)
  phs(jmn)=zhs0(jmn)+2.*(zhs1(jmn)*zc1yt+zhs2(jmn)*zs1yt &
                        +zhs3(jmn)*zc2yt+zhs4(jmn)*zs2yt)
END DO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE ozone

!==============================================================================
!+ Very Old ECMWF routine.
!------------------------------------------------------------------------------

SUBROUTINE legtri (psin, kcp, palp)

!------------------------------------------------------------------------------
!
! Description:
!   *LEGTRI* - *LEGENDRE FUNKTION FUER EINE DREIECK-ABSCHNEIDUNG
!   THIS ROUTINE COMPUTES THE VALUES *PALP* FOR THE ARGUMENT
!   *PSIN* OF THE NORMALISED *LEGENDRE ASSOCIATED FUNCTIONS IN THE
!   ORDER ((JN1=JM1,KCP),JM1=1,KCP) FOR JN=JN1-1 AND JM=JM1-1 .
!
!
! Method:
!   SIMPLE RECURRENCE FORMULA.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
REAL    (KIND=ireals),    INTENT(IN)    ::   &
  psin      ! sine of latitude

INTEGER (KIND=iintegers), INTENT(IN)    ::   &
  kcp       ! 1 + limit wave number

! Array arguments with intent(out):
REAL    (KIND=ireals),    INTENT(OUT)   ::   &
  palp(66)  ! array of results

!-----------------------------------------------------------------------
!
! Local scalars:
REAL (KIND=ireals)                      ::                       &
      zsin, zcos, zf1m, zm, z2m, zre1, ze1, zf2m, zn, zn2, ze2

INTEGER (KIND=iintegers)                ::                       &
      icp, ic, jj, jm, jm1, jm2, jn

!- End of header
!------------------------------------------------------------------------------

!
!     ------------------------------------------------------------------
!
!*         1.     PRELIMINARY SETTING.
!                 ----------- --------
!
!
ZSIN=PSIN
ICP=KCP
!
!     ------------------------------------------------------------------
!
!*         2.     COMPUTATIONS.
!                 -------------
!
!
      IC=ICP-1
      ZCOS=SQRT(1.-ZSIN**2)
      JJ=2
      PALP(1)=1.
      ZF1M=SQRT(3.)
      PALP(2)=ZF1M*ZSIN
      DO 203 JM1=1,ICP
      JM=JM1-1
      ZM=JM
      Z2M=ZM+ZM
      ZRE1=SQRT(Z2M+3.)
      ZE1=1./ZRE1
      IF(JM.EQ.0) GO TO 201
      ZF2M=ZF1M*ZCOS/SQRT(Z2M)
      ZF1M=ZF2M*ZRE1
      JJ=JJ+1
      PALP(JJ)=ZF2M
      IF(JM.EQ.IC) GO TO 203
      JJ=JJ+1
      PALP(JJ)=ZF1M*ZSIN
      IF(JM1.EQ.IC) GO TO 203
  201 CONTINUE
      JM2=JM+2
      DO 202 JN=JM2,IC
      ZN=JN
      ZN2=ZN**2
      ZE2=SQRT((4.*ZN2-1.)/(ZN2-ZM**2))
      JJ=JJ+1
      PALP(JJ)=ZE2*(ZSIN*PALP(JJ-1)-ZE1*PALP(JJ-2))
      ZE1=1./ZE2
  202 CONTINUE
  203 CONTINUE
!
!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE legtri

!==============================================================================
!+ Computation of some ground fields.
!------------------------------------------------------------------------------

SUBROUTINE ground_fields (kevert)

!------------------------------------------------------------------------------
!
! Description:
!   Compute some ground fields:
!         ground temperatures                (t_s, t_snow, t_m, t_cl)
!         spec./rel. humidity at the surface (qv_s, rh_s)
!         water content of soil layers       (w_g1, w_g2, w_g3, w_cl)
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(IN)   ::  &
  kevert               ! number of levels of target model

!------------------------------------------------------------------------------
!
! Local parameters:
REAL (KIND=ireals), PARAMETER ::                                       &
  zfacwb = 0.2,  & ! factor for splitting the 10cm ground water layer
  zd1lmw = 0.1,  & ! over 3 prognostic layers.
  zd2lmw = 0.9     ! thicknesses for 2 prog. layers ( nl_soil_lm = 2 )

! Local scalars:
INTEGER (KIND=iintegers) ::   &
  i, mt, j, l, k, s_t,        &
  mzwg1_loc_lm, mzwg2_loc_lm, & ! Locations of output variables
  mzwg3_loc_lm, mzwcl_loc_lm    ! in output variable table


REAL    (KIND=ireals)    ::   &
  zdpgl, zdplm, zdzgl, zdzlm, zfac1, zfac2, zcf_snow

! Local arrays:
REAL    (KIND=ireals)    ::   &
  zt_g_lm(ie2lm,je2lm),   & ! for computation of qv_s we need t_g_lm
  zrh_slm(ie2lm,je2lm),   & !
  zdfis  (ie2lm,je2lm)

LOGICAL                  ::   &
  lzmask (ie2lm,je2lm)      ! to call the routine tgcom

!- End of header
!------------------------------------------------------------------------------

IF (l_smi) THEN
  DO i = 1, nvar_lm
    SELECT CASE (var_lm(i)%name)
      CASE ('W_G1      '); mzwg1_loc_lm = i
      CASE ('W_G2      '); mzwg2_loc_lm = i
      CASE ('W_G3      '); mzwg3_loc_lm = i
      CASE ('W_CL      '); mzwcl_loc_lm = i
    END SELECT
  ENDDO
ENDIF

!------------------------------------------------------------------------------
! Section 1: Temperatures
!------------------------------------------------------------------------------

  ! Compute the temperature at the ground surface (t_s_lm) and
  !         the temperature at the snow surface (t_snow_lm)
  ! NOTE: for the multi-layer soil model, also t_s_lm is used here, but all
  !       input variables (t_s_gl, etc.) are derived from t_so_gme or t_so_in

  zdfis(:,:) = fis_gl(:,:) - fis_lm(:,:)

  IF (lmulti_layer_lm) THEN
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) < 0.5) THEN
          ! set freshsnw_lm to 0.0 for seapoints
          freshsnw_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  DO j = 1, je2lm
    DO i = 1, ie2lm
      IF (fr_land_lm(i,j) < 0.5) THEN
        t_s_lm(i,j)    = t_s_gl(i,j)
        t_snow_lm(i,j) = t_s_gl(i,j)

        ! set w_snow_lm to 0.0 (explicitly)
        w_snow_lm  (i,j) = 0.0_ireals
      ELSE
        ! Compute the Prandtl layer thickness (first layer above ground)
        ! for the input grid (coarse, input) and LM/HM (fine, output)
        IF (llm2lm) THEN
          ! For both grids (coarse and fine) the height between the surface 
          ! layer and the lowest HALF level (Prandtl layer thickness) is used
          ! to compute the factor for the temperature gradient.
          zdzgl = hhl_gl(i,j,ke_in) - hhl_gl(i,j,ke_in+1)
          zdzlm = hhl_lm(i,j,kelm)  - hhl_lm(i,j,kelm+1)
          t_s_gl(i,j)  = t_lm (i,j,kevert) - dtkes_gl(i,j) * zdzlm/zdzgl
        ELSE
          ! the pressure is used to compute that factor
          zdpgl = ps_gl(i,j) - (ak_in(ke_in) + bk_in(ke_in) * ps_gl(i,j))
          zdplm = ps_lm(i,j) - (ak_in(ke_in) + bk_in(ke_in) * ps_lm(i,j))
          t_s_gl(i,j)    = t_lm (i,j,kevert) - dtkes_gl(i,j) * zdplm/zdpgl
! iso code: perform height correction for water isotopes in soil. Hui Tang 2013-11-20
          risosoil_lm(i,j,1)  = riso_lm (i,j,kevert,1) - drisoke_gl(i,j,1) * zdplm/zdpgl  
          risosoil_lm(i,j,2)  = riso_lm (i,j,kevert,2) - drisoke_gl(i,j,2) * zdplm/zdpgl  
! end iso code
        ENDIF
        t_s_lm(i,j)    = t_s_gl(i,j)
        t_snow_lm(i,j) = t_s_lm(i,j)    - dtssnow_gl(i,j)
      ENDIF
    ENDDO
  ENDDO

  ! Correct temperature and height of the snow
  ! ------------------------------------------
  WHERE (w_snow_lm(:,:) == 0.0_ireals)
    t_snow_lm(:,:) = t_s_lm(:,:)
  END WHERE

  WHERE ( (w_snow_lm(:,:) /= 0.0_ireals)                               .AND. &
        ( (t_snow_lm(:,:) > T0+2.0_ireals) .OR. (zdfis(:,:) > 4000.0)) .AND. &
          ! Added dependency on the latitude
            (latlm_m(:,:) > -60._ireals) )
    w_snow_lm(:,:) = 0.0
    t_snow_lm(:,:) = t_s_lm(:,:)
  END WHERE

  WHERE ( (w_snow_lm(:,:) /= 0.0_ireals) .AND.                      &
        ((t_snow_lm(:,:) <= T0+2.0_ireals) .AND. (zdfis(:,:) <= 4000.0)) )
    t_snow_lm(:,:) = MIN(t_snow_lm(:,:), T0)
       t_s_lm(:,:) = MIN(   t_s_lm(:,:), T0 - 0.8_ireals)
  END WHERE

  ! Temperature of the ground
  ! -------------------------

  IF ( (.NOT. lmulti_layer_in) .AND.  (itype_t_cl == 0) ) THEN
    IF (linitial) THEN
      IF (lt_cl_corr) THEN
        ! A height correction of the climatological temperature is performed
        ! A fraction of the difference between the GME-height and the LM-height
        ! is added to the interpolated value. The value 0.007 is a heuristic
        ! value.
        zfac1 = 0.007_ireals / G
        zfac2 = 0.007_ireals
        WHERE (fr_land_lm(:,:) >= 0.5_ireals)
          t_cl_lm(:,:) = t_cl_lm(:,:) + fis_gl(:,:)*zfac1 - hsurf_lm(:,:)*zfac2
        END WHERE
      ELSE
        ! This is the old type of height correction dating back to the old
        ! ECMWF spectral global model (GM).
        WHERE (fr_land_lm(:,:) >= 0.5_ireals)
          t_cl_lm(:,:) = 0.5_ireals*(t_s_lm(:,:) + dtms_gl(:,:) + t_cl_lm(:,:))
        END WHERE
      ENDIF
      WHERE (fr_land_lm(:,:) < 0.5_ireals)
        t_cl_lm(:,:) = 0.0_ireals
      END WHERE
    ENDIF

    WHERE (fr_land_lm(:,:) >= 0.5_ireals)
       t_m_lm(:,:) =      t_s_lm(:,:) + dtms_gl(:,:)
    END WHERE
    WHERE (fr_land_lm(:,:) < 0.5_ireals)
       t_m_lm(:,:) = 0.0_ireals
    END WHERE
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Specific and relative humidity at the ground surface
!------------------------------------------------------------------------------

  zcf_snow    = 0.0150_ireals
  lzmask(:,:) = (fr_land_lm(:,:) >= 0.5_ireals)

  CALL tgcom ( zt_g_lm (:,:), t_snow_lm(:,:), t_s_lm (:,:), w_snow_lm(:,:),  &
               lzmask(:,:) , ie2lm, je2lm, zcf_snow, 1, ie2lm, 1, je2lm)

  DO j = 1, je2lm
    DO i = 1, ie2lm
      IF (fr_land_lm(i,j) < 0.5) THEN
        qv_s_lm(i,j) =  Rdv     * psat_w(zt_g_lm(i,j), b1, b2_w, b3, b4_w)/  &
          (ps_lm(i,j) - O_m_rdv * psat_w(zt_g_lm(i,j), b1, b2_w, b3, b4_w))
        zrh_slm(i,j)  = 1.00
      ELSE
        qv_s_lm(i,j) = rh_s_gl(i,j)  *                                       &
                        Rdv     * psat_w(zt_g_lm(i,j), b1, b2_w, b3, b4_w)/  &
          (ps_lm(i,j) - O_m_rdv * psat_w(zt_g_lm(i,j), b1, b2_w, b3, b4_w))
        zrh_slm(i,j)  = rh_s_gl(i,j)
      ENDIF

      ! when interpolating initial data to very high resolution, some problems
      ! were reported with negative qv_s values. Therefore we do a clipping here
      IF (qv_s_lm(i,j) < 0.0_ireals) THEN
        qv_s_lm(i,j) = 0.0_ireals
      ENDIF
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 3: Multiply soil moisture contents by pore volume
!------------------------------------------------------------------------------

  IF (.NOT. lmulti_layer_in) THEN
    DO j = 1, je2lm
      DO i = 1, ie2lm
        s_t = NINT(soiltyp_lm(i,j))
        IF (s_t >= 3 .AND. s_t <= 8) THEN
          IF (l_smi) THEN
            ! Scale the soil moisture with the soil moisture index SMI:
            ! smi = (sm-pwp)/(fc-pwp)
            ! Caution: sm is in volumetric units, soil depth in cm!
            w_g1_lm(i,j) = w_g1_lm(i,j) * (fcb(s_t)-pwpb(s_t)) + pwpb(s_t) * &
                (var_lm(mzwg1_loc_lm)%levbot-var_lm(mzwg1_loc_lm)%levtop)*0.01
            w_g2_lm(i,j) = w_g2_lm(i,j) * (fcb(s_t)-pwpb(s_t)) + pwpb(s_t) * &
                (var_lm(mzwg2_loc_lm)%levbot-var_lm(mzwg2_loc_lm)%levtop)*0.01
            IF (nl_soil_lm == 3) THEN
              w_g3_lm(i,j) = w_g3_lm(i,j)*(fcb(s_t)-pwpb(s_t)) + pwpb(s_t) * &
                (var_lm(mzwg3_loc_lm)%levbot-var_lm(mzwg3_loc_lm)%levtop)*0.01
            ENDIF
            IF (linitial) THEN
              w_cl_lm(i,j) = w_cl_lm(i,j)*(fcb(s_t)-pwpb(s_t)) + pwpb(s_t) * &
                (var_lm(mzwcl_loc_lm)%levbot-var_lm(mzwcl_loc_lm)%levtop)*0.01
            ENDIF
          ELSE
            w_g1_lm(i,j) = w_g1_lm(i,j) * porb(s_t)
            w_g2_lm(i,j) = w_g2_lm(i,j) * porb(s_t)
            IF (nl_soil_lm == 3) THEN
              w_g3_lm(i,j) = w_g3_lm(i,j)*porb(s_t)
            ENDIF
            IF (linitial) THEN
              w_cl_lm(i,j) = w_cl_lm(i,j)*porb(s_t)
            ENDIF
          ENDIF ! l_smi
        ELSE
          w_g1_lm(i,j) = 0.0_ireals
          w_g2_lm(i,j) = 0.0_ireals
          IF (nl_soil_lm == 3) THEN
            w_g3_lm(i,j) = 0.0_ireals
          ENDIF
          IF (linitial) THEN
            w_cl_lm(i,j) = 0.0_ireals
          ENDIF
        ENDIF ! s_t
      ENDDO
    ENDDO
  ENDIF ! .NOT. lmulti_layer_in

!------------------------------------------------------------------------------
! Section 4: New soil model
!------------------------------------------------------------------------------

  ! after introduction of new SR init_multi_layer_gme_ml, the following has only
  ! to be done, if (ke_soil_lm == ke_soil_in).
  ! But this only holds for lgme2lm; other interpolations do need that and
  ! for lgme2lm it does not harm

  IF (lmulti_layer_in .AND. (.NOT. lcomp_bound) .AND. &
     ((ke_soil_lm == ke_soil_in) .OR. lec2lm)) THEN
    ! Computation of Temperatures
    DO j = 1, je2lm
      DO i = 1, ie2lm
        ! t_so_lm(:,:,0) has to be computed for all grid points
        IF (t_s_lm(i,j) /= undef) THEN
          t_so_lm(i,j,0) = t_s_lm(i,j)
        ELSE
          t_so_lm(i,j,0) = undef
        ENDIF
      ENDDO
    ENDDO
    DO k = 1, ke_soil_lm+1
      DO j = 1, je2lm
        DO i = 1, ie2lm
          IF (fr_land_lm(i,j) >= 0.5_ireals) THEN
            IF (t_so_lm(i,j,k-1) /= undef) THEN
              t_so_lm(i,j,k) = t_so_lm(i,j,k-1) + dt_so_gl(i,j,k)
            ELSE
              t_so_lm(i,j,k) = undef
            ENDIF
          ELSE
            t_so_lm(i,j,k) = 0.0_ireals
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! Computation of Soil Moistures
    IF (.NOT. lcm2lm) THEN
      ! in case of lcm2lm, all computations are done in init_multi_layer_cm
      IF (l_smi) THEN
        ! Scale the soil moisture with the soil moisture index SMI:
        ! smi = (sm-pwp)/(fc-pwp)
        ! Caution: sm is in volumetric units, soil depth in cm!
        DO k = 1, ke_soil_lm+1
          DO j = 1, je2lm
            DO i = 1, ie2lm
              s_t = NINT(soiltyp_lm(i,j))
              IF (s_t >= 3 .AND. s_t <= 8) THEN
                IF (w_so_lm(i,j,k) /= undef) THEN
                  w_so_lm(i,j,k) = w_so_lm(i,j,k)*(fcb(s_t)-pwpb(s_t))+pwpb(s_t)
                  IF (k == 1) THEN
                    w_so_lm(i,j,k) = w_so_lm(i,j,k)*                    0.01
                  ELSE
                    w_so_lm(i,j,k) = w_so_lm(i,j,k)*(3**(k-1)-3**(k-2))*0.01
                  ENDIF
                ELSE
                  w_so_lm(i,j,k) = 0_ireals
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO k = 1, ke_soil_lm+1
          DO j = 1, je2lm
            DO i = 1, ie2lm
              s_t = NINT(soiltyp_lm(i,j))
              IF (s_t >= 3 .AND. s_t <= 8) THEN
                IF (w_so_lm(i,j,k) /= undef) THEN
                  w_so_lm(i,j,k) = w_so_lm(i,j,k) * porb(s_t)
                ELSE
                  w_so_lm(i,j,k) = undef
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF ! l_smi
    ENDIF ! lcm2lm
  ENDIF ! lmulti_layer_in

!------------------------------------------------------------------------------
! Section 5: Profiles
!------------------------------------------------------------------------------

! ! Keep the values at some gridpoints, if wanted
! IF (lprgp) THEN
!   DO l = 1, ngp
!     tsglpr (l) = t_s_gl    (igp(l), jgp(l))
!     tsnowpr(l) = t_snow_lm (igp(l), jgp(l))
!     tspr   (l) = t_s_lm    (igp(l), jgp(l))
!     tmpr   (l) = t_m_lm    (igp(l), jgp(l))
!     IF (linitial)        tclpr (l) = t_cl_lm (igp(l), jgp(l))
!     wsnowpr(l) = w_snow_lm (igp(l), jgp(l))
!     wg1pr  (l) = w_g1_lm   (igp(l), jgp(l))
!     wg2pr  (l) = w_g2_lm   (igp(l), jgp(l))
!     IF (nl_soil_lm == 3) wg3pr (l) = w_g3_lm (igp(l), jgp(l))
!     IF (linitial)        wclpr (l) = w_cl_lm (igp(l), jgp(l))
!     vio3pr (l) = vio3_lm   (igp(l), jgp(l))
!     hmo3pr (l) = hmo3_lm   (igp(l), jgp(l))
!     plcpr  (l) = plcov_lm  (igp(l), jgp(l))
!     laigp  (l) = laindex_lm(igp(l), jgp(l))
!     rootpr (l) = root_lm   (igp(l), jgp(l))
!     rh_sglpr(l)= rh_s_gl   (igp(l), jgp(l))
!     rh_spr (l) = zrh_slm   (igp(l), jgp(l))
!     frlapr (l) = fr_land_lm(igp(l), jgp(l))
!   ENDDO
! ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE ground_fields

!==============================================================================
!+ Compute plant characteristics for a given day in the year
!------------------------------------------------------------------------------

SUBROUTINE plant_characteristics (nactday, ydate, izdebug)

!------------------------------------------------------------------------------
!
! Description:
!
!------------------------------------------------------------------------------

! Parameterlist
INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  nactday,      & ! actual day in the year
  izdebug         ! for debug print out

CHARACTER (LEN=14),       INTENT(IN)  ::  &
  ydate           ! actual date in the form yyyymmddhhmmss

!------------------------------------------------------------------------------

! Local variables
INTEGER (KIND=iintegers)   ::  &
  nzday, i, j, izerror,        &
  imo1, imo2                        ! for interpolating monthly values

REAL (KIND=ireals)         ::  &
  zhred, zbvp, zdvp,           &    ! factors for vegetation period
  zwei1, zwei2                      ! weights for interpolating monthly values

REAL (KIND=ireals)         ::  &
  zvegfac(ie2lm,je2lm)              ! vegetation factor depending on grid
                                    ! point and time of year
CHARACTER (LEN=80)         ::  &
  yzerrmsg

!------------------------------------------------------------------------------

  izerror  = 0
  yzerrmsg = '    '

  !----------------------------------------------------------------------------
  ! Section 1a: Compute vegetation factor for every grid point
  !----------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, '   Compute vegetation factor for every grid point'
  ENDIF

  ! The actual state of the vegetation for a given day of the year
  ! is taken from data of the vegetation period and from the time of rest.
  DO j = 1, je2lm
    DO i = 1, ie2lm
      ! zhred:  height reduction factor
      ! zbvp:   begin of the vegetation period (in day of the year) as a
      !         function of the latitude
      ! zdvp:   length of the vegetation period (in days) as a
      !         function of the latitude
      ! zbvp & zdvp: from climatic charts of monthly mean 2m temperature
      !         > 5-10 C
      zhred = EXP(-5.0E-9_ireals * fis_lm(i,j)**2)
      zbvp  = MAX(1.0_ireals, 3.0_ireals*(ABS(latlm_m(i,j)) - 20.0_ireals))
      zdvp  = MIN(365.0_ireals,                                           &
               345.0_ireals - 4.5_ireals*(ABS(latlm_m(i,j)) - 20.0_ireals))

      IF ( latlm_m(i,j) < 0.0_ireals ) THEN
        nzday = MOD(nactday + 180_iintegers, 365_iintegers)
      ELSE
        nzday = nactday
      ENDIF

      IF ( zdvp >= 345.0_ireals ) THEN
        zvegfac(i,j) = zhred
      ELSE IF ( nzday < NINT(zbvp) ) THEN
        zvegfac(i,j) = 0.0_ireals
      ELSE IF ( nzday > NINT(zbvp+zdvp) ) THEN
        zvegfac(i,j) = 0.0_ireals
      ELSE
        zvegfac(i,j) = MAX(0.0_ireals, MIN ( 1.0_ireals, 1.12_ireals*      &
                      SIN( Pi*MAX(0.0_ireals, (nzday-zbvp))/zdvp ) ) )*zhred
      ENDIF
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 1b: Compute factors for interpolating monthly values
  !             and actual ndvi ratio, if wanted
  !----------------------------------------------------------------------------

  IF ( (itype_ndvi == 1) .OR. (itype_ndvi == 2) .OR. (itype_aerosol == 2) ) THEN

    IF (izdebug > 5) THEN
      PRINT *, '   Compute factors for monthly interpolation'
    ENDIF
    CALL month2hour (ydate, imo1, imo2, zwei1, zwei2, izdebug)

    IF (itype_ndvi == 1) THEN
      ndviratio_lm(:,:) = ndvi_mrat_lm(:,:,imo1) +                    &
                         (ndvi_mrat_lm(:,:,imo2) - ndvi_mrat_lm(:,:,imo1)) * zwei2
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Root depth
!------------------------------------------------------------------------------

  ! Choose between different treatments for rootdp

  ! rootdp_lm is always computed from the external parameter values rootdp_mx
  ! Another option tpye_rootdp=3 for the calculation of the rootdepth has been
  ! introduced by B. Rockel for the CLM

  IF     (itype_rootdp == 0) THEN
    ! Put an annual cylce to the input values
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals .AND. NINT(soiltyp_lm(i,j)) > 2) THEN
          rootdp_lm(i,j) = MIN (rootdp_mx(i,j),                             &
                              0.12_ireals + zvegfac(i,j)**2 * 0.58_ireals)
        ELSE
          rootdp_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ELSEIF (itype_rootdp == 1) THEN
    ! Just take the input values with a min of 0.12
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals .AND. NINT(soiltyp_lm(i,j)) > 2) THEN
          rootdp_lm(i,j) = MAX (rootdp_mx(i,j), 0.12_ireals)
        ELSE
          rootdp_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ELSEIF (itype_rootdp == 2) THEN
    ! Adapt the annual cylce to ECOCLIMAP niveau (was introduced by CLM)
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals .AND. NINT(soiltyp_lm(i,j)) > 2) THEN
          rootdp_lm(i,j) = MIN (rootdp_mx(i,j),                             &
                              0.12_ireals + zvegfac(i,j)**2 * 0.58_ireals)
          ! Scaling of rootdp_lm
          rootdp_lm(i,j) = MIN (rootdp_lm(i,j),1.0_ireals)
          rootdp_lm(i,j) = -0.23_ireals + 2.05_ireals* rootdp_lm(i,j)
          rootdp_lm(i,j) = MAX (rootdp_lm(i,j), 0.12_ireals)
        ELSE
          rootdp_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ELSEIF (itype_rootdp == 3) THEN
    ! Put an annual cylce to the input value without maximum cut off
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals .AND. NINT(soiltyp_lm(i,j)) > 2) THEN
          rootdp_lm(i,j) = MAX (zvegfac(i,j)**2 * rootdp_mx(i,j),0.12_ireals)
        ELSE
          rootdp_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ELSEIF (itype_rootdp == 4) THEN
    ! Just take the input values
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals .AND. NINT(soiltyp_lm(i,j)) > 2) THEN
          rootdp_lm(i,j) = rootdp_mx(i,j)
        ELSE
          rootdp_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Plant cover
!------------------------------------------------------------------------------

  IF     (itype_ndvi == 0) THEN
    ! compute the actual plant cover from plcov_mx_lm and plcov_mn_lm
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals .AND. NINT(soiltyp_lm(i,j)) > 2) THEN
          plcov_lm(i,j) = plcov_mn_lm(i,j) +                                 &
                          zvegfac(i,j) * (plcov_mx_lm(i,j) - plcov_mn_lm(i,j))
        ELSE
          plcov_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ELSEIF (itype_ndvi == 1) THEN
    ! compute the actual plant cover from plcov_mx_lm and the actual
    ! NDVI ratio (ndviratio)
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals .AND. NINT(soiltyp_lm(i,j)) > 2) THEN
          plcov_lm(i,j) = plcov_mx_lm(i,j) * ndviratio_lm (i,j)
        ELSE
          plcov_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ELSEIF (itype_ndvi == 2) THEN
    ! _dl 25.02.09  use prescribed yearly cycle from climatology data
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals) THEN
          plcov_lm(i,j) = plcov12(i,j,imo1) +                    &
                         (plcov12(i,j,imo2) - plcov12(i,j,imo1)) * zwei2
        ELSE
          plcov_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Leaf area index
!------------------------------------------------------------------------------

  IF     (itype_ndvi == 0) THEN
  ! Compute the actual leaf area index from lai_mx_lm and lai_mn_lm
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals .AND. NINT(soiltyp_lm(i,j)) > 2) THEN
          lai_lm(i,j) = lai_mn_lm(i,j) +                             &
                          zvegfac(i,j) * (lai_mx_lm(i,j) - lai_mn_lm(i,j))
        ELSE
          lai_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ELSEIF (itype_ndvi == 1) THEN
    ! compute the actual leaf area index from lai_mx_lm and the actual
    ! NDVI ratio (ndviratio)
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals .AND. NINT(soiltyp_lm(i,j)) > 2) THEN
          lai_lm(i,j) = lai_mx_lm(i,j) * ndviratio_lm(i,j)
        ELSE
          lai_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ELSEIF (itype_ndvi == 2) THEN
    ! _dl 25.02.09  use prescribed yearly cycle from climatology data
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals) THEN
          lai_lm(i,j) = lai12(i,j,imo1) +                    &
                       (lai12(i,j,imo2) - lai12(i,j,imo1)) * zwei2
        ELSE
          lai_lm(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ENDIF

!_dl 15.04.09  roughness length is variable for itype_ndvi == 2
!------------------------------------------------------------------------------
! Section 5: roughness length
!------------------------------------------------------------------------------

  IF     (itype_ndvi == 2) THEN
    DO j = 1, je2lm
      DO i = 1, ie2lm
        z0_lm(i,j) = z012(i,j,imo1) +                    &
                    (z012(i,j,imo2) - z012(i,j,imo1)) * zwei2
      ENDDO
    ENDDO
  ENDIF

!---------------------------------------------------------------------
! End of Subroutine plant_characteristics
!---------------------------------------------------------------------

END SUBROUTINE plant_characteristics

!==============================================================================
!+ Initial values of soil water and temperature on multi layer grid
!------------------------------------------------------------------------------

SUBROUTINE init_multi_layer

!------------------------------------------------------------------------------
!
! Description:
!   The module procedure init_multi_layer uses the initial values of the 
!   2/3-level version for temperature and ground water to provide initial
!   values of soil temperature and ground water for the multi layer soil model.
!   The same vertical resolution has to be taken for water and temperature 
!   calculation. ke_soil_lm is the number of active soil layers for temperature 
!   and water
!
! Method:
! a) soil water:
!    initial water distribution by interpolation of the 2/3-layer water content
!    values using mass conservation. Each layer of the 2/3-layer version must 
!    contain at least one level of the multi-layer version.
!
!    It is allowed that the lowest multi layer level is below the lowest 
!    2/3-layer version level. In this case this lowest layer uses the 2/3-layer
!    water climate value.
!
! b) soil temperature:
!    multi level initial temperature distribution is interpolated from GME 
!    temperatures using enthalpie conservation. Multi layer levels below EFR 
!    climate level are allowed. Temperatures in these levels are calculated by
!    linear extrapolation of the GME soil temperatures
!
! Code provided by: DWD, Reinhold Schrodin
!  phone:  +49  69  8062 2709
!  fax:    +49  69  8062 3721
!  email:  Reinhold.Schrodin@dwd.de
!
!------------------------------------------------------------------------------

! Local parameters:
! ----------------
  REAL    (KIND=ireals   ), PARAMETER ::  &
    zepsi  = 1.0E-6_ireals      ! security constant

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i, j   ,        & ! loop indices
    kso    ,        & ! soil loop index in z-direction
    k1     ,        & ! loop index in z-direction for soil water interpolation
    k2     ,        & ! loop index in z-direction for soil water interpolation
    kontr  ,        & ! for soil temperature profile extrapolation control
    k2_save,        & ! saves k2 for initial state of soil water interpolation
    jb     ,        & ! loop index for soil-type
    mstyp             ! soil type index

  REAL    (KIND=ireals   ) ::  &
    zwmmin, zwgmin, zwgmax, zwmmax, zx, zomb, zzdlam, zwmean, zrhoc, z4wdpv, &
    zalfbs, zbetms, zdmddb, zwqm, zeo, zeu, ztzhlto, ztzhltu, ztzhlt,        &
    zzalam(10), zdzwg(4), zwg(4), zwb(4)

! Scalars took over by LM (wsth the values as in LM)
  REAL    (KIND=ireals   ) ::  &
    ctau1  = 1.0000  , & !  first adjustment time period in EFR-method
    ctau2  = 5.0000  , & !  second adjustment time period in the EFR-method
    cdzw12 = 0.1000  , & !  thickness of upper soil water layer in two-layer model
    cdzw22 = 0.9000  , & !  thickness of lower soil water layer in two-layer model
    cdzw13 = 0.0200  , & !  thickness of upper soil water layer in three-layer model
    cdzw23 = 0.0800  , & !  thickness of middle soil water layer in three-layer model
    cdzw33 = 0.9000  , & !  thickness of lower soil water layer in three-layer model
    rho_w  = 1000.0  , & !  density of liquid water
    chc_w  = 4180.0      !  heat capacity of water

  REAL    (KIND=ireals   ) ::  &
    ! parameters for the determination of the soil heat conductivity (W/(K*m)), 
    ! pore volume, field capacity, etc.
    cala0(10) = (/2.26, 2.41, 0.30, 0.28, 0.25, 0.21, 0.18, 0.06, 1.0, 2.26/), &
    cala1(10) = (/2.26, 2.41, 2.40, 2.40, 1.58, 1.55, 1.50, 0.50, 1.0, 2.26/), &
    cporv(10) = (/1.E-10,1.E-10, 0.364, 0.445, 0.455, 0.475, 0.507, 0.863, 1.E-10, 1.E-10 /), &
    cfcap(10) = (/1.E-10,1.E-10, 0.196, 0.260, 0.340, 0.370, 0.463, 0.763, 1.E-10, 1.E-10 /), &
    cpwp (10) = (/0.0, 0.0, 0.042, 0.100, 0.110, 0.185, 0.257, 0.265, 0.0,  0.0/),            &
    crhoc(10) = (/1.92E6,2.10E6,1.28E6,1.35E6,1.42E6,1.50E6,1.63E6, 0.58E6, 4.18E6, 1.92E6/), &
    cadp (10) = (/0.0, 0.0, 0.012, 0.030, 0.035, 0.060, 0.065, 0.098, 0.0, 0.0/)            , &
    cdz1 (10)

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=ireals   ) ::  &
    zdzg   (ie2lm,je2lm)          , & !
    zdzm   (ie2lm,je2lm)          , & !
    zfcap  (ie2lm,je2lm)          , & !
    zpwp   (ie2lm,je2lm)          , & !
    zadp   (ie2lm,je2lm)          , & ! air dryness point
    zporv  (ie2lm,je2lm)          , & !
    zrocg  (ie2lm,je2lm)          , & !
    zrocm  (ie2lm,je2lm)          , & !
    zmls   (ke_soil_lm+1)         , & ! depth of soil layer full levels
    zhls   (ke_soil_lm+1)         , & ! depth of soil layer half levels
    xdzhs  (ke_soil_lm+1)         , & !
    xdzms  (ke_soil_lm+1)         , & !
    zwg_fr (ie2lm,je2lm,ke_soil_lm+1), & ! fractional water content of layers in ground
    zroc   (ie2lm,je2lm,ke_soil_lm+1)    !

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine init_multi_layer
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

! layer thickness for hydrological part
  IF(nl_soil_lm == 2) THEN    ! (2-layer model)
    zdzwg(1)  = cdzw12
    zdzwg(2)  = cdzw22
    zwg(1) = zdzwg(1)          ! depth of first half level of 2-level version
    zwg(2) = zwg(1) + zdzwg(2) ! depth of second  "   "    "      "
  ELSE                      ! (3-layer model)
    zdzwg(1)  = cdzw13
    zdzwg(2)  = cdzw23
    zdzwg(3)  = cdzw33
    zwg(1) = zdzwg(1)           ! depth of 1st half level of 3-level version
    zwg(2) = zwg(1) + zdzwg(2)  ! depth of 2d   "     "      "         "
    zwg(3) = zwg(2) + zdzwg(3)  ! depth of 3d   "     "      "         "
  ENDIF
  zdzwg(nl_soil_lm+1)= zdzwg(nl_soil_lm)
  zwg(nl_soil_lm+1) = zwg(nl_soil_lm) + zdzwg(nl_soil_lm+1) 
                                ! depth of 3d or 4th half level

  zomb    = 2.*pi/(ctau1*86400.)
  zx      = SQRT(ctau1/ctau2)
  zalfbs  = 1. + zx + zx*zx
  zbetms  = zx*SQRT(1. + zx*zx)*EXP( -zx/(1. + zx))
  zdmddb  = zalfbs/zbetms - 1.
  DO jb = 1, 10
    zzdlam   = cala1(jb) - cala0(jb)
    zwmean   = 0.5*(cfcap(jb) + cpwp(jb))
    z4wdpv   = 4.*zwmean/cporv(jb)
    zrhoc    = crhoc (jb) + rho_w*zwmean*chc_w
    zzalam(jb)   = cala0(jb) + zzdlam*(0.25 + 0.30*zzdlam   &
                 /(1.0 + 0.75*zzdlam))*MIN(z4wdpv,1.0_ireals+(z4wdpv-1.0_ireals)  &
                 *(1.0 + 0.35*zzdlam)/(1.0 + 1.95*zzdlam) )
    cdz1(jb) = SQRT (2.*zzalam(jb)/(zrhoc*zomb))/(1.+zx)
  ENDDO

  DO j = 1, je2lm
    DO i = 1, ie2lm
      IF (fr_land_lm(i,j) >= 0.5_ireals) THEN   ! for land-points only
        mstyp = NINT(soiltyp_lm(i,j))              ! soil type
        zporv(i,j)  = cporv(mstyp)              ! pore volume
        zpwp (i,j)  = cpwp (mstyp)              ! plant wilting point
        zadp (i,j)  = cadp (mstyp)              ! air dryness point
        zfcap(i,j)  = cfcap(mstyp)              ! field capacity
        zrocg(i,j)  = crhoc(mstyp)              ! heat capacity
        zwqm         = 0.5*(zfcap(i,j) + zpwp(i,j))
        zrocm(i,j)   = zrocg(i,j) + rho_w*zwqm*chc_w 
        zrocg(i,j)   = zrocm(i,j) 
        zdzg (i,j)  = cdz1 (mstyp)      ! top layer thickness
        zdzm (i,j)  = zdzg (i,j)*zdmddb ! second layer thickness
      ENDIF
    ENDDO
  ENDDO

! Definition of grids for temperature and water content

!  zhls(1)  = 2.*czmls_lm(1)    !depth of first half level
  zhls(1)  = czhls_lm(1)       !depth of first half level
  xdzhs(1) = zhls(1)           !layer thickness betw. half levels of 1st layer
  zmls(1)  = czmls_lm(1)       !depth of 1st main level
  xdzms(1) = czmls_lm(1)       !layer thickness betw. full levels of 1st layer

  IF (my_cart_id == 0) THEN
    WRITE (*,'(A)') ' kso     zhls      zmls     xdzhs     xdzms'
    kso=1
    WRITE (*,'(I3,4F10.3)') kso, zhls(kso), zmls(kso), xdzhs(kso), xdzms(kso)
  ENDIF
  DO kso = 2,ke_soil_lm+1
!    zhls(kso)  = zhls(kso-1) + 2.*(czmls_lm(kso) -zhls(kso-1))  !_br czmls -> czmls_lm
    zhls(kso)  = czhls_lm(kso)  !_br czmls -> czmls_lm
    xdzhs(kso) = zhls(kso) - zhls(kso-1)   ! layer thickness between half levels
    zmls(kso)  = czmls_lm(kso)             ! depth of main levels
    xdzms(kso) = zmls(kso) - zmls(kso-1)   ! layer thickness betw. main levels
    IF (my_cart_id == 0) THEN
      WRITE (*,'(I3,4F10.3)') kso, zhls(kso), zmls(kso), xdzhs(kso), xdzms(kso)
    ENDIF
  ENDDO

! soil temperature
  DO j = 1, je2lm
    DO i = 1, ie2lm
!US   IF (fr_land_lm(i,j) >= 0.5_ireals) THEN   ! for land-points only
      ! this is needed for all grid points
        t_so_lm(i,j,0) = t_s_lm(i,j)
!US   ENDIF
    END DO
  END DO
  DO j = 1, je2lm
    DO i = 1, ie2lm
      IF (fr_land_lm(i,j) >= 0.5_ireals) THEN   ! for land-points only
        kontr = 0
        DO kso = 1, ke_soil_lm+1
          IF(czmls_lm(kso) <= zdzg(i,j)) THEN
            t_so_lm(i,j,kso) = t_s_lm(i,j)+(t_m_lm(i,j)-t_s_lm(i,j))* &
                               zmls(kso)/zdzg(i,j)
            zroc(i,j,kso) = zrocg(i,j)
          ELSE IF (czmls_lm(kso).GT.zdzg(i,j).AND.czmls_lm(kso-1).LE.zdzg(i,j)) THEN
            zroc(i,j,kso) = (zrocg(i,j)*(zdzg(i,j) - czmls_lm(kso-1)) + &
                           zrocm(i,j)*(czmls_lm(kso)- zdzg(i,j)))/  &
                                                (czmls_lm(kso) - czmls_lm(kso-1))
            ztzhlto = t_m_lm(i,j)-(zdzg(i,j)-czmls_lm(kso-1))/zdzg(i,j)*  &
                                   (t_m_lm(i,j) -  t_s_lm(i,j))
            zeo = zrocg(i,j)*(zdzg(i,j)-czmls_lm(kso-1))*0.5*(t_m_lm(i,j) + ztzhlto)
            ztzhltu = t_m_lm(i,j) + (czmls_lm(kso) - zdzg(i,j))/zdzm(i,j)* &
                                             (t_cl_lm(i,j) - t_m_lm(i,j))
            zeu = zrocm(i,j)*(czmls_lm(kso) - zdzg(i,j))* &
                            0.5*(t_m_lm(i,j) + ztzhltu)
            t_so_lm(i,j,kso) = (zeo + zeu)/zroc(i,j,kso)/(czmls_lm(kso) - czmls_lm(kso-1))
          ELSE IF (czmls_lm(kso-1).GT.zdzg(i,j)) THEN
            zroc(i,j,kso) = zrocm(i,j)
            t_so_lm(i,j,kso) = t_m_lm(i,j) + (t_cl_lm(i,j) - t_m_lm(i,j))/ &
                                           zdzm(i,j)*(zmls(kso) - zdzg(i,j))
            IF (czmls_lm(kso).GT.(zdzg(i,j) + zdzm(i,j))) THEN
              IF (kontr.EQ.0) THEN
                ztzhlt = t_cl_lm(i,j) + (czmls_lm(kso) - (zdzg(i,j) + zdzm(i,j)))* &
                                       (t_cl_lm(i,j) - t_m_lm(i,j))/zdzm(i,j)
                kontr = 1
              ELSE
                t_so_lm(i,j,kso) = ztzhlt
              END IF

              ! below temperature climate level of 2/3-layer version: 
              ! t_so = t_cl
              t_so_lm(i,j,kso) = t_cl_lm(i,j)
            END IF
          END IF
        END DO

        zwgmin        = zadp (i,j)*zdzwg(1)     ! lower-limit for 1. layer
        zwgmax        = zporv(i,j)*zdzwg(1)     ! upper-limit for 1. layer
        zwmmin        = zadp (i,j)*zdzwg(2)     ! lower-limit for 2. layer
        zwmmax        = zporv(i,j)*zdzwg(2)     ! upper-limit for 2. layer
        w_g1_lm(i,j) = MAX( zwgmin, MIN(zwgmax,w_g1_lm(i,j)) )
        zwb(1) = w_g1_lm(i,j)
        w_g2_lm(i,j) = MAX( zwmmin, MIN(zwmmax,w_g2_lm(i,j)) )
        zwb(2) = w_g2_lm(i,j)
        IF (nl_soil_lm == 2) THEN
          w_cl_lm(i,j)     = MAX(zwmmin,MIN(zwmmax,w_cl_lm(i,j)    ))
          zwb(3) = w_cl_lm(i,j)
        ELSE
          zwmmin        = zadp (i,j)*zdzwg(3)    ! lower-limit for 3. layer
          zwmmax        = zporv(i,j)*zdzwg(3)    ! upper-limit for 3. layer
          w_g3_lm(i,j) = MAX( zwmmin, MIN(zwmmax,w_g3_lm(i,j)) )
          zwb(3) = w_g3_lm(i,j)
          w_cl_lm(i,j)       = MAX( zwmmin, MIN(zwmmax,w_cl_lm(i,j)      ) )
          zwb(4) = w_cl_lm(i,j)
        ENDIF

!       multi level water distribution derived from 2/3-level version
        k2_save = 1
        DO k1 =1,nl_soil_lm+1 ! 2/3 level index
          DO k2 = k2_save,ke_soil_lm+1 ! more layer index
            ! zwg_fr-determination if k2-layer completely in k1-layer included
            zwg_fr(i,j,k2) = zwb(k1)/zdzwg(k1)
            w_so_lm(i,j,k2)     = zwg_fr(i,j,k2)*xdzhs(k2)

            IF (czmls_lm(k2) > zwg(nl_soil_lm+1)) THEN !w_so_lm = climate value
              zwg_fr(i,j,k2) = zwb(nl_soil_lm+1)/zdzwg(nl_soil_lm+1)
              w_so_lm(i,j,k2) = zwg_fr(i,j,k2)*xdzhs(k2)
              k2_save=k2+1
              CYCLE
            ENDIF   !  (czmls_lm(k2) > zwg(nl_soil_lm+1))

            ! zwg_fr-determination if k2-layer is separated by a k1-level
            IF (czmls_lm(k2) > zwg(k1)) THEN
              zwg_fr(i,j,k2) = (zwb(k1)/zdzwg(k1)*(zwg(k1)-czmls_lm(k2-1)) +    &
                                zwb(k1+1)/zdzwg(k1+1)*(czmls_lm(k2)-zwg(k1)))/&
                               (czmls_lm(k2) - czmls_lm(k2-1))
              w_so_lm(i,j,k2) = zwg_fr(i,j,k2)*xdzhs(k2)
              k2_save = k2 + 1
              EXIT
            END IF  ! (czmls_lm(k2) > zwg(k1)....
          END DO
        END DO
      END IF  ! (landmask..
    END DO
  END DO

  ! Initialize the fresh-snow indicator with 1 
  ! (means all the snow is fresh)
  WHERE (fr_land_lm(:,:) >= 0.5_ireals)
    freshsnw_lm(:,:) = 1.0_ireals
  END WHERE

!---------------------------------------------------------------------
! End of Subroutine init_multi_layer
!---------------------------------------------------------------------

END SUBROUTINE init_multi_layer

!==============================================================================
!+ Initial values of soil water and temperature on multi layer grid for climate models
!------------------------------------------------------------------------------

SUBROUTINE init_multi_layer_cm

!------------------------------------------------------------------------------
!
! Description:
!   The module procedure init_multi_layer uses the initial values of the 
!   2/3-level version for temperature and ground water to provide initial
!   values of soil temperature and ground water for the multi layer soil model.
!   The same vertical resolution has to be taken for water and temperature 
!   calculation. ke_soil_in is the number of active soil layers for temperature 
!   and water
!
! Method:
! a) soil water:
!
! b) soil temperature:
!
! Code provided by: GKSS, Burkhardt Rockel
!
!------------------------------------------------------------------------------

! Local parameters:
! ----------------
  REAL    (KIND=ireals   ), PARAMETER ::  &
    zepsi  = 1.0E-6_ireals      ! security constant

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i, j,           & ! loop indices
    kso     ,        & ! soil loop index in z-direction
    kso_in  ,        & ! soil loop index in z-direction
    kso_in2 ,        & ! soil loop index in z-direction
    k1      ,        & ! loop index in z-direction for soil water interpolation
    k2      ,        & ! loop index in z-direction for soil water interpolation
    k2_save ,        & ! saves k2 for initial state of soil water interpolation
    jb      ,        & ! loop index for soil-type
    mstyp              ! soil type index

  REAL    (KIND=ireals   ) ::  &
    zwmmin, zwgmin, zwgmax, zwmmax, zx, zomb, zzdlam, zwmean, zrhoc, z4wdpv, &
    zalfbs, zbetms, zdmddb, zwqm,        &
    zzalam(10), zdzwg(4), zwg(4), zwb(4), &
    zvw_tmp, zweight

! Scalars took over by LM (wsth the values as in LM)
  REAL    (KIND=ireals   ) ::  &
    rho_w  = 1000.0  , & !  density of liquid water
    chc_w  = 4180.0      !  heat capacity of water

  REAL    (KIND=ireals   ) ::  &
    ! parameters for the determination of the soil heat conductivity (W/(K*m)), 
    ! pore volume, field capacity, etc.
    cala0(10) = (/2.26, 2.41, 0.30, 0.28, 0.25, 0.21, 0.18, 0.06, 1.0, 2.26/), &
    cala1(10) = (/2.26, 2.41, 2.40, 2.40, 1.58, 1.55, 1.50, 0.50, 1.0, 2.26/), &
    cporv(10) = (/1.E-10,1.E-10, 0.364, 0.445, 0.455, 0.475, 0.507, 0.863, 1.E-10, 1.E-10 /), &
    cfcap(10) = (/1.E-10,1.E-10, 0.196, 0.260, 0.340, 0.370, 0.463, 0.763, 1.E-10, 1.E-10 /), &
    cpwp (10) = (/0.0, 0.0, 0.042, 0.100, 0.110, 0.185, 0.257, 0.265, 0.0,  0.0/),            &
    crhoc(10) = (/1.92E6,2.10E6,1.28E6,1.35E6,1.42E6,1.50E6,1.63E6, 0.58E6, 4.18E6, 1.92E6/), &
    cadp (10) = (/0.0, 0.0, 0.012, 0.030, 0.035, 0.060, 0.065, 0.098, 0.0, 0.0/)            , &
    cdz1 (10)

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=ireals   ) ::  &
    zw_so_rel    (ie2lm, je2lm, ke_soil_lm+1) , & !
    zw_so_rel_in (ie_in, je_in, ke_soil_lm+1) , & !
    zadp   (ie2lm,je2lm)                 , & ! air dryness point
    zporv  (ie2lm,je2lm)                 , & !
    zhls   (ke_soil_lm+1)                , & ! depth of soil layer half levels
    zmls   (ke_soil_lm+1)                , & ! depth of soil layer full levels
    xdzhs  (ke_soil_lm+1)                    ! layer thickness between half levels

REAL    (KIND=ireals   ) ::  &
    zfcap   (ie2lm,je2lm)                    ! field capacity

! Local variables:
  INTEGER  (KIND=iintegers)  ::  &
    izerror                 ! status and error status variable

  LOGICAL                    ::  &
    lzmono, lzposdef

  CHARACTER (LEN=  1)        ::  &
    yzitype     ! interpolation type

  CHARACTER (LEN=200)        ::  &
    yzerrmsg    ! error message for error handling

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine init_multi_layer_cm
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  zw_so_rel = 0._ireals
  lzmono    = .FALSE.
  lzposdef  = .FALSE.
  yzitype   = 'M'

! Definition of grids for temperature and water content

  zmls(1)  = czmls_lm(1)                   !depth of 1st main level
  xdzhs(1) = czhls_lm(1)       !layer thickness betw. half levels of 1st layer
!  xdzhs(1) = 2._ireals * czmls_lm(1)       !layer thickness betw. half levels of 1st layer
  zhls(1)  = xdzhs(1)

! SP, 201405
!  DO kso = 2,ke_soil_lm+1
  DO kso = 1,ke_soil_lm+1
    zhls(kso)  = czhls_lm(kso)             ! depth of half levels
!    zhls(kso)  = zhls(kso-1) + 2.*(czmls_lm(kso) -zhls(kso-1))  ! depth of half levels
    zmls(kso)  = czmls_lm(kso)             ! depth of main levels
    xdzhs(kso) = zhls(kso) - zhls(kso-1)   ! layer thickness between half levels
  ENDDO

! SP, 201405: re-written this part of the code, as official version gives
! strange results (unknown reason...)
  IF (itype_w_so_rel /= 0) THEN
    zw_so_rel_in(:,:,:) = 0._ireals
    DO j = 1, je_in
      DO i = 1, ie_in
        IF (lolp_in(i,j)) THEN
          DO kso = 1, ke_soil_lm+1
            IF (czhls_lm(kso-1) > czhls_in(ke_soil_in+1)) THEN
              ! layer below lowest input level
              zw_so_rel_in(i,j,kso) = w_so_rel_in(i,j,ke_soil_in+1) * xdzhs(kso)
            ELSE
              DO kso_in = 1, ke_soil_in+1
                IF (czhls_lm(kso-1) >= czhls_in(kso_in-1) .AND. czhls_lm(kso-1) < czhls_in(kso_in)) THEN
                  IF (czhls_lm(kso) <= czhls_in(kso_in)) THEN
                    ! layer totally contained in input layer
                    zw_so_rel_in (i,j,kso) = w_so_rel_in(i,j,kso_in) * xdzhs(kso)
                    EXIT
                  ELSE
                    ! layer spread over several input layers
                    zvw_tmp = w_so_rel_in(i,j,kso_in) * (czhls_in(kso_in) - czhls_lm(kso-1))
                    DO kso_in2 = kso_in+1, ke_soil_in+1
                      IF (czhls_lm(kso) > czhls_in(kso_in2)) THEN
                        zvw_tmp = zvw_tmp + w_so_rel_in(i,j,kso_in2) * (czhls_in(kso_in2) - czhls_in(kso_in2-1))
                      ELSE
                        zvw_tmp = zvw_tmp + w_so_rel_in(i,j,kso_in2) * (czhls_lm(kso) - czhls_in(kso_in2-1))
                        EXIT
                      ENDIF
                    ENDDO
                    IF (czhls_lm(kso) > czhls_in(ke_soil_in+1))                                                      &
                      zvw_tmp = zvw_tmp + w_so_rel_in(i,j,ke_soil_in+1) * (czhls_lm(kso) - czhls_in(ke_soil_in+1))
                    zw_so_rel_in (i,j,kso) = zvw_tmp
                    EXIT
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO

!  IF (itype_w_so_rel /= 0) THEN
!    zw_so_rel_in(:,:,:) = 0._ireals
!    DO j = 1, je_in
!      DO i = 1, ie_in
!        IF (lolp_in(i,j)) THEN
!         multiplication (integration) of relative input soil moisture with layer thickness
!         i.e. transition to absolute soil moisture, scaled on input pore volume or field capacity
!          DO kso_in = 1, ke_soil_in+1
!            w_so_rel_in(i,j,kso_in)=w_so_rel_in(i,j,kso_in) * (czhls_in(kso_in) - czhls_in(kso_in-1))
!          ENDDO
!          DO kso = 1, ke_soil_lm
!          DO kso = 1, ke_soil_lm+1
!            IF (czhls_lm(kso) > czhls_in(ke_soil_in+1)) THEN
! UB
!              IF(itype_w_so_rel > 2) THEN
!                zw_so_rel_in(i,j,kso) = w_so_rel_in(i,j,ke_soil_in+1) &
!                                      * (czhls_lm(kso)-czhls_lm(kso-1))/(czhls_in(ke_soil_in+1) - czhls_in(ke_soil_in))
!              ELSE
! UB
!                zw_so_rel_in(i,j,kso) = w_so_rel_in(i,j,ke_soil_in+1)
! UB
!              ENDIF
! UB
!            ELSE
!              zvw_tmp = 0._ireals
!              zweight = 0._ireals
!              DO kso_in = 1, ke_soil_in+1
!                IF (czhls_lm(kso-1) >= czhls_in(kso_in-1) .AND. czhls_lm(kso-1) < czhls_in(kso_in)) THEN
!                  IF (czhls_lm(kso) >= czhls_in(kso_in-1) .AND. czhls_lm(kso) < czhls_in(kso_in)) THEN
! UB
!                    IF(itype_w_so_rel > 2) THEN
!                      zw_so_rel_in (i,j,kso) = w_so_rel_in(i,j,kso_in) &
!                                             * (czhls_lm(kso)-czhls_lm(kso-1))/(czhls_in(kso_in) - czhls_in(kso_in-1))
!                    ELSE
! UB
!                      zw_so_rel_in (i,j,kso) = w_so_rel_in(i,j,kso_in)
! UB
!                    ENDIF
! UB
!                  ELSE
!                    zweight = zweight + czhls_in(kso_in) -czhls_lm(kso-1)
!                    zvw_tmp = zvw_tmp + w_so_rel_in(i,j,kso_in) * (czhls_in(kso_in) - czhls_lm(kso-1))
!                    DO kso_in2 = kso_in+1, ke_soil_in+1
!                      IF (czhls_lm(kso) > czhls_in(kso_in2)) THEN
!                       zweight = zweight + czhls_in(kso_in2) - czhls_in(kso_in2-1)
!                        zvw_tmp = zvw_tmp + w_so_rel_in(i,j,kso_in2)*(czhls_in(kso_in2) - czhls_in(kso_in2-1))
!                      ELSE
!                        zweight = zweight + czhls_lm(kso) - czhls_in(kso_in2-1)
!                        zvw_tmp = zvw_tmp + w_so_rel_in(i,j,kso_in2) * (czhls_lm(kso) - czhls_in(kso_in2-1))
!                        EXIT
!                      ENDIF
!                    ENDDO
!                    zw_so_rel_in(i,j,kso) = zvw_tmp / zweight
!                    EXIT
!                  ENDIF
!                ENDIF
!              ENDDO
!            ENDIF
!          ENDDO
!        ENDIF
!      ENDDO
!    ENDDO

    ! this interpolation was also done for itype_w_so_rel=0 (artifical profile)
    ! but cannot work without bigger adaptations for GME2LM
    ! But I think it is not necessary: when using an artificial profile
    ! directly set zw_so_rel!

    DO kso = 1,ke_soil_lm+1
      CALL interp_l(zw_so_rel_in(:,:,kso), ie_in, je_in, i_index(:,:,1),    &
       j_index(:,:,1), lzmono, lzposdef, yzitype, lbd_frame_cur,            &
       lolp_in, lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),     &
       zw_so_rel(:,:,kso), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,    &
       latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,             &
       grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                      &
       yzerrmsg, izerror)
    ENDDO

! SP, 201405: skip this for this code version (important for official version!)
!      DO j = 1, je2lm
!      DO i = 1, ie2lm
!        DO kso = 1,ke_soil_lm+1
!          zw_so_rel(i,j,kso) = zw_so_rel(i,j,kso) * (czhls_lm(kso) - czhls_lm(kso-1))
!        ENDDO
!      ENDDO
!    ENDDO

  ELSE

    ! itype_w_so_rel = 0: use an artificial profile
    DO j = 1, je2lm
      DO i = 1, ie2lm
        zw_so_rel(i,j,:) = cw_so_rel_lm(:)
        DO kso = 1,ke_soil_lm+1
          zw_so_rel(i,j,kso) = zw_so_rel(i,j,kso) * (czhls_lm(kso) - czhls_lm(kso-1))
        ENDDO
      ENDDO
    ENDDO

  ENDIF

  ! here was the interpolation before

  DO j = 1, je2lm
    DO i = 1, ie2lm
      IF (lolp_lm(i,j)) THEN   ! for land-points only
        mstyp = NINT(soiltyp_lm(i,j))              ! soil type
        zporv(i,j)  = cporv(mstyp)              ! pore volume
        zadp (i,j)  = cadp (mstyp)              ! air dryness point
        zfcap(i,j)  = cfcap(mstyp)              ! field capacity
      ENDIF
    ENDDO
  ENDDO

! soil temperature and soil water
  DO j = 1, je2lm
    DO i = 1, ie2lm
      IF (lolp_lm(i,j)) THEN   ! for land-points only
!        DO kso = 1, ke_soil_lm
        DO kso = 1, ke_soil_lm+1    !_br 10.12.10
            t_so_lm(i,j,kso) = t_cl_lm(i,j) + (t_s_lm(i,j) - t_cl_lm(i,j)) * EXP(-zmls(kso)/2.8)
            zwgmin        = zadp (i,j)*xdzhs(kso)
! UB
!           IF (itype_w_so_rel < 2) THEN
!BR            IF     (itype_w_so_rel <= 1 .OR. itype_w_so_rel == 3 ) THEN
            IF     (itype_w_so_rel <= 1 .OR. itype_w_so_rel == 3 ) THEN
! UB
                zwgmax        = zporv(i,j)*xdzhs(kso)
                w_so_lm(i,j,kso) = MAX( zwgmin, MIN(zwgmax, zw_so_rel(i,j,kso) * zporv(i,j)) )
! UB
!           ELSE
            ELSEIF (itype_w_so_rel == 2 .OR. itype_w_so_rel == 4 ) THEN
! UB
                zwgmax        = zfcap(i,j)*xdzhs(kso)
                w_so_lm(i,j,kso) = MAX( zwgmin, MIN(zwgmax, zw_so_rel(i,j,kso) * zfcap(i,j)) )
            ENDIF
        END DO
        t_so_lm(i,j,ke_soil_lm+1) = t_cl_lm(i,j)
!        w_cl_lm(i,j) = w_so_lm(i,j,ke_soil_lm+1)   !_br 10.12.10
      END IF  ! (landmask)
    END DO
  END DO

  ! this is needed for GRIB format
  t_so_lm(:,:,0) = t_s_lm(:,:)

! SP, 201405
!  DO kso = 1, ke_soil_lm
  DO kso = 1, ke_soil_lm+1
    WHERE (soiltyp_lm < 1.5_ireals) w_so_lm(:,:,kso) = 0._ireals
  ENDDO

  ! US: this has also to be included for the COSMO-Model multilayer soil model
  !     also rho_snow is needed, but this can be initialized in the COSMO-Model
  !     by setting lana_rho_snow=.FALSE. in GRIBIN

  ! Initialize the fresh-snow indicator with 1
  ! (means all the snow is fresh)
  WHERE (fr_land_lm(:,:) >= 0.5_ireals)
    freshsnw_lm(:,:) = 1.0_ireals
  END WHERE

!------------------------------------------------------------------------------
! End of Subroutine init_multi_layer_cm
!------------------------------------------------------------------------------

END SUBROUTINE init_multi_layer_cm

!==============================================================================
!==============================================================================
!+ Initial values of soil water and temperature on multi layer grid for GME-ML
!------------------------------------------------------------------------------

SUBROUTINE init_multi_layer_gme_ml

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
! a) soil water:
!
! b) soil temperature:
!
! Code adapted from init_multi_layer and init_multi_layer_cm by:
!   DWD, Susanne Brienen
!
!------------------------------------------------------------------------------

! Local parameters:
! ----------------

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i, j,           & ! loop indices
    kso     ,        & ! soil loop index in z-direction
    kso_in  ,        & ! soil loop index in z-direction
    mstyp              ! soil type index

  REAL    (KIND=ireals   ) ::  &
    zwgmin, zwgmax, zeo, zeu, ztzhlto, ztzhltu

! Scalars took over by LM (wsth the values as in LM)
  REAL    (KIND=ireals   ) ::  &
    rho_w  = 1000.0  , & !  density of liquid water
    chc_w  = 4180.0      !  heat capacity of water

  REAL    (KIND=ireals   ) ::  &
    ! parameters for the determination of the soil heat conductivity (W/(K*m)),
    ! pore volume, field capacity, etc.
    cala0(10) = (/2.26, 2.41, 0.30, 0.28, 0.25, 0.21, 0.18, 0.06, 1.0, 2.26/), &
    cala1(10) = (/2.26, 2.41, 2.40, 2.40, 1.58, 1.55, 1.50, 0.50, 1.0, 2.26/), &
    cporv(10) = (/1.E-10,1.E-10, 0.364, 0.445, 0.455, 0.475, 0.507, 0.863, 1.E-10, 1.E-10 /), &
    cfcap(10) = (/1.E-10,1.E-10, 0.196, 0.260, 0.340, 0.370, 0.463, 0.763, 1.E-10, 1.E-10 /), &
    cpwp (10) = (/0.0, 0.0, 0.042, 0.100, 0.110, 0.185, 0.257, 0.265, 0.0,  0.0/),            &
    crhoc(10) = (/1.92E6,2.10E6,1.28E6,1.35E6,1.42E6,1.50E6,1.63E6, 0.58E6, 4.18E6, 1.92E6/), &
    cadp (10) = (/0.0, 0.0, 0.012, 0.030, 0.035, 0.060, 0.065, 0.098, 0.0, 0.0/)            , &
    cdz1 (10)

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=ireals   ) ::  &
    zdzg   (ie2lm,je2lm)          , & !
    zdzm   (ie2lm,je2lm)          , & !
    zfcap  (ie2lm,je2lm)          , & ! field capacity
    zpwp   (ie2lm,je2lm)          , & !
    zadp   (ie2lm,je2lm)          , & ! air dryness point
    zporv  (ie2lm,je2lm)          , & !
    zrocg  (ie2lm,je2lm)          , & ! heat capacity of bare soil   
    zmls   (ke_soil_lm+1)         , & ! depth of soil layer full levels
    zhls   (ke_soil_lm+1)         , & ! depth of soil layer half levels
    xdzhs  (ke_soil_lm+1)         , & ! layer thickness between half levels
    xdzms  (ke_soil_lm+1)         , & ! layer thickness between full levels
    zwg_fr (ie2lm,je2lm,ke_soil_lm+1),    & ! fractional water content of layers in ground
    zwg_fr_in(ie2lm,je2lm,ke_soil_in+1),  & ! fractional water content of input layers
    w_so_gl_int(ie2lm,je2lm,ke_soil_in+1),& ! intermediate soil moisture field
    t_so_int(ie2lm,je2lm,0:ke_soil_in+1), & ! intermediate soil temperature field
    zroc   (ie2lm,je2lm,ke_soil_lm+1), & ! heat capacity of LM layers
    zroc_in(ie2lm,je2lm,ke_soil_in+1)    ! heat capacity of input layers

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine init_multi_layer_gme_ml
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  DO j = 1, je2lm
    DO i = 1, ie2lm
      IF (fr_land_lm(i,j) >= 0.5_ireals) THEN   ! for land-points only
        mstyp = NINT(soiltyp_lm(i,j))              ! soil type
        zporv(i,j)  = cporv(mstyp)              ! pore volume
        zpwp (i,j)  = cpwp (mstyp)              ! plant wilting point
        zadp (i,j)  = cadp (mstyp)              ! air dryness point
        zfcap(i,j)  = cfcap(mstyp)              ! field capacity
        zrocg(i,j)  = crhoc(mstyp)              ! heat capacity
      ENDIF
    ENDDO
  ENDDO

! Definition of grids for temperature and water content

!  zhls(1)  = 2.*czmls_lm(1)    !depth of first half level
  zhls(1)  = czhls_lm(1)       !depth of first half level
  xdzhs(1) = zhls(1)           !layer thickness betw. half levels of 1st layer
  zmls(1)  = czmls_lm(1)       !depth of 1st main level
  xdzms(1) = czmls_lm(1)       !layer thickness betw. full levels of 1st layer

  IF (my_cart_id == 0) THEN
    WRITE (*,'(A)') ' kso     zhls      zmls     xdzhs     xdzms'
    kso=1
    WRITE (*,'(I3,4F10.3)') kso, zhls(kso), zmls(kso), xdzhs(kso), xdzms(kso)
  ENDIF
  DO kso = 2,ke_soil_lm+1
!    zhls(kso)  = zhls(kso-1) + 2.*(czmls_lm(kso) -zhls(kso-1))  !_br czmls -> czmls_lm
    zhls(kso)  = czhls_lm(kso)  !_br czmls -> czmls_lm
    xdzhs(kso) = zhls(kso) - zhls(kso-1)   ! layer thickness between half levels
    zmls(kso)  = czmls_lm(kso)             ! depth of main levels
    xdzms(kso) = zmls(kso) - zmls(kso-1)   ! layer thickness betw. main levels
    IF (my_cart_id == 0) THEN
      WRITE (*,'(I3,4F10.3)') kso, zhls(kso), zmls(kso), xdzhs(kso), xdzms(kso)
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
! Section 2: Interpolation of soil temperature and soil moisture
!------------------------------------------------------------------------------

! first, convert interpolated temperature differences into absolute values on input layers
    DO j = 1, je2lm
      DO i = 1, ie2lm
        ! t_so_lm(:,:,0) has to be computed for all grid points
        IF (t_s_lm(i,j) /= undef) THEN
          t_so_int(i,j,0) = t_s_lm(i,j)
        ELSE
          t_so_int(i,j,0) = undef
        ENDIF
      ENDDO
    ENDDO
    DO kso_in = 1, ke_soil_in+1
      DO j = 1, je2lm
        DO i = 1, ie2lm
          IF (fr_land_lm(i,j) >= 0.5_ireals) THEN
            IF (t_so_int(i,j,kso_in-1) /= undef) THEN
              t_so_int(i,j,kso_in) = t_so_int(i,j,kso_in-1) + dt_so_gl(i,j,kso_in)
            ELSE
              t_so_int(i,j,kso_in) = undef
            ENDIF
          ELSE
            t_so_int(i,j,kso_in) = 0.0_ireals
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  DO j = 1, je2lm
    DO i = 1, ie2lm
      t_so_lm(i,j,0) = t_so_int(i,j,0)
    END DO
  END DO
  DO j = 1, je2lm
    DO i = 1, ie2lm
      IF (fr_land_lm(i,j) >= 0.5_ireals) THEN   ! for land-points only

        ! -- soil moisture

        ! uppermost layer
! TEST: falls mit SWI:
!        w_so_lm(i,j,1) = w_so_gl(i,j,1)/czmls_in(1)*czmls_lm(1)
        ! back transformation from soil wetness index
!        w_so_lm(i,j,1) = (w_so_lm(i,j,1)+zpwp(i,j)*plcov_lm(i,j)/(zfcap(i,j)-zpwp(i,j)))&
!                          /(plcov_lm(i,j)/(zfcap(i,j)-zpwp(i,j))+(1-plcov_lm(i,j))/zfcap(i,j))
        !w_so_lm(i,j,1) = w_so_lm(i,j,1)*(zfcap(i,j)-zpwp(i,j)) + zpwp(i,j)
! sonst:
        zwgmin = zadp(i,j)*czmls_in(1)
        zwgmax = zporv(i,j)*czmls_in(1)

        w_so_gl_int(i,j,1) = MAX( zwgmin, MIN(zwgmax, w_so_gl(i,j,1)*zporv(i,j)))
        w_so_lm(i,j,1) = w_so_gl_int(i,j,1)/czhls_in(1)*czhls_lm(1)/czmls_in(1)*czmls_lm(1)

        ! other layers
        DO kso = 2, ke_soil_lm + 1

          ! new soil levels below all old levels
          IF (czmls_lm(kso) > czmls_in(ke_soil_in+1)) THEN
              w_so_lm(i,j,kso) = w_cl_lm(i,j)
! oder so?
!                zwg_fr(i,j,kso) = w_so_gl(i,j,ke_soil_in+1)/(czhls_in(kso_in+1) - czhls_in(kso_in))
!                w_so_lm(i,j,kso) = zwg_fr(i,j,kso)*xdzhs(kso)
          ELSE

            DO kso_in = 2, ke_soil_in+1

! TEST: falls mit SWI:
! in w_so_gl steht zunaechst der soil wetness index. dieser wird einfach linear
! interpoliert und dann wieder zuruecktransformiert
!              IF (czmls_lm(kso) >= czmls_in(kso_in-1) .AND. czmls_lm(kso) < czmls_in(kso_in)) THEN
!                w_so_lm(i,j,kso) = 0.5*(w_so_gl(i,j,kso_in-1)/(czhls_in(kso_in) - czhls_in(kso_in-1)) &
!                                   + w_so_gl(i,j,kso_in)/(czhls_in(kso_in+1) - czhls_in(kso_in))) &
!                                      * xdzhs(kso)
!              ENDIF


! ohne SWI:
              zwgmin = zadp(i,j)*(czhls_in(kso_in) - czhls_in(kso_in-1))
              zwgmax = zporv(i,j)*(czhls_in(kso_in) - czhls_in(kso_in-1))

              w_so_gl_int(i,j,kso_in) = MAX( zwgmin, MIN(zwgmax, w_so_gl(i,j,kso_in)*zporv(i,j)))

              ! in between two input layers
              ! Version 1.22: corrected the interpolation indices
              IF (czmls_lm(kso) >= czmls_in(kso_in-1) .AND. czmls_lm(kso) < czmls_in(kso_in)) THEN
                IF (kso_in == 2) THEN  ! uppermost input layer
                  zwg_fr(i,j,kso) = (w_so_gl_int(i,j,kso_in-1)/(czhls_in(kso_in-1)) * &
                                   (czmls_in(kso_in)-czmls_lm(kso))  + &
                                   w_so_gl_int(i,j,kso_in) / (czhls_in(kso_in)-czhls_in(kso_in-1)) * &
                                   (czmls_lm(kso)-czmls_in(kso_in-1)) )  / (czmls_in(kso_in)-czmls_in(kso_in-1))
                ELSE
                  zwg_fr(i,j,kso) = (w_so_gl_int(i,j,kso_in-1)/(czhls_in(kso_in-1)-czhls_in(kso_in-2)) * &
                                   (czmls_in(kso_in)-czmls_lm(kso))  + &
                                   w_so_gl_int(i,j,kso_in) / (czhls_in(kso_in)-czhls_in(kso_in-1)) * &
                                   (czmls_lm(kso)-czmls_in(kso_in-1)) )  / (czmls_in(kso_in)-czmls_in(kso_in-1))
                ENDIF
                w_so_lm(i,j,kso) = zwg_fr(i,j,kso) * xdzhs(kso)
              ENDIF

            ENDDO ! kso_in

! TEST: falls mit SWI:
            ! back transformation from soil wetness index
!            w_so_lm(i,j,kso) = (w_so_lm(i,j,kso)+zpwp(i,j)*plcov_lm(i,j)/(zfcap(i,j)-zpwp(i,j)))&
!                             /(plcov_lm(i,j)/(zfcap(i,j)-zpwp(i,j))+(1-plcov_lm(i,j))/zfcap(i,j))
            !w_so_lm(i,j,kso) = w_so_lm(i,j,kso)*(zfcap(i,j)-zpwp(i,j)) + zpwp(i,j)

          END IF  ! czmls_lm(kso) > czmls_in(ke_soil_in+1)

        ENDDO  ! kso

      END IF ! landmask
    ENDDO
  ENDDO


! --- soil temperature

  ! heat capacity of input layers
  DO   kso_in = 1,ke_soil_in+1
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (fr_land_lm(i,j) >= 0.5_ireals) THEN   ! for land-points only
          zwg_fr_in(i,j,kso_in) = w_so_gl(i,j,kso_in)*zporv(i,j)/(czhls_in(kso_in)-czhls_in(kso_in-1))
          zroc_in(i,j,kso_in) = zrocg(i,j) + rho_w * zwg_fr_in(i,j,kso_in) * chc_w
        END IF
      ENDDO
    ENDDO
  ENDDO

  DO j = 1, je2lm
    DO i = 1, ie2lm
      IF (fr_land_lm(i,j) >= 0.5_ireals) THEN   ! for land-points only

        ! uppermost layer
        IF (czmls_lm(1) <= czmls_in(1)) THEN
              zroc(i,j,1) = zroc_in(i,j,1)
              t_so_lm(i,j,1) = t_so_int(i,j,0) + (t_so_int(i,j,1)-t_so_int(i,j,0)) * zmls(1)/czmls_in(1)
        ELSE
              zroc(i,j,1) = (zroc_in(i,j,1) * czmls_in(1) &
                                     + zroc_in(i,j,2) * (czmls_lm(1)-czmls_in(1)) ) &
                                     / czmls_lm(1)

              ztzhlto = t_so_int(i,j,1) - czmls_in(1) / czmls_in(1) &
                        * (t_so_int(i,j,1) - t_so_int(i,j,0))
              ztzhltu = t_so_int(i,j,1) + (czmls_lm(1)-czmls_in(1))  &
                        / (czmls_in(2) - czmls_in(1) )  &
                        * (t_so_int(i,j,2) - t_so_int(i,j,1))

              zeo = zroc_in(i,j,1) * czmls_in(1)  &
                        * 0.5 * ( t_so_int(i,j,1) + ztzhlto )
              zeu = zroc_in(i,j,2) * (czmls_lm(1)-czmls_in(1))   &
                        * 0.5 * ( t_so_int(i,j,1) + ztzhltu )

              t_so_lm(i,j,1) = (zeo + zeu) / ( zroc(i,j,1) * czmls_lm(1) )
        ENDIF

        ! other layers
        DO kso = 2, ke_soil_lm+1
           ! below all input layers
           IF (czmls_lm(kso) > czmls_in(ke_soil_in+1)) THEN
              t_so_lm(i,j,kso) = t_cl_lm(i,j)
           ELSE
              DO kso_in = 1, ke_soil_in

                ! first LM layer inside first INPUT layer
                IF (czmls_lm(kso) <= czmls_in(1)) THEN  

                    zroc(i,j,kso) = zroc_in(i,j,1)

                    t_so_lm(i,j,kso) = t_so_int(i,j,0) + (t_so_int(i,j,1)-t_so_int(i,j,0)) * zmls(1)/czmls_in(1)

                ! INPUT layer boundary inside LM layer
                ELSEIF (czmls_lm(kso) > czmls_in(kso_in) .AND. czmls_lm(kso-1) <= czmls_in(kso_in) ) THEN
 
                    zroc(i,j,kso) = (zroc_in(i,j,kso_in) * (czmls_in(kso_in)-czmls_lm(kso-1)) &
                                     + zroc_in(i,j,kso_in+1) * (czmls_lm(kso)-czmls_in(kso_in)) ) &
                                     / (czmls_lm(kso)-czmls_lm(kso-1))

                    ztzhlto = t_so_int(i,j,kso_in) - (czmls_in(kso_in)-czmls_lm(kso-1)) / czmls_in(kso_in) &
                              * (t_so_int(i,j,kso_in) - t_so_int(i,j,kso_in-1))
                    ztzhltu = t_so_int(i,j,kso_in) + (czmls_lm(kso)-czmls_in(kso_in))  &
                              / (czmls_in(kso_in+1) - czmls_in(kso_in) )  &
                              * (t_so_int(i,j,kso_in+1) - t_so_int(i,j,kso_in))

                    zeo = zroc_in(i,j,kso_in) * (czmls_in(kso_in)-czmls_lm(kso-1))  &
                          * 0.5 * ( t_so_int(i,j,kso_in) + ztzhlto )
                    zeu = zroc_in(i,j,kso_in+1) * (czmls_lm(kso)-czmls_in(kso_in))   &
                          * 0.5 * ( t_so_int(i,j,kso_in) + ztzhltu )

                    t_so_lm(i,j,kso) = (zeo + zeu) / ( zroc(i,j,kso) * (czmls_lm(kso) - czmls_lm(kso-1)) )

                ! LM layer inside a further down INPUT layer
                ELSEIF (czmls_lm(kso) <= czmls_in(kso_in+1) .AND. czmls_lm(kso-1) >= czmls_in(kso_in)) THEN

                    zroc(i,j,kso) = zroc_in(i,j,kso_in+1)

                    t_so_lm(i,j,kso) = t_so_int(i,j,kso_in) &
                                       + (t_so_int(i,j,kso_in+1) - t_so_int(i,j,kso_in)) &
                                       * (czmls_lm(kso) - czmls_in(kso_in)) / (czmls_in(kso_in+1) - czmls_in(kso_in)) 


                END IF
              ENDDO
           END IF
        ENDDO ! kso
      END IF ! landmask
    ENDDO
  ENDDO

  ! Initialize the fresh-snow indicator with 1
  ! (means all the snow is fresh)
  WHERE (fr_land_lm(:,:) >= 0.5_ireals)
    freshsnw_lm(:,:) = 1.0_ireals
  END WHERE


!---------------------------------------------------------------------
! End of Subroutine init_multi_layer_gme_ml
!---------------------------------------------------------------------

END SUBROUTINE init_multi_layer_gme_ml

!<-- SB

!==============================================================================

SUBROUTINE month2hour( yactdate, m1, m2, pw1, pw2, idbg)

!------------------------------------------------------------------------------
!
!  Find the 2 nearest months m1, m2 and the weights pw1, pw2 to the actual 
!  date and time to interpolate data to the current hour from data valid in 
!  the middle of the months.
!
!------------------------------------------------------------------------------

CHARACTER(LEN=14),         INTENT(IN)  :: yactdate
INTEGER (KIND=iintegers),  INTENT(IN)  :: idbg       ! for debug output
INTEGER (KIND=iintegers),  INTENT(OUT) :: m1, m2     ! indices of nearest months
REAL    (KIND=ireals),     INTENT(OUT) :: pw1, pw2   ! weights of nearest months

!=======================================================================

INTEGER (KIND=iintegers) ::        &
  month_days(12),         & ! number of days for each month
  mmon, mday, mhour, myy, & ! month, day, hour and year of actual date
  mdayhour,               & ! actual date (in hours of month)
  mmidthhours,            & ! midth of month in hours
  i, ip1,                 & ! month indices (ip1=i+1)
  mleapy                    ! leap year (1=yes, 0=no)

REAL    (KIND=ireals)    ::        &
  zhour,      &  ! hour
  zdiff(12),  &  ! difference between midth of following months in days
  zhalf(12),  &  ! number of days for half month
  zact           ! actual time in hours of month

! Statementfunction for leap year determination (1=yes, 0=no)
! Same as in date_time.f90
      mleapy(myy) =  MAX(1-MODULO(myy,4)  ,0)     &
                   - MAX(1-MODULO(myy,100),0)     &
                   + MAX(1-MODULO(myy,400),0)

! Number of days for each month
      month_days = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

!=======================================================================

 
  ! Get values of day, month and hour out of actual date
  READ(yactdate(1:10),'(I4,3I2)') myy, mmon, mday, mhour

  ! Compute half of each month (in days)
  ! Leap Year ??
  month_days (2) = 28 + mleapy(myy)
  IF ( idbg  > 2 ) PRINT *, ' Length of February:', month_days(2)
  zhalf(:) = 0.5_ireals * REAL (month_days(:))

  ! Compute difference between the midth of actual month and the
  ! following one (in days)
  DO i = 1,12
    ip1 = MOD(i,12)+1
    zdiff(i) = zhalf(i) + zhalf(ip1)
  ENDDO

  ! Compute actual date (day and hours) and midth of actual month in hours
  mdayhour = (mday-1)*24 + mhour
  mmidthhours = NINT( zhalf(mmon)*24._ireals)
  IF ( idbg > 10 ) THEN
    PRINT *, ' *month2hour* ydate, month, days in this month, hours this month, mmidthhour:', &
               yactdate, mmon, month_days(mmon), mdayhour, mmidthhours
  ENDIF

! Determine the months needed for interpolation of current values
! Search for the position of date in relation to first of month.
! The original data are valid for the mid-month.
!
! EXAMPLE 1
!        March    !  April     !   May             X : aerosol data
!       ----X-----!-----X----o-!-----X-----        ! : first of month
!                       !    ^       !             o : current date
!                       !    ^ interpolation for that point in time
!                       !  zdiff(4)  !
!                       !zact!
!
! EXAMPLE 2
!        March    !  April     !   May             X : ndvi_ratio
!       ----X-----!-----X------!----oX-----        ! : first of month
!                       !           ^              o : current date
!                       !      interpolation for that point in time
!                       !zhalf !
!                       !  zdiff(4)  !
!                       !   zact    !
!
!

  IF( mdayhour < mmidthhours) THEN
    ! point is in first half of month (EXAMPLE 2)
    m1 = mmon - 1
    IF(mmon == 1) m1 = 12
    zact   = zhalf(m1) + REAL(mdayhour)/24._ireals
  ELSE
    ! point is in second half of month (EXAMPLE 1)
    m1 = mmon
    zact   = REAL(mdayhour-mmidthhours)/24._ireals
  ENDIF
  m2 = mod(m1,12) + 1
  pw2 =  zact / zdiff(m1)
  pw1 =  1._ireals - pw2

  IF ( idbg > 10 ) THEN
     PRINT *, '  m1, m2, pw1, pw2:', m1, m2, pw1, pw2
  ENDIF

END SUBROUTINE month2hour

!==============================================================================

END MODULE src_2d_fields
