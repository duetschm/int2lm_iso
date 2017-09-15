!+ Declaration of the COSMO-Model fields.
!==============================================================================

MODULE data_fields_lm

!==============================================================================
!
! Description:
!  This module declares all fields and the vertical coordinate of the 
!  COSMO-Model as allocatable arrays.
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
!  Introduced new fields for evergreen and deciduous forest (for_e, for_d)
! 1.3        2005/12/12 Ulrich Schaettler
!  Added field rho_snow for prognostic treatment of snow density in the LM
! V1_5         2007/07/09 Ulrich Schaettler
!  Added fields qr_lm, qs_lm, qg_lm for interpolating in LM2LM
!  Added fields fr_lake_lm, t_skin_gl (CLM)
!  Added chemistry fields cgas_lm, caero_lm
! V1_7         2007/11/26 Ulrich Schaettler
!  Added field for actual ndvi ratios on COSMO-Model grid
!  Corrected dimension for w_so_lm
! V1_8         2008/05/29 Ulrich Schaettler
!  Renamed all gz0-variables to z0 and only work with z0-values
!  Added field rootdp_mx to store external root depth
!  Added fields for subgrid scale orography and topographical corrections
! V1_9         2009/09/03 Ulrich Schaettler, Burkhardt Rockel, Daniel Luethi
!  Added fields for additional external parameters and output fields:
!  ndvi ratios, aerosol values (12 monthly and actual fields)
!  surface emissivity and stomata resistance of plants
!  Added field for depth_lk_lm  (Burkhardt Rockel)
!  Added fields for 12-monthly values for vegetation (Daniel Luethi)
! V1_10        2009/12/17 Ulrich Schaettler, Jan-Peter Schulz
!  Introduced fields hsurf_gl, p_lm (needed for lum2lm)
!  Introduced prognostic FLake variables for cold start of FLake
!  Introduced sea-ice variables for COSMO-Model (JPS)
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler, Burkhardt Rockel
!  Added new fieldsfor t_2m_gl, qv_2m_gl to process data from JMA
!  Added new field salt_lk_lm for lake water salinity (BR)
! V1_17        2011/03/11 Ulrich Schaettler
!  Added field fr_urban_lm for urban fraction data (K. Trusilova)
! V1_19        2012/06/06 Burkhardt Rockel, Juergen Helmert, Susanne Brienen
!  Introduction of prescribed surface albedo fields: alb_dry, alb_sat,
!     alb_dif, alb_dif12
!  Added new field w_so_gl, to save horizontally interpolated soil moisture fields
!  for an additional vertical interpolation.
!  Added fields for water isotope simulation (Stephan Pfahl)
!  Added fields for height correction of soil isotope values. Hui Tang 2013-11-20
!
! Code Description:
! Language:           Fortran 90.
! Software Standards: "European Standards for Writing and
!                      Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters , ONLY :   &
  ireals,    & ! KIND-type parameters for real variables
  iintegers    ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! 1. Fields for the COSMO-Model
! -----------------------------

  ! External parameters for the COSMO-Model
  REAL (KIND=ireals), TARGET, ALLOCATABLE :: &
    fis_lm     (:,:), & ! orography * G                                 (m2/s2)
    hsurf_lm   (:,:), & ! orography                                     (  m  )
    fr_land_lm (:,:), & ! land fraction of grid element                 (  1  )
    fr_lake_lm (:,:), & ! lake fraction of grid element                 (  1  )
    depth_lk_lm(:,:), & ! lake depth                                    (  m  )
    salt_lk_lm (:,:), & ! lake salinity                                 ( g/kg)
    z0_lm      (:,:), & ! roughness length                              (  m  )
    soiltyp_lm (:,:), & ! type of the soil (keys 0-9)                   (  1  )
    plcov_lm   (:,:), & ! fraction covered by plants                    (  1  )
    plcov_mx_lm(:,:), & ! plant cover during vegetation time            (  1  )
    plcov_mn_lm(:,:), & ! plant cover during time of rest               (  1  )
    lai_lm     (:,:), & ! leaf area index                               (  1  )
    lai_mx_lm  (:,:), & ! leaf area index during vegetation time        (  1  )
    lai_mn_lm  (:,:), & ! leaf area index during time of rest           (  1  )
    rootdp_lm  (:,:), & ! depth of the roots                            (  m  )
    rootdp_mx  (:,:), & ! depth of the roots from external parameters   (  m  )
    for_e_lm   (:,:), & ! ground fraction covered by evergreen forest   (  -  )
    for_d_lm   (:,:), & ! ground fraction covered by deciduous forest   (  -  )
    fr_urban_lm(:,:), & ! ground fraction covered by urban land         (  -  )
    sso_stdh_lm(:,:), & ! standard deviation of subgrid scale orography (  m  )
    sso_gamma_lm(:,:),& ! anisotropy of the orography                   (  -  )
    sso_theta_lm(:,:),& ! angle betw. principal axis of orography and E ( rad )
    sso_sigma_lm(:,:),& ! mean slope of subgrid scale orography         (  -  )
    skyview_lm (:,:), & ! sky view                                      (  1  )
    slo_asp_lm (:,:), & ! slope aspect                                  ( rad )
    slo_ang_lm (:,:), & ! slope angle                                   ( rad )
    horizon_lm(:,:,:),& ! horizon                                       ( rad )
    ps_lm      (:,:)    ! surface pressure                              ( Pa  )

  REAL (KIND=ireals), TARGET, ALLOCATABLE :: &
    alb_dry_lm (:,:), & ! surface albedo field for dry soil             (  1  )
    alb_sat_lm (:,:), & ! surface albedo field for saturated soil       (  1  )
    alb_dif_lm (:,:), & ! solar surface albedo - diffuse                (  1  )
    alb_dif12_lm(:,:,:),&!solar surface albedo - diffuse                (  1  )
    emis_rad_lm(:,:), & ! thermal radiative surface emissivity          (  1  )
    prs_min_lm (:,:), & ! minimum stomata resistance of plants          ( s/m )
    ndvi_mrat_lm(:,:,:),&!ratio of monthly mean normalized differential (  1  )
                        ! vegetation index to annual maximum for 12 months
    ndviratio_lm(:,:),& ! actual value of ndvi for a special day        (  1  )
    ! monthly mean for aerosols
    aer_su12_lm(:,:,:),&! Tegen (1997) aerosol type sulfate drops       (  -  )
    aer_du12_lm(:,:,:),&! Tegen (1997) aerosol type mineral dust coarse (  -  )
    aer_or12_lm(:,:,:),&! Tegen (1997) aerosol type organic(water solub)(  -  )
    aer_bc12_lm(:,:,:),&! Tegen (1997) aerosol type black carbon        (  -  )
    aer_ss12_lm(:,:,:),&! Tegen (1997) aerosol type sea salt            (  -  )
    ! and the actual values for a special day
    aer_su_lm  (:,:), & ! Tegen (1997) aerosol type sulfate drops       (  -  )
    aer_du_lm  (:,:), & ! Tegen (1997) aerosol type mineral dust coarse (  -  )
    aer_or_lm  (:,:), & ! Tegen (1997) aerosol type organic(water solub)(  -  )
    aer_bc_lm  (:,:), & ! Tegen (1997) aerosol type black carbon        (  -  )
    aer_ss_lm  (:,:)    ! Tegen (1997) aerosol type sea salt            (  -  )

  ! additional parameters for externally prescribed yearly
  ! cycle of lai and plcov read from external parameter file  
  REAL (KIND=ireals), TARGET, ALLOCATABLE :: &
    plcov12    (:,:,:), & ! monthly climatology of plant cover          (  1  )
    z012       (:,:,:), & ! monthly climatology of roughness length     (  m  )
    lai12      (:,:,:)    ! monthly climatology of leaf area index      (  1  )

  ! 2-dimensional fields for the COSMO-Model
  REAL (KIND=ireals), TARGET, ALLOCATABLE :: &
    t_s_lm     (:,:), & ! temperature of the ground surface             (  K  )
    t_snow_lm  (:,:), & ! temperature of the snow surface               (  K  )
    t_m_lm     (:,:), & ! temp. between upper and medium soil layer     (  K  )
    t_cl_lm    (:,:), & ! temp. between medium and lower soil layer     (  K  )
    qv_s_lm    (:,:), & ! specific water vapor content on the surface   (kg/kg)
    w_snow_lm  (:,:), & ! water content of the snow                     (m H2O)
    w_i_lm     (:,:), & ! water content of the interception storage     (m H2O)
    w_g1_lm    (:,:), & ! water content of the upper soil layer         (m H2O)
    w_g2_lm    (:,:), & ! water content of the medium soil layer        (m H2O)
    w_g3_lm    (:,:), & ! water content of the lower soil layer         (m H2O)
                        ! (if nl_soil_lm = 3, unused otherwise)
    w_cl_lm    (:,:), & ! climatological deep soil water content        (m H2O)
    freshsnw_lm(:,:), & ! weighting function indicating 'freshness' of snow
    rho_snow_lm(:,:), & ! for prognostic treatment of snow density      (kg/m3)
    hmo3_lm    (:,:), & ! height of maximum ozone concentration         ( Pa  )
    vio3_lm    (:,:), & ! total vertically integrated ozone content     (Pa O3)
    latlm_m    (:,:), & ! latitudes of the COSMO-Model grid points
    lonlm_m    (:,:), & ! longitudes of the COSMO-Model grid points
    latlm_u    (:,:), & ! latitudes of the COSMO-Model u grid points
    lonlm_u    (:,:), & ! longitudes of the COSMO-Model u grid points
    latlm_v    (:,:), & ! latitudes of the COSMO-Model v grid points
    lonlm_v    (:,:), & ! longitudes of the COSMO-Model v grid points
! iso code
    riso_lm(:,:,:,:), & ! isotope ratios in water vapor;
                        ! last index counts variables (1: 18O; 2: 2H);
                        ! unit: 1/SMOW
    risosoil_lm (:,:,:), & ! isotope ratio in soil moisture;
                           ! last index counts variables (1: 18O; 2: 2H);
                           ! unit: 1/SMOW
! Hui Tang 2013-11-20
    drisoke_gl (:,:,:)     ! riso_lm(:,:,ke,1-2) - risosoil_lm (:,:,1-2) (1: 18O; 2: 2H);
! end iso code

  ! 2-dimensional fields for the FLake model
  REAL (KIND=ireals), TARGET, ALLOCATABLE :: &
    t_mnw_lk_lm(:,:), & ! mean temperature of the water column          (  K  )
    t_wml_lk_lm(:,:), & ! mixed-layer temperature                       (  K  )
    t_bot_lk_lm(:,:), & ! temperature at the water-bottom sediment
                        ! interface                                     (  K  )
    t_b1_lk_lm (:,:), & ! temperature at the bottom of the upper layer
                        ! of the sediments                              (  K  )
    c_t_lk_lm  (:,:), & ! shape factor with respect to the
                        ! temperature profile in lake thermocline       (  -  )
    h_ml_lk_lm (:,:), & ! thickness of the mixed-layer                  (  m  )
    h_b1_lk_lm (:,:), & ! thickness of the upper layer
                        ! of bottom sediments                           (  m  )
    t_ice_lm   (:,:), & ! temperature of ice/water surface              (  K  )
    h_ice_lm   (:,:)    ! lake/sea ice thickness                        (  m  )

  LOGICAL, ALLOCATABLE :: &
    lolp_lm (:,:),    & ! Land Sea Mask of COSMO for 'M'atch Interpolation
    lmask_lm(:,:)       ! mask of points on the frame
 
  ! 2-dimensional fields for horizontal interpolated fields
  REAL (KIND=ireals), TARGET, ALLOCATABLE :: &
    ndvi_ar_gl (:,:), & ! actual ndvi ratio                             (  1  )
    hsurf_gl   (:,:), & ! interpolated orography (for lum2lm)           (  m  )
    hsurfs_gl(:,:,:), & ! interpolated splitted parts of coarse topo    (  m  )
    fis_gl     (:,:), & ! GME interpolated orography * G                (m2/s2)
    ps_gl      (:,:), & ! surface pressure on the interpol. GME orogr.  ( Pa  )
    dpsdt_gl   (:,:), & ! surface pressure tendency                     (Pa/s )
    t_s_gl     (:,:), & ! temperature of the ground surface             (  K  )
    t_skin_gl  (:,:), & ! skin temperature of the ground surface        (  K  )
    t_2m_gl    (:,:), & ! 2m temperature                                (  K  )
    qv_2m_gl   (:,:), & ! 2m humidity                                   (kg/kg)
    fic_gl     (:,:), & ! check level of geopotential                   (m2/s2)
    rh_s_gl    (:,:), & ! relative humidity at the surface              (kg/kg)
    dtms_gl    (:,:), & ! t_m_lm    - t_s_lm                            (  K  )
    dtkes_gl   (:,:), & ! t(ke)_lm  - t_s_lm                            (  K  )
    dtssnow_gl (:,:), & ! t_s_lm    - t_snow_lm                         (  K  )
    hhl_gl   (:,:,:), & ! height of half-levels on the interpolated 
                        ! COARSE COSMO-Model orography                  (  m  )
    dp0_gl   (:,:,:), & ! reference pressure thickness of layers        ( Pa  )
    rho0_gl  (:,:,:), & ! reference density at the full model levels    (kg/m3)
    p0_gl    (:,:,:)    ! reference pressure on full levels and 
                        ! interpolated COARSE COSMO-Model orography     ( Pa  )

  ! 3-dimensional fields for the COSMO-Model (u_lm, v_lm and t_lm are also
  ! used for horizontal interpolated fields)
  REAL (KIND=ireals), TARGET, ALLOCATABLE :: &
    u_lm     (:,:,:), & ! zonal wind speed                              ( m/s )
    v_lm     (:,:,:), & ! meridional wind speed                         ( m/s )
    w_lm     (:,:,:), & ! vertical wind speed (defined on half levels)  ( m/s )
    t_lm     (:,:,:), & ! temperature                                   (  K  )
    qv_lm    (:,:,:), & ! specific water vapor content                  (kg/kg)
    qc_lm    (:,:,:), & ! specific cloud water content                  (kg/kg)
    qi_lm    (:,:,:), & ! cloud ice content                             (kg/kg)
    qr_lm    (:,:,:), & ! rain content                                  (kg/kg)
    qs_lm    (:,:,:), & ! snow content                                  (kg/kg)
    qg_lm    (:,:,:), & ! graupel content                               (kg/kg)
    grh_lm   (:,:,:), & ! generalized relative humidity                 (kg/kg)
    p0_lm    (:,:,:), & ! reference pressure                            ( Pa  )
    dp0_lm   (:,:,:), & ! reference pressure thickness of layers        ( Pa  )
    rho0_lm  (:,:,:), & ! reference density at the full model levels    (kg/m3)
    p_lm     (:,:,:), & ! full pressure (needed for lum2lm)             ( Pa  )
    pp_lm    (:,:,:), & ! deviation from the reference pressure         ( Pa  )
    t0_lm    (:,:,:), & ! reference temperature                         (  K  )
    hhl_lm   (:,:,:), & ! height of half-levels of COSMO-Model          (  m  )
    t_so_lm  (:,:,:), & ! multi-layer soil temperature                  (  K  )
    dt_so_gl (:,:,:), & ! multi-layer soil temperature (for interpolation)
    w_so_lm  (:,:,:), & ! multi-layer soil moisture                     (m H2O)
    w_so_gl  (:,:,:)    ! multi-layer soil moisture (for interpolation) (m H2O)

!------------------------------------------------------------------------------

! 2. Grid correspondence GME-COSMO-Model
! --------------------------------------

! The following arrays are needed for the correspondence of the COSMO-grid
! and the GME-grid: indices, barycentric coordinates.
! These are computed in SUBROUTINE coor_gmlm
!
! All following arrays have three "subnames":
!    '_m' : for the mass points of the COSMO-grid
!    '_u' : for the u-   points of the COSMO-grid (Arakawa C-grid)
!    '_v' : for the v-   points of the COSMO-grid (Arakawa C-grid)
  
  INTEGER (KIND=iintegers), ALLOCATABLE ::    &
    index_m(:,:,:),       & !
    index_u(:,:,:),       & !
    index_v(:,:,:)          !
  
  ! index(ie2lm,je2lm, 4): array for interpolation to point (kgrid)
  ! kindex(kgrid, 1) : j1 (first  dim.) of nearest GME grid point
  ! kindex(kgrid, 2) : j2 (second dim.) of nearest GME grid point
  ! kindex(kgrid, 3) : jd (diamond) of nearest GME grid point
  ! kindex(kgrid, 4) : index of triangle containing point (kgrid)
  !
  !INTEGER (KIND=iintegers), ALLOCATABLE ::    &
  !  idxrp_m(:,:),         & !
  !  idxrp_u(:,:),         & !
  !  idxrp_v(:,:)            !
  !
  ! idxrp(ie2lm*je2lm):
  ! index of the points in the grid, which are owned
  ! by the local processor. The grid is thought as a
  ! linear array according to the FORTRAN rules, i.e. point (mgrid)
  ! gives the index kgrid
  
  REAL    (KIND=ireals),    ALLOCATABLE ::    &
    baryll_m(:,:,:),      & ! baryll(ie2lm,je2lm, 2):
    baryll_u(:,:,:),      & !   barycentric coordinates alpha and beta to
    baryll_v(:,:,:)         !   interpolate to point (kgrid)
  
  REAL    (KIND=ireals),    ALLOCATABLE ::    &
    rotang_m(:,:,:),      & ! rotang(ie2lm,je2lm, 2):
    rotang_u(:,:,:),      & !   sine and cosine of psi, the angle for
    rotang_v(:,:,:)         !   wind vector rotation
  
  REAL    (KIND=ireals),    ALLOCATABLE ::    &
    w_intpol(:,:,:)         ! interpolation weights for gme2lm

  INTEGER (KIND=iintegers), ALLOCATABLE ::    &
    n_intpol(:,:,:), & ! nearest GME gridpoint for nearest neighbor interpolation
    m_intpol(:,:,:)    ! nearest GME gridpoint with same lsm for match interpolation

  LOGICAL                 , ALLOCATABLE ::    &
    l_intpol(:,:)      ! to use a far away GME gridpoint with same lsm

  ! Global Scalars:
  INTEGER (KIND=iintegers) :: rpoints_m, rpoints_u, rpoints_v
    ! number of points in the grid, which are owned by the local processor.
  
!------------------------------------------------------------------------------

! 3. Grid correspondence: regular input grid - COSMO-Model
! --------------------------------------------------------

! The following arrays are needed for the correspondence of the COSMO-Model-grid
! and the input grid (for EC or COSMO): i- and j-indices and relative distances
! of grid points. The third dimension of all arrays is 5:
!   1: for mass grid points
!   2: for COSMO u-grid points and coarse mesh u-grid points
!   3: for COSMO u-grid points and coarse mesh v-grid points
!   4: for COSMO v-grid points and coarse mesh v-grid points
!   5: for COSMO v-grid points and coarse mesh u-grid points
! (For a non-staggered grid (like ECMWF has), 2 and 3 are the same and
!  4 and 5 are the same).

  INTEGER (KIND=iintegers), ALLOCATABLE   :: &
    i_index(:,:,:),  & ! i-index of coarse mesh grid point which is lower left
                       ! to a given COSMO (fine mesh) grid point
    j_index(:,:,:)     ! j-index of coarse mesh grid point which is lower left
                       ! to a given COSMO (fine mesh) grid point

  REAL (KIND=ireals), ALLOCATABLE         :: &
    x_wght (:,:,:),  & ! relative distance between x- (i-) coarse mesh and
                       ! fine mesh grid points
    y_wght (:,:,:)     ! relative distance between y- (j-) coarse mesh and
                       ! fine mesh grid points

!------------------------------------------------------------------------------

! 4. Fields for chemistry output
! ------------------------------

  REAL (KIND=ireals), TARGET, ALLOCATABLE ::           &
    cgas_lm  (:,:,:,:),   & ! 
    caero_lm (:,:,:,:)      !

!------------------------------------------------------------------------------
!
!- End of module header
!==============================================================================

END MODULE data_fields_lm
