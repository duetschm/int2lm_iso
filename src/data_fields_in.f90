!+ Declaration of the coarse grid input fields.
!==============================================================================

MODULE data_fields_in

!==============================================================================
!
! Description:
!  This module declares all fields and the vertical coordinates of the
!  input model(s) as allocatable arrays.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  uschaettler@dwd.d400.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2005/04/11 Ulrich Schaettler
!  Initial release for INT2LM
! 1.3        2005/12/12 Ulrich Schaettler
!  Added field rho_snow for prognostic treatment of snow density in the LM
! V1_5         2007/07/09 Ulrich Schaettler
!  Added fields qr_in, qs_in, qg_in for interpolating in LM2LM
!  Added field hsurfs_in (MCH, for SLEVE coordinate)
!  Added fields sst_in, t_skin_in, w_so_rel_in (CLM)
!  Added chemistry fields cgas_in, caero_in
! V1_7         2007/11/26 Ulrich Schaettler
!  Added fields for monthly and actual ndvi ratios on GME grid
!  Added additional fields for vegetation and rest of plcov, lai for input model
! V1_8         2008/05/29 Ulrich Schaettler
!  Renamed all gz0-variables to z0 and only work with z0-values
! V1_9         2009/09/03 Ulrich Schaettler
!  Removed fields for monthly and actual ndvi ratios on GME grid: are read
!  now from external COSMO parameters
! V1_10        2009/12/17 Oliver Fuhrer, Jan-Peter Schulz
!  Added field p_in (to read pressure or pressure deviation; also needed for um2lm)
!  Added sea-ice fields for input model (JPS)
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler
!  New fields for qv_2m_in, t_2m_in (from GSM-JMA) and rlai_in
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters, ONLY:     &
    ireals,    & ! KIND-type parameters for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! 1. Arrays for the GME as input model
! ------------------------------------

  ! Grid correspondance
  REAL (KIND=ireals), ALLOCATABLE ::           &
    grd_glob  (:,:,:,:), & ! coefficients needed for the calculation of the 
                           ! gradients in eta- and chi-direction          (1/m)
    eta_glob  (:,:,:)  , & ! 
    chi_glob  (:,:,:)  , & !
    cpsi_glob (:,:,:)  , & ! cosine of the rotation angle between the local
                           ! coordinate system of the home node and the 6 (5)
                           ! neighbouring gridpoints
    spsi_glob (:,:,:)      ! sine of the rotation angle between the local 
                           ! coordinate system of the home node and the 6 (5)
                           ! neighbouring gridpoints

  ! GME fields
  REAL (KIND=ireals), TARGET, ALLOCATABLE ::           &
    fis_gme     (:,:,:),  & ! orography * g                             (m2/s2)
    soiltyp_gme (:,:,:),  & ! type of the soil (keys 0-9)                 --
    ps_gme      (:,:,:),  & ! surface pressure                          ( Pa  )
    t_s_gme     (:,:,:),  & ! temperature of the ground surface         (  K  )
    t_snow_gme  (:,:,:),  & ! temperature of the snow-surface           (  K  )
    t_m_gme     (:,:,:),  & ! temp. between upper and medium soil layer (  K  )
    qv_s_gme    (:,:,:),  & ! specific water vapor content on surface   (kg/kg)
    w_g1_gme    (:,:,:),  & ! water content of the upper soil layer     (m H2O)
    w_g2_gme    (:,:,:),  & ! water content of the medium soil layer    (m H2O)
    w_g3_gme    (:,:,:),  & ! water content of the deepest soil layer   (m H2O)
    t_i3e_gme   (:,:,:),  & ! temperature lowest layer                  (  K  )
    grh_gme     (:,:,:)     ! generalized relative humidity at one level(  %  )

  LOGICAL, ALLOCATABLE            ::           &
    lolp_gme    (:,:,:)     ! Land Sea Mask of GME for 'M'atch Interpolation

  REAL (KIND=ireals), TARGET, ALLOCATABLE ::           &

    u_gme    (:,:,:,:),   & ! zonal wind speed                          ( m/s )
    v_gme    (:,:,:,:),   & ! meridional wind speed                     ( m/s )
    t_gme    (:,:,:,:),   & ! temperature                               (  K  )
    qv_gme   (:,:,:,:),   & ! specific water vapor content              (kg/kg)
    qc_gme   (:,:,:,:),   & ! specific cloud water content              (kg/kg)
    qi_gme   (:,:,:,:)      ! cloud ice content                         (kg/kg)

  ! and all that stuff that is perhaps not needed
  REAL (KIND=ireals), TARGET, ALLOCATABLE ::           &
    z0_gme      (:,:,:),  & ! surface roughness                         (  m  )
    fr_land_gme (:,:,:),  & ! land fraction of grid element             (  1  )
    plcmx_gme   (:,:,:),  & ! vegetation: fraction covered by plants    (  1  )
    plcmn_gme   (:,:,:),  & ! rest:       fraction covered by plants    (  1  )
    rlaimx_gme  (:,:,:),  & ! vegetation: leaf area index               (  1  )
    rlaimn_gme  (:,:,:),  & ! rest:       leaf area index               (  1  )
    root_gme    (:,:,:),  & ! depth of the roots                        (  m  )
    vio3_gme    (:,:,:),  & ! total vertically integrated ozone content (Pa O3)
    hmo3_gme    (:,:,:),  & ! height of maximum ozone concentration     ( Pa  )
    t_g_gme     (:,:,:),  & ! temperature                               (  K  )
    t_cl_gme    (:,:,:),  & ! temp.  between medium and lower soil layer(  K  )
    w_snow_gme  (:,:,:),  & ! water content of the snow                 (m H2O)
    w_i_gme     (:,:,:),  & ! water content of the interception storage (m H2O)
    w_cl_gme    (:,:,:),  & ! climatological deep soil water content    (m H2O)
    t_so_gme    (:,:,:,:),& ! temperature for new multi layer soil model(  K  )
    w_so_gme    (:,:,:,:),& ! soil moisture for multi layer soil model  (m H2O)
    freshsnw_gme(:,:,:),  & ! weighting function indicating 'freshness' of snow
    rho_snow_gme(:,:,:),  & ! for prognostic treatment of snow density  (kg/m3)
    t_ice_gme   (:,:,:),  & ! temperature of sea ice surface            (  K  )
    h_ice_gme   (:,:,:),  & ! sea ice thickness                         (  m  )
    dpsdt_gme   (:,:,:)     ! surface pressure tendency                 (Pa/s )


!------------------------------------------------------------------------------

! 2. Arrays for input models with regular grids
! ---------------------------------------------

  REAL (KIND=ireals), TARGET, ALLOCATABLE ::           &
    fis_in      (:,:),    & ! orography * g                             (m2/s2)
    hsurf_in    (:,:),    & ! orography                                 (  m  )
    hsurfs_in   (:,:,:),  & ! splitted parts of the coarse orography (SLEVE)
    soiltyp_in  (:,:),    & ! type of the soil (keys 0-9)                 --
    ps_in       (:,:),    & ! surface pressure                          ( Pa  )
    sst_in      (:,:),    & ! sea surface temperature                   (  K  )
    t_s_in      (:,:),    & ! temperature of the ground surface         (  K  )
    t_2m_in     (:,:),    & ! 2m temperature                            (  K  )
    t_skin_in   (:,:),    & ! skin temperature of the ground surface    (  K  )
    t_snow_in   (:,:),    & ! temperature of the snow-surface           (  K  )
    t_g1_in     (:,:),    & ! temperature of first soil layer           (  K  )
    t_g2_in     (:,:),    & ! temperature of second soil layer          (  K  )
    t_g3_in     (:,:),    & ! temperature of third soil layer           (  K  )
    qv_s_in     (:,:),    & ! specific water vapor content at surface   (kg/kg)
    qv_2m_in    (:,:),    & ! specific water vapor in 2m                (kg/kg)
    w_g1_in     (:,:),    & ! water content of the upper soil layer     (m H2O)
    w_g2_in     (:,:),    & ! water content of the medium soil layer    (m H2O)
    w_g3_in     (:,:),    & ! water content of the deepest soil layer   (m H2O)
    t_so_in     (:,:,:),  & ! temperature for new multi layer soil model(  K  )
    w_so_in     (:,:,:),  & ! soil moisture for multi layer soil model  (m H2O)
    w_so_rel_in (:,:,:),  & ! multi-layer volumetric soil moisture      (  1  )
    freshsnw_in (:,:),    & ! weighting function indicating 'freshness' of snow
    rho_snow_in (:,:),    & ! for prognostic treatment of snow density  (kg/m3)
    t_ice_in    (:,:),    & ! temperature of sea ice surface            (  K  )
    h_ice_in    (:,:),    & ! sea ice thickness                         (  m  )
    t_ke_in     (:,:),    & ! temperature lowest layer                  (  K  )
    grh_in      (:,:)       ! generalized relative humidity at one level(  %  )

  LOGICAL, ALLOCATABLE            ::           &
    lolp_in     (:,:)       ! Land Sea Mask of input fields for 
                            ! 'M'atch Interpolation

  REAL (KIND=ireals), TARGET, ALLOCATABLE ::           &
    u_in        (:,:,:),  & ! zonal wind speed                          ( m/s )
    v_in        (:,:,:),  & ! meridional wind speed                     ( m/s )
    w_in        (:,:,:),  & ! vertical wind speed                       ( m/s )
    t_in        (:,:,:),  & ! temperature                               (  K  )
    p_in        (:,:,:),  & ! standard pressure                         ( Pa  )
    pp_in       (:,:,:),  & ! deviation from standard pressure          ( Pa  )
    qv_in       (:,:,:),  & ! specific water vapor content              (kg/kg)
    qc_in       (:,:,:),  & ! specific cloud water content              (kg/kg)
    qi_in       (:,:,:),  & ! specific cloud ice content                (kg/kg)
    qr_in       (:,:,:),  & ! specific rain      content                (kg/kg)
    qs_in       (:,:,:),  & ! specific snow      content                (kg/kg)
    qg_in       (:,:,:),  & ! specific graupel   content                (kg/kg)
    hhl_in      (:,:,:),  & ! height of half levels of coarse grid      (  m  )
    p0_in       (:,:,:)     ! reference pressure on coarse grid

  ! and all that stuff that is perhaps not needed
  REAL (KIND=ireals), TARGET, ALLOCATABLE ::           &
    z0_in       (:,:),    & ! surface roughness                         (  m  )
    fr_land_in  (:,:),    & ! land fraction of grid element             (  1  )
    fland_in_tot(:,:),    & ! land fraction for global field            (  1  )
    plcov_in    (:,:),    & ! fraction covered by plants                (  1  )
    plcmx_in    (:,:),    & ! vegetation: fraction covered by plants    (  1  )
    plcmn_in    (:,:),    & ! rest:       fraction covered by plants    (  1  )
    rlai_in     (:,:),    & ! vegetation: leaf area index               (  1  )
    rlaimx_in   (:,:),    & ! vegetation: leaf area index               (  1  )
    rlaimn_in   (:,:),    & ! rest:       leaf area index               (  1  )
    root_in     (:,:),    & ! depth of the roots                        (  m  )
    vio3_in     (:,:),    & ! total vertically integrated ozone content (Pa O3)
    hmo3_in     (:,:),    & ! height of maximum ozone concentration     ( Pa  )
    t_g_in      (:,:),    & ! temperature at the ground surface         (  K  )
    t_m_in      (:,:),    & ! temp. between upper and medium soil layer (  K  )
    t_cl_in     (:,:),    & ! temp. between medium and lower soil layer (  K  )
    w_snow_in   (:,:),    & ! water content of the snow                 (m H2O)
    w_i_in      (:,:),    & ! water content of the interception storage (m H2O)
    w_cl_in     (:,:),    & ! climatological deep soil water content    (m H2O)
    dpsdt_in    (:,:),    & ! surface pressure tendency                 (Pa/s )
    fic_in      (:,:),    & ! control level for geopotential            (m2/s2)
! iso code
    riso_in     (:,:,:,:),& ! isotope ratios in water vapor;
                            ! last index counts variables (1: 18O; 2: 2H);
                            ! unit: 1/SMOW
    risosoil_in (:,:,:)     ! isotope ratio in soil moisture;
                            ! last index counts variables (1: 18O; 2: 2H);
                            ! unit: 1/SMOW
! end iso code

  ! and for the coordinates of fine LM grid points in the coarse grid
  REAL (KIND=ireals), TARGET, ALLOCATABLE ::           &
    lat_coarse_m   (:,:), & ! latitudes of the LM grid points
    lon_coarse_m   (:,:), & ! longitudes of the LM grid points
    lat_coarse_u   (:,:), & ! latitudes of the LM u grid points
    lon_coarse_u   (:,:), & ! longitudes of the LM u grid points
    lat_coarse_v   (:,:), & ! latitudes of the LM v grid points
    lon_coarse_v   (:,:)    ! longitudes of the LM v grid points

!------------------------------------------------------------------------------

! 3. Fields for chemistry input
! -----------------------------

  REAL (KIND=ireals), TARGET, ALLOCATABLE ::           &
    cgas_in  (:,:,:,:),   & ! 
    caero_in (:,:,:,:)      ! 

!- End of module header
!==============================================================================

END MODULE data_fields_in
