!+ External procedure for reading external parameters for LM and GME
!==============================================================================

SUBROUTINE external_data (ierror, yerror)

!==============================================================================
!
! Description:
!   This routine organizes the input of the external parameters for the fine
!   grid LM and the coarse grid fields. It also computes the reference
!   atmosphere of the LM.
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
!  Added fields for_e, for_d for reading from external parameters
! V1_5         2007/07/09 Guy de Morsier
!  Eliminated akhlm, bkhlm
!  Added switch llbc_smooth for a smooth transition from the coarse orography at
!  the lateral boundaries to the fine orography after nlbc_smooth grid points.
! V1_6         2007/09/07 Burkhardt Rockel, Uwe Boehm
!   Added: lakes, climatological deep soil temperature, climatological year
!          cressman scheme
!   Added ltcl_lm to parameter list for lcm2lm
! V1_7         2007/11/26 Ulrich Schaettler
!  Bug correction in the orography filtering when running sequential program
!  Introduced itype_t_cl to choose origin for deep soil temperature
!  Introduced itype_rootdp to choose treatment of root depths
!  Put computation of actual values for external parameters to module
!    src_2d_fields, to do that also for boundaries in case of climate mode
!  In case of interpolating external parameters from coarse grid model the
!    vegetation and rest-data set have to be read from coarse grid model
! V1_8         2008/05/29 Ulrich Schaettler
!  Renamed all gz0-variables to z0 and only work with z0-values
!  Read additional external parameters, if lsso or lradtopo
!  Do not deallocate fis_in for possible later calculations, if fic_in cannot be read
!  Save values for rootdp from external parameters to rootdp_mx
! V1_9         2009/09/03 Guenther Zaengl
!  Adaptations for new reference atmosphere
!  Moved computation of reference atmosphere for coarse COSMO grid (input model)
!  to src_read_coarse_grid.f90
!  Implemented logical flags for reading additional external parameters
! V1_10        2009/12/17 Ulrich Schaettler
!  Modifications to read UM data: interpolate coarse orography to hsurf_gl
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_12        2010/06/14 Burkhardt Rockel
!  In case DEPTH_LK is not in external data file, set depth_lk_lm to default 20 metres
! V1_14        2010/11/19 Ulrich Schaettler, Burkhardt Rockel
!  Modifications to allow processing of JMA and NCEP data
!  Initialization for field soiltyp_in must be done for the whole field
!  Added logical variable for lake salinity: lsalt_lm (BR)
! V1_15        2010/12/10 Ulrich Schaettler, Ulrich Blahak
!  Bug fix for computing the part of the coarse grid covering the fine grid
!  Bug fix when computing nearest coarse grid point for isolated points (UB)
! V1_17        2011/03/11 Ulrich Schaettler
!  Added lurban flag for reading urban fraction data fr_urban_lm (K. Trusilova)
! V1_19        2012/06/06 Ulrich Schaettler, Dmitrii Mironov,
!                         Daniel Luethi, Burkhardt Rockel
!  Added lhir2lm as internal logical flag    (US)
!  Adapted interface to SR distribute_field
!  Set lake depth consistent with fr_lake_lm (DM)
!  Initialize all variables for isolated points with -1 (for undefined)
!  CLM:
!  Introduce new logical variable lcm_hgt_coor for hybrid height coordinates
!    on input and choose proper hsurf/fis interpolation
!  Added prescribed surface albedos and new NL switch itype_albedo
!  Use longitudes_in, latitudes_in to compute isolated point instead of
!    dlon_in, dlat_in to account for non equidistant grid (e.g. gaussian)
! V1_20        2012/09/03 Michael Baldauf, Anne Roches
!  Adaptations to calls to modified SR reference_atmosphere_xx
!  Allocate zhsurfs_lm in all cases, because it is used as argument to
!   subroutine calls (Anne Roches)
! V1_21        2013/03/25 Ulrich Schaettler
!  Initialize error code with 0
! V1_22        2013/07/11 Ulrich Schaettler
!  Removed unused variables from the USE sections
!  Use routines for reference_atmospheres from module vgrid_refatm_utils
!  Introduced namelist variable for name of HHL field
! V1_23        2013/10/02 Ulrich Schaettler
!  Check hsurf and hhl-read only up to epsilon
!  Rename vcoord_out, refatm_out to vcoord, refatm (as is in COSMO-Model)
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
  p0_lm        ,   & ! reference pressure                               ( Pa  )
  dp0_lm       ,   & ! reference pressure thickness of layers           ( Pa  )
  rho0_lm      ,   & ! reference density at the full model levels       (kg/m3)
  t0_lm        ,   & ! reference temperature                            (  K  )
  hhl_lm       ,   & ! height of half-levels of LM                      (  m  )
  fis_lm       ,   & ! orography * G                                    (m2/s2)
  hsurf_lm     ,   & ! orography                                        (  m  )
  fr_land_lm   ,   & ! land fraction of grid element                    (  1  )
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
  for_e_lm     ,   & ! ground fraction covered by evergreen forest      (  1  )
  for_d_lm     ,   & ! ground fraction covered by deciduous forest      (  1  )
  fr_urban_lm  ,   & ! ground fraction covered by urban land            (  1  )
  lolp_lm      ,   & ! Land Sea Mask of LM for 'M'atch Interpolation
  lmask_lm     ,   & ! mask of points on the frame
  fis_gl       ,   & ! GME interpolated orography * G                   (m2/s2)
  fic_gl       ,   & ! check level of geopotential                      (m2/s2)
  hhl_gl       ,   & ! height of half-levels on the interpol. COARSE LM oro.(m)
  hsurf_gl     ,   & ! height of orography interpolated from coarse grid(  m  )
  hsurfs_gl    ,   & ! interpolated splitted parts of coarse topo       (  m  )
  fr_lake_lm   ,   & ! lake fraction of grid element                    (  1  )
  depth_lk_lm  ,   & ! lake depth                                       (  m  )
  salt_lk_lm         ! lake salinity                                    ( g/kg)

USE data_fields_lm, ONLY : &
  index_m      ,   & !
  baryll_m     ,   & !
  w_intpol     ,   & ! interpolation weights for gme2lm                 (  -  )
  n_intpol     ,   & ! nearest GME gridpoint for nearest neighbor interpolation
  m_intpol     ,   & ! nearest GME gridpoint with same lsm for match interpolation
  l_intpol     ,   & ! to use a far away GME gridpoint with same lsm
  latlm_m      ,   & ! latitudes of the LM grid points
  i_index      ,   & !
  j_index      ,   & !
  x_wght       ,   & !
  y_wght             !

!------------------------------------------------------------------------------

USE data_fields_in,  ONLY: &
 fis_gme     ,          & ! orography * G                               (m2/s2)
 soiltyp_gme ,          & ! type of the soil (keys 0-9)                   --
 fr_land_gme ,          & ! land fraction of grid element               (  1  )
 z0_gme      ,          & ! surface roughness                           (  m  )
 root_gme    ,          & ! depth of the roots                          (  m  )
 plcmx_gme   ,          & ! vegetation: fraction covered by plants      (  1  )
 plcmn_gme   ,          & ! rest:       fraction covered by plants      (  1  )
 rlaimx_gme  ,          & ! vegetation: leaf area index                 (  1  )
 rlaimn_gme  ,          & ! rest:       leaf area index                 (  1  )
 lolp_gme    ,          & ! Land Sea Mask of GME for 'M'atch Interpolation
 fis_in      ,          & ! orography * G                               (m2/s2)
 hsurf_in    ,          & ! orography                                   (  m  )
 hsurfs_in   ,          & ! splitted parts of the coarse orography (SLEVE)
 soiltyp_in  ,          & ! type of the soil (keys 0-9)                   --
 lolp_in     ,          & ! Land Sea Mask of input for 'M'atch Interpolation
 z0_in       ,          & ! surface roughness                           (  m  )
 fr_land_in  ,          & ! land fraction of grid element               (  1  )
 fland_in_tot,          & ! land fraction of grid element (global)      (  1  )
 root_in     ,          & ! depth of the roots                          (  m  )
 plcov_in    ,          & ! fraction covered by plants                  (  1  )
 plcmx_in    ,          & ! vegetation: fraction covered by plants      (  1  )
 plcmn_in    ,          & ! rest:       fraction covered by plants      (  1  )
 rlaimx_in   ,          & ! vegetation: leaf area index                 (  1  )
 rlaimn_in   ,          & ! rest:       leaf area index                 (  1  )
 vio3_in     ,          & ! total vertically integrated ozone content   (Pa O3)
 hmo3_in     ,          & ! height of maximum ozone concentration       ( Pa  )
 lat_coarse_m,          & ! latitudes of the LM grid points
 lon_coarse_m             ! longitudes of the LM grid points

!------------------------------------------------------------------------------

USE data_grid_lm,    ONLY: &
  ielm_tot,    & ! ie for LM, total domain
  jelm_tot,    & ! je for LM, total domain
  kelm,        & ! ke for LM
  ke1lm,       & ! ke+1
  dlat,        & ! grid point distance in zonal direction (in degrees)
  dlon,        & ! grid point distance in meridional direction (in degrees)
  pollat,      & ! latitude of the rotated north pole (in degrees, N>0)
  pollon,      & ! longitude of the rotated north pole (in degrees, E>0)
  polgam,      & ! latitude of the rotated north pole  !_br
  startlat_tot,& ! transformed latitude of the lower left grid point
                 ! of the total domain (in degrees, N>0)
  startlon_tot,& ! transformed longitude of the lower left grid point
                 ! of the total domain (in degrees, E>0)
  endlat_tot,  & ! transformed latitude of the upper right grid point
                 ! of the total domain (in degrees, N>0)
  endlon_tot,  & ! transformed longitude of the upper right grid point
                 ! of the total domain (in degrees, E>0)
  startlat,    & ! transformed latitude of the lower left grid point
                 ! of the local domain (in degrees, N>0)
  startlon,    & ! transformed longitude of the lower left grid point
                 ! of the local domain (in degrees, E>0)
  ie2lm_tot,   & ! = ielm_tot + 2
  je2lm_tot,   & ! = jelm_tot + 2
  ie2lm,       & !
  je2lm,       & !
  ie_ext,      & ! west-east size of fields with external parameters
  je_ext         ! north-south size of fields with external parameters

!------------------------------------------------------------------------------

USE data_grid_in,    ONLY: &
  ie_in_tot,   & ! ie for input grid, total domain
  je_in_tot,   & ! je for input grid, total domain
  ni_gme,      & ! resolution of GME
  j1_min,      & ! smallest index in j1-direction for LM subdomain
  j1_max,      & ! biggest index in j1-direction for LM subdomain
  j2_min,      & ! smallest index in j2-direction for LM subdomain
  j2_max,      & ! biggest index in j2-direction for LM subdomain
  jd_min,      & ! smallest index for diamonds for a LM subdomain
  jd_max,      & ! biggest index for diamonds for a LM subdomain
  nd,          & ! number of diamonds (nd = ide-ids+1 = 10)
  ni2,         & ! ni_gme=3**ni3*2**ni2 with ni3 0 or 1 and ni2 > 1
  ni3,         & !
  igg1s,       & ! start index of global array-dimension in x-direction
  igg1sm1,     & ! = igg1s - 1
  igg1sm2,     & ! = igg1s - 2
  igg1e,       & ! end index of global array-dimension in x-direction
  igg1ep1,     & ! = igg1e + 1
  igg1ep2,     & ! = igg1e + 2
  igg2s,       & ! start index of global array-dimension in y-direction
  igg2sm1,     & ! = igg2s - 1
  igg2sm2,     & ! = igg2s - 2
  igg2e,       & ! end index of global array-dimension in y-direction
  igg2ep1,     & ! = ig2e + 1
  igg2ep2,     & ! = igg2e + 2
  ispoke         ! offsets of the 6 (5) neighbouring gridpoints relative to
                 ! i1-direction use ispoke(m), m=1,6 (5), in i2-direction
                 ! use ispoke(m+6), m=1,6 (5); phys. dim. ( - )

USE data_grid_in,    ONLY: &
  ie_in,       & ! ie for input grid, local domain
  je_in,       & ! je for input grid, local domain
  ke_in,       & ! ke for input grid
  ke1in,       & !
  pollat_in,   & ! latitude of the rotated north pole (in degrees, N>0)
  pollon_in,   & ! longitude of the rotated north pole (in degrees, E>0)
  polgam_in,   & ! latitude of the rotated north pole
  startlat_in, & ! transformed latitude of the lower left grid point
                 ! of the local domain (in degrees, N>0)
  startlon_in, & ! transformed longitude of the lower left grid point
                 ! of the local domain (in degrees, E>0)
  startlat_in_tot,& ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
  startlon_in_tot,& ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
  latitudes_in,& ! latitudes of the input data
  longitudes_in,&! longitudes of the input data
  grdpt_rel_in,& ! relation between input and lm grid for cressman scheme
  lcm_hgt_coor   ! Input data has hybrid height coordinates

!------------------------------------------------------------------------------

USE data_int2lm_control,    ONLY: &
  noutput,      & ! unit number for output file
  lgme2lm,      & ! if .TRUE., gme->lm
  lgsm2lm,      & ! if .TRUE., gsm->lm
  lgfs2lm,      & ! if .TRUE., gfs->lm
  llm2lm,       & ! if .TRUE., lm ->lm
  lum2lm,       & ! if .TRUE., um ->lm
  lhir2lm,      & ! if .TRUE., hirlam ->lm
  lec2lm,       & ! if .TRUE., ec ->lm
  lcm2lm,       & ! if .TRUE., cm ->lm
  lforest,      & ! if .TRUE., run with forest (evergreen and deciduous)
  lurban,       & ! if .TRUE., run the urban module
  lsso,         & ! process parameters for sso scheme
  lradtopo,     & ! process parameters for topographic correction of radiation
  llake,        & ! if .TRUE., run with lake
  lemiss,       & ! if .TRUE., run with external parameter for surface emissivity
  lstomata,     & ! if .TRUE., run with external parameter for stomata resistance
  lfilter_oro,  & ! if .TRUE., filter the orography
  llbc_smooth,  & ! if .TRUE., run with smooth orography transition at LB
  nlbc_smooth,  & ! number of grip points for smooth orography transition at LB
  l_smi,        & ! if .TRUE., interpolate soil moisture with FC-PWP index
  l_cressman,   & ! switch for using a Cressman scheme during M-type interpolation
  l_bicub_spl,  & ! switch for using a bicubic spline interpolation
  itype_t_cl,   & ! to choose origin and treatment of deep soil temperature
  itype_rootdp, & ! to choose treatment of root depth
  itype_ndvi,   & ! to choose treatment of surface parameters (plcov, lai)
  itype_aerosol,& ! to choose treatment of surface parameters (plcov, lai)
  itype_albedo, & ! to choose treatment of solar surface albedo
  kcontrol_fi,  & ! control level for geopotential
  norder_filter,& ! order of the orography filtering
  eps_filter,   & ! parameter for orography filtering
  dt,           & ! time step used in the LM
  idbg_level,   & ! to control verbosity of output
  lprintdeb_all,& ! whether all PEs print debug output
  nl_soil_lm      ! number of soil layers in LM, resp. HM

!------------------------------------------------------------------------------

USE data_int2lm_constants, ONLY : &
    fcb_ec,     & ! field capacity of IFS soil types with CY32r3
    fcb_ec_1s,  & ! field capacity of IFS 1 soil type
    pwpb_ec,    & ! permanent wilting point of IFS soil types with CY32r3
    pwpb_ec_1s, & ! permanent wilting point of IFS 1 soil type
    Pi,         & ! circle constant
    r_d,        & ! gas constant for dry air          [J/K*kg]
    T0,         & ! 0 degree Celsius                              [Kelvin]
    G             ! gravity at sea level                          [ms-2]

!------------------------------------------------------------------------------

USE data_int2lm_io,        ONLY : &
  nuchkdat,          & ! checking the I/O data
  yuchkdat,          & ! checking the I/O data
  undefgrib,         & ! value for "undefined" in the grib routines
  undefncdf,         & ! value for "undefined" in the netcdf routines
  undef,             & ! the same with other KIND-Parameter
  ylmext_form_read,  & ! input format of external LM data
  yinext_form_read,  & ! input format of external boundary data
  ylm_form_write,    & ! output format of LM data
  ylm_hhl,           & ! name of the file LM HHL fields
  ylmext_cat,        & ! catalog-name of the file with external LM parameters
  lchkin               ! logical for print of check-values (max,min,mean)

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY :  &
    nprocx,          & ! number of processors in x-direction for LM
    nprocy,          & ! number of processors in y-direction for LM
    isubpos,         & ! positions of the subdomains in the total domain
    isubpos_coarse,  & ! the same for the coarse grid
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    num_compute,     & ! number of compute PEs
    icomm_cart,      & ! communicator that belongs to igroup_world, i.e.
    imp_integers,    & ! determines the correct INTEGER type used in the model
                       ! for MPI
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    my_cart_id         ! rank of this subdomain in the cartesian communicator

!------------------------------------------------------------------------------

USE data_profiles,             ONLY :  &
    niso_loc,              & ! number of isolated points in this subdomain
    niso_max,              & ! maximal number of isolated points in one subdomain
    niso_tot,              & ! total number of isolated points in the whole COSMO domain
    struct_for_isolated,   & ! data structure to keep isolated points
    globl_iso_points,      & ! global structure
    local_iso_points         ! local  structure

!------------------------------------------------------------------------------

USE environment,          ONLY :  get_free_unit
USE gme_utilities,        ONLY :  init_gme_interpol, pp_interp2ls
USE interp_utilities,     ONLY :  interp_l, interp_q_bs
USE parallel_utilities,   ONLY :  gather_field, distribute_field,    &
                                  i_global, j_global, global_values, gather_values
USE utilities,            ONLY :  sleve_split_oro,     &
                                  phi2phirot, rla2rlarot, rlarot2rla, phirot2phi

USE vgrid_refatm_utils,   ONLY :  reference_atmosphere, reference_atmosphere_2, &
                                  vcoord, refatm, lanalyt_calc_t0p0, vcoord_in, &
                                  lnewVGrid, nfltvc, svc1, svc2, vcoord_d,      &
                                  lhhl_lm_read

USE src_read_ext,         ONLY :  read_lm_ext, read_coarse_grid_ext
USE src_read_hhl,         ONLY :  org_read_hhl

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Parameterlist
CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

! Local arrays:
REAL (KIND=ireals)           ::   &
  zvegfac(ie2lm,je2lm),           & ! vegetation factor depending on grid 
                                    ! point and time of year
  zhhlr  (ke1lm), zhhlr_in(ke1in), zak(ke1in), zbk(ke1in)

! Local scalars:
INTEGER  (KIND=iintegers)       ::  &
  izerror, izerrstat, niostat,      & ! status and error status variables
  nzactday, nzday,                  & ! actual day of the year
  i, j, k, i_c, j_c,                & ! loop indices
  kflat,                            & !
  icn, jcn, n, n1, ncount, ii1, jj1,& !
  ilow, jlow, ihig, jhig, ix, jy,   & !
  istart, iend, jstart, jend,       & !
  istart_in, iend_in, jstart_in, jend_in,       & !
  izerrflake,                       & !
  izdebug, nzbounds                   ! for debug output

REAL (KIND=ireals)         ::  &
  zacthour,                         & ! actual hour of the day
  zhred, zbvp, zdvp,                & ! factors for vegetation period
  zwei,                             & ! weighting factor: 1=COARSE, 0=FINE
  fr_land_orig, soiltyp_orig,       & !
  zdist_min, zdist, lat_coarse, lon_coarse, lon_gpc, lat_gpc, & !
  distlat, distlon, lon_gpc2, lat_gpc2

REAL (KIND=ireals), ALLOCATABLE  ::  &
  zhsurfs_lm    (:,:,:), & ! height of splitted topography parts
  zhsurfs_lm_tot(:,:,:), & ! height of splitted topography parts of full domain
  zhsurf_lm_tot (:,:)  , & ! full topography of full domain
  zrho0         (:,:,:), & ! intermediate storage
  zdp0          (:,:,:), &
  zp0hl         (:,:,:), &
  zt0hl         (:,:,:), &
  zhhl_read     (:,:,:), & ! to read a GRIB2 HHL-file into a temporary variable
  zlsm_cosmo_tot(:,:)  , &
  zlsm_coars_tot(:,:)

LOGICAL                    ::  &
  ! logicals, to indicate whether the corresponding field has been read or
  ! interpolated (lm-fields from input-fields)
  lhsur_lm, lfis__lm, lfrla_lm, lz0___lm, lsoty_lm, lplmx_lm, lplmn_lm, &
  laimx_lm, laimn_lm, lroot_lm, lfore_lm, lford_lm, lemis_lm, lprsm_lm, &
  lurba_lm, laldr_lm, lalsa_lm,                                         &
            lfis__in, lfrla_in, lz0___in, lsoty_in, lplmx_in, lplmn_in, &
  laimx_in, laimn_in, lroot_in,                                         &
  lhsur_int2lm, lfis__int2lm, lfrla_int2lm, lsoty_int2lm, lz0___int2lm, &
  lplmx_int2lm, lplmn_int2lm, laimx_int2lm, laimn_int2lm, lroot_int2lm, &
  lflak_lm, ldept_lm, ltcl__lm, lpl12_lm, lai12_lm, lz012_lm, lsalt_lm, &
  lstdh_lm, lgamm_lm, lthet_lm, lsigm_lm,                               &
  lskyv_lm, lsang_lm, lsasp_lm, lhori_lm,                               &
  lndvi_lm, lsu12_lm, ldu12_lm, lor12_lm, lbc12_lm, lss12_lm, lal12_lm, &
  lfill_up

INTEGER(KIND=iintegers), ALLOCATABLE  :: izbufsend(:), izbufrecv(:,:)

LOGICAL                    ::  &
  lzopen,        & ! to test whether files are open
  lzframe,       & ! dummy for the interpolation routine (=.FALSE.)
  lz_lsm_missing,& ! if TRUE, a land-sea mask is missing for further calculations
  lzisolated,    & !
  lmono, lposdef   ! for interpolation types

CHARACTER (LEN=  1)        ::  &
  yitype     ! interpolation type

CHARACTER (LEN=80)         ::  &
  yzerror    ! error message for error handling

!
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierror    = 0_iintegers
  izerrstat = 0_iintegers
  izerror   = 0_iintegers
  yzerror   = '     '

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  niostat   = 0_iintegers
  IF     (yinext_form_read == 'grb1') THEN
    undef     = REAL (undefgrib, ireals)
  ELSEIF (yinext_form_read == 'ncdf') THEN
    undef     = REAL (undefncdf, ireals)
  ENDIF
  lzframe   = .FALSE.   ! no frames for the interpolation of external data

  ! get free unit number for file YUCHKDAT
  yuchkdat  = 'YUCHKDAT'
  CALL get_free_unit (nuchkdat)

  ! open file YUCHKDAT and print a headline
  IF ( (lchkin) .AND. (my_cart_id == 0) ) THEN
    OPEN(nuchkdat, FILE=yuchkdat, FORM='FORMATTED', STATUS='UNKNOWN',  &
                   POSITION='APPEND', IOSTAT=niostat)

    IF(niostat /= 0) THEN
      ierror = 1
      yerror = ' ERROR    *** Error while opening file YUCHKDAT *** '
      RETURN
    ENDIF
  ENDIF

  ! Set all logicals to .FALSE.
  lhsur_lm = .FALSE.;                              lhsur_int2lm = .FALSE.
  lfis__lm = .FALSE.;     lfis__in = .FALSE.;      lfis__int2lm = .FALSE.
  lfrla_lm = .FALSE.;     lfrla_in = .FALSE.;      lfrla_int2lm = .FALSE.
  lz0___lm = .FALSE.;     lz0___in = .FALSE.;      lz0___int2lm = .FALSE.
  lz012_lm = .FALSE.;
  lsoty_lm = .FALSE.;     lsoty_in = .FALSE.;      lsoty_int2lm = .FALSE.
  lplmx_lm = .FALSE.;     lplmx_in = .FALSE.;      lplmx_int2lm = .FALSE.
  lplmn_lm = .FALSE.;     lplmn_in = .FALSE.;      lplmn_int2lm = .FALSE.
  lpl12_lm = .FALSE.;
  laimx_lm = .FALSE.;     laimx_in = .FALSE.;      laimx_int2lm = .FALSE.
  laimn_lm = .FALSE.;     laimn_in = .FALSE.;      laimn_int2lm = .FALSE.
  lai12_lm = .FALSE.;
  lroot_lm = .FALSE.;     lroot_in = .FALSE.;      lroot_int2lm = .FALSE.

  lfore_lm = .FALSE.;     lford_lm = .FALSE.;      lemis_lm     = .FALSE.
  lurba_lm = .FALSE.;     lalsa_lm = .FALSE.;      laldr_lm     = .FALSE.
  lprsm_lm = .FALSE.;
  lstdh_lm = .FALSE.;     lgamm_lm = .FALSE.;      lthet_lm     = .FALSE.
  lsigm_lm = .FALSE.;
  lskyv_lm = .FALSE.;     lsang_lm = .FALSE.;      lsasp_lm     = .FALSE.
  lhori_lm = .FALSE.;
  lflak_lm = .FALSE.;     ldept_lm = .FALSE.;      lsalt_lm     = .FALSE.
  ltcl__lm = .FALSE.;

  lndvi_lm = .FALSE.;

  lsu12_lm = .FALSE.;     ldu12_lm = .FALSE.;      lor12_lm     = .FALSE.
  lbc12_lm = .FALSE.;     lss12_lm = .FALSE.;      lal12_lm     = .FALSE.

  ! Initialization of grdpt_rel_in to 0; 
  ! (is this really intended? how can that be controlled?)
  grdpt_rel_in = 0

!------------------------------------------------------------------------------
! Section 2: Read external parameters for output LM grid
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, ' Read external parameters for COSMO grid'
  ENDIF

  CALL read_lm_ext (lhsur_lm, lfis__lm, lfrla_lm, lz0___lm, lz012_lm,       &
                    lsoty_lm, lroot_lm, lplmx_lm, lplmn_lm, lpl12_lm,       &
                    laimx_lm, laimn_lm, lai12_lm, lfore_lm, lford_lm,       &
                    lurba_lm, lalsa_lm, laldr_lm,                           &
                    lemis_lm, lprsm_lm, lflak_lm, ldept_lm, ltcl__lm,       &
                    lndvi_lm, lstdh_lm, lgamm_lm, lthet_lm, lsigm_lm,       &
                    lskyv_lm, lsang_lm, lsasp_lm, lhori_lm, lsalt_lm,       &
                    lsu12_lm, ldu12_lm, lor12_lm, lbc12_lm, lss12_lm,       &
                    lal12_lm,                                               &
                    izerror, yzerror)

  IF (izerror /= 0_iintegers) THEN
    ierror = 2
    yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
    RETURN
  ENDIF

  ! Compute hsurf_lm from fis_lm or vice versa, if one of these fields
  ! has not been read. If both have not be read, they could be interpolated
  ! from the coarse fields later on.
  IF (.NOT. lfilter_oro) THEN
    IF (.NOT. lhsur_lm .AND. lfis__lm) THEN
      hsurf_lm(:,:) = fis_lm(:,:) / g
      lhsur_lm = .TRUE.
    ENDIF
    IF (.NOT. lfis__lm .AND. lhsur_lm) THEN
      fis_lm(:,:) = hsurf_lm(:,:) * g
      lfis__lm = .TRUE.
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Read external parameters for the coarse grid
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, ' Read external parameters for coarse grid'
  ENDIF

  CALL read_coarse_grid_ext                                                  &
            (lfis__in, lfrla_in, lsoty_in, lz0___in, lplmx_in, lplmn_in,     &
             laimx_in, laimn_in, lroot_in,                                   &
             lfis__lm, lfrla_lm, lsoty_lm, lz0___lm, lplmx_lm, lplmn_lm,     &
             laimx_lm, laimn_lm, lroot_lm, izerror, yzerror)

  IF (izerror /= 0_iintegers) THEN
    ierror = 3
    yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Do necessary calculations for INT2LM:
!            fr_land, geopotential of surface, soiltyp, surface roughness
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, ' Additional calculations for external parameters'
  ENDIF

  ! Set defaults for interpolation (which are valid for most external
  ! parameters)
  lmono          = .FALSE.
  lposdef        = .FALSE.
  lz_lsm_missing = .FALSE.

  !----------------------------------------------------------------------------
  ! Section 4.1: Fraction of Land and Land-Sea-Mask for LM and GME
  !----------------------------------------------------------------------------

  ! Compute lolp_in
  IF (lfrla_in) THEN
    IF (lgme2lm) THEN
      WHERE (fr_land_gme >= 0.5_ireals)
        lolp_gme = .TRUE.   ! land points
      ELSEWHERE
        lolp_gme = .FALSE.  ! sea points
      ENDWHERE
    ELSE !  all other models
      WHERE (fr_land_in >= 0.5_ireals)
        lolp_in  = .TRUE.   ! land points
      ELSEWHERE
        lolp_in  = .FALSE.  ! sea points
      ENDWHERE
    ENDIF
  ELSE
    ! fr_land_in is not available: set errorcode
    izerrstat = 1
    lz_lsm_missing = .TRUE.
  ENDIF

  ! Compute lolp_lm, if fr_land_lm is available,
  ! if not, interpolate fr_land_in to fr_land_lm first
  IF (lfrla_lm) THEN
    WHERE (fr_land_lm >= 0.5_ireals)
      lolp_lm = .TRUE.   ! land points
    ELSEWHERE
      lolp_lm = .FALSE.  ! sea points
    ENDWHERE
  ELSEIF (lfrla_in) THEN
    ! Try to interpolate fr_land_lm from coarse field
    ! linear interpolation: no land_sea_mask is required for that
    yitype  = 'L'

    IF (lgme2lm) THEN
!US!! ! This cannot be done with the (new) interpolation routine, because
      ! this definitely needs the weights and the land sea mask. 
      ! Implement the normal linear interpolation here directly!
      CALL pp_interp2ls(fr_land_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,   &
                        jd_min, jd_max,                                    &
                        ispoke,  baryll_m, index_m, lmono, .TRUE., yitype, &
                        lolp_gme, lolp_lm,  w_intpol, n_intpol, m_intpol,  &
                        l_intpol, undef,                                   &
                        fr_land_lm, 1, ie2lm, 1, je2lm, yzerror, izerror)
    ELSEIF (llm2lm .OR. lec2lm .OR. lcm2lm) THEN
      CALL interp_l                                                        &
              (fr_land_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),   &
               lmono, .TRUE., yitype, lzframe, lolp_in, lolp_lm,           &
               lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),              &
               fr_land_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,   &
               latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,    &
               grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,             &
               yzerror,  izerror)
    ENDIF

    IF (izerror == 0) THEN
      lfrla_int2lm = .TRUE.
    ELSE
      ierror = 4
      yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
      RETURN
    ENDIF

    WHERE (fr_land_lm >= 0.5_ireals)
      lolp_lm = .TRUE.   ! land points
    ELSEWHERE
      lolp_lm = .FALSE.  ! sea points
    ENDWHERE
  ELSE
    ! fr_land_gme/in and fr_land_lm are not available: set errorcode
    izerrstat = 1
    lz_lsm_missing = .TRUE.
  ENDIF

  IF (lgme2lm) THEN
    ! Initialize the interpolation from GME fields
    CALL init_gme_interpol                                                 &
             (lolp_gme, fr_land_gme,                                       &
              igg1sm2,  igg1ep2, igg2sm2, igg2ep2, jd_min, jd_max,         &
              lolp_lm, w_intpol, n_intpol, m_intpol, l_intpol,             &
              ispoke, baryll_m, index_m, 1, ie2lm, 1, je2lm,               &
              undef, izdebug, yzerror, izerror)
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 4.2: Geopotential of the surface
  !----------------------------------------------------------------------------

  IF (lfis__in) THEN
    IF (lfrla_in .AND. (lfrla_lm .OR. lfrla_int2lm)) THEN
      ! Interpolate fis from coarse grid to fis_gl with linear interpolation
      yitype  = 'L'
      IF (lgme2lm) THEN
        CALL pp_interp2ls(fis_gme, igg1sm2,  igg1ep2, igg2sm2, igg2ep2,      &
                          jd_min, jd_max,                                    &
                          ispoke,  baryll_m, index_m, lmono, lposdef, yitype,&
                          lolp_gme, lolp_lm,  w_intpol, n_intpol, m_intpol,  &
                          l_intpol, undef,                                   &
                          fis_gl, 1, ie2lm, 1, je2lm, yzerror, izerror)
      ELSEIF (lec2lm .OR. lcm2lm .OR. lgsm2lm .OR. lgfs2lm .OR. lhir2lm) THEN
        IF (l_bicub_spl) THEN
          CALL interp_q_bs                                                   &
             (fis_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),          &
              lmono, lposdef, lzframe, lolp_in, lolp_lm,                     &
              undef, lmask_lm, x_wght(:,:,1), y_wght(:,:,1),                 &
              fis_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,          &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,       &
              grdpt_rel_in, ie_in_tot, je_in_tot, lcm2lm, l_cressman,        &
              yzerror, izerror)
        ELSE
          CALL interp_l                                                      &
                  (fis_in,     ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
                   lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm,  &
                   lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),            &
                   fis_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,     &
                   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,  &
                   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,           &
                   yzerror,  izerror)
        ENDIF

        ! Also interpolate hsurf_in in this branch (lcm2lm)
        ! This is the same code as below for lum2lm
        IF (lcm_hgt_coor) THEN
          CALL interp_l                                                      &
                (hsurf_in,   ie_in, je_in, i_index(:,:,1), j_index(:,:,1),   &
                 lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm,    &
                 lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),              &
                 hsurf_gl(:,:), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                 latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,    &
                 grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,             &
                 yzerror,  izerror)
        ENDIF
      ELSEIF (llm2lm) THEN
        ! in case of lm2lm hsurf has to be interpolated
        CALL interp_l                                                        &
                (hsurf_in,   ie_in, je_in, i_index(:,:,1), j_index(:,:,1),   &
                 lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm,    &
                 lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),              &
                 hhl_gl(:,:,ke1in), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                 latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,    &
                 grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,             &
                 yzerror,  izerror)
        IF (vcoord_in%ivctype == 3 .OR. vcoord_in%ivctype == 4) THEN
          CALL interp_l                                                        &
                  (hsurfs_in(:,:,1),ie_in,je_in,i_index(:,:,1),j_index(:,:,1), &
                   lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm,    &
                   lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),              &
                   hsurfs_gl(:,:,1), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,    &
                   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,             &
                   yzerror,  izerror)
          CALL interp_l                                                        &
                  (hsurfs_in(:,:,2),ie_in,je_in,i_index(:,:,1),j_index(:,:,1), &
                   lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm,    &
                   lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),              &
                   hsurfs_gl(:,:,2), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,    &
                   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,             &
                   yzerror,  izerror)
        ENDIF
      ELSEIF (lum2lm) THEN
        ! in case of um2lm hsurf has to be interpolated to the field hsurf_gl
        CALL interp_l                                                        &
                (hsurf_in,   ie_in, je_in, i_index(:,:,1), j_index(:,:,1),   &
                 lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm,    &
                 lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),              &
                 hsurf_gl(:,:), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                 latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,    &
                 grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,             &
                 yzerror,  izerror)
      ENDIF

      IF ( (izerror == 0) .AND. (.NOT. lfis__lm) ) THEN
        ! fis_lm and hsurf_lm are not available: set fis_lm = fis_gl
        ! and compute hsurf_lm
        lfis__int2lm = .TRUE.
        lhsur_int2lm = .TRUE.
        IF (lgme2lm .OR. lec2lm .OR. lcm2lm) THEN
          fis_lm(:,:)   = fis_gl(:,:)
          hsurf_lm(:,:) = fis_lm(:,:) / g
        ELSEIF (llm2lm) THEN
          hsurf_lm(:,:) = hhl_gl(:,:,ke1in)
          fis_lm  (:,:) = hsurf_lm(:,:) * g
        ENDIF
      ELSEIF (izerror /= 0) THEN
        ierror = 5
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
    ENDIF

    IF (llbc_smooth) THEN
      ! At the lateral boundaries a smooth transition of nlbc_smooth grid points
      ! from the coarse to the fine orography is constructed with a linear
      ! weighting factor zwei. From 1 for coarse to 0 for fine orography.

      k = nlbc_smooth
      ! Add nboundlines grid points so that final boundary has factor 1.
      ! k - j_c + 1     (bottom, South) and for index i_c (left,  West) and
      ! j_c-je2lm_tot+k (top,    North) and for index i_c (right, East)
      DO j = 1, je2lm
        DO i = 1, ie2lm
          j_c=j_global(j)  ! global N-S position
          i_c=i_global(i)  ! global W-E position
          zwei=MAX(0.0_ireals,REAL(k - j_c + 1    +nboundlines,ireals)/REAL(k,ireals))
          zwei=MAX(      zwei,REAL(j_c-je2lm_tot+k+nboundlines,ireals)/REAL(k,ireals))
          zwei=MAX(      zwei,REAL(k - i_c + 1    +nboundlines,ireals)/REAL(k,ireals))
          zwei=MAX(      zwei,REAL(i_c-ie2lm_tot+k+nboundlines,ireals)/REAL(k,ireals))
          ! Now zwei>1.0 on boundaries and zwei=1.0 at nboundlines from boundaries
          ! Set ALL to max. 1.0 so that boundaries are now 1 + nboundlines broad
          IF ( zwei > 1 ) zwei = MIN (1.0_ireals, zwei)
          ! Use hhl_gl or fis_gl to define NEW fis_lm and hsurf_lm
          IF (llm2lm) THEN
            hsurf_lm(i,j)=zwei*hhl_gl(i,j,ke1in)+(1.0_ireals-zwei)*hsurf_lm(i,j)
            fis_lm  (i,j)=hsurf_lm(i,j) * g
          ELSE
            fis_lm  (i,j)=zwei*fis_gl(i,j)      +(1.0_ireals-zwei)*fis_lm(i,j)
            hsurf_lm(i,j)=fis_lm(i,j) / g
          ENDIF
        ENDDO
      ENDDO
    ENDIF

  ELSE
    ! fis_gme or land-sea-masks not available: set errorcode
    izerrstat = 1
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 4.3: Soiltyp
  !----------------------------------------------------------------------------

  IF (lsoty_in) THEN
    IF (.NOT. lsoty_lm) THEN
      ! Interpolate soiltyp from coarse grid to soiltyp_lm 
      ! (if no error occured above)
      IF (lfrla_in .AND. (lfrla_lm .OR. lfrla_int2lm)) THEN
        ! interpolation with 'N': nearest grid point
        yitype  = 'N'
        IF (lgme2lm) THEN
          CALL pp_interp2ls                                                 &
                 (soiltyp_gme, igg1sm2,  igg1ep2, igg2sm2, igg2ep2,         &
                  jd_min, jd_max,                                           &
                  ispoke,  baryll_m, index_m, lmono, lposdef, yitype,       &
                  lolp_gme, lolp_lm,  w_intpol, n_intpol, m_intpol,         &
                  l_intpol, undef,                                          &
                  soiltyp_lm, 1, ie2lm, 1, je2lm, yzerror, izerror)
        ELSEIF (llm2lm) THEN
          CALL interp_l                                                     &
                  (soiltyp_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),&
                   lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm, &
                   lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),           &
                   soiltyp_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
                   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,          &
                   yzerror,  izerror)
        ELSEIF (lec2lm) THEN
          ierror = 6
          yerror = 'IFS and LM soiltypes are incompatible: cannot interpolate'
          RETURN
        ENDIF
        IF (izerror == 0) THEN
          lsoty_int2lm = .TRUE.
        ELSE
          ierror = 6
          yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
          RETURN
        ENDIF
      ENDIF
    ENDIF
  ELSE ! gdm < lsoty_in = .F.
    IF (lec2lm) THEN
      IF (l_smi) THEN
        ! gdm l_smi CLEPS begin: modif. for reforecasts with ERA40 (IFS doc. CY23r4)
        soiltyp_in(:,:) = 4 ! "sandy-loam"
        fcb_ec(4) = fcb_ec_1s
        pwpb_ec(4) = pwpb_ec_1s
        WRITE (*,*) 'BEWARE: forced soiltyp_in to sandy-loam with *_ec_1s'
        ! gdm l_smi CLEPS end
      ELSE ! Force usage of single soil type (*_ec_1s)
        ! removed the WHERE lolp_in statement
        ! the whole field must be initialized with the same value, because in
        ! src_coarse_interpol it is checked, whether MINVAL /= MAXVAL!
        soiltyp_in(:,:) = 1    ! ECMWF and LM soiltypes are incompatible
      ENDIF
    ELSEIF (.NOT. lsoty_lm) THEN
      ! soiltyp from coarse model is not available: set errorcode
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 4.4: Surface roughness
  !----------------------------------------------------------------------------

  IF (.NOT. lz0___lm .AND. .NOT. lz012_lm) THEN
    IF (lz0___in) THEN
      IF (lfrla_in .AND. (lfrla_lm .OR. lfrla_int2lm)) THEN
        ! interpolate with the defaults
        yitype  = 'M'
        IF (lgme2lm) THEN
          CALL pp_interp2ls                                                 &
                 (z0_gme     , igg1sm2,  igg1ep2, igg2sm2, igg2ep2,         &
                  jd_min, jd_max,                                           &
                  ispoke,  baryll_m, index_m, lmono, lposdef, yitype,       &
                  lolp_gme, lolp_lm,  w_intpol, n_intpol, m_intpol,         &
                  l_intpol, undef,                                          &
                  z0_lm,  1, ie2lm, 1, je2lm, yzerror, izerror)
        ELSEIF (llm2lm .OR. lec2lm .OR. lcm2lm) THEN
         IF (yitype=='M' .AND. l_cressman) grdpt_rel_in=1
          CALL interp_l                                                     &
                  (z0_in,      ie_in, je_in, i_index(:,:,1), j_index(:,:,1),&
                   lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm, &
                   lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),           &
                   z0_lm,  1, ie2lm, 1, je2lm, startlat_in, startlon_in,    &
                   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
                   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,          &
                   yzerror,  izerror)
        ENDIF
        IF (izerror == 0) THEN
          lz0___int2lm = .TRUE.
        ELSE
          ierror = 7
          yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
          RETURN
        ENDIF
      ! else land-sea-masks not available
      ENDIF
    ELSE
      ! surface roughness from coarse grid model not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 4.5: Adapt fraction of land and land-sea-mask, if wanted
  !----------------------------------------------------------------------------

  ! Especially in small high resolution areas it could be that there are lakes
  ! ("water points") in the small domain without an appropriate point in the
  ! coarse domain. The interpolation would fail in these cases.
  ! As a first "quick and dirty" solution, you can fill up these points

  lfill_up  = .FALSE.    ! this can be changed on your own risk, if you have
                         ! problems with small lakes in a high resolution area

  IF (lfill_up .AND. .NOT. lz_lsm_missing) THEN
    DO j = 1, je2lm
      DO i = 1, ie2lm
        i_c = i_index(i,j,1) ! first input grid dimension
        j_c = j_index(i,j,1) ! second input grid dimension
        IF (.NOT. lolp_lm(i,j)) THEN
          IF ( (lolp_lm(i,j) .NEQV. lolp_in(i_c  ,j_c  ))  .AND. &
               (lolp_lm(i,j) .NEQV. lolp_in(i_c+1,j_c  ))  .AND. &
               (lolp_lm(i,j) .NEQV. lolp_in(i_c  ,j_c+1))  .AND. &
               (lolp_lm(i,j) .NEQV. lolp_in(i_c+1,j_c+1)) )        THEN
            fr_land_orig    = fr_land_lm(i,j)
            soiltyp_orig    = soiltyp_lm(i,j)
            lolp_lm   (i,j) = .TRUE.     ! make it a land point
            fr_land_lm(i,j) = 0.6        ! just something > 0.5
            IF (soiltyp_lm(i,j) < 3 .OR. soiltyp_lm(i,j) > 8) THEN
              soiltyp_lm(i,j) = soiltyp_in(i_c,j_c)
            ENDIF
            WRITE (*,'(A,I3,A,I3,A)') '    Filled up grid point (',        &
                       i_global(i),',', j_global(j), ') with land:' 
            WRITE (*,'(A,F6.2,A,F6.2)') '        fr_land:   orig = ',      &
                       fr_land_orig, '    now = ', fr_land_lm(i,j)
            WRITE (*,'(A,F6.2,A,F6.2)') '        soil type: orig = ',      &
                       soiltyp_orig, '    new = ', soiltyp_lm(i,j)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 4.6: Test whether isolated points are present
  !----------------------------------------------------------------------------

  IF (lz_lsm_missing) THEN
    PRINT *, ' *** WARNING: isolated points cannot be calculated due to missing land-sea-mask'
  ENDIF

  IF (.NOT. lgme2lm .AND. .NOT. lz_lsm_missing) THEN
    IF (izdebug > 5) THEN
      PRINT *, ' Check for isolated points'
    ENDIF

    nzbounds = 2 ! has to correspond to nzbounds in src_decomposition

    ! Allocate data structure for local isolated points
    ALLOCATE (local_iso_points(ie2lm*je2lm), STAT=izerror)

    ! Initialize this structure
    local_iso_points(:)%ic_glob  = -1
    local_iso_points(:)%jc_glob  = -1
    local_iso_points(:)%nc_task  = -1
    local_iso_points(:)%ic_locl  = -1
    local_iso_points(:)%jc_locl  = -1
    local_iso_points(:)%lland    = .FALSE.
    local_iso_points(:)%iin_glob = -1
    local_iso_points(:)%jin_glob = -1
    local_iso_points(:)%nin_task = -1
    local_iso_points(:)%iin_locl = -1
    local_iso_points(:)%jin_locl = -1
    local_iso_points(:)%value    = 0.0_ireals

    ! Look for isolated grid points in the local subdomain
    ! (according to the procedure in interp_l)
    niso_loc = 0 ! number of isolated points in this subdomain
    niso_max = 0 ! maximal number of isolated points in one subdomain
    niso_tot = 0 ! number of isolated points

    DO j = 1, je2lm
      DO i = 1, ie2lm
        ! Find out which is the nearest input grid point to the (i,j) COSMO grid point
        IF (x_wght(i,j,1) < 0.5_ireals) THEN
          icn = i_index(i,j,1)
        ELSE
          icn = i_index(i,j,1) + 1
        ENDIF
        IF (y_wght(i,j,1) < 0.5_ireals) THEN
          jcn = j_index(i,j,1)
        ELSE
          jcn = j_index(i,j,1) + 1
        ENDIF

        ! Check, whether any of the 25 surrounding points has the same 
        ! land sea mask
        lzisolated = .TRUE.
        DO jj1 = jcn-2, jcn+2
          DO ii1 = icn-2, icn+2
            IF (lolp_in(ii1,jj1) .EQV. lolp_lm(i,j)) THEN
              lzisolated = .FALSE.
            ENDIF
          ENDDO
        ENDDO
        
        IF (lzisolated) THEN
          niso_loc = niso_loc + 1
          local_iso_points(niso_loc)%ic_glob = i_global(i)
          local_iso_points(niso_loc)%jc_glob = j_global(j)
          local_iso_points(niso_loc)%nc_task = my_cart_id
          local_iso_points(niso_loc)%ic_locl = i
          local_iso_points(niso_loc)%jc_locl = j
          local_iso_points(niso_loc)%lland   = (fr_land_lm(i,j) >= 0.5_ireals)
        ENDIF
      ENDDO
    ENDDO

    ! Now there are niso_loc isolated points in this subdomain
    ! gather the maximum of all niso_loc to all subdomains, to allocate
    ! data structures for further gathering all these informations
    niso_max = niso_loc
    IF (num_compute > 1) THEN
      CALL global_values (niso_max, 1, 'MAX', imp_integers, icomm_cart,  &
                             -1, yzerror, izerror)

      ALLOCATE (izbufsend(6*niso_max+1),                 STAT=izerror)
      ALLOCATE (izbufrecv(6*niso_max+1,0:num_compute-1), STAT=izerror)

      izbufsend(1) = niso_loc
      DO n = 1, niso_loc
        izbufsend(1 + 6*(n-1) + 1) = local_iso_points(n)%ic_glob
        izbufsend(1 + 6*(n-1) + 2) = local_iso_points(n)%jc_glob
        izbufsend(1 + 6*(n-1) + 3) = local_iso_points(n)%nc_task
        izbufsend(1 + 6*(n-1) + 4) = local_iso_points(n)%ic_locl
        izbufsend(1 + 6*(n-1) + 5) = local_iso_points(n)%jc_locl
        IF (local_iso_points(n)%lland) THEN
          izbufsend(1 + 6*(n-1) + 6) = 1
        ELSE
          izbufsend(1 + 6*(n-1) + 6) = 0
        ENDIF
      ENDDO

      CALL gather_values(izbufsend, izbufrecv, 6*niso_max+1, num_compute, imp_integers,  &
                          -1, icomm_cart, yzerror, izerror)

      ! Compute total number of isolated points
      niso_tot = 0
      DO n = 0, num_compute-1
        niso_tot = niso_tot + izbufrecv(1,n)
      ENDDO

      ! Allocate structure to keep all isolated points
      ALLOCATE (globl_iso_points(niso_tot), STAT=izerror)

      ! Initialize this structure
      ! if there is no coarse grid point with matching land-sea-mask, 
      ! the iin_glob, jin_glob, nin_task indices will remain -1
      globl_iso_points(:)%ic_glob  = -1
      globl_iso_points(:)%jc_glob  = -1
      globl_iso_points(:)%nc_task  = -1
      globl_iso_points(:)%ic_locl  = -1
      globl_iso_points(:)%jc_locl  = -1
      globl_iso_points(:)%lland    = .FALSE.
      globl_iso_points(:)%iin_glob = -1
      globl_iso_points(:)%jin_glob = -1
      globl_iso_points(:)%nin_task = -1
      globl_iso_points(:)%iin_locl = -1
      globl_iso_points(:)%jin_locl = -1
      globl_iso_points(:)%value    = 0.0_ireals

      ! Sort all the gathered information to globl_iso_points
      ncount = 0
      DO n = 0, num_compute-1
        n1 = izbufrecv(1,n)
        IF (n1 > 0) THEN
          DO k = 1, n1
            ncount = ncount + 1
            globl_iso_points(ncount)%ic_glob = izbufrecv(1 + 6*(k-1) + 1, n)
            globl_iso_points(ncount)%jc_glob = izbufrecv(1 + 6*(k-1) + 2, n)
            globl_iso_points(ncount)%nc_task = izbufrecv(1 + 6*(k-1) + 3, n)
            globl_iso_points(ncount)%ic_locl = izbufrecv(1 + 6*(k-1) + 4, n)
            globl_iso_points(ncount)%jc_locl = izbufrecv(1 + 6*(k-1) + 5, n)
            IF (izbufrecv(1 + 6*(k-1) + 6, n) == 1) THEN
              globl_iso_points(ncount)%lland = .TRUE.
            ELSE
              globl_iso_points(ncount)%lland = .FALSE.
            ENDIF
            globl_iso_points(ncount)%value   = 0.0_ireals
          ENDDO
        ENDIF
      ENDDO

    ELSE
      niso_tot = niso_loc
      
      ! Allocate structure to keep all isolated points
      ALLOCATE (globl_iso_points(niso_tot), STAT=izerror)

      ! Initialize this structure
      globl_iso_points(:)%ic_glob  = -1
      globl_iso_points(:)%jc_glob  = -1
      globl_iso_points(:)%nc_task  = -1
      globl_iso_points(:)%ic_locl  = -1
      globl_iso_points(:)%jc_locl  = -1
      globl_iso_points(:)%lland    = .FALSE.
      globl_iso_points(:)%iin_glob = -1
      globl_iso_points(:)%jin_glob = -1
      globl_iso_points(:)%nin_task = -1
      globl_iso_points(:)%iin_locl = -1
      globl_iso_points(:)%jin_locl = -1
      globl_iso_points(:)%value    = 0.0_ireals

      DO n = 1, niso_tot
        globl_iso_points(n)%ic_glob = local_iso_points(n)%ic_glob
        globl_iso_points(n)%jc_glob = local_iso_points(n)%jc_glob
        globl_iso_points(n)%nc_task = local_iso_points(n)%nc_task
        globl_iso_points(n)%ic_locl = local_iso_points(n)%ic_locl
        globl_iso_points(n)%jc_locl = local_iso_points(n)%jc_locl
        globl_iso_points(n)%lland   = local_iso_points(n)%lland
        globl_iso_points(n)%value   = 0.0_ireals
      ENDDO
    ENDIF   

    ! Now determine for every isolated point the nearest point in the
    ! coarse grid with same land sea mask

    ! Not the full coarse grid domain can be used, but only the part
    ! covering the COSMO model domain. This part is specified in
    ! isubpos_coarse

    ! take maximum of all istart (in isubpos_coarse (j,1)
    istart_in = isubpos_coarse(0,1)
    DO j = 1, nprocy-1
      istart_in = MAX (istart_in, isubpos_coarse(j,1))
    ENDDO

    ! take minimum of all iend   (in isubpos_coarse (j,3)
    iend_in = isubpos_coarse((nprocx-1)*nprocy,3)
    DO j = 1, nprocy-1
      iend_in = MIN (iend_in, isubpos_coarse((nprocx-1)*nprocy+j,3))
    ENDDO

    ! take maximum of all jstart (in isubpos_coarse (i,2)
    jstart_in = isubpos_coarse(0,2)
    DO i = 1, nprocx-1
      jstart_in = MAX (jstart_in, isubpos_coarse(i*nprocy,2))
    ENDDO

    ! take minimum of all jend   (in isubpos_coarse (i,4)
    jend_in = isubpos_coarse(nprocy-1,4)
    DO i = 1, nprocx-1
      jend_in = MIN (jend_in, isubpos_coarse((i+1)*nprocy-1,4))
    ENDDO

    DO n = 1, niso_tot
      zdist_min = 1.0E20_ireals
      DO j = jstart_in, jend_in
        DO i = istart_in, iend_in
          ! lat and lon of COSMO grid point in coarse grid 
          IF ((pollat_in == pollat) .AND. (pollon_in == pollon)     &
                                    .AND. (polgam_in == polgam)) THEN
            lon_coarse = startlon_tot + (globl_iso_points(n)%ic_glob - 2)*dlon
            lat_coarse = startlat_tot + (globl_iso_points(n)%jc_glob - 2)*dlat
          ELSE
            IF ( ((pollat_in == 90.0_ireals) .AND. (pollon_in == 180.0_ireals)         &
                                             .AND. (polgam_in ==   0.0_ireals)) .OR.   &
                 ((pollat_in ==  0.0_ireals) .AND. (pollon_in ==   0.0_ireals)         &
                                             .AND. (polgam_in ==   0.0_ireals)) ) THEN
              ! If the coarse grid is unrotated, just transform to 
              ! geographical coordinates
              lon_gpc    = startlon_tot + (globl_iso_points(n)%ic_glob - 2)*dlon
              lat_gpc    = startlat_tot + (globl_iso_points(n)%jc_glob - 2)*dlat
              lon_coarse = rlarot2rla(lat_gpc, lon_gpc, pollat, pollon, polgam)
              lat_coarse = phirot2phi(lat_gpc, lon_gpc, pollat, pollon, polgam)
            ELSE
              ! rotate to geographical coordinates and then to the rotation of
              ! the input grid
              lon_gpc    = startlon_tot + (globl_iso_points(n)%ic_glob - 2)*dlon
              lat_gpc    = startlat_tot + (globl_iso_points(n)%jc_glob - 2)*dlat
              lon_gpc2   = rlarot2rla(lat_gpc, lon_gpc, pollat, pollon, polgam)
              lat_gpc2   = phirot2phi(lat_gpc, lon_gpc, pollat, pollon, polgam)
              lon_coarse = rla2rlarot(lat_gpc2, lon_gpc2, pollat_in, pollon_in, polgam_in)
              lat_coarse = phi2phirot(lat_gpc2, lon_gpc2, pollat_in, pollon_in)
            ENDIF
          ENDIF

          ! distance to COSMO-points
!          distlon = lon_coarse - (startlon_in_tot + (i-1)*dlon_in)
          ! should give the same results
          distlon = lon_coarse - (longitudes_in(i))
          IF (distlon <= -180.0_ireals) THEN
            ! this is most probably around the date line, because no regional
            ! domain is such big
            distlon = 360.0_ireals - ABS(distlon)
          ENDIF
!          distlat = lat_coarse - (startlat_in_tot + (j-1)*dlat_in)
          distlat = lat_coarse - (latitudes_in(j))   !_br 11.03.11
          ! should give the same results
          zdist   = SQRT(distlat*distlat + distlon*distlon)

          IF (  (zdist < zdist_min) .AND.                                       &
              ( (fland_in_tot(i,j) >= 0.5_ireals) .EQV. globl_iso_points(n)%lland) ) THEN
            zdist_min = zdist
            globl_iso_points(n)%iin_glob = i
            globl_iso_points(n)%jin_glob = j
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! Now get the subdomain and the local indices of that point
    ! for the coarse grid: for that search through all the subdomains
    DO n = 1, niso_tot
      procloop: DO n1 = 0, num_compute-1
        ilow = isubpos_coarse(n1,1)
        ihig = isubpos_coarse(n1,3)
        jlow = isubpos_coarse(n1,2)
        jhig = isubpos_coarse(n1,4)
        IF ( (ilow <= globl_iso_points(n)%iin_glob) .AND.           &
                     (globl_iso_points(n)%iin_glob <= ihig)   .AND. &
             (jlow <= globl_iso_points(n)%jin_glob) .AND.           &
                     (globl_iso_points(n)%jin_glob <= jhig) ) THEN

          globl_iso_points(n)%nin_task = n1
          globl_iso_points(n)%iin_locl = 1 + (globl_iso_points(n)%iin_glob - ilow)
          globl_iso_points(n)%jin_locl = 1 + (globl_iso_points(n)%jin_glob - jlow)

          EXIT procloop
        ENDIF
      ENDDO procloop

      IF (izdebug > 9) THEN
        WRITE (*,'(I5,7I6)') n, globl_iso_points(n)%iin_glob, globl_iso_points(n)%jin_glob,    &
                          globl_iso_points(n)%nin_task,                                  &
                          globl_iso_points(n)%iin_locl, globl_iso_points(n)%jin_locl, ie_in, je_in
      ENDIF
    ENDDO

    ! Copy the informations for the corresponding coarse grid points to the
    ! local structure local_iso_points:
    DO n = 1, niso_loc
      DO n1 = 1, niso_tot
        IF ( (my_cart_id                  ==  globl_iso_points(n1)%nc_task) .AND. &
             (local_iso_points(n)%ic_locl ==  globl_iso_points(n1)%ic_locl) .AND. &
             (local_iso_points(n)%jc_locl ==  globl_iso_points(n1)%jc_locl)) THEN 
          local_iso_points(n)%iin_glob  = globl_iso_points(n1)%iin_glob
          local_iso_points(n)%jin_glob  = globl_iso_points(n1)%jin_glob
          local_iso_points(n)%nin_task  = globl_iso_points(n1)%nin_task
          local_iso_points(n)%iin_locl  = globl_iso_points(n1)%iin_locl
          local_iso_points(n)%jin_locl  = globl_iso_points(n1)%jin_locl
          local_iso_points(n)%value     = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO

    IF (my_cart_id == 0) THEN
      WRITE (noutput, '(A)')  '  '
      WRITE (noutput, '(A,I6)')  '    There are isolated points in the COSMO Domain: ', niso_tot
      WRITE (noutput, '(A)')  '  '

      ! Some print outs:
      DO n = 1, niso_tot
        WRITE (noutput,'(6I5,L5,5I5,F5.2)')  n, globl_iso_points(n)
      ENDDO
    ENDIF

    ! Check, if there are isolated islands and only water points in the 
    ! coarse grid model: then we have to stop, because we cannot initialize 
    ! such fields
    IF ( (ALL(fland_in_tot(istart_in:iend_in,jstart_in:jend_in) < 0.5_ireals))  .AND. &
          ANY(globl_iso_points(:)%lland) ) THEN
      ierror = 12
      yerror = 'There are isolated islands, which cannot be initialized'
      RETURN
    ENDIF

    IF ( (ALL(fland_in_tot(istart_in:iend_in,jstart_in:jend_in) >= 0.5_ireals))  .AND. &
          ANY(.NOT. (globl_iso_points(:)%lland)) ) THEN
      ! Still possible: isolated water points (mainly inland water lakes), but
      !                 no water point in the coarse grid model:
      !         Then we just take the nearest grid point for initialization
      PRINT *, ' *** WARNING: There are isolated water points, that cannot be initialized  ***'
      PRINT *, ' ***          by coarse grid water points: nearest neighbour will be taken ***'
    ENDIF

  ENDIF ! not lgme2lm

!------------------------------------------------------------------------------
! Section 5: If no external parameters for root depth, plant cover and
!            leaf area index have been read for LM, try to interpolate
!            from coarse grid model.
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 5.1: Root depth
  !----------------------------------------------------------------------------

  IF (.NOT. lroot_lm) THEN
    IF (izdebug > 5) THEN
      PRINT *, ' Try to get root depth from coarse grid'
    ENDIF

    IF (lroot_in) THEN
      IF (lfrla_in .AND. (lfrla_lm .OR. lfrla_int2lm)) THEN
        ! interpolate with the defaults
        yitype  = 'M'
        IF (lgme2lm) THEN
          CALL pp_interp2ls                                                 &
                 (root_gme  ,  igg1sm2,  igg1ep2, igg2sm2, igg2ep2,         &
                  jd_min, jd_max,                                           &
                  ispoke,  baryll_m, index_m, lmono,   lposdef, yitype,     &
                  lolp_gme, lolp_lm,  w_intpol, n_intpol, m_intpol,         &
                  l_intpol, undef,                                          &
                  rootdp_lm, 1, ie2lm, 1, je2lm, yzerror, izerror)
        ELSEIF (llm2lm .OR. lec2lm .OR. lcm2lm) THEN
          IF (yitype=='M' .AND. l_cressman) grdpt_rel_in=1
          CALL interp_l                                                     &
                  (root_in,    ie_in, je_in, i_index(:,:,1), j_index(:,:,1),&
                   lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm, &
                   lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),           &
                   rootdp_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in, &
                   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
                   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,          &
                   yzerror,  izerror)
        ENDIF
        IF (izerror == 0) THEN
          lroot_int2lm = .TRUE.
        ELSE
          ierror = 8
          yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
          RETURN
        ENDIF
      ! else land-sea-masks not available
      ENDIF
    ELSE
      ! root depth from coarse grid not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.2: Plant cover
  !----------------------------------------------------------------------------

  ! Maximum
  IF (.NOT. lplmx_lm .AND. .NOT. lpl12_lm) THEN
    IF (izdebug > 5) THEN
      PRINT *, ' Try to get plant cover from coarse grid'
    ENDIF

    IF (lplmx_in) THEN
      IF (lfrla_in .AND. (lfrla_lm .OR. lfrla_int2lm)) THEN
        ! interpolate with the defaults
        yitype  = 'M'
        IF (lgme2lm) THEN
          CALL pp_interp2ls                                                   &
                 (plcmx_gme  , igg1sm2,  igg1ep2, igg2sm2, igg2ep2,           &
                  jd_min, jd_max,                                             &
                  ispoke,  baryll_m, index_m, lmono, lposdef, yitype,         &
                  lolp_gme, lolp_lm,  w_intpol, n_intpol, m_intpol,           &
                  l_intpol, undef,                                            &
                  plcov_mx_lm, 1, ie2lm, 1, je2lm, yzerror, izerror)
        ELSEIF (llm2lm .OR. lec2lm .OR. lcm2lm) THEN
          IF (yitype=='M' .AND. l_cressman) grdpt_rel_in=1
          CALL interp_l                                                       &
                  (plcmx_in,   ie_in, je_in, i_index(:,:,1), j_index(:,:,1),  &
                   lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm,   &
                   lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),             &
                   plcov_mx_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in, &
                   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                   yzerror,  izerror)
        ENDIF
        IF (izerror == 0) THEN
          lplmx_int2lm = .TRUE.
        ELSE
          ierror = 9
          yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
          RETURN
        ENDIF
      ! else land-sea-masks not available
      ENDIF
    ELSE
      ! plant cover (MX) from coarse grid not available
      izerrstat = 1
    ENDIF
  ENDIF

  ! Minimum
  IF (.NOT. lplmn_lm .AND. .NOT. lpl12_lm) THEN
    IF (lplmn_in) THEN
      IF (lfrla_in .AND. (lfrla_lm .OR. lfrla_int2lm)) THEN
        ! interpolate with the defaults
        yitype  = 'M'
        IF (lgme2lm) THEN
          CALL pp_interp2ls                                                   &
                 (plcmn_gme  , igg1sm2,  igg1ep2, igg2sm2, igg2ep2,           &
                  jd_min, jd_max,                                             &
                  ispoke,  baryll_m, index_m, lmono, lposdef, yitype,         &
                  lolp_gme, lolp_lm,  w_intpol, n_intpol, m_intpol,           &
                  l_intpol, undef,                                            &
                  plcov_mn_lm, 1, ie2lm, 1, je2lm, yzerror, izerror)
        ELSEIF (llm2lm .OR. lec2lm .OR. lcm2lm) THEN
          IF (yitype=='M' .AND. l_cressman) grdpt_rel_in=1
          CALL interp_l                                                       &
                  (plcmn_in,   ie_in, je_in, i_index(:,:,1), j_index(:,:,1),  &
                   lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm,   &
                   lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),             &
                   plcov_mn_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in, &
                   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                   yzerror,  izerror)
        ENDIF
        IF (izerror == 0) THEN
          lplmn_int2lm = .TRUE.
        ELSE
          ierror = 9
          yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
          RETURN
        ENDIF
      ! else land-sea-masks not available
      ENDIF
    ELSE
      ! plant cover (MN) from coarse grid not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.3: Leaf area index
  !----------------------------------------------------------------------------

  ! Maximum
  IF (.NOT. laimx_lm .AND. .NOT. lai12_lm) THEN
    IF (izdebug > 5) THEN
      PRINT *, ' Try to get leaf aread index from coarse grid'
    ENDIF

    IF (laimx_in) THEN
      IF (lfrla_in .AND. (lfrla_lm .OR. lfrla_int2lm)) THEN
        ! interpolate with the defaults
        yitype  = 'M'
        IF (lgme2lm) THEN
          CALL pp_interp2ls                                                   &
                 (rlaimx_gme  , igg1sm2,  igg1ep2, igg2sm2, igg2ep2,          &
                  jd_min, jd_max,                                             &
                  ispoke,  baryll_m, index_m, lmono, lposdef, yitype,         &
                  lolp_gme, lolp_lm,  w_intpol, n_intpol, m_intpol,           &
                  l_intpol, undef,                                            &
                  lai_mx_lm, 1, ie2lm, 1, je2lm, yzerror, izerror)
        ELSEIF (llm2lm .OR. lec2lm .OR. lcm2lm) THEN
          IF (yitype=='M' .AND. l_cressman) grdpt_rel_in=1
          CALL interp_l                                                       &
                  (rlaimx_in,   ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
                   lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm,   &
                   lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),             &
                   lai_mx_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,   &
                   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                   yzerror,  izerror)
        ENDIF
        IF (izerror == 0) THEN
          laimx_int2lm = .TRUE.
        ELSE
          ierror = 9
          yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
          RETURN
        ENDIF
      ! else land-sea-masks not available
      ENDIF
    ELSE
      ! leaf area index (MX) from coarse grid not available
      izerrstat = 1
    ENDIF
  ENDIF

  ! Minimum
  IF (.NOT. laimn_lm .AND. .NOT. lai12_lm) THEN
    IF (laimn_in) THEN
      IF (lfrla_in .AND. (lfrla_lm .OR. lfrla_int2lm)) THEN
        ! interpolate with the defaults
        yitype  = 'M'
        IF (lgme2lm) THEN
          CALL pp_interp2ls                                                   &
                 (rlaimn_gme  , igg1sm2,  igg1ep2, igg2sm2, igg2ep2,          &
                  jd_min, jd_max,                                             &
                  ispoke,  baryll_m, index_m, lmono, lposdef, yitype,         &
                  lolp_gme, lolp_lm,  w_intpol, n_intpol, m_intpol,           &
                  l_intpol, undef,                                            &
                  lai_mn_lm, 1, ie2lm, 1, je2lm, yzerror, izerror)
        ELSEIF (llm2lm .OR. lec2lm .OR. lcm2lm) THEN
          IF (yitype=='M' .AND. l_cressman) grdpt_rel_in=1
          CALL interp_l                                                       &
                  (rlaimn_in,   ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
                   lmono, lposdef, yitype, lzframe      , lolp_in, lolp_lm,   &
                   lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),             &
                   lai_mn_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,   &
                   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                   yzerror,  izerror)
        ENDIF
        IF (izerror == 0) THEN
          laimn_int2lm = .TRUE.
        ELSE
          ierror = 9
          yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
          RETURN
        ENDIF
      ! else land-sea-masks not available
      ENDIF
    ELSE
      ! leaf area index (MN) from coarse grid not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.4: Evergreen and deciduous forest
  !----------------------------------------------------------------------------

  IF (lforest) THEN
    IF (.NOT. lfore_lm) THEN
      ! fraction of evergreen forest not available
      izerrstat = 1
    ENDIF
    IF (.NOT. lford_lm) THEN
      ! fraction of deciduous forest not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.4.1: Urban land fraction
  !----------------------------------------------------------------------------

  IF (lurban) THEN
    IF (.NOT. lurba_lm) THEN
      ! fraction of urban land not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.5: Subgrid scale orography
  !----------------------------------------------------------------------------

  IF (lsso) THEN
    IF (.NOT. lstdh_lm) THEN
      ! sso_stdh not available
      izerrstat = 1
    ENDIF
    IF (.NOT. lgamm_lm) THEN
      ! sso_gamma
      izerrstat = 1
    ENDIF
    IF (.NOT. lthet_lm) THEN
      ! sso_theta
      izerrstat = 1
    ENDIF
    IF (.NOT. lsigm_lm) THEN
      ! sso_sigma
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.6: Topographical corrections
  !----------------------------------------------------------------------------

  IF (lradtopo) THEN
    IF (.NOT. lskyv_lm) THEN
      ! skyview not available
      izerrstat = 1
    ENDIF
    IF (.NOT. lsang_lm) THEN
      ! slo_ang
      izerrstat = 1
    ENDIF
    IF (.NOT. lsasp_lm) THEN
      ! slo_asp
      izerrstat = 1
    ENDIF
    IF (.NOT. lhori_lm) THEN
      ! horizon
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.7: Lake parameter
  !----------------------------------------------------------------------------

  ! Changes by Dmitrii:
  ! Setting the lake depth to 20 m everywhere will lead to an error in the COSMO 
  ! model. The lake depth should be set to -1 for land and sea grid boxes.
  ! Otherwise the lake parameterisation scheme FLake will be called at
  ! land/sea box, leading to an error.

  IF (llake) THEN
    IF (izdebug > 5) THEN
      PRINT *, ' Checks and calculations for FLake Model'
    ENDIF

    IF (.NOT. lflak_lm) THEN
      ! lake area fraction not available
      ! Stop the Program
      PRINT *, '*** ERROR: FR_LAKE is not found in external data file,'
      PRINT *, '***        Programm will be aborted!'
      izerrstat = 1
    ENDIF
    IF (.NOT. ldept_lm) THEN
      ! lake depth not available
      IF (lflak_lm) THEN
        ! If Lake fraction is available, set depth_lk accordingly
        IF (my_cart_id == 0) THEN
          PRINT *, 'WARNING: DEPTH_LK not found in external data file'
          PRINT *, 'set DEPTH_LK to 20 metres for lake grid boxes'
          PRINT *, 'set DEPTH_LK to -1 metre  for ocean or land grid boxes'
        ENDIF
        DO j = 1, je2lm
          DO i = 1, ie2lm
            IF(fr_lake_lm(i,j) > 0.5_ireals) THEN
              depth_lk_lm(i,j) = 20.0_ireals
            ELSE
              depth_lk_lm(i,j) = -1.0_ireals
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ELSE
      ! Check consistency of depth_lk and fr_lake
      ! NOTE: This check is only valid, if NO tile approach is used!!!
      IF (lflak_lm) THEN
        izerrflake = 0
        DO j = 1, je2lm
          DO i = 1, ie2lm
            IF ( (fr_lake_lm(i,j) < 0.5_ireals) .AND. (depth_lk_lm(i,j) >= 0.0_ireals) ) THEN
              izerrflake = 1
            ENDIF
          ENDDO
        ENDDO
        IF (izerrflake > 0) THEN
          PRINT *, ' *** ERROR *** : DEPTH_LK >= 0.0 is not consistent with FR_LAKE < 0.5 *** '
          izerrstat = 1
        ENDIF
      ENDIF
    ENDIF

    IF (.NOT. lsalt_lm) THEN
      ! lake salinity not available
      ! set to a default of 0
      IF (my_cart_id == 0) THEN
        PRINT *, 'WARNING: SALT_LK not found in external data file'
        PRINT *, 'set SALT_LK to default value 0'
      ENDIF
      salt_lk_lm(:,:) = 0.0_ireals
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.8: Climatological deep soil temperature
  !----------------------------------------------------------------------------

  IF (itype_t_cl == 1) THEN
    IF (.NOT. ltcl__lm) THEN
      ! climatological deep soil temperature not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.9: Thermal radiative surface emissivity
  !----------------------------------------------------------------------------

  IF (lemiss) THEN
    IF (.NOT. lemis_lm) THEN
      ! thermal radiative surface emissivity not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.10: Minimum stomata resistance of plants
  !----------------------------------------------------------------------------

  IF (lstomata) THEN
    IF (.NOT. lprsm_lm) THEN
      ! Minimum stomata resistance of plants not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.11.1: 12 monthly values of ndvi ratio
  !----------------------------------------------------------------------------

  IF (itype_ndvi == 1) THEN
    IF (.NOT. lndvi_lm) THEN
      ! 12 monthly values ndvi ratio not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.11.2: 12 monthly values of vegetation
  !----------------------------------------------------------------------------

  IF (itype_ndvi == 2) THEN
    IF (.NOT. lpl12_lm) THEN
      ! 12 monthly values of plant cover not available
      izerrstat = 1
    ENDIF
    IF (.NOT. lai12_lm) THEN
      ! 12 monthly values of leaf area index not available
      izerrstat = 1
    ENDIF
    IF (.NOT. lz012_lm) THEN
      ! 12 monthly values of z0 not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.12: 12 monthly values of aerosol values
  !----------------------------------------------------------------------------

  IF (itype_aerosol == 2) THEN
    IF (.NOT. lsu12_lm) THEN
      ! 12 monthly values of sulfate drops not available
      izerrstat = 1
    ENDIF
    IF (.NOT. ldu12_lm) THEN
      ! 12 monthly values of mineral dust not available
      izerrstat = 1
    ENDIF
    IF (.NOT. lor12_lm) THEN
      ! 12 monthly values of organic not available
      izerrstat = 1
    ENDIF
    IF (.NOT. lbc12_lm) THEN
      ! 12 monthly values of black carbon not available
      izerrstat = 1
    ENDIF
    IF (.NOT. lss12_lm) THEN
      ! 12 monthly values of sea salt not available
      izerrstat = 1
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5.13: Save input from external parameters for root depth
  !----------------------------------------------------------------------------

  rootdp_mx(:,:) = rootdp_lm(:,:)

  !----------------------------------------------------------------------------
  ! Section 5.14: Solar surface albedo
  !----------------------------------------------------------------------------

  IF     (itype_albedo == 2) THEN
    IF (.NOT. laldr_lm) THEN
      ! dry solar surface albedo not available
      izerrstat = 1
    ENDIF
    IF (.NOT. lalsa_lm) THEN
      ! saturated solar surface albedo not available
      izerrstat = 1
    ENDIF
  ELSEIF (itype_albedo == 3) THEN
    IF (.NOT. lal12_lm) THEN
      ! 12 monthly values of albedo not available
      izerrstat = 1
    ENDIF
  END IF

!------------------------------------------------------------------------------
! Section 6: Check and output, which fields have been read or interpolated
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, ' Checks which fields have been read / interpolated'
  ENDIF

  ! Write to noutput, which fields have been read or not read
  IF (my_cart_id == 0) THEN
    WRITE (noutput, '(A)')  '  '
    WRITE (noutput, '(A)')  '    The following external parameters have been read:'
    WRITE (noutput, '(A)')  '  '
    WRITE (noutput, '(T40,A)') ' COARSE     LM     INT2LM'
    WRITE (noutput, '(T8,A,T40,L6,T47,L6,T56,L6)')                         &
       'Height of the Orography     ',lfis__in, lhsur_lm, lhsur_int2lm
    WRITE (noutput, '(T8,A,T40,L6,T47,L6,T56,L6)')                         &
       'Geopotential of the surface ',lfis__in, lfis__lm, lfis__int2lm
    WRITE (noutput, '(T8,A,T40,L6,T47,A,T54,I2,A)')                        &
       'Orography transition at LB  ',llbc_smooth,' with ',nlbc_smooth,' grid points'
    WRITE (noutput, '(T8,A,T40,L6,T47,L6,T56,L6)')                         &
       'Fraction of land            ',lfrla_in, lfrla_lm, lfrla_int2lm
    WRITE (noutput, '(T8,A,T40,L6,T47,L6,T56,L6)')                         &
       'Soiltyp                     ',lsoty_in, lsoty_lm, lsoty_int2lm

    IF (itype_ndvi /= 2) THEN
      WRITE (noutput, '(T8,A,T40,L6,T47,L6,T56,L6)')                       &
       'Surface roughness           ',lz0___in, lz0___lm, lz0___int2lm
    ELSEIF (itype_ndvi == 2) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Surface roughness (12 month)', '---'  , lz012_lm, ' --- '       
    ENDIF

    WRITE (noutput, '(T8,A,T40,L6,T47,L6,T56,L6)')                         &
       'Root depth                  ',lroot_in, lroot_lm, lroot_int2lm

    IF (itype_ndvi /= 2) THEN
      WRITE (noutput, '(T8,A,T40,L6,T47,L6,T56,L6)')                       &
       'Plant Cover (Vegetation)    ',lplmx_in, lplmx_lm, lplmx_int2lm
      WRITE (noutput, '(T8,A,T40,L6,T47,L6,T56,L6)')                       &
       'Plant Cover (Rest)          ',lplmn_in, lplmn_lm, lplmn_int2lm
      WRITE (noutput, '(T8,A,T40,L6,T47,L6,T56,L6)')                       &
       'Leaf area index (Vegetation)',laimx_in, laimx_lm, laimx_int2lm
      WRITE (noutput, '(T8,A,T40,L6,T47,L6,T56,L6)')                       &
       'Leaf area index (Rest)      ',laimn_in, laimn_lm, laimn_int2lm
    ELSEIF (itype_ndvi == 2) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Plant Cover (12 month)      ', '---'  , lpl12_lm, ' --- '       
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Leaf area index (12 month)  ', '---'  , lai12_lm, ' --- '       
    ENDIF

    IF (itype_ndvi == 1) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'ndvi ratio                  ', ' --- ', lndvi_lm, ' --- '
    ENDIF

    IF (lforest) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Evergreen Forest            ', ' --- ', lfore_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Deciduous Forest            ', ' --- ', lford_lm, ' --- '
    ENDIF

    IF (lurban) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Urban land                  ', ' --- ', lurba_lm, ' --- '
    ENDIF

    IF (lemiss) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Surface Emissivity          ', ' --- ', lemis_lm, ' --- '
    ENDIF

    IF (lstomata) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Stomata Resistance          ', ' --- ', lprsm_lm, ' --- '
    ENDIF

    IF (lsso) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Standard deviation of SSO   ', ' --- ', lstdh_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Anisotropy of the orography ', ' --- ', lgamm_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Angle between axis          ', ' --- ', lthet_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Mean slope of SSO           ', ' --- ', lsigm_lm, ' --- '
    ENDIF

    IF (lradtopo) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Sky view                    ', ' --- ', lskyv_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Slope Angle                 ', ' --- ', lsang_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Slope Aspect                ', ' --- ', lsasp_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Horizon                     ', ' --- ', lhori_lm, ' --- '
    ENDIF

    IF (llake) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Lake area fraction          ', ' --- ', lflak_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Lake depth                  ', ' --- ', ldept_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Lake salinity               ', ' --- ', lsalt_lm, ' --- '
    ENDIF

    IF (itype_aerosol == 2) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Sulfat drops (12 month)     ', ' --- ', lsu12_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Mineral dust (12 month)     ', ' --- ', ldu12_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Organic      (12 month)     ', ' --- ', lor12_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Black carbon (12 month)     ', ' --- ', lbc12_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Sea salt     (12 month)     ', ' --- ', lss12_lm, ' --- '
    ENDIF

    IF     (itype_albedo == 2) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Dry soil surface albedo     ', ' --- ', laldr_lm, ' --- '
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Sat. soil surface albedo    ', ' --- ', lalsa_lm, ' --- '
    ELSEIF (itype_albedo == 3) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'Solar surface albedo (12 M) ', ' --- ', lal12_lm, ' --- '
    ENDIF

    IF (lcm2lm) THEN
      WRITE (noutput, '(T8,A,T43,A, T47,L6,T59,A )')                       &
       'clim. deep soil temperature ', ' --- ', ltcl__lm, ' --- '
    ENDIF
    WRITE (noutput, '(A)')  '  '
  ENDIF

  IF (izerrstat /= 0) THEN
    ! not all necessary fields could be read or interpolated: abort program
    ierror = 10
    yerror = 'Not all necessary fields could be read or interpolated'
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 7: Coordinate parameters and reference atmosphere for fine grid LM
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, ' Coordinate parameters and reference atmosphere for COSMO-Model'
  ENDIF

  ALLOCATE (zp0hl(ie2lm,je2lm,ke1lm),      zt0hl(ie2lm,je2lm,ke1lm),      &
            zhhl_read(ie2lm,je2lm,ke1lm),                   STAT=izerrstat)

  ! Allocate zhsurfs_lm in any case, because this is used as subroutine argument
  ALLOCATE (zhsurfs_lm     (ie2lm    , je2lm    , 2),    STAT=izerrstat)

  IF ((.NOT. lhhl_lm_read) .AND. (vcoord%ivctype == 3 .OR. vcoord%ivctype == 4)) THEN
    ! allocate and initialize fields for SLEVE
    ALLOCATE (zhsurf_lm_tot  (ie2lm_tot, je2lm_tot   ),                    &
              zhsurfs_lm_tot (ie2lm_tot, je2lm_tot, 2),    STAT=izerrstat)

    ! collect full topo from all PE's
    IF (num_compute == 1 ) THEN     ! we are running on one PE
       zhsurf_lm_tot(:,:) = hsurf_lm(:,:)
    ELSE                            ! we are running on multiple PE's
       CALL gather_field(hsurf_lm,      ie2lm,     je2lm,                  &
                         zhsurf_lm_tot, ie2lm_tot, je2lm_tot, 0, izerror)
    ENDIF

    IF (my_cart_id == 0) THEN
      CALL sleve_split_oro                                                  &
           (zhsurf_lm_tot, zhsurfs_lm_tot, ie2lm_tot, je2lm_tot, nfltvc,    &
            1, svc1, svc2, vcoord%vcflat, noutput, my_cart_id, izerror, yzerror)

      IF (izerror /= 0_iintegers) THEN
        ierror = 11
        yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
        RETURN
      ENDIF
    ENDIF

    ! distribute splitted topo to all PE's
    IF (num_compute == 1) THEN      ! we are running on one PE
      zhsurfs_lm(:,:,1) = zhsurfs_lm_tot(:,:,1)
      zhsurfs_lm(:,:,2) = zhsurfs_lm_tot(:,:,2)
    ELSE                            ! we are running on multiple PE's
      CALL distribute_field(zhsurfs_lm_tot(:,:,1), ie2lm_tot, je2lm_tot,     &
                            zhsurfs_lm(:,:,1),     ie2lm,     je2lm, 0, izerror)
      CALL distribute_field(zhsurfs_lm_tot(:,:,2), ie2lm_tot, je2lm_tot,     &
                            zhsurfs_lm(:,:,2),     ie2lm,     je2lm, 0, izerror)
    ENDIF
  ENDIF

  IF (lhhl_lm_read) THEN
    ! read HHL file for COSMO-Grid
    CALL org_read_hhl (zhhl_read, ie2lm, je2lm, ke1lm, ielm_tot, jelm_tot,   &
                       startlat_tot, startlon_tot, vcoord, isubpos,          &
                       ylmext_cat, ylm_hhl, 'lm')
    ! Now zhhl_read is the HHL field, with which we should go on, but it is
    ! missing the additional boundary line at each side.
    ! But in org_read_hhl we got the correct vertical coordinate parameters
    ! to compute it again: we only have to check that we use the correct
    ! orography:  hsurf_lm == zhhl_read(:,:,ke1lm)
    DO j = 2, je2lm-1
      DO i = 2, ie2lm-1
        IF ( ABS(hsurf_lm(i,j) - zhhl_read(i,j,ke1lm)) > 1E-6_ireals) THEN
          WRITE (*,'(A)') ' *** ERROR: The orography HSURF does not fit to HHL read!'
          WRITE (*,'(A,2I5,F18.8)') 'difference e.g.:  ', i,j, hsurf_lm(i,j) - zhhl_read(i,j,ke1lm)
          ierror = 12
          yerror = ' *** ERROR: The orography HSURF does not fit to HHL read!'
          RETURN
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  IF (refatm%irefatm == 1) THEN
    CALL reference_atmosphere                                                &
     ( hhl_lm, p0_lm, zp0hl, rho0_lm, t0_lm, zt0hl, dp0_lm, hsurf_lm,        &
       zhsurfs_lm, ie2lm, je2lm, kelm, refatm,                               &
       vcoord, svc1, svc2, r_d, g, lanalyt_calc_t0p0, .TRUE.)
  ELSE IF (refatm%irefatm == 2) THEN
    CALL reference_atmosphere_2                                              &
     ( hhl_lm, p0_lm, zp0hl, rho0_lm, t0_lm, zt0hl, dp0_lm, hsurf_lm,        &
       zhsurfs_lm, ie2lm, je2lm, kelm, refatm,                               &
       vcoord, svc1, svc2, r_d, g, .TRUE.)
  ENDIF

  IF ((.NOT. lhhl_lm_read) .AND. (vcoord%ivctype == 3 .OR. vcoord%ivctype == 4)) THEN
    DEALLOCATE (zhsurf_lm_tot, zhsurfs_lm_tot)
  ENDIF
  DEALLOCATE (zp0hl, zt0hl, zhsurfs_lm, zhhl_read)

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, ' Cleanup external_data'
  ENDIF

  IF (lgme2lm) THEN
    DEALLOCATE (fr_land_gme, z0_gme, root_gme, plcmx_gme, plcmn_gme,     &
                rlaimx_gme, rlaimn_gme)
  ELSE
    DEALLOCATE (fr_land_in, z0_in, plcov_in, plcmx_in, plcmn_in, &
                rlaimx_in, rlaimn_in, root_in)
  ENDIF

  IF ( (lchkin) .AND. (my_cart_id == 0) ) THEN
    CLOSE (nuchkdat, STATUS='KEEP')
  ENDIF

END SUBROUTINE external_data

!==============================================================================
