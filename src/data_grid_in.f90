!+ Variables necessary for the GME grid.
!==============================================================================

MODULE data_grid_in

!==============================================================================
!
! Description:
!  This module defines the variables necessary for the input grid.
!  Input grid can be the grid for GME, EC (European Centre), DM/SM or LM.
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
! V1_5         2007/07/09 Ulrich Schaettler
!  Added parameters necessary for non-equidistant grids, pole and date-line 
!  crossing
!  Added czmls_in, ke_soil_in for multi-layer soil model
!  Introduced nlevskip to  skip levels at the top of ECMWF model (Davide Cesari)
! V1_6         2007/09/07 Ulrich Schaettler
!  Editorial changes
! V1_9         2009/09/03 Guenther Zaengl
!  Additional parameters for new reference atmosphere of input model: 
!  irefatm_in, delta_t_in, h_scal_in
! V1_10        2009/12/17 Ulrich Schaettler
!  Added additional vectors ak_in_uv, bk_in_uv for lum2lm
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler
!  New variables for pressure levels
! V1_19        2012/06/06 Ulrich Schaettler, Daniel Luethi, Burkhardt Rockel
!  Renamed ak_in_uv, bk_in_uv to ak_in_rho, bk_in_rho according to UM conventions
!  Introduced new logical variables  (active for lcm2lm only)
!    - lcm_hgt_coor for hybrid height coordinates in coarse data input
!    - lcm_pres_coor for pressure coordinates in coarse data input
! V1_22        2013/07/11 Ulrich Schaettler
!  Eliminated values for input reference atmosphere and vertical grid
!   (moved to vgrid_refatm_utils)
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

! 1. GME Input grid
! -----------------

  ! Number of diamonds
  INTEGER (KIND=iintegers), PARAMETER      ::                      &
    ids=1 ,      & ! start index of diamonds (ids = 1)
    ide=10,      & ! end index of diamonds (ide = 10)
    nd=ide-ids+1   ! number of diamonds (nd = ide-ids+1 = 10)

  INTEGER (KIND=iintegers), PARAMETER      ::                      &
    isp11 =  1,  & ! Offsets of the 6 (5) neighbouring gridpoints relative to
    isp12 =  0,  & ! the central node for 2-dimensional array addressing. In
    isp13 = -1,  & ! i1-direction use j1 + isp1* (with * from 1 to 6), in
    isp14 = -1,  & ! i2-direction use j2 + isp2* (with * from 1 to 6) with
    isp15 =  0,  & ! (j1,j2) the indices of the central node.
    isp16 =  1,  & ! Attention:
    isp21 =  0,  & ! There are four points in the extended array which have
    isp22 =  1,  & ! other spokes ("ispokes") since they are close to the
    isp23 =  1,  & ! corners of the diamonds. For these points, the spokes 
    isp24 =  0,  & ! defined here are not valid.
    isp25 = -1,  & !
    isp26 = -1     !

  ! Resolution and vertical grid size
  INTEGER (KIND=iintegers)                 ::                      &
    ni_gme,      & ! resolution of GME
    i3e_gme,     & ! number of levels in the vertical
    ni2,         & ! ni_gme=3**ni3*2**ni2 with ni3 0 or 1 and ni2 > 1
    ni3,         & !
    i1mrp(8),    & ! i1-index of the four mirrored points of the extended array
    i2mrp(8),    & ! i2-index of the four mirrored points of the extended array
    ispoke(12),  & ! offsets of the 6 (5) neighbouring gridpoints relative to
                   ! the central node for 2-dimensional array addressing; in
                   ! i1-direction use ispoke(m), m=1,6 (5), in i2-direction
                   ! use ispoke(m+6), m=1,6 (5); phys. dim. ( - )
    ispokes(12,4)  ! offsets of the 6 neighbouring gridpoints relative
                   ! to the central node for the four special points

  REAL (KIND=ireals)                       ::                      &
    dxmin,       & ! minimum mesh width of the grid
    dhmin          ! minimum height of triangle; determines the timestep dt

  ! Decomposition for the parallel program
  INTEGER (KIND=iintegers), ALLOCATABLE    ::                      &
    ilim1(:),    & ! decomposition limits in x-direction (formerly 1)
    ilim2(:)       ! decomposition limits in y-direction (formerly 2)

  ! Variables for global horizontal dimensioning of GME fields 
  INTEGER (KIND=iintegers)                 ::                      &
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
    igg2ep2        ! = igg2e + 2

! Variables for local horizontal dimensioning of GME fields 
! and further help variables
  INTEGER (KIND=iintegers)                 ::                      &
    ig1s,        & ! start index of array-dimension in x-direction
    ig1sm1,      & ! = ig1s - 1
    ig1sm2,      & ! = ig1s - 2
    ig1e,        & ! end index of array-dimension in x-direction
    ig1ep1,      & ! = ig1e + 1
    ig1ep2,      & ! = ig1e + 2
    ig2s,        & ! start index of array-dimension in y-direction
    ig2sm1,      & ! = ig2s - 1
    ig2sm2,      & ! = ig2s - 2
    ig2e,        & ! end index of array-dimension in y-direction
    ig2ep1,      & ! = ig2e + 1
    ig2ep2         ! = ig2e + 2
 
! Start and end-indices for computations during the horizontal interpolations
! that are necessary for a LM subdomain
  INTEGER (KIND=iintegers)                 ::                      &
    j1_min,      & ! smallest index in j1-direction for LM subdomain
    j1_max,      & ! biggest index in j1-direction for LM subdomain
    j2_min,      & ! smallest index in j2-direction for LM subdomain
    j2_max,      & ! biggest index in j2-direction for LM subdomain
    jd_min,      & ! smallest index for diamonds for a LM subdomain
    jd_max         ! biggest index for diamonds for a LM subdomain

! IDs of the neighbors for GME subdomains and sizes of GME core
  INTEGER (KIND=iintegers)                 ::                      &
    nbpe,        & ! Number of neighbor poleward east
    nbaw,        & ! Number of neighbor antipoleward west
    nbpw,        & ! Number of neighbor poleward west
    nbae,        & ! Number of neighbor antipoleward east
    max_gme_core   ! size of global GME core

!==============================================================================

! 2. Other (regular) input grids
! ------------------------------

  ! Rotation and position of the grid on the earth
  REAL (KIND=ireals)          ::                                   &
    pollat_in,      & ! latitude of the rotated north pole (in degrees, N>0)
    pollon_in,      & ! longitude of the rotated north pole (in degrees, E>0)
    polgam_in,      & ! latitude of the rotated north pole
    dlat_in,        & ! grid point distance in zonal direction (in degrees)
    dlon_in,        & ! grid point distance in meridional direction
    startlat_in_tot,& ! transformed latitude of the lower left grid point
                      ! of the total domain (in degrees, N>0)
    startlon_in_tot,& ! transformed longitude of the lower left grid point
                      ! of the total domain (in degrees, E>0)
    endlat_in_tot,  & ! transformed latitude of the upper right grid point
                      ! of the total domain (in degrees, N>0)
    endlon_in_tot,  & ! transformed longitude of the upper right grid point
                      ! of the total domain (in degrees, E>0)
    startlat_in,    & ! transformed latitude of the lower left grid point
                      ! of the local domain (in degrees, N>0)
    startlon_in,    & ! transformed longitude of the lower left grid point
                      ! of the local domain (in degrees, E>0)
    endlat_in,      & ! transformed latitude of the upper right grid point
                      ! of the local domain (in degrees, N>0)
    endlon_in         ! transformed longitude of the upper right grid point
                      ! of the local domain (in degrees, E>0)

  ! horizontal size of the grid
  INTEGER (KIND=iintegers)    ::                                   &
    ie_in_tot,   & ! ie for input grid, total domain
    je_in_tot,   & ! je for input grid, total domain
    ke_in_tot,   & ! ke for input grid
    ie_in_max,   & ! Max. of ie_in (local) on all processors
    je_in_max,   & ! Max. of je_in (local) on all processors
    ie_in,       & ! ie for input grid, local domain
    je_in,       & ! je for input grid, local domain
    ke_in,       & ! ke for input grid
    ke1in,       & !
    ke_hybrid,   & ! number of hybrid levels when interpolating from pressure levels
    ke_pressure, & ! number of pressure levels when interpolating from pressure levels
    kedim_in,    & ! MAX (ke_in, ke_hybrid)
    nlevskip       ! number of levels to skip at top of input model

  ! horizontal size of the grid
  INTEGER (KIND=iintegers)    ::                                   &
    ie_ext_in,   & ! west-east size of fields with external parameters
    je_ext_in      ! north-south size of fields with external parameters
                   ! these values may be larger than the actual sizes

  ! actual vertical coordinate parameters
  REAL (KIND=ireals), ALLOCATABLE     ::       &
    ak_in(:),         & ! coefficients for half levels
    bk_in(:),         & ! coefficients for half levels
    akh_in(:),        & ! coefficients for main levels
    bkh_in(:),        & ! coefficients for main levels
    akh_in_rho(:),    & ! coefficients for main levels for u, v (lum2lm)
    bkh_in_rho(:),    & ! coefficients for main levels for u, v (lum2lm)
    dak_in(:),        & ! differences between half levels
    dbk_in(:)           !                  - " -

  INTEGER (KIND=iintegers)    ::                                   &
    klv850_in           ! approximate level where 850 hPa is reached

  REAL (KIND=ireals)          ::                                   &
    press_level(50)     ! list of available pressure levels (in Pa) for GFS
                        ! with fixed length of 50, because we can read that as
                        ! a Namelist variable then

  ! Additional parameters for CLM use
  !----------------------------------

  INTEGER (KIND=iintegers)    ::                                   &
    grdpt_rel_in   ! indicator for the relations between the coarse grid and
                   ! the LM grid points to compute correctly the distance
                   ! functions if a Cressman-type scheme is used for M-type
                   ! interpolation

  ! actual horizontal coordinate parameters
  LOGICAL                     ::                                   &
    lushift_in(2),     & ! shift of u-velocity due to grid staggering
    lvshift_in(2)        ! shift of v-velocity due to grid staggering

  REAL (KIND=ireals), POINTER ::                                   &
    latitudes_in(:),   & ! latitudes of the input data
    longitudes_in(:),  & ! longitudes of the input data
    slatitudes_in(:),  & ! shifted latitudes of the input data
    slongitudes_in(:)    ! shifted longitudes of the input data

  INTEGER (KIND=iintegers)    ::                                   &
    east_add_in,       & ! add an extra column to the East
    west_add_in,       & ! add an extra column to the West
    south_add_in,      & ! add an extra column to the South
    north_add_in         ! add an extra column to the North

  LOGICAL                     ::                                   &
    lcm_hgt_coor =.FALSE., & ! switch to declare hybrid_height_coordinates on CM input
                             ! (for HadGEM2 height_coordinates, similar to NWP lum2lm)
    lcm_pres_coor=.FALSE.    ! switch to declare pressure coordinates on CM input


  ! Variables for specifying multi-layer soil model
  ! --------------------------------------------------
  REAL      (KIND=ireals),    ALLOCATABLE       ::           &
    czmls_in(:),  & ! depth of the input soil layers in meters
    czhls_in(:)     ! depth of the input half soil layers in meters

  INTEGER   (KIND=iintegers)       ::           &
    ke_soil_in        ! number of input levels in multi-layer soil model

!==============================================================================

END MODULE data_grid_in
