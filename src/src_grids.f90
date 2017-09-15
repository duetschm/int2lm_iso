!+ Source module for generating correspondence between coarse and LM grid
!==============================================================================

MODULE src_grids

!==============================================================================
!
! Description:
!   This module determines the connection between the LM- and the GME-grid
!   or the LM and a regular grid (LM, IFS) by computing
!     in case of GME
!       - the horizontal triangular grid mesh of GME
!       - the gradient and Laplacian operator
!       - the indices and barycentric coordinates
!     in case of another regular grid
!       - the i- and j-indices of the coarse mesh grid point which is lower
!         left to a given LM (fine mesh) grid point
!       - the relative distances between x- (i-) and y- (j-) coarse mesh and
!         fine mesh grid points
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
! V1_5         2007/07/09 Friedrich Theunert
!  Take into account that the domain can cross the date line
! V1_6         2007/09/07 Burkhardt Rockel, Ulrich Schaettler
!  Added option lcm2lm and SR coor_cm_lm to interpolate data from climate model
! V1_8         2008/05/29 Ulrich Schaettler
!  Eliminated ldebug
! V1_10        2009/12/17 Ulrich Schaettler
!  Modifications for treatment of Unified Model data
!  Implemented decomposition independent computation of the interpolation weights
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_12        2010/06/14 Ulrich Schaettler
!  Bug correction in coor_coarse_grid_lm (reported by Stephan Pfahl,ETH)
! V1_14        2010/11/19 Ulrich Schaettler
!  Modifications to allow processing of JMA and NCEP data
!  Bug Fix for computation of interpolation weights
!  Do all computations of lat and lon for the coarse grid starting with 
!     startlat_tot_in and startlon_tot_in
! V1_15        2010/12/10 Ulrich Blahak
!  Bug Fix for checking the interpolation y-weights, if they are 0
! V1_17        2011/03/11 Ulrich Schaettler
!  Check correctness of INT-function when computing interpolation weights
! V1_19        2012/06/06 Ulrich Schaettler
!  Compute interpolation weights also for HIRLAM data
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
USE data_parameters , ONLY :   &
    ireals,     & ! kind-type parameter for "normal" integer variables
    iintegers     ! KIND-type parameters for real variables

!------------------------------------------------------------------------------

USE data_grid_lm,     ONLY :   &
    pollat,      & ! latitude of the rotated north pole (in degrees, N>0)
    pollon,      & ! longitude of the rotated north pole (in degrees, E>0)
    polgam,      & ! latitude of the rotated north pole   !_br
    dlat,        & ! grid point distance in zonal direction (in degrees)
    dlon,        & ! grid point distance in meridional direction (in degrees)
    startlat_tot,& ! transformed latitude of the lower left grid point
                   ! of the total domain (in degrees, N>0)
    startlon_tot,& ! transformed longitude of the lower left grid point
                   ! of the total domain (in degrees, E>0)
    startlat,    & ! transformed latitude of the lower left grid point
                   ! of the total domain (in degrees, N>0)
    startlon,    & ! transformed longitude of the lower left grid point
                   ! of the total domain (in degrees, E>0)
    ielm,        & ! ie for LM, local domain
    jelm,        & ! je for LM, local domain
    kelm,        & ! ke for LM
    ie2lm_max,   & ! Max. of ie2lm on all processors
    je2lm_max,   & ! Max. of je2lm on all processors
    ke1lm,       & !
    ie2lm,       & !
    je2lm          !

!------------------------------------------------------------------------------

USE data_grid_in,     ONLY :   &
    nd,          & ! number of diamonds
    ni_gme,      & ! resolution of GME
    ni2,         & ! ni_gme=3**ni3*2**ni2 with ni3 0 or 1 and ni2 > 1
    ni3,         & ! 
    dxmin,       & ! minimum mesh width of the grid
    dhmin,       & ! minimum height of triangle; determines the timestep dt
    i1mrp,       & ! i1-index of the four mirrored points of the extended array
    i2mrp,       & ! i2-index of the four mirrored points of the extended array
    ispoke,      & ! offsets of the 6 (5) neighbouring gridpoints relative to
                   ! the central node for 2-dimensional array addressing; in
                   ! i1-direction use ispoke(m), m=1,6 (5), in i2-direction 
                   ! use ispoke(m+6), m=1,6 (5); phys. dim. ( - )
    ispokes,     & ! offsets of the 6 neighbouring gridpoints relative
                   ! to the central node for the four special points
    isp11,       & ! Offsets of the 6 (5) neighbouring gridpoints relative to
    isp12,       & ! the central node for 2-dimensional array addressing. In
    isp13,       & ! i1-direction use j1 + isp1* (with * from 1 to 6), in
    isp14,       & ! i2-direction use j2 + isp2* (with * from 1 to 6) with
    isp15,       & ! (j1,j2) the indices of the central node.
    isp16,       & ! Attention:
    isp21,       & ! There are four points in the extended array which have
    isp22,       & ! other spokes ("ispokes") since they are close to the
    isp23,       & ! corners of the diamonds. For these points, the spokes 
    isp24,       & ! defined here are not valid.
    isp25,       & !
    isp26,       & !
    ilim1,       & ! decomposition limits in x-direction (formerly 1)
    ilim2,       & ! decomposition limits in y-direction (formerly 2)
    igg1s,       & ! start index of global array-dimension in x-direction
    igg1e,       & ! end index of global array-dimension in x-direction
    igg2s,       & ! start index of global array-dimension in y-direction
    igg2e,       & ! end index of global array-dimension in y-direction
    ig1s,        & ! start index of array-dimension in x-direction
    ig1sm1,      & ! = ig1s - 1
    ig1sm2,      & ! = ig1s - 2
    ig1e,        & ! end index of array-dimension in x-direction
    ig1ep1,      & ! = ig1e + 1
    ig1ep2         ! = ig1e + 2

USE data_grid_in,     ONLY :   &
    ig2s,        & ! start index of array-dimension in y-direction
    ig2sm1,      & ! = ig2s - 1
    ig2sm2,      & ! = ig2s - 2
    ig2e,        & ! end index of array-dimension in y-direction
    ig2ep1,      & ! = ig2e + 1
    ig2ep2,      & ! = ig2e + 2
    j1_min,      & ! smallest index in j1-direction for LM subdomain
    j1_max,      & ! biggest index in j1-direction for LM subdomain
    j2_min,      & ! smallest index in j2-direction for LM subdomain
    j2_max,      & ! biggest index in j2-direction for LM subdomain
    jd_min,      & ! smallest index for diamonds for a LM subdomain
    jd_max,      & ! biggest index for diamonds for a LM subdomain
    nbpe,        & ! Number of neighbor poleward east
    nbaw,        & ! Number of neighbor antipoleward west
    nbpw,        & ! Number of neighbor poleward west
    nbae,        & ! Number of neighbor antipoleward east
    max_gme_core,& ! size of global GME core
    pollat_in,   & ! latitude of the rotated north pole (in degrees, N>0)
    pollon_in,   & ! longitude of the rotated north pole (in degrees, E>0)
    polgam_in,   & ! latitude of the rotated north pole  !_br
    dlat_in,     & ! grid point distance in zonal direction (in degrees)
    dlon_in,     & ! grid point distance in meridional direction
    startlat_in_tot,& ! transformed latitude of the lower left grid point
                      ! of the total domain (in degrees, N>0)
    startlon_in_tot,& ! transformed longitude of the lower left grid point
                      ! of the total domain (in degrees, E>0)
    startlat_in, & ! transformed latitude of the lower left grid point
                   ! of the total domain (in degrees, N>0)
    startlon_in    ! transformed longitude of the lower left grid point
                   ! of the total domain (in degrees, E>0)

! for cm2lm:
USE data_grid_in,     ONLY :   &
    ie_in,       & ! ie for input grid, local domain
    je_in,       & ! je for input grid, local domain
    lushift_in,  & ! shift of u-velocity due to grid staggering
    lvshift_in,  & ! shift of v-velocity due to grid staggering
    latitudes_in,& ! latitudes of the input data
    longitudes_in,&! longitudes of the input data
    slatitudes_in,&! shifted latitudes of the input data
    slongitudes_in ! shifted longitudes of the input data

!------------------------------------------------------------------------------

USE data_fields_lm,   ONLY :   &
    latlm_m,     & ! latitudes of the LM grid points
    lonlm_m,     & ! longitudes of the LM grid points
    latlm_u,     & ! latitudes of the LM u grid points
    lonlm_u,     & ! longitudes of the LM u grid points
    latlm_v,     & ! latitudes of the LM v grid points
    lonlm_v,     & ! longitudes of the LM v grid points
    index_m,     & !
    index_u,     & !
    index_v,     & !
    baryll_m,    & !
    baryll_u,    & !
    baryll_v,    & !
    rotang_m,    & !
    rotang_u,    & !
    rotang_v,    & !
    rpoints_m,   & !
    rpoints_u,   & !
    rpoints_v,   & !
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
    grd_glob,    & ! coefficients needed for the calculation of the 
                   ! gradients in eta- and chi-direction               (1/m)
    eta_glob,    & !
    chi_glob,    & !
    cpsi_glob,   & !
    spsi_glob,   & !
    lat_coarse_m,& ! latitudes of the LM grid points
    lon_coarse_m,& ! longitudes of the LM grid points
    lat_coarse_u,& ! latitudes of the LM u grid points
    lon_coarse_u,& ! longitudes of the LM u grid points
    lat_coarse_v,& ! latitudes of the LM v grid points
    lon_coarse_v   ! longitudes of the LM v grid points

!------------------------------------------------------------------------------

USE data_int2lm_control,     ONLY :   &
    noutput,      & ! unit number for output file
    lgme2lm,      & ! if .TRUE., gme->lm
    lgfs2lm,      & ! if .TRUE., gfs->lm
    lgsm2lm,      & ! if .TRUE., gsm->lm
    lec2lm,       & ! if .TRUE., ec ->lm
    llm2lm,       & ! if .TRUE., lm ->lm
    lum2lm,       & ! if .TRUE., um ->lm
    lhir2lm,      & ! if .TRUE., hirlam ->lm
    lcm2lm          ! if .TRUE., cm ->lm

!------------------------------------------------------------------------------

USE data_profiles,    ONLY :   &
    lprgp,        & ! logical for print at selected grid points
    ngp,          & ! number of selected grid points in this local domain
    igp,          & ! i-indices    local domain
    jgp             ! j-indices    local domain

!------------------------------------------------------------------------------

USE data_int2lm_parallel,    ONLY :   &
    nprocx,          & ! number of processors in x-direction
    nprocy,          & ! number of processors in y-direction
    nproc1,          & ! number of processors in first direction for GME
                       ! = nprocy (!)
    nproc2,          & ! number of processors in second direction for GME
                       ! = nprocx (!)
    num_compute,     & ! number of compute PEs

    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains

    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_pos,     & ! position of this subdomain in the cartesian grid 
                       ! in x- and y-direction
    my_num1,         & ! position of this subdomain in the cartesian grid
    my_num2,         & ! in first and second direction
                       ! (my_num1 = my_cart_pos(2) and my_num2 = my_cart_pos(1))
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    isubpos,         & ! positions of the subdomains in the total domain.
    isubpos_coarse,  & ! positions of the subdomains in the total domain for 
                       ! coarse grid
    icomm_cart,      & ! communicator for the virtual cartesian topology
    igroup_world,    & ! group that belongs to MPI_COMM_WORLD, i.e. all 
                       ! processors
    icomm_world,     & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_integers,    & ! determines the correct INTEGER type used in the
                       ! model for MPI
    imp_character,   & ! determines the correct CHARACTER type used in the
                       ! model for MPI
    imp_logical        ! determines the correct LOGICAL   type used in the
                       ! model for MPI

!------------------------------------------------------------------------------

USE data_int2lm_constants,   ONLY :   &
    r_earth, & ! mean radius of the earth          [m]
    Pid5,    & !
    Omcor      !

!------------------------------------------------------------------------------

USE parallel_utilities, ONLY :  global_values
USE gme_utilities,      ONLY :  gather_gme_field, xd_p, pp_ll2gp, distance, &
                                factorni, gcpt, tricntr
USE utilities,          ONLY :  phirot2phi, rlarot2rla, phi2phirot, rla2rlarot

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Public and Private Subroutines

PUBLIC   coor_gme_lm, coor_coarse_grid_lm, compute_geo_coord, coor_cm_lm

PRIVATE  pr_genhor, glo_coor, dxdhmin, loc_coor, tri_area, hex_area,     &
         gen_grd, setspokes, setmrp, pp_horint_index_p

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure for setting up the correspondence of LM grid to coarse grid
!------------------------------------------------------------------------------

SUBROUTINE coor_coarse_grid_lm (yerror, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the i- and j-indices of the coarse mesh grid point 
!   which is lower left to a given LM (fine mesh) grid point and
!   the relative distances between x- (i-) and y- (j-) coarse mesh and
!   fine mesh grid points
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist
CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local variables
INTEGER (KIND=iintegers) :: &
  i, j, izerror, iind_total, jind_total

REAL    (KIND=ireals)    :: &
  zlats, zlons, zdifflon, zdifflat, zdifflon_tot, zlats_tot, zlons_tot, z1, z2

CHARACTER (LEN= 80)      ::  &
  yzerror       ! error message

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine coor_coarse_grid_lm
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierror  = 0_iintegers
  yerror  = '         '
  izerror = 0_iintegers
  yzerror = '         '

!------------------------------------------------------------------------------
! Section 2: Generation of coarse and fine grid correspondence
!------------------------------------------------------------------------------

  ! For every LM (fine mesh) grid point the grid point of the coarse grid
  ! is determined that is just to the lower left of the LM grid point.
  ! It has to be taken into account that the domain could cross the date line
  ! (startlon_in > endlon_in)
  ! This is done also for the staggered positions, which gives 5 possibilities:

  !  1:  mass grid points:
  DO j = 1, je2lm
    DO i = 1, ie2lm
      ! index of lower left grid point in the coarse mesh
      zdifflon = lon_coarse_m(i,j) - startlon_in
      IF (zdifflon < 0.0_ireals) zdifflon = zdifflon + 360.0_ireals

      z1 = zdifflon / dlon_in
      ! if z1 is in the proximity of an integer, the INT function might have
      ! problems to determine the exact integer number
      IF ( ABS(NINT(z1) - z1) < 1E-6_ireals) THEN
        i_index(i,j,1) = NINT(z1) + 1
      ELSE
        i_index(i,j,1) = INT (z1) + 1
      ENDIF

      z2 = (lat_coarse_m(i,j) - startlat_in) / dlat_in
      IF ( ABS(NINT(z2) - z2) < 1E-6_ireals) THEN
        j_index(i,j,1) = NINT(z2) + 1
      ELSE
        j_index(i,j,1) = INT (z2) + 1
      ENDIF

      ! weights for the interpolation: to be decomposition independent, the
      ! weights are calculated with respect to the total coarse grid

      ! first get lon-index of this grid point in the total grid
      zdifflon_tot   = lon_coarse_m(i,j) - startlon_in_tot
      IF (zdifflon_tot < 0.0_ireals) zdifflon_tot = zdifflon_tot + 360.0_ireals
      z1 = zdifflon_tot / dlon_in
      IF ( ABS(NINT(z1) - z1) < 1E-6_ireals) THEN
        iind_total = NINT(z1) + 1
      ELSE
        iind_total = INT (z1) + 1
      ENDIF

      ! compute weight for longitudes
      zdifflon       = lon_coarse_m(i,j) - longitudes_in(iind_total)
      ! it could be that these grid points have the same longitude, which
      ! could lead to a numerical small negative value instead of 0.0
      IF     (ABS(zdifflon) <   1.0E-6_ireals) THEN
        ! fine grid point and coarse grid point have the same longitude
        zdifflon = 0.0_ireals
      ELSEIF (    zdifflon  <= -1.0E-6_ireals) THEN
        ! should only be the case when crossing the date line
        zdifflon = zdifflon + 360.0
      ENDIF
      x_wght (i,j,1) = zdifflon / dlon_in

      ! get lat-index of this grid point in the total grid
      z2 = (lat_coarse_m(i,j)-startlat_in_tot)/ dlat_in
      IF ( ABS(NINT(z2) - z2) < 1E-6_ireals) THEN
        jind_total = NINT(z2) + 1
      ELSE
        jind_total = INT (z2) + 1
      ENDIF
      ! compute weight for latitudes
      zdifflat       = lat_coarse_m(i,j) - latitudes_in(jind_total)
      IF (ABS(zdifflat) < 1.0E-6_ireals) THEN
        y_wght(i,j,1) = 0.0_ireals
      ELSE
        y_wght(i,j,1) = zdifflat / dlat_in
      END IF
    ENDDO
  ENDDO
     
  !  2:  fine mesh u grid point and coarse mesh u grid point
  IF (lec2lm .OR. lgsm2lm .OR. lgfs2lm) THEN
    zlats     = startlat_in
    zlons     = startlon_in
    zlats_tot = startlat_in_tot
    zlons_tot = startlon_in_tot
  ELSEIF (llm2lm .OR. lum2lm .OR. lhir2lm) THEN
    zlats     = startlat_in
    zlons     = startlon_in + 0.5 * dlon_in
    zlats_tot = startlat_in_tot
    zlons_tot = startlon_in_tot + 0.5 * dlon_in
  ENDIF

  DO j = 1, je2lm
    DO i = 1, ie2lm
      zdifflon = lon_coarse_u(i,j) - zlons
      IF (zdifflon < 0.0_ireals) zdifflon = zdifflon + 360.0_ireals

      z1 = zdifflon / dlon_in
      ! if z1 is in the proximity of an integer, the INT function might have
      ! problems to determine the exact integer number
      IF ( ABS(NINT(z1) - z1) < 1E-6_ireals) THEN
        i_index(i,j,2) = NINT(z1) + 1
      ELSE
        i_index(i,j,2) = INT (z1) + 1
      ENDIF

      z2 = (lat_coarse_u(i,j) - zlats) / dlat_in
      IF ( ABS(NINT(z2) - z2) < 1E-6_ireals) THEN
        j_index(i,j,2) = NINT(z2) + 1
      ELSE    
        j_index(i,j,2) = INT (z2) + 1
      ENDIF

      ! weights for the interpolation: to be decomposition independent, the
      ! weights are calculated with respect to the total coarse grid

      ! first get lon-index of this grid point in the total grid
      zdifflon_tot   = lon_coarse_u(i,j) - zlons_tot
      IF (zdifflon_tot < 0.0_ireals) zdifflon_tot = zdifflon_tot + 360.0_ireals
      z1 = zdifflon_tot / dlon_in
      IF ( ABS(NINT(z1) - z1) < 1E-6_ireals) THEN
        iind_total = NINT(z1) + 1
      ELSE
        iind_total = INT (z1) + 1
      ENDIF

      ! compute weight for longitudes
      zdifflon       = lon_coarse_u(i,j) - slongitudes_in(iind_total)
      ! it could be that these grid points have the same longitude, which
      ! could lead to a numerical small negative value instead of 0.0
      IF     (ABS(zdifflon) <   1.0E-6_ireals) THEN
        ! fine grid point and coarse grid point have the same longitude
        zdifflon = 0.0_ireals
      ELSEIF (    zdifflon  <= -1.0E-6_ireals) THEN
        ! should only be the case when crossing the date line
        zdifflon = zdifflon + 360.0
      ENDIF
      x_wght (i,j,2) = zdifflon / dlon_in

      ! get lat-index of this grid point in the total grid
      z2 = (lat_coarse_u(i,j)-zlats_tot)/ dlat_in
      IF ( ABS(NINT(z2) - z2) < 1E-6_ireals) THEN
        jind_total = NINT(z2) + 1
      ELSE
        jind_total = INT (z2) + 1
      ENDIF
      ! compute weight for latitudes
      zdifflat       = lat_coarse_u(i,j) - latitudes_in(jind_total)
      IF (ABS(zdifflat) < 1.0E-6_ireals) THEN
        y_wght(i,j,2) = 0.0_ireals
      ELSE
        y_wght(i,j,2) = zdifflat / dlat_in
      END IF
    ENDDO
  ENDDO
     
  !  3:  fine mesh u grid point and coarse mesh v grid point
  IF (lec2lm .OR. lgsm2lm .OR. lgfs2lm) THEN
    zlats = startlat_in
    zlons = startlon_in
    zlats_tot = startlat_in_tot
    zlons_tot = startlon_in_tot
  ELSEIF (llm2lm .OR. lum2lm .OR. lhir2lm) THEN
    zlats = startlat_in + 0.5 * dlat_in
    zlons = startlon_in
    zlats_tot = startlat_in_tot + 0.5 * dlat_in
    zlons_tot = startlon_in_tot
  ENDIF

  DO j = 1, je2lm
    DO i = 1, ie2lm
      zdifflon = lon_coarse_u(i,j) - zlons
      IF (zdifflon < 0.0_ireals) zdifflon = zdifflon + 360.0

      z1 = zdifflon / dlon_in
      ! if z1 is in the proximity of an integer, the INT function might have
      ! problems to determine the exact integer number
      IF ( ABS(NINT(z1) - z1) < 1E-6_ireals) THEN
        i_index(i,j,3) = NINT(z1) + 1
      ELSE
        i_index(i,j,3) = INT (z1) + 1
      ENDIF

      z2 = (lat_coarse_u(i,j) - zlats) / dlat_in
      IF ( ABS(NINT(z2) - z2) < 1E-6_ireals) THEN
        j_index(i,j,3) = NINT(z2) + 1
      ELSE
        j_index(i,j,3) = INT (z2) + 1
      ENDIF

      ! weights for the interpolation: to be decomposition independent, the
      ! weights are calculated with respect to the total coarse grid

      ! first get lon-index of this grid point in the total grid
      zdifflon_tot   = lon_coarse_u(i,j) - zlons_tot
      IF (zdifflon_tot < 0.0_ireals) zdifflon_tot = zdifflon_tot + 360.0
      z1 = zdifflon_tot / dlon_in
      IF ( ABS(NINT(z1) - z1) < 1E-6_ireals) THEN
        iind_total = NINT(z1) + 1
      ELSE
        iind_total = INT (z1) + 1
      ENDIF

      ! compute weight for longitudes
      zdifflon       = lon_coarse_u(i,j) - longitudes_in(iind_total)
      ! it could be that these grid points have the same longitude, which
      ! could lead to a numerical small negative value instead of 0.0
      IF     (ABS(zdifflon) <   1.0E-6_ireals) THEN
        ! fine grid point and coarse grid point have the same longitude
        zdifflon = 0.0_ireals
      ELSEIF (    zdifflon  <= -1.0E-6_ireals) THEN
        ! should only be the case when crossing the date line
        zdifflon = zdifflon + 360.0
      ENDIF
      x_wght (i,j,3) = zdifflon / dlon_in

      ! get lat-index of this grid point in the total grid
      z2 = (lat_coarse_u(i,j)-zlats_tot)/ dlat_in
      IF ( ABS(NINT(z2) - z2) < 1E-6_ireals) THEN
        jind_total = NINT(z2) + 1
      ELSE
        jind_total = INT (z2) + 1
      ENDIF
      ! compute weight for latitudes
      zdifflat       = lat_coarse_u(i,j) - slatitudes_in(jind_total)
      IF (ABS(zdifflat) < 1.0E-6_ireals) THEN
        y_wght(i,j,3) = 0.0_ireals
      ELSE
        y_wght(i,j,3) = zdifflat / dlat_in
      END IF
    ENDDO
  ENDDO

  !  4:  fine mesh v grid point and coarse mesh v grid point
  IF (lec2lm .OR. lgsm2lm .OR. lgfs2lm) THEN
    zlats = startlat_in
    zlons = startlon_in
    zlats_tot = startlat_in_tot
    zlons_tot = startlon_in_tot
  ELSEIF (llm2lm .OR. lum2lm .OR. lhir2lm) THEN
    zlats = startlat_in + 0.5 * dlat_in
    zlons = startlon_in
    zlats_tot = startlat_in_tot + 0.5 * dlat_in
    zlons_tot = startlon_in_tot
  ENDIF

  DO j = 1, je2lm
    DO i = 1, ie2lm
      zdifflon = lon_coarse_v(i,j) - zlons
      IF (zdifflon < 0.0_ireals) zdifflon = zdifflon + 360.0

      z1 = zdifflon / dlon_in
      ! if z1 is in the proximity of an integer, the INT function might have
      ! problems to determine the exact integer number
      IF ( ABS(NINT(z1) - z1) < 1E-6_ireals) THEN
        i_index(i,j,4) = NINT(z1) + 1
      ELSE
        i_index(i,j,4) = INT (z1) + 1
      ENDIF

      z2 = (lat_coarse_v(i,j) - zlats) / dlat_in
      IF ( ABS(NINT(z2) - z2) < 1E-6_ireals) THEN
        j_index(i,j,4) = NINT(z2) + 1
      ELSE    
        j_index(i,j,4) = INT (z2) + 1
      ENDIF

      ! weights for the interpolation: to be decomposition independent, the
      ! weights are calculated with respect to the total coarse grid

      ! first get lon-index of this grid point in the total grid
      zdifflon_tot   = lon_coarse_v(i,j) - zlons_tot
      IF (zdifflon_tot < 0.0_ireals) zdifflon_tot = zdifflon_tot + 360.0
      z1 = zdifflon_tot / dlon_in
      IF ( ABS(NINT(z1) - z1) < 1E-6_ireals) THEN
        iind_total = NINT(z1) + 1
      ELSE
        iind_total = INT (z1) + 1
      ENDIF

      ! compute weight for longitudes
      zdifflon       = lon_coarse_v(i,j) - longitudes_in(iind_total)
      ! it could be that these grid points have the same longitude, which
      ! could lead to a numerical small negative value instead of 0.0
      IF     (ABS(zdifflon) <   1.0E-6_ireals) THEN
        ! fine grid point and coarse grid point have the same longitude
        zdifflon = 0.0_ireals
      ELSEIF (    zdifflon  <= -1.0E-6_ireals) THEN
        ! should only be the case when crossing the date line
        zdifflon = zdifflon + 360.0
      ENDIF
      x_wght (i,j,4) = zdifflon / dlon_in

      ! get lat-index of this grid point in the total grid
      z2 = (lat_coarse_v(i,j)-zlats_tot)/ dlat_in
      IF ( ABS(NINT(z2) - z2) < 1E-6_ireals) THEN
        jind_total = NINT(z2) + 1
      ELSE
        jind_total = INT (z2) + 1
      ENDIF
      ! compute weight for latitudes
      zdifflat       = lat_coarse_v(i,j) - slatitudes_in(jind_total)
      IF (ABS(zdifflat) < 1.0E-6_ireals) THEN
        y_wght(i,j,4) = 0.0_ireals
      ELSE
        y_wght(i,j,4) = zdifflat / dlat_in
      END IF
    ENDDO
  ENDDO

  !  5:  fine mesh v grid point and coarse mesh u grid point
  IF (lec2lm .OR. lgsm2lm .OR. lgfs2lm) THEN
    zlats = startlat_in
    zlons = startlon_in
    zlats_tot = startlat_in_tot
    zlons_tot = startlon_in_tot
  ELSEIF (llm2lm .OR. lum2lm .OR. lhir2lm) THEN
    zlats = startlat_in
    zlons = startlon_in + 0.5 * dlon_in
    zlats_tot = startlat_in_tot
    zlons_tot = startlon_in_tot + 0.5 * dlon_in
  ENDIF

  DO j = 1, je2lm
    DO i = 1, ie2lm
      zdifflon = lon_coarse_v(i,j) - zlons
      IF (zdifflon < 0.0_ireals) zdifflon = zdifflon + 360.0

      z1 = zdifflon / dlon_in
      ! if z1 is in the proximity of an integer, the INT function might have
      ! problems to determine the exact integer number
      IF ( ABS(NINT(z1) - z1) < 1E-6_ireals) THEN
        i_index(i,j,5) = NINT(z1) + 1
      ELSE
        i_index(i,j,5) = INT (z1) + 1
      ENDIF

      z2 = (lat_coarse_v(i,j) - zlats) / dlat_in
      IF ( ABS(NINT(z2) - z2) < 1E-6_ireals) THEN
        j_index(i,j,5) = NINT(z2) + 1
      ELSE
        j_index(i,j,5) = INT (z2) + 1
      ENDIF

      ! weights for the interpolation: to be decomposition independent, the
      ! weights are calculated with respect to the total coarse grid

      ! first get lon-index of this grid point in the total grid
      zdifflon_tot   = lon_coarse_v(i,j) - zlons_tot
      IF (zdifflon_tot < 0.0_ireals) zdifflon_tot = zdifflon_tot + 360.0
      z1 = zdifflon_tot / dlon_in
      IF ( ABS(NINT(z1) - z1) < 1E-6_ireals) THEN
        iind_total = NINT(z1) + 1
      ELSE
        iind_total = INT (z1) + 1
      ENDIF

      ! compute weight for longitudes
      zdifflon       = lon_coarse_v(i,j) - slongitudes_in(iind_total)
      ! it could be that these grid points have the same longitude, which
      ! could lead to a numerical small negative value instead of 0.0
      IF     (ABS(zdifflon) <   1.0E-6_ireals) THEN
        ! fine grid point and coarse grid point have the same longitude
        zdifflon = 0.0_ireals
      ELSEIF (    zdifflon  <= -1.0E-6_ireals) THEN
        ! should only be the case when crossing the date line
        zdifflon = zdifflon + 360.0
      ENDIF
      x_wght (i,j,5) = zdifflon / dlon_in

      ! get lat-index of this grid point in the total grid
      z2 = (lat_coarse_v(i,j)-zlats_tot)/ dlat_in
      IF ( ABS(NINT(z2) - z2) < 1E-6_ireals) THEN
        jind_total = NINT(z2) + 1
      ELSE
        jind_total = INT (z2) + 1
      ENDIF
      ! compute weight for latitudes
      zdifflat       = lat_coarse_v(i,j) - latitudes_in(jind_total)
      IF (ABS(zdifflat) < 1.0E-6_ireals) THEN
        y_wght(i,j,5) = 0.0_ireals
      ELSE
        y_wght(i,j,5) = zdifflat / dlat_in
      END IF
    ENDDO
  ENDDO

  ! 6: Check if all weights have reasonable values
  IF (MAXVAL(x_wght(:,:,:)) .GT. 1.001_ireals .OR. &
      MINVAL(x_wght(:,:,:)) .LT. -0.001_ireals) THEN
    ierror  = 2
    yerror = 'x-weights out of bounds'
  ENDIF
  IF (MAXVAL(y_wght(:,:,:)) .GT. 1.001_ireals .OR. &
      MINVAL(y_wght(:,:,:)) .LT. -0.001_ireals) THEN
    ierror  = 3
    yerror = 'y-weights out of bounds'
  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE coor_coarse_grid_lm

!==============================================================================
!+ Module procedure for setting up the correspondence of LM grid to coarse grid
!------------------------------------------------------------------------------

SUBROUTINE coor_cm_lm (yerror, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the i- and j-indices of the coarse mesh grid point
!   which is lower left to a given LM (fine mesh) grid point and
!   the relative distances between x- (i-) and y- (j-) coarse mesh and
!   fine mesh grid points
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist
CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local variables
INTEGER (KIND=iintegers) :: &
  izstart_in, izend_in, jzstart_in, jzend_in,               &
  i, j, iz, jz, n, nzstat, izerror

CHARACTER (LEN= 80)      ::  &
  yzerror       ! error message

REAL (KIND=ireals) , ALLOCATABLE :: &
  zlongitudes(:,:), zlatitudes(:,:)

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine coor_cm_lm
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierror  = 0_iintegers
  yerror  = '         '
  izerror = 0_iintegers
  yzerror = '         '

  ! allocate longitudes and lattitudes
  ALLOCATE (zlongitudes(1:ie_in,3), zlatitudes(1:je_in,3), STAT=ierror)
  IF (ierror /= 0) THEN
    ierror  = 1
    yzerror = 'allocation error'
    RETURN
  ENDIF

  izstart_in = isubpos_coarse (my_cart_id, 1)
  jzstart_in = isubpos_coarse (my_cart_id, 2)
  izend_in   = isubpos_coarse (my_cart_id, 3)
  jzend_in   = isubpos_coarse (my_cart_id, 4)


  !  1:  coarse mesh mass grid points:
  DO iz = 1, ie_in
    zlongitudes(iz,1) = longitudes_in(izstart_in+iz-1)
  ENDDO
  DO jz = 1, je_in
    zlatitudes(jz,1) = latitudes_in(jzstart_in+jz-1)
  ENDDO

  !  2:  coarse mesh u grid point
  IF (lushift_in(1)) THEN
    DO iz = 1, ie_in
      zlongitudes(iz,2) = slongitudes_in(izstart_in+iz-1)
    ENDDO
  ELSE
    zlongitudes(:,2) = zlongitudes(:,1)
  ENDIF
  IF (lushift_in(2)) THEN
    DO jz = 1, je_in
      zlatitudes(jz,2) = slatitudes_in(jzstart_in+jz-1)
    ENDDO
  ELSE
    zlatitudes(:,2) = zlatitudes(:,1)
  ENDIF

  !  3:  coarse mesh v grid point

  IF (lvshift_in(1)) THEN
    DO iz = 1, ie_in
      zlongitudes(iz,3) = slongitudes_in(izstart_in+iz-1)
    ENDDO
  ELSE
    zlongitudes(:,3) = zlongitudes(:,1)
  ENDIF
  IF (lvshift_in(2)) THEN
    DO jz = 1, je_in
      zlatitudes(jz,3) = slatitudes_in(jzstart_in+jz-1)
    ENDDO
  ELSE
    zlatitudes(:,3) = zlatitudes(:,1)
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Generation of coarse and fine grid correspondence
!------------------------------------------------------------------------------

  ! For every LM (fine mesh) grid point the grid point of the coarse grid
  ! is determined that is just to the lower left of the LM grid point.
  ! This is done also for the staggered positions, which gives 5 possibilities:

  !  1:  mass grid points:
  DO j = 1, je2lm
    DO i = 1, ie2lm
      DO iz = 1, ie_in
        IF (lon_coarse_m(i,j) < zlongitudes(iz,1)) THEN
          i_index(i,j,1) = iz -1
          x_wght (i,j,1) = (lon_coarse_m(i,j) - zlongitudes(iz-1,1)) /    &
                            (zlongitudes(iz,1) - zlongitudes(iz-1,1))
          EXIT
        ENDIF
      ENDDO
      DO jz = 1, je_in
        IF (lat_coarse_m(i,j) < zlatitudes(jz,1)) THEN
          j_index(i,j,1) = jz -1
          y_wght (i,j,1) = (lat_coarse_m(i,j) - zlatitudes(jz-1,1)) /     &
                            (zlatitudes(jz,1) - zlatitudes(jz-1,1))
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  !  2:  fine mesh u grid point and coarse mesh u grid point
  DO j = 1, je2lm
    DO i = 1, ie2lm
      DO iz = 1, ie_in
        IF (lon_coarse_u(i,j) < zlongitudes(iz,2)) THEN
          i_index(i,j,2) = iz -1
          x_wght (i,j,2) = (lon_coarse_u(i,j) - zlongitudes(iz-1,2)) /    &
                            (zlongitudes(iz,2) - zlongitudes(iz-1,2))
          EXIT
        ENDIF
      ENDDO
      DO jz = 1, je_in
        IF (lat_coarse_u(i,j) < zlatitudes(jz,2)) THEN
          j_index(i,j,2) = jz -1
          y_wght (i,j,2) = (lat_coarse_u(i,j) - zlatitudes(jz-1,2)) /     &
                            (zlatitudes(jz,2) - zlatitudes(jz-1,2))
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  !  3:  fine mesh u grid point and coarse mesh v grid point
  DO j = 1, je2lm
    DO i = 1, ie2lm
      DO iz = 1, ie_in
        IF (lon_coarse_u(i,j) < zlongitudes(iz,3)) THEN
          i_index(i,j,3) = iz -1
          x_wght (i,j,3) = (lon_coarse_u(i,j) - zlongitudes(iz-1,3)) /    &
                            (zlongitudes(iz,3) - zlongitudes(iz-1,3))
          EXIT
        ENDIF
      ENDDO
      DO jz = 1, je_in
        IF (lat_coarse_u(i,j) < zlatitudes(jz,3)) THEN
          j_index(i,j,3) = jz -1
          y_wght (i,j,3) = (lat_coarse_u(i,j) - zlatitudes(jz-1,3)) /     &
                            (zlatitudes(jz,3) - zlatitudes(jz-1,3))
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  !  4:  fine mesh v grid point and coarse mesh v grid point
  DO j = 1, je2lm
    DO i = 1, ie2lm
      DO iz = 1, ie_in
        IF (lon_coarse_v(i,j) < zlongitudes(iz,3)) THEN
          i_index(i,j,4) = iz -1
          x_wght (i,j,4) = (lon_coarse_v(i,j) - zlongitudes(iz-1,3)) /    &
                            (zlongitudes(iz,3) - zlongitudes(iz-1,3))
          EXIT
        ENDIF
      ENDDO
      DO jz = 1, je_in
        IF (lat_coarse_v(i,j) < zlatitudes(jz,3)) THEN
          j_index(i,j,4) = jz -1
          y_wght (i,j,4) = (lat_coarse_v(i,j) - zlatitudes(jz-1,3)) /     &
                            (zlatitudes(jz,3) - zlatitudes(jz-1,3))
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  !  5:  fine mesh v grid point and coarse mesh u grid point
  DO j = 1, je2lm
    DO i = 1, ie2lm
      DO iz = 1, ie_in
        IF (lon_coarse_v(i,j) < zlongitudes(iz,2)) THEN
          i_index(i,j,5) = iz -1
          x_wght (i,j,5) = (lon_coarse_v(i,j) - zlongitudes(iz-1,2)) /    &
                            (zlongitudes(iz,2) - zlongitudes(iz-1,2))
          EXIT
        ENDIF
      ENDDO
      DO jz = 1, je_in
        IF (lat_coarse_v(i,j) < zlatitudes(jz,2)) THEN
          j_index(i,j,5) = jz -1
          y_wght (i,j,5) = (lat_coarse_v(i,j) - zlatitudes(jz-1,2)) /     &
                            (zlatitudes(jz,2) - zlatitudes(jz-1,2))
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  ! 6: Check if all weights have reasonable values
  IF (MAXVAL(x_wght(:,:,:)) .GT. 1.001_ireals .OR. &
      MINVAL(x_wght(:,:,:)) .LT. -0.001_ireals) THEN
    ierror  = 2
    yerror = 'x-weights out of bounds'
  ENDIF
  IF (MAXVAL(y_wght(:,:,:)) .GT. 1.001_ireals .OR. &
      MINVAL(y_wght(:,:,:)) .LT. -0.001_ireals) THEN
    ierror  = 3
    yerror = 'y-weights out of bounds'
  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE coor_cm_lm

!==============================================================================
!+ Module procedure for setting up the correspondence of LM grid to coarse grid
!------------------------------------------------------------------------------

SUBROUTINE compute_geo_coord

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the latitudes and longitudes of the fine LM
!   grid points
!
! Method:
!
!------------------------------------------------------------------------------

! Local variables
INTEGER (KIND=iintegers) :: &
  i, j, itot, jtot

REAL    (KIND=ireals)    :: &
  zlats, zlons

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine compute_geo_coord
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Compute the latitudes and longitudes of the fine LM grid points
!------------------------------------------------------------------------------

  !  (Arakawa C-grid) at the locations of:
  !   - mass points (arrays with "subnames" '_m')
  !   - u-   points (arrays with "subnames" '_u')
  !   - v-   points (arrays with "subnames" '_v')

  ! geographical coordinates for the mass grid points
  DO j = 1,je2lm
    jtot = isubpos(my_cart_id,2)-2*nboundlines-1
    zlats = startlat_tot + (jtot+j-1)*dlat
    DO i = 1,ie2lm
      itot = isubpos(my_cart_id,1)-2*nboundlines-1
      zlons = startlon_tot + (itot+i-1)*dlon
      IF(zlons > 180.0_ireals) THEN
        zlons = zlons - 360.0_ireals
      ENDIF
      latlm_m(i,j)   = phirot2phi(zlats, zlons, pollat, pollon, polgam)
      lonlm_m(i,j)   = rlarot2rla(zlats, zlons, pollat, pollon, polgam)
    ENDDO
  ENDDO

  ! geographical coordinates for the u-wind grid points
  DO j = 1,je2lm
    jtot = isubpos(my_cart_id,2)-2*nboundlines-1
    zlats = startlat_tot + (jtot+j-1)*dlat
    DO i = 1,ie2lm
      itot = isubpos(my_cart_id,1)-2*nboundlines-1
      zlons = startlon_tot + (itot+i-0.5_ireals)*dlon
      IF(zlons > 180.0_ireals) THEN
        zlons = zlons - 360.0_ireals
      ENDIF
      latlm_u(i,j)   = phirot2phi(zlats, zlons, pollat, pollon, polgam)
      lonlm_u(i,j)   = rlarot2rla(zlats, zlons, pollat, pollon, polgam)
    ENDDO
  ENDDO

  ! geographical coordinates for the v-wind grid points
  DO j = 1,je2lm
    jtot = isubpos(my_cart_id,2)-2*nboundlines-1
    zlats = startlat_tot + (jtot+j-0.5_ireals)*dlat
    DO i = 1,ie2lm
      itot = isubpos(my_cart_id,1)-2*nboundlines-1
      zlons = startlon_tot + (itot+i-1)*dlon
      IF(zlons > 180.0_ireals) THEN
        zlons = zlons - 360.0_ireals
      ENDIF
      latlm_v(i,j)   = phirot2phi(zlats, zlons, pollat, pollon, polgam)
      lonlm_v(i,j)   = rlarot2rla(zlats, zlons, pollat, pollon, polgam)
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
!- Section 2: Compute the latitudes and longitudes of the fine LM grid points
!             in coordinates of the coarse grid
!------------------------------------------------------------------------------

  IF (.NOT. lgme2lm) THEN
    IF ((pollat_in == 90.0_ireals) .AND. (pollon_in == 180.0_ireals)) THEN
      ! coarse grid is unrotated
      ! lat_coarse_* = latlm_* and lon_coarse_* = lonlm_*
      lat_coarse_m(:,:) = latlm_m(:,:)
      lat_coarse_u(:,:) = latlm_u(:,:)
      lat_coarse_v(:,:) = latlm_v(:,:)
      lon_coarse_m(:,:) = lonlm_m(:,:)
      lon_coarse_u(:,:) = lonlm_u(:,:)
      lon_coarse_v(:,:) = lonlm_v(:,:)
    ELSEIF ((pollat_in == pollat) .AND. (pollon_in == pollon)     &
                                  .AND. (polgam_in == polgam)) THEN
      ! coarse grid has the same rotated pole like the fine grid
      DO j = 1, je2lm
        ! computation from startlat_tot onwards to be independent of
        ! actual domain decomposition
        jtot = isubpos(my_cart_id,2)-2*nboundlines-1
        lat_coarse_m(:,j) = startlat_tot                 + (jtot+j-1) * dlat
        lat_coarse_u(:,j) = startlat_tot                 + (jtot+j-1) * dlat
        lat_coarse_v(:,j) = startlat_tot+0.5_ireals*dlat + (jtot+j-1) * dlat
      ENDDO
      DO i = 1, ie2lm
        ! computation from startlat_tot onwards to be independent of
        ! actual domain decomposition
        itot = isubpos(my_cart_id,1)-2*nboundlines-1
        lon_coarse_m(i,:) = startlon_tot                 + (itot+i-1) * dlon
        lon_coarse_u(i,:) = startlon_tot+0.5_ireals*dlon + (itot+i-1) * dlon
        lon_coarse_v(i,:) = startlon_tot                 + (itot+i-1) * dlon
      ENDDO
    ELSE
      ! coarse grid has a different rotated pole from the fine grid
      ! rotate the geographical coordinates latlm_* and lonlm_* to the
      ! rotation of the coarse system
      DO j = 1, je2lm
        DO i = 1, ie2lm
          lat_coarse_m(i,j) = phi2phirot (latlm_m(i,j), lonlm_m(i,j), pollat_in, pollon_in)
          lat_coarse_u(i,j) = phi2phirot (latlm_u(i,j), lonlm_u(i,j), pollat_in, pollon_in)
          lat_coarse_v(i,j) = phi2phirot (latlm_v(i,j), lonlm_v(i,j), pollat_in, pollon_in)
          lon_coarse_m(i,j) = rla2rlarot (latlm_m(i,j), lonlm_m(i,j),     &
                                          pollat_in, pollon_in, polgam_in)
          lon_coarse_u(i,j) = rla2rlarot (latlm_u(i,j), lonlm_u(i,j),     &
                                          pollat_in, pollon_in, polgam_in)
          lon_coarse_v(i,j) = rla2rlarot (latlm_v(i,j), lonlm_v(i,j),     &
                                          pollat_in, pollon_in, polgam_in)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of external subroutine organize_grids
!------------------------------------------------------------------------------

END SUBROUTINE compute_geo_coord

!==============================================================================
!+ Driving routine for correspondence of LM grid to GME grid 
!------------------------------------------------------------------------------

SUBROUTINE coor_gme_lm (yerror, ierror)

!------------------------------------------------------------------------------
! Description:
! Computes the indices, barycentric coordinates, angle for wind vector rotation
! for the grid points of LM/HM (mass points, u-points and v-points)
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(OUT) ::       &
  ierror                        ! error code

CHARACTER (LEN=80)      , INTENT(OUT) ::       &
  yerror                        ! error message

!------------------------------------------------------------------------------

! Local Scalars:
INTEGER (KIND=iintegers)   :: istat, izerror,             & !
  j1_minm, j1_maxm, j2_minm, j2_maxm, jd_minm, jd_maxm,   & !
  j1_minu, j1_maxu, j2_minu, j2_maxu, jd_minu, jd_maxu,   & !
  j1_minv, j1_maxv, j2_minv, j2_maxv, jd_minv, jd_maxv,   & !
  jmin(3), jmax(3)

! The following arrays are only needed as local variables in this module.
REAL    (KIND=ireals)      ::      &
  xn     (ig1sm2:ig1ep2, ig2sm2:ig2ep2, 3, nd),    & !
  erlon  (ig1sm2:ig1ep2, ig2sm2:ig2ep2, 2, nd),    & !
  erlat  (ig1sm2:ig1ep2, ig2sm2:ig2ep2, 3, nd),    & !
  rlon   (ig1s  :ig1e  , ig2s  :ig2e  ,    nd),    & !
  rlat   (ig1s  :ig1e  , ig2s  :ig2e  ,    nd),    & !
  corio  (ig1s  :ig1e  , ig2s  :ig2e  ,    nd),    & !
  eta    (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 7),        & !
  chi    (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 7),        & !
  bary   (ig1s  :ig1e  , ig2s  :ig2e  , 6),        & !
  ceneta (ig1s  :ig1e  , ig2s  :ig2e  , 6),        & !
  cenchi (ig1s  :ig1e  , ig2s  :ig2e  , 6),        & !
  rarn   (ig1sm1:ig1ep1, ig2sm1:ig2ep1   ),        & !
  area   (ig1s  :ig1ep1, ig2sm1:ig2e  , 2),        & !
  hexwgt (ig1s  :ig1e  , ig2s  :ig2e     ),        & !
  rlap   (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 7),        & !
  cpsi   (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 6),        & !
  spsi   (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 6),        & !
  grd    (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 7, 2)

REAL    (KIND=ireals)      ::      &
  erlon_glob  (igg1s-2:igg1e+2, igg2s-2:igg2e+2, 2, nd),   & !
  erlat_glob  (igg1s-2:igg1e+2, igg2s-2:igg2e+2, 3, nd),   & !
  bary_glob   (igg1s  :igg1e  , igg2s  :igg2e  , 6),       & !
  xnglob      (igg1s  :igg1e  , igg2s  :igg2e  , 3, nd)      !

REAL    (KIND=ireals)      ::      &
  ceneta_dumm (ig1s  :ig1e  , ig2s  :ig2e , 6 ),   & !
  cenchi_dumm (ig1s  :ig1e  , ig2s  :ig2e , 6 )      !
        ! These are dummy variables for the global arrays, which are
        ! in the parameter list of *pp_horint_index_p*, but not needed

INTEGER (KIND=iintegers)   ::      &
  idxrp_m (ie2lm,je2lm),  & !
  idxrp_u (ie2lm,je2lm),  & !
  idxrp_v (ie2lm,je2lm)

CHARACTER (LEN=80)         ::      &
  yzerror 

!- End of header
!------------------------------------------------------------------------------

ierror   = 0_iintegers
yerror   = '   '

izerror  = 0_iintegers
yzerror  = '   '

!------------------------------------------------------------------------------
! Initializations: ALLOCATE and set fields to 0
!------------------------------------------------------------------------------

! Allocate the global variables needed for further processing
ALLOCATE( grd_glob    (igg1s-1:igg1e+1, igg2s-1:igg2e+1, 7, 2) )
ALLOCATE( eta_glob    (igg1s-1:igg1e+1, igg2s-1:igg2e+1, 7) )
ALLOCATE( chi_glob    (igg1s-1:igg1e+1, igg2s-1:igg2e+1, 7) )
ALLOCATE( cpsi_glob   (igg1s-1:igg1e+1, igg2s-1:igg2e+1, 6) )
ALLOCATE( spsi_glob   (igg1s-1:igg1e+1, igg2s-1:igg2e+1, 6) )

grd_glob (:,:,:,:) = 0.0_ireals
eta_glob (:,:,:)   = 0.0_ireals
chi_glob (:,:,:)   = 0.0_ireals
cpsi_glob(:,:,:)   = 0.0_ireals
spsi_glob(:,:,:)   = 0.0_ireals

!------------------------------------------------------------------------------
! Allocation of LM arrays for pp_horint_index_p for the mass points '*_m',
! the u-points '*_u' and the v-points '*_v' of the Arakawa C-grid
!------------------------------------------------------------------------------

ALLOCATE(index_m(ie2lm,je2lm,4), index_u(ie2lm,je2lm,4),       &
         index_v(ie2lm,je2lm,4), STAT=istat)

ALLOCATE(baryll_m(ie2lm,je2lm,2), baryll_u(ie2lm,je2lm,2),     &
         baryll_v(ie2lm,je2lm,2), STAT=istat)

ALLOCATE(rotang_m(ie2lm,je2lm,2), rotang_u(ie2lm,je2lm,2),     &
         rotang_v(ie2lm,je2lm,2), STAT=istat)

!------------------------------------------------------------------------------
!     Generation of the horizontal (triangular) grid mesh, computation of 
!     the gradient and Laplacian operator
!------------------------------------------------------------------------------

CALL pr_genhor (xn    , erlon , erlat , rlon  , rlat  ,           &
                corio , eta   , chi   , bary  , ceneta,           &
                cenchi, cpsi  , spsi  , grd   ,                   &
                rarn  , area  , hexwgt, rlap  , xnglob, izerror, yzerror)

IF (izerror /= 0_iintegers) THEN
  ierror = 1
  yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
  RETURN
ENDIF

!------------------------------------------------------------------------------
! Gather the global fileds needed for the interpolation
!------------------------------------------------------------------------------

CALL gather_gme_field(erlon,      ig1s , ig1e , ig2s , ig2e ,            &
                      erlon_glob, igg1s, igg1e, igg2s, igg2e, 2, 2*nd,   &
                      nproc1, nproc2, my_cart_id, my_num1, my_num2,      &
                      ilim1, ilim2, imp_reals, max_gme_core, icomm_cart, &
                      izerror)
IF (izerror /= 0) THEN
  ierror = 2
  yerror = 'gather_gme_field failed in coor_gme_lm'
  RETURN
ENDIF

CALL gather_gme_field(erlat,      ig1s , ig1e , ig2s , ig2e ,            &
                      erlat_glob, igg1s, igg1e, igg2s, igg2e, 2, 3*nd,   &
                      nproc1, nproc2, my_cart_id, my_num1, my_num2,      &
                      ilim1, ilim2, imp_reals, max_gme_core, icomm_cart, &
                      izerror)
IF (izerror /= 0) THEN
  ierror = 3
  yerror = 'gather_gme_field failed in coor_gme_lm'
  RETURN
ENDIF

CALL gather_gme_field(eta  ,      ig1s , ig1e , ig2s , ig2e ,            &
                      eta_glob,   igg1s, igg1e, igg2s, igg2e, 1, 7,      &
                      nproc1, nproc2, my_cart_id, my_num1, my_num2,      &
                      ilim1, ilim2, imp_reals, max_gme_core, icomm_cart, &
                      izerror)
IF (izerror /= 0) THEN
  ierror = 4
  yerror = 'gather_gme_field failed in coor_gme_lm'
  RETURN
ENDIF

CALL gather_gme_field(chi  ,      ig1s , ig1e , ig2s , ig2e ,            &
                      chi_glob,   igg1s, igg1e, igg2s, igg2e, 1, 7,      &
                      nproc1, nproc2, my_cart_id, my_num1, my_num2,      &
                      ilim1, ilim2, imp_reals, max_gme_core, icomm_cart, &
                      izerror)
IF (izerror /= 0) THEN
  ierror = 5
  yerror = 'gather_gme_field failed in coor_gme_lm'
  RETURN
ENDIF

CALL gather_gme_field(bary ,      ig1s , ig1e , ig2s , ig2e ,            &
                      bary_glob,  igg1s, igg1e, igg2s, igg2e, 0, 6,      &
                      nproc1, nproc2, my_cart_id, my_num1, my_num2,      &
                      ilim1, ilim2, imp_reals, max_gme_core, icomm_cart, &
                      izerror)
IF (izerror /= 0) THEN
  ierror = 6
  yerror = 'gather_gme_field failed in coor_gme_lm'
  RETURN
ENDIF

!------------------------------------------------------------------------------
! Computes the indices, barycentric coordinates, angle for wind vector rotation
! for the grid points of LM/HM (mass points, u-points and v-points)
!------------------------------------------------------------------------------

CALL  pp_horint_index_p (lonlm_m    , latlm_m    , ie2lm*je2lm,         &
                         xnglob     , erlon_glob , erlat_glob ,         &
                         eta_glob   , chi_glob   , ceneta_dumm,         &
                         cenchi_dumm, bary_glob  ,                      &
                         igg1s      , igg1e      , igg2s      , igg2e,  &
                         igg1s      , igg1e      , igg2s      , igg2e,  &
                         index_m    , baryll_m   , rotang_m,            &
                         idxrp_m    , rpoints_m)

CALL  pp_horint_index_p (lonlm_u    , latlm_u    , ie2lm*je2lm,         &
                         xnglob     , erlon_glob , erlat_glob ,         &
                         eta_glob   , chi_glob   , ceneta_dumm,         &
                         cenchi_dumm, bary_glob  ,                      &
                         igg1s      , igg1e      , igg2s      , igg2e,  &
                         igg1s      , igg1e      , igg2s      , igg2e,  &
                         index_u    , baryll_u   , rotang_u,            &
                         idxrp_u    , rpoints_u)

CALL  pp_horint_index_p (lonlm_v    , latlm_v    , ie2lm*je2lm,         &
                         xnglob     , erlon_glob , erlat_glob ,         &
                         eta_glob   , chi_glob   , ceneta_dumm,         &
                         cenchi_dumm, bary_glob  ,                      &
                         igg1s      , igg1e      , igg2s      , igg2e,  &
                         igg1s      , igg1e      , igg2s      , igg2e,  &
                         index_v    , baryll_v   , rotang_v,            &
                         idxrp_v    , rpoints_v)

!------------------------------------------------------------------------------
! Determine global min and max for j1-, j2- and jd
!------------------------------------------------------------------------------

  j1_minm  = MINVAL (index_m(:,:,1))
  j1_maxm  = MAXVAL (index_m(:,:,1))
  j2_minm  = MINVAL (index_m(:,:,2))
  j2_maxm  = MAXVAL (index_m(:,:,2))
  jd_minm  = MINVAL (index_m(:,:,3))
  jd_maxm  = MAXVAL (index_m(:,:,3))

  j1_minu  = MINVAL (index_u(:,:,1))
  j1_maxu  = MAXVAL (index_u(:,:,1))
  j2_minu  = MINVAL (index_u(:,:,2))
  j2_maxu  = MAXVAL (index_u(:,:,2))
  jd_minu  = MINVAL (index_u(:,:,3))
  jd_maxu  = MAXVAL (index_u(:,:,3))

  j1_minv  = MINVAL (index_v(:,:,1))
  j1_maxv  = MAXVAL (index_v(:,:,1))
  j2_minv  = MINVAL (index_v(:,:,2))
  j2_maxv  = MAXVAL (index_v(:,:,2))
  jd_minv  = MINVAL (index_v(:,:,3))
  jd_maxv  = MAXVAL (index_v(:,:,3))

  jmin(1)  = MIN (j1_minm, j1_minu, j1_minv)    ! j1_min
  jmax(1)  = MAX (j1_maxm, j1_maxu, j1_maxv)    ! j1_max
  jmin(2)  = MIN (j2_minm, j2_minu, j2_minv)    ! j2_min
  jmax(2)  = MAX (j2_maxm, j2_maxu, j2_maxv)    ! j2_max
  jmin(3)  = MIN (jd_minm, jd_minu, jd_minv)    ! jd_min
  jmax(3)  = MAX (jd_maxm, jd_maxu, jd_maxv)    ! jd_max

  IF (num_compute > 1) THEN
    CALL global_values (jmax, 3, 'MAX', imp_integers, icomm_cart, -1,      &
                        yzerror, izerror)
    CALL global_values (jmin, 3, 'MIN', imp_integers, icomm_cart, -1,      &
                        yzerror, izerror)
  ENDIF

  j1_min = jmin(1)
  j1_max = jmax(1)
  j2_min = jmin(2)
  j2_max = jmax(2)
  jd_min = jmin(3)
  jd_max = jmax(3)

!------------------------------------------------------------------------------
! Clean up all GME fields which are no longer needed
! The only fields which are needed any more are *grd*, *eta*, *chi*, *cpsi*
! and  *spsi* (globally)
!------------------------------------------------------------------------------

CALL gather_gme_field(grd,        ig1s , ig1e , ig2s , ig2e ,            &
                      grd_glob,   igg1s, igg1e, igg2s, igg2e, 1, 14,     &
                      nproc1, nproc2, my_cart_id, my_num1, my_num2,      &
                      ilim1, ilim2, imp_reals, max_gme_core, icomm_cart, &
                      izerror)
IF (izerror /= 0) THEN
  ierror = 7
  yerror = 'gather_gme_field failed in coor_gme_lm'
  RETURN
ENDIF

CALL gather_gme_field(cpsi ,      ig1s , ig1e , ig2s , ig2e ,            &
                      cpsi_glob,  igg1s, igg1e, igg2s, igg2e, 1, 6,      &
                      nproc1, nproc2, my_cart_id, my_num1, my_num2,      &
                      ilim1, ilim2, imp_reals, max_gme_core, icomm_cart, &
                      izerror)
IF (izerror /= 0) THEN
  ierror = 8
  yerror = 'gather_gme_field failed in coor_gme_lm'
  RETURN
ENDIF

CALL gather_gme_field(spsi ,      ig1s , ig1e , ig2s , ig2e ,            &
                      spsi_glob,  igg1s, igg1e, igg2s, igg2e, 1, 6,      &
                      nproc1, nproc2, my_cart_id, my_num1, my_num2,      &
                      ilim1, ilim2, imp_reals, max_gme_core, icomm_cart, &
                      izerror)
IF (izerror /= 0) THEN
  ierror = 9
  yerror = 'gather_gme_field failed in coor_gme_lm'
  RETURN
ENDIF

!------------------------------------------------------------------------------

END SUBROUTINE coor_gme_lm

!==============================================================================
!==============================================================================
!+ Module procedure in src_gme_lm_grid for generating GME-grid
!------------------------------------------------------------------------------

SUBROUTINE pr_genhor (xn    , erlon , erlat , rlon  , rlat  ,           &
                      corio , eta   , chi   , bary  , ceneta,           &
                      cenchi, cpsi  , spsi  , grd   ,                   &
                      rarn  , area  , hexwgt, rlap  , xnglob,           &
                      ierror, yerror)

!------------------------------------------------------------------------------
!
! Description:
!   *pr_genhor* generates the horizontal (triangular) grid mesh, 
!   calculates the gradient and Laplacian operator
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(OUT) ::       &
  ierror                        ! error code

CHARACTER (LEN=80)      , INTENT(OUT) ::       &
  yerror                        ! error message

!------------------------------------------------------------------------------

REAL    (KIND=ireals), INTENT(INOUT)       ::      &
  xn     (ig1sm2:ig1ep2, ig2sm2:ig2ep2, 3, nd),    & !
  erlon  (ig1sm2:ig1ep2, ig2sm2:ig2ep2, 2, nd),    & !
  erlat  (ig1sm2:ig1ep2, ig2sm2:ig2ep2, 3, nd),    & !
  rlon   (ig1s  :ig1e  , ig2s  :ig2e  ,    nd),    & !
  rlat   (ig1s  :ig1e  , ig2s  :ig2e  ,    nd),    & !
  corio  (ig1s  :ig1e  , ig2s  :ig2e  ,    nd),    & !
  eta    (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 7),        & !
  chi    (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 7),        & !
  bary   (ig1s  :ig1e  , ig2s  :ig2e  , 6),        & !
  ceneta (ig1s  :ig1e  , ig2s  :ig2e  , 6),        & !
  cenchi (ig1s  :ig1e  , ig2s  :ig2e  , 6),        & !
  rarn   (ig1sm1:ig1ep1, ig2sm1:ig2ep1   ),        & !
  area   (ig1s  :ig1ep1, ig2sm1:ig2e  , 2),        & !
  hexwgt (ig1s  :ig1e  , ig2s  :ig2e     ),        & !
  rlap   (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 7),        & !
  cpsi   (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 6),        & !
  spsi   (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 6),        & !
  grd    (ig1sm1:ig1ep1, ig2sm1:ig2ep1, 7, 2)


REAL    (KIND=ireals), INTENT(INOUT)       ::      &
  xnglob (igg1s:igg1e,igg2s:igg2e,3,nd)

! Local variables used
  INTEGER  (KIND=iintegers) ::    &
    mierr       ! Error flag; 0 if no error occured

  LOGICAL                   ::   lzdebug

  CHARACTER (LEN=30) yzroutine
  CHARACTER (LEN=80) yzerrmsg 

!------------------------------------------------------------------------------

  ierror = 0_iintegers
  yerror = '   '

!------------------------------------------------------------------------------
!
!     PART I  Calculate the horizontal grid
!
!     1. Factorize ni into 3**ni3 * 2**ni2 with ni3 0 or 1 and ni2 > 1,
!        if possible, otherwise abort the program
!
      mierr     = 0
      yzerrmsg  = '   '
      yzroutine = 'pr_genhor'
      lzdebug   = .FALSE.
!
      call factorni (ni_gme, lzdebug, ni2, ni3, mierr)
!
      IF (mierr.NE.0) THEN
        ierror = 1
        yerror = 'Error in *pr_genhor*: ni cannot be factorized correctly'
        RETURN
      ENDIF
      
!
!------------------------------------------------------------------------------
!
!     2. Calculate the arrays xn, rlon, rlat, and corio for all 10
!        diamonds
!
      CALL glo_coor (xn    , rlon  , rlat  , corio , xnglob, &
                     ig1sm2, ig1ep2, ig2sm2, ig2ep2, nd,     &
                     ig1sm1, ig1ep1, ig2sm1, ig2ep1,         &
                     ig1s  , ig1e  , ig2s  , ig2e  ,         &
                     ni2   , ni3   , lzdebug, mierr)
!
      IF (mierr.NE.0) THEN
        ierror = 2
        yerror = 'Error in *glo_coor*: STOP!'
        RETURN
      ENDIF
!
!------------------------------------------------------------------------------
!
!     3. Calculate the minimum mesh width and minimum height of
!        the grid
!
      CALL dxdhmin (xn   , ig1sm2, ig1ep2, ig2sm2, ig2ep2, nd,  &
                    ig1s , ig1e  , ig2s  , ig2e  , lzdebug,     &
                    dxmin, dhmin , mierr)
!
      IF (mierr.NE.0) THEN
        ierror = 3
        yerror = 'Error in *dxdhmin*: STOP!'
        RETURN
      ENDIF
!
!------------------------------------------------------------------------------
!
!     4. Calculate the arrays erlon, erlat, eta, chi, cpsi, spsi, bary,
!        ceneta, cenchi which describe the local coordinate system
!
      CALL loc_coor (xn     , ig1sm2 , ig1ep2, ig2sm2, ig2ep2, nd,      &
                     lzdebug,                                           &
                     erlon  , erlat  , eta   , chi   , cpsi  , spsi  ,  &
                     bary   , ceneta , cenchi,                          &
                     ispoke , ispokes, i1mrp , i2mrp ,                  &
                     ig1sm1 , ig1ep1 , ig2sm1, ig2ep1,                  &
                     ig1s   , ig1e   , ig2s  , ig2e  , mierr)
!
      IF (mierr.NE.0) THEN
        ierror = 4
        yerror = 'Error in *loc_coor*: STOP!'
        RETURN
      ENDIF
!
!------------------------------------------------------------------------------
!
!     5. Calculate the areas of the triangles for diamond #1
!
      CALL tri_area (xn   , ig1sm2, ig1ep2, ig2sm2, ig2ep2, nd, lzdebug, &
                     area , ig1s  , ig1ep1, ig2sm1, ig2e  , mierr) 
!
      IF (mierr.NE.0) THEN
        ierror = 5
        yerror = 'Error in *tri_area*: STOP!'
        RETURN
      ENDIF
!
!------------------------------------------------------------------------------
!
!     6. Calculate the areas of the hexagons for extended diamond #1 and
!        the weights of the gridpoints for diagnostic evaluations
!
      CALL hex_area (area  ,                                            &
                     ig1s  , ig1e  , ig2s  , ig2e  ,                    &
                     ig1sm1, ig1ep1, ig2sm1, ig2ep1, nd, lzdebug,       &
                     rarn  , hexwgt, mierr) 
!
      IF (mierr.NE.0) THEN
        ierror = 6
        yerror = 'Error in *hex_area*: STOP!'
        RETURN
      ENDIF
!
!------------------------------------------------------------------------------
!
!     7. Calculate the gradient  operator for the extended diamond #1,
!        calculate the Laplacian operator for the extended diamond #1
!
      CALL gen_grd (eta   , chi   ,                                    &
                    ig1sm1, ig1ep1, ig2sm1, ig2ep1,                    &
                    ig1s  , ig1e  , ig2s  , ig2e  ,                    &
                    lzdebug,grd   , rlap  , mierr)
!
      IF (mierr.NE.0) THEN
        ierror = 7
        yerror = 'Error in *gen_grd*: STOP!'
        RETURN
      ENDIF
!
!------------------------------------------------------------------------------

END SUBROUTINE pr_genhor

!==============================================================================
!==============================================================================
!+ Defines all subroutines for the generic routine gather_field
!------------------------------------------------------------------------------

SUBROUTINE glo_coor (pxn    , prlon  , prlat  , pcorio , pxnglob,       &
                     kig1sm2, kig1ep2, kig2sm2, kig2ep2, knd,           &
                     kig1sm1, kig1ep1, kig2sm1, kig2ep1,                &
                     kig1s  , kig1e  , kig2s  , kig2e  ,                &
                     kni2   , kni3   , lzdebug, kierr)

!------------------------------------------------------------------------------
!
! Description:
!   *glo_coor* calculates the global arrays which define the
!   triangular grid
!
! Method:
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!     Input/Output
!     pxn    (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 3, knd) REAL:  Cartesian 
!            (x,y,z) coordinates of the location vector of the grid
!            points (nodes) on the unit-sphere; phys. dim. ( - )
!            Attention: pxn is defined for the extended domain, i.e. 
!            with 2 rows and columns around the original diamond array
!            size.
!     prlon  (kig1s  :kig1e  , kig2s  :kig2e  ,    knd):  geographical
!            longitude of the gridpoints; phys. dim. (radians)
!     prlat  (kig1s  :kig1e  , kig2s  :kig2e  ,    knd):  geographical
!            latitude of the gridpoints; phys. dim. (radians)
!     pcorio (kig1s  :kig1e  , kig2s  :kig2e  ,    knd):  Coriolis 
!            parameter at the gridpoints, only the vertical component;
!            phys. dim. (1/s)
!  
!     kig1s   INTEGER  first dimension of arrays, start index (kig1s = 0)
!     kig1sm1 INTEGER  kig1sm1 = kig1s - 1
!     kig1sm2 INTEGER  kig1sm2 = kig1s - 2
!     kig1e   INTEGER  first dimension of arrays, end index (kig1e = ni)
!     kig1ep1 INTEGER  kig1ep1 = kig1e + 1
!     kig1ep2 INTEGER  kig1ep2 = kig1e + 2
!     kig2s   INTEGER  sec. dimension of arrays, start index (kig2s = 1)
!     kig2sm1 INTEGER  kig2sm1 = kig2s - 1
!     kig2sm2 INTEGER  kig2sm2 = kig2s - 2
!     kig2e   INTEGER  sec. dimension of arrays, end index (kig2e = ni+1)
!     kig2ep1 INTEGER  kig2ep1 = kig2e + 1
!     kig2ep2 INTEGER  kig2ep2 = kig2e + 2
!     knd     INTEGER  number of diamonds (knd = 10)
!     kni2    INTEGER  exponent of factor "2" in the factorization of ni,
!                      i.e. the number of bisections to perform
!     kni3    INTEGER  exponent of factor "3" in the factorization of ni,
!                      i.e. the number of trisections to perform (0 or1)
!     lzdebug LOGICAL  debug flag, if .true. print debug information
!
!------------------------------------------------------------------------------
!
!     Output
!     kierr   INTEGER  error flag, set to 0, if no error occured
!
!
!     Global header files
!      INCLUDE "comconst.h"
!     INCLUDE "globallimits.h"
!     INCLUDE "horgrid.h"  ! RJ: Only needed for xnglob
!                          ! RJ: *xnglob* should better be moved to the
!                          !     parameter list - I know ...
!
! Dummy arguments
  INTEGER  (KIND=iintegers)   ::      &
            kig1sm2, kig1ep2, kig2sm2, kig2ep2, knd , kierr,      &
            kig1sm1, kig1ep1, kig2sm1, kig2ep1, kni2, kni3 ,      &
            kig1s  , kig1e  , kig2s  , kig2e  
!
  REAL     (KIND=ireals)      ::      &
            pxn    (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 3, knd),    &
            prlon  (kig1s  :kig1e  , kig2s  :kig2e  ,    knd),    &
            prlat  (kig1s  :kig1e  , kig2s  :kig2e  ,    knd),    &
            pcorio (kig1s  :kig1e  , kig2s  :kig2e  ,    knd)

  REAL     (KIND=ireals)      ::      &
            pxnglob (igg1s:igg1e,igg2s:igg2e,3,nd)

  LOGICAL  lzdebug
!
!------------------------------------------------------------------------------
!
! Local arrays and variables
  REAL     (KIND=ireals)      ::      &
    zw   , & ! spherical angle in an icosahedron subtended by two vertices.
    zcosw, & ! cosine(zw)
    zsinw, & ! sine  (zw)
    zsgn , & ! zsgn is a hemisphere factor.
             !           zsgn =  1.0 is north  (diamonds 1- 5)
             !           zsgn = -1.0 is south  (diamonds 6-10)
    zrlon, & ! longitude of diamond vertices
    zgamma   ! fraction of great circle angle

  INTEGER  (KIND=iintegers)   ::      &
    mcosv(knd) ! meridian angle locations of the 10
               ! non-polar vertices in units of Pi/5

  INTEGER  (KIND=iintegers)   ::      &
    ml,    & ! recursive index interval
    ml2,   & ! recursive bisected index interval
    ml3,   & ! trisected index interval
    mi1,   & ! recursive row index of new node
    mi2,   & ! recursive column index of new node
    mm ,   & ! recursive number of subdivisions
    mni      ! number of subdivisions;
             ! mni=2**kni2 * 3**kni3

  INTEGER  (KIND=iintegers)   ::      &
    j1, j2, jd, jb  ! Loop indices

! Subroutines/functions used
!     EXTERNAL  gcpt,      ! computes the nodal location along a
!                            great circle arc
!    1          tricntr,   ! computes the coordinates of the 
!                            icosahedral triangle centers,; used only
!                            if a trisection is required (kni3=1)
!    2          xd_p,      ! extends the array pxn by two rows/colums
!    3          printrarr  ! prints a real 2-dimensional array
!
!
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE glo_coor
!------------------------------------------------------------------------------

! 1. Calculate the Cartesian coordinates of the gridpoints of the
!    icosahedral grid on the unit sphere.  The grid resolution
!    corresponds to a subdivision of the edges of the original
!    icosahedral triangles into mni equal parts.
!
!
! Compute angles associated with the icosahedron.
  zw      = 2.0_ireals*ACOS(1.0_ireals/(2.0_ireals*SIN(Pid5)))
  zcosw   = COS(zw)
  zsinw   = SIN(zw)
  mni     = 2**kni2 * 3**kni3

! Compute the local array mcosv, i.e. the meridian angle locations
  DO jd = 1,knd
    IF (MOD(jd,2) .eq. 1) THEN
      mcosv((jd+1)/2) = -1 + (jd - 1) - knd*((jd - 1)/7)
    ELSE   
      mcosv(jd/2+5)   = -1 + (jd - 1) - knd*((jd - 1)/7)
    ENDIF    
  ENDDO 
!
!------------------------------------------------------------------------------
!
 ! Loop over the ten diamonds computing diamond vertices (x,y,z)
 ! coordinates and then iteratively bisecting them kni2 times.
 ! First a trisection is performed, if required (kni3=1).
 DO jd = 1,knd     ! Loop over the diamonds

   ! Toggle the hemisphere
   IF (jd .GE. 6) THEN
     zsgn = -1.0_ireals    ! southern
   ELSE
     zsgn =  1.0_ireals    ! northern
   ENDIF 

   ! Compute the meridian angle for each diamond "home" vertex.
   zrlon = mcosv(jd)*Pid5

   ! Every diamond has one vertex at a pole (N or S).
   ! Label this point (0,1,,) in each diamond, and
   ! initialize it to have the (x,y,z) coordinates of
   ! the pole point on the unit sphere.
   pxnglob(  0,    1, 1, jd) =  0.0_ireals
   pxnglob(  0,    1, 2, jd) =  0.0_ireals
   pxnglob(  0,    1, 3, jd) =  zsgn

   ! Now initialize the (x,y,z) coordinates of the "home" vertex,
   ! which defines which diamond we are talking about, and label
   ! this point (mni,1,,).
   pxnglob(mni,    1, 1, jd) =  zsinw*COS(zrlon)
   pxnglob(mni,    1, 2, jd) =  zsinw*SIN(zrlon)
   pxnglob(mni,    1, 3, jd) =  zcosw*zsgn

   ! Next initialize the (x,y,z) coordinates for the corner of the
   ! diamond on the same latitude as the (mni,1,,) vertex, which
   ! is (0,mni+1,,)
   pxnglob(  0,mni+1, 1, jd) =  zsinw*COS(zrlon + 2*Pid5)
   pxnglob(  0,mni+1, 2, jd) =  zsinw*SIN(zrlon + 2*Pid5)
   pxnglob(  0,mni+1, 3, jd) =  zcosw*zsgn

   ! Initialize the last diamond vertex, which is located
   ! in the opposite hemisphere as (mni,mni+1,,)
   pxnglob(mni,mni+1, 1, jd) =  zsinw*COS(zrlon + Pid5)
   pxnglob(mni,mni+1, 2, jd) =  zsinw*SIN(zrlon + Pid5)
   pxnglob(mni,mni+1, 3, jd) = -zcosw*zsgn

!------------------------------------------------------------------------------

   ! First a trisection is performed, if required (kni3=1).
   IF (kni3 .EQ. 1) THEN
     ml3 = mni/3

     ! Trisect the rows of the diamond.
     DO j1 = 1,2
       DO j2 = 1,2
         mi1    = (j1-1)*mni
         mi2    = j2*ml3 + 1
         zgamma = REAL(j2,ireals)/3.0_ireals
         CALL gcpt (pxnglob, igg1s, igg1e, igg2s, igg2e, knd,     &
                    jd , zgamma,                                 &
                    mi1, 1,                                      &
                    mi1, mni+1,                                  &
                    mi1, mi2)
       ENDDO
     ENDDO

     ! Trisect the columns of the diamond.
     DO j1 = 1,2
       DO j2 = 1,2
         mi1    = j2*ml3
         mi2    = (j1-1)*mni + 1
         zgamma = REAL(j2,ireals)/3.0_ireals
         CALL gcpt (pxnglob, igg1s, igg1e, igg2s, igg2e, knd,     &
                    jd , zgamma,                                 &
                    0  , mi2,                                    &
                    mni, mi2,                                    &
                    mi1, mi2)
       ENDDO
     ENDDO

     ! Trisect the diagonal of the diamond.
     DO j2 = 1,2
       mi1 = mni - j2*ml3
       mi2 =   1 + j2*ml3
       zgamma = REAL(j2,ireals)/3.0_ireals
       CALL gcpt (pxnglob, igg1s, igg1e, igg2s, igg2e, knd,       &
                  jd , zgamma,                                   &
                  mni, 1,                                        &
                  0  , mni+1,                                    &
                  mi1, mi2)
     ENDDO

     ! Compute coordinates of icosahedral triangle centers.
     CALL tricntr (pxnglob, igg1s, igg1e, igg2s, igg2e, knd, jd, mni)

   ENDIF     ! End of trisection

!------------------------------------------------------------------------------

   ! Find the coordinates of the triangle nodes by iteratively
   ! bisecting the diamond intervals.
   DO jb = 0, kni2-1
     mm  = (3**kni3)*(2**jb)
     ml  = mni/mm
     ml2 = ml/2

     ! Compute the rows of the diamond.
     DO j1 = 1,mm+1
       DO j2 = 1,mm
         mi1 = (j1-1)*ml
         mi2 = (j2-1)*ml + ml2 + 1
         CALL gcpt (pxnglob, igg1s, igg1e, igg2s, igg2e, knd,     &
                    jd , 0.5_ireals,                             &
                    mi1, mi2 - ml2,                              &
                    mi1, mi2 + ml2,                              &
                    mi1, mi2)
       ENDDO
     ENDDO

     ! Compute the columns of diamond.
     DO j1 = 1,mm+1
       DO j2 = 1,mm
         mi1 = (j2-1)*ml + ml2 
         mi2 = (j1-1)*ml + 1
         CALL gcpt (pxnglob, igg1s, igg1e, igg2s, igg2e, knd,     &
                    jd , 0.5_ireals,                             &
                    mi1 - ml2,mi2,                               &
                    mi1 + ml2,mi2,                               &
                    mi1      ,mi2)
       ENDDO
     ENDDO

     ! Compute the diagonals of the diamond.
     DO j1 = 1,mm
       DO j2 = 1,mm
         mi1 = (j1-1)*ml + ml2 
         mi2 = (j2-1)*ml + ml2 + 1
         CALL gcpt (pxnglob, igg1s, igg1e, igg2s, igg2e, knd,     &
                    jd , 0.5_ireals,                             &
                    mi1 - ml2, mi2 + ml2,                        &
                    mi1 + ml2, mi2 - ml2,                        &
                    mi1      , mi2)
       ENDDO
     ENDDO
!
   ENDDO       ! end loop over bisections

!------------------------------------------------------------------------------

   ! Set pxnglob to 0 if it is less than 2.5 e-14 to avoid round-off errors
   DO j2 = igg2s,igg2e
     DO j1 = igg1s,igg1e
       IF (ABS(pxnglob(j1,j2,1,jd)) .LT. 2.5e-14_ireals)         &
               pxnglob(j1,j2,1,jd) = 0.0_ireals
       IF (ABS(pxnglob(j1,j2,2,jd)) .LT. 2.5e-14_ireals)         &
               pxnglob(j1,j2,2,jd) = 0.0_ireals
       IF (ABS(pxnglob(j1,j2,3,jd)) .LT. 2.5e-14_ireals)         &
               pxnglob(j1,j2,3,jd) = 0.0_ireals
     ENDDO
   ENDDO
!
  ENDDO         ! end loop over diamonds

!------------------------------------------------------------------------------
!
  ! Extend the nodal array by two rows/colums around
  ! This is done by the subroutine *xd* for all diamonds and the
  ! three Cartesian (x,y,z) coordinates simultaneously.

  pxn(kig1s:kig1e,kig2s:kig2e,:,:)= pxnglob(kig1s:kig1e,kig2s:kig2e,:,:)

  CALL xd_p( pxn, kig1sm2, kig1ep2, kig2sm2, kig2ep2, 1, 3, 1, 2 ,   &
             my_cart_id, num_compute, icomm_cart, imp_reals)

!------------------------------------------------------------------------------

  ! Calculate the longitude "prlon" and the latitude "prlat";
  ! only for the core of the diamonds, not the extended ones.
  DO jd = 1,knd     ! Loop over the diamonds
    DO j2 = kig2s, kig2e
      DO j1 = kig1s, kig1e
        prlon (j1,j2,jd) = ATAN2   (pxn(j1,j2,2,jd),                  &
                                    pxn(j1,j2,1,jd) + 1.e-20_ireals)
        prlat (j1,j2,jd) = ASIN    (pxn(j1,j2,3,jd))
        pcorio(j1,j2,jd) = 2.*Omcor*pxn(j1,j2,3,jd)
      ENDDO
    ENDDO
  ENDDO         ! end loop over diamonds

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE glo_coor

!==============================================================================
!==============================================================================
!+ Computes values for GME triangular grid
!------------------------------------------------------------------------------

SUBROUTINE dxdhmin (pxn    ,                                            &
                    kig1sm2, kig1ep2, kig2sm2, kig2ep2, knd,            &
                    ki1sc  , ki1ec  , ki2sc  , ki2ec  , lzdebug,        &
                    pdxmin , pdhmin , kierr)  

!------------------------------------------------------------------------------
!
! Description:
!   *dxdhmin* calculates the actual minimum mesh size "pdxmin"
!   and minimum height of the spherical triangles "pdhmin"
!   for the triangular grid. "pdhmin" determines the maximum
!   timestep allowed in connection with the maximum wind speed
!   "vmax" following dt < dhmin/(2*vmax).
!   Additionally, the local variables "zdxmax" and "zdhmax"
!   contain the maximum mesh size and triangle height.
!
! Method:
!
!------------------------------------------------------------------------------
!
!     Input
!     pxn    (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 3, knd)  REAL Cartesian 
!            (x,y,z) coordinates of the location vector of the grid
!            points (nodes) on the unit-sphere; phys. dim. ( - )
!            Attention: pxn is defined for the extended domain, i.e. 
!            with 2 rows and columns around the original diamond array
!            size.
!
!     kig1sm2 INTEGER  kig1sm2 = kig1s - 2
!                      with kig1s first dimension of arrays,
!                      start index (kig1s = 0)
!     kig1ep2 INTEGER  kig1ep2 = kig1e + 2
!                      with kig1e first dimension of arrays,
!                      end index (kig1e = ni)
!     kig2sm2 INTEGER  kig2sm2 = kig2s - 2
!                      with kig2s second dimension of arrays,
!                      start index (kig2s = 1)
!     kig2ep2 INTEGER  kig2ep2 = kig2e + 2
!                      with kig2e second dimension of arrays,
!     knd     INTEGER  number of diamonds (knd = 10)
!     ki1sc   INTEGER  start index of computation for first index
!     ki1ec   INTEGER  end   index of computation for first index
!     ki2sc   INTEGER  start index of computation for second index
!     ki2ec   INTEGER  end   index of computation for second index
!     lzdebug LOGICAL  debug flag, if .true. print debug information
!
!------------------------------------------------------------------------------
!
!     Output
!     pdxmin  REAL     minimum mesh size of triangular grid
!     pdhmin  REAL     minimum height of spherical triangle
!     kierr   INTEGER  error flag, set to 0, if no error occured
!
!------------------------------------------------------------------------------

!     Dummy arguments
  INTEGER  (KIND=iintegers)   ::                                  &
            kig1sm2, kig1ep2, kig2sm2, kig2ep2, knd, kierr,       &
            ki1sc  , ki1ec  , ki2sc  , ki2ec 

  REAL     (KIND=ireals)      ::                                  &
            pxn    (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 3, knd)
!
  LOGICAL  lzdebug
!
  REAL     (KIND=ireals)      ::     pdxmin, pdhmin
!
!     Local variables
  REAL     (KIND=ireals)      ::                                  &
            zsp1,  & ! scalar product of the location vectors of the
                     ! home node and neighbouring node # 1
            zsp2,  & ! same for # 2
            zsp3, zsp4, zsp5, zsp6, &
            zsp12, & ! scalar product of the location vectors of 
                     ! neigbouring nodes 1 and 2
            zsp23, & ! same for # 2 and # 3
            zsp34, zsp45, zsp56, zsp61, &
            zdxmax,& ! maximum mesh width of grid
            zdhmax,& ! maximum height of triangle
            zdxdhs(2)! for MPI_Allreduce
!
  INTEGER  (KIND=iintegers)   ::     j1, j2, izerror ! loop indices, etc.
!
  CHARACTER (LEN=80)          ::     yzerrmsg

!------------------------------------------------------------------------------

!     Presetting variables
      pdxmin   = 1.0E20_ireals
      pdhmin   = 1.0E20_ireals
      zdxmax   = 0.0_ireals
      zdhmax   = 0.0_ireals
      kierr    = 0_iintegers
      izerror  = 0_iintegers
      yzerrmsg = '   '
!
!------------------------------------------------------------------------------
!
!     The scalar product between the nodal vectors pxn of the home
!     node and the 6 (5) surrounding nodes is equal to the cosine of
!     the angle epsilon between the nodes and the home node. Epsilon
!     is a measure of the great circle arc distance of the gridpoints.
!     The calculation is performed only for diamond 1 since the relative
!     distances are the same in all diamonds.
!
      DO j2 = ki2sc,ki2ec
        DO j1  = ki1sc, ki1ec
        zsp1   = pxn(j1,j2,1,1)*pxn(j1+isp11,j2+isp21,1,1) +       &
                 pxn(j1,j2,2,1)*pxn(j1+isp11,j2+isp21,2,1) +       &
                 pxn(j1,j2,3,1)*pxn(j1+isp11,j2+isp21,3,1) 
        zsp1   = ACOS (zsp1)
        pdxmin = MIN (pdxmin, zsp1)
        zdxmax = MAX (zdxmax, zsp1)
        zsp2   = pxn(j1,j2,1,1)*pxn(j1+isp12,j2+isp22,1,1) +       &
                 pxn(j1,j2,2,1)*pxn(j1+isp12,j2+isp22,2,1) +       &
                 pxn(j1,j2,3,1)*pxn(j1+isp12,j2+isp22,3,1)
        zsp2   = ACOS (zsp2)
        pdxmin = MIN (pdxmin, zsp2)
        zdxmax = MAX (zdxmax, zsp2)
        zsp3   = pxn(j1,j2,1,1)*pxn(j1+isp13,j2+isp23,1,1) +       &
                 pxn(j1,j2,2,1)*pxn(j1+isp13,j2+isp23,2,1) +       &
                 pxn(j1,j2,3,1)*pxn(j1+isp13,j2+isp23,3,1)
        zsp3   = ACOS (zsp3)
        pdxmin = MIN (pdxmin, zsp3)
        zdxmax = MAX (zdxmax, zsp3)
        zsp4   = pxn(j1,j2,1,1)*pxn(j1+isp14,j2+isp24,1,1) +       &
                 pxn(j1,j2,2,1)*pxn(j1+isp14,j2+isp24,2,1) +       &
                 pxn(j1,j2,3,1)*pxn(j1+isp14,j2+isp24,3,1)
        zsp4   = ACOS (zsp4)
        pdxmin = MIN (pdxmin, zsp4)
        zdxmax = MAX (zdxmax, zsp4)
        zsp5   = pxn(j1,j2,1,1)*pxn(j1+isp15,j2+isp25,1,1) +       &
                 pxn(j1,j2,2,1)*pxn(j1+isp15,j2+isp25,2,1) +       &
                 pxn(j1,j2,3,1)*pxn(j1+isp15,j2+isp25,3,1)
        zsp5   = ACOS (zsp5)
        pdxmin = MIN (pdxmin, zsp5)
        zdxmax = MAX (zdxmax, zsp5)
        zsp6   = pxn(j1,j2,1,1)*pxn(j1+isp16,j2+isp26,1,1) +       &
                 pxn(j1,j2,2,1)*pxn(j1+isp16,j2+isp26,2,1) +       &
                 pxn(j1,j2,3,1)*pxn(j1+isp16,j2+isp26,3,1)
        zsp6   = ACOS (zsp6)
        pdxmin = MIN (pdxmin, zsp6)
        zdxmax = MAX (zdxmax, zsp6)
        ENDDO
      ENDDO
!
!------------------------------------------------------------------------------
!
!     Now the same procedure for the nodes 1 and 2, 2 and 3 and so on.
!     To avoid problems at the 4 special points (corners of the
!     diamond with only 5 neighbours) the calculation omits the
!     outer row/column.
!
      DO j2 = ki2sc+1,ki2ec-1
        DO j1  = ki1sc+1, ki1ec-1
        zsp12  = pxn(j1+isp11,j2+isp21,1,1)*pxn(j1+isp12,j2+isp22,1,1) +   &
                 pxn(j1+isp11,j2+isp21,2,1)*pxn(j1+isp12,j2+isp22,2,1) +   &
                 pxn(j1+isp11,j2+isp21,3,1)*pxn(j1+isp12,j2+isp22,3,1)  
        zsp12  = ACOS (zsp12)
        pdxmin = MIN (pdxmin, zsp12)
        zdxmax = MAX (zdxmax, zsp12)
        zsp23  = pxn(j1+isp12,j2+isp22,1,1)*pxn(j1+isp13,j2+isp23,1,1) +   &
                 pxn(j1+isp12,j2+isp22,2,1)*pxn(j1+isp13,j2+isp23,2,1) +   &
                 pxn(j1+isp12,j2+isp22,3,1)*pxn(j1+isp13,j2+isp23,3,1)  
        zsp23  = ACOS (zsp23)
        pdxmin = MIN (pdxmin, zsp23)
        zdxmax = MAX (zdxmax, zsp23)
        zsp34  = pxn(j1+isp13,j2+isp23,1,1)*pxn(j1+isp14,j2+isp24,1,1) +   &
                 pxn(j1+isp13,j2+isp23,2,1)*pxn(j1+isp14,j2+isp24,2,1) +   &
                 pxn(j1+isp13,j2+isp23,3,1)*pxn(j1+isp14,j2+isp24,3,1)  
        zsp34  = ACOS (zsp34)
        pdxmin = MIN (pdxmin, zsp34)
        zdxmax = MAX (zdxmax, zsp34)
        zsp45  = pxn(j1+isp14,j2+isp24,1,1)*pxn(j1+isp15,j2+isp25,1,1) +   &
                 pxn(j1+isp14,j2+isp24,2,1)*pxn(j1+isp15,j2+isp25,2,1) +   &
                 pxn(j1+isp14,j2+isp24,3,1)*pxn(j1+isp15,j2+isp25,3,1)  
        zsp45  = ACOS (zsp45)
        pdxmin = MIN (pdxmin, zsp45)
        zdxmax = MAX (zdxmax, zsp45)
        zsp56  = pxn(j1+isp15,j2+isp25,1,1)*pxn(j1+isp16,j2+isp26,1,1) +   &
                 pxn(j1+isp15,j2+isp25,2,1)*pxn(j1+isp16,j2+isp26,2,1) +   &
                 pxn(j1+isp15,j2+isp25,3,1)*pxn(j1+isp16,j2+isp26,3,1)  
        zsp56  = ACOS (zsp56)
        pdxmin = MIN (pdxmin, zsp56)
        zdxmax = MAX (zdxmax, zsp56)
        zsp61  = pxn(j1+isp16,j2+isp26,1,1)*pxn(j1+isp11,j2+isp21,1,1) +   &
                 pxn(j1+isp16,j2+isp26,2,1)*pxn(j1+isp11,j2+isp21,2,1) +   &
                 pxn(j1+isp16,j2+isp26,3,1)*pxn(j1+isp11,j2+isp21,3,1)  
        zsp61  = ACOS (zsp61)
        pdxmin = MIN (pdxmin, zsp61)
        zdxmax = MAX (zdxmax, zsp61)
        ENDDO
      ENDDO
!
!------------------------------------------------------------------------------

  ! Compute global maximum and minimum
  IF (num_compute > 1) THEN
    zdxdhs(1) = zdxmax
    zdxdhs(2) = zdhmax
    CALL global_values (zdxdhs, 2, 'MAX', imp_reals, icomm_cart, -1,      &
                        yzerrmsg, izerror)
    zdxmax = zdxdhs(1)
    zdhmax = zdxdhs(2)

    zdxdhs(1) = pdxmin
    zdxdhs(2) = pdhmin
    CALL global_values (zdxdhs, 2, 'MIN', imp_reals, icomm_cart, -1,      &
                        yzerrmsg, izerror)
    pdxmin = zdxdhs(1)
    pdhmin = zdxdhs(2)
  ENDIF

!------------------------------------------------------------------------------
!
!     Calculation of the minimum height of the triangles; as a good
!     estimate, we consider a spherical triangle with three sides of
!     equal length "pdxmin"; the same is done for "zdhmax".
      pdhmin = ASIN (SIN (pdxmin)*SIN (2.0_ireals*Pid5))
      zdhmax = ASIN (SIN (zdxmax)*SIN (2.0_ireals*Pid5))
!
!     Multiply by the radius of the earth r_earth to get the real distance
!     on the earth
      pdxmin = pdxmin*r_earth
      pdhmin = pdhmin*r_earth
      zdxmax = zdxmax*r_earth
      zdhmax = zdhmax*r_earth
!
!
!------------------------------------------------------------------------------
!
      IF (lzdebug) THEN
        PRINT *, '  SUBROUTINE *dxdhmin*,  minimum mesh width ',    &
                 '(dxmin, unit: m)= ', pdxmin 
        PRINT *, '  SUBROUTINE *dxdhmin*,  minimum triangle',       &
                 ' height (dhmin, unit: m)= ', pdhmin
        PRINT *, '  SUBROUTINE *dxdhmin*,  maximum mesh width ',    &
                 '(dxmax, unit: m)= ', zdxmax 
        PRINT *, '  SUBROUTINE *dxdhmin*,  maximum triangle',       &
                 ' height (dhmax, unit: m)= ', zdhmax
      ENDIF
!
!------------------------------------------------------------------------------

END SUBROUTINE dxdhmin

!==============================================================================
!==============================================================================
!+ Calculates the global arrays for local spherical coordinate system
!------------------------------------------------------------------------------

SUBROUTINE loc_coor (pxn    , kig1sm2, kig1ep2, kig2sm2, kig2ep2,       &
                     knd    , lzdebug,                                  &
                     perlon , perlat , peta   , pchi   ,                &
                     pcpsi  , pspsi  , pbary  , pceneta, pcenchi,       &
                     kispoke, kispokes,ki1mrp , ki2mrp ,                &
                     kig1sm1, kig1ep1, kig2sm1, kig2ep1,                &
                     kig1s  , kig1e  , kig2s  , kig2e  , kierr)

!------------------------------------------------------------------------------
!
! Description:
!   *loc_coor* calculates the global arrays which define the
!   the local spherical coordinate system attached at each gridpoint
!
! Method:
!
!------------------------------------------------------------------------------
!
!     Input 
!     pxn    (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 3, knd) REAL: Cartesian 
!            (x,y,z) coordinates of the location vector of the grid
!            points (nodes) on the unit-sphere; phys. dim. ( - )
!            Attention: pxn is defined for the extended domain, i.e. 
!            with 2 rows and columns around the original diamond array
!            size.
!
!     kig1s   INTEGER  first dimension of arrays, start index 
!     kig1sm1 INTEGER  kig1sm1 = kig1s - 1
!     kig1sm2 INTEGER  kig1sm2 = kig1s - 2
!     kig1e   INTEGER  first dimension of arrays, end index 
!     kig1ep1 INTEGER  kig1ep1 = kig1e + 1
!     kig1ep2 INTEGER  kig1ep2 = kig1e + 2
!     kig2s   INTEGER  second dimension of arrays, start index 
!     kig2sm1 INTEGER  kig2sm1 = kig2s - 1
!     kig2sm2 INTEGER  kig2sm2 = kig2s - 2
!     kig2e   INTEGER  second dimension of arrays, end index 
!     kig2ep1 INTEGER  kig2ep1 = kig2e + 1
!     kig2ep2 INTEGER  kig2ep2 = kig2e + 2
!     knd     INTEGER  number of diamonds 
!     lzdebug  LOGICAL  debug flag, if .true. print debug information
!
!------------------------------------------------------------------------------
!
!     Output 
!     perlon (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 2, knd):  Cartesian 
!            (x,y,z) coordinates of the unit vector of the local easter-
!            ly direction on the unit-sphere
!     perlat (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 3, knd):  Cartesian 
!            (x,y,z) coordinates of the unit vector of the local nor-   
!            therly direction on the unit-sphere
!     peta   (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7):  local longitude
!            of the 6 (5) neighbouring gridpoints relative to the home  
!            node
!     pchi   (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7):  local latitude
!            of the 6 (5) neighbouring gridpoints relative to the home  
!            node
!     pcpsi  (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 6):  sine of the 
!            rotation angle between the local spherical systems
!     pspsi  (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 6):  cosine of the
!            rotation angle between the local spherical systems
!     pbary  (kig1s  :kig1e  , kig2s  :kig2e  , 6):  factor needed  
!            for the computation of the barycentric coordinates
!     pceneta(kig1s  :kig1e  , kig2s  :kig2e  , 6):  local longitude
!            of the neighbouring triangle centers
!     pcenchi(kig1s  :kig1e  , kig2s  :kig2e  , 6):  local latitude
!            of the neighbouring triangle centers
!
!     kispoke (12):   spokes to the 6 (5) neighbours for 2-d arrays
!     kispokes(12,4): spokes to the 6 neighbours for 2-d arrays at the
!            four mirrored points of the extended array
!     ki1mrp  (8):    i1-indices of the four mirrored points
!     ki2mrp  (8):    i2-indices of the four mirrored points
!
!     kierr   INTEGER  error flag, set to 0, if no error occured
!
!
!------------------------------------------------------------------------------

!     Dummy arguments
  INTEGER  (KIND=iintegers)   ::         &
            kig1sm2, kig1ep2, kig2sm2, kig2ep2, knd , kierr,    &
            kig1sm1, kig1ep1, kig2sm1, kig2ep1,                 &
            kig1s  , kig1e  , kig2s  , kig2e
!
  REAL     (KIND=ireals)      ::         &
            pxn     (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 3, knd), &
            perlon  (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 2, knd), &
            perlat  (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 3, knd), &
            peta    (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7),      &
            pchi    (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7),      &
            pcpsi   (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 6),      &
            pspsi   (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 6),      &
            pbary   (kig1s  :kig1e  , kig2s  :kig2e  , 6),      &
            pceneta (kig1s  :kig1e  , kig2s  :kig2e  , 6),      &
            pcenchi (kig1s  :kig1e  , kig2s  :kig2e  , 6)
!
  INTEGER  (KIND=iintegers)   ::      &
            kispoke (12), kispokes(12,4), ki1mrp(8), ki2mrp(8)
!
  LOGICAL  lzdebug
!
!     Local
  REAL     (KIND=ireals)      ::      &
            zd12, & ! scalar product of two vectors
            zd13, & ! scalar product of two vectors
            zd22, & ! scalar product of two vectors
            zd23, & ! scalar product of two vectors
            z1  , & ! local temporary
            z2  , & ! local temporary
            zn  , & ! local norm
            zsn , & ! hemisphere discriminator
            zchi, & ! pi/5
            zeta, & ! meridian angle of diamond
            zd1 , & ! scalar product of two vectors
            zd2 , & ! scalar product of two vectors
            zd3     ! scalar product of two vectors
!
  INTEGER  (KIND=iintegers)   ::      &
            j1, j2, jd, jk ,jm, js  ! Loop indices
  INTEGER  (KIND=iintegers)   ::      &
            js1, js2                ! Indices of a gridpoint
!
!------------------------------------------------------------------------------
!
!     Subroutines/functions used
!     EXTERNAL  setspokes, ! computes the spokes, i.e. the offsets 
!                            for 2-d array addressing
!    1          setmrp,    ! defines the i1-, i2-indices of the four
!                            mirrored points of the extended array
!    2          printrarr  ! prints a real 2-dimensional array
!
!------------------------------------------------------------------------------

!     1. Preset "ispoke" and "ispokes", the offsets of the 6 (5) neigh-
!        bouring gridpoints relative to the central node for 2-d array
!        addressing. "ispokes" will be used at the four mirrored points
!        of arrays which have been extended by one row/column since 
!        there the normal spokes "ispoke" will take wrong values due to
!        the pentagonal structure at the corners of the diamonds.
!        Define the i1-, i2-indices of the four mirrored points of the
!        extended array.
!
      CALL setspokes (kispoke, kispokes, lzdebug, kierr)
!
      IF (kierr .NE. 0) THEN
        PRINT *,'  Error in subroutine *loc_coor* calling *setspokes*'
        RETURN
      ENDIF
!
      CALL setmrp (ki1mrp, ki2mrp, kig1s, kig1e, kig2s, kig2e,      &
                   lzdebug, kierr)
!
      IF (kierr .NE. 0) THEN
        PRINT *,'  Error in subroutine *loc_coor* calling *setmrp*'
        RETURN
      ENDIF
!
!------------------------------------------------------------------------------
!
!     2. Compute the Cartesian coordinates of the local unit vectors 
!        "prlon" (aligned to the global east direction) and "prlat"
!        (aligned to the global north direction)
!
      DO jd = 1, knd                 ! Loop over all diamonds
!
        DO j2 = kig2sm2, kig2ep2     ! Loop over the extended array
          DO j1 = kig1sm2, kig1ep2   ! Loop over the extended array
          z1    = -pxn(j1,j2,2,jd)   ! y-coordinate of node vector
          z2    =  pxn(j1,j2,1,jd)   ! x-coordinate of node vector
          zn    = 1./SQRT (z1**2 + z2**2 + 1.e-30)
!
          perlon(j1,j2,1,jd) = z1*zn ! Normalize to unit length
          perlon(j1,j2,2,jd) = z2*zn ! Normalize to unit length
!                                      z-coordinate of perlon is 0
!         perlat is the cross product of pxn and perlon
!
          perlat(j1,j2,1,jd) = -pxn(j1,j2,3,jd)*z2*zn
          perlat(j1,j2,2,jd) =  pxn(j1,j2,3,jd)*z1*zn
          perlat(j1,j2,3,jd) =  pxn(j1,j2,1,jd)*z2*zn - pxn(j1,j2,2,jd)*z1*zn
          ENDDO
        ENDDO
!
!        Define the local unit vectors for the north and south pole
!        separately (zn is undefined there). 
!        Attention:  The local unit vectors perlon, perlat for the
!                    poles differ from one diamond to the next!
        zchi = 0.8*ASIN(1.0_ireals)       ! 2*pi/5
        zsn  = 1.0_ireals                 ! hemisphere descriminator
        IF (jd .GE. 6) zsn = -1.0_ireals
        zeta = (jd-1.5)*zchi
        IF (jd .GE. 6) zeta = (jd-6)*zchi
!       
        IF(kig1s==igg1s .AND. kig2s==igg2s) THEN
          perlon(0,1,1,jd) = - SIN (zeta)
          perlon(0,1,2,jd) =   COS (zeta)
!
          perlat(0,1,1,jd) = - COS (zeta)*zsn
          perlat(0,1,2,jd) = - SIN (zeta)*zsn
          perlat(0,1,3,jd) = 0.
        ENDIF
!
      ENDDO     ! End loop over the diamonds

!------------------------------------------------------------------------------
!
!     3. Compute the local spherical coordinates "peta" (local longi-
!        tude) and "pchi" (local latitude) of the 6 (5) surrounding
!        gridpoints for each node in diamond 1 (plus one extra row and
!        column around it); compute the sine and cosine of the rotation
!        angle ("pspsi", "pcpsi") between the local coordinate systems
!        of the 6 (5) neighbouring gridpoints and the central node.
!
!        All computations are performed for diamond 1 only since the
!        relative distances are the same for the other 9 diamonds except
!        at the pole points of the diamonds (see section 8.)
!     
      DO jm = 1,6          ! Loop over the 6 (5) neighbours
        DO j2 = kig2sm1, kig2ep1
          js2 = j2 + kispoke(jm+6)
          DO j1 = kig1sm1, kig1ep1
          js1 = j1 + kispoke(jm)
!
!        Scalar product of the location vectors
          zd1 = pxn   (j1,j2,1,1)*pxn(js1,js2,1,1) +    &
                pxn   (j1,j2,2,1)*pxn(js1,js2,2,1) +    &
                pxn   (j1,j2,3,1)*pxn(js1,js2,3,1)   
!
!        Scalar product of the local longitude vector and the location
!        vector 
          zd2 = perlon(j1,j2,1,1)*pxn(js1,js2,1,1) +    &
                perlon(j1,j2,2,1)*pxn(js1,js2,2,1)  
!
!        Scalar product of the local latitude vector and the location
!        vector 
          zd3 = perlat(j1,j2,1,1)*pxn(js1,js2,1,1) +    &
                perlat(j1,j2,2,1)*pxn(js1,js2,2,1) +    &
                perlat(j1,j2,3,1)*pxn(js1,js2,3,1)  
!
!        Avoid values GT. 1 and LT. -1 which are due to round off errors
          zd3 = SIGN ( MIN (1.0_ireals, ABS (zd3)), zd3)
!
!        Local spherical coordinates "peta" and "pchi" for the diamond
!        core 
          peta (j1,j2,jm) = ATAN (zd2/(zd1 + 1.e-20_ireals))
! RJ: Abfrage auf 0 und 1 ist aber unschoen !!!!
!         IF (j1 .EQ. 0 .AND. j2 .EQ. 1) peta (j1,j2,jm) = ASIN (zd2)
          IF (j1 .EQ. igg1s .AND. j2 .EQ. igg2s)     &
            peta (j1,j2,jm) = ASIN (zd2)
          pchi (j1,j2,jm) = ASIN (zd3)
!
!        Rotation angles "pspsi" and "pcpsi"
!        Scalar product of the location vector and the local longitude
!        vector at the surrounding gridpoints
          zd12 = pxn    (j1,j2,1,1)*perlon(js1,js2,1,1) +     &
                 pxn    (j1,j2,2,1)*perlon(js1,js2,2,1) 
!
!        Scalar product of the location vector and the local latitude
!        vector at the surrounding gridpoints
          zd13 = pxn    (j1,j2,1,1)*perlat(js1,js2,1,1) +     &
                 pxn    (j1,j2,2,1)*perlat(js1,js2,2,1) +     &
                 pxn    (j1,j2,3,1)*perlat(js1,js2,3,1) 
!
!        Scalar product of the local longitude vectors at the home node
!        and at the surrounding gridpoints
          zd22 = perlon (j1,j2,1,1)*perlon(js1,js2,1,1) +     &
                 perlon (j1,j2,2,1)*perlon(js1,js2,2,1) 
!
!        Scalar product of the local longitude vector at the home node
!        and the latitude vector at the surrounding gridpoints
          zd23 = perlon (j1,j2,1,1)*perlat(js1,js2,1,1) +     &
                 perlon (j1,j2,2,1)*perlat(js1,js2,2,1) 
!
!        Rotation angles "pspsi" and "pcpsi"
          pcpsi (j1,j2,jm) = COS (peta(j1,j2,jm))*zd22 -      &
                             SIN (peta(j1,j2,jm))*zd12 
          pspsi (j1,j2,jm) = COS (peta(j1,j2,jm))*zd23 -      &
                             SIN (peta(j1,j2,jm))*zd13 
!
          ENDDO
        ENDDO
!
      ENDDO     ! End of loop over the 6 (5) neighbours

!------------------------------------------------------------------------------
!
!     4. Correct the "peta", "pchi", "pcpsi", "pspsi" at the four 
!        mirrored points of the extended array
!     
      DO js = 1,4    ! Loop over the four mirrored points
        j1 = ki1mrp(js)
        j2 = ki2mrp(js)
! RJ: Do that only if we own the point
        IF (j1 < -10) CYCLE
!
        DO jm = 1,6  ! Loop over the six neigbours
          js1 = j1 + kispokes(jm  ,js)
          js2 = j2 + kispokes(jm+6,js)
!
!        Scalar product of the location vectors
          zd1 = pxn   (j1,j2,1,1)*pxn(js1,js2,1,1) +     &
                pxn   (j1,j2,2,1)*pxn(js1,js2,2,1) +     &
                pxn   (j1,j2,3,1)*pxn(js1,js2,3,1)   
!
!        Scalar product of the local longitude vector and the location
!        vector 
          zd2 = perlon(j1,j2,1,1)*pxn(js1,js2,1,1) +     &
                perlon(j1,j2,2,1)*pxn(js1,js2,2,1)  
!
!        Scalar product of the local latitude vector and the location
!        vector 
          zd3 = perlat(j1,j2,1,1)*pxn(js1,js2,1,1) +     &
                perlat(j1,j2,2,1)*pxn(js1,js2,2,1) +     &
                perlat(j1,j2,3,1)*pxn(js1,js2,3,1)  
!
!        Avoid values GT. 1 and LT. -1 which are due to round off errors
          zd3 = SIGN ( MIN (1.0_ireals, ABS (zd3)), zd3)
!
!        Local spherical coordinates "peta" and "pchi" for the diamond
!        core 
          peta (j1,j2,jm) = ATAN (zd2/(zd1 + 1.e-20))
          pchi (j1,j2,jm) = ASIN (zd3)
!
!        Rotation angles "pspsi" and "pcpsi"
!        Scalar product of the location vector and the local longitude
!        vector at the surrounding gridpoints
          zd12 = pxn    (j1,j2,1,1)*perlon(js1,js2,1,1) +    &
                 pxn    (j1,j2,2,1)*perlon(js1,js2,2,1) 
!
!        Scalar product of the location vector and the local latitude
!        vector at the surrounding gridpoints
          zd13 = pxn    (j1,j2,1,1)*perlat(js1,js2,1,1) +    &
                 pxn    (j1,j2,2,1)*perlat(js1,js2,2,1) +    &
                 pxn    (j1,j2,3,1)*perlat(js1,js2,3,1) 
!
!        Scalar product of the local longitude vectors at the home node
!        and at the surrounding gridpoints
          zd22 = perlon (j1,j2,1,1)*perlon(js1,js2,1,1) +    &
                 perlon (j1,j2,2,1)*perlon(js1,js2,2,1) 
!
!        Scalar product of the local longitude vector at the home node
!        and the latitude vector at the surrounding gridpoints
          zd23 = perlon (j1,j2,1,1)*perlat(js1,js2,1,1) +    &
                 perlon (j1,j2,2,1)*perlat(js1,js2,2,1) 
!
!        Rotation angles "pspsi" and "pcpsi"
          pcpsi (j1,j2,jm) = COS (peta(j1,j2,jm))*zd22 -    &
                             SIN (peta(j1,j2,jm))*zd12 
          pspsi (j1,j2,jm) = COS (peta(j1,j2,jm))*zd23 -    &
                             SIN (peta(j1,j2,jm))*zd13 
!
!        Copy results to the mirror point
         peta (ki1mrp(js+4),ki2mrp(js+4),jm) = peta (j1,j2,jm)
         pchi (ki1mrp(js+4),ki2mrp(js+4),jm) = pchi (j1,j2,jm)
         pcpsi(ki1mrp(js+4),ki2mrp(js+4),jm) = pcpsi(j1,j2,jm)
         pspsi(ki1mrp(js+4),ki2mrp(js+4),jm) = pspsi(j1,j2,jm)
!
        ENDDO  ! End of loop over the six neighbours
!
      ENDDO    ! End of loop over the four mirrored points
!
!------------------------------------------------------------------------------
!
!     5. Compute the factor "pbary" which is needed for the calculation
!        of the barycentric coordinates of the gridpoints which are 
!        used for the interpolation routines
!
      DO jm = 1,6
        jk = jm + 1
        IF (jk .EQ. 7) jk = 1
        DO j2 = kig2s, kig2e
          DO j1 = kig1s, kig1e
          z1 = peta(j1,j2,jm)*pchi(j1,j2,jk) - peta(j1,j2,jk)*pchi(j1,j2,jm)
          IF (z1 .EQ. 0.) THEN
            pbary (j1,j2,jm) = 0.
          ELSE
            pbary (j1,j2,jm) = 1./z1
          ENDIF
          ENDDO
        ENDDO
      ENDDO    ! End of loop over the 6 (5) neighbours
!
!------------------------------------------------------------------------------
!
!     6. Repeat the value for the neighbour "1" (lower left) at index 
!        "7" for the arrays "peta" and "pchi"
      DO j2 = kig2sm1, kig2ep1
        DO j1 = kig1sm1, kig1ep1
        peta  (j1,j2,7) = peta  (j1,j2,1)
        pchi  (j1,j2,7) = pchi  (j1,j2,1)
        ENDDO
      ENDDO
!
!------------------------------------------------------------------------------
!
!     7. Compute the normalized spherical coordinates ("pceneta" and 
!        "pcenchi") of the triangle centers (only for the core of dia-
!        mond 1). These coordinates are needed by the subroutine
!        *findtri* to spot the triangle which contains the departure
!        point or midpoint of the particle trajectory.
!        Take care of those points with only 5 neighbours.
!
      DO jm = 1,6
        DO j2 = kig2s, kig2e
          DO j1 = kig1s, kig1e
!
          IF ( (peta (j1,j2,jm) - peta (j1,j2, jm+1) .EQ. 0.) .AND.     &
               (pchi (j1,j2,jm) - pchi (j1,j2, jm+1) .EQ. 0.) ) THEN
            pceneta (j1,j2,jm) = 0.
            pcenchi (j1,j2,jm) = 0.
          ELSE
            z1 = peta (j1,j2,jm) + peta (j1,j2, jm+1)
            z2 = pchi (j1,j2,jm) + pchi (j1,j2, jm+1)
            zn = 1./SQRT (z1**2 + z2**2)
            pceneta (j1,j2,jm) = z1*zn
            pcenchi (j1,j2,jm) = z2*zn
          ENDIF
!
          ENDDO
        ENDDO
      ENDDO
!
!------------------------------------------------------------------------------

      kierr = 0

!------------------------------------------------------------------------------

END SUBROUTINE loc_coor

!==============================================================================
!==============================================================================
!+ Calculates areas of triangles in diamond #1
!------------------------------------------------------------------------------

SUBROUTINE tri_area (pxn    , kig1sm2, kig1ep2, kig2sm2, kig2ep2,    &
                     knd    , lzdebug,                               &
                     parea  , kig1s  , kig1ep1, kig2sm1, kig2e  ,    &
                     kierr)

!------------------------------------------------------------------------------
!
! Description:
!   *tri_area* calculates the areas "area" of the triangles in diamond
!   #1 based on the formula of Hulier. The area is computed for the 
!   core of diamond #1 plus one additional row/column around. Thus 
!   there are 2*(ni+2)**2 triangles whose indices run from 0 to ni+1
!   and 0 to ni+1 and from 1 to 2. "1" are the upward pointing 
!   triangles, "2" the downward pointing ones.
!
! Method:
!
!------------------------------------------------------------------------------
!
!     Input: 
!     pxn    (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 3, knd) REAL: Cartesian 
!            (x,y,z) coordinates of the location vector of the grid
!            points (nodes) on the unit-sphere; phys. dim. ( - )
!            Attention: pxn is defined for the extended domain, i.e. 
!            with 2 rows and columns around the original diamond array
!            size.
!
!     kig1s   INTEGER  first dimension of arrays, start index (kig1s = 0)
!     kig1sm2 INTEGER  kig1sm2 = kig1s - 2
!     kig1e   INTEGER  first dimension of arrays, end index (kig1e = ni)
!     kig1ep1 INTEGER  kig1ep1 = kig1e + 1
!     kig1ep2 INTEGER  kig1ep2 = kig1e + 2
!     kig2s   INTEGER  sec. dimension of arrays, start index (kig2s = 1)
!     kig2sm1 INTEGER  kig2sm1 = kig2s - 1
!     kig2sm2 INTEGER  kig2sm2 = kig2s - 2
!     kig2e   INTEGER  sec. dimension of arrays, end index (kig2e = ni+1)
!     kig2ep2 INTEGER  kig2ep2 = kig2e + 2
!     knd     INTEGER  number of diamonds (knd = 10)
!     lzdebug  LOGICAL  debug flag, if .true. print debug information
!
!------------------------------------------------------------------------------
!
!     Output:
!     parea  (kig1s  :kig1ep1, kig2sm1:kig2e, 2):  area of the triangles
!            on the unit-sphere; phys. dim. ( - )
!
!     kierr   INTEGER  error flag, set to 0, if no error occured
!
!
!------------------------------------------------------------------------------

! Dummy arguments
  INTEGER  (KIND=iintegers)   ::                                  &
            kig1sm2, kig1ep2, kig2sm2, kig2ep2, knd , kierr,      &
            kig1s  , kig1ep1, kig2sm1, kig2e

  REAL     (KIND=ireals)      ::                                  &
            pxn    (kig1sm2:kig1ep2, kig2sm2:kig2ep2, 3, knd),    &
            parea  (kig1s  :kig1ep1, kig2sm1:kig2e  , 2)

  LOGICAL  lzdebug

! Local arrays and variables
  REAL     (KIND=ireals)   ::     &
            zxne(3,3), & ! Cartesian coordinates of a triangle
            zt1,       & !
            zt2,       & !
            zt3,       & !
            zs ,       & !
            zw           !
!
  INTEGER  (KIND=iintegers)   ::  j1 , j2      ! DO loop indices
!
!------------------------------------------------------------------------------
!
      kierr = 0
      DO j2 = kig2sm1,kig2e
         DO j1= kig1s,kig1ep1
!
!     Upward pointing triangles (index "1")
!
         zxne(1,1) = pxn(j1  ,j2  ,1,1)
         zxne(1,2) = pxn(j1  ,j2  ,2,1)
         zxne(1,3) = pxn(j1  ,j2  ,3,1)
         zxne(2,1) = pxn(j1-1,j2+1,1,1)
         zxne(2,2) = pxn(j1-1,j2+1,2,1)
         zxne(2,3) = pxn(j1-1,j2+1,3,1)
         zxne(3,1) = pxn(j1-1,j2  ,1,1)
         zxne(3,2) = pxn(j1-1,j2  ,2,1)
         zxne(3,3) = pxn(j1-1,j2  ,3,1)
         zt1 = zxne(1,1)*zxne(2,1) + zxne(1,2)*zxne(2,2)      &
             + zxne(1,3)*zxne(2,3)
         zt2 = zxne(2,1)*zxne(3,1) + zxne(2,2)*zxne(3,2)      &
             + zxne(2,3)*zxne(3,3)
         zt3 = zxne(3,1)*zxne(1,1) + zxne(3,2)*zxne(1,2)      &
             + zxne(3,3)*zxne(1,3)
!
!     Avoid values GT. 1 and LT. -1 which are due to round off errors
         zt1 = SIGN ( MIN (1.0_ireals, ABS (zt1)), zt1)
         zt2 = SIGN ( MIN (1.0_ireals, ABS (zt2)), zt2)
         zt3 = SIGN ( MIN (1.0_ireals, ABS (zt3)), zt3)

         zt1 = 0.5_ireals*ACOS(zt1)
         zt2 = 0.5_ireals*ACOS(zt2)
         zt3 = 0.5_ireals*ACOS(zt3)
         zs  = 0.5_ireals*(zt1 + zt2 + zt3)
         zw  = TAN(zs)*TAN(zs-zt1)*TAN(zs-zt2)*TAN(zs-zt3)
         parea(j1,j2,1) = 4.0_ireals*ATAN(SQRT(zw))
!
!     Downward pointing triangles (index "2")
!
         zxne(1,1) = pxn(j1  ,j2  ,1,1)
         zxne(1,2) = pxn(j1  ,j2  ,2,1)
         zxne(1,3) = pxn(j1  ,j2  ,3,1)
         zxne(2,1) = pxn(j1  ,j2+1,1,1)
         zxne(2,2) = pxn(j1  ,j2+1,2,1)
         zxne(2,3) = pxn(j1  ,j2+1,3,1)
         zxne(3,1) = pxn(j1-1,j2+1,1,1)
         zxne(3,2) = pxn(j1-1,j2+1,2,1)
         zxne(3,3) = pxn(j1-1,j2+1,3,1)
         zt1 = zxne(1,1)*zxne(2,1) + zxne(1,2)*zxne(2,2)      &
             + zxne(1,3)*zxne(2,3)
         zt2 = zxne(2,1)*zxne(3,1) + zxne(2,2)*zxne(3,2)      &
             + zxne(2,3)*zxne(3,3)
         zt3 = zxne(3,1)*zxne(1,1) + zxne(3,2)*zxne(1,2)      &
             + zxne(3,3)*zxne(1,3)
!
!     Avoid values GT. 1 and LT. -1 which are due to round off errors
         zt1 = SIGN ( MIN (1.0_ireals, ABS (zt1)), zt1)
         zt2 = SIGN ( MIN (1.0_ireals, ABS (zt2)), zt2)
         zt3 = SIGN ( MIN (1.0_ireals, ABS (zt3)), zt3)

         zt1 = 0.5_ireals*ACOS(zt1)
         zt2 = 0.5_ireals*ACOS(zt2)
         zt3 = 0.5_ireals*ACOS(zt3)
         zs  = 0.5_ireals*(zt1 + zt2 + zt3)
         zw  = TAN(zs)*TAN(zs-zt1)*TAN(zs-zt2)*TAN(zs-zt3)
         parea(j1,j2,2) = 4.0_ireals*ATAN(SQRT(zw))
 
         ENDDO 
      ENDDO 
!
!------------------------------------------------------------------------------
!
!     Due to the 4 special points at the vertices of diamond #1 with
!     only 5 neighbours, the following 6 triangles are undefined and
!     their areas set to 0.
      IF (kig1s   == igg1s   .AND. kig2sm1 == igg2s-1) THEN
        parea(kig1s  ,kig2sm1,1) = 0.
        parea(kig1s  ,kig2sm1,2) = 0.
      ENDIF
      IF (kig1ep1 == igg1e+1 .AND. kig2e   == igg2e  ) THEN
        parea(kig1ep1,kig2e  ,1) = 0.
        parea(kig1ep1,kig2e  ,2) = 0.
      ENDIF
      IF (kig1ep1 == igg1e+1 .AND. kig2sm1 == igg2s-1) THEN
        parea(kig1ep1,kig2sm1,1) = 0.
      ENDIF
      IF (kig1s   == igg1s   .AND. kig2e   == igg2e  ) THEN
        parea(kig1s  ,kig2e  ,1) = 0.
      ENDIF
!
!------------------------------------------------------------------------------

END SUBROUTINE tri_area

!==============================================================================
!==============================================================================
!+ Calculates the global arrays for local spherical coordinate system
!------------------------------------------------------------------------------

SUBROUTINE hex_area (parea  ,                                 &
                     kig1s  , kig1e  , kig2s  , kig2e  ,      &
                     kig1sm1, kig1ep1, kig2sm1, kig2ep1,      &
                     knd    , lzdebug,                        &
                     prarn  , phexwgt, kierr)

!------------------------------------------------------------------------------
!
! Description:
!   *hex_area* calculates the reciprocal "prarn" of the areas of the
!   hexagons related to the gridpoints of diamond #1 which has been 
!   extended by one row and column around, and it computes the
!   weights ("phexwgt", area of the hexagons) needed for the 
!   diagnostics.
!
! Method:
!
!------------------------------------------------------------------------------
!
!     Input
!     parea  (kig1s  :kig1ep1, kig2sm1:kig2e, 2) REAL:  area of the 
!            triangles on the unit-sphere (only diamond #1);
!            phys. dim. ( - )
!
!     kig1s   INTEGER  first dimension of arrays, start index (kig1s = 0)
!     kig1sm1 INTEGER  kig1sm1 = kig1s - 1
!     kig1e   INTEGER  first dimension of arrays, end index (kig1e = ni)
!     kig1ep1 INTEGER  kig1ep1 = kig1e + 1
!     kig2s   INTEGER  sec. dimension of arrays, start index (kig2s = 1)
!     kig2sm1 INTEGER  kig2sm1 = kig2s - 1
!     kig2e   INTEGER  sec. dimension of arrays, end index (kig2e = ni+1)
!     kig2ep1 INTEGER  kig2ep1 = kig2e + 1
!     knd     INTEGER  number of diamonds 
!     lzdebug  LOGICAL  debug flag, if .true. print debug information
!
!------------------------------------------------------------------------------
!
!     Output:
!     prarn  (kig1sm1:kig1ep1, kig2sm1:kig2ep1) REAL:   reciprocal of the
!            area of the hexagons related to the gridpoints on the 
!            unit sphere (only diamond #1); phys. dim. ( - )
!
!     phexwgt(kig1s  :kig1e  , kig2s  :kig2e  ) REAL:   weights of the
!            gridpoints for diagnostic calculations on the unit sphere
!            unit sphere (only diamond #1); phys. dim. ( - )
!
!     kierr   INTEGER  error flag, set to 0, if no error occured
!
!------------------------------------------------------------------------------

! Dummy arguments
  INTEGER  (KIND=iintegers)   ::                                &
            kig1sm1, kig1ep1, kig2sm1, kig2ep1, kierr,          &
            kig1s  , kig1e  , kig2s  , kig2e  , knd
!
  REAL     (KIND=ireals)      ::                                &
            parea  (kig1s  :kig1ep1, kig2sm1:kig2e  , 2),       &
            prarn  (kig1sm1:kig1ep1, kig2sm1:kig2ep1)   ,       &
            phexwgt(kig1s  :kig1e  , kig2s  :kig2e  )
!
  LOGICAL  lzdebug
!
!     Local arrays and variables
  REAL     (KIND=ireals)      ::                                &
            zrarn  (kig1s-2:kig1e+2, kig2s-2:kig2e+2,knd)
!
  INTEGER  (KIND=iintegers)   :: j1 , j2, jd           ! DO loop indices
!
!------------------------------------------------------------------------------
!
      kierr = 0
!
!     Compute the area of the hexagons for the core of diamond #1
      DO j2 = kig2s,kig2e
         DO j1= kig1s,kig1e
         prarn(j1,j2) = 3./( parea(j1+1,j2  ,1) + parea(j1  ,j2-1,2) +   &
                             parea(j1  ,j2  ,1) + parea(j1  ,j2  ,2) +   &
                             parea(j1+1,j2-1,1) + parea(j1+1,j2-1,2) )
!
         ENDDO 
      ENDDO 
!
!     The 4 special points at the vertices of diamond #1 have only 5
!     neighbours, therefore only 5 triangles have to be considered
      IF (kig1s == igg1s .AND. kig2s == igg2s) THEN
        prarn(kig1s, kig2s) = 3./(5.*parea(kig1s+1, kig2s  ,1))
      ENDIF
      IF (kig1e == igg1e .AND. kig2s == igg2s) THEN
        prarn(kig1e, kig2s) = 3./(5.*parea(kig1e  , kig2s  ,1))
      ENDIF
      IF (kig1e == igg1e .AND. kig2e == igg2e) THEN
        prarn(kig1e, kig2e) = 3./(5.*parea(kig1e  , kig2e-1,2))
      ENDIF
      IF (kig1s == igg1s .AND. kig2e == igg2e) THEN
        prarn(kig1s, kig2e) = 3./(5.*parea(kig1s+1, kig2e-1,1))
      ENDIF
!
!------------------------------------------------------------------------------
!
!     Now extend the array prarn by one row/column around the original
!     core of diamond #1; the subroutine *xd_p* will do this, it will be
!     provided with 10 diamonds which are all the same
      DO jd = 1,knd
        DO j2 = kig2s, kig2e
          DO j1 = kig1s, kig1e
          zrarn(j1,j2,jd) = prarn(j1,j2)
          ENDDO
        ENDDO
      ENDDO
!
      CALL xd_p( zrarn, kig1s-2, kig1e+2, kig2s-2, kig2e+2, 1, 1, 1, 2,   &
                 my_cart_id, num_compute, icomm_cart, imp_reals)
!
      prarn(:,:) = zrarn(kig1s-1:kig1e+1,kig2s-1:kig2e+1,1)
!
!------------------------------------------------------------------------------
!
!     Calculate the weight assigned to each gridpoint (which is the
!     area of the hexagon related to this point. At the pole point,
!     the area is divided by 5 since the pole is contained in 5
!     diamonds. Only the inner core of the diamond # 1 is taken to
!     avoid counting gridpoints twice, i.e. for j1=0,ni and j2=ni+1
!     as well as for j1=0 and j2=1,ni+1, the weight is set to 0.
      DO j2 = kig2s, kig2e
        DO j1 = kig1s, kig1e

        phexwgt(j1,j2) = 1./prarn(j1,j2)
        IF (j1 .EQ. igg1s .OR.  j2.EQ. igg2e) phexwgt(j1,j2) = 0.
        IF (j1 .EQ. igg1s .AND. j2.EQ. igg2s) THEN
          phexwgt(j1,j2) =  1./(5.*prarn(j1,j2))
        ENDIF

        ENDDO
      ENDDO
!
!------------------------------------------------------------------------------

END SUBROUTINE hex_area

!==============================================================================
!==============================================================================
!+ Calculates the gradient and Laplacian operators.
!------------------------------------------------------------------------------

SUBROUTINE gen_grd (peta   , pchi   ,                              &
                    kig1sm1, kig1ep1, kig2sm1, kig2ep1,            &
                    kig1s  , kig1e  , kig2s  , kig2e  ,            &
                    lzdebug, pgrd   , prlap  ,  kierr)

!------------------------------------------------------------------------------
!
! Description:
!   *gen_grd* calculates the gradient operator for diamond #1 which
!   has been extended by one row/column around. It is possible to use
!   the same operator for the other diamonds (#2 to #10), too, since
!   the relative positions are the same in all diamonds.
!   Additionally, *gen_grd* computes the Laplacian operator used for
!   the horizontal diffusion (linear forth order), again for the
!   diamond #1 only.
!
! Method:
!
!------------------------------------------------------------------------------
!
!     Input:
!     peta   (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7):  local longitude
!            of the 6 (5) neighbouring gridpoints relative to the home  
!            node
!     pchi   (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7):  local latitude
!            of the 6 (5) neighbouring gridpoints relative to the home  
!            node
!
!     kig1s   INTEGER  first dimension of arrays, start index 
!     kig1sm1 INTEGER  kig1sm1 = kig1s - 1
!     kig1e   INTEGER  first dimension of arrays, end index 
!     kig1ep1 INTEGER  kig1ep1 = kig1e + 1
!     kig2s   INTEGER  second dimension of arrays, start index 
!     kig2sm1 INTEGER  kig2sm1 = kig2s - 1
!     kig2e   INTEGER  second dimension of arrays, end index 
!     kig2ep1 INTEGER  kig2ep1 = kig2e + 1
!     lzdebug  LOGICAL  debug flag, if .true. print debug information
!
!------------------------------------------------------------------------------
!
!     Output:
!     pgrd   (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7, 2):  gradient ope-
!            rator in eta ("1") and chi ("2") direction
!     prlap  (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7   ):  Laplacian ope-
!            rator
!
!     kierr   INTEGER  error flag, set to 0, if no error occured
!
!------------------------------------------------------------------------------
!
! Dummy arguments
  INTEGER  (KIND=iintegers)   ::       &
            kig1sm1, kig1ep1, kig2sm1, kig2ep1, kierr,           &
            kig1s  , kig1e  , kig2s  , kig2e

  REAL     (KIND=ireals)      ::       &
            peta    (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7),       &
            pchi    (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7),       &
            pgrd    (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7, 2),    &
            prlap   (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7)

  LOGICAL  lzdebug
!
! Local arrays and variables
  REAL     (KIND=ireals)      ::     &
            za(5,5),     & ! matrix for the calculation of the gradient op.
            zb(5)          ! r.h.s. of equation

  REAL     (KIND=ireals)      ::     &
            zs, zt, zdot

  INTEGER  (KIND=iintegers)   ::     &
            j1, j2, jz1, jz2, jm ,   & ! Loop indices
            n1, n2, n3, n1m1

!------------------------------------------------------------------------------
!

!     Preset gradient and Laplacian operators with 0
      DO jm = 1,7
        DO j2 = kig2sm1,kig2ep1
          DO j1 = kig1sm1,kig1ep1
          pgrd (j1,j2,jm,1) = 0.
          pgrd (j1,j2,jm,2) = 0.
          prlap(j1,j2,jm)   = 0.
          ENDDO
        ENDDO
      ENDDO
!
!     Compute the gradient and Laplacian operators for the core of 
!     the diamond #1 plus one row and column around
      DO 5 j2 = kig2sm1,kig2ep1
        DO 6 j1 = kig1sm1,kig1ep1
 
!     Preset za with 0
          DO jz1 = 1,5
            DO jz2 = 1,5
            za(jz1,jz2) = 0.
            ENDDO
          ENDDO
!
!     Loop over the 6 (5) neighbours
          DO 10 jm = 1,6
!

!     Take care of the vertices of diamond #1 with only 5 neighbours
          IF (j1 .EQ. igg1s .AND. j2 .EQ. igg2s         &
                            .AND. jm .EQ. 5) GO TO 10
          IF (j1 .EQ. igg1e .AND. j2 .EQ. igg2s         &
                            .AND. jm .EQ. 6) GO TO 10
          IF (j1 .EQ. igg1s .AND. j2 .EQ. igg2e         &
                            .AND. jm .EQ. 4) GO TO 10
          IF (j1 .EQ. igg1e .AND. j2 .EQ. igg2e         &
                            .AND. jm .EQ. 1) GO TO 10
!
!     Take care of the two undefinded points in the extended array
          IF (j1 .EQ. igg1s-1 .AND. j2 .EQ. igg2s-1) GO TO 6 
          IF (j1 .EQ. igg1e+1 .AND. j2 .EQ. igg2e+1) GO TO 6 
!
!     Only the upper triangle plus the diagonal of the symmetric
!     matrix "za" has to be calculated
          za(1,1) = za(1,1) + peta(j1,j2,jm)**2
          za(1,2) = za(1,2) + peta(j1,j2,jm)   *pchi(j1,j2,jm)
          za(1,3) = za(1,3) + peta(j1,j2,jm)**3
          za(1,4) = za(1,4) + peta(j1,j2,jm)**2*pchi(j1,j2,jm)
          za(1,5) = za(1,5) + peta(j1,j2,jm)   *pchi(j1,j2,jm)**2
          za(2,2) = za(2,2) +                   pchi(j1,j2,jm)**2
          za(2,5) = za(2,5) +                   pchi(j1,j2,jm)**3
          za(3,3) = za(3,3) + peta(j1,j2,jm)**4
          za(3,4) = za(3,4) + peta(j1,j2,jm)**3*pchi(j1,j2,jm)
          za(3,5) = za(3,5) + peta(j1,j2,jm)**2*pchi(j1,j2,jm)**2
          za(4,5) = za(4,5) + peta(j1,j2,jm)   *pchi(j1,j2,jm)**3
          za(5,5) = za(5,5) +                   pchi(j1,j2,jm)**4
!
   10     CONTINUE
!
          za(2,3) = za(1,4)
          za(2,4) = za(1,5)
          za(4,4) = za(3,5)
!
!   Invert the matrix
          kierr = 0

          ! This code replaces call to linpack-routine spofa:
          ! The old matrix za is symmetric and positive definit and only the
          ! upper triangle and the diagonal have been set. After the DO n1-loop
          ! za(new) is an upper triangular matrix with 
          !                 za(old) = Transpose(za(new)) * za(new)
          DO n1 = 1,5
            kierr = n1
            zs    = 0.0_ireals
            n1m1  = n1 - 1

            IF (n1m1 > 0_iintegers) THEN
              DO n2 = 1, n1m1
                zdot      = 0.0_ireals
                DO n3 = 1, n2-1
                  zdot = zdot + za(n3,n2) * za(n3,n1)
                ENDDO
                zt        = za(n2,n1) - zdot
                zt        = zt / za(n2,n2)
                za(n2,n1) = zt
                zs        = zs + zt*zt
              ENDDO
            ENDIF

            zs    = za(n1,n1) - zs
            IF (zs <= 0.0_ireals) EXIT   ! the do n1-loop
            za(n1,n1) = SQRT (zs)
            kierr = 0
          ENDDO

          IF (kierr .NE. 0) THEN
            PRINT *,'  Error in subroutine *gen_grd* inverting matrix za*'
            RETURN
          ENDIF

!   Compute the r.h.s. of the equations
          DO 20 jm=1,6
!
!     Take care of the vertices of diamond #1 with only 5 neighbours
          IF (j1 .EQ. igg1s .AND. j2 .EQ. igg2s               &
                            .AND. jm .EQ. 5) GO TO 20
          IF (j1 .EQ. igg1e .AND. j2 .EQ. igg2s               &
                            .AND. jm .EQ. 6) GO TO 20
          IF (j1 .EQ. igg1s .AND. j2 .EQ. igg2e               &
                            .AND. jm .EQ. 4) GO TO 20
          IF (j1 .EQ. igg1e .AND. j2 .EQ. igg2e               &
                            .AND. jm .EQ. 1) GO TO 20
!
          zb(1) = peta(j1,j2,jm)
          zb(2) = pchi(j1,j2,jm)
          zb(3) = peta(j1,j2,jm)**2
          zb(4) = peta(j1,j2,jm)*pchi(j1,j2,jm)
          zb(5) = pchi(j1,j2,jm)**2
!
!     Solve linear system
          ! The following code replaces the call to linpack routine sposl.
          ! za is an upper triangular matrix which has been computed above
          ! according to linpack routine spofa.

          ! Solve system Transpose (za(new)) * y = zb
          DO n1 = 1,5
            zdot = 0.0_ireals
            DO n3 = 1, n1-1
              zdot = zdot + za(n3,n1) * zb(n3)
            ENDDO
            zb(n1) = (zb(n1) - zdot) / za(n1,n1)
          ENDDO

          ! Solve system za(new) * x = y
          DO n1 = 5, 1, -1
            zb(n1) =  zb(n1) / za(n1,n1)
            zt     = -zb(n1)
            DO n3 = 1, n1-1
              zb(n3) = zt * za(n3,n1) + zb(n3)
            ENDDO
          ENDDO

!     The gradient operator at the 6 surrounding nodes
          pgrd (j1,j2,jm+1,1) = zb(1)/r_earth
          pgrd (j1,j2,jm+1,2) = zb(2)/r_earth
!
!     The Laplacian operator at the 6 surrounding nodes
          prlap(j1,j2,jm+1  ) = 2.*(zb(3) + zb(5))/(r_earth**2)
!
   20     CONTINUE
!

    6   CONTINUE
    5 CONTINUE
!
!     Compute the gradient and Laplacian terms associated with the 
!     'home' node.
      DO jm=1,6
        DO j2 = kig2sm1,kig2ep1
          DO j1 = kig1sm1,kig1ep1
          pgrd (j1,j2,1,1) = pgrd (j1,j2,1,1) - pgrd (j1,j2,jm+1,1)
          pgrd (j1,j2,1,2) = pgrd (j1,j2,1,2) - pgrd (j1,j2,jm+1,2)
          prlap(j1,j2,1  ) = prlap(j1,j2,1  ) - prlap(j1,j2,jm+1  )
          ENDDO
        ENDDO
      ENDDO
!
!------------------------------------------------------------------------------

      kierr = 0

!------------------------------------------------------------------------------

END SUBROUTINE gen_grd

!==============================================================================
!==============================================================================
!+ Defines the spokes.
!------------------------------------------------------------------------------

      SUBROUTINE setspokes (kispoke, kispokes, lzdebug, kierr)

!------------------------------------------------------------------------------
!
! Description:
!   *setspokes* defines the spokes, i.e. the offsets of the 6 (5)
!   neighbouring gridpoints relative to the central node for 2-d 
!   array addressing.
!   See also the subroutine *setmrp* for the definition of the four
!   mirrored points.
!
! Method:
!
!------------------------------------------------------------------------------

! Input/output arrays
  INTEGER (KIND=iintegers)   ::       &
           kispoke (12),  & ! Normal spokes 
           kispokes(12,4)   ! Special spokes for the 4 mirrored nodes
                            !    of the extended array

  INTEGER (KIND=iintegers)   ::       &
           kierr            ! Error flag, set to 0, if no error

  LOGICAL lzdebug           ! Debug flag, if .true. print information

! Local variable
  INTEGER (KIND=iintegers)   ::  j   ! Do loop variable

!------------------------------------------------------------------------------
!
!     Set the normal spokes
      kispoke ( 1) =  1
      kispoke ( 2) =  0
      kispoke ( 3) = -1
      kispoke ( 4) = -1
      kispoke ( 5) =  0
      kispoke ( 6) =  1
      kispoke ( 7) =  0
      kispoke ( 8) =  1
      kispoke ( 9) =  1
      kispoke (10) =  0
      kispoke (11) = -1
      kispoke (12) = -1
!
!     Copy the normal spokes into the special ones
      DO j = 1,12
      kispokes(j,1) = kispoke(j)
      kispokes(j,2) = kispoke(j)
      kispokes(j,3) = kispoke(j)
      kispokes(j,4) = kispoke(j)
      ENDDO
!
!     Now change the special spokes at those points where they
!     deviate from the normal ones
      kispokes( 4,1) = -2
      kispokes( 9,1) =  2
      kispokes(10,1) =  2

      kispokes( 1,2) =  2
      kispokes( 7,2) = -2
      kispokes(12,2) = -2

      kispokes( 4,3) = -2
      kispokes( 5,3) = -2
      kispokes( 6,3) = -1
      kispokes(10,3) =  1
      kispokes(11,3) =  0

      kispokes( 2,4) =  1
!
!------------------------------------------------------------------------------
!
      IF (lzdebug) THEN
        WRITE (*,'(a  )') '  Subroutine *setspokes*'
        WRITE (*,'(a,a)') '  j   ispoke  ispokes(1)  ispokes(2)',   &
                                      '  ispokes(3)  ispokes(4)'
        DO j = 1,12
        WRITE (*,'(i3,i8,4i12)') j, kispoke(j), kispokes(j,1),      &
                  kispokes(j,2), kispokes(j,3), kispokes(j,4)
        ENDDO
      ENDIF
!
      kierr = 0
!
!------------------------------------------------------------------------------
!
!     Index "1" corresponds to the point (kig1s  , kig2sm1) which is 
!               identical to (kig1sm1, kig2s  )
!     Index "2" corresponds to the point (kig1e  , kig2ep1) which is 
!               identical to (kig1ep1, kig2e  )
!     Index "3" corresponds to the point (kig1ep1, kig2sm1) which is 
!               identical to (kig1e  , kig2sm1)
!     Index "4" corresponds to the point (kig1sm1, kig2e  ) which is 
!               identical to (kig1sm1, kig2ep1)
!
!------------------------------------------------------------------------------

END SUBROUTINE setspokes

!==============================================================================
!==============================================================================
!+ Defines i1- and i2-indices of the four mirrored points.
!------------------------------------------------------------------------------

SUBROUTINE setmrp (ki1mrp, ki2mrp, kig1s, kig1e, kig2s, kig2e, lzdebug, kierr)

!------------------------------------------------------------------------------
!
! Description:
!   *setmrp* defines the i1- and i2-indices of the four mirrored
!   points of arrays which has been extended by one row/column
!   around the core of the diamond.
!
! Method:
!
!------------------------------------------------------------------------------

! Input/Output
  INTEGER (KIND=iintegers)   ::     &
           ki1mrp (8),    & ! i1-index of the 4 mirrored points
           ki2mrp (8)       ! i2-index of the 4 mirrored points 
!
  INTEGER (KIND=iintegers)   ::     &
           kig1s,         & ! Start index of the core of the diamond
           kig1e,         & ! End   index of the core of the diamond
           kig2s,         & ! Start index of the core of the diamond
           kig2e            ! End   index of the core of the diamond
!
  LOGICAL lzdebug           ! Debug flag, if .true. print information
!
  INTEGER (KIND=iintegers)   ::     &
           kierr            ! Error flag, set to 0, if no error

! Local variables
  INTEGER (KIND=iintegers)   ::  j  ! Loop index

!------------------------------------------------------------------------------
!
!     RJ: First set all indices to a huge negative value,
!     so that (hopefully) every use by a processor which does not
!     own the mirror point will result in an operand range error
!
      ki1mrp = - HUGE (1_iintegers)
      ki2mrp = - HUGE (1_iintegers)

!     Define the four mirrored points (see subroutine *setspokes*)
      IF (kig1s==igg1s .AND. kig2s==igg2s) THEN
        ki1mrp (1) = kig1s
        ki2mrp (1) = kig2s - 1
        ki1mrp (5) = kig1s - 1
        ki2mrp (5) = kig2s
      ENDIF

      IF (kig1e==igg1e .AND. kig2e==igg2e) THEN
        ki1mrp (2) = kig1e
        ki2mrp (2) = kig2e + 1
        ki1mrp (6) = kig1e + 1
        ki2mrp (6) = kig2e
      ENDIF

      IF (kig1e==igg1e .AND. kig2s==igg2s) THEN
        ki1mrp (3) = kig1e + 1
        ki2mrp (3) = kig2s - 1
        ki1mrp (7) = kig1e
        ki2mrp (7) = kig2s - 1
      ENDIF

      IF (kig1s==igg1s .AND. kig2e==igg2e) THEN
        ki1mrp (4) = kig1s - 1
        ki2mrp (4) = kig2e
        ki1mrp (8) = kig1s - 1
        ki2mrp (8) = kig2e + 1
      ENDIF

!------------------------------------------------------------------------------

      IF (lzdebug) THEN
        WRITE (*,'(a)') '  Subroutine *setmrp*'
        WRITE (*,'(a)') '  j   i1mrp     i2mrp'
        DO j = 1,8
        WRITE (*,'(i3,i7,i9)') j, ki1mrp(j), ki2mrp(j)
        ENDDO
      ENDIF

      kierr = 0

!------------------------------------------------------------------------------

END SUBROUTINE setmrp

!==============================================================================
!==============================================================================
!+ Calculates indices and barycentric coordinates
!------------------------------------------------------------------------------

SUBROUTINE pp_horint_index_p ( plon    , plat   , kgrid  ,               &
                               pxnglob , perlon , perlat ,               &
                               peta    , pchi   , pceneta, pcenchi,      &
                               pbary   ,                                 &
                               kig1s   , kig1e  , kig2s  , kig2e  ,      &
                               kigg1s  , kigg1e , kigg2s , kigg2e ,      &
                               kindex  , pbaryll, protang,               &
                               kidxrp  , krpoints)

!------------------------------------------------------------------------------
!
! Description:
!   This routine is much like *pp_horint_index* but it calculates
!   the indices and barycentric coordinates only in the local domain
!   of the processor. The parameter list has also changed
!
! Method:
!
!------------------------------------------------------------------------------

! Input
  INTEGER (KIND=iintegers)   ::          &
           kgrid      ! dimension of arrays "plon" and "plat"
!
  INTEGER (KIND=iintegers)   ::          &
           kig1s  , & ! start first  dimension local  coordinates
           kig1e  , & ! end   first  dimension local  coordinates
           kig2s  , & ! start second dimension local  coordinates
           kig2e  , & ! end   second dimension local  coordinates
           kigg1s , & ! start first  dimension global coordinates
           kigg1e , & ! end   first  dimension global coordinates
           kigg2s , & ! start second dimension global coordinates
           kigg2e     ! end   second dimension global coordinates

  REAL    (KIND=ireals)      ::          &
           plon (kgrid), & ! longitude points of grid in degrees
           plat (kgrid)    ! latitude  points of grid in degrees

  REAL    (KIND=ireals)      ::          &
           pxnglob (kigg1s:kigg1e, kigg2s:kigg2e, 3, 1:10),    & !
           perlon(kig1s-2:kig1e+2, kig2s-2:kig2e+2, 2, 1:10),  & !
           perlat(kig1s-2:kig1e+2, kig2s-2:kig2e+2, 3, 1:10)     !

  REAL    (KIND=ireals)      ::          &
           peta (kig1s-1:kig1e+1, kig2s-1:kig2e+1, 7),         &
             ! eta-coordinates of the surrounding gridpoints
           pchi (kig1s-1:kig1e+1, kig2s-1:kig2e+1, 7),         &
             ! chi-coordinates of the surrounding gridpoints
           pceneta (kig1s  :kig1e  , kig2s  :kig2e  , 6),      &
             ! eta-coordinates of the surrounding triangle centres
           pcenchi (kig1s  :kig1e  , kig2s  :kig2e  , 6),      &
             ! chi-coordinates of the surrounding triangle centres
           pbary(kig1s  :kig1e  , kig2s  :kig2e  , 6)
             ! factor needed to calculate the barycentric coordinates

!------------------------------------------------------------------------------

! Output
  INTEGER (KIND=iintegers)   ::     kindex(kgrid, 4)  
     ! index array for interpolation to point (kgrid)
     ! kindex(kgrid, 1) : j1 (first  dim.) of nearest GME grid point
     ! kindex(kgrid, 2) : j2 (second dim.) of nearest GME grid point
     ! kindex(kgrid, 3) : jd (diamond) of  nearest GME grid point
     ! kindex(kgrid, 4) : index of triangle containing point (kgrid)

  REAL    (KIND=ireals)      ::     pbaryll(kgrid, 2)
     ! barycentric coordinates alpha and beta to interpolate to point
     ! (kgrid)

  REAL    (KIND=ireals)      ::     protang(kgrid, 2)
     ! sine and cosine of psi, the angle for wind vector rotation

  INTEGER (KIND=iintegers)   ::     kidxrp(*)
     ! index of the points in the grid, which are owned
     ! by the local processor. The grid is thought as a
     ! linear array according to the FORTRAN rules, i.e. point (kgrid)
     ! gives the index kgrid

  INTEGER (KIND=iintegers)   ::     krpoints
     ! number of points in the grid, which are owned
     ! by the local processor.

!------------------------------------------------------------------------------
!
!     Local variables 
!
  REAL    (KIND=ireals)      ::     petall, pchill
     ! local coordinates of point (kgrid)

  REAL    (KIND=ireals)      ::     perlonll(2)
     ! local longitude vector of point (kgrid)

  INTEGER (KIND=iintegers)   ::     mgrid, mni, mnip1, mt
  INTEGER (KIND=iintegers)   ::     j1, j2, jd, jm

  REAL    (KIND=ireals)      ::     z1, z2  

  REAL    (KIND=ireals)      ::       &
           zdet, zdetl, & ! Products of coordinates for spotting
           zdetr          ! the correct triangle

  LOGICAL lfound          ! .true., if correct triangle found

  REAL    (KIND=ireals)      ::       &
           zx, zy, zz, sp, sp_t, zd1, zd2, zd3, zd12, zd13, zd22, zd23, &
           zn, zsn

  LOGICAL lzdebug
!------------------------------------------------------------------------------
!
      lzdebug=.FALSE.
      mni = kigg1e - kigg1s
      mnip1 = mni + 1

      CALL distance (pxnglob, mni, sp, sp_t)

      krpoints = 0

      DO mgrid = 1 , kgrid

            CALL pp_ll2gp (pxnglob, plon(mgrid), plat(mgrid), mnip1,   &
                          zx, zy, zz, sp_t, jd, j1, j2, sp,            &
                          lzdebug .AND. (my_cart_id == 0))
!
!     check if this point is outside the local domain:
!
            IF ( j1<kig1s .OR. j1>kig1e .OR. j2<kig2s .OR. j2>kig2e)   &
            THEN
               kindex (mgrid, 1) = -9999
               kindex (mgrid, 2) = -9999
               kindex (mgrid, 3) = -9999
               kindex (mgrid, 4) = -9999
               pbaryll(mgrid, 1) = 0.
               pbaryll(mgrid, 2) = 0.
               protang(mgrid, 1) = 0.
               protang(mgrid, 2) = 0.
               CYCLE
            ENDIF
!
!     update krpoints, set kidxrp
!
            krpoints = krpoints+1
            kidxrp(krpoints) = mgrid

            kindex(mgrid, 1) = j1
            kindex(mgrid, 2) = j2
            kindex(mgrid, 3) = jd
!
!     local coordinates of point mgrid = zx,zy,zz
!
!!          IF ( ABS(zx) <=  1.4e-14 ) zx = 0.
!!          IF ( ABS(zy) <=  1.4e-14 ) zy = 0.
!!          IF ( ABS(zz) <=  1.4e-14 ) zz = 0.
!
            zd1 = zx * pxnglob(j1,j2,1,jd) +     &
                  zy * pxnglob(j1,j2,2,jd) +     &
                  zz * pxnglob(j1,j2,3,jd)
!    
            zd2 = zx * perlon(j1,j2,1,jd) +     &
                  zy * perlon(j1,j2,2,jd)
!
            zd3 = zx * perlat(j1,j2,1,jd) +     &
                  zy * perlat(j1,j2,2,jd) +     &
                  zz * perlat(j1,j2,3,jd)
!
            zd3 = SIGN ( MIN (1.0_ireals, ABS(zd3)),zd3)
!
            petall = ATAN(zd2/(zd1 + 1.e-6_ireals))
            IF(j1 .EQ. kigg1s .AND. j2 .EQ. kigg2s) petall = ASIN(zd2)
            pchill = ASIN(zd3)
!
!------------------------------------------------------------------------------
!
!     Calculate the index of the triangle containing the point
!     "petall", "pchill"
!     For the diamonds 1 to  5, set "zsn" to  1.
!     For the diamonds 6 to 10, set "zsn" to -1.
!
      zsn = 1.
!
      IF (jd.GE.6) zsn = -1.
!
        mt = 1
!
          lfound = .false.
!
          zdet = petall*pchi(j1,j2,mt)*zsn -               &
                 pchill*peta(j1,j2,mt)
!
          IF ((zdet .LT. 0 .AND. jd .LE. 5) .OR.           &
              (zdet .GT. 0 .AND. jd .GT. 5)) THEN
            zdetl = zdet
            zdetr = petall*pchi(j1,j2,mt+1)*zsn -          &
                    pchill*peta(j1,j2,mt+1)
            IF ((zdetr .GE. 0 .AND. jd .LE. 5) .OR.        &
                (zdetr .LE. 0 .AND. jd .GT. 5)) THEN
              lfound = .true.
               kindex(mgrid,4) = mt
            ENDIF
          ELSE
            zdetr = zdet
            mt    = mt - 1
            IF (mt .EQ. 0) mt = 6
            zdetl = petall*pchi(j1,j2,mt)*zsn -            &
                    pchill*peta(j1,j2,mt)
            IF ((zdetl .LE. 0 .AND. jd .LE. 5) .OR.        &
                (zdetl .GE. 0 .AND. jd .GT. 5)) THEN
              lfound = .true.
               kindex(mgrid,4) = mt
            ENDIF
          ENDIF
!
          IF (.NOT. lfound) THEN
            mt = mt + 1
            IF (mt .EQ. 7) mt = 1
            DO jm = 1,5
              zdetl = zdetr
              mt = mt + 1
              IF (mt .EQ. 7) mt = 1
              zdetr = petall*pchi(j1,j2,mt)*zsn -          &
                      pchill*peta(j1,j2,mt)
              IF ((zdetl .LT. 0. .AND. zdetr .GE. 0. .AND.      &
                     jd .LE. 5) .OR.                            &
                  (zdetl .GT. 0. .AND. zdetr .LE. 0. .AND.      &
                     jd .GT. 5)) THEN
                mt = mt - 1
                IF (mt .EQ. 0) mt = 6
                lfound = .true.
                 kindex(mgrid,4) = mt
                EXIT
              ENDIF
            ENDDO
          ENDIF
!
!------------------------------------------------------------------------------
!
!     Compute the barycentric coordinates 
!
       mt = kindex(mgrid,4)
       pbaryll(mgrid,1) = (    petall*pchi(j1,j2,mt+1) -        &
                           zsn*pchill*peta(j1,j2,mt+1))*        &
                           pbary(j1,j2,mt)
!
       pbaryll(mgrid,2) = (zsn*pchill*peta(j1,j2,mt) -          &
                               petall*pchi(j1,j2,mt))*          &
                           pbary(j1,j2,mt)
!
!------------------------------------------------------------------------------
!
!     Compute the rotation angles "pspsi" and "pcpsi"
!
!C       z1    = -zy                ! y-coordinate of node vector
!C       z2    =  zx                ! x-coordinate of node vector
!C       zn    = 1./SQRT (z1**2 + z2**2 + 1.e-30)
         zn    = 1./SQRT (zy**2 + zx**2 + 1.e-30)
!
         perlonll(1) = -zy*zn ! Normalize to unit length
         perlonll(2) =  zx*zn ! Normalize to unit length
!                                     z-coordinate of perlonll is 0
!
!     Take care of the polepoints
      IF(j1 == kigg1s .AND. j2 == kigg2s) THEN
          z2 = 2.*asin(1.)
          z1 = plon(mgrid) * z2/180.
          perlonll(1) = - SIN(z1)
          perlonll(2) =   COS(z1)
      ENDIF
!
!        Scalar product of the location vector and the local longitude
!        vector at the surrounding gridpoints
          zd12 = zx*perlon(j1,j2,1,jd) + zy*perlon(j1,j2,2,jd)
!
!        Scalar product of the location vector and the local latitude
!        vector at the surrounding gridpoints
          zd13 = zx*perlat(j1,j2,1,jd) + zy*perlat(j1,j2,2,jd) +     &
                 zz*perlat(j1,j2,3,jd)
!
!        Scalar product of the local longitude vectors at the interpolation
!        point and at the surrounding gridpoints
          zd22 = perlonll (1)*perlon(j1,j2,1,jd) +                   &
                 perlonll (2)*perlon(j1,j2,2,jd)
!
!        Scalar product of the local longitude vector at the interpolation
!        point and the latitude vector at the surrounding gridpoints
          zd23 = perlonll (1)*perlat(j1,j2,1,jd) +                   &
                 perlonll (2)*perlat(j1,j2,2,jd)
!
!        Rotation angles "pspsi" and "pcpsi"
          protang(mgrid,1) = COS (petall)*zd22 - SIN (petall)*zd12
          protang(mgrid,2) = COS (petall)*zd23 - SIN (petall)*zd13
!
!
      ENDDO       ! loop over grid

!------------------------------------------------------------------------------

END SUBROUTINE pp_horint_index_p

!==============================================================================

END MODULE src_grids
