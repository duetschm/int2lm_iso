!+ Variables necessary for the LM grid.
!==============================================================================

MODULE data_grid_lm

!==============================================================================
!
! Description:
!  This module defines the variables necessary for the LM grid,
!  the variables defining the size of the LM-grid 
!  and the values for the reference atmosphere.
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
!  Added polgam; czmls_lm, czhls_lm, ke_soil_lm, cvw_so_lm
!  Eliminated akhlm, bkhlm
! V1_9         2009/09/03 Guenther Zaengl
!  Additional parameters for new reference atmosphere of output model
!  irefatm, delta_t, h_scal
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_20        2012/09/03 Michael Baldauf
!  Added Namelist switch lanalyt_calc_t0p0: indicates that reference p0, t0
!   shall be computed analytically on full levels
! V1_22        2013/07/11 Ulrich Schaettler
!  Put all variables regarding the reference atmosphere  and the vertical
!  coordinate variables to new module vgrid_refatm_utils
! V1_23        2013/10/02 Ulrich Schaettler
!  Remove definition of khmax (is in vgrid_refatm_utils now)
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

! Rotation and position of the grid on the earth
  REAL (KIND=ireals)          ::                                   &
    pollat,      & ! latitude of the rotated north pole (in degrees, N>0)
    pollon,      & ! longitude of the rotated north pole (in degrees, E>0)
    polgam,      & ! latitude of the rotated north pole
    dlat,        & ! grid point distance in zonal direction (in degrees)
    dlon,        & ! grid point distance in meridional direction (in degrees)
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
    endlat,      & ! transformed latitude of the upper right grid point
                   ! of the local domain (in degrees, N>0)
    endlon         ! transformed longitude of the upper right grid point
                   ! of the local domain (in degrees, E>0)

! horizontal size of the grid
  INTEGER (KIND=iintegers)    ::                                   &
    ielm_tot,    & ! ie for LM, total domain
    jelm_tot,    & ! je for LM, total domain
    kelm_tot,    & ! ke for LM
    ielm,        & ! ie for LM, local domain
    jelm,        & ! je for LM, local domain
    kelm,        & ! ke for LM
    ie2lm_max,   & ! Max. of ie2lm on all processors
    je2lm_max,   & ! Max. of je2lm on all processors
    ie2lm_tot,   & ! = ielm_tot + 2
    je2lm_tot,   & ! = jelm_tot + 2
    ij2lm,       & ! = je2lm_tot*je2lm_tot, will be removed soon
    ke1lm,       & !
    kedim,       & !
    ie2lm,       & !
    je2lm

! Start- and end-indices of the computations
  INTEGER (KIND=iintegers)    ::                                   &
    istartpar,   & ! start index for computations in the parallel program
    iendpar,     & ! end index for computations in the parallel program
    jstartpar,   & ! start index for computations in the parallel program
    jendpar        ! end index for computations in the parallel program

! horizontal size of the grid
  INTEGER (KIND=iintegers)    ::                                   &
    ie_ext,      & ! west-east size of fields with external parameters
    je_ext         ! north-south size of fields with external parameters
                   ! these values may be larger than the actual sizes

! Variables for specifying multi-layer soil model
  REAL      (KIND=ireals),    ALLOCATABLE ::       &
    czmls_lm(:),     & ! depth of the LM soil layers in meters
    czhls_lm(:)        ! depth of the LM half soil layers in meters

  INTEGER   (KIND=iintegers)              ::       &
    ke_soil_lm         ! number of output levels in multi-layer soil model

  REAL     (KIND=ireals), ALLOCATABLE     ::       &
    cw_so_rel_lm(:)    ! artifical volumetric soil water content (0-1) profile

!==============================================================================

END MODULE data_grid_lm
