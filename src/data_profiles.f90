!+ Variables for GME and LM/HM profiles.
!==============================================================================

MODULE data_profiles

!==============================================================================
!
! Description:
!  This module declares all variables needed for storing profiles of
!  GME and LM/HM. All variables are declared as allocatable arrays. 
!  They are allocated in the routine memprofiles.
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
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler
!  Added data structures for checking isolated points
!
! Code Description:
! Language:           Fortran 90.
! Software Standards: "European Standards for Writing and
!                      Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters, ONLY : &
    ireals,        & ! KIND-type parameters for real variables
    iintegers        ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

  INTEGER (KIND=iintegers), PARAMETER ::                        &
    nmaxgp=100      ! maximal number of grid points for grid point output

  LOGICAL                          ::           &
    lprgp           ! logical for print at selected grid points
  
  INTEGER (KIND=iintegers)              :: &
    ngp_tot,      & ! number of selected grid points in the total domain
    ngp             ! number of selected grid points in this local domain
  
! Global Arrays:

  INTEGER (KIND=iintegers)              :: &
    igp    (nmaxgp), & ! i-indices    local domain
    jgp    (nmaxgp), & ! j-indices    local domain
    igp_tot(nmaxgp), & ! i-indices    total domain
    jgp_tot(nmaxgp)    ! j-indices    total domain
  
  ! Variables for GME profiles (horizontaly interpolated to LM/HM grid)
  !
  REAL (KIND=ireals), ALLOCATABLE :: &
    fisglpr  (:), & ! orography * G                                  (m2/s2)
    tsglpr   (:), & ! temperature of the ground surface              (Kelvin)
    psglpr   (:), & ! surface pressure                               ( Pa  )
    rh_sglpr (:), & ! relative humidity at the surface                 --
    uglpr  (:,:), & ! zonal wind speed                               ( m/s )
    vglpr  (:,:), & ! meridional wind speed                          ( m/s )
    wglpr  (:,:), & ! vertical wind speed (defined on half levels)   ( m/s )
    tglpr  (:,:), & ! temperature                                    (Kelvin)
    qvglpr (:,:), & ! specific water vapor content                   (kg/kg)
    qlglpr (:,:), & ! specific cloud water content                   (kg/kg)
    grhglpr(:,:)    ! generalized relative humidity                    --
  
  ! Variables for GME horizontaly AND verticaly interpolated profiles
  !
  REAL (KIND=ireals), ALLOCATABLE :: &
    upr    (:,:), & ! zonal wind speed                               ( m/s )
    vpr    (:,:), & ! meridional wind speed                          ( m/s )
    tpr    (:,:), & ! temperature                                    (Kelvin)
    qvpr   (:,:), & ! specific water vapor content                   (kg/kg)
    qlpr   (:,:), & ! specific cloud water content                   (kg/kg)
    grhpr  (:,:)    ! generalized relative humidity                    --
  
  ! Variables for LM/HM profiles (from 2D fields)
  !
  REAL (KIND=ireals), ALLOCATABLE :: &
    fispr   (:), & ! orography * G                                   (m2/s2)
    frlapr  (:), & ! part of the land                                  --
    pspr    (:), & ! final surface pressure                          ( Pa  )
    ps1pr   (:), & ! surface pressure after first correction         ( Pa  )
    rh_spr  (:), & ! relative humidity at the surface                  --
    latpr   (:), & ! geographical latitude                           ( rad )
    lonpr   (:), & ! geographical longitude                          ( rad )
    wsnowpr (:), & ! water content of snow                           (m H2O)
    wg1pr   (:), & ! water content of the upper soil layer           (m H2O)
    wg2pr   (:), & ! water content of the medium soil layer          (m H2O)
    wg3pr   (:), & ! water content of the lower soil layer           (m H2O)
                      ! (if nl_soil_lm = 3, unused otherwise)
    wclpr   (:), & ! climatological water content                    (m H2O) 
    tsnowpr (:), & ! temperature of the snow-surface                 (Kelvin)
    tspr    (:), & ! temperature of the ground surface               (Kelvin)
    tmpr    (:), & ! temperature between upper and medium soil layer (Kelvin)
    tclpr   (:), & ! temperature between medium and lower soil layer (Kelvin)
    hmo3pr  (:), & ! ozone maximum                                   ( Pa  )
    vio3pr  (:), & ! vertical integrated ozone contents              (Pa O3)
    udpr    (:), & ! u from divergent wind potential correction      ( m/s )
    vdpr    (:), & ! v from divergent wind potential correction      ( m/s )
    plcpr   (:), & ! part of plant cover                               --
    laipr   (:), & ! leaf area index of plants                         --
    rootpr  (:)    ! depth of the roots                              ( m  )
    
  ! Variables for LM/HM profiles (from 3D fields)
  !
  REAL (KIND=ireals), ALLOCATABLE :: &
    ulmpr  (:,:), & ! zonal wind speed                               ( m/s )
    vlmpr  (:,:), & ! meridional wind speed                          ( m/s )
    wlmpr  (:,:), & ! vertical wind speed (defined on half levels)   ( m/s )
    tlmpr  (:,:), & ! temperature                                    (Kelvin)
    qvlmpr (:,:), & ! specific water vapor content                   (kg/kg)
    qllmpr (:,:), & ! specific cloud water content                   (kg/kg)
    pplmpr (:,:), & ! deviation from the reference pressure          ( Pa  )
    grhlmpr(:,:)    ! generalized relative humidity                    --

!------------------------------------------------------------------------------
! Data structures for isolated points
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers)  ::        &
    niso_loc,       & ! number of isolated points in this subdomain
    niso_max,       & ! maximal number of isolated points in one subdomain
    niso_tot          ! total number of isolated points in the whole COSMO domain

  TYPE struct_for_isolated
    ! indices of isolated point in the COSMO domain
    INTEGER(KIND=iintegers) ::  ic_glob  ! i-index in global COSMO domain
    INTEGER(KIND=iintegers) ::  jc_glob  ! j-index in global COSMO domain
    INTEGER(KIND=iintegers) ::  nc_task  ! ID of COSMO subdomain with this point
    INTEGER(KIND=iintegers) ::  ic_locl  ! i-index in local COSMO subdomain
    INTEGER(KIND=iintegers) ::  jc_locl  ! j-index in local COSMO subdomain
    LOGICAL                 ::  lland    ! whether land point or not

    ! indices of the corresponding point in the coarse grid domain
    INTEGER(KIND=iintegers) ::  iin_glob ! i-index in global coarse grid domain
    INTEGER(KIND=iintegers) ::  jin_glob ! j-index in global coarse grid domain
    INTEGER(KIND=iintegers) ::  nin_task ! ID of subdomain with this point
    INTEGER(KIND=iintegers) ::  iin_locl ! i-index in local coarse grid subdomain
    INTEGER(KIND=iintegers) ::  jin_locl ! j-index in local coarse grid subdomain

    ! Value of the field under consideration (has to be specified before call
    !  to interp_l)
    REAL(KIND=ireals)       ::  value
  END TYPE struct_for_isolated

  ! data structure for isolated points in the whole global COSMO domain
  TYPE (struct_for_isolated), ALLOCATABLE ::    &
    globl_iso_points(:), & ! global structure
    local_iso_points(:)    ! local  structure

!
!- End of module header
!==============================================================================

END MODULE data_profiles
