!+ Physical and parametrization constants
!==============================================================================

MODULE data_int2lm_constants

!==============================================================================
!
! Description:
!  This module contains the definitions of physical and parameterization
!  constants.
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
! V1_9         2009/09/03 Guy deMorsier
!  New soil properties of IFS (variables with *_ec; old ones in *_ec_1s)
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:

USE data_parameters, ONLY:   ireals

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global Parameters:
! Physical and parametrization constants

  REAL (KIND=ireals)            ::           &
    t0,      & ! 0 degree Celsius                              [Kelvin]
    r_d,     & ! gas constant for dry air                      [J/K*kg]
    r_v,     & ! gas constant for water vapor                  [J/K*kg]
    rdv,     & ! = R_d/R_v, 
    rvd_m_o, & ! = R_v/R_d - 1.0,  
    o_m_rdv, & ! = 1. - Rdv
    cp_d,    & ! specific heat of dry air at constant pressure [J/K*kg]
    lh_v,    & ! Latent heat of vaporization                   [J/kg]
    lh_f,    & ! Latent heat of fusion                         [J/kg]
    lh_s,    & ! Latent heat of sublimation                    [J/kg]
    g,       & ! gravity at sea level                          [ms-2]
    r_earth, & ! mean radius of the earth                      [m]
    Day_len    ! sidereal day (Sterntag)                       [s]
           

! For the Magnus formula from JAM, 13, 1974, pp. 607
! ln e = ln a + (bT - c) / (T-d)

  REAL (KIND=ireals)            ::           &
    b1,      & !  a
    b2_w,    & !  b
    b2_i,    & !
    b3,      & !  c/b (0 degree Celsius [Kelvin])
    b4_i,    & !
    b4_w       !  d 

  REAL (KIND=ireals)            ::           &
    pi,      & ! circle constant
    raddeg,  & !
    degrad,  & !
    Pid5,    & !
    omcor      !

  REAL (KIND=ireals)            ::           &
    wimax      ! used in himbla to minimize wi

!
! Global Arrays:
  REAL (KIND=ireals)            ::           &
    porb(9),    & ! pore volume of COSMO soil types
    porb_ec_1s, & ! pore volume of IFS 1 soil type
    adpb(9),    & ! air dryness point of COSMO soil types
    pwpb(9),    & ! permanent wilting point of COSMO soil types
    pwpb_ec(6), & ! permanent wilting point of IFS soil types with CY32r3
    pwpb_ec_1s, & ! permanent wilting point of IFS 1 soil type
    fcb (9),    & ! field capacity of COSMO soil types
    fcb_ec(6),  & ! field capacity of IFS soil types with CY32r3
    fcb_ec_1s     ! field capacity of IFS 1 soil type

!- End of module header
!==============================================================================

END MODULE data_int2lm_constants
