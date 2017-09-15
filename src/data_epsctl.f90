!+ Variables necessary for ensemble mode.
!==============================================================================

MODULE data_epsctl

!==============================================================================
!
! Description:
!  This module defines the variables necessary for the ensemble mode:
!  The ensemble member ID, the ensemble type (data bas parameters) and the
!  total number of ensemble members.
!
! Method:
!
! Current Code Owner: DWD, FE15, Christoph Gebhardt
!  phone:  +49  69  8062 2689
!  email:  christoph.gebhardt@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V1_7         2007/11/26 Christoph Gebhardt
!  Initial release for INT2LM
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
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

! Global Scalars:

! ensemble definitions
  INTEGER (KIND=iintegers)    ::                                   &
    iepsmem_bc,    & ! ID of the member in the ensemble of boundary
                     ! conditions (iepsmem_bc >= 0)
    iepstyp_bc,    & ! ID of the boundary ensemble generation type
                     ! (iepstyp_bc >= 0)
    iepstot_bc       ! total number of boundary ensemble members (iepstot_bc >= 0)

  LOGICAL                     ::                                   &
    lchk_bc_typ      ! if .TRUE., check member ID of input data (if leps_bc=.TRUE.)

!==============================================================================

END MODULE data_epsctl
