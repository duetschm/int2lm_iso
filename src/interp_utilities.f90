!+ Source Module for interpolating from a regular grid
!==============================================================================

MODULE interp_utilities

!==============================================================================
!
! Description:
!   This module contains subroutines that perform linear and quadratic
!   interpolation of a field from a regular grid
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
! V1_5         2007/07/09 Ulrich Schaettler, Uwe Boehm
!  Added Cressman scheme and bicubic interpolation (Uwe Boehm)
!  Editorial changes
! V1_6         2007/09/07 Ulrich Schaettler
!  Bug correction for setting local l_cressman switch in interp_l
! V1_10        2009/12/17 Ulrich Schaettler
!  Added a debug comment for match interpolation
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler
!  Extensions for a domain decomposition independent match interpolation
! V1_15        2010/12/10 Ulrich Blahak
!  Changed INT to NINT for a better handling of rounding
! V1_19        2012/06/06 Davide Cesari, Ulrich Schaettler
!  switch off match interpolation when working with frames
!  If there is no corresponding coarse grid point for a COSMO water point,
!    then we take the linear interpolation from the surrounding coarse grid points
! V1_20        2012/09/03 Ulrich Schaettler
!  Removed obsolete Fortran features
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters , ONLY :  &
    ireals,    &! KIND-type parameters for real variables
    iintegers   ! KIND-type parameter for "normal" integer variables

USE data_int2lm_parallel   , ONLY :  &
    nprocx,          & ! number of processors in x-direction for LM
    nprocy,          & ! number of processors in y-direction for LM
    nproc,           & ! total number of processors: nprocx * nprocy
    num_compute,     & ! number of compute PEs
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    icomm_cart,          & ! communicator for the virtual cartesian topology
    imp_reals,           & ! determines the correct REAL type used in the model
                           ! for MPI
    imp_integers           ! determines the correct INTEGER type used in the
                           ! model for MPI

USE data_profiles,       ONLY :  &
    niso_loc,              & ! number of isolated points in this subdomain
    niso_max,              & ! maximal number of isolated points in one subdomain
    niso_tot,              & ! total number of isolated points in the whole COSMO domain
    struct_for_isolated,   & ! data structure to keep isolated points
    globl_iso_points,      & ! global structure
    local_iso_points         ! local  structure

USE parallel_utilities , ONLY :  &
    gather_values,   & !
    global_values      !


!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE interp_l (px, ie_in, je_in, i_index, j_index,                   &
                     lmono, lposdef, yitype, lframe,                       &
                     lolp_in ,lolp_lm, lmask_lm, undef, x_wght, y_wght,    &
                     pxi    , kilons , kilone , kilats , kilate ,          &
                     startlat_in, startlon_in, latitudes_in, longitudes_in,&
                     lat_coarse,  lon_coarse,                              &
                     grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,       &
                     yerrmsg, kierr )

!------------------------------------------------------------------------------
!
! Description:
!   *interp_l* interpolates the scalar field ("px") bi-linearly to the
!   points with the weights x_wght/y_wght using the 
!   values at the four neighbours.
!   The arrays i_index and j_index give the coordinates of the coarse
!   grid point immediately to the SW of each LM grid point
!
! Method:
!
!==============================================================================
!
! Input
INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  ie_in,    & ! first  dimension of array "px",   end   index
  je_in       ! second dimension of array "px",   end   index

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kilons ,  & ! first  dimension of array "pxi",  start index
  kilone ,  & ! first  dimension of array "pxi",  end   index
  kilats ,  & ! second dimension of array "pxi",  start index
  kilate ,  & ! second dimension of array "pxi",  end   index
  grdpt_rel_in,&! indicator for relations between the coarse and LM grid
  ie_in_tot,& ! ie for input grid, total domain
  je_in_tot   ! je for input grid, total domain

REAL    (KIND=ireals)   , INTENT(IN)   ::    &
  px   (ie_in, je_in),      &
              ! scalar field to be interpolated
  x_wght  (kilons  :kilone  , kilats  :kilate),     &
  y_wght  (kilons  :kilone  , kilats  :kilate),     &
              ! interpolation weights of points
  undef       ! value assigned to undefined grid points

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  i_index  (kilons  :kilone  , kilats  :kilate ), &
  j_index  (kilons  :kilone  , kilats  :kilate )
              ! index of the gridpoints SW of each LM gridpoint

LOGICAL, INTENT(IN)                    ::    &
  lmono, lposdef, &  ! Switch for monotonicity and positive definiteness
  lframe             ! compute boundary fields only on frame

REAL    (KIND=ireals)   , INTENT(IN)   ::    &
  startlat_in,    &  ! transformed latitude of the lower left grid point
  startlon_in,    &  ! transformed latitude of the lower left grid point
  lat_coarse(kilons:kilone, kilats:kilate ), &
  lon_coarse(kilons:kilone, kilats:kilate ), &
  latitudes_in(je_in_tot),   &  ! latitudes of the input data
  longitudes_in(ie_in_tot)      ! longitudes of the input data

LOGICAL, INTENT(IN)                    ::    &
  l_cressman,                & ! to run the Cressman scheme
  lolp_in  (ie_in, je_in),   & !
  lolp_lm  (kilons:kilone, kilats:kilate),   & !
                  ! Land Sea Mask needed for 'M'atch Interp.
  lmask_lm (kilons:kilone, kilats:kilate)      ! mask of points on the frame

CHARACTER (LEN=1), INTENT(IN) :: yitype  ! Interpolation type (L, N, M)
CHARACTER (LEN=*), INTENT(OUT) :: yerrmsg   ! for error message

!==============================================================================

! Output
REAL    (KIND=ireals),    INTENT(OUT)   ::    &
  pxi (kilons:kilone, kilats:kilate)   !  Interpolated scalar field
                !  If running on one processor and *kindex* is fully set,
                !  pxi may be defined as pxi(kilons:kilone, kilats:kilate)
                !  in the calling program
                !  scu: original DWD was pxi(*)

INTEGER (KIND=iintegers), INTENT(OUT)   ::    &
  kierr     ! error flag, kierr = 0 if no error occured

!==============================================================================

! Local variables
REAL (KIND=ireals)         ::      &
  zx0, zx1, zy0, zy1, & ! interpolation weights
  zspx,               & ! sum of surrounding values
  avg                   ! field average computed over valid points

REAL (KIND=ireals), ALLOCATABLE          ::      &
  p_cr(:), dist(:), wgt(:)  ! help arrays for Cressman sceme

INTEGER (KIND=iintegers)   ::      &
  l1, l2,             & ! Loop indices
  j1, j2,  j1u, j2u,  & ! Gridpoint indices
  j1n, j2n,           & ! Nearest gridpoint indices
  j1dn, j2dn,         & ! Relative indices of opposite-to-nearest gridpoint
  mindist, ndist,     & !
  i,j,k1,k2,          & ! Loop indices
  err, istat,         & ! Local and global error condition  
  navg                  ! number of valid points to compute field average

INTEGER (KIND = iintegers) ::      &
  n, n1, ipx ! number of surrounding values

INTEGER (KIND=iintegers), ALLOCATABLE ::      &
  cp(:), cpn(:)         ! vectors for storing neighbour indices

INTEGER (KIND=iintegers)              ::    &
  i_start_glob_in, j_start_glob_in, izerror

LOGICAL                    ::      &
  lngp,            & ! take nearest gridpoint
  lmatch_ngp,      & ! take mean over nearest grid points in case of match 
                     ! interpolation
  lund,            & ! undefined are present
  lval(2) =        & ! map integer to logical
  (/.FALSE.,.TRUE./)

LOGICAL                    ::      &
  l_switch_cr        ! for local use of Cressman scheme

CHARACTER (LEN=80)         ::  &
  yzerror    ! error message for error handling

REAL (KIND=ireals), ALLOCATABLE :: &
  rbufsend(:), rbufrecv(:,:)

REAL (KIND=ireals) :: &
   zeps = 1.E-05_ireals

!=======================================================================

  kierr   = 0
  yerrmsg = '    '

!=======================================================================

  istat      = 0
  izerror    = 0
  lngp       = .false.
  lmatch_ngp = .false.
  IF ( yitype == 'N') lngp = .true.

  IF (l_cressman .AND. .NOT.lframe) THEN
    IF ( (yitype == 'M') .AND. (grdpt_rel_in/=1) ) THEN
      ! this field is not on mass grid point and the Cressman scheme cannot
      ! be used (up to now)
      l_switch_cr = .FALSE.

      ! Print a Warning
      IF (my_cart_id == 0) THEN
        PRINT *,   'Warning from Subroutine interp_l !'
        PRINT *,   'Cressman scheme specified for M-type interpolation,'
        PRINT *,   'but field is not on mass grid point!'
        PRINT *,   'This configuration is not implemented so far'
        PRINT *,   'Cressman scheme is not used !'
      ENDIF
    ELSE
      l_switch_cr = .TRUE.
    ENDIF
  ELSE
    l_switch_cr = .FALSE.
  ENDIF

  IF (((yitype == 'M') .OR. (yitype == 'N')) .AND. .NOT.lframe) THEN
    IF (num_compute > 1) THEN
      ! gather the values of px for isolated points
      ALLOCATE (rbufsend(3*niso_tot), STAT=istat)
      rbufsend (:) = 0.0_ireals
      ALLOCATE (rbufrecv(3*niso_tot,0:num_compute-1), STAT=istat)
      rbufrecv (:,:) = 0.0_ireals

      DO n = 1, niso_tot

        ! If there is no corresponding coarse grid point for a special COSMO grid
        ! point, then nin_task = -1 and nothing will be sent
        IF (globl_iso_points(n)%nin_task == my_cart_id) THEN
          rbufsend(3*n-2) = REAL(globl_iso_points(n)%iin_locl, ireals)
          rbufsend(3*n-1) = REAL(globl_iso_points(n)%jin_locl, ireals)
          rbufsend(3*n  ) = px(globl_iso_points(n)%iin_locl,globl_iso_points(n)%jin_locl)
        ENDIF
      ENDDO

      CALL gather_values(rbufsend, rbufrecv, 3*niso_tot, num_compute, imp_reals,     &
                          -1, icomm_cart, yzerror, izerror)

      DO n1 = 1, niso_tot
        n = globl_iso_points(n1)%nin_task 
        IF (n == -1) THEN
          ! there is no corresponding coarse grid point
          globl_iso_points(n1)%value = undef
        ELSE
          IF ((NINT(rbufrecv(3*n1-2,n), iintegers) == globl_iso_points(n1)%iin_locl) .AND. &
              (NINT(rbufrecv(3*n1-1,n), iintegers) == globl_iso_points(n1)%jin_locl)) THEN
            globl_iso_points(n1)%value = rbufrecv(3*n1,n)
          ENDIF
        ENDIF
      ENDDO
      DEALLOCATE (rbufsend, rbufrecv)
    ELSE
      DO n = 1, niso_tot
        IF (globl_iso_points(n)%nin_task == -1) THEN
          globl_iso_points(n)%value = undef
        ELSE
          globl_iso_points(n)%value = px(globl_iso_points(n)%iin_locl,globl_iso_points(n)%jin_locl)
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  IF ((yitype == 'M') .AND. l_switch_cr) THEN
    ! get the global indices for the local coarse grid coordinate domain
    find_iglob: DO i = 1, ie_in_tot
      IF (ABS(longitudes_in(i) - startlon_in) < zeps) THEN
        i_start_glob_in = i 
        EXIT find_iglob
      ENDIF
    ENDDO find_iglob

    find_jglob: DO j = 1, je_in_tot
      IF (ABS(latitudes_in(j) - startlat_in) < zeps) THEN
        j_start_glob_in = j
        EXIT find_jglob
      ENDIF
    ENDDO find_jglob

    ALLOCATE (p_cr(ie_in*je_in),                      &
              dist(ie_in*je_in),                      &
              wgt(ie_in*je_in) ,  STAT=istat)
  ENDIF

  DO l2 = kilats, kilate
    DO l1 = kilons, kilone

!------------------------------------------------------------------------------

      IF (lframe) THEN
         IF (.NOT. lmask_lm(l1,l2)) CYCLE
      ENDIF

      j1 = i_index(l1,l2) ! first input grid dimension
      j2 = j_index(l1,l2) ! second input grid dimension
      ! Find out which is the nearest input grid point
      IF (x_wght(l1,l2) < 0.5_ireals) THEN
        j1n=j1
      ELSE
        j1n=j1+1
      ENDIF
      IF (y_wght(l1,l2) < 0.5_ireals) THEN
        j2n=j2
      ELSE
        j2n=j2+1
      ENDIF
      ! Compute relative indices of point diagonally opposite to the nearest one
      j1dn=-(2*(j1n-j1)-1)
      j2dn=-(2*(j2n-j2)-1)

!------------------------------------------------------------------------------

      ! Check which interpolation type is required:
      !   'L': normal Linear interpolation -- no action
      !   'N': take Nearest gridpoint (no interpolation) -- see above
      !   'M': Match interpolation,
      !        i.e. linear interpolation if all three surrounding points
      !        are of equal LSM as in LM/HM (lmatch_ngp =.FALSE.)
      !        if not (i.e. lmatch_ngp=.TRUE.): see below par 1 to 3.
      !
      IF (yitype == 'M') THEN
        IF (lframe) THEN
          lmatch_ngp = .false.
        ELSE IF ( ALL (lolp_in(j1:j1+1,j2:j2+1) .EQV. lolp_lm(l1,l2)) ) THEN
          lmatch_ngp = .false.
        ELSE
          lmatch_ngp = .true.
        ENDIF
      ENDIF

!------------------------------------------------------------------------------

      IF (lngp .AND. .NOT.lframe) THEN
        IF ( (lolp_in(j1n,j2n) .EQV. lolp_lm(l1,l2))               .AND.  &
             (px(j1n,j2n) /= undef) ) THEN
          pxi(l1,l2) = px(j1n,j2n)
        ELSEIF ( (lolp_in(j1n+j1dn,j2n) .EQV. lolp_lm(l1,l2))      .AND.  &
                 (px(j1n+j1dn,j2n) /= undef) ) THEN
          pxi(l1,l2) = px(j1n+j1dn,j2n)
        ELSEIF ( (lolp_in(j1n,j2n+j2dn) .EQV. lolp_lm(l1,l2))      .AND.  &
                 (px(j1n,j2n+j2dn) /= undef) ) THEN
          pxi(l1,l2) = px(j1n,j2n+j2dn)
        ELSEIF ( (lolp_in(j1n+j1dn,j2n+j2dn) .EQV. lolp_lm(l1,l2)) .AND.  &
                 (px(j1n+j1dn,j2n+j2dn) /= undef) ) THEN
          pxi(l1,l2) = px(j1n+j1dn,j2n+j2dn)
        ELSE
          pxi(l1,l2) = undef
          kierr = 2
          WRITE (yerrmsg, '(A,I5,A,I5)')                                  &
             'yitype = N: no input point for lat i = ',l2,' and lon j = ',l1

          ! Those who know, what they are doing, can delete the 4 lines above
          ! and use the next two ones instead. The nearest grid point will be
          ! taken then
          ! kierr = 0_iintegers
          ! pxi(l1,l2)=px(j1n,j2n)

        ENDIF

        CYCLE   ! the innermost DO-loop
      ELSE IF (lngp .AND. lframe) THEN

         kierr = 2
          WRITE (yerrmsg, '(A)')                                  &
           'nearest-grid-point interpolation and frame data incompatible'

        CYCLE   ! the innermost DO-loop
      ENDIF

!------------------------------------------------------------------------------

      IF (lmatch_ngp) THEN

        ! 1. If one or more of the 4 surrounding points are of equal
        !    LSM as in LM/HM (i.e. ipx /= 0) --> take mean of these points, 
        !    else, see below par 2. 
        zspx = 0.
        ipx  = 0
        lund  = .FALSE.
        IF (l_switch_cr) THEN
          p_cr = 0.0
          dist = 0.0
          wgt  = 0.0
        ENDIF
        DO k2=j2,j2+1
          DO k1=j1,j1+1
            IF ( (lolp_in(k1,k2) .EQV. lolp_lm(l1,l2))  .AND.    &
                 (px(k1,k2) /= undef) ) THEN
              IF (l_switch_cr) THEN
                p_cr(ipx+1) = px(k1,k2)
                dist(ipx+1) = 1.0 /           &
                   ((longitudes_in(i_start_glob_in+k1-1)-lon_coarse(l1,l2))**2 + &
                   ( latitudes_in(j_start_glob_in+k2-1)-lat_coarse(l1,l2))**2)
              ELSE
                zspx=zspx+px(k1,k2)
              ENDIF
              ipx=ipx+1
            ENDIF
          ENDDO
        ENDDO            
        IF (.NOT. l_switch_cr) THEN
          IF (ipx /= 0) THEN
            pxi(l1,l2)=zspx / REAL(ipx,ireals)
            CYCLE     ! the innermost DO-loop
          ENDIF
        ELSE
          IF (ipx >= 2) THEN
            ! Cressman scheme
            wgt(1:ipx)=   dist(1:ipx)/SUM(dist(1:ipx))
            pxi(l1,l2)=SUM(wgt(1:ipx)*    p_cr(1:ipx))
            CYCLE     ! the innermost DO-loop
          ENDIF
        ENDIF

        ! 2. If one or more of the 8 surrounding points are of equal
        !    LSM as in LM/HM (i.e. ipx /= 0) --> take mean of these points, 
        !    else, see below par 3. 
        zspx = 0.
        ipx  = 0
        lund  = .FALSE.
        IF (l_switch_cr) THEN
          p_cr=0.0
          dist=0.0
          wgt=0.0
        ENDIF
        DO k2=j2n-1,j2n+1
          DO k1=j1n-1,j1n+1
            IF ( (lolp_in(k1,k2) .EQV. lolp_lm(l1,l2))   .AND.    &
                 (px(k1,k2) /= undef) ) THEN
              IF (l_switch_cr) THEN
                p_cr(ipx+1) = px(k1,k2)
                dist(ipx+1) = 1.0/          &
                    ((longitudes_in(i_start_glob_in+k1-1)-lon_coarse(l1,l2))**2 + &
                    ( latitudes_in(j_start_glob_in+k2-1)-lat_coarse(l1,l2))**2)
              ELSE
                zspx=zspx+px(k1,k2)
              ENDIF
              ipx=ipx+1
            ENDIF
          ENDDO
        ENDDO
        IF (.NOT. l_switch_cr) THEN
          IF (ipx /= 0) THEN
            pxi(l1,l2)=zspx / REAL(ipx,ireals)
            CYCLE     ! the innermost DO-loop
          ENDIF
        ELSE
          IF (ipx >= 2) THEN
            wgt(1:ipx)=   dist(1:ipx)/SUM(dist(1:ipx))
            pxi(l1,l2)=SUM(wgt(1:ipx)*    p_cr(1:ipx))
            CYCLE     ! the innermost DO-loop
          ENDIF
        ENDIF

        ! 3. If none of the 8 surrounding points are of equal
        !    LSM as in LM/HM, then look at all 16 input points 
        !    around the 8 surrounding points
        zspx = 0.
        ipx  = 0
        lund  = .FALSE.
        IF (l_switch_cr) THEN
          p_cr=0.0
          dist=0.0
          wgt=0.0
        ENDIF

        DO k2=j2n-2,j2n+2
          DO k1=j1n-2,j1n+2
            IF ( (lolp_in(k1,k2) .EQV. lolp_lm(l1,l2))     .AND.   &
                 (px(k1,k2) /= undef) ) THEN
              IF (l_switch_cr) then
                p_cr(ipx+1) = px(k1,k2)
                dist(ipx+1) = 1.0/        &
                  ((longitudes_in(i_start_glob_in+k1-1)-lon_coarse(l1,l2))**2 + &
                  ( latitudes_in(j_start_glob_in+k2-1)-lat_coarse(l1,l2))**2)
              ELSE
                zspx=zspx+px(k1,k2)
              ENDIF
              ipx=ipx+1
            ENDIF
          ENDDO
        ENDDO
        IF (.NOT. l_switch_cr) THEN
          IF (ipx /= 0) THEN
            pxi(l1,l2)=zspx / REAL(ipx,ireals)
            CYCLE     ! the innermost DO-loop
          ENDIF
        ELSE
          IF (ipx >= 2) THEN
            wgt(1:ipx)=   dist(1:ipx)/SUM(dist(1:ipx))
            pxi(l1,l2)=SUM(wgt(1:ipx)*    p_cr(1:ipx))
            CYCLE     ! the innermost DO-loop
          ENDIF
        ENDIF

        ! 4. If this point is reached, we have to take the value from 
        !    the isolated values data structure

        DO n = 1, niso_tot
          IF (globl_iso_points(n)%nin_task == -1) THEN
            ! just do a linear interpolation of the 4 surrounding coarse grid points
            pxi(l1,l2) = 0.25_ireals * (px(j1,j2)+px(j1+1,j2)+px(j1,j2+1)+px(j1+1,j2+1))
          ELSE
            IF ( (globl_iso_points(n)%ic_locl == l1)    .AND.    &
                 (globl_iso_points(n)%jc_locl == l2)    .AND.    &
                 (globl_iso_points(n)%nc_task == my_cart_id) ) THEN
              pxi(l1,l2) = globl_iso_points(n)%value
            ENDIF
          ENDIF
        ENDDO

        CYCLE     ! the innermost DO-loop

      ENDIF       ! lmatch_ngp

!------------------------------------------------------------------------------

      zx1   = x_wght(l1,l2)
      zx0   = 1.0_ireals - zx1
      zy1   = y_wght(l1,l2)
      zy0   = 1.0_ireals - zy1

      ! Interpolate px bi-linearly
      IF (lframe) THEN
        IF (px(j1,j2  ) == undef .OR. px(j1+1,j2  ) == undef .OR. &
            px(j1,j2+1) == undef .OR. px(j1+1,j2+1) == undef) THEN
          pxi(l1,l2) = undef
          kierr = 2
          WRITE (yerrmsg, '(A,I5,A,I5,A)')                                  &
           'no input point for lat i = ',l2,' and lon j = ',l1,             &
           ' input frame is probably too narrow'
          CYCLE
        ENDIF
      ENDIF

      ! this interpolation is only done, if all surrounding points px are
      ! defined (this is also the case, if the full field is defined, but
      ! lframe=.TRUE.
      pxi(l1,l2) = (zx0*px(j1,j2  )+zx1*px(j1+1,j2  ))*zy0 + &
                   (zx0*px(j1,j2+1)+zx1*px(j1+1,j2+1))*zy1

      ! Enforce monotonicity, if required
      !IF (lmono) THEN ! Mathematically enforced
      !  pxi(l1,l2) = MIN (pxi(l1,l2), MAX (px(j1,j2),px(j1+1,j2), &
      !   px(j1,j2+1),px(j1+1,j2+1)))
      !  pxi(l1,l2) = MAX (pxi(l1,l2), MIN (px(j1,j2),px(j1+1,j2), &
      !   px(j1,j2+1),px(j1+1,j2+1)))
      !ENDIF
    ENDDO
  ENDDO

  IF (lposdef) THEN
    DO l2 = kilats, kilate
      DO l1 = kilons, kilone
        ! Enforce positive definiteness, if required
        pxi(l1,l2) = MAX (pxi(l1,l2), 0.0_ireals)
      ENDDO
    ENDDO
  ENDIF

  IF ( (yitype == 'M') .AND. l_switch_cr) THEN
    DEALLOCATE (p_cr, dist, wgt)
  ENDIF

  IF (lframe) THEN
    navg=COUNT(pxi /= undef .AND. lmask_lm)
    avg=SUM(pxi, MASK=(pxi /= undef .AND. lmask_lm))

    IF (nproc > 1) THEN
      CALL global_values (navg, 1, 'SUM', imp_integers, icomm_cart, -1, &
       yerrmsg, kierr)
      CALL global_values (avg, 1, 'SUM', imp_reals, icomm_cart, -1, &
       yerrmsg, kierr)
    ENDIF
    IF (navg > 0) THEN
      avg = avg/navg
    ELSE
      avg = 0.0_ireals
    ENDIF
    WHERE (.NOT. lmask_lm)
      pxi = avg
    ENDWHERE
  ENDIF

!------------------------------------------------------------------------------

END SUBROUTINE interp_l

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE interp_q (px, ie_in, je_in, i_index, j_index, lmono, lposdef,   &
                     lframe, undef, lmask_lm, x_wght, y_wght,              &
                     pxi, kilons , kilone , kilats , kilate ,              &
                     yerrmsg, kierr )

!------------------------------------------------------------------------------
!
! Description:
!   quadratic interpolation, not yet implemented, for now it's linear
!   *interp_l* interpolates the scalar field ("px") bi-linearly to the
!   points with the weights x_wght/y_wght using the 
!   values at the four neighbours.
!   The arrays i_index and j_index give the coordinates of the coarse
!   grid point immediately to the SW of each LM grid point
!
! Method:
!
!==============================================================================
!
! Input
INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  ie_in,    & ! first  dimension of array "px",   end   index
  je_in       ! second dimension of array "px",   end   index

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kilons ,  & ! first  dimension of array "pxi",  start index
  kilone ,  & ! first  dimension of array "pxi",  end   index
  kilats ,  & ! second dimension of array "pxi",  start index
  kilate      ! second dimension of array "pxi",  end   index

REAL    (KIND=ireals)   , INTENT(IN)   ::    &
  px   (ie_in, je_in),      &
              ! scalar field to be interpolated
  undef,    & ! value assigned to undefined grid points
  x_wght  (kilons  :kilone  , kilats  :kilate),     &
  y_wght  (kilons  :kilone  , kilats  :kilate)
              ! interpolation weights of points

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  i_index  (kilons  :kilone  , kilats  :kilate ), &
  j_index  (kilons  :kilone  , kilats  :kilate )
              ! index of the gridpoints SW of each LM gridpoint

LOGICAL, INTENT(IN)                    ::    &
  lmono, lposdef, &  ! Switch for monotonicity and positive definiteness
  lframe             ! compute boundary fields only on frame

LOGICAL, INTENT(IN)                    ::    &
  lmask_lm (kilons:kilone, kilats:kilate)   ! mask of points on the frame

CHARACTER (LEN=*), INTENT(OUT) :: yerrmsg   ! for error message

!==============================================================================

! Output
REAL    (KIND=ireals),    INTENT(OUT)   ::    &
  pxi (kilons:kilone, kilats:kilate)   !  Interpolated scalar field
                !  If running on one processor and *kindex* is fully set,
                !  pxi may be defined as pxi(kilons:kilone, kilats:kilate)
                !  in the calling program
                !  scu: original DWD was pxi(*)

INTEGER (KIND=iintegers), INTENT(OUT)   ::    &
  kierr     ! error flag, kierr = 0 if no error occured

!==============================================================================

! Local variables
REAL (KIND=ireals)         ::      &
  w, x1, x2, x3, pxj1, pxj2, pxj3, &
  avg,                & ! field average computed over valid points
  aeastnorth,         & ! statement functions for constructing Lagrange
  awestsouth            ! 2^ order polynomials from interpolation indices

INTEGER (KIND=iintegers)   ::      &
  l1, l2,             & ! Loop indices
  j1, j2, j4,         & ! Gridpoint indices
  navg                  ! number of valid points to compute field average

!=======================================================================
! Statement functions for quadratic interpolation
aeastnorth(w,x1,x2,x3) = w*(w-1.)*0.5*x3 - w*(w-2.)*x2 + (1.-1.5*w+0.5*w**2)*x1
!aeastnorth(w,x1,x2,x3) = w*(1.-w)*0.5*x3 + w**2*x2 - (1.-1.5*w+0.5*w**2)*x1
!awestsouth(w,x1,x2,x3) = w*(1.-w)*0.5*x3 + (1.5*w*x2-0.5*w**2)*x2 - (w**2-1.)*x1
awestsouth(w,x1,x2,x3) = w*(w-1.)*0.5*x3 + w*(w+1.)*0.5*x2 - (w**2-1.)*x1

!=======================================================================

  kierr   = 0
  yerrmsg = '    '

!=======================================================================

  DO l2 = kilats, kilate
    DO l1 = kilons, kilone

      IF (lframe) THEN
         IF (.NOT. lmask_lm(l1,l2)) CYCLE
      ENDIF

      j1 = i_index(l1,l2) ! first input grid dimension
!      IF (j1 == -9999) CYCLE
      j2 = j_index(l1,l2) ! second input grid dimension
      ! Find out which is the nearest input grid point

!------------------------------------------------------------------------------

      IF (y_wght(l1,l2) > 0.5_ireals) THEN
        j4 = j2 + 2
      ELSE
        j4 = j2 - 1
      ENDIF
      ! East-west interpolation
      IF (x_wght(l1,l2) > 0.5_ireals) THEN
        IF (lframe) THEN
          IF (px(j1,j2) == undef .OR. px(j1+1,j2) == undef .OR. &
           px(j1+2,j2) == undef .OR. px(j1,j2+1) == undef .OR. &
           px(j1+1,j2+1) == undef .OR. px(j1+2,j2+1) == undef &
           .OR. px(j1,j4) == undef .OR. px(j1+1,j4) == undef .OR. &
           px(j1+2,j4) == undef) THEN
            pxi(l1,l2) = undef
            kierr = 2
            WRITE (yerrmsg, '(A,I5,A,I5)')                                  &
             'yitype = N: no input point for lat i = ',l2,' and lon j = ',l1
            CYCLE
          ENDIF
        ENDIF
        pxj1 = aeastnorth(x_wght(l1,l2), px(j1,j2), px(j1+1,j2), px(j1+2,j2))
        pxj2 = aeastnorth(x_wght(l1,l2), px(j1,j2+1), px(j1+1,j2+1), px(j1+2,j2+1))
        pxj3 = aeastnorth(x_wght(l1,l2), px(j1,j4), px(j1+1,j4), px(j1+2,j4))
      ELSE
        IF (lframe) THEN
          IF (px(j1,j2) == undef .OR. px(j1+1,j2) == undef .OR. &
           px(j1-1,j2) == undef .OR. px(j1,j2+1) == undef .OR. &
           px(j1+1,j2+1) == undef .OR. px(j1-1,j2+1) == undef &
           .OR. px(j1,j4) == undef .OR. px(j1+1,j4) == undef .OR. &
           px(j1+2,j4) == undef) THEN
            pxi(l1,l2) = undef
            kierr = 2
            WRITE (yerrmsg, '(A,I5,A,I5)')                                  &
             'yitype = N: no input point for lat i = ',l2,' and lon j = ',l1
            CYCLE
          ENDIF
        ENDIF
        pxj1 = awestsouth(x_wght(l1,l2), px(j1,j2), px(j1+1,j2), px(j1-1,j2))
        pxj2 = awestsouth(x_wght(l1,l2), px(j1,j2+1), px(j1+1,j2+1), px(j1-1,j2+1))
        pxj3 = awestsouth(x_wght(l1,l2), px(j1,j4), px(j1+1,j4), px(j1-1,j4))
      ENDIF
      ! North-south interpolation
      IF (y_wght(l1,l2) > 0.5_ireals) THEN
        pxi(l1,l2) = aeastnorth(y_wght(l1,l2), pxj1, pxj2, pxj3)
      ELSE
        pxi(l1,l2) = awestsouth(y_wght(l1,l2), pxj1, pxj2, pxj3)
      ENDIF

      ! Enforce monotonicity, if required
      !IF (lmono) THEN ! Mathematically enforced
      !  pxi(l1,l2) = MIN (pxi(l1,l2), MAX (px(j1,j2),px(j1+1,j2), &
      !   px(j1,j2+1),px(j1+1,j2+1)))
      !  pxi(l1,l2) = MAX (pxi(l1,l2), MIN (px(j1,j2),px(j1+1,j2), &
      !   px(j1,j2+1),px(j1+1,j2+1)))
      !ENDIF

      ! Enforce positive definiteness, if required
      IF (lposdef) THEN
        pxi(l1,l2) = MAX (pxi(l1,l2), 0.0_ireals)
      ENDIF

    ENDDO
  ENDDO

  IF (lframe) THEN
    navg=COUNT(pxi /= undef .AND. lmask_lm)
    avg=SUM(pxi, MASK=(pxi /= undef .AND. lmask_lm))

    IF (nproc > 1) THEN
      CALL global_values (navg, 1, 'SUM', imp_integers, icomm_cart, -1, &
       yerrmsg, kierr)
      CALL global_values (avg, 1, 'SUM', imp_reals, icomm_cart, -1, &
       yerrmsg, kierr)
    ENDIF
    IF (navg > 0) THEN
      avg = avg/navg
    ELSE
      avg = 0.0_ireals
    ENDIF
    WHERE (.NOT. lmask_lm)
      pxi = avg
    ENDWHERE
!    PRINT*,'AveragedQ ',navg, avg, COUNT(.NOT. lmask_lm)
  ENDIF

END SUBROUTINE interp_q

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE interp_q_lm                                                   &
       (p_in,  ie_in, je_in,                                             &
        p_out, i_index, j_index, x_wght, y_wght, ie_lm, je_lm,           &
        lmono, lposdef, yerrmsg, kierr )

!------------------------------------------------------------------------------
!
! Description:
!   quadratic interpolation following the 16-point formula from emtodm.
!   The arrays i_index and j_index give the coordinates of the coarse
!   grid point immediately to the SW of each LM grid point;
!   x_wght and y_wght give the relative distance to this grid point
!
! Method:
!
!==============================================================================

! Input
INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  ie_in,    & ! first  dimension of array "p_in",   end   index
  je_in       ! second dimension of array "p_in",   end   index

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  ie_lm  ,  & ! first  dimension of array "p_out",
  je_lm       ! second dimension of array "p_out",

REAL    (KIND=ireals)   , INTENT(IN)   ::    &
  p_in  (ie_in, je_in)       ! scalar field to be interpolated

REAL    (KIND=ireals)   , INTENT(IN)   ::    &
  x_wght(ie_lm, je_lm),   &  ! relative distance to SW-coarse grid point
  y_wght(ie_lm, je_lm)       ! relative distance to SW-coarse grid point

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  i_index (ie_lm, je_lm), &  ! index of SW coarse grid point
  j_index (ie_lm, je_lm)     ! index of SW coarse grid point

LOGICAL, INTENT(IN)                    ::    &
  lmono, lposdef  ! Switch for monotonicity and positive definiteness


! Output
REAL    (KIND=ireals),    INTENT(OUT)   ::    &
  p_out (ie_lm, je_lm)  !  Interpolated scalar field

INTEGER (KIND=iintegers), INTENT(OUT)   ::    &
  kierr     ! error flag, kierr = 0 if no error occured

CHARACTER (LEN=*), INTENT(OUT) :: yerrmsg   ! for error message

!==============================================================================

! Local variables
INTEGER (KIND=iintegers)   ::      &
 i, j, ic, jc          ! Loop indices

REAL (KIND=ireals)         ::              &
  zdx, zdy, zmin, zmax,                    &
  a11, a12, a13, a14, a21, a22, a23, a24,  &
  a31, a32, a33, a34, a41, a42, a43, a44,  &
  t01, t02, t03, t04, t05, t06, t07, t08,  &
  t09, t10, t11, t12, t13, t14, t15, t16

!==============================================================================

  kierr   = 0
  yerrmsg = '    '

  DO j = 1, je_lm
    DO i = 1, ie_lm

      ! Compute the interpolation weights
      ic  = i_index(i,j)   ! first  dimension of coarse grid
      jc  = j_index(i,j)   ! second dimension of coarse grid
      zdx = x_wght (i,j)  
      zdy = y_wght (i,j)  

      ! Get the 16 values of the input field
      t01 = p_in(ic-1,jc-1)
      t02 = p_in(ic-1,jc  )
      t03 = p_in(ic-1,jc+1)
      t04 = p_in(ic-1,jc+2)
      t05 = p_in(ic  ,jc-1)
      t06 = p_in(ic  ,jc  )
      t07 = p_in(ic  ,jc+1)
      t08 = p_in(ic  ,jc+2)
      t09 = p_in(ic+1,jc-1)
      t10 = p_in(ic+1,jc  )
      t11 = p_in(ic+1,jc+1)
      t12 = p_in(ic+1,jc+2)
      t13 = p_in(ic+2,jc-1)
      t14 = p_in(ic+2,jc  )
      t15 = p_in(ic+2,jc+1)
      t16 = p_in(ic+2,jc+2)

      ! Determine the 16 interpolation coefficients
      a11 =                             t06
      a21 =             0.5 *t10            - 0.5 *t02
      a31 = -0.5 *t14 + 2.0 *t10 - 2.5 *t06 +      t02
      a41 =  0.5 *t14 - 1.5 *t10 + 1.5 *t06 - 0.5 *t02
      
      a12 =  0.5 *t07 - 0.5 *t05
      a22 =  0.25*t11 - 0.25*t09 - 0.25*t03 + 0.25*t01
      a32 = -0.25*t15 + 0.25*t13 +      t11 -      t09 -            &
             1.25*t07 + 1.25*t05 + 0.5 *t03 - 0.5 *t01
      a42 =  0.25*t15 - 0.25*t13 - 0.75*t11 + 0.75*t09 +            &
             0.75*t07 - 0.75*t05 - 0.25*t03 + 0.25*t01
      
      a13 =  2.0 *t07 - 2.5 *t06 +      t05 - 0.5 *t08
      a23 =       t11 - 0.25*t12 - 1.25*t10 + 0.5 *t09 +            &
             0.25*t04 -      t03 + 1.25*t02 - 0.5 *t01
      a33 =  0.25*t16 -      t15 + 1.25*t14 - 0.5 *t13 -            &
                  t12 + 4.0 *t11 - 5.0 *t10 + 2.0 *t09 +            &
             1.25*t08 - 5.0 *t07 + 6.25*t06 - 2.5 *t05 -            &
             0.5 *t04 + 2.0 *t03 - 2.5 *t02 +      t01
      a43 = -0.25*t16 +      t15 - 1.25*t14 + 0.5 *t13 +            &
             0.75*t12 - 3.0 *t11 + 3.75*t10 - 1.5 *t09 -            &
             0.75*t08 + 3.0 *t07 - 3.75*t06 + 1.5 *t05 +            &
             0.25*t04 -      t03 + 1.25*t02 - 0.5 *t01
      
      a14 =  0.5 *t08 - 1.5 *t07 + 1.5 *t06 - 0.5 *t05
      a24 =  0.25*t12 - 0.75*t11 + 0.75*t10 - 0.25*t09  -           &
             0.25*t04 + 0.75*t03 - 0.75*t02 + 0.25*t01
      a34 = -0.25*t16 + 0.75*t15 - 0.75*t14 + 0.25*t13 +            &
                  t12 - 3.0 *t11 + 3.0 *t10 -      t09 -            &
             1.25*t08 + 3.75*t07 - 3.75*t06 + 1.25*t05 +            &
             0.5 *t04 - 1.5 *t03 + 1.5 *t02 - 0.5 *t01
      a44 =  0.25*t16 - 0.75*t15 + 0.75*t14 - 0.25*t13 -            &
             0.75*t12 + 2.25*t11 - 2.25*t10 + 0.75*t09 +            &
             0.75*t08 - 2.25*t07 + 2.25*t06 - 0.75*t05 -            &
             0.25*t04 + 0.75*t03 - 0.75*t02 + 0.25*t01

      p_out(i,j) =       (a11 + a21*zdx + a31*zdx**2 + a41*zdx**3) + &
                zdy    * (a12 + a22*zdx + a32*zdx**2 + a42*zdx**3) + &
                zdy**2 * (a13 + a23*zdx + a33*zdx**2 + a43*zdx**3) + &
                zdy**3 * (a14 + a24*zdx + a34*zdx**2 + a44*zdx**3)
    ENDDO
  ENDDO

  ! Enforce monotonicity, if required
  IF (lmono) THEN
    DO j = 1, je_lm
      DO i = 1, ie_lm
        ic  = i_index(i,j)   ! first  dimension of coarse grid
        jc  = j_index(i,j)   ! second dimension of coarse grid
        zmax = MAX (p_in(ic,jc), p_in(ic+1,jc), p_in(ic,jc+1), p_in(ic+1,jc+1))
        zmin = MIN (p_in(ic,jc), p_in(ic+1,jc), p_in(ic,jc+1), p_in(ic+1,jc+1))
        p_out(i,j)   = MIN (zmax, p_out(i,j))
        p_out(i,j)   = MAX (zmin, p_out(i,j))
      ENDDO
    ENDDO
  ENDIF

  ! Enforce positive definiteness, if required
  IF (lposdef) THEN
    DO j = 1, je_lm
      DO i = 1, ie_lm
        p_out(i,j) = MAX (p_out(i,j), 0.0_ireals)
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE interp_q_lm

!==============================================================================

SUBROUTINE interp_q_bs (px, ie_in, je_in, i_index, j_index, lmono, lposdef,      &
                        lframe, lolp_in, lolp_lm, undef, lmask_lm,               &
                        x_wght, y_wght, pxi, kilons , kilone , kilats , kilate , &
                        startlat_in, startlon_in, latitudes_in, longitudes_in,   &
                        lat_coarse, lon_coarse,                                  &
                        grdpt_rel_in, ie_in_tot, je_in_tot, lcm2lm, l_cressman,  &
                        yerrmsg, kierr )

!------------------------------------------------------------------------------
!
! Description:
!   Control interface routine for performing a bicubic spline interpolation.
!   Depending on the number of available coarse-grid point rows and lines
!   surrounding the lm fine grid, either a bilinear interpolation by calling
!   subroutine *interp_l* or the spline interpolation by calling subroutine
!   *interp_q_lm* is performed at the array borders. In the interior,
!   *interp_q_lm* is called.
!
! Method:
!
!==============================================================================
!
! Input
INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  ie_in,    & ! first  dimension of array "px",   end   index
  je_in       ! second dimension of array "px",   end   index

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kilons ,  & ! first  dimension of array "pxi",  start index
  kilone ,  & ! first  dimension of array "pxi",  end   index
  kilats ,  & ! second dimension of array "pxi",  start index
  kilate ,  & ! second dimension of array "pxi",  end   index
  grdpt_rel_in,&! indicator for relations between the coarse and LM grid
  ie_in_tot,& ! ie for input grid, total domain
  je_in_tot   ! je for input grid, total domain

REAL    (KIND=ireals)   , INTENT(IN)   ::    &
  px   (ie_in, je_in),      &
              ! scalar field to be interpolated
  undef,    & ! value assigned to undefined grid points
  x_wght  (kilons  :kilone  , kilats  :kilate),     &
  y_wght  (kilons  :kilone  , kilats  :kilate)
              ! interpolation weights of points

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  i_index  (kilons  :kilone  , kilats  :kilate ), &
  j_index  (kilons  :kilone  , kilats  :kilate )
              ! index of the gridpoints SW of each LM gridpoint

REAL    (KIND=ireals)   , INTENT(IN)   ::    &
  startlat_in,    &  ! transformed latitude of the lower left grid point
  startlon_in,    &  ! transformed latitude of the lower left grid point
  lat_coarse(kilons:kilone, kilats:kilate ), &
  lon_coarse(kilons:kilone, kilats:kilate ), &
  latitudes_in(je_in_tot),   &  ! latitudes of the input data
  longitudes_in(ie_in_tot)      ! longitudes of the input data

LOGICAL, INTENT(IN)                    ::    &
  lmono, lposdef, &  ! Switch for monotonicity and positive definiteness
  lcm2lm,         &  ! if .TRUE., climate model ->lm
  l_cressman,     &  ! to run the Cressman scheme
  lframe             ! compute boundary fields only on frame

! UB
LOGICAL, INTENT(IN)                    ::    &
  lolp_in  (ie_in, je_in),   & !
  lolp_lm  (kilons:kilone, kilats:kilate)   !
                  ! Land Sea Mask needed for 'M'atch Interp.
! UB

LOGICAL, INTENT(IN)                    ::    &
  lmask_lm (kilons:kilone, kilats:kilate)   ! mask of points on the frame

CHARACTER (LEN=*), INTENT(OUT) :: yerrmsg   ! for error message

!==============================================================================

! Output
REAL    (KIND=ireals),    INTENT(OUT)   ::    &
  pxi (kilons:kilone, kilats:kilate)   !  Interpolated scalar field
                !  If running on one processor and *kindex* is fully set,
                !  pxi may be defined as pxi(kilons:kilone, kilats:kilate)
                !  in the calling program
                !  scu: original DWD was pxi(*)

INTEGER (KIND=iintegers), INTENT(OUT)   ::    &
  kierr     ! error flag, kierr = 0 if no error occured

!==============================================================================

! Local variables
REAL (KIND=ireals)         ::      &
!  zx0, zx1, zy0, zy1, & ! interpolation weights
  w, x1, x2, x3, pxj1, pxj2, pxj3, &
  avg,                & ! field average computed over valid points
  aeastnorth,         & ! statement functions for constructing Lagrange
  awestsouth            ! 2^ order polynomials from interpolation indices

! UB

REAL (KIND=ireals) ,               &
      ALLOCATABLE          ::      &
    px_bl (:,:)

REAL    (KIND=ireals)      ::      &
    zminlat,   zmaxlat,   zminlon,   zmaxlon

INTEGER (KIND=iintegers)   ::      &
  l1, l2,             & ! Loop indices
  j_s, j_n,           & ! Gridpoint indices
  i_w, i_e,           & ! in four directions
  j1,   j2,   j4,     & ! Gridpoint indices
  navg                  ! number of valid points to compute field average

INTEGER (KIND=iintegers)   ::      &
    nzbounds_e, nzbounds_w, nzbounds_s, nzbounds_n
                        ! number of boundary lines for coarse grid
! UB

!=======================================================================
! Statement functions for quadratic interpolation
aeastnorth(w,x1,x2,x3) = w*(w-1.)*0.5*x3 - w*(w-2.)*x2 + (1.-1.5*w+0.5*w**2)*x1
!aeastnorth(w,x1,x2,x3) = w*(1.-w)*0.5*x3 + w**2*x2 - (1.-1.5*w+0.5*w**2)*x1
!awestsouth(w,x1,x2,x3) = w*(1.-w)*0.5*x3 + (1.5*w*x2-0.5*w**2)*x2 - (w**2-1.)*x1
awestsouth(w,x1,x2,x3) = w*(w-1.)*0.5*x3 + w*(w+1.)*0.5*x2 - (w**2-1.)*x1

!=======================================================================

  kierr   = 0
  yerrmsg = '    '

!=======================================================================

! UB

  ! Identify the number of coarse grid points surrounding the lm
  ! grid for each processor
  ! (One more grid point at the northern and eastern boundary
  ! than in subroutine decompose_coarse_grid or decompose_cm
  ! since here, the real number of surrounding rows and lines
  ! is of interest).


  nzbounds_s =         MINVAL (j_index(kilons:kilone,kilats:kilate))
  nzbounds_n = je_in - MAXVAL (j_index(kilons:kilone,kilats:kilate)) + 1
  nzbounds_w =         MINVAL (i_index(kilons:kilone,kilats:kilate))
  nzbounds_e = ie_in - MAXVAL (i_index(kilons:kilone,kilats:kilate)) + 1

  IF (lcm2lm) THEN

    ! handle frames
    IF (lframe) THEN

      DO l2 = kilats, kilate
        DO l1 = kilons, kilone

          IF (.NOT. lmask_lm(l1,l2)) CYCLE

          j1 = i_index(l1,l2) ! first input grid dimension
!         IF (j1 == -9999) CYCLE
          j2 = j_index(l1,l2) ! second input grid dimension
          ! Find out which is the nearest input grid point

!------------------------------------------------------------------------------

          IF (y_wght(l1,l2) > 0.5_ireals) THEN
            j4 = j2 + 2
          ELSE
            j4 = j2 - 1
          ENDIF
          ! East-west interpolation
          IF (x_wght(l1,l2) > 0.5_ireals) THEN
            IF (px(j1  ,j2  ) == undef .OR. px(j1+1,j2  ) == undef .OR. &
                px(j1+2,j2  ) == undef .OR. px(j1,  j2+1) == undef .OR. &
                px(j1+1,j2+1) == undef .OR. px(j1+2,j2+1) == undef .OR. &
                px(j1,  j4)   == undef .OR. px(j1+1,j4  ) == undef .OR. &
                px(j1+2,j4  ) == undef) THEN
             pxi(l1,l2) = undef
             kierr = 2
             WRITE (yerrmsg, '(A,I5,A,I5)')                                  &
             'yitype = N: no input point for lat i = ',l2,' and lon j = ',l1
             CYCLE
            ENDIF
            pxj1 = aeastnorth(x_wght(l1,l2), px(j1,j2  ), px(j1+1,j2  ), px(j1+2,j2  ))
            pxj2 = aeastnorth(x_wght(l1,l2), px(j1,j2+1), px(j1+1,j2+1), px(j1+2,j2+1))
            pxj3 = aeastnorth(x_wght(l1,l2), px(j1,j4  ), px(j1+1,j4  ), px(j1+2,j4  ))
          ELSE
            IF (px(j1,j2    ) == undef .OR. px(j1+1,j2  ) == undef .OR. &
                px(j1-1,j2  ) == undef .OR. px(j1,  j2+1) == undef .OR. &
                px(j1+1,j2+1) == undef .OR. px(j1-1,j2+1) == undef .OR. &
                px(j1,  j4)   == undef .OR. px(j1+1,j4  ) == undef .OR. &
                px(j1+2,j4)   == undef) THEN
              pxi(l1,l2) = undef
              kierr = 2
              WRITE (yerrmsg, '(A,I5,A,I5)')                                  &
              'yitype = N: no input point for lat i = ',l2,' and lon j = ',l1
              CYCLE
            ENDIF
            pxj1 = awestsouth(x_wght(l1,l2), px(j1,j2  ), px(j1+1,j2  ), px(j1-1,j2  ))
            pxj2 = awestsouth(x_wght(l1,l2), px(j1,j2+1), px(j1+1,j2+1), px(j1-1,j2+1))
            pxj3 = awestsouth(x_wght(l1,l2), px(j1,j4  ), px(j1+1,j4  ), px(j1-1,j4  ))
          ENDIF
          ! North-south interpolation
          IF (y_wght(l1,l2) > 0.5_ireals) THEN
            pxi(l1,l2) = aeastnorth(y_wght(l1,l2), pxj1, pxj2, pxj3)
          ELSE
            pxi(l1,l2) = awestsouth(y_wght(l1,l2), pxj1, pxj2, pxj3)
          ENDIF

          ! Enforce monotonicity, if required
          !IF (lmono) THEN ! Mathematically enforced
          !  pxi(l1,l2) = MIN (pxi(l1,l2), MAX (px(j1,j2),px(j1+1,j2), &
          !   px(j1,j2+1),px(j1+1,j2+1)))
          !  pxi(l1,l2) = MAX (pxi(l1,l2), MIN (px(j1,j2),px(j1+1,j2), &
          !   px(j1,j2+1),px(j1+1,j2+1)))
          !ENDIF
          ! Enforce positive definiteness, if required
          IF (lposdef) THEN
            pxi(l1,l2) = MAX (pxi(l1,l2), 0.0_ireals)
          ENDIF
        ENDDO
      ENDDO

      navg=COUNT(pxi /= undef .AND. lmask_lm)
      avg=SUM(pxi, MASK=(pxi /= undef .AND. lmask_lm))

      IF (nproc > 1) THEN
        CALL global_values (navg, 1, 'SUM', imp_integers, icomm_cart, -1, &
         yerrmsg, kierr)
        CALL global_values (avg,  1, 'SUM', imp_reals,    icomm_cart, -1, &
         yerrmsg, kierr)
      ENDIF
      IF (navg > 0) THEN
        avg = avg/navg
      ELSE
        avg = 0.0_ireals
      ENDIF
      WHERE (.NOT. lmask_lm)
        pxi = avg
      ENDWHERE
!     PRINT*,'AveragedQ ',navg, avg, COUNT(.NOT. lmask_lm)

    ! handle full area arrays
    ELSE

    ! check for at least 2 boundary lines and rows

      IF ( nzbounds_s > 1 .AND. nzbounds_n > 1 .AND. nzbounds_w > 1 .AND. nzbounds_s > 1) THEN

        CALL interp_q_lm(px, ie_in, je_in,                                    &
                    pxi, i_index, j_index,                                    &
                    x_wght, y_wght, kilone-kilons+1, kilate-kilats+1,         &
                    lmono, lposdef, yerrmsg, kierr)
        IF (kierr /= 0) RETURN

      ELSE

        ! southern border
        IF (nzbounds_s == 1) THEN
          ! count, for which j_s the minimum value of j_index has at least
          ! 2 boundary lines
          DO j_s=1,kilate-1
              nzbounds_s = MINVAL (j_index(kilons:kilone,kilats+j_s:kilate))
              IF (nzbounds_s > 1) EXIT
          ENDDO
          ! Allocate memory for boundary line array
          ALLOCATE (px_bl (ie_in, nzbounds_s), STAT = kierr)
          IF (kierr /= 0) THEN
            yerrmsg='allocation error px_bl Subroutine interp_q_bs'
            kierr = 1
          ENDIF
          px_bl(:,1:nzbounds_s) = px(:,1:nzbounds_s)

          CALL interp_l(px_bl, ie_in, nzbounds_s,                             &
                        i_index (kilons:kilone,kilats:kilats+j_s-1),          &
                        j_index (kilons:kilone,kilats:kilats+j_s-1),          &
                        lmono, lposdef, 'L', lframe, lolp_in(:,1:nzbounds_s), &
                        lolp_lm (kilons:kilone,kilats:kilats+j_s-1),          &
                        lmask_lm(kilons:kilone,kilats:kilats+j_s-1),          &
                        undef,                                                &
                        x_wght  (kilons:kilone,kilats:kilats+j_s-1),          &
                        y_wght  (kilons:kilone,kilats:kilats+j_s-1),          &
                        pxi     (kilons:kilone,kilats:kilats+j_s-1),          &
                                 kilons,kilone,kilats,kilats+j_s-1,           &
                        startlat_in, startlon_in, latitudes_in, longitudes_in,&
                        lat_coarse, lon_coarse,                               &
                        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,       &
                        yerrmsg, kierr)

          IF (kierr /= 0) RETURN
          DEALLOCATE (px_bl)
        ENDIF ! southern border

        ! northern border
        IF (nzbounds_n == 1) THEN
          ! count, for which j_n the maximum value of j_index has at least
          ! 2 boundary lines
          DO j_n=1,kilate-1
              nzbounds_n = je_in -MAXVAL (j_index(kilons:kilone,kilats:kilate-j_n)) + 1
              IF (nzbounds_n > 1) EXIT
          ENDDO
          ! Allocate memory for boundary line array
          ALLOCATE (px_bl (ie_in, nzbounds_n), STAT = kierr)
          IF (kierr /= 0) THEN
            yerrmsg='allocation error px_bl Subroutine interp_q_bs'
            kierr = 2
          ENDIF
          px_bl(:,1:nzbounds_n) = px(:,je_in-nzbounds_n+1:je_in)

          CALL interp_l(px_bl, ie_in, nzbounds_n,                             &
                        i_index (kilons:kilone,kilate-j_n+1:kilate),          &
                        j_index (kilons:kilone,kilate-j_n+1:kilate),          &
                        lmono, lposdef, 'L', lframe, lolp_in(:,1:nzbounds_s), &
                        lolp_lm (kilons:kilone,kilate-j_n+1:kilate),          &
                        lmask_lm(kilons:kilone,kilate-j_n+1:kilate),          &
                        undef,                                                &
                        x_wght  (kilons:kilone,kilate-j_n+1:kilate),          &
                        y_wght  (kilons:kilone,kilate-j_n+1:kilate),          &
                        pxi     (kilons:kilone,kilate-j_n+1:kilate),          &
                                kilons,kilone,kilate-j_n+1,kilate,            &
                        startlat_in, startlon_in, latitudes_in, longitudes_in,&
                        lat_coarse, lon_coarse,                               &
                        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,       &
                        yerrmsg, kierr)

          IF (kierr /= 0) RETURN
          DEALLOCATE (px_bl)
        ENDIF ! northern border

        ! western border
        IF (nzbounds_w == 1) THEN
          ! count, for which i_w the minimum value of i_index has at least
          ! 2 boundary rows
          DO i_w=1,kilone-1
              nzbounds_w = MINVAL (i_index(kilons+i_w:kilone,kilats:kilate))
              IF (nzbounds_w > 1) EXIT
          ENDDO

          ! Allocate memory for boundary line array
          ALLOCATE (px_bl (nzbounds_w, je_in), STAT = kierr)
          IF (kierr /= 0) THEN
            yerrmsg='allocation error px_bl Subroutine interp_q_bs'
            kierr = 3
          ENDIF
          px_bl(1:nzbounds_w,:) = px(1:nzbounds_w,:)
          CALL interp_l(px_bl, nzbounds_w, je_in,                             &
                        i_index (kilons:kilons+i_w-1,kilats:kilate),          &
                        j_index (kilons:kilons+i_w-1,kilats:kilate),          &
                        lmono, lposdef, 'L', lframe, lolp_in(:,1:nzbounds_w), &
                        lolp_lm (kilons:kilons+i_w-1,kilats:kilate),          &
                        lmask_lm(kilons:kilons+i_w-1,kilats:kilate),          &
                        undef,                                                &
                        x_wght  (kilons:kilons+i_w-1,kilats:kilate),          &
                        y_wght  (kilons:kilons+i_w-1,kilats:kilate),          &
                        pxi     (kilons:kilons+i_w-1,kilats:kilate),          &
                                 kilons,kilons+i_w-1,kilats,kilate,           &
                        startlat_in, startlon_in, latitudes_in, longitudes_in,&
                        lat_coarse, lon_coarse,                               &
                        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,       &
                        yerrmsg, kierr)

          IF (kierr /= 0) RETURN
          DEALLOCATE (px_bl)
        ENDIF ! western border

        ! eastern border
        IF (nzbounds_e == 1) THEN
          ! count, for which i_e the maximum value of i_index has at least
          ! 2 boundary rows
          DO i_e=1,kilone-1
              nzbounds_e = ie_in -MAXVAL (i_index(kilons:kilone-i_e,kilats:kilate)) + 1
              iF (nzbounds_e > 1) EXIT
          ENDDO
          ! Allocate memory for boundary line array
          ALLOCATE (px_bl (nzbounds_e, je_in), STAT = kierr)
          IF (kierr /= 0) THEN
            yerrmsg='allocation error px_bl Subroutine interp_q_bs'
            kierr = 4
          ENDIF

          px_bl(1:nzbounds_e,:) = px(ie_in-nzbounds_e+1:ie_in,:)

          CALL interp_l(px_bl, nzbounds_e, je_in,                             &
                        i_index (kilone-i_e+1:kilone,kilats:kilate),          &
                        j_index (kilone-i_e+1:kilone,kilats:kilate),          &
                        lmono, lposdef, 'L', lframe, lolp_in(:,1:nzbounds_e), &
                        lolp_lm (kilone-i_e+1:kilone,kilats:kilate),          &
                        lmask_lm(kilone-i_e+1:kilone,kilats:kilate),          &
                        undef,                                                &
                        x_wght  (kilone-i_e+1:kilone,kilats:kilate),          &
                        y_wght  (kilone-i_e+1:kilone,kilats:kilate),          &
                        pxi     (kilone-i_e+1:kilone,kilats:kilate),          &
                                 kilone-i_e+1,kilone,kilats,kilate,           &
                        startlat_in, startlon_in, latitudes_in, longitudes_in,&
                        lat_coarse, lon_coarse,                               &
                        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,       &
                        yerrmsg, kierr)

          IF (kierr /= 0) RETURN
          DEALLOCATE (px_bl)
        ENDIF ! eastern border

        CALL interp_q_lm(px, ie_in, je_in,                                            &
                        pxi    (kilons+i_w: kilone+i_e,    kilats+j_s: kilate-j_n),   &
                        i_index(kilons+i_w: kilone+i_e,    kilats+j_s: kilate-j_n),   &
                        j_index(kilons+i_w: kilone+i_e,    kilats+j_s: kilate-j_n),   &
                        x_wght (kilons+i_w: kilone+i_e,    kilats+j_s: kilate-j_n),   &
                        y_wght (kilons+i_w: kilone+i_e,    kilats+j_s: kilate-j_n),   &
                                kilone-i_e-(kilons+i_w)+1, kilate-j_n-(kilats+j_s)+1, &
                        lmono, lposdef, yerrmsg, kierr)
        IF (kierr /= 0) RETURN


      ENDIF

    ENDIF
! ELSEIF (lec2lm) THEN
!
  ENDIF
! UB

END SUBROUTINE interp_q_bs

!==============================================================================

END MODULE interp_utilities
