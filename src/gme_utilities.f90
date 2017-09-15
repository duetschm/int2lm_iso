!+ Source module for GME utility routines.
!==============================================================================

MODULE  gme_utilities

!==============================================================================
!
! Description:
!   This module provides service utilities for the GME. It deals with
!   parallel programming and other service routines.
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
!  Bug correction in SR distance
!  Editorial changes
! V1_7         2007/11/26 Ulrich Schaettler
!  Included SR get_ndvi to compute actual ndvi ratios from GME monthly means
! V1_8         2008/05/29 Ulrich Schaettler
!  Renamed ldebug to loc_debug
!  Adapted SR xd for working with GME bitmap data
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_20        2012/09/03 Ulrich Schaettler
!  Enlarged strings for date variables to 14 characters (SR get_ndvi)
! V1_21        2013/03/25 Ulrich Schaettler
!  Enlarged debug level value for certain outputs
!
! Code Description:
! Language:           Fortran 90.
! Software Standards: "European Standards for Writing and
!                      Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters, ONLY:   &
    ireals,   &! KIND-type parameters for real variables
    iintegers  ! kind-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local data for this module that are used in setup_xd and xd_p
! (has been taken from former common block communication.h)

! Maximum number of processors:
  INTEGER (KIND=iintegers), PARAMETER   :: maxprocs = 2048

  INTEGER (KIND=iintegers), PRIVATE, ALLOCATABLE ::        &
    idx_recv(:,:), & ! Index array where to put received points
    idx_send(:,:), & ! Index array from where to get points to send
    idx_help(:,:)


  INTEGER (KIND=iintegers), PRIVATE         ::      &
    np_recv_s(0:maxprocs-1), & ! start index of recieved points
                               ! in array idx_recv
    np_recv_1(0:maxprocs-1), & ! number of points to receive from a
                               ! processor - 1 boundary line
    np_recv_2(0:maxprocs-1), & ! dto - 2 boundary lines

    np_send_s(0:maxprocs-1), & ! start index of points to send
                               ! in array idx_send
    np_send_1(0:maxprocs-1), & ! number of points to send to a
                               ! processor - 1 boundary line
    np_send_2(0:maxprocs-1), & ! dto - 2 boundary lines

    np_recv_t(0:maxprocs-1), & ! start index of recieved points
                               ! in array idx_recv at the target PE
                               ! (provided for shmem-communication)
    np_recv_tot,             & ! tot. number of points to receive
    np_send_tot,             & ! tot. number of points to send
    np_recv_max,             & ! max. number of points to receive
    np_send_max                ! max. number of points to send

! Neighborhood relationships within the diamonds
! mpw, mpe, maw, mae give the numbers of the poleward or antipoleward
! west or east neighbors of a diamond

  INTEGER (KIND=iintegers), PRIVATE        ::      &
    mpw(10) = (/ 5, 1, 2, 3, 4,10, 6, 7, 8, 9/), & !
    mpe(10) = (/ 2, 3, 4, 5, 1, 7, 8, 9,10, 6/), & !
    maw(10) = (/10, 6, 7, 8, 9, 1, 2, 3, 4, 5/), & !
    mae(10) = (/ 6, 7, 8, 9,10, 2, 3, 4, 5, 1/)

!==============================================================================

! include statements
INCLUDE "mpif.h"

!=======================================================================

CONTAINS

!==============================================================================

SUBROUTINE check_bitmap (gme_field, index, ispoke, undef,                    &
                         j1start, j1end, j2start, j2end, ie2lm, je2lm, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine checks, whether all GME points that surround LM points are
!   properly initialized or whether there are undefined points, what could
!   happen, if a GME bitmap does not fit to the LM domain.
!
! Method:
!   Check for undefined values and set the error code.
!
!------------------------------------------------------------------------------

! Input to the routine
INTEGER  (KIND=iintegers), INTENT(IN)   ::                               &
  j1start, j1end, j2start, j2end, ie2lm, je2lm   ! Field dimensions

INTEGER  (KIND=iintegers), INTENT(IN)   ::                               &
  index (ie2lm,je2lm,4), & ! Information on GME field indices for LM gridpoints
  ispoke(12)

REAL     (KIND=ireals),    INTENT(IN)   ::                               &
  undef       ! value for undefined grid points

REAL     (KIND=ireals),    INTENT(IN)   ::                               &
  gme_field (j1start:j1end,j2start:j2end,1:10)  ! field that has to be checked

! Output: the error code
INTEGER  (KIND=iintegers), INTENT(OUT)  ::                               &
  ierror

! Local variables
INTEGER  (KIND=iintegers) :: i,j, j1, j2, jd, m1, m2
REAL     (KIND=ireals)    :: gme_value_1, gme_value_2, gme_value_3

LOGICAL                   :: lzdetect

!------------------------------------------------------------------------------

  ierror   = 0_iintegers
  lzdetect = .FALSE.

  ! Loops over the LM gridpoints
  DO j = 1, je2lm
    DO i = 1, ie2lm
      ! get information on GME grid points
      j1 = index(i,j,1)
      j2 = index(i,j,2)
      jd = index(i,j,3)
      m1 = index(i,j,4)
      m2 = MOD (m1,6) + 1

      ! check three GME points needed for normal interpolation
      gme_value_1 = gme_field(j1,           j2,             jd)
      gme_value_2 = gme_field(j1+ispoke(m1),j2+ispoke(m1+6),jd)
      gme_value_3 = gme_field(j1+ispoke(m2),j2+ispoke(m2+6),jd)


      IF (gme_value_1 == undef .OR. gme_value_2 == undef            &
                               .OR. gme_value_3 == undef) THEN
        lzdetect = .TRUE.
      ENDIF
    ENDDO
  ENDDO

  IF (lzdetect) ierror = 1_iintegers

END SUBROUTINE check_bitmap

!==============================================================================

SUBROUTINE factorni (kni, loc_debug, kni2, kni3, kierr)

!------------------------------------------------------------------------------
!
! Description:
!   *factorni* computes the factors of the integer input kni, assuming
!   that kni decomposes into kni3 factors (kni3 either 0 or 1) of "3"
!   and kni2 (kni2 > 0) factors of "2". The subroutine returns the
!   number of factors of "3", kni3, number of factors of "2", kni2,
!   and sets the error flag kierr=0 if kni can be expressed this way.
!   If kni cannot be expressed in this way, the error flag kierr is 
!   set to -1. 
!   If the debug flag (loc_debug) is set to .true., all input and output
!   variables are printed.
!
! Method:
!
!=======================================================================
!
! Input
! kni         INTEGER   number to be factorized
! loc_debug   LOGICAL   debug flag; if .true. print information
!
!=======================================================================
!     
! Output
! kni3     INTEGER   exponent of "3", either 0 or 1
! kni2     INTEGER   exponent of "2", kni2 > 0
! kierr    INTEGER   error flag (0: no error, -1: error)
!
!=======================================================================
!
! Dummy arguments
  LOGICAL  loc_debug
  INTEGER  (KIND=iintegers)   :: kni, kni2, kni3, kierr
!
!=======================================================================
!
! Local variables
  INTEGER  (KIND=iintegers)   ::     &
             mx   ! auxiliary variable, set to kni initially
!
!=======================================================================

  ! Start of the factorization of kni
  mx    = kni
  kierr = 0
  kni2  = 0
  kni3  = 0
  IF (loc_debug) THEN
    PRINT *,'  SUBROUTINE *factorni*, Input:  kni= ', kni
  ENDIF

  DO WHILE (mx.GT.1) 
    IF (MOD(mx,2).EQ.0) THEN
      kni2  = kni2 + 1
      mx    = mx/2
    ELSE IF (MOD(mx,3).EQ.0) THEN
      kni3  = kni3 + 1
      mx    = mx/3
    ELSE
      kierr = -1
      RETURN
    ENDIF
  ENDDO

  ! kni3 must not be greater than 1
  IF (kni3.GT.1) kierr = -1

  IF (loc_debug) THEN
    PRINT *,'  SUBROUTINE *factorni*, Output: kni3= ', kni3,   &
            '  kni2= ', kni2, '  kierr= ', kierr
  ENDIF

!=======================================================================

END SUBROUTINE factorni

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE gcpt (pxn, kig1s, kig1e, kig2s, kig2e, knd,            &
                 kjd, pgamma, ki1, kj1, ki2, kj2, ki , kj )

!------------------------------------------------------------------------------
!
! Description:
!   *gcpt* finds the point x along the shorter great circle
!   arc between points x1 and x2 on the unit sphere, such that
!   x and x1 subtend the angle pgamma*ztheta0, where ztheta0 is the
!   angle subtended by x1 and x2. Note gamma resides on the
!   interval [0-1]. The Cartesian coordinates of the point x are
!   stored in the array pxn at the indices ki, kj. 
!
! Method:
!
!=======================================================================
!
!     Input
!     pxn     REAL      pxn(kig1s:kig1e, kig2s:kig2e, 3, knd): Carte-
!                       sian (x,y,z) coordinates of the location vector
!                       of the gridpoints (nodes) on the unit-sphere;
!                       phys. dim. ( - )
!     kig1s   INTEGER   first  dimension of array, start index
!     kig1e   INTEGER   first  dimension of array, end   index
!     kig2s   INTEGER   second dimension of array, start index
!     kig2e   INTEGER   second dimension of array, end   index
!     knd     INTEGER   number of diamond (knd = 10)
!     kjd     INTEGER   index of actual diamond (kjd = 1 to 10)
!     pgamma  REAL      the normalized great circle angular parameter
!                       interpolating between x1 and x2.
!     ki1     INTEGER   index of point x1 in i1-direction
!     kj1     INTEGER   index of point x1 in i2-direction
!     ki2     INTEGER   index of point x2 in i1-direction
!     kj2     INTEGER   index of point x2 in i2-direction
!     ki      INTEGER   index of point x  in i1-direction
!     kj      INTEGER   index of point x  in i2-direction
!
!=======================================================================
!
!     Output 
!     pxn     REAL      see above
!
! Dummy arguments
  INTEGER  (KIND=iintegers)  ::    &
            kig1s, kig1e, kig2s, kig2e, knd ,              &
            kjd  , ki1  , ki2  , kj1  , kj2 , ki , kj

  REAL     (KIND=ireals)     ::    &
            pxn    (kig1s:kig1e, kig2s:kig2e, 3, knd)

  REAL     (KIND=ireals)     ::     pgamma

!=======================================================================

! Local variables
  REAL     (KIND=ireals)     ::    &
        ztheta,          & ! great circle angle
        zchord,          & ! length of the chord connecting x1,x2
        zalpha,          & ! weighting factor for vector x1
        zbeta              ! weighting factor for vector x2


!=======================================================================

  ! Calculate "zchord", the Cartesian distance between x1 and x2
  zchord = SQRT ( (pxn(ki2,kj2,1,kjd) - pxn(ki1,kj1,1,kjd))**2 +    &
                  (pxn(ki2,kj2,2,kjd) - pxn(ki1,kj1,2,kjd))**2 +    &
                  (pxn(ki2,kj2,3,kjd) - pxn(ki1,kj1,3,kjd))**2 )

  ! Calculate "ztheta", the great circle angle between x1 and x2 
  ztheta = 2.*ASIN (0.5_ireals*zchord)

  ! Calculate the weighting factors which follow from the condition
  ! that x is a point on the unit-sphere, too.
  zbeta  = SIN (pgamma*ztheta)/SIN (ztheta)
  zalpha = SIN ((1.0_ireals-pgamma)*ztheta)/SIN (ztheta)

  ! Store the (x,y,z) coordinates of the point x into the array pxn
  pxn(ki,kj,1,kjd) = zalpha*pxn(ki1,kj1,1,kjd) + zbeta *pxn(ki2,kj2,1,kjd)
  pxn(ki,kj,2,kjd) = zalpha*pxn(ki1,kj1,2,kjd) + zbeta *pxn(ki2,kj2,2,kjd)
  pxn(ki,kj,3,kjd) = zalpha*pxn(ki1,kj1,3,kjd) + zbeta *pxn(ki2,kj2,3,kjd)

!=======================================================================

END SUBROUTINE gcpt

!==============================================================================
!==============================================================================

SUBROUTINE get_ndvi (ndvi_mr, ndvi_ar,                                     &
                     kig1sm2, kig1ep2, kig2sm2, kig2ep2, kdmin, kdmax, km, &
                     ydate_act, iprdeb)

!=======================================================================
!
! Description:
!   The ratio of actual normalized differential vegetation index to its
!   annual mean NDVIRATIO is updated according ydate_act
!
! Method:
!   A linear interpolation is done for the actual forecast time using
!   the values (ndvi_mr) of the two nearest months.
!
!=======================================================================

! Parameter list
INTEGER(KIND=iintegers), INTENT(IN)      ::           &
   kig1sm2, kig1ep2, kig2sm2, kig2ep2, & ! spatial dimensions of output
   kdmin, kdmax,                       & ! min- and max triangle number
   km,                                 & ! 12 monthly mean values
   iprdeb                                ! for debug print outs

REAL(KIND=ireals),       INTENT(IN)      ::           &
   ndvi_mr(kig1sm2:kig1ep2, kig2sm2:kig2ep2, kdmin:kdmax, km)  ! Monthly NDVI

REAL(KIND=ireals),       INTENT(OUT)     ::           &
   ndvi_ar(kig1sm2:kig1ep2, kig2sm2:kig2ep2, kdmin:kdmax)      ! Actual  NDVI

CHARACTER (LEN=14),      INTENT(IN)      ::           &
   ydate_act                             ! actual date

! Local variables
INTEGER(KIND=iintegers)                  ::           &
  mmon, mday, mhour, myy, & ! month, day, hour and year of actual date
  mdayhour,               & ! actual date (in hours of month)
  mstart, mstartp1,       & ! index of start month, next index (for interp.)
  mmidthhours,            & ! midth of month in hours
  i, ip1,                 & ! month indices (ip1=i+1)
  j1, j2, jd,             & !
  mleapy                    ! leap year (1=yes, 0=no)

REAL(KIND=ireals)                        ::           &
  zdiff(12),     & ! difference between midth of following months in days
  zhalf(12),     & ! number of days for half month
  zact,          & ! actual time in hours of month
  zdt              ! dummy time step

! Number of days for each month
INTEGER(KIND=iintegers)                  ::           &
  month_days(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

! Statementfunction for leap year determination (1=yes, 0=no)
  mleapy(myy)    =  MAX(1-MODULO(myy,4)  ,0) - MAX(1-MODULO(myy,100),0) &
                                             + MAX(1-MODULO(myy,400),0)

!=======================================================================

  READ(ydate_act(1:10),'(I4,3I2)') myy, mmon, mday, mhour

! Compute half of each month (in days)
! Leap Year ??
  month_days (2) = 28 + mleapy(myy)
  zhalf(:) = 0.5_ireals * REAL (month_days(:), ireals)
  IF (iprdeb > 19) THEN
    PRINT *, ' Length of February:', month_days(2), mleapy(myy)
  ENDIF

! Compute difference between the midth of actual month and the
! following one (in days)
  DO i = 1,12
    ip1      = MOD(i,12)+1
    zdiff(i) = zhalf(i) + zhalf(ip1)
  ENDDO

! Compute actual date (day and hours) and midth of actual month in hours
  mdayhour    = (mday-1)*24+mhour
  mmidthhours = NINT( zhalf(mmon)*24.)
  IF (iprdeb > 19) THEN
    PRINT *, ' Update actual NDVI rates for date: ', ydate_act, '  ',   &
             mdayhour, '  ', mmidthhours
  ENDIF

! Determine the months needed for interpolation of ndvi-values
! Search for the position of date in relation to first of month
! The ndvi-values are valid for the midth of month
!
! EXAMPLE 1
!        March    !  April     !   May             X : ndvi_ratio
!       ----X-----!-----X----o-!-----X-----        ! : first of month
!                       !    ^       !             o : actual date
!                       !    ^ interpolation for that point
!                       !  zdiff(4)  !
!                       !zact!
!
! EXAMPLE 2
!        March    !  April     !   May             X : ndvi_ratio
!       ----X-----!-----X------!----oX-----        ! : first of month
!                       !           ^              o : actual date
!                       !      interpolation for that point
!                       !zhalf !
!                       !  zdiff(4)  !
!                       !   zact    !
!
!

  IF (mdayhour < mmidthhours) THEN
    ! point is in first half of month (EXAMPLE 2)
    mstart = mmon - 1
    IF(mmon == 1) mstart = 12
    zact   = zhalf(mstart) + REAL(mdayhour)/24.
  ELSE
    ! point is in second half of month (EXAMPLE 1)
    mstart = mmon
    zact   = REAL(mdayhour-mmidthhours)/24.
  ENDIF

  mstartp1 = mod(mstart,12) + 1

  IF (iprdeb > 19) THEN
    PRINT *, '  Compute actual NDVI rates between months: ', mstart, mstartp1
  ENDIF

  ! Temporal interpolation of ndvi (linear) from the monthly data
  DO jd = kdmin, kdmax
    DO j2 = kig2sm2, kig2ep2
      DO j1 = kig1sm2, kig1ep2
        ndvi_ar(j1,j2,jd) = ndvi_mr(j1,j2,jd,mstart) +                           &
                           (ndvi_mr(j1,j2,jd,mstartp1)-ndvi_mr(j1,j2,jd,mstart)) &
                                                     * zact/zdiff(mstart)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE get_ndvi

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE tricntr (pxn, kig1s, kig1e, kig2s, kig2e, knd, kjd, kni)

!------------------------------------------------------------------------------
!
! Description:
!   *tricntr* finds the coordinates of the node at the center
!   of the two icosahedral triangles comprising diamond kjd. These
!   nodes are present when ni contains the factor 3, i.e. if ni3 = 1.
!
! Method:
!
!=======================================================================
!
!     Input 
!     pxn     REAL      pxn(kig1s:kig1e, kig2s:kig2e, 3, knd): Carte-
!                       sian (x,y,z) coordinates of the location vector
!                       of the gridpoints (nodes) on the unit-sphere;
!                       phys. dim. ( - )
!     kig1s   INTEGER   first  dimension of array, start index
!     kig1e   INTEGER   first  dimension of array, end   index
!     kig2s   INTEGER   second dimension of array, start index
!     kig2e   INTEGER   second dimension of array, end   index
!     knd     INTEGER   number of diamonds (knd = 10)
!     kjd     INTEGER   index of actual diamond (kjd = 1 to 10)
!     kni     INTEGER   number of intervals on a main triangle side
!
!=======================================================================
!
!     Output 
!     pxn     REAL      see above
!
! Dummy arrays and variables
  INTEGER  (KIND=iintegers)   ::     &
            kig1s, kig1e, kig2s, kig2e, knd, kjd, kni

  REAL     (KIND=ireals)      ::     &
            pxn(kig1s:kig1e, kig2s:kig2e, 3, knd)

!=======================================================================

! Local variables
  REAL    (KIND=ireals)      ::     &
    zxnorm           ! norm temporary
 
  INTEGER (KIND=iintegers)   ::     &
    j,             & ! loop index
    mi1,           & ! index of center point
    mi2              ! index of top or bottom diamond corner

!=======================================================================

  DO j = 1,2   ! Loop over the two triangles
    mi1  = j*kni/3
    mi2  = 1 + (j - 1)*kni
 
    pxn(mi1,mi1+1,1,kjd) = pxn(mi2-1,mi2  ,1,kjd) +    &
                           pxn(kni  ,1    ,1,kjd) +    &
                           pxn(0    ,kni+1,1,kjd)
    pxn(mi1,mi1+1,2,kjd) = pxn(mi2-1,mi2  ,2,kjd) +    &
                           pxn(kni  ,1    ,2,kjd) +    &
                           pxn(0    ,kni+1,2,kjd)
    pxn(mi1,mi1+1,3,kjd) = pxn(mi2-1,mi2  ,3,kjd) +    &
                           pxn(kni  ,1    ,3,kjd) +    &
                           pxn(0    ,kni+1,3,kjd)

    ! Normalize to unit-sphere

    zxnorm = 1./SQRT (pxn(mi1,mi1+1,1,kjd)**2 +        &
                      pxn(mi1,mi1+1,2,kjd)**2 +        &
                      pxn(mi1,mi1+1,3,kjd)**2) 

    pxn(mi1,mi1+1,1,kjd) = zxnorm*pxn(mi1,mi1+1,1,kjd)
    pxn(mi1,mi1+1,2,kjd) = zxnorm*pxn(mi1,mi1+1,2,kjd)
    pxn(mi1,mi1+1,3,kjd) = zxnorm*pxn(mi1,mi1+1,3,kjd)
  ENDDO    ! End of loop over triangles

!=======================================================================

END SUBROUTINE tricntr

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE pp_ll2gp (xn,slon,slat,nip1,zx,zy,zz,spd_t,kd,kj1,kj2,sp, loc_debug)

!------------------------------------------------------------------------------
!
! Description:
!   *pp_ll2gp* finds the nearest grid point on the unit sphere of
!              the triangular grid for a point of specified
!              longitude (slon) and latitude (slat)
!
! Method:
!
!=======================================================================

  INTEGER (KIND=iintegers) :: &
    nip1                    ! grid mesh dimension

! Input
  REAL (KIND=ireals)       :: &
    xn(nip1,nip1,3,10),   & ! cartesian coordinates of all nodes
    slon,                 & ! longitude (degrees) of point     
    slat,                 & ! latitude  (degrees) of point     
    spd_t                   ! threshold value for scalar product 
                            ! for initiation of additional fine
                            ! scale search with routine ZU_FUSS

  LOGICAL                  :: &
    loc_debug               !

! Output:
  REAL (KIND=ireals)       :: &
    zx,zy,zz,             & ! cartesian coordinates of point
    sp                      ! scalar product between point and
                            ! nearest GME nodal point

  INTEGER (KIND=iintegers) :: &
    kd,                   & ! diamond containing point
    kj1,kj2                 ! nodal indices of nearest grid point

! Local variables
  INTEGER (KIND=iintegers) :: &
    itest, icalls, n, ni, n2, n3, ierr, n3p1, n2pn3, jd, kpa, j1top, j2top, &
    j1ll, j2ll, j1lr, j2lr, jj1top, jj2top, jj1ll, jj2ll, jj1lr, jj2lr,     &
    j1, j2, j1_ini, j2_ini

  REAL (KIND=ireals)       :: &
    sp_min, api, zlon, zlat, smax, zxpt, zypt, zzpt, zxat, zyat, zzat, spp, &
    spa, sptop, spll, splr, sp_ini

  DATA itest /0/
  DATA icalls/0/
  DATA sp_min/1.000_ireals/

! End of header
!=======================================================================

  api  = 2.*asin(1.)
  n    = nip1
  ni   = nip1-1
  zlon = slon*api/180.
  zlat = slat*api/180.
  zx   = cos(zlon)*cos(zlat)
  zy   = sin(zlon)*cos(zlat)
  zz   =           sin(zlat)

  IF (loc_debug) WRITE (81,*) 'Suche des Punktes  ',slon,slat

  CALL factor3b(ni,n2,n3,ierr)
  IF (loc_debug) print *,' Ergebnis von FACTOR3   n2=',n2
  IF (loc_debug) print *,'                        n3=',n3
  n3p1 = n3+1
  n2pn3= n2+n3

  smax = -999.
  kd   = 0
  DO jd=1,10
    zxpt  = xn(1,1,1,jd)+xn(1,n,1,jd)+xn(n,1,1,jd)    ! Zentrum des pol-
    zypt  = xn(1,1,2,jd)+xn(1,n,2,jd)+xn(n,1,2,jd)    ! waertigen Icosa-
    zzpt  = xn(1,1,3,jd)+xn(1,n,3,jd)+xn(n,1,3,jd)    ! dreiecks
    zxat  = xn(n,n,1,jd)+xn(1,n,1,jd)+xn(n,1,1,jd)    ! Zentrum des anti-
    zyat  = xn(n,n,2,jd)+xn(1,n,2,jd)+xn(n,1,2,jd)    ! polwaertigen
    zzat  = xn(n,n,3,jd)+xn(1,n,3,jd)+xn(n,1,3,jd)    ! Icosaederdreiecks
    spp   = zx*zxpt + zy*zypt + zz*zzpt   
    spa   = zx*zxat + zy*zyat + zz*zzat   
    IF (spp.gt.smax ) THEN
      kd    = jd
      kpa   =  1
      smax  = spp
    ENDIF
    IF (spa.gt.smax ) THEN
      kd    = jd
      kpa   = -1
      smax  = spa
    ENDIF
  ENDDO

  IF (loc_debug) THEN
    print *,' grid point of longitude:',slon,'  latitude: ',slat
    print *,' belongs to diamaond       : ',kd,' !!!!'
    IF (kpa.eq.1) THEN
      print *,' in the poleward main triangle    '
    ELSE
      print *,' in anti-poleward main triangle  '
    ENDIF
    print *,'****************************************'
  ENDIF
  IF (kd.eq.0) stop 'fehler'

  IF (kpa.eq.1) THEN     ! Point im polward triangle   
    ! Binaersuche nach Unterdreieck der naechsten Stufe
    j1top = 1     ! Zeilenindex der polwaertigen Dreiecksspitze
    j2top = 1     ! Spaltenindex der polwaertigen Dreiecksspitze
    j1ll  = n     ! Zeilenindex der linken unteren Dreiecksspitze
    j2ll  = 1     ! Spaltenindex der linken unteren Dreiecksspitze
    j1lr  = 1     ! Zeilenindex der rechten unteren Dreiecksspitze
    j2lr  = n     ! Spaltenindex der rechten unteren Dreiecksspitze
  ELSE
    j1top = n     ! Zeilenindex der antipolwaertigen Dreiecksspitze
    j2top = n     ! Spaltenindex der antipolwaertigen Dreiecksspitze
    j1ll  = 1     ! Zeilenindex der linken unteren Dreiecksspitze
    j2ll  = n     ! Spaltenindex der linken unteren Dreiecksspitze
    j1lr  = n     ! Zeilenindex der rechten unteren Dreiecksspitze
    j2lr  = 1     ! Spaltenindex der rechten unteren Dreiecksspitze
  ENDIF

  IF (n3.eq.1) THEN
    call sub_t9(xn, n, kd, zx, zy, zz, j1top, j2top, j1ll, j2ll,      &
                j1lr, j2lr, loc_debug, jj1top, jj2top, jj1ll ,jj2ll, jj1lr, jj2lr)
     
    ! Swap: sub triangle corners ---> new main triangle corners
    j1top = jj1top
    j1ll  = jj1ll
    j1lr  = jj1lr
    j2top = jj2top
    j2ll  = jj2ll
    j2lr  = jj2lr
  ENDIF

  CALL sub_t4(n3p1,n2pn3,                                               &
              xn,n,kd,zx,zy,zz,                                         &
              j1top ,j2top  , j1ll  ,j2ll   , j1lr  ,j2lr   ,loc_debug, &
              jj1top,jj2top , jj1ll ,jj2ll  , jj1lr ,jj2lr )
     
  ! At end of subdivision sequence find nearest grid point:

  sptop = zx*xn(j1top,j2top,1,kd) +zy*xn(j1top,j2top,2,kd)    &
         +zz*xn(j1top,j2top,3,kd)
  spll  = zx*xn(j1ll ,j2ll ,1,kd) +zy*xn(j1ll ,j2ll ,2,kd)    &
         +zz*xn(j1ll ,j2ll ,3,kd)
  splr  = zx*xn(j1lr ,j2lr ,1,kd) +zy*xn(j1lr ,j2lr ,2,kd)    &
         +zz*xn(j1lr ,j2lr ,3,kd)

  IF (sptop.ge.spll.and.sptop.ge.splr) THEN
    j1 = j1top
    j2 = j2top
    sp = sptop
  ELSEIF (spll.ge.splr) THEN
    j1 = j1ll
    j2 = j2ll
    sp = spll 
  ELSE
    j1 = j1lr
    j2 = j2lr
    sp = splr 
  ENDIF

  IF (sp.lt.spd_t) THEN
    j1_ini = j1
    j2_ini = j2
    sp_ini = sp
    call zu_fuss   (xn,n,kd,zx,zy,zz,                          &
                    j1_ini,j2_ini,sp_ini, loc_debug, j1 ,j2  ,sp)

  ENDIF

  kj1 = j1-1  ! indices in calling routine     
  kj2 = j2    ! indices in calling routine

END SUBROUTINE pp_ll2gp

!===================================================================

SUBROUTINE factor3b (ni,n2,n3,ierr)

!------------------------------------------------------------------------------
!
! Description:
!   Subroutine factor3 computes factors of the integer input ni, 
!   assuming ni decomposes into one or zero factors of three and 
!   n2 factors of 2. The routine returns the number of factors of 
!   3, n3, the number of factors of 2, n2, and ierr=0 if ni can 
!   be expressed in this way. 
!
!   If ni cannot be expressed this way, ierr is set to -1.
!
!===================================================================
    
  INTEGER (KIND=iintegers) :: ni   ! number to be factored
      
  INTEGER (KIND=iintegers) :: &
           n2,         & !  number of factors of 2
           n3,         & !  number of factors of 3
           ierr          !  error value

  INTEGER (KIND=iintegers) :: ix
      
!===================================================================
!
!     Loop until ni has been decomposed into factors of 2 and 3.
!
!===================================================================

      ix=ni
      ierr=0
      n2=0
      n3=0
      do while( ix.gt.1)
         if ( MOD(ix,2) .eq. 0) then
             n2=n2+1
             ix=ix/2
         else if ( MOD(ix,3) .eq. 0) then
             n3=n3+1
             ix=ix/3
         else
             ierr= -1
         end if
      end do
      if (n3.gt.1) ierr=-1
      
END SUBROUTINE factor3b

!===================================================================
!
SUBROUTINE sub_t4(n3p1,n2pn3,                                   &
                  xn,n,kd,zx,zy,zz, j1top,j2top,                &
                  j1ll ,j2ll , j1lr ,j2lr ,loc_debug,           &
                  jj1top,jj2top, jj1ll ,jj2ll , jj1lr ,jj2lr )

!------------------------------------------------------------------------------
!
! Description:
!   *sub_t4* determines the grid point indices of a subtriangle
!            which contains a specified point on the unit sphere
!            The sub triangle is one of four which partition a
!            known main triangle, defined by the three pairs of
!            grid point indices of it's corners (top, lower left,
!            lower right). The main triangle must contain the
!            the point which is searched for!
!
!===================================================================
    
  INTEGER (KIND=iintegers)  ::    &
    n, kd, j1top,  j2top,  j1ll , j2ll,  j1lr,  j2lr,          &
           jj1top, jj2top, jj1ll, jj2ll, jj1lr, jj2lr, n3p1,n2pn3,ji

  real    (KIND=ireals)     ::    &
    xn(n,n,3,10)    ! cartesian coordinates of all grid points

  integer (KIND=iintegers)  ::    &
    j1(3),j2(3)     ! indices of new subtriangle corner points

  real    (KIND=ireals)     ::    &
    xtc(4),       & ! x-coordinate of subtriangle centers 
    ytc(4),       & ! y-coordinate of subtriangle centers  
    ztc(4)          ! z-coordinate of subtriangle centers 

  logical loc_debug

! Local variables
  INTEGER (KIND=iintegers) :: jd, ide

  REAL (KIND=ireals)       :: zx, zy, zz, spmin, rnorm, sp

!===================================================================
!
!     Structure of main triangle and the points defining the subtriangles
!
!
!     a) polward triangle             j1top,j2top
!                                         / \
!                                        /   \
!                                       /     \
!                                      /       \
!                                     /         \
!                                    /    D1     \
!                                   /             \
!                                  /               \
!                                 /                 \
!                       j1(1),j2(1)-----------------j1(3),j2(3)
!                               / \                 / \
!                              /   \               /   \
!                             /     \    D4       /     \
!                            /       \           /       \
!                           /         \         /         \
!                          /    D2     \       /     D3    \
!                         /             \     /             \
!                        /               \   /               \
!                       /                 \ /                 \
!                  j1ll,j2ll----------j1(2),j2(2)----------j1lr,j2lr
!
!
!     b) antipolward triangle: analogous to case a), but top is oriented
!                              away from pole and lower left refers to
!                              'eastern' corner, lower right refers to
!                              'western' corner
!
!===================================================================

 DO ji=n3p1,n2pn3

!     Subindices
      j1(1) = j1top + 0.5*(j1ll-j1top)
      j2(1) = j2top + 0.5*(j2ll-j2top) 
      j1(2) = j1top + 0.5*(j1ll-j1top)
      j2(2) = j2top + 0.5*(j2lr-j2top) 
      j1(3) = j1top + 0.5*(j1lr-j1top)
      j2(3) = j2top + 0.5*(j2lr-j2top) 

  ! Zentren der vier Subdreiecke:

  xtc(1)  = xn(j1top,j2top,1,kd) + xn(j1(1),j2(1),1,kd) + xn(j1(3),j2(3),1,kd)
  ytc(1)  = xn(j1top,j2top,2,kd) + xn(j1(1),j2(1),2,kd) + xn(j1(3),j2(3),2,kd)
  ztc(1)  = xn(j1top,j2top,3,kd) + xn(j1(1),j2(1),3,kd) + xn(j1(3),j2(3),3,kd)

  xtc(2)  = xn(j1(1),j2(1),1,kd) + xn(j1ll ,j2ll ,1,kd) + xn(j1(2),j2(2),1,kd)
  ytc(2)  = xn(j1(1),j2(1),2,kd) + xn(j1ll ,j2ll ,2,kd) + xn(j1(2),j2(2),2,kd)
  ztc(2)  = xn(j1(1),j2(1),3,kd) + xn(j1ll ,j2ll ,3,kd) + xn(j1(2),j2(2),3,kd)

  xtc(3)  = xn(j1(3),j2(3),1,kd) + xn(j1(2),j2(2),1,kd) + xn(j1lr ,j2lr ,1,kd)
  ytc(3)  = xn(j1(3),j2(3),2,kd) + xn(j1(2),j2(2),2,kd) + xn(j1lr ,j2lr ,2,kd) 
  ztc(3)  = xn(j1(3),j2(3),3,kd) + xn(j1(2),j2(2),3,kd) + xn(j1lr ,j2lr ,3,kd)

  xtc(4)  = xn(j1(2),j2(2),1,kd) + xn(j1(1),j2(1),1,kd) + xn(j1(3),j2(3),1,kd)
  ytc(4)  = xn(j1(2),j2(2),2,kd) + xn(j1(1),j2(1),2,kd) + xn(j1(3),j2(3),2,kd)
  ztc(4)  = xn(j1(2),j2(2),3,kd) + xn(j1(1),j2(1),3,kd) + xn(j1(3),j2(3),3,kd)

  spmin = 0.0
  do jd=1,4
     rnorm = 1./sqrt(xtc(jd)**2+ytc(jd)**2+ztc(jd)**2)
     xtc(jd) = rnorm*xtc(jd)     
     ytc(jd) = rnorm*ytc(jd)     
     ztc(jd) = rnorm*ztc(jd)     
     sp      = xtc(jd)*zx+ytc(jd)*zy+ztc(jd)*zz
     if (sp.gt.spmin) then
         ide     = jd
         spmin   = sp
     endif
  enddo
      
  if (ide.eq.1) then                       ! Point in D1
     jj1top = j1top
     jj2top = j2top
     jj1ll  = j1(1)
     jj2ll  = j2(1)
     jj1lr  = j1(3)
     jj2lr  = j2(3)
  else if (ide.eq.2) then  ! Point in D2
     jj1top = j1(1)
     jj2top = j2(1)
     jj1ll  = j1ll 
     jj2ll  = j2ll 
     jj1lr  = j1(2)
     jj2lr  = j2(2)
  else if (ide.eq.3) then ! Point in D3
     jj1top = j1(3)
     jj2top = j2(3)
     jj1ll  = j1(2)
     jj2ll  = j2(2)
     jj1lr  = j1lr
     jj2lr  = j2lr
  else                                     ! Point in D4
     jj1top = j1(2)
     jj2top = j2(2)
     jj1ll  = j1(3)
     jj2ll  = j2(3)
     jj1lr  = j1(1)
     jj2lr  = j2(1)
  end if

  j1top = jj1top
  j1ll  = jj1ll
  j1lr  = jj1lr
  j2top = jj2top
  j2ll  = jj2ll
  j2lr  = jj2lr
 
 ENDDO

END SUBROUTINE sub_t4

!===================================================================

SUBROUTINE sub_t9(xn,n,kd,zx,zy,zz,                                     &
                       j1top,j2top, j1ll ,j2ll , j1lr ,j2lr ,loc_debug, &
                       jj1top,jj2top, jj1ll ,jj2ll , jj1lr ,jj2lr )

!------------------------------------------------------------------------------
!
! Description:
!   *sub_t9* determines the grid point indices of a subtriangle
!            which contains a specified point on the unit sphere
!            The sub triangle is one of nine which partition a
!            known main triangle, defined by the three pairs of
!            grid point indices of it's corners (top, lower left,
!            lower right). The main triangle must contain the
!            the point which is searched for!
!
!===================================================================
    
  INTEGER (KIND=iintegers) :: n, kd, j1top, j2top, j1ll, j2ll,    &
                              j1lr, j2lr, jj1top, jj2top, jj1ll,  &
                              jj2ll, jj1lr, jj2lr

  real    (KIND=ireals)      ::    &
        xn(n,n,3,10), & ! cartesian coordinates of all grid points
        xtc(9),       & ! x-coordinate of subtriangle centers 
        ytc(9),       & ! y-coordinate of subtriangle centers  
        ztc(9)          ! z-coordinate of subtriangle centers 

  logical loc_debug

! Local variables
  INTEGER (KIND=iintegers)   ::    &
           j1(7),j2(7)  ! indices of new subtriangle corner points

  INTEGER (KIND=iintegers) :: jd, ide

  REAL (KIND=ireals)       :: zx, zy, zz, spmin, rnorm, sp

!===================================================================
!
!     Structure of main triangle and the points defining the subtriangles
!
!
!     a) polward triangle             j1top,j2top
!                                         / \
!                                        /   \
!                                       /     \
!                                      /       \
!                                     /         \
!                                    /    D1     \
!                                   /             \
!                                  /               \
!                                 /                 \
!                           j1(1),j2(1)---------j1(6),j2(6)
!                               / \                 / \
!                              /   \               /   \
!                             /     \     D7      /     \
!                            /       \           /       \
!                           /         \         /         \
!                          /    D2     \       /     D6    \
!                         /             \     /             \
!                        /               \   /               \
!                       /                 \ /                 \
!                 j1(2),j2(2)---------j1(7),j2(7)---------j1(5),j2(5)
!                     / \                 / \                 / \
!                    /   \               /   \               /   \
!                   /     \     D8      /     \      D9     /     \
!                  /       \           /       \           /       \
!                 /         \         /         \         /         \
!                /    D3     \       /    D4     \       /    D5     \
!               /             \     /             \     /             \
!              /               \   /               \   /               \
!             /                 \ /                 \ /                 \
!        j1ll,j2ll----------j1(3),j2(3)---------j1(4),j2(4)----------j1lr,j2lr
!
!
!     b) antipolward triangle: analogous to case a), but top is oriented
!                              away from pole and lower left refers to
!                              'eastern' corner, lower right refers to
!                              'western' corner
!
!===================================================================

     if (loc_debug) then
       write(*,*) ' Indizes der HauptPointe   : '
       write(*,*) ' Top    .: ',j1top ,j2top
       write(*,*) ' Links u : ',j1ll  ,j2ll  
       write(*,*) ' Rechts u: ',j1lr  ,j2lr
     end if

!     Subindices
      j1(1) = j1top +   (j1ll-j1top)/3
      j2(1) = j2top +   (j2ll-j2top)/3
      j1(2) = j1top + 2*(j1ll-j1top)/3
      j2(2) = j2top + 2*(j2ll-j2top)/3
      j1(3) = j1top + 2*(j1ll-j1top)/3
      j2(3) = j2top +   (j2lr-j2top)/3
      j1(4) = j1top +   (j1ll-j1top)/3
      j2(4) = j2top + 2*(j2lr-j2top)/3
      j1(5) = j1top + 2*(j1lr-j1top)/3
      j2(5) = j2top + 2*(j2lr-j2top)/3
      j1(6) = j1top +   (j1lr-j1top)/3
      j2(6) = j2top +   (j2lr-j2top)/3
      j1(7) = (j1(3)+j1(6))/2          
      j2(7) = (j2(1)+j2(4))/2          

       if (loc_debug) then
       write(*,*) '        Indices of subpoints:       '
       write(*,*) '        Point 1: ',j1(1),j2(1)
       write(*,*) '        Point 2: ',j1(2),j2(2)
       write(*,*) '        Point 3: ',j1(3),j2(3)
       write(*,*) '        Point 4: ',j1(4),j2(4)
       write(*,*) '        Point 5: ',j1(5),j2(5)
       write(*,*) '        Point 6: ',j1(6),j2(6)
       write(*,*) '        Point 7: ',j1(7),j2(7)
       end if

  ! Zentren der neun Subdreiecke:

  xtc(1)  = xn(j1top,j2top,1,kd) + xn(j1(1),j2(1),1,kd) + xn(j1(6),j2(6),1,kd)
  ytc(1)  = xn(j1top,j2top,2,kd) + xn(j1(1),j2(1),2,kd) + xn(j1(6),j2(6),2,kd)
  ztc(1)  = xn(j1top,j2top,3,kd) + xn(j1(1),j2(1),3,kd) + xn(j1(6),j2(6),3,kd)

  xtc(2)  = xn(j1(1),j2(1),1,kd) + xn(j1(2),j2(2),1,kd) + xn(j1(7),j2(7),1,kd)
  ytc(2)  = xn(j1(1),j2(1),2,kd) + xn(j1(2),j2(2),2,kd) + xn(j1(7),j2(7),2,kd)
  ztc(2)  = xn(j1(1),j2(1),3,kd) + xn(j1(2),j2(2),3,kd) + xn(j1(7),j2(7),3,kd)

  xtc(3)  = xn(j1(3),j2(3),1,kd) + xn(j1(2),j2(2),1,kd) + xn(j1ll ,j2ll ,1,kd)
  ytc(3)  = xn(j1(3),j2(3),2,kd) + xn(j1(2),j2(2),2,kd) + xn(j1ll ,j2ll ,2,kd)
  ztc(3)  = xn(j1(3),j2(3),3,kd) + xn(j1(2),j2(2),3,kd) + xn(j1ll ,j2ll ,3,kd)

  xtc(4)  = xn(j1(7),j2(7),1,kd) + xn(j1(3),j2(3),1,kd) + xn(j1(4),j2(4),1,kd)
  ytc(4)  = xn(j1(7),j2(7),2,kd) + xn(j1(3),j2(3),2,kd) + xn(j1(4),j2(4),2,kd)
  ztc(4)  = xn(j1(7),j2(7),3,kd) + xn(j1(3),j2(3),3,kd) + xn(j1(4),j2(4),3,kd)

  xtc(5)  = xn(j1(5),j2(5),1,kd) + xn(j1(4),j2(4),1,kd) + xn(j1lr ,j2lr ,1,kd)
  ytc(5)  = xn(j1(5),j2(5),2,kd) + xn(j1(4),j2(4),2,kd) + xn(j1lr ,j2lr ,2,kd)
  ztc(5)  = xn(j1(5),j2(5),3,kd) + xn(j1(4),j2(4),3,kd) + xn(j1lr ,j2lr ,3,kd)

  xtc(6)  = xn(j1(6),j2(6),1,kd) + xn(j1(7),j2(7),1,kd) + xn(j1(5),j2(5),1,kd)
  ytc(6)  = xn(j1(6),j2(6),2,kd) + xn(j1(7),j2(7),2,kd) + xn(j1(5),j2(5),2,kd)
  ztc(6)  = xn(j1(6),j2(6),3,kd) + xn(j1(7),j2(7),3,kd) + xn(j1(5),j2(5),3,kd)

  xtc(7)  = xn(j1(1),j2(1),1,kd) + xn(j1(7),j2(7),1,kd) + xn(j1(6),j2(6),1,kd)
  ytc(7)  = xn(j1(1),j2(1),2,kd) + xn(j1(7),j2(7),2,kd) + xn(j1(6),j2(6),2,kd)
  ztc(7)  = xn(j1(1),j2(1),3,kd) + xn(j1(7),j2(7),3,kd) + xn(j1(6),j2(6),3,kd)

  xtc(8)  = xn(j1(2),j2(2),1,kd) + xn(j1(3),j2(3),1,kd) + xn(j1(7),j2(7),1,kd)
  ytc(8)  = xn(j1(2),j2(2),2,kd) + xn(j1(3),j2(3),2,kd) + xn(j1(7),j2(7),2,kd)
  ztc(8)  = xn(j1(2),j2(2),3,kd) + xn(j1(3),j2(3),3,kd) + xn(j1(7),j2(7),3,kd)

  xtc(9)  = xn(j1(7),j2(7),1,kd) + xn(j1(4),j2(4),1,kd) + xn(j1(5),j2(5),1,kd)
  ytc(9)  = xn(j1(7),j2(7),2,kd) + xn(j1(4),j2(4),2,kd) + xn(j1(5),j2(5),2,kd)
  ztc(9)  = xn(j1(7),j2(7),3,kd) + xn(j1(4),j2(4),3,kd) + xn(j1(5),j2(5),3,kd)

         spmin = 0.0
         do jd=1,9
         rnorm = 1./sqrt(xtc(jd)**2+ytc(jd)**2+ztc(jd)**2)
         xtc(jd) = rnorm*xtc(jd)     
         ytc(jd) = rnorm*ytc(jd)     
         ztc(jd) = rnorm*ztc(jd)     
         sp      = xtc(jd)*zx+ytc(jd)*zy+ztc(jd)*zz
!       if (loc_debug) then
!       write (12,*) 'Skalarprodukt fuer Dreieck : ',jd,' = ',sp
!       end if
         if (sp.gt.spmin) then
         ide     = jd
         spmin   = sp
         end if
         end do
!       if (loc_debug) then
!        write(12,*) 'selected subtriangle    : ',ide
!       end if
      
        if (ide.eq.1) then                       ! Point in D1
!       if (loc_debug) then
!       write (12,*) 'Point in subtriangle no.1 *****'  
!       end if
        jj1top = j1top
        jj2top = j2top
        jj1ll  = j1(1)
        jj2ll  = j2(1)
        jj1lr  = j1(6)
        jj2lr  = j2(6)
        else if (ide.eq.2) then  ! Point in D2
!       if (loc_debug) then
!       write (12,*) 'Point in subtriangle no.2 *****'  
!       end if
        jj1top = j1(1)
        jj2top = j2(1)
        jj1ll  = j1(2)
        jj2ll  = j2(2)
        jj1lr  = j1(7)
        jj2lr  = j2(7)
        else if (ide.eq.3) then ! Point in D3
!       if (loc_debug) then
!       write (12,*) 'Point in subtriangle no.3 *****'  
!       end if
        jj1top = j1(2)
        jj2top = j2(2)
        jj1ll  = j1ll
        jj2ll  = j2ll
        jj1lr  = j1(3)
        jj2lr  = j2(3)
        else if(ide.eq.4) then                   ! Point in D4
!       if (loc_debug) then
!       write (12,*) 'Point in subtriangle no.4 *****'  
!       end if
        jj1top = j1(7)
        jj2top = j2(7)
        jj1ll  = j1(3)
        jj2ll  = j2(3)
        jj1lr  = j1(4)
        jj2lr  = j2(4)
        else if(ide.eq.5) then                   ! Point in D5
!       if (loc_debug) then
!       write (12,*) 'Point in subtriangle no.5 *****'  
!       end if
        jj1top = j1(5)
        jj2top = j2(5)
        jj1ll  = j1(4)
        jj2ll  = j2(4)
        jj1lr  = j1lr  
        jj2lr  = j2lr  
        else if(ide.eq.6) then                   ! Point in D6
!       if (loc_debug) then
!       write (12,*) 'Point in subtriangle no.6 *****'  
!       end if
        jj1top = j1(6)
        jj2top = j2(6)
        jj1ll  = j1(7)
        jj2ll  = j2(7)
        jj1lr  = j1(5)
        jj2lr  = j2(5)
        else if(ide.eq.7) then                   ! Point in D7
!       if (loc_debug) then
!       write (12,*) 'Point in subtriangle no.7 *****'  
!       end if
        jj1top = j1(7)
        jj2top = j2(7)
        jj1ll  = j1(6)
        jj2ll  = j2(6)
        jj1lr  = j1(1)
        jj2lr  = j2(1)
        else if(ide.eq.8) then                   ! Point in D8
!       if (loc_debug) then
!       write (12,*) 'Point in subtriangle no.8 *****'  
!       end if
        jj1top = j1(3)
        jj2top = j2(3)
        jj1ll  = j1(7)
        jj2ll  = j2(7)
        jj1lr  = j1(2)
        jj2lr  = j2(2)
        else if(ide.eq.9) then                   ! Point in D9
!       if (loc_debug) then
!       write (12,*) 'Point in subtriangle no.9 *****'  
!       end if
        jj1top = j1(4)
        jj2top = j2(4)
        jj1ll  = j1(5)
        jj2ll  = j2(5)
        jj1lr  = j1(7)
        jj2lr  = j2(7)

        else                      ! Fehler
        print *,' Fehler in SUB_t9 ide=',ide
        stop 'fehler'
        end if

END SUBROUTINE sub_t9

!===================================================================

SUBROUTINE distance (xn,ni,spd,spd_t)                     

!------------------------------------------------------------------------------
!
! Description:
!
!===================================================================

! Global and local variables
  INTEGER (KIND=iintegers)  :: ni
  REAL (KIND=ireals)        :: xh, yh, zh, spd, spd_t
  REAL (KIND=ireals)        :: xn(0:ni,1:ni+1,3,10)

  spd = xn(0,1,1,1) * xn (1,1,1,1) + xn(0,1,2,1) * xn (1,1,2,1)   &
      + xn(0,1,3,1) * xn (1,1,3,1)
      
  xh  = 0.5_ireals*(xn(0,1,1,1) + xn (1,1,1,1))
  yh  = 0.5_ireals*(xn(0,1,2,1) + xn (1,1,2,1))
  zh  = 0.5_ireals*(xn(0,1,3,1) + xn (1,1,3,1))
  spd_t = xn(0,1,1,1) * xh + xn(0,1,2,1) * yh + xn(0,1,3,1) * zh 

  spd_t = spd_t / SQRT(xh**2 + yh**2 + zh**2)

END SUBROUTINE distance

!===================================================================

SUBROUTINE xyzll(xn,ni,jd,j1,j2,x,y,z,plon,plat)

!------------------------------------------------------------------------------
!
! Description:
!
!===================================================================

! Global and local variables
  INTEGER (KIND=iintegers) :: ni, jd, j1, j2
  REAL (KIND=ireals)       :: x, y, z, plon, plat, api
  REAL (KIND=ireals)       :: xn(0:ni,1:ni+1,3,10)

      api = 2.0_ireals * ASIN(1.0_ireals)

!     Cartesian coordinates of selcted grid point
      x  = xn(j1,j2,1,jd)
      y  = xn(j1,j2,2,jd)
      z  = xn(j1,j2,3,jd)

!     Calculate the longitude "plon" and the latitude "plat";
      plon  = ATAN2 (xn(j1,j2,2,jd),xn(j1,j2,1,jd) + 1.e-20_ireals)
      plat  = ASIN  (xn(j1,j2,3,jd))

      plon = plon*180./api
      plat = plat*180./api

END SUBROUTINE xyzll

!===================================================================

SUBROUTINE zu_fuss   (xn,n,kd,zx,zy,zz, j1_ini,j2_ini,sp_ini,    &
                      loc_debug, jj1_1,jj2_1 ,sp_1)

!------------------------------------------------------------------------------
!
! Description:
!   *zu_fuss*  determines the grid point indices of a GME nodal point
!              which is closest to a specified point on the unit sphere
!===================================================================

! Global and local variables
  INTEGER (KIND=iintegers) :: n, kd, j1_ini, j2_ini, jj1_1, jj2_1,  &
                              j1_old, j2_old, jstart, jloop, jdir,  &
                              j1_new, j2_new

  REAL (KIND=ireals)       :: sp_ini, sp_1, zx, zy, zz, sp_old, sp_new

  REAL (KIND=ireals)       ::   &
           xn(n,n,3,10)   ! cartesian coordinates of all grid points

  INTEGER (KIND=iintegers) ::   &
           j1n(11),j2n(11) ! index increments for relevant neigbours
                           ! of any GME nodal point
                           ! as each node has only six neighbours,
                           ! elements 7 to 11 are just repetitions
                           ! of elements 1 to 5 to allow each neigbour
                           ! (from 1 to 6) as starting point
                           !  in a search loop, cf. Fig.1:
                           !
                           !      5---------4
                           !     / \       / \            \       /
                           !    /   \     /   \            \     /
                           !   /     \   /     \            \   /
                           !  /       \ /       \            \ /
                           ! 6---------*---------3          j1,j2
                           !  \       / \       /            / \
                           !   \     /   \     /            /   \
                           !    \   /     \   /            /     \
                           !     \ /       \ /            /       \
                           !      1---------2           j1+1     j2+1


  logical loc_debug

!===================================================================

!     Spoke No.  1  2  3  4  5  6    1  2  3  4  5
!     ==============================================
      data j1n / 1, 0,-1,-1, 0, 1,   1, 0,-1,-1, 0 / 
      data j2n / 0, 1, 1, 0,-1,-1,   0, 1, 1, 0,-1 /

      j1_old  = j1_ini                   ! 1.index of first guess
      j2_old  = j2_ini                   ! 2.index of first guess
      sp_old  = sp_ini                   ! first guess scalar product
      jstart  = 1                        ! initial search direction

!      Walking loop                                      
!      ============

       DO 2000 jloop=1,100000

        if (loc_debug) then
        write(12,*) jloop,'.Durchlauf der Suchschleife in RWALK '
        write(12,*) '-------------------------------------------'
        write(12,*) ' aktueller Gitterpunkt   : ',j1_old,j2_old    
        write(12,*) ' aktuelles Skalarprodukt : ',sp_old           
        end if

!     Directional loop
!     ----------------

      DO 1000 jdir=jstart,jstart+5
      j1_new = j1_old + j1n(jdir)
      j2_new = j2_old + j2n(jdir)
      if (j1_new.ge.1 .and. j1_new.le.n .and.     &
          j2_new.ge.1 .and. j2_new.le.n)  then
          sp_new = zx * xn(j1_new,j2_new,1,kd)    &
             + zy * xn(j1_new,j2_new,2,kd)        &
             + zz * xn(j1_new,j2_new,3,kd)

        if (loc_debug.and.jloop.ge.65) then
        write(12,*) ' Versuch in Richtung:', j1_new,j2_new,   &
                    '            SP=',sp_new
        end if

       if (sp_new.gt.sp_old) then      ! improvement !
       j1_old = j1_new
       j2_old = j2_new
       sp_old = sp_new
       jstart = jdir-((jdir-1)/6) * 6  ! try same direction first
       go to 1001
       end if                          ! improvement
      end if                           ! domain boundaries
 1000 continue                         ! loop over directions

!      If this part of the code is reached, no improvement could be
!      achieved in any direction, i.e. the fit is optimal
!      ------------------------------------------------------------

        if (loc_debug) then
        write(*,*) jloop,'.Durchlauf der Suchschleife in RWALK '
        write(*,*) '-------------------------------------------'
        write(*,*) ' Bester  Gitterpunkt   : ',j1_old,j2_old    
        write(*,*) ' Skalarprodukt         : ',sp_old           
        end if

       go to 2001


 1001 continue            ! branch adress for new, improved point
        if (loc_debug) then
        write(*,*) ' +  +  +  +  +  +  +  +  +  +  +  +  +  +  '
        write(*,*) ' neuer     Gitterpunkt   : ',j1_old,j2_old    
        write(*,*) ' neues     Skalarprodukt : ',sp_old           
        write(*,*) ' beste     Richtung      : ',jstart           
        write(*,*) ' +  +  +  +  +  +  +  +  +  +  +  +  +  +  '
        write(*,*) '                                           '
        end if
 2000 continue            ! end of walking loop           

!      If this part of the code is reached, the walk was not completed
!      whithin the walking loop 
!      ----------------------------------------------------------------

!       if (loc_debug) then
!       write(12,*) '********************************************'
!       write(12,*) 'kein optimaler Punkt wurde gefunden        '
!       write(12,*) ' aktueller Gitterpunkt   : ',j1_old,j2_old    
!       write(12,*) ' aktueller Skalarprodukt : ',sp_old           
!       end if
        stop 'kaese'

 2001 continue            ! best fit branch address      


!     Store final result in dummy arguments
!     -------------------------------------
      jj1_1 = j1_old       
      jj2_1 = j2_old       
      sp_1  = sp_old       

!       if (loc_debug) then
!       write(12,*) '     Indizes des naechsten Punktes'
!       write(12,*) ' ====================================='
!       write(12,*) '        Ecke      j1         j2      SP'
!       print *,    ' ====================================='
!       print *,    '        Ecke      j1         j2      SP'
!       write(12,*) '           ',j1_old,j2_old,sp_old
!       print *   , '           ',j1_old,j2_old,sp_old
!       end if

END SUBROUTINE zu_fuss

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE setup_xd (ig1s , ig1e , ig2s , ig2e , ids, ide, ni,             &
                     igg1s, igg1e, igg2s, igg2e,                           &
                     nproc1, nproc2, myproc, icomm, ilim1, ilim2,          &
                     yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!=======================================================================

! Parameterlist
  INTEGER (KIND=iintegers), INTENT(IN) ::    &
    ig1s , ig1e , ig2s , ig2e ,  & ! size of local field
    igg1s, igg1e, igg2s, igg2e,  & ! size of global field
    ids, ide, ni                   ! vertical size of fields

  INTEGER (KIND=iintegers), INTENT(IN) ::    &
    nproc1, nproc2, myproc, icomm  ! parallel organizational variables

  INTEGER (KIND=iintegers), INTENT(IN) ::    &
    ilim1(0:nproc1),             & !
    ilim2(0:nproc2)

  CHARACTER (LEN=25),       INTENT(OUT)  ::       &
    yerrmsg              ! for MPI error message

  INTEGER (KIND=iintegers), INTENT(OUT)::    &
    ierror                         ! error code

! Local variables
  INTEGER (KIND=iintegers)   ::     &
    marr(4,igg1s-2:igg1e+2,igg2s-2:igg2e+2,ids:ide)

  INTEGER (KIND=iintegers)   ::     &
    jd, jp1, jp2, j1, j2, jp, jb, jr, j,     &
    mbsize, mb, mi1sc, mi1ec, mi2sc, mi2ec, izmplcode

  INTEGER (KIND=iintegers)   ::     &
    mpi_status(MPI_STATUS_SIZE)

  INTEGER (KIND=iintegers)   ::     &
    idx_bound(3,(ig1e-ig1s+5)*(ig2e-ig2s+5) )

!=======================================================================

! Safety first
  IF (nproc1*nproc2 > maxprocs) THEN
    ierror  = 1 
    yerrmsg = 'maxprocs too small'
    RETURN
  ENDIF

!=======================================================================
!
!     Section 1:
!     ----------
!     Set up the array marr for all complete ( = core + extended)
!     diamonds so that it can tell us later exactly for every point,
!     from which processor and which diamond at which location
!     a particular point in the extension area has to be picked.
!
!     marr(1,.,.,.) processor number
!     marr(2,.,.,.) j1-index
!     marr(3,.,.,.) j2-index
!     marr(4,.,.,.) diamond number
!
!=======================================================================

      ierror = 0
!     Set array marr to -1 (just for debugging)
      marr = -1

!     Set up marr at the core of the diamonds

      DO jd = ids, ide
        DO jp2 = 1, nproc2
          DO j2 = ilim2(jp2-1), MIN(ilim2(jp2)-1,ni)
            DO jp1 = 1, nproc1
              DO j1 = MAX(ilim1(jp1-1),1), ilim1(jp1)-1
                jp = (jp2-1)*nproc1 + jp1-1
                marr(1,j1,j2,jd) = jp
                marr(2,j1,j2,jd) = j1
                marr(3,j1,j2,jd) = j2
                marr(4,j1,j2,jd) = jd
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!     Set  the pole

      DO jd = ids, ide
        marr(1,0,1,jd) = 0 ! Pole is allways owned by proc 0
        marr(2,0,1,jd) = 0
        marr(3,0,1,jd) = 1
        IF(jd <= 5) THEN
          marr(4,0,1,jd) = 1 ! Use allways diamond 1 for north pole
        ELSE
          marr(4,0,1,jd) = 6 ! Use allways diamond 6 for south pole
        ENDIF
      ENDDO

!     Set boundaries

      DO jd = ids, ide

!       poleward west
        DO j= 1, ni
          marr(:,j  , 0,jd) = marr(:,1,j,mpw(jd))
          marr(:,j+1,-1,jd) = marr(:,2,j,mpw(jd))
        ENDDO

!       poleward east
        DO j= 1, ni
          marr(:, 0,j+1,jd) = marr(:,j,1,mpe(jd))
          marr(:,-1,j+2,jd) = marr(:,j,2,mpe(jd))
          marr(:,-2,j+3,jd) = marr(:,j,3,mpe(jd))
        ENDDO

!       anti-poleward west
        DO j= 1, ni
          marr(:,ni+1,ni-j+1,jd) = marr(:,j,ni  ,maw(jd))
          marr(:,ni+2,ni-j+1,jd) = marr(:,j,ni-1,maw(jd))
        ENDDO

!       anti-poleward east
        DO j= 1, ni
          marr(:,ni-j+1,ni+1,jd) = marr(:,ni  ,j,mae(jd))
          marr(:,ni-j+1,ni+2,jd) = marr(:,ni-1,j,mae(jd))
          marr(:,ni-j+1,ni+3,jd) = marr(:,ni-2,j,mae(jd))
        ENDDO

      ENDDO

!     Set special points

      DO jd = ids, ide

!       Pole

        marr(:,-1, 2,jd) = marr(:,1,1,mpe(mpe(jd)))
        marr(:,-2, 3,jd) = marr(:,2,1,mpe(mpe(jd)))
        marr(:,-2, 2,jd) = marr(:,1,2,mpe(mpe(jd)))

        marr(:, 0, 0,jd) = marr(:,1,1,mpw(mpw(jd)))
        marr(:, 0,-1,jd) = marr(:,2,1,mpw(mpw(jd)))
        marr(:, 1,-1,jd) = marr(:,1,2,mpw(mpw(jd)))

        marr(:,-1, 0,jd) = marr(:, 0, 0,jd) ! Undefined
        marr(:,-1, 1,jd) = marr(:, 0, 0,jd) ! Mirror

        marr(:,-1,-1,jd) = marr(:, 0,-1,jd) ! Undefined
        marr(:,-2,-1,jd) = marr(:, 0,-1,jd) ! Undefined
        marr(:,-2, 0,jd) = marr(:, 0,-1,jd) ! Undefined
        marr(:,-2, 1,jd) = marr(:, 0,-1,jd) ! Mirror

!       West Corner

        marr(:,ni+1, 0,jd) = marr(:,ni,0,jd) ! Mirror

        marr(:,ni+2,-1,jd) = marr(:,ni+1,-1,jd) ! Undefined
        marr(:,ni+2, 0,jd) = marr(:,ni+1,-1,jd) ! Mirror

!       Antipole

        marr(:,ni+1,ni+1,jd) = marr(:,ni  ,ni+2,jd) ! Mirror
        marr(:,ni+1,ni+2,jd) = marr(:,ni  ,ni+2,jd) ! Undefined

        marr(:,ni+2,ni+1,jd) = marr(:,ni  ,ni+3,jd) ! Mirror
        marr(:,ni+2,ni+2,jd) = marr(:,ni  ,ni+3,jd) ! Undefined
        marr(:,ni+2,ni+3,jd) = marr(:,ni  ,ni+3,jd) ! Undefined
        marr(:,ni+1,ni+3,jd) = marr(:,ni  ,ni+3,jd) ! Undefined

!       East Corner

        marr(:, 0,ni+2,jd) = marr(:,-1,ni+2,jd) ! Copy
        marr(:,-1,ni+2,jd) = marr(:,-1,ni+1,jd) ! Mirror

        marr(:, 0,ni+3,jd) = marr(:,-2,ni+3,jd) ! Copy
        marr(:,-2,ni+3,jd) = marr(:,-2,ni+2,jd) ! Undefined
        marr(:,-1,ni+3,jd) = marr(:,-2,ni+2,jd) ! Mirror

      ENDDO

!=======================================================================
!
!     Section 2:
!     ----------
!     Gather all our boundary points we need during extension
!     in array idx_bound.
!
!     This is done only for 1 diamond since all other diamonds
!     are identical
!
!     idx_bound(1,.) = j1-index of boundary point
!     idx_bound(2,.) = j2-index of boundary point
!     idx_bound(3,.) == 1 if this point is needed in a 1 line exchange
!     idx_bound(3,.) == 2 if this point is needed in a 2 line exchange
!     idx_bound(3,.) == 3 these points are not really needed
!                         they are exchanged in a 2 line exchange
!                         for compatibility purposes only
!
!=======================================================================

!     Set the computational boundaries

      mi1sc = MAX(ig1s,1)
      mi1ec = ig1e
      mi2sc = ig2s
      mi2ec = MIN(ig2e,ni)

      mb = 0

      DO j1 = ig1s-2, ig1e+2
        DO j2 = ig2s-2, ig2e+2

          IF (j1<mi1sc .OR. j1>mi1ec .OR. j2<mi2sc .OR. j2>mi2ec) THEN

            ! This is a boundary point

            mb = mb+1
            idx_bound(1,mb) = j1
            idx_bound(2,mb) = j2

            IF      (j1<mi1sc-2 .OR. j1>mi1ec+2 .OR.       &
                     j2<mi2sc-2 .OR. j2>mi2ec+2) THEN
              idx_bound(3,mb) = 3
            ELSE IF (j1<mi1sc-1 .OR. j1>mi1ec+1 .OR.       &
                     j2<mi2sc-1 .OR. j2>mi2ec+1) THEN
              idx_bound(3,mb) = 2
            ELSE
              idx_bound(3,mb) = 1
            ENDIF

          ENDIF

        ENDDO
      ENDDO

      mbsize = mb

!     The processor which owns the pole (always 0) needs
!     some special treatment since it needs some more outer
!     points in order to calculate the stencils at the pole
!     correctly

      IF (myproc == 0) THEN
        DO j = 1, mbsize
          IF ( idx_bound(1,j) == -1 .AND. idx_bound(2,j) <= 2 )   &
            idx_bound(3,j) = 1 ! was previously set to 2
          IF ( idx_bound(1,j) == -2 .AND. idx_bound(2,j) <= 3 )   &
            idx_bound(3,j) = 2 ! was previously set to 3
        ENDDO
      ENDIF

!=======================================================================
!
!     Section 3:
!     ----------
!     Figure out which of our boundary points we have to receive from
!     which processor, sort them accordingly by processor and
!     set the receive arrays.
!
!=======================================================================

!     Allocate idx_recv and idx_help

      ALLOCATE ( idx_recv(3,10*mbsize) )
      ALLOCATE ( idx_help(3,10*mbsize) )

!     Sort boundary points by processor to receive them from

      mb = 0

      DO jp = 0, nproc1*nproc2-1   ! Loop over all processors
!
        np_recv_s(jp) = mb+1
!
        DO jr = 1, 2 ! Loop over the 2 boundary line cases
!
          DO jd = ids, ide ! Loop over the diamonds
!
            DO jb = 1, mbsize
!
              j1 = idx_bound(1,jb)
              j2 = idx_bound(2,jb)
!
              IF (marr(1,j1,j2,jd) /= jp) CYCLE
!
              IF ( (jr==1 .AND. idx_bound(3,jb)==1) .OR.   &
                   (jr==2 .AND. idx_bound(3,jb)/=1) ) THEN
                mb = mb+1
                idx_recv(1,mb) = j1
                idx_recv(2,mb) = j2
                idx_recv(3,mb) = jd
!
                idx_help(1,mb) = marr(2,j1,j2,jd)
                idx_help(2,mb) = marr(3,j1,j2,jd)
                idx_help(3,mb) = marr(4,j1,j2,jd)
              ENDIF
            ENDDO
          ENDDO

          IF (jr==1) THEN
            np_recv_1(jp) = mb - np_recv_s(jp) + 1
          ELSE
            np_recv_2(jp) = mb - np_recv_s(jp) + 1
          ENDIF

        ENDDO
      ENDDO

      np_recv_tot = mb

!     Internal check:

      IF (np_recv_tot /= 10*mbsize) THEN
        ierror  = 2 
        yerrmsg = 'setup_xd: Internal error'
        RETURN
      ENDIF

!=======================================================================
!
!     Section 4:
!     ----------
!     Now every processor knows, which and how many boundary points
!     it has to receive from which other processor.
!     It does not yet know, however, how many and which interior points
!     it has to send to others.
!
!     To figure that out, every processor first sends to all others
!     how many points it needs and then sends the indices of the points
!
!=======================================================================

!     Send how many points we need

      IF ( nproc1*nproc2 > 1) THEN ! no communication if only 1 processor

        DO jp = 0, nproc1*nproc2-1   ! Loop over all processors

          ! np_recv_x array -> np_send_x array-element of
          !                    the corresponding processor

          CALL MPI_Scatter(np_recv_1, 1, MPI_INTEGER, np_send_1(jp), 1,  &
                           MPI_INTEGER, jp, icomm, izmplcode)
          CALL MPI_Scatter(np_recv_2, 1, MPI_INTEGER, np_send_2(jp), 1,  &
                           MPI_INTEGER, jp, icomm, izmplcode)

          ! Provided only for the shmem communication:
          ! Tell the others where we expect our data
          ! At the moment, there is no shmem implemented in xd_p

          CALL MPI_Scatter(np_recv_s, 1, MPI_INTEGER, np_recv_t(jp), 1,  &
                           MPI_INTEGER, jp, icomm, izmplcode)

        ENDDO

      ELSE

        np_send_1(0) = np_recv_1(0)
        np_send_2(0) = np_recv_2(0)
        np_recv_t(0) = np_recv_s(0)

      ENDIF

!     Setup np_send_s, count total number of points to send

      mb = 0

      DO jp = 0, nproc1*nproc2-1
        np_send_s(jp) = mb + 1
        mb = mb + np_send_2(jp)
      ENDDO

      np_send_tot = mb

!     Allocate idx_send

      ALLOCATE ( idx_send(3,np_send_tot) )

!     Get the indices of the points we have to send

      DO jp = 0, nproc1*nproc2-1   ! Loop over all processors

        IF (jp == myproc) THEN

!         It is our turn to send which points we need

          DO jp1 = 0, nproc1*nproc2-1

            IF ( jp1 == myproc ) THEN

              ! Don't send to ourself, just copy

              DO jb = 1, np_recv_2(myproc)
                j1 = np_send_s(myproc) + jb - 1
                j2 = np_recv_s(myproc) + jb - 1
                idx_send(1,j1) = idx_help(1,j2)
                idx_send(2,j1) = idx_help(2,j2)
                idx_send(3,j1) = idx_help(3,j2)
              ENDDO

            ELSE

              ! Send what we want to receive

              IF (np_recv_2(jp1) > 0)                                 &
                CALL MPI_Send(idx_help(1,np_recv_s(jp1)),             &
                              3*np_recv_2(jp1), MPI_INTEGER, jp1, 3,  &
                              icomm, izmplcode)

            ENDIF

          ENDDO

        ELSE

!         It is our turn to receive which points are needed from
!         processor jp

          IF (np_send_2(jp) > 0)                                  &
            CALL MPI_Recv(idx_send(1,np_send_s(jp)),              &
                          3*np_send_2(jp), MPI_INTEGER, jp, 3,    &
                          icomm, mpi_status, izmplcode)

        ENDIF

        IF (nproc1*nproc2 > 1) THEN
          CALL MPI_Barrier(icomm, izmplcode)
        ENDIF
      ENDDO

!=======================================================================
!
!     Section 5:
!     ----------
!     Final work
!
!=======================================================================

!     We are all done, deallocate idx_help

      DEALLOCATE ( idx_help )

!     Calculate maximum of points to send/receive over all processors

      IF (nproc1*nproc2 > 1) THEN
        CALL MPI_Allreduce(np_recv_tot, np_recv_max, 1, MPI_INTEGER,  &
                           MPI_MAX, icomm, izmplcode)
        CALL MPI_Allreduce(np_send_tot, np_send_max, 1, MPI_INTEGER,  &
                           MPI_MAX, icomm, izmplcode)
      ELSE
        np_recv_max = np_recv_tot
        np_send_max = np_send_tot
      ENDIF

END SUBROUTINE setup_xd

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE xd (p   , kip1s , kip1e , kip2s , kip2e , ki3s , ki3e ,       &
               knt , kits  , kite  ,                                     &
               knpd, kig1s , kig1e , kig2s , kig2e , kids , kide ,       &
               kx  , lcopy , loc_debug,                                  &
               px  , kipx1s, kipx1e, kipx2s, kipx2e, knpxd, undef, kierr)       

!------------------------------------------------------------------------------
!
! Description:
!     *xd* extends the globally defined array "px" for the diamonds    
!     "kids" up to "kide" by "kx" (1 or 2) rows and columns taken from
!     the corresponding values of the appropriate neighbouring diamonds.
!     The input array "p" and the output one "px" may be identical.
!     The extension is performed for all levels/layers from "ki3s" until
!     "ki3e" and for the time levels "kits" to "kite".
!
!     Additionally, *xd* sets the values at the western boundary of the
!     diamonds to the corresponding values of the western neighbouring
!     ones. This step is necesssary since this common row/column is
!     computed independently in both diamonds. Due to round-off errors
!     the results may differ slightly which may lead to nonlinear
!     instabilities during the forecast.
!
! Method:
!
!==============================================================================

! Input
INTEGER (KIND=iintegers), INTENT(IN)  ::         &
  kip1s ,  & ! first  dimension of array "p",  start index
  kip1e ,  & ! first  dimension of array "p",  end   index
  kip2s ,  & ! second dimension of array "p",  start index
  kip2e ,  & ! second dimension of array "p",  end   index
  kipx1s,  & ! first  dimension of array "px", start index
  kipx1e,  & ! first  dimension of array "px", end   index
  kipx2s,  & ! second dimension of array "px", start index
  kipx2e,  & ! second dimension of array "px", end   index
  ki3s  ,  & ! third  dimension (# of layers/levels), start
  ki3e  ,  & ! third  dimension (# of layers/levels), end  
  knt   ,  & ! number of time levels of array "p"
  kits  ,  & ! first time levels of array "p" to be extended
  kite  ,  & ! last  time levels of array "p" to be extended
  knpd  ,  & ! number of diamonds of array "p"; must be 10!
  knpxd ,  & ! number of diamonds of array "px"
  kig1s ,  & ! first  dimension of diamond core, start index
  kig1e ,  & ! first  dimension of diamond core, end   index
  kig2s ,  & ! second dimension of diamond core, start index
  kig2e ,  & ! second dimension of diamond core, end   index
  kids  ,  & ! first  diamond to be extended
  kide  ,  & ! last   diamond to be extended
  kx         ! number of rows/columns for extension (1 or 2)

INTEGER (KIND=iintegers), INTENT(OUT) ::         &
  kierr      ! error flag, kierr = 0 if no error occured

REAL (KIND=ireals), INTENT(IN)     ::            &
  undef,   & ! for checking undefined variables
  p (kip1s :kip1e,  kip2s :kip2e,  ki3s:ki3e, knt, knpd)
             ! Global array "p", must cover all 10 diamonds

LOGICAL lcopy     ! copy  switch; if .true. the input and output
                  ! arrays are the same
LOGICAL loc_debug ! debug switch; if .true. print information

!=======================================================================
!
!     Output
REAL (KIND=ireals), INTENT(INOUT)  ::            &
  px(kipx1s:kipx1e, kipx2s:kipx2e, ki3s:ki3e, knt, knpxd)
!
!=======================================================================
!
! Local variables
! Array dimensions with 1 or 2 extension rows/columns
INTEGER (KIND=iintegers)           ::            &
  mi1sm1,  & ! = kig1s - 1 start index of extended diamond
  mi1ep1,  & ! = kig1e + 1 end   index of extended diamond
  mi2sm1,  & ! = kig2s - 1 start index of extended diamond
  mi2ep1,  & ! = kig2e + 1 end   index of extended diamond
  mi1sm2,  & ! = kig1s - 2 start index of extended diamond
  mi1ep2,  & ! = kig1e + 2 end   index of extended diamond
  mi2sm2,  & ! = kig2s - 2 start index of extended diamond
  mi2ep2     ! = kig2e + 2 end   index of extended diamond

! Index of neighbouring diamonds
INTEGER (KIND=iintegers)           ::            &
  mns  ,   & ! N/S-hemishere discriminator
  mpe  ,   & ! poleward eastern-neighbour
  mpw  ,   & ! poleward western-neighbour
  mae  ,   & ! anti-poleward eastern-neighbour
  maw  ,   & ! anti-poleward western-neighbour
  mpp  ,   & ! across pole point 'neighbour'
  mppm1      ! across pole point 'neighbour - 1'
 
! Loop indices          
INTEGER (KIND=iintegers)           :: j, j1, j2, j3, jt, jd

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine xd
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Check input variables and Initializations
!------------------------------------------------------------------------------

  kierr = 0
  IF ((kx .NE. 1) .AND. (kx .NE. 2)) THEN
    PRINT *,'  Error in subroutine *xd*, kx= ', kx, '  is not 1 or 2'
    kierr = -1
    RETURN
  ENDIF
  IF ((kip1s .GT. kig1s) .OR. (kip1e .LT. kig1e) .OR.             &
      (kip2s .GT. kig2s) .OR. (kip2e .LT. kig2e)) THEN
    PRINT *,'  Error in subroutine *xd*, dimensions of input',    &
            ' array "p" are not covering the diamond core'
    PRINT *,'  kig1s: ', kig1s,'  kip1s: ', kip1s
    PRINT *,'  kig1e: ', kig1e,'  kip1e: ', kip1e
    PRINT *,'  kig2s: ', kig2s,'  kip2s: ', kip2s
    PRINT *,'  kig2e: ', kig2e,'  kip2e: ', kip2e
    kierr = -2
    RETURN
  ENDIF
  IF (((kide-kids+1) .GT. knpxd) .OR. (knpd .NE. 10)) THEN
    PRINT *,'  Error in subroutine *xd*, number of diamonds',     &
            ' provided is insufficient'
    PRINT *,'  kids: ', kids,'  kide: ', kide 
    PRINT *,'  knpd: ', knpd
    kierr = -3
    RETURN
  ENDIF

  ! Set local array dimensions
  mi1sm1 = kig1s - 1
  mi1ep1 = kig1e + 1
  mi1sm2 = kig1s - 2
  mi1ep2 = kig1e + 2
  mi2sm1 = kig2s - 1
  mi2ep1 = kig2e + 1
  mi2sm2 = kig2s - 2
  mi2ep2 = kig2e + 2

  ! Loop over the diamonds from "kids" to "kide"
  DO jd = kids, kide

    ! Determination of appropriate neighbouring diamonds
    mns   = (jd-1)/5        
    mpe   =  jd + 1 - (jd/(5*(1+mns)))*5                   
    mpw   =  jd - 1 + ((mns*10+6-jd)/(5*(1+mns)))*5       
    mae   =  jd + 5 - 9*mns-5*(jd/10)                       
    maw   =  jd + 4 + ((6-jd)/5)*5-9*mns                  
    mpp   =  jd + 3 - ((jd+2)/5)*5+5*mns                  
    mppm1 =  jd + 2 - ((jd+1)/5)*5+5*mns                  

    ! Loop over the time levels
    DO jt = kits, kite

      ! Loop over the layers/levels
      DO j3 = ki3s, ki3e

        ! Generate extended array
        ! 1. Core of selected diamond (not necessary if the input and
        ! output arrays are the same)

        IF (.NOT.lcopy) THEN
          DO j2=kig2s,kig2e
            DO j1=kig1s,kig1e
              px (j1,j2,j3,jt,jd) = p (j1,j2,j3,jt,jd)
            ENDDO
          ENDDO
        ENDIF

!------------------------------------------------------------------------------
! Section 2: Since one row and column is the same in neighbouring diamonds
!            but are computed independently, the values of the western  
!            diamond is copied into the corresponding place of the eastern
!            neighbour to prevent instabilities which are due to round-off
!            errors.
!
!            The exchange does not include both pole points (north and
!            south) since the local coordinate system is different for each
!            diamond at the poles, thus the wind componentes have to be
!            rotated. The exchange of the pole points for the prognostic
!            variables is performed in *spreadpole*.
!------------------------------------------------------------------------------

        IF ( knpxd .EQ. 10 ) THEN
          DO j = kig2s+1,kig2e
            IF (px(kig1s, j, j3, jt, jd) /= undef) THEN
              px(j-1, kig2s, j3, jt, mpe) = px(kig1s, j, j3, jt, jd )
            ENDIF
          ENDDO

          DO j = kig1s,kig1e
            IF (px(j, kig2e, j3, jt, jd) /= undef) THEN
              px(kig1e, kig2e-j, j3, jt, mae) = px(j, kig2e, j3, jt, jd )
            ENDIF
          ENDDO
        ENDIF
      ENDDO     ! End of external loop over layers/levels
    ENDDO     ! End of external loop over time levels
  ENDDO       ! End of external loop over the diamonds


  ! Loop over the diamonds from "kids" to "kide"
  DO jd = kids, kide

    ! Determination of appropriate neighbouring diamonds
    mns   = (jd-1)/5        
    mpe   =  jd + 1 - (jd/(5*(1+mns)))*5                   
    mpw   =  jd - 1 + ((mns*10+6-jd)/(5*(1+mns)))*5       
    mae   =  jd + 5 - 9*mns-5*(jd/10)                       
    maw   =  jd + 4 + ((6-jd)/5)*5-9*mns                  
    mpp   =  jd + 3 - ((jd+2)/5)*5+5*mns                  
    mppm1 =  jd + 2 - ((jd+1)/5)*5+5*mns                  

    ! Loop over the time levels
    DO jt = kits, kite

      ! Loop over the layers/levels
      DO j3 = ki3s, ki3e  ! Loop in parallel over all layers

!------------------------------------------------------------------------------
! Section 3: Extension by 1 row/column 
!            Use data from adjacent neighbours
!------------------------------------------------------------------------------

        DO j=kig2s,kig1e
          px(j     , mi2sm1, j3, jt, jd) = p(kig1s+1  , j        , j3, jt, mpw)
          px(mi1sm1, j+1   , j3, jt, jd) = p(j-1      , kig2s+1  , j3, jt, mpe)
          px(mi1ep1, j     , j3, jt, jd) = p(kig1e+1-j, kig2e-1  , j3, jt, maw)
          px(j-1   , mi2ep1, j3, jt, jd) = p(kig1e-1  , kig2e+1-j, j3, jt, mae)
        ENDDO

        ! Special points at corner of selected diamond (note the 
        ! identical array elements)
        px(mi1sm1, kig2s , j3, jt, jd) =  p (kig1s+1, kig2s , j3, jt, mpp)
        px(kig1s , mi2sm1, j3, jt, jd) =  p (kig1s+1, kig2s , j3, jt, mpp)
        px(mi1ep1, mi2sm1, j3, jt, jd) =  px(kig1e  , mi2sm1, j3, jt, jd)
        px(mi1sm1, mi2ep1, j3, jt, jd) =  px(mi1sm1 , kig2e , j3, jt, jd)
        px(mi1ep1, kig2e , j3, jt, jd) =  p (kig1e-1, kig2s , j3, jt, mae)
        px(kig1e , mi2ep1, j3, jt, jd) =  p (kig1e-1, kig2s , j3, jt, mae)

        ! The following two corner points are undefined due to the penta-
        ! gonal structure of the pole points. To ease the vectorization,
        ! these points of extended array will be copied from neighbours
        px(mi1sm1, mi2sm1, j3, jt, jd) =  px(kig1s  , mi2sm1, j3, jt, jd)
        px(mi1ep1, mi2ep1, j3, jt, jd) =  px(kig1e+1, kig2e , j3, jt, jd)

!------------------------------------------------------------------------------
! Section 4: Second extension line, if requested
!------------------------------------------------------------------------------

        IF (kx .EQ. 2) THEN
          DO j=kig2s,kig1e
            px(j+1   , mi2sm2, j3,jt,jd) = p (kig1s+2  , j        , j3,jt,mpw)
            px(mi1sm2, j+2   , j3,jt,jd) = p (j-1      , kig2s+2  , j3,jt,mpe)
            px(mi1ep2, j     , j3,jt,jd) = p (kig1e+1-j, kig2e-2  , j3,jt,maw)
            px(j-1   , mi2ep2, j3,jt,jd) = p (kig1e-2  , kig2e+1-j, j3,jt,mae)
          ENDDO

          ! Corners of extended diamond
          ! Poleward corner (note identical array elements)
          px(mi1sm2 , kig2s  , j3,jt,jd) = p (kig1s+2, kig2s  , j3, jt, mpp)   
          px(kig1s  , mi2sm2 , j3,jt,jd) = p (kig1s+2, kig2s  , j3, jt, mpp)
          px(kig1s+1, mi2sm2 , j3,jt,jd) = p (kig1s+1, kig2s+1, j3, jt, mpp)
          px(mi1sm2 , kig2s+1, j3,jt,jd) = p (kig1s+1, kig2s+1, j3, jt, mppm1)

          ! Antipoleward corner (note identical array elements)
          px(kig1e , mi2ep2, j3, jt, jd) = p (kig1e-2, kig2s, j3, jt, mae)
          px(mi1ep2, kig2e , j3, jt, jd) = p (kig1e-2, kig2s, j3, jt, mae)

          ! Left and right corner duplication
          px(mi1ep2, mi2sm1, j3, jt, jd) = px(mi1ep1, mi2sm2, j3, jt, jd)
          px(mi1sm1, mi2ep2, j3, jt, jd) = px(mi1sm2, mi2ep1, j3, jt, jd)

!------------------------------------------------------------------------------
! Section 5. The following 8 corner points are undefined due to the penta-
!            gonal structure of the pole points. To ease the vectorization,
!            these points of extended array will be copied from neighbours
!------------------------------------------------------------------------------

          ! Set the remaining matrix elements for kx = 2
          px(mi1sm2, mi2sm2, j3, jt, jd) = px(mi1sm2, kig2s , j3, jt, jd)
          px(mi1sm1, mi2sm2, j3, jt, jd) = px(kig1s , mi2sm2, j3, jt, jd)
          px(mi1sm2, mi2sm1, j3, jt, jd) = px(mi1sm2, kig2s , j3, jt, jd)
          px(mi1ep2, mi2sm2, j3, jt, jd) = px(mi1ep1, mi2sm2, j3, jt, jd)
          px(mi1ep2, mi2ep1, j3, jt, jd) = px(mi1ep2, kig2e , j3, jt, jd)
          px(mi1sm2, mi2ep2, j3, jt, jd) = px(mi1sm2, mi2ep1, j3, jt, jd)
          px(mi1ep1, mi2ep2, j3, jt, jd) = px(kig1e , mi2ep2, j3, jt, jd)
          px(mi1ep2, mi2ep2, j3, jt, jd) = px(mi1ep2, kig2e , j3, jt, jd)
        ENDIF   ! End of extension with 2 rows/columns

      ENDDO     ! End of external loop over layers/levels
    ENDDO     ! End of external loop over time levels
  ENDDO     ! End of external loop over the diamonds

!=======================================================================
!
!     Attention:
!     For the extension by 1 row/column (kx = 1), there are four
!     identical (mirrored) and two undefined array elements, namely
!
!     Identical array elements:
!     (mi1sm1, kig2s ) = (kig1s , mi2sm1)
!     (mi1ep1, kig2e ) = (kig1e , mi2ep1)
!     (mi1ep1, mi2sm1) = (kig1e , mi2sm1)
!     (mi1sm1, mi2ep1) = (mi1sm1, kig2e )
!
!     Undefined array elements (taken from the neighbourhood to ease
!     vectorization)
!     (mi1sm1, mi2sm1)
!     (mi1ep1, mi2ep1)
!
!     For the extension by 2 rows/columns (kx = 2), there are four
!     identical (mirrored) and eight undefined array elements, additio-
!     nally to the ones for kx = 1, namely
!
!     Identical array elements:
!     (mi1sm2, kig2s ) = (kig1s , mi2sm2)
!     (kig1e , mi2ep2) = (mi1ep2, kig2e )
!     (mi1ep2, mi2sm1) = (mi1ep1, mi2sm2)
!     (mi1sm1, mi2ep2) = (mi1sm2, mi2ep1)
!
!     Undefined array elements (taken from the neighbourhood to ease
!     vectorization)
!     (mi1sm2, mi2sm2)
!     (mi1sm1, mi2sm2)
!     (mi1sm2, mi2sm1)
!     (mi1ep2, mi2sm2)
!
!     (mi1ep2, mi2ep1)
!     (mi1ep1, mi2ep2)
!     (mi1ep2, mi2ep2)
!     (mi1sm2, mi2ep2)
!
!
!=======================================================================

END SUBROUTINE xd

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE xd_p(p, kip1s, kip1e, kip2s, kip2e, kip3s, kip3e, kt, kr,   &
                myproc, nproc, icomm, imp_reals)

!------------------------------------------------------------------------------
!
! Description:
!   *xd_p* extends the subdomains of a globally defined array "p" for all
!   diamonds by "kr" (1 or 2) rows and columns taken from
!   the corresponding values of the appropriate neighbouring subdomains.
!   The input array "p" and the output one "px" may be identical.
!   The extension is performed for all levels/layers from "kip3s" until
!   "kip3e".
!
! Method:
!
!=======================================================================

! Parameterlist
  INTEGER (KIND=iintegers), INTENT(IN)     ::    &
           kip1s, kip1e, kip2s, kip2e, kip3s, kip3e, kt, kr,    &
           nproc, myproc, imp_reals, icomm

  REAL    (KIND=ireals)   , INTENT(INOUT)  ::    &
           p(kip1s:kip1e, kip2s:kip2e, kip3s:kip3e, kt, 10)

! Local variables
  REAL    (KIND=ireals)      ::    &
    zbuff_send(kip3s:kip3e,np_send_max),   & !
    zbuff_recv(kip3s:kip3e,np_recv_max)

  INTEGER (KIND=iintegers)   ::    &
    j1, j2, j3, jd, jp, jp1, j,    &
    mpr, mpl, mreq, izmplcode

  INTEGER (KIND=iintegers)   ::    &
    mpi_req(maxprocs), mpi_statuses(MPI_STATUS_SIZE,maxprocs)

!=======================================================================

      mpl = kip3e - kip3s + 1 ! number of layers

!     Put all points we have to send into the send buffer

      DO j3 = kip3s, kip3e
        DO jp = 0, nproc-1

          IF ( kr==1 ) mpr =  np_send_1(jp)
          IF ( kr==2 ) mpr =  np_send_2(jp)

          DO j = np_send_s(jp), np_send_s(jp) + mpr - 1

            j1 = idx_send(1,j)
            j2 = idx_send(2,j)
            jd = idx_send(3,j)

            zbuff_send(j3,j) = p(j1,j2,j3,1,jd)

          ENDDO
        ENDDO
      ENDDO

!     Setup non blocking receives on all processors from which
!     we expect data

      mreq = 0

      DO jp = 0, nproc-1

        IF ( jp==myproc ) CYCLE ! Don't receive from ourself

        IF ( kr==1 ) mpr =  np_recv_1(jp)
        IF ( kr==2 ) mpr =  np_recv_2(jp)

        IF ( mpr>0 ) THEN
          mreq = mreq+1
          CALL MPI_Irecv( zbuff_recv(kip3s,np_recv_s(jp)), mpr*mpl, &
                          imp_reals, jp, 1, icomm,          &
                          mpi_req(mreq), izmplcode)
        ENDIF

      ENDDO

!     Send data to all processors who want to have data from us

      DO jp1 = 0, nproc-1

        jp = MOD( jp1+myproc, nproc )

        IF ( kr==1 ) mpr =  np_send_1(jp)
        IF ( kr==2 ) mpr =  np_send_2(jp)

        IF ( jp==myproc ) THEN

!         Don't send to ourself, just copy

          DO j3 = kip3s, kip3e
            DO j = 1, mpr
              zbuff_recv(j3,np_recv_s(jp)+j-1) =    &
                 zbuff_send(j3,np_send_s(jp)+j-1)
            ENDDO
          ENDDO

        ELSE

          IF ( mpr>0 ) THEN
            CALL MPI_Send( zbuff_send(kip3s,np_send_s(jp)), mpr*mpl,  &
                            imp_reals, jp, 1, icomm,          &
                            izmplcode)
          ENDIF

        ENDIF

      ENDDO

!
!     Wait until all receives have finished
!
      IF (mreq>0) THEN
        CALL MPI_Waitall(mreq, mpi_req, mpi_statuses, izmplcode)
      ENDIF

!     Fill received data into it's location

      DO j3 = kip3s, kip3e
        DO jp = 0, nproc-1

          IF ( kr==1 ) mpr =  np_recv_1(jp)
          IF ( kr==2 ) mpr =  np_recv_2(jp)

          DO j = np_recv_s(jp), np_recv_s(jp) + mpr - 1

            j1 = idx_recv(1,j)
            j2 = idx_recv(2,j)
            jd = idx_recv(3,j)

            p(j1,j2,j3,1,jd) = zbuff_recv(j3,j)

          ENDDO
        ENDDO
      ENDDO

END SUBROUTINE xd_p

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE gather_gme_field (fld,     ig1s , ig1e , ig2s , ig2e ,          &
                             fld_tot, igg1s, igg1e, igg2s, igg2e, nb, nf,  &
                             nproc1, nproc2, myproc, my_num1, my_num2,     &
                             ilim1, ilim2, imp_reals, max_gme_core, icomm, &
                             ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Combine the partial fields fld to a total field fld_tot for GME
!
! Method:
!
!==============================================================================

! Parameterlist
  INTEGER (KIND=iintegers), INTENT(IN) ::    &
    ig1s , ig1e , ig2s , ig2e ,  & ! size of local field
    igg1s, igg1e, igg2s, igg2e,  & ! size of global field
    nb, nf                         ! variations from the size

  INTEGER (KIND=iintegers), INTENT(IN) ::    &
    nproc1, nproc2, myproc, my_num1, my_num2,  & !
    imp_reals, max_gme_core, icomm ! further organizational variables

  INTEGER (KIND=iintegers), INTENT(IN) ::    &
    ilim1(0:nproc1), ilim2(0:nproc2)

  REAL    (KIND=ireals),    INTENT(IN) ::    &
    fld(ig1s-nb:ig1e+nb,ig2s-nb:ig2e+nb,nf)

  REAL    (KIND=ireals),    INTENT(OUT)::    &
    fld_tot(igg1s-nb:igg1e+nb,igg2s-nb:igg2e+nb,nf)

  INTEGER (KIND=iintegers), INTENT(OUT)::    &
    ierror

! Local variables
  INTEGER (KIND=iintegers) :: mpi_req(nproc1*nproc2), &
                          mpi_statuses(MPI_STATUS_SIZE,nproc1*nproc2)

  REAL (KIND=ireals) :: rb( max_gme_core*nf, nproc1*nproc2) ! receive buffers

  INTEGER (KIND=iintegers) :: i1s, i1e, i2s, i2e, n, n1, n2, k, i1, i2,    &
                              count, izmplcode

!==============================================================================

! If we are not running in parallel, this is a very simple routine ...

  ierror = 0
  IF (nproc1*nproc2 == 1) THEN
    fld_tot(:,:,:) = fld(:,:,:)
    RETURN
  ENDIF

! Set up non blocking receives on all other processors
! (with the exception of ourselves)

  DO n = 0, nproc1*nproc2-1
    IF (n==myproc) CYCLE
    CALL MPI_Irecv( rb(1,n+1), max_gme_core*nf, imp_reals, n, 2, &
                    icomm, mpi_req(n+1), izmplcode)
  ENDDO

! Send our part of the field to all
! (with the exception of ourselves)
! We are sending to myproc+1, myproc+2 ....... in order to avoid
! that all procs try to send to a single one at once

  count = (ig1e-ig1s+2*nb+1)*(ig2e-ig2s+2*nb+1)*nf

  DO k = 1, nproc1*nproc2-1
    n = k + myproc
    IF ( n>nproc1*nproc2-1 ) n = n-nproc1*nproc2
    CALL MPI_Send( fld, count, imp_reals, n, 2, icomm, izmplcode)
  ENDDO

! Put our part of the field into fld_tot

  fld_tot(ig1s-nb:ig1e+nb,ig2s-nb:ig2e+nb,:) = fld

! Wait until all receives have finished

  DO n = 0, nproc1*nproc2-1
    IF (n==myproc) CYCLE
    CALL MPI_Wait( mpi_req(n+1), mpi_statuses(1,n+1), izmplcode)
  ENDDO

! Unpack the data from the receive buffers and put it into fld_glob

  DO n2 = 0, nproc2-1
    DO n1 = 0, nproc1-1

      IF ( n1==my_num1 .AND. n2==my_num2) CYCLE

      i1s = ilim1(n1)
      i1e = ilim1(n1+1)-1
      i2s = ilim2(n2)
      i2e = ilim2(n2+1)-1

!     Consistency check

      n = n2*nproc1+n1    ! Johanni used nprocx here instead of nproc1
                          ! and for him nprocx = nproc1
      CALL MPI_Get_count(  mpi_statuses(1,n+1), imp_reals, count, izmplcode)
      IF (count /= (i1e-i1s+1+2*nb)*(i2e-i2s+1+2*nb)*nf) THEN
        PRINT *,'Proc=',myproc,'n1=',n1,'n2=',n2,'received=',count, &
                'expected=',(i1e-i1s+1+2*nb)*(i2e-i2s+1+2*nb)*nf
        ierror = 1
        RETURN
      ENDIF

      n = 0

      DO k = 1, nf
        DO i2 = i2s-nb, i2e+nb
          DO i1 = i1s-nb, i1e+nb

            n = n+1

!           Take care that the boundaries are not taken if they are
!           in the inner of the core
!           If there has been a proper boundary exchange before calling
!           this routine, this shouldn't matter, but be safe!

            IF(n2>0 .AND. i2<i2s) CYCLE
            IF(n1>0 .AND. i1<i1s) CYCLE
            IF(n2<nproc2-1 .AND. i2>i2e) CYCLE
            IF(n1<nproc1-1 .AND. i1>i1e) CYCLE

            fld_tot(i1,i2,k) = rb(n,n2*nproc1+n1+1)    
                                  ! Johanni used nprocx here instead of nproc1
                                  ! and for him nprocx = nproc1

          ENDDO ! i1-Loop
        ENDDO   ! i2-Loop
      ENDDO     ! k -Loop
    ENDDO       ! n1-Loop
  ENDDO         ! n1-Loop

!==============================================================================

END SUBROUTINE gather_gme_field

!==============================================================================

SUBROUTINE init_gme_interpol                                                 &
                  (lolp_gme, field_gme, kig1sm2, kig1ep2, kig2sm2, kig2ep2,  &
                   kids, kide,                                               &
                   lolp_lm, weights_ip, nearest_ip, match_ip, log_ip,        &
                   kispoke, pbaryll, kindex, kilons, kilone, kilats, kilate, &
                   undef, idebug, yerrmsg, kierr)

!------------------------------------------------------------------------------
!
! Description:
!   init_gme_interpol initializes the horizontal interpolation of a 
!   GME field to the COSMO-model grid. For 'M'atch or 'N'earest neighbor 
!   interpolation, necessary indices are computed and stored.
!
! Method:
!
!==============================================================================
!
! Input
INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kig1sm2,  & ! first  dimension of array "px",   start index
  kig1ep2,  & ! first  dimension of array "px",   end   index
  kig2sm2,  & ! second dimension of array "px",   start index
  kig2ep2,  & ! second dimension of array "px",   end   index
  kids   ,  & ! start diamond
  kide        ! end diamond

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kilons ,  & ! first  dimension of array "pxi",  start index
  kilone ,  & ! first  dimension of array "pxi",  end   index
  kilats ,  & ! second dimension of array "pxi",  start index
  kilate      ! second dimension of array "pxi",  end   index

LOGICAL, INTENT(IN)                    ::    &
  lolp_gme (kig1sm2:kig1ep2, kig2sm2:kig2ep2,kids:kide),   & !
  lolp_lm  (kilons:kilone, kilats:kilate)

REAL    (KIND=ireals)   , INTENT(IN)   ::    &
  field_gme(kig1sm2:kig1ep2, kig2sm2:kig2ep2,kids:kide),   &
              ! to check for undefined values
  pbaryll  (kilons  :kilone  , kilats  :kilate  , 2),      &
              ! first and second barycentric coordinate of points
  undef       ! value assigned to undefined grid points

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kindex  (kilons  :kilone  , kilats  :kilate  , 4)
              ! index of the triangle containing the point: j1,j2,jd,mt

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kispoke (12)    ! Normal spokes 

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  idebug          ! for debug output

!==============================================================================

! Output
REAL    (KIND=ireals),    INTENT(OUT)   ::      &
  weights_ip (kilons:kilone, kilats:kilate,19)
            !  Weights for the interpolation, to choose, which 
            !  points have to be taken for the match interpolation

INTEGER (KIND=iintegers), INTENT(OUT)   ::      &
  nearest_ip (kilons:kilone, kilats:kilate, 3), &
            !  indices of nearest GME gridpoint for nearest neighbor interpolation
  match_ip   (kilons:kilone, kilats:kilate, 3)
            !  indices of nearest GME gridpoint for match interpolation

LOGICAL                 , INTENT(OUT)   ::      &
  log_ip     (kilons:kilone, kilats:kilate)
            !  to use a far away GME gridpoint with same lsm

CHARACTER (LEN=*),        INTENT(OUT)   ::      &
  yerrmsg   ! for error message

INTEGER (KIND=iintegers), INTENT(OUT)   ::      &
  kierr     ! error flag, kierr = 0 if no error occured

!==============================================================================

! Local variables
REAL (KIND=ireals)         ::      &
  za , zb , zc    ! Barycentric coordinates

INTEGER (KIND=iintegers)   ::      &
  l1, l2,             & ! Loop indices
  j1, j2, jd, m1, m2, & ! Gridpoint indices
  mindist, ndist,     & !
  i,k1,k2               ! Loop indices

INTEGER (KIND = iintegers) ::      &
  ipx,                & ! number of surrounding values
  i_c, i_1, i_2         ! index of GME grid points in the 19-point numbering

LOGICAL                    ::      &
  lolp3_gme(3)       ! LSM of the three surrounding grid points

!==============================================================================

  kierr   = 0
  yerrmsg(:) = ' '

  weights_ip  (:,:,:)  = 0.0_ireals
  match_ip    (:,:,:)  = -10_iintegers
  nearest_ip  (:,:,:)  = -10_iintegers
  log_ip      (:,:)    = .FALSE.

!==============================================================================

  DO l2 = kilats, kilate
    DO l1 = kilons, kilone

      !------------------------------------------------------------------------
      ! Section 1: Initializations
      !------------------------------------------------------------------------

      j1    = kindex(l1,l2,1) ! first gme grid dimension
      j2    = kindex(l1,l2,2) ! second gme grid dimension
      jd    = kindex(l1,l2,3) ! number of diamond
      m1    = kindex(l1,l2,4) ! Number of triangle for interpolation
      m2    = MOD (m1,6) + 1
      za    = pbaryll(l1,l2,1)
      zb    = pbaryll(l1,l2,2)
      zc    = 1. - za - zb

      lolp3_gme(1) = lolp_gme(j1,j2,jd)
      lolp3_gme(2) = lolp_gme(j1+kispoke(m1),j2+kispoke(m1+6),jd)
      lolp3_gme(3) = lolp_gme(j1+kispoke(m2),j2+kispoke(m2+6),jd)

      ! For the computations it is necessary to know the correspondence
      ! between the GME grid points in the usual indexing and the numbering
      ! for the 19 considered points around the centerpoint j1,j2,jd.
      ! This depends on m1:
      SELECT CASE (m1)
      CASE (1)
        i_c = 10_iintegers
        i_1 = 11_iintegers
        i_2 = 15_iintegers
      CASE (2)
        i_c = 10_iintegers
        i_1 = 15_iintegers
        i_2 = 14_iintegers
      CASE (3)
        i_c = 10_iintegers
        i_1 = 14_iintegers
        i_2 =  9_iintegers
      CASE (4)
        i_c = 10_iintegers
        i_1 =  9_iintegers
        i_2 =  5_iintegers
      CASE (5)
        i_c = 10_iintegers
        i_1 =  5_iintegers
        i_2 =  6_iintegers
      CASE (6)
        i_c = 10_iintegers
        i_1 =  6_iintegers
        i_2 = 11_iintegers
      END SELECT

      !------------------------------------------------------------------------
      ! Section 2: Nearest neighbor interpolation
      !------------------------------------------------------------------------

      ! Set up the nearest GME grid point for nearest neighbor interpolation
      ! by checking the three nearest grid points
      IF     (lolp3_gme(1) .EQV. lolp_lm(l1,l2)) THEN
        nearest_ip (l1,l2, 1) = j1
        nearest_ip (l1,l2, 2) = j2
        nearest_ip (l1,l2, 3) = jd
      ELSEIF (lolp3_gme(2) .EQV. lolp_lm(l1,l2)) THEN
        nearest_ip (l1,l2, 1) = j1 + kispoke(m1)
        nearest_ip (l1,l2, 2) = j2 + kispoke(m1+6)
        nearest_ip (l1,l2, 3) = jd
      ELSEIF (lolp3_gme(3) .EQV. lolp_lm(l1,l2)) THEN
        nearest_ip (l1,l2, 1) = j1 + kispoke(m2)
        nearest_ip (l1,l2, 2) = j2 + kispoke(m2+6)
        nearest_ip (l1,l2, 3) = jd
      ELSE
        kierr   = 1
        WRITE (yerrmsg,'(A,I4,A,I4)') 'yitype = N: no GME point for lat i = ', &
                                       l1, ' and lon j = ', l2
      ENDIF

      !------------------------------------------------------------------------
      ! Section 3: Weights for match interpolation
      !------------------------------------------------------------------------

      ! Determine the weights for the match interpolation
      ! First check the characteristics of the three surrounding points
      IF ( ALL(lolp3_gme(:) .EQV. lolp_lm(l1,l2)) ) THEN

        ! all points have the same characteristic, do the same as with 'L'
        ! but the numbering of the GME points now is according to the 19 
        ! surrounding points: where the three points are in this numbering
        ! depends on m1, m2 and was computed above:
        weights_ip (l1,l2,i_c) = 1.0_ireals - za - zb
        weights_ip (l1,l2,i_1) = za
        weights_ip (l1,l2,i_2) = zb

        IF (idebug > 50) THEN
          WRITE (*,'(A,2I5,6I6)')  'linear:  ', l1, l2, j1, j2, jd, i_c, i_1, i_2
        ENDIF
      ELSE

        ipx = 0

        ! from the 3 surrounding points take the ones with identical lsm
        IF (lolp3_gme(1) .EQV. lolp_lm(l1,l2)) THEN
          ! that is always the center point
          ipx  = ipx  + 1
          weights_ip (l1,l2,i_c) = 1.0_ireals
        ENDIF
        ! The next points depend on m1, their indices have been computed above
        IF (lolp3_gme(2) .EQV. lolp_lm(l1,l2)) THEN
          ipx  = ipx  + 1
          weights_ip (l1,l2,i_1) = 1.0_ireals
        ENDIF
        IF (lolp3_gme(3) .EQV. lolp_lm(l1,l2)) THEN
          ipx  = ipx  + 1
          weights_ip (l1,l2,i_2) = 1.0_ireals
        ENDIF

        IF (idebug > 50) THEN
          IF (ipx /= 0) THEN
            WRITE (*,'(A,2I5,6I6)')  'three:   ', l1, l2, j1, j2, jd, i_c, i_1, i_2
          ENDIF
        ENDIF

        IF (ipx == 0) THEN
          ! If no point matches up to now, look at the 6 surrounding points
          ! and take the ones with identical characteristic

          ! Point 5
          IF (lolp_gme(j1  ,j2-1,jd) .EQV. lolp_lm(l1,l2)) THEN
            ipx=ipx + 1
            weights_ip (l1,l2, 5) = 1.0_ireals
          ENDIF

          ! Point 6
          IF (lolp_gme(j1+1,j2-1,jd) .EQV. lolp_lm(l1,l2)) THEN
            ipx=ipx + 1
            weights_ip (l1,l2, 6) = 1.0_ireals
          ENDIF

          ! Point 9
          IF (lolp_gme(j1-1,j2  ,jd) .EQV. lolp_lm(l1,l2)) THEN
            ipx=ipx + 1
            weights_ip (l1,l2, 9) = 1.0_ireals
          ENDIF

          ! Point 11
          IF (lolp_gme(j1+1,j2  ,jd) .EQV. lolp_lm(l1,l2)) THEN
            ipx=ipx + 1
            weights_ip (l1,l2,11) = 1.0_ireals
          ENDIF

          ! Point 14
          IF (lolp_gme(j1-1,j2+1,jd) .EQV. lolp_lm(l1,l2)) THEN
            ipx=ipx + 1
            weights_ip (l1,l2,14) = 1.0_ireals
          ENDIF

          ! Point 15
          IF (lolp_gme(j1  ,j2+1,jd) .EQV. lolp_lm(l1,l2)) THEN
            ipx=ipx + 1
            weights_ip (l1,l2,15) = 1.0_ireals
          ENDIF

          IF (idebug > 50) THEN
            IF (ipx /= 0) THEN
              WRITE (*,'(A,2I5,3I6)')  'six:     ', l1, l2, j1, j2, jd
            ENDIF
          ENDIF

          IF (ipx == 0) THEN
            ! If no point matches up to now, look at the next 12 surrounding points
            ! and take the ones with identical characteristic

            ! Point 1
            IF (lolp_gme(j1  ,j2-2,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2, 1) = 1.0_ireals
            ENDIF

            ! Point 2
            IF (lolp_gme(j1+1,j2-2,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2, 2) = 1.0_ireals
            ENDIF

            ! Point 3
            IF (lolp_gme(j1+2,j2-2,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2, 3) = 1.0_ireals
            ENDIF

            ! Point 4
            IF (lolp_gme(j1-1,j2-1,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2, 4) = 1.0_ireals
            ENDIF

            ! Point 7
            IF (lolp_gme(j1+2,j2-1,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2, 7) = 1.0_ireals
            ENDIF

            ! Point 8
            IF (lolp_gme(j1-2,j2  ,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2, 8) = 1.0_ireals
            ENDIF

            ! Point 12
            IF (lolp_gme(j1+2,j2  ,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2,12) = 1.0_ireals
            ENDIF

            ! Point 13
            IF (lolp_gme(j1-2,j2+1,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2,13) = 1.0_ireals
            ENDIF

            ! Point 16
            IF (lolp_gme(j1+1,j2+1,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2,16) = 1.0_ireals
            ENDIF

            ! Point 17
            IF (lolp_gme(j1-2,j2+2,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2,17) = 1.0_ireals
            ENDIF

            ! Point 18
            IF (lolp_gme(j1-1,j2+2,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2,18) = 1.0_ireals
            ENDIF

            ! Point 19
            IF (lolp_gme(j1  ,j2+2,jd) .EQV. lolp_lm(l1,l2)) THEN
              ipx=ipx + 1
              weights_ip (l1,l2,19) = 1.0_ireals
            ENDIF

            IF (idebug > 50) THEN
              IF (ipx /= 0) THEN
                WRITE (*,'(A,2I5,3I6)')  'twelve:  ', l1, l2, j1, j2, jd
               ENDIF
            ENDIF

      !------------------------------------------------------------------------
      ! Section 4: Now we have to look for a far away GME point
      !------------------------------------------------------------------------

            IF (ipx == 0) THEN
              ! ipx = 0 for all surrounding points in the neighborhood
              ! As a last resort, find the next GME-Gridpoint in the
              ! same diamond which has the same LSM as in LM/HM
              ! During this search take care of undefined GME grid points 

              ! But then you have to take the value at the special grid point,
              ! so this cannot be done with weights!!!!
              ! The indices of these points have to be stored

              mindist = 10000000

              DO k2 = kig2sm2,kig2ep2
                DO k1 = kig1sm2,kig1ep2
                  IF ( (lolp_gme(k1,k2,jd) .EQV. lolp_lm(l1,l2)) .AND. &
                       (field_gme(k1,k2,jd) /= undef) ) THEN
                    ndist = ABS(k1-j1) + ABS(k2-j2)
                    IF (ndist < mindist) THEN
                      mindist = ndist
                      log_ip   (l1,l2)    = .TRUE.
                      match_ip (l1,l2, 1) = k1
                      match_ip (l1,l2, 2) = k2
                      match_ip (l1,l2, 3) = jd
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO

              IF (mindist == 10000000) THEN
                kierr = 2
                yerrmsg = '***  No point with same lsm found in neighborhood!!  '
                WRITE (yerrmsg(55:79),'(5I5)') l1, l2, j1, j2, jd
              ELSE
                IF (idebug > 50) THEN
                  WRITE (*,'(A,2I5,3I6)')  'far away:', l1, l2, k1, k2, jd
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF

        IF (ipx /= 0) THEN
          ! build the mean for all points:
!CDIR NOVECTOR
          DO i = 1, 19
            weights_ip (l1,l2,i) = weights_ip (l1,l2,i) / ipx
          ENDDO
        ENDIF
      ENDIF

    ENDDO
  ENDDO

!------------------------------------------------------------------------------

END SUBROUTINE init_gme_interpol

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE pp_interp2ls (px, kig1sm2, kig1ep2, kig2sm2, kig2ep2, kids, kide, &
                         kispoke, pbaryll, kindex, lmono, lposdef, yitype,   &
                         lolp_gme ,lolp_lm, weights_ip, nearest_ip,          &
                         match_ip, log_ip, undef,                            &
                         pxi    , kilons , kilone , kilats , kilate ,        &
                         yerrmsg, kierr  )

!------------------------------------------------------------------------------
!
! Description:
!   *pp_interp2ls* interpolates the scalar field ("px") linearly to the
!   points with the barycentric coordinates "pbaryll" using the 
!   values at the three neighbours.
!   The array "index" gives the coordinates (j1, j2, jd) of the GME
!   grid point being the nearest to the points of the regular grid and
!   the triangle containing these points
!
! Method:
!
!==============================================================================
!
! Input
INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kig1sm2,  & ! first  dimension of array "px",   start index
  kig1ep2,  & ! first  dimension of array "px",   end   index
  kig2sm2,  & ! second dimension of array "px",   start index
  kig2ep2,  & ! second dimension of array "px",   end   index
  kids   ,  & ! start diamond
  kide        ! end diamond

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kilons ,  & ! first  dimension of array "pxi",  start index
  kilone ,  & ! first  dimension of array "pxi",  end   index
  kilats ,  & ! second dimension of array "pxi",  start index
  kilate      ! second dimension of array "pxi",  end   index

REAL    (KIND=ireals)   , INTENT(IN)   ::    &
  px   (kig1sm2:kig1ep2, kig2sm2:kig2ep2,kids:kide),      &
              ! scalar field to be interpolated
  pbaryll  (kilons  :kilone  , kilats  :kilate  , 2),     &
              ! first and second barycentric coordinate of points
  weights_ip (kilons:kilone, kilats:kilate,19),           &
            !  Weights for the interpolation, to choose, which 
            !  points have to be taken for the match interpolation
  undef       ! value assigned to undefined grid points

INTEGER (KIND=iintegers), INTENT(IN)    ::      &
  nearest_ip (kilons:kilone, kilats:kilate, 3), &
            !  indices of nearest GME gridpoint for nearest neighbor interpolation
  match_ip   (kilons:kilone, kilats:kilate, 3)
            !  indices of nearest GME gridpoint for match interpolation

LOGICAL                 , INTENT(IN)    ::      &
  log_ip     (kilons:kilone, kilats:kilate)
            !  to use a far away GME gridpoint with same lsm

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kindex  (kilons  :kilone  , kilats  :kilate  , 4)
              ! index of the triangle containing the point: j1,j2,jd,mt

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kispoke (12)    ! Normal spokes 

LOGICAL, INTENT(IN)                    ::    &
  lmono, lposdef  ! Switch for monotonicity and positive definiteness

LOGICAL, INTENT(IN)                    ::    &
  lolp_gme (kig1sm2:kig1ep2, kig2sm2:kig2ep2,kids:kide),   & !
  lolp_lm  (kilons:kilone, kilats:kilate)
                  ! Land Sea Mask needed for 'M'atch Interp.

CHARACTER*1, INTENT(IN) :: yitype  ! Interpolation type (L, N, M)
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
  za , zb , zc, & ! Barycentric coordinates
  zx1, zx2, zx3   ! values at three gridpoints

INTEGER (KIND=iintegers)   ::      &
  l1, l2,             & ! Loop indices
  j1, j2, jd, m1, m2    ! Gridpoint indices

!=======================================================================

  kierr   = 0
  yerrmsg = '    '

!------------------------------------------------------------------------------
! Section 1: Nearest Neighbor interpolation
!------------------------------------------------------------------------------

  IF     (yitype == 'N') THEN

    DO l2 = kilats, kilate
      DO l1 = kilons, kilone
        j1         = nearest_ip(l1,l2, 1)
        j2         = nearest_ip(l1,l2, 2)
        jd         = nearest_ip(l1,l2, 3)
        pxi(l1,l2) = px(j1,j2,jd)
      ENDDO
    ENDDO

!------------------------------------------------------------------------------
! Section 2: Match interpolation
!------------------------------------------------------------------------------

  ELSEIF (yitype == 'M') THEN

    DO l2 = kilats, kilate
      DO l1 = kilons, kilone

        IF (.NOT. log_ip(l1,l2)) THEN
          j1 = kindex(l1,l2,1) ! first gme grid dimension
          j2 = kindex(l1,l2,2) ! second gme grid dimension
          jd = kindex(l1,l2,3) ! number of diamond

          pxi(l1,l2) =  px(j1  ,j2-2,jd) * weights_ip (l1,l2, 1) +     &
                        px(j1+1,j2-2,jd) * weights_ip (l1,l2, 2) +     &
                        px(j1+2,j2-2,jd) * weights_ip (l1,l2, 3) +     &
                        px(j1-1,j2-1,jd) * weights_ip (l1,l2, 4) +     &
                        px(j1  ,j2-1,jd) * weights_ip (l1,l2, 5) +     &
                        px(j1+1,j2-1,jd) * weights_ip (l1,l2, 6) +     &
                        px(j1+2,j2-1,jd) * weights_ip (l1,l2, 7) +     &
                        px(j1-2,j2  ,jd) * weights_ip (l1,l2, 8) +     &
                        px(j1-1,j2  ,jd) * weights_ip (l1,l2, 9) +     &
                        px(j1  ,j2  ,jd) * weights_ip (l1,l2,10) +     &
                        px(j1+1,j2  ,jd) * weights_ip (l1,l2,11) +     &
                        px(j1+2,j2  ,jd) * weights_ip (l1,l2,12) +     &
                        px(j1-2,j2+1,jd) * weights_ip (l1,l2,13) +     &
                        px(j1-1,j2+1,jd) * weights_ip (l1,l2,14) +     &
                        px(j1  ,j2+1,jd) * weights_ip (l1,l2,15) +     &
                        px(j1+1,j2+1,jd) * weights_ip (l1,l2,16) +     &
                        px(j1-2,j2+2,jd) * weights_ip (l1,l2,17) +     &
                        px(j1-1,j2+2,jd) * weights_ip (l1,l2,18) +     &
                        px(j1  ,j2+2,jd) * weights_ip (l1,l2,19)

        ELSE
          ! a far-away GME grid point has to be used
          j1 = match_ip(l1,l2,1) ! first gme grid dimension
          j2 = match_ip(l1,l2,2) ! second gme grid dimension
          jd = match_ip(l1,l2,3) ! number of diamond

          pxi (l1,l2) = px (j1,j2,jd)
        ENDIF

      ENDDO
    ENDDO

!------------------------------------------------------------------------------
! Section 3: Linear interpolation
!------------------------------------------------------------------------------

  ELSEIF (yitype == 'L') THEN

    DO l2 = kilats, kilate
      DO l1 = kilons, kilone

      j1 = kindex(l1,l2,1) ! first gme grid dimension
      j2 = kindex(l1,l2,2) ! second gme grid dimension
      jd = kindex(l1,l2,3) ! number of diamond
      m1 = kindex(l1,l2,4) ! Number of triangle for interpolation
      m2 = MOD (m1,6) + 1

      ! The 3 barycentric coordinates of the point for which the inter-
      ! polation is required
      za    = pbaryll(l1,l2,1)
      zb    = pbaryll(l1,l2,2)
      zc    = 1. - za - zb

      ! Interpolate px to the point with the barycentric coordinates
      ! (za, zb, zc)
      zx1 = px(j1,j2,jd)
      zx2 = px(j1+kispoke(m1), j2+kispoke(m1+6), jd )
      zx3 = px(j1+kispoke(m2), j2+kispoke(m2+6), jd )

      pxi(l1,l2) = zc*zx1 + za*zx2 + zb*zx3

      ! Enforce monotonicity, if required
      IF (lmono) THEN
        pxi(l1,l2) = MIN (pxi(l1,l2), MAX (zx1,zx2,zx3))
        pxi(l1,l2) = MAX (pxi(l1,l2), MIN (zx1,zx2,zx3))
      ENDIF

      ! Enforce positive definiteness, if required
      IF (lposdef) THEN
        pxi(l1,l2) = MAX (pxi(l1,l2), 0.0_ireals)
      ENDIF

      ENDDO
    ENDDO

  ENDIF

!------------------------------------------------------------------------------

END SUBROUTINE pp_interp2ls

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE pp_interp2ls_save (px, kig1sm2, kig1ep2, kig2sm2, kig2ep2, kids, kide, &
                         kispoke, pbaryll, kindex, lmono, lposdef, yitype,   &
                         lolp_gme ,lolp_lm, undef,                           &
                         pxi    , kilons , kilone , kilats , kilate ,        &
                         yerrmsg, kierr  )

!------------------------------------------------------------------------------
!
! Description:
!   *pp_interp2ls* interpolates the scalar field ("px") linearly to the
!   points with the barycentric coordinates "pbaryll" using the 
!   values at the three neighbours.
!   The array "index" gives the coordinates (j1, j2, jd) of the GME
!   grid point being the nearest to the points of the regular grid and
!   the triangle containing these points
!
! Method:
!
!==============================================================================
!
! Input
INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kig1sm2,  & ! first  dimension of array "px",   start index
  kig1ep2,  & ! first  dimension of array "px",   end   index
  kig2sm2,  & ! second dimension of array "px",   start index
  kig2ep2,  & ! second dimension of array "px",   end   index
  kids   ,  & ! start diamond
  kide        ! end diamond

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kilons ,  & ! first  dimension of array "pxi",  start index
  kilone ,  & ! first  dimension of array "pxi",  end   index
  kilats ,  & ! second dimension of array "pxi",  start index
  kilate      ! second dimension of array "pxi",  end   index

REAL    (KIND=ireals)   , INTENT(IN)   ::    &
  px   (kig1sm2:kig1ep2, kig2sm2:kig2ep2,kids:kide),      &
              ! scalar field to be interpolated
  pbaryll  (kilons  :kilone  , kilats  :kilate  , 2),     &
              ! first and second barycentric coordinate of points
  undef       ! value assigned to undefined grid points

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kindex  (kilons  :kilone  , kilats  :kilate  , 4)
              ! index of the triangle containing the point: j1,j2,jd,mt

INTEGER (KIND=iintegers), INTENT(IN)   ::    &
  kispoke (12)    ! Normal spokes 

LOGICAL, INTENT(IN)                    ::    &
  lmono, lposdef  ! Switch for monotonicity and positive definiteness

LOGICAL, INTENT(IN)                    ::    &
  lolp_gme (kig1sm2:kig1ep2, kig2sm2:kig2ep2,kids:kide),   & !
  lolp_lm  (kilons:kilone, kilats:kilate)
                  ! Land Sea Mask needed for 'M'atch Interp.

CHARACTER*1, INTENT(IN) :: yitype  ! Interpolation type (L, N, M)
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
  zsn,          & ! Hemisphere discriminator (+1 or -1)
  za , zb , zc, & ! Barycentric coordinates
  zx1, zx2, zx3,& ! values at three gridpoints
  zspx            ! sum of surrounding values

INTEGER (KIND=iintegers)   ::      &
  l1, l2,             & ! Loop indices
  j1, j2, jd, m1, m2, & ! Gridpoint indices
  mindist, ndist,     & !
  i,k1,k2               ! Loop indices

INTEGER (KIND = iintegers) ::      &
  ipx ! number of surrounding values

LOGICAL                    ::      &
  lolp3_gme(3),    & ! LSM of the three surrounding grid points
  lngp,            & ! take nearest gridpoint
  lmatch_ngp         ! take mean over nearest grid points in case of match 
                     ! interpolation

!=======================================================================

  kierr   = 0
  yerrmsg = '    '

!=======================================================================

  lngp       = .false.
  lmatch_ngp = .false.
  IF ( yitype == 'N') lngp = .true.

  DO l2 = kilats, kilate
    DO l1 = kilons, kilone

!------------------------------------------------------------------------------

      j1    = kindex(l1,l2,1) ! first gme grid dimension
      IF (j1 == -9999) CYCLE
      j2    = kindex(l1,l2,2) ! second gme grid dimension
      jd    = kindex(l1,l2,3) ! number of diamond

      !     For the diamonds 1 to  5, set "zsn" to  1.
      !     For the diamonds 6 to 10, set "zsn" to -1.
      zsn = 1.
      IF (jd .GE. 6) zsn = -1.

      m1    = kindex(l1,l2,4) ! Number of triangle for interpolation
      m2    = MOD (m1,6) + 1

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
        lolp3_gme(1) = lolp_gme(j1,j2,jd)
        lolp3_gme(2) = lolp_gme(j1+kispoke(m1),j2+kispoke(m1+6),jd)
        lolp3_gme(3) = lolp_gme(j1+kispoke(m2),j2+kispoke(m2+6),jd)

        IF ( ALL (lolp3_gme .EQV. lolp_lm(l1,l2)) ) THEN
          lmatch_ngp = .false.
        ELSE
          lmatch_ngp = .true.
        ENDIF
      ENDIF

      ! Also set lolp3 for yitype = 'N'
      IF (yitype == 'N') THEN
        lolp3_gme(1) = lolp_gme(j1,j2,jd)
        lolp3_gme(2) = lolp_gme(j1+kispoke(m1),j2+kispoke(m1+6),jd)
        lolp3_gme(3) = lolp_gme(j1+kispoke(m2),j2+kispoke(m2+6),jd)
      ENDIF

!------------------------------------------------------------------------------

      IF (lngp) THEN
        ! Checking for undefined values is not necessary here, because the
        ! search is only local around the grid point considered and these
        ! points have to be provided by the bitmap.
        IF (lolp_gme(j1,j2,jd) .EQV. lolp_lm(l1,l2)) THEN
          pxi(l1,l2) = px(j1,j2,jd)
        ELSEIF (lolp_gme(j1+kispoke(m1),j2+kispoke(m1+6),jd)    &
                                               .EQV. lolp_lm(l1,l2)) THEN
          pxi(l1,l2) = px(j1+kispoke(m1), j2+kispoke(m1+6), jd )
        ELSEIF (lolp_gme(j1+kispoke(m2),j2+kispoke(m2+6),jd)    &
                                               .EQV. lolp_lm(l1,l2)) THEN
          pxi(l1,l2) = px(j1+kispoke(m2), j2+kispoke(m2+6), jd )
        ELSE
          pxi(l1,l2) = 1E20_ireals
          kierr = 2
          WRITE (yerrmsg, '(A,I5,A,I5)')                                  &
             'yitype = N: no GME point for lat i = ',l2,' and lon j = ',l1
        ENDIF

        CYCLE   ! the innermost DO l1-loop
      ENDIF

!------------------------------------------------------------------------------

      IF (lmatch_ngp) THEN

        ! 1. If one or two of the three surrounding points are of equal
        !    LSM as in LM/HM (i.e. ipx /= 0) --> take mean of these points, 
        !    else, see below par 2. 
        zspx = 0
        ipx  = 0
        DO i=1,3
          IF (lolp3_gme(i).EQV.lolp_lm(l1,l2)) THEN
            SELECT CASE (i)
            CASE (1)
              zspx=zspx+px(j1,j2,jd)
            CASE (2)
              zspx=zspx+px(j1+kispoke(m1), j2+kispoke(m1+6), jd )
            CASE (3)   
              zspx=zspx+px(j1+kispoke(m2), j2+kispoke(m2+6), jd )
            END SELECT
            ipx=ipx+1
          ENDIF ! lolp3_gme
        END DO
            
        IF (ipx /= 0) THEN
          pxi(l1,l2)=zspx / REAL(ipx,ireals)
        ELSE
          ! 2. If one or more of the six surrounding points are of equal
          !    LSM as in LM/HM (i.e. ipx /= 0) --> take mean of these points, 
          !    else, see below par 3. 
          zspx = 0
          ipx  = 0
          DO i=1,6
            IF (lolp_gme(j1+kispoke(i),j2+kispoke(i+6),jd)    &
                                                 .EQV. lolp_lm(l1,l2)) THEN
               zspx=zspx+px(j1+kispoke(i),j2+kispoke(i+6),jd)
               ipx=ipx+1
            ENDIF
          ENDDO

          IF (ipx /= 0) THEN
            pxi(l1,l2)=zspx / REAL(ipx,ireals)
          ELSE ! ipx = 0 in the six surrounding points

            ! 3. If none of the six surrounding points are of equal
            !    LSM as in LM/HM, then look at all 12 GME-points 
            !    around the six surrounding points (ispoke)
            !    else set the value to 1.E20 
            zspx = 0
            ipx  = 0
            DO k2=-2, +2
              DO k1= -2, +2
                IF (k1+k2 > 2 .OR. k1+k2 < -2) CYCLE
                IF (lolp_gme(j1+k1,j2+k2,jd) .EQV. lolp_lm(l1,l2)) THEN
                  zspx=zspx+px(j1+k1,j2+k2,jd)
                  ipx=ipx+1
                ENDIF
              ENDDO
            ENDDO
            IF (ipx /= 0) THEN
              pxi(l1,l2)=zspx / REAL(ipx,ireals)
            ELSE ! ipx = 0 for all 12 points
            !
            ! 4. As a last resort, find the next GME-Gridpoint in the
            !    same diamond which has the same LSM as in LM/HM
            !    During this search take care of undefined GME grid points 

                  mindist = 10000000

                  DO k2 = kig2sm2,kig2ep2
                    DO k1 = kig1sm2,kig1ep2
                      IF ( (lolp_gme(k1,k2,jd) .EQV. lolp_lm(l1,l2)) .AND. &
                           (px(k1,k2,jd) /= undef) ) THEN
                        ndist = ABS(k1-j1) + ABS(k2-j2)
                        IF (ndist < mindist) THEN
                          mindist = ndist
                          zspx    = px(k1,k2,jd)
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDDO

                  ! If still no point is found, give up finally

                  IF (  mindist == 10000000) THEN
                    kierr = 1
                    WRITE (yerrmsg, '(A,T50,L6,I5,I5)')                   &
                        'ERROR in pp_interp2ls: unique lolp_lm(l1,l2)',   &
                        lolp_lm(l1,l2),l1,l2
                    RETURN
                  ENDIF

                  pxi(l1,l2) = zspx

            ENDIF  ! ipx = 0 in the twelve surrounding points
          ENDIF ! ipx = 0 in the six surrounding points
        ENDIF   ! ipx = 0 in the triangle

        CYCLE   ! the innermost DO-loop

      ENDIF      ! lmatch_ngp

!------------------------------------------------------------------------------

      ! The 3 barycentric coordinates of the point for which the inter-
      ! polation is required
      za    = pbaryll(l1,l2,1)
      zb    = pbaryll(l1,l2,2)
      zc    = 1. - za - zb

      ! Interpolate px to the point with the barycentric coordinates
      ! (za, zb, zc)
      zx1 = px(j1,j2,jd)
      zx2 = px(j1+kispoke(m1), j2+kispoke(m1+6), jd )
      zx3 = px(j1+kispoke(m2), j2+kispoke(m2+6), jd )

      pxi(l1,l2) = zc*zx1 + za*zx2 + zb*zx3

      ! Enforce monotonicity, if required
      IF (lmono) THEN
        pxi(l1,l2) = MIN (pxi(l1,l2), MAX (zx1,zx2,zx3))
        pxi(l1,l2) = MAX (pxi(l1,l2), MIN (zx1,zx2,zx3))
      ENDIF

      ! Enforce positive definiteness, if required
      IF (lposdef) THEN
        pxi(l1,l2) = MAX (pxi(l1,l2), 0.0_ireals)
      ENDIF

    ENDDO
  ENDDO

!------------------------------------------------------------------------------

END SUBROUTINE pp_interp2ls_save

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE pp_interp2qs (px     , kig1sm2, kig1ep2 , kig2sm2, kig2ep2,     &
                         peta   , pchi    , pcpsi  , pspsi  ,              &
                         pgrd   , kids    , kide   ,                       &
                         kig1sm1, kig1ep1 , kig2sm1, kig2ep1,              &
                         kispoke, kispokes, ki1mrp , ki2mrp ,              &
                         pbaryll, kindex  ,                                &
                         lmono  , lposdef ,                                &
                         pxi    , gpxeta  , gpxchi ,                       &
                         kilons , kilone  , kilats , kilate ,              &
                         j1min  , j1max   , j2min  , j2max  , jdmin, jdmax,&
                         isp11  , isp12   , isp13  , isp14, isp15, isp16,  &
                         isp21  , isp22   , isp23  , isp24, isp25, isp26,  &
                         Re     , loc_debug , kierr)

!------------------------------------------------------------------------------
!
! Description:
!   *pp_interp2qs* interpolates the scalar field "px" quadratically 
!   to the point with the barycentric coordinates "pbaryll",
!   using the values at 12 neighbours. 
!
! Method:
!
!==============================================================================
!
! Declarations:
!
! Input
  ! dimensions of icosahedral grid
  INTEGER (KIND=iintegers)   ::     &
    kig1sm2, & ! first  dimension of array "px" , start index
    kig1ep2, & ! first  dimension of array "px" , end   index
    kig2sm2, & ! second dimension of array "px" , start index
    kig2ep2, & ! second dimension of array "px" , end   index
    kig1sm1, & ! first  dimension of array "pgrd", "peta", "pchi", 
               ! "pcpsi", "pspsi",   start index
    kig1ep1, & ! first  dimension of array "pgrd", "peta", "pchi",
               ! "pcpsi", "pspsi",   end   index
    kig2sm1, & ! second dimension of array "pgrd", start index
    kig2ep1, & ! second dimension of array "pgrd", end   index
    kids   , & ! start number of diamonds
    kide       ! end number of diamonds

  ! dimensions of rectangular grid
  INTEGER (KIND=iintegers)   ::     &
    kilons , & ! first  dimension of array "pxi",  start index
    kilone , & ! first  dimension of array "pxi",  end   index
    kilats , & ! second dimension of array "pxi",  start index
    kilate     ! second dimension of array "pxi",  end   index

  ! min- and max-indices of j1-, j2- and jd-direction
  INTEGER (KIND=iintegers)   ::     &
    j1min, j1max, j2min, j2max, jdmin, jdmax

  ! fields from GME on icosahedral grid
  REAL    (KIND=ireals)      ::     &
    px   (kig1sm2:kig1ep2, kig2sm2:kig2ep2, kids:kide),               &
           ! Scalar field to be interpolated
    peta (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7),                       &
           ! eta-coordinates of the 6 (5) surrounding gridpoints  
    pchi (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7),                       &
           ! chi-coordinates of the 6 (5) surrounding gridpoints  
    pcpsi(kig1sm1:kig1ep1, kig2sm1:kig2ep1, 6),                       &
           ! cosine of the rotation angle for the wind rotation
    pspsi(kig1sm1:kig1ep1, kig2sm1:kig2ep1, 6),                       &
           ! sine   of the rotation angle for the wind rotation
    pgrd (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7, 2)
           ! coefficients for the calculation of the gradients

  ! fields on regular grid
  REAL    (KIND=ireals)      ::     &
    pbaryll (kilons  :kilone  , kilats  :kilate  , 2)
           ! first  and second barycentric coordinate of points

  ! Fields for local gradients into the eta- and chi-direction
  REAL    (KIND=ireals)      ::     &
    gpxeta (kig1sm1 :kig1ep1,  kig2sm1 :kig2ep1,  jdmin:jdmax),  &
    gpxchi (kig1sm1 :kig1ep1,  kig2sm1 :kig2ep1,  jdmin:jdmax)

  INTEGER (KIND=iintegers)   ::     &
    kispoke (12),   & ! Normal spokes 
    kispokes(12,4), & ! Spokes of the 4 mirrored points of the
                      ! extended array
    ki1mrp (8),     & ! i1-indices of the 4 mirrored points
    ki2mrp (8),     & ! i2-indices of the 4 mirrored points
    kindex (kilons  :kilone  , kilats  :kilate  , 4)

  INTEGER (KIND=iintegers)   ::     &
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
    isp26          !

  REAL (KIND=ireals)         :: Re  ! mean radius of the earth

  LOGICAL                    ::     &
    loc_debug,   & ! debug switch; if .true. print information
    lmono ,      & ! monotonic advection, if set to .true.
    lposdef        ! positive definite advection, if set to .true.

!==============================================================================
!
! Output
  INTEGER (KIND=iintegers)   ::     &
    kierr     ! error flag, kierr = 0 if no error occured

  REAL    (KIND=ireals)      ::   pxi (*)
           ! Interpolated scalar field
           ! If running on one processor and *kindex* is fully set,
           ! pxi may be defined as pxi(kilons:kilone, kilats:kilate)
           ! in the calling program
 
!==============================================================================
!
! Local arrays
!
  ! Local variables
  REAL    (KIND=ireals)      ::     &
    zsn,                            & ! Hemisphere discriminator (+1 or -1)
    zgeta1, zgeta2, zgchi1, zgchi2, & ! Rotated gradients
    zx1, zx2, zx3, zx4, zx5, zx6,   & ! Local values on triangle
    za , zb , zc,                   & ! Barycentric coordinates
    zred8,                          & ! Re/8
    zdp12 , zdp13 , zdp23

  INTEGER (KIND=iintegers)   ::     &
       mig1s, mig2s, mig1e, mig2e,  & ! 'core'-dimensions
       l1, l2, ld,                  & ! Loop indices
       j1, j2, jd, m1, m2, mi1p1,   & !
       mi2p1, mi1p2, mi2p2,         & ! Gridpoint indices
       mp                   ! counter for pxi
!
!==============================================================================
!
! Check input variables
!
! The index array "kindex" contains the dimensions j1, j2 and jd
! of the points within array "px" which are the nodes for the  
! interpolation calculation; for quadratic interpolation there are
! needed two surrounding rows.
! Check, if these are available:
!
  mig1s = kig1sm2 + 2
  mig2s = kig2sm2 + 2
  mig1e = kig1ep2 - 2
  mig2e = kig2ep2 - 2

  IF ( (j1min .LT. mig1s) .OR. (j1max .GT. mig1e) .OR.     &
       (j2min .LT. mig2s) .OR. (j2max .GT. mig2e) .OR.     &
       (jdmin .LT. kids ) .OR. (jdmax .GT. kide )  )  THEN
    PRINT *,'  Error in subroutine *pp_interp2qs*, dimensions of ',   &
         'array "px" are not covering the computation domain'
    PRINT *,'  mig1s: ', mig1s,'  j1min: ', j1min
    PRINT *,'  mig1e: ', mig1e,'  j1max: ', j1max
    PRINT *,'  mig2s: ', mig2s,'  j2min: ', j2min
    PRINT *,'  mig2e: ', mig2e,'  j2max: ', j2max
    PRINT *,'  kids : ', kids ,'  jdmin: ', jdmin
    PRINT *,'  kide : ', kide ,'  jdmax: ', jdmax
    kierr = -1
    RETURN
  ENDIF
!
!
!==============================================================================

  zred8   = 0.125*Re

!==============================================================================

  DO ld = jdmin, jdmax     ! Loop over the diamonds      
    ! Calculate the gradients of "px" into the eta- and chi-direction
    CALL grad2s (px(kig1sm2,kig2sm2,ld),                                    &
                 kig1sm2 , kig1ep2 , kig2sm2 , kig2ep2 , 1, 1, pgrd,        &
                 kig1sm1 , kig1ep1 , kig2sm1 , kig2ep1 ,                    &
                 kispokes, ki1mrp  , ki2mrp  , ld      ,                    &
                 j1min-1 , j1max+1 , j2min-1 , j2max+1 , loc_debug,         &
                 isp11   , isp12   , isp13   , isp14   , isp15 , isp16,     &
                 isp21   , isp22   , isp23   , isp24   , isp25 , isp26,     &
                 gpxeta(mig1s-1,mig2s-1,ld) ,                               &
                 gpxchi(mig1s-1,mig2s-1,ld) ,                               &
                 mig1s-1 , mig1e+1 , mig2s-1 , mig2e+1 , kierr)

    IF (kierr .NE. 0) THEN
      PRINT *,'  Error in *pp_interp2qs* calling *grad2s*'
      RETURN
    ENDIF
  ENDDO

!==============================================================================
!
! Compute the auxiliary points at the midpoints of the triangle
! sides always taking care of the proper rotation of the gradients
! into the local system of the central node
!
  mp = 0

  DO l2 = kilats, kilate
    DO l1 = kilons, kilone

      j1    = kindex(l1,l2,1)
      IF (j1 == -9999) CYCLE
      j2    = kindex(l1,l2,2)
      jd    = kindex(l1,l2,3)

      mp = mp+1

      ! For the diamonds 1 to  5, set "zsn" to  1.
      ! For the diamonds 6 to 10, set "zsn" to -1.
      zsn = 1.

      IF (jd .GE. 6) zsn = -1.

      m1    = kindex(l1,l2,4) ! Number of triangle for interpolation
      m2    = MOD (m1,6) + 1
      mi1p1 = j1 + kispoke(m1)   ! i1-index of point 1
      mi2p1 = j2 + kispoke(m1+6) ! i2-index of point 1
      mi1p2 = j1 + kispoke(m2)   ! i1-index of point 2
      mi2p2 = j2 + kispoke(m2+6) ! i2-index of point 2

      ! The 3 barycentric coordinates of the point for which the
      ! interpolation is required
      za    = pbaryll(l1,l2,1)
      zb    = pbaryll(l1,l2,2)
      zc    = 1. - za - zb

!==============================================================================

      ! Rotate the gradients at the two points (P1 and P2) into the local
      ! spherical system of point P0 (central node)

      ! Point m1, gradient in eta-direction, first comp. (i1-direction)
      zgeta1 = pcpsi(j1,j2,m1)*gpxeta(mi1p1,mi2p1,jd) + zsn*             &
               pspsi(j1,j2,m1)*gpxchi(mi1p1,mi2p1,jd)

      ! Point m1, gradient in eta-direction, secn. comp. (i2-direction)
      zgchi1 = pcpsi(j1,j2,m1)*gpxchi(mi1p1,mi2p1,jd) - zsn*             &
               pspsi(j1,j2,m1)*gpxeta(mi1p1,mi2p1,jd)

      ! Point m2, gradient in chi-direction, first comp. (i1-direction)
      zgeta2 = pcpsi(j1,j2,m2)*gpxeta(mi1p2,mi2p2,jd) + zsn*             &
               pspsi(j1,j2,m2)*gpxchi(mi1p2,mi2p2,jd)

      ! Point m2, gradient in chi-direction, secn. comp. (i2-direction)
      zgchi2 = pcpsi(j1,j2,m2)*gpxchi(mi1p2,mi2p2,jd) - zsn*             &
               pspsi(j1,j2,m2)*gpxeta(mi1p2,mi2p2,jd)

!==============================================================================

      ! Interpolate the increments at the midpoints of the triangle sides
      ! using the gradients at the three triangle vertices
      zdp12  = (gpxeta(j1,j2,jd) - zgeta1)*peta(j1,j2,m1) + zsn*         &
               (gpxchi(j1,j2,jd) - zgchi1)*pchi(j1,j2,m1)

      zdp13  = (gpxeta(j1,j2,jd) - zgeta2)*peta(j1,j2,m2) + zsn*         &
               (gpxchi(j1,j2,jd) - zgchi2)*pchi(j1,j2,m2)

      zdp23  = (zgeta1 - zgeta2)* (peta(j1,j2,m2) - peta(j1,j2,m1)) + zsn*  &
               (zgchi1 - zgchi2)* (pchi(j1,j2,m2) - pchi(j1,j2,m1))

      ! Determine the values at the vertices of the triangle (# 1, 2, 3)
      ! and at the midpoints of the triangle sides (# 4, 5, 6)
      zx1    = px(j1,j2,jd)
      zx2    = px(mi1p1,mi2p1,jd) 
      zx3    = px(mi1p2,mi2p2,jd)
      zx4    = 0.5*(zx1 + zx2) + zred8*zdp12
      zx5    = 0.5*(zx2 + zx3) + zred8*zdp23
      zx6    = 0.5*(zx1 + zx3) + zred8*zdp13

      ! Interpolated value using the barycentric coordinates
      pxi(mp) = zx1*zc*(2.*zc - 1.) + zx2*za*(2.*za - 1.) +                 &
                zx3*zb*(2.*zb - 1.) + 4.*((zx4*za*zc + zx5*za*zb) + zx6*zb*zc)

      ! Enforce monotonicity, if required
      IF (lmono) THEN
        pxi(mp) = MIN (pxi(mp), MAX (zx1,zx2,zx3))
        pxi(mp) = MAX (pxi(mp), MIN (zx1,zx2,zx3))
      ENDIF

      ! Enforce positive definiteness, if required
      IF (lposdef) THEN
        pxi(mp) = MAX (pxi(mp), 0.0_ireals)
      ENDIF

    ENDDO
  ENDDO
!
!==============================================================================

  kierr = 0

!==============================================================================

  ! Print debug information, if required
  IF (loc_debug) THEN
    PRINT *,'  *pp_interp2qs* completed '
  ENDIF

!==============================================================================

END SUBROUTINE pp_interp2qs

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE pp_interp2qv (pu     , pv      ,                                &
                         kig1sm2, kig1ep2 , kig2sm2, kig2ep2,              &
                         kids   , kide    ,                                &
                         peta   , pchi    , pcpsi  , pspsi  ,              &
                         pgrd   ,                                          &
                         kig1sm1, kig1ep1 , kig2sm1, kig2ep1,              &
                         kispoke, kispokes, ki1mrp , ki2mrp ,              &
                         pbaryll, protang , kindex ,                       &
                         pui    , pvi     ,                                &
                         gudeta , gudchi  , gvdeta , gvdchi ,              &
                         kilons , kilone  , kilats , kilate ,              &
                         j1min  , j1max   , j2min  , j2max  , jdmin, jdmax,&
                         isp11  , isp12   , isp13  , isp14, isp15, isp16,  &
                         isp21  , isp22   , isp23  , isp24, isp25, isp26,  &
                         Re     , loc_debug , kierr   )

!------------------------------------------------------------------------------
!
! Description:
!  *pp_interp2qv* interpolates the vector field ("pu", "pv") quadrati- 
!  cally to the point with the barycentric coordinates "pbaryll",
!  using the values at 12 neighbours.
!
!==============================================================================
!
! Declarations:
!
! Input
!
  INTEGER (KIND=iintegers)   ::      &
    kig1sm2, & ! first  dimension of array "px"  , start index
    kig1ep2, & ! first  dimension of array "px"  , end   index
    kig2sm2, & ! second dimension of array "px"  , start index
    kig2ep2    ! second dimension of array "px"  , end   index

  INTEGER (KIND=iintegers)   ::      &
    kig1sm1, & ! first  dimension of array "pgrd", start index
    kig1ep1, & ! first  dimension of array "pgrd", end   index
    kig2sm1, & ! second dimension of array "pgrd", start index
    kig2ep1, & ! second dimension of array "pgrd", end   index
    kids   , & ! start number of diamonds
    kide       ! end number of diamonds

  INTEGER (KIND=iintegers)   ::      &
    kilons , & ! first  dimension of array "p*i",  start index
    kilone , & ! first  dimension of array "p*i",  end   index
    kilats , & ! second dimension of array "p*i",  start index
    kilate     ! second dimension of array "p*i",  end   index

  ! min- and max-indices of j1-, j2- and jd-direction
  INTEGER (KIND=iintegers)   ::     &
    j1min, j1max, j2min, j2max, jdmin, jdmax

  REAL (KIND=ireals)         ::      &
    pu   (kig1sm2 :kig1ep2,  kig2sm2 :kig2ep2, kids:kide),                &
               ! u-component of vector field
    pv   (kig1sm2 :kig1ep2,  kig2sm2 :kig2ep2, kids:kide),                &
               ! v-component of vector field
    peta (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7),                           &
               ! eta-coordinates of the 6 (5) surrounding gridpoints  
    pchi (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7),                           &
               ! chi-coordinates of the 6 (5) surrounding gridpoints  
    pcpsi(kig1sm1:kig1ep1, kig2sm1:kig2ep1, 6),                           &
               ! cosine of the rotation angle for the wind rotation
    pspsi(kig1sm1:kig1ep1, kig2sm1:kig2ep1, 6),                           &
               ! sine   of the rotation angle for the wind rotation
    pgrd (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7, 2)
               ! coefficients for the calculation of the gradients

  REAL (KIND=ireals)         ::      &
    pbaryll (kilons  :kilone  , kilats  :kilate  , 2),                    &
               ! first and second barycentric coordinate of points
    protang (kilons  :kilone  , kilats  :kilate  , 2)
               ! rotation angels for wind

  INTEGER (KIND=iintegers)   ::      &
    kindex  (kilons  :kilone  , kilats  :kilate  , 4)
               ! index of the triangle containing the point:j1, j2, jd, mt

  ! Fields for local gradients into the eta- and chi-direction
  REAL (KIND=ireals)      ::         &
   gudeta (kig1sm1 :kig1ep1,  kig2sm1 :kig2ep1, kids:kide),     &
   gudchi (kig1sm1 :kig1ep1,  kig2sm1 :kig2ep1, kids:kide),     &
   gvdeta (kig1sm1 :kig1ep1,  kig2sm1 :kig2ep1, kids:kide),     &
   gvdchi (kig1sm1 :kig1ep1,  kig2sm1 :kig2ep1, kids:kide)

  INTEGER (KIND=iintegers)   ::      &
    kispoke (12),   & ! Normal spokes 
    kispokes(12,4), & ! Spokes of the 4 mirrored points of the extended array
    ki1mrp (8),     & ! i1-indices of the 4 mirrored points
    ki2mrp (8)        ! i2-indices of the 4 mirrored points

  INTEGER (KIND=iintegers)   ::     &
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
    isp26          !

  REAL (KIND=ireals)         :: Re  ! mean radius of the earth

  LOGICAL loc_debug   ! debug switch; if .true. print information

!==============================================================================

! Output

  REAL (KIND=ireals)   ::    pui (*)   ! Interpolated u-component of vector
  REAL (KIND=ireals)   ::    pvi (*)   ! Interpolated v-component of vector
     ! If running on one processor and *kindex* is fully set, pui and pvi 
     ! may be defined with dimensions (kilons:kilone, kilats:kilate)
     ! in the calling program

  INTEGER (KIND=iintegers)   ::      &
    kierr     ! error flag, kierr = 0 if no error occured

!==============================================================================

! Local variables
  REAL (KIND=ireals)      ::         &
    zsn,                             & ! Hemisphere discriminator (+1 or -1)
    zx1, zx2, zx3, zx4, zx5, zx6,    & ! Local values on triangle
    za , zb , zc ,                   & ! Barycentric coordinates
    zred8,                           & ! Re/8
    zgm1_11, zgm1_12, zgm1_21, zgm1_22, &
    zgm2_11, zgm2_12, zgm2_21, zgm2_22, &
    zgueta1, zguchi1, zgueta2, zguchi2, &
    zgveta1, zgvchi1, zgveta2, zgvchi2, &
    zdpu12 , zdpu13 , zdpu23 ,          &
    zdpv12 , zdpv13 , zdpv23 

  REAL (KIND=ireals)      ::  zuih, zvih

  INTEGER (KIND=iintegers)   ::  &
    l1, l2, ld,                  & ! Loop indices
    j1, j2, jd, m1, m2,          & !
    mi1p1, mi2p1, mi1p2, mi2p2,  & ! Gridpoint indices
    mig1s, mig1e, mig2s, mig2e,  & ! core dimensions
    mp                             ! counter for pui and pvi

!==============================================================================

  ! Check input variables
  ! The index array "kindex" contains the dimensions j1, j2 and jd
  ! of the points within array "px" which are the nodes for the
  ! interpolation calculation; for quadratic interpolation there are
  ! needed two surrounding rows.
  ! Check, if these are available:

  mig1s = kig1sm2 + 2
  mig2s = kig2sm2 + 2
  mig1e = kig1ep2 - 2
  mig2e = kig2ep2 - 2

  IF ( (j1min .LT. mig1s) .OR. (j1max .GT. mig1e) .OR.     &
       (j2min .LT. mig2s) .OR. (j2max .GT. mig2e) .OR.     &
       (jdmin .LT. kids ) .OR. (jdmax .GT. kide )  )  THEN
    PRINT *,'  Error in subroutine *pp_interp2qv*, dimensions of ',   &
         'array "px" are not covering the computation domain'
    PRINT *,'  mig1s: ', mig1s,'  j1min: ', j1min
    PRINT *,'  mig1e: ', mig1e,'  j1max: ', j1max
    PRINT *,'  mig2s: ', mig2s,'  j2min: ', j2min
    PRINT *,'  mig2e: ', mig2e,'  j2max: ', j2max
    PRINT *,'  kids : ', kids ,'  jdmin: ', jdmin
    PRINT *,'  kide : ', kide ,'  jdmax: ', jdmax
    kierr = -1
    RETURN
  ENDIF

!==============================================================================

  zred8   = 0.125*Re

!==============================================================================
!
  DO ld = jdmin, jdmax     ! Loop over the diamonds
    ! Calculate the gradient of "pu" and "pv" into the eta- and chi-
    ! direction
    CALL grad2v (pu(kig1sm2,kig2sm2,ld), pv(kig1sm2,kig2sm2,ld),            &
                 kig1sm2 , kig1ep2 , kig2sm2 , kig2ep2 , 1, 1, pgrd,        &
                 kig1sm1 , kig1ep1 , kig2sm1 , kig2ep1 ,                    &
                 kispokes, ki1mrp  , ki2mrp  , pcpsi   , pspsi,  ld,        &
                 j1min-1 , j1max+1 , j2min-1 , j2max+1 , loc_debug,         &
                 isp11   , isp12   , isp13   , isp14   , isp15 , isp16,     &
                 isp21   , isp22   , isp23   , isp24   , isp25 , isp26,     &
                 gudeta(kig1sm1,kig2sm1,ld) ,                               &
                 gudchi(kig1sm1,kig2sm1,ld) ,                               &
                 gvdeta(kig1sm1,kig2sm1,ld) ,                               &
                 gvdchi(kig1sm1,kig2sm1,ld) ,                               &
                 kig1sm1 , kig1ep1 , kig2sm1 , kig2ep1 , kierr)

    IF (kierr .NE. 0) THEN
      PRINT *,'  Error in *pp_interp2qv* calling *grad2v* for jd=',ld
      RETURN
    ENDIF
  ENDDO

!==============================================================================

  ! Compute the auxiliary points at the midpoints of the triangle
  ! sides always taking care of the proper rotation of the vectors
  ! into the local system of the central node

  mp = 0

  DO l2 = kilats, kilate
    DO l1 = kilons, kilone

      j1    = kindex(l1,l2,1) ! first dimension of GME grid
      IF (j1 == -9999) CYCLE
      j2    = kindex(l1,l2,2) ! second dimension of GME grid
      jd    = kindex(l1,l2,3) ! number of diamond

      mp = mp+1

!------------------------------------------------------------------------------

      ! For the diamonds 1 to  5, set "zsn" to  1.
      ! For the diamonds 6 to 10, set "zsn" to -1.

      zsn = 1.
      IF (jd .GE. 6) zsn = -1.

!------------------------------------------------------------------------------

      m1    = kindex(l1,l2,4) ! Number of triangle for interpolation
      m2    = MOD (m1,6) + 1
      mi1p1 = j1 + kispoke(m1)   ! i1-index of point 1
      mi2p1 = j2 + kispoke(m1+6) ! i2-index of point 1
      mi1p2 = j1 + kispoke(m2)   ! i1-index of point 2
      mi2p2 = j2 + kispoke(m2+6) ! i2-index of point 2

      ! The 3 barycentric coordinates of the point for which the inter-
      ! polation is required
      za    = pbaryll(l1,l2,1)
      zb    = pbaryll(l1,l2,2)
      zc    = 1. - za - zb

!==============================================================================

      ! Rotate the gradients at the two points (P1 and P2) into the local
      ! spherical system of point P0 (central node)

      ! Point m1, gradient in eta-direction, first comp. (i1-direction)
      zgm1_11 = pcpsi(j1,j2,m1)*gudeta(mi1p1,mi2p1,jd) + zsn*           &
                pspsi(j1,j2,m1)*gvdeta(mi1p1,mi2p1,jd)

      ! Point m1, gradient in eta-direction, scnd. comp. (i2-direction)
      zgm1_12 = pcpsi(j1,j2,m1)*gvdeta(mi1p1,mi2p1,jd) - zsn*           &
                pspsi(j1,j2,m1)*gudeta(mi1p1,mi2p1,jd)

      ! Point m1, gradient in chi-direction, first comp. (i1-direction)
      zgm1_21 = pcpsi(j1,j2,m1)*gudchi(mi1p1,mi2p1,jd) + zsn*           &
                pspsi(j1,j2,m1)*gvdchi(mi1p1,mi2p1,jd)

      ! Point m1, gradient in chi-direction, scnd. comp. (i2-direction)
      zgm1_22 = pcpsi(j1,j2,m1)*gvdchi(mi1p1,mi2p1,jd) - zsn*           &
                pspsi(j1,j2,m1)*gudchi(mi1p1,mi2p1,jd)

      ! Point m2, gradient in eta-direction, first comp. (i1-direction)
      zgm2_11 = pcpsi(j1,j2,m2)*gudeta(mi1p2,mi2p2,jd) + zsn*           &
                pspsi(j1,j2,m2)*gvdeta(mi1p2,mi2p2,jd)

      ! Point m2, gradient in eta-direction, scnd. comp. (i2-direction)
      zgm2_12 = pcpsi(j1,j2,m2)*gvdeta(mi1p2,mi2p2,jd) - zsn*           &
                pspsi(j1,j2,m2)*gudeta(mi1p2,mi2p2,jd)

      ! Point m2, gradient in chi-direction, first comp. (i1-direction)
      zgm2_21 = pcpsi(j1,j2,m2)*gudchi(mi1p2,mi2p2,jd) + zsn*           &
                pspsi(j1,j2,m2)*gvdchi(mi1p2,mi2p2,jd)

      ! Point m2, gradient in chi-direction, scnd. comp. (i2-direction)
      zgm2_22 = pcpsi(j1,j2,m2)*gvdchi(mi1p2,mi2p2,jd) - zsn*           &
                pspsi(j1,j2,m2)*gudchi(mi1p2,mi2p2,jd)

!==============================================================================

      ! Point m1, gradient of u in the eta-direction
      zgueta1 = pcpsi(j1,j2,m1)*zgm1_11 + zsn* pspsi(j1,j2,m1)*zgm1_21

      ! Point m1, gradient of u in the chi-direction
      zguchi1 = pcpsi(j1,j2,m1)*zgm1_21 - zsn* pspsi(j1,j2,m1)*zgm1_11

      ! Point m2, gradient of u in the eta-direction
      zgueta2 = pcpsi(j1,j2,m2)*zgm2_11 + zsn* pspsi(j1,j2,m2)*zgm2_21

      ! Point m2, gradient of u in the chi-direction
      zguchi2 = pcpsi(j1,j2,m2)*zgm2_21 - zsn* pspsi(j1,j2,m2)*zgm2_11

      ! Point m1, gradient of v in the eta-direction
      zgveta1 = pcpsi(j1,j2,m1)*zgm1_12 + zsn* pspsi(j1,j2,m1)*zgm1_22

      ! Point m1, gradient of v in the chi-direction
      zgvchi1 = pcpsi(j1,j2,m1)*zgm1_22 - zsn* pspsi(j1,j2,m1)*zgm1_12

      ! Point m2, gradient of v in the eta-direction
      zgveta2 = pcpsi(j1,j2,m2)*zgm2_12 + zsn* pspsi(j1,j2,m2)*zgm2_22

      ! Point m2, gradient of v in the chi-direction
      zgvchi2 = pcpsi(j1,j2,m2)*zgm2_22 - zsn* pspsi(j1,j2,m2)*zgm2_12

!==============================================================================

      !     Interpolate the increments at the midpoints of the triangle sides
      !     using the gradients at the three triangle vertices
      zdpu12  = (gudeta(j1,j2,jd) - zgueta1)*peta(j1,j2,m1) + zsn*      &
                (gudchi(j1,j2,jd) - zguchi1)*pchi(j1,j2,m1)

      zdpu13  = (gudeta(j1,j2,jd) - zgueta2)*peta(j1,j2,m2) + zsn*      &
                (gudchi(j1,j2,jd) - zguchi2)*pchi(j1,j2,m2)

      zdpu23  = (zgueta1-zgueta2)* (peta(j1,j2,m2) - peta(j1,j2,m1)) +   &
                (zguchi1-zguchi2)* (pchi(j1,j2,m2) - pchi(j1,j2,m1))*zsn

      zdpv12  = (gvdeta(j1,j2,jd) - zgveta1)*peta(j1,j2,m1) + zsn*      &
                (gvdchi(j1,j2,jd) - zgvchi1)*pchi(j1,j2,m1)

      zdpv13  = (gvdeta(j1,j2,jd) - zgveta2)*peta(j1,j2,m2) + zsn*      &
                (gvdchi(j1,j2,jd) - zgvchi2)*pchi(j1,j2,m2)

      zdpv23  = (zgveta1-zgveta2)* (peta(j1,j2,m2) - peta(j1,j2,m1)) +   &
                (zgvchi1-zgvchi2)* (pchi(j1,j2,m2) - pchi(j1,j2,m1))*zsn

!==============================================================================

      ! Determine the values at the vertices of the triangle (# 1, 2, 3)
      ! and at the midpoints of the triangle sides (# 4, 5, 6)
      ! First for "u"
      zx1     = pu(j1,j2,jd)
      zx2     = pcpsi(j1,j2,m1)*pu(mi1p1,mi2p1,jd) + zsn*                &
                pspsi(j1,j2,m1)*pv(mi1p1,mi2p1,jd)
      zx3     = pcpsi(j1,j2,m2)*pu(mi1p2,mi2p2,jd) + zsn*                &
                pspsi(j1,j2,m2)*pv(mi1p2,mi2p2,jd)
      zx4     = 0.5*(zx1 + zx2) + zred8*zdpu12
      zx5     = 0.5*(zx2 + zx3) + zred8*zdpu23
      zx6     = 0.5*(zx1 + zx3) + zred8*zdpu13

      ! Interpolated value using the barycentric coordinates
      zuih    = zx1*zc*(2.*zc - 1.) + zx2*za*(2.*za - 1.) +              &
                zx3*zb*(2.*zb - 1.) + 4.*((zx4*za*zc + zx5*za*zb) + zx6*zb*zc)

!==============================================================================

      ! Now for "v"
      zx1     = pv(j1,j2,jd)
      zx2     = pcpsi(j1,j2,m1)*pv(mi1p1,mi2p1,jd) - zsn*                &
                pspsi(j1,j2,m1)*pu(mi1p1,mi2p1,jd)
      zx3     = pcpsi(j1,j2,m2)*pv(mi1p2,mi2p2,jd) - zsn*                &
                pspsi(j1,j2,m2)*pu(mi1p2,mi2p2,jd)
      zx4     = 0.5*(zx1 + zx2) + zred8*zdpv12
      zx5     = 0.5*(zx2 + zx3) + zred8*zdpv23
      zx6     = 0.5*(zx1 + zx3) + zred8*zdpv13

      ! Interpolated value using the barycentric coordinates
      zvih    = zx1*zc*(2.*zc - 1.) + zx2*za*(2.*za - 1.) +              &
                zx3*zb*(2.*zb - 1.) + 4.*((zx4*za*zc + zx5*za*zb) + zx6*zb*zc)

      ! Rotate back from home node to interpolation point l1,l2
      pui(mp) = protang(l1,l2,1)*zuih + protang(l1,l2,2)*zvih
      pvi(mp) = protang(l1,l2,1)*zvih - protang(l1,l2,2)*zuih

    ENDDO
  ENDDO

!==============================================================================

  kierr = 0

!==============================================================================

  ! Print debug information, if required
  IF (loc_debug) THEN
    PRINT *,'  *interp2qv* completed '
  ENDIF

!=======================================================================

END SUBROUTINE pp_interp2qv

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE grad2s (px      , kipx1s , kipx1e , kipx2s , kipx2e ,        &
                   ki3s    , ki3e   ,                                   &
                   pgrd    , kig1sm1, kig1ep1, kig2sm1, kig2ep1,        &
                   kispokes, ki1mrp , ki2mrp , kjd    ,                 &
                   ki1sc   , ki1ec  , ki2sc  , ki2ec  , loc_debug ,     &
                   isp11   , isp12  , isp13  , isp14  , isp15  , isp16, &
                   isp21   , isp22  , isp23  , isp24  , isp25  , isp26, &
                   pgeta   , pgchi  ,                                   &
                   kipg1s  , kipg1e , kipg2s , kipg2e , kierr)

!------------------------------------------------------------------------------
!
! Description:
!   *grad2s* calculates the gradient into the eta-  and chi-
!   direction for the diamond "kjd" and the layers/levels from "ki3s"
!   up to "ki3e" for the the scalar field "px".
!   Horizontally, the calculation is performed from "ki1sc" to "ki1ec"
!   in the i1-direction and from "ki2sc" to "ki2ec" in i2-direction.
!
! Method:
!
!==============================================================================
!
! Declarations:
!
! Input
  INTEGER (KIND=iintegers)   ::        &
    kipx1s , & ! first  dimension of array "px",  start index
    kipx1e , & ! first  dimension of array "px",  end   index
    kipx2s , & ! second dimension of array "px",  start index
    kipx2e , & ! second dimension of array "px",  end   index
    ki3s   , & ! third  dimension (# of layers/levels), start
    ki3e   , & ! third  dimension (# of layers/levels), end
    kig1sm1, & ! first  dimension of array "pgrd", start index
    kig1ep1, & ! first  dimension of array "pgrd", end   index
    kig2sm1, & ! second dimension of array "pgrd", start index
    kig2ep1, & ! second dimension of array "pgrd", end   index
    kjd    , & ! index  of actual diamonds
    ki1sc  , & ! first  dimension of calculation, start index
    ki1ec  , & ! first  dimension of calculation, end   index
    ki2sc  , & ! second dimension of calculation, start index
    ki2ec      ! second dimension of calculation, end   index

  REAL    (KIND=ireals)      ::      &
    px   (kipx1s :kipx1e,  kipx2s :kipx2e,  ki3s:ki3e),                 &
               ! Field of which the gradient has to be calculated
    pgrd (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7, 2)
               ! Coefficients of the gradient operator using the 6 (5)
               ! surrounding neighbours

  INTEGER (KIND=iintegers)   ::      &
    kispokes(12,4), & ! Spokes of the 4 mirrored points of the extended array
    ki1mrp (8),     & ! i1-indices of the 4 mirrored points
    ki2mrp (8)        ! i2-indices of the 4 mirrored points

  INTEGER (KIND=iintegers)   ::     &
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
    isp26          !

  LOGICAL loc_debug      ! debug switch; if .true. print information

!==============================================================================

! Output

  INTEGER (KIND=iintegers)   ::      &
    kipg1s ,  & ! first  dimension of array "pg",  start index
    kipg1e ,  & ! first  dimension of array "pg",  end   index
    kipg2s ,  & ! second dimension of array "pg",  start index
    kipg2e      ! second dimension of array "pg",  end   index

  REAL    (KIND=ireals)      ::      &
    pgeta(kipg1s :kipg1e,  kipg2s :kipg2e,  ki3s:ki3e),                &
                ! Gradient of the field px in the eta-direction
    pgchi(kipg1s :kipg1e,  kipg2s :kipg2e,  ki3s:ki3e)
                ! Gradient of the field px in the chi-direction

  INTEGER (KIND=iintegers)   ::      &
    kierr     ! error flag, kierr = 0 if no error occured

!==============================================================================

! Local variables
  INTEGER (KIND=iintegers)   ::      &
    j1, j2, j3, js   ! Loop indices

  REAL    (KIND=ireals)      ::      &
    zsn   ! Hemisphere discriminator (+1 or -1)

!==============================================================================

  ! Check input variables
  IF ((ki1sc .LT. kipx1s) .OR. (ki1ec .GT. kipx1e) .OR.        &
      (ki2sc .LT. kipx2s) .OR. (ki2ec .GT. kipx2e)) THEN
    PRINT *,'  Error in subroutine *grad2s*, dimensions of input',     &
            ' array "px" are not covering the computation domain'
    PRINT *,'  ki1sc: ', ki1sc,'  kipx1s: ', kipx1s
    PRINT *,'  ki1ec: ', ki1ec,'  kipx1e: ', kipx1e
    PRINT *,'  ki2sc: ', ki2sc,'  kipx2s: ', kipx2s
    PRINT *,'  ki2ec: ', ki2ec,'  kipx2e: ', kipx2e
    kierr = -1
    RETURN
  ENDIF

!==============================================================================

  ! Calculate the derivative of px
  ! For the diamonds 1 to  5, set "zsn" to  1.
  ! For the diamonds 6 to 10, set "zsn" to -1.

  zsn = 1.
  IF (kjd.GE.6) zsn = -1.

  DO j3 = ki3s, ki3e   ! Loop over the layers/levels
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        pgeta(j1,j2,j3) =                                          &
            pgrd(j1,j2,1,1)*px(j1      ,j2      ,j3) +             &
            pgrd(j1,j2,2,1)*px(j1+isp11,j2+isp21,j3) +             &
            pgrd(j1,j2,3,1)*px(j1+isp12,j2+isp22,j3) +             &
            pgrd(j1,j2,4,1)*px(j1+isp13,j2+isp23,j3) +             &
            pgrd(j1,j2,5,1)*px(j1+isp14,j2+isp24,j3) +             &
            pgrd(j1,j2,6,1)*px(j1+isp15,j2+isp25,j3) +             &
            pgrd(j1,j2,7,1)*px(j1+isp16,j2+isp26,j3)

        pgchi(j1,j2,j3) = zsn*                                     &
           (pgrd(j1,j2,1,2)*px(j1      ,j2      ,j3) +             &
            pgrd(j1,j2,2,2)*px(j1+isp11,j2+isp21,j3) +             &
            pgrd(j1,j2,3,2)*px(j1+isp12,j2+isp22,j3) +             &
            pgrd(j1,j2,4,2)*px(j1+isp13,j2+isp23,j3) +             &
            pgrd(j1,j2,5,2)*px(j1+isp14,j2+isp24,j3) +             &
            pgrd(j1,j2,6,2)*px(j1+isp15,j2+isp25,j3) +             &
            pgrd(j1,j2,7,2)*px(j1+isp16,j2+isp26,j3))

      ENDDO
    ENDDO
  ENDDO   ! End of loop over the layers

!==============================================================================

  ! If the gradient is requested for the extended array, correct the
  ! four mirrored points which need special spokes
  IF ( ((ki1sc .EQ. kig1sm1) .OR. (ki1ec .EQ. kig1ep1)) .AND.         &
       ((ki2sc .EQ. kig2sm1) .OR. (ki2ec .EQ. kig2ep1)) ) THEN

    DO js = 1,4            ! Loop over the four mirrored points
      j1 = ki1mrp(js)
      j2 = ki2mrp(js)
      IF(j1<-10 .OR. j2<-10) CYCLE

      IF ( (j1 .GE. ki1sc) .AND. (j1 .LE. ki1ec) .AND.                &
           (j2 .GE. ki2sc) .AND. (j2 .LE. ki2ec) ) THEN

        DO j3 = ki3s, ki3e   ! Loop over the layers/levels
          pgeta(j1,j2,j3) =                                                 &
             pgrd(j1,j2,1,1)*px(j1               ,j2                ,j3)+   &
             pgrd(j1,j2,2,1)*px(j1+kispokes(1,js),j2+kispokes( 7,js),j3)+   &
             pgrd(j1,j2,3,1)*px(j1+kispokes(2,js),j2+kispokes( 8,js),j3)+   &
             pgrd(j1,j2,4,1)*px(j1+kispokes(3,js),j2+kispokes( 9,js),j3)+   &
             pgrd(j1,j2,5,1)*px(j1+kispokes(4,js),j2+kispokes(10,js),j3)+   &
             pgrd(j1,j2,6,1)*px(j1+kispokes(5,js),j2+kispokes(11,js),j3)+   &
             pgrd(j1,j2,7,1)*px(j1+kispokes(6,js),j2+kispokes(12,js),j3)

          pgchi(j1,j2,j3) = zsn*                                            &
            (pgrd(j1,j2,1,2)*px(j1               ,j2                ,j3)+   &
             pgrd(j1,j2,2,2)*px(j1+kispokes(1,js),j2+kispokes( 7,js),j3)+   &
             pgrd(j1,j2,3,2)*px(j1+kispokes(2,js),j2+kispokes( 8,js),j3)+   &
             pgrd(j1,j2,4,2)*px(j1+kispokes(3,js),j2+kispokes( 9,js),j3)+   &
             pgrd(j1,j2,5,2)*px(j1+kispokes(4,js),j2+kispokes(10,js),j3)+   &
             pgrd(j1,j2,6,2)*px(j1+kispokes(5,js),j2+kispokes(11,js),j3)+   &
             pgrd(j1,j2,7,2)*px(j1+kispokes(6,js),j2+kispokes(12,js),j3))

          ! Mirrored point has the same gradient
          pgeta(ki1mrp(js+4),ki2mrp(js+4),j3) = pgeta(j1,j2,j3)
          pgchi(ki1mrp(js+4),ki2mrp(js+4),j3) = pgchi(j1,j2,j3)
        ENDDO     ! End of loop over layers/levels

      ENDIF

    ENDDO       ! End of loop over the four mirrored points
  ENDIF

!==============================================================================

  ! PRINT debug information, if required
  IF (loc_debug) THEN
    PRINT *,'  *grad2s* completed for diamond #: ', kjd
  ENDIF

!==============================================================================

  kierr = 0

!==============================================================================

END SUBROUTINE grad2s

!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE grad2v (pu      , pv     ,                                   &
                   kipx1s  , kipx1e , kipx2s , kipx2e ,                 &
                   ki3s    , ki3e   ,                                   &
                   pgrd    , kig1sm1, kig1ep1, kig2sm1, kig2ep1,        &
                   kispokes, ki1mrp , ki2mrp ,                          &
                   pcpsi   , pspsi  , kjd    ,                          &
                   ki1sc   , ki1ec  , ki2sc  , ki2ec  , loc_debug ,     &
                   isp11   , isp12  , isp13  , isp14  , isp15  , isp16, &
                   isp21   , isp22  , isp23  , isp24  , isp25  , isp26, &
                   pgudeta , pgudchi, pgvdeta, pgvdchi,                 &
                   kipg1s  , kipg1e , kipg2s , kipg2e , kierr)

!------------------------------------------------------------------------------
!
! Description:
!     *grad2v* calculates the gradients into the eta- and chi-direction
!     of the wind component "pu" (zonal wind) and of "pv" (meri-
!     dional wind). The calculation is performed for the diamond "kjd"
!     and the layers/levels from "ki3s" up to "ki3e".
!
! Method:
!
!==============================================================================

! Input
  INTEGER (KIND=iintegers)   ::        &
    kipx1s , & ! first  dimension of array "px",  start index
    kipx1e , & ! first  dimension of array "px",  end   index
    kipx2s , & ! second dimension of array "px",  start index
    kipx2e , & ! second dimension of array "px",  end   index
    ki3s   , & ! third  dimension (# of layers/levels), start
    ki3e   , & ! third  dimension (# of layers/levels), end
    kig1sm1, & ! first  dimension of array "pgrd", start index
    kig1ep1, & ! first  dimension of array "pgrd", end   index
    kig2sm1, & ! second dimension of array "pgrd", start index
    kig2ep1, & ! second dimension of array "pgrd", end   index
    kjd    , & ! index  of actual diamonds
    ki1sc  , & ! first  dimension of calculation, start index
    ki1ec  , & ! first  dimension of calculation, end   index
    ki2sc  , & ! second dimension of calculation, start index
    ki2ec      ! second dimension of calculation, end   index

  REAL    (KIND=ireals)      ::        &
    pu   (kipx1s :kipx1e,  kipx2s :kipx2e,  ki3s:ki3e),                    &
               ! u-component of wind field
    pv   (kipx1s :kipx1e,  kipx2s :kipx2e,  ki3s:ki3e),                    &
               ! v-component of wind field
    pcpsi(kig1sm1:kig1ep1, kig2sm1:kig2ep1, 6),                            &
               ! cosine of the rotation angle for the wind rotation
    pspsi(kig1sm1:kig1ep1, kig2sm1:kig2ep1, 6),                            &
               ! sine   of the rotation angle for the wind rotation
    pgrd (kig1sm1:kig1ep1, kig2sm1:kig2ep1, 7, 2)
               ! coefficients of the gradient operator using the 6 (5)
               ! surrounding neighbours

  INTEGER (KIND=iintegers)   ::        &
    kispokes(12,4), & ! Spokes of the 4 mirrored points of the extended array
    ki1mrp (8),     & ! i1-indices of the 4 mirrored points
    ki2mrp (8)        ! i2-indices of the 4 mirrored points

  INTEGER (KIND=iintegers)   ::     &
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
    isp26          !

  LOGICAL loc_debug      ! debug switch; if .true. print information

!==============================================================================

! Output
  INTEGER (KIND=iintegers)   ::        &
    kipg1s , & ! first  dimension of array "pg*",  start index
    kipg1e , & ! first  dimension of array "pg*",  end   index
    kipg2s , & ! second dimension of array "pg*",  start index
    kipg2e     ! second dimension of array "pg*",  end   index

  REAL    (KIND=ireals)      ::        &
    pgudeta (kipg1s :kipg1e,  kipg2s :kipg2e,  ki3s:ki3e),   &
               ! Gradient of u in eta-direction
    pgudchi (kipg1s :kipg1e,  kipg2s :kipg2e,  ki3s:ki3e),   &
               ! Gradient of u in chi-direction
    pgvdeta (kipg1s :kipg1e,  kipg2s :kipg2e,  ki3s:ki3e),   &
               ! Gradient of v in eta-direction
    pgvdchi (kipg1s :kipg1e,  kipg2s :kipg2e,  ki3s:ki3e)
               ! Gradient of v in chi-direction

  INTEGER (KIND=iintegers)   ::        &
    kierr     ! error flag, kierr = 0 if no error occured

!==============================================================================

! Local variables
  INTEGER j1, j2, j3, js, js1, js2  ! Loop indices

  REAL    (KIND=ireals)      ::        &
    zsn,                          & ! Hemisphere discriminator (+1 or -1)
    zu1, zu2, zu3, zu4, zu5, zu6, & ! Rotated u-component
    zv1, zv2, zv3, zv4, zv5, zv6    ! Rotated v-component

!==============================================================================

  ! Check input variables
  IF ((ki1sc .LT. kipx1s) .OR. (ki1ec .GT. kipx1e) .OR.     &
      (ki2sc .LT. kipx2s) .OR. (ki2ec .GT. kipx2e)) THEN
    PRINT *,'  Error in subroutine *grad2v*, dimensions of input ',     &
            'array "pu, pv" are not covering the computation domain'
    PRINT *,'  ki1sc: ', ki1sc,'  kipx1s: ', kipx1s
    PRINT *,'  ki1ec: ', ki1ec,'  kipx1e: ', kipx1e
    PRINT *,'  ki2sc: ', ki2sc,'  kipx2s: ', kipx2s
    PRINT *,'  ki2ec: ', ki2ec,'  kipx2e: ', kipx2e
    kierr = -1
    RETURN
  ENDIF

!==============================================================================

  ! Calculate the derivative of u and v in eta- and chi-direction
  ! For the diamonds 1 to  5, set "zsn" to  1.
  ! For the diamonds 6 to 10, set "zsn" to -1.
  zsn = 1.
  IF (kjd .GE. 6) zsn = -1.

  ! Gradient of zonal wind u
  DO j3 = ki3s, ki3e   ! Loop over the layers/levels
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec

        ! The term at the central node
        pgudeta(j1,j2,j3) = pgrd(j1,j2,1,1)*pu(j1,j2,j3)
        pgudchi(j1,j2,j3) = pgrd(j1,j2,1,2)*pu(j1,j2,j3)*zsn

        ! Rotate the wind component u of the neighbours into the direction
        ! of the local spherical system
        zu1 = pcpsi(j1,j2,1)*pu(j1+isp11,j2+isp21,j3) +          &
              pspsi(j1,j2,1)*pv(j1+isp11,j2+isp21,j3)*zsn
        zu2 = pcpsi(j1,j2,2)*pu(j1+isp12,j2+isp22,j3) +          &
              pspsi(j1,j2,2)*pv(j1+isp12,j2+isp22,j3)*zsn
        zu3 = pcpsi(j1,j2,3)*pu(j1+isp13,j2+isp23,j3) +          &
              pspsi(j1,j2,3)*pv(j1+isp13,j2+isp23,j3)*zsn
        zu4 = pcpsi(j1,j2,4)*pu(j1+isp14,j2+isp24,j3) +          &
              pspsi(j1,j2,4)*pv(j1+isp14,j2+isp24,j3)*zsn
        zu5 = pcpsi(j1,j2,5)*pu(j1+isp15,j2+isp25,j3) +          &
              pspsi(j1,j2,5)*pv(j1+isp15,j2+isp25,j3)*zsn
        zu6 = pcpsi(j1,j2,6)*pu(j1+isp16,j2+isp26,j3) +          &
              pspsi(j1,j2,6)*pv(j1+isp16,j2+isp26,j3)*zsn

        ! Compute the gradient into eta- and chi-direction
        pgudeta(j1,j2,j3) = pgudeta(j1,j2,j3) +                            &
                       pgrd   (j1,j2,2,1)*zu1 + pgrd   (j1,j2,3,1)*zu2 +   &
                       pgrd   (j1,j2,4,1)*zu3 + pgrd   (j1,j2,5,1)*zu4 +   &
                       pgrd   (j1,j2,6,1)*zu5 + pgrd   (j1,j2,7,1)*zu6 

        pgudchi(j1,j2,j3) = pgudchi(j1,j2,j3) +                            &
                      (pgrd   (j1,j2,2,2)*zu1 + pgrd   (j1,j2,3,2)*zu2 +   &
                       pgrd   (j1,j2,4,2)*zu3 + pgrd   (j1,j2,5,2)*zu4 +   &
                       pgrd   (j1,j2,6,2)*zu5 + pgrd   (j1,j2,7,2)*zu6)*zsn

      ENDDO
    ENDDO
  ENDDO   ! End of loop over the layers

!==============================================================================

  !     If the gradient is requested for the extended array, correct the
  !     four mirrored points which need special spokes

  IF ( ((ki1sc .EQ. kig1sm1) .OR. (ki1ec .EQ. kig1ep1)) .AND.         &
       ((ki2sc .EQ. kig2sm1) .OR. (ki2ec .EQ. kig2ep1)) ) THEN

    DO js = 1,4            ! Loop over the four mirrored points
      j1    = ki1mrp(js)
      j2    = ki2mrp(js)

      IF ( (j1 .GE. ki1sc) .AND. (j1 .LE. ki1ec) .AND.                &
           (j2 .GE. ki2sc) .AND. (j2 .LE. ki2ec) ) THEN

        DO j3 = ki3s, ki3e   ! Loop over the layers/levels

          !  The term at the central node
          pgudeta(j1,j2,j3) = pgrd(j1,j2,1,1)*pu(j1,j2,j3)
          pgudchi(j1,j2,j3) = pgrd(j1,j2,1,2)*pu(j1,j2,j3)*zsn

          ! Rotate the wind component u of the neighbours into the direction
          ! of the local spherical system
          js1 = j1+kispokes( 1,js)
          js2 = j2+kispokes( 7,js)
          zu1 = pcpsi(j1,j2,1)*pu(js1,js2,j3) +        &
                pspsi(j1,j2,1)*pv(js1,js2,j3)*zsn
          js1 = j1+kispokes( 2,js)
          js2 = j2+kispokes( 8,js)
          zu2 = pcpsi(j1,j2,2)*pu(js1,js2,j3) +        &
                pspsi(j1,j2,2)*pv(js1,js2,j3)*zsn
          js1 = j1+kispokes( 3,js)
          js2 = j2+kispokes( 9,js)
          zu3 = pcpsi(j1,j2,3)*pu(js1,js2,j3) +        &
                pspsi(j1,j2,3)*pv(js1,js2,j3)*zsn
          js1 = j1+kispokes( 4,js)
          js2 = j2+kispokes(10,js)
          zu4 = pcpsi(j1,j2,4)*pu(js1,js2,j3) +        &
                pspsi(j1,j2,4)*pv(js1,js2,j3)*zsn
          js1 = j1+kispokes( 5,js)
          js2 = j2+kispokes(11,js)
          zu5 = pcpsi(j1,j2,5)*pu(js1,js2,j3) +        &
                pspsi(j1,j2,5)*pv(js1,js2,j3)*zsn
          js1 = j1+kispokes( 6,js)
          js2 = j2+kispokes(12,js)
          zu6 = pcpsi(j1,j2,6)*pu(js1,js2,j3) +        &
                pspsi(j1,j2,6)*pv(js1,js2,j3)*zsn

          ! Compute the gradient into eta- and chi-direction
          pgudeta(j1,j2,j3) = pgudeta(j1,j2,j3) +                         &
                         pgrd (j1,j2,2,1)*zu1 + pgrd (j1,j2,3,1)*zu2 +    &
                         pgrd (j1,j2,4,1)*zu3 + pgrd (j1,j2,5,1)*zu4 +    &
                         pgrd (j1,j2,6,1)*zu5 + pgrd (j1,j2,7,1)*zu6

          pgudchi(j1,j2,j3) = pgudchi(j1,j2,j3) +                         &
                        (pgrd (j1,j2,2,2)*zu1 + pgrd (j1,j2,3,2)*zu2 +    &
                         pgrd (j1,j2,4,2)*zu3 + pgrd (j1,j2,5,2)*zu4 +    &
                         pgrd (j1,j2,6,2)*zu5 + pgrd (j1,j2,7,2)*zu6)*zsn

          ! Mirrored point has the same gradient
          pgudeta(ki1mrp(js+4),ki2mrp(js+4),j3) = pgudeta(j1,j2,j3)
          pgudchi(ki1mrp(js+4),ki2mrp(js+4),j3) = pgudchi(j1,j2,j3)

        ENDDO  ! End of loop over the layers

      ENDIF

    ENDDO    ! End of loop over the four special points

  ENDIF

!==============================================================================

  ! Gradients of meridional wind v
  DO j3 = ki3s, ki3e   ! Loop over the layers/levels
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec

        ! The term at the central node
        pgvdeta(j1,j2,j3) = pgrd(j1,j2,1,1)*pv(j1,j2,j3)
        pgvdchi(j1,j2,j3) = pgrd(j1,j2,1,2)*pv(j1,j2,j3)*zsn

        ! Rotate the wind component v of the neighbours into the direction
        ! of the local spherical system
        zv1 = pcpsi(j1,j2,1)*pv(j1+isp11,j2+isp21,j3) -                   &
              pspsi(j1,j2,1)*pu(j1+isp11,j2+isp21,j3)*zsn
        zv2 = pcpsi(j1,j2,2)*pv(j1+isp12,j2+isp22,j3) -                   &
              pspsi(j1,j2,2)*pu(j1+isp12,j2+isp22,j3)*zsn
        zv3 = pcpsi(j1,j2,3)*pv(j1+isp13,j2+isp23,j3) -                   &
              pspsi(j1,j2,3)*pu(j1+isp13,j2+isp23,j3)*zsn
        zv4 = pcpsi(j1,j2,4)*pv(j1+isp14,j2+isp24,j3) -                   &
              pspsi(j1,j2,4)*pu(j1+isp14,j2+isp24,j3)*zsn
        zv5 = pcpsi(j1,j2,5)*pv(j1+isp15,j2+isp25,j3) -                   &
              pspsi(j1,j2,5)*pu(j1+isp15,j2+isp25,j3)*zsn
        zv6 = pcpsi(j1,j2,6)*pv(j1+isp16,j2+isp26,j3) -                   &
              pspsi(j1,j2,6)*pu(j1+isp16,j2+isp26,j3)*zsn

        ! Compute the gradient into eta- and chi-direction
        pgvdeta(j1,j2,j3) = pgvdeta(j1,j2,j3) +                           &
                         pgrd (j1,j2,2,1)*zv1 + pgrd (j1,j2,3,1)*zv2 +    &
                         pgrd (j1,j2,4,1)*zv3 + pgrd (j1,j2,5,1)*zv4 +    &
                         pgrd (j1,j2,6,1)*zv5 + pgrd (j1,j2,7,1)*zv6

        pgvdchi(j1,j2,j3) = pgvdchi(j1,j2,j3) +                           &
                        (pgrd (j1,j2,2,2)*zv1 + pgrd (j1,j2,3,2)*zv2 +    &
                         pgrd (j1,j2,4,2)*zv3 + pgrd (j1,j2,5,2)*zv4 +    &
                         pgrd (j1,j2,6,2)*zv5 + pgrd (j1,j2,7,2)*zv6)*zsn

      ENDDO
    ENDDO
  ENDDO   ! End of loop over the layers

!==============================================================================

  ! If the gradient is requested for the extended array, correct the
  ! four mirrored points which need special spokes
  IF ( ((ki1sc .EQ. kig1sm1) .OR. (ki1ec .EQ. kig1ep1)) .AND.       &
       ((ki2sc .EQ. kig2sm1) .OR. (ki2ec .EQ. kig2ep1)) ) THEN

    DO js = 1,4            ! Loop over the four mirrored points

      j1    = ki1mrp(js)
      j2    = ki2mrp(js)

      IF ( (j1 .GE. ki1sc) .AND. (j1 .LE. ki1ec) .AND.              &
           (j2 .GE. ki2sc) .AND. (j2 .LE. ki2ec) ) THEN

        DO j3 = ki3s, ki3e   ! Loop over the layers/levels

          ! The term at the central node
          pgvdeta(j1,j2,j3) = pgrd(j1,j2,1,1)*pv(j1,j2,j3)
          pgvdchi(j1,j2,j3) = pgrd(j1,j2,1,2)*pv(j1,j2,j3)*zsn

          ! Rotate the wind component v of the neighbours into the direction
          ! of the local spherical system
          js1 = j1+kispokes( 1,js)
          js2 = j2+kispokes( 7,js)
          zv1 = pcpsi(j1,j2,1)*pv(js1,js2,j3) -                     &
                pspsi(j1,j2,1)*pu(js1,js2,j3)*zsn
          js1 = j1+kispokes( 2,js)
          js2 = j2+kispokes( 8,js)
          zv2 = pcpsi(j1,j2,2)*pv(js1,js2,j3) -                     &
                pspsi(j1,j2,2)*pu(js1,js2,j3)*zsn
          js1 = j1+kispokes( 3,js)
          js2 = j2+kispokes( 9,js)
          zv3 = pcpsi(j1,j2,3)*pv(js1,js2,j3) -                     &
                pspsi(j1,j2,3)*pu(js1,js2,j3)*zsn
          js1 = j1+kispokes( 4,js)
          js2 = j2+kispokes(10,js)
          zv4 = pcpsi(j1,j2,4)*pv(js1,js2,j3) -                     &
                pspsi(j1,j2,4)*pu(js1,js2,j3)*zsn
          js1 = j1+kispokes( 5,js)
          js2 = j2+kispokes(11,js)
          zv5 = pcpsi(j1,j2,5)*pv(js1,js2,j3) -                     &
                pspsi(j1,j2,5)*pu(js1,js2,j3)*zsn
          js1 = j1+kispokes( 6,js)
          js2 = j2+kispokes(12,js)
          zv6 = pcpsi(j1,j2,6)*pv(js1,js2,j3) -                     &
                pspsi(j1,j2,6)*pu(js1,js2,j3)*zsn

          ! Compute the gradient into eta- and chi-direction
          pgvdeta(j1,j2,j3) = pgvdeta(j1,j2,j3) +                          &
                           pgrd (j1,j2,2,1)*zv1 + pgrd (j1,j2,3,1)*zv2 +   &
                           pgrd (j1,j2,4,1)*zv3 + pgrd (j1,j2,5,1)*zv4 +   &
                           pgrd (j1,j2,6,1)*zv5 + pgrd (j1,j2,7,1)*zv6

          pgvdchi(j1,j2,j3) = pgvdchi(j1,j2,j3) +                          &
                          (pgrd (j1,j2,2,2)*zv1 + pgrd (j1,j2,3,2)*zv2 +   &
                           pgrd (j1,j2,4,2)*zv3 + pgrd (j1,j2,5,2)*zv4 +   &
                           pgrd (j1,j2,6,2)*zv5 + pgrd (j1,j2,7,2)*zv6)*zsn

          ! Mirrored point has the same gradient
          pgvdeta(ki1mrp(js+4),ki2mrp(js+4),j3) = pgvdeta(j1,j2,j3)
          pgvdchi(ki1mrp(js+4),ki2mrp(js+4),j3) = pgvdchi(j1,j2,j3)

        ENDDO  ! End of loop over the layers

      ENDIF

    ENDDO    ! End of loop over the four special points

  ENDIF

!==============================================================================

      kierr = 0

!==============================================================================

  ! Print debug information, if required
  IF (loc_debug) THEN
    PRINT *,'  *grad2v* completed for diamond #: ', kjd
  ENDIF

!=======================================================================

END SUBROUTINE grad2v

!==============================================================================

END MODULE gme_utilities
