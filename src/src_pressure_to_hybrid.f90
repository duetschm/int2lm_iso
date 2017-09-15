!+ Source Module for vertical interpolation from pressure to hybrid levels
!==============================================================================

MODULE src_pressure_to_hybrid

!==============================================================================
!
! Description:
!   This module contains routines for the vertical interpolation from 
!   pressure levels to hybrid levels defined with ak, bk, similar to the GME.
!
! Current Code Owner: DWD, Ulrich Schaettler
!    phone:  +49  69  8062 2739
!    fax:    +49  69  8062 3721
!    email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V1_14        2010/11/19 Ulrich Schaettler
!  Initial release
! V1_19        2012/06/06 Ulrich Schaettler, Burkhardt Rockel
!  Modified lower boundary condition for vertical interpolation, if surface pressure
!   exceeds the "first guess" press_level(ke_in)+5000.0
!  (from CLM: Make sure that lower pressure boundary is always greater than
!    press_level(ke_in))
!  Introduced pressure level support for climate model CM input
!   (choose a 38-level version)
!  Introduced computation of control geopotential on hybrid layers, if it has not
!   been read
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

USE data_grid_lm,    ONLY: &
kedim          !

!------------------------------------------------------------------------------

USE data_grid_in,    ONLY: &
ak_in,       & ! vertical coordinate parameters for half levels
bk_in,       & !                  - " -
akh_in,      & ! vertical coordinate parameters for main levels
bkh_in,      & !                  - " -
dak_in,      & ! difference of coordinate parameters
dbk_in,      & !                  - " -
ie_in,       & ! ie for input grid, local domain
je_in,       & ! je for input grid, local domain
ke_in,       & ! ke for input grid
ke_hybrid,   & ! number of vertical hybrid levels
ke_pressure, & ! number of pressure levels when interpolating from pressure levels
kedim_in,    & ! MAX (ke_in, ke_hybrid)
press_level, & ! list of available pressure levels (in Pa) for GFS
lcm_pres_coor,&!
nlevskip,    & ! number of levels to skip at top of input model
ke1in          ! ke+1 for input grid

!------------------------------------------------------------------------------

USE data_fields_in,  ONLY: &
hsurf_in,    & ! orography                                   (  m  )
hhl_in,      & ! height of half levels of coarse grid        (  m  )
ps_in,       & ! surface pressure                            ( Pa  )
t_s_in,      & ! surface temperature                         (  K  )
t_in,        & ! temperature                                 (  K  )
u_in,        & ! u-component of wind                         ( m/s )
v_in,        & ! v-component of wind                         ( m/s )
qv_in,       & ! specific water vapor content                (kg/kg)
qc_in,       & ! specific water vapor content                (kg/kg)
fis_in,      & ! orography * g
fic_in         ! orography * g

!------------------------------------------------------------------------------

USE data_int2lm_control,    ONLY: &
idbg_level,   & ! to control verbosity of output
lprintdeb_all,& ! whether all PEs print debug output
lprog_qi,     & ! if .TRUE., interpolate qi to LM grid
lgfs2lm,      & !
lcm2lm,       & !
lcomp_bound,  & ! compute fields for boundaries
pcontrol_fi,  & ! pressure of control level for geopotential
qvmin,        & ! minimum value of water vapor (security)
qcmin,        & ! minimum value of cloud water (security)
qimin,        & ! minimum value of cloud ice content (security)
noutput         ! unit number for output file

!------------------------------------------------------------------------------

USE data_int2lm_io,         ONLY: &
  nvar_in,           & ! actual maximum number of variables in input variable table
  yin_form_read,     & ! input format of boundary data
  var_in               ! array for input model variable table

!------------------------------------------------------------------------------

USE data_int2lm_constants,     ONLY :  &
R_d,     & ! gas constant for dry air                      [J/K*kg]
Rdv,     & ! = R_d/R_v
Rvd_m_o, & ! = R_v/R_d - 1.0
O_m_rdv, & ! = 1. - Rdv
G,       & ! gravity at sea level                          [ms-2]
r_earth, & ! mean radius of the earth                      [m]
B1,      & !  a
B2_w,    & !  b
B2_i,    & !  b for ice
B3,      & !  c/b (0 degree Celsius [Kelvin])
B4_w,    & !  d
B4_i,    & !  d for ice
degrad, pi !

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY :  &
my_cart_id

!------------------------------------------------------------------------------

USE utilities,        ONLY :   &
tautsp2D             !

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================
!+ Organizes the vertical interpolation of multi-level variables
!------------------------------------------------------------------------------

SUBROUTINE org_vert_interpol_p2h (ierror, yerror)

!------------------------------------------------------------------------------
!
! Description:
!   org_vert_interpol_p2h organizes the vertical interpolations from the
!   fields on pressure levels to fields on hybrid levels similar to GME.
!
! Method:
!   First the ak, bk are set from the old GME 31-level version.
!   Then the interpolations are computed similar to the interpolations
!   in src_vert_interpol.
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
INTEGER  (KIND=iintegers), INTENT(OUT)  ::  &
ierror                  ! status and error status variable

CHARACTER (LEN=200), INTENT(OUT)        ::  &
yerror                  ! error message for error handling

! Local variables
INTEGER  (KIND=iintegers)  ::  &
izdebug,              & ! for debug output
kzgr,                 & ! top of boundary layer
i, j, k, n, istat, is1, js1, mzfi_loc_in

REAL    (KIND=ireals)      ::  &
zaq, zbq, zgqv, zpresm, zpgr, zpht

REAL    (KIND=ireals)      ::  &
zgrh_in (ie_in,je_in,kedim)

CHARACTER (LEN= 25)        ::  &
yzroutine   ! name of this routine for error handling

REAL (KIND=ireals) , ALLOCATABLE                     ::   &
  ztv    (:,:),          & ! virtual temperature
  zpo    (:,:),          & !
  zpu    (:,:),          & !
  zfiu   (:,:),          & !
  zfio   (:,:)             !

REAL (KIND=ireals) :: sf_psat_w, sf_qsat_ec, x, y, z, zi, v, w, wi

sf_psat_w  (x,y,z,v,w) = y * EXP(z*(x-v)/(x-w))
sf_qsat_ec (x,y,z,v)   = z * x / MAX( (y-v*x), 1.0_ireals)    ! EC2LM -version

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

ierror    = 0
yerror    = '                     '
yzroutine = 'org_vert_interpol_p2h'

IF (lprintdeb_all) THEN
  izdebug = idbg_level
ELSE
  IF (my_cart_id == 0) THEN
    izdebug = idbg_level
  ELSE
    izdebug = 0
  ENDIF
ENDIF

!------------------------------------------------------------------------------
! Section 2: Set ak and bk for a ke_hybrid (main) level version:
!------------------------------------------------------------------------------

IF (ke_hybrid == 31) THEN
  ! implemented by DWD for NCEP forecast data
  ak_in( 1) =  0.0000000000e+00_ireals;     bk_in( 1) =  0.0000000000e+00_ireals
  ak_in( 2) =  2.0000000000e+03_ireals;     bk_in( 2) =  0.0000000000e+00_ireals
  ak_in( 3) =  4.0000000000e+03_ireals;     bk_in( 3) =  0.0000000000e+00_ireals
  ak_in( 4) =  6.0000000000e+03_ireals;     bk_in( 4) =  0.0000000000e+00_ireals
  ak_in( 5) =  8.0000000000e+03_ireals;     bk_in( 5) =  0.0000000000e+00_ireals
  ak_in( 6) =  9.9761367188e+03_ireals;     bk_in( 6) =  3.9085815661e-04_ireals
  ak_in( 7) =  1.1820539062e+04_ireals;     bk_in( 7) =  2.9197006952e-03_ireals
  ak_in( 8) =  1.3431394531e+04_ireals;     bk_in( 8) =  9.1941319406e-03_ireals
  ak_in( 9) =  1.4736355469e+04_ireals;     bk_in( 9) =  2.0319156349e-02_ireals
  ak_in(10) =  1.5689207031e+04_ireals;     bk_in(10) =  3.6974858493e-02_ireals
  ak_in(11) =  1.6266609375e+04_ireals;     bk_in(11) =  5.9487640858e-02_ireals
  ak_in(12) =  1.6465003906e+04_ireals;     bk_in(12) =  8.7894976139e-02_ireals
  ak_in(13) =  1.6297621094e+04_ireals;     bk_in(13) =  1.2200361490e-01_ireals
  ak_in(14) =  1.5791597656e+04_ireals;     bk_in(14) =  1.6144150496e-01_ireals
  ak_in(15) =  1.4985269531e+04_ireals;     bk_in(15) =  2.0570325851e-01_ireals
  ak_in(16) =  1.3925519531e+04_ireals;     bk_in(16) =  2.5418859720e-01_ireals
  ak_in(17) =  1.2665292969e+04_ireals;     bk_in(17) =  3.0623537302e-01_ireals
  ak_in(18) =  1.1261230469e+04_ireals;     bk_in(18) =  3.6114501953e-01_ireals
  ak_in(19) =  9.7714062500e+03_ireals;     bk_in(19) =  4.1820228100e-01_ireals
  ak_in(20) =  8.2532109375e+03_ireals;     bk_in(20) =  4.7668814659e-01_ireals
  ak_in(21) =  6.7613398438e+03_ireals;     bk_in(21) =  5.3588658571e-01_ireals 
  ak_in(22) =  5.3459140625e+03_ireals;     bk_in(22) =  5.9508424997e-01_ireals
  ak_in(23) =  4.0507177734e+03_ireals;     bk_in(23) =  6.5356457233e-01_ireals
  ak_in(24) =  2.9115693359e+03_ireals;     bk_in(24) =  7.1059441566e-01_ireals
  ak_in(25) =  1.9548051758e+03_ireals;     bk_in(25) =  7.6540523767e-01_ireals
  ak_in(26) =  1.1958898926e+03_ireals;     bk_in(26) =  8.1716698408e-01_ireals
  ak_in(27) =  6.3814892578e+02_ireals;     bk_in(27) =  8.6495584249e-01_ireals
  ak_in(28) =  2.7162646484e+02_ireals;     bk_in(28) =  9.0771585703e-01_ireals
  ak_in(29) =  7.2063583374e+01_ireals;     bk_in(29) =  9.4421321154e-01_ireals
  ak_in(30) =  0.0000000000e+00_ireals;     bk_in(30) =  9.7298520803e-01_ireals
  ak_in(31) =  0.0000000000e+00_ireals;     bk_in(31) =  9.9228149652e-01_ireals
  ak_in(32) =  0.0000000000e+00_ireals;     bk_in(32) =  1.0000000000e+00_ireals
ELSEIF (ke_hybrid == 38) THEN
  ! implemented by CLM for NCEP reanalysis
  ak_in( 1) =  0.0000000000e+00_ireals;     bk_in( 1) =  0.0000000000e+00_ireals
  ak_in( 2) =  1.0000000000e+02_ireals;     bk_in( 2) =  0.0000000000e+00_ireals
  ak_in( 3) =  2.0000000000e+02_ireals;     bk_in( 3) =  0.0000000000e+00_ireals
  ak_in( 4) =  3.0000000000e+02_ireals;     bk_in( 4) =  0.0000000000e+00_ireals
  ak_in( 5) =  4.0000000000e+02_ireals;     bk_in( 5) =  0.0000000000e+00_ireals
  ak_in( 6) =  5.0000000000e+02_ireals;     bk_in( 6) =  0.0000000000e+00_ireals
  ak_in( 7) =  7.0000000000e+02_ireals;     bk_in( 7) =  0.0000000000e+00_ireals
  ak_in( 8) =  1.0000000000e+03_ireals;     bk_in( 8) =  0.0000000000e+00_ireals
  ak_in( 9) =  2.0000000000e+03_ireals;     bk_in( 9) =  0.0000000000e+00_ireals
  ak_in(10) =  4.0000000000e+03_ireals;     bk_in(10) =  0.0000000000e+00_ireals
  ak_in(11) =  6.0000000000e+03_ireals;     bk_in(11) =  0.0000000000e+00_ireals
  ak_in(12) =  8.0000000000e+03_ireals;     bk_in(12) =  0.0000000000e+00_ireals
  ak_in(13) =  9.9761367188e+03_ireals;     bk_in(13) =  3.9085815661e-04_ireals
  ak_in(14) =  1.1820539062e+04_ireals;     bk_in(14) =  2.9197006952e-03_ireals
  ak_in(15) =  1.3431394531e+04_ireals;     bk_in(15) =  9.1941319406e-03_ireals
  ak_in(16) =  1.4736355469e+04_ireals;     bk_in(16) =  2.0319156349e-02_ireals
  ak_in(17) =  1.5689207031e+04_ireals;     bk_in(17) =  3.6974858493e-02_ireals
  ak_in(18) =  1.6266609375e+04_ireals;     bk_in(18) =  5.9487640858e-02_ireals
  ak_in(19) =  1.6465003906e+04_ireals;     bk_in(19) =  8.7894976139e-02_ireals
  ak_in(20) =  1.6297621094e+04_ireals;     bk_in(20) =  1.2200361490e-01_ireals
  ak_in(21) =  1.5791597656e+04_ireals;     bk_in(21) =  1.6144150496e-01_ireals
  ak_in(22) =  1.4985269531e+04_ireals;     bk_in(22) =  2.0570325851e-01_ireals
  ak_in(23) =  1.3925519531e+04_ireals;     bk_in(23) =  2.5418859720e-01_ireals
  ak_in(24) =  1.2665292969e+04_ireals;     bk_in(24) =  3.0623537302e-01_ireals
  ak_in(25) =  1.1261230469e+04_ireals;     bk_in(25) =  3.6114501953e-01_ireals
  ak_in(26) =  9.7714062500e+03_ireals;     bk_in(26) =  4.1820228100e-01_ireals
  ak_in(27) =  8.2532109375e+03_ireals;     bk_in(27) =  4.7668814659e-01_ireals
  ak_in(28) =  6.7613398438e+03_ireals;     bk_in(28) =  5.3588658571e-01_ireals
  ak_in(29) =  5.3459140625e+03_ireals;     bk_in(29) =  5.9508424997e-01_ireals
  ak_in(30) =  4.0507177734e+03_ireals;     bk_in(30) =  6.5356457233e-01_ireals
  ak_in(31) =  2.9115693359e+03_ireals;     bk_in(31) =  7.1059441566e-01_ireals
  ak_in(32) =  1.9548051758e+03_ireals;     bk_in(32) =  7.6540523767e-01_ireals
  ak_in(33) =  1.1958898926e+03_ireals;     bk_in(33) =  8.1716698408e-01_ireals
  ak_in(34) =  6.3814892578e+02_ireals;     bk_in(34) =  8.6495584249e-01_ireals
  ak_in(35) =  2.7162646484e+02_ireals;     bk_in(35) =  9.0771585703e-01_ireals
  ak_in(36) =  7.2063583374e+01_ireals;     bk_in(36) =  9.4421321154e-01_ireals
  ak_in(37) =  0.0000000000e+00_ireals;     bk_in(37) =  9.7298520803e-01_ireals
  ak_in(38) =  0.0000000000e+00_ireals;     bk_in(38) =  9.9228149652e-01_ireals
  ak_in(39) =  0.0000000000e+00_ireals;     bk_in(39) =  1.0000000000e+00_ireals
ENDIF

DO k = 1, ke_hybrid
  akh_in(k) = (ak_in(k) + ak_in(k+1)) * 0.5_ireals
  bkh_in(k) = (bk_in(k) + bk_in(k+1)) * 0.5_ireals
  dak_in(k) =  ak_in(k+1) - ak_in(k)
  dbk_in(k) =  bk_in(k+1) - bk_in(k)
ENDDO

!------------------------------------------------------------------------------
! Section 3: Compute generalized relative humidity
!------------------------------------------------------------------------------

IF (izdebug > 10) THEN
  PRINT *, 'in vert_interpol_p2h: compute grh'
ENDIF

DO k = 1, ke_in
  DO j = 1, je_in
    DO i = 1, ie_in
      zpresm = press_level(k)
      zaq    = sf_psat_w   (t_in(i,j,k), b1, b2_w, b3, b4_w)
      zbq    = sf_qsat_ec  (zaq, zpresm, Rdv, O_m_rdv)
      zgqv   = 1.0_ireals / zbq
      zgrh_in(i,j,k) = (qv_in(i,j,k) + qc_in(i,j,k)) * zgqv
    ENDDO
  ENDDO
ENDDO

!------------------------------------------------------------------------------
! Section 4: Vertical interpolations
!------------------------------------------------------------------------------

IF (izdebug > 10) THEN
  PRINT *, 'in vert_interpol_p2h: vertical interpolations'
ENDIF

! Set boundary layer  height: kzgr
zpgr = 850.0E2_ireals
DO k = ke_in, 1, -1
  zpht = press_level(k)
  IF (zpht > zpgr) kzgr = k
ENDDO
kzgr = kzgr-1

CALL vert_interpol_p2h (t_in   , 't' , ke_hybrid, akh_in, bkh_in, ps_in,  &
                                       kzgr, ierror, yerror)

CALL vert_interpol_p2h (zgrh_in, 'rh', ke_hybrid, akh_in, bkh_in, ps_in,  &
                                       kzgr, ierror, yerror)

CALL vert_interpol_p2h (u_in   , 'u' , ke_hybrid, akh_in, bkh_in, ps_in,  &
                                       kzgr, ierror, yerror)

CALL vert_interpol_p2h (v_in   , 'v' , ke_hybrid, akh_in, bkh_in, ps_in,  &
                                       kzgr, ierror, yerror)

!------------------------------------------------------------------------------
! Section 5: Split generalized relative humidity again to qv, qc
!------------------------------------------------------------------------------

  DO k = 1, ke_hybrid
    DO j = 1, je_in
      DO i = 1, ie_in
        zpresm = akh_in(k) + bkh_in(k) * ps_in(i,j)
        zaq    = sf_psat_w  (t_in(i,j,k), b1, b2_w, b3, b4_w)
        zbq    = sf_qsat_ec (zaq, zpresm, Rdv, O_m_rdv)
        qv_in(i,j,k) =   MIN(1.0_ireals, zgrh_in(i,j,k)) * zbq
        qc_in(i,j,k) =   MAX(0.0_ireals, zgrh_in(i,j,k)-1.0_ireals) * zbq
        IF (qv_in(i,j,k) < qvmin) THEN
          qv_in(i,j,k)  = qvmin
        ENDIF
        IF (qc_in(i,j,k) < qcmin) THEN
          qc_in(i,j,k)  = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 6: Reset ke_in to ke_hybrid
!------------------------------------------------------------------------------

IF (lgfs2lm .OR. (lcm2lm .AND. lcm_pres_coor)) THEN
  ! this has to be done for ke_in,  but also in the variable tables for all
  ! 3D variables

  IF (izdebug > 10) THEN
    PRINT *, 'Reset vertical dimension to hybrid vertical levels:  ', ke_hybrid
  ENDIF

  ke_in = ke_hybrid

  DO n = 1, nvar_in
    IF (var_in(n)%rank == 3) var_in(n)%nlevels = ke_in
  ENDDO
ENDIF

!------------------------------------------------------------------------------
! Section 7: Compute control geopotential, if it has not been read
!------------------------------------------------------------------------------

  ! Find some special fields in var_in table
  mzfi_loc_in = 0;
  DO n = 1, nvar_in
    IF (TRIM(var_in(n)%name) == 'FI')   mzfi_loc_in = n
  ENDDO

  IF ( ((var_in(mzfi_loc_in)%dattyp(1:1) == 'I' .AND. .NOT. lcomp_bound)   .OR.  &
        (var_in(mzfi_loc_in)%dattyp(2:2) == 'B' .AND.       lcomp_bound) ) .AND. &
        ( .NOT. var_in(mzfi_loc_in)%lreadin ) .AND. (lcm_pres_coor)) THEN

    IF (izdebug >= 20) THEN
      PRINT *, '      Compute control geopotential: ', ke_in
    ENDIF

    ALLOCATE (ztv (ie_in,je_in), zpo (ie_in,je_in), zpu(ie_in,je_in),   &
              zfiu(ie_in,je_in), zfio(ie_in,je_in), STAT=ierror)

    zpu    (:,:) = ps_in  (:,:)
    zfiu   (:,:) = fis_in (:,:)
    fic_in (:,:) = 0.0

    DO k = ke_in, 1, - 1

      IF (k == 1) THEN
        zpo (:,:) = 0.0
        ztv (:,:) = t_in(:,:,k)*(1.0 + Rvd_m_o*qv_in(:,:,k))
        zfio(:,:) = zfiu(:,:) + R_d*ztv(:,:)*LOG(2.0_ireals)
      ELSE
        zpo (:,:) = ak_in(k)    + bk_in(k)*ps_in(:,:)
        ztv (:,:) = t_in(:,:,k)*(1.0 + Rvd_m_o *qv_in(:,:,k))
        zfio(:,:) = zfiu(:,:) + R_d*ztv(:,:)*LOG(zpu(:,:)/zpo(:,:))
      ENDIF

      WHERE (zpu(:,:) > pcontrol_fi .AND. zpo(:,:) <= pcontrol_fi)
        fic_in(:,:) = zfiu(:,:) + R_d*ztv(:,:)*LOG(zpu(:,:)/pcontrol_fi)
      END WHERE

      zpu (:,:) = zpo (:,:)
      zfiu(:,:) = zfio(:,:)

    ENDDO

    DEALLOCATE (ztv, zpo, zpu, zfiu, zfio)

    var_in(mzfi_loc_in)%lreadin = .TRUE.

    IF (my_cart_id == 0)  THEN
      PRINT *, ' *** Control Geopotential has been calculated ***'
      PRINT *, ' *** FI is marked as being read               ***'
    END IF

  ENDIF

IF (izdebug > 10) THEN
  PRINT *, 'End of vert_interpol_p2h'
ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE org_vert_interpol_p2h

!==============================================================================
!==============================================================================

!+ Vertical Interpolation of atmospheric variables.
!------------------------------------------------------------------------------

SUBROUTINE vert_interpol_p2h (xin, yfld, kevert, akmain, bkmain, ps_co, kgr, &
                              ierr_intv, yerr_intv)

!------------------------------------------------------------------------------
!
! Description:
!   Temperature, generalized relative humidity, and the wind components are
!   interpolated from pressure levels on to hybrid levels specified by the 
!   vertical coordinate parameters akmain, bkmain and the surface pressure
!   on the input coarse grid ps_co.
!
! Method:
!   The interpolation maintains the profiles of the boundary layer.
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
INTEGER (KIND=iintegers), INTENT(IN)         ::  &
  kevert,             & ! vertical dimension of coordinate parameters
                        ! (be careful when to use kevert and when ke_in!
                        !  kevert is only needed after calculating the
                        !  tension-splines)
  kgr                   ! index of boundary layer top

REAL (KIND=ireals),       INTENT (INOUT)     ::  &
  xin(ie_in, je_in, kedim_in)   ! field that has to be interpolated

REAL (KIND=ireals),       INTENT (IN)        ::  &
  akmain(kevert), bkmain(kevert), &! coordinate parameters
  ps_co (ie_in,je_in)              ! first adapted surface pressure

! Scalar arguments with intent(in):
CHARACTER (LEN= *),       INTENT(IN)         ::  &
  yfld                  ! name of the field

INTEGER (KIND=iintegers), INTENT(OUT)        ::  &
  ierr_intv

CHARACTER (LEN=80),       INTENT(OUT)        ::     &
  yerr_intv

!------------------------------------------------------------------------------
!
! Local variables
INTEGER (KIND=iintegers)       ::     &
  kzint_vec(ie_in), izln_vec(ie_in), izn(ie_in)

INTEGER (KIND=iintegers)       ::     &
  i, j, k, izerror, izdebug, idone

LOGICAL                        ::     &
  ldone(ie_in)

REAL    (KIND=ireals)          ::     &
  zgamma, zprod, zdpx

! Local arrays:
INTEGER (KIND=iintegers)       ::     &
  izindex (ie_in,kevert), &! indices used during interpolation
                           ! (vertical dimension of the target)
  kzdims(24)               ! vertical dimensions for boundary exchange

REAL (KIND=ireals)             ::     &
  zrhmax    (ie_in),             & !
  zpmain    (ie_in,kevert),      & !
  zxexp_vec (ie_in,ke_in+1),     & ! values for interpolations
  zpexp_vec (ie_in,ke_in+1),     & ! points where above values are valid
  zbreak_vec(ie_in,(ke_in+1)*3), & ! work array for tension spline routine
  zs_vec    (ie_in,(ke_in+1)*6), & ! work array for tension spline routine
  zcoef_vec (ie_in,4,3*(ke_in+4))  ! work array for tension spline routine

CHARACTER (LEN=25)             ::     &
  yzroutine

!- End of header
!------------------------------------------------------------------------------

  ierr_intv = 0_iintegers
  yerr_intv = '                 '
  yzroutine = 'vert_interpol_p2h'
  zgamma    =  5.5_ireals

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 1: Prepare input fields for tautsp
!------------------------------------------------------------------------------

  IF (izdebug > 19) THEN
    PRINT *, '      Vertical interpolation for ', yfld
  ENDIF
    
  DO j = 1, je_in

    DO k = 1, ke_in
      DO i = 1, ie_in
        zxexp_vec (i,k) = xin (i,j,k)
        zpexp_vec (i,k) = press_level(k)
      ENDDO
    ENDDO

    ! one layer around the surface as lower boundary condition
    ! (at the model top, it is hopefully o.k.)

    ! in most cases it is o.k. to set press_level(ke_in)+5000.0 as lower boundary
    ! condition. But sometimes the pressure might exceed this value.
    ! This has to be checked
    DO i = 1, ie_in
      IF (ps_co(i,j) <= press_level(ke_in)+5000.0_ireals) THEN
        zpexp_vec (i,ke_in+1) = press_level(ke_in) + 5000.0_ireals
        zxexp_vec (i,ke_in+1) = xin    (i,j,ke_in)
      ELSE
        ! then take available surface fields as lower boundary
        zpexp_vec (i,ke_in+1) = ps_co(i,j)
        IF     (yfld == 't') THEN
          zxexp_vec (i,ke_in+1) = t_s_in (i,j)
        ELSEIF (yfld == 'u') THEN
          zxexp_vec (i,ke_in+1) = 0.0_ireals
        ELSEIF (yfld == 'v') THEN
          zxexp_vec (i,ke_in+1) = 0.0_ireals
        ELSEIF (yfld == 'rh') THEN
          zxexp_vec (i,ke_in+1) = xin (i,j,ke_in)
        ENDIF
      ENDIF
    ENDDO

    kzint_vec (:) = ke_in + 1

!------------------------------------------------------------------------------
! Section 2: Vertical interpolation from pressure levels to hybrid levels
!------------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Section 2.1: Compute tension splines
    !------------------------------------------------------------------------

    ! Interpolation with Tension-Splines to hybrid levels
    DO i = 1, ie_in
      ! izln must be defined before !!!!
      kzint_vec(i) =  ke_in+1
      izln_vec (i) = (ke_in+1)*3
      ! for later limiting of rh
      zrhmax  (i) = 0.0_ireals
    ENDDO ! i

    CALL tautsp2D(zpexp_vec, zxexp_vec, kzint_vec, ie_in, 1,ie_in, ke_in+1, &
                  zgamma, zs_vec, zbreak_vec, zcoef_vec, izln_vec, izerror)

    IF (izerror == 0) THEN

      DO k = 1, kevert
        DO i = 1, ie_in
          ! Pressure on the full levels and at the right position
          zpmain(i,k) = akmain(k) + bkmain(k) * ps_co(i,j)
        ENDDO
      ENDDO

      ! check pressure at uppermost level (may occure in case of sigma coordinates )
      DO i = 1, ie_in
        IF (zpmain(i,1) < zpexp_vec(i,1)) THEN
          zpmain(i,1) = zpexp_vec(i,1) + 1.0_ireals
        ENDIF
      ENDDO

      izn(:) = 1
      DO k = 1, kevert
        ldone(:) = .FALSE.
        idone    = 0
        DO WHILE (idone < ie_in)
          DO i = 1, ie_in
            IF (.NOT. ldone(i)) THEN
              zprod = (zpmain(i,k) - zbreak_vec(i,izn(i)))*(zpmain(i,k) - zbreak_vec(i,izn(i)+1))
              IF (zprod <= 0.0) THEN
                izindex(i,k) = izn(i)
                ldone(i)     = .TRUE.
                idone        = idone + 1
              ELSE
                izn(i) = izn(i) + 1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO k = 1, kevert
        DO i = 1, ie_in
          zdpx = zpmain(i,k) - zbreak_vec(i,izindex(i,k))
          xin(i,j,k) =  zcoef_vec(i,1,izindex(i,k))+zdpx  * (zcoef_vec(i,2,izindex(i,k))+  &
              zdpx*0.5*(zcoef_vec(i,3,izindex(i,k))+zdpx/3.0*zcoef_vec(i,4,izindex(i,k))))
        ENDDO
      ENDDO

      IF ( yfld == 'rh' ) THEN
        DO k=1, kevert + 1   ! but limit to kzint_vec(i) later
          DO i = 1, ie_in
            ! Limit generalized relative humidity to 0.001 .. rhmax
            IF (k <= kzint_vec(i)) THEN
              zrhmax(i) = MAX ( zrhmax(i), zxexp_vec(i,k) )
            ENDIF
          ENDDO
        ENDDO

        DO k = 1, kevert
          DO i = 1, ie_in
            xin(i,j,k) = MIN(zrhmax(i),    xin(i,j,k))
            xin(i,j,k) = MAX(0.001_ireals, xin(i,j,k))
          ENDDO
        ENDDO
      ENDIF

    ELSE
      PRINT *, '*** ERROR in tautsp2D: while processing j-index:  ', j
      ierr_intv = 1
      yerr_intv = 'Error in tautsp'
    ENDIF
  ENDDO ! j

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE vert_interpol_p2h

!==============================================================================

END MODULE src_pressure_to_hybrid
