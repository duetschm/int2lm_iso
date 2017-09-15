!+ Source Module for computing LM fields
!==============================================================================

MODULE src_lm_fields

!==============================================================================
!
! Description:
!   This module contains routines for computing LM-fields:
!     - vertical velocity  w_lm
!     - interpolation to LM levels
!     - pressure deviation
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
! V1_5         2007/07/09 Ulrich Schaettler
!  Removed a DEALLOCATE statement and unnecessary variables
! V1_6         2007/09/07 Ulrich Schaettler
!  Editorial changes
! V1_8         2008/05/29 Ulrich Schaettler, Hans-Juergen Panitz
!  Vertical interpolation of qr, qs, qg (if present) to COSMO-Model levels
!  Adaptations, if the COSMO-Model atmosphere is higher than the input model
!  Introduction of debug output
! V1_9         2009/09/03 Ulrich Schaettler, et al.
!  Clip values of qi, qr, qs, qg after vertical interpolation, if they are
!  too small
!  Call to new SR moist_split (only commented: this does not vectorize)
!  Use qi in virtual temperature (MCH) (lmixcld)
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_12        2010/06/14 Oliver Fuhrer
!  Modifications for vectorization of routine moist_split
! V1_14        2010/11/19 Ulrich Schaettler
!  Corrected vertical dimension of qg in call to SR vert_int_lm
!  Treatment of JMA data: usage of t_2m and qv_2m for lower boundary condition
! V1_17        2011/03/11 Ulrich Schaettler
!  Prescribe the bottom boundary condition for vertical interpolation of qr, qs.
!  This was forgotten before and lead to non-reproducible results
! V1_19        2012/06/06 Ulrich Schaettler
!  Adaptations to environment.f90 (here exchg_boundaries), because of
!   unification with COSMO-Model 4.23
!  Correct computation of splines when top height of fine grid model is larger
!    than that of the coarse grid model (initialization of izln_vec)
! V1_22        2013/07/11 Ulrich Schaettler
!  Renamed lprog_qrqs to lprog_qr_qs to be consistent with other names
!  Use structure vcoord_out from vgrid_refatm_utils
! V1_23        2013/10/02 Ulrich Schaettler
!  Rename vcoord_out to vcoord (as is in COSMO-Model)
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

!------------------------------------------------------------------------------

USE data_fields_lm, ONLY : &
  fis_lm    ,      & ! orography * G                                    (m2/s2)
  hsurf_lm  ,      & ! orography                                        (  m  )
  ps_lm     ,      & ! surface pressure                                 ( Pa  )
  t_s_lm    ,      & ! temperature of the ground surface                (  K  )
  t_2m_gl   ,      & ! 2m temperature                                   (  K  )
  u_lm      ,      & ! zonal wind speed                                 ( m/s )
  v_lm      ,      & ! meridional wind speed                            ( m/s )
  w_lm      ,      & ! vertical wind speed (defined on half levels)     ( m/s )
  t_lm      ,      & ! temperature                                      (  K  )
  qv_lm     ,      & ! specific water vapor content                     (kg/kg)
  qc_lm     ,      & ! specific cloud water content                     (kg/kg)
  qi_lm     ,      & ! cloud ice content                                (kg/kg)
  qr_lm     ,      & ! rain      content                                (kg/kg)
  qs_lm     ,      & ! snow      content                                (kg/kg)
  qg_lm     ,      & ! graupel   content                                (kg/kg)
  grh_lm    ,      & ! generalized relative humidity                    (kg/kg)
  pp_lm     ,      & ! deviation from the reference pressure            ( Pa  )
  p0_lm     ,      & ! reference pressure                               ( Pa  )
  t0_lm     ,      & ! reference temperature                            (  K  )
  dp0_lm    ,      & ! reference pressure thickness of layers           ( Pa  )
  hhl_lm    ,      & ! height of half-levels of LM                      (  m  )
! iso code
  riso_lm            ! isotope ratios in water vapor
! end iso code

!------------------------------------------------------------------------------

USE data_grid_lm,   ONLY : &
  dlat,        & ! grid point distance in zonal direction (in degrees)
  dlon,        & ! grid point distance in meridional direction (in degrees)
  startlat,    & ! transformed latitude of the lower left grid point
                 ! of the local domain (in degrees, N>0)
  startlon,    & ! transformed longitude of the lower left grid point
                 ! of the local domain (in degrees, E>0)
  ie2lm,       & !
  je2lm,       & !
  kelm,        & !
  kedim,       & !
  ke1lm,       & !
  jstartpar,   & ! start index for computations in the parallel program
  jendpar        ! end index for computations in the parallel program

!------------------------------------------------------------------------------

USE data_grid_in,   ONLY : &
  ak_in ,      & ! vertical coordinate parameters for half levels
  bk_in ,      & !                  - " -
  akh_in,      & ! vertical coordinate parameters for main levels
  bkh_in,      & !                  - " -
  dak_in,      & ! differences between half levels
  dbk_in,      & !                  - " -
  ke_in ,      & ! ke for input grid
  ke1in          !

!------------------------------------------------------------------------------

USE data_int2lm_io,        ONLY : &
  var_lm         !

!------------------------------------------------------------------------------

USE data_int2lm_constants, ONLY : &
    R_d,     & ! gas constant for dry air                      [J/K*kg]
    Rdv,     & ! = R_d/R_v,
    Rvd_m_o, & ! = R_v/R_d - 1.0,
    O_m_rdv, & ! = 1. - Rdv
    G,       & ! gravity at sea level                          [ms-2]
    r_earth, & ! mean radius of the earth                      [m]
    degrad,  & !
    B1,      & !  a
    B2_w,    & !  b
    B3,      & !  c/b (0 degree Celsius [Kelvin])
    B4_w,    & !  d
    B2_i,    &
    B4_i,    &
    pi

!------------------------------------------------------------------------------

USE data_int2lm_control,       ONLY :  &
    lgsm2lm,      & ! if .TRUE., GSM => LM
    lmixcld,      & ! if .TRUE.,qi added in grh instead of being directly interp.
    lcomp_bound,  & ! compute fields for boundaries
    lvertwind_ini,& ! if .TRUE., compute vertical wind for LM for initial data
    lvertwind_bd, & ! if .TRUE., compute vertical wind for LM for boundary data
    lprog_qi,     & ! if .TRUE., interpolate qi from GME to LM grid
    lprog_qr_qs,  & ! if .TRUE., interpolate qr,qs to LM grid
    lprog_qg,     & ! if .TRUE., interpolate qg to LM grid
    idbg_level,   & ! to control verbosity of output
    lprintdeb_all,& ! whether all PEs print debug output
    qvmin,        & ! minimum value of water vapor (security)
    qcmin,        & ! minimum value of cloud water (security)
    qimin,        & ! minimum value of cloud ice content (security)
! iso code
    liso            ! if .TRUE., include variables for water isotope simulation
! end iso code

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ncomm_type,      & ! type of communication
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    icomm_cart,      & ! communicator for the virtual cartesian topology
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen        ! length of one column of sendbuf

!------------------------------------------------------------------------------

USE utilities,          ONLY :   tautsp2D
USE environment,        ONLY :   model_abort, exchg_boundaries, comm_barrier
USE parallel_utilities, ONLY :   remark, global_values
USE meteo_utilities,    ONLY :   moist_split
USE vgrid_refatm_utils, ONLY :   vcoord

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================
!+ Organizes the computation of LM-fields
!------------------------------------------------------------------------------

SUBROUTINE org_lm_fields

!------------------------------------------------------------------------------
!
! Description:
!   org_lm_fields organizes the computation of LM-fields that are necessary
!   for initial or boundary files for the nonhydrostatic LM. These
!   computations include:
!     - vertical velocity  w_lm
!     - interpolation to LM levels
!     - pressure deviation
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
! Local arrays
REAL (KIND=ireals)         ::  &
  zfi_hl(ie2lm,je2lm,ke_in+1), & ! geopotential on GME half levels
  zfi_fl(ie2lm,je2lm,ke_in  )    ! geopotential on GME full levels

INTEGER  (KIND=iintegers)  ::  &
  izerror, izdebug,            & ! status and error status variable
  i, j, k

REAL    (KIND=ireals)      ::  &
  zpu, zpo, ztv

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling

!
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  izerror   = 0
  yzerrmsg  = '  '
  yzroutine = 'org_lm_fields'

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
! Section 2: Compute the vertical velocity
!------------------------------------------------------------------------------

  ! Compute geopotential on half levels
  zfi_hl(:,:,ke_in+1) = fis_lm(:,:)

  DO k = ke_in, 1, - 1
    IF (k == 1) THEN
      zfi_hl(:,:,k) = zfi_hl(:,:,k+1)                                     &
       + R_d * t_lm(:,:,k) * (1.0 + Rvd_m_o*qv_lm(:,:,k) - qc_lm(:,:,k))  &
       * LOG( 2.0_ireals )
    ELSE
      zfi_hl(:,:,k) = zfi_hl(:,:,k+1)                                     &
       + R_d * t_lm(:,:,k) * (1.0 + Rvd_m_o*qv_lm(:,:,k) - qc_lm(:,:,k))  &
       * LOG((ak_in(k+1)+bk_in(k+1)*ps_lm(:,:))/                          &
             (ak_in(k)  +bk_in(k)  *ps_lm(:,:)))
    ENDIF
  ENDDO

  IF ( (lvertwind_ini .AND. .NOT. lcomp_bound) .OR.         &
       (lvertwind_bd  .AND.       lcomp_bound) ) THEN

    ! Compute the vertical velocity on coarse grid levels
    CALL vertical_velocity (zfi_hl)

    ! Interpolate vertical velocity to LM levels
    IF (vcoord%vcflat > 0.0_ireals) THEN
      CALL vert_int_lm (w_lm, 'w', kedim+1, zfi_hl, ke_in+1, izdebug)
    ELSE
      CALL vert_z_lm   (w_lm, 'w', kedim+1, zfi_hl, ke_in+1)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Vertical interpolations for the other variables
!------------------------------------------------------------------------------

  ! Compute geopotential on full coarse grid levels
  ! lowest full level
  k = ke_in
  DO j = 1, je2lm
    DO i = 1, ie2lm
      zpo = akh_in(k)   + bkh_in(k)   * ps_lm(i,j)
      zpu =  ak_in(k+1) +  bk_in(k+1) * ps_lm(i,j)
      ztv = t_lm(i,j,k) * (1.0 + Rvd_m_o*qv_lm(i,j,k) - qc_lm(i,j,k) )
      zfi_fl(i,j,k) = fis_lm(i,j) + R_d * ztv * LOG(zpu/zpo)
    ENDDO
  ENDDO

  ! all the other levels
  DO k = ke_in-1, 1, -1
    DO j = 1, je2lm
      DO i = 1, ie2lm
        zpo = akh_in(k)   + bkh_in(k)   * ps_lm(i,j)
        zpu =  ak_in(k+1) +  bk_in(k+1) * ps_lm(i,j)
        ztv = t_lm(i,j,k) * (1.0 + Rvd_m_o*qv_lm(i,j,k) - qc_lm(i,j,k) )
        zfi_fl(i,j,k) = zfi_hl (i,j,k+1) + R_d * ztv * LOG(zpu / zpo)
      ENDDO
    ENDDO
  ENDDO

  IF (vcoord%vcflat > 0.0_ireals) THEN
    CALL vert_int_lm (u_lm  , 'u'  , kedim, zfi_fl, ke_in, izdebug)
    CALL vert_int_lm (v_lm  , 'v'  , kedim, zfi_fl, ke_in, izdebug)
    CALL vert_int_lm (t_lm  , 't'  , kedim, zfi_fl, ke_in, izdebug)
    CALL vert_int_lm (grh_lm, 'grh', kedim, zfi_fl, ke_in, izdebug)
! iso code
    ! introduce lower limit for relative humidity
    IF (liso) THEN
      DO k = 6, ke_in
        DO j = 1, je2lm
          DO i = 1, ie2lm
            grh_lm(i,j,k) = MAX(grh_lm(i,j,k), 0.001)
          ENDDO
        ENDDO
      ENDDO

      CALL vert_int_lm (riso_lm(:,:,:,1), 'r18O', kedim, zfi_fl, ke_in,      &
                        izdebug)
      CALL vert_int_lm (riso_lm(:,:,:,2), 'r2H' , kedim, zfi_fl, ke_in,      &
                        izdebug)
    ENDIF
! end iso code
    IF (lprog_qi .AND. (.NOT. lmixcld)) THEN
      CALL vert_int_lm (qi_lm, 'qi', kedim, zfi_fl, ke_in, izdebug)
    ENDIF
    IF (lprog_qr_qs) THEN
      CALL vert_int_lm (qr_lm, 'qr', kedim, zfi_fl, ke_in, izdebug)
      CALL vert_int_lm (qs_lm, 'qs', kedim, zfi_fl, ke_in, izdebug)
    ENDIF
    IF (lprog_qg) THEN
      CALL vert_int_lm (qg_lm, 'qg', kedim, zfi_fl, ke_in, izdebug)
    ENDIF
  ELSE
    CALL vert_z_lm   (u_lm  , 'u'  , kedim, zfi_fl, ke_in)
    CALL vert_z_lm   (v_lm  , 'v'  , kedim, zfi_fl, ke_in)
    CALL vert_z_lm   (t_lm  , 't'  , kedim, zfi_fl, ke_in)
    CALL vert_z_lm   (grh_lm, 'grh', kedim, zfi_fl, ke_in)
    IF (lprog_qi .AND. (.NOT. lmixcld)) THEN
      CALL vert_z_lm   (qi_lm, 'qi', kedim, zfi_fl, ke_in)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Pressure deviation and final splitting of relative humidity
!------------------------------------------------------------------------------

  IF (num_compute > 1) THEN
    CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
  ENDIF

  CALL pressure_deviation

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE org_lm_fields

!==============================================================================
!+ Computation of the vertical velocity.
!------------------------------------------------------------------------------

SUBROUTINE vertical_velocity (fihl_gl)

!------------------------------------------------------------------------------
!
! Description:
!  Computation of the vertical velocity w on the half levels of the GME
!  vertical system. This field is interpolated to the LM vertical system
!  afterwards.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
REAL (KIND=ireals),     INTENT(IN)   ::   &
  fihl_gl (ie2lm,je2lm,ke_in+1)    ! geopotential on GME half levels

! Local scalars:
INTEGER (KIND=iintegers) ::   &
  i, j, k, k2, izerror, kzdims(24)

REAL    (KIND=ireals)    ::   &
  zdlon, zdlat, zlats, zum, zvm, zps, ztm, zfi, zfj, zpse, zpsw, zpsn, zpss

! Local arrays:
REAL    (KIND=ireals)    ::     &
  zclat   (je2lm,2),            & !
  zaclatr (je2lm,2),            & !
  zw1     (ie2lm,je2lm)           !

REAL    (KIND=ireals)    ::     &
  zdiv    (ie2lm,je2lm,ke_in),  & !
  zfeld   (ie2lm,je2lm,ke_in),  & !
  ztvlm   (ie2lm,je2lm,ke_in)     !

CHARACTER (LEN=80)       ::   &
  yzerrmsg

!- End of header
!------------------------------------------------------------------------------

  izerror  = 0
  yzerrmsg = '   '

!------------------------------------------------------------------------------
! Section 1: Computation of local fields
!------------------------------------------------------------------------------

  DO k = 1, ke_in
    ! Compute the virtual temperature
    ztvlm(:,:,k) = t_lm(:,:,k) * ( 1 + Rvd_m_o*qv_lm(:,:,k) - qc_lm(:,:,k) )
    zfeld(:,:,k) = 0.0
  ENDDO

  ! Computation of zclat and zaclatr
  ! Beware that gridpoint (1,1) is at the position (startlat,startlon)
  ! only for the local domains which are considered here
  DO j = 1, je2lm
    zlats        = startlat +  (j-1)*dlat
    zclat(j,1)   = COS(  zlats             * degrad )
    zclat(j,2)   = COS( (zlats + 0.5*dlat) * degrad )
    zaclatr(j,1) = 1.0 / (r_earth * zclat(j,1))
    zaclatr(j,2) = 1.0 / (r_earth * zclat(j,2))
  END DO

  zdlon = dlon * degrad
  zdlat = dlat * degrad

!------------------------------------------------------------------------------
! Section 2: Computation of ETA^star (Equation (2.20) of "2. ARBEITSPAPIER")
!------------------------------------------------------------------------------

  zw1(:,:) = 0.0
  DO k = 1, ke_in
    DO j = 2, je2lm-1
      DO i = 2, ie2lm-1
        ! First compute the horizontal divergence in zdiv
        zpse = dak_in(k) + dbk_in(k) * 0.5 * (ps_lm(i,  j  )+ps_lm(i+1,j  ))
        zpsw = dak_in(k) + dbk_in(k) * 0.5 * (ps_lm(i-1,j  )+ps_lm(i  ,j  ))
        zpsn = dak_in(k) + dbk_in(k) * 0.5 * (ps_lm(i  ,j  )+ps_lm(i  ,j+1))
        zpss = dak_in(k) + dbk_in(k) * 0.5 * (ps_lm(i  ,j-1)+ps_lm(i  ,j  ))

        zdiv(i,j,k) = zaclatr(j,1) *                                         &
            ( (u_lm(i,j,k) * zpse - u_lm(i-1,j,k) * zpsw) / zdlon            &
            + (v_lm(i,j  ,k) * zclat(j  ,2) * zpsn                           &
              -v_lm(i,j-1,k) * zclat(j-1,2) * zpss)       / zdlat )

        ! 1. term: surface pressure tendency: kept in zw1 (equation 2.22)
        zw1(i,j) = zw1(i,j) + zdiv(i,j,k)

        ! 2. term: integral over the divergence: kept in zfeld
        DO k2 = 1, k-1
          zfeld(i,j,k) = zfeld(i,j,k) + zdiv(i,j,k2)
        ENDDO
      ENDDO
    ENDDO
  ENDDO


  ! put ETA^star to zdiv
  DO k = 1, ke_in
    DO j = 2, je2lm-1
      DO i = 2, ie2lm-1
        zdiv(i,j,k) = bk_in(k) * zw1(i,j) - zfeld(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 3: Compute vertical velocity (equation (2.23) of "2. ARBEITSPAPIER"
!------------------------------------------------------------------------------

  ! Put w to 0.0 at top (k=1) and bottom (k=ke+1) of the atmosphere
  w_lm(:,:,1)       = 0.0
  w_lm(:,:,ke_in+1) = 0.0

  ! Compute the vertical velocity
  DO k = 2, ke_in
    DO j = 2, je2lm-1
      DO i = 2, ie2lm-1
        zum = 0.125 * (u_lm(i-1,j,k-1)+u_lm(i,j,k-1)+u_lm(i-1,j,k)+u_lm(i,j,k))
        zvm = 0.125 * (v_lm(i,j-1,k-1)+v_lm(i,j,k-1)+v_lm(i,j-1,k)+v_lm(i,j,k))
        zfi = fihl_gl(i+1,j,k)-fihl_gl(i-1,j,k)
        zfj = fihl_gl(i,j+1,k)-fihl_gl(i,j-1,k)
        ztm = 0.5*(ztvlm(i,j,k-1)+ztvlm(i,j,k))
        zps = ak_in(k) + bk_in(k)*ps_lm(i,j)

        w_lm(i,j,k) = (zaclatr(j,1) * zfi * zum / zdlon                  &
                                    + zfj * zvm / (r_earth * zdlat)      &
                     - (zdiv(i,j,k) * R_d * ztm) / zps  ) / G
      ENDDO
    ENDDO

    ! Put values at the outer boundaries with the values from the next
    ! inner gridpoint
    IF (my_cart_neigh(1) == -1) THEN
      ! no neighbor to the west
      w_lm(1,2:je2lm-1,k) = w_lm(2,2:je2lm-1,k)
    ENDIF

    IF (my_cart_neigh(3) == -1) THEN
      ! no neighbor to the east
      w_lm(ie2lm,2:je2lm-1,k) = w_lm(ie2lm-1,2:je2lm-1,k)
    ENDIF

    IF (my_cart_neigh(2) == -1) THEN
      ! no neighbor to the north
      w_lm(:,je2lm,k) = w_lm(:,je2lm-1,k)
    ENDIF

    IF (my_cart_neigh(4) == -1) THEN
      ! no neighbor to the south
      w_lm(:,1,k) = w_lm(:,2,k)
    ENDIF
  ENDDO

  IF (num_compute > 1) THEN
    ! Exchange inner boundaries for w_lm
    kzdims=(/ke1in,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                              &
       (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
        ie2lm, je2lm, kzdims, jstartpar, jendpar, 1, nboundlines,      &
        my_cart_neigh, .FALSE., .FALSE., .FALSE.,                      &
        121, ldatatypes, ncomm_type, izerror, yzerrmsg,                &
        w_lm)
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE vertical_velocity

!===============================================================================
!+ Vertical Interpolation to LM levels.
!------------------------------------------------------------------------------

SUBROUTINE vert_int_lm (xlm, yname, idim_lm, fi_gl, idim_gme, idebug)

!------------------------------------------------------------------------------
!
! Description:
!   Does the vertical interpolation on the sigma reference levels of the
!   variables: u_lm, v_lm, w_lm, t_lm and grh_lm
!   (grh = generalized relative humidity, split to qv and ql in subroutine
!   *pressure_deviation*)
!
! Method:
!   The geopotential of the GME and the LM at the full or half levels
!   are put on the abscissa (x-axis: variable pexp).
!   The variables to be interpolated are put one after each other
!   on the ordinate (y-axis: variable xexp) in the array xlm.
!   For each gridpoint the tension spline routine tautsp is called for
!   the interpolation.
!   For break(i) <= x <= break(i+1) the interpolation function has the form
!      F(X) = COEF(1,I)+DX(COEF(2,I)+DX/2(COEF(3,I)+DX/3(COEF(4,I)))
!      using DX=X-BREAK(I) and i=1,...,izl
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers),   INTENT(IN)   ::  &
  idim_lm,      & ! vertical dimension of LM fields
  idim_gme,     & ! vertical dimension of GME fields
  idebug          ! for debug output

REAL (KIND=ireals),         INTENT(IN)   ::  &
  fi_gl (ie2lm,je2lm,idim_gme)  ! geopotential on GME levels used as abscissas

REAL (KIND=ireals),         INTENT(INOUT)::  &
  xlm (ie2lm,je2lm,idim_lm)     ! field to be interpolated

CHARACTER (LEN=*),          INTENT(IN)   ::  &
  yname           ! name of the variable

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers) ::   &
  i, j, k, n, izerror,kstart, idone, &
  izln_vec(ie2lm),            & ! number of abscissas for spline
  nztau_vec(ie2lm),           & ! number of abscissas
  izind_vec(ie2lm,kedim+1),   & !
  izn      (ie2lm)

LOGICAL                        ::     &
  ldone(ie2lm)

REAL    (KIND=ireals)    ::   &
  zdx, zrfmax(ie2lm),         &
  zgamma                        ! tension-parameter = 5.5

! Local arrays:
REAL (KIND=ireals)       ::   &
  zhhl       (ie2lm,ke1lm),       & !
  zpexp_vec  (ie2lm,kedim+2),     & !
  zxexp_vec  (ie2lm,kedim+2),     & !
  zbreak_vec (ie2lm,3*(kedim+5)), & ! abscissas of spline
  zfilmk_vec(ie2lm,kedim+1)

REAL (KIND=ireals)       ::   &
  zcoef_vec (ie2lm,4,3*(kedim+5)),  & ! coefficients  of spline
  zs_vec    (ie2lm,kedim+2,6)         ! work array for tautsp

CHARACTER (LEN=80)       ::   &
  yzerrmsg        ! for error message

CHARACTER (LEN=20)       ::   &
  yzroutine       ! name of the routine

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  izerror   = 0
  yzerrmsg  = '   '
  yzroutine = 'vert_int_lm'
  zgamma    = 5.5
  IF (yname /= 'w') THEN
    nztau_vec(:) = ke_in+2  ! kedim+2 !ke_in + 2      ! u, v, t, grh, pp
  ELSE
    nztau_vec(:) = ke_in+1  ! kedim+1 !ke_in + 1      ! w
  ENDIF

  IF (idebug > 15) THEN
    PRINT *, 'vertical interpolation of:  ', yname, ke1lm, kedim
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Loop over all grid points
!------------------------------------------------------------------------------

  ! Lower boundary values
  DO j = 1, je2lm
    ! For u and v-gridpoints, the geopotential is averaged to the C-grid.
    ! On the eastern and northern boundary the values from the last mass grid
    ! point are used. For subdomains with a right or upper neighbor, this does
    ! not matter, for the rightmost subdomain, it is the only value available
    IF     (yname == 'u') THEN
      DO i = 1, ie2lm-1
        zpexp_vec(i,1) = 0.5_ireals * (fis_lm(i,j) + fis_lm(i+1,j))
      ENDDO
      zpexp_vec(ie2lm,1) = fis_lm(ie2lm,j)
    ELSEIF (yname == 'v') THEN
      IF (j < je2lm) THEN
        DO i = 1, ie2lm
          zpexp_vec(i,1) = 0.5_ireals * (fis_lm(i,j) + fis_lm(i,j+1))
        ENDDO
      ELSE
        DO i = 1, ie2lm
          zpexp_vec(i,1) = fis_lm(i,je2lm)
        ENDDO
      ENDIF
    ELSEIF (yname == 't' .OR. yname == 'grh') THEN
      IF (lgsm2lm) THEN
        ! t_2m and qv_2m given
        DO i = 1, ie2lm
          zpexp_vec(i,1) = fis_lm(i,j) + 2.0_ireals * g
        ENDDO
      ELSE
        DO i = 1, ie2lm
          zpexp_vec(i,1) = fis_lm(i,j)
        ENDDO
      ENDIF
    ELSE   ! qi, qr, qs
      DO i = 1, ie2lm
        zpexp_vec(i,1) = fis_lm(i,j)
      ENDDO
    ENDIF

    IF (yname == 't') THEN
      IF (lgsm2lm) THEN
        ! t_2m and qv_2m given
        DO i = 1, ie2lm
           zxexp_vec(i,1) = t_2m_gl(i,j)
        ENDDO
      ELSE
        DO i = 1, ie2lm
           zxexp_vec(i,1) = t_s_lm(i,j)
        ENDDO
      ENDIF
    ELSE IF (yname == 'u') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = 0.0
      ENDDO
    ELSE IF (yname == 'v') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = 0.0
      ENDDO
    ELSE IF (yname == 'grh') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = xlm(i,j,ke_in)
      ENDDO
    ELSE IF (yname=='qi' .OR. yname=='qr' .OR. yname=='qs' .OR. yname=='qg') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = xlm(i,j,ke_in)
      ENDDO
    ELSE IF (yname == 'w') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = xlm(i,j,ke_in+1)
      ENDDO
! iso code
    ELSE IF ((yname == 'r18O') .OR. (yname == 'r2H')) THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = xlm(i,j,ke_in)
      ENDDO
! end iso code
    ELSE
      ! set error code, if no bottom boundary condition has been set
      izerror  = 10
      yzerrmsg = 'Missing bottom boundary for '//yname
      CALL remark (my_cart_id, yzroutine, yzerrmsg)
      CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
    END IF

    DO k = 2, ke_in+1
      ! For u and v-gridpoints, the geopotential is averaged to the C-grid.
      ! On the eastern and northern boundary the values from the last mass grid
      ! point are used. For subdomains with a right or upper neighbor, this does
      ! not matter, for the rightmost subdomain, it is the only value available
      IF     (yname == 'u') THEN
        DO i = 1, ie2lm-1
          zpexp_vec(i,k) = 0.5_ireals *                                 &
                           (fi_gl(i,j,ke_in+2-k) + fi_gl(i+1,j,ke_in+2-k))
        ENDDO
        zpexp_vec(ie2lm,k) = fi_gl(ie2lm,j,ke_in+2-k)
      ELSEIF (yname == 'v') THEN
        IF (j < je2lm) THEN
          DO i = 1, ie2lm
            zpexp_vec(i,k) = 0.5_ireals *                               &
                           (fi_gl(i,j,ke_in+2-k) + fi_gl(i,j+1,ke_in+2-k))
          ENDDO
        ELSE
          DO i = 1, ie2lm
            zpexp_vec(i,k) = fi_gl(i,je2lm,ke_in + 2 - k)
          ENDDO
        ENDIF
      ELSE
        DO i = 1, ie2lm
          zpexp_vec(i,k) = fi_gl(i,j,ke_in + 2 - k)
        ENDDO
      ENDIF

      DO i = 1, ie2lm
        zxexp_vec(i,k) =   xlm(i,j,ke_in + 2 - k)
      ENDDO
    ENDDO

    ! Upper boundary values (set them in any case)
    DO i = 1, ie2lm
      zpexp_vec(i,ke_in+2) = zpexp_vec(i,ke_in+1) + 1000.0
      zxexp_vec(i,ke_in+2) = xlm(i,j,1)
    ENDDO

    ! Set the values for hhl at the correct grid points
    IF     (yname == 'w') THEN
      DO k = 1, kelm
        DO i = 1, ie2lm
          zhhl(i,k) = hhl_lm(i,j,k)
        ENDDO
      ENDDO
    ELSEIF (yname == 'u') THEN
      DO k = 1, kelm
        DO i = 1, ie2lm-1
          zhhl(i,k) = 0.25 * (hhl_lm(i,j,k  ) + hhl_lm(i+1,j,k  )      &
                            + hhl_lm(i,j,k+1) + hhl_lm(i+1,j,k+1))
        ENDDO
        zhhl(ie2lm,k) = 0.5 * (hhl_lm(ie2lm,j,k) + hhl_lm(ie2lm,j,k+1))
      ENDDO
    ELSEIF (yname == 'v') THEN
      DO k = 1, kelm
        IF (j < je2lm) THEN
          DO i = 1, ie2lm
            zhhl(i,k) = 0.25 * (hhl_lm(i,j,k  ) + hhl_lm(i,j+1,k  )      &
                              + hhl_lm(i,j,k+1) + hhl_lm(i,j+1,k+1))
          ENDDO
        ELSE
          DO i = 1, ie2lm
            zhhl(i,k) = 0.5 * (hhl_lm(i,je2lm,k) + hhl_lm(i,je2lm,k+1))
          ENDDO
        ENDIF
      ENDDO
    ELSE
      DO k = 1, kelm
        DO i = 1, ie2lm
          zhhl(i,k) = 0.5 * (hhl_lm(i,j,k)+hhl_lm(i,j,k+1))
        ENDDO
      ENDDO
    ENDIF

    !  If the atmosphere of the COSMO-Model (hhl) is higher than the atmosphere
    !  of the input model, the profile in the vertical interpolation is held
    !  constant up to the COSMO-model height: this is done by pushing the highest
    !  input level to the COSMO-model height
    zpexp_vec(:,ke_in+2) = MAX (zpexp_vec(:,ke_in+2), zhhl(:,1)*g)

    ! Now do the interpolation
    izln_vec(:) = (kedim+5) * 3
    CALL tautsp2D(zpexp_vec, zxexp_vec, nztau_vec, ie2lm, 1, ie2lm, kedim+2, &
                  zgamma, zs_vec, zbreak_vec, zcoef_vec, izln_vec, izerror)

    IF (izerror == 0) THEN

      IF (yname == 'w') THEN
        kstart = 2
      ELSE
        kstart = 1
      ENDIF

      ! Bug fix in 1.19: before it was izn = izln_vec
      ! but we really have to start in the interval below
      izn(:) = izln_vec(:) - 1
      DO k = kstart, kelm
        ldone(:) = .FALSE.
        idone    = 0
        DO WHILE (idone < ie2lm)
          DO i = 1, ie2lm
            IF (.NOT. ldone(i)) THEN
              zfilmk_vec(i,k) = zhhl(i,k) * G
              IF (zbreak_vec(i,izn(i)) <= zfilmk_vec(i,k) .AND. &
                  zfilmk_vec(i,k) <= zbreak_vec(i,izn(i)+1)) THEN
                izind_vec(i,k) = izn(i)
                ldone(i)       = .TRUE.
                idone          = idone + 1
              ELSE
                izn(i) = izn(i) - 1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO k = kstart, kelm
        DO i = 1, ie2lm
          zdx        = zfilmk_vec(i,k) - zbreak_vec(i,izind_vec(i,k))
          xlm(i,j,k) = zcoef_vec(i,1,izind_vec(i,k)) +          &
                       zdx*(zcoef_vec(i,2,izind_vec(i,k)) +     &
                       zdx*0.5*(zcoef_vec(i,3,izind_vec(i,k)) + &
                       zdx/3.0*zcoef_vec(i,4,izind_vec(i,k))))
        ENDDO
      ENDDO

      ! Limit generalized relative humidity to 0.001 .. rhmax
      IF (yname == 'grh') THEN
        zrfmax(:) = 0.0_ireals
        DO k = 1, kedim+2    !nztau_vec(i)
          DO i = 1, ie2lm
            IF (k <= nztau_vec(i)) THEN
              zrfmax(i) = MAX ( zrfmax(i) , zxexp_vec(i,k) )
            ENDIF
          ENDDO
        ENDDO

        DO k = 1, kelm
          DO i = 1, ie2lm
            xlm(i,j,k) = MAX ( 0.001_ireals, xlm(i,j,k) )
            xlm(i,j,k) = MIN ( xlm(i,j,k), zrfmax(i) )
          ENDDO
        ENDDO
      ENDIF

      ! Limit values of qi
      IF (yname=='qi' .OR. yname=='qr' .OR. yname=='qs' .OR. yname=='qg') THEN
        DO k = 1, kelm
          DO i = 1, ie2lm
            IF (xlm(i,j,k) < qimin) THEN
              xlm(i,j,k) = 0.0_ireals
            ENDIF
          ENDDO
        ENDDO
      ENDIF

! iso code
      ! limit isotope ratios
      IF ((yname == 'r18O') .OR. (yname == 'r2H')) THEN
        DO k = 1, kelm
          DO i = 1, ie2lm
            IF (xlm(i,j,k) < 0.0_ireals) THEN
              xlm(i,j,k) = 0.0_ireals
            ENDIF
            IF (xlm(i,j,k) > 1.5_ireals) THEN
              xlm(i,j,k) = 1.5_ireals
            ENDIF
          ENDDO
        ENDDO 
      ENDIF
! end iso code

      ! Set w_lm to 0.0 at the top and the bottom
      IF (yname == 'w') THEN
        DO i = 1, ie2lm
          xlm(i,j,    1) = 0.0_ireals
          xlm(i,j,ke1lm) = 0.0_ireals
        ENDDO
      ENDIF

    ELSE
      yzerrmsg = 'Error in tautsp'
      PRINT *, '*** ERROR in tautsp2D: while processing j-index:  ', j
      CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

  ENDDO ! j = 1, je2lm

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE vert_int_lm

!===============================================================================
!+ Vertical Interpolation to LM levels for the z-coordinate.
!------------------------------------------------------------------------------

SUBROUTINE vert_z_lm  (xlm, yname, idim_lm, fi_gl, idim_gme)

!------------------------------------------------------------------------------
!
! Description:
!   Does the vertical interpolation to the z-coordinate levels of the
!   variables: u_lm, v_lm, w_lm, t_lm and grh_lm
!   (grh = generalized relative humidity, split to qv and ql in *presdev*)
!
! Method:
!   The geopotential of the GME and the LM at the full or half levels
!   are put on the abscissa (x-axis: variable pexp). The lowest level is put
!   50 meters below the lowest LM z-level. The surface value of the
!   corresponding variable is used.
!   The variables to be interpolated are put one after each other
!   on the ordinate (y-axis: variable xexp) in the array xlm.
!   For each gridpoint the tension spline routine tautsp is called for
!   the interpolation.
!   For break(i) <= x <= break(i+1) the interpolation function has the form
!      F(X) = COEF(1,I)+DX(COEF(2,I)+DX/2(COEF(3,I)+DX/3(COEF(4,I)))
!      using DX=X-BREAK(I) and i=1,...,izl
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers),   INTENT(IN)   ::  &
  idim_lm,      & ! vertical dimension of LM fields
  idim_gme        ! vertical dimension of GME fields

REAL (KIND=ireals),         INTENT(IN)   ::  &
  fi_gl (ie2lm,je2lm,idim_gme)  ! geopotential on GME levels used as abscissas

REAL (KIND=ireals),         INTENT(INOUT)::  &
  xlm (ie2lm,je2lm,idim_lm)     ! field to be interpolated

CHARACTER (LEN=*),          INTENT(IN)   ::  &
  yname           ! name of the variable

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers) ::     &
  i, j, k, n, kstart, izerror,  &
  izln_vec(ie2lm),              & ! number of abscissas for spline
  nztau_vec(ie2lm),             & ! number of abscissas
  izind_vec(ie2lm,kedim+2)

REAL    (KIND=ireals)    ::   &
  zdx, zrfmax,                &
  zminfis, zgamma

! Local arrays:
REAL (KIND=ireals)       ::   &
  zhhl       (ie2lm,ke1lm),       & !
  zpexp_vec  (ie2lm,kedim+4),     & !
  zxexp_vec  (ie2lm,kedim+4),     & !
  zbreak_vec (ie2lm,3*(kedim+7)), & ! abscissas of spline
  zfilmk_vec(ie2lm,kedim+3)

REAL (KIND=ireals)       ::   &
  zcoef_vec (ie2lm,4,3*(kedim+7)),  & ! coefficients  of spline
  zs_vec    (ie2lm,kedim+4,6)         ! work array for tautsp

CHARACTER (LEN=80)       ::   &
  yzerrmsg        ! for error message

CHARACTER (LEN=20)       ::   &
  yzroutine       ! name of the routine

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  izerror   = 0
  yzerrmsg  = '   '
  yzroutine = 'vert_z_lm'
  zgamma    = 5.5
  izln_vec(:) = (kedim+7) * 3
  IF (yname /= 'w') THEN
    nztau_vec(:) = ke_in + 4      ! u, v, t, grh
  ELSE
    nztau_vec(:) = ke_in + 3      ! w
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Loop over all grid points
!------------------------------------------------------------------------------

  ! Determine the minimum of fis_lm:
  zminfis = MINVAL (fis_lm(:,:))

  IF (num_compute > 1) THEN
    CALL global_values (zminfis, 1, 'MIN', imp_reals, icomm_cart, -1,   &
                        yzerrmsg, izerror)
  ENDIF

  ! Lower boundary values
  DO j = 1, je2lm
    ! Set the x-axis for the lowest LM z-level:
    zpexp_vec(:,1) = (zminfis - 100.0_ireals) * g
    zpexp_vec(:,2) = (zminfis -  20.0_ireals) * g

    ! The next value that is set is the surface value
    ! For u and v-gridpoints, the geopotential is averaged to the C-grid.
    ! On the eastern and northern boundary the values from the last mass grid
    ! point are used. For subdomains with a right or upper neighbor, this does
    ! not matter, for the rightmost subdomain, it is the only value available
    IF     (yname == 'u') THEN
      DO i = 1, ie2lm-1
        zpexp_vec(i,3) = 0.5_ireals * (fis_lm(i,j) + fis_lm(i+1,j))
      ENDDO
      zpexp_vec(ie2lm,3) = fis_lm(ie2lm,j)
    ELSEIF (yname == 'v') THEN
      IF (j < je2lm) THEN
        DO i = 1, ie2lm
          zpexp_vec(i,3) = 0.5_ireals * (fis_lm(i,j) + fis_lm(i,j+1))
        ENDDO
      ELSE
        DO i = 1, ie2lm
          zpexp_vec(i,3) = fis_lm(i,je2lm)
        ENDDO
      ENDIF
    ELSE
      DO i = 1, ie2lm
        zpexp_vec(i,3) = fis_lm(i,j)
      ENDDO
    ENDIF

    IF (yname == 't') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = t_s_lm(i,j)
         zxexp_vec(i,2) = t_s_lm(i,j)
         zxexp_vec(i,3) = t_s_lm(i,j)
      ENDDO
    ELSE IF (yname == 'u') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = 0.0
         zxexp_vec(i,2) = 0.0
         zxexp_vec(i,3) = 0.0
      ENDDO
    ELSE IF (yname == 'v') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = 0.0
         zxexp_vec(i,2) = 0.0
         zxexp_vec(i,3) = 0.0
      ENDDO
    ELSE IF (yname == 'grh') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = xlm(i,j,ke_in)
         zxexp_vec(i,2) = xlm(i,j,ke_in)
         zxexp_vec(i,3) = xlm(i,j,ke_in)
      ENDDO
    ELSE IF (yname == 'qi') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = xlm(i,j,ke_in)
         zxexp_vec(i,2) = xlm(i,j,ke_in)
         zxexp_vec(i,3) = xlm(i,j,ke_in)
      ENDDO
    ELSE IF (yname == 'w') THEN
      DO i = 1, ie2lm
         zxexp_vec(i,1) = xlm(i,j,ke_in+1)
         zxexp_vec(i,2) = xlm(i,j,ke_in+1)
         zxexp_vec(i,3) = xlm(i,j,ke_in+1)
      ENDDO
    END IF

    DO k = 2, ke_in+1
      ! For u and v-gridpoints, the geopotential is averaged to the C-grid.
      ! On the eastern and northern boundary the values from the last mass grid
      ! point are used. For subdomains with a right or upper neighbor, this does
      ! not matter, for the rightmost subdomain, it is the only value available
      IF     (yname == 'u') THEN
        DO i = 1, ie2lm-1
          zpexp_vec(i,k+2) = 0.5_ireals *                                 &
                           (fi_gl(i,j,ke_in+2-k) + fi_gl(i+1,j,ke_in+2-k))
        ENDDO
        zpexp_vec(ie2lm,k+2) = fi_gl(ie2lm,j,ke_in+2-k)
      ELSEIF (yname == 'v') THEN
        IF (j < je2lm) THEN
          DO i = 1, ie2lm
            zpexp_vec(i,k+2) = 0.5_ireals *                               &
                           (fi_gl(i,j,ke_in+2-k) + fi_gl(i,j+1,ke_in+2-k))
          ENDDO
        ELSE
          DO i = 1, ie2lm
            zpexp_vec(i,k+2) = fi_gl(i,je2lm,ke_in + 2 - k)
          ENDDO
        ENDIF
      ELSE
        DO i = 1, ie2lm
          zpexp_vec(i,k+2) = fi_gl(i,j,ke_in + 2 - k)
        ENDDO
      ENDIF

      DO i = 1, ie2lm
        zxexp_vec(i,k+2) =   xlm(i,j,ke_in + 2 - k)
      ENDDO
    ENDDO
    izln_vec(:) = (kedim+7) * 3

    ! Upper boundary values
    IF (yname /= 'w') THEN
      DO i = 1, ie2lm
        zpexp_vec(i,ke_in+4) = zpexp_vec(i,ke_in  +3) + 1000.0
        zxexp_vec(i,ke_in+4) = xlm(i,j,1)
      ENDDO
    ENDIF

    ! Set the values for hhl at the correct grid points
    IF     (yname == 'w') THEN
      DO k = 1, kelm
        DO i = 1, ie2lm
          zhhl(i,k) = hhl_lm(i,j,k)
        ENDDO
      ENDDO
    ELSE
      DO k = 1, kelm
        DO i = 1, ie2lm
          zhhl(i,k) = 0.5 * (hhl_lm(i,j,k)+hhl_lm(i,j,k+1))
        ENDDO
      ENDDO
    ENDIF

    ! Now do the interpolation
    CALL tautsp2D(zpexp_vec, zxexp_vec, nztau_vec, ie2lm, 1, ie2lm, kedim+2, &
                  zgamma, zs_vec, zbreak_vec, zcoef_vec, izln_vec, izerror)

    IF (izerror == 0) THEN

      IF (yname == 'w') THEN
        kstart = 2
      ELSE
        kstart = 1
      ENDIF

        DO k = kstart, kelm ! k-loop has to start only from k=2
          DO n =  1, MAXVAL (izln_vec)    ! kedim*2
            DO i = 1, ie2lm
              IF ( n <= izln_vec(i)-1) THEN
                zfilmk_vec(i,k) = zhhl(i,k) * G
                IF (zbreak_vec(i,n) <= zfilmk_vec(i,k) .AND. &
                  zfilmk_vec(i,k) <= zbreak_vec(i,n+1)) THEN
                  izind_vec(i,k) = n
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO

        DO k = kstart, kelm
          DO i = 1, ie2lm
            zdx        = zfilmk_vec(i,k) - zbreak_vec(i,izind_vec(i,k))
            xlm(i,j,k) = zcoef_vec(i,1,izind_vec(i,k)) +          &
                         zdx*(zcoef_vec(i,2,izind_vec(i,k)) +     &
                         zdx*0.5*(zcoef_vec(i,3,izind_vec(i,k)) + &
                         zdx/3.0*zcoef_vec(i,4,izind_vec(i,k))))
          ENDDO
        ENDDO

      DO i = 1, ie2lm

        ! Limit generalized relative humidity to 0.001 .. rhmax
        IF (yname == 'grh') THEN
          zrfmax = 0.0_ireals
          DO k = 1, nztau_vec(i)
            zrfmax = MAX ( zrfmax , zxexp_vec(i,k) )
          END DO
          xlm(i,j,1:kelm) = MAX ( 0.001_ireals, xlm(i,j,1:kelm) )
          xlm(i,j,1:kelm) = MIN ( xlm(i,j,1:kelm), zrfmax )
        ENDIF

        ! Limit values of qi
        IF (yname == 'qi') THEN
          DO k = 1, kelm
            IF (xlm(i,j,k) < qimin) THEN
              xlm(i,j,k) = 0.0_ireals
            ENDIF
          ENDDO
        ENDIF

        ! Set w_lm to 0.0 at the top and the bottom
        IF (yname == 'w') THEN
          xlm(i,j,    1) = 0.0_ireals
          xlm(i,j,ke1lm) = 0.0_ireals
        ENDIF
      ENDDO

    ELSE
      yzerrmsg = 'Error in tautsp'
      PRINT *, '*** ERROR in tautsp2D: while processing j-index:  ', j
      CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

  ENDDO ! j = 1, je2lm

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE vert_z_lm

!===============================================================================
!+ Computation of the pressure deviation and splitting of grh.
!------------------------------------------------------------------------------

SUBROUTINE pressure_deviation

!------------------------------------------------------------------------------
!
! Description:
!   Computation of the pressure deviation from the reference pressure on the
!   reference sigma half levels and splitting of the generalized relative
!   humidity into the specific water vapor qv_lm and the specific cloud water
!   content qc_lm.
!
! Method:
!
!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers) ::   &
  i,j,k, kztop, kzlow, izter, k1, k2, k3, k4

REAL    (KIND=ireals)    ::   &
  z1, z2, zaq, zbq, zrhs, zleft, zmlev, dz1, dzh, rh

! Local arrays:
REAL    (KIND=ireals)    ::   &
  zt0dp0t (ie2lm,je2lm,2),    & !
  ztvdt   (ie2lm,je2lm,2),    & !
  zpm     (ie2lm,je2lm),      & !
  p1(kelm)

! Definition of statement functions
REAL    (KIND=ireals)    ::   sf_psat_w, sf_qsat, x, y, z, v, w

sf_psat_w  (x,y,z,v,w) = y * EXP(z*(x-v)/(x-w))
sf_qsat    (x,y,z,v)   = z * x / MAX( (y-v*x), 1.0_ireals)
! GME2LM before: sf_qsat    (x,y,z,v)   = z * x / (y-v*x)

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  pp_lm (:,:,:) = 0.0_ireals

!------------------------------------------------------------------------------
! Section 2: Pressure deviation
!------------------------------------------------------------------------------

  IF (vcoord%vcflat > 0.0_ireals) THEN
    ! Iteration over 3 steps
    DO izter = 1, 3

      ! indices for organization
      kztop = 2
      kzlow = 1

      DO j = 1, je2lm
        DO i = 1, ie2lm
          ! split grh_lm into qv_lm, qc_lm using actual p = pp_lm + p0_lm
          zaq = sf_psat_w (t_lm(i,j,kelm), b1, b2_w, b3, b4_w)
          zpm(i,j) = p0_lm(i,j,kelm)+pp_lm(i,j,kelm)
          zbq = sf_qsat   (zaq, zpm(i,j), Rdv, O_m_rdv)
          qv_lm(i,j,kelm) = MIN(1.0_ireals, grh_lm(i,j,kelm)) * zbq
          qc_lm(i,j,kelm) = MAX(0.0_ireals, grh_lm(i,j,kelm)-1.0_ireals) * zbq
        ENDDO
      ENDDO

      ! qi is treated in grh
      IF (lmixcld) THEN
        CAll moist_split(t_lm(:,:,kelm),zpm,grh_lm(:,:,kelm),qvmin,qcmin, &
                         qimin,pi,b1,b2_w,b2_i,b3,b4_w,b4_i,Rdv,O_m_rdv,  &
                         qv_lm(:,:,kelm),qc_lm(:,:,kelm),qi_lm(:,:,kelm), &
                         ie2lm, je2lm)
      ENDIF

      DO j = 1, je2lm
        DO i = 1, ie2lm
          ! pressure deviation on main level kelm
          z1 = G * 0.5 * (hhl_lm(i,j,kelm) - hsurf_lm(i,j))
          z2 = R_d * (1.0_ireals + Rvd_m_o*qv_lm(i,j,kelm) - qc_lm(i,j,kelm))&
                        * t_lm(i,j,kelm)
          pp_lm(i,j,kelm) = ps_lm(i,j) * EXP (-z1/z2 ) - p0_lm(i,j,kelm)

          ! virtual temperature term
          ! qi is treated in grh
          IF (lmixcld) THEN
            ztvdt(i,j,kzlow) = ((1.+Rvd_m_o*qv_lm(i,j,kelm)-qc_lm(i,j,kelm)    &
                                                           -qi_lm(i,j,kelm)) * &
                          t_lm(i,j,kelm) - t0_lm(i,j,kelm)) / t_lm(i,j,kelm)
          ELSE
            ztvdt(i,j,kzlow) = ((1.+Rvd_m_o*qv_lm(i,j,kelm)-qc_lm(i,j,kelm)) * &
                          t_lm(i,j,kelm) - t0_lm(i,j,kelm)) / t_lm(i,j,kelm)
          ENDIF

          ! coefficient of pressure term
          zt0dp0t(i,j,kzlow) = t0_lm(i,j,kelm) / p0_lm(i,j,kelm)             &
                                               / t_lm(i,j,kelm)
        ENDDO
      ENDDO

      ! Loop over levels
      DO k = kelm, 2, -1
        DO j = 1, je2lm
          DO i = 1, ie2lm
            ! split grh_lm into qv_lm, qc_lm using actual p = pp_lm + p0_lm
            zaq = sf_psat_w (t_lm(i,j,k-1), b1, b2_w, b3, b4_w)
            zpm(i,j) = p0_lm(i,j,k-1) + pp_lm(i,j,k-1)
            zbq = sf_qsat   (zaq, zpm(i,j), Rdv, O_m_rdv)
            qv_lm(i,j,k-1) = MIN(1.0_ireals, grh_lm(i,j,k-1)) * zbq
            qc_lm(i,j,k-1) = MAX(0.0_ireals, grh_lm(i,j,k-1)-1.0_ireals) * zbq
          ENDDO
        ENDDO

        ! virtual temperature term
        ! qi is treated in grh
        IF (lmixcld) THEN
          CAll moist_split(t_lm(:,:,k-1),zpm,grh_lm(:,:,k-1),qvmin,qcmin, &
                           qimin,pi,b1,b2_w,b2_i,b3,b4_w,b4_i,Rdv,O_m_rdv,&
                           qv_lm(:,:,k-1),qc_lm(:,:,k-1),qi_lm(:,:,k-1),  &
                           ie2lm,je2lm)
          DO j = 1, je2lm
            DO i = 1, ie2lm
              ztvdt(i,j,kztop) = ((1.+Rvd_m_o*qv_lm(i,j,k-1)-qc_lm(i,j,k-1)   &
                                                            -qi_lm(i,j,k-1))* &
                           t_lm(i,j,k-1) - t0_lm(i,j,k-1)) / t_lm(i,j,k-1)
            ENDDO
          ENDDO
        ELSE
          DO j = 1, je2lm
            DO i = 1, ie2lm
              ztvdt(i,j,kztop) = ((1.+Rvd_m_o*qv_lm(i,j,k-1)-qc_lm(i,j,k-1)) * &
                           t_lm(i,j,k-1) - t0_lm(i,j,k-1)) / t_lm(i,j,k-1)
            ENDDO
          ENDDO
        ENDIF

        DO j = 1, je2lm
          DO i = 1, ie2lm
            ! coefficient of pressure term
            zt0dp0t(i,j,kztop) = t0_lm(i,j,k-1) / p0_lm(i,j,k-1)             &
                                                / t_lm(i,j,k-1)

            ! pressure deviation on main level k-1
            zleft = 1.0 + 0.5 * dp0_lm(i,j,k  ) * zt0dp0t(i,j,kztop)
            zrhs = (1.0 - 0.5 * dp0_lm(i,j,k-1) * zt0dp0t(i,j,kzlow))        &
                                                               *pp_lm(i,j,k) &
                        + 0.5 * dp0_lm(i,j,k  ) * ztvdt  (i,j,kztop)         &
                        + 0.5 * dp0_lm(i,j,k-1) * ztvdt  (i,j,kzlow)
            pp_lm(i,j,k-1) = zrhs / zleft
          ENDDO
        ENDDO

        ! change organization indices
        kztop = 3 - kztop
        kzlow = 3 - kzlow
      ENDDO

    ENDDO    ! iteration loop

  ELSE

    ! In case of computing fields for the z-coordinate, the following
    ! algorithm is used to calculate the pressure deviation
    ! (by Heinz-Werner Bitzer)
    DO j = 1, je2lm
      DO i = 1, ie2lm
        k1    = 0
        p1(:) = 0.0_ireals

        DO k = 1, ke1lm
          IF (hhl_lm(i,j,k) > hsurf_lm(i,j)) THEN
            k1 = k
          ENDIF
        ENDDO
        k2 = k1 - 1
        k3 = k2 - 1

        zmlev = 0.5 * (hhl_lm(i,j,k1)+hhl_lm(i,j,k1+1))
        dz1   = hsurf_lm(i,j)-zmlev
        rh    = g / (r_d*t_lm(i,j,k2))
        p1(k1)= ps_lm(i,j) * EXP(rh*dz1)

        dzh   = 0.5 * (hhl_lm(i,j,k1-1)-hhl_lm(i,j,k1+1))
        rh    = g / (r_d*t_lm(i,j,k2))
        p1(k2)= p1(k1) * EXP(-rh*dzh)

        DO k = 1, k3
          k4    = k2-k
          rh    = g / (0.5*r_d*(t_lm(i,j,k4+1)+t_lm(i,j,k4)))
          dzh   = 0.5 * (hhl_lm(i,j,k4)-hhl_lm(i,j,k4+2))
          p1(k4)= p1(k4+1) * EXP(-rh*dzh)
        ENDDO

        DO k = 1, kelm
          pp_lm(i,j,k) = p1(k)-p0_lm(i,j,k)
        ENDDO

        DO k = 1, kelm
          IF (k > k1) THEN
            pp_lm(i,j,k) = pp_lm(i,j,k1)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  ENDIF

!------------------------------------------------------------------------------
! 3. Final splitting of relative humidity
!------------------------------------------------------------------------------

  DO k = 1, kelm
    DO j = 1, je2lm
      DO i = 1, ie2lm
        ! split grh_lm into qv_lm, qc_lm using actual p = pp_lm + p0_lm
        zaq = sf_psat_w (t_lm(i,j,k), b1, b2_w, b3, b4_w)
        zpm(i,j) = p0_lm(i,j,k) + pp_lm(i,j,k)
        zbq = sf_qsat   (zaq, zpm(i,j), Rdv, O_m_rdv)
        qv_lm(i,j,k) = MIN(1.0_ireals, grh_lm(i,j,k)) * zbq
        qc_lm(i,j,k) = MAX(0.0_ireals, grh_lm(i,j,k)-1.0_ireals) * zbq
        IF (qv_lm(i,j,k) < qvmin) THEN
          qv_lm(i,j,k)  = qvmin
        ENDIF
        IF (qc_lm(i,j,k) < qcmin) THEN
          qc_lm(i,j,k)  = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
    ! qi is treated in grh
    IF (lmixcld) THEN
      CALL moist_split(t_lm(:,:,k),zpm,grh_lm(:,:,k),qvmin,qcmin,      &
                       qimin,pi,b1,b2_w,b2_i,b3,b4_w,b4_i,Rdv,O_m_rdv, &
                       qv_lm(:,:,k),qc_lm(:,:,k),qi_lm(:,:,k),         &
                       ie2lm, je2lm)
    ENDIF
  ENDDO

! IF (lprgp) THEN
!   DO l = 1, ngp
!     IF (igp(l) == i .AND. jgp(l) == j) THEN
!       pplmpr(k,l) = pp_lm(i,j,k)
!       qvlmpr(k,l) = qv_lm(i,j,k)
!       qllmpr(k,l) = qc_lm(i,j,k)
!     END IF
!   ENDDO ! l
! END IF ! lprgp

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE pressure_deviation

!===============================================================================

END MODULE src_lm_fields
