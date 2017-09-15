!+ Source Module for vertical interpolation.
!==============================================================================

MODULE src_vert_interpol

!==============================================================================
!
! Description:
!   This module contains routines for the vertical interpolation of the
!   atmospheric fields to GME-levels which are adapted to the HM/LM-orography.
!   Also, the routines for adapting the surface pressure to the finer
!   orography are included.
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
! 1.2        2005/07/22 Ulrich Schaettler
!  Bug correction for computing the COS(lat) of v-grid points
! V1_5         2007/07/09 Ulrich Schaettler
!  Added interpolation of qr,qs,qg
!  Eliminated allocation of ps_lm, grh_lm (is done in int2lm_org)
!  Eliminated akhlm, bkhlm
! V1_6         2007/09/07 Ulrich Schaettler, Uwe Boehm, Burkhardt Rockel
!  Added option lcm2lm
! V1_7         2007/11/26 Ulrich Schaettler
!  Introduced debug printouts
! V1_8         2008/05/29 Ulrich Schaettler
!  Vectorization of vertical interpolation and necessary pre- and postprocessing
! V1_9         2009/09/03 Anne Roches, MCH
!  Call to new SR moist_split (only commented: this does not vectorize)
!  Use qi in virtual temperature (MCH) (lmixcld)
! V1_10        2009/12/17 Ulrich Schaettler
!  Repaired vectorization problems in adapt_pressure1/2
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_12        2010/06/14 Oliver Fuhrer
!  Modifications for vectorization of routine moist_split
! V1_14        2010/11/19 Ulrich Schaettler
!  Modifications to process JMA and NCEP data
! V1_16        2010/12/20 Ulrich Schaettler
!  Bug fix in SR adapt_pressure_1: There was an erroneous epsilon when
!    looking for the correct model layer
! V1_19        2012/06/06 Ulrich Schaettler, Hans-Juergen Panitz
!  Added lhir2lm processing
!  Adaptations to environment.f90 (here exchg_boundaries), because of
!   unification with COSMO-Model 4.23
!  Restructured the computations in Section 3 of SR adapt_pressure_1, to get a 
!   clear distinction of cases (and avoid numerical problems as division by zero)
! V1_22        2013/07/11 Ulrich Schaettler
!  Renamed lprog_qrqs to lprog_qr_qs to be consistent with other names
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

USE data_fields_lm, ONLY : &
  fis_lm    ,      & ! orography * G                                    (m2/s2)
  ps_lm     ,      & ! surface pressure                                 ( Pa  )
  fis_gl    ,      & ! GME interpolated orography * G                   (m2/s2)
  fic_gl    ,      & ! check level of orography                         (m2/s2)
  ps_gl     ,      & ! surface pressure on the interpol. GME orogr.     ( Pa  )
  dpsdt_gl  ,      & ! surface pressure tendency                        (Pa/s )
  u_lm      ,      & ! zonal wind speed                                 ( m/s )
  v_lm      ,      & ! meridional wind speed                            ( m/s )
  t_lm      ,      & ! temperature                                      (  K  )
  qv_lm     ,      & ! specific water vapor content                     (kg/kg)
  qc_lm     ,      & ! specific cloud water content                     (kg/kg)
  qi_lm     ,      & ! cloud ice content                                (kg/kg)
  qr_lm     ,      & ! rain      content                                (kg/kg)
  qs_lm     ,      & ! snow      content                                (kg/kg)
  qg_lm     ,      & ! graupel   content                                (kg/kg)
  grh_lm    ,      & ! generalized relative humidity                    (kg/kg)
! iso code
  riso_lm            ! isotope ratios in water vapor
! end iso code

!------------------------------------------------------------------------------

USE data_grid_lm,    ONLY: &
  dlat,        & ! grid point distance in zonal direction (in degrees)
  dlon,        & ! grid point distance in meridional direction (in degrees)
  startlat,    & ! transformed latitude of the lower left grid point
                 ! of the local domain (in degrees, N>0)
  ie2lm_tot,   & ! = ielm_tot + 2
  je2lm_tot,   & ! = jelm_tot + 2
  kelm,        & ! ke for LM
  ie2lm,       & !
  je2lm,       & !
  ke1lm,       & !
  kedim,       & !
  jstartpar,   & ! start index for computations in the parallel program
  jendpar        ! end index for computations in the parallel program

!------------------------------------------------------------------------------

USE data_grid_in,    ONLY: &
  ak_in ,      & ! vertical coordinate parameters for half levels
  bk_in ,      & !                  - " -
  akh_in,      & ! vertical coordinate parameters for main levels
  bkh_in,      & !                  - " -
  ke_in,       & ! ke for input grid
  klv850_in      ! approximate level where 850 hPa is reached

!------------------------------------------------------------------------------

USE data_int2lm_control,    ONLY: &
  lgme2lm,      & ! if .TRUE., gme->lm
  lgsm2lm,      & ! if .TRUE., gsm->lm
  lgfs2lm,      & ! if .TRUE., gfs->lm
  llm2lm,       & ! if .TRUE., lm ->lm
  lec2lm,       & ! if .TRUE., ec ->lm
  lhir2lm,      & ! if .TRUE., hirlam ->lm
  lcm2lm,       & ! if .TRUE., cm ->lm    !_br
  idbg_level,   & ! to control verbosity of output
  lprintdeb_all,& ! whether all PEs print debug output
  luvcor,       & ! if .TRUE., correct winds for given surface pressure
                  ! tendency
  lprog_qi,     & ! if .TRUE., interpolate qi to LM grid
  lmixcld,      & ! if .TRUE., qi added in grh instead of being directly interp.
  lprog_qr_qs,  & ! if .TRUE., interpolate qr,qs to LM grid
  lprog_qg,     & ! if .TRUE., interpolate qg to LM grid
  qvmin,        & ! minimum value of water vapor (security)
  qcmin,        & ! minimum value of cloud water (security)
  qimin,        & ! minimum value of cloud ice content (security)
  noutput,      & ! unit number for output file
  lprps,        & ! logical for print of different ps* und fis - fields
  lprt,         & ! logical for print of T at 2 levels (nlev1pr,nlev2)
  lpru,         & ! same as lprt but for U
  lprv,         & ! same as lprt but for V
  lprgrh,       & ! same as lprt but for General Relative Humidity
  lprqv,        & ! same as lprt but for QV
  lprqc,        & ! same as lprt but for QC
  lprud,        & ! same as lprt but for UD (divergent wind correction)
  lprvd,        & ! same as lprt but for VD (divergent wind correction)
  lprdpdt,      & ! same as lprt but for DPDT (tendency of surface pressure)
  kcontrol_fi,  & ! control level for geopotential
  pcontrol_fi,  & ! pressure of control level for geopotential
  nlev1pr,      & ! k-index for the print of first  model layer
  nlev2pr,      & ! k-index for the print of second model layer
  lbd_frame_cur,& ! if .TRUE., current boundary fields include only frames
! iso code
  liso            ! if .TRUE., include variables for water isotope simulation
! end iso code

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY :  &
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ncomm_type,      & ! type of communication
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    icomm_cart,      & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    sendbuf,         & ! sending buffer for boundary exchange
    isendbuflen        ! length of one column of sendbuf

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

USE utilities,        ONLY :   &
    tautsp2D             !

!------------------------------------------------------------------------------

USE meteo_utilities,  ONLY :   &
    moist_split

!------------------------------------------------------------------------------

USE environment,      ONLY :  &
    model_abort,       & ! aborts the program in case of errors
    exchg_boundaries     !

!------------------------------------------------------------------------------

USE parallel_utilities, ONLY:  &
    remark,            & !
    print_lm_field,    & !
    global_values,     & !
    lm_field_stats       !

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================
!+ Organizes the vertical interpolation of multi-level variables
!------------------------------------------------------------------------------

SUBROUTINE org_vert_interpol

!------------------------------------------------------------------------------
!
! Description:
!   org_vert_interpol organizes the vertical interpolations and the adaptation
!   of the surface pressure to the new LM/HM-orography.
!   For the variables qc_lm and qv_lm, the generalized relative humidity
!   is interpolated. After the interpolation, qc_lm and qv_lm are splitted
!   again.
!
! Method:
!
!------------------------------------------------------------------------------
!
! Subroutine arguments

! Local arrays
REAL (KIND=ireals)         ::  &
  zps1_lm (ie2lm,je2lm)                  ! first adaptation of surface pressure

INTEGER  (KIND=iintegers)  ::  &
  izdebug,              & ! for debug output
  izerror,              & ! status and error status variable
  kzgr,                 & ! top of boundary layer
  i, j, k, istat

REAL    (KIND=ireals)      ::  &
  zaq, zaqi, zbq, zbqi, zgqv, zpresm

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling

! Definition of statement functions
REAL (KIND=ireals) :: sf_psat_w, sf_psat_i, sf_qsat_ec, sf_qsat_gme,    &
                      x, y, z, zi, v, w, wi

sf_psat_w  (x,y,z,v,w) = y * EXP(z*(x-v)/(x-w))
sf_psat_i (x,y,zi,v,wi)= y * EXP(zi*(x-v)/(x-wi))

! ABOVE from ifs2gme file gme_comhumid.h:
!  Saturation water vapour pressure (with respect to ice, depending on the temp. "x")
! BELOW:
!  Specific humidity at saturation pressure (depending on the saturation water
!     vapour pressure "psatx" and the air pressure "px")
!  sf_qsat(psatx, px) = Rdv*MIN(psatx,0.8_ireals*px)/(px-O_m_rdv*MIN(psatx,0.8_ireals*px))
! IS DIFFERENT than the following definition:
sf_qsat_ec (x,y,z,v)   = z * x / MAX( (y-v*x), 1.0_ireals)    ! EC2LM -version
sf_qsat_gme(x,y,z,v)   = z * x / (y-v*x)                      ! GME2LM-version

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  izerror   = 0
  yzerrmsg  = '  '
  yzroutine = 'org_vert_interpol'

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
! Section 2: Compute generalized relative humidity
!------------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, 'in vert_interpol: compute grh'
  ENDIF

  IF (lgme2lm) THEN
    DO k = 1, ke_in
      DO j = 1, je2lm
        DO i = 1, ie2lm
          zpresm = akh_in(k) + bkh_in(k) * ps_gl(i,j)
          zaq    = sf_psat_w   (t_lm(i,j,k), b1, b2_w, b3, b4_w)
          zbq    = sf_qsat_gme (zaq, zpresm, Rdv, O_m_rdv)
          zgqv   = 1.0_ireals / zbq
          grh_lm(i,j,k) = (qv_lm(i,j,k) + qc_lm(i,j,k)) * zgqv
        ENDDO
      ENDDO
    ENDDO
  ELSEIF (lec2lm .OR. lcm2lm .OR. lgsm2lm .OR. lgfs2lm .OR. lhir2lm) THEN
    DO k = 1, ke_in
      DO j = 1, je2lm
        DO i = 1, ie2lm
          ! Changes by Anne Roches, MCH
          IF (qv_lm(i,j,k).LT.0.0_ireals) THEN
              qv_lm(i,j,k)=0.0_ireals
          ENDIF
          IF (qv_lm(i,j,k).GT.1.0_ireals) THEN
              qv_lm(i,j,k)=1.0_ireals
          ENDIF

          IF (qc_lm(i,j,k).LT.0.0_ireals) THEN
              qc_lm(i,j,k)=0.0_ireals
          ENDIF
          IF (qc_lm(i,j,k).GT.1.0_ireals) THEN
              qc_lm(i,j,k)=1.0_ireals
          ENDIF
          IF (lprog_qi) THEN
             IF (qi_lm(i,j,k).LT.0.0_ireals) THEN
                 qi_lm(i,j,k)=0.0_ireals
             ENDIF
             IF (qi_lm(i,j,k).GT.1.0_ireals) THEN
                 qi_lm(i,j,k)=1.0_ireals
             ENDIF
          ENDIF

          zpresm = akh_in(k) + bkh_in(k) * ps_gl(i,j)
          zaq    = sf_psat_w  (t_lm(i,j,k), B1, B2_w, B3, B4_w)
          zaqi   = sf_psat_i  (t_lm(i,j,k), B1, B2_i, B3, B4_i)
          zbq    = sf_qsat_ec (zaq, zpresm, Rdv, O_m_rdv)
          zbqi   = sf_qsat_ec (zaqi,zpresm, Rdv, O_m_rdv)
          zgqv   = 1.0_ireals / zbq

          ! Changes by Anne Roches, MCH
          IF (lmixcld) THEN
            IF (qv_lm(i,j,k) .GT. zbq) THEN
              grh_lm(i,j,k) = qv_lm(i,j,k) * zgqv
            ELSEIF ((qv_lm(i,j,k).LE. zbq) .AND. &
                    (qv_lm(i,j,k).GT. zbqi)) THEN
              grh_lm(i,j,k) = (qv_lm(i,j,k) + qc_lm(i,j,k)) * zgqv
            ELSE
              grh_lm(i,j,k) = (qv_lm(i,j,k) + qc_lm(i,j,k) + qi_lm(i,j,k)) * zgqv
            ENDIF
          ELSE IF (lprog_qi .AND. qv_lm(i,j,k) < zbqi) THEN
            grh_lm(i,j,k) = (qv_lm(i,j,k) + qc_lm(i,j,k) + qi_lm(i,j,k)) * zgqv
            qi_lm(i,j,k) = 0.0_ireals
          ELSE
            grh_lm(i,j,k) = (qv_lm(i,j,k) + qc_lm(i,j,k)) * zgqv
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Compute first adaptation of interpolated surface pressure
!------------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, 'in vert_interpol: first pressure adaptation'
  ENDIF

  CALL adapt_pressure_1 (grh_lm, zps1_lm, kzgr)

!------------------------------------------------------------------------------
! Section 4: Vertical interpolations for temperature and relative humidity
!------------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, 'in vert_interpol: vertical interpolations'
  ENDIF

  CALL vert_interpol (t_lm  , 't' , ke_in, akh_in, bkh_in, zps1_lm, kzgr)
  CALL vert_interpol (grh_lm, 'rh', ke_in, akh_in, bkh_in, zps1_lm, kzgr)
  IF (lprog_qi .AND. (.NOT. lmixcld)) THEN
    CALL vert_interpol (qi_lm, 'qi', ke_in, akh_in, bkh_in, zps1_lm, kzgr)
  ENDIF
  IF (lprog_qr_qs) THEN
    CALL vert_interpol (qr_lm, 'qr', ke_in, akh_in, bkh_in, zps1_lm, kzgr)
    CALL vert_interpol (qs_lm, 'qs', ke_in, akh_in, bkh_in, zps1_lm, kzgr)
  ENDIF
  IF (lprog_qg) THEN
    CALL vert_interpol (qg_lm, 'qg', ke_in, akh_in, bkh_in, zps1_lm, kzgr)
  ENDIF
! iso code
  IF (liso) THEN
    CALL vert_interpol (riso_lm(:,:,:,1), 'r18O', ke_in, akh_in, bkh_in,     &
                        zps1_lm, kzgr)
    CALL vert_interpol (riso_lm(:,:,:,2), 'r2H',  ke_in, akh_in, bkh_in,     &
                        zps1_lm, kzgr)
  ENDIF
! end iso code

!------------------------------------------------------------------------------
! Section 5: Compute second adaptation of surface pressure
!------------------------------------------------------------------------------

  CALL adapt_pressure_2 (grh_lm, zps1_lm, ke_in, ak_in, bk_in,          &
                                                 akh_in, bkh_in)

!------------------------------------------------------------------------------
! Section 6: Vertical interpolations for horizontal wind speeds
!------------------------------------------------------------------------------

  CALL vert_interpol (u_lm   , 'u' , ke_in, akh_in, bkh_in,   zps1_lm, kzgr)
  CALL vert_interpol (v_lm   , 'v' , ke_in, akh_in, bkh_in,   zps1_lm, kzgr)

!------------------------------------------------------------------------------
! Section 7: Correction of wind components
!------------------------------------------------------------------------------

  IF (luvcor .AND. .NOT.lbd_frame_cur) THEN
    CALL uv_correction (ke_in, ak_in, bk_in)
  ENDIF

!------------------------------------------------------------------------------
! Section 8: Computation of reduced surface pressure
!------------------------------------------------------------------------------

  IF (lprps) THEN
    CALL ps_to_sealevel (ke_in, ak_in, bk_in, kzgr)
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE org_vert_interpol

!==============================================================================

!+ First adaptation of surface pressure.
!------------------------------------------------------------------------------

SUBROUTINE adapt_pressure_1 (grh_lm, ps1_lm, kgr)

!------------------------------------------------------------------------------
!
! Description:
!  First hydrostatic adaptation of the horizontally interpolated GME surface
!  pressure to the LM/HM-orography (i.e. regarding the height difference
!  between the interpolated GME- and the LM/HM-orography).
!
! Method:
!
!------------------------------------------------------------------------------
!
! Parameterlist
REAL (KIND=ireals),       INTENT(IN)    ::   &
  grh_lm(ie2lm,je2lm,ke_in)       ! generalized relative humidity on GME-levels

REAL (KIND=ireals),       INTENT(OUT)   ::   &
  ps1_lm(ie2lm,je2lm)             ! for adapted surface pressure

INTEGER (KIND=iintegers), INTENT(OUT)   ::   &
  kgr                             ! top of boundary layer

! Local arrays:
REAL (KIND=ireals)                      ::   &
  ztv  (ie2lm,je2lm,ke_in),     & ! virtual temperature on GME-levels
  zpn  (ie2lm,je2lm,ke_in+1),   & ! pressure on GME main-levels
  zdfis(ie2lm,je2lm),           & !
  zbtv (ie2lm,je2lm),           & !
  zctv (ie2lm,je2lm),           & !
  zpq  (ie2lm,je2lm),           & !
  zpq2 (ie2lm,je2lm),           & !
  ztvq (ie2lm,je2lm),           & !
  ztvpq(ie2lm,je2lm),           & !
  zdtv (ie2lm,je2lm),           & !
  zfiu (ie2lm,je2lm),           & !
  zfio (ie2lm,je2lm)              !

! Local scalars:
INTEGER (KIND=iintegers)                ::   &
  i, j, k,                   & !
  izerror, izdebug
                            
REAL (KIND=ireals)                      ::   &
  zaq, zbq, zpdfis, zpht, zpgr, zpresm

REAL (KIND=ireals)                      ::   &
  zph(ie2lm,je2lm),   &
  zqv(ie2lm,je2lm),   &
  zqc(ie2lm,je2lm),   &
  zqi(ie2lm,je2lm)

CHARACTER (LEN=80)                         ::                              &
  yzerrmsg,            & ! for error message
  yzstring               ! for control output

CHARACTER (LEN=25)                         ::                              &
  yzroutine              ! character variable to be passed to remark

! Definition of statement functions
REAL    (KIND=ireals)    ::   sf_psat_w, sf_qsat_gme, sf_qsat_ec,          &
                              x, y, z, v, w

sf_psat_w  (x,y,z,v,w) = y * EXP(z*(x-v)/(x-w))
sf_qsat_gme(x,y,z,v)   = z * x / (y-v*x)                      ! GME2LM-version
sf_qsat_ec (x,y,z,v)   = z * x / MAX( (y-v*x), 1.0_ireals)    ! EC2LM -version

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  izerror   = 0
  yzerrmsg  = '   '
  yzroutine = 'adapt_pressure_1'

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  ! Compute virtual temperature and pressure on GME main levels

  IF (izdebug > 10) THEN
    PRINT *, 'in adapt_pressure_1: compute geopotential difference ',    &
              MINVAL(fis_gl(:,:)), MAXVAL(fis_gl(:,:)), MINVAL(fis_lm(:,:)), MAXVAL(fis_lm(:,:))
    DO k = 1, ke_in
      PRINT *, 'in adapt_pressure_1: what about inputs ', k,  &
              MINVAL(t_lm(:,:,k)), MAXVAL(t_lm(:,:,k)), akh_in(k), bkh_in(k), ak_in(k), bk_in(k)
    ENDDO
  ENDIF

  ! Keep surface pressure as lowest main level of GME pressure field
  ! and compute geopotential difference fis_gl - fis_lm
  DO j = 1, je2lm
    DO i = 1, ie2lm
      zpn    (i,j,ke_in+1) = ps_gl (i,j)
      zdfis  (i,j)         = fis_gl(i,j) - fis_lm(i,j)
    ENDDO
  ENDDO

  IF (izdebug > 10) THEN
    PRINT *, 'in adapt_pressure_1: compute pressure on half levels ',  &
              MINVAL(ps_gl(:,:)), MAXVAL(ps_gl(:,:)), ps_gl(5,5)
  ENDIF

  ! Compute the pressure on GME half levels (pn)
  ! and the virtual temperature (tv) on GME main levels
  DO k = 1, ke_in
    DO j = 1, je2lm
      DO i = 1, ie2lm
        zph(i,j)   = akh_in(k) + bkh_in(k)*ps_gl(i,j)
        zpn(i,j,k) = ak_in (k) + bk_in (k)*ps_gl(i,j)
        zaq = sf_psat_w (t_lm(i,j,k), b1, b2_w, b3, b4_w)

        ! Formulation with statement functions by Anne Roches, MCH
        IF (lgme2lm) THEN
          zbq = sf_qsat_gme (zaq, zph(i,j), Rdv, O_m_rdv)
        ELSE ! IF (lec2lm .OR. lcm2lm) THEN  ! but with ELSEIF, it does not vectorize!?
          zbq = sf_qsat_ec  (zaq, zph(i,j), Rdv, O_m_rdv)
        ENDIF

        zqv(i,j) =   MIN(1.0_ireals, grh_lm(i,j,k)) * zbq
        zqc(i,j) =   MAX(0.0_ireals, grh_lm(i,j,k)-1.0_ireals) * zbq
      ENDDO
    ENDDO

    IF ((lec2lm .OR. lcm2lm .OR. lgfs2lm .OR. lgsm2lm .OR. lhir2lm) .AND. lmixcld) THEN
      CALL moist_split(t_lm(:,:,k),zph,grh_lm(:,:,k),qvmin,qcmin,      &
                       qimin,pi,b1,b2_w,b2_i,b3,b4_w,b4_i,Rdv,O_m_rdv, &
                       zqv,zqc,zqi,ie2lm,je2lm)
      DO j = 1, je2lm
        DO i = 1, ie2lm
          ztv(i,j,k) = t_lm(i,j,k)*(1. + Rvd_m_o*zqv(i,j) - zqc(i,j) - zqi(i,j))
        ENDDO
      ENDDO
    ELSE
      DO j = 1, je2lm
        DO i = 1, ie2lm
          ztv(i,j,k) = t_lm(i,j,k)*(1. + Rvd_m_o*zqv(i,j) - zqc(i,j))
        ENDDO
      ENDDO
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
! Section 2: boundary layer and gradient of virtual temperature (zdtv)
!------------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, 'in adapt_pressure_1: find boundary layer height ',        &
               ztv (5,5,30), t_lm(5,5,30), zqv(5,5), zqc(5,5)
  ENDIF

! NOTE: the below computation is only valid, if the input vertical coordinate
!       is pressure based (which it is not if ivctype != 1)
!!$  ! Find the boundary layer height (kgr)
!!$  zpgr = 850.0E2
!!$  DO k = ke_in, 1, -1
!!$    zpht = akh_in(k) + bkh_in(k)*1.0E5_ireals
!!$    IF (zpht > zpgr) kgr = k
!!$  ENDDO
!!$  kgr = kgr - 1
  kgr = klv850_in

  ! Regression of the 4 layers above k = kgr for linear approximation of
  ! the virtual temperature (ztv(p) = zbtv*p + zctv)
  zpq  (:,:)  = 0.0_ireals
  zpq2 (:,:)  = 0.0_ireals
  ztvq (:,:)  = 0.0_ireals
  ztvpq(:,:)  = 0.0_ireals

  DO k = kgr-3, kgr
    DO j = 1, je2lm
      DO i = 1, ie2lm
        zpresm     = 0.5_ireals * (zpn(i,j,k+1) + zpn(i,j,k))
        zpq  (i,j) = zpq  (i,j) + zpresm
        zpq2 (i,j) = zpq2 (i,j) + zpresm * zpresm
        ztvq (i,j) = ztvq (i,j) + ztv(i,j,k)
        ztvpq(i,j) = ztvpq(i,j) + ztv(i,j,k) * zpresm
      ENDDO
    ENDDO
  ENDDO

  ! The gradient of the virtual temperature (zdtv) is computed from zbtv
  ! This is used for calculating the difference of the temperature
  ! regarding the change of the geopotential difference fis_gl - fis_lm.
  DO j = 1, je2lm
    DO i = 1, ie2lm
      zbtv(i,j) = (ztvpq(i,j) - ztvq(i,j) * zpq(i,j) / 4.0_ireals)  /     &
                   (zpq2(i,j) -  zpq(i,j) * zpq(i,j) / 4.0_ireals)
      zctv(i,j) = (ztvq (i,j) - zbtv(i,j) * zpq(i,j)) / 4.0_ireals
      zdtv(i,j) =  zbtv (i,j) * zpgr * zdfis  (i,j) /                     &
                     (R_d * (zbtv(i,j)*zpgr + zctv(i,j)))
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 3: Compute pressure adaptation
!------------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, 'in adapt_pressure_1: compute pressure adaptation '
  ENDIF

  zfiu (:,:)  = 0.0_ireals
  zpq  (:,:)  = 0.0_ireals
  ztvq (:,:)  = 0.0_ireals
  ztvpq(:,:)  = 0.0_ireals

  DO k = ke_in, 2, -1
    DO j = 1, je2lm
      DO i = 1, ie2lm
        zfio(i,j) = zfiu(i,j) + R_d*ztv(i,j,k)*LOG(zpn(i,j,k+1)/zpn(i,j,k))

        IF     (zdfis(i,j) < -10.0_ireals) THEN
          ! COSMO orography is higher than coarse grid orography
          ! Look for the proper layer in the coarse model, where the COSMO
          !  orography is in
          IF ( (ABS(zdfis  (i,j)) >  zfiu(i,j)) .AND.                        &
               (ABS(zdfis  (i,j)) <= zfio(i,j)) ) THEN
            ps1_lm(i,j) = zpn(i,j,k+1) *                                     &
              EXP( (ABS(zdfis(i,j)) - zfiu(i,j)) / (zfio (i,j)  - zfiu(i,j)) &
                      * LOG(zpn(i,j,k)/zpn(i,j,k+1))  )
          ENDIF

        ELSEIF (ABS(zdfis(i,j)) <= 10.0_ireals) THEN

          ! orographies are about the same, so the pressure is also
          ! this holds for all k, but ps1_lm has only to be set once
          IF (k == ke_in) THEN
            ps1_lm(i,j) = ps_gl(i,j)
          ENDIF

        ELSEIF (zdfis(i,j) > 10.0_ireals) THEN

          ! Compute ps1_lm for zdfis   > 0
          ! (coarse grid orography higher than COSMO orography)
          IF (zdfis  (i,j) > zfiu(i,j) .AND. zdfis  (i,j) <= zfio(i,j)) THEN
            ! Compute the average temperature (tvq) in the layer
            ! from fis_gl to fis_gl + zdfis   for zdfis   > 0
            zpdfis      = (zdfis  (i,j) - zfiu(i,j))/(zfio(i,j) - zfiu(i,j))*  &
                          (zpn(i,j,k) - zpn(i,j,k+1)) + zpn(i,j,k+1)
            ztvpq(i,j)  = ztv(i,j,k)*(zpn(i,j,k+1) - zpdfis)  + ztvpq(i,j)
            zpq  (i,j)  =             zpn(i,j,k+1) - zpdfis   + zpq  (i,j)
            ztvq (i,j)  = ztvpq(i,j)/zpq(i,j)
            ps1_lm(i,j) = ps_gl(i,j) *                                         &
                            EXP(zdfis  (i,j)/(R_d*(ztvq(i,j) + zdtv(i,j))))
          ELSE
            ztvpq(i,j)  = ztv(i,j,k)*(zpn(i,j,k+1) - zpn(i,j,k)) + ztvpq(i,j)
            zpq  (i,j)  =             zpn(i,j,k+1) - zpn(i,j,k)  + zpq  (i,j)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    zfiu(:,:) = zfio(:,:)
  ENDDO

!------------------------------------------------------------------------------
! Section 4: Additional control outputs
!------------------------------------------------------------------------------

  ! Control output of meanvalue, minimum and maximum of zdfis
  yzstring = '  Statistics for difference fis_gl - fis_lm (in m):'
  CALL lm_field_stats (zdfis, yzstring, 1.0_ireals/g, noutput,               &
                       yzerrmsg, izerror)

  ! Control output of meanvalue, minimum and maximum of ps_gl - ps1_lm
  zfio(:,:) = ps_gl(:,:) - ps1_lm(:,:)
  yzstring = '  Statistics for difference ps_gl - ps1_lm (in hPa):'
  CALL lm_field_stats (zfio, yzstring, 0.01_ireals, noutput, yzerrmsg,izerror)

  ! Control output of meanvalue, minimum and maximum of ps_gl
  yzstring = '  Statistics for ps_gl (in hPa):'
  CALL lm_field_stats (ps_gl, yzstring, 0.01_ireals, noutput,yzerrmsg,izerror)

  ! Control output of meanvalue, minimum and maximum of ps1_lm
  yzstring = '  Statistics for ps1_lm (in hPa):'
  CALL lm_field_stats (ps1_lm, yzstring, 0.01_ireals,noutput,yzerrmsg,izerror)

  ! Print some fields for debug/info purposes
  IF (lprps) THEN
    CALL print_lm_field (ps_gl, 'ps_gl', 'Pa ', ke1lm, 0.01_ireals,          &
         0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    CALL print_lm_field (fis_gl, 'fis_gl', 'gpm', ke1lm, 1.0_ireals/g,       &
         0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    CALL print_lm_field (fis_lm, 'fis_lm', 'gpm', ke1lm, 1.0_ireals/g,       &
         0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    CALL print_lm_field (zdfis, 'dfis_gl', 'gpm', ke1lm, 1.0_ireals/g,       &
         0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    CALL print_lm_field (ps1_lm, 'ps1_lm', 'Pa ', ke1lm, 0.01_ireals,        &
         0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE adapt_pressure_1

!==============================================================================

!+ Second adaptation of surface pressure.
!------------------------------------------------------------------------------

SUBROUTINE adapt_pressure_2 (grh_lm, ps1_lm, kevert, akhalf, bkhalf,     &
                                                     akmain, bkmain)

!------------------------------------------------------------------------------
!
! Description:
!  Second adaptation of the surface pressure using new values for temperature
!  and generalized relative humidity. The result is written to ps_lm.
!
!  Also, the generalized relative humidity is split again to the components
!  water vapor (qv_lm) and specific cloud water (qc_lm).
!
! Method:
!
!------------------------------------------------------------------------------
!
! Parameterlist
REAL (KIND=ireals),       INTENT(IN)    ::   &
  grh_lm(ie2lm,je2lm,ke_in)       ! generalized relative humidity on GME-levels

REAL (KIND=ireals),       INTENT(IN)    ::   &
  ps1_lm(ie2lm,je2lm)             ! for adapted surface pressure

INTEGER (KIND=iintegers), INTENT(IN)    ::   &
  kevert                          ! vertical dimension of coordinate parameters

REAL (KIND=ireals),       INTENT (IN)        ::  &
  akhalf(kevert+1), bkhalf(kevert+1), &! coordinate parameters for half levels
  akmain(kevert)  , bkmain(kevert)     ! coordinate parameters for main levels

! Local arrays:
REAL (KIND=ireals)                      ::   &
  ztv    (ie2lm,je2lm,kevert),   & ! virtual temperature on GME-levels
  zdfidps(ie2lm,je2lm),          & !
  zpsold (ie2lm,je2lm),          & !
  zpsnew (ie2lm,je2lm),          & !
  zficlm (ie2lm,je2lm),          & !
  zfiu   (ie2lm,je2lm),          & !
  zfio   (ie2lm,je2lm)             !

! Local scalars:
INTEGER (KIND=iintegers)                ::   &
  i, j, k,             & !
  izerror, iziter, iz1p, iz2p

REAL (KIND=ireals)                      ::   &
  zaq, zaqi, zbq, zbqi, zpno, zpnu

REAL (KIND=ireals)                      ::   &
  zph   (ie2lm,je2lm),     &
  zpc_fi(ie2lm,je2lm),     &
  zqv   (ie2lm,je2lm),     &
  zqc   (ie2lm,je2lm),     &
  zqi   (ie2lm,je2lm)

CHARACTER (LEN=80)                         ::                              &
  yzerrmsg,            & ! for error message
  yzstring               ! for control output

CHARACTER (LEN=25)                         ::                              &
  yzroutine              ! character variable to be passed to remark

! Definition of statement functions
REAL (KIND=ireals) :: sf_psat_w, sf_psat_i, sf_qsat_gme, sf_qsat_ec,       &
                      x, y, z, zi, v, w, wi

sf_psat_w  (x,y,z,v,w) = y * EXP(z*(x-v)/(x-w))
sf_psat_i (x,y,zi,v,wi)= y * EXP(zi*(x-v)/(x-wi))
sf_qsat_gme(x,y,z,v)   = z * x / (y-v*x)
sf_qsat_ec (x,y,z,v)   = z * x / MAX( (y-v*x), 1.0_ireals)

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  izerror   = 0
  yzerrmsg  = '   '
  yzroutine = 'adapt_pressure_2'

  ! Compute virtual temperature on main levels
  DO k = 2, kevert
    DO j = 1, je2lm
      DO i = 1, ie2lm
        zph(i,j) = akmain(k) + bkmain(k) * ps1_lm(i,j)
        zaq = sf_psat_w (t_lm(i,j,k), b1, b2_w, b3, b4_w)

        ! Formulation with statement functions by Anne Roches, MCH
        IF (lgme2lm) THEN
          zbq = sf_qsat_gme (zaq, zph(i,j), Rdv, O_m_rdv)
        ELSE ! IF (lec2lm .OR. lcm2lm) THEN  ! but with ELSEIF, it does not vectorize!?
          zbq = sf_qsat_ec  (zaq, zph(i,j), Rdv, O_m_rdv)
        ENDIF

        zqv(i,j) =   MIN(1.0_ireals, grh_lm(i,j,k)) * zbq
        zqc(i,j) =   MAX(0.0_ireals, grh_lm(i,j,k)-1.0_ireals) * zbq
      ENDDO
    ENDDO

    IF ((lec2lm .OR. lcm2lm .OR. lgfs2lm .OR. lgsm2lm) .AND. lmixcld) THEN
      CALL moist_split(t_lm(:,:,k),zph,grh_lm(:,:,k),qvmin,qcmin,      &
                       qimin,pi,b1,b2_w,b2_i,b3,b4_w,b4_i,Rdv,O_m_rdv, &
                       zqv,zqc,zqi,ie2lm,je2lm)
      DO j = 1, je2lm
        DO i = 1, ie2lm
          ztv(i,j,k) = t_lm(i,j,k)*(1. + Rvd_m_o*zqv(i,j) - zqc(i,j) - zqi(i,j))
        ENDDO
      ENDDO
    ELSE
      DO j = 1, je2lm
        DO i = 1, ie2lm
          ztv(i,j,k) = t_lm(i,j,k)*(1. + Rvd_m_o*zqv(i,j) - zqc(i,j))
        ENDDO
      ENDDO
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
! Section 2: Compute ps_lm with Newton iterations
!------------------------------------------------------------------------------

  iziter      = 0
  zpsold(:,:) = ps1_lm(:,:)

  iteration: DO WHILE (iziter < 5)
    iziter = iziter + 1
    zfiu   (:,:) = fis_lm(:,:)
    zficlm (:,:) = 0.0_ireals
    zdfidps(:,:) = 0.0_ireals

    IF (kcontrol_fi /= -1) THEN
      DO j = 1, je2lm
        DO i = 1, ie2lm
          zpc_fi(i,j)  = ak_in(kcontrol_fi) + bk_in(kcontrol_fi) * ps_gl(i,j)
        ENDDO
      ENDDO
    ELSE
      DO j = 1, je2lm
        DO i = 1, ie2lm
          zpc_fi(i,j)  = pcontrol_fi
        ENDDO
      ENDDO
    ENDIF

    DO k = kevert, 2, - 1
      DO j = 1, je2lm
        DO i = 1, ie2lm
          zpno         = akhalf(k)    + bkhalf(k)   * zpsold(i,j)
          zpnu         = akhalf(k+1)  + bkhalf(k+1) * zpsold(i,j)
          zfio   (i,j) = zfiu  (i,j)  + R_d * ztv(i,j,k) * LOG(zpnu/zpno)
          zdfidps(i,j) = zdfidps(i,j) + R_d * ztv(i,j,k) *             &
                         (bkhalf(k+1)*zpno - bkhalf(k)*zpnu)/(zpnu*zpno)
          IF ( zpnu > zpc_fi(i,j) .AND. zpno <= zpc_fi(i,j)) THEN
            zficlm (i,j) = zfiu   (i,j) + R_d*ztv(i,j,k)*LOG(zpnu/zpc_fi(i,j))
            zdfidps(i,j) = zdfidps(i,j) + R_d*ztv(i,j,k)*bkhalf(k+1)/zpnu
          ENDIF
        ENDDO
      ENDDO
      zfiu(:,:) = zfio(:,:)
    END DO
    zpsnew(:,:) = zpsold(:,:) - (zficlm(:,:)-fic_gl(:,:))/zdfidps(:,:)
    zpsold(:,:) = zpsnew(:,:)
  ENDDO iteration

  ps_lm(:,:)  = zpsnew(:,:)

!------------------------------------------------------------------------------
! Section 3: Split generalized relative humidity
!------------------------------------------------------------------------------

  ! Use the new surface pressure ps_lm
  IF (lgme2lm) THEN
    DO k = 1, kevert
      DO j = 1, je2lm
        DO i = 1, ie2lm
          zph(i,j) = akmain(k) + bkmain(k) * ps_lm(i,j)
          zaq = sf_psat_w   (t_lm(i,j,k), b1, b2_w, b3, b4_w)
          zbq = sf_qsat_gme (zaq, zph(i,j), Rdv, O_m_rdv)
          qv_lm(i,j,k) =   MIN(1.0_ireals, grh_lm(i,j,k)) * zbq
          qc_lm(i,j,k) =   MAX(0.0_ireals, grh_lm(i,j,k)-1.0_ireals) * zbq
          IF (qv_lm(i,j,k) < qvmin) THEN
            qv_lm(i,j,k)  = qvmin
          ENDIF
          IF (qc_lm(i,j,k) < qcmin) THEN
            qc_lm(i,j,k)  = 0.0_ireals
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ELSEIF (lec2lm .OR. lcm2lm .OR. lgsm2lm .OR. lgfs2lm .OR. lhir2lm) THEN
    DO k = 1, kevert
      DO j = 1, je2lm
        DO i = 1, ie2lm
          zph(i,j) = akmain(k) + bkmain(k) * ps_lm(i,j)
          zaq = sf_psat_w (t_lm(i,j,k), B1, B2_w, B3, B4_w)
          zaqi= sf_psat_i (t_lm(i,j,k), B1, B2_i, B3, B4_i)
          zbq = sf_qsat_ec(zaq, zph(i,j), Rdv, O_m_rdv)
          zbqi= sf_qsat_ec(zaqi, zph(i,j), Rdv, O_m_rdv)
          qv_lm(i,j,k) =   MIN(1.0_ireals, grh_lm(i,j,k)) * zbq
          qc_lm(i,j,k) =   MAX(0.0_ireals, grh_lm(i,j,k)-1.0_ireals) * zbq
          IF (qv_lm(i,j,k) < qvmin) THEN
            qv_lm(i,j,k)  = qvmin
          ENDIF
          IF (qc_lm(i,j,k) < qcmin) THEN
            qc_lm(i,j,k)  = 0.0_ireals
          ENDIF
          IF (.NOT.lmixcld .AND. lprog_qi .AND. qv_lm(i,j,k) < zbqi) THEN
            qi_lm(i,j,k)  = 0.0_ireals
          END IF
        ENDDO
      ENDDO
      IF (lmixcld) THEN
        CALL moist_split(t_lm(:,:,k),zph,grh_lm(:,:,k),qvmin,qcmin,      &
                         qimin,pi,b1,b2_w,b2_i,b3,b4_w,b4_i,Rdv,O_m_rdv, &
                         qv_lm(:,:,k),qc_lm(:,:,k),qi_lm(:,:,k),         &
                         ie2lm,je2lm)
      ENDIF
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Additional control outputs
!------------------------------------------------------------------------------

  ! Control output of meanvalue, minimum and maximum of ps1_lm - ps_lm
  yzstring = '  Statistics for difference ps1_lm - ps_lm (in hPa):'
  zpsnew(:,:) = ps1_lm(:,:) - ps_lm(:,:)
  CALL lm_field_stats (zpsnew, yzstring, 0.01_ireals,noutput,yzerrmsg,izerror)

  ! Control output of meanvalue, minimum and maximum of ps_lm
  yzstring = '  Statistics for ps_lm (in hPa):'
  CALL lm_field_stats (ps_lm, yzstring, 0.01_ireals, noutput,yzerrmsg,izerror)

  ! Print some fields for debug/info purposes
  IF (lprps) THEN
    CALL print_lm_field (zpsnew, 'ps1_lm-ps_lm', 'Pa', kevert+1,0.01_ireals, &
         0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    CALL print_lm_field (ps_lm, 'ps_lm', 'Pa ', kevert+1, 0.01_ireals,       &
         0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)
  ENDIF

  ! For t_lm, qv_lm, qc_lm and grh_lm print only 2 levels
  nlev1pr = MAX (1_iintegers, nlev1pr)
  iz1p    = MIN (     kevert, nlev1pr)
  nlev2pr = MAX (1_iintegers, nlev2pr)
  iz2p    = MIN (     kevert, nlev2pr)
  IF (lprt) THEN
    CALL print_lm_field (t_lm(:,:,iz1p), 't_lm', 'Grad C', iz1p, 10.0_ireals,&
      -273.15_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    CALL print_lm_field (t_lm(:,:,iz2p), 't_lm', 'Grad C', iz2p, 10.0_ireals,&
      -273.15_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)
  ENDIF
  IF (lprgrh) THEN
    CALL print_lm_field (grh_lm(:,:,iz1p), 'grh_lm', '%', iz1p, 100.0_ireals,&
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    CALL print_lm_field (grh_lm(:,:,iz2p), 'grh_lm', '%', iz2p, 100.0_ireals,&
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)
  ENDIF
  IF (lprqv) THEN
    CALL print_lm_field (qv_lm(:,:,iz1p),'qv_lm','kg/kg',iz1p, 1.0E4_ireals, &
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    CALL print_lm_field (qv_lm(:,:,iz2p),'qv_lm','kg/kg',iz2p, 1.0E4_ireals, &
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)
  ENDIF
  IF (lprqc) THEN
    CALL print_lm_field (qc_lm(:,:,iz1p),'qc_lm','kg/kg',iz1p, 1.0E4_ireals, &
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    CALL print_lm_field (qc_lm(:,:,iz2p),'qc_lm','kg/kg',iz2p, 1.0E4_ireals, &
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)
  END IF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE adapt_pressure_2

!==============================================================================

!+ Vertical Interpolation of atmospheric variables.
!------------------------------------------------------------------------------

SUBROUTINE vert_interpol (xlm, yfld, kevert, akmain, bkmain, ps1_lm, kgr)

!------------------------------------------------------------------------------
!
! Description:
!   Temperature, generalized relative humidity, cloud ice content and the
!   wind components are
!   interpolated from the GME-levels on to levels specified by the vertical
!   coordinate parameters akmain, bkmain and the horizontally interpolated
!   and adapted surface pressure ps_gl (temp, hum) or ps_lm (u, v) (i.e.,
!   the new levels are adapted to the HM/LM-orography).
!
!   If LM-fields are computed, the coordinate parameters from GME are still
!   used, because this is only an intermediate step and vertical interpolation
!   onto the LM-levels is done later.
!
!   If HM-fields are computed, the coordinate parameters of the regional
!   model are used.
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
  xlm(ie2lm, je2lm, kedim)   ! field that has to be interpolated

REAL (KIND=ireals),       INTENT (IN)        ::  &
  akmain(kevert), bkmain(kevert), &! coordinate parameters
  ps1_lm(ie2lm,je2lm)              ! first adapted surface pressure

! Scalar arguments with intent(in):
CHARACTER (LEN= *),       INTENT(IN)         ::  &
  yfld                  ! name of the field

!------------------------------------------------------------------------------
!
! Local variables
INTEGER (KIND=iintegers)       ::     &
  kzint_vec(ie2lm), izln_vec(ie2lm), kzgrn(ie2lm)

INTEGER (KIND=iintegers)       ::     &
  i, j, k, iz1p, iz2p, izn(ie2lm), izerror, izdebug, idone

LOGICAL                        ::     &
  ldone(ie2lm)

REAL    (KIND=ireals)          ::     &
  zgamma, zdelpchk, zdpx, zphgl, zbq, zprod

! Local arrays:
INTEGER (KIND=iintegers)       ::     &
  izindex (ie2lm,kevert), &! indices used during interpolation
                           ! (vertical dimension of the target)
  kzdims(24)               ! vertical dimensions for boundary exchange

REAL (KIND=ireals)             ::     &
  zrhmax    (ie2lm),             & !
  zdelp     (ie2lm,je2lm),       & ! difference: interpolated and computed ps
  zxexp_vec (ie2lm,ke_in+4),     & ! values for interpolations
  zpexp_vec (ie2lm,ke_in+4),     & ! points where above values are valid
  zbreak_vec(ie2lm,(ke_in+4)*3), & ! work array for tension spline routine
  zs_vec    (ie2lm,(ke_in+4)*6), & ! work array for tension spline routine
  zphgls    (ie2lm,ke_in+4),     & ! pressure on vertical levels starting from
                           ! interpolated surface pressure
                           ! (all in vertical dimension of input model)
  zpmain    (ie2lm,kevert) ! pressure on vertical levels starting from
                           ! adapted surface pressure
                           ! (in vertical dimension of the target)

REAL (KIND=ireals)             ::     &
  zcoef_vec(ie2lm,4,3*(ke_in+4)),& ! work array for tension spline routine
  zpsxgl (ie2lm,je2lm),  & ! interpolated surface pressure (for mass-, u- or)
  zpsxlm (ie2lm,je2lm),  & ! adapted surface pressure      (v-grid point    )
  zbx    (ie2lm,je2lm),  & !
  zdelx  (ie2lm,je2lm),  & !
  zpq    (ie2lm,je2lm),  & !                work arrays
  zpq2   (ie2lm,je2lm),  & !
  zxq    (ie2lm,je2lm),  & !
  zxpq   (ie2lm,je2lm)

CHARACTER (LEN=80)             ::     &
  yzerrmsg

CHARACTER (LEN=25)             ::     &
  yzroutine

!- End of header
!------------------------------------------------------------------------------

  yzroutine = 'vert_interpol'
  izerror   = 0

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
! Section 1: Prepare surface pressure fields (evtl. at u or v gridpoints)
!------------------------------------------------------------------------------
! zpsxgl   is the horizontaly interpolated surface pressure of the GME,
! zpsxlm   is the surface pressure adapted to the LM orography after
!          first (ps1_lm) or the second (ps_lm) pressure correction.
!------------------------------------------------------------------------------

  IF (izdebug > 19) THEN
    PRINT *, '      Vertical interpolation for ', yfld
  ENDIF

  IF ( yfld == 'u' ) THEN

    DO j = 1, je2lm
      DO i = 1, ie2lm-1
        zpsxgl(i,j) = 0.5*(ps_gl(i,j) + ps_gl(i+1,j))
        zpsxlm(i,j) = 0.5*(ps_lm(i,j) + ps_lm(i+1,j))
      ENDDO
    ENDDO
    IF (my_cart_neigh(3) == -1) THEN
      zpsxgl(ie2lm,:) = ps_gl(ie2lm-1,:)
      zpsxlm(ie2lm,:) = ps_lm(ie2lm-1,:)
    ENDIF

    IF (num_compute > 1) THEN
      kzdims=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                              &
         (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie2lm, je2lm, kzdims, jstartpar, jendpar, 1, nboundlines,      &
          my_cart_neigh, .FALSE., .FALSE., .FALSE.,                      &
          111, .FALSE., ncomm_type, izerror, yzerrmsg,                   &
          zpsxgl, zpsxlm)
    ENDIF

  ELSE IF ( yfld == 'v' ) THEN

    DO j = 1, je2lm-1
      DO i = 1, ie2lm
        zpsxgl(i,j) = 0.5*(ps_gl(i,j) + ps_gl(i,j+1))
        zpsxlm(i,j) = 0.5*(ps_lm(i,j) + ps_lm(i,j+1))
      ENDDO
    ENDDO
    IF (my_cart_neigh(2) == -1) THEN
      zpsxgl(:,je2lm) = ps_gl(:,je2lm-1)
      zpsxlm(:,je2lm) = ps_lm(:,je2lm-1)
    ENDIF

    IF (num_compute > 1) THEN
      kzdims=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                              &
         (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie2lm, je2lm, kzdims, jstartpar, jendpar, 1, nboundlines,      &
          my_cart_neigh, .FALSE., .FALSE., .FALSE.,                      &
          112, .FALSE., ncomm_type, izerror, yzerrmsg,                   &
          zpsxgl, zpsxlm)
    ENDIF

  ELSE    ! yfld == 't', 'rh', 'qi', 'qr', 'qs', 'qg'

    zpsxgl(:,:) = ps_gl (:,:)
    zpsxlm(:,:) = ps1_lm(:,:)

  ENDIF

!------------------------------------------------------------------------------
! Section 2: Prepare surface pressure fields (evtl. at u or v gridpoints)
!------------------------------------------------------------------------------
! For all the fields except the relative humidity a linear approximation
! of the vertical profile above the boundary layer is computed.
! This linear profile is used for the extrapolation of the fields.
!
! The index 'kgr' of the boundary layer top is computed in subroutine
! adapt_pressure_1.
!
! The linear approximation is of the form:
!  x(p) = bx*p + cx  ,  where p is the pressure at the full levels.
!
! We need only the coefficient bx and delx . We do not need cx.
!------------------------------------------------------------------------------

  IF ( (yfld == 'u') .OR. (yfld == 'v') .OR. (yfld == 't') ) THEN
    zpq  (:,:) = 0.0
    zpq2 (:,:) = 0.0
    zxq  (:,:) = 0.0
    zxpq (:,:) = 0.0
    DO k = kgr-3, kgr
      zpq (:,:) = zpq (:,:) + (akmain(k) + bkmain(k)*zpsxgl(:,:))
      zpq2(:,:) = zpq2(:,:) + (akmain(k) + bkmain(k)*zpsxgl(:,:)) &
                             *(akmain(k) + bkmain(k)*zpsxgl(:,:))
      zxq (:,:) = zxq (:,:) + xlm(:,:,k)
      zxpq(:,:) = zxpq(:,:) + xlm(:,:,k) * (akmain(k) + bkmain(k)*zpsxgl(:,:))
    ENDDO
    zbx  (:,:) = ( zxpq(:,:) - zxq(:,:)*zpq(:,:)/4.0_ireals) /    &
                  (zpq2(:,:) - zpq(:,:)**2/4.0_ireals)
    zdelx(:,:) =  zbx (:,:) * (zpsxlm(:,:) - zpsxgl(:,:))
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Keep horizontaly interpolated GME profiles at some grid points
!------------------------------------------------------------------------------

! IF (lprgp) THEN
!   DO izn = 1, ngp
!     DO k = 1, kevert
!       IF ( yfld == 't' ) THEN
!           tglpr(k,izn) = xlm(igp(izn),jgp(izn),k)
!       ELSE IF ( yfld == 'u' ) THEN
!           uglpr(k,izn) = xlm(igp(izn),jgp(izn),k)
!       ELSE IF ( yfld == 'v' ) THEN
!           vglpr(k,izn) = xlm(igp(izn),jgp(izn),k)
!       ELSE IF ( yfld == 'rh' ) THEN
!         grhglpr(k,izn) = xlm(igp(izn),jgp(izn),k)
!       ENDIF
!     ENDDO
!   ENDDO
! ENDIF

!------------------------------------------------------------------------------
! Section 4: Vertical interpolation from GME levels with interpolated surface
!            pressure to GME levels with the surface pressure adapted to the
!            LM orography.
!------------------------------------------------------------------------------

  zgamma   =  5.5_ireals
  zdelpchk = 10.0_ireals

  DO j = 1, je2lm

    IF (izdebug > 39) THEN
      PRINT *, '      Vertical interpolation for ', yfld, ';  j = ', j
    ENDIF

    !------------------------------------------------------------------------
    ! Section 4.1: Preparations
    !------------------------------------------------------------------------

    DO i = 1, ie2lm
      ! Surface pressure difference LM minus GME (horizontaly interpolated)
      zdelp(i,j) = zpsxlm(i,j) - zpsxgl(i,j)
    ENDDO

    ! Construct input model profile (the vertical coordinate parameters from
    ! input model are really needed here)
    DO k = 1, ke_in
      DO i = 1, ie2lm
        zxexp_vec (i,k) = xlm (i,j,k)
        zphgls(i,k)     = akh_in(k) + bkh_in(k)*zpsxgl(i,j)
        zpexp_vec (i,k) = zphgls(i,k)
      ENDDO
    ENDDO
    DO i = 1, ie2lm
      zxexp_vec (i,ke_in+ 1) = zxexp_vec(i,ke_in)
      zpexp_vec (i,ke_in+ 1) = zpsxgl(i,j)
      kzint_vec(i) = ke_in + 1
    ENDDO

    !  for yfld= 't', 'u', 'v' with a linear regression
    !  for yfld = 'rh'         with a constant value.
    ! For a better vectorization, the loops have been changed

    IF ( (yfld == 'u') .OR. (yfld == 'v') .OR. (yfld == 't') ) THEN
      DO k = kgr + 1, kgr + 3
        DO i = 1, ie2lm
          IF (zpsxlm(i,j) > zpsxgl(i,j) + zdelpchk) THEN
            zxexp_vec(i,k) = xlm (i,j,kgr) + zbx(i,j)*zdelp(i,j)/3.0*REAL(k-kgr,ireals)
            zpexp_vec(i,k) = zphgls(i,kgr) + zdelp(i,j)/3.0*REAL(k-kgr, ireals)
          ENDIF
        ENDDO
      ENDDO
      DO k = kgr + 4, ke_in + 3
        DO i = 1, ie2lm
          IF (zpsxlm(i,j) > zpsxgl(i,j) + zdelpchk) THEN
            zxexp_vec(i,k) = xlm (i,j,k-3) + zdelx(i,j)
            zpexp_vec(i,k) = zphgls(i,k-3) + zdelp(i,j)
          ENDIF
        ENDDO
      ENDDO
    ELSE   ! 'rh' and 'qi', 'qr', 'qs', 'qg'
      DO k = kgr + 1, kgr + 3
        DO i = 1, ie2lm
          IF (zpsxlm(i,j) > zpsxgl(i,j) + zdelpchk) THEN
            zxexp_vec(i,k) = xlm (i,j,kgr)
            zpexp_vec(i,k) = zphgls(i,kgr) + zdelp(i,j)/3.0*REAL(k-kgr, ireals)
          ENDIF
        ENDDO
      ENDDO
      DO k = kgr + 4, ke_in + 3
        DO i = 1, ie2lm
          IF (zpsxlm(i,j) > zpsxgl(i,j) + zdelpchk) THEN
            zxexp_vec(i,k) = xlm (i,j,k-3)
            zpexp_vec(i,k) = zphgls(i,k-3) + zdelp(i,j)
          ENDIF
        ENDDO
      ENDDO
    ENDIF ! yfld

    DO i = 1, ie2lm
      IF (zpsxlm(i,j) > zpsxgl(i,j) + zdelpchk) THEN
        zxexp_vec  (i,ke_in + 4) = zxexp_vec(i,ke_in+3)
        zpexp_vec  (i,ke_in + 4) = zpsxgl(i,j) + zdelp(i,j)
        kzint_vec  (i)           = ke_in + 4
      ENDIF
    ENDDO

    ! Elevate (shift) the boundary layer for zpsxlm < zpsxgl
    ! (aber so wird es immer noch nicht vektorisieren, wie lange dauert es?
    DO i = 1, ie2lm
      IF (zpsxlm(i,j) < zpsxgl(i,j) - zdelpchk) THEN
        k = kgr
        DO WHILE (zpexp_vec(i,k) > zphgls(i,kgr) + zdelp(i,j))
          k = k - 1
        ENDDO
        kzgrn(i) = k
      ENDIF
    ENDDO

    IF ( (yfld == 'u') .OR. (yfld == 'v') .OR. (yfld == 't') ) THEN
      DO k = 1, ke_in - kgr + 1
        DO i = 1, ie2lm
          IF (zpsxlm(i,j) < zpsxgl(i,j) - zdelpchk) THEN
            zxexp_vec(i,kzgrn(i) + k) = xlm (i,j,kgr + k - 1) + zbx(i,j)*zdelp(i,j)
            zpexp_vec(i,kzgrn(i) + k) = zphgls(i,kgr + k - 1) +          zdelp(i,j)
          ENDIF
        ENDDO
      ENDDO
    ELSE     ! 'rh' and 'qi', 'qr', 'qs', 'qg'
      DO k = 1, ke_in - kgr + 1
        DO i = 1, ie2lm
          IF (zpsxlm(i,j) < zpsxgl(i,j) - zdelpchk) THEN
            zxexp_vec(i,kzgrn(i) + k) = xlm (i,j,kgr + k - 1)
            zpexp_vec(i,kzgrn(i) + k) = zphgls(i,kgr + k - 1) + zdelp(i,j)
          ENDIF
        ENDDO
      ENDDO
    ENDIF ! yfld

    DO i = 1, ie2lm
      IF (zpsxlm(i,j) < zpsxgl(i,j) - zdelpchk) THEN
        kzint_vec(i)              = kzgrn(i) + ke_in - kgr + 2
        zxexp_vec(i,kzint_vec(i)) = zxexp_vec(i,kzint_vec(i)-1)
        zpexp_vec(i,kzint_vec(i)) = zpsxgl(i,j) + zdelp(i,j)
      ENDIF
    ENDDO

    !------------------------------------------------------------------------
    ! Section 4.2: Compute tension splines
    !------------------------------------------------------------------------

    ! Interpolation with Tension-Splines to COSMO Model levels
    DO i = 1, ie2lm
      ! izln must be defined before !!!!
      izln_vec(i) = (ke_in+4)*3
      ! for later limiting of rh
      zrhmax  (i) = 0.0_ireals
    ENDDO ! i

    CALL tautsp2D(zpexp_vec, zxexp_vec, kzint_vec, ie2lm, 1,ie2lm, ke_in+4, &
                  zgamma, zs_vec, zbreak_vec, zcoef_vec, izln_vec, izerror)

    IF (izerror == 0) THEN

      DO k = 1, kevert
        DO i = 1, ie2lm
          ! Pressure on the full levels and at the right position
          zpmain(i,k) = akmain(k) + bkmain(k) * zpsxlm(i,j)
        ENDDO
      ENDDO

      ! check pressure at uppermost level (may occure in case of sigma coordinates )
      DO i = 1, ie2lm
        IF (zpmain(i,1) < zpexp_vec(i,1)) THEN
          zpmain(i,1) = zpexp_vec(i,1) + 1.0_ireals
        ENDIF
      ENDDO

      izn(:) = 1
      DO k = 1, kevert
        ldone(:) = .FALSE.
        idone    = 0
        DO WHILE (idone < ie2lm)
          DO i = 1, ie2lm
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
        DO i = 1, ie2lm
          zdpx = zpmain(i,k) - zbreak_vec(i,izindex(i,k))
          xlm(i,j,k) =  zcoef_vec(i,1,izindex(i,k))+zdpx  * (zcoef_vec(i,2,izindex(i,k))+  &
              zdpx*0.5*(zcoef_vec(i,3,izindex(i,k))+zdpx/3.0*zcoef_vec(i,4,izindex(i,k))))
        ENDDO
      ENDDO

      IF ( yfld == 'rh' ) THEN
        DO k=1, ke_in + 4   ! but limit to kzint_vec(i) later
          DO i = 1, ie2lm
            ! Limit generalized relative humidity to 0.001 .. rhmax
            IF (k <= kzint_vec(i)) THEN
              zrhmax(i) = MAX ( zrhmax(i), zxexp_vec(i,k) )
            ENDIF
          ENDDO
        ENDDO

        DO k = 1, kevert
          DO i = 1, ie2lm
            xlm(i,j,k) = MIN(zrhmax(i),    xlm(i,j,k))
            xlm(i,j,k) = MAX(0.001_ireals, xlm(i,j,k))
          ENDDO
        ENDDO
      ENDIF

      ! Limit values of qi, qr, qs, qg
      IF ((yfld=='qi') .OR. (yfld=='qr') .OR. (yfld=='qs') .OR. (yfld=='qg')) THEN
        DO k = 1, kevert
          DO i = 1, ie2lm
            IF (xlm(i,j,k) < qimin) THEN
              xlm(i,j,k) = 0.0_ireals
            ENDIF
          ENDDO
        ENDDO
      ENDIF

! iso code
      ! Limit isotope ratios
      IF ((yfld=='r18O') .OR. (yfld=='r2H')) THEN
        DO k = 1, kevert
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

    ELSE
      yzerrmsg = 'Error in tautsp'
      PRINT *, '*** ERROR in tautsp2D: while processing j-index:  ', j
      CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
  ENDDO ! j

!------------------------------------------------------------------------------
! Section 5: Keep vertical interpolated HM/LM-profiles at some grid points
!------------------------------------------------------------------------------

! IF (lprgp) THEN
!
!   IF ( yfld /= 'rh' ) THEN
!     DO izn = 1, ngp
!       DO k = 1, kevert
!         IF ( yfld == 't' ) THEN
!           tpr(k,izn) = xlm(igp(izn),jgp(izn),k)
!         ELSE IF ( yfld == 'u' ) THEN
!           upr(k,izn) = xlm(igp(izn),jgp(izn),k)
!         ELSE IF ( yfld == 'v' ) THEN
!           vpr(k,izn) = xlm(igp(izn),jgp(izn),k)
!         ENDIF
!       ENDDO
!     ENDDO
!
!   ELSE IF ( yfld == 'rh' ) THEN
!
!   ! Use of the already verticaly interpolated temperature field
!     DO izn = 1, ngp
!       DO k = 1, kevert
!         zphgl = akmain(k) + bkmain(k)*ps1_lm(igp(izn),jgp(izn))
!         grhpr(k,l) =                     xlm(igp(izn),jgp(izn),k)
!         zbq = sf_qsat (sf_psat_w (t_lm(igp(izn),jgp(izn),k), &
!                                   B1, B2_w, B3, B4_w), zphgl, Rdv, O_m_rdv)
!         qvpr(k,izn) = MIN(1.0_ireals, grhpr(k,izn)) * zbq
!         qcpr(k,izn) = MAX(0.0_ireals, grhpr(k,izn)-1.0_ireals) * zbq
!       ENDDO
!     ENDDO
!
!   ENDIF ! yfld
! ENDIF ! lprgp

  ! Print the wind fields if wanted
  !--------------------------------

  IF (lpru .AND. yfld == 'u') THEN
    nlev1pr = MAX (1_iintegers, nlev1pr)
    iz1p    = MIN (     kevert, nlev1pr)
    CALL print_lm_field (u_lm(:,:,iz1p), 'u_lm', 'm/s', iz1p, 10.0_ireals,   &
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    nlev2pr = MAX (1_iintegers, nlev2pr)
    iz2p    = MIN (     kevert, nlev2pr)
    CALL print_lm_field (u_lm(:,:,iz2p), 'u_lm', 'm/s', iz2p, 10.0_ireals,&
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)
  ENDIF
  IF (lprv .AND. yfld == 'v') THEN
    nlev1pr = MAX (1_iintegers, nlev1pr)
    iz1p    = MIN (     kevert, nlev1pr)
    CALL print_lm_field (v_lm(:,:,iz1p), 'v_lm', 'm/s', iz1p, 10.0_ireals,   &
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

    nlev2pr = MAX (1_iintegers, nlev2pr)
    iz2p    = MIN (     kevert, nlev2pr)
    CALL print_lm_field (v_lm(:,:,iz2p), 'v_lm', 'm/s', iz2p, 10.0_ireals,&
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)
  ENDIF
  
!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE vert_interpol

!==============================================================================

!+ Wind correction.
!------------------------------------------------------------------------------

SUBROUTINE uv_correction (kevert, akhalf, bkhalf)

!------------------------------------------------------------------------------
!
! Description:
!   Wind correction with a divergent wind potential.
!
! Method:
!
!------------------------------------------------------------------------------
!
! Parameterlist
INTEGER (KIND=iintegers), INTENT(IN)   ::  &
  kevert               ! number of levels of target model

REAL    (KIND=ireals)   , INTENT(IN)   ::  &
  akhalf(kevert+1),  & ! coordinate parameters for half levels
  bkhalf(kevert+1)     !            --- " ---

! Local arrays:
REAL    (KIND=ireals)  ::    &
  zcoslat_m(je2lm),          & ! for mass- and u-gridpoints
  zcoslat_v(je2lm),          & ! for v-gridpoints
  zsinlat  (je2lm),          & !
  zwk      (kevert)            !

REAL    (KIND=ireals)  ::    &
  zchi_new(ie2lm,je2lm),     & !
  zchi    (ie2lm,je2lm),     & !
  zal     (ie2lm,je2lm),     & !
  zcl     (ie2lm,je2lm),     & !
  zalp    (ie2lm,je2lm),     & !
  zclp    (ie2lm,je2lm),     & !
  zbb     (ie2lm,je2lm),     & !
  zdelpq  (ie2lm,je2lm),     & !
  zdelp   (ie2lm,je2lm),     & !
  zrhs    (ie2lm,je2lm),     & !
  zud     (ie2lm,je2lm),     & !
  zvd     (ie2lm,je2lm),     & !
  zdpdtlm (ie2lm,je2lm),     & !
  zdpdtlm1(ie2lm,je2lm)        !

! Local scalars:
INTEGER (KIND=iintegers) ::  &
  i, j, k, iziter, itag,     & !
  izerror, kzdims(24)

REAL    (KIND=ireals)    ::  &
  zlat, zdlonr, zdlatr, zul, zur, zvo, zvu, zalsor, zdmax, zdmax1, zres,  &
  zdchi

CHARACTER (LEN=80)             ::     &
  yzerrmsg, yzstring

CHARACTER (LEN=25)             ::     &
  yzroutine

!
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  izerror   = 0
  yzerrmsg  = '   '
  yzroutine = 'uv_correction'

  ! Compute cos(lat) and sin(lat) terms (at mass-, u- and v-gridpoints)
  DO j = 1, je2lm
    zlat = (startlat + (j - 1)*dlat)*degrad
    zcoslat_m(j) = COS(zlat)
    zcoslat_v(j) = COS(zlat+0.5_ireals*dlat*degrad)
    zsinlat(j)   = SIN(zlat)
  ENDDO

  ! zwk(k) bestimmen
  zwk(:) = 1.0_ireals

  zdlonr = 1.0_ireals / (dlon * degrad)
  zdlatr = 1.0_ireals / (dlat * degrad)

  zdpdtlm (:,:) = 0.0_ireals
  zdpdtlm1(:,:) = 0.0_ireals
  zchi_new(:,:) = 0.0_ireals
  zchi    (:,:) = 0.0
  zdelpq  (:,:) = 0.0_ireals
  zud     (:,:) = 0.0
  zvd     (:,:) = 0.0

!------------------------------------------------------------------------------
! Section 2: Compute surface pressure tendency zdpdtlm
!------------------------------------------------------------------------------

  DO k = 1, kevert
    ! Compute layer thickness delp
    zdelp (:,:) = (akhalf(k+1)-akhalf(k)) + (bkhalf(k+1)-bkhalf(k))*ps_lm(:,:)
    zdelpq(:,:) = zdelpq(:,:) + zwk(k)*zdelp(:,:)

    DO j = 2, je2lm - 1
      DO i = 2, ie2lm - 1
        zul = u_lm(i-1,j  ,k) * 0.5 * (zdelp(i-1,j  ) + zdelp(i,j))
        zur = u_lm(i  ,j  ,k) * 0.5 * (zdelp(i+1,j  ) + zdelp(i,j))
        zvo = v_lm(i  ,j  ,k) * 0.5 * (zdelp(i  ,j+1) + zdelp(i,j))
        zvu = v_lm(i  ,j-1,k) * 0.5 * (zdelp(i  ,j-1) + zdelp(i,j))
        zdpdtlm(i,j) = zdpdtlm(i,j) - 1.0 /                                  &
          (r_earth*0.5*(zcoslat_v(j) + zcoslat_v(j-1))) * ((zur-zul)*zdlonr &
                            + (zvo*zcoslat_v(j) - zvu*zcoslat_v(j-1))*zdlatr)
      ENDDO
    ENDDO
  ENDDO

  IF (num_compute > 1) THEN
    ! Exchange one line of zdpdtlm
    kzdims=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                              &
       (1, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
        ie2lm, je2lm, kzdims, jstartpar, jendpar, 1, nboundlines,      &
        my_cart_neigh, .FALSE., .FALSE., .FALSE.,                      &
        113, ldatatypes, ncomm_type, izerror, yzerrmsg,                &
        zdpdtlm)
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Solution of the elliptic equation for the potential chi
!------------------------------------------------------------------------------

  ! rhs: right hand side
  DO j = 1, je2lm
    DO i = 1, ie2lm
      zrhs(i,j) = (zdpdtlm(i,j) - dpsdt_gl(i,j))*(r_earth*zcoslat_m(j))**2
    ENDDO
  ENDDO

  DO j = 2, je2lm - 1
    DO i = 2, ie2lm - 1
      zal (i,j) = zdlonr**2 *                                                &
                      (zdelpq(i,j) + 0.25*(zdelpq(i+1,j) - zdelpq(i-1,j)) )
      zcl (i,j) = zdlonr**2 *                                                &
                      (zdelpq(i,j) - 0.25*(zdelpq(i+1,j) - zdelpq(i-1,j)) )
      zalp(i,j) = zcoslat_m(j) * zdlatr * zdelpq(i,j) *                      &
                      (zcoslat_m(j) * zdlatr - 0.5 * zsinlat(j))             &
                   +   0.25 * (zcoslat_m(j) * zdlatr)**2 *                   &
                      (zdelpq(i,j+1) - zdelpq(i,j-1))
      zclp(i,j) = zcoslat_m(j) * zdlatr * zdelpq(i,j) *                      &
                      (zcoslat_m(j)*zdlatr + 0.5*zsinlat(j))                 &
                   -   0.25 * (zcoslat_m(j) * zdlatr)**2 *                   &
                      (zdelpq(i,j+1) - zdelpq(i,j-1))
      zbb (i,j) = - 2.0 * zdelpq(i,j) * (zdlonr**2 + (zdlatr*zcoslat_m(j))**2)
    ENDDO
  ENDDO

  ! The iteration scheme for the parallel program has changed because of the
  ! data dependencies and to avoid differences when running on a different
  ! number of processors. Therefore the number alsor has been reduced from
  ! alsor=1.70 to alsor=1.0, because the scheme would not converge otherwise.
  ! The number of maximal iterations has been set to 200 (instead of 20).
  ! (The difference is that in the sequential program chi has been overwritten
  !  at once, hence for chi(i-1,j) and chi(i,j-1) the new values were used.
  !  Now, chi_new is calculated first on the whole domain and then chi is
  !  overwritten.)

  ! SOR correction step
  zalsor = 1.0_ireals

  DO iziter = 1,200 ! RJ: Allow 200 iterations instead of 20
    zdmax = 0.0_ireals
    DO j = 2,je2lm-1
      DO i = 2,ie2lm-1
        zres  = zal (i,j) * zchi(i+1,j  ) + zcl (i,j) * zchi(i-1,j  ) + &
                zalp(i,j) * zchi(i  ,j+1) + zclp(i,j) * zchi(i  ,j-1) + &
                zbb (i,j) * zchi(i  ,j  ) - zrhs(i,j)
        zdchi = zalsor * zres / zbb(i,j)
        zdmax = MAX ( zdchi, ABS(zdmax) )
        zchi_new(i,j) = zchi(i,j) - zdchi
      ENDDO
    ENDDO

    zchi(:,:) = zchi_new(:,:)
    IF (num_compute > 1) THEN
      ! Exchange one line of zchi
      itag = 120 + iziter
      kzdims=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                              &
         (1, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie2lm, je2lm, kzdims, jstartpar, jendpar, 1, nboundlines,      &
          my_cart_neigh, .FALSE., .FALSE., .FALSE.,                      &
          itag, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
          zchi)
      CALL global_values (zdmax, 1, 'MAX', imp_reals, icomm_cart, -1,    &
                          yzerrmsg, izerror)
    ENDIF

    IF (iziter == 1) THEN
      zdmax1 = zdmax
      IF (my_cart_id == 0) THEN
        WRITE (noutput,'(A)') '     '
        WRITE (noutput,'(A)') ' SOR-Solution for Wind Potential'
      ENDIF
    ELSE
      IF (my_cart_id == 0) THEN
        WRITE (noutput,'(A,I5,A,F10.4)')                                    &
                               ' iteration: ',iziter,';   zdmax = ',zdmax
      ENDIF
      IF (zdmax / zdmax1 < 0.05) THEN
        EXIT   ! loop iter
      ENDIF
    ENDIF
  ENDDO

  IF (my_cart_id == 0) THEN
    WRITE (noutput,'(A)') '     '
    WRITE (noutput,'(A)') '     '
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Compute the wind corrections zud, zvd and correct the winds
!------------------------------------------------------------------------------

  DO j = 2, je2lm-1
    DO i = 1, ie2lm-1
      zud(i,j) = 1.0/ (r_earth*zcoslat_m(j)) * (zchi(i+1,j)-zchi(i,j))*zdlonr
    ENDDO
  ENDDO

  DO j = 1, je2lm-1
    DO i = 2, ie2lm-1
      zvd(i,j) = 1.0/   r_earth    * (zchi(i,j+1) - zchi(i,j))*zdlatr
    ENDDO
  ENDDO

  IF (num_compute > 1) THEN
    ! Exchange one line of zud and zvd
    kzdims=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                              &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
        ie2lm, je2lm, kzdims, jstartpar, jendpar, 1, nboundlines,      &
        my_cart_neigh, .FALSE., .FALSE., .FALSE.,                      &
        114, .FALSE., ncomm_type, izerror, yzerrmsg,                   &
        zud, zvd)
  ENDIF

  ! Correct LM wind
  DO k = 1,kevert
    u_lm(:,:,k) = u_lm(:,:,k) + zwk(k) * zud(:,:)
    v_lm(:,:,k) = v_lm(:,:,k) + zwk(k) * zvd(:,:)
  END DO

!------------------------------------------------------------------------------
! Section 5: Compute corrected surface pressure tendency
!------------------------------------------------------------------------------

  DO k = 1, kevert
    ! Compute layer thickness delp
    zdelp (:,:) = (akhalf(k+1)-akhalf(k)) + (bkhalf(k+1)-bkhalf(k))*ps_lm(:,:)

    DO j = 2, je2lm - 1
      DO i = 2, ie2lm - 1
        zul = u_lm(i-1,j  ,k) * 0.5 * (zdelp(i-1,j  ) + zdelp(i,j))
        zur = u_lm(i  ,j  ,k) * 0.5 * (zdelp(i+1,j  ) + zdelp(i,j))
        zvo = v_lm(i  ,j  ,k) * 0.5 * (zdelp(i  ,j+1) + zdelp(i,j))
        zvu = v_lm(i  ,j-1,k) * 0.5 * (zdelp(i  ,j-1) + zdelp(i,j))
        zdpdtlm1(i,j) = zdpdtlm1(i,j) - 1.0 /                               &
          (r_earth*0.5*(zcoslat_v(j) + zcoslat_v(j-1))) * ((zur-zul)*zdlonr &
                           + (zvo*zcoslat_v(j) - zvu*zcoslat_v(j-1))*zdlatr)
      ENDDO
    ENDDO
  ENDDO

  IF (num_compute > 1) THEN
    ! Exchange one line of zdpdtlm1
    kzdims=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                              &
       (1, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
        ie2lm, je2lm, kzdims, jstartpar, jendpar, 1, nboundlines,      &
        my_cart_neigh, .FALSE., .FALSE., .FALSE.,                      &
        115, ldatatypes, ncomm_type, izerror, yzerrmsg,                &
        zdpdtlm1)
  ENDIF

!------------------------------------------------------------------------------
! Section 6: Control outputs
!------------------------------------------------------------------------------

  ! Control output of meanvalue, minimum and maximum of zud
  yzstring = '  Statistics for absolute wind correction: ud-pot.'
  CALL lm_field_stats (ABS(zud), yzstring, 1.0_ireals,noutput,yzerrmsg,izerror)

  ! Control output of meanvalue, minimum and maximum of zvd
  yzstring = '  Statistics for absolute wind correction: vd-pot.'
  CALL lm_field_stats (ABS(zvd), yzstring, 1.0_ireals,noutput,yzerrmsg,izerror)

  ! Print zud, zvd and dpdt
  IF (lprud) THEN
    CALL print_lm_field (zud(:,:), 'ud-pot', 'm/s', kevert+1, 100.0_ireals,  &
      0.0_ireals, noutput, 2, ie2lm_tot-1, 2, je2lm_tot-1, yzerrmsg, izerror)
  ENDIF
  IF (lprvd) THEN
    CALL print_lm_field (zvd(:,:), 'vd-pot', 'm/s', kevert+1, 100.0_ireals,  &
      0.0_ireals, noutput, 2, ie2lm_tot-1, 2, je2lm_tot-1, yzerrmsg, izerror)
  ENDIF

  ! Print surface pressure tendencies if wanted
  IF (lprdpdt) THEN
    CALL print_lm_field (dpsdt_gl, 'dpsdt_gl', 'Pa/s',kevert+1,3600.0_ireals,&
      0.0_ireals, noutput, 2, ie2lm_tot-1, 2, je2lm_tot-1, yzerrmsg, izerror)

    CALL print_lm_field (zdpdtlm, 'dpsdt_gl', 'Pa/s', kevert+1,3600.0_ireals,&
      0.0_ireals, noutput, 2, ie2lm_tot-1, 2, je2lm_tot-1, yzerrmsg, izerror)

    CALL print_lm_field (zdpdtlm1, 'dpsdt_gl', 'Pa/s',kevert+1,3600.0_ireals,&
      0.0_ireals, noutput, 2, ie2lm_tot-1, 2, je2lm_tot-1, yzerrmsg, izerror)

    zrhs(:,:) = dpsdt_gl(:,:)-zdpdtlm(:,:)
    CALL print_lm_field (zrhs(:,:), 'dpsdt_gl','Pa/s',kevert+1,3600.0_ireals,&
      0.0_ireals, noutput, 2, ie2lm_tot-1, 2, je2lm_tot-1, yzerrmsg, izerror)
  ENDIF

! Keep ud and vd at some gridpoints
! IF (lprgp ) THEN
!   DO l = 1, ngp
!     udpr(l) = ud(igp(l),jgp(l))
!     vdpr(l) = vd(igp(l),jgp(l))
!   ENDDO
! ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE uv_correction

!==============================================================================

!+ Extrapolation of the surface pressure to sea level.
!------------------------------------------------------------------------------

SUBROUTINE ps_to_sealevel (kevert, akhalf, bkhalf, kgr)

!------------------------------------------------------------------------------
!
! Description: Extrapolation of the surface pressure to sea level.
!
! Method:
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(IN)    ::   &
  kevert,     & ! vertical dimension of coordinate parameters
  kgr           ! top of boundary layer

REAL (KIND=ireals),       INTENT (IN)        ::  &
  akhalf(kevert+1), bkhalf(kevert+1)   ! coordinate parameters for half levels

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers) ::   &
  i, j, k, izerror

REAL (KIND=ireals)       ::   &
  zdfis, zpgr

! Local arrays:
REAL (KIND=ireals)                      ::   &
  ztv   (ie2lm,je2lm,kevert),    & ! virtual temperature on GME-levels
  zpn   (ie2lm,je2lm,kevert+1),  & ! pressure on GME main-levels
  zpsred(ie2lm,je2lm),           & ! reduced surface pressure
  zbtv  (ie2lm,je2lm),           & !
  zctv  (ie2lm,je2lm),           & !
  zpq   (ie2lm,je2lm),           & !
  zpq2  (ie2lm,je2lm),           & !
  ztvq  (ie2lm,je2lm),           & !
  ztvpq (ie2lm,je2lm),           & !
  zdtv  (ie2lm,je2lm),           & !
  zfiu  (ie2lm,je2lm),           & !
  zfio  (ie2lm,je2lm)              !

CHARACTER (LEN=80)                         ::                              &
  yzerrmsg               ! for error message

CHARACTER (LEN=25)                         ::                              &
  yzroutine              ! character variable to be passed to remark

!
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  izerror   = 0
  yzerrmsg  = '   '
  yzroutine = 'ps_to_sealevel'

  ! Compute virtual temperature on HM/LM main levels
  DO k = 1, kevert
    ztv(:,:,k) = t_lm(:,:,k) * (1.0 + Rvd_m_o * qv_lm(:,:,k) - qc_lm(:,:,k))
  ENDDO

  ! Compute the pressure on the half levels (pn)
  DO k = 1, kevert+1
    zpn(:,:,k) = akhalf(k) + bkhalf(k)*ps_lm(:,:)
  ENDDO

!------------------------------------------------------------------------------
! Section 2: Linear approximation of the virtual temperature above kgr
!------------------------------------------------------------------------------

  ! Linear approximation of the virtual temperature above the boundary layer.
  ! Regression of 4 layers above k = kgr.
  ! The linear approximation is of the form:
  !  ztv(p) = zbtv*p + zctv  ,  where p is the pressure at the full levels.
  !
  ! The virtual temperature gradient (dtvdz) is computed from zbtv and from
  ! the temperature change we obtain the geopotential difference 0 - fis_lm

  zpgr = 850.0E2_ireals
  zpq  (:,:) = 0.0
  zpq2 (:,:) = 0.0
  ztvq (:,:) = 0.0
  ztvpq(:,:) = 0.0

  DO k = kgr-3, kgr
    zpq (:,:)  = zpq  (:,:) + 0.5 * (zpn(:,:,k+1) + zpn(:,:,k))
    zpq2(:,:)  = zpq2 (:,:) + 0.5 * (zpn(:,:,k+1) + zpn(:,:,k))* &
                              0.5 * (zpn(:,:,k+1) + zpn(:,:,k))
    ztvq (:,:) = ztvq (:,:) + ztv(:,:,k)
    ztvpq(:,:) = ztvpq(:,:) + ztv(:,:,k) * 0.5 * (zpn(:,:,k+1) + zpn(:,:,k))
  ENDDO

  zbtv (:,:) = (ztvpq(:,:) - ztvq(:,:)*zpq(:,:)/4.0) /                  &
                                       (zpq2(:,:) - zpq(:,:)**2 / 4.0)
  zctv (:,:) = (ztvq (:,:) - zbtv(:,:)*zpq(:,:)) / 4.0
  zdtv (:,:) =  zbtv (:,:) * zpgr/(R_d*(zbtv(:,:)*zpgr + zctv(:,:)))    &
                                                        * fis_lm(:,:)

!------------------------------------------------------------------------------
! Section 3: Compute pressure reduction
!------------------------------------------------------------------------------

  ! Compute the average temperature (ztvq) of the layer (0 - fis_lm)
  zfiu(:,:) = 0.0
! WHERE (fis_lm(:,:) > 0.0)
!   zfiu (:,:) = fis_lm(:,:)
! END WHERE

  DO k = kevert, 2, -1
    DO j = 1, je2lm
      DO i = 1, ie2lm
        zfio(i,j) = zfiu(i,j) + R_d * ztv(i,j,k) * LOG(zpn(i,j,k+1)/zpn(i,j,k))

        ! Compute the reduced pressure
        IF (fis_lm(i,j) >= zfiu(i,j) .AND. fis_lm(i,j) <= zfio(i,j)) THEN
          zdfis = (fis_lm(i,j) - zfiu(i,j))/(zfio(i,j) - zfiu(i,j)) * &
                       (zpn(i,j,k) - zpn(i,j,k+1)) + zpn(i,j,k+1)
          ztvpq (i,j) = ztv(i,j,k) * (zpn(i,j,k+1) - zdfis)      + ztvpq(i,j)
          zpq   (i,j) =               zpn(i,j,k+1) - zdfis       + zpq  (i,j)
          ztvq  (i,j) = ztvpq(i,j) / zpq(i,j)
          zpsred(i,j) = ps_lm(i,j)*EXP(fis_lm(i,j)/(R_d*(ztvq(i,j) + zdtv(i,j))))
        ELSEIF (fis_lm(i,j) <= 0.0) THEN
   !!!  ELSEIF (fis_lm(i,j) == 0.0) THEN   ???
          zpsred(i,j) = ps_lm(i,j)
        ELSE
          ztvpq (i,j) = ztv(i,j,k)*(zpn(i,j,k+1) - zpn(i,j,k)) + ztvpq(i,j)
          zpq   (i,j) =             zpn(i,j,k+1) - zpn(i,j,k)  + zpq  (i,j)
        ENDIF
      ENDDO
    ENDDO
    zfiu(:,:) = zfio(:,:)
  ENDDO

!------------------------------------------------------------------------------
! Section 4: Control output
!------------------------------------------------------------------------------

  CALL print_lm_field (zpsred,'psred','Pa ', kevert, 0.01_ireals,          &
      0.0_ireals, noutput, 1, ie2lm_tot, 1, je2lm_tot, yzerrmsg, izerror)

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE ps_to_sealevel

!==============================================================================

END MODULE src_vert_interpol
