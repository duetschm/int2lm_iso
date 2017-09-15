!+ Source Module for computing LM fields
!------------------------------------------------------------------------------

MODULE src_vert_inter_lm

!------------------------------------------------------------------------------
!
! Description:
!   This module contains routines for the vertical interpolation of LM levels.
!
! Current Code Owner: DWD, Ulrich Schaettler
!    phone:  +49  69  8062 2739
!    fax:    +49  69  8062 3721
!    email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2005/04/11 Guy de Morsier, MCH
!  Initial release for INT2LM
! V1_5         2007/07/09 Ulrich Schaettler
!  Added interpolation of qr,qs,qg
!  Eliminate allocation of ps_lm, grh_lm (is done in int2lm_org)
!  Added treatment of chemistry variables
! V1_7         2007/11/26 Lucio Torrisi
!  Additional splitting of relative humidity after balancing and filtering pp
! V1_8         2008/05/29 Ulrich Schaettler
!  Vectorization of vertical interpolation and necessary pre- and postprocessing
!  Eliminated ldebug
! V1_9         2009/09/03 Guenther Zaengl, Anne Roches
!  Improved initialization of perturbation pressure (Guenther Zaengl)
!  Bugfix in calculating the linear regression in vert_interp
!  Do no linear regression any more for u, v (only for t)
!  Use a value 50 m below the surface for every grid points as lower boundary
! V1_10        2009/12/17 Ulrich Schaettler
!  Adaptations to process Unified Model data
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_12        2010/06/14 Ulrich Schaettler
!  Implemented possibility to switch the reference atmosphere between incoming
!  and outgoing data
! V1_19        2012/06/06 Ulrich Schaettler, CLM
!  Renamed ak_in_uv, bk_in_uv to ak_in_rho, bk_in_rho according to UM conventions
!  Also bug fix for tautsp2D/izln_vec as in src_lm_fields
!  Some extensions for processing UM data
!  Adaptations to environment.f90, because of unification with COSMO-Model 4.23
!  CLM:
!    Several changes to interpolate hybrid height coordinates in case of lcm_hgt_coor=.TRUE.
!    Correction of calculations in the lowest model layer for lcm_hgt_coor=.TRUE.
!    Correction of if clause for lcm_hgt_coor=.TRUE. and adition of "_ireals"
!    Correction of if clauses by extension of 'pp'
! V1_22        2013/07/11 Ulrich Schaettler, KIT
!  Renamed lprog_qrqs to lprog_qr_qs to be consistent with other names
!  Use structures refatm_out, vcoord_out from vgrid_refatm_utils
!    BUG FIX: wrong vertical coordinate parameters were used to compute boundary
!      layer height for ivctype=2 and llm2lm (height instead of pressure).
!      Now always the pressure coordinates vcoord_in%sigm_coord are used.
!  Renamed grib buffers: ds_grib to ds_grib_single, ds_gribapi to ds_grib_double
!  Replaced l_chemistry by l_art, l_art_nested (KIT)
! V1_23        2013/10/02 Ulrich Schaettler
!  Initialize the boundary layer height with 0 (to detect errors early)
!  Rename vcoord_out, refatm_out to vcoord, refatm (as is in COSMO-Model)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:

USE data_parameters , ONLY :  &

! Imported Parameters:
  ireals,    &! KIND-type parameters for real variables
  iintegers   ! KIND-type parameter for "normal" integer variables

!------------------------------------------------------------------------------

USE data_fields_lm, ONLY : &
  ps_lm     ,      & ! surface pressure                                 ( Pa  )
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
  p_lm      ,      & ! full pressure (needed for lum2lm)                ( Pa  )
  pp_lm     ,      & ! deviation from the reference pressure            ( Pa  )
  p0_lm     ,      & ! reference pressure                               ( Pa  )
  p0_gl     ,      & ! ref. pres. on full levels + interpol. COARSE LM oro.(Pa)
  t0_lm     ,      & ! reference temperature                            (  K  )
  rho0_lm   ,      & ! reference density                                (kg/m^3)
  t_s_gl    ,      & ! temperature of the ground surface                (  K  )
  dp0_lm    ,      & ! reference pressure thickness of layers           ( Pa  )
  hhl_gl    ,      & ! height of half-levels on the interpol. COARSE LM oro.(m)
  hhl_lm    ,      & ! height of half-levels of LM                      (  m  )
  hsurf_gl           ! height of orography interpolated from coarse grid(  m  )

!------------------------------------------------------------------------------

USE data_grid_lm,   ONLY : &
  dlat,        & ! grid point distance in zonal direction (in degrees)
  dlon,        & ! grid point distance in meridional direction (in degrees)
  startlat,    & ! transformed latitude of the lower left grid point
                 ! of the local domain (in degrees, N>0)
  ie2lm,       & !
  je2lm,       & !
  kelm,        & !
  kedim,       & !
  ke1lm,       & !
  jstartpar,   & ! start index for computations in the parallel program
  jendpar        ! end index for computations in the parallel program

!------------------------------------------------------------------------------

USE data_grid_in,   ONLY : &
  akh_in,      & ! vertical coordinate parameters for main levels
  bkh_in,      & !                  - " -
  akh_in_rho,  & ! coefficients for main levels for u, v (lum2lm)
  bkh_in_rho,  & ! coefficients for main levels for u, v (lum2lm)
  ke1in,       & ! ke1 for input grid
  ke_in,       & ! ke for input grid
  lcm_hgt_coor,& ! Input data has hybrid height coordinates
  klv850_in      ! approximate level where 850 hPa is reached

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
    B4_w       !  d

!------------------------------------------------------------------------------

USE data_int2lm_control,       ONLY :  &
    ndebug,       & ! unit number for file with debug information
    llm2lm,       & ! if .TRUE., lm ->lm
    lum2lm,       & ! if .TRUE., um ->lm
    lcm2lm,       & ! if .TRUE., cm ->lm
    lcomp_bound,  & ! compute fields for boundaries
    lvertwind_ini,& ! if .TRUE., compute vertical wind for LM for initial data
    lvertwind_bd, & ! if .TRUE., compute vertical wind for LM for boundary data
    lprog_qi,     & ! if .TRUE., interpolate qi to LM grid
    lprog_qr_qs,  & ! if .TRUE., interpolate qr,qs to LM grid
    lprog_qg,     & ! if .TRUE., interpolate qg to LM grid
    lfilter_pp,   & ! if .TRUE., filter the pressure deviation after vertical
                    !            interpolation
    lbalance_pp,  & ! if .TRUE., compute a hydrostatic balanced pp after
                    !            vertical interpolation in LM2LM
    l_art,        & ! if .TRUE., interpolate additional int2lm art fields
    l_art_nested, & ! if .TRUE., interpolate additional lm2lm art fields
    idbg_level,   & ! to control verbosity of output
    lprintdeb_all,& ! whether all PEs print debug output
    qvmin,        & ! minimum value of water vapor (security)
    qcmin,        & ! minimum value of cloud water (security)
    qimin           ! minimum value of cloud ice content (security)

!------------------------------------------------------------------------------

USE data_int2lm_io,            ONLY :  &
  nvar_lm_norm,      & ! maximum number of variables in LM variable table
  nvar_lm_chem,      & ! maximum number of variables in LM variable table
  nvar_in_norm,      & ! maximum number of variables in input variable table
  nvar_in_chem,      & ! maximum number of variables in input variable table
  var_in,            & ! variable table for input fields
  var_lm               ! variable table for LM fields

!------------------------------------------------------------------------------

!USE data_profiles,    ONLY :   &
!    lprgp,        & ! logical for print at selected grid points
!    ngp_tot,      & ! number of selected grid points in the total domain
!    igp_tot,      & ! i-indices    total domain
!    jgp_tot         ! j-indices    total domain

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
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

USE utilities,          ONLY :   tautsp2D, horizontal_filtering

!------------------------------------------------------------------------------

USE environment,        ONLY :   model_abort, exchg_boundaries, extend_field

!------------------------------------------------------------------------------

USE parallel_utilities, ONLY :   remark, i_global, j_global, global_values

!------------------------------------------------------------------------------

USE vgrid_refatm_utils, ONLY :   refatm, refatm_in, vcoord_in

!==============================================================================

IMPLICIT NONE

INTEGER (KIND=iintegers) ::  &
  hfwidth, hfw_m_nb, ie2lm_hf, je2lm_hf


!==============================================================================

CONTAINS

!==============================================================================
!+ Organizes the interpolation of LM-fields
!------------------------------------------------------------------------------

SUBROUTINE org_vert_inter_lm

!------------------------------------------------------------------------------
!
! Description:
!   org_vert_inter_lm organizes the interpolation of LM-fields that are
!   necessary for initial or boundary files for the nonhydrostatic LM.
!   These computations include:
!     - for the variables qc_lm and qv_lm, the generalized relative humidity
!       is interpolated
!     - special vertical interpolation for pressure pp_lm
!     - vertical interpolation to LM levels
!     - pressure deviation and splitting grh_lm in qc_lm and qv_lm
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Local arrays
REAL (KIND=ireals)         ::  &
  rho0  (ie2lm,je2lm),         & ! reference density for the surface layer
  qrs   (ie2lm,je2lm),         & ! precipitation water (water loading)
  zhi_hl(ie2lm,je2lm,ke_in+1), & ! height on COARSE half levels
  zhi_fl(ie2lm,je2lm,ke_in  ), & ! height on COARSE full levels
  zhi_fl_u(ie2lm,je2lm,ke_in), & ! height on COARSE full levels for u-wind (lum2lm)
  zhi_fl_v(ie2lm,je2lm,ke_in), & ! height on COARSE full rho levels for v (lum2lm)
  zhi_fl_p(ie2lm,je2lm,ke_in)    ! height on COARSE full rho levels for p (lum2lm)

! Local allocatable arrays
REAL (KIND=ireals), ALLOCATABLE ::    &
  pp_tmp(:,:,:)

INTEGER  (KIND=iintegers)  ::  &
  istata, istatd, izerror,     & ! status and error status variable
  zpgr, kloc,                  & ! top of boundary layer for input levels
  i, j, k, kup, klow,          &
  izdebug, n, nfilt, hfjstartpar, hfjendpar, kzdims(24), nzentry, iztable,   &
  iterate

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling

! Definition of statement functions
REAL    (KIND=ireals)    ::    &
  zpm, zaq, zbq, sf_psat_w, sf_qsat, x, y, z, v, w, &
  ztvw, zphf, zfakt, zpa, zleft, zrhs, zt0, ztvw1, zgdrt, ztdbe, zbetf, zt00

REAL    (KIND=ireals)    ::    &
    zp0_lm_refmod(ie2lm,je2lm,ke1lm), &
    ztvdt   (ie2lm,je2lm,2),   &
    zt0dp0t (ie2lm,je2lm,2),   &
    zhsurf_u(ie2lm,je2lm),     &
    zhsurf_v(ie2lm,je2lm),     &
    zplm    (ie2lm,je2lm,ke1lm)
    
sf_psat_w  (x,y,z,v,w) = y * EXP(z*(x-v)/(x-w))
!!sf_qsat    (x,y,z,v)   = z * x / (y-v*x)

sf_qsat    (x,y,z,v)   = z * x / MAX( (y-v*x), 1.0_ireals)
!!sf_qsat = Rdv*psatx/MAX((px-O_m_rdv*psatx),1.0_ireals)

!
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 0: Initializations
!------------------------------------------------------------------------------

  izerror   = 0_iintegers
  kloc      = 0_iintegers
  yzerrmsg  = '  '
  yzroutine = 'org_vert_inter_lm'

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
! Section 1: Compute generalized relative humidity on horiz. interp. grid
!------------------------------------------------------------------------------

  IF     (llm2lm) THEN
    DO k = 1, ke_in
      DO j = 1, je2lm
        DO i = 1, ie2lm
          zpm = p0_gl(i,j,k) + pp_lm(i,j,k)
          zaq    = sf_psat_w (t_lm(i,j,k), B1, B2_w, B3, B4_w)
          zbq    = sf_qsat   (zaq, zpm, Rdv, O_m_rdv)
          grh_lm(i,j,k) = (qv_lm(i,j,k) + qc_lm(i,j,k)) / zbq
        ENDDO
      ENDDO
    ENDDO
  ELSEIF (lum2lm .OR. (lcm2lm .AND. lcm_hgt_coor)) THEN
    DO k = 1, ke_in
      DO j = 1, je2lm
        DO i = 1, ie2lm
          ! use full pressure horizontally interpolated to p_lm
          zaq    = sf_psat_w (t_lm(i,j,k), B1, B2_w, B3, B4_w)
          zbq    = sf_qsat   (zaq, p_lm(i,j,k), Rdv, O_m_rdv)
          grh_lm(i,j,k) = (qv_lm(i,j,k) + qc_lm(i,j,k)) / zbq
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Compute boundary layer height and interpolate pressure
!------------------------------------------------------------------------------

  ! Find the boundary layer height (kzgr)
  IF     (llm2lm) THEN
! NOTE: the below computation is only valid, if the input vertical coordinate
!       is pressure based (which it is not if ivctype != 1)
!!$    zpgr = 850.0E2
!!$    DO k = ke_in, 1, -1
!!$      zpm = akh_in(k) + bkh_in(k)*1.0E5_ireals
!!$      IF (zpm > zpgr) kz = k
!!$    ENDDO
!!$    kz = kz - 1
    kloc = klv850_in
  ELSEIF (lum2lm) THEN
    ! this is just a rough guess
    kloc = 20
  ELSEIF (lcm2lm .AND. lcm_hgt_coor) THEN
    kloc = 30
  ENDIF

  IF (my_cart_id == 0) THEN
    PRINT *, yzroutine, ": Boundary layer height in COARSE layer # :", kloc
  ENDIF

  IF     (llm2lm) THEN

    ! Compute height on full LM COARSE levels
    DO k = 1, ke_in
      zhi_fl(:,:,k) = 0.5_ireals * hhl_gl (:,:,k+1) + 0.5_ireals * hhl_gl(:,:,k)
    ENDDO

  ELSEIF (lum2lm) THEN

    ! Compute orography on u-grid points
    DO j = 1, je2lm
      DO i = 1, ie2lm-1
        zhsurf_u (i,j)   = 0.5_ireals * (hsurf_gl(i,j) + hsurf_gl(i+1,j))
      ENDDO
      zhsurf_u(ie2lm,j)  = hsurf_gl(ie2lm,j)
    ENDDO

    ! Compute orography on v-grid points
    DO j = 1, je2lm-1
      DO i = 1, ie2lm
        zhsurf_v (i,j)   = 0.5_ireals * (hsurf_gl(i,j) + hsurf_gl(i,j+1))
      ENDDO
    ENDDO
    DO i = 1, ie2lm
      zhsurf_v(i,je2lm)  = hsurf_gl(i,je2lm)
    ENDDO

    ! Compute height on UM COARSE levels (mass variables)
    ! this formula has been taken from the spanish routine "geopotencial_um"
    ! but we do compute the height, not the geopotential
    DO k = 1, ke_in
      zhi_fl   (:,:,k) = akh_in    (k) + bkh_in    (k) * hsurf_gl(:,:)  ! for T, Qx
      zhi_fl_u (:,:,k) = akh_in_rho(k) + bkh_in_rho(k) * zhsurf_u(:,:)  ! for U
      zhi_fl_v (:,:,k) = akh_in_rho(k) + bkh_in_rho(k) * zhsurf_v(:,:)  ! for V
      zhi_fl_p (:,:,k) = akh_in_rho(k) + bkh_in_rho(k) * hsurf_gl(:,:)  ! for P
    ENDDO

    ! Compute pressure deviation from (incoming) reference pressure on
    ! horizontally interpolated levels (of the incoming vertical grid)
    DO k = 1, ke_in
      pp_lm(:,:,k) = p_lm(:,:,k) - p0_gl(:,:,k)
    ENDDO

  ELSEIF (lcm2lm .AND. lcm_hgt_coor) THEN

    ! Here, all variables are defined at the same grid point
    ! and also P is defined on the same levels???
    DO k = 1, ke_in
      zhi_fl   (:,:,k) = akh_in    (k) + bkh_in    (k) * hsurf_gl(:,:)
      zhi_fl_u (:,:,k) = akh_in_rho(k) + bkh_in_rho(k) * hsurf_gl(:,:)
      zhi_fl_v (:,:,k) = akh_in_rho(k) + bkh_in_rho(k) * hsurf_gl(:,:)
    ENDDO

    ! The vertical interpolation will not be done with pressure deviation,
    ! but with LOG(p)
    DO k = 1, ke_in
      pp_lm(:,:,k) = LOG(p_lm(:,:,k))
    ENDDO

  ENDIF

  IF     (lum2lm) THEN
    ! UM pressure data are on rho-levels
    CALL vert_interp (pp_lm, 'pp',   kedim, zhi_fl_p, ke_in, kloc, izdebug)
  ELSEIF (lcm2lm .AND. lcm_hgt_coor) THEN
    ! HadGEM2 data are on usual model levels
    CALL vert_interp (pp_lm, 'logp', kedim, zhi_fl,   ke_in, ke_in,izdebug)

    ! Compute full pressure and pressure deviation
    DO k = 1, kedim
       p_lm(:,:,k) = EXP(pp_lm(:,:,k))
      pp_lm(:,:,k) = p_lm(:,:,k) - p0_lm(:,:,k)
    ENDDO
  ELSE
    CALL vert_interp (pp_lm, 'pp',   kedim, zhi_fl,   ke_in, kloc, izdebug)
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Compute or interpolate the vertical velocity
!------------------------------------------------------------------------------

  IF (llm2lm) THEN
    zhi_hl(:,:,:) = hhl_gl(:,:,:)

    IF ( lvertwind_ini .AND. .NOT. lcomp_bound ) THEN
      IF (izdebug > 5) THEN
        PRINT *, 'The vertical wind is interpolated from input'
      ENDIF
      CALL vert_interp (w_lm, 'w', kedim+1, zhi_hl, ke_in+1, kloc, izdebug)
    ENDIF

    IF ( lvertwind_bd  .AND.       lcomp_bound ) THEN
      IF (izdebug > 5) THEN
        PRINT *, 'The vertical wind is interpolated from input'
      ENDIF
      CALL vert_interp (w_lm, 'w', kedim+1, zhi_hl, ke_in+1, kloc, izdebug)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Vertical interpolations for the other variables
!------------------------------------------------------------------------------

  IF     (llm2lm) THEN
    CALL vert_interp (u_lm  , 'u'  , kedim, zhi_fl, ke_in, kloc, izdebug)
    CALL vert_interp (v_lm  , 'v'  , kedim, zhi_fl, ke_in, kloc, izdebug)
  ELSEIF (lum2lm .OR. (lcm2lm .AND. lcm_hgt_coor)) THEN
    CALL vert_interp (u_lm  , 'u'  , kedim, zhi_fl_u, ke_in, kloc, izdebug)
    CALL vert_interp (v_lm  , 'v'  , kedim, zhi_fl_v, ke_in, kloc, izdebug)
  ENDIF

  CALL vert_interp (t_lm  , 't'  , kedim, zhi_fl, ke_in, kloc, izdebug)
  CALL vert_interp (grh_lm, 'grh', kedim, zhi_fl, ke_in, kloc, izdebug)

  IF (lprog_qi) THEN
    CALL vert_interp (qi_lm, 'qi', kedim, zhi_fl, ke_in, kloc, izdebug)
  ENDIF

  IF (llm2lm) THEN
    IF (lprog_qr_qs) THEN
      CALL vert_interp (qr_lm, 'qr', kedim, zhi_fl, ke_in, kloc, izdebug)
      CALL vert_interp (qs_lm, 'qs', kedim, zhi_fl, ke_in, kloc, izdebug)
    ENDIF
    IF (lprog_qg) THEN
      CALL vert_interp (qg_lm, 'qg', kedim, zhi_fl, ke_in, kloc, izdebug)
    ENDIF

    IF (l_art .AND. l_art_nested) THEN
      DO n = nvar_in_norm+1, nvar_in_chem
        IF (var_in(n)%lreadin) THEN
          ! Search for entry in LM variable table
          nzentry = -1
          DO iztable = nvar_lm_norm+1,nvar_lm_chem
            IF (TRIM(var_lm(iztable)%name) == TRIM(var_in(n)%name)) THEN
              nzentry = iztable
              EXIT
            ENDIF
          ENDDO
          IF (nzentry > 0) THEN
            CALL vert_interp (var_lm(nzentry)%p3, var_lm(nzentry)%name,     &
                              kedim, zhi_fl, ke_in, kloc, izdebug)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 5: Switch from incoming reference atmosphere to outgoing
!------------------------------------------------------------------------------

  IF (refatm_in%irefatm /= refatm%irefatm) THEN

    IF (my_cart_id == 0) THEN
      PRINT *, ' *** Switching from incoming reference atmosphere:  ', refatm_in%irefatm
      PRINT *, ' *** to outgoing reference atmosphere:              ', refatm%irefatm
    ENDIF

    ! Compute (incoming) reference atmosphere on outgoing vertical levels
    ! (the formulas for ivctype=2 from SR reference_atmosphere (for refatm_in%irefatm=1)
    !  or reference_atmosphere_2 (for refatm_in%irefatm=2) are used)

    IF     (refatm_in%irefatm == 1) THEN

      zgdrt = g/r_d/refatm_in%t0sl
      IF (refatm_in%dt0lp /= 0.0_ireals) THEN
        ztdbe = refatm_in%t0sl/refatm_in%dt0lp
      ELSE
        ztdbe = 0.0_ireals
      ENDIF
      zbetf = 2.0_ireals*refatm_in%dt0lp*zgdrt/refatm_in%t0sl

      DO k = 1, ke1lm
        IF (refatm_in%dt0lp == 0.0_ireals) THEN
          zp0_lm_refmod(:,:,k) = refatm_in%p0sl * EXP ( - zgdrt*hhl_lm(:,:,k) )
        ELSE
          zp0_lm_refmod(:,:,k) = refatm_in%p0sl * EXP ( - ztdbe*(1.0_ireals        &
                      - SQRT(1.0_ireals - zbetf*hhl_lm(:,:,k))) )
        ENDIF

        IF (k > 1) THEN
          ! averaging to main levels
          zp0_lm_refmod(:,:,k-1) = 0.5_ireals * (zp0_lm_refmod(:,:,k-1) + zp0_lm_refmod(:,:,k))
        ENDIF
      ENDDO

    ELSEIF (refatm_in%irefatm == 2) THEN

      zt00 = refatm_in%t0sl - refatm_in%delta_t

      DO k = 1, kelm
        zp0_lm_refmod(:,:,k) = refatm_in%p0sl*EXP ( - g/r_d*refatm_in%h_scal/zt00 * LOG(  &
                  (EXP( 0.5_ireals * ( hhl_lm(:,:,k) + hhl_lm(:,:,k+1) )    &
                      /refatm_in%h_scal)*zt00 + refatm_in%delta_t)/(zt00 + refatm_in%delta_t)) )
      ENDDO

    ENDIF

    ! Compute full pressure on outgoing vertical levels
    DO k = 1, kelm
      zplm(:,:,k) = zp0_lm_refmod(:,:,k) + pp_lm(:,:,k)
    ENDDO

    ! Compute pressure deviation from (outgoing) reference atmosphere pressure
    ! on outgoing vertical levels
    DO k = 1, kelm
      pp_lm(:,:,k) = zplm(:,:,k) - p0_lm(:,:,k)
    ENDDO
! ELSE

!   ! then only for UM data, the full pressure has to be calculated
!   ! p0_lm can be used for this

!   IF (lum2lm) THEN
!     ! compute full pressure on COSMO fine levels
!     p_lm(:,:,1:kelm) = p0_lm(:,:,1:kelm) + pp_lm(:,:,1:kelm)
!   ENDIF

  ENDIF

!------------------------------------------------------------------------------
! Section 6: Splitting of relative humidity
!------------------------------------------------------------------------------

  rho0(:,:) = rho0_lm(:,:,kelm)

  DO k = 1, kelm
    DO j = 1, je2lm
      DO i = 1, ie2lm
        ! split grh_lm into qv_lm, qc_lm using actual p = pp_lm + p0_lm
        zaq = sf_psat_w (t_lm(i,j,k), b1, b2_w, b3, b4_w)
        zpm = p0_lm(i,j,k) + pp_lm(i,j,k)
        zbq = sf_qsat   (zaq, zpm, Rdv, O_m_rdv)
!USUS   zbq = sf_qsat   (zaq, p_lm(i,j,k), Rdv, O_m_rdv)
        qv_lm(i,j,k) = MIN(1.0_ireals, grh_lm(i,j,k)) * zbq
        qc_lm(i,j,k) = MAX(0.0_ireals, grh_lm(i,j,k)-1.0_ireals) * zbq
        IF (qv_lm(i,j,k) < qvmin) qv_lm(i,j,k)  = qvmin
        IF (qc_lm(i,j,k) < qcmin) qc_lm(i,j,k)  = 0.0_ireals
      ENDDO
    ENDDO
  ENDDO

  qrs(:,:) = 0.0  !  no water loading
  ps_lm(:,:) = (p0_lm(:,:,kelm)+pp_lm(:,:,kelm)) * EXP( 0.5*dp0_lm(:,:,kelm) / &
     ( t_lm(:,:,kelm)*(1.0+rvd_m_o*qv_lm(:,:,kelm)-qc_lm(:,:,kelm)-qrs(:,:)) * &
     r_d*rho0(:,:) ) )

!------------------------------------------------------------------------------
! Section 7: Recalculate a hydrostatic balanced pressure deviation
!------------------------------------------------------------------------------

  IF (lbalance_pp) THEN
   DO iterate = 1,5
    ! organizational indices
    kup  = 2
    klow = 1

    ! pressure deviation on the lowest full level
    ! and other initializations
    DO j = 1,je2lm
      DO i = 1,ie2lm

        ztvw = t_lm(i,j,kelm) * (1.0+rvd_m_o*qv_lm(i,j,kelm)-qc_lm(i,j,kelm))

        ! pressure on the full level ke
        zfakt = 0.5*(hhl_lm(i,j,kelm)-hhl_lm(i,j,kelm+1)) * g / ( ztvw*r_d )
        zphf  = ps_lm(i,j) * EXP( -zfakt )
        pp_lm(i,j,kelm) = zphf - p0_lm(i,j,kelm)

        zt0  = p0_lm(i,j,kelm) / r_d / rho0(i,j)
        ! contribution of the virtual temperature to the buoyancy term
        ztvdt(i,j,klow) = (ztvw - zt0) / t_lm(i,j,kelm)
        ! coefficient of the pressure contribution to the buoyancy term
        zt0dp0t(i,j,klow) = 0.5 / ( r_d * rho0(i,j) * t_lm(i,j,kelm) )

      ENDDO
    ENDDO

    ! pressure deviation on the full level k-1
    DO k = kelm,2, -1
      DO j = 1,je2lm
        DO i = 1,ie2lm
          IF (iterate == 1) pp_lm(i,j,k-1) = pp_lm(i,j,k)

          ! virtual temperature at levels k and k-1
          ztvw1 = t_lm(i,j,k-1) * (1.0+rvd_m_o*qv_lm(i,j,k-1)-qc_lm(i,j,k-1))
          ztvw  = t_lm(i,j,k  ) * (1.0+rvd_m_o*qv_lm(i,j,k  )-qc_lm(i,j,k  ))

          pp_lm(i,j,k-1) = pp_lm(i,j,k) + 0.125_ireals*(hhl_lm(i,j,k-1)-hhl_lm(i,j,k+1))*g*&
                           (rho0_lm(i,j,k-1)+rho0_lm(i,j,k))*( (ztvw-t0_lm(i,j,k))/&
                           (ztvw) - t0_lm(i,j,k)/ztvw*pp_lm(i,j,k)/p0_lm(i,j,k) + &
                           (ztvw1-t0_lm(i,j,k-1))/ztvw1 - t0_lm(i,j,k-1)/ztvw1* &
                           pp_lm(i,j,k-1)/p0_lm(i,j,k-1) )

          ! dieser Ausdruck kann sicher schoener hingeschrieben werden.
          ! Ausserdem werden viele Terme doppelt berechnet        
        ENDDO
      ENDDO

      ! changing the organizational indices
      klow = 3 - klow
      kup  = 3 - kup

    ENDDO
   ENDDO ! iterate
  ENDIF

!------------------------------------------------------------------------------
! Section 8: Horizontal filtering of the pressure deviation
!------------------------------------------------------------------------------

  IF (lfilter_pp) THEN
    klow = kelm
    kup  = 1
    klow = MAX( 10, klow )
    nfilt = 2

    ! width of the stencil for the horizontal filter
    hfwidth  = 4
    hfw_m_nb = hfwidth - nboundlines
    ie2lm_hf = ie2lm + 2*hfw_m_nb
    je2lm_hf = je2lm + 2*hfw_m_nb
    ALLOCATE( pp_tmp(ie2lm_hf,je2lm_hf,kelm), STAT = istata )

    CALL extend_field ( pp_lm (:,:,:), ie2lm,    je2lm,                      &
                        pp_tmp(:,:,:), ie2lm_hf, je2lm_hf, kelm,             &
                        hfw_m_nb, nboundlines, hfjstartpar, hfjendpar,       &
                        sendbuf, isendbuflen, imp_reals, icomm_cart,         &
                        my_cart_neigh, num_compute, .FALSE., .FALSE., .FALSE.)

    DO n = 1, nfilt

      CALL horizontal_filtering( pp_tmp(:,:,:), ie2lm_hf, je2lm_hf, kelm,    &
                                 nboundlines, hfwidth, 4, my_cart_neigh)

      ! exchange boundaries after filtering
      IF (num_compute > 1) THEN
      kzdims(1:24) =(/ kelm,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
      CALL exchg_boundaries                                                      &
        ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,           &
          ie2lm_hf, je2lm_hf, kzdims, hfjstartpar, hfjendpar, hfwidth, hfwidth,  &
          my_cart_neigh, .FALSE., .FALSE., .FALSE.,                              &
          2000, .FALSE., 1, izerror, yzerrmsg,                                   &
          pp_tmp(:,:,:) )
      ENDIF

    END DO

    DO k = kup, klow
      DO j = 1, je2lm
        DO i = 1, ie2lm
          pp_lm(i,j,k) = pp_tmp(i+hfw_m_nb,j+hfw_m_nb,k)
        END DO
      END DO
    END DO

    DEALLOCATE( pp_tmp, STAT = istatd )
  ENDIF

!------------------------------------------------------------------------------
! Section 9: Final splitting of relative humidity (after filter and balance)
!------------------------------------------------------------------------------

  IF (lfilter_pp .OR. lbalance_pp) THEN
    rho0(:,:) = dp0_lm(:,:,kelm) / ( hhl_lm(:,:,kelm) - hhl_lm(:,:,kelm+1) ) / g

    DO k = 1, kelm
      DO j = 1, je2lm
        DO i = 1, ie2lm
          ! split grh_lm into qv_lm, qc_lm using actual p = pp_lm + p0_lm
          zaq = sf_psat_w (t_lm(i,j,k), b1, b2_w, b3, b4_w)
          zpm = p0_lm(i,j,k) + pp_lm(i,j,k)
          zbq = sf_qsat   (zaq, zpm, Rdv, O_m_rdv)
          qv_lm(i,j,k) = MIN(1.0_ireals, grh_lm(i,j,k)) * zbq
          qc_lm(i,j,k) = MAX(0.0_ireals, grh_lm(i,j,k)-1.0_ireals) * zbq
          IF (qv_lm(i,j,k) < qvmin) qv_lm(i,j,k)  = qvmin
          IF (qc_lm(i,j,k) < qcmin) qc_lm(i,j,k)  = 0.0_ireals
        ENDDO
      ENDDO
    ENDDO

    qrs(:,:) = 0.0  !  no water loading
    ps_lm(:,:) = (p0_lm(:,:,kelm)+pp_lm(:,:,kelm)) * EXP( 0.5*dp0_lm(:,:,kelm) / &
       ( t_lm(:,:,kelm)*(1.0+rvd_m_o*qv_lm(:,:,kelm)-qc_lm(:,:,kelm)-qrs(:,:)) * &
       r_d*rho0(:,:) ) )
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE org_vert_inter_lm

!===============================================================================
!+ Vertical Interpolation to LM levels.
!------------------------------------------------------------------------------

SUBROUTINE vert_interp (xlm, yname, idim_lm, height, idim_in, kgr, idebug)

!------------------------------------------------------------------------------
!
! Description:
!   Does the vertical interpolation on the sigma reference levels of the
!   variables: u_lm, v_lm, w_lm, t_lm, pp_lm, grh_lm and qi_lm
!   (grh = generalized relative humidity)
!
! Method:
!   The routine is entered with the variable xlm on the height levels of the 
!   input profile (coarse grid). An intermediate profile is constructed by
!   shifting the input profile in order to take into account the topographical
!   difference between the coarse and the fine grid models (variable zxexp_vec
!   on the levels zpexp_vec).
!   For each gridpoint the tension spline routine tautsp is called for
!   the interpolation.
!   For break(i) <= x <= break(i+1) the interpolation function has the form
!      F(X) = COEF(1,I)+DX(COEF(2,I)+DX/2(COEF(3,I)+DX/3(COEF(4,I)))
!      using DX=X-BREAK(I) and i=1,...,izl
!   This interpolation leads to the output profile with the variable xlm on the 
!   fine grid levels zhhl.
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers),   INTENT(IN)   ::  &
  idim_lm,      & ! vertical dimension of FINE LM fields
  idim_in,      & ! vertical dimension of COARSE LM fields
  idebug,       & ! for debug output
  kgr             ! height level of boundary layer top

REAL (KIND=ireals),         INTENT(IN)   ::  &
  height (ie2lm,je2lm,idim_in)  ! height on COARSE LM levels used as abscissas

REAL (KIND=ireals),         INTENT(INOUT)::  &
  xlm (ie2lm,je2lm,idim_lm)     ! field to be interpolated

CHARACTER (LEN=*),          INTENT(IN)   ::  &
  yname           ! name of the variable

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers) ::     &
  i, j, k, kk, k1, k2, k3,      &
  n, kstart, kend, izerror, idone,    &
  izln_vec(ie2lm),              & ! number of abscissas for spline
  nztau_vec(ie2lm),             & ! number of abscissas
  kzgrn_vec(ie2lm),             & !
  izn      (ie2lm),             & !
  izind_vec(ie2lm,ke1lm)

LOGICAL                        ::     &
  ldone(ie2lm)

REAL    (KIND=ireals)    ::   &
  zdx, zrfmax(ie2lm), zp1(ie2lm), zp2(ie2lm), zp3(ie2lm), &
  zcheck(ie2lm),              & !
  zdelh(ie2lm), zdelhchk,     & ! height difference
  zgamma                        ! tension-parameter = 5.5

! Local arrays:
REAL (KIND=ireals)       ::   &
  zhhl       (ie2lm,ke1lm),   & !
  zpexp_vec  (ie2lm,idim_in+6), & !
  zxexp_vec  (ie2lm,idim_in+6), & !
  zbreak_vec (ie2lm,3*kedim)    ! abscissas of spline

REAL (KIND=ireals)       ::   &
  zcoef_vec (ie2lm,4,3*kedim),& ! coefficients  of spline
  zs_vec    (ie2lm,idim_in+6,6),& ! work array for tautsp
  zbx    (ie2lm,je2lm),       & !
  zdelx  (ie2lm,je2lm),       & !
  zpq    (ie2lm,je2lm),       & ! work arrays for regression
  zpq2   (ie2lm,je2lm),       & !
  zxq    (ie2lm,je2lm),       & !
  zxpq   (ie2lm,je2lm)

CHARACTER (LEN=80)       ::   &
  yzerrmsg        ! for error message

CHARACTER (LEN=20)       ::   &
  yzroutine       ! name of the routine

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 0: Initializations
!------------------------------------------------------------------------------

  izerror   = 0
  yzerrmsg  = '   '
  yzroutine = 'vert_interp'
  zgamma    = 5.5
  zdelhchk  = 1.0E-4_ireals

! Initialize intervals
  izind_vec(:,:) = 0.0

!roa begin                 braucht man evtl nicht unbedingt
! Initializations
! nztau_vec(:)=0.0
! zbx(:,:)=0.0
! zdelx(:,:)=0.0
! zdelh=0.0
! zpexp_vec(:,:)=-23456789
! zxexp_vec(:,:)=-23456789
!roa end

  kstart = 1
  kend = kelm
  IF (yname == 'w') kend = ke1lm

  IF (idebug > 15) THEN
    PRINT *, 'vertical interpolation of:  ', yname, kend, ke1lm, kedim, idim_in
  ENDIF

!------------------------------------------------------------------------------
! Section 1: Loop over all grid points
!------------------------------------------------------------------------------

  ! Construct intermediate profiles that can be used for interpolation:
  !   variable zxexp_vec on the levels zpexp_vec.
  ! The levels in zpexp_vec are numbered from the ground to the top of the
  ! atmosphere (just the inverse counting compared to the atmospheric 
  ! COSMO-Model variables). The first levels to be set are:
  !   zpexp_vec(:,1): under the ground
  !   zpexp_vec(:,2): at the soil surface
  !   zpexp_vec(:,3): 1st atmospheric level
  ! All other levels are taken from the coarse grid input model
  !   zpexp_vec(:,4:idim_in+2):  height(:,j;idim_in+1-2:1)

  DO j = 1, je2lm

    ! 1.1 Lower boundary values
    ! -------------------------
    ! The first value that is set is the surface value; for w, an artificial
    ! level of 20 m below the fine grid surface is introduced
    ! For u and v-gridpoints, the height is averaged to the C-grid.
    ! On the eastern and northern boundary the values from the last mass grid
    ! point are used. For subdomains with a right or upper neighbor, this does
    ! not matter, for the rightmost subdomain, it is the only value available

    IF     (yname == 'w') THEN
       zpexp_vec(:,2) = hhl_lm(:,j,ke1lm)-20.0_ireals
    ELSEIF (yname == 'u') THEN
      DO i = 1, ie2lm-1
        zpexp_vec(i,2) = 0.5_ireals * (hhl_lm(i,j,ke1lm) + hhl_lm(i+1,j,ke1lm))
      ENDDO
      zpexp_vec(ie2lm,2) = hhl_lm(ie2lm,j,ke1lm)
    ELSEIF (yname == 'v') THEN
      IF (j < je2lm) THEN
        zpexp_vec(:,2) = 0.5_ireals * (hhl_lm(:,j,ke1lm) + hhl_lm(:,j+1,ke1lm))
      ELSE
        zpexp_vec(:,2) = hhl_lm(:,je2lm,ke1lm)
      ENDIF
    ELSE   ! t, grh, qi, pp and chemistry variables
      zpexp_vec(:,2) = hhl_lm(:,j,ke1lm)
    ENDIF

    kzgrn_vec(:) = 0_iintegers

    ! For all variables, an artificial value of 50 m below the surface is
    ! introduced for the intermediate profile
    IF     (yname == 'w') THEN
      zpexp_vec(:,1) = zpexp_vec(:,2)- 30.0_ireals   ! overall: 30 m below surface
    ELSE
      zpexp_vec(:,1) = zpexp_vec(:,2)- 50.0_ireals   ! also: 50 m below surface
    ENDIF

    IF     (yname == 'u') THEN
      DO i = 1, ie2lm-1
        zpexp_vec(i,3) = 0.5_ireals * (height(i,    j,idim_in) &
                                     + height(i+1,  j,idim_in))
      ENDDO
      zpexp_vec(ie2lm,3) =             height(ie2lm,j,idim_in)
    ELSEIF (yname == 'v') THEN
      IF (j < je2lm) THEN
        zpexp_vec(:,3) = 0.5_ireals * (height(:,    j,idim_in) &
                                     + height(:,  j+1,idim_in))
      ELSE
        zpexp_vec(:,3) =               height(:,je2lm,idim_in)
      ENDIF
    ELSE   ! w, t, grh, qi, pp and chemistry variables
      zpexp_vec(:,3) =                 height(:,    j,idim_in)
    ENDIF

    IF (yname == 't') THEN
      zxexp_vec(:,1) = t_s_gl(:,j)
      zxexp_vec(:,2) = t_s_gl(:,j)
      zxexp_vec(:,3) = xlm(:,j,idim_in)
    ELSE IF (yname == 'u') THEN
      zxexp_vec(:,1) = 0.0
      zxexp_vec(:,2) = 0.0
      zxexp_vec(:,3) = xlm(:,j,idim_in)
    ELSE IF (yname == 'v') THEN
      zxexp_vec(:,1) = 0.0
      zxexp_vec(:,2) = 0.0
      zxexp_vec(:,3) = xlm(:,j,idim_in)
    ELSE IF (yname == 'grh' .OR. yname == 'qi' .OR.                    &
             yname == 'qr'  .OR. yname == 'qs' .OR. yname == 'qg' .OR. &
             yname == 'w'   .OR. yname == 'pp' .OR. yname == 'p') THEN
      zxexp_vec(:,1) = xlm(:,j,idim_in)
      zxexp_vec(:,2) = xlm(:,j,idim_in)
      zxexp_vec(:,3) = xlm(:,j,idim_in)
    ELSE IF (yname == 'logp') THEN
      zxexp_vec(:,3) = xlm(:,j,idim_in)
      zxexp_vec(:,2) = xlm(:,j,idim_in)+(zpexp_vec(:,3)- zpexp_vec(:,2))*0.000122_ireals
      zxexp_vec(:,1) = xlm(:,j,idim_in)+(zpexp_vec(:,3) + 50._ireals - zpexp_vec(:,2))*0.000122_ireals
    ELSE ! chemistry variables
      zxexp_vec(:,1) = xlm(:,j,idim_in)
      zxexp_vec(:,2) = xlm(:,j,idim_in)
      zxexp_vec(:,3) = xlm(:,j,idim_in)
    END IF

    ! 1.2 Levels in the atmosphere
    ! ----------------------------
    DO k = 2, idim_in
      ! For u and v-gridpoints, the height is averaged to the C-grid.
      ! On the eastern and northern boundary the values from the last mass grid
      ! point are used. For subdomains with a right or upper neighbor, this does
      ! not matter, for the rightmost subdomain, it is the only value available
      IF     (yname == 'u') THEN
        DO i = 1, ie2lm-1
          zpexp_vec(i,k+2) = 0.5_ireals *(height(i,    j,idim_in+1-k) &
                                        + height(i+1,  j,idim_in+1-k))
        ENDDO
        zpexp_vec(ie2lm,k+2) =            height(ie2lm,j,idim_in+1-k)
      ELSEIF (yname == 'v') THEN
        IF (j < je2lm) THEN
          zpexp_vec(:,k+2) = 0.5_ireals *(height(:,    j,idim_in+1-k) &
                                        + height(:,  j+1,idim_in+1-k))
        ELSE
          zpexp_vec(:,k+2) =              height(:,je2lm,idim_in+1-k)
        ENDIF
      ELSE
        zpexp_vec  (:,k+2) =              height(:,    j,idim_in+1-k)
      ENDIF
      zxexp_vec    (:,k+2) =                 xlm(:,    j,idim_in+1-k)
    ENDDO

    ! 1.3 Set the height vector for the fine model at the correct grid points
    !     (output profiles)
    ! -----------------------------------------------------------------------
    IF     (yname == 'w') THEN
      DO k = 1, kend
        zhhl(:,k) = hhl_lm(:,j,k)
      ENDDO
    ELSEIF (yname == 'u') THEN
      DO k = 1, kend
        DO i = 1, ie2lm-1
          zhhl(i,k) = 0.25 * (hhl_lm(i,j,k  ) + hhl_lm(i+1,j,k  )      &
                            + hhl_lm(i,j,k+1) + hhl_lm(i+1,j,k+1))
        ENDDO
        zhhl(ie2lm,k)= 0.5 * (hhl_lm(ie2lm,j,k) + hhl_lm(ie2lm,j,k+1))
      ENDDO
    ELSEIF (yname == 'v') THEN
      DO k = 1, kend
        IF (j < je2lm) THEN
          zhhl(:,k) = 0.25 * (hhl_lm(:,j,k  ) + hhl_lm(:,j+1,k  )      &
                            + hhl_lm(:,j,k+1) + hhl_lm(:,j+1,k+1))
        ELSE
          zhhl(:,k) =  0.5 * (hhl_lm(:,je2lm,k) + hhl_lm(:,je2lm,k+1))
        ENDIF
      ENDDO
    ELSE
      DO k = 1, kend
        zhhl(:,k) = 0.5 * (hhl_lm(:,j,k) + hhl_lm(:,j,k+1))
      ENDDO
    ENDIF

    nztau_vec(:) = idim_in + 2

    !--------------------------------------------------------------------------
    ! For the temperature a linear approximation of the vertical profile above 
    ! the boundary layer is computed.
    ! This linear profile is used for the extrapolation of the fields.
    !
    ! The index 'kgr' of the boundary layer top.
    !
    ! The linear approximation is of the form:
    !  x(h) = bx*h + cx  ,  where h is the height of input levels.
    !
    ! We need only the coefficient bx and delx . We do not need cx.
    !
    ! Version 1.9: 
    !  The computation of this regression had a severe bug before, where
    !  the index kgr was used for the intermediate profile, but kgr refers
    !  to the indexing in the initial profile. This has been corrected by
    !  Anne Roches, MCH.
    !
    !  Also, the regression is now computed only for temperature, not for the
    !  horizontal wind speeds any more.
    !
    !--------------------------------------------------------------------------
     IF (yname == 't') THEN
      zpq  (:,j) = 0.0
      zpq2 (:,j) = 0.0
      zxq  (:,j) = 0.0
      zxpq (:,j) = 0.0
      DO k = 3+idim_in-kgr,3+idim_in-kgr+3
        zpq (:,j) = zpq (:,j) + zpexp_vec(:,k)
        zpq2(:,j) = zpq2(:,j) + zpexp_vec(:,k) * zpexp_vec(:,k)
        zxq (:,j) = zxq (:,j) + zxexp_vec(:,k)
        zxpq(:,j) = zxpq(:,j) + zxexp_vec(:,k) * zpexp_vec(:,k)
      ENDDO
      zbx  (:,j) =( zxpq(:,j) - zxq(:,j)*zpq(:,j)/4.0_ireals) /    &
                   (zpq2(:,j) - zpq(:,j)*zpq(:,j)/4.0_ireals)
      zdelx(:,j) =  zbx (:,j) * (zhhl(:,kend) - zpexp_vec(:,3))
    ENDIF
  
    ! 1.4 Shift boundary layer depending on height differences
    ! --------------------------------------------------------

    DO i = 1, ie2lm
      ! Height difference LM minus COARSE LM (horizontaly interpolated)
      zdelh(i)  = zhhl(i,kend) - zpexp_vec(i,3)
      zcheck(i) = zpexp_vec(i,3)
    ENDDO
 
    ! Shifting within the boundary layer:
    ! Extrapolation of the profiles for hhl_gl > zhhl :
    !  for yname= 't',         with a linear regression
    !  for all other variables with a constant value.
    IF (yname == 't') THEN
      DO k = idim_in, kgr+1, -1
        kk = 3+idim_in-k
        DO i = 1, ie2lm
          IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
            zpexp_vec(i,kk) = zpexp_vec(i,kk) + zdelh(i)
            zxexp_vec(i,kk) = xlm     (i,j,k) + zdelx(i,j) ! shift all PBL
          ENDIF
        ENDDO
      ENDDO
    ELSEIF (yname /= 'logp') THEN   ! 'rh', 'qi', 'pp', 'p' and 'w'
      DO k = idim_in, kgr+1, -1
        kk = 3+idim_in-k
        DO i = 1, ie2lm
          IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
            zpexp_vec(i,kk) = zpexp_vec(i,kk) + zdelh(i)
            zxexp_vec(i,kk) = xlm     (i,j,k)        ! shift PBL height ONLY
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! Shift ground values
    IF (yname == 't') THEN
      IF (llm2lm) THEN
        DO i = 1, ie2lm
          zxexp_vec(i,2)= zxexp_vec(i,2)+zbx(i,j)*(zpexp_vec(i,2)-hhl_gl(i,j,ke_in+1))
          zxexp_vec(i,1)= zxexp_vec(i,2)
        ENDDO
       ELSEIF (lum2lm .OR. (lcm2lm .AND. lcm_hgt_coor)) THEN
        DO i = 1, ie2lm
          zxexp_vec(i,2)= zxexp_vec(i,2)+zbx(i,j)*(zpexp_vec(i,2)-hsurf_gl(i,j))
          zxexp_vec(i,1)= zxexp_vec(i,2)
        ENDDO
       ENDIF
    ENDIF

    ! At the top of the boundary layer
    IF (yname /= 'logp') THEN
      k1 = 3+idim_in-kgr+2
      k2 = 3+idim_in-kgr+1
      k3 = 3+idim_in-kgr
      DO i = 1, ie2lm
        IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
          zp1(i) = zpexp_vec(i,k1)
          zp2(i) = zpexp_vec(i,k2)
          zp3(i) = zpexp_vec(i,k3)
          zpexp_vec(i,k1) = zp1(i) + zdelh(i)/4.0_ireals * 1.0_ireals
          zpexp_vec(i,k2) = zp2(i) + zdelh(i)/4.0_ireals * 2.0_ireals
          zpexp_vec(i,k3) = zp3(i) + zdelh(i)/4.0_ireals * 3.0_ireals
        ENDIF
      ENDDO

      IF (yname == 't') THEN
        DO i = 1, ie2lm
          IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
            ! add 2 regression values
            zxexp_vec(i,k1) = xlm(i,j,kgr-2) + zbx(i,j)*zdelh(i)/4.0_ireals * 1.0_ireals
            zxexp_vec(i,k2) = xlm(i,j,kgr-1) + zbx(i,j)*zdelh(i)/4.0_ireals * 2.0_ireals
            zxexp_vec(i,k3) = xlm(i,j,kgr  ) + zbx(i,j)*zdelh(i)/4.0_ireals * 3.0_ireals
          ENDIF
        ENDDO
      ELSE   ! 'rh', 'qi', 'pp' and 'w'
        DO i = 1, ie2lm
          IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
            ! add transition zone values
            zxexp_vec(i,k1) = xlm(i,j,kgr-2)
            zxexp_vec(i,k2) = xlm(i,j,kgr-1)
            zxexp_vec(i,k3) = xlm(i,j,kgr  )
          ENDIF
        ENDDO
      ENDIF

      ! upper part of the profile

      DO k = kgr-3, 1, -1              ! TOP
        kk = 3+idim_in-k
        DO i = 1, ie2lm
          IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
            zp1(i) = zpexp_vec(i,kk)   ! no change
            zpexp_vec(i,kk) = zp1(i)
            zxexp_vec(i,kk) = xlm(i,j,k)
          ENDIF
        ENDDO
      ENDDO

    ELSE   ! yname = 'logp'
      ! at the moment, this is for HadGEM2 data
      DO i = 1, ie2lm
        zp1(i) = zpexp_vec(i,idim_in+2)
      ENDDO
    ENDIF

    DO i = 1, ie2lm
      IF (zcheck(i) > zhhl(i,kend) + zdelhchk) THEN
        zpexp_vec(i,idim_in + 4) = zp1(i) + 50.0_ireals
        zpexp_vec(i,idim_in + 3) = zp1(i) + 20.0_ireals
        zxexp_vec(i,idim_in + 4) = xlm(i,j,1)
        zxexp_vec(i,idim_in + 3) = xlm(i,j,1)
        nztau_vec(i) = idim_in + 4
      ENDIF
    ENDDO

    ! Elevate (shift) the boundary layer for hhl_gl < zhhl :
    DO i = 1, ie2lm
      IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
        k = 3+idim_in-kgr
        DO WHILE (zpexp_vec(i,k) < zpexp_vec(i,3+idim_in-kgr) + zdelh(i))
          k = k + 1
        ENDDO
        kzgrn_vec(i) = k-3-idim_in+kgr
      ENDIF
    ENDDO

    IF (yname == 't') THEN
      DO k = idim_in, kgr, -1
        kk = 3+idim_in-k
        DO i = 1, ie2lm
          IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
            zpexp_vec(i,kk) = zpexp_vec(i,kk) + zdelh(i)
            zxexp_vec(i,kk) = xlm     (i,j,k) + zdelx(i,j) ! shift all PBL
          ENDIF
        ENDDO
      ENDDO
    ELSEIF (yname == 'logp') THEN
      DO i = 1, ie2lm
        IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
          zpexp_vec(i,3) = zpexp_vec(i,3+kzgrn_vec(i))
          zxexp_vec(i,3) = xlm (i,j,idim_in-kzgrn_vec(i))
          zxexp_vec(i,2) = zxexp_vec(i,3)+(zpexp_vec(i,3) -                &
                           zpexp_vec(i,2))*0.000122_ireals
          zxexp_vec(i,1) = zxexp_vec(i,3)+(zpexp_vec(i,3) +                &
                           50._ireals - zpexp_vec(i,2))*0.000122_ireals
        ENDIF
      ENDDO
    ELSE   ! 'rh', 'qi', 'pp' and 'w'
      DO k = idim_in, kgr, -1
        kk = 3+idim_in-k
        DO i = 1, ie2lm
          IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
            zpexp_vec(i,kk) = zpexp_vec(i,kk) + zdelh(i)
            zxexp_vec(i,kk) = xlm     (i,j,k)        ! shift PBL height ONLY
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    DO k = kgr-1, 1, -1     ! LEAVE OUT zdelh (kzgrn_vec(i)) LEVELS
      DO i = 1, ie2lm
        IF (k >= kzgrn_vec(i)+1) THEN
          kk = 3+idim_in-k
          IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
            zxexp_vec(i,kk) = xlm     (i,j,k - kzgrn_vec(i))
            zpexp_vec(i,kk) = zpexp_vec(i,kk + kzgrn_vec(i))
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    DO i = 1, ie2lm
      IF (zcheck(i) < zhhl(i,kend) - zdelhchk) THEN
        nztau_vec(i)= idim_in + 2 - kzgrn_vec(i)
      ENDIF
    ENDDO

    ! 1.5. Upper boundary values
    ! --------------------------
    DO i = 1, ie2lm
      k = nztau_vec(i) + 2
      nztau_vec(i) = k
      zpexp_vec(i,k  ) = zpexp_vec(i,k-2) + 2000.0_ireals
      zpexp_vec(i,k-1) = zpexp_vec(i,k-2) + 1000.0_ireals
      zxexp_vec(i,k  ) = xlm(i,j,1)
      zxexp_vec(i,k-1) = xlm(i,j,1)
    ENDDO

    !------------------------------------------------------------------------
    ! Section 2.0: Compute tension splines
    !------------------------------------------------------------------------

    izln_vec(:) = 3*kedim

    CALL tautsp2D(zpexp_vec, zxexp_vec, nztau_vec, ie2lm, 1, ie2lm, idim_in+6, &
                  zgamma, zs_vec, zbreak_vec, zcoef_vec, izln_vec, izerror)

    IF (izerror == 0) THEN

      ! 2.1 Locate intervals and put in izind
!US   izn(:) = izln_vec(:) also bug fix as in src_lm_fields
      izn(:) = izln_vec(:) - 1
      DO k = kstart, kend
        ldone(:) = .FALSE.
        idone    = 0
        DO WHILE (idone < ie2lm)
          DO i = 1, ie2lm
            IF (.NOT. ldone(i)) THEN
              IF ( (zbreak_vec(i,izn(i)) <= zhhl(i,k))             .AND. &
                                           (zhhl(i,k) <= zbreak_vec(i,izn(i)+1)) ) THEN
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

      ! 2.2 Compute interpolated values
      !
      DO k = kstart, kend
        DO i = 1, ie2lm
          zdx        = zhhl(i,k) - zbreak_vec(i,izind_vec(i,k))
          xlm(i,j,k) = zcoef_vec(i,1,Izind_vec(i,k)) +          &
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

        DO k = 1, kend
          DO i = 1, ie2lm
            xlm(i,j,k) = MAX ( 0.001_ireals, xlm(i,j,k) )
            xlm(i,j,k) = MIN ( xlm(i,j,k), zrfmax(i) )
          ENDDO
        ENDDO
      ENDIF

      ! Limit values of qi
      IF (yname=='qi' .OR. yname=='qr' .OR. yname=='qs' .OR. yname=='qg') THEN
        DO k = 1, kend
          DO i = 1, ie2lm
            IF (xlm(i,j,k) < qimin) THEN
              xlm(i,j,k) = 0.0_ireals
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      ! Set w_lm to 0.0 at the top and the bottom
      IF (yname == 'w') THEN
        DO i = 1, ie2lm
          xlm(i,j,    1) = 0.0_ireals
          xlm(i,j,ke1lm) = 0.0_ireals
        ENDDO
      ENDIF

      ! Limit the values of the chemistry fields to 0.0
      IF (yname/='u' .AND. yname/='v'  .AND. yname/='w'   .AND. yname/='pp' &
                     .AND. yname/='t'  .AND. yname/='grh' .AND. yname/='qi' &
                     .AND. yname/='qr' .AND. yname/='qs'  .AND. yname/='qg' &
                     .AND. yname/='p') THEN
        DO k = 1, kend
          DO i = 1, ie2lm
            IF (xlm(i,j,k) < 1E-12_ireals) THEN
              xlm(i,j,k) = 0.0_ireals
            ENDIF
          ENDDO
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

END SUBROUTINE vert_interp

!===============================================================================

END MODULE src_vert_inter_lm
