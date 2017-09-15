!+ Source Module for interpolating IFS/coarser-LM-fields to LM-fields
!==============================================================================

MODULE src_coarse_interpol

!==============================================================================
!
! Description:
!   This module contains subroutines necessary for interpolating coarser
!   fields to the fine grid LM-fields.
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
!  Corrected dimensions of field zdt_so in SR interpol_coarse_special_lm
!  Interpolate T_SO(0) also for computing boundary values
! 1.3        2005/12/12 Ulrich Schaettler
!  Added fields for treatment of rho_snow
! V1_5         2007/07/09 Ulrich Schaettler
!  Replaced ke_soil to ke_soil_lm; czmls to czmls_lm
!  Eliminated un-used variables
!  Improved snow interpolation for ec2lm: if there is no snow, set, t_snow
!  to t_g1 before the horizontal interpolation (Davide Cesari)
!  Adapted debug output to idbg_level
! V1_6         2007/09/07 Ulrich Schaettler, Burkhardt Rockel, Uwe Boehm
!  Changed interfaces to interpolation routines for Cressman scheme
!  No variable loop cycling for optional variables
!  Added Subroutine interpol_coarse_special_cm
! V1_7         2007/11/26 Ulrich Schaettler
!  Corrected dimension for w_so_lm
! V1_8         2008/05/29 Hans-Juergen Panitz, Uwe Boehm
!  Interpolation of SST in case of lcm2lm not required
!  Introduced alternative calls to bicubic interpolation
! V1_9         2009/09/03 Hans-Juergen Panitz, Guy DeMorsier
!  Reset qc_lm in case lcm2lm, if it has not been read
!  Eliminated Section 8 (add qi to qc) in SR interpol_coarse_special_lm
!  Implemented option to choose new soil properties with l_smi (Guy)
!  Added additional debug output
! V1_10        2009/12/17 Ulrich Schaettler
!  Adaptations to process Unified Model data: UM data are on an unrotated grid
!  without specifying pollat, pollon; take care of not rotating them
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler
!  Adaptations to the modified vartab input table
!  Additional debug output
! V1_19        2012/06/06 Ulrich Schaettler
!  Added lhir2lm as internal logical flag
!  Added new SR interpol_coarse_special_hir
! V1_20        2012/09/03 Ulrich Schaettler
!  Removed unused variables from USE-lists and added some more debug prints
! V1_22        2013/07/11 Ulrich Schaettler, Davide Cesari
!  Removed unused variables from the USE sections
!  New branch in the interpol_coarse_special_ec routine for interpolating
!    in multi-layer style (Davide)
! V1_23        2013/10/02 Lucio Torrisi
!  IFS soiltyp has 7 soil types, not 6
!
! Added height correction for soil isotopes. Hui Tang 2013-11-20
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

USE data_fields_in, ONLY: &
 hsurf_in  ,            & ! orography                                   (  m  )
 soiltyp_in,            & ! type of the soil (keys 0-9)                   --
 t_s_in    ,            & ! temperature of the ground surface           (  K  )
 t_g1_in   ,            & ! temperature of first soil layer             (  K  )
 t_g2_in   ,            & ! temperature of second soil layer            (  K  )
 t_g3_in   ,            & ! temperature of third soil layer             (  K  )
 qv_s_in   ,            & ! specific water vapor content on the surface (kg/kg)
 w_g1_in   ,            & ! water content of the upper soil layer       (m H2O)
 w_g2_in   ,            & ! water content of the medium soil layer      (m H2O)
 w_g3_in   ,            & ! water content of the deepest soil layer     (m H2O)
 t_so_in   ,            & ! temperature for new multi layer soil model  (  K  )
 w_so_in   ,            & ! soil moisture for multi layer soil model    (m H2O)
 lolp_in   ,            & ! Land Sea Mask of input for 'M'atch Interpolation
 u_in      ,            & ! zonal wind speed                            ( m/s )
 v_in      ,            & ! meridional wind speed                       ( m/s )
 t_in      ,            & ! temperature                                 (  K  )
 qv_in     ,            & ! specific water vapor content                (kg/kg)
 qc_in     ,            & ! specific cloud water content                (kg/kg)
 qi_in                    ! specific cloud ice water content            (kg/kg)

USE data_fields_in, ONLY: &
 t_g_in    ,            & ! temperature                                 (  K  )
 t_m_in    ,            & ! temperature                                 (  K  )
 t_cl_in   ,            & ! temp.  between medium and lower soil layer  (  K  )
 t_skin_in ,            & ! skin temperature of the ground surface      (  K  )
 t_snow_in ,            & ! temperature of the snow-surface             (  K  )
 w_snow_in ,            & ! water content of the snow                   (m H2O)
 w_i_in    ,            & ! water content of the interception storage   (m H2O)
 w_cl_in   ,            & ! climatological deep soil water content      (m H2O)
 dpsdt_in  ,            & ! surface pressure tendency                   (Pa/s )
 lat_coarse_m,          & ! latitudes of the LM grid points
 lon_coarse_m,          & ! longitudes of the LM grid points
! iso code Hui Tang 2013-11-20
 riso_in                  ! isotope ratios in water vapor;
! end iso code

!------------------------------------------------------------------------------

USE data_grid_in,    ONLY: &
  pollat_in,      & ! latitude of the rotated north pole (in degrees, N>0)
  pollon_in,      & ! longitude of the rotated north pole (in degrees, E>0)
  dlat_in,        & ! grid point distance in zonal direction (in degrees)
  dlon_in,        & ! grid point distance in meridional direction
  startlat_in_tot,& ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
  startlon_in_tot,& ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
  endlat_in_tot,  & ! transformed latitude of the upper right grid point
                    ! of the total domain (in degrees, N>0)
  endlon_in_tot,  & ! transformed longitude of the upper right grid point
                    ! of the total domain (in degrees, E>0)
  startlat_in,    & ! transformed latitude of the lower left grid point
                    ! of the local domain (in degrees, N>0)
  startlon_in,    & ! transformed longitude of the lower left grid point
                    ! of the local domain (in degrees, E>0)
  endlat_in,      & ! transformed latitude of the upper right grid point
                    ! of the local domain (in degrees, N>0)
  endlon_in,      & ! transformed longitude of the upper right grid point
                    ! of the local domain (in degrees, E>0)
  ie_in,          & ! ie for input grid, local domain
  je_in,          & ! je for input grid, local domain
  ie_in_tot,      & ! ie for input grid, total domain
  je_in_tot,      & ! je for input grid, total domain
  ke_in,          & ! ke for input grid
  ke1in,          & !
  ke_soil_in,     & ! number of levels in multi-layer soil model in input
  czmls_in,       & ! depth of the coarse grid soil layers in meters
  czhls_in,       & ! depth of the coarse grid half soil layers in meters
  latitudes_in,   & ! latitudes of the input data
  longitudes_in,  & ! longitudes of the input data
                    ! of the total domain (in degrees, E>0)
  grdpt_rel_in      ! relation between input and lm grid for cressman scheme

!------------------------------------------------------------------------------

USE data_fields_lm, ONLY : &
  t_snow_lm ,      & ! temperature of the snow surface                  (  K  )
  w_i_lm    ,      & ! water content of the interception storage        (m H2O)
  w_g1_lm   ,      & ! water content of the upper soil layer            (m H2O)
  w_g2_lm   ,      & ! water content of the medium soil layer           (m H2O)
  w_g3_lm   ,      & ! water content of the lower soil layer            (m H2O)
                     ! (if nl_soil_lm = 3, unused otherwise)
  w_snow_lm ,      & ! water content of the snow                        (m H2O)
  lolp_lm   ,      & ! Land Sea Mask of LM for 'M'atch Interpolation
  lmask_lm  ,      & ! mask of points on the frame
  ps_gl     ,      & ! surface pressure on the interpol. GME orogr.     ( Pa  )
  fis_gl    ,      & ! orography * G interpolated from input grid       (m2/s2)
  p0_gl     ,      & ! ref. pres. on full levels + interpol. COARSE LM oro.(Pa)
  dp0_gl    ,      & ! reference pressure thickness of layers           ( Pa  )
  rho0_gl   ,      & ! reference density at the full model levels       (kg/m3)
  hhl_gl    ,      & ! height of half-levels on the interpol. COARSE LM oro.(m)
  fic_gl    ,      & ! check level of geopotential                      (m2/s2)
  t_2m_gl   ,      & ! 2m temperature                                   (  K  )
  qv_2m_gl  ,      & ! 2m humidity                                      (kg/kg)
  t_s_gl    ,      & ! temperature of the ground surface                (  K  )
  t_skin_gl ,      & ! skin temperature of the ground surface           (  K  )  !_br
! SP, 201405
  t_cl_lm   ,      & ! temperature between medium and lower soil layer  (  K  )
  rh_s_gl   ,      & ! relative humidity at the surface                 (kg/kg)
  dtms_gl   ,      & ! t_m_lm    - t_s_lm                               (  K  )
  dtkes_gl  ,      & ! t(ke)_lm  - t_s_lm                               (  K  )
  dtssnow_gl,      & ! t_s_lm    - t_snow_lm                            (  K  )
  t_so_lm   ,      & ! multi-layer soil temperature                     (  K  )
  dt_so_gl  ,      & ! multi-layer soil temperature (for interpolation)
  w_so_lm   ,      & ! multi-layer soil moisture                        (m H2O)
! iso code Hui Tang 2013-11-20 
  risosoil_lm,     & ! isotope ratio in soil moisture;
  drisoke_gl         ! riso_lm(:,:,ke,1-2) - risosoil_lm (:,:,1-2) (1: 18O; 2: 2H)
! iso code end

USE data_fields_lm, ONLY : &
  i_index   ,      & !
  j_index   ,      & !
  x_wght    ,      & !
  y_wght    ,      & !
  lonlm_u   ,      & ! longitudes of the LM u grid points
  latlm_u   ,      & ! latitudes of the LM u grid points
  lonlm_v   ,      & ! longitudes of the LM v grid points
  latlm_v   ,      & ! latitudes of the LM v grid points
  u_lm      ,      & ! zonal wind speed                                 ( m/s )
  v_lm      ,      & ! meridional wind speed                            ( m/s )
  t_lm      ,      & ! temperature                                      (  K  )
  p_lm      ,      & ! full pressure (needed for lum2lm)                ( Pa  )
  pp_lm     ,      & ! deviation from the reference pressure            ( Pa  )
  qv_lm     ,      & ! specific water vapor content                     (kg/kg)
  qc_lm              ! specific cloud water content                     (kg/kg)

!------------------------------------------------------------------------------

USE data_grid_lm,    ONLY: &
  pollat,      & ! latitude of the rotated north pole (in degrees, N>0)
  pollon,      & ! longitude of the rotated north pole (in degrees, E>0)
  ie2lm,       & !
  je2lm,       & !
  ke_soil_lm,  & ! number of levels in multi-layer soil model in output
  czhls_lm       ! depth of the half soil layers in meters in output

!------------------------------------------------------------------------------

USE data_int2lm_control,    ONLY: &
  luse_t_skin,  & ! if .TRUE., use ECMWF skin temperature for surface
  l_smi,        & ! if .TRUE., interpolate soil moisture with SMI
  lcomp_bound,  & ! compute fields for boundaries
  lmulti_layer_lm, & ! if .TRUE., compute soil fields for multi-layer soil model
  lmulti_layer_in, & ! if .TRUE., incoming data from new multi-layer soil model
  lbd_frame_cur,& ! if .TRUE., current boundary fields include only frames
  lgsm2lm,      & ! if .TRUE., gsm ->lm
  lgfs2lm,      & ! if .TRUE., gfs ->lm
  lec2lm,       & ! if .TRUE., ec ->lm
  llm2lm,       & ! if .TRUE., lm ->lm
  lum2lm,       & ! if .TRUE., um ->lm
  lhir2lm,      & ! if .TRUE., hirlam ->lm
  lcm2lm,       & ! if .TRUE., cm ->lm  !_br
  l_cressman,   & ! logical switch for controling the use of a Cressman scheme
  lprog_qi,     & ! if .TRUE., interpolate qi from GME to LM grid
  idbg_level,   & ! to control verbosity of output
  lprintdeb_all,& ! whether all PEs print debug output
  nl_soil_lm,   & ! number of soil layers in LM
  nl_soil_in,   & ! number of soil layers in GME
  l_bicub_spl     ! switch for using a bicubic spline interpolation

!------------------------------------------------------------------------------

USE data_int2lm_io,        ONLY : &
  nvar_lm,           & ! maximum number of variables in fine grid LM variable table
  nvar_in,           & ! actual maximum number of variables in input variable table
  undef,             & ! the same with other KIND-Parameter
  var_lm,            & ! array for fine grid LM variable table
  var_in               ! array for input model variable table

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY :  &
    my_cart_id         ! rank of this subdomain in the cartesian communicator

!------------------------------------------------------------------------------

USE data_int2lm_constants,     ONLY :  &
    r_d,        & ! gas constant for dry air                      [J/K*kg]
    rvd_m_o,    & ! = r_v/r_d - 1.0,
    rdv,        & ! = r_d/r_v,
    O_m_rdv,    & ! = 1. - rdv
    r_earth,    & ! mean radius of the earth
    b1,         & !  a
    b2_w,       & !  b
    b3,         & !  c/b (0 degree Celsius [Kelvin])
    b4_w,       & !  d
    fcb,        & ! field capacity of COSMO soil types
    fcb_ec,     & ! field capacity of IFS soil types with CY32r3
    fcb_ec_1s,  & ! field capacity of IFS 1 soil type
    pwpb,       & ! permanent wilting point of COSMO soil types
    pwpb_ec,    & ! permanent wilting point of IFS soil types with CY32r3
    pwpb_ec_1s, & ! permanent wilting point of IFS 1 soil type
    porb,       & ! pore volume of soil types
    porb_ec_1s    ! pore volume of IFS 1 soil type

!-------------------------------------------------------------------------------

USE environment,        ONLY: model_abort
USE interp_utilities,   ONLY: interp_l, interp_q, interp_q_lm, interp_q_bs
USE meteo_utilities,    ONLY: qsat, psat_w, calps
USE utilities,          ONLY: uv2uvrot_vec, uvrot2uv_vec, &
                              uv2uvrot, uvrot2uv

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE org_coarse_interpol

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

! Local variables

INTEGER (KIND=iintegers)    ::      &
 n, k, izdebug, izrank, izlevels

REAL(KIND=ireals), POINTER  :: zp2(:,:), zp3(:,:,:)

LOGICAL                     ::  lzreadin

CHARACTER (LEN =  3)        ::  yzdattyp
CHARACTER (LEN = 10)        ::  yzshname
!
!- End of header
!==============================================================================

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

var_loop: DO n = 1, nvar_in

  ! Set the necessary components from var_in to local variables
  yzshname   =  var_in(n)%name
  yzdattyp   =  var_in(n)%dattyp
  izrank     =  var_in(n)%rank
  izlevels   =  var_in(n)%nlevels
  lzreadin   =  var_in(n)%lreadin
  zp2        => var_in(n)%p2
  zp3        => var_in(n)%p3

  ! Check, whether this variable has to be interpolated; this is done by
  ! checking the data type in the control structure:
  IF (lcomp_bound) THEN
    IF ( (yzdattyp(2:2) /= 'B') .AND. (yzdattyp(3:3) /= 'O') ) THEN
      ! this variable is not needed for the boundary fields
      IF (izdebug > 10) THEN
        PRINT *, '  variable ',n, yzshname, ' not needed for boundary fields  ', nvar_in
      ENDIF
      CYCLE var_loop
    ENDIF
  ELSE
    IF ( (yzdattyp(1:1) /= 'I') .AND. (yzdattyp(3:3) /= 'O') ) THEN
      ! this variable is not needed for the initial fields
      IF (izdebug > 10) THEN
        PRINT *, '  variable ',n, yzshname, ' not needed for initial fields  ', nvar_in
      ENDIF
      CYCLE var_loop
    ENDIF
  ENDIF

  ! Postpone the interpolation of special fields
  IF (yzshname == 'U        '.OR. yzshname == 'V        '.OR. &
      yzshname == 'T_SKIN   '.OR.                             &
      yzshname == 'T_G1     '.OR. yzshname == 'T_G2     '.OR. &
      yzshname == 'T_G3     '.OR. yzshname == 'T_CL     '.OR. &
      yzshname == 'W_G1     '.OR. yzshname == 'W_G2     '.OR. &
      yzshname == 'W_G3     '.OR. yzshname == 'W_CL     ') THEN
    CYCLE var_loop
  ENDIF

  IF (lec2lm .AND. yzshname == 'T_SNOW   ') THEN
    ! also postpone interpolation of T_SNOW from ECMWF
    IF (izdebug > 10) THEN
      PRINT *, '  variable ',n, yzshname, ' is postponed for lec2lm'
    ENDIF
    CYCLE var_loop
  ENDIF

  IF (lcm2lm .AND. yzshname == 'W_SO_REL ') THEN
    ! interpolation of W_SO_REL in case of lcm2lm not required,
    ! result is not used, but overwritten
    CYCLE var_loop
  ENDIF

  ! interpolation of SST in case of lcm2lm not required
  IF (lcm2lm .AND. (yzshname == 'SST      ')) THEN
    CYCLE var_loop
  ENDIF

  IF (llm2lm .AND. lmulti_layer_in) THEN
    ! also the interpolation of t_so and w_so has to be postponed
    IF (yzshname == 'T_SO     ' .OR. yzshname == 'W_SO     ') THEN
      IF (izdebug > 10) THEN
        PRINT *, '  variable ',n, yzshname, ' is postponed for llm2lm'
      ENDIF
      CYCLE var_loop
    ENDIF
  ENDIF

  IF (llm2lm) THEN
    ! In case of lm2lm, only the pressure deviation pp is interpolated
    IF (yzshname == 'P        ') CYCLE var_loop
  ENDIF

  IF (.NOT. lzreadin) THEN
    IF (lcm2lm .AND. (yzshname == 'QC       ')) THEN
      ! reset qc_lm, if it has not been read, but only in case lcm2lm,
      ! where qc is optional
      qc_lm(:,:,:) = 0.0_ireals
      IF (izdebug > 10) THEN
        PRINT *, '  variable ',n, yzshname, ' has been reset !!!'
      ENDIF
    ELSE
      IF (izdebug > 10) THEN
        PRINT *, '  variable ',n, yzshname, ' has not been read !!!'
      ENDIF
    ENDIF
  ELSE
    IF (izrank == 2) THEN
!     IF (izdebug > 10) THEN
!       PRINT *, '  variable ',n, yzshname, ' will be  interpolated'
!     ENDIF
      CALL interpol_coarse(n, 1, zp2, izdebug)
      IF (izdebug > 10) THEN
        PRINT *, '  variable ',n, yzshname, ' has been interpolated  ', &
                               MINVAL (zp2(:,:)), MAXVAL(zp2(:,:))
      ENDIF
    ELSE
      DO k = 1, izlevels
        CALL interpol_coarse(n, k, zp3(:,:,k), izdebug)
        IF (izdebug > 10) THEN
          PRINT *, '  variable ',n, yzshname, 'lev = ', k, ' has been interpolated ', &
                               MINVAL (zp3(:,:,k)), MAXVAL(zp3(:,:,k))
        ENDIF
      ENDDO
    ENDIF
  ENDIF
ENDDO var_loop

IF (llm2lm .OR. lum2lm) THEN
  CALL interpol_coarse_uv_lm (izdebug)
ELSE IF (lec2lm .OR. lcm2lm .OR. lgsm2lm .OR. lgfs2lm .OR. lhir2lm) THEN
  CALL interpol_coarse_uv (izdebug)
ENDIF

IF      (llm2lm) THEN
  CALL interpol_coarse_special_lm (izdebug)
ELSE IF (lec2lm) THEN
  CALL interpol_coarse_special_ec (izdebug)
ELSE IF (lgfs2lm) THEN
  CALL interpol_coarse_special_gfs(izdebug)
ELSE IF (lhir2lm) THEN
  CALL interpol_coarse_special_hir(izdebug)
ELSE IF (lcm2lm) THEN
  CALL interpol_coarse_special_cm (izdebug)
ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE org_coarse_interpol

!==============================================================================
!+ Interpolates a field from EC/LM to LM-subdomains
!------------------------------------------------------------------------------

SUBROUTINE interpol_coarse (mloc_in, mlev, field_in, idbg)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
REAL     (KIND=ireals)   ,  INTENT(IN)     ::  &
  field_in (ie_in, je_in)
  ! record to be interpolated

INTEGER  (KIND=iintegers),  INTENT(IN)     ::  &
  idbg,                 & ! debug level output
  mloc_in,              & ! location of variable in input variable table
  mlev                    ! level of a multi-level field

!------------------------------------------------------------------------------

! Local variables:
INTEGER  (KIND=iintegers)  ::  &
  izerror,              & ! status and error status variable
  mloc_lm,              & ! location of variable in LM variable table
  izlen, i, j

REAL    (KIND=ireals)      ::  &
  field_lm(ie2lm,je2lm)

LOGICAL                    ::  &
  lzmono, lzposdef

CHARACTER (LEN=  1)        ::  &
  yzitype     ! interpolation type

CHARACTER (LEN= 10)        ::  &
  yzinname, & ! name of input field
  yzfuname, & ! name of input field
  yzlmname, & ! name of lm field
  yzname      ! for testing in the IF-construction

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  yzlmname   = '          '
  yzinname   = '          '
  yzfuname   = '          '
  yzname     = '          '
  izerror    = 0_iintegers
  yzroutine  = 'interpol_coarse'

  ! Set the necessary components from var_in to local variables
  yzitype    = var_in(mloc_in)%ipc(1:1)
  IF (var_in(mloc_in)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mloc_in)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF
  izlen      = LEN_TRIM(var_in(mloc_in)%name)
  yzinname   = var_in(mloc_in)%name(1:izlen)
  yzfuname   = var_in(mloc_in)%name

  ! Look for location of variable in LM variable table
  mloc_lm = 0
  DO i = 1, nvar_lm
    yzlmname(1:izlen) = var_lm(i)%name(1:izlen)
    yzname            = yzlmname(1:izlen)
    IF ( (yzlmname(1:izlen) == yzinname(1:izlen)) .AND.                &
         (yzname            == var_lm(i)%name) ) THEN
      mloc_lm = i
      EXIT
    ENDIF
  ENDDO

  IF ( (mloc_lm == 0)          .AND. (yzinname /= 'FI') .AND.  &
      (yzinname /= 'HSURF_in') .AND. (yzinname /= 'QCI').AND.  &
      (yzinname /= 'VW_SO')    .AND. (yzinname /= 'SST').AND.  &
      (yzinname /= 'T_2M')     .AND. (yzinname /= 'QV_2M')) THEN
    yzerrmsg = 'Variable '//yzinname//' not found in LM variable table after interpolation'
    CALL model_abort(my_cart_id, 9011, yzerrmsg, yzroutine)
  ELSE
    IF (idbg > 20) THEN
      IF ((yzinname == 'FI') .OR. (yzinname == 'T_2M') .OR. (yzinname == 'QV_2M')) THEN
        PRINT *, '  variable ',mloc_in, yzinname, 'lev = ', mlev,  &
                 ' needed only as intermediate field'
      ELSE
        PRINT *, '  variable ',mloc_in, yzinname, 'lev = ', mlev,  &
                 ' will be written to COSMO field ', mloc_lm, var_lm(mloc_lm)%name
      ENDIF
    ENDIF
  ENDIF

  ! initialize grdpt_rel_in  (for Cressman Scheme)
  grdpt_rel_in=0

!------------------------------------------------------------------------------
! Section 2:
!------------------------------------------------------------------------------

  SELECT CASE (yzfuname)
  CASE ('T_S      ')
    ! Interpolate t_s_in to t_s_gl
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
    CALL interp_l(field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),    &
                  lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
                  lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
                  t_s_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,      &
                  latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                  grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                  yzerrmsg, izerror)
    IF (idbg > 20) THEN
       PRINT *, '     after interpolation of T_S: t_s_gl:  ', MINVAL(t_s_gl(:,:)), MAXVAL(t_s_gl(:,:))
    ENDIF
  CASE ('T_2M     ')
    ! Interpolate t_2m_in to t_2m_gl
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
    CALL interp_l(field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),    &
                  lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
                  lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
                  t_2m_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,     &
                  latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                  grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                  yzerrmsg, izerror)
  CASE ('QV_2M     ')
    ! Interpolate qv_2m_in to qv_2m_gl
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
    CALL interp_l(field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),    &
                  lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
                  lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
                  qv_2m_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,    &
                  latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                  grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                  yzerrmsg, izerror)
  CASE ('FI       ')
    ! Interpolate control level to fic_gl
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
    CALL interp_l(field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),    &
                  lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
                  lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
                  fic_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,      &
                  latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                  grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                  yzerrmsg, izerror)
  CASE ('HSURF_in ')
    ! The input orography has already been interpolated in external_data and
    ! transfered to the hhl_gl variable.
    ! Check that the orography in initial fields is the same as in external
    ! parameter input file inext_lfn.
    ! Only extremas are checked.
    IF (MINVAL(ABS(field_in-hsurf_in)) > 0.001 .OR. &
        MAXVAL(ABS(field_in-hsurf_in)) > 0.001) THEN
      izerror = 2019
      yzerrmsg = 'Initial fields orography is not identical to ext. para. input'
    ENDIF
  CASE ('PS       ')
    ! Interpolate surface pressure to ps_gl
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
    CALL interp_l(field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),    &
                  lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
                  lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
                  ps_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,       &
                  latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                  grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                  yzerrmsg, izerror)
  CASE DEFAULT
    IF (yzitype == 'Q') THEN
      ! quadratic interpolation to LM-field
      IF (lec2lm .OR. lcm2lm .OR. lgsm2lm .OR. lgfs2lm .OR. lhir2lm) THEN
        ! Introduce option of bicubic spline-interpolation instead of quadratic
        IF (l_bicub_spl) THEN
          CALL interp_q_bs                                                   &
             (field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),        &
              lzmono, lzposdef, lbd_frame_cur, lolp_in, lolp_lm,             &
              undef, lmask_lm, x_wght(:,:,1), y_wght(:,:,1),                 &
              field_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,        &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,       &
              grdpt_rel_in, ie_in_tot, je_in_tot, lcm2lm, l_cressman,        &
              yzerrmsg, izerror)
        ELSE
          CALL interp_q(field_in, ie_in,je_in, i_index(:,:,1),j_index(:,:,1),&
                        lzmono, lzposdef, lbd_frame_cur, undef, lmask_lm,    &
                        x_wght(:,:,1), y_wght(:,:,1),                        &
                        field_lm, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
        ENDIF
      ELSEIF (llm2lm .OR. lum2lm) THEN
        CALL interp_q_lm(field_in, ie_in, je_in,                      &
                    field_lm, i_index(:,:,1), j_index(:,:,1),         &
                    x_wght(:,:,1), y_wght(:,:,1), ie2lm, je2lm,       &
                    lzmono, lzposdef, yzerrmsg, izerror)
      ENDIF
    ELSE
      ! normal linear (or match) interpolation to LM-field
      IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
      CALL interp_l(field_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
                    lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,      &
                    lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
                    field_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in, &
                    latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,&
                    grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,         &
                    yzerrmsg, izerror)
    ENDIF

    SELECT CASE (var_lm(mloc_lm)%rank)
    CASE (3)
      var_lm(mloc_lm)%p3(1:ie2lm,1:je2lm,mlev) = field_lm(1:ie2lm,1:je2lm)
    CASE (2)
      var_lm(mloc_lm)%p2(1:ie2lm,1:je2lm)      = field_lm(1:ie2lm,1:je2lm)
    END SELECT

  END SELECT ! in name

  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, 1, yzerrmsg, yzroutine)
  ENDIF


!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE interpol_coarse

!==============================================================================
!+ Interpolates a u and v from EC/LM to LM-subdomains
!------------------------------------------------------------------------------

SUBROUTINE interpol_coarse_uv (idbg)

!------------------------------------------------------------------------------
!
! Description:
!   The required local fields u and v from GME are available on each PE,
!   so that they can be interpolated to the LM subdomain.
!   With mloc_in (location in the GME
!   variable table) the interpolation code is checked and the fields are
!   interpolated using the interpolation utilities. The result is written
!   to the corresponding LM-fields.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER  (KIND=iintegers), INTENT(IN)  ::  &
  idbg                    ! debug level output

! Local variables:
INTEGER  (KIND=iintegers)  ::  &
  izerror,              & ! status and error status variable
  k

REAL    (KIND=ireals)      ::  &
  zfield_lm(ie2lm,je2lm)

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  IF (idbg > 10) THEN
    PRINT *, ' Start interpolation of U/V'
  ENDIF

  yzroutine  = 'interpol_coarse_uv'

  ! Initialization of grdpt_rel_in
  grdpt_rel_in = 0

!------------------------------------------------------------------------------
! Section 2:
!------------------------------------------------------------------------------

  DO k=1,ke_in
    IF (l_bicub_spl) THEN
      ! Introduce option of bicubic spline-interpolation instead of quadratic
      CALL interp_q_bs                                                    &
         (u_in(:,:,k), ie_in, je_in, i_index(:,:,2), j_index(:,:,2),      &
          .FALSE.,.FALSE., lbd_frame_cur, lolp_in, lolp_lm,               &
          undef, lmask_lm, x_wght(:,:,2), y_wght(:,:,2),                  &
          u_lm(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,      &
          latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,        &
          grdpt_rel_in, ie_in_tot, je_in_tot, lcm2lm, l_cressman,         &
          yzerrmsg, izerror)
      IF (pollat_in /= pollat .OR. pollon_in /= pollon) THEN
        ! I must rotate => the other component is needed
        CALL interp_q_bs                                                  &
           (v_in(:,:,k), ie_in, je_in, i_index(:,:,3), j_index(:,:,3),    &
            .FALSE.,.FALSE., lbd_frame_cur, lolp_in, lolp_lm,             &
            undef, lmask_lm, x_wght(:,:,3), y_wght(:,:,3),                &
            zfield_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,      &
            latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,      &
            grdpt_rel_in, ie_in_tot, je_in_tot, lcm2lm, l_cressman,       &
            yzerrmsg, izerror)
      ENDIF
    ELSE
      ! quadratic interpolation to LM u-points
      CALL interp_q (u_in(:,:,k), ie_in,je_in, i_index(:,:,2), j_index(:,:,2),&
                     .FALSE., .FALSE., lbd_frame_cur, undef, lmask_lm,        &
                     x_wght(:,:,2), y_wght(:,:,2),                            &
                     u_lm(:,:,k), 1, ie2lm, 1, je2lm, yzerrmsg, izerror)

      IF (pollat_in /= pollat .OR. pollon_in /= pollon) THEN
        ! I must rotate => the other component is needed
        CALL interp_q (v_in(:,:,k),ie_in,je_in, i_index(:,:,3),j_index(:,:,3),&
                       .FALSE., .FALSE., lbd_frame_cur, undef, lmask_lm,      &
                       x_wght(:,:,3), y_wght(:,:,3),                          &
                       zfield_lm, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
      ENDIF
    ENDIF

    IF ((pollat_in /= 90.0_ireals .OR. pollon_in /= 180.0_ireals) .AND. &
     (pollat_in /= pollat .OR. pollon_in /= pollon)) THEN
      ! Rotate wind back to geographic system if not yet there and
      ! a change of rotation is needed

      CALL uvrot2uv_vec (u_lm(:,:,k), zfield_lm(:,:), latlm_u(:,:),        &
                         lonlm_u(:,:), pollat_in, pollon_in, ie2lm, je2lm)
    ENDIF

    IF (pollat_in /= pollat .OR. pollon_in /= pollon) THEN
      ! Rotate wind to rotated system if a change of rotation is needed
      CALL uv2uvrot_vec (u_lm(:,:,k), zfield_lm(:,:), latlm_u(:,:),        &
                         lonlm_u(:,:), pollat, pollon, ie2lm, je2lm)
    ENDIF


    IF (l_bicub_spl) THEN
      ! Introduce option of bicubic spline-interpolation instead of quadratic
      CALL interp_q_bs                                                    &
         (v_in(:,:,k), ie_in, je_in, i_index(:,:,4), j_index(:,:,4),      &
          .FALSE.,.FALSE., lbd_frame_cur, lolp_in, lolp_lm,               &
          undef, lmask_lm, x_wght(:,:,4), y_wght(:,:,4),                  &
          v_lm(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,      &
          latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,        &
          grdpt_rel_in, ie_in_tot, je_in_tot, lcm2lm, l_cressman,         &
          yzerrmsg, izerror)
      IF (pollat_in /= pollat .OR. pollon_in /= pollon) THEN
        ! I must rotate => the other component is needed
        CALL interp_q_bs                                                  &
           (u_in(:,:,k), ie_in, je_in, i_index(:,:,5), j_index(:,:,5),    &
            .FALSE.,.FALSE., lbd_frame_cur, lolp_in, lolp_lm,             &
            undef, lmask_lm, x_wght(:,:,5), y_wght(:,:,5),                &
            zfield_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,      &
            latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,      &
            grdpt_rel_in, ie_in_tot, je_in_tot, lcm2lm, l_cressman,       &
            yzerrmsg, izerror)

      ENDIF
    ELSE
      ! quadratic interpolation to LM v-points
      CALL interp_q (v_in(:,:,k), ie_in, je_in, i_index(:,:,4), j_index(:,:,4),&
                     .FALSE., .FALSE., lbd_frame_cur, undef, lmask_lm,         &
                     x_wght(:,:,4), y_wght(:,:,4),                             &
                     v_lm(:,:,k), 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
      IF (pollat_in /= pollat .OR. pollon_in /= pollon) THEN
        ! I must rotate => the other component is needed
        CALL interp_q (u_in(:,:,k), ie_in, je_in, i_index(:,:,5),j_index(:,:,5),&
                       .FALSE., .FALSE., lbd_frame_cur, undef, lmask_lm,        &
                       x_wght(:,:,5), y_wght(:,:,5),                            &
                       zfield_lm, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
      ENDIF
    ENDIF

    IF ((pollat_in /= 90.0_ireals .OR. pollon_in /= 180.0_ireals) .AND. &
     (pollat_in /= pollat .OR. pollon_in /= pollon)) THEN
      ! Rotate wind back to geographic system if not yet there and
      ! a change of rotation is needed
      CALL uvrot2uv_vec (zfield_lm(:,:), v_lm(:,:,k), latlm_v(:,:),        &
                         lonlm_v(:,:), pollat_in, pollon_in, ie2lm, je2lm)
    ENDIF

    IF (pollat_in /= pollat .OR. pollon_in /= pollon) THEN
      ! Rotate wind to rotated system if a change of rotation is needed
      CALL uv2uvrot_vec (zfield_lm(:,:), v_lm(:,:,k), latlm_v(:,:),        &
                         lonlm_v(:,:), pollat, pollon, ie2lm, je2lm)
    ENDIF

  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE interpol_coarse_uv

!==============================================================================
!+ Interpolates a u and v from LM to LM-subdomains
!------------------------------------------------------------------------------

SUBROUTINE interpol_coarse_uv_lm (idbg)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER  (KIND=iintegers), INTENT(IN)  ::  &
  idbg                    ! debug level output

! Local variables:
INTEGER  (KIND=iintegers)  ::  &
  izerror,              & ! status and error status variable
  k

REAL    (KIND=ireals)      ::  &
  zfield_lm(ie2lm,je2lm)

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  yzroutine  = 'interpol_coarse_uv_lm'

!------------------------------------------------------------------------------
! Section 2:
!------------------------------------------------------------------------------

  DO k=1,ke_in
    ! quadratic interpolation to LM u-points
    CALL interp_q_lm(u_in(:,:,k), ie_in, je_in,                              &
                     u_lm(:,:,k), i_index(:,:,2), j_index(:,:,2),            &
                     x_wght(:,:,2), y_wght(:,:,2), ie2lm, je2lm,             &
                    .FALSE., .FALSE., yzerrmsg, izerror)

    IF (                                                                     &
      ! COSMO grids have different rotation
      (llm2lm .AND. (pollat_in   /= pollat .OR. pollon_in    /= pollon))     &
      .OR.                                                                   &
      ! UM grid and fine COSMO grid have different rotation
      (lum2lm .AND. (90.0_ireals /= pollat .OR. 180.0_ireals /= pollon)) ) THEN

      ! I must rotate => the other component is needed
      CALL interp_q_lm(v_in(:,:,k), ie_in, je_in,                            &
                       zfield_lm, i_index(:,:,3), j_index(:,:,3),            &
                       x_wght(:,:,3), y_wght(:,:,3), ie2lm, je2lm,           &
                      .FALSE., .FALSE., yzerrmsg, izerror)
    ENDIF

    IF (.NOT. lum2lm) THEN
      IF ((pollat_in /= 90.0_ireals .OR. pollon_in /= 180.0_ireals) .AND.    &
          (pollat_in /= pollat      .OR. pollon_in /= pollon)) THEN
        ! Rotate wind back to geographic system if not yet there and
        ! a change of rotation is needed
        CALL uvrot2uv_vec (u_lm(:,:,k), zfield_lm(:,:), latlm_u(:,:),        &
                           lonlm_u(:,:), pollat_in, pollon_in, ie2lm, je2lm)
      ENDIF
    ENDIF

    IF (                                                                     &
      ! COSMO grids have different rotation
      (llm2lm .AND. (pollat_in   /= pollat .OR. pollon_in    /= pollon))     &
      .OR.                                                                   &
      ! UM grid and fine COSMO grid have different rotation
      (lum2lm .AND. (90.0_ireals /= pollat .OR. 180.0_ireals /= pollon)) ) THEN

      ! Rotate wind to rotated system if a change of rotation is needed
      CALL uv2uvrot_vec (u_lm(:,:,k), zfield_lm(:,:), latlm_u(:,:),          &
                         lonlm_u(:,:), pollat, pollon, ie2lm, je2lm)
    ENDIF

    IF (idbg > 10) THEN
      PRINT *, '  variable       U         lev = ', k, ' has been interpolated'
    ENDIF

    ! quadratic interpolation to LM v-points
    CALL interp_q_lm(v_in(:,:,k), ie_in, je_in,                              &
                     v_lm(:,:,k), i_index(:,:,4), j_index(:,:,4),            &
                     x_wght(:,:,4), y_wght(:,:,4), ie2lm, je2lm,             &
                    .FALSE., .FALSE., yzerrmsg, izerror)

    IF (                                                                     &
      ! COSMO grids have different rotation
      (llm2lm .AND. (pollat_in   /= pollat .OR. pollon_in    /= pollon))     &
      .OR.                                                                   &
      ! UM grid and fine COSMO grid have different rotation
      (lum2lm .AND. (90.0_ireals /= pollat .OR. 180.0_ireals /= pollon)) ) THEN

      ! I must rotate => the other component is needed
      CALL interp_q_lm(u_in(:,:,k), ie_in, je_in,                            &
                       zfield_lm, i_index(:,:,5), j_index(:,:,5),            &
                       x_wght(:,:,5), y_wght(:,:,5), ie2lm, je2lm,           &
                      .FALSE., .FALSE., yzerrmsg, izerror)
    ENDIF

    IF (.NOT. lum2lm) THEN
      IF ((pollat_in /= 90.0_ireals .OR. pollon_in /= 180.0_ireals) .AND.    &
       (pollat_in /= pollat .OR. pollon_in /= pollon)) THEN
        ! Rotate wind back to geographic system if not yet there and
        ! a change of rotation is needed
        CALL uvrot2uv_vec (zfield_lm(:,:), v_lm(:,:,k), latlm_v(:,:),        &
                           lonlm_v(:,:), pollat_in, pollon_in, ie2lm, je2lm)
      ENDIF
    ENDIF

    IF (                                                                     &
      ! COSMO grids have different rotation
      (llm2lm .AND. (pollat_in   /= pollat .OR. pollon_in    /= pollon))     &
      .OR.                                                                   &
      ! UM grid and fine COSMO grid have different rotation
      (lum2lm .AND. (90.0_ireals /= pollat .OR. 180.0_ireals /= pollon)) ) THEN

      ! Rotate wind to rotated system if a change of rotation is needed
      CALL uv2uvrot_vec (zfield_lm(:,:), v_lm(:,:,k), latlm_v(:,:),          &
                         lonlm_v(:,:), pollat, pollon, ie2lm, je2lm)
    ENDIF

    IF (idbg > 10) THEN
      PRINT *, '  variable       V         lev = ', k, ' has been interpolated'
    ENDIF

  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE interpol_coarse_uv_lm

!==============================================================================

SUBROUTINE interpol_coarse_special_lm (idbg)

!------------------------------------------------------------------------------
!
! Description:
!   Interpolation (horizontal and vertical) of ground fields and
!   computation of nonconstant fields not available in LM COARSE input
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(IN) ::  &
  idbg                    ! debug level output

! Local variables

INTEGER (KIND=iintegers), PARAMETER :: nl_soil_in_max=4

INTEGER (KIND=iintegers) :: i,j,k,l,m,    & ! Loop indices
 izerror,                      & ! status and error status variable
 s_t,                          & ! soiltype of input
 i1,i2,                        & ! Interpolation indices
 groundl_in(nl_soil_in_max),   & ! Location of multi-level ground
 groundl_lm(nl_soil_in_max),   & !   variables in variable tables
 mzdpsdt_loc_in, mzts_loc_in,  & ! Locations of
 mztm_loc_in,                  & !   input variables
 mztso_loc_in, mzwso_loc_in,   & !   input variables
 mztcl_loc_in, mztsnow_loc_in, & !   in input
 mzqvs_loc_in, mzwg1_loc_in,   & !   variable table
 mzwg2_loc_in, mzwg3_loc_in,   & !
 mzwcl_loc_in,                 & !
 mzqi_loc_in,                  & ! ----------------
 mzdpsdt_loc_lm,               & ! Locations of
 mztm_loc_lm, mztcl_loc_lm,    & !   output variables
 mztsnow_loc_lm,               & !   in output
 mzwg1_loc_lm, mzwg2_loc_lm,   & !   variable table
 mzwg3_loc_lm, mzwcl_loc_lm      ! ----------------

REAL (KIND=ireals) :: wei,         & ! Interpolaton weights
  qrs   (ie2lm,je2lm)            , & ! precipitation water (water loading)
  zdt_so(ie_in,je_in)            , & ! base-st. pres. thickness of surface layer
 zx_g_gl(ie2lm,je2lm,nl_soil_in_max) ! Automatic temporary array

LOGICAL                    ::  &
  lzmono, lzposdef

CHARACTER (LEN=  1)        ::  &
  yzitype     ! interpolation type

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Look for locations of special variables in the variable tables
!------------------------------------------------------------------------------

yzroutine = 'interpol_coarse_special_lm'

! Initialization of grdpt_rel_in
grdpt_rel_in = 0

DO i = 1, nvar_in
  SELECT CASE (var_in(i)%name)
  CASE ('DPSDT     '); mzdpsdt_loc_in  = i
  CASE ('T_S       '); mzts_loc_in     = i
  CASE ('T_M       '); mztm_loc_in     = i
  CASE ('T_CL      '); mztcl_loc_in    = i
  CASE ('T_SNOW    '); mztsnow_loc_in  = i
  CASE ('QV_S      '); mzqvs_loc_in    = i
  CASE ('W_G1      '); mzwg1_loc_in    = i
  CASE ('W_G2      '); mzwg2_loc_in    = i
  CASE ('W_G3      '); mzwg3_loc_in    = i
  CASE ('W_CL      '); mzwcl_loc_in    = i
  CASE ('QI        '); mzqi_loc_in     = i
  CASE ('T_SO      '); mztso_loc_in    = i
  CASE ('W_SO      '); mzwso_loc_in    = i
  END SELECT
ENDDO

DO i = 1, nvar_lm
  SELECT CASE (var_lm(i)%name)
  CASE ('DPSDT     '); mzdpsdt_loc_lm  = i
  CASE ('T_M       '); mztm_loc_lm     = i
  CASE ('T_CL      '); mztcl_loc_lm    = i
  CASE ('T_SNOW    '); mztsnow_loc_lm  = i
  CASE ('W_G1      '); mzwg1_loc_lm    = i
  CASE ('W_G2      '); mzwg2_loc_lm    = i
  CASE ('W_G3      '); mzwg3_loc_lm    = i
  CASE ('W_CL      '); mzwcl_loc_lm    = i
  END SELECT
ENDDO

! Computation of surface pressure
qrs(:,:) = 0.0  !  no water loading
ps_gl(:,:) = ( p0_gl(:,:,ke_in) + pp_lm(:,:,ke_in) ) * EXP  (  0.5*dp0_gl(:,:,ke_in) / &
   ( t_lm(:,:,ke_in)*(1.0+rvd_m_o*qv_lm(:,:,ke_in)-qc_lm(:,:,ke_in)-qrs(:,:)) &
                                          *r_d*rho0_gl(:,:,ke_in) )  )

IF (.NOT. lmulti_layer_in) THEN

  !----------------------------------------------------------------------------
  ! Section 2: Soil moisture interpolation
  !----------------------------------------------------------------------------

  ! Use of indirection in order to treat multi-level ground fields as a
  ! 1d array of 2d arrays
  groundl_in(1)=mzwg1_loc_in
  groundl_in(2)=mzwg2_loc_in
  groundl_in(3)=mzwg3_loc_in
  groundl_in(4)=mzwcl_loc_in
  IF (nl_soil_in == 2) groundl_in(3) = mzwcl_loc_in

  groundl_lm(1)=mzwg1_loc_lm
  groundl_lm(2)=mzwg2_loc_lm
  groundl_lm(3)=mzwg3_loc_lm
  groundl_lm(4)=mzwcl_loc_lm
  IF (nl_soil_lm == 2) groundl_lm(3) = mzwcl_loc_lm

  ! Scale the soil moisture by the pore volume of soil
  ! for soiltypes #3 to #8; set soil moisture to 0 otherwise

! gdm 20080828 : more infos
  IF (idbg > 5) THEN
    WRITE (*,'(A)')    yzroutine//': NEW soil moisture interpolation with SMI'
    WRITE (*,'(A,I3)') yzroutine//': FOR .NOT. lmulti_layer_in WITH '// &
                                                'nl_soil_in=', nl_soil_in
    WRITE (*,'(A,L3,A)') yzroutine//': IF: l_smi = ', l_smi, &
                             ' is .TRUE. THEN w_so_lm >= pwpb'
  ENDIF

  DO k = 1, nl_soil_in + 1
    IF (l_smi) THEN
      DO i2 = 1, je_in
        DO i1 = 1, ie_in
! gdm 20071030 : new scaling of 2-3 layer soil moisture
          s_t = NINT(soiltyp_in(i1,i2))
          IF (s_t >= 3 .AND. s_t <= 8 .AND.           &
              var_in(groundl_in(k))%p2(i1,i2) /= undef) THEN
            ! Scale the soil moisture with the soil moisture index SMI:
            ! smi = (sm-pwp)/(fc-pwp)
            ! Caution: sm is in volumetric units, soil depth in cm!
            var_in(groundl_in(k))%p2(i1,i2) =                                 &
                                           var_in(groundl_in(k))%p2(i1,i2) /  &
            ((var_in(groundl_in(k))%levbot-var_in(groundl_in(k))%levtop)*0.01)
            var_in(groundl_in(k))%p2(i1,i2) = MAX (                           &
                                          (var_in(groundl_in(k))%p2(i1,i2) -  &
                                 pwpb(s_t)) / (fcb(s_t)-pwpb(s_t)), 0.0_ireals)
          ELSE
            var_in(groundl_in(k))%p2(i1,i2) = 0.0_ireals
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO i2 = 1, je_in
        DO i1 = 1, ie_in
          IF (soiltyp_in(i1,i2) >= 3 .AND. soiltyp_in(i1,i2) <= 8 .AND.  &
              var_in(groundl_in(k))%p2(i1,i2) /= undef) THEN
            var_in(groundl_in(k))%p2(i1,i2) = &
            var_in(groundl_in(k))%p2(i1,i2) / porb(NINT(soiltyp_in(i1,i2)))
          ELSE
            var_in(groundl_in(k))%p2(i1,i2) = 0.0_ireals
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! Interpolate horizontally and store in temporary arrays
    IF (var_in(groundl_in(k))%ipc(2:2) == 'T') THEN
      lzmono   = .TRUE.
    ELSE
      lzmono   = .FALSE.
    ENDIF
    IF (var_in(groundl_in(k))%ipc(3:3) == 'T') THEN
      lzposdef = .TRUE.
    ELSE
      lzposdef = .FALSE.
    ENDIF
    yzitype    = var_in(groundl_in(k))%ipc(1:1)
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
    CALL interp_l(var_in(groundl_in(k))%p2,                               &
                  ie_in, je_in, i_index(:,:,1), j_index(:,:,1),           &
                  lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,      &
                  lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
                  zx_g_gl(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                  latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,&
                  grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,         &
                  yzerrmsg, izerror)
  ENDDO

  ! Construct w_g*:
  IF ((nl_soil_in == 3) .AND. (nl_soil_lm == 2) ) THEN
    zx_g_gl(:,:,1) = zx_g_gl(:,:,1) + zx_g_gl(:,:,2)
    zx_g_gl(:,:,2) = zx_g_gl(:,:,3)
    zx_g_gl(:,:,3) = zx_g_gl(:,:,4)
  ENDIF
  IF ((nl_soil_in == 2) .AND. (nl_soil_lm == 3) ) THEN
    ! Scale w_g1_lm and w_g2_lm with thickness of soil moisture layers,
    ! if going from 2 COARSE-layers to 3 LM-layers
    zx_g_gl(:,:,4) = zx_g_gl(:,:,3)
    zx_g_gl(:,:,3) = zx_g_gl(:,:,2)
    wei = (var_lm(groundl_lm(2))%levbot - var_lm(groundl_lm(2))%levtop) / &
          (var_in(groundl_in(1))%levbot - var_in(groundl_in(1))%levtop)
    zx_g_gl(:,:,2) = zx_g_gl(:,:,1) * wei
    wei = (var_lm(groundl_lm(1))%levbot - var_lm(groundl_lm(1))%levtop) / &
          (var_in(groundl_in(1))%levbot - var_in(groundl_in(1))%levtop)
    zx_g_gl(:,:,1) = zx_g_gl(:,:,1) * wei
  ENDIF

  DO l = 1, nl_soil_lm + 1
    var_lm(groundl_lm(l))%p2(:,:) = zx_g_gl(:,:,l)
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 3: Soil temperature interpolation
  !----------------------------------------------------------------------------

  ! All temperatures have already been interpolated in interpol_coarse
  ! except T_CL .  Use multi-level formulation with k=1

  groundl_in(1)=mztcl_loc_in
  groundl_lm(1)=mztcl_loc_lm

  k = 1
    ! Interpolate horizontally and store temporarily in lm arrays
    IF (var_in(groundl_in(k))%ipc(2:2) == 'T') THEN
      lzmono   = .TRUE.
    ELSE
      lzmono   = .FALSE.
    ENDIF
    IF (var_in(groundl_in(k))%ipc(3:3) == 'T') THEN
      lzposdef = .TRUE.
    ELSE
      lzposdef = .FALSE.
    ENDIF
    yzitype    = var_in(groundl_in(k))%ipc(1:1)
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1

    IF (var_in(groundl_in(k))%lreadin) THEN
      CALL interp_l(var_in(groundl_in(k))%p2,                             &
                  ie_in, je_in, i_index(:,:,1), j_index(:,:,1),           &
                  lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,      &
                  lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
                  zx_g_gl(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                  latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,&
                  grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,         &
                  yzerrmsg, izerror)
      var_lm(groundl_lm(k))%p2(:,:)=zx_g_gl(:,:,k)
    ELSE
      ! For boundary fields var_in(tcl) not available and not needed
      var_lm(groundl_lm(k))%p2(:,:)=0.0_ireals
    ENDIF

  ! Reassign the variables for further processing
  dtms_gl(:,:) = var_lm(mztm_loc_lm)%p2(:,:) - t_s_gl(:,:)

ELSE

  IF (.NOT. lcomp_bound) THEN

    !--------------------------------------------------------------------------
    ! Section 2: Soil moisture interpolation
    !--------------------------------------------------------------------------

    ! Interpolate horizontally and store in w_so_lm
    IF (var_in(mzwso_loc_in)%ipc(2:2) == 'T') THEN
      lzmono   = .TRUE.
    ELSE
      lzmono   = .FALSE.
    ENDIF
    IF (var_in(mzwso_loc_in)%ipc(3:3) == 'T') THEN
      lzposdef = .TRUE.
    ELSE
      lzposdef = .FALSE.
    ENDIF
    yzitype    = var_in(mzwso_loc_in)%ipc(1:1)
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1

! gdm 20080828 : more infos
    IF (idbg > 5) THEN
      WRITE (*,'(A)')    yzroutine//': NEW soil moisture interpolation with SMI'
      WRITE (*,'(A,I3)') yzroutine//': FOR lmulti_layer_in WITH ke_soil_lm=',  &
                                                                      ke_soil_lm
      WRITE (*,'(A,L3,A)') yzroutine//': IF: l_smi = ', l_smi,  &
                                ' is .TRUE. THEN w_so_lm >= pwpb'
    ENDIF

    DO k = 1, ke_soil_lm + 1
      IF (l_smi) THEN
        ! Scale the soil moisture with the soil moisture index SMI:
        ! smi = (sm-pwp)/(fc-pwp)
        ! Caution: sm is in volumetric units, soil depth in cm!
        DO i2 = 1, je_in
          DO i1 = 1, ie_in
            s_t = NINT(soiltyp_in(i1,i2))
            IF (s_t >= 3 .AND. s_t <= 8 .AND.            &
                var_in(mzwso_loc_in)%p3(i1,i2,k) /= undef) THEN 
              IF (k == 1) THEN
                var_in(mzwso_loc_in)%p3(i1,i2,k) =    &
                var_in(mzwso_loc_in)%p3(i1,i2,k) / 0.01
              ELSE
                var_in(mzwso_loc_in)%p3(i1,i2,k) =                          &
                var_in(mzwso_loc_in)%p3(i1,i2,k) / ((3**(k-1)-3**(k-2))*0.01)
              ENDIF
              var_in(mzwso_loc_in)%p3(i1,i2,k) = MAX (                        &
                                           (var_in(mzwso_loc_in)%p3(i1,i2,k)- &
                                 pwpb(s_t)) / (fcb(s_t)-pwpb(s_t)), 0.0_ireals)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        ! Scale the soil moisture by the pore volume of soil
        ! for soiltypes #3 to #8; set soil moisture to 0 otherwise
        DO i2 = 1, je_in
          DO i1 = 1, ie_in
            IF (soiltyp_in(i1,i2) >= 3 .AND. soiltyp_in(i1,i2) <= 8 .AND.  &
                var_in(mzwso_loc_in)%p3(i1,i2,k) /= undef) THEN
              var_in(mzwso_loc_in)%p3(i1,i2,k) = &
              var_in(mzwso_loc_in)%p3(i1,i2,k) / porb(NINT(soiltyp_in(i1,i2)))
            ELSE
              var_in(mzwso_loc_in)%p3(i1,i2,k) = 0.0_ireals
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF (idbg > 1) THEN
        PRINT *, "Change input where soiltyp_in(i1,i2) == 1 (ice) for level ", k
      ENDIF
      DO i2 = 1, je_in
        DO i1 = 1, ie_in
          IF (soiltyp_in(i1,i2) <= 1 ) THEN
          ! search for the closest soiltyp_in NOT equal to ICE
            DO i = 1, ie_in
              DO j = 1, je_in
                l = ie_in
                m = je_in
                IF (soiltyp_in(i,j) >= 3 .AND. soiltyp_in(i,j) <= 8 ) THEN
                   l = MIN(l,ABS(i-i1))
                   m = MIN(l,ABS(j-i2))
                ENDIF
              ENDDO
            ENDDO
          ! redefine input with the maximum found in the NOT ICE vicinity
            var_in(mzwso_loc_in)%p3(i1,i2,k) = MAX(                            &
                   var_in(mzwso_loc_in)%p3(MAX(i1-l,    1),    i2,         k), &
                   var_in(mzwso_loc_in)%p3(MAX(i1-l,    1),MAX(i2-m,    1),k), &
                   var_in(mzwso_loc_in)%p3(MAX(i1-l,    1),MIN(i2+m,je_in),k), &
                   var_in(mzwso_loc_in)%p3(MIN(i1+l,ie_in),MAX(i2-m,    1),k), &
                   var_in(mzwso_loc_in)%p3(MIN(i1+l,ie_in),MIN(i2+m,je_in),k), &
                   var_in(mzwso_loc_in)%p3(    i1,         MIN(i2+m,je_in),k), &
                   var_in(mzwso_loc_in)%p3(    i1,         MAX(i2-m,je_in),k), &
                   var_in(mzwso_loc_in)%p3(MIN(i1+l,ie_in),    i2,         k)  )
          ENDIF
        ENDDO
      ENDDO

      CALL interp_l(var_in(mzwso_loc_in)%p3(:,:,k),                         &
                    ie_in, je_in, i_index(:,:,1), j_index(:,:,1),           &
                    lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,      &
                    lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
                    w_so_lm(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                    latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                    grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                    yzerrmsg, izerror)
      IF (idbg > 10) THEN
        PRINT *, '  variable ', var_in(mzwso_loc_in)%name, ' level: ', k, ' has been interpolated  ', &
                               MINVAL (w_so_lm(:,:,k)), MAXVAL(w_so_lm(:,:,k))
      ENDIF
    ENDDO

    !--------------------------------------------------------------------------
    ! Section 3: Soil temperature interpolation
    !--------------------------------------------------------------------------

    ! Interpolate horizontally and store temporarily in dt_so_gl
    IF (var_in(mztso_loc_in)%ipc(2:2) == 'T') THEN
      lzmono   = .TRUE.
    ELSE
      lzmono   = .FALSE.
    ENDIF
    IF (var_in(mztso_loc_in)%ipc(3:3) == 'T') THEN
      lzposdef = .TRUE.
    ELSE
      lzposdef = .FALSE.
    ENDIF
    yzitype    = var_in(mztso_loc_in)%ipc(1:1)
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1

    DO k = 1, ke_soil_lm+1
      DO i2 = 1, je_in
        DO i1 = 1, ie_in
          IF (t_so_in(i1,i2,k) /= undef .AND. t_so_in(i1,i2,k-1) /= undef) THEN
            zdt_so (i1,i2) = t_so_in(i1,i2,k) - t_so_in(i1,i2,k-1)
          ELSE
            zdt_so (i1,i2) = undef
          ENDIF
        ENDDO
      ENDDO

      CALL interp_l(zdt_so,                                               &
                  ie_in, je_in, i_index(:,:,1), j_index(:,:,1),           &
                  lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,      &
                  lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
                  dt_so_gl(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                  latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                  grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                  yzerrmsg, izerror)
      IF (idbg > 10) THEN
        PRINT *, '  variable ', var_in(mztso_loc_in)%name, ' level: ', k, ' has been interpolated  ', &
                               MINVAL(dt_so_gl(:,:,k)), MAXVAL(dt_so_gl(:,:,k))
      ENDIF
    ENDDO

  ENDIF

  ! Interpolate T_SO(mlev=0) to t_s_gl
  ! (in any case: for initial and boundary values)

  IF (var_in(mztso_loc_in)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mztso_loc_in)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF
  yzitype    = var_in(mztso_loc_in)%ipc(1:1)
  IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1

  CALL interp_l(t_so_in(:,:,0),                                          &
              ie_in, je_in, i_index(:,:,1), j_index(:,:,1),              &
              lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
              lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
              t_s_gl(:,:), 1, ie2lm, 1, je2lm, startlat_in, startlon_in, &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
              grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
              yzerrmsg, izerror)
  IF (idbg > 10) THEN
    PRINT *, '  variable ', var_in(mztso_loc_in)%name, ' level: ', 0, ' has been interpolated  ', &
                           MINVAL (t_s_gl(:,:)), MAXVAL(t_s_gl(:,:))
  ENDIF

ENDIF

!------------------------------------------------------------------------------
! Section 4: Set DPSDT to 0.
!------------------------------------------------------------------------------

!IF (.NOT. var_in(mzdpsdt_loc_in)%lreadin) THEN
!  ! If dpsdt is not available, and this is the case, set it to 0.
!  var_lm(mzdpsdt_loc_lm)%p2 = 0.0_ireals
!ENDIF

!------------------------------------------------------------------------------
! Section 5: Set dtkes_gl to T(lowlev) - T_S
!------------------------------------------------------------------------------

CALL interp_l(t_in(:,:,ke_in), ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
              .FALSE., .FALSE., 'L', lbd_frame_cur, lolp_in,                 &
              lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),        &
              dtkes_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,        &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,       &
              grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                &
              yzerrmsg, izerror)
dtkes_gl(:,:) = dtkes_gl(:,:) - t_s_gl(:,:)

!------------------------------------------------------------------------------
! Section 6: Interpolate qv_s and convert in rh_s
!------------------------------------------------------------------------------

IF (.NOT. var_in(mzqvs_loc_in)%lreadin) THEN
  qv_s_in = qv_in(:,:,ke_in)
ENDIF

IF (var_in(mzqvs_loc_in)%ipc(2:2) == 'T') THEN
  lzmono   = .TRUE.
ELSE
  lzmono   = .FALSE.
ENDIF
IF (var_in(mzqvs_loc_in)%ipc(3:3) == 'T') THEN
  lzposdef = .TRUE.
ELSE
  lzposdef = .FALSE.
ENDIF
yzitype    = var_in(mzqvs_loc_in)%ipc(1:1)
IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1

CALL interp_l(qv_s_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),   &
              lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,       &
              lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),  &
              rh_s_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,   &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
              grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,          &
              yzerrmsg, izerror)

DO j = 1, je2lm
  DO i = 1, ie2lm
     rh_s_gl(i,j) = rh_s_gl(i,j) /  &
     qsat(psat_w(t_s_gl(i,j),B1, B2_w, B3, B4_w),   &
     ps_gl(i,j), rdv, O_m_rdv)
  ENDDO
ENDDO

!------------------------------------------------------------------------------
! Section 7: Set dtssnow to T_S - T_SNOW
!------------------------------------------------------------------------------

IF (var_in(mztsnow_loc_in)%lreadin) THEN
  ! Use the already interpolated t_snow_lm which will be overwritten
  dtssnow_gl(:,:) = t_s_gl(:,:) - var_lm(mztsnow_loc_lm)%p2(:,:)
ELSE
  ! If tsnow is not available, set it to 0
  dtssnow_gl(:,:) = 0.0_ireals
ENDIF

END SUBROUTINE interpol_coarse_special_lm

!==============================================================================

SUBROUTINE interpol_coarse_special_ec (idbg)

!------------------------------------------------------------------------------
!
! Description:
!   Interpolation (horizontal and vertical) of ground fields and
!   computation of nonconstant fields not available in ECMWF input
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(IN) ::  &
  idbg                    ! debug level output

! Local variables

INTEGER (KIND=iintegers), PARAMETER :: nl_soil_in_max=5
INTEGER (KIND=iintegers) :: i,j,k,l,    & ! Loop indices
 izerror,                      & ! status and error status variable
 s_t, s_t_min, s_t_max,        & ! soiltypes of input
 i1,i2,                        & ! Interpolation indices
 groundl_in(nl_soil_in_max),   & ! Location of multi-level ground
 groundl_lm(nl_soil_in_max),   & ! variables in variable tables
 mzdpsdt_loc_in, mztg1_loc_in, & ! Locations of
 mztg2_loc_in, mztg3_loc_in,   & ! input variables
 mztcl_loc_in, mztsnow_loc_in, & ! in input
 mzqvs_loc_in, mzwg1_loc_in,   & ! variable table
 mzwg2_loc_in, mzwg3_loc_in,   & !
 mzwcl_loc_in, mzqc_loc_in,    & !
 mzqi_loc_in, mzts_loc_in,     & ! ----------------
 nlsoilin,                     &
 mzdpsdt_loc_lm, mzts_loc_lm,  & ! Locations of
 mztm_loc_lm, mztcl_loc_lm,    & ! output variables
 mztsnow_loc_lm,               & ! in output
 mzwg1_loc_lm, mzwg2_loc_lm,   & ! variable table
 mzwg3_loc_lm, mzwcl_loc_lm      ! ----------------

REAL (KIND=ireals) :: w1,w2,wei,           & ! Interpolaton weights
 zavg(nl_soil_in_max),                     & ! Average levels depth
 zx_g_gl(ie2lm,je2lm,0:nl_soil_in_max)       ! Automatic temporary array

REAL (KIND=ireals),ALLOCATABLE ::          &
 zczmls_in(:)                                ! Temporary for vertical levels

LOGICAL                    ::  &
  lzmono, lzposdef

CHARACTER (LEN=  1)        ::  &
  yzitype     ! interpolation type

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Look for locations of special variables in the variable tables
!------------------------------------------------------------------------------

IF (idbg > 10) THEN
  PRINT *, ' Starting interpolation of special fields'
ENDIF

yzroutine  = 'interpol_coarse_special_ec'

! Initialization of grdpt_rel_in
grdpt_rel_in = 0

DO i = 1, nvar_in
  SELECT CASE (var_in(i)%name)
  CASE ('DPSDT     '); mzdpsdt_loc_in   = i
  CASE ('T_SKIN    '); mzts_loc_in      = i
  CASE ('T_G1      '); mztg1_loc_in     = i
  CASE ('T_G2      '); mztg2_loc_in     = i
  CASE ('T_G3      '); mztg3_loc_in     = i
  CASE ('T_CL      '); mztcl_loc_in     = i
  CASE ('T_SNOW    '); mztsnow_loc_in   = i
  CASE ('QV_S      '); mzqvs_loc_in     = i
  CASE ('W_G1      '); mzwg1_loc_in     = i
  CASE ('W_G2      '); mzwg2_loc_in     = i
  CASE ('W_G3      '); mzwg3_loc_in     = i
  CASE ('W_CL      '); mzwcl_loc_in     = i
  CASE ('QC        '); mzqc_loc_in      = i
  CASE ('QI        '); mzqi_loc_in      = i
  END SELECT
ENDDO

DO i = 1, nvar_lm
  SELECT CASE (var_lm(i)%name)
  CASE ('DPSDT     '); mzdpsdt_loc_lm   = i
  CASE ('T_S       '); mzts_loc_lm      = i
  CASE ('T_M       '); mztm_loc_lm      = i
  CASE ('T_CL      '); mztcl_loc_lm     = i
  CASE ('T_SNOW    '); mztsnow_loc_lm   = i
  CASE ('W_G1      '); mzwg1_loc_lm     = i
  CASE ('W_G2      '); mzwg2_loc_lm     = i
  CASE ('W_G3      '); mzwg3_loc_lm     = i
  CASE ('W_CL      '); mzwcl_loc_lm     = i
  END SELECT
ENDDO

!------------------------------------------------------------------------------
! Input snow correction
!------------------------------------------------------------------------------

IF (var_in(mztsnow_loc_in)%lreadin) THEN

  IF (idbg > 10) THEN
    PRINT *, '   interpolation of t_snow_in'
  ENDIF

  WHERE(w_snow_in(:,:) <= 0.0_ireals)
    t_snow_in(:,:) = t_g1_in(:,:)
  END WHERE

  IF (var_in(mztsnow_loc_in)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mztsnow_loc_in)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF
  yzitype    = var_in(mzqvs_loc_in)%ipc(1:1)

  CALL interp_l(t_snow_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
              lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
              lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
              t_snow_lm, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,   &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
              grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
              yzerrmsg, izerror)
ENDIF

!------------------------------------------------------------------------------
! Section 2: Soil moisture interpolation
!------------------------------------------------------------------------------

IF (idbg > 10) THEN
  PRINT *, '   interpolation of soil moisture'
ENDIF

! Use of indirection in order to treat multi-level ground fields as a
! 1d array of 2d arrays
groundl_in(:)=0
groundl_in(1)=mzwg1_loc_in
groundl_in(2)=mzwg2_loc_in
groundl_in(3)=mzwg3_loc_in
groundl_in(4)=mzwcl_loc_in

IF (.NOT. lmulti_layer_in) THEN
! Count available soil levels
  DO k=1,nl_soil_in_max
    IF (groundl_in(k) == 0) EXIT
    IF (.NOT. var_in(groundl_in(k))%lreadin) EXIT
  ENDDO
  nlsoilin=k-1
ELSE
  nlsoilin = ke_soil_in + 1
ENDIF

s_t_min = MINVAL(NINT(soiltyp_in))
s_t_max = MAXVAL(NINT(soiltyp_in))
IF (idbg > 5) THEN
  WRITE (*,'(A)')      yzroutine//': NEW soil MOISTURE interpolation with SMI'
  WRITE (*,'(2(A,I3))')yzroutine//': WITH nlsoilin =',nlsoilin,' nl_soil_lm=',&
                                                                     nl_soil_lm
  WRITE (*,'(A,L3,A)') yzroutine//': and IF: l_smi = ',l_smi,                 &
                                           ' is .TRUE. THEN w_so_lm >= pwpb_ec'
  WRITE (*,'(A,I4,A,I4)') yzroutine//':      MIN / MAXVAL NINT(soiltyp_in):', &
                                                        s_t_min, ' / ', s_t_max
ENDIF

! LT: IFS has 7 soil types, not 6
IF (s_t_min < 0 .OR. s_t_max > 7 ) THEN
  WRITE (*,'(A,I4,A,I4)') yzroutine//': MIN / MAXVAL NINT(soiltyp_in):', &
                          s_t_min, ' / ', s_t_max
  yzerrmsg = 'ERROR: TOO MANY soil types in soiltyp_in'
  izerror = 2011
  CALL model_abort (my_cart_id, 1, yzerrmsg, yzroutine)
ELSEIF (s_t_min /= s_t_max ) THEN
  l = 0       ! there are MANY soil types
ELSE
  l = 1       ! there is ONLY ONE soil type
ENDIF

DO k = 1, nlsoilin
  DO i2 = 1, je_in
  DO i1 = 1, ie_in
    IF (l /= 0) THEN
      s_t = 1 ! there is ONLY ONE soil type
      w1 = pwpb_ec_1s
      w2 =  fcb_ec_1s - pwpb_ec_1s
    ELSE      ! there are MANY soil types
      s_t = NINT(soiltyp_in(i1,i2))
      IF (s_t == 7) s_t = 2
      IF (s_t /= 0) THEN
        w1 = pwpb_ec(s_t)
        w2 =  fcb_ec(s_t) - pwpb_ec(s_t) 
      ENDIF
    ENDIF
    wei = porb_ec_1s ! Although MANY soil types use ORIGINAL pore
                     ! volume of IFS with 1 soil type for l_smi == .F.
    IF (var_in(groundl_in(k))%p2(i1,i2) /= undef .AND. s_t /= 0) THEN
      IF (l_smi) THEN
        ! Scale the soil moisture (sm) WHICH is ALREADY multiplied by the layer
        ! depth (m) in src_gribtabs.f90 with the soil moisture index SMI:
        ! smi = (sm-pwp)/(fc-pwp)
        ! Caution: sm is in volumetric units, soil depth in cm!
        var_in(groundl_in(k))%p2(i1,i2) = var_in(groundl_in(k))%p2(i1,i2) /   &
             ((var_in(groundl_in(k))%levbot-var_in(groundl_in(k))%levtop)*0.01)
        var_in(groundl_in(k))%p2(i1,i2) =  MAX(                    &
              (var_in(groundl_in(k))%p2(i1,i2) - w1)/w2, 0.0_ireals)
      ELSE
        ! Scale the soil moisture by the pore volume of soil and
        ! by the layer depth (in m)
        var_in(groundl_in(k))%p2(i1,i2) =   var_in(groundl_in(k))%p2(i1,i2) / &
         (wei*(var_in(groundl_in(k))%levbot-var_in(groundl_in(k))%levtop)/100.)
      ENDIF ! l_smi
    ELSE
      var_in(groundl_in(k))%p2(i1,i2)=undef
    ENDIF 
  ENDDO ! i1
  ENDDO ! i2
! gdm end 

  ! Interpolate horizontally and store in temporary array
  IF (var_in(groundl_in(k))%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(groundl_in(k))%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF
  yzitype    = var_in(groundl_in(k))%ipc(1:1)
  IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
  CALL interp_l(var_in(groundl_in(k))%p2,                                  &
                ie_in, je_in, i_index(:,:,1), j_index(:,:,1),              &
                lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
                lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
                zx_g_gl(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                yzerrmsg, izerror)
ENDDO ! k

! Avoid undefined
WHERE (zx_g_gl(:,:,1:nlsoilin) == undef)
  zx_g_gl(:,:,1:nlsoilin) = 0.0_ireals
ENDWHERE

IF (.NOT. lmulti_layer_in) THEN

! Use of indirection in order to treat multi-level ground fields as a
! 1d array of 2d arrays
  groundl_lm(1)=mzwg1_loc_lm
  groundl_lm(2)=mzwg2_loc_lm
  IF (nl_soil_lm == 3) THEN
    groundl_lm(3)=mzwg3_loc_lm
    groundl_lm(4)=mzwcl_loc_lm
  ELSE
    groundl_lm(3)=mzwcl_loc_lm
  ENDIF

! Conservative vertical interpolation of soil wetness, fields are
! already scaled with the depth of the layer
  DO l=1,nl_soil_lm+1
    var_lm(groundl_lm(l))%p2(:,:)=0.0_ireals
    DO k=1,nlsoilin
      wei=MIN(var_in(groundl_in(k))%levbot, var_lm(groundl_lm(l))%levbot) - &
       MAX(var_in(groundl_in(k))%levtop, var_lm(groundl_lm(l))%levtop)
      IF (wei > 0.0_ireals) THEN  ! The intervals intersect
        wei=wei / (var_lm(groundl_lm(l))%levbot-var_lm(groundl_lm(l))%levtop)
        var_lm(groundl_lm(l))%p2(:,:) = var_lm(groundl_lm(l))%p2(:,:) &
         + wei*zx_g_gl(:,:,k)
      ENDIF
    ENDDO

! Check whether input data is not deep enough
    wei=var_lm(groundl_lm(l))%levbot - &
     MAX(var_in(groundl_in(nlsoilin))%levbot, var_lm(groundl_lm(l))%levtop)
    IF (wei > 0.0_ireals) THEN
! If so extend the lowest input field up to the output bottom
      wei=wei / (var_lm(groundl_lm(l))%levbot - var_lm(groundl_lm(l))%levtop)
      var_lm(groundl_lm(l))%p2(:,:) = var_lm(groundl_lm(l))%p2(:,:) + &
       wei*zx_g_gl(:,:,nlsoilin)
    ENDIF
  ENDDO

! integrate over lm layer depth (in m)
  DO l=1,nl_soil_lm+1
    var_lm(groundl_lm(l))%p2=var_lm(groundl_lm(l))%p2 * &
     ((var_lm(groundl_lm(l))%levbot-var_lm(groundl_lm(l))%levtop)/100.)
  ENDDO

ELSE ! lmulti_layer_in (implies lmulti_layer_lm)

  w_so_lm(:,:,:) = 0.0_ireals

! Conservative vertical interpolation of soil wetness, fields are
! already scaled with the depth of the layer
  DO l = 1, ke_soil_lm + 1

    DO k = 1, ke_soil_in + 1
      wei = MIN(czhls_in(k), czhls_lm(l)) - &
       MAX(czhls_in(k-1), czhls_lm(l-1))
!      wei=MIN(var_in(groundl_in(k))%levbot, var_lm(groundl_lm(l))%levbot) - &
!       MAX(var_in(groundl_in(k))%levtop, var_lm(groundl_lm(l))%levtop)
      IF (wei > 0.0_ireals) THEN ! the intervals intersect
        wei = wei/(czhls_lm(l) - czhls_lm(l-1))
        PRINT*,l,k,wei
!        wei=wei / (var_lm(groundl_lm(l))%levbot-var_lm(groundl_lm(l))%levtop)
        w_so_lm(:,:,l) = w_so_lm(:,:,l) + wei*zx_g_gl(:,:,k)
      ENDIF
    ENDDO

! Check whether input data is not deep enough
!    wei = var_lm(groundl_lm(l))%levbot - &
!     MAX(var_in(groundl_in(nlsoilin))%levbot, var_lm(groundl_lm(l))%levtop)
    wei = czhls_lm(l) - MAX(czhls_in(ke_soil_in+1), czhls_lm(l-1))
    IF (wei > 0.0_ireals) THEN
! If so extend the lowest input field up to the output bottom
!      wei = wei/(var_lm(groundl_lm(l))%levbot - var_lm(groundl_lm(l))%levtop)
      wei = wei/(czhls_lm(l) - czhls_lm(l-1))
      PRINT*,'bottom',l,k,wei
      w_so_lm(:,:,l) = w_so_lm(:,:,l) + wei*zx_g_gl(:,:,ke_soil_in+1)
    ENDIF
  ENDDO

ENDIF ! lmulti_layer_in

!------------------------------------------------------------------------------
! Section 3: Soil temperature interpolation
!------------------------------------------------------------------------------

IF (.NOT. lmulti_layer_in) THEN

! Use of indirection in order to treat multi-level ground fields as a
! 1d array of 2d arrays
  groundl_in(:)=0
  IF (luse_t_skin) THEN
    groundl_in(1)=mzts_loc_in
    groundl_in(2)=mztg1_loc_in
    groundl_in(3)=mztg2_loc_in
    groundl_in(4)=mztg3_loc_in
    groundl_in(5)=mztcl_loc_in
  ELSE
    groundl_in(1)=mztg1_loc_in
    groundl_in(2)=mztg2_loc_in
    groundl_in(3)=mztg3_loc_in
    groundl_in(4)=mztcl_loc_in
  ENDIF

  groundl_lm(1)=mzts_loc_lm
  groundl_lm(2)=mztm_loc_lm
  groundl_lm(3)=mztcl_loc_lm

! Count available input soil levels
  DO k=1,nl_soil_in_max
    IF (groundl_in(k) == 0) EXIT
    IF (.NOT. var_in(groundl_in(k))%lreadin) EXIT
  ENDDO
  nlsoilin=k-1
  DO k=1,nlsoilin
! Interpolate horizontally and store temporarily in lm arrays
    IF (var_in(groundl_in(k))%ipc(2:2) == 'T') THEN
      lzmono   = .TRUE.
    ELSE
      lzmono   = .FALSE.
    ENDIF
    IF (var_in(groundl_in(k))%ipc(3:3) == 'T') THEN
      lzposdef = .TRUE.
    ELSE
      lzposdef = .FALSE.
    ENDIF
    yzitype    = var_in(groundl_in(k))%ipc(1:1)
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
    CALL interp_l(var_in(groundl_in(k))%p2,                                  &
     ie_in, je_in, i_index(:,:,1), j_index(:,:,1),              &
     lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
     lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
     zx_g_gl(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
     latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
     grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
     yzerrmsg, izerror)
  ENDDO

! Avoid undefined
  WHERE (zx_g_gl(:,:,1:nlsoilin) == undef)
    zx_g_gl(:,:,1:nlsoilin) = 288.0_ireals
  ENDWHERE

! Interpolate vertically
  IF (nlsoilin == 1) THEN
! Only 1 input level available, no way to interpolate
    DO l=1,3
      var_lm(groundl_lm(l))%p2=zx_g_gl(:,:,1)
    ENDDO
  ELSE
! Compute average depths
    DO k=1,nlsoilin
      zavg(k)=(var_in(groundl_in(k))%levtop+var_in(groundl_in(k))%levbot)*0.5
    ENDDO
    DO l=1,3
      DO k=1,nlsoilin
        IF (zavg(k) > var_lm(groundl_lm(l))%levbot .OR. k == nlsoilin) THEN
! Compute indices and weights
          i1=MAX(k-1,1)
          i2=i1+1
          w1=(zavg(i2)-var_lm(groundl_lm(l))%levbot)/(zavg(i2)-zavg(i1))
          w2=(var_lm(groundl_lm(l))%levbot-zavg(i1))/(zavg(i2)-zavg(i1))
          EXIT
        ENDIF
      ENDDO
! Interpolate, remember to check undefs!!
      var_lm(groundl_lm(l))%p2(:,:)=zx_g_gl(:,:,i1)*w1 + zx_g_gl(:,:,i2)*w2
    ENDDO

  ENDIF

! Reassign the variables for further processing
  dtms_gl(:,:) = var_lm(mztm_loc_lm)%p2(:,:) - var_lm(mzts_loc_lm)%p2(:,:)
  t_s_gl(:,:)  = var_lm(mzts_loc_lm)%p2(:,:)

ELSE ! lmulti_layer_in (implies lmulti_layer_lm)

  IF (.NOT. lcomp_bound) THEN

! Interpolate horizontally and store temporarily in dt_so_gl
    lzmono = (var_in(mzts_loc_in)%ipc(2:2) == 'T')
    lzposdef = (var_in(mzts_loc_in)%ipc(3:3) == 'T')
    yzitype = var_in(mzts_loc_in)%ipc(1:1)
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1

    k = 0
!    IF (luse_t_skin) THEN
      CALL interp_l(var_in(mzts_loc_in)%p2, &
       ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
       lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in, &
       lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
       zx_g_gl(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
       latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
       grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman, &
       yzerrmsg, izerror)
      k = k + 1
!    ENDIF

    CALL interp_l(var_in(mztg1_loc_in)%p2, &
     ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
     lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in, &
     lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
     zx_g_gl(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
     latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
     grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman, &
     yzerrmsg, izerror)
    k = k + 1

    CALL interp_l(var_in(mztg2_loc_in)%p2, &
     ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
     lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in, &
     lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
     zx_g_gl(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
     latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
     grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman, &
     yzerrmsg, izerror)
    k = k + 1

    CALL interp_l(var_in(mztg3_loc_in)%p2, &
     ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
     lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in, &
     lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
     zx_g_gl(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
     latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
     grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman, &
     yzerrmsg, izerror)
    k = k + 1

    CALL interp_l(var_in(mztcl_loc_in)%p2, &
     ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
     lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in, &
     lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
     zx_g_gl(:,:,k), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
     latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
     grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman, &
     yzerrmsg, izerror)

    ALLOCATE(zczmls_in(0:ke_soil_in+1))
!    IF (luse_t_skin) THEN
      zczmls_in(0) = 0.0_ireals
      zczmls_in(1:ke_soil_in+1) = czmls_in(1:ke_soil_in+1)
!      zn_soil_in = ke_soil_in + 2
!    ELSE
!      zczmls_in(1:ke_soil_in+1) = czmls_in(1:ke_soil_in+1)
!      zn_soil_in = ke_soil_in
!    ENDIF

    IF (ke_soil_in >= 1) THEN
      DO l = 0, ke_soil_lm + 1
        k = 1 ! skip 0
        DO WHILE (zczmls_in(k) <= czhls_lm(l) .AND. k < ke_soil_in + 1)
          k = k + 1
        ENDDO
! Compute indices and weights
        i1 = k - 1
        i2 = k
! avoid extrapolation below the bottom
        w1 = MAX((zczmls_in(i2) - czhls_lm(l))/(zczmls_in(i2) - zczmls_in(i1)), 0.0_ireals)
        w2 = 1.0_ireals - w1
! Interpolate
        WHERE (zx_g_gl(:,:,i1) /= undef .AND. zx_g_gl(:,:,i2) /= undef)
          dt_so_gl(:,:,l) = zx_g_gl(:,:,i1)*w1 + zx_g_gl(:,:,i2)*w2
        ELSEWHERE
          dt_so_gl(:,:,l) = undef
        END WHERE
      ENDDO
    ELSE ! impossible case at the moment
      DO l = 0, ke_soil_lm + 1
        dt_so_gl(:,:,l) = zx_g_gl(:,:,1)
      ENDDO
    ENDIF

    DO l = ke_soil_lm + 1, 1, -1
      WHERE(dt_so_gl(:,:,l-1) /= undef .AND. dt_so_gl(:,:,l) /= undef)
        dt_so_gl(:,:,l) = dt_so_gl(:,:,l) - dt_so_gl(:,:,l-1)
      ELSEWHERE
        dt_so_gl(:,:,l) = undef
      END WHERE
    ENDDO

    DEALLOCATE(zczmls_in)

  ENDIF ! lcomp_bound

  ! Interpolate T_SO(mlev=0) to t_s_gl
  ! (in any case: for initial and boundary values)

  lzmono = (var_in(mzts_loc_in)%ipc(2:2) == 'T')
  lzposdef = (var_in(mzts_loc_in)%ipc(3:3) == 'T')
  yzitype = var_in(mzts_loc_in)%ipc(1:1)
  IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
  
  CALL interp_l(var_in(mzts_loc_in)%p2, &
   ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
   lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in, &
   lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1), &
   t_s_gl(:,:), 1, ie2lm, 1, je2lm, startlat_in, startlon_in, &
   latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
   grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman, &
   yzerrmsg, izerror)

ENDIF ! lmulti_layer_in

!------------------------------------------------------------------------------
! Section 4: Set DPSDT to 0.
!------------------------------------------------------------------------------

!IF (.NOT. var_in(mzdpsdt_loc_in)%lreadin) THEN
!  ! If dpsdt is not available, and this is the case, set it to 0.
!  var_lm(mzdpsdt_loc_lm)%p2 = 0.0_ireals
!ENDIF

!------------------------------------------------------------------------------
! Section 5: Set dtkes_gl to T(lowlev) - T_S
!------------------------------------------------------------------------------

CALL interp_l(t_in(:,:,ke_in), ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
              .FALSE., .FALSE., 'L', lbd_frame_cur, lolp_in,                 &
              lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),        &
              dtkes_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,        &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,       &
              grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                &
              yzerrmsg, izerror)
dtkes_gl(:,:) = dtkes_gl(:,:) - t_s_gl(:,:)

!------------------------------------------------------------------------------
! Section 6: Interpolate qv_s and convert in rh_s
!------------------------------------------------------------------------------

IF (.NOT. var_in(mzqvs_loc_in)%lreadin) THEN
  qv_s_in = qv_in(:,:,ke_in)
ENDIF

IF (var_in(mzqvs_loc_in)%ipc(2:2) == 'T') THEN
  lzmono   = .TRUE.
ELSE
  lzmono   = .FALSE.
ENDIF
IF (var_in(mzqvs_loc_in)%ipc(3:3) == 'T') THEN
  lzposdef = .TRUE.
ELSE
  lzposdef = .FALSE.
ENDIF
yzitype    = var_in(mzqvs_loc_in)%ipc(1:1)
  IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1

CALL interp_l(qv_s_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),     &
              lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
              lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
              rh_s_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,     &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
              grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
              yzerrmsg, izerror)

DO j = 1, je2lm
  DO i = 1, ie2lm
    rh_s_gl(i,j) = rh_s_gl(i,j) /  &
     qsat(psat_w(t_s_gl(i,j),b1, b2_w, b3, b4_w),   &
     ps_gl(i,j), rdv, O_m_rdv)
  ENDDO
ENDDO

!------------------------------------------------------------------------------
! Section 7: Set dtssnow to T_S - T_SNOW
!------------------------------------------------------------------------------

IF (var_in(mztsnow_loc_in)%lreadin) THEN
  ! Use the already interpolated t_snow_lm which will be overwritten
  dtssnow_gl(:,:) = var_lm(mzts_loc_lm)%p2(:,:) - var_lm(mztsnow_loc_lm)%p2(:,:)
ELSE
  ! If tsnow is not available, set dtssnow to 0
  dtssnow_gl(:,:) = 0.0_ireals
ENDIF

!------------------------------------------------------------------------------
! Section 8: Add QI to QC
!------------------------------------------------------------------------------

! This comment surely belongs to Section 8 of this subroutine
! (it was at the now deleted Section 8 of SR interpol_coarse_special_lm)
! roa begin
! WARNING!!!
! this part of the code is never used: it is not possible that
! var_in(mzqi_loc_in)%lreadin exists if lprog_qi=.FALSE.
! Indeed if lprog_qi=.FALSE., qi doesn't exist through the
! whole programm INT2LM
! Thus Qi is never added to Qc
! roa begin

IF (var_in(mzqi_loc_in)%lreadin .AND. .NOT.lprog_qi) THEN
  IF (var_in(mzqi_loc_in)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mzqi_loc_in)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF
  yzitype    = var_in(mzqi_loc_in)%ipc(1:1)
  IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
  DO l=1,ke_in
    CALL interp_l(qi_in(:,:,l), ie_in, je_in, i_index(:,:,1), j_index(:,:,1),&
                  lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,         &
                  lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),    &
                  zx_g_gl(:,:,1), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,&
                  latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,   &
                  grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,            &
                  yzerrmsg, izerror)
    qc_lm(:,:,l) = qc_lm(:,:,l) + zx_g_gl(:,:,1)
  ENDDO
ENDIF

IF (izerror /= 0) THEN
  CALL model_abort (my_cart_id, 1, yzerrmsg, yzroutine)
ENDIF

! Avoid undefined
WHERE (w_i_lm(:,:) == undef)
  w_i_lm(:,:) = 0.0_ireals
ENDWHERE

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE interpol_coarse_special_ec

!==============================================================================
!==============================================================================

SUBROUTINE interpol_coarse_special_gfs(idbg)

!------------------------------------------------------------------------------
!
! Description:
!   Interpolation of surface temperature to t_s_lm
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(IN) ::  &
  idbg                    ! debug level output

! Local variables

INTEGER (KIND=iintegers) :: i,j,k,l,    & ! Loop indices
 izerror,                      & ! status and error status variable
 mzdpsdt_loc_in, mzts_loc_in,  & ! Locations of
 mztm_loc_in,                  & !   input variables
 mztso_loc_in, mzwso_loc_in,   & !   input variables
 mztcl_loc_in, mztsnow_loc_in, & !   in input
 mzqvs_loc_in, mzwg1_loc_in,   & !   variable table
 mzwg2_loc_in, mzwg3_loc_in,   & !
 mzwcl_loc_in,                 & !
 mzqi_loc_in,                  & ! ----------------
 mzdpsdt_loc_lm,               & ! Locations of
 mztm_loc_lm, mztcl_loc_lm,    & !   output variables
 mztsnow_loc_lm,               & !   in output
 mzwg1_loc_lm, mzwg2_loc_lm,   & !   variable table
 mzwg3_loc_lm, mzwcl_loc_lm      ! ----------------

LOGICAL                    ::  &
  lzmono, lzposdef

CHARACTER (LEN=  1)        ::  &
  yzitype     ! interpolation type

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Look for locations of special variables in the variable tables
!------------------------------------------------------------------------------

yzroutine = 'interpol_coarse_special_gfs'

! Initialization of grdpt_rel_in
grdpt_rel_in = 0

DO i = 1, nvar_in
  SELECT CASE (var_in(i)%name)
  CASE ('DPSDT     '); mzdpsdt_loc_in  = i
  CASE ('T_S       '); mzts_loc_in     = i
  CASE ('T_M       '); mztm_loc_in     = i
  CASE ('T_CL      '); mztcl_loc_in    = i
  CASE ('T_SNOW    '); mztsnow_loc_in  = i
  CASE ('QV_S      '); mzqvs_loc_in    = i
  CASE ('W_G1      '); mzwg1_loc_in    = i
  CASE ('W_G2      '); mzwg2_loc_in    = i
  CASE ('W_G3      '); mzwg3_loc_in    = i
  CASE ('W_CL      '); mzwcl_loc_in    = i
  CASE ('QI        '); mzqi_loc_in     = i
  CASE ('T_SO      '); mztso_loc_in    = i
  CASE ('W_SO      '); mzwso_loc_in    = i
  END SELECT
ENDDO

DO i = 1, nvar_lm
  SELECT CASE (var_lm(i)%name)
  CASE ('DPSDT     '); mzdpsdt_loc_lm  = i
  CASE ('T_M       '); mztm_loc_lm     = i
  CASE ('T_CL      '); mztcl_loc_lm    = i
  CASE ('T_SNOW    '); mztsnow_loc_lm  = i
  CASE ('W_G1      '); mzwg1_loc_lm    = i
  CASE ('W_G2      '); mzwg2_loc_lm    = i
  CASE ('W_G3      '); mzwg3_loc_lm    = i
  CASE ('W_CL      '); mzwcl_loc_lm    = i
  END SELECT
ENDDO

!------------------------------------------------------------------------------
! Section 2: Set dtkes_gl to T(lowlev) - T_S
!------------------------------------------------------------------------------

CALL interp_l(t_in(:,:,ke_in), ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
              .FALSE., .FALSE., 'L', lbd_frame_cur, lolp_in,                 &
              lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),        &
              dtkes_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,        &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,       &
              grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                &
              yzerrmsg, izerror)
dtkes_gl(:,:) = dtkes_gl(:,:) - t_s_gl(:,:)

!------------------------------------------------------------------------------
! Section 3: Interpolate qv_s and convert in rh_s
!------------------------------------------------------------------------------

IF (.NOT. var_in(mzqvs_loc_in)%lreadin) THEN
  qv_s_in = qv_in(:,:,ke_in)
ENDIF

IF (var_in(mzqvs_loc_in)%ipc(2:2) == 'T') THEN
  lzmono   = .TRUE.
ELSE
  lzmono   = .FALSE.
ENDIF
IF (var_in(mzqvs_loc_in)%ipc(3:3) == 'T') THEN
  lzposdef = .TRUE.
ELSE
  lzposdef = .FALSE.
ENDIF
yzitype    = var_in(mzqvs_loc_in)%ipc(1:1)
IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1

CALL interp_l(qv_s_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),   &
              lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,       &
              lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),  &
              rh_s_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,   &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m, &
              grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,          &
              yzerrmsg, izerror)

DO j = 1, je2lm
  DO i = 1, ie2lm
     rh_s_gl(i,j) = rh_s_gl(i,j) /  &
     qsat(psat_w(t_s_gl(i,j),B1, B2_w, B3, B4_w),   &
     ps_gl(i,j), rdv, O_m_rdv)
  ENDDO
ENDDO

END SUBROUTINE interpol_coarse_special_gfs

!==============================================================================

SUBROUTINE interpol_coarse_special_cm (idbg)

!------------------------------------------------------------------------------
!
! Description:
!   Interpolation (horizontal and vertical) of ground fields and
!   computation of nonconstant fields not available in CM input
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(IN) ::  &
  idbg                    ! debug level output

! Local variables

INTEGER (KIND=iintegers), PARAMETER :: nl_soil_in_max=5
INTEGER (KIND=iintegers) :: i,j,k,l,    & ! Loop indices
 izerror,                      & ! status and error status variable
 i1,i2,                        & ! Interpolation indices
 mzdpsdt_loc_in, mztskin_loc_in, & ! Locations of input variables
 mztcl_loc_in, mztsnow_loc_in, & ! in input
 mzqvs_loc_in,                 & ! variable table
 mzwcl_loc_in, mzqc_loc_in,    & !
 mzqi_loc_in, mzts_loc_in,     & ! ----------------
 nlsoilin,                     &
 mzdpsdt_loc_lm, mzts_loc_lm,  & ! Locations of
 mztm_loc_lm, mztcl_loc_lm,    & ! output variables
 mztsnow_loc_lm, mzqvs_loc_lm, & ! in output variable table
 mzwcl_loc_lm      ! ----------------

REAL (KIND=ireals) :: w1,w2,wei,           & ! Interpolaton weights
 zavg(nl_soil_in_max),                     & ! Average levels depth
 zx_g_gl(ie2lm,je2lm,nl_soil_in_max)         ! Automatic temporary array

LOGICAL                    ::  &
  lzmono, lzposdef

CHARACTER (LEN=  1)        ::  &
  yzitype     ! interpolation type

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Look for locations of special variables in the variable tables
!------------------------------------------------------------------------------

yzroutine  = 'interpol_coarse_special_cm'

! Initialization of grdpt_rel_in
grdpt_rel_in = 0

DO i = 1, nvar_in
  SELECT CASE (var_in(i)%name)
  CASE ('DPSDT     '); mzdpsdt_loc_in   = i
  CASE ('T_S       '); mzts_loc_in      = i
  CASE ('T_SKIN    '); mztskin_loc_in   = i
! SP, 201405
  CASE ('T_CL      '); mztcl_loc_in   = i
  CASE ('T_SNOW    '); mztsnow_loc_in  = i
  CASE ('QV_S      '); mzqvs_loc_in    = i
  CASE ('W_CL      '); mzwcl_loc_in    = i
  CASE ('QC        '); mzqc_loc_in     = i
  CASE ('QI        '); mzqi_loc_in     = i
  END SELECT
ENDDO

DO i = 1, nvar_lm
  SELECT CASE (var_lm(i)%name)
  CASE ('DPSDT     '); mzdpsdt_loc_lm  = i
  CASE ('T_S       '); mzts_loc_lm     = i
  CASE ('T_M       '); mztm_loc_lm     = i
! SP, 201405
  CASE ('T_CL      '); mztcl_loc_lm    = i
  CASE ('T_SNOW    '); mztsnow_loc_lm  = i
  CASE ('W_CL      '); mzwcl_loc_lm    = i
  END SELECT
ENDDO

!------------------------------------------------------------------------------
! Section 2: Soil moisture interpolation
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Section 3: Soil temperature interpolation
!------------------------------------------------------------------------------

  lzmono   = .FALSE.
  lzposdef = .FALSE.
  yzitype    = 'M'

  IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1

  CALL interp_l(t_s_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),      &
        lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in, lolp_lm,        &
        lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),                     &
        t_s_gl(:,:), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,         &
        latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,           &
        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                    &
        yzerrmsg, izerror)

  IF (var_in(mztskin_loc_in)%lreadin) THEN
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
    CALL interp_l(t_skin_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
        lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in, lolp_lm,        &
        lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),                     &
        t_skin_gl(:,:), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,      &
        latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,           &
        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                    &
        yzerrmsg, izerror)
  ENDIF

! SP, 201405
  IF (var_in(mztcl_loc_in)%lreadin) THEN
    IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
    CALL interp_l(t_cl_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),   &
        lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in, lolp_lm,        &
        lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),                     &
        t_cl_lm(:,:), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,        &
        latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,           &
        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                    &
        yzerrmsg, izerror)
  ENDIF


!------------------------------------------------------------------------------
! Section 4: Set DPSDT to 0.
!------------------------------------------------------------------------------

!IF (.NOT. var_in(mzdpsdt_loc_in)%lreadin) THEN
!  ! If dpsdt is not available, and this is the case, set it to 0.
!  var_lm(mzdpsdt_loc_lm)%p2 = 0.0_ireals
!ENDIF

!------------------------------------------------------------------------------
! Section 5: Set dtkes_gl to T(lowlev) - T_S
!------------------------------------------------------------------------------

  IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
  CALL interp_l(t_in(:,:,ke_in), ie_in, je_in, i_index(:,:,1), j_index(:,:,1),&
        .FALSE., .FALSE., 'L', lbd_frame_cur, lolp_in, lolp_lm,            &
        lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),                     &
        dtkes_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,            &
        latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,           &
        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                    &
        yzerrmsg, izerror)

  dtkes_gl(:,:) = dtkes_gl(:,:) - t_s_gl(:,:)

! iso code Hui Tang 2013-11-20
! calculate "drisoke_gl" for height correction of R18OSOIL and R2HSOIL
  IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
  CALL interp_l(riso_in(:,:,ke_in,1), ie_in, je_in, i_index(:,:,1), j_index(:,:,1),&
        .FALSE., .TRUE., 'L', lbd_frame_cur, lolp_in, lolp_lm,                     &
        lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),                             &
        drisoke_gl(:,:,1), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,           &
        latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,                   &
        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                            &
        yzerrmsg, izerror)
  drisoke_gl(:,:,1) = drisoke_gl(:,:,1) - risosoil_lm(:,:,1)
  
  IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
  CALL interp_l(riso_in(:,:,ke_in,2), ie_in, je_in, i_index(:,:,1), j_index(:,:,1),&
        .FALSE., .TRUE., 'L', lbd_frame_cur, lolp_in, lolp_lm,                     &
        lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),                             &
        drisoke_gl(:,:,2), 1, ie2lm, 1, je2lm, startlat_in, startlon_in,           &
        latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,                   &
        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                            &
        yzerrmsg, izerror)
  drisoke_gl(:,:,2) = drisoke_gl(:,:,2) - risosoil_lm(:,:,2) 
! end iso code

!------------------------------------------------------------------------------
! Section 6: Interpolate qv_s and convert in rh_s
!------------------------------------------------------------------------------

IF (.NOT. var_in(mzqvs_loc_in)%lreadin) THEN
  qv_s_in = qv_in(:,:,ke_in)
ENDIF

IF (var_in(mzqvs_loc_in)%ipc(2:2) == 'T') THEN
  lzmono   = .TRUE.
ELSE
  lzmono   = .FALSE.
ENDIF
IF (var_in(mzqvs_loc_in)%ipc(3:3) == 'T') THEN
  lzposdef = .TRUE.
ELSE
  lzposdef = .FALSE.
ENDIF
yzitype    = var_in(mzqvs_loc_in)%ipc(1:1)

IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
CALL interp_l(qv_s_in, ie_in, je_in, i_index(:,:,1), j_index(:,:,1),       &
        lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,                 &
        lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),            &
        rh_s_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,             &
        latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,           &
        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                    &
        yzerrmsg, izerror)
DO j = 1, je2lm
  DO i = 1, ie2lm
    rh_s_gl(i,j) = rh_s_gl(i,j) /  &
    qsat(psat_w(t_s_gl(i,j),b1, b2_w, b3, b4_w),   &
    ps_gl(i,j), rdv, O_m_rdv)
  ENDDO
ENDDO

!------------------------------------------------------------------------------
! Section 7: Set dtssnow to T_S - T_SNOW
!------------------------------------------------------------------------------

IF (var_in(mztsnow_loc_in)%lreadin) THEN
  ! Use the already interpolated t_snow_lm which will be overwritten
  dtssnow_gl(:,:) = t_s_gl(:,:) - var_lm(mztsnow_loc_lm)%p2(:,:)
ELSE IF (var_in(mztskin_loc_in)%lreadin) THEN
  ! Use tskin, if available
  WHERE (w_snow_lm(:,:)/=0.0_ireals) &
    dtssnow_gl(:,:) = t_s_gl(:,:) - t_skin_gl(:,:)
ELSE
  ! If tsnow is not available, set it to 0
  dtssnow_gl(:,:) = 0.0_ireals
ENDIF

!------------------------------------------------------------------------------
! Section 8: Add QI to QC
!------------------------------------------------------------------------------

IF (var_in(mzqi_loc_in)%lreadin .AND. .NOT.lprog_qi) THEN
  IF (var_in(mzqi_loc_in)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mzqi_loc_in)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF
  yzitype    = var_in(mzqi_loc_in)%ipc(1:1)
! PIK U. Boehm - 07.09.05
  IF (yzitype=='M' .AND. l_cressman) grdpt_rel_in=1
! PIK U. Boehm - End
  DO l=1,ke_in
    CALL interp_l(qi_in(:,:,l), ie_in, je_in, i_index(:,:,1), j_index(:,:,1),&
        lzmono, lzposdef, yzitype, lbd_frame_cur, lolp_in,                 &
        lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),            &
        zx_g_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,             &
        latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,           &
        grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                    &
        yzerrmsg, izerror)
    qc_lm(:,:,l) = qc_lm(:,:,l) + zx_g_gl(:,:,1)
  ENDDO
ENDIF

IF (izerror /= 0) THEN
  CALL model_abort (my_cart_id, 1, yzerrmsg, yzroutine)
ENDIF

! Avoid undefined
WHERE (w_i_lm(:,:) == undef)
  w_i_lm(:,:) = 0.0_ireals
ENDWHERE

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE interpol_coarse_special_cm

!==============================================================================
!==============================================================================

SUBROUTINE interpol_coarse_special_hir(idbg)

!------------------------------------------------------------------------------
!
! Description:
!   Interpolation of surface temperature to t_s_lm
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(IN) ::  &
  idbg                    ! debug level output

! Local variables

INTEGER (KIND=iintegers) :: i,j,k,l,    & ! Loop indices
 izerror,                      & ! status and error status variable
 mzdpsdt_loc_in, mzts_loc_in,  & ! Locations of
 mztm_loc_in,                  & !   input variables
 mztso_loc_in, mzwso_loc_in,   & !   input variables
 mztcl_loc_in, mztsnow_loc_in, & !   in input
 mzqvs_loc_in, mzwg1_loc_in,   & !   variable table
 mzwg2_loc_in, mzwg3_loc_in,   & !
 mzwcl_loc_in,                 & !
 mzqi_loc_in,                  & ! ----------------
 mzdpsdt_loc_lm,               & ! Locations of
 mztm_loc_lm, mztcl_loc_lm,    & !   output variables
 mztsnow_loc_lm,               & !   in output
 mzwg1_loc_lm, mzwg2_loc_lm,   & !   variable table
 mzwg3_loc_lm, mzwcl_loc_lm      ! ----------------

LOGICAL                    ::  &
  lzmono, lzposdef

CHARACTER (LEN=  1)        ::  &
  yzitype     ! interpolation type

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Look for locations of special variables in the variable tables
!------------------------------------------------------------------------------

yzroutine = 'interpol_coarse_special_hir'

! Initialization of grdpt_rel_in
grdpt_rel_in = 0

DO i = 1, nvar_in
  SELECT CASE (var_in(i)%name)
  CASE ('DPSDT     '); mzdpsdt_loc_in  = i
  CASE ('T_S       '); mzts_loc_in     = i
  CASE ('T_M       '); mztm_loc_in     = i
  CASE ('T_CL      '); mztcl_loc_in    = i
  CASE ('T_SNOW    '); mztsnow_loc_in  = i
  CASE ('QV_S      '); mzqvs_loc_in    = i
  CASE ('W_G1      '); mzwg1_loc_in    = i
  CASE ('W_G2      '); mzwg2_loc_in    = i
  CASE ('W_G3      '); mzwg3_loc_in    = i
  CASE ('W_CL      '); mzwcl_loc_in    = i
  CASE ('QI        '); mzqi_loc_in     = i
  CASE ('T_SO      '); mztso_loc_in    = i
  CASE ('W_SO      '); mzwso_loc_in    = i
  END SELECT
ENDDO

DO i = 1, nvar_lm
  SELECT CASE (var_lm(i)%name)
  CASE ('DPSDT     '); mzdpsdt_loc_lm  = i
  CASE ('T_M       '); mztm_loc_lm     = i
  CASE ('T_CL      '); mztcl_loc_lm    = i
  CASE ('T_SNOW    '); mztsnow_loc_lm  = i
  CASE ('W_G1      '); mzwg1_loc_lm    = i
  CASE ('W_G2      '); mzwg2_loc_lm    = i
  CASE ('W_G3      '); mzwg3_loc_lm    = i
  CASE ('W_CL      '); mzwcl_loc_lm    = i
  END SELECT
ENDDO

!------------------------------------------------------------------------------
! Section 2: Set dtkes_gl to T(lowlev) - T_S
!------------------------------------------------------------------------------

CALL interp_l(t_in(:,:,ke_in), ie_in, je_in, i_index(:,:,1), j_index(:,:,1), &
              .FALSE., .FALSE., 'L', lbd_frame_cur, lolp_in,                 &
              lolp_lm, lmask_lm, undef, x_wght(:,:,1), y_wght(:,:,1),        &
              dtkes_gl, 1, ie2lm, 1, je2lm, startlat_in, startlon_in,        &
              latitudes_in, longitudes_in, lat_coarse_m, lon_coarse_m,       &
              grdpt_rel_in, ie_in_tot, je_in_tot, l_cressman,                &
              yzerrmsg, izerror)
dtkes_gl(:,:) = dtkes_gl(:,:) - t_s_gl(:,:)

END SUBROUTINE interpol_coarse_special_hir

!==============================================================================

END MODULE src_coarse_interpol
