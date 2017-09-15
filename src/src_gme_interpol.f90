!+ Source Module for reading GME-fields and interpolating to LM-fields
!==============================================================================

MODULE src_gme_interpol

!==============================================================================
!
! Description:
!   This module contains subroutines necessary for reading the GME-fields
!   and interpolating them to the LM-fields.
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
!  Interpolation of T_CL only for initial data
! 1.3        2005/12/12 Ulrich Schaettler
!  Added interpolation of rho_snow for prognostic treatment in the LM
! V1_5         2007/07/09 Ulrich Schaettler
!  Distribution of nunit_of_time only for num_compute larger 1
!  Eliminated 2 fields in call to interpol_gme_special
!  Corrections for bitmap checking
!  Replaced ke_soil to ke_soil_lm
!  Adaptations to changes in io_utilities
! V1_6         2007/09/07 Ulrich Schaettler
!  Introduced yinput_type and eliminated lanalysis
! V1_7         2007/11/26 Ulrich Schaettler
!  Added several debug printouts
!  Introduced itype_t_cl to choose t_cl from correct source
! V1_8         2008/05/29 Ulrich Schaettler
!  Interpolation of qr,qs also for GME variables
!  Enlarged some character strings for I/O paths
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
! V1_9         2009/09/03 Guy DeMorsier
!  Implemented option to choose new soil properties with l_smi
! V1_10        2009/12/17 Jan-Peter Schulz
!  Optional interpolation of seaice variables
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler
!  Modifications to read grib data with grib_api
! V1_19        2012/06/06 Ulrich Schaettler, Susanne Brienen
!  Modifications to read_gme_grib and sending/receiving the grib data:
!    to fix the little endian problem and to send the full message
!  Usage of ke_soil_in, msoilgrib_in for incoming data
!  Horizontal interpolation to intermediate field w_so_gl; only if 
!   ke_soil_in == ke_soil_lm, set also field w_so_lm (SB)
!  Implemented conditional compilation for GRIBDWD
! V1_20        2012/09/03 Ulrich Schaettler
!  Enlarged strings for date variables to 14 characters
!  Adapted calls to subroutine make_fn
!  Renamed 'grax' to 'apix' to be conform with other models
! V1_21        2013/03/25 Ulrich Schaettler
!  Enlarged debug level value for certain outputs
!  Modifications to read GME GRIB2 data
! V1_22        2013/07/11 Ulrich Schaettler
!  Implemented necessary ifdef GRIBAPI
!  Adapted interfaces to read_gribapi with special grib_api integer
!  Renamed lprog_qrqs to lprog_qr_qs to be consistent with other names
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

#ifdef GRIBAPI
USE grib_api
#endif

! Modules used:
USE data_parameters , ONLY :  &

! Imported Parameters:
    ireals,    & ! KIND-type parameters for real variables
    irealgrib, & ! KIND-type of the REALs in the grib library
    iintegers, & ! KIND-type parameter for standard integer variables
    intgribc,  & ! KIND-type of the c decks in the grib library
    intgribf,  & ! KIND-type of the fortran decks in the grib library
    iwlength,  & ! length of integers used in the griblib in byte
    int_ga       ! integer precision for grib_api: length of message in bytes

!------------------------------------------------------------------------------

USE data_fields_lm, ONLY : &
  t_s_lm    ,      & ! temperature of the ground surface                (  K  )
  t_snow_lm ,      & ! temperature of the snow surface                  (  K  )
  t_m_lm    ,      & ! temperature between upper and medium soil layer  (  K  )
  t_cl_lm   ,      & ! temperature between medium and lower soil layer  (  K  )
  w_snow_lm ,      & ! water content of the snow                        (m H2O)
  w_i_lm    ,      & ! water content of the interception storage        (m H2O)
  w_g1_lm   ,      & ! water content of the upper soil layer            (m H2O)
  w_g2_lm   ,      & ! water content of the medium soil layer           (m H2O)
  w_g3_lm   ,      & ! water content of the lower soil layer            (m H2O)
                     ! (if nl_soil_lm = 3, unused otherwise)
  w_cl_lm   ,      & ! climatological deep soil water content           (m H2O)
  lolp_lm   ,      & ! Land Sea Mask of LM for 'M'atch Interpolation
  fis_gl    ,      & ! GME interpolated orography * G                   (m2/s2)
  ps_gl     ,      & ! surface pressure on the interpol. GME orogr.     ( Pa  )
  fic_gl    ,      & ! check level of geopotential                      (m2/s2)
  dpsdt_gl  ,      & ! surface pressure tendency                        (Pa/s )
  t_s_gl    ,      & ! temperature of the ground surface                (  K  )
  rh_s_gl   ,      & ! relative humidity at the surface                 (kg/kg)
  dtms_gl   ,      & ! t_m_lm    - t_s_lm                               (  K  )
  dtkes_gl  ,      & ! t(ke)_lm  - t_s_lm                               (  K  )
  dtssnow_gl         ! t_s_lm    - t_snow_lm                            (  K  )

USE data_fields_lm, ONLY : &
  index_m   ,      & !
  index_u   ,      & !
  index_v   ,      & !
  baryll_m  ,      & !
  baryll_u  ,      & !
  baryll_v  ,      & !
  rotang_m  ,      & !
  rotang_u  ,      & !
  rotang_v  ,      & !
  w_intpol     ,   & ! interpolation weights for gme2lm                 (  -  )
  n_intpol     ,   & ! nearest GME gridpoint for nearest neighbor interpolation
  m_intpol     ,   & ! nearest GME gridpoint with same lsm for match interpolation
  l_intpol     ,   & ! to use a far away GME gridpoint with same lsm
  lonlm_u   ,      & ! longitudes of the LM u grid points
  latlm_u   ,      & ! latitudes of the LM u grid points
  lonlm_v   ,      & ! longitudes of the LM v grid points
  latlm_v   ,      & ! latitudes of the LM v grid points
  u_lm      ,      & ! zonal wind speed                                 ( m/s )
  v_lm      ,      & ! meridional wind speed                            ( m/s )
  t_lm      ,      & ! temperature                                      (  K  )
  qv_lm     ,      & ! specific water vapor content                     (kg/kg)
  qc_lm     ,      & ! specific cloud water content                     (kg/kg)
  qi_lm     ,      & ! cloud ice content                                (kg/kg)
  t_so_lm   ,      & ! multi-layer soil temperature                     (  K  )
  dt_so_gl  ,      & ! multi-layer soil temperature                     (  K  )
  w_so_lm   ,      & ! multi-layer soil moisture                        (m H2O)
  w_so_gl   ,      & ! multi-layer soil moisture (for interpolation)    (m H2O)
  freshsnw_lm,     & ! weighting function indicating 'freshness' of snow
  rho_snow_lm        ! for prognostic treatment of snow density         (kg/m3)

!------------------------------------------------------------------------------

USE data_fields_in,  ONLY: &
  soiltyp_gme ,         & ! type of the soil (keys 0-9)                   --
  grd_glob    ,         & ! coefficients needed for the calculation of the
                          ! gradients in eta- and chi-direction          (1/m)
  eta_glob    ,         & !
  chi_glob    ,         & !
  cpsi_glob   ,         & ! cosine of the rotation angle between the local
                          ! coordinate system of the home node and the 6 (5)
                          ! neighbouring gridpoints
  spsi_glob   ,         & ! sine of the rotation angle between the local
                          ! coordinate system of the home node and the 6 (5)
                          ! neighbouring gridpoints
  lolp_gme                ! Land Sea Mask of GME for 'M'atch Interpolation

!------------------------------------------------------------------------------

USE data_grid_lm,    ONLY: &
  pollat,      & ! latitude of the rotated north pole (in degrees, N>0)
  pollon,      & ! longitude of the rotated north pole (in degrees, E>0)
  ie2lm,       & !
  je2lm,       & !
  ke_soil_lm     ! number of levels in multi-layer soil model in output

!------------------------------------------------------------------------------

USE data_grid_in,    ONLY: &
  ak_in ,      & ! vertical coordinate parameters for half levels
  bk_in ,      & !                  - " -
  akh_in,      & ! vertical coordinate parameters for main levels
  bkh_in,      & !                  - " -
  dak_in,      & ! differences between half levels
  dbk_in,      & !                  - " -
  ids,         & ! start index of diamonds (ids = 1)
  ide,         & ! end index of diamonds (ide = 10)
  ni_gme,      & ! resolution of GME
  i3e_gme,     & ! number of levels in the vertical
  nd,          & ! number of diamonds (nd = ide-ids+1 = 10)
  ke_soil_in,  & ! number of input levels in multi-layer soil model
  czmls_in,    & ! depth of the main soil layers in meters in output
  czhls_in,    & ! depth of the half soil layers in meters in output
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
  isp26,       & !
  ni2,         & ! ni_gme=3**ni3*2**ni2 with ni3 0 or 1 and ni2 > 1
  ni3            !

USE data_grid_in,    ONLY: &
  igg1s,       & ! start index of global array-dimension in x-direction
  igg1sm1,     & ! = igg1s - 1
  igg1sm2,     & ! = igg1s - 2
  igg1e,       & ! end index of global array-dimension in x-direction
  igg1ep1,     & ! = igg1e + 1
  igg1ep2,     & ! = igg1e + 2
  igg2s,       & ! start index of global array-dimension in y-direction
  igg2sm1,     & ! = igg2s - 1
  igg2sm2,     & ! = igg2s - 2
  igg2e,       & ! end index of global array-dimension in y-direction
  igg2ep1,     & ! = ig2e + 1
  igg2ep2,     & ! = igg2e + 2
  j1_min,      & ! smallest index in j1-direction for LM subdomain
  j1_max,      & ! biggest index in j1-direction for LM subdomain
  j2_min,      & ! smallest index in j2-direction for LM subdomain
  j2_max,      & ! biggest index in j2-direction for LM subdomain
  jd_min,      & ! smallest index for diamonds for a LM subdomain
  jd_max,      & ! biggest index for diamonds for a LM subdomain
  i1mrp,       & ! i1-index of the four mirrored points of the extended array
  i2mrp,       & ! i2-index of the four mirrored points of the extended array
  ispoke,      & ! offsets of the 6 (5) neighbouring gridpoints relative to
                 ! i1-direction use ispoke(m), m=1,6 (5), in i2-direction
                 ! use ispoke(m+6), m=1,6 (5); phys. dim. ( - )
  ispokes        ! offsets of the 6 neighbouring gridpoints relative
                 ! to the central node for the four special points

!------------------------------------------------------------------------------

USE data_int2lm_control,    ONLY: &
  noutput,         & ! unit number for output file
  linitial,        & ! if .TRUE., initial data for LM
  lboundaries,     & ! if .TRUE., lateral boundaries for LM
  itype_calendar,  & ! for specifying the calendar used
  yinput_model,   & ! string to identify the input model
  lcomp_bound,     & ! compute fields for boundaries
  itype_t_cl,      & ! to choose origin and treatment of deep soil temperature
  l_smi,           & ! if .TRUE., interpolate soil moisture with SMI
  lseaice,         & ! if .TRUE., run with sea ice model
  lmulti_layer_lm, & ! if .TRUE., compute soil fields for multi-layer soil model
  lmulti_layer_in, & ! if .TRUE., incoming data from new multi-layer soil model
  msoilgrib_in,    & !
  lprog_qi,        & ! if .TRUE., interpolate qi from GME to LM grid
  lprog_qr_qs,     & ! if .TRUE., interpolate qr,qs from GME to LM grid
  lprog_rho_snow,  & ! if .TRUE., interpolate rho_snow from GME to LM grid
  ltime,           & ! detailled timings of the program are given
  nl_soil_in,      & ! number of soil layers in GME
  nl_soil_lm,      & ! number of soil layers in LM, resp. HM
  kcontrol_fi,     & ! control level for geopotential
  dt,              & ! time step used in the LM
  timings,         & ! for storing the times for different parts of the program
  idbg_level,      & ! to control verbosity of output
  lprintdeb_all      ! whether all PEs print debug output

!------------------------------------------------------------------------------

USE data_int2lm_io,        ONLY : &
  nunit_of_time,     & ! indicator for unit-of-time (1hr, 15min, 30min,...)
  yin_cat,           & ! catalog-name of the GME files
  yin_lfn,           & ! name of the file with GME input data
  nuchkdat,          & ! checking the I/O data
  ract_hour,         & ! actual hour of actual day (returned by get_utc_date)
  yuchkdat,          & ! checking the I/O data
  ymode_read,        & ! mode for opening the (read) Grib files
  yin_form_read,     & ! input format of boundary data
  npds,              & ! Dimension for product definition section (pds)
  ngds,              & ! Dimension for grid description section (gds)
  nbms,              & ! Dimension for bit map section (bms)
  nbds,              & ! Dimension for binary data section
  ndsup,             & ! Dimension for dsup
  ndims,             & ! Dimension for idims (contains all dimensions)
  lfd,               & !
  nbitmap,           & !
  nbitmappres,       & ! to check presence of bitmap
  lds,               & !
  nvar_lm,           & ! maximum number of variables in LM variable table
  nvar_in,           & ! maximum number of variables in GME variable table
  idwdednr,          & ! grib edition number for dwdlib
  igrbednr,          & ! grib edition number for grib_api
  undefgrib,         & ! value for "undefined" in the grib routines
  undef                ! the same with other KIND-Parameter

USE data_int2lm_io,        ONLY : &
  ytunit_in,         & ! time unit for input data
  yinput_type,       & ! type of input data: 'forecast', 'analysis' or 'ana_init'
  ydate_ini,         & ! start of the forecast yyyymmddhh (year,month,day,hour)
  lmmss_ini,         & ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                       ! if .FALSE. 10 digits date format (YYYYMMDDHH)
                       ! for ydate_ini and result files of INT2LM
  ydate_bd,          & ! start of the forecast from boundary model from which
                       ! data are used
  lmmss_bd,          & ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                       ! if .FALSE. 10 digits date format (YYYYMMDDHH)
                       ! for ydate_bd  and input  files of INT2LM
  iblock,            & ! array for gribed data
  idims_in,          & ! array for all dimensions
  ibmap,             & ! array for
  ipds,              & ! product definition section
  igds_in,           & ! grid description section
  ibms,              & ! bit map section
  ibds,              & ! binary data section
  dsup,              & ! Parameter for grib routines
  var_lm,            & ! array for LM variable table
  var_in,            & ! array for GME variable table
  nincwait,          & ! if ready-file is not available wait nincwait seconds
                       ! until next attempt
  nmaxwait,          & ! if ready-file is not available after nmaxwait seconds,
                       ! abort the program
  ytrans_in,         & ! directory for reading ready-files
  ytrans_out,        & ! directory for writing ready-files
  lchkin,            & ! logical for print of check-values (max,min,mean)
                       ! of GME-fields
  lcheck_bmm,        & ! if .TRUE., check bitmaps for mass grid points
  lcheck_bmu,        & ! if .TRUE., check bitmaps for u-grid points
  lcheck_bmv           ! if .TRUE., check bitmaps for v-grid points

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY :  &
    lasync_io,       & ! if .TRUE.: the model runs with extra PEs for
                       ! asynchronous IO
    num_compute,     & ! number of compute PEs
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    icomm_cart,      & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_grib,        & ! determines the REAL type used for the GRIB library
    imp_integers,    & ! determines the correct INTEGER type used in the model
                       ! for MPI
    imp_logical,     & ! determines the correct LOGICAL   type used in the
                       ! model for MPI
    imp_byte           ! determines the correct BYTE type used in the model
                       ! for MPI

!------------------------------------------------------------------------------

USE data_int2lm_constants,     ONLY :  &
    R_d,     & ! gas constant for dry air                      [J/K*kg]
    Rdv,     & ! = R_d/R_v,
    Rvd_m_o, & ! = R_v/R_d - 1.0,
    O_m_rdv, & ! = 1. - Rdv
    r_earth, & ! mean radius of the earth 
    b1,      & !  a
    b2_w,    & !  b
    b3,      & !  c/b (0 degree Celsius [Kelvin])
    b4_w,    & !  d
    porb,    & ! pore volume of COSMO soil types
    pwpb,    & ! permanent wilting point of COSMO soil types
    fcb        ! field capacity of COSMO soil types

!------------------------------------------------------------------------------

USE mpe_io,              ONLY :  mpe_io_read
USE utilities,           ONLY :  elapsed_time, uv2uvrot_vec, diff_minutes
USE meteo_utilities,     ONLY :  qsat, psat_w
USE environment,         ONLY :  model_abort, comm_barrier
USE parallel_utilities,  ONLY :  remark, distribute_values, gather_values
USE gme_utilities,       ONLY :  pp_interp2ls, pp_interp2qs, pp_interp2qv, xd,&
                                 check_bitmap, init_gme_interpol
USE io_utilities,        ONLY :  open_file, close_file, make_fn,              &
                                 check_record, check_gme_grid

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================
!+ Organizes the reading and interpolation of GME-fields
!------------------------------------------------------------------------------

SUBROUTINE org_gme_interpol (nnow, lfirst, ydate)

!------------------------------------------------------------------------------
!
! Description:
!  GME-fields are read and distributed to the available PEs. Every PE gets a
!  different record for processing. The fields with dependencies for the 
!  interpolation are stored in memory first and treated after all fields have
!  been read.
!
! Method:
!  After creating the filename and opening the grib file, all records are read
!  in a read_loop. PE 0 determines whether the records are needed or not and
!  distributes them to the PEs. 
!
! Input files:
!  Grib-file with GME-fields. 
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
INTEGER  (KIND=iintegers),  INTENT(IN)  ::  &
  nnow                    ! actual time step to be processed

LOGICAL                  ,  INTENT(INOUT)  ::  &
  lfirst                  ! for actions that have to be done only once

CHARACTER (LEN=14), INTENT(IN)    ::  &
  ydate           ! actual date in the form   yyyymmddhhmmss

!------------------------------------------------------------------------------

! Local arrays
REAL (KIND=ireals)         ::  &
  field_glob(igg1sm2:igg1ep2,igg2sm2:igg2ep2, 1:10 ),            &
  field_glob2(igg1sm2:igg1ep2,igg2sm2:igg2ep2, 1:10 ),           &
  field_interpol((igg1ep2-igg1sm2+1)*(igg2ep2-igg2sm2+1)*(jd_max-jd_min+1)), &
  field_interpol2((igg1ep2-igg1sm2+1)*(igg2ep2-igg2sm2+1)*(jd_max-jd_min+1)) 

REAL (KIND=irealgrib)      ::  &
  ds_gme    ((igg1e-igg1s+1)*(igg2e-igg2s+1)*11),                &
  field_grib(igg1s:igg1e,igg2s:igg2e, 1:11 )
                          ! field read from the GRIB-file

! Local arrays for intermediate storage of some GME fields
REAL (KIND=ireals), ALLOCATABLE ::  &
  ps_gme     (:,:,:),   & !
  ti3e_gme   (:,:,:),   & !
  t_s_gme    (:,:,:),   & !
  t_m_gme    (:,:,:),   & !
  t_snow_gme (:,:,:),   & !
  t_cl_gme   (:,:,:),   & !
  qv_s_gme   (:,:,:),   & !
  w_g1_gme   (:,:,:),   & !
  w_g2_gme   (:,:,:),   & !
  w_g3_gme   (:,:,:),   & !
  w_cl_gme   (:,:,:),   & !
  t_so_gme   (:,:,:,:), & !
  w_so_gme   (:,:,:,:), & !
  freshsnw_gme(:,:,:),  & !
  rho_snow_gme(:,:,:)     !

! Local work arrays for the interpolation routines
REAL (KIND=ireals)         ::  &
  gme_loc1 (igg1sm1:igg1ep1, igg2sm1:igg2ep1, jd_min:jd_max),   & !
  gme_loc2 (igg1sm1:igg1ep1, igg2sm1:igg2ep1, jd_min:jd_max),   & !
  gme_loc3 (igg1sm1:igg1ep1, igg2sm1:igg2ep1, jd_min:jd_max),   & !
  gme_loc4 (igg1sm1:igg1ep1, igg2sm1:igg2ep1, jd_min:jd_max)

REAL (KIND=ireals)         ::  &
  zfih (ie2lm,je2lm,i3e_gme), & ! for computing the geopotential
  zfilb(ie2lm,je2lm),         & ! 
  zpla (ie2lm,je2lm),         & ! 
  zplb (ie2lm,je2lm)            ! 

INTEGER  (KIND=intgribf), ALLOCATABLE   ::  &
  ublock (:,:), vblock(:,:)   ! for storing u and v

! Local variables:
INTEGER  (KIND=intgribf)   ::  &
  ierrf,                      & ! error code for grib routines
  izbyta,                     & ! for unpacking GDS and PDS
  mzlevf, mzlevu, mzlevv,     & ! level of u and v to be processed
  iniyy, inimm, inidd, inihh, inimin, inisec, & !
  ibdyy, ibdmm, ibddd, ibdhh, ibdmin, ibdsec, & !
  imindif

INTEGER  (KIND=intgribc)   ::  &
  izdebug,              & ! verbosity of debug output
  ierrc,                & ! Error return code
  nzincwaitc              ! as nincwait, but for calling C-routine of grib lib.

INTEGER  (KIND=iintegers)  ::  &
  nufile,               & ! unit number of opened grib file
  izerror,              & ! status and error status variable
  niostat,              & ! status and error status variable
  izlen,                & ! length of path-name for ready-files
  nzwait,               & ! seconds waited so far
  mzlocpe(0:num_compute-1),   & ! location of records of every PE in variable table
  mzlevpe(0:num_compute-1),   & ! level of a multi-level field
  mzloc, mzlev,         & ! the same for my PE
  mzlocu_gme, mzlocv_gme,&! location of u and v in GME-variable table
  mzlocu_lm, mzlocv_lm, & ! location of u and v in GME-variable table
  npr, j1, j2, j3, j123,& ! additional variables
  nzuvblock,            & ! second dimension of ublock and vblock
  nzuv,                 & ! number of u- and v-levels this PE has got
  nzuvpe(0:num_compute-1),    & ! nzuv from all PEs
  izrest, izstat,       & !
  nzlevel, i,j,k, ij, jd, nzbyte, lfdgme, ldsgme, ndiff_ini_bd, ndiff_ini_bdref, &
  isecleft, ihour, imin, isec, igribid, ireturn, nb_field_in, inrpoints

REAL    (KIND=ireals)      ::  &
  zfactor, zbias, zfactoru, zbiasu, zfactorv, zbiasv, zrealdiff, ztv, rsecleft

LOGICAL                    ::  &
  lzrequired, & ! indicates whether a record from the grib file is required
  lzreqpe(0:num_compute-1),  & ! the same for every PE
  lzeof,      & ! indicates the end of file
  lzexist,    & ! to check, whether 'gif' or 'gff' or ready-files exists
  lzcheck       ! to check, whether all data have been read

LOGICAL                    ::  &
  lu(i3e_gme), lv(i3e_gme), lt(i3e_gme), lqc(i3e_gme), lqv(i3e_gme),      &
  lqi(i3e_gme), lqr(i3e_gme), lqs(i3e_gme),                               &
  lps, lt_cl, lt_s, lt_snow, lt_m, lqv_s, lw_cl, lw_i, lw_snow, lw_g1,    &
  lw_g2, lw_g3, ldpsdt, lfic,                                             &
  lt_so(0:ke_soil_in+1), lw_so(ke_soil_in+1), lfreshsnw, lrho_snow,       &
  lt_ice, lh_ice

CHARACTER (LEN= 14)        ::  &
  yzdate      ! actual date used for make_fn

CHARACTER (LEN= 14)        ::  &
  yzfulldate  ! actual date used for the check-routine check_gme_grid,
              ! including minutes and seconds

CHARACTER (LEN=250)        ::  &
  yzpath      ! full path and name of the GME-file

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN= 21)        ::  &
  yshortname  ! short name of the variable

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling

CHARACTER (LEN=  3)        ::  &
  yzhead      ! file name header for GME-file

!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1a: Initializations
!------------------------------------------------------------------------------

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  IF (izdebug > 9) THEN
    PRINT *, '    Starting org_gme_interpol'
  ENDIF

  ! determine actual minute and second of the step that has to be processed:
  ! this is needed in the check-routines, to decide whether the data that is
  ! read is really valid for this step (time)
!US  READ (ydate(9:10),'(I2)') ihour
!US  rsecleft = (ract_hour - REAL(ihour, ireals)) * 3600.0_ireals
!US  isecleft = INT (rsecleft, iintegers)
!US  imin     = isecleft / 60
!US  isec     = isecleft - (imin*60)
!US  yzfulldate(1:10) = ydate
!US  WRITE (yzfulldate(11:14), '(2I2.2)') imin, isec

  ! With the new 14 digit dates this should be:
  yzfulldate(1:14) = ydate(1:14)

  ! initialize level counter and logical flags
  lzeof       = .FALSE.
  lzrequired  = .FALSE.
  lzexist     = .FALSE.
  ierrc       = 0_intgribc
  ierrf       = 0_intgribf
  izbyta      = 9_intgribf
  yzroutine   = 'org_gme_interpol'
  nzincwaitc  = INT (nincwait, intgribc)
  undef       = REAL (undefgrib, ireals)
  izerror     = 0_iintegers
  mzloc       = 1_iintegers ! for safety reasons
  igrbednr    = 0_iintegers ! not 1 and not 2

  ! Set dimensions for grib variables and allocate iblock, ibmap and dsup
  ! (the role of ds is taken by field_grib).
  nzbyte = 8

  ldsgme  = (ni_gme+1) * (ni_gme+1) * 11
  lfdgme = ldsgme * nzbyte / iwlength + 2000

  lfd = INT (lfdgme, intgribf)
  lds = INT (ldsgme, intgribf)

  ! Correct idims_in 
  idims_in( 5)   = nbitmap
  idims_in( 7)   = lds
  idims_in( 8)   = lfd

  ALLOCATE (iblock(lfd), ibmap(nbitmap), STAT=izstat)
  ALLOCATE (dsup(ndsup),                 STAT=izstat)

  lu(:)       = .FALSE.
  lv(:)       = .FALSE.
  lt(:)       = .FALSE.
  lqv(:)      = .FALSE.
  lqc(:)      = .FALSE.
  IF (lprog_qi) THEN
    lqi(:)      = .FALSE.
  ELSE
    lqi(:)      = .TRUE.
  ENDIF
  IF (lprog_qr_qs) THEN
    lqr(:)      = .FALSE.
    lqs(:)      = .FALSE.
  ELSE
    lqr(:)      = .TRUE.
    lqs(:)      = .TRUE.
  ENDIF
  lfic        = .FALSE.
  lps         = .FALSE.
  lt_cl       = .FALSE.
  lt_s        = .FALSE.
  lt_snow     = .FALSE.
  lt_m        = .FALSE.
  lqv_s       = .FALSE.
  lw_cl       = .FALSE.
  lw_i        = .FALSE.
  lw_snow     = .FALSE.
  lw_g1       = .FALSE.
  lw_g2       = .FALSE.
  lw_g3       = .FALSE.
  ldpsdt      = .FALSE.
  lt_so(:)    = .FALSE.
  lw_so(:)    = .FALSE.
  lfreshsnw   = .FALSE.
  lrho_snow   = .FALSE.
  lt_ice      = .FALSE.
  lh_ice      = .FALSE.

  ! Allocation of memory for the soil variables
  ALLOCATE                                                             &
     (ps_gme    (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),       &
      ti3e_gme  (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),       &
      t_snow_gme(igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),       &
      qv_s_gme  (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max), STAT=izstat)

! allocate these fields always, because they are passed as subroutine arguments
! IF (.NOT. lmulti_layer_in) THEN
    ALLOCATE                                                           &
     (t_s_gme   (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),       &
      t_m_gme   (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),       &
      t_cl_gme  (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),       &
      w_g1_gme  (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),       &
      w_g2_gme  (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),       &
      w_g3_gme  (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),       &
      w_cl_gme  (igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max), STAT=izstat)
! ELSE
!  IF (.NOT. lcomp_bound) THEN
    ALLOCATE                                                                 &
     (t_so_gme  (igg1sm2:igg1ep2,igg2sm2:igg2ep2,0:ke_soil_in+1,jd_min:jd_max), &
      w_so_gme  (igg1sm2:igg1ep2,igg2sm2:igg2ep2,1:ke_soil_in+1,jd_min:jd_max), &
      freshsnw_gme(igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),           &
      rho_snow_gme(igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),           &
                                                                 STAT=izstat)
!  ELSE
!   ALLOCATE                                                                 &
!    (t_so_gme  (igg1sm2:igg1ep2,igg2sm2:igg2ep2,0:ke_soil_in+1,jd_min:jd_max), &
!                                                                STAT=izstat)
!  ENDIF
! ENDIF

  ! Allocation of memory for ublock and vblock
  ! Compute nzuvblock
  izrest = MOD (i3e_gme, num_compute)
  IF (izrest == 0) THEN
    nzuvblock = i3e_gme / num_compute
  ELSE
    nzuvblock = i3e_gme / num_compute + 1
  ENDIF
  nzuv      = 0
  nzuvpe(:) = 0

  ALLOCATE (ublock(lfd,nzuvblock), vblock(lfd,nzuvblock), STAT=izstat)

!------------------------------------------------------------------------------
! Section 1b: Check whether old boundary data are used
!------------------------------------------------------------------------------

  IF (ydate_bd /= ydate_ini) THEN
    ! Compute the difference between start date and date of boundary data,
    ! if a difference is present
    READ( ydate_ini,'(I4,5I2)' ) iniyy, inimm, inidd, inihh, inimin, inisec
    READ( ydate_bd ,'(I4,5I2)' ) ibdyy, ibdmm, ibddd, ibdhh, ibdmin, ibdsec
    CALL diff_minutes( iniyy, inimm, inidd, inihh, inimin,                  &
                       ibdyy, ibdmm, ibddd, ibdhh, ibdmin,                  &
                       itype_calendar, imindif, ierrf )

!US    ndiff_ini_bd =  - NINT( imindif * 60 / dt , iintegers )
!US as in COSMO-Model
    ndiff_ini_bd =  - NINT( (imindif * 60 + (ibdsec - inisec)) / dt , iintegers )
    ! correct ndiff_ini_bd if ydate_ini not equal to reference time of
    ! boundary data
    ! account for seconds here!
    ndiff_ini_bdref = NINT((inimin * 60) / dt, iintegers)
    ndiff_ini_bd = ndiff_ini_bd - ndiff_ini_bdref

    yzdate = ydate
  ELSE
    ndiff_ini_bd = 0
    yzdate = ydate
  ENDIF

  IF (izdebug > 9) THEN
    PRINT *, '    Initializations ready for org_gme_interpol'
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Check if ready-files are available
!------------------------------------------------------------------------------

  ! If required, check whether ready-files are available
  IF (ytrans_in /= '   ') THEN
    IF (izdebug > 9) THEN
      PRINT *, '    Check ready files'
    ENDIF

    IF (my_cart_id == 0) THEN
      izlen   = LEN_TRIM(ytrans_in)
      nzwait  = 0
      lzexist = .FALSE.

      ! Create the file name for the GME ready-file
      IF     ( (yinput_type == 'forecast') .OR. (yinput_type == 'ana_init') ) THEN
        yzhead = 'GME'
      ELSEIF (yinput_type == 'analysis') THEN
        yzhead = 'GMA'
      ENDIF

      CALL make_fn (yzhead, yzdate, ydate_ini, 'f',' ', nnow+ndiff_ini_bd, dt, .TRUE.,    &
                    itype_calendar, ytrans_in, yzpath, lmmss_bd, izdebug, izerror)

#ifdef GRIBDWD
      ! The filename for ytunit='f' did not change, so this check has to work
      DO WHILE ( (.NOT. lzexist) .AND. (nzwait <= nmaxwait) )
        INQUIRE (FILE=yzpath, EXIST=lzexist)
        IF (.NOT. lzexist) THEN
          ! wait nincwait seconds and try again
          PRINT *, 'ready-file not available: ', yzpath(1:izlen+15)
          PRINT *, '      sleep ', nincwait,' seconds'
          CALL fsleep (nzincwaitc)
          nzwait = nzwait + nincwait
        ENDIF
      ENDDO
#endif
    ENDIF
 
    IF (num_compute > 1) THEN
      ! Synchronize the processes again
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)

      ! Distribute lzexist to all nodes
      CALL distribute_values (lzexist, 1, 0, imp_logical, icomm_cart, izerror)
    ENDIF

    IF (.NOT. lzexist) THEN
      yzerrmsg  = ' *** ERROR:  ready-file not available: '//yzpath(1:izlen+15)
      CALL model_abort (my_cart_id, 2006, yzerrmsg, yzroutine)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Create grib file name and open the file
!------------------------------------------------------------------------------

  ! Determine file header
  IF (izdebug > 9) THEN
    PRINT *, '    Create file name and open file'
  ENDIF

  IF     (yinput_type == 'forecast') THEN
    IF ((nnow == 0) .AND. (ydate_ini == ydate_bd)) THEN
      ! header could be of the form 'gif' or 'gff'
      ! Test whether a "gif" file is present, if not take a gff-file
      yzhead = 'gif'
      yzpath = yin_cat(1:LEN_TRIM(yin_cat))//'giff00000000'
      INQUIRE (FILE=yzpath, EXIST=lzexist)
      IF (.NOT. lzexist) THEN
        yzhead = 'gff'
      ENDIF
    ELSE
      ! for all other times, the header is gff
      yzhead = 'gff'
    ENDIF
  ELSEIF (yinput_type == 'analysis') THEN
    yzhead = 'gaf'
  ELSEIF (yinput_type == 'ana_init') THEN
    yzhead = 'gif'
  ENDIF

  ! Create the file name
  ! Here, lmmss_bd determines the format of the input file (depending on ytunit_in)
  CALL make_fn (yzhead, yzdate, ydate_ini, ytunit_in,' ', nnow+ndiff_ini_bd, dt, .TRUE., &
                itype_calendar, yin_cat, yzpath, lmmss_bd, izdebug, izerror)

  IF (izdebug > 5) THEN
    PRINT *, '  Call open_file for ', TRIM(yzpath)
  ENDIF

  ! All processors have to call the routine open_file. What the parallel
  ! program really does is determined in the routine.
  CALL open_file(nufile, yzpath, ymode_read, yin_form_read, icomm_cart, &
                 my_cart_id, num_compute, lasync_io, idbg_level,        &
                 yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, izerror, yzerrmsg, 'open_file')
  ELSE 
    lzeof = .FALSE.
  ENDIF

  ! Write headline in file YUCHKDAT
  IF (lchkin .AND. my_cart_id == 0) THEN
    IF (izdebug > 40) THEN
      PRINT *, '  Open file ', yuchkdat, ' again'
    ENDIF

    OPEN(nuchkdat, FILE=yuchkdat, FORM='FORMATTED', STATUS='UNKNOWN',  &
                   POSITION='APPEND', IOSTAT=niostat)
    IF(niostat /= 0) THEN
      izerror  = 1
      yzerrmsg = ' ERROR    *** Error while opening file YUCHKDAT *** '
      RETURN
    ENDIF

    ! Write a headline in YUCHKDAT for GME files     
    WRITE (nuchkdat,'(A)') '       '
    WRITE (nuchkdat,'(A)') '       '
    WRITE (nuchkdat,'(A,A)') 'Check GME file:  ', yzpath(1:LEN_TRIM(yzpath))
    WRITE (nuchkdat,'(A,I5,A,I5,A,I5)')                                  &
      '    ni_gme =', ni_gme, '  ni2 =', ni2, '  ni3 =', ni3
    WRITE (nuchkdat,'(A)') '    '
    WRITE (nuchkdat,'(A,A)')                                             &
      '   var       ee  lev        min     ',                            &
      '   location            max        location           mean'
    WRITE (nuchkdat,'(A,A)')                                             &
      '                                    ',                            &
      ' j1   j2   jd                   j1   j2   jd'
    WRITE (nuchkdat,'(A)') '  '
  ENDIF

!------------------------------------------------------------------------------
! Section 4: (Endless) loop over all records in the grib file
!------------------------------------------------------------------------------
 
  IF (izdebug > 9) THEN
    PRINT *, '    Start endless loop over all records'
  ENDIF

  read_loop: DO WHILE (.NOT. lzeof)

  !----------------------------------------------------------------------------
  ! 4.1: Get a record
  !----------------------------------------------------------------------------

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (12) = timings(12) + zrealdiff
    ENDIF

    ! Every PE gets one record from the file in an ordered manner
    ! (determined by the rank in the communicator icomm_rank). How this is
    ! done exactly is determined in the routine read_"format". This routine has
    ! to be called by every PE.

    ! The routine read_gme_grib reads and distributes the records until every
    ! PE has a gribed record in its iblock
    CALL read_gme_grib (nufile, mzlocpe, mzlevpe, lzreqpe, lzeof,          &
                        ublock, vblock, lfdgme, nzuvblock, nzuv, izerror)

    mzloc      = mzlocpe(my_cart_id)
    mzlev      = mzlevpe(my_cart_id)
    mzlevf     = INT (mzlevpe(my_cart_id), intgribf)
    lzrequired = lzreqpe(my_cart_id)

    IF (idbg_level > 15) THEN
      PRINT *, '    Read record ', var_in(mzloc)%name, '  ', mzloc, mzlev, lzrequired
    ENDIF

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (6) = timings(6) + zrealdiff
    ENDIF

  !----------------------------------------------------------------------------
  ! 4.2: Unpack and check the record
  !----------------------------------------------------------------------------

    IF (lzrequired) THEN
      IF (idbg_level > 15) THEN
        PRINT *, '    Unpack and check the record'
      ENDIF

      IF (yin_form_read == 'grb1') THEN
#ifdef GRIBDWD
        yshortname = TRIM(var_in(mzloc)%name)
        igrbednr   = 1

        ! Unpack the record
        CALL grbin1 (idwdednr, undefgrib, ndims, idims_in, iblock, ibmap, ipds,  &
                     igds_in, ibms, ibds, dsup, ds_gme, ierrf)
        IF ( ierrf /= 0 ) THEN
          yzerrmsg = 'Error in grbin1'
          CALL model_abort (my_cart_id, 4100, yzerrmsg, yzroutine)
        ENDIF

        ! get unit of time range: this is distributed after the read-loop from
        ! PE 0 to all other PEs, because it could be that not all PEs really
        ! get a record for input. Note that it is assumed that all records in
        ! this grib file do have the same unit of time.
        nunit_of_time = INT (ipds(16), iintegers)
#endif

      ELSEIF (yin_form_read == 'apix') THEN
#ifdef GRIBAPI

        ! Build the grib handle
        CALL grib_new_from_message (igribid, iblock, ireturn)
        IF (ireturn /= GRIB_SUCCESS) THEN
          PRINT *,   ' *** Error in grib_api grib_new_from_message  ', ireturn
          yzerrmsg = ' *** Error in grib_api grib_new_from_message'
          CALL model_abort (my_cart_id, 2035, yzerrmsg, yzroutine)
        ENDIF

        ! get edition number and shortname
        CALL grib_get (igribid, 'editionNumber', igrbednr,      ireturn)
        CALL grib_get (igribid, 'shortName',     yshortname,    ireturn)
        CALL grib_get (igribid, 'stepUnits',     nunit_of_time, ireturn)

        ! check bitmap and set missing value
!       CALL grib_get (igribid, 'bitmapPresent', nbitmappres,   ireturn)
        CALL grib_set (igribid, 'missingValue',  undefgrib,     ireturn)

        nb_field_in = (igg1e-igg1s+1)*(igg2e-igg2s+1)*11
        CALL grib_get_size(igribid, 'values', inrpoints, ierrf)
        IF (inrpoints > nb_field_in) THEN
          PRINT *, ' *** ERROR: size of message is too big for allocated field: ', inrpoints, nb_field_in
          izerror  = 2
          yzerrmsg = 'Error in grib_get_size'
          RETURN
        ENDIF
        CALL grib_get(igribid, 'values', ds_gme)
#endif
      ENDIF

      ! convert field to the precision used in GME2LM
      ij = 0
      DO jd = 1, 10
        DO j2 = igg2s, igg2e
          DO j1 = igg1s, igg1e
            ij = ij + 1
            IF (ds_gme(ij) /= undefgrib) THEN
              field_glob(j1,j2,jd) = REAL (ds_gme(ij), ireals)
            ELSE
              field_glob(j1,j2,jd) = undef
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ! Extend the field on the boundary lines
      CALL xd(field_glob,igg1sm2, igg1ep2, igg2sm2, igg2ep2, 1, 1, 1, 1, 1,  &
              10, igg1s, igg1e, igg2s, igg2e, 1, 10, 2, .TRUE., .FALSE.,     &
              field_glob,igg1sm2, igg1ep2, igg2sm2, igg2ep2, 10, undef, izerror)

      IF (lcheck_bmm) THEN
        ! Check whether the bitmap area is covering the LM domain
        IF (izdebug >= 2) THEN
          PRINT *, '    Check bitmap for mass grid points'
        ENDIF
        CALL check_bitmap (field_glob, index_m, ispoke, undef,               &
                           igg1sm2, igg1ep2, igg2sm2, igg2ep2, ie2lm, je2lm, &
                           izerror)
        IF (izerror /= 0) THEN
          yzerrmsg = 'Undefined GME values: Bitmap not set properly'
          CALL model_abort (my_cart_id, 4200, yzerrmsg, yzroutine)
        ENDIF
        lcheck_bmm = .FALSE.
      ENDIF

    ENDIF

    CALL check_gme_grid (igrbednr, igds_in, ngds, ipds, npds, igribid,      &
             yin_form_read, yshortname, yzfulldate,                         &
             i3e_gme, ni_gme, ni2, ni3, nd, (.NOT.lzeof) .AND. (lzrequired),&
             num_compute, icomm_cart, my_cart_id, .TRUE., itype_calendar,   &
             yinput_model, yzerrmsg, izerror)

    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 3300, yzerrmsg, yzroutine)
    ENDIF

    IF (lfirst) THEN
      CALL get_gme_vert (igribid)
    ENDIF

    ! Check the record
    IF (lchkin) THEN
      ! for safety reasons set mzloc for non required records
      IF (.NOT. lzrequired) mzloc = 1
      CALL check_record                                                    &
               (ds_gme,     igg1s, igg1e, igg2s, igg2e, ids, ide,          &
                igg1s  , igg1e  , igg2s  , igg2e  , ids, ide,  undefgrib,  &
                var_in(mzloc)%name, var_in(mzloc)%ee, mzlevf,              &
                lzrequired, nuchkdat, num_compute, icomm_cart,             &
                my_cart_id, yzerrmsg, izerror)
    ENDIF

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (7) = timings(7) + zrealdiff
    ENDIF

  !----------------------------------------------------------------------------
  ! 4.3: Scale the field
  !----------------------------------------------------------------------------

    IF (lzrequired) THEN
      zbias     = var_in(mzloc)%bias
      zfactor   = var_in(mzloc)%factor
    ENDIF
   
  !----------------------------------------------------------------------------
  ! 4.4: Distribute and interpolate the field
  !----------------------------------------------------------------------------

    DO npr = 0, num_compute-1
      IF (lzreqpe(npr)) THEN
        IF ( (idbg_level > 15) .AND. (my_cart_id == npr) ) THEN
          PRINT *, '    Distribute field ', var_in(mzlocpe(npr))%name,     &
                   ' from PE ', npr
        ENDIF

        ! Distribute the global field to all PEs, so that every PE can 
        ! interpolate this field to its LM subdomain
        IF (my_cart_id == npr) THEN
          j123 = 0
          DO j3 = jd_min, jd_max
            DO j2 = igg2sm2, igg2ep2
              DO j1 = igg1sm2, igg1ep2
                j123 = j123 + 1
                IF (field_glob(j1,j2,j3) /= undef) THEN
                  field_interpol(j123) = field_glob(j1,j2,j3) / zfactor - zbias
                ELSE
                  field_interpol(j123) = undef
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
 
        izlen = (igg1ep2-igg1sm2+1) * (igg2ep2-igg2sm2+1) * (jd_max-jd_min+1)
        IF (num_compute > 1) THEN
          CALL distribute_values (field_interpol, izlen, npr, imp_reals,    &
                                  icomm_cart, izerror)
        ENDIF
 
        IF (lfirst) THEN
          lfirst = .FALSE.
          ! Initialize the interpolation from GME fields
          CALL init_gme_interpol                                                 &
                   (lolp_gme, field_interpol,                                    &
                    igg1sm2,  igg1ep2, igg2sm2, igg2ep2, jd_min, jd_max,         &
                    lolp_lm, w_intpol, n_intpol, m_intpol, l_intpol,             &
                    ispoke, baryll_m, index_m, 1, ie2lm, 1, je2lm,               &
                    undef, izdebug, yzerrmsg, izerror)
        ENDIF

        IF (ltime) THEN
          CALL elapsed_time (zrealdiff)
          timings (9) = timings(9) + zrealdiff
        ENDIF

        ! Now do the interpolation on the LM subdomain and store LM values
        ! in memory
        CALL interpol_gme (field_interpol, mzlocpe(npr), mzlevpe(npr),       &
                           gme_loc1, gme_loc2,                               &
                           ps_gme, ti3e_gme, t_s_gme, t_m_gme, t_snow_gme,   &
                           qv_s_gme, w_g1_gme, w_g2_gme, w_g3_gme, w_cl_gme, &
                           t_cl_gme, t_so_gme, w_so_gme, freshsnw_gme,       &
                           rho_snow_gme,                                     &
                           lu, lv, lt, lqc, lqv, lqi, lqr, lqs, lps, lt_cl,  &
                           lt_s, lt_snow, lt_m, lqv_s, lw_cl, lw_i, lw_snow, &
                           lw_g1, lw_g2, lw_g3, ldpsdt, lfic, lt_so, lw_so,  &
                           lfreshsnw, lrho_snow, lt_ice, lh_ice)

        IF (ltime) THEN
          CALL elapsed_time (zrealdiff)
          timings (8) = timings(8) + zrealdiff
        ENDIF

      ENDIF
    ENDDO

#ifdef GRIBAPI
    CALL grib_release(igribid)
#endif
  ENDDO read_loop

!------------------------------------------------------------------------------
! Section 5: Interpolate special fields
!------------------------------------------------------------------------------

  CALL interpol_gme_special (ps_gme, ti3e_gme, t_s_gme, t_m_gme,             &
                        t_snow_gme, t_cl_gme, qv_s_gme, w_g1_gme, w_g2_gme,  &
                        w_g3_gme, w_cl_gme, t_so_gme, w_so_gme, field_glob,  &
                        lu, lv, lt, lqc, lqv, lqi, lps, lt_cl, lt_s, lt_snow,&
                        lt_m, lqv_s, lw_cl, lw_i, lw_snow, lw_g1, lw_g2,     &
                        lw_g3, lt_so, lw_so, ldpsdt)

!------------------------------------------------------------------------------
! Section 6: Interpolate u and v
!------------------------------------------------------------------------------

  ! Look up the location of u and v in GME- and LM-variable table
  DO i=1,nvar_lm
    IF ((var_lm(i)%name == 'U_lm       ') ) THEN
      mzlocu_lm = i
    ENDIF
    IF ((var_lm(i)%name == 'V_lm       ') ) THEN
      mzlocv_lm = i
    ENDIF
  ENDDO

  DO i=1,nvar_in
    IF ((var_in(i)%name == 'U         ') ) THEN
      mzlocu_gme = i
    ENDIF
    IF ((var_in(i)%name == 'V         ') ) THEN
      mzlocv_gme = i
    ENDIF
  ENDDO

  nzlevel = 0

  DO WHILE (nzlevel <= nzuvblock)
    nzlevel = nzlevel + 1

  !----------------------------------------------------------------------------
  ! Section 6.1: Unpack u- and v-record
  !----------------------------------------------------------------------------

    IF (nzlevel <= nzuv) THEN
      ! unpack ublock
      IF (yin_form_read == 'grb1') THEN
#ifdef GRIBDWD
        igrbednr = 1
        CALL grbin1 (idwdednr, undefgrib, ndims, idims_in, ublock(1,nzlevel),    &
                     ibmap, ipds, igds_in, ibms, ibds, dsup, ds_gme, ierrf)
        IF ( ierrf /= 0 ) THEN
          yzerrmsg = 'Error in grbin1'
          CALL model_abort (my_cart_id, 4100, yzerrmsg, yzroutine)
        ENDIF
#endif
      ELSEIF (yin_form_read == 'apix') THEN
#ifdef GRIBAPI
        ! Build the grib handle
        CALL grib_new_from_message (igribid, ublock(:,nzlevel), ireturn)
        IF (ireturn /= GRIB_SUCCESS) THEN
          PRINT *,   ' *** Error in grib_api grib_new_from_message  ', ireturn
          yzerrmsg = ' *** Error in grib_api grib_new_from_message'
          CALL model_abort (my_cart_id, 2035, yzerrmsg, yzroutine)
        ENDIF

        CALL grib_get (igribid, 'editionNumber', igrbednr,     ireturn)

        ! check bitmap and set missing value
!       CALL grib_get (igribid, 'bitmapPresent', nbitmappres,   ireturn)
        CALL grib_set (igribid, 'missingValue',  undefgrib,     ireturn)

        nb_field_in = (igg1e-igg1s+1)*(igg2e-igg2s+1)*11
        CALL grib_get_size(igribid, 'values', inrpoints, ierrf)
        IF (inrpoints > nb_field_in) THEN
          PRINT *, ' *** ERROR: size of message is too big for allocated field: ', inrpoints, nb_field_in
          izerror  = 2
          yzerrmsg = 'Error in grib_get_size'
          RETURN
        ENDIF
        CALL grib_get(igribid, 'values', ds_gme)
#endif
      ENDIF

      ! convert field to the precision used in GME2LM
      ij = 0
      DO jd = 1, 10
        DO j2 = igg2s, igg2e
          DO j1 = igg1s, igg1e
            ij = ij + 1
            IF (ds_gme(ij) /= undefgrib) THEN
              field_glob(j1,j2,jd) = REAL (ds_gme(ij), ireals)
            ELSE
              field_glob(j1,j2,jd) = undef
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ! Extend the field on the boundary lines
      CALL xd(field_glob,igg1sm2, igg1ep2, igg2sm2, igg2ep2, 1, 1, 1, 1, 1,  &
              10, igg1s, igg1e, igg2s, igg2e, 1, 10, 2, .TRUE., .FALSE.,     &
              field_glob,igg1sm2, igg1ep2, igg2sm2, igg2ep2, 10, undef, izerror)

      IF (lcheck_bmu) THEN
        ! Check whether the bitmap area is covering the LM domain
        IF (izdebug >= 2) THEN
          PRINT *, '    Check bitmap for u grid points'
        ENDIF
        CALL check_bitmap (field_glob, index_u, ispoke, undef,               &
                           igg1sm2, igg1ep2, igg2sm2, igg2ep2, ie2lm, je2lm, &
                           izerror)
        IF (izerror /= 0) THEN
          yzerrmsg = 'Undefined GME values: Bitmap not set properly'
          CALL model_abort (my_cart_id, 4200, yzerrmsg, yzroutine)
        ENDIF
        lcheck_bmu = .FALSE.
      ENDIF

    ENDIF

    CALL check_gme_grid (igrbednr, igds_in, ngds, ipds, npds, igribid,      &
             yin_form_read, 'U',        yzfulldate,                         &
             i3e_gme, ni_gme, ni2, ni3, nd, nzlevel <= nzuv,                &
             num_compute, icomm_cart, my_cart_id, .TRUE., itype_calendar,   &
             yinput_model, yzerrmsg, izerror)

    mzlevu = (nzlevel-1)*num_compute + my_cart_id + 1

    ! Check the record
    IF (lchkin) THEN
      CALL check_record                                                    &
               (ds_gme,     igg1s, igg1e, igg2s, igg2e, ids, ide,          &
                igg1s  , igg1e  , igg2s  , igg2e  , ids, ide,  undefgrib,  &
                var_in(mzlocu_gme)%name, var_in(mzlocu_gme)%ee, mzlevu,    &
                nzlevel <= nzuv, nuchkdat, num_compute, icomm_cart,        &
                my_cart_id, yzerrmsg, izerror)
    ENDIF

#ifdef GRIBAPI
    IF (nzlevel <= nzuv) THEN
      CALL grib_release (igribid)
    ENDIF
#endif

    IF (nzlevel <= nzuv) THEN
      ! unpack vblock
      IF (yin_form_read == 'grb1') THEN
#ifdef GRIBDWD
        igrbednr = 1
        CALL grbin1 (idwdednr, undefgrib, ndims, idims_in, vblock(1,nzlevel),    &
                     ibmap, ipds, igds_in, ibms, ibds, dsup, ds_gme,     ierrf)
        IF ( ierrf /= 0 ) THEN
          yzerrmsg = 'Error in grbin1'
          CALL model_abort (my_cart_id, 4100, yzerrmsg, yzroutine)
        ENDIF
#endif
      ELSEIF (yin_form_read == 'apix') THEN
#ifdef GRIBAPI
        ! Build the grib handle
        CALL grib_new_from_message (igribid, vblock(:,nzlevel), ireturn)
        IF (ireturn /= GRIB_SUCCESS) THEN
          PRINT *,   ' *** Error in grib_api grib_new_from_message  ', ireturn
          yzerrmsg = ' *** Error in grib_api grib_new_from_message'
          CALL model_abort (my_cart_id, 2035, yzerrmsg, yzroutine)
        ENDIF

        CALL grib_get (igribid, 'editionNumber', igrbednr,     ireturn)

        ! check bitmap and set missing value
!       CALL grib_get (igribid, 'bitmapPresent', nbitmappres,   ireturn)
        CALL grib_set (igribid, 'missingValue',  undefgrib,     ireturn)

        nb_field_in = (igg1e-igg1s+1)*(igg2e-igg2s+1)*11
        CALL grib_get_size(igribid, 'values', inrpoints, ierrf)
        IF (inrpoints > nb_field_in) THEN
          PRINT *, ' *** ERROR: size of message is too big for allocated field: ', inrpoints, nb_field_in
          izerror  = 2
          yzerrmsg = 'Error in grib_get_size'
          RETURN
        ENDIF
        CALL grib_get(igribid, 'values', ds_gme)
#endif
      ENDIF

      ! convert field to the precision used in GME2LM
      ij = 0
      DO jd = 1, 10
        DO j2 = igg2s, igg2e
          DO j1 = igg1s, igg1e
            ij = ij + 1
            IF (ds_gme(ij) /= undefgrib) THEN
              field_glob2(j1,j2,jd) = REAL (ds_gme(ij), ireals)
            ELSE
              field_glob2(j1,j2,jd) = undef
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL xd(field_glob2,igg1sm2, igg1ep2, igg2sm2, igg2ep2, 1, 1, 1, 1, 1, &
              10, igg1s, igg1e, igg2s, igg2e, 1, 10, 2, .TRUE., .FALSE.,     &
              field_glob2,igg1sm2, igg1ep2, igg2sm2, igg2ep2, 10, undef, izerror)

      IF (lcheck_bmv) THEN
        ! Check whether the bitmap area is covering the LM domain
        IF (izdebug >= 2) THEN
          PRINT *, '    Check bitmap for v grid points'
        ENDIF
        CALL check_bitmap (field_glob2, index_v, ispoke, undef,              &
                           igg1sm2, igg1ep2, igg2sm2, igg2ep2, ie2lm, je2lm, &
                           izerror)
        IF (izerror /= 0) THEN
          yzerrmsg = 'Undefined GME values: Bitmap not set properly'
          CALL model_abort (my_cart_id, 4200, yzerrmsg, yzroutine)
        ENDIF
        lcheck_bmv = .FALSE.
      ENDIF

    ENDIF

    CALL check_gme_grid (igrbednr, igds_in, ngds, ipds, npds, igribid,      &
             yin_form_read, 'V',        yzfulldate,                         &
             i3e_gme, ni_gme, ni2, ni3, nd, nzlevel <= nzuv,                &
             num_compute, icomm_cart, my_cart_id, .TRUE., itype_calendar,   &
             yinput_model, yzerrmsg, izerror)

    mzlevv = (nzlevel-1)*num_compute + my_cart_id + 1

    ! Check the record
    IF (lchkin) THEN
      CALL check_record                                                    &
               (ds_gme,     igg1s, igg1e, igg2s, igg2e, ids, ide,          &
                igg1s  , igg1e  , igg2s  , igg2e  , ids, ide,  undefgrib,  &
                var_in(mzlocv_gme)%name, var_in(mzlocv_gme)%ee, mzlevv,    &
                nzlevel <= nzuv, nuchkdat, num_compute, icomm_cart,        &
                my_cart_id, yzerrmsg, izerror)
    ENDIF

#ifdef GRIBAPI
    IF (nzlevel <= nzuv) THEN
      CALL grib_release (igribid)
    ENDIF
#endif

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (7) = timings(7) + zrealdiff
    ENDIF

  !----------------------------------------------------------------------------
  ! Section 6.2: get bias and factor
  !----------------------------------------------------------------------------

    zbiasu    = var_in(mzlocu_gme)%bias
    zfactoru  = var_in(mzlocu_gme)%factor
   
    zbiasv    = var_in(mzlocv_gme)%bias
    zfactorv  = var_in(mzlocv_gme)%factor

  !----------------------------------------------------------------------------
  ! Section 6.3: Distribute and interpolate the records
  !----------------------------------------------------------------------------

    ! Gather nzuv from all PEs
    IF (num_compute > 1) THEN
      CALL gather_values (nzuv, nzuvpe, 1, num_compute, imp_integers, -1,    &
                          icomm_cart, yzerrmsg, izerror)
    ELSE
      nzuvpe(0) = nzuv
    ENDIF

    ! Loop over all PEs
    DO npr = 0, num_compute-1
      IF (nzlevel <= nzuvpe(npr)) THEN
        ! Distribute the global field to all PEs, so that every PE can
        ! interpolate this field to its LM subdomain
        IF (my_cart_id == npr) THEN
          j123 = 0
          DO j3 = jd_min, jd_max
            DO j2 = igg2sm2, igg2ep2
              DO j1 = igg1sm2, igg1ep2
                j123 = j123 + 1
                IF (field_glob(j1,j2,j3) /= undef) THEN
                  field_interpol (j123) = field_glob (j1,j2,j3) / zfactoru-zbiasu
                  field_interpol2(j123) = field_glob2(j1,j2,j3) / zfactorv-zbiasv
                ELSE
                  field_interpol (j123) = undef
                  field_interpol2(j123) = undef
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        izlen = (igg1ep2-igg1sm2+1) * (igg2ep2-igg2sm2+1) * (jd_max-jd_min+1)
        IF (num_compute > 1) THEN
          CALL distribute_values (field_interpol , izlen, npr, imp_reals, &
                                  icomm_cart, izerror)
          CALL distribute_values (field_interpol2, izlen, npr, imp_reals, &
                                  icomm_cart, izerror)
        ENDIF

        IF (ltime) THEN
          CALL elapsed_time (zrealdiff)
          timings (9) = timings(9) + zrealdiff
        ENDIF

        ! Determine which level of the 3D GME-field is processed: the
        ! organization here is that the levels are processed in increasing
        ! order:
        mzlev = (nzlevel-1)*num_compute + npr + 1

        ! Now do the interpolation on the LM subdomain and store LM values
        ! in memory
        CALL interpol_gme_uv (field_interpol, field_interpol2,              &
                        mzlocu_gme, mzlocv_gme, mzlocu_lm, mzlocv_lm,       &
                        mzlev, gme_loc1, gme_loc2, gme_loc3, gme_loc4, lu, lv)

        IF (ltime) THEN
          CALL elapsed_time (zrealdiff)
          timings (8) = timings(8) + zrealdiff
        ENDIF
      ENDIF
    ENDDO
  ENDDO  ! of the DO WHILE loop

!------------------------------------------------------------------------------
! Section 7: Compute a geopotential control level
!------------------------------------------------------------------------------
 
  IF ((.NOT. lfic) .AND. (ALL(lt (:))) .AND. (ALL(lqc(:))) .AND.      &
                         (ALL(lqv(:))) .AND. lps) THEN
    ! Compute the (interpolated) geopotential on half levels
    DO k = i3e_gme, 1, -1
      DO j = 1, je2lm
        DO i = 1, ie2lm
          IF (k == i3e_gme) THEN
            zfilb(i,j) = fis_gl(i,j)
            zplb (i,j) = ps_gl (i,j)
          ENDIF

          ! Compute the virtual temperature
          ztv = t_lm(i,j,k) * (1.0_ireals + Rvd_m_o*qv_lm(i,j,k) - qc_lm(i,j,k))

          ! Compute geopotential on half levels
          IF (k == 1) THEN
            zfih (i,j,k) = zfih(i,j,k+1) + R_d * ztv * LOG(2.0_ireals)
          ELSE
            zpla (i,j)   = ak_in(k) + bk_in(k) * ps_gl(i,j)
            zfih (i,j,k) = zfilb(i,j) + R_d * ztv * LOG(zplb(i,j)/zpla(i,j))
          ENDIF
        ENDDO
      ENDDO
 
      ! Swap local arrays
      DO j = 1, je2lm
        DO i = 1, ie2lm
          zplb (i,j) = zpla(i,j)
          zfilb(i,j) = zfih(i,j,k)
        ENDDO
      ENDDO
    ENDDO

    fic_gl(:,:) = zfih(:,:,kcontrol_fi)
    lfic = .TRUE.
  ENDIF

!------------------------------------------------------------------------------
! Section 8: Cleanup
!------------------------------------------------------------------------------
 
  ! Deallocate the local variables and the grib fields
  DEALLOCATE (ublock, vblock, iblock, ibmap, dsup)
  DEALLOCATE (ps_gme, ti3e_gme, t_snow_gme, qv_s_gme)

  IF (.NOT. lmulti_layer_in) THEN
    DEALLOCATE (t_s_gme, t_m_gme, t_cl_gme, w_g1_gme, w_g2_gme, w_g3_gme, &
                w_cl_gme)
  ELSE
   IF (.NOT. lcomp_bound) THEN
    DEALLOCATE (t_so_gme, w_so_gme, freshsnw_gme)
   ELSE
    DEALLOCATE (t_so_gme)
   ENDIF
  ENDIF

  ! close file
  CALL close_file (nufile, yin_form_read, icomm_cart, my_cart_id,          &
                   num_compute, lasync_io, idbg_level, yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, izerror, yzerrmsg, 'close_file')
  ENDIF

  IF ( (lchkin) .AND. (my_cart_id == 0) ) THEN
    CLOSE (nuchkdat, STATUS='KEEP')
  ENDIF

  IF (num_compute > 1) THEN
    ! Distribute the unit of time
    CALL distribute_values(nunit_of_time, 1, 0, imp_integers, icomm_cart,  &
                           izerror)
  ENDIF

  ! Check, whether all data necessary are present
  lzcheck = .TRUE.
  DO i = 1, i3e_gme
    lzcheck = lzcheck .AND. lu (i)
    lzcheck = lzcheck .AND. lv (i)
    lzcheck = lzcheck .AND. lt (i)
    lzcheck = lzcheck .AND. lqv(i)
    lzcheck = lzcheck .AND. lqc(i)
    lzcheck = lzcheck .AND. lqi(i)
    lzcheck = lzcheck .AND. lqr(i)
    lzcheck = lzcheck .AND. lqs(i)
  ENDDO
  lzcheck = lzcheck .AND. lfic
  lzcheck = lzcheck .AND. lps
  lzcheck = lzcheck .AND. lt_snow
  lzcheck = lzcheck .AND. lw_snow
  lzcheck = lzcheck .AND. lqv_s
  IF (.NOT. lcomp_bound) THEN
    lzcheck = lzcheck .AND. lw_i
  ENDIF

  IF (.NOT. lmulti_layer_in) THEN
    lzcheck = lzcheck .AND. lt_s
    lzcheck = lzcheck .AND. lt_m
    lzcheck = lzcheck .AND. lw_g1
    lzcheck = lzcheck .AND. lw_g2
    IF (nl_soil_lm == 3) THEN
      lzcheck = lzcheck .AND. lw_g3
    ENDIF
    IF (.NOT. lcomp_bound) THEN
      lzcheck = lzcheck .AND. lt_cl
      lzcheck = lzcheck .AND. lw_cl
    ENDIF
  ELSE
    IF (.NOT. lcomp_bound) THEN
      lzcheck = lzcheck .AND. lfreshsnw
      IF (lprog_rho_snow) THEN
        lzcheck = lzcheck .AND. lrho_snow
      ENDIF
      lzcheck = lzcheck .AND. lt_so(0)
      DO i = 1, ke_soil_in+1
        lzcheck = lzcheck .AND. lt_so(i)
        lzcheck = lzcheck .AND. lw_so(i)
      ENDDO
    ELSE
      lzcheck = lzcheck .AND. lt_so(0)
    ENDIF
  ENDIF

  IF (lseaice .AND. (.NOT. lcomp_bound)) THEN
    lzcheck = lzcheck .AND. lt_ice
    lzcheck = lzcheck .AND. lh_ice
  ENDIF

  IF (.NOT. lzcheck) THEN
    ! Write check information to output file
    IF (my_cart_id == 0) THEN
      WRITE (*      , '(A)') '  '
      WRITE (*      , '(A)') '  '
      WRITE (*      , '(A)') '    Not all data necessary could be read!!!'
      WRITE (*      , '(A)') '  '
      IF (lprog_qi .AND. lprog_qr_qs) THEN
        WRITE (*      , '(A)')                                             &
              '        Level     U       V       T      QV      QC      QI      QR      QS'
        DO i = 1, i3e_gme
          WRITE (*      , '(I11,8L8)')                                     &
                             i, lu(i), lv(i), lt(i), lqv(i), lqc(i), lqi(i), lqr(i), lqs(i)
        ENDDO
      ELSEIF (lprog_qi .AND. .NOT. lprog_qr_qs) THEN
        WRITE (*      , '(A)')                                             &
              '        Level     U       V       T      QV      QC      QI'
        DO i = 1, i3e_gme
          WRITE (*      , '(I11,6L8)')                                     &
                             i, lu(i), lv(i), lt(i), lqv(i), lqc(i), lqi(i)
        ENDDO
      ELSEIF (.NOT. lprog_qi .AND. lprog_qr_qs) THEN
        WRITE (*      , '(A)')                                             &
              '        Level     U       V       T      QV      QC      QR      QS'
        DO i = 1, i3e_gme
          WRITE (*      , '(I11,7L8)')                                     &
                             i, lu(i), lv(i), lt(i), lqv(i), lqc(i), lqr(i), lqs(i)
        ENDDO
      ELSE
        WRITE (*      , '(A)')                                             &
              '        Level     U       V       T      QV      QC'
        DO i = 1, i3e_gme
          WRITE (*      , '(I11,5L8)') i, lu(i), lv(i), lt(i), lqv(i), lqc(i)
        ENDDO
      ENDIF
      WRITE (*      , '(A,L4)') '    PS:      ', lps
      WRITE (*      , '(A,L4)') '    FIC:     ', lfic
      WRITE (*      , '(A,L4)') '    T_SNOW:  ', lt_snow
      WRITE (*      , '(A,L4)') '    W_SNOW:  ', lw_snow
      WRITE (*      , '(A,L4)') '    W_I:     ', lw_i
      WRITE (*      , '(A,L4)') '    QV_S:    ', lqv_s

    IF (.NOT. lmulti_layer_in) THEN
      WRITE (*      , '(A,L4)') '    T_S:     ', lt_s
      WRITE (*      , '(A,L4)') '    T_M:     ', lt_m
      WRITE (*      , '(A,L4)') '    W_G1:    ', lw_g1
      WRITE (*      , '(A,L4)') '    W_G2:    ', lw_g2
      IF (nl_soil_lm == 3) THEN
        WRITE (*      , '(A,L4)') '    W_G3:    ', lw_g3
      ENDIF
      IF (.NOT. lcomp_bound) THEN
        WRITE (*      , '(A,L4)') '    W_CL:    ', lw_cl
        WRITE (*      , '(A,L4)') '    T_CL:    ', lt_cl
      ENDIF
    ELSE
      IF (.NOT. lcomp_bound) THEN
        WRITE (*      , '(A)') '        Level    T_SO    W_SO'
        DO i = 0, ke_soil_in+1
          IF (i==0) THEN
            WRITE (*      , '(I11,L8)')  i, lt_so(i)
          ELSE
            WRITE (*      , '(I11,2L8)') i, lt_so(i), lw_so(i)
          ENDIF
        ENDDO
        WRITE (*      , '(A,L4)') '    FRESHSNW:', lfreshsnw
        IF (lprog_rho_snow) THEN
          WRITE (*      , '(A,L4)') '    RHO_SNOW:', lrho_snow
        ENDIF
      ELSE
        WRITE (*      , '(I11,L8)')  0, lt_so(0)
      ENDIF
    ENDIF

      IF (lseaice .AND. (.NOT. lcomp_bound)) THEN
        WRITE (*      , '(A,L4)') '    T_ICE:   ', lt_ice
        WRITE (*      , '(A,L4)') '    H_ICE:   ', lh_ice
      ENDIF

      ! Abort program
      yzerrmsg  = ' *** ERROR:  Not all data available ***'
      CALL model_abort (my_cart_id, 5051, yzerrmsg, 'org_gme_interpol')
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE org_gme_interpol

!==============================================================================
!+ Reads and distributes GME grib-records to the other PEs
!------------------------------------------------------------------------------

SUBROUTINE read_gme_grib (nufile, mlocpe, mlevpe, lreqpe, leof,             &
                          ublock, vblock, nlfd, nuvblock, maxuv, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  GME-records are read and distributed to the available PEs, until every
!  PE has got a record for immediate processing. A record that requires special
!  action is also distributed after reading but kept in memory of the 
!  PE that gets this record for processing afterwards. There are the cases:
!    1) u and v are just stored (as grib-blocks):               iaction > 0
!    2) t_cl, t_snow, t_s, t_m, qv_s, w_gx are distributed, 
!          unpacked and just stored without interpolation:      iaction < 0
!    3) ps and t (lowest level) are distributed, unpacked,
!          interpolated and stored:                             iaction < 0
!    4) all other variables are distributed, unpacked,
!          and interpolated:                                    iaction = 0
!
!  The difference between 2) and 3) is handled in the interpolation routine.
!  
! Method:
!
!------------------------------------------------------------------------------

! Arguments
INTEGER (KIND=iintegers), INTENT(IN)     ::  &
  nufile      ! unit number of grib-file

LOGICAL, INTENT (OUT)      ::  &
  lreqpe(0:num_compute-1), & ! indicates whether a record is required
  leof                 ! indicates the end of file

INTEGER (KIND=iintegers), INTENT(OUT)    ::  &
  mlocpe (0:num_compute-1),&! location of records of every PE in variable table
  mlevpe (0:num_compute-1),&! level of a multi-level field
  ierror       ! error status variable

INTEGER (KIND=iintegers), INTENT(IN)     ::  &
  nlfd, nuvblock         ! dimensions of ublock and vblock

INTEGER (KIND=iintegers), INTENT(INOUT)  ::  &
  maxuv                  ! maximal number of levels this PE has got

INTEGER  (KIND=intgribf), INTENT(INOUT)   ::  &
  ublock (nlfd,nuvblock), vblock(nlfd,nuvblock)   ! for storing u and v

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER  (KIND=intgribf)   ::  &
  izbyta, izdimpds, izdimgds,  & ! for unpacking GDS and PDS
  ireadblock (nlfd),           & ! another field for reading GRIB-data
  ilenf, ierrf                   ! error status

INTEGER (KIND=int_ga)      ::  &
  ilen_ga

INTEGER  (KIND=intgribc)   ::  &
  nufilec,              & ! unit number of grib-file to pass to cuegin
  maxlenc,              & ! Maximum length of grib record
  ierrc,                & ! Error return code
  igriblenc               ! actual length of grib record read

INTEGER  (KIND=iintegers)  ::  &
  izerror,              & ! error status variable
  izednr,               & ! grib edition number
  izaction,             & ! characterizes special actions to be performed
  izloc,                & ! location of required record in GME variable table
  izlev,                & ! level of a multi-level field
  iz_info(6),           & ! characteristics of read records
  nzpr, nzprgot,        & !
  izstatus(10),         & ! for status in MPI calls
  izrest,               & ! 
  nzpe,                 & ! ID of PE that gets a special u- or v-level
  nzlevel,              & !
  maxlen,               & ! Maximum length of grib record
  igriblen                ! actual length of grib record read

LOGICAL                    ::  &
  lzreq                   ! to check, whether find_loop has to be exited

CHARACTER (LEN= 25)        ::  &
  yzroutine               ! name of this routine for error handling

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! initialize level counter and logical flags
  leof       = .FALSE.
  lreqpe(:)  = .FALSE.
  ierror     = 0_iintegers
  izloc      = 0_iintegers
  izlev      = 0_iintegers
  izaction   = 0_iintegers
  izbyta     = 9_intgribf
  ierrf      = 0_intgribf
  ierrc      = 0_intgribc
  izerror    = 0_iintegers
  yzroutine  = 'read_gme_grib'

!------------------------------------------------------------------------------
! Section 2: Actions for PE 0
!------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
 
    ! loop over all PEs for distributing the records. This loop has its 
    ! counterpart in the case my_cart_id /= 0, so that all information for
    ! a read record can be broadcasted to all PEs.

    pe_loop_0: DO nzpr = 0, num_compute-1

      ! The special_action_loop is an endless loop that is exited if a record
      ! is found that does not need a special action. Records requiring a 
      ! special action (e.g. u_gme or v_gme) are processed without increasing
      ! the loop counter nzpr of pe_loop.
      ! special_action_loop has its counterpart in the case my_cart_id /= 0

      special_action_loop_0: DO

        !----------------------------------------------------------------------
        ! 2.1: Read a record
        !----------------------------------------------------------------------

        lzreq   = .FALSE.

        ! find_loop is an endless loop that is exited, if a required record
        ! is read from the file. After exiting find_loop, it is determined
        ! whether the record just read needs a special action.

        find_loop:  DO

          ! Read in one GRIB field and get some metadata
          maxlen  = INT (idims_in(8)*iwlength, iintegers)

          SELECT CASE (yin_form_read)

          CASE ('grb1')
#ifdef GRIBDWD
!           IF (izdebug >= 20) THEN
!             PRINT *, '      Calling grib 1 reader'
!           ENDIF

            IF (lasync_io .OR. (num_compute > 1) ) THEN
              CALL mpe_io_read(nufile, ireadblock, igriblen, maxlen, izerror)
              IF ( izerror /= 0 ) THEN
                CALL model_abort (my_cart_id, izerror, 'ERROR in mpe_io_read', &
                                  yzroutine)
              ENDIF
            ELSE
              maxlenc = INT (idims_in(8)*iwlength, intgribc)
              nufilec = INT (nufile, intgribc)
              CALL cuegin (nufilec, maxlenc, ireadblock, igriblenc, ierrc)
              igriblen = INT (igriblenc, iintegers)
              izerror  = INT (ierrc, iintegers)
              IF ( izerror /= 0 ) THEN
                CALL model_abort (my_cart_id, izerror, 'ERROR in cuegin', &
                                  yzroutine)
              ENDIF
            ENDIF
#endif
          CASE ('apix')
#ifdef GRIBAPI

!           IF (izdebug >= 20) THEN
!             PRINT *, '      Calling grib 2 reader'
!           ENDIF

            ilen_ga = maxlen   ! ilenf is intent(inout)
            CALL grib_read_from_file (nufile, ireadblock, ilen_ga, ierrf)

            ! if we need the value of igriblen > HUGE(iintegers), then we are in trouble
            IF (ilen_ga < HUGE(1_iintegers)) THEN
              igriblen = INT (ilen_ga, iintegers)
            ELSE
              igriblen = 1_iintegers
            ENDIF
            IF ( ierrf /= GRIB_SUCCESS ) THEN
              IF (ierrf /= GRIB_END_OF_FILE) THEN
                izerror  = 5
                RETURN
              ELSE
                igriblen  = 0_iintegers
              ENDIF
            ENDIF
#endif

          END SELECT

          IF (igriblen > 0) THEN ! otherwise end of file reached

            IF (yin_form_read == 'grb1') THEN
#ifdef GRIBDWD
              ! Get gds and pds before further processing
              izednr = 1
              izbyta = 9_intgribf
              CALL getpd1 (idwdednr, izbyta, undefgrib, idims_in(8), idims_in(1), &
                           ireadblock, ipds, izdimpds, ierrf)
              izbyta = izbyta + ipds(1)
              CALL getgd1 (idwdednr, izbyta, undefgrib, idims_in(8), idims_in(2), &
                           ireadblock, igds_in, izdimgds, ierrf)

              ! check, whether this variable is required
              CALL check_required_dwd (ipds, izdimpds, igds_in, izdimgds,         &
                                       lzreq, izaction, izloc, izlev)
#endif
            ELSEIF (yin_form_read == 'apix') THEN
#ifdef GRIBAPI
              ! check, whether this variable is required: all meta data are retrieved
              ! within this subroutine
              CALL check_required_api (ireadblock, nlfd, lzreq, izaction,         &
                                       izloc, izlev, ierror)
#endif
            ENDIF

            IF (lzreq) THEN
              EXIT find_loop
            ENDIF

          ELSE    ! igriblen == 0:  also exit find_loop, lzreq remains false
            izaction  = 0
            EXIT find_loop
          ENDIF

        ENDDO find_loop

        ! Now PE 0 has got a record which is required for further processing
        ! or the end of file is reached (igriblen == 0)

        !----------------------------------------------------------------------
        ! 2.2: Distribute information
        !----------------------------------------------------------------------

        ! Distribute the status of record, file and action to all PEs
        IF (num_compute > 1) THEN
          iz_info(1) = igriblen                   ! length of record
          iz_info(2) = nzpr                       ! ID of PE that gets record
          iz_info(3) = izaction                   ! special action
          iz_info(4) = izloc                      ! location of record 
          iz_info(5) = izlev                      ! level of multi-level field
          IF (lzreq) THEN
            iz_info(6) = 1
          ELSE
            iz_info(6) = 0
          ENDIF
          CALL distribute_values (iz_info, 6, 0, imp_integers,        &
                                  icomm_cart, izerror)
        ENDIF

        !----------------------------------------------------------------------
        ! 2.3: Distribute the record
        !----------------------------------------------------------------------

        IF (izaction <= 0) THEN
          ! this is not a record that needs special action, so it is 
          ! just distributed to PE nzpr
          IF (igriblen > 0) THEN
            IF (nzpr == 0) THEN
              ! PE 0 keeps this record for its own
              iblock(:) = ireadblock(:)
            ELSE
              ! send record to PE nzpr
              CALL MPI_SEND (ireadblock, iz_info(1)+4, imp_byte, nzpr, 10001,  &
                             icomm_cart, izerror)
            ENDIF
          ELSE
            EXIT pe_loop_0
          ENDIF

          mlocpe(nzpr) = izloc
          mlevpe(nzpr) = izlev
          lreqpe(nzpr) = lzreq

          EXIT special_action_loop_0
        ELSE
          ! now perform the special actions: first u and v

          ! Compute ID of PE that gets this level izlev
          IF ( (izaction == 1) .OR. (izaction == 2) ) THEN
            izrest = MOD (izlev, num_compute)
            IF (izrest == 0) THEN
              nzpe    = num_compute-1
              nzlevel = izlev / num_compute
            ELSE
              nzpe    = izrest-1
              nzlevel = izlev / num_compute + 1
            ENDIF
 
            IF (nzpe == 0) THEN
              ! PE 0 keeps this record for its own
              IF (izaction == 1) THEN
                ublock(:,nzlevel) = ireadblock(:)
                maxuv = maxuv + 1 ! is only done for u: v is the same
              ELSE
                vblock(:,nzlevel) = ireadblock(:)
              ENDIF
            ELSE
              ! Send this level to PE nzpe
              CALL MPI_SEND (ireadblock, iz_info(1)+4, imp_byte, nzpe, 10002,  &
                             icomm_cart, izerror)
            ENDIF
          ENDIF
        ENDIF
 
      ENDDO special_action_loop_0

    ENDDO pe_loop_0

!------------------------------------------------------------------------------
! Section 3: Actions for all other PEs
!------------------------------------------------------------------------------

  ELSE

    ! loop over all PEs for gathering information about the read records. 
    ! This loop has its counterpart in the case my_cart_id == 0, to distribute
    ! the information about a read record to all PEs.

    pe_loop_x: DO nzpr = 0, num_compute - 1

      ! The special_action_loop is an endless loop that is exited if a record
      ! has been received that does not need a special action. Records 
      ! requiring a special action (e.g. u_gme or v_gme) are processed without 
      ! increasing the loop counter nzpr of pe_loop_0 and pe_loop_x.
      ! special_action_loop_x has its counterpart in the case my_cart_id == 0

      special_action_loop_x: DO

        !----------------------------------------------------------------------
        ! 3.1: There is no part corresponding to Section 2.1 here
        !----------------------------------------------------------------------

        !----------------------------------------------------------------------
        ! 3.2: Get information about read record
        !----------------------------------------------------------------------

        ! Get infos about the status of record, file and action from PE 0
        IF (num_compute > 1) THEN
          CALL distribute_values (iz_info, 6, 0, imp_integers,        &
                                  icomm_cart, izerror)

          igriblen   = iz_info(1)                 ! length of record
          nzprgot    = iz_info(2)                 ! ID of PE that gets record
          izaction   = iz_info(3)                 ! special action
          izloc      = iz_info(4)                 ! location in table
          izlev      = iz_info(5)                 ! level of multi-level field
          IF (iz_info(6) == 0) THEN
            lzreq = .FALSE.
          ELSE
            lzreq = .TRUE.
          ENDIF
        ENDIF

        !----------------------------------------------------------------------
        ! 3.3: Get the record
        !----------------------------------------------------------------------

        IF (izaction <= 0) THEN
          ! this is not a record that needs special action, so it is 
          ! just received from PE 0
          IF (igriblen > 0) THEN
            IF (my_cart_id == nzprgot) THEN
              ! receive record from PE 0
              CALL MPI_RECV (iblock, iz_info(1)+4, imp_byte, 0, 10001,  &
                             icomm_cart, izstatus, izerror)
            ENDIF
          ELSE
            ! end of file is reached: exit pe_loop
            EXIT pe_loop_x
          ENDIF

          mlocpe(nzprgot) = izloc
          mlevpe(nzprgot) = izlev
          lreqpe(nzprgot) = lzreq

          EXIT special_action_loop_x
        ELSE
          ! now perform the special actions: first u and v

          ! Compute ID of PE that gets this level izlev and where to store it
          IF ( (izaction == 1) .OR. (izaction == 2) ) THEN
            izrest = MOD (izlev, num_compute)
            IF (izrest == 0) THEN
              nzpe    = num_compute-1
              nzlevel = izlev / num_compute
            ELSE
              nzpe    = izrest-1
              nzlevel = izlev / num_compute + 1
            ENDIF

            ! receive this level from PE 0 in PE nzpe
            IF (my_cart_id == nzpe) THEN
              IF (izaction == 1) THEN
                ! ublock
                CALL MPI_RECV (ublock(:,nzlevel), iz_info(1)+4, imp_byte,    &
                               0, 10002, icomm_cart, izstatus, izerror)
                maxuv = maxuv + 1 ! is only done for u: v is the same
              ELSE
                ! vblock
                CALL MPI_RECV (vblock(:,nzlevel), iz_info(1)+4, imp_byte,    &
                               0, 10002, icomm_cart, izstatus, izerror)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
 
      ENDDO special_action_loop_x

    ENDDO pe_loop_x

  ENDIF  ! my_cart_id == 0 or /= 0

  IF (igriblen == 0) THEN
    leof = .TRUE.
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE read_gme_grib

!==============================================================================
!+ Checks whether a record is required
!------------------------------------------------------------------------------

SUBROUTINE check_required_dwd (mpds    , idimpds, mgds   , idimgds,         &
                               lrequire, iaction, mloc   , mlev)

!------------------------------------------------------------------------------
!
! Description:
!   A Grib-Record with given product definition section (mpds) and grid 
!   description section (mgds) is checked whether it is required for processing.
!   If a special action is necessary for this record (e.g. u or v), this is
!   indicated in iaction.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist:
INTEGER  (KIND=intgribf), INTENT(IN)   ::  &
  idimpds, idimgds               ! dimensions of PDS and GDS

INTEGER  (KIND=intgribf), INTENT(IN)   ::  &
  mpds(idimpds),               & ! product definition section
  mgds(idimgds)                  ! grid definition section

LOGICAL                 , INTENT(OUT)  ::  &
  lrequire                       ! to indicate a required record

INTEGER (KIND=iintegers), INTENT(OUT)  ::  &
  iaction,                     & ! to indicate a special action
  mloc,                        & ! location in the GME variable table name
  mlev                           ! for mlf: level in the atmosphere

!------------------------------------------------------------------------------

! Local scalars:
INTEGER  (KIND=intgribf)   ::  &
  iztabtyp, izee, izlevtyp, izlevtop, izlevbot

INTEGER (KIND=iintegers)   ::  &
  i, ii, ipos
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  lrequire    = .FALSE.
  iaction     = 0_iintegers
  mloc        = 0_iintegers
  mlev        = 0_iintegers

  iztabtyp    = mpds( 2)
  izlevtyp    = mpds( 8)
  izee        = mpds( 7)
  izlevtop    = mpds( 9)
  izlevbot    = mpds(10)

  IF (.NOT. lcomp_bound) THEN
    ipos = 1
  ELSE
    ipos = 2
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Check if it is a regular grid
!------------------------------------------------------------------------------

  IF (mgds(4) == 0) THEN
    ! this is a variable on a regular grid which is not needed
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Look up in the variable table, which variable it is
!------------------------------------------------------------------------------

  ! This search has to be different for levtyp = 109, 110, 111,112:
  ! for levtyp = 109 (multi level) only geopotential on level kcontrol_fi 
  ! for levtyp = 110 (multi level) levbot and levtop must not be compared
  ! for levtyp = 111,112 it is necessary to compare levbot and levtop

  SELECT CASE (izlevtyp)
  CASE (109)

    mloc = 0
    DO i=1,nvar_in
      IF ((var_in(i)%tabtyp == iztabtyp) .AND.  &
          (var_in(i)%levtyp == izlevtyp) .AND.  &
          (var_in(i)%ee     == izee    ) .AND.  &
          (kcontrol_fi       == izlevbot) .AND.  &
          (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
        mloc = i
        mlev = izlevbot
        EXIT
      ENDIF
    ENDDO
    IF (mloc == 0) THEN
      ! the appropriate level was not found: Return
      RETURN
    ENDIF

  CASE (1, 110)
   
    DO i=1,nvar_in
      IF ((var_in(i)%tabtyp        .EQ. iztabtyp) .AND.  &
          (var_in(i)%levtyp        .EQ. izlevtyp) .AND.  &
          (var_in(i)%ee            .EQ. izee    ) .AND.  &
          (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
        mloc = i
        mlev = izlevtop
        EXIT
      ENDIF
    ENDDO
   
  CASE (111, 112)

    IF ( (iztabtyp==201) .AND. ( (izee==197) .OR. (izee==198) ) ) THEN
      ! these are the multi-dimensional variables for the multi-layer soil
      ! model: these are needed for lmulti_layer_in = .TRUE.
      ! for these variables, levtop has to be compared to msoilgrib_in
      IF (lmulti_layer_in) THEN
        IF (.NOT. lcomp_bound) THEN
          DO i=1,nvar_in
            IF ((var_in(i)%tabtyp        .EQ. iztabtyp) .AND.  &
                (var_in(i)%levtyp        .EQ. izlevtyp) .AND.  &
                (var_in(i)%ee            .EQ. izee    ) .AND.  &
                (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mloc = i
              ! Determine mlev
              DO ii=0,ke_soil_in+1
                IF (msoilgrib_in(ii) == izlevbot) mlev = ii
              ENDDO
              EXIT
            ENDIF
          ENDDO
        ELSE
          IF (izee==197 .AND. izlevbot==0) THEN
            ! only T_SO(0) is needed (as former T_S)
            DO i=1,nvar_in
              IF ((var_in(i)%tabtyp        .EQ. iztabtyp) .AND.  &
                  (var_in(i)%levtyp        .EQ. izlevtyp) .AND.  &
                  (var_in(i)%ee            .EQ. izee    ) .AND.  &
                  (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
                mloc = i
                mlev = 0
                EXIT
              ENDIF
            ENDDO
          ELSE
            ! in this case these variables are not needed
            RETURN
          ENDIF
        ENDIF
      ELSE
        ! in this case these variables are not needed
        RETURN
      ENDIF
    ELSE
      ! these are the variables of the old soil model. These are only needed
      ! if lmulti_layer_in = .FALSE.
      IF (.NOT. lmulti_layer_in) THEN
        DO i=1,nvar_in
          IF ((var_in(i)%tabtyp        .EQ. iztabtyp) .AND.  &
              (var_in(i)%levtyp        .EQ. izlevtyp) .AND.  &
              (var_in(i)%ee            .EQ. izee    ) .AND.  &
              (var_in(i)%levtop        .EQ. izlevtop) .AND.  &
              (var_in(i)%levbot        .EQ. izlevbot) .AND.  &
              (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
            mloc = i
            mlev = 1
            EXIT
          ENDIF
        ENDDO
!US   ELSE
!US     ! in this case only T_S (ee=85,lvtyp=111,levbot=levtop=0) is needed
!US     IF (iztabtyp==2 .AND. izee==85 .AND. izlevtyp==111 .AND.     &
!US                                izlevtop==0 .AND. izlevbot==0) THEN
!US       DO i=1,nvar_in
!US         IF ((var_in(i)%tabtyp        .EQ. iztabtyp) .AND.  &
!US             (var_in(i)%levtyp        .EQ. izlevtyp) .AND.  &
!US             (var_in(i)%ee            .EQ. izee    ) .AND.  &
!US             (var_in(i)%levtop        .EQ. izlevtop) .AND.  &
!US             (var_in(i)%levbot        .EQ. izlevbot) .AND.  &
!US             (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
!US           mloc = i
!US           mlev = 1
!US           EXIT
!US         ENDIF
!US       ENDDO
!US     ELSE
!US       RETURN
!US     ENDIF
      ELSE
        ! for the multi layer soil model these variables are not needed
        RETURN
      ENDIF
    ENDIF

  CASE DEFAULT

    ! this record is not required
    RETURN
  
  END SELECT

  IF (mloc < 9) THEN
    ! This variable is not in the GME variable table (mloc == 0) 
    ! and hence not needed 
    ! OR it is a constant field (1 <= mloc <= 8) which has already been
    ! read or is not needed.
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Set the output variables
!------------------------------------------------------------------------------

  ! If this section is reached, it is a required record and mloc, mlev are
  ! already set
  lrequire = .TRUE.

  ! Determine iaction; this depends on the special variable

  ! ps and t on the lowest level have to be interpolated for themselves
  ! but also have to be stored for further interpolations of other variables.
  ! t_s, t_m, t_snow, qv_s, w_g1, w_g2 have only to be stored until all other
  ! variables are read and interpolated.
  ! t_so and w_so have to be stored for the new soil layer model
  ! iaction is set to a negative value, so that special_action_loop in
  ! read_gme_grib is exited but these variables are stored.

  IF (var_in(mloc)%name == 'PS       ') iaction = -1
  IF ((var_in(mloc)%name== 'T        ') .AND. (mlev == i3e_gme)) iaction = -2

  IF (var_in(mloc)%name == 'T_S      ') iaction = -3
  IF (var_in(mloc)%name == 'T_M      ') iaction = -4
  IF (var_in(mloc)%name == 'T_SNOW   ') iaction = -5
  IF (var_in(mloc)%name == 'QV_S     ') iaction = -6
  IF (var_in(mloc)%name == 'W_G1     ') iaction = -7
  IF (var_in(mloc)%name == 'W_G2     ') iaction = -8
  IF (var_in(mloc)%name == 'W_G3     ') iaction = -9
  IF (var_in(mloc)%name == 'W_CL     ') iaction = -10
  IF (var_in(mloc)%name == 'T_CL     ') iaction = -11
  IF (var_in(mloc)%name == 'W_SO     ') iaction = -10
  IF (var_in(mloc)%name == 'T_SO     ') iaction = -11

  ! u and v have to be stored and processed at the end of all other 
  ! interpolations: set iaction to a positive value, so that
  ! special_action_loop in read_gme_grib is not exited.

  IF (var_in(mloc)%name == 'U        ') iaction = 1
  IF (var_in(mloc)%name == 'V        ') iaction = 2

  ! For all other records iaction is set to 0 (in the initializations)
  ! so that special_action_loop in read_gme_grib is exited.

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE check_required_dwd

!==============================================================================
!==============================================================================
!+ Checks whether a record is required
!------------------------------------------------------------------------------
#ifdef GRIBAPI
SUBROUTINE check_required_api (ireadblock, ilfd, lrequire, iaction,         &
                               mloc, mlev, ierr)

!------------------------------------------------------------------------------
!
! Description:
!   A Grib-Record with given meta data characteristics is checked whether
!   it is required for processing.
!   If a special action is necessary for this record (e.g. u or v), this is
!   indicated in iaction.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist:

INTEGER (KIND=iintegers), INTENT(IN)    ::  &
  ilfd                           ! Dimension of ireadblock

INTEGER  (KIND=intgribf), INTENT(IN)    ::  &
  ireadblock (ilfd)              ! packed GRIB record

LOGICAL                 , INTENT(OUT)   ::  &
  lrequire                       ! to indicate a required record

INTEGER (KIND=iintegers), INTENT(OUT)   ::  &
  iaction,                     & ! to indicate a special action
  mloc,                        & ! location in the GME variable table name
  mlev                           ! for mlf: level in the atmosphere

INTEGER (KIND=iintegers), INTENT(INOUT) ::  &
  ierr                           ! Error code

!------------------------------------------------------------------------------

! Local scalars:
INTEGER  (KIND=intgribf)   ::  &
  ilevel, itoplevel, ibotlevel, iscalval1, iscalval2, iscalfac1, iscalfac2,  &
  izednr

REAL (KIND=ireals)         ::  &
  rfact, rdepth1, rdepth2

INTEGER (KIND=iintegers)   ::  &
  i, ii, ipos, isearchlev, igribid, ireturn

CHARACTER (LEN= 30)        ::  &
  ytypeoflevel            ! typeOfLevel

CHARACTER (LEN= 21)        ::  &
  yshortname,           & ! short name from grib_api
  ygridtyp                ! grid type

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  lrequire    = .FALSE.
  iaction     = 0_iintegers
  mloc        = 0_iintegers
  mlev        = 0_iintegers
  iscalval1   = -1
  iscalval2   = -1
  iscalfac1   = -1
  iscalfac2   = -1

  IF (.NOT. lcomp_bound) THEN
    ipos = 1
  ELSE
    ipos = 2
  ENDIF

  ! Build the grib handle
  CALL grib_new_from_message (igribid, ireadblock, ireturn)
  IF (ireturn /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_api grib_new_from_message  ', ireturn
    ierr = 6
  ENDIF

  ! get some metadata
  CALL grib_get (igribid, 'shortName',              yshortname,   ireturn)
  CALL grib_get (igribid, 'typeOfLevel',            ytypeoflevel, ireturn)
  CALL grib_get (igribid, 'typeOfGrid',             ygridtyp,     ireturn)
  CALL grib_get (igribid, 'editionNumber',          izednr,       ireturn)

  CALL grib_get (igribid, 'level',                  ilevel,       ireturn)
  CALL grib_get (igribid, 'topLevel',               itoplevel,    ireturn)
  CALL grib_get (igribid, 'bottomLevel',            ibotlevel,    ireturn)

!------------------------------------------------------------------------------
! Section 2: Check if it is a regular grid
!------------------------------------------------------------------------------

  IF (ygridtyp /= 'triangular_grid') THEN
    ! this is a variable on a regular grid which is not needed
    CALL grib_release (igribid)
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Look up in the variable table, which variable it is
!------------------------------------------------------------------------------

  ! This search has to be different for different leveltypes and sometimes
  ! also for GRIB1 and GRIB2

  ! for level type 'hybrid'          only geopotential on level kcontrol_fi
  ! for level type 'hybridLayer'     levbot and levtop must not be compared
  ! for level type 'depthBelowLand'  compare the level
  ! for level type 'depthBelowLandLayer'  compare top and bottom level

  SELECT CASE (ytypeoflevel)
  CASE ('hybrid')       ! same for GRIB1 and GRIB2

    ! only control geopotential needed as hybrid
    mloc = 0
    DO i=1,nvar_in
      IF ( (TRIM(var_in(i)%name)    == TRIM(yshortname)  ) .AND.    &
           (TRIM(var_in(i)%ylevtyp) == TRIM(ytypeoflevel)) .AND.    &
           (kcontrol_fi             == ibotlevel)               )  THEN
        mloc = i
        mlev = ibotlevel
        EXIT
      ENDIF
    ENDDO
    IF (mloc == 0) THEN
      ! the appropriate level was not found: Return
      CALL grib_release (igribid)
      RETURN
    ENDIF

  CASE ('surface', 'hybridLayer')       ! same for GRIB1 and GRIB2

    DO i=1,nvar_in
      IF ( (TRIM(var_in(i)%name)    == TRIM(yshortname ))  .AND.    &
           (TRIM(var_in(i)%ylevtyp) == TRIM(ytypeoflevel)) .AND.    &
           (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
        mloc = i
        mlev = itoplevel
        EXIT
      ENDIF
    ENDDO

  CASE ('depthBelowLand')
    ! Variables: T_S, T_M, T_CL, T_SO    from GRIB1 and GRIB2
    !            W_SO                    from GRIB1

    IF     (izednr == 1) THEN

      ! In GRIB1, T_SO and W_SO are both coded with leveltype 111 (depthBelowLand);
      ! Normally, W_SO should be coded with leveltype 112 (depthBelowLandLayer), 
      ! but the soil layers cannot properly be coded then
      ! The levels are coded in cm (integer values) in 'bottomlevel' (ipds(10))

      isearchlev = ibotlevel

      IF ((TRIM(yshortname) == 'T_SO') .OR. (TRIM(yshortname) == 'W_SO')) THEN
        ! these are the multi-dimensional variables for the multi-layer soil
        ! model: these are needed for lmulti_layer_in = .TRUE.
        ! for these variables, ibotlevel has to be compared to msoilgrib_in
        IF (lmulti_layer_in) THEN
          IF (.NOT. lcomp_bound) THEN
            DO i=1,nvar_in
              IF ( (TRIM(var_in(i)%name)    == TRIM(yshortname ))  .AND.    &
                   (TRIM(var_in(i)%ylevtyp) == TRIM(ytypeoflevel)) .AND.    &
                   (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
                mloc = i
                ! Determine mlev
                DO ii=0,ke_soil_lm+1
                  IF (msoilgrib_in(ii) == isearchlev) mlev = ii
                ENDDO
                EXIT
              ENDIF
            ENDDO
          ELSE
            IF ((TRIM(yshortname) == 'T_SO') .AND. (isearchlev == 0)) THEN
              ! only T_SO(0) is needed (as former T_S)
              DO i=1,nvar_in
                IF ( (TRIM(var_in(i)%name)    == TRIM(yshortname ))  .AND.    &
                     (TRIM(var_in(i)%ylevtyp) == TRIM(ytypeoflevel)) .AND.    &
                     (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
                  mloc = i
                  mlev = 0
                  EXIT
                ENDIF
              ENDDO
            ELSE
              ! in this case these variables are not needed
              CALL grib_release (igribid)
              RETURN
            ENDIF
          ENDIF
        ELSE
          ! in this case these variables are not needed
          CALL grib_release (igribid)
          RETURN
        ENDIF
      ELSE
        ! these are the variables of the old soil model. These are only needed
        ! if lmulti_layer_in = .FALSE.
  
        ! We assume that this will be Grib1 only!! Grib2 will not work in this way
        IF (.NOT. lmulti_layer_in) THEN
          DO i=1,nvar_in
             IF ( (TRIM(var_in(i)%name)    == TRIM(yshortname ))  .AND.    &
                  (TRIM(var_in(i)%ylevtyp) == TRIM(ytypeoflevel)) .AND.    &
                  (var_in(i)%levtop        == itoplevel)      .AND.  &
                  (var_in(i)%levbot        == ibotlevel)      .AND.  &
                  (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mloc = i
              mlev = 1
              EXIT
            ENDIF
          ENDDO
        ELSE
          ! for the multi layer soil model these variables are not needed
          CALL grib_release (igribid)
          RETURN
        ENDIF
      ENDIF

    ELSEIF (izednr == 2) THEN
      ! should only be T_SO

      IF (lmulti_layer_in) THEN
        ! in the GRIB1 to GRIB2 conversion, the first level of T_SO(0) is converted to T_S!!
        ! This has to be reset
        IF (TRIM(yshortname) == 'T_S') yshortname = 'T_SO'
      ENDIF

      ! in grib2, depth of soil layers are in meter
      ! to get the proper value in cm, scaledValueOfFirstFixedSurface is needed
      CALL grib_get (igribid, 'scaleFactorOfFirstFixedSurface',  iscalfac1,    ireturn)
      CALL grib_get (igribid, 'scaledValueOfFirstFixedSurface',  iscalval1,    ireturn)
      SELECT CASE (iscalfac1)
      CASE (0)
        rfact  =    1.0_ireals
      CASE (1)
        rfact  =    0.1_ireals
      CASE (2)
        rfact  =   0.01_ireals
      CASE (3)
        rfact  =  0.001_ireals
      CASE DEFAULT
        PRINT *,   ' *** Scaled factor for soil layers not implemented:  ', iscalfac1
        ierr = 7
        CALL grib_release (igribid)
        RETURN
      END SELECT
      rdepth1 = iscalval1 * rfact

      IF (lmulti_layer_in) THEN
        IF (.NOT. lcomp_bound) THEN
          DO i=1,nvar_in
            IF ( (TRIM(var_in(i)%name)    == TRIM(yshortname ))  .AND.    &
                 (TRIM(var_in(i)%ylevtyp) == TRIM(ytypeoflevel)) .AND.    &
                 (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mloc = i
              ! Determine mlev
              DO ii=1,ke_soil_in+1
                IF ( ABS(czmls_in(ii) - rdepth1) < 1.0E-5_ireals ) mlev = ii
              ENDDO
              EXIT
            ENDIF
          ENDDO
        ELSE
          IF ((TRIM(yshortname) == 'T_SO') .AND. (iscalval1 == 0)) THEN
            ! only T_SO(0) is needed (as former T_S)
            DO i=1,nvar_in
              IF ( (TRIM(var_in(i)%name)    == TRIM(yshortname ))  .AND.    &
                   (TRIM(var_in(i)%ylevtyp) == TRIM(ytypeoflevel)) .AND.    &
                   (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
                mloc = i
                mlev = 0
                EXIT
              ENDIF
            ENDDO
          ELSE
            ! in this case these variables are not needed
            CALL grib_release (igribid)
            RETURN
          ENDIF
        ENDIF
      ELSE
        ! in this case these variables are not needed
        CALL grib_release (igribid)
        RETURN
      ENDIF

    ENDIF ! izednr

  CASE ('depthBelowLandLayer')

    ! Variables: W_G1, W_G2, W_G3, W_CL (old soil model)  from GRIB1 (not available in GRIB2)
    !            W_SO                   (new soil model)  from GRIB2

    IF     (izednr == 1) THEN

      isearchlev = ibotlevel

      IF ( TRIM(yshortname) == 'W_SO' ) THEN
        ! these are the multi-dimensional variables for the multi-layer soil
        ! model: these are needed for lmulti_layer_in = .TRUE.
        ! for these variables, levtop has to be compared to msoilgrib_in
        IF (lmulti_layer_in) THEN
          DO i=1,nvar_in
            IF ( (TRIM(var_in(i)%name)    == TRIM(yshortname ))  .AND.    &
                 (TRIM(var_in(i)%ylevtyp) == TRIM(ytypeoflevel)) .AND.    &
                 (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mloc = i
              ! Determine mlev
              DO ii=0,ke_soil_in+1
                IF (msoilgrib_in(ii) == isearchlev) mlev = ii
              ENDDO
              EXIT
            ENDIF
          ENDDO
        ELSE
          ! in this case these variables are not needed
          CALL grib_release (igribid)
          RETURN
        ENDIF
      ELSE
        ! these are the variables of the old soil model. These are only needed
        ! if lmulti_layer_in = .FALSE.
  
        ! We assume that this will be Grib1 only!! Grib2 will not work in this way
        IF (.NOT. lmulti_layer_in) THEN
          DO i=1,nvar_in
             IF ( (TRIM(var_in(i)%name)    == TRIM(yshortname ))  .AND.    &
                  (TRIM(var_in(i)%ylevtyp) == TRIM(ytypeoflevel)) .AND.    &
                  (var_in(i)%levtop        == itoplevel)      .AND.  &
                  (var_in(i)%levbot        == ibotlevel)      .AND.  &
                  (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mloc = i
              mlev = 1
              EXIT
            ENDIF
          ENDDO
        ELSE
          ! for the multi layer soil model these variables are not needed
          CALL grib_release (igribid)
          RETURN
        ENDIF
      ENDIF

    ELSEIF (izednr == 2) THEN
      ! should only be W_SO, but when converting old files, W_Gx will also appear

      ! in grib2, depth of soil layers are in meter
      ! to get the proper value in cm, scaledValueOfFirstFixedSurface is needed
      CALL grib_get (igribid, 'scaleFactorOfFirstFixedSurface',  iscalfac1,    ireturn)
      CALL grib_get (igribid, 'scaledValueOfFirstFixedSurface',  iscalval1,    ireturn)
      SELECT CASE (iscalfac1)
      CASE (0)
        rfact =    1.0_ireals
      CASE (1)
        rfact =    0.1_ireals
      CASE (2)
        rfact =   0.01_ireals
      CASE (3)
        rfact =  0.001_ireals
      CASE DEFAULT
        PRINT *,   ' *** Scaled factor for soil layers not implemented:  ', iscalfac1
        ierr = 7
        CALL grib_release (igribid)
        RETURN
      END SELECT
      rdepth1 = iscalval1 * rfact

      ! to get the proper value in cm, scaledValueOfSecondFixedSurface is needed
      CALL grib_get (igribid, 'scaleFactorOfSecondFixedSurface',  iscalfac2,    ireturn)
      CALL grib_get (igribid, 'scaledValueOfSecondFixedSurface',  iscalval2,    ireturn)
      SELECT CASE (iscalfac2)
      CASE (0)
        rfact =    1.0_ireals
      CASE (1)
        rfact =    0.1_ireals
      CASE (2)
        rfact =   0.01_ireals
      CASE (3)
        rfact =  0.001_ireals
      CASE DEFAULT
        PRINT *,   ' *** Scaled factor for soil layers not implemented:  ', iscalfac2
        ierr = 7
        CALL grib_release (igribid)
        RETURN
      END SELECT
      rdepth2 = iscalval2 * rfact

      IF (lmulti_layer_in) THEN
        DO i=1,nvar_in
          IF (TRIM(var_in(i)%name)    == TRIM(yshortname) ) THEN
            IF (var_in(i)%name(1:4) == 'W_SO') THEN
              ! set ylevtyp for W_SO, W_SO_REL, W_SO_ICE to depthBelowLandLayer:
              var_in(i)%ylevtyp = 'depthBelowLandLayer'
            ENDIF
          ENDIF
          IF ( (TRIM(var_in(i)%name)    == TRIM(yshortname ))  .AND.    &
               (TRIM(var_in(i)%ylevtyp) == TRIM(ytypeoflevel)) .AND.    &
               (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
            mloc = i
            ! Determine mlev
            DO ii=1,ke_soil_in+1
              IF ( (ABS(czhls_in(ii-1) - rdepth1) < 1.0E-5_ireals ) .AND.  &
                   (ABS(czhls_in(ii  ) - rdepth2) < 1.0E-5_ireals ) ) mlev = ii
            ENDDO
            EXIT
          ENDIF
        ENDDO
      ELSE
        ! in this case these variables are not needed
        CALL grib_release (igribid)
        RETURN
      ENDIF

    ENDIF ! izednr

  CASE DEFAULT

    ! this record is not required
    CALL grib_release (igribid)
    RETURN

  END SELECT

  IF (mloc < 20) THEN
    ! These variables are constant fields and not needed here
    ! and hence not needed
    ! OR it is a constant field (1 <= mloc <= 8) which has already been
    ! read or is not needed.
    CALL grib_release (igribid)
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Set the output variables
!------------------------------------------------------------------------------

  ! If this section is reached, it is a required record and mloc, mlev are
  ! already set
  lrequire = .TRUE.

  ! Determine iaction; this depends on the special variable

  ! ps and t on the lowest level have to be interpolated for themselves
  ! but also have to be stored for further interpolations of other variables.
  ! t_s, t_m, t_snow, qv_s, w_g1, w_g2 have only to be stored until all other
  ! variables are read and interpolated.
  ! t_so and w_so have to be stored for the new soil layer model
  ! iaction is set to a negative value, so that special_action_loop in
  ! read_gme_grib is exited but these variables are stored.

  IF (var_in(mloc)%name == 'PS       ') iaction = -1
  IF ((var_in(mloc)%name== 'T        ') .AND. (mlev == i3e_gme)) iaction = -2

  IF (var_in(mloc)%name == 'T_S      ') iaction = -3
  IF (var_in(mloc)%name == 'T_M      ') iaction = -4
  IF (var_in(mloc)%name == 'T_SNOW   ') iaction = -5
  IF (var_in(mloc)%name == 'QV_S     ') iaction = -6
  IF (var_in(mloc)%name == 'W_G1     ') iaction = -7
  IF (var_in(mloc)%name == 'W_G2     ') iaction = -8
  IF (var_in(mloc)%name == 'W_G3     ') iaction = -9
  IF (var_in(mloc)%name == 'W_CL     ') iaction = -10
  IF (var_in(mloc)%name == 'T_CL     ') iaction = -11
  IF (var_in(mloc)%name == 'W_SO     ') iaction = -10
  IF (var_in(mloc)%name == 'T_SO     ') iaction = -11

  ! u and v have to be stored and processed at the end of all other
  ! interpolations: set iaction to a positive value, so that
  ! special_action_loop in read_gme_grib is not exited.

  IF (var_in(mloc)%name == 'U        ') iaction = 1
  IF (var_in(mloc)%name == 'V        ') iaction = 2

  ! For all other records iaction is set to 0 (in the initializations)
  ! so that special_action_loop in read_gme_grib is exited.

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

  CALL grib_release (igribid)

END SUBROUTINE check_required_api
#endif
!==============================================================================
!+ Extracts the vertical coordinate parameters for GME
!------------------------------------------------------------------------------

SUBROUTINE get_gme_vert (igribid)

!------------------------------------------------------------------------------
!
! Description:
!   PE 0 extracts the vertical coordinate parameters for GME from the
!   grid description section of a grib-record and distributes them to all
!   PEs.
!
! Method:
!
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(IN) :: igribid  ! Grib 2 handle

! Local Variables
INTEGER  (KIND=iintegers)  ::  &
  izstat, izerror, k, k1, k2

REAL (KIND=irealgrib)      :: refstf
REAL (KIND=ireals)         :: pv(2*i3e_gme+2)

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine get_gme_vert
!------------------------------------------------------------------------------

  yzroutine = 'get_gme_vert'

  ! Allocate fields for vertical coordinate parameters
  ALLOCATE (ak_in (i3e_gme+1), bk_in (i3e_gme+1),                  &
            akh_in(i3e_gme  ), bkh_in(i3e_gme  ),                  &
            dak_in(i3e_gme  ), dbk_in(i3e_gme  ), STAT = izstat)
 
  ! PE 0 extracts the vertical coordinate parameters
  IF (my_cart_id == 0) THEN
    IF (yin_form_read == 'grb1') THEN
#ifdef GRIBDWD
      k1 = 25
      k2 = 25 + i3e_gme + 1
      DO k = 1, i3e_gme+1
        ak_in(k) = REAL (refstf(igds_in(k1 + k)), ireals)
        bk_in(k) = REAL (refstf(igds_in(k2 + k)), ireals)
      ENDDO
#endif
    ELSEIF (yin_form_read == 'apix') THEN
#ifdef GRIBAPI
      CALL grib_get (igribid, 'pv', pv, izerror)
      DO k = 1, i3e_gme+1
        ak_in(k) = pv(            + k)
        bk_in(k) = pv(i3e_gme + 1 + k)
      ENDDO
#endif
    ENDIF
  ENDIF

  IF (num_compute > 1) THEN
    CALL distribute_values(ak_in, i3e_gme+1, 0, imp_reals, icomm_cart, izerror)
    CALL distribute_values(bk_in, i3e_gme+1, 0, imp_reals, icomm_cart, izerror)
  ENDIF

  ! Every PE computes the other parameters
  DO k = 1, i3e_gme
    akh_in(k) = (ak_in(k) + ak_in(k+1)) * 0.5_ireals
    bkh_in(k) = (bk_in(k) + bk_in(k+1)) * 0.5_ireals
    dak_in(k) =  ak_in(k+1) - ak_in(k)
    dbk_in(k) =  bk_in(k+1) - bk_in(k)
  ENDDO

  ! PE 0 writes the parameters to OUTPUT
  IF (my_cart_id == 0) THEN
    WRITE (noutput,'(A)') '     '
    WRITE (noutput,'(A)') '     Vertical coordinate parameters of GME:'
    WRITE (noutput,'(A)') '     '
    WRITE (noutput,'(A)') '             k        a(k) (Pa)             b(k)'
    DO k = 1,i3e_gme+1
      WRITE (noutput,'(I14,2F17.4)') k, ak_in(k), bk_in(k)
    ENDDO
    WRITE (noutput,'(A)') '     '
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE get_gme_vert

!==============================================================================
!+ Interpolates a GME-field to a LM-subdomain
!------------------------------------------------------------------------------

SUBROUTINE interpol_gme (field_gme, mloc_gme, mlev, gme_loc1, gme_loc2,    &
                         ps_gme, ti3e_gme, t_s_gme, t_m_gme, t_snow_gme,   &
                         qv_s_gme, w_g1_gme, w_g2_gme, w_g3_gme, w_cl_gme, &
                         t_cl_gme, t_so_gme, w_so_gme, freshsnw_gme,       &
                         rho_snow_gme,                                     &
                         lu, lv, lt, lqc, lqv, lqi, lqr, lqs, lps, lt_cl,  &
                         lt_s, lt_snow, lt_m, lqv_s, lw_cl, lw_i, lw_snow, &
                         lw_g1, lw_g2, lw_g3, ldpsdt, lfic, lt_so, lw_so,  &
                         lfreshsnw, lrho_snow, lt_ice, lh_ice)

!------------------------------------------------------------------------------
!
! Description:
!   The global gme field is available on every PE, so that it can be 
!   interpolated to the LM subdomain. With mloc_gme (location in the GME 
!   variable table) the interpolation code is checked and the field is
!   interpolated using the interpolation utilities. The result is written 
!   to the corresponding LM-field.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
REAL     (KIND=ireals)   ,  INTENT(IN)     ::  &
  field_gme (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max)
  ! record to be interpolated

REAL     (KIND=ireals)   ,  INTENT(OUT)    ::  &
  gme_loc1  (igg1sm1:igg1ep1, igg2sm1:igg2ep1, jd_min:jd_max),   & !
  gme_loc2  (igg1sm1:igg1ep1, igg2sm1:igg2ep1, jd_min:jd_max)
  ! these fields are only needed as working arrays for pp_interp2qs

INTEGER  (KIND=iintegers),  INTENT(IN)     ::  &
  mloc_gme,             & ! location of variable in GME-variable table
  mlev                    ! level of a multi-level field

LOGICAL                  ,  INTENT(INOUT)  ::  &
  lu(i3e_gme), lv(i3e_gme), lt(i3e_gme), lqc(i3e_gme), lqv(i3e_gme),      &
  lqi(i3e_gme), lqr(i3e_gme), lqs(i3e_gme),                               &
  lps, lt_cl, lt_s, lt_snow, lt_m, lqv_s, lw_cl, lw_i, lw_snow, lw_g1,    &
  lw_g2, lw_g3, ldpsdt, lfic, lt_so(0:ke_soil_in+1), lw_so(ke_soil_in+1),       &
  lfreshsnw, lrho_snow, lt_ice, lh_ice

! Arrays for intermediate storage of some GME fields
REAL (KIND=ireals)       ,  INTENT(INOUT)  ::  &
  ps_gme     (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  ti3e_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  t_s_gme    (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  t_m_gme    (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  t_snow_gme (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  t_cl_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  qv_s_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  w_g1_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  w_g2_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  w_g3_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  w_cl_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  t_so_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2,0:ke_soil_in+1,jd_min:jd_max), &
  w_so_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2,1:ke_soil_in+1,jd_min:jd_max), &
  freshsnw_gme(igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max),    & !
  rho_snow_gme(igg1sm2:igg1ep2,igg2sm2:igg2ep2,jd_min:jd_max)

!------------------------------------------------------------------------------
!
! Local variables:
INTEGER  (KIND=iintegers)  ::  &
  izerror,              & ! status and error status variable
  mloc_lm,              & ! location of variable in LM variable table
  izlen, i

REAL    (KIND=ireals)      ::  &
  field_lm(ie2lm,je2lm)

LOGICAL                    ::  &
  lzmono, lzposdef

CHARACTER (LEN=  1)        ::  &
  yzitype     ! interpolation type

CHARACTER (LEN= 11)        ::  &
  yzgmename, & ! name of gme field
  yzlmname,  & ! name of lm field
  yzname       ! for testing in the IF-construction

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

  yzlmname   = '           '
  yzgmename  = '           '
  yzname     = '           '
  yzroutine  = 'interpolate'
  yzitype    = var_in(mloc_gme)%ipc(1:1)
  IF (var_in(mloc_gme)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mloc_gme)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF
  izlen      = LEN_TRIM(var_in(mloc_gme)%name)
  yzgmename  = var_in(mloc_gme)%name(1:izlen)

  ! Look for location of variable in LM variable table
  mloc_lm = 0
  DO i = 1, nvar_lm
    yzlmname(1:izlen) = var_lm(i)%name(1:izlen)
    yzname            = yzlmname(1:izlen)
    IF ( (yzlmname(1:izlen) == yzgmename(1:izlen)) .AND.                &
         (yzname            == var_lm(i)%name) ) THEN
      mloc_lm = i
      EXIT
    ENDIF
  ENDDO

  IF ( (mloc_lm == 0) .AND. (yzgmename /= 'FI') ) THEN
    yzerrmsg = 'Variable not found in LM variable table after interpolation'
    CALL model_abort(my_cart_id, 9011, yzerrmsg, yzroutine)
  ENDIF

!------------------------------------------------------------------------------
! Section 2: 
!------------------------------------------------------------------------------

  SELECT CASE (var_in(mloc_gme)%name)
  CASE ('T_S      ')
    t_s_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)

    ! Interpolate t_s_gme to t_s_gl
    CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,         &
                      jd_min, jd_max, ispoke, baryll_m, index_m,             &
                      lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,          &
                      w_intpol, n_intpol, m_intpol, l_intpol,                &
                      undef, t_s_gl, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 9510_iintegers, yzerrmsg, yzroutine)
    ENDIF
    lt_s      = .TRUE.
  CASE ('T_M      ')
    t_m_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    lt_m      = .TRUE.
  CASE ('T_SNOW   ')
    t_snow_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    lt_snow   = .TRUE.
  CASE ('T_CL     ')
    t_cl_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    lt_cl     = .TRUE.
  CASE ('QV_S     ')
    qv_s_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    lqv_s     = .TRUE.
  CASE ('W_G1     ')
    w_g1_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    lw_g1     = .TRUE.
  CASE ('W_G2     ')
    w_g2_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    lw_g2     = .TRUE.
  CASE ('W_G3     ')
    w_g3_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    lw_g3     = .TRUE.
  CASE ('W_CL     ')
    w_cl_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    lw_cl     = .TRUE.
  CASE ('T_SO     ')
    t_so_gme(:,:,mlev,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    lt_so(mlev)  = .TRUE.

    IF (mlev == 0) THEN
      ! Interpolate t_so_gme(mlev=0) to t_s_gl
      CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,         &
                        jd_min, jd_max, ispoke, baryll_m, index_m,             &
                        lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,          &
                        w_intpol, n_intpol, m_intpol, l_intpol,                &
                        undef, t_s_gl, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
      IF (izerror /= 0) THEN
        CALL model_abort (my_cart_id, 9510_iintegers, yzerrmsg, yzroutine)
      ENDIF
    ENDIF
  CASE ('W_SO     ')
    w_so_gme(:,:,mlev,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    lw_so(mlev)  = .TRUE.
  CASE ('FI       ')
    ! Interpolate control level to fic_gl
    CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,         &
                      jd_min, jd_max, ispoke, baryll_m, index_m,             &
                      lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,          &
                      w_intpol, n_intpol, m_intpol, l_intpol,                &
                      undef, fic_gl, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 9510_iintegers, yzerrmsg, yzroutine)
    ENDIF

  CASE DEFAULT
    IF (yzitype == 'Q') THEN
      ! quadratic interpolation to LM-field
      CALL pp_interp2qs(field_gme, igg1sm2 , igg1ep2  , igg2sm2  , igg2ep2 ,  &
                        eta_glob , chi_glob, cpsi_glob, spsi_glob, grd_glob,  &
                        jd_min, jd_max, igg1sm1, igg1ep1, igg2sm1, igg2ep1 ,  &
                        ispoke, ispokes, i1mrp,  i2mrp, baryll_m , index_m ,  &
                        lzmono, lzposdef,                                     &
                        field_lm, gme_loc1, gme_loc2, 1, ie2lm, 1, je2lm,     &
                        j1_min , j1_max  , j2_min , j2_max , jd_min , jd_max ,&
                        isp11  , isp12   , isp13  , isp14, isp15, isp16,      &
                        isp21  , isp22   , isp23  , isp24, isp25, isp26,      &
                        r_earth, .FALSE. , izerror)
      IF (izerror /= 0) THEN
        CALL model_abort (my_cart_id, 9510_iintegers, yzerrmsg, yzroutine)
      ENDIF
    ELSE
      ! normal linear (or match) interpolation to LM-field
      CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                        jd_min, jd_max, ispoke, baryll_m, index_m,            &
                        lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                        w_intpol, n_intpol, m_intpol, l_intpol,               &
                        undef, field_lm, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
      IF (izerror /= 0) THEN
        CALL model_abort (my_cart_id, 9510_iintegers, yzerrmsg, yzroutine)
      ENDIF
    ENDIF

    SELECT CASE (var_lm(mloc_lm)%rank)
    CASE (3)
      IF (ASSOCIATED (var_lm(mloc_lm)%p3)) THEN
        var_lm(mloc_lm)%p3(1:ie2lm,1:je2lm,mlev) = field_lm(1:ie2lm,1:je2lm)
      ELSE
        print *, '*** No memory allocated for ', var_lm(mloc_lm)%name
        yzerrmsg = 'Pointer of variable table is not associated '
        CALL model_abort(my_cart_id, 9011, yzerrmsg, yzroutine)
      ENDIF
    CASE (2)
      IF (ASSOCIATED (var_lm(mloc_lm)%p2)) THEN
        var_lm(mloc_lm)%p2(1:ie2lm,1:je2lm)      = field_lm(1:ie2lm,1:je2lm)
      ELSE
        print *, '*** No memory allocated for ', var_lm(mloc_lm)%name
        yzerrmsg = 'Pointer of variable table is not associated '
        CALL model_abort(my_cart_id, 9011, yzerrmsg, yzroutine)
      ENDIF
    END SELECT

    ! ps and t on the lowest level have also to be stored
    IF (var_in(mloc_gme)%name == 'PS       ') THEN
      ps_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    ELSEIF ((var_in(mloc_gme)%name== 'T        ') .AND. (mlev == i3e_gme)) THEN
      ti3e_gme(:,:,jd_min:jd_max) = field_gme(:,:,jd_min:jd_max)
    ENDIF

  END SELECT ! gme name

! if (my_cart_id == 0) then
!   if (yzgmename /= 'FI') then
!     print *, var_lm(mloc_lm)%name, ' level ', mlev, 'has been interpolated and/or stored'
!   else
!     print *, 'FIC_gl_lm  ', ' level ', mlev, 'has been interpolated and/or stored'
!   endif
! endif

  SELECT CASE (var_in(mloc_gme)%name)
  CASE ('T        ')
    lt (mlev) = .TRUE.
  CASE ('QV       ')
    lqv(mlev) = .TRUE.
  CASE ('QC       ')
    lqc(mlev) = .TRUE.
  CASE ('QI       ')
    lqi(mlev) = .TRUE.
  CASE ('QR       ')
    lqr(mlev) = .TRUE.
  CASE ('QS       ')
    lqs(mlev) = .TRUE.
  CASE ('FI       ')
    lfic      = .TRUE.
  CASE ('PS       ')
    lps       = .TRUE.
  CASE ('W_SNOW   ')
    lw_snow   = .TRUE.
  CASE ('W_I      ')
    lw_i      = .TRUE.
  CASE ('DPSDT    ')
    ldpsdt    = .TRUE.
  CASE ('FRESHSNW ')
    lfreshsnw = .TRUE.
  CASE ('RHO_SNOW ')
    lrho_snow = .TRUE.
  CASE ('T_ICE    ')
    lt_ice    = .TRUE.
  CASE ('H_ICE    ')
    lh_ice    = .TRUE.
  END SELECT

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE interpol_gme

!==============================================================================
!+ Interpolates a u and v from GME to LM-subdomains
!------------------------------------------------------------------------------

SUBROUTINE interpol_gme_uv (u_gme, v_gme, mlocu_gme, mlocv_gme,             &
                            mlocu_lm, mlocv_lm, mlev,                       &
                            gme_loc1, gme_loc2, gme_loc3, gme_loc4, lu, lv)

!------------------------------------------------------------------------------
!
! Description:
!   The global fields u and v from GME are available on every PE, so that they
!   can be interpolated to the LM subdomain. With mloc_gme (location in the GME 
!   variable table) the interpolation code is checked and the fields are
!   interpolated using the interpolation utilities. The result is written 
!   to the corresponding LM-fields.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
REAL     (KIND=ireals)   ,  INTENT(IN)     ::  &
  u_gme (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max), & !
  v_gme (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max)

REAL     (KIND=ireals)   ,  INTENT(OUT)    ::  &
  gme_loc1  (igg1sm1:igg1ep1, igg2sm1:igg2ep1, jd_min:jd_max),   & !
  gme_loc2  (igg1sm1:igg1ep1, igg2sm1:igg2ep1, jd_min:jd_max),   & !
  gme_loc3  (igg1sm1:igg1ep1, igg2sm1:igg2ep1, jd_min:jd_max),   & !
  gme_loc4  (igg1sm1:igg1ep1, igg2sm1:igg2ep1, jd_min:jd_max)
  ! these fields are only needed as working arrays for pp_interp2qv

INTEGER  (KIND=iintegers),  INTENT(IN)     ::  &
  mlocu_gme, mlocv_gme, & ! location of variables in GME-variable table
  mlocu_lm , mlocv_lm , & ! location of variables in LM-variable table
  mlev                    ! level of a multi-level field

LOGICAL                  ,  INTENT(INOUT)  ::  &
  lu(i3e_gme), lv(i3e_gme)

!------------------------------------------------------------------------------
!
! Local variables:
INTEGER  (KIND=iintegers)  ::  &
  izerror                 ! status and error status variable

REAL    (KIND=ireals)      ::  &
  zfield_lm(ie2lm,je2lm)

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  yzroutine  = 'interpolate_uv'

!------------------------------------------------------------------------------
! Section 2: 
!------------------------------------------------------------------------------

  ! quadratic interpolation to LM u-points
  CALL pp_interp2qv(u_gme, v_gme, igg1sm2, igg1ep2, igg2sm2  , igg2ep2 ,  &
                    jd_min, jd_max,                                       &
                    eta_glob , chi_glob, cpsi_glob, spsi_glob, grd_glob,  &
                    igg1sm1 , igg1ep1  , igg2sm1  , igg2ep1 ,             &
                    ispoke, ispokes, i1mrp,  i2mrp, baryll_u , rotang_u,  &
                    index_u,                                              &
                    u_lm(:,:,mlev), zfield_lm, gme_loc1, gme_loc2,        &
                    gme_loc3, gme_loc4, 1, ie2lm, 1, je2lm,               &
                    j1_min , j1_max  , j2_min , j2_max , jd_min , jd_max ,&
                    isp11  , isp12   , isp13  , isp14, isp15, isp16,      &
                    isp21  , isp22   , isp23  , isp24, isp25, isp26,      &
                    r_earth, .FALSE. , izerror)

  CALL uv2uvrot_vec (u_lm(:,:,mlev), zfield_lm(:,:), latlm_u(:,:), lonlm_u(:,:), &
                     pollat, pollon, ie2lm, je2lm)

  ! quadratic interpolation to LM v-points
  CALL pp_interp2qv(u_gme, v_gme, igg1sm2, igg1ep2, igg2sm2  , igg2ep2 ,  &
                    jd_min, jd_max,                                       &
                    eta_glob , chi_glob, cpsi_glob, spsi_glob, grd_glob,  &
                    igg1sm1 , igg1ep1  , igg2sm1  , igg2ep1 ,             &
                    ispoke, ispokes, i1mrp,  i2mrp, baryll_v , rotang_v,  &
                    index_v,                                              &
                    zfield_lm, v_lm(:,:,mlev), gme_loc1, gme_loc2,        &
                    gme_loc3, gme_loc4, 1, ie2lm, 1, je2lm,               &
                    j1_min , j1_max  , j2_min , j2_max , jd_min , jd_max ,&
                    isp11  , isp12   , isp13  , isp14, isp15, isp16,      &
                    isp21  , isp22   , isp23  , isp24, isp25, isp26,      &
                    r_earth, .FALSE. , izerror)

  CALL uv2uvrot_vec (zfield_lm(:,:), v_lm(:,:,mlev), latlm_v(:,:), lonlm_v(:,:), &
                     pollat, pollon, ie2lm, je2lm)

  lu (mlev) = .TRUE.
  lv (mlev) = .TRUE.

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE interpol_gme_uv

!==============================================================================
!+ Interpolates the special GME-fields to the LM-subdomain
!------------------------------------------------------------------------------

SUBROUTINE interpol_gme_special (ps_gme, ti3e_gme, t_s_gme, t_m_gme,         &
                        t_snow_gme, t_cl_gme, qv_s_gme, w_g1_gme, w_g2_gme,  &
                        w_g3_gme, w_cl_gme, t_so_gme, w_so_gme, field_gme,   &
                        lu, lv, lt, lqc, lqv, lqi, lps, lt_cl, lt_s, lt_snow,&
                        lt_m, lqv_s, lw_cl, lw_i, lw_snow, lw_g1, lw_g2,     &
                        lw_g3, lt_so, lw_so, ldpsdt)

!------------------------------------------------------------------------------
!
! Description:
!   The global gme field is available on every PE, so that it can be 
!   interpolated to the LM subdomain. With mloc_gme (location in the GME 
!   variable table) the interpolation code is checked and the field is
!   interpolated using the interpolation utilities. The result is written 
!   to the corresponding LM-field.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
! Arrays for intermediate storage of some GME fields
REAL (KIND=ireals)       ,  INTENT(IN)     ::  &
  ps_gme     (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  ti3e_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  t_s_gme    (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  t_m_gme    (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  t_snow_gme (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  t_cl_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  qv_s_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max)

REAL (KIND=ireals)       ,  INTENT(INOUT)  ::  &
  w_g1_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  w_g2_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  w_g3_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  w_cl_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max),   & !
  t_so_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, 0:ke_soil_in+1, jd_min:jd_max),   & !
  w_so_gme   (igg1sm2:igg1ep2, igg2sm2:igg2ep2, 1:ke_soil_in+1, jd_min:jd_max)      !

REAL     (KIND=ireals)   ,  INTENT(OUT)    ::  &
  field_gme (igg1sm2:igg1ep2, igg2sm2:igg2ep2, jd_min:jd_max)

LOGICAL                  ,  INTENT(INOUT)  ::  &
  lu(i3e_gme), lv(i3e_gme), lt(i3e_gme), lqc(i3e_gme), lqv(i3e_gme),      &
  lqi(i3e_gme),                                                           &
  lps, lt_cl, lt_s, lt_snow, lt_m, lqv_s, lw_cl, lw_i, lw_snow, lw_g1,    &
  lw_g2, lw_g3, ldpsdt, lt_so(0:ke_soil_in+1), lw_so(1:ke_soil_in+1)

!------------------------------------------------------------------------------
!
! Local variables:
INTEGER  (KIND=iintegers)  ::  &
  izerror,              & ! status and error status variable
  i, k, j1, j2, jd, s_t

INTEGER  (KIND=iintegers)  ::  &
  mzts_loc_lm, mztm_loc_lm, mztsnow_loc_lm, mzqvs_loc_lm, mzwg1_loc_lm,      &
  mzwg2_loc_lm, mzwg3_loc_lm, mzwcl_loc_lm, mztcl_loc_lm,                    &
  mztso_loc_lm, mzwso_loc_lm,                                                &
  mzts_loc_gme, mztm_loc_gme, mztsnow_loc_gme, mzqvs_loc_gme, mzwg1_loc_gme, &
  mzwg2_loc_gme, mzwg3_loc_gme, mzwcl_loc_gme, mztcl_loc_gme,                &
  mztso_loc_gme, mzwso_loc_gme

REAL    (KIND=ireals)      ::  &
  zwg1_lm, zwg2_lm , zwg_gme

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
! Section 1: Initializations
!------------------------------------------------------------------------------

  yzroutine  = 'interpolate_special'

  ! Look for locations of special variables in the variable tables
  DO i = 1, nvar_lm
    SELECT CASE (var_lm(i)%name)
    CASE ('T_S_lm     '); mzts_loc_lm     = i
    CASE ('T_M_lm     '); mztm_loc_lm     = i
    CASE ('T_SNOW_lm  '); mztsnow_loc_lm  = i
    CASE ('T_CL_lm    '); mztcl_loc_lm    = i
    CASE ('QV_S_lm    '); mzqvs_loc_lm    = i
    CASE ('W_G1_lm    '); mzwg1_loc_lm    = i
    CASE ('W_G2_lm    '); mzwg2_loc_lm    = i
    CASE ('W_G3_lm    '); mzwg3_loc_lm    = i
    CASE ('W_CL_lm    '); mzwcl_loc_lm    = i
    CASE ('T_SO_lm    '); mztso_loc_lm    = i
    CASE ('W_SO_lm    '); mzwso_loc_lm    = i
    END SELECT
  ENDDO
  DO i = 1, nvar_in
    SELECT CASE (var_in(i)%name)
    CASE ('T_S      '); mzts_loc_gme    = i
    CASE ('T_M      '); mztm_loc_gme    = i
    CASE ('T_SNOW   '); mztsnow_loc_gme = i
    CASE ('T_CL     '); mztcl_loc_gme   = i
    CASE ('QV_S     '); mzqvs_loc_gme   = i
    CASE ('W_G1     '); mzwg1_loc_gme   = i
    CASE ('W_G2     '); mzwg2_loc_gme   = i
    CASE ('W_G3     '); mzwg3_loc_gme   = i
    CASE ('W_CL     '); mzwcl_loc_gme   = i
    CASE ('T_SO     '); mztso_loc_gme   = i
    CASE ('W_SO     '); mzwso_loc_gme   = i
    END SELECT
  ENDDO

!------------------------------------------------------------------------------
! Section 2: Interpolation of t(lowest level) - t_s --> dtkes_gl
!------------------------------------------------------------------------------

  yzitype    = var_in(mzts_loc_gme)%ipc(1:1)
  IF (var_in(mzts_loc_gme)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mzts_loc_gme)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF

  DO jd = jd_min, jd_max
    DO j2 = igg2sm2, igg2ep2
      DO j1 = igg1sm2, igg1ep2
        IF (ti3e_gme(j1,j2,jd) /= undef) THEN
          IF (.NOT. lmulti_layer_in) THEN
            field_gme(j1,j2,jd) = ti3e_gme(j1,j2,jd) - t_s_gme(j1,j2,jd)
          ELSE
            field_gme(j1,j2,jd) = ti3e_gme(j1,j2,jd) - t_so_gme(j1,j2,0,jd)
          ENDIF
        ELSE
          field_gme(j1,j2,jd) = undef
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  ! normal linear (or match) interpolation to LM-field
  CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                    jd_min, jd_max, ispoke, baryll_m, index_m,            &
                    lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                    w_intpol, n_intpol, m_intpol, l_intpol,               &
                    undef, dtkes_gl, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Interpolation of t_m - t_s --> dtms_gl
!------------------------------------------------------------------------------

IF (.NOT. lmulti_layer_in) THEN
  yzitype    = var_in(mztm_loc_gme)%ipc(1:1)
  IF (var_in(mztm_loc_gme)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mztm_loc_gme)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF

  DO jd = jd_min, jd_max
    DO j2 = igg2sm2, igg2ep2
      DO j1 = igg1sm2, igg1ep2
        IF (t_m_gme(j1,j2,jd) /= undef) THEN
          field_gme(j1,j2,jd) = t_m_gme(j1,j2,jd) - t_s_gme(j1,j2,jd)
        ELSE
          field_gme(j1,j2,jd) = undef
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  ! normal linear (or match) interpolation to LM-field
  CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                    jd_min, jd_max, ispoke, baryll_m, index_m,            &
                    lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                    w_intpol, n_intpol, m_intpol, l_intpol,               &
                    undef, dtms_gl, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
  ENDIF
ENDIF

!------------------------------------------------------------------------------
! Section 4: Interpolation of t_s - t_snow --> dtssnow_gl
!------------------------------------------------------------------------------

  yzitype    = var_in(mztsnow_loc_gme)%ipc(1:1)
  IF (var_in(mztsnow_loc_gme)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mztsnow_loc_gme)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF

  DO jd = jd_min, jd_max
    DO j2 = igg2sm2, igg2ep2
      DO j1 = igg1sm2, igg1ep2
        IF (t_snow_gme(j1,j2,jd) /= undef) THEN
          IF (.NOT. lmulti_layer_in) THEN
            field_gme(j1,j2,jd) = t_s_gme(j1,j2,jd) - t_snow_gme(j1,j2,jd)
          ELSE
            field_gme(j1,j2,jd) = t_so_gme(j1,j2,0,jd) - t_snow_gme(j1,j2,jd)
          ENDIF
        ELSE
          field_gme(j1,j2,jd) = undef
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  ! normal linear (or match) interpolation to LM-field
  CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                    jd_min, jd_max, ispoke, baryll_m, index_m,            &
                    lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                    w_intpol, n_intpol, m_intpol, l_intpol,               &
                    undef, dtssnow_gl, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
  ENDIF

!------------------------------------------------------------------------------
! Section 4a: Interpolation of t_cl
!------------------------------------------------------------------------------

IF (.NOT. lmulti_layer_in .AND. (.NOT. lcomp_bound)) THEN
  IF (itype_t_cl == 0) THEN
    ! t_cl_gme is interpolated to t_cl_lm
    yzitype    = var_in(mztcl_loc_gme)%ipc(1:1)
    IF (var_in(mztcl_loc_gme)%ipc(2:2) == 'T') THEN
      lzmono   = .TRUE.
    ELSE
      lzmono   = .FALSE.
    ENDIF
    IF (var_in(mztcl_loc_gme)%ipc(3:3) == 'T') THEN
      lzposdef = .TRUE.
    ELSE
      lzposdef = .FALSE.
    ENDIF

    DO jd = jd_min, jd_max
      DO j2 = igg2sm2, igg2ep2
        DO j1 = igg1sm2, igg1ep2
          IF (t_cl_gme(j1,j2,jd) /= undef) THEN
            field_gme(j1,j2,jd) = t_cl_gme(j1,j2,jd)
          ELSE
            field_gme(j1,j2,jd) = undef
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! normal linear (or match) interpolation to LM-field
    CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                      jd_min, jd_max, ispoke, baryll_m, index_m,            &
                      lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                      w_intpol, n_intpol, m_intpol, l_intpol,               &
                      undef, t_cl_lm,    1, ie2lm, 1, je2lm, yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
    ENDIF

!US  ELSEIF (itype_t_cl = 1) THEN
!US    ! t_cl_lm has been read with external parameters, nothing has to be done
  ENDIF
ENDIF

!------------------------------------------------------------------------------
! Section 5: Interpolation of qv_s --> rh_s_gl
!------------------------------------------------------------------------------

  yzitype    = var_in(mzqvs_loc_gme)%ipc(1:1)
  IF (var_in (mzqvs_loc_gme)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mzqvs_loc_gme)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF

  DO jd = jd_min, jd_max
    DO j2 = igg2sm2, igg2ep2
      DO j1 = igg1sm2, igg1ep2
        IF (.NOT. lmulti_layer_in) THEN
          IF (t_s_gme(j1,j2,jd) /= undef) THEN
            field_gme(j1,j2,jd) = qv_s_gme(j1,j2,jd) /                      &
                 qsat(psat_w(t_s_gme(j1,j2,jd),b1, b2_w, b3, b4_w),   &
                     ps_gme(j1,j2,jd), Rdv, O_m_rdv)
          ELSE
            field_gme(j1,j2,jd) = undef
          ENDIF
        ELSE
          IF (t_so_gme(j1,j2,0,jd) /= undef) THEN
            field_gme(j1,j2,jd) = qv_s_gme(j1,j2,jd) /                      &
                 qsat(psat_w(t_so_gme(j1,j2,0,jd),b1, b2_w, b3, b4_w),   &
                     ps_gme(j1,j2,jd), Rdv, O_m_rdv)
          ELSE
            field_gme(j1,j2,jd) = undef
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  ! normal linear (or match) interpolation to LM-field
  CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                    jd_min, jd_max, ispoke, baryll_m, index_m,            &
                    lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                    w_intpol, n_intpol, m_intpol, l_intpol,               &
                    undef, rh_s_gl, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
  ENDIF

!------------------------------------------------------------------------------
! Section 6: Interpolation of w_gx --> w_gx_lm
!------------------------------------------------------------------------------

IF (.NOT. lmulti_layer_in) THEN
  yzitype    = var_in(mzwg1_loc_gme)%ipc(1:1)
  IF (var_in(mzwg1_loc_gme)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mzwg1_loc_gme)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF

  IF (l_smi) THEN
    ! Scale the soil moisture with the soil moisture index SMI:
    ! smi = (sm-pwp)/(fc-pwp)
    ! Caution: sm is in volumetric units, soil depth in cm!
    DO jd = jd_min, jd_max
      DO j2 = igg2sm2, igg2ep2
        DO j1 = igg1sm2, igg1ep2
          IF ((soiltyp_gme(j1,j2,jd) /= undef) .AND.                       &
                 (w_g1_gme(j1,j2,jd) /= undef)) THEN
            s_t = NINT (soiltyp_gme(j1,j2,jd))
            IF (s_t >= 3 .AND. s_t <= 8) THEN
              w_g1_gme(j1,j2,jd) = w_g1_gme(j1,j2,jd) /                       &
               (var_in(mzwg1_loc_gme)%levbot-var_in(mzwg1_loc_gme)%levtop)*0.01
              w_g1_gme(j1,j2,jd) = MAX ((w_g1_gme(j1,j2,jd) - pwpb(s_t)) / &
                                          (fcb(s_t) - pwpb(s_t)),0.0_ireals)
              w_g2_gme(j1,j2,jd) = w_g2_gme(j1,j2,jd) /                       &
               (var_in(mzwg2_loc_gme)%levbot-var_in(mzwg2_loc_gme)%levtop)*0.01
              w_g2_gme(j1,j2,jd) = MAX ((w_g2_gme(j1,j2,jd) - pwpb(s_t)) / &
                                          (fcb(s_t) - pwpb(s_t)),0.0_ireals)
              IF (nl_soil_in == 3) THEN
                w_g3_gme(j1,j2,jd) = w_g3_gme(j1,j2,jd) /                      &
                (var_in(mzwg3_loc_gme)%levbot-var_in(mzwg3_loc_gme)%levtop)*0.01
                w_g3_gme(j1,j2,jd) = MAX ((w_g3_gme(j1,j2,jd) - pwpb(s_t)) / &
                                            (fcb(s_t) - pwpb(s_t)),0.0_ireals)
              ENDIF
            ELSE
              w_g1_gme(j1,j2,jd)   = 0.0_ireals
              w_g2_gme(j1,j2,jd)   = 0.0_ireals
              IF (nl_soil_in == 3) THEN
                w_g3_gme(j1,j2,jd) = 0.0_ireals
              ENDIF
            ENDIF ! s_t
          ENDIF   ! soiltyp
        ENDDO
      ENDDO
    ENDDO
  ELSE
    ! Scale the soil moisture by the pore volume of soil
    ! (for soiltypes #3 to #8; set soil moisture to 0 otherwise
    DO jd = jd_min, jd_max
      DO j2 = igg2sm2, igg2ep2
        DO j1 = igg1sm2, igg1ep2
          IF ((soiltyp_gme(j1,j2,jd) /= undef) .AND.                       &
                 (w_g1_gme(j1,j2,jd) /= undef)) THEN
            s_t   = NINT (soiltyp_gme(j1,j2,jd))
            IF (s_t >= 3 .AND. s_t <= 8) THEN
              w_g1_gme(j1,j2,jd)   = w_g1_gme(j1,j2,jd)/porb(s_t)
              w_g2_gme(j1,j2,jd)   = w_g2_gme(j1,j2,jd)/porb(s_t)
              IF (nl_soil_in == 3) THEN
                w_g3_gme(j1,j2,jd) = w_g3_gme(j1,j2,jd)/porb(s_t)
              ENDIF
            ELSE
              w_g1_gme(j1,j2,jd)   = 0.0_ireals
              w_g2_gme(j1,j2,jd)   = 0.0_ireals
              IF (nl_soil_in == 3) THEN
                w_g3_gme(j1,j2,jd) = 0.0_ireals
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! Interpolate w_g1:
  IF ((nl_soil_in == 3) .AND. (nl_soil_lm == 2) ) THEN
    WHERE (w_g1_gme(:,:,:) /= undef) 
      field_gme(:,:,:) = w_g1_gme(:,:,:) + w_g2_gme(:,:,:)
    ELSEWHERE
      field_gme(:,:,:) = undef
    ENDWHERE
  ELSE
    field_gme(:,:,:) = w_g1_gme(:,:,:)
  ENDIF

  ! normal linear (or match) interpolation to LM-field
  CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                    jd_min, jd_max, ispoke, baryll_m, index_m,            &
                    lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                    w_intpol, n_intpol, m_intpol, l_intpol,               &
                    undef, w_g1_lm, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
  ENDIF

  ! Interpolate w_g2:
  IF ((nl_soil_in == 2) .AND. (nl_soil_lm == 3) ) THEN
    field_gme(:,:,:) = w_g1_gme(:,:,:)
  ELSEIF ((nl_soil_in == 3) .AND. (nl_soil_lm == 2) ) THEN
    field_gme(:,:,:) = w_g3_gme(:,:,:)
  ELSE
    field_gme(:,:,:) = w_g2_gme(:,:,:)
  ENDIF

  ! normal linear (or match) interpolation to LM-field
  CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                    jd_min, jd_max, ispoke, baryll_m, index_m,            &
                    lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                    w_intpol, n_intpol, m_intpol, l_intpol,               &
                    undef, w_g2_lm, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
  ENDIF

  ! Interpolate w_g3, if needed:
  IF (nl_soil_lm == 3) THEN
    IF (nl_soil_in == 2) THEN
      field_gme(:,:,:) = w_g2_gme(:,:,:)
    ELSEIF (nl_soil_in == 3) THEN
      field_gme(:,:,:) = w_g3_gme(:,:,:)
    ENDIF
    CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,      &
                      jd_min, jd_max, ispoke, baryll_m, index_m,          &
                      lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,       &
                      w_intpol, n_intpol, m_intpol, l_intpol,             &
                      undef, w_g3_lm, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
    ENDIF
  ENDIF

  ! Scale w_g1_lm and w_g2_lm with thickness of soil moisture layers, 
  ! if going from 2 GME-layers to 3 LM-layers
  ! (i.e. going from w_g1_gme to w_g1_lm and w_g2_lm)
  IF ((nl_soil_lm == 3) .AND. (nl_soil_in == 2)) THEN
    ! scale factors
    zwg1_lm  = (var_lm(mzwg1_loc_lm)%levbot - var_lm(mzwg1_loc_lm)%levtop)
    zwg2_lm  = (var_lm(mzwg2_loc_lm)%levbot - var_lm(mzwg2_loc_lm)%levtop)
    zwg_gme  = 1.0_ireals /                                               &
               (var_in(mzwg1_loc_gme)%levbot - var_in(mzwg1_loc_gme)%levtop)
   
    ! do scaling
    w_g1_lm(:,:) = w_g1_lm(:,:) * zwg1_lm * zwg_gme
    w_g2_lm(:,:) = w_g2_lm(:,:) * zwg2_lm * zwg_gme
  ENDIF

!------------------------------------------------------------------------------
! Section 7: Interpolation of w_cl --> w_cl_lm
!------------------------------------------------------------------------------

  IF (.NOT. lcomp_bound) THEN
    yzitype    = var_in(mzwcl_loc_gme)%ipc(1:1)
    IF (var_in(mzwcl_loc_gme)%ipc(2:2) == 'T') THEN
      lzmono   = .TRUE.
    ELSE
      lzmono   = .FALSE.
    ENDIF
    IF (var_in(mzwcl_loc_gme)%ipc(3:3) == 'T') THEN
      lzposdef = .TRUE.
    ELSE
      lzposdef = .FALSE.
    ENDIF
  
    IF (l_smi) THEN
      ! Scale the soil moisture with the soil moisture index SMI:
      ! smi = (sm-pwp)/(fc-pwp)
      ! Caution: sm is in volumetric units, soil depth in cm!
      DO jd = jd_min, jd_max
        DO j2 = igg2sm2, igg2ep2
          DO j1 = igg1sm2, igg1ep2
            IF ((soiltyp_gme(j1,j2,jd) /= undef) .AND.                       &
                   (w_cl_gme(j1,j2,jd) /= undef)) THEN
              s_t    = NINT (soiltyp_gme(j1,j2,jd))
              IF (s_t >= 3 .AND. s_t <= 8) THEN
                w_cl_gme(j1,j2,jd) = w_cl_gme(j1,j2,jd) /                      &
                (var_in(mzwcl_loc_gme)%levbot-var_in(mzwcl_loc_gme)%levtop)*0.01
                field_gme(j1,j2,jd) = MAX((w_cl_gme(j1,j2,jd) - pwpb(s_t)) / &
                                           (fcb(s_t) - pwpb(s_t)), 0.0_ireals)
              ELSE
                field_gme(j1,j2,jd)   = 0.0_ireals 
              ENDIF
            ENDIF
          ENDDO 
        ENDDO
      ENDDO
    ELSE
      ! Scale the soil moisture by the pore volume of soil
      ! (for soiltypes #3 to #8; set soil moisture to 0 otherwise
      DO jd = jd_min, jd_max
        DO j2 = igg2sm2, igg2ep2
          DO j1 = igg1sm2, igg1ep2
            IF ((soiltyp_gme(j1,j2,jd) /= undef) .AND.                       &
                   (w_cl_gme(j1,j2,jd) /= undef)) THEN
              s_t    = NINT (soiltyp_gme(j1,j2,jd))
              IF (s_t >= 3 .AND. s_t <= 8) THEN
                field_gme(j1,j2,jd)   = w_cl_gme(j1,j2,jd)/porb(s_t)
              ELSE
                field_gme(j1,j2,jd)   = 0.0_ireals
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  
    ! normal linear (or match) interpolation to LM-field
    CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                      jd_min, jd_max, ispoke, baryll_m, index_m,            &
                      lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                      w_intpol, n_intpol, m_intpol, l_intpol,               &
                      undef, w_cl_lm, 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
    ENDIF
  ENDIF
ENDIF

IF (lmulti_layer_in .AND. .NOT. lcomp_bound) THEN

  !----------------------------------------------------------------------------
  ! Section 8: Interpolation of t_so(k) - t_so(k-1) --> dt_so_gl(k)
  !----------------------------------------------------------------------------

  yzitype    = var_in(mztso_loc_gme)%ipc(1:1)
  IF (var_in(mztso_loc_gme)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mztso_loc_gme)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF

  DO k = 1, ke_soil_in+1
    DO jd = jd_min, jd_max
      DO j2 = igg2sm2, igg2ep2
        DO j1 = igg1sm2, igg1ep2
          IF (t_so_gme(j1,j2,k,jd) /= undef) THEN
            field_gme(j1,j2,jd) = t_so_gme(j1,j2,k,jd) - t_so_gme(j1,j2,k-1,jd)
          ELSE
            field_gme(j1,j2,jd) = undef
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! normal linear (or match) interpolation to LM-field
    CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                      jd_min, jd_max, ispoke, baryll_m, index_m,            &
                      lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                      w_intpol, n_intpol, m_intpol, l_intpol,               &
                      undef, dt_so_gl(:,:,k), 1, ie2lm, 1, je2lm, yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
    ENDIF
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 9: Interpolation of w_so(k) --> w_so_gl(k)
  !----------------------------------------------------------------------------

  yzitype    = var_in(mzwso_loc_gme)%ipc(1:1)
  IF (var_in(mzwso_loc_gme)%ipc(2:2) == 'T') THEN
    lzmono   = .TRUE.
  ELSE
    lzmono   = .FALSE.
  ENDIF
  IF (var_in(mzwso_loc_gme)%ipc(3:3) == 'T') THEN
    lzposdef = .TRUE.
  ELSE
    lzposdef = .FALSE.
  ENDIF

  DO k = 1, ke_soil_in+1
    IF (l_smi) THEN
      ! Scale the soil moisture with the soil moisture index SMI:
      ! smi = (sm-pwp)/(fc-pwp)
      ! Caution: sm is in volumetric units, soil depth in cm!
      DO jd = jd_min, jd_max
        DO j2 = igg2sm2, igg2ep2
          DO j1 = igg1sm2, igg1ep2
            IF ((soiltyp_gme(j1,j2,jd) /= undef) .AND.                       &
                   (w_so_gme(j1,j2,k,jd) /= undef)) THEN
              s_t = NINT (soiltyp_gme(j1,j2,jd))
              IF (s_t >= 3 .AND. s_t <= 8) THEN
                IF (k == 1) THEN
                  w_so_gme(j1,j2,k,jd) = w_so_gme(j1,j2,k,jd) / 0.01
                ELSE
                  w_so_gme(j1,j2,k,jd) = w_so_gme(j1,j2,k,jd) /     &
                                           ((3**(k-1)-3**(k-2))*0.01)
                ENDIF
                field_gme(j1,j2,jd) = MAX((w_so_gme(j1,j2,k,jd) - pwpb(s_t))/ &
                                            (fcb(s_t) - pwpb(s_t)), 0.0_ireals)
              ELSE
                field_gme(j1,j2,jd)   = 0.0_ireals
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ELSE
      ! Scale the soil moisture by the pore volume of soil
      ! (for soiltypes #3 to #8; set soil moisture to 0 otherwise
      DO jd = jd_min, jd_max
        DO j2 = igg2sm2, igg2ep2
          DO j1 = igg1sm2, igg1ep2
            IF ((soiltyp_gme(j1,j2,jd) /= undef) .AND.                       &
                   (w_so_gme(j1,j2,k,jd) /= undef)) THEN
              s_t    = NINT (soiltyp_gme(j1,j2,jd))
              IF (s_t >= 3 .AND. s_t <= 8) THEN
                field_gme(j1,j2,jd)   = w_so_gme(j1,j2,k,jd)/porb(s_t)
              ELSE
                field_gme(j1,j2,jd)   = 0.0_ireals
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
! SB 23.3.2011
! ELSE       !USUS: but also see above l_smi ???
! TEST: if working with soil wetness index (corresponding variables are needed for that)
!              IF (lmulti_layer_in .AND. ke_soil_in /= ke_soil_lm) THEN
!                field_gme(j1,j2,jd) = (w_so_gme(j1,j2,k,jd)-cpwp(s_t))/(cfcap(s_t)-cpwp(s_t))&
!                    *plcmx_gme(j1,j2,jd) + w_so_gme(j1,j2,k,jd)/cfcap(s_t)*(1-plcmx_gme(j1,j2,jd))
!              ELSE
!               field_gme(j1,j2,jd)   = w_so_gme(j1,j2,k,jd)/porb(s_t)
!              ENDIF
!<-- SB

    ENDIF

    ! normal linear (or match) interpolation to LM-field
    CALL pp_interp2ls(field_gme, igg1sm2, igg1ep2, igg2sm2, igg2ep2,        &
                      jd_min, jd_max, ispoke, baryll_m, index_m,            &
                      lzmono, lzposdef, yzitype, lolp_gme, lolp_lm,         &
                      w_intpol, n_intpol, m_intpol, l_intpol,               &
                      undef, w_so_gl(:,:,k), 1, ie2lm, 1, je2lm,            &
                      yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 9520_iintegers, yzerrmsg, yzroutine)
    ENDIF

    ! If no vertical interpolation between soil layers is required, store horizontal 
    ! interpolated fields to w_so_lm
    IF (lmulti_layer_in .AND. ke_soil_in == ke_soil_lm) THEN
       w_so_lm(:,:,k) = w_so_gl(:,:,k)
    ENDIF

  ENDDO
ENDIF

!------------------------------------------------------------------------------
! Section 10: Check, whether dpsdt has been read and interpolated
!------------------------------------------------------------------------------

  IF (.NOT. ldpsdt) THEN
    ! set dpsdt_gl to 0.0
    dpsdt_gl(:,:) = 0.0_ireals
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE interpol_gme_special

!==============================================================================

END MODULE src_gme_interpol
