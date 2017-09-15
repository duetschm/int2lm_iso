!+ Source Module for reading input-fields for the coarse grid
!==============================================================================

MODULE src_read_coarse_grid

!==============================================================================
!
! Description:
!   This module contains subroutines necessary for reading the input-fields
!   for the coarse grid
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
!  T_SO(0) has also to be read for computing boundary values
! V1_5         2007/07/09 Ulrich Schaettler
!  Adaption of ready file name to LM (must be LMF)
!  Distribution of nunit_of_time only for num_compute larger than 1
!  Renamed ke_soil to ke_soil_in
!  Taking into account SLEVE and new type of coding vertical coordinate parameters
!    in get_vert_coord (Daniel Leuenberger)
!  Introduced nlevskip to  skip levels at the top of ECMWF model (Davide Cesari)
!  Added treatment of chemistry fields
! V1_6         2007/09/07 Ulrich Schaettler, Burkhardt Rockel
!  Introduced yinput_type and eliminated lanalysis
!  Introduced NetCDF-IO
! V1_7         2007/11/26 Ulrich Schaettler, Christoph Gebhardt
!  Renamed iw_so_rel_type to itype_w_so_rel to be consistent with other names
!  Added check for ensemble member ID in ensemble mode
! V1_8         2008/05/29 Ulrich Schaettler, Burkhardt Rockel, Hans-Juergen Panitz
!  Compute control level for geopotential, if it could not be read (for lcm2lm)
!  Corrected a bug for computing zprocarray in case of NetCDF
!  Eliminated interpolation of SST (for lcm2lm)
!  Check, whether the input data comes from the 'f'ull domain or a 's'ubdomain
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
! V1_9         2009/09/03 Guenther Zaengl, Ulrich Schaettler, Burkhardt Rockel
!  Adaptations for new reference atmosphere (Guenther)
!  Moved computation of reference atmosphere for COSMO input from external
!  data to here. (Uli)
!  Do global_values in get_vert_coord only if num_compute > 1 (Martin Suklitsch)
!  Correction due to date zone crossing
!  New reference atmosphere in netCDF
! V1_10        2009/12/17 Oliver Fuhrer, Ulrich Schaettler
!  Bug fix for reading pressure deviation and/or standard pressure
!  Modifications to allow processing of Unified Model data:
!    new SR check_um_grid; check_required_um; read vertical coordinate parameters
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_12        2010/06/14 Ulrich Schaettler
!  Computing the new reference atmosphere for UM data (instead of the old one)
!  Bug Fix in determining the names of the ready files
!  Allowing "O"ptional fields in list for boundary values (CLM)
! V1_13        2010/10/13 Christoph Gebhardt, Ulrich Schaettler
!  Modifications to recognize additional EPS types as boundaries
! V1_14        2010/11/19 Ulrich Schaettler, Anne Roches
!  Modifications to use grib_api
!  Modifications to allow processing of JMA and NCEP data
!  Adaptations in SR read_nc_gdefs for reading multi-dimensional external
!    parameter (HORIZON) (AR)
! V1_15        2010/12/10 Ulrich Schaettler
!  Added nlevskip for computing the correct vertical coordinate parameters
!  Removed some debug print outs
!  Added check for EPS members when using grib_api
! V1_17        2011/03/11 Ulrich Schaettler
!  Corrected allocation of zt0, zp0hl in SR org_read_coarse_grid, Section 6.
!  They have to be allocated with ke_in, ke1in instead of kelm, ke1lm
!  Added ECMWF fields crwc (qr) and cswc (qs) to list of possible input fields
! V1_19        2012/06/06 Davide Cesari, Susanne Brienen, Helmut Frank,
!                         Christoph Gebhardt, Burkhardt Rockel, Daniel Luethi
!  Bug correction for nlevskip /= 0 (Davide Cesari)
!  Translate GFS short names. (Helmut Frank)
!  Added ensemble type 203 for COSMO-LEPS (Christoph Gebhardt)
!  Implemented conditional compilation for NETCDF, GRIBDWD (Uli Schaettler)
!  Check whether resolution are in grib data using IBITS
!  Extensions for processing UM and HIRLAM data
!  Use msoilgrib_in (instead of msoilgrib) for incoming data (Susanne Brienen)
!  CLM:
!   Include checking of necessary input for lcm2lm
!   Several changes to read hybrid height coordinates from input data in netCDF format
!   Added 365_days support
!   Introduce pressure level support for climate model CM input
!   Allow the name "soil_depth" for dimension and name of soil layers in netCDF CM input
!     alternatively to "soil1"
!   Deactivate computation of control geopotential in case of pressure coordinates
!     Then this has to be computed in src_pressure_to_hybrid
!   Set reference atmosphere for lcm_hgt_coor: irefatm_in = 2 and compute p0_gl with
!     correct vertical coordinates (akh_in_rho for lum2lm, akh_in for lcm_hgt_coor)
!   Unified dimension IDs with COSMO (ID for topo corrections changed from 14 to 15)
! V1_20        2012/09/03 Ulrich Schaettler, Michael Baldauf,
!                         Burkhardt Rockel, Davide Cesari
!  Enlarged strings for date variables to 14 characters
!  Adapted calls to subroutine make_fn
!  Renamed 'grax' to 'apix' to be conform with other models
!  Adaptations to calls to modified SR reference_atmosphere_xx (MB)
!   Correction in case of east_add_in, west_add_in, south_add_in, north_add_in is not 0
!    (Burkhardt Rockel)
!  Added a grib_set 'missing_value' to allow for successive retrieve of a field with
!  the desired value for missing data, necessary in case of frames
! V1_21        2013/03/25 Ulrich Schaettler, Burkhardt Rockel
!  Adaptations to read W_SO, W_SO_ICE for Grib1 and for Grib2 (US)
!  Adaptation to new file name convention in netCDF I/O   (BR)
!  Make read of reference date in netCDF input more flexible   (BR)
!  rename several local variables in read_nc_gdefs and read_nc_vdefs to meet the coding standards
!  Changes in netCDF output: scalar variable instead of extra dimension of length 1 (e.g. height_2m)
!    needs changes in idims_id field indices    (BR)
! V1_22        2013/07/11 Ulrich Schaettler, Davide Cesari
!  Implemented additional "ifdef GRIBAPI", which had been forgotten before
!  Adapted interface to read_gribapi with special grib_api integer
!  Get ivctype from pv_in(4+ke1in+1) only, if list of vertical coordinate
!   parameters is long enough in SR get_vert_coord (it was not coded for old data)
!  Adapted interface to check_input_grid according to COSMO-Model (pv_in, inrvert_in)
!  Replaced zhsurfs_lm by hsurfs_in in calls to reference atmospheres for input grid
!  Use structures and routines for reference_atmospheres from module
!    vgrid_refatm_utils
!  Special treatment for IFS grib fields, which are not really multi-layer in the
!    COSMO sense (Davide)
! V1_23        2013/10/02 Ulrich Schaettler
!  Do the reference atmosphere calculations for the input grid only for the first call
!  Bug fixes in SR get_vertcoord:
!     added REFSTF when reading vertical coordinate parameters with grb1
!     added settings of vcoord_in-values for ivctype=1/2/3
! V1_24        2013/11/01 Ulrich Schaettler
!  Replaced kesoildim_in (not initialized) by ke1soildim_in-1 in SR read_nc_gdefs
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
  ireals,    & ! KIND-type parameters for real variables
  irealgrib, & ! KIND-type of the REALs in the grib library
  iintegers, & ! KIND-type parameter for "normal" integer variables
  intgribf,  & ! KIND-type of the fortran decks in the grib library
  intgribc,  & ! KIND-type of the c decks in the grib library
  iwlength,  & ! length of integers used in the griblib in byte
  int_ga       ! integer precision for grib_api: length of message in bytes

!------------------------------------------------------------------------------

USE data_grid_lm,    ONLY: &
  ie2lm_tot,   & ! = ielm_tot + 2
  je2lm_tot,   & ! = jelm_tot + 2
  ie2lm,       & !
  je2lm,       & !
  kelm,        & ! ke for LM
  ke1lm          ! ke+1

!------------------------------------------------------------------------------

USE data_grid_in,    ONLY: &
  pollat_in,      & ! latitude of the rotated north pole (in degrees, N>0)
  pollon_in,      & ! longitude of the rotated north pole (in degrees, E>0)
  polgam_in,      & ! latitude of the rotated north pole
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
  east_add_in,    & ! add an extra column to the East
  west_add_in,    & ! add an extra column to the West
  south_add_in,   & ! add an extra column to the South
  north_add_in,   & ! add an extra column to the North
  ie_in_tot,      & ! ie for input grid, total domain
  je_in_tot,      & ! je for input grid, total domain
  ie_in_max,      & ! Max. of ie_in on all processors
  je_in_max,      & ! Max. of je_in on all processors
  ie_in,          & ! ie for input grid, local domain
  je_in,          & ! je for input grid, local domain
  ke_in,          & ! ke for input grid
  ke_pressure,    & ! number of pressure levels when interpolating from pressure levels
  press_level,    & ! list of available pressure levels (in Pa) for GFS and CM
  klv850_in         ! approximate level where 850 hPa is reached

USE data_grid_in,    ONLY: &
  nlevskip,    & ! number of levels to skip at top of input model
  ke1in,       & ! ke+1 for input grid
  ke_soil_in,  & ! number of levels in multi-layer soil model
  czmls_in,    & ! depth of the input soil layers in meters
  czhls_in,    & ! depth of the half soil layers in meters in output
  ak_in ,      & ! vertical coordinate parameters for half levels
  bk_in ,      & !                  - " -
  akh_in,      & ! vertical coordinate parameters for main levels
  bkh_in,      & !                  - " -
  akh_in_rho,  & ! coefficients for main levels for u, v (lum2lm)
  bkh_in_rho,  & ! coefficients for main levels for u, v (lum2lm)
  dak_in,      & ! difference of coordinate parameters
  dbk_in,      & !                  - " -
  lcm_hgt_coor,& ! CM input data has hybrid height coordinates
  lcm_pres_coor  ! CM input data has pressure coordinates  !_br 14.03.12

!------------------------------------------------------------------------------

USE data_fields_lm,  ONLY: &
  p0_gl        ,   & ! ref. pres. on full levels + interpol. COARSE LM oro.(Pa)
  dp0_gl       ,   & ! reference pressure thickness of layers           ( Pa  )
  rho0_gl      ,   & ! reference density at the full model levels       (kg/m3)
  ps_gl        ,   & ! surface pressure on the interpol. GME orogr.     ( Pa  )
  hhl_gl       ,   & ! height of half-levels on the interpol. COARSE LM oro.(m)
  hsurf_gl     ,   & ! interpolated orography (for lum2lm)              (  m  )
  hsurfs_gl          ! interpolated splitted parts of coarse topo       (  m  )

!------------------------------------------------------------------------------

USE data_fields_in,  ONLY: &
  hsurf_in,    & ! orography                                   (  m  )
  hsurfs_in,   & ! splitted parts of the coarse orography (SLEVE)
  hhl_in,      & ! height of half levels of coarse grid        (  m  )
  p0_in,       & ! reference pressure on coarse grid
  lolp_in,     & ! Land Sea Mask of input for 'M'atch Interpolation
  qv_in,       & ! specific water vapor content              (kg/kg)
  qc_in,       & ! specific cloud water content              (kg/kg)
  ps_in,       & ! surface pressure                          ( Pa  )
  t_in,        & ! temperature                               (  K  )
  fic_in,      & ! control level for geopotential            (m2/s2)
  fis_in         ! orography * g

!------------------------------------------------------------------------------

USE data_int2lm_constants,     ONLY :  &
    R_d,             & ! gas constant for dry air                      [J/K*kg]
    Rvd_m_o,         & ! = R_v/R_d - 1.0,
    G                  ! gravity at sea level                          [ms-2]

!------------------------------------------------------------------------------

USE data_int2lm_control,    ONLY: &
  noutput,         & ! unit number for output file
  linitial,        & ! if .TRUE., initial data for LM
  lprog_qi,        & ! if .TRUE., interpolate qi from coarse model to LM grid
  nincbound,       & ! time step increment between two datasets
  lcomp_bound,     & ! compute fields for boundaries
  leps_bc,         & ! if .TRUE., ensemble mode for boundary data
  ltime,           & ! detailled timings of the program are given
  itype_calendar,  & ! for specifying the calendar used
  yinput_model,    & ! string to identify the input model
  lgsm2lm,         & ! if .TRUE., gsm ->lm
  lgfs2lm,         & ! if .TRUE., gfs ->lm
  lec2lm,          & ! if .TRUE., ec ->lm
  llm2lm,          & ! if .TRUE., lm ->lm
  lum2lm,          & ! if .TRUE., um ->lm
  lhir2lm,         & ! if .TRUE., hirlam ->lm
  lcm2lm,          & ! if .TRUE., climate model ->lm   !_br
  msoilgrib_in,    & ! grib coded depth of main soil levels in centimeters
                     ! (careful: the first level will be coded with 1,
                     !           but is in the depth of 0.5 cm!)
  lmulti_layer_lm, & ! if .TRUE., compute soil fields for multi-layer soil
                     ! model in the outgoing data
  lmulti_layer_in, & ! if .TRUE., incoming data have soil fields from the
                     ! multi-layer soil model
  luse_t_skin,     & ! if .TRUE., use ECMWF skin temperature for surface
  lradtopo,        & ! process parameters for topographic correction of radiation
  dt,              & ! time step used in the LM
  yakdat1,         & ! actual date (ydate_ini+ntstep/dt)
  yakdat2,         & ! actual date (ydate_ini+ntstep/dt) in the form
  timings,         & ! for storing the times for different parts of the program
  pcontrol_fi,     & ! pressure of control level for geopotential  !_br
  itype_w_so_rel,  & ! type of relative soil moisture input
  idbg_level,      & ! to control verbosity of output
  lprintdeb_all      ! whether all PEs print debug output

!------------------------------------------------------------------------------

USE data_int2lm_io,        ONLY : &
  nunit_of_time,& ! indicator for unit-of-time (1hr, 15min, 30min,...)
  nincwait,     & ! if ready-file is not available wait nincwait seconds
                  ! until next attempt
  nmaxwait,     & ! if ready-file is not available after nmaxwait seconds,
                  ! abort the program
  ytrans_in,    & ! directory for reading ready-files
  ytrans_out,   & ! directory for writing ready-files
  lchkin,       & ! logical for print of check-values (max,min,mean)
                  ! of input-fields
  yin_cat,      & ! catalog-name of the input files
  nuchkdat,     & ! checking the I/O data
  ract_hour,    & ! actual hour of actual day (returned by get_utc_date)
  yuchkdat,     & ! checking the I/O data
  ymode_read,   & ! mode for opening the (read) Grib files
  yin_form_read,& ! input format of boundary data
  ylmext_cat,   & ! catalog-name of the file with external LM parameters
  yin_hhl,      & ! name of the file with COSMO-input HHL fields
  ylevltypes_in,& ! to convert GRIB1 level types to grib_api string typeOfLevel
                  ! for incoming data
  rscalefac,    & ! Array to convert GRIB2 scale factors to real numbers
  npds,         & ! Dimension for product definition section (pds)
  ngds,         & ! Dimension for grid description section (gds)
  nbms,         & ! Dimension for bit map section (bms)
  nbds,         & ! Dimension for binary data section
  ndsup,        & ! Dimension for dsup
  ndims,        & ! Dimension for idims (contains all dimensions)
  lfd,          & !
  nbitmap,      & !
  lds             !

USE data_int2lm_io,        ONLY : &
  nvar_lm,      & ! maximum number of variables in fine grid LM variable table
  nvar_in,      & ! maximum number of variables in coarse grid variable table
  nvar_in_norm, & ! maximum number of variables without chemistry
  idwdednr,     & ! grib edition number for DWD grib library
  undefgrib,    & ! value for "undefined" in the grib routines
  undefncdf,    & ! value for "undefined" in the netcdf routines
  undef,        & ! the same with other KIND-Parameter
  ytunit_in,    & ! time unit for input data
  yinput_type,  & ! type of input data: 'forecast', 'analysis' or 'ana_init'
  ydate_ini,    & ! start of the forecast yyyymmddhh (year,month,day,hour)
  lmmss_ini,    & ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                  ! if .FALSE. 10 digits date format (YYYYMMDDHH)
                  ! for ydate_ini and result files of INT2LM
  ydate_bd,     & ! start of the forecast from which theboundary fields are used
  lmmss_bd,     & ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                  ! if .FALSE. 10 digits date format (YYYYMMDDHH)
                  ! for ydate_bd  and input  files of INT2LM
  iblock,       & ! array for gribed data
  idims_in,     & ! array for all dimensions
  ibmap,        & ! array for
  ipds,         & ! product definition section
  igds_in,      & ! grid description section
  ibms,         & ! bit map section
  ibds,         & ! binary data section
  dsup,         & ! Parameter for grib routines
  inrvert_in,   & ! number of vertical coordinate parameters of input data
  pv_in,        & ! array for vertical coordinate parameters of input data
  var_lm,       & ! array for fine grid LM variable table
  var_in,       & ! array for input model variable table
  ar_des_input, & ! structure for input variable table
  ndims_id,     & ! array for the IDs of the dimensions of netCDF formatted output
  idims_id        ! array for the IDs of the dimensions

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY :  &
    nproc,           & ! total number of processors: nprocx * nprocy
    lasync_io,       & ! if .TRUE.: the model runs with extra PEs for
                       ! asynchronous IO
    num_compute,     & ! number of compute PEs
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    isubpos_coarse,  & ! positions of the subdomains of the coarse grid
                       ! in the total domain.
    icomm_cart,      & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_grib,        & ! determines the REAL type used for the GRIB library
    imp_integers,    & ! determines the correct INTEGER type used in the model
                       ! for MPI
    imp_logical,     & ! determines the correct LOGICAL   type used in the
                       ! model for MPI
    imp_character,   & ! determines the correct CHARACTER type used in the model
                       ! for MPI
    imp_byte           ! determines the correct BYTE type used in the model
                       ! for MPI

!------------------------------------------------------------------------------

USE data_epsctl,     ONLY :   &
    iepsmem_bc,    & ! ID of the member in the ensemble of boundary
                     ! conditions (iepsmem_bc >= 0)
    iepstyp_bc,    & ! ID of the boundary ensemble generation type
                     ! (iepstyp_bc >= 0)
    iepstot_bc,    & ! total number of boundary ensemble members (iepstot_bc >= 0)
    lchk_bc_typ      ! if .TRUE. AND leps_bc=.TRUE., check member IDs of input data

!------------------------------------------------------------------------------

USE utilities,            ONLY : elapsed_time, sleve_split_oro, diff_minutes
USE vgrid_refatm_utils,   ONLY : reference_atmosphere, reference_atmosphere_2,  &
                                 k_index_of_pressure_levels,                    &
                                 lanalyt_calc_t0p0, vcoord_in, refatm_in,       &
                                 rundefined, nfltvc_in, svc1_in, svc2_in,       &
                                 lhhl_in_read, uuid_in_string, uuid_2char,      &
                                 lhhl_hasbeenread
USE environment,          ONLY : comm_barrier, model_abort
USE parallel_utilities,   ONLY : distribute_values, gather_values,              &
                                 scatter_values, gather_field, distribute_field,&
                                 global_values
USE io_utilities,         ONLY : open_file, read_grib, read_netcdf, close_file, &
                                 make_fn, check_record, check_input_grid,       &
                                 read_gribapi
USE src_read_hhl,         ONLY : org_read_hhl

#ifdef NETCDF
USE netcdf,           ONLY :   &
  nf90_enotvar,            &
  nf90_enotatt,            &
  nf90_get_att,            &
  nf90_get_var,            &
  nf90_inq_dimid,          &
  nf90_inq_varid,          &
! nf90_inquire,            &
  nf90_inquire_dimension,  &
  nf90_noerr,              &
  nf90_strerror
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! string variable to hold grid information
  CHARACTER (LEN=200)     grid_mapping

!==============================================================================

CONTAINS

!==============================================================================
!+ Organizes the reading and interpolation of coarse grid fields
!------------------------------------------------------------------------------

SUBROUTINE org_read_coarse_grid (nnow, lfirst, ydate)

!------------------------------------------------------------------------------
!
! Description:
!  Coarse grid fields are read and distributed to the available PEs. Every PE
!  gets a different record for grib-unpacking. The unpacked records are
!  distributed to the PEs so that every PE can get the data it needs for the
!  interpolations (the necessary values have been calculated in
!  decomp_coarse_grid in the setup.
!
! Method:
!  After creating the filename and opening the grib file, all records are read
!  in a read_loop.
!
! Input files:
!  Grib-file with input-fields.
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER  (KIND=iintegers),  INTENT(IN)    ::  &
  nnow                    ! actual time step to be processed

LOGICAL,                    INTENT(INOUT) ::  &
  lfirst                  ! for actions that have to be done only once

CHARACTER (LEN=14),         INTENT(IN)    ::  &
  ydate           ! actual date in the form   yyyymmddhhmmss

!------------------------------------------------------------------------------

! Local arrays
REAL (KIND=ireals)         ::  &
  field_glob      (ie_in_tot,je_in_tot),   &
  field_glob_aux  (ie_in_tot,je_in_tot),   &
  zprocarray      (ie_in_max*je_in_max, 0:num_compute-1), &
  zhhlr_in(ke1in), zak(ke1in), zbk(ke1in)

REAL (KIND=irealgrib)      ::  &
  field_grib((ie_in_tot+1)*(je_in_tot+1))   ! field read from the GRIB-file

REAL (KIND=irealgrib)      :: zundef

! to compute control geopotential, if not available
REAL (KIND=ireals) , ALLOCATABLE                     ::   &
  ztv    (:,:),          & ! virtual temperature
  zpo    (:,:),          & !
  zpu    (:,:),          & !
  zfiu   (:,:),          & !
  zfio   (:,:)             !

REAL (KIND=ireals), ALLOCATABLE  ::  &
  zhsurf_lm_tot (:,:)  , & ! full topography of full domain
  zrho0         (:,:,:), & ! intermediate storage
  zdp0          (:,:,:), &
  zp0hl         (:,:,:), &
  zt0hl         (:,:,:), &
  zt0           (:,:,:)

! Local variables:
INTEGER  (KIND=intgribf)   :: &
  ierrf,                      & ! error code for grib routines
  iniyy, inimm, inidd, inihh, inimin, inisec, & !
  ibdyy, ibdmm, ibddd, ibdhh, ibdmin, ibdsec, & !
  imindif, igribid, ireturn, jzscanmode, igriblen, &
  ilevel, itoplevel, ibottomlevel, &
  idisc, icatg, ipara

INTEGER  (KIND=iintegers)  :: &
  izvar_count, izlev_count,   & ! for NetCDF Organization
  myzvar, myzlev, myzlevtot,  & ! organization indices returned by write_netcdf
  k

INTEGER (KIND=iintegers)   :: &
  ivar_id(nvar_in)              ! returned list of IDs from read_nc_vdefs

INTEGER  (KIND=iintegers)  ::  &
  mztskin_loc_in, mzts_loc_in, & ! location of TSKIN and T_S in the varaible table
  mztso_loc_in,                & ! location of T_SO in the variable table
  mzt_loc_in,                  & ! location of T in the variable table
  mzsst_loc_in,                & ! location of SST in the variable table
  mzp_loc_in, mzpp_loc_in,     & ! location of P and PP in the variable table
  mzfi_loc_in                    ! location of FI (on control level) in the variable table

INTEGER  (KIND=intgribc)   :: &
  nzincwaitc                ! as nincwait, but for calling C-routine of grib lib.

INTEGER  (KIND=iintegers)  ::  &
  nufile,                 & ! unit number of opened grib file
  izerror,                & ! status and error status variable
  niostat,                & ! status and error status variable
  izdebug,                & ! for verbosity of output
  izlen,                  & ! length of path-name for ready-files
  nzwait,                 & ! seconds waited so far
  izsize,                 & ! length of grib record read
  izloc,                  & ! location of record in variable table
  izlev,                  & ! level of a mlf field
  izrank,                 & ! rang of the variable
  iz_ps_idx, iz_lnps_idx, & !
  kflat, klv850, klv950,  & !
  i,j, ij, jscan, n, nzbyte, lfdec, ldsec, izstat, iproc, ndiff_ini_bd, &
  ihour, imin, isec, isecleft, ie_p, je_p, iz1, nzlevels, nzlevrea,     &
  nzscan, igrbednr, itabtyp, iee, ilevtyp, isec1len,                    &
  ilocdef, ilocensid, iloctotnr, ilocensnr, itypeoflevel, je_in_end,    &
  ndiff_ini_bdref

INTEGER (KIND=int_ga)    ::  &
  izmaxlen, izsize_ga

REAL    (KIND=ireals)      ::  &
  zfactor, zbias, zrealdiff, rsecleft, zgdrt, ztdbe, zbetf, zt00, zpnu, zpno

LOGICAL                    ::  &
  lzconvert,  & ! for conversion from LNPS to PS
  lzfirst_uv, & ! for reading vertical coordinate parameters for UM winds
  lzgot_genvc,& ! for reading vertical coordinate parameters for general vertical coord.
  lzrefatmos, & ! to indicate whether reference atmospheres have to be computed
  lzrequired, & ! indicates whether a record from the grib file is required
  lzeof,      & ! indicates the end of file
  lzexist,    & ! to check, whether different files exist
  lzcheck,    & ! to check, whether all data have been read
  lzwarn,     & ! to warn, if chemistry variables are not available
  lzhhl_in,   & ! to indicate, if GRIB2 HHL has to be read
  lzreadin

CHARACTER (LEN=250)        ::  &
  yzpath      ! full path and name of the input-file

CHARACTER (LEN= 14)        ::  &
  yzdate      ! actual date used for make_fn

CHARACTER (LEN= 14)        ::  &
  yzfulldate  ! actual date used for the check routines, including minutes
              ! and seconds

CHARACTER (LEN= 30)        ::  &
  ytypeoflevel ! typeOfLevel

CHARACTER (LEN= 21)        ::  &
  yshortname, yzname  ! short name from grib_api

CHARACTER (LEN=  3)        ::  &
  yzdattyp    ! actual date used for make_fn

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN= 15)        ::  &
  yzfn        ! name of files to be opened

CHARACTER (LEN=200)        ::  &
  yzerrmsg, & ! error message for error handling
  yzerrms1, & ! error message for error handling
  yzerrms2    ! error message for error handling

CHARACTER (LEN=  3)        ::  &
  yzhead      ! file name header for input-file

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  igriblen = 0
  yzname = '                    '
  jzscanmode = 1
  zprocarray(:,:) = 0.0_ireals

  IF (lfirst) THEN
    ! reference atmospheres need only be computed once
    lzrefatmos = .TRUE.
  ELSE
    lzrefatmos = .FALSE.
  ENDIF

  IF (lgfs2lm .OR. (lcm2lm .AND. lcm_pres_coor)) THEN
    ! set number of vertical levels to pressure levels again
    ke_in = ke_pressure

    DO n = 1, nvar_in
      IF (var_in(n)%rank == 3) var_in(n)%nlevels = ke_in
    ENDDO

    ! Nullify humidities in top layers, because they are not read in
    qv_in(:,:,:) = 0.0_ireals
    qc_in(:,:,:) = 0.0_ireals
  ENDIF

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '  Start read_coarse_grid'
  ENDIF

  ! determine actual minute and second of the step that has to be processed:
  ! this is needed in the check-routines, to decide whether the data that is
  ! read is really valid for this step (time)
  yzfulldate(1:14) = ydate(1:14)

  ! initialize level counter and logical flags
  lzfirst_uv  = .TRUE.
  lzgot_genvc = .FALSE.
  izerror     = 0_iintegers
  yzerrmsg    = '          '
  lzeof       = .FALSE.
  lzrequired  = .FALSE.
  lzexist     = .FALSE.
  ierrf       = 0_intgribf
  yzroutine   = 'org_read_coarse_grid'
  nzincwaitc  = INT (nincwait, intgribc)
  IF ((yin_form_read == 'grb1') .OR. (yin_form_read == 'apix')) THEN
    zundef      =       undefgrib
    undef       = REAL (undefgrib, ireals)
  ELSEIF (yin_form_read == 'ncdf') THEN
    zundef      =       undefncdf
    undef       = REAL (undefncdf, ireals)
  ENDIF

  ! Set dimensions for grib variables and allocate iblock, ibmap and dsup
  ! (the role of ds is taken by field_grib).
  ! nzbyte is assumed to be 2 here: this is true in all our cases, but if
  ! if we once change the packing rate of grib code (nrbit=16) we will get
  ! in trouble here. To know this a priori, the pds of the first record has
  ! to be read and decoded to get nrbit. The 2000 are just a safety-add to
  ! take care of the definition sections that are also stored in iblock
  ! for safety: increase to 4
  nzbyte = 8

  ldsec = (ie_in_tot+1) * (je_in_tot+1)
  lfdec = ldsec * nzbyte / iwlength + 2000

  lfd     = INT (lfdec, intgribf)
  lds     = INT (ldsec, intgribf)

  ! Correct idims_in
  idims_in( 1)   = npds
  idims_in( 2)   = ngds
  idims_in( 3)   = nbms
  idims_in( 4)   = nbds
  idims_in( 5)   = nbitmap
  idims_in( 6)   = ndsup
  idims_in( 7)   = lds
  idims_in( 8)   = lfd
  idims_in(9:20) = 0

  ALLOCATE (iblock(lfd), ibmap(nbitmap), STAT=izstat)
  ALLOCATE (dsup(ndsup),                 STAT=izstat)

  ! Find some special fields in var_in table
  iz_ps_idx = 0 ; iz_lnps_idx = 0; mzfi_loc_in = 0;
  DO n = 1, nvar_in
    IF (TRIM(var_in(n)%name) == 'PS')   iz_ps_idx   = n
    IF (TRIM(var_in(n)%name) == 'LNPS') iz_lnps_idx = n
    IF (TRIM(var_in(n)%name) == 'FI')   mzfi_loc_in = n
  ENDDO

  ! Clear the data checking flags
  DO n=1,nvar_in
    IF ((var_in(n)%dattyp(1:1) == 'I' .AND. .NOT. lcomp_bound) .OR. &
        (var_in(n)%dattyp(2:2) == 'B' .AND.       lcomp_bound) ) THEN
      var_in(n)%lreadin = .FALSE.
      var_in(n)%nlevels_read  = 0
    ENDIF
  ENDDO

  ! look for some special variables
  IF (yin_form_read == 'ncdf') THEN
    DO n = 1, nvar_in
      IF (TRIM(var_in(n)%name) == 'T_SKIN') mztskin_loc_in  = n
      IF (TRIM(var_in(n)%name) == 'T_S')    mzts_loc_in  = n
      IF (TRIM(var_in(n)%name) == 'T')      mzt_loc_in   = n
      IF (lcm2lm .AND. TRIM(var_in(n)%name) == 'SST')    mzsst_loc_in = n
      IF (lcm2lm .AND. TRIM(var_in(n)%name) == 'P' )     mzp_loc_in   = n
      IF (llm2lm .AND. TRIM(var_in(n)%name) == 'P' )     mzp_loc_in   = n
      IF (llm2lm .AND. TRIM(var_in(n)%name) == 'PP')     mzpp_loc_in  = n
    ENDDO
  ELSE
    ! for grb1 we have to take care of P / PP
    DO n = 1, nvar_in
      IF (llm2lm .AND. TRIM(var_in(n)%name) == 'P' )     mzp_loc_in   = n
      IF (llm2lm .AND. TRIM(var_in(n)%name) == 'PP')     mzpp_loc_in  = n
    ENDDO
  ENDIF

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
    ndiff_ini_bd =  - NINT( (imindif * 60 + (ibdsec - inisec)) / dt , iintegers )
    ! correct ndiff_ini_bd if ydate_ini not equal to reference time of
    ! boundary data
    ! account for seconds here!
    ndiff_ini_bdref = NINT((inimin * 60) / dt, iintegers)
    ndiff_ini_bd = ndiff_ini_bd - ndiff_ini_bdref

    yzdate = ydate_bd
  ELSE
    ndiff_ini_bd = 0
    yzdate = ydate_ini
  ENDIF

!-------------------------------------------------------------------------------
! Section 2: Check if ready-files are available
!-------------------------------------------------------------------------------

  ! If required, check whether ready-files are available
  IF (ytrans_in /= '   ') THEN
    IF (izdebug > 10) THEN
      PRINT *, '  Check ready files'
    ENDIF

    IF (my_cart_id == 0) THEN
      izlen   = LEN_TRIM(ytrans_in)
      nzwait  = 0
      lzexist = .FALSE.

      ! Create the file name for the input ready-file
      IF (lec2lm) THEN
        IF   ( (yinput_type == 'forecast') .OR. (yinput_type == 'ana_init') ) THEN
          yzhead = 'EC_'
        ELSEIF (yinput_type == 'analysis') THEN
          yzhead = 'ECA'
        ENDIF
      ELSEIF (lgsm2lm) THEN
        IF   ( (yinput_type == 'forecast') .OR. (yinput_type == 'ana_init') ) THEN
          yzhead = 'JMF'
        ELSEIF (yinput_type == 'analysis') THEN
          yzhead = 'JMA'
        ENDIF
      ELSEIF (lgfs2lm) THEN
        IF   ( (yinput_type == 'forecast') .OR. (yinput_type == 'ana_init') ) THEN
          yzhead = 'NCF'
        ELSEIF (yinput_type == 'analysis') THEN
          yzhead = 'NCA'
        ENDIF
      ELSEIF (llm2lm) THEN
        IF   ( (yinput_type == 'forecast') .OR. (yinput_type == 'ana_init') ) THEN
          yzhead = 'LMF'
        ELSEIF (yinput_type == 'analysis') THEN
          yzhead = 'LMA'
        ENDIF
      ELSEIF (lum2lm) THEN
        IF   ( (yinput_type == 'forecast') .OR. (yinput_type == 'ana_init') ) THEN
          yzhead = 'UMF'
        ELSEIF (yinput_type == 'analysis') THEN
          yzhead = 'UMA'
        ENDIF
      ELSEIF (lhir2lm) THEN
        IF   ( (yinput_type == 'forecast') .OR. (yinput_type == 'ana_init') ) THEN
          yzhead = 'HMF'
        ELSEIF (yinput_type == 'analysis') THEN
          yzhead = 'HMA'
        ENDIF
      ELSEIF (lcm2lm) THEN
        IF   ( (yinput_type == 'forecast') .OR. (yinput_type == 'ana_init') ) THEN
          yzhead = 'CM_'
        ELSEIF (yinput_type == 'analysis') THEN
          yzhead = 'CMA'
        ENDIF
      ENDIF
      CALL make_fn (yzhead, yakdat1, ydate_ini, 'f',' ', nnow+ndiff_ini_bd, dt, .TRUE., &
                    itype_calendar, ytrans_in, yzpath, lmmss_bd, izdebug, izerror)

#ifdef GRIBDWD
      DO WHILE ( (.NOT. lzexist) .AND. (nzwait < nmaxwait) )
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

  ! ele begin
  ELSE

    IF (izdebug > 10) THEN
      PRINT *, '  Create file name and check availability of file'
    ENDIF

    IF (my_cart_id == 0) THEN
      nzwait  = 0
      lzexist = .FALSE.

      IF (llm2lm) THEN
        yzhead = 'l f'
      ELSEIF (lec2lm) THEN
        yzhead = 'e s'
      ELSEIF (lcm2lm) THEN
        yzhead = 'c s'
      ENDIF

      IF (yinput_type == 'forecast') THEN
        IF ((nnow == 0) .AND. (ydate_ini == ydate_bd)) THEN
          ! header could be of the form 'xif' or 'xff'
          ! Test whether a "xif" file is present, if not take a xff-file
          WRITE (yzhead(2:2), '(A1)') 'i'
          CALL make_fn (yzhead,yakdat1,ydate_ini,'f',' ',nnow+ndiff_ini_bd, dt,.TRUE., &
                        itype_calendar, yin_cat, yzpath, lmmss_ini, izdebug, izerror   )
          INQUIRE (FILE=yzpath, EXIST=lzexist)
          IF (.NOT. lzexist) THEN
            WRITE (yzhead(2:2), '(A1)') 'f'
          ENDIF
        ELSE
          ! for all other times, the header is xff
          WRITE (yzhead(2:2), '(A1)') 'f'
        ENDIF
      ELSEIF (yinput_type == 'analysis') THEN
        WRITE (yzhead(2:2), '(A1)') 'a'
      ELSEIF (yinput_type == 'ana_init') THEN
        WRITE (yzhead(2:2), '(A1)') 'i'
      ENDIF

      ! Create the file name
      CALL make_fn (yzhead,yakdat1,ydate_ini,'f',' ',nnow+ndiff_ini_bd, dt,.TRUE., &
                    itype_calendar, yin_cat, yzpath, lmmss_ini, izdebug, izerror   )

      IF (yin_form_read == 'ncdf') THEN
        yzpath = TRIM(yzpath)//'.nc'
      ENDIF

      DO WHILE ( (.NOT. lzexist) .AND. (nzwait <= nmaxwait) )
        INQUIRE (FILE=yzpath, EXIST=lzexist)
        IF (.NOT. lzexist) THEN
          ! wait nincwait seconds and try again
! gdm begin
          PRINT *, 'file not available: ', yzpath
          izlen   = LEN_TRIM(yzpath)
          PRINT *, 'file not available: ', yzpath(1:izlen)
! gdm end
          PRINT *, '      sleep ', nincwait,' seconds'
          CALL fsleep (nzincwaitc)
          nzwait = nzwait + nincwait
        ENDIF
      ENDDO
    ENDIF

    IF (num_compute > 1) THEN
      ! Synchronize the processes again
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)

      ! Distribute lzexist to all nodes
      CALL distribute_values (lzexist, 1, 0, imp_logical, icomm_cart, izerror)
    ENDIF

    IF (.NOT. lzexist) THEN
      yzerrmsg  = ' *** ERROR:  file not available: '//yzpath(1:izlen+15)
      CALL model_abort (my_cart_id, 2006, yzerrmsg, yzroutine)
    ENDIF

  ! ele end
  ENDIF

!-------------------------------------------------------------------------------
! Section 3: Create grib file name and open the file
!-------------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '  Create file name and open the file'
  ENDIF

  ! Set header (first 3 characters)
  ! 1st Character determines the model:
  IF     (llm2lm)  THEN
    yzhead(1:1) = 'l'
  ELSEIF (lec2lm)  THEN
    yzhead(1:1) = 'e'
  ELSEIF (lgfs2lm) THEN
    yzhead(1:1) = 'n'
  ELSEIF (lgsm2lm) THEN
    yzhead(1:1) = 'j'
  ELSEIF (lum2lm)  THEN
    yzhead(1:1) = 'u'
  ELSEIF (lhir2lm)  THEN
    yzhead(1:1) = 'h'
  ELSEIF (lcm2lm)  THEN
    yzhead(1:1) = 'c'
  ENDIF

  ! 2nd Character determines forecast or analyses mode:
  IF     (yinput_type == 'forecast') THEN
    IF ((nnow == 0) .AND. (ydate_ini == ydate_bd)) THEN
      ! header could be of the form 'xif' or 'xff'
      ! Test whether a "xif" file is present, if not take a xff-file
      yzhead(2:2) = 'i'
      yzhead(3:3) = 'f'
      CALL make_fn (yzhead, yakdat1, ydate_ini, ytunit_in,' ', nnow+ndiff_ini_bd, dt, &
                    .TRUE., itype_calendar, yin_cat, yzpath, lmmss_bd, izdebug, izerror)
      INQUIRE (FILE=yzpath, EXIST=lzexist)
      IF (.NOT. lzexist) THEN
        WRITE (yzhead(2:2), '(A1)') 'f'
      ENDIF
    ELSE
      ! for all other times, the header is xff
      WRITE (yzhead(2:2), '(A1)') 'f'
    ENDIF
  ELSEIF (yinput_type == 'analysis') THEN
    WRITE (yzhead(2:2), '(A1)') 'a'
  ELSEIF (yinput_type == 'ana_init') THEN
    WRITE (yzhead(2:2), '(A1)') 'i'
  ENDIF

  ! 3rd Character determines the domain: 'f'ull or 's'ub
  yzhead(3:3) = 'f'

  ! Test whether this file exists, if not, take a 's'ubdomain
  CALL make_fn (yzhead, yakdat1, ydate_ini, ytunit_in,' ', nnow+ndiff_ini_bd, dt,  &
                .TRUE., itype_calendar, yin_cat, yzpath, lmmss_bd, izdebug, izerror)
  IF (yin_form_read == 'ncdf') THEN
    yzpath = TRIM(yzpath)//'.nc'
  ENDIF

  INQUIRE (FILE=TRIM(yzpath), EXIST=lzexist)

  IF (.NOT. lzexist) THEN
    IF (my_cart_id == 0) THEN
      yzerrms1 = ' *** ATTENTION: File could not be opened:  '//TRIM(yzpath)
    ENDIF
    yzhead(3:3) = 's'
    ! Create the file name again
    CALL make_fn (yzhead, yakdat1, ydate_ini, ytunit_in,' ', nnow+ndiff_ini_bd, dt,    &
                  .TRUE., itype_calendar, yin_cat, yzpath, lmmss_bd, izdebug, izerror)
    IF (yin_form_read == 'ncdf') THEN
      yzpath = TRIM(yzpath)//'.nc'
    ENDIF
    IF (my_cart_id == 0) THEN
      yzerrms2 = ' *** ATTENTION: Try to look for a file containing a subdomain only: '//TRIM(yzpath)
    ENDIF
  ENDIF

  ! All processors have to call the routine open_file. What the parallel
  ! program really does is determined in the routine.
  CALL open_file(nufile, yzpath, ymode_read, yin_form_read, icomm_cart,  &
                 my_cart_id, num_compute, lasync_io, idbg_level,         &
                 yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    IF (my_cart_id == 0) THEN
      PRINT *, TRIM(yzerrms1)
      PRINT *, TRIM(yzerrms2)
    ENDIF
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
      izerror  = 2010
      yzerrmsg = ' ERROR  *** Error while opening file YUCHKDAT:     *** '
      WRITE (yzerrmsg(48:50),'(I3)') niostat
      CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    ! Write a headline in YUCHKDAT
    WRITE (nuchkdat,'(A)') '       '
    WRITE (nuchkdat,'(A)') '       '
    WRITE (nuchkdat,'(A,A)') 'Check input file:  ', yzpath(1:LEN_TRIM(yzpath))
    WRITE (nuchkdat,'(A,I5,A,I5)')                                       &
      '    ie_in_tot = ',ie_in_tot,'  je_in_tot = ',je_in_tot
    WRITE (nuchkdat,'(A)') '    '
    WRITE (nuchkdat,'(A,A)')                                             &
              '   var        ee  lev        min    ',                    &
              'imin jmin               max    imax jmax              mean'
    WRITE (nuchkdat,'(A)') '  '
  ENDIF

#ifdef NETCDF
  IF (yin_form_read == 'ncdf') THEN
    ! read global attributes and definitions for NetCDF

    CALL read_nc_gdefs (                                                     &
        nufile, ie_in_tot, je_in_tot, ke_in, ke_soil_in, yakdat1,            &
        startlon_in_tot, startlat_in_tot, endlon_in_tot, endlat_in_tot,      &
        pollon_in, pollat_in, polgam_in,                                     &
        east_add_in, west_add_in, south_add_in, north_add_in,                &
        icomm_cart, my_cart_id, num_compute, lfirst, .TRUE., .TRUE.,         &
        yzerrmsg, izerror)

    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 2014, yzerrmsg, yzroutine)
    ENDIF

    CALL read_nc_vdefs (                                                     &
         nufile, var_in, nvar_in, ivar_id, pollon_in, pollat_in, pcontrol_fi,&
         lzcheck, icomm_cart, my_cart_id, num_compute, yzerrmsg, izerror)

    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 2015, yzerrmsg, yzroutine)
    ENDIF

    ! the check, whether all variables have been read is done at the end
    ! of the subroutine
  ENDIF
#endif

!-------------------------------------------------------------------------------
! Section 4: (Endless) loop over all records in the grib file
!-------------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '  Start read_loop'
  ENDIF

  ! for NetCDF files this is nevertheless done in an ordered manner by going
  ! through all variables in a list. Some organizational variables have to
  ! be set for this:
  izvar_count = 1
  izlev_count = 0

  read_loop: DO WHILE (.NOT. lzeof)

    nzscan = -1 ! set scanning mode to undefined

  !-----------------------------------------------------------------------------
  ! 4.1: Get a record
  !-----------------------------------------------------------------------------

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (12) = timings(12) + zrealdiff
    ENDIF

    ! Every PE gets one record from the file in an ordered manner
    ! (determined by the rank in the communicator icomm_rank). How this is
    ! done exactly is determined in the routine read_"format". This routine has
    ! to be called by every PE.

    SELECT CASE (yin_form_read)

    CASE ('grb1')

      IF (izdebug >= 20) THEN
        PRINT *, '      Calling read_grib'
      ENDIF

      CALL read_grib (nufile, lfdec*iwlength, lfdec, icomm_cart,     &
                      iblock, izsize, num_compute, lasync_io, yzerrmsg, izerror)

      IF (idbg_level >= 20) THEN
        IF (izsize == 0) THEN
          PRINT *, '       EOF readched'
        ELSE
          PRINT *, '       Read a record with size ', izsize
        ENDIF
      ENDIF

    CASE ('apix')

      IF (izdebug >= 20) THEN
        PRINT *, '      Calling read_gribapi'
      ENDIF

      izmaxlen  = INT (lfdec*iwlength, int_ga)
      CALL read_gribapi (nufile, izmaxlen, lfdec, icomm_cart, iblock,    &
                         izsize_ga, num_compute, lasync_io, yzerrmsg, izerror)

      IF (izsize_ga == 0) THEN
        izsize = 0_iintegers
        IF (idbg_level >= 20) THEN
          PRINT *, '       EOF readched'
        ENDIF
      ELSE
        izsize = 1_iintegers ! just to indicate > 0
        IF (idbg_level >= 20) THEN
          PRINT *, '       Read a record with size ', izsize_ga
        ENDIF
      ENDIF

    CASE ('ncdf')

      CALL read_netcdf (nufile, ie_in_tot, je_in_tot, izvar_count,           &
                        izlev_count, ivar_id, nvar_in, idims_id,             &
                        ndims_id, icomm_cart, my_cart_id, num_compute,       &
                        west_add_in, east_add_in, south_add_in, north_add_in,&
                        field_grib, myzvar, myzlev, myzlevtot, izsize,       &
                        imp_grib, lasync_io, yzerrmsg, izerror)

      IF (idbg_level >= 10) THEN
        IF (izsize > 0) THEN
          PRINT *, '       Got record ', myzvar, myzlev, ivar_id(myzvar),   &
                           var_in(myzvar)%name, izsize, field_grib(1)
        ELSE
          PRINT *, '       Got no more record: EOF '
        ENDIF
      ENDIF

    END SELECT

    IF (izerror /= 0) THEN
       CALL model_abort (my_cart_id, 2013, yzerrmsg, yzroutine)
    ENDIF

    IF (izsize == 0) THEN
      ! this PE has got no more record because the end of file is reached
      lzeof = .TRUE.
    ENDIF

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (6) = timings(6) + zrealdiff
    ENDIF

    IF (lzeof .AND. (izdebug > 10)) THEN
      PRINT *, '  End of File reached'
    ENDIF

  !-----------------------------------------------------------------------------
  ! 4.2: Unpack and convert the record
  !-----------------------------------------------------------------------------

    IF (.NOT. lzeof) THEN

      IF (yin_form_read == 'grb1') THEN
#ifdef GRIBDWD
        igrbednr = 1

        ! the field has to have the same precision as the REALs in the griblib
        field_grib = 0.0_irealgrib

        ! de-grib the data with the routine grbin1 from the grib library
        CALL grbin1 (idwdednr, undefgrib, ndims, idims_in, iblock, ibmap, ipds,  &
                     igds_in, ibms, ibds, dsup, field_grib, ierrf)

        IF (ierrf /= 0) THEN
          yzerrmsg = 'Error in grbin1'
          CALL model_abort (my_cart_id, 2013, yzerrmsg, yzroutine)
        ENDIF

        IF (idbg_level >= 20) THEN
          IF (izsize == 0) THEN
            PRINT *, '       EOF readched'
          ELSE
            PRINT *, '       Got record ', izsize, ipds(2), ipds(7), ipds(8), ipds(16)
          ENDIF
        ENDIF

        isec1len   = ipds( 1)
        ilocdef    = ipds(37)
        ilocensid  = ipds(50)
        iloctotnr  = ipds(51)
        SELECT CASE (iepstyp_bc)
        CASE (0,202,203)   !COSMO-EU boundary or "DWD-SREPS (lmi)" or "COSMO-LEPS (lml)"
          ! in this case iepsmem_bc is located in pds(52)
          ! (for the COSMO-SREPS it has been moved by dwdlib to pds(53)
          !  but this will not be the case for grib_api)
          ilocensnr  = ipds(52)
        CASE (201) !COSMO-SREPS boundary
          ilocensnr  = ipds(53)
        END SELECT

        ! check for member ID if necessary
        IF (leps_bc .AND. lcomp_bound) THEN
          IF (lchk_bc_typ) THEN
            ! lchk_bc_typ is only NL input
            SELECT CASE (iepstyp_bc)
!             CASE (0)   !COSMO-EU boundary
!             add now 202, 203 for the boundaries coming from "BC-EPS" and "COSMO-LEPS"
              CASE (0,202,203)   !COSMO-EU boundary or "BC-EPS (lmi)" or "COSMO-LEPS (lml)"
                IF (isec1len == 66 .AND. ilocdef == 253) THEN
                  ! in this case iepsmem_bc is located in pds(52)
                  ! (for the COSMO-SREPS it has been moved by dwdlib to pds(53)
                  !  but this will not be the case for grib_api)
                  IF ( ilocensid /= iepstyp_bc ) THEN
                    PRINT *,'iepstyp =',iepstyp_bc
                    PRINT *,'file with wrong eps typ:',yzpath(1:LEN_TRIM(yzpath))
                    PRINT *,'eps typ in grib:', ilocensid
                    yzerrmsg = 'ERROR   *** eps typ of grib file does not match iepstyp_bc ***'
                    CALL model_abort (my_cart_id, 2014, yzerrmsg, yzroutine)
                  ENDIF
                  IF ( ilocensnr /= iepsmem_bc ) THEN
                    PRINT *,'ieps_mem=',iepsmem_bc
                    PRINT *,'file with wrong member ID:',yzpath(1:LEN_TRIM(yzpath))
                    PRINT *,'member ID in grib:', ilocensnr
                    yzerrmsg = 'ERROR   *** Member ID of grib file does not match iepsmem_bc ***'
                    CALL model_abort (my_cart_id, 2014, yzerrmsg, yzroutine)
                  ENDIF
                ENDIF
              CASE (201) !COSMO-SREPS boundary
                IF ( ilocensnr /= iepsmem_bc ) THEN
                  PRINT *,'ieps_mem=',iepsmem_bc
                  PRINT *,'file with wrong member ID:',yzpath(1:LEN_TRIM(yzpath))
                  PRINT *,'member ID in grib:', ilocensnr
                  yzerrmsg = 'ERROR   *** Member ID of grib file does not match iepsmem_bc ***'
                  CALL model_abort (my_cart_id, 2014, yzerrmsg, yzroutine)
                ENDIF
              CASE DEFAULT !unknown boundary
                PRINT *,'boundary ensemble type BCTYP=',iepstyp_bc,'NOT IMPLEMENTED'
                PRINT *,'set lchk_bc_type=.FALSE. to IGNORE this check'
                yzerrmsg = 'ERROR   *** boundary ensemble type NOT IMPLEMENTED ***'
                CALL model_abort (my_cart_id, 2014, yzerrmsg, yzroutine)
            END SELECT
          ENDIF
        ENDIF
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

        ! get edition number
        CALL grib_get (igribid, 'editionNumber', igrbednr,   ireturn)

        ! Get some keys
        CALL grib_get (igribid, 'shortName',     yshortname, ireturn)

        IF ((igrbednr == 1) .AND. (yinput_model == 'IFS')) THEN
          ! The COSMO shortnames are not available for the IFS Grib1 coding, therefore
          ! we need the Grib1 style search here: get itabtyp, iee, ilevtyp, nzscan
          CALL grib_get (igribid, 'table2Version',          itabtyp,       ireturn)
          CALL grib_get (igribid, 'indicatorOfParameter',   iee,           ireturn)
          CALL grib_get (igribid, 'indicatorOfTypeOfLevel', ilevtyp,       ireturn)

          ! Scan through the variable list to get izloc and yshortname
          get_loc_grb_ifs: DO izloc = 1, nvar_in
            IF ((var_in(izloc)%ee     == iee    ) .AND.                 &
                (var_in(izloc)%tabtyp == itabtyp) .AND.                 &
                (var_in(izloc)%levtyp == ilevtyp)) THEN
              yshortname = TRIM(var_in(izloc)%name)
              EXIT get_loc_grb_ifs
            ENDIF
          ENDDO get_loc_grb_ifs

          IF (TRIM(yshortname) == 'unknown') THEN
            PRINT *,   ' *** Error: No variable found for IFS input ', iee, itabtyp, ilevtyp
            yzerrmsg = ' *** Error: No variable found for IFS input '
            CALL model_abort (my_cart_id, 2038, yzerrmsg, yzroutine)
          ENDIF
        ENDIF

        IF ((igrbednr == 2) .AND. (yinput_model == 'IFS')) THEN
          IF (TRIM(yshortname) == 'lnsp') THEN
            yshortname = 'LNPS   '
            ytypeoflevel = 'hybrid'
            ! print *, 'lnsp renamed to LNPS'
          ENDIF
          IF (TRIM(yshortname) == 'q') THEN
            yshortname = 'QV     '
            ! print *, 'q renamed to QV'
          ENDIF
          IF (TRIM(yshortname) == 't') THEN
            yshortname = 'T      '
            ! print *, 't renamed to T'
          ENDIF
          IF (TRIM(yshortname) == 'u') THEN
            yshortname = 'U      '
            ! print *, 'u renamed to U'
          ENDIF
          IF (TRIM(yshortname) == 'v') THEN
            yshortname = 'V      '
            ! print *, 'v renamed to V'
          ENDIF

          IF (TRIM(yshortname) == 'ciwc') THEN
            yshortname = 'QI     '
            ! print *, 'ciwc renamed to QI'
          ENDIF
          IF (TRIM(yshortname) == 'clwc') THEN
            yshortname = 'QC     '
            ! print *, 'clwc renamed to QC'
          ENDIF
          IF (TRIM(yshortname) == 'crwc') THEN
            yshortname = 'QR     '
            ! print *, 'crwc renamed to QR'
          ENDIF
          IF (TRIM(yshortname) == 'cswc') THEN
            yshortname = 'QS     '
            ! print *, 'cswc renamed to QS'
          ENDIF
        ENDIF

        IF ((igrbednr == 2) .AND. (yinput_model == 'GFS')) THEN
          IF (TRIM(yshortname) == 'GH' .OR. TRIM(yshortname) == 'gh') THEN
            yshortname = 'FI'
            ! print *, 'GH renamed to FI'
          ENDIF
          IF (TRIM(yshortname) == 'sp') THEN
            yshortname = 'PS'
            ! print *, 'GH renamed to FI'
          ENDIF
          IF (TRIM(yshortname) == 'q') THEN
            yshortname = 'QV'
            ! print *, 'CLWMR renamed to QC'
          ENDIF
          IF (TRIM(yshortname) == 't') THEN
            yshortname = 'T      '
            ! print *, 't renamed to T'
          ENDIF
          IF (TRIM(yshortname) == 'u') THEN
            yshortname = 'U      '
            ! print *, 'u renamed to U'
          ENDIF
          IF (TRIM(yshortname) == 'v') THEN
            yshortname = 'V      '
            ! print *, 'v renamed to V'
          ENDIF

          IF (TRIM(yshortname) == 'CLWMR' .OR. TRIM(yshortname) == 'clwmr') THEN
            yshortname = 'QC'
            ! print *, 'CLWMR renamed to QC'
          ENDIF
        ENDIF

        IF (TRIM(yshortname) == 'unknown') THEN
          ! Test for special variables
          ! Test whether it is GFS CLWMR, which is our QC, with Grib2 ids (0,1,22)
          ! or   whether it is GFS HGT, which we need for control geopotential
          CALL grib_get (igribid, 'discipline',        idisc, ireturn)
          CALL grib_get (igribid, 'parameterCategory', icatg, ireturn)
          CALL grib_get (igribid, 'parameterNumber',   ipara, ireturn)
          IF ( ireturn /= GRIB_SUCCESS ) THEN
            PRINT *, ' *** ERROR in grib_get: dis,cat,par: ', ireturn
          ENDIF
          IF     ((idisc==0) .AND. (icatg==1) .AND. (ipara==22)) THEN
            ! it is CLWMR = QC
            yshortname = 'QC'
          ELSEIF ((idisc==0) .AND. (icatg==1) .AND. (ipara==82)) THEN
            ! it is CLIMR = QI
            yshortname = 'QI'
          ELSEIF ((idisc==0) .AND. (icatg==3) .AND. (ipara==5)) THEN
            ! it is HGT: take the level in pcontrol_fi Pa
            yshortname = 'FI'

          ! Test for special UM variables
          ELSEIF ((idisc==0) .AND. (icatg==6) .AND. (ipara==23)) THEN
            ! it is QI from UM
            yshortname = 'QI'
          ELSEIF ((idisc==0) .AND. (icatg==0) .AND. (ipara==17)) THEN
            ! it is T_SKIN from UM, take it as T_SKIN
            yshortname = 'T_SKIN'
          ENDIF
        ELSEIF (TRIM(yshortname) == 'T_G') THEN
          ! for the COSMO-Model we work with T_S
          yshortname = 'T_S'
        ENDIF

        CALL grib_get (igribid, 'jScansPositively', jzscanmode,    ireturn)
        CALL grib_get (igribid, 'Ni',               igds_in(5),    ireturn)
        CALL grib_get (igribid, 'Nj',               igds_in(6),    ireturn)
        CALL grib_get (igribid, 'typeOfLevel',      ytypeoflevel,  ireturn)
        CALL grib_get (igribid, 'level',            ilevel,        ireturn)
        CALL grib_get (igribid, 'topLevel',         itoplevel,     ireturn)
        CALL grib_get (igribid, 'bottomLevel',      ibottomlevel,  ireturn)

        IF (ireturn /= GRIB_SUCCESS) THEN
          PRINT *,   ' *** Error in grib_api grib_get keys  ', ireturn
          yzerrmsg = ' *** Error in grib_api grib_get keys  '
          CALL model_abort (my_cart_id, 2036, yzerrmsg, yzroutine)
        ELSE
          IF (izdebug > 20) THEN
            PRINT *, ' got data: ', yshortname, igds_in(5), igds_in(6), ytypeoflevel, ilevel, itoplevel, ibottomlevel
          ENDIF
        ENDIF

        ! Get size of data and data itself
        ! Set missing_value before
        CALL grib_set (igribid, 'missingValue', undefgrib)
        CALL grib_get_size (igribid, 'values', igriblen, ireturn)
        IF (igriblen > ldsec) THEN
          PRINT *, ' *** ERROR: size of message is too big for allocated field: ', igriblen, ldsec
          yzerrmsg = ' *** ERROR: Wrong size of field ***'
          CALL model_abort (my_cart_id, 2037, yzerrmsg, yzroutine)
        ENDIF

        CALL grib_get (igribid, 'values', field_grib, ireturn)

        IF (ireturn /= GRIB_SUCCESS) THEN
          PRINT *,   ' *** Error in grib_api grib_get values', ireturn
          yzerrmsg = ' *** Error in grib_api grib_get values'
          CALL model_abort (my_cart_id, 2038, yzerrmsg, yzroutine)
        ENDIF

        CALL grib_get (igribid, 'section1Length',                    isec1len,    ireturn)
        CALL grib_get (igribid, 'localDefinitionNumber',             ilocdef ,    ireturn)
        CALL grib_get (igribid, 'localTypeOfEnsembleForecast',       ilocensid,   ireturn)
        CALL grib_get (igribid, 'numberOfForecastsInEnsemble',       iloctotnr,   ireturn)
        CALL grib_get (igribid, 'perturbationNumber',                ilocensnr,   ireturn)

        ! check for member ID if necessary
        IF (leps_bc .AND. lcomp_bound) THEN
          IF (lchk_bc_typ) THEN
            ! lchk_bc_typ is only NL input
            SELECT CASE (iepstyp_bc)
!             CASE (0)   !COSMO-EU boundary
!             add now 202, 203 for the boundaries coming from "BC-EPS" or "COSMO-LEPS"
              CASE (0,202,203)   !COSMO-EU boundary or "BC-EPS (lmi)" or "COSMO-LEPS"
                IF (isec1len == 66 .AND. ilocdef == 253) THEN
                  ! in this case iepsmem_bc is located in pds(52)
                  ! (for the COSMO-SREPS it has been moved by dwdlib to pds(53)
                  !  but this will not be the case for grib_api)
                  IF ( ilocensid /= iepstyp_bc ) THEN
                    PRINT *,'iepstyp =',iepstyp_bc
                    PRINT *,'file with wrong eps typ:',yzpath(1:LEN_TRIM(yzpath))
                    PRINT *,'eps typ in grib:', ilocensid
                    yzerrmsg = 'ERROR   *** eps typ of grib file does not match iepstyp_bc ***'
                    CALL model_abort (my_cart_id, 2014, yzerrmsg, yzroutine)
                  ENDIF
                  IF ( ilocensnr /= iepsmem_bc ) THEN
                    PRINT *,'ieps_mem=',iepsmem_bc
                    PRINT *,'file with wrong member ID:',yzpath(1:LEN_TRIM(yzpath))
                    PRINT *,'member ID in grib:', ilocensnr
                    yzerrmsg = 'ERROR   *** Member ID of grib file does not match iepsmem_bc ***'
                    CALL model_abort (my_cart_id, 2014, yzerrmsg, yzroutine)
                  ENDIF
                ENDIF
              CASE (201) !COSMO-SREPS boundary
                IF ( ilocensnr /= iepsmem_bc ) THEN
                  PRINT *,'ieps_mem=',iepsmem_bc
                  PRINT *,'file with wrong member ID:',yzpath(1:LEN_TRIM(yzpath))
                  PRINT *,'member ID in grib:', ilocensnr
                  yzerrmsg = 'ERROR   *** Member ID of grib file does not match iepsmem_bc ***'
                  CALL model_abort (my_cart_id, 2014, yzerrmsg, yzroutine)
                ENDIF
              CASE DEFAULT !unknown boundary
                PRINT *,'boundary ensemble type BCTYP=',iepstyp_bc,'NOT IMPLEMENTED'
                PRINT *,'set lchk_bc_type=.FALSE. to IGNORE this check'
                yzerrmsg = 'ERROR   *** boundary ensemble type NOT IMPLEMENTED ***'
                CALL model_abort (my_cart_id, 2014, yzerrmsg, yzroutine)
            END SELECT
          ENDIF
        ENDIF
#else
        yzerrmsg = ' ***  ERROR: model not compiled for grib_api *** '
        CALL model_abort (my_cart_id, 2040, yzerrmsg, yzroutine)
#endif

      ELSEIF (yin_form_read == 'ncdf') THEN
        igrbednr = 0

        ! set nunit_of_time depending on nincbound, in case Grib output is done
        ! nincbound = 1, 2, >= 4 corresponds to hincbound = 0.25, 0.5, >= 1.0
        IF     (nincbound == 1_iintegers) THEN
          nunit_of_time = 13
        ELSEIF (nincbound == 2_iintegers) THEN
          nunit_of_time = 14
! SP, 201405
!        ELSEIF (nincbound >= 4_iintegers) THEN
        ELSE
          nunit_of_time = 1
        ENDIF

      ENDIF

      ! convert field to the precision used in INT2LM
      IF (idbg_level >= 20) THEN
        PRINT *, '       Size of field  ', igds_in(5), igds_in(6), igriblen
      ENDIF

      IF (lec2lm .AND. (yin_form_read /= 'apix')) THEN
        jzscanmode = 0
      ENDIF

      ! UM horizontal staggered V has one grid point less in north-east than all other fields
      IF (lum2lm .AND. (igrbednr == 2)) THEN
        je_in_end = je_in_tot-1
      ELSE
        je_in_end = je_in_tot
      ENDIF

      DO j = 1 , je_in_end - north_add_in - south_add_in
        ! take into account the scanning mode
        IF (jzscanmode == 0) THEN
          jscan = je_in_tot+1 - j
        ELSE
          jscan = j
        ENDIF

        DO i = 1 , ie_in_tot - east_add_in - west_add_in
          ! this is without additional lines: ij = (j-1) * ie_in_tot + i
          ij = (j-1) * (ie_in_tot-east_add_in-west_add_in) + i
          IF (field_grib(ij) /= zundef) THEN
            field_glob(i+west_add_in,jscan+south_add_in) = REAL (field_grib(ij) , ireals)
          ELSE
            field_glob(i+west_add_in,jscan+south_add_in) = undef
          ENDIF
        ENDDO
      ENDDO

      IF (lum2lm .AND. (TRIM(yshortname) == 'V') .AND. (igrbednr == 2)) THEN
        ! also set the eastern border now for V
        DO i = 1, ie_in_tot
          field_glob(i,je_in_tot) = field_glob(i,je_in_tot-1)
        ENDDO
      ENDIF

      IF (lcm2lm) THEN
        ! Set the extra lines of the field, if present
        IF (west_add_in /= 0) THEN
          ! Set extra western boundary with values from the east
          DO j = 1, je_in_tot
            field_glob(1,j) = field_glob(ie_in_tot-east_add_in,j)
          ENDDO
        ENDIF
        IF (east_add_in /= 0) THEN
          ! Set extra eastern boundary with values from the west
          DO j = 1, je_in_tot
            field_glob(ie_in_tot,j) = field_glob(1+west_add_in,j)
          ENDDO
        ENDIF
        IF (south_add_in /= 0) THEN
          ! Set extra southern boundary with values from the second southern line
          DO i = 1, ie_in_tot
            field_glob(i,1) = field_glob(i,2)
          ENDDO
        ENDIF
        IF (north_add_in /= 0) THEN
          ! Set extra northern boundary with values from the second northern line
          DO i = 1, ie_in_tot
            field_glob(i,je_in_tot) = field_glob(i,je_in_tot-1)
          ENDDO
        ENDIF
      ENDIF

      ! convert LN(ps) to ps
      IF ( lec2lm ) THEN
        lzconvert=.FALSE.
        IF (iz_lnps_idx /= 0 .AND. iz_ps_idx /= 0) THEN
          IF (yin_form_read == 'grb1') THEN
            IF (ipds(7) == var_in(iz_lnps_idx)%ee .AND.                  &
                ipds(8) == var_in(iz_lnps_idx)%levtyp) THEN
              lzconvert = .TRUE.
              ipds(7) = var_in(iz_ps_idx)%ee
              ipds(8) = var_in(iz_ps_idx)%levtyp ; ipds(9) = 0 ; ipds(10) = 0
            ENDIF
          ELSEIF (yin_form_read == 'apix') THEN
            IF (yshortname   == 'LNPS' .AND.            &
                ytypeoflevel == 'hybrid') THEN
              lzconvert = .TRUE.
              yshortname = 'PS'
              ytypeoflevel = 'surface'
            ENDIF
          ENDIF
          IF (lzconvert) THEN
            field_glob_aux(:,:) = undef
            WHERE( field_glob .NE. undef ) field_glob_aux = EXP(field_glob)
            field_glob(:,:) = field_glob_aux(:,:)
          ENDIF
        ENDIF
      ENDIF

      !USUS this can be removed once grib_api recognizes other type of levels
      ! typeOfLevel is not recognized by grib_api, so we set it here
      IF ( (igrbednr == 2) .AND. lum2lm .AND. (ytypeOfLevel == 'unknown') ) THEN
        SELECT CASE (TRIM(yshortname))
        CASE ('U','V','P')
          ytypeoflevel = 'hybridLayer'
        CASE ('QV','QI','T')
          ytypeoflevel = 'hybrid'
        END SELECT
      ENDIF

      IF (yin_form_read == 'grb1') THEN

        ! Determine whether this record is required or not (set lzrequired)
        ! and put the data to the structure procarray.
        IF (lec2lm) THEN
          CALL check_required_ec (field_glob, ipds(7), ipds(8), ipds(9),       &
                      ipds(10), zprocarray, lzrequired, izloc, izlev, izrank,  &
                      idbg_level)
        ELSEIF (llm2lm) THEN
          CALL check_required_lm (field_glob,ipds(2),ipds(7),ipds(8),ipds(9),  &
                      ipds(10), zprocarray, lzrequired, izloc, izlev, izrank,  &
                      idbg_level)
        ELSEIF (lum2lm) THEN
          CALL check_required_um (field_glob,ipds(2),ipds(7),ipds(8),ipds(9),  &
                      ipds(10), zprocarray, lzrequired, izloc, izlev, izrank,  &
                      idbg_level)
        ENDIF

        IF (lzrequired) THEN
          yshortname = TRIM(var_in(izloc)%name)

          ! get unit of time range: this is distributed after the read-loop from
          ! PE 0 to all other PEs, because it could be that not all PEs really
          ! get a record for input. Note that it is assumed that all records in
          ! this grib file do have the same unit of time.
          nunit_of_time = INT (ipds(16), iintegers)

        ELSE
          yshortname = 'unknown'
        ENDIF

      ELSEIF (yin_form_read == 'apix') THEN

#ifdef GRIBAPI
        CALL check_required_api(field_glob, igribid, yshortname, ytypeoflevel, ilevel,  &
                    itoplevel, ibottomlevel, igrbednr, zprocarray,             &
                    lzrequired, izloc, izlev, izrank, idbg_level, izerror)
        IF (izerror /= 0_iintegers) THEN
          WRITE (yzerrmsg,'(A)') 'Errors in check_required_api'
          CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF
        IF (lzrequired) THEN
          ! and get the nunit_of_time, if required
          CALL grib_get (igribid, 'stepUnits',        nunit_of_time, ireturn)
        ENDIF
#else
        yzerrmsg = ' ***  ERROR: model not compiled for grib_api *** '
        CALL model_abort (my_cart_id, 2040, yzerrmsg, yzroutine)
#endif

      ELSEIF (yin_form_read == 'ncdf') THEN

        ! set lzrequired to .TRUE., because all records are needed
        lzrequired = .TRUE.

        ! set izloc, izrank and izlev
        izloc      = myzvar
        izrank     = var_in(izloc)%rank
        izlev      = myzlev
        yshortname = TRIM(var_in(izloc)%name)

        ! put values to zprocarray
        DO n=0,num_compute-1
          ie_p = isubpos_coarse(n,3) - isubpos_coarse(n,1) + 1
          je_p = isubpos_coarse(n,4) - isubpos_coarse(n,2) + 1
          DO j = 1, je_p
            DO i = 1, ie_p
              ij = (j-1)*ie_in_max + i
              zprocarray(ij,n)      =                                     &
               field_glob(isubpos_coarse(n,1)-1 + i,isubpos_coarse(n,2)-1 + j)
            ENDDO
          ENDDO
        ENDDO

      ENDIF

    ELSE
      IF (idbg_level > 10) THEN
        PRINT *, '  End of File reached'
      ENDIF
    ENDIF   ! lzeof

  !----------------------------------------------------------------------------
  ! 4.3: Get vertical coordinate parameters
  !----------------------------------------------------------------------------

    IF ((yin_form_read == 'grb1') .OR. (yin_form_read == 'apix')) THEN
      IF ((igrbednr == 2) .AND. (llm2lm)) THEN
        ! get meta data for generalVertical coordinate
        IF (.NOT. lzgot_genvc) THEN
          CALL get_generalvertical (yshortname, ytypeoflevel, igribid,      &
                                    (.NOT.lzeof) .AND. (lzrequired), izdebug)
          lzgot_genvc = .TRUE.
        ENDIF
      ELSE
        ! all other models have the same vertical coordinates and
        ! get_vert_coord must only be called until these parameters are read
        IF (lfirst .OR. lzfirst_uv) THEN
          IF (idbg_level > 10) THEN
            PRINT *, '  Get vertical coordinate parameters'
          ENDIF
          CALL get_vert_coord (lfirst, lzfirst_uv, yshortname, igribid,      &
                               igrbednr, (.NOT.lzeof) .AND. (lzrequired),    &
                               izdebug)
        ENDIF
      ENDIF
    ENDIF

    IF (yin_form_read == 'apix') THEN
      IF (lgsm2lm .AND. lfirst) THEN
        ! set the vertical coordinate parameters
        ak_in(51) =     0.000000_ireals;   bk_in(51) = 1.000000_ireals
        ak_in(50) =     0.000000_ireals;   bk_in(50) = 0.997000_ireals
        ak_in(49) =     0.000000_ireals;   bk_in(49) = 0.994000_ireals
        ak_in(48) =     0.000000_ireals;   bk_in(48) = 0.989000_ireals
        ak_in(47) =     0.000000_ireals;   bk_in(47) = 0.982000_ireals
        ak_in(46) =     0.000000_ireals;   bk_in(46) = 0.972000_ireals
        ak_in(45) =     0.000000_ireals;   bk_in(45) = 0.960000_ireals
        ak_in(44) =     0.000000_ireals;   bk_in(44) = 0.946000_ireals
        ak_in(43) =   133.051011_ireals;   bk_in(43) = 0.926669_ireals
        ak_in(42) =   364.904149_ireals;   bk_in(42) = 0.904351_ireals
        ak_in(41) =   634.602716_ireals;   bk_in(41) = 0.879654_ireals
        ak_in(40) =   959.797167_ireals;   bk_in(40) = 0.851402_ireals
        ak_in(39) =  1347.680042_ireals;   bk_in(39) = 0.819523_ireals
        ak_in(38) =  1790.907396_ireals;   bk_in(38) = 0.785091_ireals
        ak_in(37) =  2294.841690_ireals;   bk_in(37) = 0.748052_ireals
        ak_in(36) =  2847.484778_ireals;   bk_in(36) = 0.709525_ireals
        ak_in(35) =  3468.871488_ireals;   bk_in(35) = 0.668311_ireals
        ak_in(34) =  4162.956463_ireals;   bk_in(34) = 0.624370_ireals
        ak_in(33) =  4891.880833_ireals;   bk_in(33) = 0.580081_ireals
        ak_in(32) =  5671.824240_ireals;   bk_in(32) = 0.534282_ireals
        ak_in(31) =  6476.712996_ireals;   bk_in(31) = 0.488233_ireals
        ak_in(30) =  7297.469895_ireals;   bk_in(30) = 0.442025_ireals
        ak_in(29) =  8122.159791_ireals;   bk_in(29) = 0.395778_ireals
        ak_in(28) =  8914.082201_ireals;   bk_in(28) = 0.350859_ireals
        ak_in(27) =  9656.181911_ireals;   bk_in(27) = 0.307438_ireals
        ak_in(26) = 10329.436180_ireals;   bk_in(26) = 0.265706_ireals
        ak_in(25) = 10912.638440_ireals;   bk_in(25) = 0.225874_ireals
        ak_in(24) = 11369.647830_ireals;   bk_in(24) = 0.189304_ireals
        ak_in(23) = 11695.371600_ireals;   bk_in(23) = 0.155046_ireals
        ak_in(22) = 11861.253090_ireals;   bk_in(22) = 0.124387_ireals
        ak_in(21) = 11855.434320_ireals;   bk_in(21) = 0.096446_ireals
        ak_in(20) = 11663.355370_ireals;   bk_in(20) = 0.072366_ireals
        ak_in(19) = 11285.404060_ireals;   bk_in(19) = 0.052146_ireals
        ak_in(18) = 10729.949410_ireals;   bk_in(18) = 0.035701_ireals
        ak_in(17) = 10014.615050_ireals;   bk_in(17) = 0.022854_ireals
        ak_in(16) =  9167.247036_ireals;   bk_in(16) = 0.013328_ireals
        ak_in(15) =  8226.244908_ireals;   bk_in(15) = 0.006738_ireals
        ak_in(14) =  7201.568980_ireals;   bk_in(14) = 0.002484_ireals
        ak_in(13) =  6088.673009_ireals;   bk_in(13) = 0.000113_ireals
        ak_in(12) =  4950.000000_ireals;   bk_in(12) = 0.000000_ireals
        ak_in(11) =  4000.000000_ireals;   bk_in(11) = 0.000000_ireals
        ak_in(10) =  3230.000000_ireals;   bk_in(10) = 0.000000_ireals
        ak_in( 9) =  2610.000000_ireals;   bk_in( 9) = 0.000000_ireals
        ak_in( 8) =  2105.000000_ireals;   bk_in( 8) = 0.000000_ireals
        ak_in( 7) =  1700.000000_ireals;   bk_in( 7) = 0.000000_ireals
        ak_in( 6) =  1370.000000_ireals;   bk_in( 6) = 0.000000_ireals
        ak_in( 5) =  1105.000000_ireals;   bk_in( 5) = 0.000000_ireals
        ak_in( 4) =   893.000000_ireals;   bk_in( 4) = 0.000000_ireals
        ak_in( 3) =   720.000000_ireals;   bk_in( 3) = 0.000000_ireals
        ak_in( 2) =   581.000000_ireals;   bk_in( 2) = 0.000000_ireals
        ak_in( 1) =   469.000000_ireals;   bk_in( 1) = 0.000000_ireals

        DO k = 1, ke_in
          akh_in(k) = (ak_in(k) + ak_in(k+1)) * 0.5_ireals
          bkh_in(k) = (bk_in(k) + bk_in(k+1)) * 0.5_ireals
          dak_in(k) =  ak_in(k+1) - ak_in(k)
          dbk_in(k) =  bk_in(k+1) - bk_in(k)
        ENDDO

        lfirst = .FALSE.
      ENDIF
    ENDIF

  !----------------------------------------------------------------------------
  ! 4.4: Check the record
  !----------------------------------------------------------------------------

    IF (yin_form_read == 'grb1' .OR. yin_form_read == 'apix') THEN
      CALL check_input_grid (igrbednr, igds_in, ngds, ipds, npds, igribid,    &
               yin_form_read, yshortname, yzfulldate,                         &
               ie_in_tot, je_in_tot, ke_in+nlevskip, startlat_in_tot,         &
               startlon_in_tot, dlon_in, dlat_in, pollon_in, pollat_in,       &
               inrvert_in, pv_in, (.NOT.lzeof) .AND. (lzrequired),            &
               num_compute, icomm_cart, my_cart_id, .TRUE., itype_calendar,   &
               yinput_model, yzerrmsg, izerror)

      IF (izerror /= 0) THEN
        yzerrmsg = 'wrong grid description section or Namelist parameters'
        CALL model_abort (my_cart_id, 2004, yzerrmsg, yzroutine)
      ENDIF
    ENDIF

    ! print the maximum, minimum and meanvalues of each record
    IF (lchkin) THEN
      ! just to make this call save:
      IF (.NOT. lzrequired) izloc= 1
      CALL check_record (field_grib, 1, ie_in_tot, 1, je_in_tot, 1, 1,     &
               1, ie_in_tot, 1, je_in_tot, 1, 1, zundef,                   &
               yshortname(1:10)  , var_in(izloc)%ee, izlev,                &
               (.NOT.lzeof) .AND. (lzrequired), nuchkdat, num_compute,     &
               icomm_cart, my_cart_id, yzerrmsg, izerror)
    ENDIF

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (7) = timings(7) + zrealdiff
    ENDIF

  !----------------------------------------------------------------------------
  ! 4.5: Distribute record to all PEs and put values to memory
  !----------------------------------------------------------------------------

    distribute_loop: DO iproc = 0, num_compute-1

      ! The following routine handles the distribution of one record to all
      ! others. If the record is not required, the distribute_loop is cycled
      ! (return status izstat = -2), if all records are done, the read_loop
      ! is exited (return status izstat = -1).
      ! The sender is the processor with rank=iproc, all others receive the
      ! corresponding data into the structure zlocalarray.
      CALL scatter_data (zprocarray, iproc, izloc, izrank, izlev,    &
                         lzeof, lzrequired, izstat)

      IF (izstat == -1) EXIT  read_loop          ! all records are done
      IF (izstat == -2) CYCLE distribute_loop    ! record not required

    ENDDO distribute_loop

#ifdef GRIBAPI
    ! Release the grib_api handle
    CALL grib_release(igribid)
#endif

  ENDDO read_loop

!-------------------------------------------------------------------------------
! Section 5: Closing the file
!-------------------------------------------------------------------------------

  ! Deallocate the grib fields
  DEALLOCATE (iblock, ibmap, dsup)

  ! close file
  CALL close_file (nufile, yin_form_read, icomm_cart, my_cart_id,         &
                   num_compute, lasync_io, idbg_level, yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, izerror, yzerrmsg, 'close_file')
  ENDIF

  IF ( (lchkin) .AND. (my_cart_id == 0) ) THEN
    CLOSE (nuchkdat, STATUS='KEEP')
  ENDIF

  ! Distribute the unit of time
  IF (num_compute > 1) THEN
    ! Distribute the unit of time
    CALL distribute_values(nunit_of_time, 1, 0, imp_integers, icomm_cart, &
                           izerror)
  ENDIF

!-------------------------------------------------------------------------------
! Section 6: Calculation of reference atmosphere for COSMO-Model as input model
!-------------------------------------------------------------------------------

  IF (llm2lm .AND. lzrefatmos) THEN

    ! for input model grid

    ALLOCATE (zrho0 (ie_in, je_in, ke_in), zdp0 (ie_in, je_in, ke_in),      &
              zp0hl (ie_in, je_in, ke1in), zt0  (ie_in, je_in, ke_in),      &
              zt0hl (ie_in, je_in, ke1in),                                  &
              STAT=izerror)

    ! the splitting of the coarse grid topography hsurfs_in
    ! is done in the routine read_coarse_grid_ext
    IF (refatm_in%irefatm == 1) THEN
      CALL reference_atmosphere                                              &
       ( hhl_in, p0_in, zp0hl, zrho0, zt0, zt0hl, zdp0, hsurf_in, hsurfs_in, &
         ie_in, je_in, ke_in, refatm_in,                                     &
         vcoord_in, svc1_in, svc2_in, r_d, g, lanalyt_calc_t0p0, .TRUE.)
    ELSE IF (refatm_in%irefatm == 2) THEN
      CALL reference_atmosphere_2                                            &
       ( hhl_in, p0_in, zp0hl, zrho0, zt0, zt0hl, zdp0, hsurf_in, hsurfs_in, &
         ie_in, je_in, ke_in, refatm_in,                                     &
         vcoord_in, svc1_in, svc2_in, r_d, g, .TRUE.)
    ENDIF

    DEALLOCATE (zrho0, zdp0, zp0hl, zt0, zt0hl)

    ! for input model grid hhl interpolated to fine grid  (hhl_gl)
    ! NOTE: the hhl_gl field is constructed using a linearly interpolated
    !  version of the hsurfs_in field. It would not be correct to split the
    !  hsurf_gl field, since the filter action is defined in gridpoint space!

    ALLOCATE (zt0  (ie2lm,je2lm,ke_in),  zp0hl(ie2lm,je2lm,ke1in),           &
              zt0hl(ie2lm,je2lm,ke1in), STAT=izerror)

    IF (refatm_in%irefatm == 1) THEN
      CALL reference_atmosphere                                                &
       ( hhl_gl, p0_gl, zp0hl, rho0_gl, zt0, zt0hl, dp0_gl, hhl_gl(:,:,ke1in), &
         hsurfs_gl, ie2lm, je2lm, ke_in, refatm_in,                            &
         vcoord_in, svc1_in, svc2_in, r_d, g, lanalyt_calc_t0p0, .TRUE.)
    ELSE IF (refatm_in%irefatm == 2) THEN
      CALL reference_atmosphere_2                                              &
       ( hhl_gl, p0_gl, zp0hl, rho0_gl, zt0, zt0hl, dp0_gl, hhl_gl(:,:,ke1in), &
         hsurfs_gl, ie2lm, je2lm, ke_in, refatm_in,                            &
         vcoord_in, svc1_in, svc2_in, r_d, g, .TRUE.)
    ENDIF

    DEALLOCATE (zt0, zp0hl, zt0hl)
  ENDIF

!-------------------------------------------------------------------------------
! Section 7: Calculation of reference pressure for Unified Model as input model
!-------------------------------------------------------------------------------

  IF (lum2lm .OR. (lcm2lm .AND. lcm_hgt_coor)) THEN
    ! The reference pressure for the full levels of the Unified Model
    ! (interpolated to the COSMO horizontal grid) is computed according to
    ! ivctype_in=2 and irefatm_in=2 for the COSMO model
    ! (irefatm_in=1 does not work for new UM version with levels up to 37 km
    !  in the atmosphere; therefore only irefatm_in=2 is used. For the
    !  additional variables necessary the following defaults are used):
    refatm_in%irefatm =   2
    refatm_in%delta_t =    75.0_ireals
    refatm_in%h_scal  = 10000.0_ireals
    zt00       = refatm_in%t0sl - refatm_in%delta_t

    IF (lum2lm) THEN
      ! According to UM docu, P is on rho-levels, as are U, V
      DO  k = 1, ke_in
        p0_gl(:,:,k) = refatm_in%p0sl * EXP ( - g/r_d*refatm_in%h_scal/zt00 *           &
           LOG( ( EXP((akh_in_rho(k) + bkh_in_rho(k) * hsurf_gl(:,:))/refatm_in%h_scal) &
                    *zt00 + refatm_in%delta_t) / (zt00 + refatm_in%delta_t) ))
      ENDDO
    ELSEIF (lcm2lm .AND. lcm_hgt_coor) THEN
      ! But for the clm2lm, P is on the usual levels
      DO  k = 1, ke_in
        p0_gl(:,:,k) = refatm_in%p0sl * EXP ( - g/r_d*refatm_in%h_scal/zt00 *    &
           LOG( ( EXP((akh_in(k) + bkh_in(k) * hsurf_gl(:,:))/refatm_in%h_scal)  &
                    *zt00 + refatm_in%delta_t) / (zt00 + refatm_in%delta_t) ))
      ENDDO
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Section 8: Additional calculations
!-------------------------------------------------------------------------------

  ! If control level of geopotential could not be read, it is calculated here
  IF ( ((var_in(mzfi_loc_in)%dattyp(1:1) == 'I' .AND. .NOT. lcomp_bound)   .OR.  &
        (var_in(mzfi_loc_in)%dattyp(2:2) == 'B' .AND.       lcomp_bound) ) .AND. &
       ( (.NOT. var_in(mzfi_loc_in)%lreadin)) .AND.                              &
         (.NOT. lcm_hgt_coor) .AND. (.NOT. lcm_pres_coor)                  ) THEN

    IF (izdebug >= 20) THEN
      PRINT *, '      Compute control geopotential: ', ke_in
    ENDIF

    ALLOCATE (ztv (ie_in,je_in), zpo (ie_in,je_in), zpu(ie_in,je_in),   &
              zfiu(ie_in,je_in), zfio(ie_in,je_in), STAT=izstat)

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

  ! Additional calculations for a climate model
  IF (lcm2lm) THEN

    ! If luse_t_skin==.TRUE. use T_SKIN for T_S
    IF (var_in(mzts_loc_in)%lreadin) THEN
     IF (var_in(mztskin_loc_in)%lreadin .AND. luse_t_skin) THEN
        IF (my_cart_id == 0) PRINT *, ' *** T_SKIN is used for T_S ***'
          var_in(mzts_loc_in)%p2(:,:) = var_in(mztskin_loc_in)%p2(:,:)
      ENDIF

    ! If T_S is not given, use T_SKIN for T_S or T(ke)
    ELSE
      IF (var_in(mztskin_loc_in)%lreadin .AND. luse_t_skin) THEN
        IF (my_cart_id == 0) PRINT *, ' *** T_SKIN is used for T_S ***'
        var_in(mzts_loc_in)%p2(:,:) = var_in(mztskin_loc_in)%p2(:,:)
      ELSE
         IF (my_cart_id == 0) PRINT *, ' *** T(ke) is used for T_S ***'
       var_in(mzts_loc_in)%p2(:,:) = var_in(mzt_loc_in)%p3(:,:,ke_in)
      ENDIF
      var_in(mzts_loc_in)%lreadin = .TRUE.
    ENDIF

  ENDIF ! lcm2lm

  IF (vcoord_in%ivctype == 1) THEN ! beware sigm are no SIGMA but pressure heights
    CALL k_index_of_pressure_levels (refatm_in%p0sl, vcoord_in%sigm_coord/refatm_in%p0sl, &
                                     ke_in, .FALSE., klv950, klv850        )
  ELSEIF (vcoord_in%ivctype == 2 .OR. vcoord_in%ivctype == 3) THEN
    CALL k_index_of_pressure_levels (refatm_in%p0sl, vcoord_in%sigm_coord, &
                                     ke_in, .FALSE., klv950, klv850        )
  ELSE
    ! gdm klv850_in begin
    IF (lcm2lm) THEN
      IF (my_cart_id == 0) THEN
        PRINT *, 'vcoord_in%ivctype is not defined so use akh and bkh to define klv850'
      ENDIF
      ! use definition from here:
      ! http://www.ecmwf.int/products/data/archive/descriptions/od/oper/fc/ml/index.html
      DO k = 1,ke_in-1
        zpno = akh_in(k)   + bkh_in(k)  *101325.0_ireals
        zpnu = akh_in(k+1) + bkh_in(k+1)*101325.0_ireals
        IF ( (zpno <= 850.0E2) .AND. (850.0E2 < zpnu) ) klv850 = k
      ENDDO
    ELSE
      PRINT *, 'vcoord_in%ivctype is =', vcoord_in%ivctype
      PRINT *, 'so NOT 1,2 or 3 and I do NOT know how to compute klv850 !'
      yzerrmsg = 'Can not comput klv850 !'
      CALL model_abort (my_cart_id, 2024, yzerrmsg, yzroutine)
    ENDIF
    ! gdm klv850_in end
  ENDIF
  klv850_in = klv850

  IF (my_cart_id == 0) THEN
    PRINT *, ' vcoord_in%ivctype : ', vcoord_in%ivctype
    PRINT *, 'Computation of input level corresponding to 850 hPA is ', klv850_in
  ENDIF

  ! Additional calculations for NetCDF
  IF (yin_form_read == 'ncdf') THEN
    DO n = 1, nvar_in
       IF (.NOT. var_in(n)%lreadin) CYCLE !_br 18.12.2006
       IF (TRIM(var_in(n)%name) == 'T_SO') THEN
         var_in(n)%p3(1:ie_in,1:je_in,0) = var_in(mzts_loc_in)%p2(1:ie_in,1:je_in)
       ENDIF
       IF (TRIM(var_in(n)%name) == 'W_SO_REL') THEN
         WHERE (var_in(n)%p3(1:ie_in,1:je_in,:) == undefncdf) &
             var_in(n)%p3(1:ie_in,1:je_in,:) = 0._ireals
       ENDIF
       IF (TRIM(var_in(n)%name) == 'T_SNOW') THEN
         WHERE (var_in(n)%p2(1:ie_in,1:je_in) == undefncdf) &
            var_in(n)%p2(1:ie_in,1:je_in) = var_in(mzts_loc_in)%p2(1:ie_in,1:je_in)
       ENDIF
       IF (TRIM(var_in(n)%name) == 'W_SNOW'  .OR. TRIM(var_in(n)%name) == 'W_I'     .OR. &
           TRIM(var_in(n)%name) == 'Z0'      .OR. TRIM(var_in(n)%name) == 'LAI_MX'  .OR. &
           TRIM(var_in(n)%name) == 'LAI_MN'  .OR. TRIM(var_in(n)%name) == 'PLCOV_MX'.OR. &
           TRIM(var_in(n)%name) == 'PLCOV_MN'.OR. TRIM(var_in(n)%name) == 'ROOTDP'  .OR. &
           TRIM(var_in(n)%name) == 'T_CL'    .OR. TRIM(var_in(n)%name) == 'FRESHSNW') THEN
         WHERE (var_in(n)%p2(1:ie_in,1:je_in) == undefncdf) &
            var_in(n)%p2(1:ie_in,1:je_in) = 0._ireals
       ENDIF
    ENDDO
  ENDIF

  IF (llm2lm) THEN
    ! PIK U. Boehm - 15.12.06
    ! in case P is read use it to compute PP
    IF (var_in(mzp_loc_in)%lreadin .AND. .NOT. var_in(mzpp_loc_in)%lreadin) THEN
      IF (my_cart_id == 0) PRINT *, ' *** PP is calculated from P ***'
      var_in(mzpp_loc_in)%lreadin = .TRUE.
      var_in(mzpp_loc_in)%p3(:,:,:) = var_in(mzp_loc_in)%p3(:,:,:) - p0_in(:,:,:)
    ENDIF

    ! in case PP is read set also the readin-flag for P to TRUE
    IF (var_in(mzpp_loc_in)%lreadin .AND. .NOT. var_in(mzp_loc_in)%lreadin) THEN
      var_in(mzp_loc_in)%lreadin = .TRUE.
      ! P is not computed, because we do not need it
      var_in(mzp_loc_in)%p3(:,:,:) = var_in(mzpp_loc_in)%p3(:,:,:) + p0_in(:,:,:)
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Section 9: Check that all data are read
!-------------------------------------------------------------------------------

  ! Check, whether all data necessary are present and write information to
  ! output
  IF (my_cart_id == 0) THEN
    lzcheck = .TRUE.
    lzwarn  = .FALSE.
    WRITE (noutput, '(A)')           &
     '     Variable    has been read;    number of        from total'
    WRITE (noutput, '(A)')           &
     '                                  levels read    number of levels'

    DO n = 1, nvar_in
      lzreadin     = var_in(n)%lreadin
      yzdattyp     = var_in(n)%dattyp
      yzname(1:10) = var_in(n)%name(1:10)
      nzlevrea     = var_in(n)%nlevels_read
      nzlevels     = var_in(n)%nlevels

      ! Check only necessary data
      IF ( ((yzdattyp(1:1) == 'I' .AND. .NOT. lcomp_bound) .OR.  &
            (yzdattyp(2:2) == 'B' .AND.       lcomp_bound)) .AND. &
             yzdattyp(3:3) /= 'O' .OR. TRIM(yzname) == 'T_S') THEN

        IF      (llm2lm .AND. .NOT. lzreadin) THEN
          IF (n > nvar_in_norm) THEN
            ! For l_art, only give warnings, if data is not available
            lzwarn  = .TRUE.
          ELSE
            lzcheck = .FALSE.
            lzwarn  = .FALSE.
          ENDIF
        ELSE IF (lec2lm .AND. &
          TRIM(yzname) /= 'QI'     .AND. TRIM(yzname) /= 'QC'     .AND. &
          TRIM(yzname) /= 'T_SNOW' .AND. TRIM(yzname) /= 'QV_S'   .AND. &
          TRIM(yzname) /= 'DPSDT'  .AND. TRIM(yzname) /= 'LNPS'   .AND. &
         (TRIM(yzname) /= 'T_SKIN' .OR.   luse_t_skin) .AND. &
         .NOT. lzreadin) THEN
          lzcheck = .FALSE.
        ELSE IF (lcm2lm .AND. .NOT. lzreadin) THEN 
          lzcheck = .FALSE. 
        ENDIF

        IF (lzreadin) THEN
          WRITE (noutput, '(A13,  L14,A4,I9,I17)')                    &
            yzname(1:10), lzreadin,'   ;', nzlevrea, nzlevels
        ELSE
          WRITE (noutput, '(A3,A10,L14,A4,I9,I17)') '***',            &
            yzname(1:10), lzreadin,'   ;', nzlevrea, nzlevels
        ENDIF
      ENDIF
    ENDDO

    WRITE (noutput, '(A)') '         '

    IF (lzwarn) THEN
      PRINT *, ' *** WARNING: Not all chemistry fields are read ***'
      PRINT *, ' ***          See OUTPUT for more information   ***'
    ENDIF

    IF (.NOT. lzcheck) THEN
      ! Abort program
      yzerrmsg  = ' *** ERROR:  Not all data available ***'
      CALL model_abort (my_cart_id, 5051, yzerrmsg, 'org_read_coarse_grid')
    ENDIF

  ENDIF

  ! Is there any input-model that gives negative humidities?
  DO n = 1, nvar_in
    IF ( (var_in(n)%name=='QV' .OR. var_in(n)%name=='QC') .OR.    &
                    (lprog_qi .AND. var_in(n)%name=='QI') ) THEN
      DO k = 1, ke_in
        DO j = 1, je_in
          DO i = 1, ie_in
            IF (var_in(n)%p3(i,j,k) < 0.0_ireals) THEN
                var_in(n)%p3(i,j,k) = 0.0_ireals
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE org_read_coarse_grid

!==============================================================================
!+ Checks whether a record is required from LM data
!==============================================================================

SUBROUTINE check_required_lm (field_in, itabtyp,iee,ilevtyp,ilevtop,ilevbot, &
                              field_out, lrequire, iloc, ilev, irank, idebug)

!==============================================================================
!
! Description:
!   A Grib-Record with given values for the product definition section is
!   checked whether it is required for processing.
!
! Method:
!   The variable table is searched for the position of this record.
!   If it is not found, the record is not needed.
!
!==============================================================================
!
! Parameterlist:
REAL     (KIND=ireals)  , INTENT(IN)   ::  &
  field_in (ie_in_tot, je_in_tot)   ! input with record for the total domain

INTEGER  (KIND=intgribf), INTENT(IN)   ::  &
  iee,      & ! element number
  itabtyp,  & ! table type of the record
  ilevtyp,  & ! level type of the record
  ilevtop,  & ! number of top level
  ilevbot     ! number of bottom level

INTEGER (KIND=iintegers), INTENT(IN)   ::  &
  idebug      ! for debug output

REAL     (KIND=ireals)  , INTENT(OUT)  ::  &
  field_out (ie_in_max,je_in_max, 0:num_compute-1)
     ! output field decomposed for distributing to records

LOGICAL                 , INTENT(OUT)  ::  &
  lrequire    ! to indicate a required record

INTEGER (KIND=iintegers), INTENT(OUT)  ::  &
  iloc,     & ! location in the input variable table name
  ilev,     & ! for mlf: level in the atmosphere
  irank       ! rang of the variable

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers)  ::  &
  mzloc, mzlev, i, j, ii, ie_p, je_p, iprocs

REAL     (KIND=ireals)    ::  &
  zbias, zfactor
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  lrequire    = .FALSE.
  iloc        = 0_iintegers
  ilev        = 0_iintegers
  irank       = 0_iintegers

!------------------------------------------------------------------------------
! Section 2: Look up in the variable table, which variable it is
!------------------------------------------------------------------------------

  ! This search has to be different for levtyp = 109, 110, 111,112:
  ! for levtyp = 109 (multi level)
  ! for levtyp = 110 (multi level) levbot and levtop must not be compared
  ! for levtyp = 111,112 it is necessary to compare levbot and levtop

  mzloc = 0
  SELECT CASE (ilevtyp)
  CASE (109)

    DO i=1,nvar_in
      IF ((var_in(i)%tabtyp == itabtyp) .AND.  &
          (var_in(i)%levtyp == ilevtyp) .AND.  &
          (var_in(i)%ee     == iee    ) ) THEN
        mzloc = i
        mzlev = ilevbot
        EXIT
      ENDIF
    ENDDO

  CASE (1, 110)

    DO i=1,nvar_in
      IF ((var_in(i)%tabtyp == itabtyp) .AND.  &
          (var_in(i)%levtyp == ilevtyp) .AND.  &
          (var_in(i)%ee     == iee    ) ) THEN
        mzloc = i
        mzlev = ilevtop
        EXIT
      ENDIF
    ENDDO

  CASE (111, 112)

    IF ( (itabtyp==201) .AND. ( (iee==197) .OR. (iee==198) ) ) THEN
      ! these are the multi-dimensional variables for the multi-layer soil
      ! model: these are needed for lmulti_layer_in = .TRUE.
      ! for these variables, levtop has to be compared to msoilgrib_in
      IF (lmulti_layer_in) THEN
        IF (.NOT. lcomp_bound) THEN
          DO i=1,nvar_in
            IF ((var_in(i)%tabtyp        .EQ. itabtyp) .AND.  &
                (var_in(i)%levtyp        .EQ. ilevtyp) .AND.  &
                (var_in(i)%ee            .EQ. iee    ) ) THEN
              mzloc = i
              ! Determine mzlev
              DO ii=0,ke_soil_in+1
                IF (msoilgrib_in(ii) == ilevbot) mzlev = ii
              ENDDO
              EXIT
            ENDIF
          ENDDO
        ELSE
          IF (iee==197 .AND. ilevbot==0) THEN
            ! only T_SO(0) is needed (as former T_S)
            DO i=1,nvar_in
              IF ((var_in(i)%tabtyp        .EQ. itabtyp) .AND.  &
                  (var_in(i)%levtyp        .EQ. ilevtyp) .AND.  &
                  (var_in(i)%ee            .EQ. iee    ) ) THEN
                mzloc = i
                mzlev = 0
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
          IF ((var_in(i)%tabtyp        .EQ. itabtyp) .AND.  &
              (var_in(i)%levtyp        .EQ. ilevtyp) .AND.  &
              (var_in(i)%ee            .EQ. iee    ) .AND.  &
              (var_in(i)%levtop        .EQ. ilevtop) .AND.  &
              (var_in(i)%levbot        .EQ. ilevbot) ) THEN
            mzloc = i
            mzlev = 1
            EXIT
          ENDIF
        ENDDO
      ELSE
        ! in this case only T_S (ee=85,lvtyp=111,levbot=levtop=0) is needed
        IF (itabtyp==2 .AND. iee==85 .AND. ilevtyp==111 .AND.     &
                                   ilevtop==0 .AND. ilevbot==0) THEN
          DO i=1,nvar_in
            IF ((var_in(i)%tabtyp        .EQ. itabtyp) .AND.  &
                (var_in(i)%levtyp        .EQ. ilevtyp) .AND.  &
                (var_in(i)%ee            .EQ. iee    ) .AND.  &
                (var_in(i)%levtop        .EQ. ilevtop) .AND.  &
                (var_in(i)%levbot        .EQ. ilevbot) ) THEN
              mzloc = i
              mzlev = 1
              EXIT
            ENDIF
          ENDDO
        ELSE
          RETURN
        ENDIF
      ENDIF
    ENDIF

  CASE DEFAULT

    ! this record is not required
    RETURN

  END SELECT

  IF (idebug > 10) THEN
    IF (mzloc > 0) THEN
      WRITE (*,'(A,I3,A,A,5I6)') &
         'found:  ', mzloc, '   ', var_in(mzloc)%name, iee, ilevtyp,      &
                                           itabtyp, ilevtop, ilevbot
    ELSE
      WRITE (*,'(A,5I6)') &
         'no variable found in the table:   ', iee, ilevtyp, itabtyp,     &
                                           ilevtop, ilevbot
    ENDIF
  ENDIF

  IF (mzloc == 0) THEN
    ! This variable is not in the variable table (mzloc == 0)
    ! and hence not needed
    RETURN
  ENDIF

  IF (lcomp_bound) THEN
    IF (var_in(mzloc)%dattyp(2:2) /= 'B') THEN
      ! This variable is not needed for the boundary fields
      RETURN
    ENDIF
  ELSE
    IF (var_in(mzloc)%dattyp(1:1) /= 'I') THEN
      ! This variable is not needed for the initial fields
      RETURN
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Set the output variables
!------------------------------------------------------------------------------

  ! If this section is reached, it is a required record and mzloc, mzlev are
  ! already set
  lrequire = .TRUE.
  iloc     = mzloc
  ilev     = mzlev
  irank    = var_in(iloc)%rank

  ! Scale the field with bias and factor
  zbias   = var_in(iloc)%bias
  zfactor = var_in(iloc)%factor

  DO i=0,num_compute-1
    ie_p = isubpos_coarse(i,3) - isubpos_coarse(i,1) + 1
    je_p = isubpos_coarse(i,4) - isubpos_coarse(i,2) + 1
    field_out(1:ie_p,1:je_p,i)       =                                     &
        field_in(isubpos_coarse(i,1):isubpos_coarse(i,3),                  &
                 isubpos_coarse(i,2):isubpos_coarse(i,4)) / zfactor - zbias
  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE check_required_lm

!==============================================================================
!+ Checks whether a record is required from ECMWF data
!==============================================================================

SUBROUTINE check_required_ec (field_in, iee, ilevtyp, ilevtop, ilevbot,      &
                              field_out, lrequire, iloc, ilev, irank, idebug)

!==============================================================================
!
! Description:
!   A Grib-Record with given values for the product definition section is
!   checked whether it is required for processing.
!
! Method:
!   The variable table is searched for the position of this record.
!   If it is not found, the record is not needed.
!
!==============================================================================
!
! Parameterlist:
REAL     (KIND=ireals)  , INTENT(IN)   ::  &
  field_in (ie_in_tot, je_in_tot)   ! input with record for the total domain

INTEGER  (KIND=intgribf), INTENT(IN)   ::  &
  iee,      & ! element number
  ilevtyp,  & ! level type of the record
  ilevtop,  & ! number of top level
  ilevbot     ! number of bottom level

INTEGER (KIND=iintegers), INTENT(IN)   ::  &
  idebug      ! for debug output

REAL     (KIND=ireals)  , INTENT(OUT)  ::  &
  field_out (ie_in_max,je_in_max, 0:num_compute-1)
     ! output field decomposed for distributing to records

LOGICAL                 , INTENT(OUT)  ::  &
  lrequire    ! to indicate a required record

INTEGER (KIND=iintegers), INTENT(OUT)  ::  &
  iloc,     & ! location in the input variable table name
  ilev,     & ! for mlf: level in the atmosphere
  irank       ! rang of the variable

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers)  ::  &
  mzloc, mzlev, i, ie_p, je_p

REAL     (KIND=ireals)    ::  &
  zbias, zfactor
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  lrequire    = .FALSE.
  iloc        = 0_iintegers
  ilev        = 0_iintegers
  irank       = 0_iintegers

!------------------------------------------------------------------------------
! Section 2: Look up in the variable table, which variable it is
!------------------------------------------------------------------------------

  ! This search is only done for levtyp = 1, 109, 112
  ! but ecmwf has different ee-numbers for different levels, so it is not
  ! necessary to make the same distinction as for GME
  ! For pressure level data (control level geopotential) also the
  ! value of the level is checked

  mzloc = 0
  DO i=1,nvar_in
    IF (        (var_in(i)%levtyp == ilevtyp)                               &
         .AND.  (var_in(i)%ee     == iee)                                   &
         .AND. (     (ilevtyp == 100 .AND. var_in(i)%levbot == ilevbot)     &
                .OR. (ilevtyp == 109 .AND. ilevbot > nlevskip)              &
                .OR. (ilevtyp /= 100 .AND. ilevtyp /= 109)             )    ) THEN
      mzloc = i
      mzlev = ilevbot
      IF (ilevtyp == 109) mzlev = mzlev - nlevskip
      EXIT
    ENDIF
  ENDDO

  IF (idebug > 10) THEN
    IF (mzloc > 0) THEN
      WRITE (*,'(A,I3,A,A,4I6)') &
         'found:  ', mzloc, '   ', var_in(mzloc)%name, iee, ilevtyp,      &
                                           ilevtop, ilevbot
    ELSE
      WRITE (*,'(A,4I6)') &
         'no variable found in the table:   ', iee, ilevtyp,              &
                                           ilevtop, ilevbot
    ENDIF
  ENDIF

  IF (mzloc == 0) THEN
    ! This variable is not in the ECMWF variable table (mzloc == 0)
    ! and hence not needed
    RETURN
  ENDIF

  IF (lcomp_bound) THEN
    IF (var_in(mzloc)%dattyp(2:2) /= 'B') THEN
      ! This variable is not needed for the boundary fields
      RETURN
    ENDIF
  ELSE
    IF (var_in(mzloc)%dattyp(1:1) /= 'I') THEN
      ! This variable is not needed for the initial fields
      RETURN
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Set the output variables
!------------------------------------------------------------------------------

  ! If this section is reached, it is a required record and mzloc, mzlev are
  ! already set
  lrequire = .TRUE.
  iloc     = mzloc
  ilev     = mzlev
  irank    = var_in(iloc)%rank

  ! Scale the field with bias and factor
  zbias   = var_in(iloc)%bias
  zfactor = var_in(iloc)%factor

  DO i=0,num_compute-1
    ie_p = isubpos_coarse(i,3) - isubpos_coarse(i,1) + 1
    je_p = isubpos_coarse(i,4) - isubpos_coarse(i,2) + 1
    field_out(1:ie_p,1:je_p,i)       =                                     &
        field_in(isubpos_coarse(i,1):isubpos_coarse(i,3),                  &
                 isubpos_coarse(i,2):isubpos_coarse(i,4)) / zfactor - zbias
  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE check_required_ec

!==============================================================================
!+ Checks whether a record is required from UM data
!==============================================================================

SUBROUTINE check_required_um (field_in, itabtyp,iee,ilevtyp,ilevtop,ilevbot, &
                              field_out, lrequire, iloc, ilev, irank, idebug)

!------------------------------------------------------------------------------
!
! Description:
!   A Grib-Record with given values for the product definition section is
!   checked whether it is required for processing.
!
! Method:
!   The variable table is searched for the position of this record.
!   If it is not found, the record is not needed.
!
!------------------------------------------------------------------------------
!
! Parameterlist:
REAL     (KIND=ireals)  , INTENT(IN)   ::  &
  field_in (ie_in_tot, je_in_tot)   ! input with record for the total domain

INTEGER  (KIND=intgribf), INTENT(IN)   ::  &
  itabtyp,  & ! grib table used
  iee,      & ! element number
  ilevtyp,  & ! level type of the record
  ilevtop,  & ! number of top level
  ilevbot     ! number of bottom level

INTEGER (KIND=iintegers), INTENT(IN)   ::  &
  idebug      ! for debug output

REAL     (KIND=ireals)  , INTENT(OUT)  ::  &
  field_out (ie_in_max,je_in_max, 0:num_compute-1)
     ! output field decomposed for distributing to records

LOGICAL                 , INTENT(OUT)  ::  &
  lrequire    ! to indicate a required record

INTEGER (KIND=iintegers), INTENT(OUT)  ::  &
  iloc,     & ! location in the input variable table name
  ilev,     & ! for mlf: level in the atmosphere
  irank       ! rang of the variable

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers)  ::  &
  mzloc, mzlev, i, ie_p, je_p

REAL     (KIND=ireals)    ::  &
  zbias, zfactor
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  lrequire    = .FALSE.
  iloc        = 0_iintegers
  ilev        = 0_iintegers
  irank       = 0_iintegers

!------------------------------------------------------------------------------
! Section 2: Look up in the variable table, which variable it is
!------------------------------------------------------------------------------

  ! This search is only done for levtyp = 1, 120

  mzloc = 0

  DO i=1,nvar_in
    IF ((var_in(i)%tabtyp == itabtyp) .AND.  &
        (var_in(i)%levtyp == ilevtyp) .AND.  &
        (var_in(i)%ee     == iee    ) ) THEN
      mzloc = i
      mzlev = ilevbot
      EXIT
    ENDIF
  ENDDO


  IF (idebug > 10) THEN
    IF (mzloc > 0) THEN
      WRITE (*,'(A,I3,A,A,4I6)') &
         'found:  ', mzloc, '   ', var_in(mzloc)%name, iee, ilevtyp,      &
                                           ilevtop, ilevbot
    ELSE
      WRITE (*,'(A,4I6)') &
         'no variable found in the table:   ', iee, ilevtyp,              &
                                           ilevtop, ilevbot
    ENDIF
  ENDIF

  IF (mzloc == 0) THEN
    ! This variable is not in the UM variable table (mzloc == 0)
    ! and hence not needed
    RETURN
  ENDIF

  IF (lcomp_bound) THEN
    IF (var_in(mzloc)%dattyp(2:2) /= 'B') THEN
      ! This variable is not needed for the boundary fields
      RETURN
    ENDIF
  ELSE
    IF (var_in(mzloc)%dattyp(1:1) /= 'I') THEN
      ! This variable is not needed for the initial fields
      RETURN
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Set the output variables
!------------------------------------------------------------------------------

  ! If this section is reached, it is a required record and mzloc, mzlev are
  ! already set
  lrequire = .TRUE.
  iloc     = mzloc
  ilev     = mzlev
  irank    = var_in(iloc)%rank

  ! Scale the field with bias and factor
  zbias   = var_in(iloc)%bias
  zfactor = var_in(iloc)%factor

  DO i=0,num_compute-1
    ie_p = isubpos_coarse(i,3) - isubpos_coarse(i,1) + 1
    je_p = isubpos_coarse(i,4) - isubpos_coarse(i,2) + 1
    field_out(1:ie_p,1:je_p,i)       =                                     &
        field_in(isubpos_coarse(i,1):isubpos_coarse(i,3),                  &
                 isubpos_coarse(i,2):isubpos_coarse(i,4)) / zfactor - zbias
  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE check_required_um

!==============================================================================
#ifdef GRIBAPI
SUBROUTINE check_required_api                                               &
           (field_in, igribid, yshname, ytolevel, ilevin, itople, ibotle,   &
            iednr, field_out, lrequire, iloc, ilev, irank, idebug, ierr)

!------------------------------------------------------------------------------
!
! Description:
!   A grib_api-Record with given values for the product definition section is
!   checked whether it is required for processing.
!
! Method:
!   The variable table is searched for the position of this record.
!   If it is not found, the record is not needed.
!
!------------------------------------------------------------------------------
!
! Parameterlist:
REAL     (KIND=ireals)  , INTENT(INOUT)::  &
  field_in (ie_in_tot, je_in_tot)   ! input with record for the total domain

INTEGER (KIND=iintegers), INTENT(IN)   ::  &
  igribid     ! grib handle

CHARACTER (LEN=*),        INTENT(INOUT)::  &
  yshname     ! short name

CHARACTER (LEN=*),        INTENT(IN)   ::  &
  ytolevel    ! typeOfLevel

INTEGER  (KIND=intgribf), INTENT(IN)   ::  &
  ilevin,   & ! level from grib header
  itople,   & ! top level from grib header
  ibotle,   & ! bottom level from grib header
  iednr       ! grib edition number

INTEGER (KIND=iintegers), INTENT(IN)   ::  &
  idebug      ! for debug output

REAL     (KIND=ireals)  , INTENT(OUT)  ::  &
  field_out (ie_in_max,je_in_max, 0:num_compute-1)
     ! output field decomposed for distributing to records

LOGICAL                 , INTENT(OUT)  ::  &
  lrequire    ! to indicate a required record

INTEGER (KIND=iintegers), INTENT(OUT)  ::  &
  iloc,     & ! location in the input variable table name
  ilev,     & ! for mlf: level in the atmosphere
  irank       ! rang of the variable

INTEGER (KIND=iintegers), INTENT(INOUT) ::  &
  ierr        ! Error code

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers)  ::  &
  mzloc, mzlev, i, ii, ie_p, je_p, ipos, isearchlev, iscalval1, iscalval2,   &
  iscalfac1, iscalfac2, izednr, ireturn

REAL     (KIND=ireals)    ::  &
  zbias, zfactor, rfact, rdepth1, rdepth2
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  lrequire    = .FALSE.
  iloc        = 0_iintegers
  ilev        = 0_iintegers
  irank       = 0_iintegers
  ireturn     = 0_iintegers
  ierr        = 0_iintegers
  iscalval1   = -1
  iscalval2   = -1
  iscalfac1   = -1
  iscalfac2   = -1

!------------------------------------------------------------------------------
! Section 2: Look up in the variable table, which variable it is
!------------------------------------------------------------------------------

  mzloc = 0
  mzlev = 0

  IF (.NOT. lcomp_bound) THEN
    ipos = 1
  ELSE
    ipos = 2
  ENDIF

  ! Search for the location in the variable table
  ! For W_SO and W_SO_ICE, some special things have to be done, because they
  ! have different typeoflevels in Grib1 and Grib2:
  IF (yshname(1:4) == 'W_SO') THEN
    ! only search for the string of the shortname; the level type is
    ! determined later on
    DO i=1,nvar_in
      IF (  TRIM(var_in(i)%name)    == TRIM(yshname) ) THEN
        mzloc = i
        EXIT
      ENDIF
    ENDDO
  ELSE
    DO i=1,nvar_in
      IF ( (TRIM(var_in(i)%name)                  == TRIM(yshname)) .AND.    &
           (ylevltypes_in(var_in(i)%levtyp,iednr) == TRIM(ytolevel)) ) THEN
        mzloc = i
        EXIT
      ENDIF
    ENDDO
  ENDIF

! print *, 'Looking for ', yshname, ytolevel, ':  ', mzloc

  IF (idebug > 10) THEN
    IF (mzloc > 0) THEN
      WRITE (*,'(A,I3,A,A,A,3I6)') &
         'found:  ', mzloc, '   ', var_in(mzloc)%name, TRIM(ytolevel),    &
                                    ilevin, itople, ibotle
    ELSE
      WRITE (*,'(A,A,A,A,3I6)') &
         'no variable found in the table:   ', TRIM(yshname), '  ', TRIM(ytolevel), &
                                    ilevin, itople, ibotle
    ENDIF
  ENDIF

  IF (mzloc == 0) THEN
    ! This variable is not in the variable table (mzloc == 0)
    ! and hence not needed
    RETURN
  ENDIF

  ! Determine the level in the 3D array, if applicable
  SELECT CASE (TRIM(ytolevel))

  CASE ('surface')
    mzlev = 0

  CASE ('hybrid')
    IF (lum2lm) THEN
      mzlev = ke_in+2 - ilevin
    ELSE
      mzlev = ilevin - nlevskip
    ENDIF

  CASE ('hybridLayer')
    IF (lum2lm) THEN
      mzlev = ke_in+1 - itople
    ELSE
      mzlev = itople - nlevskip
    ENDIF

  CASE ('generalVertical')

    mzlev = ilevin - nlevskip

  CASE ('generalVerticalLayer')

    mzlev = itople - nlevskip

  CASE ('depthBelowLand')
    ! Variables: T_S, T_M, T_CL, T_SO    from GRIB1 and GRIB2
    !            W_SO                    from GRIB1

    IF     (iednr == 1) THEN

      ! In GRIB1, T_SO and W_SO are both coded with leveltype 111 (depthBelowLand);
      ! Normally, W_SO should be coded with leveltype 112 (depthBelowLandLayer),
      ! but the soil layers cannot properly be coded then
      ! The levels are coded in cm (integer values) in 'bottomlevel' (ipds(10))

      isearchlev = ibotle

      IF ((TRIM(yshname) == 'T_SO') .OR. (TRIM(yshname(1:4)) == 'W_SO')) THEN
        ! these are the multi-dimensional variables for the multi-layer soil
        ! model: these are needed for lmulti_layer_in = .TRUE.
        ! for these variables, levbot has to be compared to msoilgrib_in
        IF (lmulti_layer_in) THEN
          IF (.NOT. lcomp_bound) THEN
            DO i=1,nvar_in
              IF ( (TRIM(var_in(i)%name)                  == TRIM(yshname )) .AND. &
                   (ylevltypes_in(var_in(i)%levtyp,iednr) == TRIM(ytolevel)) .AND. &
                   (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
                mzloc = i
                ! Determine mzlev
                DO ii=0,ke_soil_in+1
                  IF (msoilgrib_in(ii) == isearchlev) mzlev = ii
                ENDDO
                EXIT
              ENDIF
            ENDDO
! print *, 'found       ', yshname, ytolevel, ':  ', mzloc, mzlev, isearchlev

          ELSE
            IF ((TRIM(yshname) == 'T_SO') .AND. (isearchlev == 0)) THEN
              ! only T_SO(0) is needed (as former T_S)
              DO i=1,nvar_in
                IF ( (TRIM(var_in(i)%name)                  == TRIM(yshname )) .AND. &
                     (ylevltypes_in(var_in(i)%levtyp,iednr) == TRIM(ytolevel)) .AND. &
                     (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
                  mzloc = i
                  mzlev = 0
                  EXIT
                ENDIF
              ENDDO
! print *, 'found       ', yshname, ytolevel, ':  ', mzloc, mzlev, isearchlev
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

! print *, 'looking for soil variables:  ', TRIM(yshname ), itople, ibotle
        ! We assume that this will be Grib1 only!! Grib2 will not work in this way
        IF (.NOT. lmulti_layer_in .OR. lec2lm) THEN
          DO i=1,nvar_in
             IF ( (TRIM(var_in(i)%name)                  == TRIM(yshname )) .AND. &
                  (ylevltypes_in(var_in(i)%levtyp,iednr) == TRIM(ytolevel)) .AND. &
                  (var_in(i)%levtop                      == itople)         .AND. &
                  (var_in(i)%levbot                      == ibotle)         .AND. &
                  (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mzloc = i
              mzlev = ibotle
              EXIT
            ENDIF
          ENDDO
        ELSE
          ! for the multi layer soil model these variables are not needed
          RETURN
        ENDIF
      ENDIF

    ELSEIF (iednr == 2) THEN
      ! should only be T_SO and using grib_api

      IF (lmulti_layer_in) THEN
        ! in the GRIB1 to GRIB2 conversion, the first level of T_SO(0) is converted to T_S!!
        ! This has to be reset
        IF (TRIM(yshname) == 'T_S') yshname = 'T_SO'
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
        RETURN
      END SELECT
      rdepth1 = iscalval1 * rfact

      IF (lmulti_layer_in) THEN
        IF (.NOT. lcomp_bound) THEN
          DO i=1,nvar_in
            IF ( (TRIM(var_in(i)%name)                  == TRIM(yshname )) .AND. &
                 (ylevltypes_in(var_in(i)%levtyp,iednr) == TRIM(ytolevel)) .AND. &
                 (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mzloc = i
              ! Determine mzlev
              IF (ABS(0.0_ireals - rdepth1)  < 1.0E-5_ireals ) THEN
                mzlev = 0
              ELSE
                DO ii=1,ke_soil_in+1
                  IF ( ABS(czmls_in(ii) - rdepth1) < 1.0E-5_ireals ) mzlev = ii
                ENDDO
                EXIT
              ENDIF
            ENDIF
          ENDDO
        ELSE
          IF ((TRIM(yshname) == 'T_SO') .AND. (iscalval1 == 0)) THEN
            ! only T_SO(0) is needed (as former T_S)
            DO i=1,nvar_in
              IF ( (TRIM(var_in(i)%name)                  == TRIM(yshname )) .AND. &
                   (ylevltypes_in(var_in(i)%levtyp,iednr) == TRIM(ytolevel)) .AND. &
                   (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
                mzloc = i
                mzlev = 0
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

    ENDIF ! iednr

  CASE ('depthBelowLandLayer')

    ! Variables: W_G1, W_G2, W_G3, W_CL (old soil model)  from GRIB1 (not available in GRIB2)
    !            W_SO                   (new soil model)  from GRIB2

    IF     (iednr == 1) THEN

      isearchlev = ibotle

      IF ( TRIM(yshname(1:4)) == 'W_SO' ) THEN
        ! these are the multi-dimensional variables for the multi-layer soil
        ! model: these are needed for lmulti_layer_in = .TRUE.
        ! for these variables, levtop has to be compared to msoilgrib_in
        IF (lmulti_layer_in) THEN
          DO i=1,nvar_in
            IF ( (TRIM(var_in(i)%name)                  == TRIM(yshname )) .AND. &
                 (ylevltypes_in(var_in(i)%levtyp,iednr) == TRIM(ytolevel)) .AND. &
                 (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mzloc = i
              ! Determine mzlev
              DO ii=0,ke_soil_in+1
                IF (msoilgrib_in(ii) == isearchlev) mzlev = ii
              ENDDO
              EXIT
            ENDIF
          ENDDO
        ELSE
          ! in this case these variables are not needed
          RETURN
        ENDIF
      ELSE
        ! these are the variables of the old soil model. These are only needed
        ! if lmulti_layer_in = .FALSE.

        ! We assume that this will be Grib1 only!! Grib2 will not work in this way
        IF (.NOT. lmulti_layer_in .OR. lec2lm) THEN
          DO i=1,nvar_in
             IF ( (TRIM(var_in(i)%name)                  == TRIM(yshname )) .AND. &
                  (ylevltypes_in(var_in(i)%levtyp,iednr) == TRIM(ytolevel)) .AND. &
                  (var_in(i)%levtop                      == itople)         .AND. &
                  (var_in(i)%levbot                      == ibotle)         .AND. &
                  (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mzloc = i
              mzlev = 1
              EXIT
            ENDIF
          ENDDO
        ELSE
          ! for the multi layer soil model these variables are not needed
          RETURN
        ENDIF
      ENDIF

    ELSEIF (iednr == 2) THEN
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
        RETURN
      END SELECT
      rdepth2 = iscalval2 * rfact

      IF (lmulti_layer_in) THEN
        IF (TRIM(yshname(1:4)) == 'W_SO') THEN
          DO i=1,nvar_in
            IF ( (TRIM(var_in(i)%name)    == TRIM(yshname )) .AND.    &
                 ('depthBelowLandLayer'   == TRIM(ytolevel)) .AND.    &
                 (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mzloc = i
              ! Determine mzlev
              DO ii=1,ke_soil_in+1
                IF ( (ABS(czhls_in(ii-1) - rdepth1) < 1.0E-5_ireals ) .AND.  &
                     (ABS(czhls_in(ii  ) - rdepth2) < 1.0E-5_ireals ) ) mzlev = ii
              ENDDO
              EXIT
            ENDIF
          ENDDO
        ELSE
          DO i=1,nvar_in
            IF ( (TRIM(var_in(i)%name)                  == TRIM(yshname )) .AND. &
                 (ylevltypes_in(var_in(i)%levtyp,iednr) == TRIM(ytolevel)) .AND. &
                 (var_in(i)%dattyp(ipos:ipos) /= ' ') ) THEN
              mzloc = i
              ! Determine mzlev
              DO ii=1,ke_soil_in+1
                IF ( (ABS(czhls_in(ii-1) - rdepth1) < 1.0E-5_ireals ) .AND.  &
                     (ABS(czhls_in(ii  ) - rdepth2) < 1.0E-5_ireals ) ) mzlev = ii
              ENDDO
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ELSE
        ! in this case these variables are not needed
        RETURN
      ENDIF

    ENDIF ! iednr

  CASE ('isobaricInhPa')

    SELECT CASE (TRIM(yshname))
    CASE ('FI')
      ! control geopotential
      mzlev = ilevin
    CASE DEFAULT
      ! GFS variables
      mzlev = get_gfs_level (ilevin)
    END SELECT

  CASE DEFAULT
    ! level typ not needed for INT2LM
    RETURN
  END SELECT


  ! Check for the control geopotential
  IF (TRIM(var_in(mzloc)%name) == 'FI') THEN
    ! only the values in pcontrol_fi Pa are needed
    IF (ABS(ilevin*100.0_ireals-pcontrol_fi) > 1E-2) THEN
      RETURN
    ENDIF
    IF (lgfs2lm) THEN
      !  this is something special for GFS
      ! the values are given as height, but we want to have it as geopotential
      field_in(:,:) = field_in(:,:) * g
    ENDIF
  ENDIF

  IF (lcomp_bound) THEN
    IF (var_in(mzloc)%dattyp(2:2) /= 'B') THEN
      ! This variable is not needed for the boundary fields
      RETURN
    ENDIF
  ELSE
    IF (var_in(mzloc)%dattyp(1:1) /= 'I') THEN
      ! This variable is not needed for the initial fields
      RETURN
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Set the output variables
!------------------------------------------------------------------------------

  ! If this section is reached, it is a required record and mzloc, mzlev are
  ! already set
  lrequire = .TRUE.
  iloc     = mzloc
  ilev     = mzlev
  irank    = var_in(iloc)%rank

  ! Scale the field with bias and factor
  zbias   = var_in(iloc)%bias
  zfactor = var_in(iloc)%factor

  DO i=0,num_compute-1
    ie_p = isubpos_coarse(i,3) - isubpos_coarse(i,1) + 1
    je_p = isubpos_coarse(i,4) - isubpos_coarse(i,2) + 1
    field_out(1:ie_p,1:je_p,i)       =                                     &
        field_in(isubpos_coarse(i,1):isubpos_coarse(i,3),                  &
                 isubpos_coarse(i,2):isubpos_coarse(i,4)) / zfactor - zbias
  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE check_required_api
#endif
!==============================================================================
!==============================================================================
!+ Extracts the vertical coordinate parameters for input model
!------------------------------------------------------------------------------

SUBROUTINE get_vert_coord (lgetmass, lgetwind, yshname, igribid,          &
                           igrbednr, lrequired, idebug)

!------------------------------------------------------------------------------
!
! Description:
!   PE 0 extracts the vertical coordinate parameters for input model from the
!   grid description section of a grib-record and distributes them to all
!   PEs.
!
! Method:
!
!------------------------------------------------------------------------------
!
! Parameterlist
LOGICAL,                INTENT(INOUT) ::  lgetmass, lgetwind

INTEGER(KIND=intgribf), INTENT(IN)    :: igribid
INTEGER(KIND=iintegers),INTENT(IN)    :: igrbednr, idebug
LOGICAL,                INTENT(IN)    :: lrequired

CHARACTER (LEN=*),      INTENT(IN)    :: yshname

! Local Variables
INTEGER  (KIND=iintegers)  ::  &
  izerror, n, k, k1, k2,       &
  izvert, izsender, izrecv (0:num_compute-1), idummy, ireturn

REAL (KIND=irealgrib)      :: refstf

! Compute zhhl_in for lec2lm
REAL (KIND=ireals)         :: zhhl_in (ke_in + 1)

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling

LOGICAL                    ::  &
  lzgotmass, lzgotwind

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine get_vert_coord
!------------------------------------------------------------------------------

  yzroutine = 'get_vert_coord'
  izerror   = 0
  ireturn   = 0
  izvert    = -1
  lzgotmass = .FALSE.
  lzgotwind = .FALSE.

  ! Check whether any PE has got vertical coordinates (igds_in(2) /= 0)
  IF (lrequired) THEN
    IF     (yin_form_read == 'grb1') THEN
      izvert = INT (igds_in(2), iintegers)
    ELSEIF (yin_form_read == 'apix') THEN
#ifdef GRIBAPI
      CALL grib_get (igribid, 'numberOfVerticalCoordinateValues', izvert, ireturn)
#endif
    ENDIF
  ENDIF

  IF (idebug > 10) THEN
    PRINT *, 'Number of vertical coordinate parameters:  ', izvert
  ENDIF

!US  IF ((TRIM(yshname) == 'U') .OR. (TRIM(yshname) == 'V')) THEN
!US changed that for um2lm
  IF ((TRIM(yshname) == 'U') .OR. (TRIM(yshname) == 'V') .OR. (TRIM(yshname) == 'P')) THEN
    IF (izvert > 0) THEN
      izvert = izvert + 1000
    ENDIF
  ENDIF

  IF (num_compute > 1) THEN
    CALL gather_values (izvert, izrecv, 1, num_compute, imp_integers, -1,    &
                        icomm_cart, yzerrmsg, izerror)
  ELSE
    izrecv(0) = izvert
  ENDIF

  izsender = -1
  izvert   = -1
  search: DO n = 0, num_compute-1
    IF (izrecv(n) > 0) THEN
      izsender = n
      izvert   = izrecv(n)
      EXIT search
    ENDIF
  ENDDO search

  IF (izvert > 0 .AND. izvert < 999) THEN
    lzgotmass = .TRUE.
  ELSEIF (izvert > 1000) THEN
    lzgotwind = .TRUE.
    izvert = izvert - 1000
  ENDIF

  IF (izsender > -1) THEN

    ! Check, whether pv_in is big enough
    IF (izvert > inrvert_in) THEN
      PRINT *, 'array for vertical coordinate parameters too small:  ', inrvert_in, izvert
      yzerrmsg = 'array for vertical coordinate parameters too small:  '
      CALL model_abort (my_cart_id, 2004, yzerrmsg, yzroutine)
    ENDIF

    ! PE with rank izsender got vertical coordinates, unpacks them (for dwdgrib1)
    ! and sends them to the others
    IF (my_cart_id == izsender) THEN
      IF     (yin_form_read == 'grb1') THEN
#ifdef GRIBDWD
        IF (llm2lm) THEN
          ! check for the old or new style of coding (and how many can be decoded)
          IF (REAL(REFSTF(igds_in(26)), ireals) > 500) THEN
            ! it is the old style of coding
! print *, 'OLD', REFSTF(igds_in(26)), izvert, REFSTF(igds_in(26))
            DO n = 1, ke1in+4
              pv_in(n) = REAL (REFSTF(igds_in(25+n)), ireals)
            ENDDO
            IF (izvert > ke1in+4) THEN
              ! additional values are coded
              IF (ABS(igds_in(25+ke1in+4+1)) > 1000) THEN
                pv_in(ke1in+4+1) = REAL(REFSTF(igds_in(25+ke1in+4+1)), ireals)  ! this is ivctype
              ELSE
                pv_in(ke1in+4+1) = REAL(       igds_in(25+ke1in+4+1) , ireals)  ! this is ivctype
              ENDIF
              DO n = ke1in+4+2, izvert
                IF (ABS(igds_in(25+n)) > 1000) THEN
                  pv_in(n) = REAL(REFSTF(igds_in(25+n)), ireals)
                ELSE
                  pv_in(n) = REAL(       igds_in(25+n) , ireals)
                ENDIF
              ENDDO
            ENDIF
          ELSE
            ! it is the new style of coding
! print *, 'NEW', REFSTF(igds_in(26))
            DO n = 1, ke1in+6
              pv_in(n) = REAL (REFSTF(igds_in(25+n)), ireals)
            ENDDO
            DO n = ke1in+6+1, izvert
              IF (igds_in(25+n) /= -1) THEN
                pv_in(n) = REAL(REFSTF(igds_in(25+n)), ireals)
              ELSE
                pv_in(n) = 0.0_ireals
              ENDIF
            ENDDO
          ENDIF
        ELSE  ! not llm2lm
          DO n = 1, izvert
            pv_in(n) = REAL (REFSTF(igds_in(25+n)), ireals)
          ENDDO
        ENDIF
#endif
      ELSEIF (yin_form_read == 'apix') THEN
#ifdef GRIBAPI
        CALL grib_get (igribid, 'pv', pv_in, ireturn)
#endif
      ENDIF
    ENDIF

    IF (num_compute > 1) THEN
      CALL distribute_values(pv_in, izvert, izsender, imp_reals, icomm_cart, izerror)
    ENDIF
    inrvert_in = izvert

    IF (idebug > 10) THEN
      PRINT *, 'Vertical coordinate parameters:  ', pv_in
    ENDIF

    ! Every PE now has the real values of the vertical coordinates
    ! and has to process these

    IF (llm2lm) THEN
      lgetmass = .FALSE.
      lgetwind = .FALSE.

      ! Check for the type of the vertical coordinate parameters
      ! (In the old style of coding, this value is p0sl, so rather big,
      !  or it is the new style of coding, then it is a small integer)
      idummy = NINT(pv_in(1), iintegers)

      IF ((idummy >= 1) .AND. (idummy <= 300)) THEN

        ! This is the new grib GDS coding style introduced with INT2LM 1.5
        IF      (idummy <= 100) THEN
          vcoord_in%ivctype = idummy
          refatm_in%irefatm = 1
        ELSEIF ((idummy > 100) .AND. (idummy <= 200)) THEN
          vcoord_in%ivctype = idummy - 100
          refatm_in%irefatm = 2
        ELSEIF ((idummy > 200) .AND. (idummy <= 300)) THEN
          vcoord_in%ivctype = idummy - 200
          refatm_in%irefatm = 3
        ELSE
          PRINT *,  ' ERROR *** Type ivctype_in of vertical coordinate '// &
                                                        'not available***'
          izerror = 10
        ENDIF

        refatm_in%p0sl   = pv_in( 3)
        refatm_in%t0sl   = pv_in( 4)
        refatm_in%dt0lp  = pv_in( 5)
        vcoord_in%vcflat = pv_in( 6)

        IF (vcoord_in%ivctype == 1) THEN
          DO k = 1, ke1in
            vcoord_in%sigm_coord(k)  = pv_in(6+k)
          ENDDO
        ELSE
          DO k = 1, ke1in
            vcoord_in%vert_coord(k)  = pv_in(6+k)
          ENDDO
        ENDIF

        IF (vcoord_in%ivctype == 3 .OR. vcoord_in%ivctype == 4) THEN
          ! read three more SLEVE parameters
          svc1_in   = pv_in(6 + ke1in + 1)
          svc2_in   = pv_in(6 + ke1in + 2)
          nfltvc_in = NINT(pv_in(6 + ke1in + 3), iintegers)
        ENDIF

        IF (refatm_in%irefatm == 2) THEN
          ! read additional parameters for new reference atmosphere
          refatm_in%delta_t = pv_in(6 + ke1in + 4)
          refatm_in%h_scal  = pv_in(6 + ke1in + 5)
        ENDIF

      ELSE
        ! This is the old grib GDS coding style, kept for backwards
        ! compatibility.
        refatm_in%p0sl   = pv_in( 1)
        refatm_in%t0sl   = pv_in( 2)
        refatm_in%dt0lp  = pv_in( 3)
        vcoord_in%vcflat = pv_in( 4)

        ! Check for the type of the vertical coordinate parameters
        ! again: this was at a different location then
        !US now we have pv_in: ivctype_in = INT(igds_in(29 + ke1in + 1), iintegers)
        !   but in veeeery old data, even this was not coded and has to be checked here!
        IF (izvert > ke1in+4) THEN
          vcoord_in%ivctype = INT (pv_in(4 + ke1in +1), iintegers)
        ELSE
          vcoord_in%ivctype = -999999      ! GDS default values
        ENDIF

        IF (vcoord_in%ivctype > 100) THEN
          vcoord_in%ivctype = vcoord_in%ivctype - 100
          refatm_in%irefatm = 2
        ELSE
          refatm_in%irefatm = 1
        ENDIF

        IF (vcoord_in%ivctype == -999999) THEN
          ! This is the old type of coding
          ! The vertical coordinate type has to be determined by the
          ! ascending or descending of the vertical coordinates
          IF ( pv_in(4+2) > pv_in(4+1) ) THEN
            vcoord_in%ivctype = 1
            DO k = 1, ke1in
              vcoord_in%sigm_coord(k)  = pv_in(4 + k)
            ENDDO
          ELSEIF ( pv_in(4+2) < pv_in(4+1) ) THEN
            vcoord_in%ivctype = 2
            DO k = 1, ke1in
              vcoord_in%vert_coord(k)  = pv_in(4 + k)
            ENDDO
          ELSE
            PRINT *,  ' ERROR *** Type ivctype_in of vertical coordinate '// &
                                                        'not available***'
            izerror = 10
          ENDIF
        ELSEIF (vcoord_in%ivctype == 1) THEN
          DO k = 1, ke1in
            vcoord_in%sigm_coord(k)  = pv_in(4 + k)
          ENDDO
        ELSEIF (vcoord_in%ivctype == 2) THEN
          DO k = 1, ke1in
            vcoord_in%vert_coord(k)  = pv_in(4 + k)
          ENDDO
        ELSEIF (vcoord_in%ivctype == 3 .OR. vcoord_in%ivctype == 4) THEN
          DO k = 1, ke1in
            vcoord_in%vert_coord(k)  = pv_in(4 + k)
          ENDDO
          ! and read three more SLEVE parameters
          svc1_in   = pv_in(4 + ke1in + 2)
          svc2_in   = pv_in(4 + ke1in + 3)
          nfltvc_in = NINT(pv_in(4 + ke1in + 4), iintegers)
        ENDIF

        IF (refatm_in%irefatm == 2) THEN
          ! read additional parameters for new reference atmosphere
          refatm_in%delta_t = pv_in(4 + ke1in + 5)
          refatm_in%h_scal  = pv_in(4 + ke1in + 6)
        ENDIF

      ENDIF


      ! Now all vertical coordinate parameters and ivctype are set
      ! go on with some checks

      IF     (vcoord_in%ivctype == 1) THEN
        ! For this type the vertical coordinates should be ascending
        IF ( vcoord_in%sigm_coord(2) < vcoord_in%sigm_coord(1) ) THEN
          PRINT *,                                                         &
            ' ERROR *** Vertical coordinates not ascending for type *** ', &
             vcoord_in%ivctype, vcoord_in%sigm_coord(1), vcoord_in%sigm_coord(2)
          izerror = 5
        ENDIF
      ELSEIF (vcoord_in%ivctype == 2) THEN
        ! For this type the vertical coordinates should be descending
        IF ( vcoord_in%vert_coord(2) > vcoord_in%vert_coord(1) ) THEN
          PRINT *,                                                         &
            ' ERROR *** Vertical coordinates not descending for type ***',&
             vcoord_in%ivctype, vcoord_in%vert_coord(1), vcoord_in%vert_coord(2)
          izerror = 6
        ENDIF
      ELSEIF (vcoord_in%ivctype == 3 .OR. vcoord_in%ivctype == 4) THEN

        ! Check for meaningful values of svc1, svc2 and nfltvc
        IF ((svc1_in > vcoord_in%vert_coord(1)) .OR. (svc1_in < 0.0)) THEN
          PRINT *, ' ERROR *** svc1_in not in allowed '//                &
                          'range for ivctype_in = 3 or 4 ***'
          izerror = 7
        ENDIF

        IF ((svc2_in > vcoord_in%vert_coord(1)) .OR. (svc2_in < 0.0)) THEN
          PRINT *,  ' ERROR *** svc2_in not in allowed '//                &
                           'range for ivctype_in = 3 or 4 ***'
          izerror = 8
        ENDIF

        IF (nfltvc_in <= 0) THEN
          PRINT *,  ' ERROR *** nfltvc_in must be greater than '//        &
                                          'or equal to zero ***'
          izerror = 9
        ENDIF
      ELSE
        PRINT *,  ' ERROR *** Type ivctype_in of vertical coordinate '// &
                                                    'not available***'
        izerror = 11
      ENDIF

    ELSEIF (lec2lm) THEN

      lgetmass = .FALSE.
      lgetwind = .FALSE.

      vcoord_in%ivctype = 1
      DO k = 1, ke_in + 1
        ak_in(k) = pv_in(          nlevskip +     k)
        bk_in(k) = pv_in(ke_in + 2*nlevskip + 1 + k)
        vcoord_in%sigm_coord(k) = ak_in(k) + bk_in(k)*refatm_in%p0sl
        IF (vcoord_in%sigm_coord(k) > 0.0_ireals) THEN
          zhhl_in(k) = (r_d/g)*LOG(refatm_in%p0sl/vcoord_in%sigm_coord(k)) &
                          * (refatm_in%t0sl - 0.5 * refatm_in%dt0lp        &
                          * LOG(refatm_in%p0sl/vcoord_in%sigm_coord(k)) )
        ELSE
          zhhl_in(k) = undef
        END IF
      ENDDO

      DO k = 1, ke_in
        akh_in(k) = (ak_in(k) + ak_in(k+1)) * 0.5_ireals
        bkh_in(k) = (bk_in(k) + bk_in(k+1)) * 0.5_ireals
        dak_in(k) =  ak_in(k+1) - ak_in(k)
        dbk_in(k) =  bk_in(k+1) - bk_in(k)
      ENDDO

      ! PE 0 writes the parameters to OUTPUT
      IF (my_cart_id == 0) THEN
        WRITE (noutput,'(A)') '     '
        WRITE (noutput,'(A)') '     Vertical coordinate parameters of input:'
        WRITE (noutput,'(A)') '     '
        WRITE (noutput,'(A)') '     k        ak_in(k) [Pa]   bk_in(k)    vcoord_in(k) [Pa]  zhhl_in(k) [m]'
        DO k = 1,ke_in+1
          WRITE (noutput,'(I6,4F16.4)') k, ak_in(k), bk_in(k), vcoord_in%sigm_coord(k), zhhl_in(k)
        ENDDO
        WRITE (noutput,'(A)') '     '
      ENDIF

    ELSEIF (lum2lm) THEN

      IF     (lgetwind .AND. lzgotwind) THEN
        lgetwind = .FALSE.

        IF (igrbednr == 1) THEN
          ! for grib1, the UM data are pre-processed (by pp2grib), and the counting of
          ! levels is from top to bottom (as in all DWD models)
          DO k = 1, ke_in
            akh_in_rho(k) = pv_in(            k)
            bkh_in_rho(k) = pv_in(ke_in + 1 + k)
          ENDDO
        ELSE
          ! for grib2, the original UM data are processed, which count from bottom to top
          ! therefore the order has to be reversed
          ! Careful: in grib2, ke_in+1 values are stored for ak, bk, but the last one is nonsense
          DO k = 1, ke_in
            akh_in_rho(k) = pv_in(  ke_in+1 - k)
            bkh_in_rho(k) = pv_in(2*ke_in+2 - k)
          ENDDO
        ENDIF

        ! PE 0 writes the parameters to OUTPUT
        IF (my_cart_id == 0) THEN
          WRITE (noutput,'(A)') '     '
          WRITE (noutput,'(A)') '     Vertical coordinate parameters for wind levels of UM:'
          WRITE (noutput,'(A)') '     '
          WRITE (noutput,'(A)') '     k        akh_in_rho(k) [m]   bkh_in_rho(k) [0..1] '
          DO k = 1,ke_in
            WRITE (noutput,'(I6,2F16.4)') k, akh_in_rho(k), bkh_in_rho(k)
          ENDDO
          WRITE (noutput,'(A)') '     '
        ENDIF

      ELSEIF (lgetmass .AND. lzgotmass) THEN

        lgetmass = .FALSE.

        IF (igrbednr == 1) THEN
          ! for grib1, the UM data are pre-processed (by pp2grib), and the counting of
          ! levels is from top to bottom (as in all DWD models)
          DO k = 1, ke_in
            akh_in(k) = pv_in(            k)
            bkh_in(k) = pv_in(ke_in + 1 + k)
          ENDDO
        ELSE
          ! for grib2, the original UM data are processed, which count from bottom to top
          ! therefore the order has to be reversed
          ! Careful: in grib2, ke_in+1 values are stored for ak, bk, but the first one is not used
          DO k = 1, ke_in
            akh_in(k) = pv_in(  ke_in+2 - k)
            bkh_in(k) = pv_in(2*ke_in+3 - k)
          ENDDO
        ENDIF

        DO k = 1, ke_in
          ak_in(k)  = HUGE(1.0_ireals)   ! (ak_in(k) + ak_in(k+1)) * 0.5_ireals
          bk_in(k)  = HUGE(1.0_ireals)   ! (bk_in(k) + bk_in(k+1)) * 0.5_ireals
          dak_in(k) = HUGE(1.0_ireals)   !  ak_in(k+1) - ak_in(k)
          dbk_in(k) = HUGE(1.0_ireals)   !  bk_in(k+1) - bk_in(k)
        ENDDO

        ! PE 0 writes the parameters to OUTPUT
        IF (my_cart_id == 0) THEN
          WRITE (noutput,'(A)') '     '
          WRITE (noutput,'(A)') '     Vertical coordinate parameters for main levels of UM:'
          WRITE (noutput,'(A)') '     '
          WRITE (noutput,'(A)') '     k        akh_in(k) [m]    bkh_in(k) [0..1]  '
          DO k = 1,ke_in
            WRITE (noutput,'(I6,2F16.4)') k, akh_in(k), bkh_in(k)
          ENDDO
          WRITE (noutput,'(A)') '     '
        ENDIF
      ENDIF

    ELSEIF (lhir2lm) THEN
      lgetmass = .FALSE.
      lgetwind = .FALSE.

      ! It seems that HIRLAM provides the vertical coordinate parameters of the
      ! full levels. The coordinate parameters of the half levels are computed
      ! as arithmetic means, setting (0,0) at the top and (0,1) at the surface
      DO k = 1, ke_in
        akh_in(k) = pv_in(2*k-1)
        bkh_in(k) = pv_in(2*k  )
      ENDDO

      vcoord_in%ivctype = 1
      ak_in(   1   ) = 0.0_ireals
      bk_in(   1   ) = 0.0_ireals
      vcoord_in%sigm_coord( 1 ) = 0.0_ireals
      ak_in(ke_in+1) = 0.0_ireals
      bk_in(ke_in+1) = 1.0_ireals
      vcoord_in%sigm_coord(ke_in+1) = refatm_in%p0sl

      DO k = 2, ke_in
        ak_in (k) = (akh_in(k-1) + akh_in(k)) * 0.5_ireals
        bk_in (k) = (bkh_in(k-1) + bkh_in(k)) * 0.5_ireals
        vcoord_in%sigm_coord(k) = ak_in(k) + bk_in(k)* refatm_in%p0sl
        IF (vcoord_in%sigm_coord(k) > 0.0_ireals) THEN
          zhhl_in(k) = (r_d/g)*LOG( refatm_in%p0sl / vcoord_in%sigm_coord(k)) &
                          * (refatm_in%t0sl - 0.5 * refatm_in%dt0lp *         &
                             LOG(refatm_in%p0sl/vcoord_in%sigm_coord(k)) )
        ELSE
          zhhl_in(k) = undef
        END IF
      ENDDO

      DO k = 1, ke_in
        dak_in(k) =  ak_in(k+1) - ak_in(k)
        dbk_in(k) =  bk_in(k+1) - bk_in(k)
      ENDDO

      ! PE 0 writes the parameters to OUTPUT
      IF (my_cart_id == 0) THEN
        WRITE (noutput,'(A)') '     '
        WRITE (noutput,'(A)') '     Vertical coordinate parameters of input:'
        WRITE (noutput,'(A)') '     '
        WRITE (noutput,'(A)') '     k        ak_in(k) [Pa]   bk_in(k)    vcoord_in(k) [Pa]  zhhl_in(k) [m]'
        DO k = 1,ke_in+1
          WRITE (noutput,'(I6,4F16.4)') k, ak_in(k), bk_in(k), vcoord_in%sigm_coord(k), zhhl_in(k)
        ENDDO
        WRITE (noutput,'(A)') '     '
      ENDIF

    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE get_vert_coord

!==============================================================================
!==============================================================================
!+ Extracts the vertical coordinate parameters for general vertical coordinate
!------------------------------------------------------------------------------

SUBROUTINE get_generalvertical (yshname, ytyofle, igribid, lrequired, idebug)

!------------------------------------------------------------------------------
!
! Description:
!   One PE extracts the vertical coordinate parameters for input model from the
!   meta data of a grib-record and distributes them to all PEs.
!
! Method:
!
!------------------------------------------------------------------------------
!
! Parameterlist
INTEGER(KIND=intgribf), INTENT(IN)    :: igribid
INTEGER(KIND=iintegers),INTENT(IN)    :: idebug
LOGICAL,                INTENT(IN)    :: lrequired

CHARACTER (LEN=*),      INTENT(IN)    :: yshname, ytyofle

! Local Variables
INTEGER  (KIND=iintegers)  ::  &
  izerror, i, izexch, intbuf(3)

CHARACTER (LEN=100)        ::  &
  charbuf

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine get_vert_coord
!------------------------------------------------------------------------------

#ifdef GRIBAPI
  yzroutine = 'get_generalvertical'
  izerror   = 0

  ! Check whether any PE has got vertical coordinates (igds_in(2) /= 0)
  IF (lrequired .AND.                                                           &
     (ytyofle == 'generalVertical') .OR. (ytyofle == 'generalVerticalLayer')) THEN
    izexch = my_cart_id
  ELSE
    izexch = nproc
  ENDIF

  ! Get the minimum of all values in izexch
  IF (num_compute > 1) THEN
    CALL global_values (izexch, 1, 'MIN', imp_integers, icomm_cart,  &
                        -1, yzerrmsg, izerror)
  ENDIF

  IF (idebug > 20) THEN
    PRINT *, ' Processor ', my_cart_id, ' has got field ',           &
       yshname, ytyofle, ' with exchg status: ', izexch
  ENDIF

  ! If this is nproc, no processor has got a generalVertical field
  IF (izexch < nproc) THEN
    ! The processor with lowest id, that has got a generalVertical field
    ! distributes the meta data to all others

    IF (izexch == my_cart_id) THEN
      CALL grib_get (igribid, 'nlev',                vcoord_in%nlevels)
      CALL grib_get (igribid, 'numberOfVGridUsed',   vcoord_in%ivctype)
      CALL grib_get (igribid, 'uuidOfVGrid',         vcoord_in%vc_uuid)
      intbuf(1) = vcoord_in%nlevels
      intbuf(2) = vcoord_in%ivctype
      DO i = 1, 16
        charbuf(i:i) = vcoord_in%vc_uuid(i)
      ENDDO
    ENDIF

    ! First distribute number of vertical coordinate parameters
    IF (num_compute > 1) THEN
      CALL distribute_values (intbuf,   2, izexch, imp_integers,  icomm_cart, izerror)
      CALL distribute_values (charbuf,  1, izexch, imp_character, icomm_cart, izerror)
    ENDIF

    IF (izexch /= my_cart_id) THEN
      vcoord_in%nlevels = intbuf(1)
      vcoord_in%ivctype = intbuf(2)
      DO i = 1, 16
        vcoord_in%vc_uuid(i) = charbuf(i:i)
      ENDDO
    ENDIF
    CALL uuid_2char (vcoord_in%vc_uuid, uuid_in_string)

    ! Set lhhl_in_read to .TRUE., to indicate that a HHL-file must be read
    ! from that HHL-file we determine the necessary vcoord values later on
    lhhl_in_read = .TRUE.
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
#endif

END SUBROUTINE get_generalvertical

!==============================================================================
!==============================================================================
!+ Extracts the vertical coordinate parameters for Unified Model as input model
!==============================================================================

SUBROUTINE get_vert_coord_um (lgetmass, lgetwind)

!==============================================================================
!
! Description:
!
! Method:
!
!==============================================================================
!
! Parameterlist
LOGICAL, INTENT(INOUT)     ::  lgetmass, lgetwind

! Local Variables
INTEGER  (KIND=iintegers)  ::  &
  izerror, n, k, k1, k2,       &
  izvert, izsender_mass, izsender_wind, izrecv (0:num_compute-1), intbuf(ngds)

REAL (KIND=irealgrib)      :: refstf

! Compute zhhl_in for lec2lm
REAL (KIND=ireals)         :: zhhl_in (ke_in + 1)

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine get_vert_coord_um
!------------------------------------------------------------------------------

  yzroutine = 'get_vert_coord_um'
  izerror   = 0

  ! Check whether any PE has got vertical coordinates for mass or wind
  !  - for mass variables   (igds_in(2) /= 0, NOT (ipds(2) == 2 and ipds(7) == 33,34))
  !  - for wind variables   (igds_in(2) /= 0,      ipds(2) == 2 and ipds(7) == 33,34 )
  IF (igds_in(2) /= 0) THEN
    IF ( (ipds(2) == 2) .AND. (ipds(7) == 33 .OR. ipds(7) == 34) ) THEN
      ! this is U or V
      IF (idbg_level > 10) THEN
        PRINT *, '   Processor ', my_cart_id, ' has vertical wind coordinate parameters for UM'
      ENDIF
      izvert = 1000 + INT (igds_in(2), iintegers)
    ELSE
      ! this is another atmospheric variable
      IF (idbg_level > 10) THEN
        PRINT *, '   Processor ', my_cart_id, ' has vertical mass coordinate parameters for UM'
      ENDIF
      izvert = INT (igds_in(2), iintegers)
    ENDIF
  ELSE
    izvert = 0_iintegers
  ENDIF

  IF (num_compute > 1) THEN
    CALL gather_values (izvert, izrecv, 1, num_compute, imp_integers, -1,    &
                        icomm_cart, yzerrmsg, izerror)
  ELSE
    izrecv(0) = izvert
  ENDIF

  izsender_mass = -1
  izsender_wind = -1
  search: DO n = 0, num_compute-1
    IF     ( (izrecv(n) > 0) .AND. (izrecv(n) < 1000) ) THEN
      ! processor n has mass coordinates
      izsender_mass = n
    ELSEIF ( izrecv(n) > 1000) THEN
      ! processor n has wind coordinates
      izsender_wind = n
    ENDIF
  ENDDO search

  IF (idbg_level > 10) THEN
    PRINT *, '   Processors ', izsender_mass, izsender_wind, ' distribute vcoord for mass, wind', lgetmass, lgetwind
  ENDIF
  IF ( (izsender_mass > -1) .AND. (lgetmass) ) THEN
    ! PE with rank izsender_mass got vertical coordinates for mass variables
    ! and sends them to the others
    lgetmass = .FALSE.
    IF (num_compute > 1) THEN
      n = INT (ngds, iintegers)
      IF (my_cart_id == izsender_mass) THEN
        intbuf(1:n) = INT (igds_in(1:n), iintegers)
      ENDIF
      CALL distribute_values(intbuf, n, izsender_mass, imp_integers, &
                             icomm_cart, izerror)
      IF (my_cart_id /= izsender_mass) THEN
        igds_in(1:n) = INT (intbuf(1:n), intgribf)
      ENDIF
    ENDIF
#ifdef GRIBDWD
    ! Every PE now computes the vertical coordinate parameters for mass variables
    k1 = 25
    k2 = 25 + ke_in
    DO k = 1, ke_in
      akh_in(k) = REAL (REFSTF(igds_in(k1 + k)), ireals)
      bkh_in(k) = REAL (REFSTF(igds_in(k2 + k)), ireals)
    ENDDO
#endif
  ENDIF

  IF ( (izsender_wind > -1) .AND. (lgetwind) ) THEN
    ! PE with rank izsender_wind got vertical coordinates for wind variables
    ! and sends them to the others
    lgetwind = .FALSE.
    IF (num_compute > 1) THEN
      n = INT (ngds, iintegers)
      IF (my_cart_id == izsender_wind) THEN
        intbuf(1:n) = INT (igds_in(1:n), iintegers)
      ENDIF
      CALL distribute_values(intbuf, n, izsender_wind, imp_integers, &
                             icomm_cart, izerror)
      IF (my_cart_id /= izsender_wind) THEN
        igds_in(1:n) = INT (intbuf(1:n), intgribf)
      ENDIF
    ENDIF

#ifdef GRIBDWD
    ! Every PE now computes the vertical coordinate parameters for wind variables
    k1 = 25
    k2 = 25 + ke_in
    DO k = 1, ke_in
      akh_in_rho(k) = REAL (REFSTF(igds_in(k1 + k)), ireals)
      bkh_in_rho(k) = REAL (REFSTF(igds_in(k2 + k)), ireals)
    ENDDO
#endif
  ENDIF

  ! Set the other variables concerned with vertical coordinates to undefined
  DO k = 1, ke_in
    ak_in(k)  = HUGE(1.0_ireals)   ! (ak_in(k) + ak_in(k+1)) * 0.5_ireals
    bk_in(k)  = HUGE(1.0_ireals)   ! (bk_in(k) + bk_in(k+1)) * 0.5_ireals
    dak_in(k) = HUGE(1.0_ireals)   !  ak_in(k+1) - ak_in(k)
    dbk_in(k) = HUGE(1.0_ireals)   !  bk_in(k+1) - bk_in(k)
  ENDDO

  ! PE 0 writes the parameters to OUTPUT
  IF ( (.NOT. lgetmass) .AND. (.NOT. lgetwind) ) THEN
    IF (my_cart_id == 0) THEN
      WRITE (noutput,'(A)') '     '
      WRITE (noutput,'(A)') '     Vertical coordinate parameters for main levels of UM:'
      WRITE (noutput,'(A)') '     '
      WRITE (noutput,'(A)') '                           mass                                       wind'
      WRITE (noutput,'(A)') '     k        akh_in(k) [m]    bkh_in(k) [0..1]      akh_in_rho(k) [m]   bkh_in_rho(k) [0..1] '
      DO k = 1,ke_in
        WRITE (noutput,'(I6,4F16.4)') k, akh_in(k), bkh_in(k), akh_in_rho(k), bkh_in_rho(k)
      ENDDO
      WRITE (noutput,'(A)') '     '
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE get_vert_coord_um

!==============================================================================
!+ Subroutine for distributing the initial- and the boundary data
!------------------------------------------------------------------------------

SUBROUTINE scatter_data (procarrays, index, iloc, irank, ilev, leof,        &
                         lrequired, istat)

!------------------------------------------------------------------------------
!
! Description:
!  This subroutine is called within the read-loop over all records (with the
!  loop index "index"). In this loop up to num_compute (number of compute PE)
!  records are read and distributed to the processors. Each processor gets a
!  total record for decoding. After the decoding distribute_subarrays does the
!  distribution (with scatter_values). It works for all numbers of processors.
!  The distributed subarrays are then copied to the corresponding variables.
!
! Method:
!  Distribute the appropriate part of array procarrays to the corresponding
!  processors (in parallel mode) or just copy it into subarray (in sequential
!  mode).
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
INTEGER (KIND=iintegers), INTENT(IN)    ::  &
  iloc, irank, ilev,                 & ! characteristics of the variable
  index                                ! index of the distribution loop

LOGICAL,                  INTENT(IN)    ::  &
  leof, lrequired                      ! end of file and required field

REAL    (KIND=ireals),    INTENT(IN)    ::  &
  procarrays (ie_in_max*je_in_max, 0:num_compute-1) ! decomposed total field

INTEGER (KIND=iintegers), INTENT(OUT)   ::  &
  istat                                ! go on, cycle or exit the loop

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)   :: implcode
INTEGER (KIND=iintegers)   :: iz_info(4), i,j, izstorelev
CHARACTER (LEN=25)         :: yzroutine
CHARACTER (LEN=80)         :: yzerrmsg

! Local arrays:
REAL    (KIND=ireals)      :: zsubarray_1d(ie_in_max*je_in_max),    &
                              zsubarray_2d(ie_in_max,je_in_max)
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine scatter_data
!------------------------------------------------------------------------------

  istat    = 0
  yzroutine = 'scatter_data'

  ! check and distribute the next action
  IF (index == my_cart_id) THEN
    IF ( (.NOT.leof) .AND. (lrequired) ) THEN
      iz_info(1) = 0 ! Processor has data
      iz_info(2) = iloc
      iz_info(3) = irank
      iz_info(4) = ilev
    ELSE
      IF (leof) THEN
        iz_info(1) = -1 ! No data because EOF reached
      ELSE
        iz_info(1) = -2 ! No data because it is not required
      ENDIF
    ENDIF
  ENDIF

  IF (num_compute > 1) THEN
    CALL distribute_values (iz_info, 4, index, imp_integers, icomm_cart,   &
                            implcode)
    IF (implcode /= 0) THEN
      istat   = 1
      RETURN
    ENDIF
  ENDIF

  istat = iz_info(1)
  IF (lgsm2lm) THEN
    izstorelev = ke_in + 1 - iz_info(4)
  ELSE
    izstorelev = iz_info(4)
  ENDIF

  ! Distribute the records
  IF (istat == 0) THEN
    ! Update the values in the record var_in (ar_des_input)
    var_in(iz_info(2))%lreadin = .TRUE.
    var_in(iz_info(2))%nlevels_read = var_in(iz_info(2))%nlevels_read + 1

    IF (num_compute > 1) THEN
      CALL scatter_values (procarrays, zsubarray_1d, ie_in_max*je_in_max, &
                           num_compute, imp_reals, index, icomm_cart,     &
                           yzerrmsg, implcode)
      IF (implcode /= 0) THEN
        istat = 2
        RETURN
      ENDIF
    ELSE
      zsubarray_1d(:) = procarrays(:,0)
    ENDIF

    zsubarray_2d = RESHAPE (zsubarray_1d, (/ie_in_max, je_in_max/))

    IF ((idbg_level >= 20) .AND. (my_cart_id == 0) ) THEN
      PRINT *, '   Storing variable (name, loc, rank, lev): ',        &
                var_in(iz_info(2))%name,  iz_info(2), iz_info(3), izstorelev
    ENDIF

    ! put data into the variables
    IF(iz_info(3) == 3 ) THEN
      var_in(iz_info(2))%p3(1:ie_in,1:je_in,izstorelev)                    &
                                          = zsubarray_2d(1:ie_in,1:je_in)
    ELSE IF(iz_info(3) == 2 ) THEN
      var_in(iz_info(2))%p2(1:ie_in,1:je_in) = zsubarray_2d(1:ie_in,1:je_in)
    ENDIF

    ! this was a change by Christian Koziar
    ! lowest level of QV is also stored in QV_S
!US But there is the possibility to get QV_S from ECWMF?????
!US    IF ( (var_in(iz_info(2))%name == 'QV      ') .AND. (iz_info(4) == ke_in) ) THEN
!US      var_in(27)%p2(1:ie_in,1:je_in) = zsubarray_2d(1:ie_in,1:je_in)
!US      var_in(27)%lreadin = .TRUE.
!US      var_in(27)%nlevels_read = var_in(iz_info(2))%nlevels_read + 1
!US    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE scatter_data

!==============================================================================
!==============================================================================

FUNCTION get_gfs_level (ilevel_in_hpa)

!------------------------------------------------------------------------------
!
! Description:
!  GFS data is given on pressure levels: isobaricInhPa
!  and the level is given in hPa.
!  A priori we cannot identify, how many and which levels are given, so we
!  assume a 26-level version with the levels:
!     10,  20,  30,  50,  70, 100, 150, 200, 250, 300,
!    350, 400, 450, 500, 550, 600, 650, 700, 750, 800,
!    850, 900, 950, 975, 1000
!  and these levels are put into the variables: uppermost level = 1, etc.
!
!  Note: the humidity variables QV, QC are only given from 100 hPa on
!------------------------------------------------------------------------------

! Scalar arguments with intent(in):
INTEGER (KIND=iintegers),   INTENT(IN) :: &
  ilevel_in_hpa       ! level in hPa

! Type of function:
INTEGER (KIND=iintegers)               :: &
  get_gfs_level       ! value in the range 1 ... 26

! Local variable
INTEGER (KIND=iintegers)               :: ilev, klev, k

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

klev = -1   ! is set in the loop below
ilev = ilevel_in_hpa*100_iintegers

kloop:  DO k = 1, ke_in
  IF (NINT(press_level(k),iintegers) == ilev) THEN
    klev = k
    EXIT kloop
  ENDIF
ENDDO kloop

get_gfs_level = klev

END FUNCTION get_gfs_level

!==============================================================================
!+ Module procedure in "src_read_coarse_grid" for reading a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE read_nc_gdefs (ncid, iedim, jedim, kedim, kesoildim, ydate,     &
               startlon, startlat, endlon, endlat, pollon, pollat, polgam, &
               east_add_in, west_add_in, south_add_in, north_add_in,       &
               icomm, myid, npes, lget, lget_time, lget_soil, yerrmsg, istatus)

!------------------------------------------------------------------------------
!
! Description:
!   This routine reads global attributes from NetCDF formatted input file.
!
!------------------------------------------------------------------------------
!
!
! Subroutine / Function arguments

! Scalar arguments with intent(in):
INTEGER (kind=iintegers),   INTENT(IN) :: &
  ncid,                & ! NetCDF file ID
  icomm,               & ! MPI communicator
  myid,                & ! ID of this PE in communicator icomm
  npes,                & ! number of PEs
  iedim, jedim,        & ! dimensions of the input fields
  kedim, kesoildim,    & ! dimensions of the input fields
  east_add_in,         & ! add an extra column to the East
  west_add_in,         & ! add an extra column to the West
  south_add_in,        & ! add an extra column to the South
  north_add_in           ! add an extra column to the North

CHARACTER (LEN=*),        INTENT(IN)  ::  &
  ydate                  ! actual date from Namelist parameter

REAL    (KIND=ireals),    INTENT(IN)  ::  &
  startlat, startlon, &  ! coordinates of the lower left grid point
  endlat, endlon,     &  ! coordinates of the upper right grid point
  pollon, pollat, polgam ! coordinates of the geographical North Pole

LOGICAL, INTENT(IN)     ::  &
  lget,                     &
  lget_time,                &
  lget_soil


! Scalar arguments with intent(out):
INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
  istatus          ! error index

CHARACTER (LEN= *),        INTENT(OUT)   ::  &
  yerrmsg          ! string for error messages

!-----------------------------------------------------------------------
!
! Local scalars:

INTEGER (KIND=iintegers), PARAMETER  :: ntbds=2, ntime=1

INTEGER (KIND=iintegers) ::  &
  i, j, k,                                     & ! loop variable
  izerror, izmplcode,                          & ! error status variable
  izedim_in, jzedim_in, kzedim_in, kze1dim_in, & ! input dimensions
  kesoildim_in, ke1soildim_in,                 & !
  nztime_in, nztbds_in,                        & !
  izdb                                           ! index for finding a string

  ! local variables for replacing reference atmosphere and vertical grid
  ! parameters for putting them into new structures
INTEGER (KIND=iintegers) ::  &
  izvctype_in, izrefatm_in

REAL (KIND=ireals) :: &
  zp0sl_in, zt0sl_in, zdt0lp_in, zvcflat_in, zdelta_t_in, zh_scal_in,    &
  zvcoord_in(kedim+1)

INTEGER (kind=iintegers) ::   &
  jzgridVarID,     & ! NetCDF ID for grid mapping
  jzlonVarID,      & ! NetCDF ID for longitude
  jzlatVarID,      & ! NetCDF ID for latitude
  jzsoilVarID,     & ! NetCDF ID for soil layers
  jzvcVarID ,      & ! NetCDF ID for the vertical component
  jztimeID           ! NetCDF ID for the time

INTEGER (KIND=iintegers) :: &
  izyear_ref, izmonth_ref, izday_ref, izhour_ref, izmin_ref, izsec_ref,     &
  izyear, izmonth, izday, izhour, izmin, izsec, izref_act, izerrf, iztimepassed

REAL (KIND=ireals) :: &
  zpollat,       & ! latitude of North Pole in input data
  zpollon,       & ! longitude of North Pole in input data
  zpolgam          ! angle between the north poles of the systems

CHARACTER (LEN=80) :: &
  yzgrid_name        ! name of grid mapping

CHARACTER (LEN=40) :: &
  yzdate1           ! time unit and reference date from input

CHARACTER (LEN=14) ::  &    !_br 06.09.12 adaptation to new file name convention
  yzdate_ref        ! reference date from input

CHARACTER (LEN=20) ::  &
  yzcalendar        ! type of calendar  !_br 10.09.12

CHARACTER (LEN=4) :: &
  yzlon, yzlat       ! names for variables langitude and latitude

!  Local arrays:

REAL (kind=ireals) :: &
  ztime(ntime),   &    ! forecast time  !_br 10.09.12
  ztimemul             ! conversion factor for time in seconds  !_br 10.09.12

REAL (kind=ireals), ALLOCATABLE  :: &
  zczml_soil_in(:)

REAL (kind=ireals), ALLOCATABLE   :: &
  zlongitudes_in(:),  & ! rotated longitudes
  zlatitudes_in(:)      ! rotated latitudes

CHARACTER (LEN=3) :: &
  yzatt_value       ! attribute value

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

#ifdef NETCDF
!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

 ! initialize error status variable
  yerrmsg       = '   '
  izerror        = 0
  istatus       = 0

! Processor 0 does the job
  IF (myid == 0) THEN
    ! Get the dimension ID's and the length of the dimensions
    istatus = nf90_inq_varid(ncid, 'rotated_pole', jzgridVarID)
    IF (istatus == NF90_NOERR) THEN
      yzlon = 'rlon'
      yzlat = 'rlat'
    ELSE
      yzlon = 'lon'
      yzlat = 'lat'
    ENDIF

    istatus = nf90_inq_dimid (ncid, TRIM(yzlon), idims_id(1))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    istatus = nf90_inquire_dimension (ncid, idims_id(1), len=izedim_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF ((iedim - east_add_in - west_add_in) /= izedim_in) izerror = 1

    istatus = nf90_inq_dimid (ncid, TRIM(yzlat), idims_id(2))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    istatus = nf90_inquire_dimension (ncid, idims_id(2), len=jzedim_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF ((jedim - south_add_in - north_add_in)  /= jzedim_in) izerror = 2

    istatus = nf90_inq_dimid (ncid, 'level', idims_id(3))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    istatus = nf90_inquire_dimension (ncid, idims_id(3), len=kzedim_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF (kedim /= kzedim_in) izerror = 3

    IF ( .NOT. lcm_pres_coor) THEN
      istatus = nf90_inq_dimid (ncid, 'level1', idims_id(4))
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF

      istatus = nf90_inquire_dimension (ncid, idims_id(4), len=kze1dim_in)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      IF ((kedim+1) /= kze1dim_in) izerror = 3
    ENDIF

    ALLOCATE (zlongitudes_in(iedim), zlatitudes_in(jedim), STAT=istatus)
    IF (istatus /= 0) THEN
      yerrmsg = 'Allocation error in read_nc_gdefs'
      RETURN
    ENDIF

    istatus = nf90_inq_dimid (ncid, 'time', idims_id(5))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    istatus = nf90_inquire_dimension (ncid, idims_id(5), len=nztime_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    IF (nztime_in /= 1) izerror = 4

    istatus = nf90_inq_dimid (ncid, 'bnds', idims_id(6))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    istatus = nf90_inquire_dimension (ncid, idims_id(6), len=nztbds_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    IF (nztbds_in /= 2) izerror = 4

    IF ((llm2lm .OR. (itype_w_so_rel > 0)) .AND. lget_soil) THEN

      istatus=nf90_inq_dimid(ncid, 'soil1', idims_id(8))
      IF (istatus /= 0) istatus=nf90_inq_dimid(ncid, 'soil_depth', idims_id(8))

      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF

      istatus = nf90_inquire_dimension (ncid, idims_id(8), len=ke1soildim_in)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF

      IF (kesoildim+1 /= ke1soildim_in) izerror = 5

      ! get the values of the soil layers
      istatus = nf90_inq_varid (ncid, 'soil1', jzsoilVarID)
      IF (istatus /= 0) istatus = nf90_inq_varid (ncid, 'soil_depth', jzsoilVarID)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF

      ALLOCATE (zczml_soil_in(ke1soildim_in), STAT=istatus)
      IF (istatus /= 0) THEN
        yerrmsg = 'Allocation error in read_nc_gdefs'
        RETURN
      ENDIF

      istatus = nf90_get_var (ncid, jzsoilVarID, zczml_soil_in)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF

      istatus = 0
      DO j = 1, ke1soildim_in
        IF (ABS(zczml_soil_in(j) - czmls_in(j)) > 0.001_ireals) istatus = -1
      ENDDO
      IF (istatus /= 0) THEN
        PRINT *, 'soil layer     data file     namelist input'
        DO j = 1, ke1soildim_in
          WRITE(*,'(I8,2F14.3)') j, zczml_soil_in(j), czmls_in(j)
        ENDDO
        istatus = 30
        yerrmsg = 'Error in multi soil layer heights'
        RETURN
      ENDIF

      DEALLOCATE (zczml_soil_in)

    ENDIF ! (llm2lm .OR. (itype_w_so_rel == 1)) .AND. lget_soil


    ! Get longitude and latitude data

    istatus = nf90_inq_varid (ncid, TRIM(yzlon), jzlonVarID)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    istatus = nf90_get_var (ncid, jzlonVarID, zlongitudes_in(1+west_add_in:iedim-east_add_in))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF (west_add_in /= 0) THEN
      zlongitudes_in(1) = zlongitudes_in(iedim - east_add_in) - 360._ireals
    ENDIF
    IF (east_add_in /= 0) THEN
      zlongitudes_in(iedim) = 360._ireals + zlongitudes_in(1 + west_add_in)
    ENDIF

    istatus = nf90_inq_varid (ncid, TRIM(yzlat), jzlatVarID)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    istatus = nf90_get_var (ncid, jzlatVarID, zlatitudes_in(1+south_add_in:jedim-north_add_in))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF (south_add_in /= 0) THEN
      zlatitudes_in(1) = -90._ireals
    ENDIF
    IF (north_add_in /= 0) THEN
      zlatitudes_in(jedim) = 90._ireals
    ENDIF

    !   get grid mapping values
    istatus = nf90_inq_varid(ncid, 'rotated_pole', jzgridVarID)
    IF (istatus == NF90_NOERR) THEN   ! rotated latitude-longitude
      istatus = nf90_get_att(ncid, jzgridVarID, 'grid_mapping_name', yzgrid_name)
      IF (istatus /= NF90_NOERR) THEN
        PRINT *, TRIM(NF90_strerror(istatus))
        yerrmsg = 'Error in read_nc_gdefs / nf90_get_att - attribute "grid_mapping_name"'
        RETURN
      ENDIF
      IF (yzgrid_name(1:26) /= 'rotated_latitude_longitude') THEN
        istatus=-1
        yerrmsg = 'Error in read_nc_gdefs - Invalid value for attribute "grid_mapping_name": '//TRIM(yzgrid_name)
        RETURN
      ENDIF
      istatus = nf90_get_att(ncid, jzgridVarID, 'grid_north_pole_latitude', zpollat)
      IF (istatus /= NF90_NOERR) THEN
        PRINT *, TRIM(NF90_strerror(istatus))
        yerrmsg = 'Error in read_nc_gdefs / nf90_get_att - attribute "grid_north_pole_latitude"'
        RETURN
      ENDIF
      istatus = nf90_get_att(ncid, jzgridVarID, 'grid_north_pole_longitude', zpollon)
      IF (istatus /= NF90_NOERR) THEN
        PRINT *, TRIM(NF90_strerror(istatus))
        yerrmsg = 'Error in read_nc_gdefs / nf90_get_att - attribute "grid_north_pole_longitude"'
        RETURN
      ENDIF
      istatus = nf90_get_att(ncid, jzgridVarID, 'north_pole_grid_longitude', zpolgam)
      IF (istatus == NF90_ENOTATT) THEN
        zpolgam = 0._ireals
      ELSE IF (istatus /= NF90_NOERR) THEN
        PRINT *, TRIM(NF90_strerror(istatus))
        yerrmsg = 'Error in read_nc_gdefs / nf90_get_att - attribute "north_pole_grid_longitude"'
        RETURN
      ENDIF
    ELSE IF (istatus == NF90_enotvar) THEN  ! true latitude-longitude
      zpollon = 180.
      zpolgam =   0.
      zpollat = 90.
    ELSE
      PRINT *, 'Error in read_nc_gdefs / nf90_inq_varid'
      PRINT *, 'Variable "rotated_pole"'
      RETURN
    ENDIF

    ! compare the input values with the Namelist parameters
    IF (ABS(zpollat - pollat) > 1.0E-3_ireals) izerror = 4
    IF (ABS(zpollon - pollon) > 1.0E-3_ireals) izerror = 5
    IF (ABS(zpolgam - polgam) > 1.0E-3_ireals) izerror = 5

    IF (ABS(zlatitudes_in(1+south_add_in)  - startlat) > 1.0E-3_ireals) izerror = 6
    IF (ABS(zlongitudes_in(1+west_add_in)  - startlon) > 1.0E-3_ireals) izerror = 7

    IF (ABS(zlatitudes_in(jedim-north_add_in)  - endlat) > 1.0E-3_ireals) izerror = 6
    IF (ABS(zlongitudes_in(iedim-east_add_in)  - endlon) > 1.0E-3_ireals) izerror = 7

    ! Get the vertical co-ordinate
    IF (lget) THEN

      IF (llm2lm) THEN
        istatus = nf90_inq_varid (ncid, 'vcoord', jzvcVarID)
        IF (istatus /= NF90_NOERR) THEN
          yerrmsg = 'vcoord '// TRIM(NF90_strerror(istatus))
          RETURN
        ENDIF
        istatus = nf90_get_var (ncid, jzvcVarID, zvcoord_in)
        IF (istatus /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(istatus))
          RETURN
        ENDIF
        istatus = nf90_get_att (ncid, jzvcVarID, 'p0sl', zp0sl_in)
        IF (istatus /= NF90_NOERR) THEN
          yerrmsg = 'p0sl '//TRIM(NF90_strerror(istatus))
          RETURN
        ENDIF
        istatus = nf90_get_att (ncid, jzvcVarID, 't0sl', zt0sl_in)
        IF (istatus /= NF90_NOERR) THEN
          yerrmsg = 't0sl '//TRIM(NF90_strerror(istatus))
          RETURN
        ENDIF
        istatus = nf90_get_att (ncid, jzvcVarID, 'dt0lp', zdt0lp_in)
        IF (istatus /= NF90_NOERR) THEN
          yerrmsg = 'dt0lp '//TRIM(NF90_strerror(istatus))
          RETURN
        ENDIF
        istatus = nf90_get_att (ncid, jzvcVarID, 'vcflat', zvcflat_in)
        IF (istatus /= NF90_NOERR) THEN
          yerrmsg = 'vcflat '//TRIM(NF90_strerror(istatus))
          RETURN
        ENDIF
        istatus = nf90_get_att (ncid, jzvcVarID, 'irefatm', izrefatm_in)
        IF (istatus == NF90_NOERR) THEN
          istatus = nf90_get_att (ncid, jzvcVarID, 'ivctype', izvctype_in)
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = 'ivctype_in '//TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF
          IF (izrefatm_in == 2) THEN ! reference atmosphere type 2
            istatus = nf90_get_att (ncid, jzvcVarID, 'delta_t', zdelta_t_in)
            IF (istatus /= NF90_NOERR) THEN
              yerrmsg = 'delta_t_in '//TRIM(NF90_strerror(istatus))
              RETURN
            ENDIF
            istatus = nf90_get_att (ncid, jzvcVarID, 'h_scal', zh_scal_in)
            IF (istatus /= NF90_NOERR) THEN
              yerrmsg = 'h_scal_in '//TRIM(NF90_strerror(istatus))
              RETURN
            ENDIF
          ENDIF
          IF (izvctype_in == 3 .OR. izvctype_in == 4) THEN ! SLEVE coordinates
            istatus = nf90_get_att (ncid, jzvcVarID, 'svc1', svc1_in)
            IF (istatus /= NF90_NOERR) THEN
              yerrmsg = 'svc1_in '//TRIM(NF90_strerror(istatus))
              RETURN
            ENDIF
            istatus = nf90_get_att (ncid, jzvcVarID, 'svc2', svc2_in)
            IF (istatus /= NF90_NOERR) THEN
              yerrmsg = 'svc2_in '//TRIM(NF90_strerror(istatus))
              RETURN
            ENDIF
            istatus = nf90_get_att (ncid, jzvcVarID, 'nfltvc', nfltvc_in)
            IF (istatus /= NF90_NOERR) THEN
              yerrmsg = 'nfltvc_in '//TRIM(NF90_strerror(istatus))
              RETURN
            ENDIF
          ENDIF
        ELSE ! no attribute 'irefatm' and NF90_NOERR is false: old version
          izrefatm_in = 1
          izvctype_in = 1
        ENDIF

        ! Now all vertical coordinate parameters and ivctype_in are set
        ! go on with some checks

        IF     (izvctype_in == 1) THEN
          ! For this type the vertical coordinates should be ascending
          IF ( zvcoord_in(2) < zvcoord_in(1) ) THEN
            PRINT *,                                                         &
              ' ERROR *** Vertical coordinates not ascending for type *** ', &
               izvctype_in, zvcoord_in(1), zvcoord_in(2)
            istatus = 5
            RETURN
          ENDIF
        ELSEIF (izvctype_in == 2) THEN
          ! For this type the vertical coordinates should be descending
          IF ( zvcoord_in(2) > zvcoord_in(1) ) THEN
            PRINT *,                                                         &
               ' ERROR *** Vertical coordinates not descending for type ***',&
                izvctype_in, zvcoord_in(1), zvcoord_in(2)
            istatus = 6
            RETURN
          ENDIF
        ELSEIF (izvctype_in == 3 .OR. izvctype_in == 4) THEN
          ! Check for meaningful values of svc1, svc2 and nfltvc
          IF ((svc1_in > zvcoord_in(1)) .OR. (svc1_in < 0.0)) THEN
            PRINT *, ' ERROR *** svc1_in not in allowed '//                  &
                            'range for ivctype_in = 3 or 4 ***'
            istatus = 7
            RETURN
          ENDIF
          IF ((svc2_in > zvcoord_in(1)) .OR. (svc2_in < 0.0)) THEN
            PRINT *,  ' ERROR *** svc2_in not in allowed '//                 &
                             'range for ivctype_in = 3 or 4 ***'
            istatus = 8
            RETURN
          ENDIF
          IF (nfltvc_in <= 0) THEN
            PRINT *,  ' ERROR *** nfltvc_in must be greater than '//         &
                                            'or equal to zero ***'
            istatus = 9
            RETURN
          ENDIF
        ELSE
          PRINT *,  ' ERROR *** Type ivctype_in of vertical coordinate '//   &
                                                    'not available***'
          istatus = 11
          RETURN
        ENDIF

      ELSE ! lcm2lm

        IF (lcm_pres_coor) THEN

          ! ak, bk already defined in org_vert_interpol_p2h
          ! have to read pressure levels

          istatus = nf90_inq_varid (ncid, 'level', jzvcVarID)
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = 'level '// TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_get_var (ncid, jzvcVarID, press_level(1:kedim))
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_get_att (ncid, jzvcVarID, 'units', yzatt_value)
          IF (TRIM(yzatt_value) == "hPa") press_level(:) = press_level(:) * 100._ireals
          ! include printout of Vertical coordinate parameters
          WRITE (noutput,'(A)') '     '
          WRITE (noutput,'(A)') '     Vertical coordinate parameters of input:'
          WRITE (noutput,'(A)') '     '
          WRITE (noutput,'(A)') '             k     (Pa)'

          DO k = 1,kzedim_in
            WRITE (noutput,'(I14,2F17.4)') k, press_level(k)
          ENDDO
          WRITE (noutput,'(A)') '     '

        ELSEIF (lcm_hgt_coor) THEN

          istatus = nf90_inq_varid (ncid, 'ak', jzvcVarID)
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = 'ak '// TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_get_var (ncid, jzvcVarID, ak_in(1:kzedim_in))
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_inq_varid (ncid, 'bk', jzvcVarID)
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = 'bk '// TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_get_var (ncid, jzvcVarID, bk_in(1:kzedim_in))
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_inq_varid (ncid, 'ak_uv', jzvcVarID)
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = 'ak_uv '// TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_get_var (ncid, jzvcVarID, akh_in_rho)
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_inq_varid (ncid, 'bk_uv', jzvcVarID)
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = 'bk_uv '// TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_get_var (ncid, jzvcVarID, bkh_in_rho)
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          DO k = 1, kzedim_in
            akh_in(k) = ak_in(k)
            bkh_in(k) = bk_in(k)
            dak_in(k) = ak_in(k+1) - ak_in(k)
            dbk_in(k) = bk_in(k+1) - bk_in(k)
          ENDDO

          WRITE (noutput,'(A)') '     '
          WRITE (noutput,'(A)') '     Vertical coordinate parameters of input:'
          WRITE (noutput,'(A)') '     '
          WRITE (noutput,'(A)') '             k        a(k) (m)              b(k)'
          DO k = 1,kzedim_in
            WRITE (noutput,'(I14,2F17.4)') k, ak_in(k), bk_in(k)
          ENDDO
          WRITE (noutput,'(A)') '     '

        ELSE   ! .NOT. lcm_pres_coor   and   .NOT. lcm_hgt_coor
          istatus = nf90_inq_varid (ncid, 'ak', jzvcVarID)
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = 'ak '// TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_get_var (ncid, jzvcVarID, ak_in(1:kzedim_in+1))
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_inq_varid (ncid, 'bk', jzvcVarID)
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = 'bk '// TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF

          istatus = nf90_get_var (ncid, jzvcVarID, bk_in(1:kzedim_in+1))
          IF (istatus /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(istatus))
            RETURN
          ENDIF


          DO k = 1, kzedim_in
            akh_in(k) = (ak_in(k) + ak_in(k+1)) * 0.5_ireals
            bkh_in(k) = (bk_in(k) + bk_in(k+1)) * 0.5_ireals
            dak_in(k) =  ak_in(k+1) - ak_in(k)
            dbk_in(k) =  bk_in(k+1) - bk_in(k)
          ENDDO

          ! include printout of Vertical coordinate parameters
          WRITE (noutput,'(A)') '     '
          WRITE (noutput,'(A)') '     Vertical coordinate parameters of input:'
          WRITE (noutput,'(A)') '     '
          WRITE (noutput,'(A)') '             k        a(k) (Pa)             b(k)'
          DO k = 1,kzedim_in+1
            WRITE (noutput,'(I14,2F17.4)') k, ak_in(k), bk_in(k)
          ENDDO
          WRITE (noutput,'(A)') '     '

        ENDIF   ! lcm_pres_coor, lcm_hgt_coor or other

      ENDIF ! llm2lm

    ENDIF  ! lget

!   Get time data
    IF (lget_time) THEN

      istatus = nf90_inq_varid (ncid, 'time', jztimeID)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      istatus = nf90_get_var (ncid, jztimeID, ztime)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      istatus = nf90_get_att (ncid, jztimeID, 'calendar', yzcalendar) !_br 10.09.12
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      IF (((yzcalendar(1:8) == 'standard' .OR. yzcalendar(1:19) == 'proleptic_gregorian') &
                         .AND. (itype_calendar == 0)) .OR. &
           (yzcalendar(1:7) == '360_day' .AND. (itype_calendar == 1)) .OR. &
           (yzcalendar(1:7) == '365_day' .AND. (itype_calendar == 2))) THEN
         ! callendar attribute is valid
      ELSE
        IF (myid == 0) THEN
          PRINT *, 'calendar attribute = >',yzcalendar,'<'
          PRINT *, ' but wrong itype_calendar = ', itype_calendar
        ENDIF
        yerrmsg = ' ERROR *** Wrong calendar attribute ***'
        istatus = -1
        RETURN
      ENDIF
      istatus = nf90_get_att (ncid, jztimeID, 'units', yzdate1)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF

      ! check the date for which the product is valid
!_br 06.09.12
      !   determine time units
      izdb = index(yzdate1,'seconds')
      IF (izdb /= 0) THEN
        ztimemul = 1._ireals
      ELSE
        izdb = index(yzdate1,'minutes')
        IF (izdb /= 0) THEN
          ztimemul = 60._ireals
        ELSE
          izdb = index(yzdate1,'hours')
          IF (izdb /= 0) THEN
            ztimemul = 3600._ireals
          ELSE
            izdb = index(yzdate1,'days')
            IF (izdb /= 0) THEN
              ztimemul = 86400._ireals
            ELSE
              yerrmsg = 'No valid time units found. Must be seconds, minutes, hours or days'
              istatus = 11
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      izdb = index(yzdate1,'since') + 5
      DO i = 1,LEN_TRIM(yzdate1)
        IF(ICHAR(yzdate1(i:i)) < 48 .OR. ICHAR(yzdate1(i:i)) > 57) yzdate1(i:i)=' '
      ENDDO
        READ(yzdate1(izdb:),*,IOSTAT=istatus) izyear_ref,izmonth_ref,izday_ref,izhour_ref,izmin_ref,izsec_ref
        IF (istatus /= 0) THEN
          yerrmsg = "ERROR in units attribute of time in the netCDF input"
          RETURN
        ENDIF
      IF (lmmss_ini) THEN
        WRITE(yzdate_ref,'(I4.4,5(I2.2))') izyear_ref,izmonth_ref,izday_ref,izhour_ref,izmin_ref,izsec_ref
      ELSE
        WRITE(yzdate_ref,'(I4.4,5(I2.2))') izyear_ref,izmonth_ref,izday_ref,izhour_ref
        izmin_ref = 0
        izsec_ref = 0
      ENDIF

      ! Determine the corresponding integer values for the actual date "ydate"
      ! no minutes are given here
      IF (lmmss_ini) THEN
        READ(ydate    ,'(I4,5I2)') izyear, izmonth, izday, izhour, izmin, izsec
      ELSE
        READ(ydate    ,'(I4,3I2)') izyear, izmonth, izday, izhour
        izmin = 0
        izsec = 0
      ENDIF

      ! Determine the difference between reference date and actual date
      ! in seconds
      CALL diff_minutes(izyear_ref, izmonth_ref, izday_ref, izhour_ref, izmin_ref, &
                        izyear,     izmonth,     izday,     izhour,     izmin,     &
                        itype_calendar, izref_act,  izerrf)

      ! add difference in seconds
      izref_act = izref_act * 60._ireals + izsec - izsec_ref

      ! Now determine the time passed (since the reference date) in seconds
      iztimepassed = ztime(1) * ztimemul

      ! The timepassed and the difference of the dates must now be equal:
      IF (izref_act /= iztimepassed) THEN
        izerror = 10
      ENDIF

      IF (izerror  == 10) THEN
        WRITE(*,'(A,I4.4,5(A,I2.2))') ' The actual date           ', &
             izyear,'-',izmonth,'-',izday,' ',izhour,':',izmin,':',izsec
        WRITE(*,'(A,I4.4,5(A,I2.2))') ' and the reference date    ', &
             izyear_ref,'-',izmonth_ref,'-',izday_ref,' ',izhour_ref,':',izmin_ref,':',izsec_ref
        PRINT *, 'do not match!!'
        PRINT *, 'time difference is ', izref_act,' s, but should be ',ztime(1), ' s'
        istatus = izerror
        PRINT *, 'istatus = ',istatus
        yerrmsg = ' The actual date and the reference date do not match !!'
        RETURN
      ENDIF
!_br 06.09.12 end

    ENDIF ! lget_time

    istatus = izerror

!------------------------------------------------------------------------------
! Section 4: If necessary, print errors
!------------------------------------------------------------------------------

    IF (istatus /= 0 .and. istatus < 10) THEN
          PRINT *, '         data file          namelist input'

          PRINT *, 'ie_tot                ',izedim_in     ,'       ',iedim
          PRINT *, 'je_tot                ',jzedim_in     ,'       ',jedim
          PRINT *, 'ke_tot                ',kzedim_in     ,'       ',kedim
          PRINT *, 'ke_tot + 1            ',kze1dim_in    ,'       ',kedim+1

!     activate this when multilayer soil is valid for input
          IF (llm2lm) THEN
            PRINT *, 'ke_soil_in             ',ke1soildim_in-1 ,'       ',kesoildim
            PRINT *, 'ke_soil_in + 1         ',ke1soildim_in,   '       ',kesoildim+1
          ENDIF

          PRINT *, 'startlat_in           ',zlatitudes_in(1+south_add_in)  ,'       ',startlat
          PRINT *, 'startlon_in           ',zlongitudes_in(1+west_add_in) ,'       ',startlon

          PRINT *, 'endlat_in             ',zlatitudes_in(jedim-north_add_in)   ,'       ',endlat
          PRINT *, 'endlon_in             ',zlongitudes_in(iedim-east_add_in)  ,'       ',endlon

          PRINT *, 'pollat_in             ',zpollat,'       ', pollat
          PRINT *, 'pollon_in             ',zpollon,'       ', pollon
          PRINT *, 'polgam_in             ',zpolgam,'       ', polgam
!         PIK U. Boehm - 21.11.06
          yerrmsg = 'grid structure in input values and in Namelist parameters differs'
!         PIK U. Boehm - End
    ENDIF

    DEALLOCATE (zlongitudes_in, zlatitudes_in, STAT=istatus)
    IF (istatus /= 0) THEN
      yerrmsg = 'Allocation error in read_nc_gdefs'
      RETURN
    ENDIF

  ENDIF ! myid=0

! send information to the other processors
  IF (npes > 1) THEN
    CALL MPI_BCAST (izvctype_in,    1, imp_integers, 0, icomm, izmplcode)
    CALL MPI_BCAST (izrefatm_in,    1, imp_integers, 0, icomm, izmplcode) !_br 31.03.09
    CALL MPI_BCAST (zp0sl_in,       1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (zt0sl_in,       1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (zdt0lp_in,      1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (zvcflat_in,     1, imp_reals,    0, icomm, izmplcode)

    IF (izrefatm_in == 2) THEN
      CALL MPI_BCAST (zdelta_t_in,    1, imp_reals,    0, icomm, izmplcode)
      CALL MPI_BCAST (zh_scal_in,     1, imp_reals,    0, icomm, izmplcode)
    ENDIF

    CALL MPI_BCAST (zvcoord_in, ke1in, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (ak_in,     ke1in, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (bk_in,     ke1in, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (dak_in,    ke_in, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (dbk_in,    ke_in, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (akh_in,    ke_in, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (bkh_in,    ke_in, imp_reals,    0, icomm, izmplcode)

    IF (izvctype_in == 3 .OR. izvctype_in == 4) THEN  ! SLEVE coordinate
      CALL MPI_BCAST (svc1_in,     1, imp_reals,    0, icomm, izmplcode)
      CALL MPI_BCAST (svc2_in,     1, imp_reals,    0, icomm, izmplcode)
      CALL MPI_BCAST (nfltvc_in,   1, imp_integers, 0, icomm, izmplcode)
    ENDIF

    IF (lcm_hgt_coor) THEN
      CALL MPI_BCAST (akh_in_rho,  ke_in, imp_reals,    0, icomm, izmplcode)
      CALL MPI_BCAST (bkh_in_rho,  ke_in, imp_reals,    0, icomm, izmplcode)
    ENDIF

    IF (lcm_pres_coor) THEN
      CALL MPI_BCAST (press_level,  ke_in, imp_reals,    0, icomm, izmplcode)
    ENDIF

  ENDIF

  IF (lget) THEN
    ! transfer parameters for reference atmosphere and vertical grid to
    ! new structures refatm_in and vcoord_in
    refatm_in%irefatm    = izrefatm_in
    refatm_in%irefatm_id = 0
    refatm_in%p0sl       = zp0sl_in
    refatm_in%t0sl       = zt0sl_in
    IF     (izrefatm_in == 1) THEN
      refatm_in%dt0lp    = zdt0lp_in
      refatm_in%delta_t  = rundefined
      refatm_in%h_scal   = rundefined
    ELSEIF (izrefatm_in == 2) THEN
      refatm_in%dt0lp    = rundefined
      refatm_in%delta_t  = zdelta_t_in
      refatm_in%h_scal   = zh_scal_in
    ENDIF
    refatm_in%bvref      = rundefined
    vcoord_in%ivctype    = izvctype_in
    vcoord_in%ivcoord_id = 0
    vcoord_in%nlevels    = ke_in + 1
    IF (izvctype_in == 1) THEN
      vcoord_in%sigm_coord(1:ke1in) = zvcoord_in(1:ke1in)
      vcoord_in%vert_coord(1:ke1in) = -1.0_ireals
    ELSE
      vcoord_in%sigm_coord(1:ke1in) = -1.0_ireals
      vcoord_in%vert_coord(1:ke1in) = zvcoord_in(1:ke1in)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
#endif

END SUBROUTINE read_nc_gdefs

!==============================================================================
!==============================================================================
!+ Module procedure in "src_read_coarse_grid" for reading a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE read_nc_vdefs (ncid, listin, nvarin, var_id, pollon, pollat, &
                          pcontrol_fi, lcheck,                          &
                          icomm, myid, npes, yerrmsg, istatus)
!
!------------------------------------------------------------------------------
!
! Description:
!   This routine reads attributes for each input variable
!   from NetCDF formatted input file.
!
!------------------------------------------------------------------------------
!
!
! Subroutine / Function arguments

! Scalar arguments with intent(in):
INTEGER (kind=iintegers),   INTENT(IN) :: &
  ncid,           & ! NetCDF file ID
  nvarin,         & ! number of variables in input list
  icomm,          & ! MPI communicator
  myid,           & ! ID of this PE in communicator icomm
  npes              ! number of PEs

REAL    (KIND=ireals),    INTENT(IN)  ::  &
  pcontrol_fi,        & ! pressure of control level for geopotential  !_br061201
  pollat, pollon        ! coordinates of the rotated north pole

! Array arguments with intent(in):
TYPE(ar_des_input)  , INTENT(IN)  ::                      &
  listin(nvarin) ! List of fields for reading in

! Scalar arguments with intent(out):
INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
  istatus       ! error index

CHARACTER (LEN= *),        INTENT(OUT)    ::  &
  yerrmsg     ! string for error messages

! Array arguments with intent(out):
INTEGER (kind=iintegers), INTENT(OUT) :: &
  var_id(nvarin)       ! NetCDF ID for each variable in the input list

! Scalar arguments with intent(inout):
LOGICAL , INTENT(INOUT)                   ::  &
  lcheck     ! check whether all necessary variables has been read

!-----------------------------------------------------------------------
!
! Local scalars:

INTEGER (KIND=iintegers) ::  &
  izmplcode           ! error status variable

INTEGER (KIND=iintegers) ::  &
  jz1,               & ! loop index
  jzfi_dimid,        & ! dimension id of FI
  jzfi_dimlen,       & ! length of dimension  of FI
  jzfi_varid           ! id of FI

CHARACTER (LEN=11) ::  &
  yzcoordinates    ! string holding names of the grid coordinates"

REAL    (KIND=ireals) ::  &
  zpollat, zpollon,   &   ! coordinates of the rotated north pole
  zfi_plev                ! pressure level FI

LOGICAL                                   ::  &
  lzcheckin_r, & ! list for checking which variables should to be read
  lzcheckin_o    ! list for checking which variables may be read

LOGICAL                                   ::  &
  lzcheckin (nvarin)      ! list for checking which variables that has been read

!- End of header
!------------------------------------------------------------------------------

#ifdef NETCDF

lcheck = .TRUE.

! Processor 0 does the job
  IF (myid == 0) THEN

    var_loop: DO jz1 = 1, nvarin

      lzcheckin_r = .FALSE.
      lzcheckin_o = .FALSE.

      ! get the variable ID
      istatus = nf90_inq_varid (ncid, TRIM(listin(jz1)%name), var_id(jz1))

      IF (istatus /= NF90_NOERR) THEN
        var_id(jz1) = -1
        IF ((var_in(jz1)%dattyp(1:1) == 'I' .AND. .NOT. lcomp_bound) .OR. &
            (var_in(jz1)%dattyp(2:2) == 'B' .AND.       lcomp_bound)) THEN
          lzcheckin_r = .TRUE.
          IF (var_in(jz1)%dattyp(3:3) == 'O') THEN
            lzcheckin_o = .TRUE.
          ENDIF
        ENDIF

        IF ( lzcheckin_r .AND. .NOT. lzcheckin_o ) lcheck = .FALSE.
        istatus = 0
        CYCLE var_loop
      ELSE
        IF ((var_in(jz1)%dattyp(1:1) == 'I' .AND. .NOT. lcomp_bound) .OR. &
            (var_in(jz1)%dattyp(2:2) == 'B' .AND.       lcomp_bound)) THEN
          lzcheckin(jz1) = .TRUE.
        ELSE
            ! we do not need this variable, even if it there
          var_id(jz1) = -1
          CYCLE var_loop
        ENDIF
      ENDIF


      ! get the grid mapping
      istatus = nf90_get_att (ncid, var_id(jz1), 'grid_mapping', grid_mapping)

      IF (istatus == NF90_enotatt) THEN

        zpollat = 90._ireals
        zpollon = 180._ireals
        istatus = 0

      ELSE

        IF (istatus /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(istatus))
          RETURN
        ENDIF
        ! check, if the grid is rotated lat/lon
        IF (grid_mapping(1:12) /= 'rotated_pole') THEN
          istatus = -1
          yerrmsg = 'No rotated grid: '//TRIM(grid_mapping)
          RETURN
        ENDIF

      ENDIF

      ! check, if U coordinate is staggered
      IF ( TRIM(listin(jz1)%name) == 'U') THEN
        istatus = nf90_get_att (ncid, var_id(jz1), 'coordinates', yzcoordinates)
        IF (istatus == NF90_enotatt) THEN
          istatus = 0
        ELSEIF (istatus /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(istatus))
          RETURN
        ENDIF
        IF ((llm2lm .AND. yzcoordinates(1:11) /= 'slonu slatu')) THEN
          IF (myid == 0) THEN
            PRINT *, 'coordinates = ',yzcoordinates
          ENDIF
          yerrmsg = ' ERROR *** Wrong coordinates for U-velocity ***'
          istatus = -1
          RETURN
        ENDIF
      ENDIF

      ! check, if V coordinate is staggered
      IF ( TRIM(listin(jz1)%name) == 'V') THEN
        istatus = nf90_get_att (ncid, var_id(jz1), 'coordinates', yzcoordinates)
        IF (istatus == NF90_enotatt) THEN
          istatus = 0
        ELSEIF (istatus /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(istatus))
          RETURN
        ENDIF
        IF ((llm2lm .AND. yzcoordinates(1:11) /= 'slonv slatv')) THEN
          IF (myid == 0) THEN
            PRINT *, 'coordinates = ',yzcoordinates
          ENDIF
          yerrmsg = ' ERROR *** Wrong coordinates for V-velocity ***'
          istatus = -1
          RETURN
        ENDIF
      ENDIF

      ! check, if pcontrol_fi is the same as pressure level of FI
      ! if not calculate FI internally later
      IF ( TRIM(listin(jz1)%name) == 'FI') THEN
        istatus = nf90_inq_dimid(ncid, 'pressure', jzfi_dimid)
        IF (istatus == NF90_NOERR) THEN
          istatus = nf90_inquire_dimension(ncid, jzfi_dimid, len=jzfi_dimlen)
        ELSE
          PRINT *, '*** Cannot get pressure id for FI pressure dimension'
          PRINT *, '*** Calculate FI internally'
        ENDIF
        IF (istatus == NF90_NOERR) THEN
          IF (jzfi_dimlen /= 1) THEN
            PRINT *, '*** pressure dimension length of FI : ', jzfi_dimlen
            yerrmsg = '*** pressure dimension of FI is not equal 1'
            PRINT *, '*** Calculate FI internally'
           istatus = -1
          ENDIF
        ENDIF
        IF (istatus == NF90_NOERR) THEN
          istatus = nf90_inq_varid(ncid, 'pressure', jzfi_varid)
          IF (istatus /= NF90_NOERR) THEN
            PRINT *, '*** Cannot get pressure level value for FI'
            PRINT *, '*** Calculate FI internally'
          ENDIF
        ENDIF
        IF (istatus == NF90_NOERR) THEN
          istatus = nf90_get_var(ncid, jzfi_varid, zfi_plev)
          IF (istatus /= NF90_NOERR) THEN
            PRINT *, '*** Cannot get pressure level value for FI'
            PRINT *, '*** Calculate FI internally'
            RETURN
          ENDIF
        ENDIF
        IF (istatus == NF90_NOERR) THEN
          IF (zfi_plev /= pcontrol_fi) THEN
            PRINT *, '*** Pressure level value for FI is different from pcontrol_fi'
            PRINT *, '*** Pressure level value for FI : ', zfi_plev
            PRINT *, '***                 pcontrol_fi : ', pcontrol_fi
            PRINT *, '*** Calculate FI internally'
            istatus = -1
          ENDIF
        ENDIF
        IF (istatus == 0) THEN
           lzcheckin(jz1) = .TRUE.
        ELSE
          lzcheckin(jz1) = .FALSE.
          istatus = 0
        ENDIF
      ENDIF

    ENDDO  var_loop

  ENDIF ! myid == 0

! send information to the other processors
  IF (npes > 1) THEN
    CALL MPI_BCAST(var_id,    nvarin, imp_integers, 0, icomm, izmplcode)
    CALL MPI_BCAST(lzcheckin, nvarin, imp_logical,  0, icomm, izmplcode)
    CALL MPI_BCAST(lcheck,        1, imp_logical,  0, icomm, izmplcode)
  ENDIF

! Set the data checking flags in the structure var_in
  DO jz1=1,nvar_in
    ! lzcheckin(jz1) might be TRUE even if a variable does not exist in coarse
    ! grid input data, e.g. QV_S, but var_id is -1 in this case.
    ! To avoid problems the check is extend and only if lzcheckin(jz1)=TRUE
    ! and var_id(jz1) > 0 the variable is marked as read

    IF (lzcheckin(jz1) .AND. (var_id(jz1) > 0)) THEN
      var_in(jz1)%lreadin = .TRUE.
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
#endif

END SUBROUTINE read_nc_vdefs

!==============================================================================

END MODULE src_read_coarse_grid
