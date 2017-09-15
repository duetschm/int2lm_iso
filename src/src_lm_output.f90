!+ Module that handles output of LM fields in Grib1 or NetCDF format.
!==============================================================================

MODULE src_lm_output

!==============================================================================
!
! Description:
!   This module contains subroutines necessary for writing the result data
!   for the LM in Grib1 or NetCDF format.
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
!  Grib Adaptations for gds_out (=8) and length of Grib records in bytes
!  Adaptations to interface changes in io_utilities, parallel_utilities
!  Adaptations for output of chemistry variables
!  Replaced ke_soil to ke_soil_lm
!  Replaced gather_grib to gather_values
! V1_6         2007/09/07 Ulrich Schaettler, Burkhardt Rockel
!  Introduction of NetCDF Output
! V1_7         2007/11/26 Ulrich Schaettler, Christoph Gebhardt
!  Introduced ensemble mode for computing boundary conditions
!  Bug correction in write_nc_vdefs: kedim has to be used for 3rd dimension
!  Introduced undefsub to pass to SR check_record
! V1_8         2008/05/29 Ulrich Schaettler, Hans-Juergen Panitz
!  Moved undef-settings for sea-points in NetCDF to output field to avoid
!  inconsistencies with different real-formats
!  Adaptations for checking undef settings
! V1_9         2009/09/03 Guenther Zaengl, Burkhardt Rockel, Oli Fuhrer
!  Adaptations for new reference atmosphere
!  Replace ldwd_grib_use with l_ke_in_gds
!  Adaptations for new reference atmosphere in netCDF output (Burkhardt Rockel)
!  Changed error code for write_nc_gdefs
!  Implemented output of 3D external parameter HORIZON (Oli Fuhrer)
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler, Burkhardt Rockel, Anne Roches
!  Replaced iednr by idwdednr for grbex1 interface
!  Corrected attribute "flag_values" in netCDF output for soil type
!  Added parameter 'i' for inland water (lakes) (BR)
!  Modifications for output of external parameter HORIZON in NetCDF (AR)
! V1_15        2010/12/10 Ulrich Schaettler
!  For the new reference atmosphere, also set the SLEVE coordinate parameters
!  (=0, if SLEVE is not used, but they must not be undefined)
!  For EPS mode, ipds(51) (number of members in ensemble) was set to 0; this
!   has been removed
! V1_19        2012/06/06 Ulrich Schaettler, Dmitrii Mironov, Susanne Brienen
!                         Burkhardt Rockel
!  Implemented conditional compilation for NETCDF, GRIBDWD
!  Treatment of FLake variables, especially for GRIB output
!  Replace msoilgrib by msoilgrib_lm for outgoing data
!  Added 365_days support in NetCDF routines
!  Unified dimension IDs with COSMO (ID for topo corrections changed from 14 to 15)
! V1_20        2012/09/03 Ulrich Schaettler, Burkhardt Rockel
!  Enlarged strings for date variables to 14 characters
!  Adapted calls to subroutine make_fn
!  Burkhardt Rockel:
!  Introduction of additional global attributes in case of netCDF output
!    which can be set via namelist (analog the attributes in the namelist in COSMO)
!  Several corrections regarding netCDF output
! V1_21        2013/03/25 Ulrich Schaettler, Burkhardt Rockel
!  Implementation of grib_api for writing grib files (US)
!  Changes in netCDF output: scalar variable instead of extra dimension of length 1 
!    (e.g. height_2m) needs changes in idims_id field indices (Burkhardt Rockel)
! V1_22        2013/07/11 Ulrich Schaettler, KIT
!  Adapted settings of some GRIB2 meta data
!  Use pv_out, inrvert_out from data_int2lm_io
!  Implemented GRIB2 output for general vertical coordinate
!  Use new structures refatm_out, vcoord_out from vgrid_refatm_utils
!  Implemented possibility to write HHL file: org_lm_output is called with a 
!   special list then. Then also a UUID is created and written to this file
!  Introduced possibility to gather double precision fields for output of HHL
!  Renamed grib buffers: ds_grib to ds_grib_single, ds_gribapi to ds_grib_double
!  Use yextension and pass it to subroutine make_fn as argument (KIT)
!  Adapted interface to grib_api routines with special grib_api integer
! V1_23        2013/10/02 Ulrich Schaettler
!  Treat pv_out as double precision, just as pv_in
!  Rename vcoord_out, refatm_out to vcoord, refatm (as is in COSMO-Model)
! V1_24        2013/11/01 Ulrich Schaettler
!  Set vertical coordinate parameters in pv_out depending on vertical coordinate type
!  in init_lm_output
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Modules used:
USE data_parameters , ONLY :  &
  ireals,            & ! KIND-type parameters for real variables
  iintegers,         & ! KIND-type parameter for "normal" integer variables
  irealgrib,         & ! KIND-type of the REALs in the grib library
  iwlength,          & ! length of an integer word in byte
  intgribf,          & ! KIND-type of the fortran decks in the grib library
  intgribc,          & ! KIND-type of the c decks in the grib library
  int_ga               ! integer precision for grib_api: length of message in bytes

!------------------------------------------------------------------------------

#ifdef GRIBAPI
! grib_api interface
USE grib_api
#endif

!------------------------------------------------------------------------------

USE data_int2lm_io,            ONLY:    &
  nunit_of_time,     & ! indicator for unit-of-time (1hr, 15min, 30min,...)
  ylm_cat,           & ! catalog-name of the LM files
  ylm_lfn,           & ! name of the file with LM input data
  ylm_hhl,           & ! name of the file LM HHL fields
  nuchkdat,          & ! checking the I/O data
  yuchkdat,          & ! checking the I/O data
  ymode_write,       & ! mode for opening the (write) Grib files
  ylm_form_write,    & ! output format of LM data
  yinput_type,       & ! type of input data: 'forecast', 'analysis' or 'ana_init'
  npds,              & ! Dimension for product definition section (pds)
  ngds,              & ! Dimension for grid description section (gds)
  nbms,              & ! Dimension for bit map section (bms)
  nbds,              & ! Dimension for binary data section
  ndsup,             & ! Dimension for dsup
  ndims,             & ! Dimension for idims (contains all dimensions)
  igrib1_sample,        igrib2_sample,          & ! common samples
  igrib1_hybrid,        igrib2_hybrid,          & ! full level data (U, V, T, etc)
  igrib1_hybridlayer,   igrib2_hybridlayer,     & ! half level data (W, HHL)
  igrib1_surface,       igrib2_surface,         & ! surface data
  igrib1_depthbelow,    igrib2_depthbelow,      & ! soil data
  lfd, lfa, lds,     & !
  nbitmap,           & !
  nrbit,             & ! packrate
  nvar_lm,           & ! maximum number of variables in LM variable table
  nvar_in,           & ! maximum number of variables in input variable table
  nvar_in_norm,      & ! maximum number of variables in input variable table
  nvar_in_chem,      & ! maximum number of variables in input variable table
  nvar_lm_chem,      & ! maximum number of variables in LM variable table
  var_in,            & ! variable table for input fields
  ylevltypes_out,    & ! to convert GRIB1 level types to grib_api string typeOfLevel
                       ! for outgoing COSMO data
  idwdednr,          & ! grib edition number for dwdlib
  undefgrib,         & ! value for "undefined" in the grib routines
  undefncdf,         & ! value for "undefined" in the ncdf routines
  undef                ! value for "undefined" in the grib routines

USE data_int2lm_io,            ONLY:    &
  nlocaldefnr,       & ! local definition number for GRIB2 local section template
  nactlocdefnr,      & ! to overwrite Namelist parameter with some center default
  nprocess_ini,      & ! type of database for LM initial data
  nprocess_bd,       & ! type of database for LM boundary data
  ncenter,           & ! originating center identification
  nsubcenter,        & ! originating center identification
  ydate_ini,         & ! start of the forecast yyyymmddhh (year,month,day,hour)
  lmmss_ini,         & ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                       ! if .FALSE. 10 digits date format (YYYYMMDDHH)
                       ! for ydate_ini and result files of INT2LM
  ylm_form_write,    & ! output format of LM data
  ytunit_out,        & ! time unit for output data
  idims_id,          & ! array for the IDs of the dimensions
  iblock,            & ! array for gribed data
  ymessage,          & ! char array for grib_api
  idims_out,         & ! array for all dimensions
  ibmap,             & ! array for
  ipds,              & ! product definition section
  igds_out,          & ! grid description section
  ibms,              & ! bit map section
  ibds,              & ! binary data section
  inrvert_out,       & ! number of vertical coordinate parameters
  pv_out,            & ! array for vertical coordinate parameters
  dsup,              & ! Parameter for grib routines
  ds_grib_single,    & ! array for unpacked data
  ds_grib_double,    & ! array for unpacked grib_api data
  var_lm,            & ! array for LM variable table
  l_ke_in_gds,       & ! explicit GDS entry for number of model levels
  lchkout,           & ! logical for print of check-values (max,min,mean)
                       ! of LM/HM-fields
  ytrans_out,        & ! directory for writing ready-files
  yncglob_institution,   & ! originating center name
  yncglob_title,         & ! title string for the output
  yncglob_source,        & ! program name and version
  yncglob_project_id,    & ! identification of the project of simulation
  yncglob_experiment_id, & ! identification of the experiment of simulation
  yncglob_contact,       & ! contact e.g. email address
  yncglob_references,    & ! URL, report etc.
  ncglob_realization,    & ! number of the realization of the experiment
  yextension

!------------------------------------------------------------------------------

USE data_int2lm_control,       ONLY:    &
    nvers,           & ! for documenting purposes in Grib Code
    nl_soil_lm,      & ! number of soil layers in LM, resp. HM
    lmulti_layer_lm, & ! if .TRUE., compute soil fields for multi-layer soil
    dt,              & ! time step used in the LM
    itype_calendar,  & ! for specifying the calendar used
    lcomp_bound,     & ! compute fields for boundaries
    leps_bc,         & ! switch ensemble mode for boundary conditions on/off
    lroutine,        & ! if .TRUE., operational-job
    msoilgrib_lm,    & ! grib coded depth of main soil levels in centimeters
    idbg_level,      & ! to control verbosity of output
    nhori,           & ! number of sectors for the horizont array by the topographic
                       ! correction of the radiation
    lradtopo,        & ! if .TRUE., calculate topographic correction of radiation
    lprintdeb_all      ! whether all PEs print debug output

!------------------------------------------------------------------------------

USE data_grid_lm,       ONLY:    &
    pollat,      & ! latitude of the rotated north pole (in degrees, N>0)
    pollon,      & ! longitude of the rotated north pole (in degrees, E>0)
    polgam,      & ! latitude of the rotated north pole
    dlat,        & ! grid point distance in zonal direction (in degrees)
    dlon,        & ! grid point distance in meridional direction (in degrees)
    startlat_tot,& ! transformed latitude of the lower left grid point
                   ! of the total domain (in degrees, N>0)
    startlon_tot,& ! transformed longitude of the lower left grid point
                   ! of the total domain (in degrees, E>0)
    endlat_tot,  & ! transformed latitude of the upper right grid point
                   ! of the total domain (in degrees, N>0)
    endlon_tot,  & ! transformed longitude of the upper right grid point
                   ! of the total domain (in degrees, E>0)
    ielm_tot,    & ! ie for LM, total domain
    jelm_tot,    & ! je for LM, total domain
    kelm_tot,    & ! ke for LM
    kedim,       & !
    ke_soil_lm,  & ! number of layers in multi-layer soil model
    czmls_lm,    & ! depth of the main soil layers in meters in output
    czhls_lm,    & ! depth of the half soil layers in meters in output
    ie2lm_max,   & ! Max. of ie2lm on all processors
    je2lm_max,   & ! Max. of je2lm on all processors
    ie2lm_tot,   & ! = ielm_tot + 2
    je2lm_tot,   & ! = jelm_tot + 2
    kelm,        & !
    ke1lm,       & !
    ie2lm,       & !
    je2lm          !

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY:    &
    lasync_io,       & ! if .TRUE.: the model runs with extra PEs for
                       ! asynchronous IO
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    isubpos,         & ! positions of the subdomains in the total domain
    icomm_cart,      & ! communicator for the virtual cartesian topology
    imp_grib,        & ! determines the REAL type used for the GRIB library
    imp_reals,       & ! determines the REAL type used in the model
    imp_character,   & ! determines the correct CHARACTER type used in the model
                       ! for MPI
    imp_integers       ! determines the correct INTEGER type used in the
                       ! model for MPI

!------------------------------------------------------------------------------

USE data_epsctl,        ONLY:    &
    iepsmem_bc,    & ! ID of the member in the ensemble of boundary
                     ! conditions (iepsmem_bc >= 0)
    iepstyp_bc,    & ! ID of the boundary ensemble generation type
                     ! (iepstyp_bc >= 0)
    iepstot_bc       ! total number of boundary ensemble members (iepstot_bc>= 0)

!------------------------------------------------------------------------------

USE data_fields_lm,            ONLY:    &
    fr_lake_lm,  & ! lake fraction of grid element                    (  1  )
    depth_lk_lm, & ! depth of lakes
    lolp_lm,     & ! Land Sea Mask of LM for 'M'atch Interpolation
    w_snow_lm      ! water content of the snow

!------------------------------------------------------------------------------

USE environment,        ONLY:   model_abort

!------------------------------------------------------------------------------

USE parallel_utilities, ONLY:   gather_values, remark, distribute_values

!------------------------------------------------------------------------------

USE io_utilities,       ONLY :  open_file, write_grib, write_gribapi, write_netcdf,     &
                                close_file, check_record, make_fn

!------------------------------------------------------------------------------

USE utilities,          ONLY :  rlarot2rla, phirot2phi

!------------------------------------------------------------------------------

USE vgrid_refatm_utils, ONLY:   &
    lnewVGrid,                                                           &
    nfltvc, svc1, svc2, refatm, vcoord,                                  &
    uuid_out, uuid_out_string,                                           &
    uuid_create, uuid_2char, uuid_string_length

!------------------------------------------------------------------------------

#ifdef NETCDF
USE netcdf,           ONLY :   &
  nf90_def_dim,            &
  nf90_def_var,            &
  nf90_enddef,             &
  nf90_redef,              &
  nf90_put_att,            &
  nf90_put_var,            &
  nf90_noerr,              &
  nf90_strerror,           &
  NF90_CHAR,               &
  NF90_DOUBLE,             &
  NF90_FLOAT,              &
  NF90_GLOBAL,             &
  NF90_UNLIMITED
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Arrays used during output to store fields between subsequent calls to 
! output fields
REAL (KIND=irealgrib), ALLOCATABLE        ::   &
  procarray_single(:,:,:)

REAL (KIND=ireals),    ALLOCATABLE        ::   &
  procarray_double(:,:,:)

! string variable to hold grid information
  CHARACTER (LEN=200)     grid_mapping

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in src_lm_output for initializing the output organization
!------------------------------------------------------------------------------

SUBROUTINE init_lm_output

!------------------------------------------------------------------------------
!
! Description:
!  The routine init_lm_output initializes organizational variables of the 
!  output routines dealing with grib 1 code (such as the dimensions of
!  the different grib code sections). Also the grid description section
!  is initialized (except the location of the lower left grid point, because
!  it depends on the variable (U,V or other)).
!  
! Method:
!
!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers)                ::   &
  izerrstat, nziostat, k, nzbyte, lfdlm, ldslm, ndate, ntime, nsecond,    &
  nstatus, modnvers, izdebug

REAL (KIND=ireals)                      ::   &
  pollon_sp, pollat_sp

LOGICAL                                 ::   &
 lzopen            ! to inquire YUCHKDAT

CHARACTER  (LEN=25)                     ::   &
  yzroutine

CHARACTER  (LEN=80)                     ::   &
  yzerrmsg

#ifdef GRIBDWD
INTEGER (KIND=intgribf), EXTERNAL  :: IREFTS
#endif
!
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Section 1: Initializations
!------------------------------------------------------------------------------

  yzroutine = 'init_lm_output'
  izerrstat = 0_iintegers

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  IF (izdebug > 5) THEN
    PRINT *, ' Initialize output'
  ENDIF

  ! Set dimensions for grib variables. The allocation of iblock, ibmap, ds 
  ! and dsup is done in org_lm_output.
  ! nzbyte is assumed to be 2 here: this is true in all our cases, but if
  ! if we once change the packing rate of grib code (nrbit=16) we will get
  ! in trouble here. To know this a priori, the pds of the first record has
  ! to be read and decoded to get nrbit. The 2000 are just a safety-add to
  ! take care of the definition sections that are also stored in iblock
  nzbyte = NINT (REAL(nrbit / 8, ireals))        ! usually 2 bytes
  ldslm   = ielm_tot * jelm_tot
  lfdlm  = ldslm  * nzbyte / iwlength + 2000
  lds = INT (ldslm, intgribf)
  lfd = INT (lfdlm, intgribf)

  ! Initializations for the grib library
  !  moving arraydimensions into idims
  !  declaration dimensions
  idims_out( 1) = npds
  idims_out( 2) = ngds
  idims_out( 3) = nbms
  idims_out( 4) = nbds
  idims_out( 5) = nbitmap
  idims_out( 6) = ndsup
  idims_out( 7) = lds
  idims_out( 8) = lfd

  !  real dimensions
  idims_out(11) = 47
  IF (.NOT.l_ke_in_gds) THEN
    ! old style of coding the vertical coordinate parameters
    IF (vcoord%ivctype == 3 .OR. vcoord%ivctype == 4) THEN   ! SLEVE coordinates
      idims_out(12) = 25 + 4 + ke1lm + 4
    ELSE                     ! NON SLEVE coordinates
      idims_out(12) = 25 + 4 + ke1lm + 1
    ENDIF
    IF (refatm%irefatm == 2) idims_out(12) = 25 + 4 + ke1lm + 6
  ELSE
    ! new style of coding the vertical coordinate parameters
    IF (vcoord%ivctype == 3 .OR. vcoord%ivctype == 4) THEN   ! SLEVE coordinates
      idims_out(12) = 25 + 6 + ke1lm + 3
    ELSE                     ! NON SLEVE coordinates
      idims_out(12) = 25 + 6 + ke1lm
    ENDIF
    IF (refatm%irefatm == 2) idims_out(12) = 25 + 6 + ke1lm + 5
  ENDIF

  idims_out(13) = 3
  idims_out(14) = 5
  idims_out(15) = ielm_tot*jelm_tot
  idims_out(16) = 0
  idims_out(17) = ielm_tot*jelm_tot


  ! gridpoints, simple packing, floating point data
  ibds(2)   = 0

  ! nrbit, number of bits
  ibds(5)   = nrbit

  ! no bitmap
  ibms(3)   = -2

!------------------------------------------------------------------------------
! Section 2: Preparation for output of vertical coordinate parameters
!------------------------------------------------------------------------------

  ! number of the vertical coordinate parameters
  IF (.NOT.l_ke_in_gds) THEN
    ! old style of coding the vertical coordinate parameters
    IF (vcoord%ivctype == 3 .OR. vcoord%ivctype == 4) THEN   ! SLEVE coordinates
      inrvert_out = 8 + ke1lm
    ELSE                     ! NON SLEVE coordinates
      inrvert_out = 5 + ke1lm
    ENDIF
    IF (refatm%irefatm == 2) inrvert_out = 10 + ke1lm
  ELSE
    ! new style of coding the vertical coordinate parameters
    IF (vcoord%ivctype == 3 .OR. vcoord%ivctype == 4) THEN   ! SLEVE coordinates
      inrvert_out = 9 + ke1lm
    ELSE                     ! NON SLEVE coordinates
      inrvert_out = 6 + ke1lm
    ENDIF
    IF (refatm%irefatm == 2) inrvert_out = 11 + ke1lm
  ENDIF

  ALLOCATE (pv_out(inrvert_out),   STAT = izerrstat)

  IF (.NOT. l_ke_in_gds) THEN
    ! old style of coding: first 4 values for the reference atmosphere,
    ! then the vertical coordinate parameters and last eventually
    ! additional parameters for the SLEVE coordinate
    pv_out( 1) = refatm%p0sl
    pv_out( 2) = refatm%t0sl
    pv_out( 3) = refatm%dt0lp
    pv_out( 4) = vcoord%vcflat

    IF     (vcoord%ivctype == 1) THEN
      DO k = 1, ke1lm
        pv_out (4 + k) = vcoord%sigm_coord(k)
      ENDDO
    ELSE
      DO k = 1, ke1lm
        pv_out (4 + k) = vcoord%vert_coord(k)
      ENDDO
    ENDIF

    ! But this must be coded as integer later on because of 
    ! backwards compatibility
    pv_out (4+ke1lm+1) = vcoord%ivctype

    IF (vcoord%ivctype == 3 .OR. vcoord%ivctype == 4) THEN
      pv_out (4+ke1lm+2) = svc1
      pv_out (4+ke1lm+3) = svc2
      pv_out (4+ke1lm+4) = nfltvc
    ENDIF

    IF (refatm%irefatm == 2) THEN ! Write parameters for new reference atmosphere
      pv_out (4+ke1lm+1) = vcoord%ivctype+100
      IF (vcoord%ivctype /= 3 .AND. vcoord%ivctype /= 4) THEN 
        pv_out (4+ke1lm+2) = 0.0_ireals
        pv_out (4+ke1lm+3) = 0.0_ireals
        pv_out (4+ke1lm+4) = 0.0_ireals
      ENDIF
      pv_out (4+ke1lm+5) = refatm%delta_t
      pv_out (4+ke1lm+6) = refatm%h_scal
    ENDIF
  ELSE
    ! new style of coding: the first 2 values are the vertical coordinate
    ! type and the number of main levels used. Then 4 values for the
    ! reference atmosphere, the vertical coordinate parameters and
    ! last eventually additional parameters for the SLEVE coordinate
    IF (refatm%irefatm == 1) THEN
      pv_out ( 1) = REAL (vcoord%ivctype,     ireals)
    ELSE IF (refatm%irefatm == 2) THEN
      pv_out ( 1) = REAL (vcoord%ivctype+100, ireals)
    ENDIF

    pv_out ( 2) = REAL (kelm,   ireals)
    pv_out ( 3) = refatm%p0sl
    pv_out ( 4) = refatm%t0sl
    pv_out ( 5) = refatm%dt0lp
    pv_out ( 6) = vcoord%vcflat

    IF     (vcoord%ivctype == 1) THEN
      DO k = 1, ke1lm
        pv_out (6 + k) = vcoord%sigm_coord(k)
      ENDDO
    ELSE
      DO k = 1, ke1lm
        pv_out (6 + k) = vcoord%vert_coord(k)
      ENDDO
    ENDIF

    IF (vcoord%ivctype == 3 .OR. vcoord%ivctype ==4) THEN
      ! SLEVE coordinate: 3 new parameters are written to GDS
      pv_out (6+ke1lm+1) = svc1
      pv_out (6+ke1lm+2) = svc2
      pv_out (6+ke1lm+3) = nfltvc
    ENDIF

    IF (refatm%irefatm == 2) THEN ! Write parameters for new reference atmosphere
      IF (vcoord%ivctype /= 3 .AND. vcoord%ivctype /=4) THEN
        pv_out (6+ke1lm+1) = 0.0_ireals
        pv_out (6+ke1lm+2) = 0.0_ireals
        pv_out (6+ke1lm+3) = 0.0_ireals
      ENDIF
      pv_out (6+ke1lm+4) = refatm%delta_t
      pv_out (6+ke1lm+5) = refatm%h_scal
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Set array igds_out for the GRIB1 grid description section
!------------------------------------------------------------------------------

#if defined(GRIBDWD) || defined(GRIBAPI)
IF (ylm_form_write(4:4) == '1') THEN
  ! set all values of igds_out; they can be used for DWDLIB or grib_api
  igds_out(:) = -999999

! igds_out(1) = 42 + (inrvert_out) * 4

! ! number of the vertical coordinate parameters
! igds_out(2) = inrvert_out

! ! location of the list of vertical coordinate parameters in bytes
! igds_out(3) = 43

  ! data representation type 
  igds_out(4) = 10

  ! calculation of the left bottom and right upper corner in millidegrees:
  ! igds_out(7, 8, 10, 11)
  ! This depends on the variable (U- and V-related variables are shifted in
  ! the Arakawa C-grid) and is determined during the output step.

  ! number of gridpoints 
  igds_out( 5) = ielm_tot
  igds_out( 6) = jelm_tot
  igds_out( 9) = 8    ! that was 0 before; but with "8", also bit 5 is set to 1

  ! increments: are set to 0 explicitly, because of igds_out(9) = 0
  igds_out(12) = 0
  igds_out(13) = 0


  igds_out(14) = 64
  igds_out(15:19) = 0

  ! coordinates of the pole: in the grib code the southern pole has to be
  ! specified: for the latitude this is the negative value, for the
  ! longitude it is the value + 180.0 and then limiting the value again
  ! to -180.0 ... 180.0
  pollon_sp = pollon + 180.0_ireals
  IF (pollon_sp > 180.0_ireals) THEN
    pollon_sp = pollon_sp - 360.0_ireals
  ENDIF

  igds_out(20) = NINT(-pollat      * 1000.0_ireals)
  igds_out(21) = NINT( pollon_sp   * 1000.0_ireals)
#ifdef GRIBDWD
  !US igds_out(22) = IREFTS(0.0_irealgrib)
  igds_out(22) = IREFTS(REAL(polgam, irealgrib))
#endif

! ! vertical coordinate parameters
! DO k = 1, inrvert_out
!   igds_out(25 + k) = IREFTS( pv_out(k) )
! ENDDO

ENDIF
#endif

!------------------------------------------------------------------------------
! Section 4: Get grib samples in case of grib_api write
!------------------------------------------------------------------------------

#ifdef GRIBAPI
IF     (ylm_form_write == 'api1') THEN

  !------------------------------------------------------------------------------
  ! Section 4.1: grib samples for GRIB1
  !------------------------------------------------------------------------------

  CALL grib_new_from_samples (igrib1_hybridlayer, 'DWD_rotated_ll_7km_HI_grib1', &
                              izerrstat)
  IF (izerrstat /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_new_from_sample: hybridlayer 1 ', izerrstat
  ENDIF

  CALL grib_new_from_samples (igrib1_hybrid,      'DWD_rotated_ll_7km_H_grib1',  &
                              izerrstat)
  IF (izerrstat /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_new_from_sample: hybrid 1 ', izerrstat
  ENDIF

  CALL grib_new_from_samples (igrib1_surface,     'DWD_rotated_ll_7km_G_grib1',  &
                              izerrstat)
  IF (izerrstat /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_new_from_sample: surface 1 ', izerrstat
  ENDIF

  CALL grib_new_from_samples (igrib1_depthbelow,  'DWD_rotated_ll_7km_S_grib1',  &
                              izerrstat)
  IF (izerrstat /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_new_from_sample: depthbelowland 1 ', izerrstat
  ENDIF

  ! clone a sample to igrib1_sample
  CALL grib_clone(igrib1_hybridlayer, igrib1_sample, izerrstat)
  IF (izerrstat /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_clone: from sample 1 ', izerrstat
  ENDIF

ELSEIF (ylm_form_write == 'api2') THEN

  !------------------------------------------------------------------------------
  ! Section 4.2: grib samples for GRIB2
  !------------------------------------------------------------------------------

  CALL grib_new_from_samples (igrib2_hybridlayer, 'DWD_rotated_ll_7km_HI_grib2', &
                              izerrstat)
  IF (izerrstat /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_new_from_sample: hybridlayer 2 ', izerrstat
  ENDIF

  CALL grib_new_from_samples (igrib2_hybrid,      'DWD_rotated_ll_7km_H_grib2',  &
                              izerrstat)
  IF (izerrstat /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_new_from_sample: hybrid 2 ', izerrstat
  ENDIF

  CALL grib_new_from_samples (igrib2_surface,     'DWD_rotated_ll_7km_G_grib2',  &
                              izerrstat)
  IF (izerrstat /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_new_from_sample: surface 2 ', izerrstat
  ENDIF

  CALL grib_new_from_samples (igrib2_depthbelow,  'DWD_rotated_ll_7km_S_grib2',  &
                              izerrstat)
  IF (izerrstat /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_new_from_sample: depthbelowland 2 ', izerrstat
  ENDIF

  ! clone a sample to igrib2_sample
  CALL grib_clone(igrib2_surface, igrib2_sample, izerrstat)
  IF (izerrstat /= GRIB_SUCCESS) THEN
    PRINT *,   ' *** Error in grib_clone: from sample 2 ', izerrstat
  ENDIF

ENDIF

IF (izerrstat /= GRIB_SUCCESS) THEN
  yzerrmsg = 'ERROR in grib_new_from_sample'
  CALL model_abort (my_cart_id, izerrstat, yzerrmsg, yzroutine)
ENDIF

!------------------------------------------------------------------------------
! Section 5: Set constant GRIB meta data for grib_api
!------------------------------------------------------------------------------

IF     (ylm_form_write == 'api1') THEN

  !------------------------------------------------------------------------------
  ! Section 5.1: GRIB1 Grid Description Section
  !------------------------------------------------------------------------------

  ! here we set the constant gds- and pds-values only for one sample

  ! location of the list of vertical coordinate parameters in bytes
  ! igds_out(3) = 43         ! hopefully not with grib_api

  ! data representation type 
  ! igds_out(4) = 10
  CALL grib_set (igrib1_sample,'dataRepresentationType', igds_out(4))

  ! calculation of the left bottom and right upper corner in millidegrees:
  ! igds_out(7, 8, 10, 11)
  ! This depends on the variable (U- and V-related variables are shifted in
  ! the Arakawa C-grid) and is determined during the output step.

  ! number of gridpoints 
  CALL grib_set (igrib1_sample,'Ni', igds_out(5))            ! igds_out( 5) = ielm_tot
  CALL grib_set (igrib1_sample,'Nj', igds_out(6))            ! igds_out( 6) = jelm_tot
  CALL grib_set (igrib1_sample,'resolutionAndComponentFlags', igds_out(9))
  ! igds_out( 9) = 8    ! that was 0 before; but with "8", also bit 5 is set to 1

  ! increments: are set to 0 explicitly, because of igds_out(9) = 8
  CALL grib_set (igrib1_sample,'Di', igds_out(12))           ! igds_out(12) = 0
  CALL grib_set (igrib1_sample,'Dj', igds_out(12))           ! igds_out(13) = 0


  CALL grib_set (igrib1_sample,'scanningMode', igds_out(14)) ! igds_out(14) = 64
  ! igds_out(15:19) = 0

  ! coordinates of the pole: in the grib code the southern pole has to be
  ! specified: for the latitude this is the negative value, for the
  ! longitude it is the value + 180.0 and then limiting the value again
  ! to -180.0 ... 180.0
  pollon_sp = pollon + 180.0_ireals
  IF (pollon_sp > 180.0_ireals) THEN
    pollon_sp = pollon_sp - 360.0_ireals
  ENDIF

  CALL grib_set (igrib1_sample, 'latitudeOfSouthernPole',  igds_out(20))
  ! igds_out(20) = NINT(-pollat      * 1000.0_ireals)
  CALL grib_set (igrib1_sample, 'longitudeOfSouthernPole', igds_out(21))
  ! igds_out(21) = NINT( pollon_sp   * 1000.0_ireals)
  CALL grib_set (igrib1_sample, 'angleOfRotation', polgam)
  ! igds_out(22) = IREFTS(REAL(polgam, irealgrib))
  !US igds_out(22) = IREFTS(0.0_irealgrib)

  ! vertical coordinate parameters
  ! are set per record

ELSEIF (ylm_form_write == 'api2') THEN

  !------------------------------------------------------------------------------
  ! Section 5.2: Constant GRIB2 meta data
  !------------------------------------------------------------------------------

  ! Indicator Section
  ! -----------------

  CALL grib_set (igrib2_sample, 'centre',                      ncenter)   ! originating centre
  CALL grib_set (igrib2_sample, 'subCentre',                nsubcenter)   ! originating subcentre

  ! Identification Section
  ! ----------------------

  CALL grib_set (igrib2_sample, 'significanceOfReferenceTime', 1)         ! start of forecast

  ! Reference time of data
  READ(ydate_ini( 1: 8),'(I8)') ndate
  READ(ydate_ini( 9:12),'(I4)') ntime
  READ(ydate_ini(13:14),'(I2)') nsecond
  CALL grib_set (igrib2_sample,'dataDate',                     ndate)     ! yyyymmdd
  CALL grib_set (igrib2_sample,'dataTime',                     ntime)     ! hhmm
  CALL grib_set (igrib2_sample,'second',                     nsecond)     ! ss

  ! Production Status: depends on the center and is computed in the local use section

  ! Type of processed data
  IF (leps_bc) THEN
    ! this is set to "control and perturbed forecast products"
    CALL grib_set (igrib2_sample,'typeOfProcessedData',             5)
  ELSE
    !! then it depends on the type of the driving model, 
    !! which is specified in yinput_type:
    !IF       (yinput_type == 'analysis') THEN
    !  CALL grib_set (igrib2_sample,'typeOfProcessedData',             0)  ! Analysis
    !ELSEIF ( (yinput_type == 'forecast') .OR. (yinput_type == 'ana_init') ) THEN
    !  CALL grib_set (igrib2_sample,'typeOfProcessedData',             1)  ! Forecast
    !ENDIF
    CALL grib_set (igrib2_sample,'typeOfProcessedData',           255)
  ENDIF

  ! Local Use Section
  ! -----------------

  SELECT CASE (ncenter)
  CASE (78)
    ! Production Status: check nvers:
    modnvers = MODULO(nvers, 16384)

    SELECT CASE (modnvers)
    CASE (  1: 50)
      nstatus = 0    ! Operational
    CASE ( 51: 99)
      nstatus = 1    ! Operational Test
    CASE (100:16384)
      nstatus = 2    ! Research
    END SELECT

    CALL grib_set (igrib2_sample, 'productionStatusOfProcessedData', nstatus)

    ! local number of experiment: set "pure" experiment number from nvers
    CALL grib_set (igrib2_sample, 'localNumberOfExperiment',        modnvers)  ! GRIB1: ipds(47)

    ! local information number: not needed in INT2LM products, but for HHL-file (set later)
    CALL grib_set (igrib2_sample, 'localInformationNumber',                0)  !  ipds(41)

    SELECT CASE (nlocaldefnr)
    CASE (-1)
      ! local definition number not set per namelist: use DWD defaults:
      CALL grib_set (igrib2_sample,'localDefinitionNumber',              254)  ! deterministic system
      nactlocdefnr = 254

      IF (leps_bc) THEN
        IF (lcomp_bound .AND. iepstyp_bc >= 0) THEN
          CALL grib_set (igrib2_sample, 'localDefinitionNumber',         253)  ! ensemble BC
          CALL grib_set (igrib2_sample, 'localTypeOfEnsembleForecast', iepstyp_bc)
          nactlocdefnr = 253
        ENDIF
      ENDIF
    CASE (250,252,253,254)
      ! in principle this is the DWD default
      CALL grib_set (igrib2_sample,'localDefinitionNumber', nlocaldefnr)
      nactlocdefnr = nlocaldefnr
    CASE DEFAULT
      PRINT *, ' *** ERROR in localDefinitionNumber *** ', nlocaldefnr
      PRINT *, ' ***       This is not a valid definition number! *** '
      yzerrmsg = 'Error in settings for GRIB2 local use section'
      CALL model_abort (my_cart_id, 2022, yzerrmsg, yzroutine)
    END SELECT

  CASE DEFAULT

    ! Production Status
    IF (lroutine) THEN
      CALL grib_set (igrib2_sample,'productionStatusOfProcessedData', 0)  ! Operational
    ELSE
      CALL grib_set (igrib2_sample,'productionStatusOfProcessedData', 2)  ! Research
    ENDIF

    nactlocdefnr = nlocaldefnr ! for later use

    ! Local Use Section
    IF (nlocaldefnr /= -1) THEN
      CALL grib_set (igrib2_sample,'localDefinitionNumber', nlocaldefnr)

      SELECT CASE (nlocaldefnr)
      CASE (250,252,253,254)
        CALL grib_set (igrib2_sample, 'localNumberOfExperiment',   nvers)
      END SELECT

    ENDIF

  END SELECT


  ! Grid Definition Section
  ! -----------------------

  ! Most keys are already set in the sample grib message

  ! Number of points Ni, Nj
  CALL grib_set (igrib2_sample,'Ni',                        ielm_tot)
  CALL grib_set (igrib2_sample,'Nj',                        jelm_tot)

  ! start and end points have to be set per product because of the staggered grid

  ! resolution and component flags
  ! the data in the sample are set to:
  !   Bit 3 = 1: i direction increments given
  !   Bit 4 = 1: j direction increments given
  !   Bit 5 = 1: Resolved u- and v-components of vector quantities relative to the defined grid
  !     value of octet 55 therefore is: 56
  CALL grib_set (igrib2_sample,'ijDirectionIncrementGiven',       1)
  CALL grib_set (igrib2_sample,'uvRelativeToGrid',                1)
  CALL grib_set (igrib2_sample,'iDirectionIncrementInDegrees', dlon)
  CALL grib_set (igrib2_sample,'jDirectionIncrementInDegrees', dlat)

  ! scanning mode is set in the sample

  ! specifications of rotated pole
  ! convert COSMO north pole to rotated south pole
  pollat_sp = -pollat
  ! pollon is specified in the range -180.0...+180.0
  ! but GRIB2 gives values in the range 0.0...+360.0
  pollon_sp =  pollon + 180.0_ireals     
  ! above statement converts to south pole and automatically to range 0.0...+360.0
  ! so no additional correction as in GRIB1 is necessary

  CALL grib_set (igrib2_sample, 'latitudeOfSouthernPoleInDegrees',  pollat_sp)
  CALL grib_set (igrib2_sample, 'longitudeOfSouthernPoleInDegrees', pollon_sp)
  CALL grib_set (igrib2_sample, 'angleOfRotationInDegrees',         polgam)

ENDIF
#endif

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE init_lm_output

!==============================================================================
!+ Module procedure in src_lm_output for organizing grib 1 output
!------------------------------------------------------------------------------

SUBROUTINE org_lm_output (numlist, ylist, nnow, ydate1, lwrite_hhl)

!------------------------------------------------------------------------------
!
! Description:
!  The routine org_lm_output organizes the output of all variables contained
!  in ylist. Parallelization for the output is by layers that should be written
!  to the grib file. Every PE gets a layer and packs it into grib format.
!  If the logical lwrite_hhl is true, a special file containing the "height of
!  half" levels for GRIB2, general vertical coordinate, is written.
!
! Method:
!  - Initializations (for the grib library)
!  - Opening the output grib file
!  - Scanning through the list (loop over all variables)
!  - Closing the output grib file
!
! Output files:
!  Output grib files for initial data (laf*) or boundary files (lbff*).
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  numlist,      & ! number of elements in ylist
  nnow            ! actual time step for which the data are computed

CHARACTER (LEN=10),       INTENT(IN)     ::    &
  ylist(numlist)  ! list of variables for output

CHARACTER (LEN=14),       INTENT(IN)     ::  &
  ydate1          ! actual date in the form   ddmmyyyyhhmmss

LOGICAL,                  INTENT(IN)     ::  &
  lwrite_hhl      ! to write the file with HHL levels

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)        ::      &
  izlist, iztable, k, kztop, kzbot, nzrecords, nzbyte, lfdlm, ldslm, izstat, &
  iztabllen, izlistlen, i, j, n, kzlevtot, lfdapi

INTEGER (KIND=iintegers)        ::      &
  numloclist,        & ! number of actual elements in ylist
  nuedat, izerror,   & !
  niostat,           & !
  izdebug,           & ! for verbosity of output
  ntrans_out=19        ! unit number for ready-files

INTEGER (KIND=iintegers)        ::      &
  iorg_data(3,0:num_compute-1)      ! for controlling NetCDF IO

CHARACTER (LEN=250)             ::      &
  yzname,      & ! name of ready files
  yzpath         ! full path name of files

CHARACTER (LEN= 80)             ::      &
  yzerrmsg

CHARACTER (LEN= 25)             ::      &
  yzroutine

CHARACTER (LEN=  3)             ::      &
  yzhead   

CHARACTER (LEN=10)              ::      &
  yzloclist(nvar_lm_chem)  ! local list of variables for output

CHARACTER (LEN= 11)             ::      &
  yzlistname, yztablname

CHARACTER (LEN=100)        ::  &
  charbuf

INTEGER (KIND=iintegers)        ::      &
  iloc_table(nvar_lm_chem),             &
  ivar_id   (nvar_lm_chem)

! Local arrays:
REAL (KIND=ireals)              ::      &
  z2darray  (ie2lm_max,je2lm_max)
     ! two dimensional array for sending to another node

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Section 1: Initializations
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

  yztablname = '           '
  yzlistname = '           '
  yzroutine  = 'org_lm_output'

  ! Allocate variables
  ! lfd and lds have to be calculated again for the output
  nzbyte = NINT (REAL(nrbit / 8, ireals))        ! usually 2 bytes
  ldslm   = ielm_tot * jelm_tot
  lfdlm  = ldslm  * nzbyte / iwlength + 2000
  IF (lwrite_hhl) THEN
    lfdapi = ldslm  *  8                + 5000
  ELSE
    lfdapi = ldslm  * nzbyte            + 5000
  ENDIF
  lds = INT (ldslm,  intgribf)
  lfd = INT (lfdlm,  intgribf)
  lfa = INT (lfdapi, intgribf)

  ALLOCATE (iblock(lfd), ibmap(nbitmap), dsup(ndsup),   STAT=izstat)
  ALLOCATE (ds_grib_single(lds), ds_grib_double(lds),   STAT=izstat)
  ALLOCATE (ymessage(lfa),                              STAT=izstat)

  IF (lwrite_hhl) THEN
    ALLOCATE(procarray_double(ie2lm_max,je2lm_max,num_compute), STAT=izstat)

    IF (my_cart_id == 0) THEN
      ! Create a new UUID for the HHL file identifier
      CALL uuid_create(uuid_out)

      ! Distribute this uuid to all PEs
      DO i = 1, 16
        charbuf(i:i) = uuid_out(i)
      ENDDO
    ENDIF

    IF (num_compute > 1) THEN
      CALL distribute_values (charbuf,  1, 0, imp_character, icomm_cart, izstat)
    ENDIF

    IF (my_cart_id /= 0) THEN
      DO i = 1, 16
        uuid_out(i) = charbuf(i:i)
      ENDDO
    ENDIF

    CALL uuid_2char (uuid_out, uuid_out_string)
    vcoord%vc_uuid(:) = uuid_out(:)
    IF (my_cart_id == 0) THEN
      PRINT *, '  Created new UUID for HHL-VGrid:  ', uuid_out_string
    ENDIF

    IF (izdebug > 20) THEN
      PRINT *, '  Start output for writing HHL: ', my_cart_id, ylist(1), numlist
    ENDIF
  ELSE
    ALLOCATE(procarray_single(ie2lm_max,je2lm_max,num_compute), STAT=izstat)
    IF (izdebug > 20) THEN
      PRINT *, '  Start output for writing results: ', my_cart_id, ylist(1), ' etc.', numlist
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Open the grib file
!------------------------------------------------------------------------------

  ! creating filename

  IF (lwrite_hhl) THEN
    yzpath = TRIM(ylm_cat)//TRIM(ylm_hhl)
  ELSE
    IF (lcomp_bound) THEN
      yzhead = 'lbf'
    ELSE
      yzhead = 'laf'
    ENDIF

    CALL make_fn (yzhead, ydate1, ydate_ini, ytunit_out, yextension, nnow, dt, .TRUE., &
                  itype_calendar, ylm_cat, yzpath, lmmss_ini, izdebug, izerror)

    IF (ylm_form_write == 'ncdf') THEN
      yzpath = TRIM(yzpath)//'.nc'
    ENDIF
  ENDIF

  CALL open_file (nuedat, yzpath, ymode_write, ylm_form_write, icomm_cart, &
                  my_cart_id, num_compute, lasync_io, idbg_level,          &
                  yzerrmsg, izerror)

  ! Write the headline in YUCHKDAT for this grib file
  IF (lchkout .AND. (my_cart_id == 0)) THEN
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

    WRITE (nuchkdat,'(A)') '       '
    WRITE (nuchkdat,'(A)') '       '
    WRITE (nuchkdat,'(A,A)') 'Check LM file:  ', yzpath(1:LEN_TRIM(yzpath))
    WRITE (nuchkdat,'(A)')   '   '
    WRITE (nuchkdat,'(A,I5,A,I5,A,I5)')                                    &
     '   ielm_tot =',ielm_tot,'  jelm_tot =',jelm_tot,'  kelm_tot =',kelm_tot
    WRITE (nuchkdat,'(A)')   '   '
    WRITE (nuchkdat,'(A,A)')                                               &
          '   var        ee  lev         min   ',                          &
          'imin jmin                max   imax jmax               mean '
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Look for output variables in LM variable table
!------------------------------------------------------------------------------

  yzloclist(1:numlist) = ylist(1:numlist)
  numloclist = numlist

  ! loop over all variables that should be written and loop over all variables
  ! in the LM variable table until equal elements are found
  ! to fill up iloc_table (which is needed 2 times below)
  DO izlist = 1, numloclist
    ! indices of field in variable table
    DO iztable = 1,nvar_lm
      ! searching for the output variable in the table
      IF (var_lm(iztable)%name == yzloclist(izlist)) THEN
         iloc_table(izlist) = iztable
      ENDIF
    ENDDO
  ENDDO

#ifdef NETCDF
  ! Write the headers in the NetCDF file
  IF (ylm_form_write == 'ncdf') THEN
    CALL write_nc_gdefs (nuedat, nnow, icomm_cart, num_compute, yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 2013, yzerrmsg, yzroutine)
    ENDIF

    CALL write_nc_vdefs (nuedat, nnow, numloclist, ivar_id, iloc_table, &
                         icomm_cart, num_compute, yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 2014, yzerrmsg, yzroutine)
    ENDIF
  ENDIF
#endif

  nzrecords = 0

  ! loop over all variables that should be written and loop over all variables
  ! in the LM variable table until equal elements are found
  DO izlist = 1, numloclist
    izlistlen  = LEN_TRIM(yzloclist(izlist))
    yzlistname = yzloclist(izlist)(1:izlistlen)

    iztable = iloc_table(izlist)

    SELECT CASE (var_lm(iztable)%rank)
      ! pack data depending on the rank

    CASE(3)
      kzbot = 1      ! in LM: LBOUND(var_lm(iztable)%p3,3)
      kztop = kelm   ! in LM: UBOUND(var_lm(iztable)%p3,3)
      kzlevtot = kztop
      IF     (TRIM(var_lm(iztable)%name) == 'W' .OR. TRIM(var_lm(iztable)%name) == 'HHL') THEN
        kztop = kztop + 1
        kzlevtot = kztop
      ELSEIF (TRIM(var_lm(iztable)%name) == 'W_SO') THEN
        kztop = ke_soil_lm+1
        kzlevtot = kztop
      ELSEIF (TRIM(var_lm(iztable)%name) == 'T_SO') THEN
        IF     ( (ylm_form_write == 'grb1') .OR. (ylm_form_write(1:3) == 'api') ) THEN
          kzbot = 0
          kzlevtot = kztop + 1
        ELSEIF (ylm_form_write == 'ncdf') THEN
          kzbot = 1
          kzlevtot = kztop
        ENDIF
        kztop = ke_soil_lm+1
      ELSEIF (TRIM(var_lm(iztable)%name) == 'HORIZON') THEN
        kztop = nhori
        kzlevtot = kztop
      ENDIF
      ! Distribute the multidimensional fields to the PEs
      DO k = kzbot, kztop
        nzrecords = nzrecords + 1
        z2darray(1:ie2lm,1:je2lm) = var_lm(iztable)%p3(1:ie2lm,1:je2lm,k)

        IF (izdebug > 20) THEN
          PRINT *, '  doing output for ', yzlistname, ';  level: ', k
        ENDIF

        CALL output_field (z2darray, nzrecords, iztable, k,             &
                           nuedat, nnow, kzlevtot, ivar_id(izlist),     &
                           iorg_data, izdebug, lwrite_hhl, .FALSE.)
      ENDDO

    CASE(2)
      nzrecords   = nzrecords+1
      z2darray(1:ie2lm,1:je2lm) = var_lm(iztable)%p2(1:ie2lm,1:je2lm)

      IF (izdebug > 20) THEN
        PRINT *, '  doing output for ', yzlistname, ';  level: ', 1
      ENDIF

      CALL output_field (z2darray, nzrecords, iztable, 1,               &
                         nuedat, nnow, 1, ivar_id(izlist),              &
                         iorg_data, izdebug, lwrite_hhl, .FALSE.)
    END SELECT
  ENDDO

!------------------------------------------------------------------------------
! Section 4: Flush output buffers and close grib file 
!------------------------------------------------------------------------------

  CALL output_field (z2darray, -1, -1, -1, nuedat, nnow,                &
                     -1, -1, iorg_data, izdebug, lwrite_hhl, .TRUE.)

  CALL close_file (nuedat, ylm_form_write, icomm_cart, my_cart_id,      &
                   num_compute, lasync_io, idbg_level, yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, izerror, yzerrmsg, 'close_file')
  ENDIF

  ! Write a blank line to YUCHKDAT
  IF ( (lchkout) .AND. (my_cart_id == 0) ) THEN
    WRITE (nuchkdat,'(A)') '   '
    WRITE (nuchkdat,'(A)') '   '
  ENDIF

  IF ( (lchkout) .AND. (my_cart_id == 0) ) THEN
    CLOSE (nuchkdat, STATUS='KEEP')
  ENDIF

  ! Deallocate the grib routines
  DEALLOCATE (iblock, ibmap, dsup, ds_grib_single, ds_grib_double, ymessage)

  IF (lwrite_hhl) THEN
    DEALLOCATE(procarray_double)
  ELSE
    DEALLOCATE(procarray_single)
  ENDIF

  IF (izdebug > 20) THEN
    PRINT *, '  output finished ', my_cart_id
  ENDIF

!------------------------------------------------------------------------------
! Section 5: Write a ready-file, if required
!------------------------------------------------------------------------------

  IF ((ytrans_out /= '   ') .AND. (my_cart_id == 0) .AND. (.NOT. lwrite_hhl)) THEN
    k = LEN_TRIM(ytrans_out)
    IF (.NOT. lcomp_bound) THEN
      yzhead = 'LMA'
    ELSE
      yzhead = 'LMB'
    ENDIF

    CALL make_fn (yzhead, ydate_ini, ydate_ini, 'f', yextension, nnow, dt, .TRUE., &
                  itype_calendar, ytrans_out, yzname, lmmss_ini, izdebug, izerror)

    ! Write the file
    OPEN  (ntrans_out, FILE=yzname, FORM='FORMATTED')
    WRITE (ntrans_out, '(A)') 'ready'
    CLOSE (ntrans_out)
  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE org_lm_output

!==============================================================================
!+ distributes records to PEs and packs them into grib format
!------------------------------------------------------------------------------

SUBROUTINE output_field (fldwb, irec, itable, ilev, nudat, nnow, klevels,   &
                         incdf_var_id, my_orgdata, izdbg, lwrite_hhl, lflush)

!------------------------------------------------------------------------------
!
! Description:
!  output_field distributes records to the PEs for packing these into 
!  grib format. First, the records are only gathered from all PEs and
!  every PE stores one record into the variable procarray. Only if every
!  PE has got a record, the data are packed and written to disk (in the
!  routine "write_grib"). If some PEs have got no record because no more
!  records are left, the output buffers (variable procarray) are "flushed".
!
! Method:
!  output_field is called for every record that is processed. The PE that
!  gets a special record, saves the characteristics of this record for the
!  output step later on.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  nudat,   & ! unit number of the grib file
  nnow,    & ! actual time step to be processed
  irec,    & ! number of record to be processed
  itable,  & ! location of the variable in the LM variable table
  ilev,    & ! level of a multi-level field
  klevels, & ! number of levels this variable has
  izdbg      ! for debug output

LOGICAL                 , INTENT (IN)    ::    &
  lwrite_hhl, & ! to indicate, whether hhl is written
  lflush        ! for flushing the output buffers

! Array arguments with intent(in):
REAL (KIND=ireals)      , INTENT (IN)    ::    &
  fldwb   (ie2lm_max,je2lm_max)   ! values of the variable to be processed

! Array arguments with intent(inout):
!REAL (KIND=irealgrib)   , INTENT (INOUT) ::    &
!  procarray(ie2lm_max,je2lm_max,num_compute)

INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  incdf_var_id              ! NetCDF-ID of each variable in the output list
                            ! only PE 0 has a reasonable value here

INTEGER (KIND=iintegers), INTENT(INOUT)  ::    &
  my_orgdata(3,0:num_compute-1) ! necessary only for PE 0 to save information

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=intgribf)              ::   &
  izerrf,        & ! error code for grib routines
  iz_ps=1,       & ! for passing to grbex1
  my_ilevf         ! ilev with KIND-parameter intgribf

INTEGER (KIND=iintegers)             ::   &
  izerror,       & ! error code
  irecord_len,   & ! length of grib record 
  nzpe,          & ! ID of PE that gets the record to be processed
  iz_lfd, iz_lfa,& ! length of iblock in KIND=iintegers
  nzactlocdef,   & ! to store actual localDefinitionNumber
  izgrib1id,     & ! grib handler for this output variable
  izgrib2id        ! grib handler for this output variable

INTEGER (KIND=int_ga)                ::   &
  irecord_lga      ! length of grib record

INTEGER (KIND=iintegers)             ::   &
  ilo, iup, jlo, jup, iloclo, ilocup, jloclo, jlocup, i, j, ij, k

REAL    (KIND=ireals)                ::   &
  zstartlon_tot, zstartlat_tot,   & ! for gribing lower left corner
  zendlon_tot,   zendlat_tot,     & ! for gribing upper right corner
  zbias, zfactor                    ! for scaling the record

CHARACTER (LEN=25)                   ::   &
  yzroutine

CHARACTER (LEN=80)                   ::   &
  yzerrmsg

CHARACTER(LEN= 8)            :: yzdate
CHARACTER(LEN=10)            :: yztime

INTEGER (KIND=iintegers)    :: icy, ijj, imm, idd, ihh, imi

INTEGER (KIND=iintegers), SAVE       ::   &
  my_itable, my_ilev, my_irec   ! save these variables until the output is done

LOGICAL                 , SAVE       ::   &
  loutput

! Local arrays
REAL (KIND=irealgrib)                ::   &
  undefsub                           ! for passing to subroutines

REAL (KIND=irealgrib)                ::   &
  fld_s_max     (ie2lm_max, je2lm_max),   & ! local 2D field with max. dimensions
  fld_tot_single(ie2lm_tot, je2lm_tot)      ! global 2D field

REAL (KIND=ireals)                   ::   &
  fld_d_max     (ie2lm_max, je2lm_max),   & ! local 2D field with max. dimensions
  fld_tot_double(ie2lm_tot, je2lm_tot)      ! global 2D field

CHARACTER (LEN=5)          :: ydbtype

#ifdef GRIBDWD
INTEGER (KIND=intgribf), EXTERNAL  :: IREFTS
#endif
#ifdef GRIBAPI
INTEGER (KIND=kindOfSize) :: ibyte_size_out
#endif

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Section 1: Initializations
!------------------------------------------------------------------------------

  yzroutine = 'output_field'
  iz_lfd   = INT (lfd, iintegers)
  iz_lfa   = INT (lfa, iintegers)
  IF (irec == 1) THEN
    my_orgdata(:,:) = 0_iintegers
    loutput         = .FALSE.
  ENDIF
  izerror  = 0_iintegers
  izerrf   = 0_intgribf
  ydbtype  = '   '

  IF (ylm_form_write == 'ncdf') THEN
    undefsub  = undefncdf
  ELSE
    undefsub  = undefgrib
  ENDIF
  undef = REAL (undefncdf, ireals)

!------------------------------------------------------------------------------
!  Section 2: If this is not a call to only flush the output data,
!             gather the field on the next free PE
!------------------------------------------------------------------------------

  ! Get the number of the PE to deal with that slice
  nzpe = MOD (irec-1,num_compute)

  IF (.NOT.lflush) THEN

    !--------------------------------------------------------------------------
    !  Section 2.1: Store some organizational data
    !--------------------------------------------------------------------------

    IF (my_cart_id == nzpe ) THEN
      my_irec   = irec
      my_itable = itable
      my_ilev   = ilev
      loutput   = .TRUE.
    ENDIF

    ! Processor 0 has to save some organizational data
    IF (my_cart_id == 0) THEN
      my_orgdata(1,nzpe) = ilev
      my_orgdata(2,nzpe) = klevels
      my_orgdata(3,nzpe) = incdf_var_id
    ENDIF

    !--------------------------------------------------------------------------
    !  Section 2.2: Scale and correct data where necessary
    !--------------------------------------------------------------------------

    IF (lwrite_hhl) THEN
      DO j = 1, je2lm
        DO i = 1, ie2lm
          ! this is only HHL (up to now) and needs no scaling
          fld_d_max(i,j) = fldwb(i,j)
        ENDDO
      ENDDO
    ELSE
      IF     ( (ylm_form_write == 'grb1') .OR. (ylm_form_write(1:3) == 'api') ) THEN
        zbias   = var_lm(itable)%bias
        zfactor = var_lm(itable)%factor
        IF (TRIM(var_lm(itable)%lsm) == 'l') THEN
          ! Check for undefined values, which could come from ncdf-input
          ! These values cannot be grib-ed
          DO j = 1, je2lm
            DO i = 1, ie2lm
              IF (fldwb(i,j) == undef) THEN
                fld_s_max(i,j) =  0.0_irealgrib
              ELSE
                fld_s_max(i,j) =  REAL ((fldwb(i,j) + zbias) * zfactor, irealgrib)
              ENDIF
            ENDDO
          ENDDO
        ELSE IF (TRIM(var_lm(itable)%lsm) == 'i') THEN
          ! Check for undefined values, which could come from ncdf-input
          ! These values cannot be grib-ed
  
          ! Comments by Dmitrii:
          ! At this point, all FLake variables and external-parameter fields
          ! should be set consistently:
          ! Then it is not necessary to reset those variables to zero.
          ! Furthermore, it should not be done since
          ! - setting lake depth to zero for non-lake-type grid boxes (sea and land)
          !   is dangerous and may result in errors during the COSMO-model run,
          ! - setting lake temperatures to zero leads to a loss of accuracy
          !   due to GRIB encoding.
  
          IF (TRIM(var_lm(itable)%name) == 'DEPTH_LK' ) THEN
            ! Only if NetCDF input is used, depth_lk_lm can have -1E20 as undef-value
            ! This is reset to -1 for Grib output
            DO j = 1, je2lm
              DO i = 1, ie2lm
                IF (fldwb(i,j) < -2.0_ireals) THEN
                  ! this can happen with NetCDF input
                  fld_s_max(i,j) =  -1.0_irealgrib
                ELSE
                  fld_s_max(i,j) =  REAL ((fldwb(i,j) + zbias) * zfactor, irealgrib)
                ENDIF
              ENDDO
            ENDDO
  
          ELSEIF (TRIM(var_lm(itable)%name) == 'SALT_LK' ) THEN
            ! To be on a safe side, set the salinity to zero for non-lake-type
            ! grid boxes (sea and land) or if salt_lk_lm proved to be undefined.
            DO j = 1, je2lm
              DO i = 1, ie2lm
                IF (fldwb(i,j) == undef .OR. depth_lk_lm(i,j) < 0.0_ireals) THEN
                  fld_s_max(i,j) =  0.0_irealgrib
                ELSE
                  fld_s_max(i,j) =  REAL ((fldwb(i,j) + zbias) * zfactor, irealgrib)
                ENDIF
              ENDDO
            ENDDO
          ELSE
            ! otherwise take the field as it is
            fld_s_max(1:ie2lm,1:je2lm) =                                          &
                   REAL ((fldwb(1:ie2lm,1:je2lm) + zbias) * zfactor, irealgrib)
          ENDIF
        ELSE
          ! there should be no undefined values
          fld_s_max(1:ie2lm,1:je2lm) =                                            &
                   REAL ((fldwb(1:ie2lm,1:je2lm) + zbias) * zfactor, irealgrib)
        ENDIF
      ELSEIF (ylm_form_write == 'ncdf') THEN
        IF (TRIM(var_lm(itable)%lsm) == 'l') THEN
          DO j = 1, je2lm
            DO i = 1, ie2lm
              IF (.NOT. lolp_lm(i,j)) THEN
                fld_s_max(i,j) = undefncdf
              ELSE
                fld_s_max(i,j) = REAL (fldwb(i,j), irealgrib)
              ENDIF
            ENDDO
          ENDDO
        ELSE IF (TRIM(var_lm(itable)%lsm) == 'i') THEN
          DO j = 1, je2lm
            DO i = 1, ie2lm
              IF (fr_lake_lm(i,j) <= 0.5) THEN
                fld_s_max(i,j) = undefncdf
              ELSE
                fld_s_max(i,j) = REAL (fldwb(i,j), irealgrib)
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO j = 1, je2lm
            DO i = 1, ie2lm
                fld_s_max(i,j) = REAL (fldwb(i,j), irealgrib)
            ENDDO
          ENDDO
        ENDIF
        IF  (TRIM(var_lm(itable)%name) == 'T_SNOW' ) THEN
          ! need DO loops instead of WHERE, because w_snow_lm and fld_s_max have
          ! different dimesions
          DO j = 1, je2lm
            DO i = 1, ie2lm
              IF (w_snow_lm(i,j) == 0.0_ireals) fld_s_max(i,j) = undefncdf
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    ENDIF

    !--------------------------------------------------------------------------
    !  Section 2.3: Gather the data
    !--------------------------------------------------------------------------

    IF (num_compute > 1) THEN
      IF (lwrite_hhl) THEN
        CALL gather_values (fld_d_max, procarray_double, ie2lm_max, je2lm_max, &
                num_compute, imp_reals, nzpe, icomm_cart, yzerrmsg, izerror)
      ELSE
        CALL gather_values (fld_s_max, procarray_single, ie2lm_max, je2lm_max, &
                num_compute, imp_grib, nzpe, icomm_cart, yzerrmsg, izerror)
      ENDIF
    ELSE
      IF (lwrite_hhl) THEN
        procarray_double(:,:,1) = fld_d_max(:,:)
      ELSE
        procarray_single(:,:,1) = fld_s_max(:,:)
      ENDIF
    ENDIF

  ENDIF

  IF ( (idbg_level > 10) .AND. (my_cart_id == nzpe) ) THEN
    PRINT *, '   PE ', my_cart_id, ' deals with ', var_lm(my_itable)%name,  &
                       my_irec, my_itable, my_ilev
  ENDIF

!------------------------------------------------------------------------------
!  Section 3: If lflush is .TRUE. or all PEs have gotten data, do the output
!------------------------------------------------------------------------------

  IF ( lflush .OR. MOD(irec-1,num_compute) == num_compute-1) THEN

    IF ( loutput ) THEN

      !------------------------------------------------------------------------
      !  Section 3.1: Complete GRIB meta data
      !------------------------------------------------------------------------

#if defined(GRIBDWD) || defined(GRIBAPI)
      IF ( (ylm_form_write == 'grb1') .OR. (ylm_form_write == 'api1') ) THEN

        ! Set vertical coordinate parameters only in case of levtyp 109/110
        IF ( (var_lm(my_itable)%levtyp == 109) .OR. (var_lm(my_itable)%levtyp == 110) ) THEN
          igds_out(1) = 42 + (inrvert_out) * 4

          ! number of the vertical coordinate parameters
          igds_out(2) = inrvert_out

          ! location of the list of vertical coordinate parameters in bytes
          igds_out(3) = 43

          ! vertical coordinate parameters
#ifdef GRIBDWD
          DO k = 1, inrvert_out
            igds_out(25 + k) = IREFTS( REAL(pv_out(k), irealgrib) )
          ENDDO
#endif
        ELSE
          igds_out(1) = 42

          ! number of the vertical coordinate parameters
          igds_out(2) = 0

          ! location of the list of vertical coordinate parameters in bytes
          igds_out(3) = -999999

          ! vertical coordinate parameters
          DO k = 1, inrvert_out
            igds_out(25 + k) = -999999
          ENDDO
        ENDIF

        ! complete gds: in case of grib1 these values will also be used for grib_api
        ! calculation of the left bottom corner
        IF ( (var_lm(my_itable)%tabtyp ==   2) .AND.             &
             (var_lm(my_itable)%levtyp == 110) .AND.             &
             (var_lm(my_itable)%ee     ==  33)       ) THEN
          ! u-gridpoint
          zstartlon_tot = startlon_tot + 0.5_ireals * dlon
          zendlon_tot   = endlon_tot   + 0.5_ireals * dlon
        ELSE
          zstartlon_tot = startlon_tot
          zendlon_tot   = endlon_tot
        ENDIF
        ! All longitude values have to be limited to the range (-180.0,+180.0)
        IF (zstartlon_tot > 180.0_ireals) THEN
          zstartlon_tot = zstartlon_tot - 360.0_ireals
        ENDIF
        IF (zendlon_tot > 180.0_ireals) THEN
          zendlon_tot = zendlon_tot - 360.0_ireals
        ENDIF
  
        IF ( (var_lm(my_itable)%tabtyp ==   2) .AND.             &
             (var_lm(my_itable)%levtyp == 110) .AND.             &
             (var_lm(my_itable)%ee     ==  34)       ) THEN
          ! v-gridpoint
          zstartlat_tot = startlat_tot + 0.5_ireals * dlat
          zendlat_tot   = endlat_tot   + 0.5_ireals * dlat
        ELSE
          zstartlat_tot = startlat_tot
          zendlat_tot   = endlat_tot
        ENDIF
        igds_out( 7) = NINT (zstartlat_tot * 1000.0_ireals)
        igds_out( 8) = NINT (zstartlon_tot * 1000.0_ireals)
        igds_out(10) = NINT (zendlat_tot   * 1000.0_ireals)
        igds_out(11) = NINT (zendlon_tot   * 1000.0_ireals)
  
        ! Reset idims_out(11): length of pds
        ! This cannot be done in init_lm_output which affects ALL grib output
        ! files leading to problems if ensemble boundary (lboundaries=.TRUE.
        ! and leps_bc=.TRUE.) is combined with linitial=.TRUE. but
        ! leps_ini=.FALSE. for example)
        idims_out(11) = 47_intgribf
        IF (leps_bc) THEN
          IF (lcomp_bound .AND. iepstyp_bc >= 0) idims_out(11) = 54_intgribf
        ENDIF

        ! create pds
        CALL makedwd_pds (my_itable, my_ilev, nnow)
      ENDIF
#endif

#if defined(GRIBAPI)
      IF (ylm_form_write == 'api1') THEN

        ! clone the sample to a special grib handler
        CALL grib_clone(igrib1_sample, izgrib1id)

        CALL grib_set (izgrib1id,'La1', igds_out( 7))           ! startlat
        CALL grib_set (izgrib1id,'Lo1', igds_out( 8))           ! startlon
        CALL grib_set (izgrib1id,'La2', igds_out(10))           ! endlat
        CALL grib_set (izgrib1id,'Lo2', igds_out(11))           ! endlon

        IF ( (var_lm(my_itable)%levtyp == 109) .OR. (var_lm(my_itable)%levtyp == 110) ) THEN
          CALL grib_set (izgrib1id,'PVPresent', 1)

          ! number of the vertical coordinate parameters
          CALL grib_set (izgrib1id,'NV',        inrvert_out)

          ! vertical coordinate parameters
          CALL grib_set (izgrib1id,'pv',        pv_out)
        ELSE
          CALL grib_set (izgrib1id,'PVPresent', 0)

          ! number of the vertical coordinate parameters
          CALL grib_set (izgrib1id,'NV',        0)
        ENDIF

        ! create pds
        CALL makeapi_pds (izgrib1id, my_itable, my_ilev, nnow)

      ELSEIF (ylm_form_write == 'api2') THEN

        ! clone the sample to a special grib handler
        CALL grib_clone(igrib2_sample, izgrib2id)

        ! complete local use section (CreationDate)
        ! -----------------------------------------

        SELECT CASE (nactlocdefnr) 
        CASE (252,253,254)
          ! we are in the DWD usage

          IF (leps_bc) THEN
            IF (lcomp_bound .AND. iepstyp_bc >= 0) THEN
              ! reset the local definition number
              CALL grib_set (izgrib2id, 'localDefinitionNumber',         253)  ! ensemble BC
              CALL grib_set (izgrib2id, 'localTypeOfEnsembleForecast', iepstyp_bc)
              nactlocdefnr = 253
            ENDIF
          ENDIF

          ! Set localInformationNumber, if necessary (ipds(41))
          IF (TRIM(var_lm(my_itable)%name) == 'PP' .OR.      &
              TRIM(var_lm(my_itable)%name) == 'P') THEN
            CALL grib_set (izgrib2id, 'localInformationNumber', refatm%irefatm_id)
          ELSEIF (TRIM(var_lm(my_itable)%name) == 'HHL') THEN
            CALL grib_set (izgrib2id, 'localInformationNumber', vcoord%ivcoord_id)
          ELSE
            ! usually not necessary for INT2LM products
            CALL grib_set (izgrib2id, 'localInformationNumber',       0)    !  ipds(41)
          ENDIF

          ! Set actual date and time (for 252, 253, 254)
          CALL DATE_AND_TIME(yzdate,yztime)
          READ(yzdate,'(I4,2I2)') icy, imm, idd
          READ(yztime,'(2I2,6X)') ihh, imi
          CALL grib_set (izgrib2id, 'localCreationDateYear',      icy)
          CALL grib_set (izgrib2id, 'localCreationDateMonth',     imm)
          CALL grib_set (izgrib2id, 'localCreationDateDay',       idd)
          CALL grib_set (izgrib2id, 'localCreationDateHour',      ihh)
          CALL grib_set (izgrib2id, 'localCreationDateMinute',    imi)
        END SELECT

        ! complete grid definition section
        ! --------------------------------

        ! calculation of the left bottom corner
        IF (var_lm(my_itable)%name == 'U         ') THEN
          ! u-gridpoint
          zstartlon_tot = startlon_tot + 0.5_ireals * dlon
          zendlon_tot   = endlon_tot   + 0.5_ireals * dlon
        ELSE
          zstartlon_tot = startlon_tot
          zendlon_tot   = endlon_tot
        ENDIF

        ! All longitude values have to be transformed to the range (0.0,+360.0)
        IF (zstartlon_tot < 0.0_ireals) THEN
          zstartlon_tot = zstartlon_tot + 360.0_ireals
        ENDIF
  
        IF (var_lm(my_itable)%name == 'V         ') THEN
          ! v-gridpoint
          zstartlat_tot = startlat_tot + 0.5_ireals * dlat
          zendlat_tot   = endlat_tot   + 0.5_ireals * dlat
        ELSE
          zstartlat_tot = startlat_tot
          zendlat_tot   = endlat_tot
        ENDIF

        CALL grib_set (izgrib2id, 'latitudeOfFirstGridPointInDegrees',  zstartlat_tot )
        CALL grib_set (izgrib2id, 'longitudeOfFirstGridPointInDegrees', zstartlon_tot )
        CALL grib_set (izgrib2id, 'latitudeOfLastGridPointInDegrees',     zendlat_tot )
        CALL grib_set (izgrib2id, 'longitudeOfLastGridPointInDegrees',    zendlon_tot )

        ! product definition template
        ! ---------------------------

        CALL makeapi2_pdt (izgrib2id, my_itable, my_ilev, nnow, lwrite_hhl)

      ENDIF
#endif

      !------------------------------------------------------------------------
      !  Section 3.2: combine the subarrays in the correct order
      !------------------------------------------------------------------------

      DO i = 1, num_compute
        ! Note: the extra boundary line is not considered here!!
        ! indices for field_tot
        ilo = isubpos(i-1,1)
        iup = isubpos(i-1,3)
        jlo = isubpos(i-1,2)
        jup = isubpos(i-1,4)

        ! indices for the part of procarray_single/double
        iloclo = 1 + nboundlines
        ilocup = isubpos(i-1,3) - isubpos(i-1,1) + 1 + nboundlines
        jloclo = 1 + nboundlines
        jlocup = isubpos(i-1,4) - isubpos(i-1,2) + 1 + nboundlines

        IF (lwrite_hhl) THEN
          fld_tot_double(ilo:iup,jlo:jup) =  procarray_double (iloclo:ilocup,jloclo:jlocup,i)
          ! need that for check_record
          fld_tot_single(ilo:iup,jlo:jup) =  REAL (procarray_double (iloclo:ilocup,jloclo:jlocup,i), irealgrib)
        ELSE
          fld_tot_single(ilo:iup,jlo:jup) = procarray_single (iloclo:ilocup,jloclo:jlocup,i)
        ENDIF
      ENDDO

      ! Put the packed field in array ds_grib_single/double
      IF (lwrite_hhl) THEN
        DO j = 1,jelm_tot
          DO i = 1,ielm_tot
            ij = (j-1)*ielm_tot + i
            ds_grib_double(ij) = fld_tot_double(i+1,j+1)
          ENDDO
        ENDDO
      ELSE
        IF (ylm_form_write(1:3) == 'api') THEN
          DO j = 1,jelm_tot
            DO i = 1,ielm_tot
              ij = (j-1)*ielm_tot + i
              ! grib_api needs double fields for writing
              ds_grib_double(ij) = REAL (fld_tot_single(i+1,j+1), ireals)
            ENDDO
          ENDDO
        ELSE
          DO j = 1,jelm_tot
            DO i = 1,ielm_tot
              ij = (j-1)*ielm_tot + i
              ds_grib_single(ij) = fld_tot_single(i+1,j+1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      !------------------------------------------------------------------------
      !  Section 3.3: Pack the field to GRIB (for NetCDF only set message length)
      !------------------------------------------------------------------------

      IF     (ylm_form_write == 'grb1') THEN
#ifdef GRIBDWD
        ! Set bitmap explicit to "no" again and reduce the size of ibmap
        ! (that seems to be a dummy only)
        ibms (3)     = -2

        ! Set these values again explicit, because they could be overwritten
        ! by input data
        ibds (2)     = 0
        ibds (5)     = nrbit

        ! degrib the level
        CALL grbex1(idwdednr, iz_ps, undefgrib, ndims, idims_out, ipds, igds_out, &
                    ibms, ibds, ibmap, dsup, ds_grib_single, iblock, izerrf)
        IF (izerrf /= 0) THEN
          yzerrmsg = 'error in grbex1'
          CALL model_abort (my_cart_id, 2022, yzerrmsg, yzroutine)
        ENDIF

        ! length of GRIB record in bytes
        irecord_len = idims_out(19)
#endif

#ifdef GRIBAPI
      ELSEIF (ylm_form_write == 'api1') THEN
        ! Set bitmap explicit to "no" again and reduce the size of ibmap 
        ! (that seems to be a dummy only)
        CALL grib_set (izgrib1id, 'bitmapPresent',      0)

        ! Set these values again explicit, because they could be overwritten
        ! by input data
        !  ?? ibds (2)     = 0
        CALL grib_set (izgrib1id, 'bitsPerValue',   nrbit)

        CALL grib_get_message_size(izgrib1id, ibyte_size_out, izerrf)
        CALL grib_set (izgrib1id, 'values',  ds_grib_double(:), izerrf)
        IF (izerrf /= GRIB_SUCCESS) THEN
          yzerrmsg = 'error in grib_set: values'
          CALL model_abort (my_cart_id, 2022, yzerrmsg, yzroutine)
        ENDIF

        ! length of GRIB record in bytes
        CALL grib_get_message_size(izgrib1id, ibyte_size_out, izerrf)

        IF (ibyte_size_out <= lfa) THEN
          CALL grib_get(izgrib1id, 'totalLength', irecord_lga)
          CALL grib_copy_message(izgrib1id,  ymessage)
        ELSE
          yzerrmsg = 'error with message length: ymessage too small: '
          CALL model_abort (my_cart_id, 2022, yzerrmsg, yzroutine)
        ENDIF

      ELSEIF (ylm_form_write == 'api2') THEN

        ! Set these values explicit
        IF (lwrite_hhl) THEN
          ! HHL is written with full precision
          ! set dataRepresentationTemplateNumber to 4
          CALL grib_set (izgrib2id, 'packingType',                 'grid_ieee')
!         CALL grib_set (izgrib2id, 'dataRepresentationTemplateNumber',      4)
          CALL grib_set (izgrib2id, 'precision',                             2)   ! 64-bit
        ELSE
          CALL grib_set (izgrib2id, 'packingType',               'grid_simple')
!         CALL grib_set (izgrib2id, 'dataRepresentationTemplateNumber',      0)
          CALL grib_set (izgrib2id, 'bitsPerValue',                      nrbit)
        ENDIF

        CALL grib_get_message_size(izgrib2id, ibyte_size_out, izerrf)
        CALL grib_set (izgrib2id, 'values', ds_grib_double(:), izerrf)
        IF (izerrf /= GRIB_SUCCESS) THEN
          yzerrmsg = 'error in grib_set: values'
          CALL model_abort (my_cart_id, 2022, yzerrmsg, yzroutine)
        ENDIF

        ! length of GRIB record in bytes
        CALL grib_get_message_size(izgrib2id, ibyte_size_out, izerrf)

        IF (ibyte_size_out <= lfa) THEN
          CALL grib_get(izgrib2id, 'totalLength', irecord_lga)
          CALL grib_copy_message(izgrib2id,  ymessage)
        ELSE
          yzerrmsg = 'error with message length: ymessage too small: '
          CALL model_abort (my_cart_id, 2022, yzerrmsg, yzroutine)
        ENDIF
#endif
      ELSEIF (ylm_form_write == 'ncdf') THEN
        irecord_len = ielm_tot * jelm_tot
      ENDIF
    ELSE

      irecord_len = 0
      irecord_lga = 0

    ENDIF

    !--------------------------------------------------------------------------
    !  Section 3.4: Check data and write output to disk
    !--------------------------------------------------------------------------

    ! check the data, if wanted  (cannot be called above, because at the end
    ! not all PEs have to do output)
    IF (lchkout) THEN
      my_ilevf  = INT (my_ilev, intgribf)
      CALL check_record (fld_tot_single, 1, ie2lm_tot  , 1, je2lm_tot  , 1, 1,      &
                            2, ie2lm_tot-1, 2, je2lm_tot-1, 1, 1, undefsub,  &
               var_lm(my_itable)%name, var_lm(my_itable)%ee, my_ilevf,       &
               loutput, nuchkdat, num_compute, icomm_cart, my_cart_id,       &
               yzerrmsg, izerror)
    ENDIF

    IF     (ylm_form_write == 'grb1') THEN

      CALL write_grib (nudat, iblock, irecord_len, iz_lfd, icomm_cart,       &
                       num_compute, lflush, ydbtype, lasync_io,              &
                       yzerrmsg, izerror)

    ELSEIF (ylm_form_write(1:3) == 'api') THEN

      CALL write_gribapi (nudat, ymessage, irecord_lga, iz_lfa, icomm_cart,  &
                          my_cart_id, num_compute, lflush, lasync_io,        &
                          yzerrmsg, izerror)

#ifdef GRIBAPI
      IF     (ylm_form_write == 'api1') THEN
        CALL grib_release (izgrib1id)
      ELSEIF (ylm_form_write == 'api2') THEN
        CALL grib_release (izgrib2id)
      ENDIF
#endif

    ELSEIF (ylm_form_write == 'ncdf') THEN

     CALL write_netcdf  (nudat, ds_grib_single, ielm_tot, jelm_tot,          &
                         irecord_len, my_orgdata, icomm_cart, my_cart_id,    &
                         num_compute, imp_grib, lasync_io, yzerrmsg, izerror)
    ENDIF

    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    loutput = .FALSE.

  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE output_field

!==============================================================================
!==============================================================================
!+ Module proc in src_lm_output to create product definition block with GRIBDWD
!------------------------------------------------------------------------------

SUBROUTINE makedwd_pds ( nloc, nlevel, ntstep)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure creates the product definition block according to the 
!   description of putpd1 (dwd) and the WMO.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers) , INTENT(IN) ::    &
  nloc,       & ! location of variable in LM variable table
  nlevel,     & ! level of multi-level fields
  ntstep        ! corresponds to timestep in a model running with dt

! Local scalars:
INTEGER  (KIND=iintegers)              :: ndgb, mmgb, jjgb, nhgb, jcgb,    &
                                          nmingb, nsecgb
REAL     (KIND=ireals)                 :: hm
CHARACTER(LEN= 8)                      :: ydate
CHARACTER(LEN=10)                      :: ytime
INTEGER  (KIND=iintegers)              :: jj,mm,dd,hh,mi, timestep

!- End of header
!------------------------------------------------------------------------------

  timestep = ntstep

  ! Initialisation
  ipds(:)  = -999999
   
  ! grib table nr. of version
  ipds(2)  = var_lm(nloc)%tabtyp
   
  ! centre identification
  ipds(3)  = ncenter
   
  ! modelidentification number
  IF (.NOT. lcomp_bound) THEN
    ipds(4)  =  nprocess_ini
  ELSE
    ipds(4)  =  nprocess_bd
  ENDIF
   
  ! number of used grid
  ipds(5)  = 255
   
  ! flag, that indicates whether GDS and BMS follow
  ipds(6)  = 128
   
  ! number of element
  ipds(7)  = var_lm(nloc)%ee
   
  ! leveltyp
  ipds(8)  = var_lm(nloc)%levtyp
   
  ! upper and lower boundary of the level
  SELECT CASE(var_lm(nloc)%levtyp)
  CASE (1,100,102,103,107,113,115,117,119,125,160,200,201)
    ipds(9)  = 0
    ipds(10) = 0
  CASE (109)
    ipds(9)  = 0
    ipds(10) = nlevel
  CASE (105,111)
    IF ( (var_lm(nloc)%name == 'W_SO      ') .OR.       &
         (var_lm(nloc)%name == 'T_SO      ') ) THEN
      ipds(9)  = 0.0_ireals        ! levtop
      ipds(10) = msoilgrib_lm(nlevel) ! depth of main soil level in cm
    ELSE
      ipds(9)  = 0
      ipds(10) = var_lm(nloc)%levbot
    ENDIF
  CASE (112)
    ipds(9)  = var_lm(nloc)%levtop
    ipds(10) = var_lm(nloc)%levbot
  CASE DEFAULT
    ipds(9)  = nlevel
    ipds(10) = nlevel + 1
  END SELECT 
   
  ! Reference time of data
  READ(ydate_ini,'(7I2)') jcgb, jjgb, mmgb, ndgb, nhgb, nmingb, nsecgb

  IF(jjgb == 0) THEN
    ipds(11) = 100
    ipds(22) = jcgb
  ELSE  
    ipds(11) = jjgb
    ipds(22) = jcgb + 1
  ENDIF
  ipds(12) = mmgb
  ipds(13) = ndgb
  ipds(14) = nhgb
  ipds(15) = nmingb
   
  ! Indicator of unit of time range
  ! (hm gives the seconds for this range
  ipds(16) = INT (nunit_of_time, intgribf)

  SELECT CASE (nunit_of_time)
  CASE ( 0)                  !  1 minute
    hm  =    60.0_ireals
  CASE ( 1)                  !  1 hour
    hm  =  3600.0_ireals
  CASE ( 2)                  !  1 day
    hm  = 86400.0_ireals
  CASE (10)                  !  3 hours
    hm  = 10800.0_ireals
  CASE (11)                  !  6 hours
    hm  = 21600.0_ireals
  CASE (12)                  ! 12 hours
    hm  = 43200.0_ireals
  CASE (13)                  ! 15 minutes
    hm  =   900.0_ireals
  CASE (14)                  ! 30 minutes
    hm  =  1800.0_ireals
  CASE DEFAULT
    PRINT *, ' ERROR *** Wrong value for unit-of-time:  ', nunit_of_time
    CALL model_abort (my_cart_id, 2015, 'wrong unit-of-time', 'makedwd_pds')
  END SELECT

  ! Time p1 and p2:  
  !  p1:       actual time for which a forecast product is valid
  !            last deletion of min, max-values (tri=2)
  !       or:  begin of averaging period        (tri=3)
  !  p2:       actual time for which a forecast product is valid
  !       or:  end of averaging period          (tri=3)
  ipds(17) = NINT (timestep*dt/hm)
  ipds(18) = 0
   
  ipds(19) = 0   ! Time Range indicator is always 0 in the INT2LM
  ipds(20) = 0
  ipds(21) = 0
  ipds(24) = 0
   
  ! local use area
  ipds(37)    = 254
  ipds(38:41) = 0
  
  IF (leps_bc) THEN
    IF (lcomp_bound .AND. iepstyp_bc >= 0) THEN
      ipds( 1)    = 66_intgribf
      ipds(37)    = 253_intgribf
      ipds(48:49) = 0_intgribf
      ipds(50)    = INT (iepstyp_bc, intgribf)
      ipds(51)    = INT (iepstot_bc, intgribf)
!US   ipds(51)    = 0_intgribf                  ! this must be a mistake???
      ipds(52)    = INT (iepsmem_bc, intgribf)
      ipds(53:54) = 0_intgribf
    ENDIF
  ENDIF

  ! actual date and time
  CALL DATE_AND_TIME(ydate,ytime)
  READ(ydate,'(I4,2I2)') jj,mm,dd
  READ(ytime,'(2I2,6X)') hh,mi
  ! According to DWD standard, ipds(42) contains the offset to the year 1900
  ipds(42)    = INT (jj-1900, intgribf)
  ipds(43)    = mm
  ipds(44)    = dd
  ipds(45)    = hh
  ipds(46)    = mi
  ipds(47)    = nvers

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE makedwd_pds

!==============================================================================
#ifdef GRIBAPI
!==============================================================================
!+ Module proc in src_lm_output to create product definition block with gribapi
!------------------------------------------------------------------------------

SUBROUTINE makeapi_pds (igrbhandl, nloc, nlevel, ntstep)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure creates the product definition block using the
!   grib_api library
!   But ipds has already been set by makedwd_pds, so we only have to 
!   set these values to the grib structure
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers) , INTENT(IN) ::    &
  igrbhandl,  & ! grib-api handle
  nloc,       & ! location of variable in LM variable table
  nlevel,     & ! level of multi-level fields
  ntstep        ! corresponds to timestep in a model running with dt

! Local scalars:
INTEGER  (KIND=iintegers)              :: ndate, ntime
REAL     (KIND=ireals)                 :: hm
CHARACTER                              :: ydate*8, ytime*10
INTEGER  (KIND=iintegers)              :: cc,jj,mm,dd,hh,mi, myy_api
INTEGER  (KIND=iintegers)              :: timestep

!
!- End of header
!------------------------------------------------------------------------------

  timestep = ntstep

  ! Initialisation
   
  CALL grib_set (igrbhandl,'table2Version',               ipds( 2)) ! var_lm(nloc)%tabtyp
  CALL grib_set (igrbhandl,'centre'       ,               ipds( 3)) ! ncenter
  CALL grib_set (igrbhandl,'generatingProcessIdentifier', ipds( 4)) ! nprocess_ini/_bd
   
  ! number of used grid: should be unnnecessary, we have it in the sample hopefully
  ! CALL grib_set (igrbhandl,'gridDefinition',            ipds( 5)) ! 255
  ! flag, that indicates whether GDS and BMS follow
  ! CALL grib_set (igrbhandl,'section1Flags',             ipds( 6)) ! 128
   
  CALL grib_set (igrbhandl,'indicatorOfParameter',        ipds( 7)) ! var_lm(nloc)%ee
  CALL grib_set (igrbhandl,'indicatorOfTypeOfLevel',      ipds( 8)) ! var_lm(nloc)%levtyp
  CALL grib_set (igrbhandl,'topLevel',                    ipds( 9)) ! var_lm(nloc)%levtop
  CALL grib_set (igrbhandl,'bottomLevel',                 ipds(10)) ! var_lm(nloc)%levbot
   
  ! Reference time of data
  READ(ydate_ini( 1: 8),'(I8)') ndate
  READ(ydate_ini( 9:12),'(I4)') ntime
  CALL grib_set (igrbhandl,'dataDate',                      ndate)  ! ipds(11-13)
  CALL grib_set (igrbhandl,'dataTime',                      ntime)  ! ipds(14-15)
   
  ! Indicator of unit of time range
  CALL grib_set (igrbhandl,'indicatorOfUnitOfTimeRange',  ipds(16)) ! INT (nunit_of_time, intgribf)
  CALL grib_set (igrbhandl,'startStep',                   ipds(17)) ! NINT(timestep*dt/hm)

  ! This must not be set for pds(19)=0 (stepType='instant'), because it would
  !  also reset the startStep to 0
  !!!  CALL grib_set (igrbhandl,'endStep',                     ipds(18)) ! 0

  CALL grib_set (igrbhandl,'timeRangeIndicator',          ipds(19)) ! 0 always in INT2LM

  SELECT CASE (ncenter)
  CASE (78)   ! DWD

    ! local use area
    CALL grib_set (igrbhandl,'localDefinitionNumber',       ipds(37)) ! ipds(37) = 254 / 253
  
    IF (leps_bc) THEN
      IF (lcomp_bound .AND. iepstyp_bc >= 0) THEN
        CALL grib_set (igrbhandl,'etyp',                    ipds(50)) ! ipds(50) = INT (iepstyp_bc, intgribf)
        CALL grib_set (igrbhandl,'etot',                    ipds(51)) ! ipds(51) = INT (iepstot_bc, intgribf)
        CALL grib_set (igrbhandl,'enum',                    ipds(52)) ! ipds(52) = INT (iepsmem_bc, intgribf)
      ENDIF
    ENDIF

    ! actual date and time
    CALL DATE_AND_TIME(ydate,ytime)
    READ(ydate,'(4I2)') cc,jj,mm,dd
    READ(ytime,'(2I2,6X)') hh,mi
    ! According to DWD standard, ipds(42) contains the offset to the year 1900
    myy_api = cc*100 + jj - 1900
    CALL grib_set (igrbhandl,'localDecodeDateYear',         myy_api)
    CALL grib_set (igrbhandl,'localDecodeDateMonth',             mm)
    CALL grib_set (igrbhandl,'localDecodeDateDay',               dd)
    CALL grib_set (igrbhandl,'localDecodeDateHour',              hh)
    CALL grib_set (igrbhandl,'localDecodeDateMinute',            mi)

    CALL grib_set (igrbhandl,'localVersionNumber',          ipds(47)) !  ipds(47)    = nvers

  END SELECT

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE makeapi_pds

!==============================================================================
!==============================================================================
!+ Module proc in src_lm_output to create product definition section for GRIB2
!------------------------------------------------------------------------------

SUBROUTINE makeapi2_pdt (igrbhandl, nloc, nlevel, ntstep, lwrite_hhl)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure sets values for the GRIB2 product definition section
!   using grib_api.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers) , INTENT(IN) ::    &
  igrbhandl,  & ! grib-api handle
  nloc,       & ! location of variable in LM variable table
  nlevel,     & ! level of multi-level fields
  ntstep        ! corresponds to timestep in a model running with dt

LOGICAL,                    INTENT(IN) ::    &
  lwrite_hhl    ! to write a new HHL field and generate a new UUID

! Local scalars:
INTEGER  (KIND=iintegers)    :: ndate, ntime
REAL     (KIND=ireals)       :: hm
CHARACTER(LEN=30)            :: local_levtyp
INTEGER  (KIND=iintegers)    :: timestep, numbervc, mlev, backprocid

!- End of header
!------------------------------------------------------------------------------

  timestep = ntstep

  ! set the shortname
  CALL grib_set (igrbhandl,'shortName', TRIM(var_lm(nloc)%name))

  ! number of the vertical coordinate parameters
  ! this is overwritten later on for hybrid and hybridLayer type
! CALL grib_set (igrbhandl,'NV',        0)
! CALL grib_set (igrbhandl,'PVPresent', 0)

  ! template number
  IF (leps_bc) THEN
    CALL grib_set (igrbhandl,'productDefinitionTemplateNumber', 1) ! individual ensemble forecast
  ELSE
    CALL grib_set (igrbhandl,'productDefinitionTemplateNumber', 0) ! Standard products
  ENDIF
   
  ! generating process
  IF (leps_bc) THEN
    CALL grib_set (igrbhandl,'typeOfGeneratingProcess',         4) ! Ensemble Forecast
  ELSE
    CALL grib_set (igrbhandl,'typeOfGeneratingProcess',         2) ! Forecast
!US CALL grib_set (igrbhandl,'typeOfGeneratingProcess',         0) ! so it is in conversion
  ENDIF

  ! background process: depends on local use of center
  SELECT CASE (ncenter)
  CASE (78)  ! DWD

    backprocid = INT (nvers/16384_iintegers)
    CALL grib_set (igrbhandl,'backgroundGeneratingProcessIdentifier', backprocid)

  CASE DEFAULT

    CALL grib_set (igrbhandl,'backgroundGeneratingProcessIdentifier',        255)

  END SELECT

  ! generating process identifier
  IF (.NOT. lcomp_bound) THEN
    CALL grib_set (igrbhandl,'generatingProcessIdentifier', nprocess_ini)
  ELSE
    CALL grib_set (igrbhandl,'generatingProcessIdentifier', nprocess_bd )
  ENDIF

  ! Unit of time and Forecast time
  ! From values 0-12, the indicatorOfUnitOfTimeRange has the same meaning in GRIB1 and GRIB2
  ! But the (DWD introduced) values 13 (15 minutes) and 14 (30 minutes) from GRIB1 have no
  ! correspondence in GRIB2 and have to be transferred to GRIB2 value 0 (minutes)
  IF ( (nunit_of_time == 13) .OR. (nunit_of_time == 14) ) THEN
    nunit_of_time = 0  ! minutes
  ENDIF
  CALL grib_set (igrbhandl,'indicatorOfUnitOfTimeRange',  nunit_of_time)
   

  ! Forecast Time 
  SELECT CASE (nunit_of_time)
  CASE ( 0)                  !  1 minute
    hm  =    60.0_ireals
  CASE ( 1)                  !  1 hour
    hm  =  3600.0_ireals
  CASE ( 2)                  !  1 day
    hm  = 86400.0_ireals
  CASE (10)                  !  3 hours
    hm  = 10800.0_ireals
  CASE (11)                  !  6 hours
    hm  = 21600.0_ireals
  CASE (12)                  ! 12 hours
    hm  = 43200.0_ireals
  CASE DEFAULT 
    PRINT *, ' ERROR *** Wrong value for unit-of-time:  ', nunit_of_time
    CALL model_abort (my_cart_id, 2015, 'wrong unit-of-time', 'makeapi2_pdt')
  END SELECT
  CALL grib_set (igrbhandl,'forecastTime',         NINT(timestep*dt/hm) )

  ! level type
  ! GRIB2: soil moisture variables should/could be defined for layers!
  IF (var_lm(nloc)%name(1:4) == 'W_SO') THEN
    local_levtyp = 'depthBelowLandLayer'
  ELSE
    local_levtyp = ylevltypes_out(var_lm(nloc)%levtyp, 2)
  ENDIF

  ! the leveltype can (has to ??) be set before all other level specifications
  ! only in case of generalVertical and generalVerticalLayer it must be set
  ! after setting genVertHeightCoords; otherwise NV cannot be set
  IF (local_levtyp(1:15) == 'generalVertical') THEN
    CALL grib_set (igrbhandl,'genVertHeightCoords',        1)
    CALL grib_set (igrbhandl, 'typeOfLevel',      TRIM(local_levtyp))
  ELSE
    CALL grib_set (igrbhandl, 'typeOfLevel',      TRIM(local_levtyp))
  ENDIF

  ! level top and bottom
  SELECT CASE(TRIM(local_levtyp))

  CASE('surface')

    CALL grib_set (igrbhandl, 'topLevel',      0)
    CALL grib_set (igrbhandl, 'bottomLevel',   0)

    CALL grib_set_missing (igrbhandl, 'scaleFactorOfFirstFixedSurface')
    CALL grib_set_missing (igrbhandl, 'scaledValueOfFirstFixedSurface')

    IF     (TRIM(var_lm(nloc)%name) == 'HSURF') THEN
      ! some additional things have to be set for HSURF
      CALL grib_set (igrbhandl,'typeOfSecondFixedSurface',          101)
      CALL grib_set_missing (igrbhandl,'scaleFactorOfSecondFixedSurface')
      CALL grib_set_missing (igrbhandl,'scaledValueOfSecondFixedSurface')
!   ELSEIF (TRIM(var_lm(nloc)%name) == 'QV_S') THEN
!     CALL grib_set (igrbhandl, 'scaleFactorOfFirstFixedSurface',     0)
!     CALL grib_set_missing (igrbhandl, 'scaleFactorOfSecondFixedSurface')
!     CALL grib_set_missing (igrbhandl, 'scaledValueOfSecondFixedSurface')
    ELSE
      CALL grib_set_missing (igrbhandl, 'scaleFactorOfSecondFixedSurface')
      CALL grib_set_missing (igrbhandl, 'scaledValueOfSecondFixedSurface')
    ENDIF

  CASE('hybrid')

    CALL grib_set (igrbhandl,'typeOfFirstFixedSurface',           105)
    CALL grib_set (igrbhandl,'scaledValueOfFirstFixedSurface', nlevel)

    CALL grib_set (igrbhandl,'typeOfSecondFixedSurface',          255)
!US CALL grib_set (igrbhandl,'typeOfSecondFixedSurface',            1)   ! such it is in converted data????
    CALL grib_set_missing (igrbhandl,'scaleFactorOfSecondFixedSurface')
    CALL grib_set_missing (igrbhandl,'scaledValueOfSecondFixedSurface')

    ! number of the vertical coordinate parameters
    CALL grib_set (igrbhandl,'PVPresent', 1)
    CALL grib_set (igrbhandl,'NV',        inrvert_out)

    ! vertical coordinate parameters
    CALL grib_set (igrbhandl,'pv',        pv_out)

  CASE('hybridLayer')

    CALL grib_set (igrbhandl, 'topLevel',                   nlevel)
    CALL grib_set (igrbhandl, 'bottomLevel',                nlevel + 1)

    ! number of the vertical coordinate parameters
    CALL grib_set (igrbhandl,'PVPresent', 1)
    CALL grib_set (igrbhandl,'NV',        inrvert_out)

    ! vertical coordinate parameters
    CALL grib_set (igrbhandl,'pv',        pv_out)

  CASE('generalVertical')

    CALL grib_set (igrbhandl, 'typeOfFirstFixedSurface',           150)
    CALL grib_set (igrbhandl, 'scaleFactorOfFirstFixedSurface',      0)
    CALL grib_set (igrbhandl, 'scaledValueOfFirstFixedSurface', nlevel)
!   CALL grib_set (igrbhandl, 'level',                          nlevel)

    CALL grib_set (igrbhandl,'NV',        6)
!   CALL grib_set (igrbhandl,'numberOfVerticalGridDescriptors',        6)

    ! Only for HHL
    IF (TRIM(var_lm(nloc)%name) == 'HHL') THEN
      CALL grib_set (igrbhandl, 'typeOfSecondFixedSurface',          101)
      CALL grib_set_missing (igrbhandl,'scaleFactorOfSecondFixedSurface')
      CALL grib_set_missing (igrbhandl,'scaledValueOfSecondFixedSurface')
    ENDIF

    CALL grib_set (igrbhandl, 'nlev',               vcoord%nlevels)
    CALL grib_set (igrbhandl, 'numberOfVGridUsed',  vcoord%ivctype)
    CALL grib_set (igrbhandl, 'uuidOfVGrid',        vcoord%vc_uuid)

  CASE('generalVerticalLayer')

    CALL grib_set (igrbhandl, 'typeOfFirstFixedSurface',           150)
    CALL grib_set (igrbhandl, 'topLevel',                       nlevel)

    CALL grib_set (igrbhandl,'typeOfSecondFixedSurface',           150)
    CALL grib_set (igrbhandl, 'bottomLevel',                nlevel + 1)

    CALL grib_set (igrbhandl,'NV',        6)
!   CALL grib_set (igrbhandl,'numberOfVerticalGridDescriptors',        6)

    CALL grib_set (igrbhandl,'nlev',                vcoord%nlevels)
    CALL grib_set (igrbhandl,'numberOfVGridUsed',   vcoord%ivctype)
    CALL grib_set (igrbhandl,'uuidOfVGrid',         vcoord%vc_uuid)

  CASE('depthBelowLand')

    IF     (nlevel == 0 .AND. var_lm(nloc)%levbot == 0) THEN      !for T_S, T_SO(0)

      CALL grib_set (igrbhandl, 'level',                0)

    ELSEIF (nlevel == 0 .AND. var_lm(nloc)%levbot /= 0) THEN ! for T_M, T_CL

      CALL grib_set (igrbhandl, 'scaleFactorOfFirstFixedSurface',          2)
      CALL grib_set (igrbhandl, 'scaledValueOfFirstFixedSurface',   var_lm(nloc)%levbot)
      CALL grib_set_missing (igrbhandl, 'scaleFactorOfSecondFixedSurface')
      CALL grib_set_missing (igrbhandl, 'scaledValueOfSecondFixedSurface')

    ELSE                                         !for T_SO (1...8)

      IF(nlevel == 1) THEN
        mlev = NINT (czmls_lm(nlevel) * 1000.0_ireals, iintegers)
        CALL grib_set (igrbhandl, 'scaleFactorOfFirstFixedSurface',     3)  ! 0.005 m
      ELSE
        mlev = NINT (czmls_lm(nlevel) *  100.0_ireals, iintegers)
        CALL grib_set (igrbhandl, 'scaleFactorOfFirstFixedSurface',     2)  ! values in cm
      ENDIF
      CALL grib_set (igrbhandl, 'scaledValueOfFirstFixedSurface',    mlev)
      CALL grib_set_missing (igrbhandl, 'scaleFactorOfSecondFixedSurface')
      CALL grib_set_missing (igrbhandl, 'scaledValueOfSecondFixedSurface')

    ENDIF

  CASE('depthBelowLandLayer')

    IF ( (var_lm(nloc)%levbot == 0) .AND. (var_lm(nloc)%levtop == 0) .AND. (nlevel /= 0)) THEN ! W_SO/W_SO_ICE

      mlev = NINT (czhls_lm(nlevel-1) * 100.0_ireals, iintegers)

      CALL grib_set (igrbhandl, 'scaleFactorOfFirstFixedSurface',       2)
      CALL grib_set (igrbhandl, 'scaledValueOfFirstFixedSurface',    mlev)

      mlev = NINT (czhls_lm(nlevel  ) * 100.0_ireals, iintegers)
      CALL grib_set (igrbhandl, 'scaleFactorOfSecondFixedSurface',      2)
      CALL grib_set (igrbhandl, 'scaledValueOfSecondFixedSurface',   mlev)

    ENDIF

  CASE DEFAULT

    CALL grib_set (igrbhandl, 'topLevel',      var_lm(nloc)%levtop)
    CALL grib_set (igrbhandl, 'bottomLevel',   var_lm(nloc)%levbot)

  END SELECT

  ! if this is set before the CASE-construct above, we cannot set the
  ! NV for generalVertical ?????????
  ! CALL grib_set (igrbhandl, 'typeOfLevel',      TRIM(local_levtyp))

  ! Ensemble settings
  IF (leps_bc) THEN
    CALL grib_set (igrbhandl, 'typeOfEnsembleForecast',             192)
    CALL grib_set (igrbhandl, 'perturbationNumber',          iepsmem_bc)
    CALL grib_set (igrbhandl, 'numberOfForecastsInEnsemble', iepstot_bc)
  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE makeapi2_pdt

!==============================================================================
#endif
!==============================================================================
!+ Module procedure in "io_utilities" for writing a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE write_nc_gdefs (ncid, nstep, icomm, npes, yerrmsg, ierror)
!
!------------------------------------------------------------------------------
!
! Description:
!   This routine initializes global definitions for a NetCDF output.
!   It writes the latitude and longitudes values and the vertical coordinates.
!
!------------------------------------------------------------------------------
!
!
! Subroutine / Function arguments

! Scalar arguments with intent(in):
  INTEGER (kind=iintegers),   INTENT(IN) :: &
    nstep,        & ! actual time step
    ncid            ! NetCDF file IDr

  INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
    icomm,    & ! MPI communicator
    npes        ! number of PEs

  CHARACTER (LEN= *),        INTENT(OUT)   ::  &
    yerrmsg     ! string for error messages

! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    ierror          ! error index

!-----------------------------------------------------------------------
!
! Local scalars:

INTEGER (KIND=iintegers), PARAMETER  :: nbds=2, ntime=1

INTEGER (KIND=iintegers)  :: i, j, timestep, isect

INTEGER (kind=iintegers):: &
    jgridVarID,     & ! NetCDF ID for rotated_pole
    jlonVarID,      & ! NetCDF ID for longitude
    jlatVarID,      & ! NetCDF ID for latitude
    jrlonVarID,     & ! NetCDF ID for rotated longitude
    jrlatVarID,     & ! NetCDF ID for rotated latitude
    jslonuVarID,    & ! NetCDF ID for U-direction shifted longitude
    jslatuVarID,    & ! NetCDF ID for U-direction shifted latititude
    jslonvVarID,    & ! NetCDF ID for V-direction shifted longitude
    jslatvVarID,    & ! NetCDF ID for V-direction shifted latitude
    jsrlonVarID,    & ! NetCDF ID for shifted rotated longitude
    jsrlatVarID,    & ! NetCDF ID for shifted rotated latitude
    jvcVarID,       & ! NetCDF ID for the vertical component
    jsoilVarID,     & ! NetCDF ID for the multi soil layer component
    jsoilbdsID,     & ! NetCDF ID for the multi soil layer bounds
    jsectVarID,     & ! NetCDF ID for the sectors of the variable HORIZON
    jtimeID,        & ! NetCDF ID for the time
    jtbdsID           ! NetCDF ID for the time bounds

  CHARACTER (LEN= 40) :: ydate

  CHARACTER (LEN=1)  :: &
    cvctype           ! character version of ivctype

  CHARACTER (LEN=19) :: &
    creation_date     ! actual time when the data file is written

  CHARACTER (LEN=80)  :: &
    institution       ! name of institution where the data is produced

  CHARACTER (LEN=20)  :: &
    source            ! model name and version number

!  Local arrays:

  REAL (kind=ireals)   :: &
!   PIK U. Boehm - 20.12.06
!   longitude(ielm_tot,jelm_tot),  & ! geographic longitudes
!   latitude(ielm_tot,jelm_tot),   & ! geographic latitudes
    zlongitude(ielm_tot,jelm_tot), & ! geographic longitudes
    zlatitude(ielm_tot,jelm_tot),  & ! geographic latitudes
!   PIK U. Boehm - End
    rotlon(ielm_tot),              & ! rotated longitudes
    rotlat(jelm_tot),              & ! rotated latitudes
    zslonu(ielm_tot,jelm_tot),     & ! U-direction shifted longitudes
    zslatu(ielm_tot,jelm_tot),     & ! U-direction shifted latitudes
    zslonv(ielm_tot,jelm_tot),     & ! V-direction shifted longitudes
    zslatv(ielm_tot,jelm_tot),     & ! V-direction shifted latitudes
    srlon(ielm_tot),               & ! shifted rotated longitudes
    srlat(jelm_tot)                  ! shifted rotated latitudes

  REAL (kind=ireals)  :: &
    time(ntime), time_bnds(nbds, ntime)

  REAL (kind=ireals)  :: &
    zsoil_bnds(nbds, ke_soil_lm+1)

  REAL (kind=ireals), ALLOCATABLE   :: &
    zvcoord(:)        ! vertical coordinate

  INTEGER (kind=iintegers) :: &
    zsect(nhori),             &
    ztime_values(8)   ! holds the date and time information taken from
                      ! internal subroutine date_and_time
   INTEGER (kind=iintegers) :: my_comm_id, implcode

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

#ifdef NETCDF

IF (npes > 1) THEN
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 1
    RETURN
  ENDIF
ENDIF
IF (npes > 1) THEN
 ! Get id in communicator comm
  CALL MPI_COMM_RANK(icomm,my_comm_id,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_COMM_RANK failed'
    ierror  = 2
    RETURN
  ENDIF
ELSE
  my_comm_id = 0
ENDIF

! processor 0 does the job
IF (my_comm_id == 0) THEN

! some pre-settings
  institution   = 'GKSS'
  source = 'int2lm2 beta 1'
  ierror = 0
  timestep = nstep

  time = 0.
  time_bnds = 0.


! determine the creation date
  CALL date_and_time(values=ztime_values)
  WRITE (creation_date,'(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') &
            ztime_values(1),'-',ztime_values(2),'-',ztime_values(3),' ', &
            ztime_values(5),':',ztime_values(6),':',ztime_values(7)

! write the global attributes
  IF (TRIM(yncglob_title) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "title",    TRIM(yncglob_title))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (TRIM(yncglob_institution) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "institution",               &
                                            TRIM(yncglob_institution))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (TRIM(yncglob_source) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "source",                    &
                                            TRIM(yncglob_source))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (TRIM(yncglob_project_id) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "project_id",                &
                                            TRIM(yncglob_project_id))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (TRIM(yncglob_experiment_id) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "experiment_id",             &
                                            TRIM(yncglob_experiment_id))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (ncglob_realization /= -999) THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "realization",ncglob_realization)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "Conventions","CF-1.5")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "conventionsURL", &
                         "http://www.cfconventions.org/")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  IF (TRIM(yncglob_contact) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "contact",                  &
                                           TRIM(yncglob_contact))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (TRIM(yncglob_references) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL,"references",                &
                                           TRIM(yncglob_references))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF

  ierror=nf90_put_att(ncid, NF90_GLOBAL, "creation_date", creation_date)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

! set the values of rotated North Pole in the grid_mapping attribute
  grid_mapping = 'rotated_pole'

  ierror=nf90_def_var(ncid, "rotated_pole", NF90_CHAR, jgridVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jgridVarID, "grid_mapping_name", "rotated_latitude_longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jgridVarID, "grid_north_pole_latitude", REAL(pollat))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jgridVarID, "grid_north_pole_longitude", REAL(pollon))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF (polgam /= 0._ireals) THEN
    ierror=nf90_put_att(ncid, jgridVarID, "north_pole_grid_longitude", REAL(polgam))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
   ENDIF

! set the values for rotated and true geographic longitudes and latitudes

  DO i = 1, ielm_tot
    rotlon(i) = startlon_tot + (i-1)*dlon
  ENDDO
  DO j = 1, jelm_tot
    rotlat(j) = startlat_tot + (j-1)*dlat
  ENDDO
  srlon = rotlon + dlon*0.5
  srlat = rotlat + dlat*0.5
  DO j = 1, jelm_tot
    DO i = 1, ielm_tot
      zlongitude(i,j)  = rlarot2rla(rotlat(j), rotlon(i), pollat, pollon, polgam)
      zlatitude (i,j)  = phirot2phi(rotlat(j), rotlon(i), pollat, pollon, polgam)
      zslonu(i,j)      = rlarot2rla (rotlat(j), srlon(i), pollat, pollon, polgam)
      zslatu(i,j)      = phirot2phi (rotlat(j), srlon(i), pollat, pollon, polgam)
      zslonv(i,j)      = rlarot2rla (srlat(j), rotlon(i), pollat, pollon, polgam)
      zslatv(i,j)      = phirot2phi (srlat(j), rotlon(i), pollat, pollon, polgam)
    ENDDO
  ENDDO

! define the attributes of longitude and latitude
  ierror=nf90_def_dim(ncid,"rlon",ielm_tot, idims_id(1))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_dim(ncid,"rlat",jelm_tot, idims_id(2))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_def_dim(ncid,"srlon",ielm_tot, idims_id(9))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_dim(ncid,"srlat",jelm_tot, idims_id(10))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! rotated longitude
  ierror=nf90_def_var(ncid, "rlon", NF90_FLOAT, (/ idims_id(1) /), jrlonVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jrlonVarID, "axis", "X")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlonVarID, "standard_name", "grid_longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlonVarID, "long_name", "rotated longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlonVarID, "units", "degrees")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! rotated latitude
  ierror=nf90_def_var(ncid, "rlat", NF90_FLOAT, (/ idims_id(2) /), jrlatVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlatVarID, "axis", "Y")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlatVarID, "standard_name", "grid_latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlatVarID, "long_name", "rotated latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlatVarID, "units", "degrees")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! geographic longitude
  ierror=nf90_def_var(ncid, "lon", NF90_FLOAT, (/ idims_id(1), idims_id(2) /), jlonVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jlonVarID, "standard_name", "longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlonVarID, "long_name", "longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlonVarID, "units", "degrees_east")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! geographic latitude
  ierror=nf90_def_var(ncid, "lat", NF90_FLOAT, (/ idims_id(1), idims_id(2) /), jlatVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlatVarID, "standard_name", "latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlatVarID, "long_name", "latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlatVarID, "units", "degrees_north")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! shifted rotated longitude
  ierror=nf90_def_var(ncid, "srlon", NF90_FLOAT, (/ idims_id(9) /), jsrlonVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jsrlonVarID, "axis", "X")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlonVarID, "standard_name", "grid_longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlonVarID, "long_name", "rotated longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlonVarID, "units", "degrees")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! shifted rotated latitude
  ierror=nf90_def_var(ncid, "srlat", NF90_FLOAT, (/ idims_id(10) /), jsrlatVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlatVarID, "axis", "Y")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlatVarID, "standard_name", "grid_latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlatVarID, "long_name", "rotated latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlatVarID, "units", "degrees")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

 ! shifted geographic longitude or U wind component
  ierror=nf90_def_var(ncid, "slonu", NF90_FLOAT, (/ idims_id(9), idims_id(2) /), jslonuVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslonuVarID, "standard_name", "longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslonuVarID, "long_name", "staggered U-wind longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslonuVarID, "units", "degrees_east")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! shifted geographic latitude in U wind component
  ierror=nf90_def_var(ncid, "slatu", NF90_FLOAT, (/ idims_id(9), idims_id(2) /), jslatuVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatuVarID, "standard_name", "latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatuVarID, "long_name", "staggered U-wind latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatuVarID, "units", "degrees_north")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! shifted geographic longitude or V wind component
  ierror=nf90_def_var(ncid, "slonv", NF90_FLOAT, (/ idims_id(1), idims_id(10) /), jslonvVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslonvVarID, "standard_name", "longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslonvVarID, "long_name", "staggered V-wind longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslonvVarID, "units", "degrees_east")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! shifted geographic latitude in V wind component
  ierror=nf90_def_var(ncid, "slatv", NF90_FLOAT, (/ idims_id(1), idims_id(10) /), jslatvVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatvVarID, "standard_name", "latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatvVarID, "long_name", "staggered V-wind latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatvVarID, "units", "degrees_north")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

!  IF (vcoord%ivctype == 1) THEN    ! hybrid sigma-pressure co-ordinate  !_br 31.03.09

    ALLOCATE (zvcoord(ke1lm))
    zvcoord(1:ke1lm) = vcoord%vert_coord(1:ke1lm)

    ierror=nf90_def_dim(ncid,"level",kelm, idims_id(3))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_def_dim(ncid,"level1",ke1lm, idims_id(4))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_def_var(ncid, "vcoord", NF90_FLOAT, (/ idims_id(4) /), jvcVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    IF (vcoord%ivctype == 1) THEN
      ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                 "Reference atmosphere: Pressure based hybrid coordinate")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
!      ierror=nf90_put_att(ncid, jvcVarID, "units", "Pa")
!      IF (ierror /= NF90_NOERR) THEN
!        yerrmsg = TRIM(NF90_strerror(ierror))
!        RETURN
!      ENDIF
    ELSE IF (vcoord%ivctype == 2) THEN
      ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                "Reference atmosphere: Height-based hybrid Gal-Chen coordinate")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
!      ierror=nf90_put_att(ncid, jvcVarID, "units", "m")
!      IF (ierror /= NF90_NOERR) THEN
!        yerrmsg = TRIM(NF90_strerror(ierror))
!        RETURN
!      ENDIF
    ELSE IF (vcoord%ivctype == 3) THEN
      ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                "Reference atmosphere: SLEVE coordinate")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ELSE IF (vcoord%ivctype == 4) THEN
      ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                "Reference atmosphere: SLEVE2 coordinate")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

    ierror=nf90_put_att(ncid, jvcVarID, "ivctype", vcoord%ivctype)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "irefatm", refatm%irefatm)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "p0sl", refatm%p0sl)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "t0sl", refatm%t0sl)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "dt0lp", refatm%dt0lp)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "vcflat", vcoord%vcflat)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    IF (vcoord%ivctype == 3 .OR. vcoord%ivctype ==4) THEN
      ierror=nf90_put_att(ncid, jvcVarID, "svc1", svc1)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "svc2", svc2)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "nfltvc", nfltvc)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF
    IF (refatm%irefatm == 2) THEN
      ierror=nf90_put_att(ncid, jvcVarID, "delta_t", refatm%delta_t)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "h_scal", refatm%h_scal)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

!    ierror=nf90_put_att(ncid, jvcVarID, "formula-terms", "t.b.d.")
!    IF (ierror /= NF90_NOERR) THEN
!      yerrmsg = TRIM(NF90_strerror(ierror))
!      RETURN
!    ENDIF
!  ELSE

!    ierror = 1
!    WRITE(cvctype,'(I1)') ivctype
!    yerrmsg = 'ivctype='//cvctype//' not implemented for NetCDF output'
!    RETURN

!  ENDIF


! define fields of variable multi layer soil model

! ierror=nf90_def_dim(ncid, "soil", ke_soil_lm, idims_id(7))
! IF (ierror /= NF90_NOERR) THEN
!   yerrmsg = TRIM(NF90_strerror(ierror))
!   RETURN
! ENDIF


  ierror=nf90_def_dim(ncid, "soil1", ke_soil_lm+1, idims_id(8))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_var(ncid, "soil1", NF90_FLOAT, (/ idims_id(8) /), jsoilVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jsoilVarID, "standard_name", "depth")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsoilVarID, "long_name", "depth of soil layers")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsoilVarID, "units", "m")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsoilVarID, "positive", "down")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_dim(ncid,"bnds", nbds, idims_id(6))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_def_var(ncid, "soil1_bnds", NF90_FLOAT, (/ idims_id(6), idims_id(8) /), jsoilbdsID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsoilbdsID, "long_name", "depth of soil layer bounds")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsoilbdsID, "units", "m")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsoilbdsID, "positive", "down")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF


! define dimensions of time and time_bounds
  ierror=nf90_def_dim(ncid,"time", NF90_UNLIMITED, idims_id(5))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_def_var(ncid, "time", NF90_DOUBLE, (/ idims_id(5) /), jtimeID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_var(ncid, "time_bnds", NF90_DOUBLE, (/ idims_id(6), idims_id(5) /), jtbdsID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF(lradtopo) THEN
     ! Modifications to write multi-dimensional HORIZON field
     ierror=nf90_def_dim(ncid, "nhori", nhori, idims_id(11)) !_br 20.09.2012

     IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
     ENDIF

     ierror=nf90_def_var(ncid, "nhori", NF90_FLOAT, (/ idims_id(11) /), jsectVarID) !_br 20.09.2012
     IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
     ENDIF

!_br 07.06.12 "sectors" is no valid CF standard_name
!     ierror=nf90_put_att(ncid, jsectVarID, "standard_name", "sectors")
!     IF (ierror /= NF90_NOERR) THEN
!        yerrmsg = TRIM(NF90_strerror(ierror))
!        RETURN
!     ENDIF
!_br 07.06.12 end
     ierror=nf90_put_att(ncid, jsectVarID, "long_name", "sectors of horizontal angles around grid cells")
     IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
     ENDIF
!_br 07.06.12 "sector number starting from the North clockwise" is no valid CF unit
!     ierror=nf90_put_att(ncid, jsectVarID, "units", "sector number starting from the North clockwise")
!     IF (ierror /= NF90_NOERR) THEN
!        yerrmsg = TRIM(NF90_strerror(ierror))
!        RETURN
!     ENDIF
  ENDIF


! determine the actual forecast time

  ierror=nf90_put_att(ncid, jtimeID, "long_name", "time")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
!_br 07.09.12
  IF (lmmss_ini) THEN
    ydate = 'seconds since '// ydate_ini(1:4)//'-'//ydate_ini(5:6)//'-'//ydate_ini(7:8)//' '&
            //ydate_ini(9:10)//':00:00'
  ELSE
    ydate = 'seconds since '// ydate_ini(1:4)//'-'//ydate_ini(5:6)//'-'//ydate_ini(7:8)//' '&
           //ydate_ini(9:10)//':'//ydate_ini(11:12)//':'//ydate_ini(13:14)
  ENDIF
!_br 07.09.12 end
  ierror = nf90_put_att (ncid, jtimeID, "units", ydate)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF     (itype_calendar == 0) THEN
    ierror = nf90_put_att (ncid, jtimeID, "calendar", "standard")
  ELSEIF (itype_calendar == 1) THEN
    ierror = nf90_put_att (ncid, jtimeID, "calendar", "360_day")
  ELSEIF (itype_calendar == 2) THEN
    ierror = nf90_put_att (ncid, jtimeID, "calendar", "365_day")
  ENDIF
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror = nf90_put_att (ncid, jtimeID, "bounds", "time_bnds")

! the time will be in seconds

  time(1)        = timestep*dt
  time_bnds(1,1) = 0.
  time_bnds(2,1) = timestep*dt


! End of definition mode
!!! no more attribute definitions beyond this line !!!

  ierror=nf90_enddef(ncid)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

! write the values of longitude, latitude and vertical axis to file
  ierror=nf90_put_var(ncid, jrlonVarID, rotlon)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jrlatVarID, rotlat)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
! PIK U. Boehm - 20.12.06
! ierror=nf90_put_var(ncid, jlonVarID,  longitude)
  ierror=nf90_put_var(ncid, jlonVarID, zlongitude)
! PIK U. Boehm - End
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
! PIK U. Boehm - 20.12.06
! ierror=nf90_put_var(ncid, jlatVarID,  latitude)
  ierror=nf90_put_var(ncid, jlatVarID, zlatitude)
! PIK U. Boehm - End
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jsrlonVarID, srlon)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jsrlatVarID, srlat)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jslonuVarID, zslonu)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jslatuVarID, zslatu)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jslonvVarID, zslonv)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jslatvVarID, zslatv)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_var(ncid, jvcVarID, zvcoord)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF (lmulti_layer_lm) THEN
    ierror=nf90_put_var(ncid, jsoilVarID, czmls_lm(1:ke_soil_lm+1))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    zsoil_bnds(1, 1)       = 0._ireals
    DO i = 1, ke_soil_lm
      zsoil_bnds(2, i)    = 2._ireals*czmls_lm(i) - zsoil_bnds(1, i)
      zsoil_bnds(1, i+1)  = zsoil_bnds(2, i)
    ENDDO
    zsoil_bnds(2, ke_soil_lm+1) = 2._ireals*czmls_lm(ke_soil_lm+1) - zsoil_bnds(1, ke_soil_lm+1)
    ierror=nf90_put_var(ncid, jsoilbdsID, zsoil_bnds)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF

  IF (lradtopo) THEN
     DO isect=1,nhori
        zsect(isect)=isect
     ENDDO
     ierror=nf90_put_var(ncid, jsectVarID, zsect)
     IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
     ENDIF
  ENDIF

! write the values of time and time_bnds to file
  ierror=nf90_put_var(ncid, jtimeID, time)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_var(ncid, jtbdsID, time_bnds)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  DEALLOCATE (zvcoord)

ENDIF   ! end of IF statement for processor 0

IF (npes > 1) THEN
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 1
    RETURN
  ENDIF
ENDIF

#endif

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE write_nc_gdefs

!==============================================================================
!+ Module procedure in "io_utilities" for writing a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE write_nc_vdefs (ncid, nstep, numlist, var_id, iloc, &
                          icomm, npes, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   write NetCDF attributes for each output parameter
!
!------------------------------------------------------------------------------

! Scalar arguments with intent(in):
  INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
    nstep,         & ! actual time step
    ncid             ! NetCDF file ID

  INTEGER (KIND=iintegers), INTENT(IN)     ::  &
    numlist          ! number of variables for output

  INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
    icomm,    & ! MPI communicator
    npes        ! number of PEs

! Array arguments with intent(in):

  INTEGER (KIND=iintegers), INTENT(IN)    ::  &
    iloc(nvar_lm_chem)          ! location in the look up table

! Array arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)    ::  &
    var_id(nvar_lm)         ! NetCDF-ID of each variable 

  CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
    yerrmsg     ! string for error messages

! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    ierror          ! error index

!------------------------------------------------------------------------------

! Local scalars:
  INTEGER (KIND=iintegers)  :: iid, jid, n, timestep, iztable

  INTEGER (kind=iintegers)  :: my_comm_id, implcode

!
! Local arrays:
!
!- End of header
!==============================================================================

#ifdef NETCDF

IF (npes > 1) THEN
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 1
    RETURN
  ENDIF
ENDIF
IF (npes > 1) THEN
 ! Get id in communicator comm
  CALL MPI_COMM_RANK(icomm,my_comm_id,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_COMM_RANK failed'
    ierror  = 2
    RETURN
  ENDIF
ELSE
  my_comm_id = 0
ENDIF

var_id(:) = 0

IF (my_comm_id == 0) THEN
  ierror=0
  timestep = nstep

! Re-Enter definition mode
  ierror = NF90_redef (ncid)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! loop over all variables that should be written
  DO n = 1, numlist

    iztable = iloc(n)

    ! select netCDF dimension IDs
    IF     (TRIM(var_lm(iztable)%name) == 'U')  THEN
      iid = 9
      jid = 2
    ELSEIF (TRIM(var_lm(iztable)%name) == 'V') THEN
      iid = 1
      jid = 10
    ELSE
      iid = 1
      jid = 2
    ENDIF

!   set the dimensions of the variable regarding to its rank and other criteria
    SELECT CASE (var_lm(iztable)%rank)

    CASE(3)

      IF (UBOUND(var_lm(iztable)%p3,3) == kedim) THEN
          ierror = nf90_def_var(ncid, TRIM(var_lm(iztable)%name), NF90_FLOAT, &
                                  (/ idims_id(iid), idims_id(jid), idims_id(3), idims_id(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
      ELSEIF (UBOUND(var_lm(iztable)%p3,3) == kedim+1) THEN
          ierror = nf90_def_var(ncid, TRIM(var_lm(iztable)%name), NF90_FLOAT, &
                                  (/ idims_id(iid), idims_id(jid), idims_id(4), idims_id(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
      ELSEIF (UBOUND(var_lm(iztable)%p3,3) == ke_soil_lm+1) THEN
          ierror = nf90_def_var(ncid, TRIM(var_lm(iztable)%name), NF90_FLOAT, &
                                  (/ idims_id(iid), idims_id(jid), idims_id(8), idims_id(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
      ELSEIF (UBOUND(var_lm(iztable)%p3,3) == nhori) THEN
          ierror = nf90_def_var(ncid, TRIM(var_lm(iztable)%name), NF90_FLOAT, &
                                  (/ idims_id(iid), idims_id(jid), idims_id(11), idims_id(5) /), var_id(n)) !_br 20.09.2012
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
      ELSE
          ierror = -1
          yerrmsg = 'Invalid 3rd dimension'
          RETURN
      ENDIF

    CASE(2)

        ierror = nf90_def_var(ncid, TRIM(var_lm(iztable)%name), NF90_FLOAT, &
                                (/ idims_id(iid), idims_id(jid), idims_id(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF


    END SELECT

!   set the attributes of the variable
    IF (var_lm(iztable)%standard_name(1:1) /= '-') THEN
      ierror = nf90_put_att (ncid, var_id(n), "standard_name", TRIM(var_lm(iztable)%standard_name))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF
    ierror = nf90_put_att (ncid, var_id(n), "long_name", TRIM(var_lm(iztable)%long_name))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    IF (var_lm(iztable)%units(1:1) /= '-') THEN
      ierror = nf90_put_att (ncid, var_id(n), "units", TRIM(var_lm(iztable)%units))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF
    ierror = nf90_put_att (ncid, var_id(n), "grid_mapping", grid_mapping)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    IF     (TRIM(var_lm(iztable)%name) == 'U') THEN
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "slonu slatu")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ELSEIF (TRIM(var_lm(iztable)%name) == 'V') THEN
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "slonv slatv")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ELSE
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "lon lat")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

    IF  (TRIM(var_lm(iztable)%name) == 'HMO3') THEN
      ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "mole_fraction_of_ozone_in_air: maximum")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

    IF (TRIM(var_lm(iztable)%name) == 'SOILTYP') THEN
      ierror = nf90_put_att (ncid, var_id(n), "flag_values", (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 0. /))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror = nf90_put_att (ncid, var_id(n), "flag_meanings", &
                        "ice rock sand sandy_loam loam clay_loam clay peat sea_water sea_ice")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

    IF  (var_lm(iztable)%lsm == 'l' .OR. var_lm(iztable)%lsm == 'i') THEN  !_br 10.08.10
      ierror = nf90_put_att (ncid, var_id(n), "_FillValue", undefncdf)
      IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
      ENDIF
    ENDIF

  ENDDO     ! End of loop over all variables


! End of definition mode
!!! no more attribute definitions beyond this line !!!

  ierror=nf90_enddef(ncid)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

ENDIF


IF (npes > 1) THEN
  CALL MPI_BCAST(var_id, nvar_lm, imp_integers, 0, icomm, implcode)
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 1
    RETURN
  ENDIF
ENDIF

#endif

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE write_nc_vdefs

!==============================================================================
!==============================================================================

END MODULE src_lm_output
