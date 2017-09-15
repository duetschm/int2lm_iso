!+ Source Module for reading external parameters for LM and the coarse grid
!==============================================================================

MODULE src_read_ext

!==============================================================================
!
! Description:
!   This module contains 2 subroutines for reading the external
!   parameters for LM and the coarse grid. It uses also routines from the 
!   module "io_utilities". Besides it contains the routines for the filtering
!   of the orography.
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
!  Added fields for_e, for_d for reading from external parameters
! 1.3        2005/12/12 Ulrich Schaettler
!  Bug fix: check the table type for reading external parameters
! V1_5         2007/07/09 Ulrich Schaettler
!  Replaced scatter_grib to scatter_values
!  Introduced additional smoothing of orography for very high resolution
!  Renamed GZ0 to Z0
!  Taking into account SLEVE and new type of coding vertical coordinate parameters
!    in get_vert_coord (Daniel Leuenberger)
!  Introduced new field for accepting ECMWF orography on model levels
!  Introduced nlevskip to  skip levels at the top of ECMWF model (Davide Cesari)
! V1_6         2007/09/07 Ulrich Schaettler, Burkhardt Rockel
!  Introduced NetCDF I/O for the external parameters
! V1_7         2007/11/26 Ulrich Schaettler
!  Bug correction in the orography filtering when running sequential program
!  Bug correction in a debug printout
!  Added namelist switches itype_t_cl, itype_ndvi
!  Changed treatment of Z0 in case of netcdf: it is not scaled any more
!  Put computation of actual values for plcov, lai, rootdp to module src_2d_fields
! V1_8         2008/05/29 Ulrich Schaettler, Hans-Juergen Panitz, Burkhardt Rockel
!  Renamed all gz0-variables to z0 and only work with z0-values
!  Read additional external parameters, if lsso or lradtopo
!  Set ID of missing netcdf-variables to -1 to indicate missing parameter
!  Some adaptations in read_nc_gdefs_ext_in
!  Resolution check not done for lcm2lm, because of variable dlat, dlon
!  Compute fis_in also for NetCDF input
!  Check domain size in case of NetCDF input
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
!  Eliminated code in Section 3.4 of read_coarse_grid ext, which was duplicated
! V1_9         2009/09/03 Guenther Zaengl
!  Adaptations for new reference atmosphere
!  Values for zstartlon_ext, zstartlat_ext are now needed in any case
!  (before only for filtering); therefore remove some IFs (Lucio Torrisi)
!  Modifications to read additional external parameters for COSMO-Model
! V1_10        2009/12/17 Ulrich Schaettler
!  Modifications for treating Unified Model data
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler, Burkhardt Rockel, Anne Roches
!  Modifications to use grib_api
!  Modifications to allow processing of JMA and NCEP data
!  Broadcast of irefatm_in in SR read_nc_gdefs_ext_in
!  Added field for salt_lk_lm (BR)
!  Adaptations in SR read_nc_gdefs_ext_lm for reading multi-dimensional external
!    parameter (HORIZON)   (AR)
! V1_17        2011/03/11 Ulrich Schaettler
!  Added lurban flag for reading urban fraction data fr_urban_lm (K. Trusilova)
! V1_19        2012/06/06 Ulrich Schaettler, Helmut Frank, CLM
!  Implemented conditional compilation for NETCDF, GRIBDWD
!  Check whether resolution are in grib data using IBITS
!  Add GFS short names translation. (Helmut Frank)
!  Translate GRIB2 IFS short name FI on hybrid model level 1 to FIS_SH. (Helmut Frank)
!  Extensions for reading UM data
!  Check whether resolution are in grib data using IBITS
!  Adaptations to environment.f90, because of unification with COSMO-Model 4.23
!  CLM
!    Added prescribed surface albedos and new NL switch itype_albedo
!    Change IF clause in lfilter_oro to avoid wrong field indexing
!    Introduce pressure level support for climate model CM input
!    Unified dimension IDs with COSMO (ID for topo corrections changed from 14 to 15)
!    Added some more checks for the external parameter data set
! V1_20        2012/09/03 Ulrich Schaettler, Burkhardt Rockel
!  Renamed 'grax' to 'apix' to be conform with other models
!  Correction in case of east_add_in, west_add_in, south_add_in, north_add_in is not 0 (BR)
! V1_21        2013/03/25 Ulrich Schaettler, Burkhardt Rockel
!  Extensions to read external data sets also with grib_api
!  Changes in netCDF output: scalar variable instead of extra dimension of length 1 (e.g. height_2m)
!    needs changes in idims_id field indices
! V1_22        2013/07/11 Ulrich Schaettler
!  Implemented necessary ifdef GRIBAPI
!  Adapted interface to check_input_grid according to COSMO-Model (pv_in, inrvert_in)
!  Renamed PRS_MIN to RSMIN, which is the official shortName
!  Use structures refatm_out, vcoord_out from vgrid_refatm_utils
!  Renamed grib buffers: ds_grib to ds_grib_single, ds_gribapi to ds_grib_double
!  Adapted interface to read_gribapi with special grib_api integer
! V1_23        2013/10/02 Ulrich Schaettler
!  Rename vcoord_out, refatm_out to vcoord, refatm (as is in COSMO-Model)
! V1_24        2013/11/01 Ulrich Schaettler
!  Moved block for writing names of external parameters after NC-reading
!   (because otherwise some variables are not initialized)
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
  iintegers, & ! KIND-type parameter for standard integer variables
  irealgrib, & ! KIND-type of the REALs in the grib library
  iwlength,  & ! length of an integer word in byte
  intgribf,  & ! KIND-type of the fortran decks in the grib library
  intgribc,  & ! KIND-type of the c decks in the grib library
  int_ga       ! integer precision for grib_api: length of message in bytes

!------------------------------------------------------------------------------

USE data_int2lm_constants, ONLY : &
    G          ! gravity at sea level                          [ms-2]

!------------------------------------------------------------------------------

USE data_fields_lm, ONLY : &
  fis_lm    ,      & ! orography * G                                    (m2/s2)
  hsurf_lm           ! orography                                        (  m  )

!------------------------------------------------------------------------------

USE data_fields_in, ONLY : &
  ps_in    ,       & !
  fis_in    ,      & ! orography * G                                    (m2/s2)
  hsurf_in  ,      & ! orography                                        (  m  )
  hsurfs_in ,      & ! splitted parts of the coarse orography (SLEVE)
  fland_in_tot

!------------------------------------------------------------------------------

USE data_grid_lm,    ONLY: &
  kelm,        & ! ke for LM
  ke1lm,       & ! ke+1
  kedim,       & ! MAX (kelm, ke_in)
  pollat,      & ! latitude of the rotated north pole (in degrees, N>0)
  pollon,      & ! longitude of the rotated north pole (in degrees, E>0)
  polgam,      & ! latitude of the rotated north pole  !_br
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
  startlat,    & ! transformed latitude of the lower left grid point
                 ! of the local domain (in degrees, N>0)
  startlon,    & ! transformed longitude of the lower left grid point
                 ! of the local domain (in degrees, E>0)
  ie2lm_max,   & ! Max. of ie2lm on all processors
  je2lm_max,   & ! Max. of je2lm on all processors
  ie2lm_tot,   & ! = ielm_tot + 2
  je2lm_tot,   & ! = jelm_tot + 2
  ie2lm,       & !
  je2lm,       & !
  ie_ext,      & ! west-east size of fields with external parameters
  je_ext         ! north-south size of fields with external parameters

!------------------------------------------------------------------------------

USE data_grid_in,    ONLY: &
  nlevskip,    & ! number of levels to skip at top of input model
  lcm_pres_coor, & ! CM input data has pressure coordinates  !_br 14.03.12
  jd_min,      & ! smallest index for diamonds for a LM subdomain
  jd_max,      & ! biggest index for diamonds for a LM subdomain
  ni_gme,      & ! resolution of GME
  nd,          & ! number of diamonds (nd = ide-ids+1 = 10)
  ni2,         & ! ni_gme=3**ni3*2**ni2 with ni3 0 or 1 and ni2 > 1
  ni3,         & !
  igg1s,       & ! start index of global array-dimension in x-direction
  igg1sm2,     & ! = igg1s - 2
  igg1e,       & ! end index of global array-dimension in x-direction
  igg1ep2,     & ! = igg1e + 2
  igg2s,       & ! start index of global array-dimension in y-direction
  igg2sm2,     & ! = igg2s - 2
  igg2e,       & ! end index of global array-dimension in y-direction
  igg2ep2        ! = igg2e + 2

USE data_grid_in,    ONLY: &
  pollat_in,      & ! latitude of the rotated north pole (in degrees, N>0)
  pollon_in,      & ! longitude of the rotated north pole (in degrees, E>0)
  polgam_in,      & ! latitude of the rotated north pole
  startlat_in,    & ! transformed latitude of the lower left grid point
                    ! of the local domain (in degrees, N>0)
  startlon_in,    & ! transformed longitude of the lower left grid point
                    ! of the local domain (in degrees, E>0)
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
  ie_in_tot,      & ! ie for input grid, total domain
  je_in_tot,      & ! je for input grid, total domain
  ie_in,          & ! ie for input grid, local domain
  je_in,          & ! je for input grid, local domain
  ke_in,          & ! ke for input grid
  ke1in,          & ! ke+1 for input grid
  east_add_in,    & ! add an extra column to the East
  west_add_in,    & ! add an extra column to the West
  south_add_in,   & ! add an extra column to the South
  north_add_in,   & ! add an extra column to the North
  latitudes_in,   & ! latitudes of the input data
  longitudes_in,  & ! longitudes of the input data
  slatitudes_in,  & ! latitudes of the input data
  slongitudes_in    ! longitudes of the input data

!------------------------------------------------------------------------------

USE data_int2lm_control,    ONLY: &
  noutput,        & ! unit number for output file
  itype_calendar, & ! for specifying the calendar used
  yinput_model,   & ! string to identify the input model
  lgme2lm,        & ! if .TRUE., gme->lm
  lgfs2lm,        & ! if .TRUE., gfs->lm
  lgsm2lm,        & ! if .TRUE., gsm->lm
  lec2lm,         & ! if .TRUE., ec ->lm
  llm2lm,         & ! if .TRUE., lm ->lm
  lum2lm,         & ! if .TRUE., um ->lm
  lhir2lm,        & ! if .TRUE., hirlam ->lm
  lcm2lm,         & ! if .TRUE., cm ->lm
  llake,          & ! if .TRUE., run with lake
  lemiss,         & ! if .TRUE., run with external parameter for surface emissivity
  lstomata,       & ! if .TRUE., run with external parameter for stomata resistance
  lfilter_oro,    & ! if .TRUE., filter the orography
  lxso_first,     & ! if .TRUE., do eXtra smoothing of orography first
  lforest,        & ! if .TRUE., run with forest (evergreen and deciduous)
  lurban,         & ! if .TRUE., run the urban module
  lsso,           & ! process parameters for sso scheme
  lradtopo,       & ! process parameters for topographic correction of radiation
  nhori,          & ! number of sectors for external parameter HORIZON
  itype_t_cl,     & ! to choose origin and treatment of deep soil temperature
  itype_ndvi,     & ! to choose treatment of surface parameters (plcov, lai)
  itype_aerosol,  & ! to choose treatment of surface parameters (plcov, lai)
  itype_albedo,   & ! to choose treatment of solar surface albedo   !_br 10.02.12
  l_topo_z,       & ! if .TRUE., additional smoothing of the topography
                    ! for the LM-Z coordinate Version
  norder_filter,  & ! order of the orography filtering
  eps_filter,     & ! parameter for orography filtering
  ilow_pass_oro,  & ! type of low-pass filter for orography
  numfilt_oro,    & ! number of sequential applications of filter
  ilow_pass_xso,  & ! type of low-pass filter for extra smoothing of steep oro.
  numfilt_xso,    & ! number of sequential applications of xso-filter
  rxso_mask,      & ! mask for extra smoothing of steep oro.: dh > rxso_mask
  rfill_valley,   & ! mask for valley filling: dh > rfill_valley
  ifill_valley,   & ! type of valley filling: 1=MIN (before), 2=MIN (after)
  idbg_level,     & ! to control verbosity of output
  lprintdeb_all     ! whether all PEs print debug output

!------------------------------------------------------------------------------

USE data_int2lm_io,        ONLY : &
  ylmext_cat,        & ! directory of external data for LM
  yinext_cat,        & ! directory of external data for GME
  ylmext_lfn,        & ! file name of external data for LM
  yinext_lfn,        & ! file name of external data for GME
  nuchkdat,          & ! checking the I/O data
  yuchkdat,          & ! checking the I/O data
  ymode_read,        & ! mode for opening the (read) Grib files
  ylmext_form_read,  & ! input format of external LM data
  yinext_form_read,  & ! input format of external boundary data
  npds,              & ! Dimension for product definition section (pds)
  ngds,              & ! Dimension for grid description section (gds)
  nbms,              & ! Dimension for bit map section (bms)
  nbds,              & ! Dimension for binary data section
  ndsup,             & ! Dimension for dsup
  ndims,             & ! Dimension for idims (contains all dimensions)
  lfd,               & !
  nbitmap,           & !
  lds,               & !
  nvar_lm,           & ! maximum number of variables in LM variable table
  nvar_in,           & ! maximum number of variables in GME variable table
  idwdednr,          & ! grib edition number for DWD library
  igrbednr,          & ! grib edition number for grib_api (to be set)
  undefgrib,         & ! value for "undefined" in the grib routines
  undefncdf,         & ! value for "undefined" in the netcdf routines
  undef,             & ! the same with other KIND-Parameter
  lchkin               ! logical for print of check-values (max,min,mean)
                       ! of GME-fields

USE data_int2lm_io,        ONLY : &
  ydate_ini,         & ! start of the forecast yyyymmddhh (year,month,day,hour)
  iblock,            & ! array for gribed data
  idims_in,          & ! array for all dimensions
  ibmap,             & ! array for
  ipds,              & ! product definition section
  igds_in,           & ! grid description section
  ibms,              & ! bit map section
  ibds,              & ! binary data section
  dsup,              & ! Parameter for grib routines
  inrvert_in,        & ! number of vertical coordinate parameters of input data
  pv_in,             & ! array for vertical coordinate parameters of input data
  ds_grib_single,    & ! array for unpacked data
  ds_grib_double,    & ! array for unpacked data
  ar_des_lm,         & ! structure for LM variable table
  ar_des_input,      & ! structure for input variable table
  var_lm,            & ! array for LM variable table
  var_in,            & ! array for GME variable table
  ndims_id,          & ! array for the IDs of the dimensions of netCDF formatted output
  idims_id             ! array for the IDs of the dimensions

!-------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY :  &
    lasync_io,       & ! if .TRUE.: the model runs with extra PEs for
                       ! asynchronous IO
    num_compute,     & ! number of compute PEs
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    isubpos,         & ! positions of the subdomains in the total domain. Given
                       ! are the i- and the j-indices of the lower left and the
                       ! upper right grid point in the order
                       !                  i_ll, j_ll, i_ur, j_ur.
                       ! Only the interior of the domains are considered, not
                       ! the boundary lines.
    isubpos_coarse,  & ! positions of the subdomains in the total domain.
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    icomm_cart,      & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen,     & ! length of one column of sendbuf
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

USE environment,         ONLY :  extend_field
USE utilities,           ONLY :  horizontal_filtering, sleve_split_oro

USE mpe_io,              ONLY :  mpe_io_read
USE parallel_utilities,  ONLY :  remark, distribute_values,                   &
                                 gather_values, scatter_values
USE gme_utilities,       ONLY :  xd
USE io_utilities,        ONLY :  open_file, read_grib, read_gribapi,          &
                                 read_netcdf, close_file,                     &
                                 check_record, check_input_grid

USE vgrid_refatm_utils,  ONLY :  vcoord, refatm, refatm_in, vcoord_in,        &
                                 vcoord_d, rundefined, nfltvc_in, svc1_in,    &
                                 svc2_in

#ifdef NETCDF
! NetCDF routines and declarations
USE netcdf,           ONLY :   &
  nf90_inq_dimid,              &
  nf90_inquire_dimension,      &
  nf90_inq_varid,              &
  nf90_get_var,                &
  nf90_get_att,                &
  nf90_noerr,                  &
  nf90_strerror,               &
  nf90_enotvar,                &
  nf90_enotatt
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! string variable to hold grid information
  CHARACTER (LEN=200)     grid_mapping

!==============================================================================

PUBLIC  read_lm_ext, read_coarse_grid_ext

!==============================================================================

CONTAINS

!==============================================================================
!+ Reads the external parameters for LM
!------------------------------------------------------------------------------

SUBROUTINE read_lm_ext (lhsur_lm, lfis__lm, lfrla_lm, lz0___lm, lz012_lm,  &
                        lsoty_lm, lroot_lm, lplmx_lm, lplmn_lm, lpl12_lm,  &
                        laimx_lm, laimn_lm, lai12_lm, lfore_lm, lford_lm,  &
                        lurba_lm, laldr_lm, lalsa_lm,                      &
                        lemis_lm, lprsm_lm, lflak_lm, ldept_lm, ltcl__lm,  &
                        lndvi_lm, lstdh_lm, lgamm_lm, lthet_lm, lsigm_lm,  &
                        lskyv_lm, lsang_lm, lsasp_lm, lhori_lm, lsalt_lm,  &
                        lsu12_lm, ldu12_lm, lor12_lm, lbc12_lm, lss12_lm,  &
                        lal12_lm,                                          &
                        ierror, yerror)

!------------------------------------------------------------------------------
!
! Description:
!  Reads the fields of external parameters for the LM from a GRIB1 
!  file, unpacks and saves them in the appropriate arrays.
!  The organization of this routine is such, that it works on parallel and
!  on sequential platforms. The strategy has been taken over from the LM.
!
! Method:
!  After creating the filename and opening the grib file, all records are read
!  in a read_loop. Every processor gets a total record for de-gribing (up to 
!  the maximum number of records or processors). 
!  In `distribute_loop' (loop over processors), the records are distributed
!  according to the domain decomposition.
!
! Input files:
!  Grib-file with external parameters for LM. 
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
LOGICAL,     INTENT (INOUT)   ::                                 &
  ! logicals, to indicate whether the corresponding field has been read
  lhsur_lm, lfis__lm, lfrla_lm, lz0___lm, lz012_lm,              &
  lsoty_lm, lroot_lm, lplmx_lm, lplmn_lm, lpl12_lm,              &
  laimx_lm, laimn_lm, lai12_lm, lfore_lm, lford_lm,              &
  lurba_lm, laldr_lm, lalsa_lm,                                  &
  lemis_lm, lprsm_lm, lflak_lm, ldept_lm, ltcl__lm,              &
  lndvi_lm, lstdh_lm, lgamm_lm, lthet_lm, lsigm_lm,              &
  lskyv_lm, lsang_lm, lsasp_lm, lhori_lm, lsalt_lm,              &
  lsu12_lm, ldu12_lm, lor12_lm, lbc12_lm, lss12_lm,              &
  lal12_lm

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER  (KIND=iintegers), PARAMETER  ::  &
  nlistextparmax = 40

INTEGER  (KIND=intgribf)   ::  &
  ierrf, i_bits,        & ! error code for grib routines
  igribid, ireturn,     & ! grib and error handler for grib_api files
  ilevel, itoplevel, ibottomlevel, igriblen, ijincr, &
  ke1lmf                  ! ke1lm with KIND parameter intgribf

INTEGER  (KIND=iintegers)  ::  &
  izerror, nzbyte,                         & ! status and error status variable
  izloc, ildslm, ilfdlm,                   & ! location of field in variable table
  iz_rsize, iz_lfd, iz_info(7),            & ! characteristics of read records
  izdebug,                                 & ! for verbosity of output
  nufile,                                  & ! unit number of opened grib file
  numlistextpar,                           & ! number of external fields
  nzid_frland, nzid_hsurf, nzid_fis,       & ! ID of PEs that get hsurf or fis
  izproc, izshift, jzshift, ize, jze,      & !
  i, j, ij, k, n, nit, ndim, izstat, nfrl, & ! additional variables
  ile, iri, jlo, jup, izvar_count, izlev, izrank, izlev_count, imonth

INTEGER  (KIND=iintegers)  ::  &
  myzvar, myzlev, myzlevtot      ! organization indices returned by write_netcdf

INTEGER                    ::  &
  ivar_id(nlistextparmax)        ! ID of the external variable

LOGICAL                    :: &
  lcheckin_ext(nlistextparmax)   ! list for checking which variables have been read

LOGICAL                    :: &
  ! for checking multi-dimensional variables
  lzndvir(12), lzaersu(12), lzaerdu(12), lzaeror(12), lzaerbc(12), lzaerss(12), &
  lzplcov(12), lzlai(12), lzz0(12), lzalb(12)

! Variables for checking the grid
INTEGER  (KIND=iintegers)  ::  &
  izerrorgrid,                         & ! error status variable
  izsenderror(0:num_compute-1),        & ! for exchange of error status
  izrecverror(0:num_compute-1)           ! for exchange of error status

REAL (KIND=ireals)         ::  &
  zlatll, zlonll, zlatur, zlonur, zdlon, zdlat, zbias, zfactor,      &
  zstartlon_ext, zstartlat_ext

REAL    (KIND=ireals)      ::  &
  zmin, zmax, zsum, zmean, rbuf(2), & ! values of min, max and meanvalue
  zpollat_grib, zpollon_grib, zsouthp_lat, zsouthp_lon, zhmax_sea

CHARACTER (LEN=200)        ::   &
  yline,                        & ! line to be written for one processor
  ylines(num_compute)             ! array of lines to collect all lines from
                                  ! the other processors

INTEGER (KIND=iintegers)   ::   &
  iline (200),                  & ! the same for the integer values
  ilines(200,num_compute)

REAL (KIND=ireals), ALLOCATABLE         ::            &
  field_filt (:,:),                                   &
  field_ori  (:,:),                                   &
  frland_filt(:,:),                                   &
  xy_vec(:), a(:), b(:), c(:), d(:), e(:), f(:),      &
              ci(:), cj(:), ck(:), cl(:), cm(:)

REAL    (KIND=ireals)      ::  &
  zak(ke1lm), zbk(ke1lm), zhhlr(ke1lm), zpxx, zhh1, zhh2, zhh3, zr_d

REAL    (KIND=irealgrib)   ::  &
  zundef

REAL (KIND=irealgrib), ALLOCATABLE      ::            &
  zprocarray (:,:),                                   &
  zsubarray  (:),                                     &
  lm_field   (:,:)

LOGICAL                    ::  &
  lopen,     & ! test whether file YUCHKDAT is opened
  lexist,    & ! to check existence of LM external parameter file
  lrequired, & ! indicates whether a record from the grib file is required
  leof,      & ! indicates the end of file
  lread,     & ! go on with reading records
  linit,     & ! for the first cycle of read_loop
  lzgotfrland, lzgothsurf, lzgotfis, & ! if hsurf or fis were read 
  lzdisfrland, lzdishsurf, lzdisfis    ! if hsurf or fis were distributed

CHARACTER (LEN=10)         ::  &
  ylistextpar (nlistextparmax)     ! list of external parameters

CHARACTER (LEN=250)        ::  &
  yname      ! name of the grib file (is determined in create_filename)
CHARACTER (LEN= 25)        ::  &
  yroutine   ! name of this routine for error handling

CHARACTER (LEN= 30)        ::  &
  ytypeoflevel ! typeOfLevel

CHARACTER (LEN= 21)        ::  &
  yshortname, yzname  ! short name from grib_api

CHARACTER (LEN=14)         ::  &
  ydate        ! date in the form   yyyymmddhhmmss

! ---- Filling of deep V-valleys     and    ----
! ---- Filtering of orography (2nd version) ----
LOGICAL ::  &
  ldhsurf_xy

INTEGER (KIND=iintegers) ::  &
  hfwidth, hfw_m_nb,         &
  ie_ext_hf, je_ext_hf,      &
  hfjstartpar, hfjendpar,    &
  istata, istatd,            &
  number_mask_points,        &
  dummy_neigh(4),            &
  kzdims(24)

INTEGER (KIND=int_ga)    ::  &
  izmaxlen, iz_rsize_ga

REAL (KIND=ireals) ::   &
  zak_ke, zbk_ke,       &
  zdh_x1, zdh_x2,       &
  zdh_y1, zdh_y2,       &
  zdh_xy1, zdh_xy2,     &
  zdh_xy3, zdh_xy4,     &
  zdh_max,              &
  zff_x, zff_y, zt00

LOGICAL, ALLOCATABLE ::  &
  hfx_mask(:,:), hfy_mask(:,:)

REAL (KIND=ireals), ALLOCATABLE ::  &
  ff_tmp(:,:)

CHARACTER (LEN=200)        ::  &
  yerrmsg    ! error message for error handling

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! first initializations
  ierror      = 0_iintegers
  yerror      = '     '
  igriblen    = 0_intgribf

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  IF   ( (ylmext_form_read == 'grb1') .OR. (ylmext_form_read == 'apix') ) THEN
    zundef = undefgrib
  ELSEIF (ylmext_form_read == 'ncdf') THEN
    zundef = undefncdf
  ENDIF
  undef       = REAL (zundef, ireals)

  ierrf       = 0_intgribf
  ke1lmf      = INT (ke1lm, intgribf)

  ! initialize level counter and logical flags
  lcheckin_ext(:) = .FALSE.
  leof        = .FALSE.
  lexist      = .FALSE.
  linit       = .TRUE.
  izerrorgrid = 0_iintegers
  izstat      = 0_iintegers
  yroutine    = 'read_lm_ext'

  ! initializations for the filtering
  lzgotfrland = .FALSE.
  lzgothsurf  = .FALSE.
  lzgotfis    = .FALSE.
  lzdisfrland = .FALSE.
  lzdishsurf  = .FALSE.
  lzdisfis    = .FALSE.
  nzid_frland = -1
  nzid_hsurf  = -1
  nzid_fis    = -1

  ! Set list of external parameters that should be read
  numlistextpar   = 10
  ylistextpar( 1) = 'HSURF     '
  ylistextpar( 2) = 'FIS       '
  IF     (itype_ndvi /= 2) THEN
    ylistextpar( 3) = 'Z0        '
  ELSEIF (itype_ndvi == 2) THEN
    ylistextpar( 3) = 'Z012      '
  ENDIF
  ylistextpar( 4) = 'FR_LAND   '
  ylistextpar( 5) = 'SOILTYP   '
  ylistextpar( 6) = 'ROOTDP    '

  IF     (itype_ndvi /= 2) THEN
    ylistextpar( 7) = 'PLCOV_MX  '
    ylistextpar( 8) = 'PLCOV_MN  '
    ylistextpar( 9) = 'LAI_MX    '
    ylistextpar(10) = 'LAI_MN    '
    numlistextpar   = 10
    IF (itype_ndvi == 1) THEN
      numlistextpar   = numlistextpar + 1
      ylistextpar(numlistextpar) = 'NDVI_MRAT '
    ENDIF
  ELSEIF (itype_ndvi == 2) THEN
    ylistextpar( 7) = 'PLCOV12   '
    ylistextpar( 8) = 'LAI12     '
    numlistextpar   = 8
  ENDIF

  IF (lforest) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'FOR_E     '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'FOR_D     '
  ENDIF

  IF (lurban) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'URBAN     '
  ENDIF

  IF (lemiss) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'EMIS_RAD  '
  ENDIF

  IF     (itype_albedo == 2) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'ALB_DRY   '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'ALB_SAT   '
  ELSEIF (itype_albedo == 3) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'ALB_DIF12 '
  ENDIF

  IF (lstomata) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'RSMIN     '
  ENDIF

  IF (lsso) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'SSO_STDH  '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'SSO_GAMMA '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'SSO_THETA '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'SSO_SIGMA '
  ENDIF

  IF (lradtopo) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'SKYVIEW   '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'SLO_ANG   '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'SLO_ASP   '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'HORIZON   '
  ENDIF

  IF (llake) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'FR_LAKE   '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'DEPTH_LK  '
    ! added by CLM
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'SALT_LK   '
  ENDIF

  IF (itype_aerosol == 2) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'AER_SO412 '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'AER_DUST12'
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'AER_ORG12 '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'AER_BC12  '
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'AER_SS12  '
  ENDIF

  IF (itype_t_cl == 1) THEN
    numlistextpar   = numlistextpar + 1
    ylistextpar(numlistextpar) = 'T_CL      '
  ENDIF

  ! Grib initializations
  ! initialize variable iwlength
  i_bits      = BIT_SIZE (ierrf)
  iwlength    = INT (i_bits, iintegers) / 8

  ! Set dimensions for grib variables and allocate iblock, ibmap, ds and dsup
  ! nzbyte is assumed to be 2 here: this is true in all our cases, but if
  ! if we once change the packing rate of grib code (nrbit=16) we will get
  ! in trouble here. To know this a priori, the pds of the first record has
  ! to be read and decoded to get nrbit. The 2000 are just a safety-add to
  ! take care of the definition sections that are also stored in iblock
  nzbyte = 2

  ildslm  = ie_ext * je_ext
  lds     = INT (ildslm, intgribf)
  ilfdlm  = ildslm  * nzbyte / iwlength + 2000
  lfd     = INT (ilfdlm, intgribf)
  iz_lfd  = INT (lfd , iintegers)

  ! dimensions for grib routines
  idims_in( 1)   = npds
  idims_in( 2)   = ngds
  idims_in( 3)   = nbms
  idims_in( 4)   = nbds
  idims_in( 5)   = nbitmap
  idims_in( 6)   = ndsup
  idims_in( 7)   = lds
  idims_in( 8)   = lfd
  idims_in(9:20) = 0

  ! Allocate fields
  ALLOCATE (iblock(lfd), ibmap(nbitmap), dsup(ndsup),            STAT=izstat)
  ALLOCATE (ds_grib_single(lds), ds_grib_double(ie_ext*je_ext),  STAT=izstat)

  ALLOCATE (zprocarray (ie2lm_max* je2lm_max, num_compute),     &
            zsubarray  (ie2lm_max* je2lm_max),                  &
            lm_field   (ie2lm_tot, je2lm_tot), STAT = izstat)

  IF (lfilter_oro) THEN
    ALLOCATE (field_filt  (ie_ext,je_ext),                      &
              frland_filt (ie_ext,je_ext), STAT = izstat)
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Create file name and open the file
!------------------------------------------------------------------------------

  ! Construct the file name
  yname = TRIM(ylmext_cat)//TRIM(ylmext_lfn)

  ! because mpe_io does not return error codes, the existence of that special
  ! file has to be checked before calling open_file
  IF (my_cart_id == 0) THEN
    INQUIRE (FILE=yname, EXIST=lexist)
  ENDIF
  IF (num_compute > 1) THEN
    CALL distribute_values (lexist, 1, 0, imp_logical, icomm_cart, izerror)
  ENDIF

  ! All processors have to call the routine open_file. What the parallel
  ! program really does is determined in the routine.
  IF (lexist) THEN
    CALL open_file(nufile, yname, ymode_read, ylmext_form_read, icomm_cart, &
                   my_cart_id, num_compute, lasync_io, idbg_level,          &
                   yerrmsg, izerror)
    IF (izerror /= 0) THEN
      lread = .FALSE.
      ierror = 1
      yerror = 'ERROR: *** External parameter file for LM could not be opened ***'
      RETURN
    ELSE
      lread = .TRUE.
    ENDIF
  ELSEIF (ylmext_lfn == 'interpolate') THEN
    ! external parameters should be interpolated from coarse grid fields
    CALL remark (my_cart_id, 'read_lm_ext',                         &
                         'Get external parameters for LM from coarse grid')
    lread = .FALSE.
  ELSE
    ! file could not be opened; set error code and exit
    ierror = 1
    yerror = 'ERROR: *** External parameter file for LM could not be opened ***'
    RETURN
  ENDIF

  IF (ylmext_form_read == 'ncdf' .AND. lread) THEN
    ! read global attributes and definitions for NetCDF
    CALL read_nc_gdefs_ext_lm (nufile, ie_ext, je_ext,                       &
           startlat_tot, startlon_tot, endlat_tot, endlon_tot,               &
           dlon, dlat, pollon, pollat, polgam, zstartlon_ext, zstartlat_ext, &
           icomm_cart, my_cart_id, num_compute, yerrmsg, izerror)
    IF (izerror /= 0) THEN
      ierror = 2
      yerror = 'Error in read_nc_gdefs_ext_lm'
      PRINT *, TRIM(yerrmsg) !_br 06.04.09
      RETURN
    ENDIF
    izshift = NINT ( (startlon_tot - dlon - zstartlon_ext)/dlon )
    jzshift = NINT ( (startlat_tot - dlat - zstartlat_ext)/dlat )

    CALL read_nc_vdefs_ext_lm (nufile, var_lm, nvar_lm, ivar_id,             &
           pollon, pollat, numlistextpar, ylistextpar, lcheckin_ext,         &
           icomm_cart, my_cart_id, num_compute, yerrmsg, izerror)
    IF (izerror /= 0) THEN
      ierror = 3
      yerror = 'read_nc_vdefs_ext_lm'
      RETURN
    ENDIF

    IF (izdebug > 10) THEN
      PRINT *, 'Necessary external parameters for LM grid:'
      DO n = 1, numlistextpar
        PRINT *, ylistextpar(n), ivar_id(n), lcheckin_ext(n)
      ENDDO
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 3: (Endless) loop over all records in the grib file
!------------------------------------------------------------------------------
 
  IF (izdebug >= 10) THEN
    PRINT *, '   Start endless loop to read records: ', ylmext_form_read
  ENDIF

  ! for NetCDF files this is nevertheless done in an ordered manner by going
  ! through all variables in a list. Some organizational variables have to
  ! be set for this:
  izvar_count = 1
  izlev_count = 0

  read_loop: DO WHILE (lread)

  !----------------------------------------------------------------------------
  ! 3.1: Get a record
  !----------------------------------------------------------------------------

    ! Every PE gets one record from the file in an ordered manner
    ! (determined by the rank in the communicator icomm_rank). How this is 
    ! done exactly is determined in the routine read_"format". This routine has
    ! to be called by every PE.

    SELECT CASE (ylmext_form_read)

    CASE ('grb1')

#ifdef GRIBDWD
      IF (izdebug >= 20) THEN
        PRINT *, '      Calling read_grib'
      ENDIF

      CALL read_grib    (nufile, iz_lfd*iwlength, iz_lfd, icomm_cart,     &
                         iblock, iz_rsize, num_compute, lasync_io, yerrmsg, izerror)

      IF (idbg_level >= 20) THEN
        IF (iz_rsize == 0) THEN
          PRINT *, '       EOF reached'
        ELSE
          PRINT *, '       Read a record with size ', iz_rsize
        ENDIF
      ENDIF
#endif

    CASE ('apix')

#ifdef GRIBDWD
      IF (izdebug >= 20) THEN
        PRINT *, '      Calling read_gribapi'
      ENDIF

      izmaxlen  = INT (iz_lfd*iwlength, int_ga)
      CALL read_gribapi (nufile, izmaxlen, iz_lfd, icomm_cart,     &
                iblock, iz_rsize_ga, num_compute, lasync_io, yerrmsg, izerror)

      IF (iz_rsize_ga == 0) THEN
        iz_rsize = 0_iintegers
        IF (idbg_level >= 20) THEN
          PRINT *, '       EOF readched'
        ENDIF
      ELSE
        iz_rsize = 1_iintegers ! just to indicate > 0
        IF (idbg_level >= 20) THEN
          PRINT *, '       Read a record with size ', iz_rsize_ga
        ENDIF
      ENDIF
#endif

    CASE ('ncdf')

      CALL read_netcdf  (nufile, ie_ext, je_ext, izvar_count,                 &
                         izlev_count, ivar_id, numlistextpar, idims_id,       &
                         ndims_id, icomm_cart, my_cart_id, num_compute,       &
                         0, 0, 0, 0,                                          &
                         ds_grib_single, myzvar, myzlev, myzlevtot, iz_rsize, &
                         imp_grib, lasync_io, yerrmsg, izerror)
      IF (idbg_level >= 10) THEN
        IF (iz_rsize > 0) THEN
          PRINT *, '       Got record ', myzvar, myzlev, myzlevtot, ivar_id(myzvar),   &
                           ylistextpar(myzvar), iz_rsize, ds_grib_single(1),           &
                           izvar_count, izlev_count
        ELSE
          PRINT *, '       Got no more record: EOF '
        ENDIF
      ENDIF

    END SELECT

    IF (izerror /= 0) THEN
       ierror = 2
       yerror = 'Error while reading file'
       RETURN
    ENDIF

    IF (iz_rsize == 0) THEN
      ! this PE has got no more record because the end of file is reached
      leof = .TRUE.
    ENDIF
 
  !----------------------------------------------------------------------------
  ! 3.2: If this PE has got a record, de-grib it and put it to zprocarray
  !----------------------------------------------------------------------------

    IF (.NOT. leof) THEN 

      !------------------------------------------------------------------------
      ! 3.2.1: For GRIB data: de-grib the record
      !------------------------------------------------------------------------

      IF (ylmext_form_read == 'grb1') THEN
#ifdef GRIBDWD
        igrbednr = 1

        ! ds has to have the same precision as the REALs in the grib library
        ds_grib_single = 0.0_irealgrib

        ! de-grib the data with the routine grbin1 from the grib library
        CALL grbin1(idwdednr, undefgrib, ndims, idims_in, iblock, ibmap,     &
                    ipds, igds_in, ibms, ibds, dsup, ds_grib_single, ierrf)
        IF (ierrf /= 0) THEN
          ierror = 3
          yerror = 'Error in grbin1'
          RETURN
        ENDIF

        IF (idbg_level >= 20) THEN
          IF (iz_rsize == 0) THEN
            PRINT *, '       EOF reached'
          ELSE
            PRINT *, '       Got record ', iz_rsize, ipds(2), ipds(7), ipds(8)
          ENDIF
        ENDIF
#endif

      ELSEIF (ylmext_form_read == 'apix') THEN

#ifdef GRIBAPI
        ! Build the grib handle
        CALL grib_new_from_message (igribid, iblock, ireturn)
        IF (ireturn /= GRIB_SUCCESS) THEN
          PRINT *,  ' *** Error in grib_api grib_new_from_message  ', ireturn
          ierror = 5
          yerror = ' *** Error in grib_api grib_new_from_message'
          RETURN
        ENDIF

        ! get edition number
        CALL grib_get (igribid, 'editionNumber',   igrbednr, ireturn)

        ! Get some keys
        CALL grib_get (igribid, 'shortName',     yshortname,    ireturn)

        ! Get size of data and data itself
        ! Set missing_value before
        CALL grib_set (igribid, 'missingValue', undefgrib)
        CALL grib_get_size (igribid, 'values', igriblen, ireturn)

        IF (igriblen > lds) THEN
          PRINT *, ' *** ERROR: size of message is too big for allocated field: ', igriblen, lds
          ierror = 6
          yerror = ' *** ERROR: grib_api: Wrong size of field ***'
          RETURN
        ENDIF

        CALL grib_get (igribid, 'values', ds_grib_single, ireturn)

        IF (ireturn /= GRIB_SUCCESS) THEN
          PRINT *,   ' *** Error in grib_api grib_get values', ireturn
          ierror = 6
          yerror = ' *** Error in grib_api grib_get values'
          RETURN
        ENDIF
#endif

      ENDIF

      !------------------------------------------------------------------------
      ! 3.2.2: Preparations for filtering the orography
      !------------------------------------------------------------------------

      IF (lfilter_oro) THEN
        ! Check whether hsurf or fis are read and keep (one of) them for 
        ! filtering, if desired. In case of filtering frland has also to
        ! be kept

        ! Fraction of Land
        IF     (ylmext_form_read == 'grb1') THEN
          IF (ipds(7) == 81) THEN
            lzgotfrland = .TRUE.
          ELSE
            lzgotfrland = .FALSE.
          ENDIF
        ELSEIF (ylmext_form_read == 'apix') THEN
          IF (TRIM(yshortname) == 'FR_LAND') THEN
            lzgotfrland = .TRUE.
          ELSE
            lzgotfrland = .FALSE.
          ENDIF
        ELSEIF (ylmext_form_read == 'ncdf') THEN
          IF (ylistextpar(myzvar) == 'FR_LAND   ') THEN
            lzgotfrland = .TRUE.
          ELSE
            lzgotfrland = .FALSE.
          ENDIF
        ENDIF

        IF (lzgotfrland) THEN
          ! preparations to keep frland
          nzid_frland = my_cart_id
          DO j = 1, je_ext
            DO i = 1, ie_ext
              ij = (j-1)*ie_ext + i
              frland_filt(i,j) = REAL (ds_grib_single(ij), ireals)
            ENDDO
          ENDDO
        ELSE
          nzid_frland = -1
        ENDIF

        ! Height of orography
        IF     (ylmext_form_read == 'grb1') THEN
          IF (ipds(7) == 8) THEN
            lzgothsurf  = .TRUE.
          ELSE
            lzgothsurf  = .FALSE.
          ENDIF
        ELSEIF (ylmext_form_read == 'apix') THEN
          IF (TRIM(yshortname) == 'HSURF') THEN
            lzgothsurf  = .TRUE.
          ELSE
            lzgothsurf  = .FALSE.
          ENDIF
        ELSEIF (ylmext_form_read == 'ncdf') THEN
          IF (ylistextpar(myzvar) == 'HSURF     ') THEN
            lzgothsurf  = .TRUE.
          ELSE
            lzgothsurf  = .FALSE.
          ENDIF
        ENDIF

        IF (lzgothsurf) THEN
          ! preparations to keep hsurf
          nzid_hsurf = my_cart_id
          DO j = 1, je_ext
            DO i = 1, ie_ext
              ij = (j-1)*ie_ext + i
              field_filt(i,j) = REAL (ds_grib_single(ij), ireals)
            ENDDO
          ENDDO
        ELSE
          nzid_hsurf = -1
        ENDIF

        ! Geopotential of surface
        IF     (ylmext_form_read == 'grb1') THEN
          IF ((ipds(7) == 6) .AND. (.NOT. lzdishsurf) ) THEN
            lzgotfis    = .TRUE.
          ELSE
            lzgotfis    = .FALSE.
          ENDIF
        ELSEIF (ylmext_form_read == 'apix') THEN
          IF ( (TRIM(yshortname) == 'FIS') .AND. (.NOT. lzdishsurf) ) THEN
            lzgotfis    = .TRUE.
          ELSE
            lzgotfis    = .FALSE.
          ENDIF
        ELSEIF (ylmext_form_read == 'ncdf') THEN
          IF ((ylistextpar(myzvar) == 'FIS       ') .AND. (.NOT. lzdishsurf) ) THEN
            lzgotfis    = .TRUE.
          ELSE
            lzgotfis    = .FALSE.
          ENDIF
        ENDIF

        IF ((.NOT. lzdishsurf) .AND. (lzgotfis)) THEN
          ! If hsurf has been read already, fis is skipped, because only one 
          ! variable is filtered and the other is computed from the filtered
          ! one. But it must be notified, that fis has been read
          nzid_fis   = my_cart_id

          ! preparations to keep fis  
          DO j = 1, je_ext
            DO i = 1, ie_ext
              ij = (j-1)*ie_ext + i
              field_filt(i,j) = REAL (ds_grib_single(ij), ireals)
            ENDDO
          ENDDO
        ELSE
          nzid_fis   = -1
        ENDIF
      ENDIF

      !------------------------------------------------------------------------
      ! 3.2.3: Check the grid and get some sizes (for GRIB only)
      !------------------------------------------------------------------------

      IF ( (ylmext_form_read == 'grb1') .OR. (ylmext_form_read == 'apix') ) THEN

        ! Check size and resolution of the external parameters against the 
        ! LM grid and whether the LM domain (with 1 boundary line) is within
        ! these fields
        !-------------------------------------------------------------------
       
        IF     (ylmext_form_read == 'grb1') THEN
#ifdef GRIBDWD
          ize    = igds_in ( 5)   ! ie-size of the external parameter field
          jze    = igds_in ( 6)   ! je-size of the external parameter field
          zlatll = REAL (igds_in ( 7), ireals)*0.001_ireals ! lower left latitude
          zlonll = REAL (igds_in ( 8), ireals)*0.001_ireals ! lower left longitude
          zlatur = REAL (igds_in (10), ireals)*0.001_ireals ! upper right latitude
          zlonur = REAL (igds_in (11), ireals)*0.001_ireals ! upper right longitude
  
          ! Check the rotated pole: Instead of recomputing the given grib
          ! coordinates of the rotated south pole to the north pole, the namelist
          ! coordinates of the north pole are computed to the south pole
          zpollat_grib =   REAL(igds_in(20), ireals)*0.001_ireals
          zpollon_grib =   REAL(igds_in(21), ireals)*0.001_ireals
#endif

        ELSEIF (ylmext_form_read == 'apix') THEN

#ifdef GRIBAPI
          CALL grib_get (igribid, 'Ni',                                    ize,  ireturn)
          CALL grib_get (igribid, 'Nj',                                    jze,  ireturn)

          CALL grib_get (igribid, 'latitudeOfFirstGridPointInDegrees',  zlatll,  ireturn)
          CALL grib_get (igribid, 'longitudeOfFirstGridPointInDegrees', zlonll,  ireturn)
          CALL grib_get (igribid, 'latitudeOfLastGridPointInDegrees',   zlatur,  ireturn)
          CALL grib_get (igribid, 'longitudeOfLastGridPointInDegrees',  zlonur,  ireturn)

          CALL grib_get (igribid, 'latitudeOfSouthernPoleInDegrees',   zpollat_grib, ireturn)
          CALL grib_get (igribid, 'longitudeOfSouthernPoleInDegrees',  zpollon_grib, ireturn)
#endif

          ! Careful: for grib2, longitudes are between    0...360 degrees,
          !          for grib1, longitudes are between -180...180
          ! the COSMO convention still is -180...180
          IF (igrbednr == 2) THEN
            ! convert to COSMO convention
            IF (zlonll >= 180.0_ireals) THEN
              zlonll = zlonll - 360.0_ireals
            ENDIF
            IF (zlonur>  180.0_ireals) THEN
              zlonur = zlonur - 360.0_ireals
            ENDIF
          ENDIF

        ENDIF

        zsouthp_lat  =   - pollat
        zsouthp_lon  =   pollon + 180.0_ireals
        IF (zsouthp_lon > 180.0_ireals) THEN
          zsouthp_lon = zsouthp_lon - 360.0_ireals
        ENDIF

        ! set zstartlon_ext and zstartlat_ext
        zstartlon_ext = zlonll
        zstartlat_ext = zlatll
  
        IF     (ylmext_form_read == 'grb1') THEN
#ifdef GRIBDWD
          ijincr =  IBITS(igds_in(9),7,1)
#endif
        ELSEIF (ylmext_form_read == 'apix') THEN
#ifdef GRIBAPI
          ! for the increments check the resolution and component flag
          CALL grib_get (igribid, 'ijDirectionIncrementGiven', ijincr,  ireturn)
#endif
        ENDIF

        IF (ijincr == 0) THEN
          ! no direction increments are given
          IF (zlonur - zlonll < 0.0_ireals) THEN
            zdlon  = (zlonur - zlonll + 360.0_ireals) / (ize-1)
          ELSE
            zdlon  = (zlonur - zlonll) / (ize-1)
          ENDIF
          zdlat  = (zlatur - zlatll) / (jze-1)
        ELSE
          ! direction increments are given
          IF     (ylmext_form_read == 'grb1') THEN
#ifdef GRIBDWD
            zdlon  = REAL (igds_in(12), ireals)*0.001_ireals ! resolution for lambda
            zdlat  = REAL (igds_in(13), ireals)*0.001_ireals ! resolution for phi
#endif
          ELSEIF (ylmext_form_read == 'apix') THEN
#ifdef GRIBAPI
            CALL grib_get (igribid, 'Di', zdlon,  ireturn)
            CALL grib_get (igribid, 'Dj', zdlat,  ireturn)
#endif
          ENDIF
        ENDIF
  
        ! Do all the checks:  abort program in case of errors
        IF ( (ize /= ie_ext) .OR. (jze /= je_ext) )             izerrorgrid = 1
        IF ( (zlonll > startlon_tot - dlon) .OR.         &
             (zlatll > startlat_tot - dlat)  )                  izerrorgrid = 2
        IF ( (zlonur < endlon_tot   + dlon) .OR.         &
             (zlatur < endlat_tot   + dlat)  )                  izerrorgrid = 3
        IF ( (ABS(zdlon - dlon) > 1.0E-5_ireals) .OR.    &
             (ABS(zdlat - dlat) > 1.0E-5_ireals)      )         izerrorgrid = 4
        IF ( (ABS(zpollat_grib - zsouthp_lat) > 1.0E-5_ireals) .OR. &
             (ABS(zpollon_grib - zsouthp_lon) > 1.0E-5_ireals) )izerrorgrid = 5
        izshift = NINT((startlon_tot - zstartlon_ext)/dlon)
        jzshift = NINT((startlat_tot - zstartlat_ext)/dlat)
        IF (ABS(zstartlon_ext + izshift*zdlon - startlon_tot) > 0.0005_ireals)&
                                                                izerrorgrid = 6
        IF (ABS(zstartlat_ext + jzshift*zdlat - startlat_tot) > 0.0005_ireals)&
                                                                izerrorgrid = 7
      ENDIF
    ENDIF   ! leof

    IF ( (ylmext_form_read == 'grb1') .OR. (ylmext_form_read == 'apix') ) THEN

      !------------------------------------------------------------------------
      ! 3.2.4: Distribute zstartlat_ext and zstartlon_ext to all PEs
      !------------------------------------------------------------------------

      IF (num_compute > 1) THEN
        ! can be done by PE 0, because it got a record at least
        IF (my_cart_id == 0) THEN
          rbuf(1) = zstartlat_ext
          rbuf(2) = zstartlon_ext
        ENDIF
        CALL distribute_values (rbuf, 2, 0, imp_reals, icomm_cart, izerror)
        IF (my_cart_id /= 0) THEN
          zstartlat_ext = rbuf(1)
          zstartlon_ext = rbuf(2)
        ENDIF
      ENDIF

      !------------------------------------------------------------------------
      ! 3.2.5: Exchange izerrorgrid and print error message
      !------------------------------------------------------------------------

      izsenderror (my_cart_id) = izerrorgrid
      IF (num_compute > 1) THEN
        CALL gather_values (izsenderror(my_cart_id), izrecverror, 1, num_compute,&
                            imp_integers, -1, icomm_cart, yerrmsg, izerror)
      ELSE
        izrecverror(:) = izsenderror(:)
      ENDIF

      ! If necessary, print errors
      DO izproc = 0, num_compute-1
        IF (izrecverror(izproc)/= 0) THEN
          ! only one processor prints the error message
          IF (my_cart_id == izproc) THEN
            WRITE (*, '(     T26, A  , T50, A  )' )                              &
                            'external data file', 'namelist input'
 
            WRITE (*, '(T8,A,T30, I10, T50, I10)' )                              &
                            'ie_ext', ize   , ie_ext
            WRITE (*, '(T8,A,T30, I10, T50, I10)' )                              &
                            'je_ext', jze   , je_ext

            WRITE (*, '(T8,A,T30,F10.5,T50,F10.5)')                              &
                            'lower left lon', zlonll, startlon_tot
            WRITE (*, '(T8,A,T30,F10.5,T50,F10.5)')                              &
                            'lower left lat', zlatll, startlat_tot

            WRITE (*, '(T8,A,T30,F10.5,T50,F10.5)')                              &
                            'upper right lon', zlonur, endlon_tot
            WRITE (*, '(T8,A,T30,F10.5,T50,F10.5)')                              &
                            'upper right lat', zlatur, endlat_tot

            WRITE (*, '(T8,A,T30,F10.5,T50,F10.5)')                              &
                              'delta lon', zdlon, dlon
            WRITE (*, '(T8,A,T30,F10.5,T50,F10.5)')                              &
                          'delta lat', zdlat, dlat

            WRITE (*, '(T8,A,T30,F10.5,T50,F10.5)')                              &
                            'rot. south pole (lat)', zpollat_grib, zsouthp_lat
            WRITE (*, '(T8,A,T30,F10.5,T50,F10.5)')                              &
                            'rot. south pole (lon)', zpollon_grib, zsouthp_lon
          ENDIF

          ierror = 4
          yerror = 'wrong grid for LM external parameters'
          RETURN
          EXIT
        ENDIF
      ENDDO
    ENDIF

    !--------------------------------------------------------------------------
    ! 3.2.6: Cut out the proper domain
    !--------------------------------------------------------------------------

    ! All PEs set lrequired to .FALSE. again
    lrequired = .FALSE.

    ! Now it is going on only in the PEs with a record to process
    IF (.NOT. leof) THEN

      ! Cut out the proper global LM domain from ds_grib_double
      izshift = NINT ( (startlon_tot - dlon - zstartlon_ext)/dlon )
      jzshift = NINT ( (startlat_tot - dlat - zstartlat_ext)/dlat )
      DO j = 1, je2lm_tot
        DO i = 1, ie2lm_tot
          ij = i + izshift + (j - 1 + jzshift)*ie_ext
          lm_field(i,j) = ds_grib_single(ij)
        ENDDO
      ENDDO

      IF (ylmext_form_read == 'grb1') THEN
        ! Determine whether this record is required or not by comparing
        ! with ylistextpar (set lrequired) and put the data to the structure 
        ! procarray.

        ! This search gets more complicated as new external parameters are
        ! added, especially 3D fields with different levtyp or even with
        ! levtyp=1 but 12 monthly values. 
        ! But still, looking for ee and tabtyp should be sufficient
        search_grb1: DO izloc = 1, nvar_lm
          IF ( (var_lm(izloc)%ee     == ipds(7)) .AND. &
               (var_lm(izloc)%tabtyp == ipds(2)) )        THEN
            ! the field is in the LM variable table, check whether
            ! it is in the list of the external parameters
            DO j = 1, numlistextpar
              IF (var_lm(izloc)%name == ylistextpar(j)) THEN
                lrequired = .TRUE.
                EXIT search_grb1
              ENDIF
            ENDDO
          ENDIF
        ENDDO search_grb1
      ELSEIF (ylmext_form_read == 'apix') THEN
        ! Determine whether this record is required or not by comparing
        ! with ylistextpar (set lrequired) and put the data to the structure 
        ! procarray.
        search_apix: DO izloc = 1, nvar_lm
          IF ( TRIM(var_lm(izloc)%name) == TRIM(yshortname) ) THEN
            ! the field is in the LM variable table, check whether
            ! it is in the list of the external parameters
            DO j = 1, numlistextpar
              IF (var_lm(izloc)%name == ylistextpar(j)) THEN
                lrequired = .TRUE.
                EXIT search_apix
              ENDIF
            ENDDO
          ENDIF
        ENDDO search_apix
      ELSEIF (ylmext_form_read == 'ncdf') THEN
        search_ncdf: DO izloc = 1, nvar_lm
          IF (var_lm(izloc)%name == ylistextpar(myzvar)) THEN
            lrequired = .TRUE.
            EXIT search_ncdf
          ENDIF
        ENDDO search_ncdf
      ENDIF

      IF (lrequired) THEN

        ! Now get the rank and the level (or month) for 3D fields
        izrank = -1
        izlev  = -1
        IF (ylmext_form_read == 'grb1') THEN

          SELECT CASE (ipds(8))
          CASE (1)
            ! For 12-monthly values
            IF (  (var_lm(izloc)%name == 'NDVI_MRAT ') .OR.  &
                  (var_lm(izloc)%name == 'AER_SO412 ') .OR.  &
                  (var_lm(izloc)%name == 'AER_DUST12') .OR.  &
                  (var_lm(izloc)%name == 'AER_ORG12 ') .OR.  &
                  (var_lm(izloc)%name == 'AER_BC12  ') .OR.  &
                  (var_lm(izloc)%name == 'AER_SS12  ') .OR.  &
                  (var_lm(izloc)%name == 'ALB_DIF12 ') )        THEN
              izrank = 3
              izlev = ipds(12)    ! information about the month
            ELSE
              izrank = 2
              izlev = ipds(9)
            ENDIF
          CASE (109)
            izrank = 3
            izlev = ipds(10)
          CASE (110)
            izrank = 3
            izlev = ipds(9)
          END SELECT

#ifdef GRIBAPI
        ELSEIF (ylmext_form_read == 'apix') THEN

          CALL grib_get (igribid, 'typeOfLevel',      ytypeoflevel,  ireturn)
          CALL grib_get (igribid, 'level',            ilevel,        ireturn)
          CALL grib_get (igribid, 'topLevel',         itoplevel,     ireturn)
          CALL grib_get (igribid, 'bottomLevel',      ibottomlevel,  ireturn)

          SELECT CASE (ytypeoflevel)
          CASE ('surface')
            ! For 12-monthly values
            IF (  (var_lm(izloc)%name == 'NDVI_MRAT ') .OR.  &
                  (var_lm(izloc)%name == 'AER_SO412 ') .OR.  &
                  (var_lm(izloc)%name == 'AER_DUST12') .OR.  &
                  (var_lm(izloc)%name == 'AER_ORG12 ') .OR.  &
                  (var_lm(izloc)%name == 'AER_BC12  ') .OR.  &
                  (var_lm(izloc)%name == 'AER_SS12  ') .OR.  &
                  (var_lm(izloc)%name == 'ALB_DIF12 ') )        THEN
              CALL grib_get (igribid, 'dataDate',     ydate,  ireturn)
              READ (ydate(5:6), '(I2.2)') imonth
              izrank = 3
              izlev  = imonth      ! information about the month
            ELSE
              izrank = 2
              izlev  = ibottomlevel
            ENDIF
          CASE ('hybrid')
            izrank = 3
            izlev = itoplevel
          CASE ('hybridLayer')
            izrank = 3
            izlev = ibottomlevel
          END SELECT
#endif
        ELSEIF (ylmext_form_read == 'ncdf') THEN

          izrank = var_lm(izloc)%rank
          izlev  = myzlev

        ENDIF

        DO izproc = 0,num_compute-1
          ize = isubpos(izproc,3) - isubpos(izproc,1) + 1 + 2*nboundlines
          jze = isubpos(izproc,4) - isubpos(izproc,2) + 1 + 2*nboundlines
          DO j = 1, jze
            DO i = 1, ize
              ij = (j-1)*ize + i
              zprocarray(ij,izproc+1) =                                 &
                         lm_field(isubpos(izproc,1)-nboundlines-1+i,    &
                                  isubpos(izproc,2)-nboundlines-1+j)
            ENDDO
          ENDDO
        ENDDO

      ENDIF ! lrequired

    ENDIF   ! leof
 
  !----------------------------------------------------------------------------
  ! 3.3: Check the record
  !----------------------------------------------------------------------------

    ! print the headline for this grib file
    IF (lchkin) THEN
      IF ((linit) .AND. (my_cart_id == 0) ) THEN
        WRITE (nuchkdat,'(A)') 'Check the external LM parameters: '
        WRITE (nuchkdat,'(A,A)')                                           &
              '    File:   ',yname(1:LEN_TRIM(yname))
        WRITE (nuchkdat,'(A,I5,A,I5)')                                     &
              '    ie2lm_tot =',ie2lm_tot,'   je2lm_tot =',je2lm_tot
        WRITE (nuchkdat,'(A)') '    '
        WRITE (nuchkdat,'(A,A)')                                           &
              '   var        ee  lev        min    ',                      &
              'imin jmin               max    imax jmax              mean'
        WRITE (nuchkdat,'(A)')  '  '
        linit = .FALSE.
      ENDIF

      ! print the maximum, minimum and meanvalues of each record
      IF (.NOT. lrequired) THEN
        izloc = 1
      ENDIF
      CALL check_record (lm_field, 1, ie2lm_tot  , 1, je2lm_tot  , 1, 1,   &
                                   2, ie2lm_tot-1, 2, je2lm_tot-1, 1, 1,   &
                         zundef,    var_lm(izloc)%name, var_lm(izloc)%ee,  &
                         1     , lrequired, nuchkdat, num_compute,         &
                         icomm_cart, my_cart_id, yerrmsg, izerror)
    ENDIF

  !----------------------------------------------------------------------------
  ! 3.4: Distribute record to all PEs and put values to memory
  !----------------------------------------------------------------------------

    distribute_loop: DO izproc = 0, num_compute-1
      ! Check and distribute the next action
      IF (izproc == my_cart_id) THEN
        IF ( (.NOT. leof) .AND. (lrequired) ) THEN
          iz_info(1) = 0  ! This processor has data
          iz_info(2) = izloc
          iz_info(3) = nzid_hsurf
          iz_info(4) = nzid_fis
          iz_info(5) = nzid_frland
          iz_info(6) = izrank
          iz_info(7) = izlev
        ELSE
          IF (leof) THEN
            iz_info(1) = -1 ! No data because end of file reached
            iz_info(2) =  0
            iz_info(3) = -1
            iz_info(4) = -1
            iz_info(5) = -1
            iz_info(6) = -1
            iz_info(7) = -1
          ELSE
            iz_info(1) = -2 ! No data because data is not required
            iz_info(2) =  0
            iz_info(3) = -1
            iz_info(4) = -1
            iz_info(5) = -1
            iz_info(6) = -1
            iz_info(7) = -1
          ENDIF
        ENDIF
      ENDIF
    
      IF (num_compute > 1) THEN
        CALL distribute_values (iz_info, 7, izproc, imp_integers,      &
                                icomm_cart, izerror)
      ENDIF
      IF ( (idbg_level > 15) .AND. (my_cart_id == izproc) ) THEN
        PRINT *, '       Distribute information: ', iz_info, lzgothsurf, lzgotfis, lzgotfrland
      ENDIF

      IF (iz_info(1) == -1) EXIT  read_loop          ! all records are done
      IF (iz_info(1) == -2) CYCLE distribute_loop    ! record not required

      IF ((iz_info(3) == -1) .AND. (iz_info(4) == -1)) THEN
        ! Distribute the record to the other PEs. 
        CALL scatter_values (zprocarray, zsubarray, ie2lm_max*je2lm_max,    &
                             num_compute, imp_grib, izproc, icomm_cart,     &
                             yerrmsg, izerror)
        IF (izerror /= 0) THEN
          ierror = 5
          yerror = 'Error in scatter_values'
          RETURN
        ENDIF

        ! Put values into corresponding variables and scale the field
        IF    ( (ylmext_form_read == 'grb1') .OR. (ylmext_form_read == 'apix') ) THEN
          zbias     = var_lm(iz_info(2))%bias
          zfactor   = var_lm(iz_info(2))%factor
          IF (var_lm(iz_info(2))%rank == 2) THEN
            DO j = 1, je2lm
              DO i = 1, ie2lm
                ij = (j-1)*ie2lm + i
                var_lm(iz_info(2))%p2(i,j) =                                    &
                              REAL(zsubarray(ij), ireals) / zfactor - zbias
              ENDDO
            ENDDO
          ELSEIF (var_lm(iz_info(2))%rank == 3) THEN
            DO j = 1, je2lm
              DO i = 1, ie2lm
                ij = (j-1)*ie2lm + i
                var_lm(iz_info(2))%p3(i,j,iz_info(7)) =                         &
                              REAL(zsubarray(ij), ireals) / zfactor - zbias
              ENDDO
            ENDDO
          ENDIF
        ELSEIF (ylmext_form_read == 'ncdf') THEN
          IF (var_lm(iz_info(2))%rank == 2) THEN
            ! no scaling for NetCDF; here in INT2LM also not for Z0: values for
            ! this variable are just interpolated and given through to the output
            ! and are also not scaled there.
            DO j = 1, je2lm
              DO i = 1, ie2lm
                ij = (j-1)*ie2lm + i
                var_lm(iz_info(2))%p2(i,j) = REAL(zsubarray(ij), ireals)
              ENDDO
            ENDDO
          ELSEIF (var_lm(iz_info(2))%rank == 3) THEN
            DO j = 1, je2lm
              DO i = 1, ie2lm
                ij = (j-1)*ie2lm + i
                var_lm(iz_info(2))%p3(i,j,iz_info(7)) = REAL(zsubarray(ij), ireals)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
   
        ! Set corresponding logical flag to .TRUE.
        SELECT CASE (var_lm(iz_info(2))%name)
        CASE ('FIS       ')
          lfis__lm = .TRUE.
        CASE ('HSURF     ')
          lhsur_lm = .TRUE.
        CASE ('FR_LAND   ')
          lfrla_lm = .TRUE.
        CASE ('SOILTYP   ')
          lsoty_lm = .TRUE.
        CASE ('Z0        ')
          lz0___lm = .TRUE.
        CASE ('Z012      ')
          lzz0   (iz_info(7)) = .TRUE.
        CASE ('PLCOV_MX  ')
          lplmx_lm = .TRUE.
        CASE ('PLCOV_MN  ')
          lplmn_lm = .TRUE.
        CASE ('PLCOV12   ')
          lzplcov(iz_info(7)) = .TRUE.
        CASE ('LAI_MX    ')
          laimx_lm = .TRUE.
        CASE ('LAI_MN    ')
          laimn_lm = .TRUE.
        CASE ('LAI12     ')
          lzlai  (iz_info(7)) = .TRUE.
        CASE ('ROOTDP    ')
          lroot_lm = .TRUE.
        CASE ('FOR_E     ')
          lfore_lm = .TRUE.
        CASE ('FOR_D     ')
          lford_lm = .TRUE.
        CASE ('URBAN     ')
          lurba_lm = .TRUE.
        CASE ('EMIS_RAD  ')
          lemis_lm = .TRUE.
        CASE ('ALB_DRY   ')
          laldr_lm = .TRUE.
        CASE ('ALB_SAT   ')
          lalsa_lm = .TRUE.
        CASE ('RSMIN     ')
          lprsm_lm = .TRUE.
        CASE ('SSO_STDH  ')
          lstdh_lm = .TRUE.
        CASE ('SSO_GAMMA ')
          lgamm_lm = .TRUE.
        CASE ('SSO_THETA ')
          lthet_lm = .TRUE.
        CASE ('SSO_SIGMA ')
          lsigm_lm = .TRUE.
        CASE ('SKYVIEW   ')
          lskyv_lm = .TRUE.
        CASE ('SLO_ANG   ')
          lsang_lm = .TRUE.
        CASE ('SLO_ASP   ')
          lsasp_lm = .TRUE.
        CASE ('HORIZON   ')
          lhori_lm = .TRUE.
        CASE ('FR_LAKE   ')
          lflak_lm = .TRUE.
        CASE ('DEPTH_LK  ')
          ldept_lm= .TRUE.
        CASE ('SALT_LK   ')
          lsalt_lm= .TRUE.
        CASE ('NDVI_MRAT ')
          lzndvir(iz_info(7)) = .TRUE.
        CASE ('AER_SO412 ')
          lzaersu(iz_info(7)) = .TRUE.
        CASE ('AER_DUST12')
          lzaerdu(iz_info(7)) = .TRUE.
        CASE ('AER_ORG12 ')
          lzaeror(iz_info(7)) = .TRUE.
        CASE ('AER_BC12  ')
          lzaerbc(iz_info(7)) = .TRUE.
        CASE ('AER_SS12  ')
          lzaerss(iz_info(7)) = .TRUE.
        CASE ('ALB_DIF12 ')
          lzalb(iz_info(7)) = .TRUE.
        CASE ('T_CL      ')
          ltcl__lm  = .TRUE.
        END SELECT

        IF (lfilter_oro .AND. (iz_info(5) /= -1)) THEN
          ! distribute values for frland to all other PEs
          IF (my_cart_id == nzid_frland) THEN
            ds_grib_double = RESHAPE (frland_filt, (/ie_ext*je_ext/))
          ENDIF
          IF (num_compute > 1) THEN
            CALL distribute_values  (ds_grib_double, ie_ext*je_ext, iz_info(5),   &
                                     imp_reals, icomm_cart, izerror)
          ENDIF
          zfactor = var_lm(iz_info(2))%factor
          zbias   = var_lm(iz_info(2))%bias  
          DO j = 1, je_ext
            DO i = 1, ie_ext
              ij = (j-1)*ie_ext + i
              frland_filt(i,j) = ds_grib_double(ij) / zfactor - zbias
            ENDDO
          ENDDO
          lzdisfrland = .TRUE.
          IF ( (idbg_level > 15) .AND. (my_cart_id == iz_info(5)) ) THEN
            PRINT *, '       Distribute values for frland_filt ', lzdisfrland
          ENDIF
        ENDIF

      ELSEIF (lfilter_oro .AND. (iz_info(3) /= -1)) THEN
        ! distribute values for hsurf to all other PEs
        IF (my_cart_id == nzid_hsurf) THEN
          ds_grib_double = RESHAPE (field_filt, (/ie_ext*je_ext/))
        ENDIF
        IF (num_compute > 1) THEN
          CALL distribute_values  (ds_grib_double, ie_ext*je_ext, iz_info(3),     &
                                   imp_reals, icomm_cart, izerror)
        ENDIF
        zfactor = var_lm(iz_info(2))%factor
        zbias   = var_lm(iz_info(2))%bias  
        DO j = 1, je_ext
          DO i = 1, ie_ext
            ij = (j-1)*ie_ext + i
            field_filt(i,j) = ds_grib_double(ij) / zfactor - zbias
          ENDDO
        ENDDO
        lzdishsurf = .TRUE.
        lzdisfis   = .FALSE.   
          IF ( (idbg_level > 15) .AND. (my_cart_id == iz_info(3)) ) THEN
            PRINT *, '       Distribute values for field_filt ', lzdishsurf, lzdisfis
          ENDIF
      ELSEIF (lfilter_oro .AND. (iz_info(4) /= -1) .AND. (.NOT. lzdishsurf)) THEN
        ! distribute values for fis   to all other PEs
        IF (my_cart_id == nzid_fis) THEN
          ds_grib_double = RESHAPE (field_filt, (/ie_ext*je_ext/))
        ENDIF
        IF (num_compute > 1) THEN
          CALL distribute_values  (ds_grib_double, ie_ext*je_ext, iz_info(4),     &
                                   imp_reals, icomm_cart, izerror)
        ENDIF
        zfactor = var_lm(iz_info(2))%factor
        zbias   = var_lm(iz_info(2))%bias  
        DO j = 1, je_ext
          DO i = 1, ie_ext
            ij = (j-1)*ie_ext + i
            field_filt(i,j) = ds_grib_double(ij) / zfactor - zbias
          ENDDO
        ENDDO
        lzdisfis   = .TRUE.
          IF ( (idbg_level > 15) .AND. (my_cart_id == iz_info(4)) ) THEN
            PRINT *, '       Distribute values for field_filt ', lzdishsurf, lzdisfis
          ENDIF
      ENDIF

    ENDDO distribute_loop

    ! some values have to be re-initialized
    nzid_frland = -1
    nzid_hsurf  = -1
    nzid_fis    = -1
  ENDDO read_loop

  DEALLOCATE (iblock, ibmap, dsup, ds_grib_single)
  DEALLOCATE (ds_grib_double, lm_field, zprocarray, zsubarray)

  ! Set logical read flags for multi-dimensional variables
  IF (ALL(lzalb  (:))) lal12_lm = .TRUE.
  IF (ALL(lzndvir(:))) lndvi_lm = .TRUE.
  IF (ALL(lzaersu(:))) lsu12_lm = .TRUE.
  IF (ALL(lzaerdu(:))) ldu12_lm = .TRUE.
  IF (ALL(lzaeror(:))) lor12_lm = .TRUE.
  IF (ALL(lzaerbc(:))) lbc12_lm = .TRUE.
  IF (ALL(lzaerss(:))) lss12_lm = .TRUE.
  IF (ALL(lzplcov(:))) lpl12_lm = .TRUE.
  IF (ALL(lzlai  (:))) lai12_lm = .TRUE.
  IF (ALL(lzz0   (:))) lz012_lm = .TRUE.

!------------------------------------------------------------------------------
! Section 4: Treatment of orography
!------------------------------------------------------------------------------

 IF (lread) THEN

  !----------------------------------------------------------------------------
  ! Filling of valleys BEFORE filtering of orography
  !----------------------------------------------------------------------------

  IF ( ifill_valley == 1 ) CALL fill_valleys

  !----------------------------------------------------------------------------
  ! Filtering the orography
  !----------------------------------------------------------------------------

  IF (lfilter_oro) THEN

    ALLOCATE (field_ori   (ie_ext,je_ext), STAT = izstat)
    field_ori(:,:) = field_filt(:,:)

    SELECT CASE( ilow_pass_oro )

    CASE( 1 )

      DO nit = 1, numfilt_oro
        IF ( .NOT.lxso_first ) THEN
  
        !----------------------------------------------------------------------
        ! Standard horizontal filtering of orography (at all grid points)
        !----------------------------------------------------------------------
        ! ... BEFORE eXtra Smoothing of steep Orography
  
          !------------------------------------------------------------------------
          ! Section 4.1: filtering in x-direction 
          !------------------------------------------------------------------------
  
          ! Set the dimension
          ndim = ie_ext
    
          ! allocate the necessary fields for gaussian elimination
          ALLOCATE (xy_vec(ndim), ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim),&
                         a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim))
    
          ! Set the above variables for the gaussian elimination
          CALL set_gauss (a, b, c, d, e, f, ci, cj, ck, cl, cm, ndim, eps_filter)
    
          ! Apply the filter for every line
          DO j = 1, je_ext
            xy_vec(:) = field_filt(:,j)
            CALL low_pass_filter (xy_vec, a, b, c, d, e, f, ci, cj, ck, cl, cm,   &
                                  ndim, eps_filter)
            field_filt(:,j) = xy_vec(:)
          ENDDO
    
          ! Release memory
          DEALLOCATE (xy_vec, ci, cj, ck, cl, cm, a, b, c, d, e, f)
    
          !------------------------------------------------------------------------
          ! Section 4.2: filtering in y-direction 
          !------------------------------------------------------------------------
    
          ! Set the dimension
          ndim = je_ext
    
          ! allocate the necessary fields for gaussian elimination
          ALLOCATE (xy_vec(ndim), ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim),&
                         a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim))
    
          ! Set the above variables for the gaussian elimination
          CALL set_gauss (a, b, c, d, e, f, ci, cj, ck, cl, cm, ndim, eps_filter)
    
          ! Apply the filter for every column
          DO i = 1, ie_ext
            xy_vec(:) = field_filt(i,:)
            CALL low_pass_filter (xy_vec, a, b, c, d, e, f, ci, cj, ck, cl, cm,  &
                                  ndim, eps_filter)
            field_filt(i,:) = xy_vec(:)
          ENDDO
    
          ! Release memory
          DEALLOCATE (xy_vec, ci, cj, ck, cl, cm, a, b, c, d, e, f)
  
        ENDIF ! .NOT.lxso_first
  
        !------------------------------------------------------------------------
        ! eXtra Smoothing of steep Orography (only at masked grid points)
        !------------------------------------------------------------------------
  
        ! take into account steepness of orography in diagonal direction
        ldhsurf_xy = .TRUE.
  
        IF ( ilow_pass_xso >= ilow_pass_oro ) THEN
          IF ( rxso_mask > 0.0_ireals ) THEN
  
            IF ( lzdisfis ) rxso_mask = rxso_mask * g
  
            ! set width of the stencil for the horizontal filter
            SELECT CASE( ilow_pass_xso )
            CASE( 3, 4, 6 )
              hfwidth = 4
            CASE( 5, 8 )
              hfwidth = 6
            END SELECT
            hfw_m_nb = hfwidth - nboundlines
            ie_ext_hf = ie_ext + 2*hfw_m_nb
            je_ext_hf = je_ext + 2*hfw_m_nb
  
            ALLOCATE( hfx_mask(ie_ext_hf,je_ext_hf),  &
                      hfy_mask(ie_ext_hf,je_ext_hf),  &
                      STAT = istata )
  
            DO n = 1, numfilt_xso
  
              hfx_mask(:,:) = .FALSE.
              hfy_mask(:,:) = .FALSE.
              number_mask_points = 0
  
              ! set mask for extra smoothing
              DO j = 2, je_ext-1
                DO i = 2, ie_ext-1
  
                  zdh_x1 = ABS( field_filt(i-1,j) - field_filt(i,j) )
                  zdh_x2 = ABS( field_filt(i+1,j) - field_filt(i,j) )
                  zdh_max = MAX( zdh_x1, zdh_x2 )
                  IF ( zdh_max > rxso_mask ) THEN
                    hfx_mask(i+hfw_m_nb,j+hfw_m_nb) = .TRUE.
                    number_mask_points = number_mask_points + 1
                  ENDIF
  
                  zdh_y1 = ABS( field_filt(i,j-1) - field_filt(i,j) )
                  zdh_y2 = ABS( field_filt(i,j+1) - field_filt(i,j) )
                  zdh_max = MAX( zdh_y1, zdh_y2 )
                  IF ( zdh_max > rxso_mask ) THEN
                    hfy_mask(i+hfw_m_nb,j+hfw_m_nb) = .TRUE.
                    number_mask_points = number_mask_points + 1
                  ENDIF
  
                  IF ( ldhsurf_xy ) THEN
                    zdh_xy1 = ABS( field_filt(i-1,j-1) - field_filt(i,j) )
                    zdh_xy2 = ABS( field_filt(i+1,j+1) - field_filt(i,j) )
                    zdh_xy3 = ABS( field_filt(i-1,j+1) - field_filt(i,j) )
                    zdh_xy4 = ABS( field_filt(i+1,j-1) - field_filt(i,j) )
                    zdh_max = MAX( zdh_xy1, zdh_xy2, zdh_xy3, zdh_xy4 )
                    IF ( zdh_max > SQRT(2.0_ireals)*rxso_mask ) THEN
                      hfx_mask(i+hfw_m_nb,j+hfw_m_nb) = .TRUE.
                      hfy_mask(i+hfw_m_nb,j+hfw_m_nb) = .TRUE.
                      number_mask_points = number_mask_points + 2
                    ENDIF
                  ENDIF
  
                ENDDO
              ENDDO
  
              CALL hfilter_orography( ilow_pass_xso, .TRUE. )
              IF (idbg_level > 5) THEN
                PRINT *, " number_mask_points for ilow_pass_xso = ", number_mask_points
              ENDIF
  
            ENDDO
  
            DEALLOCATE( hfx_mask, hfy_mask, STAT = istatd )
  
          ENDIF ! ilow_pass_xso >= ilow_pass_oro
        ENDIF ! rxso_mask > 0.0_ireals
  
        IF ( lxso_first ) THEN
        !----------------------------------------------------------------------
        ! Standard horizontal filtering of orography (at all grid points)
        !----------------------------------------------------------------------
        ! ... AFTER eXtra Smoothing of steep Orography
          !------------------------------------------------------------------------
          ! Section 4.1: filtering in x-direction 
          !------------------------------------------------------------------------
    
          ! Set the dimension
          ndim = ie_ext
    
          ! allocate the necessary fields for gaussian elimination
          ALLOCATE (xy_vec(ndim), ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim),&
                         a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim))
    
          ! Set the above variables for the gaussian elimination
          CALL set_gauss (a, b, c, d, e, f, ci, cj, ck, cl, cm, ndim, eps_filter)
    
          ! Apply the filter for every line
          DO j = 1, je_ext
            xy_vec(:) = field_filt(:,j)
            CALL low_pass_filter (xy_vec, a, b, c, d, e, f, ci, cj, ck, cl, cm,   &
                                  ndim, eps_filter)
            field_filt(:,j) = xy_vec(:)
          ENDDO
    
          ! Release memory
          DEALLOCATE (xy_vec, ci, cj, ck, cl, cm, a, b, c, d, e, f)
    
          !------------------------------------------------------------------------
          ! Section 4.2: filtering in y-direction 
          !------------------------------------------------------------------------
    
          ! Set the dimension
          ndim = je_ext
    
          ! allocate the necessary fields for gaussian elimination
          ALLOCATE (xy_vec(ndim), ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim),&
                         a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim))
    
          ! Set the above variables for the gaussian elimination
          CALL set_gauss (a, b, c, d, e, f, ci, cj, ck, cl, cm, ndim, eps_filter)
    
          ! Apply the filter for every column
          DO i = 1, ie_ext
            xy_vec(:) = field_filt(i,:)
            CALL low_pass_filter (xy_vec, a, b, c, d, e, f, ci, cj, ck, cl, cm,  &
                                  ndim, eps_filter)
            field_filt(i,:) = xy_vec(:)
          ENDDO
    
          ! Release memory
          DEALLOCATE (xy_vec, ci, cj, ck, cl, cm, a, b, c, d, e, f)
  
        ENDIF ! lxso_first
      ENDDO ! numfilt_oro

    CASE( 3:6, 8 )

      !------------------------------------------------------------------------
      ! Section 4.1 + 4.2: Horizontal filtering of the orography
      !------------------------------------------------------------------------

      IF ( .NOT.lxso_first ) THEN

        !----------------------------------------------------------------------
        ! Standard horizontal filtering of orography (at all grid points)
        !----------------------------------------------------------------------
        ! ... BEFORE eXtra Smoothing of steep Orography

        ! set width of the stencil for the horizontal filter
        SELECT CASE( ilow_pass_oro )
        CASE( 3, 4, 6 )
          hfwidth = 4
        CASE( 5, 8 )
          hfwidth = 6
        END SELECT
        hfw_m_nb = hfwidth - nboundlines
        ie_ext_hf = ie_ext + 2*hfw_m_nb
        je_ext_hf = je_ext + 2*hfw_m_nb

        DO n = 1, numfilt_oro
          CALL hfilter_orography( ilow_pass_oro, .FALSE. )
        ENDDO

      ENDIF

      !------------------------------------------------------------------------
      ! eXtra Smoothing of steep Orography (only at masked grid points)
      !------------------------------------------------------------------------

      ! take into account steepness of orography in diagonal direction
      ldhsurf_xy = .TRUE.

      IF ( ilow_pass_xso >= ilow_pass_oro ) THEN
        IF ( rxso_mask > 0.0_ireals ) THEN

          IF ( lzdisfis ) rxso_mask = rxso_mask * g

          ! set width of the stencil for the horizontal filter
          SELECT CASE( ilow_pass_xso )
          CASE( 3, 4, 6 )
            hfwidth = 4
          CASE( 5, 8 )
            hfwidth = 6
          END SELECT
          hfw_m_nb = hfwidth - nboundlines
          ie_ext_hf = ie_ext + 2*hfw_m_nb
          je_ext_hf = je_ext + 2*hfw_m_nb

          ALLOCATE( hfx_mask(ie_ext_hf,je_ext_hf),  &
                    hfy_mask(ie_ext_hf,je_ext_hf),  &
                    STAT = istata )

          DO n = 1, numfilt_xso

            hfx_mask(:,:) = .FALSE.
            hfy_mask(:,:) = .FALSE.
            number_mask_points = 0

            ! set mask for extra smoothing
            DO j = 2, je_ext-1
              DO i = 2, ie_ext-1

                zdh_x1 = ABS( field_filt(i-1,j) - field_filt(i,j) )
                zdh_x2 = ABS( field_filt(i+1,j) - field_filt(i,j) )
                zdh_max = MAX( zdh_x1, zdh_x2 )
                IF ( zdh_max > rxso_mask ) THEN
                  hfx_mask(i+hfw_m_nb,j+hfw_m_nb) = .TRUE.
                  number_mask_points = number_mask_points + 1
                ENDIF

                zdh_y1 = ABS( field_filt(i,j-1) - field_filt(i,j) )
                zdh_y2 = ABS( field_filt(i,j+1) - field_filt(i,j) )
                zdh_max = MAX( zdh_y1, zdh_y2 )
                IF ( zdh_max > rxso_mask ) THEN
                  hfy_mask(i+hfw_m_nb,j+hfw_m_nb) = .TRUE.
                  number_mask_points = number_mask_points + 1
                ENDIF

                IF ( ldhsurf_xy ) THEN
                  zdh_xy1 = ABS( field_filt(i-1,j-1) - field_filt(i,j) )
                  zdh_xy2 = ABS( field_filt(i+1,j+1) - field_filt(i,j) )
                  zdh_xy3 = ABS( field_filt(i-1,j+1) - field_filt(i,j) )
                  zdh_xy4 = ABS( field_filt(i+1,j-1) - field_filt(i,j) )
                  zdh_max = MAX( zdh_xy1, zdh_xy2, zdh_xy3, zdh_xy4 )
                  IF ( zdh_max > SQRT(2.0_ireals)*rxso_mask ) THEN
                    hfx_mask(i+hfw_m_nb,j+hfw_m_nb) = .TRUE.
                    hfy_mask(i+hfw_m_nb,j+hfw_m_nb) = .TRUE.
                    number_mask_points = number_mask_points + 2
                  ENDIF
                ENDIF

              ENDDO
            ENDDO

            CALL hfilter_orography( ilow_pass_xso, .TRUE. )
            IF (idbg_level > 5) THEN
              PRINT *, " number_mask_points for ilow_pass_xso = ", number_mask_points
            ENDIF
          ENDDO

          DEALLOCATE( hfx_mask, hfy_mask, STAT = istatd )

        ENDIF
      ENDIF

      IF ( lxso_first ) THEN

        !----------------------------------------------------------------------
        ! Standard horizontal filtering of orography (at all grid points)
        !----------------------------------------------------------------------
        ! ... AFTER eXtra Smoothing of steep Orography

        ! set width of the stencil for the horizontal filter
        SELECT CASE( ilow_pass_oro )
        CASE( 3, 4, 6 )
          hfwidth = 4
        CASE( 5, 8 )
          hfwidth = 6
        END SELECT
        hfw_m_nb = hfwidth - nboundlines
        ie_ext_hf = ie_ext + 2*hfw_m_nb
        je_ext_hf = je_ext + 2*hfw_m_nb

        DO n = 1, numfilt_oro
          CALL hfilter_orography( ilow_pass_oro, .FALSE. )
        ENDDO

      ENDIF

    END SELECT

    !--------------------------------------------------------------------------
    ! Section 4.3: correct the filtered field where necessary
    !--------------------------------------------------------------------------

    ! Correct the filtered orography by setting the height to the 
    ! original value if
    !   - it is a sea-point surrounded by at least 4 other sea-points
    !   - the sign of filtered orography does not coincide with the 
    !     original sign

    zhmax_sea = 1.0_ireals
    IF ( lzdisfis ) zhmax_sea = zhmax_sea * g

    DO j = 1, je_ext
      DO i = 1, ie_ext
        IF ((frland_filt(i,j) < 0.5) .AND. (field_ori(i,j) <= zhmax_sea)) THEN
          ile = MAX (1,i-1)
          iri = MIN (i+1,ie_ext)
          jlo = MAX (1,j-1)
          jup = MIN (j+1,je_ext)
          nfrl = 0
          IF (frland_filt(ile,jlo) < 0.5) nfrl = nfrl + 1
          IF (frland_filt(i  ,jlo) < 0.5) nfrl = nfrl + 1
          IF (frland_filt(iri,jlo) < 0.5) nfrl = nfrl + 1
          IF (frland_filt(ile,j  ) < 0.5) nfrl = nfrl + 1
          IF (frland_filt(iri,j  ) < 0.5) nfrl = nfrl + 1
          IF (frland_filt(ile,jup) < 0.5) nfrl = nfrl + 1
          IF (frland_filt(i  ,jup) < 0.5) nfrl = nfrl + 1
          IF (frland_filt(iri,jup) < 0.5) nfrl = nfrl + 1
          IF (nfrl >= 4) THEN
            field_filt(i,j) = field_ori(i,j)
          ENDIF
        ENDIF

        IF (field_filt(i,j)*field_ori(i,j) < 0.0) THEN
          field_filt(i,j) = field_ori(i,j)
        ENDIF
      ENDDO
    ENDDO

    !--------------------------------------------------------------------------
    ! Filling of valleys AFTER filtering of orography
    !--------------------------------------------------------------------------

    IF ( ifill_valley == 2 ) CALL fill_valleys

    !--------------------------------------------------------------------------
    ! Section 4.4: additional smoothing for the z-coordinate version of LM
    !--------------------------------------------------------------------------
 
    IF (l_topo_z) THEN
      IF (lzdisfis) THEN
        ! this has to be done on the hsurf-values
        field_filt(:,:) = field_filt(:,:) / g
      ENDIF

      ! Computation of the height above sealevel for the half and full levels
      zr_d      =   287.05_ireals
 
      IF ( vcoord%ivctype == 1 ) THEN  ! Pressure-based vertical coordinate on input
        ! vertical coordinate parameters
        ! (this is also done in referlm, but we need it here!!!)
        DO k = 1, ke1lm
          IF (vcoord%vert_coord(k) <= vcoord%vcflat) THEN
            zak(k) = vcoord%vert_coord(k)*refatm%p0sl
            zbk(k) = 0.0
          ELSE
            zak(k) = vcoord%vcflat*refatm%p0sl*                     &
                                (1.0 - vcoord%vert_coord(k))/(1.0 - vcoord%vcflat)
            zbk(k) =  (vcoord%vert_coord(k) - vcoord%vcflat)/(1.0 - vcoord%vcflat)
          ENDIF
        ENDDO

        IF (refatm%irefatm == 1) THEN
          DO  k = 1, kelm
            zpxx = zak(k) + zbk(k)*refatm%p0sl
            zhhlr(k) = (zr_d/g)*LOG(refatm%p0sl/zpxx) &
                              *( refatm%t0sl - 0.5_ireals*refatm%dt0lp*LOG(refatm%p0sl/zpxx) )
          ENDDO
        ELSE IF (refatm%irefatm == 2) THEN
          zt00 = refatm%t0sl - refatm%delta_t
          DO  k = 1, kelm
            zpxx = zak(k) + zbk(k)*refatm%p0sl
            zhhlr(k) = refatm%h_scal*log( (exp( -zr_d*zt00/(g*refatm%h_scal) * &
                       log(zpxx/refatm%p0sl))*(zt00 + refatm%delta_t)-refatm%delta_t)/zt00)
          ENDDO
        ENDIF
        zhhlr(ke1lm) = 0.0_ireals
      ELSE                      ! vcoord%ivctype = 2, 3 or 4
        DO k = 1, ke1lm
          zhhlr(k)  = vcoord%vert_coord(k)
        ENDDO
      ENDIF

      IF (my_cart_id == 0) THEN
        WRITE (noutput,'(A)') '   Vertical Coordinate parameters and heights'
        WRITE (noutput,'(A)') '        for the Z-coordinate version of LM'
        WRITE (noutput,'(A)') '        '
        WRITE (noutput,'(A)') '     AK            BK           VCOORD         HHLR' 
        DO k = 1, ke1lm
          WRITE (noutput,'(A,4F15.5)') '          ', zak(k), zbk(k), vcoord%vert_coord(k), zhhlr(k)
        ENDDO
      ENDIF

      ! Smooth the topography in y-direction
      DO k = kelm, 2, -1
        DO j = 3, je_ext-1
          DO i = 1, ie_ext
            IF ((field_filt(i,j+1) >= zhhlr(k)) .AND. (field_filt(i,j  ) <  zhhlr(k)) ) THEN
              IF ((field_filt(i,j-1) >= zhhlr(k)) .OR. (field_filt(i,j-2) >= zhhlr(k)) ) THEN
                zhh3 = MAX (field_filt(i,j-1), field_filt(i,j-2))
                zhh1 = MAX (zhh3, field_filt(i,j+1))
                zhh2 = MIN (zhh1, zhh3)
                field_filt(i,j  ) = zhh2
                field_filt(i,j-1) = zhh2
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ! Smooth the topography in x-direction
      DO k = kelm, 2, -1
        DO j = 1, je_ext
          DO i = 3, ie_ext-1
            IF ((field_filt(i+1,j) >= zhhlr(k)) .AND. (field_filt(i,j  ) <  zhhlr(k)) ) THEN
              IF ((field_filt(i-1,j) >= zhhlr(k)) .OR. (field_filt(i-2,j) >= zhhlr(k)) ) THEN
                zhh3 = MAX (field_filt(i-1,j), field_filt(i-2,j))
                zhh1 = MAX (zhh3, field_filt(i+1,j))
                zhh2 = MIN (zhh1, zhh3)
                field_filt(i  ,j) = zhh2
                field_filt(i-1,j) = zhh2
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF (lzdisfis) THEN
        ! put values back to geopotential values
        field_filt(:,:) = field_filt(:,:) * g
      ENDIF

    ENDIF

    !--------------------------------------------------------------------------
    ! Section 4.5: Every processor has to get its domain from field_filt
    !--------------------------------------------------------------------------
 
    ! Cut out the proper global LM domain from ds_grib_double
    ! allow for a safety-epsilon when checking the longitudes
    IF (startlon - zstartlon_ext < -1.0E-5_ireals) THEN
      izshift = NINT ( (startlon - zstartlon_ext + 360.0_ireals)/dlon )
    ELSE
      izshift = NINT ( (startlon - zstartlon_ext)/dlon )
    ENDIF
    jzshift = NINT ( (startlat - zstartlat_ext)/dlat )
    IF (lzdishsurf) THEN
      DO j = 1, je2lm
        DO i = 1, ie2lm
          hsurf_lm (i,j) = field_filt(i+izshift, j+jzshift)
          fis_lm   (i,j) = field_filt(i+izshift, j+jzshift) * g
        ENDDO
      ENDDO
      lhsur_lm = .TRUE.
      lfis__lm = .TRUE.
    ELSEIF (lzdisfis) THEN
      DO j = 1, je2lm
        DO i = 1, ie2lm
          fis_lm   (i,j) = field_filt(i+izshift, j+jzshift)
          hsurf_lm (i,j) = field_filt(i+izshift, j+jzshift) / g
        ENDDO
      ENDDO
      lhsur_lm = .TRUE.
      lfis__lm = .TRUE.
    ENDIF

    DEALLOCATE (field_filt, frland_filt, field_ori)
  ENDIF

 ENDIF ! lread

!------------------------------------------------------------------------------
! Section 5: Cleanup
!------------------------------------------------------------------------------
 
  IF (lread) THEN
    ! grib file has been opened, so it is closed now
    CALL close_file (nufile, ylmext_form_read, icomm_cart, my_cart_id,    &
                     num_compute, lasync_io, idbg_level, yerrmsg, izerror)
    IF (izerror /= 0) THEN
      ierror = 6
      yerror = 'Error in close_file'
      RETURN
    ENDIF

    IF (lchkin .AND. my_cart_id == 0) THEN
      WRITE (nuchkdat,'(A)')  '  '
      WRITE (nuchkdat,'(A)')  '  '
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================

SUBROUTINE fill_valleys

  ! Check if dh <= rfill_valley in a V-valley
  ! otherwise fill valley with maximum value of neighborhood
  IF ( rfill_valley > 0.0_ireals ) THEN

    IF ( lzdisfis ) rfill_valley = rfill_valley * g

    ! Allocate and set temporary field
    ALLOCATE( ff_tmp(ie_ext,je_ext), STAT = istata )
    ff_tmp(:,:) = field_filt(:,:)

    DO j = 2, je_ext-1
      DO i = 2, ie_ext-1

        zdh_x1 = ff_tmp(i-1,j) - ff_tmp(i,j)
        zdh_x2 = ff_tmp(i+1,j) - ff_tmp(i,j)
        zdh_y1 = ff_tmp(i,j-1) - ff_tmp(i,j)
        zdh_y2 = ff_tmp(i,j+1) - ff_tmp(i,j)

        IF ( zdh_x1 > 0.0_ireals .AND. zdh_x2 > 0.0_ireals .AND.      &
             zdh_y1 > 0.0_ireals .AND. zdh_y2 > 0.0_ireals ) THEN
          zdh_max = MAX( zdh_x1, zdh_x2, zdh_y1, zdh_y2 )

          IF ( zdh_max > rfill_valley ) THEN
            ! 1-3: BEFORE  / 4-6: AFTER filtering of orography
            SELECT CASE( ifill_valley )
            CASE( 1, 2 )
              ! MIN value      of 4-point neighborhood
              field_filt(i,j) = MIN( ff_tmp(i-1,j), ff_tmp(i+1,j),    &
                                     ff_tmp(i,j-1), ff_tmp(i,j+1) )
            END SELECT
          ENDIF

        ENDIF

      ENDDO
    ENDDO

    DEALLOCATE( ff_tmp, STAT = istatd )

  ENDIF

END SUBROUTINE fill_valleys

!==============================================================================

SUBROUTINE hfilter_orography( ncutoff, lhf_mask )

INTEGER (KIND=iintegers), INTENT(in) ::  &
  ncutoff

LOGICAL, INTENT(in) ::  &
  lhf_mask

  dummy_neigh(:) = -1_iintegers

  ALLOCATE( ff_tmp(ie_ext_hf,je_ext_hf), STAT = istata )

  CALL extend_field ( field_filt(:,:), ie_ext, je_ext,                &
                      ff_tmp(:,:), ie_ext_hf, je_ext_hf, 1,           &
                      hfw_m_nb, nboundlines, hfjstartpar, hfjendpar,  &
                      sendbuf, isendbuflen, imp_reals, icomm_cart,    &
                      dummy_neigh, num_compute, .FALSE., .FALSE., .FALSE. )

  IF ( lhf_mask ) THEN
    CALL horizontal_filtering( ff_tmp(:,:), ie_ext_hf, je_ext_hf, 1,  &
                               nboundlines, hfwidth, ncutoff,         &
                               dummy_neigh, hfx_mask, hfy_mask )
  ELSE
    CALL horizontal_filtering( ff_tmp(:,:), ie_ext_hf, je_ext_hf, 1,  &
                               nboundlines, hfwidth, ncutoff,         &
                               dummy_neigh )
  ENDIF

  DO j = 1, je_ext
    DO i = 1, ie_ext
      field_filt(i,j) = ff_tmp(i+hfw_m_nb,j+hfw_m_nb)
    ENDDO
  ENDDO

  DEALLOCATE( ff_tmp, STAT = istatd )

END SUBROUTINE hfilter_orography

!==============================================================================

END SUBROUTINE read_lm_ext

!==============================================================================
!+ Reads the external parameters for the coarse grid
!------------------------------------------------------------------------------

SUBROUTINE read_coarse_grid_ext                                              &
               (lfis__in, lfrla_in, lsoty_in, lz0___in, lplmx_in, lplmn_in,  &
                laimx_in, laimn_in, lroot_in,                                &
                lfis__lm, lfrla_lm, lsoty_lm, lz0___lm, lplmx_lm, lplmn_lm,  &
                laimx_lm, laimn_lm, lroot_lm, ierror, yerror)

!------------------------------------------------------------------------------
!
! Description:
!  Reads the fields of external parameters for the coarse grid model from a 
!  GRIB1 file, unpacks and saves them in the appropriate arrays.
!  The organization of this routine is such, that it works on parallel and
!  on sequential platforms.
!
! Method:
!  After creating the filename and opening the grib file, all records are read
!  in a read_loop. Every processor gets the same total record for de-gribing,
!  because for the interpolation the global fields are needed.
!
! Input files:
!  Grib-file with external parameters for GME. 
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
LOGICAL,     INTENT (INOUT)   ::                                 &
  ! logicals, to indicate whether the corresponding field has been read
  lfis__in, lfrla_in, lsoty_in, lz0___in, lplmx_in, lplmn_in,    &
  laimx_in, laimn_in, lroot_in,                                  &
  lfis__lm, lfrla_lm, lsoty_lm, lz0___lm, lplmx_lm, lplmn_lm,    &
  laimx_lm, laimn_lm, lroot_lm

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror       ! error status

CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yerror       ! error message

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER  (KIND=intgribf)   ::  &
  ierrf, i_bits,        & ! error code for grib routines
  izbyta, izdimpds, izdimgds, idummy

INTEGER  (KIND=intgribc)   ::  &
  nufilec,              & ! unit number of opened grib file
  maxlenc,              & ! Maximum length of grib record (in intgribc)
  ierrc,                & ! Error return code
  igriblenc               ! actual length of grib record read

INTEGER  (KIND=iintegers)  ::  &
  nufile,               & ! unit number of opened grib file
  numlistextpar,        & ! number of external fields
  izerror, izstat,      & ! status and error status variable
  izdebug,              & ! for verbosity of output
  izerrorgrid, nzbyte,  & ! error status variable for grid checking
  izloc,                & ! location of field in variable table
  iz_info(3),           & ! characteristics of read records
  iz_lfd,               & !
  maxlen, igriblen,     & ! maximal and actual length of grib record read
  k, i, j, ij, n, jscan,& !
  izshift, jzshift,     & !
  izvar_count,          & !
  ildsin, ilfdin,       & !
  imonth, ilev,         & ! which month of the ndvi monthly rates have been read
  inrpoints,            & !
  iz_extfile, jz_extfile  ! field in external file can be smaller

! Variables for checking the data
INTEGER (KIND=iintegers)   ::  &
  igmemin(3), igmemax(3),      &  ! location of min and max in the GME field
  i_inmin(2), i_inmax(2)          ! location of min and max in the coarse 
                                  ! grid field

REAL    (KIND=ireals)      ::  &
  zmin, zmax, zsum, zmean,     &  ! values of min, max and meanvalue
  zfactor, zbias

REAL (KIND=ireals), ALLOCATABLE         ::  &
  ! field in the used REAL format
  field_gme(:,:,:),            & ! for GME fields
  field_in (:,:)                 ! for regular grid input fields

REAL (KIND=irealgrib), ALLOCATABLE      ::  &
  ds_gme(:),                   & !
  field_gme_file(:,:,:),       & ! field read from the file
  field_in_file (:)              !

REAL (KIND=ireals), ALLOCATABLE         ::  &
  zhsurfs_in_tot (:,:,:)         ! field for SLEVE orography splitting

REAL    (KIND=irealgrib)    :: refstf, zundef

LOGICAL                    ::  &
  lzrequired, & ! indicates whether a record from the grib file is required
  lzeof,      & ! indicates the end of file
  lzvertcoor, & ! whether vertical coordinates from LM could be read
  lread,      & ! go on with reading records
  linit         ! for the first cycle of read_loop

INTEGER                    :: &
  ivar_id(7)          ! ID of the external variable

LOGICAL                    :: &
  lcheckin_ext(7), & ! list for checking which variables has been read
  lndvi_mr(12)       ! to check whether all 12 monthly ratios have been read

CHARACTER (LEN=10)         ::  &
  ylistextpar (7)    ! list of external parameters that should be read

CHARACTER (LEN=14)         ::  &
  yzfulldate ! needed for the check-routines; including minutes and seconds

CHARACTER (LEN=250)        ::  &
  yname,   & ! name of the grib file (is determined in create_filename)
  yzerror

CHARACTER (LEN= 25)        ::  &
  yroutine   ! name of this routine for error handling

CHARACTER (LEN= 30)        ::  &
  yshortname,   & ! short name of the record read (from grib_api)
  ytypeoflevel    ! typeOfLevel

CHARACTER (LEN= 10)        ::  &
  yzlocname  ! short name from the Grib tables

CHARACTER (LEN=200)        ::  &
  yerrmsg    ! error message for error handling

! variables for grib meta data
INTEGER(KIND=intgribf)     ::  &
  iee, itabtyp, ilevtyp, byte_size, byte_size_t, byte_size_p, ireturn, &
  nb_values, nb_values_t, igribid, nb_field_in, idisc, icatg, ipara,   &
  nzscan, iz_ni_gme, iz_ni2, iz_ni3, iz_nd

LOGICAL                    ::  &
  lzgetvertcoord

REAL (KIND=ireals), ALLOCATABLE :: pv(:)

INTEGER  (KIND=iintegers)  ::  &
  izvert

!
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! for the external data, minutes and seconds should be 0
  yzfulldate(1:14) = ydate_ini(1:10)//'0000'

  ! first initializations
  ierror      = 0_iintegers
  yerror      = '    '
  igrbednr    = 0 ! neither grib1 nor grib2, then it is netcdf

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  ! initialize level counter and logical flags
  lzeof       = .FALSE.
  lzrequired  = .FALSE.
  lzvertcoor  = .FALSE.
  linit       = .TRUE.
  ierrc       = 0_intgribc
  ierrf       = 0_intgribf
  izstat      = 0_iintegers
  izerror     = 0_iintegers
  izbyta      = 9_intgribf
  yroutine    = 'read_coarse_grid_ext'
  yzerror     = '          '

  ! to check monthly mean values for gme grid
  lndvi_mr(:) = .FALSE.

  ! initialize size of field in external parameter file
  IF     (lcm2lm) THEN
    iz_extfile = ie_in_tot - east_add_in  - west_add_in
    jz_extfile = je_in_tot - south_add_in - north_add_in
  ELSE
    iz_extfile = ie_in_tot
    jz_extfile = je_in_tot
  ENDIF

  ! Put together ylistextpar (needed for NetCDF)
  ! Scan the list of variables and check the corresponding flag dattyp
  numlistextpar   = 0
  DO izloc = 1, nvar_in
    IF (var_in(izloc)%dattyp(3:3) == 'E') THEN
       numlistextpar = numlistextpar + 1
       ylistextpar(numlistextpar) =  var_in(izloc)%name
    ENDIF
  ENDDO

  ! Check for some additional parameters eventually not read for the fine grid
  IF (.NOT. lz0___lm) THEN
    numlistextpar              = numlistextpar + 1
    ylistextpar(numlistextpar) = 'Z0        '
    var_in( 4)%dattyp(3:3)     = 'E'
  ENDIF

  IF (.NOT. (lplmn_lm .AND. lplmx_lm) ) THEN
    numlistextpar              = numlistextpar + 1
    ylistextpar(numlistextpar) = 'PLCOV     '
    var_in( 7)%dattyp(3:3)     = 'E'
  ENDIF

  IF (.NOT. lroot_lm) THEN
    numlistextpar              = numlistextpar + 1
    ylistextpar(numlistextpar) = 'ROOTDP    '
    var_in(13)%dattyp(3:3)     = 'E'
  ENDIF

  ! Grib 1 initializations
  ! initialize variable iwlength
  i_bits      = BIT_SIZE (ierrf)
  iwlength    = INT (i_bits, iintegers) / 8

  ! Set dimensions for grib variables and allocate iblock, ibmap, ds and dsup
  ! nzbyte was assumed to be 2 here: but with new models and grib_api, we just
  ! change it to 8
  ! To know this a priori, the pds of the first record has
  ! to be read and decoded to get nrbit. The 2000 are just a safety-add to
  ! take care of the definition sections that are also stored in iblock
  nzbyte = 8

  ildsin  = MAX ((ni_gme+1)*(ni_gme+1)*10, ie_in_tot * je_in_tot)
  lds     = INT (ildsin, intgribf)
  ilfdin  = ildsin  * nzbyte / iwlength + 2000
  lfd     = INT (ilfdin, intgribf)
  iz_lfd  = INT (lfd , iintegers)

  ! dimensions for grib 1 routines
  idims_in( 1)   = npds
  idims_in( 2)   = ngds
  idims_in( 3)   = nbms
  idims_in( 4)   = nbds
  idims_in( 5)   = nbitmap
  idims_in( 6)   = ndsup
  idims_in( 7)   = lds
  idims_in( 8)   = lfd
  idims_in(9:20) = 0

  IF     ((yinext_form_read == 'grb1') .OR. (yinext_form_read == 'apix')) THEN
    zundef = undefgrib
  ELSEIF (yinext_form_read == 'ncdf') THEN
    zundef = undefncdf
  ENDIF
  undef       = REAL (zundef, ireals)

  ! Allocate fields
  ALLOCATE (iblock(lfd), ibmap(nbitmap), STAT=izstat)
  ALLOCATE (dsup(ndsup), ds_grib_single(lds),   STAT=izstat)
  byte_size = lfd * iwlength

  IF (lgme2lm) THEN
    nb_field_in = (igg1e-igg1s+1)*(igg2e-igg2s+1)*11
    ALLOCATE (field_gme (igg1sm2:igg1ep2,igg2sm2:igg2ep2, 1:10 ),   &
              field_gme_file(igg1s:igg1e,igg2s:igg2e, 1:11 ), STAT=izstat)
    ALLOCATE (ds_gme(nb_field_in)                           , STAT=izstat)
  ELSEIF (llm2lm .OR. lec2lm .OR. lgfs2lm .OR. lgsm2lm .OR. lum2lm .OR. lhir2lm) THEN
    ALLOCATE (field_in     (ie_in_tot, je_in_tot),                  &
              field_in_file(ie_in_tot* je_in_tot), STAT=izstat)
    nb_field_in = ie_in_tot * je_in_tot
  ELSEIF (lcm2lm) THEN
    ALLOCATE (field_in     (ie_in_tot, je_in_tot),                  &
              field_in_file(iz_extfile*jz_extfile), STAT=izstat)
    nb_field_in = iz_extfile*jz_extfile
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Create file name and open the file
!------------------------------------------------------------------------------

  ! Construct the file name
  yname = TRIM(yinext_cat)//TRIM(yinext_lfn)

  ! All processors have to call the routine open_file. What the parallel
  ! program really does is determined in the routine.
  CALL open_file(nufile, yname, ymode_read, yinext_form_read, icomm_cart, &
                 my_cart_id, num_compute, lasync_io, idbg_level,          &
                 yerrmsg, izerror)
  nufilec = INT (nufile, intgribc)
  IF (izerror /= 0) THEN
    lread = .FALSE.
    ierror = 1
    yerror = 'Error in open_file'
    RETURN
  ELSE 
    lread = .TRUE.
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Read records from the file
!------------------------------------------------------------------------------
 
  ! for NetCDF, the metadata of the file has to be scanned first
  IF (yinext_form_read == 'ncdf') THEN
    CALL read_nc_gdefs_ext_in (nufile, ie_in_tot, je_in_tot, ke_in,          &
            startlon_in_tot, startlat_in_tot, endlon_in_tot, endlat_in_tot,  &
            dlon_in, dlat_in, pollon_in, pollat_in, polgam_in,               &
            east_add_in, west_add_in, south_add_in, north_add_in,            &
            icomm_cart, my_cart_id, num_compute, yerrmsg, izerror)
    IF (izerror /= 0) THEN
      ierror = 2
      yerror = 'Error in read_nc_gdefs_ext_in'
      RETURN
    ENDIF

    ! Check, whether all data necessary are present and write information to
    ! output
    CALL read_nc_vdefs_ext_in (nufile, var_in, nvar_in, ivar_id, pollon_in,  &
            pollat_in, numlistextpar, ylistextpar, lcheckin_ext,             &
            icomm_cart, my_cart_id, num_compute, yerrmsg, izerror)
    IF (izerror /= 0) THEN
      ierror = 3
      yerror = 'Error in read_nc_vdefs_ext_in'
      RETURN
    ENDIF

    IF (izdebug > 10) THEN
      PRINT *, '  Necessary external parameters for coarse grid:'
      DO n = 1, numlistextpar
        PRINT *, '     ', ylistextpar(n), ivar_id(n), lcheckin_ext(n)
      ENDDO
    ENDIF
  ELSE
    IF (izdebug > 10) THEN
      PRINT *, '  Necessary external parameters for coarse grid:'
      DO n = 1, numlistextpar
        PRINT *, '     ', ylistextpar(n)
      ENDDO
    ENDIF
  ENDIF

  ! set logical if we've already decoded vertical coordinate parameters
  lzgetvertcoord = .TRUE.

  ! for NetCDF: initialize counter
  izvar_count = 0

  read_loop: DO WHILE (lread)

    nzscan = -1 ! set scanning mode to undefined

  !----------------------------------------------------------------------------
  ! 3.1: Get a record
  !----------------------------------------------------------------------------

    ! endless loop for reading records: this loop is exited
    ! if a record which is needed is found

    IF (my_cart_id == 0) THEN

      IF     (yinext_form_read == 'grb1') THEN
#ifdef GRIBDWD
        ! set edition number to 1
        igrbednr = 1

        ! Read in one GRIB field
        IF (lasync_io .OR. (num_compute > 1) ) THEN
          maxlen  = INT (idims_in(8)*iwlength, iintegers)
          CALL mpe_io_read(nufile, iblock, igriblen, maxlen, izerror)
          IF (izerror /= 0) THEN
            ierror = 2
            yerror = 'Error in mpe_io_read'
            RETURN
          ENDIF
        ELSE
          maxlenc = INT (idims_in(8)*iwlength, intgribc)
          CALL cuegin (nufilec, maxlenc, iblock, igriblenc, ierrc)
          igriblen = INT (igriblenc, iintegers)
          IF ( ierrc /= 0 ) THEN
            ierror = 2
            yerror = 'Error in cuegin'
            PRINT *, 'length of data:  ', maxlenc, igriblenc, ierrc
            RETURN
          ENDIF
        ENDIF
        IF (idbg_level >= 10) THEN
          PRINT *, ' Read a record with griblen = ', igriblen
        ENDIF
#endif
      ELSEIF (yinext_form_read == 'apix') THEN
#ifdef GRIBAPI

        ! Read a message and already build the grib handler here
        CALL grib_new_from_file(nufile, igribid, ireturn)

        IF ( ireturn /= GRIB_SUCCESS ) THEN
          IF (ireturn == GRIB_END_OF_FILE) THEN
            !EOF is reached
            lzeof    = .TRUE.
            igriblen = 0
          ELSE
            ierror = 2
            yerror = 'Error in grib_new_from_file'
            RETURN
          ENDIF
        ELSE
          ! get edition number
          CALL grib_get (igribid, 'editionNumber', igrbednr, ireturn)
          IF ( ireturn /= GRIB_SUCCESS ) THEN
            ierror = 2
            yerror = 'Error in grib_get:  editionNumber   '
            RETURN
          ENDIF
          ! get total length of message
          CALL grib_get (igribid, 'totalLength', igriblen, ireturn)
          IF ( ireturn /= GRIB_SUCCESS ) THEN
            ierror = 2
            yerror = 'Error in grib_get:  totalLength  '
            RETURN
          ENDIF
        ENDIF
        IF (idbg_level >= 10) THEN
          PRINT *, ' Read a grib_api record with ID = ', igribid, ' and length = ', igriblen, igrbednr
        ENDIF
#endif
      ELSEIF (yinext_form_read == 'ncdf') THEN
#ifdef NETCDF
        ! set edition number to 0 and scanning mode to 64 for NetCDF
        igrbednr = 0
        nzscan    = 64
        izvar_count = izvar_count + 1

        ! Cycle, if external parameter does not exist in coarse grid data set
        ! if parameter is necessary, program will stop at the end of external_data
        IF (ivar_id(izvar_count) == -1)  THEN
         CYCLE
        ENDIF

        IF (izvar_count > numlistextpar) THEN
          ! EOF is reached
          igriblen = 0
        ELSE
          ierror = nf90_get_var (nufile, ivar_id(izvar_count), field_in_file, &
                          start=(/1, 1, 1/), count=(/iz_extfile, jz_extfile, 1/) )
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
          igriblen = iz_extfile * jz_extfile
        ENDIF
#endif
      ENDIF

  !----------------------------------------------------------------------------
  ! 3.2: Check whether this variable is required and get location in table
  !----------------------------------------------------------------------------

      ! preset izloc; if no record is found, it will remain 0
      izloc = 0

      IF (igriblen > 0) THEN ! otherwise end of file reached

        ! Proc 0 determines the variable and whether it is required
        lzrequired  = .FALSE.

        ! get some metadata to interpret the data
        IF     (yinext_form_read=='grb1') THEN
#ifdef GRIBDWD
          ! Get gds and pds before further processing
          izbyta = 9_intgribf
          CALL getpd1 (idwdednr, izbyta, undefgrib, idims_in(8), idims_in(1),&
                       iblock, ipds, izdimpds, ierrf)
          izbyta = izbyta + ipds(1)
          CALL getgd1 (idwdednr, izbyta, undefgrib, idims_in(8), idims_in(2),&
                       iblock, igds_in, izdimgds, ierrf)

          IF (lgme2lm) THEN
            IF (igds_in(4) == 0) THEN
              ! this is a variable on a regular GME grid which is not required
              ! cycle the do loop
              CYCLE read_loop
            ENDIF

            iz_ni_gme = igds_in(8)
            iz_ni2    = igds_in(5)
            iz_ni3    = igds_in(6)
            iz_nd     = igds_in(7)
          ENDIF

          itabtyp = ipds( 2)
          iee     = ipds( 7)
          ilevtyp = ipds( 8)
          nzscan  = igds_in(14)

          ! Scan through the variable list to get izloc and yshortname
          get_loc_grb1: DO izloc = 1, nvar_in
            IF ( ((var_in(izloc)%ee     == iee    ) .AND.                   &
                  (var_in(izloc)%tabtyp == itabtyp) .AND.                   &
                  (var_in(izloc)%levtyp == ilevtyp))    .AND.               &
                  (var_in(izloc)%dattyp(3:3) == 'E') ) THEN
              yshortname = TRIM(var_in(izloc)%name)
              lzrequired = .TRUE.
              EXIT get_loc_grb1
            ENDIF
          ENDDO get_loc_grb1

          ! unpack the data, if it is required
          IF (lzrequired) THEN
            IF (lgme2lm) THEN
!             CALL grbin1 (idwdednr, undefgrib, ndims, idims_in, iblock, ibmap, ipds, &
!                          igds_in, ibms, ibds, dsup, field_gme_file, ierrf)
              CALL grbin1 (idwdednr, undefgrib, ndims, idims_in, iblock, ibmap, ipds, &
                           igds_in, ibms, ibds, dsup, ds_gme,         ierrf)
              IF (itype_ndvi == 1) THEN
                imonth = ipds(12)
              ENDIF
            ELSEIF (llm2lm .OR. lec2lm .OR. lum2lm) THEN
              CALL grbin1 (idwdednr, undefgrib, ndims, idims_in, iblock, ibmap, ipds, &
                           igds_in, ibms, ibds, dsup, field_in_file, ierrf)
            ENDIF
            IF ( ierrf /= 0 ) THEN
              ierror = 3
              yerror = 'Error in grbin1'
              RETURN
            ENDIF
          ENDIF

          IF (idbg_level >= 10) THEN
            PRINT *, ' dwdlib Record ', itabtyp, iee, ilevtyp, ' is needed:  ', lzrequired, izloc
          ENDIF
#endif
        ELSEIF (yinext_form_read=='apix') THEN
#ifdef GRIBAPI

          IF ((igrbednr == 1) .AND. (yinput_model == 'IFS')) THEN
            ! The COSMO shortnames are not available for the IFS Grib1 coding, therefore
            ! we need the Grib1 style search here: get itabtyp, iee, ilevtyp, nzscan
            CALL grib_get (igribid, 'table2Version',          itabtyp,       ireturn)
            CALL grib_get (igribid, 'indicatorOfParameter',   iee,           ireturn)
            CALL grib_get (igribid, 'indicatorOfTypeOfLevel', ilevtyp,       ireturn)
            CALL grib_get (igribid, 'scanningMode',           nzscan,        ireturn)

            ! Scan through the variable list to get izloc and yshortname
            get_loc_grb_ifs: DO izloc = 1, nvar_in
              IF ( ((var_in(izloc)%ee     == iee    ) .AND.                 &
                    (var_in(izloc)%tabtyp == itabtyp) .AND.                 &
                    (var_in(izloc)%levtyp == ilevtyp))    .AND.             &
                    (var_in(izloc)%dattyp(3:3) == 'E') ) THEN
                yshortname = TRIM(var_in(izloc)%name)
                lzrequired = .TRUE.
                EXIT get_loc_grb_ifs
              ENDIF
            ENDDO get_loc_grb_ifs
          ELSE
            ! get shortname and some other meta data
            CALL grib_get (igribid, 'shortName',   yshortname,   ireturn)
            CALL grib_get (igribid, 'typeOfLevel', ytypeoflevel, ireturn)
            IF (lgme2lm) THEN
              CALL grib_get (igribid, 'Ni', iz_ni_gme, ireturn)
              CALL grib_get (igribid, 'n2', iz_ni2   , ireturn)
              CALL grib_get (igribid, 'n3', iz_ni3   , ireturn)
              CALL grib_get (igribid, 'nd', iz_nd    , ireturn)
              IF ( ireturn /= GRIB_SUCCESS ) THEN
                PRINT *, ' *** ERROR in grib_get: GME grid values:  ', ireturn
                ierror = 4
                yerror = ' *** ERROR in grib_get: GME grid values:  '
                RETURN
              ENDIF
            ELSEIF (lhir2lm) THEN
              ! to care for some HIRLM specialities
              IF (yshortname == 'FI') yshortname = 'FIS'
              CALL grib_get (igribid, 'scanningMode',   nzscan, ireturn)
            ELSE
              CALL grib_get (igribid, 'scanningMode',   nzscan, ireturn)
              IF ( ireturn /= GRIB_SUCCESS ) THEN
                PRINT *, ' *** ERROR in grib_get: scanningMode:  ', ireturn
              ENDIF
            ENDIF

            IF (igrbednr == 2) THEN
              ! there must not be unknown names in Grib2!
              IF (TRIM(yshortname) == 'GH') THEN
                CALL grib_get (igribid, 'typeOfLevel', ytypeoflevel,  ireturn)
                IF ( ireturn /= GRIB_SUCCESS ) THEN
                  PRINT *, ' *** ERROR in grib_get: typeOfLevel: ', ireturn
                ENDIF
                IF (TRIM(ytypeOfLevel) == 'surface') THEN
                ! it is the height in gpm. As surface height it is taken as meters
                  yshortname = 'HSURF'
                END IF
              ELSEIF ( (TRIM(yshortname)== 'FI') .AND. (yinput_model == 'IFS')) THEN
                CALL grib_get (igribid, 'typeOfLevel', ytypeoflevel,  ireturn)
                IF ( ireturn /= GRIB_SUCCESS ) THEN
                  PRINT *, ' *** ERROR in grib_get: typeOfLevel: ', ireturn
                ENDIF
                IF (TRIM(ytypeOfLevel) == 'hybrid') THEN
                  CALL grib_get (igribid, 'level', ilev,  ireturn)
                  IF( ilev == 1) THEN
                    yshortname = 'FIS_SH'
                  ENDIF
                ENDIF
              ELSEIF (TRIM(yshortname) == 'unknown') THEN
                ! Test whether it is GFS HGT (surface height in gpm)
                ! with Grib2 ids (0,3,5)
                CALL grib_get (igribid, 'discipline',        idisc, ireturn)
                CALL grib_get (igribid, 'parameterCategory', icatg, ireturn)
                CALL grib_get (igribid, 'parameterNumber',   ipara, ireturn)
                IF ( ireturn /= GRIB_SUCCESS ) THEN
                  PRINT *, ' *** ERROR in grib_get: dis,cat,par: ', ireturn
                ENDIF
                IF ((idisc==0) .AND. (icatg==3) .AND. (ipara==5)) THEN
                  IF (ytypeoflevel == 'surface') THEN
                    ! it is the height in gpm. As surface height it is taken as meters
                    yshortname = 'HSURF'
                  ENDIF
                ELSEIF ((idisc==2) .AND. (icatg==0) .AND. (ipara==0) .AND. lhir2lm) THEN
                  ! it is FRLAND from HIRLAM
                  yshortname = 'FR_LAND'
                ELSE
                  PRINT *, 'UNKNOWN yshortname with dis, cat, par:  ', idisc, icatg, ipara
                ENDIF
              ENDIF
            ENDIF

            lzrequired = .FALSE.
            get_loc_apix: DO izloc = 1, nvar_in
              IF ( (TRIM(var_in(izloc)%name) == TRIM(yshortname))  .AND.   &
                  (var_in(izloc)%dattyp(3:3) == 'E') ) THEN
                lzrequired = .TRUE.
                EXIT get_loc_apix
              ENDIF
            ENDDO get_loc_apix
          ENDIF

          IF (lzrequired) THEN
            ! Get size of field and unpack data
            CALL grib_get_size(igribid, 'values', inrpoints, ierrf)
            IF (inrpoints > nb_field_in) THEN
              PRINT *, ' *** ERROR: size of message is too big for allocated field: ', inrpoints, nb_field_in
              ierror = 2
              yerror = 'Error in grib_new_from_file'
              RETURN
            ENDIF
            IF (lgme2lm) THEN
              CALL grib_get(igribid, 'values', ds_gme)
              IF (idbg_level >= 10) THEN
                PRINT *, 'grib_api record ', TRIM(yshortname), ' is needed:  ', lzrequired, izloc, &
                                             MINVAL(ds_gme(1:inrpoints)), MAXVAL(ds_gme(1:inrpoints))
              ENDIF
            ELSE
              CALL grib_get(igribid, 'values', field_in_file)
              IF (idbg_level >= 10) THEN
                PRINT *, 'grib_api record ', TRIM(yshortname), ' is needed:  ', lzrequired, izloc, &
                                             MINVAL(field_in_file), MAXVAL(field_in_file)
              ENDIF
            ENDIF
          ENDIF

#endif

        ELSEIF (yinext_form_read == 'ncdf') THEN

          get_loc_ncdf: DO izloc = 1, nvar_in
            IF (TRIM(var_in(izloc)%name) == TRIM(ylistextpar(izvar_count)) ) THEN
              lzrequired = .TRUE.
              EXIT get_loc_ncdf
            ENDIF
          ENDDO get_loc_ncdf
          yzlocname = TRIM(var_in(izloc)%name)

          IF (idbg_level >= 10) THEN
            PRINT *, 'NetCDF record ', TRIM(yzlocname), ' is needed:  ', lzrequired, izloc
          ENDIF
        ENDIF   ! check required and get_loc

        IF (.NOT. lzrequired) THEN
          ! no variable has been found in the input variable table
          CYCLE read_loop
        ENDIF
        yzlocname = TRIM(var_in(izloc)%name)

      ELSE
        lzeof = .TRUE.
      ENDIF ! igriblen

      IF (idbg_level >= 10) THEN
        IF (igriblen > 0) THEN
          PRINT *, '     Got and keep record ', yzlocname, izloc, igriblen
        ELSE
          PRINT *, '     Got no more record: EOF '
        ENDIF
      ENDIF

    ENDIF   ! my_cart_id == 0

  !----------------------------------------------------------------------------
  ! 3.3: Distribute status of reading to all PEs
  !----------------------------------------------------------------------------

    IF (num_compute > 1) THEN
      IF (my_cart_id == 0) THEN
        iz_info(1) = igriblen
        iz_info(2) = izloc
        iz_info(3) = nzscan
      ENDIF
      CALL distribute_values (iz_info, 3, 0, imp_integers, icomm_cart, izerror)
      IF (my_cart_id /= 0) THEN
        igriblen  = iz_info(1)
        izloc     = iz_info(2)
        nzscan    = iz_info(3)
        IF (izloc /= 0) yzlocname = var_in(izloc)%name
      ENDIF
      IF (igriblen == 0) THEN
        lzeof = .TRUE.
      ENDIF
    ENDIF

    IF (lzeof) THEN
      ! no more records to process
      EXIT read_loop
    ENDIF

  !----------------------------------------------------------------------------
  ! 3.4: Distribute the record to all PEs
  !----------------------------------------------------------------------------

    IF (num_compute > 1) THEN
      IF (izdebug > 10) THEN
        PRINT *, '       distribute field to all PEs'
      ENDIF
      IF (lgme2lm) THEN
        CALL distribute_values (ds_gme        , nb_field_in, 0, imp_grib, &
                                        icomm_cart,izerror)
      ELSE
        CALL distribute_values (field_in_file,  nb_field_in, 0, imp_grib, &
                                        icomm_cart,izerror)
      ENDIF
    ENDIF

    ! convert field to the precision used in INT2LM
    IF (lgme2lm) THEN
!     ij = (igg2e-igg2s+1) * (igg1e-igg1s+1)
!           ijk = (k-1)*ij + (j-1)*(igg1e-igg1s+1) + i
      ij = 0
      DO k = 1, 10
        DO j = igg2s, igg2e
          DO i = igg1s, igg1e
            ij = ij+1
            field_gme(i,j,k) = ds_gme(ij)
          ENDDO
        ENDDO
      ENDDO
!     field_gme (igg1s:igg1e,igg2s:igg2e,1:10) =                            &
!                REAL (field_gme_file(igg1s:igg1e,igg2s:igg2e,1:10), ireals)
    ELSE
      DO j = 1, jz_extfile
        ! take into account the scanning mode
        ! with grib_api this could be read out of the grib
        IF     (nzscan ==  0) THEN
          ! j scans negatively (IFS, GSM)
          jscan = jz_extfile+1 - j
        ELSEIF (nzscan == 64) THEN
          ! j scans positively (most other models)
          jscan = j
        ELSE
          ! somethings very wrong here: must abort
          PRINT *, ' *** ERROR: wrong scanning mode for the data (must be 0 or 64): ', nzscan
          ierror = 3
          yerror = 'Error in scanning mode'
          RETURN
        ENDIF
        DO i = 1, iz_extfile
          ij = (j-1) * (ie_in_tot-east_add_in-west_add_in) + i
          IF (field_in_file(ij) /= zundef) THEN
            field_in(i+west_add_in,jscan+south_add_in) = REAL (field_in_file(ij) , ireals)
          ELSE
            field_in(i+west_add_in,jscan+south_add_in) = undef
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    IF (lcm2lm) THEN
      ! Set the extra lines of the field, if present
      IF (west_add_in /= 0) THEN
        ! Set extra western boundary with values from the east
        DO j = 1, je_in_tot
          field_in(1,j) = field_in(ie_in_tot-east_add_in,j)
        ENDDO
      ENDIF
      IF (east_add_in /= 0) THEN
        ! Set extra eastern boundary with values from the west
        DO j = 1, je_in_tot
          field_in(ie_in_tot,j) = field_in(1+west_add_in,j)
        ENDDO
      ENDIF
      IF (south_add_in /= 0) THEN
        ! Set extra southern boundary with values from the second southern line
        DO i = 1, ie_in_tot
          field_in(i,1) = field_in(i,2)
        ENDDO
      ENDIF
      IF (north_add_in /= 0) THEN
        ! Set extra northern boundary with values from the second northern line
        DO i = 1, ie_in_tot
          field_in(i,je_in_tot) = field_in(i,je_in_tot-1)
        ENDDO
      ENDIF
    ENDIF

  !----------------------------------------------------------------------------
  ! 3.5: Check the record
  !----------------------------------------------------------------------------

    IF (my_cart_id == 0) THEN

      ! A check of the grid information is only done for grib. For ncdf this
      ! has been done in read_nc_gdefs_ext_in before the read-loop
      ! For all formats, min-, max- and mean-values for all records are computed

      ! This is done only in PE 0, because it is the only one who has the meta data

      IF (linit) THEN
        ! Write a headline in YUCHKDAT for external parameters
        WRITE (nuchkdat,'(A,A)') 'Check the external parameters from ', yinput_model
        WRITE (nuchkdat,'(A,A)') '    File:   ',TRIM(yname)
        IF (lgme2lm) THEN
          WRITE (nuchkdat,'(A,I5,A,I5,A,I5)')                            &
           '    ni_gme =', ni_gme,'  ni2 =', ni2,'  ni3 =', ni3
          WRITE (nuchkdat,'(A)') '    '
          WRITE (nuchkdat,'(A,A)')                                       &
           '     var         ee  lev        min     ',                   &
           ' location         max         location        mean  '
          WRITE (nuchkdat,'(A,A)')                                       &
           '                                      ',                     &
           ' j1   j2   jd                 j1   j2   jd'
        ELSE
          WRITE (nuchkdat,'(A,I5,A,I5)')                                 &
           '    ie_in_tot  =', ie_in_tot,'  je_in_tot =', je_in_tot
          WRITE (nuchkdat,'(A)') '    '
          WRITE (nuchkdat,'(A,A)')                                       &
            '   var        ee  lev        min    ',                      &
            'imin jmin               max    imax jmax              mean'
        ENDIF
        WRITE (nuchkdat,'(A)') '  '
        linit = .FALSE.
      ENDIF


      ! GME has to be treated separately due to the different grid
      IF (lgme2lm) THEN
        vcoord_in%ivctype  = 1  ! nothing else implemented here

        ! Check the grid
        IF (yinext_form_read == 'grb1' .OR. yinext_form_read == 'apix') THEN
          izerrorgrid = 0
          IF (ni_gme /= iz_ni_gme) izerrorgrid = 1
          IF (ni2    /= iz_ni2   ) izerrorgrid = 2
          IF (ni3    /= iz_ni3   ) izerrorgrid = 3
          IF (nd     /= iz_nd    ) izerrorgrid = 4

          IF (izerrorgrid /= 0) THEN
            PRINT *, '  external data file        namelist input'
            PRINT *, 'ni_gme     ',iz_ni_gme,'       ',ni_gme
            PRINT *, 'ni2        ',iz_ni2   ,'       ',ni2
            PRINT *, 'ni3        ',iz_ni3   ,'       ',ni3
            PRINT *, 'nd         ',iz_nd    ,'       ',nd
            ierror = 4
            yerror = 'wrong grid for GME external parameters'
            RETURN
          ENDIF
        ENDIF

        ! Check the record
        IF (lchkin) THEN
          ! print the maximum, minimum and meanvalues of each record
          igmemin(:) =   MINLOC ( field_gme (igg1s:igg1e,igg2s:igg2e,1:10) )
          igmemax(:) =   MAXLOC ( field_gme (igg1s:igg1e,igg2s:igg2e,1:10) )
          zmin       =   MINVAL ( field_gme (igg1s:igg1e,igg2s:igg2e,1:10) )
          zmax       =   MAXVAL ( field_gme (igg1s:igg1e,igg2s:igg2e,1:10) )
          zsum       =   SUM    ( field_gme (igg1s:igg1e,igg2s:igg2e,1:10) )
          zmean      =   zsum / REAL((igg1e-igg1s+1)*(igg2e-igg2s+1)*10,ireals)

          ! print the values
          IF (TRIM(var_in(izloc)%name) == 'NDVI_MR') THEN
            ilev = imonth
          ELSE
            ilev = 1
          ENDIF
          WRITE(nuchkdat,'(A,A,   I4,F14.6,3I5,F14.6,3I5,F14.7)')             &
            '    ', yzlocname,                      ilev,                     &
            zmin, igmemin(1), igmemin(2), igmemin(3),                         &
            zmax, igmemax(1), igmemax(2), igmemax(3), zmean
        ENDIF ! lchkin

      ELSE  ! all other models have regular (rotated) lat-lon grid

        ! Check the grid
        IF (yinext_form_read == 'grb1' .OR. yinext_form_read == 'apix') THEN
          CALL check_input_grid (igrbednr, igds_in, ngds, ipds, npds, igribid,    &
                   yinext_form_read, yzlocname, yzfulldate,                       &
                   ie_in_tot, je_in_tot, ke_in+nlevskip, startlat_in_tot,         &
                   startlon_in_tot, dlon_in, dlat_in, pollon_in, pollat_in,       &
                   inrvert_in, pv_in, (.NOT.lzeof) .AND. (lzrequired), 1, -1, 0,  &
                   .FALSE., itype_calendar, yinput_model, yzerror, izerror)
        ENDIF

        ! Check the record
        IF (lchkin) THEN
          ! print the maximum, minimum and meanvalues of each record
          i_inmin(:) =   MINLOC ( field_in (:,:) )
          i_inmax(:) =   MAXLOC ( field_in (:,:) )
          zmin       =   MINVAL ( field_in (:,:) )
          zmax       =   MAXVAL ( field_in (:,:) )
          zsum       =   SUM    ( field_in (:,:) )
          zmean      =   zsum / REAL (ie_in_tot * je_in_tot, ireals)

          ! print the values

          WRITE(nuchkdat,'(A,A,I4,I4,F15.6,2I5,A,F16.6,2I5,A,F16.7)')      &
            '  ', yzlocname, var_in(izloc)%ee,          1,                 &
                  zmin, i_inmin(1), i_inmin(2), '     ',                   &
                  zmax, i_inmax(1), i_inmax(2), '     ', zmean
        ENDIF ! lchkin

      ENDIF   ! lgme2lm or other models

    ENDIF     ! my_cart_id == 0

  !----------------------------------------------------------------------------
  ! 3.6: Extract vertical coordinate parameters
  !----------------------------------------------------------------------------

    ! NOTE: some of the vertical coordinate parameters are required for lm2lm
    !       since we need to do some processing below (e.g. vcoord_in%ivctype for SLEVE input)

    IF (llm2lm .AND. lzgetvertcoord) THEN

      IF (my_cart_id == 0) THEN

        ! Check whether we have got vertical coordinates (igds_in(2) /= 0)
        IF     (yinext_form_read == 'grb1') THEN
          izvert = INT (igds_in(2), iintegers)
#ifdef GRIBAPI
        ELSEIF (yinext_form_read == 'apix') THEN
          CALL grib_get (igribid, 'numberOfVerticalCoordinateValues', izvert, ireturn)
#endif
        ENDIF
  
        IF (izvert > 0) THEN
    
          ALLOCATE(pv(izvert), STAT=ireturn)
      
          ! PE with rank izsender got vertical coordinates, unpacks them (for dwdgrib1)
          ! and sends them to the others
          IF     (yinext_form_read == 'grb1') THEN
            ! check for the old or new style of coding (and how many can be decoded)
            IF (REAL(REFSTF(igds_in(26)), ireals) > 500) THEN
              ! it is the old style of coding
              DO n = 1, ke1in+4
                pv(n) = REAL (REFSTF(igds_in(25+n)), ireals)
              ENDDO
              IF (izvert > ke1in+4) THEN
                ! additional values are coded
                pv(ke1in+4+1) = REAL(igds_in(25+ke1in+4+1), ireals)  ! this is vcoord_in%ivctype
                DO n = ke1in+4+2, izvert
                  IF (igds_in(25+n) /= -1) pv(n) = REAL(REFSTF(igds_in(25+n)), ireals)
                ENDDO
              ENDIF
            ELSE
              ! it is the new style of coding
              DO n = 1, ke1in+6
                pv(n) = REAL (REFSTF(igds_in(25+n)), ireals)
              ENDDO
              DO n = ke1in+6+1, izvert
                IF (igds_in(25+n) /= -1) THEN
                  pv(n) = REAL(REFSTF(igds_in(25+n)), ireals)
                ELSE
                  pv(n) = 0.0_ireals
                ENDIF
              ENDDO
            ENDIF

#ifdef GRIBAPI
          ELSEIF (yinext_form_read == 'apix') THEN
            CALL grib_get (igribid, 'pv', pv, ireturn)
#endif
          ENDIF

          ! process the vertical coordinate parameters stored in pv
  
          ! Check for the type of the vertical coordinate parameters
          ! (In the old style of coding, this value is p0sl, so rather big,
          !  or it is the new style of coding, then it is a small integer)
          idummy = NINT(pv(1), iintegers)
  
          IF ((idummy >= 1) .AND. (idummy <= 300)) THEN
  
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
              PRINT *,  ' ERROR *** Type vcoord_in%ivctype of vertical coordinate '// &
                                                            'not available***'
              izerror = 10
            ENDIF
  
            ! This is the new grib GDS coding style introduced with INT2LM 1.5
            refatm_in%p0sl   = pv( 3)
            refatm_in%t0sl   = pv( 4)
            refatm_in%dt0lp  = pv( 5)
            vcoord_in%vcflat = pv( 6)
  
            DO k = 1, ke1in
              vcoord_in%vert_coord(k)  = pv(6+k)
            ENDDO
  
            IF (vcoord_in%ivctype == 3 .OR. vcoord_in%ivctype == 4) THEN
              ! read three more SLEVE parameters
              svc1_in   = pv(6 + ke1in + 1)
              svc2_in   = pv(6 + ke1in + 2)
              nfltvc_in = NINT(pv(6 + ke1in + 3), iintegers)
            ENDIF
  
            IF (refatm_in%irefatm == 2) THEN
              ! read additional parameters for new reference atmosphere
              refatm_in%delta_t = pv(6 + ke1in + 4)
              refatm_in%h_scal  = pv(6 + ke1in + 5)
            ENDIF
  
          ELSE
            ! This is the old grib GDS coding style, kept for backwards
            ! compatibility.
            refatm_in%p0sl   = pv( 1)
            refatm_in%t0sl   = pv( 2)
            refatm_in%dt0lp  = pv( 3)
            vcoord_in%vcflat = pv( 4)
  
            DO k = 1, ke1in
              vcoord_in%vert_coord(k)  = pv(4 + k)
            ENDDO
  
            ! Check for the type of the vertical coordinate parameters
            ! again: this was at a different location then
            vcoord_in%ivctype = INT(igds_in(29 + ke1in + 1), iintegers)
  
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
              IF ( vcoord_in%vert_coord(2) > vcoord_in%vert_coord(1) ) THEN
                vcoord_in%ivctype = 1
              ELSEIF ( vcoord_in%vert_coord(2) < vcoord_in%vert_coord(1) ) THEN
                vcoord_in%ivctype = 2
              ELSE
                PRINT *,  ' ERROR *** Type vcoord_in%ivctype of vertical coordinate '// &
                                                            'not available***'
                izerror = 10
              ENDIF
            ELSEIF (vcoord_in%ivctype == 3 .OR. vcoord_in%ivctype == 4 ) THEN
              ! read three more SLEVE parameters
              svc1_in   = pv(4 + ke1in + 2)
              svc2_in   = pv(4 + ke1in + 3)
              nfltvc_in = NINT(pv(4 + ke1in + 4), iintegers)
            ENDIF
  
            IF (refatm_in%irefatm == 2) THEN
              ! read additional parameters for new reference atmosphere
              refatm_in%delta_t = pv(4 + ke1in + 5)
              refatm_in%h_scal  = pv(4 + ke1in + 6)
            ENDIF
  
          ENDIF
  
          ! Now all vertical coordinate parameters and vcoord_in%ivctype are set
          ! go on with some checks
  
          IF     (vcoord_in%ivctype == 1) THEN
            ! For this type the vertical coordinates should be ascending
            IF ( vcoord_in%vert_coord(2) < vcoord_in%vert_coord(1) ) THEN
              PRINT *,                                                         &
                ' ERROR *** Vertical coordinates not ascending for type *** ', &
                 vcoord_in%ivctype, vcoord_in%vert_coord(1), vcoord_in%vert_coord(2)
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
                              'range for vcoord_in%ivctype = 3 or 4 ***'
              izerror = 7
            ENDIF
  
            IF ((svc2_in > vcoord_in%vert_coord(1)) .OR. (svc2_in < 0.0)) THEN
              PRINT *,  ' ERROR *** svc2_in not in allowed '//                &
                               'range for vcoord_in%ivctype = 3 or 4 ***'
              izerror = 8
            ENDIF
  
            IF (nfltvc_in <= 0) THEN
              PRINT *,  ' ERROR *** nfltvc_in must be greater than '//        &
                                              'or equal to zero ***'
              izerror = 9
            ENDIF
          ELSE
            PRINT *,  ' ERROR *** Type vcoord_in%ivctype of vertical coordinate '// &
                                                        'not available***'
            izerror = 11
          ENDIF
  
          DEALLOCATE(pv)
  
        END IF
  
        WRITE (noutput,'(A)')         '     '
        WRITE (noutput,'(A)')         '     Vertical coordinate parameters of input:'
        WRITE (noutput,'(A)')         '     '
        WRITE (noutput,'(A,I10)')     '     vcoord_in%ivctype = ', vcoord_in%ivctype
        WRITE (noutput,'(A,I10)')     '     refatm_in%irefatm = ', refatm_in%irefatm
        IF (vcoord_in%ivctype == 3 .OR. vcoord_in%ivctype == 4) THEN
          WRITE (noutput,'(A,F13.2)')   '     svc1_in   = ', svc1_in
          WRITE (noutput,'(A,F13.2)')   '     svc2_in   = ', svc2_in
          WRITE (noutput,'(A,I11)')     '     nfltvc_in = ', nfltvc_in
        ENDIF
        WRITE (noutput,'(A)')         '     '
        IF (refatm_in%irefatm == 2) THEN
          WRITE (noutput,'(A,F10.2)') '     refatm_in%delta_t = ', refatm_in%delta_t
          WRITE (noutput,'(A,F10.2)') '     refatm_in%h_scal  = ', refatm_in%h_scal
          WRITE (noutput,'(A)')       '     '
        ENDIF
        WRITE (noutput,'(A)')         '     k        vcoord_in(k)'
        DO k = 1,ke_in+1
          WRITE (noutput,'(I6,F16.4)') k, vcoord_in%vert_coord(k)
        ENDDO
        WRITE (noutput,'(A)')         '     '

      END IF ! my_cart_id == 0

      lzgetvertcoord = .FALSE.
      CALL distribute_values (lzgetvertcoord, 1, 0, imp_logical, icomm_cart,izerror)
      CALL distribute_values (vcoord_in%ivctype, 1, 0, imp_integers, icomm_cart,izerror)
      CALL distribute_values (refatm_in%irefatm, 1, 0, imp_integers, icomm_cart,izerror)
      CALL distribute_values (refatm_in%p0sl, 1, 0, imp_reals, icomm_cart,izerror)
      CALL distribute_values (refatm_in%t0sl, 1, 0, imp_reals, icomm_cart,izerror)
      CALL distribute_values (refatm_in%dt0lp, 1, 0, imp_reals, icomm_cart,izerror)
      CALL distribute_values (vcoord_in%vcflat, 1, 0, imp_reals, icomm_cart,izerror)
      CALL distribute_values (vcoord_in%vert_coord, ke1in, 0, imp_reals, icomm_cart,izerror)
      IF (vcoord_in%ivctype == 3 .OR. vcoord_in%ivctype == 4) THEN
        CALL distribute_values (svc1_in, 1, 0, imp_reals, icomm_cart,izerror)
        CALL distribute_values (svc2_in, 1, 0, imp_reals, icomm_cart,izerror)
        CALL distribute_values (nfltvc_in, 1, 0, imp_integers, icomm_cart,izerror)
      ENDIF
      IF (refatm_in%irefatm == 2) THEN
        CALL distribute_values (refatm_in%delta_t, 1, 0, imp_reals, icomm_cart,izerror)
        CALL distribute_values (refatm_in%h_scal, 1, 0, imp_reals, icomm_cart,izerror)
      ENDIF

    END IF ! lm2lm .AND. (.NOT. lzgetvertcoord)

  !----------------------------------------------------------------------------
  ! 3.7: Put values to memory and set logical flags
  !----------------------------------------------------------------------------

    IF (izdebug > 10) THEN
      PRINT *, '       put field to memory', izloc
    ENDIF

    IF (lgme2lm) THEN
      ! Extend the field on the boundary lines, scale it with factor and bias
      ! and put it into memory
      CALL xd(field_gme ,igg1sm2, igg1ep2, igg2sm2, igg2ep2, 1, 1, 1, 1, 1,  &
              10, igg1s, igg1e, igg2s, igg2e, 1, 10, 2, .TRUE., .FALSE.,     &
              field_gme ,igg1sm2, igg1ep2, igg2sm2, igg2ep2, 10, undef, izerror)

      ! this is only for DWD Grib1 table!!!
      zbias     = var_in (izloc)%bias
      zfactor   = var_in (izloc)%factor
      IF (TRIM(var_in(izloc)%name) == 'NDVI_MR') THEN
        ! the only external GME field with rank 4 but without bias and factor
        var_in(izloc)%p4(:,:,jd_min:jd_max,imonth) = field_gme (:,:,jd_min:jd_max) &
                                                         / zfactor - zbias
      ELSE
        var_in(izloc)%p3(:,:,jd_min:jd_max) = field_gme (:,:,jd_min:jd_max)  &
                                                         / zfactor - zbias
      ENDIF
    ELSE
      IF     (yinext_form_read == 'grb1') THEN
        zbias     = var_in (izloc)%bias
        zfactor   = var_in (izloc)%factor
      ELSEIF (yinext_form_read == 'apix') THEN
        zbias     = var_in (izloc)%bias
        zfactor   = var_in (izloc)%factor
      ELSE
        zbias     = 0.0_ireals
        zfactor   = 1.0_ireals
      ENDIF

      ! Cut out the proper global coarse grid LM domain from field_in
      IF (lcm2lm) THEN
        ! determine the shift with isubpos_coarse
        izshift = isubpos_coarse(my_cart_id, 1) - 1
        jzshift = isubpos_coarse(my_cart_id, 2) - 1
      ELSE
        IF (startlon_in - startlon_in_tot < 0.0) THEN
          izshift = NINT ( (startlon_in - startlon_in_tot + 360.0_ireals)/dlon_in )
        ELSE
          izshift = NINT ( (startlon_in - startlon_in_tot)/dlon_in )
        ENDIF
        jzshift = NINT ( (startlat_in - startlat_in_tot)/dlat_in )
      ENDIF

      IF     (yinext_form_read == 'grb1') THEN
        DO j = 1, je_in
          DO i = 1, ie_in
            var_in(izloc)%p2(i,j) = field_in(i+izshift,j+jzshift)
          ENDDO
        ENDDO
      ELSEIF (yinext_form_read == 'apix') THEN
        DO j = 1, je_in
          DO i = 1, ie_in
            var_in(izloc)%p2(i,j) = field_in(i+izshift,j+jzshift)
          ENDDO
        ENDDO
      ELSEIF (yinext_form_read == 'ncdf') THEN
        DO j = 1, je_in
          DO i = 1, ie_in
            var_in(izloc)%p2(i,j) = field_in(i+izshift,j+jzshift)
          ENDDO
        ENDDO
      ENDIF

      IF (vcoord_in%ivctype == 3 .OR. vcoord_in%ivctype == 4) THEN
        IF (llm2lm .AND. (TRIM(yzlocname) == 'HSURF')) THEN
          ! we split the coarse grid topography for the SLEVE coordinate here,
          ! because here we have access to it on the WHOLE coarse grid area

          ! allocate splitted topography parts
          ALLOCATE (hsurfs_in  (ie_in,je_in,2), STAT=izstat)
          hsurfs_in   = 0.0_ireals
          izerror = izstat

          ALLOCATE (zhsurfs_in_tot (ie_in_tot, je_in_tot, 2),    STAT=izstat)
          zhsurfs_in_tot = 0.0_ireals
          izerror = izerror + izstat

          IF (izerror /= 0) THEN
            ierror  = 2
            yerror  = 'Memory allocation failed for coarse grid SLEVE field(s)'
            RETURN
          ENDIF

          ! write Info to OUTPUT
          IF (my_cart_id == 0) THEN
            write(noutput,'(A,A)') '   Splitting of coarse LM grid ',          &
                                 'topography for SLEVE coordinate'
          ENDIF

          ! split coarse grid topography (field_in) to large scale part 
          ! (zhsurfs_in_tot(:,:,1) and small scale part (zhsurfs_in_tot(:,:,2).
          ! This is done on every PE, so no gathering/distributing is necessary
          CALL sleve_split_oro                                                        &
               (field_in, zhsurfs_in_tot, ie_in_tot, je_in_tot, nfltvc_in,            &
               0_iintegers, svc1_in, svc2_in, vcoord_in%vcflat, noutput, my_cart_id,  &
               izerror, yzerror)

          IF (izerror /= 0_iintegers) THEN
            ierror = 12
            yerror(1:LEN_TRIM(yzerror)) = yzerror(1:LEN_TRIM(yzerror))
            RETURN
          ENDIF

          ! cut out appropriate part of hsurfs_in (like above for field_in)
          DO j = 1, je_in
            DO i = 1, ie_in
              hsurfs_in(i,j,:) = zhsurfs_in_tot(i+izshift,j+jzshift,:)
            ENDDO
          ENDDO

          ! deallocate temporary field
          DEALLOCATE(zhsurfs_in_tot)
        ENDIF
      ELSE
        ! allocate a dummy variable, so that when passed to reference_atmosphere() not NULL-pointer
        ! error will be triggered
! SP, 201405: this line produces a "double allocate" error, reason unclear
!        ALLOCATE(hsurfs_in(1,1,1))
      ENDIF

    ENDIF
   
    ! Set corresponding logical flag to .TRUE.
    SELECT CASE (TRIM(yzlocname))
    CASE ('HSURF')
      lfis__in = .TRUE.
      fis_in(:,:) = hsurf_in(:,:) * g
    CASE ('FIS','FIS_SH')
      lfis__in = .TRUE.
    CASE ('FR_LAND') 
      lfrla_in = .TRUE.
      ! save the global field for computation of isolated points
      IF (.NOT. lgme2lm) THEN
        ALLOCATE (fland_in_tot(ie_in_tot,je_in_tot), STAT=izerror)
        fland_in_tot(:,:) = field_in(:,:)
      ENDIF
    CASE ('SOILTYP') 
      lsoty_in = .TRUE.
    CASE ('NDVI_MR') 
      lndvi_mr(imonth) = .TRUE.
    CASE ('Z0') 
      lz0___in = .TRUE.
    CASE ('PLCOV_MX') 
      lplmx_in = .TRUE.
    CASE ('PLCOV_MN') 
      lplmn_in = .TRUE.
    CASE ('LAI_MX') 
      laimx_in = .TRUE.
    CASE ('LAI_MN') 
      laimn_in = .TRUE.
    CASE ('ROOTDP') 
      lroot_in = .TRUE.
    END SELECT

#ifdef GRIBAPI
    IF (yinext_form_read=='apix') THEN
      ! The handle can be released now
      CALL grib_release(igribid)
    ENDIF
#endif

  ENDDO read_loop

!------------------------------------------------------------------------------
! Section 4: Cleanup
!------------------------------------------------------------------------------
 
  ! Deallocate fields
  DEALLOCATE (iblock, ibmap, dsup, ds_grib_single)
  IF (lgme2lm) THEN
    DEALLOCATE (field_gme, field_gme_file, ds_gme)
  ELSEIF (llm2lm .OR. lec2lm) THEN
    DEALLOCATE (field_in, field_in_file)
  ENDIF

  IF (lread) THEN
    ! grib file has been opened, so it is closed now
    CALL close_file (nufile, yinext_form_read, icomm_cart, my_cart_id,    &
                     num_compute, lasync_io, idbg_level, yerrmsg, izerror)
    IF (izerror /= 0) THEN
      ierror = 12
      yerror = 'Error in close_file'
      RETURN
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE read_coarse_grid_ext

!==============================================================================
!+ Module procedure in "src_read_ext" for filtering the orography
!------------------------------------------------------------------------------

SUBROUTINE low_pass_filter (xy_vec, a, b, c, d, e, f, ci, cj, ck, cl, cm,   &
                            ndim, eps)

!------------------------------------------------------------------------------
!
! Description:
!   The filter of Raymond (Raymond, W.H., 1988: High-Order Low-Pass Implicit
!   Tangent Filters for Use in Finite Area Calculations, MWR 116, 2132-2141)
!   is used for filtering orography. We use tenth order filter with p=5.
!   Computations are as proposed in the appendix of the paper.
!   Some corrections to the filtered field are done afterwards.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist

INTEGER (KIND=iintegers), INTENT(IN)    ::   &
  ndim         ! dimension of xy_vec

REAL    (KIND=ireals)   , INTENT(INOUT) ::   &
  xy_vec(ndim) ! one dimensional vector to be filtered

REAL    (KIND=ireals)   , INTENT(IN)    ::   &
  eps,              &   ! filter parameter
  ! variables for gaussian elimination
  a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim),    &
           ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim)

! Local variables
REAL    (KIND=ireals)                   ::   &
  rhs(ndim),   & ! RHS of equation
  phi(ndim),   & ! filter increment
    h(ndim)      ! for intermediate storage

INTEGER (KIND=iintegers)                :: n

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE low_pass_filter
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations and checks
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 2: Specify RHS of equation
!------------------------------------------------------------------------------

  rhs(     1) = 0
  rhs(ndim  ) = 0

  rhs(     2) = eps*(xy_vec(     1)-2.0*xy_vec(     2)+xy_vec(     3))
  rhs(ndim-1) = eps*(xy_vec(ndim-2)-2.0*xy_vec(ndim-1)+xy_vec(ndim  )) 

  rhs(     3) = eps*(-1.0*(xy_vec(     1)+xy_vec(     5))          & 
                     +4.0*(xy_vec(     2)+xy_vec(     4))          &
                     -6.0* xy_vec(     3)               )
  rhs(ndim-2) = eps*(-1.0*(xy_vec(ndim  )+xy_vec(ndim-4))          &
                     +4.0*(xy_vec(ndim-1)+xy_vec(ndim-3))          &
                     -6.0* xy_vec(ndim-2)         )

  rhs(     4) = eps*(      xy_vec(   1)+xy_vec(   7)           &
                     -6.0*(xy_vec(   2)+xy_vec(   6))          &
                    +15.0*(xy_vec(   3)+xy_vec(   5))          &
                    -20.0* xy_vec(   4)             )
  rhs(ndim-3) = eps*(      xy_vec(ndim-6)+xy_vec(  ndim)           &
                     -6.0*(xy_vec(ndim-5)+xy_vec(ndim-1))          &
                    +15.0*(xy_vec(ndim-4)+xy_vec(ndim-2))          &
                    -20.0* xy_vec(ndim-3)               )

  rhs(     5) = eps*(-1.0*(xy_vec(   1)+xy_vec(   9))          &
                     +8.0*(xy_vec(   2)+xy_vec(   8))          &
                    -28.0*(xy_vec(   3)+xy_vec(   7))          &
                    +56.0*(xy_vec(   4)+xy_vec(   6))          &
                    -70.0* xy_vec(   5)             )
  rhs(ndim-4) = eps*(-1.0*(xy_vec(ndim-8)+xy_vec(  ndim))          &
                     +8.0*(xy_vec(ndim-7)+xy_vec(ndim-1))          &
                    -28.0*(xy_vec(ndim-6)+xy_vec(ndim-2))          &
                    +56.0*(xy_vec(ndim-5)+xy_vec(ndim-3))          &
                    -70.0* xy_vec(ndim-4)               )
  
  DO n = 6, ndim-5
    rhs(n) = eps*(         xy_vec( n-5)+xy_vec( n+5)           &
                    -10.0*(xy_vec( n-4)+xy_vec( n+4))          &
                    +45.0*(xy_vec( n-3)+xy_vec( n+3))          &
                   -120.0*(xy_vec( n-2)+xy_vec( n+2))          &
                   +210.0*(xy_vec( n-1)+xy_vec( n+1))          &
                   -252.0* xy_vec(   n)             )
  ENDDO
  
!------------------------------------------------------------------------------
!- Section 3: Compute filter increment and filtered variables
!------------------------------------------------------------------------------

  ! Step III
  
  h(1) =  rhs(1)/f(1)
  h(2) = (rhs(2)-h(1)*e(2))/f(2)
  h(3) = (rhs(3)-h(1)*d(3)-h(2)*e(3))/f(3)
  h(4) = (rhs(4)-h(1)*c(4)-h(2)*d(4)-h(3)*e(4))/f(4)
  h(5) = (rhs(5)-h(1)*b(5)-h(2)*c(5)-h(3)*d(5)-h(4)*e(5))/f(5)
  
  DO n = 6, ndim
    h(n) = (rhs(n)-h(n-5)*a(n)-h(n-4)*b(n)-h(n-3)*c(n)-h(n-2)*d(n)- &
            h(n-1)*e(n))/f(n)
  ENDDO
  
  ! Step IV
  ! Determine filter increment phi 
  phi(ndim  ) = h(ndim)
  phi(ndim-1) = h(ndim-1)-ci(ndim-1)*phi(ndim)
  phi(ndim-2) = h(ndim-2)-ci(ndim-2)*phi(ndim-1)-cj(ndim-2)*phi(ndim)
  phi(ndim-3) = h(ndim-3)-ci(ndim-3)*phi(ndim-2)-cj(ndim-3)*phi(ndim-1)-   &
                          ck(ndim-3)*phi(ndim)
  phi(ndim-4) = h(ndim-4)-ci(ndim-4)*phi(ndim-3)-cj(ndim-4)*phi(ndim-2)-   &
                          ck(ndim-4)*phi(ndim-1)-cl(ndim-4)*phi(ndim)
  
  DO n = ndim-5,1,-1
    phi(n) = h(n)-ci(n)*phi(n+1)-cj(n)*phi(n+2)-ck(n)*phi(n+3)-      &
             cl(n)*phi(n+4)-cm(n)*phi(n+5)
  ENDDO

  DO n = 1, ndim
    xy_vec(n) = xy_vec(n) + phi(n)
  ENDDO
   
END SUBROUTINE low_pass_filter

!==============================================================================
!+ Module procedure in "src_read_ext" for initializing the gauss elimination
!------------------------------------------------------------------------------

SUBROUTINE set_gauss (a, b, c, d, e, f, ci, cj, ck, cl, cm, ndim, eps)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist

INTEGER (KIND=iintegers), INTENT(IN)    ::   &
  ndim         ! dimension of xy_vec

REAL    (KIND=ireals)   , INTENT(IN)    ::   &
  eps                   ! filter parameter

REAL    (KIND=ireals)   , INTENT(OUT)   ::   &
  a(ndim),  b(ndim),  c(ndim),  d(ndim),  e(ndim),  f(ndim),    &
           ci(ndim), cj(ndim), ck(ndim), cl(ndim), cm(ndim)

! Local variables
REAL    (KIND=ireals)                   ::   &
  o(ndim),  p(ndim),  q(ndim),  r(ndim),  s(ndim),  t(ndim),    &
  u(ndim),  v(ndim),  w(ndim),  x(ndim),  h(ndim)

INTEGER (KIND=iintegers)                :: n

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE set_gauss
!------------------------------------------------------------------------------

  ! Setup up the gaussian elimination step
  ! coefficients of matrix
  
   a(:) = 0.0;  b(:) = 0.0;  c(:) = 0.0;  d(:) = 0.0;  e(:) = 0.0; f(:) = 0.0;
  ci(:) = 0.0; cj(:) = 0.0; ck(:) = 0.0; cl(:) = 0.0; cm(:) = 0.0;
   o(:) = 0.0;  p(:) = 0.0;  q(:) = 0.0;  r(:) = 0.0;  s(:) = 0.0;  
   t(:) = 0.0;  u(:) = 0.0;  v(:) = 0.0;  w(:) = 0.0;  x(:) = 0.0;
  
  s(   1) = 1.0
  s(ndim  ) = s(1)
  s(   2) = 2.0*(1.0+eps)
  s(ndim-1) = s(2)
  s(   3) = 6.0*(1.0+eps)
  s(ndim-2) = s(3)
  s(   4) =20.0*(1.0+eps)
  s(ndim-3) = s(4)
  s(   5) =70.0*(1.0+eps)
  s(ndim-4) = s(5)
  
  t(   2) = 1.0-eps
  t(ndim-1) = t(2)
  t(   3) = 4.0*(1.0-eps)
  t(ndim-2) = t(3)
  t(   4) =15.0*(1.0-eps)
  t(ndim-3) = t(4)
  t(   5) =56.0*(1.0-eps)
  t(ndim-4) = t(5)
  
  u(   3) = 1.0+eps
  u(ndim-2) = u(3)
  u(   4) = 6.0*(1.0+eps)
  u(ndim-3) = u(4)
  u(   5) =28.0*(1.0+eps)
  u(ndim-4) = u(5)
  
  v(   4) = 1.0-eps
  v(ndim-3) = v(4)
  v(   5) = 8.0*(1.0-eps)
  v(ndim-4) = v(5)
  
  w(   5) = 1.0+eps
  w(ndim-4) = w(5)
  
  DO n = 6, ndim-5
    x(n) =  1.0-eps
    w(n) = 10.0*(1.0+eps)
    v(n) = 45.0*(1.0-eps)
    u(n) =120.0*(1.0+eps)
    t(n) =210.0*(1.0-eps)
    s(n) =252.0*(1.0+eps)
  ENDDO
  
  r(:) = t(:)
  q(:) = u(:)
  p(:) = v(:)
  o(:) = w(:)
  a(:) = x(:)
  
  !-------------Matrix inversion--------------
  
  ! Step I
  
  f (1) = s(1)
  ci(1) = t(1)/f(1)
  e (2) = r(2)
  f (2) = s(2)-e(2)*ci(1)
  cj(1) = u(1)/f(1)
  ci(2) = (t(2)-e(2)*cj(1))/f(2)
  d (3) = q(3)
  e (3) = r(3)-d(3)*ci(1)
  f (3) = s(3)-e(3)*ci(2)-d(3)*cj(1)
  ck(1) = v(1)/f(1)
  cj(2) = (u(2)-e(2)*ck(1))/f(2)
  ci(3) = (t(3)-d(3)*ck(1)-e(3)*cj(2))/f(3)
  c (4) = p(4)
  d (4) = q(4)-c(4)*ci(1)
  e (4) = r(4)-d(4)*ci(2)-c(4)*cj(1)
  f (4) = s(4)-e(4)*ci(3)-d(4)*cj(2)-c(4)*ck(1)
  cl(1) = w(1)/f(1)
  ck(2) = (v(2)-e(2)*cl(1))/f(2)
  cj(3) = (u(3)-d(3)*cl(1)-e(3)*ck(2))/f(3)
  ci(4) = (t(4)-c(4)*cl(1)-d(4)*ck(2)-e(4)*cj(3))/f(4)
  b (5) = o(5)
  c (5) = p(5)-b(5)*ci(1)
  d (5) = q(5)-c(5)*ci(2)-b(5)*cj(1)
  e (5) = r(5)-d(5)*ci(3)-c(5)*cj(2)-b(5)*ck(1)
  f (5) = s(5)-e(5)*ci(4)-d(5)*cj(3)-c(5)*ck(2)-b(5)*cl(1) 
  
  ! Step II
  DO n = 6, ndim
    cm(n-5) = x(n-5)/f(n-5)
    cl(n-4) = (w(n-4)-e(n-4)*cm(n-5))/f(n-4)
    ck(n-3) = (v(n-3)-d(n-3)*cm(n-5)-e(n-3)*cl(n-4))/f(n-3)
    cj(n-2) = (u(n-2)-c(n-2)*cm(n-5)-d(n-2)*cl(n-4)-e(n-2)*ck(n-3))/f(n-2)
    ci(n-1) = (t(n-1)-b(n-1)*cm(n-5)-c(n-1)*cl(n-4)-d(n-1)*ck(n-3)-  &
              e(n-1)*cj(n-2))/f(n-1)
    b (n  ) = o(n)-a(n)*ci(n-5)
    c (n  ) = p(n)-b(n)*ci(n-4)-a(n)*cj(n-5)
    d (n  ) = q(n)-c(n)*ci(n-3)-b(n)*cj(n-4)-a(n)*ck(n-5)
    e (n  ) = r(n)-d(n)*ci(n-2)-c(n)*cj(n-3)-b(n)*ck(n-4)-a(n)*cl(n-5)
    f (n  ) = s(n)-e(n)*ci(n-1)-d(n)*cj(n-2)-c(n)*ck(n-3)-b(n)*cl(n-4)- &
              a(n)*cm(n-5)
  ENDDO

END SUBROUTINE set_gauss

!==============================================================================
!==============================================================================
!+ Module procedure in "src_read_ext" for reading a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE read_nc_axis   (iedim, jedim, kedim,&
                           startlon, startlat, endlon, endlat, pollon, pollat, &
                           lushift_in, lvshift_in,   &
                           east_add_in, west_add_in, south_add_in, north_add_in, &
                           icomm, myid, npes, yerrmsg, istatus)
!
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
    icomm,             & ! MPI communicator
    myid,              & ! ID of this PE in communicator icomm
    npes,              & ! number of PEs
    iedim, jedim, kedim, &  ! dimensions of the input fields
    east_add_in,   & ! add an extra column to the East
    west_add_in,   & ! add an extra column to the West
    south_add_in,  & ! add an extra column to the South
    north_add_in     ! add an extra column to the North

  REAL    (KIND=ireals),    INTENT(IN)  ::  &
    startlat, startlon, & ! coordinates of the lower left grid point
    endlat, endlon,     & ! coordinates of the upper right grid point
    pollon, pollat        ! coordinates of the geographical North Pole

  LOGICAL, INTENT(IN) :: &
    lushift_in(2),  & ! shift of u-velocity due to grid staggering
    lvshift_in(2)     ! shift of v-velocity due to grid staggering

! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    istatus          ! error index

  CHARACTER (LEN= *),        INTENT(OUT)   ::  &
    yerrmsg     ! string for error messages

!-----------------------------------------------------------------------
!
! Local scalars:

INTEGER (KIND=iintegers), PARAMETER  :: ntbds=2, ntime=1

  INTEGER (KIND=iintegers) ::  &
    ierror, izmplcode,                       & ! error status variable
    iedim_in, jedim_in, kedim_in, ke1dim_in, & ! input dimensions
    ncid                                       ! NetCDF file ID

  INTEGER (KIND=iintegers) ::   &
    jgridVarID,     & ! NetCDF ID for grid mapping !_br 31.03.09
    jlonVarID,      & ! NetCDF ID for longitude
    jlatVarID,      & ! NetCDF ID for latitude
    jslonVarID,     & ! NetCDF ID for slongitude
    jslatVarID        ! NetCDF ID for slatitude

  CHARACTER (LEN=5) :: &
    ylon, ylat,     & ! names for variables langitude and latitude
    yslon, yslat      ! names for variables langitude and latitude

  CHARACTER (LEN=250)        ::  &
    yzname      ! full path and name of the input-file

  CHARACTER (LEN=80)         ::   yzerror

  CHARACTER (LEN=1)         ::   yzerror_star(4)

!  Local arrays:

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

#ifdef NETCDF
!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! initialize error status variable
  yerrmsg         = '   '
  yzerror_star(:) = ' '
  ierror          = 0
  istatus         = 0

  yzname = TRIM(yinext_cat)//TRIM(yinext_lfn)

  CALL open_file(ncid, yzname, ymode_read, yinext_form_read, icomm_cart, &
                   my_cart_id, num_compute, lasync_io, idbg_level,         &
                   yzerror, ierror)
  IF (ierror /= 0) THEN
    PRINT *, 'Error in open_nc : ', TRIM(yzerror)
    RETURN
  ENDIF

  ! Processor 0 does the job
  IF (myid == 0) THEN
    ! Get the dimension ID's and the length of the dimensions
    istatus = nf90_inq_varid(ncid, 'rotated_pole', jgridVarid)
    IF (istatus == NF90_NOERR) THEN
      ylon = 'rlon'
      ylat = 'rlat'
      yslon = 'srlon'
      yslat = 'srlat'
    ELSE
      ylon = 'lon'
      ylat = 'lat'
      yslon = 'slon'
      yslat = 'slat'
    ENDIF

    istatus = nf90_inq_dimid (ncid, TRIM(ylon), idims_id(1))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_inquire_dimension (ncid, idims_id(1), len=iedim_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF ((iedim - east_add_in - west_add_in) /= iedim_in) istatus= 1
    IF (istatus /= 0 ) THEN
      yerrmsg = 'Number of west-east grid points in GRID_IN and yinext_lfn differs'
      RETURN
    ENDIF
    istatus = nf90_inq_dimid (ncid, TRIM(ylat), idims_id(2))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_inquire_dimension (ncid, idims_id(2), len=jedim_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF ((jedim - south_add_in - north_add_in) /= jedim_in) istatus= 2
    IF (istatus /= 0 ) THEN
      yerrmsg = 'Number of south-north grid points in GRID_IN and yinext_lfn differs'
      RETURN
    ENDIF

    IF ( .NOT. lcm_pres_coor) THEN
      istatus = nf90_inq_dimid (ncid, 'level', idims_id(3))
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      istatus = nf90_inquire_dimension (ncid, idims_id(3), len=kedim_in)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      IF (kedim /= kedim_in) istatus= 3
      IF (istatus /= 0 ) THEN
        yerrmsg = 'Number of vertical levels in GRID_IN and yinext_lfn differs'
        RETURN
      ENDIF

      istatus = nf90_inq_dimid (ncid, 'level1', idims_id(4))
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      istatus = nf90_inquire_dimension (ncid, idims_id(4), len=ke1dim_in)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      IF ((kedim+1) /= ke1dim_in) istatus= 3
      IF (istatus /= 0 ) THEN
        yerrmsg = 'Number of vertical levels+1 in GRID_IN and yinext_lfn differs'
      RETURN
      ENDIF
    ENDIF


    istatus = nf90_inq_varid (ncid, TRIM(ylon), jlonVarID)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_get_var (ncid, jlonVarID, longitudes_in(1+west_add_in:iedim-east_add_in))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF (west_add_in /= 0) THEN
      longitudes_in(1) = longitudes_in(iedim - east_add_in) - 360._ireals
    ENDIF
    IF (east_add_in /= 0) THEN
      longitudes_in(iedim) = 360._ireals + longitudes_in(1 + west_add_in)
    ENDIF

    istatus = nf90_inq_varid (ncid, TRIM(ylat), jlatVarID)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_get_var (ncid, jlatVarID, latitudes_in(1+south_add_in:jedim-north_add_in))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF (south_add_in /= 0) THEN
      latitudes_in(1) = -90._ireals
    ENDIF
    IF (north_add_in /= 0) THEN
      latitudes_in(jedim) = 90._ireals
    ENDIF

    IF (ABS(latitudes_in(1+south_add_in)  - startlat) > 1.0E-3_ireals) THEN
      ierror = 1
      yzerror_star(1) = '*'
    ENDIF

    IF (ABS(longitudes_in(1+west_add_in)  - startlon) > 1.0E-3_ireals) THEN
      ierror = 1
      yzerror_star(2) = '*'
    ENDIF

    IF (ABS(latitudes_in(jedim-north_add_in)  - endlat) > 1.0E-3_ireals) THEN
      ierror = 1
      yzerror_star(3) = '*'
    ENDIF

    IF (ABS(longitudes_in(iedim-east_add_in)  - endlon) > 1.0E-3_ireals) THEN
      ierror = 1
      yzerror_star(4) = '*'
    ENDIF

    istatus=ierror

    IF (istatus /= 0 ) THEN
          PRINT *, '                data file                   namelist input'
          PRINT *, yzerror_star(1),'startlat_in ',latitudes_in(1+south_add_in) ,'   ',startlat
          PRINT *, yzerror_star(2),'startlon_in ',longitudes_in(1+west_add_in) ,'   ',startlon

          PRINT *, yzerror_star(3),'endlat_in   ',latitudes_in(jedim-north_add_in),'   ',endlat
          PRINT *, yzerror_star(4),'endlon_in   ',longitudes_in(iedim-east_add_in),'   ',endlon
          yerrmsg = 'axis boundaries in namelist GRID_IN and input file differ'
          RETURN
    ENDIF

    IF (lushift_in(1) .OR. lvshift_in(1)) THEN

      istatus = nf90_inq_dimid (ncid, TRIM(yslon), idims_id(9))
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      istatus = nf90_inq_varid (ncid, TRIM(yslon), jslonVarID)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      istatus = nf90_get_var (ncid, jslonVarID, slongitudes_in(1+west_add_in:iedim-east_add_in))
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        istatus= -1
        RETURN
      ENDIF
      IF (west_add_in /= 0) THEN
        slongitudes_in(1) = slongitudes_in(iedim - east_add_in) - 360._ireals
      ENDIF
      IF (east_add_in /= 0) THEN
        slongitudes_in(iedim) = 360._ireals + slongitudes_in(1 + west_add_in)
      ENDIF

    ELSE

      slongitudes_in(:) = longitudes_in(:)

    ENDIF

    IF (lushift_in(2) .OR. lvshift_in(2)) THEN

      istatus = nf90_inq_dimid (ncid, TRIM(yslat), idims_id(10))
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      istatus = nf90_inq_varid (ncid, TRIM(yslat), jslatVarID)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      istatus = nf90_get_var (ncid, jslatVarID, slatitudes_in(1+south_add_in:jedim-north_add_in))
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(istatus))
        istatus= -1
        RETURN
      ENDIF
      IF (south_add_in /= 0) THEN
        slatitudes_in(1) = -90._ireals
      ENDIF
      IF (north_add_in /= 0) THEN
        slatitudes_in(jedim) = 90._ireals
      ENDIF

    ELSE

      slatitudes_in(:) = latitudes_in(:)

    ENDIF

  ENDIF ! myid == 0

! send information to the other processors
  IF (npes > 1) THEN
    CALL MPI_BCAST (longitudes_in,  iedim, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (latitudes_in,   jedim, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (slongitudes_in, iedim, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (slatitudes_in,  jedim, imp_reals,    0, icomm, izmplcode)
  ENDIF
 
  CALL close_file (ncid, yinext_form_read, icomm_cart, my_cart_id,    &
                   num_compute, lasync_io, idbg_level, yzerror, ierror)
  IF (ierror /= 0) THEN
    PRINT *, 'Error in close_nc : ', TRIM(yzerror)
    RETURN
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------
#endif

END SUBROUTINE read_nc_axis

!==============================================================================
!==============================================================================
!+ Module procedure in "src_read_ext" for reading a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE read_nc_gdefs_ext_lm (ncid, iedim, jedim,  &
                           startlat, startlon, endlat, endlon,  &
                           dlon, dlat, pollon, pollat, polgam,  &
                           startlon_ext, startlat_ext,  &
                           icomm, myid, npes, yerrmsg, istatus)
!
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
    ncid,              & ! NetCDF file ID
    icomm,             & ! MPI communicator
    myid,              & ! ID of this PE in communicator icomm
    npes,              & ! number of PEs
    iedim, jedim         ! dimensions of the input fields


  REAL    (KIND=ireals),    INTENT(IN)  ::  &
    startlat, startlon, & ! coordinates of the lower left grid point
    endlat, endlon,     & ! coordinates of the upper right grid point
    pollon, pollat, polgam,  & ! coordinates of the geographical North Pole
    dlon,  dlat           ! grid resolution in lambda and phi direction

! Scalar arguments with intent(out):
  REAL    (KIND=ireals),    INTENT(OUT)  ::  &
    startlat_ext, startlon_ext    ! coordinates of the lower left grid point 

  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    istatus          ! error index

  CHARACTER (LEN= *),        INTENT(OUT)   ::  &
    yerrmsg     ! string for error messages

!-----------------------------------------------------------------------
!
! Local scalars:

INTEGER (KIND=iintegers), PARAMETER  :: ntbds=2, ntime=1

  INTEGER (KIND=iintegers) ::  &
    izerrorgrid, izmplcode,                  & ! error status variable
    nhori_in,                  &
    ierror,                    &
    iedim_in, jedim_in

  INTEGER (kind=iintegers) ::   &
    jgridVarID,     & ! NetCDF ID for grid mapping
    jsectVarID,     & ! NetCDF ID for sections for HORIZON
    jlonVarID,      & ! NetCDF ID for longitude
    jlatVarID         ! NetCDF ID for latitude

  REAL (KIND=ireals) :: &
    zpollat,       & ! latitude of North Pole in input data
    zpollon,       & ! longitude of North Pole in input data
    zpolgam          ! latitude of the rotated north pole

  CHARACTER (LEN=80) :: &
    grid_name        ! name of grid mapping

  REAL (KIND=ireals)         ::  &
    zlatll, zlonll, zlatur, zlonur   ! longitudes and latitudes of lower left and
                                     ! upper right corner

  REAL    (KIND=ireals)     ::  &
    zdlon, zdlat             ! grid resolution in latitude and longitude

  CHARACTER (LEN=4) :: &
    ylon, ylat       ! names for variables langitude and latitude

! Local arrays:
  REAL (kind=ireals), ALLOCATABLE   :: &
    zlongitude(:),  & ! rotated longitudes
    zlatitude(:)      ! rotated latitudes

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

#ifdef NETCDF
!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! initialize error status variable
  yerrmsg         = '   '
  grid_name(1:80) = ' '
  istatus         = 0

! Processor 0 does the job
  IF (myid == 0) THEN
    ! Get the dimension ID's and the length of the dimensions
!_br 06.04.09
    istatus = nf90_inq_varid(ncid, 'rotated_pole', jgridVarid)
    IF (istatus == NF90_NOERR) THEN
      ylon = 'rlon'
      ylat = 'rlat'
    ELSE
      ylon = 'lon'
      ylat = 'lat'
    ENDIF
!_br 06.04.09 end

    istatus = nf90_inq_dimid (ncid, TRIM(ylon), idims_id(1))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_inquire_dimension (ncid, idims_id(1), len=iedim_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    istatus = nf90_inq_dimid (ncid, TRIM(ylat), idims_id(2))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_inquire_dimension (ncid, idims_id(2), len=jedim_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    ALLOCATE (zlongitude(iedim_in), zlatitude(jedim_in), STAT=istatus)
    IF (istatus /= 0) THEN
      yerrmsg = 'Allocation error in read_nc_gdefs_ext_lm'
      RETURN
    ENDIF

    istatus = nf90_inq_varid (ncid, TRIM(ylon), jlonVarID)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_get_var (ncid, jlonVarID, zlongitude)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    istatus = nf90_inq_varid (ncid, TRIM(ylat), jlatVarID)
     IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_get_var (ncid, jlatVarID, zlatitude)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    IF(lradtopo) THEN
       ! Read multi-dimensional external parameter HORIZON
       istatus = nf90_inq_dimid (ncid, "nhori", idims_id(14)) ! gdm 2.12.2013
       IF (istatus /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(istatus))
          RETURN
       ENDIF
       istatus = nf90_inquire_dimension (ncid, idims_id(14), len=nhori_in) ! gdm 2.12.2013
       IF (istatus /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(istatus))
          RETURN
       ENDIF

       IF ((nhori) /= nhori_in) THEN
          ierror = 14
          PRINT *, '  external data file        namelist input'
          PRINT *, 'nhori_in    ',nhori_in  ,'nhori     ',nhori
       ENDIF

       IF (istatus /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(istatus))
          RETURN
       ENDIF
    ENDIF

 !  get grid mapping values

    istatus = nf90_inq_varid(ncid, 'rotated_pole', jgridVarid)

    IF (istatus == NF90_NOERR) THEN   ! rotated latitude-longitude
      istatus = nf90_get_att(ncid, jgridVarid, 'grid_mapping_name', grid_name)

      IF (istatus /= NF90_NOERR) THEN
        PRINT *, TRIM(NF90_strerror(istatus))
        yerrmsg = 'Error in read_nc_gdefs_ext_lm / nf90_get_att - attribute "grid_mapping_name"'
        RETURN
      ENDIF

      IF (grid_name(1:26) /= 'rotated_latitude_longitude') THEN
        yerrmsg = 'Error in read_nc_gdefs_ext_lm - invalid value for attribute ' &
                  //'"grid_mapping_name": '//TRIM(grid_name)
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF

      istatus = nf90_get_att(ncid, jgridVarid, 'grid_north_pole_latitude', zpollat)
      IF (istatus /= NF90_NOERR) THEN
        PRINT *, TRIM(NF90_strerror(istatus))
        yerrmsg = 'Error in read_nc_gdefs_ext_lm / nf90_get_att - '     &
                     //'attribute "grid_north_pole_latitude"'
        RETURN
      ENDIF

      istatus = nf90_get_att(ncid, jgridVarid, 'grid_north_pole_longitude', zpollon)
      IF (istatus /= NF90_NOERR) THEN
        PRINT *, TRIM(NF90_strerror(istatus))
        yerrmsg = 'Error in read_nc_gdefs_ext_lm / nf90_get_att - '     &
                      //'attribute "grid_north_pole_longitude"'
        RETURN
      ENDIF

      IF (polgam /= 0.0_ireals) THEN
        istatus = nf90_get_att(ncid, jgridVarid, 'north_pole_grid_longitude', zpolgam)
        IF (istatus /= NF90_NOERR) THEN
          PRINT *, TRIM(NF90_strerror(istatus))
          yerrmsg = 'Error in read_nc_gdefs / nf90_get_att - '          &
                        //'attribute "north_pole_grid_longitude"'
          RETURN
        ENDIF
      ELSE
        zpolgam = 0.0_ireals
      ENDIF

    ELSE IF (istatus == NF90_enotvar) THEN  ! true latitude-longitude
      zpollat =  90.0_ireals
      zpollon = 180.0_ireals
      zpolgam =   0.0_ireals
      istatus = 0_iintegers
    ELSE
       PRINT *, 'Error in read_nc_gdefs_ext_lm / nf90_inq_varid'
       PRINT *, 'Variable "rotated_pole"'
       RETURN
    ENDIF

    ! Check size and resolution of the external parameters against the
    ! LM grid and whether the LM domain (with 1 boundary line) is within
    ! these fields

    zlatll = zlatitude(1)
    zlonll = zlongitude(1)
    zlatur = zlatitude(jedim_in)
    zlonur = zlongitude(iedim_in)

    startlon_ext = zlonll
    startlat_ext = zlatll

    zdlat   = (zlatitude(jedim_in) - zlatitude(1))/ REAL(jedim_in-1)
    IF (zlongitude(iedim_in) - zlongitude(1) < 0.0_ireals) THEN
      ! If the area is located around the 180-Meridian, zlongitude(iedim_in) - zlongitude(1)
      ! will be negative and 360 degrees have to be added
      zdlon   = (zlongitude(iedim_in) - zlongitude(1) + 360.0_ireals) / (iedim_in-1)
    ELSE
      zdlon   = (zlongitude(iedim_in) - zlongitude(1))/ (iedim_in-1)
    ENDIF

    ! Do all the checks:  abort program in case of errors
    izerrorgrid = 0
    IF ( (iedim_in /= iedim) .OR. (jedim_in /= jedim) )     izerrorgrid = 1

    IF ( ((zlonll > startlon - dlon) .AND.                          &
             (ABS(zlonll-(startlon - dlon)) > 1.0E-3_ireals)) .OR.  &
         ((zlatll > startlat - dlat) .AND.                          &
             (ABS(zlonll-(startlon - dlon)) > 1.0E-3_ireals)))      &
                                                            izerrorgrid = 2
    IF (((zlonur < endlon) .AND.                                    &
             (ABS(zlonur - endlon) > 1.0E-3_ireals)) .OR.           &
        ((zlatur < endlat) .AND.                                    &
             (ABS(zlatur - endlat) > 1.0E-3_ireals)))               &
                                                            izerrorgrid = 3

    IF ( (ABS(zdlon - dlon) > 1.0E-5_ireals) .OR.    &
         (ABS(zdlat - dlat) > 1.0E-5_ireals)      )         izerrorgrid = 4
    IF (ABS(zpollat - pollat) > 1.0E-3_ireals) izerrorgrid = 5
    IF (ABS(zpollon - pollon) > 1.0E-3_ireals) izerrorgrid = 6
    IF (ABS(zpolgam - polgam) > 1.0E-3_ireals) izerrorgrid = 6
    IF (MINVAL(ABS(zlongitude-startlon)) > 1.0E-3_ireals*dlon) izerrorgrid = 7
    IF (MINVAL(ABS(zlatitude -startlat)) > 1.0E-3_ireals*dlat) izerrorgrid = 8

!------------------------------------------------------------------------------
! Section 4: If necessary, print errors
!------------------------------------------------------------------------------

    IF ( izerrorgrid /= 0 ) THEN
      PRINT *, '  external data file        namelist input'

      PRINT *, 'ie_ext     ',iedim_in  ,'       ',iedim
      PRINT *, 'je_ext     ',jedim_in  ,'       ',jedim

      PRINT *, 'lonll      ',zlonll ,'       ',startlon
      PRINT *, 'latll      ',zlatll ,'       ',startlat

      PRINT *, 'lonur      ',zlonur ,'       ',endlon
      PRINT *, 'latur      ',zlatur ,'       ',endlat

      PRINT *, 'dlon       ',zdlon  ,'       ',dlon
      PRINT *, 'dlat       ',zdlat  ,'       ',dlat

      PRINT *, 'pollat     ',zpollat,'       ', pollat
      PRINT *, 'pollon     ',zpollon,'       ', pollon
      PRINT *, 'polgam     ',zpolgam,'       ', polgam

      yerrmsg = 'wrong grid for LM external parameters'
      istatus= -1

      RETURN
    ENDIF

    DEALLOCATE (zlatitude, zlongitude, STAT=istatus)
    IF (istatus /= 0) THEN
      yerrmsg = 'Deallocation error in read_nc_gdefs_ext_lm'
      RETURN
    ENDIF

  ENDIF ! (myid == 0)

! send information to the other processors
  IF (npes > 1) THEN
    CALL MPI_BCAST (startlat_ext,  1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (startlon_ext,  1, imp_reals,    0, icomm, izmplcode)
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
#endif

END SUBROUTINE read_nc_gdefs_ext_lm

!==============================================================================
!==============================================================================
!+ Module procedure in "src_read_ext" for reading a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE read_nc_gdefs_ext_in (ncid, iedim, jedim, kedim,           &
                startlon, startlat, endlon, endlat, dlon, dlat,       &
                pollon, pollat, polgam,                               &
                east_add_in, west_add_in, south_add_in, north_add_in, &
                icomm, myid, npes, yerrmsg, istatus)
!
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
  INTEGER (KIND=iintegers),   INTENT(IN) :: &
    ncid,              & ! NetCDF file ID
    icomm,             & ! MPI communicator
    myid,              & ! ID of this PE in communicator icomm
    npes,              & ! number of PEs
    iedim, jedim, kedim, &  ! dimensions of the input fields
    east_add_in,   & ! add an extra column to the East
    west_add_in,   & ! add an extra column to the West
    south_add_in,  & ! add an extra column to the South
    north_add_in     ! add an extra column to the North

  REAL    (KIND=ireals),    INTENT(IN)  ::  &
    startlat, startlon, & ! coordinates of the lower left grid point
    endlat, endlon,     & ! coordinates of the upper right grid point
    dlon,   dlat,       & ! coordinate increments
    pollon, pollat,     & ! coordinates of the geographical North Pole
    polgam                !

! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    istatus          ! error index

  CHARACTER (LEN= *),        INTENT(OUT)   ::  &
    yerrmsg     ! string for error messages

!-----------------------------------------------------------------------
!
! Local scalars:

INTEGER (KIND=iintegers), PARAMETER  :: ntbds=2, ntime=1

INTEGER (KIND=iintegers) ::  &
  ierror, izmplcode,                       & ! error status variable
  iedim_in, jedim_in, kedim_in, ke1dim_in, & ! input dimensions
  ivctype_in, irefatm_in                     ! just for reading

INTEGER (KIND=iintegers) ::  &
  jgridVarID,      & ! NetCDF ID for grid mapping
  jvcVarID,        & ! NetCDF ID for the vertical component
  jlonVarID,       & ! NetCDF ID for longitude
  jlatVarID          ! NetCDF ID for latitude

REAL    (KIND=irealgrib) ::  &
  zpollat_in,      & ! latitude of North Pole in input data
  zpollon_in,      & ! longitude of North Pole in input data
  zpolgam_in         ! angle between the north poles of the systems

CHARACTER (LEN=80)       ::  &
  grid_name_in       ! name of grid mapping

REAL    (KIND=irealgrib)     ::  &
  zdlon_in, zdlat_in ! grid resolution in latitude and longitude

REAL (KIND=irealgrib), ALLOCATABLE   :: &
  vcoordgrib(:),   & ! vertical coordinate
  longitude(:),    & ! rotated longitudes
  latitude(:)        ! rotated latitudes

REAL (KIND=irealgrib)                :: &
  p0slgrib, t0slgrib, dt0lpgrib, vcflatgrib, &
  deltatgrib, hscalgrib, svc1grib, svc2grib !_br 06.04.09

CHARACTER (LEN=4) :: &
  ylon, ylat       ! names for variables langitude and latitude

CHARACTER (LEN=80)       ::   yzerror

REAL (KIND=ireals)                   :: &
  ! these values are later transferred to the new structures refatm_in, vcoord_in
  p0sl_in, t0sl_in, dt0lp_in, delta_t_in, h_scal_in, vcflat_in, zvcoord_in(kedim+1)

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

#ifdef NETCDF
!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! initialize error status variable
  yerrmsg       = '   '
  ierror        = 0
  istatus       = 0

! Processor 0 does the job
  IF (myid == 0) THEN

!_br 06.04.09
    istatus = nf90_inq_varid(ncid, 'rotated_pole', jgridVarid)
    IF (istatus == NF90_NOERR) THEN
      ylon = 'rlon'
      ylat = 'rlat'
    ELSE
      ylon = 'lon'
      ylat = 'lat'
    ENDIF
!_br 06.04.09 end

    ! Get the spatial dimensions

    istatus = nf90_inq_dimid (ncid, TRIM(ylon) , idims_id(1))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'Error 1 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
      PRINT *, TRIM(yerrmsg)
      RETURN
    ENDIF
    istatus = nf90_inquire_dimension (ncid, idims_id(1), len=iedim_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'Error 2 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
      PRINT *, TRIM(yerrmsg)
      RETURN
    ENDIF

    istatus = nf90_inq_dimid (ncid, TRIM(ylat) , idims_id(2))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'Error 3 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
      PRINT *, TRIM(yerrmsg)
      RETURN
    ENDIF
    istatus = nf90_inquire_dimension (ncid, idims_id(2), len=jedim_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'Error 4 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
      PRINT *, TRIM(yerrmsg)
      RETURN
    ENDIF

    IF (lcm_pres_coor) THEN
      kedim_in  = kedim
      ke1dim_in = kedim + 1
    ELSE
      istatus = nf90_inq_dimid (ncid, 'level', idims_id(3))
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error 5 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF
      istatus = nf90_inquire_dimension (ncid, idims_id(3), len=kedim_in)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error 6 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF

      istatus = nf90_inq_dimid (ncid, 'level1', idims_id(4))
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error 7 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF
      istatus = nf90_inquire_dimension (ncid, idims_id(4), len=ke1dim_in)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error 8 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF
    ENDIF

    ! Get longitude and latitude data

    ALLOCATE (longitude(iedim_in), latitude(jedim_in), vcoordgrib(ke1dim_in), &
              STAT=istatus)
    IF (istatus /= 0) THEN
      yerrmsg = 'Allocation error in read_nc_gdefs_ext_in'
      PRINT *, TRIM(yerrmsg)
      RETURN
    ENDIF

    istatus = nf90_inq_varid (ncid, TRIM(ylon), jlonVarID)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'Error 9 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
      PRINT *, TRIM(yerrmsg)
      RETURN
    ENDIF
    istatus = nf90_get_var (ncid, jlonVarID, longitude)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'Error 10 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
      PRINT *, TRIM(yerrmsg)
      RETURN
    ENDIF

    istatus = nf90_inq_varid (ncid, TRIM(ylat), jlatVarID)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'Error 11 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
      PRINT *, TRIM(yerrmsg)
      RETURN
    ENDIF
    istatus = nf90_get_var (ncid, jlatVarID, latitude)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'Error 12 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
      PRINT *, TRIM(yerrmsg)
      RETURN
    ENDIF

    zdlat_in   = (latitude(jedim_in) - latitude(1))/ REAL(jedim_in-1)
    IF (longitude(iedim_in) - longitude(1) < 0.0) THEN
      ! If the area is located around the 180-Meridian, longitude(iedim_in) - longitude(1)
      ! will be negative and 360 degrees have to be added
      zdlon_in   = (longitude(iedim_in) - longitude(1) + 360.0) / (iedim_in-1)
    ELSE
      zdlon_in   = (longitude(iedim_in) - longitude(1))/ (iedim_in-1)
    ENDIF

    ! get grid mapping values

    istatus = nf90_inq_varid(ncid, 'rotated_pole', jgridVarid)
    IF (istatus == NF90_NOERR) THEN   ! rotated latitude-longitude
      istatus = nf90_get_att(ncid, jgridVarid, 'grid_mapping_name', grid_name_in)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error 13 in read_nc_gdefs_ext_in:  '//TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF
      IF (grid_name_in(1:26) /= 'rotated_latitude_longitude') THEN
        PRINT *, 'Error 14 in read_nc_gdefs_ext_in'
        PRINT *, 'Invalid value for attribute "grid_mapping_name"'
        yerrmsg = TRIM(grid_name_in)
      ENDIF
      istatus = nf90_get_att(ncid, jgridVarid, 'grid_north_pole_latitude',  &
                             zpollat_in)
      IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Error 15 in read_nc_gdefs_ext_in / nf90_get_att'
        PRINT *, 'Attribute "grid_north_pole_latitude"'
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      istatus = nf90_get_att(ncid, jgridVarid, 'grid_north_pole_longitude',  &
                             zpollon_in)
      IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Error 16 in read_nc_gdefs_ext_in / nf90_get_att'
        PRINT *, 'Attribute "grid_north_pole_longitude"'
        yerrmsg = TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      istatus = nf90_get_att(ncid, jgridVarid, 'north_pole_grid_longitude',  &
                             zpolgam_in)
      IF (istatus == NF90_ENOTATT) THEN
        zpolgam_in = 0.0_ireals
      ELSE IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Error 17 in read_nc_gdefs_ext_in / nf90_get_att'
        PRINT *, 'Attribute "north_pole_grid_longitude"'
        PRINT *, TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
    ELSE IF (istatus == NF90_ENOTVAR) THEN  ! true latitude-longitude
      zpollon_in =180.0_ireals
      zpollat_in = 90.0_ireals
      zpolgam_in =  0.0_ireals
    ELSE
       PRINT *, 'Error 18 in read_nc_gdefs_ext_in / nf90_inq_varid'
       yerrmsg = 'Variable "rotated_pole"'
       RETURN
    ENDIF

    ! Get the vertical co-ordinate parameters in case of llm2lm

    IF (llm2lm) THEN
      istatus = nf90_inq_varid (ncid, 'vcoord', jvcVarID)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error in read_nc_gdefs_ext_in: vcoord '// TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF
      istatus = nf90_get_var (ncid, jvcVarID, vcoordgrib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error in read_nc_gdefs_ext_in: '//TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF
      zvcoord_in = REAL(vcoordgrib, ireals)

      istatus = nf90_get_att (ncid, jvcVarID, 'p0sl', p0slgrib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error in read_nc_gdefs_ext_in: p0sl '//TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF
      p0sl_in = REAL (p0slgrib, ireals)

      istatus = nf90_get_att (ncid, jvcVarID, 't0sl', t0slgrib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error in read_nc_gdefs_ext_in: t0sl '//TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF
      t0sl_in = REAL (t0slgrib, ireals)

      istatus = nf90_get_att (ncid, jvcVarID, 'dt0lp', dt0lpgrib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error in read_nc_gdefs_ext_in: dt0lp '//TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF
      dt0lp_in = REAL (dt0lpgrib, ireals)

      istatus = nf90_get_att (ncid, jvcVarID, 'vcflat', vcflatgrib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'Error in read_nc_gdefs_ext_in: vcflat '//TRIM(NF90_strerror(istatus))
        PRINT *, TRIM(yerrmsg)
        RETURN
      ENDIF
      vcflat_in = REAL (vcflatgrib, ireals)

!_br 06.04.09
  istatus = nf90_get_att (ncid, jvcVarID, 'irefatm', irefatm_in)
  IF (istatus == NF90_NOERR) THEN
    istatus = nf90_get_att (ncid, jvcVarID, 'ivctype', ivctype_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'ivctype '//TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF (irefatm_in == 2) THEN ! reference atmosphere type 2
      istatus = nf90_get_att (ncid, jvcVarID, 'delta_t', deltatgrib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'delta_t '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      delta_t_in = REAL(deltatgrib, ireals)
      istatus = nf90_get_att (ncid, jvcVarID, 'h_scal', hscalgrib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'h_scal '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      h_scal_in = REAL(hscalgrib, ireals)
    ENDIF
    IF (ivctype_in == 3 .OR. ivctype_in == 4) THEN ! SLEVE coordinates
      istatus = nf90_get_att (ncid, jvcVarID, 'svc1', svc1grib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'svc1 '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      svc1_in = REAL(svc1grib, ireals)
      istatus = nf90_get_att (ncid, jvcVarID, 'svc2', svc2grib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'svc2 '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      svc1_in = REAL(svc1grib, ireals)
      istatus = nf90_get_att (ncid, jvcVarID, 'nfltvc', nfltvc_in)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'nfltvc '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
    ENDIF
  ELSE ! old version
    irefatm_in = 1
    ivctype_in = 1
  ENDIF
!_br 06.04.09 end
    ENDIF

    ! compare the input values with the Namelist parameters
    IF ((iedim - east_add_in - west_add_in) /= iedim_in)           ierror = 1
    IF ((jedim - south_add_in - north_add_in) /= jedim_in)         ierror = 2
    IF (kedim /= kedim_in)                                         ierror = 3
    IF ((kedim+1) /= ke1dim_in)                                    ierror = 4

    IF (ABS(REAL(zpollat_in,ireals) - pollat) > 1.0E-3_ireals)     ierror = 5
    IF (ABS(REAL(zpollon_in,ireals) - pollon) > 1.0E-3_ireals)     ierror = 6
    IF (ABS(REAL(zpolgam_in,ireals) - polgam) > 1.0E-3_ireals)     ierror = 7

    IF (ABS(REAL(latitude(1),ireals)  - startlat) > 1.0E-3_ireals) ierror = 8
    IF (ABS(REAL(longitude(1),ireals) - startlon) > 1.0E-3_ireals) ierror = 9

    ! In case of lcm2lm, dlat and dlon could be variable. This check makes
    ! no sense then
    IF (.NOT. lcm2lm) THEN
      IF (ABS(REAL(zdlat_in,ireals) - dlat) > 1.0E-3_ireals)       ierror = 10
      IF (ABS(REAL(zdlon_in,ireals) - dlon) > 1.0E-3_ireals)       ierror = 11
    ENDIF

    IF (ierror /= 0) THEN
      PRINT *, '         data file          namelist input'

      PRINT *, 'ie_in_tot             ',iedim_in     ,'       ', iedim-east_add_in-west_add_in
      PRINT *, 'je_in_tot             ',jedim_in     ,'       ', jedim-south_add_in-north_add_in
      PRINT *, 'ke_in_tot             ',kedim_in     ,'       ', kedim
      PRINT *, 'ke_in_tot + 1         ',ke1dim_in    ,'       ', kedim+1

      PRINT *, 'startlat_in_tot       ',latitude(1)  ,'       ', startlat
      PRINT *, 'startlon_in_tot       ',longitude(1) ,'       ', startlon

      IF (.NOT. lcm2lm) THEN
        PRINT *, 'dlat_in               ',zdlat_in   ,'       ', dlat
        PRINT *, 'dlon_in               ',zdlon_in   ,'       ', dlon
      ENDIF

      PRINT *, 'pollat_in             ',zpollat_in   ,'       ', pollat
      PRINT *, 'pollon_in             ',zpollon_in   ,'       ', pollon
      PRINT *, 'polgam_in             ',zpolgam_in   ,'       ', polgam

      istatus = ierror
      RETURN
    ENDIF

    DEALLOCATE (latitude, longitude, STAT=istatus)
    IF (istatus /= 0) THEN
      PRINT *, 'Error 21 in read_nc_gdefs_ext_in / deallocation error'
      yerrmsg = 'Deallocation error in read_nc_gdefs_ext_in'
      istatus = 21
      RETURN
    ENDIF

  ENDIF ! myid == 0

! send information to the other processors
  IF (npes > 1) THEN
    CALL MPI_BCAST (irefatm_in,      1, imp_integers, 0, icomm, izmplcode)
    CALL MPI_BCAST (ivctype_in,      1, imp_integers, 0, icomm, izmplcode)
    CALL MPI_BCAST (p0sl_in,         1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (t0sl_in,         1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (dt0lp_in,        1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (vcflat_in,       1, imp_reals,    0, icomm, izmplcode)
!_br 06.04.09
    IF (irefatm_in == 2) THEN
      CALL MPI_BCAST (delta_t_in,  1, imp_reals,    0, icomm, izmplcode)
      CALL MPI_BCAST (h_scal_in,   1, imp_reals,    0, icomm, izmplcode)
    ENDIF
!_br 06.04.09 end
    CALL MPI_BCAST (zvcoord_in, kedim+1, imp_reals,    0, icomm, izmplcode)
!_br 06.04.09
    IF (ivctype_in == 3 .OR. ivctype_in == 4) THEN  ! SLEVE coordinate
    CALL MPI_BCAST (svc1_in,     1, imp_reals,    0, icomm, izmplcode) 
    CALL MPI_BCAST (svc2_in,     1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (nfltvc_in,   1, imp_integers, 0, icomm, izmplcode)
    ENDIF
!_br 06.04.09 end
  ENDIF

! transfer these values to the new structures refatm_in, vcoord_in
refatm_in%irefatm     = irefatm_in
refatm_in%irefatm_id  = 0
refatm_in%p0sl        = p0sl_in
refatm_in%t0sl        = t0sl_in
IF (irefatm_in == 1) THEN
  refatm_in%dt0lp       = dt0lp_in
  refatm_in%delta_t     = rundefined
  refatm_in%h_scal      = rundefined
ELSE
  refatm_in%dt0lp       = rundefined
  refatm_in%delta_t     = delta_t_in
  refatm_in%h_scal      = h_scal_in
ENDIF
refatm_in%bvref       = rundefined

vcoord_in%ivctype     = ivctype_in
vcoord_in%ivcoord_id  = 0
vcoord_in%nlevels     = kedim+1 
vcoord_in%kflat       = -1
vcoord_in%vc_uuid(:)  = 'x'
vcoord_in%vcflat      = vcflat_in
IF (ivctype_in == 1) THEN
  vcoord_in%sigm_coord(1:kedim+1) = zvcoord_in(1:kedim+1)
  vcoord_in%vert_coord(1:kedim+1) = -1.0_ireals
ELSE
  vcoord_in%sigm_coord(1:kedim+1) = -1.0_ireals
  vcoord_in%vert_coord(1:kedim+1) = zvcoord_in(1:kedim+1)
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------
#endif

END SUBROUTINE read_nc_gdefs_ext_in

!==============================================================================
!==============================================================================
!+ Module procedure in "src_read_ext" for reading a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE read_nc_vdefs_ext_lm (ncid, listin, nvarin, var_id, pollon, pollat, &
                             numlistextpar, ylistextpar, lcheckin_ext,  &
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
    numlistextpar,   & ! number of external fields
    icomm,          & ! MPI communicator
    myid,           & ! ID of this PE in communicator icomm
    npes              ! number of PEs

  REAL    (KIND=ireals),    INTENT(IN)  ::  &
    pollat, pollon        ! coordinates of the rotated north pole

! Array arguments with intent(in):
  TYPE(ar_des_lm)  , INTENT(IN)  ::                      &
    listin(nvarin) ! List of fields for reading in

  CHARACTER (LEN=*)         ::  &
    ylistextpar (numlistextpar)     ! list of external parameters that should be read


! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    istatus       ! error index

  CHARACTER (LEN= *),        INTENT(OUT)    ::  &
    yerrmsg     ! string for error messages

! Array arguments with intent(out):
  INTEGER (kind=iintegers), INTENT(OUT) :: &
    var_id(numlistextpar)       ! NetCDF ID for each variable in the input list

  LOGICAL                    ::  &
    lcheckin_ext(numlistextpar)     ! list for checking which variables has been read


!-----------------------------------------------------------------------
!
! Local scalars:

  INTEGER (KIND=iintegers) ::  &
    jvar_id,          & ! ID for the external variable
    izmplcode           ! error status variable

  INTEGER (KIND=iintegers) ::  &
    j1, j2,           & ! loop index
    idx_lm              ! position of string "_lm"

  CHARACTER (LEN=10) :: &
    varname

  REAL    (KIND=ireals) ::  &
    zpollat, zpollon    ! coordinates of the rotated north pole

  LOGICAL    ::      &
    lzfound             ! true, if variable has been found


! Local arrays:

!- End of header
!------------------------------------------------------------------------------

#ifdef NETCDF
! Processor 0 does the job
  IF (myid == 0) THEN
DO j1 = 1, nvarin
ENDDO
      DO j2 = 1, numlistextpar
ENDDO

    var_loop: DO j1 = 1, nvarin

      idx_lm = index(TRIM(listin(j1)%name),'_lm')
      IF (idx_lm == 0) THEN
        varname = TRIM(listin(j1)%name)
      ELSE
        varname = listin(j1)%name(1:idx_lm-1)
      ENDIF
      lzfound = .FALSE.

      DO j2 = 1, numlistextpar
        IF (TRIM(listin(j1)%name) == TRIM(ylistextpar(j2)) ) THEN
          ! get the variable ID
           istatus = nf90_inq_varid (ncid, TRIM(varname), jvar_id)
           IF (istatus /= NF90_NOERR) THEN
             istatus = 0
             ! Set ID to -1, to indicate missing parameter
             var_id(j2)       = -1
             CYCLE var_loop
           ELSE
             lcheckin_ext(j2) = .TRUE.
             lzfound          = .TRUE.
             var_id(j2)       = jvar_id
             EXIT
           ENDIF
        ENDIF
      ENDDO
      IF (.NOT. lzfound ) CYCLE

      ! get the grid mapping
      istatus = nf90_get_att (ncid, jvar_id, 'grid_mapping', grid_mapping)

      IF (istatus == NF90_enotatt) THEN

        zpollat = 90._ireals
!       zpollon = 180._ireals
!       PIK U. Boehm - 27.11.06
!       zpollon = 0._ireals
        zpollon = 180._ireals
!       PIK U. Boehm - End
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

    ENDDO  var_loop

  ENDIF ! (myid == 0)

! send information to the other processors
  IF (npes > 1) THEN
    CALL MPI_BCAST(var_id, numlistextpar, imp_integers, 0, icomm, izmplcode)
    CALL MPI_BCAST(lcheckin_ext, numlistextpar, imp_logical, 0, icomm, izmplcode)
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
#endif

END SUBROUTINE read_nc_vdefs_ext_lm

!==============================================================================
!==============================================================================
!+ Module procedure in "src_read_ext" for reading a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE read_nc_vdefs_ext_in (ncid, listin, nvarin, var_id, pollon, pollat, &
                                 numlistextpar, ylistextpar, lcheckin_ext,     &
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
    numlistextpar,   & ! number of external fields
    icomm,          & ! MPI communicator
    myid,           & ! ID of this PE in communicator icomm
    npes              ! number of PEs

  REAL    (KIND=ireals),    INTENT(IN)  ::  &
    pollat, pollon        ! coordinates of the rotated north pole

! Array arguments with intent(in):
  TYPE(ar_des_input)  , INTENT(IN)  ::                      &
    listin(nvarin) ! List of fields for reading in

  CHARACTER (LEN=*)         ::  &
    ylistextpar (7)     ! list of external parameters that should be read

! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    istatus       ! error index

  CHARACTER (LEN= *),        INTENT(OUT)    ::  &
    yerrmsg     ! string for error messages

! Array arguments with intent(out):
  INTEGER (kind=iintegers), INTENT(OUT) :: &
    var_id(numlistextpar)       ! NetCDF ID for each variable in the input list

  LOGICAL                    ::  &
    lcheckin_ext(numlistextpar)     ! list for checking which variables has been read

!-----------------------------------------------------------------------
!
! Local scalars:

  INTEGER (KIND=iintegers) ::  &
    jvar_id,          & ! ID for the external variable
    izmplcode           ! error status variable

  INTEGER (KIND=iintegers) ::  &
    j1, j2

  CHARACTER (LEN=10) :: &
    varname

  REAL    (KIND=ireals) ::  &
    zpollat, zpollon    ! coordinates of the rotated north pole

  LOGICAL    ::      &
    lzfound             ! true, if variable has been found

!- End of header
!------------------------------------------------------------------------------

#ifdef NETCDF
! Processor 0 does the job
  IF (myid == 0) THEN

    var_loop: DO j1 = 1, nvarin

      varname = TRIM(listin(j1)%name)
      lzfound = .FALSE.
      DO j2 = 1, numlistextpar
        IF (TRIM(listin(j1)%name) == TRIM(ylistextpar(j2)) ) THEN
          ! get the variable ID
          istatus = nf90_inq_varid (ncid, TRIM(varname), jvar_id)
          IF (istatus /= NF90_NOERR) THEN
            istatus = 0
            ! Set ID to -1, to indicate missing parameter
            var_id(j2)       = -1
            CYCLE var_loop
          ELSE
            lcheckin_ext(j2) = .TRUE.
            lzfound          = .TRUE.
            var_id(j2)       = jvar_id
            EXIT
          ENDIF
        ENDIF
      ENDDO
      IF (.NOT. lzfound ) CYCLE

      ! get the grid mapping
      istatus = nf90_get_att (ncid, jvar_id, 'grid_mapping', grid_mapping)

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

    ENDDO  var_loop

  ENDIF ! myid == 0

! send information to the other processors
  IF (npes > 1) THEN
    CALL MPI_BCAST(var_id, numlistextpar, imp_integers, 0, icomm, izmplcode)
    CALL MPI_BCAST(lcheckin_ext, numlistextpar, imp_logical, 0, icomm, izmplcode)
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
#endif

END SUBROUTINE read_nc_vdefs_ext_in

!==============================================================================

END MODULE src_read_ext
