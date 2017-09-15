!+ Data module for variables concerned with I/O
!==============================================================================

MODULE data_int2lm_io

!==============================================================================
!
! Description:
!  This data module contains all data necessary for input and output of grib
!  files.
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
!  Increased nmaxlist to 50
!  Added nvar_lm_chem, nvar_in_chem; nvar_lm_norm, nvar_in_norm
!      set nvar_in, nvar_lm according to l_chemistry
!  Additional variables for NetCDF IO
! V1_6         2007/09/07 Ulrich Schaettler
!  Splitted ytunitbd to ytunit_in, ytunit_out
!  Introduced character variable for type of input data: yinput_type
! V1_7         2007/11/26 Ulrich Schaettler
!  Increased nvar_in to 60
! V1_8         2008/05/29 Ulrich Schaettler
!  Increased nvar_lm_norm to 70
! V1_9         2009/09/03 Guenther Zaengl
!  Replaced ldwd_grib_use by l_ke_in_gds
!  Increased nvar_lm_norm to 90
! V1_10        2009/12/17 Ulrich Schaettler
!  Increased nvar_lm_norm to 100
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_14        2010/11/19 Ulrich Schaettler, Anne Roches
!  Introduce idwdednr (old iednr) and igrbednr (for input with grib_api)
!  Modified type definition of ar_des_input for working with grib_api
!  Introduced additional NetCDF dimension for external parameter for topographical
!     correction
! V1_19        2012/06/06 Burkhardt Rockel
!  Unified dimension IDs with COSMO (ID for topo corrections changed from 14 to 15)
! V1_20        2012/09/03 Ulrich Schaettler, Burkhardt Rockel, Davide Cesari
!  Enlarged strings for date variables to 14 characters
!  Introduced (internal) variable lmmss to indicate whether the 14 digits or the
!    10 digits format is used
!  Introduction of additional global attributes in case of netCDF output
!    which can be set via namelist (analog the attributes in the namelist in COSMO)
!    (Burkhardt Rockel)
!  Increased nmaxlist to 60 to be able to allow all possible new interpolation features
!    (Davide Cesari)
! V1_21        2013/03/25 Ulrich Schaettler, Burkhardt Rockel
!  Introduced global variables for the handles of grib sample data (US)
!  Added string for level type in variable table var_lm for COSMO-Model (US)
!  Added new NL variables for GRIB2 handling: nsubcenter, nlocaldefnr (US)
!  Changes in netCDF output: scalar variable instead of extra dimension of length 1
!    (e.g. height_2m) needs less definitions of idims_id -> ndims_id is changed. (BR)
! V1_22        2013/07/11 Ulrich Schaettler, KIT
!  Implemented arrays for vertical coordinate parameters (pv_in, pv_out)
!  Introduced new organizational variables for Grib handling:
!    ylevltypes_in, ylevltypes_out, ysteptypes, rscalefac
!  New Namelist variables ylm_hhl, yin_hhl for HHL file names
!  New list variables youtlist_hhl, numlist_hhl for writing HHL
!  Renamed grib buffers: ds_grib to ds_grib_single, ds_gribapi to ds_grib_double
!  Increased number of variables in variable table for INT2LM_ART (KIT)
!  Introduced a possible extension to files: yextension (KIT)
! V1_23        2013/10/02 Ulrich Schaettler
!  Define pv_in as double precision, because it is treated as such in the INT2LM
!
!  Increased nvar_lm_norm to 110 (Stephan Pfahl)
!  Increased nvar_in_norm to 76 (Hui Tang 2013-11-20)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters, ONLY : &
  ireals,    & ! KIND-type parameters for real variables
  iintegers, & ! KIND-type parameter for standard integer variables
  irealgrib, & ! kind-type parameter for the real variables in the grib library
  intgribf     ! kind-type parameter for the fortran files of the grib library

!==============================================================================

IMPLICIT NONE

!==============================================================================

! 1.1. Global parameters, dimensions and arrays for the GRIB-library
! ------------------------------------------------------------------

  INTEGER (KIND=intgribf)             :: &
    idwdednr  =  1,    & ! grib edition number for DWD library
    igrbednr             ! grib edition number for grib_api (to be set)

  INTEGER (KIND=intgribf),  PARAMETER :: &
    npds   =    321_intgribf, & ! Dimension for product definition section(pds)
    ngds   =    626_intgribf, & ! Dimension for grid description section (gds)
    nbms   =      3_intgribf, & ! Dimension for bit map section (bms)
    nbds   =     11_intgribf, & ! Dimension for binary data section
    ndsup  =     73_intgribf, & ! Dimension for dsup
    ndims  =     20_intgribf    ! Dimension for idims (contains all dimensions)

  ! Grib handles for the grib_api sample data
  INTEGER (KIND=intgribf)             :: &
    igrib1_sample,        igrib2_sample,          & ! common samples
    igrib1_hybrid,        igrib2_hybrid,          & ! full level data (U, V, T, etc)
    igrib1_hybridlayer,   igrib2_hybridlayer,     & ! half level data (W, HHL)
    igrib1_surface,       igrib2_surface,         & ! surface data
    igrib1_depthbelow,    igrib2_depthbelow         ! soil data

  ! The following dimensions are set during program execution
  INTEGER (KIND=intgribf)             :: &
    lfd,                      & ! Dimension for iblock
    lfa,                      & ! Dimension for grib_api message in bytes
    nbitmap,                  & ! Dimension for bitmap: this value has to be at 
                                ! least the same as lbmax in subroutine grbin1
    lds,                      & ! Dimension for unpacked data
    inrvert_in,               & ! number of vertical coordinate parameters of input data
    inrvert_out                 ! number of vertical coordinate parameters of output data

  INTEGER (KIND=iintegers)            :: &
    nvar_lm,         & ! maximum number of variables in LM variable table
    nvar_in            ! maximum number of variables in input variable table

  INTEGER (KIND=iintegers), PARAMETER :: &
! iso code: increase value
    nvar_lm_norm = 110, & ! maximum number of variables in LM variable table
    nvar_in_norm =  76, & ! maximum number of variables in input variable table
                          ! without chemistry fields
    nvar_lm_chem = 650, & ! maximum number of variables in LM variable table
    nvar_in_chem = 310    ! maximum number of variables in input variable table
                          ! using chemistry fields

  REAL    (KIND=irealgrib), PARAMETER :: &
    undefgrib =  -1.E7_irealgrib, & ! value for "undefined" in the grib routines
    undefncdf =  -1.E20_irealgrib
 
  REAL    (KIND=ireals)               :: &
    undef             ! the same as undefgrib but with other KIND-Parameter

  INTEGER (KIND=iintegers)            :: &
    nrbit             ! packrate

  INTEGER (KIND=intgribf)              :: &
    idims_in  (ndims),  & ! array for all dimensions (input)
    idims_out (ndims),  & ! array for all dimensions (output)
    ipds      (npds),   & ! product definition section
    igds_in   (ngds),   & ! grid description section for input
    igds_out  (ngds),   & ! grid description section for output
    ibms      (nbms),   & ! bit map section
    ibds      (nbds)      ! binary data section

  ! packed grib_field -> grib_api
  CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: ymessage

  ! Arrays that are allocated during program execution
  INTEGER (KIND=intgribf), ALLOCATABLE :: &
    iblock    (:),      & ! array for gribed data
    ibmap     (:)         ! array for bit map

  REAL   (KIND=irealgrib), ALLOCATABLE :: &
    dsup   (:),       & ! array for special data
    ds_grib_single(:)   ! array for unpacked data

  REAL   (KIND=ireals),    ALLOCATABLE :: &
    pv_in  (:),       & ! array for vertical coordinate parameters for input data
    pv_out (:),       & ! array for vertical coordinate parameters for output data
    ds_grib_double(:)   ! array for unpacked data

! Arrays to convert GRIB1 level types and time range indicators to grib_api strings
! for GRIB1 (ylevltypes(:,1)) and GRIB2 (ylevltypes(:,2)), resp.
CHARACTER(LEN=30)  :: ylevltypes_in(0:255,2), ylevltypes_out(0:255,2)
CHARACTER(LEN=10)  :: ysteptypes(0:10)

! Array to convert GRIB2 scale factors to real numbers
REAL (KIND=ireals) :: rscalefac(0:9)

!------------------------------------------------------------------------------

! 1.2. Global parameters, dimensions and arrays for the NetCDF-library
! --------------------------------------------------------------------

  INTEGER (KIND=iintegers), PARAMETER :: &
    ndims_id=14  ! array for the IDs of the dimensions of netCDF formatted output
                 ! The different dimensions are:
                 !   1: ID for ie,                                 "rlon"
                 !   2: ID for je,                                 "rlat"
                 !   3: ID for ke,                                 "level" ("pressure", "altitude")
                 !   4: ID for ke1,                                "level1"
                 !   5: ID for ntime,                              "time"
                 !   6: ID for nbnds,                              "bnds"
                 !   7: ID for ke_soil,                            "soil"
                 !   8: ID for ke_soil+1,                          "soil1"
                 !   9: ID for staggered ie,                       "srlon"
                 !  10: ID for staggered je,                       "srlat"
                 !  11: ID for sections of topo ext. parameter     "nhori"
                 !  12: ID for all products of syn. sat. data      "nsynmsg"
                 !  13: ID for one group of products of MSG        "msgchan"
                 !  14: ID for ke_snow                             "snow_layer"

  INTEGER (KIND=iintegers)   :: &
    idims_id (ndims_id) ! array for the IDs of the dimensions

! NetCDF global attributes
CHARACTER (LEN=100)  :: &  ! length should be <= length of charbuf
  yncglob_institution,   & ! originating center name
  yncglob_title,         & ! title string for the output
  yncglob_source,        & ! program name and version
  yncglob_project_id,    & ! identification of the project of simulation
  yncglob_experiment_id, & ! identification of the experiment of simulation
  yncglob_contact,       & ! contact e.g. email address
  yncglob_references       ! URL, report etc.

INTEGER   (KIND=iintegers)          :: &
  ncglob_realization       ! number of the realization of the experiment

!------------------------------------------------------------------------------

! 2. Variables for handling the Gribfile and NetCDF I/O:
! ------------------------------------------------------

  INTEGER (KIND=iintegers)            :: &
    nlocaldefnr,    & ! local definition number for GRIB local section (Namelist parameter)
    nactlocdefnr,   & ! to overwrite Namelist parameter with some center default
    nbitmappres,    & ! to check presence of bitmap
    nunit_of_time,  & ! indicator for unit-of-time (1hr, 15min, 30min,...)
    nprocess_ini,   & ! type of database for initial LM data (analysis)
    nprocess_bd,    & ! type of database for boundary LM data (forecasts)
    ncenter,        & ! originating center identification
    nsubcenter        ! originating subcenter identification

  INTEGER   (KIND=iintegers)       ::           &
    nincwait,     & ! if ready-file is not available wait nincwait seconds
                    ! until next attempt
    nmaxwait        ! if ready-file is not available after nmaxwait seconds,
                    ! abort the program

  LOGICAL                          ::           &
    l_ke_in_gds , & ! explicit GDS entry for number of model levels
    lchkin,       & ! logical for print of check-values (max,min,mean)
                    ! of GME-fields
    lchkout,      & ! logical for print of check-values (max,min,mean)
                    ! of LM/HM-fields
    lcheck_bmm,   & ! if .TRUE., check bitmaps for mass grid points
    lcheck_bmu,   & ! if .TRUE., check bitmaps for u-grid points
    lcheck_bmv      ! if .TRUE., check bitmaps for v-grid points

  CHARACTER (LEN=250)              ::           &
    ytrans_in,    & ! directory for reading ready-files
    ytrans_out      ! directory for writing ready-files

  INTEGER (KIND=iintegers), PARAMETER       ::           &
    nmaxlist=60      ! maximal number of elements in output lists

  INTEGER (KIND=iintegers)         ::           &
    numlist_ini,   & ! number of elements in youtlist_ini
    numlist_bd,    & ! number of elements in youtlist_bd
    numlist_hhl      ! number of elements in youtlist_hhl

  CHARACTER (LEN=10)               ::           &
    youtlist_ini (nmaxlist),   &  ! list of output variables for initial data
    youtlist_bd  (nmaxlist),   &  ! list of output variables for boundary data
    youtlist_hhl (nmaxlist)       ! list of output variables for boundary data

  CHARACTER (LEN=250) ::    &
    ylmext_cat,  & ! catalog of external LM parameters
    yinext_cat,  & ! catalog of external parameters from input model
    ylm_cat,     & ! catalog of the files where LM fields are written
    yin_cat        ! catalog of the files where input fields are read

  CHARACTER (LEN= 50) ::    &
    ylmext_lfn,  & ! name of the file with external LM parameters
    yinext_lfn,  & ! name of the file with external GME parameters
    ylm_lfn,     & ! name of the file with LM output data
    yin_lfn,     & ! name of the file with GME input data
    ylm_hhl,     & ! name of the file with output (GRIB2) HHL fields
    yin_hhl        ! name of the file with input (GRIB2) HHL fields
 
  CHARACTER (LEN=4)                ::           &
    ylmext_form_read,   & ! input format of external LM data
    yinext_form_read,   & ! input format of external boundary data
    yin_form_read,      & ! input format of boundary data
    ylm_form_write        ! output format of LM data

  CHARACTER (LEN=8)                ::           &
    yinput_type           ! type of input data: 'forecast', 'analysis' or 'ana_init'

  CHARACTER (LEN= 1)  :: ytunit_in   ! time unit for input data
  CHARACTER (LEN= 1)  :: ytunit_out  ! time unit for output data
                                     ! only for boundaries
  CHARACTER (LEN=14)  :: ydate_ini   ! start of the forecast
                                     ! yyyymmddhh (year, month, day, hour)
  CHARACTER (LEN=14)  :: ydate_bd    ! start of the forecast from which the
                                     ! boundary fields are used

  LOGICAL                     ::    &
    lmmss_ini,  & ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                  ! if .FALSE. 10 digits date format (YYYYMMDDHH)
                  ! for ydate_ini and result files of INT2LM
    lmmss_bd      ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                  ! if .FALSE. 10 digits date format (YYYYMMDDHH)
                  ! for ydate_bd  and input  files of INT2LM

  CHARACTER (LEN=  8) :: yuchkdat    ! checking the I/O data
  CHARACTER (LEN=  3) :: ymode_read  ! mode for opening the (read) Grib files
  CHARACTER (LEN=  3) :: ymode_write ! mode for opening the (write) Grib files

  INTEGER (KIND=iintegers)     ::    &
    nuchkdat,   & ! logical UNIT-number for checking the I/O data
    njulianday    ! julian day in the year (returned by get_utc_date)

  REAL (KIND=ireals)           ::    &
    ract_hour     ! actual hour of actual day (returned by get_utc_date)

  CHARACTER(LEN=1)    ::  yextension = ' '! file extension necessary for int2lm-art

! 3. Defined data types for I/O
! -----------------------------

  ! Type for input model variable table
  TYPE ar_des_input
    CHARACTER (LEN=10)             :: name           ! name of variable
    CHARACTER (LEN=30)             :: ylevtyp        ! name of level type
    INTEGER(KIND=intgribf)         :: tabtyp         ! grib table number
    INTEGER(KIND=intgribf)         :: levtyp         ! code for leveltype
    INTEGER(KIND=intgribf)         :: ee             ! element code
    INTEGER(KIND=intgribf)         :: levtop         ! top of layer
    INTEGER(KIND=intgribf)         :: levbot         ! bottom of layer
    REAL(KIND=ireals)              :: factor         ! factor
    REAL(KIND=ireals)              :: bias           ! bias
    INTEGER(KIND=intgribf)         :: rank           ! number of dimensions
    REAL(KIND=ireals), POINTER     :: p4(:,:,:,:)    ! pointer for rank4 var
    REAL(KIND=ireals), POINTER     :: p3(:,:,:)      ! pointer for rank3 var
    REAL(KIND=ireals), POINTER     :: p2(:,:)        ! pointer for rank2 var
    CHARACTER (LEN=3)              :: dattyp         ! kind of data
    CHARACTER (LEN=3)              :: ipc            ! Interpolation Code:
    LOGICAL                        :: lreadin        ! variable has been read
    INTEGER(KIND=iintegers)        :: nlevels        ! number of levels
    INTEGER(KIND=iintegers)        :: nlevels_read   ! number of levels read
  END TYPE ar_des_input

  ! Type for LM variable table
  TYPE ar_des_lm
    CHARACTER (LEN=10)             :: name           ! name of variable
    CHARACTER (LEN=30)             :: ylevtyp        ! name of level type
    INTEGER(KIND=intgribf)         :: tabtyp         ! grib table number
    INTEGER(KIND=intgribf)         :: levtyp         ! code for leveltype
    INTEGER(KIND=intgribf)         :: ee             ! element code
    INTEGER(KIND=intgribf)         :: levtop         ! top of layer
    INTEGER(KIND=intgribf)         :: levbot         ! bottom of layer
    REAL(KIND=ireals)              :: factor         ! factor
    REAL(KIND=ireals)              :: bias           ! bias
    INTEGER(KIND=intgribf)         :: rank           ! number of dimensions
    REAL(KIND=ireals), POINTER     :: p3(:,:,:)      ! pointer for rank3 var
    REAL(KIND=ireals), POINTER     :: p2(:,:)        ! pointer for rank2 var
    CHARACTER (LEN=20)             :: units          ! unit of variable
    CHARACTER (LEN=80)             :: standard_name  ! standard name of variable
    CHARACTER (LEN=80)             :: long_name      ! description of variable
    CHARACTER (LEN=1)              :: lsm            ! land/sea mask flag
  END TYPE ar_des_lm

  TYPE (ar_des_lm),    ALLOCATABLE   :: var_lm (:)
  TYPE (ar_des_input), ALLOCATABLE   :: var_in (:)

!===============================================================================

END MODULE data_int2lm_io
