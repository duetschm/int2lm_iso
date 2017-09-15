!+ Data module for organizational control variables
!==============================================================================

MODULE data_int2lm_control

!==============================================================================
!
! Description:
!  This module contains variables for controlling the run of the interpolation
!  program.
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
! 1.2        2005/07/22 Ulrich Schaettler
!  New Namelist parameter lforest for use of external fields for_e, for_d
! 1.3        2005/12/12 Ulrich Schaettler
!  New Namelist parameter lprog_rho_snow for prognostic treatment of rho_snow
! V1_5         2007/07/09 Ulrich Schaettler
!  Added parameters for treatment of lateral boundaries for qr, qs, qg
!  Added parameters for extra smoothing of topography
!  Added parameters for smooth orography transition at lateral boundaries
!  Added additional parameters for CLM
!  Added debug variable idbg_level
!  Added parameter for treatment of chemistry variables
! V1_6         2007/09/07 Ulrich Schaettler
!  Eliminated logical switch lanalysis
!  Introduced real parameter dt=900
!  Introduced variables for actual data
! V1_7         2007/11/26 Ulrich Schaettler, Christoph Gebhardt
!  Added switch leps_bc to compute boundary conditions for ensemble mode
!  Introduced new namelist parameters for treatment of soil and surface variables
!     itype_rootdp, itype_ndvi, itype_t_cl
!  Renamed switch iw_so_rel_type to itype_w_so_rel to be consistent in the names
! V1_8         2008/05/29 Ulrich Schaettler
!  New Namelist parameter lsso for subgrib scale orography and 
!  lradtopo for topographical corrections for radiation scheme
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
!  Eliminated ldebug
! V1_9         2009/09/03 Ulrich Schaettler, et al.
!  Introduced new NL parameters for reading new external data
!   (lemiss, lstomata, itype_aerosol)
!  Introduced new NL parameters for treatment of humidity (MCH)
!  Eliminated variable lanafg (has been replaced by yinput_type)
! V1_10        2009/12/17 Ulrich Schaettler, Jan-Peter Schulz
!  Introduced new Namelist switch lum2lm to process UKMO data
!  Introduced new Namelist switch llake_coldstart to initialize FLake variables
!  Introduced new Namelist switch lseaice to initialize sea ice variables (JPS)
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_12        2010/06/14 Ulrich Schaettler
!  Updated comments for new types for itype_w_so_rel (from Uwe Boehm, CLM)
! V1_14        2010/11/19 Ulrich Schaettler
!  New NL variable yinput_model (replaces lgme2lm, llm2lm, etc. as NL input)
! V1_17        2011/03/11 Ulrich Schaettler
!  Added lurban flag for reading urban fraction data fr_urban (K. Trusilova)
! V1_19        2012/06/06 Ulrich Schaettler, Burkhardt Rockel, Susanne Brienen
!  Added lhir2lm as internal logical flag
!  Introduction of new NL variable itype_albedo, for choosing different types
!    of albedo
!  Comment line for 365-day year added
!  Splitted msoilgrib (which was used for input and output) to
!    msoilgrib_in (for input) and msoilgrib_lm (for output) (SB)
! V1_20        2012/09/03 Ulrich Schaettler
!  Enlarged strings for date variables to 14 characters
! V1_22        2013/07/11 Ulrich Schaettler, KIT
!  Renamed lprog_qrqs to lprog_qr_qs to be consistent with other names
!  Additional logical variables to indicate INT2LM_ART usage
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters , ONLY :   &
  ireals,    & ! KIND-type parameters for real variables
  iintegers    ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

! 1. Start and end of the run
! ---------------------------

  INTEGER   (KIND=iintegers)       ::           &
    nstart,       & ! start time (in time steps)
    nstop,        & ! end time (in time steps)
    nincbound       ! time step increment between two datasets

! 2. Files for ASCII I/O
! ----------------------

  INTEGER   (KIND=iintegers)       ::           &
    ndebug,       & ! unit number for file with debug information
    nhori,        & ! number of sectors for the horizon array by the 
                    ! topographic correction of the radiation
    ninput,       & ! unit number for NAMELIST input file
    noutput         ! unit number for output file

! 3. Controlling the run
! ----------------------

  LOGICAL                          ::           &
    linitial,     & ! if .TRUE., initial data for LM
    lboundaries,  & ! if .TRUE., lateral boundaries for LM
    lcomp_bound,  & ! compute fields for boundaries
    lgme2lm,      & ! if .TRUE., gme->lm
    lgfs2lm,      & ! if .TRUE., gfs->lm
    lgsm2lm,      & ! if .TRUE., gsm->lm
    lec2lm,       & ! if .TRUE., ec ->lm
    llm2lm,       & ! if .TRUE., lm ->lm
    lum2lm,       & ! if .TRUE., um ->lm
    lhir2lm,      & ! if .TRUE., hirlam ->lm
    lhm2lm,       & ! if .TRUE., hm ->lm
    lcm2lm,       & ! if .TRUE., climate model ->lm   !_br
    l_smi,        & ! if .TRUE., interpolate soil moisture with SMI 
                    ! (Soil Moisture Index; by MCH)
    luse_t_skin,  & ! if .TRUE., use ECMWF skin temperature for surface
    lante_0006,   & ! if .TRUE., force to use ECMWF dataset before 27 June 2000
    lpost_0006,   & ! if .TRUE., force to use ECMWF dataset after 27 June 2000
    luvcor,       & ! if .TRUE., correct winds for given surface pressure  
                    ! tendency
    lvertwind_ini,& ! if .TRUE., compute vertical wind for LM for initial data
    lvertwind_bd, & ! if .TRUE., compute vertical wind for LM for boundary data
    lprog_qi,     & ! if .TRUE., interpolate qi from GME to LM grid
    lmixcld,      & ! if .TRUE., qi added in grh instead of being interpolated straightforward
    lprog_qr_qs,  & ! if .TRUE., interpolate qr,qs from LM to LM grid
    lprog_qg,     & ! if .TRUE., interpolate qg from LM to LM grid
    lprog_rho_snow,&! if .TRUE., interpolate rho_snow from GME to LM grid
    lfilter_oro,  & ! if .TRUE., filter the orography
    lxso_first,   & ! if .TRUE., do eXtra smoothing of orography first
    lfilter_pp,   & ! if .TRUE., filter the pressure deviation after vertical
                    !            interpolation
    lbalance_pp,  & ! if .TRUE., compute a hydrostatic balanced pp after
                    !            vertical interpolation in LM2LM
    llbc_smooth     ! if .TRUE., run with smooth orography transition at LB

  LOGICAL                          ::           &
    lt_cl_corr,   & ! if .TRUE., a height-correction of T_CL is performed
    lbdclim,      & ! if .TRUE., special boundary data for climate mode
    lforest,      & ! if .TRUE., run with forest (evergreen and deciduous)
    lurban,       & ! if .TRUE., use urban fraction data
    lsso,         & ! process parameters for sso scheme
    lradtopo,     & ! process parameters for topographic correction of radiation
    lseaice,      & ! if .TRUE., run with sea ice model
    llake,        & ! if .TRUE., run with lake parameters    !_br
    llake_coldstart,& ! if .TRUE., initialize prognostic lake variables for cold start
    lemiss,       & ! if .TRUE., run with external parameter for surface emissivity
    lstomata,     & ! if .TRUE., run with external parameter for stomata resistance
    lroutine,     & ! if .TRUE., routine-job
    lclock,       & ! if .TRUE., system clock is present
    ltime,        & ! detailed timings of the program are given
    ltime_mean,   & ! if .TRUE., mean values of the timings are printed
    ltime_proc,   & ! if .TRUE., timings for each processor are printed
! iso code
    liso            ! if .TRUE., include variables for water isotope simulation
! end iso code

  REAL (KIND=ireals)               ::           &
    qvmin,        & ! minimum value of water vapor (security)
    qcmin,        & ! minimum value of cloud water (security)
    qimin           ! minimum value of cloud ice content (security)

  CHARACTER (LEN=5)                ::           &
    yinput_model    ! string to identify the input model
                    ! valid options: 'GME', 'COSMO', 'IFS', 'UM', 'GSM', 'GFS'

! 4. Controlling print/debug information
! --------------------------------------

  LOGICAL                          ::           &
    lprps,        & ! logical for print of different ps* und fis - fields
    lprt,         & ! logical for print of T at 2 levels (nlev1pr,nlev2)
    lpru,         & ! same as lprt but for U
    lprv,         & ! same as lprt but for V
    lprgrh,       & ! same as lprt but for General Relative Humidity
    lprqv,        & ! same as lprt but for QV
    lprqc,        & ! same as lprt but for QC
    lprqi,        & ! same as lprt but for QI
    lprud,        & ! same as lprt but for UD (divergent wind correction)
    lprvd,        & ! same as lprt but for VD (divergent wind correction)
    lprdpdt         ! same as lprt but for DPDT (tendency of surface pressure)

  LOGICAL                          ::           &
    lprintdeb_all   ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output

  INTEGER   (KIND=iintegers)       ::           &
    idbg_level,   & ! to control the verbosity of debug output
    nlev1pr,      & ! k-index for the print of first  model layer
    nlev2pr         ! k-index for the print of second model layer

! 5. Additional control variables
! -------------------------------

  LOGICAL                          ::           &
    lbd_frame,    & ! if .TRUE., boundary fields include only frames
    lbd_frame_cur   ! if .TRUE., current boundary fields include only frames

  INTEGER   (KIND=iintegers)       ::           &
    nvers,        & ! for documenting purposes in Grib-Code
    nl_soil_in,   & ! number of soil layers in input model
    nl_soil_lm,   & ! number of soil layers in LM, resp. HM
    norder_filter,& ! order of the orography filtering
    ilow_pass_oro,& ! type of low-pass filter for orography
    numfilt_oro,  & ! number of sequential applications of filter
    ilow_pass_xso,& ! type of low-pass filter for extra smoothing of steep oro.
    numfilt_xso,  & ! number of sequential applications of xso-filter
    ifill_valley, & ! type of valley filling: 1=MIN (before), 2=MIN (after)
    nlbc_smooth,  & ! number of grid points for smooth orography transition at LB
    kcontrol_fi,  & ! control level for geopotential
    npstrframe      ! thickness of output frames

  REAL (KIND=ireals), PARAMETER    ::           &
    dt=900.0_ireals ! time step used in the LM

  REAL (KIND=ireals)               ::           &
    eps_filter,   & ! parameter for orography filtering
    rxso_mask,    & ! mask for extra smoothing of steep oro.: dh > rxso_mask
    rfill_valley, & ! mask for valley filling: dh > rfill_valley
    pcontrol_fi     ! pressure of control level for geopotential

  LOGICAL                          ::           &
    l_cressman,  &  ! logical switch for controling the use of a Cressmann 
                    ! scheme for horizontal interpolation in case of a 'Match'
                    ! type interpolation when no bilinear interpolation is 
                    ! possible
    l_bicub_spl     ! logical switch for controling the use of a bicubic
                    ! spline interpolation (16-points formula)

! 6. Variables for date and timings
! ---------------------------------

  CHARACTER (LEN=14)               ::           &
    yakdat1   ! actual date (ydate_ini+ntstep/dt) in the form
              ! yyyymmddhhmmss (year, month, day, hour, minutes, seconds)
  CHARACTER (LEN=28)               ::           &
    yakdat2   ! actual date (ydate_ini+ntstep/dt) in the form
              ! wd dd.mm.yyyy  hh mm ss UTC  (weekday, ...)

  REAL (KIND=ireals), ALLOCATABLE  ::           &
    timings (:)   ! for storing the times for different parts of the program

  INTEGER   (KIND=iintegers)       ::           &
    itype_calendar   ! for specifying the calendar used
                     !  = 0: gregorian calendar (default)
                     !    (but this needs a bug fix in get_utc_date,
                     !    because up to now we only have julian calendar)
                     !  = 1: every year has 360 days
                     !  = 2: every year has 365 days

! 7. Variables for soil and surface variables
! -------------------------------------------

  INTEGER   (KIND=iintegers), ALLOCATABLE       ::           &
    msoilgrib_lm(:), & ! grib coded depth of main soil levels in centimeters for
                       ! multi-layer soil model in the outgoing data
    msoilgrib_in(:)    ! grib coded depth of main soil levels in centimeters for
                       ! multi-layer soil model in the incoming data
                       ! (careful: the first level will be coded with 1,
                       !           but is in the depth of 0.5 cm!)

  LOGICAL                          ::           &
    lmulti_layer_lm, & ! if .TRUE., compute soil fields for multi-layer soil
                       ! model in the outgoing data
    lmulti_layer_in    ! if .TRUE., incoming data have soil fields from the
                       ! multi-layer soil model

  INTEGER   (KIND=iintegers)       ::           &
    itype_w_so_rel     ! type of relative soil moisture input
                       ! 0 = artificial profile relative to pore volume
                       ! 1 = relative to pore volume (read from coarse grid data)
                       ! 2 = relative to field capacity (read from coarse grid data)
                       ! 3 = as 1, but constant relative soil moisture below deepest input soil layer
                       ! 4 = as 2, but constant relative soil moisture below deepest input soil layer

  INTEGER   (KIND=iintegers)       ::           &
    itype_t_cl,      & ! to choose origin and treatment of deep soil temperature
                       !   0: take t_cl from driving coarse grid model
                       !   1: take t_cl from external parameters
                       !
    itype_rootdp,    & ! to choose treatment of root depth
                       !   0: take input from external parameter and modify
                       !      with a yearly cycle
                       !   1: just take the input from the external parameters
                       !   2: like 0, but with adapting to ECOCLIMAP niveau
                       !
    itype_ndvi,      & ! to choose treatment of surface parameters (plcov, lai)
                       !   0: no monthly ndvi rates are taken, plcov and lai are
                       !      computed from min- and max-values by a yearly cycle
                       !   1: take monthly ndvi rates and max-values from plcov
                       !      and lai to get an actual value for every day
                       !      additional external parameter are read for that
                       !   2: take monthly climatological values for plcov, lai
                       !      (and also z0) and interpolate to a special day
                       !
    itype_aerosol,   & ! to choose treatment of surface parameters (plcov, lai)
                       !   0: not implemented yet
                       !   1: take default aerosol background concentration
                       !   2: take external parameters for different aerosol types
    itype_albedo       ! type of surface albedo treatment
                       ! 1: surface albedo is a function of soiltype
                       !    (default: no actions needed in INT2LM)
                       ! 2: surface albedo is prescribed by external fields
                       !    for dry and saturated soil
                       ! 3: background albedo is prescribed by external fields
                       ! 4: vegetation albedo is modified by forest fraction
                       !    (done in the COSMO-Model, no actions needed in INT2LM)

! 8. Variables for the LM-Z coordinate Version
! --------------------------------------------

  LOGICAL                          ::           &
    l_topo_z           ! if .TRUE., additional smoothing of the topography
                       ! for the LM-Z coordinate Version

! 9. controlling ensemble mode (EPS)
! -----------------------------------

  LOGICAL                          ::           &
    leps_bc            ! switch ensemble mode on/off for boundary data

! 10. Variables for Chemistry
! ---------------------------

  LOGICAL                          ::           &
    l_art_nested,                               &
    l_art

!==============================================================================

END MODULE data_int2lm_control
