!+ Source module for the namelist input of the interpolation program
!==============================================================================

MODULE src_namelists

!==============================================================================
!
! Description:
!   This module performs the input of the namelist groups for the 
!   interpolation program.
!   - /CONTRL/
!   - /GRID_IN/
!   - /LMGRID/
!   - /DBASE/
!   - /DATA/
!   - /PRICTR/
!   - /EPSCTL/
!   the input is organized by the routine "read_namelists"
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
!  Adapt list of output fields accordingly
!  Check that kcontrol_fi (only for GME) and pcontrol_fi (only for IFS) are
!  not both set.
! 1.3        2005/12/12 Ulrich Schaettler
!  New Namelist parameter lprog_rho_snow for prognostic treatment of 
!  snow density
! V1_5         2007/07/09 Ulrich Schaettler
!  Introduced new namelist variables for 
!   - lateral boundary treatment for qr, qs, qg
!   - extra smoothing of topography
!   - smooth transition of topography at lateral boundaries
!   - treatment for CLM: several
!   - treatment of I/O (grib, netcdf)
!   - to  skip levels at the top of ECMWF model (Davide Cesari)
!   - treatment of chemistry variables
!  Limited values of endlon_in_tot to less than 180.0
!  Adjusted lengths of output strings
! V1_6         2007/09/07 Burkhardt Rockel, Ulrich Schaettler
!  Eliminated logical switch lanalysis
!  Introduced yinput_type
!  If lbdclim is set, change default of yvarbd
! V1_7         2007/11/26 Christoph Gebhardt, Ulrich Schaettler
!  Renamed iw_so_rel_type to itype_w_so_rel to be consistent with other names
!  Introduced new namelist group EPSCTL for ensemble mode
!  Introduced new namelist parameters for treatment of soil and surface variables
!     itype_rootdp, itype_ndvi, itype_t_cl
!  Bug correction for setting nstart, nstop, nincbound
! V1_8         2008/05/29 Ulrich Schaettler, Hans-Juergen Panitz
!  New Namelist parameters lsso and lradtopo
!  Bug fix: lcm2lm has been forgotten for the distribution
!  Allow max value of 3 for itype_rootdp (another option added)
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
!  Enlarged size of some character strings
!  Included 32-layer version again
! V1_9         2009/09/03 Guenther Zaengl, Ulrich Schaettler, et al.
!  Added namelist parameters for new reference atmosphere
!  Allow bicubic spline interpolation only for lcm2lm=.TRUE. (Uli)
!  Rename ldwd_grib_use to l_ke_in_gds
!  New NL parameters for additional external parameters: 
!    lstomata, lemiss, itype_aerosol
!  New NL parameter for treatment of soil properties: l_smi  (Guy deMorsier)
!  New NL parameter to add qi to generalized relative humidity (MCH)
!  New NL parameter to specify  number of sectors for the horizon array (MCH)
!  Increase itype_ndvi   to 2 (Daniel Luethi)
!  Increase itype_rootdp to 4 (Burkhardt Rockel)
!  Eliminated lanafg
! V1_10        2009/12/17 Oliver Fuhrer, Heike Vogel, Ulrich Schaettler
!  Call to distribute_values only if number of processors is greater than 1
!  delete (NOT lec2lm) for lmulti_layer_in=FALSE to get variables for old soil model
!  Read new Namelist switches lum2lm, llake_coldstart, lseaice
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_12        2010/06/14 Oliver Fuhrer
!  Some technical adaptations
! V1_14        2010/11/19 Ulrich Schaettler, Burkhardt Rockel
!  Read new NL switch yinput_model and eliminate lgme2lm, lec2lm etc. as NL switches
!  Read new NL variable press_level for available pressure levels
!  Added field for salt_lk_lm (BR)
! V1_15        2010/12/10 Ulrich Schaettler
!  Increased default value for nbitmap_d for higher resolution GME
! V1_17        2011/03/11 Ulrich Schaettler
!  Add a namelist variable lurban for reading urban fraction data from the
!  external parameter dataset (K. Trusilova)
! V1_19        2012/06/06 Ulrich Schaettler, Dmitrii Mironov, Burkhardt Rockel
!  Implemented conditional compilation for NETCDF
!  Check consistency between lmixcld and lprog_qi
!  Split long directory names for distribution using charbuf (len=100)
!  Distribute values of czmls_in for Namelist distribution instead of czml_soil_in
!    (by Susanne Brienen)
!  Set internal logical flag lhir2lm, depending on yinput_model
!  Allow UMR (for regional) and UMG (for global) model
!  Allocate t_ice_lm, h_ice_lm also for lake cold starts (DM)
!  Use of msoilgrib_in, msoilgrib_lm instead of msoilgrib
!  CLM:
!    Added 365_days support
!    Added new NL switch itype_albedo, to choose surface albedo treatment
!    Deactivated l_cressman until it is tested
! V1_20        2012/09/03 Michael Baldauf, Burkhardt Rockel
!  Implemented namelist switch lanalyt_calc_t0p0 in Namelist group /LMGRID/
!   (only for irefatm=1)
!     =.FALSE.: current method (is default)
!     =.TRUE. : compute t0, p0 on full levels analytically
!               (needed for new fast waves solver)
!  Implemented 14-digit format for date variables;
!  Initialize lmmss_[ini,bd] according to COSMO-Model variable lmmss
!  Renamed 'grax' to 'apix' to be conform with other models (US)
!  Burkhardt Rockel:
!  Introduction of additional global attributes in case of netCDF output
!    which can be set via namelist (analog the attributes in the namelist in COSMO)
!  Change print statements to account for new namelist parameter yinput_model
!  Changes to write T_S for both climate and forecast version in case of netCDF output
! V1_21        2013/03/25 Ulrich Schaettler
!  Added new NL variables for GRIB2 handling: nsubcenter, nlocaldefnr
!  Set and check proper relation between lroutine and nvers for Grib2 meta data handling
!  Set ydate_bd in all cases with 14 digits
! V1_22        2013/07/11 Ulrich Schaettler, Davide Cesari, KIT
!  Put Salinity to list of initial fields only in case of climate runs
!  Renamed lprog_qrqs to lprog_qr_qs to be consistent with other names
!  Renamed PRS_MIN to RSMIN, which is the official shortName
!  Removed p0sl_in, t0sl_in, dt0lp_in from Namelist GRIDIN
!  Grouped reference atmosphere parameters to new type refatm_out
!  Grouped vertical coordinate  parameters to new type vcoord_out
!     and removed some (old) defaults for vcoord
!  New Namelist variable :  lnewVGrid for GRIB2 general vertical coordinate
!  New Namelist variables: ylm_hhl, yin_hhl with names of HHL-files for in- and output
!  Allow new namelist settings, allocate and compute multi-layer depths for IFS
!  Replaced namelist variable l_chemistry by l_art, l_art_nested (KIT)
! V1_23        2013/10/02 Ulrich Schaettler
!  If a new set of vcoord is used, only abort in case of GRIB2 output
!  Set values of press_level_d which are not used to 0
!  Rename vcoord_out, refatm_out to vcoord, refatm (as is in COSMO-Model)
!  Set default for l_ke_in_gds to .TRUE.
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

!------------------------------------------------------------------------------

USE data_grid_lm,     ONLY :   &
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
    startlat,    & ! transformed latitude of the lower left grid point
                   ! of the local domain (in degrees, N>0)
    startlon,    & ! transformed longitude of the lower left grid point
                   ! of the local domain (in degrees, E>0)
    endlat,      & ! transformed latitude of the upper right grid point
                   ! of the local domain (in degrees, N>0)
    endlon,      & ! transformed longitude of the upper right grid point
                   ! of the local domain (in degrees, E>0)
    ielm_tot,    & ! ie for LM, whole area
    jelm_tot,    & ! je for LM, whole area
    kelm_tot,    & ! number of vertical levels for LM
    ielm,        & ! ie for LM, local domain
    jelm,        & ! je for LM, local domain
    kelm,        & ! ke for LM
    ie2lm_max,   & ! Max. of ie2lm on all processors
    je2lm_max,   & ! Max. of je2lm on all processors
    ie2lm_tot,   & ! = ielm_tot + 2
    je2lm_tot,   & ! = jelm_tot + 2
    ij2lm,       & ! = je2lm_tot*je2lm_tot, will be removed soon
    kedim,       & !
    ie2lm,       & !
    je2lm

USE data_grid_lm,     ONLY :   &
    istartpar,   & ! start index for computations in the parallel program
    iendpar,     & ! end index for computations in the parallel program
    jstartpar,   & ! start index for computations in the parallel program
    jendpar,     & ! end index for computations in the parallel program
    ie_ext,      & ! west-east size of fields with external parameters
    je_ext,      & ! north-south size of fields with external parameters
                   ! these values may be larger than the actual sizes
    cw_so_rel_lm,& ! artifical volumetric soil water content (0-1) profile
    ke_soil_lm,  & ! number of levels in multi-layer soil model in output
    czmls_lm,    & ! depth of the soil layers in meters in output
    czhls_lm       ! depth of the half soil layers in meters in output

!------------------------------------------------------------------------------

USE data_grid_in,     ONLY :   &
    ids,         & ! start index of diamonds (ids = 1)
    ide,         & ! end index of diamonds (ide = 10)
    ni_gme,      & ! resolution of GME
    i3e_gme,     & ! number of levels in the vertical
    ilim1,       & ! decomposition limits in x-direction (formerly 1)
    ilim2,       & ! decomposition limits in y-direction (formerly 2)
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
    ig1s,        & ! start index of array-dimension in x-direction
    ig1sm1,      & ! = ig1s - 1
    ig1sm2,      & ! = ig1s - 2
    ig1e,        & ! end index of array-dimension in x-direction
    ig1ep1,      & ! = ig1e + 1
    ig1ep2,      & ! = ig1e + 2
    ig2s,        & ! start index of array-dimension in y-direction
    ig2sm1,      & ! = ig2s - 1
    ig2sm2,      & ! = ig2s - 2
    ig2e,        & ! end index of array-dimension in y-direction
    ig2ep1,      & ! = ig2e + 1
    ig2ep2,      & ! = ig2e + 2
    nbpe,        & ! Number of neighbor poleward east
    nbaw,        & ! Number of neighbor antipoleward west
    nbpw,        & ! Number of neighbor poleward west
    nbae,        & ! Number of neighbor antipoleward east
    max_gme_core   ! size of global GME core

USE data_grid_in,     ONLY :   &
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
    startlat_in,    & ! transformed latitude of the lower left grid point
                      ! of the local domain (in degrees, N>0)
    startlon_in,    & ! transformed longitude of the lower left grid point
                      ! of the local domain (in degrees, E>0)
    endlat_in,      & ! transformed latitude of the upper right grid point
                      ! of the local domain (in degrees, N>0)
    endlon_in,      & ! transformed longitude of the upper right grid point
                      ! of the local domain (in degrees, E>0)
    latitudes_in,   & ! latitudes of the input data
    longitudes_in,  & ! longitudes of the input data
    slatitudes_in,  & ! latitudes of the input data
    slongitudes_in, & ! longitudes of the input data
    lushift_in,     & ! shift of u-velocity due to grid staggering
    lvshift_in,     & ! shift of v-velocity due to grid staggering
    east_add_in,    & ! add an extra column to the East
    west_add_in,    & ! add an extra column to the West
    south_add_in,   & ! add an extra column to the South
    north_add_in      ! add an extra column to the North

USE data_grid_in,     ONLY :   &
    pollat_in,      & ! latitude of the rotated north pole (in degrees, N>0)
    ie_in_tot,      & ! ie for input grid, total domain
    je_in_tot,      & ! je for input grid, total domain
    ke_in_tot,      & ! ke for input grid
    ie_in_max,      & ! Max. of ie_in (local) on all processors
    je_in_max,      & ! Max. of je_in (local) on all processors
    ie_in,          & ! ie for input grid, local domain
    je_in,          & ! je for input grid, local domain
    ke_in,          & ! ke for input grid
    ke_hybrid,      & ! number of hybrid levels when interpolating from pressure levels
    ke_pressure,    & ! number of pressure levels when interpolating from pressure levels
    kedim_in,       & ! MAX (ke_in, ke_hybrid)
    press_level,    & ! list of available pressure levels (in Pa) for GFS
    nlevskip,       & ! number of levels to skip at top of input model
    ke_soil_in,     & ! number of input levels in multi-layer soil model
    czmls_in,       & ! depth of the coarse grid soil layers in meters
    czhls_in,       & ! depth of the coarse grid half soil layers in meters
    ie_ext_in,      & ! west-east size of fields with external parameters
    je_ext_in         ! north-south size of fields with external parameters

!------------------------------------------------------------------------------

USE data_int2lm_io,          ONLY :   &
    nlocaldefnr,  & ! local definition number for GRIB2 local section template
    nprocess_ini, & ! type of database for LM initial data
    nprocess_bd,  & ! type of database for LM boundary data
    ncenter,      & ! originating center identification
    nsubcenter,   & ! originating subcenter identification
    nincwait,     & ! if ready-file is not available wait nincwait seconds
                    ! until next attempt
    nmaxwait,     & ! if ready-file is not available after nmaxwait seconds,
                    ! abort the program
    l_ke_in_gds,  & ! explicit GDS entry for number of model levels
    lchkin,       & ! logical for print of check-values (max,min,mean)
                    ! of GME-fields
    lchkout,      & ! logical for print of check-values (max,min,mean)
                    ! of LM/HM-fields
    lcheck_bmm,   & ! if .TRUE., check bitmaps for mass grid points
    lcheck_bmu,   & ! if .TRUE., check bitmaps for u-grid points
    lcheck_bmv,   & ! if .TRUE., check bitmaps for v-grid points
    ytrans_in,    & ! directory for reading ready-files
    ytrans_out,   & ! directory for writing ready-files
    nrbit,        & ! packrate for grib records
    nbitmap,      & ! user dimension for bitmap
    ydate_ini,    & ! start of the forecast
    lmmss_ini,    & ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                    ! if .FALSE. 10 digits date format (YYYYMMDDHH)
                    ! for ydate_ini and result files of INT2LM
    ydate_bd,     & ! start of the forecast from the boundary model
    lmmss_bd        ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                    ! if .FALSE. 10 digits date format (YYYYMMDDHH)
                    ! for ydate_bd  and input  files of INT2LM

USE data_int2lm_io,          ONLY :   &
    ytunit_in,    & ! time unit for input data
    ytunit_out,   & ! time unit for output data
    ylmext_cat,   & ! catalog-name of the file with external LM parameters
    yinext_cat,   & ! catalog-name of the file with external GME parameters
    ylm_cat,      & ! catalog-name of the LM files
    yin_cat,      & ! catalog-name of the GME files
    ylm_hhl,      & ! name of the file with output (GRIB2) HHL fields
    yin_hhl,      & ! name of the file with input (GRIB2) HHL fields
    ylmext_lfn,   & ! name of the file with external LM parameters
    yinext_lfn,   & ! name of the file with external GME parameters
    ymode_read,   & ! mode for opening the (read) Grib files
    ymode_write,  & ! mode for opening the (write) Grib files
    nmaxlist,     & ! maximal number of elements in output lists
    numlist_ini,  & ! number of elements in youtlist_ini
    numlist_bd,   & ! number of elements in youtlist_bd
    numlist_hhl,  & ! number of elements in youtlist_hhl
    youtlist_ini, & ! list of output variables for initial data
    youtlist_bd,  & ! list of output variables for boundary data
    youtlist_hhl, & ! list of output variables for HHL file
    yinput_type,      & ! type of input data: 'forecast', 'analysis' or 'ana_init'
    ylmext_form_read, & ! input format of external LM data
    yinext_form_read, & ! input format of external boundary data
    yin_form_read,    & ! input format of boundary data
    ylm_form_write,   & ! output format of LM data
    yncglob_institution,   & ! originating center name
    yncglob_title,         & ! title string for the output
    yncglob_source,        & ! program name and version
    yncglob_project_id,    & ! identification of the project of simulation
    yncglob_experiment_id, & ! identification of the experiment of simulation
    yncglob_contact,       & ! contact e.g. email address
    yncglob_references,    & ! URL, report etc.
    ncglob_realization       ! number of the realization of the experiment

!------------------------------------------------------------------------------

USE data_int2lm_control,     ONLY :   &
    nstart,         & ! start time (in time steps)
    nstop,          & ! end time (in time steps)
    nhori,          & ! number of sectors for the horizont array by the topographic
                      ! correction of the radiation
    nincbound,      & ! time step increment between two datasets
    ninput,         & ! unit number for NAMELIST input file
    noutput,        & ! unit number for output file
    linitial,       & ! if .TRUE., initial data for LM
    lboundaries,    & ! if .TRUE., lateral boundaries for LM
    leps_bc,        & ! switch ensemble mode on/off for boundary data
    itype_calendar, & ! for specifying the calendar used
    lgme2lm,        & ! if .TRUE., gme->lm
    lgfs2lm,        & ! if .TRUE., gfs->lm
    lgsm2lm,        & ! if .TRUE., gsm->lm
    lec2lm,         & ! if .TRUE., ec ->lm
    llm2lm,         & ! if .TRUE., lm ->lm
    lum2lm,         & ! if .TRUE., um ->lm
    lhir2lm,        & ! if .TRUE., hirlam ->lm
    lhm2lm,         & ! if .TRUE., hm ->lm
    lcm2lm,         & ! if .TRUE., cm ->lm
    yinput_model,   & ! string to identify the input model
    luse_t_skin,    & ! if .TRUE., use ECMWF skin temperature for surface
    lante_0006,     & ! if .TRUE., force to use ECMWF dataset before 27 June 2000
    lpost_0006,     & ! if .TRUE., force to use ECMWF dataset after 27 June 2000
    luvcor,         & ! if .TRUE., correct winds for given surface pressure  
                      ! tendency
    lprog_qi,       & ! if .TRUE., interpolate qi from GME to LM grid
    lmixcld,        & ! if .TRUE., qi added in grh instead of being directly interp.
    l_smi,          & ! if .TRUE., interpolate soil moisture with SMI
    lprog_qr_qs,    & ! if .TRUE., interpolate qr,qs from LM to LM grid
    lprog_qg,       & ! if .TRUE., interpolate qg from LM to LM grid
    lprog_rho_snow, & ! if .TRUE., interpolate rho_snow from GME to LM grid
    l_topo_z,       & ! if .TRUE., additional smoothing of the topography
    qvmin,          & ! minimum value of water vapor (security)
    qcmin,          & ! minimum value of cloud water (security)
    qimin,          & ! minimum value of cloud ice content (security)
    lvertwind_ini,  & ! if .TRUE., compute vertical wind for LM for initial data
    lvertwind_bd,   & ! if .TRUE., compute vertical wind for LM for boundary data
    lbdclim           ! if .TRUE., special boundary data for climate mode

USE data_int2lm_control,     ONLY :   &
    lforest,      & ! if .TRUE., run with forest (evergreen and deciduous)
    lurban,       & ! if .TRUE., use data on urban fraction
    lsso,         & ! process parameters for sso scheme
    lradtopo,     & ! process parameters for topographic correction of radiation
    lseaice,      & ! if .TRUE., run with sea ice model
    llake,        & ! if .TRUE., run with lake
    llake_coldstart,& ! if .TRUE., initialize prognostic lake variables for cold start
    lemiss,       & ! if .TRUE., run with external parameter for surface emissivity
    itype_albedo, & ! choose treatment of solar surface albedo
    lstomata,     & ! if .TRUE., run with external parameter for stomata resistance
    llbc_smooth,  & ! if .TRUE., run with smooth orography transition at LB
    nlbc_smooth,  & ! number of grip points for smooth orography transition at LB
    lroutine,     & ! if .TRUE., routine-job
    lclock,       & ! if .TRUE., system clock is present
    ltime,        & ! detailled timings of the program are given
    ltime_mean,   & ! if .TRUE., mean values of the timings are printed
    ltime_proc,   & ! if .TRUE., timings for each processor are printed
    idbg_level,   & ! to control verbosity of output
    lprintdeb_all,& ! whether all PEs print debug output
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
    lfilter_oro,  & ! if .TRUE., filter the orography
    lxso_first,   & ! if .TRUE., do eXtra smoothing of orography first
    lfilter_pp,   & ! if .TRUE., filter the pressure deviation after vertical
                    !            interpolation
    lbalance_pp     ! if .TRUE., compute a hydrostatic balanced pp after
                    !            vertical interpolation in LM2LM

USE data_int2lm_control,     ONLY :   &
    norder_filter,& ! order of the orography filtering
    eps_filter,   & ! parameter for orography filtering
    ilow_pass_oro,& ! type of low-pass filter for orography
    numfilt_oro,  & ! number of sequential applications of filter
    ilow_pass_xso,& ! type of low-pass filter for extra smoothing of steep oro.
    numfilt_xso,  & ! number of sequential applications of xso-filter
    rxso_mask,    & ! mask for extra smoothing of steep oro.: dh > rxso_mask
    rfill_valley, & ! mask for valley filling: dh > rfill_valley
    ifill_valley, & ! type of valley filling: 1=MIN (before), 2=MIN (after)
    lt_cl_corr,   & ! if .TRUE., a height-correction of T_CL is performed
    nlev1pr,      & ! k-index for the print of first  model layer
    nlev2pr,      & ! k-index for the print of second model layer
    lbd_frame,    & ! if .TRUE., boundary fields include only frames
    npstrframe,   & ! thickness of output frames
    nvers,        & ! for documenting purposes in Grib Code
    nl_soil_in,   & ! number of soil layers in GME
    nl_soil_lm,   & ! number of soil layers in LM, resp. HM

    kcontrol_fi,  & ! control level for geopotential
    dt,           & ! time step used in the LM
    pcontrol_fi,  & ! pressure of control level for geopotential
    timings,      & ! for storing the times for different parts of the program
    msoilgrib_lm, & ! grib coded depth of main soil levels in centimeters
    msoilgrib_in, & ! grib coded depth of main soil levels in centimeters
                    ! (careful: the first level will be coded with 1,
                    !           but is in the depth of 0.5 cm!)
    itype_w_so_rel,&! type of relative soil moisture input (0,1,2,3,4)
    itype_aerosol,& ! to choose treatment of aerosol parameters
    lmulti_layer_lm,&! if .TRUE., compute soil fields for multi-layer soil
                    ! model in the outgoing data
    lmulti_layer_in,&! if .TRUE., incoming data have soil fields from the
                    ! multi-layer soil model
    l_cressman,   & ! indicator for using a cressmann scheme during M interpolation
    l_bicub_spl,  & ! indicator for using a bicubic spline interpolation
                    ! for variables with interpolation code ipc(1:1)='Q'
    l_art      ,  & ! if .TRUE., interpolate additional cosmo-art fields
    l_art_nested, & ! if true nesting lm-art2lm-art
    itype_t_cl,   & ! to choose origin and treatment of deep soil temperature
    itype_rootdp, & ! to choose treatment of root depth
    itype_ndvi,   & ! to choose treatment of surface parameters (plcov, lai)
! iso code
    liso            ! if .TRUE., include variables for water isotope simulation
! end iso code

!------------------------------------------------------------------------------

USE data_epsctl,     ONLY :   &
    iepsmem_bc,    & ! ID of the member in the ensemble of boundary
                     ! conditions (iepsmem_bc >= 0)
    iepstyp_bc,    & ! ID of the boundary ensemble generation type
                     ! (iepstyp_bc >= 0)
    iepstot_bc,    & ! total number of boundary ensemble members (iepstot_bc >= 0)
    lchk_bc_typ      ! check member ID of input data (if leps_bc=.TRUE.)

!------------------------------------------------------------------------------

USE data_profiles,    ONLY :   &
    nmaxgp,       & ! maximal number of grid points for grid point output
    lprgp,        & ! logical for print at selected grid points
    ngp_tot,      & ! number of selected grid points in the total domain
    ngp,          & ! number of selected grid points in this local domain
    igp,          & ! i-indices    local domain
    jgp,          & ! j-indices    local domain
    igp_tot,      & ! i-indices    total domain
    jgp_tot         ! j-indices    total domain

!------------------------------------------------------------------------------

USE data_int2lm_parallel,    ONLY :   &
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    lasync_io,       & ! if .TRUE.: the model runs with extra PEs for
                       ! asynchronous IO
    nprocx,          & ! number of processors in x-direction
    nprocy,          & ! number of processors in y-direction
    nprocio,         & ! number of extra processors for asynchronous IO
    nproc,           & ! total number of processors: nprocx * nprocy + nprocio
    num_compute,     & ! number of compute PEs
    num_io,          & ! number of IO PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ncomm_type,      & ! type of communication for boundary exchange
    my_world_id,     & ! rank of this subdomain in the global communicator
    icomm_world,     & ! communicator that belongs to igroup_world
    lreorder,        & ! during the creation of the virtual topology the
                       ! ranking of the processors may be reordered
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_integers,    & ! determines the correct INTEGER type used in the
                       ! model for MPI
    imp_character,   & ! determines the correct CHARACTER type used in the
                       ! model for MPI
    imp_logical        ! determines the correct LOGICAL   type used in the
                       ! model for MPI

!------------------------------------------------------------------------------

USE environment,        ONLY:  get_free_unit, release_unit
USE mpe_io,             ONLY:  mpe_readnldb
USE parallel_utilities, ONLY:  distribute_values

USE vgrid_refatm_utils, ONLY:                                               &
    imax_vcoordtype, vcoord_defaults, vcoord,                               &
    imax_refatmtype, refatm_defaults, refatm, rundefined,                   &
    set_vcoord_defaults, set_refatm_defaults,                               &
    lanalyt_calc_t0p0, nfltvc, svc1,svc2, khmax, vcoord_d,                  &
    lnewVGrid, lhhl_lm_read, lhhl_in_read

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Public and Private Subroutines

PUBLIC   read_namelists

PRIVATE  input_contrl, input_grid_in, input_lmgrid, input_data,             &
         input_prictr, input_epsctl,  input_dbase

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in src_namelists for the organization
!------------------------------------------------------------------------------

SUBROUTINE read_namelists (iglob_error, yglob_error)

!------------------------------------------------------------------------------
!
! Description:
!  read_namelists is the driver routine for the input of all namelist groups
!  for the interpolation program. The different groups are:
!   - /CONTRL/   basic control: date, start, end, specification of tasks
!   - /GRID_IN/  specification of input grid: gme or other regular grids
!   - /LMGRID/   specification of (fine) LM grid
!   - /DBASE/    specification of data base job
!   - /DATA/     directories, filenames, ...
!   - /PRICTR/   for debug- and control-output
!   - /EPSCTL/   for ensemble specifications
!
! Method:
!  MPI-Task 0 sets defaults for all namelist variables and then reads the
!  input. After an error and consistency check the variables are distributed
!  to all other MPI tasks.
!
! Input files:
!  File INPUT containing all Namelist-groups.
!
! Output files:
!  File OUTPUT containing the variables of every Namelist-group with default
!  and actual value and further information gathered during the program
!  run.
!
!==============================================================================
!
! Parameterlist
CHARACTER (LEN= 80),      INTENT(OUT) ::  &
  yglob_error        ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  iglob_error        ! error status

! Local scalars:
INTEGER (KIND=iintegers)   ::       &
  izerrstat,       & ! error status of Namelist input
  izerrstat1,      & ! error status of Namelist input
  izbuflen,        & ! length of buffer distributing the NAMELIST
  nziostat,        & ! for error-code of I/O
  nzstat             ! for error-code on allocation

CHARACTER (LEN= 8)         ::       &
  youtput,   & ! output file for protocolling the program run
  yinput       ! Namelist INPUT file

!
! Local arrays:
! Buffers for distributing the Namelists
INTEGER (KIND=iintegers), ALLOCATABLE ::   intbuf  (:)
REAL    (KIND=ireals)   , ALLOCATABLE ::   realbuf (:)
LOGICAL                 , ALLOCATABLE ::   logbuf  (:)
CHARACTER (LEN=100)     , ALLOCATABLE ::   charbuf (:)

!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!- Begin SUBROUTINE read_namelists
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  iglob_error = 0_iintegers
  yglob_error = '     '

  izerrstat   = 0_iintegers
  izerrstat1  = 0_iintegers

  ! Allocate space for sending buffers
  izbuflen  = 1000
  ALLOCATE ( intbuf(izbuflen)   , STAT=nzstat )
  intbuf (:) = 0
  ALLOCATE ( realbuf(izbuflen)  , STAT=nzstat )
  realbuf(:) = 0.0_ireals
  ALLOCATE ( logbuf(izbuflen)   , STAT=nzstat )
  logbuf (:) = .FALSE.
  ALLOCATE ( charbuf(izbuflen)  , STAT=nzstat )
  charbuf(:) = ' '

!------------------------------------------------------------------------------
! Section 2: Prepare NAMELIST-input
!------------------------------------------------------------------------------

  IF (my_world_id == 0) THEN
    PRINT *,'    INPUT OF THE NAMELISTS'

    ! get free unit numbers for the files INPUT and OUTPUT
    CALL get_free_unit (ninput)
    CALL get_free_unit (noutput)

    ! Open files for input of the NAMELISTS and control output
    yinput   = 'INPUT   '
    youtput  = 'OUTPUT  '

    OPEN(ninput, FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
         IOSTAT=nziostat)
    IF(nziostat /= 0) THEN
      yglob_error = ' ERROR    *** Error while opening file INPUT *** '
      iglob_error = 1001
      RETURN
    ENDIF

    OPEN(noutput, FILE=youtput, FORM=  'FORMATTED', STATUS='UNKNOWN',  &
         IOSTAT=nziostat)
    IF(nziostat /= 0) THEN
      yglob_error = ' ERROR    *** Error while opening file OUTPUT *** '
      iglob_error = 1001
      RETURN
    ENDIF
    REWIND noutput

    ! Print a headline in file OUTPUT
    WRITE (noutput, '(A2)')  '  '
    WRITE (noutput, '(A55)')                                                 &
                   '  The NAMELIST for INT2LM contains the following values:'
    WRITE (noutput, '(A55)')                                                 &
                   '  ======================================================'
    WRITE (noutput, '(A2)')  '  '
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Read the NAMELIST-variables
!------------------------------------------------------------------------------

  ! Read all NAMELIST-groups
  CALL input_contrl (realbuf, intbuf, logbuf, charbuf, izbuflen,          &
                     noutput, ninput, izerrstat1)
  IF (izerrstat1 < 0) THEN
    izerrstat  = 1000
    izerrstat1 = 0
    IF (my_world_id == 0) THEN
      PRINT *, ' ERROR    *** while reading NAMELIST Group /CONTRL/ ***'
    ENDIF
  ELSE
    izerrstat = izerrstat + izerrstat1
  ENDIF

  CALL input_grid_in(realbuf, intbuf, logbuf, charbuf, izbuflen,          &
                     noutput, ninput, izerrstat1)
  IF (izerrstat1 < 0) THEN
    izerrstat  = 1000
    izerrstat1 = 0
    IF (my_world_id == 0) THEN
      PRINT *, ' ERROR    *** while reading NAMELIST Group /GRID_IN/ ***'
    ENDIF
  ELSE
    izerrstat = izerrstat + izerrstat1
  ENDIF

  CALL input_lmgrid (realbuf, intbuf, logbuf, charbuf, izbuflen,          &
                     noutput, ninput, izerrstat1)
  IF (izerrstat1 < 0) THEN
    izerrstat  = 1000
    izerrstat1 = 0
    IF (my_world_id == 0) THEN
      PRINT *, ' ERROR    *** while reading NAMELIST Group /LMGRID/ ***'
    ENDIF
  ELSE
    izerrstat = izerrstat + izerrstat1
  ENDIF

  CALL input_dbase  (realbuf, intbuf, logbuf, charbuf, izbuflen,          &
                     noutput, ninput, izerrstat1)
  IF (izerrstat1 < 0) THEN
    izerrstat  = 1000
    izerrstat1 = 0
    IF (my_world_id == 0) THEN
      PRINT *, ' ERROR    *** while reading NAMELIST Group /DBASE/ ***'
    ENDIF
  ELSE
    izerrstat = izerrstat + izerrstat1
  ENDIF

  CALL input_data   (realbuf, intbuf, logbuf, charbuf, izbuflen,          &
                     noutput, ninput, izerrstat1)
  IF (izerrstat1 < 0) THEN
    izerrstat  = 1000
    izerrstat1 = 0
    IF (my_world_id == 0) THEN
      PRINT *, ' ERROR    *** while reading NAMELIST Group /DATA/ ***'
    ENDIF
  ELSE
    izerrstat = izerrstat + izerrstat1
  ENDIF

  CALL input_prictr (realbuf, intbuf, logbuf, charbuf, izbuflen,          &
                     noutput, ninput, izerrstat1)
  IF (izerrstat1 < 0) THEN
    izerrstat  = 1000
    izerrstat1 = 0
    IF (my_world_id == 0) THEN
      PRINT *, ' ERROR    *** while reading NAMELIST Group /PRICTR/ ***'
    ENDIF
  ELSE
    izerrstat = izerrstat + izerrstat1
  ENDIF

  IF (leps_bc) THEN
    CALL input_epsctl (realbuf, intbuf, logbuf, charbuf, izbuflen, &
                       noutput, ninput, izerrstat1)
    IF (izerrstat1 < 0) THEN
      izerrstat  = 1000
      izerrstat1 = 0
      IF (my_world_id == 0) THEN
        PRINT *, ' ERROR    *** while reading NAMELIST Group /EPSCTL/ ***'
      ENDIF
    ELSE
      izerrstat = izerrstat + izerrstat1
    ENDIF
  ENDIF

  IF (my_world_id == 0) THEN
    ! Close file for input of the NAMELISTS
    CLOSE (ninput  , STATUS='KEEP')
  ENDIF

  IF (izerrstat /= 0) THEN
    iglob_error = 1002
    yglob_error = ' ERROR    *** Wrong values occured in NAMELIST input ***'
    RETURN
  ENDIF

  ! Deallocate space for sending buffers
  DEALLOCATE ( intbuf , STAT=nzstat )
  DEALLOCATE ( realbuf, STAT=nzstat )
  DEALLOCATE ( logbuf , STAT=nzstat )
  DEALLOCATE ( charbuf, STAT=nzstat )

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE read_namelists

!==============================================================================
!+ Module procedure in "setup" for the input of NAMELIST contrl
!------------------------------------------------------------------------------

SUBROUTINE input_contrl (realbuf, intbuf, logbuf, charbuf, ibuflen,          &
                         nuout, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group contrl. 
!   The group contrl contains variables defining the date, start- and end
!   of the run and the specifications of what has to be done.
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for 
!   consistency. If wrong input values are detected the program prints 
!   an error message. The program is not stopped in this routine but an 
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists. 
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.    
!   Both, default and input values are written to the file OUTPUT.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    ibuflen,      & ! dimension of the buffers
    nuout,        & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Sending buffers (are given to routine distribute_values)
  INTEGER (KIND=iintegers), INTENT(INOUT)    ::  intbuf  (ibuflen)
  REAL    (KIND=ireals)   , INTENT(INOUT)    ::  realbuf (ibuflen)
  LOGICAL                 , INTENT(INOUT)    ::  logbuf  (ibuflen)
  CHARACTER (LEN=100)     , INTENT(INOUT)    ::  charbuf (ibuflen)

! Local variables
  INTEGER (KIND=iintegers)   ::  i, invar

! Variables for default values
  INTEGER (KIND=iintegers)   ::       &
    nprocx_d,        & ! number of processors in east-west direction
    nprocy_d,        & ! number of processors in north-south direction
    nprocio_d,       & ! number of extra processors for asynchronous IO
    nboundlines_d,   & ! number of boundary lines at each side of a subdomain
    ncomm_type_d,    & ! type of communication
    nincwait_d,      & ! seconds to wait until next attempt if a ready file is
                       ! not available
    nmaxwait_d,      & ! maximum seconds to wait until abort if a ready file is
                       ! not available
    norder_filter_d, & ! order of the orography filtering
    ilow_pass_oro_d, & ! type of low-pass filter for orography
    numfilt_oro_d,   & ! number of sequential applications of filter
    ilow_pass_xso_d ,& ! type of low-pass filter for extra smoothing steep oro.
    numfilt_xso_d,   & ! number of sequential applications of xso-filter
    ifill_valley_d,  & ! type of valley filling: 1=MIN (before), 2=MIN (after)
    nlbc_smooth_d,   & ! number of grip points for smooth orography transition at LB
    nhori_d,         & ! number of sectors for the horizont array by the topographic
                       ! correction of the radiation
    itype_calendar_d,& ! for specifying the calendar used
    itype_w_so_rel_d,& ! type of relative soil moisture input (0,1,2)
    itype_t_cl_d,    & ! to choose origin and treatment of deep soil temperature
    itype_rootdp_d,  & ! to choose treatment of root depth
    itype_ndvi_d,    & ! to choose treatment of surface parameters (plcov, lai)
    itype_aerosol_d, & ! to choose treatment of aerosol parameters
    itype_albedo_d,  & ! choose treatment of solar surface albedo
    idbg_level_d       ! to control verbosity of debug output

  REAL (KIND=ireals)         ::       &
    hstart_d,        & ! first hour to be processed  (default)
    hstart,          & ! first hour to be processed 
    hstop_d,         & ! last hour to be processed (default)
    hstop,           & ! last hour to be processed
    htest,           & ! just for testing
    hincbound_d,     & ! increment in hours for the processing (default)
    hincbound,       & ! increment in hours for the processing
    eps_filter_d,    & ! parameter for orography filtering
    rxso_mask_d,     & ! mask for extra smoothing of steep oro.: dh > rxso_mask
    rfill_valley_d,  & ! mask for valley filling: dh > rfill_valley
    qvmin_d,         & ! minimum value of water vapor (security)
    qcmin_d,         & ! minimum value of cloud water (security)
    qimin_d            ! minimum value of cloud ice content (security)

  CHARACTER (LEN=14)         ::       &
    ydate_ini_d,  & ! start of the forecast
    ydate_bd_d      ! start of the forecast from the boundary model from which
                    ! data are used

  CHARACTER (LEN=250)        ::       &
    ytrans_in_d,  & ! directory for reading ready-files
    ytrans_out_d    ! directory for writing ready-files

  LOGICAL                    ::       &
    ldatatypes_d, & ! if .TRUE.: use MPI-Datatypes for some communications
    luse_t_skin_d,& ! if .TRUE., use ECMWF skin temperature for surface
    lante_0006_d, & ! if .TRUE., force to use ECMWF dataset before 27 June 2000
    lpost_0006_d, & ! if .TRUE., force to use ECMWF dataset after 27 June 2000
    luvcor_d,     & ! if .TRUE., correct winds for given surface pres. tendency
    lvertwind_ini_d,&!if .TRUE., compute vertical wind for LM for initial data
    lvertwind_bd_d,&! if .TRUE., compute vertical wind for LM for boundary data
    lprog_qi_d,   & ! if .TRUE., interpolate qi to LM grid
    lprog_qr_qs_d,& ! if .TRUE., interpolate qr,qs to LM grid
    lprog_qg_d,   & ! if .TRUE., interpolate qg to LM grid
    lprog_rho_snow_d, & ! if .TRUE., interpolate rho_snow to LM grid
    linitial_d,   & ! if .TRUE., initial data for LM
    lboundaries_d,& ! if .TRUE., lateral boundaries for LM
    leps_bc_d,    & ! if .TRUE., ensemble mode for boundary data
    ltime_mean_d, & ! if .TRUE., mean values of the timings are printed 
    ltime_proc_d, & ! if .TRUE., timings for each processor are printed
    lmulti_layer_lm_d,&! compute data for multi-layer soil model
    lmulti_layer_in_d,&! data from multi-layer soil model in the incoming data
    lfilter_oro_d,& ! if .TRUE., filter the orography
    lxso_first_d, & ! if .TRUE., do eXtra smoothing of orography first
    lfilter_pp_d, & ! if .TRUE., filter the pressure deviation after vertical
    lbalance_pp_d,& ! if .TRUE., compute a hydrostatic balanced pp after
                    !            vertical interpolation in LM2LM
    lt_cl_corr_d, & ! if .TRUE., a height-correction of T_CL is performed
    lbdclim_d,    & ! if .TRUE., special boundary data for climate mode
    lforest_d,    & ! if .TRUE., run with forest (evergreen and deciduous)
    lurban_d,     & ! if .TRUE., use data on urban fraction
    lsso_d,       & ! process parameters for sso scheme
    lradtopo_d,   & ! process parameters for topographic correction of radiation
    lseaice_d,    & ! if .TRUE., run with sea ice model
    llake_d,      & ! if .TRUE., run with lake
    llake_coldstart_d,& ! if .TRUE., initialize prognostic lake variables for cold start
    lemiss_d,     & ! if .TRUE., run with external parameter for surface emissivity
    lstomata_d,   & ! if .TRUE., run with external parameter for stomata resistance
    llbc_smooth_d,& ! if .TRUE., run with smooth orography transition to LB
    lroutine_d      ! if .TRUE., routine-job

  LOGICAL                    ::       &
    l_cressman_d, & ! switch for using a cressmann scheme during M interpolation
    l_bicub_spl_d,& ! switch for using a bicubic spline interpolation
    l_art_d,      & ! switch for additional cosmo-art fields
    l_art_nested_d,&! switch for cosmo-art2cosmo-art
    l_smi_d,      & ! if .TRUE., interpolate soil moisture with SMI
    lmixcld_d,    & ! if .TRUE., qi added in grh instead of being directly interp.
    lasync_io_d,  & ! if .TRUE.: the model runs with extra PEs for asynchr. IO
    lreorder_d,   & ! during the creation of the virtual topology the
                    ! ranking of the processors may be reordered
    l_topo_z_d,   & ! if .TRUE., additional smoothing of the topography
                    ! for the LM-Z coordinate Version
! iso code
    liso_d,       & ! if .TRUE., include variables for water isotope simulation
! end iso code
    lprintdeb_all_d ! whether all PEs print debug output

  CHARACTER (LEN=5)                ::           &
    yinput_model_d  ! string to identify the input model

  INTEGER (KIND=iintegers)   :: ierr, nzylen, ierror, k
  LOGICAL                    :: lzequiv
  CHARACTER (LEN= 8)         :: ydatearg, ytime
  CHARACTER (LEN=10)         :: ydate, ytimearg
  CHARACTER (LEN=15)         :: yroutine

! Define the namelist group
  NAMELIST /contrl/ &
    ydate_ini,     ydate_bd,      hstart,        hstop,         hincbound,     &
    nincwait,      nmaxwait,      ytrans_in,     ytrans_out,    nprocx,        &
    nprocy,        nprocio,       nboundlines,   luse_t_skin,   lante_0006,    &
    lpost_0006,    luvcor,        lvertwind_ini, lvertwind_bd,  lprog_qi,      &
    lprog_qr_qs,   lprog_qg,      lprog_rho_snow,linitial,      lboundaries,   &
    leps_bc,       itype_calendar,lfilter_oro,   lxso_first,    lfilter_pp,    &
    lbalance_pp,   norder_filter, eps_filter,    ilow_pass_oro, numfilt_oro,   &
    ilow_pass_xso, numfilt_xso,   rxso_mask,     rfill_valley,  ifill_valley,  &
    qvmin,         qcmin,         qimin,         ltime_mean,    ltime_proc,    &
    yinput_model,  lbdclim,       lforest,       lsso,          lradtopo,      &
    llake,         llbc_smooth,   nlbc_smooth,   lroutine,      lasync_io,     &
    lreorder,      itype_w_so_rel,lmulti_layer_lm, lmulti_layer_in,            &
    lt_cl_corr,    ldatatypes,    ncomm_type,    l_topo_z,      idbg_level,    &
    lprintdeb_all, l_cressman,    l_bicub_spl,   l_art,         itype_t_cl,    &
    itype_rootdp,  itype_ndvi,    itype_aerosol, lemiss,        lstomata,      &
    l_smi,         lmixcld,       nhori,         lseaice,       llake_coldstart,&
    lurban,        itype_albedo,  l_art_nested,                                &
! iso code
    liso
! end iso code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_contrl
!------------------------------------------------------------------------------

yroutine = 'input_contrl'
ierrstat = 0_iintegers
ierror   = 0_iintegers

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  ydate_ini_d         =  '              '
  ydate_bd_d          =  '              '
  hstart_d            =   0.0_ireals
  hstop_d             =   0.0_ireals
  hincbound_d         =   0.0_ireals
  eps_filter_d        =  10.0_ireals
  rxso_mask_d         =   0.0_ireals
  rfill_valley_d      =   0.0_ireals
  nincwait_d          =     0_iintegers
  nmaxwait_d          =     0_iintegers
  ytrans_in_d         =  '   '
  ytrans_out_d        =  '   '
  nprocx_d            =     1_iintegers
  nprocy_d            =     1_iintegers
  nprocio_d           =     0_iintegers
  ncomm_type_d        =     1_iintegers
  nboundlines_d       =     1_iintegers
  norder_filter_d     =     5_iintegers
  ilow_pass_oro_d     =     1_iintegers
  numfilt_oro_d       =     1_iintegers
  ilow_pass_xso_d     =     0_iintegers
  numfilt_xso_d       =     1_iintegers
  ifill_valley_d      =     0_iintegers
  luse_t_skin_d       =  .FALSE.
  lante_0006_d        =  .FALSE.
  lpost_0006_d        =  .FALSE.
  luvcor_d            =  .TRUE.
  lvertwind_ini_d     =  .TRUE.
  lvertwind_bd_d      =  .FALSE.
  lprog_qi_d          =  .FALSE.
  l_smi_d             =  .FALSE.
  lmixcld_d           =  .FALSE.
  lprog_qr_qs_d       =  .FALSE.
  lprog_qg_d          =  .FALSE.
  lprog_rho_snow_d    =  .FALSE.
  qvmin_d             =  1.0E-12_ireals
  qcmin_d             =  1.0E-12_ireals
  qimin_d             =  1.0E-12_ireals
  linitial_d          =  .FALSE.
  lboundaries_d       =  .TRUE.
  leps_bc_d           =  .FALSE.
  lfilter_oro_d       =  .FALSE.
  lxso_first_d        =  .FALSE.
  lfilter_pp_d        =  .FALSE.
  lbalance_pp_d       =  .FALSE.
  lt_cl_corr_d        =  .FALSE.
  ltime_mean_d        =  .FALSE.
  ltime_proc_d        =  .FALSE.
  lbdclim_d           =  .FALSE.
  lforest_d           =  .FALSE.
  lurban_d            =  .FALSE.
  lsso_d              =  .FALSE.
  lradtopo_d          =  .FALSE.
  nhori_d             =  24
  lseaice_d           =  .FALSE.
  llake_d             =  .FALSE.
  llake_coldstart_d   =  .FALSE.
  lemiss_d            =  .FALSE.
  lstomata_d          =  .FALSE.
  llbc_smooth_d       =  .FALSE.
  nlbc_smooth_d       =  20
  lroutine_d          =  .FALSE.
  yinput_model_d      =  '     '
  lasync_io_d         =  .FALSE.
  lreorder_d          =  .TRUE.
  ldatatypes_d        =  .FALSE.
  l_topo_z_d          =  .FALSE.
  itype_calendar_d    =  0
  itype_w_so_rel_d    =  1
  itype_t_cl_d        =  0 
  itype_rootdp_d      =  0 
  itype_ndvi_d        =  0 
  itype_aerosol_d     =  1 ! (constant values are used in COSMO Model)
  itype_albedo_d      =  1
  lmulti_layer_lm_d   =  .FALSE.
  lmulti_layer_in_d   =  .FALSE.
  idbg_level_d        =  2
  lprintdeb_all_d     = .FALSE.
  l_cressman_d        =  .FALSE.
  l_bicub_spl_d       =  .FALSE.
  l_art_d             =  .FALSE.
  l_art_nested_d      =  .FALSE.
! iso code
  liso_d              =  .FALSE.
! end iso code

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with default
!------------------------------------------------------------------------------
 
  ydate_ini           = ydate_ini_d
  ydate_bd            = ydate_bd_d
  hstart              = hstart_d
  hstop               = hstop_d
  hincbound           = hincbound_d
  eps_filter          = eps_filter_d
  rxso_mask           = rxso_mask_d
  rfill_valley        = rfill_valley_d
  nincwait            = nincwait_d
  nmaxwait            = nmaxwait_d
  ytrans_in           = ytrans_in_d
  ytrans_out          = ytrans_out_d
  nprocx              = nprocx_d
  nprocy              = nprocy_d
  nprocio             = nprocio_d
  nboundlines         = nboundlines_d
  ncomm_type          = ncomm_type_d
  norder_filter       = norder_filter_d
  ilow_pass_oro       = ilow_pass_oro_d
  numfilt_oro         = numfilt_oro_d
  ilow_pass_xso       = ilow_pass_xso_d
  numfilt_xso         = numfilt_xso_d
  ifill_valley        = ifill_valley_d
  luse_t_skin         = luse_t_skin_d
  lante_0006          = lante_0006_d
  lpost_0006          = lpost_0006_d
  luvcor              = luvcor_d
  lvertwind_ini       = lvertwind_ini_d
  lvertwind_bd        = lvertwind_bd_d
  lprog_qi            = lprog_qi_d
  l_smi               = l_smi_d
  lmixcld             = lmixcld_d
  lprog_qr_qs         = lprog_qr_qs_d
  lprog_qg            = lprog_qg_d
  lprog_rho_snow      = lprog_rho_snow_d
  qvmin               = qvmin_d
  qcmin               = qcmin_d
  qimin               = qimin_d
  linitial            = linitial_d
  lboundaries         = lboundaries_d
  leps_bc             = leps_bc_d
  lfilter_oro         = lfilter_oro_d
  lxso_first          = lxso_first_d
  lfilter_pp          = lfilter_pp_d
  lbalance_pp         = lbalance_pp_d
  lt_cl_corr          = lt_cl_corr_d
  ltime_mean          = ltime_mean_d
  ltime_proc          = ltime_proc_d
  lbdclim             = lbdclim_d
  lforest             = lforest_d
  lurban              = lurban_d
  lsso                = lsso_d
  lradtopo            = lradtopo_d
  nhori               = nhori_d
  lseaice             = lseaice_d
  llake               = llake_d
  llake_coldstart     = llake_coldstart_d
  lemiss              = lemiss_d
  lstomata            = lstomata_d
  llbc_smooth         = llbc_smooth_d
  nlbc_smooth         = nlbc_smooth_d
  lroutine            = lroutine_d
  yinput_model        = yinput_model_d
  lasync_io           = lasync_io_d
  lreorder            = lreorder_d
  ldatatypes          = ldatatypes_d
  l_topo_z            = l_topo_z_d
  itype_calendar      = itype_calendar_d
  itype_w_so_rel      = itype_w_so_rel_d
  itype_t_cl          = itype_t_cl_d
  itype_rootdp        = itype_rootdp_d
  itype_ndvi          = itype_ndvi_d
  itype_aerosol       = itype_aerosol_d
  itype_albedo        = itype_albedo_d
  lmulti_layer_lm     = lmulti_layer_lm_d
  lmulti_layer_in     = lmulti_layer_in_d
  idbg_level          = idbg_level_d
  lprintdeb_all       = lprintdeb_all_d
  l_cressman          = l_cressman_d
  l_bicub_spl         = l_bicub_spl_d
  l_art               = l_art_d
  l_art_nested        = l_art_nested_d
! iso code
  liso                = liso_d
! end iso code

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  READ (nuin, contrl, IOSTAT=ierror)
ENDIF

! distribute error status to all processors
IF (nproc > 1) THEN
  CALL distribute_values  (ierror, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF
IF (ierror /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

  ! Check whether an input model has been specified and set lxxx2lm
  ! (this is done here in proc 0, to check other settings; must be done
  !  in all other procs later on)
  lgme2lm = .FALSE.
  lgfs2lm = .FALSE.
  lgsm2lm = .FALSE.
  llm2lm  = .FALSE.
  lcm2lm  = .FALSE.
  lec2lm  = .FALSE.
  lum2lm  = .FALSE.
  lhir2lm = .FALSE.

  SELECT CASE (TRIM(yinput_model))
  CASE ('GME')
    lgme2lm = .TRUE.
  CASE ('GFS')
    lgfs2lm  = .TRUE.
  CASE ('GSM')
    lgsm2lm  = .TRUE.
  CASE ('COSMO')
    llm2lm  = .TRUE.
  CASE ('IFS')
    lec2lm  = .TRUE.
  CASE ('UMR', 'UMG')
    lum2lm  = .TRUE.
  CASE ('HIRLM')
    lhir2lm = .TRUE.
  CASE ('CM')
    lcm2lm  = .TRUE.
  CASE DEFAULT
    PRINT *, 'ERROR: *** No or wrong input model has been specified', yinput_model, ' ***'
    ierrstat = 1002
  END SELECT

  ! Check whether the values for start and end of the forecast are
  ! given in hours and calculate the values in time steps
  IF ( hstart /= hstart_d ) THEN
    ! only allow multiples of 0.25
    htest = hstart / 0.25_ireals
    IF ( ABS(NINT(htest) - htest) > 1E-5) THEN
      ! then it is not a multiple of 0.25
      PRINT *, 'ERROR: *** This is not a valid hstart: ', hstart, ' ***'
      PRINT *, '       *** only multiples of 0.25 are allowed       ***'
      ierrstat = 1002
    ENDIF
  ENDIF
  nstart    = NINT(hstart * 3600.0_ireals/dt)

  IF ( hstop /= hstop_d ) THEN
    ! only allow multiples of 0.25
    htest = hstop  / 0.25_ireals
    IF ( ABS(NINT(htest) - htest) > 1E-5) THEN
      ! then it is not a multiple of 0.25
      PRINT *, 'ERROR: *** This is not a valid hstop : ', hstop , ' ***'
      PRINT *, '       *** only multiples of 0.25 are allowed       ***'
      ierrstat = 1002
    ENDIF
  ENDIF
  nstop     = NINT(hstop * 3600.0_ireals/dt)

  IF ( hincbound /= hincbound_d) THEN
    IF ( ABS(NINT(hincbound) - hincbound) > 1E-5) THEN
      ! then it is not a full hour, only allow 0.25 and 0.5
      IF ( (hincbound /= 0.50_ireals) .AND. (hincbound /= 0.25_ireals) ) THEN
        PRINT *, 'ERROR: *** This is not a valid hincbound: ', hincbound, ' ***'
        PRINT *, '       *** only values = n.0 / 0.5 / 0.25 are allowed   ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF
  nincbound = NINT (hincbound * 3600.0_ireals / dt)

  ! Check whether the start date has been set, because this
  ! is needed:
  IF ( ydate_ini == ydate_ini_d ) THEN
    PRINT *,' ERROR   ***  ydate_ini not set ***'
    PRINT *,'         ***  Please specify ydate_ini in the format YYYYMMDDHH[MMSS]  ***'
    ierrstat = 1025
  ELSE
    ! Check, whether 10 or 14 digits are used for the date format
    ! NOTE: lmmss_ini must be distributed to other PEs
    IF     (LEN_TRIM(ydate_ini) == 10) THEN
      ydate_ini(11:14) = '0000'
      lmmss_ini = .FALSE.
      PRINT *, ' *** NOTE: Old 10 digit date format is used for output files of INT2LM'
    ELSEIF (LEN_TRIM(ydate_ini) == 14) THEN
      lmmss_ini = .TRUE.
      PRINT *, ' *** NOTE: New 14 digit date format is used for output files of INT2LM'
    ELSE
      PRINT *, ' *** ERROR: Wrong number of digits for ydate_ini! *** '
      PRINT *, ' ***        Must be 10 or 14, but are ', LEN_TRIM(ydate_ini)
      ierrstat = 1025
    ENDIF
  ENDIF

  ! Check whether a date is given for the start of the forecast for the
  ! boundary fields
  IF ( ydate_bd == ydate_bd_d ) THEN
    ! Use the same number of digits as ydate_ini
    IF (lmmss_ini) THEN
      ydate_bd(1:14) = ydate_ini(1:14)
      lmmss_bd  = .TRUE.
    ELSE
      ydate_bd(1:14) = ydate_ini(1:10)//'0000'
      lmmss_bd  = .FALSE.
    ENDIF
  ELSE
    ! Check, whether 10 or 14 digits are used for the boundary date format
    ! NOTE: lmmss_bd  must be distributed to other PEs
    IF     (LEN_TRIM(ydate_bd) == 10) THEN
      ydate_bd(11:14) = '0000'
      lmmss_bd = .FALSE.
      PRINT *, ' *** NOTE: Old 10 digit date format is used for input files of INT2LM'
    ELSEIF (LEN_TRIM(ydate_bd) == 14) THEN
      lmmss_bd = .TRUE.
      PRINT *, ' *** NOTE: New 14 digit date format is used for input files of INT2LM'
    ELSE
      PRINT *, ' *** ERROR: Wrong number of digits for ydate_bd! *** '
      PRINT *, ' ***        Must be 10 or 14, but are ', LEN_TRIM(ydate_bd)
      ierrstat = 1025
    ENDIF
  ENDIF

  ! Check whether linitial and lboundaries are not both .FALSE.
  IF ( (.NOT. linitial) .AND. (.NOT. lboundaries) ) THEN
    PRINT *, 'linitial and lboundaries = .FALSE. ?'
    ierrstat = 1001
  ENDIF

  ! Check nincwait and nmaxwait
  IF (nincwait < 0) THEN
    PRINT *, ' ERROR    *** nincwait < 0:  ', nincwait, '  *** '
    ierrstat = 1001
  ENDIF
  IF (nmaxwait < 0) THEN
    PRINT *, ' ERROR    *** nmaxwait < 0:  ', nmaxwait, '  *** '
    ierrstat = 1001
  ENDIF

  ! Check length of the directory-names
  nzylen = LEN_TRIM(ytrans_in)
  IF( nzylen > 0 ) THEN
    IF( ytrans_in(nzylen:nzylen) /= '/') THEN
      IF( nzylen < LEN(ytrans_in) ) THEN
        ytrans_in = ytrans_in(1:nzylen)//'/'
      ELSE
        PRINT *,' ERROR    *** ytrans_in is too long *** '
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF

  nzylen = LEN_TRIM(ytrans_out)
  IF( nzylen > 0 ) THEN
    IF( ytrans_out(nzylen:nzylen) /= '/') THEN
      IF( nzylen < LEN(ytrans_out) ) THEN
        ytrans_out = ytrans_out (1:nzylen)//'/'
      ELSE
        PRINT *,' ERROR    *** ytrans_out is too long *** '
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF

  ! Check whether nprocio and lasync_io do fit
  IF ( (nprocio /= 0) .AND. (lasync_io .EQV. .FALSE.) ) THEN
    PRINT *,' WARNING  *** For non-asynchronous IO nprocio = 0 is needed ***'
    nprocio   = 0
    lasync_io = .FALSE.
  ENDIF
  IF ( (nprocio == 0) .AND. (lasync_io .EQV. .TRUE.) ) THEN
    PRINT *,' WARNING  *** For nprocio == 0 synchronous IO is needed ***'
    lasync_io = .FALSE.
  ENDIF

  ! Check whether nprocx*nprocy + nprocio = nproc
  IF ( nprocx * nprocy + nprocio /= nproc ) THEN
    PRINT *,' ERROR    *** Wrong size of the processor-grid for nproc *** '
    PRINT *,'          *** ',nprocx,' * ',nprocy,' + ',nprocio,' /= ',nproc
    ierrstat = 1002
  ENDIF

  IF ( rfill_valley > 0.0_ireals ) THEN
    SELECT CASE( ifill_valley )
    CASE( 1:2 )
      ! okay
    CASE default
      PRINT *,' WARNING  *** To be relevant / do the filling of valleys *** '
      PRINT *,'          *** ifill_valley has to be between  1 and 2:   *** '
      PRINT *,'          *** 1 = MIN    filling BEFORE filter           *** '
      PRINT *,'          *** 2 = MIN    filling  AFTER filter           *** '
      PRINT *,'          *** set ifill_valley = 0   (default value) and *** '
      PRINT *,'          ***     rfill_valley = 0.0 (default value) !!! *** '
      ifill_valley = 0_iintegers
      rfill_valley = 0.0_ireals
    END SELECT
  END IF

  IF ( lfilter_oro ) THEN
    ! Check filter order
    IF ( (norder_filter < 2) .OR. (norder_filter > 5)) THEN
      PRINT *,' ERROR    *** Order of filter has to be  2, 3, 4 or 5 *** '
      ierrstat = 1002
    ENDIF
    ! Check setting of low-pass filter
    SELECT CASE( ilow_pass_oro )
    CASE( 1 )
      ! okay
      ! Default Raymond filter
    CASE( 3:6, 8 )
      ! okay
      ! Low-pass filter with predefined filter weights
      IF ( (numfilt_oro < 1) .OR. (numfilt_oro > 10) ) THEN
        PRINT *,' WARNING  *** numfilt_oro has to be between  1 and 10 *** '
        PRINT *,'          *** set numfilt_oro = 1 (default value)!    *** '
        numfilt_oro = 1_iintegers
      ENDIF
    CASE default
      PRINT *,' WARNING  *** ilow_pass_oro has to be 1, 3, 4, 5, 6 or 8 *** '
      PRINT *,'          *** set ilow_pass_oro = 1 (default value)!   *** '
      ilow_pass_oro = 1_iintegers
    END SELECT
    ! Check setting of low-pass filter for eXtra smoothing of steep orography
    IF ( ilow_pass_xso >= ilow_pass_oro ) THEN
      SELECT CASE( ilow_pass_xso )
      CASE( 3:6, 8 )
        ! okay
        ! Low-pass filter with predefined filter weights
        IF ( (numfilt_xso < 1) .OR. (numfilt_xso > 10) ) THEN
          PRINT *,' WARNING  *** numfilt_xso has to be between  1 and 10 *** '
          PRINT *,'          *** set numfilt_xso = 1 (default value)!    *** '
          numfilt_xso = 1_iintegers
        ENDIF
      CASE default
        PRINT *,' WARNING  *** ilow_pass_xso has to be 3, 4, 5, 6 or 8 *** '
        PRINT *,'          *** and to be  greater/equal ilow_pass_oro *** '
        PRINT *,'          *** set ilow_pass_xso = 0 (default value)! *** '
        ilow_pass_xso = 0_iintegers
      END SELECT
    ELSE
      IF ( ilow_pass_xso /= 0_iintegers ) THEN
        PRINT *,' WARNING  *** ilow_pass_xso has to be 3, 4, 5, 6 or 8 *** '
        PRINT *,'          *** and to be  greater/equal ilow_pass_oro *** '
        PRINT *,'          *** set ilow_pass_xso = 0 (default value)! *** '
        ilow_pass_xso = 0_iintegers
      ENDIF
    ENDIF
  ENDIF


  ! Check settings of soil model
  IF ((lmulti_layer_in) .AND. (.NOT. lmulti_layer_lm)) THEN
    PRINT *,' WARNING: *** if   lmulti_layer_in is .TRUE.,        ***'
    PRINT *,'        : *** also lmulti_layer_lm has to be .TRUE., ***'
    lmulti_layer_lm = .TRUE.
  ENDIF

  IF (lec2lm .AND. lmulti_layer_in .AND. .NOT.l_smi) THEN
    PRINT *," WARNING: *** if   yinput_model='IFS' and lmulti_layer_in is .TRUE.,        ***"
    PRINT *,"        : *** also l_smi has to be .TRUE.,                                  ***"
    l_smi = .TRUE.
  ENDIF

  IF (lcm2lm .AND. .NOT. (lmulti_layer_lm .AND. lmulti_layer_in)) THEN
    PRINT *," ERROR  *** soil layer specification for yinput_model='CM' "
    PRINT *,'  lmulti_layer_lm and lmulti_layer_in must be set .TRUE. *** '
    ierrstat = 1003
  ENDIF

  ! Check whether type of communication is in the correct range
  IF ( (ncomm_type < 1) .OR. (ncomm_type > 3) ) THEN
    PRINT *,' ERROR    *** unknown type of communication type ***'
    PRINT *,'          ***   (1 <= ncomm_type <= 3), but is ', ncomm_type,' ***'
    ierrstat = 1004
  ENDIF

  ! Check whether lfilter_oro and l_topo_z do match (l_topo_z may only be true,
  ! if also lfilter_oro is set)
  IF ( (.NOT. lfilter_oro) .AND. (l_topo_z) ) THEN
    PRINT *,' ERROR *** l_topo_z is .TRUE., but lfilter_oro is .FALSE.  ***'
    PRINT *,'       ***   but lfilter_oro has also to be .TRUE. then    ***'
    ierrstat = 1005
  ENDIF

  ! Check lprog_rho_snow (only possible for lgme2lm and llm2lm)
  IF ( lprog_rho_snow .AND. .NOT. (lgme2lm .OR. llm2lm)) THEN
    PRINT *, "ERROR *** lprog_rho_snow=.TRUE. is only possible for yinput_model='GME' ***"
    ierrstat = 1007
  ENDIF

  ! Check itype_calendar
  IF ( itype_calendar < 0 .OR. itype_calendar > 2) THEN
    PRINT *, 'ERROR *** itype_calendar has invalid value'
    PRINT *,'       ***   (0 <= itype_calendar <= 2), but is ', itype_calendar,' ***'
    ierrstat = 1007
  ENDIF

  ! Check itype_w_so_rel
  IF ( itype_w_so_rel < 0 .OR. itype_w_so_rel > 4) THEN
    PRINT *, 'ERROR *** itype_w_so_rel has invalid value'
    PRINT *,'       ***   (0 <= itype_w_so_rel <= 4), but is ', itype_w_so_rel,' ***'
    ierrstat = 1007
  ENDIF
! SP, 201405
  IF ( lcm2lm .AND. ((itype_w_so_rel == 1) .OR. (itype_w_so_rel == 2))) THEN
    PRINT *, 'ERROR *** itype_w_so_rel has invalid value'
    PRINT *,'       *** (1 or 2 not implemented for lcm2lm = .TRUE.)'
    ierrstat = 1007
  ENDIF

  ! Check itype_t_cl
  IF ( itype_t_cl < 0 .OR. itype_t_cl > 1) THEN
    PRINT *, 'ERROR *** itype_t_cl has invalid value'
    PRINT *,'       ***   (0 <= itype_t_cl     <= 1), but is ', itype_t_cl    ,' ***'
    ierrstat = 1007
  ENDIF
  IF (itype_t_cl == 1 .AND. llm2lm) THEN
    itype_t_cl = 0
    PRINT *, "WARNING *** itype_t_cl is ',itype_t_cl,' but has to be 0 in case yinput_model='COSMO'"
    PRINT *, '            itype_t_cl is set to 0  ***'
  ENDIF

  ! Check itype_rootdp 
  IF ( itype_rootdp  < 0 .OR. itype_rootdp  > 4) THEN
    PRINT *, 'ERROR *** itype_rootdp   has invalid value'
    PRINT *,'       ***   (0 <= itype_rootdp   <= 4), but is ', itype_rootdp  ,' ***'
    ierrstat = 1007
  ENDIF

  ! Check itype_ndvi
  IF ( itype_ndvi < 0 .OR. itype_ndvi > 2) THEN
    PRINT *, 'ERROR *** itype_ndvi has invalid value'
    PRINT *,'       ***   (0 <= itype_ndvi     <= 2), but is ', itype_ndvi    ,' ***'
    ierrstat = 1007
  ENDIF

  ! Check itype_aerosol
  IF ( itype_aerosol < 1 .OR. itype_aerosol > 2) THEN
    PRINT *, 'ERROR *** itype_aerosol has invalid value'
    PRINT *,'       ***   (1 <= itype_aerosol <= 2), but is ', itype_aerosol    ,' ***'
    ierrstat = 1007
  ENDIF

  ! Check itype_albedo
  IF ( itype_albedo  < 1 .OR. itype_albedo  > 4) THEN
    PRINT *, 'ERROR *** itype_albedo  has invalid value'
    PRINT *,'       ***   (1 <= itype_albedo  <= 4), but is ', itype_albedo     ,' ***'
    ierrstat = 1007
  ENDIF

  ! check consistency between input model and l_bicub_spl
  ! (bicubic spline interpolation only possible for lcm2lm=.TRUE., because of
  !  an IF (lcm2lm) statement in the routine interp_q_bs.
  !  Can that be eliminated????)
  IF ( l_bicub_spl .AND. .NOT. lcm2lm) THEN
    PRINT *, 'ERROR *** bicubic spline interpolation (l_bicub_spl=.true.) ***'
    PRINT *,"      *** not possible for yinput_model='CM')  ***"
    ierrstat = 1008
  ENDIF

  ! Check llake and llake_coldstart
  IF (.NOT. llake .AND. llake_coldstart) THEN
    PRINT *, 'ERROR *** if llake_coldstart is set, also llake must be set ***'
    ierrstat = 1008
  ENDIF

  ! Cressman not yet compatible
  IF (l_cressman) THEN
    PRINT *, 'WARNING *** Cressman interpolation not compatible with this version ***'
    PRINT *, '        *** l_cressman is reset to .FALSE.                          ***'
    l_cressman = .FALSE.
  ENDIF

  ! Check consistency of lmixcld and lprog_qi
  IF (lmixcld .AND. (.NOT. lprog_qi)) THEN
    PRINT *, 'ERROR *** if lmixcld is .TRUE., then also lprog_qi has to be .TRUE.'
    ierrstat = 1007
  ENDIF

! iso code
  ! Checks for 	compatibility of liso
  IF ( liso .AND. (.NOT. (lec2lm .OR. lcm2lm)) ) THEN
    PRINT *, 'ERROR  *** Isotope data processing only implemented for ***'
    PRINT *, '       *** lec2lm = .TRUE. or lcm2lm = .TRUE.           ***'
    ierrstat = 1007
  ENDIF
  IF ( liso .AND. (ylm_form_write == 'ncdf')) THEN
    PRINT *, 'ERROR  *** Isotope data processing not implemented for  ***'
    PRINT *, '       *** ylm_form_write = ncdf                        ***'
    ierrstat = 1007
  ENDIF
! end iso code

ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf  ( 1) = nstart
    intbuf  ( 2) = nstop
    intbuf  ( 3) = nincbound
    intbuf  ( 4) = nprocx
    intbuf  ( 5) = nprocy
    intbuf  ( 6) = nprocio
    intbuf  ( 7) = nboundlines
    intbuf  ( 8) = nincwait
    intbuf  ( 9) = nmaxwait
    intbuf  (10) = norder_filter
    intbuf  (12) = ncomm_type
    intbuf  (13) = ilow_pass_oro
    intbuf  (14) = numfilt_oro
    intbuf  (15) = ilow_pass_xso
    intbuf  (16) = numfilt_xso
    intbuf  (17) = ifill_valley
    intbuf  (18) = nlbc_smooth
    intbuf  (19) = idbg_level
    intbuf  (20) = itype_w_so_rel
    intbuf  (21) = itype_t_cl
    intbuf  (22) = itype_rootdp
    intbuf  (23) = itype_ndvi
    intbuf  (24) = itype_calendar
    intbuf  (25) = nhori
    intbuf  (26) = itype_aerosol
    intbuf  (27) = itype_albedo
    realbuf ( 1) = hstart
    realbuf ( 2) = hstop
    realbuf ( 3) = hincbound
    realbuf ( 5) = eps_filter
    realbuf ( 6) = qvmin
    realbuf ( 7) = qcmin
    realbuf ( 8) = qimin
    realbuf ( 9) = rxso_mask
    realbuf (10) = rfill_valley
    logbuf  ( 1) = luvcor
    logbuf  ( 2) = lvertwind_ini
    logbuf  ( 3) = lvertwind_bd
    logbuf  ( 4) = linitial
    logbuf  ( 5) = lboundaries
    logbuf  ( 6) = ltime_proc
    logbuf  ( 7) = ltime_mean
    logbuf  ( 9) = lroutine
    logbuf  (16) = lasync_io
    logbuf  (17) = lfilter_oro
    logbuf  (18) = lxso_first
    logbuf  (19) = lfilter_pp
    logbuf  (20) = lbalance_pp
    logbuf  (21) = lreorder
    logbuf  (22) = lmulti_layer_lm
    logbuf  (23) = lmulti_layer_in
    logbuf  (24) = lprog_qi
    logbuf  (25) = lt_cl_corr
    logbuf  (26) = luse_t_skin
    logbuf  (27) = lante_0006
    logbuf  (28) = lpost_0006
    logbuf  (29) = ldatatypes
    logbuf  (30) = l_topo_z
    logbuf  (31) = lforest
    logbuf  (32) = lprog_rho_snow
    logbuf  (33) = lprog_qr_qs
    logbuf  (34) = lprog_qg
    logbuf  (35) = llbc_smooth
    logbuf  (36) = llake
    logbuf  (38) = lbdclim
    logbuf  (39) = lprintdeb_all
    logbuf  (40) = l_cressman
    logbuf  (41) = l_bicub_spl
    logbuf  (42) = l_art
    logbuf  (43) = leps_bc
    logbuf  (44) = lsso
    logbuf  (45) = lradtopo
    logbuf  (46) = lstomata
    logbuf  (47) = lemiss
    logbuf  (48) = l_smi
    logbuf  (49) = lmixcld
    logbuf  (50) = lum2lm
    logbuf  (51) = lseaice
    logbuf  (52) = llake_coldstart
    logbuf  (53) = lurban
    logbuf  (54) = lmmss_ini
    logbuf  (55) = lmmss_bd
    logbuf  (56) = l_art_nested
! iso code
    logbuf  (57) = liso
! end iso code
    charbuf ( 1) = ydate_ini
    charbuf ( 2) = ydate_bd
    charbuf ( 3)(  1:100) = ytrans_in (  1:100)
    charbuf ( 4)(  1:100) = ytrans_in (101:200)
    charbuf ( 5)(  1: 50) = ytrans_in (201:250)
    charbuf ( 6)(  1:100) = ytrans_out(  1:100)
    charbuf ( 7)(  1:100) = ytrans_out(101:200)
    charbuf ( 8)(  1: 50) = ytrans_out(201:250)
    charbuf ( 9) = yinput_model
  ENDIF

  CALL distribute_values  (intbuf ,27, 0, imp_integers,  icomm_world, ierr)
  CALL distribute_values  (realbuf,10, 0, imp_reals,     icomm_world, ierr)
! iso code: changed value
  CALL distribute_values  (logbuf ,57, 0, imp_logical,   icomm_world, ierr)
  CALL distribute_values  (charbuf, 9, 0, imp_character, icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    nstart        = intbuf  ( 1)
    nstop         = intbuf  ( 2)
    nincbound     = intbuf  ( 3)
    nprocx        = intbuf  ( 4)
    nprocy        = intbuf  ( 5)
    nprocio       = intbuf  ( 6)
    nboundlines   = intbuf  ( 7)
    nincwait      = intbuf  ( 8)
    nmaxwait      = intbuf  ( 9)
    norder_filter = intbuf  (10)
    ncomm_type    = intbuf  (12)
    ilow_pass_oro = intbuf  (13)
    numfilt_oro   = intbuf  (14)
    ilow_pass_xso = intbuf  (15)
    numfilt_xso   = intbuf  (16)
    ifill_valley  = intbuf  (17)
    nlbc_smooth   = intbuf  (18)
    idbg_level    = intbuf  (19)
    itype_w_so_rel= intbuf  (20)
    itype_t_cl    = intbuf  (21)
    itype_rootdp  = intbuf  (22)
    itype_ndvi    = intbuf  (23)
    itype_calendar= intbuf  (24)
    nhori         = intbuf  (25)
    itype_aerosol = intbuf  (26)
    itype_albedo  = intbuf  (27)
    hstart        = realbuf ( 1)
    hstop         = realbuf ( 2)
    hincbound     = realbuf ( 3)
    eps_filter    = realbuf ( 5)
    qvmin         = realbuf ( 6)
    qcmin         = realbuf ( 7)
    qimin         = realbuf ( 8)
    rxso_mask     = realbuf ( 9)
    rfill_valley  = realbuf (10)
    luvcor        = logbuf  ( 1)
    lvertwind_ini = logbuf  ( 2)
    lvertwind_bd  = logbuf  ( 3)
    linitial      = logbuf  ( 4)
    lboundaries   = logbuf  ( 5)
    ltime_proc    = logbuf  ( 6)
    ltime_mean    = logbuf  ( 7)
    lroutine      = logbuf  ( 9)
    lasync_io     = logbuf  (16)
    lfilter_oro   = logbuf  (17)
    lxso_first    = logbuf  (18)
    lfilter_pp    = logbuf  (19)
    lbalance_pp   = logbuf  (20)
    lreorder      = logbuf  (21)
    lmulti_layer_lm = logbuf(22)
    lmulti_layer_in = logbuf(23)
    lprog_qi      = logbuf  (24)
    lt_cl_corr    = logbuf  (25)
    luse_t_skin   = logbuf  (26)
    lante_0006    = logbuf  (27)
    lpost_0006    = logbuf  (28)
    ldatatypes    = logbuf  (29)
    l_topo_z      = logbuf  (30)
    lforest       = logbuf  (31)
    lprog_rho_snow= logbuf  (32)
    lprog_qr_qs   = logbuf  (33)
    lprog_qg      = logbuf  (34)
    llbc_smooth   = logbuf  (35)
    llake         = logbuf  (36)
    lbdclim       = logbuf  (38)
    lprintdeb_all = logbuf  (39)
    l_cressman    = logbuf  (40)
    l_bicub_spl   = logbuf  (41)
    l_art         = logbuf  (42)
    leps_bc       = logbuf  (43)
    lsso          = logbuf  (44)
    lradtopo      = logbuf  (45)
    lstomata      = logbuf  (46)
    lemiss        = logbuf  (47)
    l_smi         = logbuf  (48)
    lmixcld       = logbuf  (49)
    lum2lm        = logbuf  (50)
    lseaice       = logbuf  (51)
    llake_coldstart=logbuf  (52)
    lurban        = logbuf  (53)
    lmmss_ini     = logbuf  (54)
    lmmss_bd      = logbuf  (55)
    l_art_nested  = logbuf  (56)
! iso code
    liso          = logbuf  (57)
! end iso code
    ydate_ini     = charbuf ( 1)
    ydate_bd      = charbuf ( 2)
    ytrans_in (  1:100) = charbuf ( 3)(  1:100)
    ytrans_in (101:200) = charbuf ( 4)(  1:100)
    ytrans_in (201:250) = charbuf ( 5)(  1:100)
    ytrans_out(  1:100) = charbuf ( 6)(  1:100)
    ytrans_out(101:200) = charbuf ( 7)(  1:100)
    ytrans_out(201:250) = charbuf ( 8)(  1:100)
    yinput_model  = charbuf ( 9)
  ENDIF

  ! Check whether an input model has been specified and set lxxx2lm
  lgme2lm = .FALSE.
  lgfs2lm = .FALSE.
  lgsm2lm = .FALSE.
  llm2lm  = .FALSE.
  lcm2lm  = .FALSE.
  lec2lm  = .FALSE.
  lum2lm  = .FALSE.
  lhir2lm = .FALSE.

  SELECT CASE (yinput_model)
  CASE ('GME')
    lgme2lm = .TRUE.
  CASE ('GFS')
    lgfs2lm  = .TRUE.
  CASE ('GSM')
    lgsm2lm  = .TRUE.
  CASE ('COSMO')
    llm2lm  = .TRUE.
  CASE ('IFS')
    lec2lm  = .TRUE.
  CASE ('UMR', 'UMG')
    lum2lm  = .TRUE.
  CASE ('HIRLM')
    lhir2lm = .TRUE.
  CASE ('CM')
    lcm2lm  = .TRUE.
  CASE DEFAULT
    PRINT *, 'ERROR: *** No or wrong input model has been specified', yinput_model, ' ***'
    ierrstat = 1002
  END SELECT

ENDIF

  ! Set ltime depending on ltime_mean and ltime_proc
  IF (lclock .EQV. .FALSE.) THEN
    ! If no system clock is present, ltime has to be .FALSE.
    ltime = .FALSE.
    IF (my_world_id == 0) THEN
      PRINT *,'  WARNING  ***  NO SYSTEM CLOCK PRESENT ***'
      PRINT *,'           ***  ltime IS SET TO .FALSE. ***'
    ENDIF
  ELSE
    ltime = ltime_mean .OR. ltime_proc
    IF (ltime_mean .AND. ltime_proc) THEN
      ltime_proc = .FALSE.
      IF (my_world_id == 0) THEN
        PRINT *,'WARNING  *** ltime_mean and ltime_proc cannot both be set ***'
        PRINT *,'         *** ltime_mean is set ***'
      ENDIF
    ENDIF
  ENDIF

! Allocate space for timings
  IF (ltime .EQV. .TRUE.) THEN
    ALLOCATE ( timings (40), STAT=ierr )
    timings (:) = 0.0_ireals
  ENDIF

  ! Determine number of compute PEs and number of IO PEs
  num_compute = nprocx * nprocy
  IF (lasync_io) THEN
    num_io = nprocio
  ELSE
    num_io = 1
  ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN
  ! define time and date
  CALL DATE_AND_TIME (ydatearg, ytimearg)
  ydate = ydatearg(7:8)//'.'//ydatearg(5:6)//'.'//ydatearg(1:4)
  ytime = ytimearg(1:2)//'.'//ytimearg(3:4)//'.'//ytimearg(5:6)

  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(A23)') '      NAMELIST:  contrl'
  WRITE (nuout, '(A23)') '      -----------------'
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value',       &
                                               'Default Value', 'Format'

  WRITE (nuout, '(T8,A,T21,  A  ,T40,  A  ,T59,A3)')                          &
                          'yinput_model', yinput_model, yinput_model_d,'C* 5'
  IF (lmmss_ini) THEN
    WRITE (nuout, '(T8,A,T21,  A  ,T40,  A  ,T59,A3)')                          &
                     'ydate_ini ' ,ydate_ini(1:14), ydate_ini_d(1:14), 'C*14'
  ELSE
    WRITE (nuout, '(T8,A,T21,  A  ,T40,  A  ,T59,A3)')                          &
                     'ydate_ini ' ,ydate_ini(1:10), ydate_ini_d(1:10), 'C*10'
  ENDIF
  IF (lmmss_bd) THEN
    WRITE (nuout, '(T8,A,T21,  A  ,T40,  A  ,T59,A3)')                          &
                     'ydate_bd  ' ,ydate_bd (1:14), ydate_bd_d (1:14), 'C*14'
  ELSE
    WRITE (nuout, '(T8,A,T21,  A  ,T40,  A  ,T59,A3)')                          &
                     'ydate_bd  ' ,ydate_bd (1:10), ydate_bd_d (1:10), 'C*10'
  ENDIF
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                           'nprocx', nprocx ,nprocx_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                           'nprocy', nprocy ,nprocy_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                        'nprocio', nprocio, nprocio_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                            'nboundlines', nboundlines, nboundlines_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                               'ncomm_type', ncomm_type, ncomm_type_d ,' I '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                           'hstart ',hstart ,hstart_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                           'hstop  ',hstop  ,hstop_d  ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                'hincbound', hincbound, hincbound_d   ,' R '
  WRITE (nuout, '(T8,A,T25,I8   ,T40,I12  ,T59,A3)')                          &
                   'itype_calendar', itype_calendar, itype_calendar_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                     'nincwait', nincwait, nincwait_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                     'nmaxwait', nmaxwait, nmaxwait_d ,' I '
  WRITE (nuout, '(T8,A,T21,  A  ,          T59,A5)')                          &
                                 'ytrans_in ', TRIM(ytrans_in),        'C*250'
  WRITE (nuout, '(T8,A,T21,  A  ,          T59,A5)')                          &
                                 'ytrans_out', TRIM(ytrans_out),       'C*250'
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                              'lante_0006',lante_0006, lante_0006_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                              'lpost_0006',lpost_0006, lpost_0006_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                           'luse_t_skin',luse_t_skin, luse_t_skin_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                          'luvcor',luvcor, luvcor_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                     'lvertwind_ini',lvertwind_ini, lvertwind_ini_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                        'lvertwind_bd',lvertwind_bd, lvertwind_bd_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprog_qi ',lprog_qi , lprog_qi_d    ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lmixcld ' , lmixcld , lmixcld_d     ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                              'l_smi', l_smi, l_smi_d                 ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                        'lprog_qr_qs ',lprog_qr_qs , lprog_qr_qs_d    ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprog_qg ',lprog_qg , lprog_qg_d    ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
               'lprog_rho_snow ',lprog_rho_snow , lprog_rho_snow_d    ,' L '
  WRITE (nuout, '(T8,A,T21,G12.4,T40,G12.4,T59,A3)')                          &
                                           'qvmin  ',qvmin  ,qvmin_d  ,' R '
  WRITE (nuout, '(T8,A,T21,G12.4,T40,G12.4,T59,A3)')                          &
                                           'qcmin  ',qcmin  ,qcmin_d  ,' R '
  WRITE (nuout, '(T8,A,T21,G12.4,T40,G12.4,T59,A3)')                          &
                                           'qimin  ',qimin  ,qimin_d  ,' R '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                    'linitial',linitial, linitial_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                           'lboundaries',lboundaries, lboundaries_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'leps_bc',leps_bc, leps_bc_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                           'ltime_mean ',ltime_mean , ltime_mean_d    ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                           'ltime_proc ',ltime_proc , ltime_proc_d    ,' L '
  WRITE (nuout, '(T8,A,T25,I12  ,T40,I12  ,T59,A3)')                          &
                  'itype_w_so_rel',itype_w_so_rel, itype_w_so_rel_d   ,' I '
  WRITE (nuout, '(T8,A,T25,I12  ,T40,I12  ,T59,A3)')                          &
                  'itype_t_cl    ',itype_t_cl    , itype_t_cl_d       ,' I '
  WRITE (nuout, '(T8,A,T25,I12  ,T40,I12  ,T59,A3)')                          &
                  'itype_rootdp  ',itype_rootdp  , itype_rootdp_d     ,' I '
  WRITE (nuout, '(T8,A,T25,I12  ,T40,I12  ,T59,A3)')                          &
                  'itype_ndvi    ',itype_ndvi    , itype_ndvi_d       ,' I '
  WRITE (nuout, '(T8,A,T25,I12  ,T40,I12  ,T59,A3)')                          &
            'itype_aerosol', itype_aerosol    , itype_aerosol_d       ,' I '
  WRITE (nuout, '(T8,A,T25,I12  ,T40,I12  ,T59,A3)')                          &
            'itype_albedo',  itype_albedo    , itype_albedo_d         ,' I '   
  WRITE (nuout, '(T8,A,T23,L10  ,T40,L12  ,T59,A3)')                          &
                'lmulti_layer_lm ',lmulti_layer_lm,lmulti_layer_lm_d  ,' L '
  WRITE (nuout, '(T8,A,T23,L10  ,T40,L12  ,T59,A3)')                          &
                'lmulti_layer_in ',lmulti_layer_in,lmulti_layer_in_d  ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                    'lbdclim' ,lbdclim,  lbdclim_d    ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                    'lforest' ,lforest,  lforest_d    ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                    'lurban' ,lurban,  lurban_d       ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                    'lsso',    lsso,        lsso_d    ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lradtopo' ,lradtopo,  lradtopo_d    ,' L '
  WRITE (nuout, '(T8,A,T25,I12  ,T40,I12  ,T59,A3)')                          &
                  'nhori',nhori    , nhori_d       ,' I '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                    'lseaice' ,lseaice,  lseaice_d    ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                      'llake' ,llake,  llake_d        ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
        'llake_coldstart' ,llake_coldstart,  llake_coldstart_d        ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                             'lstomata' ,lstomata,  lstomata_d        ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                   'lemiss' ,lemiss,  lemiss_d        ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                               'llbc_smooth',llbc_smooth,llbc_smooth_d,' L '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                             'nlbc_smooth', nlbc_smooth, nlbc_smooth_d,' I '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                    'lroutine',lroutine, lroutine_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lasync_io',lasync_io, lasync_io_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lreorder ',lreorder , lreorder_d    ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                             'ldatatypes',ldatatypes, ldatatypes_d    ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                             'lt_cl_corr', lt_cl_corr, lt_cl_corr_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                           'lfilter_oro',lfilter_oro, lfilter_oro_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                             'lxso_first', lxso_first, lxso_first_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                              'lfilter_pp',lfilter_pp, lfilter_pp_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                           'lbalance_pp',lbalance_pp, lbalance_pp_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                    'l_topo_z',l_topo_z, l_topo_z_d   ,' L '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                       'norder_filter', norder_filter,norder_filter_d ,' I '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                             'eps_filter', eps_filter, eps_filter_d   ,' R '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                       'ilow_pass_oro', ilow_pass_oro,ilow_pass_oro_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                          'numfilt_oro', numfilt_oro, numfilt_oro_d   ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                       'ilow_pass_xso', ilow_pass_xso,ilow_pass_xso_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                          'numfilt_xso', numfilt_xso, numfilt_xso_d   ,' I '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                'rxso_mask', rxso_mask, rxso_mask_d   ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                       'rfill_valley', rfill_valley, rfill_valley_d   ,' R '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                       'ifill_valley', ifill_valley, ifill_valley_d   ,' I '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                              'l_cressman',l_cressman, l_cressman_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                           'l_bicub_spl',l_bicub_spl, l_bicub_spl_d   ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                           'l_art',    l_art,     l_art_d             ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                       'l_art_nested', l_art_nested, l_art_nested_d,   ' L '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                 'idbg_level', idbg_level, idbg_level, ' I '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                       'lprintdeb_all', lprintdeb_all, lprintdeb_all_d,' L '
! iso code
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'liso', liso, liso_d                 ,' L '
! end iso code
  WRITE (nuout, '(A2)')  '  '
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_contrl

!==============================================================================
!+ Module procedure in "setup" for the input of NAMELIST grid_in
!------------------------------------------------------------------------------

SUBROUTINE input_grid_in (realbuf, intbuf, logbuf, charbuf, ibuflen,          &
                          nuout, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group grid_in. 
!   The group grid_in contains variables defining the resolution and vertical
!   size of the coarse input grid (GME, LM or IFS).
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for 
!   consistency. If wrong input values are detected the program prints 
!   an error message. The program is not stopped in this routine but an 
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists. 
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.    
!   Both, default and input values are written to the file OUTPUT.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    ibuflen,      & ! dimension of the buffers
    nuout,        & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Sending buffers (are given to routine distribute_values)
  INTEGER (KIND=iintegers), INTENT(INOUT)    ::  intbuf  (ibuflen)
  REAL    (KIND=ireals)   , INTENT(INOUT)    ::  realbuf (ibuflen)
  LOGICAL                 , INTENT(INOUT)    ::  logbuf  (ibuflen)
  CHARACTER (LEN=100)     , INTENT(INOUT)    ::  charbuf (ibuflen)

! Local variables

! Variables for default values
  INTEGER (KIND=iintegers)   ::       &
    ni_gme_d,   & ! resolution of GME
    i3e_gme_d,  & ! number of levels in the vertical
    kcontrol_fi_d ! control level for geopotential

  INTEGER (KIND=iintegers)   ::       &
    i, invar, ierr, ierror

  INTEGER (KIND=iintegers)    ::                                   &
    ie_in_tot_d,    & ! ie for input grid, total domain
    je_in_tot_d,    & ! je for input grid, total domain
    ke_in_tot_d,    & ! ke for input grid, total domain
    ke_hybrid_d,    & ! number of hybrid levels when interpolating from pressure levels
    nlevskip_d,     & ! default number of missing levels in input grid
    ke_soil_in_d,   & ! number of layers in multi-layer soil model in input
    east_add_in_d,  & ! add an extra column to the East
    west_add_in_d,  & ! add an extra column to the West
    south_add_in_d, & ! add an extra column to the South
    north_add_in_d    ! add an extra column to the North

  REAL (KIND=ireals)          ::                                   &
    pollat_in_d,      & ! latitude of the rotated north pole (in degrees, N>0)
    pollon_in_d,      & ! longitude of the rotated north pole (in degrees, E>0)
    polgam_in_d,      & ! latitude of the rotated north pole
    dlat_in_d,        & ! grid point distance in zonal direction (in degrees)
    dlon_in_d,        & ! grid point distance in meridional direction
    startlat_in_tot_d,& ! transformed latitude of the lower left grid point
                        ! of the total domain (in degrees, N>0)
    startlon_in_tot_d,& ! transformed longitude of the lower left grid point
                        ! of the total domain (in degrees, E>0)
    endlat_in_tot_d,  & ! transformed latitude of the upper right grid point
                        ! of the total domain (in degrees, N>0)
    endlon_in_tot_d,  & ! transformed longitude of the upper right grid point
                        ! of the total domain (in degrees, E>0)
    pcontrol_fi_d,    & ! pressure of control level for geopotential
    press_level_d(50)   ! list of available pressure levels (in Pa) for GFS

  REAL (KIND=ireals)          ::                                   &
    czml_soil_in_d(20), & ! depth of bottom level of soil layers (defaults)
    czml_soil_in  (20)    ! depth of bottom level of soil layers (read in)

  LOGICAL   :: &
    lushift_in_d(2),  & ! shift of u-velocity due to grid staggering
    lvshift_in_d(2)     ! shift of v-velocity due to grid staggering

  CHARACTER (LEN=15)         :: yroutine

! Define the namelist group
  NAMELIST /grid_in/ &
    ni_gme,         & ! resolution of GME
    i3e_gme,        & ! number of levels in the vertical
    kcontrol_fi,    & ! control level for geopotential
    ie_in_tot,      & ! ie for input grid, total domain
    je_in_tot,      & ! je for input grid, total domain
    ke_in_tot,      & ! ke for input grid, total domain
    ke_hybrid,      & ! number of pressure levels (for GFS)
    press_level,    & ! list of available pressure levels (in Pa) for GFS
    nlevskip,       & ! number of missing levels in input grid
    pcontrol_fi,    & ! pressure of control level for geopotential
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
    ke_soil_in,     & ! number of layers in multi-layer soil model in output
    czml_soil_in,   & ! depth of bottom level of soil layers (read in)
    lushift_in,     & ! shift of u-velocity due to grid staggering
    lvshift_in,     & ! shift of v-velocity due to grid staggering
    east_add_in,    & ! add an extra column to the East
    west_add_in,    & ! add an extra column to the West
    south_add_in,   & ! add an extra column to the South
    north_add_in      ! add an extra column to the North

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_grid_in
!------------------------------------------------------------------------------

yroutine = 'input_grid_in'
ierror   = 0

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  ni_gme_d          = 128_iintegers
  i3e_gme_d         =  31_iintegers
  kcontrol_fi_d     =  -1_iintegers
  ie_in_tot_d       = 141_iintegers
  je_in_tot_d       =  71_iintegers
  ke_in_tot_d       =  60_iintegers
  ke_hybrid_d       =   0_iintegers
  nlevskip_d        =   0_iintegers
  pcontrol_fi_d     =  -1.0_ireals
  pollat_in_d       =  90.0_ireals
  pollon_in_d       = 180.0_ireals
  polgam_in_d       =   0.0_ireals
  dlat_in_d         =   0.5_ireals
  dlon_in_d         =   0.5_ireals
  startlat_in_tot_d =  35.0_ireals
  startlon_in_tot_d = -30.0_ireals
  endlat_in_tot_d   =  0.0_ireals
  endlon_in_tot_d   = -40.0_ireals
  lushift_in_d(:)   = .FALSE.
  lvshift_in_d(:)   = .FALSE.
  east_add_in_d     = 0_iintegers
  west_add_in_d     = 0_iintegers
  south_add_in_d    = 0_iintegers
  north_add_in_d    = 0_iintegers
  ke_soil_in_d        = 7_iintegers
  czml_soil_in_d(1) =  0.005_ireals
  czml_soil_in_d(2) =  0.02_ireals
  czml_soil_in_d(3) =  0.06_ireals
  czml_soil_in_d(4) =  0.18_ireals
  czml_soil_in_d(5) =  0.54_ireals
  czml_soil_in_d(6) =  1.62_ireals
  czml_soil_in_d(7) =  4.86_ireals
  czml_soil_in_d(8) = 14.58_ireals

  press_level_d( 1) =   1000.0_ireals
  press_level_d( 2) =   2000.0_ireals
  press_level_d( 3) =   3000.0_ireals
  press_level_d( 4) =   5000.0_ireals
  press_level_d( 5) =   7000.0_ireals
  press_level_d( 6) =  10000.0_ireals
  press_level_d( 7) =  15000.0_ireals
  press_level_d( 8) =  20000.0_ireals
  press_level_d( 9) =  25000.0_ireals
  press_level_d(10) =  30000.0_ireals
  press_level_d(11) =  35000.0_ireals
  press_level_d(12) =  40000.0_ireals
  press_level_d(13) =  45000.0_ireals
  press_level_d(14) =  50000.0_ireals
  press_level_d(15) =  55000.0_ireals
  press_level_d(16) =  60000.0_ireals
  press_level_d(17) =  65000.0_ireals
  press_level_d(18) =  70000.0_ireals
  press_level_d(19) =  75000.0_ireals
  press_level_d(20) =  80000.0_ireals
  press_level_d(21) =  85000.0_ireals
  press_level_d(22) =  90000.0_ireals
  press_level_d(23) =  92500.0_ireals
  press_level_d(24) =  95000.0_ireals
  press_level_d(25) =  97500.0_ireals
  press_level_d(26) = 100000.0_ireals
  press_level_d(27:50) =   0.0_ireals

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with default
!------------------------------------------------------------------------------

  ni_gme            = ni_gme_d
  i3e_gme           = i3e_gme_d
  kcontrol_fi       = kcontrol_fi_d
  ie_in_tot         = ie_in_tot_d
  je_in_tot         = je_in_tot_d
  ke_in_tot         = ke_in_tot_d
  ke_hybrid         = ke_hybrid_d
  nlevskip          = nlevskip_d
  pcontrol_fi       = pcontrol_fi_d
  pollat_in         = pollat_in_d
  pollon_in         = pollon_in_d
  polgam_in         = polgam_in_d
  dlat_in           = dlat_in_d
  dlon_in           = dlon_in_d
  startlat_in_tot   = startlat_in_tot_d
  startlon_in_tot   = startlon_in_tot_d
  endlat_in_tot     = endlat_in_tot_d
  endlon_in_tot     = endlon_in_tot_d
  lushift_in(:)     = lushift_in_d(:)
  lvshift_in(:)     = lvshift_in_d(:)
  east_add_in       = east_add_in_d
  west_add_in       = west_add_in_d
  south_add_in      = south_add_in_d
  north_add_in      = north_add_in_d
  ke_soil_in        = ke_soil_in_d
  czml_soil_in(:)   = -1.0_ireals
  press_level(:)    = press_level_d(:)

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  READ (nuin, grid_in, IOSTAT=ierror)
ENDIF

! distribute error status to all processors
IF (nproc > 1) THEN
  CALL distribute_values  (ierror, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF
IF (ierror /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

  ! Check the control level
  IF (lec2lm) THEN
    IF (kcontrol_fi > 0_iintegers) THEN
      ! control level for GME must not be set for IFS2LM
      PRINT *, ' WARNING  *** kcontrol_fi must not be set for IFS2LM ***'
      kcontrol_fi = -1_iintegers
    ENDIF
    IF     (pcontrol_fi < 0.0_ireals) THEN
      ! control level pressure value is not given for IFS2LM
      PRINT *, ' ERROR    *** pcontrol_fi is not set for IFS2LM ***'
      ierrstat = 2001
    ELSEIF ((pcontrol_fi <  1000.0_ireals)                          &
       .OR. (pcontrol_fi > 99000.0_ireals)) THEN
      ! wrong values for pcontrol_fi
      PRINT *, ' ERROR    *** pcontrol_fi (',pcontrol_fi,           &
                                ') not in 1000.0....99000.0 *** '
      ierrstat = 2001
    ENDIF
  ENDIF

  IF (lgme2lm) THEN
    IF (pcontrol_fi > 0.0_ireals) THEN
      ! control level for IFS must not be set for GME2LM
      PRINT *, " WARNING  *** pcontrol_fi must not be set for yinput_model='GME' ***"
      pcontrol_fi = -1.0_ireals
    ENDIF
    IF     (kcontrol_fi < 0_iintegers) THEN
      ! pressure control level is not given for GME2LM
      PRINT *, ' ERROR    *** kcontrol_fi is not set for GME2LM ***'
      ierrstat = 2001
    ELSEIF (kcontrol_fi < 1 .OR. kcontrol_fi > i3e_gme) THEN
      ! wrong values for kcontrol_fi
      PRINT *, ' ERROR    *** kcontrol_fi',kcontrol_fi,             &
                                ' not in 1..',i3e_gme,' *** '
      ierrstat = 2001
    ENDIF
  ENDIF

  IF (lcm2lm) THEN
    IF (kcontrol_fi > 0_iintegers) THEN
      ! control level for GME must not be set for CM input
      PRINT *, " WARNING  *** kcontrol_fi must not be set for CM input (yinput_model='CM') ***"
      kcontrol_fi = -1_iintegers
    ENDIF

    IF(east_add_in  > 1_iintegers) THEN
      ! only one additional eastern boundary column for global input grids allowed so far
      PRINT *, ' WARNING  *** east_add_in  = ',east_add_in,            &
                                ' restricted to 1 ***'
      east_add_in  = 1_iintegers
    ENDIF
    IF(west_add_in  > 1_iintegers) THEN
      ! only one additional western boundary column for global input grids allowed so far
      PRINT *, ' WARNING  *** west_add_in  = ',west_add_in,            &
                                ' restricted to 1 ***'
      west_add_in  = 1_iintegers
    ENDIF
    IF(south_add_in > 1_iintegers) THEN
      ! only one additional southern boundary row for global input grids allowed so far
      PRINT *, ' WARNING  *** south_add_in = ',south_add_in,           &
                                ' restricted to 1 ***'
      south_add_in = 1_iintegers
    ENDIF
    IF(north_add_in > 1_iintegers) THEN
      ! only one additional northern boundary row for global input grids allowed so far
      PRINT *, ' WARNING  *** north_add_in = ',north_add_in,           &
                                ' restricted to 1 ***'
      north_add_in = 1_iintegers
    ENDIF

    ie_in_tot = ie_in_tot + east_add_in
    ie_in_tot = ie_in_tot + west_add_in
    je_in_tot = je_in_tot + south_add_in
    je_in_tot = je_in_tot + north_add_in

  ENDIF

  IF (lmulti_layer_in) THEN
    IF (lec2lm) THEN
      IF (ke_soil_in /= 3) THEN
        ke_soil_in = 3
        PRINT *, ' WARNING  *** ke_soil_in set unconditionally = 3 for IFS input'
      ENDIF
      ALLOCATE (czmls_in(1:ke_soil_in+1), STAT=ierr)
      ALLOCATE (czhls_in(0:ke_soil_in+1), STAT=ierr)
      ALLOCATE (msoilgrib_in(0:ke_soil_in+1), STAT=ierr)

      ! fixed IFS soil depths
      czmls_in(:) = (/0.035_ireals, 0.175_ireals, 0.64_ireals, 2.5_ireals/)

      ! this is probably useless, but we do it for completeness
      DO i = 1, ke_soil_in+1
        msoilgrib_in(i) = NINT (100.0_ireals * czmls_in(i)+1.0E-7_ireals)
      ENDDO

    ELSE ! all the other models

      ALLOCATE (czmls_in(1:ke_soil_in+1), STAT=ierr)
      ALLOCATE (czhls_in(0:ke_soil_in+1), STAT=ierr)  !_br
      ALLOCATE (msoilgrib_in(0:ke_soil_in+1), STAT=ierr)
       ! Check, how many soil levels have been specified
      invar = COUNT(czml_soil_in(:) /= -1.0_ireals)
      IF (invar == 0) THEN
        ! no level specifications have been read
        IF (ke_soil_in == ke_soil_in_d) THEN
          ! use the default
          PRINT *,' *** Default specifications of input soil main levels are used *** '
          czmls_in(1:ke_soil_in+1) = czml_soil_in_d(1:ke_soil_in+1)
        ELSEIF (ke_soil_in /= ke_soil_in_d) THEN
          PRINT *,' ERROR  *** no specifications of input soil levels, *** '
          PRINT *,' ERROR  *** but using wrong default                              *** ', &
                    ke_soil_in_d, ke_soil_in
          ierrstat = 1002
        ENDIF
      ELSE
        IF (ke_soil_in+1 == invar) THEN
          ! use specifications read in
          PRINT *,'  *** specifications of input soil main levels *** '
          PRINT *,'  *** from Namelist INPUT are used       *** '
          czmls_in(1:ke_soil_in+1) = czml_soil_in(1:ke_soil_in+1)
        ELSE
          PRINT *,' ERROR  *** wrong number of specifications ',           &
                  'for input soil levels  *** ', ke_soil_in, invar
  
          ierrstat = 1002
        ENDIF
      ENDIF
  
      IF (ierrstat == 0) THEN
        ! compute grib coded values of the depth of main soil levels
        msoilgrib_in(0) = 0_iintegers
        DO i = 1, ke_soil_in+1
          msoilgrib_in(i) = NINT (100.0_ireals * czmls_in(i)+1.0E-7_ireals)
        ENDDO
      ENDIF

    ENDIF
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf  ( 1) = ni_gme
    intbuf  ( 2) = i3e_gme
    intbuf  ( 3) = kcontrol_fi
    intbuf  ( 4) = ie_in_tot
    intbuf  ( 5) = je_in_tot
    intbuf  ( 6) = ke_in_tot
    intbuf  ( 7) = east_add_in
    intbuf  ( 8) = west_add_in
    intbuf  ( 9) = south_add_in
    intbuf  (10) = north_add_in
    intbuf  (11) = ke_soil_in
    intbuf  (12) = nlevskip
    intbuf  (13) = ke_hybrid
    realbuf ( 1) = pollat_in
    realbuf ( 2) = pollon_in
    realbuf ( 3) = polgam_in
    realbuf ( 4) = dlat_in
    realbuf ( 5) = dlon_in
    realbuf ( 9) = startlat_in_tot
    realbuf (10) = startlon_in_tot
    realbuf (11) = endlat_in_tot
    realbuf (12) = endlon_in_tot
    realbuf (13) = pcontrol_fi
    IF (lmulti_layer_in) THEN
      DO i = 1, ke_soil_in+1
        realbuf(13+i) = czmls_in(i)
      ENDDO
    ELSE
      realbuf(14:33) = 0.0_ireals
    ENDIF
    DO i = 1, 50
      realbuf(33+i) = press_level(i)
    ENDDO
    logbuf  (1: 2) = lushift_in
    logbuf  (3: 4) = lvshift_in
  ENDIF

  CALL distribute_values  (intbuf ,13, 0, imp_integers, icomm_world, ierr)
  CALL distribute_values  (realbuf,83, 0, imp_reals,    icomm_world, ierr)
  CALL distribute_values  (logbuf,  4, 0, imp_logical,  icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    ni_gme          = intbuf  ( 1)
    i3e_gme         = intbuf  ( 2)
    kcontrol_fi     = intbuf  ( 3)
    ie_in_tot       = intbuf  ( 4)
    je_in_tot       = intbuf  ( 5)
    ke_in_tot       = intbuf  ( 6)
    east_add_in     = intbuf  ( 7)
    west_add_in     = intbuf  ( 8)
    south_add_in    = intbuf  ( 9)
    north_add_in    = intbuf  (10)
    ke_soil_in      = intbuf  (11)
    nlevskip        = intbuf  (12)
    ke_hybrid       = intbuf  (13)
    pollat_in       = realbuf ( 1)
    pollon_in       = realbuf ( 2)
    polgam_in       = realbuf ( 3)
    dlat_in         = realbuf ( 4)
    dlon_in         = realbuf ( 5)
    startlat_in_tot = realbuf ( 9)
    startlon_in_tot = realbuf (10)
    endlat_in_tot   = realbuf (11)
    endlon_in_tot   = realbuf (12)
    pcontrol_fi     = realbuf (13)
    IF (lmulti_layer_in) THEN
      ALLOCATE (czmls_in(1:ke_soil_in+1), STAT=ierr)
      ALLOCATE (czhls_in(0:ke_soil_in+1), STAT=ierr)
      ALLOCATE (msoilgrib_in(0:ke_soil_in+1), STAT=ierr)
      msoilgrib_in(0) = 0_iintegers
      DO i = 1, ke_soil_in+1
        czmls_in(i)     = realbuf(13+i)
        msoilgrib_in(i) = NINT (100.0_ireals * czmls_in(i)+1.0E-7_ireals)
      ENDDO
    ENDIF
    DO i = 1, 50
      press_level(i) = realbuf(33+i)
    ENDDO
    lushift_in      = logbuf  (1: 2)
    lvshift_in      = logbuf  (3: 4)
  ENDIF

ENDIF

IF (lgme2lm) THEN
  ke_in = i3e_gme
ELSE   !  all other models
  ke_in = ke_in_tot - nlevskip
ENDIF

ALLOCATE ( longitudes_in(ie_in_tot),  latitudes_in(je_in_tot), &
          slongitudes_in(ie_in_tot), slatitudes_in(je_in_tot), STAT=ierr)

IF (.NOT. lcm2lm) THEN
  ! for cm it is calculated in read_nc_axis called in src_decomposition
  endlat_in_tot = startlat_in_tot + (je_in_tot-1) * dlat_in
  endlon_in_tot = startlon_in_tot + (ie_in_tot-1) * dlon_in

  ! Limit values of endlon_in_tot to -180.0 ... +180.0
  IF (endlon_in_tot > 180.0_ireals) THEN
    endlon_in_tot = endlon_in_tot - 360.0_ireals
  ENDIF

  DO i = 1+west_add_in, ie_in_tot-east_add_in
    longitudes_in(i) = startlon_in_tot + (i-1) * dlon_in
!   IF (lushift_in(1) .OR. lvshift_in(1)) THEN
    IF (llm2lm .OR. lum2lm .OR.  lhir2lm) THEN
      slongitudes_in(i) = longitudes_in(i) + 0.5_ireals * dlon_in
    ELSE
      slongitudes_in(i) = longitudes_in(i)
    ENDIF
  ENDDO

  DO i = 1+south_add_in, je_in_tot-north_add_in
    latitudes_in(i) = startlat_in_tot + (i-1) * dlat_in
!   IF (lushift_in(2) .OR. lvshift_in(2)) THEN
    IF (llm2lm .OR. lum2lm .OR.  lhir2lm) THEN
      slatitudes_in(i) = latitudes_in(i) + 0.5 * dlat_in
    ELSE
      slatitudes_in(i) = latitudes_in(i)
    ENDIF
  ENDDO

  IF (west_add_in /= 0) THEN
    longitudes_in (1) = longitudes_in (ie_in_tot - east_add_in) - 360.0_ireals
    slongitudes_in(1) = slongitudes_in(ie_in_tot - east_add_in) - 360.0_ireals
  ENDIF
  IF (east_add_in /= 0) THEN
    longitudes_in (ie_in_tot) = 360.0_ireals + longitudes_in (1 + west_add_in)
    slongitudes_in(ie_in_tot) = 360.0_ireals + slongitudes_in(1 + west_add_in)
  ENDIF
  IF (south_add_in /= 0) THEN
    latitudes_in (1) = -90._ireals
    slatitudes_in(1) = -90._ireals
  ENDIF
  IF (north_add_in /= 0) THEN
    latitudes_in (je_in_tot) = 90._ireals
    slatitudes_in(je_in_tot) = 90._ireals
  ENDIF
ENDIF

IF (lmulti_layer_in) THEN
  ! determine depth of half levels out of czmls
  IF (lec2lm) THEN
    ! fixed IFS soil depths
    czhls_in(:) = (/0.0_ireals, 0.07_ireals, 0.28_ireals, 1.0_ireals, 4.0_ireals/)
  ELSE
    czhls_in(0) = 0.0_ireals
    DO i = 1, ke_soil_in+1
      czhls_in(i) = 2.0_ireals * czmls_in(i) - czhls_in(i-1)
    ENDDO
  ENDIF
ENDIF

ke_pressure = ke_in  ! to save this value throughout the program
kedim_in    = MAX (ke_in, ke_hybrid)

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(A23)') '      NAMELIST: grid_in'
  WRITE (nuout, '(A23)') '      -----------------'
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value',       &
                                               'Default Value', 'Format'

  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                           'ni_gme', ni_gme ,ni_gme_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                        'i3e_gme', i3e_gme, i3e_gme_d, ' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                            'kcontrol_fi', kcontrol_fi, kcontrol_fi_d, ' I '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                            'pcontrol_fi', pcontrol_fi, pcontrol_fi_d, ' R '

  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                  'ie_in_tot', ie_in_tot, ie_in_tot_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                  'je_in_tot', je_in_tot, je_in_tot_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                  'ke_in_tot', ke_in_tot, ke_in_tot_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                  'ke_hybrid', ke_hybrid, ke_hybrid_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                     'nlevskip', nlevskip, nlevskip_d ,' I '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                  'pollat_in', pollat_in, pollat_in_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                  'pollon_in', pollon_in, pollon_in_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                  'polgam_in', polgam_in, polgam_in_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                        'dlat_in', dlat_in, dlat_in_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                        'dlon_in', dlon_in, dlon_in_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                'startlat_in_tot', startlat_in_tot, startlat_in_tot_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                'startlon_in_tot', startlon_in_tot, startlon_in_tot_d ,' R '
  IF (.NOT. lcm2lm) THEN
    WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                        &
                      'endlat_in_tot', endlat_in_tot, endlat_in_tot_d ,' R '
    WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                        &
                      'endlon_in_tot', endlon_in_tot, endlon_in_tot_d ,' R '
  ENDIF
  WRITE (nuout, '(T8,A,T21,2L6,T40,2L6,T59,A3)')                              &
                               'lushift_in', lushift_in, lushift_in_d ,' L '
  WRITE (nuout, '(T8,A,T21,2L6,T40,2L6,T59,A3)')                              &
                               'lvshift_in', lvshift_in, lvshift_in_d ,' L '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                            'east_add_in', east_add_in, east_add_in_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                            'west_add_in', west_add_in, west_add_in_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                         'south_add_in', south_add_in, south_add_in_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                         'north_add_in', north_add_in, north_add_in_d ,' I '

  IF (lmulti_layer_in) THEN
    ! Write specification for soil levels
    WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                        &
                                        'ke_soil_in ',ke_soil_in ,ke_soil_in_d, ' I '
    WRITE (nuout, '(T10,A)') 'Main levels of input soil layers'
    WRITE (nuout, '(T10,A)') '                  (m)           (cm)'
    WRITE (nuout, '(A,I12)') '              0:               ', msoilgrib_in(0)
    DO i = 1, ke_soil_in+1
      WRITE (nuout, '(I15,A,F12.4,I12)') i, ':   ',czmls_in(i), msoilgrib_in(i)
    ENDDO
  ENDIF
  WRITE (nuout, '(A2)')  '  '

  IF (lgfs2lm) THEN
  WRITE (nuout, '(T8,A)') 'Pressure Levels for GFS: '
    DO i = 1, ke_in
      WRITE (nuout, '(I15,A,F12.4)') i, ':   ',press_level(i)
    ENDDO
  ENDIF

  WRITE (nuout, '(A2)')  '  '
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_grid_in

!==============================================================================
!+ Module procedure in "setup" for the input of NAMELIST lmgrid
!------------------------------------------------------------------------------

SUBROUTINE input_lmgrid (realbuf, intbuf, logbuf, charbuf, ibuflen,          &
                         nuout, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group lmgrid. 
!   The group lmgrid contains variables defining the size and resolution of
!   the LM grid (horizontal and vertical).
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for 
!   consistency. If wrong input values are detected the program prints 
!   an error message. The program is not stopped in this routine but an 
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists. 
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.    
!   Both, default and input values are written to the file OUTPUT.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    ibuflen,      & ! dimension of the buffers
    nuout,        & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Sending buffers (are given to routine distribute_values)
  INTEGER (KIND=iintegers), INTENT(INOUT)    ::  intbuf  (ibuflen)
  REAL    (KIND=ireals)   , INTENT(INOUT)    ::  realbuf (ibuflen)
  LOGICAL                 , INTENT(INOUT)    ::  logbuf  (ibuflen)
  CHARACTER (LEN=100)     , INTENT(INOUT)    ::  charbuf (ibuflen)

! Local variables

! Variables for default values
  INTEGER (KIND=iintegers)   ::       &
    ielm_tot_d,  & ! number of gridpoints in "west-east"   direction
    jelm_tot_d,  & ! number of gridpoints in "south-north" direction
    kelm_tot_d,  & ! number of levels in the vertical
    ke_soil_lm_d,& ! number of layers in multi-layer soil model in output
    ivctype, ivctype_d, & ! vertical coordinate type (for LM)
    irefatm, irefatm_d    ! type of the reference atmosphere

  LOGICAL                    ::       &
    lanalyt_calc_t0p0_d,&! to compute values p0, t0 of reference atmosphere
                         ! analytically on full levels (only for irefatm=1)
    lnewVGrid_d          ! to write a new vertical grid HHL

  REAL (KIND=ireals)         ::       &
    ! horizontal grid
    pollat_d,    & ! geographical latitude  of the rotated north-pole
    pollon_d,    & ! geographical longitude of the rotated north-pole
    polgam_d,    & ! latitude of the rotated north pole
    dlon_d,      & ! grid mesh (in degrees) in "west-east"   direction
    dlat_d,      & ! grid mesh (in degrees) in "south-north" direction
    startlat_tot_d,  & ! latitude  (in degrees, rotated) of lower left corner
    startlon_tot_d     ! longitude (in degrees, rotated) of lower left corner

  REAL (KIND=ireals)         ::       &
    ! vertical grid for LM
    vcflat, vcflat_d,   & ! coordinate value where system changes back to z-system
    p0sl, p0sl_d,       & ! constant reference pressure on sea-level
    t0sl, t0sl_d,       & ! constant reference temperature on sea-level
    dt0lp, dt0lp_d,     & ! d (t0) / d (ln p0)
    delta_t, delta_t_d, & ! temp. difference between sea level and stratosphere (for irefatm = 2)
    h_scal, h_scal_d      ! scale height for irefatm = 2

  REAL (KIND=ireals)         ::       &
    svc1_d, svc2_d ! decay rates for large-scale and small-scale parts
                   ! of the topography for SLEVE coordinate

  INTEGER (KIND=iintegers)   ::       &
    nfltvc_d       ! number of applications of filter for
                   ! decomposition of topography in large-scale
                   ! and small-scale part for SLEVE coordinate

  REAL (KIND=ireals)         ::       &
    czvw_so_lm_d(20), & ! artificial volumetric soil water content (0-1) profile (defaults)
    czvw_so_lm  (20), & ! artificial volumetric soil water content (0-1) profile (read in)
    czml_soil_lm_d(20),&! depth of bottom level of soil layers (defaults)
    czml_soil_lm  (20)  ! depth of bottom level of soil layers (read in)

  INTEGER (KIND=iintegers)   :: i, k, n, invar, ierr, ierror, nzmaxrealbuf
  CHARACTER (LEN=120)        :: yerrmsg
  CHARACTER (LEN=15)         :: yroutine

! Define the namelist group
  NAMELIST /lmgrid/ &
    ielm_tot,    & ! number of gridpoints in "west-east"   direction
    jelm_tot,    & ! number of gridpoints in "south-north" direction
    kelm_tot,    & ! number of levels in the vertical
    ke_soil_lm,  & ! number of layers in multi-layer soil model in output
    pollat,      & ! geographical latitude  of the rotated north-pole
    pollon,      & ! geographical longitude of the rotated north-pole
    polgam,      & ! latitude of the rotated north pole
    dlon,        & ! grid mesh (in degrees) in "west-east"   direction
    dlat,        & ! grid mesh (in degrees) in "south-north" direction
    startlat_tot,& ! latitude  (in degrees, rotated) of lower left corner
    startlon_tot,& ! longitude (in degrees, rotated) of lower left corner
    ivctype,     & ! reference-pressure based hybrid system
    lnewVGrid,   & ! to write a new vertical grid HHL
    irefatm,     & ! type of the reference atmosphere
    lanalyt_calc_t0p0,&! to compute values p0, t0 of reference atmosphere
                       ! analytically on full levels (only for irefatm=1)
    svc1,        & ! decay rate for large-scale part of topography
    svc2,        & ! decay rate for small-scale part of topography
    nfltvc,      & ! number of filter applications for topo decomposition
    vcflat,      & ! coordinate value where system changes back to z-system
    vcoord_d,    & ! default vertical coordinate parameters
    p0sl,        & ! constant reference pressure on sea-level
    t0sl,        & ! constant reference temperature on sea-level
    dt0lp,       & ! d (t0) / d (ln p0)
    delta_t,     & ! temp. difference between sea level and stratosphere (for irefatm = 2)
    h_scal,      & ! scale height for irefatm = 2
    czml_soil_lm,& ! depth of bottom level of soil layers (read in)
    czvw_so_lm     ! artificial volumetric soil water content profile

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_lmgrid
!------------------------------------------------------------------------------

yroutine = 'input_lmgrid'
ierror   = 0

! set defaults for vertical coordinates and reference atmospheres
CALL set_vcoord_defaults
CALL set_refatm_defaults

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  ielm_tot_d          =  213_iintegers
  jelm_tot_d          =  213_iintegers
  kelm_tot_d          =   20_iintegers
  ivctype_d           =    2_iintegers
  lnewVGrid_d         = .FALSE.           ! only if a new HHL is constructed
  irefatm_d           =    1_iintegers
  lanalyt_calc_t0p0_d = .FALSE.   ! to be compatible with older versions
  svc1_d              =  1.0E4_ireals
  svc2_d              =  1.0E4_ireals
  nfltvc_d            =  100_iintegers
  vcflat_d            =    0.220_ireals  ! where levels change back to z-system
  pollat_d            =   32.5_ireals
  pollon_d            = -170.0_ireals
  polgam_d            =    0.0_ireals
  dlon_d              =    0.0625_ireals
  dlat_d              =    0.0625_ireals
  startlat_tot_d      =  -14.375_ireals
  startlon_tot_d      =   -6.875_ireals

  ! default values of the reference-atmosphere of the nonhydrostatic COSMO-Model
  p0sl_d              =   1.0E5_ireals
  t0sl_d              =  288.15_ireals
  dt0lp_d             =   42.0_ireals
  delta_t_d           =   75.0_ireals
  h_scal_d            = 10000.0_ireals
  vcoord_d(:)         = HUGE(0.0_ireals) ! defaults are set after reading the namelist group

  ke_soil_lm_d        = 7_iintegers
  czml_soil_lm_d(1)   =  0.005_ireals
  czml_soil_lm_d(2)   =  0.02_ireals
  czml_soil_lm_d(3)   =  0.06_ireals
  czml_soil_lm_d(4)   =  0.18_ireals
  czml_soil_lm_d(5)   =  0.54_ireals
  czml_soil_lm_d(6)   =  1.62_ireals
  czml_soil_lm_d(7)   =  4.86_ireals
  czml_soil_lm_d(8)   = 14.58_ireals
  czvw_so_lm_d(1)     =  0.75_ireals
  czvw_so_lm_d(2)     =  0.75_ireals
  czvw_so_lm_d(3)     =  0.75_ireals
  czvw_so_lm_d(4)     =  0.75_ireals
  czvw_so_lm_d(5)     =  0.75_ireals
  czvw_so_lm_d(6)     =  0.75_ireals
  czvw_so_lm_d(7)     =  0.75_ireals
  czvw_so_lm_d(8)     =  0.75_ireals

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with default
!------------------------------------------------------------------------------

  ielm_tot          = ielm_tot_d
  jelm_tot          = jelm_tot_d
  kelm_tot          = kelm_tot_d
  ivctype           = ivctype_d
  lnewVGrid         = lnewVGrid_d       ! only if a new HHL is constructed
  irefatm           = irefatm_d
  lanalyt_calc_t0p0 = lanalyt_calc_t0p0_d
  svc1              = svc1_d
  svc2              = svc2_d
  nfltvc            = nfltvc_d
  vcflat            = vcflat_d
  pollat            = pollat_d
  pollon            = pollon_d
  polgam            = polgam_d
  dlon              = dlon_d
  dlat              = dlat_d
  startlat_tot      = startlat_tot_d
  startlon_tot      = startlon_tot_d
  p0sl              = p0sl_d
  t0sl              = t0sl_d
  dt0lp             = dt0lp_d
  delta_t           = delta_t_d
  h_scal            = h_scal_d
  ke_soil_lm        = ke_soil_lm_d
  czvw_so_lm(:)     = -1.0_ireals
  czml_soil_lm(:)   = -1.0_ireals

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  READ (nuin, lmgrid, IOSTAT=ierror)
ENDIF

! distribute error status to all processors
IF (nproc > 1) THEN
  CALL distribute_values  (ierror, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF
IF (ierror /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4.1: Set defaults for vertical coordinate parameters and check values
!------------------------------------------------------------------------------

  ! initialize the ID with -1; if a default set is used, it will be set to
  ! the ID of this default set and later distributed to all other PEs
  vcoord%ivcoord_id = -1

  IF (vcoord_d(1) == HUGE(0.0_ireals)) THEN
    ! No values have been read for the vertical coordinate parameters
    ! check, if appropriate defaults are available in the vcoord_type-defaults:

    search_ivc_def: DO n = 1, imax_vcoordtype
      IF      ( ivctype   == vcoord_defaults(n)%ivctype .AND.        &
               kelm_tot+1 == vcoord_defaults(n)%nlevels)        THEN
        vcflat                 =  vcoord_defaults(n)%vcflat
        IF (ivctype == 1) THEN
          vcoord_d(1:kelm_tot+1) =  vcoord_defaults(n)%sigm_coord(1:kelm_tot+1)
        ELSE
          vcoord_d(1:kelm_tot+1) =  vcoord_defaults(n)%vert_coord(1:kelm_tot+1)
        ENDIF

        ! set the ID of vcoord to the ID of the default set
        vcoord%ivcoord_id = n

        EXIT search_ivc_def
      ENDIF
    ENDDO search_ivc_def

    IF (lnewVGrid .AND. (vcoord_d(1) == HUGE(0.0_ireals)) ) THEN
      ! no appropriate defaults could be found:
      WRITE (yerrmsg,'(2(A,I4),A)')                                          &
          ' *** ERROR: For kelm_tot= ',kelm_tot,', ivctype= ',ivctype,       &
          ' NO default vertical coordinate parameters'
      PRINT *, yerrmsg
      ierrstat = 2031
    ENDIF
  ELSE
    ! check whether the vertical coordinate parameters from Namelist are already
    ! coded as a default.
    search_ivc_def2: DO n = 1, imax_vcoordtype
      IF    ( ( ivctype   == vcoord_defaults(n)%ivctype) .AND.        &
              (kelm_tot+1 == vcoord_defaults(n)%nlevels) .AND.        &
              (ABS (vcflat - vcoord_defaults(n)%vcflat) < 0.001_ireals) ) THEN

        IF (ivctype == 1) THEN
          ! check whether vcoord_d and sigm_coord are identical
          IF (ALL (ABS(vcoord_d(1:kelm_tot+1) -                             &
                       vcoord_defaults(n)%sigm_coord(1:kelm_tot+1)) < 0.001_ireals)) THEN
            ! set the ID of vcoord to the ID of the default set
            vcoord%ivcoord_id = n
          ENDIF
          EXIT search_ivc_def2
        ELSE
          ! check whether vcoord_d and vert_coord are identical
          IF (ALL (ABS(vcoord_d(1:kelm_tot+1) -                             &
                       vcoord_defaults(n)%vert_coord(1:kelm_tot+1)) < 0.001_ireals)) THEN
            ! set the ID of vcoord to the ID of the default set
            vcoord%ivcoord_id = n
          ENDIF
          EXIT search_ivc_def2
        ENDIF

      ENDIF
    ENDDO search_ivc_def2
  ENDIF

  ! check of the vertical coordinate parameters
  IF (ivctype == 1) THEN
    IF (vcflat == 0.0 .OR. vcflat > 1.0 ) THEN
      PRINT *, ' *** ERROR:  vcflat not in allowed range ***'
      ierrstat = 2035
    ENDIF
    IF ( (vcoord_d(1) > vcflat) .OR. (vcoord_d(kelm_tot+1) /= 1.0) ) THEN
      PRINT *, ' *** ERROR:  wrong vertical coordinates for ivctype = 1 ***'
      ierrstat = 2036
    ENDIF
    DO k = 1, kelm_tot
      IF ( vcoord_d(k+1) < vcoord_d(k) ) THEN
        PRINT *, ' *** ERROR:  wrong vertical coordinates for ivctype = 1 ***'
        ierrstat = 2037
      ENDIF
    END DO
  ELSE IF (ivctype == 2) THEN
    IF (vcflat > 20000.0 .OR. vcflat < 1.0 ) THEN
      PRINT *, ' *** ERROR:  vcflat not in allowed range ***'
      ierrstat = 2038
    ENDIF
    IF ( (vcoord_d(1) < vcflat) .OR. (vcoord_d(kelm_tot+1) /= 0.0) ) THEN
      PRINT *, ' *** ERROR:  wrong vertical coordinates for ivctype = 2 ***'
      ierrstat = 2039
    ENDIF
    DO k= 1, kelm_tot
      IF ( vcoord_d(k+1) > vcoord_d(k) ) THEN
        PRINT *, ' *** ERROR:  wrong vertical coordinates for ivctype = 2 ***'
        ierrstat = 2040
      ENDIF
    END DO
  ELSE IF (ivctype == 3 .OR. ivctype == 4) THEN
    IF (vcflat < 0.0 ) THEN
      PRINT *, ' *** ERROR:  vcflat not in allowed range ***'
      ierrstat = 2041
    ENDIF

    IF ( (vcoord_d(1) < vcflat) .OR. (vcoord_d(kelm_tot+1) /= 0.0) ) THEN
      PRINT *, ' *** ERROR:  wrong vertical coordinates for ivctype = 3 or 4 ***'
      ierrstat = 2042
    ENDIF

    DO k= 1, kelm_tot
      IF ( vcoord_d(k+1) > vcoord_d(k) ) THEN
        PRINT *, ' *** ERROR:  wrong vertical coordinates for ivctype = 3 or 4 ***'
        ierrstat = 2043
      ENDIF
    END DO

    ! check for reasonable values of svc1, svc2 and nfltvc

    IF (svc1 > vcoord_d(1) .OR. svc1 < 0.0) THEN
      PRINT *, ' *** ERROR:  svc1 not in allowed range ***'
      ierrstat = 2044
    ENDIF

    IF (svc2 > vcoord_d(1) .OR. svc2 < 0.0) THEN
      PRINT *, ' *** ERROR:  svc2 not in allowed range ***'
      ierrstat = 2045
    ENDIF

    IF (nfltvc <= 0) THEN
      PRINT *, ' *** ERROR:  nfltvc must be greater than or equal to zero ***'
      ierrstat = 2046
    ENDIF
  ELSE
    PRINT *, ' *** ERROR:  ivctype /= 1, 2, 3 or 4. This option does not exist ***'
    ierrstat = 2047
  ENDIF ! ivctype

!------------------------------------------------------------------------------
!- Section 4.2: Check values for the reference atmosphere
!------------------------------------------------------------------------------

  ! test if p0sl is meaningful
  IF ((p0sl < 90000.0_ireals) .OR. (p0sl > 110000.0_ireals)) THEN
    PRINT *, ' *** ERROR:  p0sl must be in the range 90000.0 ... 110000.0 *** '
    ierrstat = 2033
  ENDIF

  ! test if t0sl is meaningful
  IF ((t0sl < 270.0_ireals) .OR. (t0sl > 300.0_ireals)) THEN
    PRINT *, ' *** ERROR:  t0sl must be in the range 270.0 ... 300.0 *** '
    ierrstat = 2033
  ENDIF

  IF     (irefatm == 1) THEN
    ! test if dt0lp for LM is non-negative
    IF (dt0lp  < 0.0) THEN
      PRINT *, ' *** ERROR:  dt0lp < 0.0 is not allowed ***'
      ierrstat = 2033
    ENDIF
  ELSEIF (irefatm == 2) THEN
    ! test for reasonable values of delta_t
    IF ( (delta_t  < 50.0_ireals) .OR. (delta_t > 100.0_ireals) ) THEN
      PRINT *, ' *** ERROR:  delta_t must be in the range 50.0 ... 100.0 *** '
      PRINT *, ' ***         but is delta_t = ', delta_t
      ierrstat = 2032
    ENDIF

    ! test for reasonable values of h_scal
    IF ( (h_scal  < 7000.0_ireals) .OR. (h_scal > 12000.0_ireals) ) THEN
      PRINT *, ' *** ERROR:  h_scal must be in the range 7000.0 ... 12000.0 *** '
      PRINT *, ' ***         but is h_scal  = ', h_scal
      ierrstat = 2032
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- Section 4.2: Check other values
!------------------------------------------------------------------------------

  ! test for the maximum number of layers for LM
  IF (kelm_tot+1 > khmax) THEN
    PRINT *, ' *** ERROR:  too many layers in LM kelm_tot+1 > khmax ***'
    ierrstat = 2034
  ENDIF

! -180.0 <= pollon, startlon_tot <= 180.0
  IF ( (-180.0_ireals > pollon) .OR. (pollon > 180.0_ireals) ) THEN
    PRINT *,' *** ERROR: WRONG VALUE OF VARIABLE pollon: ',pollon,'  ***'
    ierrstat = 2048
  ENDIF
  IF ((-180.0_ireals > startlon_tot) .OR. (startlon_tot > 180.0_ireals)) THEN
    PRINT *,' *** ERROR: WRONG VALUE OF VARIABLE startlon_tot: ',startlon_tot,'  ***'
    ierrstat = 2049
  ENDIF

! -90.0 <= pollat, startlat_tot <= 90.0
  IF ((-90.0_ireals > pollat) .OR. (pollat > 90.0_ireals)) THEN
    PRINT *,' *** ERROR: WRONG VALUE OF VARIABLE pollat: ',pollat,'  ***'
    ierrstat = 2050
  ENDIF
  IF ( (-90.0_ireals > startlat_tot) .OR. (startlat_tot > 90.0_ireals) ) THEN
    PRINT *,' *** ERROR: WRONG VALUE OF VARIABLE startlat_tot: ',startlat_tot,'  ***'
    ierrstat = 2051
  ENDIF

! dlon, dlat > epsilon
  IF (dlon < 1E-6_ireals) THEN
    PRINT *,' *** ERROR: WRONG VALUE OF VARIABLE dlon:  ',dlon,'  *** '
    ierrstat = 2052
  ENDIF
  IF (dlat < 1E-6_ireals) THEN
    PRINT *,' *** ERROR: WRONG VALUE OF VARIABLE dlat:  ',dlat,'  *** '
    ierrstat = 2053
  ENDIF

! ie_tot, je_tot >= 3, ke >= 1
  IF (ielm_tot < 3) THEN
    PRINT *,' *** ERROR: WRONG VALUE OF VARIABLE ielm_tot:  ',ielm_tot,'  *** '
    ierrstat = 2054
  ENDIF
  IF (jelm_tot < 3) THEN
    PRINT *,' *** ERROR: WRONG VALUE OF VARIABLE jelm_tot:  ',jelm_tot,'  *** '
    ierrstat = 2055
  ENDIF
  IF (kelm_tot < 1) THEN
    PRINT *,' *** ERROR: WRONG VALUE OF VARIABLE kelm_tot:  ',kelm_tot,'  *** '
    ierrstat = 2056
  ENDIF

  IF (lmulti_layer_lm) THEN
    ALLOCATE (czmls_lm(1:ke_soil_lm+1), STAT=ierr)
    ALLOCATE (czhls_lm(0:ke_soil_lm+1), STAT=ierr)
    ALLOCATE (msoilgrib_lm(0:ke_soil_lm+1), STAT=ierr)
                    ! (careful: the level k=1 will be coded with 1,
                    !           but is in the depth of 0.5 cm!)
      ! Check, how many soil levels have been specified
    invar = COUNT(czml_soil_lm(:) /= -1.0_ireals)
    IF (invar == 0) THEN
      ! no level specifications have been read
      IF (ke_soil_lm == ke_soil_lm_d) THEN
        ! use the default
        PRINT *,' *** Default specifications of LM soil main levels are used *** '
        czmls_lm(1:ke_soil_lm+1) = czml_soil_lm_d(1:ke_soil_lm+1)
      ELSEIF (ke_soil_lm /= ke_soil_lm_d) THEN
        PRINT *,' ERROR  *** no specifications of LM soil levels, *** '
        PRINT *,' ERROR  *** but using wrong default              *** ', &
                  ke_soil_lm_d, ke_soil_lm
        ierrstat = 1002
      ENDIF
    ELSE
      IF (ke_soil_lm+1 == invar) THEN
        ! use specifications read in
        PRINT *,'  *** specifications of LM soil main levels *** '
        PRINT *,'  *** from Namelist INPUT are used          *** '
        czmls_lm(1:ke_soil_lm+1) = czml_soil_lm(1:ke_soil_lm+1)
      ELSE
        PRINT *,' ERROR  *** wrong number of specifications ',           &
                'for LM soil levels  *** ', ke_soil_lm, invar
        ierrstat = 1002
      ENDIF
    ENDIF

    IF (lmulti_layer_in .AND. llm2lm) THEN 
      ! check if input and output soil levels are identical

      IF (ke_soil_lm /= ke_soil_in) THEN
        PRINT *,' ERROR  *** number of soil levels in input and output differ '
        PRINT *,'ke_soil_in = ',ke_soil_in, '   ke_soil_lm = ',ke_soil_lm
        ierrstat = 1002
      ENDIF

      IF (ke_soil_lm == ke_soil_in) THEN
        DO i = 1, ke_soil_lm
          IF (czmls_lm(i) /= czmls_in(i)) THEN
            PRINT *,' ERROR  *** soil levels in input and output differ '
            WRITE(*,'(A)') ' level  czmls_in      czmls_in'
            DO k = 1, ke_soil_lm
              WRITE (*,'(1X,I4,A3,2F12.4)') k, ' : ', czmls_in(k), czmls_lm(k)
            ENDDO
           ierrstat = 1002
            EXIT
          ENDIF
        ENDDO
      ENDIF

    ENDIF

    IF (ierrstat == 0) THEN
      ! compute grib coded values of the depth of main soil levels
      msoilgrib_lm(0) = 0_iintegers
      DO i = 1, ke_soil_lm+1
        msoilgrib_lm(i) = NINT (100.0_ireals * czmls_lm(i)+1.0E-7_ireals)
      ENDDO
    ENDIF

  ! Check constant volumetric moisture content
    IF (itype_w_so_rel == 0) THEN
      ALLOCATE (cw_so_rel_lm(1:ke_soil_lm+1), STAT=ierr)
      invar = COUNT(czvw_so_lm(:) /= -1.0_ireals)
      IF (invar == 0) THEN
        ! no level specifications have been read
        IF (ke_soil_lm == ke_soil_lm_d) THEN
        ! use the default
          PRINT *,' *** Default specifications of artificial vw_so are used  *** '
          cw_so_rel_lm(1:ke_soil_lm+1) = czvw_so_lm_d(1:ke_soil_lm+1)
        ELSEIF (ke_soil_lm /= ke_soil_lm_d) THEN
          PRINT *,' ERROR  *** in defining cw_so_rel_lm*** '
          PRINT *,' ERROR  *** no specifications of LM soil levels, *** '
          PRINT *,' ERROR  *** but using wrong default              *** ', &
                  ke_soil_lm_d, ke_soil_lm
          ierrstat = 1004
        ENDIF
      ELSE
        cw_so_rel_lm(1:ke_soil_lm+1) = czvw_so_lm(1:ke_soil_lm+1)
        DO k = 1, invar
          IF (czvw_so_lm(k) < 0._ireals .OR. czvw_so_lm(k) > 1._ireals) THEN
            PRINT *,'ERROR     *** cw_so_rel_lm is out of range (0-1) : ', k, czvw_so_lm(k)
            ierrstat = 1004
          ENDIF
        ENDDO
      ENDIF
    ENDIF

  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf  ( 1) = ielm_tot
    intbuf  ( 2) = jelm_tot
    intbuf  ( 3) = kelm_tot
    intbuf  ( 4) = ivctype
    intbuf  ( 5) = nfltvc
    intbuf  ( 6) = ke_soil_lm
    intbuf  ( 7) = irefatm

    logbuf  ( 1) = lanalyt_calc_t0p0
    logbuf  ( 2) = lnewVGrid

    realbuf ( 1) = vcflat
    realbuf ( 2) = pollat
    realbuf ( 3) = pollon
    realbuf ( 4) = polgam
    realbuf ( 5) = dlat
    realbuf ( 6) = dlon
    realbuf ( 7) = startlat_tot
    realbuf ( 8) = startlon_tot
    realbuf ( 9) = p0sl
    realbuf (10) = t0sl
    realbuf (11) = dt0lp
    realbuf (12) = svc1
    realbuf (13) = svc2
    IF (lmulti_layer_lm) THEN
      nzmaxrealbuf = 13
      DO i = 1, ke_soil_lm+1
        realbuf(nzmaxrealbuf+i) = czmls_lm(i)
      ENDDO
      nzmaxrealbuf = nzmaxrealbuf + ke_soil_lm+1
      IF (itype_w_so_rel == 0) THEN
        DO i = 1, ke_soil_lm+1
          realbuf(nzmaxrealbuf+i) = czvw_so_lm(i)
        ENDDO
        nzmaxrealbuf = nzmaxrealbuf + ke_soil_lm+1
      ENDIF
    ELSE
      realbuf(13:52) = 0.0_ireals
    ENDIF
    realbuf (53) = delta_t
    realbuf (54) = h_scal
  ENDIF

  CALL distribute_values  (intbuf ,      7, 0, imp_integers,  icomm_world, ierr)
  CALL distribute_values  (logbuf ,      2, 0, imp_logical,   icomm_world, ierr)
  CALL distribute_values  (realbuf,     54, 0, imp_reals,     icomm_world, ierr)
  CALL distribute_values  (vcoord_d, khmax, 0, imp_reals,     icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    ielm_tot    = intbuf  ( 1)
    jelm_tot    = intbuf  ( 2)
    kelm_tot    = intbuf  ( 3)
    ivctype     = intbuf  ( 4)
    nfltvc      = intbuf  ( 5)
    ke_soil_lm  = intbuf  ( 6)
    irefatm     = intbuf  ( 7)

    lanalyt_calc_t0p0 = logbuf  ( 1)
    lnewVGrid         = logbuf  ( 2)

    vcflat      = realbuf ( 1)
    pollat      = realbuf ( 2)
    pollon      = realbuf ( 3)
    polgam      = realbuf ( 4)
    dlat        = realbuf ( 5)
    dlon        = realbuf ( 6)
    startlat_tot= realbuf ( 7)
    startlon_tot= realbuf ( 8)
    p0sl        = realbuf ( 9)
    t0sl        = realbuf (10)
    dt0lp       = realbuf (11)
    svc1        = realbuf (12)
    svc2        = realbuf (13)
    nzmaxrealbuf = 13
    IF (lmulti_layer_lm) THEN
      ALLOCATE (czmls_lm(1:ke_soil_lm+1), STAT=ierr)
      ALLOCATE (czhls_lm(0:ke_soil_lm+1), STAT=ierr)  !_br
      ALLOCATE (cw_so_rel_lm(1:ke_soil_lm+1), STAT=ierr)
      ALLOCATE (msoilgrib_lm(0:ke_soil_lm+1), STAT=ierr)

      msoilgrib_lm(0) = 0_iintegers
      DO i = 1, ke_soil_lm+1
        czmls_lm(i)  = realbuf(nzmaxrealbuf+i)
        msoilgrib_lm(i) = NINT (100.0_ireals * czmls_lm(i)+1.0E-7_ireals)
      ENDDO
      nzmaxrealbuf = nzmaxrealbuf + ke_soil_lm+1
      IF (itype_w_so_rel == 0) THEN
        DO i = 1, ke_soil_lm+1
          cw_so_rel_lm(i)  = realbuf(nzmaxrealbuf+i)
        ENDDO
      ENDIF
    ENDIF
    delta_t     = realbuf (53)
    h_scal      = realbuf (54)
  ENDIF

ENDIF

! variables related to the horizontal sizes of the fields 
ie2lm_tot  = ielm_tot + 2*nboundlines
je2lm_tot  = jelm_tot + 2*nboundlines
endlon_tot = startlon_tot + (ielm_tot-1)*dlon
IF (endlon_tot > 180.0_ireals) THEN
  endlon_tot = endlon_tot - 360.0_ireals
ENDIF
endlat_tot = startlat_tot + (jelm_tot-1)*dlat

! maximal vertical dimension
kedim = MAX (kelm_tot, kedim_in)

IF (lmulti_layer_lm) THEN
  ! determine depth of half levels out of czmls
  czhls_lm(0) = 0.0_ireals
  DO i = 1, ke_soil_lm+1
    czhls_lm(i) = 2.0_ireals * czmls_lm(i) - czhls_lm(i-1)
  ENDDO
ENDIF

! Set type vcoord
vcoord%ivctype    = ivctype
vcoord%nlevels    = kelm_tot+1
vcoord%vc_uuid(:) = 'x'
vcoord%vcflat     = vcflat

IF (ivctype == 1) THEN
  vcoord%sigm_coord(1:vcoord%nlevels) = vcoord_d(1:vcoord%nlevels)
ELSE
  vcoord%vert_coord(1:vcoord%nlevels) = vcoord_d(1:vcoord%nlevels)
ENDIF
    
IF (num_compute > 1) THEN
  ! to set the ID, the value from PE 0, which has been set above, has to be 
  ! distributed to all other PEs
  CALL distribute_values  (vcoord%ivcoord_id, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF

IF (my_world_id == 0) THEN
  IF (vcoord%ivcoord_id > 0) THEN
    print *, ' *** A default set for vcoord parameters is used: ', vcoord%ivcoord_id
  ELSE
    print *, ' *** A   new   set for vcoord parameters is used: ', vcoord%ivcoord_id
  ENDIF
ENDIF

! Set type refatm
refatm%irefatm     = irefatm
refatm%irefatm_id  = 0
refatm%p0sl        = p0sl
refatm%t0sl        = t0sl
IF     (irefatm == 1) THEN
  refatm%dt0lp     = dt0lp
  refatm%delta_t   = rundefined
  refatm%h_scal    = rundefined
ELSEIF (irefatm == 2) THEN
  refatm%dt0lp     = dt0lp     ! not really necessary, but for backwards compatibility
  refatm%delta_t   = delta_t
  refatm%h_scal    = h_scal
ENDIF

refatm%irefatm_id = 0
search_refatm: DO n = 1, imax_refatmtype
  IF     (irefatm == 1) THEN
    IF    ( (irefatm    == refatm_defaults(n)%irefatm)                .AND.      &
            (ABS(p0sl    - refatm_defaults(n)%p0sl   ) < 1E-5_ireals)    .AND.   &
            (ABS(t0sl    - refatm_defaults(n)%t0sl   ) < 1E-5_ireals)    .AND.   &
            (ABS(dt0lp   - refatm_defaults(n)%dt0lp  ) < 1E-5_ireals) ) THEN
      refatm%irefatm_id = n
      EXIT search_refatm
    ENDIF
  ELSEIF (irefatm == 2) THEN
    IF    ( (irefatm    == refatm_defaults(n)%irefatm)                .AND.      &
            (ABS(p0sl    - refatm_defaults(n)%p0sl   ) < 1E-5_ireals)    .AND.   &
            (ABS(t0sl    - refatm_defaults(n)%t0sl   ) < 1E-5_ireals)    .AND.   &
            (ABS(delta_t - refatm_defaults(n)%delta_t) < 1E-5_ireals)    .AND.   &
            (ABS(h_scal  - refatm_defaults(n)%h_scal ) < 1E-5_ireals) ) THEN
      refatm%irefatm_id = n
      EXIT search_refatm
    ENDIF
  ENDIF
ENDDO search_refatm


IF (my_world_id == 0) THEN
  IF (refatm%irefatm_id > 0) THEN
    print *, ' *** A default set for refatm parameters is used: ', refatm%irefatm_id
  ELSE
    print *, ' *** A   new   set for refatm parameters is used: ', refatm%irefatm_id
  ENDIF
ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(A23)') '      NAMELIST:  lmgrid'
  WRITE (nuout, '(A23)') '      -----------------'
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value',       &
                                               'Default Value', 'Format'

  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                     'ielm_tot', ielm_tot, ielm_tot_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                     'jelm_tot', jelm_tot, jelm_tot_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                     'kelm_tot', kelm_tot, kelm_tot_d ,' I '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                           'pollat', pollat, pollat_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                           'pollon', pollon, pollon_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                           'polgam', polgam, polgam_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                                 'dlat', dlat, dlat_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                                 'dlon', dlon, dlon_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                         'startlat_tot', startlat_tot, startlat_tot_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                         'startlon_tot', startlon_tot, startlon_tot_d ,' R '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                        'ivctype', ivctype, ivctype_d ,' I '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                  'lnewVGrid', lnewVGrid, lnewVGrid_d, ' L '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                        'irefatm', irefatm, irefatm_d ,' I '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
          'lanalyt_calc_t0p0', lanalyt_calc_t0p0, lanalyt_calc_t0p0_d, ' L '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                           'vcflat', vcflat, vcflat_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                                 'p0sl', p0sl, p0sl_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                                 't0sl', t0sl, t0sl_d ,' R '
  WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                              'dt0lp', dt0lp, dt0lp_d ,' R '
  IF (irefatm == 2) THEN
     WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                       &
                                        'delta_t', delta_t, delta_t_d ,' R '
     WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                       &
                                           'h_scal', h_scal, h_scal_d ,' R '
  ENDIF

  IF (ivctype == 3 .OR. ivctype == 4) THEN
     WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                        'nfltvc', nfltvc, nfltvc_d ,' I '
     WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                              'svc1', svc1, svc1_d ,' R '
     WRITE (nuout, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                          &
                                              'svc2', svc2, svc2_d ,' R '
  ENDIF

  WRITE (nuout, '(A2)')  '  Vertical Coordinate Parameters:  '
  DO k = 1, kelm_tot+1
      WRITE(nuout   ,'(10X,I5,F12.4)' ) k, vcoord_d(k)
  ENDDO
  WRITE (nuout, '(A2)')  '  '

  IF (lmulti_layer_lm) THEN
    ! Write specification for soil levels
    IF (itype_w_so_rel == 0) THEN
      WRITE (nuout, '(T10,A)') 'artificial volumetric soil water content profile'
      WRITE (nuout, '(T10,A)') '                  (0-1)'
      DO i = 1, ke_soil_lm+1
        WRITE (nuout, '(I15,A,F12.4)') i, ':   ',cw_so_rel_lm(i)
      ENDDO
    ENDIF
    WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                        'ke_soil_lm ',ke_soil_lm ,ke_soil_lm_d, ' I '
    WRITE (nuout, '(T10,A)') 'Main levels of the LM soil layers'
    WRITE (nuout, '(T10,A)') '                  (m)           (cm)'
    WRITE (nuout, '(A,I12)') '              0:               ', msoilgrib_lm(0)
    DO i = 1, ke_soil_lm+1
      WRITE (nuout, '(I15,A,F12.4,I12)') i, ':   ',czmls_lm(i), msoilgrib_lm(i)
    ENDDO
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_lmgrid

!==============================================================================
!+ Module procedure in "src_setup_gme2lm" for the input of NAMELIST dbase
!------------------------------------------------------------------------------

SUBROUTINE input_dbase  (realbuf, intbuf, logbuf, charbuf, ibuflen,          &
                         nuout, nuin, ierrstat )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group database.
!   Unlike the other routines that read a NAMELIST-group, this routine just
!   calls a routine of MPE_IO for reading the group dbase, because these
!   values must only be known by the IO-PEs.
!
! Method:
!   Directly call module procedure mpe_readnldb of => MPE_IO
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    ibuflen,      & ! dimension of the buffers
    nuout,        & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Sending buffers (are given to routine distribute_values)
  INTEGER (KIND=iintegers), INTENT(INOUT)    ::  intbuf  (ibuflen)
  REAL    (KIND=ireals)   , INTENT(INOUT)    ::  realbuf (ibuflen)
  LOGICAL                 , INTENT(INOUT)    ::  logbuf  (ibuflen)
  CHARACTER (LEN=100)     , INTENT(INOUT)    ::  charbuf (ibuflen)

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_dbase
!------------------------------------------------------------------------------

  IF (my_world_id == 0) THEN
    CALL mpe_readnldb(nuin, nuout, ierrstat)
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_dbase

!==============================================================================
!+ Module procedure in "setup" for the input of NAMELIST data
!------------------------------------------------------------------------------

SUBROUTINE input_data (realbuf, intbuf, logbuf, charbuf, ibuflen,          &
                       nuout, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group data. 
!   The group data contains variables defining the directories and names
!   of files and/or databases where data are read from or written to.
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for 
!   consistency. If wrong input values are detected the program prints 
!   an error message. The program is not stopped in this routine but an 
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists. 
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.    
!   Both, default and input values are written to the file OUTPUT.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    ibuflen,      & ! dimension of the buffers
    nuout,        & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Sending buffers (are given to routine distribute_values)
  INTEGER (KIND=iintegers), INTENT(INOUT)    ::  intbuf  (ibuflen)
  REAL    (KIND=ireals)   , INTENT(INOUT)    ::  realbuf (ibuflen)
  LOGICAL                 , INTENT(INOUT)    ::  logbuf  (ibuflen)
  CHARACTER (LEN=100)     , INTENT(INOUT)    ::  charbuf (ibuflen)

! Local variables

! Variables for default values
  LOGICAL                          ::           &
    lbd_frame_d,   & ! if .TRUE., boundary fields include only frames
    l_ke_in_gds_d    ! if .FALSE., use former DWD specific Grib settings

  INTEGER (KIND=iintegers)   ::       &
    nlocaldefnr_d, & ! local definition number for GRIB2 local section template
    nprocess_ini_d,& ! type number of initial LM/HM-data
    nprocess_bd_d, & ! type number of boundary LM/HM-data
    ncenter_d,     & ! originating center identification
    nsubcenter_d,  & ! originating subcenter identification
    nrbit_d,       & ! packrate for grib records
    nbitmap_d,     & ! user dimension for bitmaps
    nvers_d,       & ! for documenting purposes in Grib Code
    nl_soil_in_d,  & ! number of soil layers in GME
    nl_soil_lm_d,  & ! number of soil layers in LM, resp. HM
    ie_ext_d,      & ! west-east size of fields with external parameters
    je_ext_d,      & ! north-south size of fields with external parameters
                     ! these values may be larger than the actual sizes
    npstrframe_d     ! thickness of output frames

  INTEGER (KIND=iintegers)   ::       &
    nzylen, invar, i ! length of catalog-names

  INTEGER (KIND=iintegers)   ::       &
    ncglob_realization_d  ! nr. of realization of the experiment

  CHARACTER (LEN=250)        ::       &
    ylmext_cat_d,  & ! directory of the external fields for LM/HM
    yinext_cat_d,  & ! directory of the external fields for GME
    ylm_cat_d,     & ! directory of the LM/HM-fields
    yin_cat_d        ! directory of the GME-fields

  CHARACTER (LEN= 50)        ::       &
    ylm_hhl_d,     & ! name of the file with output (GRIB2) HHL fields
    yin_hhl_d,     & ! name of the file with input (GRIB2) HHL fields
    ylmext_lfn_d,  & ! name of the file with the external fields for LM/HM
    yinext_lfn_d     ! name of the file with the external fields for GME

  CHARACTER (LEN=  3)        ::       &
    ymode_read_d,  & ! mode for opening the (read) Grib files
    ymode_write_d    ! mode for opening the (write) Grib files

  CHARACTER (LEN= 10)        ::       &
    yvarini_d (nmaxlist), & ! list of initial fields for the LM
    yvarini   (nmaxlist), & ! list of initial fields for the LM
    yvarbd_d  (nmaxlist), & ! list of boundary fields for the LM
    yvarbd    (nmaxlist)    ! list of boundary fields for the LM

  CHARACTER (LEN= 80)        ::       &
    ybitmap_cat,   & ! directory of an optional bitmap for GME data
    ybitmap_cat_d    ! for default value

  CHARACTER (LEN= 80)        ::       &
    ybitmap_lfn,   & ! name of the file with an optional bitmap for GME data
    ybitmap_lfn_d    ! for default value

  INTEGER (KIND=iintegers)   :: ierr, ierror, nzmxini_d, nzmxbd_d
  CHARACTER (LEN=120)        :: yerrmsg
  CHARACTER (LEN=15)         :: yroutine

  CHARACTER (LEN= 1)  :: ytunit_in_d, ytunit_out_d    ! time units

  CHARACTER (LEN=4)                ::           &
    ylmext_form_read_d,   & ! input format of external LM data
    yinext_form_read_d,   & ! input format of external boundary data
    yin_form_read_d,      & ! input format of boundary data
    ylm_form_write_d        ! output format of LM data

  CHARACTER (LEN=8)                ::           &
    yinput_type_d           ! type of input data: 'forecast', 'analysis' or 'ana_init'

  CHARACTER (LEN=100)  ::    &
    yncglob_institution_d,   & ! originating center name
    yncglob_title_d,         & ! title string for the output
    yncglob_source_d,        & ! program name and version
    yncglob_project_id_d,    & ! identification of the project of simulation
    yncglob_experiment_id_d, & ! identification of the experiment of simulation
    yncglob_contact_d,       & ! contact e.g. email address
    yncglob_references_d       ! URL, report etc.

! Define the namelist group
  NAMELIST /data/ &
    ylmext_cat,  & ! directory of the external fields for LM/HM
    ylmext_lfn,  & ! name of the file with the external fields for LM/HM
    yinext_cat,  & ! directory of the external fields for GME
    yinext_lfn,  & ! name of the file with the external fields for GME
    ybitmap_cat, & ! directory of an optional bitmap for GME data
    ybitmap_lfn, & ! name of the file with an optional bitmap for GME data
    yin_cat,     & ! directory of the GME-fields
    ylm_cat,     & ! directory of the LM/HM-fields
    ylm_hhl,     & ! name of the file with output (GRIB2) HHL fields
    yin_hhl,     & ! name of the file with input (GRIB2) HHL fields
    ymode_read,  & ! mode for opening the (read) Grib files
    ymode_write, & ! mode for opening the (write) Grib files
    nvers,       & ! for documenting purposes in Grib Code
    yvarini,     & ! list of initial fields for LM
    yvarbd,      & ! list of boundary fields for LM
    nlocaldefnr, & ! local definition number for GRIB2 local section template
    nprocess_ini,& ! type number of initial LM/HM-data
    nprocess_bd, & ! type number of boundary LM/HM-data
    ncenter,     & ! originating center identification
    nsubcenter,  & ! originating subcenter identification
    nrbit,       & ! packrate for grib records
    nbitmap,     & ! user dimension for bitmaps
    nl_soil_in,  & ! number of soil layers in GME
    nl_soil_lm,  & ! number of soil layers in LM, resp. HM
    ie_ext,      & ! west-east size of fields with external parameters
    je_ext,      & ! north-south size of fields with external parameters
                   ! these values may be larger than the actual sizes
    npstrframe,  & ! thickness of output frames
    lbd_frame,   & ! if .TRUE., boundary fields include only frames
    l_ke_in_gds, & ! if .FALSE., use former DWD specific Grib settings
    ytunit_in,   & ! time unit for input data
    ytunit_out,  & ! time unit for output data
    yinput_type,        & ! type of input data: 'forecast', 'analysis' or 'ana_init'
    ylmext_form_read,   & ! input format of external LM data
    yinext_form_read,   & ! input format of external boundary data
    yin_form_read,      & ! input format of boundary data
    ylm_form_write        ! output format of LM data

  NAMELIST /data/ &
    yncglob_institution,   & ! originating center name
    yncglob_title,         & ! title string for the output
    yncglob_source,        & ! program name and version
    yncglob_project_id,    & ! identification of the project of simulation
    yncglob_experiment_id, & ! identification of the experiment of simulation
    yncglob_contact,       & ! contact e.g. email address
    yncglob_references,    & ! URL, report etc.
    ncglob_realization       ! nr. of realization of the experiment

! This COMMON-Block is necessary to get the name and directory of the optional
! GME-bitmap into the supplib because the parameterlists of the subroutines 
! in these libs must not be changed. This COMMON-Block is also contained in 
! the routine bitmin of supplib.

  COMMON / gme_bitmap / ybitmap_cat, ybitmap_lfn

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_data
!------------------------------------------------------------------------------

yroutine = 'input_data'
ierror   = 0

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  ylmext_cat_d   = '/rhome/routfor/routfox/lm/const/'
  ylmext_lfn_d   = 'lm7boden.stf'
  yinext_cat_d   = '/rhome/routfor/routfox/gme/const/'
  yinext_lfn_d   = 'invar.i128a'
  ylm_hhl_d      = 'COSMO_HHL_xxx_yyy.g1'
  yin_hhl_d      = 'COSMO_HHL_xxx_yyy.g1'
  yin_cat_d      = '   '
  ylm_cat_d      = '   '
  ymode_read_d   = 'r  '
  ymode_write_d  = 'w  '
  nvers_d        =    1_iintegers
  nlocaldefnr_d  =   -1_iintegers     ! means: not defined
  nprocess_ini_d =  -999999_iintegers
  nprocess_bd_d  =  -999999_iintegers
  ncenter_d      =   78_iintegers
  nsubcenter_d   =  255_iintegers
  nrbit_d        =   16_iintegers
  ybitmap_cat_d  = '   '
  ybitmap_lfn_d  = '   '
  nbitmap_d      = 12000*2
  nl_soil_in_d   =    2_iintegers
  nl_soil_lm_d   =    2_iintegers
  ie_ext_d       = 1081_iintegers
  je_ext_d       = 1081_iintegers
  npstrframe_d   =    8_iintegers
  lbd_frame_d    = .FALSE.
  l_ke_in_gds_d  = .TRUE.

  ytunit_in_d         = 'f'
  ytunit_out_d        = 'f'
  ylmext_form_read_d  = 'grb1'
  yinext_form_read_d  = 'grb1'
  yin_form_read_d     = 'grb1'
  ylm_form_write_d    = 'grb1'
  yinput_type_d       = 'forecast'

! Variables for initial data
  yvarini_d(:)   = '          '
  yvarini_d( 1)  = 'HSURF     '; yvarini_d( 2)  = 'FR_LAND   '
  yvarini_d( 3)  = 'Z0        '; yvarini_d( 4)  = 'SOILTYP   '
  yvarini_d( 5)  = 'PLCOV     '; yvarini_d( 6)  = 'LAI       '
  yvarini_d( 7)  = 'ROOTDP    '; yvarini_d( 8)  = 'VIO3      '
  yvarini_d( 9)  = 'HMO3      '; yvarini_d(10)  = 'U         '
  yvarini_d(11)  = 'V         '; yvarini_d(12)  = 'T         '
  yvarini_d(13)  = 'QV        '; yvarini_d(14)  = 'QC        '
  yvarini_d(15)  = 'PP        '; yvarini_d(16)  = 'T_SNOW    '
  yvarini_d(17)  = 'QV_S      '; yvarini_d(18)  = 'W_I       '
  yvarini_d(19)  = 'W_SNOW    ';
  nzmxini_d      = 19

  IF (lforest) THEN
    yvarini_d(nzmxini_d+1)  = 'FOR_E     '
    yvarini_d(nzmxini_d+2)  = 'FOR_D     '
    nzmxini_d      = nzmxini_d+2
  ENDIF

  IF (lurban) THEN
    yvarini_d(nzmxini_d+1)  = 'URBAN     '
    nzmxini_d      = nzmxini_d+1
  ENDIF

  IF (lstomata) THEN
    yvarini_d(nzmxini_d+1)  = 'RSMIN     '
    nzmxini_d      = nzmxini_d+1
  ENDIF

  IF (lemiss) THEN
    yvarini_d(nzmxini_d+1)  = 'EMIS_RAD  '
    nzmxini_d      = nzmxini_d+1
  ENDIF

  IF     (itype_albedo == 2) THEN
    yvarini_d(nzmxini_d+1)  = 'ALB_DRY   '
    yvarini_d(nzmxini_d+2)  = 'ALB_SAT   '
    nzmxini_d      = nzmxini_d+2
  ELSEIF (itype_albedo == 3) THEN
    yvarini_d(nzmxini_d+1)  = 'ALB_DIF   '
    nzmxini_d      = nzmxini_d+1
  ENDIF

  IF (lsso) THEN
    yvarini_d(nzmxini_d+1)  = 'SSO_STDH  '
    yvarini_d(nzmxini_d+2)  = 'SSO_GAMMA '
    yvarini_d(nzmxini_d+3)  = 'SSO_THETA '
    yvarini_d(nzmxini_d+4)  = 'SSO_SIGMA '
    nzmxini_d      = nzmxini_d+4
  ENDIF

  IF (lradtopo) THEN
    yvarini_d(nzmxini_d+1)  = 'HORIZON   '
    yvarini_d(nzmxini_d+2)  = 'SLO_ANG   '
    yvarini_d(nzmxini_d+3)  = 'SLO_ASP   '
    yvarini_d(nzmxini_d+4)  = 'SKYVIEW   '
    nzmxini_d      = nzmxini_d+4
  ENDIF

  IF (llake) THEN
    yvarini_d(nzmxini_d+1)  = 'FR_LAKE   '
    yvarini_d(nzmxini_d+2)  = 'DEPTH_LK  '
    nzmxini_d      = nzmxini_d+2
    IF (lbdclim) THEN
      yvarini_d(nzmxini_d+3)  = 'SALT_LK   '
      nzmxini_d      = nzmxini_d+3
    ENDIF
  ENDIF

  IF (llake_coldstart) THEN
    yvarini_d(nzmxini_d+1)  = 'T_B1_LK   '
    yvarini_d(nzmxini_d+2)  = 'H_B1_LK   '
    yvarini_d(nzmxini_d+3)  = 'T_WML_LK  '
    yvarini_d(nzmxini_d+4)  = 'T_MNW_LK  '
    yvarini_d(nzmxini_d+5)  = 'T_BOT_LK  '
    yvarini_d(nzmxini_d+6)  = 'C_T_LK    '
    yvarini_d(nzmxini_d+7)  = 'H_ML_LK   '
    nzmxini_d      = nzmxini_d+7
  ENDIF

  ! If a FLake coldstart is done without using DWD external analyis,
  ! then these fields should also be added
  IF (lseaice .OR. llake_coldstart) THEN
    yvarini_d(nzmxini_d+1)  = 'T_ICE     '
    yvarini_d(nzmxini_d+2)  = 'H_ICE     '
    nzmxini_d      = nzmxini_d+2
  ENDIF

  IF (itype_aerosol == 2) THEN
    yvarini_d(nzmxini_d+1)  = 'AER_SO4   '
    yvarini_d(nzmxini_d+2)  = 'AER_DUST  '
    yvarini_d(nzmxini_d+3)  = 'AER_ORG   '
    yvarini_d(nzmxini_d+4)  = 'AER_BC    '
    yvarini_d(nzmxini_d+5)  = 'AER_SS    '
    nzmxini_d      = nzmxini_d+5
  ENDIF

  IF (lvertwind_ini) THEN
    yvarini_d(nzmxini_d+1)  = 'W         '
    nzmxini_d      = nzmxini_d+1
  ENDIF

  IF (lprog_qi) THEN
    yvarini_d(nzmxini_d+1)  = 'QI        '
    nzmxini_d      = nzmxini_d+1
  ENDIF

  IF (lprog_qr_qs) THEN
    yvarini_d(nzmxini_d+1)  = 'QR        '
    yvarini_d(nzmxini_d+2)  = 'QS        '
    nzmxini_d      = nzmxini_d+2
  ENDIF

  IF (lprog_qg) THEN
    yvarini_d(nzmxini_d+1)  = 'QG        '
    nzmxini_d      = nzmxini_d+1
  ENDIF

! the "Not lec2lm" has been removed again, because otherwise you could not
! run the old soil model with IFS data any more
! IF (.NOT. lmulti_layer_in .AND. .NOT. lec2lm) THEN
  IF (.NOT. lmulti_layer_in .AND. .NOT. lmulti_layer_lm) THEN
    ! add variables for old soil model
    yvarini_d(nzmxini_d+1)  = 'T_S       '
    yvarini_d(nzmxini_d+2)  = 'T_M       '
    yvarini_d(nzmxini_d+3)  = 'T_CL      '
    yvarini_d(nzmxini_d+4)  = 'W_G1      '
    yvarini_d(nzmxini_d+5)  = 'W_G2      '
    yvarini_d(nzmxini_d+6)  = 'W_CL      '
    nzmxini_d = nzmxini_d+6
  ENDIF
  IF (lmulti_layer_lm) THEN
    ! add variables for new soil model
    yvarini_d(nzmxini_d+1)  = 'T_SO      ';
    yvarini_d(nzmxini_d+2)  = 'W_SO      '
    yvarini_d(nzmxini_d+3)  = 'FRESHSNW  '
    nzmxini_d      = nzmxini_d+3
    IF (lprog_rho_snow) THEN
      yvarini_d(nzmxini_d+1)  = 'RHO_SNOW  '
      nzmxini_d      = nzmxini_d+1
    ENDIF
!_br 02.08.12 the following is moved down after the namelist has been read
!             because it depends also on the namelist parameter ymode_write
!    IF (lbdclim .AND. lmulti_layer_in) THEN
!      yvarini_d(nzmxini_d+1)  = 'T_S       '
!      nzmxini_d      = nzmxini_d+1
!    ENDIF
! SP, 201405
    IF (lcm2lm .AND. lmulti_layer_in) THEN
      yvarini_d(nzmxini_d+1)  = 'T_S       '
      nzmxini_d      = nzmxini_d+1
    ENDIF
  ENDIF

! iso code
  IF (liso) THEN
    yvarini_d(nzmxini_d+1)  = 'QV18O   '
    yvarini_d(nzmxini_d+2)  = 'QV2H    '
    yvarini_d(nzmxini_d+3)  = 'R18OSOIL  '
    yvarini_d(nzmxini_d+4)  = 'R2HSOIL   '
    nzmxini_d      = nzmxini_d+4
  ENDIF
! end iso code

  IF ( nzmxini_d > nmaxlist) THEN
    ierrstat = -3
    RETURN
  ENDIF

! Variables for boundary data
  yvarbd_d (:)   = '          '
  yvarbd_d ( 1)  = 'U         '; yvarbd_d ( 2)  = 'V         '
  yvarbd_d ( 3)  = 'T         '; yvarbd_d ( 4)  = 'QV        '
  yvarbd_d ( 5)  = 'QC        '; yvarbd_d ( 6)  = 'PP        '
  yvarbd_d ( 7)  = 'T_SNOW    '; yvarbd_d ( 8)  = 'QV_S      '
  yvarbd_d ( 9)  = 'W_SNOW    '
  nzmxbd_d       = 9

! the "Not lec2lm" has been removed again, because otherwise you could not
! run the old soil model with IFS data any more
! IF (.NOT. lmulti_layer_in .AND. .NOT. lec2lm) THEN
  IF (.NOT. lmulti_layer_in .AND. .NOT. lmulti_layer_lm) THEN
    ! add variables for old soil model
    yvarbd_d (10)  = 'T_S       '
    yvarbd_d (11)  = 'T_M       '
    yvarbd_d (12)  = 'W_G1      '
    yvarbd_d (13)  = 'W_G2      '
    nzmxbd_d       = 13
  ENDIF

  IF (lbdclim) THEN
    yvarbd_d(nzmxbd_d+1)  = 'PLCOV     '
    yvarbd_d(nzmxbd_d+2)  = 'LAI       '
    yvarbd_d(nzmxbd_d+3)  = 'ROOTDP    '
    yvarbd_d(nzmxbd_d+4)  = 'VIO3      '
    yvarbd_d(nzmxbd_d+5)  = 'HMO3      '
    nzmxbd_d      = nzmxbd_d+5
    IF (.NOT. lmulti_layer_in) THEN
      yvarbd_d (nzmxbd_d+1)  = 'T_CL      '
      yvarbd_d (nzmxbd_d+2)  = 'W_CL      '
      nzmxbd_d      = nzmxbd_d+2
!_br 02.08.12 the following is moved down after the namelist has been read
!             because it depends also on the namelist parameter ymode_write
!
!    ELSE
!      yvarbd_d (nzmxbd_d+1)  = 'T_S       '
!      nzmxbd_d      = nzmxbd_d+1
    ENDIF
    IF (itype_aerosol == 2) THEN
      yvarbd_d(nzmxbd_d+1)  = 'AER_SO4   '
      yvarbd_d(nzmxbd_d+2)  = 'AER_DUST  '
      yvarbd_d(nzmxbd_d+3)  = 'AER_ORG   '
      yvarbd_d(nzmxbd_d+4)  = 'AER_BC    '
      yvarbd_d(nzmxbd_d+5)  = 'AER_SS    '
      nzmxbd_d      = nzmxbd_d+5
    ENDIF
  ENDIF

  IF (lvertwind_bd) THEN
    yvarbd_d(nzmxbd_d+1)  = 'W         '
    nzmxbd_d      = nzmxbd_d+1
  ENDIF

  IF (lprog_qi) THEN
    yvarbd_d(nzmxbd_d+1)  = 'QI        '
    nzmxbd_d      = nzmxbd_d+1
  ENDIF

  IF (lprog_qr_qs) THEN
    yvarbd_d(nzmxbd_d+1)  = 'QR        '
    yvarbd_d(nzmxbd_d+2)  = 'QS        '
    nzmxbd_d      = nzmxbd_d+2
  ENDIF

  IF (lprog_qg) THEN
    yvarbd_d(nzmxbd_d+1)  = 'QG        '
    nzmxbd_d      = nzmxbd_d+1
  ENDIF

! iso code
  IF (liso) THEN
    yvarbd_d(nzmxbd_d+1)  = 'QV18O   '
    yvarbd_d(nzmxbd_d+2)  = 'QV2H    '
    yvarbd_d(nzmxbd_d+3)  = 'R18OSOIL  '
    yvarbd_d(nzmxbd_d+4)  = 'R2HSOIL   '
    nzmxbd_d      = nzmxbd_d+4
  ENDIF
! end iso code

  IF ( nzmxbd_d > nmaxlist) THEN
    ierrstat = -2
    RETURN
  ENDIF

  yncglob_institution_d   = '-'
  yncglob_title_d         = '-'
  yncglob_source_d        = '-'
  yncglob_project_id_d    = '-'
  yncglob_experiment_id_d = '-'
  yncglob_contact_d       = '-'
  yncglob_references_d    = '-'
  ncglob_realization_d    = -999

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with default
!------------------------------------------------------------------------------

  yvarini(:)    = ''
  yvarbd (:)    = ''
  ylmext_cat    = ylmext_cat_d
  ylmext_lfn    = ylmext_lfn_d
  yinext_cat    = yinext_cat_d
  yinext_lfn    = yinext_lfn_d
  ybitmap_cat   = ybitmap_cat_d
  ybitmap_lfn   = ybitmap_lfn_d
  ylm_hhl       = ylm_hhl_d
  yin_hhl       = yin_hhl_d
  yin_cat       = yin_cat_d
  ylm_cat       = ylm_cat_d
  ymode_read    = ymode_read_d
  ymode_write   = ymode_write_d
  nvers         = nvers_d
  nlocaldefnr   = nlocaldefnr_d
  nprocess_ini  = nprocess_ini_d
  nprocess_bd   = nprocess_bd_d
  ncenter       = ncenter_d
  nsubcenter    = nsubcenter_d
  nrbit         = nrbit_d
  nbitmap       = nbitmap_d
  nl_soil_in    = nl_soil_in_d
  nl_soil_lm    = nl_soil_lm_d
  ie_ext        = ie_ext_d
  je_ext        = je_ext_d
  npstrframe    = npstrframe_d
  lbd_frame     = lbd_frame_d
  l_ke_in_gds   = l_ke_in_gds_d
  ytunit_in     = ytunit_in_d
  ytunit_out    = ytunit_out_d
  ylmext_form_read  = ylmext_form_read_d
  yinext_form_read  = yinext_form_read_d
  yin_form_read     = yin_form_read_d
  ylm_form_write    = ylm_form_write_d
  yinput_type       = yinput_type_d

  yncglob_institution     = yncglob_institution_d
  yncglob_title           = yncglob_title_d
  yncglob_source          = yncglob_source_d
  yncglob_project_id      = yncglob_project_id_d
  yncglob_experiment_id   = yncglob_experiment_id_d
  yncglob_contact         = yncglob_contact_d
  yncglob_references      = yncglob_references_d
  ncglob_realization      = ncglob_realization_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  READ (nuin, data, IOSTAT=ierror)

  IF ((lbdclim .OR. ylm_form_write=='ncdf') .AND. lmulti_layer_lm) THEN
    IF (.NOT.lec2lm) THEN 
      ! DR, KIT: Otherwise the combination of ec and netcdf leads to 
      !          two variables T_S causing an error
      yvarini_d(nzmxini_d+1)  = 'T_S       '
      nzmxini_d      = nzmxini_d+1
      yvarbd_d(nzmxbd_d+1)  = 'T_S       '
      nzmxbd_d      = nzmxbd_d+1
    ENDIF
  ENDIF

ENDIF

! distribute error status to all processors
IF (nproc > 1) THEN
  CALL distribute_values  (ierror, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF
IF (ierror /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

  ! determination of the initial input variables
  invar = COUNT(yvarini(:) /= '')
  IF( invar == 0) THEN
    ! nothing has been read and the default list is used
    youtlist_ini(1:nzmxini_d)  = yvarini_d(1:nzmxini_d)
    youtlist_ini(nzmxini_d+1:) = ''
    numlist_ini         = nzmxini_d
  ELSEIF ( (invar >= 1) .AND. (invar <= nmaxlist) ) THEN
    youtlist_ini(1:invar)  = yvarini(1:invar)
    youtlist_ini(invar+1:) = ''
    numlist_ini     = invar
  ELSEIF (invar > nmaxlist) THEN
    PRINT *,' ERROR  *** number of variables for initial input is too big ***'
    ierrstat = 1002
  ELSE
    ! something went wrong: cannot determine number of variables for input
    PRINT *,' ERROR    *** cannot determine number of variables for input *** '
    ierrstat = 1002
  ENDIF

  ! determination of the boundary input variables
  invar = COUNT(yvarbd(:) /= '')
  IF( invar == 0) THEN
    ! nothing has been read and the default list is used
    youtlist_bd(1:nzmxini_d)  = yvarbd_d(1:nzmxini_d)
    youtlist_bd(nzmxini_d+1:) = ''
    numlist_bd          = nzmxbd_d
  ELSEIF ( (invar >= 1) .AND. (invar <= nmaxlist) ) THEN
    youtlist_bd(1:invar)  = yvarbd(1:invar)
    youtlist_bd(invar+1:) = ''
    numlist_bd      = invar
  ELSEIF (invar > nmaxlist) THEN
    PRINT *,' ERROR  *** number of variables for boundary input is too big ***'
    ierrstat = 1002
  ELSE
    ! something went wrong: cannot determine number of variables for input
    PRINT *,' ERROR  *** cannot determine number of variables for input ***'
    ierrstat = 1002
  ENDIF

  ! setting of the hhl output list
  ! there just is a default list with one entry
  numlist_hhl     = 1
  youtlist_hhl(1) = 'HHL       '

  ! Check length of the directory-names
  nzylen = LEN_TRIM(ylmext_cat)
  IF( nzylen > 0 ) THEN
    IF( ylmext_cat(nzylen:nzylen) /= '/') THEN
      IF( nzylen < LEN(ylmext_cat) ) THEN
        ylmext_cat = ylmext_cat (1:nzylen)//'/'
      ELSE
        PRINT *,' ERROR    *** ylmext_cat is too long *** '
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF

  nzylen = LEN_TRIM(yinext_cat)
  IF( nzylen > 0 ) THEN
    IF( yinext_cat(nzylen:nzylen) /= '/') THEN
      IF( nzylen < LEN(yinext_cat) ) THEN
        yinext_cat = yinext_cat (1:nzylen)//'/'
      ELSE
        PRINT *,' ERROR    *** yinext_cat is too long *** '
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF

  nzylen = LEN_TRIM(ylm_cat)
  IF( nzylen > 0 ) THEN
    IF( ylm_cat(nzylen:nzylen) /= '/') THEN
      IF( nzylen < LEN(ylm_cat) ) THEN
        ylm_cat = ylm_cat (1:nzylen)//'/'
      ELSE
        PRINT *,' ERROR    *** ylm_cat is too long *** '
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF

  nzylen = LEN_TRIM(yin_cat)
  IF( nzylen > 0 ) THEN
    IF( yin_cat(nzylen:nzylen) /= '/') THEN
      IF( nzylen < LEN(yin_cat) ) THEN
        yin_cat = yin_cat (1:nzylen)//'/'
      ELSE
        PRINT *,' ERROR    *** yin_cat is too long *** '
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF

  nzylen = LEN_TRIM(ybitmap_cat)
  IF ((nzylen > 0) .AND. (ybitmap_cat /= '   ')) THEN
    IF( ybitmap_cat(nzylen:nzylen) /= '/') THEN
      IF( nzylen < LEN(ybitmap_cat) ) THEN
        ybitmap_cat = ybitmap_cat (1:nzylen)//'/'
      ELSE
        PRINT *,' ERROR    *** ybitmap_cat is too long *** '
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF

  ! Check the range of the packrate
  IF (nrbit <= 0 .OR. nrbit > 24) THEN
    PRINT *,' ERROR    *** packrate nrbit is out of range: 0 < nrbit <= 24 ***'
    ierrstat = 1002
  ENDIF

! Check time unit for boundary data
  IF ((ytunit_in /= 'f') .AND. (ytunit_in /= 'd')) THEN
    PRINT *,' ERROR    *** ytunit_in not valid ', ytunit_in
    ierrstat = 1002
  ENDIF
  IF ((ytunit_out /= 'f') .AND. (ytunit_out /= 'd')) THEN
    PRINT *,' ERROR    *** ytunit_out not valid ', ytunit_out
    ierrstat = 1002
  ENDIF

 ! Check I/O format
  IF ((TRIM(ylmext_form_read) /= 'grb1') .AND. (TRIM(ylmext_form_read) /= 'ncdf') &
                                         .AND. (TRIM(ylmext_form_read) /= 'apix')) THEN
    PRINT *,' ERROR    *** ylmext_form_read not valid ', TRIM(ylmext_form_read)
    ierrstat = 1002
#ifndef GRIBDWD
  ELSEIF (TRIM(ylmext_form_read) == 'grb1') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of DWD Grib1 library *** '
    ierrstat = 1002
#endif
#ifndef GRIBAPI
  ELSEIF (TRIM(ylmext_form_read) == 'apix') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of grib_api *** '
    ierrstat = 1002
#endif
#ifndef NETCDF
  ELSEIF (TRIM(ylmext_form_read) == 'ncdf') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of NETCDF   *** '
    ierrstat = 1002
#endif
  ENDIF
  IF ((TRIM(yinext_form_read) /= 'grb1') .AND. (TRIM(yinext_form_read) /= 'ncdf') &
                                         .AND. (TRIM(yinext_form_read) /= 'apix')) THEN
    PRINT *,' ERROR    *** yinext_form_read not valid ', TRIM(yinext_form_read)
    ierrstat = 1002
#ifndef GRIBDWD
  ELSEIF (TRIM(yinext_form_read) == 'grb1') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of DWD Grib1 library *** '
    ierrstat = 1002
#endif
#ifndef GRIBAPI
  ELSEIF (TRIM(yinext_form_read) == 'apix') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of grib_api *** '
    ierrstat = 1002
#endif
#ifndef NETCDF
  ELSEIF (TRIM(yinext_form_read) == 'ncdf') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of NETCDF   *** '
    ierrstat = 1002
#endif
  ENDIF
  IF ((TRIM(yin_form_read) /= 'grb1') .AND. (TRIM(yin_form_read) /= 'ncdf')       &
                                      .AND. (TRIM(yin_form_read) /= 'apix')) THEN
    PRINT *,' ERROR    *** yin_form_read not valid ', TRIM(yin_form_read)
    ierrstat = 1002
#ifndef GRIBDWD
  ELSEIF (TRIM(yin_form_read) == 'grb1') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of DWD Grib1 library *** '
    ierrstat = 1002
#endif
#ifndef GRIBAPI
  ELSEIF (TRIM(yin_form_read) == 'apix') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of grib_api *** '
    ierrstat = 1002
#endif
#ifndef NETCDF
  ELSEIF (TRIM(yin_form_read) == 'ncdf') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of NETCDF   *** '
    ierrstat = 1002
#endif
  ENDIF
  IF ((TRIM(ylm_form_write) /= 'grb1') .AND. (TRIM(ylm_form_write) /= 'ncdf') .AND. &
      (TRIM(ylm_form_write) /= 'api1') .AND. (TRIM(ylm_form_write) /= 'api2') ) THEN
    PRINT *,' ERROR    *** ylm_form_write not valid ', TRIM(ylm_form_write)
    ierrstat = 1002
#ifndef GRIBDWD
  ELSEIF (TRIM(ylm_form_write) == 'grb1') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of DWD Grib1 library *** '
    ierrstat = 1002
#endif
#ifndef NETCDF
  ELSEIF (TRIM(ylm_form_write) == 'ncdf') THEN
    PRINT *,' ERROR    *** INT2LM not compiled for use of NETCDF   *** '
    ierrstat = 1002
#endif
  ENDIF

  IF (TRIM(ylm_form_write) == 'api2') THEN
    IF ((vcoord%ivcoord_id == -1) .AND. (.NOT. lnewVGrid)) THEN
      ! given values do not correspond to a default: this is not possible at the moment
      WRITE (yerrmsg,'(A)')                                          &
          ' *** ERROR: HHL-file will be read, but values for vcoord do not correspond to a default list '
      PRINT *, yerrmsg
      ierrstat = 2031
    ENDIF
  ENDIF

  IF (.NOT. (llm2lm .OR. lcm2lm)) THEN
    IF (TRIM(yinext_form_read) == 'ncdf') THEN
      PRINT *,' ERROR    *** yinext_form_read '
      PRINT *,'       netCDF format not implemented for the chosen input model'
      ierrstat = 1002
    ENDIF
    IF (TRIM(yin_form_read) == 'ncdf') THEN
      PRINT *,' ERROR    *** yin_form_read '
      PRINT *,'       netCDF format not implemented for the chosen input model'
     ierrstat = 1002
    ENDIF
  ENDIF

  IF (TRIM(yinext_form_read) == 'grb1' .AND. lcm2lm ) THEN
    ! GRIB input of external data does not work for this configuration
    PRINT *,' ERROR    *** yinext_form_read '
    PRINT *," GRIB input of coarse-grid external data not possible (so far) for yinput_model='CM'"
    ierrstat = 1002
  ENDIF

  IF (TRIM(yinext_form_read) == 'ncdf'   .AND. lgme2lm) THEN
    ! NetCDF input of external data does not work for this configuration
    PRINT *,' ERROR    *** yinext_form_read '
    PRINT *," NetCDF input of coarse-grid external data not possible (so far) for yinput_model='GME'"
    ierrstat = 1002
  ENDIF

  IF (TRIM(yin_form_read) == 'grb1'    .AND. lcm2lm ) THEN
    ! GRIB input of coarse grid data does not work for this configuration
    PRINT *,' ERROR    *** yin_form_read '
    PRINT *," GRIB input of coarse-grid data not possible (so far) for yinput_model='CM'"
    ierrstat = 1002
  ENDIF

  IF (TRIM(yin_form_read) == 'ncdf'      .AND. lgme2lm) THEN
    ! NeCDF input of coarse grid data does not work for this configuration
    PRINT *,' ERROR    *** yin_form_read '
    PRINT *," NetCDF input of coarse-grid data not possible (so far) for yinput_model='GME'"
    ierrstat = 1002
  ENDIF

  IF (TRIM(yin_form_read) == 'ncdf'      .AND. lec2lm ) THEN
    ! NeCDF input of ECMWF data does not work for this configuration
    PRINT *,' ERROR    *** yin_form_read '
    PRINT *," NetCDF input of coarse-grid data not possible (so far) for yinput_model='IFS'"
    ierrstat = 1002
  ENDIF

  IF ((yinput_type == 'ana_init') .AND. (ytunit_in /= 'd') ) THEN
    PRINT *,' ERROR    *** for initialized analysis, the time unit d is needed'
    ierrstat = 1002
  ENDIF

  IF (itype_ndvi == 2 .AND. TRIM(ylmext_form_read) /= 'ncdf' ) THEN
    ! GRIB input of external data does not work for this configuration
    PRINT *,' ERROR    *** ylmext_form_read '
    PRINT *,' GRIB input of external data not possible (so far) for itype_ndvi=2'
    ierrstat = 1002
  ENDIF

  ! Check consistency between lroutine and nvers for GRIB2:

  IF ( (ncenter == 78) .AND. (ylm_form_write == 'api2') ) THEN
    IF (lroutine) THEN
      ! for lroutine:  nvers=1..50 (Operations) or nvers=51..99 (Operations test)
      IF ( (nvers < 0) .OR. (nvers > 100) ) THEN
        PRINT *,' ERROR    *** wrong settings of lroutine=.TRUE. and nvers = ', nvers, ' *** '
        PRINT *,'          *** nvers must be in the range 1..99 for correct grib meta data settings'
        ierrstat = 1002
      ENDIF
    ELSE
      ! for not lroutine: nvers does not matter
      IF ( (nvers < 0) .OR. (nvers > 16384) ) THEN
        PRINT *,' ERROR    *** wrong settings of lroutine=.FALSE. and nvers = ', nvers, ' *** '
        PRINT *,'          *** nvers must be in the range 1..16384 for correct grib meta data settings'
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF
ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf  ( 1) = nprocess_ini
    intbuf  ( 2) = nprocess_bd
    intbuf  ( 3) = ncenter
    intbuf  ( 4) = nl_soil_in 
    intbuf  ( 5) = nl_soil_lm
    intbuf  ( 6) = ie_ext
    intbuf  ( 7) = je_ext
    intbuf  ( 8) = nrbit
    intbuf  ( 9) = npstrframe
    intbuf  (10) = nvers
    intbuf  (11) = numlist_ini
    intbuf  (12) = numlist_bd
    intbuf  (13) = nbitmap
    intbuf  (14) = ncglob_realization  !_br 05.06.12
    intbuf  (15) = nlocaldefnr
    intbuf  (16) = nsubcenter
    intbuf  (17) = numlist_hhl

    logbuf  ( 1) = lbd_frame
    logbuf  ( 2) = l_ke_in_gds
    charbuf ( 1)(  1:100) = ylmext_cat(  1:100)
    charbuf ( 2)(  1:100) = ylmext_cat(101:200)
    charbuf ( 3)(  1: 50) = ylmext_cat(201:250)
    charbuf ( 4) = ylmext_lfn
    charbuf ( 5)(  1:100) = yinext_cat(  1:100)
    charbuf ( 6)(  1:100) = yinext_cat(101:200)
    charbuf ( 7)(  1: 50) = yinext_cat(201:250)
    charbuf ( 8) = yinext_lfn
    charbuf ( 9) = ybitmap_cat
    charbuf (10) = ybitmap_lfn
    charbuf (11)(  1:100) = yin_cat(  1:100)
    charbuf (12)(  1:100) = yin_cat(101:200)
    charbuf (13)(  1: 50) = yin_cat(201:250)
    charbuf (14)(  1:100) = ylm_cat(  1:100)
    charbuf (15)(  1:100) = ylm_cat(101:200)
    charbuf (16)(  1: 50) = ylm_cat(201:250)
    charbuf (17) = ymode_read
    charbuf (18) = ymode_write
    charbuf (19) = ylmext_form_read
    charbuf (20) = yinext_form_read
    charbuf (21) = yin_form_read
    charbuf (22) = ylm_form_write
    charbuf (23) = ytunit_in
    charbuf (24) = ytunit_out
    charbuf (25) = yinput_type
    charbuf (26) = yncglob_institution
    charbuf (27) = yncglob_title
    charbuf (28) = yncglob_source
    charbuf (29) = yncglob_project_id
    charbuf (30) = yncglob_experiment_id
    charbuf (31) = yncglob_contact
    charbuf (32) = yncglob_references
    charbuf (33) = ylm_hhl
    charbuf (34) = yin_hhl
    DO i=1,nmaxlist
      charbuf(34+           i) = youtlist_ini(i)
    ENDDO
    DO i=1,nmaxlist
      charbuf(34+  nmaxlist+i) = youtlist_bd (i)
    ENDDO
    DO i=1,nmaxlist
      charbuf(34+2*nmaxlist+i) = youtlist_hhl(i)
    ENDDO
  ENDIF

  CALL distribute_values  (intbuf ,17, 0, imp_integers,  icomm_world, ierr)
  CALL distribute_values  (charbuf,34+3*nmaxlist, 0,                        &
                                          imp_character, icomm_world, ierr)
  CALL distribute_values  (logbuf , 2, 0, imp_logical,   icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    nprocess_ini = intbuf  ( 1)
    nprocess_bd  = intbuf  ( 2)
    ncenter      = intbuf  ( 3)
    nl_soil_in   = intbuf  ( 4)
    nl_soil_lm   = intbuf  ( 5)
    ie_ext       = intbuf  ( 6)
    je_ext       = intbuf  ( 7)
    nrbit        = intbuf  ( 8)
    npstrframe   = intbuf  ( 9)
    nvers        = intbuf  (10)
    numlist_ini  = intbuf  (11)
    numlist_bd   = intbuf  (12)
    nbitmap      = intbuf  (13)
    ncglob_realization = intbuf (14) !_br 05.06.12
    nlocaldefnr  = intbuf  (15)
    nsubcenter   = intbuf  (16)
    numlist_hhl  = intbuf  (17)
    lbd_frame    = logbuf  ( 1)
    l_ke_in_gds  = logbuf  ( 2)
    ylmext_cat(  1:100) = charbuf ( 1)(  1:100)
    ylmext_cat(101:200) = charbuf ( 2)(  1:100)
    ylmext_cat(201:250) = charbuf ( 3)(  1:100)
    ylmext_lfn          = charbuf ( 4)
    yinext_cat(  1:100) = charbuf ( 5)(  1:100)
    yinext_cat(101:200) = charbuf ( 6)(  1:100)
    yinext_cat(201:250) = charbuf ( 7)(  1:100)
    yinext_lfn          = charbuf ( 8)
    ybitmap_cat         = charbuf ( 9)
    ybitmap_lfn         = charbuf (10)
    yin_cat   (  1:100) = charbuf (11)(  1:100)
    yin_cat   (101:200) = charbuf (12)(  1:100)
    yin_cat   (201:250) = charbuf (13)(  1:100)
    ylm_cat   (  1:100) = charbuf (14)(  1:100)
    ylm_cat   (101:200) = charbuf (15)(  1:100)
    ylm_cat   (201:250) = charbuf (16)(  1:100)
    ymode_read          = charbuf (17)
    ymode_write         = charbuf (18)
    ylmext_form_read    = charbuf (19)
    yinext_form_read    = charbuf (20)
    yin_form_read       = charbuf (21)
    ylm_form_write      = charbuf (22)
    ytunit_in           = charbuf (23)
    ytunit_out          = charbuf (24)
    yinput_type         = charbuf (25)
    yncglob_institution = charbuf (26)
    yncglob_title       = charbuf (27)
    yncglob_source      = charbuf (28)
    yncglob_project_id  = charbuf (29)
    yncglob_experiment_id = charbuf (30)
    yncglob_contact     = charbuf (31)
    yncglob_references  = charbuf (32)
    ylm_hhl             = charbuf (33)
    yin_hhl             = charbuf (34)
    DO i=1,nmaxlist
      youtlist_ini(i) = charbuf(34+           i)
    ENDDO
    DO i=1,nmaxlist
      youtlist_bd (i) = charbuf(34+  nmaxlist+i)
    ENDDO
    DO i=1,nmaxlist
      youtlist_hhl(i) = charbuf(34+2*nmaxlist+i)
    ENDDO
  ENDIF

ENDIF

! Additional settings for all processors
IF (ybitmap_lfn == '   ') THEN
  ! no bitmap has been specified, so no bitmap has to be checked
  lcheck_bmm = .FALSE.
  lcheck_bmu = .FALSE.
  lcheck_bmv = .FALSE.
ELSE
  ! the bitmap has to be checked
  lcheck_bmm = .TRUE.
  lcheck_bmu = .TRUE.
  lcheck_bmv = .TRUE.
ENDIF

! set the logicals, to indicate whether HHL-files have to be read
! (depend on lnewVGrid and ylm_form_write
IF (ylm_form_write == 'api2') THEN
  IF (lnewVGrid) THEN
    lhhl_lm_read = .FALSE.
  ELSE
    lhhl_lm_read = .TRUE.
  ENDIF
ELSE
  lhhl_lm_read = .FALSE.
ENDIF

! lhhl_in_read is set in src_read_coarse_grid, if COSMO GRIB2 data are read
lhhl_in_read = .FALSE.

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(A23)') '      NAMELIST:    data'
  WRITE (nuout, '(A23)') '      -----------------'
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value',       &
                                               'Default Value', 'Format'

  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                            'nlocaldefnr', nlocaldefnr, nlocaldefnr_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                         'nprocess_ini', nprocess_ini, nprocess_ini_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                            'nprocess_bd', nprocess_bd, nprocess_bd_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                        'ncenter', ncenter, ncenter_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                               'nsubcenter', nsubcenter, nsubcenter_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                              'nrbit', nrbit, nrbit_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                        'nbitmap', nbitmap, nbitmap_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                               'nvers     ', nvers     , nvers_d      ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                               'nl_soil_lm', nl_soil_lm, nl_soil_lm_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                              'nl_soil_in', nl_soil_in, nl_soil_in_d  ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                           'ie_ext', ie_ext, ie_ext_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                           'je_ext', je_ext, je_ext_d ,' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                               'npstrframe', npstrframe, npstrframe_d ,' I '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                  'lbd_frame', lbd_frame, lbd_frame_d ,' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                       'l_ke_in_gds', l_ke_in_gds, l_ke_in_gds_d,' L '
  WRITE (nuout, '(T8,A,T21,  A  ,T40,A    ,T59,A)')                           &
                              'ymode_read ',ymode_read, ymode_read_d  ,'C*3'
  WRITE (nuout, '(T8,A,T21,  A  ,T40,A    ,T59,A)')                           &
                              'ymode_write',ymode_write, ymode_write_d,'C*3'
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                              'ylmext_cat',     TRIM(ylmext_cat),      ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                              'ylmext_lfn',     TRIM(ylmext_lfn),      ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                              'yinext_cat',     TRIM(yinext_cat),      ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                              'yinext_lfn',     TRIM(yinext_lfn),      ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                              'yin_hhl',        TRIM(yin_hhl),         ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                              'ybitmap_cat',    TRIM(ybitmap_cat),     ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                              'ybitmap_lfn',    TRIM(ybitmap_lfn),     ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                               'yin_cat',       TRIM(yin_cat),         ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                              'ylm_cat',        TRIM(ylm_cat),         ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                              'ylm_hhl',        TRIM(ylm_hhl),         ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                            'yinput_type',      TRIM(yinput_type),     ' C '
  WRITE (nuout, '(T8,A,T31,  A  ,T50,A    ,T59,A)')                           &
          'ylmext_form_read',TRIM(ylmext_form_read),ylmext_form_read_d,' C '
  WRITE (nuout, '(T8,A,T31,  A  ,T50,A    ,T59,A)')                           &
          'yinext_form_read',TRIM(yinext_form_read),yinext_form_read_d,' C '
  WRITE (nuout, '(T8,A,T31,  A  ,T50,A    ,T59,A)')                           &
                'yin_form_read',TRIM(yin_form_read), yin_form_read_d,  ' C '
  WRITE (nuout, '(T8,A,T31,  A  ,T50,A    ,T59,A)')                           &
              'ylm_form_write',TRIM(ylm_form_write), ylm_form_write_d, ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                            'ytunit_in',                 ytunit_in,    ' C '
  WRITE (nuout, '(T8,A,T21,  A            ,T59,A)')                           &
                            'ytunit_out',                ytunit_out,   ' C '
  WRITE (nuout, '(T8,A,T30,A)')                                           &
                         'yncglob_institution',   TRIM(yncglob_institution)
  WRITE (nuout, '(T8,A,T30,A)')                                           &
                         'yncglob_title',         TRIM(yncglob_title)
  WRITE (nuout, '(T8,A,T30,A)')                                           &
                         'yncglob_source',        TRIM(yncglob_source)
  WRITE (nuout, '(T8,A,T30,A)')                                           &
                         'yncglob_project_id',    TRIM(yncglob_project_id)
  WRITE (nuout, '(T8,A,T30,A)')                                           &
                         'yncglob_experiment_id', TRIM(yncglob_experiment_id)
  WRITE (nuout, '(T8,A,T30,A)')                                           &
                         'yncglob_contact',       TRIM(yncglob_contact)
  WRITE (nuout, '(T8,A,T30,A)')                                           &
                         'yncglob_references',    TRIM(yncglob_references)
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
            'ncglob_realization',ncglob_realization, ncglob_realization,' I '

  WRITE (nuout, '(A26)')  'Variables for initial data'
  DO i=1,numlist_ini
     WRITE (nuout, '(T7,A,I3,A,T21,A,T58,A)')                                &
                                     'yvarini(',i,')', youtlist_ini(i),'C10'
  ENDDO
  WRITE (nuout, '(A2)')  '  '

  WRITE (nuout, '(A27)')  'Variables for boundary data'
  DO i=1,numlist_bd
     WRITE (nuout, '(T7,A,I3,A,T21,A,T58,A)')                                &
                                     'yvarbd(',i,')', youtlist_bd(i),'C10'
  ENDDO
  WRITE (nuout, '(A2)')  '  '
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_data

!==============================================================================
!------------------------------------------------------------------------------
!+ Module procedure in "setup" for the input of NAMELIST prictr
!------------------------------------------------------------------------------

SUBROUTINE input_prictr (realbuf, intbuf, logbuf, charbuf, ibuflen,          &
                         nuout, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group prictr. 
!   The group prictr contains variables which indicate whether some debug or
!   control output should be done. Also it contains the variables for the
!   output of profiles.
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for 
!   consistency. If wrong input values are detected the program prints 
!   an error message. The program is not stopped in this routine but an 
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists. 
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.    
!   Both, default and input values are written to the file OUTPUT.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    ibuflen,      & ! dimension of the buffers
    nuout,        & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Sending buffers (are given to routine distribute_values)
  INTEGER (KIND=iintegers), INTENT(INOUT)    ::  intbuf  (ibuflen)
  REAL    (KIND=ireals)   , INTENT(INOUT)    ::  realbuf (ibuflen)
  LOGICAL                 , INTENT(INOUT)    ::  logbuf  (ibuflen)
  CHARACTER (LEN=100)     , INTENT(INOUT)    ::  charbuf (ibuflen)

! Local variables

! Variables for default values
  INTEGER (KIND=iintegers)   ::       &
    nlev1pr_d,          & ! k-index for printing the first  model layer
    nlev2pr_d,          & ! k-index for printing the second model layer
    igp_tot_d(nmaxgp),  & ! i-index for printing selected grid points
    jgp_tot_d(nmaxgp)     ! j-index for printing selected grid points

  LOGICAL                    ::       &
    lprps_d,       & ! if .TRUE., print some ps- and fis-fields
    lprt_d,        & ! if .TRUE., print t at 2 levels (nlev1pr,nlev2)
    lpru_d,        & ! if .TRUE., print u at 2 levels
    lprv_d,        & ! if .TRUE., print v at 2 levels
    lprgrh_d,      & ! if .TRUE., print general relative humidity  at 2 levels
    lprqv_d,       & ! if .TRUE., print qv at 2 levels
    lprqc_d,       & ! if .TRUE., print qc at 2 levels
    lprud_d,       & ! if .TRUE., print ud (divergent wind correction)
    lprvd_d,       & ! if .TRUE., print vd (divergent wind correction)
    lprdpdt_d,     & ! if .TRUE., print dpdt (tendency of surface pressure)
    lprgp_d,       & ! if .TRUE., print profiles at selected grid points
    lchkin_d,      & ! if .TRUE., print check-values of input-fields
    lchkout_d        ! if .TRUE., print check-values of output-fields

  INTEGER (KIND=iintegers)   ::       &
    nzgp, nztest, nsend, n, ierr, ierror

  CHARACTER (LEN=15)         :: yroutine

! Define the namelist group
  NAMELIST /prictr/ &
    nlev1pr,     & ! k-index for printing the first  model layer
    nlev2pr,     & ! k-index for printing the second model layer
    igp_tot,     & ! i-index for printing selected grid points
    jgp_tot,     & ! j-index for printing selected grid points
    lprps,       & ! if .TRUE., print some ps- and fis-fields
    lprt,        & ! if .TRUE., print t at 2 levels (nlev1pr,nlev2)
    lpru,        & ! if .TRUE., print u at 2 levels
    lprv,        & ! if .TRUE., print v at 2 levels
    lprgrh,      & ! if .TRUE., print general relative humidity  at 2 levels
    lprqv,       & ! if .TRUE., print qv at 2 levels
    lprqc,       & ! if .TRUE., print qc at 2 levels
    lprud,       & ! if .TRUE., print ud (divergent wind correction)
    lprvd,       & ! if .TRUE., print vd (divergent wind correction)
    lprdpdt,     & ! if .TRUE., print dpdt (tendency of surface pressure)
    lprgp,       & ! if .TRUE., print profiles at selected grid points
    lchkin,      & ! if .TRUE., print check-values of input-fields
    lchkout        ! if .TRUE., print check-values of output-fields

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_prictr
!------------------------------------------------------------------------------

yroutine = 'input_prictr'
ierror = 0

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  nlev1pr_d = 10
  nlev2pr_d = 20
  lprps_d   = .FALSE.
  lprt_d    = .FALSE.
  lpru_d    = .FALSE.
  lprv_d    = .FALSE.
  lprgrh_d  = .FALSE.
  lprqv_d   = .FALSE.
  lprqc_d   = .FALSE.
  lprud_d   = .FALSE.
  lprvd_d   = .FALSE.
  lprdpdt_d = .FALSE.
  lprgp_d   = .FALSE.
  lchkin_d  = .FALSE.
  lchkout_d = .FALSE.

  ngp_tot   = 7
  igp_tot_d(1) = 10
  jgp_tot_d(1) = 10
  igp_tot_d(2) = 10
  jgp_tot_d(2) = 30
  igp_tot_d(3) = 26
  jgp_tot_d(3) =  8
  igp_tot_d(4) = 17
  jgp_tot_d(4) = 36
  jgp_tot_d(5) = 29
  igp_tot_d(5) = 20
  igp_tot_d(6) = 17
  jgp_tot_d(6) = 10
  igp_tot_d(7) = 40
  jgp_tot_d(7) = 20

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with default
!------------------------------------------------------------------------------

  nlev1pr   = nlev1pr_d
  nlev2pr   = nlev2pr_d
  lprps     = lprps_d
  lprt      = lprt_d
  lpru      = lpru_d 
  lprv      = lprv_d 
  lprgrh    = lprgrh_d 
  lprqv     = lprqv_d 
  lprqc     = lprqc_d 
  lprud     = lprud_d 
  lprvd     = lprvd_d 
  lprdpdt   = lprdpdt_d 
  lprgp     = lprgp_d 
  lchkin    = lchkin_d
  lchkout   = lchkout_d

  DO n = 1,nmaxgp
    igp_tot(n) = 0
    jgp_tot(n) = 0
  ENDDO

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  READ (nuin, prictr, IOSTAT=ierror)
ENDIF

! distribute error status to all processors
IF (nproc > 1) THEN
  CALL distribute_values  (ierror, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF
IF (ierror /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

  ! Check whether grid point indices have been read; otherwise they are
  ! initialized with the defaults
  nztest = 0
  nzgp   = 0

  DO n = 1,nmaxgp
    nztest = nztest + igp_tot(n)
    IF ( igp_tot(n) /= 0 ) THEN
      nzgp = nzgp+1
    ENDIF
  ENDDO

  IF (nztest == 0) THEN
    ! no grid point indices have been read; use default values
    DO n = 1,ngp_tot
      igp_tot(n) = igp_tot_d(n)
      jgp_tot(n) = jgp_tot_d(n)
    ENDDO
  ELSE
    ! number of grid points
    ngp_tot   = nzgp
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf  ( 1) = nlev1pr
    intbuf  ( 2) = nlev2pr
    intbuf  ( 3) = ngp_tot

    ! indices of selected grid points
    DO n = 1,ngp_tot
      intbuf (3 + n)           = igp_tot (n)
      intbuf (3 + ngp_tot + n) = jgp_tot (n)
    ENDDO

    logbuf  ( 1) = lprps
    logbuf  ( 2) = lprt
    logbuf  ( 3) = lpru
    logbuf  ( 4) = lprv
    logbuf  ( 5) = lprgrh
    logbuf  ( 6) = lprqv
    logbuf  ( 7) = lprqc
    logbuf  ( 8) = lprud
    logbuf  ( 9) = lprvd
    logbuf  (10) = lprdpdt
    logbuf  (11) = lprgp
    logbuf  (12) = lchkin
    logbuf  (13) = lchkout
  ENDIF
  nsend = 3 + 2 * nmaxgp

  CALL distribute_values  (intbuf , nsend, 0, imp_integers, icomm_world, ierr)
  CALL distribute_values  (logbuf ,    13, 0, imp_logical,  icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    nlev1pr     = intbuf  ( 1)
    nlev2pr     = intbuf  ( 2)
    ngp_tot     = intbuf  ( 3)

    ! indices of selected grid points
    DO n = 1,ngp_tot
      igp_tot (n) = intbuf (3 + n)
      jgp_tot (n) = intbuf (3 + ngp_tot + n)
    ENDDO

    lprps       = logbuf  ( 1)
    lprt        = logbuf  ( 2)
    lpru        = logbuf  ( 3)
    lprv        = logbuf  ( 4)
    lprgrh      = logbuf  ( 5)
    lprqv       = logbuf  ( 6)
    lprqc       = logbuf  ( 7)
    lprud       = logbuf  ( 8)
    lprvd       = logbuf  ( 9)
    lprdpdt     = logbuf  (10)
    lprgp       = logbuf  (11)
    lchkin      = logbuf  (12)
    lchkout     = logbuf  (13)
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(A23)') '      NAMELIST:  prictr'
  WRITE (nuout, '(A23)') '      -----------------'
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value',       &
                                               'Default Value', 'Format'

  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                 'nlev1pr',   nlev1pr,    nlev1pr_d ,  ' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                 'nlev2pr',   nlev2pr,    nlev2pr_d ,  ' I '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprps',     lprps,      lprps_d,     ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprt' ,     lprt,       lprt_d,      ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lpru' ,     lpru,       lpru_d,      ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprv' ,     lprv,       lprv_d,      ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprgrh',    lprgrh,     lprgrh_d,    ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprqv',     lprqv,      lprqv_d,     ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprqc',     lprqc,      lprqc_d,     ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprud',     lprud,      lprud_d,     ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprvd',     lprvd,      lprvd_d,     ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprdpdt',   lprdpdt,    lprdpdt_d,   ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lprgp',     lprgp,      lprgp_d,     ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lchkin',    lchkin,     lchkin_d,    ' L '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lchkout',   lchkout,    lchkout_d,   ' L '
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(A38)') '      Indices of selected grid points:'
  WRITE (nuout, '(T9,A,T18,A,T28,A,T56,A,T66,A,T76,A)')                       &
                 'Nr.: ','i-index','j-index','Nr.:','i-index', 'j-index'
  DO n = 0,ngp_tot-1,2
    IF((ngp_tot-n) == 1) THEN
      WRITE (nuout,'(T10,I2,T19,I3,T29,I3)') ngp_tot, igp_tot(ngp_tot),       &
                                                         jgp_tot(ngp_tot)
    ELSE
      WRITE (nuout,'(T10,I2,T19,I3,T29,I3,T57,I2,T67,I3,T77,I3)')             &
       n+1, igp_tot(n+1), jgp_tot(n+1), n+2, igp_tot(n+2), jgp_tot(n+2)
    ENDIF
  ENDDO

  WRITE (nuout, '(A2)')  '  '
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_prictr

!==============================================================================
!------------------------------------------------------------------------------
!+ Module procedure in "setup" for the input of NAMELIST epsctl
!------------------------------------------------------------------------------

SUBROUTINE input_epsctl (realbuf, intbuf, logbuf, charbuf, ibuflen,          &
                         nuout, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group epsctl. 
!   The group epsctl contains variables which charcterize the ensemble of
!   boundary data (ensemble typ, member ID, total number of members).
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for 
!   consistency. If wrong input values are detected the program prints 
!   an error message. The program is not stopped in this routine but an 
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists. 
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.    
!   Both, default and input values are written to the file OUTPUT.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    ibuflen,      & ! dimension of the buffers
    nuout,        & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Sending buffers (are given to routine distribute_values)
  INTEGER (KIND=iintegers), INTENT(INOUT)    ::  intbuf  (ibuflen)
  REAL    (KIND=ireals)   , INTENT(INOUT)    ::  realbuf (ibuflen)
  LOGICAL                 , INTENT(INOUT)    ::  logbuf  (ibuflen)
  CHARACTER (LEN=100)     , INTENT(INOUT)    ::  charbuf (ibuflen)

! Local variables

! Variables for default values
  INTEGER (KIND=iintegers)   ::       &
    iepsmem_bc_d,       & ! ID of the member in the ensemble of boundary
                          ! conditions (iepsmem_bc >= 0)
    iepstyp_bc_d,         & ! ID of the boundary ensemble generation type
                          ! (iepstyp_bc >= 0)
    iepstot_bc_d            ! total number of boundary ensemble members (iepstot_bc >= 0)
   
  LOGICAL                    ::       &
    lchk_bc_typ_d         ! if .TRUE., check member ID of input data (if leps_bc=.TRUE.)

  INTEGER (KIND=iintegers)   ::       &
    ierr, ierror

  CHARACTER (LEN=15)         :: yroutine

! Define the namelist group
  NAMELIST /epsctl/ &
    iepsmem_bc,    & ! ID of the member in the ensemble of boundary
                     ! conditions (iepsmem_bc >= 0)
    iepstyp_bc,    & ! ID of the boundary ensemble generation type
                     ! (iepstyp_bc >= 0)
    iepstot_bc,    & ! total number of boundary ensemble members (iepstot_bc >= 0)
    lchk_bc_typ      !  if .TRUE., check member ID of input data (if leps_bc=.TRUE.)

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_epsctl
!------------------------------------------------------------------------------

yroutine = 'input_epsctl'
ierror = 0

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  iepsmem_bc_d  = -(1_iintegers)
  iepstyp_bc_d  = -(1_iintegers)
  iepstot_bc_d  = 0_iintegers
  lchk_bc_typ_d = .FALSE.

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with default
!------------------------------------------------------------------------------

  iepsmem_bc  = iepsmem_bc_d
  iepstyp_bc  = iepstyp_bc_d
  iepstot_bc  = iepstot_bc_d
  lchk_bc_typ = lchk_bc_typ_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  READ (nuin, epsctl, IOSTAT=ierror)
ENDIF

! distribute error status to all processors
IF (nproc > 1) THEN
  CALL distribute_values  (ierror, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF
IF (ierror /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

! Check whether values are OK
  IF ( iepsmem_bc < 0 .OR. iepstot_bc < 0 .OR. &
       iepstyp_bc < 0 ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** iepsmem_bc, iepstot_bc, iepstyp_bc *** ', &
       iepsmem_bc,iepstot_bc,iepstyp_bc
    ierrstat = 1002
    RETURN
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf  ( 1) = iepsmem_bc
    intbuf  ( 2) = iepstyp_bc
    intbuf  ( 3) = iepstot_bc
    logbuf  ( 1) = lchk_bc_typ
  ENDIF

  CALL distribute_values  (intbuf , 3, 0, imp_integers, icomm_world, ierr)
  CALL distribute_values  (logbuf ,    1, 0, imp_logical,  icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    iepsmem_bc     = intbuf  ( 1)
    iepstyp_bc     = intbuf  ( 2)
    iepstot_bc     = intbuf  ( 3)
    lchk_bc_typ    = logbuf  ( 1) 
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(A23)') '      NAMELIST:  epsctl'
  WRITE (nuout, '(A23)') '      -----------------'
  WRITE (nuout, '(A2)')  '  '
  WRITE (nuout, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value',       &
                                               'Default Value', 'Format'

  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                 'iepsmem_bc',iepsmem_bc, iepsmem_bc_d  ,  ' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                 'iepstot_bc',iepstot_bc, iepstot_bc_d  ,  ' I '
  WRITE (nuout, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                          &
                                 'iepstyp_bc',iepstyp_bc, iepstyp_bc_d  ,  ' I '
  WRITE (nuout, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                          &
                                 'lchk_bc_typ',lchk_bc_typ, lchk_bc_typ_d   ,' L '
  WRITE (nuout, '(A2)')  '  '
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_epsctl

!==============================================================================

END MODULE src_namelists
