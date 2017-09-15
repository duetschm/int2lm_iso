!+ Main organization program of the interpolation to LM fields.
!==============================================================================

PROGRAM int2lm_org

!==============================================================================
!
! Description:
!   Main program: organizes the start up for the transformation
!   from coarse model fields to LM-fields fields.
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
!  Initial release
! V1_5         2007/07/09 Ulrich Schaettler
!  Added switch lyear_360
! V1_6         2007/09/07 Ulrich Schaettler
!  Introduced variables for actual date
! V1_8         2008/05/29 Ulrich Schaettler, Hans-Jürgen Panitz
!  Call to org_read_coarse grid and org_coarse_interpol also for lcm2lm
!  Deallocate memory for external parameters that are only written to initial file
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
! V1_9         2009/09/03 Ulrich Schaettler
!  Deallocate memory for external parameters after computing initial fields
! V1_10        2009/12/17 Ulrich Schaettler
!  Modifications to allow processing of Unified Model Data
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_12        2010/06/14 Ulrich Schaettler
!  Putting the call for the info-module to the right place
! V1_14        2010/11/19 Ulrich Schaettler
!  Modifications to allow processing of JMA and NCEP data
! V1_17        2011/03/11 Ulrich Schaettler
!  Added field fr_urban for urban fraction data (K. Trusilova)
! V1_19        2012/06/06 Ulrich Schaettler, Burkhardt Rockel
!  Added lhir2lm as internal logical flag
!  Call org_vert_interpol_p2h also for lcm_pres_coor
!  Call org_vert_inter_lm and do not call org_lm_fields in case of lcm_hgt_coor=.TRUE.
!  Deallocate eventual surface albedo fields after initial data processing
! V1_22        2013/07/11 Ulrich Schaettler, KIT
!  Implemented additional call to org_lm_output for writing HHL-file
!  Implemented conditional compilation for ART usage (KIT)
!  Allow call to org_2d_fields also for lgsm2lm
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Modules used:

! Load the library information data:
USE info_int2lm, ONLY: info_define, info_readnl, info_print

USE data_int2lm_control,  ONLY: nstart, nstop, nincbound, linitial,           &
                                lboundaries, lcomp_bound, dt, timings,        &
                                lgme2lm, lec2lm, llm2lm, lcm2lm, ndebug,      &
                                ltime, noutput, lbd_frame, lbd_frame_cur,     &
                                itype_calendar, yakdat1, yakdat2, lforest,    &
                                llake, lsso, lradtopo, lemiss, lstomata,      &
                                lum2lm, lgsm2lm, lgfs2lm, lhir2lm, lurban,    &
                                itype_albedo, l_art
USE data_fields_lm,       ONLY: grh_lm, ps_lm, for_e_lm, for_d_lm,            &
                                skyview_lm, slo_asp_lm, slo_ang_lm,           &
                                horizon_lm, sso_stdh_lm, sso_gamma_lm,        &
                                sso_theta_lm, sso_sigma_lm, fr_lake_lm,       &
                                depth_lk_lm, emis_rad_lm, prs_min_lm,         &
                                fr_urban_lm, alb_dry_lm, alb_sat_lm
USE data_grid_lm,         ONLY: ie2lm, je2lm, kedim
USE data_grid_in,         ONLY: lcm_hgt_coor, lcm_pres_coor
USE data_int2lm_io,       ONLY: numlist_ini, numlist_bd, numlist_hhl,         &
                                youtlist_ini, youtlist_bd, youtlist_hhl,      &
                                nmaxlist, ydate_ini, njulianday, ract_hour,   &
                                nunit_of_time, ylm_form_write
USE data_int2lm_parallel, ONLY: lcompute_pe, my_cart_id, my_world_id
USE data_parameters,      ONLY: ireals, iintegers


USE mpe_io,               ONLY: mpe_io_node, mpe_io_shutdown
USE environment,          ONLY: get_free_unit, final_environment, model_abort
USE utilities,            ONLY: get_utc_date, elapsed_time
USE vgrid_refatm_utils,   ONLY: lnewVGrid


USE src_2d_fields,          ONLY: org_2d_fields
USE src_cleanup,            ONLY: org_cleanup
USE src_coarse_interpol,    ONLY: org_coarse_interpol
USE src_gme_interpol,       ONLY: org_gme_interpol
USE src_lm_fields,          ONLY: org_lm_fields
USE src_lm_output,          ONLY: init_lm_output, org_lm_output
USE src_pressure_to_hybrid, ONLY: org_vert_interpol_p2h
USE src_read_coarse_grid,   ONLY: org_read_coarse_grid
USE src_vert_inter_lm,      ONLY: org_vert_inter_lm
USE src_vert_interpol,      ONLY: org_vert_interpol

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local Scalars:
INTEGER (KIND=iintegers) ::       &
  izerror, istat,         & ! error flag 
  izfiles,                & ! boundary file to be processed
  nnow                      ! for loop over all boundary fields
                            
REAL    (KIND=ireals)    ::       &
  zrealdiff,                & ! for time measuring
  gme_core

LOGICAL                    ::  &
  lzfirst                     ! for actions that have to be done only once

CHARACTER (LEN=200)      ::  yzerror    ! error message

!- End of header
!==============================================================================

izerror  = 0
yzerror  = '   '

!------------------------------------------------------------------------------
!- Section 1: Setup of the model
!------------------------------------------------------------------------------

CALL setup_int2lm (izerror, yzerror)

IF (my_cart_id == 0) THEN
  ! Print the default information to stdout:
  CALL info_define ('int2lm_org')        ! Pre-define the program name as binary name
  CALL info_readnl ('INPUT')             ! Read additional information from namelist file
  CALL info_print ('?')                  ! Print ALL the information to stdout
ENDIF

IF (izerror /= 0) THEN
  CALL model_abort (my_world_id, izerror, yzerror, 'int2lm_org')
ENDIF

#ifdef ART
IF (l_art) THEN
  CALL organize_art    ("setup     ", 0_iintegers, izerror)

  IF (izerror /= 0) THEN
    yzerror  = ' *** art module error *** '
    CALL model_abort (my_cart_id, 9999, yzerror , 'int2lm_org')
  ENDIF
ENDIF
#endif

#ifdef FPEABORT
! Floating point exception trapping
  CALL initialize_fpe_trap(.TRUE., istat)
  IF ( istat /= 0 ) THEN
    yzerror = 'Error initializing ftp_trap'
    CALL model_abort (my_world_id, istat, yzerror, 'int2lm_org')
  ENDIF
#endif

istat    = 0
IF (ltime) THEN
  CALL elapsed_time (zrealdiff)
  timings (2) = timings(2) + zrealdiff
ENDIF

!------------------------------------------------------------------------------
!- Section 2: Reading external parameters and setting LM reference atmosphere
!------------------------------------------------------------------------------

! Now comes the part for the compute PEs. The part for the IO-PEs is in
! the ELSE-part at the end of the program.

comp_pe: IF (lcompute_pe) THEN

  ! Open file OUTPUT for compute PEs
  IF (my_cart_id == 0) THEN
    ! get again a unit number for OUTPUT
    CALL get_free_unit (noutput)

    OPEN(noutput, FILE='OUTPUT', FORM=  'FORMATTED', STATUS='OLD',      &
         POSITION='APPEND', IOSTAT=istat)
    IF(istat /= 0) THEN
      yzerror  = ' ERROR    *** Error while opening file OUTPUT *** '
      CALL model_abort (my_cart_id, 1001, yzerror , 'intorg')
    ENDIF
  ENDIF

  CALL external_data (izerror, yzerror)

  IF (izerror /= 0) THEN
    CALL model_abort (my_world_id, izerror, yzerror , 'int2lm_org')
  ENDIF

  IF (ltime) THEN
    CALL elapsed_time (zrealdiff)
    timings (3) = timings(3) + zrealdiff
  ENDIF

!------------------------------------------------------------------------------
!- Section 3: Initialize output and do HHL-output, if necessary
!------------------------------------------------------------------------------

  ! Initialize the output
  CALL init_lm_output

  IF (lnewVGrid .AND. ylm_form_write == 'api2') THEN
    CALL org_lm_output (numlist_hhl, youtlist_hhl, 0, ydate_ini, .TRUE.)
  ENDIF

!------------------------------------------------------------------------------
!- Section 4: Loop over all boundary files
!------------------------------------------------------------------------------

  nnow          = nstart
  izfiles       = 0
  lzfirst       = .TRUE.

  loop_files: DO WHILE (nnow <= nstop)

  !----------------------------------------------------------------------------
  !- Section 4.1: Initialization of the loop
  !----------------------------------------------------------------------------

    izfiles = izfiles + 1

    ! check whether initial fields have to be computed first
    IF (linitial) THEN
      lcomp_bound = .FALSE.
      izfiles = 0
    ELSE
      lcomp_bound = .TRUE.
    ENDIF

    ! get new date
    CALL get_utc_date (nnow, ydate_ini, dt, itype_calendar, yakdat1, yakdat2, &
                       njulianday, ract_hour)

  !----------------------------------------------------------------------------
  !- Section 4.2: Reading and interpolation of coarse grid fields
  !----------------------------------------------------------------------------

    ! Set unit-of-time to an "undefined" value:
    ! it has to be read from the input fields
    nunit_of_time = -1

    IF (lgme2lm) THEN
      CALL org_gme_interpol (nnow, lzfirst, yakdat1)
    ELSE   ! all other models
      ! check if current output should include only frames
      IF (lbd_frame .AND. lcomp_bound) THEN
        lbd_frame_cur=.TRUE.
      ELSE
        lbd_frame_cur=.FALSE.
      ENDIF

      ! the fields for EC2LM and LM2LM are first all read
      CALL org_read_coarse_grid (nnow, lzfirst, yakdat1)

      IF (lgfs2lm .OR. (lcm2lm .AND. lcm_pres_coor)) THEN
        ! interpolate pressure levels to hybrid levels
        CALL org_vert_interpol_p2h (izerror, yzerror)
        IF (izerror /= 0) THEN
          yzerror  = ' ERROR    *** Error interpolating pressure to hybrid levels *** '
          CALL model_abort (my_cart_id, 1005, yzerror , 'intorg')
        ENDIF
      ENDIF

      ! and then interpolated
      CALL org_coarse_interpol
    ENDIF

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (8) = timings(8) + zrealdiff
    ENDIF

  !----------------------------------------------------------------------------
  !- Section 4.3: Vertical interpolation of atmospheric variables
  !----------------------------------------------------------------------------

    ! Allocation of additional LM fields:
    ALLOCATE (ps_lm  (ie2lm,je2lm),                                  &
              grh_lm (ie2lm,je2lm,kedim),            STAT = istat)

    IF (llm2lm .OR. lum2lm .OR. lcm_hgt_coor) THEN
      CALL org_vert_inter_lm
    ELSE
      CALL org_vert_interpol
    ENDIF

  !----------------------------------------------------------------------------
  !- Section 4.4: Additional two-dimensional fields
  !----------------------------------------------------------------------------

    ! For some input models this is skipped up to now
    IF ((.NOT. lum2lm) .AND. (.NOT. lhir2lm)) THEN
      CALL org_2d_fields (yakdat1, njulianday, ract_hour)
    ENDIF

  !----------------------------------------------------------------------------
  !- Section 4.5: Computation of LM-fields
  !----------------------------------------------------------------------------

    IF ((.NOT. llm2lm) .AND. (.NOT. lum2lm) .AND. (.NOT. lcm_hgt_coor)) THEN
      CALL org_lm_fields
    ENDIF

#ifdef ART
    ! Treatment of tracers
    IF (l_art) THEN
      CALL organize_art("loop      ", nnow, izerror)
      IF (izerror /= 0) THEN
        yzerror  = ' *** art module error *** '
        CALL model_abort (my_cart_id, 9999, yzerror , 'int2lm_org')
      ENDIF
    ENDIF
#endif

    DEALLOCATE (ps_lm, grh_lm)

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (10) = timings(10) + zrealdiff
    ENDIF

  !----------------------------------------------------------------------------
  !- Section 4.6: Writing HM/LM-results to disk
  !----------------------------------------------------------------------------

    IF (.NOT. lcomp_bound) THEN
      CALL org_lm_output (numlist_ini, youtlist_ini, nnow, yakdat1, .FALSE.)

      ! Deallocate some memory
      IF (llake) THEN
        DEALLOCATE (fr_lake_lm, depth_lk_lm)
      ENDIF
      IF (lforest) THEN
        DEALLOCATE  (for_e_lm, for_d_lm)
      ENDIF
      IF (lurban) THEN
        DEALLOCATE  (fr_urban_lm)
      ENDIF
      IF (lsso)    THEN
        DEALLOCATE  (sso_stdh_lm, sso_gamma_lm, sso_theta_lm, sso_sigma_lm)
      ENDIF
      IF (lradtopo) THEN
        DEALLOCATE (skyview_lm, slo_asp_lm, slo_ang_lm, horizon_lm)
      ENDIF
      IF (lemiss) THEN
        DEALLOCATE (emis_rad_lm)
      ENDIF
      IF (lstomata) THEN
        DEALLOCATE (prs_min_lm)
      ENDIF
      IF     (itype_albedo == 2) THEN
        DEALLOCATE (alb_dry_lm, alb_sat_lm)
        ! the fields from option itype_albedo=3 could also be used for boundary
        ! fields (in principle; not done yet); they are deallocated in src_cleanup
      ENDIF
    ELSE
      CALL org_lm_output (numlist_bd,  youtlist_bd,  nnow, yakdat1, .FALSE.)
    ENDIF

    ! Time measurement for output
    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (11) = timings(11) + zrealdiff
    ENDIF

  !----------------------------------------------------------------------------
  !- Section 4.7: Preparing the next loop cycle or exiting
  !----------------------------------------------------------------------------

    lzfirst = .FALSE.
    IF (.NOT. linitial) THEN
      nnow = nnow + nincbound
    ELSE
      ! initial fields have been computed in the first cycle;
      ! check, whether boundary fields have to be computed now:
      linitial = .FALSE.
      IF (.NOT. lboundaries) THEN
        EXIT loop_files
      ENDIF
    ENDIF

  ENDDO loop_files  ! loop over all boundary files

!------------------------------------------------------------------------------
!- Section 5: Cleanup
!------------------------------------------------------------------------------

  ! Cleanup of the compute PEs
  IF (my_cart_id == 0) THEN
    PRINT *, 'CLEANUP'
  ENDIF
  CALL org_cleanup

#ifdef ART
  IF (l_art) THEN
    CALL organize_art("cleanup   ", nnow, izerror)
    IF (izerror /= 0) THEN
      yzerror  = ' *** art module error *** '
      CALL model_abort (my_cart_id, 9999, yzerror , 'int2lm_org')
    ENDIF
  ENDIF
#endif

  ! Finish tasks of IO PEs
  CALL mpe_io_shutdown()

!------------------------------------------------------------------------------
!- Section 6: Part of the IO-PEs
!------------------------------------------------------------------------------

ELSE comp_pe

   CALL mpe_io_node()

ENDIF comp_pe

!------------------------------------------------------------------------------
! Section 7: Finalize the environment
!------------------------------------------------------------------------------

  CALL final_environment (izerror, yzerror)

!=======================================================================

END PROGRAM int2lm_org
