!+ Source Module for reading input-fields for the coarse grid
!==============================================================================

MODULE src_read_hhl

!==============================================================================
!
! Description:
!   This module contains subroutines necessary to read the HHL file for the
!   COSMO and for the input coarse grid
!
! Current Code Owner: DWD, Ulrich Schaettler
!    phone:  +49  69  8062 2739
!    fax:    +49  69  8062 3721
!    email:  ulrich.schaettler@dwd.de
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V1_22        2013/07/11 Ulrich Schaettler
!  Initial release for INT2LM
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
  iwlength,  & ! length of an integer word in byte
  iintegers, & ! KIND-type parameter for "normal" integer variables
  intgribf,  & ! KIND-type of the fortran decks in the grib library
  int_ga       ! integer precision for grib_api: length of message in bytes

!------------------------------------------------------------------------------

USE data_int2lm_control,    ONLY: &
  ltime,        & ! detailed timings of the program are given
  timings,      & ! for storing the times for different parts of the program
  idbg_level,   & ! to control verbosity of output
  lprintdeb_all   ! whether all PEs print debug output

!------------------------------------------------------------------------------

USE data_int2lm_io,        ONLY : &
  lchkin,       & ! logical for print of check-values (max,min,mean)
                  ! of input-fields
  nuchkdat,     & ! checking the I/O data
  yuchkdat,     & ! checking the I/O data
  lfd,          & !
  lds,          & !
  undefgrib,    & ! value for "undefined" in the grib routines
  iblock          ! array for gribed data

!------------------------------------------------------------------------------

USE data_int2lm_parallel,      ONLY :  &
    nproc,           & ! total number of processors: nprocx * nprocy
    lasync_io,       & ! if .TRUE.: the model runs with extra PEs for
                       ! asynchronous IO
    num_compute,     & ! number of compute PEs
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    nboundlines,     & ! number of boundary lines of the domain for which
    icomm_cart,      & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_character,   & ! determines the correct CHARACTER type used in the model
                       ! for MPI
    imp_integers       ! determines the correct INTEGER type used in the model
                       ! for MPI

!------------------------------------------------------------------------------

USE utilities,            ONLY : elapsed_time
USE vgrid_refatm_utils,   ONLY : vcoord_type, imax_vcoordtype, vcoord_defaults, &
                                 uuid_2char, uuid_in_string, uuid_out_string

USE environment,          ONLY : comm_barrier, model_abort
USE parallel_utilities,   ONLY : distribute_values, global_values,              &
                                 scatter_values
USE io_utilities,         ONLY : open_file, close_file, check_record,           &
                                 read_gribapi

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================
!+ Organizes the reading and interpolation of coarse grid fields
!------------------------------------------------------------------------------

SUBROUTINE org_read_hhl  (hhl, idim, jdim, kdim, idim_tot, jdim_tot,        &
                          startlat_tot, startlon_tot, vc_type,              &
                          nsubpos, ydir, yfname, yfile)

!------------------------------------------------------------------------------
!
! Description:
!  HHL for the COSMO-grid (argument yfile='lm') or for the input-grid (yfile='in')
!  is read and distributed to the available PEs. In case of COSMO-grid there is 
!  the speciality, that only the COSMO-domain without boundary lines is read.
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
INTEGER(KIND=iintegers),  INTENT(IN)    ::   &
  idim, jdim, kdim,      & ! dimensions for local hhl field
  idim_tot, jdim_tot,    & ! dimensions for global hhl field
  nsubpos(0:num_compute-1,4) ! decomposition of grid

TYPE(vcoord_type),        INTENT(INOUT) ::   &
  vc_type

REAL(KIND=ireals),        INTENT(IN)    ::   &
  startlat_tot, startlon_tot ! to check horizontal grid

REAL(KIND=ireals),        INTENT(OUT)   ::   &
  hhl(idim,jdim,kdim)      ! height of half levels

CHARACTER(LEN=*),         INTENT(IN)    ::   &
  ydir,                  & ! directory
  yfname,                & ! name of the file
  yfile                    ! to indicate whether COSMO or Input-files are read

!------------------------------------------------------------------------------

! Local arrays
REAL (KIND=ireals)         ::                           &
  field_local     (idim    ,jdim    ),                  &
  field_api       (idim_tot*jdim_tot),                  &
  field_glob      (idim_tot,jdim_tot),                  &
  field_temp      (idim_tot+2,jdim_tot+2)

REAL (KIND=ireals), ALLOCATABLE         ::                           &
  zprocarray      (:,:)

REAL (KIND=irealgrib)      ::                           &
  field_grib      (idim_tot,jdim_tot)

REAL (KIND=ireals)         :: zundef, rla1, rlo1, zrealdiff

! Local variables:
INTEGER  (KIND=intgribf)   :: &
  igribid, ireturn, igriblen, ilevel, itoplevel, ibottomlevel

INTEGER  (KIND=iintegers)  ::  &
  nufile,                 & ! unit number of opened grib file
  izerror, izstat,        & ! status and error status variable
  niostat,                & ! status and error status variable
  izdebug,                & ! for verbosity of output
  iztestdim, jztestdim,   & ! to check meta data
  izee_hhl,               & ! location of record in variable table
  i,j,k,n, ij, lfdec, ldsec, iproc, ie_p, je_p, igrbednr, &
  izexch, izaee, nzlevels, izvctype, idim_max, jdim_max

INTEGER (KIND=int_ga)      ::  &
  izmaxlen, izsize_ga

CHARACTER(LEN=1) :: yvc_uuid(16)

LOGICAL                    ::  &
  lzrequired, & ! indicates whether a record from the grib file is required
  lzeof,      & ! indicates the end of file
  lzexist,    & ! to check, whether different files exist
  lzall,      & ! whether all levels have been read
  lzcheck(kdim) ! to check, whether all data have been read

CHARACTER (LEN=250)        ::  &
  yzpath      ! full path and name of the input-file

CHARACTER (LEN= 30)        ::  &
  ytypeoflevel ! typeOfLevel

CHARACTER (LEN= 21)        ::  &
  yshortname, yzname  ! short name from grib_api

CHARACTER (LEN= 25)        ::  &
  yzroutine   ! name of this routine for error handling

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling

INTEGER (KIND=iintegers)   ::  &
  intbuf      (3) ! an integer buffer for sending

CHARACTER (LEN=100)        ::  &
  charbuf

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  yzroutine   = 'org_read_hhl'
  yzerrmsg(:) = ' '
  izerror     = 0_iintegers
  izstat      = 0_iintegers
  lzcheck(:)  = .FALSE.

  igriblen = 0
  yzname   = '                    '
  izee_hhl = 8
  zundef   = REAL(undefgrib, ireals)

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
    PRINT *, '  Start to read HHL file '
  ENDIF

  ! initialize level counter and logical flags
  lzeof       = .FALSE.
  lzrequired  = .FALSE.
  lzexist     = .FALSE.

  ! Determine maximum for all idim,jdim -values from all tasks
  intbuf(1) = idim
  intbuf(2) = jdim
  CALL global_values (intbuf, 2, 'MAX', imp_integers, icomm_cart, -1,    &
                        yzerrmsg, izerror)
  idim_max = intbuf(1)
  jdim_max = intbuf(2)

  ALLOCATE (zprocarray((idim_max+1)*(jdim_max+1), 0:num_compute-1), STAT=izerror)


  ldsec = idim_tot * jdim_tot
  lfdec = (ldsec *8) / iwlength  + 5000
  lfd     = INT (lfdec, intgribf)
  lds     = INT (ldsec, intgribf)

  ALLOCATE (iblock(lfd),    STAT=izstat)

!-------------------------------------------------------------------------------
! Section 2: Create grib file name and open the file
!-------------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '  Open the HHL file'
  ENDIF

  ! Test whether the file exists
  yzpath = TRIM(ydir)//TRIM(yfname)
  INQUIRE (FILE=TRIM(yzpath), EXIST=lzexist)

  IF (.NOT. lzexist) THEN
    PRINT *,   ' HHL File does not exist:  ', TRIM(yzpath)
    yzerrmsg = ' *** ERROR:  HHL file for COSMO input does not exist *** '
    izerror  = 5
    CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
  ELSE
    ! All processors have to call the routine open_file. What the parallel
    ! program really does is determined in the routine.
    CALL open_file(nufile, yzpath, 'r  ', 'apix', icomm_cart,         &
                   my_cart_id, num_compute, lasync_io, idbg_level,    &
                   yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, izerror, yzerrmsg, 'open_file')
    ELSE 
      lzeof = .FALSE.
    ENDIF
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
    WRITE (nuchkdat,'(A,A)') 'Check HHL input file:  ', TRIM(yzpath)
    WRITE (nuchkdat,'(A,I5,A,I5)')                                       &
      '    idim_tot = ',idim_tot,'  jdim_tot = ',jdim_tot
    WRITE (nuchkdat,'(A)') '    '
    WRITE (nuchkdat,'(A,A)')                                             &
              '   var        ee  lev        min    ',                    &
              'imin jmin               max    imax jmax              mean'
    WRITE (nuchkdat,'(A)') '  '
  ENDIF

!-------------------------------------------------------------------------------
! Section 3: (Endless) loop over all records in the grib file
!-------------------------------------------------------------------------------
 
  IF (izdebug > 10) THEN
    PRINT *, '  Start read_loop'
  ENDIF

  read_loop: DO WHILE (.NOT. lzeof)

  !-----------------------------------------------------------------------------
  ! 3.1: Get a record
  !-----------------------------------------------------------------------------

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (12) = timings(12) + zrealdiff
    ENDIF

    ! Every PE gets one record from the file in an ordered manner
    ! (determined by the rank in the communicator icomm_rank). How this is
    ! done exactly is determined in the routine read_"format". This routine has
    ! to be called by every PE.

    IF (izdebug >= 20) THEN
      PRINT *, '      Calling read_gribapi'
    ENDIF

    izmaxlen  = INT (lfdec*iwlength, int_ga)
    CALL read_gribapi (nufile, izmaxlen, lfdec, icomm_cart, iblock, izsize_ga, &
                       num_compute, lasync_io, yzerrmsg, izerror)
    IF (izerror /= 0) THEN
       CALL model_abort (my_cart_id, 2013, yzerrmsg, yzroutine)
    ENDIF

    IF (idbg_level >= 20) THEN
      IF (izsize_ga == 0) THEN
        PRINT *, '       EOF readched'
      ELSE
        PRINT *, '       Read a record with size ', izsize_ga
      ENDIF
    ENDIF

    IF (izsize_ga == 0) THEN
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
  ! 3.2: Unpack and convert the record
  !-----------------------------------------------------------------------------

    IF (.NOT. lzeof) THEN
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

      ! Get necessary meta data
      CALL grib_get (igribid, 'shortName',     yshortname, ireturn)
      CALL grib_get (igribid, 'typeOfLevel',   ytypeoflevel,  ireturn)
      CALL grib_get (igribid, 'level',         ilevel,        ireturn)
      CALL grib_get (igribid, 'topLevel',      itoplevel,     ireturn)
      CALL grib_get (igribid, 'bottomLevel',   ibottomlevel,  ireturn)

      IF (ireturn /= GRIB_SUCCESS) THEN
        PRINT *,   ' *** Error in grib_api grib_get keys  ', ireturn
        yzerrmsg = ' *** Error in grib_api grib_get keys  '
        CALL model_abort (my_cart_id, 2036, yzerrmsg, yzroutine)
      ELSE
        IF (izdebug > 20) THEN
          PRINT *, ' got data: ', yshortname, iztestdim, jztestdim, ytypeoflevel, ilevel, itoplevel, ibottomlevel
        ENDIF
      ENDIF

      ! Check the meta data
      IF (TRIM(yshortname) == 'HHL' .AND. TRIM(ytypeoflevel) == 'generalVertical') THEN
        lzrequired = .TRUE.

        ! check the horizontal grid
        CALL grib_get (igribid, 'Ni',            iztestdim,     ireturn)
        CALL grib_get (igribid, 'Nj',            jztestdim,     ireturn)
        CALL grib_get (igribid, 'latitudeOfFirstGridPointInDegrees',  rla1,  ireturn)
        CALL grib_get (igribid, 'longitudeOfFirstGridPointInDegrees', rlo1,  ireturn)
        IF (rlo1 >= 180.0_ireals) THEN
          rlo1 = REAL(rlo1, ireals) - 360.0_ireals
        ENDIF

        IF (iztestdim /= idim_tot) izerror = 1
        IF (jztestdim /= jdim_tot) izerror = 2
        IF (ABS(startlat_tot - rla1) > 1E-5_ireals) izerror = 3
        IF (ABS(startlon_tot - rlo1) > 1E-5_ireals) izerror = 4
        IF (izerror > 0) THEN
          WRITE (*,'(A,A)')      ' Error in the horizontal grid of HHL:  ', TRIM(yfname)
          WRITE (*,'(A,A)')      '             Model values    File values'
          WRITE (*,'(A,2I13)')   ' iedim:    ', idim_tot, iztestdim
          WRITE (*,'(A,2I13)')   ' jedim:    ', jdim_tot, jztestdim
          WRITE (*,'(A,2F13.4)') ' startlat: ', startlat_tot, rla1
          WRITE (*,'(A,2F13.4)') ' startlon: ', startlon_tot, rlo1
          CALL model_abort (my_cart_id, 2036, 'wrong grid for HHL', yzroutine)
        ENDIF

        ! get the PE with smallest ID to exchange necessary meta data
        izexch = my_cart_id
      ELSE
        lzrequired = .FALSE.
        izexch = nproc
      ENDIF
    ELSE   ! leof
      lzrequired = .FALSE.
      izexch = nproc

      IF (idbg_level > 10) THEN
        PRINT *, '  End of File reached'
      ENDIF
    ENDIF  ! leof

    ! Get the minimum of all values in izexch
    IF (num_compute > 1) THEN
      CALL global_values (izexch, 1, 'MIN', imp_integers, icomm_cart,  &
                          -1, yzerrmsg, ireturn)
    ENDIF

    IF (izexch < nproc) THEN
      ! The processor with lowest id, that has got a HHL field
      ! and distributes the meta data to all others

      IF (izexch == my_cart_id) THEN
        CALL grib_get (igribid, 'nlev',                nzlevels)
        CALL grib_get (igribid, 'numberOfVGridUsed',   izvctype)
        CALL grib_get (igribid, 'uuidOfVGrid',         yvc_uuid)

        ! get localInformationNumber => izaee
        CALL grib_get (igribid, 'localInformationNumber', izaee)

        intbuf(1) = nzlevels
        intbuf(2) = izvctype
        intbuf(3) = izaee
        DO i = 1, 16
          charbuf(i:i) = yvc_uuid(i)
        ENDDO
      ENDIF

      ! First distribute number of vertical coordinate parameters
      IF (num_compute > 1) THEN
        CALL distribute_values (intbuf,  3, izexch, imp_integers,  icomm_cart, ireturn)
        CALL distribute_values (charbuf, 1, izexch, imp_character, icomm_cart, ireturn)
      ENDIF

      IF (izexch /= my_cart_id) THEN
        nzlevels = intbuf(1)
        izvctype = intbuf(2)
        izaee    = intbuf(3)
        DO i = 1, 16
          yvc_uuid(i) = charbuf(i:i)
        ENDDO
      ENDIF
    ENDIF

    IF     (yfile == 'lm') THEN
      CALL uuid_2char (yvc_uuid, uuid_out_string)
      ! print *, 'got meta data:  ', nzlevels, izvctype, izaee, uuid_out_string
    ELSEIF (yfile == 'in') THEN
      CALL uuid_2char (yvc_uuid, uuid_in_string)
      ! print *, 'got meta data:  ', nzlevels, izvctype, izaee, uuid_in_string
    ENDIF


    ! Now check the vertical coordinate parameters
    IF (vc_type%ivcoord_id == -1) THEN
      ! no values set for vertical coordinates and up to now no default set could 
      ! be found: check whether there is a default set for izaee
      IF (izaee <= imax_vcoordtype) THEN
        ! set the defaults to vc_type
        vc_type%ivctype    = vcoord_defaults(izaee)%ivctype
        vc_type%ivcoord_id = vcoord_defaults(izaee)%ivcoord_id
        vc_type%nlevels    = vcoord_defaults(izaee)%nlevels
        vc_type%vc_uuid(:) = yvc_uuid(:)
        vc_type%vcflat     = vcoord_defaults(izaee)%vcflat
        IF      (vc_type%ivctype == 1) THEN
          vc_type%sigm_coord(1:vc_type%nlevels) = vcoord_defaults(izaee)%sigm_coord(1:vc_type%nlevels)
        ELSEIF ((vc_type%ivctype == 2) .OR. (vc_type%ivctype == 3)) THEN
          vc_type%vert_coord(1:vc_type%nlevels) = vcoord_defaults(izaee)%vert_coord(1:vc_type%nlevels)
        ENDIF
      ELSE
        WRITE (*,'(A)') ' *** ERROR: NO default vertical coordinate parameters available: ', izaee
        CALL model_abort (my_cart_id, 2036, 'No default vertical parameters available', yzroutine)
      ENDIF
    ELSE
      ! the default set vc_type%ivcoord_id is used for the vertical
      ! coordinate: check whether this is compatible with izaee
      IF (izaee /= vc_type%ivcoord_id) THEN
        WRITE (*,'(A)') ' *** ERROR: The set for vertical coordinates used in HHL-file'
        WRITE (*,'(A)') ' ***        and in the namelist INPUT do not match'
        CALL model_abort (my_cart_id, 2036, 'wrong meta data for HHL-file', yzroutine)
      ELSE
        ! set the uuid to vc_type
        vc_type%vc_uuid(:) = yvc_uuid(:)
      ENDIF
    ENDIF

    ! Perhaps this is not really necessary:
    IF ((nzlevels /= vc_type%nlevels) .OR. (izvctype /= vc_type%ivctype)) THEN
      WRITE (*,'(A,A)')      '             Namelist Input    File values'
      WRITE (*,'(A,2I13)')   ' nlevels   ', nzlevels, vc_type%nlevels
      WRITE (*,'(A,2I13)')   ' ivctype   ', izvctype, vc_type%ivctype
      CALL model_abort (my_cart_id, 2036, 'wrong meta data for HHL-file', yzroutine)
    ENDIF

    IF (lzrequired) THEN
      ! Get size of data and data itself
      ! Set missing_value before
      CALL grib_set (igribid, 'missingValue', zundef)
      CALL grib_get_size (igribid, 'values', igriblen, ireturn)
      IF (igriblen > ldsec) THEN
        PRINT *, ' *** ERROR: size of message is too big for allocated field: ', igriblen, ldsec
        yzerrmsg = ' *** ERROR: Wrong size of field ***'
        CALL model_abort (my_cart_id, 2037, yzerrmsg, yzroutine)
      ENDIF

      CALL grib_get (igribid, 'values', field_api, ireturn)
      IF (ireturn /= GRIB_SUCCESS) THEN
        PRINT *,   ' *** Error in grib_api grib_get values', ireturn
        yzerrmsg = ' *** Error in grib_api grib_get values'
        CALL model_abort (my_cart_id, 2038, yzerrmsg, yzroutine)
      ENDIF

      DO j = 1, jdim_tot
        DO i = 1, idim_tot
          ij = (j-1)*idim_tot + i
          field_glob(i,j) =      field_api(ij)
          ! we need this for check_record
          field_grib(i,j) = REAL(field_api(ij), irealgrib)
        ENDDO
      ENDDO

      IF (yfile == 'lm') THEN
        field_temp(:,:) = 0.0_ireals
        ! put a boundary line at each side
        ! this field has to be used for distribution
        DO j = 1, jdim_tot
          DO i = 1, idim_tot
            field_temp(i+1,j+1) = field_glob(i,j)
          ENDDO
        ENDDO
      ENDIF

      ! put values to zprocarray
      DO n=0,num_compute-1
        IF     (yfile == 'lm') THEN
          ie_p = nsubpos(n,3) - nsubpos(n,1) + 1 + 2*nboundlines
          je_p = nsubpos(n,4) - nsubpos(n,2) + 1 + 2*nboundlines
          DO j = 1, je_p
            DO i = 1, ie_p
              ij = (j-1)*(idim_max+1) + i
              zprocarray(ij,n) = field_temp(nsubpos(n,1)-nboundlines-1+i,   &
                                            nsubpos(n,2)-nboundlines-1+j)
            ENDDO
          ENDDO
        ELSEIF (yfile == 'in') THEN
          ie_p = nsubpos(n,3) - nsubpos(n,1) + 1
          je_p = nsubpos(n,4) - nsubpos(n,2) + 1
          DO j = 1, je_p
            DO i = 1, ie_p
              ij = (j-1)*(idim_max+1) + i
              zprocarray(ij,n) = field_glob(nsubpos(n,1)-1 + i,nsubpos(n,2)-1 + j)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
#endif
    ENDIF

  !----------------------------------------------------------------------------
  ! 3.3: Check the record
  !----------------------------------------------------------------------------

    ! print the maximum, minimum and meanvalues of each record
    IF (lchkin) THEN
      ! just to make this call save:
      CALL check_record (field_grib, 1, idim_tot, 1, jdim_tot, 1, 1,       &
               1, idim_tot, 1, jdim_tot, 1, 1, undefgrib,                  &
               yshortname(1:10)  , izee_hhl,         ilevel,               &
               (.NOT.lzeof) .AND. (lzrequired), nuchkdat, num_compute,     &
               icomm_cart, my_cart_id, yzerrmsg, izerror)
    ENDIF

    IF (ltime) THEN
      CALL elapsed_time (zrealdiff)
      timings (7) = timings(7) + zrealdiff
    ENDIF

  !----------------------------------------------------------------------------
  ! 3.4: Distribute record to all PEs and put values to memory
  !----------------------------------------------------------------------------

    distribute_loop: DO iproc = 0, num_compute-1

      ! The following routine handles the distribution of one record to all
      ! others. If the record is not required, the distribute_loop is cycled
      ! (return status izstat = -2), if all records are done, the read_loop
      ! is exited (return status izstat = -1).
      ! The sender is the processor with rank=iproc, all others receive the
      ! corresponding data into the structure zlocalarray.
      CALL scatter_hhl  (zprocarray, hhl, lzcheck, idim, jdim, kdim,       &
                         idim_max, jdim_max, iproc,                        &
                         yshortname, ilevel, lzeof, lzrequired, izstat)

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
  DEALLOCATE (iblock)

  ! close file
  CALL close_file (nufile, 'apix', icomm_cart, my_cart_id,         &
                   num_compute, lasync_io, idbg_level, yzerrmsg, izerror)
  IF (izerror /= 0) THEN
    CALL model_abort (my_cart_id, izerror, yzerrmsg, 'close_file')
  ENDIF

  IF ( (lchkin) .AND. (my_cart_id == 0) ) THEN
    CLOSE (nuchkdat, STATUS='KEEP')
  ENDIF

!-------------------------------------------------------------------------------
! Section 9: Check that all data are read
!-------------------------------------------------------------------------------

  ! Check, whether all data necessary are present and write information to 
  ! output
  IF (my_cart_id == 0) THEN
    IF (ALL(lzcheck(:)) ) THEN
      WRITE (*,'(A)') '  All levels of HHL have been read'
    ELSE
      WRITE (*,'(A)') '  NOT All levels of HHL could be read'
      DO k = 1, kdim
        IF (lzcheck(k)) THEN
          WRITE (*,'(A,I5,L5)') 'HHL: ****    ', k, lzcheck
        ELSE
          WRITE (*,'(A,I5,L5)') 'HHL:         ', k, lzcheck
        ENDIF
      ENDDO
      CALL model_abort (my_cart_id, 2038, 'HHL could not be read', yzroutine)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE org_read_hhl

!==============================================================================
!==============================================================================
!+ Subroutine for distributing the initial- and the boundary data
!------------------------------------------------------------------------------

SUBROUTINE scatter_hhl  (procarrays, hhl_field, ltest, idim, jdim, kdim,     &
                         idim_max, jdim_max, itask,                          &
                         yshname, ilev, leof, lrequired, istat)

!------------------------------------------------------------------------------
!
! Description:
!  This subroutine is called within the read-loop over all records (with the
!  loop index "itask"). In this loop up to num_compute (number of compute PE)
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
  idim, jdim,  kdim, ilev,           & ! characteristics of the variable
  idim_max, jdim_max,                & ! max of idim,jdim for all tasks
  itask                                ! index of the distribution loop

CHARACTER(LEN=*),         INTENT(IN)    ::  &
  yshname

LOGICAL,                  INTENT(INOUT) ::  &
  ltest(kdim)

LOGICAL,                  INTENT(IN)    ::  &
  leof, lrequired                      ! end of file and required field

REAL    (KIND=ireals),    INTENT(IN)    ::  &
  procarrays ((idim_max+1)*(jdim_max+1), 0:num_compute-1) ! decomposed total field

REAL    (KIND=ireals),    INTENT(INOUT) ::  &
  hhl_field(idim,jdim,kdim)

INTEGER (KIND=iintegers), INTENT(OUT)   ::  &
  istat                                ! go on, cycle or exit the loop

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)   :: implcode
INTEGER (KIND=iintegers)   :: iz_info(2), i,j, izstorelev
CHARACTER (LEN=25)         :: yzroutine
CHARACTER (LEN=80)         :: yzerrmsg

! Local arrays:
REAL    (KIND=ireals)      :: zsubarray_1d((idim_max+1)*(jdim_max+1)),    &
                              zsubarray_2d((idim_max+1),(jdim_max+1))
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine scatter_hhl
!------------------------------------------------------------------------------

  istat    = 0
  yzroutine = 'scatter_hhl'

  ! check and distribute the next action
  IF (itask == my_cart_id) THEN
    IF ( (.NOT.leof) .AND. (lrequired) ) THEN
      iz_info(1) = 0 ! Processor has data
      iz_info(2) = ilev
    ELSE
      IF (leof) THEN
        iz_info(1) = -1 ! No data because EOF reached
      ELSE
        iz_info(1) = -2 ! No data because it is not required
      ENDIF
    ENDIF
  ENDIF

  IF (num_compute > 1) THEN
    CALL distribute_values (iz_info, 2, itask, imp_integers, icomm_cart,   &
                            implcode)
    IF (implcode /= 0) THEN
      istat   = 1
      RETURN
    ENDIF
  ENDIF

  istat = iz_info(1)

  ! Distribute the records
  IF (istat == 0) THEN
    ! Update the values in the record var_in (ar_des_input)
    ltest(iz_info(2)) = .TRUE.

    IF (num_compute > 1) THEN
      CALL scatter_values (procarrays, zsubarray_1d, (idim_max+1)*(jdim_max+1), &
                           num_compute, imp_reals, itask, icomm_cart,           &
                           yzerrmsg, implcode)
      IF (implcode /= 0) THEN
        istat = 2
        RETURN
      ENDIF
    ELSE
      zsubarray_1d(:) = procarrays(:,0)
    ENDIF

    zsubarray_2d = RESHAPE (zsubarray_1d, (/idim_max+1, jdim_max+1/))

    IF ((idbg_level >= 20) .AND. (my_cart_id == 0) ) THEN
      PRINT *, '   Storing variable: ', yshname, iz_info(2)
    ENDIF

    ! put data into the variables
    hhl_field(1:idim,1:jdim,iz_info(2)) = zsubarray_2d(1:idim,1:jdim)
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE scatter_hhl

!==============================================================================
!==============================================================================

END MODULE src_read_hhl
