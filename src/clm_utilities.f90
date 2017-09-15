!+ Source module for clm_utility routines
!==============================================================================

MODULE  clm_utilities

!==============================================================================
!
! Description:
!   This module provides service utilities for the climate mode version. 
!
! Current Code Owner: Helmholtz-Zentrum Geesthacht
!  phone:  +49  4152 87 1803
!  phone:  +49  4152 87 8 1803
!  email:  burkhardt.rockel@hzg.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V1_19        2012/06/06 Burkhardt Rockel
!  Initial release for INT2LM
! V1_20        2012/09/03 Ulrich Schaettler
!  Enlarged strings for date variables to 14 characters: Adapted check for
!   existence of files to length required according to lmmss_bd
!
! Code Description:
! Language:           Fortran 90.
! Software Standards: "European Standards for Writing and
!                      Documenting Exchangeable Fortran 90 Code".
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
#ifdef NETCDF
USE netcdf
#endif

USE data_parameters , ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers, & ! KIND-type parameter for standard integer variables
    irealgrib, & ! KIND-type parameter for real variables in the grib library
    idouble,   & ! KIND-type parameter for double precision real variables
    isingle      ! KIND-type parameter for single precision real variables

USE data_int2lm_parallel,    ONLY :   &
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    nproc,           & ! total number of processors: nprocx * nprocy + nprocio
    icomm_world,     & ! communicator that belongs to igroup_world
    imp_logical        ! determines the correct LOGICAL   type used in the
                       ! model for MPI

USE data_int2lm_control,  ONLY: &
  dt,              & ! time step used in the COSMO
  itype_calendar,  & ! for specifying the calendar used
  yakdat1,         & ! actual date (ydate_ini+ntstep/dt)
  yakdat2,         & ! actual date (ydate_ini+ntstep/dt) in the form
  idbg_level,      & ! to control verbosity of output
  nstart             ! start time (in time steps)

USE data_int2lm_io,        ONLY : &
  njulianday,      & ! julian day in the year (returned by get_utc_date)
  ract_hour,       & ! actual hour of actual day (returned by get_utc_date)
  ydate_ini,       & ! start of the forecast yyyymmddhh (year,month,day,hour)
  lmmss_bd,        & ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                     ! if .FALSE. 10 digits date format (YYYYMMDDHH)
                     ! for ydate_bd  and input  files of INT2LM
  yin_cat,         & ! catalog-name of the input files
  ytunit_in          ! time unit for input data

USE data_grid_in,       ONLY: &
    lcm_hgt_coor,   & ! Input data has hybrid height coordinates
    lcm_pres_coor     ! Input data has pressure coordinates  !_br 14.03.12

USE environment,         ONLY :  &
    model_abort    ! one process stops the whole parallel program

USE utilities,            ONLY: &
    get_utc_date   ! determines the actual date

USE parallel_utilities, ONLY:  &
    distribute_values

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Interface Blocks

#ifdef NETCDF
INTERFACE nc_get_att
  MODULE PROCEDURE                                &
     nc_get_att_FourByteReals,                    &
     nc_get_att_Char
END INTERFACE
#endif

!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE setup_clm
  
#ifdef NETCDF

  CHARACTER (LEN=200) :: &
    att_value,  &    ! attribute value
    yerror           ! error message

  CHARACTER (LEN=250)        ::  &
    yzpath           ! full path and name of the input-file

  INTEGER (KIND=iintegers)   ::       &
    izbuflen,        & ! length of buffer distributing the NAMELIST
    izerror,         & ! error value
    izdebug,         & ! for verbosity of output
    izcharlen,       & ! length of character variable for date string 
                       ! written to the file name
    nzstat             ! for error-code on allocation

  REAL      (KIND=ireals)     ::             &
    racthour    ! actual hour of the day

  LOGICAL :: &
    lzexist     ! to check, whether different files exist

  LOGICAL  , ALLOCATABLE ::   &
    logbuf  (:)

  ! Allocate space for sending buffers
  izbuflen  = 1000
  ALLOCATE ( logbuf(izbuflen)   , STAT=nzstat )
  logbuf (:) = .FALSE.

  IF (my_cart_id == 0) THEN

    izdebug = idbg_level

    CALL get_utc_date (nstart, ydate_ini, dt, itype_calendar, yakdat1, yakdat2, &
                       njulianday, ract_hour)

    ! Since INT2LM Version 1.20, yakdat1 is given in the form yyyymmddhhmmss with
    ! with minutes and seconds. Whether the fields from the driving model are given
    ! in the same form or not, is specified by the form of ydate_bd:
    ! If ydate_bd is given with minutes and seconds, lmmss_bd is .TRUE., and then
    !    also the file names should correspond to that form
    ! If ydate_bd is given without minutes and seconds, lmmss_bd is .FALSE., and then
    !    the file names should be given as it was before.
    ! This is accounted for now by adapting the length of yakdat1, which is written
    !   to the file name

    IF (lmmss_bd) THEN
      izcharlen = 14
    ELSE
      izcharlen = 10
    ENDIF

    yzpath=''
    yzpath=TRIM(yin_cat)//'caf'//yakdat1(1:izcharlen)//'.nc'
    INQUIRE (FILE=TRIM(yzpath), EXIST=lzexist)
    IF (.NOT. lzexist) THEN      
      ! Create the file name again, trying cas as prefix
       yzpath=TRIM(yin_cat)//'cas'//yakdat1(1:izcharlen)//'.nc'
      INQUIRE (FILE=TRIM(yzpath), EXIST=lzexist)
      IF (.NOT. lzexist) THEN      
        ! Create the file name again, trying cffd as prefix
         yzpath=TRIM(yin_cat)//'cffd'//yakdat1(1:izcharlen)//'.nc'
        INQUIRE (FILE=TRIM(yzpath), EXIST=lzexist)
        IF (.NOT. lzexist) THEN      
          ! Create the file name again, trying cfsd as prefix
           yzpath=TRIM(yin_cat)//'cfsd'//yakdat1(1:izcharlen)//'.nc'
           INQUIRE (FILE=TRIM(yzpath), EXIST=lzexist)
           IF (.NOT. lzexist) THEN      
             izerror = 1
             yerror='None of the possible coarse grid output files caf, cas, cffd, cfsd exist'
             PRINT *, 'None of the coarse grid output files exist: '
             PRINT *, '  ',TRIM(yin_cat)//'caf'//yakdat1(1:izcharlen)//'.nc'
             PRINT *, '  ',TRIM(yin_cat)//'cas'//yakdat1(1:izcharlen)//'.nc'
             PRINT *, '  ',TRIM(yin_cat)//'cffd'//yakdat1(1:izcharlen)//'.nc'
             PRINT *, '  ',TRIM(yin_cat)//'cfsd'//yakdat1(1:izcharlen)//'.nc'
             CALL model_abort (my_cart_id, izerror, yerror, 'setup_clm')
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    izerror = nc_get_att (yzpath, 'level', 'standard_name', att_value, yerror)
    IF (izerror == NF90_NOERR .AND. TRIM(att_value) == "atmosphere_hybrid_height_coordinate") THEN
      lcm_hgt_coor = .TRUE.
    ELSE
      lcm_hgt_coor = .FALSE.
    ENDIF
!_br 14.03.12
     izerror = nc_get_att (yzpath, 'level', 'units', att_value, yerror)
    IF (izerror == NF90_NOERR .AND. (TRIM(att_value) == "hPa" .OR. TRIM(att_value) == "Pa")) THEN
      lcm_pres_coor = .TRUE.
    ELSE
      lcm_pres_coor = .FALSE.
    ENDIF
!_br 14.03.12 end
    
  ENDIF

  IF (nproc > 1) THEN

    IF (my_cart_id == 0) THEN
      logbuf  ( 1) = lcm_hgt_coor
      logbuf  ( 2) = lcm_pres_coor !_br 14.03.12
    ENDIF

    CALL distribute_values  (logbuf ,2, 0, imp_logical,   icomm_world, izerror)

    IF (my_cart_id /= 0) THEN
      lcm_hgt_coor  = logbuf  ( 1)
      lcm_pres_coor = logbuf  ( 2) !_br 14.03.12
    ENDIF

  ENDIF

  DEALLOCATE ( logbuf , STAT=nzstat )

#endif 
  
END SUBROUTINE setup_clm

!==============================================================================
!==============================================================================
#ifdef NETCDF

INTEGER FUNCTION nc_get_att_FourByteReals (ncfile, var_name, att_name, att_value, &
                                           yerrmsg)  

!------------------------------------------------------------------------------
! Description:
!   Retrieve a float attribute value from a netCDF file
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
  CHARACTER (LEN=*), INTENT(IN) :: &
    ncfile,             & ! name of NetCDF file
    var_name,           & ! name of variable to be read
    att_name              ! name of attribute to be read
  REAL,  INTENT(OUT) :: &
    att_value            ! array receiving the input data
  CHARACTER (LEN=200), OPTIONAL, INTENT(OUT) :: &
    yerrmsg              ! error message

!------------------------------------------------------------------------------
! local variables

  INTEGER :: &
    istatus,            & ! status variable
    ncid,               & ! ID of the NetCDF file
    varid                 ! ID of the variable "var_name"

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine nc_get_att_FourByteReals
!------------------------------------------------------------------------------

  att_value = -1.E20


! open the data file
  nc_get_att_FourByteReals = nf90_open(ncfile, NF90_NOWRITE, ncid)
  IF (nc_get_att_FourByteReals /= NF90_NOERR) THEN
!    PRINT *, 'Error in nc_get_att / nf90_open' 
!    PRINT *, TRIM(NF90_strerror(nc_get_att_FourByteReals))
    IF (PRESENT(yerrmsg)) THEN
      yerrmsg = 'Error in nc_get_att / nf90_open -- '// &
                TRIM(NF90_strerror(nc_get_att_FourByteReals))
    ENDIF
    istatus = nf90_close(ncid)
    RETURN
  ENDIF
  

  IF (TRIM(var_name) == 'GLOBAL' .OR. TRIM(var_name) == 'global') THEN
    varid = NF90_GLOBAL
  ELSE
!   set the ID for the quantity "label"
    nc_get_att_FourByteReals = nf90_inq_varid(ncid, TRIM(var_name), varid)
    IF (nc_get_att_FourByteReals /= NF90_NOERR) THEN
!      PRINT *, 'Error in nc_get_att / nf90_inq_varid' 
!      PRINT *, TRIM(NF90_strerror(nc_get_att_FourByteReals))
      IF (PRESENT(yerrmsg)) THEN
        yerrmsg = 'Error in nc_get_att / nf90_inq_varid -- '// &
                  TRIM(NF90_strerror(nc_get_att_FourByteReals))
      ENDIF
      istatus = nf90_close(ncid)
      RETURN
    ENDIF
  ENDIF

! get the attribute value
  nc_get_att_FourByteReals = nf90_get_att (ncid, varid, att_name, att_value)
  IF (nc_get_att_FourByteReals /= NF90_NOERR) THEN
!    PRINT *, 'Error in nc_get_att / nf90_get_att'
!    PRINT *, TRIM(NF90_strerror(nc_get_att_FourByteReals))
    IF (PRESENT(yerrmsg)) THEN
      yerrmsg = 'Error in nc_get_att / nf90_get_att -- ' // &
                TRIM(NF90_strerror(nc_get_att_FourByteReals))
    ENDIF
    istatus = nf90_close(ncid)
    RETURN
  ENDIF
  
! close the data file
  nc_get_att_FourByteReals = nf90_close(ncid)
  IF (nc_get_att_FourByteReals /= NF90_NOERR) THEN
!    PRINT *, 'Error in nc_get_att / nf90_close' 
!    PRINT *, TRIM(NF90_strerror(nc_get_att_FourByteReals))
    IF (PRESENT(yerrmsg)) THEN
      yerrmsg = 'Error in nc_get_att / nf90_close -- '// &
                TRIM(NF90_strerror(nc_get_att_FourByteReals))
    ENDIF
    RETURN
  ENDIF


END FUNCTION nc_get_att_FourByteReals

!------------------------------------------------------------------------------

!+ Function for getting an attribute

INTEGER FUNCTION nc_get_att_Char (ncfile, var_name, att_name, att_value, yerrmsg)  

!------------------------------------------------------------------------------
! Description:
!   Retrieve a character attribute value from a netCDF file
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
  CHARACTER (LEN=*), INTENT(IN) :: &
    ncfile,             & ! name of NetCDF file
    var_name,           & ! name of variable to be read
    att_name              ! name of attribute to be read
  CHARACTER (LEN=*),  INTENT(OUT) :: &
    att_value            ! array receiving the input data
  CHARACTER (LEN=200), OPTIONAL, INTENT(OUT) :: &
    yerrmsg              ! error message

!------------------------------------------------------------------------------
! local variables

  INTEGER :: &
    istatus,            & ! status variable
    ncid,               & ! ID of the NetCDF file
    varid                 ! ID of the variable "var_name"

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine nc_get_att_Char
!------------------------------------------------------------------------------

  att_value = 'NONE'

! open the data file
  nc_get_att_Char = nf90_open(ncfile, NF90_NOWRITE, ncid)
  IF (nc_get_att_Char /= NF90_NOERR) THEN
!    PRINT *, 'Error in nc_get_att / nf90_open' 
!    PRINT *, TRIM(NF90_strerror(nc_get_att_Char))
    IF (PRESENT(yerrmsg)) THEN
      yerrmsg = 'Error in nc_get_att / nf90_open -- ' // &
                 TRIM(NF90_strerror(nc_get_att_Char))
    ENDIF
    istatus = nf90_close(ncid)
    RETURN
  ENDIF

  IF (TRIM(var_name) == 'GLOBAL' .OR. TRIM(var_name) == 'global') THEN
    varid = NF90_GLOBAL
  ELSE
!   set the ID for the quantity "label"
    nc_get_att_Char = nf90_inq_varid(ncid, TRIM(var_name), varid)
    IF (nc_get_att_Char /= NF90_NOERR) THEN
!      PRINT *, 'Error in nc_get_att / nf90_inq_varid' 
!      PRINT *, TRIM(NF90_strerror(nc_get_att_Char))
      IF (PRESENT(yerrmsg)) THEN
        yerrmsg = 'Error in nc_get_att / nf90_inq_varid -- '// &
                  TRIM(NF90_strerror(nc_get_att_Char))
      ENDIF
      istatus = nf90_close(ncid)
      RETURN
    ENDIF
  ENDIF

! read the attribute value
  nc_get_att_Char = nf90_get_att (ncid, varid, att_name, att_value)
  IF (nc_get_att_Char /= NF90_NOERR) THEN
!    PRINT *, 'Error in nc_get_att / nf90_get_att' 
!    PRINT *, TRIM(NF90_strerror(nc_get_att_Char))
    IF (PRESENT(yerrmsg)) THEN
      yerrmsg = 'Error in nc_get_att / nf90_get_att -- '// &
                TRIM(NF90_strerror(nc_get_att_Char))
    ENDIF
    istatus = nf90_close(ncid)
    RETURN
  ENDIF
  
! close the data file
  nc_get_att_Char = nf90_close(ncid)
  IF (nc_get_att_Char /= NF90_NOERR) THEN
!    PRINT *, 'Error in nc_get_att / nf90_close' 
!    PRINT *, TRIM(NF90_strerror(nc_get_att_Char))
    IF (PRESENT(yerrmsg)) THEN
      yerrmsg = 'Error in nc_get_att / nf90_close -- '// &
                TRIM(NF90_strerror(nc_get_att_Char))
    ENDIF
    RETURN
  ENDIF


END FUNCTION nc_get_att_Char

#endif
!==============================================================================
!==============================================================================

END MODULE clm_utilities
