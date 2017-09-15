! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!+ Dummy interfaces for DWD database routines
!------------------------------------------------------------------------------
!
! Description:
!   This file provides dummy interfaces for the calls to the DWD
!   database system, which is not available on other than DWD machines.
!   The meaning of this file is to avoid complaints about "unsatisfied
!   external references".
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! !VERSION!  !DATE!     Ulrich Schaettler
!  Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

SUBROUTINE csodb_init( DB_init_order, DB_OUT_SOCKETS, DB_IN_SOCKETS, &
                       DB_retry, cso_dbg_level, error)

    INTEGER   DB_OUT_SOCKETS, DB_IN_SOCKETS,DB_retry, cso_dbg_level, error
    CHARACTER (LEN=*) DB_init_order

END SUBROUTINE csodb_init

!==============================================================================

SUBROUTINE csodb_arr_wait_out(ihits, err_close)

    INTEGER   ihits, err_close

PRINT *, '*** The routine csodb_arr_wait_out used is only a dummy ***'
PRINT *, '***          and should have never been called!!        ***'
err_close = 9999

END SUBROUTINE csodb_arr_wait_out

!==============================================================================

SUBROUTINE csodb_arrout_req( DB_order, request,                   &
                             DB_hits,error, buff_len, ilenc, data)

    INTEGER   request, DB_hits,error, ilenc, data, buff_len (*)
    CHARACTER (LEN=*) DB_order

PRINT *, '*** The routine csodb_arrout_req   used is only a dummy ***'
PRINT *, '***          and should have never been called!!        ***'
error = 9999

END SUBROUTINE csodb_arrout_req

!==============================================================================

SUBROUTINE csodb_arrout_end( request, DB_hits, error )

    INTEGER request, DB_hits, error

PRINT *, '*** The routine csodb_arrout_end   used is only a dummy ***'
PRINT *, '***          and should have never been called!!        ***'
error = 9999

END SUBROUTINE csodb_arrout_end

!==============================================================================

SUBROUTINE csodb_arrin_req(DB_order, DB_request, DB_err )
  
    INTEGER DB_err, DB_request(*)
    CHARACTER (LEN=*) DB_order

PRINT *, '*** The routine csodb_arrin_req    used is only a dummy ***'
PRINT *, '***          and should have never been called!!        ***'
DB_err = 9999

END SUBROUTINE csodb_arrin_req

!==============================================================================

SUBROUTINE csodb_arrin_end(DB_request, DB_hits, DB_err, ilfd, ilen, data)
  
    INTEGER DB_hits, DB_err, ilfd, ilen, data, DB_request(*)

PRINT *, '*** The routine csodb_arrin_end    used is only a dummy ***'
PRINT *, '***          and should have never been called!!        ***'
DB_err = 9999

END SUBROUTINE csodb_arrin_end

!==============================================================================

SUBROUTINE csodb_close(DB_err)

    INTEGER DB_err

PRINT *, '*** The routine csodb_close used is only a dummy ***'
PRINT *, '***     and should have never been called!!      ***'
DB_err = 9999

END SUBROUTINE csodb_close

!==============================================================================

SUBROUTINE DB_stat_out                                               &
           (input_time, IMBytes, inp_calls, inp_succ,                &
            output_time, OMBytes, outp_calls, outp_succ,             &
            env_time, env_calls, env_succ)

   REAL  input_time, IMBytes, output_time, OMBytes, env_time
   INTEGER  inp_calls, inp_succ, outp_calls, outp_succ, env_calls, env_succ

PRINT *, '*** The routine DB_stat_out used is only a dummy ***'
PRINT *, '***     and should have never been called!!      ***'

END SUBROUTINE DB_stat_out

!==============================================================================

SUBROUTINE grib_save_init (yname, myid, isize, ierror)

    CHARACTER (LEN=*)    yname
    INTEGER              myid, isize, ierror

PRINT *, '*** The routine grib_save_init used is only a dummy ***'
PRINT *, '***       and should have never been called!!       ***'
ierror = 9999

END SUBROUTINE grib_save_init

!==============================================================================

SUBROUTINE grib_save_open (yname, myid, ierror)

    CHARACTER (LEN=*)    yname
    INTEGER              myid, ierror

PRINT *, '*** The routine grib_save_open used is only a dummy ***'
PRINT *, '***       and should have never been called!!       ***'
ierror = 9999

END SUBROUTINE grib_save_open

!==============================================================================

SUBROUTINE grib_save_delete (myid, ierror)

    INTEGER              myid, ierror

PRINT *, '*** The routine grib_save_delete used is only a dummy ***'
PRINT *, '***        and should have never been called!!        ***'
ierror = 9999

END SUBROUTINE grib_save_delete

!==============================================================================

SUBROUTINE grib_save (yname, myid, data, ilen, ierror)

    CHARACTER (LEN=*)    yname
    INTEGER              myid, ilen, ierror, data(*)

PRINT *, '*** The routine grib_save used is only a dummy ***'
PRINT *, '***    and should have never been called!!     ***'
ierror = 9999

END SUBROUTINE grib_save

!==============================================================================

SUBROUTINE grib_save_close (myid, ierror)

    INTEGER              myid, ierror

PRINT *, '*** The routine grib_save_close used is only a dummy ***'
PRINT *, '***       and should have never been called!!        ***'
ierror = 9999

END SUBROUTINE grib_save_close

!==============================================================================

SUBROUTINE grib_save_filename (yname, myid, index, ierror)

    CHARACTER (LEN=*)    yname
    INTEGER              myid, index, ierror

PRINT *, '*** The routine grib_save_filename used is only a dummy ***'
PRINT *, '***         and should have never been called!!         ***'
ierror = 9999

END SUBROUTINE grib_save_filename

!==============================================================================

SUBROUTINE grib_retrieve (yname, myid, ibuffer, ibuflen, igriblen, ierror)

    CHARACTER (LEN=*)    yname
    INTEGER              myid, ibuflen, igriblen, ierror, ibuffer(*)

PRINT *, '*** The routine grib_retrieve used is only a dummy ***'
PRINT *, '***       and should have never been called!!      ***'
ierror = 9999

END SUBROUTINE grib_retrieve

!==============================================================================

SUBROUTINE io_stat (isend, rsend)

    INTEGER       isend
    REAL          rsend

PRINT *, '*** The routine io_stat used is only a dummy ***'
PRINT *, '***    and should have never been called!!   ***'

END SUBROUTINE io_stat

!==============================================================================
