!+ Data module for organizational variables for the parallel GME2LM
!==============================================================================

MODULE data_int2lm_parallel

!==============================================================================
!
! Description:
!  This module contains all variables that are necessary for the organization
!  of the parallel program (number of processors, virtual topology, neighbors,
!  etc.). In a sequential environment some of these variables are also
!  necessary and are set accordingly (e.g. nproc = 1 and my_*_id = 0).
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
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_22        2013/07/11 Ulrich Schaettler
!  Introduced new MPI data type for special grib_api integers kindOfSize
!   (which could be a 4 or a 8-byte integer)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters , ONLY :   &
           ireals,    & ! KIND-type parameters for real variables
           iintegers    ! kind-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! 1. Information about the processor grid
! ---------------------------------------

  LOGICAL                           ::           &
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    lasync_io          ! if .TRUE.: the model runs with extra PEs for
                       ! asynchronous IO

  INTEGER   (KIND=iintegers)       ::           &
    nprocx,          & ! number of processors in x-direction for LM
    nprocy,          & ! number of processors in y-direction for LM
    nprocio,         & ! number of extra processors for doing asynchronous IO
    nproc,           & ! total number of processors: nprocx * nprocy + nprocio
    num_compute,     & ! number of compute PEs
    num_io,          & ! number of IO PEs
    nproc1,          & ! number of processors in first direction for GME
                       ! = nprocy (!)
    nproc2,          & ! number of processors in second direction for GME
                       ! = nprocx (!)

    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ncomm_type,      & ! type of communication for boundary exchange

    my_world_id,     & ! rank of this subdomain in the global communicator
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_pos(2),  & ! position of this subdomain in the cartesian grid 
                       ! in x- and y-direction
    my_num1,         & ! position of this subdomain in the cartesian grid 
    my_num2,         & ! in first and second direction 
                       ! (my_num1 = my_cart_pos(2) and my_num2 = my_cart_pos(1))
    my_cart_neigh(4)   ! neighbors of this subdomain in the cartesian grid

  INTEGER   (KIND=iintegers), ALLOCATABLE       ::           &
    isubpos(:,:),    & ! positions of the subdomains in the total domain. Given
                       ! are the i- and the j-indices of the lower left and the
                       ! upper right grid point in the order 
                       !                  i_ll, j_ll, i_ur, j_ur.
                       ! Only the interior of the domains are considered, not 
                       ! the boundary lines.
    isubpos_coarse(:,:)! the same for the coarse grid

! 2. Further information for MPI
! ------------------------------

  INTEGER   (KIND=iintegers)       ::           &
    igroup_world,    & ! group that belongs to MPI_COMM_WORLD, i.e. all 
                       ! processors
    icomm_world,     & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    icomm_compute,   & ! communicator for the group of compute PEs
    igroup_cart,     & ! group of the compute PEs
    icomm_cart,      & ! communicator for the virtual cartesian topology
    icomm_row,       & ! communicator for a east-west row of processors
    iexch_req(4),    & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_grib,        & ! determines the REAL type used for the GRIB library
    imp_integers,    & ! determines the correct INTEGER type used in the
                       ! model for MPI
    imp_integ_ga,    & ! determines the correct INTEGER type used for grib_api
    imp_byte,        & ! determines the correct BYTE type used in the model
                       ! for MPI
    imp_character,   & ! determines the correct CHARACTER type used in the
                       ! model for MPI
    imp_logical        ! determines the correct LOGICAL   type used in the
                       ! model for MPI

  LOGICAL                                       &
    lcompute_pe,     & ! indicates whether this is a compute PE or not
    lreorder           ! during the creation of the virtual topology the
                       ! ranking of the processors may be reordered

! 3. Static Send buffers
! ----------------------

  REAL (KIND=ireals), ALLOCATABLE  ::           &
    sendbuf(:,:)       ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving

  INTEGER (KIND=iintegers)         ::           &
    isendbuflen        ! length of one column of sendbuf 

!==============================================================================

END MODULE data_int2lm_parallel
