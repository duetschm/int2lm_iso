!+ Module for the decomposition of the different grids on parallel computers
!==============================================================================

MODULE src_decomposition

!==============================================================================
!
! Description:
!   This module provides routines for the decomposition of the LM output grid
!   and for the different possible input grids. In Debug mode, information
!   on the decomposition is written to file YUDEBUG.
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
!  Eliminated unused variables
!  Replaced nproc by num_compute, because only compute PEs are involved
!  Introduced nlevskip to  skip levels at the top of ECMWF model (Davide Cesari)
!  Adapted coarse grid check to domains that are located around the date line
!    (by Friedrich Theunert)
! V1_6         2007/09/07 Burkhardt Rockel, Ulrich Schaettler
!  Added lcm2lm and SR decompose_cm for processing data from climate models
! V1_8         2008/05/29 Hans-Juergen Panitz
!  Bug fix in check_decomposition for gathering information from the nodes
! V1_9         2009/09/03 Ulrich Schaettler, Burkhardt Rockel
!  Corrected check of coarse domain in decompose_cm (minlon belongs to west
!  and maxlon to east)
!  Scale endlon from the total domain to 0...360 in SR decompose_coarse_grid
!  Added an error message for field boundaries
! V1_11        2010/04/23 Michael Gertz
!  Adaptions to SVN
! V1_17        2011/03/11 Ulrich Schaettler
!  Bug fix for isolated points-reproducibility for coarse grids that have a
!  different rotation than the fine COSMO grid
! V1_18        2011/03/11 Ulrich Schaettler
!  ANSI fix in MIN-,MAX-arguments
! V1_19        2012/06/06 Ulrich Schaettler, Daniel Luethi
!  Added initialization of ke1in also for GME data in decompose_gme
!  Adaptations of decompose_cm when domain crosses the date-line
!  Make sure that the same part of coarse grid domain is considered for 
!    sequential and for parallel runs for checking the isolated grid points
!    (Adaptation in Section 3 of decompose_cm)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Modules used:
USE data_parameters , ONLY :   &
    ireals,     & ! kind-type parameter for "normal" integer variables
    iintegers     ! KIND-type parameters for real variables

!------------------------------------------------------------------------------

USE data_int2lm_control,     ONLY :   &
    idbg_level,   & ! to control verbosity of output
    lprintdeb_all,& ! whether all PEs print debug output
    lgme2lm,      & ! if .TRUE., gme->lm, if .FALSE. gme->hm
    lec2lm,       & ! if .TRUE., ec ->lm
    llm2lm,       & ! if .TRUE., lm ->lm
    lcm2lm          ! if .TRUE., cm ->lm

!------------------------------------------------------------------------------

USE data_grid_lm,     ONLY :   &
    pollat,      & ! latitude of the rotated north pole (in degrees, N>0)
    pollon,      & ! longitude of the rotated north pole (in degrees, E>0)
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
    ke1lm,       & !
    ie2lm,       & !
    je2lm,       & !
    istartpar,   & ! start index for computations in the parallel program
    iendpar,     & ! end index for computations in the parallel program
    jstartpar,   & ! start index for computations in the parallel program
    jendpar        ! end index for computations in the parallel program

!------------------------------------------------------------------------------

USE data_grid_in,     ONLY :   &
    ids,         & ! start index of diamonds (ids = 1)
    ide,         & ! end index of diamonds (ide = 10)
    ni_gme,      & ! resolution of GME
    i3e_gme,     & ! number of vertical layers in GME
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
    ie_in_tot,   & ! ie for input grid, total domain
    je_in_tot,   & ! je for input grid, total domain
    ke_in_tot,   & ! ke for input grid
    ie_in_max,   & ! Max. of ie_in (local) on all processors
    je_in_max,   & ! Max. of je_in (local) on all processors
    ie_in,       & ! ie for input grid, local domain
    je_in,       & ! je for input grid, local domain
    ke_in,       & ! ke for input grid
    nlevskip,    & ! number of levels to skip at top of input model
    ke1in          ! ke_in+1

USE data_grid_in,     ONLY :   &
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

!------------------------------------------------------------------------------

USE data_fields_lm,     ONLY:                              &
  latlm_m   ,      & ! latitudes of the LM grid points
  lonlm_m   ,      & ! longitudes of the LM grid points
  latlm_u   ,      & ! latitudes of the LM u grid points
  lonlm_u   ,      & ! longitudes of the LM u grid points
  latlm_v   ,      & ! latitudes of the LM v grid points
  lonlm_v            ! longitudes of the LM v grid points

!------------------------------------------------------------------------------

USE data_fields_in,     ONLY:                              &
    lat_coarse_m,& ! latitudes of the LM grid points
    lon_coarse_m,& ! longitudes of the LM grid points
    lat_coarse_u,& ! latitudes of the LM u grid points
    lon_coarse_u,& ! longitudes of the LM u grid points
    lat_coarse_v,& ! latitudes of the LM v grid points
    lon_coarse_v   ! longitudes of the LM v grid points

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
    nprocx,          & ! number of processors in x-direction
    nprocy,          & ! number of processors in y-direction
    nprocio,         & ! number of extra processors for asynchronous IO
    num_compute,     & ! number of compute PEs
    num_io,          & ! number of IO PEs
    nproc1,          & ! number of processors in first direction for GME
                       ! = nprocy (!)
    nproc2,          & ! number of processors in second direction for GME
                       ! = nprocx (!)

    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains

    my_world_id,     & ! rank of this subdomain in the global communicator
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_pos,     & ! position of this subdomain in the cartesian grid 
                       ! in x- and y-direction
    my_num1,         & ! position of this subdomain in the cartesian grid
    my_num2,         & ! in first and second direction
                       ! (my_num1 = my_cart_pos(2) and my_num2 = my_cart_pos(1))
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    isubpos,         & ! positions of the subdomains in the total domain.
    isubpos_coarse,  & ! positions of the subdomains in the total domain.
    icomm_compute,   & ! communicator for the group of compute PEs
    igroup_cart,     & ! group of the compute PEs
    icomm_cart,      & ! communicator for the virtual cartesian topology
    icomm_row,       & ! communicator for a east-west row of processors
    igroup_world,    & ! group that belongs to MPI_COMM_WORLD, i.e. all 
                       ! processors
    icomm_world,     & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_integers,    & ! determines the correct INTEGER type used in the
                       ! model for MPI
    imp_byte,        & ! determines the correct BYTE type used in the model
                       ! for MPI
    imp_character,   & ! determines the correct CHARACTER type used in the
                       ! model for MPI
    imp_logical,     & ! determines the correct LOGICAL   type used in the
                       ! model for MPI
    lcompute_pe,     & ! indicates whether this is a compute PE or not
    lreorder,        & ! during the creation of the virtual topology the
                       ! ranking of the processors may be reordered
    sendbuf,         & ! sending buffer for boundary exchange
    isendbuflen        ! length of one column of sendbuf.

!------------------------------------------------------------------------------

USE data_int2lm_io,     ONLY :   &
    yinext_cat,   & ! catalog-name of the file with external GME parameters
    yinext_lfn      ! name of the file with external GME parameters

!------------------------------------------------------------------------------

USE parallel_utilities, ONLY:  global_values, gather_values
USE utilities,          ONLY:  phi2phirot, rla2rlarot

USE src_read_ext,       ONLY:  read_nc_axis

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Public and Private Subroutines

PUBLIC  decompose_lm, decompose_gme, decompose_coarse_grid, decompose_cm,  &
        check_decomposition

!==============================================================================

CONTAINS

!==============================================================================
!+ Subroutine that decomposes the LM domain for distributed memory computers
!------------------------------------------------------------------------------

SUBROUTINE decompose_lm

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine computes the decomposition of the total LM domain. For 
!   dealing with the cartesian grid, it uses the following organizational
!   variables which are determined in init_procgrid.
!    - my_cart_id:       rank and id of this processor in the virtual topology
!    - my_cart_pos(2):   position in the cartesian processor grid in 
!                        x- and y-direction
!    - my_cart_neigh(4): neighbours of this processor in the order west, north,
!                        east, south
!    - icomm_cart:       MPI-communicator for the cartesian grid
!
!   With the above information, the decomposition of the total domain is 
!   computed. For every subdomain the indices of the lower left and the
!   upper right corner in the total domain are computed and stored in the
!   variable isubpos (num_compute,4) in the order (i_ll,j_ll,i_ur,j_ur).
!   Note that only the interior of the subdomains are considered and
!   boundary lines are neglected. That means
!
!         total domain                    subdomain
!
!      (i_ll,j_ll) corresponds to  (1+nboundlines , 1+nboundlines)
!      (i_ur,j_ur) corresponds to  (ie-nboundlines,je-nboundlines)
!
!   For the diagnostic computations, the selected grid points for output 
!   of profiles are distributed to the appropriate nodes.
!
!   For the sequential program the variables are set accordingly.
!
! Method:
!
!------------------------------------------------------------------------------

! Local variables

  INTEGER (KIND=iintegers)   ::       &
    nz1d, nzsubi, nzsubj, nzcompi, nzcompj, nzix1, nzix2, nzjy1, nzjy2,  &
    nzix2right, nzix2left, nzjy2lower, nzjy2upper, ix, iy,               &
    iup, ilo, jup, jlo, n, izerror

  INTEGER (KIND=iintegers)   ::       &
    intvec(2)       ! for collecting values from all nodes

  CHARACTER (LEN=80)         ::       &
    yzerrmsg

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE decompose_lm
!------------------------------------------------------------------------------

  izerror  = 0
  yzerrmsg = '   '

!------------------------------------------------------------------------------
!- Section 1: Compute the domain decomposition
!------------------------------------------------------------------------------

  IF (num_compute > 1) THEN

  !----------------------------------------------------------------------------
  !- Section 1.1: Sizes and distribution of the subdomains
  !----------------------------------------------------------------------------

    ! Number of grid points that have to be distributed in each direction.
    ! The first nboundlines boundary lines of the total domain are not
    ! considered.
    nzcompi = ie2lm_tot - 2*nboundlines   ! = ielm_tot
    nzcompj = je2lm_tot - 2*nboundlines   ! = jelm_tot

    ! Number of grid points a subdomain gets at least: nzsubi, nzsubj
    nzsubi  = nzcompi / nprocx
    nzsubj  = nzcompj / nprocy

    ! Determine how many subdomains will get nzsubi (nzix1) and how many will
    ! get nzsubi+1 (nzix2) grid points: nzix1, nzix2
    nzix2   = nzcompi - nprocx * nzsubi
    nzix1   = nprocx - nzix2

    ! Determine how many subdomains will get nzsubj (nzjy1) and how many will
    ! get nzsubj+1 (nzjy2) grid points
    nzjy2   = nzcompj - nprocy * nzsubj
    nzjy1   = nprocy - nzjy2

    ! Determine the distribution of the subdomains with different sizes.
    ! The ones with more grid points are placed to the interior, the ones
    ! with less grid points to the boundary of the processor grid.
    nzix2left  = nzix1 / 2
    nzix2right = nzix1 - nzix2left
    nzjy2lower = nzjy1 / 2
    nzjy2upper = nzjy1 - nzjy2lower
   
  !----------------------------------------------------------------------------
  !- Section 1.2: Position of the subdomains in the total domain
  !----------------------------------------------------------------------------

    DO ix = 0,nprocx-1
      DO iy = 0,nprocy-1
        ! 1D numbering of the processors: rank
        nz1d = ix * nprocy + iy

        IF ( (0 <= iy) .AND. (iy <= nzjy2lower-1) ) THEN
          isubpos (nz1d,2) =  iy    *  nzsubj + nboundlines + 1
          isubpos (nz1d,4) = (iy+1) *  nzsubj + nboundlines
        ELSEIF ( (nzjy2lower <= iy) .AND. (iy <= nzjy2lower+nzjy2-1) ) THEN
          isubpos (nz1d,2) =  iy    * (nzsubj+1) - nzjy2lower + nboundlines + 1
          isubpos (nz1d,4) = (iy+1) * (nzsubj+1) - nzjy2lower + nboundlines
        ELSEIF ( (nzjy2lower+nzjy2 <= iy) .AND. (iy <= nprocy-1) ) THEN
          isubpos (nz1d,2) =  iy    *  nzsubj + nzjy2 + nboundlines + 1
          isubpos (nz1d,4) = (iy+1) *  nzsubj + nzjy2 + nboundlines
        ENDIF

        IF ( (0 <= ix) .AND. (ix <= nzix2left-1) ) THEN
          isubpos (nz1d,1) =  ix    *  nzsubi + nboundlines + 1
          isubpos (nz1d,3) = (ix+1) *  nzsubi + nboundlines
        ELSEIF ( (nzix2left <= ix) .AND. (ix <= nzix2left+nzix2-1) ) THEN
          isubpos (nz1d,1) =  ix    * (nzsubi+1) - nzix2left + nboundlines + 1
          isubpos (nz1d,3) = (ix+1) * (nzsubi+1) - nzix2left + nboundlines
        ELSEIF ( (nzix2left+nzix2 <= ix) .AND. (ix <= nprocx-1) ) THEN
          isubpos (nz1d,1) =  ix    *  nzsubi + nzix2 + nboundlines + 1
          isubpos (nz1d,3) = (ix+1) *  nzsubi + nzix2 + nboundlines
        ENDIF

      ENDDO
    ENDDO

  !----------------------------------------------------------------------------
  !- Section 2.3: Compute lmgrid variables for this subdomain
  !----------------------------------------------------------------------------

    ! Attention: in contrast to the LM, ielm and jelm are the sizes of a 
    ! subdomain without (!) the boundary lines. That means that the LM fields
    ! have to be allocated with ie2lm and je2lm, resp. 
    ielm = isubpos (my_cart_id,3) - isubpos (my_cart_id,1) + 1
    jelm = isubpos (my_cart_id,4) - isubpos (my_cart_id,2) + 1
    kelm = kelm_tot

    ie2lm = ielm + 2*nboundlines
    je2lm = jelm + 2*nboundlines
    ke1lm = kelm_tot + 1

    ! Calculate ie2lm_max and je2lm_max
    intvec(1) = ie2lm
    intvec(2) = je2lm
    CALL global_values (intvec, 2, 'MAX', imp_integers, icomm_cart, -1,    &
                        yzerrmsg, izerror)
    ie2lm_max = intvec(1)
    je2lm_max = intvec(2)

    ! These values give the start- and end-latitude or longitude, resp.
    ! of a subdomain with boundary lines. That means, gridpoint (1,1)
    ! of the subdomain has (startlon,startlat) as longitude or latitude, resp.
    startlon = startlon_tot + (isubpos(my_cart_id,1)-2*nboundlines-1) * dlon
    startlat = startlat_tot + (isubpos(my_cart_id,2)-2*nboundlines-1) * dlat
    endlon   = startlon_tot + (isubpos(my_cart_id,3)              -1) * dlon
    endlat   = startlat_tot + (isubpos(my_cart_id,4)              -1) * dlat

    ! Determine start- and end-indices for computations
    IF (my_cart_neigh(1) == -1) THEN
      istartpar = 1
    ELSE
      istartpar = 1 + nboundlines
    ENDIF

    IF (my_cart_neigh(4) == -1) THEN
      jstartpar = 1
    ELSE
      jstartpar = 1 + nboundlines
    ENDIF

    IF (my_cart_neigh(3) == -1) THEN
      iendpar   = ie2lm
    ELSE
      iendpar   = ie2lm - nboundlines
    ENDIF

    IF (my_cart_neigh(2) == -1) THEN
      jendpar   = je2lm
    ELSE
      jendpar   = je2lm - nboundlines
    ENDIF

  ELSE   ! num_compute==1, ie. sequential version
    ! set the variables accordingly for the sequential program
    isubpos(0,1) = 1 + nboundlines
    isubpos(0,2) = 1 + nboundlines
    isubpos(0,3) = ie2lm_tot - nboundlines
    isubpos(0,4) = je2lm_tot - nboundlines
    ielm         = ielm_tot
    jelm         = jelm_tot
    kelm         = kelm_tot
    ie2lm        = ielm + 2*nboundlines
    je2lm        = jelm + 2*nboundlines
    ke1lm        = kelm_tot + 1
    ie2lm_max    = ie2lm
    je2lm_max    = je2lm
    startlon     = startlon_tot - dlon
    startlat     = startlat_tot - dlat
    endlon       = endlon_tot   + dlon
    endlat       = endlat_tot   + dlat
    istartpar    = 1
    jstartpar    = 1
    iendpar      = ie2lm
    jendpar      = je2lm
  ENDIF

!------------------------------------------------------------------------------
!- Section 3: Decomposition of grid points for profiles
!------------------------------------------------------------------------------

  ! Grid points can be chosen on the boundary lines of the total domain.
  ! Therefore determine ilo, iup, jlo, jup to look in which subdomain a 
  ! certain grid point is located.

  IF (num_compute > 1) THEN

    IF (my_cart_neigh(1) == -1) THEN
      ilo = 1
    ELSE
      ilo = isubpos(my_cart_id,1)
    ENDIF

    IF (my_cart_neigh(2) == -1) THEN
      jup = jelm_tot
    ELSE
      jup = isubpos(my_cart_id,4)
    ENDIF

    IF (my_cart_neigh(3) == -1) THEN
      iup = ielm_tot
    ELSE
      iup = isubpos(my_cart_id,3)
    ENDIF

    IF (my_cart_neigh(4) == -1) THEN
      jlo = 1
    ELSE
      jlo = isubpos(my_cart_id,2)
    ENDIF

    ngp = 0
    
    DO n=1,ngp_tot
      IF ( (ilo <= igp_tot(n)) .AND. (igp_tot(n) <= iup) ) THEN
        IF ( (jlo <= jgp_tot(n)) .AND. (jgp_tot(n) <= jup) ) THEN

          ! grid point n is located in this subdomain (my_cart_id)
          ngp = ngp + 1
          igp(ngp) = nboundlines +                                         &
                        (igp_tot(n) - isubpos(my_cart_id,1) + 1)
          jgp(ngp) = nboundlines +                                         &
                        (jgp_tot(n) - isubpos(my_cart_id,2) + 1)
        ENDIF
      ENDIF
    ENDDO

  ELSE    ! num_compute==1, i.e. sequential version

    ! set the *_tot variables to the local variables
    ngp          = ngp_tot
    igp(:)       = igp_tot(:)
    jgp(:)       = jgp_tot(:)

  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE decompose_lm

!==============================================================================
!+ Subroutine that decomposes the GME domain for distributed memory computers
!------------------------------------------------------------------------------

SUBROUTINE decompose_gme

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine decomposes the GME domain on a 2-dimensional processor
!   grid. The following variables are computed:
!    - ilim1, ilim2:
!        decomposition limits in x- and y-direction
!    - ig1s, ig1e, ig2s, ig2e:
!        local start and end-indices of array-dimension in x- and y-direction
!    - igg1s, igg1e, igg2s, igg2e:
!        global start and end-indices of array-dimension in x- and y-direction
!    - ig1sm1, ig1sm2, ig1ep1, ig1ep2: 
!    - ig2sm1, ig2sm2, ig2ep1, ig2ep2:
!        further help variables
!    - nbpe, nbdw, nbpw, nbae:
!        neighbors for GME subdomains
!
! Method:
!
!------------------------------------------------------------------------------

! Local variables
  INTEGER (KIND=iintegers)     ::   nrows, nrem, i, izerror

  CHARACTER (LEN=80)           ::   yzerrmsg

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

  izerror  = 0
  yzerrmsg = '  '

!------------------------------------------------------------------------------
!- Begin SUBROUTINE decompose_gme
!------------------------------------------------------------------------------

  ! Set number of PEs in the two GME-directions
  nproc1 = nprocy
  nproc2 = nprocx

  IF (num_compute > 1) THEN

    ! Set decomposition limits in first direction (y-direction in LM)
    nrows = (ni_gme+1)/nproc1
    nrem  = MOD(ni_gme+1,nproc1)
    ilim1(0) = 0
    DO i=1,nrem
      ilim1(i) = ilim1(i-1) + (nrows+1)
    ENDDO
    DO i=nrem+1,nproc1
      ilim1(i) = ilim1(i-1) + nrows
    ENDDO

    ! Set decomposition limits in second direction (x-direction in LM)
    nrows = (ni_gme+1)/nproc2
    nrem  = MOD(ni_gme+1,nproc2)
    ilim2(0) = 1
    DO i=1,nrem
      ilim2(i) = ilim2(i-1) + (nrows+1)
    ENDDO
    DO i=nrem+1,nproc2
      ilim2(i) = ilim2(i-1) + nrows
    ENDDO

    ! Set ig1s ... ig2e
    my_num1 = MOD (my_cart_id, nproc1)   ! = my_cart_pos(2)
    my_num2 = my_cart_id/nproc1          ! = my_cart_pos(1)

    ig1s = ilim1(my_num1  )
    ig1e = ilim1(my_num1+1) - 1

    ig2s = ilim2(my_num2  )
    ig2e = ilim2(my_num2+1) - 1

    ! Set neighbors
    IF (my_num1 == 0) THEN
      nbpe = -1
    ELSE
      nbpe = my_cart_id-1
    ENDIF

    IF (my_num1 == nproc1-1) THEN
      nbaw = -1
    ELSE
      nbaw = my_cart_id+1
    ENDIF

    IF (my_num2 == 0) THEN
      nbpw = -1
    ELSE
      nbpw = my_cart_id-nproc1
    ENDIF

    IF (my_num2 == nproc2-1) THEN
      nbae = -1
    ELSE
      nbae = my_cart_id+nproc1
    ENDIF

  ELSE    ! sequential version: num_compute = 1

    ! Set decomposition limits
    ilim1(0) = 0
    ilim1(1) = ni_gme+1
    ilim2(0) = 1
    ilim2(1) = ni_gme+2

    nrows = (ni_gme+1)/num_compute
    nrem  = MOD(ni_gme+1,num_compute)

    ! Set ig1s ... ig2e
    ig1s = 0
    ig1e = ni_gme

    ig2s = 1
    ig2e = ni_gme+1

    ! Set neighbors
    nbpe = -1
    nbaw = -1
    nbpw = -1
    nbae = -1

  ENDIF

  ! Parameters for both Versions:
  ! Size of gme core
  max_gme_core = (ig1e-ig1s+5)*(ig2e-ig2s+5)
  IF (num_compute > 1) THEN
    CALL global_values (max_gme_core, 1, 'MAX', imp_integers, icomm_cart,  &
                        -1, yzerrmsg, izerror)
  ENDIF

  ig1sm1 = ig1s - 1
  ig1sm2 = ig1s - 2
  ig1ep1 = ig1e + 1
  ig1ep2 = ig1e + 2

  ig2sm1 = ig2s - 1
  ig2sm2 = ig2s - 2
  ig2ep1 = ig2e + 1
  ig2ep2 = ig2e + 2

  igg1s   = 0
  igg1e   = ni_gme
  igg2s   = 1
  igg2e   = ni_gme+1

  igg1sm1 = igg1s - 1
  igg1sm2 = igg1s - 2
  igg1ep1 = igg1e + 1
  igg1ep2 = igg1e + 2

  igg2sm1 = igg2s - 1
  igg2sm2 = igg2s - 2
  igg2ep1 = igg2e + 1
  igg2ep2 = igg2e + 2

! This was missing before, but ke1in is used for exchanging w in src_lm_fields
  ke1in = i3e_gme + 1

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE decompose_gme

!==============================================================================
!+ Subroutine that decomposes the coarse grid domain for distributed memory
!------------------------------------------------------------------------------

SUBROUTINE decompose_coarse_grid (yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine computes the start- and end-indices of the grid points 
!   from the coarse grid that are necessary to interpolate to the fine grid.
!   We add 2 additional rows and columns for the interpolation.
!   The data from the coarse grid then is decomposed accordingly.
!   The start- and end-indices of the coarse grid are stored in isubpos_coarse
!   with the same convention used for isubpos (of the fine grid): For every
!   processor the order is: i_ll, j_ll, i_ur, j_ur (ll: lower left; ur: upper
!   right).
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT(OUT) ::       &
    ierror                        ! error code

  CHARACTER (LEN=80)      , INTENT(OUT) ::       &
    yerrmsg                       ! error message

! Local variables
  INTEGER (KIND=iintegers) :: &
    nzbounds,                 & ! number of boundary lines for coarse grid
    izstart_in, izend_in, jzstart_in, jzend_in, nzy, nzx, nzproc,            &
    i, j, n, izerror, isendbuf(6), irecvbuf(6,0:num_compute-1)

  REAL    (KIND=ireals)    :: &
    zlats, zlons, zendlonintot,                  &
    zminlat_m, zmaxlat_m, zminlon_m, zmaxlon_m,  &
    zminlat_u, zmaxlat_u, zminlon_u, zmaxlon_u,  &
    zminlat_v, zmaxlat_v, zminlon_v, zmaxlon_v,  &
    zminlat,   zmaxlat,   zminlon,   zmaxlon  ,  &
    zminlon_tst, zmaxlon_tst,                    &
    zsouthminlat, znorthmaxlat,                  &
    zwest_minlon, zeast_maxlon,                  &
    rsendbuf(4), rrecvbuf(4,0:num_compute-1)

  CHARACTER (LEN=80)           ::   yzerror

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

  ierror   = 0_iintegers
  yerrmsg  = '   '

  izerror  = 0_iintegers
  yzerror = '  '

!------------------------------------------------------------------------------
!- Begin SUBROUTINE decompose_coarse_grid
!------------------------------------------------------------------------------

  nzbounds = 2     ! this number can be changed if it turns out that
                   ! the subdomains of the coarse grid are not big enough

!------------------------------------------------------------------------------
!- Section 2: Determine most northern, southern, western, eastern part of the 
!             domain and check whether coarse grid is big enough
!------------------------------------------------------------------------------

  ! For latitudes this can be done with MIN and MAX.
  ! For longitudes we have to be careful, if the domain is located around
  ! the date line (startlon_in_tot > endlon_in_tot)

  IF (startlon_in_tot > endlon_in_tot) THEN
    ! normalize longitude values to )0 ... 360)
    ! then MINVAL, MAXVAL can be used
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (lon_coarse_m(i,j) <= 0.0_ireals) THEN
          lon_coarse_m(i,j) = lon_coarse_m(i,j) + 360.0_ireals
        ENDIF
        IF (lon_coarse_u(i,j) <= 0.0_ireals) THEN
          lon_coarse_u(i,j) = lon_coarse_u(i,j) + 360.0_ireals
        ENDIF
        IF (lon_coarse_v(i,j) <= 0.0_ireals) THEN
          lon_coarse_v(i,j) = lon_coarse_v(i,j) + 360.0_ireals
        ENDIF
      ENDDO
    ENDDO

    ! the endlon from the total domain also has to be normalized
    zendlonintot = endlon_in_tot + 360.0_ireals
  ELSE
    zendlonintot = endlon_in_tot
  ENDIF

  ! for the mass grid points
  zminlat_m = MINVAL (lat_coarse_m(1:ie2lm,1:je2lm))
  zmaxlat_m = MAXVAL (lat_coarse_m(1:ie2lm,1:je2lm))
  zminlon_m = MINVAL (lon_coarse_m(1:ie2lm,1:je2lm))
  zmaxlon_m = MAXVAL (lon_coarse_m(1:ie2lm,1:je2lm))

  ! for the u-grid points
  zminlat_u = MINVAL (lat_coarse_u(1:ie2lm,1:je2lm))
  zmaxlat_u = MAXVAL (lat_coarse_u(1:ie2lm,1:je2lm))
  zminlon_u = MINVAL (lon_coarse_u(1:ie2lm,1:je2lm))
  zmaxlon_u = MAXVAL (lon_coarse_u(1:ie2lm,1:je2lm))

  ! for the v-grid points
  zminlat_v = MINVAL (lat_coarse_v(1:ie2lm,1:je2lm))
  zmaxlat_v = MAXVAL (lat_coarse_v(1:ie2lm,1:je2lm))
  zminlon_v = MINVAL (lon_coarse_v(1:ie2lm,1:je2lm))
  zmaxlon_v = MAXVAL (lon_coarse_v(1:ie2lm,1:je2lm))

  ! overall
  zminlat   = MIN (zminlat_m, zminlat_u, zminlat_v)
  zmaxlat   = MAX (zmaxlat_m, zmaxlat_u, zmaxlat_v)
  zminlon   = MIN (zminlon_m, zminlon_u, zminlon_v)
  zmaxlon   = MAX (zmaxlon_m, zmaxlon_u, zmaxlon_v)

  ! Check whether coarse grid is big enough
  IF (zminlat-nzbounds*dlat_in < startlat_in_tot) THEN
    ierror  = 2
    WRITE (*,'(A)')      &
       ' *** ERROR: Coarse domain is not big enough in the south:'
    WRITE (*,'(A,F9.3)') &
       '            startlat for coarse domain   : ',  startlat_in_tot
    WRITE (*,'(A,F9.3)') &
       '            necessary lat for fine domain: ',  zminlat-nzbounds*dlat_in
  ENDIF
  IF (zmaxlat+nzbounds*dlat_in > endlat_in_tot) THEN
    ierror  = 3
    WRITE (*,'(A)')      &
       ' *** ERROR: Coarse domain is not big enough in the north:'
    WRITE (*,'(A,F9.3)') &
       '            endlat for coarse domain     : ',  endlat_in_tot
    WRITE (*,'(A,F9.3)') &
       '            necessary lat for fine domain: ',  zmaxlat+nzbounds*dlat_in
  ENDIF
  IF (zminlon-nzbounds*dlon_in < startlon_in_tot) THEN
    ierror  = 4
    WRITE (*,'(A)')      &
       ' *** ERROR: Coarse domain is not big enough in the west:'
    WRITE (*,'(A,F9.3)') &
       '            startlon for coarse domain   : ',  startlon_in_tot
    WRITE (*,'(A,F9.3)') &
       '            necessary lon for fine domain: ',  zminlon-nzbounds*dlon_in
  ENDIF
  IF (zmaxlon+nzbounds*dlon_in > zendlonintot) THEN
    ierror  = 5
    WRITE (*,'(A)')      &
       ' *** ERROR: Coarse domain is not big enough in the east:'
    WRITE (*,'(A,F9.3)') &
       '            endlon for coarse domain     : ',  zendlonintot
    WRITE (*,'(A,F9.3)') &
       '            necessary lon for fine domain: ',  zmaxlon+nzbounds*dlon_in
  ENDIF

  IF (startlon_in_tot > endlon_in_tot) THEN
    ! re-calculate the original values
    DO j = 1, je2lm
      DO i = 1, ie2lm
        IF (lon_coarse_m(i,j) >= 180.0_ireals) THEN
          lon_coarse_m(i,j) = lon_coarse_m(i,j) - 360.0_ireals
        ENDIF
        IF (lon_coarse_u(i,j) >= 180.0_ireals) THEN
          lon_coarse_u(i,j) = lon_coarse_u(i,j) - 360.0_ireals
        ENDIF
        IF (lon_coarse_v(i,j) >= 180.0_ireals) THEN
          lon_coarse_v(i,j) = lon_coarse_v(i,j) - 360.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
!- Section 3: Indices of coarse grid covering this part of fine grid
!------------------------------------------------------------------------------

  IF (num_compute > 1) THEN
    ! to get reproducible results also for the Match-interpolation if isolated
    ! points have to be treated, the local domains at the boundaries have to
    ! take care that the same coarse grid part is used as in the sequential case.
    ! Therefore we have to gather the global min- and max-values from all boundary
    ! domains. For technical reasons we just gather all these values to every
    ! subdomain and the subdomains compute, what they need.
    rsendbuf(1) = zminlat
    rsendbuf(2) = zmaxlat
    rsendbuf(3) = zminlon
    rsendbuf(4) = zmaxlon

    ! Collect values to all nodes
    IF (num_compute > 1) THEN
      CALL gather_values (rsendbuf, rrecvbuf, 4, num_compute, imp_reals,  &
                          -1, icomm_cart, yzerror, izerror)

      zsouthminlat = zminlat
      znorthmaxlat = zmaxlat
      zwest_minlon = zminlon
      zeast_maxlon = zmaxlon

      DO nzx = 0, nprocx-1
        DO nzy = 0 , nprocy-1
          nzproc = nzx * nprocy + nzy

          ! get the minimal start latitude for all subdomains at the southern border
          ! (nzy = 0, nzx = 0,...,nprocx-1)
          IF (nzy == 0) THEN
            zsouthminlat = MIN(zsouthminlat, rrecvbuf(1,nzproc))
          ENDIF

          ! get the maximal end latitude for all subdomains at the northern border
          ! (nzy = nprocy-1, nzx = 0,...,nprocx-1)
          IF (nzy == nprocy-1) THEN
            znorthmaxlat = MAX(znorthmaxlat, rrecvbuf(2,nzproc))
          ENDIF

          ! get the minimal start longitude for all subdomains at the western border
          ! (nzx = 0, nzy = 0,...,nprocy-1)
          IF (nzx == 0) THEN
            zwest_minlon = MIN(zwest_minlon, rrecvbuf(3,nzproc))
          ENDIF

          ! get the maximal end longitude for all subdomains at the eastern border
          ! (nzx = nprocx-1, nzy = 0,...,nprocx-1)
          IF (nzx == nprocx-1) THEN
            zeast_maxlon = MAX(zeast_maxlon, rrecvbuf(4,nzproc))
          ENDIF

        ENDDO
      ENDDO

      ! now every subdomain knows all min and max values
      ! the subdomains at the boundaries have to store these values now
      DO nzx = 0, nprocx-1
        DO nzy = 0 , nprocy-1
          nzproc = nzx * nprocy + nzy
          IF ((nzy == 0) .AND. (my_cart_id == nzproc)) THEN
            zminlat = zsouthminlat
          ENDIF
          IF ((nzy == nprocy-1) .AND. (my_cart_id == nzproc)) THEN
            zmaxlat = znorthmaxlat
          ENDIF
          IF ((nzx == 0) .AND. (my_cart_id == nzproc)) THEN
            zminlon = zwest_minlon
          ENDIF
          IF ((nzx == nprocx-1) .AND. (my_cart_id == nzproc)) THEN
            zmaxlon = zeast_maxlon
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! looking for southern grid point
  DO j = 1, je_in_tot
    zlats  = startlat_in_tot + (j-1) * dlat_in
    IF (zlats >= zminlat) THEN
      jzstart_in = j-nzbounds-1
      EXIT
    ENDIF
  ENDDO

  ! looking for northern grid point
  DO j = 1, je_in_tot
    zlats  = startlat_in_tot + (j-1) * dlat_in
    IF (zlats >  zmaxlat) THEN
      jzend_in   = j+nzbounds
      EXIT
    ENDIF
  ENDDO

  ! looking for western grid point
  DO i = 1, ie_in_tot
    zlons  = startlon_in_tot + (i-1) * dlon_in
    IF (zlons >= zminlon) THEN
      izstart_in = i-nzbounds-1
      EXIT
    ENDIF
  ENDDO

  ! looking for eastern grid point
  DO i = 1, ie_in_tot
    zlons  = startlon_in_tot + (i-1) * dlon_in
    IF (zlons >  zmaxlon) THEN
      izend_in   = i+nzbounds
      EXIT
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- Section 4: Compute isubpos_coarse
!------------------------------------------------------------------------------

  ! Store values in sendbuf
  isendbuf (1) = izstart_in
  isendbuf (2) = jzstart_in
  isendbuf (3) = izend_in
  isendbuf (4) = jzend_in
  isendbuf (5) = 0
  isendbuf (6) = 0

  ! Collect values of isubpos_coarse to all nodes
  IF (num_compute > 1) THEN
    CALL gather_values (isendbuf, irecvbuf, 6, num_compute, imp_integers,  &
                        -1, icomm_cart, yzerror, izerror)
  ELSE
    irecvbuf(1:4,0) = isendbuf(1:4)
  ENDIF

  DO i = 0, num_compute-1
    isubpos_coarse (i, 1) = irecvbuf(1,i)
    isubpos_coarse (i, 2) = irecvbuf(2,i)
    isubpos_coarse (i, 3) = irecvbuf(3,i)
    isubpos_coarse (i, 4) = irecvbuf(4,i)
  ENDDO

!--------------------------------------------------------------------------------
!- Section 5: Determine size of this part of coarse grid
!--------------------------------------------------------------------------------

  ! size
  ie_in = izend_in - izstart_in + 1
  je_in = jzend_in - jzstart_in + 1
  ke_in = ke_in_tot - nlevskip
  ke1in = ke_in + 1

  ! Calculate ie_in_max and je_in_max
  isendbuf(1) = ie_in
  isendbuf(2) = je_in
  CALL global_values (isendbuf, 2, 'MAX', imp_integers, icomm_cart, -1,   &
                      yzerror, izerror)
  ie_in_max = isendbuf(1)
  je_in_max = isendbuf(2)


  ! start- and end-latitude and -longitude
  startlat_in = startlat_in_tot + (jzstart_in-1) * dlat_in
  endlat_in   = startlat_in_tot + (jzend_in  -1) * dlat_in
  startlon_in = startlon_in_tot + (izstart_in-1) * dlon_in
  endlon_in   = startlon_in_tot + (izend_in  -1) * dlon_in

  IF (idbg_level > 10) THEN
    WRITE (*,'(A,3I4,2F9.3)') 'LATs for subdomain:', jzstart_in, jzend_in, &
                                   je_in, startlat_in, endlat_in
    WRITE (*,'(A,3I4,2F9.3)') 'LONs for subdomain:', izstart_in, izend_in, &
                                   ie_in, startlon_in, endlon_in
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE decompose_coarse_grid

!==============================================================================
!==============================================================================
!+ Subroutine that decomposes the CM domain for distributed memory computers
!------------------------------------------------------------------------------

SUBROUTINE decompose_cm (yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine computes the start- and end-indices of the grid points
!   from the CM grid that are necessary to interpolate to the fine grid.
!   We add 2 additional rows and columns for the interpolation.
!   The data from the coarse grid then is decomposed accordingly.
!   The start- and end-indices of the coarse grid are stored in isubpos_coarse
!   with the same convention used for isubpos (of the fine grid): For every
!   processor the order is: i_ll, j_ll, i_ur, j_ur (ll: lower left; ur: upper
!   right).
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT(OUT) ::       &
    ierror                        ! error code

  CHARACTER (LEN=80)      , INTENT(OUT) ::       &
    yerrmsg                       ! error message

! Local variables
  INTEGER (KIND=iintegers) :: &
    nzbounds_e, nzbounds_w, nzbounds_s, nzbounds_n,       & 
               ! number of boundary lines for coarse grid
    izstart_in, izend_in, jzstart_in, jzend_in,  nzy, nzx, nzproc,          &
    i, j, n, izerror, isendbuf(6), irecvbuf(6,0:num_compute-1), &
    nufile     ! unit number of opened netCDF file

  REAL    (KIND=ireals)    :: &
    zlats, zlons,             &
    zminlat_m, zmaxlat_m, zminlon_m, zmaxlon_m,  &
    zminlat_u, zmaxlat_u, zminlon_u, zmaxlon_u,  &
    zminlat_v, zmaxlat_v, zminlon_v, zmaxlon_v,  &
    zminlat,   zmaxlat,   zminlon,   zmaxlon,    &
    zsouthminlat, znorthmaxlat,                  &
    zwest_minlon, zeast_maxlon,                  &
    rsendbuf(4), rrecvbuf(4,0:num_compute-1)

  CHARACTER (LEN=80)           ::   yzerror

  CHARACTER (LEN=100)        ::  &
    yzname      ! full path and name of the input-file

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

  ierror   = 0_iintegers
  yerrmsg  = '   '

  izerror  = 0_iintegers
  yzerror = '  '

!------------------------------------------------------------------------------
!- Begin SUBROUTINE decompose_cm
!------------------------------------------------------------------------------

  IF ((east_add_in + west_add_in + north_add_in + south_add_in) == 0) THEN
    nzbounds_e = 2     ! these numbers can be changed if it turns out that
    nzbounds_w = 2     ! the subdomains of the coarse grid are not big enough
    nzbounds_s = 2     !
    nzbounds_n = 2     !
  ELSE
    IF ((east_add_in + west_add_in) == 0) THEN
      nzbounds_e = 2
      nzbounds_w = 2
    ELSE
      nzbounds_e = 0
      nzbounds_w = 0
    ENDIF
    IF (south_add_in == 0) THEN
      nzbounds_s = 2
    ELSE
      nzbounds_s = 0
    ENDIF
    IF (north_add_in == 0) THEN
      nzbounds_n = 2
    ELSE
      nzbounds_n = 0
    ENDIF
  ENDIF
!------------------------------------------------------------------------------
!- Section 1: get latitudes and longitudes from netCDF Input file
!
!------------------------------------------------------------------------------

!  yzname = TRIM(yinext_cat)//TRIM(yinext_lfn)
!  CALL open_nc(nufile, yzname, 'r  ', icomm_cart, num_compute, &
!                 yzerror, ierror)
!  IF (ierror /= 0) THEN
!    PRINT *, 'Error in open_nc : ', TRIM(yzerror)
!    RETURN
!  ENDIF
  CALL read_nc_axis  (ie_in_tot, je_in_tot, ke_in, &
                      startlon_in_tot, startlat_in_tot, endlon_in_tot, endlat_in_tot,  &
                      pollon_in, pollat_in, lushift_in, lvshift_in,  &
                      east_add_in, west_add_in, south_add_in, north_add_in, &
                      icomm_cart, my_cart_id, num_compute,yzerror, ierror)
  IF (ierror /= 0) THEN
    PRINT *, 'Error in read_nc_axis : ', TRIM(yzerror)
    RETURN
  ENDIF
!  CALL close_nc (nufile, icomm_cart, num_compute, yzerror, ierror)
!  IF (ierror /= 0) THEN
!    PRINT *, 'Error in close_nc : ', TRIM(yzerror)
!    RETURN
!  ENDIF

  !DL: treat case where domain crosses date line
  DO i=2,ie_in_tot
    IF (longitudes_in(i) < longitudes_in(i-1)) longitudes_in(i)=longitudes_in(i) + 360.
  ENDDO

  WHERE ( lon_coarse_m < startlon_in_tot ) lon_coarse_m=lon_coarse_m+360.
  WHERE ( lon_coarse_u < startlon_in_tot ) lon_coarse_u=lon_coarse_u+360.
  WHERE ( lon_coarse_v < startlon_in_tot ) lon_coarse_v=lon_coarse_v+360.

!------------------------------------------------------------------------------
!- Section 2: Determine MIN and MAX of latitude and longitude
!             and check whether coarse grid is big enough
!------------------------------------------------------------------------------

  ! for the mass grid points
  zminlat_m = MINVAL (lat_coarse_m(1:ie2lm,1:je2lm))
  zmaxlat_m = MAXVAL (lat_coarse_m(1:ie2lm,1:je2lm))
  zminlon_m = MINVAL (lon_coarse_m(1:ie2lm,1:je2lm))
  zmaxlon_m = MAXVAL (lon_coarse_m(1:ie2lm,1:je2lm))

  ! for the u-grid points
  zminlat_u = MINVAL (lat_coarse_u(1:ie2lm,1:je2lm))
  zmaxlat_u = MAXVAL (lat_coarse_u(1:ie2lm,1:je2lm))
  zminlon_u = MINVAL (lon_coarse_u(1:ie2lm,1:je2lm))
  zmaxlon_u = MAXVAL (lon_coarse_u(1:ie2lm,1:je2lm))

  ! for the v-grid points
  zminlat_v = MINVAL (lat_coarse_v(1:ie2lm,1:je2lm))
  zmaxlat_v = MAXVAL (lat_coarse_v(1:ie2lm,1:je2lm))
  zminlon_v = MINVAL (lon_coarse_v(1:ie2lm,1:je2lm))
  zmaxlon_v = MAXVAL (lon_coarse_v(1:ie2lm,1:je2lm))

  ! overall
  zminlat   = MIN (zminlat_m, zminlat_u, zminlat_v)
  zmaxlat   = MAX (zmaxlat_m, zmaxlat_u, zmaxlat_v)
  zminlon   = MIN (zminlon_m, zminlon_u, zminlon_v)
  zmaxlon   = MAX (zmaxlon_m, zmaxlon_u, zmaxlon_v)

  ! Check whether coarse grid is big enough
  IF (zminlat < latitudes_in(1+nzbounds_s) .OR.             &
      zminlat < slatitudes_in(1+nzbounds_s)) THEN
    ierror  = 2
    yerrmsg = 'Coarse domain is not big enough'
    PRINT *, 'zminlat ', zminlat, ' < ', & 
      MAX(latitudes_in(1+nzbounds_s), slatitudes_in(1+nzbounds_s))
  ENDIF
  IF (zmaxlat > latitudes_in(je_in_tot-nzbounds_n) .OR.     &
      zmaxlat > slatitudes_in(je_in_tot-nzbounds_n)) THEN
    ierror  = 3
    yerrmsg = 'Coarse domain is not big enough'
    PRINT *, 'zmaxlat ', zmaxlat, ' > ', &
      MAX(latitudes_in(je_in_tot-nzbounds_n), slatitudes_in(je_in_tot-nzbounds_n))
  ENDIF
!US  IF (zminlon < longitudes_in(1+nzbounds_e) .OR.            &
!US      zminlon < slongitudes_in(1+nzbounds_e)) THEN
  IF (zminlon < longitudes_in(1+nzbounds_w) .OR.            &
      zminlon < slongitudes_in(1+nzbounds_w)) THEN
    ierror  = 4
    yerrmsg = 'Coarse domain is not big enough'
    PRINT *, 'zminlon ', zminlon, ' < ', &
      MAX(longitudes_in(1+nzbounds_w), slongitudes_in(1+nzbounds_w))
  ENDIF
!US  IF (zmaxlon > longitudes_in(ie_in_tot-nzbounds_w) .OR.    &
!US      zmaxlon > slongitudes_in(ie_in_tot-nzbounds_w)) THEN
  IF (zmaxlon > longitudes_in(ie_in_tot-nzbounds_e) .OR.    &
      zmaxlon > slongitudes_in(ie_in_tot-nzbounds_e)) THEN
   ierror  = 5
    yerrmsg = 'Coarse domain is not big enough'
    PRINT *, 'zmaxlon ', zmaxlon, ' > ', &
      MAX(longitudes_in(ie_in_tot-nzbounds_e), slongitudes_in(ie_in_tot-nzbounds_e))
  ENDIF

!------------------------------------------------------------------------------
!- Section 3: Indices of coarse grid covering this part of fine grid
!------------------------------------------------------------------------------

  IF (num_compute > 1) THEN
    ! to get reproducible results also for the Match-interpolation if isolated
    ! points have to be treated, the local domains at the boundaries have to
    ! take care that the same coarse grid part is used as in the sequential case.
    ! Therefore we have to gather the global min- and max-values from all boundary
    ! domains. For technical reasons we just gather all these values to every
    ! subdomain and the subdomains compute, what they need.
    rsendbuf(1) = zminlat
    rsendbuf(2) = zmaxlat
    rsendbuf(3) = zminlon
    rsendbuf(4) = zmaxlon

    ! Collect values to all nodes
    IF (num_compute > 1) THEN
      CALL gather_values (rsendbuf, rrecvbuf, 4, num_compute, imp_reals,  &
                          -1, icomm_cart, yzerror, izerror)
  
      zsouthminlat = zminlat
      znorthmaxlat = zmaxlat
      zwest_minlon = zminlon
      zeast_maxlon = zmaxlon 

      DO nzx = 0, nprocx-1
        DO nzy = 0 , nprocy-1
          nzproc = nzx * nprocy + nzy
  
          ! get the minimal start latitude for all subdomains at the southern border
          ! (nzy = 0, nzx = 0,...,nprocx-1)
          IF (nzy == 0) THEN
            zsouthminlat = MIN(zsouthminlat, rrecvbuf(1,nzproc))
          ENDIF  

          ! get the maximal end latitude for all subdomains at the northern border
          ! (nzy = nprocy-1, nzx = 0,...,nprocx-1)
          IF (nzy == nprocy-1) THEN
            znorthmaxlat = MAX(znorthmaxlat, rrecvbuf(2,nzproc))
          ENDIF 
  
          ! get the minimal start longitude for all subdomains at the western border
          ! (nzx = 0, nzy = 0,...,nprocy-1)
          IF (nzx == 0) THEN
            zwest_minlon = MIN(zwest_minlon, rrecvbuf(3,nzproc))
          ENDIF
  
          ! get the maximal end longitude for all subdomains at the eastern border
          ! (nzx = nprocx-1, nzy = 0,...,nprocx-1)
          IF (nzx == nprocx-1) THEN
            zeast_maxlon = MAX(zeast_maxlon, rrecvbuf(4,nzproc))
          ENDIF  

        ENDDO
      ENDDO
  
      ! now every subdomain knows all min and max values
      ! the subdomains at the boundaries have to store these values now
      DO nzx = 0, nprocx-1
        DO nzy = 0 , nprocy-1
          nzproc = nzx * nprocy + nzy
          IF ((nzy == 0) .AND. (my_cart_id == nzproc)) THEN
            zminlat = zsouthminlat
          ENDIF
          IF ((nzy == nprocy-1) .AND. (my_cart_id == nzproc)) THEN
            zmaxlat = znorthmaxlat
          ENDIF
          IF ((nzx == 0) .AND. (my_cart_id == nzproc)) THEN
            zminlon = zwest_minlon
          ENDIF
          IF ((nzx == nprocx-1) .AND. (my_cart_id == nzproc)) THEN
            zmaxlon = zeast_maxlon
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! looking for southern grid point
  DO j = 1, je_in_tot
    IF (latitudes_in(j) >= zminlat) THEN
      jzstart_in = j-nzbounds_s-1
      EXIT
   ENDIF
  ENDDO

  ! looking for northern grid point
  DO j = 1, je_in_tot
    IF (latitudes_in(j) >  zmaxlat) THEN
      jzend_in   = j+nzbounds_n
      EXIT
    ENDIF
  ENDDO

  ! looking for western grid point
  DO i = 1, ie_in_tot
    IF (longitudes_in(i) >= zminlon) THEN
      izstart_in = i-nzbounds_e-1
      EXIT
    ENDIF
  ENDDO

  ! looking for eastern grid point
  DO i = 1, ie_in_tot
    IF (longitudes_in(i) >  zmaxlon) THEN
      izend_in   = i+nzbounds_w
      EXIT
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- Section 4: Compute isubpos_coarse
!------------------------------------------------------------------------------

  ! Store values in sendbuf
  isendbuf (1) = izstart_in
  isendbuf (2) = jzstart_in
  isendbuf (3) = izend_in
  isendbuf (4) = jzend_in
  isendbuf (5) = 0
  isendbuf (6) = 0

  ! Collect values of isubpos_coarse to all nodes
  IF (num_compute > 1) THEN
    CALL gather_values (isendbuf, irecvbuf, 6, num_compute, imp_integers,  &
                        -1, icomm_cart, yzerror, izerror)
  ELSE
    irecvbuf(1:4,0) = isendbuf(1:4)
  ENDIF

  DO i = 0, num_compute-1
    isubpos_coarse (i, 1) = irecvbuf(1,i)
    isubpos_coarse (i, 2) = irecvbuf(2,i)
    isubpos_coarse (i, 3) = irecvbuf(3,i)
    isubpos_coarse (i, 4) = irecvbuf(4,i)
  ENDDO

!--------------------------------------------------------------------------------
!- Section 5: Determine size of this part of coarse grid
!--------------------------------------------------------------------------------

  ! size
  ie_in = izend_in - izstart_in + 1
  je_in = jzend_in - jzstart_in + 1
  ke_in = ke_in_tot
  ke1in = ke_in + 1

  ! Calculate ie_in_max and je_in_max
  isendbuf(1) = ie_in
  isendbuf(2) = je_in
  CALL global_values (isendbuf, 2, 'MAX', imp_integers, icomm_cart, -1,   &
                      yzerror, izerror)
  ie_in_max = isendbuf(1)
  je_in_max = isendbuf(2)


  ! start- and end-latitude and -longitude
  startlat_in = latitudes_in(jzstart_in)
  endlat_in   = latitudes_in(jzend_in)
  startlon_in = longitudes_in(izstart_in)
  endlon_in   = longitudes_in(izend_in)

  IF (idbg_level > 10) THEN
    WRITE (*,'(A,3I4,2F9.3)') 'LATs for subdomain:', jzstart_in, jzend_in, &
                                   je_in, startlat_in, endlat_in
    WRITE (*,'(A,3I4,2F9.3)') 'LONs for subdomain:', izstart_in, izend_in, &
                                   ie_in, startlon_in, endlon_in
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE decompose_cm

!==============================================================================
!==============================================================================
!+ Send control information about the decomposition to the root
!------------------------------------------------------------------------------

SUBROUTINE check_decomposition (ibuflen, ndebug, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine puts the variables computed in decompose_*** into
!   a buffer and sends it to the root-process where it is printed in a file
!   for debugging purposes.
!
! Method:
!   With routine gather_values.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT(IN)  ::       &
    ibuflen,                    & ! length of sending buffer
    ndebug                        ! unit number of debug output file

  INTEGER (KIND=iintegers), INTENT(OUT) ::       &
    ierror                        ! error code

  CHARACTER (LEN=40)      , INTENT(OUT) ::       &
    yerrmsg                       ! error message 

! Local variables
  INTEGER (KIND=iintegers)               ::      &
    isendbuf(ibuflen),   &        ! buffer for sending and receiving
    irecvbuf(ibuflen,0:num_compute-1)    

  REAL    (KIND=ireals)      :: &
    rsendbuf(4), rrecvbuf(4,0:num_compute-1)

  INTEGER (KIND=iintegers)   ::       &
    implcode,                   & ! for MPI error code
    nzsendcount, nzrecvcount, nzroot, noffset, n, l,    &
    itot, jtot

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE check_decomposition
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ! Initializations
  ierror      = 0
  nzroot      = 0
  nzsendcount = ibuflen
  nzrecvcount = nzsendcount

!------------------------------------------------------------------------------
!- Section 2:  Gather and print values for fine grid LM
!------------------------------------------------------------------------------

  ! Put data for LM decomposition into isendbuf
  isendbuf ( 1) = my_cart_id
  isendbuf ( 2) = my_world_id
  isendbuf ( 3) = my_cart_pos(1)
  isendbuf ( 4) = my_cart_pos(2)
  isendbuf ( 5) = my_cart_neigh(1)
  isendbuf ( 6) = my_cart_neigh(2)
  isendbuf ( 7) = my_cart_neigh(3)
  isendbuf ( 8) = my_cart_neigh(4)
  isendbuf ( 9) = isubpos(my_cart_id,1)
  isendbuf (10) = isubpos(my_cart_id,2)
  isendbuf (11) = isubpos(my_cart_id,3)
  isendbuf (12) = isubpos(my_cart_id,4)
  isendbuf (13) = ielm
  isendbuf (14) = jelm
  isendbuf (15) = kelm
  isendbuf (16) = ngp

  noffset = 16
  DO n = 1,ngp
    noffset = noffset + 1
    isendbuf (noffset) = igp(n)
    noffset = noffset + 1
    isendbuf (noffset) = jgp(n)
  ENDDO

  ! Gather the values
  CALL gather_values (isendbuf, irecvbuf, ibuflen    , num_compute,         &
                      imp_integers, nzroot, icomm_cart, yerrmsg, implcode )

  IF (my_cart_id == 0) THEN
    ! Print a headline in file YUDEBUG
    WRITE (ndebug, '(A2)')  '  '
    WRITE (ndebug, '(A)')                                                    &
          '      The decomposition for the LM grid was calculated as follows:'
    WRITE (ndebug, '(A)')                                                    &
          '      ============================================================'
    WRITE (ndebug, '(A2)')  '  '

    ! Print the information from all processes
    DO n = 0,num_compute-1
      WRITE (ndebug, '(A2)')  '  '
      WRITE (ndebug, '(A,I10,A,I10)')                                        &
      '       my_cart_id:  ',irecvbuf(1,n),'    my_world_id:  ',irecvbuf(2,n)
      WRITE (ndebug, '(A41,I2,A1,I2,A1)')                                    &
      '       Position in the cartesian grid:  (',                           &
                                          irecvbuf(3,n),',',irecvbuf(4,n),')'

      WRITE (ndebug, '(A)') '       LM Neighbors '
      WRITE (ndebug, '(T13,I5)')        irecvbuf(6,n)
      WRITE (ndebug, '(T8, I5,T18,I5)') irecvbuf(5,n), irecvbuf(7,n)
      WRITE (ndebug, '(T13,I5)')        irecvbuf(8,n)

      WRITE (ndebug, '(A57)')                                                &
      '       LM Location of this subdomain in the total domain:'
      WRITE (ndebug, '(A13,I5,A5,I5,A13,I5,A5,I5)')                          &
      '         i = ',irecvbuf( 9,n),',...,',irecvbuf(11,n),                 &
      '         j = ',irecvbuf(10,n),',...,',irecvbuf(12,n)

      WRITE (ndebug, '(A36)') '       Dimensions of this subdomain:'
      WRITE (ndebug, '(A16,I4,A12,I4,A12,I4)')                               &
      '         ielm = ',irecvbuf(13,n),'     jelm = ',irecvbuf(14,n),       &
      '     kelm = ',irecvbuf(15,n)

      WRITE (ndebug, '(A61,I4)')                                             &
      '       Number of selected grid points in this LM subdomain:  ',       &
                                               irecvbuf(16,n)
      noffset = 26 + nprocx + nprocy
      WRITE (ndebug, '(A)')                                                  &
      '          no.    i-local   j-local          i_total   j_total'

      noffset = 16
      DO l = 1,irecvbuf(16,n)
        noffset = noffset + 2

        ! Compute total indices
        itot = isubpos(n,1) - nboundlines - 1 + irecvbuf(noffset-1,n)
        jtot = isubpos(n,2) - nboundlines - 1 + irecvbuf(noffset  ,n)
        WRITE (ndebug, '(T10,I5,T18,I5,T28,I5,T45,I5,T55,I5)')               &
                 l, irecvbuf(noffset-1,n),irecvbuf(noffset,n), itot, jtot
      ENDDO

      WRITE (ndebug, '(A2)')  '  '

    ENDDO
  ENDIF

!------------------------------------------------------------------------------
!- Section 3:  Gather and print values for the coarse grids
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !- Section 3.1:  GME
  !----------------------------------------------------------------------------

  IF (lgme2lm) THEN
    ! Put data for GME decomposition into isendbuf
    isendbuf (1) = ig1s
    isendbuf (2) = ig1e
    isendbuf (3) = ig2s
    isendbuf (4) = ig2e
    isendbuf (5) = nbpe
    isendbuf (6) = nbaw
    isendbuf (7) = nbpw
    isendbuf (8) = nbae
  
    noffset = 8
    DO n = 0, nprocy
      noffset = noffset + 1
      isendbuf (noffset) = ilim1(n)
    ENDDO
    DO n = 0, nprocx
      noffset = noffset + 1
      isendbuf (noffset) = ilim2(n)
    ENDDO

    ! Gather the values
    CALL gather_values (isendbuf, irecvbuf, ibuflen    , num_compute,         &
                      imp_integers, nzroot, icomm_cart, yerrmsg, implcode )

    IF (my_cart_id == 0) THEN
      WRITE (ndebug, '(A2)')  '  '
      WRITE (ndebug, '(A2)')  '  '
      WRITE (ndebug, '(A2)')  '  '
      WRITE (ndebug, '(A)')                                                  &
        '      The decomposition for the GME grid was calculated as follows:'
      WRITE (ndebug, '(A)')                                                  &
        '      ============================================================='
      WRITE (ndebug, '(A2)')  '  '

      ! Print the information from all processes for GME decomposition
      DO n = 0,num_compute-1
        WRITE (ndebug, '(A,I5)') '  Subdomain:  ', n
        WRITE (ndebug, '(T8,A,T35,A)')                                     &
                       'GME Neighbors:', 'Location in the total domain:'
        WRITE (ndebug, '(T8,A)')  '        /\       '
        WRITE (ndebug, '(T8,I5,A,I5,T40,A,I5,A,I5)')                       &
                        irecvbuf(7,n), '  /  \ ', irecvbuf(5,n),            &
                      'j1 = ',irecvbuf(1,n),',...,',irecvbuf(2,n)
        WRITE (ndebug, '(T8,A)')  '      /    \     '
        WRITE (ndebug, '(T8,A)')  '      \    /     '
        WRITE (ndebug, '(T8,I5,A,I5,T40,A,I5,A,I5)')                       &
                        irecvbuf(6,n), '  \  / ', irecvbuf(8,n),            &
                      'j2 = ',irecvbuf(3,n),',...,',irecvbuf(4,n)
        WRITE (ndebug, '(T8,A)')  '        \/       '
        WRITE (ndebug, '(A2)')  '  '
      ENDDO

      ! Only from PE 0:
      DO n = 0, nprocy
        WRITE (ndebug, '(A,I3,A,I3)') '      ilim1 (',n,') = ',ilim1(n)
      ENDDO
      DO n = 0, nprocx
        WRITE (ndebug, '(A,I3,A,I3)') '      ilim2 (',n,') = ',ilim2(n)
      ENDDO
      WRITE (ndebug, '(A2)')  '  '
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  !- Section 3.2:  EC or LM
  !----------------------------------------------------------------------------

  IF (llm2lm .OR. lec2lm .OR. lcm2lm) THEN
    ! Put data for coarse grid decomposition into sendbuf
    isendbuf(1) = isubpos_coarse(my_cart_id, 1)
    isendbuf(2) = isubpos_coarse(my_cart_id, 2)
    isendbuf(3) = ie_in
    isendbuf(4) = isubpos_coarse(my_cart_id, 3)
    isendbuf(5) = isubpos_coarse(my_cart_id, 4)
    isendbuf(6) = je_in

    rsendbuf(1) = startlat_in
    rsendbuf(2) = endlat_in
    rsendbuf(3) = startlon_in
    rsendbuf(4) = endlon_in

    ! Gather the values
    CALL gather_values (isendbuf, irecvbuf, ibuflen, num_compute,         &
                     imp_integers, nzroot, icomm_cart, yerrmsg, implcode )

    CALL gather_values (rsendbuf, rrecvbuf, 4, num_compute,         &
                        imp_reals, nzroot, icomm_cart, yerrmsg, implcode )

    IF (my_cart_id == 0) THEN
      ! Loop over all processors
      WRITE (ndebug, '(A)')    '   '
      WRITE (ndebug, '(A)')    '   '
      WRITE (ndebug, '(A)')    '      Data decomposition of coarse grid'
      WRITE (ndebug, '(A)')    '      ================================='
      WRITE (ndebug, '(A,I5)')                                        &
        '                   start       end       start     end    size'
      DO n = 0, num_compute-1
        WRITE (ndebug, '(A)')    '     '
        WRITE (ndebug, '(A,I5)') '  Subdomain:  ', n
        WRITE (ndebug, '(A,F10.2,F10.2,A,I5,I8,I8)')                  &
          '    latitude: ', rrecvbuf(1,n), rrecvbuf(2,n),             &
          '    j: ', irecvbuf(2,n), irecvbuf(5,n), irecvbuf(6,n)
        WRITE (ndebug, '(A,F10.2,F10.2,A,I5,I8,I8)')                  &
          '    longitude:', rrecvbuf(3,n), rrecvbuf(4,n),             &
          '    i: ', irecvbuf(1,n), irecvbuf(4,n), irecvbuf(3,n)
      ENDDO
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE check_decomposition

!==============================================================================

END MODULE src_decomposition
