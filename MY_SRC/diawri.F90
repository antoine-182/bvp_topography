MODULE diawri
   !!======================================================================
   !!                     ***  MODULE  diawri  ***
   !! Ocean diagnostics :  write ocean output files
   !!=====================================================================
   !! History :  OPA  ! 1991-03  (M.-A. Foujols)  Original code
   !!            4.0  ! 1991-11  (G. Madec)
   !!                 ! 1992-06  (M. Imbard)  correction restart file
   !!                 ! 1992-07  (M. Imbard)  split into diawri and rstwri
   !!                 ! 1993-03  (M. Imbard)  suppress writibm
   !!                 ! 1998-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
   !!                 ! 1999-02  (E. Guilyardi)  name of netCDF files + variables
   !!            8.2  ! 2000-06  (M. Imbard)  Original code (diabort.F)
   !!   NEMO     1.0  ! 2002-06  (A.Bozec, E. Durand)  Original code (diainit.F)
   !!             -   ! 2002-09  (G. Madec)  F90: Free form and module
   !!             -   ! 2002-12  (G. Madec)  merge of diabort and diainit, F90
   !!                 ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2008-11  (B. Lemaire) creation from old diawri
   !!            3.7  ! 2014-01  (G. Madec) remove eddy induced velocity from no-IOM output
   !!                 !                     change name of output variables in dia_wri_state
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_wri       : create the standart output files
   !!   dia_wri_state : create an output NetCDF file for a single instantaeous ocean state and forcing fields
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE isf_oce
   USE isfcpl
   USE abl            ! abl variables in case ln_abl = .true.
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE dianam         ! build name of file (routine)
   USE diahth         ! thermocline diagnostics
   USE dynadv   , ONLY: ln_dynadv_vec
   USE icb_oce        ! Icebergs
   USE icbdia         ! Iceberg budgets
   USE ldftra         ! lateral physics: eddy diffusivity coef.
   USE ldfdyn         ! lateral physics: eddy viscosity   coef.
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE sbcssr         ! restoring term toward SST/SSS climatology
   USE sbcwave        ! wave parameters
   USE wet_dry        ! wetting and drying
   USE zdf_oce        ! ocean vertical physics
   USE zdfdrg         ! ocean vertical physics: top/bottom friction
   USE zdfmxl         ! mixed layer
   !
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager ! I/O manager
   USE dia25h         ! 25h Mean output
   USE iom            !
   USE ioipsl         !

#if defined key_si3
   USE ice
   USE icewri
#endif
   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary
   USE diu_bulk        ! diurnal warm layer
   USE diu_coolskin    ! Cool skin

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_wri                 ! routines called by step.F90
   PUBLIC   dia_wri_state
   PUBLIC   dia_wri_alloc           ! Called by nemogcm module
#if ! defined key_iomput
   PUBLIC   dia_wri_alloc_abl       ! Called by sbcabl  module (if ln_abl = .true.)
#endif
   INTEGER ::   nid_T, nz_T, nh_T, ndim_T, ndim_hT   ! grid_T file
   INTEGER ::          nb_T              , ndim_bT   ! grid_T file
   INTEGER ::   nid_U, nz_U, nh_U, ndim_U, ndim_hU   ! grid_U file
   INTEGER ::   nid_V, nz_V, nh_V, ndim_V, ndim_hV   ! grid_V file
   INTEGER ::   nid_F, nz_F, nh_F, ndim_F, ndim_hF   ! grid_F file
   INTEGER ::   nid_W, nz_W, nh_W                    ! grid_W file
   INTEGER ::   nid_A, nz_A, nh_A, ndim_A, ndim_hA   ! grid_ABL file
   INTEGER ::   ndex(1)                              ! ???
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hT, ndex_hU, ndex_hV, ndex_hF
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hA, ndex_A ! ABL
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_T, ndex_U, ndex_V, ndex_F
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_bT

   !! * Substitutions
#  include "do_loop_substitute.h90"
!!st12
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: diawri.F90 13416 2020-08-20 10:10:55Z gm $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_iomput'                                        use IOM library
   !!----------------------------------------------------------------------
   INTEGER FUNCTION dia_wri_alloc()
      !
      dia_wri_alloc = 0
      !
   END FUNCTION dia_wri_alloc


   SUBROUTINE dia_wri( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri  ***
      !!
      !! ** Purpose :   Standard output of opa: dynamics and tracer fields
      !!      NETCDF format is used by default
      !!
      !! ** Method  :  use iom_put
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kmm     ! ocean time level index
      !!
      INTEGER ::   ji, jj, jk       ! dummy loop indices
      INTEGER ::   ikbot            ! local integer
      REAL(wp)::   zztmp , zztmpx   ! local scalar
      REAL(wp)::   zztmp2, zztmpy   !   -      -
      REAL(wp)::   ze3              !   -      -
      REAL(wp), DIMENSION(jpi,jpj)     ::   z2d   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   z3d   ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dia_wri')
      !
      ! Output the initial state and forcings
      IF( ninist == 1 ) THEN
         CALL dia_wri_state( Kmm, 'output.init' )
         ninist = 0
      ENDIF

      ! Output of initial vertical scale factor
      CALL iom_put("e3t_0", e3t_0(:,:,:) )
      CALL iom_put("e3u_0", e3u_0(:,:,:) )
      CALL iom_put("e3v_0", e3v_0(:,:,:) )
      CALL iom_put("e3f_0", e3f_0(:,:,:) )
      !
      IF ( iom_use("ahmts") ) CALL iom_put("ahmts", ahmt(:,:,1) )
      IF ( iom_use("ahmfs") ) CALL iom_put("ahmfs", ahmf(:,:,1) )
      !
!!st13
#if ! defined key_qco
      IF ( iom_use("e3t") .OR. iom_use("e3tdef") ) THEN  ! time-varying e3t
         DO jk = 1, jpk
            z3d(:,:,jk) =  e3t(:,:,jk,Kmm)
         END DO
         CALL iom_put( "e3t"     ,     z3d(:,:,:) )
         CALL iom_put( "e3tdef"  , ( ( z3d(:,:,:) - e3t_0(:,:,:) ) / e3t_0(:,:,:) * 100 * tmask(:,:,:) ) ** 2 )
      ENDIF
      IF ( iom_use("e3u") ) THEN                         ! time-varying e3u
         DO jk = 1, jpk
            z3d(:,:,jk) =  e3u(:,:,jk,Kmm)
         END DO
         CALL iom_put( "e3u" , z3d(:,:,:) )
      ENDIF
      IF ( iom_use("e3v") ) THEN                         ! time-varying e3v
         DO jk = 1, jpk
            z3d(:,:,jk) =  e3v(:,:,jk,Kmm)
         END DO
         CALL iom_put( "e3v" , z3d(:,:,:) )
      ENDIF
      IF ( iom_use("e3w") ) THEN                         ! time-varying e3w
         DO jk = 1, jpk
            z3d(:,:,jk) =  e3w(:,:,jk,Kmm)
         END DO
         CALL iom_put( "e3w" , z3d(:,:,:) )
      ENDIF
#endif
!!st13

# if defined key_bvp
      IF( iom_use("rpo" ) )   CALL iom_put( "rpo",   rpo(:,:,1) )
      IF( iom_use("rpou") )   CALL iom_put( "rpou", rpou(:,:,1) )
      IF( iom_use("rpov") )   CALL iom_put( "rpov", rpov(:,:,1) )
      IF( iom_use("rpof") )   CALL iom_put( "rpof", rpof(:,:,1) )
      IF( iom_use("bmpt") )   CALL iom_put( "bmpt", bmpt(:,:,1) )
      IF( iom_use("bmpu") )   CALL iom_put( "bmpu", bmpu(:,:,1) )
      IF( iom_use("bmpv") )   CALL iom_put( "bmpv", bmpv(:,:,1) )
#endif

# if defined key_bvp_bath
      IF( iom_use("batht") )   CALL iom_put( "batht", batht(:,:) )
#endif

# if ! defined key_bvp
      IF( ll_wd ) THEN
         CALL iom_put( "ssh" , (ssh(:,:,Kmm)+ssh_ref)*tmask(:,:,1) )   ! sea surface height (brought back to the reference used for wetting and drying)
      ELSE
         CALL iom_put( "ssh" , ssh(:,:,Kmm) )                          ! sea surface height
      ENDIF
# else
      IF( ll_wd ) THEN
         CALL iom_put( "ssh" , (r1_rpo(:,:,1) * ssh(:,:,Kmm)+ssh_ref)*tmask(:,:,1) )   ! sea surface height (brought back to the reference used for wetting and drying)
      ELSE
         CALL iom_put( "ssh" , r1_rpo(:,:,1) * ssh(:,:,Kmm) )                          ! sea surface height
      ENDIF

#endif

      !
      IF( iom_use("ht") )   &                  ! water column at t-point
         CALL iom_put( "ht" , ht_0(:,:) + ssh(:,:,Kmm) )
      !
      IF( iom_use("hu") )   THEN                   ! water column at u-point
         z2d(:,:) = 0._wp
         DO_2D_10_10
            z2d(ji,jj) = 0.5_wp * ( e3t(ji  ,jj,1,Kmm) * e1e2t(ji  ,jj)  &
               &                  + e3t(ji+1,jj,1,Kmm) * e1e2t(ji+1,jj)  ) * r1_e1e2u(ji,jj)
         END_2D
         CALL lbc_lnk( 'diawri', z2d, 'U', 1._wp )
         CALL iom_put( "hu", z2d )
      ENDIF
     !
      IF( iom_use("hv") )   THEN                  ! water column at v-point
         z2d(:,:) = 0._wp
         DO_2D_10_10
            z2d(ji,jj) = 0.5_wp * (  e3t(ji,jj+1,1,Kmm) * e1e2t(ji,jj+1)   &
              &                    + e3t(ji,jj  ,1,Kmm) * e1e2t(ji,jj  )   ) * r1_e1e2v(ji,jj)
         END_2D
         CALL lbc_lnk( 'diawri', z2d, 'V', 1._wp )
         CALL iom_put( "hv", z2d )
      ENDIF
      !
      IF( iom_use("hf") )  THEN                    ! water column at f-point
         z2d(:,:) = 0._wp
         DO_2D_10_10
            z2d(ji,jj) = 0.25_wp * (  e3t(ji,jj+1,1,Kmm) * e1e2t(ji,jj+1) + e3t(ji+1,jj+1,1,Kmm) * e1e2t(ji+1,jj+1)    &
              &                     + e3t(ji,jj  ,1,Kmm) * e1e2t(ji,jj  ) + e3t(ji+1,jj  ,1,Kmm) * e1e2t(ji+1,jj  )  ) * r1_e1e2f(ji,jj)
         END_2D
         z2d(:,:) = z2d(:,:) * ssfmask(:,:)
         CALL lbc_lnk( 'diawri', z2d, 'F', 1._wp )
         CALL iom_put( "hf", z2d )
      ENDIF
!!an

      IF( iom_use("wetdep") )   &                  ! wet depth
         CALL iom_put( "wetdep" , ht_0(:,:) + ssh(:,:,Kmm) )

      IF ( iom_use("taubot") ) THEN                ! bottom stress
         zztmp = rho0 * 0.25
         z2d(:,:) = 0._wp
         DO_2D_00_00
            zztmp2 = (  ( rCdU_bot(ji+1,jj)+rCdU_bot(ji  ,jj) ) * uu(ji  ,jj,mbku(ji  ,jj),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji  ,jj)+rCdU_bot(ji-1,jj) ) * uu(ji-1,jj,mbku(ji-1,jj),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj  ) ) * vv(ji,jj  ,mbkv(ji,jj  ),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji,jj  )+rCdU_bot(ji,jj-1) ) * vv(ji,jj-1,mbkv(ji,jj-1),Kmm)  )**2
            z2d(ji,jj) = zztmp * SQRT( zztmp2 ) * tmask(ji,jj,1)
            !
         END_2D
         CALL lbc_lnk( 'diawri', z2d, 'T', 1. )
         CALL iom_put( "taubot", z2d )
      ENDIF
      !
      CALL iom_put(  "ssu" , uu(:,:,1,Kmm) )            ! surface i-current
      CALL iom_put(  "ssv" , vv(:,:,1,Kmm) )            ! surface j-current
      CALL iom_put(  "woce", ww )                       ! vertical velocity
      !

      IF ( iom_use("sKE") ) THEN                        ! surface kinetic energy at T point
         z2d(:,:) = 0._wp
         DO_2D_00_00
            z2d(ji,jj) =    0.25_wp * ( uu(ji  ,jj,1,Kmm) * uu(ji  ,jj,1,Kmm) * e1e2u(ji  ,jj) * e3u(ji  ,jj,1,Kmm)  &
               &                      + uu(ji-1,jj,1,Kmm) * uu(ji-1,jj,1,Kmm) * e1e2u(ji-1,jj) * e3u(ji-1,jj,1,Kmm)  &
               &                      + vv(ji,jj  ,1,Kmm) * vv(ji,jj  ,1,Kmm) * e1e2v(ji,jj  ) * e3v(ji,jj  ,1,Kmm)  &
               &                      + vv(ji,jj-1,1,Kmm) * vv(ji,jj-1,1,Kmm) * e1e2v(ji,jj-1) * e3v(ji,jj-1,1,Kmm)  )  &
               &                    * r1_e1e2t(ji,jj) / e3t(ji,jj,1,Kmm) * tmask(ji,jj,1)
         END_2D
         !
         CALL lbc_lnk( 'diawri', z2d, 'T', 1. )
         IF ( iom_use("sKE" ) )  CALL iom_put( "sKE" , z2d )

      ENDIF

!!an sKEf
      IF ( iom_use("sKEf") ) THEN                        ! surface kinetic energy at F point
         z2d(:,:) = 0._wp                                ! CAUTION : only valid in SWE, not with bathymetry
         DO_2D_00_00
            z2d(ji,jj) =    0.25_wp * ( uu(ji,jj  ,1,Kmm) * uu(ji,jj  ,1,Kmm) * e1e2u(ji,jj  ) * e3u(ji,jj  ,1,Kmm)  &
               &                      + uu(ji,jj+1,1,Kmm) * uu(ji,jj+1,1,Kmm) * e1e2u(ji,jj+1) * e3u(ji,jj+1,1,Kmm)  &
               &                      + vv(ji  ,jj,1,Kmm) * vv(ji,jj  ,1,Kmm) * e1e2v(ji  ,jj) * e3v(ji  ,jj,1,Kmm)  &
               &                      + vv(ji+1,jj,1,Kmm) * vv(ji+1,jj,1,Kmm) * e1e2v(ji+1,jj) * e3v(ji+1,jj,1,Kmm)  )  &
               &                    * r1_e1e2f(ji,jj) / e3f(ji,jj,1) * ssfmask(ji,jj)
         END_2D
         !
         CALL lbc_lnk( 'diawri', z2d, 'F', 1. )
         CALL iom_put( "sKEf", z2d )
      ENDIF
!!an
      !
      CALL iom_put( "hdiv", hdiv )                  ! Horizontal divergence
      !
      ! Output of vorticity terms
      IF ( iom_use("relvor")    .OR. iom_use("plavor")    .OR.   &
         & iom_use("relpotvor") .OR. iom_use("abspotvor") .OR.   &
         & iom_use("enstro")                                        ) THEN
         !
         z2d(:,:) = 0._wp
         ze3 = 0._wp
         DO_2D_10_10
            z2d(ji,jj) = (   e2v(ji+1,jj  ) * vv(ji+1,jj  ,1,Kmm) - e2v(ji,jj) * vv(ji,jj,1,Kmm)    &
            &              - e1u(ji  ,jj+1) * uu(ji  ,jj+1,1,Kmm) + e1u(ji,jj) * uu(ji,jj,1,Kmm)  ) * r1_e1e2f(ji,jj)
         END_2D
         CALL lbc_lnk( 'diawri', z2d, 'F', 1. )
         CALL iom_put( "relvor", z2d )                  ! relative vorticity ( zeta )
         !
         CALL iom_put( "plavor", ff_f )                 ! planetary vorticity ( f )
         !
         DO_2D_10_10
            ze3 = (  e3t(ji,jj+1,1,Kmm) * e1e2t(ji,jj+1) + e3t(ji+1,jj+1,1,Kmm) * e1e2t(ji+1,jj+1)    &
              &    + e3t(ji,jj  ,1,Kmm) * e1e2t(ji,jj  ) + e3t(ji+1,jj  ,1,Kmm) * e1e2t(ji+1,jj  )  ) * r1_e1e2f(ji,jj)
            IF( ze3 /= 0._wp ) THEN   ;   ze3 = 4._wp / ze3
            ELSE                      ;   ze3 = 0._wp
            ENDIF
            z2d(ji,jj) = ze3 * z2d(ji,jj)
         END_2D
         CALL lbc_lnk( 'diawri', z2d, 'F', 1. )
         CALL iom_put( "relpotvor", z2d )                  ! relative potential vorticity (zeta/h)
         !
         DO_2D_10_10
            ze3 = (  e3t(ji,jj+1,1,Kmm) * e1e2t(ji,jj+1) + e3t(ji+1,jj+1,1,Kmm) * e1e2t(ji+1,jj+1)    &
              &    + e3t(ji,jj  ,1,Kmm) * e1e2t(ji,jj  ) + e3t(ji+1,jj  ,1,Kmm) * e1e2t(ji+1,jj  )  ) * r1_e1e2f(ji,jj)
            IF( ze3 /= 0._wp ) THEN   ;   ze3 = 4._wp / ze3
            ELSE                      ;   ze3 = 0._wp
            ENDIF
            z2d(ji,jj) = ze3 * ff_f(ji,jj) + z2d(ji,jj)
         END_2D
         CALL lbc_lnk( 'diawri', z2d, 'F', 1. )
         CALL iom_put( "abspotvor", z2d )                  ! absolute potential vorticity ( q )
         !
         DO_2D_10_10
            z2d(ji,jj) = 0.5_wp * z2d(ji,jj)  * z2d(ji,jj)
         END_2D
         CALL lbc_lnk( 'diawri', z2d, 'F', 1. )
         CALL iom_put( "enstro", z2d )                        ! potential enstrophy ( 1/2*q2 )
         !
      ENDIF

      !
      IF( ln_timing )   CALL timing_stop('dia_wri')
      !
   END SUBROUTINE dia_wri

#else
   !!----------------------------------------------------------------------
   !!   Default option                                  use IOIPSL  library
   !!----------------------------------------------------------------------

   INTEGER FUNCTION dia_wri_alloc()
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ierr
      !!----------------------------------------------------------------------
      IF( nn_write == -1 ) THEN
         dia_wri_alloc = 0
      ELSE
         ierr = 0
         ALLOCATE( ndex_hT(jpi*jpj) , ndex_T(jpi*jpj*jpk) ,     &
            &      ndex_hU(jpi*jpj) , ndex_U(jpi*jpj*jpk) ,     &
            &      ndex_hV(jpi*jpj) , ndex_V(jpi*jpj*jpk) , STAT=ierr(1) )
         !
         dia_wri_alloc = MAXVAL(ierr)
         CALL mpp_sum( 'diawri', dia_wri_alloc )
         !
      ENDIF
      !
   END FUNCTION dia_wri_alloc

   INTEGER FUNCTION dia_wri_alloc_abl()
      !!----------------------------------------------------------------------
	  ALLOCATE(   ndex_hA(jpi*jpj), ndex_A (jpi*jpj*jpkam1), STAT=dia_wri_alloc_abl)
      CALL mpp_sum( 'diawri', dia_wri_alloc_abl )
      !
   END FUNCTION dia_wri_alloc_abl


   SUBROUTINE dia_wri( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri  ***
      !!
      !! ** Purpose :   Standard output of opa: dynamics and tracer fields
      !!      NETCDF format is used by default
      !!
      !! ** Method  :   At the beginning of the first time step (nit000),
      !!      define all the NETCDF files and fields
      !!      At each time step call histdef to compute the mean if ncessary
      !!      Each nn_write time step, output the instantaneous or mean fields
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kmm  ! ocean time level index
      !
      LOGICAL ::   ll_print = .FALSE.                        ! =T print and flush numout
      CHARACTER (len=40) ::   clhstnam, clop, clmx           ! local names
      INTEGER  ::   inum = 11                                ! temporary logical unit
      INTEGER  ::   ji, jj, jk                               ! dummy loop indices
      INTEGER  ::   ierr                                     ! error code return from allocation
      INTEGER  ::   iimi, iima, ipk, it, itmod, ijmi, ijma   ! local integers
      INTEGER  ::   ipka                                     ! ABL
      INTEGER  ::   jn, ierror                               ! local integers
      REAL(wp) ::   zsto, zout, zmax, zjulian                ! local scalars
      !
      REAL(wp), DIMENSION(jpi,jpj)   :: zw2d       ! 2D workspace
!!st14
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zw3d, ze3t, zgdept       ! 3D workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zw3d_abl   ! ABL 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( ninist == 1 ) THEN     !==  Output the initial state and forcings  ==!
         CALL dia_wri_state( Kmm, 'output.init' )
         ninist = 0
      ENDIF
      !
      IF( nn_write == -1 )   RETURN   ! we will never do any output
      !
      IF( ln_timing )   CALL timing_start('dia_wri')
      !
      ! 0. Initialisation
      ! -----------------

      ll_print = .FALSE.                  ! local variable for debugging
      ll_print = ll_print .AND. lwp

      ! Define frequency of output and means
      clop = "x"         ! no use of the mask value (require less cpu time and otherwise the model crashes)
#if defined key_diainstant
      zsto = nn_write * rn_Dt
      clop = "inst("//TRIM(clop)//")"
#else
      zsto=rn_Dt
      clop = "ave("//TRIM(clop)//")"
#endif
      zout = nn_write * rn_Dt
      zmax = ( nitend - nit000 + 1 ) * rn_Dt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk
      IF(ln_abl) ipka = jpkam1

      ! define time axis
      it = kt
      itmod = kt - nit000 + 1
!!st15
      ! store e3t for subsitute
      DO jk = 1, jpk
         ze3t  (:,:,jk) =  e3t  (:,:,jk,Kmm)
         zgdept(:,:,jk) =  gdept(:,:,jk,Kmm)
      END DO
!!st15

      ! 1. Define NETCDF files and fields at beginning of first time step
      ! -----------------------------------------------------------------

      IF( kt == nit000 ) THEN

         ! Define the NETCDF files (one per grid)

         ! Compute julian date from starting date of the run
         CALL ymds2ju( nyear, nmonth, nday, rn_Dt, zjulian )
         zjulian = zjulian - adatrj   !   set calendar origin to the beginning of the experiment
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'Date 0 used :', nit000, ' YEAR ', nyear,   &
            &                    ' MONTH ', nmonth, ' DAY ', nday, 'Julian day : ', zjulian
         IF(lwp)WRITE(numout,*) ' indexes of zoom = ', iimi, iima, ijmi, ijma,   &
                                 ' limit storage in depth = ', ipk

         ! WRITE root name in date.file for use by postpro
         IF(lwp) THEN
            CALL dia_nam( clhstnam, nn_write,' ' )
            CALL ctl_opn( inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
            WRITE(inum,*) clhstnam
            CLOSE(inum)
         ENDIF

         ! Define the T grid FILE ( nid_T )

         CALL dia_nam( clhstnam, nn_write, 'grid_T' )
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam    ! filename
         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           &  ! Horizontal grid: glamt and gphit
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, rn_Dt, nh_T, nid_T, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_T, "deptht", "Vertical T levels",      &  ! Vertical grid: gdept
            &           "m", ipk, gdept_1d, nz_T, "down" )
         !                                                            ! Index of ocean points
         CALL wheneq( jpi*jpj*ipk, tmask, 1, 1., ndex_T , ndim_T  )      ! volume
         CALL wheneq( jpi*jpj    , tmask, 1, 1., ndex_hT, ndim_hT )      ! surface
         !
         IF( ln_icebergs ) THEN
            !
            !! allocation cant go in dia_wri_alloc because ln_icebergs is only set after
            !! that routine is called from nemogcm, so do it here immediately before its needed
            ALLOCATE( ndex_bT(jpi*jpj*nclasses), STAT=ierror )
            CALL mpp_sum( 'diawri', ierror )
            IF( ierror /= 0 ) THEN
               CALL ctl_stop('dia_wri: failed to allocate iceberg diagnostic array')
               RETURN
            ENDIF
            !
            !! iceberg vertical coordinate is class number
            CALL histvert( nid_T, "class", "Iceberg class",      &  ! Vertical grid: class
               &           "number", nclasses, class_num, nb_T )
            !
            !! each class just needs the surface index pattern
            ndim_bT = 3
            DO jn = 1,nclasses
               ndex_bT((jn-1)*jpi*jpj+1:jn*jpi*jpj) = ndex_hT(1:jpi*jpj)
            ENDDO
            !
         ENDIF

         ! Define the U grid FILE ( nid_U )

         CALL dia_nam( clhstnam, nn_write, 'grid_U' )
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam    ! filename
         CALL histbeg( clhstnam, jpi, glamu, jpj, gphiu,           &  ! Horizontal grid: glamu and gphiu
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, rn_Dt, nh_U, nid_U, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_U, "depthu", "Vertical U levels",      &  ! Vertical grid: gdept
            &           "m", ipk, gdept_1d, nz_U, "down" )
         !                                                            ! Index of ocean points
         CALL wheneq( jpi*jpj*ipk, umask, 1, 1., ndex_U , ndim_U  )      ! volume
         CALL wheneq( jpi*jpj    , umask, 1, 1., ndex_hU, ndim_hU )      ! surface

         ! Define the V grid FILE ( nid_V )

         CALL dia_nam( clhstnam, nn_write, 'grid_V' )                   ! filename
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam
         CALL histbeg( clhstnam, jpi, glamv, jpj, gphiv,           &  ! Horizontal grid: glamv and gphiv
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, rn_Dt, nh_V, nid_V, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_V, "depthv", "Vertical V levels",      &  ! Vertical grid : gdept
            &          "m", ipk, gdept_1d, nz_V, "down" )
         !                                                            ! Index of ocean points
         CALL wheneq( jpi*jpj*ipk, vmask, 1, 1., ndex_V , ndim_V  )      ! volume
         CALL wheneq( jpi*jpj    , vmask, 1, 1., ndex_hV, ndim_hV )      ! surface

         ! Define the W grid FILE ( nid_W )

         CALL dia_nam( clhstnam, nn_write, 'grid_W' )                   ! filename
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam
         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           &  ! Horizontal grid: glamt and gphit
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, rn_Dt, nh_W, nid_W, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_W, "depthw", "Vertical W levels",      &  ! Vertical grid: gdepw
            &          "m", ipk, gdepw_1d, nz_W, "down" )

         IF( ln_abl ) THEN
         ! Define the ABL grid FILE ( nid_A )
            CALL dia_nam( clhstnam, nn_write, 'grid_ABL' )
            IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam    ! filename
            CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           &  ! Horizontal grid: glamt and gphit
               &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
               &          nit000-1, zjulian, rn_Dt, nh_A, nid_A, domain_id=nidom, snc4chunks=snc4set )
            CALL histvert( nid_A, "ght_abl", "Vertical T levels",      &  ! Vertical grid: gdept
               &           "m", ipka, ght_abl(2:jpka), nz_A, "up" )
            !                                                            ! Index of ocean points
			ALLOCATE( zw3d_abl(jpi,jpj,ipka) )
			zw3d_abl(:,:,:) = 1._wp
			CALL wheneq( jpi*jpj*ipka, zw3d_abl, 1, 1., ndex_A , ndim_A  )      ! volume
            CALL wheneq( jpi*jpj     , zw3d_abl, 1, 1., ndex_hA, ndim_hA )      ! surface
			DEALLOCATE(zw3d_abl)
         ENDIF

         ! Declare all the output fields as NETCDF variables

         !                                                                                      !!! nid_T : 3D
         CALL histdef( nid_T, "votemper", "Temperature"                        , "C"      ,   &  ! tn
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         CALL histdef( nid_T, "vosaline", "Salinity"                           , "PSU"    ,   &  ! sn
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         IF(  .NOT.ln_linssh  ) THEN
            CALL histdef( nid_T, "vovvle3t", "Level thickness"                    , "m"      ,&  ! ze3t(:,:,:,Kmm)
            &             jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
            CALL histdef( nid_T, "vovvldep", "T point depth"                      , "m"      ,&  ! ze3t(:,:,:,Kmm)
            &             jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
            CALL histdef( nid_T, "vovvldef", "Squared level deformation"          , "%^2"    ,&  ! ze3t(:,:,:,Kmm)
            &             jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         ENDIF
         !                                                                                      !!! nid_T : 2D
         CALL histdef( nid_T, "sosstsst", "Sea Surface temperature"            , "C"      ,   &  ! sst
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosaline", "Sea Surface Salinity"               , "PSU"    ,   &  ! sss
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sossheig", "Sea Surface Height"                 , "m"      ,   &  ! ssh
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowaflup", "Net Upward Water Flux"              , "Kg/m2/s",   &  ! (emp-rnf)
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sorunoff", "River runoffs"                      , "Kg/m2/s",   &  ! runoffs
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosfldow", "downward salt flux"                 , "PSU/m2/s",  &  ! sfx
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         IF(  ln_linssh  ) THEN
            CALL histdef( nid_T, "sosst_cd", "Concentration/Dilution term on temperature"     &  ! emp * ts(:,:,1,jp_tem,Kmm)
            &                                                                  , "KgC/m2/s",  &  ! sosst_cd
            &             jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "sosss_cd", "Concentration/Dilution term on salinity"        &  ! emp * ts(:,:,1,jp_sal,Kmm)
            &                                                                  , "KgPSU/m2/s",&  ! sosss_cd
            &             jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         ENDIF
         CALL histdef( nid_T, "sohefldo", "Net Downward Heat Flux"             , "W/m2"   ,   &  ! qns + qsr
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "soshfldo", "Shortwave Radiation"                , "W/m2"   ,   &  ! qsr
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "somixhgt", "Turbocline Depth"                   , "m"      ,   &  ! hmld
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "somxl010", "Mixed Layer Depth 0.01"             , "m"      ,   &  ! hmlp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "soicecov", "Ice fraction"                       , "[0,1]"  ,   &  ! fr_i
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowindsp", "wind speed at 10m"                  , "m/s"    ,   &  ! wndm
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         !
         IF( ln_abl ) THEN
            CALL histdef( nid_A, "t_abl", "Potential Temperature"     , "K"        ,       &  ! t_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout )
            CALL histdef( nid_A, "q_abl", "Humidity"                  , "kg/kg"    ,       &  ! q_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout )
            CALL histdef( nid_A, "u_abl", "Atmospheric U-wind   "     , "m/s"        ,     &  ! u_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout )
            CALL histdef( nid_A, "v_abl", "Atmospheric V-wind   "     , "m/s"    ,         &  ! v_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout )
            CALL histdef( nid_A, "tke_abl", "Atmospheric TKE   "     , "m2/s2"    ,        &  ! tke_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout )
            CALL histdef( nid_A, "avm_abl", "Atmospheric turbulent viscosity", "m2/s"   ,  &  ! avm_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout )
            CALL histdef( nid_A, "avt_abl", "Atmospheric turbulent diffusivity", "m2/s2",  &  ! avt_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout )
            CALL histdef( nid_A, "pblh", "Atmospheric boundary layer height "  , "m",      &  ! pblh
               &          jpi, jpj, nh_A,  1  , 1, 1   , -99 , 32, clop, zsto, zout )
#if defined key_si3
            CALL histdef( nid_A, "oce_frac", "Fraction of open ocean"  , " ",      &  ! ato_i
               &          jpi, jpj, nh_A,  1  , 1, 1   , -99 , 32, clop, zsto, zout )
#endif
            CALL histend( nid_A, snc4chunks=snc4set )
         ENDIF
         !
         IF( ln_icebergs ) THEN
            CALL histdef( nid_T, "calving"             , "calving mass input"                       , "kg/s"   , &
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "calving_heat"        , "calving heat flux"                        , "XXXX"   , &
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "berg_floating_melt"  , "Melt rate of icebergs + bits"             , "kg/m2/s", &
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "berg_stored_ice"     , "Accumulated ice mass by class"            , "kg"     , &
               &          jpi, jpj, nh_T, nclasses  , 1, nclasses  , nb_T , 32, clop, zsto, zout )
            IF( ln_bergdia ) THEN
               CALL histdef( nid_T, "berg_melt"           , "Melt rate of icebergs"                    , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_buoy_melt"      , "Buoyancy component of iceberg melt rate"  , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_eros_melt"      , "Erosion component of iceberg melt rate"   , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_conv_melt"      , "Convective component of iceberg melt rate", "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_virtual_area"   , "Virtual coverage by icebergs"             , "m2"     , &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "bits_src"           , "Mass source of bergy bits"                , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "bits_melt"          , "Melt rate of bergy bits"                  , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "bits_mass"          , "Bergy bit density field"                  , "kg/m2"  , &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_mass"           , "Iceberg density field"                    , "kg/m2"  , &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_real_calving"   , "Calving into iceberg class"               , "kg/s"   , &
                  &          jpi, jpj, nh_T, nclasses  , 1, nclasses  , nb_T , 32, clop, zsto, zout )
            ENDIF
         ENDIF

         IF( ln_ssr ) THEN
            CALL histdef( nid_T, "sohefldp", "Surface Heat Flux: Damping"         , "W/m2"   ,   &  ! qrp
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "sowafldp", "Surface Water Flux: Damping"        , "Kg/m2/s",   &  ! erp
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "sosafldp", "Surface salt flux: damping"         , "Kg/m2/s",   &  ! erp * sn
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         ENDIF

         clmx ="l_max(only(x))"    ! max index on a period
!         CALL histdef( nid_T, "sobowlin", "Bowl Index"                         , "W-point",   &  ! bowl INDEX
!            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clmx, zsto, zout )
#if defined key_diahth
         CALL histdef( nid_T, "sothedep", "Thermocline Depth"                  , "m"      ,   & ! hth
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "so20chgt", "Depth of 20C isotherm"              , "m"      ,   & ! hd20
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "so28chgt", "Depth of 28C isotherm"              , "m"      ,   & ! hd28
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sohtc300", "Heat content 300 m"                 , "J/m2"   ,   & ! htc3
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif

         CALL histend( nid_T, snc4chunks=snc4set )

         !                                                                                      !!! nid_U : 3D
         CALL histdef( nid_U, "vozocrtx", "Zonal Current"                      , "m/s"    ,   &  ! uu(:,:,:,Kmm)
            &          jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout )
         IF( ln_wave .AND. ln_sdw) THEN
            CALL histdef( nid_U, "sdzocrtx", "Stokes Drift Zonal Current"         , "m/s"    ,   &  ! usd
               &          jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout )
         ENDIF
         !                                                                                      !!! nid_U : 2D
         CALL histdef( nid_U, "sozotaux", "Wind Stress along i-axis"           , "N/m2"   ,   &  ! utau
            &          jpi, jpj, nh_U, 1  , 1, 1  , - 99, 32, clop, zsto, zout )

         CALL histend( nid_U, snc4chunks=snc4set )

         !                                                                                      !!! nid_V : 3D
         CALL histdef( nid_V, "vomecrty", "Meridional Current"                 , "m/s"    ,   &  ! vv(:,:,:,Kmm)
            &          jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout )
         IF( ln_wave .AND. ln_sdw) THEN
            CALL histdef( nid_V, "sdmecrty", "Stokes Drift Meridional Current"    , "m/s"    ,   &  ! vsd
               &          jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout )
         ENDIF
         !                                                                                      !!! nid_V : 2D
         CALL histdef( nid_V, "sometauy", "Wind Stress along j-axis"           , "N/m2"   ,   &  ! vtau
            &          jpi, jpj, nh_V, 1  , 1, 1  , - 99, 32, clop, zsto, zout )

         CALL histend( nid_V, snc4chunks=snc4set )

         !                                                                                      !!! nid_W : 3D
         CALL histdef( nid_W, "vovecrtz", "Vertical Velocity"                  , "m/s"    ,   &  ! ww
            &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         CALL histdef( nid_W, "votkeavt", "Vertical Eddy Diffusivity"          , "m2/s"   ,   &  ! avt
            &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         CALL histdef( nid_W, "votkeavm", "Vertical Eddy Viscosity"             , "m2/s"  ,   &  ! avm
            &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )

         IF( ln_zdfddm ) THEN
            CALL histdef( nid_W,"voddmavs","Salt Vertical Eddy Diffusivity"    , "m2/s"   ,   &  ! avs
               &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         ENDIF

         IF( ln_wave .AND. ln_sdw) THEN
            CALL histdef( nid_W, "sdvecrtz", "Stokes Drift Vertical Current"   , "m/s"    ,   &  ! wsd
               &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         ENDIF
         !                                                                                      !!! nid_W : 2D
         CALL histend( nid_W, snc4chunks=snc4set )

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization'
         IF(ll_print) CALL FLUSH(numout )

      ENDIF

      ! 2. Start writing data
      ! ---------------------

      ! ndex(1) est utilise ssi l'avant dernier argument est different de
      ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
      ! donne le nombre d'elements, et ndex la liste des indices a sortir

      IF( lwp .AND. MOD( itmod, nn_write ) == 0 ) THEN
         WRITE(numout,*) 'dia_wri : write model outputs in NetCDF files at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~ '
      ENDIF
!!st16
      IF( .NOT.ln_linssh ) THEN
         CALL histwrite( nid_T, "votemper", it, ts(:,:,:,jp_tem,Kmm) * ze3t(:,:,:) , ndim_T , ndex_T  )   ! heat content
         CALL histwrite( nid_T, "vosaline", it, ts(:,:,:,jp_sal,Kmm) * ze3t(:,:,:) , ndim_T , ndex_T  )   ! salt content
         CALL histwrite( nid_T, "sosstsst", it, ts(:,:,1,jp_tem,Kmm) * ze3t(:,:,1) , ndim_hT, ndex_hT )   ! sea surface heat content
         CALL histwrite( nid_T, "sosaline", it, ts(:,:,1,jp_sal,Kmm) * ze3t(:,:,1) , ndim_hT, ndex_hT )   ! sea surface salinity content
!!st16
      ELSE
         CALL histwrite( nid_T, "votemper", it, ts(:,:,:,jp_tem,Kmm) , ndim_T , ndex_T  )   ! temperature
         CALL histwrite( nid_T, "vosaline", it, ts(:,:,:,jp_sal,Kmm) , ndim_T , ndex_T  )   ! salinity
         CALL histwrite( nid_T, "sosstsst", it, ts(:,:,1,jp_tem,Kmm) , ndim_hT, ndex_hT )   ! sea surface temperature
         CALL histwrite( nid_T, "sosaline", it, ts(:,:,1,jp_sal,Kmm) , ndim_hT, ndex_hT )   ! sea surface salinity
      ENDIF
      IF( .NOT.ln_linssh ) THEN
!!st17 if ! defined key_qco
         zw3d(:,:,:) = ( ( ze3t(:,:,:) - e3t_0(:,:,:) ) / e3t_0(:,:,:) * 100 * tmask(:,:,:) ) ** 2
         CALL histwrite( nid_T, "vovvle3t", it, ze3t (:,:,:)     , ndim_T , ndex_T  )   ! level thickness
         CALL histwrite( nid_T, "vovvldep", it, zgdept , ndim_T , ndex_T  )   ! t-point depth
!!st17
         CALL histwrite( nid_T, "vovvldef", it, zw3d             , ndim_T , ndex_T  )   ! level thickness deformation
      ENDIF
      CALL histwrite( nid_T, "sossheig", it, ssh(:,:,Kmm)          , ndim_hT, ndex_hT )   ! sea surface height
      CALL histwrite( nid_T, "sowaflup", it, ( emp-rnf )   , ndim_hT, ndex_hT )   ! upward water flux
      CALL histwrite( nid_T, "sorunoff", it, rnf           , ndim_hT, ndex_hT )   ! river runoffs
      CALL histwrite( nid_T, "sosfldow", it, sfx           , ndim_hT, ndex_hT )   ! downward salt flux
                                                                                  ! (includes virtual salt flux beneath ice
                                                                                  ! in linear free surface case)
      IF( ln_linssh ) THEN
         zw2d(:,:) = emp (:,:) * ts(:,:,1,jp_tem,Kmm)
         CALL histwrite( nid_T, "sosst_cd", it, zw2d, ndim_hT, ndex_hT )          ! c/d term on sst
         zw2d(:,:) = emp (:,:) * ts(:,:,1,jp_sal,Kmm)
         CALL histwrite( nid_T, "sosss_cd", it, zw2d, ndim_hT, ndex_hT )          ! c/d term on sss
      ENDIF
      CALL histwrite( nid_T, "sohefldo", it, qns + qsr     , ndim_hT, ndex_hT )   ! total heat flux
      CALL histwrite( nid_T, "soshfldo", it, qsr           , ndim_hT, ndex_hT )   ! solar heat flux
      CALL histwrite( nid_T, "somixhgt", it, hmld          , ndim_hT, ndex_hT )   ! turbocline depth
      CALL histwrite( nid_T, "somxl010", it, hmlp          , ndim_hT, ndex_hT )   ! mixed layer depth
      CALL histwrite( nid_T, "soicecov", it, fr_i          , ndim_hT, ndex_hT )   ! ice fraction
      CALL histwrite( nid_T, "sowindsp", it, wndm          , ndim_hT, ndex_hT )   ! wind speed
      !
      IF( ln_abl ) THEN
         ALLOCATE( zw3d_abl(jpi,jpj,jpka) )
         IF( ln_mskland )   THEN
            DO jk=1,jpka
               zw3d_abl(:,:,jk) = tmask(:,:,1)
            END DO
         ELSE
            zw3d_abl(:,:,:) = 1._wp
         ENDIF
         CALL histwrite( nid_A,  "pblh"   , it, pblh(:,:)                  *zw3d_abl(:,:,1     ), ndim_hA, ndex_hA )   ! pblh
         CALL histwrite( nid_A,  "u_abl"  , it, u_abl   (:,:,2:jpka,nt_n  )*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! u_abl
         CALL histwrite( nid_A,  "v_abl"  , it, v_abl   (:,:,2:jpka,nt_n  )*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! v_abl
         CALL histwrite( nid_A,  "t_abl"  , it, tq_abl  (:,:,2:jpka,nt_n,1)*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! t_abl
         CALL histwrite( nid_A,  "q_abl"  , it, tq_abl  (:,:,2:jpka,nt_n,2)*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! q_abl
         CALL histwrite( nid_A,  "tke_abl", it, tke_abl (:,:,2:jpka,nt_n  )*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! tke_abl
         CALL histwrite( nid_A,  "avm_abl", it, avm_abl (:,:,2:jpka       )*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! avm_abl
         CALL histwrite( nid_A,  "avt_abl", it, avt_abl (:,:,2:jpka       )*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! avt_abl
#if defined key_si3
         CALL histwrite( nid_A,  "oce_frac"   , it, ato_i(:,:)                                  , ndim_hA, ndex_hA )   ! ato_i
#endif
         DEALLOCATE(zw3d_abl)
      ENDIF
      !
      IF( ln_icebergs ) THEN
         !
         CALL histwrite( nid_T, "calving"             , it, berg_grid%calving      , ndim_hT, ndex_hT )
         CALL histwrite( nid_T, "calving_heat"        , it, berg_grid%calving_hflx , ndim_hT, ndex_hT )
         CALL histwrite( nid_T, "berg_floating_melt"  , it, berg_grid%floating_melt, ndim_hT, ndex_hT )
         !
         CALL histwrite( nid_T, "berg_stored_ice"     , it, berg_grid%stored_ice   , ndim_bT, ndex_bT )
         !
         IF( ln_bergdia ) THEN
            CALL histwrite( nid_T, "berg_melt"           , it, berg_melt        , ndim_hT, ndex_hT   )
            CALL histwrite( nid_T, "berg_buoy_melt"      , it, buoy_melt        , ndim_hT, ndex_hT   )
            CALL histwrite( nid_T, "berg_eros_melt"      , it, eros_melt        , ndim_hT, ndex_hT   )
            CALL histwrite( nid_T, "berg_conv_melt"      , it, conv_melt        , ndim_hT, ndex_hT   )
            CALL histwrite( nid_T, "berg_virtual_area"   , it, virtual_area     , ndim_hT, ndex_hT   )
            CALL histwrite( nid_T, "bits_src"            , it, bits_src         , ndim_hT, ndex_hT   )
            CALL histwrite( nid_T, "bits_melt"           , it, bits_melt        , ndim_hT, ndex_hT   )
            CALL histwrite( nid_T, "bits_mass"           , it, bits_mass        , ndim_hT, ndex_hT   )
            CALL histwrite( nid_T, "berg_mass"           , it, berg_mass        , ndim_hT, ndex_hT   )
            !
            CALL histwrite( nid_T, "berg_real_calving"   , it, real_calving     , ndim_bT, ndex_bT   )
         ENDIF
      ENDIF

      IF( ln_ssr ) THEN
         CALL histwrite( nid_T, "sohefldp", it, qrp           , ndim_hT, ndex_hT )   ! heat flux damping
         CALL histwrite( nid_T, "sowafldp", it, erp           , ndim_hT, ndex_hT )   ! freshwater flux damping
         zw2d(:,:) = erp(:,:) * ts(:,:,1,jp_sal,Kmm) * tmask(:,:,1)
         CALL histwrite( nid_T, "sosafldp", it, zw2d          , ndim_hT, ndex_hT )   ! salt flux damping
      ENDIF
!      zw2d(:,:) = FLOAT( nmln(:,:) ) * tmask(:,:,1)
!      CALL histwrite( nid_T, "sobowlin", it, zw2d          , ndim_hT, ndex_hT )   ! ???

#if defined key_diahth
      CALL histwrite( nid_T, "sothedep", it, hth           , ndim_hT, ndex_hT )   ! depth of the thermocline
      CALL histwrite( nid_T, "so20chgt", it, hd20          , ndim_hT, ndex_hT )   ! depth of the 20 isotherm
      CALL histwrite( nid_T, "so28chgt", it, hd28          , ndim_hT, ndex_hT )   ! depth of the 28 isotherm
      CALL histwrite( nid_T, "sohtc300", it, htc3          , ndim_hT, ndex_hT )   ! first 300m heaat content
#endif

      CALL histwrite( nid_U, "vozocrtx", it, uu(:,:,:,Kmm)            , ndim_U , ndex_U )    ! i-current
      CALL histwrite( nid_U, "sozotaux", it, utau          , ndim_hU, ndex_hU )   ! i-wind stress

      CALL histwrite( nid_V, "vomecrty", it, vv(:,:,:,Kmm)            , ndim_V , ndex_V  )   ! j-current
      CALL histwrite( nid_V, "sometauy", it, vtau          , ndim_hV, ndex_hV )   ! j-wind stress

      IF( ln_zad_Aimp ) THEN
         CALL histwrite( nid_W, "vovecrtz", it, ww + wi     , ndim_T, ndex_T )    ! vert. current
      ELSE
         CALL histwrite( nid_W, "vovecrtz", it, ww          , ndim_T, ndex_T )    ! vert. current
      ENDIF
      CALL histwrite( nid_W, "votkeavt", it, avt            , ndim_T, ndex_T )    ! T vert. eddy diff. coef.
      CALL histwrite( nid_W, "votkeavm", it, avm            , ndim_T, ndex_T )    ! T vert. eddy visc. coef.
      IF( ln_zdfddm ) THEN
         CALL histwrite( nid_W, "voddmavs", it, avs         , ndim_T, ndex_T )    ! S vert. eddy diff. coef.
      ENDIF

      IF( ln_wave .AND. ln_sdw ) THEN
         CALL histwrite( nid_U, "sdzocrtx", it, usd         , ndim_U , ndex_U )    ! i-StokesDrift-current
         CALL histwrite( nid_V, "sdmecrty", it, vsd         , ndim_V , ndex_V )    ! j-StokesDrift-current
         CALL histwrite( nid_W, "sdvecrtz", it, wsd         , ndim_T , ndex_T )    ! StokesDrift vert. current
      ENDIF

      ! 3. Close all files
      ! ---------------------------------------
      IF( kt == nitend ) THEN
         CALL histclo( nid_T )
         CALL histclo( nid_U )
         CALL histclo( nid_V )
         CALL histclo( nid_W )
         IF(ln_abl) CALL histclo( nid_A )
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('dia_wri')
      !
   END SUBROUTINE dia_wri
#endif

   SUBROUTINE dia_wri_state( Kmm, cdfile_name )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_state  ***
      !!
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains
      !!      the instantaneous ocean state and forcing fields.
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! ** Method  :   NetCDF files using ioipsl
      !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
      !!      File 'output.abort.nc' is created in case of abnormal job end
      !!----------------------------------------------------------------------
      INTEGER           , INTENT( in ) ::   Kmm              ! time level index
      CHARACTER (len=* ), INTENT( in ) ::   cdfile_name      ! name of the file created
      !!
      INTEGER :: inum, jk
!!st18  TBR
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ze3t, zgdept      ! 3D workspace !!st patch to use substitution
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
      IF(lwp) WRITE(numout,*) '                and named :', cdfile_name, '...nc'
!!st19 TBR
      DO jk = 1, jpk
         ze3t(:,:,jk) =  e3t(:,:,jk,Kmm)
         zgdept(:,:,jk) =  gdept(:,:,jk,Kmm)
      END DO
!!st19
#if defined key_si3
     CALL iom_open( TRIM(cdfile_name), inum, ldwrt = .TRUE., kdlev = jpl )
#else
     CALL iom_open( TRIM(cdfile_name), inum, ldwrt = .TRUE. )
#endif

      CALL iom_rstput( 0, 0, inum, 'votemper', ts(:,:,:,jp_tem,Kmm) )    ! now temperature
      CALL iom_rstput( 0, 0, inum, 'vosaline', ts(:,:,:,jp_sal,Kmm) )    ! now salinity
      CALL iom_rstput( 0, 0, inum, 'sossheig', ssh(:,:,Kmm)              )    ! sea surface height
      CALL iom_rstput( 0, 0, inum, 'vozocrtx', uu(:,:,:,Kmm)                )    ! now i-velocity
      CALL iom_rstput( 0, 0, inum, 'vomecrty', vv(:,:,:,Kmm)                )    ! now j-velocity
      IF( ln_zad_Aimp ) THEN
         CALL iom_rstput( 0, 0, inum, 'vovecrtz', ww + wi        )    ! now k-velocity
      ELSE
         CALL iom_rstput( 0, 0, inum, 'vovecrtz', ww             )    ! now k-velocity
      ENDIF
      CALL iom_rstput( 0, 0, inum, 'risfdep', risfdep            )    ! now k-velocity
      CALL iom_rstput( 0, 0, inum, 'ht'     , ht(:,:)            )    ! now water column height

      IF ( ln_isf ) THEN
         IF (ln_isfcav_mlt) THEN
            CALL iom_rstput( 0, 0, inum, 'fwfisf_cav', fwfisf_cav          )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'rhisf_cav_tbl', rhisf_tbl_cav    )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'rfrac_cav_tbl', rfrac_tbl_cav    )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'misfkb_cav', REAL(misfkb_cav,wp) )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'misfkt_cav', REAL(misfkt_cav,wp) )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'mskisf_cav', REAL(mskisf_cav,wp), ktype = jp_i1 )
         END IF
         IF (ln_isfpar_mlt) THEN
            CALL iom_rstput( 0, 0, inum, 'isfmsk_par', REAL(mskisf_par,wp) )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'fwfisf_par', fwfisf_par          )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'rhisf_par_tbl', rhisf_tbl_par    )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'rfrac_par_tbl', rfrac_tbl_par    )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'misfkb_par', REAL(misfkb_par,wp) )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'misfkt_par', REAL(misfkt_par,wp) )    ! now k-velocity
            CALL iom_rstput( 0, 0, inum, 'mskisf_par', REAL(mskisf_par,wp), ktype = jp_i1 )
         END IF
      END IF
      !
      IF( ALLOCATED(ahtu) ) THEN
         CALL iom_rstput( 0, 0, inum,  'ahtu', ahtu              )    ! aht at u-point
         CALL iom_rstput( 0, 0, inum,  'ahtv', ahtv              )    ! aht at v-point
      ENDIF
      IF( ALLOCATED(ahmt) ) THEN
         CALL iom_rstput( 0, 0, inum,  'ahmt', ahmt              )    ! ahmt at u-point
         CALL iom_rstput( 0, 0, inum,  'ahmf', ahmf              )    ! ahmf at v-point
      ENDIF
      CALL iom_rstput( 0, 0, inum, 'sowaflup', emp - rnf         )    ! freshwater budget
      CALL iom_rstput( 0, 0, inum, 'sohefldo', qsr + qns         )    ! total heat flux
      CALL iom_rstput( 0, 0, inum, 'soshfldo', qsr               )    ! solar heat flux
      CALL iom_rstput( 0, 0, inum, 'soicecov', fr_i              )    ! ice fraction
      CALL iom_rstput( 0, 0, inum, 'sozotaux', utau              )    ! i-wind stress
      CALL iom_rstput( 0, 0, inum, 'sometauy', vtau              )    ! j-wind stress
!!st20 TBR
      IF(  .NOT.ln_linssh  ) THEN
         CALL iom_rstput( 0, 0, inum, 'vovvldep', zgdept        )    !  T-cell depth
         CALL iom_rstput( 0, 0, inum, 'vovvle3t', ze3t          )    !  T-cell thickness
      END IF
      IF( ln_wave .AND. ln_sdw ) THEN
         CALL iom_rstput( 0, 0, inum, 'sdzocrtx', usd            )    ! now StokesDrift i-velocity
         CALL iom_rstput( 0, 0, inum, 'sdmecrty', vsd            )    ! now StokesDrift j-velocity
         CALL iom_rstput( 0, 0, inum, 'sdvecrtz', wsd            )    ! now StokesDrift k-velocity
      ENDIF
      IF ( ln_abl ) THEN
         CALL iom_rstput ( 0, 0, inum, "uz1_abl",   u_abl(:,:,2,nt_a  ) )   ! now first level i-wind
         CALL iom_rstput ( 0, 0, inum, "vz1_abl",   v_abl(:,:,2,nt_a  ) )   ! now first level j-wind
         CALL iom_rstput ( 0, 0, inum, "tz1_abl",  tq_abl(:,:,2,nt_a,1) )   ! now first level temperature
         CALL iom_rstput ( 0, 0, inum, "qz1_abl",  tq_abl(:,:,2,nt_a,2) )   ! now first level humidity
      ENDIF

#if defined key_si3
      IF( nn_ice == 2 ) THEN   ! condition needed in case agrif + ice-model but no-ice in child grid
         CALL ice_wri_state( inum )
      ENDIF
#endif
      !
      CALL iom_close( inum )
      !
   END SUBROUTINE dia_wri_state

   !!======================================================================
END MODULE diawri
