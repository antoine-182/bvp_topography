MODULE step
   !!======================================================================
   !!                       ***  MODULE step  ***
   !! Time-stepping   : manager of the shallow water equation time stepping
   !!======================================================================
   !! History :  NEMO !  2020-03  (A. Nasser, G. Madec)  Original code from  4.0.2
   !!----------------------------------------------------------------------
#if defined key_qco
   !!----------------------------------------------------------------------
   !!   'key_qco'      EMPTY MODULE      Quasi-Eulerian vertical coordonate
   !!----------------------------------------------------------------------
#else
   !!----------------------------------------------------------------------
   !!   stp             : Shallow Water time-stepping
   !!----------------------------------------------------------------------
   USE step_oce         ! time stepping definition modules
   USE phycst           ! physical constants
   USE usrdef_nam
   USE lib_mpp        ! MPP library
   USE iom              ! xIOs server

   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp   ! called by nemogcm.F90

   !!----------------------------------------------------------------------
   !! time level indices
   !!----------------------------------------------------------------------
   INTEGER, PUBLIC ::   Nbb, Nnn, Naa, Nrhs   !! used by nemo_init
   ! INTEGER, PUBLIC, PARAMETER ::   np_rfr = 0   ! ENS scheme
   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: step.F90 13154 2020-06-24 13:31:32Z gm $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_agrif
   RECURSIVE SUBROUTINE stp( )
      INTEGER             ::   kstp   ! ocean time-step index
#else
   SUBROUTINE stp( kstp )
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
#endif
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp  ***
      !!
      !! ** Purpose : - Time stepping of shallow water (SHW) (momentum and ssh eqs.)
      !!
      !! ** Method  : -1- Update forcings
      !!              -2- Update the ssh at Naa
      !!              -3- Compute the momentum trends (Nrhs)
      !!              -4- Update the horizontal velocity
      !!              -5- Apply Asselin time filter to uu,vv,ssh
      !!              -6- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk   ! dummy loop indice
      INTEGER ::   indic        ! error indicator if < 0
!!gm kcall can be removed, I guess
      INTEGER ::   kcall        ! optional integer argument (dom_vvl_sf_nxt)
      REAL(wp)::   z1_2rho0     ! local scalars

      REAL(wp) ::   zue3a, zue3n, zue3b    ! local scalars
      REAL(wp) ::   zve3a, zve3n, zve3b    !   -      -
      REAL(wp) ::   ze3t_tf, ze3u_tf, ze3v_tf, zua, zva
      REAL(wp), DIMENSION(jpi,jpj)    ::   zssh, zFu, zFv, zuu, zvv
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ze3u, ze3v, z3du, z3dv
      !! ---------------------------------------------------------------------
#if defined key_agrif
      kstp = nit000 + Agrif_Nb_Step()
      Kbb_a = Nbb; Kmm_a = Nnn; Krhs_a = Nrhs   ! agrif_oce module copies of time level indices
      IF( lk_agrif_debug ) THEN
         IF( Agrif_Root() .and. lwp)   WRITE(*,*) '---'
         IF(lwp)   WRITE(*,*) 'Grid Number', Agrif_Fixed(),' time step ', kstp, 'int tstep', Agrif_NbStepint()
      ENDIF
      IF( kstp == nit000 + 1 )   lk_agrif_fstep = .FALSE.
# if defined key_iomput
      IF( Agrif_Nbstepint() == 0 )   CALL iom_swap( cxios_context )
# endif
#endif
      !
      IF( ln_timing )   CALL timing_start('stp')
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! model timestep
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      IF( l_1st_euler ) THEN     ! start or restart with Euler 1st time-step
         rDt   =  rn_Dt
         r1_Dt = 1._wp / rDt
      ENDIF

      IF ( kstp == nit000 )   ww(:,:,:) = 0._wp   ! initialize vertical velocity one for all to zero

      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! update I/O and calendar
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                             indic = 0                ! reset to no error condition

      IF( kstp == nit000 ) THEN                       ! initialize IOM context (must be done after nemo_init for AGRIF+XIOS+OASIS)
                             CALL iom_init( cxios_context, ld_closedef=.FALSE. )   ! for model grid (including passible AGRIF zoom)
         IF( lk_diamlr   )   CALL dia_mlr_iom_init    ! with additional setup for multiple-linear-regression analysis
                             CALL iom_init_closedef
         IF( ln_crs      )   CALL iom_init( TRIM(cxios_context)//"_crs" )  ! for coarse grid
      ENDIF
      IF( kstp /= nit000 )   CALL day( kstp )         ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp - nit000 + 1,      cxios_context          )   ! tell IOM we are at time step kstp
      IF( ln_crs         )   CALL iom_setkt( kstp - nit000 + 1, TRIM(cxios_context)//"_crs" )   ! tell IOM we are at time step kstp

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update external forcing (tides, open boundaries, ice shelf interaction and surface boundary condition (including sea-ice)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_tide    )   CALL tide_update( kstp )                     ! update tide potential
      IF( ln_apr_dyn )   CALL sbc_apr ( kstp )                        ! atmospheric pressure (NB: call before bdy_dta which needs ssh_ib)
      IF( ln_bdy     )   CALL bdy_dta ( kstp, Nnn )                   ! update dynamic & tracer data at open boundaries
                         CALL sbc     ( kstp, Nbb, Nnn )              ! Sea Boundary Condition

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Ocean physics update
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  LATERAL  PHYSICS
      !                                                                        ! eddy diffusivity coeff.
      IF( l_ldfdyn_time )   CALL ldf_dyn( kstp, Nbb )                          ! eddy viscosity coeff.

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  Ocean dynamics : hdiv, ssh, e3, u, v, w
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                            CALL ssh_nxt       ( kstp, Nbb, Nnn, ssh, Naa )    ! after ssh (includes call to div_hor)
                         uu(:,:,:,Nrhs) = 0._wp            ! set dynamics trends to zero
                         vv(:,:,:,Nrhs) = 0._wp

      !!an compute the after scale factors from the ssh after
      !!an nothing is done on e3f because no need of a trend for it
      IF( .NOT.ln_linssh )  CALL dom_vvl_sf_nxt( kstp, Nbb, Nnn,      Naa )    ! after vertical scale factors

      IF( ln_bdy     )      CALL bdy_dyn3d_dmp ( kstp, Nbb,      uu, vv, Nrhs )  ! bdy damping trends

#if defined key_agrif
      IF(.NOT. Agrif_Root())  &
               &            CALL Agrif_Sponge_dyn        ! momentum sponge
#endif

!!an - calcul du gradient de pression horizontal (explicit)
      zssh(:,:) = ssh(:,:,Nnn)
# if defined key_bvp
      IF ( kstp == nit000 ) WRITE(numout,*) 'grad(ssh/rpo) in the RHS (step) !'
      zssh(:,:) = zssh(:,:) * r1_rpo(:,:,1)

      !!!! Useful to modelise bathymetry
      ! IF ( kstp == nit000 ) WRITE(numout,*) 'grad(ssh) penalised in the RHS (step) !'
#endif
      !
      DO_3D_00_00( 1, jpkm1 )
         uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - grav * ( zssh(ji+1,jj) - zssh(ji,jj) ) * r1_e1u(ji,jj)
         vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - grav * ( zssh(ji,jj+1) - zssh(ji,jj) ) * r1_e2v(ji,jj)
      END_3D
      !

!!an
! Tentative - Echec
      ! IF ( kstp == nit000 ) WRITE(numout,*) 'gradient rigolo'
      ! DO_2D_00_00
      !    zuu(ji,jj) = ( zssh(ji+1,jj  ) - zssh(ji,jj) ) * r1_e1u(ji,jj)
      !    zvv(ji,jj) = ( zssh(ji  ,jj+1) - zssh(ji,jj) ) * r1_e2v(ji,jj)
      ! END_2D
      ! CALL lbc_lnk_multi( 'stp', zuu, 'U', 1. , zvv, 'V', 1. )
      ! DO_2D_00_00
      !    zFu(ji,jj) = 0.25_wp * ( zuu(ji  ,jj+1) + zuu(ji,jj) )
      !    zFv(ji,jj) = 0.25_wp * ( zvv(ji+1,jj  ) + zvv(ji,jj) )
      ! END_2D
      ! CALL lbc_lnk_multi( 'stp', zFu, 'F', 1. , zFu, 'F', 1. )
      ! DO_3D_00_00( 1, jpkm1 )
      !     ! à l'intérieur du bassin
      !    IF (fmask(ji,jj,1)*fmask(ji,jj+1,1) == 1._wp) THEN
      !      uu(ji,jj+1,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - grav * ( zFu(ji,jj+1) + zFu(ji,jj))
      !    ELSE
      !      uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - grav * ( zssh(ji+1,jj) - zssh(ji,jj) ) * r1_e1u(ji,jj)
      !    ENDIF
      !    IF (fmask(ji,jj,1)*fmask(ji+1,jj,1) == 1._wp) THEN
      !      vv(ji+1,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - grav * ( zFv(ji+1,jj) + zFv(ji,jj))
      !   ELSE
      !     vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - grav * ( zssh(ji,jj+1) - zssh(ji,jj) ) * r1_e2v(ji,jj)
      !    ENDIF
      ! END_3D
!
!!an
!      IF( kstp == nit000 .AND. lwp ) THEN
!         WRITE(numout,*)
!         WRITE(numout,*) 'step.F90 : classic script used'
!         WRITE(numout,*) '~~~~~~~'
!         IF(       ln_dynvor_ens_adVO .OR. ln_dynvor_ens_adKE .OR. ln_dynvor_ens_adKEVO   &
!         &    .OR. ln_dynvor_ene_adVO .OR. ln_dynvor_ene_adKE .OR. ln_dynvor_ene_adKEVO   ) THEN
!            CALL ctl_stop('STOP','step : alternative direction asked but classis step'  )
!         ENDIF
!      ENDIF
!!an

# if defined key_bath
      DO_3D_00_00(1,jpkm1)
         uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) + grav * ( batht(ji+1,jj  ) - batht(ji,jj) ) * r1_e1u(ji,jj)
         vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) + grav * ( batht(ji  ,jj+1) - batht(ji,jj) ) * r1_e1v(ji,jj)
      END_3D
#endif

                         CALL dyn_adv( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! advection (VF or FF)	==> RHS

                         CALL dyn_vor( kstp, Nbb, Nnn      , uu, vv, Nrhs)   ! vorticity           	==> RHS

                         CALL dyn_ldf( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! lateral mixing

      ! add wind stress forcing and layer linear friction to the RHS
      ze3u(:,:,:) = e3u(:,:,:,Nnn)
      ze3v(:,:,:) = e3v(:,:,:,Nnn)
# if defined key_bvp
      IF ( kstp == nit000 ) WRITE(numout,*) 'Wind stress divided by H (not H_tilde)'
      !!!! if e3 not penalised, no need to divide h by rpo
      ze3u(:,:,:) = ze3u(:,:,:) * r1_rpou(:,:,:)
      ze3v(:,:,:) = ze3v(:,:,:) * r1_rpov(:,:,:)

      ! IF ( kstp == nit000 ) WRITE(numout,*) 'Wind stress divided by H_tilde (not H)'
      ! !!!! if e3 penalised, no need to multiply h by rpo
      ! ze3u(:,:,:) = ze3u(:,:,:) * rpou(:,:,:)
      ! ze3v(:,:,:) = ze3v(:,:,:) * rpov(:,:,:)
#endif
      z1_2rho0 = 0.5_wp * r1_rho0    ! == Wind stress
      !IF( kstp == 10 ) WRITE(*,*) utau_b(1,1), utau(1,1)
      DO_3D_00_00(1,jpkm1)
         uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) + z1_2rho0 * ( utau_b(ji,jj) + utau(ji,jj) ) / ze3u(ji,jj,jk)
         vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) + z1_2rho0 * ( vtau_b(ji,jj) + vtau(ji,jj) ) / ze3v(ji,jj,jk)
      END_3D

# if defined key_bvp
      IF ( kstp == nit000 ) write(numout,*) 'step : layer drag penalisation used nn_rfr=',nn_rfr
      SELECT CASE( nn_rfr )           ! == layer drag formulation
      CASE ( -2 )
        DO_3D_00_00(1,jpkm1)
           uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - rn_rfr * uu(ji,jj,jk,Nbb) * rpou(ji,jj,jk) * rpou(ji,jj,jk)
           vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - rn_rfr * vv(ji,jj,jk,Nbb) * rpov(ji,jj,jk) * rpov(ji,jj,jk)
        END_3D
      CASE ( -1 )
        DO_3D_00_00(1,jpkm1)
           uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - rn_rfr * uu(ji,jj,jk,Nbb) * rpou(ji,jj,jk)
           vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - rn_rfr * vv(ji,jj,jk,Nbb) * rpov(ji,jj,jk)
        END_3D
      CASE ( 0 )
        DO_3D_00_00(1,jpkm1)
           uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - rn_rfr * uu(ji,jj,jk,Nbb)
           vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - rn_rfr * vv(ji,jj,jk,Nbb)
        END_3D
      CASE ( 1 )
        DO_3D_00_00(1,jpkm1)
           uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - rn_rfr * uu(ji,jj,jk,Nbb) * r1_rpou(ji,jj,jk)
           vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - rn_rfr * vv(ji,jj,jk,Nbb) * r1_rpov(ji,jj,jk)
        END_3D
      CASE ( 2 )
        DO_3D_00_00(1,jpkm1)
           uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - rn_rfr * uu(ji,jj,jk,Nbb) * SQRT(r1_rpou(ji,jj,jk))
           vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - rn_rfr * vv(ji,jj,jk,Nbb) * SQRT(r1_rpov(ji,jj,jk))
        END_3D
      CASE ( 3 )
        DO_3D_00_00(1,jpkm1)
           uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - rn_rfr * uu(ji,jj,jk,Nbb) * ( r1_rpou(ji,jj,jk) ** ( REAL(nn_rfr, wp) ) )
           vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - rn_rfr * vv(ji,jj,jk,Nbb) * ( r1_rpov(ji,jj,jk) ** ( REAL(nn_rfr, wp) ) )
        END_3D
      CASE DEFAULT                                             ! error
         CALL ctl_stop('STOP','step: wrong value for nn_rfr'  )
      END SELECT
! # elif defined key_bath
!       z3du(:,:,:) = 0._wp ; z3dv(:,:,:) = 0._wp
!       WHERE(bathu(:,:) <= 4750._wp)
!         z3du(:,:,1) = 1._wp
!         z3du(:,:,2) = 1._wp
!       END WHERE
!       WHERE(bathv(:,:) <= 4750._wp)
!         z3dv(:,:,1) = 1._wp
!         z3dv(:,:,2) = 1._wp
!       END WHERE
!       DO_3D_00_00(1,jpkm1)
!          uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - ( rn_rfr + rn_fsp*z3du(ji,jj,jk) ) * uu(ji,jj,jk,Nbb)
!          vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - ( rn_rfr + rn_fsp*z3dv(ji,jj,jk) ) * vv(ji,jj,jk,Nbb)
!       END_3D
# else
      DO_3D_00_00(1,jpkm1)
        uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - rn_rfr * uu(ji,jj,jk,Nbb)
        vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - rn_rfr * vv(ji,jj,jk,Nbb)
      END_3D
#endif


# if defined key_bvp
      !  Add frictionnal term   - sigma * u
      !
      DO_3D_00_00(1,jpkm1)
         uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - bmpu(ji,jj,jk) * uu(ji,jj,jk,Nbb)
         vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - bmpv(ji,jj,jk) * vv(ji,jj,jk,Nbb)
      END_3D
#endif
      !

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Leap-Frog time splitting + Robert-Asselin time filter on u,v,e3
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!! what about  IF( .NOT.ln_linssh )  ?
!!an futur module dyn_nxt (a la place de dyn_atf)

!!an pour tester le vector form en linear
     ! IF( ln_dynadv_vec ) THEN      ! vector invariant form : applied on velocity
     IF( ln_dynadv_vec .OR. ln_dynadv_OFF ) THEN      ! vector invariant form : applied on velocity
         IF( l_1st_euler ) THEN           ! Euler time stepping (no Asselin filter)
            IF(lwp) write(numout,*) 'step : timestepping in Vect Form'
            DO_3D_00_00(1,jpkm1)
               uu(ji,jj,jk,Naa) = uu(ji,jj,jk,Nbb) + rDt * uu(ji,jj,jk,Nrhs) * umask(ji,jj,jk)
               vv(ji,jj,jk,Naa) = vv(ji,jj,jk,Nbb) + rDt * vv(ji,jj,jk,Nrhs) * vmask(ji,jj,jk)
             END_3D
         ELSE                             ! Leap Frog time stepping + Asselin filter
            DO_3D_11_11(1,jpkm1)
               zua = uu(ji,jj,jk,Nbb) + rDt * uu(ji,jj,jk,Nrhs) * umask(ji,jj,jk)
               zva = vv(ji,jj,jk,Nbb) + rDt * vv(ji,jj,jk,Nrhs) * vmask(ji,jj,jk)
               !                                                                  ! Asselin time filter on u,v (Nnn)
               uu(ji,jj,jk,Nnn) = uu(ji,jj,jk,Nnn) + rn_atfp * (uu(ji,jj,jk,Nbb) - 2._wp * uu(ji,jj,jk,Nnn) + zua)
               vv(ji,jj,jk,Nnn) = vv(ji,jj,jk,Nnn) + rn_atfp * (vv(ji,jj,jk,Nbb) - 2._wp * vv(ji,jj,jk,Nnn) + zva)
               !
               ze3u_tf = e3u(ji,jj,jk,Nnn) + rn_atfp * ( e3u(ji,jj,jk,Nbb) - 2._wp * e3u(ji,jj,jk,Nnn)  + e3u(ji,jj,jk,Naa) )
               ze3v_tf = e3v(ji,jj,jk,Nnn) + rn_atfp * ( e3v(ji,jj,jk,Nbb) - 2._wp * e3v(ji,jj,jk,Nnn)  + e3v(ji,jj,jk,Naa) )
               ze3t_tf = e3t(ji,jj,jk,Nnn) + rn_atfp * ( e3t(ji,jj,jk,Nbb) - 2._wp * e3t(ji,jj,jk,Nnn)  + e3t(ji,jj,jk,Naa) )
               !
               e3u(ji,jj,jk,Nnn) = ze3u_tf
               e3v(ji,jj,jk,Nnn) = ze3v_tf
               e3t(ji,jj,jk,Nnn) = ze3t_tf
               !
               uu(ji,jj,jk,Naa) = zua
               vv(ji,jj,jk,Naa) = zva
            END_3D
         ENDIF
         !
      ELSE                          ! flux form : applied on thickness weighted velocity
         IF( l_1st_euler ) THEN           ! Euler time stepping (no Asselin filter)
            IF(lwp) write(numout,*) 'step : timestepping in Flux Form'
            DO_3D_00_00(1,jpkm1)
               zue3b = e3u(ji,jj,jk,Nbb) * uu(ji,jj,jk,Nbb)
               zve3b = e3v(ji,jj,jk,Nbb) * vv(ji,jj,jk,Nbb)
               !                                                ! LF time stepping
               zue3a = zue3b + rDt * e3u(ji,jj,jk,Nrhs) * uu(ji,jj,jk,Nrhs) * umask(ji,jj,jk)
               zve3a = zve3b + rDt * e3v(ji,jj,jk,Nrhs) * vv(ji,jj,jk,Nrhs) * vmask(ji,jj,jk)
               !
               uu(ji,jj,jk,Naa) = zue3a / e3u(ji,jj,jk,Naa)
               vv(ji,jj,jk,Naa) = zve3a / e3v(ji,jj,jk,Naa)
            END_3D
         ELSE                             ! Leap Frog time stepping + Asselin filter
            DO_3D_11_11(1,jpkm1)
               zue3n = e3u(ji,jj,jk,Nnn) * uu(ji,jj,jk,Nnn)
               zve3n = e3v(ji,jj,jk,Nnn) * vv(ji,jj,jk,Nnn)
               zue3b = e3u(ji,jj,jk,Nbb) * uu(ji,jj,jk,Nbb)
               zve3b = e3v(ji,jj,jk,Nbb) * vv(ji,jj,jk,Nbb)
               !                                                ! LF time stepping
!!an from sibylle, correction e3u(Nrhs)->Nnn
              zue3a = zue3b + rDt * e3u(ji,jj,jk,Nnn) * uu(ji,jj,jk,Nrhs) * umask(ji,jj,jk)
              zve3a = zve3b + rDt * e3v(ji,jj,jk,Nnn) * vv(ji,jj,jk,Nrhs) * vmask(ji,jj,jk)
!!an  Nrhs = wrong !
              ! zue3a = zue3b + rDt * e3u(ji,jj,jk,Nrhs) * uu(ji,jj,jk,Nrhs) * umask(ji,jj,jk)
              ! zve3a = zve3b + rDt * e3v(ji,jj,jk,Nrhs) * vv(ji,jj,jk,Nrhs) * vmask(ji,jj,jk)
               !                                                ! Asselin time filter on e3u/v/t
               ze3u_tf = e3u(ji,jj,jk,Nnn) + rn_atfp * ( e3u(ji,jj,jk,Nbb) - 2._wp * e3u(ji,jj,jk,Nnn)  + e3u(ji,jj,jk,Naa) )
               ze3v_tf = e3v(ji,jj,jk,Nnn) + rn_atfp * ( e3v(ji,jj,jk,Nbb) - 2._wp * e3v(ji,jj,jk,Nnn)  + e3v(ji,jj,jk,Naa) )
               ze3t_tf = e3t(ji,jj,jk,Nnn) + rn_atfp * ( e3t(ji,jj,jk,Nbb) - 2._wp * e3t(ji,jj,jk,Nnn)  + e3t(ji,jj,jk,Naa) )
               !                                                ! Asselin time filter on u,v (Nnn)
               uu(ji,jj,jk,Nnn) = ( zue3n + rn_atfp * ( zue3b - 2._wp * zue3n  + zue3a ) ) / ze3u_tf
               vv(ji,jj,jk,Nnn) = ( zve3n + rn_atfp * ( zve3b - 2._wp * zve3n  + zve3a ) ) / ze3v_tf
               !
               e3u(ji,jj,jk,Nnn) = ze3u_tf
               e3v(ji,jj,jk,Nnn) = ze3v_tf
               e3t(ji,jj,jk,Nnn) = ze3t_tf
               !
               uu(ji,jj,jk,Naa) = zue3a / e3u(ji,jj,jk,Naa)
               vv(ji,jj,jk,Naa) = zve3a / e3v(ji,jj,jk,Naa)
            END_3D
         ENDIF
      ENDIF


      CALL lbc_lnk_multi( 'stp', uu(:,:,:,Nnn), 'U', -1., vv(:,:,:,Nnn), 'V', -1.,   &   !* local domain boundaries
         &                       uu(:,:,:,Naa), 'U', -1., vv(:,:,:,Naa), 'V', -1.    )

!!an shapiro filter 2D
!Falissard, Fabrice. (2013).
! DO_3D_11_11(1,jpkm1)
!     uu (ji,jj,jk,Nnn) = 0.0625_wp * (8*uu(ji,jj,jk,Nnn) + 2*uu(ji+1,jj,jk,Nnn) - uu(ji+2,jj,jk,Nnn) + uu(ji+1,jj+1,jk,Nnn) )
!     vv (ji,jj,jk,Nnn) = 0.0625_wp * (8*vv(ji,jj,jk,Nnn) + 2*vv(ji+1,jj,jk,Nnn) - vv(ji+2,jj,jk,Nnn) + vv(ji+1,jj+1,jk,Nnn) )
!     e3t(ji,jj,jk,Nnn) = 0.0625_wp * (8*e3t(ji,jj,jk,Nnn) + 2*e3t(ji+1,jj,jk,Nnn) - e3t(ji+2,jj,jk,Nnn) + e3t(ji+1,jj+1,jk,Nnn) )
!     e3u(ji,jj,jk,Nnn) = 0.0625_wp * (8*e3u(ji,jj,jk,Nnn) + 2*e3u(ji+1,jj,jk,Nnn) - e3u(ji+2,jj,jk,Nnn) + e3u(ji+1,jj+1,jk,Nnn) )
!     e3v(ji,jj,jk,Nnn) = 0.0625_wp * (8*e3v(ji,jj,jk,Nnn) + 2*e3v(ji+1,jj,jk,Nnn) - e3v(ji+2,jj,jk,Nnn) + e3v(ji+1,jj+1,jk,Nnn) )
!     !
!     uu (ji,jj,jk,Naa) = 0.0625_wp * (8*uu(ji,jj,jk,Naa) + 2*uu(ji+1,jj,jk,Naa) - uu(ji+2,jj,jk,Naa) + uu(ji+1,jj+1,jk,Naa) )
!     vv (ji,jj,jk,Naa) = 0.0625_wp * (8*vv(ji,jj,jk,Naa) + 2*vv(ji+1,jj,jk,Naa) - vv(ji+2,jj,jk,Naa) + vv(ji+1,jj+1,jk,Naa) )
!     e3t(ji,jj,jk,Naa) = 0.0625_wp * (8*e3t(ji,jj,jk,Naa) + 2*e3t(ji+1,jj,jk,Naa) - e3t(ji+2,jj,jk,Naa) + e3t(ji+1,jj+1,jk,Naa) )
!     e3u(ji,jj,jk,Naa) = 0.0625_wp * (8*e3u(ji,jj,jk,Naa) + 2*e3u(ji+1,jj,jk,Naa) - e3u(ji+2,jj,jk,Naa) + e3u(ji+1,jj+1,jk,Naa) )
!     e3v(ji,jj,jk,Naa) = 0.0625_wp * (8*e3v(ji,jj,jk,Naa) + 2*e3v(ji+1,jj,jk,Naa) - e3v(ji+2,jj,jk,Naa) + e3v(ji+1,jj+1,jk,Naa) )
! END_3D
! CALL lbc_lnk_multi( 'stp', uu(:,:,:,Nnn), 'U', -1., vv(:,:,:,Nnn), 'V', -1.,   &   !* local domain boundaries
!    &                       uu(:,:,:,Naa), 'U', -1., vv(:,:,:,Naa), 'V', -1.    )
! CALL lbc_lnk_multi( 'stp', e3t(:,:,:,Nnn), 'T', 1., e3u(:,:,:,Nnn), 'U', -1.,e3v(:,:,:,Nnn), 'V', -1.,   &   !* local domain boundaries
!   &                        e3t(:,:,:,Naa), 'T', 1., e3u(:,:,:,Naa), 'U', -1.,e3v(:,:,:,Naa), 'V', -1.    )
!DO_3D_11_11(1,jpkm1)
!     uu (ji,jj,jk,Nnn) = 0.25 * (2*uu(ji,jj,jk,Nnn) + uu(ji+1,jj,jk,Nnn) - uu(ji-1,jj,jk,Nnn) )
!     vv (ji,jj,jk,Nnn) = 0.25 * (2*uu(ji,jj,jk,Nnn) + uu(ji+1,jj,jk,Nnn) - uu(ji-1,jj,jk,Nnn) )
!     e3t(ji,jj,jk,Nnn) = 0.25 * (2*uu(ji,jj,jk,Nnn) + uu(ji+1,jj,jk,Nnn) - uu(ji-1,jj,jk,Nnn) )
!     e3u(ji,jj,jk,Nnn) = 0.25 * (2*uu(ji,jj,jk,Nnn) + uu(ji+1,jj,jk,Nnn) - uu(ji-1,jj,jk,Nnn) )
!     e3v(ji,jj,jk,Nnn) = 0.25 * (2*uu(ji,jj,jk,Nnn) + uu(ji+1,jj,jk,Nnn) - uu(ji-1,jj,jk,Nnn) )
! END_3D
! CALL lbc_lnk_multi( 'stp', uu(:,:,:,Nnn), 'U', -1., vv(:,:,:,Nnn), 'V', -1.,   &   !* local domain boundaries
!    &                       uu(:,:,:,Naa), 'U', -1., vv(:,:,:,Naa), 'V', -1.    )

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Set boundary conditions, time filter and swap time levels
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!!an TO BE ADDED : dyn_nxt
!!                         CALL dyn_atf       ( kstp, Nbb, Nnn, Naa, uu, vv, e3t, e3u, e3v  )  ! time filtering of "now" velocities and scale factors
!!an TO BE ADDED : a simplifier
!!                         CALL ssh_atf       ( kstp, Nbb, Nnn, Naa, ssh )                     ! time filtering of "now" sea surface height
      IF ( .NOT.( l_1st_euler ) ) THEN   ! Only do time filtering for leapfrog timesteps
         !                                                  ! filtering "now" field
         ssh(:,:,Nnn) = ssh(:,:,Nnn) + rn_atfp * ( ssh(:,:,Nbb) - 2 * ssh(:,:,Nnn) + ssh(:,:,Naa) )
      ENDIF
!!an


      ! Swap time levels
      Nrhs = Nbb
      Nbb = Nnn
      Nnn = Naa
      Naa = Nrhs
      !
      !!an recompute scale factor for ?
                         CALL dom_vvl_sf_update( kstp, Nbb, Nnn, Naa )  ! recompute vertical scale factors
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! diagnostics and outputs
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_floats  )   CALL flo_stp   ( kstp, Nbb, Nnn )      ! drifting Floats
      IF( ln_diacfl  )   CALL dia_cfl   ( kstp,      Nnn )      ! Courant number diagnostics

                         CALL dia_wri   ( kstp,      Nnn )      ! ocean model: outputs
      !
      IF( lrst_oce   )   CALL rst_write    ( kstp, Nbb, Nnn )   ! write output ocean restart file


#if defined key_agrif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! AGRIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         Kbb_a = Nbb; Kmm_a = Nnn; Krhs_a = Nrhs      ! agrif_oce module copies of time level indices
                         CALL Agrif_Integrate_ChildGrids( stp )       ! allows to finish all the Child Grids before updating

                         IF( Agrif_NbStepint() == 0 ) THEN
                            CALL Agrif_update_all( )                  ! Update all components
                         ENDIF
#endif
      IF( ln_diaobs  )   CALL dia_obs      ( kstp, Nnn )      ! obs-minus-model (assimilation) diagnostics (call after dynamics update)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL stp_ctl      ( kstp, Nbb, Nnn, indic )


      IF( kstp == nit000 ) THEN                          ! 1st time step only
                                        CALL iom_close( numror )   ! close input  ocean restart file
         IF(lwm)                        CALL FLUSH    ( numond )   ! flush output namelist oce
         IF(lwm .AND. numoni /= -1 )    CALL FLUSH    ( numoni )   ! flush output namelist ice (if exist)
      ENDIF

      !
#if defined key_iomput
      IF( kstp == nitend .OR. indic < 0 ) THEN
                      CALL iom_context_finalize(      cxios_context          ) ! needed for XIOS+AGRIF
                      IF(lrxios) CALL iom_context_finalize(      crxios_context          )
      ENDIF
#endif
      !
      IF( l_1st_euler ) THEN         ! recover Leap-frog timestep
         rDt = 2._wp * rn_Dt
         r1_Dt = 1._wp / rDt
         l_1st_euler = .FALSE.
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('stp')
      !
   END SUBROUTINE stp
#endif
   !
   !!======================================================================
END MODULE step
