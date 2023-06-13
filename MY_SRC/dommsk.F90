MODULE dommsk
   !!======================================================================
   !!                       ***  MODULE dommsk   ***
   !! Ocean initialization : domain land/sea mask
   !!======================================================================
   !! History :  OPA  ! 1987-07  (G. Madec)  Original code
   !!            6.0  ! 1993-03  (M. Guyon)  symetrical conditions (M. Guyon)
   !!            7.0  ! 1996-01  (G. Madec)  suppression of common work arrays
   !!             -   ! 1996-05  (G. Madec)  mask computed from tmask
   !!            8.0  ! 1997-02  (G. Madec)  mesh information put in domhgr.F
   !!            8.1  ! 1997-07  (G. Madec)  modification of kbat and fmask
   !!             -   ! 1998-05  (G. Roullet)  free surface
   !!            8.2  ! 2000-03  (G. Madec)  no slip accurate
   !!             -   ! 2001-09  (J.-M. Molines)  Open boundaries
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and module
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!            3.6  ! 2015-05  (P. Mathiot) ISF: add wmask,wumask and wvmask
   !!            4.0  ! 2016-06  (G. Madec, S. Flavoni)  domain configuration / user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_msk       : compute land/ocean mask
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE usrdef_fmask   ! user defined fmask
   USE bdy_oce        ! open boundary
   !
   USE in_out_manager ! I/O manager
   USE iom            ! IOM library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! Massively Parallel Processing library
   !
   USE usrdef_nam , ONLY : rn_abp, r1_abp, rn_fsp, rn_cnp, rn_dx, nn_AM98, nn_fsp,nn_smo, ln_hdiv_AD

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_msk    ! routine called by inidom.F90

   !                            !!* Namelist namlbc : lateral boundary condition *
!!an rn_shlat changed to PUBLIC
   REAL(wp), PUBLIC :: rn_shlat   ! type of lateral boundary condition on velocity
   LOGICAL , PUBLIC :: ln_vorlat  !  consistency of vorticity boundary condition
   !                                            with analytical eqs.

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dommsk.F90 13416 2020-08-20 10:10:55Z gm $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_msk( k_top, k_bot )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-
      !!      zontal velocity points (u & v), vorticity points (f) points.
      !!
      !! ** Method  :   The ocean/land mask  at t-point is deduced from ko_top
      !!      and ko_bot, the indices of the fist and last ocean t-levels which
      !!      are either defined in usrdef_zgr or read in zgr_read.
      !!                The velocity masks (umask, vmask, wmask, wumask, wvmask)
      !!      are deduced from a product of the two neighboring tmask.
      !!                The vorticity mask (fmask) is deduced from tmask taking
      !!      into account the choice of lateral boundary condition (rn_shlat) :
      !!         rn_shlat = 0, free slip  (no shear along the coast)
      !!         rn_shlat = 2, no slip  (specified zero velocity at the coast)
      !!         0 < rn_shlat < 2, partial slip   | non-linear velocity profile
      !!         2 < rn_shlat, strong slip        | in the lateral boundary layer
      !!
      !!      tmask_i : interior ocean mask at t-point, i.e. excluding duplicated
      !!                rows/lines due to cyclic or North Fold boundaries as well
      !!                as MPP halos.
      !!      tmask_h : halo mask at t-point, i.e. excluding duplicated rows/lines
      !!                due to cyclic or North Fold boundaries as well as MPP halos.
      !!
      !! ** Action :   tmask, umask, vmask, wmask, wumask, wvmask : land/ocean mask
      !!                         at t-, u-, v- w, wu-, and wv-points (=0. or 1.)
      !!               fmask   : land/ocean mask at f-point (=0., or =1., or
      !!                         =rn_shlat along lateral boundaries)
      !!               tmask_i : interior ocean mask
      !!               tmask_h : halo mask
      !!               ssmask , ssumask, ssvmask, ssfmask : 2D ocean mask
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), INTENT(in) ::   k_top, k_bot   ! first and last ocean level
      !
      INTEGER  ::   ji, jj, jk, jl     ! dummy loop indices
      INTEGER  ::   iif, iil       ! local integers
      INTEGER  ::   ijf, ijl       !   -       -
      INTEGER  ::   iktop, ikbot   !   -       -
      INTEGER  ::   ios, inum
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zwf   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   z3d   ! 3D workspace
      REAL(wp) :: zES,zEN,zWS,zWN,zl,z0,ze1, ze2, ze3,z1d, z1x, z1y, z1x1, z1x2, z1y1, z1y2 ! local real
      !!
      NAMELIST/namlbc/ rn_shlat, ln_vorlat
      NAMELIST/nambdy/ ln_bdy ,nb_bdy, ln_coords_file, cn_coords_file,         &
         &             ln_mask_file, cn_mask_file, cn_dyn2d, nn_dyn2d_dta,     &
         &             cn_dyn3d, nn_dyn3d_dta, cn_tra, nn_tra_dta,             &
         &             ln_tra_dmp, ln_dyn3d_dmp, rn_time_dmp, rn_time_dmp_out, &
         &             cn_ice, nn_ice_dta,                                     &
         &             ln_vol, nn_volctl, nn_rimwidth
      !!---------------------------------------------------------------------
      !
      READ  ( numnam_ref, namlbc, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namlbc in reference namelist' )
      READ  ( numnam_cfg, namlbc, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namlbc in configuration namelist' )
      IF(lwm) WRITE ( numond, namlbc )

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dommsk : ocean mask '
         WRITE(numout,*) '~~~~~~'
         WRITE(numout,*) '   Namelist namlbc'
         WRITE(numout,*) '      lateral momentum boundary cond.    rn_shlat  = ',rn_shlat
         WRITE(numout,*) '      consistency with analytical form   ln_vorlat = ',ln_vorlat
      ENDIF
      !
      IF(lwp) WRITE(numout,*)
      IF     (      rn_shlat == 0.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ==>>>   ocean lateral  free-slip'
      ELSEIF (      rn_shlat == 2.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ==>>>   ocean lateral  no-slip'
      ELSEIF ( 0. < rn_shlat .AND. rn_shlat < 2. ) THEN   ;   IF(lwp) WRITE(numout,*) '   ==>>>   ocean lateral  partial-slip'
      ELSEIF ( 2. < rn_shlat                     ) THEN   ;   IF(lwp) WRITE(numout,*) '   ==>>>   ocean lateral  strong-slip'
      ELSE
         CALL ctl_stop( 'dom_msk: wrong value for rn_shlat (i.e. a negalive value). We stop.' )
      ENDIF

      !  Ocean/land mask at t-point  (computed from ko_top and ko_bot)
      ! ----------------------------
      !
      tmask(:,:,:) = 0._wp
      DO_2D_11_11
         iktop = k_top(ji,jj)
         ikbot = k_bot(ji,jj)
         IF( iktop /= 0 ) THEN       ! water in the column
            tmask(ji,jj,iktop:ikbot  ) = 1._wp
         ENDIF
      END_2D
      !
      ! the following call is mandatory
      ! it masks boundaries (bathy=0) where needed depending on the configuration (closed, periodic...)
      CALL lbc_lnk( 'dommsk', tmask  , 'T', 1._wp )      ! Lateral boundary conditions
      !
     ! Mask corrections for bdy (read in mppini2)
      READ  ( numnam_ref, nambdy, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambdy in reference namelist' )
      READ  ( numnam_cfg, nambdy, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'nambdy in configuration namelist' )
      ! ------------------------
      IF ( ln_bdy .AND. ln_mask_file ) THEN
         CALL iom_open( cn_mask_file, inum )
         CALL iom_get ( inum, jpdom_data, 'bdy_msk', bdytmask(:,:) )
         CALL iom_close( inum )
         DO_3D_11_11( 1, jpkm1 )
            tmask(ji,jj,jk) = tmask(ji,jj,jk) * bdytmask(ji,jj)
         END_3D
      ENDIF
      !
      !  Ocean/land mask at u-, v-, and f-points   (computed from tmask)
      ! ----------------------------------------
      ! NB: at this point, fmask is designed for free slip lateral boundary condition
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector loop
               umask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)
               vmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji  ,jj+1,jk)
            END DO
            DO ji = 1, jpim1      ! NO vector opt.
               fmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
                  &            * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
            END DO
         END DO
      END DO
      CALL lbc_lnk_multi( 'dommsk', umask, 'U', 1., vmask, 'V', 1., fmask, 'F', 1. )      ! Lateral boundary conditions

      ! Ocean/land mask at wu-, wv- and w points    (computed from tmask)
      !-----------------------------------------
      wmask (:,:,1) = tmask(:,:,1)     ! surface
      wumask(:,:,1) = umask(:,:,1)
      wvmask(:,:,1) = vmask(:,:,1)
      DO jk = 2, jpk                   ! interior values
         wmask (:,:,jk) = tmask(:,:,jk) * tmask(:,:,jk-1)
         wumask(:,:,jk) = umask(:,:,jk) * umask(:,:,jk-1)
         wvmask(:,:,jk) = vmask(:,:,jk) * vmask(:,:,jk-1)
      END DO


      ! Ocean/land column mask at t-, u-, and v-points   (i.e. at least 1 wet cell in the vertical)
      ! ----------------------------------------------
      ssmask (:,:) = MAXVAL( tmask(:,:,:), DIM=3 )
      ssumask(:,:) = MAXVAL( umask(:,:,:), DIM=3 )
      ssvmask(:,:) = MAXVAL( vmask(:,:,:), DIM=3 )
      ! ssfmask(:,:) = MAXVAL( fmask(:,:,:), DIM=3 )
      DO_2D_10_10
         ssfmask(ji,jj) = MAX(  tmask(ji,jj+1,1), tmask(ji+1,jj+1,1),  &
            &                   tmask(ji,jj  ,1), tmask(ji+1,jj  ,1)   )
      END_2D
      CALL lbc_lnk( 'dommsk', ssfmask, 'F', 1._wp )

# if defined key_bvp
      !  Penalize the first inner layer of fluid
      ! -----------------------------------------------------
!!an to be added : lissage sur rn_cnp dx
      rpo    (:,:,:) = 1._wp
      rpou   (:,:,:) = 1._wp
      rpov   (:,:,:) = 1._wp
      rpof   (:,:,:) = 1._wp
      !
      bmpt(:,:,:) = 1._wp
      !      !
      WHERE(tmask(:,:,:) == 0._wp)
        rpo   (:,:,:) = rn_abp             ! used on T points
        bmpt  (:,:,:) = 0._wp              ! used on T points
      END WHERE
      WHERE(umask(:,:,:) == 0._wp)
        rpou  (:,:,:) = rn_abp             ! used on U points
      END WHERE
      WHERE(vmask(:,:,:) == 0._wp)
        rpov  (:,:,:) = rn_abp             ! used on V points
      END WHERE
      ! CAUTION
      WHERE(fmask(:,:,:) == 0._wp)
        rpof (:,:,:) = rn_abp              ! used on F points
      END WHERE
      !
      CALL lbc_lnk_multi( 'dommsk', rpo ,  'T', 1._wp,                      &
              &                     bmpt,  'T', 1._wp,                      &
              &                     rpou,  'U', 1._wp,                      &
              &                     rpov,  'V', 1._wp,                      &
              &                     rpof,  'F', 1._wp, kfillmode=jpfillcopy )
      !  tanh penalisation (T point) coast = middle of rpo
      ! ----------------------------------------
      z1d =  rn_dx / REAL(nn_AM98, wp)
      DO jk = 1, jpkm1
        DO_2D_00_00
          z1x =    tanh( (glamt(ji,jj) - 0._wp       )/(rn_cnp * z1d * 0.25_wp) ) +       &
              &  - tanh( (glamt(ji,jj) - 2000000._wp )/(rn_cnp * z1d * 0.25_wp) )
          z1x = z1x * 0.5_wp
          !
          z1y =    tanh( (gphit(ji,jj) - 0._wp       )/(rn_cnp * z1d * 0.25_wp) ) +       &
              &  - tanh( (gphit(ji,jj) - 2000000._wp )/(rn_cnp * z1d * 0.25_wp) )
          z1y = z1y * 0.5_wp
          !
          rpo(ji,jj,jk) = MIN(z1x,z1y)
          !
        END_2D
      END DO
      CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      WHERE(rpo(:,:,:) < rn_abp)
         rpo   (:,:,:) = rn_abp             ! used on T points
      END WHERE
      CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      !!an
      ! pour le fun
      ! WRITE(numout,*) 'dommsk : Pour le fun...'
      ! rpo(:,:,:) = 1._wp + rn_abp - rpo(:,:,:)
      !!an
      !
      !  linear penalisation (speed points) coast = middle of porosity
      ! -------------------------------------------------------------
      ! z1d =  rn_dx / REAL(nn_AM98, wp)
      ! DO jk = 1, jpkm1
      !   DO_2D_00_00
      !     z1x1 =  (1._wp - rn_abp)*(glamt(ji,jj) -       0._wp)/(rn_cnp * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     z1x1 = MAX(MIN(z1x1,1._wp), rn_abp)
      !     z1x2 = -(1._wp - rn_abp)*(glamt(ji,jj) - 2000000._wp)/(rn_cnp * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     z1x2 = MAX(MIN(z1x2,1._wp), rn_abp)
      !     z1x  = MIN(z1x1,z1x2)
      !     !
      !     z1y1 =  (1._wp - rn_abp)*(gphit(ji,jj) -       0._wp)/(rn_cnp * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     z1y1 = MAX(MIN(z1y1,1._wp), rn_abp)
      !     z1y2 = -(1._wp - rn_abp)*(gphit(ji,jj) - 2000000._wp)/(rn_cnp * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     z1y2 = MAX(MIN(z1y2,1._wp), rn_abp)
      !     z1y  = MIN(z1y1,z1y2)
      !     !
      !     rpo(ji,jj,jk) = MIN(z1x,z1y)
      !     ! rpo(ji,jj,jk) = 1._wp
      !     !
      !   END_2D
      ! END DO
      ! CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      ! WHERE(tmask(:,:,:) == 0._wp)
      !   rpo   (:,:,:) = rn_abp             ! used on T points
      ! END WHERE
      ! ! !!an
      ! ! rpo(:,:,:) = rn_abp
      ! ! !!an
      ! CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      !
      !  linear penalisation only on NS or WS coast
      ! ---------------------------------------------------------
      ! z1d =  rn_dx / REAL(nn_AM98, wp)
      ! DO jk = 1, jpkm1
      !   DO_2D_00_00
      !     !
      !     !                                !== South West Coast light  ==!
      !     ! W
      !     ! WRITE(numout,*) 'dommsk : West coast lightly penalised'
      !     ! z1x1 =  (1._wp - rn_abp)*(glamt(ji,jj) -       0._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     ! ze2 = MAX(MIN(z1x1,1._wp), rn_abp)
      !     ! !
      !     ! z1y1 = -(1._wp - ze2)*(gphit(ji,jj) -  15._wp*z1d)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + ze2)
      !     ! ze1 = MAX(MIN(z1y1,1._wp), ze2)
      !     ! !
      !     ! ! bmpt(ji,jj,jk) = ze1
      !     ! ! bmpt(ji,jj,jk) = MIN(ze1,1._wp)
      !     ! ! bmpt(ji,jj,jk) = MAX(MIN(ze1,1._wp), rn_abp)
      !     ! rpo(ji,jj,jk) = ze1
      !     !
      !     !                                !== South West Coast + SW corner  ==!
      !     ! W
      !     ! WRITE(numout,*) 'dommsk : West coast penalised'
      !     ! z1x1 =  (1._wp - rn_abp)*(glamt(ji,jj) -       0._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     ! ze2 = MAX(MIN(z1x1,1._wp), rn_abp)
      !     ! !
      !     ! zl = 1._wp * REAL(15,wp)*z1d ! radius   !   so rpo=1/2 take off the real coastline (3 - sqrt(2))
      !     ! z0 = zl/3._wp                    ! center   !   so the mean take off the real coastline (1 - 1/3)
      !     ! WRITE(numout,*) 'dommsk : West-South corner circle-penalised'
      !     ! !                          so rpo=rn_abp take off the real coastline
      !     ! z1x =  ( glamt(ji,jj) - (               z0 - REAL(rn_cnp,wp)*z1d/2 ) )/zl
      !     ! z1y =  ( gphit(ji,jj) - (               z0 - REAL(rn_cnp,wp)*z1d/2 ) )/zl
      !     ! zWS   = - (z1x*z1x + z1y*z1y - 1)  !circle
      !     ! !
      !     ! zWS =  MAX(MIN(zWS,1._wp), rn_abp)
      !     ! ! IF ( z1y > 0._wp .OR. z1x > 0._wp ) THEN   ! quart de cercle coin
      !     ! !   zWS = 1._wp
      !     ! ! ENDIF
      !     ! !
      !     ! ze1 = zWS
      !     ! IF (     z1x > 0._wp ) THEN   ! quart de cercle coin
      !     !   ze1 = 1._wp   ! in the ocean
      !     ! ELSEIF ( z1y > 0._wp )THEN
      !     !   ze1 = ze2    ! in the West coast penalised
      !     ! ENDIF          ! in the quarter circle
      !     ! ! bmpt(ji,jj,jk) = ze1
      !     ! ! bmpt(ji,jj,jk) = MIN(ze1,1._wp)
      !     ! ! bmpt(ji,jj,jk) = MAX(MIN(ze1,1._wp), rn_abp)
      !     ! rpo(ji,jj,jk) = ze1
      !     ! !
      !     !                                           !   tought for d=15 and gc = 8
      !     ! zl = 3._wp * REAL(rn_cnp,wp)*z1d ! radius   !   so rpo=1/2 take off the real coastline (3 - sqrt(2))
      !     ! z0 = zl/3._wp                    ! center   !   so the mean take off the real coastline (1 - 1/3)
      !     ! !
      !     ! !                               !==  West South Circle Corner  ==!
      !     ! WRITE(numout,*) 'dommsk : West-South corner circle-penalised'
      !     ! !
      !     ! z1x =  ( glamt(ji,jj) - (               z0 - REAL(rn_cnp,wp)*z1d/2 ) )/zl
      !     ! z1y =  ( gphit(ji,jj) - (               z0 - REAL(rn_cnp,wp)*z1d/2 ) )/zl
      !     ! !
      !     ! zWS   = - (z1x*z1x + z1y*z1y - 1)  !circle
      !     ! !
      !     ! IF ( z1y > 0._wp .OR. z1x > 0._wp ) THEN   ! quart de cercle coin
      !     !   zWS = 1._wp
      !     ! ENDIF
      !     ! zWS =  MAX(MIN(zWS,1._wp), rn_abp)
      !     ! ! rpo(ji,jj,jk) = zWS
      !     ! !
      !     ! !                               !==  East South Circle Corner  ==!
      !     ! WRITE(numout,*) 'dommsk : East-South corner circle-penalised'
      !     ! !
      !     ! z1x =  ( glamt(ji,jj) - ( 2000000._wp - z0 + REAL(rn_cnp,wp)*z1d/2 ) )/zl
      !     ! z1y =  ( gphit(ji,jj) - (               z0 - REAL(rn_cnp,wp)*z1d/2 ) )/zl
      !     ! !
      !     ! zES   = - (z1x*z1x + z1y*z1y - 1)  !circle
      !     ! !
      !     ! IF ( z1y > 0._wp .OR. z1x < 0._wp ) THEN   ! quart de cercle coin
      !     ! zES = 1._wp
      !     ! ENDIF
      !     ! zES =  MAX(MIN(zES,1._wp), rn_abp)
      !     ! rpo(ji,jj,jk) = zES
      !     ! !
      !     ! !                               !==  West North Circle Corner  ==!
      !     ! WRITE(numout,*) 'dommsk : West-North corner circle-penalised'
      !     ! !
      !     ! z1x =  ( glamt(ji,jj) - (               z0 - REAL(rn_cnp,wp)*z1d/2 ) )/zl
      !     ! z1y =  ( gphit(ji,jj) - ( 2000000._wp - z0 + REAL(rn_cnp,wp)*z1d/2 ) )/zl
      !     ! !
      !     ! zWN   = - (z1x*z1x + z1y*z1y - 1)  !circle
      !     ! !
      !     ! IF ( z1y < 0._wp .OR. z1x > 0._wp ) THEN   ! quart de cercle coin
      !     !   zWN = 1._wp
      !     ! ENDIF
      !     ! !
      !     ! zWN = MAX(MIN(zWN,1._wp), rn_abp)
      !     ! ! rpo(ji,jj,jk) = zWN
      !     ! !
      !     ! !                             !==  East North Circle Corner  ==!
      !     ! WRITE(numout,*) 'dommsk : East-North corner circle-penalised'
      !     ! !
      !     ! z1x =  ( glamt(ji,jj) - ( 2000000._wp - z0 + REAL(rn_cnp,wp)*z1d/2 ) )/zl
      !     ! z1y =  ( gphit(ji,jj) - ( 2000000._wp - z0 + REAL(rn_cnp,wp)*z1d/2 ) )/zl
      !     ! !
      !     ! zEN   = - (z1x*z1x + z1y*z1y - 1)  !circle
      !     ! !
      !     ! IF ( z1y < 0._wp .OR. z1x < 0._wp ) THEN   ! quart de cercle coin
      !     ! zEN = 1._wp
      !     ! ENDIF
      !     ! !
      !     ! zEN = MAX(MIN(zEN,1._wp), rn_abp)
      !     ! ! rpo(ji,jj,jk) = zEN
      !     ! !
      !     ! !          !==  West (North+South) Corners  ==!
      !     ! ! WRITE(numout,*) 'dommsk : North corners circle-penalised'
      !     ! ! rpo(ji,jj,jk) = MIN(zWS,zWN)
      !     ! !                   !==  All Corners  ==!
      !     ! WRITE(numout,*) 'dommsk : All four corners circle-penalised'
      !     ! rpo(ji,jj,jk) = MIN(zWS,zWN,zES,zEN)   ! ils sont chacun indépendanant
      !     !
      !     !                                !==  South West coastlines  ==!
      !     !! SW coast
      !     ! WRITE(numout,*) 'dommsk : West-South coast penalised'
      !     ! z1x1   =  (1._wp - rn_abp)*(glamt(ji,jj) -       0._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     ! z1x1 = MAX(MIN(z1x1,1._wp), rn_abp)
      !     ! !
      !     ! z1y1 =  (1._wp - rn_abp)*(gphit(ji,jj) -       0._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     ! z1y1 = MAX(MIN(z1y1,1._wp), rn_abp)
      !     ! z1y  = MIN(z1y1,z1y2)
      !     ! !
      !     ! rpo(ji,jj,jk) = MIN(z1x1, z1y)
      !     !
      !     ! W
      !     WRITE(numout,*) 'dommsk : West coast penalised'
      !     z1x1 =  (1._wp - rn_abp)*(glamt(ji,jj) -       0._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     z1x1 = MAX(MIN(z1x1,1._wp), rn_abp)
      !     !
      !     ! ! EW
      !     ! WRITE(numout,*) 'dommsk : West-East coast penalised'
      !     ! z1x1 =  (1._wp - rn_abp)*(glamt(ji,jj) -       0._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     ! z1x1 = MAX(MIN(z1x1,1._wp), rn_abp)
      !     ! z1x2 = -(1._wp - rn_abp)*(glamt(ji,jj) - 2000000._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     ! z1x2 = MAX(MIN(z1x2,1._wp), rn_abp)
      !     ! z1x  = MIN(z1x1,z1x2)
      !     !
      !     !NS
      !     ! WRITE(numout,*) 'dommsk : North-South coast penalised'
      !     ! z1y1 =  (1._wp - rn_abp)*(gphit(ji,jj) -       0._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     ! z1y1 = MAX(MIN(z1y1,1._wp), rn_abp)
      !     ! z1y2 = -(1._wp - rn_abp)*(gphit(ji,jj) - 2000000._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
      !     ! z1y2 = MAX(MIN(z1y2,1._wp), rn_abp)
      !     ! z1y  = MIN(z1y1,z1y2)
      !     !
      !     ! rpo(ji,jj,jk) = MIN(z1x,z1y)
      !     rpo(ji,jj,jk) = z1x1   ! W
      !     !
      !   END_2D
      ! END DO
      ! CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      ! WHERE(tmask(:,:,:) == 0._wp)
      !   rpo   (:,:,:) = rn_abp             ! used on T points
      ! END WHERE
      ! ! !!an
      ! ! rpo(:,:,:) = rn_abp
      ! ! !!an
      ! CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )

      !
      !  tanh penalisation (T point) coast = min of rpo
      ! ------------------------------------------------
      ! z1d =  rn_dx / REAL(nn_AM98, wp)
      ! DO jk = 1, jpkm1
      !   DO_2D_00_00
      !     z1x =    tanh( (glamt(ji,jj) - 0._wp       - rn_cnp*z1d )/(rn_cnp * z1d) )         &
      !         &  - tanh( (glamt(ji,jj) - 2000000._wp + rn_cnp*z1d )/(rn_cnp * z1d) )
      !     z1x = z1x * 0.5_wp
      !     z1x = (1._wp - rn_abp) * z1x + rn_abp
      !     !
      !     z1y =    tanh( (gphit(ji,jj) - 0._wp       - rn_cnp*z1d )/(rn_cnp * z1d) )         &
      !         &  - tanh( (gphit(ji,jj) - 2000000._wp + rn_cnp*z1d )/(rn_cnp * z1d) )
      !     z1y = z1y * 0.5_wp
      !     z1y = (1._wp - rn_abp) * z1y + rn_abp
      !     !
      !     rpo(ji,jj,jk) = MIN(z1x,z1y)
      !     !
      !   END_2D
      ! END DO
      ! CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      ! WHERE(tmask(:,:,:) == 0._wp)
      !   rpo   (:,:,:) = rn_abp             ! used on T points
      ! END WHERE
      ! CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      !  tanh penalisation (T point) coast = middle of rpo
      ! ----------------------------------------
      ! z1d =  rn_dx / REAL(nn_AM98, wp)
      ! DO jk = 1, jpkm1
      !   DO_2D_00_00
      !     z1x =    tanh( (glamt(ji,jj) - 0._wp       - z1d*0.5_wp )/(rn_cnp * z1d) ) +       &
      !         &  - tanh( (glamt(ji,jj) - 2000000._wp + z1d*0.5_wp )/(rn_cnp * z1d) )
      !     z1x = z1x * 0.5_wp
      !     !
      !     z1y =    tanh( (gphit(ji,jj) - 0._wp       - z1d*0.5_wp )/(rn_cnp * z1d) ) +       &
      !         &  - tanh( (gphit(ji,jj) - 2000000._wp + z1d*0.5_wp )/(rn_cnp * z1d) )
      !     z1y = z1y * 0.5_wp
      !     !
      !     rpo(ji,jj,jk) = MIN(z1x,z1y)
      !     !
      !   END_2D
      ! END DO
      ! CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      ! WHERE(rpo(:,:,:) < rn_abp)
      !    rpo   (:,:,:) = rn_abp             ! used on T points
      ! END WHERE
      ! CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      !
      ! rn_abp = 1.e-2    EXP06 (with smoothing)
      ! goes 1 cell off shore at 45 degree
      ! first ocean cell at 0 degree
      !                                                                    !== T-points ==!
      ! DO jk = 1, jpkm1
      !   DO_2D_00_00
      !     IF ( tmask(ji,jj,1) == 1._wp) THEN
      !        IF  ( (      fmask(ji-1,jj  ,1) + fmask(ji,jj  ,1)           +             &
      !         &           fmask(ji-1,jj-1,1) + fmask(ji,jj-1,1)           )  < 4._wp    )   THEN
      !           !
      !           rpo   (ji,jj,1:jpkm1) = rn_abp             ! T points
      !           !
      !        ENDIF
      !     ENDIF
      !   END_2D
      ! END DO
      ! CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      ! !                                                                    !== F-points ==!
      ! DO jk = 1, jpkm1
      !   DO_2D_00_00
      !     IF ( fmask(ji,jj,1) == 1._wp) THEN
      !       IF             ( (      fmask(ji-1,jj+1,1) + fmask(ji+1,jj+1,1)           +             &
      !           &                   fmask(ji-1,jj-1,1) + fmask(ji+1,jj-1,1)           )  < 4._wp    )   THEN
      !           ! side
      !           rpof(ji,jj,1) = rn_abp
      !           !
      !       ELSE   IF      ( (      fmask(ji-1,jj+1,1) + fmask(ji+1,jj+1,1)           +             &
      !             &                 fmask(ji-1,jj-1,1) + fmask(ji+1,jj-1,1)           ) == 1._wp    )   THEN
      !           ! corners
      !           rpof(ji,jj,1) = rn_abp/2._wp
      !           !
      !       ENDIF
      !     ENDIF
      !   END_2D
      ! END DO
      ! CALL lbc_lnk( 'dommsk', rpof,  'F', 1._wp, kfillmode=jpfillcopy )
      !
      ! rn_abp = 0.5      EXP03 (without smoothing) and EXP04 (with smoothing) or EXP05 (without smoothing)
      !                   EXP07
      !                   EXP11 (non rotated basin) : rpo and rpov/rpou
      !                !== T-points ==!
      ! Nearest neighbourg
      ! DO jk = 1, jpkm1
      !   DO_2D_00_00
      !      IF( ( (                    vmask(ji,jj  ,jk)                     +      &
      !       &      umask(ji-1,jj,jk) +                  umask(ji,jj,jk)     +      &
      !       &                         vmask(ji,jj-1,jk)        )  < 3._wp  ) .AND. &
      !       &                         tmask(ji,jj,jk) == 1._wp                 )   THEN
      !         !
      !         rpo   (ji,jj,jk) = rn_abp             ! T points
      !         !
      !      ENDIF
      !      IF( ( (                    vmask(ji,jj  ,jk)                     +      &
      !       &      umask(ji-1,jj,jk) +                  umask(ji,jj,jk)     +      &
      !       &                         vmask(ji,jj-1,jk)        )  == 3._wp  ) .AND. &
      !       &                         tmask(ji,jj,jk) == 1._wp                 )   THEN
      !         !
      !         rpo   (ji,jj,jk) = 0.5_wp * rn_abp             ! T points
      !         !
      !      ENDIF
      !   END_2D
      ! END DO
      ! CALL lbc_lnk( 'dommsk', rpo ,  'T', 1._wp, kfillmode=jpfillcopy )
      !
      ! utile pour le 0°
      ! DO jk = 1,jpkm1
      !   DO_2D_10_10
      !      rpou(ji,jj,jk) = MAX( rpo(ji,jj,jk), rpo(ji+1,jj  ,jk) )
      !      rpov(ji,jj,jk) = MAX( rpo(ji,jj,jk), rpo(ji  ,jj+1,jk) )
      !   END_2D
      ! END DO
      ! CALL lbc_lnk_multi( 'dommsk', rpov  , 'V', 1._wp, rpou  , 'U', 1._wp, kfillmode=jpfillcopy )
      !
      !
      !              !== F-points ==!
      ! half e1e2f (45 rotated)
      ! DO jk = 1, jpkm1
      !   DO_2D_00_00
      !       IF      (       (    umask(ji,jj,jk) + umask(ji  ,jj+1,jk)    &
      !           &           +    vmask(ji,jj,jk) + vmask(ji+1,jj  ,jk)    )  == 2._wp   ) THEN
      !       ! IF      ( ssfmask(ji,jj) + fmask(ji,jj,jk) == 1._wp   ) THEN
      !         ! side
      !         rpof   (ji,jj,jk) = rn_abp
      !         !
      !     ELSE IF (       (    umask(ji,jj,jk) + umask(ji  ,jj+1,jk)    &
      !         &           +    vmask(ji,jj,jk) + vmask(ji+1,jj  ,jk)    )  == 1._wp   ) THEN
      !         ! corners
      !         rpof   (ji,jj,jk) = 0.5_wp * rn_abp
      !         !
      !     ENDIF
      !   END_2D
      ! END DO
      ! CALL lbc_lnk( 'dommsk', rpof,  'F', 1._wp, kfillmode=jpfillcopy )
      !
      !
      ! !
      !
      ! Spontaneous penalisation (mean)
      ! DO_3D_00_00(1,jpkm1)
      !    IF( ( (       vmask(ji-1,jj,1) + vmask(ji+1,jj,1) ) < 2._wp   )  .AND. &
      !     &                         tmask(ji,jj,1) == 1._wp                   )   THEN
      !       rpov   (ji,jj,jk) = 0.5_wp  * (    rpo(ji,jj  ,jk) +    rpo(ji  ,jj+1,jk)   )
      !       r1_rpov(ji,jj,jk) = 1._wp / rpov(ji,jj,jk)
      !    ENDIF
      !    !
      !    IF( ( (       umask(ji,jj-1,1) + umask(ji,jj+1,1) ) < 2._wp   )  .AND. &
      !     &                         tmask(ji,jj,1) == 1._wp                   )   THEN
      !       rpou   (ji,jj,jk) = 0.5_wp  * (    rpo(ji,jj  ,jk) +    rpo(ji  ,jj+1,jk)   )
      !       r1_rpou(ji,jj,jk) = 1._wp / rpou(ji,jj,jk)
      !    ENDIF
      !    !
      !    rpof   (ji,jj,jk) = 0.25_wp * (    rpo(ji,jj+1,jk) +    rpo(ji+1,jj+1,jk)   &
      !        &                         +    rpo(ji,jj  ,jk) +    rpo(ji+1,jj  ,jk)   )
      !    r1_rpof(ji,jj,jk) = 1._wp / rpof(ji,jj,jk)
      ! END_3D
      !
      !! "Partial cells" (not rotated) (exp02)
      ! DO_3D_00_00(1,jpkm1)
      !   rpou   (ji,jj,jk) = MAX( rpo(ji,jj,jk), rpo(ji+1,jj  ,jk))
      !   r1_rpou(ji,jj,jk) = 1._wp / rpou(ji,jj,jk)
      !   rpov   (ji,jj,jk) = MAX( rpo(ji,jj,jk), rpo(ji  ,jj+1,jk) )
      !   r1_rpov(ji,jj,jk) = 1._wp / rpov(ji,jj,jk)
      !   rpof   (ji,jj,jk) = 1._wp
      !   r1_rpof(ji,jj,jk) = 1._wp / rpof(ji,jj,jk)
      ! END_3D
      !
      !
      ! CALL lbc_lnk_multi( 'dommsk', rpo ,  'T', 1._wp,                     &
      !         &                     rpou,  'U', 1._wp,                     &
      !         &                     rpov,  'V', 1._wp,                     &
      !         &                     rpof,  'F', 1._wp, kfillmode=jpfillcopy )
      !
      !         !== SMOOTHING ==!
      ! Shapiro(1/2) on i and j components separatly
      WRITE(numout,*) 'dommsk : shapiro filter applied on rpo (',nn_smo,' times)'
      DO jl = 1, nn_smo
        z3d(:,:,:) = rpo(:,:,:)
        DO jk = 1,jpkm1
          DO jj = 2, jpjm1
             DO ji = 2, jpim1
               IF (tmask(ji-1,jj,jk)*tmask(ji,jj,jk)*tmask(ji+1,jj,jk) == 1._wp) THEN
                 z3d(ji,jj,jk) = 0.5_wp * rpo(ji,jj,jk) + 0.25_wp * (rpo(ji-1,jj,jk) + rpo(ji+1,jj,jk))
               ENDIF
             END DO
          END DO
        END DO
        CALL lbc_lnk( 'dommsk', z3d,  'T', 1._wp, kfillmode=jpfillcopy )
        rpo(:,:,:) = z3d(:,:,:)
        CALL lbc_lnk( 'dommsk', rpo,  'T', 1._wp, kfillmode=jpfillcopy )
        !
        z3d(:,:,:) = rpo(:,:,:)
        DO jk = 1,jpkm1
          DO jj = 2, jpjm1
             DO ji = 2, jpim1
               IF (tmask(ji,jj-1,jk)*tmask(ji,jj,jk)*tmask(ji,jj+1,jk) == 1._wp) THEN
                  z3d(ji,jj,jk) = 0.5_wp * rpo(ji,jj,jk) + 0.25_wp * (rpo(ji,jj-1,jk) + rpo(ji,jj+1,jk))
               ENDIF
             END DO
          END DO
        END DO
        CALL lbc_lnk( 'dommsk', z3d,  'T', 1._wp, kfillmode=jpfillcopy )
        rpo(:,:,:) = z3d(:,:,:)
        !
        WHERE(tmask(:,:,:) == 0._wp)
          rpo   (:,:,:) = rn_abp             ! used on T points
        END WHERE
      END DO
      !
      !! classic penalisation
      ! DO_3D_00_00(1,jpkm1) // END_3D
      ! rpou   (ji,jj,jk) = 0.5_wp  * (    rpo(ji,jj  ,jk) +    rpo(ji+1,jj  ,jk)   )
      ! rpov   (ji,jj,jk) = 0.5_wp  * (    rpo(ji,jj  ,jk) +    rpo(ji  ,jj+1,jk)   )
      ! rpof   (ji,jj,jk) = 0.25_wp * (    rpo(ji,jj+1,jk) +    rpo(ji+1,jj+1,jk)   &
      ! &                         +    rpo(ji,jj  ,jk) +    rpo(ji+1,jj  ,jk)   )
      DO jk = 1,jpk
        DO jj = 1, jpjm1
           DO ji = 1, jpim1
               IF (umask(ji,jj,jk) == 1._wp ) THEN
                 rpou(ji,jj,jk) = 0.5_wp  * (    rpo(ji,jj  ,jk) +    rpo(ji+1,jj  ,jk)   )
               ENDIF
               IF (vmask(ji,jj,jk) == 1._wp ) THEN
                 rpov(ji,jj,jk) = 0.5_wp  * (    rpo(ji,jj  ,jk) +    rpo(ji  ,jj+1,jk)   )
              ENDIF
               IF (fmask(ji,jj,jk) == 1._wp ) THEN
                 rpof(ji,jj,jk) = 0.25_wp * (    rpo(ji,jj+1,jk) +    rpo(ji+1,jj+1,jk)   &
                     &                      +    rpo(ji,jj  ,jk) +    rpo(ji+1,jj  ,jk)   )
              ENDIF
          END DO
       END DO
     END DO
      !
      CALL lbc_lnk_multi( 'dommsk', rpov  , 'V', 1._wp,     &
          &                         rpou  , 'U', 1._wp,     &
          &                         rpof  , 'F', 1._wp,     kfillmode=jpfillcopy )
      !
      !         !== INVERSE ==!
      r1_rpo (:,:,:) = 1._wp
      r1_rpou(:,:,:) = 1._wp
      r1_rpov(:,:,:) = 1._wp
      r1_rpof(:,:,:) = 1._wp
      DO jk = 1,jpk
        DO jj = 1, jpj
           DO ji = 1, jpi
             r1_rpo (ji,jj,jk) = 1._wp / rpo (ji,jj,jk)
             r1_rpou(ji,jj,jk) = 1._wp / rpou(ji,jj,jk)
             r1_rpov(ji,jj,jk) = 1._wp / rpov(ji,jj,jk)
             r1_rpof(ji,jj,jk) = 1._wp / rpof(ji,jj,jk)
           END DO
        END DO
      END DO
      !
      CALL lbc_lnk_multi( 'dommsk', r1_rpo   , 'T', 1._wp,                     &
              &                     r1_rpou  , 'U', 1._wp,                     &
              &                     r1_rpov  , 'V', 1._wp,                     &
              &                     r1_rpof  , 'F', 1._wp, kfillmode=jpfillcopy )
      !

      !
      !
      !  Friction applied on the inner boundary layer
      ! -------------------------------------------------
      bmpu(:,:,:) = 0._wp
      bmpv(:,:,:) = 0._wp
      !
      ! especially for upwellings
      z1d =  rn_dx / REAL(nn_AM98, wp)
      DO_2D_00_00
        ! normal to North and South coast
        ! z1y =    tanh( (gphit(ji,jj) - 0._wp       - z1d*0.5_wp )/(rn_cnp * z1d) ) +       &
        !     &  - tanh( (gphit(ji,jj) - 2000000._wp + z1d*0.5_wp )/(rn_cnp * z1d) ) - 1._wp
        ! z1y = (1._wp - rn_abp) * z1y + rn_abp
        ! z1x =    tanh( (gphit(ji,jj+1) - 0._wp       - z1d*0.5_wp )/(rn_cnp * z1d) ) +       &
        !     &  - tanh( (gphit(ji,jj+1) - 2000000._wp + z1d*0.5_wp )/(rn_cnp * z1d) ) - 1._wp
        ! z1x = (1._wp - rn_abp) * z1x + rn_abp
        ! z1d = 0.5_wp * ( z1y + z1x )
        ! ! bmpv(ji,jj, 1:jpkm1) = rn_fsp * 0.5_wp * ( ff_t(ji,jj) + ff_t(ji  ,jj+1) )   &
        ! !    &                                   * (1._wp - z1d)
        ! bmpv(ji,jj, 1:jpkm1) = rn_fsp * (1._wp - z1d)
        !
        ! normal to the East and West coast
        ! z1y =    tanh( (glamt(ji,jj) - 0._wp    v   - z1d*0.5_wp )/(rn_cnp * z1d) ) +       &
        !     &  - tanh( (glamt(ji,jj) - 2000000._wp + z1d*0.5_wp )/(rn_cnp * z1d) ) - 1._wp
        ! z1y = (1._wp - rn_abp) * z1y + rn_abp
        ! z1x =    tanh( (glamt(ji+1,jj) - 0._wp       - z1d*0.5_wp )/(rn_cnp * z1d) ) +       &
        !     &  - tanh( (glamt(ji+1,jj) - 2000000._wp + z1d*0.5_wp )/(rn_cnp * z1d) ) - 1._wp
        ! z1x = (1._wp - rn_abp) * z1x + rn_abp
        ! z1d = 0.5_wp * ( z1y + z1x )
        !
        ! North and South paralell to the coast (0°)
        ! z1y =    tanh( (gphit(ji,jj) - 0._wp       - z1d*0.5_wp )/(rn_cnp * z1d) ) +       &
        !     &  - tanh( (gphit(ji,jj) - 2000000._wp + z1d*0.5_wp )/(rn_cnp * z1d) ) - 1._wp
        ! z1y = (1._wp - rn_abp) * z1y + rn_abp
        ! !
        ! (0°)
        ! bmpu(ji,jj, 1:jpkm1) = rn_fsp * (1._wp - z1y)
        !(45°)
        ! bmpu(ji,jj, 1:jpkm1) = rn_fsp * (1._wp - z1y) / SQRT(2._wp)
        ! bmpv(ji,jj, 1:jpkm1) = rn_fsp * (1._wp - z1y) / SQRT(2._wp)
        ! bmpu(ji,jj, 1:jpkm1) = rn_fsp * MIN(1._wp - z1y, 1._wp - z1d)
        !
      END_2D
      CALL lbc_lnk_multi( 'dommsk', bmpu,  'U', 1._wp, bmpv  , 'V', 1._wp, kfillmode=jpfillcopy )
      WHERE(vmask(:,:,:) == 0._wp)
        bmpv  (:,:,:) = 0             ! used on V points
      END WHERE
      WHERE(umask(:,:,:) == 0._wp)
        bmpu  (:,:,:) = 0             ! used on U points
      END WHERE
      !
      WRITE(numout,*) 'dommsk : impermeability bmp (sigma) used nn_fsp=',nn_fsp
      SELECT CASE( nn_fsp )           ! == layer drag formulation
      CASE ( -5 )
        WHERE(rpou(:,:,:) > 0.5_wp)
          bmpu  (:,:,:) = rn_fsp           ! used on V points
        END WHERE
        WHERE(rpov(:,:,:) > 0.5_wp)
          bmpv  (:,:,:) = rn_fsp             ! used on U points
        END WHERE
      CASE ( -4 )
        WHERE(rpou(:,:,:) <= 1._wp)
          bmpu  (:,:,:) = rn_fsp * rpou(:,:,:)             ! used on V points
        END WHERE
        WHERE(rpov(:,:,:) <= 1._wp)
          bmpv  (:,:,:) = rn_fsp * rpov(:,:,:)            ! used on U points
        END WHERE
      CASE ( -1 )
        bmpu  (:,:,:) = 0._wp
        bmpv  (:,:,:) = 0._wp
      CASE ( -2 )
        bmpu  (:,:,:) = 0._wp
        bmpv  (:,:,:) = 0._wp
      CASE ( 0 )
        WHERE(rpou(:,:,:) <= 0.5_wp)
          bmpu  (:,:,:) = rn_fsp             ! used on V points
        END WHERE
        WHERE(rpov(:,:,:) <= 0.5_wp)
          bmpv  (:,:,:) = rn_fsp             ! used on U points
        END WHERE
      CASE ( 1 )
        WHERE(rpou(:,:,:) <= 1._wp)
          bmpu  (:,:,:) = rn_fsp * r1_rpou(:,:,:)             ! used on V points
        END WHERE
        WHERE(rpov(:,:,:) <= 1._wp)
          bmpv  (:,:,:) = rn_fsp * r1_rpov(:,:,:)
        END WHERE
      CASE ( 2 )
        WHERE(rpou(:,:,:) <= 1._wp)
          bmpu  (:,:,:) = rn_fsp * SQRT(r1_rpou(:,:,:))             ! used on V points
        END WHERE
        WHERE(rpov(:,:,:) <= 1._wp)
          bmpv  (:,:,:) = rn_fsp * SQRT(r1_rpov(:,:,:))
        END WHERE
      CASE DEFAULT                                             ! error
         CALL ctl_stop('STOP','dommsk: wrong value for nn_fsp'  )
      END SELECT
      !
      ! DO_2D_00_00
      !    bmpu(ji,jj, 1:jpkm1) = rn_fsp * 0.5_wp * ( ff_t(ji,jj) + ff_t(ji+1,jj  ) )   &
      !       &                          * (1._wp - rpou(ji,jj,1:jpkm1)) * r1_rpou(ji,jj,1:jpkm1)
      !    bmpv(ji,jj, 1:jpkm1) = rn_fsp * 0.5_wp * ( ff_t(ji,jj) + ff_t(ji  ,jj+1) )   &
      !       &                          * (1._wp - rpov(ji,jj,1:jpkm1)) * r1_rpov(ji,jj,1:jpkm1)
      ! END_2D

      ! DO_2D_00_00
      !    bmpu(ji,jj, 1:jpkm1) = rn_fsp * ( 1._wp - rpou(ji,jj,1:jpkm1) )
      !    bmpv(ji,jj, 1:jpkm1) = rn_fsp * ( 1._wp - rpov(ji,jj,1:jpkm1) )
      ! END_2D
      !
      ! constant
      ! DO_2D_00_00
      !    IF( rpou(ji,jj,1) < 1._wp )   bmpu(ji,jj, 1:jpkm1) = rn_fsp
      !    IF( rpov(ji,jj,1) < 1._wp )   bmpv(ji,jj, 1:jpkm1) = rn_fsp
      ! END_2D
      CALL lbc_lnk_multi( 'dommsk', bmpu,  'U', 1._wp, bmpv  , 'V', 1._wp, kfillmode=jpfillcopy )
#endif

# if defined key_bvp_bath
      !  Add bathymetry from the porosity field
      ! --------------------------------------
      ! H = 5000m ; zalpha = (b-a)/2 ; zbeta = (b+a)/2
      batht(:,:) = 4500._wp ! meters
      !
      z1d =  rn_dx / REAL(nn_AM98, wp) ; ze1 = 0._wp
      DO jk = 1, jpkm1
        DO_2D_00_00
          !  tanh penalisation (T point) coast = middle of rpo
          ! ----------------------------------------
          z1x =    tanh( (glamt(ji,jj) - 0._wp       )/(rn_cnp * z1d * 0.25_wp) ) +       &
            &  - tanh( (glamt(ji,jj) - 2000000._wp )/(rn_cnp * z1d * 0.25_wp) )
          z1x = z1x * 0.5_wp
          !
          z1y =    tanh( (gphit(ji,jj) - 0._wp       )/(rn_cnp * z1d * 0.25_wp) ) +       &
            &  - tanh( (gphit(ji,jj) - 2000000._wp )/(rn_cnp * z1d * 0.25_wp) )
          z1y = z1y * 0.5_wp
          !
          !  linear penalisation (speed points) coast = middle of porosity
          ! -------------------------------------------------------------
          !
          ! z1x1 =  (1._wp - rn_abp)*(glamt(ji,jj) -       0._wp - (rn_cnp*ze1/4._wp))/(rn_cnp * z1d) + 0.5_wp*(1._wp + rn_abp)
          ! z1x1 = MAX(MIN(z1x1,1._wp), rn_abp)
          ! z1x2 = -(1._wp - rn_abp)*(glamt(ji,jj) - 2000000._wp + (rn_cnp*ze1/4._wp))/(rn_cnp * z1d) + 0.5_wp*(1._wp + rn_abp)
          ! z1x2 = MAX(MIN(z1x2,1._wp), rn_abp)
          ! z1x  = MIN(z1x1,z1x2)
          ! !
          ! z1y1 =  (1._wp - rn_abp)*(gphit(ji,jj) -       0._wp - (rn_cnp*ze1/4._wp))/(rn_cnp * z1d) + 0.5_wp*(1._wp + rn_abp)
          ! z1y1 = MAX(MIN(z1y1,1._wp), rn_abp)
          ! z1y2 = -(1._wp - rn_abp)*(gphit(ji,jj) - 2000000._wp + (rn_cnp*ze1/4._wp))/(rn_cnp * z1d) + 0.5_wp*(1._wp + rn_abp)
          ! z1y2 = MAX(MIN(z1y2,1._wp), rn_abp)
          ! z1y  = MIN(z1y1,z1y2)
          !
          !
          ! W
          ! WRITE(numout,*) 'dommsk : West coast bathy'
          ! z1x1 =  (1._wp - rn_abp)*(glamt(ji,jj) -       0._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
          ! z1x1 = MAX(MIN(z1x1,1._wp), rn_abp)
          !
          ! ! EW
          ! WRITE(numout,*) 'dommsk : West-East coast penalised'
          ! z1x1 =  (1._wp - rn_abp)*(glamt(ji,jj) -       0._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
          ! z1x1 = MAX(MIN(z1x1,1._wp), rn_abp)
          ! z1x2 = -(1._wp - rn_abp)*(glamt(ji,jj) - 2000000._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
          ! z1x2 = MAX(MIN(z1x2,1._wp), rn_abp)
          ! z1x  = MIN(z1x1,z1x2)
          !
          !NS
          ! WRITE(numout,*) 'dommsk : North-South coast penalised'
          ! z1y1 =  (1._wp - rn_abp)*(gphit(ji,jj) -       0._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
          ! z1y1 = MAX(MIN(z1y1,1._wp), rn_abp)
          ! z1y2 = -(1._wp - rn_abp)*(gphit(ji,jj) - 2000000._wp)/(REAL(rn_cnp,wp) * z1d) + 0.5_wp*(1._wp + rn_abp)
          ! z1y2 = MAX(MIN(z1y2,1._wp), rn_abp)
          ! z1y  = MIN(z1y1,z1y2)
          !
          !
          !
          ! ze1 = z1x1
          ze1 = MIN(z1x,z1y)
          batht(ji,jj) = 4500._wp * ( 1 - ze1 ) + ze1 * 5000._wp
          !
        END_2D
      END DO
      !
      CALL lbc_lnk( 'dommsk', batht,  'T', 1._wp)
      !
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
             IF (umask(ji,jj,1) == 1._wp ) THEN
               bathu(ji,jj) = 0.5_wp  * (    batht(ji,jj  ) +    batht(ji+1,jj  )   )
             ENDIF
             IF (vmask(ji,jj,1) == 1._wp ) THEN
               bathv(ji,jj) = 0.5_wp  * (    batht(ji,jj  ) +    batht(ji  ,jj+1)   )
            ENDIF
        END DO
     END DO
      !
      CALL lbc_lnk_multi( 'dommsk', bathv  , 'V', 1._wp,     &
          &                         bathu  , 'U', 1._wp,     kfillmode=jpfillcopy )
#endif

  IF (ln_hdiv_AD) THEN
    ! ONLY TRUE IN FREESLIP
      write(numout,*) 'dommsk : r1_e1e2t half at the coast (ln_hdiv_AD=T)'
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
            IF  ( (     fmask(ji-1,jj  ,1) + fmask(ji  ,jj  ,1)  &
                  &   + fmask(ji-1,jj-1,1) + fmask(ji  ,jj-1,1)  ) == 1._wp) THEN
                 r1_e1e2t(ji,jj) = 2._wp*r1_e1e2t(ji,jj)   ! dont touch e1e2t otherwise e3u,v will be badly impacted
            ENDIF
          END DO
       END DO
       CALL lbc_lnk( 'dommsk', r1_e1e2t  , 'T', 1._wp )
  ENDIF

      ! Interior domain mask  (used for global sum)
      ! --------------------
      !
      iif = nn_hls   ;   iil = nlci - nn_hls + 1
      ijf = nn_hls   ;   ijl = nlcj - nn_hls + 1
      !
      !                          ! halo mask : 0 on the halo and 1 elsewhere
      tmask_h(:,:) = 1._wp
      tmask_h( 1 :iif,   :   ) = 0._wp      ! first columns
      tmask_h(iil:jpi,   :   ) = 0._wp      ! last  columns (including mpp extra columns)
      tmask_h(   :   , 1 :ijf) = 0._wp      ! first rows
      tmask_h(   :   ,ijl:jpj) = 0._wp      ! last  rows (including mpp extra rows)
      !
      !                          ! north fold mask
      tpol(1:jpiglo) = 1._wp
      fpol(1:jpiglo) = 1._wp
      IF( jperio == 3 .OR. jperio == 4 ) THEN      ! T-point pivot
         tpol(jpiglo/2+1:jpiglo) = 0._wp
         fpol(     1    :jpiglo) = 0._wp
         IF( mjg(nlej) == jpjglo ) THEN                  ! only half of the nlcj-1 row for tmask_h
            DO ji = iif+1, iil-1
               tmask_h(ji,nlej-1) = tmask_h(ji,nlej-1) * tpol(mig(ji))
            END DO
         ENDIF
      ENDIF
      !
      IF( jperio == 5 .OR. jperio == 6 ) THEN      ! F-point pivot
         tpol(     1    :jpiglo) = 0._wp
         fpol(jpiglo/2+1:jpiglo) = 0._wp
      ENDIF
      !
      !                          ! interior mask : 2D ocean mask x halo mask
      tmask_i(:,:) = ssmask(:,:) * tmask_h(:,:)


      ! Lateral boundary conditions on velocity (modify fmask)
      ! ---------------------------------------
      IF( rn_shlat /= 0 ) THEN      ! Not free-slip lateral boundary condition
         !
         ALLOCATE( zwf(jpi,jpj) )
         !
         DO jk = 1, jpk
            zwf(:,:) = fmask(:,:,jk)
            DO_2D_00_00
               IF( fmask(ji,jj,jk) == 0._wp ) THEN
                  fmask(ji,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
                     &                                           zwf(ji-1,jj), zwf(ji,jj-1)  )  )
               ENDIF
            END_2D
            DO jj = 2, jpjm1
               IF( fmask(1,jj,jk) == 0._wp ) THEN
                  fmask(1  ,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(2,jj), zwf(1,jj+1), zwf(1,jj-1) ) )
               ENDIF
               IF( fmask(jpi,jj,jk) == 0._wp ) THEN
                  fmask(jpi,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(jpi,jj+1), zwf(jpim1,jj), zwf(jpi,jj-1) ) )
               ENDIF
            END DO
            DO ji = 2, jpim1
               IF( fmask(ji,1,jk) == 0._wp ) THEN
                  fmask(ji, 1 ,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,1), zwf(ji,2), zwf(ji-1,1) ) )
               ENDIF
               IF( fmask(ji,jpj,jk) == 0._wp ) THEN
                  fmask(ji,jpj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,jpj), zwf(ji-1,jpj), zwf(ji,jpjm1) ) )
               ENDIF
            END DO
#if defined key_agrif
            IF( .NOT. AGRIF_Root() ) THEN
               IF ((nbondi ==  1).OR.(nbondi == 2)) fmask(nlci-1 , :     ,jk) = 0.e0      ! east
               IF ((nbondi == -1).OR.(nbondi == 2)) fmask(1      , :     ,jk) = 0.e0      ! west
               IF ((nbondj ==  1).OR.(nbondj == 2)) fmask(:      ,nlcj-1 ,jk) = 0.e0      ! north
               IF ((nbondj == -1).OR.(nbondj == 2)) fmask(:      ,1      ,jk) = 0.e0      ! south
            ENDIF
#endif
         END DO
         !
         DEALLOCATE( zwf )
         !
         CALL lbc_lnk( 'dommsk', fmask, 'F', 1._wp )      ! Lateral boundary conditions on fmask
         !
         ! CAUTION : The fmask may be further modified in dyn_vor_init ( dynvor.F90 ) depending on ln_vorlat
         !
      ENDIF

      ! User defined alteration of fmask (use to reduce ocean transport in specified straits)
      ! --------------------------------
      !
      CALL usr_def_fmask( cn_cfg, nn_cfg, fmask )
      !
   END SUBROUTINE dom_msk

   !!======================================================================
END MODULE dommsk
