MODULE usrdef_zgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_zgr  ***
   !!
   !!                       ===  AM98 configuration  ===
   !!
   !! User defined : vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-06  (G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system
   !!      zgr_z      : reference 1D z-coordinate
   !!      zgr_top_bot: ocean top and bottom level indices
   !!      zgr_zco    : 3D verticl coordinate in pure z-coordinate case
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce        ! ocean domain
   USE depth_e3       ! depth <=> e3
   USE usrdef_nam
   USE phycst         , ONLY : rad      ! physical constants for rad

   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr        ! called by domzgr.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_zgr( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
      &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,    &   ! 1D reference vertical coordinate
      &                    pdept , pdepw ,                             &   ! 3D t & w-points depth
      &                    pe3t  , pe3u  , pe3v   , pe3f ,             &   ! vertical scale factors
      &                    pe3w  , pe3uw , pe3vw         ,             &   !     -      -      -
      &                    k_top , k_bot                 )                 ! top & bottom ocean level
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(out) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags
      LOGICAL                   , INTENT(out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pdept, pdepw                ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3w , pe3uw, pe3vw         ! i-scale factors
      INTEGER , DIMENSION(:,:)  , INTENT(out) ::   k_top, k_bot                ! first & last ocean level
      !
      INTEGER  ::   inum   ! local logical unit
      REAL(WP) ::   z_zco, z_zps, z_sco, z_cav
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : AM98 configuration (z-coordinate closed flat box ocean without cavities)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      ! type of vertical coordinate
      ! ---------------------------
      ld_zco    = .FALSE.         ! AM98 case:  z-coordinate without ocean cavities
      ld_zps    = .FALSE.
      ld_sco    = .TRUE.
      ld_isfcav = .FALSE.
      !
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      CALL zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
      CALL zgr_msk_top_bot( k_top , k_bot)                 ! masked top and bottom ocean t-level indices
      !
      !                                                     ! z-coordinate (3D arrays) from the 1D z-coord.
      CALL zgr_zco( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   ! in  : 1D reference vertical coordinate
         &          pdept   , pdepw   ,                     &   ! out : 3D t & w-points depth
         &          pe3t    , pe3u    , pe3v   , pe3f   ,   &   !       vertical scale factors
         &          pe3w    , pe3uw   , pe3vw             )     !           -      -      -
      !
   END SUBROUTINE usr_def_zgr


   SUBROUTINE zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! 1D reference vertical coordinate
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z  ***
      !!
      !! ** Purpose :   set the 1D depth of model levels and the resulting
      !!              vertical scale factors.
      !!
      !! ** Method  :   1D z-coordinate system (use in all type of coordinate)
      !!       The depth of model levels is set from dep(k), an analytical function:
      !!                   w-level: depw_1d  = dep(k)
      !!                   t-level: dept_1d  = dep(k+0.5)
      !!       The scale factors are the discrete derivative of the depth:
      !!                   e3w_1d(jk) = dk[ dept_1d ]
      !!                   e3t_1d(jk) = dk[ depw_1d ]
      !!           with at top and bottom :
      !!                   e3w_1d( 1 ) = 2 * ( dept_1d( 1 ) - depw_1d( 1 ) )
      !!                   e3t_1d(jpk) = 2 * ( dept_1d(jpk) - depw_1d(jpk) )
      !!       The depth are then re-computed from the sum of e3. This ensures
      !!    that depths are identical when reading domain configuration file.
      !!    Indeed, only e3. are saved in this file, depth are compute by a call
      !!    to the e3_to_depth subroutine.
      !!
      !!       Here the Madec & Imbard (1996) function is used.
      !!
      !! ** Action  : - pdept_1d, pdepw_1d : depth of T- and W-point (m)
      !!              - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!
      !! Reference : Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
      !!             Madec and Imbard, 1996, Clim. Dyn.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d   ! 1D grid-point depth        [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d    ! 1D vertical scale factors  [m]
      !
      INTEGER  ::   jk       ! dummy loop indices
        !!----------------------------------------------------------------------
      !
        !
      IF(lwp) THEN            ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates '
         WRITE(numout,*) '    ~~~~~~~'
      ENDIF

      !
      ! 1D Reference z-coordinate    (using Madec & Imbard 1996 function)
      ! -------------------------
      !
      ! depth at T and W-points
      pdepw_1d(1) =   0._wp
      pdept_1d(1) = 250._wp
      !
      pdepw_1d(2) = 500._wp
      pdept_1d(2) = 750._wp
      !
      !                       ! e3t and e3w from depth
      CALL depth_to_e3( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d )
      !
      !                       ! recompute depths from SUM(e3)  <== needed
      CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d )
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
   END SUBROUTINE zgr_z




   SUBROUTINE zgr_msk_top_bot( k_top , k_bot)
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_msk_top_bot  ***
      !!
      !! ** Purpose :   set the masked top and bottom ocean t-levels
      !!
      !! ** Method  :   AM98 case = closed flat box ocean without ocean cavities
      !!                   k_top = 1     except along north, south, east and west boundaries
      !!                   k_bot = jpk-1 except along north, south, east and west boundaries
      !!
      !! ** Action  : - k_top : first wet ocean level index
      !!              - k_bot : last  wet ocean level index
      !!----------------------------------------------------------------------
      INTEGER , DIMENSION(:,:), INTENT(out) ::   k_top , k_bot   ! first & last wet ocean level
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D local workspace
      INTEGER  ::   ji, jj, jzo                    ! dummy loop indices
      REAL(wp) ::   zylim0, zylim1, zxlim0, zxlim1, ze1, zex, zey ! limit of the domain [m]
      REAL(wp) ::   zcoeff, ze0, zey_do, zey_up, zex_le, zex_ri    ! local scalar
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_top_bot : defines the top and bottom wet ocean levels.'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '       AM98 case : closed flat box ocean without ocean cavities'
      !
      z2d(:,:) = REAL( jpkm1 , wp )          ! flat bottom
      !
      CALL lbc_lnk( 'usrdef_zgr', z2d, 'T', 1. )           ! set surrounding land to zero (here jperio=0 ==>> closed)
      !
# if defined key_bvp_bath
      WRITE(numout,*) 'usrdef_zgr : (bath) mask enlarged of nn_gc'
      IF ( lwp ) WRITE(numout,*) 'usrdef_zgr : when nn_gc!=0, the basin is widened on each side of nn_gc (rn_cnp/2<nn_gc)     '
      ze1 = REAL(nn_gc,wp) * ( rn_dx / REAL(nn_AM98, wp) )      ! [m] gridspacing used
      !
      zylim0 =       0._wp - ze1   ! +0km
      zylim1 = 2000000._wp + ze1   ! 2000km
      zxlim0 =       0._wp - ze1   ! +0km
      zxlim1 = 2000000._wp + ze1   ! 2000km
#else
      IF ( ln_mir_adv .OR. ln_mir_coriolis ) THEN
        WRITE(numout,*) 'usrdef_zgr (ln_mir_adv or ln_mir_coriolis =T) bassin 45° squared'
        ze1 =  rn_dx / REAL(nn_AM98, wp)      ! [m] gridspacing used
        ! at 45° T points are centers sqrt(2) = 1.414
        ! des coins 3->2 et 2->(1)->3
        IF      (nn_AM98==4) THEN
          zylim0 =       0._wp - ze1 ! S
          zylim1 = 2000000._wp ! N
          zxlim0 =       0._wp ! W
          zxlim1 = 2000000._wp + ze1 ! E
        ELSE IF (nn_AM98 == 16) THEN ! no need apparently
          zylim0 =       0._wp ! S
          zylim1 = 2000000._wp ! N
          zxlim0 =       0._wp ! W
          zxlim1 = 2000000._wp ! E
        ENDIF
      ELSE
        zylim0 =       0._wp
        zylim1 = 2000000._wp
        zxlim0 =       0._wp
        zxlim1 = 2000000._wp
      ENDIF
#endif
      !
      IF ( ln_obstacle ) THEN
        WRITE(numout,*) 'usrdef_zgr : Obstacle kept at 1°/4 resolution (ln_obstacle=T)'
        WRITE(numout,*) '             solved at 1°/',nn_AM98,' resolution           '
          DO jj = 1, jpj
             DO ji = 1, jpi
              IF ( k_top0(ji,jj) == 1 )  THEN
                k_top(ji,jj) = 1    ! = ocean
                k_bot(ji,jj) = NINT( z2d(ji,jj) )
              ELSE
                k_top(ji,jj) = 0    ! = land
                k_bot(ji,jj) = 0
              END IF
             END DO
          END DO
      ELSE
          WRITE(numout,*) 'usrdef_zgr : Solved at 1°/',nn_AM98,' resolution (ln_obstacle=F)'
          ! for limit
          ze1 =  0.1_wp * rn_dx / REAL(nn_AM98, wp)                   ! [m] gridspacing used -10%
          zylim0 = zylim0 - ze1
          zylim1 = zylim1 + ze1
          zxlim0 = zxlim0 - ze1
          zxlim1 = zxlim1 + ze1
          !
          k_top(:,:) = 0    ! = land
          k_bot(:,:) = 0
          !
          DO jj = 1, jpj
             DO ji = 1, jpi
               ! au bit près peut faire des glitches, mieux vaut rallonger de 10%
               IF ( gphit(ji,jj) > zylim0 .AND. gphit(ji,jj) < zylim1 .AND. &
                  & glamt(ji,jj) > zxlim0 .AND. glamt(ji,jj) < zxlim1       )  THEN
                 ! on élargit de 10% les frontières, et V et U à droite et en haut de T (0.9 = -0.1 + 1) (0.1 = +0.1)
                 ! IF ( gphiv(ji,jj) > (zylim0 + 0.9_wp*zey_do) .AND. gphiv(ji,jj) < (zylim1 + 0.1_wp*zey_up) .AND. &
                 !    & glamu(ji,jj) > (zxlim0 + 0.9_wp*zex_le) .AND. glamu(ji,jj) < (zxlim1 + 0.1_wp*zex_ri)       )  THEN
                 ! pour 45° on élargit (-0.1 = -0.1 + 0)
                 ! IF ( gphiv(ji,jj) > (zylim0 - 0.1_wp*zey) .AND. gphiv(ji,jj) < (zylim1 + 1.1_wp*zey) .AND. &
                 !    & glamu(ji,jj) > (zxlim0 - 0.1_wp*zex) .AND. glamu(ji,jj) < (zxlim1 + 1.1_wp*zex)       )  THEN
                 k_top(ji,jj) = 1    ! = ocean (T point)
                 k_bot(ji,jj) = NINT( z2d(ji,jj) )
               END IF
             END DO
          END DO
          ! mask the lonely corners
          DO jj = 2, jpjm1
             DO ji = 2, jpim1
             zcoeff = k_top(ji+1,jj) + k_top(ji,jj+1)   &
                +     k_top(ji-1,jj) + k_top(ji,jj-1)
             IF ( zcoeff <= 1._wp )   THEN
                k_top(ji,jj) = 0    ! = land
                k_bot(ji,jj) = 0
             END IF
             END DO
          END DO
      ENDIF
      !
      !
   END SUBROUTINE zgr_msk_top_bot


   SUBROUTINE zgr_zco( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   ! in : 1D reference vertical coordinate
      &                pdept   , pdepw   ,                     &   ! out: 3D t & w-points depth
      &                pe3t    , pe3u    , pe3v   , pe3f   ,   &   !      vertical scale factors
      &                pe3w    , pe3uw   , pe3vw             )     !          -      -      -
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_zco  ***
      !!
      !! ** Purpose :   define the reference z-coordinate system
      !!
      !! ** Method  :   set 3D coord. arrays to reference 1D array
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth       [m]
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   pe3t_1d , pe3w_1d           ! 1D vertical scale factors [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! grid-point depth          [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors    [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         !    -       -      -
      !
      INTEGER  ::   jk
      !!----------------------------------------------------------------------
      !
      DO jk = 1, jpk
         pdept(:,:,jk) = pdept_1d(jk)
         pdepw(:,:,jk) = pdepw_1d(jk)
         pe3t (:,:,jk) = pe3t_1d (jk)
         pe3u (:,:,jk) = pe3t_1d (jk)
         pe3v (:,:,jk) = pe3t_1d (jk)
         pe3f (:,:,jk) = pe3t_1d (jk)
         pe3w (:,:,jk) = pe3w_1d (jk)
         pe3uw(:,:,jk) = pe3w_1d (jk)
         pe3vw(:,:,jk) = pe3w_1d (jk)
      END DO
      !
   END SUBROUTINE zgr_zco

   !!======================================================================
END MODULE usrdef_zgr
