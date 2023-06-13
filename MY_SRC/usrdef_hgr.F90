MODULE usrdef_hgr
   !!======================================================================
   !!                     ***  MODULE usrdef_hgr   ***
   !!
   !!                     ===  AM98 configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  4.0 ! 2016-03  (S. Flavoni)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr   : initialize the horizontal mesh
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp, njmpp       ! ocean space and time domain
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE usrdef_nam     !
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_hgr   ! called in domhgr.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_hgr.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_hgr( plamt , plamu , plamv  , plamf  ,   &   ! geographic position (required)
      &                    pphit , pphiu , pphiv  , pphif  ,   &   !
      &                    kff   , pff_f , pff_t  ,            &   ! Coriolis parameter  (if domain not on the sphere)
      &                    pe1t  , pe1u  , pe1v   , pe1f   ,   &   ! scale factors       (required)
      &                    pe2t  , pe2u  , pe2v   , pe2f   ,   &   !
      &                    ke1e2u_v      , pe1e2u , pe1e2v ,   &   ! u- & v-surfaces (if gridsize reduction is used in strait(s))
      &                    plamt0, pphit0, pk_top0              )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE usr_def_hgr  ***
      !!
      !! ** Purpose :   user defined mesh and Coriolis parameter
      !!
      !! ** Method  :   set all intent(out) argument to a proper value
      !!
      !!                Here AM98 configuration :
      !!          Rectangular mid-latitude domain
      !!          - with axes rotated by 45 degrees
      !!          - a constant horizontal resolution of 106 km
      !!          - on a beta-plane
      !!
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in degrees)
      !!              - define coriolis parameter at f-point if the domain in not on the sphere (on beta-plane)
      !!              - define i- & j-scale factors at t-, u-, v- and f-points (in meters)
      !!              - define u- & v-surfaces (if gridsize reduction is used in some straits) (in m2)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   plamt, plamu, plamv, plamf   ! longitude outputs                     [degrees]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pphit, pphiu, pphiv, pphif   ! latitude outputs                      [degrees]
      INTEGER                 , INTENT(out) ::   kff                          ! =1 Coriolis parameter computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pff_f, pff_t                 ! Coriolis factor at f-point                [1/s]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1t, pe1u, pe1v, pe1f       ! i-scale factors                             [m]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe2t, pe2u, pe2v, pe2f       ! j-scale factors                             [m]
      INTEGER                 , INTENT(out) ::   ke1e2u_v                     ! =1 u- & v-surfaces computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1e2u, pe1e2v               ! u- & v-surfaces (if reduction in strait)   [m2]
      !
      INTEGER  ::   ji, jj, jk, ji0, jj0, jk0, jei0               ! dummy loop indices
      REAL(wp) ::   zlam1, zlam0, zcos_theta, zim1 , zjm1 , ze1, ze0  , ze1deg ! local scalars
      REAL(wp) ::   zphi1, zphi0, zsin_theta, zim05, zjm05, znorme        !   -      -
      REAL(wp) ::   zgl, zbl, zgl0, zbl0       !   -      -
      REAL(wp) ::   ze2, zelrefine, ztm, ztmp, zel
      INTEGER  ::   jflag
      REAL(wp) ::   z1x, z1y, z1x1, z1x2, z1y1, z1y2 ! local real
      !
      REAL(wp) ::   zylim0, zylim1, zxlim0, zxlim1, zex, zey ! limit of the domain [m]
      REAL(wp) ::   zcoeff      ! local scalar

      REAL(wp), DIMENSION(:,:), INTENT(out) ::   plamt0, pphit0, pk_top0  ! longitude / latitude bigger grid                     [degrees]
      !!-------------------------------------------------------------------------------
      !
      !     !==  beta-plane with regular grid-spacing and rotated domain ==!  (AM98 configuration)
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : AM98 configuration (beta-plane with rotated regular grid-spacing)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      zcos_theta = COS( rn_theta * rad)
      zsin_theta = SIN( rn_theta * rad)

      !
      IF (ln_obstacle ) THEN
        ze0 =  rn_dx / 4._wp   ! obstacle are at 1°/4
        zgl0 = rn_domsiz + (4._wp + 2._wp * REAL(nn_gc, wp) ) * ze0       ! [m] length of the square with ghostcells
        zbl0 = zgl0 * ( zcos_theta + zsin_theta )   ! length side bigger domain [m]
      ENDIF
      !
      !                       !==  grid point position  ==!
      !
      ze1 =  rn_dx / REAL(nn_AM98, wp)                   ! [m] gridspacing used
      zgl =  rn_domsiz + (4._wp + 2._wp * REAL(nn_gc, wp) ) * ze1   ! [m] length of the square with ghostcells
      ! fit the best square around the square + ghost cells
      zbl = zgl * ( COS( rn_theta * rad ) + SIN( rn_theta * rad ) )   ! length side bigger domain [m]
      !
      ! Translation vers le coin bas-gauche du carré tourné
      zlam1 =  zbl * COS((rn_theta + 45. )* rad ) / SQRT( 2._wp ) - rn_domsiz/2._wp
      zphi1 =  zbl * SIN((rn_theta + 45. )* rad ) / SQRT( 2._wp ) - rn_domsiz/2._wp
      ! select the nearest integer coordonate point
      zlam0 = REAL( anint( zlam1 / (ze1 * zcos_theta)  ), wp ) * ze1 * zcos_theta
      zphi0 = REAl( anint( zphi1 / (ze1 * zcos_theta)  ), wp ) * ze1 * zcos_theta
      !
      IF (ln_obstacle ) THEN
        ! exact origin in meters
        zlam1 =  zbl0 * COS((rn_theta + 45 )* rad ) / SQRT( 2._wp )  - rn_domsiz/2._wp
        zphi1 =  zbl0 * SIN((rn_theta + 45 )* rad ) / SQRT( 2._wp )  - rn_domsiz/2._wp
        ! origin put in the true corner of a cell so there will be no cropping
        ! of the edge cells
        zlam0 = REAL( anint( zlam1 / (ze0 * zcos_theta)  ), wp ) * ze0 * zcos_theta
        zphi0 = REAl( anint( zphi1 / (ze0 * zcos_theta)  ), wp ) * ze0 * zcos_theta
      ENDIF
      ! zlam0 = REAL( anint( zlam1 / ze1 ), wp ) * ze1
      ! zphi0 = REAl( anint( zphi1 / ze1 ), wp ) * ze1
      !
      IF(lwp) WRITE(numout,*) '                           origin position    zlam0   = ', zlam0/1000,   ' km'
      IF(lwp) WRITE(numout,*) '                           origin position    zphi0   = ', zphi0/1000,   ' km'
      ! O1M = OM x rotation_theta - OO1
      ! zim1, zim05, zjm1, zjm05 fit for 2 ghost cells on each side
      DO jj = 1, jpj
         DO ji = 1, jpi
            zim1 = REAL( ji + nimpp - 1 )   ;   zim05 = REAL( ji + nimpp - 1 ) - 0.5
            zjm1 = REAL( jj + njmpp - 1 )   ;   zjm05 = REAL( jj + njmpp - 1 ) - 0.5
            !
        !
        !glamt(i,j) position (meters) at T-point
        !gphit(i,j) position (meters) at T-point
        plamt(ji,jj) =  zim05 * ze1 * zcos_theta - zjm05 * ze1 * zsin_theta - zlam0
        pphit(ji,jj) = +zim05 * ze1 * zsin_theta + zjm05 * ze1 * zcos_theta - zphi0
        !
        !glamu(i,j) position (meters) at U-point
        !gphiu(i,j) position (meters) at U-point
        plamu(ji,jj) =  zim1  * ze1 * zcos_theta - zjm05 * ze1 * zsin_theta - zlam0
        pphiu(ji,jj) = +zim1  * ze1 * zsin_theta + zjm05 * ze1 * zcos_theta - zphi0
        !
        !glamv(i,j) position (meters) at V-point
        !gphiv(i,j) position (meters) at V-point
        plamv(ji,jj) =  zim05 * ze1 * zcos_theta - zjm1  * ze1 * zsin_theta - zlam0
        pphiv(ji,jj) = +zim05 * ze1 * zsin_theta + zjm1  * ze1 * zcos_theta - zphi0
        !
        !glamf(i,j) position (meters) at F-point
        !gphif(i,j) position (meters) at F-point
        plamf(ji,jj) =  zim1  * ze1 * zcos_theta - zjm1  * ze1 * zsin_theta - zlam0
        pphif(ji,jj) = +zim1  * ze1 * zsin_theta + zjm1  * ze1 * zcos_theta - zphi0
     END DO
  END DO
      !
      !
      IF (ln_obstacle ) THEN
        DO jj = 1, jpj
           DO ji = 1, jpi
              ! nn*nn-1 trop d'informations mais parallélisable
              ! ji : indice local dans l'array partionné
              ! ji+nimpp-1 : indice global dans le domaine
              ji0 = INT((ji+nimpp-1)*4/nn_AM98) ; jj0 = INT((jj+njmpp-1)*4/nn_AM98)
              ! IF(lwp) WRITE(numout,*) 'ji0=',ji0,'jj0=',jj0
              zim05 = REAL(ji0, wp) - 0.5 ; zjm05 = REAL(jj0, wp) - 0.5
              !
              plamt0(ji,jj) =   zim05 * ze0 * zcos_theta - zjm05 * ze0 * zsin_theta - zlam0
              pphit0(ji,jj) = + zim05 * ze0 * zsin_theta + zjm05 * ze0 * zcos_theta - zphi0
            END DO
         END DO
         !
         !! False ktop - for the corner (especially for the lonely corners)
         !ze0 =  rn_dx / 4._wp                   ! [m] gridspacing used -10%
         ! zex = ze0 * COS( rn_theta * rad) ; zey = ze0 * COS( rn_theta * rad)
         ! !
         ! zylim0 =       0._wp + 0.9_wp*zey
         ! zylim1 = 2000000._wp + 0.1_wp*zey
         ! zxlim0 =       0._wp + 0.9_wp*zex
         ! zxlim1 = 2000000._wp + 0.1_wp*zex
         ! !
         ! for limit
         zylim0 =       0._wp - 0.1_wp * ze0
         zylim1 = 2000000._wp + 0.1_wp * ze0
         zxlim0 =       0._wp - 0.1_wp * ze0
         zxlim1 = 2000000._wp + 0.1_wp * ze0
         !
         pk_top0(:,:) = 0
         !
         WRITE(numout,*) 'usrdef_hgr : lonely 1°/4 corners (ln_obstacle=T) removed'
           ! for limit
           DO jj = 1, jpj
              DO ji = 1, jpi
               IF ( pphit0(ji,jj) > zylim0 .AND. pphit0(ji,jj) < zylim1 .AND. &
                  & plamt0(ji,jj) > zxlim0 .AND. plamt0(ji,jj) < zxlim1       )  THEN
                 pk_top0(ji,jj) = 1    ! = ocean
               ELSE
                 pk_top0(ji,jj) = 0    ! = land
               END IF
              END DO
           END DO
           ! mask the lonely corners
           ! problematic in usrdef_zgr because done in parallel (cannot use leap indices)
           zex = ze0 * COS( rn_theta * rad) ; zey = ze0 * COS( rn_theta * rad)
           DO jj = 1, jpj
              DO ji = 1, jpi
                ! corners are the closest to both boundaries
                IF ( ( SIGN(1.,pphit0(ji,jj) + zey - zylim0)*SIGN(1.,pphit0(ji,jj) - zey - zylim0) < 0. .AND.    &
                   &   SIGN(1.,plamt0(ji,jj) + zex - zxlim0)*SIGN(1.,plamt0(ji,jj) - zex - zxlim0) < 0. )   .OR. &
                   !
                     ( SIGN(1.,pphit0(ji,jj) + zey - zylim1)*SIGN(1.,pphit0(ji,jj) - zey - zylim1) < 0. .AND.    &
                   &   SIGN(1.,plamt0(ji,jj) + zex - zxlim0)*SIGN(1.,plamt0(ji,jj) - zex - zxlim0) < 0. )   .OR. &
                   !
                     ( SIGN(1.,pphit0(ji,jj) + zey - zylim0)*SIGN(1.,pphit0(ji,jj) - zey - zylim0) < 0. .AND.    &
                   &   SIGN(1.,plamt0(ji,jj) + zex - zxlim1)*SIGN(1.,plamt0(ji,jj) - zex - zxlim1) < 0. )   .OR. &
                   !
                     ( SIGN(1.,pphit0(ji,jj) + zey - zylim1)*SIGN(1.,pphit0(ji,jj) - zey - zylim1) < 0. .AND.    &
                   &   SIGN(1.,plamt0(ji,jj) + zex - zxlim1)*SIGN(1.,plamt0(ji,jj) - zex - zxlim1) < 0. )        )  THEN
                   !
                   pk_top0(ji,jj) = 0    ! = land
                   !
                END IF
              END DO
            END DO
           ! pphit(:,:) = pphit0(:,:) ; plamt(:,:) = plamt0(:,:)
      ENDIF
      !
      !                       !== Horizontal scale factors ==! (in meters)
      !
      !                                         ! constant grid spacing
      pe1t(:,:) =  ze1
      pe1u(:,:) =  ze1
      pe1v(:,:) =  ze1
      pe1f(:,:) =  ze1

      pe2t(:,:) = ze1
      pe2u(:,:) = ze1
      pe2v(:,:) = ze1
      pe2f(:,:) = ze1

      IF (ln_bnd_refine) THEN
        WRITE(numout,*) 'North boundary locally refine'
        !
      ze2 = ze1 / rn_bnd_ref ; ztm = rpi / REAL(nn_bnd_ntr,wp)
        zelrefine = 2000000_wp - rn_bnd_len*1000._wp - rn_bnd_eqi * ze1
        WRITE(numout,*) 'zelrefine = ', zelrefine
        jflag = 0
        !
        ! 1/ gphiv
        DO jj = 1, jpj
            IF ( jflag == 0 ) THEN
              !
              zjm1 = REAL( jj + njmpp - 1 ) * ze1 - zphi0
              !
              IF ( zjm1 >= zelrefine ) jflag = 1 ; ztmp = 0._wp
            ELSE
              IF (ztmp <= nn_bnd_ntr) THEN
                WRITE(numout,*) 'COS( ztm *',ztmp,' )=',0.5_wp * ( (ze1 - ze2)*COS( ztm * ztmp ) + ze1 + ze2)
                zjm1 = zjm1 + 0.5_wp * ( (ze1 - ze2)*COS( ztm * ztmp ) + ze1 + ze2)
                ztmp = ztmp + 1._wp
              ELSE
                zjm1 = zjm1 + ze2
              ENDIF
            ENDIF
            !
            DO ji = 1, jpi
              zim05 = REAL( ji + nimpp - 1 ) - 0.5
              pphiv(ji,jj) = +zim05 * ze1 * zsin_theta + zjm1
            END DO
          END DO
          !
          !
          ! ze1 est déjà 100km/nn_AM98
          ! 1/ gphiv
          ! ratio = 2
          ! x=0.3, Nb_cells 3, EquivLength 2.0
          ! ['1.000', '1.300', '1.690']
          ! ratio = 4
          ! x=0.3, Nb_cells 6, EquivLength 3.2
          ! ['1.000', '1.300', '1.690', '2.197', '2.856', '3.713']
          !
          ! znI = 5._wp ; znq = 0.125_wp ! nn_smo = znI / rn_abp = znq
          !
          !znI = nn_refn ;  znq = rn_refr ! (actual Nb of needed cells and 1/ratio), rn_cnp=variation (x)
          ! znI = 6 ;  znq = 0.25_wp ! 1/16 (1/r)
          ! znI = 3 ;  znq = 0.5_wp ! 1/8 (1/r)
          !
          !jei0 = jpj - (nn_gc + 2 + nn_refn_equiv ) ! taille hypothètique (nn_smo = number of equiv cells)
          !
          ! DO jj = 1, jpj
          !    DO ji = 1, jpi
          !      !
          !      zim05 = REAL( ji + nimpp - 1 ) - 0.5
          !      !
          !      zel = 0._wp
          !      DO jk = 0,(jj + njmpp - 1)
          !        zek = REAL(jk - jei0, wp)   !
          !        IF (zek <  0._wp ) zfac = 1._wp                                      ! resolution nn_AM98
          !        IF (zek >= 0._wp ) zfac = MAX( znq*(1+rn_refx)**(znI-1 - zek) ,znq)   !
          !        ! IF (zek >= 0._wp ) zfac = znq   !
          !        zel = zel + zfac
          !      END DO
          !      zjm1 = zel
          !      !
          !      pphiv(ji,jj) = +zim05 * ze1 * zsin_theta + zjm1 * ze1 * zcos_theta - zphi0
          !   END DO
          ! END DO
          ! 2/ gphit
          DO jj = 2, jpj
             DO ji = 1, jpi
               pphit(ji,jj) = 0.5_wp * (pphiv(ji,jj) + pphiv(ji,jj-1))
            END DO
          END DO
          ! 3/ gphiu and f
          pphiu(:,:) = pphit(:,:) ; pphif(:,:) = pphiv(:,:)
          ! 3/ e2
          DO jj = 2, jpj
             DO ji = 1, jpi
               pe2t(ji,jj  ) = pphiv(ji,jj) - pphiv(ji,jj-1)
               pe2u(ji,jj  ) = pphiv(ji,jj) - pphiv(ji,jj-1)
               !
               pe2v(ji,jj-1) = pphit(ji,jj) - pphit(ji,jj-1)
               pe2f(ji,jj-1) = pphit(ji,jj) - pphit(ji,jj-1)
            END DO
          END DO
      ENDIF
      !
      !                                         ! NO reduction of grid size in some straits
      ke1e2u_v = 0                              !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp                       !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp                       !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                                            !  indicate not to compute ff afterward
      !
      pff_f(:,:) =  REAL( rn_f0, wp ) + REAL( rn_beta, wp ) * pphif(:,:)  ! f = f0 +beta* y
      pff_t(:,:) =  REAL( rn_f0, wp ) + REAL( rn_beta, wp ) * pphit(:,:)  ! f = f0 +beta* y
      !
      IF(lwp) WRITE(numout,*) '                           beta-plane used. f0   = ', rn_f0 ,  ' 1/s'
      IF(lwp) WRITE(numout,*) '                           beta-plane used. beta = ', rn_beta, ' 1/(s.m)'
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr

! #else
!       pe1t(:,:) =  ze1     ;      pe2t(:,:) = ze1
!       pe1u(:,:) =  ze1     ;      pe2u(:,:) = ze1
!       pe1v(:,:) =  ze1     ;      pe2v(:,:) = ze1
!       pe1f(:,:) =  ze1     ;      pe2f(:,:) = ze1
! #endif
