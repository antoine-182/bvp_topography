MODULE dynvor
   !!======================================================================
   !!                       ***  MODULE  dynvor  ***
   !! Ocean dynamics: Update the momentum trend with the relative and
   !!                 planetary vorticity trends
   !!======================================================================
   !! History :  OPA  ! 1989-12  (P. Andrich)  vor_ens: Original code
   !!            5.0  ! 1991-11  (G. Madec)  vor_ene, vor_mix: Original code
   !!            6.0  ! 1996-01  (G. Madec)  s-coord, suppress work arrays
   !!   NEMO     0.5  ! 2002-08  (G. Madec)  F90: Free form and module
   !!            1.0  ! 2004-02  (G. Madec)  vor_een: Original code
   !!             -   ! 2003-08  (G. Madec)  add vor_ctl
   !!             -   ! 2005-11  (G. Madec)  add dyn_vor (new step architecture)
   !!            2.0  ! 2006-11  (G. Madec)  flux form advection: add metric term
   !!            3.2  ! 2009-04  (R. Benshila)  vvl: correction of een scheme
   !!            3.3  ! 2010-10  (C. Ethe, G. Madec)  reorganisation of initialisation phase
   !!            3.7  ! 2014-04  (G. Madec)  trend simplification: suppress jpdyn_trd_dat vorticity
   !!             -   ! 2014-06  (G. Madec)  suppression of velocity curl from in-core memory
   !!             -   ! 2016-12  (G. Madec, E. Clementi) add Stokes-Coriolis trends (ln_stcor=T)
   !!            4.0  ! 2017-07  (G. Madec)  linear dynamics + trends diag. with Stokes-Coriolis
   !!             -   ! 2018-03  (G. Madec)  add two new schemes (ln_dynvor_enT and ln_dynvor_eet)
   !!             -   ! 2018-04  (G. Madec)  add pre-computed gradient for metric term calculation
   !!            4.x  ! 2020-03  (G. Madec, A. Nasser)  make ln_dynvor_msk truly efficient on relative vorticity
   !!            4.x  ! 2020-03  (G. Madec, A. Nasser)  alternate direction computation of vorticity tendancy
   !!                 !                                 for ENS, ENE
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_vor       : Update the momentum trend with the vorticity trend
   !!       vor_ens   : enstrophy conserving scheme       (ln_dynvor_ens=T)
   !!       vor_ene   : energy conserving scheme          (ln_dynvor_ene=T)
   !!       vor_een   : energy and enstrophy conserving   (ln_dynvor_een=T)
   !!   dyn_vor_init  : set and control of the different vorticity option
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE dommsk         ! ocean mask
   USE dynadv         ! momentum advection
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   USE sbcwave        ! Surface Waves (add Stokes-Coriolis force)
   USE sbc_oce , ONLY : ln_stcor    ! use Stoke-Coriolis force
   !
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
!!an
   USE dynkeg, ONLY : dyn_kegAD
   USE dynadv, ONLY : nn_dynkeg
   USE usrdef_nam , ONLY : ln_mir_coriolis !, ln_mir_e3f
!!an
   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_vor        ! routine called by step.F90
   PUBLIC   dyn_vor_init   ! routine called by nemogcm.F90

   !                                   !!* Namelist namdyn_vor: vorticity term
   LOGICAL, PUBLIC ::   ln_dynvor_ens      !: enstrophy conserving scheme          (ENS)
   LOGICAL, PUBLIC ::   ln_dynvor_ens_adVO   = .FALSE.   !: AD enstrophy conserving scheme       (ENS_ad)
   LOGICAL, PUBLIC ::   ln_dynvor_ens_adKE   = .FALSE.   !: AD enstrophy conserving scheme       (ENS_ad)
   LOGICAL, PUBLIC ::   ln_dynvor_ens_adKEVO = .FALSE.   !: AD enstrophy conserving scheme       (ENS_ad)
   LOGICAL, PUBLIC ::   ln_dynvor_ene      !: f-point energy conserving scheme     (ENE)
   LOGICAL, PUBLIC ::   ln_dynvor_ene_adVO   = .FALSE.   !: f-point AD energy conserving scheme  (ENE_ad)
   LOGICAL, PUBLIC ::   ln_dynvor_ene_adKE   = .FALSE.   !: f-point AD energy conserving scheme  (ENE_ad)
   LOGICAL, PUBLIC ::   ln_dynvor_ene_adKEVO = .FALSE.   !: f-point AD energy conserving scheme  (ENE_ad)
   LOGICAL, PUBLIC ::   ln_dynvor_enT      !: t-point energy conserving scheme     (ENT)
   LOGICAL, PUBLIC ::   ln_dynvor_eeT      !: t-point energy conserving scheme     (EET)
   LOGICAL, PUBLIC ::   ln_dynvor_een      !: energy & enstrophy conserving scheme (EEN)
   INTEGER, PUBLIC ::      nn_een_e3f         !: e3f=masked averaging of e3t divided by 4 (=0) or by the sum of mask (=1)
   LOGICAL, PUBLIC ::   ln_dynvor_mix      !: mixed scheme                         (MIX)
!!an
   LOGICAL, PUBLIC ::   ln_dynvor_mix_een_ens      !: mixed scheme                         (MIX)
   LOGICAL, PUBLIC ::   ln_dynvor_mix_ens_een      !: mixed scheme                         (MIX)
!!an
   LOGICAL, PUBLIC ::   ln_dynvor_msk      !: vorticity multiplied by fmask (=T) or not (=F) (all vorticity schemes)

   INTEGER, PUBLIC ::   nvor_scheme     !: choice of the type of advection scheme
   !                                       ! associated indices:
   INTEGER, PUBLIC, PARAMETER ::   np_ENS = 0   ! ENS scheme
   INTEGER, PUBLIC, PARAMETER ::   np_ENE = 1   ! ENE scheme
   INTEGER, PUBLIC, PARAMETER ::   np_ENT = 2   ! ENT scheme (t-point vorticity)
   INTEGER, PUBLIC, PARAMETER ::   np_EET = 3   ! EET scheme (EEN using e3t)
   INTEGER, PUBLIC, PARAMETER ::   np_EEN = 4   ! EEN scheme
   INTEGER, PUBLIC, PARAMETER ::   np_MIX = 5   ! MIX scheme
!!an
   INTEGER, PUBLIC, PARAMETER ::   np_MIX1 = 6   ! MIX scheme f-ENS dyn-EEN
   INTEGER, PUBLIC, PARAMETER ::   np_MIX2 = 7   ! MIX scheme f-EEN dyn-ENS
!!an
   INTEGER, PUBLIC, PARAMETER ::   np_ENS_adKE   = 11   ! ENS scheme - AD scheme (KE only)
   INTEGER, PUBLIC, PARAMETER ::   np_ENS_adVO   = 12   ! ENS scheme - AD scheme (VOR only)
   INTEGER, PUBLIC, PARAMETER ::   np_ENS_adKEVO = 13   ! ENS scheme - AD scheme (KE+VOR)
   INTEGER, PUBLIC, PARAMETER ::   np_ENE_adKE   = 21   ! ENE scheme - AD scheme (KE only)
   INTEGER, PUBLIC, PARAMETER ::   np_ENE_adVO   = 22   ! ENE scheme - AD scheme (VOR only)
   INTEGER, PUBLIC, PARAMETER ::   np_ENE_adKEVO = 23   ! ENE scheme - AD scheme (KE+VOR)
!!an

!!an ds step on pourra spécifier la valeur de ntot = np_COR ou np_COR + np_RVO
   INTEGER, PUBLIC ::   ncor, nrvm, ntot   ! choice of calculated vorticity
   !                               ! associated indices:
   INTEGER, PUBLIC, PARAMETER ::   np_COR = 1         ! Coriolis (planetary)
   INTEGER, PUBLIC, PARAMETER ::   np_RVO = 2         ! relative vorticity
   INTEGER, PUBLIC, PARAMETER ::   np_MET = 3         ! metric term
   INTEGER, PUBLIC, PARAMETER ::   np_CRV = 4         ! relative + planetary (total vorticity)
   INTEGER, PUBLIC, PARAMETER ::   np_CME = 5         ! Coriolis + metric term

   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   di_e2u_2        ! = di(e2u)/2          used in T-point metric term calculation
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   dj_e1v_2        ! = dj(e1v)/2           -        -      -       -
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   di_e2v_2e1e2f   ! = di(e2u)/(2*e1e2f)  used in F-point metric term calculation
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   dj_e1u_2e1e2f   ! = dj(e1v)/(2*e1e2f)   -        -      -       -

   REAL(wp) ::   r1_4  = 0.250_wp         ! =1/4
   REAL(wp) ::   r1_8  = 0.125_wp         ! =1/8
   REAL(wp) ::   r1_12 = 1._wp / 12._wp   ! 1/12

   !! * Substitutions
#  include "do_loop_substitute.h90"
!!st23
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynvor.F90 13154 2020-06-24 13:31:32Z gm $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_vor( kt, Kbb, Kmm, puu, pvv, Krhs)
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!
      !! ** Action : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the now vorticity term trend
      !!             - save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) and send them to trd_dyn
      !!               for futher diagnostics (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk   ! dummy loop indice
      INTEGER                             , INTENT( in  ) ::   kt               ! ocean time-step index
      INTEGER                             , INTENT( in  ) ::   Kmm, Krhs, Kbb   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::   puu, pvv         ! ocean velocity field and RHS of momentum equation
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  zuu, zvv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_vor')
      !
      IF( l_trddyn ) THEN     !==  trend diagnostics case : split the added trend in two parts  ==!
         !
         ALLOCATE( ztrdu(jpi,jpj,jpk), ztrdv(jpi,jpj,jpk) )
         !
         ztrdu(:,:,:) = puu(:,:,:,Krhs)            !* planetary vorticity trend (including Stokes-Coriolis force)
         ztrdv(:,:,:) = pvv(:,:,:,Krhs)
         SELECT CASE( nvor_scheme )
         CASE( np_ENS )           ;   CALL vor_ens( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! enstrophy conserving scheme
            IF( ln_stcor )            CALL vor_ens( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
         CASE( np_ENE, np_MIX )   ;   CALL vor_ene( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! energy conserving scheme
            IF( ln_stcor )            CALL vor_ene( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
         CASE( np_ENT )           ;   CALL vor_enT( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! energy conserving scheme (T-pts)
            IF( ln_stcor )            CALL vor_enT( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
         CASE( np_EET )           ;   CALL vor_eeT( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! energy conserving scheme (een with e3t)
            IF( ln_stcor )            CALL vor_eeT( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
         CASE( np_EEN )           ;   CALL vor_een( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! energy & enstrophy scheme
            IF( ln_stcor )            CALL vor_een( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
         END SELECT
         ztrdu(:,:,:) = puu(:,:,:,Krhs) - ztrdu(:,:,:)
         ztrdv(:,:,:) = pvv(:,:,:,Krhs) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_pvo, kt, Kmm )
         !
         IF( n_dynadv /= np_LIN_dyn ) THEN   !* relative vorticity or metric trend (only in non-linear case)
            ztrdu(:,:,:) = puu(:,:,:,Krhs)
            ztrdv(:,:,:) = pvv(:,:,:,Krhs)
            SELECT CASE( nvor_scheme )
            CASE( np_ENT )           ;   CALL vor_enT( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! energy conserving scheme (T-pts)
            CASE( np_EET )           ;   CALL vor_eeT( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! energy conserving scheme (een with e3t)
            CASE( np_ENE )           ;   CALL vor_ene( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! energy conserving scheme
            CASE( np_ENS, np_MIX )   ;   CALL vor_ens( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! enstrophy conserving scheme
            CASE( np_EEN )           ;   CALL vor_een( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! energy & enstrophy scheme
            END SELECT
            ztrdu(:,:,:) = puu(:,:,:,Krhs) - ztrdu(:,:,:)
            ztrdv(:,:,:) = pvv(:,:,:,Krhs) - ztrdv(:,:,:)
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_rvo, kt, Kmm )
         ENDIF
         !
         DEALLOCATE( ztrdu, ztrdv )
         !
      ELSE              !==  total vorticity trend added to the general trend  ==!
         !
         SELECT CASE ( nvor_scheme )      !==  vorticity trend added to the general trend  ==!
         CASE( np_ENT )                        !* energy conserving scheme  (T-pts)
                             CALL vor_enT( kt, Kmm, ntot, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! total vorticity trend
            IF( ln_stcor )   CALL vor_enT( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
         CASE( np_EET )                        !* energy conserving scheme (een scheme using e3t)
                             CALL vor_eeT( kt, Kmm, ntot, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! total vorticity trend
            IF( ln_stcor )   CALL vor_eeT( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
         CASE( np_ENE )                        !* energy conserving scheme
                             CALL vor_ene( kt, Kmm, ntot, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! total vorticity trend
            IF( ln_stcor )   CALL vor_ene( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
         CASE( np_ENS )                        !* enstrophy conserving scheme
                             CALL vor_ens( kt, Kmm, ntot, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! total vorticity trend
            IF( ln_stcor )   CALL vor_ens( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! add the Stokes-Coriolis trend
         CASE( np_MIX )                        !* mixed ene-ens scheme
                             CALL vor_ens( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! relative vorticity or metric trend (ens)
                             CALL vor_ene( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! planetary vorticity trend (ene)
            IF( ln_stcor )   CALL vor_ene( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
          CASE( np_MIX1)                        !* mixed f-ens dyn-een scheme
                              CALL vor_een( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )
                              CALL vor_ens( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )
          CASE( np_MIX2)                        !* mixed f-een dyn-ens scheme
                              CALL vor_ens( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )
                              CALL vor_een( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )
         CASE( np_EEN )                        !* energy and enstrophy conserving scheme
                             CALL vor_een( kt, Kmm, ntot, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! total vorticity trend
            IF( ln_stcor )   CALL vor_een( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
         END SELECT
         !
      ENDIF
      !
      !                       ! print sum trends (used for debugging)
      IF(sn_cfctl%l_prtctl) CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' vor  - Ua: ', mask1=umask,               &
         &                                tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_vor')
      !
   END SUBROUTINE dyn_vor


   SUBROUTINE vor_enT( kt, Kmm, kvor, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE vor_enT  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!       and t-point evaluation of vorticity (planetary and relative).
      !!       conserves the horizontal kinetic energy.
      !!         The general trend of momentum is increased due to the vorticity
      !!       term which is given by:
      !!          voru = 1/bu  mj[ ( mi(mj(bf*rvor))+bt*f_t)/e3t  mj[vn] ]
      !!          vorv = 1/bv  mi[ ( mi(mj(bf*rvor))+bt*f_t)/e3f  mj[un] ]
      !!       where rvor is the relative vorticity at f-point
      !!
      !! ** Action : - Update (pu_rhs,pv_rhs) with the now vorticity term trend
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt               ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kmm              ! ocean time level index
      INTEGER                         , INTENT(in   ) ::   kvor             ! total, planetary, relative, or metric
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu, pv           ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs   ! total v-trend
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zx1, zy1, zx2, zy2   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj)     ::   zwx, zwy, zwt   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwz             ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj)     ::   zwtf, zwtz   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj)     ::   zr1_e1e2f        ! 2D workspace penalised
      REAL(wp), DIMENSION(jpi,jpj,jpk,3) ::   ze3t, ze3u, ze3v             ! 4D workspace penalised
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_enT : vorticity term: t-point energy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      !
      zr1_e1e2f(:,:) = r1_e1e2f(:,:)
      ze3t(:,:,:,:) = e3t(:,:,:,:)
      ze3u(:,:,:,:) = e3u(:,:,:,:)
      ze3v(:,:,:,:) = e3v(:,:,:,:)
!
      SELECT CASE( kvor )                 !== relative vorticity considered  ==!
      CASE ( np_RVO , np_CRV )                  !* relative vorticity at f-point is used
         DO jk = 1, jpkm1                                ! Horizontal slab
            DO_2D_10_10
               zwz(ji,jj,jk) = (  e2v(ji+1,jj) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk)  &
                  &             - e1u(ji,jj+1) * pu(ji,jj+1,jk) + e1u(ji,jj) * pu(ji,jj,jk)  ) * zr1_e1e2f(ji,jj)
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask relative vorticity
               DO_2D_10_10
                  zwz(ji,jj,jk) = zwz(ji,jj,jk) * fmask(ji,jj,jk)
               END_2D
            ENDIF
         END DO
         CALL lbc_lnk( 'dynvor', zwz, 'F', 1. )
         !
      END SELECT
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         SELECT CASE( kvor )                 !==  volume weighted vorticity considered  ==!
         !
         CASE ( np_COR )                           !* Coriolis (planetary vorticity)
            zwt(:,:) = ff_t(:,:) * e1e2t(:,:) * ze3t(:,:,jk,Kmm)
         CASE ( np_RVO )                           !* relative vorticity
            DO_2D_01_01
               zwt(ji,jj) = r1_4 * (   zwz(ji-1,jj  ,jk) + zwz(ji,jj  ,jk)   &
                  &                  + zwz(ji-1,jj-1,jk) + zwz(ji,jj-1,jk) ) &
                  &                * e1e2t(ji,jj)*e3t(ji,jj,jk,Kmm)
            END_2D
         CASE ( np_MET )                           !* metric term
            DO_2D_01_01
               zwt(ji,jj) = (   ( pv(ji,jj,jk) + pv(ji,jj-1,jk) ) * di_e2u_2(ji,jj)     &
                  &           - ( pu(ji,jj,jk) + pu(ji-1,jj,jk) ) * dj_e1v_2(ji,jj)   ) &
                  &         * e3t(ji,jj,jk,Kmm)
            END_2D
         CASE ( np_CRV )                           !* Coriolis + relative vorticity
            DO_2D_01_01
               zwt(ji,jj) = (  ff_t(ji,jj) + r1_4 * ( zwz(ji-1,jj  ,jk) + zwz(ji,jj  ,jk)      &
                  &                                 + zwz(ji-1,jj-1,jk) + zwz(ji,jj-1,jk) )  ) &
                  &         * e1e2t(ji,jj)*ze3t(ji,jj,jk,Kmm)
            END_2D
         CASE ( np_CME )                           !* Coriolis + metric
            DO_2D_01_01
               zwt(ji,jj) = (  ff_t(ji,jj) * e1e2t(ji,jj)                             &
                    &        + ( pv(ji,jj,jk) + pv(ji,jj-1,jk) ) * di_e2u_2(ji,jj)    &
                    &        - ( pu(ji,jj,jk) + pu(ji-1,jj,jk) ) * dj_e1v_2(ji,jj)  ) &
                    &       * ze3t(ji,jj,jk,Kmm)
            END_2D
         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','dyn_vor: wrong value for kvor'  )
         END SELECT
         DO_2D_00_00
          pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + r1_4 * r1_e1e2u(ji,jj) / ze3u(ji,jj,jk,Kmm)                &
             &                                * (  zwt(ji+1,jj) * ( pv(ji+1,jj,jk) + pv(ji+1,jj-1,jk) )   &
             &                                   + zwt(ji  ,jj) * ( pv(ji  ,jj,jk) + pv(ji  ,jj-1,jk) )   )
             !
          pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) - r1_4 * r1_e1e2v(ji,jj) / ze3v(ji,jj,jk,Kmm)                &
             &                                * (  zwt(ji,jj+1) * ( pu(ji,jj+1,jk) + pu(ji-1,jj+1,jk) )   &
             &                                   + zwt(ji,jj  ) * ( pu(ji,jj  ,jk) + pu(ji-1,jj  ,jk) )   )
         END_2D
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
   END SUBROUTINE vor_enT

   SUBROUTINE vor_ene( kt, Kmm, kvor, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE vor_ene  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!       and the Sadourny (1975) flux form formulation : conserves the
      !!       horizontal kinetic energy.
      !!         The general trend of momentum is increased due to the vorticity
      !!       term which is given by:
      !!          voru = 1/e1u  mj-1[ (rvor+f)/e3f  mi(e1v*e3v pvv(:,:,:,Kmm)) ]
      !!          vorv = 1/e2v  mi-1[ (rvor+f)/e3f  mj(e2u*e3u puu(:,:,:,Kmm)) ]
      !!       where rvor is the relative vorticity
      !!
      !! ** Action : - Update (pu_rhs,pv_rhs) with the now vorticity term trend
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      INTEGER                                  , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                                  , INTENT(in   ) ::   Kmm              ! ocean time level index
      INTEGER                                  , INTENT(in   ) ::   kvor        ! total, planetary, relative, or metric
      REAL(wp), DIMENSION(jpi,jpj,jpk)         , INTENT(inout) ::   pu, pv    ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk),OPTIONAL, INTENT(inout) ::   pu_rhs, pv_rhs    ! total v-trend
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zx1, zy1, zx2, zy2,z1d,zmsk   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zwx, zwy, zwz   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj) ::   zr1_e1e2f            ! 2D workspace penalised
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ze3f             ! 3D workspace penalised
      REAL(wp), DIMENSION(jpi,jpj,jpk,3) ::   ze3u, ze3v             ! 4D workspace penalised
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_ene : vorticity term: energy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
!!an in ENE (or ENS), zeta is divided by the penalised volume bf (curl(u)/hf tilde)
!!an and multiplied by the penalised flux. Coriolis is untouched
      zr1_e1e2f(:,:) = r1_e1e2f(:,:)
      ze3f(:,:,:) = e3f(:,:,:)
      ze3u(:,:,:,:) = e3u(:,:,:,:)
      ze3v(:,:,:,:) = e3v(:,:,:,:)
    !
      SELECT CASE( nn_een_e3f )           ! == reciprocal of e3 at F-point
        CASE ( 0 )                                   ! original formulation  (masked averaging of e3t divided by 4)
           ! rien
           ! à comparer
           ! DO_3D_10_10(1,jpkm1)
           !    z1d = (   e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)   &
           !       &    + e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk)   &
           !       &    + e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)   &
           !       &    + e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk)  )
           !    IF( z1d /= 0._wp ) THEN   ;   ze3f(ji,jj,jk) = z1d * 0.25_wp
           !    ELSE                      ;   ze3f(ji,jj,jk) = 1._wp
           !    ENDIF
           ! END_3D
        CASE ( 1 )                                   ! new formulation  (masked averaging of e3t divided by the sum of mask)
           DO_3D_10_10(1,jpkm1)
              z1d = (   e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)   &
                 &    + e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk)   &
                 &    + e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)   &
                 &    + e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk)  )
              zmsk = (                    tmask(ji,jj+1,jk) +                     tmask(ji+1,jj+1,jk)   &
                 &                      + tmask(ji,jj  ,jk) +                     tmask(ji+1,jj  ,jk)  )
              IF( z1d /= 0._wp ) THEN   ;   ze3f(ji,jj,jk) = z1d / zmsk
              ELSE                      ;   ze3f(ji,jj,jk) = 1._wp
              ENDIF
           END_3D
      END SELECT
    !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         SELECT CASE( kvor )                 !==  vorticity considered  ==!
         CASE ( np_COR )                           !* Coriolis (planetary vorticity)
            zwz(:,:) = ff_f(:,:)
         CASE ( np_RVO )                           !* relative vorticity
            DO_2D_10_10
               zwz(ji,jj) = (  e2v(ji+1,jj  ) * pv(ji+1,jj  ,jk) - e2v(ji,jj) * pv(ji,jj,jk)    &
                  &          - e1u(ji  ,jj+1) * pu(ji  ,jj+1,jk) + e1u(ji,jj) * pu(ji,jj,jk)  ) * zr1_e1e2f(ji,jj)
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask the relative vorticity
               DO_2D_10_10
                  zwz(ji,jj) = zwz(ji,jj) * fmask(ji,jj,jk)
               END_2D
            ENDIF
         CASE ( np_MET )                           !* metric term
            DO_2D_10_10
               zwz(ji,jj) = ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &       - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)
            END_2D
         CASE ( np_CRV )                           !* Coriolis + relative vorticity
            DO_2D_10_10
               zwz(ji,jj) = ff_f(ji,jj) + (  e2v(ji+1,jj) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk)      &
                  &                        - e1u(ji,jj+1) * pu(ji,jj+1,jk) + e1u(ji,jj) * pu(ji,jj,jk)  ) * zr1_e1e2f(ji,jj)
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask the relative vorticity (NOT the Coriolis term)
               DO_2D_10_10
                  zwz(ji,jj) = ( zwz(ji,jj) - ff_f(ji,jj) ) * fmask(ji,jj,jk) + ff_f(ji,jj)
               END_2D
            ENDIF
         CASE ( np_CME )                           !* Coriolis + metric
            DO_2D_10_10
               zwz(ji,jj) = ff_f(ji,jj) + ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &                     - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)
            END_2D
         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','dyn_vor: wrong value for kvor'  )
         END SELECT
         !
         IF( PRESENT( pu_rhs ) .AND. PRESENT( pv_rhs ) ) THEN     !***  NO alternating direction  ***!
            !
            !                                   !==  horizontal fluxes and potential vorticity ==!
            zwx(:,:) = e2u(:,:) * ze3u(:,:,jk,Kmm) * pu(:,:,jk)
            zwy(:,:) = e1v(:,:) * ze3v(:,:,jk,Kmm) * pv(:,:,jk)
            zwz(:,:) = zwz(:,:) / ze3f(:,:,jk)
            !
            !                                   !==  compute and add the vorticity term trend  =!
            DO_2D_00_00
               zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
               zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
               zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
               zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
               pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + r1_4 * r1_e1u(ji,jj) * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
               pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) - r1_4 * r1_e2v(ji,jj) * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
            END_2D
            !
            !
         ELSEIF(       PRESENT( pu_rhs ) .AND. .NOT. PRESENT( pv_rhs ) ) THEN            !***  Alternating direction : i-component  ***!
            !
            !
            !                                   !==  horizontal fluxes and potential vorticity ==!
!!anAD remplacer pu par uu(Nnn) et on fait correctement la direction alternée
!      on transportera (en AD) de la bonne manière
            !zwy(:,:) = e1v(:,:) * e3v(:,:,jk,Kmm) * pv(:,:,jk)
            zwy(:,:) = e1v(:,:) * e3v(:,:,jk,Kmm) * vv(:,:,jk,Kmm)
            zwz(:,:) = zwz(:,:) / e3f(:,:,jk)
            !
            !                                   !==  compute and add the vorticity term trend  =!
            DO_2D_00_00
               zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
               zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
               pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + r1_4 * r1_e1u(ji,jj) * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
            END_2D
            !
         ELSEIF( .NOT. PRESENT( pu_rhs ) .AND.       PRESENT( pv_rhs ) ) THEN            !***  Alternating direction : j-component  ***!
            !
            !
            !                                   !==  horizontal fluxes and potential vorticity ==!
!!anAD remplacer pu par uu(Nnn) et on fait correctement la direction alternée
            !zwx(:,:) = e2u(:,:) * e3u(:,:,jk,Kmm) * pu(:,:,jk)
            zwx(:,:) = e2u(:,:) * e3u(:,:,jk,Kmm) * uu(:,:,jk,Kmm)
            zwz(:,:) = zwz(:,:) / e3f(:,:,jk)
            !
            !                                   !==  compute and add the vorticity term trend  =!
            DO_2D_00_00
               zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
               zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
               pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) - r1_4 * r1_e2v(ji,jj) * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
            END_2D
            !
         ENDIF
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
   END SUBROUTINE vor_ene


   SUBROUTINE vor_ens( kt, Kmm, kvor, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_ens  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      and the Sadourny (1975) flux FORM formulation : conserves the
      !!      potential enstrophy of a horizontally non-divergent flow. the
      !!      trend of the vorticity term is given by:
      !!          voru = 1/e1u  mj-1[ (rvor+f)/e3f ]  mj-1[ mi(e1v*e3v pvv(:,:,:,Kmm)) ]
      !!          vorv = 1/e2v  mi-1[ (rvor+f)/e3f ]  mi-1[ mj(e2u*e3u puu(:,:,:,Kmm)) ]
      !!      Add this trend to the general momentum trend:
      !!          (u(rhs),v(Krhs)) = (u(rhs),v(Krhs)) + ( voru , vorv )
      !!
      !! ** Action : - Update (pu_rhs,pv_rhs)) arrays with the now vorticity term trend
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      INTEGER                                  , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                                  , INTENT(in   ) ::   Kmm              ! ocean time level index
      INTEGER                                  , INTENT(in   ) ::   kvor        ! total, planetary, relative, or metric
      REAL(wp), DIMENSION(jpi,jpj,jpk)         , INTENT(inout) ::   pu, pv    ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk),OPTIONAL, INTENT(inout) ::   pu_rhs, pv_rhs    ! total v-trend
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zuav, zvau, z1d, zmsk   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj)     ::   zwx, zwy, zwz, zww   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj)     ::   zr1_e1e2f            ! 2D workspace penalised
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ze3f, z3d, zuu, zvv                 ! 3D workspace penalised
      REAL(wp), DIMENSION(jpi,jpj,jpk,3) ::   ze3u, ze3v         ! 4D workspace penalised
    !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_ens : vorticity term: enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      zr1_e1e2f(:,:) = r1_e1e2f(:,:)
      ze3f(:,:,:) = e3f(:,:,:)
      ze3u(:,:,:,:) = e3u(:,:,:,:)
      ze3v(:,:,:,:) = e3v(:,:,:,:)
      z3d(:,:,:) = 1._wp
      !
      !
      !
      zuu(:,:,:) = pu(:,:,:) ; zvv = pv(:,:,:)
      !
      ! IF ( ln_mir_coriolis ) THEN
      !   DO jj = 1,jpjm1
      !     DO ji = 1,jpim1
      !       !
      !       ! wetmean
      !     END DO
      !   END DO
      ! ENDIF
      !
      IF ( ln_mir_coriolis ) THEN
        DO jj = 1,jpjm1
          DO ji = 1,jpim1
              IF ( (     tmask(ji  ,jj+1,1) + tmask(ji+1,jj+1,1)  &
                   &   + tmask(ji  ,jj  ,1) + tmask(ji+1,jj  ,1)  ) == 3._wp) THEN ! salient angle
                  !
                  IF      ( tmask(ji  ,jj  ,1)==0._wp ) THEN   ! South
                   zuu(ji,jj,1) = - pv(ji+1,jj  ,1)
                   zvv(ji,jj,1) = - pu(ji  ,jj+1,1)
                  ELSE IF ( tmask(ji+1,jj  ,1)==0._wp ) THEN   ! East
                    zuu(ji  ,jj,1) = + pv(ji,jj  ,1)
                    zvv(ji+1,jj,1) = + pu(ji,jj+1,1)
                  ELSE IF ( tmask(ji+1,jj+1,1)==0._wp ) THEN   ! North
                    zuu(ji  ,jj+1,1) = - pv(ji,jj,1)
                    zvv(ji+1,jj  ,1) = - pu(ji,jj,1)
                  ELSE IF ( tmask(ji  ,jj+1,1)==0._wp ) THEN   ! West
                    zuu(ji,jj+1,1) = + pv(ji+1,jj,1)
                    zvv(ji,jj  ,1) = + pu(ji  ,jj,1)
                  END IF
                  !
                END IF
              !
            END DO
          END DO
        ENDIF
      !
      SELECT CASE( nn_een_e3f )           ! == reciprocal of e3 at F-point
        CASE ( 0 )
          ! rien                          ! original formulation  (masked averaging of e3t divided by 4)
          !
        CASE ( 1 )                                   ! new formulation  (masked averaging of e3t divided by the sum of mask)
           DO_3D_10_10(1,jpkm1)
              z1d = (   e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)   &
                 &    + e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk)   &
                 &    + e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)   &
                 &    + e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk)  )
              zmsk = (                    tmask(ji,jj+1,jk) +                     tmask(ji+1,jj+1,jk)   &
                 &                      + tmask(ji,jj  ,jk) +                     tmask(ji+1,jj  ,jk)  )
              IF( z1d /= 0._wp ) THEN   ;   ze3f(ji,jj,jk) = z1d / zmsk
              ELSE                      ;   ze3f(ji,jj,jk) = 1._wp
              ENDIF
           END_3D
      END SELECT
      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         SELECT CASE( kvor )                 !==  vorticity considered  ==!
         CASE ( np_COR )                           !* Coriolis (planetary vorticity)
            zwz(:,:) = ff_f(:,:)
         CASE ( np_RVO )                           !* relative vorticity
            DO_2D_10_10
               zwz(ji,jj) = (  e2v(ji+1,jj  ) * pv(ji+1,jj  ,jk) - e2v(ji,jj) * pv(ji,jj,jk)    &
                  &          - e1u(ji  ,jj+1) * pu(ji  ,jj+1,jk) + e1u(ji,jj) * pu(ji,jj,jk)  ) * zr1_e1e2f(ji,jj)
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask the relative vorticity
               DO_2D_10_10
                  zwz(ji,jj) = zwz(ji,jj) * fmask(ji,jj,jk)
               END_2D
            ENDIF
         CASE ( np_MET )                           !* metric term
            DO_2D_10_10
               zwz(ji,jj) = ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &       - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)
            END_2D
         CASE ( np_CRV )                           !* Coriolis + relative vorticity
            DO_2D_10_10
               zwz(ji,jj) = ff_f(ji,jj) + (  e2v(ji+1,jj  ) * pv(ji+1,jj  ,jk) - e2v(ji,jj) * pv(ji,jj,jk)  &
                  &                        - e1u(ji  ,jj+1) * pu(ji  ,jj+1,jk) + e1u(ji,jj) * pu(ji,jj,jk)  ) * zr1_e1e2f(ji,jj)
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask the relative vorticity (NOT the Coriolis term)
               DO_2D_10_10
                  zwz(ji,jj) = ( zwz(ji,jj) - ff_f(ji,jj) ) * fmask(ji,jj,jk) * z3d(ji,jj,jk) + ff_f(ji,jj)
               END_2D
            ENDIF
         CASE ( np_CME )                           !* Coriolis + metric
            DO_2D_10_10
               zwz(ji,jj) = ff_f(ji,jj) + ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &                     - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)
            END_2D
         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','dyn_vor: wrong value for kvor'  )
         END SELECT
         !
         !
!!an wut ? v et u
         IF( PRESENT( pu_rhs ) .AND. PRESENT( pv_rhs ) ) THEN     !***  NO alternating direction  ***!
        !
#if defined key_am98h
              ! vertical scale factor not taken into account in the ENS scheme (like AM98)
              !                                   !==  horizontal fluxes and potential vorticity ==!
              zwx(:,:) = e2u(:,:) * pu(:,:,jk)
              zwy(:,:) = e1v(:,:) * pv(:,:,jk)
              zwz(:,:) = zwz(:,:)
              !
#else

              !                                   !==  horizontal fluxes and potential vorticity ==!
              zwx(:,:) = e2u(:,:) * e3u(:,:,jk,Kmm) * zuu(:,:,jk)
              zwy(:,:) = e1v(:,:) * e3v(:,:,jk,Kmm) * zvv(:,:,jk)
              zwz(:,:) = zwz(:,:) / ze3f(:,:,jk)
#endif
              !
              !                                   !==  compute and add the vorticity term trend  =!
              DO_2D_00_00
                 zuav = r1_8 * r1_e1u(ji,jj) * (  zwy(ji  ,jj-1) + zwy(ji+1,jj-1)  &
                    &                           + zwy(ji  ,jj  ) + zwy(ji+1,jj  )  )
                 zvau =-r1_8 * r1_e2v(ji,jj) * (  zwx(ji-1,jj  ) + zwx(ji-1,jj+1)  &
                    &                           + zwx(ji  ,jj  ) + zwx(ji  ,jj+1)  )
                 pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zuav * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
                 pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zvau * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
              END_2D
              !

         ELSEIF(       PRESENT( pu_rhs ) .AND. .NOT. PRESENT( pv_rhs ) ) THEN            !***  Alternating direction : i-component  ***!
            !
            !
            !                                   !==  horizontal fluxes and potential vorticity ==!
            zwy(:,:) = e1v(:,:) * e3v(:,:,jk,Kmm) * pv(:,:,jk)
            zwz(:,:) = zwz(:,:) / e3f(:,:,jk)
            !
            !                                   !==  compute and add the vorticity term trend  =!
            DO_2D_00_00
               zuav = r1_8 * r1_e1u(ji,jj) * (  zwy(ji  ,jj-1) + zwy(ji+1,jj-1)  &
                  &                           + zwy(ji  ,jj  ) + zwy(ji+1,jj  )  )
               pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zuav * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
            END_2D
            !
         ELSEIF( .NOT. PRESENT( pu_rhs ) .AND.       PRESENT( pv_rhs ) ) THEN            !***  Alternating direction : j-component  ***!
            !
            !
            !                                   !==  horizontal fluxes and potential vorticity ==!
            zwx(:,:) = e2u(:,:) * e3u(:,:,jk,Kmm) * pu(:,:,jk)
            zwz(:,:) = zwz(:,:) / e3f(:,:,jk)
            !
            !                                   !==  compute and add the vorticity term trend  =!
            DO_2D_00_00
               zvau =-r1_8 * r1_e2v(ji,jj) * (  zwx(ji-1,jj  ) + zwx(ji-1,jj+1)  &
                  &                           + zwx(ji  ,jj  ) + zwx(ji  ,jj+1)  )
               pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zvau * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
            END_2D
            !
         ENDIF
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
   END SUBROUTINE vor_ens


   SUBROUTINE vor_een( kt, Kmm, kvor, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_een  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      and the Arakawa and Lamb (1980) flux form formulation : conserves
      !!      both the horizontal kinetic energy and the potential enstrophy
      !!      when horizontal divergence is zero (see the NEMO documentation)
      !!      Add this trend to the general momentum trend (pu_rhs,pv_rhs).
      !!
      !! ** Action : - Update (pu_rhs,pv_rhs) with the now vorticity term trend
      !!
      !! References : Arakawa and Lamb 1980, Mon. Wea. Rev., 109, 18-36
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kmm              ! ocean time level index
      INTEGER                         , INTENT(in   ) ::   kvor        ! total, planetary, relative, or metric
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu, pv    ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs    ! total v-trend
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ierr         ! local integer
      REAL(wp) ::   zua, zva     ! local scalars
      REAL(wp) ::   zmsk, ze3f   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj)     ::   zwx , zwy , z1_e3f
      REAL(wp), DIMENSION(jpi,jpj)     ::   ztnw, ztne, ztsw, ztse
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwz
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_een : vorticity term: energy and enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         SELECT CASE( nn_een_e3f )           ! == reciprocal of e3 at F-point
         CASE ( 0 )                                   ! original formulation  (masked averaging of e3t divided by 4)
            ! DO_2D_10_10
            !    ze3f = (  e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)   &
            !       &    + e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk)   &
            !       &    + e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)   &
            !       &    + e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk)  )
            !    IF( ze3f /= 0._wp ) THEN   ;   z1_e3f(ji,jj) = 4._wp / ze3f
            !    ELSE                       ;   z1_e3f(ji,jj) = 0._wp
            !    ENDIF
            ! END_2D
            z1_e3f(:,:) = 1._wp / e3f(:,:,jk)
         CASE ( 1 )                                   ! new formulation  (masked averaging of e3t divided by the sum of mask)
            DO_2D_10_10
               ze3f = (  e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)   &
                  &    + e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk)   &
                  &    + e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)   &
                  &    + e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk)  )
               zmsk = (                    tmask(ji,jj+1,jk) +                     tmask(ji+1,jj+1,jk)   &
                  &                      + tmask(ji,jj  ,jk) +                     tmask(ji+1,jj  ,jk)  )
               IF( ze3f /= 0._wp ) THEN   ;   z1_e3f(ji,jj) = zmsk / ze3f
               ELSE                       ;   z1_e3f(ji,jj) = 0._wp
               ENDIF
            END_2D
         END SELECT
         !
         SELECT CASE( kvor )                 !==  vorticity considered  ==!
         CASE ( np_COR )                           !* Coriolis (planetary vorticity)
            DO_2D_10_10
               zwz(ji,jj,jk) = ff_f(ji,jj) * z1_e3f(ji,jj)
            END_2D
         CASE ( np_RVO )                           !* relative vorticity
            DO_2D_10_10
               zwz(ji,jj,jk) = ( e2v(ji+1,jj  ) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk)  &
                  &            - e1u(ji  ,jj+1) * pu(ji,jj+1,jk) + e1u(ji,jj) * pu(ji,jj,jk)  ) * r1_e1e2f(ji,jj)*z1_e3f(ji,jj)
            END_2D
            !
            IF( ln_dynvor_msk ) THEN          !==  mask/unmask vorticity ==!
               DO_2D_10_10
                  zwz(ji,jj,jk) = zwz(ji,jj,jk) * fmask(ji,jj,jk)
               END_2D
            ENDIF
            !
         CASE ( np_MET )                           !* metric term
            DO_2D_10_10
               zwz(ji,jj,jk) = (   ( pv(ji+1,jj,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &              - ( pu(ji,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)   ) * z1_e3f(ji,jj)
            END_2D
         CASE ( np_CRV )                           !* Coriolis + relative vorticity
            DO_2D_10_10
               zwz(ji,jj,jk) =    e2v(ji+1,jj  ) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk)      &
                  &             - e1u(ji  ,jj+1) * pu(ji,jj+1,jk) + e1u(ji,jj) * pu(ji,jj,jk)
            END_2D
            !
            IF( ln_dynvor_msk ) THEN          !==  mask/unmask vorticity ==!
               DO_2D_10_10
                  zwz(ji,jj,jk) = zwz(ji,jj,jk) * fmask(ji,jj,jk)
               END_2D
            ENDIF
            !
            zwz(:,:,jk) = (  ff_f(:,:) + zwz(:,:,jk) * r1_e1e2f(:,:) ) * z1_e3f(:,:)
            !
         CASE ( np_CME )                           !* Coriolis + metric
            DO_2D_10_10
               zwz(ji,jj,jk) = (   ff_f(ji,jj) + ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &                            - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)   ) * z1_e3f(ji,jj)
            END_2D
         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','dyn_vor: wrong value for kvor'  )
         END SELECT
         !
      END DO                                           !   End of slab
         !
      CALL lbc_lnk( 'dynvor', zwz, 'F', 1. )

      DO jk = 1, jpkm1                                 ! Horizontal slab
         !
         !                                   !==  horizontal fluxes  ==!
         zwx(:,:) = e2u(:,:) * e3u(:,:,jk,Kmm) * pu(:,:,jk)
         zwy(:,:) = e1v(:,:) * e3v(:,:,jk,Kmm) * pv(:,:,jk)

         !                                   !==  compute and add the vorticity term trend  =!
         jj = 2
         ztne(1,:) = 0   ;   ztnw(1,:) = 0   ;   ztse(1,:) = 0   ;   ztsw(1,:) = 0
         DO ji = 2, jpi          ! split in 2 parts due to vector opt.
               ztne(ji,jj) = zwz(ji-1,jj  ,jk) + zwz(ji  ,jj  ,jk) + zwz(ji  ,jj-1,jk)
               ztnw(ji,jj) = zwz(ji-1,jj-1,jk) + zwz(ji-1,jj  ,jk) + zwz(ji  ,jj  ,jk)
               ztse(ji,jj) = zwz(ji  ,jj  ,jk) + zwz(ji  ,jj-1,jk) + zwz(ji-1,jj-1,jk)
               ztsw(ji,jj) = zwz(ji  ,jj-1,jk) + zwz(ji-1,jj-1,jk) + zwz(ji-1,jj  ,jk)
         END DO
         DO jj = 3, jpj
            DO ji = 2, jpi   ! vector opt. ok because we start at jj = 3
               ztne(ji,jj) = zwz(ji-1,jj  ,jk) + zwz(ji  ,jj  ,jk) + zwz(ji  ,jj-1,jk)
               ztnw(ji,jj) = zwz(ji-1,jj-1,jk) + zwz(ji-1,jj  ,jk) + zwz(ji  ,jj  ,jk)
               ztse(ji,jj) = zwz(ji  ,jj  ,jk) + zwz(ji  ,jj-1,jk) + zwz(ji-1,jj-1,jk)
               ztsw(ji,jj) = zwz(ji  ,jj-1,jk) + zwz(ji-1,jj-1,jk) + zwz(ji-1,jj  ,jk)
            END DO
         END DO
         DO_2D_00_00
            zua = + r1_12 * r1_e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
               &                             + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
            zva = - r1_12 * r1_e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
               &                             + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
            pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zua
            pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zva
         END_2D
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
   END SUBROUTINE vor_een

   SUBROUTINE vor_eeT( kt, Kmm, kvor, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_eeT  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      and the Arakawa and Lamb (1980) vector form formulation using
      !!      a modified version of Arakawa and Lamb (1980) scheme (see vor_een).
      !!      The change consists in
      !!      Add this trend to the general momentum trend (pu_rhs,pv_rhs).
      !!
      !! ** Action : - Update (pu_rhs,pv_rhs) with the now vorticity term trend
      !!
      !! References : Arakawa and Lamb 1980, Mon. Wea. Rev., 109, 18-36
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kmm              ! ocean time level index
      INTEGER                         , INTENT(in   ) ::   kvor        ! total, planetary, relative, or metric
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu, pv    ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs    ! total v-trend
      !
      INTEGER  ::   ji, jj, jk     ! dummy loop indices
      INTEGER  ::   ierr           ! local integer
      REAL(wp) ::   zua, zva       ! local scalars
      REAL(wp) ::   zmsk, z1_e3t   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj)     ::   zwx , zwy
      REAL(wp), DIMENSION(jpi,jpj)     ::   ztnw, ztne, ztsw, ztse
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwz
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_een : vorticity term: energy and enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         !
         SELECT CASE( kvor )                 !==  vorticity considered  ==!
         CASE ( np_COR )                           !* Coriolis (planetary vorticity)
            DO_2D_10_10
               zwz(ji,jj,jk) = ff_f(ji,jj)
            END_2D
         CASE ( np_RVO )                           !* relative vorticity
            DO_2D_10_10
               zwz(ji,jj,jk) = (  e2v(ji+1,jj  ) * pv(ji+1,jj  ,jk) - e2v(ji,jj) * pv(ji,jj,jk)    &
                  &             - e1u(ji  ,jj+1) * pu(ji  ,jj+1,jk) + e1u(ji,jj) * pu(ji,jj,jk)  ) &
                  &          * r1_e1e2f(ji,jj)
            END_2D
         CASE ( np_MET )                           !* metric term
            DO_2D_10_10
               zwz(ji,jj,jk) = ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &          - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)
            END_2D
         CASE ( np_CRV )                           !* Coriolis + relative vorticity
            DO_2D_10_10
               zwz(ji,jj,jk) = (  ff_f(ji,jj) + (  e2v(ji+1,jj  ) * pv(ji+1,jj  ,jk) - e2v(ji,jj) * pv(ji,jj,jk)    &
                  &                              - e1u(ji  ,jj+1) * pu(ji  ,jj+1,jk) + e1u(ji,jj) * pu(ji,jj,jk)  ) &
                  &                         * r1_e1e2f(ji,jj)    )
            END_2D
         CASE ( np_CME )                           !* Coriolis + metric
            DO_2D_10_10
               zwz(ji,jj,jk) = ff_f(ji,jj) + ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &                        - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)
            END_2D
         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','dyn_vor: wrong value for kvor'  )
         END SELECT
         !
         IF( ln_dynvor_msk ) THEN          !==  mask/unmask vorticity ==!
            DO_2D_10_10
               zwz(ji,jj,jk) = zwz(ji,jj,jk) * fmask(ji,jj,jk)
            END_2D
         ENDIF
      END DO
      !
      CALL lbc_lnk( 'dynvor', zwz, 'F', 1. )
      !
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !
         !                                   !==  horizontal fluxes  ==!
         zwx(:,:) = e2u(:,:) * e3u(:,:,jk,Kmm) * pu(:,:,jk)
         zwy(:,:) = e1v(:,:) * e3v(:,:,jk,Kmm) * pv(:,:,jk)

         !                                   !==  compute and add the vorticity term trend  =!
         jj = 2
         ztne(1,:) = 0   ;   ztnw(1,:) = 0   ;   ztse(1,:) = 0   ;   ztsw(1,:) = 0
         DO ji = 2, jpi          ! split in 2 parts due to vector opt.
               z1_e3t = 1._wp / e3t(ji,jj,jk,Kmm)
               ztne(ji,jj) = ( zwz(ji-1,jj  ,jk) + zwz(ji  ,jj  ,jk) + zwz(ji  ,jj-1,jk) ) * z1_e3t
               ztnw(ji,jj) = ( zwz(ji-1,jj-1,jk) + zwz(ji-1,jj  ,jk) + zwz(ji  ,jj  ,jk) ) * z1_e3t
               ztse(ji,jj) = ( zwz(ji  ,jj  ,jk) + zwz(ji  ,jj-1,jk) + zwz(ji-1,jj-1,jk) ) * z1_e3t
               ztsw(ji,jj) = ( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj-1,jk) + zwz(ji-1,jj  ,jk) ) * z1_e3t
         END DO
         DO jj = 3, jpj
            DO ji = 2, jpi   ! vector opt. ok because we start at jj = 3
               z1_e3t = 1._wp / e3t(ji,jj,jk,Kmm)
               ztne(ji,jj) = ( zwz(ji-1,jj  ,jk) + zwz(ji  ,jj  ,jk) + zwz(ji  ,jj-1,jk) ) * z1_e3t
               ztnw(ji,jj) = ( zwz(ji-1,jj-1,jk) + zwz(ji-1,jj  ,jk) + zwz(ji  ,jj  ,jk) ) * z1_e3t
               ztse(ji,jj) = ( zwz(ji  ,jj  ,jk) + zwz(ji  ,jj-1,jk) + zwz(ji-1,jj-1,jk) ) * z1_e3t
               ztsw(ji,jj) = ( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj-1,jk) + zwz(ji-1,jj  ,jk) ) * z1_e3t
            END DO
         END DO
         DO_2D_00_00
            zua = + r1_12 * r1_e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
               &                             + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
            zva = - r1_12 * r1_e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
               &                             + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
            pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zua
            pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zva
         END_2D
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
   END SUBROUTINE vor_eeT


   SUBROUTINE dyn_vor_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_vor_init  ***
      !!
      !! ** Purpose :   Control the consistency between cpp options for
      !!              tracer advection schemes
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk    ! dummy loop indices
      INTEGER ::   ioptio, ios   ! local integer
      !!
      NAMELIST/namdyn_vor/ ln_dynvor_ens, ln_dynvor_ene, ln_dynvor_enT, ln_dynvor_eeT,   &
         &                 ln_dynvor_een, nn_een_e3f   , ln_dynvor_mix, ln_dynvor_msk,   &
         &                 ln_dynvor_mix_ens_een, ln_dynvor_mix_een_ens,                 &   ! mix ENS/EEN
         &                 ln_dynvor_ens_adVO, ln_dynvor_ens_adKE, ln_dynvor_ens_adKEVO, &   ! Alternative direction parameters
         &                 ln_dynvor_ene_adVO, ln_dynvor_ene_adKE, ln_dynvor_ene_adKEVO
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_vor_init : vorticity term : read namelist and control the consistency'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      READ  ( numnam_ref, namdyn_vor, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_vor in reference namelist' )
      READ  ( numnam_cfg, namdyn_vor, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn_vor in configuration namelist' )
      IF(lwm) WRITE ( numond, namdyn_vor )
      !
      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*) '   Namelist namdyn_vor : choice of the vorticity term scheme'
         WRITE(numout,*) '      enstrophy conserving scheme                    ln_dynvor_ens = ', ln_dynvor_ens
         WRITE(numout,*) '      f-point energy conserving scheme               ln_dynvor_ene = ', ln_dynvor_ene
         WRITE(numout,*) '      t-point energy conserving scheme               ln_dynvor_enT = ', ln_dynvor_enT
         WRITE(numout,*) '      energy conserving scheme  (een using e3t)      ln_dynvor_eeT = ', ln_dynvor_eeT
         WRITE(numout,*) '      enstrophy and energy conserving scheme         ln_dynvor_een = ', ln_dynvor_een
         WRITE(numout,*) '         e3f = averaging /4 (=0) or /sum(tmask) (=1)    nn_een_e3f = ', nn_een_e3f
         WRITE(numout,*) '      mixed enstrophy/energy conserving scheme       ln_dynvor_mix = ', ln_dynvor_mix
         WRITE(numout,*) '      masked (=T) or unmasked(=F) vorticity          ln_dynvor_msk = ', ln_dynvor_msk
      ENDIF

!!an      IF( ln_dynvor_msk )   CALL ctl_stop( 'dyn_vor_init:   masked vorticity is not currently not available')

!!gm  this should be removed when choosing a unique strategy for fmask at the coast
      ! If energy, enstrophy or mixed advection of momentum in vector form change the value for masks
      ! at angles with three ocean points and one land point
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '      change fmask value in the angles (T)           ln_vorlat = ', ln_vorlat
      IF( ln_vorlat .AND. ( ln_dynvor_ene .OR. ln_dynvor_ens .OR. ln_dynvor_mix ) ) THEN
         DO_3D_10_10( 1, jpk )
            IF(    tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk)              &
               & + tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk) == 3._wp )   fmask(ji,jj,jk) = 1._wp
         END_3D
         !
         CALL lbc_lnk( 'dynvor', fmask, 'F', 1._wp )      ! Lateral boundary conditions on fmask
         !
      ENDIF
!!gm end

      ioptio = 0                     ! type of scheme for vorticity (set nvor_scheme)
      IF( ln_dynvor_ens )        THEN   ;   ioptio = ioptio + 1     ;   nvor_scheme = np_ENS        ;   ENDIF
      IF( ln_dynvor_ens_adVO )   THEN   ;   ioptio = ioptio + 1     ;   nvor_scheme = np_ENS_adVO   ;   ENDIF
      IF( ln_dynvor_ens_adKE )   THEN   ;   ioptio = ioptio + 1     ;   nvor_scheme = np_ENS_adKE   ;   ENDIF
      IF( ln_dynvor_ens_adKEVO ) THEN   ;   ioptio = ioptio + 1     ;   nvor_scheme = np_ENS_adKEVO ;   ENDIF
      IF( ln_dynvor_ene )        THEN   ;   ioptio = ioptio + 1     ;   nvor_scheme = np_ENE        ;   ENDIF
      IF( ln_dynvor_ene_adVO )   THEN   ;   ioptio = ioptio + 1     ;   nvor_scheme = np_ENE_adVO   ;   ENDIF
      IF( ln_dynvor_ene_adKE )   THEN   ;   ioptio = ioptio + 1     ;   nvor_scheme = np_ENE_adKE   ;   ENDIF
      IF( ln_dynvor_ene_adKEVO ) THEN   ;   ioptio = ioptio + 1     ;   nvor_scheme = np_ENE_adKEVO ;   ENDIF
      IF( ln_dynvor_enT ) THEN   ;   ioptio = ioptio + 1   ;   nvor_scheme = np_ENT   ;   ENDIF
      IF( ln_dynvor_eeT ) THEN   ;   ioptio = ioptio + 1   ;   nvor_scheme = np_EET   ;   ENDIF
      IF( ln_dynvor_een ) THEN   ;   ioptio = ioptio + 1   ;   nvor_scheme = np_EEN   ;   ENDIF
      IF( ln_dynvor_mix ) THEN   ;   ioptio = ioptio + 1   ;   nvor_scheme = np_MIX   ;   ENDIF
      IF( ln_dynvor_mix_ens_een ) THEN   ;   ioptio = ioptio + 1   ;   nvor_scheme = np_MIX1   ;   ENDIF
      IF( ln_dynvor_mix_een_ens ) THEN   ;   ioptio = ioptio + 1   ;   nvor_scheme = np_MIX2   ;   ENDIF
      !
      IF( ioptio /= 1 ) CALL ctl_stop( ' use ONE and ONLY one vorticity scheme' )
      !
      IF(lwp) WRITE(numout,*)        ! type of calculated vorticity (set ncor, nrvm, ntot)
      ncor = np_COR                       ! planetary vorticity
      SELECT CASE( n_dynadv )
      CASE( np_LIN_dyn )
         IF(lwp) WRITE(numout,*) '   ==>>>   linear dynamics : total vorticity = Coriolis'
         nrvm = np_COR        ! planetary vorticity
         ntot = np_COR        !     -         -
      CASE( np_VEC_c2  )
         IF(lwp) WRITE(numout,*) '   ==>>>   vector form dynamics : total vorticity = Coriolis + relative vorticity'
         nrvm = np_RVO        ! relative vorticity
         ntot = np_CRV        ! relative + planetary vorticity
      CASE( np_FLX_c2 , np_FLX_ubs, np_FLX_c4 )
         IF(lwp) WRITE(numout,*) '   ==>>>   flux form dynamics : total vorticity = Coriolis + metric term'
         nrvm = np_MET        ! metric term
         ntot = np_CME        ! Coriolis + metric term
         !
         SELECT CASE( nvor_scheme )    ! pre-computed gradients for the metric term:
         CASE( np_ENT )                      !* T-point metric term :   pre-compute di(e2u)/2 and dj(e1v)/2
            ALLOCATE( di_e2u_2(jpi,jpj), dj_e1v_2(jpi,jpj) )
            DO_2D_00_00
               di_e2u_2(ji,jj) = ( e2u(ji,jj) - e2u(ji-1,jj  ) ) * 0.5_wp
               dj_e1v_2(ji,jj) = ( e1v(ji,jj) - e1v(ji  ,jj-1) ) * 0.5_wp
            END_2D
            CALL lbc_lnk_multi( 'dynvor', di_e2u_2, 'T', -1. , dj_e1v_2, 'T', -1. )   ! Lateral boundary conditions
            !
         CASE DEFAULT                        !* F-point metric term :   pre-compute di(e2u)/(2*e1e2f) and dj(e1v)/(2*e1e2f)
            ALLOCATE( di_e2v_2e1e2f(jpi,jpj), dj_e1u_2e1e2f(jpi,jpj) )
            DO_2D_10_10
               di_e2v_2e1e2f(ji,jj) = ( e2v(ji+1,jj  ) - e2v(ji,jj) )  * 0.5 * r1_e1e2f(ji,jj)
               dj_e1u_2e1e2f(ji,jj) = ( e1u(ji  ,jj+1) - e1u(ji,jj) )  * 0.5 * r1_e1e2f(ji,jj)
            END_2D
            CALL lbc_lnk_multi( 'dynvor', di_e2v_2e1e2f, 'F', -1. , dj_e1u_2e1e2f, 'F', -1. )   ! Lateral boundary conditions
         END SELECT
         !
      END SELECT

      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         SELECT CASE( nvor_scheme )

         CASE( np_ENS    )   ;   WRITE(numout,*) '   ==>>>   enstrophy conserving scheme (ENS)'
         CASE( np_ENS_adVO ) ;   WRITE(numout,*) '   ==>>>   AD enstrophy conserving scheme (ENS_adVO) on vorticity only'
         CASE( np_ENS_adKE ) ;   WRITE(numout,*) '   ==>>>   AD enstrophy conserving scheme (ENS_adKE) on kinetic energy only'
         CASE( np_ENS_adKEVO ) ;   WRITE(numout,*) '   ==>>>   AD enstrophy conserving scheme (ENS_adKEVO) on kinetic energy and vorticity'

         CASE( np_ENE    )   ;   WRITE(numout,*) '   ==>>>   energy conserving scheme (Coriolis at F-points) (ENE)'
         CASE( np_ENE_adVO ) ;   WRITE(numout,*) '   ==>>>   AD energy conserving scheme (Coriolis at F-points) (ENE_adVO) on vorticity only'
         CASE( np_ENE_adKE ) ;   WRITE(numout,*) '   ==>>>   AD energy conserving scheme (Coriolis at F-points) (ENE_adKE) on kinetic energy only'
         CASE( np_ENE_adKEVO ) ;   WRITE(numout,*) '   ==>>>   AD energy conserving scheme (Coriolis at F-points) (ENE_adKEVO) on kinetic energy and vorticity'

         CASE( np_ENT    )   ;   WRITE(numout,*) '   ==>>>   energy conserving scheme (Coriolis at T-points) (ENT)'
         CASE( np_EET    )   ;   WRITE(numout,*) '   ==>>>   energy conserving scheme (EEN scheme using e3t) (EET)'
         CASE( np_EEN    )   ;   WRITE(numout,*) '   ==>>>   energy and enstrophy conserving scheme (EEN)'
         CASE( np_MIX    )   ;   WRITE(numout,*) '   ==>>>   mixed enstrophy/energy conserving scheme (MIX)'
         CASE( np_MIX1    )   ;   WRITE(numout,*) '   ==>>>   mixed enstrophy/energy conserving scheme (ENS/EEN)'
         CASE( np_MIX2    )   ;   WRITE(numout,*) '   ==>>>   mixed enstrophy/energy conserving scheme (EEN/ENS)'
         END SELECT
      ENDIF
      !
   END SUBROUTINE dyn_vor_init

   !!==============================================================================
END MODULE dynvor
