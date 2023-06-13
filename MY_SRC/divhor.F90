MODULE divhor
   !!==============================================================================
   !!                       ***  MODULE  divhor  ***
   !! Ocean diagnostic variable : now horizontal divergence
   !!==============================================================================
   !! History :  1.0  ! 2002-09  (G. Madec, E. Durand)  Free form, F90
   !!             -   ! 2005-01  (J. Chanut) Unstructured open boundaries
   !!             -   ! 2003-08  (G. Madec)  merged of cur and div, free form, F90
   !!             -   ! 2005-01  (J. Chanut, A. Sellar) unstructured open boundaries
   !!            3.3  ! 2010-09  (D.Storkey and E.O'Dea) bug fixes for BDY module
   !!             -   ! 2010-10  (R. Furner, G. Madec) runoff and cla added directly here
   !!            3.7  ! 2014-01  (G. Madec) suppression of velocity curl from in-core memory
   !!             -   ! 2014-12  (G. Madec) suppression of cross land advection option
   !!             -   ! 2015-10  (G. Madec) add velocity and rnf flag in argument of div_hor
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   div_hor    : Compute the horizontal divergence field
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce, ONLY : ln_rnf      ! river runoff
   USE sbcrnf , ONLY : sbc_rnf_div ! river runoff
   USE isf_oce, ONLY : ln_isf      ! ice shelf
   USE isfhdiv, ONLY : isf_hdiv    ! ice shelf
#if defined key_asminc
   USE asminc          ! Assimilation increment
#endif
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   div_hor    ! routine called by step.F90 and istate.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: divhor.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE div_hor( kt, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE div_hor  ***
      !!
      !! ** Purpose :   compute the horizontal divergence at now time-step
      !!
      !! ** Method  :   the now divergence is computed as :
      !!         hdiv = 1/(e1e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
      !!      and correct with runoff inflow (div_rnf) and cross land flow (div_cla)
      !!
      !! ** Action  : - update hdiv, the now horizontal divergence
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt        ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm  ! ocean time level indices
      !
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::   zraur, zdep   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) :: ztmp, z2d

      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('div_hor')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'div_hor : horizontal velocity divergence '
         IF(lwp) WRITE(numout,*) '~~~~~~~   '
         hdiv(:,:,:) = 0._wp    ! initialize hdiv for the halos at the first time step
      ENDIF

      ! e3 stays penalised for consistency in continuity equation [1/s]
      ! ssh then becomes penalised implicitly
      DO_3D_00_00( 1, jpkm1 )
         hdiv(ji,jj,jk) = (   e2u(ji  ,jj) * e3u(ji  ,jj  ,jk, Kmm ) * uu(ji  ,jj,jk,Kmm)      &
            &               - e2u(ji-1,jj) * e3u(ji-1,jj  ,jk, Kmm ) * uu(ji-1,jj,jk,Kmm)      &
            &               + e1v(ji,jj  ) * e3v(ji  ,jj  ,jk, Kmm ) * vv(ji,jj  ,jk,Kmm)      &
            &               - e1v(ji,jj-1) * e3v(ji  ,jj-1,jk, Kmm ) * vv(ji,jj-1,jk,Kmm)  )   &
            &            * r1_e1e2t(ji,jj) / e3t(ji  ,jj  ,jk, Kmm)
      END_3D
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN
         IF( nbondi == -1 .OR. nbondi == 2 )   hdiv(   2   ,  :   ,:) = 0._wp      ! west
         IF( nbondi ==  1 .OR. nbondi == 2 )   hdiv( nlci-1,  :   ,:) = 0._wp      ! east
         IF( nbondj == -1 .OR. nbondj == 2 )   hdiv(   :   ,  2   ,:) = 0._wp      ! south
         IF( nbondj ==  1 .OR. nbondj == 2 )   hdiv(   :   ,nlcj-1,:) = 0._wp      ! north
      ENDIF
#endif
      !
      IF( ln_rnf )   CALL sbc_rnf_div( hdiv, Kmm )                     !==  runoffs    ==!   (update hdiv field)
      !
#if defined key_asminc
      IF( ln_sshinc .AND. ln_asmiau )   CALL ssh_asm_div( kt, Kbb, Kmm, hdiv )   !==  SSH assimilation  ==!   (update hdiv field)
      !
#endif
      !
      IF( ln_isf )                      CALL isf_hdiv( kt, Kmm, hdiv )           !==  ice shelf         ==!   (update hdiv field)
      !
      CALL lbc_lnk( 'divhor', hdiv, 'T', 1. )   !   (no sign change)
      !
      IF( ln_timing )   CALL timing_stop('div_hor')
      !
   END SUBROUTINE div_hor

   !!======================================================================
END MODULE divhor
