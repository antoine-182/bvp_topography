MODULE dynldf_lap_blp
   !!======================================================================
   !!                   ***  MODULE  dynldf_lap_blp  ***
   !! Ocean dynamics:  lateral viscosity trend (laplacian and bilaplacian)
   !!======================================================================
   !! History : 3.7  ! 2014-01  (G. Madec, S. Masson)  Original code, re-entrant laplacian
   !!           4.0  ! 2020-04  (A. Nasser, G. Madec)  Add symmetric mixing tensor
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf_lap   : update the momentum trend with the lateral viscosity using an iso-level   laplacian operator
   !!   dyn_ldf_blp   : update the momentum trend with the lateral viscosity using an iso-level bilaplacian operator
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE ldfslp         ! iso-neutral slopes
   USE zdf_oce        ! ocean vertical physics
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library
   !
!!an
  USE usrdef_nam , ONLY : nn_dynldf_lap_typ, nn_dynldf_flx, ln_mir_ldf    ! use laplacian parameter
   !
   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_ldf_lap  ! called by dynldf.F90
   PUBLIC dyn_ldf_blp  ! called by dynldf.F90
!!anSYM
   INTEGER, PUBLIC, PARAMETER ::   np_dynldf_lap_rot       = 1         ! div-rot   laplacian
   INTEGER, PUBLIC, PARAMETER ::   np_dynldf_lap_rot_noh   = 10         ! div-rot   laplacian
   INTEGER, PUBLIC, PARAMETER ::   np_dynldf_lap_sym       = 2         ! symmetric laplacian (Griffies&Hallberg 2000)
   INTEGER, PUBLIC, PARAMETER ::   np_dynldf_lap_sym_noh   = 20         ! symmetric laplacian (Griffies&Hallberg 2000)
   INTEGER, PUBLIC, PARAMETER ::   np_dynldf_lap_con       = 3         ! conventional laplacian (cartesian)
   INTEGER, PUBLIC, PARAMETER ::   np_dynldf_lap_con_noh   = 30         ! conventional laplacian (cartesian)

   !INTEGER, PUBLIC, PARAMETER ::   nn_dynldf_lap_typ = 1         ! choose type of laplacian (ideally from namelist)
!!anSYM
   !! * Substitutions
#  include "do_loop_substitute.h90"
!!st21
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynldf_lap_blp.F90 13416 2020-08-20 10:10:55Z gm $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf_lap( kt, Kbb, Kmm, pu, pv, pu_rhs, pv_rhs, kpass )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_lap  ***
      !!
      !! ** Purpose :   Compute the before horizontal momentum diffusive
      !!      trend and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The Laplacian operator apply on horizontal velocity is
      !!      writen as :   grad_h( ahmt div_h(U )) - curl_h( ahmf curl_z(U) )
      !!      writen as :   grad_h( ahmt div_h(U )) - curl_h( ahmf curl_z(U) )
      !!
      !! ** Action : - pu_rhs, pv_rhs increased by the harmonic operator applied on pu, pv.
      !!
      !! Reference : S.Griffies, R.Hallberg 2000 Mon.Wea.Rev., DOI:/
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kbb, Kmm   ! ocean time level indices
      INTEGER                         , INTENT(in   ) ::   kpass      ! =1/2 first or second passage
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pu, pv     ! before velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs   ! velocity trend   [m/s2]
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zsign        ! local scalars
      REAL(wp) ::   zua, zva, ze1     ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zcur, zdiv
      REAL(wp), DIMENSION(jpi,jpj) ::   zten, zshe   ! tension (diagonal) and shearing (anti-diagonal) terms
      REAL(wp), DIMENSION(jpi,jpj) ::   zxu, zyu, zxv, zyv  ! shears and tension (Conventionnal)
      !
      REAL(wp), DIMENSION(jpi,jpj)     ::   ze1e2f, zr1_e1e2f, zr1_e1e2t     ! penalised scale factors
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ze3f                   !
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zu   , zv      ! temporary
      REAL(wp), DIMENSION(jpi,jpj,jpk,3) :: ze3t, ze3u, ze3v             ! 4D workspace penalised
      !!----------------------------------------------------------------------
      !
!!anSYM TO BE ADDED : reading of laplacian operator (ln_dynldf_lap_typ -> to be written nn_) shall be added in dyn_ldf_init
!!                 as the writing
!!                 and an integer as np_dynldf_lap for instance taken as argument by dyn_ldf_lap call in dyn_ldf
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ldf : iso-level harmonic (laplacian) operator, pass=', kpass
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '                  nn_dynldf_lap_typ = ', nn_dynldf_lap_typ
         SELECT CASE( nn_dynldf_lap_typ )             ! print the choice of operator
         CASE( np_dynldf_lap_rot )   ;   WRITE(numout,*) '   ==>>>   div-rot   laplacian'
         CASE( np_dynldf_lap_sym )   ;   WRITE(numout,*) '   ==>>>   symmetric laplacian (covariant form)'
         CASE( np_dynldf_lap_con)    ;   WRITE(numout,*) '   ==>>>   conventional laplacian '
         CASE( np_dynldf_lap_rot_noh )   ;   WRITE(numout,*) '   ==>>>   div-rot   laplacian (no h)'
         CASE( np_dynldf_lap_sym_noh )   ;   WRITE(numout,*) '   ==>>>   symmetric laplacian (no h)'
         CASE( np_dynldf_lap_con_noh )   ;   WRITE(numout,*) '   ==>>>   conventional laplacian (no h)'
         END SELECT
      ENDIF
      !
      IF( kpass == 1 ) THEN   ;   zsign =  1._wp      ! bilaplacian operator require a minus sign
      ELSE                    ;   zsign = -1._wp      !  (eddy viscosity coef. >0)
      ENDIF
      !
      zu(:,:,:) = pu(:,:,:)*umask(:,:,:) ; zv(:,:,:) = pv(:,:,:)*vmask(:,:,:) ! masked arrays (pu et pv déjà maské en fait)
      !
      IF ( ln_mir_ldf ) THEN   ! free-slip mirror boundary condition, works only at 45°
          !
          DO jk = 1, jpkm1
            DO_2D_10_10
                  IF ( (     tmask(ji  ,jj+1,1) + tmask(ji+1,jj+1,1)  &
                       &   + tmask(ji  ,jj  ,1) + tmask(ji+1,jj  ,1)  ) == 3._wp) THEN ! salient angle (45°)
                      !
                      IF      ( tmask(ji  ,jj  ,1)==0._wp ) THEN   ! South
                       zu(ji,jj,jk) = - pv(ji+1,jj  ,jk)
                       zv(ji,jj,jk) = - pu(ji  ,jj+1,jk)
                      ELSE IF ( tmask(ji+1,jj  ,1)==0._wp ) THEN   ! East
                        zu(ji  ,jj,jk) = + pv(ji,jj  ,jk)
                        zv(ji+1,jj,jk) = + pu(ji,jj+1,jk)
                      ELSE IF ( tmask(ji+1,jj+1,1)==0._wp ) THEN   ! North
                        zu(ji  ,jj+1,jk) = - pv(ji,jj,jk)
                        zv(ji+1,jj  ,jk) = - pu(ji,jj,jk)
                      ELSE IF ( tmask(ji  ,jj+1,1)==0._wp ) THEN   ! West
                        zu(ji,jj+1,jk) = + pv(ji+1,jj  ,jk)
                        zv(ji,jj  ,jk) = + pu(ji  ,jj,jk)
                      END IF
                      !
                  END IF
                !
            END_2D
          END DO
          CALL lbc_lnk_multi( 'dynldf_lap_blp', zu , 'U', 1._wp,   &
                    &                           zv , 'V', 1._wp    )
          !
      END IF
      !
      SELECT CASE( nn_dynldf_lap_typ )
         !
         CASE ( np_dynldf_lap_rot )       !==  Vorticity-Divergence form  ==!
           !
            DO jk = 1, jpkm1                                 ! Horizontal slab
               !
               ! !!an 01 : 2 -> jpi ; 00 : 2 -> jpim1 ; 0 = interieur, 1 = vers l'exterieur
               !
               DO_2D_01_01
               !                                      ! ahm * e3 * curl  (computed from 1 to jpim1/jpjm1)
!!gm note that ahmf has already been multiplied by fmask
!!an en discret, pour l'energie, pour démontrer la conservation il faut ces e3f/e3u/e3v (voir annex)
!!an zcurl se définit en ji=1
!!an zdiv se définit en ji=2
                zcur(ji-1,jj-1) =  &
                   &       ahmf(ji-1,jj-1,jk) * e3f(ji-1,jj-1,jk) * r1_e1e2f(ji-1,jj-1)      &
                   &  * (  e2v(ji  ,jj-1) * zv(ji  ,jj-1,jk) - e2v(ji-1,jj-1) * zv(ji-1,jj-1,jk)  &
                   &     - e1u(ji-1,jj  ) * zu(ji-1,jj  ,jk) + e1u(ji-1,jj-1) * zu(ji-1,jj-1,jk)  )
            !                                      ! ahm * div        (computed from 2 to jpi/jpj)
!!gm note that ahmt has already been multiplied by tmask
                zdiv(ji,jj)     = ahmt(ji,jj,jk) * r1_e1e2t(ji,jj) / e3t(ji,jj,jk,Kbb)                                         &
                   &     * (  e2u(ji,jj)*e3u(ji,jj,jk,Kbb) * zu(ji,jj,jk) - e2u(ji-1,jj)*e3u(ji-1,jj,jk,Kbb) * zu(ji-1,jj,jk)  &
                   &        + e1v(ji,jj)*e3v(ji,jj,jk,Kbb) * zv(ji,jj,jk) - e1v(ji,jj-1)*e3v(ji,jj-1,jk,Kbb) * zv(ji,jj-1,jk)  )
                   END_2D
               !
               DO_2D_00_00
                  pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign * (                                             &
                     &              - ( zcur(ji  ,jj) - zcur(ji,jj-1) ) * r1_e2u(ji,jj) / e3u(ji,jj,jk,Kmm)   &
                     &              + ( zdiv(ji+1,jj) - zdiv(ji,jj  ) ) * r1_e1u(ji,jj)                       )
                     !
                  pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign * (                                             &
                     &                ( zcur(ji,jj  ) - zcur(ji-1,jj) ) * r1_e1v(ji,jj) / e3v(ji,jj,jk,Kmm)   &
                     &              + ( zdiv(ji,jj+1) - zdiv(ji  ,jj) ) * r1_e2v(ji,jj)                       )
               END_2D
               !
            END DO                                           !   End of slab
            !
          CASE ( np_dynldf_lap_rot_noh )       !==  Vorticity-Divergence form  ==!
            !
             DO jk = 1, jpkm1                                 ! Horizontal slab
                !
                DO_2D_01_01
                 !                                      !-- without e3
                 zcur(ji-1,jj-1) =  &
                    &       ahmf(ji-1,jj-1,jk) * r1_e1e2f(ji-1,jj-1)      &
                    &  * (  e2v(ji  ,jj-1) * zv(ji  ,jj-1,jk) - e2v(ji-1,jj-1) * zv(ji-1,jj-1,jk)  &
                    &     - e1u(ji-1,jj  ) * zu(ji-1,jj  ,jk) + e1u(ji-1,jj-1) * zu(ji-1,jj-1,jk)  )
                 !
                 zdiv(ji,jj)     = ahmt(ji,jj,jk) * r1_e1e2t(ji,jj)                                         &
                    &     * (  e2u(ji,jj) * zu(ji,jj,jk) - e2u(ji-1,jj) * zu(ji-1,jj,jk)  &
                    &        + e1v(ji,jj) * zv(ji,jj,jk) - e1v(ji,jj-1) * zv(ji,jj-1,jk)  )
                END_2D
                !
                DO_2D_00_00
                   pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign * (                                             &
                      &              - ( zcur(ji  ,jj) - zcur(ji,jj-1) ) * r1_e2u(ji,jj)    &
                      &              + ( zdiv(ji+1,jj) - zdiv(ji,jj  ) ) * r1_e1u(ji,jj)                       )
                      !
                   pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign * (                                             &
                      &                ( zcur(ji,jj  ) - zcur(ji-1,jj) ) * r1_e1v(ji,jj)    &
                      &              + ( zdiv(ji,jj+1) - zdiv(ji  ,jj) ) * r1_e2v(ji,jj)                       )
                END_2D

                !
             END DO                                           !   End of slab
             !
         CASE ( np_dynldf_lap_sym )       !==  Symmetric form  ==!   (Griffies&Hallberg 2000)
            !            !
            !
            DO jk = 1, jpkm1                                 ! Horizontal slab
               !
               DO_2D_01_01
                  !                                      ! shearing stress component (F-point)   NB : ahmf has already been multiplied by fmask
                  zshe(ji-1,jj-1) = ahmf(ji-1,jj-1,jk) * e3f(ji-1,jj-1,jk)                                          &
                     &     * (    e1f(ji-1,jj-1)    * r1_e2f(ji-1,jj-1)                                             &
                     &         * ( zu(ji-1,jj  ,jk) * r1_e1u(ji-1,jj  )  - zu(ji-1,jj-1,jk) * r1_e1u(ji-1,jj-1) )   &
                     &         +  e2f(ji-1,jj-1)    * r1_e1f(ji-1,jj-1)                                             &
                     &         * ( zv(ji  ,jj-1,jk) * r1_e2v(ji  ,jj-1)  - zv(ji-1,jj-1,jk) * r1_e2v(ji-1,jj-1) )   )
                  !                                      ! tension stress component (T-point)   NB : ahmt has already been multiplied by tmask
                  zten(ji,jj)    = ahmt(ji,jj,jk) *  e3t(ji,jj,jk,Kmm)                                                      &
                     &     * (    e2t(ji,jj)    * r1_e1t(ji,jj)                                         &
                     &         * ( zu(ji,jj,jk) * r1_e2u(ji,jj)  - zu(ji-1,jj,jk) * r1_e2u(ji-1,jj) )   &
                     &         -  e1t(ji,jj)    * r1_e2t(ji,jj)                                         &
                     &         * ( zv(ji,jj,jk) * r1_e1v(ji,jj)  - zv(ji,jj-1,jk) * r1_e1v(ji,jj-1) )   )
               END_2D
               !
               DO_2D_00_00
                  pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)        &
                     &    * (   (   zten(ji+1,jj  ) * e2t(ji+1,jj  )*e2t(ji+1,jj  )                        &
                     &            - zten(ji  ,jj  ) * e2t(ji  ,jj  )*e2t(ji  ,jj  )    ) * r1_e2u(ji,jj)   &
                     &        + (   zshe(ji  ,jj  ) * e1f(ji  ,jj  )*e1f(ji  ,jj  )                        &
                     &            - zshe(ji  ,jj-1) * e1f(ji  ,jj-1)*e1f(ji  ,jj-1)    ) * r1_e1u(ji,jj) )
                  !
                  pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)        &
                     &    * (   (   zshe(ji  ,jj  ) * e2f(ji  ,jj  )*e2f(ji  ,jj  )                        &
                     &            - zshe(ji-1,jj  ) * e2f(ji-1,jj  )*e2f(ji-1,jj  )  ) * r1_e2v(ji,jj)     &
                     &        - (   zten(ji  ,jj+1) * e1t(ji  ,jj+1)*e1t(ji  ,jj+1)                        &
                     &            - zten(ji  ,jj  ) * e1t(ji  ,jj  )*e1t(ji  ,jj  )  ) * r1_e1v(ji,jj) )
                   !
               END_2D
               !
            END DO                                           !   End of slab
            !
          CASE ( np_dynldf_lap_sym_noh )       !==  Symmetric form  ==!
             !            !
             !
             DO jk = 1, jpkm1                                 ! Horizontal slab
                !
                DO_2D_01_01
                   !                                      ! Ds - shearing stress component (F-point)   NB : ahmf has already been multiplied by fmask
                   zshe(ji-1,jj-1) = ahmf(ji-1,jj-1,jk)                           &
                      &         *  (  ( zu(ji-1,jj  ,jk) - zu(ji-1,jj-1,jk) ) * r1_e2f(ji-1,jj-1)   &
                      &            +  ( zv(ji  ,jj-1,jk) - zv(ji-1,jj-1,jk) ) * r1_e1f(ji-1,jj-1)   )
                   !                                      ! Dt - tension stress component (T-point)   NB : ahmt has already been multiplied by tmask
                   zten(ji,jj)    = ahmt(ji,jj,jk)                          &
                      &         *  ( ( zu(ji,jj,jk) - zu(ji-1,jj  ,jk)  ) * r1_e1t(ji,jj)   &
                      &            - ( zv(ji,jj,jk) - zv(ji  ,jj-1,jk)  ) * r1_e2t(ji,jj)   )
                END_2D
                !
                DO_2D_00_00
                   pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign   &
                      &    * (   (   zten(ji+1,jj  ) - zten(ji  ,jj  )   ) * r1_e1u(ji,jj)  &
                      &        + (   zshe(ji  ,jj  ) - zshe(ji  ,jj-1)   ) * r1_e2u(ji,jj)  )
                   !
                   pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign   &
                      &    * (   (   zshe(ji  ,jj  ) - zshe(ji-1,jj  ) ) * r1_e1v(ji,jj) &
                      &        - (   zten(ji  ,jj+1) - zten(ji  ,jj  ) ) * r1_e2v(ji,jj) )
                    !
                    !
                END_2D
                !
                !
             END DO                                           !   End of slab
             !
         CASE ( np_dynldf_lap_con )       !==  Conventionnal form  ==!
            !
            !
            DO jk = 1, jpkm1                                 ! Horizontal slab
               !
               DO_2D_01_01
                  !                                      ! shearing stress component (F-point)   NB : ahmf has already been multiplied by fmask
                 zyu(ji-1,jj-1) = ahmf(ji-1,jj-1,jk) * e1f(ji-1,jj-1) * e3f(ji-1,jj-1,jk)     &
                    &           * r1_e2f(ji-1,jj-1) * ( pu(ji-1,jj  ,jk) - pu (ji-1,jj-1,jk) )
                 zxv(ji-1,jj-1) = ahmf(ji-1,jj-1,jk) * e2f(ji-1,jj-1) * e3f(ji-1,jj-1,jk)     &
                    &           * r1_e1f(ji-1,jj-1) * ( pv(ji  ,jj-1,jk) - pv (ji-1,jj-1,jk) )
                  !                                      ! tension stress component (T-point)   NB : ahmt has already been multiplied by tmask
                 zxu(ji,jj)     = ahmt(ji,jj,jk) *  e2t(ji,jj   ) * e3t(ji  ,jj  ,jk,Kmm)          &
                    &           * r1_e1t(ji,jj)  * ( pu(ji,jj,jk) - pu (ji-1,jj  ,jk) )
                 zyv(ji,jj)     = ahmt(ji,jj,jk) *  e1t(ji,jj   ) * e3t(ji  ,jj,  jk,Kmm)          &
                    &           * r1_e2t(ji,jj)  * ( pv(ji,jj,jk) - pv (ji  ,jj-1,jk) )
               END_2D
               !
               DO_2D_00_00
                  pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)  &
                     &    * (   zxu(ji+1,jj  ) - zxu(ji  ,jj  )                                      &
                     &        + zyu(ji  ,jj  ) - zyu(ji  ,jj-1)                                      )
                  !
                  pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)  &
                     &    * (   zxv(ji  ,jj  )  - zxv(ji-1,jj  )                                     &
                     &        + zyv(ji  ,jj+1)  - zyv(ji  ,jj  )                                     )
                   !
               END_2D
               !
               !
            END DO                                           !   End of slab
            !
          CASE ( np_dynldf_lap_con_noh )       !==  Conventionnal form  ==!
             !
             !
             DO jk = 1, jpkm1                                 ! Horizontal slab
                !
                DO_2D_01_01
                   !                                      ! shearing stress component (F-point)   NB : ahmf has already been multiplied by fmask
                  zyu(ji-1,jj-1) = ahmf(ji-1,jj-1,jk)                                            &
                     &           * r1_e2f(ji-1,jj-1) * ( pu(ji-1,jj  ,jk) - pu(ji-1,jj-1,jk) )
                  zxv(ji-1,jj-1) = ahmf(ji-1,jj-1,jk)                                            &
                     &           * r1_e1f(ji-1,jj-1) * ( pv(ji  ,jj-1,jk) - pv(ji-1,jj-1,jk) )
                   !                                      ! tension stress component (T-point)   NB : ahmt has already been multiplied by tmask
                  zxu(ji,jj)     = ahmt(ji,jj,jk)                                         &
                     &           * r1_e1t(ji,jj) * ( pu(ji,jj,jk) - pu(ji-1,jj  ,jk) )
                  zyv(ji,jj)     = ahmt(ji,jj,jk)                                         &
                     &           * r1_e2t(ji,jj) * ( pv(ji,jj,jk) - pv(ji  ,jj-1,jk) )
                END_2D
                !
                DO_2D_00_00
                   pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign * r1_e1e2u(ji,jj)   &
                      &    * (   zxu(ji+1,jj  ) * e2t(ji+1,jj  )         &  ! multiplier dans la zxu (gagne 2*)
                      &        - zxu(ji  ,jj  ) * e2t(ji  ,jj  )         &
                      &        + zyu(ji  ,jj  ) * e1f(ji  ,jj  )         &
                      &        - zyu(ji  ,jj-1) * e1f(ji  ,jj-1)         )
                   !
                   pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign * r1_e1e2v(ji,jj)  &
                      &    * (   zxv(ji  ,jj  ) * e2f(ji  ,jj  )     &
                      &        - zxv(ji-1,jj  ) * e2f(ji-1,jj  )     &
                      &        + zyv(ji  ,jj+1) * e1t(ji  ,jj+1)     &
                      &        - zyv(ji  ,jj  ) * e1t(ji  ,jj  )     )
                    !
                END_2D
                !
             END DO                                           !   End of slab
             !
         CASE DEFAULT                                     ! error
            CALL ctl_stop('STOP','dyn_ldf_lap: wrong value for nn_dynldf_lap_typ'  )
         END SELECT
         !
       IF ( nn_dynldf_flx == 1 ) THEN
         SELECT CASE( nn_dynldf_lap_typ )
            !
            CASE ( np_dynldf_lap_sym     )
              DO_2D_00_00
              ! rajouter un paramètre pour la compensation, pour choisir avec quel flux diffusif on compense
              ! devrait converger vers le même flux car e3u/e3t -> 1
                  IF ( (     tmask(ji  ,jj+1,1) + tmask(ji+1,jj+1,1)  &
                       &   + tmask(ji  ,jj  ,1) + tmask(ji+1,jj  ,1)  ) == 3._wp) THEN ! salient angle
                      !
                      IF      ( tmask(ji  ,jj  ,1)==0._wp ) THEN   ! South
                        !! compensation
                        pu_rhs(ji,jj+1,jk)  = pu_rhs(ji,jj+1,jk) - 1._wp * zsign * r1_e1e2u(ji,jj+1) / e3u(ji,jj+1,jk,Kmm)  &  ! ++--
                             &            * e2t(ji+1,jj+1) *e2t(ji+1,jj+1) * e3t(ji+1,jj+1,jk,Kmm)  * r1_e2u(ji,jj+1)       &
                             &            * ahmt(ji+1,jj+1,jk)   * e1t(ji+1,jj+1) * r1_e2t(ji+1,jj+1) *  zv(ji+1,jj,jk) * r1_e1v(ji+1,jj) ! -Dt(i+1,j+1)
                        pv_rhs(ji+1,jj,jk)= pv_rhs(ji+1,jj,jk) - 1._wp * zsign * r1_e1e2v(ji+1,jj) / e3v(ji+1,jj,jk,Kmm)    & ! -++-
                             &            * e1t(ji+1,jj+1) * e1t(ji+1,jj+1) * e3t(ji+1,jj+1,jk,Kmm) * r1_e1v(ji+1,jj)       &
                             &            * ahmt(ji+1,jj+1,jk) * e2t(ji+1,jj+1) * r1_e1t(ji+1,jj+1) * zu(ji,jj+1,jk) * r1_e2u(ji,jj+1)   ! +Dt(i+1,j+1)
                       !! divrot flux
                         pu_rhs(ji,jj+1,jk)  = pu_rhs(ji,jj+1,jk) - 1._wp * zsign     * r1_e1u(ji,jj+1)                    &
                              &            * ahmt(ji+1,jj+1,jk)   * r1_e1e2t(ji+1,jj+1) / e3t(ji+1,jj+1,jk, Kbb)           &
                              &                                 * zv(ji+1,jj,jk)    * e1v(ji,jj+1) * e3v(ji+1,jj,jk, Kbb)
                         pv_rhs(ji+1,jj,jk)= pv_rhs(ji+1,jj,jk) - 1._wp * zsign     * r1_e2v(ji+1,jj)                       &
                              &            * ahmt(ji+1,jj+1,jk)   * r1_e1e2t(ji+1,jj+1) / e3t(ji+1,jj+1,jk, Kbb)            &
                              &                                 * zu(ji,jj+1,jk)      * e2u(ji,jj+1) * e3u(ji,jj+1,jk, Kbb)
                      ELSE IF ( tmask(ji+1,jj  ,1)==0._wp ) THEN   ! East
                        !! compensation
                        pu_rhs(ji,jj+1,jk)  = pu_rhs(ji,jj+1,jk) + 1._wp * zsign * r1_e1e2u(ji,jj+1) / e3u(ji,jj+1,jk,Kmm)    &  ! +---
                             &            * e2t(ji,jj+1) *e2t(ji,jj+1) * e3t(ji,jj+1,jk,Kmm)  * r1_e2u(ji,jj+1)               &
                             &            * ahmt(ji,jj+1,jk)   * e1t(ji,jj+1) * r1_e2t(ji,jj+1) *  zv(ji,jj,jk) * r1_e1v(ji,jj) ! -Dt(i,j+1)
                        pv_rhs(ji,jj,jk)= pv_rhs(ji,jj,jk) + 1._wp * zsign * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)              & ! -+++
                             &            * e1t(ji,jj+1) * e1t(ji,jj+1) * e3t(ji,jj+1,jk,Kmm) * r1_e1v(ji,jj)                 &
                             &            * ahmt(ji,jj+1,jk) * e2t(ji,jj+1) * r1_e1t(ji,jj+1) * zu(ji,jj+1,jk) * r1_e2u(ji,jj+1)   ! +Dt(i,j+1)
                        !! divrot flux
                        pu_rhs(ji,jj+1,jk)  = pu_rhs(ji,jj+1,jk) + 1._wp * zsign     * r1_e1u(ji,jj+1)                         &
                             &            * ahmt(ji,jj+1,jk)   * r1_e1e2t(ji,jj+1) / e3t(ji,jj+1,jk, Kbb)                  &
                             &                                 * zv(ji,jj,jk)    * e1v(ji,jj) * e3v(ji,jj,jk, Kbb)
                        pv_rhs(ji,jj,jk)= pv_rhs(ji,jj,jk) + 1._wp * zsign     * r1_e2v(ji,jj)                       &
                             &            * ahmt(ji,jj+1,jk)   * r1_e1e2t(ji,jj+1) / e3t(ji,jj+1,jk, Kbb)                  &
                             &                                 * zu(ji,jj+1,jk)      * e2u(ji,jj+1) * e3u(ji,jj+1,jk, Kbb)
                      ELSE IF ( tmask(ji+1,jj+1,1)==0._wp ) THEN   ! North
                        !! compensation
                        pu_rhs(ji,jj,jk)  = pu_rhs(ji  ,jj,jk) - 1._wp * zsign * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)        &  ! +--+
                             &            * e2t(ji,jj  ) *e2t(ji,jj  ) * e3t(ji,jj  ,jk,Kmm)  * r1_e2u(ji,jj)               &
                             &            * ahmt(ji,jj,jk)   * e1t(ji,jj) * r1_e2t(ji,jj) *  zv(ji,jj,jk) * r1_e1v(ji,jj) ! -Dt(i,j)
                        pv_rhs(ji,jj,jk)= pv_rhs(ji,jj,jk) - 1._wp * zsign * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)            & ! --++
                             &            * e1t(ji,jj) * e1t(ji,jj) * e3t(ji,jj,jk,Kmm) * r1_e1v(ji,jj)                     &
                             &            * ahmt(ji,jj,jk) * e2t(ji,jj) * r1_e1t(ji,jj) * zu(ji,jj,jk) * r1_e2u(ji,jj)   ! -Dt(i,j)
                       !! divrot flux
                       pu_rhs(ji,jj,jk)  = pu_rhs(ji  ,jj,jk) - 1._wp * zsign     * r1_e1u(ji,jj)                    &
                            &            * ahmt(ji,jj,jk)   * r1_e1e2t(ji,jj) / e3t(ji,jj,jk, Kbb)                   &
                            &                                 * zv(ji,jj,jk)    * e1v(ji,jj) * e3v(ji,jj,jk, Kbb)
                       pv_rhs(ji,jj,jk)= pv_rhs(ji,jj,jk) - 1._wp * zsign     * r1_e2v(ji,jj)                    &
                            &            * ahmt(ji,jj,jk)   * r1_e1e2t(ji,jj) / e3t(ji,jj,jk, Kbb)                   &
                            &                                 * zu(ji,jj,jk)      * e2u(ji,jj) * e3u(ji,jj,jk, Kbb)
                      ELSE IF ( tmask(ji  ,jj+1,1)==0._wp ) THEN   ! West
                        !! compensation
                        pu_rhs(ji,jj,jk)  = pu_rhs(ji  ,jj,jk) + 1._wp * zsign * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)            &  ! ++-+
                             &            * e2t(ji+1,jj  ) *e2t(ji+1,jj  ) * e3t(ji+1,jj,jk,Kmm)  * r1_e2u(ji,jj)               &
                             &            * ahmt(ji+1,jj,jk)   * e1t(ji+1,jj) * r1_e2t(ji+1,jj) *  zv(ji+1,jj,jk) * r1_e1v(ji+1,jj) ! +Dt(i+1,j)
                        pv_rhs(ji+1,jj,jk)= pv_rhs(ji+1,jj,jk) + 1._wp * zsign * r1_e1e2v(ji+1,jj) / e3v(ji+1,jj,jk,Kmm)        & ! --+-
                             &            * e1t(ji+1,jj) * e1t(ji+1,jj) * e3t(ji+1,jj,jk,Kmm) * r1_e1v(ji+1,jj)                 &
                             &            * ahmt(ji+1,jj,jk) * e2t(ji+1,jj) * r1_e1t(ji+1,jj) * zu(ji,jj,jk) * r1_e2u(ji,jj)   ! -Dt(i+1,j)
                        !! divrot flux
                        pu_rhs(ji,jj,jk)  = pu_rhs(ji  ,jj,jk) + 1._wp * zsign     * r1_e1u(ji,jj)                         &
                             &            * ahmt(ji+1,jj,jk)   * r1_e1e2t(ji+1,jj) / e3t(ji+1,jj,jk, Kbb)                  &
                             &                                 * zv(ji+1,jj,jk)    * e1v(ji+1,jj) * e3v(ji+1,jj,jk, Kbb)
                        pv_rhs(ji+1,jj,jk)= pv_rhs(ji+1,jj,jk) + 1._wp * zsign     * r1_e2v(ji+1,jj)                       &
                             &            * ahmt(ji+1,jj,jk)   * r1_e1e2t(ji+1,jj) / e3t(ji+1,jj,jk, Kbb)                  &
                             &                                 * zu(ji,jj,jk)      * e2u(ji,jj) * e3u(ji,jj,jk, Kbb)
                      END IF
                    END IF
                      !
              END_2D
              !
            CASE ( np_dynldf_lap_sym_noh )
                DO_2D_00_00
                  IF ( (     tmask(ji  ,jj+1,1) + tmask(ji+1,jj+1,1)  &
                       &   + tmask(ji  ,jj  ,1) + tmask(ji+1,jj  ,1)  ) == 3._wp) THEN ! salient angle
                      !
                      IF      ( tmask(ji  ,jj  ,1)==0._wp ) THEN   ! South
                       ! ahmt(ji+1,jj+1,jk) = -1
                       pu_rhs(ji  ,jj+1,jk) = pu_rhs(ji  ,jj+1,jk) - 2._wp * zsign * r1_e1e2u(ji  ,jj+1) * ahmt(ji+1,jj+1,jk) * pv(ji+1,jj  ,jk)
                       pv_rhs(ji+1,jj  ,jk) = pv_rhs(ji+1,jj  ,jk) - 2._wp * zsign * r1_e1e2v(ji+1,jj  ) * ahmt(ji+1,jj+1,jk) * pu(ji  ,jj+1,jk)
                      ELSE IF ( tmask(ji+1,jj  ,1)==0._wp ) THEN   ! East
                       ! ahmt(ji,jj+1,jk) = -1
                       pu_rhs(ji,jj+1,jk) = pu_rhs(ji,jj+1,jk) + 2._wp * zsign * r1_e1e2u(ji,jj+1) * ahmt(ji,jj+1,jk) * pv(ji,jj  ,jk)
                       pv_rhs(ji,jj  ,jk) = pv_rhs(ji,jj  ,jk) + 2._wp * zsign * r1_e1e2v(ji,jj  ) * ahmt(ji,jj+1,jk) * pu(ji,jj+1,jk)
                      ELSE IF ( tmask(ji+1,jj+1,1)==0._wp ) THEN   ! North
                       ! ahmt(ji,jj,jk) = -1
                       pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) - 2._wp * zsign * r1_e1e2u(ji,jj) * ahmt(ji,jj,jk) * pv(ji,jj,jk)
                       pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) - 2._wp * zsign * r1_e1e2v(ji,jj) * ahmt(ji,jj,jk) * pu(ji,jj,jk)
                      ELSE IF ( tmask(ji  ,jj+1,1)==0._wp ) THEN   ! West
                        ! ahmt(ji+1,jj,jk) = -1
                        pu_rhs(ji  ,jj,jk) = pu_rhs(ji  ,jj,jk) + 2._wp * zsign * r1_e1e2u(ji,jj  ) * ahmt(ji+1,jj,jk) * pv(ji+1,jj,jk)
                        pv_rhs(ji+1,jj,jk) = pv_rhs(ji+1,jj,jk) + 2._wp * zsign * r1_e1e2v(ji,jj+1) * ahmt(ji+1,jj,jk) * pu(ji  ,jj,jk)
                      END IF
                      !
                   END IF
                END_2D
            CASE ( np_dynldf_lap_con     )
                DO_2D_00_00
                IF ( (     tmask(ji  ,jj+1,1) + tmask(ji+1,jj+1,1)  &
                     &   + tmask(ji  ,jj  ,1) + tmask(ji+1,jj  ,1)  ) == 3._wp) THEN ! salient angle
                    !
                    IF      ( tmask(ji  ,jj  ,1)==0._wp ) THEN   ! South
                      pu_rhs(ji,jj+1,jk)  = pu_rhs(ji,jj+1,jk) - 1._wp * zsign * r1_e1e2u(ji,jj+1)                          &
                           &             * ahmt(ji+1,jj+1,jk) *  pv(ji+1,jj,jk)                               ! -Dt(i+1,j+1)
                      pv_rhs(ji+1,jj,jk)= pv_rhs(ji+1,jj,jk) - 1._wp * zsign * r1_e1e2v(ji+1,jj)                          &
                           &            * ahmt(ji+1,jj+1,jk) * pu(ji,jj+1,jk)                                  ! +Dt(i+1,j+1)
                    ELSE IF ( tmask(ji+1,jj  ,1)==0._wp ) THEN   ! East
                      pu_rhs(ji,jj+1,jk)  = pu_rhs(ji,jj+1,jk) + 1._wp * zsign * r1_e1e2u(ji,jj+1)                            &
                           &            * ahmt(ji,jj+1,jk)    *  pv(ji,jj,jk)                                 ! -Dt(i,j+1)
                      pv_rhs(ji,jj,jk)= pv_rhs(ji,jj,jk) + 1._wp * zsign * r1_e1e2v(ji,jj)                                    &
                           &            * ahmt(ji,jj+1,jk)  * pu(ji,jj+1,jk)                                   ! +Dt(i,j+1)
                    ELSE IF ( tmask(ji+1,jj+1,1)==0._wp ) THEN   ! North
                      pu_rhs(ji,jj,jk)  = pu_rhs(ji  ,jj,jk) - 1._wp * zsign * r1_e1e2u(ji,jj)                            &
                           &            * ahmt(ji,jj,jk) *  pv(ji,jj,jk)                                      ! -Dt(i,j)
                      pv_rhs(ji,jj,jk)= pv_rhs(ji,jj,jk) - 1._wp * zsign * r1_e1e2v(ji,jj)                                       &
                           &            * ahmt(ji,jj,jk)  * pu(ji,jj,jk)                                   ! -Dt(i,j)
                    ELSE IF ( tmask(ji  ,jj+1,1)==0._wp ) THEN   ! West
                      pu_rhs(ji,jj,jk)  = pu_rhs(ji  ,jj,jk) + 1._wp * zsign * r1_e1e2u(ji,jj)                                &
                           &            * ahmt(ji+1,jj,jk) *  pv(ji+1,jj,jk)                                   ! +Dt(i+1,j)
                      pv_rhs(ji+1,jj,jk)= pv_rhs(ji+1,jj,jk) + 1._wp * zsign * r1_e1e2v(ji+1,jj)                              &
                           &            * ahmt(ji+1,jj,jk)  * pu(ji,jj,jk)         ! -Dt(i+1,j)
                    END IF
                    !
                  END IF
                END_2D
            CASE DEFAULT                                     ! error
               CALL ctl_stop('STOP','dyn_ldf_lap: wrong value for nn_dynldf_lap_typ'  )
            END SELECT
            !
          ELSEIF ( nn_dynldf_flx == 2 ) THEN
              SELECT CASE( nn_dynldf_lap_typ )
                 !
                 CASE ( np_dynldf_lap_sym     )
                   DO_2D_00_00
                    IF ( (     tmask(ji  ,jj+1,1) + tmask(ji+1,jj+1,1)  &
                         &   + tmask(ji  ,jj  ,1) + tmask(ji+1,jj  ,1)  ) == 3._wp) THEN ! salient angle
                        !
                        IF      ( tmask(ji  ,jj  ,1)==0._wp ) THEN   ! South
                          !! compensation + symetric flux
                          pu_rhs(ji,jj+1,jk)  = pu_rhs(ji,jj+1,jk) - 2._wp * zsign * r1_e1e2u(ji,jj+1) / e3u(ji,jj+1,jk,Kmm)  &  ! ++--
                               &            * e2t(ji+1,jj+1) *e2t(ji+1,jj+1) * e3t(ji+1,jj+1,jk,Kmm)  * r1_e2u(ji,jj+1)       &
                               &            * ahmt(ji+1,jj+1,jk)   * e1t(ji+1,jj+1) * r1_e2t(ji+1,jj+1) *  zv(ji+1,jj,jk) * r1_e1v(ji+1,jj) ! -Dt(i+1,j+1)
                          pv_rhs(ji+1,jj,jk)= pv_rhs(ji+1,jj,jk) - 2._wp * zsign * r1_e1e2v(ji+1,jj) / e3v(ji+1,jj,jk,Kmm)    & ! -++-
                               &            * e1t(ji+1,jj+1) * e1t(ji+1,jj+1) * e3t(ji+1,jj+1,jk,Kmm) * r1_e1v(ji+1,jj)       &
                               &            * ahmt(ji+1,jj+1,jk) * e2t(ji+1,jj+1) * r1_e1t(ji+1,jj+1) * zu(ji,jj+1,jk) * r1_e2u(ji,jj+1)   ! +Dt(i+1,j+1)
                        ELSE IF ( tmask(ji+1,jj  ,1)==0._wp ) THEN   ! East
                          !! compensation + symetric flux
                          pu_rhs(ji,jj+1,jk)  = pu_rhs(ji,jj+1,jk) + 2._wp * zsign * r1_e1e2u(ji,jj+1) / e3u(ji,jj+1,jk,Kmm)    &  ! +---
                               &            * e2t(ji,jj+1) *e2t(ji,jj+1) * e3t(ji,jj+1,jk,Kmm)  * r1_e2u(ji,jj+1)               &
                               &            * ahmt(ji,jj+1,jk)   * e1t(ji,jj+1) * r1_e2t(ji,jj+1) *  zv(ji,jj,jk) * r1_e1v(ji,jj) ! -Dt(i,j+1)
                          pv_rhs(ji,jj,jk)= pv_rhs(ji,jj,jk) + 2._wp * zsign * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)              & ! -+++
                               &            * e1t(ji,jj+1) * e1t(ji,jj+1) * e3t(ji,jj+1,jk,Kmm) * r1_e1v(ji,jj)                 &
                               &            * ahmt(ji,jj+1,jk) * e2t(ji,jj+1) * r1_e1t(ji,jj+1) * zu(ji,jj+1,jk) * r1_e2u(ji,jj+1)   ! +Dt(i,j+1)
                        ELSE IF ( tmask(ji+1,jj+1,1)==0._wp ) THEN   ! North
                          !! compensation + symetric flux
                          pu_rhs(ji,jj,jk)  = pu_rhs(ji  ,jj,jk) - 2._wp * zsign * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)        &  ! +--+
                               &            * e2t(ji,jj  ) *e2t(ji,jj  ) * e3t(ji,jj  ,jk,Kmm)  * r1_e2u(ji,jj)               &
                               &            * ahmt(ji,jj,jk)   * e1t(ji,jj) * r1_e2t(ji,jj) *  zv(ji,jj,jk) * r1_e1v(ji,jj) ! -Dt(i,j)
                          pv_rhs(ji,jj,jk)= pv_rhs(ji,jj,jk) - 2._wp * zsign * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)            & ! --++
                               &            * e1t(ji,jj) * e1t(ji,jj) * e3t(ji,jj,jk,Kmm) * r1_e1v(ji,jj)                     &
                               &            * ahmt(ji,jj,jk) * e2t(ji,jj) * r1_e1t(ji,jj) * zu(ji,jj,jk) * r1_e2u(ji,jj)   ! -Dt(i,j)
                        ELSE IF ( tmask(ji  ,jj+1,1)==0._wp ) THEN   ! West
                          !! compensation + symetric flux
                          pu_rhs(ji,jj,jk)  = pu_rhs(ji  ,jj,jk) + 2._wp * zsign * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)            &  ! ++-+
                               &            * e2t(ji+1,jj  ) *e2t(ji+1,jj  ) * e3t(ji+1,jj,jk,Kmm)  * r1_e2u(ji,jj)               &
                               &            * ahmt(ji+1,jj,jk)   * e1t(ji+1,jj) * r1_e2t(ji+1,jj) *  zv(ji+1,jj,jk) * r1_e1v(ji+1,jj) ! +Dt(i+1,j)
                          pv_rhs(ji+1,jj,jk)= pv_rhs(ji+1,jj,jk) + 2._wp * zsign * r1_e1e2v(ji+1,jj) / e3v(ji+1,jj,jk,Kmm)        & ! --+-
                               &            * e1t(ji+1,jj) * e1t(ji+1,jj) * e3t(ji+1,jj,jk,Kmm) * r1_e1v(ji+1,jj)                 &
                               &            * ahmt(ji+1,jj,jk) * e2t(ji+1,jj) * r1_e1t(ji+1,jj) * zu(ji,jj,jk) * r1_e2u(ji,jj)   ! -Dt(i+1,j)
                        END IF
                      END IF
                    END_2D
                CASE DEFAULT                                     ! error
                   CALL ctl_stop('STOP','dyn_ldf_lap: wrong value for nn_dynldf_lap_typ'  )
                END SELECT
            END IF
      !
   END SUBROUTINE dyn_ldf_lap


   SUBROUTINE dyn_ldf_blp( kt, Kbb, Kmm, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_ldf_blp  ***
      !!
      !! ** Purpose :   Compute the before lateral momentum viscous trend
      !!              and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The lateral viscous trends is provided by a bilaplacian
      !!      operator applied to before field (forward in time).
      !!      It is computed by two successive calls to dyn_ldf_lap routine
      !!
      !! ** Action :   pt(:,:,:,:,Krhs)   updated with the before rotated bilaplacian diffusion
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kbb, Kmm   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pu, pv     ! before velocity fields
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs   ! momentum trend
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zulap, zvlap   ! laplacian at u- and v-point
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf_blp : bilaplacian operator momentum '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      zulap(:,:,:) = 0._wp
      zvlap(:,:,:) = 0._wp
      !
      CALL dyn_ldf_lap( kt, Kbb, Kmm, pu, pv, zulap, zvlap, 1 )   ! rotated laplacian applied to pt (output in zlap,Kbb)
      !
      CALL lbc_lnk_multi( 'dynldf_lap_blp', zulap, 'U', -1., zvlap, 'V', -1. )             ! Lateral boundary conditions
      !
      CALL dyn_ldf_lap( kt, Kbb, Kmm, zulap, zvlap, pu_rhs, pv_rhs, 2 )   ! rotated laplacian applied to zlap (output in pt(:,:,:,:,Krhs))
      !
   END SUBROUTINE dyn_ldf_blp

   !!======================================================================
END MODULE dynldf_lap_blp
