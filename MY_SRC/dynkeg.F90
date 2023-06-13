MODULE dynkeg
   !!======================================================================
   !!                       ***  MODULE  dynkeg  ***
   !! Ocean dynamics:  kinetic energy gradient trend
   !!======================================================================
   !! History :  1.0  !  1987-09  (P. Andrich, M.-A. Foujols)  Original code
   !!            7.0  !  1997-05  (G. Madec)  Split dynber into dynkeg and dynhpg
   !!  NEMO      1.0  !  2002-07  (G. Madec)  F90: Free form and module
   !!            3.6  !  2015-05  (N. Ducousso, G. Madec)  add Hollingsworth scheme as an option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_keg      : update the momentum trend with the horizontal tke
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE trd_oce         ! trends: ocean variables
   USE trddyn          ! trend manager: dynamics
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE prtctl          ! Print control
   USE timing          ! Timing
   USE bdy_oce         ! ocean open boundary conditions

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_keg    ! routine called by step module
   PUBLIC   dyn_kegAD   ! routine called by step module

   INTEGER, PARAMETER, PUBLIC  ::   nkeg_C2     = 0   !: 2nd order centered scheme (standard scheme)
   INTEGER, PARAMETER, PUBLIC  ::   nkeg_HW     = 1   !: Hollingsworth et al., QJRMS, 1983
   INTEGER, PARAMETER, PUBLIC  ::   nkeg_RN     = 2   !: Renner, 1981
   !!an45
   INTEGER, PARAMETER, PUBLIC  ::   nkeg_C2_mirror = 3   !: 2nd order centered scheme (wet point only : coastline at 45 degrees)
   !
   REAL(wp) ::   r1_48 = 1._wp / 48._wp   !: =1/(2*3*8)   for nkeg_HW
   REAL(wp) ::   r1_24 = 1._wp / 24._wp   !: =1/(2*3*4)   for nkeg_RN

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynkeg.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_keg( kt, kscheme, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_keg  ***
      !!
      !! ** Purpose :   Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy and add it to the
      !!      general momentum trend.
      !!
      !! ** Method  : * kscheme = nkeg_C2 : 2nd order centered scheme that
      !!      conserve kinetic energy. Compute the now horizontal kinetic energy
      !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
      !!              * kscheme = nkeg_HW : Hollingsworth correction following
      !!      Arakawa (2001). The now horizontal kinetic energy is given by:
      !!         zhke = 1/6 [ mi-1(  2 * un^2 + ((u(j+1)+u(j-1))/2)^2  )
      !!                    + mj-1(  2 * vn^2 + ((v(i+1)+v(i-1))/2)^2  ) ]
      !!
      !!      Take its horizontal gradient and add it to the general momentum
      !!      trend.
      !!         u(rhs) = u(rhs) - 1/e1u di[ zhke ]
      !!         v(rhs) = v(rhs) - 1/e2v dj[ zhke ]
      !!
      !! ** Action : - Update the (puu(:,:,:,Krhs), pvv(:,:,:,Krhs)) with the hor. ke gradient trend
      !!             - send this trends to trd_dyn (l_trddyn=T) for post-processing
      !!
      !! ** References : Arakawa, A., International Geophysics 2001.
      !!                 Hollingsworth et al., Quart. J. Roy. Meteor. Soc., 1983.
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt               ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  kscheme          ! =0/1/2   type of KEG scheme
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs        ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv         ! ocean velocities and RHS of momentum equation
      !
      INTEGER  ::   ji, jj, jk             ! dummy loop indices
      REAL(wp) ::   zu, zv                   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk)        ::   zhke
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_keg')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_keg : kinetic energy gradient trend, scheme number=', kscheme
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      IF( l_trddyn ) THEN           ! Save the input trends
         ALLOCATE( ztrdu(jpi,jpj,jpk) , ztrdv(jpi,jpj,jpk) )
         ztrdu(:,:,:) = puu(:,:,:,Krhs)
         ztrdv(:,:,:) = pvv(:,:,:,Krhs)
      ENDIF

      zhke(:,:,jpk) = 0._wp

      SELECT CASE ( kscheme )             !== Horizontal kinetic energy at T-point  ==!
!!an45 : "wet points only" dans x2 le terme quadratic a la coast (nn_dynkeg_adv = 3)
      CASE ( nkeg_C2_mirror )                          !-- Mirrored Standard scheme  --!
         DO_3D_01_01( 1, jpkm1 )
            zu =  (   puu(ji-1,jj  ,jk,Kmm) * puu(ji-1,jj  ,jk,Kmm)   &
               &    + puu(ji  ,jj  ,jk,Kmm) * puu(ji  ,jj  ,jk,Kmm)   ) * ( 2._wp - umask(ji-1,jj,jk) * umask(ji,jj,jk) )
            zv =  (   pvv(ji  ,jj-1,jk,Kmm) * pvv(ji  ,jj-1,jk,Kmm)   &
               &    + pvv(ji  ,jj  ,jk,Kmm) * pvv(ji  ,jj  ,jk,Kmm)   ) * ( 2._wp - vmask(ji,jj-1,jk) * vmask(ji,jj,jk) )
            zhke(ji,jj,jk) = 0.25_wp * ( zv + zu )
         END_3D
!!an45
      !
      CASE ( nkeg_C2 )                          !--  Standard scheme  --!
         DO_3D_01_01( 1, jpkm1 )
            zu =    puu(ji-1,jj  ,jk,Kmm) * puu(ji-1,jj  ,jk,Kmm)   &
               &  + puu(ji  ,jj  ,jk,Kmm) * puu(ji  ,jj  ,jk,Kmm)
            zv =    pvv(ji  ,jj-1,jk,Kmm) * pvv(ji  ,jj-1,jk,Kmm)   &
               &  + pvv(ji  ,jj  ,jk,Kmm) * pvv(ji  ,jj  ,jk,Kmm)
            zhke(ji,jj,jk) = 0.25_wp * ( zv + zu )
         END_3D
      CASE ( nkeg_HW )                          !--  Hollingsworth scheme  --!
         DO_3D_00_00( 1, jpkm1 )
            zu = 8._wp * ( puu(ji-1,jj  ,jk,Kmm) * puu(ji-1,jj  ,jk,Kmm)    &
               &         + puu(ji  ,jj  ,jk,Kmm) * puu(ji  ,jj  ,jk,Kmm) )  &
               &   +     ( puu(ji-1,jj-1,jk,Kmm) + puu(ji-1,jj+1,jk,Kmm) ) * ( puu(ji-1,jj-1,jk,Kmm) + puu(ji-1,jj+1,jk,Kmm) )   &
               &   +     ( puu(ji  ,jj-1,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) ) * ( puu(ji  ,jj-1,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) )
               !
            zv = 8._wp * ( pvv(ji  ,jj-1,jk,Kmm) * pvv(ji  ,jj-1,jk,Kmm)    &
               &         + pvv(ji  ,jj  ,jk,Kmm) * pvv(ji  ,jj  ,jk,Kmm) )  &
               &  +      ( pvv(ji-1,jj-1,jk,Kmm) + pvv(ji+1,jj-1,jk,Kmm) ) * ( pvv(ji-1,jj-1,jk,Kmm) + pvv(ji+1,jj-1,jk,Kmm) )   &
               &  +      ( pvv(ji-1,jj  ,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) ) * ( pvv(ji-1,jj  ,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) )
            zhke(ji,jj,jk) = r1_48 * ( zv + zu )
         END_3D
      CASE ( nkeg_RN )                          !--  Renner 1981 scheme  --!
         DO_3D_00_00( 1, jpkm1 )
            zu =               puu(ji-1,jj+1,jk,Kmm) * puu(ji-1,jj+1,jk,Kmm) + puu(ji,jj+1,jk,Kmm) * puu(ji,jj+1,jk,Kmm)    &
               &  + 4._wp * (  puu(ji-1,jj  ,jk,Kmm) * puu(ji-1,jj  ,jk,Kmm) + puu(ji,jj  ,jk,Kmm) * puu(ji,jj  ,jk,Kmm)  ) &
               &  +            puu(ji-1,jj-1,jk,Kmm) * puu(ji-1,jj-1,jk,Kmm) + puu(ji,jj-1,jk,Kmm) * puu(ji,jj-1,jk,Kmm)
               !
            zv =               pvv(ji+1,jj-1,jk,Kmm) + pvv(ji+1,jj-1,jk,Kmm) + pvv(ji+1,jj,jk,Kmm) * pvv(ji+1,jj,jk,Kmm)    &
               &  + 4._wp * (  pvv(ji  ,jj-1,jk,Kmm) * pvv(ji  ,jj-1,jk,Kmm) + pvv(ji  ,jj,jk,Kmm) * pvv(ji  ,jj,jk,Kmm)  ) &
               &  +            pvv(ji-1,jj-1,jk,Kmm) + pvv(ji-1,jj-1,jk,Kmm) + pvv(ji-1,jj,jk,Kmm) * pvv(ji-1,jj,jk,Kmm)
            zhke(ji,jj,jk) = r1_24 * ( zv + zu )
         END_3D
         CALL lbc_lnk( 'dynkeg', zhke, 'T', 1. )
         !
      END SELECT
      !
      DO_3D_00_00( 1, jpkm1 )
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - ( zhke(ji+1,jj  ,jk) - zhke(ji,jj,jk) ) / e1u(ji,jj)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - ( zhke(ji  ,jj+1,jk) - zhke(ji,jj,jk) ) / e2v(ji,jj)
      END_3D
      !
      IF( l_trddyn ) THEN                 ! save the Kinetic Energy trends for diagnostic
         ztrdu(:,:,:) = puu(:,:,:,Krhs) - ztrdu(:,:,:)
         ztrdv(:,:,:) = pvv(:,:,:,Krhs) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_keg, kt, Kmm )
         DEALLOCATE( ztrdu , ztrdv )
      ENDIF
      !
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' keg  - Ua: ', mask1=umask,   &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_keg')
      !
   END SUBROUTINE dyn_keg


   SUBROUTINE dyn_kegAD( kt, kscheme, puu, pvv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_kegAD  ***
      !!
      !! ** Purpose :   Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy and add it to the
      !!      general momentum trend.
      !!
      !! ** Method  : * kscheme = nkeg_C2 : 2nd order centered scheme that
      !!      conserve kinetic energy. Compute the now horizontal kinetic energy
      !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
      !!              * kscheme = nkeg_HW : Hollingsworth correction following
      !!      Arakawa (2001). The now horizontal kinetic energy is given by:
      !!         zhke = 1/6 [ mi-1(  2 * un^2 + ((u(j+1)+u(j-1))/2)^2  )
      !!                    + mj-1(  2 * vn^2 + ((v(i+1)+v(i-1))/2)^2  ) ]
      !!
      !!      Take its horizontal gradient and add it to the general momentum
      !!      trend.
      !!         u(rhs) = u(rhs) - 1/e1u di[ zhke ]
      !!         v(rhs) = v(rhs) - 1/e2v dj[ zhke ]
      !!
      !! ** Action : - Update the (puu(:,:,:,Krhs), pvv(:,:,:,Krhs)) with the hor. ke gradient trend
      !!             - send this trends to trd_dyn (l_trddyn=T) for post-processing
      !!
      !! ** References : Arakawa, A., International Geophysics 2001.
      !!                 Hollingsworth et al., Quart. J. Roy. Meteor. Soc., 1983.
      !!----------------------------------------------------------------------
      !
      INTEGER                                  , INTENT( in )  ::  kt               ! ocean time-step index
      INTEGER                                  , INTENT( in )  ::  kscheme          ! =0/1/2   type of KEG scheme
      REAL(wp), DIMENSION(jpi,jpj,jpk)         , INTENT(inout) ::  puu, pvv         ! ocean velocities at Kmm
      REAL(wp), DIMENSION(jpi,jpj,jpk),OPTIONAL, INTENT(inout) ::  pu_rhs, pv_rhs   ! RHS
      !
      INTEGER  ::   ji, jj, jk             ! dummy loop indices
      REAL(wp) ::   zu, zv                   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk)        ::   zhke
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_keg')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_kegAD : kinetic energy gradient trend, scheme number=', kscheme
         IF(lwp) WRITE(numout,*) '~~~~~~~~~'
      ENDIF

      zhke(:,:,jpk) = 0._wp

      SELECT CASE ( kscheme )             !== Horizontal kinetic energy at T-point  ==!
!!an45 to be ADDED : que cas C2 - "wet points only" il suffit de x2 le terme quadratic a la coast (nn_dynkeg_adv = 2)
CASE ( nkeg_C2_mirror )                          !--  Standard scheme  --!
         DO_3D_01_01( 1, jpkm1 )
            zu =  (   puu(ji-1,jj  ,jk) * puu(ji-1,jj  ,jk)   &
               &    + puu(ji  ,jj  ,jk) * puu(ji  ,jj  ,jk)   ) * ( 2._wp - umask(ji-1,jj,jk) * umask(ji,jj,jk) )
            zv =  (   pvv(ji  ,jj-1,jk) * pvv(ji  ,jj-1,jk)   &
               &    + pvv(ji  ,jj  ,jk) * pvv(ji  ,jj  ,jk)   ) * ( 2._wp - vmask(ji,jj-1,jk) * vmask(ji,jj,jk) )
            zhke(ji,jj,jk) = 0.25_wp * ( zv + zu )
         END_3D
!!an45
      !
      CASE ( nkeg_C2 )                          !--  Standard scheme  --!
         DO_3D_01_01( 1, jpkm1 )
            zu =    puu(ji-1,jj  ,jk) * puu(ji-1,jj  ,jk)   &
               &  + puu(ji  ,jj  ,jk) * puu(ji  ,jj  ,jk)
            zv =    pvv(ji  ,jj-1,jk) * pvv(ji  ,jj-1,jk)   &
               &  + pvv(ji  ,jj  ,jk) * pvv(ji  ,jj  ,jk)
            zhke(ji,jj,jk) = 0.25_wp * ( zv + zu )
         END_3D
!!an 00_00 ?
      CASE ( nkeg_HW )                          !--  Hollingsworth scheme  --!
         DO_3D_00_00( 1, jpkm1 )
            zu = 8._wp * ( puu(ji-1,jj  ,jk) * puu(ji-1,jj  ,jk)    &
               &         + puu(ji  ,jj  ,jk) * puu(ji  ,jj  ,jk) )  &
               &   +     ( puu(ji-1,jj-1,jk) + puu(ji-1,jj+1,jk) ) * ( puu(ji-1,jj-1,jk) + puu(ji-1,jj+1,jk) )   &
               &   +     ( puu(ji  ,jj-1,jk) + puu(ji  ,jj+1,jk) ) * ( puu(ji  ,jj-1,jk) + puu(ji  ,jj+1,jk) )
               !
            zv = 8._wp * ( pvv(ji  ,jj-1,jk) * pvv(ji  ,jj-1,jk)    &
               &         + pvv(ji  ,jj  ,jk) * pvv(ji  ,jj  ,jk) )  &
               &  +      ( pvv(ji-1,jj-1,jk) + pvv(ji+1,jj-1,jk) ) * ( pvv(ji-1,jj-1,jk) + pvv(ji+1,jj-1,jk) )   &
               &  +      ( pvv(ji-1,jj  ,jk) + pvv(ji+1,jj  ,jk) ) * ( pvv(ji-1,jj  ,jk) + pvv(ji+1,jj  ,jk) )
            zhke(ji,jj,jk) = r1_48 * ( zv + zu )
         END_3D
         CALL lbc_lnk( 'dynkeg', zhke, 'T', 1. )
         !
      END SELECT
      !
         IF( PRESENT( pu_rhs ) .AND. PRESENT( pv_rhs ) ) THEN     !***  NO alternating direction  ***!
            !
            DO_3D_00_00( 1, jpkm1 )
               pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) - ( zhke(ji+1,jj  ,jk) - zhke(ji,jj,jk) ) / e1u(ji,jj)
               pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) - ( zhke(ji  ,jj+1,jk) - zhke(ji,jj,jk) ) / e2v(ji,jj)
            END_3D
            !
         ELSEIF(       PRESENT( pu_rhs ) .AND. .NOT. PRESENT( pv_rhs ) ) THEN            !***  Alternating direction : i-component  ***!
            !
            DO_3D_00_00( 1, jpkm1 )
               pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) - ( zhke(ji+1,jj  ,jk) - zhke(ji,jj,jk) ) / e1u(ji,jj)
            END_3D
            !
         ELSEIF( .NOT. PRESENT( pu_rhs ) .AND.       PRESENT( pv_rhs ) ) THEN            !***  Alternating direction : j-component  ***!
            !
            DO_3D_00_00( 1, jpkm1 )
               pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) - ( zhke(ji  ,jj+1,jk) - zhke(ji,jj,jk) ) / e2v(ji,jj)
            END_3D
            !
         ENDIF
      IF( ln_timing )   CALL timing_stop('dyn_kegAD')
      !
   END SUBROUTINE dyn_kegAD
   !!======================================================================
END MODULE dynkeg
