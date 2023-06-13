MODULE dom_oce
   !!======================================================================
   !!                       ***  MODULE dom_oce  ***
   !!
   !! ** Purpose :   Define in memory all the ocean space domain variables
   !!======================================================================
   !! History :  1.0  ! 2005-10  (A. Beckmann, G. Madec)  reactivate s-coordinate
   !!            3.3  ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!            3.4  ! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.5  ! 2012     (S. Mocavero, I. Epicoco) Add arrays associated
   !!                             to the optimization of BDY communications
   !!            3.7  ! 2015-11  (G. Madec) introduce surface and scale factor ratio
   !!             -   ! 2015-11  (G. Madec, A. Coward)  time varying zgr by default
   !!            4.1  ! 2019-08  (A. Coward, D. Storkey) rename prognostic variables in preparation for new time scheme.
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   Agrif_Root    : dummy function used when lk_agrif=F
   !!   Agrif_CFixed  : dummy function used when lk_agrif=F
   !!   dom_oce_alloc : dynamical allocation of dom_oce arrays
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters

   IMPLICIT NONE
   PUBLIC             ! allows the acces to par_oce when dom_oce is used (exception to coding rules)

   PUBLIC dom_oce_alloc  ! Called from nemogcm.F90

   !!----------------------------------------------------------------------
   !! time & space domain namelist
   !! ----------------------------
   !                                   !!* Namelist namdom : time & space domain *
   LOGICAL , PUBLIC ::   ln_linssh      !: =T  linear free surface ==>> model level are fixed in time
   LOGICAL , PUBLIC ::   ln_meshmask    !: =T  create a mesh-mask file (mesh_mask.nc)
   REAL(wp), PUBLIC ::   rn_Dt          !: time step for the dynamics and tracer
   REAL(wp), PUBLIC ::   rn_atfp        !: asselin time filter parameter
   LOGICAL , PUBLIC ::   ln_1st_euler   !: =T start with forward time step or not (=F)
   LOGICAL , PUBLIC ::   ln_crs         !: Apply grid coarsening to dynamical model output or online passive tracers

   !! Free surface parameters
   !! =======================
   LOGICAL , PUBLIC :: ln_dynspg_exp    !: Explicit free surface flag
   LOGICAL , PUBLIC :: ln_dynspg_ts     !: Split-Explicit free surface flag

   !! Time splitting parameters
   !! =========================
   LOGICAL,  PUBLIC :: ln_bt_fw         !: Forward integration of barotropic sub-stepping
   LOGICAL,  PUBLIC :: ln_bt_av         !: Time averaging of barotropic variables
   LOGICAL,  PUBLIC :: ln_bt_auto       !: Set number of barotropic iterations automatically
   INTEGER,  PUBLIC :: nn_bt_flt        !: Filter choice
   INTEGER,  PUBLIC :: nn_e          !: Number of barotropic iterations during one baroclinic step (rn_Dt)
   REAL(wp), PUBLIC :: rn_bt_cmax       !: Maximum allowed courant number (used if ln_bt_auto=T)
   REAL(wp), PUBLIC :: rn_bt_alpha      !: Time stepping diffusion parameter


   !                                   !!! associated variables
   LOGICAL , PUBLIC ::   l_1st_euler    !: Euler 1st time-step flag (=T if ln_restart=F or ln_1st_euler=T)
   REAL(wp), PUBLIC ::   rDt, r1_Dt     !: Current model timestep and reciprocal
                                        !: rDt = 2 * rn_Dt if leapfrog and l_1st_euler = F
                                        !:     =     rn_Dt if leapfrog and l_1st_euler = T
                                        !:     =     rn_Dt if RK3

   !!----------------------------------------------------------------------
   !! space domain parameters
   !!----------------------------------------------------------------------
   INTEGER, PUBLIC ::   jperio   !: Global domain lateral boundary type (between 0 and 7)
   !                                !  = 0 closed                 ;   = 1 cyclic East-West
   !                                !  = 2 cyclic North-South     ;   = 3 North fold T-point pivot
   !                                !  = 4 cyclic East-West AND North fold T-point pivot
   !                                !  = 5 North fold F-point pivot
   !                                !  = 6 cyclic East-West AND North fold F-point pivot
   !                                !  = 7 bi-cyclic East-West AND North-South
   LOGICAL, PUBLIC ::   l_Iperio, l_Jperio   !   should we explicitely take care I/J periodicity

   !                             !: domain MPP decomposition parameters
   INTEGER             , PUBLIC ::   nimpp, njmpp     !: i- & j-indexes for mpp-subdomain left bottom
   INTEGER             , PUBLIC ::   nreci, nrecj     !: overlap region in i and j
   INTEGER             , PUBLIC ::   nproc            !: number for local processor
   INTEGER             , PUBLIC ::   narea            !: number for local area
   INTEGER             , PUBLIC ::   nbondi, nbondj   !: mark of i- and j-direction local boundaries
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondi_bdy(:)    !: mark i-direction local boundaries for BDY open boundaries
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondj_bdy(:)    !: mark j-direction local boundaries for BDY open boundaries
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondi_bdy_b(:)  !: mark i-direction of neighbours local boundaries for BDY open boundaries
   INTEGER, ALLOCATABLE, PUBLIC ::   nbondj_bdy_b(:)  !: mark j-direction of neighbours local boundaries for BDY open boundaries

   INTEGER, PUBLIC ::   npolj             !: north fold mark (0, 3 or 4)
   INTEGER, PUBLIC ::   nlci, nldi, nlei  !: i-dimensions of the local subdomain and its first and last indoor indices
   INTEGER, PUBLIC ::   nlcj, nldj, nlej  !: i-dimensions of the local subdomain and its first and last indoor indices
   INTEGER, PUBLIC ::   noea, nowe        !: index of the local neighboring processors in
   INTEGER, PUBLIC ::   noso, nono        !: east, west, south and north directions
   INTEGER, PUBLIC ::   nidom             !: ???

   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mig        !: local  ==> global domain i-index
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mjg        !: local  ==> global domain j-index
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mi0, mi1   !: global ==> local  domain i-index (mi0=1 and mi1=0 if the global index
   !                                                                !                                             is not in the local domain)
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mj0, mj1   !: global ==> local  domain j-index (mj0=1 and mj1=0 if the global index
   !                                                                !                                             is not in the local domain)
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nimppt, njmppt   !: i-, j-indexes for each processor
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ibonit, ibonjt   !: i-, j- processor neighbour existence
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nlcit , nlcjt    !: dimensions of every subdomain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nldit , nldjt    !: first, last indoor index for each i-domain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nleit , nlejt    !: first, last indoor index for each j-domain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: nfiimpp, nfipproc, nfilcit

   !!----------------------------------------------------------------------
   !! horizontal curvilinear coordinate and scale factors
   !! ---------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   glamt , glamu, glamv , glamf    !: longitude at t, u, v, f-points [degree]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   gphit , gphiu, gphiv , gphif    !: latitude  at t, u, v, f-points [degree]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::   e1t   , e2t  , r1_e1t, r1_e2t   !: t-point horizontal scale factors    [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::   e1u   , e2u  , r1_e1u, r1_e2u   !: horizontal scale factors at u-point [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::   e1v   , e2v  , r1_e1v, r1_e2v   !: horizontal scale factors at v-point [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) ::   e1f   , e2f  , r1_e1f, r1_e2f   !: horizontal scale factors at f-point [m]
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e1e2t , r1_e1e2t                !: associated metrics at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e1e2u , r1_e1e2u , e2_e1u       !: associated metrics at u-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e1e2v , r1_e1e2v , e1_e2v       !: associated metrics at v-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e1e2f , r1_e1e2f                !: associated metrics at f-point
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ff_f  , ff_t                    !: Coriolis factor at f- & t-points  [1/s]

!!an
   ! https://stackoverflow.com/questions/2893097/fortran-save-statement
   ! SAVE statement
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,  DIMENSION(:,:) ::   glamt0, gphit0  ! longitude / latitude bigger grid                     [degrees]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,  DIMENSION(:,:) ::   k_top0           ! False mask (oops)
!!an


   !!----------------------------------------------------------------------
   !! vertical coordinate and scale factors
   !! ---------------------------------------------------------------------
   LOGICAL, PUBLIC ::   ln_zco       !: z-coordinate - full step
   LOGICAL, PUBLIC ::   ln_zps       !: z-coordinate - partial step
   LOGICAL, PUBLIC ::   ln_sco       !: s-coordinate or hybrid z-s coordinate
   LOGICAL, PUBLIC ::   ln_isfcav    !: presence of ISF
   !                                                        !  reference scale factors
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::     e3t_0   !: t- vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::     e3u_0   !: u- vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::     e3v_0   !: v- vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::     e3f_0   !: f- vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::     e3w_0   !: w- vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::    e3uw_0   !: uw-vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::    e3vw_0   !: vw-vert. scale factor [m]
   !                                                        !  time-dependent scale factors
!!st1
#if ! defined key_qco
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   e3t, e3u, e3v, e3w, e3uw, e3vw  !: vert. scale factor [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   e3f                             !: F-point vert. scale factor [m]
#endif
   !                                                        !  time-dependent ratio ssh / h_0
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   r3t, r3u, r3v                   !: time-dependent    ratio at t-, u- and v-point [-]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   r3f                             !: mid-time-level    ratio at f-point            [-]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   r3t_f, r3u_f, r3v_f             !: now time-filtered ratio at t-, u- and v-point [-]

   !                                                        !  reference depths of cells
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   gdept_0  !: t- depth              [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   gdepw_0  !: w- depth              [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   gde3w_0  !: w- depth (sum of e3w) [m]
   !                                                        !  time-dependent depths of cells
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   gdept, gdepw
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   gde3w
!!st2
   !                                                        !  reference heights of ocean water column and its inverse
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   ht_0, r1_ht_0   !: t-depth        [m] and [1/m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   hu_0, r1_hu_0   !: u-depth        [m] and [1/m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   hv_0, r1_hv_0   !: v-depth        [m] and [1/m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   hf_0, r1_hf_0   !: f-depth        [m] and [1/m]
   !                                                        ! time-dependent heights of ocean water column
#if ! defined key_qco
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   ht          !: t-points           [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hu, r1_hu   !: u-depth            [m] and [1/m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hv, r1_hv   !: v-depth            [m] and [1/m]
#endif
!!st2

# if defined key_bath
 !! bathymetry field
 !! ------------
 REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   batht             !: bathymetry     [-]
 REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   bathu             !: bathymetry     [-]
 REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   bathv             !: bathymetry     [-]
#endif

   INTEGER, PUBLIC ::   nla10              !: deepest    W level Above  ~10m (nlb10 - 1)
   INTEGER, PUBLIC ::   nlb10              !: shallowest W level Bellow ~10m (nla10 + 1)

   !! 1D reference  vertical coordinate
   !! =-----------------====------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   gdept_1d, gdepw_1d !: reference depth of t- and w-points (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   e3t_1d  , e3w_1d   !: reference vertical scale factors at T- and W-pts (m)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   risfdep, bathy

   !!----------------------------------------------------------------------
   !! masks, top and bottom ocean point position
   !! ---------------------------------------------------------------------
!!gm Proposition of new name for top/bottom vertical indices
!   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mtk_t, mtk_u, mtk_v   !: top    first wet T-, U-, and V-level (ISF)
!   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mbk_t, mbk_u, mbk_v   !: bottom last  wet T-, U-, and V-level
!!gm
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mbkt, mbku, mbkv   !: bottom last wet T-, U- and V-level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tmask_i            !: interior domain T-point mask
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tmask_h            !: internal domain T-point mask (Figure 8.5 NEMO book)

   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mikt, miku, mikv, mikf  !: top first wet T-, U-, V-, F-level           (ISF)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssmask, ssumask, ssvmask, ssfmask    !: surface mask at T-,U-, V- and F-pts
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:), TARGET :: tmask, umask, vmask, fmask   !: land/ocean mask at T-, U-, V- and F-pts
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:), TARGET :: wmask, wumask, wvmask        !: land/ocean mask at WT-, WU- and WV-pts

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   tpol, fpol          !: north fold mask (jperio= 3 or 4)

   !!----------------------------------------------------------------------
   !! calendar variables
   !! ---------------------------------------------------------------------
   INTEGER , PUBLIC ::   nyear         !: current year
   INTEGER , PUBLIC ::   nmonth        !: current month
   INTEGER , PUBLIC ::   nday          !: current day of the month
   INTEGER , PUBLIC ::   nhour         !: current hour
   INTEGER , PUBLIC ::   nminute       !: current minute
   INTEGER , PUBLIC ::   ndastp        !: time step date in yyyymmdd format
   INTEGER , PUBLIC ::   nday_year     !: current day counted from jan 1st of the current year
   INTEGER , PUBLIC ::   nsec_year     !: seconds between 00h jan 1st of the current  year and half of the current time step
   INTEGER , PUBLIC ::   nsec_month    !: seconds between 00h 1st day of the current month and half of the current time step
   INTEGER , PUBLIC ::   nsec_monday   !: seconds between 00h         of the last Monday   and half of the current time step
   INTEGER , PUBLIC ::   nsec_day      !: seconds between 00h         of the current   day and half of the current time step
   REAL(wp), PUBLIC ::   fjulday       !: current julian day
   REAL(wp), PUBLIC ::   fjulstartyear !: first day of the current year in julian days
   REAL(wp), PUBLIC ::   adatrj        !: number of elapsed days since the begining of the whole simulation
   !                                   !: (cumulative duration of previous runs that may have used different time-step size)
   INTEGER , PUBLIC, DIMENSION(  0: 2) ::   nyear_len     !: length in days of the previous/current/next year
   INTEGER , PUBLIC, DIMENSION(-11:25) ::   nmonth_len    !: length in days of the months of the current year
   INTEGER , PUBLIC, DIMENSION(-11:25) ::   nmonth_beg    !: second since Jan 1st 0h of the current year and the half of the months
   INTEGER , PUBLIC                  ::   nsec1jan000     !: second since Jan 1st 0h of nit000 year and Jan 1st 0h the current year
   INTEGER , PUBLIC                  ::   nsec000_1jan000   !: second since Jan 1st 0h of nit000 year and nit000
   INTEGER , PUBLIC                  ::   nsecend_1jan000   !: second since Jan 1st 0h of nit000 year and nitend

   !!----------------------------------------------------------------------
   !! agrif domain
   !!----------------------------------------------------------------------
#if defined key_agrif
   LOGICAL, PUBLIC, PARAMETER ::   lk_agrif = .TRUE.    !: agrif flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_agrif = .FALSE.   !: agrif flag
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dom_oce.F90 12983 2020-05-27 15:46:51Z techene $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if ! defined key_agrif
   !!----------------------------------------------------------------------
   !! NOT 'key_agrif'      dummy function                     No AGRIF zoom
   !!----------------------------------------------------------------------
   LOGICAL FUNCTION Agrif_Root()
      Agrif_Root = .TRUE.
   END FUNCTION Agrif_Root

   CHARACTER(len=3) FUNCTION Agrif_CFixed()
      Agrif_CFixed = '0'
   END FUNCTION Agrif_CFixed
#endif
!!st3: dom_oce_alloc modified to ease the ifdef if necessary (gm stuff)
   INTEGER FUNCTION dom_oce_alloc()
      !!----------------------------------------------------------------------
      INTEGER                ::   ii
      INTEGER, DIMENSION(30) :: ierr
      !!----------------------------------------------------------------------
      ii = 0   ;   ierr(:) = 0
      !
      ii = ii+1
      ALLOCATE( mig(jpi), mjg(jpj), STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( mi0 (jpiglo) , mi1 (jpiglo),  mj0(jpjglo) , mj1 (jpjglo) ,     &
         &      tpol(jpiglo) , fpol(jpiglo)                              , STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( glamt(jpi,jpj) ,    glamu(jpi,jpj) ,  glamv(jpi,jpj) ,  glamf(jpi,jpj) ,     &
         &      gphit(jpi,jpj) ,    gphiu(jpi,jpj) ,  gphiv(jpi,jpj) ,  gphif(jpi,jpj) ,     &
         &       e1t (jpi,jpj) ,     e2t (jpi,jpj) , r1_e1t(jpi,jpj) , r1_e2t(jpi,jpj) ,     &
         &       e1u (jpi,jpj) ,     e2u (jpi,jpj) , r1_e1u(jpi,jpj) , r1_e2u(jpi,jpj) ,     &
         &       e1v (jpi,jpj) ,     e2v (jpi,jpj) , r1_e1v(jpi,jpj) , r1_e2v(jpi,jpj) ,     &
         &       e1f (jpi,jpj) ,     e2f (jpi,jpj) , r1_e1f(jpi,jpj) , r1_e2f(jpi,jpj) ,     &
         &      e1e2t(jpi,jpj) , r1_e1e2t(jpi,jpj)                                     ,     &
         &      e1e2u(jpi,jpj) , r1_e1e2u(jpi,jpj) , e2_e1u(jpi,jpj)                   ,     &
         &      e1e2v(jpi,jpj) , r1_e1e2v(jpi,jpj) , e1_e2v(jpi,jpj)                   ,     &
         &      e1e2f(jpi,jpj) , r1_e1e2f(jpi,jpj)                                     ,     &
         &      ff_f (jpi,jpj) ,    ff_t (jpi,jpj)                                     , STAT=ierr(ii) )
         !
         ALLOCATE(glamt0(jpi,jpj) ,     &
            &     gphit0(jpi,jpj) ,     &
            &     k_top0(jpi,jpj) , STAT=ierr(ii) )
         !
#if defined key_bath
    ALLOCATE( batht(jpi,jpj), bathu(jpi,jpj), bathv(jpi,jpj), STAT=ierr(ii) )
#endif
         !
      ii = ii+1
      ALLOCATE(  e3t_0(jpi,jpj,jpk) , e3u_0 (jpi,jpj,jpk) , e3v_0 (jpi,jpj,jpk) , e3f_0(jpi,jpj,jpk) ,      &
         &       e3w_0(jpi,jpj,jpk) , e3uw_0(jpi,jpj,jpk) , e3vw_0(jpi,jpj,jpk)                      ,  STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( gdept_0(jpi,jpj,jpk)     , gdepw_0(jpi,jpj,jpk)     , gde3w_0(jpi,jpj,jpk) ,      &
         &      gdept  (jpi,jpj,jpk,jpt) , gdepw  (jpi,jpj,jpk,jpt) , gde3w  (jpi,jpj,jpk) , STAT=ierr(ii) )
         !
!!st4
#if ! defined key_qco
      ii = ii+1
      ALLOCATE( e3t(jpi,jpj,jpk,jpt) , e3u (jpi,jpj,jpk,jpt) , e3v (jpi,jpj,jpk,jpt) , e3f(jpi,jpj,jpk) ,      &
         &      e3w(jpi,jpj,jpk,jpt) , e3uw(jpi,jpj,jpk,jpt) , e3vw(jpi,jpj,jpk,jpt)                    ,  STAT=ierr(ii) )
#endif
!!st4
         !
      ii = ii+1
      ALLOCATE( r3t  (jpi,jpj,jpt)   , r3u  (jpi,jpj,jpt)    , r3v  (jpi,jpj,jpt)    , r3f  (jpi,jpj) ,  &
         &      r3t_f(jpi,jpj)       , r3u_f(jpi,jpj)        , r3v_f(jpi,jpj)                         ,  STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( ht_0(jpi,jpj) ,    hu_0(jpi,jpj)    ,    hv_0(jpi,jpj)     , hf_0(jpi,jpj) ,       &
         &   r1_ht_0(jpi,jpj) , r1_hu_0(jpi,jpj) ,    r1_hv_0(jpi,jpj),   r1_hf_0(jpi,jpj) ,   STAT=ierr(ii)  )
         !
#if ! defined key_qco
      ii = ii+1
      ALLOCATE( ht  (jpi,jpj) ,    hu  (jpi,jpj,jpt),    hv  (jpi,jpj,jpt)                 ,       &
         &                      r1_hu  (jpi,jpj,jpt), r1_hv  (jpi,jpj,jpt)                 ,   STAT=ierr(ii)  )
#endif
         !
      ii = ii+1
      ALLOCATE( risfdep(jpi,jpj) , bathy(jpi,jpj) , STAT=ierr(ii)  )
         !
      ii = ii+1
      ALLOCATE( gdept_1d(jpk) , gdepw_1d(jpk) , e3t_1d(jpk) , e3w_1d(jpk) , STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( tmask_i(jpi,jpj) , tmask_h(jpi,jpj) ,                                           &
         &      ssmask (jpi,jpj) , ssumask(jpi,jpj) , ssvmask(jpi,jpj) , ssfmask(jpi,jpj) ,     &
         &      mbkt   (jpi,jpj) , mbku   (jpi,jpj) , mbkv   (jpi,jpj) ,                    STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( mikt(jpi,jpj), miku(jpi,jpj), mikv(jpi,jpj), mikf(jpi,jpj), STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( tmask(jpi,jpj,jpk) , umask(jpi,jpj,jpk) ,     &
         &      vmask(jpi,jpj,jpk) , fmask(jpi,jpj,jpk) , STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( wmask(jpi,jpj,jpk) , wumask(jpi,jpj,jpk), wvmask(jpi,jpj,jpk) , STAT=ierr(ii) )
      !
      dom_oce_alloc = MAXVAL(ierr)
      !
   END FUNCTION dom_oce_alloc

   !!======================================================================
END MODULE dom_oce
