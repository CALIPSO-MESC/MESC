! this version is designed to work with monthly focings from ORCHIDEE model
!+==========github checked version 0.0
! tasks to do
! (6) add comments and tidy up the codes
!
! tasks done
! (1) modify "bgc_fractions" to include woody litter
! (1) write/read RESTART file
! (3) add dimension (mpft) to "mic_param_xscale"
! (5) assigm the default "micpxdef" values in "function funct"
! (2) write output file
! (4) read in a PFT-dependent parameter table
! (7) check the CUE _T depenendence (switch it off)
! (8) check the soil moisture: MIMICS used MILLENNIAL2 function; the combined model uses Yan et al. (2018)1G
! (9) microbial turnover rates: varying with NPP (MIMICS); not varying with NPP (the combined model)
! (10) use bgctype based cluster analysis to group soils "micglobal%bgctype(mbgc)"
! (11) need to "functn_frc" and "functn_c14" working with new approach
!
! this version run the model for global sensitivity analysis for
! for kinetcs=3 (the combined model)                        [defualt values] (scaling factor range)
!  1: xav:      scaling factor for V                            [1]              (0-30)               8.0e-6
!  2: xak:      scaling factor for K                            [1]              (0-30)               10.0
!  3: xfm:      scaling factor for fm                           [1]              (0.1-5.0)            0.05
!  4: xfs:      scaling factor for fs                           [1]              (0.1-5.0)            0.05
!  5: xtvmic:   scaling factor for tvmicR (0-10)                [1]              (0.1,10)             sqrt(NPPP)
!  5: xtvmic :  scaling factor for tvmicK (0-10)                [1]              =xtvmicR             sqrt(NPPP)
!  6: xtvp:     scaling factor for tvppool(0-10)                [1]            (0.1,10)               1/25  year-1    ! rate of disaggregation
!  7: xtvc:     scaling factor for tvcpool(0-10)                [1]            (0.1,10)               1/100 year-1    ! rate of MAOC breakdown
!  8: xtvac:    scaling factor for tvac   (0-10)                [1]            (0.1,10)               1/2.0 year-1    ! leaching rate
!  9: xkba:     scaling factor for kba    (0.2-5)               [1]            (0.5,10)               2.0             ! ratio of adsorption/desoprtion
! 10: xqmaxcoeff:coefficient of Qmax on clay+silt (0.4-0.8)     [1]            (0.5,5.0)              0.6
! 11: xdiffsoc:  SOC diffusion/bioturbation rate                [1]            (0.1,10.0)             (1.0/24.0)* 2.74e-3 (1/day)
! 12: xNPP:      carbon input                                   [1]            (0.5,2.0)              NPP
! 13: xrootbeta: scaling for depth-dependent of root C input    [1]            (0.5,5.0)              2.0
! 14: xvmaxbeta: scaling for depth-dependent of vmax            [1]            (0.5,5.0)              2.0
! 15: xfp2ax:    scaling factor for fp2ax                       [1]            (0.5-2.0)                              ! not used
! 16: xdesorp:   desorption coefficient (kinetics=1,2)          [1]            (0.1,10.0)                             ! not used
! 17: xbeta:     beta parameter                                 [1]            (0.55,1.0)                             ! fixed to 2
module mic_constant
  IMPLICIT NONE

  ! All module members are public by default
  public

  integer,  parameter  :: r_2 = SELECTED_REAL_KIND(12, 60)
  integer,  parameter  :: diag=0       ! =1 for printout 0 no prinout
  integer,  parameter  :: outp=1       ! output site
  !integer,  parameter  :: msite=213   ! number of sites
  integer                 :: mp           ! number of site the model runs for
  integer                 :: ntime        ! 365  !12 * 4 ! 4 year's monthly global forcings
  integer                 :: mpft         ! =17 !15      ! number of PFTs =17 FOR cable AND =19 FOR orchidee
  integer                 :: mbgc         ! number of soil categories
  integer,  parameter  :: ms= 7        !15       ! soil layers
  real(r_2) :: zse(ms)
  data zse/0.2,0.2,0.2,0.2,0.2,0.5,0.5/
  integer,  parameter  :: nlon =180
  integer,  parameter  :: nlat =90
  integer,  parameter  :: mcpool=10    ! number of C pools
  integer,  parameter  :: nfvar=22     ! number of data input variables
  real(r_2),parameter  :: delt= 1.0    ! one hour
  real(r_2),PARAMETER  :: tvc14 = (1.0/(24.0*365.0))* log(2.0)/5730.0    ! 1/hour
  integer,  parameter  :: nyic14=1940  ! year 0 of 14C record
  integer,  parameter  :: nyec14=2020  ! last yr of 14C calculation
  real(r_2),parameter  :: thresh_patchfrac=1.0e-6   ! minimial patch area fraction
!  real(r_2),PARAMETER  :: diffsoc  =(1.0/24.0)* 2.74e-3  !cm2/hour
!                                       ! m2/hour  ! see Table 1,Camino-Serrano et al. (2018)

end module mic_constant

module mic_variable
  use mic_constant
  IMPLICIT NONE

  ! All module members are public by default
  public

  SAVE

  TYPE mic_param_xscale
    ! parameter scaling
     real(r_2),dimension(:), allocatable  :: &
      xav,         &
      xak,         &
      xfp2ax,      &
      xfm,         &
      xfs,         &
      xtvmic,      &
      xtvp,        &
      xtvc,        &
      xtvac,       &
      xkba,        &
      xqmaxcoeff,  &
      xbeta,       &
      xdiffsoc,    &
      xnpp,        &
      xdesorp,     &
      xrootbeta,   &
      xvmaxbeta
  END TYPE mic_param_xscale

  TYPE mic_param_default
     !default values for Michaelis-Menten K
     real(r_2)  ::   &
     sk =0.017,      &
     skx=0.027,      &
     ak = 10.0,      &
     bk = 3.19,      &
     xk1 =8.0,       &
     xk2 =2.0,       &
     xk3 =4.0,       &
     xj1 =2.0,       &
     xj2 =4.0,       &
     xj3 =6.0
     !default values for Michaelis-Menten Vmax
     real(r_2)  ::   &
     sv = 0.063,     &
     av = 10.0*8.0e-6,    &
     bv = 5.47,      &
     xv1= 10.0,      &
     xv2= 2.0,       &
     xv3= 10.0,      &
     xw1= 3.0,       &
     xw2= 3.0,       &
     xw3= 2.0
     ! default values of MM kinetics (Wieder et al. 2015)
     real(r_2)  ::   &
     Q1=4.0,         &
     Q2=4.0,         &
     fm=0.5,         &
     fs=0.5
     ! microbial turnover rate parameter values (Wieder et al. (2015))
     real(r_2)  ::        &
     xtv      = 100.0,    &
     betamic  = 2.0,      &
     tvmicR   = 0.00052,  &
     tvmicK   = 0.00024
     !dependence on the partitioning of necromass on soil clay and substrate quality
     real(r_2)  ::   &
     fmicsom1=0.432, &
     fmicsom2=0.098, &
     fmicsom3=10.56, &
     fmicsom4=29.78, &
     fmicsom5=2.61
     !microbial carbon use efficiency
     real(r_2)  ::          &
     cuemax    = 0.80,      &
     cue_coef1 = 0.66,      &
     cue_coef2 = 1.23,      &
     epislon1 = 0.5,        &
     epislon2 = 0.25,       &
     epislon3 = 0.7,        &
     epislon4 = 0.35
     !adsorption dependence on soil pH (!Table A1 Abramoff et al. (2022))
     real(r_2)  ::          &
     phcoeff1 = 0.2429,      &
     phcoeff2 = -0.3632
!     phcoeff1 = 0.186,      &
!     phcoeff2 = 0.216
     !dependence on soil moisture  Yan et al (2018)
     real(r_2)  ::          &
     smkdesorp = 0.1,       &
     smexpns   = 2.0,       &
     smexpb    = 0.75
     ! dependence of Qmax on soil texture
     real(r_2) :: qmaxcoeff = 0.4    ! Georgiou et al. (2022)
     ! SOC diffusion coefficient. see Table 1,Camino-Serrano et al. (2018)
     real(r_2):: diffsoc  =(1.0/24.0)* 2.74e-3  !cm2/hour
     ! kinetic parameter values of rkinetics=3. Table A1  Abramoff2022
     real(r_2) ::           &
     kadsorpx = 0.001,      &
     kbax     = 6.0,        &
     fp2ax    = 1.143 * 0.33,       &
     tvcpoolx = 0.02*(1.0/1.60)*1.0/(365.0*24.0), &
     tvppoolx = 0.10*(1.0/0.44)*1.0/(365.0*24.0), &
     tvacx = 0.0015/24.0,        &
     rootbeta = 6.1,          &
     vmaxbeta = 0.5

   ! previous values
   !   rootbeta = 2.0,          &
   !   tvcpoolx = 0.102 * 0.02/24.0/2.0,  &
   !   tvppoolx = 4.705 * 0.019/24./10.00, &
   !   tvacx    = 0.1   * 0.0015/24.0,&
  END TYPE mic_param_default

  TYPE mic_parameter
  real(r_2), dimension(:,:),    allocatable  :: K1,K2,K3,J1,J2,J3
  real(r_2), dimension(:,:),    allocatable  :: V1,V2,V3,W1,W2,W3
  real(r_2), dimension(:,:),    allocatable  :: desorp
  real(r_2), dimension(:,:),    allocatable  :: Q1,Q2,fm,fs
  real(r_2), dimension(:,:),    allocatable  :: mgeR1,mgeR2,mgeR3,mgeK1,mgeK2,mgeK3
  real(r_2), dimension(:,:),    allocatable  :: tvmicR,tvmicK,betamicR,betamicK
  real(r_2), dimension(:,:),    allocatable  :: fmetave
  real(r_2), dimension(:,:,:),  allocatable  :: cn_r
  real(r_2), dimension(:,:),    allocatable  :: fr2p,fk2p,fr2c,fk2c,fr2a,fk2a
  real(r_2), dimension(:),      allocatable  :: xcnleaf,xcnroot,xcnwood,fligleaf,fligroot,fligwood
  real(r_2), dimension(:),      allocatable  :: diffsocx
  ! additional parameters for kinetics3
  real(r_2), dimension(:,:),    allocatable  :: kdesorp   !mg C cm-3 hour-1
  real(r_2), dimension(:,:),    allocatable  :: kadsorp   !1/hour
  real(r_2), dimension(:,:),    allocatable  :: fp2a
  real(r_2), dimension(:,:),    allocatable  :: tvcpool   !1/hour
  real(r_2), dimension(:,:),    allocatable  :: tvppool   !1/hour
  real(r_2), dimension(:,:),    allocatable  :: tvac      !1/hour (leaching rate coefficient)
  real(r_2), dimension(:,:),    allocatable  :: qmaxcoeff !coefficient relating qmax to soil clay+silt

  ! the following are alrealy available in CABLE
  integer,   dimension(:),      allocatable  :: pft,bgctype,isoil,sorder,region,siteid,dataid
  real(r_2), dimension(:,:),    allocatable  :: sdepth,fracroot
  real(r_2), dimension(:,:),    allocatable  :: csoilobs,csoilobsp,csoilobsm
  real(r_2), dimension(:),      allocatable  :: c14soilobsp,c14soilobsm     !14C obs
  real(r_2), dimension(:,:,:),  allocatable  :: c14atm         !atmospheric 14C
  integer,   dimension(:),      allocatable  :: nyc14obs       !year at which 14C was observed
  integer,   dimension(:),      allocatable  :: top,bot

  END TYPE mic_parameter

  TYPE mic_input
  real(r_2), dimension(:,:),    allocatable  :: tavg,wavg,tair,ph,clay,silt,porosity,bulkd,matpot
  real(r_2), dimension(:),      allocatable  :: dleaf,dwood,droot
  real(r_2), dimension(:,:),    allocatable  :: cinputm
  real(r_2), dimension(:,:),    allocatable  :: cinputs
  real(r_2), dimensioN(:),      allocatable  :: fcnpp

  END TYPE mic_input

  TYPE mic_global_input
    real(r_2), dimension(:),      allocatable  :: lon,lat,time
    integer,   dimension(:),      allocatable  :: pft,bgctype,isoil,sorder,siteid
    real(r_2), dimension(:),      allocatable  :: area,npp,ph,clay,silt,poros,bulkd,avgts,avgms
    real(r_2), dimension(:,:,:),  allocatable  :: tsoil,moist,matpot
    real(r_2), dimension(:),      allocatable  :: ligleaf,ligwood,ligroot
    real(r_2), dimension(:,:),    allocatable  :: dleaf,dwood,droot,cnleaf,cnwood,cnroot
  END TYPE mic_global_input

  TYPE mic_output
  real(r_2), dimension(:),    allocatable  :: fluxcinput
  real(r_2), dimension(:),    allocatable  :: fluxrsoil
  real(r_2), dimension(:),    allocatable  :: fluxcleach
  END TYPE mic_output

  TYPE mic_cpool
  real(r_2), dimension(:,:,:),  allocatable  :: cpool
  real(r_2), dimension(:,:,:),  allocatable  :: cpooleq
  real(r_2), dimension(:),      allocatable  :: cpooleqp,cpooleqm,c12pooleqp,c12pooleqm
  END TYPE mic_cpool

  TYPE mic_npool
  real(r_2), dimension(:,:),    allocatable  :: mineralN
  END TYPE mic_npool


 CONTAINS

  SUBROUTINE mic_allocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
   IMPLICIT NONE
   TYPE(mic_parameter),    INTENT(INOUT)  :: micparam
   TYPE(mic_param_xscale), INTENT(INOUT)  :: micpxdef
   integer  :: mpft,mbgc,mp,ms

    allocate(micpxdef%xav(mbgc),  &
      micpxdef%xak(mbgc),         &
      micpxdef%xfp2ax(mbgc),      &
      micpxdef%xfm(mbgc),         &
      micpxdef%xfs(mbgc),         &
      micpxdef%xtvmic(mbgc),      &
      micpxdef%xtvp(mbgc),        &
      micpxdef%xtvc(mbgc),        &
      micpxdef%xtvac(mbgc),       &
      micpxdef%xkba(mbgc),        &
      micpxdef%xqmaxcoeff(mbgc),  &
      micpxdef%xbeta(mbgc),       &
      micpxdef%xdiffsoc(mbgc),    &
      micpxdef%xnpp(mpft),        &
      micpxdef%xdesorp(mbgc),     &
      micpxdef%xrootbeta(mpft),   &
      micpxdef%xvmaxbeta(mbgc))

    allocate(micparam%K1(mp,ms),  &
             micparam%K2(mp,ms),  &
             micparam%K3(mp,ms),  &
             micparam%J1(mp,ms),  &
             micparam%J2(mp,ms),  &
             micparam%J3(mp,ms),  &
             micparam%V1(mp,ms),  &
             micparam%V2(mp,ms),  &
             micparam%V3(mp,ms),  &
             micparam%W1(mp,ms),  &
             micparam%W2(mp,ms),  &
             micparam%W3(mp,ms),  &
             micparam%desorp(mp,ms),  &
             micparam%Q1(mp,ms),      &
             micparam%Q2(mp,ms),      &
             micparam%fm(mp,ms),      &
             micparam%fs(mp,ms),      &
             micparam%mgeR1(mp,ms),   &
             micparam%mgeR2(mp,ms),   &
             micparam%mgeR3(mp,ms),   &
             micparam%mgeK1(mp,ms),   &
             micparam%mgeK2(mp,ms),   &
             micparam%mgeK3(mp,ms),   &
             micparam%fmetave(mp,ms), &
             micparam%tvmicR(mp,ms),  &
             micparam%tvmicK(mp,ms),  &
             micparam%betamicR(mp,ms),     &
             micparam%betamicK(mp,ms),     &
             micparam%cn_r(mp,ms,mcpool),  &
             micparam%fr2p(mp,ms),   &
             micparam%fk2p(mp,ms),   &
             micparam%fr2c(mp,ms),   &
             micparam%fk2c(mp,ms),   &
             micparam%fr2a(mp,ms),   &
             micparam%fk2a(mp,ms))

    allocate(micparam%xcnleaf(mp),   &
             micparam%xcnroot(mp),   &
             micparam%xcnwood(mp),   &
             micparam%fligleaf(mp),  &
             micparam%fligroot(mp),  &
             micparam%fligwood(mp),  &
             micparam%diffsocx(mp))

    allocate(micparam%pft(mp),       &
             micparam%bgctype(mp),   &
             micparam%isoil(mp),     &
             micparam%sorder(mp),    &
             micparam%region(mp),    &
             micparam%siteid(mp),    &
             micparam%dataid(mp))

    allocate(micparam%sdepth(mp,ms),   &
             micparam%fracroot(mp,ms), &
             micparam%csoilobs(mp,ms), &
             micparam%csoilobsp(mp,ms), &
             micparam%csoilobsm(mp,ms), &
             micparam%c14soilobsp(mp), &
             micparam%c14soilobsm(mp), &
             micparam%c14atm(79,5,2),   &
             micparam%nyc14obs(mp),    &
             micparam%top(mp),        &
             micparam%bot(mp))

! additional variables for kinetics3
    allocate(micparam%kdesorp(mp,ms), &
             micparam%kadsorp(mp,ms), &
             micparam%fp2a(mp,ms),    &
             micparam%tvcpool(mp,ms), &
             micparam%tvppool(mp,ms), &
             micparam%tvac(mp,ms),    &
             micparam%qmaxcoeff(mp,ms))
  END SUBROUTINE mic_allocate_parameter

  SUBROUTINE mic_allocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
   IMPLICIT NONE
   integer :: mp,ms,nlon,nlat,ntime
   TYPE(mic_input),        INTENT(INOUT)  :: micinput
   TYPE(mic_global_input), INTENT(INOUT)  :: micglobal

    allocate(micinput%tavg(mp,ms),    &
             micinput%wavg(mp,ms),    &
             micinput%ph(mp,ms),      &
             micinput%clay(mp,ms),    &
             micinput%silt(mp,ms),    &
             micinput%bulkd(mp,ms),   &
             micinput%porosity(mp,ms),&
             micinput%matpot(mp,ms),  &
             micinput%tair(mp,365),   &
             micinput%fcnpp(mp),      &
             micinput%dleaf(mp),      &
             micinput%dwood(mp),      &
             micinput%droot(mp),      &
             micinput%cinputm(mp,ms), &
             micinput%cinputs(mp,ms) )

    allocate(micglobal%lon(mp),             &
             micglobal%lat(mp),             &
             micglobal%time(ntime),         &
             micglobal%pft(mp),             &
             micglobal%bgctype(mp),         &
             micglobal%isoil(mp),           &
             micglobal%sorder(mp),          &
             micglobal%siteid(mp),          &
             micglobal%area(mp),            &
             micglobal%npp(mp),             &
             micglobal%ph(mp),              &
             micglobal%clay(mp),            &
             micglobal%silt(mp),            &
             micglobal%poros(mp),           &
             micglobal%bulkd(mp),           &
             micglobal%avgts(mp),           &
             micglobal%avgms(mp),           &
             micglobal%tsoil(mp,ms,ntime),  &
             micglobal%moist(mp,ms,ntime),  &
             micglobal%matpot(mp,ms,ntime), &
             micglobal%ligleaf(mp),         &
             micglobal%ligwood(mp),         &
             micglobal%ligroot(mp),         &
             micglobal%dleaf(mp,ntime),     &
             micglobal%dwood(mp,ntime),     &
             micglobal%droot(mp,ntime),     &
             micglobal%cnleaf(mp,ntime),    &
             micglobal%cnwood(mp,ntime),    &
             micglobal%cnroot(mp,ntime))
  END SUBROUTINE mic_allocate_input

  SUBROUTINE mic_allocate_output(mp,micoutput)
   IMPLICIT NONE
   TYPE(mic_output), INTENT(INOUT)  :: micoutput
   integer  :: mp

   allocate(micoutput%fluxcinput(mp))
   allocate(micoutput%fluxrsoil(mp))
   allocate(micoutput%fluxcleach(mp))

 END SUBROUTINE mic_allocate_output

 SUBROUTINE mic_allocate_cpool(mp,ms,miccpool)
   IMPLICIT NONE
   integer :: mp,ms
   TYPE(mic_cpool), INTENT(INOUT)  :: miccpool
   allocate(miccpool%cpool(mp,ms,mcpool), &
            miccpool%cpooleq(mp,ms,mcpool), &
            miccpool%cpooleqp(mp),   &
            miccpool%cpooleqm(mp),  &
            miccpool%c12pooleqp(mp), &
            miccpool%c12pooleqm(mp))
 END SUBROUTINE mic_allocate_cpool


  SUBROUTINE mic_allocate_npool(mp,ms,micnpool)
   IMPLICIT NONE
   integer :: mp,ms
   TYPE(mic_npool), INTENT(INOUT)  :: micnpool

   ALLOCATE(micnpool%mineralN(mp,ms))

  END SUBROUTINE mic_allocate_npool

  ! deallocate to free up storage

  SUBROUTINE mic_deallocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
   IMPLICIT NONE
   TYPE(mic_parameter), INTENT(INOUT)     :: micparam
   TYPE(mic_param_xscale), INTENT(INOUT)  :: micpxdef
   integer  :: mpft,mbgc,mp,ms

    deallocate(micpxdef%xav,  &
      micpxdef%xak,         &
      micpxdef%xfp2ax,      &
      micpxdef%xfm,         &
      micpxdef%xfs,         &
      micpxdef%xtvmic,      &
      micpxdef%xtvp,        &
      micpxdef%xtvc,        &
      micpxdef%xtvac,       &
      micpxdef%xkba,        &
      micpxdef%xqmaxcoeff,  &
      micpxdef%xbeta,       &
      micpxdef%xdiffsoc,    &
      micpxdef%xnpp,        &
      micpxdef%xdesorp,     &
      micpxdef%xrootbeta,   &
      micpxdef%xvmaxbeta)


    deallocate(micparam%K1,  &
             micparam%K2,  &
             micparam%K3,  &
             micparam%J1,  &
             micparam%J2,  &
             micparam%J3,  &
             micparam%V1,  &
             micparam%V2,  &
             micparam%V3,  &
             micparam%W1,  &
             micparam%W2,  &
             micparam%W3,  &
             micparam%desorp,  &
             micparam%Q1,      &
             micparam%Q2,      &
             micparam%fm,      &
             micparam%fs,      &
             micparam%mgeR1,   &
             micparam%mgeR2,   &
             micparam%mgeR3,   &
             micparam%mgeK1,   &
             micparam%mgeK2,   &
             micparam%mgeK3,   &
             micparam%fmetave, &
             micparam%tvmicR,  &
             micparam%tvmicK,  &
             micparam%betamicR,     &
             micparam%betamicK,     &
             micparam%cn_r,   &
             micparam%fr2p,   &
             micparam%fk2p,   &
             micparam%fr2c,   &
             micparam%fk2c,   &
             micparam%fr2a,   &
             micparam%fk2a)

    deallocate(micparam%xcnleaf,   &
             micparam%xcnroot,   &
             micparam%xcnwood,   &
             micparam%fligleaf,  &
             micparam%fligroot,  &
             micparam%fligwood,  &
             micparam%diffsocx)

    deallocate(micparam%pft,     &
             micparam%bgctype,     &
             micparam%isoil,     &
             micparam%sorder,    &
             micparam%region,    &
             micparam%siteid)

    deallocate(micparam%sdepth,   &
             micparam%fracroot,   &
             micparam%csoilobs,   &
             micparam%csoilobsp,  &
             micparam%csoilobsm,  &
             micparam%c14soilobsp,&
             micparam%c14soilobsm,&
             micparam%c14atm,     &
             micparam%nyc14obs,   &
             micparam%top,        &
             micparam%bot)

! additional variables for kinetics3
    deallocate(micparam%kdesorp,  &
             micparam%kadsorp,  &
             micparam%fp2a,     &
             micparam%tvcpool,  &
             micparam%tvppool,  &
             micparam%tvac,     &
             micparam%qmaxcoeff)
  END SUBROUTINE mic_deallocate_parameter

  SUBROUTINE mic_deallocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
   IMPLICIT NONE
   integer :: mp,ms,nlon,nlat,ntime
   TYPE(mic_input), INTENT(INOUT)        :: micinput
   TYPE(mic_global_input), INTENT(INOUT) :: micglobal

    deallocate(micinput%tavg,    &
               micinput%wavg,    &
               micinput%ph,      &
               micinput%clay,    &
               micinput%silt,    &
               micinput%bulkd,   &
               micinput%porosity,&
               micinput%tair,    &
               micinput%fcnpp,   &
               micinput%dleaf,   &
               micinput%dwood,   &
               micinput%droot,   &
               micinput%cinputm, &
               micinput%cinputs)


    deallocate(micglobal%lon,     &
               micglobal%lat,     &
               micglobal%time,    &
               micglobal%pft,     &
               micglobal%bgctype, &
               micglobal%isoil,   &
               micglobal%sorder,  &
               micglobal%siteid,  &
               micglobal%area,    &
               micglobal%npp,     &
               micglobal%ph,      &
               micglobal%clay,    &
               micglobal%silt,    &
               micglobal%poros,   &
               micglobal%bulkd,   &
               micglobal%avgts,   &
               micglobal%avgms,   &
               micglobal%tsoil,   &
               micglobal%moist,   &
               micglobal%matpot,  &
               micglobal%ligleaf, &
               micglobal%ligwood, &
               micglobal%ligroot, &
               micglobal%dleaf,   &
               micglobal%dwood,   &
               micglobal%droot,   &
               micglobal%cnleaf,  &
               micglobal%cnwood,  &
               micglobal%cnroot)

  END SUBROUTINE mic_deallocate_input

  SUBROUTINE mic_deallocate_output(mp,micoutput)
   IMPLICIT NONE
   TYPE(mic_output), INTENT(INOUT)  :: micoutput
   integer  :: mp
    deallocate(micoutput%fluxcinput)
    deallocate(micoutput%fluxrsoil)
    deallocate(micoutput%fluxcleach)

 END SUBROUTINE mic_deallocate_output

 SUBROUTINE mic_deallocate_cpool(mp,ms,miccpool)
   IMPLICIT NONE
   integer :: mp,ms
   TYPE(mic_cpool), INTENT(INOUT)  :: miccpool
   deallocate(miccpool%cpool,  &
              miccpool%cpooleq, &
              miccpool%cpooleqp, &
              miccpool%cpooleqm, &
              miccpool%c12pooleqp,&
              miccpool%c12pooleqm)

 END SUBROUTINE mic_deallocate_cpool


  SUBROUTINE mic_deallocate_npool(mp,ms,micnpool)
   IMPLICIT NONE
   integer :: mp,ms
   TYPE(mic_npool), INTENT(INOUT)  :: micnpool

   DEALLOCATE(micnpool%mineralN)

  END SUBROUTINE mic_deallocate_npool

end module mic_variable

 real*8 function functn_c14(nx,xparam16)
   use mic_constant
   use mic_variable
   implicit none
    TYPE(mic_param_xscale)    :: micpxdef
    TYPE(mic_param_default)   :: micpdef
    TYPE(mic_parameter)       :: micparam
    TYPE(mic_input)           :: micinput
    TYPE(mic_global_input)    :: micglobal
    TYPE(mic_cpool)           :: miccpool
    TYPE(mic_npool)           :: micnpool
    TYPE(mic_output)          :: micoutput

    !local variables
    real*8,    dimension(16)           :: xparam16
    integer    :: nx
    integer,   dimension(16)           :: nxopt
    real*8,    dimension(16)           :: xopt
    real*8     :: totcost1,totcost2
    integer    :: ifsoc14,kinetics,bgcopt,jopt,nyeqpool,isoc14,jglobal,jmodel
    integer :: jrestart,nparam
    character(len=99) :: frestart_in,frestart_out,foutput
    character(len=120) :: frac14c,f14c(5)

      jrestart=0;xopt(:)=1.0
      do nparam=1,16
         nxopt(nparam) = nparam
      enddo

      frestart_in='miccpool_in.nc'
      frestart_out='miccpool_out.nc'
      foutput='vmic_output.nc'

      open(1,file='params1.txt')
!      open(91,file='modobs.txt')
!      open(92,file='modobs2.txt')
!      open(93,file='modobs_c14.txt')
!      open(94,file='modobs2_c14.txt')
      read(1,*)
      read(1,*) jglobal,ifsoc14,kinetics,bgcopt,jopt,jrestart
      read(1,11) frac14c
      read(1,11) f14c(1)
      read(1,11) f14c(2)
      read(1,11) f14c(3)
      read(1,11) f14c(4)
      read(1,11) f14c(5)
11    format(a120)
      read(1,*) xopt(1:14)

      if(jopt==0) then
         read(1,*) nxopt(1:nx)
         do nparam=1,nx
            xopt(nxopt(nparam)) = xparam16(nparam)
         enddo
      endif
      print*, xopt

      mp = 213


      totcost1 = 0.0; totcost2=0.0
      nyeqpool= 500;jmodel=1;mpft=17;mbgc=12

      call mic_allocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_allocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_allocate_output(mp,micoutput)
      call mic_allocate_cpool(mp,ms,miccpool)
      call mic_allocate_npool(mp,ms,micnpool)

          isoc14 = 0
          print *, 'isoc14 =',isoc14,'--getdata_c14'
          call getdata_c14(frac14c,f14c,micinput,micparam,micnpool)
          call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)
          print *, 'vmicsoil_c14'
          call vmicsoil_c14(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,ifsoc14,bgcopt,nyeqpool, &
                        micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

          print *, 'calcost_c14'
          call calcost_c14(nx,isoc14,bgcopt,xopt,micparam,miccpool,micinput,totcost1)

          miccpool%c12pooleqp(:) = miccpool%cpooleqp(:)
          miccpool%c12pooleqm(:) = miccpool%cpooleqm(:)

          isoc14 = 1
          print *, 'isoc14 =',isoc14,'--getdata_c14'
          call getdata_c14(frac14c,f14c,micinput,micparam,micnpool)
          call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)
          print *, 'vmicsoil_c14'
          call vmicsoil_c14(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,ifsoc14,bgcopt,nyeqpool+2000, &
                        micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

          print *, 'calcost_c14'
          call calcost_c14(nx,isoc14,bgcopt,xopt,micparam,miccpool,micinput,totcost2)
          functn_c14 = totcost1+totcost2
          print *,'tot1 = ',totcost1
          print *,'tot2 = ',totcost2

      close(1)
!      close(91)
!      close(92)
!      close(93)
!      close(94)

!      functn = totcost

      call mic_deallocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_deallocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_deallocate_output(mp,micoutput)
      call mic_deallocate_cpool(mp,ms,miccpool)
      call mic_deallocate_npool(mp,ms,micnpool)

END function functn_c14

real*8 function functn_frc(nx,xparam16)
   use mic_constant
   use mic_variable
   implicit none
    TYPE(mic_param_xscale)    :: micpxdef
    TYPE(mic_param_default)   :: micpdef
    TYPE(mic_parameter)       :: micparam
    TYPE(mic_input)           :: micinput
    TYPE(mic_global_input)    :: micglobal
    TYPE(mic_cpool)           :: miccpool
    TYPE(mic_npool)           :: micnpool
    TYPE(mic_output)          :: micoutput

    !local variables
    integer    :: nx
    integer,   dimension(16)           :: nxopt
    real*8,    dimension(16)           :: xparam16
    real*8,    dimension(16)           :: xopt
    real*8     :: totcost1
    integer    :: ifsoc14,kinetics,bgcopt,jopt,nyeqpool,isoc14,jglobal,jmodel
    integer :: jrestart,nparam
    character(len=99) :: frestart_in,frestart_out,foutput
    character(len=120) :: Cfraction


      jrestart=0;xopt(:)=1.0
      do nparam=1,16
         nxopt(nparam) = nparam
      enddo

      frestart_in='miccpool_in.nc'
      frestart_out='miccpool_out.nc'
      foutput='vmic_output.nc'

      open(1,file='params1.txt')
!      open(91,file='modobs.txt')
!      open(92,file='modobs2.txt')

      read(1,*)
      read(1,*) jglobal,ifsoc14,kinetics,bgcopt,jopt,jrestart
      read(1,11) Cfraction
11    format(a120)
      read(1,*) xopt(1:14)

      if(jopt==0) then
         read(1,*) nxopt(1:nx)
         do nparam=1,nx
            xopt(nxopt(nparam)) = xparam16(nparam)
         enddo
      endif
      print*, xopt


      !mp = 2210
      mp = 2206

      totcost1 = 0.0
      nyeqpool= 1000
      isoc14 = 0
      jmodel=1;mpft=17;mbgc=12

      call mic_allocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_allocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_allocate_output(mp,micoutput)
      call mic_allocate_cpool(mp,ms,miccpool)
      call mic_allocate_npool(mp,ms,micnpool)


          print *, 'isoc14 =',isoc14,'--getdata_frc'
          call getdata_frc(Cfraction,micinput,micparam,micnpool)
          call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)
          print *, 'vmicsoil_frc'

          call vmicsoil_c14(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,ifsoc14,bgcopt,nyeqpool, &
                        micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

          print *, 'calcost_frc'
          call calcost_frc(nx,bgcopt,xopt,micpxdef,micparam,miccpool,micinput,totcost1)

      close(1)
!      close(91)
!      close(92)

     functn_frc = totcost1
!      functn     = totcost1

      call mic_deallocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_deallocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_deallocate_output(mp,micoutput)
      call mic_deallocate_cpool(mp,ms,miccpool)
      call mic_deallocate_npool(mp,ms,micnpool)

END function functn_frc

real*8 function functn_soc_wosis(nx,xparam16)
   use mic_constant
   use mic_variable
   implicit none
    TYPE(mic_param_xscale)    :: micpxdef
    TYPE(mic_param_default)   :: micpdef
    TYPE(mic_parameter)       :: micparam
    TYPE(mic_input)           :: micinput
    TYPE(mic_global_input)    :: micglobal
    TYPE(mic_cpool)           :: miccpool
    TYPE(mic_npool)           :: micnpool
    TYPE(mic_output)          :: micoutput

    !local variables
    integer    :: nx
    integer,   dimension(16)           :: nxopt
    real*8,    dimension(16)           :: xparam16
    real*8,    dimension(16)           :: xopt
    real*8     :: totcost1
    integer    :: ifsoc14,kinetics,bgcopt,jopt,nyeqpool,isoc14,jglobal,jmodel
    integer :: jrestart,nparam
    character(len=99) :: frestart_in,frestart_out,foutput
    character(len=120) :: finputsoc


      jrestart=0;xopt(:)=1.0
      do nparam=1,16
         nxopt(nparam) = nparam
      enddo

      frestart_in='miccpool_in.nc'
      frestart_out='miccpool_out.nc'
      foutput='vmic_output.nc'

      open(1,file='params1.txt')
!      open(91,file='modobs.txt')
!      open(92,file='modobs2.txt')
      read(1,*)
      read(1,*) jglobal,ifsoc14,kinetics,bgcopt,jopt,jrestart
      read(1,11) finputsoc
11    format(a120)
      read(1,*) xopt(1:14)

      if(jopt==0) then
         read(1,*) nxopt(1:nx)
         do nparam=1,nx
            xopt(nxopt(nparam)) = xparam16(nparam)
         enddo
      endif
      print*, xopt

      mp = 4058  !wosis
      jmodel=1; mpft=17;mbgc=12
      nyeqpool= 1000

      isoc14 = 0

      totcost1 = 0.0

      call mic_allocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_allocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_allocate_output(mp,micoutput)
      call mic_allocate_cpool(mp,ms,miccpool)
      call mic_allocate_npool(mp,ms,micnpool)


 !         print *, "isoc14 =",isoc14,'--getdata_soc'
          call getdata_soc(finputsoc,micinput,micparam,micnpool)
          call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)
 !         print *, 'vmicsoil_soc'
          call vmicsoil_soc(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,bgcopt,nyeqpool, &
                        micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

 !         print *, 'calcost_soc'
          call calcost_soc(nx,isoc14,bgcopt,xopt,micpxdef,micparam,miccpool,micinput,totcost1)
          functn_soc_wosis = totcost1

      close(1)
!      close(91)
!      close(92)

!      functn = totcost

      call mic_deallocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_deallocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_deallocate_output(mp,micoutput)
      call mic_deallocate_cpool(mp,ms,miccpool)
      call mic_deallocate_npool(mp,ms,micnpool)

END function functn_soc_wosis


real*8 function functn(nx,xparam16)
 ! real*8 function functn_soc_hwsd(nx,xparam16)
 ! this function is yet to bet set up for running with SCE_UA optimization
   use mic_constant
   use mic_variable
   implicit none
    !local variables
    integer    :: nx
    integer,   dimension(16)           :: nxopt
    real*8,    dimension(16)           :: xparam16
    real*8,    dimension(16)           :: xopt
    TYPE(mic_param_xscale)    :: micpxdef
    TYPE(mic_param_default)   :: micpdef
    TYPE(mic_parameter)       :: micparam
    TYPE(mic_input)           :: micinput
    TYPE(mic_global_input)    :: micglobal
    TYPE(mic_cpool)           :: miccpool
    TYPE(mic_npool)           :: micnpool
    TYPE(mic_output)          :: micoutput

    !local variables
    integer    :: ifsoc14,kinetics,bgcopt,jopt,nyeqpool,isoc14,jglobal,jmodel
    integer :: jrestart,nf,ok,nparam,mpx,timex
    character(len=99)  :: frestart_in,frestart_out,fparam_global,foutput
    character(len=120) :: fhwsdsoc
    real(r_2)     :: totcost1

      isoc14=0;nyeqpool = 500;ok=0;totcost1=0.0

      jrestart=0;xopt(:)=1.0
      do nparam=1,16
         nxopt(nparam) = nparam
      enddo

      frestart_in='miccpool_in.nc'
      frestart_out='miccpool_out.nc'
      foutput='vmic_output.nc'

    !  open(91,file='modobs.txt')
    !  open(92,file='modobs2.txt')

      open(1,file='params1.txt')
      read(1,*)
      read(1,*) jglobal,ifsoc14,kinetics,bgcopt,jopt,jrestart,jmodel
      read(1,101) fhwsdsoc
      read(1,*)   xopt(1:14)
      read(1,*) nxopt(1:nx)
      do nparam=1,nx
         xopt(nxopt(nparam)) = xparam16(nparam)
      enddo
      close(1)
      print*, xopt
101   format(a120)

      ! get dimensions
      call getdata_hwsd_dim(fhwsdsoc,mpx,timex)
      mp=mpx
      ntime=timex
      if(jmodel==1) mpft=17   !CABLE
      if(jmodel==2) mpft=18   !ORCHIDEE
      mbgc=12

    !  print *, 'dimensions: mp ms ntime = ', mp,ms,ntime

      call mic_allocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_allocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_allocate_output(mp,micoutput)
      call mic_allocate_cpool(mp,ms,miccpool)
      call mic_allocate_npool(mp,ms,micnpool)
    !  print *, ' all  arrays are allocated!'

      call getdata_hwsd(fhwsdsoc,jglobal,bgcopt,jopt,jmodel,micparam,micglobal)
    !  print *, 'global input data are read in'
    !  print *, 'pft=', micglobal%pft

      call profile_hwsd(jmodel,micparam)
!      call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)

    !  print *, 'vmicsoil_hwsd'
!      call vmicsoil_hwsd_cpu(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,bgcopt,nyeqpool, &
!                         micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

    !  print *, 'calcost_hwsd'

!      call calcost_hwsd(nx,bgcopt,xopt,micpxdef,micparam,miccpool,micinput,micglobal,totcost1)

      call mic_deallocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_deallocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_deallocate_output(mp,micoutput)
      call mic_deallocate_cpool(mp,ms,miccpool)
      call mic_deallocate_npool(mp,ms,micnpool)

      close(1)
      print *, 'cost and parameter values',totcost1, xopt(nxopt(1:nx))
    !  close(91)
    !  close(92)
      functn = totcost1

    !  functn_soc_hwsd=real(ok)
!END function functn_soc_hwsd
END function functn

  subroutine vmic_restart_read(miccpool,micnpool,frestart_in)
  ! read soil carbon pool sizes "miccpool%cpool(mp,ms,mcpool)"
    use netcdf
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_cpool),              INTENT(INOUT)   :: miccpool
    TYPE(mic_npool),              INTENT(INOUT)   :: micnpool
    character :: frestart_in*99
    ! local variables
    integer :: mpx,msx,mcpoolx
    integer :: status,ncid,varid
    real(r_2), dimension(mp,ms,mcpool)  :: fcpool
    real(r_2), dimension(mp,ms)         :: fnpool

   ! open restart file
    status = nf90_open(frestart_in,nf90_nowrite,ncid)
    if(status /= nf90_noerr) CALL nc_abort(STATUS, 'Error opening '//frestart_in)

    ! get dimensions
    status = nf90_inq_dimid(ncid,'mp',varid)
    if(status /= nf90_noerr) CALL nc_abort(STATUS, 'Error inquiring dimensions mp_id')
    status = nf90_inquire_dimension(ncid,varid,len=mpx)
    if(status /= nf90_noerr) CALL nc_abort(STATUS, 'Error reading mp')

    status = nf90_inq_dimid(ncid,'ms',varid)
    if(status /= nf90_noerr) CALL nc_abort(STATUS,'Error inquiring dimensions ms_id')
    status = nf90_inquire_dimension(ncid,varid,len=msx)
    if(status /= nf90_noerr) CALL nc_abort(STATUS, 'Error reading ms')

    status = nf90_inq_dimid(ncid,'mcpool',varid)
    if(status /= nf90_noerr) CALL nc_abort(STATUS, 'Error inquiring dimensions mccpool_id')
    status = nf90_inquire_dimension(ncid,varid,len=mcpoolx)
    if(status /= nf90_noerr) CALL nc_abort(STATUS,'Error reading mcpool')

    ! get variables
    status = nf90_inq_varid(ncid,'mic_cpool',varid)
    if(status /= nf90_noerr) CALL nc_abort(STATUS, 'Error inquiring miccpoolc')
    status = nf90_get_var(ncid,varid,fcpool)
    if(status /= nf90_noerr) CALL nc_abort(STATUS,'Error reading fcpool')

    status = nf90_inq_varid(ncid,'mic_npool',varid)
    if(status /= nf90_noerr) CALL nc_abort(STATUS, 'Error inquiring micnpoolc')
    status = nf90_get_var(ncid,varid,fnpool)
    if(status /= nf90_noerr) CALL nc_abort(STATUS,'Error reading fnpool')

    ! close the file
    status = NF90_close(ncid)
    if(status /= nf90_noerr) call nc_abort(status, 'Error in clsoing netCDF input file')

    ! assign the values from the restart file
    if(mpx/=mp .or. msx/=ms .or. mcpoolx/=mcpool) then
       print *, 'dimensions do not match! ', mp,mpx,ms,msx,mcpool,mcpoolx
       STOP
    endif
    miccpool%cpool    = fcpool
    micnpool%mineralN = fnpool


  end subroutine vmic_restart_read


  subroutine vmic_restart_write(frestart_out,miccpool,micnpool)
  ! write out soil carbon pool sizes "miccpool%cpool(mp,ms,mcpool)"
    use netcdf
    use mic_constant
    use mic_variable
    implicit None
    TYPE(mic_cpool),              INTENT(INOUT)   :: miccpool
    TYPE(mic_npool),              INTENT(INOUT)   :: micnpool
    INTEGER*4                :: STATUS
    INTEGER*4                :: FILE_ID, mp_ID, miccarb_ID, soil_ID
    CHARACTER                :: CDATE*10,frestart_out*99
    INTEGER*4                :: cmic_ID, nmic_ID
    integer :: values(10)
    real(r_2)  :: missreal

    missreal=-1.0e10
    call date_and_time(values=values)
    WRITE(CDATE, '(I4.4,"-",I2.2,"-",I2.2)') values(1),values(2),values(3)

    ! Create NetCDF file:
    STATUS = NF90_create(frestart_out, NF90_CLOBBER, FILE_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error creating restart file ')

    WRITE(*,*) 'writing mic restart', frestart_out
    ! Put the file in define mode:
    STATUS = NF90_redef(FILE_ID)

    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, 'Valid restart date', CDATE )

    ! Define dimensions:
    ! mp (number of patches)
    STATUS = NF90_def_dim(FILE_ID, 'mp'   , mp     , mp_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error defining mp dimension ')

    ! ms: number of soil layers
    STATUS = NF90_DEF_DIM(FILE_ID, 'ms', ms, soil_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error defining soil dimension ' )

    ! mcpool: number of soil carbon pools
    STATUS = NF90_def_dim(FILE_ID, 'mcpool', mcpool, miccarb_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error defining mic_carbon_pools dimension ' )

    STATUS = NF90_def_var(FILE_ID,'mic_cpool',NF90_FLOAT,[mp_ID,soil_ID,miccarb_ID],cmic_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error defining mic_cpool variable ' )

    STATUS = NF90_def_var(FILE_ID,'mic_npool',NF90_FLOAT,[mp_ID,soil_ID],nmic_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error defining mic_npool variable ' )

    ! End define mode:
    STATUS = NF90_enddef(FILE_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error ending define mode ' )

    ! PUT VARS
    STATUS = NF90_PUT_VAR(FILE_ID, cmic_ID, REAL(miccpool%cpool, 4) )
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error writing mic_cpool variable ' )

    STATUS = NF90_PUT_VAR(FILE_ID, nmic_ID, REAL(micnpool%mineralN, 4) )
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error writing mic_npool variable ')

    ! Close NetCDF file:
    STATUS = NF90_close(FILE_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error closing restart file '  )

    write(*, *) 'restart file written to ', frestart_out

  end subroutine vmic_restart_write

  SUBROUTINE nc_abort( ok, message )
    USE netcdf
    implicit none
    ! Input arguments
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER, INTENT(IN) :: ok

    WRITE(*,*) message ! error from subroutine
    WRITE(*,*) NF90_STRERROR(ok) ! netcdf error details

    STOP

  END SUBROUTINE nc_abort

  subroutine vmic_output_write(foutput,micinput,micoutput)
    ! fNPP is not quite right yet. It shoudl be the sump of "cinputm+cinputs"
    use netcdf
    use mic_constant
    use mic_variable
    implicit None
    TYPE(mic_input),         INTENT(INout)   :: micinput
    TYPE(mic_output),        INTENT(INout)   :: micoutput
    real(r_2)     :: missreal
    INTEGER*4                :: STATUS
    INTEGER*4                :: FILE_ID, mp_ID
    CHARACTER                :: CDATE*10,foutput*99
    INTEGER*4                :: cinput_ID, rsoil_ID, cleach_ID
    integer :: values(10)

    missreal=-1.0e10
    call date_and_time(values=values)
    WRITE(CDATE, '(I4.4,"-",I2.2,"-",I2.2)') values(1),values(2),values(3)
    ! Create NetCDF file:
    STATUS = NF90_create(foutput, NF90_CLOBBER, FILE_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error creating output file ')

    WRITE(*,*) 'writing output file', foutput
    print *, CDATE

    ! Put the file in define mode:
    STATUS = NF90_redef(FILE_ID)

    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, 'Valid output date', CDATE  )

    ! Define dimensions:
    ! mp (number of patches)
    STATUS = NF90_def_dim(FILE_ID, 'mp'   , mp     , mp_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error defining mp dimension ')

    STATUS = NF90_def_var(FILE_ID,'Cinput',NF90_FLOAT,[mp_ID],cinput_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error defining NPP ' )


    STATUS = NF90_def_var(FILE_ID,'rsoil',NF90_FLOAT,[mp_ID],rsoil_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error defining rsoil ' )


    STATUS = NF90_def_var(FILE_ID,'Cleach',NF90_FLOAT,[mp_ID],cleach_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error defining cleach ' )

    ! End define mode:
    STATUS = NF90_enddef(FILE_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error ending define mode ' )

    ! put attributes
    STATUS = NF90_PUT_ATT(FILE_ID,cinput_ID,'unit','g C m-2 year-1')
    STATUS = NF90_PUT_ATT(FILE_ID,cinput_ID,'missing_value', real(missreal,4))

    STATUS = NF90_PUT_ATT(FILE_ID,rsoil_ID,'unit','g C m-2 year-1')
    STATUS = NF90_PUT_ATT(FILE_ID,rsoil_ID,'missing_value', real(missreal,4))

    STATUS = NF90_PUT_ATT(FILE_ID,cleach_ID,'unit','g C m-2 year-1')
    STATUS = NF90_PUT_ATT(FILE_ID,cleach_ID,'missing_value', real(missreal,4))

    ! PUT VARS
    STATUS = NF90_PUT_VAR(FILE_ID, cinput_ID, REAL(micoutput%fluxcinput,4) )
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error writing NPP ' )

    STATUS = NF90_PUT_VAR(FILE_ID, rsoil_ID, REAL(micoutput%fluxrsoil,4) )
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error writing Rsoil ')

    STATUS = NF90_PUT_VAR(FILE_ID, cleach_ID, REAL(micoutput%fluxcleach,4) )
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error writing Cleach ')

    ! Close NetCDF file:
    STATUS = NF90_close(FILE_ID)
    IF(STATUS /= NF90_NOERR) CALL nc_abort(STATUS, 'Error closing restart file '  )

    write(*, *) 'output written to ', foutput

  end subroutine vmic_output_write

  subroutine getparam_global(micpxdef)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_param_xscale)    :: micpxdef
    integer :: npft,ipft,n
    real(r_2), dimension(20)    :: x

    open(100,file='parameters_global.csv')
    read(100,*)
    do npft=1,mpft
       read(100,*) ipft, (x(n),n=1,17)


       micpxdef%xav(npft)        = x(1)
       micpxdef%xak(npft)        = x(2)
       micpxdef%xfm(npft)        = x(3)
       micpxdef%xfs(npft)        = x(4)
       micpxdef%xtvmic(npft)     = x(5)
       micpxdef%xtvp(npft)       = x(6)
       micpxdef%xtvc(npft)       = x(7)
       micpxdef%xtvac(npft)      = x(8)
       micpxdef%xkba(npft)       = x(9)
       micpxdef%xqmaxcoeff(npft) = x(10)
       micpxdef%xdiffsoc(npft)   = x(11)
       micpxdef%xnpp(npft)       = x(12)
       micpxdef%xrootbeta(npft)  = x(13)
       micpxdef%xvmaxbeta(npft)  = x(14)
       ! the following parameters are fixed to 1.0
       micpxdef%xfp2ax(npft)     = x(15)
       micpxdef%xbeta(npft)      = x(16)
       micpxdef%xdesorp(npft)    = x(17)
    enddo
    close(100)

!    print *, 'xdesorp=', micpxdef%xdesorp(:)
!    print *, 'xtvc=',    micpxdef%xtvc(npft)

  end subroutine getparam_global

  subroutine getdata_c14(frac14c,f14c,micinput,micparam,micnpool)
    use netcdf
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(INout)   :: micparam
    TYPE(mic_input),     INTENT(INout)   :: micinput
    TYPE(mic_npool),     INTENT(INOUT)   :: micnpool

    integer:: ncid,varid,status
    integer:: np,ns,i,j
    integer:: nz
    character(len=120) :: frac14c,f14c(5)

    character(len = nf90_max_name):: name
    real(r_2),dimension(:,:),allocatable:: fclay,fsilt,fph,ftemp,fmoist,fporosity,fmatpot
    real(r_2),dimension(:),allocatable:: fsoc,fpoc,fmaoc,ffmpoc,ffmmaoc,fbulkd
    real(r_2),dimension(:),allocatable:: fnpp,fanpp,fbnpp,flignin,fcna,fcnb
    integer,dimension(:),allocatable:: fid,fpft,ftop,fbot,fyear,fregion

    allocate(fsoc(mp))

    allocate(fclay(mp,ms))
    allocate(fsilt(mp,ms))
    allocate(fph(mp,ms))
    allocate(ftemp(mp,ms))
    allocate(fmoist(mp,ms))
    allocate(fporosity(mp,ms))
    allocate(fmatpot(mp,ms))

    allocate(fnpp(mp))
    allocate(fanpp(mp))
    allocate(fbnpp(mp))
    allocate(flignin(mp))
    allocate(fcna(mp))
    allocate(fcnb(mp))
    allocate(fid(mp))
    allocate(fpft(mp))

    ! inputdata for 14C
    allocate(fpoc(mp))
    allocate(fmaoc(mp))
    allocate(ffmpoc(mp))
    allocate(ffmmaoc(mp))
    allocate(fbulkd(mp))
    allocate(ftop(mp)) !! upper depth of observed soil layer
    allocate(fbot(mp)) !! lower depth of observed soil layer
    allocate(fyear(mp)) !! year at which c14 was observed
    allocate(fregion(mp)) !! north/south hemisphere zone of c14

   ! open .nc file
    status = nf90_open(frac14c,nf90_nowrite,ncid)
    if(status /= nf90_noerr) print*, 'Error opening frc_c14.nc'

    ! get dimensions/profile_id
    status = nf90_inq_varid(ncid,'nsite',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring dimensions/profile_id'
    status = nf90_get_var(ncid,varid,fid)
    if(status /= nf90_noerr) print*,'Error reading profile_id'

    ! get variables
    status = nf90_inq_varid(ncid,'SOC',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soc'
    status = nf90_get_var(ncid,varid,fsoc)
    if(status /= nf90_noerr) print*,'Error reading soc'

    status = nf90_inq_varid(ncid,'bulkd',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring bulk density'
    status = nf90_get_var(ncid,varid,fbulkd)
    if(status /= nf90_noerr) print*,'Error reading bulk density'

    status = nf90_inq_varid(ncid,'clay',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring clay'
    status = nf90_get_var(ncid,varid,fclay)
    if(status /= nf90_noerr) print*,'Error reading clay'

    status = nf90_inq_varid(ncid,'silt',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring silt'
    status = nf90_get_var(ncid,varid,fsilt)
    if(status /= nf90_noerr) print*,'Error reading silt'

    status = nf90_inq_varid(ncid,'ph',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring ph'
    status = nf90_get_var(ncid,varid,fph)
    if(status /= nf90_noerr) print*,'Error reading ph'

    status = nf90_inq_varid(ncid,'temp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil temperature'
    status = nf90_get_var(ncid,varid,ftemp)
    if(status /= nf90_noerr) print*,'Error reading soil temperature'

    status = nf90_inq_varid(ncid,'moist',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil moisture'
    status = nf90_get_var(ncid,varid,fmoist)
    if(status /= nf90_noerr) print*,'Error reading soil moisture'

    status = nf90_inq_varid(ncid,'porosity',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil porosity'
    status = nf90_get_var(ncid,varid,fporosity)
    if(status /= nf90_noerr) print*,'Error reading soil porosity'

    status = nf90_inq_varid(ncid,'matpot',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil matric potential'
    status = nf90_get_var(ncid,varid,fmatpot)
    if(status /= nf90_noerr) print*,'Error reading soil matric potential'

    status = nf90_inq_varid(ncid,'npp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring npp'
    status = nf90_get_var(ncid,varid,fnpp)
    if(status /= nf90_noerr) print*,'Error reading npp'

    status = nf90_inq_varid(ncid,'anpp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring anpp'
    status = nf90_get_var(ncid,varid,fanpp)
    if(status /= nf90_noerr) print*,'Error reading anpp'

    status = nf90_inq_varid(ncid,'bnpp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring bnpp'
    status = nf90_get_var(ncid,varid,fbnpp)
    if(status /= nf90_noerr) print*,'Error reading bnpp'

    status = nf90_inq_varid(ncid,'lignin_C',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring lignin/C'
    status = nf90_get_var(ncid,varid,flignin)
    if(status /= nf90_noerr) print*,'Error reading lignin/C'

    status = nf90_inq_varid(ncid,'cna',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring C/N aboveground'
    status = nf90_get_var(ncid,varid,fcna)
    if(status /= nf90_noerr) print*,'Error reading C/N aboveground'

    status = nf90_inq_varid(ncid,'cnb',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring C/N belowground'
    status = nf90_get_var(ncid,varid,fcnb)
    if(status /= nf90_noerr) print*,'Error reading C/N belowground'

    status = nf90_inq_varid(ncid,'pft',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring plant functional type'
    status = nf90_get_var(ncid,varid,fpft)
    if(status /= nf90_noerr) print*,'Error reading plant functional type'

      status = nf90_inq_varid(ncid,'POC',varid)
      if(status /= nf90_noerr) print*, 'Error inquiring POC'
      status = nf90_get_var(ncid,varid,fpoc)
      if(status /= nf90_noerr) print*,'Error reading POC'

      status = nf90_inq_varid(ncid,'MAOC',varid)
      if(status /= nf90_noerr) print*, 'Error inquiring MAOC'
      status = nf90_get_var(ncid,varid,fmaoc)
      if(status /= nf90_noerr) print*,'Error reading MAOC'

      status = nf90_inq_varid(ncid,'fm_poc',varid)
      if(status /= nf90_noerr) print*, 'Error inquiring fm_poc'
      status = nf90_get_var(ncid,varid,ffmpoc)
      if(status /= nf90_noerr) print*,'Error reading fm_poc'

      status = nf90_inq_varid(ncid,'fm_maoc',varid)
      if(status /= nf90_noerr) print*, 'Error inquiring fm_maoc'
      status = nf90_get_var(ncid,varid,ffmmaoc)
      if(status /= nf90_noerr) print*,'Error reading fm_maoc'

      status = nf90_inq_varid(ncid,'top_depth',varid)
      if(status /= nf90_noerr) print*, 'Error inquiring top depth'
      status = nf90_get_var(ncid,varid,ftop)
      if(status /= nf90_noerr) print*,'Error reading top depth'

      status = nf90_inq_varid(ncid,'bot_depth',varid)
      if(status /= nf90_noerr) print*, 'Error inquiring bottom depth'
      status = nf90_get_var(ncid,varid,fbot)
      if(status /= nf90_noerr) print*,'Error reading bottom depth'

      status = nf90_inq_varid(ncid,'c14_year',varid)
      if(status /= nf90_noerr) print*, 'Error inquiring c14 year'
      status = nf90_get_var(ncid,varid,fyear)
      if(status /= nf90_noerr) print*,'Error reading c14 year'

      status = nf90_inq_varid(ncid,'c14_region',varid)
      if(status /= nf90_noerr) print*, 'Error inquiring c14 region'
      status = nf90_get_var(ncid,varid,fregion)
      if(status /= nf90_noerr) print*,'Error reading c14 region'

    ! Close netcdf file
    status = NF90_CLOSE(ncid)

      ! we need to include additional data for kinetics3

      micparam%csoilobs(:,:) = -999.0
      do np=1, mp

         micparam%pft(np)    = int(fpft(np))
         micparam%siteid(np) = int(fid(np))

            micparam%top(np)         = int(ftop(np))
            micparam%bot(np)         = int(fbot(np))
            micparam%nyc14obs(np)    = int(fyear(np)) !! year when c14 is observed
            micparam%region(np)      = int(fregion(np)) !! south/north hemisphere zone of c14
            micparam%c14soilobsp(np) = ffmpoc(np) !! poc c14 fraction modern
            micparam%c14soilobsm(np) = ffmmaoc(np) !! maoc c14 fraction modern

         ! make sure "*delt" is not repeated in the model called by rk4
          micinput%fcnpp(np)      = fnpp(np)
          micinput%Dleaf(np)      = fanpp(np)/(24.0*365.0)*delt    !gc/m2/delt
          micinput%Droot(np)      = fbnpp(np)/(24.0*365.0)*delt     !gc/m2/delt
          !micinput%Dwood(np)      = forcdata(np,17)/(24.0*365.0)*delt     !gc/m2/delt

          micparam%xcnleaf(np)    = fcna(np)
          micparam%xcnroot(np)    = fcnb(np)
          !micparam%xcnwood(np)    = forcdata(np,20)
          micparam%fligleaf(np)   = flignin(np)
          micparam%fligroot(np)   = flignin(np)
          !micparam%fligwood(np)   = forcdata(np,23)

         do ns=1,ms
            micinput%tavg(np,ns)     = ftemp(np,ns)  ! average temperature in deg C
            micinput%wavg(np,ns)     = fmoist(np,ns)  ! average soil water content mm3/mm3
            micinput%clay(np,ns)     = fclay(np,ns)  ! clay content (fraction)
            micinput%silt(np,ns)     = fsilt(np,ns)  ! silt content (fraction)
            micinput%ph(np,ns)       = fph(np,ns)
            micinput%porosity(np,ns) = fporosity(np,ns) !porosity mm3/mm3
            micinput%matpot(np,ns)   = fmatpot(np,ns)  ! soil matric potential -kPa

            micparam%csoilobs(np,ns)    = fsoc(np)
            micinput%bulkd(np,ns)       = fbulkd(np)

            micparam%csoilobsp(np,ns)   = fpoc(np)
            micparam%csoilobsm(np,ns)   = fmaoc(np)

            !micnpool%mineralN(np,ns) = forcdata(np,7)*0.001 ! mineral N: "0.001" mg N /kg soil --> g N /kg soil
         enddo !"ns"
      enddo    ! "np=1,mp"

         ! read in the standard 14C atmospheric data for five zones
!         f14c(1) ='/g/data/w97/lw9370/combined-model/c14/code-structure/data/NH1-C14.csv'
!         f14c(2) ='/g/data/w97/lw9370/combined-model/c14/code-structure/data/NH2-C14.csv'
!         f14c(3) ='/g/data/w97/lw9370/combined-model/c14/code-structure/data/NH3-C14.csv'
!         f14c(4) ='/g/data/w97/lw9370/combined-model/c14/code-structure/data/SH12-C14.csv'
!         f14c(5) ='/g/data/w97/lw9370/combined-model/c14/code-structure/data/SH3-C14.csv'
         do nz=1,5
             call get14catm(nz,f14c(nz),micparam)
         enddo

    deallocate(fsoc)
    deallocate(fbulkd)
    deallocate(fclay)
    deallocate(fsilt)
    deallocate(fph)
    deallocate(ftemp)
    deallocate(fmoist)
    deallocate(fporosity)
    deallocate(fmatpot)

    deallocate(fnpp)
    deallocate(fanpp)
    deallocate(fbnpp)
    deallocate(flignin)
    deallocate(fcna)
    deallocate(fcnb)
    deallocate(fid)
    deallocate(fpft)

    deallocate(fpoc)
    deallocate(fmaoc)
    deallocate(ffmpoc)
    deallocate(ffmmaoc)
    deallocate(ftop) !! upper depth of observed soil layer
    deallocate(fbot) !! bottom depth of observed soil layer
    deallocate(fyear) !! year at which c14 was observed
    deallocate(fregion) !! north/south hemisphere zone of c14

   end subroutine getdata_c14

   subroutine getdata_frc(Cfraction,micinput,micparam,micnpool)
    use netcdf
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(INout)   :: micparam
    TYPE(mic_input),     INTENT(INout)   :: micinput
    TYPE(mic_npool),     INTENT(INOUT)   :: micnpool

    integer:: ncid,varid,status
    integer:: np,ns,i,j
    integer:: nz
    character(len=120) :: Cfraction

    character(len = nf90_max_name):: name
    real(r_2),dimension(:),allocatable:: fclay,fsilt,fph,ftemp,fmoist,fporosity,fmatpot
    real(r_2),dimension(:),allocatable:: fsoc,fpoc,fmaoc,fbulkd
    real(r_2),dimension(:),allocatable:: fnpp,fanpp,fbnpp,flignin,fcna,fcnb
    integer,dimension(:),allocatable:: fid,fpft,ftop,fbot,fdataid

    allocate(fsoc(mp))

    allocate(fclay(mp))
    allocate(fsilt(mp))
    allocate(fph(mp))
    allocate(ftemp(mp))
    allocate(fmoist(mp))
    allocate(fporosity(mp))
    allocate(fmatpot(mp))

    allocate(fnpp(mp))
    allocate(fanpp(mp))
    allocate(fbnpp(mp))
    allocate(flignin(mp))
    allocate(fcna(mp))
    allocate(fcnb(mp))
    allocate(fid(mp))
    allocate(fpft(mp))

    ! inputdata for 14C
    allocate(fpoc(mp))
    allocate(fmaoc(mp))
    allocate(fbulkd(mp))
    allocate(ftop(mp)) !! upper depth of observed soil layer
    allocate(fbot(mp)) !! lower depth of observed soil layer
    allocate(fdataid(mp)) !! 1 for LUCAS; 2 for AUS; 3 for KG

   ! open .nc file
    status = nf90_open(Cfraction,nf90_nowrite,ncid)
    if(status /= nf90_noerr) print*, 'Error opening c_fraction.nc'

    ! get dimensions/profile_id
    status = nf90_inq_varid(ncid,'nsite',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring dimensions/profile_id'
    status = nf90_get_var(ncid,varid,fid)
    if(status /= nf90_noerr) print*,'Error reading profile_id'

    ! get variables
    status = nf90_inq_varid(ncid,'dataid',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring data ID'
    status = nf90_get_var(ncid,varid,fdataid)
    if(status /= nf90_noerr) print*,'Error reading data ID'

    status = nf90_inq_varid(ncid,'SOC',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soc'
    status = nf90_get_var(ncid,varid,fsoc)
    if(status /= nf90_noerr) print*,'Error reading soc'

    status = nf90_inq_varid(ncid,'bulkd',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring bulk density'
    status = nf90_get_var(ncid,varid,fbulkd)
    if(status /= nf90_noerr) print*,'Error reading bulk density'

    status = nf90_inq_varid(ncid,'clay',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring clay'
    status = nf90_get_var(ncid,varid,fclay)
    if(status /= nf90_noerr) print*,'Error reading clay'

    status = nf90_inq_varid(ncid,'silt',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring silt'
    status = nf90_get_var(ncid,varid,fsilt)
    if(status /= nf90_noerr) print*,'Error reading silt'

    status = nf90_inq_varid(ncid,'ph',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring ph'
    status = nf90_get_var(ncid,varid,fph)
    if(status /= nf90_noerr) print*,'Error reading ph'

    status = nf90_inq_varid(ncid,'temp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil temperature'
    status = nf90_get_var(ncid,varid,ftemp)
    if(status /= nf90_noerr) print*,'Error reading soil temperature'

    status = nf90_inq_varid(ncid,'moist',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil moisture'
    status = nf90_get_var(ncid,varid,fmoist)
    if(status /= nf90_noerr) print*,'Error reading soil moisture'

    status = nf90_inq_varid(ncid,'porosity',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil porosity'
    status = nf90_get_var(ncid,varid,fporosity)
    if(status /= nf90_noerr) print*,'Error reading soil porosity'

    status = nf90_inq_varid(ncid,'matpot',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil matric potential'
    status = nf90_get_var(ncid,varid,fmatpot)
    if(status /= nf90_noerr) print*,'Error reading soil matric potential'

    status = nf90_inq_varid(ncid,'npp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring npp'
    status = nf90_get_var(ncid,varid,fnpp)
    if(status /= nf90_noerr) print*,'Error reading npp'

    status = nf90_inq_varid(ncid,'anpp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring anpp'
    status = nf90_get_var(ncid,varid,fanpp)
    if(status /= nf90_noerr) print*,'Error reading anpp'

    status = nf90_inq_varid(ncid,'bnpp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring bnpp'
    status = nf90_get_var(ncid,varid,fbnpp)
    if(status /= nf90_noerr) print*,'Error reading bnpp'

    status = nf90_inq_varid(ncid,'lignin_C',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring lignin/C'
    status = nf90_get_var(ncid,varid,flignin)
    if(status /= nf90_noerr) print*,'Error reading lignin/C'

    status = nf90_inq_varid(ncid,'cna',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring C/N aboveground'
    status = nf90_get_var(ncid,varid,fcna)
    if(status /= nf90_noerr) print*,'Error reading C/N aboveground'

    status = nf90_inq_varid(ncid,'cnb',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring C/N belowground'
    status = nf90_get_var(ncid,varid,fcnb)
    if(status /= nf90_noerr) print*,'Error reading C/N belowground'

    status = nf90_inq_varid(ncid,'pft',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring plant functional type'
    status = nf90_get_var(ncid,varid,fpft)
    if(status /= nf90_noerr) print*,'Error reading plant functional type'

    status = nf90_inq_varid(ncid,'POC',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring POC'
    status = nf90_get_var(ncid,varid,fpoc)
    if(status /= nf90_noerr) print*,'Error reading POC'

    status = nf90_inq_varid(ncid,'MAOC',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring MAOC'
    status = nf90_get_var(ncid,varid,fmaoc)
    if(status /= nf90_noerr) print*,'Error reading MAOC'

    status = nf90_inq_varid(ncid,'top_depth',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring top depth'
    status = nf90_get_var(ncid,varid,ftop)
    if(status /= nf90_noerr) print*,'Error reading top depth'

    status = nf90_inq_varid(ncid,'bot_depth',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring bottom depth'
    status = nf90_get_var(ncid,varid,fbot)
    if(status /= nf90_noerr) print*,'Error reading bottom depth'

    ! Close netcdf file
    status = NF90_CLOSE(ncid)

      ! we need to include additional data for kinetics3

      micparam%csoilobs(:,:) = -999.0
      do np=1, mp

         micparam%pft(np)    = int(fpft(np))
         micparam%siteid(np) = int(fid(np))
         micparam%dataid(np) = int(fdataid(np))
         micparam%top(np)         = int(ftop(np))
         micparam%bot(np)         = int(fbot(np))

         ! make sure "*delt" is not repeated in the model called by rk4
          micinput%fcnpp(np)      = fnpp(np)
          micinput%Dleaf(np)      = fanpp(np)/(24.0*365.0)*delt    !gc/m2/delt
          micinput%Droot(np)      = fbnpp(np)/(24.0*365.0)*delt     !gc/m2/delt
          !micinput%Dwood(np)      = forcdata(np,17)/(24.0*365.0)*delt     !gc/m2/delt

          micparam%xcnleaf(np)    = fcna(np)
          micparam%xcnroot(np)    = fcnb(np)
          !micparam%xcnwood(np)    = forcdata(np,20)
          micparam%fligleaf(np)   = flignin(np)
          micparam%fligroot(np)   = flignin(np)
          !micparam%fligwood(np)   = forcdata(np,23)

         do ns=1,ms
            micinput%tavg(np,ns)     = ftemp(np)  ! average temperature in deg C
            micinput%wavg(np,ns)     = fmoist(np)  ! average soil water content mm3/mm3
            micinput%clay(np,ns)     = fclay(np)  ! clay content (fraction)
            micinput%silt(np,ns)     = fsilt(np)  ! silt content (fraction)
            micinput%ph(np,ns)       = fph(np)
            micinput%porosity(np,ns) = fporosity(np) !porosity mm3/mm3
            micinput%matpot(np,ns)   = fmatpot(np)  ! soil matric potential -kPa

            micparam%csoilobs(np,ns)    = fsoc(np)
            micinput%bulkd(np,ns)       = fbulkd(np)

            micparam%csoilobsp(np,ns)   = fpoc(np)
            micparam%csoilobsm(np,ns)   = fmaoc(np)

            !micnpool%mineralN(np,ns) = forcdata(np,7)*0.001 ! mineral N: "0.001" mg N /kg soil --> g N /kg soil
         enddo !"ns"
      enddo    ! "np=1,mp"

    deallocate(fsoc)
    deallocate(fbulkd)
    deallocate(fclay)
    deallocate(fsilt)
    deallocate(fph)
    deallocate(ftemp)
    deallocate(fmoist)
    deallocate(fporosity)
    deallocate(fmatpot)

    deallocate(fnpp)
    deallocate(fanpp)
    deallocate(fbnpp)
    deallocate(flignin)
    deallocate(fcna)
    deallocate(fcnb)
    deallocate(fid)
    deallocate(fpft)

    deallocate(fpoc)
    deallocate(fmaoc)
    deallocate(ftop) !! upper depth of observed soil layer
    deallocate(fbot) !! bottom depth of observed soil layer
    deallocate(fdataid)

   end subroutine getdata_frc

   subroutine getdata_soc(finputsoc,micinput,micparam,micnpool)
    use netcdf
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(INout)   :: micparam
    TYPE(mic_input),     INTENT(INout)   :: micinput
    TYPE(mic_npool),     INTENT(INOUT)   :: micnpool

    character(len=120) :: finputsoc
    integer:: ncid,varid,status
    integer:: np,ns,i,j
    character(len = nf90_max_name):: name
    real(r_2),dimension(:,:),allocatable:: fclay,fsilt,fph,ftemp,fmoist,fbulkd,fmatpot,fporosity,fsoc
    real(r_2),dimension(:),allocatable:: fnpp,fanpp,fbnpp,flignin,fcna,fcnb
    integer,dimension(:),allocatable:: fid,fpft,fcluster

    allocate(fsoc(mp,ms))
    allocate(fbulkd(mp,ms))
    allocate(fclay(mp,ms))
    allocate(fsilt(mp,ms))
    allocate(fph(mp,ms))
    allocate(ftemp(mp,ms))
    allocate(fmoist(mp,ms))
    allocate(fmatpot(mp,ms))
    allocate(fporosity(mp,ms))
    allocate(fnpp(mp))
    allocate(fanpp(mp))
    allocate(fbnpp(mp))
    allocate(flignin(mp))
    allocate(fcna(mp))
    allocate(fcnb(mp))
    allocate(fid(mp))
    allocate(fpft(mp))
    allocate(fcluster(mp))

 !   print *, 'mp ms',mp,ms

    ! open .nc file
    status = nf90_open(finputsoc,nf90_nowrite,ncid)
    if(status /= nf90_noerr) print*, 'Error opening wosis_input.nc'

    ! get dimensions/profile_id
    status = nf90_inq_varid(ncid,'nsite',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring dimensions/profile_id'
    status = nf90_get_var(ncid,varid,fid)
    if(status /= nf90_noerr) print*,'Error reading profile_id'

    ! get variables
    status = nf90_inq_varid(ncid,'soc',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soc'
    status = nf90_get_var(ncid,varid,fsoc)
    if(status /= nf90_noerr) print*,'Error reading soc'

    status = nf90_inq_varid(ncid,'bulkd',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring bulk density'
    status = nf90_get_var(ncid,varid,fbulkd)
    if(status /= nf90_noerr) print*,'Error reading bulk density'

    status = nf90_inq_varid(ncid,'clay',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring clay'
    status = nf90_get_var(ncid,varid,fclay)
    if(status /= nf90_noerr) print*,'Error reading clay'

    status = nf90_inq_varid(ncid,'silt',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring silt'
    status = nf90_get_var(ncid,varid,fsilt)
    if(status /= nf90_noerr) print*,'Error reading silt'

    status = nf90_inq_varid(ncid,'ph',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring ph'
    status = nf90_get_var(ncid,varid,fph)
    if(status /= nf90_noerr) print*,'Error reading ph'

    status = nf90_inq_varid(ncid,'temp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil temperature'
    status = nf90_get_var(ncid,varid,ftemp)
    if(status /= nf90_noerr) print*,'Error reading soil temperature'

    status = nf90_inq_varid(ncid,'moist',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil moisture'
    status = nf90_get_var(ncid,varid,fmoist)
    if(status /= nf90_noerr) print*,'Error reading soil moisture'

    status = nf90_inq_varid(ncid,'matpot',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil matric potential'
    status = nf90_get_var(ncid,varid,fmatpot)
    if(status /= nf90_noerr) print*,'Error reading soil matric potential'

    status = nf90_inq_varid(ncid,'porosity',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil porosity'
    status = nf90_get_var(ncid,varid,fporosity)
    if(status /= nf90_noerr) print*,'Error reading soil porosity'

    status = nf90_inq_varid(ncid,'npp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring npp'
    status = nf90_get_var(ncid,varid,fnpp)
    if(status /= nf90_noerr) print*,'Error reading npp'

    status = nf90_inq_varid(ncid,'anpp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring anpp'
    status = nf90_get_var(ncid,varid,fanpp)
    if(status /= nf90_noerr) print*,'Error reading anpp'

    status = nf90_inq_varid(ncid,'bnpp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring bnpp'
    status = nf90_get_var(ncid,varid,fbnpp)
    if(status /= nf90_noerr) print*,'Error reading bnpp'

    status = nf90_inq_varid(ncid,'lignin_C',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring lignin/C'
    status = nf90_get_var(ncid,varid,flignin)
    if(status /= nf90_noerr) print*,'Error reading lignin/C'

    status = nf90_inq_varid(ncid,'cna',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring C/N aboveground'
    status = nf90_get_var(ncid,varid,fcna)
    if(status /= nf90_noerr) print*,'Error reading C/N aboveground'

    status = nf90_inq_varid(ncid,'cnb',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring C/N belowground'
    status = nf90_get_var(ncid,varid,fcnb)
    if(status /= nf90_noerr) print*,'Error reading C/N belowground'

    status = nf90_inq_varid(ncid,'pft',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring plant functional type'
    status = nf90_get_var(ncid,varid,fpft)
    if(status /= nf90_noerr) print*,'Error reading plant functional type'

    status = nf90_inq_varid(ncid,'cluster',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring environmental cluster'
    status = nf90_get_var(ncid,varid,fcluster)
    if(status /= nf90_noerr) print*,'Error reading environmental cluster'

   !  print*,'matpot = ',fmatpot(1:4,1)
   !  print*,'siteid = ',fid(1:4)
   !  print*,'soc = ',fsoc(1:4,1)

   !  Close netcdf file

    status = NF90_CLOSE(ncid)

      micparam%csoilobs(:,:) = -999.0
 !     print *, 'size and pft', size(fcluster), fcluster(:)

      do np=1, mp

         ! use cluster
         micparam%pft(np)    = int(fcluster(np))
         ! use PFT
         ! micparam%pft(np)    = int(fpft(np))
         !
         micparam%siteid(np) = int(fid(np))

         ! make sure "*delt" is not repeated in the model called by rk4
          micinput%fcnpp(np)      = fnpp(np)
          micinput%Dleaf(np)      = fanpp(np)/(24.0*365.0)*delt    !gc/m2/delt
          micinput%Droot(np)      = fbnpp(np)/(24.0*365.0)*delt     !gc/m2/delt
          !micinput%Dwood(np)      = forcdata(np,17)/(24.0*365.0)*delt     !gc/m2/delt

          micparam%xcnleaf(np)    = fcna(np)
          micparam%xcnroot(np)    = fcnb(np)
          !micparam%xcnwood(np)    = forcdata(np,20)
          micparam%fligleaf(np)   = flignin(np)
          micparam%fligroot(np)   = flignin(np)
          !micparam%fligwood(np)   = forcdata(np,23)

         do ns=1,ms
            micinput%tavg(np,ns)     = ftemp(np,ns)  ! average temperature in deg C
            micinput%wavg(np,ns)     = fmoist(np,ns)  ! average soil water content mm3/mm3
            micinput%clay(np,ns)     = fclay(np,ns)  ! clay content (fraction)
            micinput%silt(np,ns)     = fsilt(np,ns)  ! silt content (fraction)
            micinput%ph(np,ns)       = fph(np,ns)
            micinput%bulkd(np,ns)    = fbulkd(np,ns)
            micinput%porosity(np,ns) = fporosity(np,ns) !porosity mm3/mm3
            micinput%matpot(np,ns)   = fmatpot(np,ns) ! -kPa
            micparam%csoilobs(np,ns) = fsoc(np,ns)

            !micnpool%mineralN(np,ns) = forcdata(np,7)*0.001 ! mineral N: "0.001" mg N /kg soil --> g N /kg soil
         enddo !"ns"
      enddo    ! "np=1,mp"

    deallocate(fsoc)
    deallocate(fbulkd)
    deallocate(fclay)
    deallocate(fsilt)
    deallocate(fph)
    deallocate(ftemp)
    deallocate(fmoist)
    deallocate(fporosity)
    deallocate(fmatpot)

    deallocate(fnpp)
    deallocate(fanpp)
    deallocate(fbnpp)
    deallocate(flignin)
    deallocate(fcna)
    deallocate(fcnb)
    deallocate(fid)
    deallocate(fpft)

   end subroutine getdata_soc

   subroutine get14catm(nz,f14cz,micparam)
   ! get the atmospheric 14C data 1941-2019 (inclusive, Hua et al. 2020)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(INout)   :: micparam
    integer :: i, nz, ny, nc14atm(100,5)
    real(r_2)  :: year,c14del,sdx1,c14fm,sdx2
    character(len=80) :: f14cz
    ! give 14C zones globally
    ! 14C zone        region code
    ! NH zone 1       11
    ! NH zone 2       12
    ! NH zone 3       13
    ! SH zone 1,2     14
    ! SH zone 3       15

      micparam%c14atm(:,nz,:) = 0.0
      open(13,file=f14cz)
      do i=1,4
          read(13,*)
      enddo

      do i=1,79 !! 1941-2019
        read(13,*,end=91) year,c14del,sdx1,c14fm,sdx2
        ny = year - 1940
         if(ny<1 .or. ny>79) then
            print *, 'year', year, 'outside the range'
            stop
         else
            micparam%c14atm(ny,nz,1) = c14del !!! delta c14
            micparam%c14atm(ny,nz,2) = c14fm
         endif
      enddo
91    close(13)
   end subroutine get14catm

   subroutine getdata_hwsd_dim(fhwsdsoc,mpx,timex)
    use netcdf
    use mic_constant
    use mic_variable
    implicit none
    character(len=120) :: fhwsdsoc
    integer :: mpx,timex
    integer:: ncid,varid,status
   ! open .nc file
    status = nf90_open(fhwsdsoc,nf90_nowrite,ncid)
    if(status /= nf90_noerr) print*, 'Error opening c_fraction.nc'

    ! get dimensions/profile_id
    status = nf90_inq_dimid(ncid,'nsite',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring dimensions/nsite'
    status = nf90_inquire_dimension(ncid,varid,len=mpx)
    if(status /= nf90_noerr) print*,'Error reading profile_id'

    !
    status = nf90_inq_dimid(ncid,'time',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring dimensions/ntime'
    status = nf90_inquire_dimension(ncid,varid,len=timex)
    if(status /= nf90_noerr) print*,'Error reading profile_id'

    ! Close netcdf file
    status = NF90_CLOSE(ncid)
   end subroutine  getdata_hwsd_dim

   subroutine getdata_hwsd(fhwsdsoc,jglobal,bgcopt,jopt,jmodel,micparam,micglobal)
    !use micglobal%area (area fraction) as a switch to run for selected sites during parameter optimization (jopt==0)
    !model only runs for those sites with micglobal%area(np) > 0.0
    use netcdf
    use mic_constant
    use mic_variable
    implicit none
    character(len=120) :: fhwsdsoc
    integer :: jglobal,bgcopt,jopt,jmodel
    TYPE(mic_parameter),          INTENT(INout) :: micparam
    TYPE(mic_global_input),       INTENT(INout) :: micglobal
    ! local variables
    integer:: ncid,varid,status
    integer:: np,ns,k,ipft,nsocobs
    integer:: intval,msite,isite
    integer,           dimension(:),     allocatable     :: ivarx1
    real,              dimension(:),     allocatable     :: varx1float
    double precision,  dimension(:),     allocatable     :: varx1db,avgts,avgms
    double precision,  dimension(:,:),   allocatable     :: varx2db,fsoc
    double precision,  dimension(:,:,:), allocatable     :: varx3db

    !
    ! CABLE PFT-dependent parameter values
    real(r_2), dimension(17)         :: cnleaf1,cnroot1,cnwood1,ligleaf1,ligroot1,ligwood1
    data cnLeaf1/99.60,46.20,118.60,62.80,75.20,69.60,88.00,98.40,43.20,50.00,99.60,46.20,62.80,100.00,80.00,80.00,80.00/
    data cnwood1/250.63,142.00,256.63,164.42,149.58,157.89,157.89,155.05,157.89,131.58,250.63,142.00,164.42,157.89,157.89,142.11,157.89/
    data cnroot1/81.89,68.00,83.33,70.22,74.56,71.67,69.67,76.67,67.44,78.89,81.89,68.00,70.22,78.89,78.89,78.89,78.89/
    data ligleaf1/0.25,0.20,0.20,0.20,0.20,0.10,0.10,0.10,0.10,0.10,0.25,0.20,0.20,0.15,0.15,0.25,0.10/
    data ligwood1/0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40,0.40/
    data ligroot1/0.25,0.20,0.20,0.20,0.20,0.10,0.10,0.10,0.10,0.10,0.25,0.20,0.20,0.15,0.15,0.25,0.10/

    ! ORCHIDEE PFT-dependent parameter values
    real(r_2), dimension(18)         :: cnleaf2,cnroot2,cnwood2,ligleaf2,ligroot2,ligwood2
    data cnleaf2/47.12,51.42,94.78,49.99,54.55,94.78,54.55,94.78,75.29,101.21,75.29,93.80,75.29,101.21,75.29,101.21,75.29,101.21/
    data cnwood2/115.74,113.64,208.96,125.00,145.45,208.96,145.45,36.36,16.00,22.00,16.00,16.00,16.00,22.00,16.00,22.00,16.00,22.00/
    data cnroot2/81.43,82.53,151.87,87.49,100.00,151.87,100.00,65.57,45.64,61.60,45.64,54.90,45.64,61.60,45.64,61.60,45.64,61.60/
    data ligleaf2/0.20,0.20,0.25,0.20,0.20,0.25,0.20,0.25,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10/
    data ligwood2/0.25,0.25,0.30,0.25,0.25,0.30,0.25,0.30,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10/
    data ligroot2/0.20,0.20,0.25,0.20,0.20,0.25,0.20,0.25,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10/


    allocate(ivarx1(mp))
    allocate(varx1float(mp),varx1db(mp),avgts(mp),avgms(mp))
    allocate(varx2db(mp,ntime),fsoc(mp,ms))
    allocate(varx3db(mp,ms,ntime))


   ! open .nc file
    print *, ' calling getdata_hwsd'
    print *,'input file', fhwsdsoc
    print *,'mp ms bgcopt=',    mp,ms,bgcopt

    status = nf90_open(fhwsdsoc,nf90_nowrite,ncid)
    if(status /= nf90_noerr) print*, 'Error opening c_fraction.nc'

    ! get variables
    status = nf90_inq_varid(ncid,'lat',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring data lat'
    status = nf90_get_var(ncid,varid,varx1db)
    if(status /= nf90_noerr) print*,'Error reading data lat'
    micglobal%lat = real(varx1db,kind=r_2)

    status = nf90_inq_varid(ncid,'lon',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring data lont'
    status = nf90_get_var(ncid,varid,varx1db)
    if(status /= nf90_noerr) print*,'Error reading data lon'
    micglobal%lon=real(varx1db,kind=r_2)

    status = nf90_inq_varid(ncid,'max_PFT',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring data PFT'
    status = nf90_get_var(ncid,varid,ivarx1)
    if(status /= nf90_noerr) print*,'Error reading data PFT'
    micglobal%pft = ivarx1

    status = nf90_inq_varid(ncid,'USDA_SoilSuborder',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil order'
    status = nf90_get_var(ncid,varid,ivarx1)
    if(status /= nf90_noerr) print*,'Error reading soil order'
    micglobal%sorder = ivarx1

    if(jmodel==1) then
       status = nf90_inq_varid(ncid,'isoil',varid)
       if(status /= nf90_noerr) print*, 'Error inquiring soil texturep'
       status = nf90_get_var(ncid,varid,ivarx1)
       if(status /= nf90_noerr) print*,'Error reading soil texure'
       micglobal%isoil = ivarx1
    endif
    if(jmodel==2) then
       status = nf90_inq_varid(ncid,'USDA_Soil_texture_class',varid)
       if(status /= nf90_noerr) print*, 'Error inquiring soil texturep'
       status = nf90_get_var(ncid,varid,ivarx1)
       if(status /= nf90_noerr) print*,'Error reading soil texure'
       micglobal%isoil = ivarx1
    endif

    status = nf90_inq_varid(ncid,'max_PFTfrac',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring max_PFTfrac'
    status = nf90_get_var(ncid,varid,varx1float)
    if(status /= nf90_noerr) print*,'Error reading max_PFTfrac'
    micglobal%area = real(varx1float,kind=r_2)

    status = nf90_inq_varid(ncid,'npp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring npp'
    status = nf90_get_var(ncid,varid,varx1db)
    if(status /= nf90_noerr) print*,'Error reading npp'
    micglobal%npp = real(varx1db,kind=r_2)

    status = nf90_inq_varid(ncid,'pH',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring ph'
    status = nf90_get_var(ncid,varid,varx1float)
    if(status /= nf90_noerr) print*,'Error reading ph'
    micglobal%ph = real(varx1float,kind=r_2)

    status = nf90_inq_varid(ncid,'clay',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring clay'
    status = nf90_get_var(ncid,varid,varx1float)
    if(status /= nf90_noerr) print*,'Error reading clay'
    micglobal%clay = real(varx1float,kind=r_2)

    status = nf90_inq_varid(ncid,'silt',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring silt'
    status = nf90_get_var(ncid,varid,varx1float)
    if(status /= nf90_noerr) print*,'Error reading silt'
    micglobal%silt = real(varx1float,kind=r_2)

    status = nf90_inq_varid(ncid,'bulk_density',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil bulk density'
    status = nf90_get_var(ncid,varid,varx1float)
    if(status /= nf90_noerr) print*,'Error reading bulk density'
    micglobal%bulkd= real(varx1float,kind=r_2)

    status = nf90_inq_varid(ncid,'HWSD_SOC',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil carbon'
    status = nf90_get_var(ncid,varid,fsoc)
    if(status /= nf90_noerr) print*,'Error reading soil carbon'

    status = nf90_inq_varid(ncid,'SoilTemp',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil temperature'
    status = nf90_get_var(ncid,varid,varx3db)
    if(status /= nf90_noerr) print*,'Error reading soil temperature'
    micglobal%tsoil=real(varx3db,kind=r_2)

    status = nf90_inq_varid(ncid,'SoilMoist',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil moisture'
    status = nf90_get_var(ncid,varid,varx3db)
    if(status /= nf90_noerr) print*,'Error reading soil moisture'
    micglobal%moist=real(varx3db,kind=r_2)

    status = nf90_inq_varid(ncid,'water_potential',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring soil matric potential'
    status = nf90_get_var(ncid,varid,varx3db)
    if(status /= nf90_noerr) print*,'Error reading soil matric potential'
    micglobal%matpot=real(varx3db,kind=r_2)

    status = nf90_inq_varid(ncid,'Leaf_fall',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring Leaf_fall'
    status = nf90_get_var(ncid,varid,varx2db)
    if(status /= nf90_noerr) print*,'Error reading Leaf_fall'
    micglobal%dleaf=real(varx2db,kind=r_2)

    status = nf90_inq_varid(ncid,'Belowground_litter_fall',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring Belowground_litter_fall'
    status = nf90_get_var(ncid,varid,varx2db)
    if(status /= nf90_noerr) print*,'Error reading Belowground_litter_fall'
    micglobal%droot=real(varx2db,kind=r_2)

    status = nf90_inq_varid(ncid,'non_leaf_aboveground_litterfall',varid)
    if(status /= nf90_noerr) print*, 'Error inquiring non_leaf_aboveground_litterfall'
    status = nf90_get_var(ncid,varid,varx2db)
    if(status /= nf90_noerr) print*,'Error reading non_leaf_aboveground_litterfall'
    micglobal%dwood =real(varx2db,kind=r_2)

    ! Close netcdf file
    status = NF90_CLOSE(ncid)

    do k=1,ntime
       micglobal%time(k)= real(k*1.0,kind=r_2)
    enddo

    ! print *, 'PFT=', micglobal%pft
    msite = 0
    do np=1, mp
       micglobal%siteid(np)  = np
       micglobal%bgctype(np) = micglobal%sorder(np)
       micglobal%poros(np)   = 1.0 - micglobal%bulkd(np)/2650.0
       micparam%siteid(np)   = micglobal%siteid(np)
       micparam%pft(np)      = micglobal%pft(np)
       micparam%bgctype(np)  = micglobal%bgctype(np)
       micparam%isoil(np)    = micglobal%isoil(np)
       micparam%sorder(np)   = micglobal%sorder(np)
       if(jmodel==1) then      !CABLE
          ipft =  micglobal%pft(np)
          micparam%xcnleaf(np)  = cnleaf1(ipft)
          micparam%xcnroot(np)  = cnroot1(ipft)
          micparam%xcnwood(np)  = cnwood1(ipft)
          micparam%fligleaf(np) = ligleaf1(ipft)
          micparam%fligroot(np) = ligroot1(ipft)
          micparam%fligwood(np) = ligwood1(ipft)
       endif
       if(jmodel==2) then      !ORCHIDEE
          ipft =  micglobal%pft(np)
          micparam%xcnleaf(np)  = cnleaf2(ipft)
          micparam%xcnroot(np)  = cnroot2(ipft)
          micparam%xcnwood(np)  = cnwood2(ipft)
          micparam%fligleaf(np) = ligleaf2(ipft)
          micparam%fligroot(np) = ligroot2(ipft)
          micparam%fligwood(np) = ligwood2(ipft)
       endif

       nsocobs=0
       do ns=1,ms
          micparam%csoilobs(np,ns)=real(fsoc(np,ns),kind=r_2)
          if(micparam%csoilobs(np,ns) >0.0 .and. micparam%csoilobs(np,ns) < 1000.0) nsocobs = nsocobs + 1
       enddo

       ! using "micglobal%area" to filter out some sites
       micglobal%npp(np) = sum(micglobal%dleaf(np,:) + micglobal%dwood(np,:) + micglobal%droot(np,:))
       if(micglobal%npp(np)<100.0 .or. micglobal%ph(np)<3.0 .or. nsocobs==0) micglobal%area(np) = -1.0

       if(micglobal%bgctype(np) ==bgcopt .and. micglobal%area(np) >0) msite = msite + 1
    enddo    ! "np=1,mp"

    if(msite>50) then

       intval = msite/50; isite=0
       do np=1,mp
          if(micglobal%bgctype(np) == bgcopt .and.micglobal%area(np) > 0.0) then
             isite = isite +1
             if(int(isite/intval)*intval /= isite.or. isite/intval >50) micglobal%area(np) = -1.0
          endif
          if(micglobal%area(np) > 0.0 .and. micglobal%bgctype(np) == bgcopt) then
             write(*,103) isite,np, micglobal%bgctype(np), micglobal%area(np),micglobal%npp(np),micglobal%ph(np)
          endif
       enddo
    else

      isite=0
      do np=1,mp
         if(micglobal%area(np) > 0.0 .and. micglobal%bgctype(np) == bgcopt) then
            isite=isite+1
            write(*,103) isite,np,micglobal%bgctype(np),micglobal%area(np),micglobal%npp(np),micglobal%ph(np)
         endif
      enddo
      if(isite<10) print *, 'too few sites ', isite

    endif

    micglobal%avgts(:) = sum(sum(micglobal%tsoil(:,:,:),dim=3),dim=2)/real(ms*ntime)
    micglobal%avgms(:) = sum(sum(micglobal%moist(:,:,:),dim=3),dim=2)/real(ms*ntime)

       open(100,file='inputdata.txt')
       do np=1,mp
          write(100,101) micparam%siteid(np),micglobal%area(np),micparam%pft(np),micparam%isoil(np),micparam%sorder(np),micparam%bgctype(np),   &
          micglobal%npp(np),micglobal%ph(np),micglobal%clay(np)+micglobal%silt(np),micglobal%bulkd(np), &
          micglobal%avgts(np),micglobal%avgms(np),micparam%csoilobs(np,:)
       enddo
       close(100)
101 format(i5,1x,f8.4,1x,4(i3,1x),30(f10.4,1x))
103 format(' run site', 2(i6,1x),10(f10.3,1x))

    deallocate(ivarx1)
    deallocate(varx1float,varx1db)
    deallocate(varx2db,fsoc)
    deallocate(varx3db)
    deallocate(avgts,avgms)

!    print *, 'exit getdata_hwsd'

   end subroutine getdata_hwsd


SUBROUTINE profile_hwsd(jmodel,micparam)
    ! compute the mean soil C profile for each PFT
    use mic_constant
    use mic_variable
    implicit none
    integer :: jmodel
    TYPE(mic_parameter),          INTENT(INout) :: micparam
    ! local variables
    real(r_2)  :: sdepth1,sdepth2
    real(r_2), dimension(ms)  :: sdepth,y
    real(r_2), dimension(mp)  :: y0,slope,corr
    real(r_2), dimension(mpft):: avgslope,sdslope
    integer,   dimension(mpft):: nslope
    integer :: ns,p,np

        open(21,file='profile_hwsd.txt')
        sdepth1=0.0; sdepth2=0.0
        y0(:)=0.0;slope(:)=0.0;corr(:)=0.0
        avgslope(:) =0.0; sdslope(:) = 0.0; nslope(:) =0
        do ns=1,ms
           sdepth2    = sdepth2 + zse(ns)
           if(ns>1)   sdepth(ns) = 0.5 * (sdepth1+sdepth2)
           sdepth1    = sdepth2
        enddo
        do np=1,mp
           do ns=2,ms
              y(ns-1) = log(micparam%csoilobs(np,ns))
          enddo
          call linreg(ms-1,sdepth,y,y0(np),slope(np),corr(np))
          write(21,201) micparam%siteid(np),micparam%pft(np),micparam%isoil(np), &
                        micparam%sorder(np),jmodel,y0(np),slope(np),corr(np)
        enddo

        do p=1,mpft
           do np=1,mp
              if(slope(np) /= slope(np) .or. slope(np) > 0.0) then
                 print *, 'invalid profile', np, micparam%csoilobs(np,:)
              else
                 if(micparam%pft(np)==p) then
                    avgslope(p)= avgslope(p) + slope(np)
                    sdslope(p) = sdslope(p)  + slope(np)*slope(np)
                    nslope(p)  = nslope(p)   + 1
                 endif
              endif
           enddo
           if(nslope(p) >1) then
              avgslope(p)= avgslope(p)/real(nslope(p))
              sdslope(p) = (sdslope(p) - nslope(p) *avgslope(p)*avgslope(p))/real(nslope(p))
              sdslope(p) = sqrt(sdslope(p))
           else
              sdslope(p) = 0.0
           endif
           write(21,202) jmodel,p,avgslope(p),sdslope(p),nslope(p)
        enddo
       close(21)
201 format(5(i4,1xs),3(f10.4,1x))
202 format(2(i4,2x),2(f10.4,2x),i4)

END SUBROUTINE profile_hwsd

subroutine linreg(ndata,x,y,a,b,r)
   use mic_constant
   implicit none                                                                    ! no default data types
   ! input
   integer :: ndata
   real(r_2), dimension(ndata)       ::  x, y
   ! output
   real(r_2) :: a,b,r
   !local variables
   integer :: n
   real(r_2) :: sumx,sumy,sumx2, sumy2, sumxy

   sumx =0.0; sumx2=0.0; sumxy=0.0;sumy=0.0;sumy2=0.0
   a=0.0;b=0.0;r=0.0
   do n=1,ndata                                                                     ! loop for all data points
      sumx  = sumx + x(n)                                                           ! compute sum of x
      sumx2 = sumx2 + x(n) * x(n)                                                   ! compute sum of x**2
      sumxy = sumxy + x(n) * y(n)                                                   ! compute sum of x * y
      sumy  = sumy + y(n)                                                           ! compute sum of y
      sumy2 = sumy2 + y(n) * y(n)                                                   ! compute sum of y**2
   end do

   b = (ndata * sumxy  -  sumx * sumy) / (ndata * sumx2 - sumx**2)                  ! compute slope
   a = (sumy * sumx2  -  sumx * sumxy) / (ndata * sumx2  -  sumx**2)                ! compute y-intercept
   r = (sumxy - sumx * sumy /real(ndata)) /                 &
       sqrt((sumx2 - sumx**2/real(ndata)) * (sumy2 - sumy**2/real(ndata)))  ! compute correlation coefficient

   print *, 'a b r', a,b,r
end subroutine linreg

SUBROUTINE vmic_param_constant(kinetics,micpxdef,micpdef,micparam)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_param_xscale),       INTENT(IN)    :: micpxdef
    TYPE(mic_param_default),      INTENT(IN)    :: micpdef
    TYPE(mic_parameter),          INTENT(INout) :: micparam
    !local variables
    real(r_2), dimension(:,:), allocatable      :: froot
    real(r_2), dimension(:),   allocatable      :: totroot
    integer    :: nopt,npft,np,ns,kinetics
    real(r_2)  :: depths1,depths2,krootx


      allocate(froot(mpft,ms))
      allocate(totroot(mpft))

      do np=1,mp
!         print *, 'np pft', np,micparam%pft(np)
          nopt=micparam%bgctype(np)
          do ns=1,ms
            micparam%Q1(np,ns)=micpdef%Q1
            micparam%Q2(np,ns)=micpdef%Q2
            micparam%fm(np,ns)=micpdef%fm * micpxdef%xfm(nopt)
            micparam%fs(np,ns)=micpdef%fs * micpxdef%xfs(nopt)
         enddo  !ns
      enddo     ! np

      depths1=0.0;depths2=0.0
      do ns=1,ms
          depths2 = depths2 + zse(ns)
          do npft=1,mpft
              krootx = micpdef%rootbeta * micpxdef%xrootbeta(npft)
              froot(npft,ns) = (1.0/krootx) *( exp(-krootx*depths1)-exp(-krootx*depths2))
          enddo  !npft
          depths1=depths2
      enddo   !ns


      do npft=1,mpft
         totroot(npft) =sum(froot(npft,1:ms))
      enddo

      ! !normalizing
      do ns=1,ms
         do npft=1,mpft
            froot(npft,ns) = froot(npft,ns)/totroot(npft)
         enddo
      enddo

      ! calculate mp by ms all parameter values
      do np=1, mp
         npft=min(mpft,max(1,micparam%pft(np)))
         nopt=micparam%bgctype(np)
         do ns=1,ms
            micparam%sdepth(np,ns)   = zse(ns)
            micparam%fracroot(np,ns) = froot(npft,ns)
         enddo !"ns"
          micparam%diffsocx(np) = micpxdef%xdiffsoc(nopt) * micpdef%diffsoc  !"diffsoc" from mic_constant
      enddo    ! "np=1,mp"

      if(diag==1) then
         print *, micparam%fracroot(outp,:)
         print *, micparam%sdepth(outp,:)
         print *, micparam%diffsocx(outp)
      endif

      ! the following parameters are specific to kinetics3
      if(kinetics==3) then
         do np=1,mp
            do ns=1,ms
               nopt=micparam%bgctype(np)
               micparam%kadsorp(np,ns)  = micpdef%kadsorpx
               micparam%kdesorp(np,ns)  = micparam%kadsorp(np,ns) /(micpdef%kbax * micpxdef%xkba(nopt))
               micparam%fp2a(np,ns)     = micpdef%fp2ax     * micpxdef%xfp2ax(nopt)
               micparam%tvcpool(np,ns)  = micpdef%tvcpoolx  * micpxdef%xtvc(nopt)
               micparam%tvppool(np,ns)  = micpdef%tvppoolx  * micpxdef%xtvp(nopt)
               micparam%tvac(np,ns)     = micpdef%tvacx     * micpxdef%xtvac(nopt)
               micparam%qmaxcoeff(np,ns)= micpdef%qmaxcoeff * micpxdef%xqmaxcoeff(nopt)
            enddo  !ns
         enddo  !np
      endif ! kinetics


    deallocate(froot)
    deallocate(totroot)
END SUBROUTINE vmic_param_constant

subroutine vmic_param_time(kinetics,micpxdef,micpdef,micparam,micinput,micnpool)
    ! time-dependent model parameters, called every time step if the forcing, such air temperature
    ! varies every time step
    ! otherwise only called at the start the integration
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_param_xscale),       INTENT(IN)      :: micpxdef
    TYPE(mic_param_default),      INTENT(IN)      :: micpdef
    TYPE(mic_parameter),          INTENT(INout)   :: micparam
    TYPE(mic_input),              INTENT(INout)   :: micinput
    TYPE(mic_npool),              INTENT(INOUT)   :: micnpool

    integer  :: kinetics
      ! compute fractions
      call bgc_fractions(micpxdef,micpdef,micparam,micinput)
      ! compute microbial growth efficiency
      call mget(micpdef,micparam,micinput,micnpool)
      ! compute microbial turnover rates
      call turnovert(kinetics,micpxdef,micpdef,micparam,micinput)
      if(kinetics/=3) call Desorpt(micpxdef,micparam,micinput)
      call Vmaxt(micpxdef,micpdef,micparam,micinput)
      call Kmt(micpxdef,micpdef,micparam,micinput)

end subroutine vmic_param_time

subroutine vmic_init(miccpool,micnpool)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_cpool),              INTENT(INOUT)   :: miccpool
    TYPE(mic_npool),              INTENT(INOUT)   :: micnpool
    integer :: ip
    real(r_2), dimension(:), allocatable    :: cpooldef

      allocate(cpooldef(mcpool))
!	  print *, 'calling vmic_init'

      cpooldef(1) = 16.5*0.1;     cpooldef(2) = 16.5*0.1
      cpooldef(3) = 16.5*0.025;   cpooldef(4) = 16.5*0.025
      cpooldef(5) = 16.5*0.1125;  cpooldef(6) = 16.5*0.375;  cpooldef(7) = 16.5*0.2625
      cpooldef(8) = 0.0;          cpooldef(9) = 0.0;         cpooldef(10)= 0.0

      do ip=1,mcpool
         miccpool%cpool(:,:,ip) = cpooldef(ip)
      enddo
!    print *, 'at np=1 ns=1 cpool', miccpool%cpool(1,1,1:mcpool)
    deallocate(cpooldef)
end subroutine vmic_init


subroutine vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_param_xscale),  INTENT(INOUT)   :: micpxdef
    integer :: bgcopt,jmodel
    real*8, dimension(16)                    :: xopt
    real(r_2), dimension(17)                 :: xrootcable
    real(r_2), dimension(18)                 :: xrootorchidee
    data xrootcable/1.43,0.94,1.43,1.04,0.77,0.85,0.62,1.77,0.94,0.94,1.43,0.94,1.04,0.53,1.00,1.00,1.00/
    data xrootorchidee/0.94,0.94,1.04,1.04,1.04,1.43,1.43,1.43,0.85,0.62,0.94,0.94,0.85,0.85,0.85,0.85,0.85,0.85/
    integer :: i

!  1: xav:      scaling factor for V                            [1]                 (0-30)
!  2: xak:      scaling factor for K                            [1]                 (0-30)
!  3: xfp2ax:   scaling factor for fp2ax                        [0.33]              (0.5-2.0)
!  4: xfm:      scaling factor for fm                           [0.05]              (0.1-2.0)
!  5: xfs:      scaling factor for fs                           [0.05]              (0.1-2.0)
!  6: xtvmic:   scaling factor for tvmicR (0-10)                [0.00052]           (0.1,10)
!  6: xtvmic :  scaling factor for tvmicK (0-10)                [0.00024]           =xtvmicR
!  7: xtvp:     scaling factor for tvppool(0-10)                [0.019]             (0.1,10)
!  8: xtvc:     scaling factor for tvcpool(0-10)                [0.02]              (0.1,10)
!  9: xtvac:    scaling factor for tvac   (0-10)                [0.0015]            (0.1,10)
! 10: xkba:     scaling factor for kba    (0.2-5)               [8]                 (0.1,2)
! 11: xqmaxcoeff:coefficient of Qmax on clay+silt (0.4-0.8)     [0.6]               (0.5,1.5)
! 12: xbeta:     beta parameter            (1.01-2.0)           [2.0]               (1.05,2.0)
! 13: xdiffsoc:  SOC diffusion/bioturbation rate                [1.0]               (0.1,10.0)
! 14: xNPP:      carbon input                                   [1.0]               (0.5,2.0)
! 15: xdesorp:   desorption coefficient (kinetics=1,2)          [1.0]               (0.1,10.0)

     ! assign the default values
     ! this shoudl be replaced by a parameter lookup tables for gloabl simulations
      micpxdef%xav       = 1.0
      micpxdef%xak       = 1.0
      micpxdef%xfm       = 1.0
      micpxdef%xfs       = 1.0
      micpxdef%xtvmic    = 1.0
      micpxdef%xtvp      = 1.0
      micpxdef%xtvc      = 1.0   ! unstop the backflow from MAOC to aggregate pool
      micpxdef%xtvac     = 1.0
      micpxdef%xkba      = 1.0
      micpxdef%xqmaxcoeff= 1.0
      micpxdef%xdiffsoc  = 1.0
      micpxdef%xNPP      = 1.0
!      micpxdef%xrootbeta = 1.0
      micpxdef%xvmaxbeta = 1.0

      micpxdef%xfp2ax    = 1.0
      micpxdef%xbeta     = 1.0
      micpxdef%xdesorp   = 1.0

      do i=1,mpft
         if(jmodel==1) then
            micpxdef%xrootbeta(i) = xrootcable(i)
         endif
         if(jmodel==2) then
            micpxdef%xrootbeta(i) = xrootorchidee(i)
         endif
      enddo

      ! assign the values to the optimized parameters

      micpxdef%xav(bgcopt)        = xopt(1)
      micpxdef%xak(bgcopt)        = xopt(2)
      micpxdef%xfm(bgcopt)        = xopt(3)
      micpxdef%xfs(bgcopt)        = xopt(4)
      micpxdef%xtvmic(bgcopt)    = xopt(5)
      micpxdef%xtvp(bgcopt)       = xopt(6)
      micpxdef%xtvc(bgcopt)       = xopt(7)
      micpxdef%xtvac(bgcopt)      = xopt(8)
      micpxdef%xkba(bgcopt)       = xopt(9)
      micpxdef%xqmaxcoeff(bgcopt) = xopt(10)
      micpxdef%xdiffsoc(bgcopt)   = xopt(11)
      micpxdef%xnpp(:)            = xopt(12)
!      micpxdef%xrootbeta(pftopt)  = xopt(13)
      micpxdef%xvmaxbeta(bgcopt)  = xopt(14)

end subroutine vmic_param_xscale


subroutine variable_time(year,doy,micglobal,micinput,micnpool)
    use mic_constant
    use mic_variable
    implicit none
    integer :: year,doy
    TYPE(mic_global_input),  INTENT(INout)   :: micglobal
    TYPE(mic_input),         INTENT(INout)   :: micinput
    TYPE(mic_npool),         INTENT(INOUT)   :: micnpool
    integer :: np,ns


!        print *, 'calling global2np- ntime', ntime
     do np=1,mp
        micinput%fcnpp(np)      = max(0.0,micglobal%npp(np))              !gc/m2/year
        micinput%Dleaf(np)      = (micglobal%dleaf(np,doy)/24.0)*delt     !gc/m2/delt
        micinput%Droot(np)      = (micglobal%droot(np,doy)/24.0)*delt     !gc/m2/delt
        micinput%Dwood(np)      = (micglobal%dwood(np,doy)/24.0)*delt     !gc/m2/delt

        do ns=1,ms
           micinput%tavg(np,ns)     = micglobal%tsoil(np,ns,doy)  ! average temperature in deg C
           micinput%wavg(np,ns)     = micglobal%moist(np,ns,doy)  ! average soil water content mm3/mm3
           micinput%matpot(np,ns)   = micglobal%matpot(np,ns,doy)
           micinput%clay(np,ns)     = micglobal%clay(np)     ! clay content (fraction)
           micinput%silt(np,ns)     = micglobal%silt(np)     ! silt content (fraction)
           micinput%ph(np,ns)       = micglobal%ph(np)
           micinput%porosity(np,ns) = micglobal%poros(np)    ! porosity mm3/mm3
           micinput%bulkd(np,ns)    = micglobal%bulkd(np)
           micnpool%mineralN(np,ns) = 0.1                    ! g N /kg soil
        enddo !"ns"
     enddo !"np"

end subroutine variable_time

    subroutine vmicsoil_c14(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,ifsoc14,bgcopt,nyeqpool, &
                    micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)
    use mic_constant
    use mic_variable
   !  use omp_lib
    implicit none
    TYPE(mic_param_xscale),  INTENT(INOUT)   :: micpxdef
    TYPE(mic_param_default), INTENT(IN)      :: micpdef
    TYPE(mic_parameter),     INTENT(INout)   :: micparam
    TYPE(mic_input),         INTENT(INout)   :: micinput
    TYPE(mic_global_input),  INTENT(INout)   :: micglobal
    TYPE(mic_cpool),         INTENT(INOUT)   :: miccpool
    TYPE(mic_npool),         INTENT(INOUT)   :: micnpool
    TYPE(mic_output),        INTENT(INout)   :: micoutput

    integer :: ifsoc14,isoc14,kinetics,bgcopt,nyeqpool

    ! local variables
    real(r_2),    dimension(mcpool)    :: xpool0,xpool1
    real(r_2),    dimension(ms)        :: ypooli,ypoole,fluxsoc
    real(r_2),    dimension(ms)        :: cfluxa

!    integer       ndelt,n1,n2,i,j,year,ip,np,ns,ny,nyrun
    integer       :: ndelt,i,j,year,ip,np,ns,ny,nyrun
    real(r_2)     :: timex,delty,fluxdocsx,diffsocxx

    integer    :: jrestart
    character(len=99) :: frestart_in,frestart_out,foutput
    real(r_2)  :: cpool0, cpool1, totcinput

   ! local variables
   real(r_2)                      :: deltD !,tot0,tot1,totflux
   real(r_2), dimension(:), allocatable    :: xzse
   real(r_2), dimension(:), allocatable    :: sdepthx
   real(r_2)                                :: coeffA, coeffB
   real(r_2), dimension(:), allocatable    :: at,bt,ct,rt
   real(r_2), dimension(:), allocatable    :: xpool
   real(r_2)  :: cleachloss


    allocate(xzse(ms))
    allocate(sdepthx(ms+1))
    allocate(at(ms),bt(ms),ct(ms),rt(ms))
    allocate(xpool(ms))

      call vmic_param_constant(kinetics,micpxdef,micpdef,micparam)
      call vmic_init(miccpool,micnpool)
      call vmic_param_time(kinetics,micpxdef,micpdef,micparam,micinput,micnpool)

    !  print *, 'initial pool size np=1 ns=1', miccpool%cpool(1,1,:)
    !  print *, 'xav=', bgcopt,micpxdef%xav(:)

      if(jrestart==1) call vmic_restart_read(miccpool,micnpool,frestart_in)

      ndelt   = int(24*365/delt) ! number of time step per year in "delt" unit

! "data copyin" : all data accessed by GPU
! "create": intermediate variables used by GPU
! "copyout" copy data out from GPU
! "private": every cell-dependent variables used in paralelling computing

!$acc data copyin(micpdef,micparam,miccpool,micinput,micoutput,  &
!$acc nyrun,bgcopt,ndelt,zse,kinetics,nyeqpool,isoc14)      &
!$acc create(delty,timex,fluxsoc,year,ny,i,ns,np,ip,xpool0,xpool1,ypooli,ypoole,diffsocxx,cfluxa, &
!$acc j,deltD,xzse,sdepthx,coeffA,coeffB,at,bt,ct,rt,xpool,cleachloss)              &
!$acc copyout(miccpool%cpool,micoutput%fluxcinput,micoutput%fluxrsoil,micoutput%fluxcleach)
!$acc PARALLEL LOOP                                                      &
!$acc private(delty,timex,fluxsoc,year,ny,i,ns,ip,np,xpool0,xpool1,ypooli,ypoole,diffsocxx,cfluxa,&
!$acc j,deltD,xzse,sdepthx,coeffA,coeffB,at,bt,ct,rt,xpool,cleachloss)

      do np=1,mp

      if(micparam%bgctype(np)==bgcopt) then
         if (ifsoc14 == 1) then
             nyrun = micparam%nyc14obs(np) - 1940 + nyeqpool !! how many years to run to get equilibrium
         else
             nyrun = nyeqpool
         endif

    !     print *,'np pft npp anpp bnpp = ',np,micparam%pft(np),micinput%fcnpp(np), micinput%dleaf(np)*365.0*24.0, micinput%droot(np)*365.0*24.0

         do year=1,nyrun
            ny = year-nyrun
            micoutput%fluxcinput(np)=0.0; micoutput%fluxrsoil(np) = 0.0; micoutput%fluxcleach(np)= 0.0    ! yearly fluxes

            do i=1,365

               ! for each soil layer
               ! sum last all C pools of all layers for compute the soil respiration = input - sum(delCpool)
               ! before leaching is computed
                cpool0 =0.0; cpool1 =0.0; totcinput = 0.0
               do ns=1,ms
                 ! micinput%cinputm(np,ns)+micinput%cinputs(np,ns) in mg C/cm3/delt
                  totcinput =totcinput + (micinput%cinputm(np,ns)+micinput%cinputs(np,ns)) *1000.0 * zse(ns)   ! convert to g C/m2/delt/zse

                  do ip=1,mcpool
                     xpool0(ip) = miccpool%cpool(np,ns,ip)
                     cpool0     = cpool0  + xpool0(ip) * zse(ns) * 1000.0  ! 1000 for mg C/cm3 to g C/m2/zse
                  enddo

                 ! here the integration step is "delty" in rk4 and "ndelt" is number of "delt (hour) per year
                  timex=real(i*delt)
                  delty = real(ndelt)/(365.0*delt)  ! time step in rk4 in "24 * delt (or daily)", all C input are in " per delt"
                  call rk4modelx(timex,delty,ny,isoc14,np,ns,kinetics,micpdef,micparam,micinput,xpool0,xpool1)

                  do ip=1,mcpool
                     miccpool%cpool(np,ns,ip) = max(xpool1(ip),1.0e-8)
                     cpool1 = cpool1 + miccpool%cpool(np,ns,ip) * zse(ns) * 1000.0  ! 1000 for mg C/cm3 to g C/m2/zse
                  enddo

               enddo    ! "ns"

               micoutput%fluxcinput(np)= micoutput%fluxcinput(np) + totcinput * real(delty)
               micoutput%fluxrsoil(np) = micoutput%fluxrsoil(np)  + totcinput * real(delty) + (cpool1 - cpool0)

                if(diag==1) then
                   print *, 'year day site np1', year, i, outp,micparam%diffsocx(outp)
                    do ns=1,ms
                       print *, ns, miccpool%cpool(outp,ns,:)
                    enddo
                endif

               do ip=1,mcpool
                  do ns=1,ms
                     ypooli(ns) = miccpool%cpool(np,ns,ip)      ! in mg c/cm3
                  enddo  !"ns"

                  fluxsoc(:) = 0.0  ! This flux is added in "modelx"
                  diffsocxx= micparam%diffsocx(np)


                 !Move bioturb here to work around memory auto allocation failure
                 !call bioturb(int(delty/delty),ms,zse,delty,diffsocxx,fluxsoc,ypooli,ypoole)  ! only do every 24*delt
                 !subroutine bioturb(ndelt,ms,zse,delt,diffsocxx,fluxsoc,xpooli,xpoole)

                  sdepthx(1) = 0.0          ! depth of a layer from the top (x_0.5=0.0 eg soil surface)
                  do j=2,ms+1
                     sdepthx(j) = sdepthx(j-1) + zse(j-1)*100.0     ! depth of the bottom of each layer (eg x_j+0.5)
                                                                    !*100 to convert from m to cm
                  enddo

                  do j=1,ms
                     xzse(j) = 0.5 * (sdepthx(j) + sdepthx(j+1))    ! depth of midpoint of a layer j  (x_j)
                  enddo

                  deltD = diffsocxx * delty
                  xpool = ypooli

                  !do i=1,1 ( int(delty/delty) == 1 )
                     do j=1,ms
                        if(j==1) then
                           coeffB = 1.0/(sdepthx(2)-sdepthx(1))
                           coeffA = deltD*coeffB/(xzse(2)-xzse(1))
                           ! Crank-Nicholson
                           at(1) = 0.0
                           bt(1) = 1.0 + 0.5 * coeffA
                           ct(1) =     - 0.5 * coeffA
                           rt(1) = (1.0-0.5*coeffA) * xpool(1) + 0.5 * coeffA * xpool(2) &
                                 +  fluxsoc(1) * delt
                        endif
                        if(j>1.and.j<ms) then
                          coeffA = deltD/((xzse(j+1)-xzse(j))*(sdepthx(j+1)-sdepthx(j)))
                          coeffB = (xzse(j+1)-xzse(j))/(xzse(j)-xzse(j-1))
                          ! Crank-Nicholson
                          at(j) =    -0.5 * coeffA * coeffB
                          bt(j) = 1.0+0.5 * coeffA *(1.0+coeffB)
                          ct(j) =    -0.5 * coeffA
                          rt(j) = 0.5 * coeffA * coeffB * xpool(j-1)        &
                                  +(1.0-0.5* coeffA*(1.0+coeffB))*xpool(j)  &
                                  + 0.5* coeffA * xpool(j+1)                &
                                  + fluxsoc(j) *delt
                        endif
                        if(j==ms) then
                            coeffA = deltD/((xzse(ms)-xzse(ms-1))*(sdepthx(ms+1) - sdepthx(ms)))
                          ! Crank-Nicholson
                            at(ms) = -0.5 * coeffA
                            bt(ms) = 1.0 + 0.5 * coeffA
                            ct(ms) = 0.0
                            rt(ms) = 0.5* coeffA  * xpool(ms-1) + (1.0-0.5 * coeffA) * xpool(ms) &
                                   + fluxsoc(ms) * delt
                        endif
                     enddo
                     call tridag(at,bt,ct,rt,xpool,ms)
                  !enddo
                  ypoole = xpool

                  !!! end bioturb

                  do ns=1,ms
                     miccpool%cpool(np,ns,ip) = ypoole(ns)
                  enddo
               enddo ! "ip=1,mcpool"

               ! computing daily leaching loss from bottom-layer LWMC
               cleachloss = micparam%tvac(np,ms) * sqrt(micinput%wavg(np,ms)/micinput%porosity(np,ms)) *  miccpool%cpool(np,ms,7)  * 24.0
               cleachloss = max(0.0,min(cleachloss,miccpool%cpool(np,ms,7)))
               micoutput%fluxcleach(np) = micoutput%fluxcleach(np) + cleachloss
               miccpool%cpool(np,ms,7)  = miccpool%cpool(np,ms,7)  - cleachloss

            enddo   ! "day"

         enddo !"year"

      endif   !bgctype(np)=bgcopt
   enddo !"mp"

!$ACC END PARALLEL
!$ACC END DATA


    miccpool%cpooleq(:,:,:) = miccpool%cpool(:,:,:)

   !  call vmic_output_write(foutput,micinput,micoutput)
   !  call vmic_restart_write(frestart_out,miccpool,micnpool)

    deallocate(xzse)
    deallocate(sdepthx)
    deallocate(at,bt,ct,rt)
    deallocate(xpool)

    end subroutine vmicsoil_c14

    subroutine vmicsoil_soc(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,bgcopt,nyeqpool, &
                    micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)
    use mic_constant
    use mic_variable

    implicit none
    TYPE(mic_param_xscale),  INTENT(INOUT)   :: micpxdef
    TYPE(mic_param_default), INTENT(IN)      :: micpdef
    TYPE(mic_parameter),     INTENT(INout)   :: micparam
    TYPE(mic_input),         INTENT(INout)   :: micinput
    TYPE(mic_global_input),  INTENT(INout)   :: micglobal
    TYPE(mic_cpool),         INTENT(INOUT)   :: miccpool
    TYPE(mic_npool),         INTENT(INOUT)   :: micnpool
    TYPE(mic_output),        INTENT(INout)   :: micoutput

    integer :: isoc14,kinetics,bgcopt,nyeqpool
    ! local variables
    real(r_2),    dimension(mcpool)    :: xpool0,xpool1
    real(r_2),    dimension(ms)        :: ypooli,ypoole,fluxsoc
    real(r_2),    dimension(ms)        :: cfluxa

    integer       :: ndelt,i,j,year,ip,np,ns,ny
    real(r_2)     :: timex,delty,fluxdocsx,diffsocxx

    integer    :: jrestart
    character(len=99) :: frestart_in,frestart_out,foutput
    real(r_2)  :: cpool0, cpool1, totcinput

       ! local variables
   ! for numerical solution
   real(r_2),    parameter            :: tol = 1.0E-04
   real(r_2),    parameter            :: tolx = 0.0001
   real(r_2),    parameter            :: tolf = 0.000001
   integer,      parameter            :: ntrial = 100

   ! local variables
   real(r_2)                      :: deltD !,tot0,tot1,totflux
   real(r_2), dimension(:), allocatable    :: xzse
   real(r_2), dimension(:), allocatable    :: sdepthx
   real(r_2)                      :: coeffA, coeffB
   real(r_2), dimension(:), allocatable    :: at,bt,ct,rt
   real(r_2), dimension(:), allocatable    :: xpool
   real(r_2)  :: cleachloss


    allocate(xzse(ms))
    allocate(sdepthx(ms+1))
    allocate(at(ms),bt(ms),ct(ms),rt(ms))
    allocate(xpool(ms))

      call vmic_param_constant(kinetics,micpxdef,micpdef,micparam)
      call vmic_init(miccpool,micnpool)
      call vmic_param_time(kinetics,micpxdef,micpdef,micparam,micinput,micnpool)

      if(jrestart==1) call vmic_restart_read(miccpool,micnpool,frestart_in)

      ndelt   = int(24*365/delt) ! number of time step per year in "delt" unit

! "data copyin" : all data accessed by GPU
! "create": intermediate variables used by GPU
! "copyout" copy data out from GPU
! "private": every cell-dependent variables used in paralelling computing

!$acc data copyin(micpdef,micparam,miccpool,micinput,micoutput,  &
!$acc bgcopt,ndelt,zse,kinetics,nyeqpool,isoc14)      &
!$acc create(delty,timex,fluxsoc,year,ny,i,ns,np,ip,xpool0,xpool1,ypooli,ypoole,diffsocxx,cfluxa, &
!$acc j,deltD,xzse,sdepthx,coeffA,coeffB,at,bt,ct,rt,xpool,cleachloss)              &
!$acc copyout(miccpool%cpool,micoutput%fluxcinput,micoutput%fluxrsoil,micoutput%fluxcleach)
!$acc PARALLEL LOOP                                                      &
!$acc private(delty,timex,fluxsoc,year,ny,i,ns,ip,np,xpool0,xpool1,ypooli,ypoole,diffsocxx,cfluxa,&
!$acc j,deltD,xzse,sdepthx,coeffA,coeffB,at,bt,ct,rt,xpool,cleachloss)

      do np=1,mp

      if(micparam%bgctype(np)==bgcopt) then

         do year=1,nyeqpool
            ny = year-nyeqpool
            micoutput%fluxcinput(np)=0.0; micoutput%fluxrsoil(np) = 0.0; micoutput%fluxcleach(np)= 0.0    ! yearly fluxes

            do i=1,365

               ! for each soil layer
               ! sum last all C pools of all layers for compute the soil respiration = input - sum(delCpool)
               ! before leaching is computed
                cpool0 =0.0; cpool1 =0.0; totcinput = 0.0
               do ns=1,ms
                 ! micinput%cinputm(np,ns)+micinput%cinputs(np,ns) in mg C/cm3/delt
                  totcinput =totcinput + (micinput%cinputm(np,ns)+micinput%cinputs(np,ns)) *1000.0 * zse(ns)   ! convert to g C/m2/delt/zse

                  do ip=1,mcpool
                     xpool0(ip) = miccpool%cpool(np,ns,ip)
                     cpool0     = cpool0  + xpool0(ip) * zse(ns) * 1000.0  ! 1000 for mg C/cm3 to g C/m2/zse
                  enddo

                 ! here the integration step is "delty" in rk4 and "ndelt" is number of "delt (hour) per year
                  timex=real(i*delt)
                  delty = real(ndelt)/(365.0*delt)  ! time step in rk4 in "24 * delt (or daily)", all C input are in " per delt"

                  call rk4modelx(timex,delty,ny,isoc14,np,ns,kinetics,micpdef,micparam,micinput,xpool0,xpool1)

                  do ip=1,mcpool
                     miccpool%cpool(np,ns,ip) = max(xpool1(ip),1.0e-8)
                     cpool1 = cpool1 + miccpool%cpool(np,ns,ip) * zse(ns) * 1000.0  ! 1000 for mg C/cm3 to g C/m2/zse
                  enddo

               enddo    ! "ns"

               micoutput%fluxcinput(np)= micoutput%fluxcinput(np) + totcinput * real(delty)
               micoutput%fluxrsoil(np) = micoutput%fluxrsoil(np)  + totcinput * real(delty) + (cpool1 - cpool0)

               ! do labile carbon leaching only for kinetics=3
               ! the following leachate transport calculations caused mass imbalance: disabled temporarily
            !   if(kinetics==3) then
            !     cfluxa(:)=0.0
            !     do ns=1,ms
            !        cfluxa(ns) = sqrt(micinput%wavg(np,ns)/micinput%porosity(np,ns)) * micparam%tvac(np,ns) * miccpool%cpool(np,ns,7) * delty
            !        if(ns==1) then
            !           miccpool%cpool(np,ns,7) = miccpool%cpool(np,ns,7) - cfluxa(ns)
            !        else
            !           miccpool%cpool(np,ns,7) = miccpool%cpool(np,ns,7) + cfluxa(ns-1) -cfluxa(ns)
            !        endif
            !     enddo
            !    ! converting flux from mg C cm-3 delty-1 to g C m-2 delty-1
            !    micoutput%fluxcleach(np) = micoutput%fluxcleach(np) + cfluxa(ms) * zse(ms) * 1000.0
            !  endif

                if(diag==1) then
                   print *, 'year day site np1', year, i, outp,micparam%diffsocx(outp)
                    do ns=1,ms
                       print *, ns, miccpool%cpool(outp,ns,:)
                    enddo
                endif

               do ip=1,mcpool
                  do ns=1,ms
                     ypooli(ns) = miccpool%cpool(np,ns,ip)      ! in mg c/cm3
                  enddo  !"ns"

                  fluxsoc(:) = 0.0  ! This flux is added in "modelx"
                  diffsocxx= micparam%diffsocx(np)

                 !Move bioturb here to work around memory auto allocation failure
                 !call bioturb(int(delty/delty),ms,zse,delty,diffsocxx,fluxsoc,ypooli,ypoole)  ! only do every 24*delt
                 !subroutine bioturb(ndelt,ms,zse,delt,diffsocxx,fluxsoc,xpooli,xpoole)

                  sdepthx(1) = 0.0          ! depth of a layer from the top (x_0.5=0.0 eg soil surface)
                  do j=2,ms+1
                     sdepthx(j) = sdepthx(j-1) + zse(j-1)*100.0     ! depth of the bottom of each layer (eg x_j+0.5)
                                                                    !*100 to convert from m to cm
                  enddo

                  do j=1,ms
                     xzse(j) = 0.5 * (sdepthx(j) + sdepthx(j+1))    ! depth of midpoint of a layer j  (x_j)
                  enddo

                  deltD = diffsocxx * delty
                  xpool = ypooli

                  !do i=1,1 ( int(delty/delty) == 1 )
                     do j=1,ms
                        if(j==1) then
                           coeffB = 1.0/(sdepthx(2)-sdepthx(1))
                           coeffA = deltD*coeffB/(xzse(2)-xzse(1))
                           ! Crank-Nicholson
                           at(1) = 0.0
                           bt(1) = 1.0 + 0.5 * coeffA
                           ct(1) =     - 0.5 * coeffA
                           rt(1) = (1.0-0.5*coeffA) * xpool(1) + 0.5 * coeffA * xpool(2) &
                                 +  fluxsoc(1) * delt
                        endif
                        if(j>1.and.j<ms) then
                          coeffA = deltD/((xzse(j+1)-xzse(j))*(sdepthx(j+1)-sdepthx(j)))
                          coeffB = (xzse(j+1)-xzse(j))/(xzse(j)-xzse(j-1))
                          ! Crank-Nicholson
                          at(j) =    -0.5 * coeffA * coeffB
                          bt(j) = 1.0+0.5 * coeffA *(1.0+coeffB)
                          ct(j) =    -0.5 * coeffA
                          rt(j) = 0.5 * coeffA * coeffB * xpool(j-1)        &
                                  +(1.0-0.5* coeffA*(1.0+coeffB))*xpool(j)  &
                                  + 0.5* coeffA * xpool(j+1)                &
                                  + fluxsoc(j) *delt
                        endif
                        if(j==ms) then
                            coeffA = deltD/((xzse(ms)-xzse(ms-1))*(sdepthx(ms+1) - sdepthx(ms)))
                          ! Crank-Nicholson
                            at(ms) = -0.5 * coeffA
                            bt(ms) = 1.0 + 0.5 * coeffA
                            ct(ms) = 0.0
                            rt(ms) = 0.5* coeffA  * xpool(ms-1) + (1.0-0.5 * coeffA) * xpool(ms) &
                                   + fluxsoc(ms) * delt
                        endif
                     enddo
                     call tridag(at,bt,ct,rt,xpool,ms)
                  !enddo
                  ypoole = xpool

                  !!! end bioturb

                  do ns=1,ms
                     miccpool%cpool(np,ns,ip) = ypoole(ns)
                  enddo
               enddo ! "ip=1,mcpool"

               ! computing daily leaching loss from bottom-layer LWMC
               cleachloss = micparam%tvac(np,ms) * sqrt(micinput%wavg(np,ms)/micinput%porosity(np,ms)) *  miccpool%cpool(np,ms,7)  * 24.0
               cleachloss = max(0.0,min(cleachloss,miccpool%cpool(np,ms,7)))
               micoutput%fluxcleach(np) = micoutput%fluxcleach(np) + cleachloss
               miccpool%cpool(np,ms,7)  = miccpool%cpool(np,ms,7)  - cleachloss

            enddo   ! "day"
         enddo !"year"
      endif   !bgctype(np)=bgcopt
   enddo !"mp"
!$ACC END PARALLEL
!$ACC END DATA

    miccpool%cpooleq(:,:,:) = miccpool%cpool(:,:,:)

   !  call vmic_output_write(foutput,micinput,micoutput)
   !  call vmic_restart_write(frestart_out,miccpool,micnpool)

    deallocate(xzse)
    deallocate(sdepthx)
    deallocate(at,bt,ct,rt)
    deallocate(xpool)

    end subroutine vmicsoil_soc

subroutine vmicsoil_hwsd_cpu(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,bgcopt,nyeqpool, &
                    micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)
    use mic_constant
    use mic_variable
   !  use omp_lib
    implicit none
    TYPE(mic_param_xscale),  INTENT(INOUT)   :: micpxdef
    TYPE(mic_param_default), INTENT(IN)      :: micpdef
    TYPE(mic_parameter),     INTENT(INout)   :: micparam
    TYPE(mic_input),         INTENT(INout)   :: micinput
    TYPE(mic_global_input),  INTENT(INout)   :: micglobal
    TYPE(mic_cpool),         INTENT(INOUT)   :: miccpool
    TYPE(mic_npool),         INTENT(INOUT)   :: micnpool
    TYPE(mic_output),        INTENT(INout)   :: micoutput

    integer :: jrestart,isoc14,kinetics,bgcopt,nyeqpool

    ! local variables
    real(r_2),    dimension(:), allocatable  :: xpool0,xpool1
    real(r_2),    dimension(:), allocatable  :: ypooli,ypoole,fluxsoc,cfluxa

    integer       :: ndelt,i,j,year,ip,np,ns,ny
    integer       :: nyrun,ip5
    real(r_2)     :: timex,delty,fluxdocsx,diffsocxx

    character(len=99)  :: frestart_in,frestart_out,foutput
    real(r_2)     :: cpool0, cpool1, totcinput


       allocate(xpool0(mcpool),xpool1(mcpool))
       allocate(ypooli(ms),ypoole(ms),fluxsoc(ms),cfluxa(ms))

    !   print *, 'calling vmic_param_constant'
       call vmic_param_constant(kinetics,micpxdef,micpdef,micparam)

    !   print *, 'calling vmic_init'
       call vmic_init(miccpool,micnpool)

       if(jrestart==1) call vmic_restart_read(miccpool,micnpool,frestart_in)

       ndelt   = int(24*365/delt) ! number of time step per year in "delt" unit

   do year=1,nyeqpool
      ny = year-nyeqpool
       do i=1,ntime
         call variable_time(year,i,micglobal,micinput,micnpool)

        ! calculate parameter values that depend on soil temperature or moisture (varying with time)
         call vmic_param_time(kinetics,micpxdef,micpdef,micparam,micinput,micnpool)

!$OMP PARALLEL DEFAULT(NONE) SHARED (micparam,micpxdef,micnpool,micinput,micglobal,miccpool,micoutput,micpdef,&
!$OMP kinetics,isoc14,nyeqpool,bgcopt,ndelt,zse,mp,ny,i,year) &
!$OMP PRIVATE (np,timex,delty,ns,ip,&
!$OMP xpool0,xpool1,fluxsoc,diffsocxx,ypooli,ypoole,cpool0,cpool1,totcinput,cfluxa)
!$OMP DO

        do np=1,mp

         if(micparam%bgctype(np)==bgcopt .and. micglobal%area(np)>0.0) then   ! optimizing parameters for each soil order
            micoutput%fluxcinput(np)=0.0; micoutput%fluxrsoil(np) = 0.0; micoutput%fluxcleach(np)= 0.0    ! yearly fluxes

               ! for each soil layer
               ! sum last all C pools of all layers for compute the soil respiration = input - sum(delCpool)
               ! before leaching is computed
                cpool0 =0.0; cpool1 =0.0; totcinput = 0.0
               do ns=1,ms
                 ! micinput%cinputm(np,ns)+micinput%cinputs(np,ns) in mg C/cm3/delt
                  totcinput =totcinput + (micinput%cinputm(np,ns)+micinput%cinputs(np,ns)) *1000.0 * zse(ns)   ! convert to g C/m2/delt/zse

                  do ip=1,mcpool
                     xpool0(ip) = miccpool%cpool(np,ns,ip)
                     cpool0     = cpool0  + xpool0(ip) * zse(ns) * 1000.0  ! 1000 for mg C/cm3 to g C/m2/zse
                  enddo

                 ! here the integration step is "delty" in rk4 and "ndelt" is number of "delt (hour) per year
                  timex=real(i*delt)
                  delty = real(ndelt)/(365.0*delt)  ! time step in rk4 in "24 * delt (or daily)", all C input are in " per delt"
                  call rk4modelx(timex,delty,ny,isoc14,np,ns,kinetics,micpdef,micparam,micinput,xpool0,xpool1)

                  do ip=1,mcpool
                     miccpool%cpool(np,ns,ip) = max(xpool1(ip),1.0e-8)
                     cpool1 = cpool1 + miccpool%cpool(np,ns,ip) * zse(ns) * 1000.0  ! 1000 for mg C/cm3 to g C/m2/zse
                  enddo

                ! for checking mass balance
        !          write(*,101) np,ns, micinput%cinputm(np,ns)+micinput%cinputs(np,ns),sum(xpool1(1:7)-xpool0(1:7))/real(delty), &
        !                              micinput%cinputm(np,ns)+micinput%cinputs(np,ns)-sum(xpool1(1:7)-xpool0(1:7))/real(delty)
101 format('vmicsoil input sumdelC rsoil',2(i3,1x),3(f10.6,1x))

               enddo    ! "ns"

               micoutput%fluxcinput(np)= micoutput%fluxcinput(np) + totcinput * real(delty)
               micoutput%fluxrsoil(np) = micoutput%fluxrsoil(np)  + totcinput * real(delty) + (cpool1 - cpool0)

            !! do labile carbon leaching only for kinetics=3
            !! the following leachate transport calculations caused mass imbalance: disabled temporarily
            !   if(kinetics==3) then
            !      cfluxa(:)=0.0
            !      do ns=1,ms
            !         cfluxa(ns) = sqrt(micinput%wavg(np,ns)/micinput%porosity(np,ns)) * micparam%tvac(np,ns) * miccpool%cpool(np,ns,7) * delty
            !         cfluxa(ns) = 0.0
            !         miccpool%cpool(np,ns,7) = miccpool%cpool(np,ns,7) - cfluxa(ns)
            !         if(ns==1) then
            !            miccpool%cpool(np,ns,7) = miccpool%cpool(np,ns,7) - cfluxa(ns)
            !         else
            !            miccpool%cpool(np,ns,7) = miccpool%cpool(np,ns,7) + cfluxa(ns-1) -cfluxa(ns)
            !         endif
            !      enddo
            !     ! converting flux from mg C cm-3 delty-1 to g C m-2 delty-1
            !      micoutput%fluxcleach(np) = micoutput%fluxcleach(np) + cfluxa(ms) * zse(ms) * 1000.0
            !   endif

                if(diag==1) then
                   print *, 'year day site np1', year, i, outp,micparam%diffsocx(outp)
                    do ns=1,ms
                       print *, ns, miccpool%cpool(outp,ns,:)
                    enddo
                endif

               do ip=1,mcpool
                  do ns=1,ms
                     ypooli(ns) = miccpool%cpool(np,ns,ip)      ! in mg c/cm3
                  enddo  !"ns"

                  fluxsoc(:) = 0.0  ! This flux is added in "modelx"
                  diffsocxx= micparam%diffsocx(np)

                  call bioturb(int(delty/delty),ms,zse,delty,diffsocxx,fluxsoc,ypooli,ypoole)  ! only do every 24*delt

                  do ns=1,ms
                     miccpool%cpool(np,ns,ip) = ypoole(ns)
                  enddo
               enddo ! "ip=1,mcpool"

            ! print out the time series of pool sizes
            ! if(micglobal%bgctype(np)==bgcopt) then
            !   write(*,201) year, np, miccpool%cpool(np,1,:),miccpool%cpool(np,ms,:)
201             format('vmicsoil:cpool',2(i5,1x),30(f7.4,1x))
            ! endif

            endif   !bgctype(np) = bgcopt
          enddo !"mp"
!$OMP END DO
!$OMP END PARALLEL
         enddo   !"i: day of year"
      enddo !"year"

     miccpool%cpooleq(:,:,:) = miccpool%cpool(:,:,:)

    ! call vmic_output_write(foutput,micinput,micoutput)
    ! call vmic_restart_write(frestart_out,miccpool,micnpool)

     deallocate(xpool0,xpool1)
     deallocate(ypooli,ypoole,fluxsoc,cfluxa)

    end subroutine vmicsoil_hwsd_cpu

   subroutine vmicsoil_hwsd_gpu(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,bgcopt,nyeqpool, &
                    micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)
    use mic_constant
    use mic_variable

    implicit none
    TYPE(mic_param_xscale),  INTENT(INOUT)   :: micpxdef
    TYPE(mic_param_default), INTENT(IN)      :: micpdef
    TYPE(mic_parameter),     INTENT(INout)   :: micparam
    TYPE(mic_input),         INTENT(INout)   :: micinput
    TYPE(mic_global_input),  INTENT(INout)   :: micglobal
    TYPE(mic_cpool),         INTENT(INOUT)   :: miccpool
    TYPE(mic_npool),         INTENT(INOUT)   :: micnpool
    TYPE(mic_output),        INTENT(INout)   :: micoutput

    integer :: isoc14,kinetics,bgcopt,nyeqpool
    ! local variables
    real(r_2),    dimension(:), allocatable  :: xpool0,xpool1
    real(r_2),    dimension(:), allocatable  :: ypooli,ypoole,fluxsoc,cfluxa


    integer       :: ndelt,i,j,year,ip,np,ns,ny
    integer       :: nyrun
    real(r_2)     :: timex,delty,fluxdocsx,diffsocxx

    integer    :: jrestart
    character(len=99) :: frestart_in,frestart_out,foutput
    real(r_2)  :: cpool0, cpool1, totcinput

   ! local variables
   ! for numerical solution
   real(r_2),    parameter            :: tol = 1.0E-04
   real(r_2),    parameter            :: tolx = 0.0001
   real(r_2),    parameter            :: tolf = 0.000001
   integer,      parameter            :: ntrial = 100

   ! local variables
   real(r_2)                      :: deltD !,tot0,tot1,totflux
   real(r_2), dimension(:), allocatable    :: xzse
   real(r_2), dimension(:), allocatable    :: sdepthx
   real(r_2)                      :: coeffA, coeffB
   real(r_2), dimension(:), allocatable    :: at,bt,ct,rt
   real(r_2), dimension(:), allocatable    :: xpool
   real(r_2)  :: cleachloss

    allocate(xpool0(mcpool),xpool1(mcpool))
    allocate(ypooli(ms),ypoole(ms),fluxsoc(ms),cfluxa(ms))
    allocate(xzse(ms))
    allocate(sdepthx(ms+1))
    allocate(at(ms),bt(ms),ct(ms),rt(ms))
    allocate(xpool(ms))

      call vmic_param_constant(kinetics,micpxdef,micpdef,micparam)
      call vmic_init(miccpool,micnpool)
      call vmic_param_time(kinetics,micpxdef,micpdef,micparam,micinput,micnpool)

      if(jrestart==1) call vmic_restart_read(miccpool,micnpool,frestart_in)

      ndelt   = int(24*365/delt) ! number of time step per year in "delt" unit

! "data copyin" : all data accessed by GPU
! "create": intermediate variables used by GPU
! "copyout" copy data out from GPU
! "private": every cell-dependent variables used in paralelling computing

!$acc data copyin(micpdef,micparam,miccpool,micinput,micoutput,micglobal,  &
!$acc bgcopt,ndelt,zse,kinetics,nyeqpool,isoc14)      &
!$acc create(delty,timex,fluxsoc,year,ny,i,ns,np,ip,xpool0,xpool1,ypooli,ypoole,diffsocxx,cfluxa, &
!$acc j,deltD,xzse,sdepthx,coeffA,coeffB,at,bt,ct,rt,xpool,cleachloss)              &
!$acc copyout(miccpool%cpool,micoutput%fluxcinput,micoutput%fluxrsoil,micoutput%fluxcleach)
!$acc PARALLEL LOOP                                                      &
!$acc private(delty,timex,fluxsoc,year,ny,i,ns,ip,np,xpool0,xpool1,ypooli,ypoole,diffsocxx,cfluxa,&
!$acc j,deltD,xzse,sdepthx,coeffA,coeffB,at,bt,ct,rt,xpool,cleachloss)

      do np=1,mp

         if(micparam%bgctype(np)==bgcopt .and. micglobal%area(np)>0.0) then   ! optimizing parameters for each soil order

         do year=1,nyeqpool
            ny = year-nyeqpool
            micoutput%fluxcinput(np)=0.0; micoutput%fluxrsoil(np) = 0.0; micoutput%fluxcleach(np)= 0.0    ! yearly fluxes

            do i=1,365

               call variable_time(year,i,micglobal,micinput,micnpool)

               ! calculate parameter values that depend on soil temperature or moisture (varying with time)
               call vmic_param_time(kinetics,micpxdef,micpdef,micparam,micinput,micnpool)

               ! for each soil layer
               ! sum last all C pools of all layers for compute the soil respiration = input - sum(delCpool)
               ! before leaching is computed
                cpool0 =0.0; cpool1 =0.0; totcinput = 0.0
               do ns=1,ms
                 ! micinput%cinputm(np,ns)+micinput%cinputs(np,ns) in mg C/cm3/delt
                  totcinput =totcinput + (micinput%cinputm(np,ns)+micinput%cinputs(np,ns)) *1000.0 * zse(ns)   ! convert to g C/m2/delt/zse

                  do ip=1,mcpool
                     xpool0(ip) = miccpool%cpool(np,ns,ip)
                     cpool0     = cpool0  + xpool0(ip) * zse(ns) * 1000.0  ! 1000 for mg C/cm3 to g C/m2/zse
                  enddo

                 ! here the integration step is "delty" in rk4 and "ndelt" is number of "delt (hour) per year
                  timex=real(i*delt)
                  delty = real(ndelt)/(365.0*delt)  ! time step in rk4 in "24 * delt (or daily)", all C input are in " per delt"

                  call rk4modelx(timex,delty,ny,isoc14,np,ns,kinetics,micpdef,micparam,micinput,xpool0,xpool1)

                  do ip=1,mcpool
                     miccpool%cpool(np,ns,ip) = max(xpool1(ip),1.0e-8)
                     cpool1 = cpool1 + miccpool%cpool(np,ns,ip) * zse(ns) * 1000.0  ! 1000 for mg C/cm3 to g C/m2/zse
                  enddo

               enddo    ! "ns"

               micoutput%fluxcinput(np)= micoutput%fluxcinput(np) + totcinput * real(delty)
               micoutput%fluxrsoil(np) = micoutput%fluxrsoil(np)  + totcinput * real(delty) + (cpool1 - cpool0)

               ! do labile carbon leaching only for kinetics=3
               ! the following leachate transport calculations caused mass imbalance: disabled temporarily
            !   if(kinetics==3) then
            !     cfluxa(:)=0.0
            !     do ns=1,ms
            !        cfluxa(ns) = sqrt(micinput%wavg(np,ns)/micinput%porosity(np,ns)) * micparam%tvac(np,ns) * miccpool%cpool(np,ns,7) * delty
            !        if(ns==1) then
            !           miccpool%cpool(np,ns,7) = miccpool%cpool(np,ns,7) - cfluxa(ns)
            !        else
            !           miccpool%cpool(np,ns,7) = miccpool%cpool(np,ns,7) + cfluxa(ns-1) -cfluxa(ns)
            !        endif
            !     enddo
            !    ! converting flux from mg C cm-3 delty-1 to g C m-2 delty-1
            !    micoutput%fluxcleach(np) = micoutput%fluxcleach(np) + cfluxa(ms) * zse(ms) * 1000.0
            !  endif

                if(diag==1) then
                   print *, 'year day site np1', year, i, outp,micparam%diffsocx(outp)
                    do ns=1,ms
                       print *, ns, miccpool%cpool(outp,ns,:)
                    enddo
                endif

               do ip=1,mcpool
                  do ns=1,ms
                     ypooli(ns) = miccpool%cpool(np,ns,ip)      ! in mg c/cm3
                  enddo  !"ns"

                  fluxsoc(:) = 0.0  ! This flux is added in "modelx"
                  diffsocxx= micparam%diffsocx(np)

                 !Move bioturb here to work around memory auto allocation failure
                 !call bioturb(int(delty/delty),ms,zse,delty,diffsocxx,fluxsoc,ypooli,ypoole)  ! only do every 24*delt
                 !subroutine bioturb(ndelt,ms,zse,delt,diffsocxx,fluxsoc,xpooli,xpoole)

                  sdepthx(1) = 0.0          ! depth of a layer from the top (x_0.5=0.0 eg soil surface)
                  do j=2,ms+1
                     sdepthx(j) = sdepthx(j-1) + zse(j-1)*100.0     ! depth of the bottom of each layer (eg x_j+0.5)
                                                                    !*100 to convert from m to cm
                  enddo

                  do j=1,ms
                     xzse(j) = 0.5 * (sdepthx(j) + sdepthx(j+1))    ! depth of midpoint of a layer j  (x_j)
                  enddo

                  deltD = diffsocxx * delty
                  xpool = ypooli

                  !do i=1,1 ( int(delty/delty) == 1 )
                     do j=1,ms
                        if(j==1) then
                           coeffB = 1.0/(sdepthx(2)-sdepthx(1))
                           coeffA = deltD*coeffB/(xzse(2)-xzse(1))
                           ! Crank-Nicholson
                           at(1) = 0.0
                           bt(1) = 1.0 + 0.5 * coeffA
                           ct(1) =     - 0.5 * coeffA
                           rt(1) = (1.0-0.5*coeffA) * xpool(1) + 0.5 * coeffA * xpool(2) &
                                 +  fluxsoc(1) * delt
                        endif
                        if(j>1.and.j<ms) then
                          coeffA = deltD/((xzse(j+1)-xzse(j))*(sdepthx(j+1)-sdepthx(j)))
                          coeffB = (xzse(j+1)-xzse(j))/(xzse(j)-xzse(j-1))
                          ! Crank-Nicholson
                          at(j) =    -0.5 * coeffA * coeffB
                          bt(j) = 1.0+0.5 * coeffA *(1.0+coeffB)
                          ct(j) =    -0.5 * coeffA
                          rt(j) = 0.5 * coeffA * coeffB * xpool(j-1)        &
                                  +(1.0-0.5* coeffA*(1.0+coeffB))*xpool(j)  &
                                  + 0.5* coeffA * xpool(j+1)                &
                                  + fluxsoc(j) *delt
                        endif
                        if(j==ms) then
                            coeffA = deltD/((xzse(ms)-xzse(ms-1))*(sdepthx(ms+1) - sdepthx(ms)))
                          ! Crank-Nicholson
                            at(ms) = -0.5 * coeffA
                            bt(ms) = 1.0 + 0.5 * coeffA
                            ct(ms) = 0.0
                            rt(ms) = 0.5* coeffA  * xpool(ms-1) + (1.0-0.5 * coeffA) * xpool(ms) &
                                   + fluxsoc(ms) * delt
                        endif
                     enddo
                     call tridag(at,bt,ct,rt,xpool,ms)
                  !enddo
                  ypoole = xpool

                  !!! end bioturb

                  do ns=1,ms
                     miccpool%cpool(np,ns,ip) = ypoole(ns)
                  enddo
               enddo ! "ip=1,mcpool"

               ! computing daily leaching loss from bottom-layer LWMC
               cleachloss = micparam%tvac(np,ms) * sqrt(micinput%wavg(np,ms)/micinput%porosity(np,ms)) *  miccpool%cpool(np,ms,7)  * 24.0
               cleachloss = max(0.0,min(cleachloss,miccpool%cpool(np,ms,7)))
               micoutput%fluxcleach(np) = micoutput%fluxcleach(np) + cleachloss
               miccpool%cpool(np,ms,7)  = miccpool%cpool(np,ms,7)  - cleachloss

            enddo   ! "day"
         enddo !"year"
      endif   !bgctype(np)=bgcopt
   enddo !"mp"
!$ACC END PARALLEL
!$ACC END DATA

    miccpool%cpooleq(:,:,:) = miccpool%cpool(:,:,:)

   !  call vmic_output_write(foutput,micinput,micoutput)
   !  call vmic_restart_write(frestart_out,miccpool,micnpool)
    deallocate(xpool0,xpool1)
    deallocate(ypooli,ypoole,fluxsoc,cfluxa)
    deallocate(xzse)
    deallocate(sdepthx)
    deallocate(at,bt,ct,rt)
    deallocate(xpool)

    end subroutine vmicsoil_hwsd_gpu


    subroutine calcost_c14(nx,isoc14,bgcopt,xopt,micparam,miccpool,micinput,totcost)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(INOUT) :: micparam
    TYPE(mic_cpool),     INTENT(INOUT) :: miccpool
    TYPE(mic_input),     INTENT(IN)    :: micinput
    integer :: nx,isoc14,bgcopt
    real*8  :: totcost
    real*8, dimension(16)              :: xopt

    ! cost function
    real(r_2), dimension(:), allocatable           :: xcost,xobs,xobsp,xobsm
    real(r_2), dimension(:), allocatable           :: xmod,xmodp,xmodm !! weighted modelled SOC, POC and MAOC
    real(r_2), dimension(:), allocatable           :: xfracpmod,xfracmmod,xfracpobs,xfracmobs
    integer   :: np,ns,ip
    real(r_2)  :: xbdz,small

    real(r_2) :: xtop,xbot,x1,x2 !! cm
    real(r_2) :: weight



      allocate(xcost(mp),xobs(mp),xobsp(mp),xobsm(mp))
      allocate(xmod(mp),xmodp(mp),xmodm(mp))
      allocate(xfracpmod(mp),xfracmmod(mp),xfracpobs(mp),xfracmobs(mp))

      small = 1.0e-6
      xcost(:)=0.0
      xobs(:)=0.0; xobsp(:)=0.0;xobsm(:)=0.0

      do np=1,mp
         if(micparam%bgctype(np)==bgcopt) then
            xobs(np)  = micparam%csoilobs(np,1)   !! only one observation for each site
            xobsp(np) = micparam%csoilobsp(np,1)
            xobsm(np) = micparam%csoilobsm(np,1)
            xfracpobs(np) = min(1.0,max(0.0,xobsp(np)/(xobsp(np)+xobsm(np))))
            xfracmobs(np) = 1.0 - xfracpobs(np)

           !! calculated weighted average of modelled values to correspond to observations
            xmod(np)=0.0; xmodp(np)=0.0; xmodm(np)=0.0
            xtop = real(micparam%top(np))*0.01   ! convert from cm to m
            xbot = real(micparam%bot(np))*0.01

            x1 =0.0; x2 =0.0
            do ns=1,ms
               x2=x1+zse(ns)
               if (x1 <= xbot .and. x2 >= xtop) then
                  weight = min(xbot,x2) - max(xtop,x1)
                  xmod(np)  = xmod(np)  + weight * sum(miccpool%cpooleq(np,ns,3:9)) !! unit: mg C cm-3
                  xmodp(np) = xmodp(np) + weight * (sum(miccpool%cpooleq(np,ns,3:5))+sum(miccpool%cpooleq(np,ns,7:8))) !! observed poc = modelled (aggregate c + lwmc)
                  xmodm(np) = xmodm(np) + weight * (miccpool%cpooleq(np,ns,6) + miccpool%cpooleq(np,ns,9))             !! observed maoc = modelled maoc
               end if
               x1=x2
            enddo
           ! this assumes bulk density constant in the vertical
           xmod(np)   = 1000.0 * xmod(np) /((xbot-xtop) * micinput%bulkd(np,1))   !! unit: g C/kg soil
           xmodp(np)  = 1000.0 * xmodp(np)/((xbot-xtop) * micinput%bulkd(np,1))
           xmodm(np)  = 1000.0 * xmodm(np)/((xbot-xtop) * micinput%bulkd(np,1))

           xfracpmod(np) = min(1.0,max(0.0,xmodp(np)/(xmodp(np)+xmodm(np))))
           xfracmmod(np) = 1.0 - xfracpmod(np)

           miccpool%cpooleqp(np) = xmodp(np)
           miccpool%cpooleqm(np) = xmodm(np)
           ! avoid dividing by zero
           miccpool%c12pooleqp(np) = max(small,miccpool%c12pooleqp(np))
           micparam%c14soilobsp(np)= max(small,micparam%c14soilobsp(np))
           miccpool%c12pooleqm(np) = max(small,miccpool%c12pooleqm(np))
           micparam%c14soilobsm(np)= max(small,micparam%c14soilobsm(np))

           if (xmod(np) >= 1000.0 .or. xmodp(np) >= 1000.0 .or. xmodm(np) >= 1000.0) then
               print *, 'abnormal value of model simulation site=', micparam%siteid(np)
               print *, 'parameter values = ',  xopt(1:nx)
            !   stop
           endif

           if(isoc14 == 0) then
             if((xobsp(np)+xobsm(np))>0.0 .and. (xobsp(np)+xobsm(np))<1000.0) then
        !     xcost(np) = xcost(np) + ((xmodp(np) - xobsp(np))/xobsp(np))**2 +((xmodm(np) - xobsm(np))/xobsm(np))**2
              xcost(np) = xcost(np) + (xfracpmod(np)-xfracpobs(np))** 2 + (xfracmmod(np) - xfracmobs(np))**2
             endif
             write(91,901) micparam%siteid(np),micparam%pft(np),micparam%top(np),micparam%bot(np),xobs(np),xmod(np),xobsp(np),xmodp(np),xobsm(np),xmodm(np)

             do ns = 1,ms
                write(92,*) micparam%siteid(np),micparam%pft(np), ns, (1000.0*miccpool%cpooleq(np,ns,ip)/micinput%bulkd(np,ns),ip=1,mcpool)
             enddo

           else
             xcost(np) = xcost(np) +(micparam%c14soilobsp(np) - xmodp(np)/miccpool%c12pooleqp(np))**2  &
                                   +(micparam%c14soilobsm(np) - xmodm(np)/miccpool%c12pooleqm(np))**2

             write(93,901) micparam%siteid(np),micparam%pft(np),micparam%top(np),micparam%bot(np), &
                           xfracpobs(np),xfracpmod(np),xfracmobs(np),xfracmmod(np),                &
                           micparam%c14soilobsp(np),xmodp(np)/miccpool%c12pooleqp(np),             &
                           micparam%c14soilobsm(np),xmodm(np)/miccpool%c12pooleqm(np)

             do ns = 1,ms
                write(94,*) micparam%siteid(np),micparam%pft(np), ns, (1000.0*miccpool%cpooleq(np,ns,ip)/micinput%bulkd(np,ns),ip=1,mcpool)
             enddo
           endif

        endif  !bgctype(np)=bgcopt
      enddo  !"np"
      totcost = sum(xcost(1:mp))

      deallocate(xcost,xobs,xobsp,xobsm)
      deallocate(xmod,xmodp,xmodm)
      deallocate(xfracpmod,xfracmmod,xfracpobs,xfracmobs)

901   format(i5,2x,3(i4,2x),12(f10.4,2x))
921   format(i7,2x,2(i4,2x),10(f12.4,1x))
    end subroutine calcost_c14



    subroutine calcost_frc(nx,bgcopt,xopt,micpxdef,micparam,miccpool,micinput,totcost)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_param_xscale),INTENT(IN)  :: micpxdef
    TYPE(mic_parameter), INTENT(IN)    :: micparam
    TYPE(mic_cpool),     INTENT(INOUT) :: miccpool
    TYPE(mic_input),     INTENT(IN)    :: micinput
    integer :: nx,bgcopt
    real*8  :: totcost
    real*8, dimension(16)              :: xopt

    ! cost function
    real(r_2), dimension(:), allocatable        :: xcost,xobs,xobsp,xobsm
    real(r_2), dimension(:), allocatable        :: xmod,xmodp,xmodm !! weighted modelled SOC, POC and MAOC
    integer   :: np,ns,ip,ipsite
    real(r_2)  :: xbdz

    real(r_2) :: xtop,xbot,x1,x2 !! cm
    real(r_2) :: weight
    real(r_2),dimension(:), allocatable         :: xmodfracp,xmodfracm,xobsfracp,xobsfracm


      allocate(xcost(mp),xobs(mp),xobsp(mp),xobsm(mp))
      allocate(xmod(mp),xmodp(mp),xmodm(mp))
      allocate(xmodfracp(mp),xmodfracm(mp),xobsfracp(mp),xobsfracm(mp))

      ipsite=-999
      xcost(:)=0.0
      xobs(:)=0.0; xobsp(:)=0.0;xobsm(:)=0.0;xmodfracp(:)=0.0;xmodfracm(:)=0.0;xobsfracp(:)=0.0;xobsfracm(:)=0.0

      do np=1,mp
         if(micparam%bgctype(np)==bgcopt) then
            ipsite=np
            xobs(np)  = micparam%csoilobs(np,1) !! only one observation for each site
            xobsp(np) = micparam%csoilobsp(np,1)
            xobsm(np) = micparam%csoilobsm(np,1)

           !! calculated weighted average of modelled values to correspond to observations
            xmod(np)=0.0; xmodp(np)=0.0; xmodm(np)=0.0
            xtop = real(micparam%top(np))*0.01   ! convert from cm to m
            xbot = real(micparam%bot(np))*0.01

            x1 =0.0; x2 =0.0
            do ns=1,ms
               x2=x1+zse(ns)
               if (x1 <= xbot .and. x2 >= xtop) then
                  weight = min(xbot,x2) - max(xtop,x1)
                  xmod(np)  = xmod(np)  + weight * sum(miccpool%cpooleq(np,ns,3:9)) !! unit: mg C cm-3
                  xmodp(np) = xmodp(np) + weight * (sum(miccpool%cpooleq(np,ns,3:5))+sum(miccpool%cpooleq(np,ns,7:8))) !! observed poc = modelled (aggregate c + lwmc)
                  xmodm(np) = xmodm(np) + weight * (miccpool%cpooleq(np,ns,6) + miccpool%cpooleq(np,ns,9))             !! observed maoc = modelled maoc
               end if
               x1=x2
            enddo

           xmod(np)   = 1000.0 * xmod(np) /((xbot-xtop)*micinput%bulkd(np,1))   !! unit: g C/kg soil
           xmodp(np)  = 1000.0 * xmodp(np)/((xbot-xtop)*micinput%bulkd(np,1))
           xmodm(np)  = 1000.0 * xmodm(np)/((xbot-xtop)*micinput%bulkd(np,1))

           xmodfracp(np) = xmodp(np)/(xmodp(np)+xmodm(np)+1.0e-6)
           xmodfracm(np) = xmodm(np)/(xmodp(np)+xmodm(np)+1.0e-6)
           xobsfracp(np) = xobsp(np)/(xobsp(np)+xobsm(np)+1.0e-6)
           xobsfracm(np) = xobsm(np)/(xobsp(np)+xobsm(np)+1.0e-6)

           if (xmod(np) >= 1000.0 .or. xmodp(np) >= 1000.0 .or. xmodm(np) >= 1000.0) then
               print *, 'abnormal value of model simulation site=', micparam%siteid(np)
               print *, 'parameter values = ',  xopt(1:nx)
               print *, xmodp(np),xmodm(np)
            !   stop
            endif


            if(xobsp(np) >0.0 .and. xobsm(np)>0.0) then  !! POC is not available for some sites
               xcost(np) = xcost(np) + 200.0*((xmodfracp(np) - xobsfracp(np))**2 + (xmodfracm(np) - xobsfracm(np))**2) &
                                     + 0.01 *(xmodm(np)+xmodp(np) - xobsm(np)-xobsp(np))**2
            endif

            write(91,901) micparam%dataid(np),micparam%siteid(np),micparam%pft(np),micparam%top(np),micparam%bot(np), &
                          xobsp(np),xmodp(np),xobsm(np),xmodm(np),xobsfracp(np),xmodfracp(np),xobsfracm(np),xmodfracm(np)

            do ns = 1,ms
               write(92,*) micparam%dataid(np),micparam%siteid(np),micparam%pft(np), ns, &
                           (1000.0*miccpool%cpooleq(np,ns,ip)/micinput%bulkd(np,ns),ip=1,mcpool)
            enddo


        endif  !bgctype(np)=bgcopt
      enddo  !"np"
      totcost = sum(xcost(1:mp))




!      write(93,911) bgcopt,ipsite,  &
!                    xmodp(np),xobsp(np),xmodm(np),xobsm(np),xmodfracp(np),xobsfracp(np),xmodfracm(np),xobsfracm(np),totcost


      deallocate(xcost,xobs,xobsp,xobsm)
      deallocate(xmod,xmodp,xmodm)
      deallocate(xmodfracp,xmodfracm,xobsfracp,xobsfracm)

901   format(5(i6,1x),10(f12.4,1x))
921   format(4(i6,1x),14(f12.4,1x))
911   format(2(i6,1x),4(f12.4,1x),4(f9.5,1x),e15.6)

    end subroutine calcost_frc

    subroutine calcost_soc(nx,isoc14,bgcopt,xopt,micpxdef,micparam,miccpool,micinput,totcost)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_param_xscale), INTENT(IN) :: micpxdef
    TYPE(mic_parameter), INTENT(IN)    :: micparam
    TYPE(mic_cpool),     INTENT(INOUT) :: miccpool
    TYPE(mic_input),     INTENT(IN)    :: micinput
    integer :: nx,isoc14,bgcopt
    real*8  :: totcost
    real*8, dimension(16)              :: xopt
    ! cost function
    real(r_2), dimension(:),   allocatable        :: xcost
    real(r_2), dimension(:,:), allocatable        :: xobs, xmod
    real(r_2), dimension(5)            :: xobsv, xmodv
    integer   :: np,ns,ipsite,v,ip
    real(r_2)  :: xbdz

    allocate(xcost(mp))
    allocate(xobs(mp,ms), xmod(mp,ms))

      ! write out the observed and modelled SOC (gC/kg soil) for each site
      xobs(:,:)=0.0; xmod(:,:)=0.0;xcost(:)=0.0
      ipsite=-999
   do np=1,mp
      if(micparam%bgctype(np)==bgcopt) then
         ipsite=np
         xcost(np) = 0.0

         do ns = 1,ms
            xobs(np,ns) = micparam%csoilobs(np,ns)
            xmod(np,ns) = 1000.0 * sum(miccpool%cpooleq(np,ns,3:mcpool))/micinput%bulkd(np,ns)

               if(xmod(np,ns) <0.0.or.xmod(np,ns) >1.0e3) then
                  print *, 'abnormal value of model simulation site=', micparam%siteid(np)
                  print *, 'parameter values = ',  xopt(1:nx)
                  print *,  xobs(np,ns),xmod(np,ns),xobs(np,ns)-xmod(np,ns)
                  print *, ' modelled pool size= ', ns,miccpool%cpooleq(np,ns,:)/micinput%bulkd(np,ns)
                  !stop
               endif
         enddo !"ns"

         xobsv(1) = sum(xobs(np,1:3))/3.0
         xobsv(2) = sum(xobs(np,4:6))/3.0
         xobsv(3) = sum(xobs(np,7:9))/3.0
         xobsv(4) = sum(xobs(np,10:12))/3.0
         xobsv(5) = sum(xobs(np,13:15))/3.0

         xmodv(1) = sum(xmod(np,1:3))/3.0
         xmodv(2) = sum(xmod(np,4:6))/3.0
         xmodv(3) = sum(xmod(np,7:9))/3.0
         xmodv(4) = sum(xmod(np,10:12))/3.0
         xmodv(5) = sum(xmod(np,13:15))/3.0

         do v = 1,5
             if(xobsv(v) /= -999.0 .and. xmodv(v) <1.0e3) then
                xcost(np) = xcost(np) + (xmodv(v) - xobsv(v))**2
             else
                if(xmodv(v) >1.0e3)  xcost(np) = xcost(np) + 9.0e9
                if(xobsv(v)==-999.0) xcost(np) = xcost(np)
             endif
             write(91,901) micparam%siteid(np),micparam%pft(np),v,xobsv(v),xmodv(v),xmodv(v)-xobsv(v),xcost(np)
         enddo

         do ns = 1,ms
            write(92,921) micparam%siteid(np),micparam%pft(np), ns, (1000.0*miccpool%cpooleq(np,ns,ip)/micinput%bulkd(np,ns),ip=1,mcpool)
         enddo
      endif !"pft"
   enddo  !"np"
   totcost = sum(xcost(1:mp))

      write(93,911) bgcopt,ipsite,  &
                    micparam%V1(ipsite,1),micparam%K1(ipsite,1),micparam%fm(ipsite,1), &
                    micparam%fs(ipsite,1),micparam%tvmicR(ipsite,1)*365.0*24.0, &
                    micparam%tvppool(ipsite,1)*365.0*24.0,micparam%tvcpool(ipsite,1)*365.0*24.0, &
                    micparam%tvac(ipsite,1)*365.0*24.0, micparam%kdesorp(ipsite,1), &
                    micparam%qmaxcoeff(ipsite,1), micparam%diffsocx(ipsite)*24.0,micpxdef%xNPP(1),micparam%fracroot(ipsite,1),  &
                    micparam%V1(ipsite,1), totcost

901   format(i7,2x,2(i4,2x),3(f10.4,2x),f13.1)
921   format(i7,2x,2(i4,2x),10(f12.4,1x))
911   format(2(i6,1x),14(f12.4,1x),e15.6)
    deallocate(xcost)
    deallocate(xobs, xmod)

    end subroutine calcost_soc


    subroutine calcost_hwsd(nx,bgcopt,xopt,micpxdef,micparam,miccpool,micinput,micglobal,totcost)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_param_xscale), INTENT(IN)    :: micpxdef
    TYPE(mic_parameter),    INTENT(IN)    :: micparam
    TYPE(mic_cpool),        INTENT(INOUT) :: miccpool
    TYPE(mic_input),        INTENT(IN)    :: micinput
    TYPE(mic_global_input), INTENT(IN)    :: micglobal

    integer :: nx,bgcopt
    real*8  :: totcost
    real*8, dimension(16)              :: xopt
    ! cost function
    real(r_2), dimension(:),   allocatable        :: xcost
    real(r_2), dimension(:,:), allocatable        :: xobs, xmod
    real(r_2)                                     :: fracpocm,fracmaocm,fracmicm,fraclabm
    integer   :: np,ns,ipsite,v,ip
    real(r_2)  :: xbdz

    allocate(xcost(mp))
    allocate(xobs(mp,ms), xmod(mp,ms))

      xobs(:,:)=0.0; xmod(:,:)=0.0;xcost(:)=0.0
      ipsite=-999
    do np=1,mp
       if(micparam%bgctype(np)==bgcopt .and. micglobal%area(np) > 0.0) then
         ipsite=np
         xcost(np) = 0.0
         v = 0
         do ns = 1,ms
            xobs(np,ns) = micparam%csoilobs(np,ns)
            xmod(np,ns) = 1000.0 * sum(miccpool%cpooleq(np,ns,3:mcpool))/micinput%bulkd(np,ns)   ! gC/kg soil

            if(xmod(np,ns) <0.0.or.xmod(np,ns) >1.0e3) then
               print *, 'abnormal value of model simulation site=', micparam%siteid(np)
               print *, 'parameter values = ',  xopt(1:nx)
               print *,  xobs(np,ns),xmod(np,ns),xobs(np,ns)-xmod(np,ns)
               print *, ' modelled pool size= ', ns,miccpool%cpooleq(np,ns,:)/micinput%bulkd(np,ns)
               !stop
            endif

            if(xobs(np,ns) > 0.0 .and. xobs(np,ns) <1.0e3) then
               xcost(np) = xcost(np) + (xobs(np,ns)-xmod(np,ns))**2
               fracpocm = (sum(miccpool%cpooleq(np,ns,3:8))-miccpool%cpooleq(np,ns,6))/(sum(miccpool%cpooleq(np,ns,3:mcpool))+1.0e-6)
               fracmaocm = (miccpool%cpooleq(np,ns,6)+miccpool%cpooleq(np,ns,9))/(sum(miccpool%cpooleq(np,ns,3:mcpool))+1.0e-6)
               fracmicm = (miccpool%cpooleq(np,ns,3)+miccpool%cpooleq(np,ns,4))/(sum(miccpool%cpooleq(np,ns,3:mcpool))+1.0e-6)
               fraclabm = miccpool%cpooleq(np,ns,7)/(sum(miccpool%cpooleq(np,ns,3:mcpool))+1.0e-6)

               write(91,901) micparam%siteid(np),micparam%pft(np),micparam%isoil(np),micparam%sorder(np),ns, &
                             xobs(np,ns),xmod(np,ns),fracpocm,fracmaocm,fracmicm,fraclabm
            endif
         enddo !"ns"

         do ns = 1,ms
            write(92,921) micparam%siteid(np),micparam%sorder(np), ns, (1000.0*miccpool%cpooleq(np,ns,ip)/micinput%bulkd(np,ns),ip=1,mcpool)
         enddo

      endif !"pft"
   enddo  !"np"
   totcost = sum(xcost(1:mp))

    deallocate(xcost)
    deallocate(xobs, xmod)

901   format(i6,1x,4(i3,1x),10(f12.4,1x))
921   format(i7,2x,2(i4,2x),10(f12.4,1x))
    end subroutine calcost_hwsd

    subroutine rk4modelx(timex,delty,ny,isoc14,np,ns,kinetics,micpdef,micparam,micinput,xpool0,xpool1)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_param_default), INTENT(IN)  :: micpdef
    TYPE(mic_parameter),     INTENT(IN)  :: micparam
    TYPE(mic_input),         INTENT(IN)  :: micinput
    integer      :: np,ns, kinetics,ny,isoc14
    real(r_2)    :: timex,delty,h
    real(r_2),   dimension(mcpool),intent(inout)     :: xpool0,xpool1
    real(r_2),   dimension(mcpool)                   :: y1,y2,y3,y4,dy1dt,dy2dt,dy3dt,dy4dt

     h=delty
     y1(:) = xpool0(:)

     call vmic_c(ny,isoc14,np,ns,kinetics,micpdef,micparam,micinput,y1,dy1dt)
     y2(:) = y1(:) + 0.5 * h * dy1dt(:)
     call vmic_c(ny,isoc14,np,ns,kinetics,micpdef,micparam,micinput,y2,dy2dt)
     y3(:) = y1(:) + 0.5 * h * dy2dt(:)
     call vmic_c(ny,isoc14,np,ns,kinetics,micpdef,micparam,micinput,y3,dy3dt)
     y4(:) = y1(:) +       h * dy3dt(:)
     call vmic_c(ny,isoc14,np,ns,kinetics,micpdef,micparam,micinput,y4,dy4dt)
    ! RK4
     xpool1(:) = xpool0(:) + (dy1dt(:)/6.0 + dy2dt(:)/3.0 + dy3dt(:)/3.0 + dy4dt(:)/6.0) * h

    ! Euler
    ! xpool1(:) = xpool0(:) + dy1dt(:) * h


!    write(*,101) np,ns,delty,micinput%cinputm(np,ns)+micinput%cinputs(np,ns),sum(dy1dt(1:7)), &
!                 micinput%cinputm(np,ns)+micinput%cinputs(np,ns)-sum(dy1dt(1:7))
101 format('rk4: input, sumdelc rsoil',2(i3,1x),f6.2,1x,3(f10.6,1x))

    end subroutine rk4modelx

    subroutine Kmt(micpxdef,micpdef,micparam,micinput)
      ! unit: mg Mic C/cm3
      use mic_constant
      use mic_variable
      implicit none
      TYPE(mic_param_xscale), INTENT(IN)       :: micpxdef
      TYPE(mic_param_default), INTENT(IN)      :: micpdef
      TYPE(mic_parameter),     INTENT(INOUT)   :: micparam
      TYPE(mic_input),         INTENT(IN)      :: micinput


      ! local variable
      real(r_2), dimension(:,:), allocatable   :: xkclay,km,kmx
      integer :: nopt,np,ns

     allocate(xkclay(mp,ms),km(mp,ms),kmx(mp,ms))
     do np=1,mp
      do ns=1,ms
         nopt=micparam%bgctype(np)
         xkclay(np,ns) = 1.0/(2.0*exp(-2.0*sqrt(micinput%clay(np,ns))))
         km(np,ns) =  micpxdef%xak(nopt) * micpdef%ak * exp(micpdef%sk * micinput%tavg(np,ns) + micpdef%bk)
         micparam%K1(np,ns) =  km(np,ns)/micpdef%xk1
         micparam%K3(np,ns) =  km(np,ns) * xkclay(np,ns)/micpdef%xk3
         micparam%J1(np,ns) =  km(np,ns)/micpdef%xj1
         micparam%J3(np,ns) =  km(np,ns) * xkclay(np,ns)/micpdef%xj3

         kmx(np,ns) =  micpxdef%xak(nopt) * micpdef%ak * exp(micpdef%skx * micinput%tavg(np,ns) + micpdef%bk)
         micparam%K2(np,ns) =  kmx(np,ns)/micpdef%xk2
         micparam%J2(np,ns) =  kmx(np,ns)/micpdef%xj2
       enddo
      enddo

       if(diag==1.and.np==outp) then
          print *, 'Kmt',micinput%clay(outp,1),micinput%tavg(outp,1),km(outp,1),kmx(outp,1)
          print *, 'K1=',micparam%K1(outp,1)
          print *, 'K2=',micparam%K2(outp,1)
          print *, 'K3=',micparam%K3(outp,1)
          print *, 'J1=',micparam%J1(outp,1)
          print *, 'J2=',micparam%J2(outp,1)
          print *, 'J3=',micparam%J3(outp,1)
       endif
      deallocate(xkclay,km,kmx)

    end subroutine Kmt

    subroutine Vmaxt(micpxdef,micpdef,micparam,micinput)
      ! mg Cs per mg mic C per hour
      use mic_constant
      use mic_variable
      implicit none
      TYPE(mic_param_xscale),  INTENT(IN)     :: micpxdef
      TYPE(mic_param_default), INTENT(IN)     :: micpdef
      TYPE(mic_parameter),     INTENT(INOUT)  :: micparam
      TYPE(mic_input),         INTENT(IN)     :: micinput


      ! local variables
      real(r_2),dimension(:,:), allocatable :: vmax
      integer :: nopt,np,ns
      real(r_2), dimension(:), allocatable   :: sdepthz

      allocate(vmax(mp,ms))
      allocate(sdepthz(ms))

      sdepthz=0.0

      do np=1,mp
           sdepthz=0.0
           nopt=micparam%bgctype(np)
          do ns=1,ms
            if(ns==1) then
               sdepthz(ns) = 0.5 * micparam%sdepth(np,ns)
            else
                sdepthz(ns) = sdepthz(ns-1) + micparam%sdepth(np,ns)
            endif
        !   vmax(np,ns) =  micpxdef%xav * micpdef%av * exp(micpdef%sv*micinput%tavg(np,ns) + micpdef%bv) * delt
        !   vmax(np,ns) =  exp(-2.0* sdepthz(ns)) * micpxdef%xav(npft) * micpdef%av * exp(micpdef%sv*micinput%tavg(np,ns) + micpdef%bv) * delt
           vmax(np,ns) =  exp(-micpdef%vmaxbeta * micpxdef%xvmaxbeta(nopt) * sdepthz(ns))     &
                                 * micpxdef%xav(nopt) * micpdef%av * exp(micpdef%sv*micinput%tavg(np,ns) + micpdef%bv)  * delt

           micparam%V1(np,ns)   =  micpdef%xv1 * vmax(np,ns)
           micparam%V2(np,ns)   =  micpdef%xv2 * vmax(np,ns)
           micparam%V3(np,ns)   =  micpdef%xv3 * vmax(np,ns)

           micparam%W1(np,ns)   =  micpdef%xw1 * vmax(np,ns)
           micparam%W2(np,ns)   =  micpdef%xw2 * vmax(np,ns)
           micparam%W3(np,ns)   =  micpdef%xw3 * vmax(np,ns)
          enddo
       enddo

        if(diag==1.and.np==outp) then
           print *, 'Vmaxt',micinput%tavg(outp,1),vmax(outp,1)
           print *, 'V1=',micparam%V1(outp,1)
           print *, 'V2=',micparam%V2(outp,1)
           print *, 'V3=',micparam%V3(outp,1)
           print *, 'W1=',micparam%W1(outp,1)
           print *, 'W2=',micparam%W2(outp,1)
           print *, 'W3=',micparam%W3(outp,1)
        endif

      deallocate(vmax)
      deallocate(sdepthz)

    end subroutine Vmaxt

    subroutine Desorpt(micpxdef,micpdef,micparam,micinput)
      use mic_constant
      use mic_variable
      implicit none
!      real(r_2)              xdesorp
      TYPE(mic_param_xscale),  INTENT(IN)     :: micpxdef
      TYPE(mic_param_default), INTENT(IN)     :: micpdef
      TYPE(mic_parameter), INTENT(INOUT)      :: micparam
      TYPE(mic_input), INTENT(IN)             :: micinput
      integer :: nopt,np,ns

     do np=1,mp
      do ns=1,ms
         nopt=micparam%bgctype(np)
         micparam%desorp(np,ns) = micpxdef%xdesorp(nopt) * (1.5e-5) * exp(-1.5*micinput%clay(np,ns))
      enddo
     enddo

      if(diag==1.and. np==outp) then
         print *, 'Desorpt'
         print *, 'desorpt=',micparam%desorp(outp,:)
      endif

    end subroutine Desorpt

  subroutine mget(micpdef,micparam,micinput,micnpool)
     use mic_constant
     use mic_variable
     implicit none
     TYPE(mic_param_default), INTENT(IN)     :: micpdef
     TYPE(mic_parameter),     INTENT(INOUT)  :: micparam
     TYPE(mic_input),         INTENT(IN)     :: micinput
     TYPE(mic_npool),         INTENT(IN)     :: micnpool

     ! local variables
     integer :: np,ns

      do np=1,mp
       do ns=1,ms
          ! variable mge

         !  micparam%mgeR1(np,ns) = micpdef%cuemax*min(1.0,(micparam%cn_r(np,ns,1)/micparam%cn_r(np,ns,3)) &
         !                          **(micpdef%cue_coef1*(micnpool%mineralN(np,ns)-micpdef%cue_coef2)))

         !  micparam%mgeR2(np,ns) = micpdef%cuemax*min(1.0,(micparam%cn_r(np,ns,2)/micparam%cn_r(np,ns,3)) &
         !                          **(micpdef%cue_coef1*(micnpool%mineralN(np,ns)-micpdef%cue_coef2)))

         !  micparam%mgeR3(np,ns) = micpdef%cuemax*min(1.0,(micparam%cn_r(np,ns,7)/micparam%cn_r(np,ns,3)) &
         !                          **(micpdef%cue_coef1*(micnpool%mineralN(np,ns)-micpdef%cue_coef2)))

         !  micparam%mgeK1(np,ns) = micpdef%cuemax*min(1.0,(micparam%cn_r(np,ns,1)/micparam%cn_r(np,ns,4)) &
         !                          **(micpdef%cue_coef1*(micnpool%mineralN(np,ns)-micpdef%cue_coef2)))

         !  micparam%mgeK2(np,ns) = micpdef%cuemax*min(1.0,(micparam%cn_r(np,ns,2)/micparam%cn_r(np,ns,4)) &
         !                          **(micpdef%cue_coef1*(micnpool%mineralN(np,ns)-micpdef%cue_coef2)))

         !  micparam%mgeK3(np,ns) = micpdef%cuemax*min(1.0,(micparam%cn_r(np,ns,7)/micparam%cn_r(np,ns,4)) &
         !                          **(micpdef%cue_coef1*(micnpool%mineralN(np,ns)-micpdef%cue_coef2)))
         ! fixed mge
          micparam%mgeR1(np,ns) = micpdef%epislon1 * exp(-0.015 *micinput%tavg(np,ns))
          micparam%mgeR2(np,ns) = micpdef%epislon2 * exp(-0.015 *micinput%tavg(np,ns))
          micparam%mgeR3(np,ns) = micpdef%epislon1 * exp(-0.015 *micinput%tavg(np,ns))
          micparam%mgeK1(np,ns) = micpdef%epislon3 * exp(-0.015 *micinput%tavg(np,ns))
          micparam%mgeK2(np,ns) = micpdef%epislon4 * exp(-0.015 *micinput%tavg(np,ns))
          micparam%mgeK3(np,ns) = micpdef%epislon3 * exp(-0.015 *micinput%tavg(np,ns))
       enddo
      enddo

       if(diag==1.and.np==outp) then
          print *, 'mget'
          print *, 'epislon1-4=',micpdef%epislon1,micpdef%epislon2,micpdef%epislon3,micpdef%epislon4
       endif

  end subroutine mget

  subroutine turnovert(kinetics,micpxdef,micpdef,micparam,micinput)
      use mic_constant
      use mic_variable
      implicit none
      TYPE(mic_param_xscale),  INTENT(IN)      :: micpxdef
      TYPE(mic_param_default), INTENT(IN)      :: micpdef
      TYPE(mic_parameter),     INTENT(INOUT)   :: micparam
      TYPE(mic_input),         INTENT(IN)      :: micinput

      integer :: nx,kinetics
      real(r_2)  :: xbeta

      ! local variable
      integer :: nopt,np,ns
      real(r_2), dimension(:), allocatable    :: tvref

      allocate(tvref(mp))
       do np=1,mp
           nopt=micparam%bgctype(np)
           tvref(np) = sqrt(micinput%fcnpp(np)/micpdef%xtv)
           tvref(np) = max(0.6,min(1.3,tvref(np)))          ! 0.8-1.2 based on Wieder et al., 2015

  !         if(kinetics==3) then
  !            tvref(np) = 1.0
  !            tvref(np) = 1.0
  !         endif

           do ns=1,ms
              micparam%tvmicR(np,ns)   = micpxdef%xtvmic(nopt) * micpdef%tvmicR * tvref(np) * exp(0.3 * micparam%fmetave(np,ns)) * delt
              micparam%tvmicK(np,ns)   = micpxdef%xtvmic(nopt) * micpdef%tvmicK * tvref(np) * exp(0.1 * micparam%fmetave(np,ns)) * delt
              micparam%betamicR(np,ns) = micpdef%betamic * micpxdef%xbeta(nopt)
              micparam%betamicK(np,ns) = micpdef%betamic * micpxdef%xbeta(nopt)
           enddo
       enddo

        if(diag==1.and.np==outp) then
          print *, 'turnovert'
          print *, 'tvmicR=',micparam%tvmicR(outp,:)
          print *, 'tvmicR=',micparam%tvmicR(outp,:)
        endif
      deallocate(tvref)
  end subroutine turnovert


    subroutine bgc_fractions(micpxdef,micpdef,micparam,micinput)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_param_xscale), INTENT(IN)      :: micpxdef
    TYPE(mic_param_default), INTENT(IN)     :: micpdef
    TYPE(mic_parameter),     INTENT(INOUT)  :: micparam
    TYPE(mic_input),         INTENT(INOUT)  :: micinput
    !local variables
    integer :: npft,np,ns
    real(r_2), dimension(:),     allocatable :: fmetleaf,fmetroot,fmetwood
    real(r_2), dimension(:,:),   allocatable :: dleafx,drootx,dwoodx
    real(r_2), dimension(:,:),   allocatable :: cinputm,cinputs
    real(r_2), dimension(:,:,:), allocatable :: cninp


     allocate(fmetleaf(mp),fmetroot(mp),fmetwood(mp))
     allocate(dleafx(mp,ms),drootx(mp,ms),dwoodx(mp,ms))
     allocate(cinputm(mp,ms),cinputs(mp,ms))
     allocate(cninp(mp,ms,2))

     do np=1,mp
          npft=micparam%pft(np)
          fmetleaf(np) = max(0.0, 1.0 * (0.85 - 0.013 * micparam%fligleaf(np) * micparam%xcnleaf(np)))
          fmetroot(np) = max(0.0, 1.0 * (0.85 - 0.013 * micparam%fligroot(np) * micparam%xcnroot(np)))
          fmetwood(np) = max(0.0, 1.0 * (0.85 - 0.013 * micparam%fligwood(np) * micparam%xcnwood(np)))

          ! Initial C:N ratio of each C pool
          do ns=1,ms

             ! **this is a temporary solution, to be modified after N cycle is included
             micparam%cn_r(np,ns,1) = max( 5.0,0.5*(micparam%xcnleaf(np)+micparam%xcnroot(np)))
             micparam%cn_r(np,ns,2) = max(10.0,0.5*micparam%xcnleaf(np))
             micparam%cn_r(np,ns,3) =  7.4
             micparam%cn_r(np,ns,4) = 13.4
             micparam%cn_r(np,ns,5) = 12.0
             micparam%cn_r(np,ns,6) = 16.0
             micparam%cn_r(np,ns,7) = 10.0


       ! here zse in m, litter input in g/m2/delt, *0.001 to mgc/cm3/delt and "zse" in m.
             if(ns==1) then
                dleafx(np,ns) = micpxdef%xNPP(npft) * 0.001* micinput%dleaf(np)/micparam%sdepth(np,1)                               ! mgc/cm3/delt
                drootx(np,ns) = micpxdef%xNPP(npft) * 0.001* micparam%fracroot(np,1) * micinput%droot(np)/micparam%sdepth(np,1)     ! mgc/cm3/delt
                dwoodx(np,ns) = micpxdef%xNPP(npft) * 0.001* micinput%dwood(np)/micparam%sdepth(np,1)                               ! mgc/cm3/delt
             else
                dleafx(np,ns) = 0.0
                drootx(np,ns) = micpxdef%xNPP(npft) * 0.001 * micparam%fracroot(np,ns) * micinput%droot(np)/micparam%sdepth(np,ns)  ! mgc/cm3/delt
                dwoodx(np,ns) = 0.0
             endif

          !! calculate soil texture and litter quaility dependent parameter values
          ! C input to metabolic litter
             micinput%cinputm(np,ns) = dleafx(np,ns)*fmetleaf(np)        &
                                     + drootx(np,ns)*fmetroot(np)        &
                                     + dwoodx(np,ns)*fmetwood(np)
          ! C input to structural litter
             micinput%cinputs(np,ns) = dleafx(np,ns)*(1.0-fmetleaf(np))  &
                                     + drootx(np,ns)*(1.0-fmetroot(np))  &
                                     + dwoodx(np,ns)*(1.0-fmetwood(np))

          ! if((dleafx(np,ns)+drootx(np,ns))>0.0) then
          ! C:N input of litter input to the metabolic pool
                cninp(np,ns,1) = micinput%cinputm(np,ns)                          &
                               /(dleafx(np,ns)*fmetleaf(np)/micparam%xcnleaf(np)  &
                               +drootx(np,ns)*fmetroot(np)/micparam%xcnroot(np)   &
                               +dwoodx(np,ns)*fmetwood(np)/micparam%xcnwood(np))
          ! C:N input of litter input to the structural pool
                cninp(np,ns,2) = micinput%cinputs(np,ns)                               &
                               /(dleafx(np,ns)*(1.0-fmetleaf(np))/micparam%xcnleaf(np) &
                               +drootx(np,ns)*(1.0-fmetroot(np))/micparam%xcnroot(np)  &
                               +dwoodx(np,ns)*(1.0-fmetwood(np))/micparam%xcnwood(np))

                micparam%fmetave(np,ns) = (dleafx(np,ns)*fmetleaf(np) + drootx(np,ns)*fmetroot(np) + dwoodx(np,ns) * fmetwood(np))  &
                                        /(dleafx(np,ns) + drootx(np,ns) + dwoodx(np,ns) + 1.0e-10)

            !  else
            !    if(ns==1) then
            !       cninp(np,ns,1)          = micparam%xcnleaf(np)
            !       cninp(np,ns,2)          = micparam%xcnleaf(np)
            !       micparam%fmetave(np,ns) = fmetleaf(np)
            !    else
            !       cninp(np,ns,1)          = micparam%xcnroot(np)
            !       cninp(np,ns,2)          = micparam%xcnroot(np)
            !       micparam%fmetave(np,ns) = fmetroot(np)
            !    endif
            !  endif

             micparam%cn_r(np,ns,1) = cninp(np,ns,1); micparam%cn_r(np,ns,2)=cninp(np,ns,2)

            ! micparam%fr2p(np,ns) = micpdef%fmicsom1 * 0.30 * exp(1.3*micinput%clay(np,ns)) *1.0                   ! 3.0
            ! micparam%fk2p(np,ns) = micpdef%fmicsom2 * 0.20 * exp(0.8*micinput%clay(np,ns)) *1.0                   ! 3.0
            ! micparam%fr2c(np,ns) = min(1.0-micparam%fr2p(np,ns), micpdef%fmicsom3 * 0.10 * exp(-micpdef%fmicsom5 * micparam%fmetave(np,ns))*1.0 )    ! 9.0   to invoid a negative value of fr2a  ZHC
            ! micparam%fk2c(np,ns) = min(1.0-micparam%fk2p(np,ns), micpdef%fmicsom4 * 0.30 * exp(-micpdef%fmicsom5 * micparam%fmetave(np,ns))*1.0)     ! 9.0   to invoid a negative value of fk2a ZHC
            ! micparam%fr2a(np,ns) = 1.00 - micparam%fr2p(np,ns) - micparam%fr2c(np,ns)
            ! micparam%fk2a(np,ns) = 1.00 - micparam%fk2p(np,ns) - micparam%fk2c(np,ns)
            ! changes made to accommodate added aggregated pools

             micparam%fr2p(np,ns) =  0.30 * exp(1.3*micinput%clay(np,ns))                    ! 3.0
             micparam%fk2p(np,ns) =  0.20 * exp(0.8*micinput%clay(np,ns))                    ! 3.0
             micparam%fr2c(np,ns) = min(1.0, micparam%fr2p(np,ns) + 0.10 * exp(-3.0 * micparam%fmetave(np,ns)))     ! 9.0   to invoid a negative value of fr2a  ZHC
             micparam%fk2c(np,ns) = min(1.0, micparam%fr2p(np,ns) + 0.30 * exp(-3.0 * micparam%fmetave(np,ns)))     ! 9.0   to invoid a negative value of fk2a ZHC
             micparam%fr2p(np,ns) =  0.0
             micparam%fk2p(np,ns) =  0.0
             micparam%fr2a(np,ns) = max(0.0,1.00 - micparam%fr2c(np,ns))
             micparam%fk2a(np,ns) = max(0.0,1.00 - micparam%fk2c(np,ns))
          enddo   !"ns"
     enddo       !"np"

      if(diag==1.and.np ==outp) then
         print *,'bgc_fraction parameters'
         print *, 'empirical params1-4=', micpdef%fmicsom1,micpdef%fmicsom2,micpdef%fmicsom3,micpdef%fmicsom4
         print *, 'clay=', micinput%clay(np,:)
         print *, 'cinputm=', micinput%cinputm(outp,:)
         print *, 'cinputs=',micinput%cinputs(outp,:)
         print *, 'fmetave=',micparam%fmetave(outp,:)
         print *, 'cn_r1=',micparam%cn_r(outp,:,1)
         print *, 'cn_r2=',micparam%cn_r(outp,:,2)
         print *, 'fr2p=',micparam%fr2p(outp,:)
         print *, 'fk2p=',micparam%fk2p(outp,:)
         print *, 'fr2c=',micparam%fr2c(outp,:)
         print *, 'fk2c=',micparam%fk2c(outp,:)
         print *, 'fr2a=',micparam%fr2a(outp,:)
         print *, 'fk2a=',micparam%fk2a(outp,:)
      endif
     deallocate(fmetleaf,fmetroot,fmetwood)
     deallocate(dleafx,drootx,dwoodx)
     deallocate(cinputm,cinputs)
     deallocate(cninp)

   end subroutine bgc_fractions

   subroutine bioturb(ndelt,ms,zse,delt,diffsocxx,fluxsoc,xpooli,xpoole)
   ! multi-layered soil BGC including DOC and bioturbation using microbially-based BGC modeling
   ! step 1: litter-C and SOC bioturbation treated as a diffusion process
   ! step 2: advection of DOC along with water flux
   ! solve dc/dt=Dd2c/dx +F(z) where c is total SOC concentration in each soil layer
   ! bioturbation diffusion rate
   ! boundary conditions: at the top     -Ddc/dx = F0+F(1)  at x=0
   !                      at the bottom: dC/dx=0            at x=h
   ! using the fully implicit method together with Thomas method
   ! unit for pool:                 mgc/cm3      (=kg C/m3)
   !      for flux:                 mgc/cm3/delt (=kg c/m3/delt): g/m2/delt = 0.1 mg/cm2/delt
   !      for length:               cm
   !      for diffsion coefficient: cm2/delt
   use mic_constant,  ONLY : r_2
   implicit none
   integer                        :: ndelt,ms
   real(r_2), dimension(ms)    :: zse
   real(r_2)                      :: delt,diffsocxx
   real(r_2), dimension(ms)    :: xpooli,xpoole,xpool,fluxsoc
   ! local variables
   integer                        :: i,j
   real(r_2)                      :: deltD,tot0, tot1, totflux
   real(r_2), dimension(ms)    :: xzse
   real(r_2), dimension(ms+1)  :: sdepthx
   real(r_2)                      :: coeffA, coeffB
   real(r_2), dimension(ms)    :: at,bt,ct,rt

      ! calculate the mid-point of each layer
     sdepthx(1) = 0.0          ! depth of a layer from the top (x_0.5=0.0 eg soil surface)
     do j=2,ms+1
        sdepthx(j) = sdepthx(j-1) + zse(j-1)*100.0     ! depth of the bottom of each layer (eg x_j+0.5)
                                                       !*100 to convert from m to cm
     enddo

     do j=1,ms
        xzse(j) = 0.5 * (sdepthx(j) + sdepthx(j+1))    ! depth of midpoint of a layer j  (x_j)
     enddo

     deltD = diffsocxx * delt

      xpool = xpooli
      tot0 = 0.0
     do j=1,ms
        tot0 = tot0 + xpool(j) * zse(j)*100.0         !*100 convert m to cm
     enddo

     do i=1,ndelt
        do j=1,ms
           if(j==1) then
              coeffB = 1.0/(sdepthx(2)-sdepthx(1))
              coeffA = deltD*coeffB/(xzse(2)-xzse(1))
              ! Crank-Nicholson
              at(1) = 0.0
              bt(1) = 1.0 + 0.5 * coeffA
              ct(1) =     - 0.5 * coeffA
              rt(1) = (1.0-0.5*coeffA) * xpool(1) + 0.5 * coeffA * xpool(2) &
                    +  fluxsoc(1) * delt
           endif
           if(j>1.and.j<ms) then
             coeffA = deltD/((xzse(j+1)-xzse(j))*(sdepthx(j+1)-sdepthx(j)))
             coeffB = (xzse(j+1)-xzse(j))/(xzse(j)-xzse(j-1))
             ! Crank-Nicholson
             at(j) =    -0.5 * coeffA * coeffB
             bt(j) = 1.0+0.5 * coeffA *(1.0+coeffB)
             ct(j) =    -0.5 * coeffA
             rt(j) = 0.5 * coeffA * coeffB * xpool(j-1)        &
                     +(1.0-0.5* coeffA*(1.0+coeffB))*xpool(j)  &
                     + 0.5* coeffA * xpool(j+1)                &
                     + fluxsoc(j) *delt
           endif
           if(j==ms) then
               coeffA = deltD/((xzse(ms)-xzse(ms-1))*(sdepthx(ms+1) - sdepthx(ms)))
             ! Crank-Nicholson
               at(ms) = -0.5 * coeffA
               bt(ms) = 1.0 + 0.5 * coeffA
               ct(ms) = 0.0
               rt(ms) = 0.5* coeffA  * xpool(ms-1) + (1.0-0.5 * coeffA) * xpool(ms) &
                    + fluxsoc(ms) * delt
            endif
        enddo

        call tridag(at,bt,ct,rt,xpool,ms)
     enddo
     xpoole = xpool

     tot1 = 0.0
     totflux=0.0
     do j=1,ms
        tot1 = tot1 + xpool(j) * zse(j) *100.0
        totflux = totflux + fluxsoc(j) * zse(j) *100.0
     enddo

end subroutine bioturb

   subroutine tridag(at,bt,ct,rt,u,ms)
   ! solving the triadigonal matrix (numerical recipes, p43)
   ! linear equation: A* u(i-1) + B *u(i) + C * u(i+1) = R,
   ! where i is soil layer, u(i-1), u(i) and u(i+1) are at time step t
   ! NOTE: bt(1) should not be equal to 0.0, otherwise rewrite the equation
    use mic_constant,  ONLY : r_2
    implicit none
    integer, parameter    :: nmax=500
    integer :: ms
    real(r_2), dimension(ms)    :: at,bt,ct,rt,u
    integer :: j
    real(r_2) :: bet
    real(r_2), dimension(nmax) :: gam

      bet  = bt(1)
      u(1) = rt(1)/bet
      do j=2,ms
         gam(j) = ct(j-1)/bet
         bet = bt(j)-at(j)* gam(j)
         if(bet ==0) then
            print *, 'triag failed'
            stop
         endif
         u(j) = (rt(j) - at(j) * u(j-1))/bet
      enddo
      do j=ms-1,1,-1
         u(j) = u(j) -gam(j+1) * u(j+1)
      enddo
    end subroutine tridag


    subroutine advecdoc(deltx,fluxsoilwx,fluxdocsx,vsoilwx,ypool)
    ! to be modified using an implicit solver to ensure mass conservation
    !
    use mic_constant
    implicit none
    real(r_2)                          :: deltx
    real(r_2), dimension(ms)        :: fluxsoilwx,vsoilwx,ypool
    real(r_2), dimension(ms)        :: dypool,ypool1
    real(r_2)                          :: fluxdocsx,totdoc0,totdoc1,fluxdocbot
    integer :: ns,iter

     ypool1= ypool
     fluxdocbot = 0.0
     do iter=1,100
      do ns=1,ms
        if(ns==1) then
           dypool(1)  = (fluxdocsx - fluxsoilwx(1)*ypool1(1)/vsoilwx(1))*deltx*0.01/zse(1)
        else
           dypool(ns) = (fluxsoilwx(ns-1)*ypool1(ns-1)/vsoilwx(ns-1) &
                        -fluxsoilwx(ns)  *ypool1(ns)/vsoilwx(ns))       *deltx*0.01/zse(ns)
        endif
        if(ns==ms) then
           fluxdocbot = fluxdocbot + fluxsoilwx(ns)  *ypool1(ns)/vsoilwx(ns) *deltx* 0.01
        endif
      enddo
      ypool1 = max(0.0,ypool1+dypool)
     enddo
     ! check mass conservation
     totdoc0=0.0; totdoc1=0.0
     do ns=1,ms
        totdoc0 = totdoc0 + ypool(ns)  *zse(ns)
        totdoc1 = totdoc1 + ypool1(ns) *zse(ns)
     enddo
    ! print *, 'mass cons DOC', totdoc0,totdoc1,(totdoc1-totdoc0)-(fluxdocsx - fluxdocbot)*deltx

     ypool = ypool1

    end subroutine advecdoc


    subroutine vmic_c(ny,isoc14,np,ns,kinetics,micpdef,micparam,micinput,xpool,y)
    ! MIMICS as modified by Zhang et al. (2019, GCB).
    ! Seven pools: metabolic litter (1), Structural litter, microbe-R (3), microbe-K(4),
    !              Physical protected (5), chemically-protected (6), active (7)
    ! for each layer
    ! input: aboveground litter and belowground in the surface layer
    !        belowground input in other layers
    ! kinetics: Michaelis-Mennten
    ! unit:
    ! all carbon pools : mg C/cm3
    ! time step :        one hour
    !
     use mic_constant
     use mic_variable
     implicit none

     TYPE(mic_param_default), INTENT(IN)     :: micpdef
     TYPE(mic_parameter),     INTENT(IN)     :: micparam
     TYPE(mic_input),         INTENT(IN)     :: micinput

     real(r_2),  parameter                           ::  kamin = 0.2      !          Abramoff et al. (2022)
     real(r_2),  parameter                           ::  lamda = 2.01e-4  !1/kPa     Abramoof et al. (2022)

     real(r_2),  dimension(mcpool),  INTENT(IN)        :: xpool
     real(r_2),  dimension(mcpool),  INTENT(INOUT)     :: y   !=dxpool/dt     ! local variables
     ! local variables
     integer     :: np,ns,kinetics,ny,isoc14,ip

     real(r_2)  :: betamicR,betamicK,                 &
                cinputmx,cinputsx,fmx,fsx,         &
                fr2px,fr2cx,fr2ax,                 &
                fk2px,fk2cx,fk2ax,                 &
                mgeRx1,mgeRx2,mgeRx3,              &
                mgeKx1,mgeKx2,mgeKx3,              &
                tvmicRx,tvmicKx,                   &
                tavgx,clayx,                       &
                desorpx,                           &
                V1x,V2x,V3x,W1x,W2x,W3x,           &
                J1x,J2x,J3x,K1x,K2x,K3x,           &
                Q1x,Q2x


     real(r_2) :: cfluxm2r, cfluxm2k, cfluxs2r, cfluxs2k, cfluxr,   cfluxk
     real(r_2) :: cfluxr2p, cfluxk2p, cfluxp2a, cfluxr2c, cfluxk2c
     real(r_2) :: cfluxc2a, cfluxr2a, cfluxk2a, cfluxa2r, cfluxa2k

     ! additional variables for kinetics3
     real(r_2)  :: cfluxa,cfluxp,cfluxc2p,cfluxa2c,cfluxp2c
     real(r_2)  :: kadsorpx,kdesorpx,fp2ax,moistx,soilphx,porex,xwater,phx1,phx2,phx3,siltx,tvcpoolx,tvppoolx,tvacx
     real(r_2)  :: smexpa,smopt,qmaxcoeffx,qmax
     real(r_2)  :: swbx,swdx,matpotx,xwater1,xwater2
     real(r_2)  :: rsoil
     real(r_2)  :: cfluxp2m,cfluxp2s

      ! matpotx =-15.0 !dummy value for the time being

      if(isoc14==1) then

         if(ny<1) then
            cinputmx = micinput%cinputm(np,ns) *  micparam%c14atm(1,micparam%region(np),2)        ! using fraction modern before 1941
            cinputsx = micinput%cinputs(np,ns) *  micparam%c14atm(1,micparam%region(np),2)
         else
            cinputmx = micinput%cinputm(np,ns) *  micparam%c14atm(ny,micparam%region(np),2)       ! using fraction modern after 1941
            cinputsx = micinput%cinputs(np,ns) *  micparam%c14atm(ny,micparam%region(np),2)
         endif
      else
         cinputmx = micinput%cinputm(np,ns)
         cinputsx = micinput%cinputs(np,ns)
      endif

      tavgx    = micinput%tavg(np,ns);      clayx    = micinput%clay(np,ns);    siltx    = micinput%silt(np,ns)

      fmx      = micparam%fm(np,ns);        fsx      = micparam%fs(np,ns)
      fr2px    = micparam%fr2p(np,ns);      fr2cx    = micparam%fr2c(np,ns)
      fr2ax    = micparam%fr2a(np,ns);      fk2px    = micparam%fk2p(np,ns)
      fk2cx    = micparam%fk2c(np,ns);      fk2ax    = micparam%fk2a(np,ns)
      mgeRx1   = micparam%mgeR1(np,ns);     mgeRx2   = micparam%mgeR2(np,ns);   mgeRx3 = micparam%mgeR3(np,ns)
      mgeKx1   = micparam%mgeK1(np,ns);     mgeKx2   = micparam%mgeK2(np,ns);   mgeKx3 = micparam%mgeK3(np,ns)
      tvmicRx  = micparam%tvmicR(np,ns);    tvmicKx  = micparam%tvmicK(np,ns)
      desorpx  = micparam%desorp(np,ns)
      V1x      = micparam%V1(np,ns);        V2x      = micparam%V2(np,ns);      V3x    = micparam%V3(np,ns)
      W1x      = micparam%W1(np,ns);        W2x      = micparam%W2(np,ns);      W3x    = micparam%W3(np,ns)
      K1x      = micparam%K1(np,ns);        K2x      = micparam%K2(np,ns);      K3x    = micparam%K3(np,ns)
      J1x      = micparam%J1(np,ns);        J2x      = micparam%J2(np,ns);      J3x    = micparam%J3(np,ns)
      Q1x      = micparam%Q1(np,ns);        Q2x      = micparam%Q2(np,ns)
      betamicR = micparam%betamicR(np,ns);  betamicK = micparam%betamicK(np,ns)


!      print *, 'p1=',  fmx, fsx
!      print *, 'p2=',  fr2px,fr2cx
!      print *, 'p3=',  fr2ax,fk2px
!      print *, 'p4=',  fk2cx,fk2ax
!      print *, 'p5=',  mgeRx1,mgeRx2,mgeRx3
!      print *, 'p6=',  mgeKx1,mgeKx2,mgeKx3
!      print *, 'p7=',  tvmicRx,tvmicKx
!      print *, 'p8=',  desorpx
!      print *, 'p9=',  V1x,V2x,V3x
!      print *, 'p10=', W1x,W2x,W3x
!      print *, 'p11=', K1x,K2x,K3x
!      print *, 'p12=', J1x,J2x,J3x
!      print *, 'p13=', Q1x,Q2x
!      print *, 'p14=', betamicR,betamicK

      ! additional parameters and input for kinetics3
      if(kinetics==3) then
         moistx   = micinput%wavg(np,ns);          matpotx    = micinput%matpot(np,ns)
         soilphx  = micinput%ph(np,ns);            porex      = micinput%porosity(np,ns)
         kadsorpx = micparam%kadsorp(np,ns);       tvcpoolx   = micparam%tvcpool(np,ns);   tvppoolx =micparam%tvppool(np,ns)
         tvacx    = micparam%tvac(np,ns);          fp2ax      = micparam%fp2a(np,ns)
         kdesorpx = micparam%kdesorp(np,ns);       qmaxcoeffx = micparam%qmaxcoeff(np,ns)

         ! we applied a single water-limiting function from Yan et al. 2018, Nature Coomunitation. eqn1)
         if(clayx <= 0.016) then
            smexpa=0.0
         else if(clayx >0.016 .and. clayx <=0.37) then
            smexpa=2.8* clayx-0.046
         else
            smexpa=1.0
         endif

         smopt     = 0.65 * porex

         if(moistx < smopt) then
            xwater = ((micpdef%smkdesorp+smopt)/(micpdef%smkdesorp + moistx))     &
                   * (moistx/smopt)**(1.0+smexpa *micpdef%smexpns)

         else
            xwater = ((porex - moistx)/(porex-smopt)) **micpdef%smexpb
         endif


         ! soil water limitation from Abramoff et al. (2022) eqn (4) and eqn (15)
         xwater1 = sqrt(moistx/porex)                                                                ! eqn (4)
         xwater2 = exp(lamda * matpotx) *(kamin + (1.0- kamin) *sqrt(1.0-moistx/porex))*xwater1      ! eqn (15)

         phx1      = exp(-micpdef%phcoeff1* soilphx - micpdef%phcoeff2)                          ! eqn(10) Abramoff2022
         phx2     = 1.0/(1.0+exp(-(soilphx-4.798)/0.4246))        !bacteria
         phx3     = 1.0/(1.0+exp(-(soilphx-3.022)/0.428))         !fungi
         qmax     = qmaxcoeffx * (clayx + siltx)*100.0  ! mg C/g soil, based on Georgiou et al. (2022), media value (their  Fig 1)
         ! unit conversion qmax to mg C/cm3
         qmax     = qmax * micinput%bulkd(np,ns) *0.001                        ! bulkd in kg/m3 multiply by 0.001 into g/cm3

         !    xwater = sqrt(micinput%moist(np,ns)/micinput%poros(np,ns))      ! eqn (4) Abramoff2022
         !    xwmic    = xwdecomp * exp(lambda * (-micinput%matpot(np,ns)) * (swmin + (1.0-swmin(np,ns)) * &
         !               ((micinput%poros(np) - micinput%moist(np,ns))/micinput%poros(np,ns)) **0.5)    !eqn(15) Abramoff2022

      endif
      ! carbon fluxes
      if(kinetics==1) then
        ! forward Michaelis-Menten
        cfluxm2r = xpool(3) * V1x * xpool(1)/(K1x + xpool(1))
        cfluxs2r = xpool(3) * V2x * xpool(2)/(K2x + xpool(2))
        cfluxa2r = xpool(3) * V3x * xpool(7)/(K3x + xpool(7))

        cfluxm2k = xpool(4) * W1x * xpool(1)/(J1x + xpool(1))
        cfluxs2k = xpool(4) * W2x * xpool(2)/(J2x + xpool(2))
        cfluxa2k = xpool(4) * W3x * xpool(7)/(J3x + xpool(7))

        cfluxr   = tvmicRx * xpool(3) ** betamicR
        cfluxk   = tvmicKx * xpool(4) ** betamicK

        cfluxr2p = fr2px * cfluxr
        cfluxk2p = fk2px * cfluxk

        cfluxr2c = fr2cx   * cfluxr
        cfluxk2c = fk2cx   * cfluxk

        cfluxp2a = desorpx * xpool(5)
        cfluxr2a = fr2ax   * cfluxr
        cfluxk2a = fk2ax   * cfluxk
        cfluxc2a = xpool(3)* V2x * xpool(6)/(Q1x*K2x + xpool(6))   &
                 + xpool(4)* W2x * xpool(6)/(Q2x*J2x + xpool(6))
      endif
      if(kinetics ==2 )then
        !=======================================================
        ! reverse Michaelis-Menten
        cfluxm2r = xpool(1) * V1x * xpool(3)/(K1x + xpool(3))
        cfluxs2r = xpool(2) * V2x * xpool(3)/(K2x + xpool(3))
        cfluxa2r = xpool(7) * V3x * xpool(3)/(K3x + xpool(3))

        cfluxm2k = xpool(1) * W1x * xpool(4)/(J1x + xpool(4))
        cfluxs2k = xpool(2) * W2x * xpool(4)/(J2x + xpool(4))
        cfluxa2k = xpool(7) * W3x * xpool(4)/(J3x + xpool(4))

        cfluxr   = tvmicRx * xpool(3) ** betamicR
        cfluxk   = tvmicKx * xpool(4) ** betamicK


        cfluxr2p = fr2px   * cfluxr
        cfluxk2p = fk2px   * cfluxk

        cfluxr2c = fr2cx * cfluxr
        cfluxk2c = fk2cx * cfluxk

        cfluxp2a = desorpx * xpool(5)
        cfluxr2a = fr2ax * cfluxr
        cfluxk2a = fk2ax * cfluxk
        cfluxc2a = xpool(6) * V2x * xpool(3)/(Q1x*K2x + xpool(3))   &
                 + xpool(6) * W2x * xpool(4)/(Q2x*J2x + xpool(4))
      endif

      !===================================================
      !
      if(kinetics ==1 .or. kinetics==2) then
         ! metabolic litter  [=Im*(1-fm)-A1-A5]
         y(1) = cinputmx * (1.0-fmx) - cfluxm2r - cfluxm2k

         ! structural litter [=Is*(1-fs)-A2-A6]
         y(2) = cinputsx * (1.0-fsx) - cfluxs2r - cfluxs2k

        ! these two are incorrect
        ! !microbe R          [mge1*A1+mge2*A2+mge3*A3-A4]
        ! y(3) = mgeRx1 * cfluxm2r + mgeRx2 * cfluxs2r + mgeRx3 * cfluxa2r - cfluxr

        ! !microbe K          [mge3*A5+mge4*A6+mge3*A7-A8]
        ! y(4) = mgeKx1 * cfluxm2k + mgeKx2 * cfluxs2k + mgeKx2 * cfluxa2k - cfluxk

        ! !microbe R          [mge1*A1+mge2*A2+mge1*A3-A4]
         y(3) = mgeRx1 * cfluxm2r + mgeRx2 * cfluxs2r + mgeRx1 * cfluxa2r - cfluxr

        ! !microbe K          [mge3*A5+mge4*A6+mge3*A7-A8]
         y(4) = mgeKx1 * cfluxm2k + mgeKx2 * cfluxs2k + mgeKx1 * cfluxa2k - cfluxk

         !physically protected SOM: [Lm*fm+fpr*A4+fpk*A8-A9]
         y(5) = cinputmx * fmx + cfluxr2p + cfluxk2p - cfluxp2a

         ! chemically protected SOM: [Is*fs+fcr*A4+fck*A8-A10]
         y(6) = cinputsx * fsx + cfluxr2c + cfluxk2c - cfluxc2a

         !active SOM: [far*A4+fak*A8+A9+A10-A3-A7]
         y(7) = cfluxr2a + cfluxk2a + cfluxp2a + cfluxc2a - cfluxa2r - cfluxa2k
         ! additional dummy pools
         y(8) = 0.0
         y(9) = 0.0
         y(10)= 0.0

      endif

      ! the new soil carbon model combining MIMICS and MILLENNIAL2
      ! we use two litter pools (m,s) and two microbial pool (r,k) and LWC (pool a), aggregate C (pool p) amd MAOC (pool C)
      ! see documentation on the combined model
      if(kinetics==3) then
        ! reverse Michaelis-Menten for litter and forward MM for pool 7
        cfluxm2r = xpool(1) * V1x * phx2 * xwater2 * xpool(3)/(K1x + xpool(3))    ! eqn 2 Abramoff2022
        cfluxs2r = xpool(2) * V2x * phx2 * xwater2 * xpool(3)/(K2x + xpool(3))    ! eqn 2 Abramoff2022
        cfluxa2r = xpool(3) * V3x * phx2 * xwater2 * xpool(7)/(K3x + xpool(7))
!        cfluxa2r = xpool(7) * V3x * xwater2 * xpool(3)/(K3x + xpool(3))    ! eqn 2 Abramoff2022

        cfluxm2k = xpool(1) * W1x * phx3 * xwater2 * xpool(4)/(J1x + xpool(4))    ! eqn 2 Abramoff2022
        cfluxs2k = xpool(2) * W2x * phx3 * xwater2 * xpool(4)/(J2x + xpool(4))    ! eqn 2 Abramoff2022
        cfluxa2k = xpool(4) * W3x * phx3 * xwater2 * xpool(7)/(J3x + xpool(7))
!        cfluxa2k = xpool(7) * W3x * xwater * xpool(4)/(J3x + xpool(4))    ! eqn 2 Abramoff2022
       ! forward Michaelis-Menten
!        cfluxm2r = xpool(3) * V1x * xpool(1)/(K1x + xpool(1))
!        cfluxs2r = xpool(3) * V2x * xpool(2)/(K2x + xpool(2))
!        cfluxa2r = xpool(3) * V3x * xpool(7)/(K3x + xpool(7))

!        cfluxm2k = xpool(4) * W1x * xpool(1)/(J1x + xpool(1))
!        cfluxs2k = xpool(4) * W2x * xpool(2)/(J2x + xpool(2))
!        cfluxa2k = xpool(4) * W3x * xpool(7)/(J3x + xpool(7))

        cfluxa   = tvacx    * xwater1 * xpool(7)                           ! eqn 8 Abramoff2022 (leaching)
        cfluxa   = 0.0                                                    ! labile C leaching is done separately
        cfluxr   = tvmicRx  * xpool(3) ** betamicR                        ! eqv. eqn(16) Abramoff2022
        cfluxk   = tvmicKx  * xpool(4) ** betamicK                        ! eqv. eqn(16) Abramoff2022
!        cfluxp   = tvppoolx * xwater1 * xpool(5)                          ! eqv. eqn (6) abramoff2022 (aggregate breakdown)
        cfluxc2p = tvcpoolx * xwater1 * xpool(6)                          ! eqn(18) Abramoff2022  , all flux to aggregate C

        ! flux to MAOC (pool c)
        cfluxa2c = kadsorpx * phx1 * xwater1 * xpool(7) * (1.0 -  xpool(6)/qmax)   ! eqn(9)  abramoff2022
        cfluxr2c = fr2cx * cfluxr                                                ! p_b*F_bm in eqn(19) Abramoff2022
        cfluxk2c = fk2cx * cfluxk                                                ! p_b*F_bm in eqn(19) Abramoff2022
!        cfluxp2c = (1.0-fp2ax) * cfluxp                                          !(1-p_a)*F_a  in eqn(19) of Abramoff2022

        ! flux to low weight mass C (pool a)
        cfluxr2a = (1.0-fr2cx) * cfluxr                                          !(1-p_b)*F_bm in eqn(7) Abramoff2022
        cfluxk2a = (1.0-fk2cx) * cfluxk                                          !(1-p_b)*F_bm in eqn(7) Abramoff2022
!        cfluxp2a = fp2ax * cfluxp                                                ! p_a*F_a  different from Abramoff2022, Aggregate -> active, not litter
!        ! this equation is wrong  (see Wang et al. 2022, their GCB papes, Tables S3 and S5)
!        cfluxc2a = kadsorpx * xpool(6)/qmax                                      ! eqn(12) Abramoff2022
        ! based on Wang et al. (2022)
        cfluxc2a = kdesorpx * xpool(6)/qmax

!       disaggregation fluxes
        cfluxp2m = tvppoolx * xwater1 * xpool(5)   !metabolic pool
        cfluxp2s = tvppoolx * xwater1 * xpool(8)   !structural litter pool
        cfluxp2c = tvppoolx * xwater1 * xpool(9)   !MAOC  assuming disaggregation of aggregagated MAOC is slower than POC


         y(1) = cinputmx * (1.0-fmx) + cfluxp2m - cfluxm2r - cfluxm2k                       ! same as kinetics=2

         ! structural litter [=Is*(1-fs)-A2-A6]
         y(2) = cinputsx * (1.0-fsx) + cfluxp2s - cfluxs2r - cfluxs2k                       ! same as kinetics=2

         !microbe R          [mge1*A1+mge2*A2+mge3*A3-A4]
         y(3) = mgeRx1 * cfluxm2r + mgeRx2 * cfluxs2r + mgeRx1 * cfluxa2r - cfluxr ! same as kinetics=2

         !microbe K          [mge3*A5+mge4*A6+mge3*A7-A8]
         y(4) = mgeKx1 * cfluxm2k + mgeKx2 * cfluxs2k + mgeKx1 * cfluxa2k - cfluxk ! same as kinetics =2

         !Aggregate metabolic C (pool p)
         y(5) = cinputmx * fmx   - cfluxp2m

         !MAOC (pool c)
         y(6) = cfluxr2c + cfluxk2c + cfluxp2c + cfluxa2c  - cfluxc2a - cfluxc2p    ! eqn(19) of abramoff2022
                                                                                    ! cfluxa2c<->F_lm (sorption);   cfluxc2a<->F_ld (desorption)
                                                                                    ! cinputsx * fsx<-> no match;   (cfluxr2c + cfluxk2c)<->p_b * F_bm
                                                                                    ! cfluxp2c<->(1-p_a)*F_a;        cfluxc2p<->F_ma
                                                                                    !
         !LWC
         y(7) = cfluxr2a + cfluxk2a + cfluxc2a - cfluxa - cfluxa2r - cfluxa2k - cfluxa2c              ! eqn(7) Abramoff2022
                                                                                    ! no litter litter input. None <-> F_i
                                                                                    ! cfluxa: F_l (leaching)
                                                                                    ! no depolymeration <->F_pl :: litter input does not enter this pool directly
                                                                                    ! cfluxa2c<->F_lm (adsorption)
                                                                                    ! (cfluxa2r + cfluxa2k)<-> F_lb
                                                                                    ! (cfluxr2a + cfluxk2a)<->(1-p_b)*F_bm,  necromass input
                                                                                    ! cfluxc2a<->F_ld, (desorption)
                                                                                    ! clfuxp2a: de-aggregation       ! different from Abramoff2022
        ! aggregated structural C
        y(8) = cinputsx * fsx - cfluxp2s
        ! aggregated MAOC
        y(9) = cfluxc2p  - cfluxp2c
        ! dummy pool
        y(10) = 0.0
        ! check mass balance
        rsoil = (1.0-mgeRx1) * (cfluxm2r+cfluxa2r) + (1.0-mgeKx1)* (cfluxm2k + cfluxa2k)  &
              + (1.0-mgeRx2) * cfluxs2r            + (1.0-mgeKx2)* cfluxs2k  - cfluxa

!       write(*,101) np,ns, cinputmx+cinputsx,sum(y(1:7)),rsoil, cinputmx+cinputsx-sum(y(1:7))-rsoil
101 format('vmic_c: input, sumdelc rsoil',2(i3,1x),10(f10.6,1x))
      endif

!      print *, ' @ vmic_c xpool =', xpool(:)
!      print *, ' @ vmic_c y =', y(:)

      if(isoc14==1) then

         do ip=1,mcpool
            y(ip) = y(ip) - tvc14 * max(0.0,xpool(ip))
         enddo
      endif

   end subroutine vmic_c

