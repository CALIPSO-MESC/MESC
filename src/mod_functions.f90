!> Per-mode orchestrator functions for MESC calibration.
!>
!> Each function allocates arrays, ingests observed data, runs the soil
!> carbon model via an appropriate driver, computes a cost against
!> observations, and returns the scalar cost value.  The dispatcher
!> [[functn]] selects the active mode based on `case.txt`.
module function_module
  use precision_module, only: dp, r_2
  use mic_constant, only: mp, mpft, mbgc, ntime, nlon, nlat, ms
  use mic_variable, only: mic_param_xscale, mic_param_default, mic_parameter, &
                          mic_input, mic_global_input, mic_cpool, mic_npool, mic_output, &
                          mic_allocate_parameter, mic_allocate_input, mic_allocate_output, &
                          mic_allocate_cpool, mic_allocate_npool, &
                          mic_deallocate_parameter, mic_deallocate_input, mic_deallocate_output, &
                          mic_deallocate_cpool, mic_deallocate_npool
  use mesc_inout_module, only: getdata_c14, getdata_frc_dim, getdata_frc, &
                               getdata_hwsd_dim, getdata_hwsd, screenout, &
                               getparam_global,getpatch_global, &
                               getdata_global4_cable,getdata_global4_orchidee, getdata_aust_dim,getdata_aust
  use mesc_interface_module, only: vmic_param_xscale, vmic_param_time, vmic_param_time_single, vmicsoil_c14, &
                                   vmicsoil_frc1_cpu, vmicsoil_hwsd_cpu, vmicsoil_hwsd_gpu
  use calcost_module, only: calcost_c14, calcost_frc1, calcost_hwsd3, calcost_global_hwsd, calcost_aust
  implicit none

 Contains

 !> Dispatcher that selects the appropriate per-mode orchestrator.
 !!
 !! Reads `case.txt` to determine the run mode, then delegates to one
 !! of the mode-specific functions:
 !!
 !! * 0 -> [[functn_soc_aust]]
 !! * 1 -> [[functn_c14]]
 !! * 2 -> [[functn_frc1]]
 !! * 3 -> [[functn_soc_hwsd]]
 !! * 4 -> [[functn_global4]]
 real(dp) function functn(nx,xparam16)

     integer, intent(in) :: nx
         !! Number of optimized parameters.
     real(dp), dimension(16), intent(in) :: xparam16
         !! Values of the `nx` optimized parameters.

     integer :: runcase

     open(1,file='case.txt')
     read(1,*) runcase
     close(1)
     SELECT CASE (runcase)
       CASE(0)
         functn = functn_soc_aust(nx,xparam16)
       CASE (1)    ! run model for 14C
         functn = functn_c14(nx,xparam16)
       CASE (2)    ! run model for POC/MAOC fractions
         functn = functn_frc1(nx,xparam16)
       CASE (3)    ! run model for site HWSD SOC profile
         functn = functn_soc_hwsd(nx,xparam16)
       CASE (4)    ! run model globally with CABLE/ORCHIDEE spatial resolution
         functn = functn_global4(nx,xparam16)
     !  CASE (5)    ! run model with prescribed forcing (not developed yet for 1% per year offline)
     !    functn = functn_offline(nx,xparam16)
     END SELECT

 END function functn

 real(dp) function functn_c14(nx,xparam16)
    !! Orchestrator for 14C calibration mode.
    !!
    !! Runs the model twice: once for stable C (12C) and once for
    !! radiocarbon (14C), then returns the combined cost.  Reads
    !! configuration and default parameters from `params1.txt`.

    integer, intent(in) :: nx
        !! Number of optimized parameters.
    real(dp), dimension(16), intent(in) :: xparam16
        !! Values of the `nx` optimized parameters.

    TYPE(mic_param_xscale)    :: micpxdef
    TYPE(mic_param_default)   :: micpdef
    TYPE(mic_parameter)       :: micparam
    TYPE(mic_input)           :: micinput
    TYPE(mic_global_input)    :: micglobal
    TYPE(mic_cpool)           :: miccpool
    TYPE(mic_npool)           :: micnpool
    TYPE(mic_output)          :: micoutput

    integer,   dimension(16)           :: nxopt
    real(dp),    dimension(16)           :: xopt
    real(dp)     :: totcost1,totcost2
    integer    :: ifsoc14,kinetics,bgcopt,jopt,nyeqpool,isoc14,jglobal,jmodel
    integer :: jrestart,nparam
    character(len=140) :: frestart_in,frestart_out,foutput
    character(len=140) :: frac14c,f14c(5),filecluster
    real(r_2), dimension(:), allocatable :: zse

      jrestart=0;xopt(:)=1.0
      do nparam=1,16
         nxopt(nparam) = nparam
      enddo

      frestart_in='miccpool_in.nc'
      frestart_out='miccpool_out.nc'
      foutput='vmic_output.nc'

      open(1,file='params1.txt')
      read(1,*)
      read(1,*) jglobal,ifsoc14,kinetics,bgcopt,jopt,jrestart
      read(1,11) frac14c
      read(1,11) f14c(1)
      read(1,11) f14c(2)
      read(1,11) f14c(3)
      read(1,11) f14c(4)
      read(1,11) f14c(5)
      read(1,11) filecluster
11    format(a140)
      read(1,*) xopt(1:14)

      if(jopt==0) then
         read(1,*) nxopt(1:nx)
         do nparam=1,nx
            xopt(nxopt(nparam)) = xparam16(nparam)
         enddo
      endif
    !  print*, xopt

      mp = 213

      totcost1 = 0.0; totcost2=0.0
      nyeqpool= 500;jmodel=1;mpft=17;mbgc=12;ntime=1;nlon=1;nlat=1
      ms=15
      allocate(zse(ms))
      zse(1:ms)=0.1

      call mic_allocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_allocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_allocate_output(mp,micoutput)
      call mic_allocate_cpool(mp,ms,miccpool)
      call mic_allocate_npool(mp,ms,micnpool)

          isoc14 = 0
      !    print *, "isoc14 =",isoc14,'--getdata_c14'
          call getdata_c14(frac14c,f14c,filecluster,micinput,micparam,micnpool,zse)
          call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)
      !    print *, 'vmicsoil_c14'
          call vmicsoil_c14(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,ifsoc14,bgcopt,nyeqpool, &
                        zse,micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

       !   print *, 'calcost_c14'
          call calcost_c14(nx,isoc14,bgcopt,xopt,micparam,miccpool,micinput,zse,totcost1)

          miccpool%c12pooleqp(:) = miccpool%cpooleqp(:)
          miccpool%c12pooleqm(:) = miccpool%cpooleqm(:)

          isoc14 = 1
       !   print *, "isoc14 =",isoc14,'--getdata_c14'
          call getdata_c14(frac14c,f14c,filecluster,micinput,micparam,micnpool,zse)
          call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)
       !   print *, 'vmicsoil_c14'
          call vmicsoil_c14(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,ifsoc14,bgcopt,nyeqpool+2000, &
                        zse,micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

        !  print *, 'calcost_c14'
          call calcost_c14(nx,isoc14,bgcopt,xopt,micparam,miccpool,micinput,zse,totcost2)
          functn_c14 = totcost1+totcost2
        !  print *,"tot1 = ",totcost1
        !  print *,"tot2 = ",totcost2
           call screenout('c14run    ',jmodel,bgcopt,xopt,functn_c14)
      close(1)

!      functn = totcost

      call mic_deallocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_deallocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_deallocate_output(mp,micoutput)
      call mic_deallocate_cpool(mp,ms,miccpool)
      call mic_deallocate_npool(mp,ms,micnpool)
      deallocate(zse)
END function functn_c14


real(dp) function functn_frc1(nx,xparam16)
    !! Orchestrator for POC/MAOC fraction calibration mode.
    !!
    !! Fits simulated soil organic carbon fractions to observed
    !! fractionation data.  Reads configuration from `params1.txt`.

    integer, intent(in) :: nx
        !! Number of optimized parameters.
    real(dp), dimension(16), intent(in) :: xparam16
        !! Values of the `nx` optimized parameters.

     TYPE(mic_param_xscale)    :: micpxdef
     TYPE(mic_param_default)   :: micpdef
     TYPE(mic_parameter)       :: micparam
     TYPE(mic_input)           :: micinput
     TYPE(mic_global_input)    :: micglobal
     TYPE(mic_cpool)           :: miccpool
     TYPE(mic_npool)           :: micnpool
     TYPE(mic_output)          :: micoutput

    integer,   dimension(16)           :: nxopt
    real(dp),    dimension(16)           :: xopt
    real(dp)     :: totcost1
    integer    :: ifsoc14,kinetics,bgcopt,jopt,nyeqpool,isoc14,jglobal,jmodel
    integer :: jrestart,nparam
    character(len=140) :: frestart_in,frestart_out,foutput
    character(len=140) :: cfraction
    real(r_2), dimension(:), allocatable :: zse
    integer :: mpx

      jrestart=0;xopt(:)=1.0
      do nparam=1,16
         nxopt(nparam) = nparam
      enddo

      frestart_in='miccpool_in.nc'
      frestart_out='miccpool_out.nc'
      foutput='vmic_output.nc'

      open(1,file='params1.txt')
      read(1,*)
      read(1,*) jglobal,ifsoc14,kinetics,bgcopt,jopt,jrestart
      read(1,11) cfraction
11    format(a140)
      read(1,*) xopt(1:14)

      if(jopt==0) then
         read(1,*) nxopt(1:nx)
         do nparam=1,nx
            xopt(nxopt(nparam)) = xparam16(nparam)
         enddo
      endif
    !  print *, xopt

      close(1)
      !mp = 2206
      ntime=365

      totcost1 = 0.0
      nyeqpool= 1000
      isoc14 = 0
      jmodel=1;mpft=17;mbgc=12;nlon=1;nlat=1
      ms = 10
      allocate(zse(ms))
      zse(1) =0.02;zse(2)=0.04;zse(3)=0.06;zse(4)=0.08
      zse(5:8)=0.2;zse(9:10)=0.5
      call getdata_frc_dim(cfraction,mpx)
      mp = mpx
      call mic_allocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_allocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_allocate_output(mp,micoutput)
      call mic_allocate_cpool(mp,ms,miccpool)
      call mic_allocate_npool(mp,ms,micnpool)


    !  print *, "isoc14 =",isoc14,'--getdata_frc'
      call getdata_frc(cfraction,jglobal,bgcopt,micinput,micparam,micnpool,micglobal,zse)
      call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)

    !  print *, 'vmicsoil_frc1_cpu'
    !  call vmicsoil_frc1_cpu(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,ifsoc14,bgcopt,nyeqpool, &
    !                         zse,micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

      call vmicsoil_hwsd_cpu(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,bgcopt,nyeqpool, &
                         zse,micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

    !  print *, 'calcost_frc1'
      call calcost_frc1(nx,bgcopt,xopt,micpxdef,micparam,miccpool,micinput,micglobal,zse,totcost1)

      close(1)

      functn_frc1    = totcost1
      call screenout('fraction  ',jmodel,bgcopt,xopt,functn_frc1)
      call mic_deallocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_deallocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_deallocate_output(mp,micoutput)
      call mic_deallocate_cpool(mp,ms,miccpool)
      call mic_deallocate_npool(mp,ms,micnpool)
      deallocate(zse)
END function functn_frc1


  real(dp) function functn_soc_hwsd(nx,xparam16)
    !! Orchestrator for HWSD SOC profile calibration mode.
    !!
    !! Fits simulated soil carbon profiles to Harmonized World Soil
    !! Database observations.  Reads configuration from `params1.txt`.

    integer, intent(in) :: nx
        !! Number of optimized parameters.
    real(dp), dimension(16), intent(in) :: xparam16
        !! Values of the `nx` optimized parameters.

    integer,   dimension(16)           :: nxopt
    real(dp),    dimension(16)           :: xopt
    TYPE(mic_param_xscale)    :: micpxdef
    TYPE(mic_param_default)   :: micpdef
    TYPE(mic_parameter)       :: micparam
    TYPE(mic_input)           :: micinput
    TYPE(mic_global_input)    :: micglobal
    TYPE(mic_cpool)           :: miccpool
    TYPE(mic_npool)           :: micnpool
    TYPE(mic_output)          :: micoutput

    integer    :: ifsoc14,kinetics,bgcopt,jopt,nyeqpool,isoc14,jglobal,jmodel
    integer :: jrestart,nf,ok,nparam,mpx,timex
    character(len=140)  :: frestart_in,frestart_out,fparam_global,foutput
    character(len=140) :: fhwsdsoc,fmodis,fanoc
    real(r_2)     :: totcost1
    integer       :: ns
    real(r_2), dimension(:), allocatable :: zse


      isoc14=0;nyeqpool = 500;ok=0;totcost1=0.0

      jrestart=0;xopt(:)=1.0
      do nparam=1,16
         nxopt(nparam) = nparam
      enddo

      frestart_in='miccpool_in.nc'
      frestart_out='miccpool_out.nc'
      foutput='vmic_output.nc'

      open(1,file='params1.txt')
      read(1,*)
      read(1,*) jglobal,ifsoc14,kinetics,bgcopt,jopt,jrestart,jmodel
      read(1,101) fhwsdsoc
      read(1,101) fmodis
      read(1,101) fanoc
      read(1,*)   xopt(1:14)
      read(1,*) nxopt(1:nx)
      do nparam=1,nx
         xopt(nxopt(nparam)) = xparam16(nparam)
      enddo
      close(1)

!      print *, 'nx xparam16 =', nx, nxopt(1:nx),xparam16(1:nx)
!      print *, 'parameter values used= ', xopt
!      print *, 'ms zse', ms, zse(:)
101   format(a140)

      ! get dimensions
      call getdata_hwsd_dim(fhwsdsoc,mpx,timex)
      mp=mpx
      ntime=timex
      if(jmodel==1)                mpft=17   !CABLE
      if(jmodel==2 .or. jmodel==3) mpft=19   !ORCHIDEE
      mbgc=12;nlon=1;nlat=1
      ms=7
      allocate(zse(ms))
      zse(1:5)=0.2;zse(6:7)=0.5

      call mic_allocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_allocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_allocate_output(mp,micoutput)
      call mic_allocate_cpool(mp,ms,miccpool)
      call mic_allocate_npool(mp,ms,micnpool)

      call getdata_hwsd(fhwsdsoc,fmodis,fanoc,jglobal,bgcopt,jopt,jmodel,micparam,micglobal,zse)

      !  call profile()
      call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)

      call vmicsoil_hwsd_cpu(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,bgcopt,nyeqpool, &
                         zse,micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

      call calcost_hwsd3(nx,bgcopt,xopt,micpxdef,micparam,miccpool,micinput,micglobal,zse,totcost1)
      call mic_deallocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_deallocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_deallocate_output(mp,micoutput)
      call mic_deallocate_cpool(mp,ms,miccpool)
      call mic_deallocate_npool(mp,ms,micnpool)

      call screenout('hwsd_soc  ',jmodel,bgcopt,xopt,totcost1)


      functn_soc_hwsd = totcost1

      deallocate(zse)
END function functn_soc_hwsd

 real(dp) function functn_global4(nx,xparam16)
    !! Orchestrator for global-scale calibration with CABLE/ORCHIDEE forcing.
    !!
    !! Runs the model at global spatial resolution.  Reads configuration
    !! and grid patch definitions from `params1.txt`.  Not yet set up for
    !! SCE_UA optimization.

    integer, intent(in) :: nx
        !! Number of optimized parameters.
    real(dp), dimension(16), intent(in) :: xparam16
        !! Values of the `nx` optimized parameters.

    integer,   dimension(16)  :: nxopt
    real(dp),    dimension(16)  :: xopt
    TYPE(mic_param_xscale)    :: micpxdef
    TYPE(mic_param_default)   :: micpdef
    TYPE(mic_parameter)       :: micparam
    TYPE(mic_input)           :: micinput
    TYPE(mic_global_input)    :: micglobal
    TYPE(mic_cpool)           :: miccpool
    TYPE(mic_npool)           :: micnpool
    TYPE(mic_output)          :: micoutput

    integer    :: ifsoc14,kinetics,bgcopt,jopt,nyeqpool,isoc14,jglobal,jmodel
    integer    :: jrestart,nf,ok,nparam
    character(len=140) :: frestart_in,frestart_out,fparam_global,foutput
    character(len=140) :: fglobal(10)
    real(r_2) :: totcost1
    real(r_2), dimension(:), allocatable :: zse

      isoc14=0
      nyeqpool = 500
      ok=0
      totcost1=0.0
      jmodel=1;mpft=17;mbgc=10;ntime=365;nlon=192;nlat=112
      ms=7
      allocate(zse(ms))
      zse(1:5)=0.2;zse(6:7)=0.5

      frestart_in='miccpool_in.nc'
      frestart_out='miccpool_out.nc'
      foutput='vmic_output.nc'

      jrestart=0;xopt(:)=1.0
      do nparam=1,16
         nxopt(nparam) = nparam
      enddo
      xopt = 1.0

!      open(91,file='modobs.txt')
!      open(92,file='modobs2.txt')

      open(1,file='params1.txt')
      read(1,*)
      read(1,*) jglobal,ifsoc14,kinetics,bgcopt,jopt,jrestart,jmodel
      do nf=1,7
         read(1,101) fglobal(nf)
      enddo
      read(1,*)   xopt(1:14)
      read(1,*)   nxopt(1:nx)
      do nparam=1,nx
         xopt(nxopt(nparam)) = xparam16(nparam)
      enddo
      close(1)

      close(1)
101   format(a140)
      print *, xopt

      if(jmodel==2 .or. jmodel==3) then
         mpft=19; nlon=720; nlat=360
      endif

      ! reading global parameter values here      xopt =xparam16(1:nx)
      call getpatch_global(fglobal(1),jmodel,mp)
      print *, 'total number of patches= ', mp

      call mic_allocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_allocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_allocate_output(mp,micoutput)
      call mic_allocate_cpool(mp,ms,miccpool)
      call mic_allocate_npool(mp,ms,micnpool)

      print *, ' all  arrays are allocated!'

      if(jmodel==1)                call getdata_global4_cable(fglobal,jglobal,bgcopt,jopt,jmodel,micglobal,micparam,zse)
      if(jmodel==2 .or. jmodel==3) call getdata_global4_orchidee(fglobal,jglobal,bgcopt,jopt,jmodel,micglobal,micparam,zse)
      print *, 'global input data are read in'

      if(jopt==0) call getparam_global(fglobal(4),jmodel,micpxdef)     ! reading global parameter lookup table
      if(jopt==1) call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)  ! parameter optimization

      print *, 'vmicsoil_global'
      call vmicsoil_hwsd_cpu(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,bgcopt,nyeqpool, &
                         zse,micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

      call calcost_hwsd3(nx,bgcopt,xopt,micpxdef,micparam,miccpool,micinput,micglobal,zse,totcost1)

      call mic_deallocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_deallocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_deallocate_output(mp,micoutput)
      call mic_deallocate_cpool(mp,ms,miccpool)
      call mic_deallocate_npool(mp,ms,micnpool)

!      close(91)
!      close(92)

      functn_global4=totcost1
      print *, 'total cost =', totcost1
      deallocate(zse)

END function functn_global4

  real(dp) function functn_soc_aust(nx,xparam16)
    !! Orchestrator for Australian SOC profile calibration.
    !!
    !! Fits simulated soil carbon profiles to Australian soil data.
    !! Reads configuration from `params1.txt`.

    integer, intent(in) :: nx
        !! Number of optimized parameters.
    real(dp), dimension(16), intent(in) :: xparam16
        !! Values of the `nx` optimized parameters.

    integer,   dimension(16)           :: nxopt
    real(dp),    dimension(16)           :: xopt
    TYPE(mic_param_xscale)    :: micpxdef
    TYPE(mic_param_default)   :: micpdef
    TYPE(mic_parameter)       :: micparam
    TYPE(mic_input)           :: micinput
    TYPE(mic_global_input)    :: micglobal
    TYPE(mic_cpool)           :: miccpool
    TYPE(mic_npool)           :: micnpool
    TYPE(mic_output)          :: micoutput

    integer    :: ifsoc14,kinetics,bgcopt,jopt,nyeqpool,isoc14,jglobal,jmodel
    integer    :: jrestart,nf,ok,nparam,mpx,timex
    character(len=140) :: frestart_in,frestart_out,foutput
    character(len=140) :: faustsoc
    real(r_2)     :: totcost1
    integer       :: ns
    real(r_2), dimension(:), allocatable :: zse


      isoc14=0;nyeqpool = 500;ok=0;totcost1=0.0

      jrestart=0;xopt(:)=1.0
      do nparam=1,16
         nxopt(nparam) = nparam
      enddo

      frestart_in='miccpool_in.nc'
      frestart_out='miccpool_out.nc'
      foutput='vmic_output.nc'

      open(1,file='params1.txt')
      read(1,*)
      read(1,*) jglobal,ifsoc14,kinetics,bgcopt,jopt,jrestart,jmodel
      read(1,101) faustsoc
      read(1,*)   xopt(1:14)
      read(1,*) nxopt(1:nx)
      do nparam=1,nx
         xopt(nxopt(nparam)) = xparam16(nparam)
      enddo
      close(1)

101   format(a140)

      ! get dimensions
      call getdata_aust_dim(faustsoc,mpx,timex)
      mp=mpx
      ntime=timex
      mpft=9;mbgc=9;nlon=1;nlat=1
      ms=5
      allocate(zse(ms))
      zse(1:5)=0.1


      call mic_allocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_allocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_allocate_output(mp,micoutput)
      call mic_allocate_cpool(mp,ms,miccpool)
      call mic_allocate_npool(mp,ms,micnpool)

      call getdata_aust(faustsoc,jglobal,bgcopt,jopt,jmodel,micparam,micglobal,zse)

      !  call profile()
      call vmic_param_xscale(xopt,bgcopt,jmodel,micpxdef)

      call vmicsoil_hwsd_cpu(jrestart,frestart_in,frestart_out,foutput,kinetics,isoc14,bgcopt,nyeqpool, &
                         zse,micpxdef,micpdef,micparam,micinput,micglobal,miccpool,micnpool,micoutput)

      call calcost_aust(nx,bgcopt,xopt,micpxdef,micparam,miccpool,micinput,micglobal,zse,totcost1)
      call mic_deallocate_parameter(mpft,mbgc,mp,ms,micpxdef,micparam)
      call mic_deallocate_input(mp,ms,nlon,nlat,ntime,micinput,micglobal)
      call mic_deallocate_output(mp,micoutput)
      call mic_deallocate_cpool(mp,ms,miccpool)
      call mic_deallocate_npool(mp,ms,micnpool)

      functn_soc_aust = totcost1

      deallocate(zse)
END function functn_soc_aust

end module function_module
