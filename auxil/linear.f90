program lreg
  ! linear regression of the modelled(y) an observed(x) soil carbon
  implicit none
  real, dimension(10000,4) :: varx,vary
  real, dimension(10000)   :: x,y
  real :: a,b,r,a1,a2,b1,b2,r1,r2
  integer :: ndata,n,n1,n2,n3,n4,n5,i,ndataset,idata,ms,is
  real :: x1,x2,x3,x4,y1,y2,y3,y4,obs1,mod1,area
  ! idata=3
  integer :: siteid,pft,isoil,sorder,ns,bgctype
  real    :: xobs,xmod,fracpocm,fracmaocm,fracmicm,fraclabm
  real, dimension(7)  :: zsoil,profobs,profmod
  character(len=20) :: dataset(4)

    open(1,file="input.txt")
    open(2,file="output.txt")
    open(3,file="linefit.txt")

    open(10,file="param_lreg.txt")
    read(10,*) idata
    close(10)
    if(idata==2) ndataset=4
    if(idata==3) then
       ndataset=1
       ms =7
       zsoil(1) = 0.1
       zsoil(2) = 0.3
       zsoil(3) = 0.5
       zsoil(4) = 0.7
       zsoil(5) = 0.9
       zsoil(6) = 1.25
       zsoil(7) = 1.75
       open(20,file="profobs.txt")
       open(21,file="profmod.txt")
    end if
    ndata=0
    ! for idata==2
    dataset(1) ="total POC"
    dataset(2) ="total MAOC"
    dataset(3) ="POC fraction"
    dataset(4) ="MAOC fraction"
    ! for idata==3
    dataset(1) ="total SOC"

    do n=1,10000
       if(idata==1) then   ! 14C data
        print *, "yet to be done"
        stop
       end if !idata=1

       if(idata==2) then  ! soil C fraciton data
          !  write(91,901) micparam%dataid(np),micparam%siteid(np),micparam%bgctype(np),micparam%top(np),micparam%bot(np),
          !   xobsp(np),xmodp(np),xobsm(np),xmodm(np),xobsfracp(np),xmodfracp(np),xobsfracm(np),xmodfracm(np)
          read(1,*,end=15) n1,n2,n3,n4,n5,y1,x1,y2,x2,y3,x3,y4,x4
          if(min(y1,y2,y3,y4)>0.0) then
             ndata=ndata +1
             varx(ndata,1) = x1
             vary(ndata,1) = y1
             varx(ndata,2) = x2
             vary(ndata,2) = y2
             varx(ndata,3) = x3
             vary(ndata,3) = y3
             varx(ndata,4) = x4
             vary(ndata,4) = y4
             write(2,11) n1,n2,n3,n4,n5,y1,x1,y2,x2,y3,x3,y4,x4
          end if
       end if !idata=2

       if(idata==3) then   !HWSD data
          profmod(:) = 1.0e-6; profobs(:) = 1.0e-6
          do ns=1,ms
             read(1,*,end=15) siteid,pft,isoil,sorder,bgctype,area,is,&
                              xobs,xmod,fracpocm,fracmaocm,fracmicm,fraclabm
             if(min(xobs,xmod)>0.0) then
                ndata=ndata +1
                varx(ndata,1) = log(xobs)
                vary(ndata,1) = log(xmod)
                write(2,31) siteid,pft,isoil,sorder,bgctype,area,ns,xobs,xmod
             end if
             if(ns==1) then
                obs1=xobs;  mod1=xmod
                profobs(1) = 1.0; profmod(1) = 1.0
             else
                profobs(ns) = xobs/obs1
                profmod(ns) = xmod/mod1
             end if
           end do
           do ns=1,ms
              profobs(ns) = log(profobs(ns))
              profmod(ns) = log(profmod(ns))
              print *, ns, profobs(ns),profmod(ns)
           end do
           ! computing the intercept and slope of the normalized exponetial profiles
           call linreg(ms,zsoil,profobs,a1,b1,r1)
           call linreg(ms,zsoil,profmod,a2,b2,r2)
           ! write the profile results
           write(20,201) siteid,pft,isoil,sorder,a1,b1,r1,exp(profobs(1:ms))
           write(21,211) siteid,pft,isoil,sorder,a2,b2,r2,exp(profmod(1:ms))

        end if  !idata=3

    end do   !enddo
15  print *, "total number datapoints = ", ndata

    do n=1,ndataset
       x(:)=varx(:,n)
       y(:)=vary(:,n)
       call linreg(ndata,x,y,a,b,r)
       write(3,21) dataset(n),n,ndata,a,b,r*r
    end do
    close(1)
    close(2)
    close(3)
    if(idata==3) then
       close(20)
       close(21)
    end if
11  format(5(i6,1x),10(f10.5,1x))
21  format(a20,1x,2(i6,1x),3(f10.4,1x))
31  format(5(i6,1x),f7.3,1x,i3,1x,8(f10.5,1x))
201 format("obs: ",4(i4,1x),20(f9.4,1x))
211 format("mod: ",4(i4,1x),20(f9.4,1x))
    end program lreg


subroutine linreg(ndata,x,y,a,b,r)

   implicit none                                                                    ! no default data types
   ! input
   integer :: ndata
   real, dimension(ndata)       ::  x, y
   ! output
   real :: a,b,r
   !local variables
   integer :: n
   real :: sumx,sumy,sumx2, sumy2, sumxy

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

end subroutine linreg
