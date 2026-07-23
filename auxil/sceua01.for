      program scemain
c$debug
c
c  THIS IS THE MAIN PROGRAM CALLING SUBROUTINES SCEIN AND SCEUA
      implicit real*8 (a-h,o-z)
      dimension a(16), bl(16), bu(16), jseed(10)
      data jseed/2,3,5,7,11,13,17,19,23,29/

      write (*,*) ' ENTER THE MAIN PROGRAM --- '
c
      call scein(a,bl,bu,nopt,maxn,kstop,pcento,iseed,
     &           ngs,npg,nps,nspl,mings,iniflg,iprint)
c
      if (iseed .gt. 0) then
        nrun = min(iseed,10)
      else
        nrun = 1
      end if

      do i=1,nrun
        if (nrun .ne. 1) iseed = jseed(i)
        write (*,*) '@ SCE-UA Run Number',i,' Random Seed Value',iseed
        call sceua(a,bl,bu,nopt,maxn,kstop,pcento,iseed,
     &           ngs,npg,nps,nspl,mings,iniflg,iprint)
      end do
c
c  END OF THE MAIN PROGRAM
      stop
      end


      subroutine scein(a,bl,bu,nopt,maxn,kstop,pcento,iseed,
     &                 ngs,npg,nps,nspl,mings,iniflg,iprint)
c$debug
c
c   THIS SUBROUTINE READS AND PRINTS THE INPUT VARIABLES FOR
c   SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
c     -- Version 2.1
c
c   WRITTEN BY QINGYUN DUAN - UNIVERSITY OF ARIZONA, APRIL 1992
c
c
      implicit real*8 (a-h,o-z)
      dimension a(16),bl(16),bu(16)
      common /iopar/ in,ipr
      character*10 pcntrl,deflt,usrsp
      character*4 reduc,initl,ysflg,noflg,xname(16)
c
      data deflt/' DEFAULT  '/
      data usrsp/'USER SPEC.'/
      data ysflg/'YES '/
      data noflg/'NO  '/
      data xname /'  X1','  X2','  X3','  X4','  X5','  X6','  X7',
     &'  X8','  X9',' X10',' X11',' X12',' X13',' X14',' X15',' X16'/
c
      write (*,*) ' ENTER THE SCEIN SUBROUTINE --- '
c
c
c  INITIALIZE I/O VARIABLES
      in = 2
      ipr = 5
      open(unit=in,file='scein.dat',status='old')
      open(unit=ipr,file='sceout.dat',status='unknown')
c
      ierror = 0
      iwarn = 0
      write(ipr,700)
  700 format(10x,'SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION',
     &       /,10x,46(1h=))
c
c
c  READ THE SCE CONTROL PARAMETERS
      ideflt = 0
Cypw      read(in,800) maxn,kstop,pcento,ngs,iseed,ideflt
C read using free format: ypw
      read(in,*) maxn,kstop,pcento,ngs,iseed,ideflt

  800 format(2i5,f5.2,3i5)
      if (iseed .eq. 0) iseed = 1969
c
c
c  IF ideflt IS EQUAL TO 1, READ THE SCE CONTROL PARAMETERS
      if (ideflt .eq. 1) Then
C        read(in,810) npg,nps,nspl,mings,iniflg,iprint
C read using free format 
        read(in,*) npg,nps,nspl,mings,iniflg,iprint
  810   format(6i5)
        pcntrl = usrsp
      else
        read(in,*)
        pcntrl = deflt
      end if
c
c
c  READ THE INITIAL PARAMETER VALUES AND THE PARAMETER BOUNDS
      iopt = 0
  820 iopt = iopt + 1
C read using free format
C      read(in,830,end=840) a(iopt), bl(iopt), bu(iopt)
      read(in,*,end=840) a(iopt), bl(iopt), bu(iopt)
  830 format(3f10.3)
      go to 820
  840 nopt = iopt - 1
c
c
c  IF ideflt IS EQUAL TO 0, SET THE SCE CONTROL PARAMETERS TO
c  THE DEFAULT VALUES
      if (ideflt .eq. 0) then
        npg = 2*nopt + 1
        nps = nopt + 1
        nspl = npg
        mings = ngs
        iniflg = 0
        iprint = 0
      end if
c
c
c  CHECK IF THE SCE CONTROL PARAMETERS ARE VALID
      if (ngs .lt. 1 .or. ngs .ge. 1320) then
        write(ipr,900) ngs
  900   format(//,1x,'**ERROR** NUMBER OF COMPLEXES IN INITIAL ',
     *         ' POPULATION ',i5,' IS NOT A VALID CHOICE')
        ierror = ierror + 1
      end if
c
      if (kstop .lt. 0 .or. kstop .ge. 20) then
        write(ipr,901) kstop
  901   format(//,1x,'**WARNING** THE NUMBER OF SHUFFLING LOOPS IN',
     *  ' WHICH THE CRITERION VALUE MUST CHANGE ',/,13x,'SHOULD BE',
     *  ' GREATER THAN 0 AND LESS THAN 10.  ','kstop = ',i2,
     *  ' WAS SPECIFIED.'/,13x,'BUT kstop = 5 WILL BE USED INSTEAD.')
        iwarn = iwarn + 1
        kstop=5
      end if
c
      if (mings .lt. 1 .or. mings .gt. ngs) then
        write(ipr,902) mings
  902   format(//,1x,'**WARNING** THE MINIMUM NUMBER OF COMPLEXES ',
     *         i2,' IS NOT A VALID CHOICE. SET IT TO DEFAULT')
        iwarn = iwarn + 1
        mings = ngs
      end if
c
      if (npg .lt. 2 .or. npg .gt. 1320/max(ngs,1)) then
        write(ipr,903) npg
  903   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN A COMPLEX ',
     *         I4,' IS NOT A VALID CHOICE, SET IT TO DEFAULT')
        iwarn = iwarn + 1
        npg = 2*nopt+1
      end if
c
      if (nps.lt.2 .or. nps.gt.npg .or. nps.gt.50) then
        write(ipr,904) nps
  904   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN A SUB-',
     *  'COMPLEX ',i4,' IS NOT A VALID CHOICE, SET IT TO DEFAULT')
        iwarn = iwarn + 1
        nps = nopt + 1
      end if
c
      if (nspl .lt. 1) then
        write(ipr,905) nspl
  905   format(//,1x,'**WARNING** THE NUMBER OF EVOLUTION STEPS ',
     *         'TAKEN IN EACH COMPLEX BEFORE SHUFFLING ',I4,/,13x,
     *         'IS NOT A VALID CHOICE, SET IT TO DEFAULT')
        iwarn = iwarn + 1
        nspl = npg
      end if
c
c  COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPULATION
      npt = ngs * npg
c
      if (npt .gt. 1320) then
        write(ipr,906) npt
  906   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN INITIAL ',
     *         'POPULATION ',i5,' EXCEED THE POPULATION LIMIT,',/,13x,
     *         'SET NGS TO 2, AND NPG, NPS AND NSPL TO DEFAULTS')
        iwarn = iwarn + 1
        ngs = 2
        npg = 2*nopt + 1
        nps = nopt + 1
        nspl = npg
      end if
c
c  PRINT OUT THE TOTAL NUMBER OF ERROR AND WARNING MESSAGES
      if (ierror .ge. 1) write(ipr,907) ierror
  907 format(//,1x,'*** TOTAL NUMBER OF ERROR MESSAGES IS ',i2)
c
      if (iwarn .ge. 1) write(ipr,908) iwarn
  908 format(//,1x,'*** TOTAL NUMBER OF WARNING MESSAGES IS ',i2)
c
      if (mings .lt. ngs) then
        reduc = ysflg
      else
        reduc = noflg
      end if
c
      if (iniflg .ne. 0) then
        initl = ysflg
      else
        initl = noflg
      end if
c
c
c  PRINT SHUFFLED COMPLEX EVOLUTION OPTIMIZATION OPTIONS
  104 write(ipr,910)
  910 format(//,2x,'SCE CONTROL',5x,'MAX TRIALS',5x,
     &'REQUIRED IMPROVEMENT',5x,'RANDOM',/,3x,'PARAMETER',8x,
     &'ALLOWED',6x,'PERCENT',4x,'NO. LOOPS',6x,'SEED',/,
     &2x,11(1h-),5x,10(1H-),5x,7(1h-),4x,9(1h-),5x,6(1h-))
c
      pcenta=pcento*100.
      write(ipr,912) pcntrl,maxn,pcenta,kstop,iseed
  912 format(3x,a10,7x,i5,10x,f3.1,9x,i2,9x,i5)
      write(ipr,914) ngs,npg,npt,nps,nspl
  914 format(//,18x,'SCE ALGORITHM CONTROL PARAMETERS',/,18x,32(1H=),
     &//,2x,'NUMBER OF',5x,'POINTS PER',5x,'POINTS IN',6x,'POINTS PER',
     &4x,'EVOL. STEPS',/,2x,'COMPLEXES',6X,'COMPLEX',6x,'INI. POPUL.',
     &5x,'SUB-COMPLX',4x,'PER COMPLEX',/,2x,9(1h-),5x,10(1h-),4x,
     &11(1h-),5x,10(1h-),4x,11(1h-),5x,/,2x,5(i5,10x))
      write(ipr,915) reduc,mings,initl
  915 format(//,15x,'COMPLX NO.',5x,'MIN COMPLEX',5x,'INI. POINT',/,
     &15x,'REDUCTION',6x,'NO. ALLOWED',6x,'INCLUDED',/,
     &15x,10(1h-),5x,11(1h-),5x,10(1h-),/,18x,a4,6x,i8,13x,a4)
      write(ipr,916)
  916 format(//,8x,'INITIAL PARAMETER VALUES AND PARAMETER BOUNDS',/,
     &       8x,45(1h=),//,2x,'PARAMETER',5x,'INITIAL VALUE',5x,
     &       'LOWER BOUND',5x,'UPPER BOUND',/,2x,9(1h-),5x,13(1h-),5x,
     &       11(1h-),5x,11(1h-))
      do 920 i = 1, nopt
        write(ipr,918) xname(i),a(i),bl(i),bu(i)
  918   format(4x,a4,4x,3(6x,f10.3))
  920 continue
      if (ierror .ge. 1) then
      write(ipr,922)
  922 format(//,'*** THE OPTIMIZATION SEARCH IS NOT CONDUCTED BECAUSE',
     &       ' OF INPUT DATA ERROR ***')
      stop
      end if
c
C  END OF SUBROUTINE SCEIN
      return
      end
      subroutine sceua(a,bl,bu,nopt,maxn,kstop,pcento,iseed,
     &                 ngs,npg,nps,nspl,mings,iniflg,iprint)
c$debug
c
c
c  SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
c     -- Version 2.1
c
c  by QINGYUN DUAN
c  DEPARTMENT OF HYDROLOGY & WATER RESOURCES
c  UNIVERSITY OF ARIZONA, TUCSON, AZ 85721
c  (602) 621-9360, email: duan@hwr.arizona.edu
c
c  WRITTEN IN OCTOBER 1990.
c  REVISED IN AUGUST 1991
c  REVISED IN APRIL 1992
c
c  STATEMENT BY AUTHOR:
c  --------------------
c
c     This general purpose global optimization program is developed at
c     the Department of Hydrology & Water Resources of the University
c     of Arizona.  Further information regarding the SCE-UA method can
c     be obtained from Dr. Q. Duan, Dr. S. Sorooshian or Dr. V.K. Gupta
c     at the address and phone number listed above.  We request all
c     users of this program make proper reference to the paper entitled
c     'Effective and Efficient Global Optimization for Conceptual
c     Rainfall-runoff Models' by Duan, Q., S. Sorooshian, and V.K. Gupta,
c     Water Resources Research, Vol 28(4), pp.1015-1031, 1992.
c
c
c  LIST OF INPUT ARGUEMENT VARIABLES
c
c     a(.) = initial parameter set
c     bl(.) = lower bound on parameters
c     bu(.) = upper bound on parameters
c     nopt = number of parameters to be optimized
c
c
c  LIST OF SCE ALGORITHMIC CONTROL PARAMETERS:
c
c     ngs = number of complexes in the initial population
c     npg = number of points in each complex
c     npt = total number of points in initial population (npt=ngs*npg)
c     nps = number of points in a sub-complex
c     nspl = number of evolution steps allowed for each complex before
c         complex shuffling
c     mings = minimum number of complexes required, if the number of
c         complexes is allowed to reduce as the optimization proceeds
c     iseed = initial random seed
c     iniflg = flag on whether to include the initial point in population
c         = 0, not included
c         = 1, included
c     iprint = flag for controlling print-out after each shuffling loop
c         = 0, print information on the best point of the population
c         = 1, print information on every point of the population
c
c
c  CONVERGENCE CHECK PARAMETERS
c
c     maxn = max no. of trials allowed before optimization is terminated
c     kstop = number of shuffling loops in which the criterion value must
c         chang by the given percentage before optimization is terminated
c     pcento = percentage by which the criterion value must change in
c         given number of shuffling loops
c     ipcnvg = flag indicating whether parameter convergence is reached
c         (i.e., check if gnrng is less than 0.001)
c         = 0, parameter convergence not satisfied
c         = 1, parameter convergence satisfied
c
c
c  LIST OF LOCAL VARIABLES
c     x(.,.) = coordinates of points in the population
c     xf(.) = function values of x(.,.)
c     xx(.) = coordinates of a single point in x
c     cx(.,.) = coordinates of points in a complex
c     cf(.) = function values of cx(.,.)
c     s(.,.) = coordinates of points in the current simplex
c     sf(.) = function values of s(.,.)
c     bestx(.) = best point at current shuffling loop
c     bestf = function value of bestx(.)
c     worstx(.) = worst point at current shuffling loop
c     worstf = function value of worstx(.)
c     xnstd(.) = standard deviation of parameters in the population
c     gnrng = normalized geometric mean of parameter ranges
c     lcs(.) = indices locating position of s(.,.) in x(.,.)
c     bound(.) = bound on ith variable being optimized
c     ngs1 = number of complexes in current population
c     ngs2 = number of complexes in last population
c     iseed1 = current random seed
c     criter(.) = vector containing the best criterion values of the last
c         10 shuffling loops
c
      implicit real*8 (a-h,o-z)
c
c  ARRAYS FROM THE INPUT DATA
      dimension a(16),bl(16),bu(16)
c
c  LOCAL ARRAYS
      dimension x(2000,16),xx(16),bestx(16),worstx(16),xf(2000)
      dimension s(50,16),sf(50),lcs(50),cx(2000,16),cf(2000)
      dimension xnstd(16),bound(16),criter(20),unit(16)
c
      common /iopar/ in,ipr
c
      character*4 xname(16)
      data xname /'  X1','  X2','  X3','  X4','  X5','  X6','  X7',
     &'  X8','  X9',' X10',' X11',' X12',' X13',' X14',' X15',' X16'/
c
      write (*,*) ' ENTER THE SCEUA SUBROUTINE --- '     
c
c  INITIALIZE VARIABLES
      nloop = 0
      loop = 0
      igs = 0
      nopt1 = 8
      if (nopt.lt.8) nopt1 = nopt
      nopt2 = 12
      if (nopt.lt.12) nopt2 = nopt
c
c  INITIALIZE RANDOM SEED TO A NEGATIVE INTEGER
      iseed1 = -abs(iseed)
c
c  COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPUALTION
      npt = ngs * npg
      ngs1 = ngs
      npt1 = npt
c
      write(ipr,400)
      write (*,*) ' ***  Evolution Loop Number ',nloop
c
c  COMPUTE THE BOUND FOR PARAMETERS BEING OPTIMIZED
      do j = 1, nopt
        bound(j) = bu(j) - bl(j)
        unit(j) = 1.0
      end do
c
c  COMPUTE THE FUNCTION VALUE OF THE INITIAL POINT
      fa = functn(nopt,a)
c
c  PRINT THE INITIAL POINT AND ITS CRITERION VALUE
      write(ipr,500)
      write(ipr,510) (xname(j),j=1,nopt2)
      write(ipr,520) fa,(a(j),j=1,nopt2)
      if (nopt.lt.13) go to 101
      write(ipr,530) (xname(j),j=13,nopt)
      write(ipr,540) (a(j),j=13,nopt)
  101 continue
c
c  GENERATE AN INITIAL SET OF npt1 POINTS IN THE PARAMETER SPACE
c  IF iniflg IS EQUAL TO 1, SET x(1,.) TO INITIAL POINT a(.)
      if (iniflg .eq. 1) then
        do j = 1, nopt
          x(1,j) = a(j)
        end do
        xf(1) = fa
c
c  ELSE, GENERATE A POINT RANDOMLY AND SET IT EQUAL TO x(1,.)
      else
        call getpnt(nopt,1,iseed1,xx,bl,bu,unit,bl)
        do j=1, nopt
          x(1,j) = xx(j)
        end do
        xf(1) = functn(nopt,xx)
      end if
      icall = 1
      if (icall .ge. maxn) go to 9000
c
c  GENERATE npt1-1 RANDOM POINTS DISTRIBUTED UNIFORMLY IN THE PARAMETER
c  SPACE, AND COMPUTE THE CORRESPONDING FUNCTION VALUES
      do i = 2, npt1
        call getpnt(nopt,1,iseed1,xx,bl,bu,unit,bl)

        do j = 1, nopt
          x(i,j) = xx(j)
        end do
        xf(i) = functn(nopt,xx)
        icall = icall + 1
        if (icall .ge. maxn) then
          npt1 = i
          go to 45
        end if
      end do
c
c  ARRANGE THE POINTS IN ORDER OF INCREASING FUNCTION VALUE
   45 call sort(npt1,nopt,x,xf)
c
c  RECORD THE BEST AND WORST POINTS
      do j = 1, nopt
        bestx(j) = x(1,j)
        worstx(j) = x(npt1,j)
      end do
      bestf = xf(1)
      worstf = xf(npt1)
c
c  COMPUTE THE PARAMETER RANGE FOR THE INITIAL POPULATION
      call parstt(npt1,nopt,x,xnstd,bound,gnrng,ipcnvg)
c
c  PRINT THE RESULTS FOR THE INITIAL POPULATION
      write(ipr,600)
      write(ipr,610) (xname(j),j=1,nopt1)
      if (nopt .lt. 9) go to 201
      write(ipr,620) (xname(j),j=9,nopt)
  201 continue
      write(ipr,630) nloop,icall,ngs1,bestf,worstf,gnrng,
     &               (bestx(j),j=1,nopt1)
      if (nopt .lt. 9) go to 301
      write(ipr,640) (bestx(j),j=9,nopt)
  301 continue
      if (iprint .eq. 1) then
        write(ipr,650) nloop
        do i = 1, npt1
          write(ipr,660) xf(i),(x(i,j),j=1,nopt1)
          if (nopt .lt. 9) go to 401
          write(ipr,640) (x(i,j),j=9,nopt)
  401   end do
      end if
c
      if (icall .ge. maxn) go to 9000
      if (ipcnvg .eq. 1) go to 9200
c
c  BEGIN THE MAIN LOOP ----------------
 1000 continue
      nloop = nloop + 1
c
      write (*,*) ' ***  Evolution Loop Number ',nloop
c
c  BEGIN LOOP ON COMPLEXES
      do 3000 igs = 1, ngs1
c
c  ASSIGN POINTS INTO COMPLEXES
      do k1 = 1, npg
        k2 = (k1-1) * ngs1 + igs
        do j = 1, nopt
          cx(k1,j) = x(k2,j)
        end do
        cf(k1) = xf(k2)
      end do
c
c  BEGIN INNER LOOP - RANDOM SELECTION OF SUB-COMPLEXES ---------------
      do 2000 loop = 1, nspl
c
c  CHOOSE A SUB-COMPLEX (nps points) ACCORDING TO A LINEAR
c  PROBABILITY DISTRIBUTION
      if (nps .eq. npg) then
        do k = 1, nps
          lcs(k) = k
        end do
        go to 85
      end if
c
      rand = ran1(iseed1)
      lcs(1) = 1 + dint(npg + 0.5 - dsqrt( (npg+.5)**2 -
     &         npg * (npg+1) * rand ))
      do k = 2, nps
   60   rand = ran1(iseed1)
        lpos = 1 + dint(npg + 0.5 - dsqrt((npg+.5)**2 -
     &         npg * (npg+1) * rand ))
        do k1 = 1, k-1
          if (lpos .eq. lcs(k1)) go to 60
        end do
        lcs(k) = lpos
      end do
c
c  ARRANGE THE SUB-COMPLEX IN ORDER OF INCEASING FUNCTION VALUE
      call sort1(nps,lcs)
c
c  CREATE THE SUB-COMPLEX ARRAYS
   85 do k = 1, nps
        do j = 1, nopt
          s(k,j) = cx(lcs(k),j)
        end do
        sf(k) = cf(lcs(k))
      end do
c
c  USE THE SUB-COMPLEX TO GENERATE NEW POINT(S)
      call cce(nopt,nps,s,sf,bl,bu,xnstd,icall,maxn,iseed1)
c
c  IF THE SUB-COMPLEX IS ACCEPTED, REPLACE THE NEW SUB-COMPLEX
c  INTO THE COMPLEX
      do k = 1, nps
        do j = 1, nopt
          cx(lcs(k),j) = s(k,j)
        end do
        cf(lcs(k)) = sf(k)
      end do
c
c  SORT THE POINTS
      call sort(npg,nopt,cx,cf)
c
c  IF MAXIMUM NUMBER OF RUNS EXCEEDED, BREAK OUT OF THE LOOP
      if (icall .ge. maxn) go to 2222
c
c  END OF INNER LOOP ------------
 2000 continue
 2222 continue
c
c  REPLACE THE NEW COMPLEX INTO ORIGINAL ARRAY x(.,.)
      do k1 = 1, npg
        k2 = (k1-1) * ngs1 + igs
        do j = 1, nopt
          x(k2,j) = cx(k1,j)
        end do
        xf(k2) = cf(k1)
      end do
      if (icall .ge. maxn) go to 3333
c
c  END LOOP ON COMPLEXES
 3000 continue
c
c  RE-SORT THE POINTS
 3333 call sort(npt1,nopt,x,xf)
c
c  RECORD THE BEST AND WORST POINTS
      do j = 1, nopt
        bestx(j) = x(1,j)
        worstx(j) = x(npt1,j)
      end do
      bestf = xf(1)
      worstf = xf(npt1)
c
c  TEST THE POPULATION FOR PARAMETER CONVERGENCE
      call parstt(npt1,nopt,x,xnstd,bound,gnrng,ipcnvg)
c
c  PRINT THE RESULTS FOR CURRENT POPULATION
      if (mod(nloop,5) .ne. 0) go to 501
      write(ipr,610) (xname(j),j=1,nopt1)
      if (nopt .lt. 9) go to 501
      write(ipr,620) (xname(j),j=9,nopt)
  501 continue
      write(ipr,630) nloop,icall,ngs1,bestf,worstf,gnrng,
     &               (bestx(j),j=1,nopt1)
      if (nopt.lt.9) go to 601
      write(ipr,640) (bestx(j),j=9,nopt)
  601 continue
      if (iprint .eq. 1) then
        write(ipr,650) nloop
        do i = 1, npt1
          write(ipr,660) xf(i),(x(i,j),j=1,nopt1)
          if (nopt .lt. 9) go to 701
          write(ipr,640) (x(i,j),j=9,nopt)
  701   end do
      end if
c
c  TEST IF MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED
      if (icall .ge. maxn) go to 9000
c
c  COMPUTE THE COUNT ON SUCCESSIVE LOOPS W/O FUNCTION IMPROVEMENT
      criter(20) = bestf
      if (nloop .lt. (kstop+1)) go to 132
      denomi = dabs(criter(20-kstop) + criter(20)) / 2.
      timeou = dabs(criter(20-kstop) - criter(20)) / denomi
      if (timeou .lt. pcento) go to 9100
  132 continue
      do l = 1, 19
        criter(l) = criter(l+1)
      end do
c
c  IF POPULATION IS CONVERGED INTO A SUFFICIENTLY SMALL SPACE
      if (ipcnvg .eq. 1) go to 9200
c
c  NONE OF THE STOPPING CRITERIA IS SATISFIED, CONTINUE SEARCH
c
c  CHECK FOR COMPLEX NUMBER REDUCTION
      if (ngs1 .gt .mings) then
        ngs2 = ngs1
        ngs1 = ngs1 - 1
        npt1 = ngs1 * npg
        call comp(nopt,npt1,ngs1,ngs2,npg,x,xf,cx,cf)
      end if
c
c  END OF MAIN LOOP -----------
      go to 1000
c
c  SEARCH TERMINATED
 9000 continue
      write(ipr,800) maxn,loop,igs,nloop
      go to 9999
 9100 continue
      write(ipr,810) pcento*100.,kstop
      go to 9999
 9200 write(ipr,820) gnrng*100.
 9999 continue
c
c  PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
      write(ipr,830)
      write(ipr,510) (xname(j),j=1,nopt2)
      write(ipr,520) bestf,(bestx(j),j=1,nopt2)
      if (nopt.lt.13) go to 801
      write(ipr,530) (xname(j),j=13,nopt)
      write(ipr,540) (bestx(j),j=13,nopt)
  801 continue
c
c  END OF SUBROUTINE SCEUA
      return
  400 format(//,2x,50(1h=),/,2x,'ENTER THE SHUFFLED COMPLEX EVOLUTION',
     &       ' GLOBAL SEARCH',/,2x,50(1h=))
  500 format(//,'*** PRINT THE INITIAL POINT AND ITS CRITERION ',
     &       'VALUE ***')
  510 format(/,' CRITERION',12(6x,a4),/1x,60(1h-))
  520 format(g10.3,12f10.3)
  530 format(10x,12(6x,a4))
  540 format(10x,12f10.3)
  600 format(//,1x,'*** PRINT THE RESULTS OF THE SCE SEARCH ***')
  610 format(/,1x,'LOOP',1x,'TRIALS',1x,'COMPLXS',2x,'BEST F',3x,
     &       'WORST F',3x,'PAR RNG',1x,8(6x,a4))
  620 format(49x,8(6x,a4))
  630 format(i5,1x,i5,3x,i5,3g10.3,8(f10.3))
  640 format(49x,8(f10.3))
  650 format(/,1x,'POPULATION AT LOOP ',i3,/,1x,22(1h-))
  660 format(15x,g10.3,20x,8(f10.3))
  800 format(//,1x,'*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE',
     &       ' LIMIT ON THE MAXIMUM',/,5x,'NUMBER OF TRIALS ',i5,
     &       ' EXCEEDED.  SEARCH WAS STOPPED AT',/,5x,'SUB-COMPLEX ',
     &       i3,' OF COMPLEX ',i3,' IN SHUFFLING LOOP ',i3,' ***')
  810 format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE CRITERION',
     &       ' VALUE HAS NOT CHANGED ',/,5x,f5.2,' PERCENT IN',i3,
     &       ' SHUFFLING LOOPS ***')
  820 format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE POPULATION',
     &       ' HAS CONVERGED INTO ',/,4x,f5.2,' PERCENT OF THE',
     &       ' FEASIBLE SPACE ***')
  830 format(//,'*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS',
     &       ' CRITERION VALUE ***')
      end
c
c
c
c====================================================================
      subroutine cce(nopt,nps,s,sf,bl,bu,xnstd,icall,maxn,iseed)
c$debug
c
c  ALGORITHM GENERATE A NEW POINT(S) FROM A SUB-COMPLEX
c
c  SUB-COMPLEX VARIABLES
      implicit real*8 (a-h,o-z)
      parameter (c1=0.8,c2=0.4)
      dimension s(50,16),sf(50),bu(16),bl(16),xnstd(16)
c
c  LIST OF LOCAL VARIABLES
c    sb(.) = the best point of the simplex
c    sw(.) = the worst point of the simplex
c    w2(.) = the second worst point of the simplex
c    fw = function value of the worst point
c    ce(.) = the centroid of the simplex excluding wo
c    snew(.) = new point generated from the simplex
c    iviol = flag indicating if constraints are violated
c          = 1 , yes
c          = 0 , no
c
      dimension sw(16),sb(16),ce(16),snew(16)
c
c  EQUIVALENCE OF VARIABLES FOR READABILTY OF CODE
      n = nps
      m = nopt
      alpha = 1.0
      beta = 0.5
c
c  IDENTIFY THE WORST POINT wo OF THE SUB-COMPLEX s
c  COMPUTE THE CENTROID ce OF THE REMAINING POINTS
c  COMPUTE step, THE VECTOR BETWEEN wo AND ce
c  IDENTIFY THE WORST FUNCTION VALUE fw
      do j = 1, m
        sb(j) = s(1,j)
        sw(j) = s(n,j)
        ce(j) = 0.0
        do i = 1, n-1
          ce(j) = ce(j) + s(i,j)
        end do
        ce(j) = ce(j)/dble(n-1)
      end do
      fw = sf(n)
c
c  COMPUTE THE NEW POINT snew
c
c  FIRST TRY A REFLECTION STEP
      do j = 1, m
        snew(j) = ce(j) + alpha * (ce(j) - sw(j))
      end do
c
c  CHECK IF snew SATISFIES ALL CONSTRAINTS
      call chkcst(nopt,snew,bl,bu,ibound)
c
c
c  snew IS OUTSIDE THE BOUND,
c  CHOOSE A POINT AT RANDOM WITHIN FEASIBLE REGION ACCORDING TO
c  A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
c  AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
      if (ibound .ge. 1) call getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb)
c
c
c  COMPUTE THE FUNCTION VALUE AT snew
      fnew = functn(nopt,snew)
      icall = icall + 1
c
c  COMPARE fnew WITH THE WORST FUNCTION VALUE fw
c
c  fnew IS LESS THAN fw, ACCEPT THE NEW POINT snew AND RETURN
      if (fnew .le. fw) go to 2000
      if (icall .ge. maxn) go to 3000
c
c
c  fnew IS GREATER THAN fw, SO TRY A CONTRACTION STEP
      do j = 1, m
        snew(j) = ce(j) - beta * (ce(j) - sw(j))
      end do
c
c  COMPUTE THE FUNCTION VALUE OF THE CONTRACTED POINT
      fnew = functn(nopt,snew)
      icall = icall + 1
c
c  COMPARE fnew TO THE WORST VALUE fw
c  IF fnew IS LESS THAN OR EQUAL TO fw, THEN ACCEPT THE POINT AND RETURN
      if (fnew .le. fw) go to 2000
      if (icall .ge. maxn) go to 3000
c
c
c  IF BOTH REFLECTION AND CONTRACTION FAIL, CHOOSE ANOTHER POINT
c  ACCORDING TO A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
c  AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
 1000 call getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb)
c
c  COMPUTE THE FUNCTION VALUE AT THE RANDOM POINT
      fnew = functn(nopt,snew)
      icall = icall + 1
c
c
c  REPLACE THE WORST POINT BY THE NEW POINT
 2000 continue
      do j = 1, m
        s(n,j) = snew(j)
      end do
      sf(n) = fnew
 3000 continue
c
c  END OF SUBROUTINE CCE
      return
      end
c
c
c
c===================================================================
      subroutine getpnt(nopt,idist,iseed,x,bl,bu,std,xi)
c$debug      
c
c     This subroutine generates a new point within feasible region
c
c     x(.) = new point
c     xi(.) = focal point
c     bl(.) = lower bound
c     bu(.) = upper bound
c     std(.) = standard deviation of probability distribution
c     idist = probability flag
c           = 1 - uniform distribution
c           = 2 - Gaussian distribution
c
      implicit real*8 (a-h,o-z)
      dimension x(16),bl(16),bu(16),std(16),xi(16)
c
    1 do j=1, nopt
    2   if (idist .eq. 1) rand = ran1(iseed)
        if (idist .eq. 2) rand = gasdev(iseed)
        x(j) = xi(j) + std(j) * rand * (bu(j) - bl(j))
c
c     Check explicit constraints
c        
        call chkcst(1,x(j),bl(j),bu(j),ibound)
        if (ibound .ge. 1) go to 2
      end do
c
c     Check implicit constraints
c      
      call chkcst(nopt,x,bl,bu,ibound)
      if (ibound .ge. 1) go to 1
c
      return
      end
c
c
c
c===================================================================
      subroutine parstt(npt,nopt,x,xnstd,bound,gnrng,ipcnvg)
c$debug      
c
c  SUBROUTINE CHECKING FOR PARAMETER CONVERGENCE
      implicit real*8 (a-h,o-z)
      dimension x(2000,16),xmax(16),xmin(16)
      dimension xmean(16),xnstd(16),bound(16)
      parameter (delta = 1.0d-20,peps=1.0d-3)
c
c  COMPUTE MAXIMUM, MINIMUM AND STANDARD DEVIATION OF PARAMETER VALUES
      gsum = 0.d0
      do k = 1, nopt
        xmax(k) = -1.0d+20
        xmin(k) = 1.0d+20
        xsum1 = 0.d0
        xsum2 = 0.d0
        do i = 1, npt
          xmax(k) = dmax1(x(i,k), xmax(k))
          xmin(k) = dmin1(x(i,k), xmin(k))
          xsum1 = xsum1 + x(i,k)
          xsum2 = xsum2 + x(i,k)*x(i,k)
        end do
        xmean(k) = xsum1 / dble(npt)
        xnstd(k) = (xsum2 / dble(npt) - xmean(k)*xmean(k))
        if (xnstd(k) .le. delta) xnstd(k) = delta
        xnstd(k) = dsqrt(xnstd(k))
        xnstd(k) = xnstd(k) / bound(k)
        gsum = gsum + dlog( delta + (xmax(k)-xmin(k))/bound(k) )
      end do
      gnrng = dexp(gsum/dble(nopt))
c
c  CHECK IF NORMALIZED STANDARD DEVIATION OF PARAMETER IS <= eps
      ipcnvg = 0
      if (gnrng .le. peps) then
        ipcnvg = 1
      end if
c
c  END OF SUBROUTINE PARSTT
      return
      end
c
c
c
c====================================================================
      subroutine comp(n,npt,ngs1,ngs2,npg,a,af,b,bf)
c$debug      
c
c
c  THIS SUBROUTINE REDUCE INPUT MATRIX a(n,ngs2*npg) TO MATRIX
c  b(n,ngs1*npg) AND VECTOR af(ngs2*npg) TO VECTOR bf(ngs1*npg)
      implicit real*8 (a-h,o-z)
      dimension a(2000,16),af(2000),b(2000,16),bf(2000)
      do igs=1, ngs1
        do ipg=1, npg
          k1=(ipg-1)*ngs2 + igs
          k2=(ipg-1)*ngs1 + igs
          do i=1, n
            b(k2,i) = a(k1,i)
          end do
          bf(k2) = af(k1)
        end do
      end do
c
      do j=1, npt
        do i=1, n
          a(j,i) = b(j,i)
        end do
        af(j) = bf(j)
      end do
c
c  END OF SUBROUTINE COMP
      return
      end
c
c
c
c===================================================================
      subroutine sort(n,m,rb,ra)
c$debug      
c
c
c  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
c  BY W.H. PRESS ET AL., pp. 233-234
c
c  LIST OF VARIABLES
c     ra(.) = array to be sorted
c     rb(.,.) = arrays ordered corresponding to rearrangement of ra(.)
c     wk(.,.), iwk(.) = local varibles
c
      implicit real*8 (a-h,o-z)
      dimension ra(2000),rb(2000,16),wk(2000,16),iwk(2000)
c
      call indexx(n, ra, iwk)
      do 11 i = 1, n
      wk(i,1) = ra(i)
   11 continue
      do 12 i = 1, n
      ra(i) = wk(iwk(i),1)
   12 continue
      do 14 j = 1, m
      do 13 i = 1, n
      wk(i,j) = rb(i,j)
   13 continue
   14 continue
      do 16 j = 1, m
      do 15 i = 1, n
      rb(i,j) = wk(iwk(i),j)
   15 continue
   16 continue
c
c  END OF SUBROUTINE SORT
      return
      end
c
c
c===========================================================
      subroutine sort1(n,ra)
c$debug      
c
c
c  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
c  BY W.H. PRESS ET AL., pp. 231
c
c  LIST OF VARIABLES
c     ra(.) = integer array to be sorted
c
      implicit real*8 (a-h,o-z)
      dimension ra(n)
c
      integer ra, rra
c
      l = (n / 2) + 1
      ir = n
   10 continue
      if (l .gt. 1) then
      l = l - 1
      rra = ra(l)
      else
      rra = ra(ir)
      ra(ir) = ra(1)
      ir = ir - 1
      if (ir .eq. 1) then
      ra(1) = rra
      return
      end if
      end if
      i = l
      j = l + l
   20 if (j .le. ir) then
      if (j .lt. ir) then
      if (ra(j) .lt. ra(j + 1)) j = j + 1
      end if
      if (rra .lt. ra(j)) then
      ra(i) = ra(j)
      i = j
      j = j + j
      else
      j = ir + 1
      end if
      goto 20
      end if
      ra(i) = rra
      goto 10
c
c  END OF SUBROUTINE SORT1
      end
c
c
c
c=======================================================
      subroutine indexx(n, arrin, indx)
c$debug      
c
c
c  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
      implicit real*8 (a-h,o-z)
      dimension arrin(n), indx(n)
c
      do 11 j = 1, n
      indx(j) = j
   11 continue
      l = (n / 2) + 1
      ir = n
   10 continue
      if (l .gt. 1) then
      l = l - 1
      indxt = indx(l)
      q = arrin(indxt)
      else
      indxt = indx(ir)
      q = arrin(indxt)
      indx(ir) = indx(1)
      ir = ir - 1
      if (ir .eq. 1) then
      indx(1) = indxt
      return
      end if
      end if
      i = l
      j = l + l
   20 if (j .le. ir) then
      if (j .lt. ir) then
      if (arrin(indx(j)) .lt. arrin(indx(j + 1))) j = j + 1
      end if
      if (q .lt. arrin(indx(j))) then
      indx(i) = indx(j)
      i = j
      j = j + j
      else
      j = ir + 1
      end if
      goto 20
      end if
      indx(i) = indxt
      goto 10
c
c  END OF SUBROUTINE INDEXX
      end
c
c
c
c==============================================================
      real*8 function ran1(idum)
c$debug
c
c
c  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
      implicit real*8 (a-h,o-z)
      dimension r(97)
      parameter (m1 = 259200, ia1 = 7141, ic1 = 54773, rm1 =
     &3.8580247e-6)
      parameter (m2 = 134456, ia2 = 8121, ic2 = 28411, rm2 =
     &7.4373773e-6)
      parameter (m3 = 243000, ia3 = 4561, ic3 = 51349)
      save
      data iff / 0 /
      if ((idum .lt. 0) .or. (iff .eq. 0)) then
      iff = 1
      ix1 = mod(ic1 - idum,m1)
      ix1 = mod((ia1 * ix1) + ic1,m1)
      ix2 = mod(ix1,m2)
      ix1 = mod((ia1 * ix1) + ic1,m1)
      ix3 = mod(ix1,m3)
      do 11 j = 1, 97
      ix1 = mod((ia1 * ix1) + ic1,m1)
      ix2 = mod((ia2 * ix2) + ic2,m2)
      r(j) = (dble(ix1) + (dble(ix2) * rm2)) * rm1
   11 continue
      idum = 1
      end if
      ix1 = mod((ia1 * ix1) + ic1,m1)
      ix2 = mod((ia2 * ix2) + ic2,m2)
      ix3 = mod((ia3 * ix3) + ic3,m3)
      j = 1 + ((97 * ix3) / m3)
      if ((j .gt. 97) .or. (j .lt. 1)) pause
      ran1 = r(j)
      r(j) = (dble(ix1) + (dble(ix2) * rm2)) * rm1
c
c  END OF SUBROUTINE RAN1
      return
      end
c
c
c
c===============================================================
      real*8 function gasdev(idum)
c$debug
c
c
c  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
      implicit real*8 (a-h,o-z)
      common /gasblk/ iset
      data iset / 0 /
      if (iset .eq. 0) then
    1 v1 = (2. * ran1(idum)) - 1.
      v2 = (2. * ran1(idum)) - 1.
      r = (v1 ** 2) + (v2 ** 2)
      if (r .ge. 1.) goto 1
      fac = sqrt(- ((2. * log(r)) / r))
      gset = v1 * fac
      gasdev = v2 * fac
      iset = 1
      else
      gasdev = gset
      iset = 0
      end if
c
c  END OF SUBROUTINE GASDEV
      return
      end

c
c==================================================================
      subroutine chkcst(nopt,x,bl,bu,ibound)
c
c     This subroutine check if the trial point satisfies all
c     constraints.
c
c     ibound - violation indicator
c            = -1 initial value
c            = 0  no violation
c            = 1  violation
c     nopt = number of optimizing variables
c     ii = the ii'th variable of the arrays x, bl, and bu
c
      implicit real*8 (a-h,o-z)
      dimension x(nopt),bl(nopt),bu(nopt)
c
      ibound = -1
c
c     Check if explicit constraints are violated
c
      do ii=1, nopt
        if (x(ii) .lt. bl(ii) .or. x(ii) .gt. bu(ii)) go to 10
      end do
      if (nopt .eq. 1) go to 9
c
c     Check if implicit constraints are violated
c     (Note - Implicit constraints are commented out, 
c      you can activate these constraints by removing
c      "c" in the first columns)
c      c1 = 5.0*x(1) + 2.0*x(3) - 100.0
c      c2 = x(3) + 50.*x(4) - 50.0
c      c3 = 10.*x(5) - x(6) - 2.0
c      c4 = -10.*x(5) + x(6) - 2.0
c      if (c1 .gt. 0.) go to 10
c      if (c2 .gt. 0.) go to 10
c      if (c3 .gt. 0.) go to 10
c      if (c4 .gt. 0.) go to 10
c
c     No constraints are violated
c
    9 ibound = 0
      return
   10 ibound = 1
      return
      end
c
c
c
c=======================================================================
      real*8 function hmle()
      implicit real*8 (a-h,o-z)
      common /fnblk/ rlamda, ad
      common /block1/ ndata, ns, iobj
      common /block2/ p(1000), obsq(1000), s0(2)
      common /block3/ q(1000), r(1000), b(1000), s(1000)
      dimension ra(2)
      parameter (eps=5.d-02, del=5.d-02)
      data iflag /0/
c
c  COMPUTE THE MEAN OF LOGARITHM OF OBSERVED FLOWS
      if (iflag .eq. 0) then
        ad = 0.d0
        do 10 i = 1, ndata
          ad = ad + dlog(obsq(i))
   10   continue
        ad = ad / dble(ndata)
        rlamda = 1.d0
        iflag = 1
      end if
c
c  ESTIMATE THE LAMDA VALUE
      lcount = 0
      ict = 1
      ra(1) = 0.d0
      ra(2) = 0.d0
   25 continue
      lcount = lcount + 1
      if(lcount .gt. 40) then
        write(*,*) 'LAMDA ITERATION GO OVER 40', rlamda, ra(1), ra(2)
        go to 50
      end if
      rd = 0.d0
      rn = 0.d0
      do 30 i = 1, ndata
        a = dlog(obsq(i)) / ad
        w = obsq(i)**(2*(rlamda-1.d0))
        rd = rd + w*(obsq(i) - q(i))**2
        rn = rn + w*(obsq(i) - q(i))**2 * a
   30 continue
      ra(ict) = rn / rd - 1.d0
      if (dabs(ra(ict)) .le. eps) go to 50
      isign = -1
      if (ra(ict) .lt. 0.d0) isign = 1
      rlamda = rlamda + isign * del
      if (ict .eq. 2) go to 35
      ict = 2
      go to 25
c
   35 continue
      if (ra(1)*ra(2) .lt. 0.d0) go to 40
      ra(1) = ra(2)
      go to 25
c
   40 continue
      rlamda = rlamda - isign * del / 2.d0
c
c  COMPUTE HMLE
   50 continue
      hmle = 0.d0
      ex = 2. * (rlamda - 1.)
      do 60 i = 1, ndata
        hmle = hmle + obsq(i)**ex * (obsq(i) -q(i))**2
   60 continue
      hmle = hmle / dble(ndata)
      hmle = hmle / dexp(ex * ad)
c
      return
      end
