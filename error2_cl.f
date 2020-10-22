      PROGRAM dlsmax2
      PARAMETER (NMAX=1000)

      REAL sigma(nmax),X(NMAX),Y(NMAX),SD(NMAX),ysUM(NMAX),sd2(nmax)
      real xi(nmax),yi(nmax),y2i(nmax),y2(nmax)
      real xx(nmax),yy(nmax),sdsd(nmax),yy2(nmax),sdsd2(nmax)
      character*34 bname,cname
      character*32 aname
      integer jmtot(nmax)
C********************************************************************
C     READ Y,SD FROM FILE IN.DAT
C********************************************************************
      character(len=32) arg
      real pi
      PI=3.14159

C********************************************************************

      call getarg(1,arg)
      aname=arg

      bname='sd'//aname
      cname='er'//aname
      OPEN(1,FILE=ANAME,STATUS='UNKNOWN')

      do 10 i=1,100000
      read(1,*,end=11,err=12)xi(i),yi(i),y2i(i)
c      xi(i)=0.000001*xi(i) Eliminated by Peter Ulvskov Oct 2020. 
c      yi(i)=0.000001*yi(i) The reason for dividing by a million is unknown
c      y2i(i)=0.000001*y2i(i)
   10 continue

   12 OPEN(2,FILE=BNAME,STATUS='UNKNOWN')
      write(2,*)'error'
      close(2)
      goto 999

   11 ntot=i-1
   14 close(1)

      j=0
      k=0
      xtest=xi(1)
  13  k=k+1
      jmtot(k)=0

      do 18 i=1,ntot
      if(xi(i).eq.xtest) then
      j=j+1
      x(j)=xi(i)
      y(j)=yi(i)
      y2(j)=y2i(i)
      xi(i)=0
      yi(i)=0
      y2i(i)=0
      jmtot(k)=jmtot(k)+1
      endif
  18  continue
   
      do 19 i=k,ntot
      if(xi(i).ne.0) then 
      xtest=xi(i)
      goto 13  
      endif
  19  continue
      
  20  continue

      do 30 i=1,ntot
      write(55,*)x(i),y(i),y2(i),jmtot(i)
   30 continue
      itot=k

      OPEN(2,FILE=BNAME,STATUS='UNKNOWN')

      k=0
      do 596 i=1,itot
      ym=0
      ym2=0
      sdx=0
      sd2x=0
      do 597 j=1,jmtot(i)  
      k=k+1
      xi(j)=x(k)
      yi(j)=y(k)
      y2i(j)=y2(k)    
      ym=ym+yi(j)
      ym2=ym2+y2i(j)
      write(56,*)xi(j),yi(j),y2i(j),jmtot(i)
  597 continue
      ym=ym/jmtot(i)
      ym2=ym2/jmtot(i)
      do 598 j=1,jmtot(i)
      sdx=sdx+(yi(j)-ym)**2
      sd2x=sd2x+(y2i(j)-ym2)**2
  598 continue
      sdx=sdx/jmtot(i)+0.01*ym**2/(jmtot(i)-0.98)
      sd2x=sd2x/jmtot(i)+0.01*ym2**2/(jmtot(i)-0.98)
      if(ym.lt.0) sdx=sdx*100
      if(ym2.lt.0) sd2x=sd2x*100
      sdx=abs(sdx)
      sd2x=abs(sd2x)
c      sdx=sdx/jmtot(i)+ym**2/(jmtot(i)-0.98)
c      sd2x=sd2x/jmtot(i)+ym2**2/(jmtot(i)-0.98)
      write(2,*)xi(1),ym,sqrt(sdx),ym2,sqrt(sd2x)
      xx(i)=xi(1)
      yy(i)=ym
      sdsd(i)=sqrt(sdx)
      yy2(i)=ym2
      sdsd2(i)=sqrt(sd2x)
  596 continue
      close(2)
C*************
      do 700 i=1,itot
      x(i)=xx(i)
      y(i)=yy(i)
      sd(i)=sdsd(i)
      y2(i)=yy2(i)
      sd2(i)=sdsd2(i)
      xx(i)=0
      yy(i)=0
      sdsd(i)=0
      yy2(i)=0
      sdsd2(i)=0
  700 continue

      do 850 j=1,itot
      xmin=1.e7
      do 800 i=1,itot
      if(x(i).le.xmin) then
      xmin=x(i)
      imin=i
      endif
 800  continue
      xx(j)=x(imin)
      yy(j)=y(imin)
      yy2(j)=y2(imin)
      sdsd(j)=sd(imin)
      sdsd2(j)=sd2(imin)
      x(imin)=1.e10
 850  continue
      
C***************
      OPEN(3,FILE=cNAME,STATUS='UNKNOWN')
      do 600 i=1,itot
      if(xx(i).eq.1.e10) goto 600
      if(jmtot(i).gt.1) write(3,*)xx(i),yy(i),sdsd(i)
      if(jmtot(i).eq.1) write(3,*)xx(i),yy(i),sdsd(i)/4
  600 continue
      do 601 i=1,itot
      if(xx(i).eq.1.e10) goto 601
      if(jmtot(i).gt.1) write(3,*)xx(i),yy2(i),sdsd2(i)
      if(jmtot(i).eq.1) write(3,*)xx(i),yy2(i),sdsd2(i)/4
  601 continue
      close(3)
  999 stop
      end

