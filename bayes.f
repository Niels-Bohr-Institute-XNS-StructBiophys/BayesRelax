      PROGRAM BAYES
C***********************************************************************
C     Details in Hansen (2008)
C***********************************************************************
      PARAMETER (NMAX=1000)
      PARAMETER (NDIST=5000)

C     maximum 1000 points used for data or relaxation spectrum
C     maximum 5000 distributions

      REAL A(NMAX,NMAX),B(NMAX,NMAX),FT(NMAX,NMAX),SMEAR(NMAX,NMAX)
      REAL X(NMAX),Y(NMAX),SD(NMAX),XF(NMAX),F(NMAX),u(nmax,nmax)
      REAL M(NMAX),YSUM(NMAX),FM(NMAX),w(nmax),w2(nmax)
      real sigma(nmax),sigf(nmax)
      real ftot(nmax,ndist),prob(ndist),xmean(ndist),qmin,dot(ndist)
      integer clock,clockold

      common r11,r12,r21,r22,rl11,rl12,rl21,rl22
      common rlam1,rlam2,dlam1,dlam2,ntot1,ntot2
      common /cputime/ clockold,cpumax


      CHARACTER*28 ANAME
      character*22 redaname
      chARACTER*30 BNAME,aaname
      character*32 fit1name,fit2name
      CHARACTER*33 hxname
      CHARACTER*34 KNAME,LNAME
      CHARACTER*35 PNAME,QNAME,in1name,in2name,dataname,diamname
      character*40 xin1name,xin2name
      CHARACTER*36 FNAME
      CHARACTER*6 CHAR
      CHARACTER*1 ANSWER,test
      character*40 dummy
      character(len=32) arg
      real pi
      PI=3.14159

C********************************************************************

c max cpu time in seconds
      cpumax=30
      call system_clock(clock)
      clockold=clock

      call getarg(1,arg)
      aname=arg

      write(6,*)'Name of data file read from command-line'
      write(6,*)'Value for input parameters read from file inputfile.d'

c    1 WRITE(6,3)
c    3 FORMAT(1X,'Write name of the inputfile                     => ',$)
    4 FORMAT(A)
      bname='me'//ANAME

      FNAME='trans_'//BNAME        

      fit1name='fit1'//aname
      fit2name='fit2'//aname
      LNAME='tot_'//BNAME   
      hxname='in_'//bname
      PNAME='prob_'//BNAME   
      OPEN(1,FILE=ANAME,STATUS='UNKNOWN')
      OPEN(40,FILE=PNAME,STATUS='UNKNOWN')
   10 format(1X,'#',3x'alpha',8X,'Chi-square',5X,'   N_g  ',
     -8X,'Evidence',6X,'Check')
      write(40,10)

      OPEN(2,FILE=BNAME,STATUS='UNKNOWN')
      write(6,4)aname
      open(55,file='inputfile.d',status='unknown')

      qmin=0
      WRITE(6,131)
  131 FORMAT(1X,'wmin <enter>                        => ',$)
      read(55,4)dummy
      if(dummy.ne.' ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)qmin
      close(50)
      endif

      qmax=1.e10
      WRITE(6,132)
  132 FORMAT(1X,'wmax or <enter>                        => ',$)
      read(55,4)dummy
      if(dummy.ne.' ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)qmax
      close(50)
      endif

      mtot=0     

      OPEN(1,FILE=aname,STATUS='UNKNOWN')
c************************************************************
c skipping text lines in data file and writing data to dummy
c*************************************************************
      nofm=0
      open(111,file='dummy.d',status='unknown')
      ndata=0
   2  if(ndata.ge.1) goto 8
  22  continue
      read(1,*,err=2,end=8)x1,y1,sd1
      ndata=ndata+1
      write(111,*)x1,y1,sd1
      goto 22
      
   8  rewind(111)
      close(111)
      rewind(111)
       
      close(1)

      if(ndata.lt.350) then
      OPEN(1,FILE='dummy.d',STATUS='unknown')
      x1sum=x1
      goto 7777
      endif
c****************************************************************
c if ndata > 350 reading data from dummy and rebinning into dummy2
c****************************************************************
      nr=ndata/350
      nr2=ndata/nr
      nrest=ndata-nr2*nr

      write(6,*)nr,nr2,nrest

c     nr: rebin no.
c     nr2(+1): number of data points used (after rebinning)
      open(111,file='dummy.d',status='old')
      open(1,file='dummy2.d',status='unknown')

      if(nrest.eq.0) goto 1310
      x1sum=0
      y1sum=0
      sd1sum=0
      do 1311 i=1,nrest
      read(111,*)x1,y1,sd1
      x1sum=x1sum+x1
      y1sum=y1sum+y1
      sd1sum=sd1sum+sd1**2
 1311 continue
      x1sum=x1sum/nrest
      y1sum=y1sum/nrest
      sd1sum=sqrt(sd1sum)/nrest
      write(1,*)x1sum,y1sum,sd1sum

 1310 do 1313 j=1,nr2
      x1sum=0
      y1sum=0
      sd1sum=0
      do 1312 i=1,nr
      read(111,*)x1,y1,sd1
      x1sum=x1sum+x1
      y1sum=y1sum+y1
      sd1sum=sd1sum+sd1**2
 1312 continue
      x1sum=x1sum/nr
      y1sum=y1sum/nr
      sd1sum=sqrt(sd1sum)/nr
      write(1,*)x1sum,y1sum,sd1sum
 1313 continue

      close(1)
      OPEN(1,FILE='dummy2.d',STATUS='old')
c****************************************************************
c Reading data from dummy3 using only points within [q_min;q_max]
c Adding points to constrain S(q) -> 1 for q > qmax/1.6
c The number 1.6 should probably be changed to 2.0 or even larger.
c The error of the added points may also be changed to ensure that
c S(q) -> 1 
c Change the constraints for S(q) here:
c****************************************************************
 7777 close(111)
      

      mtot=0
      sdmax=0.
      do 7 i=1,10000
c    5 READ(1,err=652,end=7)Xi,Yi,sdi   problem here!!!!!!!!!!!!!!
    5 READ(1,*,err=7,end=7)Xi,Yi,sdi
c        write(6,*)i,xi,yi,sdi
        if((xi.ge.qmin).and.(xi.le.qmax)) then
        mtot=mtot+1
        x(mtot)=xi
        Y(mtot)=yi
        sd(mtot)=sdi
        if(sdi.gt.sdmax) sdmax=sdi
c        write(6,*)xi,yi,sdi,mtot
        endif
  798 continue
    7 continue
      close(1)

 652  OPEN(24,FILE=HxNAME,STATUS='UNKNOWN')
      mtot7=mtot
      write(24,*)'#--------------------------------------------------'
      write(24,*)'#                     Data'
      write(24,*)'# Format: (frequency, modulus, standard deviation)'
      write(24,*)'# First part of file: elastic modulus'
      write(24,*)'# Second part of file: loss modulus'
      write(24,*)'# Standard deviation was calculated from raw data'
      write(24,*)'#--------------------------------------------------'
      DO 201 I=1,MTOT
      if(sd(i).eq.0) sd(i)=sdmax
      WRITE(24,*)X(I),Y(I),sd(i)
      if(x(i).le.0.0) mtot7=-1
      if(sd(i).le.0.0) mtot7=-2
  201 CONTINUE
      close(24)

      write(6,*)'Number of points in data file      = ',ndata
      write(6,*)'Number of points used for calc.    = ',mtot
      if(mtot7.le.0) then
      write(6,*)'******************************'
      write(6,*)'****  ERROR IN INPUTFILE  ****'
      write(6,*)'******************************'
      if(mtot7.eq.0) write(6,*) 'No points read'
      if(mtot7.eq.-1) write(6,*) 'Non positive frequencies'
      if(mtot7.eq.-2) write(6,*) 'Non positive standard dev'
      goto 999
      endif
c********************************************************************
c     End of input of data
c********************************************************************
      ALPHA=1.
      CAIM=1.00
c      MAXIT=1000
c      OMEGA=.1
      MAXIT=2000
      OMEGA=.05
   61 CAIMOLD=CAIM
      xmin=1.e9
      xmax=0                  
      DO 600 I=1,MTOT
      SD(I)=SD(I)**2
      if(x(i).le.xmin) xmin=x(i)
      if(x(i).ge.xmax) xmax=x(i)
  600 CONTINUE
  560 ALPHASTART=ALPHA

c********** choice of input parameters

      dmax=(1/xmin)*10
      dmin=(1/xmax)*0.1
c  478 WRITE(6,1113)
c 1113 FORMAT(1X,'tau_min, tau_max, or <enter>                   => ',$)
c      read(55,4)dummy
c      write(6,4)dummy
c      if(dummy.ne.'    ') then
c      open(unit=50,file='dummy.d',status='unknown')
c      write(50,4)dummy
c      close(50)
c      open(unit=50,file='dummy.d',status='unknown')
c      read(50,*)dmin,dmax
c      close(50)
c      endif

      alphamin=log(1.e-4)
      WRITE(6,1114)
 1114 FORMAT(1X,'Alpha_min, or <enter>  => ',$)
      read(55,4)dummy
c      write(6,4)dummy
      if(dummy.ne.' ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      write(6,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)alphamin
      close(50)
      endif

      alphamax=log(1.e2)
      WRITE(6,1115)
 1115 FORMAT(1X,'Alpha_max  <enter>  => ',$)
      read(55,4)dummy
c     write(6,4)dummy
      if(dummy.ne.' ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      write(6,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)alphamax
      close(50)
      endif

      alphamin=exp(alphamin)
      alphamax=exp(alphamax)

      ncalc=30
      WRITE(6,1116)
 1116 FORMAT(1X,'Number of calculations <enter>  => ',$)
      read(55,4)dummy
c     write(6,4)dummy
      if(dummy.ne.' ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
c      write(6,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)ncalc
      close(50)
      endif
c     write(6,4)dummy

      alpharel=exp(log(alphamin/alphamax)/ncalc)

      ntot0=100      
      ntot=30
      WRITE(6,113)
  113 FORMAT(1X,'Minimum number of points in p(r) or <enter>     => ',$)
c      write(6,4)dummy
      read(55,4)dummy
      if(dummy.ne.' ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)ntot
      close(50)
      endif
      close(66)

      WRITE(6,114)
  114 FORMAT(1X,'Free end points enter Yes (default No)          => ',$)
c      write(6,4)dummy
      read(55,4)dummy
      nfree=0
      if(dummy.ne.' ') then
      nfree=1
      endif

   40 FORMAT(1X,' No ',1X,'  N_g   ',2X,'Alpha',3X,'Chi-square',
     -2X,'Constr.',4X,'Check',5X,'Sum',5x,'-Evidence')
C*****************************************
c     Start variation of diameter
c*****************************************

      nof=0
 555  continue

      CALL PRIOR(M,NTOT,X,Y,XF,MTOT,DFX,NMAX,dmin,dmax)
      CALL TRANS(A,FT,SMEAR,X,XF,NTOT,MTOT,DFX,NMAX)

      DO 48 I=1,NTOT
      ysum(i)=0.
      DO 47 K=1,MTOT
      YSUM(I)=YSUM(I)+Y(K)*A(K,I)/SD(K)
   47 CONTINUE
   48 CONTINUE
 
      DO 14 I=1,NTOT
      DO 13 J=1,NTOT
        BIJ=0
          DO 12 K=1,MTOT
          BIJ=BIJ+A(K,I)*A(K,J)/SD(K)
   12     CONTINUE
        B(I,J)=BIJ
   13 CONTINUE
   14 CONTINUE

 1399 C1=0
      C2=0
      SUMM=0
      DO 1411 I=1,3 
      FM(I)=0
      DO 1400 J=1,NTOT
      FM(I)=FM(I)+A(I,J)*M(J)
 1400 CONTINUE
      C1=C1+FM(I)/SD(I)
      C2=C2+Y(I)/SD(I)
 1411 CONTINUE
      SUMM=C2/C1
 
      DO 35 j=1,NTOT
c      M(j)=SUMM*M(j)
      F(j)=M(j)*1.0111                     
      sigma(j)=1
c      sigma(j)=dfx**2
   35 CONTINUE

      evmax=-1.e10

      alpha=alphamax
 3333 continue

c**************************************************************
c     Start variation of alpha
c**************************************************************

      nof=nof+1
      alpha=alpha*alpharel

      ITE=1
      if(nof.eq.1) WRITE(6,40)
c      goto 110

C*************************************************************************
C     Calculation of p(r) for given diameter and alpha
C*************************************************************************

   50 ITE=ITE+1

      SUMM=0.                              
      SUMF=0.
      DO 63 J=1,NTOT
      SUMM=SUMM+M(J)
      SUMF=SUMF+F(J)
      sigma(j)=m(j)
   63 CONTINUE
      SUMFACT=SUMF/SUMM
      m(1)=f(2)/2
      m(ntot)=f(ntot-1)/2
c    Below: BSW-spectra      Above: end points = 0
      if(nfree.eq.1) then
      m(1)=f(2)*f(3)/f(4)
      m(ntot)=f(ntot-1)*f(ntot-2)/f(ntot-3)
      endif
      DO 64 J=2,NTOT-1
      m(j)=(f(j-1)+f(j+1))/2
   64 CONTINUE
      
   70 DO 100 I=1,NTOT
      FSUMI=0
      DO 75 J=1,NTOT
      FSUMI=FSUMI+B(I,J)*F(J)
   75 CONTINUE
      FSUMI=FSUMI-B(I,I)*F(I)      


      EE=2*(YSUM(I)-FSUMI)/ALPHA
      BB=M(I)*2*B(I,I)/ALPHA
      FX=ALPHA/(2*B(I,I))
      fx=fx*xx(bb,ee)

ce      fx=(alpha*m(i)/sigma(i)+ysum(i)-fsumi)/(alpha/sigma(i)+b(i,i))
      F(I)=(1.-OMEGA)*F(I)+OMEGA*FX
  100 CONTINUE

c**************************************************************************

  110 S=0
      C=0
      GRADSI=0
      WGRADS=0
      GRADCI=0
      WGRADC=0
      DOTSP=0

      DO 150 I=1,NTOT
ce      SADD=-(f(i)-m(i))**2/sigma(i)      
      SADD=-F(I)*LOG(F(I)/M(I))+F(I)-M(I)
      S=S+SADD      
      GRADSI=-LOG(F(I)/M(I))
ce      GRADSI=-2*(f(i)-m(i))/sigma(i)

      WGRADS=WGRADS+GRADSI**2
      FSUMI=0
      DO 130 J=1,NTOT
      FSUMI=FSUMI+B(I,J)*F(J) 
  130 CONTINUE
      GRADCI=2*(FSUMI-YSUM(I))
      WGRADC=WGRADC+GRADCI**2
      DOTSP=DOTSP+GRADSI*GRADCI
  150 CONTINUE

      DO 141 I=1,MTOT
      FM(I)=0
      DO 140 J=1,NTOT
      FM(I)=FM(I)+A(I,J)*F(J)
  140 CONTINUE
      CADD=(Y(I)-FM(I))**2/SD(I)
      C=C+CADD
  141 CONTINUE
      C=C/MTOT

      call system_clock(clock)
      cpu=(clock-clockold)*0.001

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if((cpu.ge.cpumax).or.(ndata.eq.0)) then

 777   open(unit=88,file='plot1.pl',status='unknown')
	write(88,*)'set term gif'
	write(88,*)'set nokey'
	write(88,*)'set output "fig1.gif"'
       if(ndata.eq.0) then
       write(88,*)
     -'set label "Error in inputfile - no data read. " at -6,0.55'
       else
       write(88,*)
     -'set label "Too many iterations or max cpu time." at -6,0.55'
       write(88,*)
     -'set label "Try reducing the number of points." at -4.5,0.45'
       endif
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88)

      open(unit=88,file='plot2.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
       if(ndata.eq.0) then
       write(88,*)
     -'set label "Error in inputfile - no data read. " at -6,0.55'
       else
	write(88,*)'set output "fig2.gif"'
       write(88,*)
     -'set label "Too many iterations or max cpu time." at -6,0.55'
       write(88,*)
     -'set label "Try reducing the number of points." at -4.5,0.45'
       endif
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88) 
      goto 9999

      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      IF((WGRADS.EQ.0).OR.(WGRADC.EQ.0)) THEN
      WRITE(6,*)'ITE    = ',ITE
      WRITE(6,*)'WGRADS = ',WGRADS
      WRITE(6,*)'WGRADC = ',WGRADC
      GOTO 50
      eNDIF
      WGRADS=SQRT(WGRADS)
      WGRADC=SQRT(WGRADC)
      DOTSP=DOTSP/(WGRADS*WGRADC)
 1798 continue
        if(ite.ge.20) then
        if(mod(ite,100).eq.0) goto 160
        goto 162
        endif
  160  continue
  162 continue 
c 1621 IF((ITE.EQ.MAXIT).or.(dotsp.ge.0.9999)) GOTO 1799
 1621 IF(ITE.EQ.MAXIT) GOTO 1799
      goto 50

C************************************************************************
C     Estimate of p(r) written to file
C     and the probability for this particular solution is calculated
C************************************************************************

 1799 continue
      dot(nof)=dotsp

      alpha=alpha/(2*4)
      S=2*S*4
      do 1801 i=1,ntot
      do 1797 j=1,ntot
      u(i,j)=sqrt(f(i)*f(j))*b(i,j)/alpha
      if(abs(i-j).eq.0) u(i,j)=u(i,j)+2
      if(abs(i-j).eq.1) u(i,j)=u(i,j)-1
 1797 continue
 1801 continue 
      call SVDCMP(u,ntot,Ntot,nmax,Nmax,W,smear)
 1813 rlogdet=0
      det=1.
      do 1802 i=1,ntot
      w2(i)=w(i)*alpha
      det=det*w(i)
      rlogdet=rlogdet+log(abs(w(i)))  
 1802 continue
      rnorm=0.5*(log(ntot+1.)+ntot*log(0.5))

      evidence=rnorm+(alpha*S-0.5*c*mtot*caimold)-0.5*rlogdet-log(alpha)

      if(dotsp.lt.0.9) evidence=-1.e10

      if(evidence.gt.evmax) then
      nofm=nof
      evmax=evidence
      chi2=c*caimold
      endif

c      write(35,995)
c     -alpha*s,-c*mtot*caimold/2.,-rlogdet/2.,
c     -rnorm,-log(alpha),evidence
c  995 format(7(1x,e10.4)) 
      do 997 j=1,ntot
      ftot(j,nof)=f(j)
  997 continue
      do 9971 j=ntot+1,nmax-4
      ftot(j,nof)=0.
 9971 continue
      prob(nof)=evidence
      ftot(nmax-1,nof)=dmax
      ftot(nmax-2,nof)=c*caimold
      ftot(nmax-3,nof)=alpha*8

      do 1777 i=1,ntot
      do 1778 j=1,ntot
      u(i,j)=sqrt(f(i)*f(j))*b(i,j)
 1778 continue
 1777 continue

      call SVDCMP(u,ntot,Ntot,nmax,Nmax,W,smear)
      sum=0
      sumx=0
      do 850 i=1,ntot
      sum=sum+w(i)/(w(i)+alpha)
      sumx=sumx+w(i)/w2(i)
  850 continue
      ftot(nmax-5,nof)=0
      ftot(nmax-6,nof)=-2*alpha*S
      ftot(nmax-4,nof)=sumx
   
      fmax=0
      do 852 i=1,ntot
      if(f(i).ge.fmax) fmax=f(i)
  852 continue 
      nstart=ntot/2
      do 853 i=nstart,ntot
      if(f(i).le.(0.01*fmax)) then
      ftot(nmax-7,nof)=xf(i)
      df=(f(i-1)-f(i)) / (xf(i)-xf(i-1))
      ftot(nmax-8,nof)=xf(i)+f(i)/df
      goto 854
      endif
  853 continue
      ftot(nmax-7,nof)=xf(ntot)
      ftot(nmax-8,nof)=xf(ntot)
  854 continue
     
changed *8******************************''
      alpha=alpha*2.*4
      S=S/(2.*4)

  161 FORMAT(1X,I4,1x,f6.2,1X,5(1X,E9.4),1x,e10.4)
  180 WRITE(6,161)nof,sum,ALPHA,c,-S*alpha*2,DOTSP,
     -ftot(nmax-4,nof-1),-evidence

      if(alpha.lt.alphamin) goto 998
      goto 3333
c************************************
c     End of variation of alpha
c************************************
      
  998 continue
  300 continue

      write(6,*)sum,ftot(nmax-2,nof),ftot(nmax-3,nof)      
      write(6,*)

      evimax=-1.e20
      do 9991 k=1,nof
      if(prob(k).gt.evimax) evimax=prob(k)
 9991 continue

      do 9992 k=1,nof
      prob(k)=exp(prob(k)-evimax)
 9992 continue

ccccccccccccccccccccccc

      summodel=0
      do 9993 k=1,nof
      write(40,*)
     -ftot(nmax-3,k),ftot(nmax-2,k),ftot(nmax-4,k),prob(k),dot(k)
      summodel=summodel+prob(k)
 9993 continue


      sand=prob(1)
      do 5541 k=1,nof-1
      if(ftot(nmax-1,k).eq.ftot(nmax-1,k+1)) sand=sand+prob(k+1)
      if(ftot(nmax-1,k).ne.ftot(nmax-1,k+1)) then
      sand=prob(k+1)
      endif
 5541 continue

      totprob=0.
      do 9118 k=1,nof
      totprob=totprob+prob(k)
 9118 continue

C*****************************************************************************
C     Output estimate of H(tau) - including errors
C*****************************************************************************
      write(2,*)'#-----------------------------------------------'
      write(2,*)'#                 Relaxation function'
      write(2,*)'# Format: (tau, H(tau), standard deviation)'
      write(2,*)'#-----------------------------------------------'
      do 9995 j=1,ntot

      xmean(j)=0
      do 9994 k=1,nof
      xmean(j)=xmean(j)+ftot(j,k)*prob(k)
9994  continue
      xmean(j)=xmean(j)/totprob

      std=0.
      do 9321 k=1,nof
      std=std+prob(k)*(ftot(j,k)-xmean(j))**2
 9321 continue
      std=sqrt(std/totprob)

      f(j)=xmean(j)
      sigf(j)=std
      write(2,*) xf(j),f(j),sigf(j)
 9995 continue
      close(2)
c******************************************************
c     Output
c******************************************************
      alphamean=0
      sumng=0
      sumng2=0
      sng=0
      probtot=0
      cav=0
      cmax=0
      do 1132 k=1,nof
      alphamean=alphamean+log(ftot(nmax-3,k))*prob(k)
      cav=cav+ftot(nmax-2,k)*prob(k)
      if(prob(k).eq.1) cmax=ftot(nmax-2,k)
      sumng=sumng+ftot(nmax-4,k)*prob(k)
      sumng2=sumng2+ftot(nmax-6,k)*prob(k)
      sng=sng+ftot(nmax-5,k)*prob(k)
      probtot=probtot+prob(k)
 1132 continue
      sumng=sumng/probtot
      sumng2=sumng2/probtot
      cav=cav/probtot
      alphamean=exp(alphamean/probtot)

      chisd=0
      do 1333 k=1,nof
      chisd=chisd+prob(k)*(cav-ftot(nmax-2,k))**2
 1333 continue
      chisd=sqrt(chisd/probtot)
c***********************************************************
c     Number of good parameters calculated for optimum alpha
c***********************************************************
      do 7775 i=1,ntot
      do 7778 j=1,ntot
      u(i,j)=sqrt(abs(f(i)*f(j)))*b(i,j)
 7778 continue
 7775 continue
      call SVDCMP(u,ntot,Ntot,nmax,Nmax,W,smear)
 
      do 8777 i=1,ntot
      do 8778 j=1,ntot
      u(i,j)=b(i,j)
 8778 continue
 8777 continue
  
      call SVDCMP(u,ntot,Ntot,nmax,Nmax,W2,smear)
    
      write(6,*)'<Ng> , SNg= ', sumng,sumng2
      write(6,*)'<Chi> , Chi_max = ', cav,cmax
      write(6,*)'Dotsp = ',dot(nofm)
 1133 format(a,f6.3,1x,f5.3,x,f5.5,x,f5.2,
     -x,f5.2,1x,f6.3,1x,f6.3,x,f6.2,x,f6.2)
      write(6,*)

c******************************************************
c     Output fit of data
c******************************************************

      OPEN(11,FILE=FNAME,STATUS='UNKNOWN')
      OPEN(51,FILE=Fit1NAME,STATUS='UNKNOWN')
      OPEN(52,FILE=Fit2NAME,STATUS='UNKNOWN')
      write(11,*)'#-------------------------------------------------'
      write(11,*)'#                     Fit of data'
      write(11,*)'# Format: (frequency, modulus)'
      write(11,*)'# First part of file: elastic modulus'
      write(11,*)'# Second part of file: loss modulus'
      write(11,*)'#-------------------------------------------------'
      DO 8141 I=1,MTOT
      FM(I)=0
      DO 8140 J=1,NTOT
      FM(I)=FM(I)+A(I,J)*F(J)
 8140 CONTINUE
      WRITE(11,*)X(I),FM(I),0.
      if(i.le.mtot/2) write(51,*)x(i),fm(i)
      if(i.gt.mtot/2) write(52,*)x(i),fm(i)
 8141 CONTINUE
      close(51)
      close(52)
      close(11)

  251 continue
      close(43)
      close(14)
C***************************************************************************** 

  899 WRITE(6,900)BNAME
  900 FORMAT(1X,'Estimated function in :       ',A)
      WRITE(6,901)FNAME
  901 FORMAT(1X,'Fit of data in        : ',A)
      WRITE(6,910)PNAME
  910 FORMAT(1X,'Prob.dist. in         : ',A)
      CLOSE(11)
      CLOSE(14)
      CLOSE(20) 
      CLOSE(36) 
      CLOSE(37) 
  999 call plot(dmin,dmax,mtot,f,fm,aname,cav,alphamean,x,dot(nofm))
 9999 call system_clock(clock)
      cpu=(clock-clockold)*0.001
      write(6,955)cpu,cpumax
  955 format(1x,'Cpu time used: ',f4.1,
     -' seconds of maximum ',f4.1,' seconds')
      stop
      END

C*************************************************************************
C     Solution of the equation XX*EXP(XX)=BB*EXP(EE)=AA
C*************************************************************************
      FUNCTION XX(BB,EE)
      IF(EE.GT.50) GOTO 85
      if(ee.lt.-50) then
      xx=0
      goto 90
      endif
      AA=BB*EXP(EE)
      IF(AA.LE.0.05) XX=-.5+SQRT(.25+AA)
      IF(AA.LE.0.05) GOTO 90
      XO=LOG(AA+1.)
   80 FXO=XO*EXP(XO)
      XX=XO*(1-(FXO-AA)/(FXO*XO+1))
      IF(XX.LT.0) XX=0
      FX=XX*EXP(XX)
      DF=ABS(FX-AA)/AA
      IF(DF.LT.0.001) GOTO 90
      XO=XX
      GOTO 80
 
   85 CORR=EE+LOG(BB)
      XX=CORR
   86 CORX=XX+LOG(XX)
      TEST=ABS((CORR-CORX)/CORR)
      IF(TEST.LT.0.0005) GOTO 90
      XO=XX
      XX=XO+(CORR-CORX)*(1+1./XO)
      GOTO 86
 
   90 RETURN
      END
**************************************************************************
********************************************************************
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      PARAMETER (NMAX=1000)
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      G=0.0
      SCALE=0.0
      ANORM=0.0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0
        S=0.0
        SCALE=0.0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0
        S=0.0
        SCALE=0.0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0
            V(J,I)=0.0
31        CONTINUE
        ENDIF
        V(I,I)=1.0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0) THEN
          G=1.0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
C          IF (ITS.EQ.300) PAUSE 'No convergence in 300 iterations'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=SQRT(F*F+1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0) THEN
              Z=1.0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END
**************************************************************************
*     Definition of the prior
**************************************************************************
      SUBROUTINE PRIOR(M,NTOT,X,Y,XF,MTOT,DFX,NMAX,dmin,dmax)

      real X(NMAX),Y(NMAX),XF(NMAX),M(NMAX)
      REAL PMAX,PPMIN
      CHARACTER*20 NAME
      CHARACTER*1 ANS

      DFX=DMAX/NTOT
      w=(dmin/dmax)**(1./(ntot-1))                     

      DO 12 J=1,NTOT 
      XF(J)=DMAX*w**(J-1) 
  12  CONTINUE 

      do 444 j=1,ntot
      m(j)=1+100*(j-j/2)**2/ntot**2
  444 continue

 999  RETURN
      END
C*********************************************************************
C     SETTING UP OF TRANSFORMATION MATRIX A    ( A=Smear*FT )
C     (Calculation of smearing and/or Fourier transformation)
C*********************************************************************
      SUBROUTINE TRANS(A,FT,SMEAR,x,xf,NTOT,MTOT,DFX,NMAX)
      REAL A(NMAX,NMAX),SMEAR(NMAX,NMAX),FT(NMAX,NMAX)
      REAL r(NMAX),q(NMAX),x(nmax),xf(nmax),xj(nmax)
      common r11,r12,r21,r22,rl11,rl12,rl21,rl22
      common rlam1,rlam2,dlam1,dlam2,ntot1,ntot2
      real lambda, omega

      itot=mtot/2
      dr=log(abs((xf(1)/xf(2))))
      
      DO 803 I=1,itot
      DO 802 J=1,NTOT
      lambda=xf(j)
      omega=x(i)
     
      a(i,j)=(omega*lambda)**2/(1+(omega*lambda)**2)*dr
      a(i+itot,j)=(omega*lambda)/(1+(omega*lambda)**2)*dr

  802 CONTINUE
  803 CONTINUE


  999 return
      end
c**********************************************************************
c     This makes a plot file to be used by wgnuplot
c**********************************************************************
      subroutine plot(dmin,dmax,mtot,f,fm,oldfile2,chi,alpha,x,dotspm)      
      PARAMETER (NMAX=1000)

      real xf(nmax),f(nmax),fm(nmax),x(nmax)
      
c      CHARACTER*18 ANAME
      character*18 oldfile2
      character*12 oldfile
c      CHARACTER*20 BNAME
      CHARACTER*6 ANAME
      CHARACTER*10 BNAME
      character*20 comment
 
      CHaRACTER*8,pl1name,pl2name
c      CHaRACTER*21,pl1name,pl2name
c      CHaRACTER*22,fit1name,fit2name
      character*6 fit1name, fit2name
      CHARACTER*23 hxname
      CHARACTER*24 KNAME
      CHARACTER*25 PNAME
      CHARACTER*26 FNAME

      oldfile=oldfile2

c      aname=oldfile


c      bname='me'//ANAME
      bname='estimate.d'
      aname='data.d'

      if(dotspm.ge.0.9) comment='conv. ok'

      if(dotspm.le.0.9) comment='conv. NOT ok'

c      FNAME='trans_'//BNAME        
c      hxname='in_'//bname
c      PNAME='prob_'//BNAME   
c      pl1name='pl1'//aname
c      pl2name='pl2'//aname
c      fit1name='fit1'//aname
c      fit2name='fit2'//aname
      fit1name='fit1.d'
      fit2name='fit2.d'
      pl1name='plot1.pl'
      pl2name='plot2.pl'

      open(unit=88,file=pl1name,status='unknown')

	write(88,*)'set noclip points'
	write(88,*)'set clip one'
	write(88,*)'set noclip two'
	write(88,*)'set bar 1.000000'
	write(88,*)'set border 31 lt -1 lw 1.000'
	write(88,*)'set xdata'
	write(88,*)'set ydata'
	write(88,*)'set zdata'
	write(88,*)'set x2data'
	write(88,*)'set y2data'
	write(88,*)'set boxwidth'
	write(88,*)'set dummy x,y'
	write(88,*)'set format x "%g"'
	write(88,*)'set format y "%g"'
	write(88,*)'set format x2 "%g"'
	write(88,*)'set format y2 "%g"'
	write(88,*)'set format z "%g"'
	write(88,*)'set angles radians'
	write(88,*)'set nogrid'
	write(88,*)'set key title ""'
c	write(88,*)'set key right top Right noreverse box linetype -2 linewidth 1.000          c     -samplen 4 spacing 1 width 0'
	write(88,*)'set nolabel'
	write(88,*)'set noarrow'
c	write(88,*)'unset style line'
	write(88,*)'set nologscale'
	write(88,*)'set logscale x 10'
	write(88,*)'set offsets 0, 0, 0, 0'
	write(88,*)'set pointsize 1'
	write(88,*)'set encoding default'
	write(88,*)'set nopolar'
	write(88,*)'set noparametric'
	write(88,*)'set view 60, 30, 1, 1'
	write(88,*)'set samples 100, 100'
	write(88,*)'set isosamples 10, 10'
	write(88,*)'set surface'
	write(88,*)'set nocontour'
	write(88,*)'set clabel "%8.3g"'
	write(88,*)'set mapping cartesian'
	write(88,*)'set nohidden3d'
	write(88,*)'set cntrparam order 4'
	write(88,*)'set cntrparam linear'
	write(88,*)'set cntrparam levels auto 5'
	write(88,*)'set cntrparam points 5'
	write(88,*)'set size ratio 0 1,1'
	write(88,*)'set origin 0,0'
	write(88,*)'set data style points'
	write(88,*)'set function style lines'
	write(88,*)'set xzeroaxis lt -2 lw 1.000'
	write(88,*)'set x2zeroaxis lt -2 lw 1.000'
	write(88,*)'set yzeroaxis lt -2 lw 1.000'
	write(88,*)'set y2zeroaxis lt -2 lw 1.000'
	write(88,*)'set tics in'
	write(88,*)'set ticslevel 0.5'
	write(88,*)'set ticscale 1 0.5'
	write(88,*)'set mxtics default'
	write(88,*)'set mytics default'
	write(88,*)'set mx2tics default'
	write(88,*)'set my2tics default'
	write(88,*)'set xtics border mirror norotate autofreq '
	write(88,*)'set ytics border mirror norotate autofreq '
	write(88,*)'set ztics border nomirror norotate autofreq' 
	write(88,*)'set nox2tics'
	write(88,*)'set noy2tics'
c	write(88,*)'set title "" 0.000000,0.000000  ""'
C	write(88,*)'set timestamp "" bottom norotate 0.000000,0.000000  ""'
	write(88,*)'set rrange [ * : * ] noreverse nowriteback  '
	write(88,*)'set trange [ * : * ] noreverse nowriteback  '
	write(88,*)'set urange [ * : * ] noreverse nowriteback  '
	write(88,*)'set vrange [ * : * ] noreverse nowriteback  '
c	write(88,*)'set xlabel "" 0.000000,0.000000  ""'
c	write(88,*)'set x2label "" 0.000000,0.000000  ""'
c	write(88,*)'set timefmt "%d/%m/%y\n%H:%M"'
	write(88,*)'set xrange [ * : * ] noreverse nowriteback  '
	write(88,*)'set x2range [ * : * ] noreverse nowriteback '
c	write(88,*)'set ylabel "" 0.000000,0.000000  '
c	write(88,*)'set y2label "" 0.000000,0.000000  ""'
	write(88,*)'set yrange [ * : * ] noreverse nowriteback  '
	write(88,*)'set y2range [ * : * ] noreverse nowriteback  '
c	write(88,*)'set zlabel "" 0.000000,0.000000  ""'
	write(88,*)'set zrange [ * : * ] noreverse nowriteback '
	write(88,*)'set zero 1e-008'
	write(88,*)'set lmargin -1'
	write(88,*)'set bmargin -1'
	write(88,*)'set rmargin -1'
	write(88,*)'set tmargin -1'
	write(88,*)'set locale "C"'
       write(88,*)'set nokey'
c       write(88,*)'set logscale y'
       xmin=x(1)
       xmax=x(mtot)
       if(x(1).ge.x(mtot))then
       xmax=x(1)
       xmin=x(mtot)
       endif
        write(88,444)oldfile,xmin,xmax
 444  format(1x,'set title "Data and model fit for ',a,'
     -from ',f10.4,' to ',f10.4,'"')

       if(mtot.eq.0) then
       write(88,*)'set title "ERROR IN INPUTFILE"'
       write(88,*)'set label "ERROR IN INPUTFILE" at graph 0.41,0.5'
       endif

        write(88,*)'set xlabel "omega [Hz]"'
        write(88,*)'set ylabel "G'' and G''''"'
        write(88,*)'set term gif'
        write(88,*)'set output "fig1.gif"'
        write(88,888)aname,fit1name,fit2name
 888  format(1x,'plot ''',a,''' w e 1,''',a,''' w l 3,''',a,''' w l 3')        
 999  write(88,*)'exit gnuplot'
      close(88)


      open(unit=88,file=pl2name,status='unknown')

	write(88,*)'set noclip points'
	write(88,*)'set clip one'
	write(88,*)'set noclip two'
	write(88,*)'set bar 1.000000'
	write(88,*)'set border 31 lt -1 lw 1.000'
	write(88,*)'set xdata'
	write(88,*)'set ydata'
	write(88,*)'set zdata'
	write(88,*)'set x2data'
	write(88,*)'set y2data'
	write(88,*)'set boxwidth'
	write(88,*)'set dummy x,y'
	write(88,*)'set format x "%g"'
	write(88,*)'set format y "%g"'
	write(88,*)'set format x2 "%g"'
	write(88,*)'set format y2 "%g"'
	write(88,*)'set format z "%g"'
	write(88,*)'set angles radians'
	write(88,*)'set nogrid'
	write(88,*)'set key title ""'
c	write(88,*)'set key right top Right noreverse box linetype -2 linewidth 1.000          c     -samplen 4 spacing 1 width 0'
	write(88,*)'set nolabel'
	write(88,*)'set noarrow'
c	write(88,*)'unset style line'
	write(88,*)'set nologscale'
	write(88,*)'set logscale x 10'
	write(88,*)'set offsets 0, 0, 0, 0'
	write(88,*)'set pointsize 1'
	write(88,*)'set encoding default'
	write(88,*)'set nopolar'
	write(88,*)'set noparametric'
	write(88,*)'set view 60, 30, 1, 1'
	write(88,*)'set samples 100, 100'
	write(88,*)'set isosamples 10, 10'
	write(88,*)'set surface'
	write(88,*)'set nocontour'
	write(88,*)'set clabel "%8.3g"'
	write(88,*)'set mapping cartesian'
	write(88,*)'set nohidden3d'
	write(88,*)'set cntrparam order 4'
	write(88,*)'set cntrparam linear'
	write(88,*)'set cntrparam levels auto 5'
	write(88,*)'set cntrparam points 5'
	write(88,*)'set size ratio 0 1,1'
	write(88,*)'set origin 0,0'
	write(88,*)'set data style points'
	write(88,*)'set function style lines'
	write(88,*)'set xzeroaxis lt -2 lw 1.000'
	write(88,*)'set x2zeroaxis lt -2 lw 1.000'
	write(88,*)'set yzeroaxis lt -2 lw 1.000'
	write(88,*)'set y2zeroaxis lt -2 lw 1.000'
	write(88,*)'set tics in'
	write(88,*)'set ticslevel 0.5'
	write(88,*)'set ticscale 1 0.5'
	write(88,*)'set mxtics default'
	write(88,*)'set mytics default'
	write(88,*)'set mx2tics default'
	write(88,*)'set my2tics default'
	write(88,*)'set xtics border mirror norotate autofreq '
	write(88,*)'set ytics border mirror norotate autofreq '
	write(88,*)'set ztics border nomirror norotate autofreq' 
	write(88,*)'set nox2tics'
	write(88,*)'set noy2tics'
c	write(88,880)aname
c  880 format(1x,'set label ''',a,''')
C at graph 1.3,1.4')
c	write(88,*)'set title "" 0.000000,0.000000  ""'
C	write(88,*)'set timestamp "" bottom norotate 0.00,0.00  ""'
	write(88,*)'set rrange [ * : * ] noreverse nowriteback  '
	write(88,*)'set trange [ * : * ] noreverse nowriteback  '
	write(88,*)'set urange [ * : * ] noreverse nowriteback  '
	write(88,*)'set vrange [ * : * ] noreverse nowriteback  '
c	write(88,*)'set xlabel "" 0.000000,0.000000  ""'
c	write(88,*)'set x2label "" 0.000000,0.000000  ""'
c	write(88,*)'set timefmt "%d/%m/%y\n%H:%M"'
	write(88,*)'set xrange [ * : * ] noreverse nowriteback  '
	write(88,*)'set x2range [ * : * ] noreverse nowriteback '
c	write(88,*)'set ylabel "" 0.000000,0.000000  '
c	write(88,*)'set y2label "" 0.000000,0.000000  ""'
	write(88,*)'set yrange [ * : * ] noreverse nowriteback  '
	write(88,*)'set y2range [ * : * ] noreverse nowriteback  '
c	write(88,*)'set zlabel "" 0.000000,0.000000  ""'
	write(88,*)'set zrange [ * : * ] noreverse nowriteback '
	write(88,*)'set zero 1e-008'
	write(88,*)'set lmargin -1'
	write(88,*)'set bmargin -1'
	write(88,*)'set rmargin -1'
	write(88,*)'set tmargin -1'
       write(88,*)'set nokey'
c       write(88,*)'set logscale y'
c       write(88,222)
c 222   format(1x,'set label "chi" at 0.5,0.5')
	write(88,*)'set locale "C"'
        write(88,544)oldfile,chi,log(alpha),comment
 544  format(1x,'set title "  Relax. ',
     -'spectrum for ',a,'  chi = ',f6.3,' (log(alpha) =',f6.1,
     -') ',a,'"')

       if(mtot.eq.0) then
       write(88,*)'set title "ERROR IN INPUTFILE"'
       write(88,*)'set label "ERROR IN INPUTFILE" at graph 0.41,0.5'
       endif

        write(88,*)'set xlabel "tau [s]"'
        write(88,*)'set ylabel "H(tau)"'
        write(88,*)'set term gif'
        write(88,*)'set output "fig2.gif"'
	write(88,588)bname,bname
 588  format(1x,'plot ''',a,''' w e 1,''',a,''' w l 3')
 9991 write(88,*)'exit gnuplot'
      close(88)

      write(6,*)'Gnuplot files (=> gif): plot1.pl and plot2.pl'


      RETURN
      END
