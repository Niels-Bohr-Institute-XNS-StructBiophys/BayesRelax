      program conv
      real x(5000),y(5000),sd(5000)
      character*1 tab
      character*40 aname,bname
      character*60 array(7)
      tab=char(9)
      aname='data.d'
      bname='data2.d'
      open(unit=1,file=aname,status='old')
  4   format(a)
      do 5 k=1,7
      read(1,4)array(k)
  5   continue
      do 10 i=1,5000
      read(1,*,end=11) x(i),y(i),sd(i)
  10  continue
  11  ntot=i-1
      close(1)
      write(6,*)'ntot = ',ntot
      open(unit=2,file=bname,status='unknown')
      do 15 k=1,7
      write(2,4)array(k)
  15  continue
      do 20 j=1,ntot
      write(2,21)x(j),tab,y(j),tab,sd(j)
  20  continue
  21  format(1x,e12.6,a,e12.6,a,e12.6)
      close(2)
      stop 
      end