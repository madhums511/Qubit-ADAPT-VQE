C*
      PROGRAM vorh2u
C*
C*======================================================================
C*
C*    This program changes the parameter for the H_2 program according
C*    to the actual parameters.
C*
C*----------------------------------------------------------------------
C*
      IMPLICIT NONE
C*
      INTEGER i,j,nrmax,nblmax,kmx,lmx,pmx,iin,iout,icheck,icheck1,
     +    icheck2,ismold1,iend,iout2
C*
      CHARACTER*256 filenm
      CHARACTER*80 line
C*
C*---------------------------------------------------------------------
C*

      iin = 1
      iout = 2
      iout2 = 3
C*
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) '               ********************'
      WRITE(*,*) '               *  Program VORH2U  *'
      WRITE(*,*) '               ********************'
C*
      CALL getenv('INPDAT',filenm)
      OPEN ( UNIT=iin, FILE=filenm, STATUS='OLD',
     +       ACCESS='SEQUENTIAL', FORM='FORMATTED' )
C*
      CALL getenv('OUTDAT',filenm)
      OPEN ( UNIT=iout, FILE=filenm, STATUS='NEW',
     +       ACCESS='SEQUENTIAL', FORM='FORMATTED' )
C*
      CALL getenv('OUT2DAT',filenm)
      OPEN ( UNIT=iout2, FILE=filenm, STATUS='NEW',
     +       ACCESS='SEQUENTIAL', FORM='FORMATTED' )
C*

      DO 100 i=1,3
        READ(iin,'(A)') line
  100 CONTINUE
      READ(iin,*) nrmax
      j = MOD(nrmax,10)
      IF (j.NE.0) j=1
      iend = ( nrmax / 10 ) + j + 2
      DO 120 i=1,iend
        READ(iin,'(A)') line
  120 CONTINUE
C*
      READ(iin,140) lmx,pmx,kmx
  140 FORMAT(3I3)
      CLOSE(iin)
C*
      icheck = (nrmax*(nrmax+1)/2) + 800
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) '  icheck = ',icheck
c      IF (icheck.GT.20213) THEN
        ismold1 = 2*(lmx+2)*(pmx+1)**2 + (6*(kmx+1)*(lmx+1)) + 1
        WRITE(*,*) ' ismold1 = ',ismold1
        icheck1 = (2*nrmax*(nrmax+1)) + 2*nrmax
        WRITE(*,*) ' icheck1 = ',icheck1
        icheck2 = ismold1 + ( ( 3*nrmax*(nrmax+1) ) / 2 )
        WRITE(*,*) ' icheck2 = ',icheck2
        IF (icheck1.GT.icheck2) THEN
          nblmax = icheck1
        ELSE
          nblmax = icheck2
        ENDIF
c      ELSE
c        nblmax = 30113 + ( (3*nrmax*(nrmax+1)) / 2 )
c      ENDIF
      WRITE(*,*) ' nblmax = ',nblmax
C*
      WRITE(iout,160)
  160 FORMAT('C*')
      WRITE(iout,170)
  170 FORMAT('C*    Parameter File for the H_2 program:')
      WRITE(iout,160)
      WRITE(iout,200) nrmax,nblmax
  200 FORMAT('      PARAMETER ( nrmax = ',I4,', nblmax = ',I12,
     +   ' )')
      WRITE(iout,160) 
      CLOSE(iout)
C*
      WRITE(iout2,160)
      WRITE(iout2,220) 
  220 FORMAT('C*    2nd Parameter File for the H_2 program:')
      WRITE(iout2,160)
      WRITE(iout2,300) nrmax
  300 FORMAT('      PARAMETER ( nrdim = ',I4,')')
      WRITE(iout2,160)
      CLOSE(iout2)
C*
      WRITE(*,*) ' '
      WRITE(*,*) ' Program VORH2U finished normally. '
      WRITE(*,*) ' '
      END
