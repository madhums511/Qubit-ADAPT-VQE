C*
C*==================================================================
C*
      PROGRAM hhsp02
C*
C*------------------------------------------------------------------
C*
C*    This program performs a calculation of the eigenstates of an 
C*    exotic diatomic molecule as the hydrogen/antihydrogen system 
C*    containing two leptons (more specifically: one electron and 
C*    one positron). The two baryons are supposed to have also 
C*    oppositely signed charges, i.e. one has a positive and one 
C*    a negative charge. 
C*
C*    In contrast to the Hamiltonian of normal diatomic molecules
C*    (H_2, HeH^+, etc.) the interleptonic repulsion is exchanged
C*    to an attraction (opposite sign). In addition, the interaction
C*    between one baryon and the two leptons is supposed to be once 
C*    attractive and once repulsive.
C* 
C*    This is a modified version of the H_2 program written by 
C*    Kolos (1985-87) and already earlier modified by his 
C*    collaborators. However, it is not based on Kolos' special 
C*    version for the hydrogen/antihydrogen interaction from 1975, 
C*    but on the H_2 version.
C*
C*    This program allows calculation for two-electron diatomic 
C*    systems. However, this version considers only \Sigma symmetry.
C*    Explicitly correlated two-electron basis functions expressed 
C*    in prolate-spheroidal coordinates are used in this program.
C*
C*    The basis functions used are not symmetry adapted to the 
C*    hydrogen/antihydrogen problem, i.e. they are not symmetric
C*    or antisymmetric with respect to a combined exchange of the 
C*    leptons and a reflection on a mirror plane perpendicular to
C*    the internuclear axis (so-called Q symmetry).
C*
C*    ++++  In contrast to HHSP01.F the output is minimised!  ++++
C*
C*    The program can be installed in DOUBLE or QUADRUPLE precision.
C*    Both has been tested by A. Saenz. The quadruple precision should
C*    be helpful in checking for numerical inaccuracies, but is 
C*    extremely slow. In addition, it has been used (and thus tested)
C*    only very few times. The quadruple-precision version has been 
C*    implemented only on a CRAY computer, while the double-precision
C*    version has been implemented on DEC risc machines (running 
C*    ULTRIX), IBM risc machines (running AIX), and a CRAY (running
C*    UNICOS). The quadruple version is activated by uncommenting 
C*    the lines following CQP and commenting out the lines following 
C*    CDP (where usually the number of lines to uncomment or to 
C*    comment out is given after CQP or CDP). HOWEVER, on e.g. a 
C*    CRAY computer the statements of the type DLOG(X)=LOG(X),  
C*    DEXP(X)=EXP(X), DSIN(X)=SIN(X) etc. should NOT be uncommented
C*    when using quadruple precision, although this is stated in the
C*    program!
C*
C*    This program has been written using the FORTRAN speciality 
C*    of IMPLICIT statements. To make things even worse, DIFFERENT 
C*    IMPLICITs are used IN DIFFERENT SUBROUTINES!!! So be careful 
C*    when doing modifications, since otherwise the constant pi 
C*    may have the value 3.0, as it happened earlier ...
C*
C*    Many modifications have been done by A. Saenz from 1992 on,
C*    and not all of them are explicitly stated. Especially, a 
C*    number of comments was added to the original puristic code
C*    (`no comments', even not to separate subroutines, logical 
C*    fragments etc., since comments vaste space...).
C*
C*    The most important modification was to increase the size of 
C*    the basis set which can be used in the program, since the 
C*    version inherited did only work properly for basis sets 
C*    in the range of about 10 to 30 basis functions. The increase 
C*    was done in steps, first to 100, then to 200 and finally to 
C*    400 basis functions. The problem of a lower limit (that was 
C*    already present in the inherited version) is still existent.
C*    Although in most applications this should not be a problem, 
C*    this is annoying when trying to perform small test runs. 
C*
C*    In contrast to the original program, an additional storage 
C*    optimisation was performed by first calculating the actually
C*    needed dimensions of the largest arrays using a preprogram, 
C*    and then to produce an executable prior to program execution. 
C*    In this way the run-time storage space is optimised.
C*
C*    P. Froelich has introduced the possibility of writing out 
C*    two matrices which are needed for performing complex-scaling
C*    calculations. There was however some confusion regarding the 
C*    correct matrices to be stored, which has been later resolved 
C*    by A. Saenz, who also introduced the possibility to extract 
C*    some other matrices which may be useful for different purposes.
C*
C*    Note, that the subroutines are simply ordered in alphabetical
C*    order and not by any program logics.
C*
C*    Since the code has gone through some hands, and since especially 
C*    A. Saenz did not have any direct contact with any of the 
C*    father(s) of the code when doing the modifications, there is 
C*    an obvious risk of errors.
C*
C*      ---> Therefore, do not take any of the results (and even 
C*           not the comments) for granted!!! 
C*
C*
C*------------------------------------------------------------------- 
C*
CQP      1
C     IMPLICIT REAL*16 (A-H,O-Z)
CDP      1
      IMPLICIT REAL*8 (A-H,O-Z)
C
C*    The idea of the program is to store nearly all of the arrays
C*    in one big blank common with the name F. Since this usually huge 
C*    one-dimensional array contains auxiliary functions for 
C*    calculating matrix elements as well as the different matrices
C*    of overlap, kinetic energy, etc., the dimensioning of that 
C*    array depends on a number of different parameters. In the 
C*    case of small basis sets and thus small matrices, the storage
C*    requirement is mainly determined by the number of auxiliary
C*    functions, while in the case of large basis sets the storage 
C*    requirements for the matrices becomes dominant.
C*
C*    The version inherited by A. Saenz contained unfortunately at 
C*    different places contradictory statements regarding the 
C*    dimensioning requirements. The most useful was the following:
C*  
C     MINIMAL SIZE OF THE BLANK COMMON IS DECIDED BY THE SUBROUTINES 
C     INTFI AND ASYFI WHERE IT IS EQUAL 20213. IT IS SET TO A NUMBER
C     DIFFERENT FROM 1 ALSO IN DIAG (601), DIAGDR (301) AND IN
C     EINT (101) BECAUSE OF SOME EQIVALENCING.
C     FOR SMALL BASIS SETS THE SIZE OF BLANK COMMON SHOULD BE
C     30113 + 3*NR*(NR+1)/2
C     FOR LARGER BASIS SETS, I.E. IF (NR*(NR+1)/2 + 600) .GT. 20213
C     THE SIZE OF THE BLANK COMMON SHOULD BE 2*NR*(NR+1)+600
C     IN EACH CASE THE SIZE OF THE BALNK COMMON MUST BE AT
C     LEAST ISMOLD1 + 3*NR*(NR+1)/2, WHERE
C     ISMOLD1 = 2*(LMX+2)*(PMX+1)**2 + 6*(KMX+1)*(LMX+1)
C*
C*    As stated above, the storage requirements have been linked 
C*    to the actual problem by evaluating the dimensioning parameters
C*    using a preprogram and recompiling this program for every 
C*    calculation.
C*
      INTEGER nblmax,nrmax
C*
C*    The parameter nrmax and nblmax are contained in the 
C*    file h2upar.for:
C*
      INCLUDE 'h2upar.for'
C*
      COMMON F(nblmax)
C*
      DIMENSION HS(nrmax*nrmax)
      DATA ICORE,NR/nblmax,nrmax/
C*
C*------------------------------------------------------------------
C*
C*    The real main program is contained in subroutine HHMOL which is 
C*    the core of the program. In that subroutine the different steps
C*    (integral calculation, matrix formation, matrix diagonalisation
C*     etc.) are initiated by calling the corresponding subroutines.
C*
      CALL HHMOL(HS,NR,ICORE)
C*
      STOP
      END
C*
C*======================================================
C*
      BLOCK DATA
C*
C*------------------------------------------------------
C*
CQP      3
C     IMPLICIT REAL*16(A-H),INTEGER(P)
C     REAL*16 S,V
C     REAL*16 R,SEPAT,IZ1,IZ2,X
CDP      3
      IMPLICIT REAL*8 (A-H),INTEGER(P)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
      PARAMETER ( nrdim2 = 2 * nrdim )
C*
      REAL*8 S,V
      REAL*8 R,SEPAT,IZ1,IZ2,X
      INTEGER OUT
C
      COMMON/AREAR/A1,A2,B1,B2
C*ASB  1  21.07.93
c      COMMON/AREAI/LMX,PMX,MAXPN,NOFTE,ITRIP,IPAR,I(300)
      COMMON/AREAI/LMX,PMX,MAXPN,NOFTE,ITRIP,IPAR,I(nrdim2)
C*
      COMMON/AREA1R/A(36)
      COMMON/AREA1I/N
      COMMON/AREA2/J(30)
      COMMON/AREA3/JA(20)
      COMMON/AREA4R/C(110)
      COMMON/AREA4I/K1
      COMMON/AREA5R/D(7)
      COMMON/AREA5I/M(5)
      COMMON/AREA6R/X(2)
      COMMON/AREA6I/LTE(10),KMX
      COMMON/AREA7R/R,IZ1,IZ2,SEPAT
      COMMON/AREA7I/IPNCH,IPRFUN
      COMMON/AREA8R/CONVE,E
      COMMON/AREA8I/ITER,NIN,IASY,NOFT1,ITEST,NOFT2,NOFT3,
     *ISTP,ISTP1,IDGPR,IO(8)
C
      DATA A1,A2,B1,B2/4*0./
C*ASB 1  21.07.93
c      DATA LMX,PMX,MAXPN,NOFTE,ITRIP,IPAR,I/306*0/
      DATA LMX,PMX,MAXPN,NOFTE,ITRIP,IPAR/6*0/
      DATA I/nrdim2*0/
C*
      DATA A/36*0./,N/0/
      DATA J/30*0/
      DATA JA/20*0/
      DATA C/110*0./
      DATA D/7*0./,M/5*0/
      DATA X/0.,0./,LTE,KMX/11*0/
      DATA R/0./,IZ1,IZ2/2*0./,SEPAT/0./,IPNCH,IPRFUN/2*0/
      DATA CONVE,E/2*0./,ITER,NIN,IASY,NOFT1,ITEST,NOFT2/6*0/,
     $NOFT3,ISTP,ISTP1,IDGPR/4*0/
      END
C*
C*=================================================================
C*
      SUBROUTINE CCI
C*
C*-----------------------------------------------------------------
C*
CQP      3
C     IMPLICIT INTEGER(R,S),REAL*16(A-C)
C     GENERIC
C     REAL*16 IZ1,IZ2
CDP      2
      IMPLICIT INTEGER(R,S),REAL*8(A-C)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      REAL*8 IZ1,IZ2
      COMMON /AREA3/S(3),MPN,MMN,RPR,RMR,SPS,SMS
      COMMON /AREA4R/C(24),APA,AMA,BPB,BMB
      COMMON /AREA5R/A(2)
      COMMON /AREA5I/IT(4),J2
      COMMON /AREA7R/AR,IZ1,IZ2
C*
C*-----------------------------------------------------------------
C*
      C(6) = 2.0D+00
      C(7) = 2.0D+00
C*
C*
C*       Set C(8) = - 4 * ( Z_1 + Z_2 )  and
C*           C(9) = + 4 * ( Z_1 - Z_2 )
C*       ( where Z_1 (=IZ1) is the nuclear charge of nucleus 1 and
C*               Z_2 (=IZ2) is the nuclear charge of nucleus 2. )
C* 
      C(8) = C(23)
      C(9) = C(24)
C*
C*          If \mu_j .NEQ. \mu_i  go to label 6:
      IF (MMN.NE.0) GO TO 6
C*
C*
C*               ***  If \mu_j = \mu_i ***
C* 
C*          then a number of terms simplify, i.e. 
C*          C(1) = \mu_j + \mu_i   and
C*          C(2) = C(3) = C(4) = C(5) = C(10) = C(11) = C(12) = C(13) = 0:
C*
      C(1) = MPN
      DO 5 I=2,5
        C(I+8) = 0.0D+00
        C(I) = 0.0D+00
    5 CONTINUE
      GO TO 7
C*
C*
C*               ***  If \mu_1 .NEQ. \mu_2  ***
C*
C*          set:   
C*          C(1) = (\mu_j-\mu_i) * [ (\mu_j-\mu_i) + (r_j-r_i) 
C*                       + (s_j-s_i) ] + (\mu_j + \mu_i)  ,
C*          C(3) = (\mu_j-\mu_i) * (\alpha - \bar{\alpha}) ,
C*          C(5) = (\mu_j-\mu_i) * ( (s_j-s_i) - (r_j-r_i) )
C*          C(10) =   2 * C(3)
C*          C(12) = - 2 * (s_j-s_i) * (\mu_j-\mu_i)
C*          C(13) =   2 * (r_j-r_i) * (\mu_j-\mu_i)   :  
C*
C*          NOTE: s_i or s_j may also be \bar{s}_i or \bar{s}_j and
C*                r_i or r_j may also be \bar{r}_i or \bar{r}_j, 
C*                since those quantities are sometimes exchanged in 
C*                subroutine MX(ni) (after label 8 of that routine).
C*
    6 CONTINUE
      C(1) = MMN * ( MMN + RMR + SMS ) + MPN
      C(3) = MMN * AMA
      C(5) = MMN * ( SMS - RMR )
      C(10) = 2 * C(3)
      C(13) = 2 * RMR * MMN
      C(12) = - 2 * SMS * MMN
C*
C*
C*               *** In all cases (\mu_j .EQ. or .NEQ. \mu_i) ***
C*
C*         C(17) to C(22) are defined by:
C*
    7 CONTINUE
      C(17) = SMS**2 - SPS
      C(19) = AMA**2
      C(21) = RPR - RMR**2
      C(20) = - 2 * AMA * RMR
      C(18) = C(3) - C(20) - 2 * APA
      C(22) = - C(5) - C(21) - C(17) + 2 * ( RPR - SPS ) - C(19)
      GO TO 9
C*
      ENTRY CCHI
C*    ==========
C*
      C(22) = C(22) + C(15)
C*
C*
C*       Exchange the values of BPB and BMB:
C* 
      C(15) = BPB
      BPB = BMB
      BMB = C(15)
C*
      IF (C(6).GT.0) GO TO 9
C*
      DO 8 I=1,22
        C(I) = - C(I)
    8 CONTINUE
C*
    9 CONTINUE
C*
      C(15) = - BMB**2
      C(22) = C(22) - C(15)
      C(16) = - 2 * BMB * SMS
      C(2) = - MMN * BMB
      C(11) = - 2 * C(2)
      C(14) = - C(2) - C(16) - 2 * BPB
C*
      GO TO(12,10,12),J2
C*
   10 CONTINUE
C*
      DO 11 I=1,22
        C(I) = - C(I)
   11 CONTINUE
C*
   12 CONTINUE
C*
      RETURN
      END
C*
C*================================================================
C*
      SUBROUTINE CNTR(MAT,VEC,VE,N,J)
C*
C*----------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(M,V)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(M,V)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      DIMENSION MAT(J),VEC(N),VE(N)
C*
C*----------------------------------------------------------------
C*
      DO 3 I=1,N
    3 VE(I)=0
      K=N
    4 I=J
      I1=N
    5 VE(K)=VE(K)+MAT(I)*VEC(I1)
      I1=I1-1
      IF(K-I1)6,6,7
    6 I=I-I1
      GO TO 5
    7 I=I-1
      IF(I1)8,8,5
    8 K=K-1
      J=J-1
      IF(K)9,9,4
    9 GO TO 13
C*
      ENTRY CNTRA(MAT,VEC,VE,N0,J,N)
C*    ==============================
      DO 12 I=1,N
   12 VE(I)=0.
      I=J-(N-N0)*N
      K=N0
      DO 11 I3=1,N0
      I=I-N+N0
      I1=N0
      DO 10 I2=1,N0
      VE(K)=VE(K)+MAT(I)*VEC(I1)
      I=I-1
   10 I1=I1-1
   11 K=K-1
C*
      ENTRY CNTRC(MAT,VEC,VE,N,J,NOF)
C*    ===============================
      DO 21 I=1,N
   21 VE(I)=0.0
      K=N
   22 I=J
      I1=N
   23 VE(K)=VE(K)+MAT(I)*VEC(I1)
      I1=I1-1
      I=I-1
      IF(I1)24,24,23
   24 K=K-1
      J=J-NOF
      IF(K)13,13,22
   13 RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE DIAG(N)
C*
C*------------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(A-F,R-X)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-F,R-X)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
C*ASB 1  21.07.93
c      COMMON F(601)
      COMMON F((2*nrdim)+1)
C*
      COMMON /AREAR/A(4)
      COMMON /AREAI/L(6),IPOWR(50)
      COMMON /AREA2/IT(17),ISM,IVM,IKM
      COMMON /AREA7R/R
      COMMON /AREA8R/ CONVE,E
      COMMON /AREA8I/ ITER,INP(9)
      COMMON /AREA11/X
      DIMENSION AMAT(1),C(1),VC(1),EI(11)
C*ASB 1  21.07.93 (2 changes!!!)
c      EQUIVALENCE (C(1),F(1)),(VC(1),F(301)),(AMAT(1),F(601))
      EQUIVALENCE (C(1),F(1)),(VC(1),F(nrdim+1)),
     +   (AMAT(1),F((2*nrdim)+1))
C*
C*-----------------------------------------------------------------
C*
CQP      2
C     DABS(X)=ABS(X)
C     DSQRT(X)=SQRT(X)
C
      IDGP=INP(9)
      EI(1)=E
      N1=(N*(N+1))/2
      III=0
C
  100 FORMAT(/1X,'SUBROUTINE DIAG'/)
  101 FORMAT(1X,'ENERGY E=',E14.7)
  200 FORMAT(/7E14.7)
  201 FORMAT(//)
C
      DO 8 K=1,10
    4 DO 5 I=1,N1
      J=IKM+I-1
      J1=ISM+I-1
    5 AMAT(I)=F(J)-E*F(J1)
      J=ISM+N1-1
      CALL CNTR(F,C,VC,N,J)
      DO 51 I=1,N
   51 C(I)=VC(I)
      CALL SEQU(AMAT,C,N,N1)
      J=ISM+N1-1
      CALL CNTR(F,C,VC,N,J)
      CALL SCPR(C,VC,N)
      X=DABS(X)
      X=1./DSQRT(X)
      DO 6 I=1,N
    6 C(I)=X*C(I)
      IF(ITER.LE.0) GO TO 7
      ITER=ITER-1
      GO TO 4
    7 J=IKM+N1-1
C*
C     WRITE(6,1000)
C1000 FORMAT(10X,' SUBROUTINE DIAG TABLICE F C VC X  ')
C     WRITE(6,1001) J,N,IKM,N1
C1001 FORMAT(1X,10I6)
C      WRITE(6,200) (F(KKK),KKK=IKM,J)
C     WRITE(6,201)
C*
      CALL CNTR(F,C,VC,N,J)
C*
C      WRITE(6,200) (F(KKK),KKK=IKM,J)
C     WRITE(6,201)
C     WRITE(6,200)(C(KKK),KKK=1,N)
C     WRITE(6,201)
C     WRITE(6,200)(VC(KKK),KKK=1,N)
C     WRITE(6,201)
C     WRITE(6,1001) J,N,IKM,N1
C     WRITE(6,201)
C*
      CALL SCPR(C,VC,N)
C*
C     WRITE(6,200)(C(KKK),KKK=1,N)
C     WRITE(6,201)
C     WRITE(6,200)(VC(KKK),KKK=1,N)
C     WRITE(6,201)
C     WRITE(6,1001) J,N,IKM,N1
C     WRITE(6,1002)X
C1002 FORMAT(10X,' X  W DIAG',E14.7)
C*
      E=X
      EI(K+1)=X
      X=1.-EI(K)/X
      X=DABS(X)
      IF(X.LE.CONVE) GO TO 9
    8 CONTINUE
      WRITE(6,41)
      K=11
      GO TO 90
    9 IF(IDGP.EQ.0)GO TO 91
   90 WRITE(6,42)(EI(I),I=1,K)
   91 X=2./R
      X=X**3
      DO 10 I=1,N
   10 C(I)=X*C(I)
   41 FORMAT(1H0,16HPOOR CONVERGENCE)
   42 FORMAT(1H ,5D20.10)
      RETURN
      END
C*
C*================================================================
C*
      SUBROUTINE DIAGDR(H,NR,N1,ITEST)
C*
C*----------------------------------------------------------------
C*
CQP      3
C     IMPLICIT REAL*16(A-H,O-X)
C     GENERIC
C     REAL*16 YSEP,IZ1,IZ2
CDP      2
      IMPLICIT REAL*8(A-H,O-X)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      REAL*8 YSEP,IZ1,IZ2
C*ASB  1  21.07.93
c      COMMON F(301)
      COMMON F(nrdim+1)
C*
      COMMON /AREAR/A(4)
C*ASB  1  21.07.93  
c      COMMON /AREAI/L(6),IPOWR(300)
      COMMON /AREAI/L(6),IPOWR(nrdim)
C*
      COMMON /AREA2/IT(17),ISM,IVM,IKM
      COMMON /AREA6I/LTE(10),KMX
      COMMON /AREA7R/R,IZ1,IZ2,YSEP
      COMMON /AREA7I/IPNCH,IPRFUN
      COMMON /AREA8R/CONVE,E
      COMMON /AREA8I/ITER,NIN,IP(8)
      COMMON /AREA11/X
      COMMON /S/IO,IO2,NS,NEGV,NDEL,IHESHE,INPVC
      COMMON /S2/CINP(10)
      COMMON /FPS/MXRN,NAP,IAP
      common /teta/teta
      DIMENSION H(NR,1),ENER(nrdim),V(nrdim)
      DIMENSION C(1),VC(nrdim)
C*ASB 1  21.07.93
c      EQUIVALENCE (C(1),F(1)),(VC(1),F(301))
      EQUIVALENCE (C(1),F(1)),(VC(1),F(nrdim+1))
C*
      EQUIVALENCE (L(4),IN)
C*
      CHARACTER*256 filnam
C*
C*-----------------------------------------------------------------
C*
CQP      2
C     DABS(X)=ABS(X)
C     DSQRT(X)=SQRT(X)
C*
      REWIND 9
      JJ=1
      LC=0
      IF(NIN.LT.0)GO TO 4
C*
      J=IP(6)/100
      I=IP(6)-100*J
      STEP=-J*(0.1**I)
      J=IP(7)/100
      I=IP(7)-100*J
      STEP1=-J*(0.1**I)
      IF(ITEST.EQ.1) WRITE(6,98) IP(6),IP(7),STEP,STEP1
      I1=(L(4)*(L(4)+1))/2
      X=1./R
      X1=-2*X*X
      X=2*X
      I2=(N1*(N1-1))/2
      J=IKM+I2-1
      J1=IVM+I2-1
      I2=I2+1
      DO 5 I=I2,I1
      J=J+1
      J1=J1+1
    5 F(J)=X1*F(J)+X*F(J1)
    4 NIN=IABS(NIN)
      IF(NIN.EQ.0.OR.NIN.GT.IN) NIN=IN
      IF(E.EQ.0.0) E=F(IKM)/F(ISM)
      EN=E
      IF(ITER.NE.0)EN=1.
      J1=L(4)
      J=NIN
      IV=ISM-J1-1
      N=J
      C(1)=1.
      DO 6 I=2,IN
    6 C(I)=0.0
      IF(INPVC.EQ.0)GO TO 64
C*
      I=ISM+I1-1
      CALL CNTR(F,CINP,VC,IN,I)
      CALL SCPR(CINP,VC,IN)
      X=DABS(X)
      X=1./DSQRT(X)
C*
      DO 61 KK=1,IN
   61 C(KK)=X*CINP(KK)
C*
      I=IKM+I1-1
C*
      CALL CNTR(F,C,VC,IN,I)
      CALL SCPR(C,VC,IN)
C*
      E=X
      GO TO 3001
C*
   64 KK=(N*(N+1))/2-1
       KIVM=IVM+KK
       KIKM=IKM+KK
       KISM=ISM+KK
C*
      IPR5=0
      IF(IPNCH.EQ.0)GO TO 119
C
C      ***** Save results and parameters *****
C
CKSZ92
      WRITE(10,90)IN,R,(A(KJ),KJ=1,4),L(1),L(2),KMX,L(3),teta
      WRITE(10,79)(IPOWR(KJ),KJ=1,IN)
  119 CONTINUE
       ILE=0
    7 CONTINUE
      IF(IO.EQ.0) GO TO 3000
C
C      ***** Compute single eigenvalue (if chosen) *****
C*
      CALL DIAG(N)
      GO TO 3001
C
C       ***** Complete diagonalisation *****
C*
 3000 CONTINUE
      NNH=IKM-1
      NNS=ISM-1
 2050 CALL DIAGON(H,NR,N,V,ENER,DIAGG,ISM,IVM,IKM)
C*
 2056 E=ENER(ITER)
      ILE=ILE+1
      WRITE(9) ILE,N,(ENER(IIE),IIE=1,N)
C*
      DO 2002 IIE=1,N
 2002 C(IIE)=H(IIE,ITER)
C
C      ***** TEST OF TERMS *****
C*
 3001 IF(ITEST.EQ.0)GO TO 10
      DE=E-EN
      J=J+1
      EN=E
      WRITE(6,99) J,J1,DE,STEP
   99 FORMAT(1X,2I5,2D10.2)
      IF(DE.LE.STEP)GO TO 8
      IF(DE.LT.STEP1.AND.J.GT.J1)GO TO 8
      C(N)=0.
      WRITE(6,41)IPOWR(N),E,DE
      E=E-DE
      EN=E
      IF(N.GE.IN)GO TO 72
      CALL RELOC(ISM,IV,IN,N)
      CALL RELOC(IVM,IV,IN,N)
      CALL RELOC(IKM,IV,IN,N)
      J2=IPOWR(N)
      J3=IN-1
      DO 71 I=N,J3
   71 IPOWR(I)=IPOWR(I+1)
      IPOWR(IN)=J2
   72 IF(DE.LE.STEP1)GO TO 9
      IN=IN-1
      GO TO 9
    8 WRITE(6,40)IPOWR(N),E,DE
      N=N+1
    9 IF(N.LE.IN)GO TO 7
      IF(IN.LT.IP(4).AND.IP(5).GT.0)GO TO 15
      N=IN
C
C      ***** Print output *****
C
   10 K=IVM+(N*(N+1))/2-1
C*
 201  FORMAT(/'MATRIX IN DIAGDR BEFORE DIAGONALIZATION')
 200  FORMAT(6E14.7)
C*
      CALL CNTR(F,C,VC,N,K)
C*
      CALL SCPR(C,VC,N)
C*
      X=X*(R/2)**5
      IF(IO.EQ.0)X=X*(2./R)**6
      VC(2)=219474.62D00*(YSEP-E)
      VC(1)=(X-E*2.)/R
      IF(IPR5.EQ.1)GO TO 105
      IF(MXRN.GT.0)GO TO 100
      WRITE(6,42) 
  100 IF(IZ1.EQ.DABS(IZ2))GO TO 101
      WRITE(6,75)
      GO TO 102
C*
  101 IF(IP(1).EQ.0)GO TO 103
      WRITE(6,80) 
      GO TO 102
C*
  103 WRITE(6,43)
  102 WRITE(6,48)IZ1,IZ2,L(3),L(2),L(1)
      WRITE(6,50)R,(A(I),I=1,4),teta
      VC(1)=(X-E*2.)/R
      WRITE(6,52)E,VC(2),X,VC(1),YSEP
C*ASB   7  17.09.1998
      CALL getenv('VIRIAL',filnam)
      OPEN (96, FILE=filnam, FORM='FORMATTED', 
     +      ACCESS='SEQUENTIAL')
      WRITE(96,*) E
      WRITE(96,*) X
      WRITE(96,*) VC(1)
      WRITE(96,*) R
      CLOSE(96)
C***      
      IF(IPRFUN.EQ.1) GO TO 106
      IF(MXRN.GT.0)GO TO 106
CC      WRITE(6,76)
      NP=N+1
      NP=NP/2.
      ID=1
      IC=1
   11 IP1=IPOWR(IC)/10000
      LP=IPOWR(IC)-10000*IP1
      IP2=LP/1000
      LP=LP-1000*IP2
      IP3=LP/100
      LP=LP-100*IP3
      IP4=LP/10
      IP5=LP-10*IP4
      IIC=IC+NP
      ID=ID+1
      IF(ID.LE.N)GO TO 1100
CC      WRITE(6,44)IC,IP1,IP2,IP3,IP4,IP5,C(IC)
      GO TO 13
 1100 JP1=IPOWR(IIC)/10000
      LP=IPOWR(IIC)-10000*JP1
      JP2=LP/1000
      LP=LP-1000*JP2
      JP3=LP/100
      LP=LP-100*JP3
      JP4=LP/10
      JP5=LP-10*JP4
CC      WRITE(6,45)IC,IP1,IP2,IP3,IP4,IP5,C(IC),IIC,JP1,JP2,JP3,JP4,JP5,C(
CC     CIIC)
      ID=ID+1
      IC=IC+1
      IF(ID.LE.N)GO TO 11
   13 GO TO 106
  105 IP1=IPOWR(N)/10000
      LP=IPOWR(N)-10000*IP1
      IP2=LP/1000
      LP=LP-1000*IP2
      IP3=LP/100
      LP=LP-100*IP3
      IP4=LP/10
      IP5=LP-10*IP4
      DELD=DELD-VC(2)
CC      WRITE(6,74)N,IP1,IP2,IP3,IP4,IP5,E,VC(2),DELD,VC(1)
  106 CONTINUE
      IF(IPNCH.EQ.0) GO TO 117
C
C     ***   PUNCH OUTPUT TAPE   ***
C                                  *** IMPORTANT.
C                                      OLD OUTPUT TAPES MUST HAVE
C                                      A LINE WITH NUMBER OF FUNC
C                                      ADDED BETWEEN POWERS AND E
C                                                    AUGUST 84 KS ***
C     *** TO AVOID DOUBLE PUNCH WHEN N=IN ***
C
      IF(LC.NE.0) GO TO 108
C
C     *** THIS LINE INTRODUCES INCOMPATIBILITY WITH OLD TAPE FORMAT
C
      WRITE(10,94) N
      IF(IO.NE.0.OR.IO2.NE.0)GO TO 110
      WRITE(10,91)(ENER(J),J=1,NEGV)
C*
C*    Instead of writing the eigenvectors into a FORMATTED file,
C*    they are now written (together with the eigenvalues) to an 
C*    UNFORMATTED (and thus smaller) file:
C***
      IF (NEGV.NE.0) THEN
        CALL getenv('EIGVEC',filnam)
        OPEN ( UNIT=35 , FILE =filnam, STATUS='UNKNOWN',
     +         FORM='UNFORMATTED', ACCESS='SEQUENTIAL' )
        WRITE(35) N,NEGV 
        WRITE(35) (ENER(J),J=1,NEGV)
        WRITE(35) ((H(J,JV),J=1,N),JV=1,NEGV)
        CLOSE(35)
      ENDIF
C*
c      WRITE(10,78)((H(J,JV),J=1,N),JV=1,NEGV)
C***
      GO TO 121
  110 WRITE(10,92)E
      WRITE(10,78)(C(J),J=1,N)
  121 CONTINUE
      INEXT=1
      IF(N.GE.IN) INEXT=0
      WRITE(10,93) INEXT
  117 N=N+NDEL
      IF(N.GT.IN)GO TO 108
      IPR5=1
      DELD=VC(2)
      IF(JJ.EQ.0)GO TO 107
      WRITE(6,69)
      WRITE(6,70)
      JJ=0
C*
C     LOOP OVER ENERGY CASCADE 
C*
  107 GO TO 7
  108 IF(LC.NE.0)GO TO 15
  109 LC=1
      N=IN
      IF(IPR5.EQ.1)GO TO 100
C
   31 FORMAT (2F10.5,2I5)
   40 FORMAT(1H ,4X,I5,10X,2F15.10)
   41 FORMAT(1H ,50X,I5,10X,2F15.10)
   42 FORMAT(1H1)
   43 FORMAT('0',40X,'H anti{H} clamped nuclei computation',18X,'
     CMunich, 1998'//)
   44 FORMAT(4X,I3,I4,2(I2,I1),D19.9)
   45 FORMAT(4X,I3,I4,2(I2,I1),D19.9,6X,I3,I4,2(I2,I1),D19.9)
   46 FORMAT(1H ,15X,I5,10X,D18.8)
   48 FORMAT(25X,'Z1 =',F7.4,3X,'Z2 =',F7.4,3X,
     C           'Int.P. =',I4,3X,'Pmax =',I3,3X,'Lmax =',I3/)
   50 FORMAT(28X,'R =',F8.4,3X,'a1 =',F12.8,3X,'a2 =',F12.8/
     C       40X,'b1 =',F12.8,3X,'b2 =',F12.8/
     c       40x,' theta = ',f12.6/)
CQP      2
C  52 FORMAT(18X,'E =',F30.25,3X,'D =',F30.20,' 1/CM   V =',F30.25/
C    C       18X,'dE/dR =',F30.25,3X,'Sep.At. =',F30.25//)
CDP      2
   52 FORMAT(18X,'E =',F20.15,3X,'D =',F20.10,' 1/CM   V =',F20.15/
     C       18X,'dE/dR =',F20.15,3X,'Sep.At. =',F20.15//)
   69 FORMAT(1H0,64X,13HTEST OF TERMS)
   70 FORMAT(1H0,35X,3H  N,3X,8HFUNCTION,6X,11HENERGY (AU),7X,8HD (1/CM)
     C,7X,7HDELTA D,15X,5HDE/DR)
   74 FORMAT(1H ,36X,I3,3X,I1,2(I2,I1),5X,F13.10,3X,F11.4,4X,F11.4,5X,F1
     C5.10)
   75 FORMAT('0',40X,'Exotic molecule clamped nuclei computation',18X,'
     CMunich, 1998'//)
   76 FORMAT(1H0,5X,1HN,10H  Function,8X,6HVector,13X,1HN,10H  Function,
     C       7X,6HVector/)
   77 FORMAT(2(I2),F8.4,I4,4(F7.3),F13.9,F10.2,F13.9)
   78 FORMAT(3D25.14)
   79 FORMAT(10I6)
   80 FORMAT(' This is NOT a Pi-state calculation...')
CKSZ92
   90 FORMAT(I3,F10.4,4F10.6,4I3,f6.3)
   91 FORMAT(4D20.9)
   92 FORMAT(D20.9)
   93 FORMAT(10I1)
   94 FORMAT(10I5)
   98 FORMAT('1SELECTION PARM. IP(6),IP(7),STEP,STEP1 ',2I3,2D10.2)
C
   15 IF(IO.NE.0.OR.ITEST.NE.0)GO TO 3002
CC      WRITE(6,2019)
      N1=1
      N10=5
C      *********** PRINT ALL EIGENVALUES ****************
C      ***** COMPUTED FOR SEVERAL EXPANSION LENGTHS *****
      NN1=N1
      REWIND 9
 2008 READ(9,END=2010) ILE,NNE,(H(ILE,IIE),IIE=1,NNE)
      GO TO 2008
 2010 CONTINUE
 2023 IF(NN1.GT.IN) GO TO 2024
      IF(N10.GT.IN) N10=IN
CC      WRITE(6,2020) (MJ,MJ=NN1,N10)
      III=0
      DO 2021 II=NIN,IN,NDEL
      III=III+1
      N11=N10
C     IF(NN1.GT.II) WRITE(6,2022) II
      IF(NN1.LE.II) GO TO 2025
      GO TO 2021
 2025 IF(N10.GT.II) N11=II
CC      WRITE(6,2022) II,(H(III,JJ),JJ=NN1,N11)
 2021 CONTINUE
      NN1=N10+1
      N10=N10+5
      GO TO 2023
 2024 CONTINUE
 2019 FORMAT(/,11X,' Number of functions / eigenvalues',/)
 2020 FORMAT(15X,I4,4(17X,I4))
 2022 FORMAT(1X,I4,3X,5(F20.15))
 3002 RETURN
      END
C*
C*===================================================================
C*
      SUBROUTINE DIAGON(HS,NR,N,E,V,DIA,ISM,IVM,IKM)
C*
C*-------------------------------------------------------------------
C*
C     SPECIAL VERSION FOR H2 PROGRAM. 9/6/85 KSZ.
C*
CQP      2
C     IMPLICIT REAL*16(A-H,O-Z)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-H,O-Z)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar.for'
      INCLUDE 'h2upar2.for'
C*  
c      COMMON F(1)
C*ASB 1   02.09.1998
      COMMON F(nblmax)
C*
      DIMENSION HS(NR,N)
      DIMENSION E(1),V(1)
C*
      CHARACTER*256 filnam
C*
C*------------------------------------------------------------------
C*
      DIA=16.0D0**(-13)
  602 FORMAT('  EIGENVALUES',3D25.15/(12X,3D25.15))
C                           S --> HS
      ISC=ISM-1
      DO 10 I=1,N
      DO 10 J=1,I
      ISC=ISC+1
C***
c          IF ((isc.EQ.ism).OR.(isc.EQ.ism+1).OR.(isc.EQ.ism+2)) THEN
c            WRITE(*,*) ' Checkpt.5a: ISC, F(ISC) = ',isc,F(isc)
c          ENDIF
C***
   10 HS(I,J)=F(ISC)
C*
C*    Store the overlap matrix:
C*
C*ASB 4  20.04.1994
      CALL getenv('OVLAP',filnam)
      OPEN (32, FILE=filnam, STATUS='UNKNOWN',
     +      FORM='UNFORMATTED', POSITION='APPEND')
      WRITE(32) ((HS(IPIO,JPIO),JPIO=1,IPIO),IPIO=1,N)
      CLOSE(32)
C*
      CALL SPOL(HS,NR,N)
      CALL REV(HS,NR,N)
C*
C*                        Q --> F(2*nrdim+1,...)
C*
C*ASB 1  21.07.93
c      I1=600
      I1=2*nrdim
C*
      DO 20 I=1,N
      DO 20 J=1,I
      I1=I1+1
   20 F(I1)=HS(I,J)
C*
C*
C*    Store the Q matrix that transforms the basis into an 
C*    orthogonal basis:
C*
      CALL getenv('QMAT',filnam)
      OPEN (32, FILE=filnam, STATUS='UNKNOWN',
     +          FORM='UNFORMATTED', POSITION='APPEND')
C*
      WRITE(32) ((HS(IPIO,JPIO),JPIO=1,IPIO),IPIO=1,N)
C*
C*
C*    Store the matrix of potential energy (V)
C*
C*                       V --> HS
C*
      ivc = ivm - 1
      DO 28 i=1,n
        DO 26 j=1,i
          ivc = ivc + 1
C***
c          IF ((ivc.EQ.ivm).OR.(ivc.EQ.ivm+1).OR.(ivc.EQ.ivm+2)) THEN
c            WRITE(*,*) ' Checkpt.5b: IVC, F(IVC) = ',ivc,F(ivc)
c          ENDIF
C***
          hs(i,j) = f(ivc)
   26   CONTINUE
   28 CONTINUE
C*
C*
C*    Now the UNTRANSFORMED potential energy matrix V is stored:
C*
      CALL getenv('VMAT',filnam)
      OPEN (UNIT=38, FILE=filnam, STATUS='UNKNOWN',
     +      FORM='UNFORMATTED', POSITION='APPEND')
      WRITE(38) ((hs(i,j),j=1,i),i=1,n)
      CLOSE(38)
C*
C*    Transform V by multiplication with the Q matrix and its transpose:
C*
c      CALL QFQ(HS,F(601),NR,N,E)
      CALL QFQ(HS,F((2*nrdim)+1),NR,N,E)
C*
C*    Store the TRANSFORMED potential energy matrix V:
C*
C      OPEN (UNIT=38, FILE='qvq.out', STATUS='UNKNOWN',
C     +      FORM='UNFORMATTED', POSITION='APPEND')
C      WRITE(38) ((hs(i,j),j=1,i),i=1,n)
C      CLOSE(38)
C
C*      
C*                       H --> HS
C*
      IHC=IKM-1
      DO 30 I=1,N
      DO 30 J=1,I
      IHC=IHC+1
C***
c          IF ((ihc.EQ.ikm).OR.(ihc.EQ.ikm+1).OR.(ihc.EQ.ikm+2)) THEN
c            WRITE(*,*) ' Checkpt.5c: IHC, F(IHC) = ',ihc,F(ihc)
c          ENDIF
C***
   30 HS(I,J)=F(IHC)
C*
C*ASB 1  21.07.93
c      CALL QFQ(HS,F(601),NR,N,E)
      CALL QFQ(HS,F((2*nrdim)+1),NR,N,E)
C*
C*
C*    Store the TRANSFORMED Hamiltonian matrix:
C*
      CALL getenv('QHQMAT',filnam)
      OPEN (UNIT=33, FILE=filnam, STATUS='UNKNOWN',
     +      FORM='UNFORMATTED', POSITION='APPEND')
      WRITE(33) ((HS(IPIO,JPIO),JPIO=1,IPIO),IPIO=1,N)
C*
      CALL TRED2(NR,N,HS,V,E,HS)
      CALL TQL2(NR,N,V,E,HS,IERR)
C*ASB 1  21.07.93
c      CALL TRANS(HS,F(601),NR,N)
      CALL TRANS(HS,F((2*nrdim)+1),NR,N)
C*
      RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE EINT(BET,LMX,PMX,IA)
C*
C*------------------------------------------------------------------
C*
CQP      2
C     IMPLICIT INTEGER(P),REAL*16(A-F)
C     GENERIC
CDP      1
      IMPLICIT INTEGER(P),REAL*8(A-F)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON F(nrdim+1)
      DIMENSION A(1),C(2)
      EQUIVALENCE (F(1),A(1)),(F(nrdim+1),C(1))
C*
C*------------------------------------------------------------------
C*
CQP      1
C     DSQRT(X)=SQRT(X)
C*
      KMX1=PMX+1
      D1=1.
      I1=PMX+LMX
      I2=I1+LMX
      DO 6 I=1,I1
      A(I)=2./D1
    6 D1=D1+2.
      C(1)=1.
      D1=1.
      DO 7 I=1,I2
C***CHECK:
      IF (DABS(C(I)).LT.1.0D-300) THEN
        C(I+1) = 0.0D+00
      ELSE
        C(I+1)=C(I)*BET/D1
      ENDIF
C***
    7 D1=D1+1.
      L=0
      L1=-1
      D1=0.5
    8 D2=DSQRT(D1)
      L1=-L1
      N3=-1
      DO 13 K=1,KMX1
      N3=-N3
      N1=L+2-K
      N2=1
      B=0.
      IF(N1.GT.0) GO TO 10
      N1=1
      N2=(K-L+2)/2
      IF(L1.NE.N3)N1=2
   10 DO 11 N=N2,I1
      AINCR=C(N1)*A(N)
      B1=B+AINCR
      IF(B1.EQ.B) GO TO 12
      B=B1
   11 N1=N1+2
C*
CQP      1
C     IF(AINCR/B1.GT.1.E-26) WRITE(6,900) L,K,IA,BET,AINCR,B 
CDP      1
      IF(AINCR/B1.GT.1.E-13) WRITE(6,900) L,K,IA,BET,AINCR,B
C*
  900 FORMAT(' POOR CONV IN EINT: L,K,IA,BET,INCR,B ',3I5,F10.5,
     $       2E20.10)
   12 F(IA)=B*D2
   13 IA=IA+1
      D1=D1+1.
      L=L+1
      IF(L.GT.LMX) GO TO 15
      B=L
      B1=2.*B+1.
      DO 14 I=1,I1
      A(I)=A(I)*B/B1
      B=B+2.
   14 B1=B1+2.
      GO TO 8
   15 RETURN
      END
C*
C*================================================================
C*
      SUBROUTINE FIEV(ALF,DALF,IPMX,IF)
C*
C*----------------------------------------------------------------
C*
C*    In the case that no explicit 1/r_{12} term is present, 
C*    the so-called \Phi_{n,\bar{n}}^l integral introduced by 
C*    Ruedenberg (J. Chem. Phys. vol.19, p.1459(1951)) reduces 
C*    to a simple product of two simple one-dimensional integrals 
C*    that can be solved analytically by a stable recursion 
C*    formula. This subroutine calculates those simple integrals
C*    and multiplies them in the end in order to obtain \Phi.
C*
C*---------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(A-F)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-F)
      COMMON F(1)
      DIMENSION A(20),B(20)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
C*----------------------------------------------------------------
C*
CQP      1
C     DEXP(X)=EXP(X)
C
      L1=IPMX+1
      A(1)=DEXP(-ALF)/ALF
      DO 2 L=2,L1
    2 A(L)=A(L-1)*(L-1)/ALF+A(1)
      DO 3 L=1,L1
    3 B(L)=A(L)
C*
C*    If DALF=0.0 then only \Phi(\alpha,\alpha) has to be calculated,
C*    i.e. both arguments are identical. This occurs either when the
C*    user chose \alpha=\bar{\alpha}, or when the mixed term 
C*    \Phi(\alpha+\bar{\alpha},\alpha+\bar{\alpha_2}) has to be 
C*    calculated.
C*    (The sense of all that is to prevent the recursion to be done 
C*     twice for the same argument...)
C*
      IF(DALF.EQ.0.0)GO TO 5
C*
      ALF=ALF+2*DALF
      A(1)=DEXP(-ALF)/ALF
      DO 4 L=2,L1
    4 A(L)=A(L-1)*(L-1)/ALF+A(1)
C*
    5 J=IF
      DO 7 L=1,L1
      DO 6 I=1,L1
      F(J)=A(I)*B(L)
    6 J=J+1
    7 CONTINUE
C*
      RETURN
      END
C*
C*================================================================
C*
      SUBROUTINE FLM(IF)
C*
C*----------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(A-F),INTEGER(P)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-F),INTEGER(P)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON F(1)
      COMMON/AREAR/A(4)
      COMMON/AREAI/LMX,PMX
      COMMON/AREA5R/AA(2)
      COMMON/AREA5I/M
      COMMON/AREA6I/J(2),IDEL,IDEL1
C*
C*---------------------------------------------------------------
C*
CQP      1
C     DSQRT(X)=SQRT(X)
C
      IT=IF
      IF=IF+IDEL
      K=PMX+1
      K1=PMX+1-M
      L1=LMX-M
      DO 10 L=M,L1
      B2=2*L+1.
      B1=-(L+M-1)/B2
      B2=-(L+2-M)/B2
      DO 9 J1=1,K1
      I=IF+K*(J1-1)
      I1=I+IDEL+K+1
      I2=I+2*IDEL
      DO 8 J2=1,K1
      F(I)=F(I1)+B1*F(I)+B2*F(I2)
      F(I)=-F(I)
      IF(M.EQ.1)F(I)=2*F(I)
      I=I+1
      I1=I1+1
    8 I2=I2+1
    9 CONTINUE
   10 IF=IF+IDEL
      GO TO 15
C*
      ENTRY BLM(IF)
C*    =============
      IT=IF
      K1=IDEL1-M
      L1=LMX-M
      DO 14 L=M,L1
      B1=L+M
      B2=(L-M+1)/B1
      B1=B1*(2*L-1)
      B1=(2*L+1)*(L+M-1)/B1
      B1=-DSQRT(B1)
      B2=DSQRT(B2)
      I=IF
      I2=I+IDEL1+1
      DO 13 J1=1,K1
      F(I)=B1*F(I)+B2*F(I2)
      I=I+1
   13 I2=I2+1
   14 IF=IF+IDEL1
   15 IF=IT
      RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE HHMOL(HEV,NR,ICORE)
C*
C*------------------------------------------------------------------
C*
CQP      3
C     IMPLICIT REAL*16(A-H),INTEGER(P)
C     GENERIC
C     REAL*16 R,RH,RHE,SEPAT,IZ1,IZ2,X
C*ASB
C     REAL*16 teta,tp
CDP      2
      IMPLICIT REAL*8(A-H),INTEGER(P)
      REAL*8 R,RH,RHE,SEPAT,IZ1,IZ2,X
C*ASB 1   15.09.1993
      INCLUDE 'h2upar2.for'
C*ASB
      REAL*8 teta,tp,pi
      CHARACTER*256 filenm
C*
      INTEGER OUT
      COMMON F(1)
      COMMON /AREAR/A1,A2,B1,B2
C*ASB 1  21.07.93
c      COMMON /AREAI/LMX,PMX,MAXPN,NOFTE,ITRIP,IPAR,I(600)
      COMMON /AREAI/LMX,PMX,MAXPN,NOFTE,ITRIP,IPAR,I(2*nrdim)
C*
      COMMON /AREA1R/A(36)
      COMMON /AREA1I/N
      COMMON /AREA2/J(30)
      COMMON /AREA3/JA(20)
      COMMON /AREA4R/C(110)
      COMMON /AREA5R/D(7)
      COMMON /AREA5I/M
      COMMON /AREA6R/X(2)
      COMMON /AREA6I/LTE(10),KMX
      COMMON /AREA7R/R,IZ1,IZ2,SEPAT
      COMMON /AREA7I/IPNCH,IPRFUN
      COMMON /AREA8R/CONVE,E
      COMMON /AREA8I/ITER,NIN,IASY,NOFT1,ITEST,NOFT2,NOFT3,
     CISTP,ISTP1,IDGPR,IO(7),ITRANS
      COMMON/AREA9/ IBI(8)
      COMMON/S/IOO,IO2,NS,NEGV,NDEL,IHESHE,INPVC
      COMMON/S2/CINP(10)
      COMMON/FPS/MXRN,NAP,IAP
      common/teta/teta
      DIMENSION HEV(NR,1),EP(100),TP(100)
      DIMENSION INPUT(10)
      DIMENSION AB1(4)
C*ASB  1  21.07.93
c      DIMENSION IPOW(300),AB(4)
      DIMENSION IPOW(nrdim),AB(4)
C*
      DIMENSION EDP(3),ABOP(3)
C*ASB      DIMENSION OUT(400),PHIA(100),PHIB(100)
      DIMENSION CF1(nrdim),CF2(nrdim),EHE(nrdim)
      EQUIVALENCE (INPUT(1),ITER)
C*ASB 1  21.07.93
c      EQUIVALENCE (IPOW(1),I(301)),(AB(1),A1)
      EQUIVALENCE (IPOW(1),I(nrdim+1)),(AB(1),A1)
C*
      DATA NPRIM,NEMAT,NCOEFF,NMAT,NINT /8,9,10,18,23/
      DATA IPR1,IPR2,IPR3,IPR4/4*0/
      NAMELIST /NAM1/X,AB,R,E,CONVE,I,IPOW,LMX,PMX,MAXPN,NOFTE,ITRIP,
     CIPAR,IZ1,IZ2,ITER,NIN,IASY,NOFT1,ITEST,NOFT2,NOFT3,ISTP,ISTP1,
     CKMX,IDGPR,IPR1,IPR2,IPR3,IPR4,SEPAT,MXRN,PRM,DPAR
       NAMELIST/NAM2/ MAXPN,ITRIP,IPAR,IO,ITRANS,IPR1,IPR3,NOFTE,
     $NOFT1,IPR4,IOO,IO2
C*
      CHARACTER*256 filnam
C*
C*----------------------------------------------------------------------
C*
C*    X(1) is the convergence parameter for the matrix element
C*    calculation. The convergence is checked in subroutine LTEST.
C*
CQP      1
C2600 X(1)=1.E26
C*
CDP      1
 2600 X(1)=1.E13
C*
C*
C*
C*    Open the scratch files 9, 18, and 23 and result file 10:
C*
      CALL getenv('F10DAT',filnam)      
      OPEN ( UNIT=10 , FILE = filnam, STATUS='UNKNOWN',
     +       FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C*
      CALL getenv('F09DAT',filnam)      
      OPEN ( UNIT=9 , FILE = filnam, STATUS='UNKNOWN',
     +       FORM='UNFORMATTED', ACCESS='SEQUENTIAL' )
C*
      CALL getenv('F18DAT',filnam)      
      OPEN ( UNIT=18 , FILE = filnam, STATUS='UNKNOWN',
     +       FORM='UNFORMATTED', ACCESS='SEQUENTIAL' )
C*
      CALL getenv('F23DAT',filnam)      
      OPEN ( UNIT=23 , FILE = filnam, STATUS='UNKNOWN',
     +       FORM='UNFORMATTED', ACCESS='SEQUENTIAL' )
C*
C*    The following part is needed in (automatic) parameter 
C*    optimisation. Thus it is skipped in the first call of the
C*    program by setting IPR4=0:
C*
      IPR4=0
    1 IF(IPR4.NE.2)GO TO 102
      READ 14,E,ITER
      READ 11,SEPAT
      READ 9,IPR4
      GO TO 30
C*
C*    Starting point for a standard calculation (no parameter 
C*    optimisation):
C*
C*    ***** First, the parameters are read from the input file: *****
C*
  102 READ 9,ITRANS
C
C       ITRANS=0 for standard energy computation.
C       ( ITRANS is the first parameter in the input file. )
C
      IF (ITRANS.NE.0) THEN
        WRITE(*,*) ' '
        WRITE(*,*) '   ITRANS unequal 0, this is not implemented !!!'
        WRITE(*,*) ' '
        RETURN
      ENDIF
C*
C*
C*    Read the non-linear parameters A1=\alpha, A2=\bar{\alpha}, 
C*    B1=\beta, B2=\bar{\beta}, and the scaling parameter teta=\theta.
C*    (For usual runs choose \theta=1.0, to which it is set by 
C*     default, if it is found to be zero).
C*
  104 READ 1978,A1,A2,B1,B2,teta
CASB 8  29.03.1995
      IF (teta.EQ.0.0D+00) THEN
        WRITE(*,*) ' '
        WRITE(*,*) '  ***** Theta = 0.0 !?!     *****'
        WRITE(*,*) '  ***** Theta is set to 1.0 *****'
        teta = 1.0D+00
      ENDIF
      a1=a1*teta
      a2=a2*teta
CASB
C*
C*
C*    Read the value of the internuclear distance R:
C*
      READ 11,R
C*
C*
C*    Read NOFTE=number of basis functions, NIN=?, NDEL=?
C*       (If NIN.NE.0 then No INtegrals are calculated).
C*
      READ 52,NOFTE,NIN,NDEL
      IF(NDEL.EQ.0)NDEL=1
C*
C*
C*    Read the NOFTE basis functions defined by the 5 integer 
C*    exponents \mu, r, s, \bar{r}, \bar{s} that are stored 
C*    as a five-digit integer in array I.
C*    (since the non-linear parameters \alpha, \beta etc. are 
C*     the same for all basis functions, two basis functions 
C*     having the same 5 integer exponents are identical).
C*
      READ 13,(I(K),K=1,NOFTE)
      DO 105 K=1,NOFTE
      KK=K+1
      DO 105 KKK=KK,NOFTE
      IF(I(K)-I(KKK).NE.0)GO TO 109
  108 WRITE(6,106)I(K),K,KKK
  106 FORMAT(///10X,'IDENTICAL TERMS IN WAVEFUNCTION:',3I7)
      RETURN
C
  109 CONTINUE
  105 CONTINUE
C*
C*
C*    Read E (Energy for initial guess when only one eigenvalue 
C*            is determined which is achieved by setting IOO=1)
C*         ITER = ? (probably number of iterations in iterative
C*                   determination of one eigenvalue).
C*        [ ICC = 0(1) chooses real (complex) diagonalisation 
C*                                       (probably not active) ]
C*         IWRITE = 0: usual calculation
C*                = 1: store overlap, kinetic-energy, and 
C*                     potential-energy matrices.
C*                = 2: do not calculate, but read overlap, kinetic-
C*                     energy, and potential-energy matrix (useful 
C*                     e.g. for restarts after using IWRITE=1 in 
C*                     the first (interrupted) calculation). 
C*         IPRKS = ?  (some additional output is produced)
C*
      READ 14,E,ITER,ICC,IWRITE,IPRKS
C*
C*
C*    Read ITRIP and IPAR, but both parameters have no meaning for 
C*         this program, and thus they are set to zero.
C*
      READ 15,ITRIP,IPAR
      ITRIP = 0
      IPAR  = 0
C*
C*
C*    Read LMX, PMX, KMX, MAXPN which are the parameters for the 
C*    integrations, i.e. number of steps in numerical 1D-integration, 
C*    highest exponents of \xi and \eta in the calcuation of the basic 
C*    integrals which is the upper limit of values used in the 
C*    von-Neumann expansion:
C*    ( LMX  : maximum number of terms used in the von-Neumann expansion,
C*             i.e. the \Phi integrals are calculated in the range 
C*             \Phi_{i,j}^{0..LMX} and the B integrals in the range 
C*             B_{j}^{0..LMX}.
C*      PMX  : maximum index of \Phi_{i,j}^l, i.e. the \Phi integrals 
C*             are calculated in the range \Phi_{0..PMX,0..PMX}^{0..LMX}. 
C*      KMX  : maximum index of the B integral to be calculated, i.e. 
C*             the B integrals are calculated in the range 
C*             B_{0..KMX}^{0..LMX}.
C*      MAXPN: number of grid points in numerical integration)
C*
      READ 16,LMX,PMX,KMX,MAXPN
C*
C*
C*    Read CONVE the energy convergence criterion:
C*
      READ 11,CONVE
C*
C*
C*    Read SEPAT:  energy in the separated atom limit (ground state). 
C*                 In the case of H_2 (the ground state dissociating 
C*                 into H(1s)+H(1s)) its value should be -1.00 (= two
C*                 times the energy of H(1s)). In the case of HeH^+
C*                 (the ground state dissociating into He + H^+) its
C*                 values should be the one of He(1s^2)=-2.903...
C*                 This parameter is only used for the evaluation of
C*                 the dissociation energy D given in the output and 
C*                 thus it is not really of importance, since it does
C*                 not affect the calculation itself.
C* 
C*         INPVC (a parameter that defines whether CI coefficients 
C*                are given as input or not; here INPVC=0 means that 
C*                no vectors are given which is the default. Probably
C*                input vectors are given for parameter optimisation 
C*                or to speed up the calculation by using the vectors
C*                e.g. obtained for one value of R for another (nearby)
C*                value of R).
C*
      READ 14,SEPAT,INPVC
C*
C*
C*    Read IPNCH = 1 write informations (in DIAGDR) into xxx.f10 file,
C*         IPR1 > 0 prints the (basic) integrals, (IPR1=2 only a 
C*                  subgroup of the integrals).
C*         IPR2 > 0 skips the matrix evaluation,
C*         IPR3 > 0 skips the matrix diagonalisation,
C          IPRFUN = 1(0) linear coefficients are printed (or not).
C*
      READ 19,IPNCH,IPR1,IPR2,IPR3,IPRFUN
      WRITE(6,*) 'Punch parameter =',ipnch
C*
C*
C*    Read IZ1 = nuclear charge of first nucleus (nucleus A),
C*         IZ2 = nuclear charge of second nucleus (nucleus B).
C*
      READ 57,IZ1,IZ2
C*
C*
C*    Read IOO = 0 complete diagonalisation of the Hamiltonian
C*             = 1 only one eigenvalue is determined,
C*         IO2,NAP,IAP = ? (If NAP unequal 0 or 1, diagonalisation 
C*                          is skipped, but it has additional use, since
C*                          certain integral calculations are skipped 
C*                          depending on the values NAP and IAP).
C*         NEGV = number of eigenvectors to be written to the xxx.f10 
C*                file (should be set to ONE for usual runs, since 
C*                otherwise up to 8 vectors are printed).
C*
C*ASB 5  03.09.1998
      IF (NOFTE.LT.1000) THEN
        READ 5200,IOO,IO2,NEGV,NAP,IAP
      ELSE
        READ 5202,IOO,IO2,NEGV,NAP,IAP
      ENDIF
C*
c      READ 5200,IOO,IO2,NEGV,NAP,IAP
C***
      IF(NAP.EQ.0)NAP=1
      IF(NOFTE.LE.8.AND.NEGV.EQ.0)NEGV=NOFTE
      IF(NOFTE.GT.8.AND.NEGV.EQ.0)NEGV=8
      IF(IOO.NE.0)NEGV=1
C*
C*
C*    Read MXRN, PRM, DPAR = ? (These parameters are used in the 
C*    parameter optimisation. If MXRN<>0, then parameter optimisation
C*    is performed).
C*
      READ 18,MXRN,PRM,DPAR
C*
C*
C*    Read INPUT(3): switch between Sigma (0) and Pi (1) symmetry,
C*                   leads however in this program version only
C*                   to the exciting statement that Pi states are
C*                   calculated, but it remains a Sigma-state 
C*                   calculation... (thus it is set to zero).
C*         INPUT(4..10) = ? 
C*    (This array is identical (via equivalencing) to:
C*     INPUT(1)=ITER, (2)=NIN, (3)=IASY, (4)=NOFT1, (5)=ITEST,
C*     INPUT(6)=NOFT2, (7)=NOFT3, (8)=ISTP, (9)=ISTP1, (10)=IDGPR)
C*
      READ 17,(INPUT(K),K=3,10)
      INPUT(3) = 0
C*
C*
C*    Read the input vectors, if requested by INPVC.NE.0 (and given 
C*    in the input file):
C*
      IF(INPVC.EQ.0)GO TO 2
      READ 905,(CINP(K),K=1,NOFTE)
C*
C*
C*    Read IPR4 (Another parameter related to parameter optimisation.
C*    It seems as this parameter indicates whether the parameter 
C*    optimisation is ready or not).
C*
    2 READ 9,IPR4
  101 IF(NIN.LT.0)GO TO 30
      X(2)=0.
      N1=1
      K6=1
      EDP(3)=0.
      K5=121
      K4=NOFT1+120
      IF(IWRITE.EQ.2) GO TO 122
C
C       ***** Computation of the necessary integrals *****
C
    3 CONTINUE
      ALF=2*A1
      DALF=A2-A1
      IF(IPR1.LT.0)GO TO 4
C*
C*       Calculate the \Phi integrals (defined by Ruedenberg) for 
C*       the possible combinations of \alpha and the case that 
C*       a 1/r_{12} term (requiring a von-Neumann expansion) is 
C*       present (however, the summation over the von-Neumann 
C*       terms is done when forming the matrix elements):
C*
      CALL INTFI(MAXPN,PMX,LMX,ALF,DALF)
C*
C*       Calculate the \Phi(\alpha1,\alpha2) integrals (for the trivial 
C*       case that no 1/r_{12} term is present) and also the B(\beta) 
C*       integrals (but them for all cases):
C*
    4 CALL INTPAK
C*
      K2=NOFT2
      IF(K2.LT.NOFTE)K2=NOFTE
      K=(K2*(K2+1))/2
      ISMOLD=J(18)
      DO 42 L=19,20
   42 J(L)=J(L-1)+K
C*
C*    J(18), J(19), J(20) contain the addresses of the first elements of
C*    S,     V, and K, respectively (which are, however, not yet 
C*    calculated).
C*
      J(20)=ICORE-K+1
      J(19)=J(20)-K
      J(18)=J(19)-K
C
C*    ICOEF is the size of the space in the blank COMMON for 
C*    vectors used in the Ostrowski diagonalisation:
C*
C*ASB 1  21.07.93
c      ICOEF=600
      ICOEF=2*nrdim
C*
C*    Check whether there is enough space for the integrals:
C*
      IF(ISMOLD.GT.J(18)) THEN
        WRITE(6,*) ' ***** ERROR: ISMOLD .GT. J(18) ***** ' 
        GO TO 6
      ENDIF
C
C     Check whether there is also enough space for the Ostrowski
C*    diagonalisation:
C
      IF(J(18)-K.GT.ICOEF) THEN 
        GO TO 5
      ELSE
        WRITE(6,*) ' ***** ERROR: J(18)-K .LT. ICOEF ***** ' 
      ENDIF
    6 CONTINUE
C*
      WRITE(6,928) J(18),K,ICORE,ICOEF,NNN,ISMOLD
      WRITE(6,56)
      RETURN
C*
    5 IF(IPR1.LE.0) GO TO 20
C*
C*    *** Print integrals, if required ***
C*
      WRITE(6,930) ISMOLD
      K=1
      K1=ISMOLD-1
      IF(IPR1.EQ.2)K=J(12)
      WRITE(6,50)(F(L),L=K,K1)
   20 IF(IPR2.LT.0)GO TO 30
C
C       ***** Form the matrices *****
C
      CALL MDRIV(N1,AB1)
C
      K2=J(18)
      K3=K2+J(19)-J(18)-1
      IF(IPR2.EQ.0)GO TO 555
      WRITE(6,10) (F(K),K=K2,K3)
C
  555 CONTINUE
      IF(ITRANS.NE.0) GO TO 31
      IF(F(1).NE.0.)GO TO 22
      WRITE(6,*) ' ***  F(1) = 0 ERROR  *** '
      WRITE(6,56)
      GO TO 31
   22 IF(IPR2.EQ.0)GO TO 30
      DO 29 L=18,20
      WRITE(6,54)
      K2=J(L)
      K3=K2+(NOFTE*(NOFTE+1))/2
      K3=K3-1
   29 WRITE(6,50)(F(K),K=K2,K3)
      WRITE(6,51)X(2)
   30 IF(IPR3.LT.0)GO TO 31
      IF(X(2).NE.0.0)WRITE(6,55)X(2)
C
C        ***** Diagonalise the Hamiltonian *****
C
      IKM=J(20)
      IVM=J(19)
      ISM=J(18)
      NN2=(NOFTE*(NOFTE+1))/2
      LKM=IKM+NN2-1
      LVM=IVM+NN2-1
      LSM=ISM+NN2-1
      IF(IWRITE.EQ.0) GO TO 111
      WRITE(18) IKM,IVM,ISM,LKM,LVM,LSM
      WRITE(18) (F(II),II=IKM,LKM)
      WRITE(18) (F(II),II=IVM,LVM)
      WRITE(18) (F(II),II=ISM,LSM)
  111 CONTINUE
      GO TO 133
  122 READ(18) IKM,IVM,ISM,LKM,LVM,LSM
      READ(18) (F(II),II=IKM,LKM)
      READ(18) (F(II),II=IVM,LVM)
      READ(18) (F(II),II=ISM,LSM)
      J(20)=IKM
      J(19)=IVM
      J(18)=ISM
  133 CONTINUE
      IF(NAP.NE.1)GO TO 31
C*
      CALL DIAGDR(HEV,NR,N1,ITEST)
C*
      IF(MXRN.LE.0)GO TO 39
C
C        ***** Optimise exponents, if required *****
C
      MXRN=MXRN-1
      EDP(K6)=E
      ABOP(K6)=AB(PRM)
      GO TO (34,35,37),K6
   34 AB(PRM)=ABOP(1)+DPAR
      GO TO 40
   35 IF(E.GT.EDP(1))GO TO 36
      K6=1
      EDP(3)=EDP(1)
      ABOP(1)=ABOP(2)
      EDP(1)=E
      GO TO 34
   36 AB(PRM)=ABOP(1)-DPAR
      IF(EDP(3).NE.0)GO TO 38
      GO TO 40
   37 IF(E.GT.EDP(1))GO TO 38
      EDP(2)=EDP(1)
      ABOP(2)=ABOP(1)
      EDP(1)=E
      ABOP(1)=AB(PRM)
      EDP(3)=0.
      K6=2
      GO TO 36
   40 K6=K6+1
      IF(DALF.EQ.0.)AB(PRM+1)=AB(PRM)
      IF(DALF.EQ.0.AND.PRM.EQ.3) AB(4)=-AB(4)
      GO TO (3,3,4,4),PRM
   38 AB(PRM)=ABOP(1)+(EDP(3)-EDP(2))*DPAR/2/(EDP(3)+EDP(2)-2*EDP(1))
      IF(DALF.EQ.0)AB(PRM+1)=AB(PRM)
      IF(DALF.EQ.0.AND.PRM.EQ.3) AB(4)=-AB(4)
      IF(IPR4.NE.0)GO TO 381
      WRITE(6,11)AB(PRM)
      GO TO 31
  381 READ 18,MXRN,PRM,DPAR
      READ 9,IPR4
      GO TO 101
C
   39 CONTINUE
      IF(K4.LE.120.OR.ITEST.EQ.0) GO TO 31
      IF(NOFTE.GE.NOFT2.OR.NOFT3.EQ.0)GO TO 31
      N1=NOFTE+1
      NIN=N1
      DO 32 K=1,NOFT3
      NOFTE=NOFTE+1
      I(NOFTE)=I(K5)
      K5=K5+1
      K4=K4-1
C     IF(K4.EQ.120)NOFT3=0
      IF(NOFTE.EQ.NOFT2.OR.K4.LE.120)GO TO 33
   32 CONTINUE
C     CHANGED 2/18/85 KSZ
   33 IF(K4.EQ.120) NOFT3=0
      GO TO 4
   31 IF(IPR4.NE.0)GO TO 1
      RETURN
C
    9 FORMAT(I1)
   10 FORMAT(4D15.6)
 1978 FORMAT(4D15.6,D12.6)
   11 FORMAT(D20.10)
   12 FORMAT(2I3)
   13 FORMAT(10I6)
   14 FORMAT(D15.6,4I1)
   15 FORMAT(2I1)
   16 FORMAT(4I3)
   17 FORMAT(8I3)
   18 FORMAT(2I1,D15.5)
   19 FORMAT(5I1)
   56 FORMAT(1H1,11HINPUT ERROR)
   50 FORMAT(1H ,6D19.10)
   51 FORMAT(1H ,E14.7)
   52 FORMAT(3I4)
 5200 FORMAT(5I3)
C*ASB 1  03.09.1998
 5202 FORMAT(2I3,I4,2I3)
C***
   54 FORMAT(1H0)
   55 FORMAT(///'0LMX REACHED X=',E10.2////)
   57 FORMAT(2D15.8)
  900 FORMAT(3D25.14)
 8900 FORMAT(2D25.14)
  901 FORMAT(9X,I3,3F15.8,F10.3) 
 8901 FORMAT(1X,I3,3(F12.5,F12.5),F8.2,F8.2) 
 8903 FORMAT(1X,I3,D24.17,2X,D24.17,4X,D24.17,2X,D24.17) 
  902 FORMAT(3I1)
  903 FORMAT(10I6)
  904 FORMAT(I3,F10.4,4F10.6,4I3)
 9044 FORMAT(I3,F10.4,4F10.6,4I3,F10.6)
  905 FORMAT(4D20.9)
 8905 FORMAT(4D20.9)
  906 FORMAT(D20.9)
 8906 FORMAT(2D20.9)
  908 FORMAT(5X,'GAINESVILLE 1984'///)
  910 FORMAT(32X,'INPUT INFORMATION')
  911 FORMAT(32X,'-----------------'//)
  912 FORMAT(26X,'INTERNUCLEAR DISTANCE =',F8.4//)
  914 FORMAT(31X,'ENERGY =',F11.7/)
  916 FORMAT(10X,'EXPONENTS :',4F9.4,' TETA = ',F9.4/)
  917 FORMAT(29X,I3,'-TERM WAVEFUNCTION')
  918 FORMAT(9X,10I6)
  919 FORMAT(/31X,'PARAMETERS FOR HEH+'/)
  920 FORMAT(//37X,'OUTPUT')
  921 FORMAT(37X,'------'/)
  923 FORMAT(1H+,29X,I3)
  924 FORMAT(/)
  926 FORMAT(10I5)
  927 FORMAT(10I1)
  928 FORMAT(' ISM, K, ICORE, ICOEF, NNN, ISMOLD ',6I10)
  930 FORMAT(' BASIC INTEGRALS'/' NUMBER OF INTS =',I10/)
      END
C*
C*============================================================
C*
      SUBROUTINE ILOCZ
C*
C*------------------------------------------------------------
C*
C*    Fill the array A(1..36) which is transferred via COMMON 
C*    AREA1R into e.g. I00,I02,I11,...I100 with values by 
C*    forming up to 36 products of two F values (i.e. two 
C*    integrals).
C*
C*------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(A-F)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-F)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON F(1)
      COMMON /AREA1R/A(36)
      COMMON /AREA1I/N
      COMMON /AREA2/IA,IB
C*
C*----------------------------------------------------------
C*
      K=1
      L=0
    2 J=IB+L
      I=IA
    3 A(K)=F(I)*F(J)
      K=K+1
      I=I+1
      J=J-1
      IF(J-IB)4,3,3
    4 L=L+2
      IF(L-N)2,2,5
    5 RETURN
      END
C*
C*================================================================
C*
      SUBROUTINE INTFI(M1,P,LMX,ALF,DALF)
C*
C*----------------------------------------------------------------
C*
C*    Calculation of the \Phi_{n,\bar{n}}^l function (defined by 
C*    Ruedenberg) for M=0 (\Sigma symmetry) and the case where a 
C*    1/r_{12} term has to be considered which couples \xi_1 and 
C*    \xi_2 and thus requires to split the double integration into 
C*    two, depending on whether \xi_1<\xi_2 or vice versa.
C*
C*    ***  The maximum index n or \bar{n} is defined by P, the *** 
C*    ***  maximum value of l is defined by LMX.               ***
C*    ***  Thus \Phi_{1..P,1..P}^{1..LMX} is calculated!!!     ***
C*
C*    The simple case of no 1/r_{1,2} term is considered in FIEV.
C*
C*    The resulting integrals are written to file fort.23 and 
C*    later on read into the blank COMMON (and thus into array F)
C*    in subroutine INTPAK.
C*
C*--------------------------------------------------------------- 
C*
CQP      2
C     IMPLICIT  INTEGER(P),REAL*16(A-H,S-Z)
C     GENERIC
CDP      1
      IMPLICIT  INTEGER(P),REAL*8(A-H,S-Z)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON FA(20,20),B(901,2),SN(902),SI(900,21),TS(901,9)
C*
C*----------------------------------------------------------------
C*
CQP      3
C     DLOG(X)=LOG(X)
C     DEXP(X)=EXP(X)
C     DABS(X)=ABS(X)
C*
      REWIND 23
      IRUN=1
      IRUN1=1
      IZERO=0
      IONE=1
      M=M1+1
      TS(1,1)=1/DFLOAT(M1)
      L=M-2
      DO 2 I=2,L,1
    2 TS(I,1)=TS(I-1,1)+TS(1,1)
      TS(L+1,1)=TS(L,1)+TS(1,1)/2.
      TS(M,1)=TS(L,1)+TS(1,1)
      DO 3 I=1,M,1
      TS(I,6)=1.
      TS(I,7)=1.
      TS(I,2)=1/TS(I,1)
      TS(I,3)=1-TS(I,1)**2
      T1=(P+2)*DLOG(TS(I,2))-ALF*TS(I,2)
      TS(I,4)=DEXP(T1)
      IF(DALF.EQ.0)GO TO 3
      T2=DALF*TS(I,2)
      TS(I,5)=DEXP(T1-T2)
      TS(I,8)=DEXP(T1-2*T2)
    3 CONTINUE
C     WRITE(6,1000)  ((TS(I,J),J=1,8),I=1,M1)
 1000 FORMAT(1X,8E14.7)
C
    6 H1=(TS(1,1)/3)**3
C*
C*        Begin FI
C*
      P=P+1
      L1=LMX+1
      DO 45 L2=1,L1,1
C*
C*       Compute BL
C*
      L=L2-1
      T1=2*L
      T2=T1
      T3=0
      DO 12 I=1,M,1
      B(I,1)=1
   12 B(I,2)=1.
      I2=2
      I1=1
      IF(L-1)17,151,13
   13 I1=1
      DO 15 J=1,L,1
      DO 14 I=1,M,1
   14 B(I,I1)=(T2*B(I,I2)+T3*TS(I,3)*B(I,I1))/T1
      T1=T1-1.
      T2=T2-2.
      T3=T3+1.
      I3=I1
      I1=I2
   15 I2=I3
  151 TS(M,7)=TS(M,7)*TS(M-3,1)
      IF(DABS(B(M-3,I2)).LE.1.E-20)WRITE(6,1001) B(M-3,I2),M,I2,L2,L,
     $I1,I3
 1001 FORMAT(1X,'151 B(M-3,I2),M,I2,L2,L,I1,I3',E14.7,1X,6I5)
      T1=B(M,I2)/B(M-3,I2)*TS(M,7)
      DO 16 I=3,M,1
      TS(I-2,7)=TS(I-2,7)*TS(I-2,1)*TS(I,2)
      TS(I-1,6)=TS(I-1,6)*TS(I-1,1)*TS(I,2)
      IF(DABS(B(I-1,I2)).LE.1.E-20)WRITE(6,1002) B(I-1,I2),M,I2,L2,L,
     $I1,I3,I
      IF(DABS(B(I-2,I2)).LE.1.E-20)WRITE(6,1003) B(I-2,I2),M,I2,L2,L,
     $I1,I3,I
 1002 FORMAT(1X,'B(I-1,I2),M,I2,L2,L,I1,I3,I',E14.7,7I5)
 1003 FORMAT(1X,'B(I-2,I2),M,I2,L2,L,I1,I3,I',E14.7,7I5)
      B(I-1,I1)=TS(I-1,6)*B(I,I2)/B(I-1,I2)
   16 B(I-2,I2)=TS(I-2,7)*B(I,I2)/B(I-2,I2)
      B(M,I2)=T1
C*
C*         BL ready
C*
   17 J2=3
      J3=1
      T1=0.
      IF(DALF.EQ.T1) J2=2
      SN(1)=0.
      DO 18 I=1,M,1
   18 SN(I+1)=TS(I,4)
      GO TO 24
C*
   20 DO 23 I=1,M,1
      IF(J2-2)22,21,21
   21 SN(I+1)=TS(I,5)
      GO TO 23
C*
   22 SN(I+1)=TS(I,8)
   23 CONTINUE
   24 P1=P
   25 P4=P1
      IF(J2.EQ.IONE)P4=P+1
      SI(M-1,P4)=(SN(M-1)+4*SN(M)*B(M-2,I1)+SN(M+1)*B(M-2,I2))/2.
      SI(M-2,P4)=SN(M-2)+4*SN(M-1)*B(M-3,I1)+SN(M+1)*B(M,I2)
      K=M-3
      GO TO 28
C*
   27 K=M-4
   28 SI(K,P4)=SN(K)+4*SN(K+1)*B(K-1,I1)+(SN(K+2)+SI(K+2,P4))*B(K-1,I2)
      K=K-2
      IF(K-1)27,291,28
  291 IF(L.NE.IZERO)GO TO 29
      SI(1,P4)=4*SN(2)+SN(3)+SI(3,P4)
   29 IF(J2-2)30,32,32
   30 J3=P1
      GO TO 37
C*
   32 P1=P1-1
      IF(P1.EQ.IZERO)GO TO 34
      DO 33 I=1,M,1
   33 SN(I+1)=SN(I+1)*TS(I,1)
      GO TO 25
C*
   34 IF(J2-2)35,37,36
C*
C*    Some of the \Phi functions are now written to file fort.23:
C*
   35 WRITE(23)((FA(I,J),J=1,P),I=1,P)
      IRUN1=IRUN1+1
      IRUN=IRUN1
      J2=2
      J3=1
      GO TO 20
C*
   36 J2=1
      GO TO 20
C*
   37 DO 44 P2=1,P,1
      P3=J3
  370 P4=P3
      IF(J2.EQ.1) P4=P+1
   38 T1=0
C
      K=M-1
   39 T1=T1+SI(K,P2)*SI(K,P4)/TS(K-1,3)
      K=K-2
      IF(K-1)40,41,39
   40 T1=2*T1
      K=M-2
      GO TO 39
C*
   41 T1=2*T1
      IF(L.NE.IZERO)GO TO 42
      T1=T1+SI(1,P2)*SI(1,P4)
   42 FA(P2,P3)=T1*H1
      IF(J2.EQ.1)GO TO 44
   43 FA(P3,P2)=FA(P2,P3)
      P3=P3+1
      IF(P3.LE.P2) GO TO 370
   44 CONTINUE
      IF(J2.EQ.1)GO TO 32
C*
C*    More of the \Phi functions are written to file fort.23:
C*
      WRITE(23)((FA(I,J),J=1,P),I=1,P)
      IRUN1=IRUN1+1
      IRUN=IRUN1
   45 CONTINUE
      P=P-1
      RETURN
      END
C*
C*=================================================================
C*
      SUBROUTINE INTPAK
C*
C*-----------------------------------------------------------------
C*
CQP      2
C     IMPLICIT INTEGER(P),REAL*16(A-F)
C     GENERIC
CDP      1
      IMPLICIT INTEGER(P),REAL*8(A-F)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      INTEGER COSH
      COMMON F(1)
      COMMON /AREAR/A1,A2,B1,B2
      COMMON /AREAI/LMX,PMX
      COMMON /AREA6I/LTES,L,IDEL,IDEL1,LT(6),KMX
      COMMON /AREA3/IT(19),COSH
      COMMON /AREA2/IA,IB,I,P1,P2,P1M1,P01,P02,IF,IFA,IFS,IBI(6)
     C   ,ISM
C*
C*----------------------------------------------------------------
C*
      REWIND 23
      IRUN=1
      IRUN1=1
C*
C*
C*           Determine starting addresses for the integrals: 
C*           ===============================================
C*
      IDEL1 = PMX + 1
      IDEL = IDEL1**2
      I = ( LMX + 2 ) * IDEL
C*
C* 
C*           ###  \Phi_{0,0} (2*\alpha, 2*\bar{\alpha}) at IFA = 1  ###
C*
      IFA = 1
C*
C*
C*           ###  \Phi_{0,0} ( \alpha+\bar{\alpha}, \alpha+\bar{\alpha} )
C*                at IFS = IFA = 1,             if \alpha = \bar{\alpha}
C*                at IFS = IFA + ( LMX + 2 ) * ( PMX + 1 )^2   otherwise ###
C*
      IFS = 1
      DALF = A2 - A1
      IF (DALF.NE.0.0) IFS = IFS + I
C*
      IF = 0
C*
C*
C*          ###  B_0^0 ( 2 * \beta )  
C*               at IBI(1) = IFS + ( LMX + 2 ) * ( PMX + 1 )^2  ### 
C*
      IBI(1) = IFS + I
C*
      IDEL1 = KMX + 1
C*
C*
C*       Check whether \beta = \bar{\beta} = 0.0. 
C*          If this is the case, set COSH=0, otherwise set COSH=1. 
C*          (If COSH=0, only one type of B integrals occur, i.e. 
C*           B_{l,k} (\beta=0.0) ):
C*
      COSH = 0
      IF (B1) 6,5,6
    5 IF (B2) 6,7,6
C*
    6 COSH = 1
      IF = IDEL1 * ( LMX + 1 )
C*
C*          If COSH = 1:
C*          ------------
C*
C*          ###  B_0^0 ( 2 * \bar{\beta} )  
C*                   at IBI(2) = IBI(1) + [ (KMX+1) * (LMX+1) ]  ### 
C*          ###  B_0^0 ( \beta + \bar{\beta} ) 
C*                   at IBI(3) = IBI(2) + [ (KMX+1) * (LMX+1) ]  ###
C*          ###  B_0^0 ( 0.000 ) 
C*                   at IBI(4) = IBI(3) + [ (KMX+1) * (LMX+1) ]  ###
C*          ###  B_0^0 ( \beta - \bar{\beta} ) 
C*                   at IBI(5) = IBI(4) + [ (KMX+1) * (LMX+1) ]  ###
C*          ###  B_0^0 ( \bar{\beta} - \beta ) 
C*                   at IBI(6) = IBI(5) + [ (KMX+1) * (LMX+1) ]  ###
C*
C*          If COSH = 0:
C*          ------------
C*
C*          ###  IBI(2) = IBI(3) = IBI(4) = IBI(5) = IBI(6) = IBI(1) ###
C*
    7 DO 8 J=2,6
    8 IBI(J) = IBI(J-1) + IF
C*
C*
C*         The starting point for the overlap matrix is finally 
C*         determined to be ISM = IBI(6) + [ (KMX+1) * (LMX+1) ]
C*  
      ISM = IBI(6) + IDEL1 * ( LMX + 1 )
C*

      IF = IBI(1) - 6
      F(IF) = 2 * B1
      F(IF+1) = 2 * B2
      F(IF+2) = B1 + B2
      F(IF+3) = 0.
      F(IF+4) = B1 - B2
      F(IF+5) = B2 - B1
C*
C*    Calculation of the integrals over \eta, i.e.  
C*    B_{0..KMX}^{0..LMX} (\eta), for all possible combinations of 
C*    the combined \beta exponents, i.e. 
C*    2*\beta, 2*\bar{\beta}, \beta+\bar{\beta}, 
C*    \beta-\beta=\bar{\beta}-\bar{\beta}=0.0, \beta-\bar{\beta}, 
C*    and \bar{\beta}-\beta. (The fact is used that the values for 
C*    the integral with exponent \beta is +1 or -1 times the 
C*    integral with exponent -\beta.)
C*    If \beta=\bar{\beta}=0 (COSH=0), then the DO loop 
C*    is left before termination (using GOTO 10). 
C*    The integrals are stored starting at F(IBI(1..6)), i.e. 
C*    the integrals for 2*\beta are stored at F(IBI(1)...), 
C*    the ones for 2*\bar{\beta} at F(IBI(2)...) etc. 
C*
      DO 9 J=1,6
        IA = IBI(J)
        BET = F(IF)
        CALL EINT(BET,LMX,KMX,IA)
        IF (COSH.EQ.0) GO TO 10
    9 IF=IF+1
C*
C*
C*    Calculate the simple integrals over \xi_1 and \xi_2 that 
C*    occur in cases where no 1/r_{12} term is present. In that 
C*    case the double integral splits into a simple product of 
C*    two 1D integrals which is calculated in routine FIEV.
C* 
C*       1st call: Calculate \Phi( 2*\alpha, 2*\bar{\alpha} ):
C*       --------- (and store these Asymmetric functions at 
C*                  F(IFA...))
C*
   10 ALF = 2 * A1
      CALL FIEV(ALF,DALF,PMX,IFA)
C*
      IF (DALF.EQ.0.0) GO TO 11
C*
C*       2nd call: Calculate 
C*       ---------   \Phi( \alpha + \bar{\alpha}, \bar{\alpha} + \alpha )
C*                 (only needed, if \alpha <> \bar{\alpha}):
C*                 (and store these Symmetric functions at
C*                  F(IFS...), while, if \alpha = \bar{\alpha}, 
C*                  IFS=IFA)
C*
      ALF=A1+A2
      DALF=0.
      CALL FIEV(ALF,DALF,PMX,IFS)
C*
   11 IF=IFA+IDEL
      IF1=IFS+IDEL
   12 L=0
      I=IF+IDEL-1
      I1=IF1+IDEL-1
   13 READ(23)(F(J),J=IF,I)
      IRUN1=IRUN1+1
      IRUN=IRUN1
      L=L+1
      IF=IF+IDEL
      I=I+IDEL
      IF(IFA.EQ.IFS)GO TO 14
      READ(23)(F(J),J=IF1,I1)
      IRUN1=IRUN1+1
      IRUN=IRUN1
      IF1=IF1+IDEL
      I1=I1+IDEL
   14 IF(L-LMX)13,13,16
   16 RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE JA3B
C*
C*------------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(A-F),INTEGER(P)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-F),INTEGER(P)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON F(1)
      COMMON /AREA2/IA,IB,I1,P1,P2,P1M1,P01,P02
CQP      1
C     REAL*16 I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,I33,I42,
CDP      1
      REAL*8 I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,I33,I42,
     CI51,I60,I08,I17,I26,I35,I44,I53,I62,I71,I80,I010,I19,I28,
     CI37,I46,I55,I64,I73,I82,I91,I100
      COMMON /AREA1R/I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,
     CI33,I42,I51,I60,I08,I17,I26,I35,I44,I53,I62,I71,I80,I010,I19,I28,
     CI37,I46,I55,I64,I73,I82,I91,I100
      COMMON /AREA1I/N
C*
C*---------------------------------------------------------------------
C*
      I=I1
      N=10
C*
C*    Fill the variables I00,I02,I11,...I100 with values by calling
C*    subroutine ILOCZ in which up to 36 products of two F values 
C*    (i.e. two integrals) are formed and stored in the variables 
C*    I00 etc.:
C*
      CALL ILOCZ
C*
      B=(30*(I42+I24)-20*I22-12*(I62+I26)-36*I44+I82+9*(I64+I46)+I28)
     C*F(I)
      I=I+P2
      B=B+(20*I02-30*I04-24*I42+12*I06+8*I62-I08+27*I44-6*I64-6*I46)
     C*F(I)
      I=I-P1M1
      B=B+(36*(I53-I33+I35)-6*(I73+I37)-24*I55)*F(I)
      I=I-P1M1
      B=B+(20*I20-30*I40-24*I24+27*I44+12*I60+8*I26-6*(I64+I46)-I80)*
     CF(I)
      I=I+P02
      B=B+(36*I40-30*I20+24*I22-27*I42-9*I60+6*I62)*F(I)
      I=I+P1M1
      B=B+(36*(I31-I51)-18*I35+6*I71+12*I55)*F(I)
      I=I+P1M1
      B=B+(24*(I04+I40)-20*I00-8*I06-48*I44-8*I60+18*(I64+I46))*F(I)
      I=I+P1M1
      B=B+(36*(I13-I15)-18*I53+6*I17+12*I55)*F(I)
      I=I+P1M1
      B=B+(24*I22-30*I02+36*I04-27*I24-9*I06+6*I26)*F(I)
      I=I+P2
      B=B+(12*I02-8*I22-9*I04+6*I24)*F(I)
      I=I-P1M1
      B=B+(18*I33-36*I13+24*I15-12*I35)*F(I)
      I=I-P1M1
      B=B+(30*I00-24*I20-27*I04+48*I24+6*I06-18*I26)*F(I)
      I=I-P1M1
      B=B+(18*(I15+I51)-36*I11-20*I55)*F(I)
      I=I-P1M1
      B=B+(30*I00-24*I02-27*I40+48*I42+6*I60-18*I62)*F(I)
      I=I-P1M1
      B=B+(24*I51-36*I31+18*I33-12*I53)*F(I)
      I=I-P1M1
      B=B+(12*I20-9*I40+6*I42-8*I22)*F(I)
      I=I+P02
      B=B-I20*F(I)
      I=I+P1M1
      B=B+6*I31*F(I)
      I=I+P1M1
      B=B+(8*I02-12*I00+6*I40-18*I42)*F(I)
      I=I+P1M1
      B=B+(36*I11-18*I13-12*I51+20*I53)*F(I)
      I=I+P1M1
      B=B+(27*(I20+I02)-36*I00-48*I22)*F(I)
      I=I+P1M1
      B=B+(36*I11-18*I31-12*I15+20*I35)*F(I)
      I=I+P1M1
      B=B+(8*I20-12*I00+6*I04-18*I24)*F(I)
      I=I+P1M1
      B=B+6*I13*F(I)
      I=I+P1M1
      B=B-I02*F(I)
      I=I+P02
      B=B+I00*F(I)
      I=I-P1M1
      B=B-6*I11*F(I)
      I=I-P1M1
      B=B+(9*I00-6*(I20+I02)+18*I22)*F(I)
      I=I-P1M1
      B=B+(12*(I31+I13)-24*I11-20*I33)*F(I)
      I=I-P1M1
      B=B+(9*I00-6*(I02+I20)+18*I22)*F(I)
      I=I-P1M1
      B=B-6*I11*F(I)
      I=I-P1M1
      B=B+I00*F(I)
      I00=B
      RETURN
      END
C*
C*====================================================================
C*
      SUBROUTINE KA
C*
C*--------------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(A-F),INTEGER(P)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-F),INTEGER(P)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON F(10000)
      COMMON /AREA2/IA,IB,I1,P1,P2,P1M1,P01,P02
CQP      1
C     REAL*16 I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,I33,I42,
CDP      1
      REAL*8 I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,I33,I42,
     CI51,I60,I08,I17,I26,I35,I44,I53,I62,I71,I80,I010,I19,I28,
     CI37,I46,I55,I64,I73,I82,I91,I100
      COMMON /AREA1R/I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,
     CI33,I42,I51,I60,I08,I17,I26,I35,I44,I53,I62,I71,I80,I010,I19,I28,
     CI37,I46,I55,I64,I73,I82,I91,I100
      COMMON /AREA1I/N
C*
C*---------------------------------------------------------------------
C*
      I=I1
      N=4
C*
C*    Fill the variables I00,I02,I11,...I100 with values by calling
C*    subroutine ILOCZ in which up to 36 products of two F values 
C*    (i.e. two integrals) are formed and stored in the variables 
C*    I00 etc.:
C*
      CALL ILOCZ
C*
      B=(2.0*I02-I22-I04)*F(I)
      I=I+P2
      B=B-I02*F(I)
      I=I-P1M1
      B=B+2.0*I13*F(I)
      I=I-P1M1
      B=B+(I20-2.0*I00)*F(I)
      I=I+P2
      B=B+I00*F(I)
      I=I-P1M1
      B=B-2.0*I11*F(I)
      I=I-P1M1
      B=B+I00*F(I)
      GO TO 1
C*
      ENTRY JA
C*    ========
      I=I1
      N=6
C*
C*    Fill the variables I00,I02,I11,...I100 with values by calling
C*    subroutine ILOCZ in which up to 36 products of two F values 
C*    (i.e. two integrals) are formed and stored in the variables 
C*    I00 etc.:
C*
      CALL ILOCZ
C*
      B=(I42-2.0*I22+I24)*F(I)
      I=I+P2
      B=B+(2.0*I02-I04)*F(I)
      I=I-P1M1
      B=B-2.0*I33*F(I)
      I=I-P1M1
      B=B+(2.0*I20-I40)*F(I)
      I=I+P02
      B=B-I20*F(I)
      I=I+P1M1
      B=B+2.0*I31*F(I)
      I=I+P1M1
      B=B-2.*I00*F(I)
      I=I+P1M1
      B=B+2.*I13*F(I)
      I=I+P1M1
      B=B-I02*F(I)
      I=I+P02
      B=B+I00*F(I)
      I=I-P1M1
      B=B-2.*I11*F(I)
      I=I-P1M1
      B=B+I00*F(I)
      GO TO 1
C*
      ENTRY KA2B
C*    ==========
      I=I1
      N=6
C*
C*    Fill the variables I00,I02,I11,...I100 with values by calling
C*    subroutine ILOCZ in which up to 36 products of two F values 
C*    (i.e. two integrals) are formed and stored in the variables 
C*    I00 etc.:
C*
      CALL ILOCZ
C*
      B=(6.*(-I02+I22+I04)-I42-4.*I24-I06)*F(I)
      I=I+P2
      B=B+(6.*I02-4.*I22-4.*I04+2.*I24)*F(I)
      I=I-P1M1
      B=B+4.*(-2.*I13+I33+I15)*F(I)
      I=I-P1M1
      B=B+(6.*I00-6.*I20+I40-3.*I04+2.*I24)*F(I)
      I=I+P02
      B=B+(-6.*I00+4.*I20+3.*I02-2.*I22)*F(I)
      I=I+P1M1
      B=B+(8.*I11-4.*I31)*F(I)
      I=I+P1M1
      B=B+(-6.*I00+4.*I20+2.*I04-6.*I24)*F(I)
      I=I+P1M1
      B=B+4.*I13*F(I)
      I=I+P1M1
      B=B-I02*F(I)
      I=I+P02
      B=B+I00*F(I)
      I=I-P1M1
      B=B-4.*I11*F(I)
      I=I-P1M1
      B=B+(4.*I00-2.*I20-2.*I02+6.*I22)*F(I)
      I=I-P1M1
      B=B-4.*I11*F(I)
      I=I-P1M1
      B=B+I00*F(I)
      GO TO 1
C*
      ENTRY JA2B
C*    ==========
      I=I1
      N=8
C*
C*    Fill the variables I00,I02,I11,...I100 with values by calling
C*    subroutine ILOCZ in which up to 36 products of two F values 
C*    (i.e. two integrals) are formed and stored in the variables 
C*    I00 etc.:
C*
      CALL ILOCZ
C*
      B=(6.*(I22-I42-I24)+I62+4.*I44+I26)*F(I)
      I=I+P2
      B=B+(-6.*I02+6.*I04+3.*I42-I06-2.*I44)*F(I)
      I=I-P1M1
      B=B+4.*(2.*I33-I53-I35)*F(I)
      I=I-P1M1
      B=B+(-6.*I20+6.*I40-I60+3.*I24-2.*I44)*F(I)
      I=I+P02
      B=B+(6.*I20-4.*I40-3.*I22+2.*I42)*F(I)
      I=I+P1M1
      B=B+(-8.*I31+4.*I51)*F(I)
      I=I+P1M1
      B=B+(6.*I00-3.*I40-3.*I04+6.*I44)*F(I)
      I=I+P1M1
      B=B+(-8.*I13+4.*I15)*F(I)
      I=I+P1M1
      B=B+(6.*I02+2.*I24-3.*I22-4.*I04)*F(I)
      I=I+P2
      B=B-I02*F(I)
      I=I-P1M1
      B=B+4.*I13*F(I)
      I=I-P1M1
      B=B+(-6.*I00+3.*I20+2.*I04-6.*I24)*F(I)
      I=I-P1M1
      B=B+8.*I11*F(I)
      I=I-P1M1
      B=B+(-6.*I00+3.*I02+2.*I40-6.*I42)*F(I)
      I=I-P1M1
      B=B+4.*I31*F(I)
      I=I-P1M1
      B=B-I20*F(I)
      I=I+P2
      B=B+I00*F(I)
      I=I+P1M1
      B=B-4.*I11*F(I)
      I=I+P1M1
      B=B+(4.*I00-2.*I20-2.*I02+6.*I22)*F(I)
      I=I+P1M1
      B=B-4.*I11*F(I)
      I=I+P1M1
      B=B+I00*F(I)
    1 I00=B
      RETURN
      END
C*
C*===================================================================
C*
      SUBROUTINE KZW
C*
C*-------------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(A-F),INTEGER(P)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-F),INTEGER(P)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON F(10000)
      COMMON /AREA2/IA,IB,I1,P1,P2,P1M1,P01,P02
CQP      1
C     REAL*16 I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,I33,I42,
CDP      1
      REAL*8 I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,I33,I42,
     CI51,I60,I08,I17,I26,I35,I44,I53,I62,I71,I80,I010,I19,I28,
     CI37,I46,I55,I64,I73,I82,I91,I100
      COMMON /AREA1R/I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,
     CI33,I42,I51,I60,I08,I17,I26,I35,I44,I53,I62,I71,I80,I010,I19,I28,
     CI37,I46,I55,I64,I73,I82,I91,I100
      COMMON /AREA1I/N
C*
C*--------------------------------------------------------------------
C*
      I = I1
      N = I + P02
      B = ( F(N) * F(IB) - F(I) * F(IB+2) ) * F(IA)
      GO TO 1
C*
      ENTRY JZW
C*    =========
      I=I1
      N=I+P02
      B=(F(I)*F(IB+2)-F(N)*F(IB))*F(IA+2)
      I=I+P2
      N=N+P2
      B=B-(F(I)*F(IB+2)-F(N)*F(IB))*F(IA)
      GO TO 1
C*
      ENTRY KAS1
C*    ==========
      I=I1
      N=6
C*
C*    Fill the variables I00,I02,I11,...I100 with values by calling
C*    subroutine ILOCZ in which up to 36 products of two F values 
C*    (i.e. two integrals) are formed and stored in the variables 
C*    I00 etc.:
C*
      CALL ILOCZ
C*
      B=(5*(-I02+I22+I04)-I42-3*I24-I06)*F(I)
      I=I+P2
      B=B+(5*I02-3*(I22+I04)+I24)*F(I)
      I=I-P1M1
      B=B+4*(I33+I15-2*I13)*F(I)
      I=I-P1M1
      B=B+(5*(I00-I20)-2*I04+I40+I24)*F(I)
      I=I+P02
      B=B+(-5*I00+3*I20+2*I02-I22)*F(I)
      I=I+P1M1
      B=B+(8*I11-4*I31)*F(I)
      I=I+P1M1
      B=B+(3*I20-5*(I00+I24)+I04)*F(I)
      I=I+P1M1
      B=B+4*I13*F(I)
      I=I+P1M1
      B=B-I02*F(I)
      I=I+P02
      B=B+I00*F(I)
      I=I-P1M1
      B=B-4*I11*F(I)
      I=I-P1M1
      B=B+(3*I00+5*I22-I20-I02)*F(I)
      I=I-P1M1
      B=B-4*I11*F(I)
      I=I-P1M1
      B=B+I00*F(I)
      B=3.0*B
      GO TO 1
C*
      ENTRY JAS1
C*    ==========
      I=I1
      N=8
C*
C*    Fill the variables I00,I02,I11,...I100 with values by calling
C*    subroutine ILOCZ in which up to 36 products of two F values 
C*    (i.e. two integrals) are formed and stored in the variables 
C*    I00 etc.:
C*
      CALL ILOCZ
C*
      B=(5*(I22-I42-I24)+I62+3*I44+I26)*F(I)
      I=I+P2
      B=B+(2*I42-5*(I02-I04)-I06-I44)*F(I)
      I=I-P1M1
      B=B+4*(2*I33-I53-I35)*F(I)
      I=I-P1M1
      B=B+(5*(I40-I20)+2*I24-I60-I44)*F(I)
      I=I+P02
      B=B+(5*I20-3*I40-2*I22+I42)*F(I)
      I=I+P1M1
      B=B+(4*I51-8*I31)*F(I)
      I=I+P1M1
      B=B+(5*(I00+I44)-2*(I04+I40))*F(I)
      I=I+P1M1
      B=B+(4*I15-8*I13)*F(I)
      I=I+P1M1
      B=B+(5*I02-2*I22-3*I04+I24)*F(I)
      I=I+P02
      B=B+(2*I20-5*(I00+I24)+I04)*F(I)
      I=I-P1M1
      B=B+8*I11*F(I)
      I=I-P1M1
      B=B+(2*I02-5*(I00+I42)+I40)*F(I)
      I=I-P1M1
      B=B+4*I31*F(I)
      I=I-P1M1
      B=B-I20*F(I)
      I=I+P2
      B=B+I00*F(I)
      I=I+P1M1
      B=B-4*I11*F(I)
      I=I+P1M1
      B=B+(3*I00+5*I22-I20-I02)*F(I)
      I=I+P1M1
      B=B-4*I11*F(I)
      I=I+P1M1
      B=B+I00*F(I)
      I=I-P02
      B=B-I02*F(I)
      I=I-P1M1
      B=B+4*I13*F(I)
      B=3.0*B
    1 I00=B
      RETURN
      END
C*
C*===================================================================
C*
      SUBROUTINE KA3B
C*
C*-------------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(A-F),INTEGER(P)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-F),INTEGER(P)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON F(10000)
      COMMON /AREA2/IA,IB,I1,P1,P2,P1M1,P01,P02
CQP      1
C     REAL*16 I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,I33,I42,
CDP      1
      REAL*8 I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,I33,I42,
     CI51,I60,I08,I17,I26,I35,I44,I53,I62,I71,I80,I010,I19,I28,
     CI37,I46,I55,I64,I73,I82,I91,I100
      COMMON /AREA1R/I00,I02,I11,I20,I04,I13,I22,I31,I40,I06,I15,I24,
     CI33,I42,I51,I60,I08,I17,I26,I35,I44,I53,I62,I71,I80,I010,I19,I28,
     CI37,I46,I55,I64,I73,I82,I91,I100
      COMMON /AREA1I/N
C*
C*---------------------------------------------------------------------
C*
      I=I1
      N=8
C*
C*    Fill the variables I00,I02,I11,...I100 with values by calling
C*    subroutine ILOCZ in which up to 36 products of two F values 
C*    (i.e. two integrals) are formed and stored in the variables 
C*    I00 etc.:
C*
      CALL ILOCZ
C*
      B=(20*I02-30*(I22+I04)+12*(I42+I06)+36*I24-I62-9*(I44+I26)-I08)
     C*F(I)
      I=I+P2
      B=B+(36*(I22+I04-I24)-30*I02-9*(I42+I06)+6*(I44+I26))*F(I)
      I=I-P1M1
      B=B+(36*(I13-I33-I15)+6*I53+24*I35+6*I17)*F(I)
      I=I-P1M1
      B=B+(30*I20-20*I00+24*I04-27*I24-12*I40-8*I06+6*(I44+I26)+I60)*
     CF(I)
      I=I+P02
      B=B+(30*I00-36*I20-24*I02+27*I22+9*I40-6*I42)*F(I)
      I=I+P1M1
      B=B+(36*(I31-I11)+18*I15-6*I51-12*I35)*F(I)
      I=I+P1M1
      B=B+(30*I00-36*I20-27*I04+54*I24+9*I40-18*I44+6*I06-18*I26)*F(I)
      I=I+P1M1
      B=B+(24*(I33+I15)-36*I13-12*I35)*F(I)
      I=I+P1M1
      B=B+(12*I02-9*(I22+I04)+6*I24)*F(I)
      I=I+P2
      B=B-I02*F(I)
      I=I-P1M1
      B=B+6*I13*F(I)
      I=I-P1M1
      B=B+(9*I20-12*I00+6*I04-18*I24)*F(I)
      I=I-P1M1
      B=B+(36*I11-24*I31-12*I15+20*I35)*F(I)
      I=I-P1M1
      B=B+(36*(I20-I00)+27*I02-54*I22-6*I40+18*I42)*F(I)
      I=I-P1M1
      B=B+(36*I11-24*I31-18*I13+12*I33)*F(I)
      I=I-P1M1
      B=B+(9*I20-12*I00+8*I02-6*I22)*F(I)
      I=I+P02
      B=B+I00*F(I)
      I=I+P1M1
      B=B-6*I11*F(I)
      I=I+P1M1
      B=B+(9*I00-6*(I20+I02)+18*I22)*F(I)
      I=I+P1M1
      B=B+(12*(I31+I13)-24*I11-20*I33)*F(I)
      I=I+P1M1
      B=B+(9*I00-6*(I20+I02)+18*I22)*F(I)
      I=I+P1M1
      B=B-6*I11*F(I)
      I=I+P1M1
      B=B+I00*F(I)
      GO TO 1
    1 I00=B
      RETURN
      END
C*
C*===================================================================
C*
      SUBROUTINE LTEST(LL)
C*
C*-------------------------------------------------------------------
C*
C*     Check whether the matrix elements are converged with respect to
C*     the von-Neumann expansion.
C*
C*-------------------------------------------------------------------
C*
CQP      4
C     IMPLICIT REAL*16(X,Y)
C     GENERIC
C     REAL*16 CONV,ALAR1,A,B
C     REAL*16 DABS
CDP      3
      IMPLICIT REAL*8(X,Y)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      REAL*8 CONV,ALAR1,A,B
      REAL*8 DABS
      COMMON /AREA2/IA,IB,I
      COMMON /AREA5R/X,Y
      COMMON /AREA5I/M
      COMMON /AREAR/XA(4)
      COMMON /AREAI/LMAX
      COMMON /AREA6R/CONV,ALAR1
      COMMON /AREA6I/LTES,L,IDEL,IDEL1
C*
C*-----------------------------------------------------------------
C*
CQP      1
C     DABS(X)=ABS(X)
      LL=0
      LMX=LMAX-M
      A= X*CONV
      A=DABS(A)
      B=DABS(Y)
C                        CHECK IF ABS(DEL(Y)/Y) GT 1/CONV
C                        NO MEANS CONVERGED
C                        YES MEANS NOT CONVERGED
      IF(A.GT.B)GO TO 1
C                        THRESHOLD MET. MUST BE MET TWICE.
      LTES=LTES+1
C                        IF THE THRESHOLD MET FOR THE FIRST
C                        TIME, PERFORM ONE MORE STEP
      GO TO (2,5),LTES
C                        THRESHOLD NOT MET
    1 LTES=0
C                        CHECK FOR THE MAX NUMBER OF STEPS
    2 IF(L.GT.LMX)GO TO 3
      L=L+1
      I=I+IDEL
      IA=IA+IDEL1
      IB=IB+IDEL1
      GO TO 4
    3 A=A/B
C                       STORE LARGEST VALUE OF THE ERROR AND RETURN
      IF(A.LE.ALAR1)GO TO 5
      ALAR1=A
      GO TO 5
    4 RETURN
C                       THIS RETURN MEANS END OF SUMMATION
C                       CAN BE CONVERGED OR NOT
    5 LL=1
      RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE MDRIV(N1,AB)
C*
C*------------------------------------------------------------------
C*
CQP      3
C     IMPLICIT REAL*16(A-F),INTEGER(P,Z)
C     GENERIC
C     REAL*16 R,Z1,Z2
CDP      
      IMPLICIT REAL*8(A-F),INTEGER(P,Z)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      REAL*8 R,Z1,Z2
      INTEGER COSH
      COMMON F(1)
      COMMON /AREAR/B(4)
      COMMON /AREAI/L(3),NOFTE,N(2),POWR(60)
      COMMON /AREA2/IT(8),IF,IFA,IFS,IB(6),ISM,IVM,IKM 
      COMMON /AREA7R/ R,Z1,Z2
      COMMON /AREA3/ II(19),COSH
      COMMON /AREA4R/ C(24)
      COMMON /AREA5R/ AA(2)
      COMMON /AREA5I/ M
      COMMON /AREA6I/ IIA(7)
      COMMON /AREA8R/ DT(2)
      COMMON /AREA8I/ ITTE(14),IADIAB,IT2(2),ITRANS
      COMMON /AREA9/IBI(8)
      DIMENSION AB(4)
C*
C*-----------------------------------------------------------------
C*
      N1=1
      IF(IIA(5).EQ.-1) GO TO 13
      C(23)=-4*(Z1+Z2)
      C(24)=4*(Z1-Z2)
      I=NOFTE*(NOFTE+1)/2
      I1=((N1-1)*N1)/2+1
      I2=IVM-ISM
      K=ISM-1
C*     WRITE(6,200)I1,I2,K,ISM,L(1),L(2),L(3),NOFTE,N(1),N(2)
      DO 2 I3=1,3
      DO 1 K1=I1,I
      J=K+K1
    1 F(J)=0.
    2 K=K+I2
      DO 3 I=6,7
   3  IIA(I)=0
  200 FORMAT(10I5)
      KIVM=IVM+I
      KIKM=IKM+I
      KISM=ISM+I
C*
C      WRITE(6,202)
C      WRITE(6,201)(F(KKK),KKK=IVM,KIVM)
C      WRITE(6,202)
C      WRITE(6,201)(F(KKK),KKK=IKM,KIKM)
C      WRITE(6,202)
C      WRITE(6,201)(F(KKK),KKK=ISM,KISM)
C*
  201 FORMAT(1H ,7E14.7)
  202 FORMAT(/'  MATRIX IN MDRIV')
   31 CONTINUE
C*
C*
C*       Check for the largest exponent \mu (= exponent of r_{12}).
C*       The largest value of \mu is assigned to M1. 
C*       If M1 is found to be larger than the maximum value of \mu that    
C*         can be handled by the program (which seems to be 3 and not 2, 
C*         as I thought???), the program stops with an error message (by 
C*         setting F(1)=0.0).
C*       Otherwise, M1 defines the number of calls of subroutine MX which 
C*         actually calculates the integrals.
C*
      K=POWR(1)
      DO 5 I=2,NOFTE
        IF (K-POWR(I)) 4,5,5
    4   K=POWR(I)
    5 CONTINUE
      M1=K/10000
      IF(M1.LE.3) GO TO 51
      F(1)=0.
      GO TO 12
C*
C*
C*       Call subroutine MX(n1) (M1-2) times, but at least once
C*       (which is done by initialising M=0, and after calling 
C*        MX(n1) increasing M by 1; if M1-M is then smaller or equal 
C*        1, no further calls of MX(n1) are done). (The actual value
C*        of M is transfered into MX(n1) via COMMON block AREA5I):
C*
   51 M=0
      GO TO 9
C*
    6 CALL FLM(IFA)
C*
      IF(IFS-IFA)71,71,7
    7 CALL FLM(IFS)
   71 DO 8 J=1,6
        IF=IB(J)
C*
        CALL BLM(IF)
C*
        IF (COSH) 8,9,8
    8 CONTINUE
C*
    9 CALL MX(N1)
C*
      M = M + 1
      IF (M1-M) 10,6,6
C*
C*       Add the nuclear repulsion term to the potential
C*       energy matrix elements:
C*       (due to the non-orthogonality of the basis functions,
C*        the nuclear repulsion term has to be added to all 
C*        matrix elements (weighted by the overlap matrix element)).
C*
   10 A=Z1*Z2
      IF(IIA(5).EQ.-1) GO TO 12
      A=A/2
      I=NOFTE*(NOFTE+1)/2
      DO 11 J=I1,I
      K1=ISM+J-1
      K=IVM+J-1
C***
c      IF ((K.EQ.IVM).OR.(K.EQ.IVM+1).OR.(K.EQ.IVM+2)) THEN
c        WRITE(*,*) ' Checkpt.4a: K, F(K) = ',K,F(K)
c      ENDIF
C***
      F(K)=F(K)+A*F(K1)
C***
c      IF ((K.EQ.IVM).OR.(K.EQ.IVM+1).OR.(K.EQ.IVM+2)) THEN
c        WRITE(*,*) ' Checkpt.4b: K, F(K) = ',K,F(K)
c      ENDIF
C***
   11 CONTINUE    
      GO TO 12
C*
   13 A1=AB(1)
      A2=AB(2)
      A3=AB(3)
      A4=AB(4)
      K=NOFTE*ITTE(4)
      I=K+ISM-1
      DO 101 J=ISM,I
  101 F(J)=0.
      K=POWR(1)
      I1=2
      J=NOFTE
      DO 150 K1=1,2
      IF(J.EQ.1) I1=1
      DO 105 I=I1,J
      IF(K-POWR(I))104,105,105
  104 K=POWR(I)
  105 CONTINUE
C*ASB 2  21.07.93
c      I1=301
      I1=nrdim+1
c  150 J=ITTE(4)+300
  150 J=ITTE(4)+nrdim
C*
      M1=K/10000
      IF(M1.LE.3) GO TO 151
      F(1)=0.
      GO TO 12
C*
  151 M=0
      GO TO 109
C*
  106 CALL FLM(IFA)
C*
      IF(IFS-IFA)171,171,107
C*
  107 CALL FLM(IFS)
C*
  171 DO 108 J=1,8
      IF(COSH.EQ.0.AND.J.GT.2) GO TO 108
      IF=IBI(J)
C*
      CALL BLM(IF)
C*
  108 CONTINUE
C*
  109 CONTINUE
      WRITE(*,*) ' '
      WRITE(*,*) '    *****   ERROR *****'
      WRITE(*,*) '    This place of the program should not '
      WRITE(*,*) '    be reached!!! '
      WRITE(*,*) '    Change parameter ITRANS to zero!!! '
      WRITE(*,*) ' '
      STOP
C*
      M=M+1
      IF(M.LE.M1) GO TO 106
   12 RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE MTRI
C*
C*------------------------------------------------------------------
C*
CQP      1
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON /AREA2/IA,IB
      COMMON /AREA6I/I(4),NOS,NOT,NOV
C*
C*-----------------------------------------------------------------
C*
C*       Save IA in K and IB in K1 for later use in this routine:
C*       (probably IA and IB are modified in the subroutines KIN, 
C*        SMAT and VMAT):
C*
      K=IA
      K1=IB
C*
C*        Evidently, there had been a program version that allowed 
C*        (via parameters NOT, NOS, and NOV) to selectively omit the 
C*        calculation of the kinetic-energy matrix T (NOT= no T matrix),
C*        of the overlap matrix S (NOS), and the matrix of the potential
C*        energy V (NOV):
C*
C      IF(NOT.NE.0) GO TO 1
C*
C*    Calculate a fragment of the kinetic energy matrix element:
C*
      CALL KIN
C*
      IA=K
      IB=K1
C*
C    1 IF(NOS.NE.0) GO TO 2
C*
C*    Calculate a fragment of the overlap matrix element:
C*
      CALL SMAT
C*
      IA=K
      IB=K1
C*
C    2 IF(NOV.NE.0) GO TO 3
C*
C*    Calculate a fragment of the potential energy matrix element:
C*
      CALL VMAT
C*
C*         Make sure to leave the values IA and IB unchanged after
C*         this subroutine:
C* 
      IA = K
      IB = K1
C*
    3 RETURN
      END
C*
C*================================================================
C*
      SUBROUTINE MX(N1)
C*
C*----------------------------------------------------------------
C*
CQP      3
C     IMPLICIT REAL*16(A-F,X),INTEGER(P,R,S,T)
C     GENERIC
C     REAL*16 IZ1,IZ2
C*
C*ASB 1   15.09.1993
      INCLUDE 'h2upar2.for'
C*
CDP      2
      IMPLICIT REAL*8(A-F,X),INTEGER(P,R,S,T)
      REAL*8 IZ1,IZ2
C*
C*ASB 2   09.06.1998
      INTEGER AD(nrdim)     
c      INTEGER AD(100)
C*
      INTEGER COSH
      COMMON /AREAR/A1,A2,B1,B2
C*
C*ASB 2   09.06.1998
      COMMON /AREAI/LMX,PMX,MAXPN,NOFTE,TRIPL,PARIT,POWR(nrdim)
c      COMMON /AREAI/LMX,PMX,MAXPN,NOFTE,TRIPL,PARIT,POWR(60)
C*
      COMMON/AREA3/TEST,TEST1,TEST2,MPN,MMN,RPR,RMR,SPS,SMS,
     C             MU(2),R(2),S(2),RB(2),SB(2),COSH
      COMMON /AREA2/IA,IB,T,P10,P20,P1M1,P01,P02,IF,IFA,IFB,
     C              IB1,IB2,IB3,IB4,IB5,IB6,ISM,IVM,IKM
      COMMON /AREA4R/C(24),APA,AMA,BPB,BMB
      COMMON /AREA5R/X(2)
      COMMON /AREA5I/M,TES1,TES2,J1,J2
      COMMON /AREA7R/AR,IZ1,IZ2
      COMMON /AREA8R/D(2)
      COMMON /AREA8I/III(2),IASY
      COMMON /FPS/  MXRN,NAP,IAP
      COMMON /CALLS/ICALL
C*
      EQUIVALENCE (AD(1),POWR(1)),(P,PMX)
      EQUIVALENCE (IFS,IFB)
      EQUIVALENCE (NPM,MPN)
C*
C*-----------------------------------------------------------------
C*
      AA1=2*A1
      AA2=2*A2
      BB1=2*B1
      BB2=2*B2
      J1=((N1-1)*N1)/2-1
C*
      DO 19 I=N1,NOFTE
        K1=2
        K=I
C*
C*
C*         Analyse the integer tupels of the two basis functions
C*         (note that AD(i)=POWR(i) due to equivalencing!).
C*         The integer tupel of one basis function is defined via
C*         the integers  MU(=\mu), R(=r), S(=s), RB(=\bar{r}), and 
C*         SB(=\bar{s}):
C*
    3   MU(K1)=AD(K)/10000
        L=AD(K)-10000*MU(K1)
        R(K1)=L/1000
        L=L-1000*R(K1)
        S(K1)=L/100
        L=L-100*S(K1)
        RB(K1)=L/10
        SB(K1)=L-10*RB(K1)
C*
C*         If only the integer tupel of the second basis function i 
C*         (K1=2) is read in yet, go to 4 and read the integer tupel 
C*         defining the first basis function j (K1=1). 
C*         Otherwise, if already the tupels for both basis functions 
C*         have been read in (K1=1), continue with the integral 
C*         calculation by going to label 5: 
        GO TO(5,4),K1
C*
    4   DO 18 J=1,I
          J1=J1+1
          JS=(J1+1)/NAP
          JS=(J1+1)-JS*NAP
          IF(JS.NE.IAP)GO TO 18
C*
C*           Now read the integer tupel defining the first basis 
C*           function j (K1=1, GOTO 3):
          K1=1
          K=J
          GO TO 3
C*
C*
C*             Calculate MPN = NPM (equivalence!) = \mu_j + \mu_i,
C*                       MMN = \mu_j - \mu_i, and
C*
C*                              /  0  ,  if s_j + \bar{s}_j + PARIT odd  
C*                       K2  =  |
C*                              \  1  ,  if s_j + \bar{s}_j + PARIT even
    5     MPN = MU(1) + MU(2)
          MMN = MU(1) - MU(2)
          K1 = S(1) + SB(1) + PARIT
          K2 = K1 / 2
          K2 = K1 - 2 * K2
C*
C*
C*           In subroutine STODC the three variables TEST, TEST1, 
C*           and TEST2 are determined based on the actual values of 
C*           MPN(=\mu_j+\mu_i) and M(=number of calls of this subroutine,
C*           i.e. MX(ni)). The result is transfered back via COMMON block 
C*           AREA3):
C* 
          CALL STODC
C*
          M1=(NPM+1)/2
          IF (M1.LT.M) GO TO 18
C*
          TES1 = 2 * M1 - NPM
          TES2 = TES1 + 1
C*
C*
C*           Via equivalencing, P = PMX where PMX is an input parameter
C*           specifying the maximum argument of the \Phi integral, i.e. 
C*           \Phi(0..PMX,0..PMX). In order to avoid arrays starting with
C*           argument 0, this is mapped to \Phi(1..PMX+1,1..PMX+1). The 
C*           corresponding quantity of interest, i.e. PMX+1, is saved in
C*           P10  (while P20 = 2 * P10). The real range of arguments is 
C*           however - for the reasons stated above - P1M1 = PMX:
C*  
          P10 = P + 1
          P20 = 2 * P10
          P01 = 1
          P02 = 2
          P1M1 = P10 - P01
C*
          K1=0
C*
C*            IFA points to the first \Phi integral with doubled exponents,
C*                i.e. to \Phi_{0,0} ( 2*\alpha , 2*\bar{\alpha} ):
          IF = IFA
C*
C*            IB1 points to B_0^0 (2*\beta):
          IA = IB1
C*
C*            IB2 points to B_0^0 (2*\bar{\beta}):
          IB = IB2
C*
C*
C*             *****   First integral fragment  *****
C*
C*             (i.e. first call of CCI and MTRI [and CCHI and MTRI, 
C*              if needed]; K1=1)
C*
C*                APA = \alpha + \alpha = 2.0 * \alpha
C*                AMA = \alpha - \alpha = 0.0
C*                BPB = \beta  + \beta  = 2.0 * \beta
C*                BMB = \beta  - \beta  = 0.0
C*
C*                When calling CCI/MTRI:
C*                ----------------------
C*
C*                  J2 = 1
C*                  IF pointing on 
C*                      \Phi_{r_j+r_i,\bar{r}_j+\bar{r}_i}^0 
C*                                 ( 2 * \alpha , 2 * \bar{\alpha} )
C*                  IA pointing on 
C*                      B_{s_j+s_i}^0 ( 2 * \beta )
C*                  IB pointing on 
C*                      B_{\bar{s}_j+\bar{s}_i}^0 ( 2 * \bar{\beta} ) 
C*
C*                When calling CCHI/MTRI:
C*                -----------------------
C*
C*                  J2 = 1 + K2   (where K2 = 0 if s_j+\bar{s}_j+PARIT odd
C*                                   or  K2 = 1 if       "             even)
C*                  IF pointing on 
C*                      \Phi_{r_j+r_i,\bar{r}_j+\bar{r}_i}^0 
C*                                 ( 2 * \alpha , 2 * \bar{\alpha} )
C*                  IA pointing on 
C*                      B_{s_j+s_i}^0 ( 0.0D+00 )
C*                  IB pointing on 
C*                      B_{\bar{s}_j+\bar{s}_i}^0 ( 0.0D+00 ) 
C*
          AMA = 0.
          APA = AA1
          BMB = 0.
          BPB = BB1
          J2 = 1
          GO TO 11
C*
C*
C*             *****   Second integral fragment  *****
C*
C*             (i.e. second call of CCI and MTRI [and CCHI and MTRI, 
C*              if needed]; K1=2)
C*
C*                APA = \bar{\alpha} + \alpha
C*                AMA = \bar{\alpha} - \alpha
C*                BPB = \bar{\beta}  + \beta  
C*                BMB = \bar{\beta}  - \beta 
C*
C*                [ R(2) = r_i <--> \bar{r}_i = RB(2) ]
C*                [ S(2) = s_i <--> \bar{s}_i = SB(2) ]
C*                
C*                When calling CCI/MTRI:
C*                ----------------------
C*
C*                  J2 = TRIPL + 1
C*                  IF pointing on 
C*                      \Phi_{r_j+\bar{r}_i,\bar{r}_j+r_i}^0 
C*                              (\alpha+\bar{\alpha},\alpha+\bar{\alpha} )
C*                  IA pointing on 
C*                      B_{s_j+\bar{s}_i}^0 ( \beta + \bar{\beta} )
C*                  IB pointing on 
C*                      B_{\bar{s}_j+s_i}^0 ( \beta + \bar{\beta} ) 
C*                
C*                When calling CCHI/MTRI:
C*                -----------------------
C*
C*                  J2 = TRIPL + K2 + 1     (K2 = 0 or 1, see above)
C*                  IF pointing on 
C*                      \Phi_{r_j+\bar{r}_i,\bar{r}_j+r_i}^0 
C*                               (\alpha+\bar{\alpha},\alpha+\bar{\alpha} )
C*                  IA pointing on 
C*                      B_{s_j+\bar{s}_i}^0 ( \bar{\beta} - \beta )
C*                  IB pointing on 
C*                      B_{\bar{s}_j+s_i}^0 ( \beta - \bar{\beta} ) 
C*
C*
    6     AMA = A2 - A1
          APA = A2 + A1
          BMB = B2 - B1
          BPB = B2 + B1
C*
    7     CONTINUE
C*
C*            IFS points to the first \Phi integral with mixed alpha's, 
C*                i.e. to  
C*                \Phi_{0,0} (\alpha+\bar{\alpha},\alpha+\bar{\alpha}))
          IF=IFS
C*
C*            IB3 points to B_0^0 ( \beta + \bar{\beta} )
          IA=IB3
          IB=IB3
C*
          J2=TRIPL+1
          GO TO 72
C*
   72     K = 2
C*
C*           Exchange r and \bar{r} as well as s and \bar{s} for 
C*           basis function K (where K can be j(\equiv 1) or 
C*           i(\equiv 2)): 
    8     L=RB(K)
          RB(K)=R(K)
          R(K)=L
          L=SB(K)
          SB(K)=S(K)
          S(K)=L
          GO TO 11
C*
C*
C*             *****   Third integral fragment  *****
C*
C*             (i.e. third call of CCI and MTRI [and CCHI and MTRI, 
C*              if needed]; K1=3)
C*
C*                APA = \bar{\alpha} + \bar{\alpha} = 2.0 * \bar{\alpha}
C*                AMA = \bar{\alpha} - \bar{\alpha} = 0.0
C*                BPB = \bar{\beta}  + \bar{\beta}  = 2.0 * \bar{\beta}
C*                BMB = \bar{\beta}  - \bar{\beta}  = 0.0 
C*
C*                [ R(1) = r_j <--> \bar{r}_j = RB(1) ]
C*                [ S(1) = s_j <--> \bar{s}_j = SB(1) ]
C*                
C*                When calling CCI/MTRI:
C*                ----------------------
C*
C*                  J2 = 1
C*                  IF pointing on 
C*                      \Phi_{\bar{r}_j+\bar{r}_i,r_j+r_i}^0 
C*                                  ( 2 * \alpha, 2 * \bar{\alpha} )
C*                  IA pointing on 
C*                      B_{\bar{s}_j+\bar{s}_i}^0 ( 2 * \bar{\beta} )
C*                  IB pointing on 
C*                      B_{s_j+s_i}^0 ( 2 * \beta ) 
C*                
C*                When calling CCHI/MTRI:
C*                -----------------------
C*
C*                  J2 = 1 + K2     (where K2 = 0 or 1, see above)
C*                  IF pointing on 
C*                      \Phi_{\bar{r}_j+\bar{r}_i,r_j+r_i}^0 
C*                                  ( 2 * \alpha, 2 * \bar{\alpha} )
C*                  IA pointing on 
C*                      B_{\bar{s}_j+\bar{s}_i}^0 ( 0.0D+00 )
C*                  IB pointing on 
C*                      B_{s_j+s_i}^0 ( 0.0D+00 ) 
C*
    9     AMA = 0.
          APA = AA2
          BMB = 0.
          BPB = BB2
C*
          K = 1
C*
C*            IFA points to the first \Phi integral with doubled exponents,
C*                i.e. to \Phi_{0,0} ( 2*\alpha , 2*\bar{\alpha} ):
          IF = IFA
C*
C*            IB2 points to B_0^0 (2*\bar{\beta}):
          IA = IB2
C*
C*            IB1 points to B_0^0 (2*\beta):
          IB = IB1
C*
C*            Exchange P01 <--> P10 and P02 <--> P20:
          P01 = P10
          P02 = P20
          P10 = 1
          P20 = 2
          P1M1 = P10 - P01
          J2 = 1
          GO TO 8
C*
C*
C*             *****   Fourth integral fragment  *****
C*
C*             (i.e. fourth call of CCI and MTRI [and CCHI and MTRI, 
C*              if needed]; K1=4)
C*
C*                APA = \alpha + \bar{\alpha} 
C*                AMA = \alpha - \bar{\alpha} 
C*                BPB = \beta  + \bar{\beta}  
C*                BMB = \beta  - \bar{\beta}  
C*
C*                [ R(2) = r_i <--> \bar{r}_i = RB(2) ]
C*                [ S(2) = s_i <--> \bar{s}_i = SB(2) ]
C*
C*                When calling CCI/MTRI:
C*                ----------------------
C*
C*                  J2 = TRIPL + 1
C*                  IF pointing on 
C*                      \Phi_{\bar{r}_j+r_i,r_j+\bar{r}_i}^0 
C*                             ( \alpha+\bar{\alpha},\alpha+\bar{\alpha} )
C*                  IA pointing on 
C*                      B_{\bar{s}_j+s_i}^0 ( \beta + \bar{\beta} )
C*                  IB pointing on 
C*                      B_{s_j+\bar{s}_i}^0 ( \beta + \bar{\beta} ) 
C*
C*                When calling CCHI/MTRI:
C*                -----------------------
C*
C*                  J2 = TRIPL + 1 + K2   (where K2 = 0 or 1, see above)
C*                  IF pointing on 
C*                      \Phi_{\bar{r}_j+r_i,r_j+\bar{r}_i}^0 
C*                             ( \alpha+\bar{\alpha},\alpha+\bar{\alpha} )
C*                  IA pointing on 
C*                      B_{\bar{s}_j+s_i}^0 ( \beta - \bar{\beta} )
C*                  IB pointing on 
C*                      B_{s_j+\bar{s}_i}^0 ( \bar{\beta} - \beta ) 
C*
   10     AMA = A1 - A2
          APA = A1 + A2
          BMB = B1 - B2
          BPB = B1 + B2
          GO TO 7
C*
   11     K1 = K1 + 1
C*
          SMS = S(1) - S(2)
          RMR = R(1) - R(2)
          SPS = S(1) + S(2)
          RPR = R(1) + R(2)
C*
          IF = IF + P10 * RPR + P01 * ( RB(1) + RB(2) )
          IA = IA + SPS
          IB = IB + SB(1) + SB(2)
C*
          IF ((K1.EQ.1).OR.(K1.EQ.3)) THEN 
            ICALL = K1
            CALL CCI
            CALL MTRI
          ENDIF
C*
C*
C*           If it is a HETEROATOMIC molecule (non-equal nuclear charges),
C*           the heteroatomic type of wavefunction is used, i.e. 
C*           exp(\beta \eta_1 + \bar{\beta} \eta_2), and not cosh(...).
C*           in this case the current fragment of the matrix element is 
C*           ready:
C*
          IF (IZ1.NE.IZ2) GO TO 17
C*
C*           If it is a HOMOATOMIC molecule (equal nuclear charges), but
C*           \beta_1 = \beta_2 = 0.0 has been chosen for the basis set,
C*           the cosh(...) reduces also to the exp(...), and thus the 
C*           current fragment of the matrix element is ready:
C*
c          IF (COSH.EQ.0) GO TO 17
C*
          GO TO(13,14,13,15),K1
C*
C*               IB4 points to B_0^0 (0.0D+00)
   13     IA = IB4
          IB = IB4
          GO TO 16
C*
C*
C*               IB6 points to B_0^0 ( \bar{\beta} - \beta )
   14     IA = IB6
C*
C*               IB5 points to B_0^0 ( \beta - \bar{\beta} )
          IB = IB5
          GO TO 16
C*
C*
C*               IB5 points to B_0^0 ( \beta - \bar{\beta} )
   15     IA = IB5
C*
C*               IB6 points to B_0^0 ( \bar{\beta} - \beta )
          IB = IB6
C*
   16     IA = IA + SPS
          IB = IB + SB(1) + SB(2)
          J2 = J2 + K2
C*
          IF ((K1.EQ.1).OR.(K1.EQ.3)) THEN
            ICALL = K1
            CALL CCHI
            CALL MTRI
          ENDIF
C*
   17     GO TO(6,9,10,18),K1
   18   CONTINUE
   19 CONTINUE
C*
      RETURN
      END
C*
C*===================================================================
C*
      SUBROUTINE QFQ(F,S,NR,IOR,HULP)
C*
C*-------------------------------------------------------------------
C*
C     MODIFIED FOR ONE DIMENSIONAL S. 9/6/85
CQP      2
C     IMPLICIT REAL*16(A-H,O-Z)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-H,O-Z)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      DIMENSION HULP(1)
      DIMENSION F(NR,NR),S(1)
C*
C*------------------------------------------------------------------
C*
      DO 801 I=1,IOR
      II=IOR-I+1
      DO 202 K=1,II
      SUM=0D0
      DO 204 L=1,II
      FKL=F(K,L)
      IF(K.LT.L) FKL=F(L,K)
  204 SUM=SUM+FKL*S((II*(II-1))/2+L)
  202 HULP(K)=SUM
      DO 802 J=1,II
      SUM=0D0
      DO 206 K=1,J
  206 SUM=SUM+S((J*(J-1))/2+K)*HULP(K)
      F(II,J)=SUM
  802 F(J,II)=SUM
  801 CONTINUE
      RETURN
      END
C*
C*=================================================================
C*
      SUBROUTINE RELOC(IM,IV,N,S)
C*
C*-----------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(F)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(F)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      INTEGER S
      COMMON F(1)
C*
C*----------------------------------------------------------------
C*
      J=IM+(S*(S-1))/2
      I=J
      K1=S
      K2=IV+N-1
      K3=0
      DO 7 K=IV,K2
      K3=K3+1
      F(K)=F(J)
      IF(K3-S)5,6,6
    5 J=J+1
      GO TO 7
    6 J=J+K1
      K1=K1+1
    7 CONTINUE
      K=IV+S-1
      F(K2+1)=F(K)
      DO 8 J=K,K2
    8 F(J)=F(J+1)
      DO 12 K=S,N
      K1=K-1
      DO 10 L=1,K1
      J=I+K
      F(I)=F(J)
   10 I=I+1
   12 CONTINUE
      J=IM+(N*(N-1))/2
      DO 13 K=IV,K2
      F(J)=F(K)
   13 J=J+1
      RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE REV(Q,NR,N)
C*
C*------------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(A-H,O-Z)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-H,O-Z)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      DIMENSION Q(NR,NR)
C*
C*------------------------------------------------------------------
C*
      DO 1 J=1,N
      DO 1 I=J,N
      IF(I-J) 10,11,10
   11 II=1
      GO TO 12
   10 II=0
   12 S=0.0D0
      KK=I-1
      IF(KK.LT.J) GO TO 200
      DO 2 K=J,KK
      S=S-Q(I,K)*Q(K,J)
    2 CONTINUE
  200 CONTINUE
      Q(I,J)=(II+S)/Q(I,I)
    1 CONTINUE
      RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE SCPR(VEC,VE,N)
C*
C*------------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(X,V)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(X,V)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON /AREA11/X
      DIMENSION VEC(N),VE(N)
C*
C*------------------------------------------------------------------
C*
      X=0.
      DO 1 I=1,N
    1 X=X+VEC(I)*VE(I)
      RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE SEQU(MAT,C,N,N1)
C*
C*------------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(C,M,X)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(C,M,X)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      DIMENSION MAT(N1),C(N)
C*
C*-----------------------------------------------------------------
C*
      K=N
      J=N1
    3 J2=J
      J1=J
      I=K
    4 J1=J1-I
      I=I-1
      J2=J2-1
      X=MAT(J2)/MAT(J)
      C(I)=C(I)-X*C(K)
      L1=J1
      L2=J2
      DO 5 L=1,I
      MAT(L1)=MAT(L1)-X*MAT(L2)
      L1=L1-1
    5 L2=L2-1
      IF(I-1)6,6,4
    6 J=J-K
      K=K-1
      IF(K-1)7,7,3
    7 C(1)=C(1)/MAT(J)
      L=J
      L2=N-1
      DO 9 K=1,L2
      X=0
      DO 8 I=1,K
      L1=L+I
    8 X=X+C(I)*MAT(L1)
      L=L+K+1
    9 C(K+1)=(C(K+1)-X)/MAT(L)
      RETURN
      END
C*
C*=================================================================
C*
      SUBROUTINE SM2
C*
C*-----------------------------------------------------------------
C*
CQP      2
C     IMPLICIT INTEGER(T),REAL*16(X-Z)
C     GENERIC
CDP      1
      IMPLICIT INTEGER(T),REAL*8(X-Z)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON /AREA1R/Y
      COMMON /AREA5R/X,Z
      COMMON /AREA3/TEST,TEST1,TEST2
C*
C*-----------------------------------------------------------------
C*
      IF(X.EQ.0.0) GO TO 18
      GO TO(1,2,3,4,5),TEST2
C*
      ENTRY KM2
C*    =========
      IF(X.EQ.0.0) GO TO 18
      GO TO(1,12,13,14,15),TEST2
C*
      ENTRY SM1
C*    =========
      IF(X.EQ.0.0) GO TO 18
      GO TO(1,2,3,4,5,20,22,24),TEST1
C*
      ENTRY KM1
C*    =========
      IF(X.EQ.0.0) GO TO 18
      GO TO(1,12,13,14,15,21,23,25),TEST1
C*
      ENTRY SZE
C*    =========
      IF(X.EQ.0.0) GO TO 18
      GO TO(1,2,3,4,5,20,22,24),TEST
C*
      ENTRY KZE
C*    =========
      IF(X.EQ.0.0) GO TO 18
      GO TO(1,12,13,14,15,21,23,25),TEST
    1 GO TO 18
C*
    2 CALL JZW
      GO TO 17
    3 CALL JA
      GO TO 17
    4 CALL JA2B
      GO TO 17
    5 CALL JA
      GO TO 16
   12 CALL KZW
      GO TO 17
   13 CALL KA
      GO TO 17
   14 CALL KA2B
      GO TO 17
   20 CALL JA3B
      GO TO 17
   21 CALL KA3B
      GO TO 17
   22 CALL JAS1
      GO TO 17
   23 CALL KAS1
      GO TO 17
   24 CALL JA
      Y=3*Y
      GO TO 17
   25 CALL KA
      Y=3*Y
      GO TO 17
   15 CALL KA
   16 Y=2*Y
   17 X=X*Y
      Z=Z+X
   18 RETURN
      END
C*
C*=====================================================================
C*
      SUBROUTINE SMAT
C*
C*---------------------------------------------------------------------
C*
CQP      3
C     IMPLICIT REAL*16(A-F,X-Z),INTEGER(P,T)
C     GENERIC
C     REAL*16 IZ1,IZ2,HCN,HCN2
CDP      2
      IMPLICIT REAL*8(A-F,X-Z),INTEGER(P,T)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      REAL*8 IZ1,IZ2,HCN,HCN2
      COMMON F(1)
      COMMON /AREA2/IA,IB,I,P10,P20,P1M1,P01,P02,IF,IT(8),ISM,
     CIVM,IKM,IVE
      COMMON /AREA3/ITE(4),MMN
      COMMON /AREA4R/C(24)
      COMMON /AREA5R/X,Y
      COMMON /AREA5I/M,TES1,TES2,J1
      COMMON /AREA6R/HCN,HCN2
      COMMON /AREA6I/LTES,L,IDEL
      COMMON /AREA7R/AR,IZ1,IZ2
      COMMON /CALLS/ICALL
C*
C*------------------------------------------------------------------
C*
C*
C*    Calculation of the overlap matrix elements:
C*    ===========================================
C*    (Note, that this part of the program is called more than once 
C*     for every matrix element!) 
C*
      LTES = 0
      L = TES1
      J = ISM + J1
      Y = F(J)
      IF (TES1.EQ.0) GO TO 1
      I = IF + IDEL
      GO TO 2
    1 I = IF
    2 X = C(6)
C*
      CALL SZE
C*
      IF(TES1.EQ.0)GO TO 10
C*
      CALL LTEST(LL)
C*
      IF(LL.EQ.1) GO TO 10
      GO TO 2
C*
C*
C*    Calculation of the potential-energy matrix elements:
C*    ====================================================
C*    (Note, that this part of the program is called more than once 
C*     for every matrix element!) 
C*
      ENTRY VMAT
C*    ==========
C*
      LTES=0
      K=IA
      K1=IB
      J=IVM+J1
C*
C***
c      IF ((J.EQ.IVM).OR.(J.EQ.IVM+1).OR.(J.EQ.IVM+2)) THEN
c        WRITE(*,*) ' Checkpt.1: J, F(J) = ',J,F(J)
c      ENDIF
C***
      Y=F(J)
C*
C*
C*    Calculation of the electron-electron repulsion matrix elements:
C*    ---------------------------------------------------------------
C*
      IF (TES1.EQ.0) GO TO 5
      L = 0
      L1 = 1
      I = IF
      GO TO 6
C*
    5 L = 1
      L1 = 0
      I = IF + IDEL
C***
C*    Check, whether this is really the calculation of the electron-
C*    repulsion integrals. If this is the case, setting the coefficient
C*    (X=C(7)=2.0) to zero should produce the pure nuclear-attraction
C*    integrals in the place where the potential-energy matrix is build.
C*    This has been checked to be really the case (at least for \mu=0,2)
C*    basis functions.
C*
c    6 X = C(7)
    6 X = - C(7)
c    6 X = 0.0D+00
C***
C*
      CALL SM1
C*
      IF (TES1.NE.0) GO TO 7
C*
      CALL LTEST(LL)
C*
      IF (LL.EQ.1) GO TO 7
      GO TO 6
C*
    7 CONTINUE
C*
C*
C*    Calculation of the nuclear-attraction matrix elements:
C*    ------------------------------------------------------
C*
C*      Calculation of that part of the nuclear-attraction 
C*      matrix elements that is multiplied by - 4 * (Z_1 + Z_2) = C(8), 
C*      i.e. the terms \xi_1 * (\xi_2^2 - \eta_2^2) and 
C*      \xi_2 * (\xi_1^2 - \eta_1^2):
C*
C***
c      IF ((J.EQ.IVM).OR.(J.EQ.IVM+1).OR.(J.EQ.IVM+2)) THEN
c        WRITE(*,*) ' Checkpt.2: J, Y = ',J,Y
c      ENDIF
C***
      L = L1
      IA = K
      IB = K1
      I = IF + P10
      IF (TES1.NE.0) I = I + IDEL
    4 CONTINUE
      IF ((ICALL.EQ.2).OR.(ICALL.EQ.3)) THEN
        X = - C(8)
      ELSE
        X = C(8)
      ENDIF
C*
      CALL KZE
C*
C*      Calculation of that part of the nuclear-attraction 
C*      matrix elements that is multiplied by 4 * (Z_1 - Z_2) = C(9), 
C*      i.e. the terms \eta_1 * (\xi_2^2 - \eta_2^2) and 
C*      \eta_2 * (\xi_1^2 - \eta_1^2):
C*
      IF (C(9).EQ.0.) GO TO 3
      IA = IA + 1
      I = I - P10
      IF ((ICALL.EQ.2).OR.(ICALL.EQ.3)) THEN
        X = - C(9)
      ELSE
        X = C(9)
      ENDIF
C*
      CALL KZE
C*
      IA = IA - 1
      I = I + P10
    3 CONTINUE
      IF (TES1.EQ.0) GO TO 10
C*
      CALL LTEST(LL)
C*
      IF(LL.EQ.1) GO TO 10
      GO TO 4
C*
C*
C*    Calculation of the kinetic-energy matrix elements:
C*    ==================================================
C*    (Note, that this part of the program is called more than once 
C*     for every matrix element!) 
C*
      ENTRY KIN
C*    =========
C*
      J=J1+IKM
      L=TES1
      LTES=0
      I=IF
      IF(TES1.NE.0)I=I+IDEL
    8 Y=0.
      X=C(1)
C*
      CALL SM2
C*
      IF(MMN.EQ.0)GO TO 9
      IA=IA+1
      X=C(2)
C*
      CALL SM2
C*
      IA=IA-1
      I=I+P10
      X=C(3)
C*
      CALL SM2
C*
      IA=IA+1
      IB=IB+2
      I=I-P10
      X=C(2)
C*
      CALL KM2
C*
      IA=IA-1
      X=C(5)
C*
      CALL KM2
C*
      I=I+P10
      X=-C(3)
C*
      CALL KM2
C*
      X=-C(3)
      I=I+P02
      IB=IB-2
C*
      CALL KM2
C*
      I=I-P10
      IA=IA+1
      X=C(2)
C*
      CALL KM2
C*
      IA=IA-1
      X=C(5)
C*
      CALL KM2
C*
      I=I-P01
      IA=IA+1
      IB=IB+1
      X=C(10)
C*
      CALL KM2
C*
      IA=IA-1
      I=I+P10
      X=C(11)
C*
      CALL KM2
C*
      IA=IA-1
      X=C(12)
C*
      CALL KM2
C*
      IA=IA+2
      I=I-P20
      X=C(13)
C*
      CALL KM2
C*
      I=I+P1M1
      IA=IA-1
      IB=IB-1
    9 X=C(22)
C*
      CALL KZE
C*
      IA=IA+1
      X=C(14)
C*
      CALL KZE
C*
      IA=IA+1
      X=C(15)
C*
      CALL KZE
C*
      IA=IA-3
      X=C(16)
C*
      CALL KZE
C*
      IA=IA-1
      X=C(17)
C*
      CALL KZE
C*
      IA=IA+2
      I=I+P10
      X=C(18)
C*
      CALL KZE
C*
      I=I+P10
      X=C(19)
C*
      CALL KZE
C*
      I=I-3*P10
      X=C(20)
C*
      CALL KZE
C*
      I=I-P10
      X=C(21)
C*
      CALL KZE
C*
      F(J)=F(J)+Y
      IF(TES1.EQ.0)GO TO 12
      I=I+P20
      X=Y
      Y=F(J)
C*
      CALL LTEST(LL)
C*
      IF(LL.EQ.1) GO TO 12
      GO TO 8
   10 CONTINUE
C***
c      IF ((J.EQ.IVM).OR.(J.EQ.IVM+1).OR.(J.EQ.IVM+2)) THEN
c        WRITE(*,*) ' Checkpt.3: j = ',j
c        WRITE(*,*) ' y, y-F(j): ',Y,Y-F(J)
c        WRITE(*,*) ' '
c      ENDIF
C***
      F(J)=Y
C*
   12 RETURN
      END
C*
C*===================================================================
C*
      SUBROUTINE SPOL(S,NR,N)
C*
C*-------------------------------------------------------------------
C*
CQP      2
C     IMPLICIT REAL*16(A-H,O-Z)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-H,O-Z)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      DIMENSION S(NR,NR)
C*
      CHARACTER*256 filnam
C*
C*------------------------------------------------------------------
C*
CQP      2
C     DABS(X)=ABS(X)
C     DSQRT(X)=SQRT(X)
C
  100 FORMAT(' MATRIX DECOMPOSED IS NOT POSITIVE DEFINITE')
  101 FORMAT(D23.15)
      DET=1D0
      DO 1 K=1,N
      W=S(K,K)
      KK=K-1
      IF(KK.LT.1) GO TO 200
      DO 2 J=1,KK
      W=W-S(K,J)**2
    2 CONTINUE
  200 CONTINUE
      IF(W) 11,11,10
   11 WRITE(6,100)
C*
C*ASB 2  07.02.1995
      CALL getenv('H2ERR',filnam)      
      OPEN(UNIT=99, FILE=filnam, STATUS='NEW')
      CLOSE(99)
C***
      RETURN
C*
   10 S(K,K)=DSQRT(DABS(W))
C     DET=DET*W
      KL=K+1
      IF(N.LT.KL) GO TO 1
      DO 20 L=KL,N
      W=S(L,K)
      KU=K-1
      IF(KU.LT.1) GO TO 300
      DO 30 J=1,KU
      W=W-S(K,J)*S(L,J)
   30 CONTINUE
  300 CONTINUE
      S(L,K)=W/S(K,K)
   20 CONTINUE
    1 CONTINUE
  103 FORMAT(' DETERMINANT=',D23.15)
      RETURN
      END
C*
C*=================================================================
C*
      SUBROUTINE STODC
C*
C*-----------------------------------------------------------------
C*
      IMPLICIT INTEGER(T)
CQP      1
C     GENERIC
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      COMMON /AREA3/TEST,TEST1,TEST2,MPN
      COMMON /AREA5I/M
C*
C*----------------------------------------------------------------
C*
      MPNP1=MPN+1
      M1=M+1
      GO TO(1,2,3,31),M1
    1 GO TO(4,5,6,7,8,15,16),MPNP1
    2 GO TO(9,10,11,12,13,17,18),MPNP1
    3 GO TO(9,9,9,10,11,19,20),MPNP1
   31 GO TO(9,9,9,9,9,10,11),MPNP1
    4 TEST2=1
      TEST1=2
      TEST=2
      GO TO 14
    5 TEST2=2
      TEST1=2
      TEST=3
      GO TO 14
    6 TEST2=2
      TEST1=3
      TEST=3
      GO TO 14
    7 TEST2=3
      TEST1=3
      TEST=4
      GO TO 14
    8 TEST2=3
      TEST1=4
      TEST=4
      GO TO 14
    9 TEST2=1
      TEST1=1
      TEST=1
      GO TO 14
   10 TEST2=1
      TEST1=1
      TEST=2
      GO TO 14
   11 TEST2=1
      TEST1=2
      TEST=1
      GO TO 14
   12 TEST2=2
      TEST1=1
      TEST=5
      GO TO 14
   13 TEST2=1
      TEST1=5
      TEST=1
      GO TO 14
   15 TEST2=4
      TEST1=4
      TEST=6
      GO TO 14
   16 TEST2=4
      TEST1=6
      TEST=6
      GO TO 14
   17 TEST2=5
      TEST1=1
      TEST=7
      GO TO 14
   18 TEST2=1
      TEST1=7
      TEST=1
      GO TO 14
   19 TEST2=2
      TEST1=1
      TEST=8
      GO TO 14
   20 TEST2=1
      TEST1=8
      TEST=1
   14 RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C*
C*------------------------------------------------------------------
C*
CQP      1
C     GENERIC
C
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      INTEGER I,J,K,L,M,N,L1,NM,MML,IERR,II
CQP      4
C     REAL*16 D(N),E(N),Z(NM,N)
C     REAL*16 B,C,F,G,H,P,R,S
C     REAL*16 MACHEP
C     REAL*16 DABS,DSQRT,DSIGN
CDP      3
      REAL*8 D(N),E(N),Z(NM,N)
      REAL*8 B,C,F,G,H,P,R,S
      REAL*8 MACHEP
C
C      THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C      A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD
C      THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C      BE FOUND IF TRED2 HAS BEEN USED TO REDUCED THIS
C      FULL MATRIX TO TRIDIAGONAL FORM.
C
C
C     ON INPUT
C
C     NM MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C     ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C
C     N IS THE ORDER OF THE MATIRX
C
C     D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX
C
C     E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C      IN ITS LAST N-1 POSITIONS. E(1) IS ARBITRARY
C
C     Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C       REDUCTION BY TRED2, IF PERFORMED. IF THE EIGENVECTORS
C       OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C     THE IDENTITY MATRIX
C
C     ON OUTPUT
C
C     D CONTAINS THE EIGENVALUES IN ASCENDING ORDER. IF AN
C       ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C     UNORDERED FOR INDICES 1,2,....,IERR-1.
C
C     E IS DESTROYED
C
C     Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C     TRIDIAGONAL (OR FULL) MATRIX. IF AN ERROR EXIT IS MADE,
C     Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C     EIGENVALUES
C
C     IERR IS SET TO
C      ZERO    FOR NORMAL RETURN
C      J       IF THE J-TH EIGENVALUES HAS NOT BEEN
C              DETERMINED AFTER 30 ITERATIONS.
C
C----------------------------------------------------------------------
C
C                  MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                  THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                  MACHEP=16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                  ON S360
CVAX     1
      DATA  MACHEP/1.D-16/
CQP      3
C     DABS(X)=ABS(X)
C     DSQRT(X)=SQRT(X)
C     DSIGN(X,Y)=SIGN(X,Y)
CIBM     1
C     MACHEP=16.D0**(-13)
CQP      1
C     MACHEP=16D0**(-26)
C
      IERR=0
      IF(N.EQ.1) GO TO 1001
C
      DO 100 I=2,N
  100 E(I-1)=E(I)
C
      F=0.
      B=0.
      E(N)=0.
C
      DO 240 L=1,N
      J=0
      H=MACHEP*(DABS(D(L))+DABS(E(L)))
      IF(B.LT.H) B=H
C
C            LOOK FOR SMALL SUB-DIAGONAL ELEMENT
C
      DO 110 M=L,N
      IF(DABS(E(M)).LE.B) GO TO 120
C         E(N) -IS ALWAYS ZERO, SO THERE IS NO EXIT
C         THROUGH THE BOTTOM OF THE LOOP
  110 CONTINUE
C
  120 IF(M.EQ.L) GO TO 220
  130 IF(J.EQ.30) GO TO 1000
C                    FORM SHIFT
      J=J+1
      L1=L+1
      G=D(L)
      P=(D(L1)-G)/(2.*E(L))
      R=DSQRT(P*P+1.)
      D(L)=E(L)/(P+DSIGN(R,P))
      H=G-D(L)
C
      DO 140 I=L1,N
  140 D(I)=D(I)-H
      F=F+H
C                    QL TRANSFORMATION
      P=D(M)
      C=1.
      S=0.
      MML=M-L
      DO 200 II=1,MML
      I=M-II
      G=C*E(I)
      H=C*P
      IF(DABS(P).LT.DABS(E(I))) GO TO 150
      C=E(I)/P
      R=DSQRT(C*C+1.)
      E(I+1)=S*P*R
      S=C/R
      C=1./R
      GO TO 160
  150 C=P/E(I)
      R=DSQRT(C*C+1)
      E(I+1)=S*E(I)*R
      S=1./R
      C=C*S
  160 P=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
C                  FORM VECTOR
      DO 180 K=1,N
      H=Z(K,I+1)
      Z(K,I+1)=S*Z(K,I)+C*H
      Z(K,I)=C*Z(K,I)-S*H
  180 CONTINUE
C
  200 CONTINUE
C
      E(L)=S*P
      D(L)=C*P
      IF(DABS(E(L)).GT.B) GO TO 130
  220 D(L)=D(L)+F
  240 CONTINUE
C                  ORDER EIGENVALUES AND IEGENVECTORS
      DO 300 II=2,N
      I=II-1
      K=I
      P=D(I)
C
      DO 260 J=II,N
      IF(D(J).GE.P) GO TO 260
      K=J
      P=D(J)
  260 CONTINUE
C
      IF(K.EQ.I) GO TO 300
      D(K)=D(I)
      D(I)=P
C
      DO 280 J=1,N
      P=Z(J,I)
      Z(J,I)=Z(J,K)
      Z(J,K)=P
  280 CONTINUE
C
  300 CONTINUE
      GO TO 1001
C                    SET ERROR -- NO CONVERGENCE TO AN
C                    EIGENVALUE AFTER 30 ITERATIONS
 1000 IERR=L
 1001 RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE TRANS(F,P,NR,IOR)
C*
C*------------------------------------------------------------------
C*
C     MODIFIED FOR USE OF ONE DIMENSIONAL P. 9/6/85 KSZ.
CQP      2
C     IMPLICIT REAL*16(A-H,O-Z)
C     GENERIC
CDP      1
      IMPLICIT REAL*8(A-H,O-Z)
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      DIMENSION F(NR,NR),P(1)
      DO 550 K=1,IOR
      DO 550 I=1,IOR
      SS=0.0D0
      DO 551 J=I,IOR
      PP=P((J*(J-1))/2+I)
      IF(J.LT.I) PP=P((I*(I-1))/2+J)
      SS=SS+PP*F(J,K)
  551 CONTINUE
      F(I,K)=SS
  550 CONTINUE
      RETURN
      END
C*
C*==================================================================
C*
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
C*
C*------------------------------------------------------------------
C*
CQP      1
C     GENERIC
C*ASB 2   15.09.1993
      INCLUDE 'h2upar2.for'
C*
      INTEGER I,J,K,L,N,II,NM,JP1
CQP      3
C     REAL*16 Z(NM,N),D(N),E(N),A(NM,N)
C     REAL*16 F,G,H,HH,SCALE
C     REAL*16 X,Y,DABS,DSQRT,DSIGN
CDP      3
      REAL*8 Z(NM,N),D(N),E(N),A(NM,N)
      REAL*8 F,G,H,HH,SCALE
      REAL*8 X,Y,DABS,DSQRT,DSIGN
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
C     ORTHOGONAL SIMILARITY  TRANSFORMATIONS.
C
C       ON INPUT
C     NM MUST BE SET THE ROW DIMENSION OF TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C       DIMENSION STATEMENT
C
C     N IS THE ORDER OF THE MATRIX
C
C     Z CONTAINS THE REAL SYMMETRIC INPUT MATRIX. ONLY THE
C       LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C     D CONTAINS THE DIAGONAL ELEMENST OF THE TRIDIAGONAL MATRIX.
C
C     E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C      MATRIX IN ITS LAST N-1 POSITIONS. E(1) IS SET TO ZERO
C
C     Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
C       PRODUCED IN THE REDUCTION.
C
C---------------------------------------------------------------------
C
CQP       3
C     DSQRT(X)=SQRT(X)
C     DSIGN(X,Y)=SIGN(X,Y)
C     DABS(X)=ABS(X)
      DO 100 I=1,N
      DO 100 J=1,I
      Z(I,J)=A(I,J)
  100 CONTINUE
      IF(N.EQ.1) GO TO 320
      DO 300 II=2,N
C     IF(II/10*10.EQ.II) WRITE(6,5) II
    5 FORMAT(10X,'II=  ',I6)
         I=N+2-II
        L=I-1
        H=0.
        SCALE=0.
        IF(L.LT.2) GO TO 130
C               SCALE ROW
        DO 120 K=1,L
  120   SCALE=SCALE+DABS(Z(I,K))
C
      IF(SCALE.NE.0.) GO TO 140
  130 E(I)=Z(I,L)
      GO TO 290
C
  140 DO 150 K=1,L
        Z(I,K)=Z(I,K)/SCALE
        H=H+Z(I,K)*Z(I,K)
  150 CONTINUE
C
      F=Z(I,L)
      G=-DSIGN(DSQRT(H),F)
      E(I)=SCALE*G
      H=H-F*G
      Z(I,L)=F-G
      F=0.
C
      DO 240 J=1,L
      Z(J,I)=Z(I,J)/H
      G=0.
C       FORM ELEMENT OF A*U
      DO 180 K=1,J
  180 G=G+Z(J,K)*Z(I,K)
C
      JP1=J+1
      IF(L.LT.JP1) GO TO 220
C
      DO 200 K=JP1,L
  200 G=G+Z(K,J)*Z(I,K)
C        FORM ELEMENT OF P
  220 E(J)=G/H
      F=F+E(J)*Z(I,J)
  240 CONTINUE
C
      HH=F/(H+H)
C        FORM REDUCED A
      DO 260 J=1,L
      F=Z(I,J)
      G=E(J)-HH*F
      E(J)=G
C
      DO 260 K=1,J
      Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  260 CONTINUE
C
  290 D(I)=H
  300 CONTINUE
C
  320 D(1)=0.
      E(1)=0.
C        ACCUMULATION OF TRANSFORMATION MATRICES
      DO 500 I=1,N
      L=I-1
      IF(D(I).EQ.0.) GO TO 380
C
      DO 360 J=1,L
      G=0.
C
      DO 340 K=1,L
  340 G=G+Z(I,K)*Z(K,J)
C
      DO 360 K=1,L
      Z(K,J)=Z(K,J)-G*Z(K,I)
  360 CONTINUE
C
  380 D(I)=Z(I,I)
      Z(I,I)=1.
      IF(L.LT.1) GO TO 500
C
      DO 400 J=1,L
      Z(I,J)=0.
      Z(J,I)=0.
  400 CONTINUE
C
  500 CONTINUE
C
      RETURN
      END
