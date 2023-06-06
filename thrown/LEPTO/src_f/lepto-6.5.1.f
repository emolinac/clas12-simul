C######################################################################C
C                                                                      C
C                          L E P T O                                   C
C                                                                      C
C                  A Monte Carlo Generator for                         C
C            Deep Inelastic Lepton-Nucleon Scattering                  C
C                                                                      C
C                  version 6.5.1,  October 31, 1996                    C
C                                                                      C
C   Authors:                                                           C
C   Gunnar Ingelman, DESY theory group     Dept. of Radiation Sciences C
C                    (room 202  bldg 2a)   Uppsala University          C
C                    D-22603 Hamburg, FRG  S-751 21 Uppsala, Sweden    C
C             phone: +49(40)8998-2795      +46(18)18-3884              C
C               fax:            -2777               -3833              C
C            e-mail: ingelman@desy.de      (ingelman@tsl.uu.se)        C
C   Anders Edin, Dept. of Radiation Sciences, edin@tsl.uu.se           C
C   Johan Rathsman, Dept. of Radiation Sciences, rathsman@tsl.uu.se    C
C                                                                      C
C   Contributions on parton cascades: M. Bengtsson, T. Sjostrand       C
C                                                                      C
C   Home page: http://www3.tsl.uu.se/thep/lepto with code and manual   C
C                                                                      C
C   Manual: G.Ingelman, A.Edin, J.Rathsman, DESY 96-057                C
C                                                                      C
C   Please report any problems or suggestions for improvements.        C
C                                                                      C
C######################################################################C

      SUBROUTINE LTIMEX(TIME)
C...Interface routine to transfer a call to some machine-dependent
C...routine to get the execution time used since job started.
C...Nice, but not necessary information. Can also be called by user.

      TIME=0.
C...Use of CERN library routine Z007, replace/delete if not available.
      CALL TIMEX(TIME)
      RETURN
      END

C **********************************************************************

      BLOCK DATA LEPTOD

C...Give sensible default values to switches and parameters.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LFLMIX/ CABIBO(4,4)
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON /LGRID/ NXX,NWW,XX(31),WW(21),PQG(31,21,3),PQQB(31,21,2),
     &QGMAX(31,21,3),QQBMAX(31,21,2),YCUT(31,21),XTOT(31,21),NP
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     &FLMT(41,16)
      COMMON /LYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4),
     &MAXFIN,RELUP,RELERR,RELER2,FCNMAX
      COMMON /LMINUC/ NAMKIN(4),NAM(30)
      CHARACTER*10 NAMKIN,NAM

C...LEPTOU: Cuts, basic switches and parameters.
      DATA CUT/1.E-04,1.,0.,1.,4.,1.E+08,5.,1.E+08,1.,1.E+08,1.,1.E+08,
     &0.,3.1416/
C...            0    1    2    3    4    5    6    7    8    9
      DATA LST/      0,   1,   5,   1,   3,   1,   1,  12,   5,
     1          1,   0,   4,   5,   4,   9,   1,   0,   2, -10,
     2          5,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     3          0,   0,   0,   0,   2,   1,   5*0/
C...            0    1    2    3    4    5    6    7    8    9
      DATA PARL/     1.,  1.,  0.44,0.75,.2319,0., 0.5, .04, 4.,
     1         0.44, 0.01,0.01,0.1, 0.35,0.01,7.29735E-03,
     &                                   1.16639E-05,0.044,0.03,
     2         0.1,10*0./
C...Internally used variables.
      DATA PARI/40*0./
      DATA QC/-.33333,.66667,-.33333,.66667,-.33333,.66667,
     &        -.33333,.66667/
      DATA CABIBO/.95,.05,2*0.,.05,.948,.002,2*0.,.002,.998,4*0.,1./
      DATA OPTX/1.,3*0./,OPTY/1.,3*0./,OPTQ2/1.,3*0./,OPTW2/1.,3*0./
      DATA NXX,NWW/31,21/
      DATA PQG,PQQB,QGMAX,QQBMAX/6510*0./,YCUT/651*0./,XTOT/651*0./
      DATA NFX,NFQ/41,16/,FLQT,FLGT,FLMT/1968*0./
      DATA XKIN/1.,2.,3.,4./,UKIN,WKIN,AIN,BIN/16*0./,MAXFIN/2000/
      DATA RELUP,RELERR,RELER2/0.1,0.05,0.05/
      DATA NAMKIN/'         x','          ','          ','          '/
      DATA IPY/
     1 0,     0,     2,     2,     6,     1,     1,     6,     3,     1,
     2 3,     1,     1,     2,     1,     1,     4,     1,     1,     1,
     3 0,     1,     1,     1,     1,     1,     1,     0,     0,     0,
     4 1,     2,     1,     1,    30,    33,     1,     1,     7,     0,
     5 0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     6 0,     0,     0,     1,   100,     0,     0,     0,     0,     0,
     7 0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     8 0,     0,     0,     0,     0,     0,     0,     0,     0,     0/
      DATA (PYPAR(I),I=1,40)/
     1   7.299E-03,   2.290E-01,   2.000E-01,   2.500E-01,   4.000E+00,
     1   1.000E+00,   4.400E-01,   4.400E-01,   7.500E-02,   0.000E+00,
     2   2.000E+00,   2.000E+00,   1.000E+00,   0.000E+00,   3.000E+00,
     2   1.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   1.000E+00,
     3   2.500E-01,   1.000E+00,   2.000E+00,   1.000E-03,   1.000E+00,
     3   1.000E+00,   1.000E+00,  -2.000E-02,  -1.000E-02,   0.000E+00,
     4   0.000E+00,   1.600E+00,   0.500E+00,   0.200E+00,   3.894E-01,
     4   1.000E+00,   3.300E-01,   6.600E-01,   0.000E+00,   1.000E+00/
      DATA (PYPAR(I),I=41,80)/
     5   2.260E+00,   1.000E+04,   1.000E-04,   0.000E+00,   0.000E+00,
     5   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     6   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     6   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     7   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     7   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     8   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     8   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00/
      DATA PYVAR/80*0./
      END

C **********************************************************************


      SUBROUTINE LINIT(LFILE,LEPIN,PLZ,PPZ,INTER)

C...Initialize for an incoming lepton (type LEPIN, momentum pz=PLZ)
C...and target nucleon (momentum pz=PPZ) to interact via INTER.
C...Find maximum of differential cross section, calculate QCD event
C...probabilities or read them from logical file LFILE (if >0).
C...Numerical integration to obtain total cross-section.

      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LGRID/ NXX,NWW,XX(31),WW(21),PQG(31,21,3),PQQB(31,21,2),
     &QGMAX(31,21,3),QQBMAX(31,21,2),YCUT(31,21),XTOT(31,21),NP
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4),
     &MAXFIN,RELUP,RELERR,RELER2,FCNMAX
      COMMON /LMINUC/ NAMKIN(4),NAM(30)
      COMMON /LPFLAG/ LST3
      COMMON /LYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARIPY(200)
      CHARACTER*10 NAMKIN,NAM
      DIMENSION LSTW(40),PARLW(30)
      DOUBLE PRECISION DTHETA,DPHI,DBETA
      DATA PI/3.1415927/,NCALL/0/
      
      DOUBLE PRECISION Target
      DOUBLE PRECISION KFermi

      NCALL=NCALL+1
      LST3=LST(3)
      IF(LST(18).GE.1) THEN
C...W, Z masses from theta-Weinberg, Fermi constant GF and rad. corr.
        PMAS(24,1)=SQRT(PI*PARL(16)/(SQRT(2.)*PARL(17)*PARL(5)*
     &  (1.-PARL(18))))
        PMAS(23,1)=PMAS(24,1)/SQRT(1.-PARL(5))
      ENDIF
C...Couplings between Z0 and left/right-handed leptons and quarks.
      ZL(1,1)=-.5+PARL(5)
      ZL(1,2)=PARL(5)
      ZL(2,1)=ZL(1,2)
      ZL(2,2)=ZL(1,1)
      ZL(1,3)=0.5
      ZL(2,3)=0.
      ZL(1,4)=0.
      ZL(2,4)=0.5
      DO 10 IFL=1,8
      ZQ(1,IFL)=SIGN(0.5,QC(IFL))-QC(IFL)*PARL(5)
   10 ZQ(2,IFL)=-QC(IFL)*PARL(5)

C...Set initial state.
      LST(23)=INTER
      KSAVE(1)=LEPIN
      KSAVE(2)=2212
      K(1,1)=21
      K(1,2)=KSAVE(1)
      K(1,3)=0
      K(1,4)=0
      K(1,5)=0
      K(2,1)=21
      K(2,2)=KSAVE(2)
      K(2,3)=0
      K(2,4)=0
      K(2,5)=0
      P(1,1)=0.
      P(1,2)=0.
      P(1,3)=PLZ
      P(1,5)=ULMASS(KSAVE(1))
      P(1,4)=SQRT(P(1,3)**2+P(1,5)**2)

      P(2,1)=0.
      P(2,2)=0.
      P(2,3)=PPZ
C      write (*,*)'After c++ call, fortran last value of P(2,3) = ',P(2,3)
      P(2,5)=ULMASS(KSAVE(2))
      P(2,4)=SQRT(P(2,3)**2+P(2,5)**2)
      N=2
      LST(28)=3
C...Save momentum vectors of incoming particles
      DO 20 I=1,2
      DO 20 J=1,5
   20 PSAVE(3,I,J)=P(I,J)
C...Dot-product of initial particles, cms energy
      PARL(21)=2.*(DBLE(P(1,4))*DBLE(P(2,4))-DBLE(P(1,3))*DBLE(P(2,3)))
      ROOTS=SQRT((DBLE(P(1,4))+DBLE(P(2,4)))**2
     &          -(DBLE(P(1,3))+DBLE(P(2,3)))**2)
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,1000)
     &LEPIN,(P(1,J),J=1,3),PARL(1),PARL(2),(P(2,J),J=1,3),INTER,ROOTS
      IF(PLZ*PPZ.GT.0.1) THEN 
        WRITE(6,1010)
        STOP
      ENDIF



C...Reduced header for Jetset/Pythia
      MSTU(12)=0
      MSTP(122)=0
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) 
     &WRITE(6,1020) MSTU(181),MSTU(182),MSTP(181),MSTP(182)
C...If JETSET version before 7.402, problem with azimuthal dependence
C...in LUSHOW solved by chosing flat azimuthal dependence.
      IF(MSTU(181).LE.7.AND.MSTU(182).LT.402) THEN
        MSTJ(46)=0
        WRITE(6,1030) MSTJ(46)
      ENDIF
C...Initialize PYTHIA for parton densities.
      IF(LST(15).GT.0) THEN
C...Set switches and parameters for parton densities in PYSTFU.
        MSTP(51)=LST(15)
        MSTP(52)=LST(16)
        MSTP(58)=LST(12)
      ENDIF
      CALL PYINIT('NONE','e-','p',ROOTS)
      PARL(26)=PARP(1)
CAE--     use Lambda from parton densities in initial cascade
      PYPAR(21)=PARP(1)
C...Reset PYTHIA 4.8 parameters from LEPTO parameters.
      IF(MOD(LST(8),10).EQ.3.OR.MOD(LST(8),10).EQ.5) IPY(13)=0
      IF(LST(35).EQ.0.AND.
     &(MOD(LST(8),10).EQ.4.OR.MOD(LST(8),10).EQ.5)) IPY(14)=0
      IPY(8)=LST(12)

      IF(PSAVE(3,1,3).LT.0.) THEN
C...Flip event to have initial lepton along +z axis
        P(1,3)=-P(1,3)
        P(2,3)=-P(2,3)
      ENDIF
C...Boost parameters to cms of incoming particles
      DBETA(1,1)=0.D0
      DBETA(1,2)=0.D0
      DBETA(1,3)=(DBLE(P(1,3))+DBLE(P(2,3)))/(DBLE(P(1,4))+DBLE(P(2,4)))
      SPHI(1)=0.D0
      STHETA(1)=0.D0
      IF(LST(17).NE.0) THEN
C...For varying beam energies, transform to cms, lepton along +z axis.
        CALL LUDBRB(0,0,0.,0.,0.D0,0.D0,-DBETA(1,3))
        SPHI(1)=ULANGL(P(1,1),P(1,2))
        CALL LUDBRB(0,0,0.,-SPHI(1),0.D0,0.D0,0.D0)
        STHETA(1)=ULANGL(P(1,3),P(1,1))
        CALL LUDBRB(0,0,-STHETA(1),0.,0.D0,0.D0,0.D0)
        LST(28)=2
      ENDIF

C...Effective limits on kinematic variables x, y, Q**2, W**2
      PM2=P(2,5)**2
      S=PARL(21)
      XMIN=MAX(CUT(1),0.)
      XMAX=MIN(CUT(2),1.)
      YMIN=MAX(CUT(3),0.)
      YMAX=MIN(CUT(4),1.)
      Q2MIN=MAX(CUT(5),0.)
      Q2MAX=MIN(CUT(6),S)
      W2MIN=MAX(CUT(7),0.)
      W2MAX=MIN(CUT(8),S+PM2)
      UMIN=MAX(CUT(9),0.)
      UMAX=MIN(CUT(10),S/(2.*P(2,5)))
      DO 40 I=1,2
      XMIN=MAX(XMIN,Q2MIN/(S*YMAX),Q2MIN/(2.*P(2,5)*UMAX),
     &1.-(W2MAX-PM2)/MAX(S*YMIN,1.E-22),
     &1.-(W2MAX-PM2)/MAX(2.*P(2,5)*UMIN,1.E-22))
      XMAX=MIN(XMAX,Q2MAX/MAX(S*YMIN,1.E-22),
     &Q2MAX/MAX(2.*P(2,5)*UMIN,1.E-22),
     &1.-(W2MIN-PM2)/(S*YMAX),1.-(W2MIN-PM2)/(2.*P(2,5)*UMAX))
      YMIN=MAX(YMIN,Q2MIN/(S*XMAX),(W2MIN-PM2)/(S*(1.-XMIN)),
     &(W2MIN-PM2+Q2MIN)/S,2.*P(2,5)*UMIN/S)
      YMAX=MIN(YMAX,Q2MAX/MAX(S*XMIN,1.E-22),
     &(W2MAX-PM2)/MAX(S*(1.-XMAX),1.E-22),
     &(W2MAX-PM2+Q2MAX)/S,2.*P(2,5)*UMAX/S)
      Q2MIN=MAX(Q2MIN,S*XMIN*YMIN,S*YMIN-W2MAX+PM2,
     &2.*P(2,5)*UMIN*XMIN,(W2MIN-PM2)*XMIN/(1.-XMIN))
      Q2MAX=MIN(Q2MAX,S*XMAX*YMAX,S*YMAX-W2MIN+PM2,
     &2.*P(2,5)*UMAX*XMAX,(W2MAX-PM2)*XMAX/MAX(1.-XMAX,1.E-22))
      W2MIN=MAX(W2MIN,S*(1.-XMAX)*YMIN+PM2,Q2MIN*(1.-XMAX)/XMAX+PM2,
     &S*YMIN-Q2MAX+PM2,2.*P(2,5)*UMIN*(1.-XMAX)+PM2)
      W2MAX=MIN(W2MAX,S*(1.-XMIN)*YMAX+PM2,
     &Q2MAX*(1.-XMIN)/MAX(XMIN,1.E-22)+PM2,
     &S*YMAX-Q2MIN+PM2,2.*P(2,5)*UMAX*(1.-XMIN)+PM2)
C     UMIN=MAX(UMIN,....)
C     UMAX=MIN(UMAX,....)
   40 CONTINUE
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,1050)
     &CUT,XMIN,XMAX,YMIN,YMAX,Q2MIN,Q2MAX,W2MIN,W2MAX,UMIN,UMAX
      IF(XMAX.LT.XMIN.OR.YMAX.LT.YMIN.OR.Q2MAX.LT.Q2MIN.OR.
     &W2MAX.LT.W2MIN) THEN
        IF(LST(3).GE.1) WRITE(6,1100)
        IF(LST(3).GE.2) THEN
          WRITE(6,1900)
          STOP
        ENDIF
      ENDIF
      IF(XMIN.LT.1.E-10.OR.Q2MIN.LT.1.E-01) THEN
        IF(LST(3).GE.1) WRITE(6,1110)
        IF(LST(3).GE.2) THEN
          WRITE(6,1900)
          STOP
        ENDIF
      ENDIF

      PARI(11)=(PARL(1)-PARL(2))/PARL(1)
      KSAVE(4)=LEPIN
      ILEP=1
      IF(LEPIN.LT.0) ILEP=2
      INU=0
      IF(IABS(LEPIN).EQ.12.OR.IABS(LEPIN).EQ.14
     &.OR.IABS(LEPIN).EQ.16) INU=1
      IF(INU.EQ.1) THEN
C...Set full polarisation for incoming neutrino.
        PARL(6)=-1.
        IF(LEPIN.LT.0) PARL(6)=1.
      ENDIF
      IF(LST(23).EQ.1.AND.INU.EQ.0) THEN
C...Electromagnetic interaction.
        KSAVE(3)=22
        IG=1
        IZ=0
      ELSEIF(LST(23).EQ.2) THEN
C...Weak charged current, only one helicity state contributes.
        IF(KSAVE(1).LT.0.AND.PARL(6).LT.-0.99
     &  .OR.KSAVE(1).GT.0.AND.PARL(6).GT.0.99) THEN
          IF(LST(3).GE.1) WRITE(6,1150) LEPIN,PARL(6)
          IF(LST(3).GE.2) THEN
            WRITE(6,1900)
            STOP
          ENDIF
        ENDIF
        IF(MOD(IABS(LEPIN),2).EQ.0) THEN
          KSAVE(3)=ISIGN(24,LEPIN)
          KSAVE(4)=ISIGN(IABS(LEPIN)-1,LEPIN)
        ELSE
          KSAVE(3)=ISIGN(24,-LEPIN)
          KSAVE(4)=ISIGN(IABS(LEPIN)+1,LEPIN)
        ENDIF
      ELSEIF(LST(23).EQ.3.OR.(LST(23).EQ.4.AND.INU.EQ.1)) THEN
C...Weak neutral current.
        KSAVE(3)=23
        IG=0
        IZ=1
      ELSEIF(LST(23).EQ.4.AND.INU.EQ.0) THEN
C...Neutral current, electromagnetic and weak with interference.
        KSAVE(3)=23
        IG=1
        IZ=1
      ELSE
        IF(LST(3).GE.1) WRITE(6,1200) INTER,LEPIN
        IF(LST(3).GE.2) THEN
          WRITE(6,1900)
          STOP
        ENDIF
      ENDIF

C...Choice of independent variables.
      IF(LST(1).EQ.0) THEN
        LST(31)=1
        IF(INTER.EQ.2.OR.INTER.EQ.3) LST(31)=2
      ELSE
        LST(31)=IABS(LST(1))
      ENDIF
      IF(LST(31).LT.1.OR.LST(31).GT.3) THEN
        IF(LST(3).GE.1) WRITE(6,1210) LST(1),LST(31)
        IF(LST(3).GE.2) THEN
          WRITE(6,1900)
          STOP
        ENDIF
      ENDIF
      IF(LST(1).LT.0) THEN
C...User-defined optimization parameters.
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &  WRITE(6,1220) OPTX,OPTY,OPTQ2,OPTW2
      ELSE
C...Set optimization parameters.
        DO 50 I=1,4
        OPTX(I)=0.
        OPTY(I)=0.
        OPTQ2(I)=0.
   50   OPTW2(I)=0.
        IF(INTER.EQ.1) THEN
          OPTX(2)=1.
          OPTY(1)=1.
          OPTQ2(3)=1.
          OPTW2(3)=1.
        ELSEIF(INTER.EQ.4) THEN
          OPTX(1)=0.1
          OPTX(2)=1.
          OPTY(1)=1.
          OPTQ2(1)=0.5
          OPTQ2(2)=0.5
          OPTQ2(3)=1.
          OPTW2(1)=0.5
          OPTW2(2)=0.5
          OPTW2(3)=1.
        ELSE
          OPTX(1)=1.
          OPTY(1)=1.
          OPTQ2(1)=1.
          OPTW2(1)=1.
        ENDIF
      ENDIF

C...Initialize Monte Carlo estimate of cross section.
      PARL(24)=0.
      PARI(27)=0.
      PARI(28)=0.
      PARI(29)=0.
      PARI(30)=0.
      PARI(32)=0.
      IF(LST(23).EQ.2) THEN
C...Constant factor GF**2/pi for CC, transformation to picobarn.
        PARI(31)=PARL(17)**2/PI*0.39E+09
      ELSE
C...Constant factor 2*pi*alpha**2 for NC, transformation to picobarn.
        PARI(31)=2.*PI*PARL(16)**2*0.39E+09
      ENDIF
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &WRITE(6,1250) (I,LST(I),LST(I+10),PARL(I),PARL(I+10),I=1,10)

C...Set up grid with longitudinal structure function, QCD & target mass;
C...only when photon exchange is included
      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      IF(LST(11).NE.0.AND.(INTER.EQ.1.OR.INTER.EQ.4)) CALL FLTABL

C...Get integrated cross-section.
      PARL(23)=0.
      IF(LST(10).GT.0) CALL LXSECT
      IF(LQCD.EQ.2.OR.LTM.EQ.2) THEN
        WRITE(6,1300)
        IF(LQCD.EQ.2) WRITE(6,1310)
        IF(LTM .EQ.2) WRITE(6,1320)
        WRITE(6,1330)
      ENDIF

      IF(LST(2).EQ.1) THEN
C...Find max value of differential cross section for rejection.
        UKIN(1)=(XMAX+XMIN)/2.
        WKIN(1)=0.8*(XMAX-XMIN)/2.
        AIN(1)=XMIN
        BIN(1)=XMAX
        IF(LST(31).EQ.1) THEN
          UKIN(2)=(Q2MAX+Q2MIN)/2.
          WKIN(2)=0.8*(Q2MAX-Q2MIN)/2.
          AIN(2)=Q2MIN
          BIN(2)=Q2MAX
          NAMKIN(2)='      Q**2'
        ELSEIF(LST(31).EQ.2) THEN
          UKIN(2)=(YMAX+YMIN)/2.
          WKIN(2)=0.8*(YMAX-YMIN)/2.
          AIN(2)=YMIN
          BIN(2)=YMAX
          NAMKIN(2)='         y'
        ELSEIF(LST(31).EQ.3) THEN
          UKIN(2)=(W2MAX+W2MIN)/2.
          WKIN(2)=0.8*(W2MAX-W2MIN)/2.
          AIN(2)=W2MIN
          BIN(2)=W2MAX
          NAMKIN(2)='      W**2'
        ENDIF
C...Maximum obtained by minimizing -(diff. x-section).
        CALL LTIMEX(TI1)
        CALL LMINEW
        CALL LTIMEX(TI2)
        PARI(LST(23))=FCNMAX*1.1
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &  WRITE(6,1400) PARI(LST(23)),TI2-TI1
      ENDIF

      IF(LFILE.GT.0.AND.LST(19).GE.0) THEN
C...Read QCD weights from file.
        READ(LFILE) LSTW,PARLW,NXX,NWW,NP,XX,WW
        IPMAX=2
        IF(LSTW(17).NE.0) IPMAX=3
        READ(LFILE) (((PQG(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,NP),
     &  (((PQQB(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,NP),
     &  (((QGMAX(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,IPMAX),
     &  (((QQBMAX(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,MIN(2,IPMAX)),
     &  YCUT
        IF(NP.NE.1) READ(LFILE) XTOT
        CLOSE(LFILE)
C...Reset parameters for matrix element integration.
        PARL(8)=PARLW(8)
        PARL(9)=PARLW(9)
        PARL(11)=PARLW(11)
        PARL(12)=PARLW(12)
        PARL(13)=PARLW(13)
C...Check current parameter values against those used when
C...calculating weights.
        IF(LST(12).NE.LSTW(12).OR.LST(13).NE.LSTW(13)
     &  .OR.LST(15).NE.LSTW(15).OR.LST(16).NE.LSTW(16)
     &  .OR.LST(17).NE.LSTW(17).OR.LST(23).NE.LSTW(23)
     &  .OR.ABS(PARL(1)-PARLW(1)).GT.0.1.OR.ABS(PARL(2)-PARLW(2)).GT.0.1
     &  .OR.ABS(PARL(5)-PARLW(5)).GT.0.01
     &  .OR.ABS(PARL(6)-PARLW(6)).GT.0.1) THEN
         IF(LST(3).GE.1)
     &    WRITE(6,1500) LST(12),LSTW(12),LST(13),LSTW(13),LST(15),
     &    LSTW(15),LST(16),LSTW(16),LST(17),LSTW(17),LST(23),LSTW(23),
     &    PARL(1),PARLW(1),PARL(2),PARLW(2),PARL(5),PARLW(5),PARL(6),
     &    PARLW(6)
          IF(LST(3).GE.2) THEN
            WRITE(6,1900)
            STOP
          ENDIF
        ENDIF
      ELSEIF((LST(19).GE.0.OR.LST(19).EQ.-10).AND.
     &(LST(8).EQ.1.OR.LST(8)/10.EQ.1.OR.MOD(LST(8),10).EQ.9)) THEN
C...Calculate weights if 1st order QCD from grid is requested.
        CALL LTIMEX(TI1)
        CALL LWEITS(LFILE)
        CALL LTIMEX(TI2)
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &  WRITE(6,1510) TI2-TI1
      ENDIF

C...Reset counters to zero for Monte Carlo estimate of cross section.
      PARI(27)=0.
      PARI(28)=0.
      PARI(29)=0.
      PARI(30)=0.
      LST(32)=0
      RETURN

 1000 FORMAT(' ',//,5X,
     &'A MONTE CARLO GENERATOR FOR DEEP INELASTIC LEPTON-'
     &,'NUCLEON SCATTERING',/,5X,68('='),//,
     &25X,'LEPTO version 6.5.1, October 31, 1996',//,
C    &25X,'PRELIMINARY VERSION, DO NOT CIRCULATE',//,
     &' Lepton: type =',I3,5X,'momentum (px,py,pz) =',3F8.4,
     &' GeV',//,' Target: A, Z =',2F3.0,2X,
     &'momentum (px,py,pz) =',3F8.4,' GeV',//,
     &' Interaction :',I3,14X,' CMS energy =',1PG12.4,' GeV',/)
 1010 FORMAT(' Warning: lepton and nucleon momenta in same direction',
     &' not allowed.',/,10X,'Execution stopped.')
 1020 FORMAT(/,' JETSET version ',I3,'.',I3,' is used.',/,
     &' Parton densities in PYTHIA version ',I3,'.',I3,' are used.',/)
 1030   FORMAT(' Warning (LINIT): JETSET version before 7.402, MSTJ(46)'
     &  ,' set to',I4,/,18X,'to avoid mismatch LEPTO<-->LUSHOW.',/)
 1050 FORMAT(/,' User applied cuts (+ phase space) : ',1P,
     &      G12.4,' <   x   < ',G12.4,
     &/,37X,G12.4,' <   y   < ',G12.4,
     &/,37X,G12.4,' < Q**2  < ',G12.4,
     &/,37X,G12.4,' < W**2  < ',G12.4,
     &/,37X,G12.4,' <  nu   < ',G12.4,
     &/,37X,G12.4,' <  E''   < ',G12.4,
     &/,37X,G12.4,' < theta < ',G12.4,/,
     &/,       ' Effective ranges (from above cuts): ',
     &      G12.4,' <   x   < ',G12.4,
     &/,37X,G12.4,' <   y   < ',G12.4,
     &/,37X,G12.4,' < Q**2  < ',G12.4,
     &/,37X,G12.4,' < W**2  < ',G12.4,
     &/,37X,G12.4,' <  nu   < ',G12.4)
 1100 FORMAT(' Warning: effective upper limit of kinematical ',
     &'variable(s) smaller than corresponding lower limit.')
 1110 FORMAT(' Warning: lower limit in x and/or Q2 too small for ',
     &'DIS formalism.')
 1150 FORMAT(' Warning: weak charged current cross section zero for ',
     &'specified lepton helicity; LEPIN, PARL(6) =',I3,F5.2)
 1200 FORMAT(' Warning: unrecognized interaction in LINIT call: ',
     &'INTER = ',I5,'  for lepton LEPIN =',I5)
 1210 FORMAT(' Warning: unallowed value of LST(1) =',I3,
     &' and/or LST(31) =',I3)
 1220 FORMAT(/,' User-defined optimization parameters:',
     &/,5X,'OPTX(1...4)  =',4G11.3,/,5X,'OPTY(1...4)  =',4G11.3,
     &/,5X,'OPYQ2(1...4) =',4G11.3,/,5X,'OPTW2(1...4) =',4G11.3,/)
 1250 FORMAT(/,' Parameter values:',//,9X,'I',4X,'LST(I)',1X,
     &'LST(I+10)',8X,'PARL(I)',5X,'PARL(I+10)',1P,
     &/,5X,55('-'),10(/,3I10,2G15.4),/)
 1300 FORMAT(' Warning: cross section, PARL(23), excludes FL (see ',
     &'LST(11)) from:')
 1310 FORMAT(10X,'QCD, since evaluated event by event for LQCD=2')
 1320 FORMAT(10X,'TM , since evaluated event by event for LTM =2')
 1330 FORMAT(' Cross section in PARL(24) includes these contributions.')
 1400 FORMAT(' Max of differential cross section (for weighting) =',
     &E12.4,/,' obtained in ',F7.2,' seconds.',/)
 1500 FORMAT(//,' Warning: current parameter values do not match ',
     &'with those used when calculating QCD weights.',//,15X,
     &'current value     value for weights',/,
     &/,'     LST(12)   ',I12,10X,I12,
     &/,'     LST(13)   ',I12,10X,I12,
     &/,'     LST(15)   ',I12,10X,I12,
     &/,'     LST(16)   ',I12,10X,I12,
     &/,'     LST(17)   ',I12,10X,I12,
     &/,'     LST(23)   ',I12,10X,I12,
     &/,'     PARL(1)   ',E12.4,10X,E12.4,
     &/,'     PARL(2)   ',E12.4,10X,E12.4,
     &/,'     PARL(5)   ',E12.4,10X,E12.4,
     &/,'     PARL(6)   ',E12.4,10X,E12.4)
 1510 FORMAT(/,' Time for calculating QCD weights =',F5.1,' seconds',/)
 1900 FORMAT(' Execution stopped ',/)
      END

C **********************************************************************

      SUBROUTINE LEPTO

C...Administer the generation of an event.
C...Note: if error flag LST(21) is non-zero, no proper event generated.

      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
      COMMON /ARDAT1/ PARA(40),MSTA(40)
      DOUBLE PRECISION DTHETA,DPHI,DBETA,DETOT,DARI29,DARI30
      DIMENSION SPQ(17)
      DATA NUMMIS,NWARN/0,10/,DARI29,DARI30/2*0.D0/

      DOUBLE PRECISION Target
      DOUBLE PRECISION KFermi
      DOUBLE PRECISION CosTheta
      DOUBLE PRECISION Phi
      DOUBLE PRECISION KFx
      DOUBLE PRECISION KFy
      DOUBLE PRECISION KFz
      
      L17=0
    1 LST(21)=0
      DO 10 I=1,10
      DO 10 J=1,5
      K(I,J)=0
   10 V(I,J)=0.
      DO 15 I=1,4
      K(I,1)=21
   15 K(I,2)=KSAVE(I)
      K(4,1)=1
      N=2
      IF(LST(17).NE.0.AND.LST(2).GT.0) THEN
C...Lepton and/or nucleon energy may vary from event to event,
        IF(L17.EQ.0) THEN

C...Calling for csparam to be taken as vaiables nucleon init mom
      Target = PARL(1) 
      CAll csparam(Target, KFermi, CosTheta, Phi, KFx, KFy, KFz);
c.....Target = PARL(1) 
      P(2,1)  = KFx
      P(2,2)  = KFy
      P(2,3)  = KFz
      P(2,4)=SQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)
c.....P(2,4)=SQRT(P(2,3)**2+P(2,5)**2) for lst(17) = 0
c.....P(2,3)  = KFermi
            
C...Momentum vectors from P(i,j) i=1,2 j=1,2,3 on entry in LEPTO
          DO 20 I=1,2
          P(I,5)=ULMASS(K(I,2))
          P(I,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2)
          DO 20 J=1,5
   20     PSAVE(3,I,J)=P(I,J)
        ELSE
C...Momentum vectors from PSAVE if new try, i.e. jump back to 1
          DO 25 I=1,2
          DO 25 J=1,5
   25     P(I,J)=PSAVE(3,I,J)
        ENDIF
        L17=1
C...Transform to cms of incoming particles, lepton along +z axis.
        DO 30 J=1,3
   30   DBETA(1,J)=(DBLE(P(1,J))+DBLE(P(2,J)))/
     &             (DBLE(P(1,4))+DBLE(P(2,4)))
        CALL LUDBRB(0,0,0.,0.,-DBETA(1,1),-DBETA(1,2),-DBETA(1,3))
        SPHI(1)=ULANGL(P(1,1),P(1,2))
        CALL LUDBRB(0,0,0.,-SPHI(1),0.D0,0.D0,0.D0)
        STHETA(1)=ULANGL(P(1,3),P(1,1))
        CALL LUDBRB(0,0,-STHETA(1),0.,0.D0,0.D0,0.D0)
        LST(28)=2
        PARL(21)=2.*(P(1,4)*P(2,4)-P(1,3)*P(2,3))
      ELSE
C...Initial state momenta fixed from LINIT call.
        DO 42 I=1,2
        DO 40 J=1,5
   40   P(I,J)=PSAVE(3,I,J)
   42   IF(PSAVE(3,1,3).LT.0.) P(I,3)=-PSAVE(3,I,3)
        LST(28)=3
      ENDIF

      CALL LEPTOX
C...Return if error or if no event to be generated.
      IF(LST(21).NE.0.OR.LST(2).LE.0.OR.LST(7).EQ.-1) RETURN

      IF(PARI(29).LT.0.5) THEN
C...For first call, reset double precision counters.
        DARI29=0.D0
        DARI30=0.D0
      ENDIF
      DARI29=DARI29+1.D0
      PARI(29)=DARI29

C     CALL GULIST(-3,2)
C...Scattered lepton and exchanged boson added to event record in LKINEM
C...Transform to lepton-nucleon cms if not made earlier
      IF(LST(17).EQ.0) THEN
        DO 46 I=3,4
        DO 45 J=1,5
   45   PSAVE(3,I,J)=P(I,J)
   46   IF(PSAVE(3,1,3).LT.0.) PSAVE(3,I,3)=-P(I,3)
        CALL LUDBRB(0,0,0.,0.,0.D0,0.D0,-DBETA(1,3))
        LST(28)=2
      ENDIF
      DO 50 I=1,4
      DO 50 J=1,5
   50 PSAVE(2,I,J)=P(I,J)
C     CALL GULIST(-2,2)

C...Prepare for parton cascade.
      IF(LST(8).GE.2.AND.MOD(LST(8),10).NE.9) CALL LSHOWR(0)

C...Transform to hadronic cms, boost parameters in double precision.
      DETOT=DBLE(P(1,4))-DBLE(P(4,4))+DBLE(P(2,4))
      DBETA(2,1)=-DBLE(P(4,1))/DETOT
      DBETA(2,2)=-DBLE(P(4,2))/DETOT
      DBETA(2,3)=(DBLE(P(1,3))-DBLE(P(4,3))+DBLE(P(2,3)))/DETOT
      CALL LUDBRB(0,0,0.,0.,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
      SPHI(2)=0.
      STHETA(2)=ULANGL(P(3,3),P(3,1))
      CALL LUDBRB(0,0,-STHETA(2),0.,0.D0,0.D0,0.D0)
      LST(28)=1
      DO 60 I=1,4
      DO 60 J=1,5
   60 PSAVE(1,I,J)=P(I,J)
C...Save momentum of exchanged boson (used in subroutine LFRAME).
      DO 70 J=1,5
   70 PB(J)=P(3,J)
C     CALL GULIST(-1,2)

   90 N=4
      MSTU(1)=N+1
      LST(26)=N+1
      LST(27)=0
      PARL(25)=ULALPS(Q2)
      IF(LST(8).EQ.1.OR.LST(8)/10.EQ.1.OR.MOD(LST(8),10).EQ.9) THEN
C...Probabilities for hard, first order QCD events.
CAE...Corrected what to do when LQGEV or LQQBEV fail. Now make LQEV.
        CALL LQCDPR(QG,QQB)
        
CAE        WRITE(6,*) 'Lepto:',X,Q2,PARL(27)
CAE        WRITE(6,*) 'weights:',1-QG-QQB,QG,QQB
        
        DO 100 I=1,17
  100   SPQ(I)=PQ(I)
  200   SRLU=RLU(0)
        IF(SRLU.GT.QQB+QG) THEN
          CALL LQEV
        ELSEIF(SRLU.GT.QQB) THEN
          IF(LST(8).EQ.9) THEN
            DO 211 I=1,17
  211         PQ(I)=SPQ(I)
            CALL LQEV
          ELSE
            CALL LQGEV
            DO 212 I=1,17
              PQ(I)=SPQ(I)
  212       CONTINUE
          ENDIF
        ELSE
          CALL LQQBEV
          DO 213 I=1,17
            PQ(I)=SPQ(I)
  213     CONTINUE
          IF(LST(8).EQ.9.AND.LST(21).EQ.0) THEN
            IF(PLU(5,11).LT.Q2*PARA(20)) THEN
              DO 220 I=1,17
  220         PQ(I)=SPQ(I)
              CALL LQEVAR(K(5,2),K(7,2))
            ENDIF
          ENDIF
        ENDIF
        IF(LST(21).NE.0) THEN
CAE          WRITE(6,*) 'error:',LST(21),LST(24)
  230     CALL LQEV
          IF(LST(21).NE.0) GOTO 230
        ENDIF   
      ELSE
C...QPM model without QCD corrections (cascade applied later).
  300   CALL LQEV
        IF(LST(21).NE.0) GOTO 300
      ENDIF

      NS=MSTU(1)
      MSTU(1)=0
C     CALL GULIST(-3,2)
C     WRITE(6,*) ' LST(24)=',LST(24)
CJR--       no preclustering of small systems
          MSTJ(14)=-1
CJR--            
      IF(LST(8).LE.1.OR.MOD(LST(8),10).EQ.9) THEN
C...No parton cascade, introduce primordial kt.
        IF(PARL(3).GT.1.E-03) THEN
          CALL LPRIKT(PARL(3),PT,PHI)
          CALL LUDBRB(NS,N,0.,-PHI,0.D0,0.D0,0.D0)
          CALL LUDBRB(NS,N,ATAN(2.*PT/SQRT(W2)),PHI,0.D0,0.D0,0.D0)
        ENDIF
        IF(MOD(LST(8),10).NE.9) THEN
C...Check system against fragmentation cuts.
          MSTU(24)=0
          CALL LUPREP(0)
          IF(MSTU(24).NE.0) THEN
            IF(LST(3).GE.1) WRITE(6,*)'LUPREP error MSTU(24)=',MSTU(24),
     &                                ', New event generated'
            LST(21)=1
            GOTO 1
          ENDIF
        ENDIF
      ELSEIF(LST(24).EQ.1) THEN
C...Include parton cascades (+ remnant & kt) on q-event
        CALL LSHOWR(1)
      ELSE
C...Include parton cascades (+ remnant & kt) on qg- or qqbar-event
        CALL LMEPS
      ENDIF
      IF(LST(21).NE.0) THEN
C        IF(LST(3).GE.1)
C     &     WRITE(6,*)'Cascade error LST(21)= ',LST(21),
C     &               ', New event generated'
        GOTO 1
      ENDIF

CJR--       Soft colour interactions
      IF(LST(34).EQ.1 .OR. LST(34).EQ.2) CALL LSCI(PARL(7))
      IF(LST(21).NE.0) GOTO 1
CJR--       take care of small systems
      CALL LSMALL
      IF(LST(21).NE.0) THEN
        IF(LST(3).GE.1) WRITE(6,*)' LSMALL error LST(21)= ',LST(21),
     &                            ', New event generated'
        GOTO 1
      ENDIF
      MSTJ(14)=1
      CALL LUPREP(0)
      IF(MSTU(24).NE.0) THEN
         IF(LST(3).GE.1) WRITE(6,*)' LUPREP error MSTU(24)= ',MSTU(24),
     &                             ', New event generated'
         LST(21)=1
      ENDIF
CJR--            
      IF(LST(21).NE.0) GOTO 1

      DO 400 I=1,N
C...Correct energy-momentum-mass mismatch for real particle
      IF(P(I,5).LT.0.) GOTO 400
      ENERGY=SQRT(DBLE(P(I,5))**2+DBLE(P(I,1))**2+DBLE(P(I,2))**2+
     &DBLE(P(I,3))**2)
      P2=DBLE(P(I,4))**2-DBLE(P(I,1))**2-DBLE(P(I,2))**2-DBLE(P(I,3))**2
      IF(ABS(ENERGY-P(I,4))/(PSAVE(3,1,4)+PSAVE(3,2,4)).GT.PARU(11))THEN
        NUMMIS=NUMMIS+1
C...For testing purposes
C       IF(LST(3).GE.1.AND.NUMMIS.LE.NWARN) THEN
C         WRITE(6,1000) I,(K(I,J),J=1,2),(P(I,J),J=1,5),
C    &    SIGN(SQRT(ABS(P2)),P2),ENERGY,INT(DARI29),NWARN
C         IF(ABS(P2-P(I,5)**2).GT.400.) CALL LULIST(2)
C       ENDIF
CAE        WRITE(6,*) 'Energy mismatch',LST(24),PARL(28),PARL(29),NUMMIS
        GOTO 90
      ENDIF
      P(I,4)=ENERGY
  400 CONTINUE

      DARI30=DARI30+1.D0
      PARI(30)=DARI30
Ctest IF(LST(23).EQ.2) PARL(24)=PARL(24)*DARI30/DARI29

      DO 500 I=1,N
      DO 500 J=1,5
  500 V(I,J)=0.
      IF(LST(7).EQ.1) THEN
        CALL LUEXEC
        IF(MSTU(24).NE.0) THEN
          WRITE(6,*) ' Error from JETSET, new event made'
          GOTO 90
        ENDIF
      ENDIF

C     CALL GULIST(-1,2)
C...Transform to desired frame
C     LST(28)=1
      LST(29)=0
      PHIR=6.2832*RLU(0)
      IF(LST(17).EQ.0) THEN
        IF(LST(5).GE.2) CALL LFRAME(LST(5),0)
C...Restore momenta (e,p,boson,L) due to numerical errors from boosts
        DO 600 I=1,4
        DO 600 J=1,5
  600   P(I,J)=PSAVE(LST(28),I,J)
        IF(LST(6).EQ.1.AND.LST(28).GE.2) THEN
C...Random rotation in azimuthal angle
          CALL LUDBRB(0,0,0.,PHIR,0.D0,0.D0,0.D0)
          LST(29)=1
        ENDIF
      ELSE
        IF(LST(5).GE.2) CALL LFRAME(LST(5),LST(6))
      ENDIF
C...Deactivate scattered lepton
      IF(MOD(LST(4),10).EQ.0) K(4,1)=21
C     CALL GULIST(0,2)

      RETURN
 1000 FORMAT(' Warning: too large numerical mismatch in ',
     &'particle energy-momentum-mass',
     &/,3X,'I K(I,1) ..2)  P(I,1)  P(I,2)  P(I,3)',
     &'  P(I,4)  P(I,5)    mass  energy',/,I4,2I6,7F8.3,/,
     &' Event no.',I8,' regenerated. Only first',I5,' warnings printed')
      END

C **********************************************************************

      SUBROUTINE LEPTOX

C...Select process and choose kinematical variables (x,y; x,Q2; x,W2)
C...according to the differential cross section.

      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON /FLINFO/ RFLQ,RFLG,RFLM,RFLT
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      DIMENSION PQH(17,2),PNT(2,2),XPQ(-6:6)
      DOUBLE PRECISION DARI27,DARI28
      DATA DARI27,DARI28/2*0.D0/
      DATA W2LOW,W2UPP,YLOW,YUPP,Q2LOW,Q2UPP/6*0./

      DO 10 IH=1,2
      DO 5 I=1,2
    5 PNT(I,IH)=0.
      DO 6 I=1,8
      EWQC(1,IH,I)=0.
    6 EWQC(2,IH,I)=0.
      DO 10 I=1,17
   10 PQH(I,IH)=0.
      DO 20 I=1,17
   20 PQ(I)=0.

      LST(21)=0
      NCUT=0
      S=PARL(21)
      PM2=PSAVE(3,2,5)**2
      IF(LST(2).NE.1) THEN
       Q2LOW=MAX(Q2MIN,X*YMIN*S,(W2MIN-PM2)*X/MAX(1.-X,1.E-22))
       Q2UPP=MIN(Q2MAX,X*YMAX*S,(W2MAX-PM2)*X/MAX(1.-X,1.E-22))
       YLOW=MAX(YMIN,Q2MIN/MAX(S*X,1.E-22),
     & (W2MIN-PM2)/MAX(S*(1.-X),1.E-22))
       YUPP=MIN(YMAX,Q2MAX/MAX(S*X,1.E-22),
     & (W2MAX-PM2)/MAX(S*(1.-X),1.E-22))
       W2LOW=MAX(W2MIN,(1.-X)*YMIN*S+PM2,Q2MIN*(1.-X)/MAX(X,1.E-22)+PM2)
       W2UPP=MIN(W2MAX,(1.-X)*YMAX*S+PM2,Q2MAX*(1.-X)/MAX(X,1.E-22)+PM2)
       GOTO 110
      ENDIF

      IF(PARI(28).LT.0.5) THEN
C...For first call, reset double precision counters.
        DARI27=0.D0
        DARI28=0.D0
      ENDIF
  100 DARI28=DARI28+1.D0
      PARI(28)=DARI28
  101 CONTINUE
C...Choose x according to the distribution
C...hx(x) =  a + b/x + c/x**2 + d/x**3. In detail
C...hq=OPTX(1)/(XMAX-XMIN) + 1/ln(XMAX/XMIN)*OPTX(2)/X
C...   +XMIN*XMAX/(XMAX-XMIN)*OPTX(3)/X**2
C...   +2*(XMIN*XMAX)**2/(XMAX**2-XMIN**2)*OPTX(4)/X**3
      WHICH=(OPTX(1)+OPTX(2)+OPTX(3)+OPTX(4))*RLU(0)
      IF(WHICH.LE.OPTX(1)) THEN
        X=XMIN+RLU(0)*(XMAX-XMIN)
      ELSEIF(WHICH.LE.(OPTX(1)+OPTX(2))) THEN
        X=XMIN*(XMAX/XMIN)**RLU(0)
      ELSEIF(WHICH.LE.(OPTX(1)+OPTX(2)+OPTX(3))) THEN
        X=XMIN*XMAX/(XMAX+RLU(0)*(XMIN-XMAX))
      ELSE
        X=SQRT((XMIN*XMAX)**2/(XMAX**2+RLU(0)*(XMIN**2-XMAX**2)))
      ENDIF
      IF(LST(31).EQ.1) THEN
C...Choose Q**2 according to the distribution
C...hq(Q2) =  a + b/(Q2) + c/(Q2)**2 + d/(Q2)**3. In detail
C...hq=OPTQ2(1)/(Q2MAX-Q2MIN) + 1/ln(Q2MAX/Q2MIN)*OPTQ2(2)/Q2
C...   +Q2MIN*Q2MAX/(Q2MAX-Q2MIN)*OPTQ2(3)/Q2**2
C...   +2*(Q2MIN*Q2MAX)**2/(Q2MAX**2-Q2MIN**2)*OPTQ2(4)/Q2**3
      

C.......TEST 
c        Q2LOW=MAX(0.5*Q2MIN,0.5*X*YMIN*S,0.5*(W2MIN-PM2)*X/(1.-X))
c        Q2UPP=MIN(2*Q2MAX,2*X*YMAX*S,2*(W2MAX-PM2)*X/(1.-X))
        Q2LOW=MAX(Q2MIN,X*YMIN*S,(W2MIN-PM2)*X/(1.-X))
        Q2UPP=MIN(Q2MAX,X*YMAX*S,(W2MAX-PM2)*X/(1.-X))
       
        IF(Q2UPP.LT.Q2LOW) GOTO 101
        WHICH=(OPTQ2(1)+OPTQ2(2)+OPTQ2(3)+OPTQ2(4))*RLU(0)
        IF(WHICH.LE.OPTQ2(1)) THEN
          Q2=Q2LOW+RLU(0)*(Q2UPP-Q2LOW)
        ELSEIF(WHICH.LE.(OPTQ2(1)+OPTQ2(2))) THEN
          Q2=Q2LOW*(Q2UPP/Q2LOW)**RLU(0)
        ELSEIF(WHICH.LE.(OPTQ2(1)+OPTQ2(2)+OPTQ2(3))) THEN
          Q2=Q2LOW*Q2UPP/(Q2UPP+RLU(0)*(Q2LOW-Q2UPP))
        ELSE
         Q2=SQRT((Q2LOW*Q2UPP)**2/(Q2UPP**2+RLU(0)*(Q2LOW**2-Q2UPP**2)))
        ENDIF
        Y=Q2/(PARL(21)*X)
        IF(Y.LT.YMIN.OR.Y.GT.YMAX) GOTO 100
      ELSEIF(LST(31).EQ.2) THEN
C...Choose y according to the distribution
C...hy(y) =  a + b/y + c/y**2 + d/y**3. In detail
C...hy=OPTY(1)/(YMAX-YMIN) + 1/ln(YMAX/YMIN)*OPTY(2)/Y
C...   +YMIN*YMAX/(YMAX-YMIN)*OPTY(3)/Y**2
C...   +2*(YMIN*YMAX)**2/(YMAX**2-YMIN**2)*OPTY(4)/Y**3
        YLOW=MAX(YMIN,Q2MIN/(S*X),(W2MIN-PM2)/(S*(1.-X)))
        YUPP=MIN(YMAX,Q2MAX/(S*X),(W2MAX-PM2)/(S*(1.-X)))
        IF(YUPP.LT.YLOW) GOTO 101
        WHICH=(OPTY(1)+OPTY(2)+OPTY(3)+OPTY(4))*RLU(0)
        IF(WHICH.LE.OPTY(1)) THEN
          Y=YLOW+RLU(0)*(YUPP-YLOW)
        ELSEIF(WHICH.LE.(OPTY(1)+OPTY(2))) THEN
          Y=YLOW*(YUPP/YLOW)**RLU(0)
        ELSEIF(WHICH.LE.(OPTY(1)+OPTY(2)+OPTY(3))) THEN
          Y=YLOW*YUPP/(YUPP+RLU(0)*(YUPP-YLOW))
        ELSE
          Y=SQRT((YLOW*YUPP)**2/(YUPP**2+RLU(0)*(YLOW**2-YUPP**2)))
        ENDIF
        Q2=X*Y*PARL(21)
        IF(Q2.LT.Q2MIN.OR.Q2.GT.Q2MAX) GOTO 100
      ELSEIF(LST(31).EQ.3) THEN
C...Choose W**2 according to the distribution
C...hw(W2) =  a + b/(W2) + c/(W2)**2 + d/(W2)**3. In detail
C...hw=OPTW2(1)/(W2MAX-W2MIN) + 1/ln(W2MAX/W2MIN)*OPTW2(2)/W2
C...   +W2MIN*W2MAX/(W2MAX-W2MIN)*OPTW2(3)/W2**2
C...   +2*(W2MIN*W2MAX)**2/(W2MAX**2-W2MIN**2)*OPTW2(4)/W2**3
        W2LOW=MAX(W2MIN,(1.-X)*YMIN*S+PM2,Q2MIN*(1.-X)/X+PM2)
        W2UPP=MIN(W2MAX,(1.-X)*YMAX*S+PM2,Q2MAX*(1.-X)/X+PM2)
        IF(W2UPP.LT.W2LOW) GOTO 101
        WHICH=(OPTW2(1)+OPTW2(2)+OPTW2(3)+OPTW2(4))*RLU(0)
        IF(WHICH.LE.OPTW2(1)) THEN
          W2=W2LOW+RLU(0)*(W2UPP-W2LOW)
        ELSEIF(WHICH.LE.(OPTW2(1)+OPTW2(2))) THEN
          W2=W2LOW*(W2UPP/W2LOW)**RLU(0)
        ELSEIF(WHICH.LE.(OPTW2(1)+OPTW2(2)+OPTW2(3))) THEN
          W2=W2LOW*W2UPP/(W2UPP+RLU(0)*(W2LOW-W2UPP))
        ELSE
         W2=SQRT((W2LOW*W2UPP)**2/(W2UPP**2+RLU(0)*(W2LOW**2-W2UPP**2)))
        ENDIF
        Y=(W2-PM2)/((1.-X)*PARL(21))
        Q2=X*Y*PARL(21)
        IF(Y.LT.YMIN.OR.Y.GT.YMAX) GOTO 100
        IF(Q2.LT.Q2MIN.OR.Q2.GT.Q2MAX) GOTO 100
      ENDIF
  110 IF(LKINEM(LST(2)).NE.0) THEN
        NCUT=NCUT+1
        IF(LST(2).EQ.1) THEN
          IF(NCUT.LE.9999) GOTO 100
          IF(LST(3).GE.1) WRITE(6,1200)
        ENDIF
        LST(21)=2
        RETURN
      ENDIF

      PARI(24)=(1.+(1.-Y)**2)/2.
      PARI(25)=1.-Y
      PARI(26)=(1.-(1.-Y)**2)/2.
      CALL LNSTRF(X,Q2,XPQ)
C...Lepton helicity state, only one contributes in some cases.
      IH=1
      IF(PARL(6).GT.+0.99) IH=2
  200 LST(30)=SIGN(1.,IH-1.5)
      PQH(17,IH)=0.
      PNT(1,IH)=0.
      PNT(2,IH)=0.
      IF(LST(23).EQ.2) THEN
C...Charged current: zero cross-section for one helicity state.
        IF(KSAVE(1).LT.0.AND.IH.EQ.1
     &  .OR.KSAVE(1).GT.0.AND.IH.EQ.2) GOTO 240
        YQ=PARI(24)-LST(30)*PARI(26)
        YQB=PARI(24)+LST(30)*PARI(26)
        IF(PARI(11).GT.1.E-06) THEN
          IF(K(3,2).LT.0) THEN
            PNT(1,IH)=(1.-PARI(11))*PARI(13)*YQ
            PNT(2,IH)=PARI(11)*PARI(12)*YQ
          ELSE
            PNT(1,IH)=(1.-PARI(11))*PARI(12)*YQ
            PNT(2,IH)=PARI(11)*PARI(13)*YQ
          ENDIF
        ENDIF
        DO 220 I=1,LST(12)
        IF(K(3,2)*QC(I).LT.0) THEN
          PQH(I,IH)=XPQ(I)*YQ
        ELSE
          PQH(I+LST(12),IH)=XPQ(-I)*YQB
        ENDIF
  220   CONTINUE
      ELSE
C...Neutral current: electromagnetic or weak or both with interference.
        GFQ2=Q2/(PMAS(23,1)**2+Q2)*SQRT(2.)*PARL(17)*PMAS(23,1)**2/
     &  (3.1415927*PARL(16))
C...Correction to obtain Q**2 dependent alpha-em, if desired.
        AEMCOR=1.
        IF(LST(18).GE.2) AEMCOR=ULALEM(Q2)/PARL(16)
        II=3-IH
        ZLEP=ZL(IH,ILEP+2*INU)
        DO 230 I=1,MAX(LST(12),LST(13))
        A=(-IG*QC(I)*AEMCOR+IZ*GFQ2*ZLEP*ZQ(IH,I))**2
        B=(-IG*QC(I)*AEMCOR+IZ*GFQ2*ZLEP*ZQ(II,I))**2
C...Save helicity-dependent electroweak quark couplings for later use.
        EWQC(1,IH,I)=A
        EWQC(2,IH,I)=B
        IF(I.GT.LST(12)) GOTO 230
        FYQ=(A+B)*PARI(24)+(A-B)*PARI(26)
        PQH(I,IH)=XPQ(I)*FYQ
        IF(I.LE.2.AND.PARI(11).GT.1.E-06) THEN
          PNT(1,IH)=PNT(1,IH)+(1.-PARI(11))*PARI(11+I)*FYQ
          PNT(2,IH)=PNT(2,IH)+PARI(11)*PARI(14-I)*FYQ
        ENDIF
        PQH(I+LST(12),IH)=XPQ(-I)*((A+B)*PARI(24)-(A-B)*PARI(26))
  230   CONTINUE
      ENDIF
  240 CONTINUE
      DO 300 I=1,LST(12)
  300 PQH(17,IH)=PQH(17,IH)+PQH(I,IH)+PQH(I+LST(12),IH)

      IF(ABS(PARL(6)).LT.0.99.AND.IH.EQ.1) THEN
        IH=2
        GOTO 200
      ENDIF

      FLQ=0.
      FLG=0.
      FLM=0.
      FLT=0.
      IF(LST(11).NE.0.AND.(LST(23).EQ.1.OR.LST(23).EQ.4)
     &.AND.LST(2).NE.-3) THEN
C...Include F_L for photon exchange (unless QCD grid being set up)
        LQCD=MOD(LST(11),10)
        LTM=MOD(LST(11)/10,10)
        LHT=LST(11)/100
C...Include QCD, target mass and/or higher twist contr. to long. str fcn
C...FL from interpolation.
        IF(LQCD.EQ.1.OR.LTM.EQ.1) CALL FLIPOL(FLQ,FLG,FLM)
C...Event simulation: if requested, get FL by event-by-event integration
        IF(LST(2).GT.0.AND.
     &  (LQCD.EQ.2.OR.LTM.EQ.2)) CALL FLINTG(FLQ,FLG,FLM)
        IF(LTM.GE.1.OR.LHT.GE.1) THEN
          F2EM=0.
          DO 301 I=1,LST(12)
  301     F2EM=F2EM+QC(I)**2*(XPQ(I)+XPQ(-I))
          IF(LTM.GE.1) FLM=FLM-2.*X**2*PSAVE(3,2,5)**2/Q2*F2EM
          IF(LHT.GE.1) FLT=8.*PARL(19)/Q2*F2EM
        ENDIF
        DO 305 IH=1,2
        PQH17=PQH(17,IH)
C...Note factor 2 at the end, since PQH(IH,17) contains overall factor 2
        PQH(17,IH)=PQH(17,IH)-Y**2*(FLQ+FLG+FLM+FLT)
        DO 305 I=1,16
  305   PQH(I,IH)=PQH(I,IH)*PQH(17,IH)/PQH17
      ENDIF

      DO 310 I=1,17
  310 PQ(I)=(1.-PARL(6))/2.*PQH(I,1)+(1.+PARL(6))/2.*PQH(I,2)

C...Relative contribution from longitudinal str. fcn. and higher twist.
       RFLQ=-Y**2*FLQ/MAX(PQ(17),1.E-33)
       RFLG=-Y**2*FLG/MAX(PQ(17),1.E-33)
       RFLM=-Y**2*FLM/MAX(PQ(17),1.E-33)
       RFLT=-Y**2*FLT/MAX(PQ(17),1.E-33)

C...Common factor for matrix elements.
      IF(LST(31).EQ.1) THEN
        IF(LST(23).EQ.2) THEN
          COMFAC=1./X/(1.+Q2/PMAS(24,1)**2)**2
        ELSE
          COMFAC=1./X/Q2**2
        ENDIF
      ELSEIF(LST(31).EQ.2) THEN
        IF(LST(23).EQ.2) THEN
          COMFAC=1./(1.+Q2/PMAS(24,1)**2)**2*PARL(21)
        ELSE
          COMFAC=1./Q2**2*PARL(21)
        ENDIF
      ELSEIF(LST(31).EQ.3) THEN
        IF(LST(23).EQ.2) THEN
          COMFAC=1./X/(1.+Q2/PMAS(24,1)**2)**2  * X/(1.-X)
        ELSE
          COMFAC=1./X/Q2**2 * X/(1.-X)
        ENDIF
      ENDIF
C-check: Move change of COMFAC to below??
C...Prepare for Q2 weighting.


C TEST WEIGHT=1/Q2**2
      WEIGHT=1.D0
      COMFAC=COMFAC/WEIGHT
      IF(LST(2).LE.-2) RETURN
      HX=OPTX(1)/(XMAX-XMIN) + 1./ALOG(XMAX/XMIN)*OPTX(2)/X
     &+XMIN*XMAX/(XMAX-XMIN)*OPTX(3)/X**2
     &+2*(XMIN*XMAX)**2/(XMAX**2-XMIN**2)*OPTX(4)/X**3
      XFACT=OPTX(1)+OPTX(2)+OPTX(3)+OPTX(4)
      IF(LST(31).EQ.1) THEN
        HQ2=OPTQ2(1)/(Q2UPP-Q2LOW)
     &  +1./ALOG(Q2UPP/Q2LOW)*OPTQ2(2)/Q2
     &  +Q2LOW*Q2UPP/(Q2UPP-Q2LOW)*OPTQ2(3)/Q2**2
     &  +2*(Q2LOW*Q2UPP)**2/(Q2UPP**2-Q2LOW**2)*OPTQ2(4)/Q2**3
        Q2FACT=OPTQ2(1)+OPTQ2(2)+OPTQ2(3)+OPTQ2(4)
        COMFAC=COMFAC*XFACT*Q2FACT/HX/HQ2
      ELSEIF(LST(31).EQ.2) THEN
        HY=OPTY(1)/(YUPP-YLOW)+1./ALOG(YUPP/YLOW)*OPTY(2)/Y
     &  +YLOW*YUPP/(YUPP-YLOW)*OPTY(3)/Y**2
     &  +2*(YLOW*YUPP)**2/(YUPP**2-YLOW**2)*OPTY(4)/Y**3
        YFACT=OPTY(1)+OPTY(2)+OPTY(3)+OPTY(4)
        COMFAC=COMFAC*XFACT*YFACT/HX/HY
      ELSEIF(LST(31).EQ.3) THEN
        HW2=OPTW2(1)/(W2UPP-W2LOW)
     &  +1./ALOG(W2UPP/W2LOW)*OPTW2(2)/W2
     &  +W2LOW*W2UPP/(W2UPP-W2LOW)*OPTW2(3)/W2**2
     &  +2*(W2LOW*W2UPP)**2/(W2UPP**2-W2LOW**2)*OPTW2(4)/W2**3
        W2FACT=OPTW2(1)+OPTW2(2)+OPTW2(3)+OPTW2(4)
        COMFAC=COMFAC*XFACT*W2FACT/HX/HW2
      ENDIF
      IF(LST(2).LE.0) RETURN

C-check: Move change of COMFAC to here?
      SIGL=(1.-PARL(6))/2.*PQH(17,1)
      SIGR=(1.+PARL(6))/2.*PQH(17,2)
      SIGMA=SIGL+SIGR
      IF(LST(2).EQ.1) THEN
C...When chosing (x,y), reject according to maximum of "cross-section",
C...update cross-section estimate.
        DARI27=DARI27+DBLE(SIGMA)*DBLE(COMFAC)*WEIGHT
        PARI(27)=DARI27
        VIOL=SIGMA*COMFAC/PARI(LST(23))
        IF(VIOL.GT.PARI(32)) THEN
          PARI(32)=VIOL
          IF(PARI(32).GT.1.) THEN
            PARI(LST(23))=PARI(LST(23))*PARI(32)
            IF(LST(3).GE.1) WRITE(6,1300) PARI(32),INT(PARI(30)+1),
     &      PARI(LST(23)),X,Y,Q2,W2
            PARI(32)=1.
          ENDIF
        ENDIF
        IF(VIOL.LT.RLU(0)) GOTO 100
        PARL(24)=PARI(31)*DARI27/DARI28
      ENDIF

      IF(ABS(PARL(6)).LT.0.99) THEN
C...Choose helicity of incoming lepton.
        IH=1
        IF(RLU(0)*SIGMA.GT.SIGL) IH=2
      ENDIF
      LST(30)=SIGN(1.,IH-1.5)

C...Choose target nucleon, proton or neutron.
      LST(22)=1
      K(2,2)=2212
      IF(PARI(11).GT.1.E-06) THEN
        IF(RLU(0).LT.(PARI(11)*(PQH(17,IH)-PNT(1,IH)-PNT(2,IH))+
     &  PNT(2,IH))/PQH(17,IH)) THEN
          LST(22)=2
          K(2,2)=2112
        ENDIF
      ENDIF

      RETURN
 1200 FORMAT(' Warning: LEPTOX is looping, cannot find allowed ',
     &'phase space point due to cuts,',/,
     &10X,'check, in particular, CUT(11) to CUT(14)')
 1300 FORMAT(' Warning: maximum violated by a factor ',F7.3,
     &' in event ',I7,/,' maximum increased by this factor to ',E12.3,
     &/,' Point of violation: x, y, Q**2, W**2 = ',4G10.3)
      END

C **********************************************************************

      FUNCTION LKINEM(L)

C...Calculate kinematical variables and reject (optionally) if outside
C...required limits.

      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
      DOUBLE PRECISION DTHETA,DPHI,DBETA,DE,DPZ,DPT,DETOT

      LKINEM=1
      IF(L.EQ.-3) THEN
C...x,W known from LWEITS, no cuts applied.
        U=(W2-P(2,5)**2)/(2.*P(2,5)*(1.-X))
        Q2=2.*P(2,5)*U*X
        Y=Q2/(PARL(21)*X)
        GOTO 200
      ENDIF
C...x,y given.
      PARL(22)=Y*PARL(21)
      Q2=X*PARL(22)
      U=PARL(22)/(2.*P(2,5))
      W2=PARL(22)*(1.-X)+P(2,5)**2
      P(4,5)=ULMASS(K(4,2))
      IF(P(4,5)/SQRT(PARL(21)).LT.0.001) THEN
C...Simpler formulae for effectively massless scattered lepton.
        DE=DBLE(P(1,4))*(1.-DBLE(Y))+DBLE(X)*DBLE(Y)*DBLE(ABS(P(2,3)))
        DPZ=DE-DBLE(X)*DBLE(Y)*(DBLE(P(2,4))+DBLE(ABS(P(2,3))))
      ELSE
C...Formulae for massive scattered lepton.
        DE=DBLE(P(1,4))+(DBLE(ABS(P(2,3)))*(DBLE(Q2)+DBLE(P(1,5))**2+
     &  DBLE(P(4,5))**2)/(2.D0*DBLE(P(1,4)))-DBLE(PARL(22))/2.D0)/
     &  (DBLE(P(2,4))+DBLE(ABS(P(2,3))))
        DPZ=DBLE(P(1,4))-(DBLE(P(2,4))*(DBLE(Q2)+DBLE(P(1,5))**2+
     &  DBLE(P(4,5))**2)/(2.D0*DBLE(P(1,4)))+DBLE(PARL(22))/2.D0)/
     &  (DBLE(P(2,4))+DBLE(ABS(P(2,3))))
      ENDIF
      DPT=DE**2-DPZ**2-DBLE(P(4,5))**2
      IF(DPT.LT.0.D0) RETURN
      DPT=SQRT(DPT)
      P(4,1)=DPT
      P(4,2)=0.
      P(4,3)=DPZ
      P(4,4)=DE
      P(3,1)=-DPT
      P(3,2)=0.
      P(3,3)=DBLE(P(1,3))-DPZ
      P(3,4)=DBLE(P(1,4))-DE
      P(3,5)=-SQRT(Q2)
      K(3,3)=1
      K(4,3)=1
      N=4
      IF(L.EQ.3) GOTO 200

      IF(X.LT.XMIN.OR.X.GT.XMAX) RETURN
      IF(Y.LT.YMIN.OR.Y.GT.YMAX) RETURN
      IF(Q2.LT.Q2MIN.OR.Q2.GT.Q2MAX) RETURN
      IF(W2.LT.W2MIN.OR.W2.GT.W2MAX) RETURN
C-check: CUT(9),CUT(10) --> UMIN,UMAX needs change in /LINTRL/ --> next update 
      IF(U.LT.CUT(9).OR.U.GT.CUT(10)) THEN
      WRITE(*,*) 'Nu is out of range!'
      RETURN
      ENDIF
      IF(LST(17).EQ.0) THEN
        IF(P(4,4).LT.CUT(11).OR.P(4,4).GT.CUT(12))  RETURN
        THETAL=PLU(4,13)
C       THETAL=ACOS((P(1,1)*P(4,1)+P(1,2)*P(4,2)+P(1,3)*P(4,3))
C    &  /SQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2)/
C    &  SQRT(P(4,1)**2+P(4,2)**2+P(4,3)**2))
      ELSE
C...No cuts on energy, angle for initialisation of varying energy mode
        IF(LST(32).NE.0) GOTO 200
C...Transform scattered lepton back to lab system to make cut
C...in energy and angle (defined as space angle to incoming lepton).
        DO 110 J=1,5
        K(6,J)=K(4,J)
  110   P(6,J)=P(4,J)
        CALL LUDBRB(6,6,STHETA(1),SPHI(1),0.D0,0.D0,0.D0)
        CALL LUDBRB(6,6,0.,0.,DBETA(1,1),DBETA(1,2),DBETA(1,3))
        IF(P(6,4).LT.CUT(11).OR.P(6,4).GT.CUT(12))  RETURN
        THETAL=ACOS((PSAVE(3,1,1)*P(6,1)+PSAVE(3,1,2)*P(6,2)+
     &  PSAVE(3,1,3)*P(6,3))
     &  /SQRT(PSAVE(3,1,1)**2+PSAVE(3,1,2)**2+PSAVE(3,1,3)**2)/
     &  SQRT(P(6,1)**2+P(6,2)**2+P(6,3)**2))
      ENDIF
      IF(THETAL.LT.CUT(13).OR.THETAL.GT.CUT(14))  RETURN
  200 LKINEM=0
      RETURN
      END
C **********************************************************************

      SUBROUTINE LQCDPR(QG,QQB)

C...Probabilities for hard QCD events, qg or qqb, from integration of
C...QCD matrix elements event-by event or interpolation on x-W grid.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LGRID/ NXX,NWW,XX(31),WW(21),PQG(31,21,3),PQQB(31,21,2),
     &QGMAX(31,21,3),QQBMAX(31,21,2),YCUT(31,21),XTOT(31,21),NP
      DIMENSION PQSAVE(17)
      EXTERNAL DSIGMA,DSIGM2
      DATA NOUT,NABOVE/2*0/,NWARN/10/
      
      LOGICAL ZOOM

C...Get ycut from grid
      IF(LST(19).GE.0.OR.LST(19).EQ.-10) THEN
C...
C...qg and qqb event probabilities from interpolation on grid
      QG=0.
      QQB=0.
C...QCD weight zero for x->1 above grid and W small below grid
      IF(X.GT.XX(NXX).AND.X.GT.0.999) RETURN
      IF(LST(19).LT.10.AND.SQRT(W2).LT.WW(1).AND.WW(1).LT.6.) RETURN

      XP=X
C...Local variable W is W or y
      W=SQRT(W2)
      IF(LST(19).GE.10.OR.LST(19).EQ.-10) W=Y
      IF(X.LT.XX(1).OR.X.GT.XX(NXX).OR.
     &W.LT.WW(1).OR.W.GT.WW(NWW)) THEN
C...x and/or W/y outside limits of grid, write warning NWARN first times
        NOUT=NOUT+1
        IF(LST(3).GE.1.AND.NOUT.LE.NWARN)
     &  WRITE(6,1000) X,W,INT(PARI(29)),NWARN
        IF(X.LT.XX(1)) XP=XX(1)
        IF(X.GT.XX(NXX)) XP=XX(NXX)
        IF(W.LT.WW(1)) W=WW(1)
        IF(W.GT.WW(NWW)) W=WW(NWW)
      ENDIF

      IH=1
      IF(LST(30).EQ.1) IH=2
      IX=0
  100 IX=IX+1
      IF(XP.GT.XX(IX+1)) GOTO 100
      IW=0
  200 IW=IW+1
      IF(W.GT.WW(IW+1)) GOTO 200
      WD=(W-WW(IW))/(WW(IW+1)-WW(IW))
      XD=(XP-XX(IX))/(XX(IX+1)-XX(IX))

      DO 500 IP=1,NP
      X1P=(PQG(IX+1,IW,IP)-PQG(IX,IW,IP))*XD+PQG(IX,IW,IP)
      X2P=(PQG(IX+1,IW+1,IP)-PQG(IX,IW+1,IP))*XD+PQG(IX,IW+1,IP)
      QGIP=(X2P-X1P)*WD+X1P
      IF(NP.EQ.1) THEN
        QG=QGIP
        PARI(15)=MAX(QGMAX(IX,IW,IH),QGMAX(IX+1,IW+1,IH),
     &  QGMAX(IX+1,IW,IH),QGMAX(IX,IW+1,IH))
      ELSE
        QG=QG+PARI(23+IP)*QGIP
        PARI(14+IP)=MAX(QGMAX(IX,IW,IP),QGMAX(IX+1,IW+1,IP),
     &  QGMAX(IX+1,IW,IP),QGMAX(IX,IW+1,IP))
      ENDIF
      IF(IP.EQ.3) GOTO 500
      X1P=(PQQB(IX+1,IW,IP)-PQQB(IX,IW,IP))*XD+PQQB(IX,IW,IP)
      X2P=(PQQB(IX+1,IW+1,IP)-PQQB(IX,IW+1,IP))*XD+PQQB(IX,IW+1,IP)
      QQBIP=(X2P-X1P)*WD+X1P
      IF(NP.EQ.1) THEN
        QQB=QQBIP
        PARI(18)=MAX(QQBMAX(IX,IW,IH),QQBMAX(IX+1,IW+1,IH),
     &  QQBMAX(IX+1,IW,IH),QQBMAX(IX,IW+1,IH))
      ELSE
        QQB=QQB+PARI(23+IP)*QQBIP
        PARI(17+IP)=MAX(QQBMAX(IX,IW,IP),QQBMAX(IX+1,IW+1,IP),
     &  QQBMAX(IX+1,IW,IP),QQBMAX(IX,IW+1,IP))
      ENDIF
  500 CONTINUE

      IF(NP.NE.1) THEN
C...Get total x-section from interpolation to be used for normalization.
        X1P=(XTOT(IX+1,IW)-XTOT(IX,IW))*XD+XTOT(IX,IW)
        X2P=(XTOT(IX+1,IW+1)-XTOT(IX,IW+1))*XD+XTOT(IX,IW+1)
        PQ17=(X2P-X1P)*WD+X1P
        QG=QG/PQ17
        QQB=QQB/PQ17
      ENDIF

C..Interpolate in the grid
      X1P=(YCUT(IX+1,IW)-YCUT(IX,IW))*XD+YCUT(IX,IW)
      X2P=(YCUT(IX+1,IW+1)-YCUT(IX,IW+1))*XD+YCUT(IX,IW+1)
      PARL(27)=(X2P-X1P)*WD+X1P
C...Include alpha-strong in weight.
      QG=QG*PARL(25)
      QQB=QQB*PARL(25)
C...Get value of y-cut,
      IF(LST(19).GE.0) THEN
        IF(LST(33).EQ.-91) THEN
C...Include 3-jet cross section in denominator
          QTOT=1.+QG+QQB
          QG =QG/QTOT
          QQB=QQB/QTOT
        ENDIF
        IF(QG+QQB.GT.1) THEN
C...Sum of QCD event probabilities larger than unity, rescale to unity
C...and print warning for first NWARN cases.
          NABOVE=NABOVE+1
          IF(LST(3).GE.1.AND.NABOVE.LE.NWARN)
     &    WRITE(6,1100) QG,QQB,X,W,INT(PARI(29)),NWARN
          QGQQB=QG+QQB
          QG=QG/QGQQB
          QQB=QQB/QGQQB
        ENDIF
      ELSE
        IF(MAX(YCUT(IX,IW),YCUT(IX+1,IW+1),
     &         YCUT(IX+1,IW),YCUT(IX,IW+1))-
     &     MIN(YCUT(IX,IW),YCUT(IX+1,IW+1),
     &         YCUT(IX+1,IW),YCUT(IX,IW+1)).EQ.0.0) THEN
          RETURN
        ELSE
C...Get the minimum from the grid
          PARL(27)=MIN(YCUT(IX,IW),YCUT(IX+1,IW+1),
     &                 YCUT(IX+1,IW),YCUT(IX,IW+1))
        ENDIF
      ENDIF

C...Grid
      ENDIF
      
C...Calculate probabilities directly or refine value from grid
      IF(LST(19).LE.0) THEN
      
C...qg and qqbar event probabilities (and max values for simulation)
C...obtained by integrating QCD matrix elements for each event.
C     LST2=LST(2)
C     LST(2)=-3
C     NP=1
      LST(32)=1

      DO 1 I=1,17
    1 PQSAVE(I)=PQ(I)

      PARL(25)=ULALPS(Q2)
      PARI(20)=PQ(17)
      IF(LST(19).GT.-10) THEN
        IF(LST(20).LE.1) THEN
          PARL(27)=MAX(PARL(9)**2/W2,PARL(8))
          P27MAX=1.0
        ELSEIF(LST(20).EQ.2) THEN
          PARL(27)=MAX(PARL(9)**2/Q2,PARL(8))
          P27MAX=W2/Q2
        ELSEIF(LST(20).GE.3.AND.LST(20).LE.5) THEN
          PARL(27)=PARL(8)
          P27MAX=0.5
        ELSEIF(LST(20).EQ.6) THEN
          PARL(27)=PARL(9)
          P27MAX=W2
        ENDIF
      ELSE
        IF(LST(20).LE.1) THEN
          P27MAX=1.0
        ELSEIF(LST(20).EQ.2) THEN
          P27MAX=W2/Q2
        ELSEIF(LST(20).GE.3.AND.LST(20).LE.5) THEN
          P27MAX=0.5
        ELSEIF(LST(20).EQ.6) THEN
          P27MAX=W2
        ENDIF
      ENDIF

      ZOOM=.FALSE.
      IYCUT=0
      YCMIN=PARL(27)
      YCMAX=PARL(27)
   10 IYCUT=IYCUT+1
      RQG=0.
      RQQB=0.
CAE.Scheme for ME cutoff: W2, Q2, mixed
      IF(LST(20).LE.1) THEN
        XPMIN=DBLE(X)/(1.D0-2.D0*(1.D0-DBLE(X))*DBLE(PARL(27)))
        XPMAX=DBLE(X)/(DBLE(X)+(1.D0-DBLE(X))*DBLE(PARL(27)))
      ELSEIF(LST(20).EQ.2) THEN
        XPMIN=DBLE(X)/(1.D0-2.D0*DBLE(X)*DBLE(PARL(27)))
        XPMAX=1.D0/(1.D0+DBLE(PARL(27)))
      ELSEIF(LST(20).EQ.3.OR.LST(20).EQ.4) THEN
        XPMIN=X
        XPMAX=1./(1.+PARL(9))
      ELSEIF(LST(20).EQ.5) THEN
        XPMIN=X
        XPMAX=Q2/(Q2+PARL(9))
      ELSEIF(LST(20).EQ.6) THEN
        XPMIN=X
        XPMAX=Q2/(Q2+PARL(27))
      ELSE
        WRITE(6,*) 'LQCDPR: No such jet scheme!'
      ENDIF
CAE
      IF(XPMIN.LT.X.OR.XPMIN.GT.1.) GOTO 40
      IF(XPMIN.GE.XPMAX) GOTO 40

      PARI(15)=0.
      PARI(16)=0.
      PARI(18)=0.
      PARI(19)=0.
C...QCD-Compton -> qg-event
      LST(24)=2
      EPS=PARL(11)
CAE      CALL GADAP(XPMIN,XPMAX,DSIGMA,EPS,RQG)
      CALL GADAP(LOG(1.0-XPMAX),LOG(1.0-XPMIN),DSIGM2,EPS,RQG)
C...QCD-fusion  -> qq-event
      LST(24)=3
      EPS=PARL(11)
CAE      CALL GADAP(XPMIN,XPMAX,DSIGMA,EPS,RQQB)
      CALL GADAP(LOG(1.0-XPMAX),LOG(1.0-XPMIN),DSIGM2,EPS,RQQB)
C...q-event
      RQ=1.-RQG-RQQB
CAE      WRITE(6,*) IYCUT,RQ,PARL(27),YCMIN,YCMAX
      IF(.NOT.ZOOM) THEN
CAE.First find interval so that RQ>0
        IF(RQ.LT.0.AND.IYCUT.LT.10) THEN
          PARL(27)=MIN(1.1*EXP(-2.0*RQ)*PARL(27),P27MAX)
          YCMIN=YCMAX
          YCMAX=PARL(27)
        ELSEIF(RQ.LT.0.AND.IYCUT.GE.10) THEN
C...Terminate procedure after some iterations
          RTOT=(RQG+RQQB)*1.05
          RQG=RQG/RTOT
          RQQB=RQQB/RTOT
          RQ=1.-RQG-RQQB
C          IF(LST(3).GE.1) THEN
C            WRITE(6,*) 'Warning! sigma>tot for x,q2,cut=',X,Q2,PARL(27)
C            WRITE(6,*) 'Weights set to=',RQ,RQG,RQQB
C          ENDIF
C...Break loop
          GOTO 40
        ELSEIF(IYCUT.GE.2.AND.RQ.GT.PARL(13)) THEN
C...If RQ>PARL(13), then ycut was increased to much
          ZOOM=.TRUE.
          PARL(27)=(YCMIN+YCMAX)/2.
        ELSE
C...correct ycut found
          GOTO 40
        ENDIF
      ELSE
C...Zoom in on ycut so that 0<RQ<PARL(13)
        IF(RQ.LT.0.AND.IYCUT.LT.40) THEN
          YCMIN=PARL(27)
          PARL(27)=(YCMIN+YCMAX)/2.
        ELSEIF(RQ.GT.PARL(13).AND.IYCUT.LT.40) THEN
          YCMAX=PARL(27)
          PARL(27)=(YCMIN+YCMAX)/2.
C...Catch infinite loop
        ELSEIF(IYCUT.GE.40) THEN
          IF(LST(3).GE.1) THEN
            WRITE(6,*) 'LQCDPR: Warning, PARL(27) not found.'
          ENDIF
          RTOT=(RQG+RQQB)*1.05
          RQG=RQG/MAX(1.0,RTOT)
          RQQB=RQQB/MAX(1.0,RTOT)
          RQ=1.-RQG-RQQB
C...Break loop
          GOTO 40
        ELSE
C...ycut found, break loop
          GOTO 40
        ENDIF
      ENDIF
C...Loop until correct weights found
      GOTO 10

   40 CONTINUE
CAE
      IF(LST(33).EQ.-91) THEN
C...Include 3-jet cross section in denominator
        QTOT=1.+RQG+RQQB
        RQG =RQG/QTOT
        RQQB=RQQB/QTOT
        RQ=1.-RQG-RQQB
      ENDIF

C     LST(2)=LST2
      LST(32)=0
      DO 90 I=1,17
   90 PQ(I)=PQSAVE(I)
      QG=RQG
      QQB=RQQB

C...Refine      
      ENDIF
      
 1000 FORMAT(' Warning: x=',F7.4,' or W/y=',F10.4,' outside QCD grid',
     &' in event no.',I8,/,10X,
     &'weight on limit of grid used. Only first',I5,' warnings printed')
 1100 FORMAT(' Warning: Sum of QCD probabilities larger than unity ',
     &' QG, QQB =',2F8.4,/10X,'occurs at x, W/y =',2F10.4,
     &' in event no.',I8,/,10X,
     &'Weights rescaled to unit sum. Only first',I5,' warnings printed')
      RETURN
      END

C **********************************************************************

      SUBROUTINE LQEV

C...Generate an ordinary 2-jet event, q-event.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)

      LST(24)=1
      W=SQRT(W2)

C...Choose flavour of scattered quark and target remnant.
  200 CALL LFLAV(IFL,IFLR)
      IF(LST(21).NE.0) GOTO 200

      GOTO 210
C...Entry used for Ariadne
      ENTRY LQEVAR(IFLAR,IFLRAR)
      IFL=IFLAR
      IFLR=IFLRAR
      LST(24)=1
      W=SQRT(W2)

  210 CONTINUE
      MSTJ(93)=1
      AMIFL=ULMASS(IFL)
      MSTJ(93)=1
      AMIFLR=ULMASS(IFLR)
      IF(LST(14).EQ.0.OR.IFLR.GT.10
     &.OR.(LST(8).GE.2.AND.MOD(LST(8),10).NE.9)) THEN
C...Check if energy in jet system is enough for fragmentation.
        IF(W.LT.AMIFL+AMIFLR+PARJ(32)) GOTO 200
        CALL LU2ENT(MSTU(1),IFL,IFLR,W)
        K(MSTU(1)+1,3)=2
      ELSE
C...Target remnant is not a simple diquark, special treatment needed.
        IF(W.LT.AMIFL+AMIFLR+0.9+PARJ(32)) GOTO 200
        IFLRO=IFLR
        NREMH=0
  300   NREMH=NREMH+1
        IF(NREMH.GT.100) GOTO 999
C...Give balancing pt to IFLQ and IFLQQ.
        CALL LPRIKT(PARL(14),PT,PHI)
        CALL LREMH(IFLRO,PT,IFLR,K2,XT)
        MSTJ(93)=1
        AMIFLR=ULMASS(IFLR)
CJR--
        KFIFLR=LUCOMP(IFLR)
        IF (KFIFLR.EQ.90) THEN
          AMIFLR=AMIFLR-2.*PARL(20)
        ELSEIF (1.LE.KFIFLR .AND. KFIFLR.LE.6) THEN
          AMIFLR=AMIFLR-PARL(20)
        ENDIF
        MSTJ(93)=1
        AMK2=ULMASS(K2)
        KFK2=LUCOMP(K2)
        IF (KFK2.EQ.90) THEN
          AMK2=AMK2-2.*PARL(20)
        ELSEIF (1.LE.KFK2 .AND. KFK2.LE.6) THEN
          AMK2=AMK2-PARL(20)
        ENDIF
CJR--
        PT2=PT**2
        TM2K2=AMK2**2+PT2
        EK2=.5*(XT*W+TM2K2/XT/W)
        PZK2=-.5*(XT*W-TM2K2/XT/W)
        EPZ=W-TM2K2/XT/W
        WT=(1.-XT)*W*EPZ-PT2
C...Check if energy in jet system is enough for fragmentation.
        IF(WT.LT.(AMIFL+AMIFLR+PARJ(32))**2) GOTO 300
        WT=SQRT(WT+PT2)
        TMIFLR=AMIFLR**2+PT2
        EIFL=.5*(WT+(AMIFL**2-TMIFLR)/WT)
        EIFLR=.5*(WT+(-AMIFL**2+TMIFLR)/WT)
        THER=ULANGL(-SQRT(EIFLR**2-TMIFLR),PT)
C...Form jet system.
        CALL LU1ENT(-MSTU(1),IFL,EIFL,0.,0.)
        CALL LU1ENT(MSTU(1)+1,IFLR,EIFLR,THER,PHI)
        CALL LUDBRB(MSTU(1),0,0.,0.,0.D0,0.D0,
     &  (DBLE(EPZ)-(1.D0-DBLE(XT))*DBLE(W))/
     &  (DBLE(EPZ)+(1.D0-DBLE(XT))*DBLE(W)))
        THEK2=ULANGL(PZK2,PT)
C...Add formed "target" particle.
        MSTU(10)=1
        P(MSTU(1)+2,5)=AMK2
        CALL LU1ENT(MSTU(1)+2,K2,EK2,THEK2,PHI+3.1415927)
        MSTU(10)=2
        K(MSTU(1)+1,3)=2
        K(MSTU(1)+2,3)=2
CIC...Target remnants required to go backwards in hadronic cms
        IF(P(MSTU(1)+1,3).GT.0..OR.P(MSTU(1)+2,3).GT.0.) GOTO 300
      ENDIF
      
CAE...Set reasonable values to the ME variables xp,zq and phi
      PARL(28)=1.0
      PARL(29)=1.0
      PARL(30)=0.0

      LST(21)=0
      RETURN

  999 LST(21)=3
      RETURN
      END

C **********************************************************************

      SUBROUTINE LQGEV

C...Generate quark-gluon jet event, choose xp and zp according to QCD
C...matrix elements and apply cuts for soft and collinear gluons.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)

      LST(24)=2
      W=SQRT(W2)
      J1=MSTU(1)
      J2=MSTU(1)+1
      J3=MSTU(1)+2
      J4=MSTU(1)+3

      CALL LXP(XP,IFAIL)
      IF(IFAIL.NE.0) GOTO 999

C...Choose flavour of scattered quark and target remnant.
  200 CALL LFLAV(IFL,IFLR)
      IF(LST(21).NE.0) RETURN
      CALL LZP(XP,ZP,IFAIL)
      IF(IFAIL.NE.0) GOTO 999
      MSTJ(93)=1
      AMIFL=ULMASS(IFL)
      MSTJ(93)=1
      AMIFLR=ULMASS(IFLR)
      
      IF(LST(14).EQ.0.OR.IFLR.GT.10
     &.OR.(LST(8).GE.2.AND.MOD(LST(8),10).NE.9)) THEN
        IF(W.LT.AMIFL+AMIFLR+PARJ(32)) GOTO 999
        IF(LQMCUT(XP,ZP,AMIFL,0.,AMIFLR).NE.0) GOTO 999
        CALL LU3ENT(J1,IFL,21,IFLR,W,PARI(21),PARI(23))
        K(MSTU(1)+2,3)=2
        CALL LUROBO(ACOS(-P(J3,3)/SQRT(P(J3,3)**2+P(J3,1)**2)),
     &  0.,0.,0.,0.)
      ELSE
C...Target remnant is not a simple diquark, special treatment needed.
        IF(W.LT.AMIFL+AMIFLR+1.+PARJ(32)) GOTO 999
        IF(LQMCUT(XP,ZP,AMIFL,0.,1.).NE.0) GOTO 999
        IFLRO=IFLR
        NREMH=0
  300   NREMH=NREMH+1
        IF(NREMH.GT.100) GOTO 999
        CALL LPRIKT(PARL(14),PT,PHI)
        CALL LREMH(IFLRO,PT,IFLR,K2,XT)
        MSTJ(93)=1
        AMIFLR=ULMASS(IFLR)
CJR--
        KFIFLR=LUCOMP(IFLR)
        IF (KFIFLR.EQ.90) THEN
          AMIFLR=AMIFLR-2.*PARL(20)
        ELSEIF (1.LE.KFIFLR .AND. KFIFLR.LE.6) THEN
          AMIFLR=AMIFLR-PARL(20)
        ENDIF
        MSTJ(93)=1
        AMK2=ULMASS(K2)
        KFK2=LUCOMP(K2)
        IF (KFK2.EQ.90) THEN
          AMK2=AMK2-2.*PARL(20)
        ELSEIF (1.LE.KFK2 .AND. KFK2.LE.6) THEN
          AMK2=AMK2-PARL(20)
        ENDIF
CJR--
        P(J1,5)=AMIFL
        P(J2,5)=0.
        PT2=PT**2
        TM2K2=AMK2**2+PT2
        TMIFLR=AMIFLR**2+PT2
        P(J3,5)=SQRT(TM2K2/XT+TMIFLR/(1.-XT))
        IF(LQMCUT(XP,ZP,AMIFL,0.,P(J3,5)).NE.0) GOTO 300
        MSTU(10)=1
        CALL LU3ENT(J1,IFL,21,IFLR,W,PARI(21),PARI(23))
        K(MSTU(1)+2,3)=2
        MSTU(10)=2
        CALL LUROBO(ACOS(-P(J3,3)/SQRT(P(J3,3)**2+P(J3,1)**2)),
     &  0.,0.,0.,0.)
        EPZ=P(J3,4)-P(J3,3)
        P(J3,1)=PT*COS(PHI)
        P(J3,2)=PT*SIN(PHI)
        P(J3,3)=-0.5*((1.-XT)*EPZ-TMIFLR/(1.-XT)/EPZ)
        P(J3,4)= 0.5*((1.-XT)*EPZ+TMIFLR/(1.-XT)/EPZ)
        P(J3,5)=AMIFLR
        P(J4,1)=-P(J3,1)
        P(J4,2)=-P(J3,2)
        P(J4,3)=-0.5*(XT*EPZ-TM2K2/XT/EPZ)
        P(J4,4)= 0.5*(XT*EPZ+TM2K2/XT/EPZ)
        P(J4,5)=AMK2
        K(J4,1)=1
        K(J4,2)=K2
        K(J4,3)=2
        K(J4,4)=0
        K(J4,5)=0
        N=J4
        IF((P(J3,4)+P(J2,4)/2.)**2-(P(J3,1)+P(J2,1)/2.)**2-P(J3,2)**2
     &  -(P(J3,3)+P(J2,3)/2.)**2.LT.(AMIFLR+2.5*PARJ(32))**2) GOTO 300
      ENDIF

      CALL LAZIMU(XP,ZP)
      LST(21)=0
      RETURN

  999 LST(21)=4
      RETURN
      END

C **********************************************************************

      SUBROUTINE LQQBEV

C...Generate boson-gluon fusion event, choose xp and zp according to
C...QCD matrix elements and apply cuts for softness and collinearness.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)

      LST(24)=3
      W=SQRT(W2)
      J1=MSTU(1)
      J2=MSTU(1)+1
      J3=MSTU(1)+2
      J4=MSTU(1)+3

      CALL LXP(XP,IFAIL)
      IF(IFAIL.NE.0) GOTO 999

C...Choose flavour of produced quark-antiquark pair.
  200 CALL LFLAV(IFL1,IFL3)
      IF(LST(21).NE.0) RETURN
      IF(IFL1.LT.0) THEN
C...Put quark in first position
        IFL1S=IFL1
        IFL1=IFL3
        IFL3=IFL1S
      ENDIF
      CALL LZP(XP,ZP,IFAIL)
      IF(IFAIL.NE.0) GOTO 999
      IFL1A=IABS(IFL1)
      IFL3A=IABS(IFL3)
      MSTJ(93)=1
      AMIFL1=ULMASS(IFL1)
      MSTJ(93)=1
      AMIFL3=ULMASS(IFL3)

      IF(LST(14).EQ.0.OR.(LST(8).GE.2.AND.MOD(LST(8),10).NE.9)) THEN
C...If baryon production from target remnant is neglected the
C...target remnant is approximated by a gluon.
        IF(W.LT.AMIFL1+AMIFL3+PARJ(32)) GOTO 999
        IF(LQMCUT(XP,ZP,AMIFL1,0.,AMIFL3).NE.0) GOTO 999
        CALL LU3ENT(J1,IFL1,21,IFL3,W,PARI(21),PARI(23))
        K(MSTU(1)+1,3)=2
C...Align target remnant (gluon) along -z axis
        CALL LUROBO(-ACOS(-P(J2,3)/SQRT(P(J2,3)**2+P(J2,1)**2)),
     &  0.,0.,0.,0.)
C...Phi-rotation to bring quark to phi=0.
        CALL LUROBO(0.,-PLU(J1,15),0.,0.,0.)
      ELSE

        IF(W.LT.AMIFL1+AMIFL3+0.9+2.*PARJ(32)) GOTO 999
        IF(LQMCUT(XP,ZP,AMIFL1,1.,AMIFL3).NE.0) GOTO 999
        P(J1,5)=AMIFL1
        P(J3,5)=AMIFL3
C...Choose target valence quark/diquark to form jet system with
C...produced antiquark/quark.
        IFLR2=INT(1.+LST(22)/3.+RLU(0))
        IF(IFLR2.EQ.LST(22)) THEN
          IFLR1=2101
          IF(RLU(0).GT.PARL(4)) IFLR1=2103
        ELSE
          IFLR1=1000*IFLR2+100*IFLR2+3
        ENDIF
        IFLR2=3-IFLR2
        MSTJ(93)=1
        AMR1=ULMASS(IFLR1)
CJR--
        KFIFL1=LUCOMP(IFLR1)
        IF (KFIFL1.EQ.90) THEN
          AMR1=AMR1-2*PARL(20)
        ELSEIF (1.LE.KFIFL1 .AND. KFIFL1.LE.6) THEN
          AMR1=AMR1-PARL(20)
        ENDIF
        MSTJ(93)=1
        AMR2=ULMASS(IFLR2)
        KFIFL2=LUCOMP(IFLR2)
        IF (KFIFL2.EQ.90) THEN
          AMR2=AMR2-2.*PARL(20)
        ELSEIF (1.LE.KFIFL2 .AND. KFIFL2.LE.6) THEN
          AMR2=AMR2-PARL(20)
        ENDIF
CJR--
        NREMH=0
  310   NREMH=NREMH+1
        IF(NREMH.GT.100) GOTO 999
        CALL LPRIKT(PARL(14),PT,PHI)
        CALL LREMH(0,PT,IFLR1,IFLR2,XT)
        PT2=PT**2
        TM2R1=AMR1**2+PT2
        TM2R2=AMR2**2+PT2
        P(J2,5)=SQRT(TM2R1/(1.-XT)+TM2R2/XT)
        IF(LQMCUT(XP,ZP,AMIFL1,P(J2,5),AMIFL3).NE.0) GOTO 310
        MSTU(10)=1
        CALL LU3ENT(J1,IFL1,21,IFL3,W,PARI(21),PARI(23))
        MSTU(10)=2
C...Align target remnant (effective gluon) along -z axis
        CALL LUROBO(-ACOS(-P(J2,3)/SQRT(P(J2,3)**2+P(J2,1)**2)),
     &0.,0.,0.,0.)
C...Phi-rotation to bring quark to phi=0.
        CALL LUROBO(0.,-PLU(J1,15),0.,0.,0.)
        EPZ=P(J2,4)-P(J2,3)
        IF(IFL1.GT.0) THEN
          IR1=J2
          IR2=J4
        ELSE
          IR1=J4
          IR2=J2
        ENDIF
        P(IR1,1)=PT*COS(PHI)
        P(IR1,2)=PT*SIN(PHI)
        P(IR1,3)=-0.5*((1.-XT)*EPZ-TM2R1/(1.-XT)/EPZ)
        P(IR1,4)= 0.5*((1.-XT)*EPZ+TM2R1/(1.-XT)/EPZ)
        P(IR1,5)=AMR1
        P(IR2,1)=-P(IR1,1)
        P(IR2,2)=-P(IR1,2)
        P(IR2,3)=-0.5*(XT*EPZ-TM2R2/XT/EPZ)
        P(IR2,4)= 0.5*(XT*EPZ+TM2R2/XT/EPZ)
        P(IR2,5)=AMR2
        K(IR1,1)=1
        K(IR1,2)=IFLR1
        K(IR2,1)=1
        K(IR2,2)=IFLR2
        K(J3,1)=2
        DO 320 I=J1,J4
          DO 320 J=3,5
  320       K(I,J)=0
        N=J4
        K(IR1,3)=2
        K(IR2,3)=2
       IF((P(J1,4)+P(J2,4))**2-(P(J1,1)+P(J2,1))**2-(P(J1,3)+P(J2,3))**2
     &  -P(J2,2)**2.LT.(P(J1,5)+P(J2,5)+PARJ(32))**2) GOTO 310
       IF((P(J3,4)+P(J4,4))**2-(P(J3,1)+P(J4,1))**2-(P(J3,3)+P(J4,3))**2
     &  -P(J4,2)**2.LT.(P(J3,5)+P(J4,5)+PARJ(32))**2) GOTO 310

      ENDIF

      CALL LAZIMU(XP,ZP)
      LST(21)=0
      RETURN

  999 LST(21)=5
      RETURN
      END

C **********************************************************************

      SUBROUTINE LXP(XP,IFAIL)

C...Choose value of XP according to QCD matrix elements weighted by
C...structure functions.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      DOUBLE PRECISION DXPMAX

      IFAIL=1
CAE.Scheme for ME cutoff: W2, Q2, mixed, z-shat
      IF(LST(20).LE.1) THEN
        XPMIN=DBLE(X)/(1.D0-2.D0*(1.D0-DBLE(X))*DBLE(PARL(27)))
        DXPMAX=DBLE(X)/(DBLE(X)+(1.D0-DBLE(X))*DBLE(PARL(27)))
        XPMAX=SNGL(DXPMAX)
      ELSEIF(LST(20).EQ.2) THEN
        XPMIN=DBLE(X)/(1.D0-2.D0*DBLE(X)*DBLE(PARL(27)))
        DXPMAX=1.D0/(1.D0+DBLE(PARL(27)))
        XPMAX=SNGL(DXPMAX)
      ELSEIF(LST(20).EQ.3.OR.LST(20).EQ.4) THEN
        XPMIN=X
        DXPMAX=1.D0/(1.D0+DBLE(PARL(9)))
        XPMAX=DXPMAX
      ELSEIF(LST(20).EQ.5) THEN
        XPMIN=X
        DXPMAX=DBLE(Q2)/(DBLE(Q2)+DBLE(PARL(9)))
        XPMAX=DXPMAX
      ELSEIF(LST(20).EQ.6) THEN
        XPMIN=X
        DXPMAX=DBLE(Q2)/(DBLE(Q2)+DBLE(PARL(27)))
        XPMAX=DXPMAX
      ELSE
        WRITE(6,*) 'LXP: No such jet scheme!'
      ENDIF
CAE
      IF(XPMIN.LT.X.OR.XPMIN.GT.1.) RETURN
      IF(XPMIN.GE.XPMAX) RETURN
      AP=1.-XPMIN
      BP=(1.D0-DXPMAX)/AP
      IF(LST(24).EQ.2) THEN
        QXPMAX=PARI(15)
        IF(LST(17).NE.0.AND.LST(19).GE.0) QXPMAX=
     &  PARI(24)*PARI(15)+PARI(25)*PARI(16)+PARI(26)*PARI(17)
      ELSE
        QXPMAX=PARI(18)
        IF(LST(17).NE.0.AND.LST(19).GE.0) QXPMAX=
     &  PARI(24)*PARI(18)+PARI(25)*PARI(19)
      ENDIF
C...Safety factor on max value.
      QXPMAX=QXPMAX*1.05
      LOOP=0
  100 LOOP=LOOP+1
      IF(LOOP.GT.1000) RETURN
      XP=1.-AP*BP**RLU(0)
      XPWEIT=DSIGMA(XP)/QXPMAX
      IF(XPWEIT.LT.RLU(0)) GOTO 100
      IFAIL=0
      RETURN
      END

C **********************************************************************

      SUBROUTINE LZP(XP,ZP,IFAIL)

C...Choose value of ZP according to QCD matrix elements.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      DATA C1,C2/0.2122066,0.0795775/,DZPMAX,SZP,CP/3*0./
      FQG(DZ,DX,DA,DB,DC)=DA*(DZ**2+DX**2)/(1.-DX)+2.*DA*DX*DZ*(1.-DZ)
     &+2.*DA*(1.-DZ)+4.*DB*DX*DZ*(1.-DZ)+DC*(DZ**2+DX**2)/(1.-DX)+
     &2.*DC*(DX+DZ)*(1.-DZ)
      FQQ(DZ,DX,DA,DB,DC,DD,DE)=DA*DD*(DZ**2+(1.-DZ)**2)+DB*DE*DZ*
     &(1.-DZ)+DC*DD*(2.*DZ-1.)

      IFAIL=1
      IH=1
      IF(LST(30).EQ.1) IH=2
CAE.Scheme for ME cutoff: W2, Q2, mixed, z-shat
      IF(LST(20).LE.1) THEN
        ZPMIN=(1.-X)*XP/(XP-X)*PARL(27)
      ELSEIF(LST(20).EQ.2) THEN
        ZPMIN=X*XP/(XP-X)*PARL(27)
      ELSEIF(LST(20).GE.3.AND.LST(20).LE.5) THEN
        ZPMIN=PARL(27)
      ELSEIF(LST(20).EQ.6) THEN
        ZPMIN=PARL(8)
      ELSE
        WRITE(6,*) 'LZP: No such jet scheme!'
      ENDIF
CAE
      IF(ZPMIN.LE.0..OR.ZPMIN.GE.0.5) RETURN
      ZPMAX=1.-ZPMIN
      I=IABS(LST(25))
      AP=1.-ZPMIN
      BP=ZPMIN/AP
      IF(LST(23).EQ.2) THEN
        A=PARI(24)
        B=PARI(25)
        CSIGN=-LST(30)*ISIGN(1,LST(25))*PARI(26)
      ELSE
        A=(EWQC(1,IH,I)+EWQC(2,IH,I))*PARI(24)
        B=(EWQC(1,IH,I)+EWQC(2,IH,I))*PARI(25)
        C=(EWQC(1,IH,I)-EWQC(2,IH,I))*PARI(26)
        CSIGN=-C*LST(30)*ISIGN(1,LST(25))
      ENDIF
      IF(LST(24).EQ.2) THEN
        DZPMAX=MAX(FQG(ZPMIN,XP,A,B,CSIGN),FQG(ZPMAX,XP,A,B,CSIGN))
        AA=2.*(A+CSIGN)/(1.-XP)-4.*A*XP-8.*B*XP-4.*CSIGN
        IF(ABS(AA).GT.1.E-20) THEN
          BB=2.*A*(XP-1.)+4.*B*XP+2.*CSIGN*(1.-XP)
          Z1=-BB/AA
          IF(Z1.GT.ZPMIN.AND.Z1.LT.ZPMAX)
     &    DZPMAX=MAX(DZPMAX,FQG(Z1,XP,A,B,CSIGN))
        ENDIF
        DZPMAX=DZPMAX*C1*1.05
      ELSEIF(LST(24).EQ.3) THEN
        CP=1./BP**2
        D=XP**2+(1.-XP)**2
        E=8.*XP*(1-XP)
        DZPMAX=MAX(FQQ(ZPMIN,XP,A,B,CSIGN,D,E),
     &  FQQ(ZPMAX,XP,A,B,CSIGN,D,E))
        AA=4.*A*D-2.*B*E
        IF(ABS(AA).GT.1.E-20) THEN
          BB=B*E-2.*A*D+2.*CSIGN*D
          Z1=-BB/AA
          IF(Z1.GT.ZPMIN.AND.Z1.LT.ZPMAX)
     &    DZPMAX=MAX(DZPMAX,FQQ(Z1,XP,A,B,CSIGN,D,E))
        ENDIF
        DZPMAX=DZPMAX*C2*1.05
      ENDIF
      IPART=LST(24)-1
      LOOP=0
  100 LOOP=LOOP+1
      IF(LOOP.GT.1000) RETURN
      IF(LST(24).EQ.2) THEN
        ZP=1.-AP*BP**RLU(0)
        SZP=1.-ZP
      ELSEIF(LST(24).EQ.3) THEN
        DP=BP*CP**RLU(0)
        ZP=DP/(1.+DP)
        SZP=ZP*(1.-ZP)
      ENDIF
      ZPWEIT=SZP*(A*DQCD(0,IPART,1,XP,ZP,0.)+B*DQCD(0,IPART,2,XP,ZP,0.)
     &+CSIGN*DQCD(0,IPART,3,XP,ZP,0.))/DZPMAX
     
      IF(ZPWEIT.LT.RLU(0)) GOTO 100
      IFAIL=0
      RETURN
      END

C **********************************************************************

      FUNCTION LQMCUT(XP,ZP,AM1,AM2,AM3)

C...Apply cuts, if necessary, on the event configuration
C...obtained from QCD matrix elements.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      DATA S12,S23,S13/3*0./

      IF(LST(24).EQ.2) THEN
        S12=Q2*(1.-XP)/XP
        S23=Q2*(XP-X)*(1.-ZP)/X/XP+AM2**2+AM3**2
        S13=Q2*(XP-X)*ZP/X/XP+AM1**2+AM3**2
      ELSEIF(LST(24).EQ.3) THEN
        S13=Q2*(1.-XP)/XP
        S23=Q2*(XP-X)*(1.-ZP)/X/XP+AM2**2+AM3**2
        S12=Q2*(XP-X)*ZP/X/XP+AM1**2+AM2**2
        IF(S13.LT.(AM1+AM3)**2) GOTO 900
      ENDIF

      W=SQRT(W2)
      X1=1.-(S23-AM1**2)/W2
      X3=1.-(S12-AM3**2)/W2
      X2=2.-X1-X3
      PARI(21)=X1
      PARI(22)=X2
      PARI(23)=X3
      
           
      IF(X1.GT.1..OR.X2.GT.1..OR.X3.GT.1.) GOTO 900
      IF(X1*W/2..LT.AM1.OR.X2*W/2..LT.AM2.OR.X3*W/2..LT.AM3) GOTO 900
      PA1=SQRT((0.5*X1*W)**2-AM1**2)
      PA2=SQRT((0.5*X2*W)**2-AM2**2)
      PA3=SQRT((0.5*X3*W)**2-AM3**2)
      IF(ABS((PA3**2-PA1**2-PA2**2)/(2.*PA1*PA2)).GE.1.) GOTO 900
      IF(ABS((PA2**2-PA1**2-PA3**2)/(2.*PA1*PA3)).GE.1.) GOTO 900
      LQMCUT=0
      RETURN

  900 LQMCUT=1
      RETURN
      END

C **********************************************************************

      SUBROUTINE LAZIMU(XP,ZP)

C...Choose azimuthal angle (PHI) according to QCD matrix elements.
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)

      J=LST(24)-1
      SGN=SIGN(1.,2.5-LST(24))
      IFL=LST(25)
      I=IABS(IFL)
      IH=1
      IF(LST(30).EQ.1) IH=2

      IF(LST(23).EQ.2) THEN
        A=PARI(24)*DQCD(0,J,1,XP,ZP,Y)+PARI(25)*DQCD(0,J,2,XP,ZP,Y)
     &  -LST(30)*ISIGN(1,IFL)*PARI(26)*DQCD(0,J,3,XP,ZP,Y)
        B=DQCD(1,J,1,XP,ZP,Y)
     &  +SGN*LST(30)*ISIGN(1,IFL)*DQCD(1,J,3,XP,ZP,Y)
        C=DQCD(2,J,1,XP,ZP,Y)
      ELSE
        A=(EWQC(1,IH,I)+EWQC(2,IH,I))*(PARI(24)*DQCD(0,J,1,XP,ZP,Y)+
     &    PARI(25)*DQCD(0,J,2,XP,ZP,Y))
     &    -LST(30)*ISIGN(1,IFL)*(EWQC(1,IH,I)-EWQC(2,IH,I))
     &    *PARI(26)*DQCD(0,J,3,XP,ZP,Y)
        B=(EWQC(1,IH,I)+EWQC(2,IH,I))*DQCD(1,J,1,XP,ZP,Y)
     &    +SGN*LST(30)*ISIGN(1,IFL)*(EWQC(1,IH,I)-EWQC(2,IH,I))
     &    *DQCD(1,J,3,XP,ZP,Y)
        C=(EWQC(1,IH,I)+EWQC(2,IH,I))*DQCD(2,J,1,XP,ZP,Y)
      ENDIF

      PHIMAX=ABS(A)+ABS(B)+ABS(C)
  100 PHI=6.2832*RLU(0)
      IF(A+B*COS(PHI)+C*COS(2.*PHI).LT.RLU(0)*PHIMAX) GOTO 100
      CALL LUROBO(0.,PHI,0.,0.,0.)

CAE.Store ME variables
      PARL(28)=XP
      PARL(29)=ZP
      PARL(30)=PHI
CAE
      RETURN
      END

C **********************************************************************

      FUNCTION DSIGMA(XP)

C...Differential cross section for first order QCD processes.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      DIMENSION XPQ(-6:6),PQH(17,2)

      DSIGMA=0.
      DO 10 I=1,17
        PQH(I,1)=0.
        PQH(I,2)=0.
        PQ(I)=0.
 10   CONTINUE
      MSTJ(93)=1
      AMU=ULMASS(2)
      IF(LST(32).EQ.1.AND.LST(19).GE.0.AND.LST(17).EQ.1) THEN
        IL=LST(32)
        IU=LST(32)
      ELSE
        IL=1
        IU=3
        IF(LST(23).EQ.1.OR.LST(24).EQ.3) IU=2
      ENDIF
      XI=X/XP
C...Scheme for ME cutoff: W2, Q2, mixed
      IF(LST(20).LE.1) THEN
        ZPMIN=(1.-X)*XP/(XP-X)*PARL(27)
      ELSEIF(LST(20).EQ.2) THEN
        ZPMIN=X*XP/(XP-X)*PARL(27)
      ELSEIF(LST(20).GE.3.AND.LST(20).LE.5) THEN
        ZPMIN=PARL(27)
      ELSEIF(LST(20).GE.6) THEN
        ZPMIN=PARL(8)
      ENDIF

      IF(ZPMIN.LE.0..OR.ZPMIN.GE.0.5) RETURN
      ZPMAX=1.D0-DBLE(ZPMIN)
      CALL LNSTRF(XI,Q2,XPQ)
      IF(LST(24).EQ.3) GOTO 3000

C...Gluon bremsstrahlung process, i.e. qg-event.
 2000 DO 2500 IP=IL,IU
      SIG=DQCDI(1,IP,XP,ZPMIN,ZPMAX)
      SGN=SIGN(1.,5.-2.*IP)
      DO 2300 IH=1,2
      IF(IH.EQ.1) THEN
        IF(PARL(6).GT.0.99) GOTO 2300
        IF(LST(32).EQ.0.AND.LST(30).NE.-1) GOTO 2300
      ELSEIF(IH.EQ.2) THEN
        IF(PARL(6).LT.-0.99) GOTO 2300
        IF(LST(32).EQ.0.AND.LST(30).NE.1) GOTO 2300
      ENDIF
      IF(LST(32).NE.0) LST(30)=SIGN(1.,IH-1.5)
      IF(LST(23).NE.2) THEN
        DO 2100 I=1,LST(12)
        WQ=XPQ(I)*SIG*(EWQC(1,IH,I)+SGN*EWQC(2,IH,I))
        WQB=XPQ(-I)*SIG*SGN*(EWQC(1,IH,I)+SGN*EWQC(2,IH,I))
C...Include y-dependence.
        WQ=WQ*PARI(23+IP)
        WQB=WQB*PARI(23+IP)
        PQH(I,IH)=PQH(I,IH)+WQ
        PQH(I+LST(12),IH)=PQH(I+LST(12),IH)+WQB
        PQH(17,IH)=PQH(17,IH)+WQ+WQB
 2100   CONTINUE
      ELSEIF(LST(23).EQ.2) THEN
C...Zero CC cross-section for one helicity state.
        IF(KSAVE(1).LT.0.AND.IH.EQ.1
     &  .OR.KSAVE(1).GT.0.AND.IH.EQ.2) GOTO 2300
        IF(IP.EQ.3) THEN
          TQ=-LST(30)
          TQB=-TQ
        ELSE
          TQ=1.
          TQB=1.
        ENDIF
        DO 2200 I=1,LST(12)
        T1=-K(3,2)*QC(I)
        IF(T1.GT.0) THEN
          WQ=XPQ(I)*SIG*TQ
          WQB=0.
        ELSE
          WQB=XPQ(-I)*SIG*TQB
          WQ=0.
        ENDIF
C...Include y-dependence.
        WQ=WQ*PARI(23+IP)
        WQB=WQB*PARI(23+IP)
        PQH(I,IH)=PQH(I,IH)+WQ
        PQH(I+LST(12),IH)=PQH(I+LST(12),IH)+WQB
        PQH(17,IH)=PQH(17,IH)+WQ+WQB
 2200   CONTINUE
      ENDIF
 2300 CONTINUE
 2500 CONTINUE
      DO 2600 I=1,17
 2600 PQ(I)=(1.-PARL(6))/2.*PQH(I,1)+(1.+PARL(6))/2.*PQH(I,2)
      IH=1
      IF(LST(30).EQ.1) IH=2
      IF(LST(32).EQ.0) THEN
C...Simulation: cross section for chosen helicity state only.
        DSIGMA=PQH(17,IH)
      ELSEIF(LST(19).EQ.-1) THEN
C...Integration event-by-event: normalize and include alpha_s*1/(1-xp)
        DSIGMA=PQH(17,IH)/PARI(20)*PARL(25)/(1.-XP)
C...Max of dsigma/dxp for L- and R-handed lepton.
        IF(PQH(17,1).GT.PARI(15)) PARI(15)=PQH(17,1)
        IF(PQH(17,2).GT.PARI(16)) PARI(16)=PQH(17,2)
      ELSE
C...Integration for grid: normalize and include alpha_s*1/(1-xp)
        DSIGMA=PQ(17)/PARI(20)*PARL(25)/(1.-XP)
        IF(LST(17).EQ.0) THEN
C...Fixed beam energy, max of dsigma/dxp for L- and R-handed lepton.
          IF(PQH(17,1).GT.PARI(15)) PARI(15)=PQH(17,1)
          IF(PQH(17,2).GT.PARI(16)) PARI(16)=PQH(17,2)
        ELSE
C...Variable beam energy, max of dsigma/dxp for S, T, I contributions.
          IF(PQ(17)/PARI(23+LST(32)).GT.PARI(14+LST(32)))
     &    PARI(14+LST(32))=PQ(17)/PARI(23+LST(32))
        ENDIF
      ENDIF
      RETURN

C...Boson-gluon fusion, i.e. qq-event.
 3000 S13=Q2*(1.-XP)/XP
      IF(S13.LT.4.*AMU**2) RETURN
      DO 3500 IP=IL,IU
      SIG=XPQ(0)*DQCDI(2,IP,XP,ZPMIN,ZPMAX)
      DO 3300 IH=1,2
      IF(IH.EQ.1) THEN
        IF(PARL(6).GT.0.99) GOTO 3300
        IF(LST(32).EQ.0.AND.LST(30).NE.-1) GOTO 3300
      ELSEIF(IH.EQ.2) THEN
        IF(PARL(6).LT.-0.99) GOTO 3300
        IF(LST(32).EQ.0.AND.LST(30).NE.1) GOTO 3300
      ENDIF
      IF(LST(32).NE.0) LST(30)=SIGN(1.,IH-1.5)
      IF(LST(23).NE.2) THEN
        DO 3100 I=1,LST(13)
        MSTJ(93)=1
        IF(S13.LT.4.*ULMASS(I)**2) GOTO 3100
        WQ=SIG/2.*(EWQC(1,IH,I)+EWQC(2,IH,I))
        WQB=WQ
C...Include y-dependence.
        WQ=WQ*PARI(23+IP)
        WQB=WQB*PARI(23+IP)
        PQH(I,IH)=PQH(I,IH)+WQ
        PQH(I+LST(13),IH)=PQH(I+LST(13),IH)+WQB
        PQH(17,IH)=PQH(17,IH)+WQ+WQB
 3100   CONTINUE
      ELSEIF(LST(23).EQ.2) THEN
C...Zero CC cross-section for one helicity state.
        IF(KSAVE(1).LT.0.AND.IH.EQ.1
     &  .OR.KSAVE(1).GT.0.AND.IH.EQ.2) GOTO 3300
        DO 3200 I=1,LST(13)
        MSTJ(93)=1
        IF(S13.LT.(AMU+ULMASS(I))**2) GOTO 3200
        IF(K(3,2)*QC(I).LT.0) THEN
          WQ=SIG
          WQB=0.
        ELSE
          WQB=SIG
          WQ=0.
        ENDIF
C...Include y-dependence.
        WQ=WQ*PARI(23+IP)
        WQB=WQB*PARI(23+IP)
        PQH(I,IH)=PQH(I,IH)+WQ
        PQH(I+LST(13),IH)=PQH(I+LST(13),IH)+WQB
        PQH(17,IH)=PQH(17,IH)+WQ+WQB
 3200   CONTINUE
      ENDIF
 3300 CONTINUE
 3500 CONTINUE
      DO 3600 I=1,17
 3600 PQ(I)=(1.-PARL(6))/2.*PQH(I,1)+(1.+PARL(6))/2.*PQH(I,2)
      IH=1
      IF(LST(30).EQ.1) IH=2
      IF(LST(32).EQ.0) THEN
C...Simulation: cross section for chosen helicity state only.
        DSIGMA=PQH(17,IH)
      ELSEIF(LST(19).EQ.-1) THEN
C...Integration event-by-event: normalize and include alpha_s*1/(1-xp)
        DSIGMA=PQH(17,IH)/PARI(20)*PARL(25)/(1.-XP)
C...Max of dsigma/dxp for L- and R-handed lepton.
        IF(PQH(17,1).GT.PARI(18)) PARI(18)=PQH(17,1)
        IF(PQH(17,2).GT.PARI(19)) PARI(19)=PQH(17,2)
      ELSE
C...Integration for grid: normalize and include alpha_s*1/(1-xp)
        DSIGMA=PQ(17)/PARI(20)*PARL(25)/(1.-XP)
        IF(LST(17).EQ.0) THEN
C...Fixed beam energy, max of dsigma/dxp for L- and R-handed lepton.
          IF(PQH(17,1).GT.PARI(18)) PARI(18)=PQH(17,1)
          IF(PQH(17,2).GT.PARI(19)) PARI(19)=PQH(17,2)
        ELSE
C...Variable beam energy, max of dsigma/dxp for S, T, I contributions.
          IF(PQ(17)/PARI(23+LST(32)).GT.PARI(17+LST(32)))
     &    PARI(17+LST(32))=PQ(17)/PARI(23+LST(32))
        ENDIF
      ENDIF
      RETURN
      END

      FUNCTION DSIGM2(X)
      IMPLICIT NONE
      REAL DSIGM2,X

CAE...Modified version of DSIGMA with a factor 1/(1-XP) is removed
C...by a variable transformation. This makes the integration of the 
C...cross section faster.

      REAL XP,DSIGMA
      
      XP=1.0-EXP(X)
      DSIGM2=(1.0-XP)*DSIGMA(XP)

      RETURN
      END

C **********************************************************************

      FUNCTION DQCD(ICOSFI,IPART,IP,XP,ZP,Y)

C...First order QCD matrix elements from R.D. Peccei and R. Ruckl:
C...Nucl. Phys. B162 (1980) 125

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,YY,W2,Q2,U
C...Constants C1 to C5 are resp. 2/3/pi, 1/4/pi, 4/3/pi, 1/2/pi, 1/pi
      DATA C1,C2,C3,C4,C5/0.2122066,0.0795775,0.4244132,0.1591549,
     &     0.3183099/

      IF(LST(8).EQ.19.AND.IPART.EQ.1) THEN
C...No QCD Compton for Ariadne
        DQCD=0.0
        RETURN
      ENDIF

      IF(ICOSFI.EQ.0) THEN
        IF(IPART.EQ.1) THEN
          IF(IP.EQ.1) THEN
            DQCD=C1*((ZP**2+XP**2)/(1.-XP)/(1.-ZP)+2.*(XP*ZP+1.))
          ELSEIF(IP.EQ.2) THEN
            DQCD=C1*4.*XP*ZP
          ELSEIF(IP.EQ.3) THEN
            DQCD=C1*((ZP**2+XP**2)/(1.-XP)/(1.-ZP)+2.*(XP+ZP))
          ELSE
            WRITE(6,1000) ICOSFI,IPART,IP
          ENDIF
        ELSEIF(IPART.EQ.2) THEN
          IF(IP.EQ.1) THEN
            DQCD=C2*(XP**2+(1.-XP)**2)*(ZP**2+(1.-ZP)**2)/(1.-ZP)/ZP
          ELSEIF(IP.EQ.2) THEN
            DQCD=C2*8.*XP*(1.-XP)
          ELSEIF(IP.EQ.3) THEN
            DQCD=C2*(XP**2+(1.-XP)**2)*(ZP-(1.-ZP))/(1.-ZP)/ZP
          ELSE
            WRITE(6,1000) ICOSFI,IPART,IP
          ENDIF
        ELSE
          WRITE(6,1000) ICOSFI,IPART,IP
        ENDIF

      ELSEIF(ICOSFI.EQ.1) THEN
        IF(IPART.EQ.1) THEN
          IF(IP.EQ.1) THEN
            DQCD=C3*Y*SQRT((1.-Y)*XP*ZP/(1.-XP)/(1.-ZP))*
     &      (1.-2./Y)*(1.-ZP-XP+2.*XP*ZP)
          ELSEIF(IP.EQ.3) THEN
            DQCD=C3*Y*SQRT((1.-Y)*XP*ZP/(1.-XP)/(1.-ZP))*
     &      (1.-XP-ZP)
          ELSE
            WRITE(6,1000) ICOSFI,IPART,IP
          ENDIF
        ELSEIF(IPART.EQ.2) THEN
          IF(IP.EQ.1) THEN
            DQCD=C4*Y*SQRT((1.-Y)*XP*(1.-XP)/ZP/(1.-ZP))*
     &      (1.-2./Y)*(1.-2.*ZP)*(1.-2.*XP)
          ELSEIF(IP.EQ.3) THEN
            DQCD=C4*Y*SQRT((1.-Y)*XP*(1.-XP)/ZP/(1.-ZP))*
     &      (1.-2.*XP)
          ELSE
            WRITE(6,1000) ICOSFI,IPART,IP
          ENDIF
        ENDIF

      ELSEIF(ICOSFI.EQ.2) THEN
        IF(IPART.EQ.1) THEN
          DQCD=C3*(1.-Y)*XP*ZP
        ELSEIF(IPART.EQ.2) THEN
          DQCD=C5*(1.-Y)*XP*(1.-XP)
        ELSE
          WRITE(6,1000) ICOSFI,IPART,IP
        ENDIF

      ELSE
        WRITE(6,1000) ICOSFI,IPART,IP
      ENDIF
      RETURN

 1000 FORMAT(' Error in routine DQCD     ',
     &' ICOSFI, IPART, IP = ',3I10)
      END

C **********************************************************************

      FUNCTION DQCDI(IPART,IP,XP,ZPMIN,ZPMAX)

C...First order QCD matrix elements as in function DQCD but analytically
C...integrated over ZP from ZPMIN to ZPMAX, also a factor 1/(1-XP) is
C...factored out (since XP chosen randomly according to 1/(1-XP) distr.)

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      DATA C1,C2/0.2122066,0.0795775/

      IF(LST(8).EQ.19.AND.IPART.EQ.1) THEN
C...No QCD Compton for Ariadne
        DQCDI=0.0
        RETURN
      ENDIF

      IF(IPART.EQ.1) THEN
        IF(IP.EQ.1) THEN
          ZLOG=ALOG(ZPMAX/ZPMIN)
          DQCDI=C1*(XP**2*ZLOG+ZPMIN-ZPMAX+(ZPMIN**2-ZPMAX**2)/2.+ZLOG+
     &    XP*(1.-XP)*(ZPMAX**2-ZPMIN**2)+2.*(1.-XP)*(ZPMAX-ZPMIN))
        ELSEIF(IP.EQ.2) THEN
          DQCDI=C1*2.*XP*(1.-XP)*(ZPMAX**2-ZPMIN**2)
        ELSEIF(IP.EQ.3) THEN
          ZLOG=ALOG(ZPMAX/ZPMIN)
          DQCDI=C1*(XP**2*ZLOG+ZPMIN-ZPMAX+(ZPMIN**2-ZPMAX**2)/2.+ZLOG+
     &    2.*XP*(1.-XP)*(ZPMAX-ZPMIN)+(1.-XP)*(ZPMAX**2-ZPMIN**2))
        ELSE
          WRITE(6,1000) IPART,IP
        ENDIF

      ELSEIF(IPART.EQ.2) THEN
        IF(IP.EQ.1) THEN
          DQCDI=C2*(1.-XP)*(XP**2+(1.-XP)**2)*(2.*(ZPMIN-ZPMAX)+
     &    ALOG(ZPMAX**2/ZPMIN**2))
        ELSEIF(IP.EQ.2) THEN
          DQCDI=C2*8.*XP*(1.-XP)**2*(ZPMAX-ZPMIN)
        ELSEIF(IP.EQ.3) THEN
          DQCDI=0.
        ELSE
          WRITE(6,1000) IPART,IP
        ENDIF

      ELSE
        WRITE(6,1000) IPART,IP
      ENDIF
      RETURN

 1000 FORMAT(' Error in routine DQCDI     ',
     &' IPART, IP = ',2I10)
      END

C **********************************************************************

      SUBROUTINE LFLAV(IFL,IFLR)

C...Choose flavour of struck quark and the
C...corresponding flavour of the target remnant jet.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LFLMIX/ CABIBO(4,4)

      LST(21)=0
      IF(LST(24).EQ.3) THEN
        NFL=LST(13)
      ELSE
        NFL=LST(12)
      ENDIF

   20 R=RLU(0)*PQ(17)
      PSUB=0.
      DO 30 I=1,2*NFL
      IFL=I
      PSUB=PSUB+PQ(I)
      IF(R.LE.PSUB) GOTO 40
   30 CONTINUE
   40 CONTINUE
      IF(IFL.GT.NFL) IFL=NFL-IFL
      LST(25)=IFL
      IFLR=-IFL

      IF(LST(23).EQ.2) THEN
C...Weak charged current, change the flavour of the struck
C...quark using generalized Cabibbo mixing matrix.
        IFLA=IABS(IFL)
        J1=(IFLA+1)/2
        M1=MOD(IFLA,2)
        M2=MOD(IFLA+1,2)
        R=RLU(0)
        PSUB=0.
        DO 100 J=1,4
        J2=J
        PSUB=PSUB+CABIBO(M1*J2+M2*J1,M2*J2+M1*J1)
        IF(R.LT.PSUB) GOTO 200
  100   CONTINUE
  200   IFL=2*J2-M2
        IF(LST(25).LT.0) IFL=-IFL
      ENDIF

      IFLA=IABS(IFL)
      IFLRA=IABS(IFLR)
      IF(IFLA.GE.4.OR.IFLRA.GE.4) THEN
C...Threshold function for heavy quarks of flavour IFLA and IFLRA.
        MSTJ(93)=1
        AMU=ULMASS(1)
        MSTJ(93)=1
        AMIFL=ULMASS(IFLA)
        MSTJ(93)=1
        AMIFLR=ULMASS(IFLRA)
        IF(1.-(.938+AMIFL+AMIFLR+2.*AMU)**2/W2.LT.RLU(0))
     &  GOTO(20,999,999) LST(24)
      ENDIF

C...Remnant flavour taken care of later for qqbar event and ME+PS case
      IF(LST(24).EQ.3) RETURN
      IF(LST(8).GT.10.AND.MOD(LST(8),10).NE.9) RETURN

C...With LST(14)=0/1(default) baryon production from the target remnant
C...is excluded/included.
      IF(LST(14).EQ.0) RETURN
      IF(IFLR.EQ.-2) THEN
        IF(LST(22).EQ.1) THEN
          IFLR=2101
          IF(RLU(0).GT.PARL(4)) IFLR=2103
        ELSE
          IFLR=1103
        ENDIF
      ELSEIF(IFLR.EQ.-1) THEN
        IF(LST(22).EQ.1) THEN
          IFLR=2203
        ELSE
          IFLR=2101
          IF(RLU(0).GT.PARL(4)) IFLR=2103
        ENDIF
      ENDIF
      RETURN

  999 LST(21)=6
      RETURN
      END

C **********************************************************************

      SUBROUTINE LREMH(IFLRO,PT,IFLR,K2,Z)

C...Gives flavour and energy-momentum fraction Z for the particle
C...to be produced out of the target remnant when that is not a
C...simple diquark.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)

C...Flavours fixed when calling from LQQBEV OR LYREMN
      IF(IFLRO.EQ.0) GOTO 200

C...Split target remnant qqqQ -> qqQ + q or qqqQbar -> qQbar + qq
C...Q (Qbar) is the partner to the struck sea quark
C...qqq are the nucleon valence quarks from which a quark q or a
C...diquark qq is chosen at random to form a jet system with the
C...scattered sea antiquark or quark, respectively, the other parton
C...forms a baryon qqQ or meson qQbar, respectively.
  100 IFLQ=INT(1.+LST(22)/3.+RLU(0))
      IF(IFLQ.EQ.LST(22)) THEN
        IFLQQ=2101
        IF(RLU(0).GT.PARL(4)) IFLQQ =2103
      ELSE
        IFLQQ=1000*IFLQ+100*IFLQ+3
      ENDIF
      IFLQ=3-IFLQ

C...Choose flavour of hadron and parton for jet system
      IF(IFLRO.GT.0) THEN
        CALL LUKFDI(IFLQQ,IFLRO,IDUM,K2)
        IF(K2.EQ.0) GOTO 100
        IFLR=IFLQ
      ELSE
        CALL LUKFDI(IFLQ,IFLRO,IDUM,K2)
        IF(K2.EQ.0) GOTO 100
        IFLR=IFLQQ
      ENDIF

C...Entry from LQQBEV & LYREMN with flavours given, choose E-p fraction
  200 KSP=IFLR
C...Split energy-momentum of target remnant according to functions P(z)
C...z=E-pz fraction for qq (q) forming jet-system with struck Q (Qbar)
C...1-z=E-pz fraction for qQbar (qqQ) hadron
C...mq=mass of (light) parton remnant q (qq) in jet system
C...mQ=mass of produced (heavy flavour) hadron
      MSTJ(93)=1
      AMSP=ULMASS(KSP)
      MSTJ(93)=1
      AMK2=ULMASS(K2)
      IF(LST(14).EQ.1) THEN
C...P(z)=(a+1)(1-z)**a with <z>=1/(a+2)=1/3 since a=1 fixed
        Z=1.-SQRT(RLU(0))
C...Flip if baryon produced
        KC2=IABS(LUCOMP(K2))
        IF(KC2.GE.301.AND.KC2.LE.400) Z=1.-Z
      ELSEIF(LST(14).EQ.2) THEN
C...P(z)=(a+1)(1-z)**a with <z>=1/(a+2)=mq/(mq+mQ) --> a=a(mq,mQ)
        A=(AMSP+AMK2)/AMSP - 2.
        Z=RLU(0)**(1./(A+1.))
      ELSEIF(LST(14).EQ.3) THEN
C...Using Peterson fragmentation function
C...P(z)=N/(z(1-1/z-c/(1-z))**2)  where c=(mq/mQ)**2  (FC=-c)
        FC=-(AMSP/AMK2)**2
  300   Z=RLU(0)
        IF(-4.*FC*Z*(1.-Z)**2.LT.RLU(0)*((1.-Z)**2-FC*Z)**2) GOTO 300
      ELSEIF(LST(14).EQ.4) THEN
C...Using chosen fragmentation function in JETSET
        TM2=AMK2**2+PT**2
        CALL LUZDIS(1,0,TM2,Z)
      ENDIF
      LST(27)=1
      K2A=IABS(K2)
      IF((K2A.GE.1.AND.K2A.LE.8).OR.K2A.EQ.21.OR.LUCOMP(K2A).EQ.90)
     &LST(27)=2

      RETURN
      END

C **********************************************************************

      SUBROUTINE LPRIKT(S,PT,PHI)

C...Size (PT) and azimuthal angle (PHI) of primordial kt according
C...to a Gaussian distribution.

      PT=S*SQRT(-ALOG(RLU(0)))
      PHI=6.2832*RLU(0)
      RETURN
      END

C **********************************************************************

      SUBROUTINE LFRAME(IFR,IPH)

C...Make transformation from hadronic CM frame to lab frame.

      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
      DOUBLE PRECISION DTHETA,DPHI,DBETA
      INTEGER IFR,IPH,IFRAME,IPHI

      IFRAME=IFR
      IPHI=IPH
      IF(IFRAME.LT.1.OR.IFRAME.GT.4.OR.IPHI.LT.0.OR.IPHI.GT.1) GOTO 999
      IF(IFRAME.EQ.1) IPHI=0

CGI -- Also boost lines up to N+MSTU(3)
      N=N+1+MSTU(3)
      DO 5 J=1,5
    5 P(N,J)=PB(J)

   10 CONTINUE
      IF(IPHI.NE.LST(29)) THEN
         IFRAME=2
      ELSE
         IFRAME=IFR
      ENDIF
      IF((IFRAME.EQ.LST(28)).AND.(IPHI.EQ.LST(29))) THEN
        DO 15 J=1,5
   15   PB(J)=P(N,J)
        N=N-1-MSTU(3)
        RETURN
      ENDIF

      GOTO(100,200,300,400), LST(28)
      GOTO 999

  100 IF(IFRAME.GE.2) THEN
        CALL LUDBRB(0,0,STHETA(2),SPHI(2),0.D0,0.D0,0.D0)
        CALL LUDBRB(0,0,0.,0.,DBETA(2,1),DBETA(2,2),DBETA(2,3))
        LST(28)=2
      ELSE
        GOTO 999
      ENDIF
      GOTO 10

  200 IF(IPHI.NE.LST(29)) THEN
       CALL LUDBRB(0,0,0.,SIGN(PHIR,FLOAT(IPHI-LST(29))),0.D0,0.D0,0.D0)
       LST(29)=IPHI
      ENDIF

      IF(IFRAME.EQ.1) THEN
        CALL LUDBRB(0,0,0.,0.,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
        CALL LUDBRB(0,0,-STHETA(2),0.,0.D0,0.D0,0.D0)
        LST(28)=1
      ELSEIF(IFRAME.GE.3) THEN
        IF(LST(17).EQ.0) THEN
          CALL LUDBRB(0,0,0.,0.,0.D0,0.D0,DBETA(1,3))
          IF(PSAVE(3,1,3).LT.0.) THEN
            DO 210 I=1,N
            V(I,3)=-V(I,3)
  210       P(I,3)=-P(I,3)
          ENDIF
        ELSE
          CALL LUDBRB(0,0,STHETA(1),SPHI(1),0.D0,0.D0,0.D0)
          CALL LUDBRB(0,0,0.,0.,DBETA(1,1),DBETA(1,2),DBETA(1,3))
        ENDIF
        LST(28)=3
      ENDIF
      GOTO 10

  300 IF(IFRAME.LE.2) THEN
        IF(LST(17).EQ.0) THEN
          IF(PSAVE(3,1,3).LT.0.) THEN
            DO 310 I=1,N
            V(I,3)=-V(I,3)
  310       P(I,3)=-P(I,3)
          ENDIF
          CALL LUDBRB(0,0,0.,0.,0.D0,0.D0,-DBETA(1,3))
        ELSE
          CALL LUDBRB(0,0,0.,0.,-DBETA(1,1),-DBETA(1,2),-DBETA(1,3))
          CALL LUDBRB(0,0,0.,-SPHI(1),0.D0,0.D0,0.D0)
          CALL LUDBRB(0,0,-STHETA(1),0.,0.D0,0.D0,0.D0)
        ENDIF
        LST(28)=2
      ELSEIF(IFRAME.EQ.4) THEN
        THEBOS=PLU(N,13)
        PHIBOS=PLU(N,15)
        CALL LUDBRB(0,0,0.,-PHIBOS,0.D0,0.D0,0.D0)
        CALL LUDBRB(0,0,-THEBOS,0.,0.D0,0.D0,0.D0)
        LST(28)=4
      ENDIF
      GOTO 10

  400 IF(IFRAME.LE.3) THEN
        CALL LUDBRB(0,0,THEBOS,PHIBOS,0.D0,0.D0,0.D0)
        LST(28)=3
      ENDIF
      GOTO 10

  999 WRITE(6,1000) IFRAME,IPHI,LST(28),LST(29)
 1000 FORMAT(' BAD VARIABLES IN SUBROUTINE LFRAME: IFRAME,IPHI,',
     &'LST(28),LST(29) =',4I5)
      RETURN
      END

C ********************************************************************

      SUBROUTINE LWBB(ENU)

C...Gives energy (ENU) of a (anti-)neutrino chosen from a simple
C...parametrization of a wide band beam.

      DATA EMEAN,SLOPE,EMIN,EMAX/30.,0.02,12.,300./
      A1=1./(EMEAN-12.)
      A2=EXP(EMEAN*SLOPE)
  100 ENU=300.*RLU(0)
      IF(ENU.LT.EMEAN)THEN
      E=A1*(ENU-12.)
      ELSE
      E=A2*EXP(-ENU*SLOPE)
      ENDIF
      IF(ENU.LT.EMIN.OR.ENU.GT.EMAX) GOTO 100
      IF(E.LT.RLU(0)) GOTO 100
      RETURN
      END

C **********************************************************************

      SUBROUTINE LWEITS(LFILE)

C...Integrates the QCD matrix elements to obtain probabilities for
C...qg- and qq-events as a function of (x,W). Also finds various
C...maximum values to be used for the QCD simulation. Results stored
C...in common LGRID and optionally written to logical file LFILE.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON /LGRID/ NXX,NWW,XX(31),WW(21),PQG(31,21,3),PQQB(31,21,2),
     &QGMAX(31,21,3),QQBMAX(31,21,2),YCUT(31,21),XTOT(31,21),NP
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      DIMENSION WWI(21,4),XXI(31,4)
      EXTERNAL DSIGMA,DSIGM2
      LOGICAL ZOOM
      DATA WWI/
     1  5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,
     1  20.,21.,22.,23.,24.,25.,
     2  5.,7.,9.,11.,13.,15.,17.,19.,21.,23.,25.,27.,29.,31.,33.,
     2  35.,37.,39.,41.,43.,45.,
     3  5.,7.,10.,15.,20.,25.,30.,40.,50.,75.,100.,125.,150.,175.,
     3  200.,225.,250.,275.,300.,325.,350.,
     4  5.,10.,15.,20.,30.,50.,75.,100.,150.,200.,250.,300.,400.,
     4  500.,700.,900.,1200.,1500.,1800.,2100.,2500./
      DATA XXI/
     1  1.E-3,2.E-3,3.E-3,4.E-3,5.E-3,6.E-3,8.E-3,
     1  1.E-2,2.E-2,3.E-2,4.E-2,5.E-2,6.E-2,8.E-2,
     1  .1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.72,.8,.87,.94,.999,
     2  1.E-3,2.E-3,3.E-3,4.E-3,5.E-3,6.E-3,8.E-3,
     2  1.E-2,2.E-2,3.E-2,4.E-2,5.E-2,6.E-2,8.E-2,
     2  .1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.72,.8,.87,.94,.999,
     3  1.E-5,2.E-5,4.E-5,6.E-5,8.E-5,1.E-4,2.E-4,4.E-4,6.E-4,8.E-4,
     3  1.E-3,2.E-3,4.E-3,6.E-3,8.E-3,1.E-2,2.E-2,4.E-2,6.E-2,8.E-2,
     3  .1,.2,.3,.4,.5,.6,.7,.8,.87,.94,.999,
     4  1.E-5,2.E-5,4.E-5,6.E-5,8.E-5,1.E-4,2.E-4,4.E-4,6.E-4,8.E-4,
     4  1.E-3,2.E-3,4.E-3,6.E-3,8.E-3,1.E-2,2.E-2,4.E-2,6.E-2,8.E-2,
     4  .1,.2,.3,.4,.5,.6,.7,.8,.87,.94,.999/
      DATA NCALL/0/,XSPLIT/0.1/,YSPLIT/2./

      NCALL=NCALL+1
      LST2=LST(2)
      LST(2)=-3
      WMAX=SQRT(PARL(21)) 

      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &WRITE(6,1000) PARL(11),LST(13),MSTU(112),PARU(112),
     &PARL(8),PARL(9),PARL(12),PARL(13)
      IF(LST(17).EQ.0) THEN
        NP=1
        IPMAX=2
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,1010)
      ELSE
        NP=3
        IF(LST(23).EQ.1) NP=2
        IPMAX=3
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) THEN
          IF(LST(19).GE.0.AND.LST(19).LT.10) WRITE(6,1020)
          IF(LST(19).GE.10.OR.LST(19).EQ.-10) WRITE(6,2020)
        ENDIF
      ENDIF

      IF(LST(19).GE.0.AND.LST(19).LT.10) THEN
C...Fixed grid in x and W
        IF(LST(19).EQ.0)THEN
C...Grid specified by user.
          READ(5,*) NWW,NXX
          READ(5,*) (WW(IW),IW=1,NWW)
          READ(5,*) (XX(IX),IX=1,NXX)
          IF(XX(NXX).GT..99) XX(NXX)=.99
        ELSEIF(LST(19).GE.1.AND.LST(19).LE.4) THEN
C...Grid taken from data in arrays WWI, XXI.
          DO 10 IW=1,NWW
   10     WW(IW)=WWI(IW,LST(19))
          DO 20 IX=1,NXX
   20     XX(IX)=XXI(IX,LST(19))
        ENDIF
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &  WRITE(6,1030) LST(19),NXX,NWW,XX,WW
        IF(WMAX.GT.WW(NWW)) THEN
          IF(LST(3).GE.1) WRITE(6,1040) WMAX,WW(NWW)
          IF(LST(3).GE.2) THEN
            WRITE(6,1900)
            STOP
          ENDIF
        ENDIF
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,1100)
      ELSEIF(LST(19).GE.10.OR.LST(19).EQ.-10) THEN
C...Grid in x,y automatically defined from available x,y-ranges
        NX=NXX
        NY=NWW
        IF(XMIN.GE.XSPLIT) THEN
          DO 30 I=1,NX
   30     XX(I)=MIN(0.999,XMIN+(XMAX-XMIN)*(I-1)/FLOAT(NX-1))
        ELSEIF(XMAX.GT.XSPLIT) THEN
          NL=MIN(2.*NX/3.,MAX(NX/3.,NX*LOG(XSPLIT/XMIN)/LOG(XMAX/XMIN)))
          DO 40 I=1,NL
   40     XX(I)=MIN(0.999,
     &   10.**(LOG10(XMIN)+(LOG10(XSPLIT)-LOG10(XMIN))*(I-1)/FLOAT(NL)))
          DO 41 I=NL+1,NX
   41     XX(I)=MIN(0.999,XSPLIT+(XMAX-XSPLIT)*(I-NL-1)/FLOAT(NX-NL-1))
        ELSE
          DO 50 I=1,NX
   50     XX(I)=MIN(0.999,
     &   10.**(LOG10(XMIN)+(LOG10(XMAX)-LOG10(XMIN))*(I-1)/FLOAT(NX-1)))
        ENDIF
C...y-variable stored in same array as W
        IF(YMIN.GE.YSPLIT) THEN
          DO 60 I=1,NY
   60     WW(I)=MIN(0.999,YMIN+(YMAX-YMIN)*(I-1)/FLOAT(NY-1))
        ELSEIF(YMAX.GT.YSPLIT) THEN
          NL=MIN(2.*NY/3.,MAX(NY/3.,NY*LOG(YSPLIT/YMIN)/LOG(YMAX/YMIN)))
          DO 70 I=1,NL
   70     WW(I)=MIN(0.999,
     &   10.**(LOG10(YMIN)+(LOG10(YSPLIT)-LOG10(YMIN))*(I-1)/FLOAT(NL)))
          DO 71 I=NL+1,NY
   71     WW(I)=MIN(0.999,YSPLIT+(YMAX-YSPLIT)*(I-NL-1)/FLOAT(NY-NL-1))
        ELSE
          DO 80 I=1,NY
   80     WW(I)=MIN(0.999,
     &   10.**(LOG10(YMIN)+(LOG10(YMAX)-LOG10(YMIN))*(I-1)/FLOAT(NY-1)))
        ENDIF
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &  WRITE(6,2030) LST(19),NXX,NWW,XX,WW
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,2100)
      ENDIF

      LW=0
      DO 500 IW=1,NWW
      W=WW(IW)
      Y=WW(IW)
      IF(LW.GT.0) GOTO 600
      IF(LST(19).LT.10.AND.W.GT.WMAX) LW=LW+1
      LX=0
      DO 400 IX=1,NXX
      X=XX(IX)
      IF(LST(19).GE.0.AND.LST(19).LT.10) THEN
C...x,W given.
        W2=W**2
        U=(W2-PSAVE(3,2,5)**2)/(2.*PSAVE(3,2,5)*(1.-X))
        Q2=2.*PSAVE(3,2,5)*U*X
        Y=Q2/(PARL(21)*X)
      ELSEIF(LST(19).GE.10.OR.LST(19).EQ.-10) THEN
C...x,y given.
        PARL(22)=Y*PARL(21)
        Q2=X*PARL(22)
        U=PARL(22)/(2.*PSAVE(3,2,5))
        W2=PARL(22)*(1.-X)+PSAVE(3,2,5)**2
        W=SQRT(W2)
      ENDIF
C...Protection against too small Q2 in structure functions
Ctest IF(Q2.LT.2.098) GOTO 400
      IF(LX.GT.0) GOTO 500
      IF(LST(19).GE.0.AND.LST(19).LT.10.AND.X.GT.1.-W2/WMAX**2) LX=LX+1
      CALL LEPTO
      PQCOM=PARI(31)*PQ(17)*COMFAC
      PARL(25)=ULALPS(Q2)
      PARI(20)=PQ(17)
      XTOT(IX,IW)=PQ(17)
      IF(LST(20).LE.1) THEN
        PARL(27)=MAX(PARL(9)**2/W2,PARL(8))
        P27MAX=1.0
      ELSEIF(LST(20).EQ.2) THEN
        PARL(27)=MAX(PARL(9)**2/Q2,PARL(8))
        P27MAX=W2/Q2
      ELSEIF(LST(20).EQ.3.OR.LST(20).EQ.4) THEN
        PARL(27)=PARL(8)
        P27MAX=0.5
      ELSEIF(LST(20).EQ.5) THEN
        PARL(27)=PARL(8)
        P27MAX=0.5
      ELSEIF(LST(20).EQ.6) THEN
        PARL(27)=PARL(9)
        P27MAX=W2
      ENDIF

      ZOOM=.FALSE.
      YCMIN=PARL(27)
      YCMAX=PARL(27)
      IYCUT=0
  100 IYCUT=IYCUT+1
      RQG=0.
      RQQB=0.
CAE.Scheme for ME cutoff: W2, Q2, mixed
      IF(LST(20).LE.1) THEN
        XPMIN=DBLE(X)/(1.D0-2.D0*(1.D0-DBLE(X))*DBLE(PARL(27)))
        XPMAX=DBLE(X)/(DBLE(X)+(1.D0-DBLE(X))*DBLE(PARL(27)))
      ELSEIF(LST(20).EQ.2) THEN
        XPMIN=DBLE(X)/(1.D0-2.D0*DBLE(X)*DBLE(PARL(27)))
        XPMAX=1.D0/(1.D0+DBLE(PARL(27)))
      ELSEIF(LST(20).EQ.3.OR.LST(20).EQ.4) THEN
        XPMIN=X
        XPMAX=1.D0/(1.D0+DBLE(PARL(9)))
      ELSEIF(LST(20).EQ.5) THEN
        XPMIN=X
        XPMAX=Q2/(Q2+PARL(9))
      ELSEIF(LST(20).EQ.6) THEN
        XPMIN=X
        XPMAX=Q2/(Q2+PARL(27))
      ELSE
        WRITE(6,*) 'LWEITS: No such jet scheme!'
      ENDIF
CAE
      IF(XPMIN.GE.XPMAX.OR.XPMIN.LE.0.) GOTO 210
      DO 200 IP=1,NP
      IF(LST(17).EQ.0) THEN
        PARI(15)=0.
        PARI(16)=0.
        PARI(18)=0.
        PARI(19)=0.
      ELSE
        PARI(14+IP)=0.
        IF(IP.LE.2) PARI(17+IP)=0.
      ENDIF
      LST(32)=IP
      LST(24)=2
      EPS=PARL(11)
CAE      CALL GADAP(XPMIN,XPMAX,DSIGMA,EPS,RESULT)
      CALL GADAP(LOG(1.0-XPMAX),LOG(1.0-XPMIN),DSIGM2,EPS,RESULT)
      RQG=RQG+RESULT
      PQG(IX,IW,IP)=RESULT/PARL(25)
      IF(LST(17).EQ.0) THEN
        QGMAX(IX,IW,1)=PARI(15)
        QGMAX(IX,IW,2)=PARI(16)
      ELSE
        PQG(IX,IW,IP)=RESULT*PARI(20)/PARI(23+IP)/PARL(25)
        QGMAX(IX,IW,IP)=PARI(14+IP)
      ENDIF
      IF(IP.EQ.3) GOTO 200
      LST(24)=3
      EPS=PARL(11)
CAE      CALL GADAP(XPMIN,XPMAX,DSIGMA,EPS,RESULT)
      CALL GADAP(LOG(1.0-XPMAX),LOG(1.0-XPMIN),DSIGM2,EPS,RESULT)
      RQQB=RQQB+RESULT
      PQQB(IX,IW,IP)=RESULT/PARL(25)
      IF(LST(17).EQ.0) THEN
        QQBMAX(IX,IW,1)=PARI(18)
        QQBMAX(IX,IW,2)=PARI(19)
      ELSE
        PQQB(IX,IW,IP)=RESULT*PARI(20)/PARI(23+IP)/PARL(25)
        QQBMAX(IX,IW,IP)=PARI(17+IP)
      ENDIF
  200 CONTINUE
  210 CONTINUE
      RQ=1.-RQG-RQQB
      IF(.NOT.ZOOM) THEN
CAE.First find interval so that RQ>0
        IF(RQ.LT.0.AND.IYCUT.LT.10) THEN
          PARL(27)=MIN(1.1*EXP(-2.0*RQ)*PARL(27),P27MAX)
          YCMIN=YCMAX
          YCMAX=PARL(27)
        ELSEIF(RQ.LT.0.AND.IYCUT.GE.10) THEN
C...Terminate procedure after some iterations
          WRITE(6,*) 'Warning! sigma>tot for x,q2,cut=',X,Q2,PARL(27)
          WRITE(6,*) 'Weights=',RQ,RQG,RQQB
          RTOT=(RQG+RQQB)*1.05
          RQG=RQG/RTOT
          RQQB=RQQB/RTOT
          RQ=1.-RQG-RQQB
          DO 220 IP=1,3
            PQG(IX,IW,IP)=PQG(IX,IW,IP)/RTOT
            QGMAX(IX,IW,IP)=QGMAX(IX,IW,IP)/RTOT
220       CONTINUE
          DO 230 IP=1,2
            PQQB(IX,IW,IP)=PQQB(IX,IW,IP)/RTOT
            QQBMAX(IX,IW,IP)=QQBMAX(IX,IW,IP)/RTOT
230       CONTINUE
C...Break loop
          GOTO 250
        ELSEIF(IYCUT.GE.2.AND.RQ.GT.PARL(13)) THEN
C...If RQ>PARL(13), then ycut was increased to much
          ZOOM=.TRUE.
          PARL(27)=(YCMIN+YCMAX)/2.
        ELSE
C...correct ycut found
          GOTO 250
        ENDIF
      ELSE
C...Zoom in on ycut so that 0<RQ<PARL(13)
        IF(RQ.LT.0.AND.IYCUT.LT.100) THEN
          YCMIN=PARL(27)
          PARL(27)=(YCMIN+YCMAX)/2.
        ELSEIF(RQ.GT.PARL(13).AND.IYCUT.LT.100) THEN
          YCMAX=PARL(27)
          PARL(27)=(YCMIN+YCMAX)/2.
        ELSEIF(IYCUT.GE.100) THEN
          IF(LST(3).GE.1) THEN
            WRITE(6,*) 'LWEITS: Warning, PARL(27) not found.'
          ENDIF
          RTOT=(RQG+RQQB)*1.05
          RQG=RQG/MAX(1.0,RTOT)
          RQQB=RQQB/MAX(1.0,RTOT)
          RQ=1.-RQG-RQQB
C...Break loop
          GOTO 250        
        ELSE
C...ycut found, break loop
          GOTO 250
        ENDIF
      ENDIF
C...Loop until correct weights found
      GOTO 100

 250  CONTINUE
CAE
      YCUT(IX,IW)=PARL(27)
      IF(LST(33).EQ.-91) THEN
Ctest...Include 3-jet cross section in denominator
        QTOT=1.+RQG+RQQB
        RQG =RQG/QTOT
        RQQB=RQQB/QTOT
        RQ=1.-RQG-RQQB
      ENDIF
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) THEN
        IF(LST(19).LT.10) WRITE(6,1200) W,X,Y,Q2,PARL(25),PQCOM,
     &  PARL(27),IYCUT,RQ,RQG,RQQB,(QGMAX(IX,IW,IP),IP=1,IPMAX),
     &  (QQBMAX(IX,IW,IP),IP=1,MIN(2,IPMAX))
        IF(LST(19).GE.10) WRITE(6,2200) X,Y,Q2,W,PARL(25),PQCOM,
     &  PARL(27),IYCUT,RQ,RQG,RQQB,(QGMAX(IX,IW,IP),IP=1,IPMAX),
     &  (QQBMAX(IX,IW,IP),IP=1,MIN(2,IPMAX))
      ENDIF
  400 CONTINUE
  500 CONTINUE
  600 CONTINUE

      LST(2)=LST2
      IF(LFILE.LT.0) THEN
C...Write results on logical file number IABS(LFILE)
        WRITE(IABS(LFILE)) LST,PARL,NXX,NWW,NP,XX,WW
        WRITE(IABS(LFILE))(((PQG(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,NP),
     &  (((PQQB(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,NP),
     &  (((QGMAX(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,IPMAX),
     &  (((QQBMAX(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,MIN(2,IPMAX)),
     &  YCUT
        IF(NP.NE.1) WRITE(IABS(LFILE)) XTOT
        CLOSE(IABS(LFILE))
      ENDIF
      RETURN

 1000 FORMAT('1',/,5X,'Integration of 1st order QCD matrix elements',
     &           /,5X,'============================================',
     &/,' for gluon radiation (qg-event) and boson-gluon fusion ',
     &'(qq-event) probability.',
     &//,' Required precision in integration, PARL(11) =',F8.4,
     &//,' Heaviest flavour produced in boson-gluon fusion, LST(13) =',
     &I5,//,' Alpha-strong parameters: # flavours, MSTU(112) =',I3,
     &/,25X,' QCD lambda, PARU(112) =',F6.3,' GeV',
     &//,' Cuts on matrix elements:',
     &/,' PARL(8), PARL(9), PARL(12), PARL(13) =',4F8.4,/)
 1010 FORMAT(' Lepton energy not allowed to vary in simulation.',/)
 1020 FORMAT(' Lepton energy allowed to vary in simulation, ',/,
     &' y in table below calculated assuming max energy.',/)
 1030 FORMAT(' Grid choice, LST(19) =',I3,5X,'# grid points in x, W =',
     &2I5,/,' x-values in array XX:',/,10F8.5,/,10F8.5,/,11F8.5,
     &    /,' W-values in array WW:',/,10F7.1,/,11F7.1,/)
 1040 FORMAT(' Warning: max W outside grid, Wmax, grid-max =',2F12.1)
 1100 FORMAT(//,6X,'W',7X,'x',7X,'y',6X,'Q**2',1X,'alpha',1X,'dsigma',
     &9X,'cut',' it',2X,'q-event',1X,'qg-event',
     &1X,'qq-event',' max of matrix elements qg & qq; L,R or T,S,I',
     &/,1X,132(1H-),/)
 1200 FORMAT(F7.1,2F8.4,1PG10.3,0PF6.2,1PG11.3,0PF8.4,I3,3F9.4,1P,5E9.2)
 1900 FORMAT(' Execution stopped ',/)
 2020 FORMAT(' Lepton energy allowed to vary in simulation, ',/,
     &' W in table below calculated assuming max energy.',/)
 2030 FORMAT(' Grid choice, LST(19) =',I3,5X,'# grid points in x, y =',
     &2I5,/,' x-values in array XX:',/,10F8.5,/,10F8.5,/,11F8.5,
     &    /,' y-values in array WW:',/,10F7.4,/,11F7.4,/)
 2100 FORMAT(//,7X,'x',7X,'y',6X,'Q**2',6X,'W',1X,'alpha',1X,'dsigma',
     &9X,'cut',' it',2X,'q-event',1X,'qg-event',
     &1X,'qq-event',' max of matrix elements qg & qq; L,R or T,S,I',
     &/,1X,132(1H-),/)
 2200 FORMAT(2F8.5,1PG10.3,0PF7.1,F6.2,1PG11.3,0PF8.4,I3,3F9.4,1P,5E9.2)
      END

C **********************************************************************

      SUBROUTINE LPRWTS(NSTEP)

C...Prints probabilities for q-, qg- and qqbar-events using the present
C...QCD weights stored in common block LGRID.

      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LGRID/ NXX,NWW,XX(31),WW(21),PQG(31,21,3),PQQB(31,21,2),
     &QGMAX(31,21,3),QQBMAX(31,21,2),YCUT(31,21),XTOT(31,21),NP
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ

      WMAX=SQRT(PARL(21)+PSAVE(3,1,5)**2+PSAVE(3,2,5)**2)
      WRITE(6,1000) PARL(11),LST(13),MSTU(112),PARU(112),
     &PARL(8),PARL(9),PARL(12),PARL(13)
      IF(NP.EQ.1) THEN
        WRITE(6,1010)
      ELSE
        IF(LST(19).LT.10) WRITE(6,1020)
        IF(LST(19).GE.10) WRITE(6,2020)
      ENDIF
      IF(LST(19).LT.10) THEN
        WRITE(6,1030) LST(19),NXX,NWW,XX,WW
        IF(WMAX.GT.WW(NWW)) WRITE(6,1040) WMAX,WW(NWW)
        WRITE(6,1100)
      ELSEIF(LST(19).GE.10) THEN
        WRITE(6,2030) LST(19),NXX,NWW,XX,WW
        WRITE(6,2100)
      ENDIF

      LW=0
      DO 500 IW=1,NWW,MAX(1,NSTEP)
      W=WW(IW)
      Y=WW(IW)
      IF(LW.GT.0) GOTO 600
      IF(LST(19).LT.10.AND.W.GT.WMAX) LW=LW+1
      W2=W**2
      LX=0
      DO 400 IX=1,NXX,MAX(1,NSTEP)
      X=XX(IX)
      IF(LX.GT.0) GOTO 500
      IF(LST(19).LT.10) THEN
C...x,W2 given.
        U=(W2-PSAVE(3,2,5)**2)/(2.*PSAVE(3,2,5)*(1.-X))
        Q2=2.*PSAVE(3,2,5)*U*X
        Y=Q2/(PARL(21)*X)
      ELSEIF(LST(19).GE.10) THEN
C...x,y given.
        PARL(22)=Y*PARL(21)
        Q2=X*PARL(22)
        U=PARL(22)/(2.*PSAVE(3,2,5))
        W2=PARL(22)*(1.-X)+PSAVE(3,2,5)**2
      ENDIF
      PARI(24)=(1.+(1.-Y)**2)/2.
      PARI(25)=1.-Y
      PARI(26)=(1.-(1.-Y)**2)/2.
      PARL(25)=ULALPS(Q2)
      IF(Y.GT.1.) LX=LX+1
      RQG=0.
      RQQB=0.
      DO 200 IP=1,NP
      IF(NP.EQ.1) THEN
        RQG=PQG(IX,IW,IP)
        RQQB=PQQB(IX,IW,IP)
      ELSE
        RQG=RQG+PQG(IX,IW,IP)*PARI(23+IP)/XTOT(IX,IW)
        IF(IP.LT.3) RQQB=RQQB+PQQB(IX,IW,IP)*PARI(23+IP)/XTOT(IX,IW)
      ENDIF
  200 CONTINUE
C...Include alpha-strong in weight.
      RQG=RQG*PARL(25)
      RQQB=RQQB*PARL(25)
      IF(LST(33).EQ.-91) THEN
C...Include 3-jet cross section in denominator
        QTOT=1.+RQG+RQQB
        RQG =RQG/QTOT
        RQQB=RQQB/QTOT
      ENDIF
      RQ=1.-RQG-RQQB
      IF(LST(19).LT.10) THEN
        WRITE(6,1200) W,X,Y,Q2,PARL(25),YCUT(IX,IW),RQ,RQG,RQQB
      ELSEIF(LST(19).GE.10) THEN
        WRITE(6,2200) X,Y,Q2,W,PARL(25),YCUT(IX,IW),RQ,RQG,RQQB
      ENDIF
  400 CONTINUE
  500 CONTINUE
  600 CONTINUE
      RETURN

 1000 FORMAT('1',/,5X,'SUMMARY OF QCD MATRIX ELEMENT INTEGRATION',
     &           /,5X,'-----------------------------------------',//,
     &/,' for gluon radiation (qg-event) and boson-gluon fusion ',
     &'(qq-event) probability.',
     &//,' Required precision in integration, PARL(11) =',F8.4,
     &//,' Heaviest flavour produced in boson-gluon fusion, LST(13) =',
     &I5,//,' Alpha-strong parameters: # flavours, MSTU(112) =',I3,
     &'  QCD lambda, PARU(112) =',F6.3,' GeV',
     &//,' Cuts on matrix elements:',
     &/,' PARL(8), PARL(9), PARL(12), PARL(13) =',4F8.4,/)
 1010 FORMAT(' Lepton energy not allowed to vary in simulation.',/)
 1020 FORMAT(' Lepton energy allowed to vary in simulation, ',/,
     &' y in table below calculated assuming max energy.',/)
 1030 FORMAT(' Grid choice, LST(19) =',I3,5X,'# grid points in x, W =',
     &2I5,/,' x-values in array XX:',/,10F8.5,/,10F8.5,/,11F8.5,
     &    /,' W-values in array WW:',/,10F7.1,/,11F7.1,/)
 1040 FORMAT(' Max W outside grid, execution stopped ! Wmax, grid-max ='
     &,2F12.1)
 1100 FORMAT(//,6X,'W',7X,'x',7X,'y',8X,'Q**2',3X,'alpha',
     &5X,'cut',2X,'q-event',1X,'qg-event',1X,'qq-event',
     &/,1X,77(1H-),/)
 1200 FORMAT(F7.1,2F8.4,1PG12.3,0PF8.2,F8.4,3F9.4)
 2020 FORMAT(' Lepton energy allowed to vary in simulation, ',/,
     &' W in table below calculated assuming max energy.',/)
 2030 FORMAT(' Grid choice, LST(19) =',I3,5X,'# grid points in x, y =',
     &2I5,/,' x-values in array XX:',/,10F8.5,/,10F8.5,/,11F8.5,
     &    /,' y-values in array WW:',/,10F7.4,/,11F7.4,/)
 2100 FORMAT(//,7X,'x',7X,'y',8X,'Q**2',6X,'W',3X,'alpha',
     &5X,'cut',2X,'q-event',1X,'qg-event',1X,'qq-event',
     &/,1X,77(1H-),/)
 2200 FORMAT(2F8.5,1PG12.3,0PF7.1,F8.2,F8.4,3F9.4)
      END
C **********************************************************************

      SUBROUTINE LSIGMX(NPAR,DERIV,DIFSIG,XF,IFLAG)

C...Calculates the negative of the differential cross-section.
C...In the generation procedure the maximum of the differential cross-
C...section is needed for weighting purposes. This maximum is found by
C...minimizing the negative differential cross-section using the MINUIT
C...routines which are then calling this routine.
C...More precisly, only the part of the cross-section formula which is
C...needed for the weighting procedure is included here.
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      DIMENSION DERIV(30),XF(30)
      DATA NCALLS/0/

      DUMMY=NPAR+DERIV(1)
      IF(IFLAG.EQ.1) NCALLS=0
      IF(IFLAG.EQ.2) WRITE(6,1000)

      DIFSIG=1.E+12
      NCALLS=NCALLS+1
      X=XF(1)
      IF(X.LT.XMIN) THEN
        DIFSIG=(XMIN-X)**2*1.E+11 
        RETURN
      ELSEIF(X.GT.XMAX) THEN 
        DIFSIG=(X-XMAX)**2*1.E+11 
        RETURN
      ENDIF
      S=PARL(21)
      PM2=PSAVE(3,2,5)**2
      Q2LOW=MAX(Q2MIN,X*YMIN*S,(W2MIN-PM2)*X/(1.-X))
      Q2UPP=MIN(Q2MAX,X*YMAX*S,(W2MAX-PM2)*X/(1.-X))
      YLOW=MAX(YMIN,Q2MIN/(S*X),(W2MIN-PM2)/(S*(1.-X)))
      YUPP=MIN(YMAX,Q2MAX/(S*X),(W2MAX-PM2)/(S*(1.-X)))
      W2LOW=MAX(W2MIN,(1.-X)*YMIN*S+PM2,Q2MIN*(1.-X)/X+PM2)
      W2UPP=MIN(W2MAX,(1.-X)*YMAX*S+PM2,Q2MAX*(1.-X)/X+PM2)
      IF(LST(31).EQ.1) THEN
        Q2=XF(2)
        IF(Q2.LT.Q2LOW) THEN
          DIFSIG=(Q2LOW-Q2)**2*1.E+11 
          RETURN
        ELSEIF(Q2.GT.Q2UPP) THEN 
          DIFSIG=(Q2-Q2UPP)**2*1.E+11 
          RETURN
        ENDIF
        Y=Q2/(PARL(21)*X)
        W2=(1.-X)*Y*PARL(21)+PM2
      ELSEIF(LST(31).EQ.2) THEN
        Y=XF(2)
        IF(Y.LT.YLOW) THEN
          DIFSIG=(YLOW-Y)**2*1.E+11 
          RETURN
        ELSEIF(Y.GT.YUPP) THEN 
          DIFSIG=(Y-YUPP)**2*1.E+11 
          RETURN
        ENDIF
        Q2=Y*X*PARL(21)
        W2=(1.-X)*Y*PARL(21)+PM2
      ELSEIF(LST(31).EQ.3) THEN
        W2=XF(2)
        IF(W2.LT.W2LOW) THEN
          DIFSIG=(W2LOW-W2)**2*1.E+11 
          RETURN
        ELSEIF(W2.GT.W2UPP) THEN 
          DIFSIG=(W2-W2UPP)**2*1.E+11 
          RETURN
        ENDIF
        Y=(W2-PM2)/((1.-X)*PARL(21))
        Q2=X*Y*PARL(21)
      ENDIF
      IF(Q2.LT.Q2LOW) THEN
        DIFSIG=(Q2LOW-Q2)**2*1.E+11 
        RETURN
      ELSEIF(Q2.GT.Q2UPP) THEN 
        DIFSIG=(Q2-Q2UPP)**2*1.E+11 
        RETURN
      ENDIF
      IF(Y.LT.YLOW) THEN
        DIFSIG=(YLOW-Y)**2*1.E+11 
        RETURN
      ELSEIF(Y.GT.YUPP) THEN 
        DIFSIG=(Y-YUPP)**2*1.E+11 
        RETURN
      ENDIF
      IF(W2.LT.W2LOW) THEN
        DIFSIG=(W2LOW-W2)**2*1.E+11 
        RETURN
      ELSEIF(W2.GT.W2UPP) THEN 
        DIFSIG=(W2-W2UPP)**2*1.E+11 
        RETURN
      ENDIF
      LST2=LST(2)
      LST(2)=-1
      CALL LEPTO
      LST(2)=LST2
      DIFSIG=0. 
      IF(LST(21).NE.0) RETURN
      DIFSIG=-PQ(17)*COMFAC

      IF(LST(3).GE.4.AND.IFLAG.EQ.3)
     &WRITE(6,1100) NCALLS,DIFSIG,X,Y,Q2,W2
      RETURN

 1000 FORMAT(' Warning: IFLAG = 2 in call to LSIGMX, which does not '
     &,'calculate derivatives.')
 1100 FORMAT(/,5X,'Terminating entry in LSIGMX after ',I5,' calls.',/,
     &5X,'Best estimate of minimum found to be ',E12.4,/,
     &5X,'located at x, y, Q**2, W**2 = ',4G10.3,/)

      END

C **********************************************************************

      SUBROUTINE LXSECT

C...Integrate differential cross-section using GADAP, RIWIAD or DIVONNE

      COMMON /LINTEG/ NTOT,NPASS
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      DOUBLE PRECISION ACC,VALUE,ERRIW,FLOW,FHIGH
      COMMON /PARAMS/ ACC,NDIM,NSUB,ITER
      COMMON /ANSWER/ VALUE,ERRIW
      COMMON /BNDLMT/ FLOW,FHIGH
      COMMON /SAMPLE/ NPOINT
      DIMENSION XMINUS(2),XPLUS(2)
      EXTERNAL DCROSS,DLOWER,DUPPER,RIWFUN
      DATA NCALL/0/

      NCALL=NCALL+1
      CALL LTIMEX(TI1)
      NTOT=0
      NPASS=0
      SIGMA=0.
      ERREST=0.
      NDIM=2
C...Parameters for RIWIAD integration.
      ACC=PARL(15)
      NSUB=100
      ITER=100
C...Parameters for DIVON integration.
      DO 20 I=1,2
      XMINUS(I)=0.
   20 XPLUS(I)=1.
      EPS=PARL(15)
      MAXNUM=50000
      FLOW=-1.D0
      FHIGH=1.D+20
      NPOINT=100
C...Additional parameters for detailed DIVON integration.
      SPRDMX=2.
      MAXPTS=50000
      JDEG=0
      NPT=1000

      SIGMA=0.
      ERREST=0.
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,1000)
      IF(LST(10).EQ.1) THEN
C...Integration using GADAP.
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,1100)
        ACCUR=PARL(15)
        IT=0
  100   IT=IT+1
        ERREST=ACCUR
        CALL GADAP2(XMIN,XMAX,DLOWER,DUPPER,DCROSS,ERREST,SIGMA)
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &  WRITE(6,1110) IT,NTOT,NPASS,SIGMA
        IF(SIGMA.GT.1.) THEN
          IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &    WRITE(6,1120) ACCUR
        ELSE
          IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &    WRITE(6,1130) ACCUR,ACCUR/MAX(1.E-22,SIGMA),PARL(15)
          ACCUR=MAX(1.E-22,SIGMA*PARL(15))
          IF(IT.LT.2) GOTO 100
        ENDIF
      ELSEIF(LST(10).EQ.2) THEN
C...Integration using RIWIAD. When RIWIAD cannot be loaded:
C...activate next two lines and deactivate RIWIAD call.
C       WRITE(6,*) ' RIWIAD not available, execution stopped.'
C       STOP
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &  WRITE(6,1200) SNGL(ACC),NSUB,ITER
        CALL RIWIAD(RIWFUN)
        SIGMA=SNGL(VALUE)
        ERREST=SNGL(ERRIW)
      ELSEIF(LST(10).EQ.3) THEN
C...Integration using simple DIVONNE. When DIVONNE cannot be loaded:
C...activate next two lines and deactivate DIVONNE call.
C       WRITE(6,*) ' DIVONNE not available, execution stopped.'
C       STOP
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &  WRITE(6,1300) EPS,MAXNUM,SNGL(FLOW),SNGL(FHIGH),NPOINT
        CALL DIVON(NDIM,XMINUS,XPLUS,EPS,MAXNUM,SIGMA,ERREST)
      ELSEIF(LST(10).EQ.4) THEN
C...Integration using detailed DIVONNE. When DIVONNE cannot be loaded:
C...activate next two lines and deactivate PARTN and INTGRL calls.
C       WRITE(6,*) ' DIVONNE not available, execution stopped.'
C       STOP
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1))
     &  WRITE(6,1400) EPS,MAXNUM,
     &  SNGL(FLOW),SNGL(FHIGH),NPOINT,SPRDMX,MAXPTS,JDEG,NPT
        CALL PARTN(NDIM,XMINUS,XPLUS,SPRDMX,MAXPTS)
        CALL INTGRL(NDIM,JDEG,NPT,SIGMA,ERREST)
      ELSE
        IF(LST(3).GE.1) WRITE(6,*) ' Warning: LST(10) = ',LST(10),
     &  ' not allowed.'
      ENDIF
      CALL LTIMEX(TI2)
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)
     &.OR.(LST(3).GE.1.AND.NPASS.EQ.0)) THEN
        WRITE(6,1500) SIGMA,ERREST,NTOT,NPASS,TI2-TI1
        IF(LST(3).GE.1.AND.NPASS.EQ.0) WRITE(6,1600)
      ENDIF
      PARL(23)=SIGMA

      RETURN
 1000 FORMAT(/,' Integration of cross section:',/,1X,28('-'))
 1100 FORMAT(5X,'using GADAP = adaptive Gaussian integration')
 1110 FORMAT(5X,'Iteration #',I3,/,
     &10X,'# function evaluations; total & non-zero =',2I8,/,
     &10X,'sigma =',G10.2,' pb')
 1120 FORMAT(10X,'required  relative error = ',G10.2)
 1130 FORMAT(10X,'effective absolute error = ',G10.2,/,
     &       10X,'effective relative error = ',G10.2,/,
     &       10X,'required  relative error = ',G10.2)
 1200 FORMAT(5X,'using RIWIAD with parameters: rel. acc. = ',F10.4,
     &/,5X,'# of subvolumes = ',I5,5X,'max # iterations = ',I5)
 1300 FORMAT(5X,'using automatic DIVONNE with parameters: ',
     &'rel. acc. = ',F10.4,/,5X,'max # function calls = ',I5,
     &/,5X,'lower and upper bound on integrand =',2E12.4,
     &/,5X,'# sample points/region =',I5)
 1400 FORMAT(5X,'using detailed DIVONNE with parameters: ',
     &'rel. acc. = ',F10.4,/,5X,'max # function calls = ',I5,
     &/,5X,'lower and upper bound on integrand =',2E12.4,
     &/,5X,'# sample points/region =',I5,
     &/,5X,'SPRDMX, MAXPTS, JDEG, NPT =',F5.2,3I10)
 1500 FORMAT(/,' ===> Cross-section =',1P,G12.3,
     &' pb,  error estimate = ',G12.3,/,
     &6X,'# of integrand evaluations; total & non-zero =',2I8,/,
     &6X,'cpu time for integration =',G12.3,' seconds',/)
 1600 FORMAT(' Warning: integrand always zero, probably no allowed',
     &' phase space due to cuts',/,
     &10X,'check, in particular, CUT(11) to CUT(14)')
      END

C **********************************************************************

      SUBROUTINE RIWIBD
C   BLOCK DATA SUBSTITUTE from RIWIAD
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/STORE/XA(11),XB(11),XC(11),XD(11),MA(11),MB(11),MC(11)
      COMMON/STORE1/R(10000),LR
      COMMON/OPTION/IPRRIW,ICONV,IRESET
      COMMON/RANDOM/NSHOTS
      COMMON/INTERN/FACTOR,ALFA,BETA,GAMMA,DELTA,LEVEL,NMIN
      COMMON /LPFLAG/ LST3
      DATA INIT/0/
      IF(INIT.EQ.1) RETURN
      INIT=1
      MA(1)=0
      LR=10000
      ICONV=1
      IRESET=0
      NSHOTS=2
      FACTOR=1.65
      LEVEL=90
      ALFA=0.3
      BETA=0.3
      GAMMA=0.3
      DELTA=.7
      NMIN=2
C...Print flag to be changed here.
      IPRRIW=0
      IF(LST3.GE.4) WRITE(6,1000) IPRRIW
      RETURN
 1000 FORMAT(5X,'RIWIAD print flag changed: IPRRIW =',I5)
      END

C **********************************************************************

      SUBROUTINE DVNOPT

C...Change of default options in DIVONNE

      COMMON /PRINT/ IPRDIV
      COMMON /LPFLAG/ LST3
      IPRDIV=0
      IF(LST3.GE.2) IPRDIV=1000
      IF(LST3.GE.4) IPRDIV=10
      IF(LST3.GE.4) WRITE(6,1000) IPRDIV
      RETURN
 1000 FORMAT(5X,'DIVON4 print flag changed: IPRDIV =',I5)
      END

C **********************************************************************

      DOUBLE PRECISION FUNCTION DFUN(NDIM,X)
      INTEGER NDIM
      DOUBLE PRECISION X(NDIM)
      DOUBLE PRECISION RIWFUN !Ahmed

      DFUN=RIWFUN(X)
      RETURN
      END

C **********************************************************************

      DOUBLE PRECISION FUNCTION RIWFUN(V)
      DOUBLE PRECISION V(2)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      DATA V2MIN,V2MAX/2*0./

      RIWFUN=0.D0
      V1MIN=XMIN
      V1MAX=XMAX
      IF(LST(31).EQ.1) THEN
        V2MIN=Q2MIN
        V2MAX=Q2MAX
      ELSEIF(LST(31).EQ.2) THEN
        V2MIN=YMIN
        V2MAX=YMAX
      ELSEIF(LST(31).EQ.3) THEN
        V2MIN=W2MIN
        V2MAX=W2MAX
      ENDIF
      V1=V1MIN+V(1)*(V1MAX-V1MIN)
      V2=V2MIN+V(2)*(V2MAX-V2MIN)
      RIWFUN=DCROSS(V1,V2)*(V1MAX-V1MIN)*(V2MAX-V2MIN)

      RETURN
      END

C **********************************************************************

      FUNCTION DCROSS(V1,V2)

C...Differential cross-section dsigma/dv1dv2; v1=x, v2=Q2 or y or W2.
C...Used for numerical integration etc.
C...Note, non-zero result only for region defined by cuts through CUT.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON /LINTEG/ NTOT,NPASS
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)

      DCROSS=0.
      NTOT=NTOT+1
C...Variable V1 is x, variable V2 is either Q**2, y or W**2
      X=V1
      IF(X.LT.XMIN.OR.X.GT.XMAX) RETURN
      S=PARL(21)
      PM2=PSAVE(3,2,5)**2
      IF(LST(31).EQ.1) THEN
        Q2=V2
        Y=Q2/(PARL(21)*X)
        W2=(1.-X)*Y*PARL(21)+PSAVE(3,2,5)**2
      ELSEIF(LST(31).EQ.2) THEN
        Y=V2
        Q2=Y*X*PARL(21)
        W2=(1.-X)*Y*PARL(21)+PSAVE(3,2,5)**2
      ELSEIF(LST(31).EQ.3) THEN
        W2=V2
        Y=(W2-PSAVE(3,2,5)**2)/((1.-X)*PARL(21))
        Q2=X*Y*PARL(21)
      ENDIF
      Q2LOW=MAX(Q2MIN,X*YMIN*S,(W2MIN-PM2)*X/(1.-X))
      Q2UPP=MIN(Q2MAX,X*YMAX*S,(W2MAX-PM2)*X/(1.-X))
      YLOW=MAX(YMIN,Q2MIN/(S*X),(W2MIN-PM2)/(S*(1.-X)))
      YUPP=MIN(YMAX,Q2MAX/(S*X),(W2MAX-PM2)/(S*(1.-X)))
      W2LOW=MAX(W2MIN,(1.-X)*YMIN*S+PM2,Q2MIN*(1.-X)/X+PM2)
      W2UPP=MIN(W2MAX,(1.-X)*YMAX*S+PM2,Q2MAX*(1.-X)/X+PM2)
      IF(Q2.LT.Q2LOW.OR.Q2.GT.Q2UPP) RETURN
      IF(Y.LT.YLOW.OR.Y.GT.YUPP) RETURN
      IF(W2.LT.W2LOW.OR.W2.GT.W2UPP) RETURN
      LST2=LST(2)
      LST(2)=-2
      CALL LEPTO
      LST(2)=LST2
      IF(LST(21).NE.0) RETURN
      NPASS=NPASS+1
      DCROSS=PARI(31)*PQ(17)*COMFAC

      RETURN
      END

C **********************************************************************

      FUNCTION DLOWER(V1)

C...Lower limit on second variable (y, Q**2 or W**2) depending on first
C...variable x=V1. Used for integrating differential cross-section.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
C...CMS energy squared and target nucleon mass.
      S=PARL(21)
      PM2=PSAVE(3,2,5)**2
      IF(LST(31).EQ.1) THEN
        DLOWER=MAX(Q2MIN,V1*YMIN*S,(W2MIN-PM2)*V1/MAX(1.-V1,1.E-22))
      ELSEIF(LST(31).EQ.2) THEN
        DLOWER=MAX(YMIN,Q2MIN/(S*V1),(W2MIN-PM2)/MAX(S*(1.-V1),1.E-22))
      ELSEIF(LST(31).EQ.3) THEN
        DLOWER=MAX(W2MIN,(1.-V1)*YMIN*S+PM2,
     &  Q2MIN*(1.-V1)/MAX(V1,1.E-22)+PM2)
      ENDIF
      RETURN
      END

C **********************************************************************

      FUNCTION DUPPER(V1)

C...Upper limit on second variable (y, Q**2 or W**2) depending on first
C...variable x=V1. Used for integrating differential cross-section.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
C...CMS energy squared and target nucleon mass.
      S=PARL(21)
      PM2=PSAVE(3,2,5)**2
      IF(LST(31).EQ.1) THEN
        DUPPER=MIN(Q2MAX,V1*YMAX*S,(W2MAX-PM2)*V1/MAX(1.-V1,1.E-22))
      ELSEIF(LST(31).EQ.2) THEN
        DUPPER=MIN(YMAX,Q2MAX/(S*V1),(W2MAX-PM2)/MAX(S*(1.-V1),1.E-22))
      ELSEIF(LST(31).EQ.3) THEN
        DUPPER=MIN(W2MAX,(1.-V1)*YMAX*S+PM2,
     &  Q2MAX*(1.-V1)/MAX(V1,1.E-22)+PM2)
      ENDIF
      RETURN
      END

C **********************************************************************

      SUBROUTINE FLTABL

C...Integrates the longitudinal structure function, store on grid
C...in x, Q**2.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LINTEG/ NTOT,NPASS
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     &FLMT(41,16)
      EXTERNAL FLQINT,FLGINT,FLTINT

      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      LHT=LST(11)/100
      IF(LST(3).GE.3) WRITE(6,1000) LST(11),LQCD,LTM,LHT
      IF(LQCD.LT.1.AND.LTM.LT.1) RETURN
      CALL LTIMEX(T1)
      DO 10 IX=1,NFX
      DO 10 IQ=1,NFQ
      FLQT(IX,IQ)=0.
      FLGT(IX,IQ)=0.
   10 FLMT(IX,IQ)=0.
      QR(1)=Q2MIN
      XR(1)=XMIN
      XR(2)=XMAX
      DO 500 IX=1,NFX
      X=10**(ALOG10(XR(1))+(ALOG10(XR(2))-ALOG10(XR(1)))*(IX-1)/(NFX-1))
      QR(2)=X*PARL(21)
      IF(QR(1).GT.QR(2)) GOTO 500
      LQ=0
      DO 400 IQ=1,NFQ
      Q2=10**(ALOG10(QR(1))+(ALOG10(QR(2))-ALOG10(QR(1)))*
     &(IQ-1)/(NFQ-1))
Ctest IF(LQ.GT.0) GOTO 500
      IF(Q2.GT.PARL(21)) LQ=LQ+1
      Y=Q2/(PARL(21)*X)
      IF(Y.LT.0.0.OR.Y.GT.1.0) LQ=LQ+1
      PARL(25)=ULALPS(Q2)
      IF(LQCD.EQ.1) THEN
C...Quark part.
        ACCUR=PARL(11)
        IT=0
  100   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLQINT,EPS,FLQ)
        IF(FLQ.LT.1) THEN
          ACCUR=FLQ*PARL(11)
          IF(IT.LT.2) GOTO 100
        ENDIF
        FLQT(IX,IQ)=FLQ
C...Gluon part.
        ACCUR=PARL(11)
        IT=0
  200   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLGINT,EPS,FLG)
        IF(FLG.LT.1.) THEN
          ACCUR=FLG*PARL(11)
          IF(IT.LT.2) GOTO 200
        ENDIF
        FLGT(IX,IQ)=FLG
      ENDIF
      IF(LTM.EQ.1) THEN
C...Target mass  part.
        ACCUR=PARL(11)
        IT=0
  300   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLTINT,EPS,FLM)
        IF(FLM.LT.1) THEN
          ACCUR=FLM*PARL(11)
          IF(IT.LT.2) GOTO 300
        ENDIF
        FLMT(IX,IQ)=FLM
      ENDIF
  400 CONTINUE
  500 CONTINUE
  600 CONTINUE
      CALL LTIMEX(T2)
      IF(LST(3).GE.3) WRITE(6,1100) T2-T1
      RETURN

 1000 FORMAT(' Initialisation for FL; QCD, target mass, higher twist: ',
     &/,' LST(11) =',I5,' --> LQCD, LTM, LHT =',3I3)
 1100 FORMAT(' FL integrations performed if LQCD=1 and/or LTM=1, ',
     &'results on grid.'/,' Time for FL integrations is ',F7.1,' sec.')
      END

C **********************************************************************

      SUBROUTINE FLIPOL(FLQ,FLG,FLM)

C...QCD and target mass contributions to longitudinal structure function
C...from interpolation on x,Q2 grid.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     &FLMT(41,16)
      DATA NOUT/0/,NWARN/10/

      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      LHT=LST(11)/100
      XP=X
      Q2P=Q2
C...NOTE: tiny mismatch between present x-value and those on grid.
      QR(2)=X*PARL(21)
      IF(QR(1).GT.QR(2)) RETURN
      IF(X.LT.XR(1).OR.X.GT.XR(2).OR.
     &Q2.LT.QR(1).OR.Q2.GT.QR(2)) THEN
C...x and/or Q2 outside grid limits, write warning for first NWARN cases
        IF(LST(2).GE.0) THEN
          NOUT=NOUT+1
          IF(LST(3).GE.1.AND.NOUT.LE.NWARN) WRITE(6,1000) X,Q2,NWARN
        ENDIF
        IF(X.LT.XR(1)) XP=XR(1)
        IF(X.GT.XR(2)) XP=XR(2)
        IF(Q2.LT.QR(1)) Q2P=QR(1)
        IF(Q2.GT.QR(2)) Q2P=QR(2)
      ENDIF

      IX=(ALOG10(XP)-ALOG10(XR(1)))/
     &(ALOG10(XR(2))-ALOG10(XR(1)))*(NFX-1)+1
      IQ=(ALOG10(Q2P)-ALOG10(QR(1)))/
     &(ALOG10(QR(2))-ALOG10(QR(1)))*(NFQ-1)+1
      IX=MIN(IX,NFX-1)
      IQ=MIN(IQ,NFQ-1)
      Q2L=10**(ALOG10(QR(1))+(ALOG10(QR(2))-ALOG10(QR(1)))*
     &(IQ-1)/(NFQ-1))
      Q2H=10**(ALOG10(QR(1))+(ALOG10(QR(2))-ALOG10(QR(1)))*
     &(IQ  )/(NFQ-1))
      XL=10**(ALOG10(XR(1))+(ALOG10(XR(2))-ALOG10(XR(1)))*
     &(IX-1)/(NFX-1))
      XH=10**(ALOG10(XR(1))+(ALOG10(XR(2))-ALOG10(XR(1)))*
     &(IX  )/(NFX-1))
      QD=(Q2P-Q2L)/(Q2H-Q2L)
      XD=(XP-XL)/(XH-XL)

      IF(LQCD.EQ.1) THEN
        X1P=(FLQT(IX+1,IQ)-FLQT(IX,IQ))*XD+FLQT(IX,IQ)
        X2P=(FLQT(IX+1,IQ+1)-FLQT(IX,IQ+1))*XD+FLQT(IX,IQ+1)
        FLQ=(X2P-X1P)*QD+X1P
        X1P=(FLGT(IX+1,IQ)-FLGT(IX,IQ))*XD+FLGT(IX,IQ)
        X2P=(FLGT(IX+1,IQ+1)-FLGT(IX,IQ+1))*XD+FLGT(IX,IQ+1)
        FLG=(X2P-X1P)*QD+X1P
      ENDIF
      IF(LTM.EQ.1) THEN
        X1P=(FLMT(IX+1,IQ)-FLMT(IX,IQ))*XD+FLMT(IX,IQ)
        X2P=(FLMT(IX+1,IQ+1)-FLMT(IX,IQ+1))*XD+FLMT(IX,IQ+1)
        FLM=(X2P-X1P)*QD+X1P
      ENDIF

      RETURN
 1000 FORMAT(' Warning: x=',F7.4,' or Q2=',F6.1,' outside grid,',
     &' for FL interpolation',/,10X,'value on grid limit used.',
     &' Only first',I5,' warnings printed.',/)
      END

C **********************************************************************

      SUBROUTINE FLINTG(CFLQ,CFLG,CFLM)

C...Event-by-event calculation of contribution to longitudinal
C...structure function from QCD and target mass effects.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTEG/ NTOT,NPASS
      EXTERNAL FLQINT,FLGINT,FLTINT

      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      LHT=LST(11)/100
      PARL(25)=ULALPS(Q2)
      IF(LQCD.EQ.2) THEN
C...FL from QCD, quark and gluon contributions.
        ACCUR=PARL(11)
        IT=0
  100   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLQINT,EPS,CFLQ)
        IF(CFLQ.LT.1) THEN
          ACCUR=CFLQ*PARL(11)
          IF(IT.LT.2) GOTO 100
        ENDIF
        ACCUR=PARL(11)
        IT=0
  200   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLGINT,EPS,CFLG)
        IF(CFLG.LT.1.) THEN
          ACCUR=CFLG*PARL(11)
          IF(IT.LT.2) GOTO 200
        ENDIF
      ENDIF
      IF(LTM.EQ.2) THEN
        ACCUR=PARL(11)
        IT=0
  300   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLTINT,EPS,CFLM)
        IF(CFLM.LT.1.) THEN
          ACCUR=CFLM*PARL(11)
          IF(IT.LT.2) GOTO 300
        ENDIF
      ENDIF

      RETURN
      END

C **********************************************************************

      FUNCTION FLQINT(Z)

C...Quark contribution integrand to QCD longitudinal structure function.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
      DIMENSION XPQ(-6:6)
      DATA PI/3.14159/
      NTOT=NTOT+1
      CALL LNSTRF(Z,Q2,XPQ)
      FLQINT=0.
      DO 100 I=-LST(12),LST(12)
      IF(I.EQ.0) GOTO 100
      FLQINT=FLQINT+QC(IABS(I))**2*XPQ(I)
  100 CONTINUE
      FLQINT=4./3.*PARL(25)/PI*(X/Z)**2*FLQINT/Z
      NPASS=NPASS+1

      RETURN
      END

C **********************************************************************

      FUNCTION FLGINT(Z)

C...Gluon contribution integrand to QCD longitudinal structure function.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
      DIMENSION XPQ(-6:6)
      DATA PI/3.14159/
      NTOT=NTOT+1
      CALL LNSTRF(Z,Q2,XPQ)
      FLGINT=20./9.*PARL(25)/PI*(X/Z)**2*(1.-X/Z)/Z*XPQ(0)
      NPASS=NPASS+1

      RETURN
      END

C **********************************************************************

      FUNCTION FLTINT(Z)

C...Integrand for target mass correction contribution to
C...quark longitudinal structure function

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
      DIMENSION XPQ(-6:6)
      DATA PM2/0.8804/
      NTOT=NTOT+1
      CALL LNSTRF(Z,Q2,XPQ)
      FLTINT=0.
      DO 100 I=-LST(12),LST(12)
      IF(I.EQ.0) GOTO 100
      FLTINT=FLTINT+QC(IABS(I))**2*XPQ(I)
  100 CONTINUE
      FLTINT=4.*PM2/Q2*(X/Z)**2*X*FLTINT
      NPASS=NPASS+1

      RETURN
      END
C **********************************************************************

      SUBROUTINE LSCI(PROB)
C--                                                                   --C
C--   Created:        950319                                          --C
C--   Last update:    960730                                          --C
C--   Purpose:        to generate random switches of parton           --C
C--                   colours in the partonic final state             --C
C--                   Two versions: LST(34).EQ.1 original             --C
C--                                 LST(34).NE.1 no switch between    --C
C--                                              perturbative partons --C

      IMPLICIT NONE

C--       global variables
      COMMON /LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      INTEGER N,K
      REAL P,V
      COMMON/LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      REAL CUT,PARL,X,Y,W2,Q2,U
      INTEGER LST

C--       functions
      REAL RLU
C--       local variables
      INTEGER I,J,LUCOMP,NS,NEXT,THIS,INIT
      LOGICAL QUARK,QUARK1,QUARK2,AQUARK1,AQUARK2,GLUON1,GLUON2,FIRST
      REAL PROB

C--       Assign colour and anticolour pointers to all partons. Colour
C--       pointers are in K(I,4) and anticolour pointers are in K(I,5).
C--       The pointer points to the row where the respective anticolour
C--       and colour is.

      FIRST=.TRUE.
      DO 10 I=5,N
         IF (K(I,1).LT.10 .AND. K(I,1).GT.0) THEN
C--       check if parton is a quark, antiquark or diquark
            IF (ABS(K(I,2)).LT.10 .OR. LUCOMP(K(I,2)).EQ.90) THEN
               IF (K(I,2).LT.10 .AND. K(I,2).GT.0 .OR.
     &             K(I,2).LT.-1000) THEN
        	  QUARK=.TRUE.
               ELSE
        	  QUARK=.FALSE.
               ENDIF
C--       reset pointers
              K(I,4)=0
              K(I,5)=0
C--       the first quark, antiquark or diquark in a string points
C--       to the parton in the next line
               IF (FIRST) THEN
        	  IF (QUARK) THEN
                     K(I,4)=(I+1)
        	  ELSE
                     K(I,5)=(I+1)
        	  ENDIF
        	  FIRST=.FALSE.
C--       the last quark, antiquark or diquark in a string points
C--       to the parton in the previous line
               ELSE
        	  IF (QUARK) THEN
                     K(I,4)=(I-1)
        	  ELSE
                     K(I,5)=(I-1)
        	  ENDIF
        	  FIRST=.TRUE.
               ENDIF
               K(I,1)=3
C--       check if parton gluon 
            ELSEIF (K(I,2).EQ.21) THEN
C--       if the previous colour points to this gluon then its anticolour
C--       should point back and its colour should point to the next line 
               IF(K(I-1,4).EQ.I) THEN
        	  K(I,4)=(I+1)
        	  K(I,5)=(I-1)
               ELSE
        	  K(I,4)=(I-1)
        	  K(I,5)=(I+1)
               ENDIF
               K(I,1)=3
            ENDIF
         ENDIF
10    CONTINUE

C--       find first parton in colour switch
      DO 20 I=5,N
         QUARK1=.FALSE.
         AQUARK1=.FALSE.
         GLUON1=.FALSE.
         IF (K(I,1).EQ.3) THEN 
C--       check if parton quark or antidiquark
            IF (K(I,4).NE.0 .AND. K(I,5).EQ.0) THEN
               QUARK1=.TRUE.
C--       check if parton antiquark or diquark
            ELSEIF (K(I,4).EQ.0 .AND. K(I,5).NE.0) THEN
               AQUARK1=.TRUE.
C--       check if parton gluon 
            ELSEIF (K(I,2).EQ.21) THEN
               GLUON1=.TRUE.
            ENDIF
C--       find second parton in colour switch
            DO 30 J=I+1,N
               QUARK2=.FALSE.
               AQUARK2=.FALSE.
               GLUON2=.FALSE.
               IF (K(J,1).EQ.3 .AND.
C--       at least one remnant parton if LST(34).NE.1
     &            (K(J,3).EQ.2 .OR. K(I,3).EQ.2 .OR. LST(34).EQ.1))THEN
C--       check if second parton quark or antidiquark
                  IF (K(J,4).NE.0 .AND. K(J,5).EQ.0) THEN
                     QUARK2=.TRUE.
C--       check if second parton antquark or diquark
                  ELSEIF (K(J,4).EQ.0 .AND. K(J,5).NE.0) THEN
                     AQUARK2=.TRUE.
C--       check if second parton gluon 
                  ELSEIF (K(J,2).EQ.21) THEN
        	     GLUON2=.TRUE.
                  ENDIF
C--       switch colour pointers
                  IF (QUARK1.AND.QUARK2) THEN
                     IF (RLU(0).LT.PROB) CALL LECSWI(I,J)
                  ELSEIF (K(I,4).NE.J .AND. K(J,4).NE.I .AND. 
     &            (QUARK1.AND.GLUON2 .OR. GLUON1.AND.QUARK2)) THEN
                     IF (RLU(0).LT.PROB) CALL LECSWI(I,J)
                  ELSEIF (AQUARK1.AND.AQUARK2) THEN
                     IF (RLU(0).LT.PROB) CALL LEASWI(I,J)
                  ELSEIF (K(I,5).NE.J .AND. K(J,5).NE.I .AND. 
     &            (AQUARK1.AND.GLUON2 .OR. GLUON1.AND.AQUARK2)) THEN
                     IF (RLU(0).LT.PROB) CALL LEASWI(I,J)
                  ELSEIF (K(I,4).NE.J .AND. K(J,4).NE.I .AND.
     &                 GLUON1.AND.GLUON2) THEN
                     IF (RLU(0).LT.PROB) CALL LECSWI(I,J)
                     IF (RLU(0).LT.PROB) CALL LEASWI(I,J)
                  ENDIF
               ENDIF
30          CONTINUE
         ENDIF
20    CONTINUE

C--       restore colour order in strings ready for hadronisation
      NS=N
      DO 40 I=5,NS
C--       find first quark (or anti diquark) string end
         IF (K(I,1).EQ.3 .AND. K(I,4).NE.0 .AND. K(I,5).EQ.0 )THEN
            NEXT=I
50          CONTINUE
            N=N+1
            IF(N.GT.4000) THEN 
              IF(LST(3).GE.1) WRITE(6,*) 'LSCI: N>4000!'
              LST(21)=101
              RETURN
            ENDIF
            THIS=NEXT
C--       copy to last row in event-record and update K-vector
            DO 60 J=1,5
               P(N,J)=P(THIS,J)
               V(N,J)=V(THIS,J)
               K(N,J)=K(THIS,J)
60          CONTINUE
            K(THIS,1)=13
            K(N,1)=2
            K(N,3)=THIS
            K(N,4)=0
            K(N,5)=0
C--       find next parton in string in row K(THIS,4)
            NEXT=K(THIS,4)
            IF (NEXT.NE.0) GOTO 50
C--       this is the last parton in string
            K(N,1)=1    
         ENDIF
40    CONTINUE
      DO 70 I=5,NS
C--       find first gluon string end
         IF (K(I,1).EQ.3 .AND. K(I,2).EQ.21) THEN
            INIT=I
            NEXT=I
80          CONTINUE
            N=N+1
            IF(N.GT.4000) THEN 
              IF(LST(3).GE.1) WRITE(6,*) 'LSCI: N>4000!'
              LST(21)=101
              RETURN
            ENDIF
            THIS=NEXT
C--       copy to last row in event-record and update K-vector
            DO 90 J=1,5
               P(N,J)=P(THIS,J)
               V(N,J)=V(THIS,J)
               K(N,J)=K(THIS,J)
90          CONTINUE
            K(THIS,1)=13
            K(N,1)=2
            K(N,3)=THIS
            K(N,4)=0
            K(N,5)=0
C--       find next parton in string in row K(THIS,4)
            NEXT=K(THIS,4)
            IF (NEXT.NE.INIT) GOTO 80
C--       this is the last parton in string
            K(N,1)=1    
         ENDIF
70    CONTINUE

      END


C ********************************************************************

      SUBROUTINE LEASWI(I,J)

C--       switch anticolour pointers for partons in rows I and J
C--       and colour pointers for the partons pointing back

      IMPLICIT NONE

C--       global variables
      INTEGER N,K
      REAL P,V
      COMMON /LUJETS/N,K(4000,5),P(4000,5),V(4000,5)

C--       local variables
      INTEGER I,J,KI5,KJ5

      KI5=K(I,5)
      KJ5=K(J,5)
      K(KI5,4)=J
      K(KJ5,4)=I
      K(I,5)=KJ5
      K(J,5)=KI5         

      END

C ********************************************************************

      SUBROUTINE LECSWI(I,J)

C--       switch colour pointers for partons in rows I and J 
C--       and anticolour pointers for the partons pointing back

      IMPLICIT NONE

C--       global variables
      INTEGER N,K
      REAL P,V
      COMMON /LUJETS/N,K(4000,5),P(4000,5),V(4000,5)

C--       local variables
      INTEGER I,J,KI4,KJ4

      KI4=K(I,4)
      KJ4=K(J,4)
      K(KI4,5)=J
      K(KJ4,5)=I
      K(I,4)=KJ4
      K(J,4)=KI4

      END
C ********************************************************************

      SUBROUTINE LSMALL
C--                                                                  --C
C--   Created:     951031                                            --C
C--   Last update: 960109                                            --C
C--   Purpose:     Take care of small mass systems with one diquark  --C
C--                                                                  --C

      IMPLICIT NONE

C--       global variables
      INTEGER N,K
      REAL P,V
      COMMON /LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      INTEGER LST
      REAL CUT,PARL,X,Y,W2,Q2,U
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      INTEGER MSTU,MSTJ
      REAL PARU,PARJ
      COMMON /LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      INTEGER KCHG
      REAL PMAS,PARF,VCKM
      COMMON /LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 

C--       functions
      REAL RLU,ULMASS
      INTEGER LUCOMP
C--       local variables
      INTEGER I,IFIRST,I1,I2,IDQ,IQ,NTRY1,NTRY2,KFB,KFH1,KFH2
      INTEGER KIDQ,KIQ,KDUMMY
      REAL*8 ENERGY,MINENERGY,MAXENERGY,INVMASS,MININVMASS,MAXINVMASS
      REAL*8 PLIGHT(5),PSUM(5)
      REAL*8 TOT2,M1,M2,ROTARG,PABS,COSTHE,PTEMP,PHI,PI
      REAL*8 PCPS,PC2,PN2,PS2,A,B,C,EPS2,EPS1
      LOGICAL FIRST
      DATA PI/3.14159265359D0/

C--       find lightest singlet system
      FIRST=.TRUE.
      MINENERGY=PARL(21)
      DO 20 I=1,N
         IF (K(I,1).EQ.2) THEN
            IF (FIRST) THEN
               PSUM(1)=P(I,1)
               PSUM(2)=P(I,2)
               PSUM(3)=P(I,3)
               PSUM(4)=P(I,4)
               PSUM(5)=P(I,5)
               MSTJ(93)=1
               PSUM(5)=ULMASS(K(I,2)) 
               FIRST=.FALSE.
               IFIRST=I
            ELSE
               PSUM(1)=PSUM(1)+P(I,1)
               PSUM(2)=PSUM(2)+P(I,2)
               PSUM(3)=PSUM(3)+P(I,3)
               PSUM(4)=PSUM(4)+P(I,4)
               PSUM(5)=PSUM(5)+P(I,5)
            ENDIF
         ELSEIF (K(I,1).EQ.1 .AND. .NOT. FIRST) THEN
            PSUM(1)=PSUM(1)+P(I,1)
            PSUM(2)=PSUM(2)+P(I,2)
            PSUM(3)=PSUM(3)+P(I,3)
            PSUM(4)=PSUM(4)+P(I,4)
            PSUM(5)=PSUM(5)+P(I,5)
            MSTJ(93)=1
            PSUM(5)=PSUM(5)+ULMASS(K(I,2)) 
            INVMASS=PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2
            ENERGY=SQRT(MAX(0.D0,INVMASS))-PSUM(5)
            FIRST=.TRUE.
C--       only system with a diquark are of interest
            IF (ENERGY.LT.MINENERGY .AND. 
     &         LUCOMP(K(IFIRST,2)).EQ.90 .OR. LUCOMP(K(I,2)).EQ.90) THEN
               MINENERGY=ENERGY
               MININVMASS=INVMASS
               I1=IFIRST
               I2=I
               PLIGHT(1)=PSUM(1)
               PLIGHT(2)=PSUM(2)
               PLIGHT(3)=PSUM(3)
               PLIGHT(4)=PSUM(4)
            ENDIF
         ENDIF
20    CONTINUE

C--       lightest singlet system is in row I1 to I2
      IF (MINENERGY.LT.DBLE(PARJ(32))) THEN

C--       add system as cluster
         N=N+1
         DO 30 I=1,4
            P(N,I)=PLIGHT(I)
30       CONTINUE
         P(N,5)=SQRT(MAX(0.D0,MININVMASS))
         K(N,1)=11
         K(N,2)=91
         K(N,3)=I1
         K(N,4)=N+1
         K(N,5)=N+2

C--       inactivate old system
         DO 40 I=I1,I2
            K(I,1)=K(I,1)+10
            K(I,4)=N
40       CONTINUE
         
C--       try to make two particles
         NTRY1=0
         NTRY2=0
50       CONTINUE
C--       take diquark end first
         IF (LUCOMP(K(I1,2)).EQ.90) THEN
            IDQ=I1
            IQ=I2
         ELSE
            IDQ=I2
            IQ=I1
         ENDIF
         KIDQ=K(IDQ,2)
         KIQ=K(IQ,2)
         KDUMMY=0
         KFB=0
         KFH1=0
         KFH2=0
         CALL LUKFDI(KIDQ,KDUMMY,KFB,KFH1)
         KDUMMY=0
         CALL LUKFDI(KIQ,-KFB,KDUMMY,KFH2)
         IF (KFH1.EQ.0 .OR. KFH2.EQ.0) THEN
            NTRY1=NTRY1+1
            IF (NTRY1.GE.100) THEN
               LST(21)=200
               RETURN
            ENDIF
            GOTO 50
         ENDIF
C--       consistency checks
         IF (LUCOMP(KFH1).EQ.0 .OR. LUCOMP(KFH2).EQ.0) THEN
            LST(21)=201
            RETURN
         ENDIF
         IF (KCHG(LUCOMP(KFH1),2)*ISIGN(1,KFH1)+
     +       KCHG(LUCOMP(KFH2),2)*ISIGN(1,KFH2) .NE.0) THEN
            LST(21)=202
            RETURN
         ENDIF
         P(N+1,5)=ULMASS(KFH1)
         P(N+2,5)=ULMASS(KFH2)
         IF (P(N,5).LE.P(N+1,5)+P(N+2,5)+PARJ(64) .AND.
     &       P(N,5).GE.ULMASS(2212)+ULMASS(111)+PARJ(64) .AND.
     &       NTRY2.LE.100) THEN
            NTRY2=NTRY2+1
            GOTO 50
         ENDIF
         IF (P(N,5).GE.P(N+1,5)+P(N+2,5)+PARJ(64)) THEN
C--       make two particles
C--       isotropic decay in cms (dcostheta*dphi)
            TOT2=MININVMASS
            M1=DBLE(P(N+1,5))
            M2=DBLE(P(N+2,5))            
            ROTARG=(TOT2-M1**2-M2**2)**2-4.D0*M1**2*M2**2
            IF (ROTARG.LT.0.D0) THEN
               LST(21)=203
               RETURN
            ENDIF
            PABS=0.5D0*SQRT(ROTARG/TOT2)
            COSTHE=-1.D0+2.D0*RLU(0)
            PTEMP=PABS*SQRT(1.D0-COSTHE**2)
            PHI=2.D0*PI*RLU(0)
            P(N+1,4)=SQRT(PABS**2+M1**2)
            P(N+1,3)=PABS*COSTHE
            P(N+1,2)=PTEMP*COS(PHI)
            P(N+1,1)=PTEMP*SIN(PHI)
            P(N+2,4)=SQRT(PABS**2+M2**2)
            P(N+2,3)=-P(N+1,3)
            P(N+2,2)=-P(N+1,2)
            P(N+2,1)=-P(N+1,1)
C--       K-vector
            K(N+1,1)=1
            K(N+1,2)=KFH1
            K(N+1,3)=N
            K(N+1,4)=0
            K(N+1,5)=0
            K(N+2,1)=1
            K(N+2,2)=KFH2
            K(N+2,3)=N
            K(N+2,4)=0
            K(N+2,5)=0
C--       boost to lab
            MSTU(33)=1
            CALL LUDBRB(N+1,N+2,0.,0.,PLIGHT(1)/PLIGHT(4),
     &               PLIGHT(2)/PLIGHT(4),PLIGHT(3)/PLIGHT(4)) 
C--       V-vector
            DO 60 I=1,5 
               V(N,I)=V(IDQ,I) 
               V(N+1,I)=V(IDQ,I) 
               V(N+2,I)=V(IQ,I) 
60          CONTINUE 
            V(N,5)=0. 
            V(N+1,5)=0. 
            V(N+2,5)=0.
            N=N+2 
         ELSE
C--       make one particle instead
70          CONTINUE
            KIDQ=K(IDQ,2)
            KIQ=K(IQ,2)
            KDUMMY=0
            KFH1=0
            CALL LUKFDI(KIDQ,KIQ,KDUMMY,KFH1)
            IF (KFH1.EQ.0) GOTO 70
C--       isospin conservation
            IF (KFH1.EQ.2214) KFH1=2212
            IF (KFH1.EQ.2114) KFH1=2112
            IF (KFH1.EQ.-2214) KFH1=-2212
            IF (KFH1.EQ.-2114) KFH1=-2112
            K(N+1,1)=1
            K(N+1,2)=KFH1
            K(N+1,3)=N
            K(N+1,4)=0
            K(N+1,5)=0
            P(N+1,5)=ULMASS(KFH1)
C--       find particle to shuffle energy & momentum to and from
            MAXENERGY=0.D0
            I1=0
            DO 80 I=1,N-1
               IF (0.LT.K(I,1) .AND. K(I,1).LT.10 .AND.
     &             LUCOMP(K(I,2)).NE.0 .AND.
     &             (K(I,2).EQ.21 .OR. ABS(K(I,2)).LT.10 .OR.
     &              ABS(K(I,2)).GT.100) ) THEN
                  INVMASS=(P(N,4)+P(I,4))**2-(P(N,1)+P(I,1))**2-
     -                    (P(N,2)+P(I,2))**2-(P(N,3)+P(I,3))**2
                  ENERGY=SQRT(MAX(0.D0,INVMASS))-P(N,5)-P(I,5)
                  IF (ENERGY.GT.MAXENERGY ) THEN
                     I1=I
                     MAXENERGY=ENERGY
                  ENDIF
               ENDIF
80          CONTINUE
C--       shuffle energy & momentum
            IF (I1.NE.0) THEN
              PCPS=DBLE(P(N,4))*DBLE(P(I1,4))-DBLE(P(N,1))*DBLE(P(I1,1))
     -            -DBLE(P(N,2))*DBLE(P(I1,2))-DBLE(P(N,3))*DBLE(P(I1,3))
              PC2=DBLE(P(N,5))**2
              PN2=DBLE(P(N+1,5))**2
              PS2=DBLE(P(I1,5))**2
              A=PC2+PS2+2.D0*PCPS
              B=PC2+PN2+2.D0*PCPS
              C=(PN2-PC2)*(4.D0*PCPS*(PCPS+PC2)-PC2*(PN2-PC2))/
     /          4.D0/(PCPS**2-PC2*PS2)
              IF (B**2-4.D0*C*A.LT.0.D0) THEN
                 LST(21)=204
                 RETURN
              ENDIF
              EPS2=(-B+SQRT(MAX(0.D0,B**2-4.D0*C*A)))/2.D0/A
              EPS1=(PN2-PC2+2.D0*EPS2*(PS2+PCPS))/2.D0/(PC2+PCPS)
              DO 90 I=1,4 
                 P(N+1,I)=(1.+EPS1)*P(N,I)-EPS2*P(I1,I) 
                 P(I1,I)=(1.+EPS2)*P(I1,I)-EPS1*P(N,I) 
                 V(N,I)=V(I1,I) 
                 V(N+1,I)=V(I1,I) 
90            CONTINUE 
              V(N,5)=0. 
              V(N+1,5)=0. 
              N=N+1
            ELSE
              LST(21)=205
              RETURN
            ENDIF
         ENDIF
      ENDIF
      
      RETURN
      
      END

C#######################################################################
C
C  The following routines for parton cascades were made together
C  with M. Bengtsson and T. Sjostrand (Z. Phys. C37 (1988) 465,
C  Nucl. Phys. B301 (1988) 554). Contain modifications of
C  routines in PYTHIA 4.8 (Sjostrand, Bengtsson, CPC 46 (1987) 43).
C
C **********************************************************************

      SUBROUTINE LSHOWR(ICALL)

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
      COMMON /LYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON /LYPROC/ ISUB,KFL(3,2),XPY(2),SH,TH,UH,Q2PY,XSEC(0:40)
      COMMON /LYINT1/ XQPY(2,-6:6),DSIG(-6:6,-6:6,5),FSIG(10,10,3)
      DOUBLE PRECISION DTHETA,DPHI,DBETA
      DOUBLE PRECISION DPQ2,DPB(3),DPA(3),DCTHET,DROBO(5)
      DIMENSION KS(9,5),PS(9,5),ROBO(5),XPQ(-6:6)
      SAVE KS,PS

      IF(ICALL.EQ.0) THEN
C...Initialize cascade for each event, save event record in overall cms.
        DO 10 I=1,9
        DO 10 J=1,5
        KS(I,J)=0
   10   PS(I,J)=0.
        DO 20 J=1,5
        KS(1,J)=K(1,J)
        PS(1,J)=P(1,J)
        KS(2,J)=K(2,J)
        PS(2,J)=P(2,J)
        KS(5,J)=K(3,J)
        PS(5,J)=P(3,J)
        KS(7,J)=K(4,J)
   20   PS(7,J)=P(4,J)
        KS(5,3)=3
        KS(7,1)=21
        KS(7,3)=5
C       CALL GULIST(0,2)
        RETURN
      ENDIF

C     CALL GULIST(1,2)
C...Apply parton cascade on QPM event.
C...Save incoming and outgoing quark as well as scattered lepton.
      KS(6,1)=21
      KS(6,2)=LST(25)
      KS(6,3)=4
      KS(8,1)=21
      KS(8,2)=K(5,2)
      KS(8,3)=6
      KS(9,1)=0
      KS(9,2)=K(4,2)
      KS(9,3)=5
      DO 110 J=1,5
      PS(6,J)=0.
      PS(8,J)=P(5,J)
  110 PS(9,J)=P(4,J)
      XR=X
      DPQ2=DBLE(Q2)
      PMA1=0.
      PS(6,5)=PMA1
      PMA2=PS(8,5)
      DPB(1)=0.5D0*(DPQ2*(1D0/XR-1D0)+DBLE(PS(1,5))**2-
     &ULMASS(IABS(KS(7,2)))**2)/(PS(1,4)+PS(2,4))
      DPB(2)=DSQRT(DPB(1)**2+DPQ2)
      DCTHET=(DBLE(PS(2,4))*DPB(1)-DPQ2/(2D0*XR))/(DBLE(PS(2,3))*
     &DPB(2))
      DPA(1)=(DPB(2)*DCTHET)**2-DPB(1)**2
      DPA(2)=DPQ2-DBLE(PMA1)**2+DBLE(PMA2)**2
      PS(6,4)=-(DPA(2)*DPB(1)-DPB(2)*DCTHET*DSQRT(DPA(2)**2+4D0*
     &DBLE(PMA1)**2*DPA(1)))/(2D0*DPA(1))
      PS(6,3)=-SQRT((PS(6,4)+PMA1)*(PS(6,4)-PMA1))
C...Partons with colour information in hadronic cms frame.
      DO 120 I=10,26
      DO 120 J=1,5
      K(I,J)=0
      P(I,J)=0.
  120 V(I,J)=0.
      NS=20
      K(NS+1,1)=21
      K(NS+1,2)=K(3,2)
      K(NS+1,3)=3
      K(NS+2,1)=-1
      K(NS+2,3)=NS+1
      K(NS+3,2)=KS(6,2)
      DO 130 J=1,5
  130 P(NS+1,J)=P(3,J)
      K(NS+3,1)=13
      K(NS+3,3)=2
      P(NS+3,5)=0.
      K(NS+4,1)=-1
      K(NS+4,3)=NS+3
      K(NS+3,4)=NS+5
      K(NS+3,5)=NS+5
      P(NS+4,3)=NS+5
      P(NS+4,4)=NS+5
      K(NS+5,1)=3
      K(NS+5,3)=8
      K(NS+5,2)=KS(8,2)
      K(NS+6,1)=-1
      K(NS+6,3)=NS+5
      DO 140 J=1,4
      P(NS+5,J)=P(5,J)
      P(NS+3,J)=P(NS+5,J)-P(NS+1,J)
  140 CONTINUE
      P(NS+5,5)=PMA2
      P(NS+6,1)=NS+3
      P(NS+6,2)=NS+3
      K(NS+5,4)=(NS+3)*MSTU(5)
      K(NS+5,5)=(NS+3)*MSTU(5)
      N=NS+6
C     CALL GULIST(2,2)
C...Copy saved record in overall cms to line 1 through 9.
C...Lines 1,2,5,6,7 in ep cms, 8,9 in hadronic cms
      DO 150 I=1,9
      DO 150 J=1,5
      K(I,J)=KS(I,J)
  150 P(I,J)=PS(I,J)
C     CALL GULIST(3,2)
C...Scale for bremsstrahlung etc.
      Q2PY=Q2
      IPY(40)=8
      IPY(47)=N
C...Save quantities for later use.
      XPY(1)=1.
      XPY(2)=XR
      CALL LYSTFU(K(2,2),XR,Q2,XPQ)
      DO 160 IFL=-6,6
  160 XQPY(2,IFL)=XPQ(IFL)
      IF(LST(23).EQ.1) THEN
        ISUB=39
        IPY(11)=1
      ELSEIF(LST(23).EQ.3) THEN
        ISUB=39
        IPY(11)=2
      ELSEIF(LST(23).EQ.4) THEN
        ISUB=39
        IPY(11)=3
      ELSEIF(LST(23).EQ.2) THEN
        ISUB=40
      ENDIF
      KFL(2,1)=K(5,2)
      IFL1=K(6,2)
      IFL2=K(8,2)
      KFL(2,2)=IFL1
      KFL(1,1)=KFL(2,1)
      KFL(1,2)=KFL(2,2)
      IF(ISUB.EQ.39) KFL(3,1)=K(1,2)
      IF(ISUB.EQ.40) KFL(3,1)=K(1,2)+ISIGN(1,K(1,2))
      KFL(3,2)=IFL2
      PYVAR(2)=(P(1,4)+P(2,4))**2
      PYVAR(1)=SQRT(PYVAR(2))
      PYVAR(3)=P(1,5)
      PYVAR(4)=P(2,5)
      PYVAR(5)=P(1,3)
      IPY(41)=K(1,2)
      IPY(42)=K(2,2)
      IPY(48)=0

C...Generate timelike parton shower (if required)
      IF(IPY(13).EQ.1) THEN
        CALL LSCALE(1,QMAX)
        QMAX=MIN(QMAX,P(25,4))
        CALL LUSHOW(25,0,QMAX)
      ENDIF
      IT=25
      IF(N.GE.27) IT=27
      NS=N
C     CALL GULIST(4,2)

C...Generate spacelike parton shower (if required)
      IPU1=0
      IPU2=23
      IF(XPY(2)*(1.+(P(IT,5)**2+PYPAR(22))/P(21,5)**2).GT.0.999) THEN
        LST(21)=7
        RETURN
      ENDIF
      IF(IPY(14).GE.1) THEN
        CALL LYSSPA(IPU1,IPU2)
        IF(LST(21).NE.0) RETURN
      ENDIF
      IF(N.EQ.NS) THEN 
        DO 220 I=NS+1,NS+4
        DO 220 J=1,5
        K(I,J)=0
        P(I,J)=0.
  220   V(I,J)=0.
        K(NS+1,1)=11
        K(NS+1,2)=KFL(2,1)
        K(NS+1,3)=21
        DO 230 J=1,5
  230   P(NS+1,J)=P(21,J)
        K(NS+2,1)=-1
        K(NS+2,3)=NS+1
        K(NS+3,1)=13
        K(NS+3,2)=KFL(2,2)
        K(NS+3,3)=23
        K(NS+3,4)=23
        K(NS+3,5)=23
        P(NS+3,3)=(P(IT,5)**2+Q2)*(P(21,4)-P(21,3))/(2.*Q2)
        P(NS+3,4)=-P(NS+3,3)
        K(NS+4,1)=-1
        K(NS+4,3)=NS+3
        P(NS+4,3)=23
        P(NS+4,4)=23
        P(24,1)=NS+3
        P(24,2)=NS+3
        K(23,4)=K(23,4)+(NS+3)*MSTU(5)
        K(23,5)=K(23,5)+(NS+3)*MSTU(5)
        IPU1=0
        IPU2=NS+3
        N=N+4
      ENDIF
C     CALL GULIST(5,2)

      IF(ABS(P(IT,1)).GT.0.1.OR.ABS(P(IT,2)).GT.0.1) THEN
C       WRITE(6,*) 'Warning: non-zero pt on final shower initiator'
C       WRITE(6,*) '0:',IT,K(IT,2),P(IT,1),P(IT,2),P(IT,3),P(IT,4),P(IT,5)
        LST(21)=8
        RETURN
      ENDIF
      P(IT,1)=0.
      P(IT,2)=0.

C...Rotate and boost outgoing parton shower
      IF(N.GT.30) THEN
        K(N+1,1)=0
        DO 210 J=1,4
  210   P(N+1,J)=P(NS+1,J)+P(NS+3,J)
        IF(P(N+1,4).LE.1.01*P(IT,5)) THEN
          LST(21)=9
          RETURN
        ENDIF
        ROBO(1)=ULANGL(P(IT,3),SQRT(P(IT,1)**2+P(IT,2)**2))
        ROBO(2)=ULANGL(P(IT,1),P(IT,2))
      IF(ABS(ROBO(1)).GT.0.001.OR.ABS(ROBO(2)).GT.0.001) THEN
      WRITE(6,*) '0:',IT,K(IT,2),P(IT,1),P(IT,2),P(IT,3),P(IT,4),P(IT,5)
      WRITE(6,*) '   ROBO(1-2)=',ROBO(1),ROBO(2)
      ENDIF
        CALL LUDBRB(25,NS,0.,-ROBO(2),0.D0,0.D0,0.D0)
        CALL LUDBRB(25,NS,-ROBO(1),0.,0.D0,0.D0,0.D0)
        DROBO(5)=-(P(IT,3)*P(IT,4)-P(N+1,4)*SQRT(P(N+1,4)**2-
     &  P(IT,4)**2+P(IT,3)**2))/(P(IT,3)**2+P(N+1,4)**2)
        CALL LUDBRB(25,NS,0.,0.,0.D0,0.D0,DROBO(5))
        ROBO(1)=ULANGL(P(N+1,3),SQRT(P(N+1,1)**2+P(N+1,2)**2))
        ROBO(2)=ULANGL(P(N+1,1),P(N+1,2))
        CALL LUDBRB(25,NS,ROBO(1),ROBO(2),0.D0,0.D0,0.D0)
      ENDIF
C     CALL GULIST(6,2)

      Q2PY=Q2
C...Hadron remnant and primordial kt
      IPY(47)=N
      CALL LYREMN(IPU1,IPU2)
      IF(IPY(48).EQ.1) THEN
        LST(21)=10
        RETURN
      ENDIF
C     CALL GULIST(7,2)

C...Transform line 1,2 and 5-7 to hadronic cms frame.
      CALL LUDBRB(1,2,0.,0.,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
      CALL LUDBRB(1,2,-STHETA(2),0.,0.D0,0.D0,0.D0)
      CALL LUDBRB(5,7,0.,0.,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
      CALL LUDBRB(5,7,-STHETA(2),0.,0.D0,0.D0,0.D0)
C     CALL GULIST(8,2)

C...Rearrange partons along strings
      MSTU(24)=0
      CALL LUPREP(0)
      IF(MSTU(24).NE.0) THEN
C       CALL GULIST(88,2)
        IF(LST(3).GE.1) WRITE(6,*) ' LUPREP error MSTU(24)= ',MSTU(24)
        LST(21)=11
        RETURN
      ENDIF
C     CALL GULIST(9,2)

C...Clean up event record -> order:
C...1=inc. lepton; 2=inc. nucleon; 3=exch boson; 4=scat. lepton;
C...5=inc. parton before initial shower; 6=inc. quark at boson vertex
C...after shower; 7=scat. quark at boson vertex before final shower
      LST(26)=7
      DO 510 J=1,5
      K(N+1,J)=K(4,J)
  510 P(N+1,J)=P(4,J)
      DO 520 J=1,5
      K(3,J)=K(5,J)
      P(3,J)=P(5,J)
      K(4,J)=K(9,J)
      P(4,J)=P(9,J)
      K(5,J)=K(N+1,J)
      P(5,J)=P(N+1,J)
C     K(7,J)=K(8,J)
C     P(7,J)=P(8,J)
      K(6,J)=K(NS+3,J)
      P(6,J)=P(NS+3,J)
      K(7,J)=K(IT,J)
      P(7,J)=P(IT,J)
  520 CONTINUE
      K(3,3)=1
      K(4,3)=1
      K(6,1)=21
      K(6,3)=5
      K(6,4)=0
      K(6,5)=0
      K(7,1)=21
      K(7,3)=6
      K(7,4)=0
      K(7,5)=0
C...Activate line with scattered lepton.
      K(4,1)=1
C...Deactivate obsolete lines 8, 9 and 21, NS+1 (extra lines with boson)
      K(8,1)=0
      K(9,1)=0
      K(21,1)=0
      IF(K(NS+1,2).EQ.K(3,2)) K(NS+1,1)=0
C...Zero irrelevant lines with K(I,1)<0
      DO 540 I=1,N
      IF(K(I,1).LT.0) THEN
        DO 530 J=1,5
        K(I,J)=0
  530   P(I,J)=0.
      ENDIF
  540 CONTINUE
C     CALL GULIST(10,2)
C...Delete internal parton lines, i.e. with K(I,1)=13,14
      IF(MOD(LST(4)/10,10).EQ.0) THEN
        CALL LTIMEX(T1)
        CALL LUEDIT(14)
        CALL LTIMEX(T2)
C       CALL GULIST(11,2)
      ENDIF
C...Delete empty lines
      CALL LTIMEX(T1)
      CALL LUEDIT(12)
      CALL LTIMEX(T2)
C     CALL GULIST(12,2)

      RETURN
      END

C **********************************************************************

      SUBROUTINE LMEPS

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
      COMMON /LYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON /LYPROC/ ISUB,KFL(3,2),XPY(2),SH,TH,UH,Q2PY,XSEC(0:40)
      COMMON /LYINT1/ XQPY(2,-6:6),DSIG(-6:6,-6:6,5),FSIG(10,10,3)
      DOUBLE PRECISION DTHETA,DPHI,DBETA
      DOUBLE PRECISION DPQ2,DPB(3),DPA(3),DCTHET,DROBO(5)
      DOUBLE PRECISION DELTAP(4),DPLONG,DBTOT,DGAMMA,DROOT
      DIMENSION KS(9,5),PS(9,5),ROBO(5),XPQ(-6:6)
      SAVE KS,PS

C     CALL GULIST(100,2)
C...Save event record in hadronic cms
      DO 10 I=1,7
      DO 10 J=1,5
      KS(I,J)=K(I,J)
   10 PS(I,J)=P(I,J)
C...Rearrange event record to PYSSPA standard
      IP2=6
      IF(LST(24).EQ.3) IP2=7
      DO 20 J=1,5
      K(3,J)=0.
      P(3,J)=0.
      K(4,J)=0
      P(4,J)=0.
      K(5,J)=KS(3,J)
      P(5,J)=PS(3,J)
      K(7,J)=KS(4,J)
      P(7,J)=PS(4,J)
      K(8,J)=KS(5,J)
      P(8,J)=PS(5,J)
      K(9,J)=KS(4,J)
      P(9,J)=PS(4,J)
      K(10,J)=KS(IP2,J)
   20 P(10,J)=PS(IP2,J)
      K(5,3)=3
      K(6,3)=4
      K(7,3)=5
      K(8,3)=6
      K(9,3)=5
      K(10,3)=6
      DO 30 I=5,10
   30 K(I,1)=21
      K(9,1)=0
C...Incoming parton = outgoing 2 parton - boson fourvectors
      DO 40 J=1,4
   40 P(6,J)=P(8,J)+P(10,J)-P(5,J)
      P(6,5)=0.
      K(6,2)=LST(25)
      IF(LST(24).EQ.3) K(6,2)=21
      N=10
C     CALL GULIST(101,2)

      XR=X
      DPQ2=DBLE(Q2)
C...Partons with colour information in hadronic cms frame.
      DO 120 I=11,27
      DO 120 J=1,5
      K(I,J)=0
      P(I,J)=0.
  120 V(I,J)=0.
      NS=20
      DO 130 J=1,5
      K(NS+1,J)=K(5,J)
      P(NS+1,J)=P(5,J)
      K(NS+3,J)=K(6,J)
      P(NS+3,J)=P(6,J)
      K(NS+5,J)=K(8,J)
      P(NS+5,J)=P(8,J)
      K(NS+6,J)=K(10,J)
  130 P(NS+6,J)=P(10,J)
C...Old standard continuation lines
      K(NS+2,1)=-1
      K(NS+2,3)=NS+1
      K(NS+4,1)=-1
      K(NS+4,3)=NS+3
      P(NS+4,3)=27
      P(NS+4,4)=27
C...Origin and colour info for incoming parton
      K(NS+3,1)=13
      K(NS+3,3)=2
      K(NS+3,4)=27
      K(NS+3,5)=27
C...Colour info for two outgoing partons
      K(NS+5,1)=3
      K(NS+6,1)=3
      IF(K(NS+6,2).EQ.21) THEN
C...qg-event
        IF(K(NS+5,2).GT.0) THEN
          K(NS+5,4)=(NS+6)*MSTU(5)
          K(NS+5,5)=(NS+7)*MSTU(5)
          K(NS+6,4)=(NS+7)*MSTU(5)
          K(NS+6,5)=(NS+5)*MSTU(5)
        ELSE
          K(NS+5,4)=(NS+7)*MSTU(5)
          K(NS+5,5)=(NS+6)*MSTU(5)
          K(NS+6,4)=(NS+5)*MSTU(5)
          K(NS+6,5)=(NS+7)*MSTU(5)
        ENDIF
      ELSE
C...qqbar-event
        K(NS+5,4)=(NS+7)*MSTU(5)
        K(NS+5,5)=(NS+7)*MSTU(5)
        K(NS+6,4)=(NS+7)*MSTU(5)
        K(NS+6,5)=(NS+7)*MSTU(5)
      ENDIF
C...Effective outgoing parton = sum of both outgoing partons
      K(NS+7,1)=14
      K(NS+7,3)=3
      IF(LST(24).EQ.2) THEN
        K(NS+7,2)=K(NS+5,2)
        IF(K(NS+7,2).EQ.21) WRITE(6,*) ' Warning: K(NS+7,2)=',K(NS+7,2)
        IF(K(NS+7,2).GT.0) THEN
          K(NS+7,4)=(NS+3)*MSTU(5)+26
          K(NS+7,5)=(NS+3)*MSTU(5)+25
        ELSE
          K(NS+7,4)=(NS+3)*MSTU(5)+25
          K(NS+7,5)=(NS+3)*MSTU(5)+26
        ENDIF
      ELSE
        K(NS+7,2)=21
        IF(K(NS+5,2).GT.0) THEN
          K(NS+7,4)=(NS+3)*MSTU(5)+25
          K(NS+7,5)=(NS+3)*MSTU(5)+26
        ELSE
          K(NS+7,4)=(NS+3)*MSTU(5)+26
          K(NS+7,5)=(NS+3)*MSTU(5)+25
        ENDIF
      ENDIF
      DO 140 J=1,4
  140 P(NS+7,J)=P(8,J)+P(10,J)

      IT=NS+7
      IF(ABS(P(IT,1)).GT.0.1.OR.ABS(P(IT,2)).GT.0.1) THEN
C       WRITE(6,*) 'Warning: non-zero pt on final shower initiator'
C       WRITE(6,*) '1:',IT,K(IT,2),P(IT,1),P(IT,2),P(IT,3),P(IT,4),P(IT,5)
C       WRITE(6,*) '1:',8 ,K( 8,2),P( 8,1),P( 8,2),P( 8,3),P( 8,4),P( 8,5)
C       WRITE(6,*) '1:',10,K(10,2),P(10,1),P(10,2),P(10,3),P(10,4),P(10,5)
        LST(21)=12
        RETURN
      ENDIF
      P(IT,1)=0.
      P(IT,2)=0.

      P(NS+7,5)=SQRT(MAX(0.,P(NS+7,4)**2-P(NS+7,1)**2-P(NS+7,2)**2-
     &P(NS+7,3)**2))
      N=NS+7
C     CALL GULIST(103,2)

C...Scale for bremsstrahlung etc.
      Q2PY=Q2
      IPY(40)=10
      IPY(47)=N
C...Save quantities for later use.
      XPY(1)=1.
      XPY(2)=XR
      CALL LYSTFU(K(2,2),XR,Q2,XPQ)
      DO 160 IFL=-6,6
  160 XQPY(2,IFL)=XPQ(IFL)
      IF(LST(23).EQ.1) THEN
        ISUB=39
        IPY(11)=1
      ELSEIF(LST(23).EQ.3) THEN
        ISUB=39
        IPY(11)=2
      ELSEIF(LST(23).EQ.4) THEN
        ISUB=39
        IPY(11)=3
      ELSEIF(LST(23).EQ.2) THEN
        ISUB=40
      ENDIF
      KFL(2,1)=K(5,2)
      KFL(2,2)=K(6,2)
      KFL(1,1)=KFL(2,1)
      KFL(1,2)=KFL(2,2)
      IF(ISUB.EQ.39) KFL(3,1)=K(1,2)
      IF(ISUB.EQ.40) KFL(3,1)=K(1,2)+ISIGN(1,K(1,2))
      KFL(3,2)=K(27,2)
      PYVAR(2)=PARL(21)
      PYVAR(1)=SQRT(PYVAR(2))
      PYVAR(3)=P(1,5)
      PYVAR(4)=P(2,5)
      PYVAR(5)=PYVAR(1)/2.
      IPY(41)=K(1,2)
      IPY(42)=K(2,2)
      IPY(48)=0

C...Generate timelike parton shower (if required)
      IF(IPY(13).EQ.1) THEN
        CALL LSCALE(1,QMAX)
        CALL LUSHOW(25,26,QMAX)
      ENDIF
      IT=25
      IF(N.GE.27) IT=27
      NS=N
C     CALL GULIST(104,2)

C...Generate spacelike parton shower (if required)
      IPU1=0
      IPU2=23
      IF(XPY(2)*(1.+(P(IT,5)**2+PYPAR(22))/P(21,5)**2).GT.0.999) THEN
        LST(21)=13
        RETURN
      ENDIF
      IF(IPY(14).GE.1) THEN
        CALL LYSSPA(IPU1,IPU2)
        IF(LST(21).NE.0) RETURN
      ENDIF
      IF(N.EQ.NS) THEN 
        DO 220 I=NS+1,NS+4
        DO 220 J=1,5
        K(I,J)=0
        P(I,J)=0.
  220   V(I,J)=0.
        K(NS+1,1)=11
        K(NS+1,2)=KFL(2,1)
        K(NS+1,3)=21
        DO 230 J=1,5
  230   P(NS+1,J)=P(21,J)
        K(NS+2,1)=-1
        K(NS+2,3)=NS+1
        K(NS+3,1)=13
        K(NS+3,2)=KFL(2,2)
        K(NS+3,3)=23
        K(NS+3,4)=23
        K(NS+3,5)=23
        P(NS+3,3)=(P(IT,5)**2+Q2)*(P(21,4)-P(21,3))/(2.*Q2)
        P(NS+3,4)=-P(NS+3,3)
        K(NS+4,1)=-1
        K(NS+4,3)=NS+3
        P(NS+4,3)=23
        P(NS+4,4)=23
        P(24,1)=NS+3
        P(24,2)=NS+3
        K(23,4)=K(23,4)+(NS+3)*MSTU(5)
        K(23,5)=K(23,5)+(NS+3)*MSTU(5)
        IPU1=0
        IPU2=NS+3
        N=N+4
      ENDIF
C     CALL GULIST(105,2)

C...Rotate and boost outgoing parton shower
      IF(N.GT.31) THEN
        K(N+1,1)=0
        DO 210 J=1,4
  210   P(N+1,J)=P(NS+1,J)+P(NS+3,J)
        IF(P(N+1,4).LE.1.01*P(IT,5)) THEN
          LST(21)=14
          RETURN
        ENDIF
        ROBO(1)=ULANGL(P(IT,3),SQRT(P(IT,1)**2+P(IT,2)**2))
        ROBO(2)=ULANGL(P(IT,1),P(IT,2))
      IF(ABS(ROBO(1)).GT.0.001.OR.ABS(ROBO(2)).GT.0.001) THEN
      WRITE(6,*) '1:',IT,K(IT,2),P(IT,1),P(IT,2),P(IT,3),P(IT,4),P(IT,5)
      WRITE(6,*) '   ROBO(1-2)=',ROBO(1),ROBO(2)
      ENDIF
        CALL LUDBRB(25,NS,0.,-ROBO(2),0.D0,0.D0,0.D0)
        CALL LUDBRB(25,NS,-ROBO(1),0.,0.D0,0.D0,0.D0)
C...Replace old rotation method with x,y,z-boost to preserve QCD phi-dep
        DELTAP(1)=DBLE(P(N+1,1))
        DELTAP(2)=DBLE(P(N+1,2))
        DELTAP(3)=DBLE(P(N+1,3)) - DBLE(P(IT,3))
        DELTAP(4)=SQRT(DELTAP(1)**2+DELTAP(2)**2+DELTAP(3)**2)
        IF(DELTAP(4).LT.1.D-11) GOTO 410
        DPLONG=-(DBLE(P(IT,3))*DELTAP(3))/DELTAP(4)
        DROOT=MAX(0.D0,DBLE(P(N+1,4))**2-DBLE(P(IT,4))**2+DPLONG**2)
        DBTOT=-(DPLONG*DBLE(P(IT,4))-DBLE(P(N+1,4))*SQRT(DROOT))/
C    &  SQRT(DBLE(P(N+1,4))**2-DBLE(P(IT,4))**2+DPLONG**2))/
     &  (DPLONG**2+DBLE(P(N+1,4))**2)
        DGAMMA=1.D0/SQRT(1.D0-DBTOT**2)
        DO 400 I = 1,3
  400   DROBO(I+2)=DELTAP(I)/(DGAMMA/(DGAMMA+1.D0)*
     &  (DBLE(P(N+1,4))-DGAMMA*DBLE(P(IT,4)))+DGAMMA*DBLE(P(IT,4)))
        CALL LUDBRB(25,NS,0.,0.,DROBO(3),DROBO(4),DROBO(5))
  410   CONTINUE
C...End phi-correction
      ENDIF
C     CALL GULIST(106,2)

      Q2PY=Q2
C...Hadron remnant and primordial kt
      IPY(47)=N
      CALL LYREMN(IPU1,IPU2)
      IF(IPY(48).EQ.1) THEN
        LST(21)=15
        RETURN
      ENDIF
C     CALL GULIST(107,2)

C...Rearrange partons along strings
      MSTU(24)=0
      CALL LUPREP(0)
      IF(MSTU(24).NE.0) THEN
C       CALL GULIST(188,2)
        IF(LST(3).GE.1) WRITE(6,*) ' LUPREP error MSTU(24)= ',MSTU(24)
        LST(21)=16
        RETURN
      ENDIF
C     CALL GULIST(109,2)

C...Clean up event record -> order:
C...1=inc. lepton; 2=inc. nucleon; 3=exch boson; 4=scat. lepton;
C...5=inc. parton before initial shower; 6=inc. parton at hard scattering
C...after shower; 7,8=first,second parton from hard scattering
C...before final shower
      LST(26)=7
      DO 510 J=1,5
      K(N+1,J)=K(4,J)
  510 P(N+1,J)=P(4,J)
      DO 520 J=1,5
      K(3,J)=K(5,J)
      P(3,J)=P(5,J)
      K(4,J)=K(9,J)
      P(4,J)=P(9,J)
      K(5,J)=K(N+1,J)
      P(5,J)=P(N+1,J)
      K(6,J)=K(NS+3,J)
      P(6,J)=P(NS+3,J)
C     K(7,J)=K(IT,J)
C     P(7,J)=P(IT,J)
      K(7,J)=K(25,J)
      P(7,J)=P(25,J)
      K(8,J)=K(26,J)
      P(8,J)=P(26,J)
  520 CONTINUE
      K(3,3)=1
      K(4,3)=1
      K(6,1)=21
      K(6,3)=5
      K(6,4)=0
      K(6,5)=0
      K(7,1)=21
      K(7,3)=6
      K(7,4)=0
      K(7,5)=0
      K(8,1)=21
      K(8,3)=6
      K(8,4)=0
      K(8,5)=0
C...Activate line with scattered lepton.
      K(4,1)=1
C...Deactivate obsolete lines 9, 10, 21, NS+1 (extra lines with boson)
      K(9,1)=0
      K(10,1)=0
      K(21,1)=0
      IF(K(NS+1,2).EQ.K(3,2)) K(NS+1,1)=0
C...Zero irrelevant lines with K(I,1)<0
      DO 540 I=1,N
      IF(K(I,1).LT.0) THEN
        DO 530 J=1,5
        K(I,J)=0
  530   P(I,J)=0.
      ENDIF
  540 CONTINUE
C     CALL GULIST(110,2)
C...Delete internal parton lines, i.e. with K(I,1)=13,14
      IF(MOD(LST(4)/10,10).EQ.0) THEN
        CALL LTIMEX(T1)
        CALL LUEDIT(14)
        CALL LTIMEX(T2)
C       CALL GULIST(111,2)
      ENDIF
C...Delete empty lines
      CALL LTIMEX(T1)
      CALL LUEDIT(12)
      CALL LTIMEX(T2)
C     CALL GULIST(112,2)

      RETURN
      END

C **********************************************************************

      SUBROUTINE LSCALE(INFIN,QMAX)

C...Give maximum virtuality of partons in parton showers.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
C...Power in f(x0)=(1-x0)**power used for scale x0*W2
      DATA POWER/3./
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)

      QMAX=0.1
      IF(LST(8).GE.2.AND.LST(8).LE.5) THEN
C...Parton showers without matrix elements matching
        IF(LST(9).EQ.1) THEN
          QMAX=Q2
        ELSEIF(LST(9).EQ.2) THEN
          QMAX=W2
        ELSEIF(LST(9).EQ.3) THEN
          QMAX=SQRT(W2*Q2)
        ELSEIF(LST(9).EQ.4) THEN
          QMAX=Q2*(1.-X)
        ELSEIF(LST(9).EQ.5) THEN
          QMAX=Q2*(1.-X)*MAX(1.,LOG(1./MAX(1.E-06,X)))
        ELSEIF(LST(9).EQ.6) THEN
          X0=1.D0-(1.D0-DBLE(X))*RLU(0)**(1./(POWER+1.))
          QMAX=X0*W2
        ELSEIF(LST(9).EQ.9) THEN
          QMAX=W2**(2./3.)
        ELSE
          WRITE(6,*) ' Warning, LSCALE: LST(9)=',LST(9),' not allowed'
        ENDIF
      ELSEIF(LST(8).GT.10.AND.LST(24).EQ.1.AND.MOD(LST(8),10).NE.9) THEN
C...Parton showers added to q-event from 1st order matrix elements
        IF(LST(20).LE.1) THEN
          QMAX=PARL(27)*W2
        ELSEIF(LST(20).EQ.2) THEN
          QMAX=PARL(27)*Q2
        ELSEIF(LST(20).EQ.3) THEN
          QMAX=PARL(9)*Q2
        ELSEIF(LST(20).EQ.4) THEN
          QMAX=PARL(9)*Q2
          IF(INFIN.LT.0) QMAX=PARL(27)*Q2/X
        ELSEIF(LST(20).EQ.5) THEN
          QMAX=PARL(9)
          IF(INFIN.LT.0) QMAX=PARL(27)*Q2/X
        ELSEIF(LST(20).EQ.6) THEN
          QMAX=PARL(27)
          IF(INFIN.LT.0) QMAX=PARL(8)*Q2/X
        ELSE
          WRITE(6,*) 'LSCALE: No such jet scheme!'
        ENDIF
      ELSEIF(LST(8).GT.10.AND.MOD(LST(8),10).NE.9) THEN
C...Parton showers added to qg-/qqbar-event from 1st order matrix elements
C...Scale given by invariant mass of final parton pair
        QMAX=P(27,5)**2
        IF(INFIN.LT.0)
     &  QMAX=MAX(ABS(-Q2-2.*FOUR(25,21)),ABS(-Q2-2.*FOUR(26,21)))
      ENDIF
      QMAX=SQRT(QMAX)

      RETURN
      END
C ********************************************************************

      SUBROUTINE LYSSPA(IPU1,IPU2)
C--                                                                   --C
C--   Created:        -                                               --C
C--   Last update:    960801                                          --C

C...NEW X REDEFINITION
C...GENERATES SPACELIKE PARTON SHOWERS
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /LYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON /LYPROC/ ISUB,KFL(3,2),X(2),SH,TH,UH,Q2,XSEC(0:40)
      COMMON /LYINT1/ XQ(2,-6:6),DSIG(-6:6,-6:6,5),FSIG(10,10,3)
      DIMENSION IFLS(4),IS(2),XS(2),ZS(2),Q2S(2),TEVS(2),ROBO(5),
     &XFS(2,-6:6),XFA(-6:6),XFB(-6:6),WTAP(-6:6),WTSF(-6:6)
      DOUBLE PRECISION DQ2(3),DSH,DSHZ,DSHR,DPLCM,DPC(3),DPD(4),DMS,
     &DMSMA,DPT2,DPB(4),DBE1(4),DBE2(4),DBEP,DGABEP,DPQ(4),DPQS(2),
     &DM2,DQ2B,DROBO(5),DBEZ,DTEMP
C-GI &DQ23,DPH(4),DM2,DQ2B,DQM2
CJR--begin
      LOGICAL SLAST
      LOGICAL SEAQUARK,SPLIT
      REAL XFT(-6:6)
      REAL XQUARK,XGLUON,XSEA,ZSPLIT,ZSOFT,ZMAX
CTEST      COMMON /SEAQTE/ XQUARK,XGLUON,XSEA,ZSPLIT,ZSOFT,ZMAX,SPLIT
      INTEGER LASTFL,SEAFL
CTEST      COMMON /FLAVOR/ LASTFL,SEAFL
      REAL LEXSEA
CJR--end
      DATA IFLA,NQ/0,0/,Z,XE0,XA/3*0./,DSHZ,DMSMA,DPT2,DSHR/4*0.D0/

C...COMMON CONSTANTS, SET UP INITIAL VALUES
      ILEP=0
      IF(IPU1.EQ.0) ILEP=1
      IF(IPU2.EQ.0) ILEP=2
      Q2E=Q2
C-GI  IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.3) Q2E=Q2E/PYPAR(26)
      IF(ISUB.EQ.27) Q2E=PMAS(23,1)**2
      IF(ISUB.EQ.28) Q2E=PMAS(24,1)**2
      TMAX=ALOG(PYPAR(26)*PYPAR(27)*Q2E/PYPAR(21)**2)
      IF(ILEP.GE.1) THEN
        SH=P(25,5)**2
        IF(N.GE.27) SH=P(27,5)**2
        CALL LSCALE(-1,QMAX)
        Q2E=QMAX**2
        Q2E=MAX(PYPAR(21)**2,MIN(Q2E,(0.95/X(3-ILEP)-1.)*Q2-SH,
     &  Q2/2.+SH))
        TMAX=ALOG(Q2E/PYPAR(21)**2)
      ENDIF
CJR--begin
CJR-- varying cut-off in virtuality of incoming seaquarks
      IF (LST(35).EQ.2) THEN
         CALL LPRIKT(PARL(10),PT,PHI)
         PYPAR22=MAX(PYPAR(21)**2,PT**2)
      ELSE
         PYPAR22=PYPAR(22)
      ENDIF
CJR--end
      IF (MOD(LST(8),10).EQ.4 .OR. MOD(LST(8),10).EQ.5) THEN
         Q2E=PYPAR(22)
         TMAX=ALOG(Q2E/PYPAR(21)**2)
      ELSEIF(PYPAR(26)*Q2E.LT.MAX(PYPAR(22),2.*PYPAR(21)**2).OR.
     &TMAX.LT.0.2) THEN
         RETURN
      ENDIF
CJR--end
      IF(ILEP.EQ.0) XE0=2.*PYPAR(23)/PYVAR(1)
      B0=(33.-2.*IPY(8))/6.
      NS=N
      MSTU(2)=0
CJR--begin
      NTRY=0
  100 N=NS
      NTRY=NTRY+1
      IF (NTRY.GT.100) THEN
         LST(21)=17
         RETURN
      ENDIF
CJR--end
CJR  100 N=NS
CJR--begin
      SEAQUARK=.FALSE.
      SPLIT=.FALSE.
      SLAST=.FALSE.
CJR--end
      IF(ILEP.GE.1) THEN
        NQ=IPU2-2
        IF(ILEP.EQ.2) NQ=IPU1+2
        DPQS(1)=DBLE(P(NQ,3))
        DPQS(2)=DBLE(P(NQ,4))
        XBMIN=X(3-ILEP)*MAX(0.5,SH/Q2)
        CALL LYSTFU(IPY(43-ILEP),XBMIN,Q2,XFB)
        DO 110 IFL=-6,6
  110   XQ(3-ILEP,IFL)=XFB(IFL)
      ENDIF
      DO 120 JT=1,2
      IFLS(JT)=KFL(2,JT)
      IF(KFL(2,JT).EQ.21) IFLS(JT)=0
      IFLS(JT+2)=IFLS(JT)
      XS(JT)=X(JT)
      ZS(JT)=1.
      IF(ILEP.EQ.0) Q2S(JT)=PYPAR(26)*Q2E
      TEVS(JT)=TMAX
      DO 120 IFL=-6,6
  120 XFS(JT,IFL)=XQ(JT,IFL)
      IF(ILEP.GE.1) THEN
        Q2S(ILEP)=P(NQ,5)**2
        DQ2(ILEP)=Q2S(ILEP)
        Q2S(3-ILEP)=Q2E
      ENDIF
      DSH=SH
      IHFC=0
      IHFX=0

C...PICK UP LEG WITH HIGHEST VIRTUALITY
CJR
      NTURN=0
CJR
  130 CONTINUE
CJR
      NTURN=NTURN+1
CJR
      IF(N.GT.MSTU(4)-10) THEN
        WRITE(6,*) ' LYSSPA: no more memory in LUJETS'
        LST(21)=18
        RETURN
      ENDIF
      DO 133 I=N+1,N+8
      DO 133 J=1,5
      K(I,J)=0
  133 P(I,J)=0.
C     CALL GULIST(21,2)
      N=N+2
      JT=1
      IF((N.GT.NS+2.AND.Q2S(2).GT.Q2S(1).AND.ILEP.EQ.0).OR.ILEP.EQ.1)
     &JT=2
      JR=3-JT
      IFLB=IFLS(JT)
      XB=XS(JT)
      IF(ILEP.GE.1.AND.N.EQ.NS+2) XB=XS(JT)*MAX(SH/Q2,0.5)
      DO 140 IFL=-6,6
  140 XFB(IFL)=XFS(JT,IFL)
      Q2B=Q2S(JT)
      TEVB=TEVS(JT)
      IF(IPY(14).GE.9.AND.N.GT.NS+4) THEN
        Q2B=0.5*(1./ZS(JT)+1.)*Q2S(JT)+0.5*(1./ZS(JT)-1.)*(Q2S(3-JT)-
     &  SNGL(DSH)+SQRT((SNGL(DSH)+Q2S(1)+Q2S(2))**2+8.*Q2S(1)*Q2S(2)*
     &  ZS(JT)/(1.-ZS(JT))))
        TEVB=ALOG(PYPAR(27)*Q2B/PYPAR(21)**2)
      ENDIF
      IF(ILEP.EQ.0) THEN
        DSHR=2.*DSQRT(DSH)
        DSHZ=DSH/DBLE(ZS(JT))
      ELSEIF(ILEP.GE.1) THEN
        DSHZ=DSH
        IF(N.GT.NS+4) DSHZ=(DSH+DQ2(JR)-DQ2(JT))/ZS(JT)-DQ2(JR)+
     &  PYPAR22
        DPD(2)=DSHZ+DQ2(JR)+DBLE(PYPAR22)

        MSTJ(93)=1
        QMASS=ULMASS(IABS(IFLB))
        IF(IABS(IFLB).EQ.0) QMASS=ULMASS(21)
C...CHECK IF QUARK PAIR CREATION ONLY POSSIBILITY
        IF(DQ2(JR).LT.4.*QMASS**2) THEN
          DM2=QMASS**2
          DPC(1)=DQ2(JR)*(DBLE(PYPAR22)+DM2)**2
          DPC(2)=DPD(2)*(DPD(2)-2D0*PYPAR22)*(PYPAR22+DM2)
          DPC(3)=PYPAR22*(DPD(2)-2D0*PYPAR22)**2
          XE0=1D0-(DPC(2)-DSQRT(DPC(2)**2-4D0*DPC(1)*DPC(3)))/
     &    (2D0*DPC(1))
CJR--begin
          ZMAX=(DPC(2)-DSQRT(DPC(2)**2-4D0*DPC(1)*DPC(3)))/(2D0*DPC(1))
          XE0=XB*(1./ZMAX-1.)
CJR--end
        ELSE
          XE0=1D0-(DPD(2)-2D0*DBLE(PYPAR22))*(DPD(2)-DSQRT(DPD(2)**2-
     &    4D0*DQ2(JR)*DBLE(PYPAR22)))/(2D0*DQ2(JR)*DBLE(PYPAR22))
CJR--begin
          ZMAX=(DPD(2)-2D0*DBLE(PYPAR22))*(DPD(2)-DSQRT(DPD(2)**2-
     &    4D0*DQ2(JR)*DBLE(PYPAR22)))/(2D0*DQ2(JR)*DBLE(PYPAR22))
          XE0=XB*(1./ZMAX-1.)
CJR--end
        ENDIF
CJR--begin
CJR--          radiated parton energy cut
C        XE0=MAX(XE0,2.*PYPAR(23)/SQRT(W2LP))
CJR--end
      ENDIF
CJR  145 XE=MAX(XE0,XB*(1./(1.-PYPAR(24))-1.))
CJR--begin
      N145=0
145   CONTINUE
      N145=N145+1
      IF (N145.GT.100) THEN
CJR         WRITE(*,*) '145'
         GOTO 100
      ENDIF
      XE=MAX(XE0,XB*(1./(1.-PYPAR(24))-1.))
        ZMAX=XB/(XB+XE)
CJR--end
      IF(XB+XE.GE.0.999) THEN
        Q2B=0.
        GOTO 210
      ENDIF

C...CALCULATE ALTARELLI-PARISI AND STRUCTURE FUNCTION WEIGHTS
      DO 150 IFL=-6,6
      WTAP(IFL)=0.
  150 WTSF(IFL)=0.
      IF(IFLB.EQ.0) THEN
        WTAPQ=16.*(1.-SQRT(XB+XE))/(3.*SQRT(XB))
        DO 160 IFL=-IPY(8),IPY(8)
        IF(IFL.EQ.0) WTAP(IFL)=6.*ALOG((1.-XB)/XE)
  160   IF(IFL.NE.0) WTAP(IFL)=WTAPQ
      ELSE
        WTAP(0)=0.5*XB*(1./(XB+XE)-1.)
        WTAP(IFLB)=8.*ALOG((1.-XB)*(XB+XE)/XE)/3.
      ENDIF
  170 WTSUM=0.
      IF(IHFC.EQ.0) THEN
        DO 180 IFL=-IPY(8),IPY(8)
        WTSF(IFL)=XFB(IFL)/MAX(1E-10,XFB(IFLB))
  180   WTSUM=WTSUM+WTAP(IFL)*WTSF(IFL)
        IF(IABS(IFLB).GE.4.AND.WTSUM.GT.1E3) THEN
          IHFX=1
          DO 185 IFL=-IPY(8),IPY(8)
  185     WTSF(IFL)=WTSF(IFL)*1E3/WTSUM
          WTSUM=1E3
        ENDIF
      ENDIF

C...CHOOSE NEW T AND FLAVOUR
CJR  190 IF(IPY(14).LE.6.OR.IPY(14).GE.9) THEN
CJR--begin
      NJR=0
190   CONTINUE
      SEAQUARK=.FALSE.
      NJR=NJR+1
      IF (NJR.GT.100) THEN
CJR         WRITE(*,*) '190'
         GOTO 100
      ENDIF
CJR--end
      IF(IPY(14).LE.6.OR.IPY(14).GE.9) THEN
        TEVXP=B0/MAX(0.0001,WTSUM)
      ELSE
        TEVXP=B0/MAX(0.0001,5.*WTSUM)
      ENDIF
      TEVB=TEVB*EXP(MAX(-100.,ALOG(RLU(0))*TEVXP))
      Q2REF=PYPAR(21)**2*EXP(TEVB)/PYPAR(27)
      Q2B=Q2REF/PYPAR(27)
CJR--begin -- 
      Q2BOLD=Q2B
CAE--seaquarks up to LST(12), if the quark density != 0
      IF( LST(35).EQ.1 .OR. LST(35).EQ.2) THEN
         IF( Q2B.LT.PYPAR(22) .AND.
     &       (ABS(IFLB).LE.LST(12).AND.ABS(IFLB).GE.1)) THEN
            Q2REF=MIN(PYPAR(22),Q2S(JT))
            Q2B=MIN(PYPAR22,Q2S(JT))
            IF(ILEP.GE.1.AND.N.EQ.NS+2) THEN
               XT=X(JT)*(1.+(DSH-Q2B)/DQ2(JR))
            ELSE
               XT=XB
            ENDIF
            CALL LYSTFU(IPY(40+JT),XT,Q2REF,XFT)
            IF(XFT(IFLB).EQ.0.0.AND.XFT(-IFLB).EQ.0.0) THEN
              SEARATIO=1.0
            ELSEIF(XFT(ABS(IFLB)).EQ.0.0) THEN
              SEARATIO=1.0
            ELSE
              SEARATIO=XFT(-ABS(IFLB))/XFT(ABS(IFLB))
            ENDIF
CJR-- (protons only)
            IF (RLU(0).LT.SEARATIO) THEN
               SEAQUARK=.TRUE.
               XQUARK=XT
            ELSE
               Q2B=Q2BOLD
               SEAQUARK=.FALSE.
            ENDIF
         ENDIF
      ENDIF
CJR--end
      DQ2B=Q2B
      IF(ILEP.GE.1) THEN
        DSHZ=DSH
        IF(N.GT.NS+4) DSHZ=(DSH+DQ2(JR)-DQ2(JT))/DBLE(ZS(JT))-DQ2(JR)+
     &  DQ2B
      ENDIF
CJR      IF(Q2B.LT.PYPAR(22)) THEN
      IF(Q2B.LT.PYPAR(22).AND.(.NOT.SEAQUARK)) THEN
         Q2B=0.
      ELSE
        WTRAN=RLU(0)*WTSUM
        IFLA=-IPY(8)-1
  200   IFLA=IFLA+1
        WTRAN=WTRAN-WTAP(IFLA)*WTSF(IFLA)
        IF(IFLA.LT.IPY(8).AND.WTRAN.GT.0.) GOTO 200

CJR--begin
        IF (SEAQUARK) THEN
           SEAFL=-IFLB
           IFLA=0
CT           XE=XB*(1./(1.-0.001)-1.)
        ENDIF
CJR--end

C...CHOOSE Z VALUE AND CORRECTIVE WEIGHT
        IF(IFLB.EQ.0.AND.IFLA.EQ.0) THEN
          Z=1./(1.+((1.-XB)/XB)*(XE/(1.-XB))**RLU(0))
          WTZ=(1.-Z*(1.-Z))**2
        ELSEIF(IFLB.EQ.0) THEN
          Z=XB/(1.-RLU(0)*(1.-SQRT(XB+XE)))**2
          WTZ=0.5*(1.+(1.-Z)**2)*SQRT(Z)
        ELSEIF(IFLA.EQ.0) THEN
          Z=XB*(1.+RLU(0)*(1./(XB+XE)-1.))
          WTZ=1.-2.*Z*(1.-Z)
        ELSE
          Z=1.-(1.-XB)*(XE/((XB+XE)*(1.-XB)))**RLU(0)
          WTZ=0.5*(1.+Z**2)
        ENDIF
CJR--begin
C        IF (SEAQUARK) THEN
C           XSEA=LEXSEA(0.15*XT,Q2B)
C           XE=MIN(XE,XSEA)
C           Z=XT/(XSEA+XT)
C        ENDIF
CJR--end
C...REWEIGHT FIRST LEG BECAUSE OF MODIFIED XB OR CHECK PHASE SPACE
        IF(ILEP.GE.1.AND.N.EQ.NS+2) THEN
          XBNEW=X(JT)*(1.+(DSH-Q2B)/DQ2(JR))
          IF(XBNEW.GT.MIN(Z,0.999)) GOTO 190
          XB=XBNEW
        ENDIF
C...SUM UP SOFT GLUON EMISSION AS EFFECTIVE Z SHIFT
CJR--       should this realy always be done ??
        IF(IPY(15).GE.1) THEN
          RSOFT=6.
          IF(IFLB.NE.0) RSOFT=8./3.
          Z=Z*(TEVB/TEVS(JT))**(RSOFT*XE/((XB+XE)*B0))
          IF(Z.LE.XB) GOTO 190
CJR--begin
          ZSOFT=(TEVB/TEVS(JT))**(RSOFT*XE/((XB+XE)*B0))
          ZMAX=XB/(XB+XE)
CJR--end
        ENDIF
C...CHECK IF HEAVY FLAVOUR BELOW THRESHOLD
        IHFT=0
CIC...Skip for intrinsic charm/bottom simulation, charm quark should
CIC...not come from gluon but is non-perturbative part of proton.
        IF(LST(15).EQ.-4.OR.LST(15).EQ.-5) GOTO 205
        MSTJ(93)=1
        IF(ILEP.GE.1.AND.IABS(IFLB).GE.4.AND.(XFB(IFLB).LT.1E-10.OR.
     &    Q2B.LT.5.*ULMASS(IABS(IFLB))**2)) THEN
          IHFT=1
          IFLA=0
        ENDIF
  205   CONTINUE

C...FOR LEPTOPRODUCTION, CHECK Z AGAINST NEW LIMIT
        IF(ILEP.GE.1) THEN
          DPD(2)=DSHZ+DQ2(JR)+DQ2B
          MSTJ(93)=1
          DM2=ULMASS(IABS(IFLA-IFLB))**2
          IF(IABS(IFLA-IFLB).EQ.0) DM2=ULMASS(21)**2
          DPC(1)=DQ2(JR)*(DQ2B+DM2)**2
          DPC(2)=DPD(2)*(DPD(2)-2D0*DQ2B)*(DQ2B+DM2)
          DPC(3)=DQ2B*(DPD(2)-2D0*DQ2B)**2
          ZU=(DPC(2)-DSQRT(DPC(2)**2-4D0*DPC(1)*DPC(3)))/(2D0*DPC(1))
          IF(Z.GE.ZU) GOTO 190
        ENDIF

C...OPTION WITH EVOLUTION IN KT2=(1-Z)Q2:
        IF(IPY(14).GE.5.AND.IPY(14).LE.6.AND.N.LE.NS+4) THEN
C...CHECK THAT (Q2)LAST BRANCHING < (Q2)HARD
          IF(Q2B/(1.-Z).GT.PYPAR(26)*Q2) GOTO 190
        ELSEIF(IPY(14).GE.3.AND.IPY(14).LE.6.AND.N.GE.NS+6) THEN
C...CHECK THAT Z,Q2 COMBINATION IS KINEMATICALLY ALLOWED
          Q2MAX=0.5*(1./ZS(JT)+1.)*DQ2(JT)+0.5*(1./ZS(JT)-1.)*
     &    (DQ2(3-JT)-DSH+SQRT((DSH+DQ2(1)+DQ2(2))**2+8.*DQ2(1)*DQ2(2)*
     &    ZS(JT)/(1.-ZS(JT))))
          IF(Q2B/(1.-Z).GE.Q2MAX) GOTO 190

        ELSEIF(IPY(14).EQ.7.OR.IPY(14).EQ.8) THEN
C...OPTION WITH ALPHAS((1-Z)Q2): DEMAND KT2 > CUTOFF, REWEIGHT
          IF((1.-Z)*Q2B.LT.PYPAR22) GOTO 190
          ALPRAT=TEVB/(TEVB+ALOG(1.-Z))
          IF(ALPRAT.LT.5.*RLU(0)) GOTO 190
          IF(ALPRAT.GT.5.) WTZ=WTZ*ALPRAT/5.
        ENDIF

C...WEIGHTING WITH NEW STRUCTURE FUNCTIONS
        CALL LYSTFU(IPY(40+JT),XB,Q2REF,XFB)
        XA=XB/Z
        CALL LYSTFU(IPY(40+JT),XA,Q2REF,XFA)
        IF(IHFT.EQ.1.OR.IHFX.EQ.1) THEN
           IF(XFA(IFLA).LT.1E-10) IHFC=1
           GOTO 210
        ELSEIF(XFB(IFLB).LT.1E-20) THEN
          GOTO 100
        ENDIF
        IF(WTZ*XFA(IFLA)/XFB(IFLB).LT.RLU(0)*WTSF(IFLA)) THEN
          IF(ILEP.GE.1.AND.N.EQ.NS+2) GOTO 145
          GOTO 170
        ENDIF
CJR--begin
        IF (SEAQUARK) THEN
           SPLIT=.TRUE.
           XGLUON=XA
           XSEA=XA-XB
           ZSPLIT=Z
           SEAQUARK=.FALSE.
        ENDIF
CJR--end
      ENDIF

210   CONTINUE
      IF(N.EQ.NS+4-2*MIN(1,ILEP)) THEN
C...DEFINE TWO HARD SCATTERERS IN THEIR CM-FRAME
        DQ2(JT)=Q2B
        IF(IPY(14).GE.3.AND.IPY(14).LE.6) DQ2(JT)=Q2B/(1.-Z)
        IF(ILEP.EQ.0) THEN
          DPLCM=DSQRT((DSH+DQ2(1)+DQ2(2))**2-4.*DQ2(1)*DQ2(2))/DSHR
          DO 220 JR=1,2
          I=NS+2*JR-1
          IPO=19+2*JR
          K(I,1)=14
          K(I,2)=IFLS(JR+2)
          IF(IFLS(JR+2).EQ.0) K(I,2)=21
          K(I,3)=0
          K(I,4)=IPO
          K(I,5)=IPO
          P(I,1)=0.
          P(I,2)=0.
          P(I,3)=DPLCM*(-1)**(JR+1)
          P(I,4)=(DSH+DQ2(3-JR)-DQ2(JR))/DSHR
          P(I,5)=-SQRT(SNGL(DQ2(JR)))
          K(I+1,1)=-1
          K(I+1,2)=K(IPO+1,2)
          K(I+1,3)=I
          K(I+1,4)=0
          K(I+1,5)=0
          P(I+1,1)=0.
          P(I+1,2)=0.
          P(I+1,3)=IPO
          P(I+1,4)=IPO
          P(I+1,5)=0.
          P(IPO+1,1)=I
          P(IPO+1,2)=I
          K(IPO,4)=MOD(K(IPO,4),MSTU(5))+I*MSTU(5)
          K(IPO,5)=MOD(K(IPO,5),MSTU(5))+I*MSTU(5)
  220     CONTINUE
        ELSE
C..LEPTOPRODUCTION EVENTS: BOSON AND HADRON REST FRAME
          I1=NS+2*ILEP-1
          I2=NS-2*ILEP+5
          DO 225 ITEMP=NS+1,NS+4
          DO 225 J=1,5
          K(ITEMP,J)=0
  225     P(ITEMP,J)=0.
          DO 230 J=1,5
  230     P(I1,J)=P(NQ,J)
          K(NS+1,1)=11
          K(NS+3,1)=14
          IF(ILEP.EQ.2) THEN
            K(NS+1,1)=14
            K(NS+3,1)=11
          ENDIF
          K(NS+2,1)=-1
          K(NS+4,1)=-1
          K(NS+1,3)=0
          K(NS+2,3)=NS+1
          K(NS+3,3)=0
          K(NS+4,3)=NS+3
          K(I1,2)=KFL(2,ILEP)
          K(I2,2)=KFL(2,3-ILEP)
          DPD(1)=DSH+DQ2(1)+DQ2(2)
          DPD(3)=(3-2*ILEP)*DSQRT(DPD(1)**2-4D0*DQ2(1)*DQ2(2))
          P(I2,3)=(DPQS(2)*DPD(3)-DPQS(1)*DPD(1))/
     &    (2D0*DQ2(JR))
          P(I2,4)=(DPQS(1)*DPD(3)-DPQS(2)*DPD(1))/
     &    (2D0*DQ2(JR))
          P(I2,5)=-SQRT(SNGL(DQ2(3-ILEP)))
          P(I2+1,3)=MAX(IPU1,IPU2)
          P(I2+1,4)=MAX(IPU1,IPU2)
          K(I2,4)=K(I2,4)-MOD(K(I2,4),MSTU(5))+MAX(IPU1,IPU2)
          K(I2,5)=K(I2,5)-MOD(K(I2,5),MSTU(5))+MAX(IPU1,IPU2)
          P(26-2*ILEP,1)=I2
          P(26-2*ILEP,2)=I2
          K(25-2*ILEP,4)=MOD(K(25-2*ILEP,4),MSTU(5))+I2*MSTU(5)
          K(25-2*ILEP,5)=MOD(K(25-2*ILEP,5),MSTU(5))+I2*MSTU(5)
          N=N+2
        ENDIF

      ELSEIF(N.GT.NS+4) THEN
C...FIND MAXIMUM ALLOWED MASS OF TIMELIKE PARTON
        DQ2(3)=Q2B
        IF(IPY(14).GE.3.AND.IPY(14).LE.6) DQ2(3)=Q2B/(1.-Z)
        IF(IS(1).GE.1.AND.IS(1).LE.MSTU(4)) THEN
          DPC(1)=P(IS(1),4)
          DPC(3)=0.5*(ABS(P(IS(1),3))+ABS(P(IS(2),3)))
        ELSE
C...IS(1) not initialized
          DPC(1)=0.
          DPC(3)=0.5*(       0.      +ABS(P(IS(2),3)))
        ENDIF
        DPC(2)=P(IS(2),4)
        DPD(1)=DSH+DQ2(JR)+DQ2(JT)
        DPD(2)=DSHZ+DQ2(JR)+DQ2(3)
        DPD(3)=DSQRT(DPD(1)**2-4.*DQ2(JR)*DQ2(JT))
        DPD(4)=DSQRT(DPD(2)**2-4.*DQ2(JR)*DQ2(3))
        IKIN=0
        IF((Q2S(JR).GE.0.5*PYPAR22.AND.DPD(1)-DPD(3).GE.1D-10*DPD(1))
     &  .OR.ILEP.GE.1) IKIN=1
        IF(IKIN.EQ.0) DMSMA=(DQ2(JT)/DBLE(ZS(JT))-DQ2(3))*(DSH/
     &  (DSH+DQ2(JT))-DSH/(DSHZ+DQ2(3)))
        IF(IKIN.EQ.1) DMSMA=(DPD(1)*DPD(2)-DPD(3)*DPD(4))/(2.*
     &  DQ2(JR))-DQ2(JT)-DQ2(3)

C...GENERATE TIMELIKE PARTON SHOWER (IF REQUIRED)
        IT=N-1
        K(IT,1)=3
        K(IT,2)=IFLB-IFLS(JT+2)
        IF(IFLB-IFLS(JT+2).EQ.0) K(IT,2)=21
        MSTJ(93)=1
        P(IT,5)=ULMASS(K(IT,2))
        IF (SLAST) P(IT,5)=P(IT,5)-PARL(20)
        IF(SNGL(DMSMA).LE.P(IT,5)**2) GOTO 100
        P(IT,2)=0.
        DO 240 J=1,5
        K(IT+1,J)=0
  240   P(IT+1,J)=0.
        K(IT+1,1)=-1
        K(IT+1,2)=K(IS(JT)+1,2)
        K(IT+1,3)=IT
        IF(MOD(IPY(14),2).EQ.0) THEN
          P(IT,1)=0.
          IF(ILEP.EQ.0) P(IT,4)=(DSHZ-DSH-P(IT,5)**2)/DSHR
          IF(ILEP.GE.1) P(IT,4)=0.5*(P(IS(JT),3)*DPD(2)+
     &    DPQS(1)*(DQ2(JT)+DQ2(3)+P(IT,5)**2))/(P(IS(JT),3)*DPQS(2)-
     &    P(IS(JT),4)*DPQS(1))-DPC(JT)
          P(IT,3)=SQRT(MAX(0.,P(IT,4)**2-P(IT,5)**2))
CJR--begin
          IF (SLAST) THEN
             CALL LUSHOW(IT,0,P(IT,5))
          ELSE
             CALL LUSHOW(IT,0,SQRT(MIN(SNGL(DMSMA),PYPAR(25)*Q2)))
          ENDIF
CJR--end
          IF(N.GE.IT+2) P(IT,5)=P(IT+2,5)
          IF(N.GT.MSTU(4)-10) THEN
            WRITE(6,*) ' LYSSPA: no more memory in LUJETS'
            LST(21)=19
            RETURN
          ENDIF
          DO 243 I=N+1,N+8
          DO 243 J=1,5
          K(I,J)=0
  243     P(I,J)=0.
        ENDIF

C...RECONSTRUCT KINEMATICS OF BRANCHING: TIMELIKE PARTON SHOWER
        DMS=P(IT,5)**2
        IF(IKIN.EQ.0.AND.ILEP.EQ.0) DPT2=(DMSMA-DMS)*(DSHZ+DQ2(3))/
     &  (DSH+DQ2(JT))
        IF(IKIN.EQ.1.AND.ILEP.EQ.0) DPT2=(DMSMA-DMS)*(0.5*DPD(1)*
     &  DPD(2)+0.5*DPD(3)*DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)+DMS))/
     &  (4.*DSH*DPC(3)**2)
        IF(IKIN.EQ.1.AND.ILEP.GE.1) DPT2=(DMSMA-DMS)*(0.5*DPD(1)*
     &  DPD(2)+0.5*DPD(3)*DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)+DMS))/
     &  DPD(3)**2
        IF(DPT2.LT.0.) GOTO 100
        K(IT,3)=N+1
        IF (SLAST) K(IT,3)=2
        P(IT,1)=SQRT(SNGL(DPT2))
        IF(ILEP.EQ.0) THEN
          DPB(1)=(0.5*DPD(2)-DPC(JR)*(DSHZ+DQ2(JR)-DQ2(JT)-DMS)/
     &    DSHR)/DPC(3)-DPC(3)
          P(IT,3)=DPB(1)*(-1)**(JT+1)
          P(IT,4)=(DSHZ-DSH-DMS)/DSHR
        ELSE
          DPC(3)=DQ2(JT)+DQ2(3)+DMS
          DPB(2)=DPQS(2)*DBLE(P(IS(JT),3))-DPQS(1)*DPC(JT)
          DPB(1)=0.5D0*(DPC(JT)*DPD(2)+DPQS(2)*DPC(3))/DPB(2)-
     &    DBLE(P(IS(JT),3))
          P(IT,3)=DPB(1)
          P(IT,4)=0.5D0*(DBLE(P(IS(JT),3))*DPD(2)+
     &    DPQS(1)*DPC(3))/DPB(2)-DPC(JT)
        ENDIF
        IF(N.GE.IT+2) THEN
          MSTU(1)=IT+2
          DPB(1)=DSQRT(DPB(1)**2+DPT2)
          DPB(2)=DSQRT(DPB(1)**2+DMS)
          DPB(3)=P(IT+2,3)
          DPB(4)=DSQRT(DPB(3)**2+DMS)
          DBEZ=(DPB(4)*DPB(1)-DPB(3)*DPB(2))/(DPB(4)*DPB(2)-DPB(3)*
     &    DPB(1))
          CALL LUDBRB(MSTU(1),MSTU(2),0.,0.,0.D0,0.D0,DBEZ)
          THE=ULANGL(P(IT,3),P(IT,1))
          CALL LUDBRB(MSTU(1),MSTU(2),THE,0.,0.D0,0.D0,0.D0)
        ENDIF

C...RECONSTRUCT KINEMATICS OF BRANCHING: SPACELIKE PARTON
        K(N+1,1)=14
        K(N+1,2)=IFLB
        IF(IFLB.EQ.0) K(N+1,2)=21
        K(N+1,3)=0
CJR--begin 
CJR--       give all radiated partons 5 as mother particle 
        K(N+1,3)=5
CJR--end
        P(N+1,1)=P(IT,1)
        P(N+1,2)=0.
        P(N+1,3)=P(IT,3)+P(IS(JT),3)
        P(N+1,4)=P(IT,4)+P(IS(JT),4)
        P(N+1,5)=-SQRT(SNGL(DQ2(3)))
        DO 250 J=1,5
        K(N+2,J)=0
  250   P(N+2,J)=0.
        K(N+2,1)=-1
        K(N+2,2)=K(IS(JT)+1,2)
        K(N+2,3)=N+1

C...DEFINE COLOUR FLOW OF BRANCHING
        K(IS(JT),1)=14
        K(IS(JT),3)=N+1
        ID1=IT
      KN1=ISIGN(500+IABS(K(N+1,2)),2*K(N+1,2)+1)
      KD1=ISIGN(500+IABS(K(ID1,2)),2*K(ID1,2)+1)
      IF(K(N+1,2).EQ.21) KN1=500
      IF(K(ID1,2).EQ.21) KD1=500
        IF((KN1.GE.501.AND.KD1.GE.501).OR.(KN1.LT.0.AND.
     &  KD1.EQ.500).OR.(KN1.EQ.500.AND.KD1.EQ.500.AND.
     &  RLU(0).GT.0.5).OR.(KN1.EQ.500.AND.KD1.LT.0))
     &  ID1=IS(JT)
        ID2=IT+IS(JT)-ID1
        P(N+2,3)=ID1
        P(N+2,4)=ID2
        P(ID1+1,1)=N+1
        P(ID1+1,2)=ID2
        P(ID2+1,1)=ID1
        P(ID2+1,2)=N+1
        K(N+1,4)=K(N+1,4)-MOD(K(N+1,4),MSTU(5))+ID1
        K(N+1,5)=K(N+1,5)-MOD(K(N+1,5),MSTU(5))+ID2
        K(ID1,4)=MOD(K(ID1,4),MSTU(5))+(N+1)*MSTU(5)
        K(ID1,5)=MOD(K(ID1,5),MSTU(5))+ID2*MSTU(5)
        K(ID2,4)=MOD(K(ID2,4),MSTU(5))+ID1*MSTU(5)
        K(ID2,5)=MOD(K(ID2,5),MSTU(5))+(N+1)*MSTU(5)
        N=N+2
C     CALL GULIST(22,2)

C...BOOST TO NEW CM-FRAME
        MSTU(1)=NS+1
        IF(ILEP.EQ.0) THEN
          CALL LUDBRB(MSTU(1),MSTU(2),0.,0.,
     &    -DBLE(P(N-1,1)+P(IS(JR),1))/DBLE(P(N-1,4)+P(IS(JR),4)),
     &    0.D0,-DBLE(P(N-1,3)+P(IS(JR),3))/DBLE(P(N-1,4)+P(IS(JR),4)))
          IR=N-1+(JT-1)*(IS(1)-N+1)
          CALL LUDBRB(MSTU(1),MSTU(2),
     &    -ULANGL(P(IR,3),P(IR,1)),PARU(2)*RLU(0),0.D0,0.D0,0.D0)
        ELSE
C...REORIENTATE EVENT WITHOUT CHANGING THE BOSON FOUR MOMENTUM
          DO 260 J=1,4
  260     DPQ(J)=P(NQ,J)
          DBE1(4)=DPQ(4)+DBLE(P(N-1,4))
          DO 270 J=1,3,2
  270     DBE1(J)=-(DPQ(J)+DBLE(P(N-1,J)))/DBE1(4)
          DTEMP=1.D0-DBE1(1)**2-DBE1(3)**2
          IF(DTEMP.LE.0.D0) THEN
            LST(21)=20
            IF(LST(3).GE.1) WRITE(6,*) ' Warning from LYSSPA: sqrt of',
     &      DTEMP,' New event generated.'
            RETURN
          ENDIF
          DBE1(4)=1.D0/DSQRT(DTEMP)
          DBEP=DBE1(1)*DPQ(1)+DBE1(3)*DPQ(3)
          DGABEP=DBE1(4)*(DBE1(4)*DBEP/(1D0+DBE1(4))+DPQ(4))
          DO 280 J=1,3,2
  280     DPQ(J)=DPQ(J)+DGABEP*DBE1(J)
          DPQ(4)=DBE1(4)*(DPQ(4)+DBEP)
          DPC(1)=DSQRT(DPQ(1)**2+DPQ(3)**2)
          DBE2(4)=-(DPQ(4)*DPC(1)-DPQS(2)*DSQRT(DPQS(2)**2+DPC(1)**2-
     &    DPQ(4)**2))/(DPC(1)**2+DPQS(2)**2)
          THE=ULANGL(SNGL(DPQ(3)),SNGL(DPQ(1)))
          DBE2(1)=DBE2(4)*DSIN(DBLE(THE))
          DBE2(3)=DBE2(4)*DCOS(DBLE(THE))
          DBE2(4)=1D0/(1D0-DBE2(1)**2-DBE2(3)*2)

C...CONSTRUCT THE COMBINED BOOST
          DPB(1)=DBE1(4)**2*DBE2(4)/(1D0+DBE1(4))
          DPB(2)=DBE1(1)*DBE2(1)+DBE1(3)*DBE2(3)
          DPB(3)=DBE1(4)*DBE2(4)*(1D0+DPB(2))
          DO 290 JB=1,3,2
  290     DROBO(JB+2)=(DBE1(4)*DBE2(4)*DBE1(JB)+DBE2(4)*DBE2(JB)+
     &    DPB(1)*DBE1(JB)*DPB(2))/DPB(3)
          CALL LUDBRB(MSTU(1),MSTU(2),0.,0.,DROBO(3),0.D0,DROBO(5))
          IF(ILEP.EQ.1) THE=ULANGL(P(NS+1,3),P(NS+1,1))
          IF(ILEP.EQ.2) THE=PARU(1)+ULANGL(P(NS+3,3),P(NS+3,1))
          CALL LUDBRB(MSTU(1),MSTU(2),-THE,PARU(2)*RLU(0),0D0,0D0,0D0)
        ENDIF
        MSTU(1)=0
      ENDIF

C...SAVE QUANTITIES, LOOP BACK
      IS(JT)=N-1
      IF(ILEP.EQ.2.AND.N.EQ.NS+4) IS(JT)=N-3
      Q2S(JT)=Q2B
      DQ2(JT)=Q2B
      IF(IPY(14).GE.3.AND.IPY(14).LE.6) DQ2(JT)=Q2B/(1.-Z)
      DSH=DSHZ
      IF(Q2B.GE.0.5*PYPAR22) THEN
        IFLS(JT+2)=IFLS(JT)
        IFLS(JT)=IFLA
        XS(JT)=XA
        ZS(JT)=Z
        DO 300 IFL=-6,6
  300   XFS(JT,IFL)=XFA(IFL)
        TEVS(JT)=TEVB
      ELSE
        IF(JT.EQ.1) IPU1=N-1
        IF(JT.EQ.2) IPU2=N-1
      ENDIF
      IF((MAX(IABS(1-ILEP)*Q2S(1),MIN(1,2-ILEP)*Q2S(2))
     &.GE.0.5*PYPAR22.OR.N.LE.NS+2) .AND. SPLIT) SLAST=.TRUE.
      IF(MAX(IABS(1-ILEP)*Q2S(1),MIN(1,2-ILEP)*Q2S(2)).GE.0.5*PYPAR22
     &.OR.N.LE.NS+2) GOTO 130

      IF(ILEP.EQ.0) THEN
C...BOOST HARD SCATTERING PARTONS TO FRAME OF SHOWER INITIATORS
        DO 310 J=1,3
  310   DROBO(J+2)=(P(NS+1,J)+P(NS+3,J))/(P(NS+1,4)+P(NS+3,4))
        DO 320 J=1,5
  320   P(N+2,J)=P(NS+1,J)
        MSTU(1)=N+2
        MSTU(2)=N+2
        CALL LUDBRB(N+2,N+2,0.,0.,-DROBO(3),-DROBO(4),-DROBO(5))
        ROBO(2)=ULANGL(P(N+2,1),P(N+2,2))
        ROBO(1)=ULANGL(P(N+2,3),SQRT(P(N+2,1)**2+P(N+2,2)**2))
        MSTU(1)=4
        MSTU(2)=NS
        CALL LUDBRB(4,NS,ROBO(1),ROBO(2),DROBO(3),DROBO(4),DROBO(5))
        MSTU(1)=0
        MSTU(2)=0
      ENDIF

C...STORE USER INFORMATION
      K(21,1)=14
      IF(ILEP.NE.0) K(21,1)=11
      K(23,1)=14
      K(21,3)=NS+1
      K(23,3)=NS+3
      DO 330 JT=1,2
      KFL(1,JT)=IFLS(JT)
      IF(IFLS(JT).EQ.0) KFL(1,JT)=21
  330 PYVAR(30+JT)=XS(JT)

      DO 340 I=NS+1,N
      DO 340 J=1,5
  340 V(I,J)=0.
  
CJR--begin
      LASTFL=IFLA
CJR--end

      RETURN
      END

C **********************************************************************

      SUBROUTINE LYREMN(IPU1,IPU2)

C...ADDS ON TARGET REMNANTS (ONE OR TWO FROM EACH SIDE) AND
C...INCLUDES PRIMORDIAL KT.
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON /LYPROC/ ISUB,KFL(3,2),X(2),SH,TH,UH,Q2,XSEC(0:40)
      DIMENSION KFLCH(2),KFLSP(2),CHI(2),PMS(6),IS(2),ROBO(5)
      DOUBLE PRECISION DBETAX,DBETAZ,DROBO(5)
      DATA IPU,IQ/0,0/,PEI,PE,PZI,PZ,SHS,PZH,PEH/7*0./

C...FIND EVENT TYPE, SET POINTERS
      IF(IPU1.EQ.0.AND.IPU2.EQ.0) RETURN
      ILEP=0
      IF(IPU1.EQ.0) ILEP=1
      IF(IPU2.EQ.0) ILEP=2
      IF(ISUB.EQ.7) ILEP=-1
      IF(ILEP.EQ.1) IQ=21
      IF(ILEP.EQ.2) IQ=23
      IP=MAX(IPU1,IPU2)
      NS=N

C...DEFINE INITIAL PARTONS, INCLUDING PRIMORDIAL KT
  100 DO 120 I=3,4
      IF(I.EQ.3) IPU=IPU1
      IF(I.EQ.4) IPU=IPU2
      K(I,1)=21
      K(I,3)=I-2
      DO 110 J=1,5
  110 P(I,J)=0.
      IF(ISUB.EQ.7) THEN
        K(I,2)=21
        SHS=0.
      ELSEIF(IPU.NE.0) THEN
        K(I,2)=K(IPU,2)
        P(I,5)=P(IPU,5)
        CALL LPRIKT(PARL(3),PTSPL,PHISPL)
        P(I,1)=PTSPL*COS(PHISPL)
        P(I,2)=PTSPL*SIN(PHISPL)
        PMS(I-2)=P(I,5)**2+P(I,1)**2+P(I,2)**2
      ELSE
        K(I,2)=K(IQ,2)
        P(I,5)=-SQRT(Q2)
        PMS(I-2)=-Q2
        SHS=(1.-X(5-I))*Q2/X(5-I)+PYVAR(7-I)**2
      ENDIF
  120 CONTINUE

C...KINEMATICS CONSTRUCTION FOR INITIAL PARTONS
      IF(ILEP.EQ.0) SHS=PYVAR(31)*PYVAR(32)*PYVAR(2)+
     &(P(3,1)+P(4,1))**2+(P(3,2)+P(4,2))**2
      SHR=SQRT(MAX(0.,SHS))
      IF(ILEP.EQ.0) THEN
        IF((SHS-PMS(1)-PMS(2))**2-4.*PMS(1)*PMS(2).LE.0.) GOTO 100
        P(3,4)=0.5*(SHR+(PMS(1)-PMS(2))/SHR)
        P(3,3)=SQRT(MAX(0.,P(3,4)**2-PMS(1)))
        P(4,4)=SHR-P(3,4)
        P(4,3)=-P(3,3)
      ELSEIF(ILEP.EQ.1) THEN
        P(3,4)=P(IQ,4)
        P(3,3)=P(IQ,3)
        P(4,4)=P(IP,4)
        P(4,3)=P(IP,3)
      ELSEIF(ILEP.EQ.2) THEN
        P(3,4)=P(IP,4)
        P(3,3)=P(IP,3)
        P(4,4)=P(IQ,4)
        P(4,3)=P(IQ,3)
      ENDIF

C...TRANSFORM PARTONS TO OVERALL CM-FRAME (NOT FOR LEPTOPRODUCTION)
      IF(ILEP.EQ.0) THEN
        MSTU(1)=3
        MSTU(2)=4
        DROBO(3)=(P(3,1)+P(4,1))/SHR
        DROBO(4)=(P(3,2)+P(4,2))/SHR
        CALL LUDBRB(MSTU(1),MSTU(2),0.,0.,-DROBO(3),-DROBO(4),0.D0)
        ROBO(2)=ULANGL(P(3,1),P(3,2))
        CALL LUDBRB(MSTU(1),MSTU(2),0.,-ROBO(2),0.D0,0.D0,0.D0)
        ROBO(1)=ULANGL(P(3,3),P(3,1))
        CALL LUDBRB(MSTU(1),MSTU(2),-ROBO(1),0.,0.D0,0.D0,0.D0)
        MSTU(2)=MAX(IPY(47),IPU1,IPU2)
        CALL LUDBRB(MSTU(1),MSTU(2),
     &  ROBO(1),ROBO(2),DROBO(3),DROBO(4),0.D0)
        DROBO(5)=MAX(-0.999999,MIN(0.999999,(PYVAR(31)-PYVAR(32))/
     &  (PYVAR(31)+PYVAR(32))))
        CALL LUDBRB(MSTU(1),MSTU(2),0.,0.,0.D0,0.D0,DROBO(5))
        MSTU(1)=0
        MSTU(2)=0
      ENDIF

C...CHECK INVARIANT MASS OF REMNANT SYSTEM:
C...HADRONIC EVENTS OR LEPTOPRODUCTION
      IF(ILEP.LE.0) THEN
        IF(IPY(12).LE.0.OR.ISUB.EQ.7) PYVAR(33)=0.
        IF(IPY(12).LE.0.OR.ISUB.EQ.7) PYVAR(34)=0.
        PEH=P(3,4)+P(4,4)+0.5*PYVAR(1)*(PYVAR(33)+PYVAR(34))
        PZH=P(3,3)+P(4,3)+0.5*PYVAR(1)*(PYVAR(33)-PYVAR(34))
        SHH=(PYVAR(1)-PEH)**2-(P(3,1)+P(4,1))**2-(P(3,2)+P(4,2))**2-
     &  PZH**2
        MSTJ(93)=1
        AMK32=ULMASS(K(3,2))
        MSTJ(93)=1
        AMK42=ULMASS(K(4,2))
        PMMIN=P(1,5)+P(2,5)+AMK32+AMK42
        IF(SHR.GE.PYVAR(1).OR.SHH.LE.(PMMIN+PYPAR(12))**2) THEN
          IPY(48)=1
          RETURN
        ENDIF
        SHR=SQRT(SHH+(P(3,1)+P(4,1))**2+(P(3,2)+P(4,2))**2)
      ELSE
        PEI=P(IQ,4)+P(IP,4)
        PZI=P(IQ,3)+P(IP,3)
        PMS(ILEP)=MAX(0.,PEI**2-PZI**2+P(5-ILEP,1)**2+P(5-ILEP,2)**2)
        MSTJ(93)=1
        PMMIN=P(3-ILEP,5)+ULMASS(K(5-ILEP,2))+SQRT(PMS(ILEP))
        IF(SHR.LE.PMMIN+PYPAR(12)) THEN
          IPY(48)=1
          RETURN
        ENDIF
      ENDIF

C...SUBDIVIDE REMNANT IF NECESSARY, STORE FIRST PARTON
  130 I=NS-1
      DO 160 JT=1,2
      IF(JT.EQ.ILEP) GOTO 160
      IF(JT.EQ.1) IPU=IPU1
      IF(JT.EQ.2) IPU=IPU2
      CALL LYSPLI(IPY(40+JT),KFL(1,JT),KFLCH(JT),KFLSP(JT))
      I=I+2
      IS(JT)=I
      K(I,1)=3
      K(I,2)=KFLSP(JT)
      K(I,3)=JT
      MSTJ(93)=1
      P(I,5)=ULMASS(K(I,2))
CJR--
      KFI2=LUCOMP(K(I,2))
      IF (KFI2.EQ.90) THEN
         P(I,5)=P(I,5)-2.*PARL(20)
      ELSEIF (1.LE.KFI2 .AND. KFI2.LE.6) THEN
         P(I,5)=P(I,5)-PARL(20)
      ENDIF
CJR--
C...FIRST PARTON COLOUR CONNECTIONS AND TRANSVERSE MASS
      K(I+1,1)=-1
      K(I+1,3)=I
      K(I+1,2)=1000
      IF(IPY(34).GE.1) K(I+1,2)=1000+JT
      DO 140 J=1,5
  140 P(I+1,J)=0.
      IF(KFLSP(JT).EQ.21) THEN
        P(I+1,3)=IPU
        P(I+1,4)=IPU
        P(IPU+1,1)=I
        P(IPU+1,2)=I
        K(I,4)=IPU+IPU*MSTU(5)
        K(I,5)=IPU+IPU*MSTU(5)
        K(IPU,4)=MOD(K(IPU,4),MSTU(5))+I*MSTU(5)
        K(IPU,5)=MOD(K(IPU,5),MSTU(5))+I*MSTU(5)
      ELSE
        IFLS=(3-ISIGN(1,KFLSP(JT)*(1102-IABS(KFLSP(JT)))))/2
        P(I+1,IFLS+2)=IPU
        P(IPU+1,3-IFLS)=I
        K(I,IFLS+3)=IPU
        K(IPU,6-IFLS)=MOD(K(IPU,6-IFLS),MSTU(5))+I*MSTU(5)
      ENDIF
      IF(KFLCH(JT).EQ.0) THEN
        P(I,1)=-P(JT+2,1)
        P(I,2)=-P(JT+2,2)
        PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2
      ELSE
C...WHEN EXTRA REMNANT PARTON OR HADRON: FIND RELATIVE PT, STORE
C...PRIMORDIAL KT SPLIT SHARED BETWEEN REMNANTS
        CALL LPRIKT(PARL(14),PTSPL,PHISPL)
C...RELATIVE DISTRIBUTION OF ENERGY; EXTRA PARTON COLOUR CONNECTION
        CALL LREMH(0,PTSPL,KFLSP(JT),KFLCH(JT),CHI(JT))
        P(I,1)=-P(JT+2,1)*(1.-CHI(JT))+PTSPL*COS(PHISPL)
        P(I,2)=-P(JT+2,2)*(1.-CHI(JT))+PTSPL*SIN(PHISPL)
        PMS(JT+2)=P(I,5)**2+P(I,1)**2+P(I,2)**2
        I=I+2
        DO 150 J=1,5
        K(I,J)=0
        K(I+1,J)=0
        P(I,J)=0.
  150   P(I+1,J)=0.
        K(I,1)=1
        K(I,2)=KFLCH(JT)
        K(I,3)=JT
        MSTJ(93)=1
        P(I,5)=ULMASS(K(I,2))
CJR--
        KFI2=LUCOMP(K(I,2))
        IF (KFI2.EQ.90) THEN
           P(I,5)=P(I,5)-2.*PARL(20)
        ELSEIF (1.LE.KFI2 .AND. KFI2.LE.6) THEN
           P(I,5)=P(I,5)-PARL(20)
        ENDIF
CJR--
        P(I,1)=-P(JT+2,1)*CHI(JT)-PTSPL*COS(PHISPL)
        P(I,2)=-P(JT+2,2)*CHI(JT)-PTSPL*SIN(PHISPL)
        PMS(JT+4)=P(I,5)**2+P(I,1)**2+P(I,2)**2
C...end of update
        PMS(JT)=PMS(JT+4)/CHI(JT)+PMS(JT+2)/(1.-CHI(JT))
        K(I+1,1)=-1
        K(I+1,3)=I
        K(I+1,2)=1000
        IF(IPY(34).GE.1) K(I+1,2)=1000+JT
        IF((IABS(KFLCH(JT)).GE.1.AND.IABS(KFLCH(JT)).LE.8).OR.
     &  IABS(KFLCH(JT)).EQ.21.OR.LUCOMP(IABS(KFLCH(JT))).EQ.90) THEN
          IFLS=(3-ISIGN(1,KFLCH(JT)*(1102-IABS(KFLCH(JT)))))/2
          P(I+1,IFLS+2)=IPU
          P(IPU+1,3-IFLS)=I
          K(I,1)=3
          K(I,IFLS+3)=IPU
          K(IPU,6-IFLS)=MOD(K(IPU,6-IFLS),MSTU(5))+I*MSTU(5)
        ELSE
          IF(IPY(34).GE.1) THEN
            K(I,1)=1
            K(I,3)=JT
          ENDIF
        ENDIF
      ENDIF
  160 CONTINUE
      IF(SHR.LE.SQRT(PMS(1))+SQRT(PMS(2))) GOTO 130
      N=I+1

C...RECONSTRUCT KINEMATICS OF REMNANTS
      DO 170 JT=1,2
      IF(JT.EQ.ILEP) GOTO 170
      PE=0.5*(SHR+(PMS(JT)-PMS(3-JT))/SHR)
      PZ=SQRT(PE**2-PMS(JT))
      IF(KFLCH(JT).EQ.0) THEN
        P(IS(JT),4)=PE
        P(IS(JT),3)=PZ*(-1)**(JT-1)
      ELSE
        PW1=CHI(JT)*(PE+PZ)
        P(IS(JT)+2,4)=0.5*(PW1+PMS(JT+4)/PW1)
        P(IS(JT)+2,3)=0.5*(PW1-PMS(JT+4)/PW1)*(-1)**(JT-1)
        P(IS(JT),4)=PE-P(IS(JT)+2,4)
        P(IS(JT),3)=PZ*(-1)**(JT-1)-P(IS(JT)+2,3)
      ENDIF
  170 CONTINUE

C     CALL GULIST(31,2)
C...HADRONIC EVENTS: BOOST REMNANTS TO CORRECT LONGITUDINAL FRAME
      IF(ILEP.LE.0) THEN
        MSTU(1)=NS+1
        CALL LUDBRB(MSTU(1),MSTU(2),
     &  0.,0.,0.D0,0.D0,-DBLE(PZH)/(DBLE(PYVAR(1))-DBLE(PEH)))
        MSTU(1)=0
C...LEPTOPRODUCTION EVENTS: BOOST COLLIDING SUBSYSTEM
      ELSE
        IMIN=21
        IMAX=MAX(IP,IPY(47))
        PEF=SHR-PE
        PZF=PZ*(-1)**(ILEP-1)
        PT2=P(5-ILEP,1)**2+P(5-ILEP,2)**2
        PHIPT=ULANGL(P(5-ILEP,1),P(5-ILEP,2))
        CALL LUDBRB(IMIN,IMAX,0.,-PHIPT,0.D0,0.D0,0.D0)
        RQP=P(IQ,3)*(PT2+PEI**2)-P(IQ,4)*PEI*PZI
        SINTH=P(IQ,4)*SQRT(PT2*(PT2+PEI**2)/(RQP**2+PT2*
     &  P(IQ,4)**2*PZI**2))*SIGN(1.,-RQP)
        CALL LUDBRB(IMIN,IMAX,ASIN(SINTH),0.,0.D0,0.D0,0.D0)
        DBETAX=(-DBLE(PEI)*PZI*SINTH+
     &  SQRT(DBLE(PT2)*(PT2+PEI**2-(PZI*SINTH)**2)))/
     &  (DBLE(PT2)+PEI**2)
        CALL LUDBRB(IMIN,IMAX,0.,0.,DBETAX,0.D0,0.D0)
        CALL LUDBRB(IMIN,IMAX,0.,PHIPT,0.D0,0.D0,0.D0)
        PEM=P(IQ,4)+P(IP,4)
        PZM=P(IQ,3)+P(IP,3)
        DBETAZ=(-DBLE(PEM)*PZM+
     &  PZF*SQRT(DBLE(PZF)**2+PEM**2-PZM**2))/(DBLE(PZF)**2+PEM**2)
        CALL LUDBRB(IMIN,IMAX,0.,0.,0.D0,0.D0,DBETAZ)
C...Avoid double application of kt
        P(4,1)=0.
        P(4,2)=0.
        CALL LUDBRB(3,4,ASIN(SINTH),0.,DBETAX,0.D0,0.D0)
        CALL LUDBRB(3,4,0.,PHIPT,0.D0,0.D0,DBETAZ)
      ENDIF

      RETURN
      END

C **********************************************************************

      SUBROUTINE LYSPLI(KPART,KFLIN,KFLCH,KFLSP)

C...IN CASE OF A HADRON REMNANT WHICH IS MORE COMPLICATED THAN JUST A
C...QUARK OR A DIQUARK, SPLIT IT INTO TWO (PARTONS OR HADRON + PARTON).
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP

      IFLIN=KFLIN
      KSIGN=ISIGN(1,KPART)
      IFL=KFLIN*KSIGN
      KFLCH=0
      IDUM=0

      IF(LST(14).EQ.0) THEN
C...If baryon production from remnant excluded, remnant is antiflavour
        KFLSP=-KFLIN
        IF(KFLIN.EQ.21) KFLSP=21
        RETURN
      ENDIF

      IF(IABS(KPART).EQ.211) THEN
C...DECOMPOSE PI+ (PI-).
        IF(IFL.EQ.2) THEN
C...VALENCE U (UBAR) REMOVED.
          KFLSP=-1*KSIGN
        ELSEIF(IFL.EQ.-1) THEN
C...VALENCE D (DBAR) REMOVED.
          KFLSP=2*KSIGN
        ELSEIF(KFLIN.EQ.21) THEN
C...GLUON REMOVED.
          R=2.*RLU(0)
          IF(R.LT.1.) THEN
            KFLCH=2*KSIGN
            KFLSP=-1*KSIGN
          ELSE
            KFLCH=-1*KSIGN
            KFLSP=2*KSIGN
          ENDIF
        ELSEIF((IFL.GE.1.AND.IFL.LE.8).AND.IFL.NE.2) THEN
C...SEA QUARK (ANTIQUARK) REMOVED.
          CALL LUKFDI(-IFLIN,2*KSIGN,IDUM,KFLCH)
          KFLSP=-1*KSIGN
        ELSEIF((IFL.GE.-8.AND.IFL.LE.-1).AND.IFL.NE.-1) THEN
C...SEA ANTIQUARK (QUARK) REMOVED.
          CALL LUKFDI(-IFLIN,-1*KSIGN,IDUM,KFLCH)
          KFLSP=2*KSIGN
        ENDIF

      ELSEIF(IABS(KPART).EQ.2212) THEN
C...DECOMPOSE PROTON (ANTIPROTON).
        IF(IFL.EQ.2) THEN
C...VALENCE U (UBAR) REMOVED.
          R=4.*RLU(0)
          IF(R.LT.3.) THEN
            KFLSP=2101*KSIGN
          ELSE
            KFLSP=2103*KSIGN
          ENDIF
        ELSEIF(IFL.EQ.1) THEN
C...VALENCE D (DBAR) REMOVED.
          KFLSP=2203*KSIGN
        ELSEIF(KFLIN.EQ.21) THEN
C...GLUON REMOVED.
          R=6.*RLU(0)
          IF(R.LT.3.) THEN
            KFLCH=2*KSIGN
            KFLSP=2101*KSIGN
          ELSEIF(R.LT.4.) THEN
            KFLCH=2*KSIGN
            KFLSP=2103*KSIGN
          ELSE
            KFLCH=1*KSIGN
            KFLSP=2203*KSIGN
          ENDIF
        ELSEIF(IFL.GT.2) THEN
C...SEA QUARK (ANTIQUARK) REMOVED.
          R=6*RLU(0)
          IF(R.LT.3.) THEN
            CALL LUKFDI(-IFLIN,2*KSIGN,IDUM,KFLCH)
            KFLSP=2101*KSIGN
          ELSEIF(R.LT.4.) THEN
            CALL LUKFDI(-IFLIN,2*KSIGN,IDUM,KFLCH)
            KFLSP=2103*KSIGN
          ELSE
            CALL LUKFDI(-IFLIN,1*KSIGN,IDUM,KFLCH)
            KFLSP=2203*KSIGN
          ENDIF
        ELSEIF(IFL.LT.0) THEN
C...SEA ANTIQUARK (QUARK) REMOVED.
  100     R=6*RLU(0)
          IF(R.LT.3.) THEN
            CALL LUKFDI(2101*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=2*KSIGN
          ELSEIF(R.LT.4.) THEN
            CALL LUKFDI(2103*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=2*KSIGN
          ELSE
            CALL LUKFDI(2203*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=1*KSIGN
          ENDIF
          IF(KFLCH.EQ.0) GOTO 100
        ENDIF

      ELSEIF(IABS(KPART).EQ.2112) THEN
C...DECOMPOSE NEUTRON (ANTINEUTRON).
        IF(IFL.EQ.2) THEN
C...VALENCE U (UBAR) REMOVED.
          KFLSP=1103*KSIGN
        ELSEIF(IFL.EQ.1) THEN
C...VALENCE D (DBAR) REMOVED.
          R=4.*RLU(0)
          IF(R.LT.3.) THEN
            KFLSP=2101*KSIGN
          ELSE
            KFLSP=2103*KSIGN
          ENDIF
        ELSEIF(KFLIN.EQ.21) THEN
C...GLUON REMOVED.
          R=6.*RLU(0)
          IF(R.LT.2.) THEN
            KFLCH=2*KSIGN
            KFLSP=1103*KSIGN
          ELSEIF(R.LT.5.) THEN
            KFLCH=1*KSIGN
            KFLSP=2101*KSIGN
          ELSE
            KFLCH=1*KSIGN
            KFLSP=2103*KSIGN
          ENDIF
        ELSEIF(IFL.GT.2) THEN
C...SEA QUARK (ANTIQUARK) REMOVED.
          R=6*RLU(0)
          IF(R.LT.2.) THEN
            CALL LUKFDI(-IFLIN,2*KSIGN,IDUM,KFLCH)
            KFLSP=1103*KSIGN
          ELSEIF(R.LT.5.) THEN
            CALL LUKFDI(-IFLIN,1*KSIGN,IDUM,KFLCH)
            KFLSP=2101*KSIGN
          ELSE
            CALL LUKFDI(-IFLIN,1*KSIGN,IDUM,KFLCH)
            KFLSP=2103*KSIGN
          ENDIF
        ELSEIF(IFL.LT.0) THEN
C...SEA ANTIQUARK (QUARK) REMOVED.
  110     R=6*RLU(0)
          IF(R.LT.2.) THEN
            CALL LUKFDI(1103*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=2*KSIGN
          ELSEIF(R.LT.5.) THEN
            CALL LUKFDI(2101*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=1*KSIGN
          ELSE
            CALL LUKFDI(2103*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=1*KSIGN
          ENDIF
          IF(KFLCH.EQ.0) GOTO 110
        ENDIF
      ENDIF

      RETURN
      END

C#######################################################################
C
C   The following routines are slightly modified minimization routines
C   from the MINUIT program package.
C
C **********************************************************************

      SUBROUTINE LMCMND

C...This is the MINUIT routine COMAND.
CC        GETS IFORMATION FROM /LMINUI/  AND TAKES APPROPRIATE ACTION,
CC        EITHER DIRECTLY BY SKIPPING TO THE CORRESPONDING CODE IN
CC        LMCMND, OR BY SETTING UP A CALL TO A SUBROUTINE
CC
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4),
     &MAXFIN,RELUP,RELERR,RELER2,FCNMAX
      COMMON /LPFLAG/ LST3

      COMMON
     1/LMMINE/ ERP(30)  ,ERN(30)
     2/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     3/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     4/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     5/LMVARI/ V(15,15)
     7/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     7/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     C/LMCASC/ Y(16)    ,JH       ,JL
     F/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     G/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     J/LMVART/ VT(15,15)
      COMMON
     6/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     8/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     9/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     A/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     B/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC
      FVAL3 = 2.0*AMIN+1.0
C                                        . . . . . . . . . . ERROR DEF
      WORD7(1)=RELUP*ABS(AMIN)
      UP = WORD7(1)
      IF (UP .LE. 0.)  UP = 1.0
      IF (ISW(2) .GE. 1)  CALL LMPRIN(1,AMIN)
      WORD7(1)=MAXFIN
      WORD7(2)=RELERR*UP
      NFCNMX = WORD7(1) + 0.5
      IF (NFCNMX .LE. 0)  NFCNMX = 1000
      EPSI = WORD7(2)
      IF (EPSI .LE. 0.)  EPSI = 0.1 * UP
      NEWMIN = 0
      ITAUR = 0
      ISW(1) = 0
      CALL LMSIMP
      IF(ABS(DIRIN(1)).LE.ABS(EPSMAC*X(1)).AND.
     &   ABS(DIRIN(2)).LE.ABS(EPSMAC*X(2))) THEN
        IF(LST3.GE.1) WRITE(6,2100)
        GOTO 500
      ENDIF
      WORD7(1)=MAXFIN
      RELERR=RELER2*RELERR
      WORD7(2)=RELERR*UP
      NFCNMX = WORD7(1) + 0.5
      IF (NFCNMX .LE. 0)  NFCNMX = 1000
      EPSI = WORD7(2)
      IF (EPSI .LE. 0.)  EPSI = 0.1 * UP
      CALL LMSIMP
  500 FCNMAX=ABS(AMIN)
      IF(ISW(1).GE.1) THEN
        IF(LST3.GE.1) WRITE(6,2200)
        FCNMAX=FCNMAX*1.25
      ENDIF
      FMAX=ABS(AMIN)
C                                        . . . . . . . . . . END, EXIT
      WORD7(1)=0.
 1100 IT = WORD7(1) + 0.5
      IF (FVAL3 .EQ. AMIN .OR. IT .GT. 0) RETURN
      IFLAG = 3
      CALL LSIGMX(NPAR,GIN,F,U,IFLAG)
      NFCN = NFCN + 1
      IF(LST3.GE.1.AND.ABS(F).GT.FMAX) WRITE(6,2300) F
      RETURN

 2100 FORMAT(' Warning: Stepsizes are less than machine accuracy ',
     &'times variable values. No further minimization attempted.')
 2200 FORMAT(' Warning: Simplex minimization has not converged ',
     &'properly.',/,10X,'Returned maximum increased by a factor 1.25.')
 2300 FORMAT(' Warning from LMCMND: function at minimum, ',E12.4,
     &'  is smaller than stored minimum.')

      END

C **********************************************************************

      SUBROUTINE LMINTO(PINT)

C...This is the MINUIT routine INTOEX.
CC        TRANSFORMS FROM INTERNAL COORDINATES (PINT) TO EXTERNAL
CC        PARAMETERS (U).   THE MINIMIZING ROUTINES WHICH WORK IN
CC        INTERNAL COORDINATES CALL THIS ROUTINE BEFORE CALLING FCN.

      COMMON
     1/LMMINE/ ERP(30)  ,ERN(30)
     2/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     3/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     4/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     5/LMVARI/ V(15,15)
     7/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     7/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     C/LMCASC/ Y(16)    ,JH       ,JL
     F/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     G/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     J/LMVART/ VT(15,15)
      COMMON
     6/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     8/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     9/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     A/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     B/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC

      DIMENSION PINT(2)
      DO 100 I= 1, NU
      J = LCORSP(I)
      IF ( J )  100,100,50
   50 CONTINUE
      IF (LCODE(I) .EQ. 1)  GO TO 80
      AL = ALIM(I)
      U(I) = AL + 0.5 *(SIN(PINT(J)) +1.0) * (BLIM(I) -AL)
      GO TO 100
   80 U(I) = PINT(J)
  100 CONTINUE
      RETURN
      END

C **********************************************************************

      SUBROUTINE LMIDAT

C...This is the MINUIT routine MIDATA.
CC        GETS PARAMETERS FROM /LMINUI/  AND /LMINUC/
CC        AND SETS UP THE STARTING PARAMETER LISTS.
CC        CONTROL THEN PASSES TO LMCMND FOR READING THE COMMAND "CARDS".
CC

      COMMON
     1/LMMINE/ ERP(30)  ,ERN(30)
     2/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     3/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     4/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     5/LMVARI/ V(15,15)
     7/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     7/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     C/LMCASC/ Y(16)    ,JH       ,JL
     F/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     G/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     J/LMVART/ VT(15,15)
      COMMON
     6/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     8/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     9/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     A/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     B/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4),
     &MAXFIN,RELUP,RELERR,RELER2,FCNMAX
      COMMON /LMINUC/ NAMKIN(4),NAM(30)
      COMMON /LPFLAG/ LST3
      CHARACTER*10 NAMKIN,NAM,NAMK,BLANK
      CHARACTER XTITLE*60
      REAL LMPINT
      DATA BLANK/'          '/
      DATA XTITLE/' FIND MINIMUM OF -(DIFFERENTIAL CROSS SECTION)'/
      DATA MNINIT/0/,IFATAL,NINT/0,0/
C                                        . INITIALIZE NEW DATA BLOCK . .
      IF (MNINIT .EQ. 0)  NBLOCK=0
      MNINIT = 1
      NBLOCK = NBLOCK + 1
      VERSN = 11.79
      IF(LST3.GE.5) THEN
        WRITE (ISYSWR,1004)  MAXINT,MAXEXT,VERSN,NBLOCK
        WRITE (ISYSWR,1005)
      ENDIF
      DO 50 I= 1, 7
   50 ISW(I) = 0
      SIGMA = 0.
      CALL LTIMEX(TIME)
      IF(LST3.GE.5) THEN
        WRITE (ISYSWR,1110)  XTITLE,TIME,EPSMAC
        WRITE (ISYSWR,1005)
      ENDIF
      NPFIX = 0
      NINT = 0
      NU = 0
      NPAR = 0
      IFATAL = 0
      IF(LST3.GE.5) WRITE (ISYSWR,1005)
      DO 100 I= 1, MAXEXT
      U(I) = 0.0
      NAM(I) = BLANK
      ERP(I) = 0.0
      ERN(I) = 0.0
      LCODE(I) = 0
  100 LCORSP (I) = 0
      UP = 1.0
      ISW(5) = 1
      IUNIT = ISYSRD
C                                        . . . READ PARAMETER CARDS . .
      ENTRY LMIDA2
      DO 200 I= 1, 200
      IF(I.GE.5) GOTO 250
      XK=XKIN(I)
      NAMK=NAMKIN(I)
      UK=UKIN(I)
      WK=WKIN(I)
      A=AIN(I)
      B=BIN(I)
      K = XK + 0.1
      NU = MAX0(NU,K)
      IF (K .LE. 0)  GO TO 250
      IF (K .LE. MAXEXT)  GO TO 115
      IFATAL = IFATAL + 1
      IF(LST3.GE.1) THEN
        WRITE (ISYSWR,1009)  K,MAXEXT
        WRITE (ISYSWR,1002)  K,NAMK,UK,WK,A,B
      ENDIF
      GO TO 200
  115 CONTINUE
      IF(NAM(K).EQ.BLANK) GO TO 120
C         PREVIOUSLY DEFINED PARAMETER IS BEING REDEFINED
      IF(LST3.GE.1) WRITE(ISYSWR,1007)
      IF(WERR(K).GT..0) NINT=NINT-1
  120 CONTINUE
      NAM(K) = NAMK
      U(K) = UK
      WERR(K) = WK
      IF (WK .GT. 0.0)  GO TO 122
C                                        . . . FIXED PARAMETER . . . .
      IF(LST3.GE.5) WRITE (ISYSWR, 1002)  K,NAMK,UK
      LCODE(K) = 0
      GO TO 160
C                                        . . . VARIABLE PARAMETER . . .
  122 IF(LST3.GE.5) WRITE (ISYSWR, 1002)  K,NAMK,UK,WK,A,B
      NINT = NINT + 1
      ISW(2) = 0
      IF  (A)  140,130,140
  130 IF  (B)  140,135,140
  135 LCODE(K) = 1
      GO TO 160
  140 IF (B-A)  145,142,150
  142 IFATAL = IFATAL + 1
      IF(LST3.GE.1) WRITE (ISYSWR,1010)
      GO TO 150
  145 SAV = B
      B = A
      A = SAV
      IF(LST3.GE.1) WRITE (ISYSWR,1003)
  150 ALIM(K) = A
      BLIM(K) = B
      LCODE(K) = 4
      IF ((B-U(K))*(U(K)-A))  153,155,160
  153 IFATAL = IFATAL + 1
      IF(LST3.GE.1) WRITE (ISYSWR,1011)
      GO TO 160
  155 IF(LST3.GE.1) WRITE (ISYSWR,1006)
  160 CONTINUE
  200 CONTINUE
      IFATAL = IFATAL + 1
      IF(LST3.GE.1) WRITE (ISYSWR,1012)
C                                       . . . END PARAMETER CARDS
C                                       . .   . STOP IF FATAL ERROR
  250 IF(LST3.GE.5) WRITE (ISYSWR,1005)
      IF (NINT .LE. MAXINT)  GO TO 253
      IF(LST3.GE.1) WRITE (ISYSWR,1008)  NINT,MAXINT
      IFATAL = IFATAL + 1
  253 IF (IFATAL .LE. 0)  GO TO 280
      IF(LST3.GE.1) WRITE (ISYSWR,1013)  IFATAL
      IF(LST3.GE.2) STOP
C                                       CALCULATE STEP SIZES DIRIN
  280 NPAR = 0
      DO 300 K= 1, NU
      IF (LCODE(K) .LE. 0)  GO TO 300
      NPAR = NPAR + 1
      LCORSP(K) = NPAR
      SAV = U(K)
      X(NPAR) = LMPINT(SAV,K)
      XT(NPAR) = X(NPAR)
      SAV2 = SAV + WERR(K)
      VPLU = LMPINT(SAV2,K) - X(NPAR)
      SAV2 = SAV - WERR(K)
      VMINU = LMPINT(SAV2,K) - X(NPAR)
      DIRIN(NPAR) = 0.5 * (ABS(VPLU) +ABS(VMINU))
      G2(NPAR) = 2.0 / DIRIN(NPAR)**2
      GSTEP(NPAR) = DIRIN(NPAR)
      IF (LCODE(K) .GT. 1)  GSTEP(NPAR) = -GSTEP(NPAR)
  300 CONTINUE
      SIGMA = 1.0E10
      IUNIT = ISYSRD
      RETURN
C...     THE FORMAT BELOW IS MACHINE-DEPENDENT. (A10) , (A6,4X) , ETC.
 1002 FORMAT (I10,2X,A10,2X,2G12.6,2X,2G12.6)
 1003 FORMAT(' WARNING           - ABOVE LIMITS HAVE BEEN REVERSED.')
 1004 FORMAT (1H1/42X,21(1H*)/42X,21H*   D506   MINUIT   */42X,
     112H* DIMENSIONS, I3, 1H/, I3, 2H */   42X,
     1'*  MODIFICATION OF  *',/,42X,
     111H*   VERSION ,F6.2,4H   */42X,16H* DATA BLOCK NO. ,I3,2H *)
 1005 FORMAT (4X,96(1H*))
 1006 FORMAT(' WARNING           - ABOVE PARAMETER IS AT LIMIT ')
 1007 FORMAT(' WARNING  *******  - PARAMETER REQUESTED ON FOLLOWING',
     1' CARD HAS ALREADY APPEARED.  PREVIOUS VALUES IGNORED.')
 1008 FORMAT('0   TOO MANY VARIABLE PARAMETERS.  YOU REQUEST',I5/,
     +'   THIS VERSION OF MINUIT IS ONLY DIMENSIONED FOR',I4//)
 1009 FORMAT('0FATAL ERROR. PARAMETER NUMBER',I11,' GREATER THAN ',
     +'ALLOWED MAXIMUM',I4)
 1010 FORMAT(' FATAL ERROR. UPPER AND LOWER LIMITS ARE EQUAL.')
 1011 FORMAT(' FATAL ERROR. PARAMETER OUTSIDE LIMITS',/)
 1012 FORMAT('0FATAL ERROR. MORE THAN 200 PARAMETER CARDS',/)
 1013 FORMAT(/I5,' FATAL ERRORS ON PARAMETER CARDS.  ABORT.',//)
 1110 FORMAT(5X,A60,5X,'TIME',F8.3,' SECONDS',/,70X,'MACH. PREC.=',
     &E10.2)
      END

C **********************************************************************

      SUBROUTINE LMINEW

C...This is the MINUIT routine MINNEW.
CC        THIS IS THE MAIN PROGRAM, DISGUISED AS A SUBROUTINE FOR
CC        REASONS OF COMPATIBILITY BETWEEN SYSTEMS.    IT INITIALIZES
CC        SOME CONSTANTS IN COMMON (INCLUDING THE LOGICAL I/O UNIT NOS.)
CC        THEN VERIFIES THAT FCN GIVES THE SAME VALUE WHEN CALLED
CC        TWICE WITH THE SAME ARGUMENTS, AND PASSES CONTROL TO LMCMND.
CC

      COMMON /LPFLAG/ LST3
      COMMON
     1/LMMINE/ ERP(30)  ,ERN(30)
     2/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     3/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     4/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     5/LMVARI/ V(15,15)
     7/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     7/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     C/LMCASC/ Y(16)    ,JH       ,JL
     F/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     G/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     J/LMVART/ VT(15,15)
      COMMON
     6/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     8/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     9/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     A/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     B/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC

C         UNIT NUMBERS FOR CARD READER, PRINTER, PUNCH
C
      ISYSRD = 5
      ISYSWR = 6
      ISYSPU = 7
      MAXINT=15
      MAXEXT=30
C                   DETERMINE MACHINE ACCURACY EPSMAC
      EPSMAC = 0.5
      DO 33 I= 1, 100
      EPSMAC = EPSMAC * 0.5
      IF ((1.0+EPSMAC) .EQ. 1.0)  GO TO 35
   33 CONTINUE
      EPSMAC = 1.0E-6
   35 EPSMAC = 2.0 * EPSMAC
C                             . . . . . . . . .
  110 CONTINUE
      NFCN = 1
      CALL LMIDAT
      CALL LMINTO(X)
      IF(LST3.GE.5) WRITE (ISYSWR,120)
  120 FORMAT (/,'0FIRST ENTRY TO FCN ')
      CALL LSIGMX(NPAR,GIN,AMIN,U,1)
      CALL LSIGMX(NPAR,GIN,AMIN,U,4)
      CALL LMPRIN(1,AMIN)
      CALL LSIGMX(NPAR,GIN,F   ,U,4)
      IF  (F .NE. AMIN)  GO TO 160
      NFCN = 3
      CALL LMCMND
      RETURN
  160 CONTINUE
      IF(LST3.GE.1) WRITE (ISYSWR,880) AMIN, F
      IF(LST3.GE.2) STOP
  880 FORMAT('0FOR THE ABOVE VALUES OF THE PARAMETERS, FCN IS TIME-',
     +'DEPENDENT',/,'0F = ',E22.14,' FOR FIRST CALL',/,' F =',E22.14,
     +' FOR SECOND')
      END

C **********************************************************************

      SUBROUTINE LMPRIN  (IKODE,FVAL)

C...This is the MINUIT routine MPRINT.
CC        PRINTS THE VALUES OF THE PARAMETERS AT THE TIME OF THE CALL.
CC        ALSO PRINTS OTHER RELEVANT INFORMATION SUCH AS FUNCTION VALUE,
CC        ESTIMATED DISTANCE TO MINIMUM, PARAMETER ERRORS, STEP SIZES.
CC        ACCORDING TO THE VALUE OF IKODE,THE PRINTOUT IS LONG FORMAT,
CC        SHORT FORMAT, OR MINOS FORMAT (0,1,2)
CC

      COMMON
     1/LMMINE/ ERP(30)  ,ERN(30)
     2/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     3/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     4/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     5/LMVARI/ V(15,15)
     7/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     7/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     C/LMCASC/ Y(16)    ,JH       ,JL
     F/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     G/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     J/LMVART/ VT(15,15)
      COMMON
     6/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     8/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     9/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     A/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     B/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC
      COMMON /LMINUC/ NAMKIN(4),NAM(30)
      COMMON /LPFLAG/ LST3
      CHARACTER*10 NAMKIN,NAM
C                                        . GET TIME AND PRINT HEADINGS .
      CALL LTIMEX(TI)
      IF(LST3.GE.5) WRITE (ISYSWR,1000)
      E = SIGMA
      KOUNT = 0
C                                        . . . LOOP OVER PARAMETERS . .
      DO 200 I= 1, NU
      IF(NAM(I).EQ.'          ') GOTO 200
   20 L = LCORSP(I)
      IF (L .EQ. 0)  GO TO 55
C              VARIABLE PARAMETER.  CALCULATE EXTERNAL ERROR IF V EXISTS
      IF (ISW(2) .LT. 1)  GO TO 27
      DX = SQRT(ABS(V(L,L)*UP))
      IF (LCODE(I) .LE. 1)  GO TO 26
      AL = ALIM(I)
      BA = BLIM(I) - AL
      DU1 = AL + 0.5 *(SIN(X(L)+DX) +1.0) * BA - U(I)
      DU2 = AL + 0.5 *(SIN(X(L)-DX) +1.0) * BA - U(I)
      IF (DX .GT. 1.0)  DU1 = BA
      DX = 0.5 * (ABS(DU1) + ABS(DU2))
   26 WERR(I) = DX
   27 X1 = X(L)
      X2 = DIRIN(L)
      IF (IKODE .LT. 2)  GO TO 29
      X1 = ERP(I)
      X2 = ERN(I)
   29 IF (KOUNT)  30,30,40
   30 KOUNT = 1
      IF(LST3.GE.5)
     &WRITE (ISYSWR,1001) FVAL,NFCN,TI,E, L,I,NAM(I),U(I),WERR(I),X1,X2
      GO TO 45
   40 IF(LST3.GE.5) WRITE (ISYSWR,1002)   L,I,NAM(I),U(I),WERR(I),X1,X2
   45 IF (LCODE(I) .LE. 1)  GO TO 200
      IF(LST3.GE.1.AND. ABS(COS(X(L))) .LT. 0.001)  WRITE (ISYSWR,1004)
      GO TO 200
C                           FIXED PARAMETER.  PRINT ONLY IF IKODE .GT.0
   55 IF (IKODE .EQ. 0)  GO TO 200
      IF (KOUNT)  60,60,70
   60 KOUNT = 1
      IF(LST3.GE.5) WRITE (ISYSWR,1001) FVAL,NFCN,TI,E, L,I,NAM(I),U(I)
      GO TO 200
   70 IF(LST3.GE.5) WRITE (ISYSWR,1003)                   I,NAM(I),U(I)
  200 CONTINUE
      IF(LST3.GE.5.AND.
     &IKODE.GE.1 .AND.ISW(2).GE.1)  WRITE (ISYSWR,1005)  UP
      RETURN
 1000 FORMAT(/ 4X,'FCN VALUE',5X,'CALLS',4X,'TIME',4X,' EDM  ',4X ,
     +'INT.EXT. PARAMETER     VALUE         ERROR      INTERN.VALUE  ',
     +'INT.STEP SIZE')
 1001 FORMAT(E15.7,I7,F9.2,E11.2,I6,I4,1X,A10,4E14.5)
 1002 FORMAT(1H ,41X,I6,I4,1X,A10,4E14.5)
 1003 FORMAT(1H  ,47X  ,I4,1X,A10,4E14.5)
 1004 FORMAT(1H ,52X  ,'WARNING -   - ABOVE PARAMETER IS AT LIMIT.')
 1005 FORMAT(/45X,'ERRORS CORRESPOND TO FUNCTION CHANGE OF ',E12.4)
      END

C **********************************************************************

      REAL FUNCTION LMPINT(PEXTI,I)

C...This is the MINUIT routine PINTF.
CC        CALCULATES THE INTERNAL PARAMETER VALUE LMPINT CORRESPONDING
CC        TO THE EXTERNAL VALUE PEXTI FOR PARAMETER I.
CC
      COMMON
     1/LMMINE/ ERP(30)  ,ERN(30)
     2/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     3/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     4/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     5/LMVARI/ V(15,15)
     7/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     7/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     C/LMCASC/ Y(16)    ,JH       ,JL
     F/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     G/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     J/LMVART/ VT(15,15)
      COMMON
     6/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     8/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     9/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     A/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     B/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC
      COMMON /LPFLAG/ LST3
      DATA BIG, SMALL / 1.570796326795 , -1.570796326795 /
      IGO = LCODE(I)
      GO TO (100,200,300,400),IGO
C--       IGO = 1  MEANS NO LIMITS
  100 LMPINT = PEXTI
      GO TO 800
  200 CONTINUE
  300 CONTINUE
C--       IGO = 4  MEANS THERE ARE TWO LIMITS
  400 ALIMI = ALIM(I)
      BLIMI = BLIM(I)
      IF (PEXTI-ALIMI)  440,500,460
  440 A = SMALL
  450 LMPINT = A
      PEXTI = ALIMI + 0.5* (BLIMI-ALIMI) *(SIN(A) +1.0)
      LIMSET=1
      IF(LST3.GE.1) WRITE (ISYSWR,241)  I
      GO TO 800
  460 IF (BLIMI-PEXTI)  470,520,480
  470 A = BIG
      GO TO 450
  480 YY=2.0*(PEXTI-ALIMI)/(BLIMI-ALIMI) - 1.0
      LMPINT = ATAN(YY/SQRT(1.0- YY**2) )
      GO TO 800
  500 LMPINT = SMALL
      GO TO 800
  520 LMPINT = BIG
  800 RETURN
  241 FORMAT(' WARNING - VARIABLE',I3,' HAS BEEN BROUGHT BACK IN',
     +'SIDE LIMITS BY LMPINT.')
      END

C **********************************************************************

      SUBROUTINE LMRAZZ(YNEW,PNEW)

C...This is the MINUIT routine RAZZIA.
CC        CALLED ONLY BY SIMPLEX (AND IMPROV) TO ADD A NEW POINT
CC        AND REMOVE AN OLD ONE FROM THE CURRENT SIMPLEX, AND GET THE
CC        ESTIMATED DISTANCE TO MINIMUM.
CC
      COMMON
     1/LMMINE/ ERP(30)  ,ERN(30)
     2/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     3/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     4/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     5/LMVARI/ V(15,15)
     7/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     7/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     C/LMCASC/ Y(16)    ,JH       ,JL
     F/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     G/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     J/LMVART/ VT(15,15)
      COMMON
     6/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     8/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     9/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     A/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     B/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC
      COMMON /LPFLAG/ LST3
      DIMENSION PNEW(15)
      DO 10 I=1,NPAR
10    P(I,JH)=PNEW(I)
      Y(JH)=YNEW
      IF(YNEW.GE.AMIN) GO TO 18
      DO 15 I=1,NPAR
15    X(I)=PNEW(I)
      CALL LMINTO(X)
      AMIN=YNEW
      JL=JH
18    CONTINUE
      JH=1
      NPARP1=NPAR+1
20    DO 25 J=2,NPARP1
      IF (Y(J) .GT. Y(JH))  JH = J
25    CONTINUE
      SIGMA = Y(JH) - Y(JL)
      IF (SIGMA .LE. 0.)  GO TO 45
      US = 1.0/SIGMA
      DO 35 I= 1, NPAR
      PBIG = P(I,1)
      PLIT = PBIG
      DO 30 J= 2, NPARP1
      IF (P(I,J) .GT. PBIG)  PBIG = P(I,J)
      IF (P(I,J) .LT. PLIT)  PLIT = P(I,J)
   30 CONTINUE
      DIRIN(I) = PBIG - PLIT
      IF (ITAUR .LT. 1 )  V(I,I) = 0.5*(V(I,I) +US*DIRIN(I)**2)
   35 CONTINUE
   40 RETURN
   45 IF(LST3.GE.1) WRITE (ISYSWR, 1000)  NPAR
      GO TO 40
 1000 FORMAT('0***** FUNCTION VALUE DOES NOT SEEM TO DEPEND ON ANY ',
     +'OF THE',I3,' VARIABLE PARAMETERS',/15X ,'VERIFY THAT STEP SIZES',
     +' ARE BIG ENOUGH AND CHECK FCN LOGIC.',/1X,81(1H*)/1X,81(1H*)//)
      END

C **********************************************************************

      SUBROUTINE LMSIMP

C...This is the MINUIT routine SIMPLEX.
CC        PERFORMS A MINIMIZATION USING THE SIMPLEX METHOD OF NELDER
CC        AND MEAD (REF. -- COMP. J. 7,308 (1965)).
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4),
     &MAXFIN,RELUP,RELERR,RELER2,FCNMAX
      COMMON /LPFLAG/ LST3
      COMMON
     1/LMMINE/ ERP(30)  ,ERN(30)
     2/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     3/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     4/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     5/LMVARI/ V(15,15)
     7/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     7/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     C/LMCASC/ Y(16)    ,JH       ,JL
     F/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     G/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     J/LMVART/ VT(15,15)
      COMMON
     6/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     8/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     9/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     A/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     B/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC

      DATA ALPHA,BETA,GAMMA,RHOMIN,RHOMAX / 1.0, 0.5, 2.0, 4.0, 8.0/
      IF (NPAR .LE. 0)  RETURN
      NPFN=NFCN
      NPARP1=NPAR+1
      RHO1 = 1.0 + ALPHA
      RHO2 = RHO1 + ALPHA*GAMMA
      WG = 1.0/FLOAT(NPAR)
      IFLAG=4
      IF(LST3.GE.5) WRITE(ISYSWR,100) EPSI
      DO 2 I= 1, NPAR
      IF (ISW(2) .GE. 1)  DIRIN(I) = SQRT(V(I,I)*UP)
      IF (ABS(DIRIN(I)) .LT. 1.0E-10*ABS(X(I))) DIRIN(I)=1.0E-8*X(I)
      IF(ITAUR.LT. 1)  V(I,I) =    DIRIN(I)**2/UP
    2 CONTINUE
      IF (ITAUR .LT. 1)  ISW(2) = 1
C**       CHOOSE THE INITIAL SIMPLEX USING SINGLE-PARAMETER SEARCHES
    1 CONTINUE
      YNPP1 = AMIN
      JL = NPARP1
      Y(NPARP1) = AMIN
      ABSMIN = AMIN
      DO 10 I= 1, NPAR
      AMING = AMIN
      PBAR(I) = X(I)
      BESTX = X(I)
      KG = 0
      NS = 0
      NF = 0
    4 X(I) = BESTX + DIRIN(I)
      CALL LMINTO(X)
      CALL LSIGMX(NPAR,GIN, F, U, 4)
      NFCN = NFCN + 1
      IF (F .LE. AMING)  GO TO 6
C         FAILURE
      IF (KG .EQ. 1)  GO TO 8
      KG = -1
      NF = NF + 1
      DIRIN(I) = DIRIN(I) * (-0.4)
      IF (NF .LT. 3)  GO TO 4
      NS = 6
C         SUCCESS
    6 BESTX = X(I)
      DIRIN(I) = DIRIN(I) * 3.0
      AMING = F
      KG = 1
      NS = NS + 1
      IF (NS .LT. 6)  GO TO 4
C         LOCAL MINIMUM FOUND IN ITH DIRECTION
    8 Y(I) = AMING
      IF (AMING .LT. ABSMIN)  JL = I
      IF (AMING .LT. ABSMIN)  ABSMIN = AMING
      X(I) = BESTX
      DO 9 K= 1, NPAR
    9 P(K,I) = X(K)
   10 CONTINUE
      JH = NPARP1
      AMIN=Y(JL)
      CALL LMRAZZ(YNPP1,PBAR)
      DO 20 I= 1, NPAR
   20 X(I) = P(I,JL)
      CALL LMINTO(X)
      DO 30 I=1,NPAR
   30 IF(ABS(DIRIN(I)).LE.ABS(EPSMAC*X(I))) DIRIN(I)=4.*EPSMAC*X(I)
      IF (ISW(5) .GE. 1)  CALL LMPRIN(0,AMIN)
      SIGMA = SIGMA * 10.
      SIG2 = SIGMA
      NCYCL=0
C                                        . . . . .  START MAIN LOOP
   50 CONTINUE
C...Change in SIMPLX; error redefined for second call to LMSIMP.
      UP=RELUP*ABS(AMIN)
      EPSI=RELERR*UP
      IF (SIG2 .LT. EPSI .AND. SIGMA.LT.EPSI)     GO TO 76
      SIG2 = SIGMA
      IF ((NFCN-NPFN) .GT. NFCNMX)  GO TO 78
C         CALCULATE NEW POINT * BY REFLECTION
      DO 60 I= 1, NPAR
      PB = 0.
      DO 59 J= 1, NPARP1
   59 PB = PB + WG * P(I,J)
      PBAR(I) = PB - WG * P(I,JH)
   60 PSTAR(I)=(1.+ALPHA)*PBAR(I)-ALPHA*P(I,JH)
      CALL LMINTO(PSTAR)
      CALL LSIGMX(NPAR,GIN,YSTAR,U,4)
      NFCN=NFCN+1
      IF(YSTAR.GE.AMIN) GO TO 70
C         POINT * BETTER THAN JL, CALCULATE NEW POINT **
      DO 61 I=1,NPAR
   61 PSTST(I)=GAMMA*PSTAR(I)+(1.-GAMMA)*PBAR(I)
      CALL LMINTO(PSTST)
      CALL LSIGMX(NPAR,GIN,YSTST,U,4)
      NFCN=NFCN+1
C         TRY A PARABOLA THROUGH PH, PSTAR, PSTST.  MIN = PRHO
      Y1 = (YSTAR-Y(JH)) * RHO2
      Y2 = (YSTST-Y(JH)) * RHO1
      RHO = 0.5 * (RHO2*Y1 -RHO1*Y2) / (Y1 -Y2)
      IF (RHO .LT. RHOMIN)  GO TO 66
      IF (RHO .GT. RHOMAX)  RHO = RHOMAX
      DO 64 I= 1, NPAR
   64 PRHO(I) = RHO*PBAR(I) + (1.0-RHO)*P(I,JH)
      CALL LMINTO(PRHO)
      CALL LSIGMX(NPAR,GIN,YRHO, U,4)
      NFCN = NFCN + 1
      IF (YRHO .LT. Y(JL) .AND. YRHO .LT. YSTST)  GO TO 65
      IF (YSTST .LT. Y(JL))  GO TO 67
      IF (YRHO .GT. Y(JL))  GO TO 66
C         ACCEPT MINIMUM POINT OF PARABOLA, PRHO
   65 CALL LMRAZZ (YRHO,PRHO)
      GO TO 68
   66 IF (YSTST .LT. Y(JL))  GO TO 67
      CALL LMRAZZ(YSTAR,PSTAR)
      GO TO 68
   67 CALL LMRAZZ(YSTST,PSTST)
   68 NCYCL=NCYCL+1
      IF (ISW(5) .LT. 2)  GO TO 50
      IF (ISW(5) .GE. 3 .OR. MOD(NCYCL, 10) .EQ. 0) CALL LMPRIN(0,AMIN)
      GO TO 50
C         POINT * IS NOT AS GOOD AS JL
   70 IF (YSTAR .GE. Y(JH))  GO TO 73
      JHOLD = JH
      CALL LMRAZZ(YSTAR,PSTAR)
      IF (JHOLD .NE. JH)  GO TO 50
C         CALCULATE NEW POINT **
   73 DO 74 I=1,NPAR
   74 PSTST(I)=BETA*P(I,JH)+(1.-BETA)*PBAR(I)
      CALL LMINTO (PSTST)
      CALL LSIGMX(NPAR,GIN,YSTST,U,4)
      NFCN=NFCN+1
      IF(YSTST.GT.Y(JH)) GO TO 1
C     POINT ** IS BETTER THAN JH
      IF (YSTST .LT. AMIN)  GO TO 67
      CALL LMRAZZ(YSTST,PSTST)
      GO TO 50
C                                        . . . . . .  END MAIN LOOP
   76 IF(LST3.GE.5) WRITE(ISYSWR,120)
      GO TO 80
   78 IF(LST3.GE.5) WRITE(ISYSWR,130)
      ISW(1) = 1
   80 DO 82 I=1,NPAR
      PB = 0.
      DO 81 J=1,NPARP1
   81 PB = PB + WG * P(I,J)
   82 PBAR(I) = PB - WG * P(I,JH)
      CALL LMINTO(PBAR)
      CALL LSIGMX(NPAR,GIN,YPBAR,U,IFLAG)
      NFCN=NFCN+1
      IF (YPBAR .LT. AMIN)  CALL LMRAZZ(YPBAR,PBAR)
      CALL LMINTO(X)
      IF (NFCNMX+NPFN-NFCN .LT. 3*NPAR)  GO TO 90
      IF (SIGMA .GT. 2.0*EPSI)  GO TO 1
   90 CALL LMPRIN(1-ITAUR, AMIN)
      RETURN
  100 FORMAT(' START SIMPLEX MINIMIZATION          ',8X   ,'CON',
     +'VERGENCE CRITERION -- ESTIMATED DISTANCE TO MINIMUM (EDM) .LT.',
     +E10.2 )
120   FORMAT(1H ,'SIMPLEX MINIMIZATION HAS CONVERGED')
130   FORMAT(1H ,'SIMPLEX TERMINATES WITHOUT CONVERGENCE')
      END

C#######################################################################
C
C   One- and two-dimensional adaptive Gaussian integration routines.
C
C **********************************************************************

      SUBROUTINE GADAP(A0,B0,F,EPS,SUM)
C
C   PURPOSE           - INTEGRATE A FUNCTION F(X)
C   METHOD            - ADAPTIVE GAUSSIAN
C   USAGE             - CALL GADAP(A0,B0,F,EPS,SUM)
C   PARAMETERS  A0    - LOWER LIMIT (INPUT,REAL)
C               B0    - UPPER LIMIT (INPUT,REAL)
C               F     - FUNCTION F(X) TO BE INTEGRATED. MUST BE
C                       SUPPLIED BY THE USER. (INPUT,REAL FUNCTION)
C               EPS   - DESIRED RELATIVE ACCURACY. IF SUM IS SMALL EPS
C                       WILL BE ABSOLUTE ACCURACY INSTEAD. (INPUT,REAL)
C               SUM   - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL)
C   PRECISION         - SINGLE
C   REQ'D PROG'S      - F
C   AUTHOR            - T. JOHANSSON, LUND UNIV. COMPUTER CENTER, 1973
C   REFERENCE(S)      - THE AUSTRALIAN COMPUTER JOURNAL,3 P.126 AUG. -71
C
      COMMON/GADAP1/ NUM,IFU
      EXTERNAL F
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)
    1 FORMAT(16H GADAP:I TOO BIG)
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.3
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)
      F2(1)=F(0.5*(A0+B0))
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)
      IFU=3
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)
  100 CONTINUE
      L=L+1
      N(L)=3
      EPS=EPS*RED
      A(I+1)=A(I)+C*(B(I)-A(I))
      B(I+1)=B(I)
      A(I+2)=A(I)+B(I)-A(I+1)
      B(I+2)=A(I+1)
      A(I+3)=A(I)
      B(I+3)=A(I+2)
      W1=A(I)+(B(I)-A(I))/5.
      U2=2.*W1-(A(I)+A(I+2))/2.
      F1(I+1)=F(A(I)+B(I)-W1)
      F2(I+1)=F3(I)
      F3(I+1)=F(B(I)-A(I+2)+W1)
      F1(I+2)=F(U2)
      F2(I+2)=F2(I)
      F3(I+2)=F(B(I+2)+A(I+2)-U2)
      F1(I+3)=F(A(I)+A(I+2)-W1)
      F2(I+3)=F1(I)
      F3(I+3)=F(W1)
      IFU=IFU+6
      IF(IFU.GT.5000) GOTO 130
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))
      SS=S(I+1)+S(I+2)+S(I+3)
      I=I+3
      IF(I.GT.300)GOTO 120
      SOLD=S(I-3)
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100
      SUM=SUM+SS
      I=I-4
      N(L)=0
      L=L-1
  110 CONTINUE
      IF(L.EQ.1) GOTO 130
      N(L)=N(L)-1
      EPS=EPS/RED
      IF(N(L).NE.0) GOTO 100
      I=I-1
      L=L-1
      GOTO 110
  120 WRITE(6,1)
  130 RETURN
      END

C **********************************************************************

      SUBROUTINE GADAP2(A0,B0,FL,FU,F,EPS,SUM)
C
C   PURPOSE           - INTEGRATE A FUNCTION F(X,Y) OF TWO VARIABLES
C   METHOD            - ADAPTIVE GAUSSIAN IN BOTH DIRECTIONS
C   USAGE             - CALL GADAP2(A0,B0,FL,FU,F,EPS,SUM)
C   PARAMETERS  A0    - LOWER X-LIMIT (INPUT,REAL)
C               B0    - UPPER X-LIMIT (INPUT,REAL)
C               FL    - USER SUPPLIED FUNCTION FL(X) GIVING THE LOWER
C                       Y-LIMIT FOR A GIVEN X-VALUE
C                       (INPUT,REAL FUNCTION)
C               FU    - USER SUPPLIED FUNCTION FU(X) GIVING THE UPPER
C                       Y-LIMIT FOR A GIVEN X-VALUE
C                       (INPUT,REAL FUNCTION)
C               F     - USER SUPPLIED FUNCTION F(X,Y) TO BE INTEGRATED
C                       (INPUT,REAL FUNCTION)
C               EPS   - DESIRED ACCURACY (INPUT,REAL)
C               SUM   - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL)
C   PRECISION         - SINGLE
C   REQ'D PROG'S      - FL,FU,F,GADAPF
C   AUTHOR            - THOMAS JOHANSSON, LDC,1973
C
      COMMON/GADAP1/ NUM,IFU
      EXTERNAL F,FL,FU
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)
    1 FORMAT(16H GADAP:I TOO BIG)
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.4
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      X=0.5*(1+C)*A0+0.5*(1-C)*B0
      AY=FL(X)
      BY=FU(X)
      F1(1)=GADAPF(X,AY,BY,F,EPS)
      X=0.5*(A0+B0)
      AY=FL(X)
      BY=FU(X)
      F2(1)=GADAPF(X,AY,BY,F,EPS)
      X=0.5*(1-C)*A0+0.5*(1+C)*B0
      AY=FL(X)
      BY=FU(X)
      F3(1)=GADAPF(X,AY,BY,F,EPS)
      IFU=3
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)
  100 CONTINUE
      L=L+1
      N(L)=3
      EPS=EPS*RED
      A(I+1)=A(I)+C*(B(I)-A(I))
      B(I+1)=B(I)
      A(I+2)=A(I)+B(I)-A(I+1)
      B(I+2)=A(I+1)
      A(I+3)=A(I)
      B(I+3)=A(I+2)
      W1=A(I)+(B(I)-A(I))/5.
      U2=2.*W1-(A(I)+A(I+2))/2.
      X=A(I)+B(I)-W1
      AY=FL(X)
      BY=FU(X)
      F1(I+1)=GADAPF(X,AY,BY,F,EPS)
      F2(I+1)=F3(I)
      X=B(I)-A(I+2)+W1
      AY=FL(X)
      BY=FU(X)
      F3(I+1)=GADAPF(X,AY,BY,F,EPS)
      X=U2
      AY=FL(X)
      BY=FU(X)
      F1(I+2)=GADAPF(X,AY,BY,F,EPS)
      F2(I+2)=F2(I)
      X=B(I+2)+A(I+2)-U2
      AY=FL(X)
      BY=FU(X)
      F3(I+2)=GADAPF(X,AY,BY,F,EPS)
      X=A(I)+A(I+2)-W1
      AY=FL(X)
      BY=FU(X)
      F1(I+3)=GADAPF(X,AY,BY,F,EPS)
      F2(I+3)=F1(I)
      X=W1
      AY=FL(X)
      BY=FU(X)
      F3(I+3)=GADAPF(X,AY,BY,F,EPS)
      IFU=IFU+6
      IF(IFU.GT.5000) GOTO 130
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))
      SS=S(I+1)+S(I+2)+S(I+3)
      I=I+3
      IF(I.GT.300)GOTO 120
      SOLD=S(I-3)
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100
      SUM=SUM+SS
      I=I-4
      N(L)=0
      L=L-1
  110 CONTINUE
      IF(L.EQ.1) GOTO 130
      N(L)=N(L)-1
      EPS=EPS/RED
      IF(N(L).NE.0) GOTO 100
      I=I-1
      L=L-1
      GOTO 110
  120 WRITE(6,1)
 130  RETURN
      END

C **********************************************************************


      FUNCTION GADAPF(X,A0,B0,F,EPS)
      COMMON/GADAP1/ NUM,IFU
      EXTERNAL F
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)
    1 FORMAT(16H GADAP:I TOO BIG)
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.4
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      F1(1)=F(X,0.5*(1+C)*A0+0.5*(1-C)*B0)
      F2(1)=F(X,0.5*(A0+B0))
      F3(1)=F(X,0.5*(1-C)*A0+0.5*(1+C)*B0)
      IFU=3
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)
  100 CONTINUE
      L=L+1
      N(L)=3
      EPS=EPS*RED
      A(I+1)=A(I)+C*(B(I)-A(I))
      B(I+1)=B(I)
      A(I+2)=A(I)+B(I)-A(I+1)
      B(I+2)=A(I+1)
      A(I+3)=A(I)
      B(I+3)=A(I+2)
      W1=A(I)+(B(I)-A(I))/5.
      U2=2.*W1-(A(I)+A(I+2))/2.
      F1(I+1)=F(X,A(I)+B(I)-W1)
      F2(I+1)=F3(I)
      F3(I+1)=F(X,B(I)-A(I+2)+W1)
      F1(I+2)=F(X,U2)
      F2(I+2)=F2(I)
      F3(I+2)=F(X,B(I+2)+A(I+2)-U2)
      F1(I+3)=F(X,A(I)+A(I+2)-W1)
      F2(I+3)=F1(I)
      F3(I+3)=F(X,W1)
      IFU=IFU+6
      IF(IFU.GT.5000) GOTO 130
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))
      SS=S(I+1)+S(I+2)+S(I+3)
      I=I+3
      IF(I.GT.300)GOTO 120
      SOLD=S(I-3)
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100
      SUM=SUM+SS
      I=I-4
      N(L)=0
      L=L-1
  110 CONTINUE
      IF(L.EQ.1) GOTO 130
      N(L)=N(L)-1
      EPS=EPS/RED
      IF(N(L).NE.0) GOTO 100
      I=I-1
      L=L-1
      GOTO 110
  120 WRITE(6,1)
  130 GADAPF=SUM
      EPS=EPS/RED
      RETURN
      END

C ********************************************************************

      SUBROUTINE LNSTRF(X,Q2,XPQ)

C...Structure function per nucleon for a proton/neutron mixture
C...according to defined nucleus.

      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      DIMENSION XPQ(-6:6)

      CALL LYSTFU(2212,X,Q2,XPQ)

      IF(PARI(11).LE.1.E-06) RETURN
      XDV=XPQ(1)-XPQ(-1)
      XUV=XPQ(2)-XPQ(-2)
C...For nuclear target, mix u- and d-valence distributions.
      XPQ(1)=(1.-PARI(11))*XDV+PARI(11)*XUV + XPQ(-1)
      XPQ(2)=(1.-PARI(11))*XUV+PARI(11)*XDV + XPQ(-2)
C...Save d and u valence in proton
      PARI(12)=XDV
      PARI(13)=XUV

      RETURN
      END

C **********************************************************************

      SUBROUTINE LYSTFU(KF,X,Q2,XPQ)

C...Interface to PYSTFU in PYTHIA 5.7 to get parton density distributions,
C...i.e. momentum weighted probability distributions xq(x,Q2), xg(x,Q2).
C...Also gives intrinsic charm and beauty distributions.
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON /ARSTRF/ KFSAVE(2),XSAVE(2),XQ2SAV(2),XPQSAV(2,-6:6)
      DIMENSION XPQ(-6:6),XPYST(-25:25)

C...Reset arrays etc.
      DO 100 KFL=-6,6
      XPQ(KFL)=0.0
  100 XPQSAV(1,KFL)=0.
      XSAVE(1)=X
      XQ2SAV(1)=Q2
      KFSAVE(1)=KF
C...Check x and particle species.
      IF(X.LE.0..OR.X.GE.1.) THEN
        WRITE(MSTU(11),5000) X
        RETURN
      ENDIF

      IF(LST(15).EQ.-4.OR.LST(15).EQ.-5) THEN
C...Intrinsic charm/bottom quark distribution in the proton...
        IF(Q2.LT.1.) RETURN
C...from Phys. Lett 93B (1980) 451, 
C...Amount of intrinsic charm PARL(12)=BETA^2
        XPQ(4)=X**3*1800.*PARL(12)*
     &         ((1.-X)/3.*(1.+10.*X+X**2)+2.*X*(1.+X)*LOG(X))
C...plus first order QCD-correction parametrized with polynomia
        IF(X.LT.0.9) THEN
          XCORR=0.22024E-1*X-0.77833E-1*X**2-0.47292*X**3+
     &          2.104*X**4-2.1698*X**5-0.84891*X**6+1.8882*X**7+
     &          0.8989*X**8-2.1072*X**9+0.76351*X**10
        ELSE
          XCORR=-1.
        ENDIF
C...and a Q2 dependence on that
CJR        XCORR=1.125*XCORR*0.190424*EXP(1.15*LOG(LOG(Q2)))
        IF(Q2.GT.1) THEN
           XCORR=1.125*XCORR*0.190424*EXP(1.15*LOG(LOG(Q2)))
        ELSE
           XCORR=1.125*XCORR*0.190424
        ENDIF
C...smooth cut-off of the structure function !
        XPQ(4)=MAX(XPQ(4)+XCORR,XPQ(4)/Q2)
        XPQ(-4)=XPQ(4)
        IF(LST(15).EQ.-5) THEN
C...Intrinsic bottom assumed to have the same shape as zeroth
C...approximation but suppressed by (mc/mb)**2=0.1 approximately
          XPQ(5)=XPQ(4)*0.1
          XPQ(-5)=XPQ(5)
          XPQ(4)=0.
          XPQ(-4)=0.
        ENDIF
      ELSE
C...Parton densities from PYSTFU in PYTHIA 5.7
        CALL PYSTFU(KF,X,Q2,XPYST)
        DO 110 KFL=-6,6
  110   XPQ(KFL)=XPYST(KFL)
      ENDIF

      DO 120 KFL=-6,6
  120 XPQSAV(1,KFL)=XPQ(KFL)

C...Formats for error printouts.
 5000 FORMAT(' Error in LYSTFU: x =',1P,E12.4,' outside physical range')

      RETURN
      END
