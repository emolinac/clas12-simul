C**********************************************************************
C...  
C...  PROGRAM QP1
      SUBROUTINE QP1
c     IMPLICIT NONE
      INTEGER NSET
      PARAMETER(NSET=20)
C...Demonstration job for LEPTO
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      
      REAL CUT,PARL,X,Y,W2,Q2,U

c      INTEGER LST
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDATR/MRLU(6),RRLU(100)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
c     COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
c     COMMON/PYSUBS/MSEL, MSLEPD, MSUB(500), KFIN(2,-40:40), CKIN(200)

c      REAL MRLU, RRLU
      integer*4 today(3), now(3)
c      REAL P,V
      DOUBLE PRECISION P,V
      integer partpid
      INTEGER N,K
      INTEGER NEVENT
      REAL TarA,TarZ

      CHARACTER OUTNAME*20

      INTEGER NE,SET,Npid,Neta
      REAL T1,T2,T3
      real plz,ppz
      logical check_pid

c     Beam energy
      plz = 10.6

      open(unit=31,file="lepto_input.txt", status="old")
      read(31,*) NEVENT,TarA,TarZ
      close(31)

c      mstj(21) = 0  !all particle decays are inhibited (PYTHIA, PYDAT1)

c       MDCY(LUCOMP(111),1)=0     !inhibit pi0 decay (JETSET, LUDAT3)
c       MDCY(LUCOMP(221),1)=0     !inhibit eta decay
c       MDCY(LUCOMP(331),1)=0     !inhibit eta' decay
c       MDCY(LUCOMP(321),1)=0     !inhibit K+ decay
c       MDCY(LUCOMP(310),1)=0     !inhibit K0s decay
c       MDCY(LUCOMP(-323),1)=0    !inhibit K*- decay
c       MDCY(LUCOMP(3122),1)=0     !inhibit Lambda
c       MDCY(LUCOMP(223),1)=0     !inhibit Omega
     
      lst(1)  =  1   ! choice of indep. variables; x and Q2
      lst(2)  =  1   ! D1; kin variables chosen fron diff cross section
      lst(3)  =  5   ! D5; verbosity of output; full output(5)
      lst(5)  =  3   ! which reference frame; lab system(3), CM(1)
c     lst(7)  =  0   !!!D = 1 hadron decays in Lepto; 0 no h decay/shower(later in GSIM) 
      lst(8)  =  9   ! (1) + D12 does not work; simulation of QCD effects in hadronic f.s.
      lst(9)  =  5   ! D5 scale in hadronic showers
      lst(15) = 10   ! (10) parton distribution functions, from CTEQ2D.
      lst(16) = 1    ! + use pdflib rather than internal pythia selection of pdf's 
   


c     FERMI :
      lst(17) = 1     ! D0; energy allow to vary(1). vary energies when using Fermi
  
      lst(19) = -1    !D-10 choice of the grid
      lst(23) = 1    ! which process simulated, e.g. electromagnetic
      lst(25) = 1    ! flavor of struck quark, 1=d, 2=u 
      lst(26) = 5    ! entry line in event record of struck quark
      lst(27) = 1    ! split of non-trivial nucleon remnant
      lst(31) = 1    ! vairiables of integration are (x,Q2)

c     EXTRA Modifications (Saclay)
c     parl(3) = 0.44
      parl(8) = 0.02 ! cut-off against divergencies D=0.04 (for D lst(20)=5)
      parl(9) = 3.5  ! cut-off against divergencies D = 4  (for D lst(20)=5)  

c     JetSet LUDAT1 switch
      parj(21)= 0.3  !width of the primary gauss pT in fragmentation (D = 0.35)

c     Pythia PYDAT1 : COMMON/PYDAT1/MSTU(200), PARU(200),MSTJ(200),PARJ(200)
      parj(23) = 0.02   !0.02 !tail of the small gauss for the pT frag (D = 0.01) 
      parj(24) = 3.5   !3.5  !width of the small gauss (D=2.0)
      parj(41) = 0.6  !param 'a' of sym Lund fragmentation fnc
      parj(42) = 0.1  !param 'b' -- with D choice of frag fnc mstj(11) = 4
   
c     These are for the general survey plots
      cut(1)  = 0.01   ! 0.02 ! lower limit Bj x
      cut(2)  = 0.9    ! 0.90 !upper limit Bj x
      cut(3)  = 0.     ! lower limit y
      cut(4)  = 1.     ! upper limit y
      cut(5)  = 0.8    ! 0.1 ! 0.8    ! lower limit Q2
      cut(6)  = plz    ! 10upper limit Q2 
c      cut(6)  = 22.    ! 10upper limit Q2 
      cut(7)  = 4.     ! 4 ! lower limit W2
      cut(8)  = 20.    ! upper limit W2 !lepto.f l.277 
c      cut(8)  = 22.    ! upper limit W2 !lepto.f l.277 
      cut(9)  = 0.     ! lower limit nu
      cut(10) = plz - 0.01    ! upper limit nu !previously = beam_e - .3
      cut(11) = 0.01   ! 0.03 ! 0.3    ! lower limit E'
      cut(12) = plz    ! upper limit E' in given ref frame
c     cut(13) =  0.08  ! lower limit electron scat. angle, rad. (0.08 ~ 5 degrees)
c      cut(14) =  1.05  ! upper limit electron scat. angle, rad. (1.05 ~ 60 degrees)


c      parj(33) = 0.00001; !string tension
c      parj(34) = 0.0001;

c     Hayk's
c      parj(41)= 1.13
c      parj(42)= 0.37
c      mstp(82)= 1



c     Pythia PYPARS : COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)     
c     lst(15) =  0  !parton densities controlled through Pythia MSTP 
c     mstp(52)= 2    ! analog of lst(16)   
c     mstp(51)= 5012  !5012 !choice of pdf, GRV98LO lst(15)
c     pdf = 5012    


c      parl(14) = 0.35    !+ D0.35 rT  width of G dist in pT, i.e width in the remnant split (lst27)  
c      parl(20) = 0.1     !+ D0.1 mass reduction for the complex remnant 
c      parl(26) = 0.22   !+ D?! value of Lambda QCD from parl(25)
c      parl(18) = 0.044   !+ D0.044 delta r from rad.corrections


c     TARGET

      parl(1) =  tarA
      parl(2) =  tarZ

      ppz = 0 
      CALL TIMEX(T1)
      print*,'Calling LINIT'
      CALL LINIT(0,11,plz,ppz,lst(23))
      print*,'Calling compleeteeeee'
      CALL TIMEX(T2) 
      call itime(now)
      call idate(today)
      print*,'Call to LINIT is complete'
      print*,'Default random seed : ',MRLU(1)
      print*, rand(0),today(1),today(2),today(3),now(1),now(2),now(3)
      if(today(1).eq.0)today(1)=1
      if(today(2).eq.0)today(2)=1
      if(today(3).eq.0)today(3)=1
      if(now(1).eq.0)now(1)=1
      if(now(2).eq.0)now(2)=1
      if(now(3).eq.0)now(3)=1
c     DO NOT USE THE YEAR IN THE PRODUCT !

      MRLU(1)=today(1)*today(2)*now(1)*now(2)*now(3)      

      print*,'>>>>>>>>>>>> RANDOM SEED USED FOR THIS RUN : ',MRLU(1)
      print*,'>>>>>>>>>>>> TARGET CHOOSEN A = ', PARL(1)
              
      WRITE(*, 950) "#" ,"NEvents", "Xb", "y", "W", "Q2", "Nu"
 950  FORMAT(A8,X,A8,X,A8,X,A8,X,A8,X,A8,X,A8)
      
      DO 500 NE=1,NEVENT

c     Checking the loop
c      print*,"processing event: ",ne

  100 CALL LEPTO

      IF(LST(21).NE.0) THEN
        WRITE(6,*) ' Warning: LST(21)=',LST(21),' event skipped'
        GOTO 100
      ENDIF

      IF(MOD(NE,1).EQ.0) THEN
        WRITE(1, 951) N, X, Y, SQRT(W2), Q2, U
 951    FORMAT(I4,X ,F10.7,X ,F10.7,X ,F10.7,X ,F10.7,X ,F10.7)
        CALL LULIST(1)
      ENDIF
  500 CONTINUE
      CALL TIMEX(T3)
      WRITE(6,1000) PARL(23),PARL(24),NEVENT,
     &              T2-T1,T3-T2,(T3-T2)/MAX(1,NEVENT)

 1000 FORMAT(/,' Cross section from numerical integration:',G12.3,' pb',
     &       /,' Cross section from MonteCarlo simulation:',G12.3,' pb',
     &       /,' Number of events generated:              ',I12,
     &       /,' Time for initialization (LINIT):         ',F12.3,' s',
     &       /,' Time for event generation,     total:    ',F12.3,' s',
     &       /,'                            per event:    ',F12.3,' s')
      END

      logical function check_pid(pid)
      integer pid
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      DOUBLE PRECISION P,V
      INTEGER N,K,Neta,Npid
      check_pid=.FALSE.
      Neta = 0
      DO 200 Npid=1,N
      IF(K(NPID,2).EQ.pid .AND. K(NPID,3).NE.0) THEN ! check if particle selected is in final state
         check_pid=.TRUE.
      ENDIF
 200  CONTINUE

      return
      end
