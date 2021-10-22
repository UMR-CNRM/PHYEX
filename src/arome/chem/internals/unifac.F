!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/tulet/MESONH/prepsource/test/unifac.f,v $ $Revision: 1.1 $ $Date: 2005/06/02 10:27:34 $
!-----------------------------------------------------------------

      SUBROUTINE UNIFAC (NMOL,NFUNC,NU,X,A,RG,QG,Z,TEMP,                &
     & GAMA)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK

!********************************************************
! Function: unifac
!           SUBROUTINE TO ESTIMATE ACTIVITY COEFFICIENTS &
!           AW OF MULTICOMPONENT MIXTURES FROM UNIFAC
!
! Preconditions required: called by unidriver, the C
!                         driver that replaces the Fortran
!                         codes above
!
! Return values:  GAMA(I): ACTIVITY COEFFICIENT OF COMPONENT I ON
!                 MOLE FRACTION SCALE
!
!
! Notes: 1. FOR DIMENSIONAL PURPOSES, WE HAVE ASSUMED THAT
!           MAXIMUM NUMBER OF FUNCTIONAL GROUPS =50 AND
!           MAXIMUM NUMBER OF MOLECULAR ENTITIES =50.
!
!        2. INPUT SUMMARY: 1) GROUP VOLUME (RG), GROUP SURFACE
!           AREA (QG) OF EACH FUNCTIONAL GROUP IN EACH
!           COMPONENT; 2) MOLE FRACTION OF ALL COMPONENTS
!           3) INTER-FUNCTIONAL GROUP INTERACTION PARAMETERS
!           A(J1,J2) FOR ALL FUNCTIONAL GROUP COMBOS.
!           "COMPONENT" MEANS A MOLECULAR LEVEL ENTITY
!
!        3. AVOID SENDING ZERO MOLE FRACTIONS IN ARRAY X(I).  THE
!           PROGRAM CAN HANDLE X(I) = 0 (IT SETS ACTIVITY COEFF TO 1)
!           BUT BETTER TO AVOID IT. (note from P.Saxena.)
!           (BKP NOTE: activity coefficient is not set to 1 at x = 0)
!
!        4. 3-17-98 CHANGED ALOG FUNCTION CALLS TO DLOG FOR TYPE
!           COMPATABILITY
!
! Revision History: Developed by Pradeep Saxena, EPRI, 95
!                   Revised by Betty Pun, AER, Nov 99 Under EPRI/CARB
!                   funding to comply with models-3 standards
! *****************************************************************

      IMPLICIT NONE
!      IMPLICIT REAL*8 (A-H,O-Z)

!........Arguments and their description

!     INTEGER*4 NMOL [VALUE]
!     INTEGER*4 NFUNC [VALUE]

      INTEGER*4 NMOL         ! TOTAL NO. OF MOLECULAR ENTITIES IN THE MIXTURE
      INTEGER*4 NFUNC        ! TOTAL NO. OF FUNCTIONAL GROUPS IN THE MIXTURE

!     INTEGER*4 NU [REFERENCE] (50, 50)
!     REAL*8 X [REFERENCE] (50)
!     REAL*8 A [REFERENCE] (50, 50)
!     REAL*8 RG [REFERENCE] (50)
!     REAL*8 QG [REFERENCE] (50)
!     REAL*8 Z [VALUE]
!     REAL*8 TEMP [VALUE]
!     REAL*8 GAMA [REFERENCE] (50)

      INTEGER*4 NU(NMOL,NFUNC)    ! NU(I,J): VECTOR OF NO. OF A PARTICULAR
                             ! FUNCTIONAL GROUP J IN MOLECULE I
                             ! (FROM STOICHIOMETRY)
      REAL*8 X(NMOl)           ! MOLE FRACTION (or amount) OF MOLECULE I
      REAL*8 A(NFUNC,NFUNC)        ! INTERACTION PARAMETER FOR GROUPS J1 & J2
      REAL*8 RG(NFUNC)          ! VAN DER VAAL VOLUME FOR GROUP J
      REAL*8 QG(NFUNC)          ! VAN DER VAAL SURFACE AREA FOR GROUP J
      REAL*8 GAMA(NMOL)        ! calculated activity coefficients


      REAL*8 Z               ! COORDINATION NUMBER FOR THE SOLVENT = 10
      REAL*8 TEMP            ! TEMPERATURE IN DEGREES KELVIN

! ....  Parameters and their descriptions
!      DIMENSION R(50),Q(50),RL(50),RX(50),QX(50),XL(50),PHI(50)
!      DIMENSION THETA(50),XGM(50),XGP(50,50),THTAGP(50,50),
!     $     THTAGM(50),TTSIM(50),TTSIP(50),SI(50,50),GAMMLN(50),
!     $     GAMPLN(50,50),YCLN(50),YRLN(50)

      REAL*8 R(50)           ! total volume (R) for each molecule
      REAL*8 Q(50)           ! total surface area (Q) for each molecule
      REAL*8 RL(50)
      REAL*8 RX(50)          ! R(I) * X(I)
      REAL*8 QX(50)          ! Q(I) * X(I)
      REAL*8 XL(50)          ! RL(I) * X(I)
      REAL*8 PHI(50)
      REAL*8 THETA(50)
      REAL*8 XGM(50)
      REAL*8 XGP(50,50)
      REAL*8 THTAGP(50,50)
      REAL*8 THTAGM(50)
      REAL*8 TTSIM(50)
      REAL*8 TTSIP(50)
      REAL*8 SI(50,50)
      REAL*8 GAMMLN(50)
      REAL*8 GAMPLN(50,50)
      REAL*8 YCLN(50)
      REAL*8 YRLN(50)

      REAL*8 ZHALF           ! 0.5 * Z
      REAL*8 ZOT1
      REAL*8 ZOT2
      REAL*8 ZOT4
      REAL*8 ZOT5
      REAL*8 ZOT6
      REAL*8 ZOT7
      REAL*8 ZOT8
      REAL*8 ZOT9
      REAL*8 ZOT10
      REAL*8 ZOT11
      REAL*8 SUMXGP
      REAL*8 SUMXGM
      REAL*8 SUMXL          ! SUM(RL(I) * X(I))
      REAL*8 SUMJX
      REAL*8 CTEMP
      REAL*8 SUMTGM
      REAL*8 SUMTGP
      REAL*8 SUMQX          ! SUM (Q(I) * X(I))
      REAL*8 SUMRX          ! SUM (R(I) * X(I))

      INTEGER I, J, J1, J2  ! loop counters


!    Subroutine writes to error file, need to open
!swu      OPEN(7,FILE='unifac.err')

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('UNIFAC',0,ZHOOK_HANDLE)
      ZHALF = 0.5 * Z

!     **************************************************************
!     COMBINATORIAL PART--THIS IS FOR ORIGINAL UNIFAC
!     NEED TO REVISE IT FOR MODIFIED UNIFAC
!     **************************************************************
!     COMPUTE R AND Q FOR EACH MOLECULE
      DO 20 I = 1,NMOL
      R(I) = 0.0
      Q(I) = 0.0
!
      DO 10 J = 1, NFUNC
      R(I) = R(I) + FLOAT(NU(I,J)) * RG(J)
      Q(I) = Q(I) + FLOAT(NU(I,J)) * QG(J)
10    CONTINUE
!
20    CONTINUE
!     WRITE(6,9900)Q,R
9900  FORMAT('Q AND R',4F12.5)
!
!
      DO 30 I = 1,NMOL
      RL(I) = ZHALF * (R(I) - Q(I)) - R(I) + 1.
      RX(I) = R(I) * X(I)
      QX(I) = Q(I) * X(I)
      XL(I) = X(I) * RL(I)
30    CONTINUE
!
      SUMRX = 0.
      SUMQX = 0.
      SUMXL = 0.
      DO 50 I = 1,NMOL
      SUMRX = SUMRX + RX(I)
      SUMQX = SUMQX + QX(I)
      SUMXL = SUMXL + XL(I)
50    CONTINUE
!
      ZOT1 = 0.0
      IF (SUMRX .NE. 0.) THEN
      ZOT1 = 1./SUMRX
      ELSE
      WRITE (7,5000)
5000  FORMAT('ZERO DIVISION IN COMBINATORIAL--SUMRX=0.')
      STOP
      END IF
      ZOT2 = 0.0
      IF (SUMQX .NE. 0.) THEN
      ZOT2 = 1./ SUMQX
      ELSE
      WRITE(7,5010)
5010  FORMAT('ZERO DIVISION IN COMBINATORIAL-SUMQX=0.')
      STOP
      END IF
!
      DO 60 I = 1,NMOL
      PHI(I) = RX(I) * ZOT1
      THETA(I) = QX(I) * ZOT2
60    CONTINUE
!
!     SOLVE FOR THE COMBINATORIAL LN(GAMAC)
      DO 70 I = 1,NMOL
      YCLN(I) = 0.0
      IF (X(I) .EQ. 0.0) THEN
      GO TO 70
      ELSE
      END IF
      IF((PHI(I)) .NE. 0. .AND. (THETA(I)) .NE. 0.) THEN
      ZOT4 = PHI(I)/X(I)
      ZOT5 = THETA(I)/PHI(I)
      YCLN(I) = LOG(ZOT4) + ZHALF * Q(I) * LOG(ZOT5) + RL(I) -          &
     & ZOT4 * SUMXL
      ELSE
      WRITE(7,5015)
5015  FORMAT('THETA(I) ,OR. PHI(I) .EQ. 0.0')
      STOP
      END IF
70    CONTINUE
!     WRITE(6,9991)YCLN(1),YCLN(2)
9991  FORMAT('YCLN',2E12.5)
!
!     *****************************************************************
!     NOW COMPUTE THE RESIDUAL PART--THIS IS FOR ORIGINAL
!     UNIFAC.  NEED TO REVISE IT FOR MODIFIED UNIFAC
!     *****************************************************************
!
!       COMPUTE GROUP MOLE FRACTIONS
!
!       COMPUTE TOTAL GROUP MOLES AND GROUP MOLE FRACTIONS
!       IN PURE COMPOUNDS AND MIXTURE
!
      SUMXGM = 0.0
      DO 150 I = 1, NMOL
!
      SUMXGP = 0.0
      DO 110 J = 1, NFUNC
      SUMXGM = SUMXGM + FLOAT(NU(I,J)) * X(I)
      SUMXGP =SUMXGP + FLOAT(NU(I,J))
110   CONTINUE
!
!     NOW LOOP THRU TO GET MOLE FRACTIONS FOR PURE C COMPOUNDS
      IF (SUMXGP .NE. 0.) THEN
      ZOT6 = 1./ SUMXGP
      ELSE
      WRITE(7,5020)
5020  FORMAT('ZERO DIVISION IN RESIDUAL PART-SUMXGP=0.')
      STOP
      END IF
      DO 120 J = 1, NFUNC
      XGP(I,J) = FLOAT(NU(I,J)) * ZOT6
!     WRITE(6,7770)I,J,XGP(I,J)
7770  FORMAT('XGP',2I10,2E12.5)
120   CONTINUE
!
!
150   CONTINUE
!
!
!     NOW WE HAVE TOTAL MOLES OF ALL GROUPS IN THE MIX
!     LOOP ONCE THRU TO COMPUTE GROUP MOLE FRACTION
!     IN THE MIXTURE
      IF (SUMXGM .LE. 0.) THEN
      WRITE(7,5030)
5030  FORMAT('ZERO DIVISION IN RESIDUAL PART-SUMXGM=0.')
      STOP
      ELSE
      ZOT7 = 1./SUMXGM
      END IF
      DO 180 J = 1,NFUNC
      SUMJX = 0.0
      DO 170 I=1,NMOL
      SUMJX = FLOAT(NU(I,J)) * X(I) + SUMJX
170   CONTINUE
      XGM(J) = SUMJX * ZOT7
!     WRITE(6,7771)J,XGM(J)
7771  FORMAT('XGM',I10,E12.5)
180   CONTINUE
!
!     COMPUTE TOTAL GROUP AREA & GROUP AREA FRACTIONS
!     IN PURE COMPOUNDS
!
      DO 250 I = 1, NMOL
      SUMTGP = 0.0
!
      DO 210 J = 1, NFUNC
      SUMTGP =SUMTGP + XGP(I,J) * QG(J)
210   CONTINUE
!
!     NOW LOOP THRU TO GET AREA FRACTIONS FOR PURE
!     COMPOUNDS
      IF (SUMTGP .NE. 0.) THEN
      ZOT8 = 1./ SUMTGP
      ELSE
      WRITE(7,5040)
5040  FORMAT('ZERO DIVISION IN RESIDUAL PART-SUMTGP=0.')
      STOP
      END IF
!
      DO 220 J = 1, NFUNC
      THTAGP(I,J) = QG(J) * XGP(I,J)* ZOT8
!     WRITE(6,7772)I,J,THTAGP(I,J)
7772  FORMAT('THTAGP',2I10,E12.5)
220   CONTINUE
!
!
250   CONTINUE
!
!     COMPUTE TOTAL GROUP AREA IN MIXTURE
      SUMTGM = 0.0
      DO 260 J = 1, NFUNC
      SUMTGM = SUMTGM + XGM(J) * QG(J)
260   CONTINUE
!
!     NOW WE HAVE TOTAL AREA OF ALL GROUPS IN THE MIX
!     LOOP ONCE THRU TO COMPUTE GROUP AREA FRACTION
!     IN THE MIXTURE
      IF (SUMTGM .LE. 0.) THEN
      WRITE(7,5050)
5050  FORMAT('ZERO DIVISION IN RESIDUAL PART-SUMTGM=0.')
      STOP
      ELSE
      ZOT9 = 1./SUMTGM
      END IF
      DO 280 J = 1,NFUNC
      THTAGM(J) = QG(J) * XGM(J) * ZOT9
!     WRITE(6,7773)J,THTAGM(J)
7773  FORMAT('THTAGM',I10,E12.5)
280   CONTINUE
!
!
!
!     COMPUTE SI VALUES FROM A(J1,J2) AND TEMPERATURE
      IF (TEMP .NE. 0.)THEN
      CTEMP = 1./TEMP
      ELSE
      WRITE(7,5060)
5060  FORMAT('ZERO DIVIDE IN RESIDUAL-TEMP=0.0')
      STOP
      END IF
      DO 300 J1 = 1,NFUNC
      DO 300 J2 = 1,NFUNC
      SI(J1,J2) = EXP (-A(J1,J2)*CTEMP)
!     WRITE(6,7774)J1,J2,SI(J1,J2)
7774  FORMAT('SI',2I10,E12.5)
300   CONTINUE
!
!
!
!     NOW COMPUTE CAPITAL GAMMAS FOR PURE COMPOUNDS
!     AND MIXTURES
!
!     CAPITAL GAMAS FOR FUNCTIONAL GROUPS IN MIXTURES
!     ARE INDEPENDENT OF MOLECULAR ENTITIES
!
      DO 500 J1 = 1, NFUNC
      TTSIM(J1) = 0.0
      DO 400 J2 = 1,NFUNC
      TTSIM(J1) = TTSIM(J1) + THTAGM(J2) * SI(J2,J1)
400   CONTINUE
500   CONTINUE
!
!
      DO 550 J1 = 1,NFUNC
      ZOT10 = 0.0
      DO 520 J2 = 1,NFUNC
      IF (TTSIM(J2) .NE. 0.0) THEN
      ZOT10 = ZOT10 + THTAGM(J2) * SI(J1,J2)/TTSIM(J2)
      ELSE
      WRITE(7,5070)
5070  FORMAT('ZERO DIVIDE IN RESIDUAL TTSIM = 0.0')
      STOP
      END IF
520   CONTINUE
      IF (TTSIM(J1) .NE. 0.0) THEN
      GAMMLN(J1) = QG(J1) * (1. - LOG(TTSIM(J1)) - ZOT10)

!     WRITE(6,7775)J1,GAMMLN(J1)
7775  FORMAT('GAMMLN',I10,E12.5)
      ELSE
      WRITE(7,5080)
5080  FORMAT('ZERO ARGUMENT OF LOG IN RES TTSIM')
      STOP
      END IF
550   CONTINUE
!
!
!     CAPITAL GAMAS FOR FUNCTIONAL GROUP IN PURE
!     COMPOUNDS ARE DEPENDENT ON MOLECULE--SO
!     WE LOOP THRU EACH MOLECULE-FUNCTIONAL GROUP
!     PAIR
!
      DO 900 I = 1, NMOL
!
!
      DO 700 J1 = 1, NFUNC
      TTSIP(J1) = 0.0
      DO 600 J2 = 1,NFUNC
      TTSIP(J1) = TTSIP(J1) + THTAGP(I,J2) * SI(J2,J1)
600   CONTINUE
700   CONTINUE
!
!
      DO 750 J1 = 1,NFUNC
      ZOT11 = 0.0
      DO 720 J2 = 1,NFUNC
      IF (TTSIP(J2) .NE. 0.0) THEN
      ZOT11 = ZOT11 + THTAGP(I,J2) * SI(J1,J2)/TTSIP(J2)
      ELSE
      WRITE(7,5090)
5090  FORMAT('ZERO DIVIDE IN RESIDUAL TTSIP = 0.0')
      STOP
      END IF
720   CONTINUE
      IF (TTSIP(J1) .NE. 0.0) THEN
      GAMPLN(I,J1) = QG(J1) * (1. - LOG(TTSIP(J1)) - ZOT11)
!     WRITE(6,7776)I,J1,GAMPLN(I,J1)
7776  FORMAT('GAMPLN',2I10,E12.5)
      ELSE
      WRITE(7,6000)
6000  FORMAT('ZERO ARGUMENT OF LOG IN RES TTSIP')
      STOP
      END IF
750   CONTINUE
!
!
900   CONTINUE
!
!     COMPUTE RESIDUAL PART AND ADD TO
!     COMBINATORIAL PART
!
      DO 1100 I = 1, NMOL
!
      YRLN(I) = 0.0
      DO 1000 J = 1,NFUNC
      YRLN(I) = YRLN(I) + FLOAT(NU(I,J)) * (GAMMLN(J) - GAMPLN(I,J))
1000  CONTINUE
!
      GAMA(I) = EXP ( YCLN(I) + YRLN(I))
!
1100  CONTINUE
!
      IF (LHOOK) CALL DR_HOOK('UNIFAC',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE UNIFAC
