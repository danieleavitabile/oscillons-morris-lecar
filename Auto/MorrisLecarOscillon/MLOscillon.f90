!----------------------------------------------------------------------
!   ML Bursting Network: Synaptically Connected Network of Moris-Lecar Neurons
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION ::IAPP,RHO,GCA,SPEED,BETA,VT,KAPPA,A,B,T
      DOUBLE PRECISION ::V2,V3,V4,GKCA,KCA,EPSI,MU
      INTEGER :: I,J,IV,IIN,IC,IS,ROW,COL,JV,JN,JC,JS

      INTEGER, PARAMETER :: NX = 460
      INTEGER, PARAMETER :: MAX_ROWS = NX
      DOUBLE PRECISION, PARAMETER :: DX = 60.0d0/(NX-1)
      DOUBLE PRECISION MINF,NINF,FR,ICA,HILL,R,DMINF,DNINF,DR,DFR,DHILL,V,N,CA
      DOUBLE PRECISION DFV,DFN,DFCA,DGV,DGN,DGCA,DHV,DHN,DHCA

      DOUBLE PRECISION, DIMENSION(NX,NX) :: W
      COMMON W

      !AUX FUNCTIONS
      MINF(V)         = (1.0+TANH((V+1.2)/18.0))/2.0
      NINF(V)         = (1.0+TANH((V-12.0)/30.0))/2.0
      R(V,RHO)        = RHO*COSH((V-12.0)/(2.0*30.0))
      DMINF(V,V2)     = (COSH((5.0*V+6.0)/90.0)**2)/(9.0*(1.0+COSH((5.0*V+6.0)/45.0)**2))
      DNINF(V,V3,V4)  = (COSH((V-12.0)/30.0)**2)/(15.0*(1.0+COSH((V-12.0)/15.0)**2))
      DR(V,V3,V4,RHO) = (RHO/60.0)*SINH((V-12.0)/60.0)
      FR(V,KAPPA,VT)  = 1.0/(1.0 + EXP(-KAPPA*(V-VT)))
      DFR(V,KAPPA,VT) = KAPPA*FR(V,KAPPA,VT)*(1.0-FR(V,KAPPA,VT))
      HILL(CA)        = CA/(CA+10.0)
      DHILL(CA)       = 10.0/(CA+10.0)**2

      !JACOBIAN COMPONENT FUNCTIONS
      DFV(V,N,CA,GKCA,GCA,V2)    = (-8.0*N + GKCA*HILL(CA) - 2.0 -GCA*(DMINF(V,V2)*(V-120.0) + MINF(V)))/20.0
      DFN(V,N,CA)                = -8.0*(V+84.0)/20.0
      DFCA(V,N,CA,GKCA)          = -GKCA*DHILL(CA)*(V+84.0)

      DGV(V,N,CA,V3,V4,RHO)      = DR(V,V3,V4,RHO)*(NINF(V)-N) + R(V,RHO)*DNINF(V,V3,V4)
      DGN(V,N,CA,RHO)            = -R(V,RHO)
      DGCA(V,N,CA)               = 0.0

      DHV(V,N,CA,EPSI,MU,GCA,V2) = -EPSI*MU*GCA*(DMINF(V,V2)*(V-120.0)+MINF(V))
      DHN(V,N,CA)                = 0.0
      DHCA(V,N,CA,EPSI)          = -EPSI


      !ASSIGN PARAMETERS
      IAPP  = PAR(1)
      RHO   = PAR(2)
      GCa   = PAR(3)
      SPEED = PAR(4)
      BETA  = PAR(5)
      VT    = PAR(6)
      KAPPA = PAR(7)
      A     = PAR(8)
      B     = PAR(9)
      V4    = PAR(10)
      T     = PAR(11)
      GKCA  = PAR(12)
      KCA   = PAR(13)
      EPSI  = PAR(14)
      MU    = PAR(15)
      V2    = PAR(16)
      V3    = PAR(17)

      !Right-hand side
      PRINT *, "Evaluating RHS... with IJAC=", IJAC
      DO  I=1,NX

          IV  = I
          IIN = NX + I
          IC  = 2*NX + I
          IS  = 3*NX + I

          F(IV)  = U(IS) + (-(8.0*U(IIN) + GKCA*HILL(U(IC)))*(U(IV)+84.0)-2.0*(U(IV)+60.0) &
                   -GCA*MINF(U(IV))*(U(IV)-120.0) + IAPP)/20.0
          F(IIN) = R(U(IV),RHO)*(NINF(U(IV))-U(IIN))
          F(IC)  = EPSI*(-MU*(GCA*MINF(U(IV))*(U(IV)-120.0))-U(IC))
          F(IS)  = -BETA*U(IS)

           DO J=1,NX
              F(IS) = F(IS) + W(I,J)*FR(U(J),KAPPA,VT)*DX
           END DO

      ENDDO
      PRINT *, "DONE"

      PRINT *, "Evaluating RHS..."
      !Jacobian
      IF (IJAC .EQ. 1) THEN

        ! Initialise to 0
        DFDU = 0.0d0

        DO  I=1,NX

          IV  = I
          IIN = NX + I
          IC  = 2*NX + I
          IS  = 3*NX + I

          DFDU(IV,IV)   = DFV(U(IV),U(IIN),U(IC),GKCA,GCA,V2)
          DFDU(IV,IIN)  = DFN(U(IV),U(IIN),U(IC))
          DFDU(IV,IC)   = DFCA(U(IV),U(IIN),U(IC),GKCA)
          DFDU(IV,IS)   = 1.0

          DFDU(IIN,IV)  = DGV(U(IV),U(IIN),U(IC),V3,V4,RHO)
          DFDU(IIN,IIN) = DGN(U(IV),U(IIN),U(IC),RHO)
          DFDU(IIN,IC)  = 0.0
          DFDU(IIN,IS)  = 0.0

          DFDU(IC,IV)   = DHV(U(IV),U(IIN),U(IC),EPSI,MU,GCA,V2)
          DFDU(IC,IIN)  = 0.0
          DFDU(IC,IC)   = -EPSI
          DFDU(IC,IS)   = 0.0

          DFDU(IS,IS)   = -BETA

          DO  J=1,NX

            DFDU(IS,J)   = DX*W(IV,J)*DFR(U(J),KAPPA,VT)

          END DO

        END DO

      END IF
      PRINT *, "DONE"

      !RESCALE
      F = T*F
      IF (IJAC .EQ. 1) THEN
        DFDU = T*DFDU
      ENDIF

!     -----------------------------------------------------------------
      END SUBROUTINE FUNC


      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ----------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBlE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION IAPP,GCA,SPEED,BETA,RHO,VT,KAPPA,A,B
      DOUBLE PRECISION V2,V3,V4,GKCA,KCA,EPSI,MU

      INTEGER, PARAMETER :: NX = 460
      INTEGER :: I,J

      DOUBLE PRECISION, DIMENSION(NX,NX) :: W
      COMMON W

      IAPP   = 150
      RHO   = 0.13
      GCa   = 4.0
      SPEED = 0.175
      BETA  = 0.50
      VT    = 15.0
      KAPPA = 10.0
      A     = 2.0
      B     = 1.0
      V4    = 30.0
      GKCA  = 0.25
      KCA   = 1.0
      EPSI  = 0.003
      MU    = 0.2
      V2    = 18.0
      V3    = 12.0

      PAR(1)  = IAPP
      PAR(2)  = RHO
      PAR(3)  = GCa
      PAR(4)  = SPEED
      PAR(5)  = BETA
      PAR(6)  = VT
      PAR(7)  = KAPPA
      PAR(8)  = A
      PAR(9)  =  B
      PAR(10) = V4
      PAR(11) = 3.6500000e+01
      PAR(12) = GKCA
      PAR(13) = KCA
      PAR(14) = EPSI
      PAR(15) = MU
      PAR(16) = V2
      PAR(17) = V3

      OPEN(UNIT=1,FILE="Weights.dat")
      PRINT *, "Reading weights..."
      DO I = 1,NX
         READ(1,*) (W(I,J),J=1,NX)
      END DO
      CLOSE(1)
      PRINT *, "DONE"

      !U = 0.0d0
      ! STOP

!     ----------------------------------------------------------------
      END SUBROUTINE STPNT

      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!     ----------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
      INTEGER :: I

      DO I = 1,NDIM
        FB(I) = U0(I) - U1(I)
      ENDDO

!     ----------------------------------------------------------------
      END SUBROUTINE BCND


      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
!     ----------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
      DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

      INTEGER :: I

      FI(1) = 0
      DO I = 1,NDIM
        FI(1) = FI(1) + UPOLD(I) * ( U(I) - UOLD(I) )
      ENDDO

      ! PRINT*,FI(1)

      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
