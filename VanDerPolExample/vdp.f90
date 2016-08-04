!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!   frc :      A periodically forced system
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
        DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
        DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

        DOUBLE PRECISION :: MU,T
        INTEGER :: I

        MU  = PAR(1) 
        T   = PAR(11) 

        F(1) = U(2)
        F(2) = MU*U(2)*(1.0d0-U(1)**2.0d0)-U(1)

        DO I=1,NDIM
          F(I) = T*F(I)
        ENDDO

        ! PRINT*, OMEGA, A, I0, DELTA, J
        !PRINT*, P, F(3)
        ! PRINT*, F
        ! STOP

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)  
!     ---------- -----

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
        DOUBLE PRECISION, INTENT(IN) :: T

        DOUBLE PRECISION MU

        MU  = 10.0d0

        PAR(1)  = MU
        PAR(11) = 1.9072769e+01

        ! U(1)=
        ! U(2)=
        ! U(3)=
        ! U(4)=

      END SUBROUTINE STPNT

      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
      !--------- ---- 
      
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
        DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
        DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
        DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

        INTEGER :: I
      
        DO I=1,NDIM
          FB(I) = U1(I)- U0(I)
        ENDDO
      
      END SUBROUTINE BCND

      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
      !--------- ----

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
        DOUBLE PRECISION, INTENT(IN) :: PAR(*)
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
        DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
        DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

        INTEGER :: I
        
        FI(1) = 0.0d0
        DO I=1,NDIM
          FI(1) = FI(1) + UPOLD(I)*( U(I) - UOLD(I) )
        ENDDO

      END SUBROUTINE ICND



      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
