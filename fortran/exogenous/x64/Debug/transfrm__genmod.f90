        !COMPILER-GENERATED INTERFACE MODULE: Fri Oct  2 13:23:05 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRANSFRM__genmod
          INTERFACE 
            SUBROUTINE TRANSFRM(DIM,DTR,X,TX)
              INTEGER(KIND=4), INTENT(IN) :: DIM
              INTEGER(KIND=4), INTENT(IN) :: DTR
              REAL(KIND=8), INTENT(IN) :: X(DIM)
              REAL(KIND=8), INTENT(OUT) :: TX(DIM)
            END SUBROUTINE TRANSFRM
          END INTERFACE 
        END MODULE TRANSFRM__genmod
