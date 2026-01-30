
MODULE tweedie_params_mod
  ! Save commonly-used parameters, for access throughout functions and subroutines
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE , C_BOOL
  
  IMPLICIT NONE
  SAVE
  
  ! Global Parameters
  REAL(KIND=C_DOUBLE), ALLOCATABLE  :: Cmu(:), Cphi(:), Cy(:) 
  REAL(KIND=C_DOUBLE)               :: Cp
  LOGICAL(C_BOOL)                   :: CpSmall, Cverbose, Cpdf
  INTEGER(C_INT)                    :: CN
  
END MODULE tweedie_params_mod

