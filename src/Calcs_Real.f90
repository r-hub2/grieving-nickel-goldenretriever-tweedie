MODULE Calcs_Real

  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  PUBLIC :: evaluateRek, evaluateRekd
  
CONTAINS

  SUBROUTINE evaluateRek(i, t, Rek)
    ! Find the value of Re k(t)
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)    :: t
    INTEGER(C_INT), INTENT(IN)         :: i
    REAL(KIND=C_DOUBLE), INTENT(OUT)   :: Rek
  
    REAL(KIND=C_DOUBLE) :: current_mu, current_phi
    REAL(KIND=C_DOUBLE) :: pi, omega, pindex, front, alpha, tanArg
  
    
    ! Grab the relevant scalar values for this iteration:
    current_mu   = Cmu(i)   ! Access mu value for index i
    current_phi  = Cphi(i)  ! Access phi value for index i
    
  
    pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
    pindex = (2.0E0_C_DOUBLE - Cp)
    front = current_mu ** pindex  / ( current_phi * pindex)
    tanArg = (1.0E0_C_DOUBLE - Cp) * t * current_phi / (current_mu ** (1.0E0_C_DOUBLE - Cp) )
    omega = DATAN( tanArg )

    ! Safety check
    IF ((omega .GT. 0.0E0_C_DOUBLE ) .OR. (omega .LT. (-pi/2.0E0_C_DOUBLE))) THEN
       ! Error!
        CALL DBLEPR("ERROR (evaluateRek): Omega out of range:", -1, omega, 1)
        CALL DBLEPR("ERROR (evaluateRek): t:", -1, t, 1)
        CALL DBLEPR("         p:", -1, Cp, 1)
        CALL DBLEPR("        mu:", -1, current_mu, 1)
        CALL DBLEPR("       phi:", -1, current_phi, 1)
       RETURN
    END IF
    
    alpha = (2.0E0_C_DOUBLE - Cp)/(1.0E0_C_DOUBLE - Cp)
    Rek = front * ( DCOS(omega * alpha)/(DCOS(omega)**alpha) - 1.0E0_C_DOUBLE )

  END SUBROUTINE evaluateRek



  SUBROUTINE evaluateRekd(i, t, Redk)
    ! Find the value of Re k'(t)
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)    :: t
    INTEGER(C_INT), INTENT(IN)         :: i
    REAL(KIND=C_DOUBLE), INTENT(OUT)   :: Redk
    
    REAL(KIND=C_DOUBLE)       :: omega, pindex
    REAL(KIND=C_DOUBLE)       :: current_mu, current_phi
  
  
    ! Grab the relevant scalar values for this iteration:
    current_mu   = Cmu(i)   ! Access mu value for index i
    current_phi  = Cphi(i)  ! Access phi value for index i
  
    
    pindex = 1.0E0_C_DOUBLE / (1.0E0_C_DOUBLE - Cp)
    omega = DATAN( ( (1.0E0_C_DOUBLE - Cp) * t * current_phi) / &
                   (current_mu ** (1.0E0_C_DOUBLE - Cp) ) )
    
    Redk = current_mu * &
           DSIN( omega * pindex ) / &
           (DCOS(omega)**pindex)
    
  END SUBROUTINE evaluateRekd
  
  
  
  
    
  SUBROUTINE evaluateLambda(i, lambda)
    ! Find lambda, such that P(Y = 0) = exp( -lambda ) when 1 < p < 2 
    
    USE tweedie_params_mod, ONLY: Cmu, Cphi, Cp, CpSmall
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
    IMPLICIT NONE
  
    INTEGER(C_INT), INTENT(IN)        :: i
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: lambda 
    
    REAL(KIND=C_DOUBLE)               :: current_mu, current_phi
  
  
    ! Grab the relevant scalar values for this iteration:
    current_mu   = Cmu(i)   ! Access mu value for index i
    current_phi  = Cphi(i)  ! Access phi value for index i
    
  
    lambda = 0.0E0_C_DOUBLE
    IF (CpSmall) THEN
      ! The calculation for lambda (used in P(Y=0) = exp(-lambda))
      lambda = (current_mu ** (2.0E0_C_DOUBLE - Cp) ) / &
               (current_phi * (2.0E0_C_DOUBLE - Cp) )
      ! NOTE: No negative sign in front
    END IF
    
  END SUBROUTINE evaluateLambda
     
      
      
      

END MODULE Calcs_Real