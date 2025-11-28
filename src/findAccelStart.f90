SUBROUTINE findAccelStart(i, tmax, tStartAcc) 
  ! Find the value of t for when to invoke the acceleration algorithm 
  
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE), INTENT(IN)   :: tmax     
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: tStartAcc     ! The output starting point for acceleration
  INTEGER(C_INT), INTENT(IN)        :: i

  REAL(KIND=C_DOUBLE)               :: pi, t_nmax
  INTEGER(C_INT)                    :: n_max
  REAL(KIND=C_DOUBLE)               :: current_y, current_mu, current_phi
  
  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
  
  ! Start acceleration once the largest t of the region to integrate
  ! (a) exceeds tmax;
  ! (b) exceeds the final TP (inflection?) of Im k(t); AND
  ! (c) exceeds the final TP (inflection?) of Re k(t).
  !
  ! The number of infltecions points of Im k(t) is the same, or one larger, tha for Re k(t).
  ! So find the vale of t correspinding toi the number of inflections in Im k(t) only,
  ! and tRightMost is the maximum of t_max and this inflection value.

  ! Check the Im k(t) criterion
  n_max = FLOOR( Cp / (2.0_C_DOUBLE * (Cp - 1.0_C_DOUBLE) ) ) ! The maximum number of TPs
  t_nmax = current_mu ** (1.0_C_DOUBLE - Cp) / ( current_phi * (1.0_C_DOUBLE - Cp)) *  &
           DTAN( 1.0_C_DOUBLE - Cp * n_max * pi / Cp )
  ! This is the value of t corresponding to the largest inflection point in Im k(t),
  ! which is never smaller than the largest inflection pint of Re k(t).
  
  tStartAcc = MAX(tmax, t_nmax)

  
END SUBROUTINE findAccelStart
