
SUBROUTINE TweedieIntegration(i, funvalueI, exitstatus, relerr, count_Integration_Regions) 
  ! Compute the value of the integrals in the Fourier-inversion expressions for the PDF and CDF

  USE Integrands_MOD, ONLY: Integrands
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE rprintf_mod
  USE Calcs_Imag
  USE Calcs_Real
  USE Calcs_K
  
  IMPLICIT NONE
  
  INTEGER(C_INT), INTENT(IN)        :: i              ! Observation index
  INTEGER(C_INT), INTENT(OUT)       :: exitstatus     ! Output status
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: funvalueI      ! Computed result
  REAL(KIND=C_DOUBLE), INTENT(OUT)  :: relerr         ! Estimate of relative error
  INTEGER(C_INT), INTENT(OUT)       :: count_Integration_Regions ! Num int regions


  ! Local Variables: All local variables defined here
  INTEGER(C_INT)    :: mmax, mfirst, mOld, accMax
  INTEGER(C_INT)    :: count_PreAcc_Regions, count_Acc_Regions
  INTEGER(C_INT)    :: m, min_Acc_Regions
  LOGICAL(C_BOOL)   :: converged_Accelerating, converged_Pre, convergence_Acc
  REAL(KIND=C_DOUBLE)   :: kmax, tmax, aimrerr
  REAL(KIND=C_DOUBLE)   :: epsilon, areaT, pi, psi, zero
  REAL(KIND=C_DOUBLE)   :: zeroL, zeroR, area0, area1, areaA, sumA
  REAL(KIND=C_DOUBLE)   :: current_y, current_mu, current_phi
  REAL(KIND=C_DOUBLE), ALLOCATABLE   :: Mmatrix(:, :), Nmatrix(:, :), xvec(:), wvec(:)
  REAL(KIND=C_DOUBLE)   :: TMP,   leftPreAccZero, leftAccZero
  REAL(KIND=C_DOUBLE)   :: zeroStartPoint
  LOGICAL(C_BOOL)       :: leftOfMax, flip_To_Other_Side
  

  INTERFACE

    SUBROUTINE GaussQuadrature(i, a, b, area)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      REAL(KIND=C_DOUBLE), INTENT(OUT)      :: area
      REAL(KIND=C_DOUBLE), INTENT(IN)       :: a, b
      INTEGER(C_INT), INTENT(IN)            :: i
    END SUBROUTINE GaussQuadrature


    SUBROUTINE accelerate(xvec, wvec, nzeros, Mmatrix, NMatrix, West)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER, INTENT(IN)                 :: nzeros
      REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: xvec(:), wvec(:), Mmatrix(:, :), Nmatrix(:, :)
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: West
    END SUBROUTINE accelerate

  END INTERFACE


  ! Grab the relevant scalar values for this iteration:
  current_y    = Cy(i)    ! Access y value for index i
  current_mu   = Cmu(i)   ! Access mu value for index i
  current_phi  = Cphi(i)  ! Access phi value for index i
  
  ! ALLOCATE these arrays onto the HEAP
  ALLOCATE(Mmatrix(2, 501))
  ALLOCATE(Nmatrix(2, 501))
  ALLOCATE(xvec(501))
  ALLOCATE(wvec(501))
  
  IF (Cverbose) THEN
    CALL DBLEPR("*** Computing for p =", -1, Cp, 1)
    CALL DBLEPR("*** Computing for y =", -1, current_y, 1)
    CALL DBLEPR("*** Computing for mu =", -1, current_mu, 1)
    CALL DBLEPR("*** Computing for phi =", -1, current_phi, 1)
  END IF


  ! --- Initialization ---
  pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
  aimrerr = 1.0E-12_C_DOUBLE ! 1.0E-12_C_DOUBLE
  mOld = 0
  m = 0
  exitstatus = 0
  relerr = 1.0_C_DOUBLE
  convergence_Acc = .FALSE.
  epsilon = 1.0E-12_C_DOUBLE
  Mmatrix = 0.0_C_DOUBLE
  Nmatrix = 0.0_C_DOUBLE
  xvec = 0.0_C_DOUBLE
  wvec = 0.0_C_DOUBLE
  mmax = 0  
  count_Integration_Regions = 0_C_INT
  zeroStartPoint = 0.0_C_DOUBLE
  flip_To_Other_Side = .FALSE.
  zeroR = 0.0_c_DOUBLE
  zeroL = 0._C_DOUBLE

  ! --- Integration initialization ---
  area0 = 0.0_C_DOUBLE
  area1 = 0.0_C_DOUBLE
  zero  = 0.0_C_DOUBLE

  ! --- Find kmax, tmax, mmax ---
  CALL findKmax(i, kmax, tmax, mmax, mfirst, leftOfMax)
  IF (Cverbose) THEN
    CALL DBLEPR("  -            kmax:", -1, kmax, 1 )
    CALL DBLEPR("  -            tmax:", -1, tmax, 1 )
    CALL INTPR( "  -            mmax:", -1, mmax, 1 )
    CALL INTPR( "  - first zero at m:", -1, mfirst, 1 )
  END IF
  
  m = mfirst
  
  ! INTEGRATION
  ! Three integration zones:
  !   1. The *initial* area, which is not between zeros of Im{k(t)}: area0
  !   2. The initial area *before* Sidi acceleration is invoked: area1
  !      (For instance, wait until Im{k(t)} is on the downturn.)
  !   3. The area thereafter, upon which Sidi acceleration is
  !      applied; the area returned by acceleration is areaA


  ! ----------------------------------------------------------------------------
  ! --- 1. INTEGRATE FIRST (sometimes non-standard) REGION: area0 ---
!WRITE(*,*) "Initial area"
  CALL integrateFirstRegion(i, mfirst, leftOfMax,         & ! INPUTS
                            area0, zeroR)                   ! OUTPUTS

  IF (Cverbose) THEN
    CALL DBLEPR("Initial region area:", -1, area0, 1)
    CALL DBLEPR("      between 0 and:", -1, zeroR, 1)
    CALL INTPR( "      using right m:", -1, m, 1)
  END IF




  ! ----------------------------------------------------------------------------
  ! --- 2. INTEGRATE: the PRE-ACCELERATION regions: area1 ---
!WRITE(*,*) "Starting pre-acceleration"
  zeroL = zeroR  ! The last region's right-side zero is next region's left-side zero
  leftPreAccZero = zeroL
  CALL integratePreAccRegions(m, mfirst, leftOfMax, zeroL,  tmax, mmax,             & ! INPUTS
                              area1, zeroR, count_PreAcc_regions,  converged_Pre)     ! OUTPUTS
  count_Integration_Regions = count_Integration_Regions + count_PreAcc_regions

  IF (Cverbose) THEN
    CALL DBLEPR(" Pre-acc area:", -1, area1, 1)
    CALL DBLEPR("      between:", -1, leftPreAccZero, 1)
    CALL DBLEPR("          and:", -1, zeroR, 1)
    CALL INTPR( "using right m:", -1, m,     1)
  
    IF (converged_Pre) THEN
      CALL DBLEPR(" Accelerating not needed; convergence by t =", -1, zeroR, 1)
    ELSE
      CALL DBLEPR(" Accelerating starting after t =", -1, zeroR, 1)
    END IF
  END IF



  ! ----------------------------------------------------------------------------
  ! --- 3. INTEGRATE: the ACCELERATION regions: areaA ---
!WRITE(*,*) "Starting acceleration"
  IF ( converged_Pre) THEN
    areaA = 0.0_C_DOUBLE
    count_Acc_Regions = 0
  ELSE
    leftAccZero = zeroL
    zeroL = zeroR  ! The last region's right-side zero is next region's left-side zero
    CALL integrateAccelerationRegions(m, leftOfMax, zeroL,  tmax,                               & ! INPUTS
                                      areaA, zeroR, count_Acc_Regions, converged_Accelerating)    ! OUTPUTS
    IF (Cverbose) THEN
      IF ( .NOT.(converged_Accelerating) ) THEN
        CALL DBLEPR(" Accelerating dod not converged by t =", -1, zeroR, 1)
      END IF
      
      CALL DBLEPR("     Acc area:", -1, areaA, 1)
      CALL DBLEPR("      between:", -1, leftAccZero, 1)
      CALL DBLEPR("          and:", -1, zeroR, 1)
      CALL INTPR( "      up to m:", -1, m,     1)
    END IF  
  END IF
  

  ! --- WIND THINGS UP ---
  count_Integration_Regions = 1_C_INT  +              &   ! Initial zone has one integration region
                              count_PreAcc_Regions +  &   ! Pre-acc regions
                              count_Acc_Regions           ! Acc regions
  areaT = area0 + area1 + areaA

  IF (Cverbose) THEN
    CALL DBLEPR("* Initial area0: ", -1, area0, 1)
    CALL DBLEPR("* Pre-acc area1: ", -1, area1, 1)
    CALL DBLEPR("*     Acc area!: ", -1, areaA, 1)
    CALL DBLEPR("***       TOTAL: ", -1, areaT, 1)
    CALL INTPR( "   over regions: ", -1, count_Integration_Regions, 1)
  END IF

  ! We have the value of the integral in the PDF/CDF calculation.
  ! So now work out the actualy PDF/CDF
  
  IF (Cpdf) THEN
    funvalueI = areaT/pi 
  ELSE
    funvalueI = -areaT/pi + 0.5E0_C_DOUBLE
  END IF  
  IF (Cverbose) CALL DBLEPR("***    Fun. value:", -1, funvalueI, 1)
  

  CONTAINS
  
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE integrateFirstRegion(i, mfirst, leftOfMax,           & ! INPUTS
                                    area0, zeroR)                     ! OUTPUT
      ! Integrates the initial region of the Fourier integral

      IMPLICIT NONE
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: area0, zeroR
      INTEGER(C_INT), INTENT(IN)          :: mfirst, i
      LOGICAL(C_BOOL), INTENT(INOUT)      :: leftOfMax
      
      REAL(KIND=C_DOUBLE)                 :: t_Start_Point, zeroL, zeroBoundL, zeroBoundR
      REAL(KIND=C_DOUBLE)                 :: pi
      LOGICAL(C_BOOL)                     :: error


      pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)

      ! Grab the relevant scalar values for this iteration:
      current_y    = Cy(i)    ! Access y value for index i
      current_mu   = Cmu(i)   ! Access mu value for index i
      current_phi  = Cphi(i)  ! Access phi value for index i
  
      
      ! Initialisation
      t_Start_Point = 0.0_C_DOUBLE
      zeroL = 0.0_C_DOUBLE
      zeroBoundL = 0.0_C_DOUBLE
      zeroBoundR = 0.0_C_DOUBLE

!WRITE(*,*) "leftOfMax:", leftOfmax
      ! Find starting point for the first zero
      IF (leftOfMax) THEN
        t_Start_Point = pi / current_y  
        zeroBoundL = 0
        zeroBoundR = tmax   ! WAS: t_Start_Point * 2.0_C_DOUBLE
      ELSE
        ! Searching to the right of tmax
        t_Start_Point = tmax + pi / current_y  
        zeroBoundL = tmax
        zeroBoundR = t_Start_Point * 2.0_C_DOUBLE
      END IF
!WRITE(*,*) "BOUNDS", zeroBoundL, t_Start_Point, zeroBoundR
      IF ( (t_Start_Point .GT. zeroBoundR) .OR. (t_Start_Point .LT. zeroBoundL) ) Then
        t_Start_Point = (zeroBoundL + zeroBoundR) / 2.0_c_DOUBLE
      END IF
    
      ! Find the zero
      zeroL = 0.0_C_DOUBLE
!WRITE(*,*) "BOUNDS---about to call find exact zeros", zeroBoundL, t_Start_Point, zeroBoundR

      CALL findExactZeros(i, mfirst, zeroBoundL, zeroBoundR, t_Start_Point, zeroR, leftOfMax)
!WRITE(*,*) "FOUND---zeroR", zeroR
      CALL evaluateImk(i, zeroR, TMP, error)
      IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, zeroR, 1)

      ! Find the area
      CALL GaussQuadrature(i, zeroL, zeroR, area0)
!WRITE(*,*) "From", zeroL, "to", zeroR, "(m=", m, "), area=", area0
      IF (Cverbose) THEN
        CALL DBLEPR("  *** INITIAL area:", -1, area0, 1 )
        CALL DBLEPR("         between t:", -1, 0.0_C_DOUBLE, 1 )
        CALL DBLEPR("             and t:", -1, zeroR, 1 )
        CALL INTPR( "   using (right) m:", -1, mfirst, 1)
      END IF

      END SUBROUTINE integrateFirstRegion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE integratePreAccRegions(m, mfirst, leftOfMax, zeroL,  tmax, mmax,                     & ! INPUTS
                                     area1, zeroR, count_PreAcc_regions, converged_Pre)               ! OUTPUTS
      ! Integration of the region *after* the initial, but *before* acceleration is invoked.
      ! Potentially, everything converges in this step (without ever needing acceleration).
      
      IMPLICIT NONE
      
      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: area1, zeroR
      REAL(KIND=C_DOUBLE), INTENT(IN)     :: tmax
      REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: zeroL
      INTEGER(C_INT), INTENT(IN)          :: mfirst, mmax
      INTEGER(C_INT), INTENT(INOUT)       :: m
      LOGICAL(C_BOOL), INTENT(INOUT)      :: leftOfMax
      INTEGER(C_INT), INTENT(OUT)         :: count_PreAcc_Regions
      LOGICAL(C_BOOL), INTENT(OUT)        :: converged_Pre

      INTEGER(C_INT)                      :: mOld
      REAL(KIND=C_DOUBLE)                 :: zeroBoundL, zeroBoundR
      REAL(KIND=C_DOUBLE)                 :: area1Old, tolerance, sumAOld
      LOGICAL(C_BOOL)                     :: stop_PreAccelerate
  
  
!      pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
      converged_Pre = .FALSE.           ! .TRUE. if it seems the integration has converged
      tolerance = 1.0E-12_C_DOUBLE      ! tolerance for concluding the area is not vhanging
      sumAOld = 0.0_C_DOUBLE            ! Previous estimation of the area
      count_PreAcc_Regions = 0_C_INT    ! Count how many pre-acc regions are evaluated

      IF (mfirst .EQ. -1) THEN
        area1 = 0.0_C_DOUBLE
        mOld = m
        zeroR = zeroL 

        IF (Cverbose) CALL DBLEPR("  - No PRE-ACC area for y:", -1, current_y, 1 )
      ELSE
        CALL advanceM(i, m, mmax, mOld, leftOfMax, flip_To_Other_Side)
        area1 = 0.0_C_DOUBLE
        area1Old = 10.0_C_DOUBLE
        mOld = m
    
        stop_PreAccelerate = .FALSE.
        DO WHILE ( .NOT.(stop_PreAccelerate ) )
          IF (leftOfMax) THEN
            zeroBoundL = zeroR
            zeroBoundR = tmax
            IF (flip_To_Other_Side) THEN 
              leftOfMax = .FALSE.
              flip_To_Other_Side = .FALSE.
            END IF
          ELSE
            zeroBoundL = tmax
            zeroBoundR = zeroBoundL * 20.0_C_DOUBLE
          END IF
          zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE
          zeroL = zeroR
          CALL improveKZeroBounds(i, m, leftOfMax, zeroStartPoint, &
                                  zeroBoundL, zeroBoundR)
          zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE
          CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, &
                              zeroStartPoint, zero, leftOfMax)
          zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE
          zeroR = zero
    
          IF (count_PreAcc_Regions .GT. 1) THEN
            area1Old = area1
            sumAOld = sumA
          END IF

          CALL GaussQuadrature(i, zeroL, zeroR, sumA)
!WRITE(*,*) "From", zeroL, "to", zeroR, "(m=", m, "), area=", sumA

          area1 = area1 + sumA
          count_PreAcc_Regions = count_PreAcc_Regions + 1

          IF (Cverbose) THEN 
            CALL DBLEPR("   Pre-acc:  Left t:", -1, zeroL, 1)
            CALL DBLEPR("   Pre-acc: Right t:", -1, zeroR, 1)
            CALL DBLEPR("               sumA:", -1, sumA, 1)
            CALL DBLEPR("              area1:", -1, area1, 1)
            CALL INTPR( "            with m::", -1, m, 1)
          END IF          
!CALL INTPR( "       >>>>>     with m::", -1, m, 1)
          
          ! Determine if we can stop pre-accelerating yet
          CALL stopPreAcc(tmax, zeroR, stop_PreAccelerate)

          mOld = m
          CALL advanceM(i, m, mmax, mOld, leftOfMax, flip_To_Other_Side)
        END DO 
    
        IF (Cverbose) THEN
          CALL DBLEPR("  *** PRE-ACC area:", -1, sumA, 1 )
          CALL DBLEPR("         between t:", -1, zeroL, 1 )
          CALL DBLEPR("             and t:", -1, zeroR, 1 )
          CALL INTPR( "           using m:", -1, m, 1 )
        END IF
      END IF

    END SUBROUTINE integratePreAccRegions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE integrateAccelerationRegions(m, leftOfMax, zeroL,  tmax, & ! INPUTS
                                            West, zeroR, its_Acceleration, converged_Accelerating)
      ! Integration of the region *after* the initial and pre-acceleration regions, the tail area.
      ! Potentially, everything converges in pre-acceleration (without needing acceleration).


      IMPLICIT NONE

      REAL(KIND=C_DOUBLE), INTENT(OUT)    :: zeroR
      REAL(KIND=C_DOUBLE)                 :: tmax, zeroL
      INTEGER(C_INT), INTENT(INOUT)       :: m
      LOGICAL(C_BOOL), INTENT(INOUT)      :: leftOfMax
      INTEGER(C_INT), INTENT(OUT)         :: its_Acceleration

      INTEGER(C_INT)                      :: mOld
      REAL(KIND=C_DOUBLE)                 :: zeroBoundL, zeroBoundR
      REAL(KIND=C_DOUBLE)                 :: West, Wold, Wold2
      LOGICAL(C_BOOL)                     :: keep_Accelerating, converged_Accelerating


      ! Initialisation of acceleration 
      West = 3.0_C_DOUBLE                 ! The current estimate of the tail area
      Wold = 2.0_C_DOUBLE                 ! Wold and Wold2 are the previous two estimates of the tail area,
      Wold2 = 1.0_C_DOUBLE                !   to allow testing for convergence
      its_Acceleration = 0                ! Number of integration regions used
      keep_Accelerating = .TRUE.          ! .TRUE. if we need more integration regions
      converged_Accelerating = .FALSE.    ! .TRUE. if it seems the integration has converged
      its_Acceleration = 0_C_INT          ! Number of acceleration integration regions
      
      accMax = 100_C_INT                  ! Maximum number of regions in acceleration; arbitrary
      min_Acc_Regions = 3_C_INT           ! Minimum number of acceleration regions to use; need at least three

      ! This will be the very first, left-most value of t used in acceleration
      xvec(1) = zeroL



!      CALL advanceM(i, m, mmax, mOld, leftOfMax, flip_To_Other_Side)

      DO WHILE (keep_Accelerating )
        its_Acceleration = its_Acceleration + 1_C_INT
        mOld = m
        CALL advanceM(i, m, mmax, mOld, leftOfMax, flip_To_Other_Side)
        zeroBoundL = zeroL
    
        IF (leftOfMax) THEN
          zeroBoundR     = tmax
          zeroStartPoint = (zeroBoundL + zeroBoundR) / 2.0_C_DOUBLE
        ELSE
          ! Searching to the right of the maximum of Im k(t)
          IF (flip_To_Other_Side) THEN 
            ! Then this is the first time on the right side of the maximum
            leftOfMax = .FALSE.
            flip_To_Other_Side = .FALSE.
            zeroBoundL = tmax
            zeroBoundR = tmax * 2.0_C_DOUBLE
          ELSE
            zeroStartPoint = zeroL + (pi / current_y)
            zeroBoundR = zeroStartPoint * 5.0_C_DOUBLE
          END IF
        END IF

        ! Find the exact zero
        CALL findExactZeros(i, m, zeroBoundL, zeroBoundR, zeroStartPoint, zeroR, leftOfMax)
        IF (its_Acceleration .GE. accMax) THEN
          IF (Cverbose) CALL INTPR("Max acceleration regions reached. Stopping acceleration at m:", -1, m, 1)
          convergence_Acc = .TRUE.
          EXIT
        END IF
    
        xvec(its_Acceleration + 1) = zeroR
        
        CALL GaussQuadrature(i, zeroL, zeroR, psi)           ! psi: area of the latest region
!WRITE(*,*) "From", zeroL, "to", zeroR, "(m=", m, "), area=", psi
        wvec(its_Acceleration) = psi 
          ! wvec contains the sequence of integration areas, starting with the first and up to the limit
    
        ! Update past estimates
        Wold2 = Wold
        Wold  = West    
        ! Estimate the tail area using acceleration
        CALL accelerate(xvec, wvec, its_Acceleration, Mmatrix, Nmatrix, West)

        ! Check for convergence
        relerr = ( DABS(West - Wold) + DABS(West - Wold2)) / (DABS(West) + epsilon)
        IF (Cverbose) THEN
          CALL DBLEPR("  *** This acceleration area:", -1, psi, 1)
          CALL DBLEPR("                   between t:", -1, zeroL, 1)
          CALL DBLEPR("                       and t:", -1, zeroR, 1)
          CALL INTPR( "                     using m:", -1, m, 1)
          CALL DBLEPR("                with rel err:", -1, relerr, 1)
        END IF
    
        ! Declare convergence_Acc of we have sufficient regions, and relerr estimate is small
        IF ( (its_Acceleration .GE. min_Acc_Regions) .AND. &
             (relerr .LT. aimrerr) ) THEN
          keep_Accelerating = .FALSE.
          converged_Accelerating = .TRUE.
          convergence_Acc = .TRUE.
        END IF
        
        IF (its_Acceleration .GT. accMax) THEN
             keep_Accelerating = .FALSE.
             converged_Accelerating = .FALSE.
        END IF
!IF (m .LT. -13) STOP

        zeroL = zeroR ! Next left zero, is previous region's right zero

      END DO  
    
    END SUBROUTINE integrateAccelerationRegions
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    SUBROUTINE stopPreAcc(tmax, zero, stop_PreAccelerate)
      ! Determine if it is OK to stop pre-accelerating, and start using acceleration
      
      IMPLICIT NONE
    
      REAL(KIND=C_DOUBLE), INTENT(IN) :: tmax, zero
      LOGICAL(C_BOOL), INTENT(OUT)    :: stop_PreAccelerate
      
      INTEGER (C_INT)       :: nmax
      REAL(KIND=C_DOUBLE)   :: MM, Rek, Rekd, tstop
      
      ! Stop condition for pre-acceleration.
      ! Ensure that we have passed the peak of Im k(t), so that acceleration can be used
      stop_PreAccelerate = .FALSE.
      IF ( Cp .GT. 2) THEN
        IF (current_y .GT. current_mu) THEN
          stop_PreAccelerate = .TRUE.
        ELSE
          IF (zero .GT. tmax) THEN
            stop_PreAccelerate = .TRUE.
          END IF
        END IF
      ELSE
        IF (current_y .GT. current_mu) THEN
          stop_PreAccelerate = .TRUE.
        ELSE
          MM = 1.0_C_DOUBLE / (2.0_C_DOUBLE  * (Cp - 1.0_C_DOUBLE))
          IF ( MM .LE. 1.0_C_DOUBLE ) THEN
            stop_PreAccelerate = .TRUE.
          ELSE
          ! Check when t is larger than the last turning point of Re k(t) 
          ! Or when exp{Re(k)/t} is so small that it makes no difference... but 
          ! care is needed: Re k(t) is not necessarily convex here
             IF ( ABS( MM - FLOOR(MM) ) < 1.0E-09_C_DOUBLE ) then
                nmax = FLOOR(MM) - 1
             ELSE
                nmax = FLOOR(MM)
             END IF
             tstop = current_mu**(1.0_C_DOUBLE - Cp) / ((1.0_C_DOUBLE - Cp) * current_phi) *   & 
                     DTAN( DBLE(nmax) * pi * (1.0_C_DOUBLE - Cp) )
             IF (zeroL .GT. tstop) stop_PreAccelerate = .TRUE.
          END IF
        END IF
      END IF
      
      ! Sometimes this takes forever to flag stop_preacc as TRUE,
      ! so also check if exp{Re k(t)/t} is small
      CALL evaluateRek( i, zero, Rek)
      CALL evaluateRekd(i, zero, Rekd)
      IF ( ( (DEXP(Rek)/zeroL) .LT. 1.0E-07_C_DOUBLE) .AND.          & 
           (Rekd .LT. 0.0_C_DOUBLE) ) then
        stop_PreAccelerate = .TRUE.
        IF (Cverbose) THEN 
          CALL DBLEPR("Pre-accelerating stopping. Rek(t) small:", -1, (DEXP(Rek)/zeroL), 1)
        END IF
      END IF

      IF ( ( (DEXP(Rek)/zeroL) .LT. 1.0E-15_C_DOUBLE)  ) then
        stop_PreAccelerate = .TRUE.
        IF (Cverbose) THEN 
          CALL DBLEPR("Pre-accelerating stopping. Rek(t) small:", -1, (DEXP(Rek)/zeroL), 1)
        END IF
      END IF


    END SUBROUTINE stopPreAcc
    
END SUBROUTINE TweedieIntegration

