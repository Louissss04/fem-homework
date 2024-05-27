      DOUBLE PRECISION FUNCTION DPMEPS()
C     **********
C
C     SUBROUTINE DPMEPS
C
C     THIS SUBROUTINE COMPUTES THE MACHINE PRECISION PARAMETER
C     DPMEPS AS THE SMALLEST FLOATING POINT NUMBER SUCH THAT
C     1 + DPMEPS DIFFERS FROM 1.
C
C     THIS SUBROUTINE IS BASED ON THE SUBROUTINE MACHAR DESCRIBED IN
C
C     W. J. CODY, MACHAR: A SUBROUTINE TO DYNAMICALLY DETERMINE
C     MACHINE PARAMETERS, ACM TRANS. MATH. SOFTWARE 14 (1988), 303-311.
C
C     THE FUNCTION STATEMENT IS
C
C       DOUBLE PRECISION DPMEPS()
C
C     MINPACK-2 PROJECT. FEBRUARY 1991.
C     ARGONNE NATIONAL LABORATORY AND UNIVERSITY OF MINNESOTA.
C     BRETT M. AVERICK.
C
C     **********
      INTEGER I, IBETA, IRND, IT, ITEMP, NEGEP
      DOUBLE PRECISION A, B, BETA, BETAIN, BETAH, TEMP, TEMPA, TEMP1
      DOUBLE PRECISION ZERO, ONE, TWO
      DATA ZERO, ONE, TWO/0.0D0, 1.0D0, 2.0D0/

C     DETERMINE IBETA, BETA ALA MALCOLM.

      A = ONE
      B = ONE
   10 CONTINUE
      A = A + A
      TEMP = A + ONE
      TEMP1 = TEMP - A
      IF (TEMP1-ONE .EQ. ZERO) GO TO 10
   20 CONTINUE
      B = B + B
      TEMP = A + B
      ITEMP = INT(TEMP-A)
      IF (ITEMP .EQ. 0) GO TO 20
      IBETA = ITEMP
      BETA = DBLE(IBETA)

C     DETERMINE IT, IRND.

      IT = 0
      B = ONE
   30 CONTINUE
      IT = IT + 1
      B = B*BETA
      TEMP = B + ONE
      TEMP1 = TEMP - B
      IF (TEMP1-ONE .EQ. ZERO) GO TO 30
      IRND = 0
      BETAH = BETA/TWO
      TEMP = A + BETAH
      IF (TEMP-A .NE. ZERO) IRND = 1
      TEMPA = A + BETA
      TEMP = TEMPA + BETAH
      IF ((IRND .EQ. 0) .AND. (TEMP-TEMPA .NE. ZERO)) IRND = 2

C     DETERMINE DPMEPS.

      NEGEP = IT + 3
      BETAIN = ONE/BETA
      A = ONE
      DO 40 I = 1, NEGEP
         A = A*BETAIN
   40 CONTINUE
   50 CONTINUE
      TEMP = ONE + A
      IF (TEMP-ONE .NE. ZERO) GO TO 60
      A = A*BETA
      GO TO 50
   60 CONTINUE
      DPMEPS = A
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 70
      A = (A*(ONE+A))/TWO
      TEMP = ONE + A
      IF (TEMP-ONE .NE. ZERO) DPMEPS = A

   70 CONTINUE

      END
