// *DECK DUMACH
// C***BEGIN PROLOGUE  DUMACH
// C***PURPOSE  Compute the unit roundoff of the machine.
// C***CATEGORY  R1
// C***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
// C***KEYWORDS  MACHINE CONSTANTS
// C***AUTHOR  Hindmarsh, Alan C., (LLNL)
// C***DESCRIPTION
// C *Usage:
// C        DOUBLE PRECISION  A, DUMACH
// C        A = DUMACH()
// C
// C *Function Return Values:
// C     A : the unit roundoff of the machine.
// C
// C *Description:
// C     The unit roundoff is defined as the smallest positive machine
// C     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
// C     in a machine-independent manner.
// C
// C***REFERENCES  (NONE)
// C***ROUTINES CALLED  DUMSUM
// C***REVISION HISTORY  (YYYYMMDD)
// C   19930216  DATE WRITTEN
// C   19930818  Added SLATEC-format prologue.  (FNF)
// C   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
// C***END PROLOGUE  DUMACH
// C
let DUMACH = function () {
  let U, COMP
  U = 1;
  do{
    U = U * 0.5;
    DUMSUM(1, U, COMP);
    } while(COMP !== 1)

  return U*2;
}