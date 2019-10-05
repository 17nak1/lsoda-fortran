// *DECK DINTDY

// C***BEGIN PROLOGUE  DINTDY
// C***SUBSIDIARY
// C***PURPOSE  Interpolate solution derivatives.
// C***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
// C***AUTHOR  Hindmarsh, Alan C., (LLNL)
// C***DESCRIPTION
// C
// C  DINTDY computes interpolated values of the K-th derivative of the
// C  dependent variable vector y, and stores it in DKY.  This routine
// C  is called within the package with K = 0 and T = TOUT, but may
// C  also be called by the user for any K up to the current order.
// C  (See detailed instructions in the usage documentation.)
// C
// C  The computed values in DKY are gotten by interpolation using the
// C  Nordsieck history array YH.  This array corresponds uniquely to a
// C  vector-valued polynomial of degree NQCUR or less, and DKY is set
// C  to the K-th derivative of this polynomial at T.
// C  The formula for DKY is:
// C               q
// C   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
// C              j=K
// C  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
// C  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
// C  communicated by COMMON.  The above sum is done in reverse order.
// C  IFLAG is returned negative if either K or T is out of bounds.
// C
// C***SEE ALSO  DLSODE
// C***ROUTINES CALLED  XERRWD
// C***COMMON BLOCKS    DLS001
// C***REVISION HISTORY  (YYMMDD)
// C   791129  DATE WRITTEN
// C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
// C   890503  Minor cosmetic changes.  (FNF)
// C   930809  Renamed to allow single/double precision versions. (ACH)
// C   010418  Reduced size of Common block /DLS001/. (ACH)
// C   031105  Restored 'own' variables to Common block /DLS001/, to
// C           enable interrupt/restart feature. (ACH)
// C   050427  Corrected roundoff decrement in TP. (ACH)
// C***END PROLOGUE  DINTDY
// C**End
let DINTDY = function (T, K, YH, NYH, DKY, IFLAG) {
      let K, NYH, IFLAG;
      let T, YH, DKY;
      let YH = Array(NYH).fill([]), DKY = [];
      let IOWND, IOWNS,ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU;
      let ROWNS,CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND;
     //  COMMON /DLS001/ ROWNS(209),
     // 1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     // 2   IOWND(6), IOWNS(6),
     // 3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     // 4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     // 5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      let I, IC, J, JB, JB2, JJ, JJ1, JP1;
      let C, R, S, TP;
      let MSG // CHARACTER(LEN=80) 
      // FIRST EXECUTABLE STATEMENT  DINTDY
      let goto_variable = 0;
      while (true){
        switch (goto_variable){
          case 0:// FIRST EXECUTABLE STATEMENT  DCFODE
            IFLAG = 0
            if (K < 0 || K > NQ) {
              goto_variable = 80;
            }
            TP = TN - HU -  100*UROUND*SIGN(ABS(TN) + ABS(HU), HU)
            if ((T-TP)*(T-TN) > 0) {
              goto_variable = 90;
            }
            S = (T - TN)/H
            IC = 1
            if (K === 0) {
              goto_variable = 15;
            }
            JJ1 = L - K
            for (let JJ = JJ1; JJ <= NQ; JJ++) {
              IC = IC*JJ
            }
            break;
          case 15:



 15   C = IC
      DO 20 I = 1,N
        DKY(I) = C*YH(I,L)
 20   CONTINUE
      if (K === NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1,JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        if (K === 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1,J
          IC = IC*JJ
 30     CONTINUE
 35     C = IC
        DO 40 I = 1,N
          DKY(I) = C*YH(I,JP1) + S*DKY(I)
 40     CONTINUE
 50     CONTINUE
      if (K === 0) RETURN
 55   R = H**(-K)
      DO 60 I = 1,N
        DKY(I) = R*DKY(I)
 60   CONTINUE
      RETURN
C
 80   MSG = 'DINTDY-  K (=I1) illegal      '
      CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0, 0)
      IFLAG = -1
      RETURN
 90   MSG = 'DINTDY-  T (=R1) illegal      '
      CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      RETURN
C----------------------- END OF SUBROUTINE DINTDY ----------------------
      END