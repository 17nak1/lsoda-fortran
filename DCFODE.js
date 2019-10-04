  
// C ***BEGIN PROLOGUE  DCFODE
// C ***SUBSIDIARY
// C ***PURPOSE  Set ODE integrator coefficients.
// C ***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
// C ***AUTHOR  Hindmarsh, Alan C., (LLNL)
// C ***DESCRIPTION

// C  DCFODE is called by the integrator routine to set coefficients
// C  needed there.  The coefficients for the current method, as
// C  given by the value of METH, are set for all orders and saved.
// C  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
// C  (A smaller value of the maximum order is also allowed.)
// C  DCFODE is called once at the beginning of the problem,
// C  and is not called again unless and until METH is changed.

// C  The ELCO array contains the basic method coefficients.
// C  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
// C  order nq are stored in ELCO(i,nq).  They are given by a genetrating
// C  polynomial, i.e.,
// C      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
// C  For the implicit Adams methods, l(x) is given by
// C      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
// C  For the BDF methods, l(x) is given by
// C      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
// C  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).

// C  The TESCO array contains test constants used for the
// C  local error test and the selection of step size and/or order.
// C  At order nq, TESCO(k,nq) is used for the selection of step
// C  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
// C  nq + 1 if k = 3.

// C ***SEE ALSO  DLSODE
// C ***ROUTINES CALLED  (NONE)
// C ***REVISION HISTORY  (YYMMDD)
// C   791129  DATE WRITTEN
// C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
// C   890503  Minor cosmetic changes.  (FNF)
// C   930809  Renamed to allow single/double precision versions. (ACH)
// C ***END PROLOGUE  DCFODE
// C **End
let DCFODE = function (METH, ELCO, TESCO) {
  let METH;
  let I, IB, NQ, NQM1, NQP1;
  let ELCO, TESCO;
  let AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ,RQFAC, RQ1FAC, TSIGN, XPIN;
  let ELCO = Array(13).fill(Array(12)), TESCO = Array(3).fill(Array(12))
  let PC = Array(12);
  let goto_variable = 0;
  while (true){
    switch (goto_variable){
      case 0:// FIRST EXECUTABLE STATEMENT  DCFODE
        if (METH === 1){
          goto_variable = 100;
        } else if (METH === 2) {
          goto_variable = 200;
        }
        break;
      case 100:
        ELCO[1][1] = 1;
        ELCO[2][1] = 1;
        TESCO[1][1] = 0;
        TESCO[2][1] = 2;
        TESCO[1][2] = 1;
        TESCO[3][1] = 0;
        PC[1] = 1;
        RQFAC = 1;
        for(let NQ = 2; NQ <= 12; NQ++) {
          RQ1FAC = RQFAC;
          RQFAC = RQFAC / NQ;
          NQM1 = NQ - 1;
          FNQM1 = NQM1;
          NQP1 = NQ + 1;
          // Form coefficients of p(x)*(x+nq-1). ----------------------------------
          PC(NQ) = 0;
          for(let IB = 1; IB <= NQM1; IB++) {
            I = NQP1 - IB;
            PC[I] = PC[I-1] + FNQM1 * PC[I];
          }
          PC[1] = FNQM1 * PC[1];
          // Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
          PINT = PC[1];
          XPIN = PC[1] / 2;
          TSIGN = 1;
          for(let I = 2; I <= NQ; I++) {
            TSIGN = -TSIGN;
            PINT = PINT + TSIGN * PC[I]/I;
            XPIN = XPIN + TSIGN * PC[I]/(I+1);
          }
          // Store coefficients in ELCO and TESCO. --------------------------------
          ELCO[1][NQ] = PINT * RQ1FAC;
          ELCO[2][NQ] = 1;
          for (let I = 2; I <= NQ; I++) {
            ELCO[I+1][NQ] = RQ1FAC * PC[I] / I;
          }
          AGAMQ = RQFAC * XPIN;
          RAGQ = 1/AGAMQ;
          TESCO[2][NQ] = RAGQ;
          if (NQ < 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1;
          TESCO[3][NQM1] = RAGQ;
        }
        return;
      case 200:
        PC[1] = 1;
        RQ1FAC = 1;
        for (let NQ = 1; NQ <= 5; NQ++) {
          FNQ = NQ;
          NQP1 = NQ + 1;
          PC(NQP1) = 0;
          for (let IB = 1; IB <= NQ; IB++) {
          I = NQ + 2 - IB;
          PC(I) = PC(I-1) + FNQ*PC(I);
          }
          PC(1) = FNQ*PC(1);
          for(let I = 1; I <= NQP1; I++) {
            ELCO(I,NQ) = PC(I)/PC(2);
          }
          ELCO(2,NQ) = 1;
          TESCO(1,NQ) = RQ1FAC;
          TESCO(2,NQ) = NQP1/ELCO(1,NQ);
          TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ);
          RQ1FAC = RQ1FAC/FNQ;
        }
        return;
    }
  }
}
