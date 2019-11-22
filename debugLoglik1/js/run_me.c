int run_me(int lengthBuffer, double covarLength, double * yout, double * times,  double * covar1_p, double * covar2_p, double y1, double y2, double y3, double y4, double y5, double y6, double y7, double y8, double y9, double y10, double y11, double y12,
           double y13, double y14 , double data0, double data1, double data2, double data3, double data4 , double data5 , double data6 , double data7 , double data8 ,
            double data9 , double data10 , double data11 , double data12 , double data13 , double data14 , double data15, double data16, double data17 , double data18 , double data19 , double data20 , double data21 ,
             double data22 , double data23 , double data24 , double data25 , double data26)
{
  double          rwork1, rwork5, rwork6, rwork7;
  double          atol[15], rtol[15], t, y[15];
  int             iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9;
  int             neq = 14;
  int             itol, itask, istate, iopt, jt, iout, ijs;
  double          data[28];
  double          tout;
  iwork1 = iwork2 = iwork5 = iwork6 = iwork7 = iwork8 = iwork9 = 0;
  rwork1 = rwork5 = rwork6 = rwork7 = 0.0;
  y[1] = y1;
  y[2] = y2;
  y[3] = y3;
  y[4] = y4;
  y[5] = y5;
  y[6] = y6;
  y[7] = y7;
  y[8] = y8;
  y[9] = y9;
  y[10] = y10;
  y[11] = y11;
  y[12] = y12;
  y[13] = y13;
  y[14] = y14;
  data[0] = data0;                                
  data[1] = data1;
  data[2] = data2;
  data[3] = data3;
  data[4] = data4;
  data[5] = data5;
  data[6] = data6;
  data[7] = data7;
  data[8] = data8;
  data[9] = data9;
  data[10] = data10;
  data[11] = data11;
  data[12] = data12;
  data[13] = data13;
  data[14] = data14;
  data[15] = data15;
  data[16] = data16;
  data[17] = data17;
  data[18] = data18;
  data[19] = data19;
  data[20] = data20;
  data[21] = data21;
  data[22] = data22;
  data[23] = data23;
  data[24] = data24;
  data[25] = data25;
  data[26] = data26;
  data[27] = covarLength;
 
 
  itol = 2;
  rtol[0] = 0.0;
  atol[0] = 0.0;
  rtol[1] = rtol[3] = 1.0E-6;
  rtol[2] = rtol[4] = rtol[5] = rtol[6] = rtol[7] = rtol[8] = rtol[9] = rtol[10] = rtol[11] =rtol[12] = rtol[13] = rtol[14] = 1.0E-6;
  atol[1] = 1.0E-6;
  atol[2] = 1.0E-6;
  atol[3] = 1.0E-6;
  atol[4] = 1.0E-6;
  atol[5] = 1.0E-6;
  atol[6] = 1.0E-6;
  atol[7] = 1.0E-6;
  atol[8] = 1.0E-6;
  atol[9] = 1.0E-6;
  atol[10] = 1.0E-6;
  atol[11] = 1.0E-6;
  atol[12] = 1.0E-6;
  atol[13] = 1.0E-6;
  atol[14] = 1.0E-6;
  
  itask = 1;
  istate = 1;
  iopt = 0;
  jt = 2;
  for (iout = 0; iout < lengthBuffer; iout++) {
    t = times[iout];
    tout = times[iout + 1];
    lsoda(skel, neq, y, &t, tout, itol, rtol, atol, itask, &istate, iopt, jt,
          iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9,
          rwork1, rwork5, rwork6, rwork7, data, covar1_p, covar2_p);
    yout[iout] = y[14];
    // printf("%14.6e \n", yout[iout]);
    if (istate <= 0) {
      printf("error istate = %d\n", istate);
      // exit(0);
      n_lsoda_terminate();
      return -1;
    }
  } 
  n_lsoda_terminate();

  return 0;
}
