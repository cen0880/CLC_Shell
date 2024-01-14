
/**
 * @file   calculate_Fel.h
 * @author Yucen Han
 * @date   Jan 14 2024
 * @brief  modified from 1Dimension's nematic shell code
 * 
 * 
 */
#ifndef _CALCULATE_FEL_H
#define _CALCULATE_FEL_H

using namespace std;

#include <iostream>
#include "global.h"
#include "Basis.h"

double calc_F_el(double *qnlm, double *qijk)  
{

  double a,b,c;
  double trQ2,trQ3,trQ4;
  double L1,L2,L3;

  //Dimensionless parameters in LdG functional
  a = 0.5 * landau_t;
  b = - sqrt(6);
  c = 0.5;
  L1 = 1.0/Rad/Rad;
  L2 = etaL*L1;
  L3 = 1.0; 
  int i, j, k;
  double intG;
  double temp_x,temp_y,temp_z;

  double b1, b2, b3;

  double c11, c12, c13, c21, c22, c23, c31, c32, c33;

  double q1, q2, q3, q4, q5;

  double s0 = sqrt(1.5) * (3.0 + sqrt(9.0 - 8 * landau_t))/4;
  double f_bulkmin = s0 * s0 / 54*(9 * landau_t - 3*sqrt(6)*s0);
  
  double gamma_dr = 1;

  if(landau_t > 9.0/8)
  {
    s0 = 0;
  }

  for(i = 0; i < 5; i++)
  {
    calc_fijk(qnlm + i * Basis, qijk + i * Point);
    calc_dr_fijk(qnlm + i * Basis, dr_qijk + i * innerPoint);
    calc_dt_fijk(qnlm + i * Basis, dt_qijk + i * innerPoint);
    calc_dp_fijk(qnlm + i * Basis, dp_qijk + i * innerPoint);
  }

  for(int ix = 0; ix < innerPoint; ix++)
  {
    q1 = qijk[ix];
    q2 = qijk[Point + ix];
    q3 = qijk[2 * Point + ix];
    q4 = qijk[3 * Point + ix];
    q5 = qijk[4 * Point + ix];

    Q1x[ix] = dr_qijk[ix] * drdx[ix] + dt_qijk[ix] * dtdx[ix] + dp_qijk[ix] * dpdx[ix];
    Q1y[ix] = dr_qijk[ix] * drdy[ix] + dt_qijk[ix] * dtdy[ix] + dp_qijk[ix] * dpdy[ix];
    Q1z[ix] = dr_qijk[ix] * drdz[ix] + dt_qijk[ix] * dtdz[ix] + dp_qijk[ix] * dpdz[ix];

    Q2x[ix] = dr_qijk[ix + innerPoint] * drdx[ix] + dt_qijk[ix + innerPoint] * dtdx[ix] + dp_qijk[ix + innerPoint] * dpdx[ix];
    Q2y[ix] = dr_qijk[ix + innerPoint] * drdy[ix] + dt_qijk[ix + innerPoint] * dtdy[ix] + dp_qijk[ix + innerPoint] * dpdy[ix];
    Q2z[ix] = dr_qijk[ix + innerPoint] * drdz[ix] + dt_qijk[ix + innerPoint] * dtdz[ix] + dp_qijk[ix + innerPoint] * dpdz[ix];

    Q3x[ix] = dr_qijk[ix + 2 * innerPoint] * drdx[ix] + dt_qijk[ix + 2 * innerPoint]*dtdx[ix] + dp_qijk[ix + 2 * innerPoint]*dpdx[ix];
    Q3y[ix] = dr_qijk[ix + 2 * innerPoint] * drdy[ix] + dt_qijk[ix + 2 * innerPoint]*dtdy[ix] + dp_qijk[ix + 2 * innerPoint]*dpdy[ix];
    Q3z[ix] = dr_qijk[ix + 2 * innerPoint] * drdz[ix] + dt_qijk[ix + 2 * innerPoint]*dtdz[ix] + dp_qijk[ix + 2 * innerPoint]*dpdz[ix];

    Q4x[ix] = dr_qijk[ix + 3 * innerPoint] * drdx[ix] + dt_qijk[ix + 3 * innerPoint]*dtdx[ix] + dp_qijk[ix + 3 * innerPoint]*dpdx[ix];
    Q4y[ix] = dr_qijk[ix + 3 * innerPoint] * drdy[ix] + dt_qijk[ix + 3 * innerPoint]*dtdy[ix] + dp_qijk[ix + 3 * innerPoint]*dpdy[ix];
    Q4z[ix] = dr_qijk[ix + 3 * innerPoint] * drdz[ix] + dt_qijk[ix + 3 * innerPoint]*dtdz[ix] + dp_qijk[ix + 3 * innerPoint]*dpdz[ix];

    Q5x[ix] = dr_qijk[ix + 4 * innerPoint] * drdx[ix] + dt_qijk[ix + 4 * innerPoint]*dtdx[ix] + dp_qijk[ix + 4 * innerPoint]*dpdx[ix];	
    Q5y[ix] = dr_qijk[ix + 4 * innerPoint] * drdy[ix] + dt_qijk[ix + 4 * innerPoint]*dtdy[ix] + dp_qijk[ix + 4 * innerPoint]*dpdy[ix];
    Q5z[ix] = dr_qijk[ix + 4 * innerPoint] * drdz[ix] + dt_qijk[ix + 4 * innerPoint]*dtdz[ix] + dp_qijk[ix + 4 * innerPoint]*dpdz[ix];

    trQ2 = 2*(q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4);
    trQ3 = 3*(2*q2*q3*q5 - (q1*q1-q2*q2+q3*q3)*q4 + q1*(q2*q2-q4*q4-q5*q5));
    trQ4 = trQ2*trQ2;

    fbulk[ix] = a*trQ2 + b*trQ3 + c*trQ4 - f_bulkmin;

    grad_fbulk[ix] = a * (4*q1+2*q4) + b*3*(q2*q2-q4*q4-q5*q5-2*q1*q4) + c*2*trQ2*(4*q1+2*q4);
    grad_fbulk[innerPoint + ix] = a * 4*q2 + b*6*(q3*q5 + q2*q4 + q1*q2) + c*2*trQ2*4*q2;
    grad_fbulk[2 * innerPoint + ix] = a * 4*q3 + b*6*(q2*q5 - q3*q4) + c*2*trQ2*4*q3;
    grad_fbulk[3 * innerPoint + ix] = a * (4*q4+2*q1) + b*3*(q2*q2-q1*q1-q3*q3-2*q1*q4) + c*2*trQ2*(4*q4+2*q1);
    grad_fbulk[4 * innerPoint + ix] = a * 4*q5 + b*6*(q2*q3 - q1*q5) + c*2*trQ2*4*q5;
     
    b1 = Q1x[ix] + Q2y[ix] + Q3z[ix];
    b2 = Q2x[ix] + Q4y[ix] + Q5z[ix];
    b3 = Q3x[ix] + Q5y[ix] + (-Q1z[ix]-Q4z[ix]);

    c11 = Q3y[ix] - Q2z[ix] + 2*q0*q1;
    c12 = Q5y[ix] - Q4z[ix] + 2*q0*q2;
    c13 = (-Q1y[ix]-Q4y[ix]) - Q5z[ix] + 2*q0*q3;
    c21 = -Q3x[ix] + Q1z[ix] + 2*q0*q2;
    c22 = -Q5x[ix] + Q2z[ix] + 2*q0*q4;
    c23 = -(-Q1x[ix]-Q4x[ix]) + Q3z[ix] + 2*q0*q5;
    c31 = Q2x[ix] - Q1y[ix] + 2*q0*q3;
    c32 = Q4x[ix] - Q2y[ix] + 2*q0*q5;
    c33 = Q5x[ix] - Q3y[ix] + 2*q0*(-q1-q4);

    QL1[ix] = 0.5 * L2 * (b1*b1 + b2*b2 + b3*b3) + 0.5 * L1 * (c11*c11 + c12*c12 + c13*c13 + c21*c21 + c22*c22 + c23*c23 + c31*c31 + c32*c32 + c33*c33) + L3 * fbulk[ix];

    temp_x = L2 * b1 + L1 * c23;
    temp_y = - L1 * c13 - L1 * c31;
    temp_z = - L2 * b3 + L1 * c21;

    grad_QL1[ix] = 2 * L1 * q0 * c11 - 2 * L1 * q0 * c33;
    grad_QL1_dr[ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];
    
    temp_x = L2 * b2 + L1 * c31;
    temp_y = L2 * b1 - L1 * c32;
    temp_z = - L1 * c11 + L1 * c22;

    grad_QL1[innerPoint + ix] = 2 * L1 * q0 * c12 + 2 * L1 * q0 * c21;
    grad_QL1_dr[innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];
	
    temp_x = L2 * b3 - L1 * c21;
    temp_y = L1 * c11 - L1 * c33;
    temp_z = L2 * b1 + L1 * c23;

    grad_QL1[2 * innerPoint + ix] = 2 * L1 * q0 * c13 + 2 * L1 * q0 * c31;
    grad_QL1_dr[2 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[2 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[2 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];

    temp_x = L1 * c23 + L1 * c32;
    temp_y = L2 * b2 - L1 * c13;
    temp_z = - L2 * b3 - L1 * c12;
    
    grad_QL1[3 * innerPoint + ix] = 2 * L1 * q0 * c22 - 2 * L1 * q0 * c33;
    grad_QL1_dr[3 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[3 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[3 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];

    temp_x = L1 * c33 - L1 * c22; 
    temp_y = L2 * b3 + L1 * c12;
    temp_z = L2 * b2 - L1 * c13;

    grad_QL1[4 * innerPoint + ix] = 2 * L1 * q0 * c23 + 2 * L1 * q0 * c32;
    grad_QL1_dr[4 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[4 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[4 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];

  }

  intG = IntBisphereJacobi(QL1);


  for (int i = 0; i < 5; i++) 
  {
    calc_gnlm(grad_fbulk + i * innerPoint, grad_Fb + i * Basis);
    calc_gnlm(grad_QL1 + i * innerPoint, grad_FeL + i * Basis);
    calc_dr_gnlm(grad_QL1_dr + i * innerPoint, grad_FeL1_r + i * Basis);
    calc_dt_gnlm(grad_QL1_dt + i * innerPoint, grad_FeL1_t + i * Basis);
    calc_dp_gnlm(grad_QL1_dp + i * innerPoint, grad_FeL1_p + i * Basis);
  }

  for (int i = 0; i < 5 * Basis; i++) 
  {
    grad_FeL1[i] = grad_FeL1_r[i] + grad_FeL1_t[i] + grad_FeL1_p[i];
  }
  // exit();

  return intG;

}


#endif
