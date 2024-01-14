/**
 * @file   calculate_Fs.h
 * @author Yucen Han
 * @date   Jan 14 2024
 * 
 * @brief  calculate surface Energy, modified from 1Dimension's nematic shell code
 * 
 * 
 */
#include "global.h"
#include "Basis.h"
#include "initial.h"
#include <cmath>

void calc_fs(double *Qb, double *x_b, double *y_b, double *z_b, double *fs, double *grad_fs, int Ind) /*计算surf energy density*/
{
  // int n, j, k;
  double s0 = sqrt(1.5) * (3.0 + sqrt(9.0 - 8 * landau_t))/4;
  if(landau_t > 9.0/8)
  {
    s0 = 0;
  }

  double q1,q2,q3,q4,q5;

	if(strcmp(boundary, "arvtan") == 0)
	{
      for(int ix = 0; ix < I * K; ix++)
      {

	    q1 = (Qb[ix * 5 + 0] + s0/3) * x_b[ix] + Qb[ix * 5 + 1] * y_b[ix] + Qb[ix * 5 + 2] * z_b[ix];
	    q2 = Qb[ix * 5 + 1] * x_b[ix] + (Qb[ix * 5 + 3] + s0/3) * y_b[ix] + Qb[ix * 5 + 4] * z_b[ix];
	    q3 = Qb[ix * 5 + 2] * x_b[ix] + Qb[ix * 5 + 4] * y_b[ix] - (Qb[ix * 5 + 0] + Qb[ix * 5 + 3] - s0/3) * z_b[ix];
      
	    fs[ix] = q1 * q1 + q2 * q2 + q3 * q3 ; 
      
	    grad_fs[ix] = 2 * q1 * x_b[ix] - 2 * q3 * z_b[ix];
	    grad_fs[I * K + ix] = 2 * q1 * y_b[ix] + 2 * q2 * x_b[ix];
	    grad_fs[2 * I * K + ix] = 2 * q1 * z_b[ix]  + 2 * q3 * x_b[ix];
	    grad_fs[3 * I * K + ix] = 2 * q2 * y_b[ix] - 2 * q3 * z_b[ix];
		  grad_fs[4 * I * K + ix] = 2 * q2 * z_b[ix] + 2 * q3 * y_b[ix];
	    }
	}
  else if(strcmp(boundary, "vertical") == 0)
{
    
      for(int ix = 0; ix < I * K; ix++)
      {

        q1 = Qb[ix * 5 + 0] - s0 * (x_b[ix] * x_b[ix] - 1.0/3);
        q2 = Qb[ix * 5 + 1] - s0 * x_b[ix] * y_b[ix];
        q3 = Qb[ix * 5 + 2] - s0 * x_b[ix] * z_b[ix];
        q4 = Qb[ix * 5 + 3] - s0 * (y_b[ix] * y_b[ix] - 1.0/3);
        q5 = Qb[ix * 5 + 4] - s0 * y_b[ix] * z_b[ix];

        fs[ix] = 2 * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4);
    
        grad_fs[ix] = 2 * (2 * q1 + q4);
        grad_fs[I * K + ix] = 2*(2 * q2);
        grad_fs[2 * I * K + ix] = 2*(2 * q3);
        grad_fs[3 * I * K + ix] = 2*(2 * q4 + q1);
        grad_fs[4 * I * K + ix] = 2*(2 * q5);  
      } 
}
else if(strcmp(boundary, "mix") == 0) 
{
    if(Ind == 2)
    {
      for(int ix = 0; ix < I * K; ix++)
      {
    
  q1 = Qb[ix * 5 + 0] - s0 * (x_b[ix] * x_b[ix] - 1.0/3);
  q2 = Qb[ix * 5 + 1] - s0 * x_b[ix] * y_b[ix];
  q3 = Qb[ix * 5 + 2] - s0 * x_b[ix] * z_b[ix];
  q4 = Qb[ix * 5 + 3] - s0 * (y_b[ix] * y_b[ix] - 1.0/3);
  q5 = Qb[ix * 5 + 4] - s0 * y_b[ix] * z_b[ix];

  fs[ix] = 2 * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4);
    
  grad_fs[ix] = 2 * (2 * q1 + q4);
  grad_fs[I * K + ix] = 2*(2 * q2);
  grad_fs[2 * I * K + ix] = 2*(2 * q3);
  grad_fs[3 * I * K + ix] = 2*(2 * q4 + q1);
  grad_fs[4 * I * K + ix] = 2*(2 * q5);  

      }
    }
    else
    {
      for(int ix = 0; ix < I * K; ix++)
      {
  q1 = (Qb[ix * 5 + 0] + s0/3) * x_b[ix] + Qb[ix * 5 + 1] * y_b[ix] + Qb[ix * 5 + 2] * z_b[ix];
  q2 = Qb[ix * 5 + 1] * x_b[ix] + (Qb[ix * 5 + 3] + s0/3) * y_b[ix] + Qb[ix * 5 + 4] * z_b[ix];
  q3 = Qb[ix * 5 + 2] * x_b[ix] + Qb[ix * 5 + 4] * y_b[ix] - (Qb[ix * 5 + 0] + Qb[ix * 5 + 3] - s0/3) * z_b[ix];

  fs[ix] = (q1 * q1 + q2 * q2 + q3 * q3); // / (R1 * R1);

  grad_fs[ix] = (2 * q1 * x_b[ix] - 2 * q3 * z_b[ix]); // / (R1 * R1);
  grad_fs[I * K + ix] = (2 * q1 * y_b[ix] + 2 * q2 * x_b[ix]); // / (R1 * R1);
  grad_fs[2 * I * K + ix] = (2 * q1 * z_b[ix]  + 2 * q3 * x_b[ix]); // / (R1 * R1);
  grad_fs[3 * I * K + ix] = (2 * q2 * y_b[ix] - 2 * q3 * z_b[ix]); // / (R1 * R1);
  grad_fs[4 * I * K + ix] = (2 * q2 * z_b[ix] + 2 * q3 * y_b[ix]); // / (R1 * R1);
      }
    }
}
else if(strcmp(boundary, "mix2") == 0) 
{
    if(Ind == 1)
    {
      for(int ix = 0; ix < I * K; ix++)
      {
    
  q1 = Qb[ix * 5 + 0] - s0 * (x_b[ix] * x_b[ix] - 1.0/3);
  q2 = Qb[ix * 5 + 1] - s0 * x_b[ix] * y_b[ix];
  q3 = Qb[ix * 5 + 2] - s0 * x_b[ix] * z_b[ix];
  q4 = Qb[ix * 5 + 3] - s0 * (y_b[ix] * y_b[ix] - 1.0/3);
  q5 = Qb[ix * 5 + 4] - s0 * y_b[ix] * z_b[ix];

  fs[ix] = 2 * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4);
    
  grad_fs[ix] = 2 * (2 * q1 + q4);
  grad_fs[I * K + ix] = 2*(2 * q2);
  grad_fs[2 * I * K + ix] = 2*(2 * q3);
  grad_fs[3 * I * K + ix] = 2*(2 * q4 + q1);
  grad_fs[4 * I * K + ix] = 2*(2 * q5);  

      }
    }
    else
    {
      for(int ix = 0; ix < I * K; ix++)
      {
  q1 = (Qb[ix * 5 + 0] + s0/3) * x_b[ix] + Qb[ix * 5 + 1] * y_b[ix] + Qb[ix * 5 + 2] * z_b[ix];
  q2 = Qb[ix * 5 + 1] * x_b[ix] + (Qb[ix * 5 + 3] + s0/3) * y_b[ix] + Qb[ix * 5 + 4] * z_b[ix];
  q3 = Qb[ix * 5 + 2] * x_b[ix] + Qb[ix * 5 + 4] * y_b[ix] - (Qb[ix * 5 + 0] + Qb[ix * 5 + 3] - s0/3) * z_b[ix];

  fs[ix] = (q1 * q1 + q2 * q2 + q3 * q3); // / (R1 * R1);

  grad_fs[ix] = (2 * q1 * x_b[ix] - 2 * q3 * z_b[ix]); // / (R1 * R1);
  grad_fs[I * K + ix] = (2 * q1 * y_b[ix] + 2 * q2 * x_b[ix]); // / (R1 * R1);
  grad_fs[2 * I * K + ix] = (2 * q1 * z_b[ix]  + 2 * q3 * x_b[ix]); // / (R1 * R1);
  grad_fs[3 * I * K + ix] = (2 * q2 * y_b[ix] - 2 * q3 * z_b[ix]); // / (R1 * R1);
  grad_fs[4 * I * K + ix] = (2 * q2 * z_b[ix] + 2 * q3 * y_b[ix]); // / (R1 * R1);
      }
    }
}

}

double calc_energy_Fs(double *Q)
{
  double Fs = 0;
  double *fs1 = new double[I * K];
  double *fs2 = new double[I * K];
  double *Qb = new double[2 * I * K * 5];
  int ix, ix1;
  
  for(int bi = 0; bi < 2; bi++)
  {
    for(int i = 0; i < I; i++)
    {
      for(int k = 0; k < K; k++)
      {
	for (int n = 0; n < 5; n++)
	{
 	  Qb[bi * 5 * I * K + i * K * 5 + k * 5 + n] = Q[I * J * K + bi * I * K + i * K + k + n * Point];
	}
      }
    }
  }

  int Ind = 1;
  calc_fs(Qb, xb, yb, zb, fs1, grad_fs1, Ind);
  Ind  = 2;
  calc_fs(Qb + 5 * I * K, xb + I * K, yb + I * K, zb + I * K, fs2, grad_fs2, Ind);

  ix = 0;
  double fs1_tmp = 0, fs2_tmp = 0; 
  for (int i = 0; i < I; i++)
  {
    fs1_tmp = 0;
    fs2_tmp = 0;
    for (int k = 0; k < K; k++)
    {
      fs1_tmp += fs1[i * K + k];
      fs2_tmp += fs2[i * K + k];
    }     
    Fs += fs1_tmp * surface_jacobi[i] * coe_mu[i] + fs2_tmp * surface_jacobi[i + I] * coe_mu[i];
    // Fs += fs1_tmp * surface_jacobi[i] * coe_r[i];
  }

  Fs = Fs * 2 * PI / K;

  // cout << Fs << endl;

  double *fnm1, *fnm2;
  fnm1 = new double[(2 * M - 1) * N - M * (M - 1)]();
  fnm2 = new double[(2 * M - 1) * N - M * (M - 1)]();

  int sig1 = 0, sig2 = 1;

  for(int s = 0; s < 5; s++)
  {
    calc_gnm(grad_fs1 + s * I * K, fnm1, sig1);
    calc_gnm(grad_fs2 + s * I * K, fnm2, sig2);

    ix = 0;
    ix1 = 0;
      for(int l = 0; l < L; l++)
	{
	  ix1 = 0;
	  for(int m = 1 - M; m < M; m++) 
	  {
	  for(int n = abs(m); n < N; n++)
		{
		  grad_Fs[s * Basis + ix] = fnm1[ix1] * Plp[J * L + l] + fnm2[ix1] * Plp[(J + 1) * L + l];
		  // grad_Fs[s * Basis + ix] = fnm1[ix1] * Plp[J * L + l];
		  ix++;
		  ix1++;
		}
	  }
	}
  }

  delete[] Qb;
  delete[] fs1;
  delete[] fs2;
  delete[] fnm1;
  delete[] fnm2;

  return Fs;
}

