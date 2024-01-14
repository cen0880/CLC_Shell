/**
 * @file   main.cpp
 * @author Yucen Han
 * @date   Jan 14 2024
 * @brief main program modified from 1Dimension's nematic shell code
 *
 *
 */
using namespace std;

#include <iomanip>

#include <cstdio>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "lbfgs.h"

#include "global.h"
#include "Basis.h"
#include "initial.h"
#include "calculate_Fel.h"
#include "calculate_Fs.h"

#include "EigInfoTwo.h"
#include "in_and_out.h"
#include "calc_vector.h"
#include "point.h"

/*********************************************************************/
/**
 *
 * @param A  Zernike coefficients
 * @param V
 * @param Q
 * @param Felas Bulk_energy + Elastic_energy
 * @param Fs Surface Energy
 *
 *
 */
double cal_F(double *A, double *Q, double& Felas, double& Fs)
{

  Felas = calc_F_el(A, Q);
  Fs = calc_energy_Fs(Q);

  return (xi1 - xi0)/2 * Felas + eta * Fs;

}

/*********************************************************************/

void cal_dF(double *A, double *Q,double *grad_Energy)
{

  for(int i = 0; i < 5 * Basis; i++)
  {
    grad_Energy[i] = (xi1 - xi0)/2 * (grad_Fb[i] + grad_FeL[i] + grad_FeL1[i]) + eta * grad_Fs[i];
  }

}

static lbfgsfloatval_t evaluate(void *instance,
				const double *Bnlm,
				double *grad_Energy,
				const int n,
				const lbfgsfloatval_t step)
{
  double Energy;
  Energy = cal_F(Anlm, Qijk, Felas, Fs);
  cal_dF(Anlm, Qijk, grad_Energy);
  return Energy;
}

static int progress(void *instance,
		    double *Anlm,
 		    double *grad_Energy,
		    double Energy,
		    const lbfgsfloatval_t xnorm,
		    const lbfgsfloatval_t gnorm,
		    const lbfgsfloatval_t step,
		    int n,
		    int k,
		    int ls)
{

  if(k % 100 == 0)
  {
    printf("Iteration %d:  ",k);
    printf("Energy = %16.15f  ",Energy);
    printf("normdF = %16.15f  step = %16.15f \n Fs = %16.15f, Felas = %16.15f\n",gnorm,step, Fs, Felas);
  }



  if(k % 1000 == 0)
  {

    // sprintf(fname_loc, "%s_T_%d", fname_arv, k);

    // cout << fname_loc << endl;


    oput(Anlm, landau_t, Rad, eta, N, L, M, fname_arv);
    cout << "Update done!" << endl;

    // delete[] filename;
  }


  return 0;
}

/*********************************************************************/
int main(int argc,char *argv[])
{
  int i, j, k, n, l, s, m, ret = 0, con = 0;
  double Energy,normdF;
  double fr_Rad,fr_t,fr_eta,st_Rad,st_t,st_eta,ed_Rad,ed_t,ed_eta;
  double S0;
  lbfgs_parameter_t param;
  FILE *fp = fopen("Log.txt","a+");


  fname_arv = argv[1];
  fname_loc = new char[30];

  //设定物理参数

  double A = - 0.172 * 1e6, B = -2.12 * 1e6, C = 1.73 * 1e6, period = 5, K1 = 4 * 1e-11;
  R1 = 0.7;
  //double ep =1e-4; //symmetric shell
  double ep = 0.1; // asymmetric shell


  a = sqrt((pow(1.0 - R1*R1 - ep * ep, 2.0) - 4 * ep * ep * R1 * R1))/(2.0 * ep); //parameter in bispherical coordinate
  xi0 = asinh(a);
  xi1 = asinh(a/R1);

  cout << "rho = " <<  R1 << " c = " << ep << " a = " << a << " xi0 = " << xi0 << " xi1 = " << xi1 << endl;
 
  landau_t = 27 * A * C / (B * B);
  Rad = 50;
  etaL = 1.0; 
  q0 = 2.0*PI*period;
  eta = 1e-2; //penalty constant in surface energy
  cout << "t = " << landau_t << " Rad = " << Rad << " etaL = " <<  etaL << " q0 = " << q0 << " eta = " << eta << endl;
  cout << "R0 = " << sqrt(27*C*K1*Rad*Rad/B/B) << "m W = " << eta*B*B*sqrt(27*C*K1*Rad*Rad/B/B)/0.5/27/C << "J/m^3" << endl;
  
  Qscale = - 3.0*sqrt(6)/2 * C / B;
  beig = Qscale/3;


  S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;


  Basis_init(64, 32, 64, 256, 65, 256);
  malloc_variables_init();  //memory allocation and variable-definition

  cout << "Done!" << endl;

if(! thisAnderson.initialized)
  {
    thisAnderson.m_max = 30;
    thisAnderson.N = 5 * Basis;
    thisAnderson.Initialize();
  }

  iput(Anlm, landau_t, Rad, eta, 64, 32, 64, argv[1], 0, 3);  //set an initial condition in term of Zernike coefficient Anlm

  lbfgs_parameter_init(&param);
  param.m = 50;   ///the number of step of Hessen Matrix updating
  param.epsilon = 1e-15;
  param.delta = 1e-10;
  param.max_iterations = 100000;
  param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
  // param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
  // param.linesearch = LBFGS_LINESEARCH_MORETHUENTE;
  ret = lbfgs(5 * Basis, Anlm, &Energy, evaluate, progress, NULL, &param);   ///BFGS method to find minimizers


  oput(Anlm, landau_t, Rad, eta, N, L, M, fname_arv);
  Energy = cal_F(Anlm, Qijk, Felas, Fs);  ///calculating total energy
  cal_dF(Anlm, Qijk, grad_Energy);
  printf("R = %.2f t = %.2f Energy = %16.15f normdF = %16.15f ", Rad, landau_t, Energy, Norm(grad_Energy, 5*Basis));
  fprintf(fp,"%.2f %.2f %.6f, %.2f, %.2f, %.2f, %.10f, %s\n",landau_t, ep, Rad, eta, etaL, q0, Energy, argv[1]);

  fclose(fp);
  return 0;
}
