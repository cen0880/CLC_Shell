/*
 * @file   in_and_out.h
 * @author Yucen Han
 * @date   Jan 14 2024
 * 
 * @brief initial conditions, modified from 1Dimension's nematic shell code
 * case 1: n = ez
 * case 2: n = ex z>0 n = ey z<0
 * case 3: out homeotropic in planar
 * case 4: in homeotropic out planar
 * case 5: in and out homeotropic
 * case 6: planar state
 */
using namespace std;

#include <iostream>
#include <fstream>
#include "Basis.h"
#include "global.h"

void initial(double *A,int mode) 
{
  int i,j,k,n,ix;

  double Theta = PI/2.0; 

  double theta_b, phi_b, r_b;

  double x, y, z;
  double xn, yn, zn;

  int Boundary_index = I * J * K;
 
  double *Q = new double[5 * Point]();
  double S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;




  if(landau_t > 9.0/8)
  {  
    S0 = 0;
  }     

  cout << "S0 = "<< S0 << endl;
  
  switch (mode) 
  {
  case 0:

     
    srand((unsigned)time(NULL));
    for (i = 0; i < 5 * Basis; i++)
      A[i] = 0.2*(2.0*rand()/RAND_MAX - 1);
   

    break;

  case 1:  ///n = (0,0,1),Q = s0(n\otimes n - I/3)

       
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
	for (k = 0; k < K; k++) 
	{

	  Q[i * J * K + j * K + k] = S0 * (-1.0/3.0) ;
	  Q[Point + i * J * K + j * K + k] = 0;
	  Q[2 * Point + i * J * K + j * K + k] = 0;
	  Q[3 * Point + i * J * K + j * K + k] =  S0 * (- 1.0/3.0);
	  Q[4 * Point + i * J * K + j * K + k] = 0;

	  
	}
      }
    }
  
    for (n = 0; n < 5; n++)
    {
      calc_fnlm(Q + n * Point, A + n * Basis);
    }

    break;

  case 2://n = ex for z>0 and n = ey for z<0

  
    for (j = 0; j < J; j++)
    { 
      ix = 0;

      for (i = 0; i < I; i++)
      {
	for (k = 0; k < K; k++)
	{

	  x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
	  y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
	  z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i]))  - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));

	  xn = x;
	  yn = y;
	  zn = z;

	  if(z > 0)
	  {
	    Q[i * J * K + j * K + k] = S0 * (1 * 1 - 1.0/3);
	    Q[Point + i * J * K + j * K + k] = 0;
	    Q[2 * Point + i * J * K + j * K + k] = 0;
	    Q[3 * Point + i * J * K + j * K + k] = S0 * (- 1.0/3);
	    Q[4 * Point + i * J * K + j * K + k] = 0;
	    
	    
	  }
	  else
	  {
	    Q[i * J * K + j * K + k] = S0 * (0 - 1.0/3);
	    Q[Point + i * J * K + j * K + k] = 0;
	    Q[2 * Point + i * J * K + j * K + k] = 0;
	    Q[3 * Point + i * J * K + j * K + k] = S0 * (1 - 1.0/3);
	    Q[4 * Point + i * J * K + j * K + k] = 0;
	   }
	 

	  ix++; 
	}
      }
    }    
       
    for (n = 0; n < 5; n++)
    {
      calc_fnlm(Q + n * Point, A + n * Basis);
    }

    break;

  case 3://symmetric shell focal conic out homeotropic in planar
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
	for (k = 0; k < K; k++) 
	{
	    Q[i * J * K + j * K + k] = 0;
	    Q[Point + i * J * K + j * K + k] = 0;
	    Q[2 * Point + i * J * K + j * K + k] = 0;
	    Q[3 * Point + i * J * K + j * K + k] = 0;
	    Q[4 * Point + i * J * K + j * K + k] = 0;
    }
	}
      }
    int center_number;
    double center_phi[100], center_theta[100], theta_array[10];
    int center_shape[50];
    double x_u, y_u, z_u, r_u, phi_u, theta_u, theta0, phi0, theta1, phi1, nx, ny, nz, nx_tmp, ny_tmp, nz_tmp, dx_tmp, dy_tmp, dz_tmp, range, tmp0;
    center_number = 0;
    center_theta[center_number] = 0.0; center_phi[center_number] = 0.0; center_shape[center_number] = 5;
    center_number++;
    center_theta[center_number] = PI; center_phi[center_number] = 0.0; center_shape[center_number] = 5;
    center_number++;
    double len,alpha;
    len = PI/(3.0+2.0*sqrt(3.0) + 2.0*sin(54.0/180.0*PI));
    theta_array[1] = len*sin(54.0/180.0*PI) + len*sqrt(3.0)/2.0; theta_array[2] = 3.0*len-0.05; theta_array[3] = len*sin(54.0/180.0*PI) + len*sqrt(3.0)/2.0*3.0+0.05; 
    theta_array[4] = 3*len + len*sin(54.0/180.0*PI) + len*sqrt(3.0)/2.0-0.05; theta_array[5] = 2.0*len*sin(54.0/180.0*PI) + len*sqrt(3.0)*2.0+0.05; theta_array[6] = len*sin(54.0/180.0*PI) + len*sqrt(3.0)*3.0/2.0 + 3.0*len;
    for (i = 1; i < 7; i++) 
    {
        for (j = 0; j < 5; j++) 
        {
            center_theta[center_number] = theta_array[i];
            if (i % 2 == 0)
              {center_phi[center_number] = j*PI/5.0*2.0;}
            else
              {center_phi[center_number] = (j+0.5)*PI/5.0*2.0;}
            if (i == 1 || i == 3 || i == 4 || i == 6)
              {center_shape[center_number] = 6;}
            else
              {center_shape[center_number] = 5;}
            cout << "i = " << i << " j = " << j << " center number" << center_number << "theta = " << center_theta[center_number] << " phi = " << center_phi[center_number] << " shape = " << center_shape[center_number] << endl;
            center_number++;
        }
    }

    cout << center_number << endl;
    for (int u = 0; u<center_number; u++){
      if (center_shape[u] == 5)
      {
        //x_u = sin(center_theta[u])*cos(center_phi[u]);
        //y_u = sin(center_theta[u])*sin(center_phi[u]);
        //z_u = cos(center_theta[u]);
       // cout << center_theta[u] << " " << center_phi[u] << endl;
        x_u = a * sin(center_theta[u])/(cosh(realxi(p[J-1])) - cos(center_theta[u])) * cos(center_phi[u]);
	      y_u = a * sin(center_theta[u])/(cosh(realxi(p[J-1])) - cos(center_theta[u])) * sin(center_phi[u]);
	      z_u = a * sinh(realxi(p[J-1]))/(cosh(realxi(p[J-1])) - cos(center_theta[u])) - a * cosh(realxi(p[J-1]))/sinh(realxi(p[J-1]));
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
	for (k = 0; k < K; k++) 
	{
        x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
	      y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
	      z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i])) - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));
        r_u = sqrt((x-x_u)*(x-x_u)+(y-y_u)*(y-y_u)+(z-z_u)*(z-z_u));
        dx_tmp = (x-x_u)*cos(center_theta[u])*cos(center_phi[u]) + (y-y_u)*cos(center_theta[u])*sin(center_phi[u]) - (z-z_u)*sin(center_theta[u]);
        dy_tmp = -(x-x_u)*sin(center_phi[u]) + (y-y_u)*cos(center_phi[u]);
        dz_tmp = (x-x_u)*sin(center_theta[u])*cos(center_phi[u]) + (y-y_u)*sin(center_theta[u])*sin(center_phi[u]) + (z-z_u)*cos(center_theta[u]);
        phi_u = atan2(dy_tmp,dx_tmp);
        theta_u = abs(atan2(sqrt((dx_tmp)*(dx_tmp)+(dy_tmp)*(dy_tmp)),dz_tmp));
        phi0 = atan2(y,x);
        theta0 = atan2(sqrt(x*x+y*y),z);
        tmp0 = (x-x_u)*sin(theta0)*cos(phi0) + (y-y_u)*sin(theta0)*sin(phi0) + (z-z_u)*cos(theta0);
        range = sqrt(pow(x-x_u-tmp0*sin(theta0)*cos(phi0),2) + pow(y-y_u-tmp0*sin(theta0)*sin(phi0),2) + pow(z-z_u-tmp0*cos(theta0),2));
        //if (range < len && r_u < len*sqrt(2.0))
        if (range < len*sqrt(3.0)/2.0 && r_u < len*sqrt(3.0)/2.0*sqrt(2.0))
        //if (range < 0.3 && r_u < 0.3*sqrt(2.0))
        {
         //nx_tmp = cos(q0*r_u)*cos(theta_u)*cos(phi_u)-sin(q0*r_u)*sin(phi_u), ny_tmp = cos(q0*r_u)*cos(theta_u)*sin(phi_u)+sin(q0*r_u)*cos(phi_u), nz_tmp = -cos(q0*r_u)*sin(theta_u);
         nx_tmp = cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(theta_u)*cos(phi_u)-sin(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*sin(phi_u), ny_tmp = cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(theta_u)*sin(phi_u)+sin(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(phi_u), nz_tmp = -cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*sin(theta_u);
         //nx_tmp = sin(q0*r_u)*cos(theta_u)*cos(phi_u)-cos(q0*r_u)*sin(phi_u), ny_tmp = sin(q0*r_u)*cos(theta_u)*sin(phi_u)+cos(q0*r_u)*cos(phi_u), nz_tmp = -sin(q0*r_u)*sin(theta_u);
         nx = nx_tmp*cos(center_theta[u])*cos(center_phi[u]) - ny_tmp*sin(center_phi[u]) + nz_tmp*sin(center_theta[u])*cos(center_phi[u]);
         ny = nx_tmp*cos(center_theta[u])*sin(center_phi[u]) + ny_tmp*cos(center_phi[u]) + nz_tmp*sin(center_theta[u])*sin(center_phi[u]);
         nz = -nx_tmp*sin(center_theta[u]) + nz_tmp*cos(center_theta[u]);
         Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
         Q[Point + i * J * K + j * K + k] = S0 * (nx * ny);
         Q[2 * Point + i * J * K + j * K + k] = S0 * nx * nz;
         Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
         Q[4 * Point + i * J * K + j * K + k] = S0 * (ny * nz);
        }

	}
      }
    }
  }
    }
    for (int u = 0; u<center_number; u++){
    //for (int u = 16; u<17; u++){
      if (center_shape[u] == 6)
      {
        //x_u = sin(center_theta[u])*cos(center_phi[u]);
        //y_u = sin(center_theta[u])*sin(center_phi[u]);
        //z_u = cos(center_theta[u]);
       // cout << center_theta[u] << " " << center_phi[u] << endl;
        x_u = a * sin(center_theta[u])/(cosh(realxi(p[J-1])) - cos(center_theta[u])) * cos(center_phi[u]);
	    y_u = a * sin(center_theta[u])/(cosh(realxi(p[J-1])) - cos(center_theta[u])) * sin(center_phi[u]);
	    z_u = a * sinh(realxi(p[J-1]))/(cosh(realxi(p[J-1])) - cos(center_theta[u])) - a * cosh(realxi(p[J-1]))/sinh(realxi(p[J-1]));
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
  for (k = 0; k < K; k++) 
  {
        x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
        y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
        z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i])) - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));
        r_u = sqrt((x-x_u)*(x-x_u)+(y-y_u)*(y-y_u)+(z-z_u)*(z-z_u));
        dx_tmp = (x-x_u)*cos(center_theta[u])*cos(center_phi[u]) + (y-y_u)*cos(center_theta[u])*sin(center_phi[u]) - (z-z_u)*sin(center_theta[u]);
        dy_tmp = -(x-x_u)*sin(center_phi[u]) + (y-y_u)*cos(center_phi[u]);
        dz_tmp = (x-x_u)*sin(center_theta[u])*cos(center_phi[u]) + (y-y_u)*sin(center_theta[u])*sin(center_phi[u]) + (z-z_u)*cos(center_theta[u]);
        phi_u = atan2(dy_tmp,dx_tmp);
        theta_u = atan2(sqrt((dx_tmp)*(dx_tmp)+(dy_tmp)*(dy_tmp)),dz_tmp);
        phi0 = atan2(y,x);
        theta0 = atan2(sqrt(x*x+y*y),z);
        tmp0 = (x-x_u)*sin(theta0)*cos(phi0) + (y-y_u)*sin(theta0)*sin(phi0) + (z-z_u)*cos(theta0);
        range = sqrt(pow(x-x_u-tmp0*sin(theta0)*cos(phi0),2) + pow(y-y_u-tmp0*sin(theta0)*sin(phi0),2) + pow(z-z_u-tmp0*cos(theta0),2));
        
        //if (range < len && r_u < len*sqrt(2.0))
        if (range < len*sqrt(3.0)/2.0 && r_u < len*sqrt(3.0)/2.0*sqrt(2.0))
        //if (range < 0.3 && r_u < 0.3*sqrt(2.0))
        {
         //nx = cos(alpha + q0*r_u)*cos(theta_u)*cos(phi_u)-sin(alpha + q0*r_u)*sin(phi_u), ny = cos(alpha + q0*r_u)*cos(theta_u)*sin(phi_u)+sin(alpha + q0*r_u)*cos(phi_u), nz = -cos(alpha + q0*r_u)*sin(theta_u);
         nx_tmp = cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(theta_u)*cos(phi_u)-sin(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*sin(phi_u), ny_tmp = cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(theta_u)*sin(phi_u)+sin(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(phi_u), nz_tmp = -cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*sin(theta_u);
         //nx_tmp = sin(q0*r_u)*cos(theta_u)*cos(phi_u)-cos(q0*r_u)*sin(phi_u), ny_tmp = sin(q0*r_u)*cos(theta_u)*sin(phi_u)+cos(q0*r_u)*cos(phi_u), nz_tmp = -sin(q0*r_u)*sin(theta_u);
         nx = nx_tmp*cos(center_theta[u])*cos(center_phi[u]) - ny_tmp*sin(center_phi[u]) + nz_tmp*sin(center_theta[u])*cos(center_phi[u]);
         ny = nx_tmp*cos(center_theta[u])*sin(center_phi[u]) + ny_tmp*cos(center_phi[u]) + nz_tmp*sin(center_theta[u])*sin(center_phi[u]);
         nz = -nx_tmp*sin(center_theta[u]) + nz_tmp*cos(center_theta[u]);
         Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
         Q[Point + i * J * K + j * K + k] = S0 * (nx * ny);
         Q[2 * Point + i * J * K + j * K + k] = S0 * nx * nz;
         Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
         Q[4 * Point + i * J * K + j * K + k] = S0 * (ny * nz);
        }

  }
      }
    }
  }
    }
    //delete [] center_theta;
    //delete [] center_phi;
    for (n = 0; n < 5; n++)
    {
      calc_fnlm(Q + n * Point, A + n * Basis);
    }
    break;

    case 4://symmetric shell focal conic in homeotropic out planar
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
	for (k = 0; k < K; k++) 
	{
	    Q[i * J * K + j * K + k] = 0;
	    Q[Point + i * J * K + j * K + k] = 0;
	    Q[2 * Point + i * J * K + j * K + k] = 0;
	    Q[3 * Point + i * J * K + j * K + k] = 0;
	    Q[4 * Point + i * J * K + j * K + k] = 0;
    }
	}
      }
    center_number = 0;
    center_theta[center_number] = 0.0; center_phi[center_number] = 0.0; center_shape[center_number] = 5;
    center_number++;
    center_theta[center_number] = PI; center_phi[center_number] = 0.0; center_shape[center_number] = 5;
    center_number++;
    len = PI/(3.0+2.0*sqrt(3.0) + 2.0*sin(54.0/180.0*PI));
    theta_array[1] = len*sin(54.0/180.0*PI) + len*sqrt(3.0)/2.0; theta_array[2] = 3.0*len-0.05; theta_array[3] = len*sin(54.0/180.0*PI) + len*sqrt(3.0)/2.0*3.0+0.05; 
    theta_array[4] = 3*len + len*sin(54.0/180.0*PI) + len*sqrt(3.0)/2.0-0.05; theta_array[5] = 2.0*len*sin(54.0/180.0*PI) + len*sqrt(3.0)*2.0+0.05; theta_array[6] = len*sin(54.0/180.0*PI) + len*sqrt(3.0)*3.0/2.0 + 3.0*len;
    for (i = 1; i < 7; i++) 
    {
        for (j = 0; j < 5; j++) 
        {
            center_theta[center_number] = theta_array[i];
            if (i % 2 == 0)
              {center_phi[center_number] = j*PI/5.0*2.0;}
            else
              {center_phi[center_number] = (j+0.5)*PI/5.0*2.0;}
            if (i == 1 || i == 3 || i == 4 || i == 6)
              {center_shape[center_number] = 6;}
            else
              {center_shape[center_number] = 5;}
            cout << "i = " << i << " j = " << j << " center number" << center_number << "theta = " << center_theta[center_number] << " phi = " << center_phi[center_number] << " shape = " << center_shape[center_number] << endl;
            center_number++;
        }
    }

    cout << center_number << endl;
    for (int u = 0; u<center_number; u++){
      if (center_shape[u] == 5)
      {
        //x_u = sin(center_theta[u])*cos(center_phi[u]);
        //y_u = sin(center_theta[u])*sin(center_phi[u]);
        //z_u = cos(center_theta[u]);
       // cout << center_theta[u] << " " << center_phi[u] << endl;
        x_u = a * sin(center_theta[u])/(cosh(realxi(p[1])) - cos(center_theta[u])) * cos(center_phi[u]);
	    y_u = a * sin(center_theta[u])/(cosh(realxi(p[1])) - cos(center_theta[u])) * sin(center_phi[u]);
	    z_u = a * sinh(realxi(p[1]))/(cosh(realxi(p[1])) - cos(center_theta[u])) - a * cosh(realxi(p[1]))/sinh(realxi(p[1]));
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
	for (k = 0; k < K; k++) 
	{
        x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
	    y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
	    z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i])) - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));
        r_u = sqrt((x-x_u)*(x-x_u)+(y-y_u)*(y-y_u)+(z-z_u)*(z-z_u));
        dx_tmp = (x-x_u)*cos(center_theta[u])*cos(center_phi[u]) + (y-y_u)*cos(center_theta[u])*sin(center_phi[u]) - (z-z_u)*sin(center_theta[u]);
        dy_tmp = -(x-x_u)*sin(center_phi[u]) + (y-y_u)*cos(center_phi[u]);
        dz_tmp = (x-x_u)*sin(center_theta[u])*cos(center_phi[u]) + (y-y_u)*sin(center_theta[u])*sin(center_phi[u]) + (z-z_u)*cos(center_theta[u]);
        phi_u = atan2(dy_tmp,dx_tmp);
        theta_u = abs(atan2(sqrt((dx_tmp)*(dx_tmp)+(dy_tmp)*(dy_tmp)),dz_tmp));
        phi0 = atan2(y,x);
        theta0 = atan2(sqrt(x*x+y*y),z);
        tmp0 = (x-x_u)*sin(theta0)*cos(phi0) + (y-y_u)*sin(theta0)*sin(phi0) + (z-z_u)*cos(theta0);
        range = sqrt(pow(x-x_u-tmp0*sin(theta0)*cos(phi0),2) + pow(y-y_u-tmp0*sin(theta0)*sin(phi0),2) + pow(z-z_u-tmp0*cos(theta0),2));
        //if (range < len*R1 && r_u < len*sqrt(2.0)*R1)
        if (range < len*sqrt(3.0)/2.0*R1 && r_u < len*sqrt(3.0)/2.0*sqrt(2.0)*R1)
        //if (range < 0.3 && r_u < 0.3*sqrt(2.0))
        {
         nx_tmp = cos(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*cos(theta_u)*cos(phi_u)-sin(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*sin(phi_u), ny_tmp = cos(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*cos(theta_u)*sin(phi_u)+sin(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*cos(phi_u), nz_tmp = -cos(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*sin(theta_u);
         //nx_tmp = cos(q0*r_u)*cos(theta_u)*cos(phi_u)-sin(q0*r_u)*sin(phi_u), ny_tmp = cos(q0*r_u)*cos(theta_u)*sin(phi_u)+sin(q0*r_u)*cos(phi_u), nz_tmp = -cos(q0*r_u)*sin(theta_u);
         //nx_tmp = sin(q0*r_u)*cos(theta_u)*cos(phi_u)-cos(q0*r_u)*sin(phi_u), ny_tmp = sin(q0*r_u)*cos(theta_u)*sin(phi_u)+cos(q0*r_u)*cos(phi_u), nz_tmp = -sin(q0*r_u)*sin(theta_u);
         nx = nx_tmp*cos(center_theta[u])*cos(center_phi[u]) - ny_tmp*sin(center_phi[u]) + nz_tmp*sin(center_theta[u])*cos(center_phi[u]);
         ny = nx_tmp*cos(center_theta[u])*sin(center_phi[u]) + ny_tmp*cos(center_phi[u]) + nz_tmp*sin(center_theta[u])*sin(center_phi[u]);
         nz = -nx_tmp*sin(center_theta[u]) + nz_tmp*cos(center_theta[u]);
         Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
         Q[Point + i * J * K + j * K + k] = S0 * (nx * ny);
         Q[2 * Point + i * J * K + j * K + k] = S0 * nx * nz;
         Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
         Q[4 * Point + i * J * K + j * K + k] = S0 * (ny * nz);
        }

	}
      }
    }
  }
    }
    for (int u = 0; u<center_number; u++){
    //for (int u = 16; u<17; u++){
      if (center_shape[u] == 6)
      {
        //x_u = sin(center_theta[u])*cos(center_phi[u]);
        //y_u = sin(center_theta[u])*sin(center_phi[u]);
        //z_u = cos(center_theta[u]);
       // cout << center_theta[u] << " " << center_phi[u] << endl;
        x_u = a * sin(center_theta[u])/(cosh(realxi(p[1])) - cos(center_theta[u])) * cos(center_phi[u]);
	    y_u = a * sin(center_theta[u])/(cosh(realxi(p[1])) - cos(center_theta[u])) * sin(center_phi[u]);
	    z_u = a * sinh(realxi(p[1]))/(cosh(realxi(p[1])) - cos(center_theta[u])) - a * cosh(realxi(p[1]))/sinh(realxi(p[1]));
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
  for (k = 0; k < K; k++) 
  {
        x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
        y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
        z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i])) - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));
        r_u = sqrt((x-x_u)*(x-x_u)+(y-y_u)*(y-y_u)+(z-z_u)*(z-z_u));
        dx_tmp = (x-x_u)*cos(center_theta[u])*cos(center_phi[u]) + (y-y_u)*cos(center_theta[u])*sin(center_phi[u]) - (z-z_u)*sin(center_theta[u]);
        dy_tmp = -(x-x_u)*sin(center_phi[u]) + (y-y_u)*cos(center_phi[u]);
        dz_tmp = (x-x_u)*sin(center_theta[u])*cos(center_phi[u]) + (y-y_u)*sin(center_theta[u])*sin(center_phi[u]) + (z-z_u)*cos(center_theta[u]);
        phi_u = atan2(dy_tmp,dx_tmp);
        theta_u = atan2(sqrt((dx_tmp)*(dx_tmp)+(dy_tmp)*(dy_tmp)),dz_tmp);
        phi0 = atan2(y,x);
        theta0 = atan2(sqrt(x*x+y*y),z);
        tmp0 = (x-x_u)*sin(theta0)*cos(phi0) + (y-y_u)*sin(theta0)*sin(phi0) + (z-z_u)*cos(theta0);
        range = sqrt(pow(x-x_u-tmp0*sin(theta0)*cos(phi0),2) + pow(y-y_u-tmp0*sin(theta0)*sin(phi0),2) + pow(z-z_u-tmp0*cos(theta0),2));
        
        //if (range < len*R1 && r_u < len*sqrt(2.0)*R1)
        if (range < len*sqrt(3.0)/2.0*R1 && r_u < len*sqrt(3.0)/2.0*sqrt(2.0)*R1)
        //if (range < 0.3 && r_u < 0.3*sqrt(2.0))
        {
         //nx = cos(alpha + q0*r_u)*cos(theta_u)*cos(phi_u)-sin(alpha + q0*r_u)*sin(phi_u), ny = cos(alpha + q0*r_u)*cos(theta_u)*sin(phi_u)+sin(alpha + q0*r_u)*cos(phi_u), nz = -cos(alpha + q0*r_u)*sin(theta_u);
         nx_tmp = cos(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*cos(theta_u)*cos(phi_u)-sin(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*sin(phi_u), ny_tmp = cos(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*cos(theta_u)*sin(phi_u)+sin(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*cos(phi_u), nz_tmp = -cos(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*sin(theta_u);
         //nx_tmp = cos(q0*r_u)*cos(theta_u)*cos(phi_u)-sin(q0*r_u)*sin(phi_u), ny_tmp = cos(q0*r_u)*cos(theta_u)*sin(phi_u)+sin(q0*r_u)*cos(phi_u), nz_tmp = -cos(q0*r_u)*sin(theta_u);
         //nx_tmp = sin(q0*r_u)*cos(theta_u)*cos(phi_u)-cos(q0*r_u)*sin(phi_u), ny_tmp = sin(q0*r_u)*cos(theta_u)*sin(phi_u)+cos(q0*r_u)*cos(phi_u), nz_tmp = -sin(q0*r_u)*sin(theta_u);
         nx = nx_tmp*cos(center_theta[u])*cos(center_phi[u]) - ny_tmp*sin(center_phi[u]) + nz_tmp*sin(center_theta[u])*cos(center_phi[u]);
         ny = nx_tmp*cos(center_theta[u])*sin(center_phi[u]) + ny_tmp*cos(center_phi[u]) + nz_tmp*sin(center_theta[u])*sin(center_phi[u]);
         nz = -nx_tmp*sin(center_theta[u]) + nz_tmp*cos(center_theta[u]);
         Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
         Q[Point + i * J * K + j * K + k] = S0 * (nx * ny);
         Q[2 * Point + i * J * K + j * K + k] = S0 * nx * nz;
         Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
         Q[4 * Point + i * J * K + j * K + k] = S0 * (ny * nz);
        }

  }
      }
    }
  }
    }
    //delete [] center_theta;
    //delete [] center_phi;
    for (n = 0; n < 5; n++)
    {
      calc_fnlm(Q + n * Point, A + n * Basis);
    }
    break;

    case 5://symmetric shell focal conic in and out homeotropic

    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
  for (k = 0; k < K; k++) 
  {
      Q[i * J * K + j * K + k] = 0;
      Q[Point + i * J * K + j * K + k] = 0;
      Q[2 * Point + i * J * K + j * K + k] = 0;
      Q[3 * Point + i * J * K + j * K + k] = 0;
      Q[4 * Point + i * J * K + j * K + k] = 0;
    }
  }
      }
  center_number = 0;
    center_theta[center_number] = 0.0; center_phi[center_number] = 0.0; center_shape[center_number] = 5;
    center_number++;
    center_theta[center_number] = PI; center_phi[center_number] = 0.0; center_shape[center_number] = 5;
    center_number++;
    len = PI/(3.0+2.0*sqrt(3.0) + 2.0*sin(54.0/180.0*PI));
    theta_array[1] = len*sin(54.0/180.0*PI) + len*sqrt(3.0)/2.0; theta_array[2] = 3.0*len-0.05; theta_array[3] = len*sin(54.0/180.0*PI) + len*sqrt(3.0)/2.0*3.0+0.05; 
    theta_array[4] = 3*len + len*sin(54.0/180.0*PI) + len*sqrt(3.0)/2.0-0.05; theta_array[5] = 2.0*len*sin(54.0/180.0*PI) + len*sqrt(3.0)*2.0+0.05; theta_array[6] = len*sin(54.0/180.0*PI) + len*sqrt(3.0)*3.0/2.0 + 3.0*len;
    for (i = 1; i < 7; i++) 
    {
        for (j = 0; j < 5; j++) 
        {
            center_theta[center_number] = theta_array[i];
            if (i % 2 == 0)
              {center_phi[center_number] = j*PI/5.0*2.0;}
            else
              {center_phi[center_number] = (j+0.5)*PI/5.0*2.0;}
            if (i == 1 || i == 3 || i == 4 || i == 6)
              {center_shape[center_number] = 6;}
            else
              {center_shape[center_number] = 5;}
            cout << "i = " << i << " j = " << j << " center number" << center_number << "theta = " << center_theta[center_number] << " phi = " << center_phi[center_number] << " shape = " << center_shape[center_number] << endl;
            center_number++;
        }
    }

    cout << center_number << endl;
    for (int u = 0; u<center_number; u++){
      if (center_shape[u] == 5)
      {
        x_u = a * sin(center_theta[u])/(cosh(realxi(p[J-1])) - cos(center_theta[u])) * cos(center_phi[u]);
      y_u = a * sin(center_theta[u])/(cosh(realxi(p[J-1])) - cos(center_theta[u])) * sin(center_phi[u]);
      z_u = a * sinh(realxi(p[J-1]))/(cosh(realxi(p[J-1])) - cos(center_theta[u])) - a * cosh(realxi(p[J-1]))/sinh(realxi(p[J-1]));
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
  for (k = 0; k < K; k++) 
  {
        x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
      y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
      z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i])) - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));
        r_u = sqrt((x-x_u)*(x-x_u)+(y-y_u)*(y-y_u)+(z-z_u)*(z-z_u));
        dx_tmp = (x-x_u)*cos(center_theta[u])*cos(center_phi[u]) + (y-y_u)*cos(center_theta[u])*sin(center_phi[u]) - (z-z_u)*sin(center_theta[u]);
        dy_tmp = -(x-x_u)*sin(center_phi[u]) + (y-y_u)*cos(center_phi[u]);
        dz_tmp = (x-x_u)*sin(center_theta[u])*cos(center_phi[u]) + (y-y_u)*sin(center_theta[u])*sin(center_phi[u]) + (z-z_u)*cos(center_theta[u]);
        phi_u = atan2(dy_tmp,dx_tmp);
        theta_u = abs(atan2(sqrt((dx_tmp)*(dx_tmp)+(dy_tmp)*(dy_tmp)),dz_tmp));
        phi0 = atan2(y,x);
        theta0 = atan2(sqrt(x*x+y*y),z);
        tmp0 = (x-x_u)*sin(theta0)*cos(phi0) + (y-y_u)*sin(theta0)*sin(phi0) + (z-z_u)*cos(theta0);
        range = sqrt(pow(x-x_u-tmp0*sin(theta0)*cos(phi0),2) + pow(y-y_u-tmp0*sin(theta0)*sin(phi0),2) + pow(z-z_u-tmp0*cos(theta0),2));
        //if (r_u < len*sqrt(3.0)/2.0) 
        if (range < len/2.0*sqrt(3.0)*0.8 && r_u < len*sqrt(2.0)/2.0*sqrt(3.0)*0.8)
        {
         //nx_tmp = cos(q0*r_u)*cos(theta_u)*cos(phi_u)-sin(q0*r_u)*sin(phi_u), ny_tmp = cos(q0*r_u)*cos(theta_u)*sin(phi_u)+sin(q0*r_u)*cos(phi_u), nz_tmp = -cos(q0*r_u)*sin(theta_u);
         nx_tmp = cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(theta_u)*cos(phi_u)-sin(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*sin(phi_u), ny_tmp = cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(theta_u)*sin(phi_u)+sin(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(phi_u), nz_tmp = -cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*sin(theta_u);
         //nx_tmp = sin(q0*r_u)*cos(theta_u)*cos(phi_u)-cos(q0*r_u)*sin(phi_u), ny_tmp = sin(q0*r_u)*cos(theta_u)*sin(phi_u)+cos(q0*r_u)*cos(phi_u), nz_tmp = -sin(q0*r_u)*sin(theta_u);
         nx = nx_tmp*cos(center_theta[u])*cos(center_phi[u]) - ny_tmp*sin(center_phi[u]) + nz_tmp*sin(center_theta[u])*cos(center_phi[u]);
         ny = nx_tmp*cos(center_theta[u])*sin(center_phi[u]) + ny_tmp*cos(center_phi[u]) + nz_tmp*sin(center_theta[u])*sin(center_phi[u]);
         nz = -nx_tmp*sin(center_theta[u]) + nz_tmp*cos(center_theta[u]);
         Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
         Q[Point + i * J * K + j * K + k] = S0 * (nx * ny);
         Q[2 * Point + i * J * K + j * K + k] = S0 * nx * nz;
         Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
         Q[4 * Point + i * J * K + j * K + k] = S0 * (ny * nz);
        }

  }
      }
    }
  }
    }
    for (int u = 0; u<center_number; u++){
      if (center_shape[u] == 6)
      {
        x_u = a * sin(center_theta[u])/(cosh(realxi(p[J-1])) - cos(center_theta[u])) * cos(center_phi[u]);
      y_u = a * sin(center_theta[u])/(cosh(realxi(p[J-1])) - cos(center_theta[u])) * sin(center_phi[u]);
      z_u = a * sinh(realxi(p[J-1]))/(cosh(realxi(p[J-1])) - cos(center_theta[u])) - a * cosh(realxi(p[J-1]))/sinh(realxi(p[J-1]));
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
  for (k = 0; k < K; k++) 
  {
        x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
        y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
        z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i])) - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));
        r_u = sqrt((x-x_u)*(x-x_u)+(y-y_u)*(y-y_u)+(z-z_u)*(z-z_u));
        dx_tmp = (x-x_u)*cos(center_theta[u])*cos(center_phi[u]) + (y-y_u)*cos(center_theta[u])*sin(center_phi[u]) - (z-z_u)*sin(center_theta[u]);
        dy_tmp = -(x-x_u)*sin(center_phi[u]) + (y-y_u)*cos(center_phi[u]);
        dz_tmp = (x-x_u)*sin(center_theta[u])*cos(center_phi[u]) + (y-y_u)*sin(center_theta[u])*sin(center_phi[u]) + (z-z_u)*cos(center_theta[u]);
        phi_u = atan2(dy_tmp,dx_tmp);
        theta_u = atan2(sqrt((dx_tmp)*(dx_tmp)+(dy_tmp)*(dy_tmp)),dz_tmp);
        phi0 = atan2(y,x);
        theta0 = atan2(sqrt(x*x+y*y),z);
        tmp0 = (x-x_u)*sin(theta0)*cos(phi0) + (y-y_u)*sin(theta0)*sin(phi0) + (z-z_u)*cos(theta0);
        range = sqrt(pow(x-x_u-tmp0*sin(theta0)*cos(phi0),2) + pow(y-y_u-tmp0*sin(theta0)*sin(phi0),2) + pow(z-z_u-tmp0*cos(theta0),2));
        
        //if (r_u < len)
        //if (r_u < len/2.0*sqrt(3.0))
        if (range < len/2.0*sqrt(3.0)*0.8 && r_u < len*sqrt(2.0)/2.0*sqrt(3.0)*0.8)
        {
         //nx = cos(alpha + q0*r_u)*cos(theta_u)*cos(phi_u)-sin(alpha + q0*r_u)*sin(phi_u), ny = cos(alpha + q0*r_u)*cos(theta_u)*sin(phi_u)+sin(alpha + q0*r_u)*cos(phi_u), nz = -cos(alpha + q0*r_u)*sin(theta_u);
         //nx_tmp = cos(q0*r_u)*cos(theta_u)*cos(phi_u)-sin(q0*r_u)*sin(phi_u), ny_tmp = cos(q0*r_u)*cos(theta_u)*sin(phi_u)+sin(q0*r_u)*cos(phi_u), nz_tmp = -cos(q0*r_u)*sin(theta_u);
         nx_tmp = cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(theta_u)*cos(phi_u)-sin(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*sin(phi_u), ny_tmp = cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(theta_u)*sin(phi_u)+sin(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*cos(phi_u), nz_tmp = -cos(q0*r_u*2.4/(cosh(1.60944)-cos(center_theta[u])))*sin(theta_u);
         //nx_tmp = sin(q0*r_u)*cos(theta_u)*cos(phi_u)-cos(q0*r_u)*sin(phi_u), ny_tmp = sin(q0*r_u)*cos(theta_u)*sin(phi_u)+cos(q0*r_u)*cos(phi_u), nz_tmp = -sin(q0*r_u)*sin(theta_u);
         nx = nx_tmp*cos(center_theta[u])*cos(center_phi[u]) - ny_tmp*sin(center_phi[u]) + nz_tmp*sin(center_theta[u])*cos(center_phi[u]);
         ny = nx_tmp*cos(center_theta[u])*sin(center_phi[u]) + ny_tmp*cos(center_phi[u]) + nz_tmp*sin(center_theta[u])*sin(center_phi[u]);
         nz = -nx_tmp*sin(center_theta[u]) + nz_tmp*cos(center_theta[u]);
         Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
         Q[Point + i * J * K + j * K + k] = S0 * (nx * ny);
         Q[2 * Point + i * J * K + j * K + k] = S0 * nx * nz;
         Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
         Q[4 * Point + i * J * K + j * K + k] = S0 * (ny * nz);
        }

  }
      }
    }
  }
    }
  center_number = 0;
  len_i = PI/(1.0/sin(36.0/180.0*PI) + 1.0/sin(36.0/180.0*PI)*sin(54.0/180.0*PI) + 1.0 + 2.0*sqrt(3.0));
    theta_array[1] = len_i/2.0/sin(36.0/180.0*PI); theta_array[2] = len_i/2.0/sin(36.0/180.0*PI) + len_i; theta_array[3] = len_i/sin(36.0/180.0*PI) + len_i-0.1; 
    theta_array[4] = len_i/sin(36.0/180.0*PI) + len_i + len_i/2.0/sin(36.0/180.0*PI)*sin(54.0/180.0*PI); theta_array[5] = len_i/sin(36.0/180.0*PI) + len_i + len_i/2.0/sin(36.0/180.0*PI)*sin(54.0/180.0*PI) + sqrt(3.0)/2.0*len_i; 
    theta_array[6] = len_i/sin(36.0/180.0*PI)*sin(54.0/180.0*PI) + sqrt(3.0)*2.0*len_i+0.1;
    theta_array[7] = len_i/2.0/sin(36.0/180.0*PI) + len_i/sin(36.0/180.0*PI)*sin(54.0/180.0*PI) + sqrt(3.0)*2.0*len_i;
    theta_array[8] = len_i/2.0/sin(36.0/180.0*PI) + len_i + len_i/sin(36.0/180.0*PI)*sin(54.0/180.0*PI) + sqrt(3.0)*2.0*len_i;
    inv = 2*PI/(5*(len_i+2*cos(54.0/180.0*PI)*len_i));
    phi_array[0] = 1.5*PI/20.0; phi_array[1] = 6.5*PI/20.0; phi_array[2] = 9.5*PI/20.0; phi_array[3] = 14.5*PI/20.0; phi_array[4] = 17.5*PI/20.0;
    phi_array[5] = 22.5*PI/20.0; phi_array[6] = 25.5*PI/20.0; phi_array[7] = 30.5*PI/20.0; phi_array[8] = 33.5*PI/20.0; phi_array[9] = 38.4*PI/20.0;
    
    phi_array_4[0] = 1.5*PI/20.0-0.03; phi_array_4[1] = 6.5*PI/20.0+0.03; phi_array_4[2] = 9.5*PI/20.0-0.03; phi_array_4[3] = 14.5*PI/20.0+0.03; phi_array_4[4] = 17.5*PI/20.0-0.03;
    phi_array_4[5] = 22.5*PI/20.0+0.03; phi_array_4[6] = 25.5*PI/20.0-0.03; phi_array_4[7] = 30.5*PI/20.0+0.03; phi_array_4[8] = 33.5*PI/20.0-0.03; phi_array_4[9] = 38.4*PI/20.0+0.03;
    for (i = 1; i < 9; i++) 
    {
      if (i == 1 || i == 2)
      {
        for (j = 0; j < 5; j++) 
        {
            center_theta[center_number] = theta_array[i];
            center_phi[center_number] = j*PI*2.0/5.0;
            cout << " center number" << center_number << "theta = " << center_theta[center_number] << " phi = " << center_phi[center_number] << endl;
            center_number++;
        }
      }
      if (i == 3)
      {
        for (j = 0; j < 10; j++) 
        {
            center_theta[center_number] = theta_array[i];
            center_phi[center_number] = phi_array[j] + PI/5.0;
            cout << " center number" << center_number << "theta = " << center_theta[center_number] << " phi = " << center_phi[center_number] << endl;
            center_number++;
        }
      }
      if (i == 4)
      {
        for (j = 0; j < 10; j++) 
        {
            center_theta[center_number] = theta_array[i];
            center_phi[center_number] = phi_array_4[j];
            cout << " center number" << center_number << "theta = " << center_theta[center_number] << " phi = " << center_phi[center_number] << endl;
            center_number++;
        }
      }
      if (i == 5)
      {
        for (j = 0; j < 10; j++) 
        {
            center_theta[center_number] = theta_array[i];
            center_phi[center_number] = phi_array_4[j]-PI/5.0;
            cout << " center number" << center_number << "theta = " << center_theta[center_number] << " phi = " << center_phi[center_number] << endl;
            center_number++;
        }
      }
      if (i == 6)
      {
        for (j = 0; j < 10; j++) 
        {
            center_theta[center_number] = theta_array[i];
            center_phi[center_number] = phi_array[j];
            cout << " center number" << center_number << "theta = " << center_theta[center_number] << " phi = " << center_phi[center_number] << endl;
            center_number++;
        }
      }
      if (i == 7 || i == 8)
      {
        for (j = 0; j < 5; j++) 
        {
            center_theta[center_number] = theta_array[i];
            center_phi[center_number] = (j + 0.5)*PI*2.0/5.0;
            cout << " center number" << center_number << "theta = " << center_theta[center_number] << " phi = " << center_phi[center_number] << endl;
            center_number++;
        }
      }
    }

    cout << center_number << endl;
    for (int u = 0; u<center_number; u++)
    //for (int u = 0; u<1; u++)
      {
      x_u = a * sin(center_theta[u])/(cosh(realxi(p[0])) - cos(center_theta[u])) * cos(center_phi[u]);
      y_u = a * sin(center_theta[u])/(cosh(realxi(p[0])) - cos(center_theta[u])) * sin(center_phi[u]);
      z_u = a * sinh(realxi(p[0]))/(cosh(realxi(p[0])) - cos(center_theta[u])) - a * cosh(realxi(p[0]))/sinh(realxi(p[0]));
    //cout << " x = " << x_u << " y = " << y_u << " z = " << z_u << endl;
    
      for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
        for (k = 0; k < K; k++) 
        {
        x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
        y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
        z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i])) - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));
        r_u = sqrt((x-x_u)*(x-x_u)+(y-y_u)*(y-y_u)+(z-z_u)*(z-z_u));
        dx_tmp = (x-x_u)*cos(center_theta[u])*cos(center_phi[u]) + (y-y_u)*cos(center_theta[u])*sin(center_phi[u]) - (z-z_u)*sin(center_theta[u]);
        dy_tmp = -(x-x_u)*sin(center_phi[u]) + (y-y_u)*cos(center_phi[u]);
        dz_tmp = (x-x_u)*sin(center_theta[u])*cos(center_phi[u]) + (y-y_u)*sin(center_theta[u])*sin(center_phi[u]) + (z-z_u)*cos(center_theta[u]);
        phi_u = atan2(dy_tmp,dx_tmp);
        theta_u = atan2(sqrt((dx_tmp)*(dx_tmp)+(dy_tmp)*(dy_tmp)),dz_tmp);
        phi0 = atan2(y,x);
        theta0 = atan2(sqrt(x*x+y*y),z);
        //if (r_u < len_i*3.0/8.0)
        if (r_u < len_i*1.0/4.0)
        {
         //nx_tmp = cos(q0*r_u)*cos(theta_u)*cos(phi_u)-sin(q0*r_u)*sin(phi_u), ny_tmp = cos(q0*r_u)*cos(theta_u)*sin(phi_u)+sin(q0*r_u)*cos(phi_u), nz_tmp = -cos(q0*r_u)*sin(theta_u);
         nx_tmp = cos(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*cos(theta_u)*cos(phi_u)-sin(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*sin(phi_u), ny_tmp = cos(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*cos(theta_u)*sin(phi_u)+sin(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*cos(phi_u), nz_tmp = -cos(q0*r_u*2.4/(cosh(1.94591)-cos(center_theta[u])))*sin(theta_u);
         nx = nx_tmp*cos(center_theta[u])*cos(center_phi[u]) - ny_tmp*sin(center_phi[u]) + nz_tmp*sin(center_theta[u])*cos(center_phi[u]);
         ny = nx_tmp*cos(center_theta[u])*sin(center_phi[u]) + ny_tmp*cos(center_phi[u]) + nz_tmp*sin(center_theta[u])*sin(center_phi[u]);
         nz = -nx_tmp*sin(center_theta[u]) + nz_tmp*cos(center_theta[u]);
         Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
         Q[Point + i * J * K + j * K + k] = S0 * (nx * ny);
         Q[2 * Point + i * J * K + j * K + k] = S0 * nx * nz;
         Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
         Q[4 * Point + i * J * K + j * K + k] = S0 * (ny * nz);
        }
        }
      }
    }
      }
    
    for (n = 0; n < 5; n++)
    {
      calc_fnlm(Q + n * Point, A + n * Basis);
    }
    break;

    
    case 6://planar state
    double r0;
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
  for (k = 0; k < K; k++) 
  {
        x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
        y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
        z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i])) - a * cosh(realxi(p[J-1]))/sinh(realxi(p[J-1]));
        r0 = sqrt(x*x + y*y + z*z);
        phi0 = atan2(y,x);
        theta0 = atan2(sqrt(x*x+y*y),z);
        nx = sin(phi0+q0*(1.0-r0))*(cos(theta0))*cos(phi0)-cos(phi0+q0*(1.0-r0))*sin(phi0), ny = sin(phi0+q0*(1.0-r0))*(cos(theta0))*sin(phi0)+cos(phi0+q0*(1.0-r0))*cos(phi0), nz = -sin(phi0+q0*(1.0-r0))*sin(theta0);
         Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
         Q[Point + i * J * K + j * K + k] = S0 * (nx * ny);
         Q[2 * Point + i * J * K + j * K + k] = S0 * nx * nz;
         Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
         Q[4 * Point + i * J * K + j * K + k] = S0 * (ny * nz);
        }
      }
  }
     
        
    for (n = 0; n < 5; n++)
    {
      calc_fnlm(Q + n * Point, A + n * Basis);
    }
    break;
    }
  //   delete [] s;
  delete [] Q;
}

/** 
 * @param fname 
 * @param t 
 * @param R 
 * @param eta 
 * @param N 
 * @param L 
 * @param M 
 * @param num 
 */
void getfname(char fname[],double t,double R,double eta,int N,int L,int M,char num[]) 
{
  char strR[20];
  char strN[20];
  char strL[20];
  char strM[20];
  
  if (R < 10)
    sprintf(strR,"0%.2f",R);
  else
    sprintf(strR,"%.2f",R);

  if (N < 10)
    sprintf(strN,"0%d",N);
  else
    sprintf(strN,"%d",N);

  if (L < 10)
    sprintf(strL,"0%d",L);
  else
    sprintf(strL,"%d",L);
  
  if (M < 10)
    sprintf(strM,"0%d",M);
  else
    sprintf(strM,"%d",M);
  
  sprintf(fname,"%s/Result/%s%s_t_%+.2f_R_%s_N_%s_L_%s_M_%s_eta_%.2f_etaL_%+.1f_q0_%+.1f_s0_%.2f_%s.txt",DIR,func_type,boundary,t,strR,strN,strL,strM,eta,etaL,q0,beig,num);

  //  sprintf(fname,"%s/Result/%s%s_t_%+.2f_R_%s_N_%s_L_%s_M_%s_eta_%.0e_%s.txt",DIR,func_type,boundary,t,strR,strN,strL,strM,eta,num);

}

/***************************************************************************/
/** 
 * 
 * 
 * @param A expansion coefficient
 * @param t reduced temperature
 * @param R radius
 * @param eta parameters
 * @param FN 
 * @param FL 
 * @param FM 
 * @param num 
 * @param firstvalue 0 creat an initial condition, otherwise read an initial condition from a file
 * @param mode 
 */
void iput(double *A, double t, double R,double eta,int FN,int FL,int FM,char num[],int firstvalue,int mode) 
{
  int i,n,l,m,ix,err;
  int Max;
  FILE *fp;
  char path[200];

  if (firstvalue == 0)
  {
    initial(A, mode);
  }

  else  //read expansion coefficient from file
  {
    
    if (M >= FM)
      Max = M;
    else
      Max = FM;
    
    Max = 128;

    double *DataAnlm = new double[5 * Max * Max * (2 * Max - 1)]();
    getfname(path,t,R,eta,FN,FL,FM,num);

    cout << path << endl;

    if ((fp = fopen(path,"r")) == NULL) printf("Open Anlm Error!\n");
    /**
    ix = 0;
    for(i = 0; i < 5 * Basis; i++)
    {
      err = fscanf(fp,"%lf\n", &A[ix]);
      ix++;
    }
    fclose(fp);
    */

    
    for(i = 0; i < 5; i++)
    {
      for(l = 0; l < FL; l++)
      {
	for(m = 1 - FM; m <= FM - 1; m++)
	{
	  for(n = abs(m); n < FN; n++)
	  {
	    err = fscanf(fp,"%lf\n",&DataAnlm[i * Max * Max * (2 * Max - 1) + l * (2 * Max -1) * Max + (m + Max - 1) * Max + n]);
	  }
	}
      }
    }
	 
   fclose(fp);

    ix = 0;
    for (i = 0; i < 5; i++)
    {
      for(l = 0; l < L; l++)
      {
	for(m = 1 - M; m <= M - 1; m++)
	{
	  for(n = abs(m); n < N; n++)
	  {
	    A[ix++] = DataAnlm[i * Max * Max * (2 * Max - 1) + l * (2 * Max -1) * Max + (m + Max - 1) * Max + n];
	  }
	}
      }
    }

    delete [] DataAnlm;
  }

}

/***************************************************************************/


void oput(double *A,double t,double R,double eta,int FN,int FL,int FM,char num[]) 
{
  int i,n,l,m,ix,err;
  FILE* fp;
  char path[200];
  
  getfname(path,t,R,eta,FN,FL,FM,num);
  fp = fopen(path,"w");  //path = ./Result/......，Directory “Result” needed

  for (int i = 0; i < 5 * Basis; i++)
  {
    fprintf(fp,"%16.15e\n",A[i]);
  }

  fclose(fp);
}



