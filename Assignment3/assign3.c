#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double gam = 1.;
double omega = 30.0;
int N_x = 2;
int N_y = 2;

////////////Physics///////
double  exact(double t){
  double om = sqrt(omega*omega - 0.25*gam*gam);
  double ans = exp(-0.5*gam*t)*(cos(om*t) + 0.5*gam/om*sin(om*t));
  //printf("%e %e\n", t, ans);
  return(ans);
}
void osc(double * F){
  F[0] = 0; //F00
  F[2] = 1; //F01
  F[1] = -pow(omega,2); //F10
  F[3] = -gam; //F11
}
////////////local functions/////////////
void mult2x2(double * F, double* a){
  double a_n[2] = {a[0], a[1]};
  a_n[0] = F[0]*a[0] + F[2]*a[1];
  a_n[1] = F[1]*a[0] + F[3]*a[1];
  a[0] = a_n[0];
  a[1] = a_n[1];
}
void invert2x2( double *F, double * Fi){
  double det = F[0]*F[3] - F[1]*F[2];
  Fi[0] =  F[3]/det;
  Fi[1] = -F[1]/det;
  Fi[2] = -F[2]/det;
  Fi[3] =  F[0]/det;
}
void inverse33( double ** a , double ** inv ){

   int i;
   double det = 0.0; 
   for( i=0 ; i<3 ; ++i )
      det += a[0][i]*(a[1][(i+1)%3]*a[2][(i+2)%3] - a[1][(i+2)%3]*a[2][(i+1)%3]) ;
 
   int j;
   for(i=0;i<3;i++){
      for(j=0;j<3;j++)
           inv[j][i] = ((a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]) - (a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3]))/det ;
   }
}
void Feuler2x2( double dt, double * z, double * F){
  double z_n[2] = {z[0], z[1]}; 
  mult2x2(F, z_n);
  z[0] += dt*z_n[0];
  z[1] += dt*z_n[1];
}
void Mid2x2(double dt, double * z, double * F){
  double z_half[2] = {z[0], z[1]};
  Feuler2x2(0.5*dt, z_half, F);
  mult2x2(F, z_half);
  z[0] += dt*z_half[0];
  z[1] += dt*z_half[1];
  }
void Beuler2x2(double dt, double *z, double * F){
  double T[4]= {0}; //* T = (double*)malloc(N_x*N_y*sizeof(double));
  double Ti[4]= {0};// = (double*)malloc(N_x*N_y*sizeof(double));
  double y = 0.0;
  int i;
  int j;
  for( i=0; i<2; ++i){
    for( j=0; j<2; ++j){
      if(i==j){
	y = 1.0;
      }
      else{
	y = 0.0;
      }
      T[i + 2*j] = y - dt*F[i + 2*j]; 
    }
  }
  invert2x2(T, Ti);
  double z_n[2] = {z[0], z[1]};
  mult2x2(Ti, z_n);
  z[0] = z_n[0];
  z[1] = z_n[1];
  //free(T);
  //free(Ti);
}
void Crank(double dt, double *z, double * F){
  double T[4] = {0};// = (double*)malloc(N_x*N_y*sizeof(double));
  double Ti[4]= {0};// = (double*)malloc(N_x*N_y*sizeof(double));
  double  z_n[2] = {z[0], z[1]};
  double y;
  int i, j;
  for( i=0; i<2; ++i){
    for( j=0; j<2; ++j){
      if(i==j){
	y = 1.0;
      }
      else{
	y = 0.0;
      }
      T[i + 2*j] = y + 0.5*dt*F[i + 2*j]; 
    }
  }
  mult2x2(T, z_n);
  for( i=0; i<2; ++i){
    for( j=0; j<2; ++j){
      if(i==j){
	y = 1.0;
      }
      else{
	y = 0.0;
      }
      T[i + 2*j] = y - 0.5*dt*F[i + 2*j]; 
    }
  }
  invert2x2(T, Ti);
  mult2x2(Ti, z_n);
  z[0] = z_n[0];
  z[1] = z_n[1];
  //free(T);
  //free(Ti);
}
void advance(double dt, double * z, int method){
  double F[4] = {0};// = (double*)malloc(N_x*N_y*sizeof(double));//F[i + j*N_y]
  osc(F);
  //printf("%e   %e\n%e    %e\n\n", F[0], F[2], F[1], F[3]);
  if( method == 1){
    Feuler2x2(dt, z, F);  
  }
  else if(method ==2){
    Mid2x2(dt, z, F);
    }
  else if(method == 3){
    Beuler2x2(dt, z, F);
  }
  else if(method ==4){
    Crank(dt, z, F);
  }
  //free(F);
}

double error( z0, t, h ){
  double diff = pow((z0-exact(t)),2);
  return(h*diff);
}
////////////////Main///////////////////
int main(int argc, char **argv){
  // /a.out T N method

  FILE* fid,* fa;
  fid = fopen("oscdata.dat", "w");
  fa = fopen("L21.dat","a");
  double T = atof(argv[1]);
  int N = atoi(argv[2]);
  int method = atoi(argv[3]);
  double dt = T/(double)N;
  double t = 0.0;
  double z[2] = {1., 0.}; 
  double sum = 0.;
  
  int i;
  for(i = 0; i<N; ++i){
    //printf("%e %e %e\n", t, z[0], exact(t));
    fprintf(fid, "%e %e %e\n", t, z[0], exact(t));
    sum += dt*pow(fabs(z[0]-exact(t)),2);
    //printf("%e\n", dt*pow(z[0]-exact(t),2));
    //printf("%e %e %e\n", t, z[0], z[1]);
    advance(dt, z, method);
    t += dt;
  }
  double L2 = sqrt(sum);
  fprintf(fa, "%d %e\n", N, L2);
}
