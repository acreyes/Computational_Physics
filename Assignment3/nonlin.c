#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void advance(double  h, double * z, int method);
double TOL = 1.e-3;
double h0 = 1.e-5;
double K[3] = {1. , 100. , 1000.};
////////////PHYSICS///////////////
void F(double * z, double * z_n){
  double An = -K[0]*z[0]*z[1];
  double Xn = K[0]*z[0]*z[1] - K[1]*z[1]*z[2];
  double Yn = K[1]*z[1]*z[2] - K[2]*z[2];
  z_n[0] = An;
  z_n[1] = Xn;
  z_n[2] = Yn;
}
void DF(double * z, double ** dF){
  dF[0][0] = -K[0]*z[1];
  dF[0][1] = -K[0]*z[0];
  dF[0][2] = 0.0;
  dF[1][0] = K[0]*z[0];
  dF[1][1] = K[0]*z[0] - K[1]*z[2];
  dF[1][2] = -K[1]*z[1];
  dF[2][0] = 0.0;
  dF[2][1] = K[1]*z[2];
  dF[2][2] = K[1]*z[1] - K[2];
  //this is a -dF thats why theres a plus later
}
////////////MATH/////////////////
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
void mult33(double ** A, double ** B, double ** C){
  int i, j, k;
  for(i=0;i<3;++i){
    for(j=0;j<3;++j){
      C[i][j] = 0.;
      for(k=0;k<3;++k){
	C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
}
void dot33(double ** A, double * z, double * z_n){
  int i,j;
  double zn[3] = {0,0,0};
  for(i=0;i<3;++i){
    for(j=0;j<3;++j){
      zn[i] += A[i][j]*z[j];
    }
  }
  z_n[0] = zn[0];
  z_n[1] = zn[1];
  z_n[2] = zn[2];
}
void init33(double *** A){
  *A = (double**)malloc(3*sizeof(double*));
  int i, j;
  for (i=0;i<3;++i){
    (*A)[i] = (double*)malloc(3*sizeof(double));
    for (j=0; j<3;++j){
      (*A)[i][j] = 0.;
    }
  }
}
double getstep(double * z, int method){
  double z1[3];
  double z2[3];
  int i;
  for(i=0;i<3;++i){
    z1[i] = z[i];
    z2[i] = z[i];
  }
  double h1 = h0;
  double h2 = 2.*h0;
  double order = 1.0;

  advance( h1, z1, method);
  advance( h1, z1, method);
  advance( h2, z2, method);
 
  double delta = 0.;
  double del;
  for(i=0;i<3;++i){
    del = fabs(z1[i] - z2[i]);
    if(del>delta){
      delta = del;
    }
  }
  double h_n = pow(TOL/delta, 1./order)*h0;
  if(h_n > 0.1){
    h_n = 0.1;
  }
  else if(h_n < 4.e-5){
    h_n = 4.e-5;
  }
  return(h_n);
}
///////////INTEGRATION FUNCTIONS////////////
void Feuler(double * z, double h){
  double zn[3];
  F(z, zn);
  z[0] += h*zn[0];
  z[1] += h*zn[1];
  z[2] += h*zn[2];
}
void Beuler(double * z, double h){
  double ** dF, **Ti, **T;
  init33(&dF);
  init33(&Ti);
  init33(&T);
  DF(z, dF);
  double y = 0.0;
  int i, j;
  for( i=0; i<3; ++i){
    for( j=0; j<3; ++j){
      if(i==j){
	y = 1.0;
      }
      else{
	y = 0.0;
      }
      T[i][j] = y - h*dF[i][j]; 
    }
  }
  inverse33(T, Ti);
  double zo[3], Fz[3];
  F(z, Fz);
  dot33(Ti, Fz, zo);
  for ( i = 0;i<3;++i){
    z[i]+= h*zo[i];
  }
}
void test(double * z, double  h){
  int i;
  for (i=0;i<3;++i){
    z[i] = z[i] - 2*h*z[i];
  }
}

void advance(double  h, double * z, int method){
  if (method == 1){
    Feuler(z, h);
  }
  else if (method == 2){
    Beuler(z, h);
  }
  else if (method == 7){
    test(z, h);
  }
}
//////////////MAIN//////////////
int main(int argc, char **argv){
  // ./a.out method h
  double T = 3.0;
  double A = 150.;
  double X = 10.;
  double Y = 10.;
  double z[3] = {A, X, Y};
  double t = 0.;
  FILE * fid;
  double h;
  //h = atof(argv[2]);
  int method = atoi(argv[1]);
  fid = fopen("best.dat", "w");
  while(t<T){
    fprintf(fid, "%e %e %e %e\n", t, z[0], z[1], z[2]);
    h = getstep(z, method);
    //printf("%e\n", h);
    advance(h, z, method);
    printf("%e\n", t);
    t+= h;
    }
}
