#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// I want to use all matrices of from A[Nx*i + j] linear arrays
// for a Nx by Ny matrices
void print_array(double * a, int N){
  int i;for(i=0;i<N;++i){
    printf("%e  ", a[i]);
  }
  printf("\n\n");
}
void PHI(double * A, double * phi, int N){
  phi[N] = 1;
  phi[N-1] = A[(N-1)*N + N-1];
  int i;
  for(i = N-2; i>=0; --i){
    phi[i] = A[(i-1)*N + i-1]*phi[i+1] - A[(i+1 - 1)*N + i-1]*A[(i-1)*N + (i+1)-1]*phi[i+2];
  }
}
void THETA(double * A, double * theta, int N){
  theta[0] = 1.; theta[1] = A[(1-1)*N + 1-1];
  int i;
  for(i = 2; i<=N; ++i){
    theta[i] = A[(i-1)*N + i-1] * theta[i - 1] - A[(i-1 - 1)*N + i-1]*A[(i-1)*N + (i-1)-1]*theta[i-2];
  }
}
void inverse_triD(double * A, double * Ai, int N){
  double theta[N+1], phi[N+1];
  THETA(A, theta, N);
  PHI(A, phi, N);
  int i, j;
  for(j=0;j<N;++j){
    for(i=0;i<N;++i){
      if(i <= j){
	double min = pow(-1., i+j);
	Ai[i*N + j] = min*theta[i-1 + 1]*phi[j+1]/theta[N];
	//printf("%e\n", Ai[i*N+j]);
	int k;for(k = i;k <= j-1; ++k){
	  Ai[i*N + j] *= A[k*N + k+1]; 
	}
      }
      else if( i > j){
	double minone = pow(-1., i+j);
	Ai[i*N + j] = minone*theta[j-1 + 1]*phi[i+1]/theta[N];
	int l;for(l = j; l<= i-1; ++l){
	  Ai[i*N + j] *= A[(l+1)*N + l];
	}
      }
    }
  }
}
void identity(double * A, int N){
  int i, j;
  for(i=0;i<N;++i){
    for(j = 0;j<N;++j){
      if(i == j){
	A[i*N + j] = 1.;
      }
      else{
	A[i*N + j] = 0.;
      }
    }
  }
}
void maketriD(double * A, int N){
  int i, j;
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      if( i == j){
	A[i*N + j] = 2.;
      }
      else if(i == j+1){
	A[i*N + j] = 1.;
      }
      else if( i+1 == j){
	A[i*N + j] = 3.;
      }
      else{
	A[i*N + j] = 0.;
      }
  }
}
}
void printmatrix(double * A, int N){
  int i, j;
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      printf("%e  ", A[i*N + j]);
    }
    printf("\n");
  }
}
make_zero(double * A, int N){
  int i, j;
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      A[i*N + j] = 0.;
    }
  }
}
void multsq_mat(double * A, double * B, double * AB, int N){
  make_zero(AB, N);
  int i, j, k;
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      for(k=0;k<N;++k){
	AB[i*N + j] += A[i*N + k]*B[k*N + j];
      }
    }
  }
}
void matrix_dot(double * A, double * u, double * p, int N){
  int i,j; for(i=0;i<N;++i){
    p[i] = 0.;
    for(j=0;j<N;++j){
      p[i] += A[i*N + j]*u[j];
    }
  }
}
double D2( double uipo, double ui, double uimo, double dx){
  return((uipo - 2.*ui + uimo)/pow(dx, 2));
}
int condition(int i, int Nh, int range){
  return(i> Nh-range && i<Nh+range);
}
void initialize_system(double * U, int N){
  int Nh = N/2;int range = N/10;
  int i, j;
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      if( condition(i, Nh, range) && condition(j, Nh, range)){
	U[i*N + j] = 10.;
      }
      else{
	U[i*N + j] = 0.;
      }
    }
  }
}
void initialize_laplace(double * U, int N, double L){
  double dx = 2.*3.14156*L/(double)N;
  int i,j; for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      if(i == 0){
	U[i*N + j] = sin(dx*j);
      }
      else{
	U[i*N + j] = 0.;
      }
    }
  }
}
void enforce_boundary_laplace(double * U, int N, double L){
  double dx = 2.*3.14156*L/(double)N;
  int i,j; for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      if(i == 0){
	U[i*N + j] = sin(dx*j);
      }
      else if(i == N-1 || j == 0 || j == N-1 ) U[i*N + j] = 0.;
    }
  }
  //printf("enforced\n");
}
void write_U(FILE * fid, double * U, int N, double ds){
  int i,j;
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      fprintf(fid, "%e   %e   %e\n",(double)i *ds, (double)j *ds, U[i*N+j]);
    }
    //fprintf(fid, "\n");
  }
}
void getrow(double * U, double * u, int i, int N){ //note that this does not get the boundary terms
  int j;for(j=0;j<N;++j){                          
    u[j] = U[i*N + j];
  }
}
void getcol(double * U, double * u, int j, int N){ //note that this does not get the boundary terms
  int i;for(i=0;i<N;++i){                          
    u[i] = U[i*N + j]; 
  }
}
void setcol(double * U, double *u, int j, int N){
  int i;for(i=0;i<N;++i){
    U[(i)*N + j] = u[i];
  }
}
void setrow(double * U, double * u, int i, int N){
  int j;for(j=0;j<N;++j){
    U[i*N + j] = u[j];
  }
}
void copy_matrix(double * U, double * Uc, int N){
  int i,j;for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      Uc[i*N + j] = U[i*N + j];
    }
  }
}
void make_D(double * D, int N, double dx, double dt, double pm){
  double diffusion = 1.;
  double coeff = diffusion*pm*dt/(2.*pow(dx,2));
  int l,m;
  for(l=0;l<N;++l){
    for(m=0;m<N;++m){
      if(l == m) D[l*N + m] = 1. - coeff*2.;
      else if( l+1 == m) D[l*N + m] = coeff;
      else if( l == m+1) D[l*N + m] = coeff;
      else D[l*N + m] = 0.; 
    }
  }
}
void advance_yexp(double * U, int N, double dt, double dx){
  int k; 
  for(k=0;k<N;++k){ //loop over columns of U
    double Dy[(N)*(N)]; make_D(Dy, N, dx, dt, 1.);  //make update matrix
    double uy[N]; getcol(U, uy, k, N);  
    double uy_new[N]; matrix_dot(Dy, uy, uy_new, N);//update column
    setcol(U, uy_new, k, N);
  }
}
void advance_xexp(double * U, int N, double dt, double dx){
  int k;
  for(k=0;k<N;++k){
    double Dx[(N)*(N)]; make_D(Dx, N, dx, dt, 1.);//make derivative matrix
    double ux[N], ux_new[N]; getrow(U, ux, k, N);   //get row to update
    matrix_dot(Dx, ux, ux_new, N); setrow(U, ux_new, k, N);
  }
}
void advance_ximp(double * U, int N, double dt, double dx){
  int k;
  for(k=0;k<N;++k){
    double Dx[(N)*(N)]; make_D(Dx, N, dx, dt, -1.);//make derivative matrix
    double Dxi[(N)*(N)]; inverse_triD(Dx, Dxi, N);//invert it
    double ux[N], ux_new[N]; getrow(U, ux, k, N);   //get row to update
    matrix_dot(Dxi, ux, ux_new, N); setrow(U, ux_new, k, N);
    //    printf("\n"); print_array(ux, N-2);
    //printf("\n"); print_array(ux_new, N-2);
  }
}
void advance_yimp(double * U, int N, double dt, double dx){
  int k; 
  for(k=0;k<N;++k){ //loop over columns of U
    double Dy[(N)*(N)]; make_D(Dy, N, dx, dt, -1.);  //make update matrix
    double Dyi[(N)*(N)]; inverse_triD(Dy, Dyi, N);
    double uy[N]; getcol(U, uy, k, N);  
    double uy_new[N]; matrix_dot(Dyi, uy, uy_new, N);//update column
    setcol(U, uy_new, k, N);
  }
}
double calc_convergence(double * U, double * U_old, int N){
  double con = 0.;
  int i,j;for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      con+= fabs(U[i*N+j]-U_old[i*N+j]);
    }
  }
  return(con);
}
double advance_system(double * U, int N, double dt, double dx){
  double U_old[N*N]; copy_matrix(U, U_old, N);
  advance_yexp(U, N, dt, dx);
  advance_ximp(U, N, dt, dx);
  advance_xexp(U, N, dt, dx);
  advance_yimp(U, N, dt, dx);
  enforce_boundary_laplace(U, N, 1.);
  double con = calc_convergence(U, U_old, N);
  return(con);
}
/////////////////////////MAIN/////////////////////////////
int main(void){
  // ./a.out T dt
  FILE * fid, * fnew;
  fnew = fopen("relax.dat","w");
  fid = fopen("init.dat", "w");
  int N = 100;
  //double T = atof(argv[1]); 
  double dt = 0.00001; 
  double L = 1.; double ds = L/((double) N);
  double U[N*N];
  double TOL = pow(10., -9); double con = 1; int COUNT = 0;
  initialize_laplace(U, N, L);
  //make_zero(U, N); //enforce_boundary_laplace(U, N, 1.);
  write_U(fid, U, N, ds);
  while(con > TOL){
    con = advance_system(U, N, dt, ds);
    if(COUNT > 300){
      printf("broke\n");
      break;
    }
    COUNT += 1;
    //break;
  }
  write_U(fnew, U, N, ds);
}
