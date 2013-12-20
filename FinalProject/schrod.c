#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <math.h>

// I want to use all matrices of from A[Nx*i + j] linear arrays
// for a Nx by Ny matrices
void print_array(complex * a, int N){
  int i;for(i=0;i<N;++i){
    printf("%e  ", cabs(a[i]));
  }
  printf("\n\n");
}
void PHI(complex * A, complex * phi, int N){
  phi[N] = 1;
  phi[N-1] = A[(N-1)*N + N-1];
  int i;
  for(i = N-2; i>=0; --i){
    phi[i] = A[(i-1)*N + i-1]*phi[i+1] - A[(i+1 - 1)*N + i-1]*A[(i-1)*N + (i+1)-1]*phi[i+2];
  }
}
void THETA(complex * A, complex * theta, int N){
  theta[0] = 1.; theta[1] = A[(1-1)*N + 1-1];
  int i;
  for(i = 2; i<=N; ++i){
    theta[i] = A[(i-1)*N + i-1] * theta[i - 1] - A[(i-1 - 1)*N + i-1]*A[(i-1)*N + (i-1)-1]*theta[i-2];
  }
}
void inverse_triD(complex * A, complex * Ai, int N){
  complex theta[N+1], phi[N+1];
  THETA(A, theta, N);
  PHI(A, phi, N);
  int i, j;
  for(j=0;j<N;++j){
    for(i=0;i<N;++i){
      if(i <= j){
	complex min = cpow(-1., i+j);
	Ai[i*N + j] = min*theta[i-1 + 1]*phi[j+1]/theta[N];
	//printf("%e\n", Ai[i*N+j]);
	int k;for(k = i;k <= j-1; ++k){
	  Ai[i*N + j] *= A[k*N + k+1]; 
	}
      }
      else if( i > j){
	complex minone = cpow(-1., i+j);
	Ai[i*N + j] = minone*theta[j-1 + 1]*phi[i+1]/theta[N];
	int l;for(l = j; l<= i-1; ++l){
	  Ai[i*N + j] *= A[(l+1)*N + l];
	}
      }
    }
  }
}
void identity(complex * A, int N){
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
void maketriD(complex * A, int N){
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
void printmatrix(complex * A, int N){
  int i, j;
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      printf("%e  ", creal(A[i*N + j]));
    }
    printf("\n");
  }
}
void make_zero(complex * A, int N){
  int i, j;
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      A[i*N + j] = 0.;
    }
  }
}
void multsq_mat(complex * A, complex * B, complex * AB, int N){
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
void matrix_dot(complex * A, complex * u, complex * p, int N){
  int i,j; for(i=0;i<N;++i){
    p[i] = 0.;
    for(j=0;j<N;++j){
      p[i] += A[i*N + j]*u[j];
    }
  }
}
complex D2( complex uipo, complex ui, complex uimo, double dx){
  return((uipo - 2.*ui + uimo)/cpow((complex)dx, 2));
}
int condition(int i, int Nh, int range){
  return(i> Nh-range && i<Nh+range);
}
void initialize_system1(complex * U, int N, double dx){
  int Nh = N/2;
  int i, j;
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      complex x = (i-Nh)*dx;
      complex y = (j-Nh)*dx;
      complex r2 = cpow(x,2)+cpow(y,2);
      U[i*N + j] = cexp(-(r2));
      U[i*N +j] *= sin(3.14*y*2./6.);
    }
  }
}
void initialize_system2(complex * U, int N, double dx){
  int Nh = N/2;
  int i, j;
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      complex x = (i)*dx;
      complex y = (j)*dx;
      complex r2 = cpow(x,2)+cpow(y,2);
      U[i*N + j] = cpow(x,2)*cpow(y,2);
    }
  }
}
void write_U(FILE * fid, complex * U, int N, double ds){
  int i,j;
  for(i=1;i<N-1;++i){
    for(j=1;j<N-1;++j){
      fprintf(fid, "%e   %e   %e\n",(double)i *ds, (double)j *ds, creal(conj(U[i*N+j])*U[i*N+j]));
    }
    //fprintf(fid, "\n");
  }
}
void getrow(complex * U, complex * u, int i, int N){ //note that this does not get the boundboury terms
  int j;for(j=0;j<N;++j){                          
    u[j] = U[i*N + j];
  }
}
void getcol(complex * U, complex * u, int j, int N){ //note that this does not get the boundary terms
  int i;for(i=0;i<N;++i){                          
    u[i] = U[i*N + j]; 
  }
}
void setcol(complex * U, complex *u, int j, int N){
  int i;for(i=0;i<N;++i){
    U[(i)*N + j] = u[i];
  }
}
void setrow(complex * U, complex * u, int i, int N){
  int j;for(j=0;j<N;++j){
    U[(i)*N + j] = u[j];
  }
}
void make_D(complex * D, int N, double dx, double dt, complex pm){
  complex coeff = I*pm*dt/(2.*cpow((complex)dx,2));
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
void add_V(complex * ui, complex * vi, complex * uio, int N, double dt){
  complex coeff = dt*I/2.;
  int i;for(i=0;i<N;++i){
    ui[i] -= coeff*vi[i]*uio[i];
  }
}
void advance_yexp(complex * U, int N, double dt, double dx, complex * V){
  int k; 
  for(k=0;k<N;++k){ //loop over columns of U
    complex Dy[(N)*(N)]; make_D(Dy, N, dx, dt, 1.);  //make update matrix
    complex uy[N]; getcol(U, uy, k, N);
    complex vy[N]; getcol(V, vy, k, N);  
    complex uy_new[N]; matrix_dot(Dy, uy, uy_new, N);//update column
    add_V(uy_new, vy, uy, N, dt);
    setcol(U, uy_new, k, N);
  }
}
void advance_xexp(complex * U, int N, double dt, double dx, complex * V){
  int k;
  for(k=0;k<N;++k){
    complex Dx[(N)*(N)]; make_D(Dx, N, dx, dt, 1.);//make derivative matrix
    complex ux[N], ux_new[N]; getrow(U, ux, k, N);   //get row to update
    complex vx[N]; getrow(V, vx, k, N);
    matrix_dot(Dx, ux, ux_new, N); add_V(ux_new, vx, ux, N, dt);
    setrow(U, ux_new, k, N);
  }
}
void advance_ximp(complex * U, int N, double dt, double dx){
  int k;
  for(k=0;k<N;++k){
    complex Dx[(N)*(N)]; make_D(Dx, N, dx, dt, -1.);//make derivative matrix
    complex Dxi[(N)*(N)]; inverse_triD(Dx, Dxi, N);//invert it
    complex ux[N], ux_new[N]; getrow(U, ux, k, N);   //get row to update
    matrix_dot(Dxi, ux, ux_new, N); setrow(U, ux_new, k, N);
    //    printf("\n"); print_array(ux, N-2);
    //printf("\n"); print_array(ux_new, N-2);
  }
}
void advance_yimp(complex * U, int N, double dt, double dx){
  int k; 
  for(k=0;k<N;++k){ //loop over columns of U
    complex Dy[(N)*(N)]; make_D(Dy, N, dx, dt, -1.);  //make update matrix
    complex Dyi[(N)*(N)]; inverse_triD(Dy, Dyi, N);
    complex uy[N]; getcol(U, uy, k, N);  
    complex uy_new[N]; matrix_dot(Dyi, uy, uy_new, N);//update column
    setcol(U, uy_new, k, N);
  }
}
void gaussbound(complex * U, int N, double dx){
  int Nh = N/2;
  int i, j;for(j=0;j<N;++j){
    i=0;
    complex x = (i)*dx;
    complex y = (j-Nh)*dx;
    complex r2 = cpow(x,2)+cpow(y,2);
    U[i*N + j] = cexp(-r2);  
  }
}
void make_V0(complex * V, int N, double dx){
  int i, j;for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      V[i*N+j] = 0.;
    }
  }
  //printmatrix(V, N);
}
void make_V1(complex * V, int N, double dx){
  int i, j;for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      complex x = (i-N/2)*dx;
      complex y = (j-N/2)*dx;
      if( i == 0. || j == 0 || i==N-1 || j == N-1) V[i*N+j] = 0.;//dont touch boundary
      else V[i*N+j] = cpow(x,2) + cpow(y,2);
    }
  }
  //printmatrix(V, N);
}
void boundary_V1(complex * U, int N, double dx, double t){
  int i, j;for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      if( i == 0 || j == 0) U[i*N+j] = 0.;
      /*else if( i == N-1 || j == N-1){
	complex x = i*dx; complex y = j*dx;
	U[i*N+j] = cpow(x,2)*cpow(y,2)*cexp(I*t);
	}*/
    }
  }
}
void boundary0(complex * U, int N){
  int i, j; complex Utop[N], Ubot[N], Uright[N], Uleft[N];
  getrow(U, Ubot, 0, N); getrow(U, Utop, N-1, N);getcol(U, Uleft, 0, N); getcol(U, Uright, N-1, N);
  setrow(U, Utop, 0, N); setrow(U, Ubot, N-1, N);setcol(U, Uright, N-1, N); setcol(U, Uleft, 0, N);
}
void advance_system(complex * U, int N, double dt, double dx, double t){
  t += dt;
  complex V[N*N]; make_V1(V, N, dx);
  advance_yexp(U, N, dt, dx, V);
  //boundary_V1(U, N, dx, dt);
  advance_ximp(U, N, dt, dx);
  //boundary_V1(U, N, dx, dt);
  advance_xexp(U, N, dt, dx, V);
  //boundary_V1(U, N, dx, dt);
  advance_yimp(U, N, dt, dx);
  //boundary_V1(U, N, dx, dt);
  boundary_V1(U, N, dx, t);
  // gaussbound(U, N, dx);
}

/////////////////////////MAIN/////////////////////////////
int main(int argc, char ** argv){
  // ./a.out T dt
  FILE * fid, * fnew, * f2;
  f2 = fopen("schrod25.dat","w");
  fnew = fopen("schrod5.dat","w");
  fid = fopen("initiali.dat", "w");
  int N = 50;
  double T = atof(argv[1]); double dt = atof(argv[2]); double t = 0.;
  double L = 6.; double ds = L/((double) N);
  complex U[N*N];
  initialize_system1(U, N, ds);
  write_U(fid, U, N, ds);
  int thing = 0;
  while(t<T){
    advance_system(U, N, dt, ds, t);
    t+= dt;
    if( t > 0.25 && thing == 0){
      write_U(f2, U, N, ds);
      thing = 1;
    }
    //break;
  }
  write_U(fnew, U, N, ds);
}
