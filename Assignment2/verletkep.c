#include <stdio.h>
#include <math.h>
double a2(double x){
  return(pow(x, 2));
}
double h(double M, double E, double e){
  return(M - E + e*sin(E));
}
double dfdE( double E, double e){
  return(-1.0 + e*cos(E));
}
double del( double M, double E, double e){
  return(-h(M, E, e)/dfdE(E, e));
}
double Energy( double x, double y, double v_x, double v_y){
  double E = 0.5*(a2(v_x) + a2(v_y)) - 1/(sqrt( a2(x) + a2(y) ));
  return(E);
}
double v(double v_x, double v_y){
  return(sqrt(pow(v_x, 2) + pow(v_y, 2)));
}
double dv_idt( double t, double i, double j, double v){
  double denom = pow( a2(i) + a2(j), 1.5);
  return( -i/denom );
}
void verlet( double t, double dt, double x, double y, double v_x, double v_y, double* x_n, double* y_n, double* v_nx, double* v_ny ){
  /////solves for x(t) and v(t+dt/2) using x(t-dt) and v(t-dt/2)//////
  *v_nx = v_x + dt*dv_idt( t, x, y, v_x );
  *v_ny = v_y + dt*dv_idt( t, y, x, v_y );

  *x_n = x + dt*(*v_nx);
  *y_n = y + dt*(*v_ny); 
}

double dr2( double dt, double x, double y, double x_t, double y_t){
  double dx = x - x_t;
  double dy = y - y_t;
  return( dt*(pow(dx, 2) + pow(dy, 2)));   
}
int main(void){
  printf("%f\n", dr2(1.0, 0.0, 1, 0.5, 2.0));
  FILE *fid, *fl2;
  fl2 = fopen("l2.dat", "w");
  fid = fopen("kepler.dat", "w");

  double x, y, a, b, E, E_n, e, f, M, l, phi, v_r, v_p, v_x, v_y, U, T, dt, t, tol, sum;
  double xv, yv, v_xv, v_yv, x_n, y_n, v_nx, v_ny; 
  int N = 10;
  tol = 10e-11;
  T = 6.0;
  dt = T/(double)N;
  t = 0.0;
  x = 1.2;
  y = 0.0;
  v_x = -0.0;
  v_y = 1.0;
  xv = x;
  yv = y;
  v_xv = v_x;
  v_yv = v_y;
  double r = sqrt(pow(x,2) + pow(y,2));

  U = Energy(x, y, v_x, v_y);
  l = (v_y*x + v_x*y);

  a = 0.5*1/(fabs(U));
  b = l/(sqrt(2.0*fabs(U)));
  f = sqrt(1.0 - 2.0*pow(l,2)*fabs(U))/(2.0*fabs(U));
  e = f/a;

  int i, j;
  for (j = 0; j< 1000; j++){
    int count;
    
    sum = 0;
    for(i = 0; i<N; i++){
      if (j == 0){
	fprintf(fid, "%e %e %e\n", t, x, y);
      }
      //do newton rapheson//
      count = 0;
      while(count < 30){
	M = t*l/(a*b);
	E_n = E + del(M, E, e);
	E = E_n;
	if ( h( M, E, e) < tol){
	  {break;}
	}
	count += 1;
      }
      //////calculate x and y////////
      x = a*cos(E) - f;
      y = b*sin(E);
      //////////do verlet///////////
      verlet( t, dt, xv, yv, v_xv, v_yv, &x_n, &y_n, &v_nx, &v_ny);
      xv = x_n;
      yv = y_n;
      v_xv = v_nx;
      v_yv = v_ny;

      /////Integrate L2 Error////////
      sum += dr2(dt, xv, yv, x, y);
      t += dt;
    }
    fprintf(fl2, "%d %e\n", N, sqrt(fabs(sum)));
    printf("%d %e\n", N, fabs(sum) ) ;
    N += 100;
    dt = T/(double)N;
    t = 0;
    x = 1.2;
    y = 0.0;
    v_x = -0.0;
    v_y = 1.0;
    xv = x;
    yv = y;
    v_xv = v_x;
    v_yv = v_y;
  }
}
