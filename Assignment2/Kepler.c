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
int main(void){
  FILE *fid;
  fid = fopen("kepler.dat", "w");

  double x, y, a, b, E, E_n, e, f, M, l, phi, v_r, v_p, v_x, v_y, U, T, dt, t, tol;
  int N = 5000;
  tol = 10e-11;
  T = 200.0;
  dt = T/(double)N;
  t = 0.0;
  x = 1.2;
  y = 0;
  v_x = -0.0;
  v_y = 1.0;
  double r = sqrt(pow(x,2) + pow(y,2));

  U = Energy(x, y, v_x, v_y);
  l = (v_y*x + v_x*y);

  a = 0.5*1/(fabs(U));
  b = l/(sqrt(2.0*fabs(U)));
  f = sqrt(1.0 - 2.0*pow(l,2)*fabs(U))/(2.0*fabs(U));
  e = f/a;
  printf("%f %f %f %f %f %f\n", a, b, f, e, l, U);
  int i;
  int count;
  for(i = 0; i<N; i++){
    fprintf(fid, "%e %e %e\n", t, x, y);
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
    t += dt;
  }
}
