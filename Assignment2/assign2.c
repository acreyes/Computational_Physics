#include <stdio.h>
#include <math.h>
double a2(double x){
  double temp = pow(x,2);
  return( temp );
}
//integration functions
double dv_idt( double t, double i, double j, double v){
  double denom = pow( a2(i) + a2(j), 1.5);
  return( -i/denom );
}
double didt( double t, double i, double v_i){
  return(v_i);
}
//end integration functions
//solutions
double rk4final( double x, double h1, double h2, double h3, double h4){
  return(x + (h1 + 2*h2 + 2*h3 + h4)/6);
}
void rk4( double t, double dt, double x, double y, double v_x, double v_y, double* x_n, double* y_n, double* v_nx, double* v_ny){
  double hvx1 = dv_idt(t, x, y, v_x)*dt;
  double hvy1 = dv_idt(t, y, x, v_y)*dt;
  double hx1 = didt(t, x, v_x)*dt;
  double hy1 = didt(t, y, v_y)*dt;

  double hvx2 = dv_idt(t + 0.5*dt, x + 0.5*hx1, y + 0.5*hy1, v_x + 0.5*hvx1)*dt;
  double hvy2 = dv_idt(t + 0.5*dt, y + 0.5*hy1, x + 0.5*hx1, v_y+0.5*hvy1)*dt;
  double hx2 = didt(t+0.5*dt, x+0.5*hx1, v_x+0.5*hvx1)*dt;
  double hy2 = didt(t+0.5*dt, y+0.5*hy1, v_y+0.5*hvy1)*dt;

  double hvx3 = dv_idt(t + 0.5*dt, x + 0.5*hx2, y + 0.5*hy2, v_x + 0.5*hvx2)*dt;
  double hvy3 = dv_idt(t + 0.5*dt, y + 0.5*hy2, x + 0.5*hx2, v_y + 0.5*hvy2)*dt;
  double hx3 = didt(t+0.5*dt, x + 0.5*hx2, v_x + 0.5*hvx2)*dt;
  double hy3 = didt(t + 0.5*dt, y + 0.5*hy2, v_y + 0.5*hvy2)*dt;

  double hvx4 = dv_idt(t + dt, x + hx3, y + hy3, v_x + hvx3)*dt;
  double hvy4 = dv_idt(t + dt, y + hy3, x + hx3, v_y + hvy3)*dt;
  double hx4 = didt(t + dt, x + hx3, v_x + hvx3)*dt;
  double hy4 = didt(t+dt, y + hy3, v_y + hvy3)*dt;

  *v_nx = rk4final(v_x, hvx1, hvx2, hvx3, hvx4);
  *v_ny = rk4final(v_y, hvy1, hvy2, hvy3, hvy4);
  *x_n = rk4final(x, hx1, hx2, hx3, hx4);
  *y_n = rk4final(y, hy1, hy2, hy3, hy4); 
}

void verlet( double t, double dt, double x, double y, double v_x, double v_y, double* x_n, double* y_n, double* v_nx, double* v_ny ){
  /////solves for x(t) and v(t+dt/2) using x(t-dt) and v(t-dt/2)//////
  *v_nx = v_x + dt*dv_idt( t, x, y, v_x );
  *v_ny = v_y + dt*dv_idt( t, y, x, v_y );

  *x_n = x + dt*(*v_nx);
  *y_n = y + dt*(*v_ny); 
}

double Energy( double x, double y, double v_x, double v_y){
  double E = 0.5*(a2(v_x) + a2(v_y)) - 1/(sqrt( a2(x) + a2(y) ));
  return(E);
}

///////////////Begin Code///////////////////

int main(void){
  /////////////Set Initial Conditions///////
  double T = 200.0;
  double t, xv, yv, v_xv, v_yv, x4, y4, v_x4, v_y4, v_nx, v_ny, x_n, y_n, x_old, y_old, v_yold, v_xold, E4, Ev, v_x, v_y, E0;
  t = 0;
  int N = 5000; /////////////////////
  int i;
  double dt = T/(double)N;
  FILE *fid;
  fid = fopen("assign2.dat", "w");
  
  x4 = 1.2;
  y4 = 0.0;
  xv = x4;
  yv = y4;
  v_x4 = -0.0;
  v_y4 = 1.0;
  //////////RK4 to get next value/////////////
  rk4(t, dt, x4, y4, v_x4, v_y4, &x_n, &y_n, &v_nx, &v_ny);
  //////////Average v's to get v at the half step//////
  v_xv = 0.5*(v_x4 + v_nx);
  v_yv = 0.5*(v_y4 + v_ny); 
  E0 = Energy(x4, y4, v_x4, v_y4);
  E4 = fabs( (Energy(x4, y4, v_x4, v_y4) - E0)/E0);
  Ev = E4; 
  printf("%f %f\n", E4, Ev);

  ///////loop over everything//////////
  for(i=0; i<N; ++i){
    fprintf(fid, "%e %e %e %e %e %e %e\n", t, x4, y4, xv, yv, E4, Ev );
    /////////Run RK4/////////////
    /*
    rk4(t, dt, x4, y4, v_x4, v_y4, &x_n, &y_n, &v_nx, &v_ny);
    x4 = x_n;
    y4 = y_n;
    v_x4 = v_nx;
    v_y4 = v_ny;
    E4 = fabs( (Energy(x4, y4, v_x4, v_y4) - E0)/E0);
    */
    ////////Do Verlet///////////
    verlet( t, dt, xv, yv, v_xv, v_yv, &x_n, &y_n, &v_nx, &v_ny);
    xv = x_n;
    yv = y_n;
    v_x = 0.5*(v_xv + v_nx);
    v_y = 0.5*(v_yv + v_ny);
    v_xv = v_nx;
    v_yv = v_ny;
    Ev = fabs( (Energy(xv, yv, v_xv, v_yv) - E0)/E0 );
    t += dt;
  }

}
