#include <stdio.h>
#include <math.h>



double a2(double x){
  double temp = pow(x,2);
  return( temp );
}

double dv_idt( double t, double i, double j, double v){
  double denom = pow( a2(i) + a2(j), 2.5);
  return( -i/denom );
}

double didt( double t, double i, double v_i){
  return(v_i);
}

double rk4( double x, double h1, double h2, double h3, double h4){
  return(x + (h1 + 2*h2 + 2*h3 + h4)/6);
}

    


double U(double x, double y){
  double denom = pow( a2(x) + a2(y), 0.5);
    return(1.0/denom);
}

double Edif(double E0, double Et){
  double dif = Et - E0;
  return(log10(fabs(dif/E0)));
}

double KE(double v_x, double v_y){
  double v2 = a2(v_x) + a2(v_y);
  return(0.5*pow(v2, 0.5));
}
	 

int main(void){

  double T = 100.0;
  double t = 0.0;
  double y = 0.0;
  double x = 1.1;
  double v_x = 0.0;
  double v_y = 1.0;
  int N = 1000;
  double dt = T/(double)N;
  double U0 = U(x, y);
  double Ut = U(x, y);
  double KEt = KE(v_x, v_y);
  double KE0 = KE(v_x, v_y);
  double Udif = Edif(U0, Ut);
  double KEdif = Edif(KE0, KEt);
  double ETdif = Edif(KE0+U0, KEt+Ut);

  FILE *fp, *fe;
  fp = fopen("asign1.dat", "w");
  fe = fopen("Energy.dat", "w");

  int i;
  for(i = 0; i<N; i++){
    fprintf(fp, "%f %f %f\n", t, x, y);
    fprintf(fe, "%f %f %f %f\n", t, Udif, KEdif, ETdif);
    /*
    //euler
    double v_nx = v_x + dv_idt(t, x, y, v_x)*dt;
    double x_n = x + didt(t, x, v_x)*dt;

    double v_ny = v_y + dv_idt(t, y, x, v_y)*dt;
    double y_n = y + didt(t, y, v_y)*dt;
    //end euler
    */
    /*
    //Midpoint method
    double kvx1 = dv_idt(t, x, y, v_x)*dt;
    double kvy1 = dv_idt(t, y, x, v_y)*dt;
    double kx1 = didt(t, x, v_x)*dt;
    double ky1 = didt(t, y, v_y)*dt;

    double kvx2 = dv_idt(t+0.5*dt, x+0.5*kx1, y + 0.5*ky1, v_x+0.5*kvx1)*dt;
    double kvy2 = dv_idt(t+0.5*dt, y+0.5*ky1, x + 0.5*kx1, v_y+0.5*kvy1)*dt;
    double kx2 = didt(t+0.5*dt, x+0.5*kx1, v_x+0.5*kvx1)*dt;
    double ky2 = didt(t+0.5*dt, y+0.5*ky1, v_y+0.5*kvy1)*dt;


    double x_n = x + kx2;
    double y_n = y + ky2;

    double v_nx = v_x + kvx2;
    double v_ny = v_y + kvy2;
    //end midpoint
    */
    //Start RK4
    
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

    double v_nx = rk4(v_x, hvx1, hvx2, hvx3, hvx4);
    double v_ny = rk4(v_y, hvy1, hvy2, hvy3, hvy4);

    double x_n = rk4(x, hx1, hx2, hx3, hx4);
    double y_n = rk4(y, hy1, hy2, hy3, hy4);


    //end RK4
    

    //make full step
    t = t + dt;
    x = x_n;
    y = y_n;
    v_x = v_nx;
    v_y = v_ny;    

    Ut = U(x, y);
    KEt = KE(v_x, v_y);
    Udif = Edif(U0, Ut);
    KEdif = Edif(KE0, KEt);
    ETdif = Edif(KE0+U0, KEt+Ut);
    
  }
  return(0);
}

  





