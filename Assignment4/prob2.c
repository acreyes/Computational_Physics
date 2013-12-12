#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double gam = 1.4;
double alph = 1.3;

void initnxN(double *** A, int N, int n){
  *A = (double**)malloc(n*sizeof(double*));
  int i, j;
  for (i=0;i<3;++i){
    (*A)[i] = (double*)malloc(N*sizeof(double));
    for (j=0; j<3;++j){
      (*A)[i][j] = 0.;
    }
  }
}
double f(double x){
  return(exp(-pow((x+0.4)*50,2)));
}
double density(double x){
  return(1+alph*f(x));
}
double Pressurei(double rho){
  double P = pow(rho, gam);
  return(P);
}
double Velocityi(double  rho){
  double gama = 2./(gam-1.);
  double v = gama*sqrt(gam) * (pow(rho,1/gama) - 1);
  return(v);
}
void initialize_system(double ** U, double * rho, double * P, double * V, int NT){
  double v = 0.; double dx = 1./(double)NT; double x = -0.5-dx; int i = 0;
  while(x<.5+dx){
    rho[i] = density(x);
    P[i] = Pressurei(rho[i]);
    V[i] = Velocityi(rho[i]);
    x += dx; ++i;
  }
  for(i=0;i<NT;++i){
    U[0][i] = rho[i];
    U[1][i] = rho[i]*V[i];
    U[2][i] = pow(U[1][i],2)/(2*U[0][i]) + P[i]/(gam - 1);
  }
  //printf("last guy=%f\n", U[0][NT-2]);
}
void Reset_BC(double ** U, int NT){
  double Ul[3], Ur[3];
  int i;for(i=0;i<3;++i){
    Ul[i] = U[i][1];
    Ur[i] = U[i][NT-2];
    U[i][0] = Ur[i];
    U[i][NT-1] = Ul[i];
  }
}
void getcol(double **U, int j, double * col){
  int i;for(i=0;i<3;++i){
    col[i] = U[i][j];
  }
  //printf("j %d column %f  %f   %f\n", j,U[0][j],U[1][j],U[2][j]);
}
double Pressure(double* Ui){
  double Et = Ui[2] - pow(Ui[1],2)/(2*Ui[0]);
  double P = Et*(gam-1);
  //printf("%f  %f  %f  %f\n", Ui[0],Ui[1],Ui[2], P);
  return(P);
}
void make_P(double ** U, int NT, double * Press){
  Press[0] = 1.;
  int i;for(i=1;i<NT-1;++i){
    double Ui[3]; getcol(U, i, Ui);
    Press[i] = Pressure(Ui);
    //printf("Ui_%d %f  %f P=%f\n", i, Ui[0], Ui[1], Press[i]);
  }
  Press[NT] = 0.1;
  //printf("%f\n", Press[NT]);
}
double Velocity(double * Ui){
  double v = Ui[1]/Ui[0];
  //printf("%f   %f   %f\n", Ui[1], Ui[0], v);
  return(v);
}
double C_s(double P, double rho){
  double c = sqrt(gam*P/rho);
  //printf("%f  %f  %f\n", P, rho, c);
  return(c);
}
double lambdapm(double v, double c_s, int pm){
  double lam = v + (double)pm *c_s;
  //printf("%f  %f\n",v, c_s);
  return(lam);
}
double Max(double* a, int N){
  double max = a[0];
  int i; for(i=0;i<N;++i){
    if(max<a[i]) max = a[i];
  }
  return(max); 
}
double alphapm(double* Ui, double* Uimo, int pm, double  Pi, double Pimo){
  double vi = Velocity(Ui); double ci = C_s(Pi, Ui[0]);
  double lami = lambdapm(vi, ci, pm);
  double vimo = Velocity(Uimo); double cimo = C_s(Pimo, Uimo[0]);
  double lamimo = lambdapm(vimo, cimo, pm);
  double tobemaxed[3] = {0, (double)pm*lami, (double)pm*lamimo};
  double alpha = Max(tobemaxed, 3);
  //printf("%f  %f\n", lami, lamimo);
  return(alpha);
}
void Flux(double * Ui, double * Fi, double P){
  Fi[0] = Ui[1]; 
  Fi[1] = pow(Ui[1], 2)/Ui[0] + P;
  Fi[2] = (Ui[2] + P)*Velocity(Ui);
}
void Fi_mh(double * Ui, double * Uimo, double * Fimh, double * alphas, double Pi, double Pimo){
  double Fi[3], Fimo[3];//double  Fhll[3];
  Flux(Ui, Fi, Pi); Flux(Uimo, Fimo, Pimo);
  double alphap = alphapm(Ui, Uimo, 1, Pi, Pimo);
  double alpham = alphapm(Ui, Uimo, -1, Pi, Pimo);
  alphas[0] = alphap; alphas[1] = alpham;   
  int i; for(i=0;i<3;++i){
    Fimh[i] = (alphap*Fimo[i] + alpham*Fi[i] - alphap*alpham*(Ui[i]-Uimo[i]))/(alphap+alpham);
  }
}
void make_F(double ** U, double ** F, double **alpha, int NT, double * P){
  double alphas[2];

  int i,j; for(i=1;i<NT;++i){
    double Ui[3], Uimo[3]; getcol(U, i, Ui); getcol(U, i-1, Uimo);
    double Fimh[3]; Fi_mh(Ui, Uimo, Fimh, alphas, P[i], P[i-1]);
    for(j=0;j<2;++j){
      if(i>0 && i<NT-1){
	alpha[j][i-1]=alphas[j];
      }
    }
    for(j=0;j<3;++j){
      if(i>0 && i<NT){
      F[j][i] = Fimh[j];
      }      
    }
  }
  for(j=0;j<3;++j){
    F[j][0] = F[j][NT-1];
    F[j][NT] = F[j][1];
  }
}
double get_step(double maxalpha, double dx){
  double dt = dx/(10*maxalpha);
  //dt = 0.00001;
  return(dt);
}
double get_maxalpha(double ** alpha, int NT){
  double max = fabs(alpha[0][0]); //printf("here is the first one%f\n", max);
  int i,j; for(i=0;i<NT-2;++i){
    for(j=0;j<2;++j){
      if(max<fabs(alpha[j][i])) max = fabs(alpha[j][i]);
      //printf("other alphas %f\n", alpha[j][i]);
    }
  }
  return(max);
}
double advance_system(double ** U, double dx , double t, int NT, double * P){
 
  double ** alpha; initnxN(&alpha, NT-2, 2);
  make_P(U, NT, P);
  double ** F; initnxN(&F, NT+1, 3); make_F(U, F, alpha, NT, P);
  double maxalpha = get_maxalpha(alpha, NT);
  double dt = get_step(maxalpha, dx);
  double Un[3];
  double L;
  int i,j; for(i=0;i<NT-1;++i){
    for(j=0;j<3;++j){
      L = - dt*(F[j][i+1]-F[j][i])/dx;
      Un[j] = U[j][i] + L;
      //printf("j=%d  L = %f\n", j, L);
      U[j][i] = Un[j];
    }
    //printf("i %d column %f  %f   %f\n", i,U[0][i],U[1][i],U[2][i]);
  }
  Reset_BC(U, NT);
  double tn;
  return(t+dt);
}
void write(double * P, FILE * fid, int NT, double dx){
  double x = -0.5;
  int i;for(i=0;i<NT-2;++i){
    fprintf(fid, "%e    %e    %d\n", x, P[i], i);
    //printf("%d %f\n", i,  P[i]);
    x+= dx;
  }
}
void getRow(double ** U, int NT, double * thing, int row){
  int count = 1;
  int i;for(i=1;i<NT-1;++i){
    thing[i-1] = U[row][i];
    //printf("%d %f\n", i-1, U[row][i]);
    count+=1;
  }
  //printf("before i write=%f NT-2=%d\n", thing[NT-3], NT-3);
}
double Entropy(double ** U, int NT, double dx, double * sumv, double * sump){
  double vsum = 0;
  double psum = 0;
  int i;for(i=1;i<NT-1;++i){
    double Ui[3]; getcol(U, i, Ui);
    double P = Pressure(Ui);
    double V = Velocity(Ui);
    double Pi = Pressurei(Ui[0]);
    double Vi = Velocityi(Ui[0]);
    vsum += dx*fabs(V - Vi);
    psum += dx*fabs(P - Pi);
  }
  *sumv = vsum; *sump = psum;
}

int main(int argc, char ** argv){
  // ./a.out Nx
  FILE * fid1, * fid2, * fid3;
  fid1 = fopen("initial2.dat", "w");
  fid3 = fopen("entropy.dat", "w");
  double T = atof(argv[1]); double t = 0; 
  double L = 1.;   int Nx = 5000; double dx = L/(double)Nx;int NT = Nx + 2;
  double **U; double rho[Nx+2]; double P[Nx+2];double V[Nx+2];
  initnxN(&U, Nx+2, 3);
  initialize_system(U, rho, P, V, NT);
  //printf("i write=%f\n", U[0][NT-2]);
  double dense[Nx]; getRow(U, NT, dense, 0);
  write(dense, fid1, NT, dx);


  //advance_system( U , dx , t , NT , P );
  //printf("t=%f\n", t); 
  //sleep(1);
  
  while(t<T){
    double tn = advance_system(U, dx, t, NT, P); 
    double sumv, sump;
    Entropy(U, NT, dx, &sumv, &sump);
    fprintf(fid3, "%e   %e   %e\n", t, sumv, sump);
    //printf("t=%f\n", tn); 
    //sleep(1);
    t = tn;
    //printf("i 21 column %f  %f  %f", U[0][21], U[1][21], U[2][21]);
    //break;

    }
  fid2 = fopen("final2.dat", "w");
  getRow(U, NT, dense, 0);
  write(dense, fid2, NT, dx);
  
  return(0);
}

