#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double gam = 1.4;

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
void initialize_system(double ** U, double * rho, double * P, double * V, int N){
  double v = 0.; double dx = 1./(double)N; double x = -0.5-dx; int i = 0;
  int NT = N + 2;
  while(x<0){
    rho[i] = 1.;
    P[i] = 1.;
    x += dx; ++i;
  }
  while(x<.5+dx){
    //printf("%f\n", x);
    rho[i] = 0.125;
    P[i] = 0.1;
    x+= dx; ++i;
  }
  int count = 0;
  for(i=0;i<NT;++i){
    V[i] = 0;
    U[0][i] = rho[i];
    U[1][i] = rho[i]*v;
    U[2][i] = pow(U[1][i],2)/(2*U[0][i]) + P[i]/(gam - 1);
    //printf("i %d start %f  %f   %f\n", i,U[0][i],U[1][i],U[2][i]);
  }
}
void Reset_BC(double ** U, int NT){
  U[0][0] = 1.;
  U[1][0] = 0.;
  U[2][0] = 2.5;

  U[0][NT-1] = .125;
  U[1][NT-1] = 0.;
  U[2][NT-1] = 0.25;
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
  Fi[1] = pow(Ui[1], 2)/Ui[0] + Pressure(Ui);
  Fi[2] = (Ui[2] + Pressure(Ui))*Velocity(Ui);
}
void Fi_mh(double * Ui, double * Uimo, double * Fimh, double * alphas, double Pi, double Pimo){
  double Fi[3], Fimo[3];//double  Fhll[3];
  Flux(Ui, Fi, Pi); Flux(Uimo, Fimo, Pimo);
  double alphap = alphapm(Ui, Uimo, 1, Pi, Pimo);
  double alpham = alphapm(Ui, Uimo, -1, Pi, Pimo);
  //printf("alpha=%f  %f  Pi = %f\n", alphap, alpham, Pi);
  alphas[0] = alphap; alphas[1] = alpham; 
  int i; for(i=0;i<3;++i){
    Fimh[i] = (alphap*Fimo[i] + alpham*Fi[i] - alphap*alpham*(Ui[i]-Uimo[i]))/(alphap+alpham);
  }
  //printf("Fimh = %f\n", Fimh[0], Fimh[1], Fimh[2]); 
}
void make_F(double ** U, double ** F, double **alpha, int NT, double * P){
  double alphas[2];
  F[0][0] = 0.; F[1][0] = 1.; F[2][0] = 0.;
  int i,j; for(i=1;i<NT;++i){
    double Ui[3], Uimo[3]; getcol(U, i, Ui); getcol(U, i-1, Uimo);
    double Fimh[3]; Fi_mh(Ui, Uimo, Fimh, alphas, P[i], P[i-1]);
    for(j=0;j<2;++j){
      if(i>0 && i<NT){
	alpha[j][i-1]=alphas[j];
      }
      //printf("alphas %f   %f\n", alpha[j][i], alphas[j]); 
    }
    for(j=0;j<3;++j){
      if(i>0 && i<NT){
      F[j][i] = Fimh[j];
      }      
    }
    // printf("Fi_%d = %f  %f  %f\n",i, F[0][i],F[1][i],F[2][i]);
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
}
double advance_system(double ** U, double dx , double t, int NT, double * P){
  double ** alpha; initnxN(&alpha, NT-2, 2);
  make_P(U, NT, P);
  double ** F; initnxN(&F, NT, 3); make_F(U, F, alpha, NT, P);
  double maxalpha = get_maxalpha(alpha, NT);
  double dt = get_step(maxalpha, dx);
  double Un[3];
  //printf(" column %f  %f   %f  F[1][1] = %f\n", U[0][0],U[1][0],U[2][0], F[1][1]);
  double L;
  int i,j; for(i=0;i<NT;++i){
    for(j=0;j<3;++j){
      L = - dt*(F[j][i+1]-F[j][i])/dx;
      Un[j] = U[j][i] + L;
      //printf("j=%d  L = %f\n", j, L);
      U[j][i] = Un[j];
    }
    //printf("i %d column %f  %f   %f\n", i,U[0][i],U[1][i],U[2][i]);
  }
  //printf("%f\n", dt);
  Reset_BC(U, NT);
  for(i=0;i<NT;++i){
    //printf("i %d column %f  %f   %f\n", i,U[0][i],U[1][i],U[2][i]);
  }
  t += dt;
  return(t);
}
void write(double ** U, FILE * fid, int NT, double dx){
  double x = -0.5;
  int i;for(i=1;i<NT-1;++i){
    fprintf(fid, "%e    %e   %e   %e\n", x, U[0][i], U[1][i], U[2][i]);
    x+= dx;
  }
}
void getRow(double ** U, int NT, double * thing, int row){
  int count = 1;
  int i;for(i=1;i<NT-1;++i){
    thing[i-1] = U[row][i];
    count+=1;
  }
}

int main(int argc, char ** argv){
  // ./a.out Nx
  FILE * fid1, * fid2;
  fid1 = fopen("initial.dat", "w");
  double T = 0.25; double t = 0; double tn; 
  double L = 1.;   int Nx = atoi(argv[1]); double dx = L/(double)Nx;int NT = Nx + 2;
  double **U; double rho[Nx+2]; double P[Nx+2];double V[Nx+2];
  initnxN(&U, Nx+2, 3);
  initialize_system(U, rho, P, V, Nx);
  write(U, fid1, NT, dx);
  while(t<T){
    tn = advance_system(U, dx, t, NT, P); t = tn; //printf("t=%f\n", t);
    //printf("i 21 column %f  %f  %f", U[0][21], U[1][21], U[2][21]);
    //break;
    //printf("E = %f\n", U[2][1]);
     }
  getRow(U, NT, rho, 0);
  fid2 = fopen("final.dat", "w");
  write(U, fid2, NT, dx);
}
