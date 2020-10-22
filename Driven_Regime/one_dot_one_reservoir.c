#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  double k;
  double T;
  double mu;
  double up;
  double down;
  double fermi;
} RES;


void set_to_zero(int *, double *, double *, int );
double energy(double,  double, double, double , char *, int);
double p_theory(double , double, double, double, char *);
double rate(double,double , RES, char * );
double power(double, int);
double absolute(double);
void input_dealer(char **, int , int *, double *, int *, double *, double *, double *, double *);

int main(int argc, char *argv[])
{
  // ITERATION VARIABLES
  int n_dt, MAXSTRING = 4, MAXSTEPS, N;
  double time[MAXSTEPS], t;
  double MAXTIME;
  double dt;
  double T,mu,k;

  char *type;
  type = (char* )malloc(5*sizeof(char));
  if (argc > 1 )
  {
    input_dealer(argv, argc, &N, &MAXTIME, &MAXSTEPS, &dt, &T, &mu, &k);
    type = argv[8];
  }
  else
  {
    printf("%s\n","HERE");
    N = 1000;
    MAXTIME = 100.0;
    MAXSTEPS = 1000;
    dt = MAXTIME / MAXSTEPS;
    T = 1.0;
    mu = 1.0;
    k = 1.0;
  }
  for (int i=0; i<MAXSTEPS ; ++i)
  {
    time[i] = i*dt;
  }

  // RESERVOIR VALUES
  RES res;
  res.T = T;
  res.k = k;
  res.mu = mu;

  // ENERGY Boundary conditions
  double e0, e1;
  e0 = -10*T ;
  e1 = 10*T;

  // ARRRAYS
  int *x[N];
  double *dq[N], *dw[N];
  double *ddQ[N], *ddW[N];
  double *p, *dQ, *dW, *en, *sigma, *pTheory;
  p = (double *)malloc(MAXSTEPS* sizeof(double));
  dW = (double *)malloc(MAXSTEPS* sizeof(double));
  dQ = (double *)malloc(MAXSTEPS* sizeof(double));
  en = (double *)malloc(MAXSTEPS* sizeof(double));
  pTheory = (double *)malloc(MAXSTEPS* sizeof(double));
  for (int i=0; i<N; i++)
  {
    x[i] = (int *)malloc(MAXSTEPS * sizeof(int));
    dq[i] = (double *)malloc(MAXSTEPS* sizeof(double));
    dw[i] = (double *)malloc(MAXSTEPS* sizeof(double));
    ddQ[i] = (double *)malloc(MAXSTEPS* sizeof(double));
    ddW[i] = (double *)malloc(MAXSTEPS* sizeof(double));
  }

  for (int i=0; i<N ; ++i)
  {
    ///// Main Algorythm

    double t = 0, old_e = 0, r,e;
    int n;
    // Set arrays to zero
    set_to_zero(x[i], dq[i], dw[i], MAXSTEPS );
    if (strcmp(type, "up") == 0)
    {
      n = 1;
    }
    else
    {
      n = 0;
    }
    for (int j=0; j<MAXSTEPS; j++)
    {
      if (j != 0)
      {
        old_e = energy(time[j-1], e0,e1,MAXTIME, type, MAXSTEPS );
      }
      else
      {
        old_e = energy(time[j], e0, e1, MAXTIME, type, MAXSTEPS);
      }
      e = energy(time[j], e0,e1,MAXTIME, type, MAXSTEPS);
      t = time[j];
      res.up = rate(t, e,res, "up");
      res.down = rate(t, e,res, "down");
      r = (double) rand() / ((double) RAND_MAX);
      if (n == 0)
      {
        if (r < res.up * dt)
        {
          dw[i][j] = 0. ;
          dq[i][j] = e;
          n = 1;
        }
        else
        {
          dq[i][j] = 0.;
          dw[i][j] = 0.;
        }
      }
      else
      {
        if (r < res.down * dt)
        {
          dq[i][j] = -e;
          dw[i][j] = 0.;
          n = 0;
        }
        else
        {
          dw[i][j] = e - old_e;
          dq[i][j] = 0.;
        }
      }
      x[i][j] = n;
    }

  }

  //Time deribvative of the work and heat
  for (int i=0 ; i< N; i++)
  {
    for (int j=0; j<MAXSTEPS; j++)
    {
      ddQ[i][j] = (dq[i][j]) / dt;
      ddW[i][j] = (dw[i][j]) / dt;
    }
  }
  // Mean values
  double sumx=0, sumdq=0, sumdw=0;
  for (int i = 0; i<MAXSTEPS ; i++)
  {
    sumx = 0;
    sumdq = 0;
    sumdw = 0;
    for (int j=0; j<N; j++)
    {
      sumx += x[j][i];
      sumdq += ddQ[j][i];
      sumdw += ddW[j][i];
    }
    p[i] = sumx / N;
    dQ[i] = sumdq / N;
    dW[i] = sumdw / N;
  }

  // Write data to File
  FILE *f;
  f = fopen("data.txt","w");
  fprintf(f, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","time", "p","x", "dw","dw_0", "dQ" ,"dQ_0", "sigma","e", "p_theory");
  for (int i=0; i < MAXSTEPS; i++)
  {
    en[i] = energy(time[i], e0,e1,MAXTIME, type, MAXSTEPS);
    pTheory[i] = p_theory(time[i], e0,e1,MAXTIME, type);
    fprintf(f,"%f,%f,%d,%f,%f,%f,%f,%f,%f,%f\n", time[i], p[i], x[0][i], dW[i],dw[0][i], dQ[i],dq[0][i], sigma[i],en[i], pTheory[i]);
  }
  fclose(f);

  return 0;
}

double rate(double t,double e, RES res, char *flag)
{
  /*
    Calculates the rates
  */
  double d;
  if (strcmp(flag, "up") == 0)
  {
    d = exp( e  / res.T)  + 1;
  }
  else
  {
    d = exp( - e  / res.T)  + 1;
  }
  return res.k / d;
}

void set_to_zero(int x[], double dq[] , double dw[] , int MAXSTEPS )
{
  for (int i=0; i<MAXSTEPS; i++)
  {
    x[i] = 0;
    dq[i] = 0.;
    dw[i] = 0.;
  }
}


void input_dealer(char *a[], int len_a, int *N, double *MAXTIME, int *MAXSTEPS, double *dt, double *T, double *mu, double *k)
{
  /*
  Takes the command line inputs and updates the corresponding variables
  */

  char *b, *new_flag = "false";
  if (len_a == 9)
  {
    *N = strtof(a[1],&b) ;
    *MAXTIME = strtof(a[2],&b);
    *MAXSTEPS = strtof(a[3],&b);
    *dt = strtof(a[4],&b);
    *T = strtof(a[5],&b);
    *mu = strtof(a[6],&b);
    *k = strtof(a[7],&b);
  }
}

double energy(double t,  double e0, double e1, double MAXTIME, char *type,int MAXSTEPS)
{
  /*
    Calculate the energy
  */
  double K, k;
  if (strcmp(type, "up") == 0)
  {
    k = (e1 - e0) / (MAXTIME + 2);
    K = (1. / 4) * k*k;
    return e0 + 2*(t + 1) * sqrt(K);
  }
  else if (strcmp(type, "down") == 0)
  {
    double tmp;
    tmp =  e1;
    e1 = e0;
    e0 = tmp;
    k = (e1 - e0) / (MAXTIME + 2);
    K = (1. / 4) * k*k;
    return e0 - 2*(t + 1) * sqrt(K);
  }
  else
  {
    double gamma, omega, T, value, test;
    e0 = 0.1;
    T = 10;
    gamma = 10;
    omega = 2*M_PI / T;
    test = sin(omega * t);
    value = test;
    if (test < 0)
    {
      value = - test;
    }
    return e0 * (1 + gamma * value);
  }
}

double p_theory(double t, double e0, double e1, double MAXTIME, char *type)
{
  /*
    Hight temperature Regime Optimal robability
  */
  double K, k;
  if (strcmp(type, "up") == 0 || strcmp(type, "down") == 0)
  {
    k = (e1 - e0) / (MAXTIME + 2);
    K = (1. / 4) * k*k;
    return 0.5 - (1. / 4) * (e0 + 2*sqrt(K)*t);
  }
  else
  {
    return 0;
  }
}

double power(double a, int b)
{
  int res=1, abso;
  abso = (int) absolute(b);
  for (int i; i<abso; i++)
  {
    if (b < 0)
    {
      res *= 1/a;
    }
    else
    {
      res *= a;
    }
  }
  return res;
}

double absolute(double a)
{
  if (a > 0)
  {
    return a;
  }
  else
  {
    return -a;
  }
}
