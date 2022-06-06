#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_sf_gamma.h>

#include "lattice_generator.h"
#include "alias_table.h"
#include "find_root.h"

//=============================================================================

//#define L 64
// #define sg 0.1

// #define ti 10
// #define tf ti
// #define dt 0.1

 #define tp 'C'

// #define hi 0.00
// #define hf hi
// #define dh 0.1

#define pi 4*atan(1.0)

gsl_rng *gslran;

//=============================================================================
double energy (int L, int Spin[L],int M,double h,double dist[L])
{
  int i,j;
  double E = 0;
  for (i = 0; i < L; ++i) for (j = i+1; j < L; ++j ) {
    E += -Spin[i]*Spin[j]*dist[j-i-1];
  } 
  E += -h*M;
  return E;
}
//=============================================================================

void metro (int L, int Spin[L], double beta, double h, double dist[L])
{
  int MCini = 3*L;
  int mea   = 100000;
  int MCd   = 20;
  int MCtot = MCini+mea*MCd;
  int i,j;
  
  double sum = 0,sum_E = 0,sum_M = 0,sum_M2 = 0;
  
  int M = 0;
  for (i = 0; i < L; ++i) M += Spin[i];
  double E = energy(L, Spin,M,h,dist);
   
  int imc;
  for (imc = 0; imc < MCtot; ++imc)  {
    int ipas;
    for (ipas = 0; ipas < L; ++ipas) {
      int I = gsl_rng_uniform_int(gslran,L);

      double E_1 = 0;
      for (i = 0; i < L; ++i) {
        if ( I != i) {
          E_1 += -Spin[I]*Spin[i]*dist[abs(I-i)-1];
        }
      }
      E_1 += -h*Spin[I];
      
      double delta = -2*E_1;
      
      if (delta <= 0) {
        E += delta;
        Spin[I] = -Spin[I];
      } else {
        double x = gsl_rng_uniform(gslran);
        if (x < exp(-beta*delta)) {
          E += delta;
          Spin[I] = -Spin[I];
        } 
      }
    } 
    
    if ((imc > MCini) && (MCd*(imc/MCd)) == imc) {
      
      M = 0;
      for (i = 0; i < L; ++i) M += Spin[i];
      E = energy(L, Spin,M,h,dist);

      sum += 1.;
      sum_E += E;
      sum_M += abs(M);
      sum_M2 += M*M;
          
    } 
    
  }

  double e  = sum_E/(sum*L);
  double m  = sum_M/(sum*L);
  double m2 = sum_M2/(sum*L*L);
  double x  = L*(m2-m*m);
  double T  = 1.0/beta;
  
  printf("METRO:   L= %d | T= %2.4f | h= %2.2f | e= %1.5f | m= %1.5f | x= %4.5f \n", L,T,h,e,m,x);
}
//=============================================================================

void wolff (int L, int Spin[L], double beta, double h,double dist[L])
{
  int MCini = 3*L;
  int mea   = 100000;
  int MCd   = 20;
  int MCtot = MCini+mea*MCd;
  int i,j;
  
//   int sigma = 10*sg;
//   char out_data[60];
//   FILE *file_data;
//   sprintf(out_data,"WOLFF_1D_%05d_%02d_PBC_TL.dat",L,sigma);
//   file_data = fopen(out_data,"a+");
  
  double sum = 0,sum_E = 0,sum_M = 0,sum_M2 = 0;
  
  int imc;
  for (imc = 0; imc < MCtot; ++imc)  {
    int I_clu[L];
    for ( i = 0; i < L; ++i) I_clu[i] = 0;
    int I = gsl_rng_uniform_int(gslran,L);
    int S = Spin[I];
    I_clu[0] = I;
    Spin[I] = - S;
    int N_clu = 0;
    int B_clu = 0;

    while (B_clu < N_clu+1) {
      int k;
      for ( k = B_clu ; k < N_clu+1; ++k) {
        I = I_clu[k];
	int ai;
        for ( ai = 0; ai < L; ++ai) {
	  if ( I != ai) {
            int In = ai;   
            if (Spin[In] == S) {
	      double prob = 1.0-exp(-2.0*beta*dist[abs(I-ai)-1]);
              if (gsl_rng_uniform(gslran) < prob) {
                Spin[In] = -Spin[In];
                N_clu += 1;
                I_clu[N_clu] = In;
        } } } }
        B_clu += 1;
    } }
    
    //external source energy
    //double delta = 0;
    //for (i = 0; i < N_clu+1; i++) delta += Spin[I_clu[i]];
//    double delta = 2*h*beta*(N_clu+1)*S;
   // printf("delta= %4.4f N_clu= %d \n",delta,S*(N_clu+1));
    
  //  if (delta > 0) {
    //  if (gsl_rng_uniform(gslran) > exp(-delta)) {
        //printf("second delta= %4.4f \n",delta);
       // for (i = 0; i < N_clu+1; i++) {
      //    Spin[I_clu[i]] = S;
       //   printf("unchanged \n");
       // }
     // }
   // }
    
    if ((imc > MCini) && (MCd*(imc/MCd)) == imc) {
      
      int M = 0;
      for (i = 0; i < L; ++i) M += Spin[i];
      double E = energy(L, Spin,M,h,dist);
      //double E = 0;
      
      sum += 1.;
      sum_E += E;
      sum_M += abs(M);
      sum_M2 += M*M;
       
    } 
    
  }

  double e  = sum_E/(sum*L);
  double m  = sum_M/(sum*L);
  double m2 = sum_M2/(sum*L*L);
  double x  = L*(m2-m*m);
  double T  = 1.0/beta;
  
  printf("WOLFF:   L= %d | T= %2.4f | h= %2.2f | e= %1.5f | m= %1.5f | x= %4.5f \n", L,T,h,e,m,x);
  //fprintf(file_data,"L= %d T= %2.4f h= %2.2f e= %1.5f m= %1.5f x= %4.5f \n", L,T,h,e,m,x);
  //fclose(file_data); 
}
//=============================================================================
double modes (int L, int *Spin, int m)
{
  int i;
  double chi,mag_cos,mag_sin,kk;
   mag_cos = 0;
   mag_sin = 0;
   kk = 2*pi/(L);
   for (i = 0; i < L; ++i) {
     mag_sin += sin((double) m*kk*(i+1.0) )*Spin[i];
     mag_cos += cos((double) m*kk*(i+1.0) )*Spin[i];
   }
      chi = (mag_cos*mag_cos+mag_sin*mag_sin)/(L);
  return chi;
}//=============================================================================

void luijten (int L, int *Spin, double *bp, int *a, double beta, double Jtot, double *dist,double h)
{
  int MCini = 3*L;
  int mea   = 100000;
  int MCd   = 20;
  int MCtot = MCini+mea*MCd;
  int i,j,k;
  
//   int sigma = 10*sg;
//   char out_data[60];
//   FILE *file_data;
//   sprintf(out_data,"Luijten_1D_%05d_%02d_PBC_T%1c.dat",L,sigma,tp);
//   file_data = fopen(out_data,"a+");
  
  double sum = 0,sum_E = 0,sum_E2 = 0,sum_M = 0,sum_M2 = 0;
  
  //setting prob
  double p[L];
  for (i = 1; i < L; ++i) p[i] = 1.0-exp(-2.0*beta*dist[i-1]);
  p[0] = 0;
  double P[L];
  for (i = 0; i < L; ++i) P[i] = p[i];  
  for (i = 0; i < L; ++i) {
    for (k = 0; k < i; ++k)  P[i] *= (1.0-p[k]);
  }
  double C[L];
  for (i = 0; i < L; ++i) C[i] = 0;
  for (i = 0; i < L; ++i) {
    for (k = 0; k < i+1; ++k)  C[i] += P[k];
  }
  
 // for (i = 0; i < L; ++i) printf("i= %d, p= %6.6f, P= %6.6f, C= %6.6f \n",i,p[i],P[i],C[i]);
  //starting MC
  int imc;
  for (imc = 0; imc < MCtot; ++imc)  {
   // printf("sweep %d \n",imc);
    int I_clu[L];
    for ( i = 0; i < L; ++i) I_clu[i] = 0;
    int I = gsl_rng_uniform_int(gslran,L);
    int S = Spin[I];
    I_clu[0] = I;
    Spin[I] = - S;
    int N_clu = 0;
    int B_clu = 0;
    //printf(" I= %d \n",I);
    while (B_clu < N_clu+1) {
      int k;
      for ( k = B_clu ; k < N_clu+1; ++k) {
        I = I_clu[k];
	
	
	
	
	//printf(" I= %d \n",I);
	//int cs = gsl_rng_uniform_int(gslran,L-1)+1;
	//int I2 = I+cs;
	//if (I2>L-1) I2 = I2-L;
 	//printf(" cs= %d, I2= %d \n",cs,I2);
	
	//if (Spin[I2] == S) {
         //     Spin[I2] = -Spin[I2];
          //    N_clu += 1;
           //   I_clu[N_clu] = I2;
 	  //   printf(" %d added \n",I2);
	   // } 
 	//getchar(); 
	double g = 0; 
	//double g = gsl_rng_uniform(gslran); 
        // g = C[cs]+(1.0-C[cs])*g;	
//  	printf(" second g= %4.4f, previous cs= %d \n",g,cs);
//  	getchar(); 
	
	//while ((cs < L)&&(g<C[L-1])) {
	  while (g<C[L-1]) {
	  int mm = 0;
	  while (g > C[mm]) mm += 1;
	  int cs = mm;
//  	 printf(" next cs= %d \n",cs);
 	// getchar(); 
	  int I2 = I+cs;
	  if (I2>L-1) I2 = I2-L;
 	// printf(" I2= %d, cs= %d \n",I2,cs);
            if (Spin[I2] == S) {
              Spin[I2] = -Spin[I2];
              N_clu += 1;
              I_clu[N_clu] = I2;
	    //  printf(" %d added \n",I2);
	      
	    } 
	g = gsl_rng_uniform(gslran);     
	g = C[cs]+(1.0-C[cs])*g;	
//  	printf("  next g= %4.4f \n",g);
//         getchar(); 
	}
	//getchar();
        B_clu += 1;
    } }
    
    
    
    if ((imc > MCini) && (MCd*(imc/MCd)) == imc) {
        
      int M = 0;
      for (i = 0; i < L; ++i) M += Spin[i];
     double E = energy(L, Spin,M,h,dist);
      //double E = 0;
      
      
      
      sum += 1.;
      //sum_E += sum_mu;
      //sum_E2 += sum_mu*sum_mu;
      sum_E += E;
      sum_E2 += E*E;
      sum_M += abs(M);
      sum_M2 += M*M;
      
     // fprintf(file_data,"L= %d h= %2.2f E= %6.6f M= %d\n", L,h,E,M);     
    } 
    
  }
  
  sum_E /= sum; sum_E2 /= sum; sum_M /= sum; sum_M2 /= sum; 
 // double e = (Jtot-sum_E/(beta))/(L);
 // double e2 = (Jtot-sum_E/(beta))/(L)-h*m;
//  double c = (sum_E2-sum_E*sum_E-sum_E)/(L);
  double e = sum_E/(L);
//  double c = beta*beta*(sum_E2-sum_E*sum_E)/(L*L);
  double m  = sum_M/(L);
  double m2 = sum_M2/(L*L);
  double x  = L*(m2-m*m);
  double T  = 1.0/beta;
  
  printf("LUIJTEN: L= %d | T= %2.4f | h= %2.2f | e= %1.5f | m= %1.5f | x= %4.5f \n", L,T,h,e,m,x);
//   fprintf(file_data,"%d %2.4f %2.4f %1.5f %4.5f %1.5f %4.5f %4.5f %4.5f %16.5f\n", L,T,h,e,c,m,x,mk1,xk1,xi);
//   fclose(file_data); 
}
//=============================================================================

void todo (int L, int *Spin, double *bp, int *a, double beta, double Jtot, double *dist,double h)
{
  int MCini = 3*L;
  int mea   = 10000;
  int MCd   = 20;
  int MCtot = MCini+mea*MCd;
  int i,j,k;
  
   int sigma = 10*0.1;
   char out_data[60];
   FILE *file_data;
   sprintf(out_data,"data/TODO_1D_%05d_%02d_PBC_%6.6f.csv",L,sigma,1.0/beta);
   file_data = fopen(out_data,"a+");
   fprintf(file_data,"K,E,M\n");
  
  double sum = 0,sum_E = 0, sum_E_k = 0, sum_E2 = 0,sum_M = 0,sum_M2 = 0;
  
  int imc;
  for (imc = 0; imc < MCtot; ++imc)  {
    int clu[L];
    for (k = 0; k < L; ++k) clu[k] = -1;
    int sum_mu = 0;
    double lam = 2*Jtot*beta;
    int K = gsl_ran_poisson(gslran,lam);
    int ms;
    for (ms = 0; ms < K; ++ms) {
      int I1 = gsl_rng_uniform_int(gslran,L);
      int cs = gsl_rng_uniform_int(gslran,L-1);
      if (gsl_rng_uniform(gslran) > bp[cs]) cs = a[cs];
      int I2 = I1+cs+1;
      if (I2 > L-1) I2 = I2-L;
      if (Spin[I1]*Spin[I2] == 1) {
        sum_mu += 1;
        int R1 = find_root (I1,clu);
        int R2 = find_root (I2,clu);
        if (R1 != R2) {
          if (clu[R1] >= clu[R2]) {
            clu[R2] = clu[R2] + clu[R1];
            clu[R1] = R2;
          } else {
            clu[R1] = clu[R1] + clu[R2];
            clu[R2] = R1;
          }
        }
      }
    }  
    for (k = 0; k < L; ++k) {
      if (clu[k] < 0) {
        if (gsl_rng_uniform(gslran) < 0.5) {
          Spin[k] = -Spin[k];
        }
      }
    }

     //printf("MCS= %d | \n",imc);
     //for (k = 0; k < L; ++k) printf("Spin= %d, Clus= %d,\n",k,clu[k]);
//      for (k = 0; k < L; ++k) {
//        if (clu[k] < 0) {
	// printf("Spin= %d, Clus= %d,\n",k,clu[k]);
//          if (gsl_rng_uniform(gslran) < 0.5) {
// 	   Spin[k] = -Spin[k];
	   //double delta = 2*h*beta*clu[k]*Spin[k];
           //if (delta > 0) {
             //if (gsl_rng_uniform(gslran) > exp(-delta)) {
               //Spin[k] = -Spin[k];
             //}   
           //}
//          }
//        }
//      }
    for (k = 0; k < L; ++k) {
      if (clu[k] >= 0) {
        int R = clu[k];
        while (clu[R] >= 0) R = clu[R];
        Spin[k] = Spin[R];
      }
    }
    if ((imc > MCini) && (MCd*(imc/MCd)) == imc) {
        
      int M = 0;
      for (i = 0; i < L; ++i) M += Spin[i];
      double E = energy(L, Spin,M,h,dist);
      int E_k = sum_mu;
      
      
      
      sum += 1.;
      sum_E += E;
      sum_E_k += E_k;
      //sum_E2 += sum_mu*sum_mu;
     // sum_E += E;
      sum_E2 += E*E;
      sum_M += abs(M);
      sum_M2 += M*M;
      
      
//      fprintf(file_data,"L= %d h= %2.2f E= %6.6f M= %d\n", L,h,E,M);
      fprintf(file_data,"%6.6f,%d,%d\n",E,M,E_k);
    }
    
  }
  
  sum_E /= sum; sum_E2 /= sum; sum_M /= sum; sum_M2 /= sum; 
  double e = (Jtot-sum_E/(beta))/(L);
 // double e2 = (Jtot-sum_E/(beta))/(L)-h*m;
//  double c = (sum_E2-sum_E*sum_E-sum_E)/(L);
 // double e = sum_E/(L);
  //double c = beta*beta*(sum_E2-sum_E*sum_E)/(L*L);
  double m  = sum_M/(L);
  double m2 = sum_M2/(L*L);
  double x  = L*(m2-m*m);
  double T  = 1.0/beta;
  
  printf("TODO:    L= %d | T= %2.4f | h= %2.2f | e= %1.5f | m= %1.5f | x= %4.5f \n", L,T,h,e,m,x);
//  fprintf(file_data,"%d %2.4f %2.4f %1.5f %4.5f %1.5f %4.5f %4.5f %4.5f %16.5f\n", L,T,h,e,c,m,x,mk1,xk1,xi);
  fclose(file_data);
}
//=============================================================================


int simulation(int L, double sg, double T_ini, double T_fin, double dT, double h_ini, double h_fin, double dh)
{
//   int sigma = 10*sg;
//   char out_data[60];
//   FILE *file_data;
//   sprintf(out_data,"TODO_1D_%05d_%02d_PBC_T%1c.dat",L,sigma,tp);
//   file_data = fopen(out_data,"w+");
//   fprintf(file_data,"# L    T    h    e    m    x \n");
//   fclose(file_data); 
  
  gslran = gsl_rng_alloc(gsl_rng_mt19937);
//  srand(time(NULL));
  int rseed = 1234567890 +1;
  gsl_rng_set(gslran, rseed);
  int i,j;
  
  double Jtot = 0;
  double dist[L];
  double pow(double x, double y);
  for (i = 0; i < L-1; ++i) {
   // dist[i] = pow((double) i+1, -(1+sg));
    double hz1 =  gsl_sf_hzeta((double)1+sg,(double)(i+1)/(L));
    double hz2 =  gsl_sf_hzeta((double)1+sg,(double)(L-i-1)/(L));
    dist [i] = (hz1+hz2)*pow((double)L,-(1+sg));
    Jtot += dist[i];
  }
  Jtot = Jtot*L/2;
  
  int a[L];
  for (i = 0; i < L; ++i) a[i] = -1;
  double bp[L];
  alias_table(L, a,bp,Jtot,dist);
  double hii = h_ini;
  int hde = (int)((h_fin-h_ini)/dh)+1;
  //for (j = 0; j < 20; ++j) {
    
    
    double h = hii;//+j*dh;
    hii *= 1.5;
    int tde = (int)((T_fin-T_ini)/dT)+1;
    for (i = 0; i < tde; ++i) {
      double T = T_ini+i*dT;
      double beta = 1.0/T;
      
      int Spin[L];
      
     //lattice (Spin);
     // metro (Spin,beta,h,dist);
    //    clock_t start1 = clock();
     // lattice (Spin);
     // wolff (Spin,beta,h,dist);
      //   double u1 = (((double)clock() - start1) / CLOCKS_PER_SEC);
      // printf("Time elapsed for L = %d: %5.9f s\n", L,u1);
    //  lattice (Spin);
     // luijten (Spin,bp,a,beta,Jtot,dist,h);
      
    //    clock_t start2 = clock();

//      lattice (L, Spin);
      lattice_generator(L, Spin, gslran);
      todo (L, Spin,bp,a,beta,Jtot,dist,h);
    //   double u2 = (((double)clock() - start2) / CLOCKS_PER_SEC);
      //  printf("Time elapsed for L = %d: %5.9f s\n", L,u2);
      
   // }
  }
  return L;
}
