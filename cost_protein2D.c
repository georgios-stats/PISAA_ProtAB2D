

/*
 * Georgios Karagiannis
 * Postdoctoral research associate
 * Department of Mathematics, Purdue University
 * 150 N. University Street
 * West Lafayette, IN 47907-2067, USA
 *
 * Telephone: +1 (765) 496-1007
 *
 * Email: gkaragia@purdue.edu
 *
 * Contact email: georgios.stats@gmail.com
*/


#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include <stdio.h>

#ifndef INFINITY
	#include <float.h>
	#define INFINITY DBL_MAX
#endif

static int *seq=NULL ;

void get_data(char file_name[], int N_monomer, int *N_dimension){

	int i ;
	FILE *ins;

	*N_dimension = N_monomer -2 ;

	seq = ivector(1,N_monomer) ;

	ins = fopen(file_name, "r") ;

	if (ins==NULL) {
		printf("No data set\n") ;
		abort() ;
	}

	for(i=1; i<=N_monomer; i++)
		fscanf(ins, " %d", &seq[i]);

	fclose(ins);

}

/*void cost_bounds(double *z_min, double *z_max, int i){
	if (i > 0){
		*z_min = (double) -INFINITY ;
		*z_max = (double) INFINITY ;

	}
}*/

double cost(double *z, int N_dimension){

	int N_monomer = N_dimension+2 ;
	/*double alpha[N_monomer+1],D[N_monomer+1][3], sum1,sum2,r1,r2, a;*/
	double sum1,sum2,r1,r2, a;
	double *alpha, **D ;
	int i, j, k, l;

	double pi = 3.14159265358979 ;

	alpha = dvector(0,N_monomer) ;
	D = dmatrix(0,N_monomer,0,2) ;

/*	double z_min, z_max ;
	for (i=1; i<=N_dimension; i++){
		cost_bounds(&z_min, &z_max, i) ;
		if ( z[i]<z_min || z[i]>z_max ) return (double) INFINITY ;
	}*/

	D[1][1]=D[1][2]=0.0; D[2][1]=1.0; D[2][2]=0.0;
	alpha[1]=0.0;
	for(i=2; i<=N_monomer-1; i++){
		a=z[i-1]/(2.0*pi);
		k=floor(a+0.5);
		alpha[i]=z[i-1]-k*2.0*pi;
		z[i-1]=alpha[i];
		D[i+1][1]=D[i][1]+cos(alpha[i]);
		D[i+1][2]=D[i][2]+sin(alpha[i]);
	}

	for(sum1=0.0,i=1; i<=N_monomer-2; i++){
		r1=(D[i+1][1]-D[i][1])*(D[i+2][1]-D[i+1][1]);
		r2=(D[i+1][2]-D[i][2])*(D[i+2][2]-D[i+1][2]);
		sum1+=0.25*(1-r1-r2);
	  }

	sum2=0.0;
	for(i=1; i<=N_monomer-2; i++)
	  for(j=i+2; j<=N_monomer; j++){
		  r2=(D[i][1]-D[j][1])*(D[i][1]-D[j][1])+(D[i][2]-D[j][2])*(D[i][2]-D[j][2]);
		  a=(1+seq[i]+seq[j]+5*seq[i]*seq[j])/8.0;
		  sum2+=4.0*(pow(r2,-6)-a*pow(r2,-3));
		}

	free_dvector(alpha,0,N_monomer) ;
	free_dmatrix(D,0,N_monomer,0,2) ;

	return sum1+sum2;
}
