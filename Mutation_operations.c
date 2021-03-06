/*
 * Copyrigtht 2014 Georgios Karagiannis
 *
 * This file is part of PISAA_ProtAB2D.
 *
 * PISAA_ProtAB2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * PISAA_ProtAB2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISAA_ProtAB2D.  If not, see <http://www.gnu.org/licenses/>.
*/

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

#include <math.h>

#include "RNG.h"
#include "cost_protein2D.h"
#include "Self_adjastment_prosedure.h"

#define ES_MH 0.5

void Mutation_HitAndRun(double *z, double *fz,
				int N_dimension,
				double *theta, double *grid_points, int grid_size,
				double temp, double scl, double *accpr_pop,
				double *z_new){

	/* THIS IS THE HIT AND RUN METROPOLIS HASTINGS UPDATE */

	/* Smith, Robert L. "Efficient Monte Carlo procedures for generating points
	 * uniformly distributed over bounded regions." Operations Research 32.6
	 * (1984): 1296-1308. */

	/*Chen, Ming-Hui, and Bruce Schmeiser. "Performance of the Gibbs, hit-and-run,
	 * and Metropolis samplers." Journal of computational and graphical
	 * statistics 2.3 (1993): 251-272.*/

	int k_old ;
	int k_new ;
	int i;

	double fz_new ;
	double rat ;
	double un ;

	self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

	/* propose a new solution */
	uniformdirectionrng( z_new , N_dimension ) ;
	un = normalrng()*scl *pow(k_old/(grid_size+1.0),ES_MH) ;
	for ( i = 1 ; i <= N_dimension ; i++ )
		z_new[i] = z[i] + z_new[i]*un ;

	fz_new = cost( z_new, N_dimension ) ;

	self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

	rat = -theta[k_new] -fz_new/temp +theta[k_old] + *fz/temp  ;

	/* Accept reject */

	*accpr_pop = ( (rat>0.0) ? 1.0 : exp( rat ) ) ;

	un = uniformrng() ;
	if ( *accpr_pop>un ){
		*fz = fz_new;
		for ( i = 1 ; i <= N_dimension ; i++) z[i] = z_new[i] ;
	}

}

/*K POINT OPERATION*/

void Mutation_Kpoint(double *z, double *fz,
				int N_dimension,
				double *theta, double *grid_points, int grid_size,
				double temp, double scl, double *accpr_pop,
				double *z_new){

	int k_old ;
	int k_new ;
	int i;
	int lab_1 ;
	int lab_2 ;

	double fz_new ;
	double rat ;
	double un ;

	self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

	/* propose a new solution */

	for ( i=1 ; i<=N_dimension ; i++) z_new[i] = z[i] ;

	/*1st dimension*/
	lab_1 = integerrng(1,N_dimension) ;
	un = normalrng()*scl *pow(k_old/(grid_size+1.0),ES_MH) ;
	z_new[lab_1] += un ;

	/*2nd dimension*/
	if ( (uniformrng()<0.5) && (N_dimension>=2) ){
		lab_2 = integerrng(1,N_dimension-1) ;
		if ( lab_2 >= lab_1 ) lab_2++ ;
		un = normalrng()*scl ;
		z_new[lab_2] += un ;
	}

	fz_new = cost( z_new, N_dimension ) ;

	self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

	rat = -theta[k_new] -fz_new/temp +theta[k_old] + *fz/temp  ;

	/* Accept reject */

	*accpr_pop = ( (rat>0.0) ? 1.0 : exp( rat ) ) ;

	un = uniformrng() ;
	if ( *accpr_pop>un ){
		*fz = fz_new;
		for ( i = 1 ; i <= N_dimension ; i++) z[i] = z_new[i] ;
	}

}

/*METROPOLIS RANDOM WALK OPERATION*/

void Mutation_Metropolis(double *z, double *fz,
				int N_dimension,
				double *theta, double *grid_points, int grid_size,
				double temp, double scl, double *accpr_pop,
				double *z_new){

	int k_old ;
	int k_new ;
	int i;

	double fz_new ;
	double rat ;
	double un ;

	self_adj_index_search(&k_old, *fz, grid_points, grid_size) ;

	/* propose a new solution */

	for ( i = 1 ; i <= N_dimension ; i++ )
		z_new[i] = z[i] +normalrng()*scl *pow(k_old/(grid_size+1.0),ES_MH) ;

	fz_new = cost( z_new, N_dimension ) ;

	self_adj_index_search(&k_new, fz_new, grid_points, grid_size) ;

	rat = -theta[k_new] -fz_new/temp +theta[k_old] + *fz/temp  ;

	/* Accept reject */

	*accpr_pop = ( (rat>0.0) ? 1.0 : exp( rat ) ) ;

	un = uniformrng() ;
	if ( *accpr_pop>un ){
		*fz = fz_new;
		for ( i = 1 ; i <= N_dimension ; i++) z[i] = z_new[i] ;
	}

}









