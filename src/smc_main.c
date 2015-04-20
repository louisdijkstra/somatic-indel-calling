/*
 * Copyright (C) 2015 Louis Dijkstra
 * 
 * This file is part of somatic-indel-calling
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * smc_main.c
 * 
 * Takes in an observations file (see bin/extract-observations.py) and 
 * determines for every indel in the file the maximum a posteriori estimates
 * of the variant allele frequency (vaf) of the healthy (h) and cancer (c) 
 * cells. 
 * 
 * Author: Louis Dijkstra
 * E-mail: dijkstra@cwi.nl
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"

#include "smc_likelihood.h"

int main(int argc, char *argv[])
{
	// parse the command line   
	char input_filename[256]; 
	parameters p ; 
	parse_arguments(&p, input_filename, argc, argv) ;

	if (p.verbose == 1) {
		print_parameters(&p) ;
	} 

	FILE *fp = fopen(input_filename, "r") ; 
	if (fp == NULL) {
		printf("Could not open file %s\n", input_filename) ; 
		exit(EXIT_FAILURE) ; 
	}

	// set the global parameters used for the likelihood function (see sm_likelihood.h)
	ALPHA 		= p.alpha; 
	MU_H 		= p.mu_h; 
	SIGMA_H 	= p.sigma_h; 
	RATE_H		= p.rate_h; 
	MU_C 		= p.mu_c; 
	SIGMA_C 	= p.sigma_c; 
	RATE_C		= p.rate_c; 
	EPS_A 		= p.eps_a; 
	EPS_P 		= p.eps_p_del; 
	MAX_ITER 	= p.max_iter; 
	EPSABS 		= p.epsabs; 
	N_PANELS 	= p.n_panels; 

	/* Allocate memory */
	variant v ;
	v.chromosome = malloc(LENGTH_CHROM_BUFFER * sizeof(char)) ; 
	data D ; 
	
	double mle_h_vaf, mle_c_vaf, max_logl; 
	double p_somatic, p_germline, p_not_present ; 
	
	size_t status ; 

	while(1) {
		// read in the data 
		status = obtainVariant(fp, &v, &p) ; 
		if (status == END_OF_FILE_REACHED) {
			break ; 
		}
	
		// print info on the VCF record
		printf("%c\t%s\t%zd\t%zd", v.type, v.chromosome, v.position, v.length);		
	
		if (status == VALID_VARIANT) { // variant is of the right type
			D = obtainData(fp, &p) ; 
			if (v.type == '+') {
				D.delta = -1.0*v.length ; 
				EPS_P = p.eps_p_ins; 
			} else {
				D.delta = v.length ;
				EPS_P = p.eps_p_del; 
			}

			if (unlikelyInsertSize(D)) {
				D = removeUnlikelyInsertSizes(D) ; 
			}
	
			if (uniqueGlobalMaximumExists(D) == 1) {
				max_logl = computeMLE(&mle_h_vaf, &mle_c_vaf, D) ; 
				determinePosteriorProbabilities(&p_somatic, &p_germline, &p_not_present, D) ; 
				printf("\t%f\t%f\t%f\t%f\t%f\t%f\n", mle_h_vaf, mle_c_vaf, max_logl, p_somatic, p_germline, p_not_present) ; 
			} else {
				printf("\t.\t.\t.\t.\t.\t.\n") ; 
			}
		} else {
			skipDataEntry(fp) ; 
			printf("\t.\t.\t.\t.\t.\t.\n") ; 
		}	
	}
	
    	fclose(fp);
    	exit(EXIT_SUCCESS) ; 
}
