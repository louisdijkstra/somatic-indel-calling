#!/usr/bin/env python

"""
Copyright (C) 2015 Louis Dijkstra

This file is part of somatic-indel-calling

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function, division
from optparse import OptionParser
import os
import sys

__author__ = "Louis Dijkstra"

usage = """%prog <vcf-file> <healthy-bam> <tumour-bam> 

	<vcf-file> 	Sorted VCF file with potential somatic variants
	<healthy-bam>	Sorted and indexed BAM file of the healthy/
					control sample
	<tumour-bam>	Sorted and indexed BAM file of the tumour/
					case sample. 

(See 'sort-vcf-file.py', 'sort-bam-file.py' and 'index-bam-file.py' 
for sorting and indexing the VCF and BAM files.)

This program runs the entire pipeline to determine for all indels 
in the given VCF, whether or not they are somatic or not somatic 
(i.e., germline/not present) on the basis of the provided BAM files 
of the healthy/control and tumour/case sample. 

NOTE: this program assumes that it is placed in the main directory
of the repository.
"""

def execute(command, description = None): 
	if description != None: 
		barrier = '*' * 50 
		print("\n%s\n%s\n%s"%(barrier, description, barrier))
	print("EXECUTING ---> \t%s\n"%command)
	os.system(command)

def readParameters():
		"""Reads in the parameters from parameters.txt"""
		parameters_dict = dict() ; # allocate memory 
		for line in open('parameters.txt', 'r'):
			if line[0] == '$': # line in which parameter is set
				parameter_label, value = line.split('=')
				parameters_dict[parameter_label.strip()[1:]] = value.strip() 
		return parameters_dict 

def correctInputAndRepository(input_vcf_filename, bam_healthy_filename, bam_tumour_filename): 
	"""Checks the correctness of the input and repository. Return True when correct, False otherwise"""

	# check whether the input files exist
	if not os.path.exists(input_vcf_filename): 
		print("ERROR: the VCF file %s cannot be found"%input_vcf_filename) 
		return False
	
	if not os.path.exists(bam_healthy_filename): 
		print("ERROR: the BAM file %s cannot be found"%bam_healthy_filename) 
		return False
	
	if not os.path.exists(bam_tumour_filename): 
		print("ERROR: the BAM file %s cannot be found"%bam_tumour_filename) 
		return False
	

	# check whether the intermediate-results and results directories exist
	if not os.path.isdir(os.getcwd() + '/intermediate-results/'):
		print("ERROR: repository must contain a folder named 'intermediate-results'") 
		return False
	
	if not os.path.isdir(os.getcwd() + '/results/'):
		print("ERROR: repository must contain a folder named 'results'") 
		return False

	# read parameters from parameters.txt
	if not os.path.exists('parameters.txt'): 
		print("ERROR: parameter file 'parameters.txt' must be in main directory") 
		return False

	return True 

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-v", action="store_true", dest="verbose", default=False, 
				  		help="Verbose.")
	(options, args) = parser.parse_args()
	
	if (len(args)!=3):
		parser.print_help()
		return 1	
	
	input_vcf_filename 	= os.path.abspath(args[0])
	bam_healthy_filename 	= os.path.abspath(args[1])
	bam_tumour_filename 	= os.path.abspath(args[2])
	
	# check whether the input and repository are correct:
	if not correctInputAndRepository(input_vcf_filename, bam_healthy_filename, bam_tumour_filename):
		return 1 

	parameters = readParameters() # read the parameters from parameters.txt
	
	# files: 
	basename 		= os.path.basename(input_vcf_filename)[:-4] 	
	observations_filename	= os.getcwd() + '/intermediate-results/' + basename + '.raw-observations'
	calls_filename		= os.getcwd() + '/intermediate-results/' + basename + '.calls'
	output_vcf_filename	= os.getcwd() + '/results/' + basename + '.somatic.vcf'

	if options.verbose:
		print("Parameter settings")
		print("==================\n")
		print("label\tvalue")
		print("-----\t-----") 
		for key, value in sorted(parameters.items()): 
			print('%s\t%s'%(key, value))
		print("-----\t-----\n") 
	
	### EXTRACTING OBSERVATIONS FROM VCF AND BAMs ###
	command = "python bin/extract-observations.py %s %s %s %s > %s"%(
				parameters['ALIGNER'],
				input_vcf_filename,
				bam_healthy_filename, 
				bam_tumour_filename,
				observations_filename,
			)
	
	if not os.path.exists(observations_filename): 		
		execute(command, description = "Extracting raw observations")
	else: 
		print("\nThe raw observations file %s already exists.\n"%observations_filename) 
	
	### CALLING THE VARIANTS ###
	command = "bin/sm_caller -B -a %s -e %s -E %s -y %s %s > %s"%(
				parameters['ALPHA'],
				parameters['EPS_P_DEL'],
				parameters['EPS_P_INS'],
				parameters['EPS_A'], 
				observations_filename,
				calls_filename
			)
	
	if not os.path.exists(calls_filename): 
		execute(command, description = "Calling")
	else: 
		print("\nThe calls file %s already exists.\n"%calls_filename) 
	
	### CREATING A VCF FILE WITH THE SOMATIC CALLS ###
	command = "python bin/calls-to-vcf.py -b %s -s DreamChallenge5.0 %s %s %s"%(
				parameters['BETA'], 
				input_vcf_filename,
				calls_filename,
				output_vcf_filename
			) 

	if not os.path.exists(output_vcf_filename): 
		execute(command, description = "Creating the VCF file with somatic calls.")
	else: 
		print("\nThe final VCF file %s already exists.\n"%output_vcf_filename) 
	

if __name__ == '__main__':
	sys.exit(main())

