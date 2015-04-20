#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import vcf

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')
from Indel import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <vcf-file> <calls-file> <new-vcf-file> 

	<vcf-file> 	tabix-indexed VCF file 
	<calls-file>	output file generated by sm_caller 
	<new-vcf-file> 	annoted VCF file

Annotates the given VCF file with the results generated by the sm_caller 
application. Outputs somatic calls only! 
"""

header = """##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

def printHeader(outputfile, source="POSOM"):
	print("##fileformat=VCFv4.1", file=outputfile)
	print("##source=%s"%source, file=outputfile)
	print("%s"%header, file=outputfile)

def determinePosteriorProbabilityThreshold(calls_filename, beta, verbose=False): 
	"""Determines the posterior probability threshold for calling
	   something somatic by maximizing the posterior expected
	   F_beta score for classifying items as somatic/not somatic."""

	if verbose: print("Reading in all somatic posterior probabilities") 

	calls_file = open(calls_filename, 'r')

	n = 0 
	somatic_posterior_probabilities = [] 
	for line in calls_file: # get the posterior probabilities 
		values 	= line.split('\t')
		if len(values) != 10:
			break 
		if values[7].strip() != '.':	
			somatic_posterior_probabilities.append(float(values[7]))
			n += 1

	if verbose: 
		print("DONE Reading in all somatic posterior probabilities") 
		print("Sorting all probabilities...")

	somatic_posterior_probabilities.sort(reverse=True)

	if verbose: 
		print("DONE Sorting all probabilities...")
		print("Maximizing expected posterior F-score (beta is set to %f)"%beta)

	# initialization
	beta2		= beta**2		
	S_n 		= beta2*sum(somatic_posterior_probabilities)
	S_k 		= somatic_posterior_probabilities[0]
	max_f_score 	= (1.0 + beta2) * S_k / (S_n + 1.0) # F-beta score for k = 1
	threshold	= S_k

	for k in range(2,n+1):
		S_k += somatic_posterior_probabilities[k-1] 
		f_score = (1.0 + beta2) * S_k / (S_n + k)
		#print("%d\t%f\t%f"%(k, f_score, somatic_posterior_probabilities[k-1]))
		if f_score > max_f_score:
			max_f_score = f_score
			threshold = somatic_posterior_probabilities[k-1] 

	if verbose: print("Threshold is set to %f"%threshold) 
	
	calls_file.close()

	return threshold 

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-b", action="store", dest="beta", type=float, default=1.0,
                      		help="Beta value of the F-score metric. (Default = 1.0)")
	parser.add_option("-s", action="store", dest="source", default="POSOM",
                      		help="Source used for the VCF file header. (Default = POSOM)")
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Prints regularly how many variants have been processed.")
	(options, args) = parser.parse_args()
	
	if (len(args)!=3):
		parser.print_help()
		return 1

	vcf_reader		= vcf.Reader(open(os.path.abspath(args[0])))
	calls_filename 		= os.path.abspath(args[1])
	new_vcf_file 		= open(os.path.abspath(args[2]), 'w')

	printHeader(new_vcf_file, source=options.source)

	posterior_threshold = determinePosteriorProbabilityThreshold(calls_filename, options.beta, verbose=options.verbose)

	calls_file = open(calls_filename, 'r') 
	variant_description = '' 
	n = 0 
	for line in calls_file: 
		n += 1 
		values 	= line.split('\t')
		
		if len(values) != 10: 
			print("ERROR: line number %d does not contain 10 values as required. The line in question is:\n\t%s\n"%(n, line))
			return 1 

		if options.verbose and n % 1000 == 0: 
			print("Processed %d variants"%n)

		position 	= int(values[2])

		call = 'UNKNOWN'
		post_prob_call = 'NONE'
		h_vaf, c_vaf, max_logl = 'NONE', 'NONE', 'NONE' 
		p_somatic, p_germline, p_not_present = 'NONE', 'NONE', 'NONE'
	
		if values[4].strip() != '.':	h_vaf 		= float(values[4])
		if values[5].strip() != '.': 	c_vaf 		= float(values[5])
		if values[6].strip() != '.':	max_logl	= float(values[6])
		if values[7].strip() != '.':	p_somatic 	= float(values[7])
		if values[8].strip() != '.':	p_germline 	= float(values[8])
		if values[9].strip() != '.':	p_not_present 	= float(values[9])
		
		if not (p_somatic is 'NONE' and p_germline is 'NONE' and p_not_present is 'NONE'):
			if p_somatic >= posterior_threshold: 
				call = 'SOMATIC'
			else: 
				call = 'NOT-SOMATIC'  

		vcf_record = vcf_reader.next() 
		while (vcf_record.POS != position):
			vcf_record = vcf_reader.next()

		if call == 'SOMATIC': 
			svlen = returnIndelLength(vcf_record)
			if isDeletion(vcf_record): 
				end = vcf_record.POS + svlen 
				variant_description = "%s\t%s\t.\t%s\t%s\t.\tPASS\tSOMATIC;SVTYPE=DEL;END=%d;SVLEN=%d"%(
					vcf_record.CHROM,
					vcf_record.POS,
					vcf_record.REF,
					vcf_record.ALT[0],
					vcf_record.POS + svlen,
					-1 * svlen)
			elif isInsertion(vcf_record):
				variant_description =  "%s\t%s\t.\t%s\t%s\t.\tPASS\tSOMATIC;SVTYPE=INS;END=%d;SVLEN=%d"%(
					vcf_record.CHROM,
					vcf_record.POS,
					vcf_record.REF,
					vcf_record.ALT[0],
					vcf_record.POS,
					svlen)
			print("%s"%variant_description, file=new_vcf_file)

	calls_file.close()
	new_vcf_file.close()	 

if __name__ == '__main__':
	sys.exit(main())
