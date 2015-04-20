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
import vcf
import sys
import operator

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')
from Indel import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <vcf-file> <output-vcf-file>
 
	<vcf-file> 			original VCF file
	<output-vcf-file>	filtered VCF records will be stored here
	
Filters a given VCF file. Type '%prog -h' for the filter options. 	
"""

def fallsInInterval(value, min_value, max_value): 
	if min_value != None and value < min_value:
		return False
	if max_value != None and value > max_value:
		return False
	return True

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("--deletions", action="store_true", dest="deletions", default=False, 
				help = "Filter for deletions")
	parser.add_option("--indels", action="store_true", dest="indels", default=False, 
				help = "Filter for indels")
	parser.add_option("--insertions", action="store_true", dest="insertions", default=False, 
				help = "Filter for insertions")
	parser.add_option("--NOT", action="store_true", dest="negation", default=False, 
				help = "Outputs those records that do NOT satisfy the given constraints. E.g., '--deletions --NOT' returns all records that are not deletions.") 
	parser.add_option("--snps", action="store_true", dest="snps", default=False, 
				help = "Filter for SNPs")
	parser.add_option("-k", action="store", dest="min_length", default=None, type=int, 
				help="Minimal length of the indels (Default = all indels are outputed)")
	parser.add_option("-l", action="store", dest="max_length", default=None, type=int, 
				help="Maximal length of the indels (Default = all indels are outputed)")
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
				help="Verbose. Prints regularly how many variants have been processed.")
	parser.add_option("-x", action="store", dest="chromosome", default=None, 
				help="Processes only this chromosome. (Default = all are processed)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=2):
		parser.print_help()
		return 1

	if options.indels: 
		options.deletions = True
		options.insertions = True 
	
	vcf_filename 		= os.path.abspath(args[0])
	vcf_output_filename	= os.path.abspath(args[1])

	vcf_reader 		= vcf.Reader(open(vcf_filename))
	vcf_writer 		= vcf.Writer(open(vcf_output_filename, 'w'), vcf_reader)

	n, n_passed_filter = 0, 0 
	satisfies_filter = True
	for vcf_record in vcf_reader:
		
		n += 1
		if options.verbose and n % 10000 == 0: 
			print('Processed %d variants'%n)
		
		if options.chromosome == None or options.chromosome == vcf_record.CHROM: 
			satisfies_filter = False
			is_deletion 	= isDeletion(vcf_record) 
			is_insertion 	= isInsertion(vcf_record) 
			is_snp 		= isSNP(vcf_record) 

			if options.snps and is_snp: 
				satisfies_filter = True
			if options.deletions and is_deletion: 
				if fallsInInterval(returnIndelLength(vcf_record), options.min_length, options.max_length):
					satisfies_filter = True
			if options.insertions and is_insertion: 
				if fallsInInterval(returnIndelLength(vcf_record), options.min_length, options.max_length):
					satisfies_filter = True

		if (not options.negation) and satisfies_filter:
			n_passed_filter += 1
			vcf_writer.write_record(vcf_record)
		if options.negation and (not satisfies_filter):
			n_passed_filter += 1
			vcf_writer.write_record(vcf_record)	
		
	vcf_writer.close()

	print("Processed %d variants in total. %d variants passed the filter."%(n, n_passed_filter))

if __name__ == '__main__':
	sys.exit(main())

