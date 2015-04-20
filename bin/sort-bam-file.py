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
import pysam

__author__ = "Louis Dijkstra"

usage = """%prog <bam-file> 

	<bam-file> 		original BAM file

Sorts the BAM file using PySAM (which employs SAMTools). 
If the given file is '<name>.bam', the output is stored under
'<name>.sorted.bam'. 
"""

def main():

	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1
		
	filename = args[0]
	
	if filename[-4:] != '.bam': 
		print("ERROR: file extension should be '.bam'. File is not sorted.")
		return 1 
		
	print("BAM file %s is sorted. Output is stored at %s"%(filename, filename[:-3] + 'sorted.bam'))
	pysam.sort(filename, filename[:-3] + 'sorted')
# 	pysam.index(args[1])

if __name__ == '__main__':
	sys.exit(main())

