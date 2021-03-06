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
from BAMProcessor import * # import the super class 

__author__ = "Louis Dijkstra"

"""
	DefaultBAMProcessor.py contains the functionality to process a given unspecified BAM file. 
	(This used to be the way we processed a BAM file)
"""

class DefaultBAMProcessor(BAMProcessor): 
	"""Class for processing BAM files. Contains all functionality needed for retrieving alignments."""
	def __init__(self, bam_filename, search_range = 5000, primary_alignments_only = False, centerpoints_thres_del = 50, centerpoints_thres_ins = 50, length_thres_del = 20, length_thres_ins = 20):
		BAMProcessor.__init__(self, bam_filename, search_range = search_range, primary_alignments_only = primary_alignments_only)
		self.centerpoints_thres_del 	= centerpoints_thres_del	
		self.centerpoints_thres_ins 	= centerpoints_thres_ins
		self.length_thres_del		= length_thres_del 
		self.length_thres_ins		= length_thres_ins 		

	def determineSupportDeletionSingleAlignment(self, deletion, alignment): 
		i = alignment.pos 
		for (cigar_type, cigar_length) in alignment.cigar: # walk through the cigar string
			if cigar_type == 2: # deletion
				if abs(cigar_length - deletion.length) > self.length_thres_del:
					continue
				centerpoints_split = None
				if cigar_length % 2 == 0: 
					centerpoints_split = [i + cigar_length / 2 - 1, i + cigar_length / 2]
				else:
					centerpoints_split = [i + cigar_length / 2]
				if returnMinimumDifference(centerpoints_split, deletion.centerpoints) <= self.centerpoints_thres_del:
					return 1 
			i += cigar_length
		return 0 

	def determineSupportInsertionSingleAlignment(self, insertion, alignment):
		i = alignment.pos 
		for (cigar_type, cigar_length) in alignment.cigar: # walk through the cigar string
			if cigar_type == 1: # insertion
				if abs(cigar_length - insertion.length) > self.length_thres_ins:
					continue
				centerpoints_split = None
				if cigar_length % 2 == 0: 
					centerpoints_split = [i + cigar_length / 2 - 1, i + cigar_length / 2]
				else:
					centerpoints_split = [i + cigar_length / 2]
				if returnMinimumDifference(centerpoints_split, [insertion.position, insertion.position + 1]) <= self.centerpoints_thres_ins:
					return 1 
			i += cigar_length
		return 0
	

