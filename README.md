Somatic Indel Calling
=====================

This repository contains the code/scripts for discovering somatic mutations by using a Bayesian latent variable model that links (next-generation sequencing) reads and their alignments to the (latent) allele frequencies of indels in both the healthy/control as the cancer/disease sample. The pipeline is written in both C and Python. 

***

## Installation 

### Dependencies 

The compilation of the C code requires the following libraries to be installed:  

* The _GNU scientific library_ (GSL - see http://www.gnu.org/software/gsl), in particular:
	* `gsl_math.h`,
	* `gsl_min.h` and 
	* `gsl_errno.h`. 

* The _GMP library_ for arbitrary precision arithmetic (see https://gmplib.org). 

and

* _CMake_ (see http://www.cmake.org)

The project depends for Python on the following packages: 

* _PySAM_ (see https://code.google.com/p/pysam/) for working with BAM/SAM files and
* _PyVCF_ (see https://github.com/jamescasbon/PyVCF) for working with VCF files. 

PySAM requires the installation of 

* _SAMtools_ (see http://samtools.sourceforge.net)

### Installation instructions 

In order to compile the C code in the folder `src/`, type in the main directory: 

```
	$ cmake . 
	$ make
	$ make install 
```

The executable `sm_caller` will be placed in the `bin/` folder together with the Python scripts. CMake will check automatically whether the GSL and GMP libraries are installed. 

***

## Usage 

First, make sure to set the parameters in `parameters.txt` (in the main directory). The parameters are currently set to their default values. 

In order to use the somatic indel calling pipeline, type in the main directory: 

	$ python call-somatic-variants.py <vcf-file> <healthy-bam> <cancer-bam>  

where 

* `<vcf-file>` - sorted VCF file with potential somatic variants.
	
* `<healthy-bam>` - sorted and indexed BAM file of the healthy/control sample.
	
* `<cancer-bam>` - sorted and indexed BAM file of the cancer/disease sample. 

(One could use the option `-v` for a more elaborate output). A VCF file containing all indels deemed somatic is stored in the `results/` directory. See for a more elaborate description of the pipeline the following section. 

## Pipeline

The main script `call-somatic-mutations.py` executes the following steps: 

1. It reads in the parameters used throughout the entire pipeline from the `parameters.txt` file. 

2. It runs the `extract-observations.py` script; for every indel in the given VCF, the script collects all relevant alignments from the healthy/control and cancer/disease BAMs. (Relevant meaning that the alignment overlaps with the indel in question.) Both insert sizes and splits are taken into consideration.
The output is stored as a `.raw-observations` file (see the File formats section) in the `intermediate-results/` folder. 

3. After having collected all the relevant data, the program `bin/sm_caller` is called. This program makes the actual calls and determines the posterior probabilities of the variants to be somatic, germline or absent. The results are stored as a `.calls` file (see the File formats section) in the `intermediate-results/` folder. 

4. Finally, the script `bin/calls-to-vcf.py` takes in both the `.calls`-file of the previous step and the original VCF file and outputs a VCF file (in format v4.1) in the `results/` directory. Note that it only outputs the variants deemed somatic. 

*** 

## Directory structure 

The repository consists of the following directories: 

* `bin/` - contains the Python scripts and executables used for calling the somatic mutations. 

* `include/` - contains the header files for the C-code.

* `intermediate-results/` - will contain all the intermediate results, i.e., the so-called 'raw observations' and the 'calls' file. See for a description of their format Section File formats.

* `results/` - will contain the VCF file with the somatic calls.

* `python/` - contains Python code that is re-used throughout some of the scripts in the `bin`-folder.

* `src/` - contains the C-code for the somatic mutation caller. 

In addition, the repository contains in the main directory the file `parameters.txt`. In this file, one can set the main parameters that play a role for calling the somatic mutations. 

***

## File formats

This section contains a description of the two (novel) file formats used in this projects: `.raw-observations` and `.calls`. 

### .raw-observations

This file format is used as a summary of the relevant data
from a VCF file and two BAM files; one from the healthy and one
from the cancer sample. Every indel in the given VCF file is 
represented by 9 lines. The first line contains information
on the variant and is formatted as follows (tab-delimited): 

	<type>	<chr>	<pos>	<length>

where `<type>` is '+' in case of an insertion and '-' for a 
deletion. The chromosome and the variant's postion are given by 
`<chr>` and `<pos>`. The indel's length is given by `<length>`. 

The next four lines represent the data of the healthy sample. The 
later four are for the cancer sample. The first line represents in
both cases the observed insert sizes. The second the associated 
alignment probabilities. The third line contains the split read 
observations, i.e., an observation is 1 when a read contained a split 
that supported the presence of the indel and a 0 otherwise. The 
fourth line contains the alignment probabilities associated with 
the split read observations. 

### .calls

Every line represents one variant and consists of 10 columns in total (tab-delimited):

1. _type_ - the symbol `+` represents an insertion and `-` a deletion;

2. _chromosome_ - the chromosome the variant is on;

3. _position_ - its position (the same as in the original VCF file);

4. _length_ - the indel's length;

5. estimate of the healthy variant allele frequency. This value can either be 0.0 (the variant is absent on both chromosomes), 0.5 (heterozygous for the variant) and 1.0 (homozygous for the variant); 

6. estimate of the cancer variant allele frequency. This value lies between 0 and 1 and reflects the proportion of cancer haplotypes that harbour the variant in question;

7. the maximum loglikelihood associated with the estimates; 

8. the posterior probability of the variant to be SOMATIC; 

9. the posterior probability of the variant to be GERMLINE;  

10. the posterior probability of the variant to be NOT PRESENT.    

***

## Contact

Louis Dijkstra

__E-mail__: louisdijkstra (at) gmail.com
