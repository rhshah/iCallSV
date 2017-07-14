iCallSV: Structural Aberration Detection from NGS datasets
================================================================

:Author: Ronak H Shah
:Contact: rons.shah@gmail.com
:Source code: http://github.com/rhshah/iCallSV
:Wiki: http://icallsv.readthedocs.io/en/latest/
:License: `Apache License 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_

.. image:: https://landscape.io/github/rhshah/iCallSV/master/landscape.svg?style=flat
   :target: https://landscape.io/github/rhshah/iCallSV/master
   :alt: Code Health
   

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.184864.svg
	:target: https://doi.org/10.5281/zenodo.184864
   

.. image:: https://codecov.io/gh/rhshah/iCallSV/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/rhshah/iCallSV
   

iCallSV is a Python library and command-line software toolkit to call structural aberrations from Next Generation DNA sequencing data. Behind the scenes it uses Delly2 to do structural variant calling. It is designed for use with hybrid capture, including both whole-exome and custom target panels, and
short-read sequencing platforms such as Illumina.

The filtering process can be observed here:  `Workflow <https://www.draw.io/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1&title=iCallSV_Filters.html#Uhttps%3A%2F%2Fdrive.google.com%2Fuc%3Fid%3D0Bwn1ij1qNCi_cE4xOW1NS0JJaTg%26export%3Ddownload>`_

Citation
========

We are in the process of publishing a manuscript describing iCallSV as part of the Structural Variant Detection framework.
If you use this software in a publication, for now, please cite our website `iCallSV <http://github.com/rhshah/iCallSV>`_.

Note
====

For some reason docstrings is not shown by `Read The Docs <https:read-the-docs.readthedocs.io>`_

So please use these url from `Github <https:github.com>`_ with `Html Preview <https://htmlpreview.github.io/>`_ for each module information:

`Per Module Info <https://htmlpreview.github.io/?https://raw.githubusercontent.com/rhshah/iCallSV/master/docs/_build/html/iCallSV.html>`_ 


Required Packages
=================
We require that you install:

:pandas: `v0.16.2 <http://pandas.pydata.org/>`_
:biopython: `v1.65 <http://biopython.org/wiki/Main_Page>`_
:pysam: `v0.8.4 <https://pypi.python.org/pypi/pysam>`_
:pyvcf: `0.6.7 <https://pypi.python.org/pypi/PyVCF>`_
:Delly: `v0.7.5 <https://github.com/tobiasrausch/delly>`_
:targetSeqView: `master <https://github.com/Eitan177/targetSeqView>`_
:iAnnotateSV: `v1.0.6 <https://github.com/rhshah/iAnnotateSV/tree/1.0.5>`_
:coloredlogs: `v5.2 <https://pypi.python.org/pypi/coloredlogs>`_

Required Data Files
===================

This files are given in the ``data`` folder inside iCallSV.


:BlackListFile: (blacklist.txt) Tab-delimited file wihout header having black listed regions in order:
				
				**chromosome 1, breakpoint 1, chromosome 2, breakpoint 2**
				

	:Example:
	
		7	140498077	5	175998094
		

:BlackListGenes: (blacklistgenes.txt) Gene listed one per line wihout header that are to be removed 

	:Example:
	
		LINC00486
		
		CNOT4
		
:HotspotFile: (hotspotgenes.txt) Tab-delimited file wihout header having hotspot regions in order:
			  
			  **chromosome, start, end, name**
	
	:Example:
	
		2	29416089	30143525	ALK

:GenesToKeep: (genesToInclude.txt) Gene listed one per line wihout header that are to be kept
	
	:Example:
	
		ALK
		
		BRAF
		

Configuration File Format
=========================

.. code-block:: ini
	
	#~~~Template configuration file to run iCallSV~~~#
	#### Path to python executable ###
	[Python]
	PYTHON:
	#### Path to R executable and R Lib ###
	[R]
	RHOME: 
	RLIB: 
	#### Path to delly, bcftools executables and Version of delly (supports only 0.7.3)###
	[SVcaller]
	DELLY:
	DellyVersion:
	BCFTOOLS:
	#### Path to hg19 Referece Fasta file ###
	[ReferenceFasta]
	REFFASTA:
	#### Path to file containing regions to exclude please follow Delly documentation for this ###
	[ExcludeRegion]
	EXREGIONS:
	#### Path to file containing regions to where lenient threshold will be used; and file containing genes to keep ###
	[HotSpotRegions]
	HotspotFile:
	GenesToKeep:
	#### Path to file containing regions/genes to filter ###
	[BlackListRegions]
	BlackListFile:
	BlackListGenes:
	#### Path to samtools executable ###
	[SAMTOOLS]
	SAMTOOLS:
	#### Path to iAnnotateSV.py and all its required files, please follow iAnnotateSV documentation ###
	[iAnnotateSV]
	ANNOSV:
	GENOMEBUILD:
	DISTANCE:
	CANONICALTRANSCRIPTFILE:
	UNIPROTFILE:
	CosmicCensus:
	CosmicFusionCounts:
	RepeatRegionAnnotation:
	DGvAnnotations:
	#### TargetSeqView Parameters ###
	[TargetSeqView]
	CalculateConfidenceScore:
	GENOMEBUILD:
	ReadLength:
	#### Parameters to run Delly ###
	[ParametersToRunDelly]
	MAPQ: 20
	NumberOfProcessors: 4
	[ParametersToFilterDellyResults]
	####Case Allele Fraction Hotspot####
	CaseAltFreqHotspot: 0.05
	####Total Case Coverage Hotspot#####
	CaseCoverageHotspot = 5
	####Control Allele Fraction Hotspot####
	ControlAltFreqHotspot = 0
	####Case Allele Fraction####
	CaseAltFreq: 0.10
	####Total Case Coverage#####
	CaseCoverage = 10
	####Control Allele Fraction####
	ControlAltFreq = 0
	###Overall Supporting Read-pairs ###
	OverallSupportingReads: 5
	###Overall Supporting Read-pairs Hotspot ###
	OverallSupportingReadsHotspot: 3
	###Overall Supporting splitreads ###
	OverallSupportingSplitReads: 0
	###Overall Supporting splitreads Hotspot ###
	OverallSupportingSplitReadsHotspot: 0
	###Case Supporting Read-pairs ###
	CaseSupportingReads: 2
	###Case Supporting splitreads ###
	CaseSupportingSplitReads: 0
	###Case Supporting Read-pairs Hotspot ###
	CaseSupportingReadsHotspot: 1
	###Case Supporting splitreads Hotspot ###
	CaseSupportingSplitReadsHotspot: 0
	###Control Supporting Read-pairs ###
	ControlSupportingReads: 3
	###Control Supporting Read-pairs Hotspot ###
	ControlSupportingReadsHotspot: 3
	###Control Supporting splitreads ###
	ControlSupportingSplitReads: 3
	###Control Supporting splitreads Hotspot ###
	ControlSupportingSplitReadsHotspot: 3
	###Length of Structural Variant###
	LengthOfSV: 500
	###Overall Mapping Quality Threshold###
	OverallMapq: 20
	###Overall Mapping Quality Threshold Hotspot###
	OverallMapqHotspot: 5
	


Quick Usage
===========

.. code-block:: sh

	python iCallSV.py -sc /path/to/template.ini -abam /path/to/casebamFile -bbam /path/to/controlbamFile -aId caseID -bId controlId -o /path/to/output/directory -op prefix_for_the_output_files


.. code-block:: sh
	
	> python iCallSV.py -h
	
	usage: iCallSV.py [-h] [-v] [-V] -sc config.ini -abam caseBAMFile.bam -bbam
	                  controlBAMFile.bam -aId caseID -bId controlID -o
	                  /somepath/output -op TumorID

	iCallSV.iCallSV -- wrapper to run iCallSV package

	  Created by Ronak H Shah on 2015-03-30.
	  Copyright 2015-2016 Ronak H Shah. All rights reserved.

	  Licensed under the Apache License 2.0
	  http://www.apache.org/licenses/LICENSE-2.0

	  Distributed on an "AS IS" basis without warranties
	  or conditions of any kind, either express or implied.

	USAGE

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --verbose         set verbosity level [default: True]
	  -V, --version         show program's version number and exit
	  -sc config.ini, --svConfig config.ini
	                        Full path to the structural variant configuration
	  -abam caseBAMFile.bam, --caseBam caseBAMFile.bam
	                        Full path to the case bam file
	  -bbam controlBAMFile.bam, --controlBam controlBAMFile.bam
	                        Full path to the control bam file
	  -aId caseID, --caseId caseID
	                        Id of the case to be analyzed, this will be the sub-
	                        folder
	  -bId controlID, --controlId controlID
	                        Id of the control to be used, this will be used for
	                        filtering variants
	  -o /somepath/output, --outDir /somepath/output
	                        Full Path to the output dir.
	  -op TumorID, --outPrefix TumorID
	                        Id of the Tumor bam file which will be used as the
	                        prefix for output files

Running on SGE or LSF
=====================

.. sidebar:: Note:
			
	For both SGE and LSF you need to provide total number of cores based on the number of threads you have assinged to delly installation using **OMP_NUM_THREADS**. 

.. sidebar:: Note:

	For example: if you set **OMP_NUM_THREADS** as `export OMP_NUM_THREADS=3` then you need to set total number of cores to be 13 (12 + 1 extra as buffer) so for each of the Delly program it utilizes 3 cores. Here I use pythons multiprocessing module to launch delly, so all four programs would be launch as seprate process utilizing number of threads given to them but setting the **OMP_NUM_THREADS**
	
SGE
---

.. code-block:: sh
	
	qsub -q some.q -N iCallSV_JobName -o iCallSV.stdout -e iCallSV.stderr -V -l h_vmem=6G,virtual_free=6G -pe smp 13 -wd /some/path/to/working/dir -sync y  -b y python iCallSV.py -sc template.ini -bbam control.bam -abam case.bam -aId caseID -bId controlID -op outputPrefix -o  /some/path/to/output/dir -v 

LSF
---

.. code-block:: sh

	bsub -q some.q -J iCallSV_JobName -o iCallSV.stdout -e iCallSV.stderr -We 24:00 -R "rusage[mem=20]" -M 30 -n 13 -cwd /some/path/to/working/dir "python iCallSV.py -sc template.ini -bbam control.bam -abam case.bam -aId caseID -bId controlID -op outputPrefix -o  /some/path/to/output/dir -v"

						
Utilities
=========

Running iCallSV on MSK-IMPACT Pools
-----------------------------------

**This is only for MSK-IMPACT internal samples**

.. code-block:: sh
	
	> python iCallSV_dmp_wrapper.py -h
	
	usage: iCallSV_dmp_wrapper.py [options]

	Run iCallSV on selected pools using MSK data

	optional arguments:
	  -h, --help            show this help message and exit
	  -fl folders.fof, --folderList folders.fof
	                        Full path folders file of files.
	  -qc /some/path/qcLocation, --qcLocation /some/path/qcLocation
	                        Full path qc files.
	  -b /some/path/bamlocation, --bamLocation /some/path/bamlocation
	                        Full path bam files.
	  -P /somepath/python, --python /somepath/python
	                        Full path Pyhton executables.
	  -icsv /somepath/iCallSV.py, --iCallSV /somepath/iCallSV.py
	                        Full path iCallSV.py executables.
	  -conf /somepath/template.ini, --iCallSVconf /somepath/template.ini
	                        Full path configuration file to run iCallSV
	  -q all.q or clin.q, --queue all.q or clin.q
	                        Name of the SGE queue
	  -qsub /somepath/qsub, --qsubPath /somepath/qsub
	                        Full Path to the qsub executables of SGE.
	  -t 5, --threads 5     Number of Threads to be used to run iCallSV
	  -v, --verbose         make lots of noise [default]
	  -o /somepath/output, --outDir /somepath/output
	                        Full Path to the output dir.
	  -of outputfile.txt, --outDir outputfile.txt
					  	    Name of the final output file.


Taking the iCallSV and chechking for processed transcript/cDNA in samples
-------------------------------------------------------------------------

.. code-block:: sh
	
	> python check_cDNA_contamination.py -h
	usage: check_cDNA_contamination.py [options]

	Calculate cDNA contamination per sample based of the Structural Variants
	Pipeline result

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --verbose         make lots of noise [default]
	  -s SVfile.txt, --svFile SVfile.txt
	                        Location of the structural variant file to be used
	  -o cDNA_contamination, --outputFileName cDNA_contamination
	                        Full path name for the output file
	

	
	