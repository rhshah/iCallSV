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
   
   
iCallSV is a Python library and command-line software toolkit to call structural aberrations from Next Generation DNA sequencing data. Behind the scenes it uses Delly2 to do structural variant calling. It is designed for use with hybrid capture, including both whole-exome and custom target panels, and
short-read sequencing platforms such as Illumina.

Citation
========

We are in the process of publishing a manuscript describing iCallSV as part of the Structural Variant Detection framework.
If you use this software in a publication, for now, please cite our website `iCallSV <http://github.com/rhshah/iCallSV>`_.

Required Packages
=================
We require that you install:

:pandas: `v0.16.2 <http://pandas.pydata.org/>`_
:biopython: `v1.65 <http://biopython.org/wiki/Main_Page>`_
:pysam: `v0.8.4 <https://pypi.python.org/pypi/pysam>`_
:pyvcf: `0.6.7 <https://pypi.python.org/pypi/PyVCF>`_
:Delly: `v0.7.3 <https://github.com/tobiasrausch/delly>`_
:targetSeqView: `master <https://github.com/Eitan177/targetSeqView>`_
:iAnnotateSV: `v1.0.5 <https://github.com/rhshah/iAnnotateSV/tree/1.0.5>`_


Required Data Files
===================

This files are given in the ``data`` folder inside iCallSV, they are uploaded using `git-lfs <https://git-lfs.github.com/>`_ and need to be downloaded with `git-lfs <https://git-lfs.github.com/>`_

:blacklistRegionsFile: Tab-delimited file wihout header having black listed regions.

	:Example:
	
		7	140498077	5	175998094
		

:blacklistGenes: Gene listed one per line wihout header that are to be removed 

	:Example:
	
		LINC00486
		
		CNOT4
		
		
:genesToInclude: Gene listed one per line wihout header that are to be kept
	
	:Example:
	
		ALK
		
		BRAF

Configuration File Format
=========================

.. code-block:: sh

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
	####Case Allele Fraction Hotspot ###
	CaseAltFreqHotspot: 0.05
	####Total Case Coverage Hotspot ###
	CaseCoverageHotspot = 5
	####Control Allele Fraction Hotspot###
	ControlAltFreqHotspot = 0
	#### Case Allele Fraction ###
	CaseAltFreq: 0.08
	#### Total Case Coverage ###
	CaseCoverage = 8
	#### Control Allele Fraction ###
	ControlAltFreq = 0
	### Overall Supporting Read-pairs ###
	OverallSupportingReads: 5
	### Overall Supporting Read-pairs Hotspot ###
	OverallSupportingReadsHotspot: 3
	### Overall Supporting splitreads ###
	OverallSupportingSplitReads: 0
	### Overall Supporting splitreads Hotspot ###
	OverallSupportingSplitReadsHotspot: 0
	### Case Supporting Read-pairs ###
	CaseSupportingReads: 2
	### Case Supporting splitreads ###
	CaseSupportingSplitReads: 0
	### Case Supporting Read-pairs Hotspot ###
	CaseSupportingReadsHotspot: 1
	### Case Supporting splitreads Hotspot ###
	CaseSupportingSplitReadsHotspot: 0
	### Control Supporting Read-pairs ###
	ControlSupportingReads: 5
	### Control Supporting Read-pairs Hotspot ###
	ControlSupportingReadsHotspot: 5
	### Control Supporting splitreads ###
	ControlSupportingSplitReads: 5
	### Control Supporting splitreads Hotspot ###
	ControlSupportingSplitReadsHotspot: 5
	### Length of Structural Variant ###
	LengthOfSV: 500
	### Overall Mapping Quality Threshold###
	OverallMapq: 20
	### Overall Mapping Quality Threshold Hotspot ###
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
	



Submodules
==========

iCallSV.FilterDellyCalls module
-------------------------------

.. automodule:: iCallSV.FilterDellyCalls
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.Run_Delly module
------------------------

.. automodule:: iCallSV.Run_Delly
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.Run_iAnnotateSV module
------------------------------

.. automodule:: iCallSV.Run_iAnnotateSV
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.Run_samblaster module
-----------------------------

.. automodule:: iCallSV.Run_samblaster
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.Run_targetSeqView module
--------------------------------

.. automodule:: iCallSV.Run_targetSeqView
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.checkBlackList module
-----------------------------

.. automodule:: iCallSV.checkBlackList
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.checkHotSpotList module
-------------------------------

.. automodule:: iCallSV.checkHotSpotList
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.checkparameters module
------------------------------

.. automodule:: iCallSV.checkparameters
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.combineVCF module
-------------------------

.. automodule:: iCallSV.combineVCF
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.dellyVcf2Tab module
---------------------------

.. automodule:: iCallSV.dellyVcf2Tab
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.dellyVcf2targetSeqView module
-------------------------------------

.. automodule:: iCallSV.dellyVcf2targetSeqView
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.filterAnnotatedSV module
--------------------------------

.. automodule:: iCallSV.filterAnnotatedSV
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.iCallSV module
----------------------

.. automodule:: iCallSV.iCallSV
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.iCallSV_dmp_wrapper module
----------------------------------

.. automodule:: iCallSV.iCallSV_dmp_wrapper
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.launchThreads module
----------------------------

.. automodule:: iCallSV.launchThreads
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.launch_FilterDellyCalls module
--------------------------------------

.. automodule:: iCallSV.launch_FilterDellyCalls
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.launch_Run_Delly module
-------------------------------

.. automodule:: iCallSV.launch_Run_Delly
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.make_analysis_dir module
--------------------------------

.. automodule:: iCallSV.make_analysis_dir
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.makebamindex module
---------------------------

.. automodule:: iCallSV.makebamindex
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.mergeFinalFiles module
------------------------------

.. automodule:: iCallSV.mergeFinalFiles
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.sortbamByCoordinate module
----------------------------------

.. automodule:: iCallSV.sortbamByCoordinate
    :members:
    :undoc-members:
    :show-inheritance:

iCallSV.sortbamByReadName module
--------------------------------

.. automodule:: iCallSV.sortbamByReadName
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: iCallSV
    :members:
    :undoc-members:
    :show-inheritance:
