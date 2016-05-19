Note:This Is Work In Progress. Use this at your own RISK.

iCallSV: Structural Aberration Detection from NGS datasets
================================================================

:Author: Ronak H Shah
:Contact: rons.shah@gmail.com
:Source code: http://github.com/rhshah/iCallSV
:License: `Apache License 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_

.. image:: https://landscape.io/github/rhshah/iCallSV/master/landscape.svg?style=flat
   :target: https://landscape.io/github/rhshah/iCallSV/master
   :alt: Code Health
   
   
iCallSV is a Python library and command-line software toolkit to call structural aberrations from Next Generation DNA sequencing data. It is designed for use with hybrid capture, including both whole-exome and custom target panels, and
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
:Delly: `v0.6.1 and 0.7.3 <https://github.com/tobiasrausch/delly>`_


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
	