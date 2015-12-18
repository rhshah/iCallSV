'''
Created on December 18, 2015
Description: This module will merge multiple vcf file with same headers
@author: Ronak H Shah
'''
'''
::Input::
vcfFiles : List of VCF Files to merge in list data structure
mergedVCF: Name of the merged vcf to output

::Output::
It is a merged vcf file 
'''
import sys
import vcf
vcf_header = vcf.Reader(filename=vcfFiles[1])
vcf_output = vcf.Writer(mergedVCF,vcf_header)
for vcffile in vcfFiles:
	vcf_reader = vcf.Reader(open(vcffile, 'r'))
	for each in vcf_reader:
		vcf_output.write_record(each)
return()