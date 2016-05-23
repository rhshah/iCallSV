"""
Created on Mar 19, 2015
Description: Convert VCF to targetSeqView
#Example:
SampleDesc    Chr1    Start1    End1    LeftSideSegDup    Chr2    Start2    End2    RightSideSeqDup    ValidationStatus    Sample    SplitsSample
Ramos    15    22462315    22462465    TRUE    14    106467050    106467150    TRUE    Failed PCR    1320KB0009MultipleAlnsort.bam    1320KB0009.bam
::Input::
sampleName: Name of the sample that has the structural abberations
sampleBamName: Name of the bam file.
sampleSplitBaName: Name of the split bam file (Use bam file if you dont have split bam file)
vcfFile: Input Delly VCF file for the conversion
outputDir: Directory to write the output file
outputFileName: Name of the output File
::Output::
outputFile: TargetSeqView format text file for a given vcf file.
@author: Ronak H Shah
"""

import vcf
import checkparameters as cp
import logging

logger = logging.getLogger(__name__)

def Convert2targetSeqView(
        sampleName,
        sampleBamName,
        sampleSplitBamName,
        vcfFile,
        outputDir,
        outputFileName):
    
    logger.info("Convert2targetSeqView: Will convert vcf to targetSeqVie format")
    cp.checkFile(vcfFile)
    cp.checkDir(outputDir)
    logger.info(
        "Convert2targetSeqView: All Input Parameters look good. Lets convert to tab-delimited file")
    vcf_reader = vcf.Reader(open(vcfFile, 'r'))
    outputFile = outputDir + "/" + outputFileName
    outputHandle = open(outputFile, "w")
    outputHandle.write(
        "SampleDesc\tChr1\tStart1\tEnd1\tLeftSideSegDup\tChr2\tStart2\tEnd2\tRightSideSeqDup\tValidationStatus\tSample\tSplitsSample\n")
    for record in vcf_reader:
        (chrom1,
         start1,
         start2,
         chrom2,
         contype,
         str1,
         str2) = (None for i in range(7))
        chrom1 = record.CHROM
        start1 = record.POS
        if("END" in record.INFO):
            start2 = record.INFO['END']
        if("CHR2" in record.INFO):
            chrom2 = record.INFO['CHR2']
        if("CT" in record.INFO):
            contype = record.INFO['CT']
        (startCT, endCT) = contype.split("to")
        outputHandle.write(
            sampleName +
            "\t" +
            str(chrom1) +
            "\t" +
            str(int(start1)-50) +
            "\t" +
            str(int(start1)+50) +
            "\tFALSE\t" +
            str(chrom2) +
            "\t" +
            str(int(start2)-50) +
            "\t" +
            str(int(start2)+50) +
            "\tFALSE\tFailed PCR\t" +
            str(sampleBamName) +
            "\t" +
            str(sampleSplitBamName) +
            "\n")
    outputHandle.close()
    logger.info("Convert2targetSeqView: Finished conversion of Vcf file to targetSeqView file format.")
    logger.info("Convert2targetSeqView: Output can be found: %s", outputFile)
    return(outputFile)
# Test module
