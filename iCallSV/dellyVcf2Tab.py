'''
Created on Mar 18, 2015
Description: This module converts the Delly Vcf file having tumor normal, to tab-delimited format for input to iAnnotateSV
@author: Ronak H Shah
'''
'''
::Input::
vcfFile: Input vcf file to convert
outputFileName: Name of the output file
OutputDir: Directory for output file
::Output::
outputFile: Tab-delimited file containing:
chr1: Its the chromosome name for first break point [1,2,3,4,5,6,7 etc..],
pos1: Its the chromosome loaction for first break point [1-based],
str1: Its the read direction for the first break point [0=top/plus/reference, 1=bottom/minus/complement],
chr2: Its the chromosome name for second break point [1,2,3,4,5,6,7 etc..],
pos2: Its the chromosome loaction for second break point [1-based],
str2: Its the read direction for the second break point [0=top/plus/reference, 1=bottom/minus/complement],
'''
import vcf
import checkparameters as cp


def vcf2tab(vcfFile, outputFileName, outputDir):
    cp.checkFile(vcfFile)
    cp.checkDir(outputDir)
    print "dellyVcf2Tab: All Input Parameters look good. Lets convert to tab-delimited file\n"
    vcf_reader = vcf.Reader(open(vcfFile, 'r'))
    outputFile = outputDir + "/" + outputFileName
    outputHandle = open(outputFile, "w")
    outputHandle.write("chr1\tpos1\tstr1\tchr2\tpos2\tstr2\n")
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
        if((int(startCT) == 3) and (int(endCT) == 3)):
            str1 = 0
            str2 = 0
        elif((int(startCT) == 3) and (int(endCT) == 5)):
            str1 = 0
            str2 = 1
        elif((int(startCT) == 5) and (int(endCT) == 3)):
            str1 = 1
            str2 = 0
        elif((int(startCT) == 5) and (int(endCT) == 5)):
            str1 = 1
            str2 = 1
        else:
            print "dellyVcf2Tab: The connection type (CT) given in the vcf file is incorrect.CT:", contype, "\n"
        outputHandle.write(
            str(chrom1) +
            "\t" +
            str(start1) +
            "\t" +
            str(str1) +
            "\t" +
            str(chrom2) +
            "\t" +
            str(start2) +
            "\t" +
            str(str2) +
            "\n")
    outputHandle.close()
    print "dellyVcf2Tab: Finished conversion of Vcf file to tab-delimited file.\n"
    print "dellyVcf2Tab: Output can be found: ", outputFile, "\n"
    return(outputFile)

#Test Module
#vcf2tab("/dmp/hot/shahr2/IMPACT/Test/SVtest//35462375-T_bc44_jmp.stdfilter.vcf","35462375-T_bc44_jmp.stdfilter.tab" ,"/dmp/hot/shahr2/IMPACT/Test/SVtest/" )
