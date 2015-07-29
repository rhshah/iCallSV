'''
Created on Mar 18, 2015
Description: This module will sort bam file by name
@author: Ronak H Shah
'''
import os
import sys
import pysam


def sortBam(inputBam, outputBamName, outputDir):
    print "sortbamByReadName: Trying to sort BAM file by Read Name"
    if(os.path.isdir(outputDir)):
        print "sortbamByReadName: The output directory ", outputDir, " exists"
    outputFile = outputDir + "/" + outputBamName

    if(os.path.isfile(inputBam)):
        try:
            pysam.sort("-n", inputBam, outputFile)
        except IndexError as err:
            print "Index error({0}): {1}".format(err.errno, err.strerror)
        except IOError as err:
            print "I/O error({0}): {1}".format(err.errno, err.strerror)
    else:
        print inputBam, " File doesnot exists !!"
        sys.exit()
    print "sortbamByReadName: Finished sorting BAM file by Read Name.\n"
    return(outputFile)

#sortBam('/home/shahr2/M15-2555.recal.bam', "M15-2555.recal.NSORT","/home/shahr2/")
