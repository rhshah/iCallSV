'''
Created on Mar 19, 2015
Description: This module will sort bam file by coordinate
@author: Ronak H Shah
'''
import os
import sys
import pysam


def sortBam(inputBam, outputBamName, outputDir):
    print "sortbamByCoordinate: Trying to sort BAM file by Coordinate"
    if(os.path.isdir(outputDir)):
        print "sortbamByCoordinate: The output directory ", outputDir, " exists"
    outputFile = outputDir + "/" + outputBamName

    if(os.path.isfile(inputBam)):
        try:
            pysam.sort(inputBam, outputFile)
        except IndexError as err:
            print "Index error({0}): {1}".format(err.errno, err.strerror)
        except IOError as err:
            print "I/O error({0}): {1}".format(err.errno, err.strerror)
    else:
        print inputBam, " File doesnot exists !!"
        sys.exit()
    print "sortbamByCoordinate: Finished sorting BAM file by Coordinate.\n"
    return(outputFile)

#sortBam('/home/shahr2/M15-2555.recal.bam', "M15-2555.recal.NSORT","/home/shahr2/")
