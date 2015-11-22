"""
Created on Mar 16, 2015.

@author: Ronak H Shah

"""
# Use PySAM to make bam index
import pysam
import os
import sys


def MakeIndex(bamFile):
    print "makebamindex: Trying to make index for bam file\n"
    if(os.path.isfile(bamFile)):
        try:
            pysam.index(bamFile)
        except IndexError as err:
            print "Index error({0}): {1}".format(err.errno, err.strerror)
        except IOError as err:
            print "I/O error({0}): {1}".format(err.errno, err.strerror)
    else:
        print bamFile, " File doesnot exists !!"
        sys.exit()

# MakeIndex('/home/shahr2/M15-2555.recal.bam')
