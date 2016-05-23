"""
Created on Mar 16, 2015.

@author: Ronak H Shah

"""
# Use PySAM to make bam index
import pysam
import os
import sys
import logging


def MakeIndex(bamFile):
    logger = logging.getLogger(__name__)
    logger.info("makebamindex: Trying to make index for bam file")
    if(os.path.isfile(bamFile)):
        try:
            pysam.index(bamFile)
        except IndexError as err:
            exception = "Index error({0}): {1}".format(err.errno, err.strerror)
            logger.info("%s", exception)
        except IOError as err:
            exception = "I/O error({0}): {1}".format(err.errno, err.strerror)
            logger.info("%s", exception)
    else:
        logger.info("Bam File %s does not exists", bamFile)
        sys.exit()

# MakeIndex('/home/shahr2/M15-2555.recal.bam')
