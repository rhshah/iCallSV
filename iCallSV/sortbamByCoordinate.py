"""
sortbamByCoordinate
~~~~~~~~~~~~~~~~~~~

:Description: This module will sort bam file by coordinate
"""
'''
Created on Mar 19, 2015
Description: This module will sort bam file by coordinate
@author: Ronak H Shah
'''

import os
import sys
import pysam
import logging
import coloredlogs

logger = logging.getLogger('iCallSV.sortbamByCoordinate')
coloredlogs.install(level='DEBUG')

def sortBam(inputBam, outputBamName, outputDir):
    logger.info("sortbamByCoordinate: Trying to sort BAM file by Coordinate")
    if(os.path.isdir(outputDir)):
        logger.info("sortbamByCoordinate: The output directory %s exists", outputDir)
        outputFile = outputDir + "/" + outputBamName
    else:
        logger.info("sortbamByCoordinate:The output directory %s does not exists !!", outputDir)
        sys.exit()

    if(os.path.isfile(inputBam)):
        try:
            pysam.sort(inputBam, outputFile)
        except IndexError as err:
            exception = "Index error({0}): {1}".format(err.errno, err.strerror)
            logger.info("%s", exception)
        except IOError as err:
            exception = "I/O error({0}): {1}".format(err.errno, err.strerror)
            logger.info("%s", exception)
    else:
        logger.info("sortbamByCoordinate: bam File %s does not exists !!", inputBam)
        sys.exit()
    logger.info("sortbamByCoordinate: Finished sorting BAM file by Coordinate.")
    return(outputFile)
