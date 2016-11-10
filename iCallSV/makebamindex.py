"""
makebamindex
~~~~~~~~~~~~

:Description: Use PySAM to make bam index

"""
# Use PySAM to make bam index
import pysam
import os
import sys
import logging
import coloredlogs

logger = logging.getLogger('iCallSV.makebamindex')
coloredlogs.install(level='DEBUG')

def MakeIndex(bamFile):
    """
    This will make bam index if not there for a bam file using pysam.

    :param str bamFile: Path to bam file
    :return: None
    :rtype: None

    """
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
