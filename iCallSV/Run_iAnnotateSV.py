"""
Run_iAnnotate
~~~~~~~~~~~~~

:Description: This module will run iAnnotateSV package

"""
'''
Created on Mar 18, 2015
Description: This module will run iAnnotateSV package
@author: Ronak H Shah

::Input::
python : Location for the python executable.
iAnnotateSV : Location of the wrapper iAnnotateSV package (iAnnotateSV.py)
build : Which human reference file to be used, hg18,hg19 or hg38
inputTabFile : Tab-Delimited Input FIle compatible with iAnnotateSV package.
outputPrefix: Prefix of the output files/DIR with Annotations and images
outputDir : Name of the output directory where the outputPrefix will be written
uniprotFile: Location for ucsc uniprot file
cosmicFile: Location for cosmic census file
repeatregionFile: Location for repeat region file
dgvFile: Location for database of Genomic Variants file
::Output::
It is a tab-delimited file with annotation in them
'''

import os
import sys
from subprocess import Popen
import shlex
import time
from datetime import date, timedelta
import checkparameters as cp
import logging
import coloredlogs

logger = logging.getLogger('iCallSV.Run_iAnnotateSV')
coloredlogs.install(level='DEBUG')

def run(
        python,
        iAnnotateSV,
        build,
        distance,
        canonicalTranscriptFile,
        uniprotFile,
        cosmicFile,
        cosmicCountsFile,
        repeatregionFile,
        dgvFile,
        inputTabFile,
        outputPrefix,
        outputDir):
    """
    This module will run iAnnotateSV package.

    :param str python: Location for the python executable.
    :param str iAnnotateSV: Location of the wrapper iAnnotateSV package (iAnnotateSV.py)
    :param str build: Which human reference file to be used, hg18,hg19 or hg38
    :param str inputTabFile: Tab-Delimited Input FIle compatible with iAnnotateSV package.
    :param str outputPrefix: Prefix of the output files/DIR with Annotations and images
    :param str outputDir: Name of the output directory where the outputPrefix will be written
    :param str uniprotFile: Location for ucsc uniprot file
    :param str cosmicFile: Location for cosmic census file
    :param str cosmicCountsFile: Location for cosmic fusion counts file
    :param str repeatregionFile: Location for repeat region file
    :param str dgvFile: Location for database of Genomic Variants file
    :return: str of the output file
    :rtype: str

    """

    start_time = time.time()
    cp.checkDir(outputDir)
    cp.checkFile(iAnnotateSV)
    cp.checkFile(inputTabFile)
    cp.checkFile(python)
    cp.checkInt(distance, "Distance for extending the promoter region")
    cp.checkEmpty(build, "Which human reference file to be used, hg18,hg19 or hg38")
    cp.checkFile(canonicalTranscriptFile)
    cp.checkFile(uniprotFile)
    cp.checkFile(cosmicFile)
    cp.checkFile(cosmicCountsFile)
    cp.checkFile(repeatregionFile)
    cp.checkFile(dgvFile)
    logger.info("Run_iAnnotateSV: All input parameters look good. Lets run the package.")
    myPid = os.getpid()
    day = date.today()
    today = day.isoformat()
    logger.info("Run_iAnnotateSV: ProcessID:%s, Date:%s", myPid, today)
    outputFile = outputDir + "/" + outputPrefix + "_Annotated.txt"
    cmd = python + " " + iAnnotateSV + " -r " + build + " -i " + inputTabFile + " -o " + outputDir + " -ofp " + outputPrefix + " -d " + str(
        distance) + " -c " + canonicalTranscriptFile + " -rr " + repeatregionFile + " -cc " + cosmicFile + " -cct " + cosmicCountsFile + " -dgv " + dgvFile + " -v -p -u " + uniprotFile
    args = shlex.split(cmd)
    logger.info("Run_iAnnotateSV: Command that will be run: %s", cmd)
    # Remove if the file exists
    if(os.path.isfile(outputFile)):
        os.remove(outputFile)
    proc = Popen(args)
    proc.wait()
    retcode = proc.returncode
    if(retcode >= 0):
        end_time = time.time()
        totaltime = str(timedelta(seconds=end_time - start_time))
        logger.info(
            "Run_iAnnotateSV: We have finished running iAnnotateSV for %s using local machine",
            inputTabFile)
        logger.info("Run_iAnnotateSV Duration: %s", totaltime)
    else:
        logger.info(
            "Run_iAnnotateSV: iAnnotateSV is either still running on local machine or it errored out with return code %d for %s",
            retcode,
            inputTabFile)
        sys.exit()
    return(outputFile)
