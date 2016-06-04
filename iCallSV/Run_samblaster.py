"""
Run_samblaster
~~~~~~~~~~~~~~

:Description: This module will run samblaster for extracting discordant and spit reads in sam format

"""
'''
Created on Mar 20, 2015
Description: This module will run samblaster for extracting discordant and spit reads in sam format
@author: Ronak Shah
::Input::
samtools: To stream bam to sam into samblaster
samblaster: Location of the samblaster executables
bamFile: ReadName sorted bamFile
discordantFileName : Name of the discordant SAM file
splitFileName : Name of the split SAM file
outputDir : name the output directory
::Output::
discordantFile : SAM file containing discordant read entries
splitFile : SAM file containing split read entries
'''

import os
import sys
from subprocess import Popen
import shlex
import time
from datetime import date, timedelta
import checkparameters as cp
import logging

logger = logging.getLogger(__name__)

def run(samtools, samblaster, bamFile, discordantFileName, splitFileName, outputDir):
    start_time = time.time()
    print "We will now run samblaster to extract from BAM file discordant and spit reads\n"
    myPid = os.getpid()
    day = date.today()
    today = day.isoformat()
    cp.checkDir(outputDir)
    cp.checkFile(samblaster)
    cp.checkFile(bamFile)
    cp.checkEmpty(discordantFileName, "Name of the Discordant Read SAM File output by samblaster")
    cp.checkEmpty(splitFileName, "Name of the split Read SAM File output by samblaster")
    print "Run_samblaster: All the input parameters look good for running samblaster\n"
    print "Run_samblaster: ProcessID:", myPid, " Date:", today, "\n"
    discordantFile = outputDir + "/" + discordantFileName
    splitFile = outputDir + "/" + splitFileName
    cmd = samtools + " view -h " + bamFile + " | " + samblaster + \
        " -a -e --maxSplitCount 10 -d " + discordantFile + " -s " + splitFile + " -o /dev/null"
    print "Run_samblaster: Command that will be run\n", cmd, "\n"
    #args = shlex.split(cmd)
    proc = Popen(cmd, shell=True)
    proc.wait()
    retcode = proc.returncode
    if(retcode >= 0):
        end_time = time.time()
        print "Run_samblaster: We have finished running samblaster for ", bamFile, " using local machine.\n"
        print "Run_samblaster Duration:", str(timedelta(seconds=end_time - start_time)), "\n"
    else:
        print "Run_samblaster: samblaster is either still running on local machine or it errored out with return code", retcode, " for", bamFile, "\n"
        sys.exit()
    return(discordantFile, splitFile)
