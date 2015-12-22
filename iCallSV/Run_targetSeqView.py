'''
Created on Mar 19, 2015
Description: This module will run targetSeqView
@author: Ronak H Shah
'''
'''
::Input::
RLocation : Location of the R executable (>3.1.2).
targetSeqView: Location of R script that will run tragetSeqView
nodes : Number of parallel nodes for running targetSeqView
bamFile : Location of the bamFile which has the  structural variant events.
svFile : targetSeqView compatible input structural variant file.
build : Which human reference file to be used, hg18,hg19 or hg38
readLength : Sequencing Read Length (101)
outputDir : Directory for output files
outsvFile : Name of the output structural variant file that has added confidence score to it.
::Output::
Tab-delimited File with added confidence score and Image of each SV
'''
import os
import sys
from subprocess import Popen
import time
from datetime import date, timedelta
import checkparameters as cp
import logging

def run(
        RLocation,
        targetSeqView,
        nodes,
        bamFilePath,
        svFile,
        build,
        readLength,
        outputDir,
        outsvFileName):
    start_time = time.time()
    logging.info("We will now be running targetSeqView. Hope fully the R package targetSeqView is installed.")
    cp.checkFile(targetSeqView)
    cp.checkInt(nodes, "Number of nodes to run targetSeqView")
    cp.checkDir(bamFilePath)
    cp.checkFile(svFile)
    cp.checkEmpty(build, "Genome build to be used for targetSeqView")
    cp.checkInt(readLength, "Sequencing Read Length")
    cp.checkDir(outputDir)
    logging.info("All Input Parameters look good. Lets Run targetSeqView")
    RLocation = RLocation + "/bin/R"
    myPid = os.getpid()
    day = date.today()
    today = day.isoformat()
    logging.info("Run_targetSeqView: ProcessID: %s, Date: %s", myPid, today )
    outputFile = outputDir + "/" + outsvFileName
    stdoutFile = outputDir + "/" + outsvFileName[:-4] + "_" + str(myPid) + ".stdout"
    stderrFile = outputDir + "/" + outsvFileName[:-4] + "_" + str(myPid) + ".stderr"
    cmd = RLocation + " --slave --vanilla --args " + str(nodes) + " " + bamFilePath + " " + svFile + " " + build + " " + str(
        readLength) + " " + outputDir + " " + outsvFileName + " < " + targetSeqView + " > " + stdoutFile + " 2> " + stderrFile
    logging.info("Run_targetSeqView: Command that will be run %s", cmd)
    # Remove if the file exists
    if(os.path.isfile(outputFile)):
        os.remove(outputFile)
    proc = Popen(cmd, shell=True)
    proc.wait()
    retcode = proc.returncode
    if(retcode >= 0):
        end_time = time.time()
        totaltime = str(timedelta(seconds=end_time - start_time))
        logging.info("Run_targetSeqView: We have finished running targetSeqView for %s using local machine.", svFile)
        logging.info("Run_targetSeqView Duration: %s", totaltime )
    else:
        logging.info("Run_targetSeqView: targetSeqView is either still running on local machine or it errored out with return code %d for %s", retcode, svFile)
        sys.exit()
    return(outputFile)
# # Test module
# run("/home/shahr2/.Renv/versions/3.1.2/bin/R",
#     "/home/shahr2/workspace/iCallSV/iCallSV/R/Rscripts/calculateConfidenceScore.R", 5,
#     "/home/shahr2/.Renv/versions/3.1.2/lib64/R/library/targetSeqView/extdata",
#     "/home/shahr2/.Renv/versions/3.1.2/lib64/R/library/targetSeqView/extdata/targetCaptureSVs.txt",
#     "hg19", 100, "/dmp/hot/shahr2/IMPACT/Test/SVtest/", "testTargetSeqView.txt")
