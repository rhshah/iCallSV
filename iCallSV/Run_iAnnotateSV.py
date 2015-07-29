'''
Created on Mar 18, 2015
Description: This module will run iAnnotateSV package
@author: Ronak H Shah
'''
'''
::Input::
python : Location for the python executable.
iAnnotateSV : Location of the wrapper iAnnotateSV package (iAnnotateSV.py)
build : Which human reference file to be used, hg18,hg19 or hg38
inputTabFile : Tab-Delimited Input FIle compatible with iAnnotateSV package.
outputTabFile: Name of the output Tab-Delimited file with Annotations
outputDir : Name of the output directory where the outputTabFile will be written

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


def run(
        python,
        iAnnotateSV,
        build,
        distance,
        canonicalTranscriptFile,
        inputTabFile,
        outputTabFile,
        outputDir):
    start_time = time.time()
    cp.checkDir(outputDir)
    cp.checkFile(iAnnotateSV)
    cp.checkFile(inputTabFile)
    cp.checkFile(python)
    cp.checkInt(distance, "Distance for extending the promoter region")
    cp.checkEmpty(build, "Which human reference file to be used, hg18,hg19 or hg38")
    cp.checkFile(canonicalTranscriptFile)
    print "Run_iAnnotateSV: All input parameters look good. Lets run the package.\n"
    myPid = os.getpid()
    day = date.today()
    today = day.isoformat()
    print "Run_iAnnotateSV: ProcessID:", myPid, " Date:", today, "\n"
    outputFile = outputDir + "/" + outputTabFile
    cmd = python + " " + iAnnotateSV + " -r " + build + " -i " + inputTabFile + \
        " -o " + outputFile + " -d " + str(distance) + " -c " + canonicalTranscriptFile
    args = shlex.split(cmd)
    print "Run_iAnnotateSV: Command that will be run\n", cmd, "\n"
    # Reomove if the file exists
    if(os.path.isfile(outputFile)):
        os.remove(outputFile)
    proc = Popen(args)
    proc.wait()
    retcode = proc.returncode
    if(retcode >= 0):
        end_time = time.time()
        print "Run_iAnnotateSV: We have finished running iAnnotateSV for ", inputTabFile, " using local machine.\n"
        print "Run_iAnnotateSV Duration:", str(timedelta(seconds=end_time - start_time)), "\n"
    else:
        print "Run_iAnnotateSV: iAnnotateSV is either still running on local machine or it errored out with return code", retcode, " for", inputTabFile, "\n"
        sys.exit()
    return(outputFile)
# # test moudule
# run("/dmp/resources/prod/tools/system/python/production/bin/python",
#     "/home/shahr2/workspace/iAnnotateSV/iAnnotateSV/iAnnotateSV.py", "hg19", 3000,
#     "/home/shahr2/workspace/iAnnotateSV/iAnnotateSV/data/canonicalInfo/cannonical_transcripts_cv5.txt",
#     "/dmp/hot/shahr2/IMPACT/Test/SVtest//35462375-T_bc44_jmp.stdfilter.tab",
#     "35462375-T_bc44_jmp.stdfilter.anno.tab", "/dmp/hot/shahr2/IMPACT/Test/SVtest/")
