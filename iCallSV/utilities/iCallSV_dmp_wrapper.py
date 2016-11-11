"""
iCallSV_dmp_wrapper
~~~~~~~~~~~~~~~~~~~

:Description: iCallSV is a wrapper to run the iCallSV package on MSKCC data
:author:     Ronak H Shah
:copyright:  (c) 2015-2016 by Ronak H Shah for Memorial Sloan Kettering Cancer Center. All rights reserved.
:license:    Apache License 2.0
:contact:    rons.shah@gmail.com
:deffield    updated: Updated

"""


import argparse
import os
import sys
import time
import glob
import pandas as pd
from subprocess import Popen, PIPE
import shlex
import re
# from threading import Thread, Lock
from gridmap import Job, process_jobs


def main():
    
    parser = argparse.ArgumentParser(
        prog='iCallSV_dmp_wrapper.py',
        description='Run iCallSV on selected pools using MSK data',
        usage='%(prog)s [options]')
    parser.add_argument(
        "-fl",
        "--folderList",
        action="store",
        dest="folderList",
        required=True,
        metavar='folders.fof',
        help="Full path folders file of files.")
    parser.add_argument(
        "-qc",
        "--qcLocation",
        action="store",
        dest="qcLocation",
        required=True,
        metavar='/some/path/qcLocation',
        help="Full path qc files.")
    parser.add_argument(
        "-b",
        "--bamLocation",
        action="store",
        dest="bamLocation",
        required=True,
        metavar='/some/path/bamlocation',
        help="Full path bam files.")
    parser.add_argument(
        "-P",
        "--python",
        action="store",
        dest="python",
        required=False,
        default='/dmp/resources/prod/tools/system/python/production/bin/python',
        metavar='/somepath/python',
        help="Full path Pyhton executables.")
    parser.add_argument(
        "-icsv",
        "--iCallSV",
        action="store",
        dest="icsv",
        required=False,
        default='/home/shahr2/git/iCallSV-master/iCallSV/iCallSV.py',
        metavar='/somepath/iCallSV.py',
        help="Full path iCallSV.py executables.")
    parser.add_argument(
        "-conf",
        "--iCallSVconf",
        action="store",
        dest="conf",
        required=False,
        default='/home/shahr2/git/iCallSV-master/iCallSV/configuration/template_dmp.ini',
        metavar='/somepath/template.ini',
        help="Full path configuration file to run iCallSV")
    parser.add_argument(
        "-q",
        "--queue",
        action="store",
        dest="queue",
        required=False,
        default='test.q',
        metavar='all.q or clin.q',
        help="Name of the SGE queue")
    parser.add_argument(
        "-qsub",
        "--qsubPath",
        action="store",
        dest="qsub",
        required=False,
        default='/common/sge/bin/lx-amd64/qsub',
        metavar='/somepath/qsub',
        help="Full Path to the qsub executables of SGE.")
    parser.add_argument(
        "-t",
        "--threads",
        action="store",
        dest="threads",
        required=False,
        default='5',
        metavar='5',
        help="Number of Threads to be used to run iCallSV")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        help="make lots of noise [default]")
    parser.add_argument(
        "-o",
        "--outDir",
        action="store",
        dest="outDir",
        required=True,
        metavar='/somepath/output',
        help="Full Path to the output dir.")
    parser.add_argument(
        "-of",
        "--outFile",
        action="store",
        dest="outFile",
        required=True,
        metavar='outputfile.txt',
        help="Name of the final output file.")

    args = parser.parse_args()
    if(args.verbose):
        print "Starting the Process to Run iCallSV."
    with open(args.folderList, 'r') as filecontent:
        for line in filecontent:
            poolName = line.rstrip()
            # print poolName
            if poolName:
                if(args.verbose):
                    print "####Processing data for:", poolName
                (titleFile,
                 outdir,
                 HSmetricsFileList,
                 bamFileList) = SetupRun(poolName, args)
                RunPerPool(
                    titleFile,
                    outdir,
                    HSmetricsFileList,
                    bamFileList,
                    args)
    interesting_files = glob.glob(
        os.path.join(
            args.outDir,
            "IMPACT*",
            "StructuralVariantAnalysis","DellyDir","*","*final.txt"))
    outfilename = os.path.join(args.outDir, args.outFile)
    header_saved = False
    with open(outfilename,'wb') as fout:
        for filename in interesting_files:
            with open(filename) as fin:
                header = next(fin)
                if not header_saved:
                    fout.write(header)
                    header_saved = True
                for line in fin:
                    fout.write(line)
    
    if(args.verbose):
        print "Finished the Process to Run iCallSV."


def SetupRun(poolName, args):
    """This will setup the run to be analyzed.


    :param str poolName: str of pool to be analyzed
    :param Namespace args: Namespace of args to get other variables
    :return: Multiple objects
    :rtype: list

    """
    bamlocation = args.bamLocation + "/" + poolName
    # Get the qc data location
    baseqclocation = args.qcLocation + "/" + poolName + "/*"
    # print bamlocation,baseqclocation
    all_qc_subdirs = getSubDirs(baseqclocation)
    # print all_qc_subdirs
    qclocation = max(all_qc_subdirs, key=os.path.getmtime)
    if(os.path.isdir(bamlocation)):
        if(args.verbose):
            print "\tBAM Location:", bamlocation, "\n"
    else:
        if(args.verbose):
            print "\tBAM LOCATION", bamlocation, " DOES NOT EXISTS!!!, Please Review you bamLocation INPUT\n"
        sys.exit()
    if(os.path.isdir(qclocation)):
        if(args.verbose):
            print "\tQC Location:", qclocation, "\n"
    else:
        if(args.verbose):
            print "\tQC LOCATION", qclocation, " DOES NOT EXISTS!!!, Please Review you qcLocation INPUT\n"
    titleFile = poolName.replace("/", "")
    titleFile = qclocation + "/" + titleFile + "_title.txt"
    if(os.path.isfile(titleFile)):
        if(args.verbose):
            print "\tTitleFile:", titleFile, "\n"
    else:
        if(args.verbose):
            print "\tTitleFile:", titleFile, " DOSE NOT EXISITS!!!, Please make sure that all qcLocation have title file\n"

    HSmetricsFile = qclocation + "/*.HSmetrics.txt"
    HSmetricsFileList = glob.glob(HSmetricsFile)
    if(args.verbose):
        print "\tAllHSmetricsFiles:\n"
    for file in HSmetricsFileList:
        if(os.path.isfile(file)):
            if(args.verbose):
                print"\t\t", file, "\n"
        else:
            if(args.verbose):
                print "\t\t", file, "Does Not Exist!!!, Please make sure qcLocation have all HSmetrics.txt files\n"

    bamFile = bamlocation + "/*.bam"
    bamFileList = glob.glob(bamFile)
    if(args.verbose):
        print "\tAllBamFiles:\n"
    for file in bamFileList:
        if(os.path.isfile(file)):
            if(args.verbose):
                print"\t\t", file, "\n"
        else:
            if(args.verbose):
                print "\t\t", file, "Does Not Exist!!!, Please make sure bamLocation have all bam files\n"
    outdir = poolName.replace("/", "")
    outdir = args.outDir + "/" + outdir
    return(titleFile, outdir, HSmetricsFileList, bamFileList)


def RunPerPool(titleFile, outdir, HSmetricsFileList, bamFileList, args):
    """This will run the pool to be analyzed.


    :param str titleFile: str of meta information file
    :param str outdir: str of output directory
    :param list HSmetricsFileList: list of picard hsmetrics files
    :param list bamFileList: list of bam files
    :param Namespace args: Namespace of args to get other variables
    :return: None
    :rtype: None

    """
    # Run Preprocess
    titleFileDF = pd.read_csv(titleFile, sep='\t', header=0, keep_default_na='True')
    groupByPatientId = titleFileDF.groupby('Patient_ID')
    baseNames = {}
    jobs = []
    poolidRegXcompile = re.compile('.*[PoolNormal|PooledNormal].*')
    poolHsmetricsFile = filter(poolidRegXcompile.match, HSmetricsFileList).pop()
    poolbamFile = filter(poolidRegXcompile.match, bamFileList).pop()
    for patientID, group in groupByPatientId:
        print patientID, ":"
        tsampleId = ''
        tBamFile = ''
        nBamFile = ''
        basename = ''
        nsampleId = ''
        if(os.path.isdir(outdir)):
            if(args.verbose):
                print "Pool Output Dir:", outdir, "exists!!!"
        else:
            os.mkdir(outdir)
            os.chmod(outdir, 0o755)
        for count, row in group.iterrows():
            sampleId = row.loc['Sample_ID']
            patientId = row.loc['Patient_ID']
            sampleClass = row.loc['Class']
            idRegXcompile = re.compile('.*' + sampleId + '.*')
            if(sampleClass == "Tumor"):
                basename = sampleId
                tBamFile = filter(idRegXcompile.match, bamFileList).pop()
                os.symlink(tBamFile, os.path.join(outdir, os.path.basename(tBamFile)))
                tBamFile = os.path.join(outdir, os.path.basename(tBamFile))
                tsampleId = sampleId
            if(sampleClass == "Normal"):
                nBamFile = filter(idRegXcompile.match, bamFileList).pop()
                os.symlink(nBamFile, os.path.join(outdir, os.path.basename(nBamFile)))
                nBamFile = os.path.join(outdir, os.path.basename(nBamFile))
                nsampleId = sampleId
                nHSmetricsFile = filter(idRegXcompile.match, HSmetricsFileList).pop()
                (decision) = SelectNormal(nHSmetricsFile, poolHsmetricsFile)
                if(decision == 'UnMatched'):
                    nBamFile = poolbamFile
                else:
                    if(args.verbose):
                        print "Matched Sample\n"
        if(os.path.isfile(tBamFile) and (os.path.isfile(nBamFile))):
            jobId = "iCallSV_" + str(count) + "_" + str(basename)
            cmdList = []
            cmd = args.python + " " + args.icsv + " -sc " + args.conf + " -bbam " + nBamFile + " -abam " + \
                tBamFile + " -aId " + tsampleId + " -bId " + nsampleId + " -op " + tsampleId + " -o " + outdir + " -v"
            # cmd = str(cmd)
            threads = int(args.threads)
            threads = threads + 1
            qsub_cmd = args.qsub + " -q " + args.queue + " -N " + jobId + " -o " + jobId + ".stdout" + " -e " + jobId + ".stderr" + \
                " -V -l h_vmem=6G,virtual_free=6G -pe smp " + str(threads) + " -wd " + outdir + " -sync y " + " -b y " + cmd
            print "qsub_cmd:", qsub_cmd, "\n"
            cmdList.append(qsub_cmd)
            job = Job(
                RunJob,
                cmdList,
                kwlist=None,
                cleanup=True,
                mem_free="2G",
                name=jobId,
                num_slots=1,
                queue=args.queue)
            jobs.append(job)
    print("sending function jobs to cluster")
    print("")

    job_outputs = process_jobs(
        jobs,
        max_processes=10,
        temp_dir='/dmp/analysis/SCRATCH/',
        white_list=None,
        quiet=False,
        local=False)

    print("results from each job")
    for (i, result) in enumerate(job_outputs):
        print("Job {0}- result: {1}".format(i, result))

    return


def RunJob(cmd):
    """Given a command run the job.


    :param str cmd: str of command to be run on the local machine
    :return: None
    :rtype: None

    """
    args = shlex.split(cmd)
    proc = Popen(args)
    proc.wait()
    retcode = proc.returncode
    if(retcode >= 0):
        print "I have finished running process using SGE"

'''
    #iterate over jobs and put each into the queue in sequence
    for job in jobs:
        print "inserting job into the queue: %s"%(job)
        jobqueue.put(job)
    #start some threads, each one will process one job from the queue#
    for i in range(mp.cpu_count()-1):
        th = Thread(target=processor,args = (i,jobqueue))
        th.setDaemon(True)
        th.start()
    #wait until all jobs are processed before quitting
    jobqueue.join()
    with jobqueue.mutex:
        jobqueue.queue.clear()
'''


def SelectNormal(normal, poolnormal):
    """Select the best possible normal.

    :param str normal: str of match normal
    :param str poolnormal: str of pool normal
    :return: str with decision whether to run matched or unmatched
    :rtype: str

    """
    filesToProcess = [normal, poolnormal]
    fcount = 0
    coveragePerFile = []
    covgIndex = 0
    decision = ''
    for file in filesToProcess:

        with open(file, 'r') as filecontent:
            coverageForSample = 0
            linecount = 0
            for line in filecontent:
                line = line.rstrip()
                if line:
                    if line.startswith("#"):
                        linecount = linecount + 1
                        continue
                    # print linecount,line
                    data = line.rstrip('\n').split('\t')
                    if(fcount == 0):
                        if(data[0] == "BAIT_SET"):
                            for content in data:
                                # print linecount,content ,"C\n"
                                if(content == 'MEAN_TARGET_COVERAGE'):
                                    covgIndex = covgIndex
                                    # print covgIndex,"CI\n"
                                    break
                                else:
                                    covgIndex = covgIndex + 1
                        else:
                            # print fcount,covgIndex,data,"Inside\n"
                            coverageForSample = data[covgIndex]
                    if(fcount > 0):
                        if(data[0] == "BAIT_SET"):
                            continue
                        # print fcount,covgIndex,data,"All\n"
                        coverageForSample = data[covgIndex]
                    linecount = linecount + 1
        fcount = fcount + 1
        coveragePerFile.append(coverageForSample)
    mnormalCovg = coveragePerFile[0]
    umnormalCovg = coveragePerFile[1]
    if(float(mnormalCovg) > 50.0):
        decision = 'Matched'
    else:
        decision = 'UnMatched'
    return(decision)


def getSubDirs(dirLocation):
    """
    Get all sub directories.


    :param str dirLocation: str of directory location
    :return: list of all sub directories
    :rtype: list

    """
    dirs = []
    for d in glob.glob(dirLocation):
        if os.path.isdir(d):
            dirs.append(d)
    return(dirs)


def processor(i, jobqueue):
    """Operate on a jobqueue.


    :param int i: count of the job
    :param Namespace jobqueue: Namespace for jobqueue
    :return: None
    :rtype: None

    """
    devnull = open('/dev/null', 'w')
    if jobqueue.empty():
        print "the Queue is empty!"
        sys.exit(1)
    try:
        job = jobqueue.get()
        print "I'm operating on job item: %s\n" % i, job
        qsub_args = shlex.split(job)
        p = Popen(qsub_args, stdout=PIPE, stderr=devnull)
        jobqueue.task_done()
    except:
        print "Failed to operate on job\n"

if __name__ == "__main__":
    start_time = time.time()
    # mp.freeze_support()
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))
