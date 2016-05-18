"""
Created on 10/28/2014.

@Ronak Shah

"""
import argparse
import os
import sys
import time
import stat
import glob
import pandas as pd
from subprocess import Popen, PIPE
import shlex
import shutil
from pprint import pprint
import re
# import multiprocessing as mp
from Queue import *
# from threading import Thread, Lock
from gridmap import Job, process_jobs


def main():
    parser = argparse.ArgumentParser(
        prog='Run_PyLoh.py',
        description='Run_PyLOH on selected pools using MSK data',
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
        "-pl",
        "--pyloh",
        action="store",
        dest="pyloh",
        required=False,
        default='/home/shahr2/.local/bin/PyLOH.py',
        metavar='/somepath/PyLOH.py',
        help="Full path PyLOH.py executables.")
    parser.add_argument(
        "-r",
        "--referenceFile",
        action="store",
        dest="ref",
        required=False,
        default='/dmp/data/pubdata/hg-fasta/production/Homo_sapiens_assembly19.fasta',
        metavar='/somepath/Homo_Sapeins_hg19.fasta',
        help="Full Path to the reference file with the bwa index.")
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
        "-md",
        "--minimum_depth",
        action="store",
        dest="minDepth",
        required=False,
        default='20',
        metavar='20',
        help="Minimum depth in both normal and tumor sample required to use a site in the analysis.")
    parser.add_argument(
        "-mbq",
        "--minimum_base_qual",
        action="store",
        dest="minBQ",
        required=False,
        default='20',
        metavar='20',
        help="Minimum base quality required for each base.")
    parser.add_argument(
        "-mmq",
        "--minimum_mapping_qual",
        action="store",
        dest="minMQ",
        required=False,
        default='0',
        metavar='0',
        help="Minimum mapping quality required for each base.")
    parser.add_argument(
        "-t",
        "--threads",
        action="store",
        dest="threads",
        required=False,
        default='5',
        metavar='5',
        help="Number of Threads to be used to run Pindel")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=True,
        help="make lots of noise [default]")
    parser.add_argument(
        "-o",
        "--outDir",
        action="store",
        dest="outDir",
        required=True,
        metavar='/somepath/output',
        help="Full Path to the output dir.")

    args = parser.parse_args()
    jobqueue = Queue()

    if(args.verbose):
        print "Starting the Process to Run PyLOH."
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
                     bamFileList,
                     segmentFileList) = SetupRun(poolName,
                                                 args)
                    RunPerPool(
                        titleFile,
                        outdir,
                        HSmetricsFileList,
                        bamFileList,
                        segmentFileList,
                        args,
                        jobqueue)

    if(args.verbose):
        print "Finished the Process to Run PyLOH."


def SetupRun(poolName, args):

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

    segFile = qclocation + "/*_copynumber.seg"
    segFileList = glob.glob(segFile)
    if(args.verbose):
        print "\tAllSegFiles:\n"
    for file in segFileList:
        if(os.path.isfile(file)):
            if(args.verbose):
                print"\t\t", file, "\n"
        else:
            if(args.verbose):
                print "\t\t", file, "Does Not Exist!!!, Please make sure qcLocation have all segments files\n"

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
    return(titleFile, outdir, HSmetricsFileList, bamFileList, segFileList)


def RunPerPool(titleFile, outdir, HSmetricsFileList, bamFileList, segmentFileList, args, jobqueue):
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
        outTargetFile = ''
        tBamFile = ''
        nBamFile = ''
        basename = ''
        toutdir = ''
        if(os.path.isdir(outdir)):
            if(args.verbose):
                print "Pool Output Dir:", outdir, "exists!!!"
        else:
            os.mkdir(outdir)
            os.chmod(outdir, 0o755)
        for count, row in group.iterrows():
            bcId = row.loc['Barcode']
            poolId = row.loc['Pool']
            sampleId = row.loc['Sample_ID']
            patientId = row.loc['Patient_ID']
            sampleClass = row.loc['Class']
            idRegXcompile = re.compile('.*' + sampleId + '.*')
            if(sampleClass == "Tumor"):
                toutdir = outdir + "/" + sampleId
                if(os.path.isdir(toutdir)):
                    if(args.verbose):
                        print "Output Dir:", toutdir, "exists!!!"
                else:
                    os.mkdir(toutdir)
                    os.chmod(toutdir, 0o755)
                outTargetFile = toutdir + "/" + sampleId + "_targetRegion.bed"
                txt_fh = open(outTargetFile, "wb")
                txt_fh.write("chrom\tloc.start\tloc.end\n")
                basename = sampleId
                tBamFile = filter(idRegXcompile.match, bamFileList).pop()
                segfile = filter(idRegXcompile.match, segmentFileList).pop()
                if(segfile):
                    segFileDF = pd.read_csv(segfile, sep=' ', header=0, keep_default_na='True')
                    for segcount, segrow in segFileDF.iterrows():
                        chr = segrow.loc['chrom']
                        start = segrow.loc['loc.start']
                        end = segrow.loc['loc.end']
                        txt_fh.write(str(chr) + "\t" + str(start) + "\t" + str(end) + "\n")
                    txt_fh.close()
            if(sampleClass == "Normal"):
                nBamFile = filter(idRegXcompile.match, bamFileList).pop()
                nHSmetricsFile = filter(idRegXcompile.match, HSmetricsFileList).pop()
                (decision) = SelectNormal(nHSmetricsFile, poolHsmetricsFile)
                if(decision == 'UnMatched'):
                    nBamFile = poolbamFile
                else:
                    if(args.verbose):
                        print "Matched Sample\n"
        if(os.path.isfile(tBamFile) and (os.path.isfile(nBamFile))and (os.path.isfile(outTargetFile))):
            jobId_preprocess = "Preprocess_" + str(count) + "_" + str(basename)
            baseNames[basename] = toutdir + "#" + jobId_preprocess
            outsegFile = toutdir + "/" + basename + '.PyLOH.segments'
            if(os.path.isfile(outsegFile)):
                continue
            else:
                cmdList = []
                cmd = args.python + " " + args.pyloh + " preprocess " + args.ref + " " + nBamFile + " " + tBamFile + " " + basename + " --segments_bed " + outTargetFile + \
                    " --min_depth " + args.minDepth + " --min_base_qual " + args.minBQ + " --min_map_qual " + args.minMQ + " --process_num " + args.threads + " --WES"
                # cmd = str(cmd)
                threads = int(args.threads)
                threads = threads + 1
                qsub_cmd = args.qsub + " -q " + args.queue + " -N " + jobId_preprocess + " -o " + jobId_preprocess + ".stdout" + " -e " + \
                    jobId_preprocess + ".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp " + str(threads) + " -wd " + toutdir + " -sync y " + " -b y " + cmd
                print "qsub_cmd:", qsub_cmd, "\n"
                cmdList.append(qsub_cmd)
                job = Job(
                    RunJob,
                    cmdList,
                    kwlist=None,
                    cleanup=True,
                    mem_free="2G",
                    name=jobId_preprocess,
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

    # RunModel
    count = 0
    jobs = []
    for basename, jdata in baseNames.iteritems():
        (toutdir, jobId_preprocess) = jdata.split('#', 1)
        jobId_runmodel = "RunModel_" + str(count) + "_" + str(basename)
        outsegFile = toutdir + "/" + basename + '.PyLOH.segments.extended'
        if(os.path.isfile(outsegFile)):
            continue
        else:
            cmdList = []
            # + " --allelenumber_max 2 --max_iters 100 --stop_value 1e-7"
            cmd = args.python + " " + args.pyloh + " run_model " + basename
            qsub_cmd = args.qsub + " -q " + args.queue + " -N " + jobId_runmodel + " -o " + jobId_runmodel + ".stdout" + " -e " + \
                jobId_runmodel + ".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp 1" + " -wd " + toutdir + " -sync y " + "-b y " + cmd
            print "qsub_cmd:", qsub_cmd, "\n"
            cmdList.append(qsub_cmd)
            job = Job(
                RunJob,
                cmdList,
                kwlist=None,
                cleanup=True,
                mem_free="2G",
                name=jobId_runmodel,
                num_slots=1,
                queue=args.queue)
            jobs.append(job)
            count = count + 1

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

    # Run BAF
    count = 0
    jobs = []
    for basename, jdata in baseNames.iteritems():
        (toutdir, jobId_preprocess) = jdata.split('#', 1)
        jobId_BAF = "BAF_" + str(count) + "_" + str(basename)
        cmdList = []
        cmd = args.python + " " + args.pyloh + " BAF_heatmap " + basename
        qsub_cmd = args.qsub + " -q " + args.queue + " -N " + jobId_BAF + " -o " + jobId_BAF + ".stdout" + " -e " + \
            jobId_BAF + ".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp 1" + " -wd " + toutdir + " -sync y " + " -b y " + cmd
        print "qsub_cmd:", qsub_cmd, "\n"
        cmdList.append(qsub_cmd)
        job = Job(
            RunJob,
            cmdList,
            kwlist=None,
            cleanup=True,
            mem_free="2G",
            name=jobId_BAF,
            num_slots=1,
            queue=args.queue)
        jobs.append(job)
        count = count + 1

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
    dirs = []
    for d in glob.glob(dirLocation):
        if os.path.isdir(d):
            dirs.append(d)
    return(dirs)


def processor(i, jobqueue):
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
