'''
Created on March 17, 2015
Description: Runs the delly program on case and control bam file to give its results
@author: Ronak H Shah
'''
'''
::Inputs::
delly:Path to delly executables (0.6.1 or above)
type: What ot run in delly, DEL:Deletion, DUP: Duplication,TRA:Translocation, INV:Inversion
reference: Reference Genome that was used to align the reads.
controlBam: Path to control/normal bam file
caseBam: Path to case/tumor bam file
controlID: Id of the control/normal sample
caseID: Id of the case/tumor sample
mapq: mapping quality cutoff for delly
excludeRegions: Regions to be excluded for calling structural variation.
outputdir: directory for the output of delly
debug: If you just wish to test what we will do
::Output::
VCF file having all the structural variants called
'''
import os
import sys
from subprocess import Popen
import shlex
import time
from datetime import date, timedelta
import checkparameters as cp
import makebamindex as mbi

# This function will run delly based on given inputs


def run(
        delly,
        analysisType,
        reference,
        controlBam,
        caseBam,
        controlId,
        caseId,
        mapq,
        excludeRegions,
        outputdir,
        debug):
    start_time = time.time()
    print "Run_Delly: We are now going to run Delly for you. It going to be exciting time.\n"
    myPid = os.getpid()
    day = date.today()
    today = day.isoformat()
    tag = analysisType.lower()
    outputVcf = outputdir + "/" + caseId + "_" + tag + ".vcf"
    # Check all input variables
    cp.checkFile(controlBam)
    cp.checkFile(caseBam)
    cp.checkFile(delly)
    cp.checkFile(reference)
    cp.checkFile(excludeRegions)
    cp.checkDir(outputdir)
    cp.checkInt(mapq, "Delly MAPQ")
    cp.checkEmpty(controlId, "Delly Control BAM ID")
    cp.checkEmpty(caseId, "Delly Case BAM ID")
    cp.checkDellyAnalysisType(analysisType)
    print "Run_Delly: All the input parameters look good for running delly\n"
    print "Run_Delly: ProcessID:", myPid, " Date:", today, "\n"
    if(debug):
        cmd = delly + " -t " + analysisType + " -g " + reference + " -x " + excludeRegions + " -q " + str(
            mapq) + " -o " + outputVcf + " " + caseBam + " " + controlBam
        print "Run_Delly: Command that will be run\n", cmd, "\n"
    else:
        # Check if bam index files are there else make them
        controlBai = controlBam + ".bai"
        if(os.path.isfile(controlBai)):
            print "Run_Delly: Bam Index file is present for ", controlBai, " \n"
        else:
            print "Run_Delly: Bam Index file is not present and we will make it for ", controlBai, " \n"
            mbi.MakeIndex(controlBam)
        caseBai = caseBam + ".bai"
        if(os.path.isfile(caseBai)):
            print "Run_Delly: Bam Index file is present for ", caseBai, " \n"
        else:
            print "Run_Delly: Bam Index file is not present and we will make it for ", caseBai, " \n"
            mbi.MakeIndex(caseBam)
        cmd = delly + " -t " + analysisType + " -g " + reference + " -x " + excludeRegions + " -q " + str(
            mapq) + " -o " + outputVcf + " " + caseBam + " " + controlBam
        print "Run_Delly: Command that will be run\n", cmd, "\n"
        args = shlex.split(cmd)
        proc = Popen(args)
        proc.wait()
        retcode = proc.returncode
        if(retcode >= 0):
            end_time = time.time()
            print "Run_Delly: We have finished running Delly for ", caseId, " using local machine.\n"
            print "Run_Delly Duration:", str(timedelta(seconds=end_time - start_time)), "\n"
        else:
            print "Run_Delly: Delly is either still running on local machine or it errored out with return code", retcode, " for", caseId, "\n"
            sys.exit()
    return(outputVcf)

# Testing the module
# run("/dmp/resources/dev2/bin/delly", "TRA", "/dmp/data/pubdata/hg-fasta/production/Homo_sapiens_assembly19.fasta", "/dmp/hot/shahr2/IMPACT/Test/SVtest/35462375-N_bc45_IMPACTv5-CLIN-20150050_L000_mrg_cl_aln_srt_MD_IR_BR.bam", "/dmp/hot/shahr2/IMPACT/Test/SVtest/35462375-T_bc44_IMPACTv5-CLIN-20150050_L000_mrg_cl_aln_srt_MD_IR_BR.bam", "35462375-N", "35462375-T", 20, "/dmp/data/mskdata/sv-files/production/human.hg19.excl.tsv", "/dmp/hot/shahr2/IMPACT/Test/SVtest/", False)
