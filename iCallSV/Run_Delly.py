'''
Created on March 17, 2015
Description: Runs the delly program on case and control bam file to give its results
@author: Ronak H Shah
'''
'''
::Inputs::
delly:Path to delly executables (0.7.3 or above)
bcftools:Path to bcftools executables (1.3.1 or above)
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
import logging
# This function will run delly based on given inputs


def run(
        delly,
        bcftools,
        analysisType,
        reference,
        controlBam,
        caseBam,
        caseId,
        mapq,
        excludeRegions,
        outputdir,
        verbose,
        debug):
    start_time = time.time()
    if(verbose):
        logging.info("Run_Delly: We are now going to run Delly for you. It going to be exciting time.")
    myPid = os.getpid()
    day = date.today()
    today = day.isoformat()
    tag = analysisType.lower()
    outputBcf = outputdir + "/" + caseId + "_" + tag + ".bcf"
    outputVcf = outputdir + "/" + caseId + "_" + tag + ".vcf"
    # Check all input variables
    cp.checkFile(controlBam)
    cp.checkFile(caseBam)
    cp.checkFile(delly)
    cp.checkFile(bcftools)
    cp.checkFile(reference)
    cp.checkFile(excludeRegions)
    cp.checkDir(outputdir)
    cp.checkInt(mapq, "Delly MAPQ")
    cp.checkEmpty(caseId, "Delly Case BAM ID")
    cp.checkDellyAnalysisType(analysisType)
    if(verbose):
        logging.info("Run_Delly: All the input parameters look good for running delly")
        logging.info("Run_Delly: ProcessID:%s,Date:%s", myPid, today)
    if(debug):
        cmd = delly + " -t " + analysisType + " -g " + reference + " -x " + excludeRegions + \
            " -q " + str(mapq) + " -o " + outputBcf + " " + caseBam + " " + controlBam
        logging.debug("Run_Delly: Command that will be run %s", cmd)
    else:
        # Check if bam index files are there else make them
        controlBai = controlBam + ".bai"
        if(os.path.isfile(controlBai)):
            if(verbose):
                logging.info("Run_Delly: Bam Index file is present for %s ", controlBai)
        else:
            if(verbose):
                logging.warn(
                    "Run_Delly: Bam Index file is not present and we will make it for %s ",
                    controlBai)
            mbi.MakeIndex(controlBam)
        caseBai = caseBam + ".bai"
        if(os.path.isfile(caseBai)):
            if(verbose):
                logging.info("Run_Delly: Bam Index file is present for %s ", caseBai)
        else:
            if(verbose):
                logging.warn(
                    "Run_Delly: Bam Index file is not present and we will make it for %s ",
                    caseBai)
            mbi.MakeIndex(caseBam)
        cmd = delly + "call -t " + analysisType + " -g " + reference + " -x " + excludeRegions + \
            " -q " + str(mapq) + " -o " + outputBcf + " " + caseBam + " " + controlBam
        if(verbose):
            logging.info("Run_Delly: Command that will be run:%s", cmd)
        args = shlex.split(cmd)
        proc = Popen(args)
        proc.wait()
        retcode = proc.returncode
        if(retcode >= 0):
            end_time = time.time()
            totaltime = str(timedelta(seconds=end_time - start_time))
            if(verbose):
                logging.info("Run_Delly: We have finished running Delly for %s using local machine", caseId)
                logging.info("Run_Delly Duration: %s", totaltime)
            if(os.path.isfile(outputBcf)):
                retcode = 1
            else:
                cmd = bcftools  + " view "  + outputBcf  + " > " + outputVcf
                if(verbose):
                    logging.info("Run_Delly_bcf2vcf: Command that will be run:%s", cmd)
                args = shlex.split(cmd)
                proc = Popen(args)
                proc.wait()
                retcode = proc.returncode
                if(retcode >= 0):  
                    end_time = time.time()
                    totaltime = str(timedelta(seconds=end_time - start_time))
                    if(verbose):
                        logging.info(
                            "Run_Delly_bcf2vcf: We have finished running bcftools for %s using local machine", caseId)
                        logging.info("Run_Delly_bcf2vcf Duration: %s", totaltime)
                else:
                    logging.fatal(
                    "Run_Delly_bcf2vcf: bcftools is either still running on local machine or it errored out with return code %d for %s",
                    retcode,
                    caseId)
                    sys.exit(1)
        else:
            logging.fatal(
                "Run_Delly: Delly is either still running on local machine or it errored out with return code %d for %s",
                retcode,
                caseId)
            sys.exit(1)
    return(outputVcf)


# Testing the module
# run("/dmp/resources/dev2/bin/delly","/dmp/resources/dev2/bin/bcftools", "TRA", "/dmp/data/pubdata/hg-fasta/production/Homo_sapiens_assembly19.fasta", "/dmp/hot/shahr2/IMPACT/Test/SVtest/35462375-N_bc45_IMPACTv5-CLIN-20150050_L000_mrg_cl_aln_srt_MD_IR_BR.bam", "/dmp/hot/shahr2/IMPACT/Test/SVtest/35462375-T_bc44_IMPACTv5-CLIN-20150050_L000_mrg_cl_aln_srt_MD_IR_BR.bam","35462375-T", 20, "/dmp/data/mskdata/sv-files/production/human.hg19.excl.tsv", "/dmp/hot/shahr2/IMPACT/Test/SVtest/", False)
