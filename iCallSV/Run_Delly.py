"""
Created on March 17, 2015
Description: Runs the delly program on case and control bam file to give its results
@author: Ronak H Shah
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
"""

import os
import sys
from subprocess import Popen
import shlex
import time
from datetime import date, timedelta
import checkparameters as cp
import logging
import makebamindex as mbi
from distutils.version import LooseVersion, StrictVersion
# This function will run delly based on given inputs


def run(
        delly,
        version,
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
        debug,
        loggeroutput):
    """
    Created on March 17, 2015
    Description: Runs the delly program on case and control bam file to give its results
    @author: Ronak H Shah
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
    
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(loggeroutput)
    handler.setLevel(logging.INFO)
    formatter='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    start_time = time.time()
    if(verbose):
        logger.info("Run_Delly: We are now going to run Delly for you. It going to be exciting time.")
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
    cp.checkEmpty(version,"Delly Version")
    cp.checkFile(bcftools)
    cp.checkFile(reference)
    cp.checkFile(excludeRegions)
    cp.checkDir(outputdir)
    cp.checkInt(mapq, "Delly MAPQ")
    cp.checkEmpty(caseId, "Delly Case BAM ID")
    cp.checkDellyAnalysisType(analysisType)
    if(verbose):
        logger.info("Run_Delly: All the input parameters look good for running delly")
        logger.info("Run_Delly: ProcessID:%s,Date:%s", myPid, today)
    if(debug):
        if(version >= StrictVersion('0.7.3')):
            cmd = delly + " -t " + analysisType + " -g " + reference + " -x " + excludeRegions + \
                " -q " + str(mapq) + " -o " + outputBcf + " " + caseBam + " " + controlBam
            logger.debug("Run_Delly: Command that will be run %s", cmd)
        else:
            cmd = delly + " -t " + analysisType + " -g " + reference + " -x " + excludeRegions + \
                " -q " + str(mapq) + " -o " + outputVcf + " " + caseBam + " " + controlBam
            logger.debug("Run_Delly: Command that will be run %s", cmd)
    else:
        # Check if bam index files are there else make them
        controlBai = controlBam + ".bai"
        if(os.path.isfile(controlBai)):
            if(verbose):
                logger.info("Run_Delly: Bam Index file is present for %s ", controlBai)
        else:
            if(verbose):
                logger.warn(
                    "Run_Delly: Bam Index file is not present and we will make it for %s ",
                    controlBai)
            mbi.MakeIndex(controlBam)
        caseBai = caseBam + ".bai"
        if(os.path.isfile(caseBai)):
            if(verbose):
                logger.info("Run_Delly: Bam Index file is present for %s ", caseBai)
        else:
            if(verbose):
                logger.warn(
                    "Run_Delly: Bam Index file is not present and we will make it for %s ",
                    caseBai)
            mbi.MakeIndex(caseBam)
        if(version >= StrictVersion('0.7.3')):
            cmd = delly + " call -t " + analysisType + " -g " + reference + " -x " + excludeRegions + \
                " -q " + str(mapq) + " -o " + outputBcf + " " + caseBam + " " + controlBam
        else:
            cmd = delly + " -t " + analysisType + " -g " + reference + " -x " + excludeRegions + \
                " -q " + str(mapq) + " -o " + outputVcf + " " + caseBam + " " + controlBam
        if(verbose):
            logger.info("Run_Delly: Command that will be run:%s", cmd)
        args = shlex.split(cmd)
        proc = Popen(args)
        proc.wait()
        retcode = proc.returncode
        if(retcode >= 0):
            end_time = time.time()
            totaltime = str(timedelta(seconds=end_time - start_time))
            if(verbose):
                logger.info(
                    "Run_Delly: We have finished running Delly for %s using local machine", caseId)
                logger.info("Run_Delly Duration: %s", totaltime)
            if(version >= StrictVersion('0.7.3')):
                if(os.path.isfile(outputBcf)):
                    cmd = bcftools + " view " + outputBcf + " -O v -o " + outputVcf
                    if(verbose):
                        logger.info("Run_Delly_bcf2vcf: Command that will be run:%s", cmd)
                    args = shlex.split(cmd)
                    proc = Popen(args)
                    proc.wait()
                    retcode = proc.returncode
                    if(retcode >= 0):
                        end_time = time.time()
                        totaltime = str(timedelta(seconds=end_time - start_time))
                        if(verbose):
                            logger.info(
                                "Run_Delly_bcf2vcf: We have finished running bcftools for %s using local machine", caseId)
                            logger.info("Run_Delly_bcf2vcf Duration: %s", totaltime)
                    else:
                        if(verbose):
                            logger.fatal(
                                "Run_Delly_bcf2vcf: bcftools is either still running on local machine or it errored out with return code %d for %s",
                                retcode,
                                caseId)
                        sys.exit(1)
                else:
                    if(verbose):
                        logger.fatal(
                                      "Run_Delly_bcf2vcf: bcftools is either still running on local machine or it errored out with return code %d for %s",
                                      retcode,
                                      caseId)
                        sys.exit(1)
            else:
                if(os.path.isfile(outputVcf)):
                    return(outputVcf)
                else:
                    if(verbose):
                        logger.fatal(
                                      "Run_Delly: Delly is either still running on local machine or it errored out with return code %d for %s",
                                      retcode,
                                      caseId)
                    sys.exit(1)
        else:
            if(verbose):
                logger.fatal(
                    "Run_Delly: Delly is either still running on local machine or it errored out with return code %d for %s",
                    retcode,
                    caseId)
            sys.exit(1)
    return(outputVcf)


# Testing the module
# run("/dmp/resources/dev2/bin/delly","/dmp/resources/dev2/bin/bcftools", "TRA", "/dmp/data/pubdata/hg-fasta/production/Homo_sapiens_assembly19.fasta", "/dmp/hot/shahr2/IMPACT/Test/SVtest/35462375-N_bc45_IMPACTv5-CLIN-20150050_L000_mrg_cl_aln_srt_MD_IR_BR.bam", "/dmp/hot/shahr2/IMPACT/Test/SVtest/35462375-T_bc44_IMPACTv5-CLIN-20150050_L000_mrg_cl_aln_srt_MD_IR_BR.bam","35462375-T", 20, "/dmp/data/mskdata/sv-files/production/human.hg19.excl.tsv", "/dmp/hot/shahr2/IMPACT/Test/SVtest/", False)
