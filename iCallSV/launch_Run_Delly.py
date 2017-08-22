"""
launch_Run_Delly
~~~~~~~~~~~~~~~~

:Description: This module will be launching delly using Run_Delly
"""

'''
Created on November 19, 2015
Description: This module will be launching delly using Run_Delly
@author: Ronak H Shah
::Inputs::
args: Arguments passed to iCallSV
config: configuration file passed to iCallSV
sampleOutdirForDelly: Output directory for delly vcf files.
'''
import os
import logging
import Run_Delly as rd
import multiprocessing as mp
import makebamindex as mbi
import coloredlogs

logger = logging.getLogger('iCallSV.launch_Run_Delly')
coloredlogs.install(level='DEBUG')

def launch_delly_for_different_analysis_type(args, config, sampleOutdirForDelly):
    """
    This will launch delly calls in parallel.

    :param Namespace args: Namespace of args to get other variables
    :param Namespace config: configuration file passed to iCallSV
    :param str sampleOutdirForDelly: Output directory for delly vcf files.
    :return: Multiple objects
    :rtype: list

    """
    verbose = args.verbose
    pool = mp.Pool(processes=4)
    analyisisType = ["DEL", "DUP", "INV", "TRA"]
    if(verbose):
        logger.info(
            "launch_Run_Delly: Launched Delly for Deletion, Duplication, Inversion and Translocation Events")
    # check of index file before run
    controlBai = args.controlBam + ".bai"
    if(os.path.isfile(controlBai)):
        if(verbose):
            logger.info("Run_Delly: Bam Index file is present for %s ", controlBai)
    else:
        if(verbose):
            logger.warn(
                "Run_Delly: Bam Index file is not present and we will make it for %s ",
                controlBai)
        mbi.MakeIndex(args.controlBam)
    caseBai = args.caseBam + ".bai"
    if(os.path.isfile(caseBai)):
        if(verbose):
            logger.info("Run_Delly: Bam Index file is present for %s ", caseBai)
    else:
        if(verbose):
            logger.warn(
                "Run_Delly: Bam Index file is not present and we will make it for %s ",
                caseBai)
        mbi.MakeIndex(args.caseBam)
    # launch commands
    results = [pool.apply_async(rd.run, args=(
        config.get("SVcaller", "DELLY"),
        config.get("SVcaller", "DellyVersion"),
        config.get("SVcaller", "BCFTOOLS"),
        x,
        config.get("ReferenceFasta", "REFFASTA"),
        args.controlBam,
        args.caseBam,
        args.caseId,
        int(config.get("ParametersToRunDelly", "MAPQ")),
        config.get("ExcludeRegion", "EXREGIONS"),
        sampleOutdirForDelly,
        verbose,
        False)) for x in analyisisType]
    output = [p.get() for p in results]
    del_vcf, dup_vcf, inv_vcf, tra_vcf = output
    return(del_vcf, dup_vcf, inv_vcf, tra_vcf)

'''
    # Run Delly for Deletion
    if(verbose):
        logger.info("launch_Run_Delly: Launched Delly for Deletion Events")

    del_vcf = rd.run(
        delly=config.get("SVcaller", "DELLY"),
        analysisType="DEL",
        reference=config.get("ReferenceFasta", "REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.caseId,
        mapq=int(config.get("ParametersToRunDelly", "MAPQ")),
        excludeRegions=config.get("ExcludeRegion", "EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=False)

# Run Delly for duplication
    if(verbose):
        logger.info("launch_Run_Delly: Launched Delly for Duplication Events")

    dup_vcf = rd.run(
        delly=config.get("SVcaller", "DELLY"),
        analysisType="DUP",
        reference=config.get("ReferenceFasta", "REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.caseId,
        mapq=int(config.get("ParametersToRunDelly", "MAPQ")),
        excludeRegions=config.get("ExcludeRegion", "EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=False)

# Run Delly for inversion
    if(verbose):
        logger.info("launch_Run_Delly: Launched Delly for Inversion Events")

    inv_vcf = rd.run(
        delly=config.get("SVcaller", "DELLY"),
        analysisType="INV",
        reference=config.get("ReferenceFasta", "REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.caseId,
        mapq=int(config.get("ParametersToRunDelly", "MAPQ")),
        excludeRegions=config.get("ExcludeRegion", "EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=False)

# Run Delly for Translocation
    if(verbose):
        logger.info("launch_Run_Delly: Launched Delly for Translocation Envents")

    tra_vcf = rd.run(
        delly=config.get("SVcaller", "DELLY"),
        analysisType="TRA",
        reference=config.get("ReferenceFasta", "REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.caseId,
        mapq=int(config.get("ParametersToRunDelly", "MAPQ")),
        excludeRegions=config.get("ExcludeRegion", "EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=False)

    return(del_vcf, dup_vcf, inv_vcf, tra_vcf)

# Run Delly for Insertion
    if(verbose):
        logger.info("launch_Run_Delly: Launched Delly for Insertion Events")

    ins_vcf = rd.run(
        delly=config.get("SVcaller", "DELLY"),
        analysisType="INS",
        reference=config.get("ReferenceFasta", "REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.caseId,
        mapq=int(config.get("ParametersToRunDelly", "MAPQ")),
        excludeRegions=config.get("ExcludeRegion", "EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=False)
'''
