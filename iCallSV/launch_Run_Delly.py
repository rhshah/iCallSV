'''
Created on November 19, 2015
Description: This module will be launching delly using Run_Delly
@author: Ronak H Shah
'''
'''
::Inputs::
args: Arguments passed to iCallSV
config: configuration file passed to iCallSV
sampleOutdirForDelly: Output directory for delly vcf files.
'''
import logging
import Run_Delly as rd


def launch_delly_for_different_analysis_type(args, config, sampleOutdirForDelly):
    verbose = args.verbose

    # Run Delly for Deletion
    if(verbose):
        logging.info("launch_Run_Delly: Launched Delly for Deletion Events")

    del_vcf = rd.run(
        delly=config.get("SVcaller", "DELLY"),
        analysisType="DEL",
        reference=config.get("ReferenceFasta", "REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.patientId,
        mapq=config.get("ParametersToRunDelly", "MAPQ"),
        excludeRegions=config.get("ExcludeRegion", "EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=False)

# Run Delly for duplication
    if(verbose):
        logging.info("launch_Run_Delly: Launched Delly for Duplication Events")

    dup_vcf = rd.run(
        delly=config.get("SVcaller", "DELLY"),
        analysisType="DUP",
        reference=config.get("ReferenceFasta", "REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.patientId,
        mapq=config.get("ParametersToRunDelly", "MAPQ"),
        excludeRegions=config.get("ExcludeRegion", "EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=False)

# Run Delly for inversion
    if(verbose):
        logging.info("launch_Run_Delly: Launched Delly for Inversion Events")

    inv_vcf = rd.run(
        delly=config.get("SVcaller", "DELLY"),
        analysisType="INV",
        reference=config.get("ReferenceFasta", "REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.patientId,
        mapq=config.get("ParametersToRunDelly", "MAPQ"),
        excludeRegions=config.get("ExcludeRegion", "EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=False)

# Run Delly for Translocation
    if(verbose):
        logging.info("launch_Run_Delly: Launched Delly for Translocation Envents")

    tra_vcf = rd.run(
        delly=config.get("SVcaller", "DELLY"),
        analysisType="TRA",
        reference=config.get("ReferenceFasta", "REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.patientId,
        mapq=config.get("ParametersToRunDelly", "MAPQ"),
        excludeRegions=config.get("ExcludeRegion", "EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=False)
'''
# Run Delly for Insertion
    if(verbose):
        logging.info("launch_Run_Delly: Launched Delly for Insertion Events")

    ins_vcf = rd.run(
        delly=config.get("SVcaller", "DELLY"),
        analysisType="INS",
        reference=config.get("ReferenceFasta", "REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.patientId,
        mapq=config.get("ParametersToRunDelly", "MAPQ"),
        excludeRegions=config.get("ExcludeRegion", "EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=False)
'''
    return(del_vcf, dup_vcf, inv_vcf, tra_vcf)
