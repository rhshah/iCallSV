"""
Created on November 20, 2015
Description: This module will filter delly results and create filtered delly vcf files
@author: Ronak H Shah
"""
"""
::Inputs::
args: Arguments passed to iCallSV
config: configuration file passed to iCallSV
sampleOutdirForDelly: Output directory for delly vcf files.
del_vcf: Path to deletion based vcf file
dup_vcf: Path to duplication based vcf file
inv_vcf: Path to inversion based vcf file
tra_vcf: Path to translocation based vcf file
ins_vcf: Path to insertion based vcf file
"""

import os
import logging
import FilterDellyCalls as fdc
import multiprocessing as mp


def launch_filterdellycalls_for_different_analysis_type(
        args, config, sampleOutdirForDelly, del_vcf, dup_vcf, inv_vcf, tra_vcf):
    logger = logging.getLogger(__name__)
    verbose = args.verbose
    fileType = [del_vcf, dup_vcf, inv_vcf, tra_vcf]
    pool = mp.Pool(processes=4)
    if(verbose):
        logger.info(
            "launch_FilterDellyCalls: Launched FilterDellyCalls for Deletion, Duplication, Inversion and Translocation Events")
    results = [pool.apply_async(fdc.run, args=(
        x,
        sampleOutdirForDelly,
        args.controlId,
        args.caseId,
        config.get("HotSpotRegions", "HotspotFile"),
        config.get("BlackListRegions", "BlackListFile"),
        int(config.get("ParametersToFilterDellyResults", "LengthOfSV")),
        int(config.get("ParametersToFilterDellyResults", "OverallMapq")),
        int(config.get("ParametersToFilterDellyResults", "OverallMapqHotspot")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingSplitReadsHotspot")),
        verbose)) for x in fileType]
    output = [p.get() for p in results]
    filter_del_vcf, filter_dup_vcf, filter_inv_vcf, filter_tra_vcf = output
    return(filter_del_vcf, filter_dup_vcf, filter_inv_vcf, filter_tra_vcf)

'''
    # Run Delly for Deletion
    if(verbose):
        logger.info("launch_Run_Delly: Launched Delly for Deletion Events")
    filter_del_vcf = os.path.splitext(os.path.basename(del_vcf))[0] + "_filtered.vcf"
    filter_dup_vcf = os.path.splitext(os.path.basename(dup_vcf))[0] + "_filtered.vcf"
    filter_del_vcf = fdc.run(
        del_vcf,
        sampleOutdirForDelly,
        args.controlId,
        args.caseId,
        config.get("HotSpotRegions", "HotspotFile"),
        config.get("BlackListRegions", "BlackListFile"),
        int(config.get("ParametersToFilterDellyResults", "LengthOfSV")),
        int(config.get("ParametersToFilterDellyResults", "OverallMapq")),
        int(config.get("ParametersToFilterDellyResults", "OverallMapqHotspot")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingSplitReadsHotspot")),
        verbose)
# Run Delly for duplication
    if(verbose):
        logger.info("launch_Run_Delly: Launched Delly for Duplication Events")
    filter_dup_vcf = os.path.splitext(os.path.basename(dup_vcf))[0] + "_filtered.vcf"
    filter_dup_vcf = fdc.run(
        dup_vcf,
        sampleOutdirForDelly,
        args.controlId,
        args.caseId,
        config.get("HotSpotRegions", "HotspotFile"),
        config.get("BlackListRegions", "BlackListFile"),
        int(config.get("ParametersToFilterDellyResults", "LengthOfSV")),
        int(config.get("ParametersToFilterDellyResults", "OverallMapq")),
        int(config.get("ParametersToFilterDellyResults", "OverallMapqHotspot")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingSplitReadsHotspot")),
        verbose)
# Run Delly for inversion
    if(verbose):
        logger.info("launch_Run_Delly: Launched Delly for Inversion Events")
    filter_inv_vcf = os.path.splitext(os.path.basename(inv_vcf))[0] + "_filtered.vcf"
    filter_inv_vcf = fdc.run(
        inv_vcf,
        sampleOutdirForDelly,
        args.controlId,
        args.caseId,
        config.get("HotSpotRegions", "HotspotFile"),
        config.get("BlackListRegions", "BlackListFile"),
        int(config.get("ParametersToFilterDellyResults", "LengthOfSV")),
        int(config.get("ParametersToFilterDellyResults", "OverallMapq")),
        int(config.get("ParametersToFilterDellyResults", "OverallMapqHotspot")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingSplitReadsHotspot")),
        verbose)
# Run Delly for Translocation
    if(verbose):
        logger.info("launch_Run_Delly: Launched Delly for Translocation Envents")
    filter_tra_vcf=os.path.splitext(os.path.basename(tra_vcf))[0] + "_filtered.vcf"
    filter_tra_vcf=fdc.run(
        tra_vcf,
        sampleOutdirForDelly,
        args.controlId,
        args.caseId,
        config.get("HotSpotRegions", "HotspotFile"),
        config.get("BlackListRegions", "BlackListFile"),
        int(config.get("ParametersToFilterDellyResults", "LengthOfSV")),
        int(config.get("ParametersToFilterDellyResults", "OverallMapq")),
        int(config.get("ParametersToFilterDellyResults", "OverallMapqHotspot")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "OverallSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "CaseSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingReads")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults", "ControlSupportingSplitReadsHotspot")),
        verbose)
    return(filter_del_vcf, filter_dup_vcf, filter_inv_vcf, filter_tra_vcf)
# Run Delly for Insertion
    if(verbose):
        logger.info("launch_Run_Delly: Launched Delly for Insertion Events")
    filter_ins_vcf = os.path.splitext(os.path.basename(tra_vcf))[0] + "_filtered.vcf"
    ins_vcf = fdc.run(
        ins_vcf,
        sampleOutdirForDelly,
        args.controlId,
        args.caseId,
        config.get("HotSpotRegions","HotspotFile"),
        config.get("BlackListRegions","BlackListFile"),
        int(config.get("ParametersToFilterDellyResults","LengthOfSV")),
        int(config.get("ParametersToFilterDellyResults","OverallMapq")),
        int(config.get("ParametersToFilterDellyResults","OverallMapqHotspot")),
        int(config.get("ParametersToFilterDellyResults","OverallSupportingReads")),
        int(config.get("ParametersToFilterDellyResults","OverallSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults","OverallSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults","OverallSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults","CaseSupportingReads")),
        int(config.get("ParametersToFilterDellyResults","CaseSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults","CaseSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults","CaseSupportingSplitReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults","ControlSupportingReads")),
        int(config.get("ParametersToFilterDellyResults","ControlSupportingSplitReads")),
        int(config.get("ParametersToFilterDellyResults","ControlSupportingReadsHotspot")),
        int(config.get("ParametersToFilterDellyResults","ControlSupportingSplitReadsHotspot")),
        verbose)
'''