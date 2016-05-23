"""
Created on November 19, 2015
Description: This module will make directory structure for running analysis
@author: Ronak H Shah
"""
"""
::Inputs::
args: Arguments passed to iCallSV
"""

import os
import re
import logging


def makeOutputDir(args, tool):
    logger = logging.getLogger(__name__)
    SampleDirName = args.caseId
    static_SV_Dir = "StructuralVariantAnalysis"
    static_tool_Dir = tool
    AnalysisDir = os.path.join(args.outdir, static_SV_Dir)
    ToolDir = os.path.join(AnalysisDir, static_tool_Dir)
    SampleAnalysisDir = os.path.join(ToolDir, SampleDirName)
    try:
        os.mkdir(AnalysisDir)
    except OSError:
        if(args.verbose):
            logger.warn("make_output_dir:Dir=>%s exists thus we wont be making it", AnalysisDir)
    try:
        os.mkdir(ToolDir)
    except OSError:
        if(args.verbose):
            logger.warn("make_output_dir:Dir=>%s exists thus we wont be making it", ToolDir)

    if os.path.isdir(SampleAnalysisDir):
        if(args.verbose):
            logger.fatal(
                "make_output_dir:Dir=>%s exists and we wont run the analysis",
                SampleAnalysisDir)
            logger.info("make_output_dir:Please delete this directory and rerun the program")
        tag = False
    else:
        os.mkdir(SampleAnalysisDir)
        tag = True

    return(tag, SampleAnalysisDir)
