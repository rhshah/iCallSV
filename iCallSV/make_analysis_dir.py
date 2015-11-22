'''
Created on November 19, 2015
Description: This module will make directory structure for running analysis
@author: Ronak H Shah
'''
'''
::Inputs::
args: Arguments passsed to iCallSV
'''
import os
import re
import logging

def makeOutputDir(args,tool):
    SampleDirName = args.patientId
    static_SV_Dir = "StructuralVariantAnalysis"
    static_analysis_Dir = tool
    AnalysisDir = os.path.join(args.outdir, static_SV_Dir)
    ToolDir = os.path.join(AnalysisDir, static_Delly_Dir)
    SampleAnalysisDir = os.path.join(ToolDir, SampleDirName)
    try:
        os.mkdir(AnalysisDir)
    except OSError:
        if(args.verbose):
            logging.warn("Dir:", AnalysisDir, " exists thus we wont be making it")
    try:
        os.mkdir(ToolDir)
    except OSError:
        if(args.verbose):
             logging.warn("Dir:", ToolDir, " exists thus we wont be making it.")
    
    if os.path.isdir(SampleAnalysisDir):
            if(args.verbose):
                logging.info("Dir:", SampleAnalysisDir, " exists and we wont run the analysis")
            tag = 0
    else:
        os.mkdir(SampleAnalysisDir)
        tag = 1
return(tag,SampleAnalysisDir)