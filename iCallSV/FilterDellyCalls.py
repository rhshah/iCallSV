"""
Created on Mar 17, 2015
Description: This module will filter calls made by Delly which are in a VCF format
@author: Ronak H Shah

::Inputs::
inputVcf: Input VCF file name with path
outputDir: Output directory
controlId: Control Sample ID (Should be part of Sample Name in VCF)
caseID: Case Sample ID (Should be part of Sample Name in VCF)
hospotFile: List of Genes that have Hotspot Structural Variants (Tab-delimited Format without header:chr    start    end    geneName).
blacklistFile: List of Genes that have blacklist of Structural Variants (Tab-delimited Format without header:chr    start1    chr2     start2; where chr1==chr2, end==start2).
peSupport: overall pair-end read support threshold for the event
srSupport: overall split-reads support threshold for the event
peSupportHotspot: overall pair-end read support threshold for the event in hot-spot region
srSupportHotspot: overall split-reads support threshold for the event in hot-spot region
peSupportCase: pair-end read support threshold for the event in the Case sample
srSupportCase: split-reads support threshold for the event in the Case sample
peSupportHotspotCase: pair-end read support threshold for the event in hot-spot region for the Case sample
srSupportHotspotCase: split-reads support threshold for the event in hot-spot region for the Case sample
peSupportControl: pair-end read support threshold for the event in the Control sample
srSupportControl: split-reads support threshold for the event in the Control sample
peSupportHotspotControl: pair-end read support threshold for the event in hot-spot region for the Control sample
srSupportHotspotControl: split-reads support threshold for the event in hot-spot region for the Control sample
svlength: length of the structural variants
mapq: overall mapping quality
mapqHotspot: mapping quality for hot-spots

::Output::
Filtered VCF files
"""
import os
import sys
import vcf
import re
import checkparameters as cp
import checkHotSpotList as chl
import checkBlackList as cbl
import logging


def run(
        inputVcf,
        outputDir,
        controlId,
        caseID,
        hotspotFile,
        blacklistFile,
        verbose,
        svlength,
        mapq,
        mapqHotspot,
        peSupport,
        srSupport,
        peSupportHotspot,
        srSupportHotspot,
        peSupportCase,
        srSupportCase,
        peSupportHotspotCase,
        srSupportHotspotCase,
        peSupportControl,
        srSupportControl,
        peSupportHotspotControl,
        srSupportHotspotControl):
    if(verbose):
        logging.info("FilterDellyCalls: We will now check all the input parameters")
    # Check input parameters
    cp.checkDir(outputDir)
    cp.checkFile(inputVcf)
    cp.checkFile(hotspotFile)
    cp.checkFile(blacklistFile)
    cp.checkEmpty(controlId, "Control Bam ID")
    cp.checkEmpty(caseID, "Case Bam ID")
    cp.checkInt(svlength, "Length of Structural Variant Threshold")
    cp.checkInt(mapq, "Mapping quality of Reads threshold")
    cp.checkInt(mapqHotspot, "Mapping quality of Reads threshold for hotspot events ")
    cp.checkInt(peSupport, "overall pair-end read support threshold for the event")
    cp.checkInt(srSupport, "overall split-reads support threshold for the event")
    cp.checkInt(
        peSupportHotspot,
        "overall pair-end read support threshold for the event in hot-spot region")
    cp.checkInt(
        srSupportHotspot,
        "overall split-reads support threshold for the event in hot-spot region")
    cp.checkInt(
        peSupportCase,
        "overall pair-end read support threshold for the event for the Case sample")
    cp.checkInt(
        srSupportCase,
        "overall split-reads support threshold for the event for the Case sample")
    cp.checkInt(
        peSupportHotspotCase,
        "overall pair-end read support threshold for the event in hot-spot region for the Case sample")
    cp.checkInt(
        srSupportHotspotCase,
        "overall split-reads support threshold for the event in hot-spot region for the Case sample")
    cp.checkInt(
        peSupportControl,
        "overall pair-end read support threshold for the event for the Control sample")
    cp.checkInt(
        srSupportControl,
        "overall split-reads support threshold for the event for the Control sample")
    cp.checkInt(
        peSupportHotspotControl,
        "overall pair-end read support threshold for the event in hot-spot region for the Control sample")
    cp.checkInt(
        srSupportHotspotControl,
        "overall split-reads support threshold for the event in hot-spot region for the Control sample")
    if(verbose):
        logging.info("FilterDellyCalls: All Input Parameters look good for filtering these VCF file.")
        logging.info("FilterDellyCalls: We will filter the given VCF file now.")
    # Make a string of all the variables
    thresholdVariablesList = [svlength,
                              mapq,
                              mapqHotspot,
                              peSupport,
                              srSupport,
                              peSupportHotspot,
                              srSupportHotspot,
                              peSupportCase,
                              srSupportCase,
                              peSupportHotspotCase,
                              srSupportHotspotCase,
                              peSupportControl,
                              srSupportControl,
                              peSupportHotspotControl,
                              srSupportHotspotControl]
    thresholdVariables = ",".join(str(v) for v in thresholdVariablesList)

    hotspotDict = chl.ReadHotSpotFile(hotspotFile)
    blacklist = cbl.ReadBlackListFile(blacklistFile)
    outputVcf = os.path.splitext(os.path.basename(inputVcf))[0] + "_filtered.vcf"
    vcf_reader = vcf.Reader(open(inputVcf, 'r'))
    outputFile = outputDir + "/" + outputVcf
    vcf_writer = vcf.Writer(open(outputFile, 'w'), vcf_reader)
    samples = vcf_reader.samples
    pattern = re.compile(caseID)
    # Get the case and control id
    caseIDinVcf = None
    controlIDinVcf = None
    for sample in samples:
        match = re.search(pattern, sample)
        if(match):
            caseIDinVcf = sample
        else:
            controlIDinVcf = sample
    #
    for record in vcf_reader:
        # Define all variables:
        (chrom1,
         start1,
         start2,
         chrom2,
         filter,
         svlengthFromDelly,
         mapqFromDelly,
         svtype,
         preciseFlag,
         peSupportFromDelly,
         srSupportFromDelly,
         contype,
         caseDR,
         caseDV,
         caseRR,
         caseRV,
         caseFT,
         contrilDR,
         controlDV,
         controlRR,
         controlRV,
         controlFT) = (None for i in range(22))
        chrom1 = record.CHROM
        start1 = record.POS
        filter = record.FILTER
        preciseFlag = record.is_sv_precise
        if("END" in record.INFO):
            start2 = record.INFO['END']
        if("CHR2" in record.INFO):
            chrom2 = record.INFO['CHR2']
        if("SVTYPE" in record.INFO):
            svtype = record.INFO['SVTYPE']
        if("SVLEN" in record.INFO):
            svlengthFromDelly = record.INFO['SVLEN']
        else:
            if(svtype is not "TRA"):
                svlengthFromDelly = abs(start2 - start1)
            else:
                svlengthFromDelly = None
        if("MAPQ" in record.INFO):
            mapqFromDelly = record.INFO['MAPQ']

        if("PE" in record.INFO):
            peSupportFromDelly = record.INFO['PE']
        if("SR" in record.INFO):
            srSupportFromDelly = record.INFO['SR']
        if("CT" in record.INFO):
            contype = record.INFO['CT']
        caseCalls = record.genotype(caseIDinVcf)
        controlCalls = record.genotype(controlIDinVcf)
        if(hasattr(caseCalls.data, "FT")):
            caseFT = caseCalls.data.FT
        if(hasattr(caseCalls.data, "DR")):
            caseDR = caseCalls.data.DR
        if(hasattr(caseCalls.data, "DV")):
            caseDV = caseCalls.data.DV
        if(hasattr(caseCalls.data, "RR")):
            caseRR = caseCalls.data.RR
        if(hasattr(caseCalls.data, "RV")):
            caseRV = caseCalls.data.RV

        if(hasattr(controlCalls.data, "FT")):
            controlFT = controlCalls.data.FT
        if(hasattr(controlCalls.data, "DR")):
            controlDR = controlCalls.data.DR
        if(hasattr(controlCalls.data, "DV")):
            controlDV = controlCalls.data.DV
        if(hasattr(controlCalls.data, "RR")):
            controlRR = controlCalls.data.RR
        if(hasattr(controlCalls.data, "RV")):
            controlRV = controlCalls.data.RV
        # Make a string of all the variables
        dellyVariablesList = [chrom1,
                              start1,
                              start2,
                              chrom2,
                              filter,
                              svlengthFromDelly,
                              mapqFromDelly,
                              svtype,
                              preciseFlag,
                              peSupportFromDelly,
                              srSupportFromDelly,
                              contype,
                              caseFT,
                              caseDR,
                              caseDV,
                              caseRR,
                              caseRV,
                              controlFT,
                              controlDR,
                              controlDV,
                              controlRR,
                              controlRV]
        dellyVariables = ",".join(str(v) for v in dellyVariablesList)
        # print chrom1, start1, start2, chrom2, svlengthFromDelly, mapqFromDelly,
        # svtype, peSupportFromDelly, srSupportFromDelly, contype, caseDV, caseRV,
        # controlDV, controlRV
        filterFlag = GetFilteredRecords(dellyVariables, thresholdVariables, hotspotDict, blacklist)
        if(filterFlag):
            vcf_writer.write_record(record)
    vcf_writer.close()
    if(verbose):
        logging.info("FilterDellyCalls: We have finished filtering: %s file", inputVcf)
        logging.info("FilterFellyCalls: Output hass been written in: %s file", outputFile)
    return(outputFile)


def GetFilteredRecords(dellyVarialbles, thresholdVariables, hotspotDict, blacklist):
    (svlength,
     mapq,
     mapqHotspot,
     peSupport,
     srSupport,
     peSupportHotspot,
     srSupportHotspot,
     peSupportCase,
     srSupportCase,
     peSupportHotspotCase,
     srSupportHotspotCase,
     peSupportControl,
     srSupportControl,
     peSupportHotspotControl,
     srSupportHotspotControl) = thresholdVariables.split(",")
    (chrom1,
     start1,
     start2,
     chrom2,
     filter,
     svlengthFromDelly,
     mapqFromDelly,
     svtype,
     preciseFlag,
     peSupportFromDelly,
     srSupportFromDelly,
     contype,
     caseFT,
     caseDR,
     caseDV,
     caseRR,
     caseRV,
     controlFT,
     controlDR,
     controlDV,
     controlRR,
     controlRV) = dellyVarialbles.split(",")
    # Get if its a hotspot or not
    hotspotTag = chl.CheckIfItIsHotspot(chrom1, start1, chrom2, start2, hotspotDict)
    # Get if its a blacklist or not
    blacklistTag = cbl.CheckIfItIsBlacklisted(chrom1, start1, chrom2, start2, blacklist, 20)
    # Get the flag for pass and fail for tumor and normal
    casePassFlag = GetCaseFlag(caseDR, caseDV, preciseFlag, caseRR, caseRV)
    controlPassFlag = GetControlFlag(controlDR, controlDV, preciseFlag, controlRR, controlRV)
    filterFlag = False
    if(casePassFlag and controlPassFlag):
        if(hotspotTag):
            if(filter == "PASS" and controlFT == "LowQual"):
                if(svlengthFromDelly != "None"):
                    if(int(svlengthFromDelly) >= int(svlength)):
                        filterFlag = True
                    else:
                        filterFlag = False
                else:
                    filterFlag = True
            if(not filterFlag):
                if(svlengthFromDelly != "None"):
                    if((int(svlengthFromDelly) >= int(svlength)) and (int(mapqFromDelly) >= int(mapqHotspot)) and (int(peSupportFromDelly) >= int(peSupportHotspot)) and (int(caseDV) > int(peSupportHotspotCase)) and (int(controlDV) <= int(peSupportHotspotControl)) and (int(controlDV) < int(caseDV))):
                        if(preciseFlag):
                            if((int(srSupportFromDelly) >= int(srSupportHotspot)) and (int(caseRV) >= int(srSupportHotspotCase)) and (int(controlRV) <= int(srSupportHotspotControl)) and (int(controlRV) < int(caseRV))):
                                filterFlag = True
                            else:
                                filterFlag = False
                        else:
                            filterFlag = True
                    else:
                        filterFlag = False
                else:
                    if((int(mapqFromDelly) >= int(mapqHotspot)) and (int(peSupportFromDelly) >= int(peSupportHotspot)) and (int(caseDV) >= int(peSupportHotspotCase)) and (int(controlDV) <= int(peSupportHotspotControl)) and (int(controlDV) < int(caseDV))):
                        if(preciseFlag):
                            if((int(srSupportFromDelly) >= int(srSupportHotspot)) and (int(caseRV) >= int(srSupportHotspotCase)) and (int(controlRV) <= int(srSupportHotspotControl)) and (int(controlRV) < int(caseRV))):
                                filterFlag = True
                            else:
                                filterFlag = False
                        else:
                            filterFlag = True
                    else:
                        filterFlag = False
        else:
            if(svlengthFromDelly != "None"):
                # print svlengthFromDelly, svlength, mapqFromDelly, mapq,
                # peSupportFromDelly, peSupport, caseDV, peSupportCase, controlDV,
                # peSupportControl
                if((int(svlengthFromDelly) >= int(svlength)) and (int(mapqFromDelly) >= int(mapq)) and (int(peSupportFromDelly) >= int(peSupport)) and (int(caseDV) >= int(peSupportCase)) and (int(controlDV) <= int(peSupportControl)) and (int(controlDV) < int(caseDV))):
                    if(preciseFlag):
                        if((int(srSupportFromDelly) >= int(srSupport)) and (int(caseRV) >= int(srSupportCase)) and (int(controlRV) <= int(srSupportControl)) and (int(controlRV) < int(caseRV))):
                            filterFlag = True
                        else:
                            filterFlag = False
                    else:
                        filterFlag = True
                else:
                    filterFlag = False
            else:
                if((int(mapqFromDelly) >= int(mapq)) and (int(peSupportFromDelly) >= int(peSupport)) and (int(caseDV) >= int(peSupportCase)) and (int(controlDV) <= int(peSupportControl)) and (int(controlDV) < int(caseDV))):
                    if(preciseFlag):
                        if((int(srSupportFromDelly) >= int(srSupport)) and (int(caseRV) >= int(srSupportCase)) and (int(controlRV) <= int(srSupportControl)) and (int(controlRV) < int(caseRV))):
                            filterFlag = True
                        else:
                            filterFlag = False
                    else:
                        filterFlag = True
                else:
                    filterFlag = False
    else:
        filterFlag = True

    if(blacklistTag):
        filterFlag = True
    else:
        filterFlag = filterFlag

    return(filterFlag)


def GetCaseFlag(caseDR, caseDV, preciseFlag, caseRR, caseRV):
    caseAltAf = 0.0
    caseCovg = 0
    caseFlag = False
    if(preciseFlag):
        caseCovg = int(caseRR) + int(caseRV)
        if((caseRR is not None or caseRR != 0) and (caseRV is not None or caseRV != 0)):
            caseAltAf = caseRV / float(caseRR + caseRV)

    else:
        caseCovg = int(caseDR) + int(caseDV)
        if((caseDR is not None or caseDR != 0) and (caseDV is not None or caseDV != 0)):
            caseAltAf = caseDV / float(caseDR + caseDV)

    if(caseAltAf >= 0.2 and caseCovg >= 10):
        caseFlag = True
    else:
        caseFlag = False
    return(caseFlag)


def GetControlFlag(controlDR, controlDV, preciseFlag, controlRR, controlRV):
    controlAltAf = 0.0
    controlCovg = 0
    controlFlag = False
    if(preciseFlag):
        if((controlRR is not None or controlRR != 0) and (controlRV is not None or controlRV != 0)):
            controlAltAf = controlRV / float(controlRR + controlRV)

    else:
        if((controlDR is not None or controlDR != 0) and (controlDV is not None or controlDV != 0)):
            controlAltAf = controlDV / float(controlDR + controlDV)

    if(controlAltAf <= 0.0):
        controlFlag = True
    else:
        controlFlag = False
    return(controlFlag)

# # Test the module
# run(
#     inputVcf="/ifs/e63data/bergerm1/Analysis/Projects/Test/Ronak/Test/iCallSV/StructuralVariantAnalysis/DellyDir/s-EV-crc-001-M2/s-EV-crc-001-M2_del.vcf",
#     outputDir="/ifs/e63data/bergerm1/Analysis/Projects/Test/Ronak/Test/iCallSV/StructuralVariantAnalysis/DellyDir/s-EV-crc-001-M2/",
#     controlId="s-EV-crc-001-N2",
#     caseID="s-EV-crc-001-M2",
#     hotspotFile="/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv5/structuralvariants_geneInterval.txt",
#     blacklistFile = "/home/shahr2/git/iCallSV/iCallSV/data/blacklist.txt",
#     svlength=500,
#     mapq=20,
#     mapqHotspot=5,
#     peSupport=5,
#     srSupport=0,
#     peSupportHotspot=3,
#     srSupportHotspot=0,
#     peSupportCase=1,
#     srSupportCase=0,
#     peSupportHotspotCase=0,
#     srSupportHotspotCase=0,
#     peSupportControl=5,
#     srSupportControl=5,
#     peSupportHotspotControl=5,
#     srSupportHotspotControl=5,
#     verbose=True)
