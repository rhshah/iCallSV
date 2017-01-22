"""
FilterDellyCalls
~~~~~~~~~~~~~~~~

:Description: This module will filter calls made by Delly which are in a VCF format
"""
'''
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
caseAltFreq: Alternate Allele Frequency threshold for case
caseTotalCount: Total ReadCount threshold for case
ccontrolAltFreq: Alternate Allele Frequency threshold for control
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
'''
import os
import sys
import vcf
import re
import checkparameters as cp
import checkHotSpotList as chl
import checkBlackList as cbl
import logging
import coloredlogs

logger = logging.getLogger('iCallSV.FilterDellyCalls')
coloredlogs.install(level='DEBUG')

def run(
        inputVcf,
        outputDir,
        controlId,
        caseID,
        hotspotFile,
        blacklistFile,
        svlength,
        mapq,
        mapqHotspot,
        caseAltFreqHotspot,
        caseTotalCountHotspot,
        controlAltFreqHotspot,
        caseAltFreq,
        caseTotalCount,
        controlAltFreq,
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
        srSupportHotspotControl,
        verbose):
    """``main:``Filter calls made by Delly which are in a VCF format


    :param str inputVcf: Input VCF file name with path
    :param str outputDir: Output directory
    :param str controlId: Control Sample ID (Should be part of Sample Name in VCF)
    :param str caseID: Case Sample ID (Should be part of Sample Name in VCF)
    :param str hospotFile: List of Genes that have Hotspot Structural Variants (Tab-delimited Format without header:chr    start    end    geneName).
    :param str blacklistFile: List of Genes that have blacklist of Structural Variants (Tab-delimited Format without header:chr    start1    chr2     start2; where chr1==chr2, end==start2).
    :param float caseAltFreq: Alternate Allele Frequency threshold for case
    :param int caseTotalCount: Total ReadCount threshold for case
    :param flaot ccontrolAltFreq: Alternate Allele Frequency threshold for control
    :param int peSupport: overall pair-end read support threshold for the event
    :param int srSupport: overall split-reads support threshold for the event
    :param int peSupportHotspot: overall pair-end read support threshold for the event in hot-spot region
    :param int srSupportHotspot: overall split-reads support threshold for the event in hot-spot region
    :param int peSupportCase: pair-end read support threshold for the event in the Case sample
    :param int srSupportCase: split-reads support threshold for the event in the Case sample
    :param int peSupportHotspotCase: pair-end read support threshold for the event in hot-spot region for the Case sample
    :param int srSupportHotspotCase: split-reads support threshold for the event in hot-spot region for the Case sample
    :param int peSupportControl: pair-end read support threshold for the event in the Control sample
    :param int srSupportControl: split-reads support threshold for the event in the Control sample
    :param int peSupportHotspotControl: pair-end read support threshold for the event in hot-spot region for the Control sample
    :param int srSupportHotspotControl: split-reads support threshold for the event in hot-spot region for the Control sample
    :param int svlength: length of the structural variants
    :param int mapq: overall mapping quality
    :param int mapqHotspot: mapping quality for hot-spots
    :return: A str name of filtered vcf file
    :rtype: str
    """
    if(verbose):
        logger.info("FilterDellyCalls: We will now check all the input parameters")
    # Check input parameters
    cp.checkDir(outputDir)
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
        logger.info("FilterDellyCalls: All Input Parameters look good for filtering these VCF file.")
        logger.info("FilterDellyCalls: We will filter the given VCF file now.")
    # Make a string of all the variables
    thresholdVariablesList = [svlength,
                              mapq,
                              mapqHotspot,
                              caseAltFreqHotspot,
                              caseTotalCountHotspot,
                              controlAltFreqHotspot,
                              caseAltFreq,
                              caseTotalCount,
                              controlAltFreq,
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
    outputFile = os.path.join(outputDir, outputVcf)
    if(not os.path.isfile(inputVcf)):
        if(verbose):
            logger.warning("VCF file %s does not exists.", inputVcf)
        return(outputFile)
    vcf_reader = vcf.Reader(open(inputVcf, 'r'))
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
    # Check if ID are assigned properly or not
    if(caseIDinVcf is None):
        logger.error("FilterDellyCalls: caseID was not assigned properly, please make sure that the vcf case id and the provided case id match")
        sys.exit(1)
    else:
        if(verbose):
            logger.info("FilterDellyCalls:Case ID is: %s file", caseIDinVcf)
    if(controlIDinVcf is None):
        logger.error("FilterDellyCalls: controlID was not assigned properly, please make sure that the vcf control id and the provided control id match")
        sys.exit(1)
    else:
        if(verbose):
            logger.info("FilterDellyCalls:Control ID is: %s file", controlIDinVcf)
    
    # Traversing the VCF
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
         controlDR,
         controlDV,
         controlRR,
         controlRV,
         controlFT) = (None for i in range(22))
        chrom1 = record.CHROM
        start1 = record.POS
        filter = record.FILTER
        if(len(filter) < 1):
            filter = "PASS"
        else:
            filter = filter[0]
        preciseFlag = record.is_sv_precise
        # print "Precise:",preciseFlag,":",type(preciseFlag)
        if("END" in record.INFO):
            start2 = record.INFO['END']
        if("CHR2" in record.INFO):
            chrom2 = record.INFO['CHR2']
        if("SVTYPE" in record.INFO):
            svtype = record.INFO['SVTYPE']
        if("SVLEN" in record.INFO):
            svlengthFromDelly = record.INFO['SVLEN']
        else:
            if(svtype == "TRA"):
                svlengthFromDelly = None
            else:
                svlengthFromDelly = abs(start2 - start1)
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
        # print chrom1, start1, start2, chrom2, "Coordinate"
        # print svlengthFromDelly, mapqFromDelly, svtype, peSupportFromDelly, srSupportFromDelly, contype, "Overall"
        # print caseDR, caseDV, caseRR, caseRV, "Case"
        # print controlDR, controlDV, controlRR, controlRV, "Control"
        filterFlag = GetFilteredRecords(dellyVariables, thresholdVariables, hotspotDict, blacklist)
        if(filterFlag):
            # print "Passs"
            vcf_writer.write_record(record)
    vcf_writer.close()
    if(verbose):
        logger.info("FilterDellyCalls: We have finished filtering: %s file", inputVcf)
        logger.info("FilterFellyCalls: Output hass been written in: %s file", outputFile)
    return(outputFile)


def GetFilteredRecords(dellyVarialbles, thresholdVariables, hotspotDict, blacklist):
    """
    This will ``Filter one record at a time``

    :param str dellyVariables: str having all delly variables separated by ","
    :param str thresholdVariables: str having all delly threshold variables separated by ","
    :param dict hotspotDict: A dict containing hotspot regions
    :param list blacklist: A list containing blacklist regions
    :return: A boolean tag indicating True or False
    :rtype: bool

    """
    (svlength,
     mapq,
     mapqHotspot,
     caseAltFreqHotspot,
     caseTotalCountHotspot,
     controlAltFreqHotspot,
     caseAltFreq,
     caseTotalCount,
     controlAltFreq,
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
    # print "PreciseFlagInFilter:",preciseFlag,":",type(preciseFlag)
    # Get if its a hotspot or not
    hotspotTag = chl.CheckIfItIsHotspot(chrom1, start1, chrom2, start2, hotspotDict)
    # Get if its a blacklist or not
    blacklistTag = cbl.CheckIfItIsBlacklisted(chrom1, start1, chrom2, start2, blacklist, 20)
    # Get the flag for pass and fail for tumor and normal
    filterFlag = False
    # print "CaseControlPassFlag:", casePassFlag, " : ", controlPassFlag

    if(hotspotTag):
        casePassFlag = GetCaseFlag(
            caseDR,
            caseDV,
            preciseFlag,
            caseRR,
            caseRV,
            caseAltFreqHotspot,
            caseTotalCountHotspot)
        controlPassFlag = GetControlFlag(
            controlDR,
            controlDV,
            preciseFlag,
            controlRR,
            controlRV,
            controlAltFreqHotspot)
        if(filter == "PASS" and controlPassFlag and casePassFlag):
            if(svlengthFromDelly != "None"):
                if(int(svlengthFromDelly) >= int(svlength)):
                    filterFlag = True
                    return(filterFlag)
                else:
                    filterFlag = False
            else:
                filterFlag = True
                return(filterFlag)

            if(filterFlag is False):
                if(svlengthFromDelly != "None"):
                    if((int(svlengthFromDelly) >= int(svlength)) and (int(mapqFromDelly) >= int(mapqHotspot)) and (int(peSupportFromDelly) >= int(peSupportHotspot)) and (int(caseDV) > int(peSupportHotspotCase)) and (int(controlDV) <= int(peSupportHotspotControl)) and (int(controlDV) < int(caseDV))):
                        if(preciseFlag == "True"):
                            if((int(srSupportFromDelly) >= int(srSupportHotspot)) and (int(caseRV) >= int(srSupportHotspotCase)) and (int(controlRV) <= int(srSupportHotspotControl)) and (int(controlRV) < int(caseRV))):
                                filterFlag = True
                                return(filterFlag)
                            else:
                                filterFlag = False
                        else:
                            filterFlag = True
                            return(filterFlag)
                    else:
                        filterFlag = False
                else:
                    if((int(mapqFromDelly) >= int(mapqHotspot)) and (int(peSupportFromDelly) >= int(peSupportHotspot)) and (int(caseDV) >= int(peSupportHotspotCase)) and (int(controlDV) <= int(peSupportHotspotControl)) and (int(controlDV) < int(caseDV))):
                        if(preciseFlag == "True"):
                            if((int(srSupportFromDelly) >= int(srSupportHotspot)) and (int(caseRV) >= int(srSupportHotspotCase)) and (int(controlRV) <= int(srSupportHotspotControl)) and (int(controlRV) < int(caseRV))):
                                filterFlag = True
                                return(filterFlag)
                            else:
                                filterFlag = False
                        else:
                            filterFlag = True
                            return(filterFlag)
                    else:
                        filterFlag = False
            else:
                filterFlag = False
        else:
            filterFlag = False
    else:
        casePassFlag = GetCaseFlag(caseDR, caseDV, preciseFlag, caseRR,
                                   caseRV, caseAltFreq, caseTotalCount)
        controlPassFlag = GetControlFlag(
            controlDR,
            controlDV,
            preciseFlag,
            controlRR,
            controlRV,
            controlAltFreq)
        if(filter == "PASS" and controlPassFlag and casePassFlag):
            if(svlengthFromDelly != "None"):
                if(int(svlengthFromDelly) >= int(svlength)):
                    filterFlag = True
                    return(filterFlag)
                else:
                    filterFlag = False
            else:
                filterFlag = True
                return(filterFlag)
        if(filterFlag is False):
            if(svlengthFromDelly != "None"):
                # print svlengthFromDelly, svlength, mapqFromDelly, mapq,
                # peSupportFromDelly, peSupport, caseDV, peSupportCase, controlDV,
                # peSupportControl
                if((int(svlengthFromDelly) >= int(svlength)) and (int(mapqFromDelly) >= int(mapq)) and (int(peSupportFromDelly) >= int(peSupport)) and (int(caseDV) >= int(peSupportCase)) and (int(controlDV) <= int(peSupportControl)) and (int(controlDV) < int(caseDV))):
                    if(preciseFlag == "True"):
                        if((int(srSupportFromDelly) >= int(srSupport)) and (int(caseRV) >= int(srSupportCase)) and (int(controlRV) <= int(srSupportControl)) and (int(controlRV) < int(caseRV))):
                            filterFlag = True
                            return(filterFlag)
                        else:
                            filterFlag = False
                    else:
                        filterFlag = True
                        return(filterFlag)
                else:
                    filterFlag = False
            else:
                if((int(mapqFromDelly) >= int(mapq)) and (int(peSupportFromDelly) >= int(peSupport)) and (int(caseDV) >= int(peSupportCase)) and (int(controlDV) <= int(peSupportControl)) and (int(controlDV) < int(caseDV))):
                    if(preciseFlag == "True"):
                        if((int(srSupportFromDelly) >= int(srSupport)) and (int(caseRV) >= int(srSupportCase)) and (int(controlRV) <= int(srSupportControl)) and (int(controlRV) < int(caseRV))):
                            filterFlag = True
                            return(filterFlag)
                        else:
                            filterFlag = False
                    else:
                        filterFlag = True
                        return(filterFlag)
                else:
                    filterFlag = False
        else:
            filterFlag = False

        if(blacklistTag):
            filterFlag = True
            return(filterFlag)
        else:
            filterFlag = filterFlag

    # return(filterFlag)


def GetCaseFlag(caseDR, caseDV, preciseFlag, caseRR, caseRV, caseAltFreq, caseTotalCount):
    """
    This will ``check if the case sample passes or not``

    :param int caseDR: int representing number of reference reads for case reported by delly
    :param int caseDV: int representing number of variant reads for case reported by delly
    :param str preciseFlag: str representing if an event is precise or imprecise
    :param int caseRR: int representing number of split reference reads for case reported by delly
    :param int caseRV: int representing number of split variant reads for case reported by delly
    :param float caseAltFreq: float representing altratio threshold for case
    :param int caseTotalCount: int repeseting readcount threshold for case
    :return: A boolean tag indicating True or False
    :rtype: bool

    """
    caseAltAf = 0.0
    caseCovg = 0
    caseFlag = False
    if(caseDR is None):
        caseDR = 0
    if(caseDV is None):
        caseDV = 0
    if(caseRR is None):
        caseRR = 0
    if(caseRV is None):
        caseRR = 0
    if(preciseFlag == "True"):
        caseCovg = int(caseRR) + int(caseRV)
        if((float(caseRR) != 0.0) or (float(caseRV) != 0.0)):
            caseAltAf = float(caseRV) / float(int(caseRR) + int(caseRV))

    else:
        caseCovg = int(caseDR) + int(caseDV)
        if((float(caseDR) != 0.0) or (float(caseDV) != 0.0)):
            caseAltAf = float(caseDV) / float(int(caseDR) + int(caseDV))

    if(caseAltAf >= float(caseAltFreq) and caseCovg >= int(caseTotalCount)):
        caseFlag = True
    elif((preciseFlag == "True") and (caseRV >= 5) and (caseCovg >= int(caseTotalCount))):
        caseFlag = True
    else:
        caseFlag = False
    return(caseFlag)


def GetControlFlag(controlDR, controlDV, preciseFlag, controlRR, controlRV, controlAltFreq):
    """
    This will ``check if the control sample passes or not``

    :param int controlDR: int representing number of reference reads for control reported by delly
    :param int controlDV: int representing number of variant reads for control reported by delly
    :param str preciseFlag: str representing if an event is precise or imprecise
    :param int controlRR: int representing number of split reference reads for control reported by delly
    :param int controlRV: int representing number of split variant reads for control reported by delly
    :param float controlAltFreq: float representing altratio threshold for control
    :return: A boolean tag indicating True or False
    :rtype: bool

    """
    controlAltAf = 0.0
    controlCovg = 0
    controlFlag = False
    if(controlDR is None):
        controlDR = 0
    if(controlDV is None):
        controlDV = 0
    if(controlRR is None):
        controlRR = 0
    if(controlRV is None):
        controlRR = 0
    if(preciseFlag == "True"):
        if((float(controlRR) != 0.0) or (float(controlRV) != 0.0)):
            controlAltAf = float(controlRV) / float(int(controlRR) + int(controlRV))

    else:
        if((float(controlDR) != 0.0) or (float(controlDV) != 0.0)):
            controlAltAf = float(controlDV) / float(int(controlDR) + int(controlDV))

    if(controlAltAf <= float(controlAltFreq)):
        controlFlag = True
    else:
        controlFlag = False
    return(controlFlag)
