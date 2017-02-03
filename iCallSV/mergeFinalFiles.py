"""
mergeFinalFiles
~~~~~~~~~~~~~~~

:Description: Merge VCF, iAnnotateSV tab and targetSeqView tab file into a single tab-delimited file
"""
'''
Created on May 17, 2015
Description: Merge VCF, iAnnotateSV tab and targetSeqView tab file into a single tab-delimited file
@author: Ronak H Shah
::Input::
aId: Sample ID for case that has the structural abberations
bId: Sample ID for control
vcfFile: Delly filtered and merged VCF file
annoTab: iAnnotateSV tab-delimited file with annotations
confTab: targetSeqView tab-delimited file with probability score
outputDir: Directory to write the output file
outputPrefix: Output File Prefix
::Output::
outputFile: File with following header
"TumorId\tNormalId\tChr1\tPos1\tChr2\tPos2\tSV_Type\tGene1\tGene2\tTranscript1\tTranscript2\tSite1Description\tSite2Description\tFusion\tProbabilityScore\tConfidence\tComments\tConnection_Type\tSV_LENGTH\tMAPQ\tPairEndReadSupport\tSplitReadSupport\tBrkptType\tConsensusSequence\tTumorVariantCount\tTumorSplitVariantCount\tTumorReadCount\tTumorGenotypeQScore\tNormalVariantCount\tNormalSplitVariantCount\tNormalReadCount\tNormalGenotypeQScorerepName-repClass-repFamily:-site1\trepName-repClass-repFamily:-site2\tCC_Chr_Band\tCC_Tumour_Types(Somatic)\tCC_Cancer_Syndrome\tCC_Mutation_Type\tCC_Translocation_Partner\tDGv_Name-DGv_VarType-site1\tDGv_Name-DGv_VarType-site2\n";
'''
import sys
import os
import logging
import vcf
import checkparameters as cp
import re
import coloredlogs
import numpy as np
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp() #So that matplotlib doesnot complain stale file handle
try:
    import pandas as pd
except ImportError, e:
    print "mergeFinalFiles: pandas is not installed, please install pandas as it is required to run the mapping."
    sys.exit(1)
    
logger = logging.getLogger('iCallSV.mergeFinalFiles')
coloredlogs.install(level='DEBUG')


def run(aId, bId, vcfFile, annoTab, confTab, outDir, outputPrefix, verbose):
    """
    This will Merge VCF, iAnnotateSV tab and targetSeqView tab file into a single tab-delimited file

    :param str aId: Sample ID for case that has the structural abberations
    :param str bId: Sample ID for control
    :param str vcfFile: Delly filtered and merged VCF file
    :param str annoTab: iAnnotateSV tab-delimited file with annotations
    :param str confTab: targetSeqView tab-delimited file with probability score
    :param str outputDir: Directory to write the output file
    :param str outputPrefix: Output File Prefix
    :return: str of the tab-delimited file
    :rtype: str

    """
    if(verbose):
        logger.info(
            "iCallSV::MergeFinalFile: Merging Delly Filtered VCF, iAnnotateSV tab and targetSeqView tab file into a single tab-delimited file")
    cp.checkFile(vcfFile)
    cp.checkFile(annoTab)
    # cp.checkFile(confTab)
    cp.checkDir(outDir)
    outDF = pd.DataFrame(
        columns=[
            "TumorId",
            "NormalId",
            "Chr1",
            "Pos1",
            "Chr2",
            "Pos2",
            "SV_Type",
            "Gene1",
            "Gene2",
            "Transcript1",
            "Transcript2",
            "Site1Description",
            "Site2Description",
            "Fusion",
            "ProbabilityScore",
            "Confidence",
            "Comments",
            "Connection_Type",
            "SV_LENGTH",
            "MAPQ",
            "PairEndReadSupport",
            "SplitReadSupport",
            "BrkptType",
            "ConsensusSequence",
            "TumorReferenceCount",
            "TumorSplitReferenceCount",
            "TumorVariantCount",
            "TumorSplitVariantCount",
            "TumorReadCount",
            "TumorGenotypeQScore",
            "NormalReferenceCount",
            "NormalSplitReferenceCount",
            "NormalVariantCount",
            "NormalSplitVariantCount",
            "NormalReadCount",
            "NormalGenotypeQScore",
            "Cosmic_Fusion_Counts",
            "repName-repClass-repFamily:-site1",
            "repName-repClass-repFamily:-site2",
            "CC_Chr_Band",
            "CC_Tumour_Types(Somatic)",
            "CC_Cancer_Syndrome",
            "CC_Mutation_Type",
            "CC_Translocation_Partner",
            "DGv_Name-DGv_VarType-site1",
            "DGv_Name-DGv_VarType-site2"])
    annoDF = pd.read_csv(annoTab, sep="\t", header=0, keep_default_na='True')
    if(os.path.isfile(confTab)):
        confDF = pd.read_csv(confTab, sep="\t", header=0, keep_default_na='True')
    else:
        confDF = None
    # Read VCF and Traverse through it
    vcf_reader = vcf.Reader(open(vcfFile, 'r'))
    samples = vcf_reader.samples
    pattern = re.compile(aId)
    # Get the case and control id
    caseIDinVcf = None
    controlIDinVcf = None
    for sample in samples:
        match = re.search(pattern, sample)
        if(match):
            caseIDinVcf = sample
        else:
            controlIDinVcf = sample
    # traverse through the vcf
    count = 0
    for record in vcf_reader:
        # Define all variables:
        (chrom1,
         start1,
         start2,
         chrom2,
         filter,
         svtype,
         brktype,
         contype,
         conseq) = (None for i in range(9))
         
        (startCT,
         endCT,
         str1,
         str2,
         svlengthFromDelly,
         mapqFromDelly,
         peSupportFromDelly,
         srSupportFromDelly,
         ciEndNeg,
         ciEndPos,
         ciPosNeg,
         ciPosPos,
         caseRC,
         caseGQ,
         caseDR,
         caseDV,
         caseRR,
         caseRV,
         controlGQ,
         controlRC,
         controlDR,
         controlDV,
         controlRR,
         controlRV) = (0 for i in range(24))
        chrom1 = str(record.CHROM)
        start1 = record.POS
        filter = record.FILTER
        if(len(filter) < 1):
            filter = None
        else:
            filter = filter[0]
        preciseFlag = record.is_sv_precise
        if("END" in record.INFO):
            start2 = record.INFO['END']
        if("CHR2" in record.INFO):
            chrom2 = str(record.INFO['CHR2'])
        if("SVTYPE" in record.INFO):
            svtype = record.INFO['SVTYPE']
        if("SVLEN" in record.INFO):
            svlengthFromDelly = np.int(record.INFO['SVLEN'])
        else:
            if(svtype == "TRA"):
                svlengthFromDelly = 0
            else:
                svlengthFromDelly = np.int(abs(start2 - start1))
        if("MAPQ" in record.INFO):
            mapqFromDelly =  np.int(record.INFO['MAPQ'])
        if("PE" in record.INFO):
            peSupportFromDelly = np.int(record.INFO['PE'])
        if("SR" in record.INFO):
            srSupportFromDelly = np.int(record.INFO['SR'])
        if("CT" in record.INFO):
            contype = record.INFO['CT']
        (startCT, endCT) = contype.split("to")
        if((int(startCT) == 3) and (int(endCT) == 3)):
            str1 = 0
            str2 = 0
        elif((int(startCT) == 3) and (int(endCT) == 5)):
            str1 = 0
            str2 = 1
        elif((int(startCT) == 5) and (int(endCT) == 3)):
            str1 = 1
            str2 = 0
        elif((int(startCT) == 5) and (int(endCT) == 5)):
            str1 = 1
            str2 = 1
        else:
            if(verbose):
                logger.warning(
                    "mergeFinalFiles: The connection type (CT) given in the vcf file is incorrect.CT: %s",
                    contype)
        if("CONSENSUS" in record.INFO):
            conseq = record.INFO['CONSENSUS']
        if(record.is_sv_precise):
            brktype = "PRECISE"
        else:
            brktype = "IMPPRECISE"
        if("CIEND" in record.INFO):
            ciEndNeg, ciEndPos = record.INFO['CIEND']
        if(abs(ciEndNeg) < 50):
            ciEndNeg = 50
        if(abs(ciEndPos) < 50):
            ciEndNeg = 50
        if("CIPOS" in record.INFO):
            ciPosNeg, ciPosPos = record.INFO['CIPOS']
        if(abs(ciPosNeg) < 50):
            ciPosNeg = 50
        if(abs(ciPosPos) < 50):
            ciPosNeg = 50

        caseCalls = record.genotype(caseIDinVcf)
        controlCalls = record.genotype(controlIDinVcf)

        if(hasattr(caseCalls.data, "GQ")):
            caseGQ = np.int(caseCalls.data.GQ)
        if(hasattr(caseCalls.data, "RC")):
            caseRC = np.int(caseCalls.data.RC)
        if(hasattr(caseCalls.data, "DR")):
            caseDR = np.int(caseCalls.data.DR)
        if(hasattr(caseCalls.data, "DV")):
            caseDV = np.int(caseCalls.data.DV)
        if(hasattr(caseCalls.data, "RR")):
            caseRR = np.int(caseCalls.data.RR)
        if(hasattr(caseCalls.data, "RV")):
            caseRV = np.int(caseCalls.data.RV)

        if(hasattr(controlCalls.data, "GQ")):
            controlGQ = np.int(controlCalls.data.GQ)
        if(hasattr(controlCalls.data, "RC")):
            controlRC = np.int(controlCalls.data.RC)
        if(hasattr(controlCalls.data, "DR")):
            controlDR = np.int(controlCalls.data.DR)
        if(hasattr(controlCalls.data, "DV")):
            controlDV = np.int(controlCalls.data.DV)
        if(hasattr(controlCalls.data, "RR")):
            controlRR = np.int(controlCalls.data.RR)
        if(hasattr(controlCalls.data, "RV")):
            controlRV = np.int(controlCalls.data.RV)

        # Get data from annotation file
        (indexList,
         annoIndex,
         gene1,
         gene2,
         transcript1,
         transcript2,
         site1,
         site2,
         fusion,
         rr_site1,
         rr_site2,
         cc_chr_band,
         cc_t_t,
         cc_c_s,
         cc_m_t,
         cc_t_p,
         dgv_site1,
         dgv_site2
         ) = (None for i in range(18))
        cosmic_fusion_counts = 0

        annoDF[['chr1', 'chr2']] = annoDF[['chr1', 'chr2']].astype(str)
        annoDF['Cosmic_Fusion_Counts'].fillna(0, inplace=True)
        annoDF[['Cosmic_Fusion_Counts']] = annoDF[['Cosmic_Fusion_Counts']].astype(int)
        indexList = annoDF.loc[annoDF['chr1'].isin([chrom1]) &
                               annoDF['pos1'].isin([int(start1)]) &
                               annoDF['chr2'].isin([chrom2]) &
                               annoDF['pos2'].isin([int(start2)]) & 
                               annoDF['str1'].isin([str1]) &
                               annoDF['str2'].isin([str2])].index.tolist()
        if(len(indexList) > 1):
            if(verbose):
                logger.fatal(
                    "iCallSV::MergeFinalFile: More then one sv have same coordinate in same sample for annotated file. Please check and rerun")
            sys.exit(1)
        else:
            annoIndex = indexList[0]

        gene1 = annoDF.iloc[annoIndex]['gene1']
        gene2 = annoDF.iloc[annoIndex]['gene2']
        transcript1 = annoDF.iloc[annoIndex]['transcript1']
        transcript2 = annoDF.iloc[annoIndex]['transcript2']
        site1 = annoDF.iloc[annoIndex]['site1']
        site2 = annoDF.iloc[annoIndex]['site2']
        fusion = annoDF.iloc[annoIndex]['fusion']
        rr_site1 = annoDF.iloc[annoIndex]['repName-repClass-repFamily:-site1']
        rr_site2 = annoDF.iloc[annoIndex]['repName-repClass-repFamily:-site2']
        cosmic_fusion_counts = int(annoDF.iloc[annoIndex]['Cosmic_Fusion_Counts'])
        cc_chr_band = annoDF.iloc[annoIndex]['CC_Chr_Band']
        cc_t_t = annoDF.iloc[annoIndex]['CC_Tumour_Types(Somatic)']
        cc_c_s = annoDF.iloc[annoIndex]['CC_Cancer_Syndrome']
        cc_m_t = annoDF.iloc[annoIndex]['CC_Mutation_Type']
        cc_t_p = annoDF.iloc[annoIndex]['CC_Translocation_Partner']
        dgv_site1 = annoDF.iloc[annoIndex]['DGv_Name-DGv_VarType-site1']
        dgv_site2 = annoDF.iloc[annoIndex]['DGv_Name-DGv_VarType-site2']

        if(confDF is None):
            confidenceScore = None
        else:
            # Get information for confidence score
            confIndex = None
            confidenceScore = None
            confDF[['Chr1', 'Chr2']] = confDF[['Chr1', 'Chr2']].astype(str)
            indexList = confDF.loc[
                confDF['Chr1'].isin([chrom1]) & confDF['Start1'].isin(
                    [int(start1 - abs(ciPosNeg))]) & confDF['Chr2'].isin([chrom2]) &
                confDF['Start2'].isin([int(start2 - abs(ciEndNeg))])].index.tolist()

            if(len(indexList) > 1):
                if(verbose):
                    logger.fatal(
                        "iCallSV::MergeFinalFile: More then one sv have same coordinate in same sample for confidence score. Please check and rerun")
                sys.exit(1)
            else:
                confIndex = indexList[0]
            confidenceScore = np.float(confDF.iloc[confIndex]['ProbabilityScore'])

        # populate final dataframe
        outDF.loc[count,
                  ["TumorId", "NormalId", "Chr1", "Pos1", "Chr2", "Pos2", "SV_Type", "Gene1",
                   "Gene2", "Transcript1", "Transcript2", "Site1Description", "Site2Description",
                   "Fusion", "ProbabilityScore", "Confidence", "Comments", "Connection_Type",
                   "SV_LENGTH", "MAPQ", "PairEndReadSupport", "SplitReadSupport", "BrkptType",
                   "ConsensusSequence", "TumorReferenceCount", "TumorSplitReferenceCount",
                   "TumorVariantCount", "TumorSplitVariantCount", "TumorReadCount",
                   "TumorGenotypeQScore", "NormalReferenceCount", "NormalSplitReferenceCount",
                   "NormalVariantCount", "NormalSplitVariantCount", "NormalReadCount",
                   "NormalGenotypeQScore", "Cosmic_Fusion_Counts", "repName-repClass-repFamily:-site1",
                   "repName-repClass-repFamily:-site2", "CC_Chr_Band", "CC_Tumour_Types(Somatic)",
                   "CC_Cancer_Syndrome", "CC_Mutation_Type", "CC_Translocation_Partner",
                   "DGv_Name-DGv_VarType-site1", "DGv_Name-DGv_VarType-site2"]] = [aId, bId, chrom1,
                                                                                   start1, chrom2, start2, svtype, gene1, gene2, transcript1, transcript2, site1, site2,
                                                                                   fusion, confidenceScore, None, None, contype, svlengthFromDelly, mapqFromDelly,
                                                                                   peSupportFromDelly, srSupportFromDelly, brktype, conseq, caseDR, caseRR, caseDV, caseRV,
                                                                                   caseRC, caseGQ, controlDR, controlRR, controlDV, controlRV, controlRC, controlGQ, cosmic_fusion_counts,
                                                                                   rr_site1, rr_site2, cc_chr_band, cc_t_t, cc_c_s, cc_m_t, cc_t_p, dgv_site1, dgv_site2]

        count = count + 1

    # Write Output
    outFile = outDir + "/" + outputPrefix + "_merged.txt"
    outDF.to_csv(outFile, sep='\t', index=False)
    if(verbose):
        logger.info("iCallSV::MergeFinalFile: Finished merging, Final data written in %s", outFile)
    return(outFile)
