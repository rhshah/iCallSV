'''
Created on May 17, 2015
Description: Merge VCF, iAnnotateSV tab and targetSeqView tab file into a single tab-delimited file

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

@author: Ronak H Shah
'''
import sys
import logging
import vcf
import checkparameters as cp
import pandas as pd
import re


def run(aId, bId, vcfFile, annoTab, confTab, outDir, outputPrefix, verbose):
    if(verbose):
        logging.info(
            "iCallSV::MergeFinalFile: Merging Delly Filtered VCF, iAnnotateSV tab and targetSeqView tab file into a single tab-delimited file")
    cp.checkFile(vcfFile)
    cp.checkFile(annoTab)
    cp.checkFile(confTab)
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
            "TumorVariantCount",
            "TumorSplitVariantCount",
            "TumorReadCount",
            "TumorGenotypeQScore",
            "NormalVariantCount",
            "NormalSplitVariantCount",
            "NormalReadCount",
            "NormalGenotypeQScorerepName-repClass-repFamily:-site1",
            "repName-repClass-repFamily:-site2",
            "CC_Chr_Band",
            "CC_Tumour_Types(Somatic)",
            "CC_Cancer_Syndrome",
            "CC_Mutation_Type",
            "CC_Translocation_Partner",
            "DGv_Name-DGv_VarType-site1",
            "DGv_Name-DGv_VarType-site2"])
    annoDF = pd.read_csv(annoTab, sep="\t", header=0, keep_default_na='True')
    confDF = pd.read_csv(confTab, sep="\t", header=0, keep_default_na='True')
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
         svlengthFromDelly,
         mapqFromDelly,
         svtype,
         brktype,
         peSupportFromDelly,
         srSupportFromDelly,
         contype,
         conseq,
         caseRC,
         caseGQ,
         caseDR,
         caseDV,
         caseRR,
         caseRV,
         caseFT,
         controlGQ,
         controlRC,
         controlDR,
         controlDV,
         controlRR,
         controlRV,
         controlFT) = (None for i in range(27))
        chrom1 = record.CHROM
        start1 = record.POS
        filter = record.FILTER
        if("END" in record.INFO):
            start2 = record.INFO['END']
        if("CHR2" in record.INFO):
            chrom2 = record.INFO['CHR2']
        if("SVLEN" in record.INFO):
            svlengthFromDelly = record.INFO['SVLEN']
        else:
            svlengthFromDelly = abs(start2 - start1)
        if("MAPQ" in record.INFO):
            mapqFromDelly = record.INFO['MAPQ']
        if("SVTYPE" in record.INFO):
            svtype = record.INFO['SVTYPE']
        if("PE" in record.INFO):
            peSupportFromDelly = record.INFO['PE']
        if("SR" in record.INFO):
            srSupportFromDelly = record.INFO['SR']
        if("CT" in record.INFO):
            contype = record.INFO['CT']
        if("CONSENSUS" in record.INFO):
            conseq = record.INFO['CONSENSUS']
        if(record.INFO.is_sv_precise()):
            brktype = "PRECISE"
        else:
            brktype = "IMPPRECISE"

        caseCalls = record.genotype(caseIDinVcf)
        controlCalls = record.genotype(controlIDinVcf)

        if(hasattr(caseCalls.data, "FT")):
            caseFT = caseCalls.data.FT
        if(hasattr(caseCalls.data, "GQ")):
            caseGQ = caseCalls.data.GQ
        if(hasattr(caseCalls.data, "RC")):
            caseRC = caseCalls.data.RC
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
        if(hasattr(controlCalls.data, "GQ")):
            controlGQ = controlCalls.data.GQ
        if(hasattr(controlCalls.data, "RC")):
            controlRC = controlCalls.data.RC
        if(hasattr(controlCalls.data, "DR")):
            controlDR = controlCalls.data.DR
        if(hasattr(controlCalls.data, "DV")):
            controlDV = controlCalls.data.DV
        if(hasattr(controlCalls.data, "RR")):
            controlRR = controlCalls.data.RR
        if(hasattr(controlCalls.data, "RV")):
            controlRV = controlCalls.data.RV

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
        indexList = ((str(annoDF['chr1']) == str(chrom1)) &
                     (str(annoDF['pos1']) == str(start1)) &
                     (str(annoDF['chr2']) == str(chrom2)) &
                     (str(annoDF['pos2']) == str(start2))).index.tolist()
        if(len(indexList) > 1):
            if(verbose):
                logging.fatal(
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
        cc_chr_band = annoDF.iloc[annoIndex]['CC_Chr_Band']
        cc_t_t = annoDF.iloc[annoIndex]['CC_Tumour_Types(Somatic)']
        cc_c_s = annoDF.iloc[annoIndex]['CC_Cancer_Syndrome']
        cc_m_t = annoDF.iloc[annoIndex]['CC_Mutation_Type']
        cc_t_p = annoDF.iloc[annoIndex]['CC_Translocation_Partner']
        dgv_site1 = annoDF.iloc[annoIndex]['DGv_Name-DGv_VarType-site1']
        dgv_site2 = annoDF.iloc[annoIndex]['DGv_Name-DGv_VarType-site2']

        # Get information for confidence score
        confIndex = None
        confidenceScore = None
        indexList = ((str(confDF['Chr1']) == str(chrom1)) &
                     (str(confDF['Start1']) == str(start1)) &
                     (str(confDF['Chr2']) == str(chrom2)) &
                     (str(confDF['Start2']) == str(start2))).index.tolist()
        if(len(indexList) > 1):
            if(verbose):
                logging.fatal(
                    "iCallSV::MergeFinalFile: More then one sv have same coordinate in same sample for confidence score. Please check and rerun")
            sys.exit(1)
        else:
            confIndex = indexList[0]
        confidenceScore = confDF.iloc[confIndex]['ProbabilityScore']

        # populate final dataframe
        outDF.loc[count,
                  ["TumorId", "NormalId", "Chr1", "Pos1", "Chr2", "Pos2", "SV_Type", "Gene1",
                   "Gene2", "Transcript1", "Transcript2", "Site1Description", "Site2Description",
                   "Fusion", "ProbabilityScore", "Confidence", "Comments", "Connection_Type",
                   "SV_LENGTH", "MAPQ", "PairEndReadSupport", "SplitReadSupport", "BrkptType",
                   "ConsensusSequence", "TumorVariantCount", "TumorSplitVariantCount",
                   "TumorReadCount", "TumorGenotypeQScore", "NormalVariantCount",
                   "NormalSplitVariantCount", "NormalReadCount", "NormalGenotypeQScore",
                   "repName-repClass-repFamily:-site1", "repName-repClass-repFamily:-site2",
                   "CC_Chr_Band", "CC_Tumour_Types(Somatic)", "CC_Cancer_Syndrome",
                   "CC_Mutation_Type", "CC_Translocation_Partner", "DGv_Name-DGv_VarType-site1",
                   "DGv_Name-DGv_VarType-site2"]] = [aId, bId, chrom1, start1, chrom2, start2,
                                                     svtype, gene1, gene2, transcript1, transcript2, site1, site2, fusion, confidenceScore,
                                                     None, None, contype, svlengthFromDelly, mapqFromDelly, peSupportFromDelly,
                                                     srSupportFromDelly, brktype, conseq, caseDV, caseRV, caseDR, caseRC, caseGQ, controlDV,
                                                     controlRV, controlDR, controlRC, controlGQ, rr_site1, rr_site2, cc_chr_band, cc_t_t,
                                                     cc_c_s, cc_m_t, cc_t_p, dgv_site1, dgv_site2]

        count = count + 1

    # Write Output
    outFile = outDir + "/" + outputPrefix + "_final.txt"
    outDF.to_csv(outFile, sep='\t', index=False)
    if(verbose):
        logging.info("iCallSV::MergeFinalFile: Finished merging, Final data written in %s", outFile)
    return(outFile)
