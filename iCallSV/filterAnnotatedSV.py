"""
filterAnnotatedSV
~~~~~~~~~~~~~~~~~

:Description: This module will filter calls from the merged file
"""

'''
Created on Mar 17, 2015
Description: This module will filter calls from the merged file
@author: Ronak H Shah

::Inputs::
inputTxt: Filter Text File
outputDir: Output directory
outPrefix: Prefix of the output file
blacklistGenesFile: List of genes that should be eliminated
genesToKeepFile: List of genes that should be kept
verbose: Mode

::Output::
Filtered Output files
'''
import os
import pandas as pd
import logging
import checkparameters as cp
import re

# Initiate logger
logger = logging.getLogger('iCallSV.FilterDellyCalls')


def run(inputTxt, outputDir, outPrefix, blacklistGenesFile, genesToKeepFile=None, verbose):
    """
    This will ``filter sv calls`` from the final merged file.

    :param str inputTxt: str for the txt file to be filtered
    :param str outputDir: str for the output directory
    :param str outputPrefix: str prefix for the output File
    :param str blacklistGenesFile: str for the txt file containing blacklisted genes
    :param str genesToKeepFile: str for the txt file containing genes to keep
    :param bool verbose: a boolean
    :return: A str name of final sv file
    :rtype: str

    """
    cp.checkFile(inputTxt)
    cp.checkFile(blacklistGenesFile)
    cp.checkDir(outputDir)
    cp.checkEmpty(outPrefix, "Prefix for the output file")
    if(os.path.isfile(genesToKeepFile)):
        keepGenes = [line.strip() for line in open(genesToKeepFile, 'r')]
    else:
        keepGenes = None
    inputDF = pd.read_csv(inputTxt, sep="\t", header=0, keep_default_na='True')
    outputDF = pd.DataFrame(columns=inputDF.columns)
    outputFile = os.path.join(outputDir, outPrefix + "_final.txt")
    count = 0
    for index, row in inputDF.iterrows():
        gene1 = row.loc['Gene1']
        gene2 = row.loc['Gene2']
        transcript1 = row.loc['Transcript1']
        transcript2 = row.loc['Transcript2']
        site1 = row.loc['Site1Description']
        site2 = row.loc['Site2Description']
        fusion = row.loc['Fusion']
        # skip IGR records
        if("IGR" in site1 and "IGR" in site2):
            igrFlag = True
        else:
            igrFlag = False

        # check records from these gene
        if(keepGenes):
            keepGeneFlag = checkBlackListGene(gene1, gene2, keepGenes)
        else:
            keepGeneFlag = True
        # check records from these gene
        blacklistGenes = [line.strip() for line in open(blacklistGenesFile, 'r')]
        blacklistGeneFlag = checkBlackListGene(gene1, gene2, blacklistGenes)
        # skip record occurring within intron
        eventInIntronFlag = False
        if((gene1 == gene2) and ((igrFlag is False) or (blacklistGeneFlag is False)) and ("Intron" in site1 and "Intron" in site2)):
            eventInIntronFlag = checkEventInIntronFlag(gene1, gene2, site1, site2)
        else:
            pass

        if((keepGeneFlag is False) or (igrFlag) or (blacklistGeneFlag) or (eventInIntronFlag)):
            if(verbose):
                logger.warn(
                    "iCallSV::FilterFinalFile: Record will be Filtered as keepGeneFlag:%s, IGR:%s, blackListGene:%s, Intronic Event:%s",
                    keepGeneFlag,
                    igrFlag,
                    blacklistGeneFlag,
                    eventInIntronFlag)
            continue
        else:
            outputDF.loc[count] = row
            count = count + 1

    outputDF.to_csv(outputFile, sep='\t', index=False)
    if(verbose):
        logger.info(
            "iCallSV::FilterFinalFile: Finished Filtering, Final data written in %s",
            outputFile)

    return(outputFile)

# Check if the gene is a blacklist gene


def checkGeneListToKeep(gene1, gene2, keepGenes):
    if((gene1 in keepGenes) or (gene2 in keepGenes)):
        kgFlag = True
    else:
        kgFlag = False
    return(kgFlag)

# Check if the gene is a blacklist gene


def checkBlackListGene(gene1, gene2, blacklistGenes):
    """
    This will ``check for blacklisted genes``

    :param str gene1: str for the name of gene at breakpoint 1
    :param str gene2: str for the name of gene at breakpoint 2
    :param list blacklistGenes:  list containing blacklisted genes
    :param str genesToKeepFile: str for the txt file containing genes to keep
    :return: A boolean tag indicating True or False
    :rtype: bool

    """
    if((gene1 in blacklistGenes) or (gene2 in blacklistGenes)):
        bgFlag = True
    else:
        bgFlag = False
    return(bgFlag)


# Check if the event is in the intron only and not affecting splicing
def checkEventInIntronFlag(gene1, gene2, site1, site2):
    """
    This will ``Check if the event is in the intron only and not affecting
    splicing``


    :param str gene1: str for the name of gene at breakpoint 1
    :param str gene2: str for the name of gene at breakpoint 2
    :param str site1: str for the description of site in breakpoint 1
    :param str site2: str for the description of site in breakpoint 2
    :return: A boolean tag indicating True or False
    :rtype: bool

    """
    if(gene1 == gene2):
        (s1A, s1B) = site1.split(":")
        (s2A, s2B) = site2.split(":")
        (s1a, s1b, s1c, s1d) = s1B.split(" ")
        (s2a, s2b, s2c, s2d) = s2B.split(" ")
        if(("before" in site1 and "before" in site2) or ("after" in site1 and "after" in site2)):
            if(int(s1d) == int(s2d)):
                if("bp" in s1a):
                    s1location = re.findall(r'\d+', s1a)[0]
                    s2location = re.findall(r'\d+', s2a)[0]
                    if(int(s1location) < 5 or int(s2location) < 5):
                        eviFlag = False
                    else:
                        eviFlag = True
                else:
                    eviFlag = False
            else:
                eviFlag = False
        else:
            eviFlag = False
    else:
        eviFlag = False
    return(eviFlag)
