"""
Created on Mar 17, 2015
Description: This module will read the hotspot file and tell if it is a hotspot or not
@author: Ronak H Shah
"""
"""
::Inputs::
HotSpotFile: List of Genes that have Hotspot Structural Variants (Tab-delimited Format without header:chr    start    end    geneName).
chr1: Chromosome location for 1st breakpoint
start1: Start location of the breakpoint
chr2: Chromosome location for 2nd breakpoint
start2: Start Location of the second breakpoint
"""
import os
import re


# Read the hotspot file and make a dictionary of it
def ReadHotSpotFile(HotSpotFile):
    hotspotDict = {}
    if os.path.isfile(HotSpotFile):
        with open(HotSpotFile, 'r') as filecontent:
            for line in filecontent:
                (chrom, start, end, gene) = line.rstrip().split("\t")
                if chrom in hotspotDict:
                    val = hotspotDict[chrom]
                    hotspotDict[chrom] = val + "#" + start + ":" + end + ":" + gene
                else:
                    hotspotDict[chrom] = start + ":" + end + ":" + gene
    return(hotspotDict)


# Read the hotspot dictionary and tell if the event occurs in hotspot or not
def CheckIfItIsHotspot(chr1, start1, chr2, start2, hotspotDict):
    hotspotTag = False
    chr1 = str(chr1)
    chr2 = str(chr2)
    start1 = int(start1)
    start2 = int(start2)
    if chr1 == chr2:
        start = start1
        end = start2
        if chr1 in hotspotDict:
            allCordinates = hotspotDict[chr1]
            pattern = re.compile("#")
            match = re.search(pattern, allCordinates)
            if(match):
                coordList = re.split(pattern, allCordinates)
                for item in coordList:
                    pattern = re.compile(":")
                    (hsStart, hsEnd, gene) = re.split(pattern, item)
                    if((start >= int(hsStart) and start <= int(hsEnd)) or (end >= int(hsStart) and end <= int(hsEnd))):
                        hotspotTag = True
                        return(hotspotTag)
                    else:
                        hotspotTag = False
            else:
                pattern = re.compile(":")
                (hsStart, hsEnd, gene) = re.split(pattern, allCordinates)
                if((start >= int(hsStart) and start <= int(hsEnd)) or (end >= int(hsStart) and end <= int(hsEnd))):
                    hotspotTag = True
                    return(hotspotTag)
                else:
                    hotspotTag = False
        else:
            hotspotTag = False
    else:
        if chr1 in hotspotDict:
            allCordinates = hotspotDict[chr1]
            pattern = re.compile("#")
            match = re.search(pattern, allCordinates)
            if(match):
                coordList = re.split(pattern, allCordinates)
                for item in coordList:
                    pattern = re.compile(":")
                    (hsStart, hsEnd, gene) = re.split(pattern, item)
                    if(start1 >= int(hsStart) and start1 <= int(hsEnd)):
                        hotspotTag = True
                        return(hotspotTag)
                    else:
                        hotspotTag = False
            else:
                pattern = re.compile(":")
                (hsStart, hsEnd, gene) = re.split(pattern, allCordinates)
                if(start1 >= int(hsStart) and start1 <= int(hsEnd)):
                    hotspotTag = True
                    return(hotspotTag)
                else:
                    hotspotTag = False
        if ((chr2 in hotspotDict) and (hotspotTag is False)):
            allCordinates = hotspotDict[chr2]
            pattern = re.compile("#")
            match = re.search(pattern, allCordinates)
            if(match):
                coordList = re.split(pattern, allCordinates)
                for item in coordList:
                    pattern = re.compile(":")
                    (hsStart, hsEnd, gene) = re.split(pattern, item)
                    if(start2 >= int(hsStart) and start2 <= int(hsEnd)):
                        hotspotTag = True
                        return(hotspotTag)
                    else:
                        hotspotTag = False
            else:
                pattern = re.compile(":")
                (hsStart, hsEnd, gene) = re.split(pattern, allCordinates)
                if(start2 >= int(hsStart) and start2 <= int(hsEnd)):
                    hotspotTag = True
                    return(hotspotTag)
                else:
                    hotspotTag = False
    return(hotspotTag)
    
# # Test module
# hotspotDict = ReadHotSpotFile(
#     "/home/shahr2/workspace/dmp-data/mskdata/interval-lists/structuralvariants_geneInterval.txt")
# hotspotTag = CheckIfItIsHotspot(2, 29417089, 2, 30143125, hotspotDict)
# print hotspotTag
