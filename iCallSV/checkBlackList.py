'''
Created on Nov 20, 2015
Description: This module will read the Black List file and tell if and event is blacklisted or not
@author: Ronak H Shah
'''
'''
::Inputs::
BlackListFile: List of Position that have Black List Structural Variants (Tab-delimited Format without header:chr1    start1    end1	chr2    start2    end2).
chr1: Chromosome location for 1st breakpoint
start1: Start location of the 1st breakpoint
chr2: Chromosome location for 2nd breakpoint
start2: Start Location of the second breakpoint

'''
import os
import re


# Read the hotspot file and make a dictionary of it
def ReadBlackListFile(BlackListFile):
    blacklist = ()
    if os.path.isfile(BlackListFile):
        with open(BlackListFile, 'r') as filecontent:
            for line in filecontent:
                #(chrom1, start1, chrom2, start2) = line.rstrip().split("\t")
                blacklist.append(line)
    return(blacklist)


# Read the hotspot dictionary and tell if the event occurs in blacklist or not
def CheckIfItIsBlacklisted(chr1, start1, chr2, start2, blacklist, range):
    blacklistTag = None
    chr1 = str(chr1)
    chr2 = str(chr2)
    start1 = int(start1)
    start2 = int(start2)
    if chr1 == chr2:
        start = start1
        end = start2
        for entries in blacklist:
            (bchr1, bstart1, bchr2, bstart2) = entries.rstrip().split("\t")
            if(bchr1 == bchr2):
                if(chr1 == bchr1):
                    if(((start <= int(bstart1 - range)) and (start >= int(bstart1 + range))) 
                       and ((end <= int(bstart2 - range)) and (end <= int(bstart2 + range)))):
                        blacklistTag=True
                    else:
                        blacklistTag=False
                    if(((end <= int(bstart1 - range)) and (end >= int(bstart1 + range))) 
                        and ((start <= int(bstart2 - range)) and (start <= int(bstart2 + range)))):
                        blacklistTag=True
        			else:
        				blacklistTag=False
                else:
                    continue
            else:
                continue

    else:
         for entries in blacklist:
             (bchr1, bstart1, bchr2, bstart2)=entries.rstrip().split("\t")
        		if(((chr1 == bchr1) and (chr2 == bchr2)) or ((chr1 == bchr2) and (chr2 == bchr1))):
        			if(((start <= int(bstart1 - range)) and (start >= int(bstart1 + range))) and ((end <= int(bstart2 - range)) and (end <= int(bstart2 + range)))):
        					blacklistTag=True
        			else:
        					blacklistTag=False
        			if(((end <= int(bstart1 - range)) and (end >= int(bstart1 + range))) and ((start <= int(bstart2 - range)) and (start <= int(bstart2 + range)))):
        					blacklistTag=True
        			else:
        					blacklistTag=False
                else:
                    continue

    return(blacklistTag)
