"""
checkBlackList
~~~~~~~~~~~~~~

:Description: This module will read the Black List file and tell if and event is blacklisted or not

"""
'''
Created on Nov 20, 2015
@author: Ronak H Shah
::Inputs::
BlackListFile: List of Position that have Black List Structural Variants (Tab-delimited Format without header:chr1    start1    end1	chr2    start2    end2).
chr1: Chromosome location for 1st breakpoint
start1: Start location of the 1st breakpoint
chr2: Chromosome location for 2nd breakpoint
start2: Start Location of the second breakpoint
extention: How much should the intervals be extended in positive and negative directions
'''
import os


# Read the hotspot file and make a dictionary of it
def ReadBlackListFile(BlackListFile):
    """
    Read the ``blacklist region file``

    :param str BlackListFile: str of file to be read.
    :return: A list containing black listed regions.
    :rtype: list.

    """
    blacklist = []
    if os.path.isfile(BlackListFile):
        with open(BlackListFile, 'r') as filecontent:
            for line in filecontent:
                #(chrom1, start1, chrom2, start2) = line.rstrip().split("\t")
                blacklist.append(line)
    return(blacklist)


# Read the hotspot dictionary and tell if the event occurs in blacklist or not
def CheckIfItIsBlacklisted(chr1, start1, chr2, start2, blacklist, extention):
    """
    Check if coordinate are present in the ``blacklist region file``

    :param str chr1: str of the breakpoint in first chromosome
    :param int start1: int of the start location of the breakpoint in first chromosome
    :param str chr2: str of the breakpoint in second chromosome
    :param int start2: int of the start location of the breakpoint in second chromosome
    :param list blacklist: A list containing black listed regions
    :param int extension: an value for search in positive and negative direction of the start1 and start2 location
    :return: A boolean tag indicating True or False
    :rtype: bool

    """
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
            bstart1 = int(bstart1)
            bstart2 = int(bstart2)
            if(bchr1 == bchr2):
                if(chr1 == bchr1):
                    if(((start <= int(bstart1 - extention)) and (start >= int(bstart1 + extention)))
                       and ((end <= int(bstart2 - extention)) and (end <= int(bstart2 + extention)))):
                        blacklistTag = True
                    else:
                        blacklistTag = False
                    if(((end <= int(bstart1 - extention)) and (end >= int(bstart1 + extention)))
                            and ((start <= int(bstart2 - extention)) and (start <= int(bstart2 + extention)))):
                        blacklistTag = True
                    else:
                        blacklistTag = False
                else:
                    continue
            else:
                continue

    else:
        for entries in blacklist:
            (bchr1, bstart1, bchr2, bstart2) = entries.rstrip().split("\t")
            bstart1 = int(bstart1)
            bstart2 = int(bstart2)
            if(((chr1 == bchr1) and (chr2 == bchr2)) or ((chr1 == bchr2) and (chr2 == bchr1))):
                if(((start1 <= int(bstart1 - extention)) and (start1 >= int(bstart1 + extention)))
                   and ((start2 <= int(bstart2 - extention)) and (start2 <= int(bstart2 + extention)))):
                    blacklistTag = True
                else:
                    blacklistTag = False
                if(((start2 <= int(bstart1 - extention)) and (start2 >= int(bstart1 + extention)))
                   and ((start1 <= int(bstart2 - extention)) and (start1 <= int(bstart2 + extention)))):
                    blacklistTag = True
                else:
                    blacklistTag = False
            else:
                continue

    return(blacklistTag)
