'''
Created on Mar 16, 2015
Description: This modules checks the parameters for various type of inputs.
@author: Ronak H Shah
'''

import os
import sys

# Check if the file exists


def checkFile(fileToCheck):
    if(os.path.isfile(fileToCheck)):
        print "Given File: ", fileToCheck, "exists.\n"
    else:
        print "Given File: ", fileToCheck, "does not exists. Sorry please check this Input and Run Again.\n"
        sys.exit()
# Check if the Dierectory exists


def checkDir(folderToCheck):
    if(os.path.isdir(folderToCheck)):
        print "Given Directory: ", folderToCheck, "exists.\n"
    else:
        print "Given Directory: ", folderToCheck, "does not exists. Sorry please check this Input and Run Again.\n"
        sys.exit()
# Check if the variable is and Integer


def checkInt(variableToCheck, variableName):
    if(isinstance(variableToCheck, int)):
        print variableName, " Variable:", variableToCheck, " is an Integer.\n"
    else:
        print variableName, " Variable:", variableToCheck, " is not an Integer. Sorry please check this Input and Run Again.\n"
        sys.exit()
# Check if the given variable is not Empty


def checkEmpty(variableToCheck, variableName):
    if(variableToCheck):
        print variableName, " Variable:", variableToCheck, " is not empty.\n"
    else:
        print variableName, " Variable:", variableToCheck, " is empty. Sorry please check this Input and Run Again.\n"
        sys.exit()
# Check the Delly Analysis Type is Valid or Not


def checkDellyAnalysisType(varaibleToCheck):
    if(varaibleToCheck == "DEL" or varaibleToCheck == "DUP" or varaibleToCheck == "INV" or varaibleToCheck == "TRA"):
        print "Given Delly Type:", varaibleToCheck, " is valid.\n"
    else:
        print "Given Delly Type:", varaibleToCheck, " is not valid. Sorry please check this Input and Run Again.\n"
        sys.exit()
