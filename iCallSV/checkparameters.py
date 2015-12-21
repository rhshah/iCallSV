'''
Created on Mar 16, 2015
Description: This modules checks the parameters for various type of inputs.
@author: Ronak H Shah
'''

import os
import sys
import logging


# Check if the file exist

def checkFile(fileToCheck):
    if(os.path.isfile(fileToCheck)):
        logging.info("checkparameters:Given File: %s exists.", fileToCheck)
    else:
        logging.fatal(
            "checkparameters:Given File: %s does not exists. Sorry please check this Input and Run Again.",
            fileToCheck)
        sys.exit(1)


# Check if the Directory exists

def checkDir(folderToCheck):
    if(os.path.isdir(folderToCheck)):
        logging.info("checkparameters:Given Directory: %s exists.", folderToCheck)
    else:
        logging.fatal(
            "checkparameters:Given Directory: %s does not exists. Sorry please check this Input and Run Again.",
            folderToCheck)
        sys.exit(1)


# Check if the variable is and Integer

def checkInt(variableToCheck, variableName):
    if(isinstance(variableToCheck, int)):
        logging.info(
            "checkparameters: %s Variable: %d is an Integer.",
            variableName,
            variableToCheck)
    else:
        logging.fatal(
            "checkparameters: %s Variable: %d is not an Integer. Sorry please check this Input and Run Again.",
            variableName,
            variableToCheck)
        sys.exit(1)


# Check if the given variable is not Empty

def checkEmpty(variableToCheck, variableName):
    if(variableToCheck):
        logging.info("checkparameters: %s Variable:%s is not empty.", variableName, variableToCheck)
    else:
        logging.fatal(
            "checkparameters: %s Variable:%s is empty. Sorry please check this Input and Run Again.",
            variableName,
            variableToCheck)
        sys.exit(1)


# Check the Delly Analysis Type is Valid or Not

def checkDellyAnalysisType(varaibleToCheck):
    if(varaibleToCheck == "DEL" or varaibleToCheck == "DUP" or varaibleToCheck == "INV" or varaibleToCheck == "TRA"):
        logging.info("checkparameters:Given Delly Analysis Type:%s is valid.", varaibleToCheck)
    else:
        logging.fatal(
            "checkparameters:Given Delly Analysis Type: %s is not valid. Sorry please check this Input and Run Again.",
            varaibleToCheck)
        sys.exit(1)
