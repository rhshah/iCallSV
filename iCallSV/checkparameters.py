"""
Created on Mar 16, 2015
Description: This modules checks the parameters for various type of inputs.
@author: Ronak H Shah
"""

import os
import sys
import logging


# Check if the file exist

def checkFile(fileToCheck):
    logger = logging.getLogger(__name__)
    if(os.path.isfile(fileToCheck)):
        logger.info("checkparameters:Given File: %s exists.", fileToCheck)
    else:
        logger.fatal(
            "checkparameters:Given File: %s does not exists. Sorry please check this Input and Run Again.",
            fileToCheck)
        sys.exit(1)


# Check if the Directory exists

def checkDir(folderToCheck):
    logger = logging.getLogger(__name__)
    if(os.path.isdir(folderToCheck)):
        logger.info("checkparameters:Given Directory: %s exists.", folderToCheck)
    else:
        logger.fatal(
            "checkparameters:Given Directory: %s does not exists. Sorry please check this Input and Run Again.",
            folderToCheck)
        sys.exit(1)


# Check if the variable is and Integer

def checkInt(variableToCheck, variableName):
    logger = logging.getLogger(__name__)
    if(isinstance(variableToCheck, int)):
        logger.info(
            "checkparameters: %s Variable: %d is an Integer.",
            variableName,
            variableToCheck)
    else:
        logger.fatal(
            "checkparameters: %s Variable: %d is not an Integer. Sorry please check this Input and Run Again.",
            variableName,
            variableToCheck)
        sys.exit(1)


# Check if the given variable is not Empty

def checkEmpty(variableToCheck, variableName):
    logger = logging.getLogger(__name__)
    if(variableToCheck):
        logger.info("checkparameters: %s Variable:%s is not empty.", variableName, variableToCheck)
    else:
        logger.fatal(
            "checkparameters: %s Variable:%s is empty. Sorry please check this Input and Run Again.",
            variableName,
            variableToCheck)
        sys.exit(1)


# Check the Delly Analysis Type is Valid or Not

def checkDellyAnalysisType(varaibleToCheck):
    logger = logging.getLogger(__name__)
    if(varaibleToCheck == "DEL" or varaibleToCheck == "DUP" or varaibleToCheck == "INV" or varaibleToCheck == "TRA"):
        logger.info("checkparameters:Given Delly Analysis Type:%s is valid.", varaibleToCheck)
    else:
        logger.fatal(
            "checkparameters:Given Delly Analysis Type: %s is not valid. Sorry please check this Input and Run Again.",
            varaibleToCheck)
        sys.exit(1)
