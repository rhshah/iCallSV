"""
checkparameters
~~~~~~~~~~~~~~~

:Description: This modules checks the parameters for various type of inputs.

"""
'''
Created on Mar 16, 2015
Description: This modules checks the parameters for various type of inputs.
@author: Ronak H Shah
'''

import os
import sys
import logging
import coloredlogs

# initiate Logger
logger = logging.getLogger('iCallSV.checkparameters')
coloredlogs.install(level='DEBUG')
# Check if the file exist


def checkFile(fileToCheck):
    """
    Check `if the file exists or not``

    :param str fileToCheck: Name of the file to be checked.
    :return: None
    :rtype: None

    """
    if(os.path.isfile(fileToCheck)):
        logger.info("checkparameters:Given File: %s exists.", fileToCheck)
    else:
        logger.fatal(
            "checkparameters:Given File: %s does not exists. Sorry please check this Input and Run Again.",
            fileToCheck)
        sys.exit(1)


# Check if the Directory exists

def checkDir(folderToCheck):
    """
    Check `if the folder exists or not``

    :class:`str`.

    :param str folderToCheck: Name of the folder to be checked.
    :return: None
    :rtype: None

    """
    if(os.path.isdir(folderToCheck)):
        logger.info("checkparameters:Given Directory: %s exists.", folderToCheck)
    else:
        logger.fatal(
            "checkparameters:Given Directory: %s does not exists. Sorry please check this Input and Run Again.",
            folderToCheck)
        sys.exit(1)


# Check if the variable is and Integer

def checkInt(variableToCheck, variableName):
    """
    Check `if the variable is int or not``

    :param int variableToCheck: Check if it is int or not
    :param str variableName: Name of the int object to be verified
    :return: None
    :rtype: None

    """
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
    """
    Check `if the variable is None or not``

    :param str variableToCheck: check if str is None or not
    :param str variableName: Name of the None object to be verified
    :return: None
    :rtype: None

    """
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
    """
    Check `if the variable for Delly analysis exists or not``

    :param str variableToCheck: check if str is DEL|DUP|INV|TRA
    :return: None
    :rtype: None

    """
    if(varaibleToCheck == "DEL" or varaibleToCheck == "DUP" or varaibleToCheck == "INV" or varaibleToCheck == "TRA"):
        logger.info("checkparameters:Given Delly Analysis Type:%s is valid.", varaibleToCheck)
    else:
        logger.fatal(
            "checkparameters:Given Delly Analysis Type: %s is not valid. Sorry please check this Input and Run Again.",
            varaibleToCheck)
        sys.exit(1)
