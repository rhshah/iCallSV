#!/usr/local/bin/python2.7
# encoding: utf-8
"""
iCallSV
~~~~~~~

:Description: iCallSV is a wrapper to the iCallSV package which facilitates calling structural variants from Next Generation Sequencing methods such as Illumina
:author:     Ronak H Shah
:copyright:  (c) 2015-2017 by Ronak H Shah for Memorial Sloan Kettering Cancer Center. All rights reserved.
:license:    Apache License 2.0
:contact:    rons.shah@gmail.com

"""

import sys
import os
import time
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import ConfigParser as configparser
import logging

try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    print "iCallSV: coloredlogs is not installed, please install it if you wish to see color in logs on standard out."
    pass
try:
    import make_analysis_dir as mad
    import checkparameters as cp
    import launch_Run_Delly as lrd
    import launch_FilterDellyCalls as lfd
    import dellyVcf2Tab as dvcf2tab
    import combineVCF as cvcf
    import Run_iAnnotateSV as annSV
    import dellyVcf2targetSeqView as dvcf2tsv
    import Run_targetSeqView as rtsv
    import mergeFinalFiles as mff
    import filterAnnotatedSV as fas
    import helper as hp
except ImportError, e:
    print "iCallSV: sub python were not imported, please make sure that sub python scripts are in same folder as iCallSV.py."
    sys.exit(1)


__all__ = []
__version_info__ = ('0', '0', '7')
__version__ = '.'.join(__version_info__)
__date__ = '2015-03-30'
__updated__ = '2017-01-24'


def main(argv=None):  # IGNORE:C0111
    """
    Command Line Usage.

    .. code-block:: sh

    > python iCallSV.py -h

    usage: iCallSV.py [-h] [-v] [-V] -sc config.ini -abam caseBAMFile.bam -bbam
                      controlBAMFile.bam -aId caseID -bId controlID -o
                      /somepath/output -op TumorID

    iCallSV.iCallSV -- wrapper to run iCallSV package

      Created by Ronak H Shah on 2015-03-30.
      Copyright 2016-2017 Ronak H Shah. All rights reserved.

      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0

      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.

    USAGE

    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         set verbosity level [default: True]
      -V, --version         show program's version number and exit
      -sc config.ini, --svConfig config.ini
                            Full path to the structural variant configuration
      -abam caseBAMFile.bam, --caseBam caseBAMFile.bam
                            Full path to the case bam file
      -bbam controlBAMFile.bam, --controlBam controlBAMFile.bam
                            Full path to the control bam file
      -aId caseID, --caseId caseID
                            Id of the case to be analyzed, this will be the sub-
                            folder
      -bId controlID, --controlId controlID
                            Id of the control to be used, this will be used for
                            filtering variants
      -o /somepath/output, --outDir /somepath/output
                            Full Path to the output dir.
      -op TumorID, --outPrefix TumorID
                            Id of the Tumor bam file which will be used as the
                            prefix for output files


    :Example:
    .. code-block:: sh

    python iCallSV.py -sc /path/to/template.ini -abam /path/to/casebamFile -bbam /path/to/controlbamFile -aId caseID -bId controlId -o /path/to/output/directory -op prefix_for_the_output_files

    """

    # Command line options.

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by Ronak H Shah on %s.
  Copyright 2015-2016 Ronak H Shah. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    # Setup argument parser
    parser = ArgumentParser(
        description=program_license,
        formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="set verbosity level [default: %(default)s]")
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument(
        "-sc",
        "--svConfig",
        action="store",
        dest="config_file",
        required=True,
        metavar='config.ini',
        help="Full path to the structural variant configuration")
    parser.add_argument(
        "-abam",
        "--caseBam",
        action="store",
        dest="caseBam",
        required=True,
        metavar='caseBAMFile.bam',
        help="Full path to the case bam file")
    parser.add_argument(
        "-bbam",
        "--controlBam",
        action="store",
        dest="controlBam",
        required=True,
        metavar='controlBAMFile.bam',
        help="Full path to the control bam file")
    #parser.add_argument("-afastq", "--caseFastq", action="store", dest="caseFastq", required=False, metavar='caseFastqFile.fastq', help="Full path to the case fastq file")
    #parser.add_argument("-bfastq", "--controlFastq", action="store", dest="controlFastq", required=False, metavar='controlFastqFile.fastq', help="Full path to the control fastq file")
    parser.add_argument(
        "-aId",
        "--caseId",
        action="store",
        dest="caseId",
        required=True,
        metavar='caseID',
        help="Id of the case to be analyzed, this will be the sub-folder")
    parser.add_argument(
        "-bId",
        "--controlId",
        action="store",
        dest="controlId",
        required=True,
        metavar='controlID',
        help="Id of the control to be used, this will be used for filtering variants")
    #parser.add_argument("-t", "--threads", action="store", dest="threads", required=True, metavar='5', help="Number of Threads to be used to run tools")
    parser.add_argument(
        "-o",
        "--outDir",
        action="store",
        dest="outdir",
        required=True,
        metavar='/somepath/output',
        help="Full Path to the output dir.")
    parser.add_argument(
        "-op",
        "--outPrefix",
        action="store",
        dest="outprefix",
        required=True,
        metavar='TumorID',
        help="Id of the Tumor bam file which will be used as the prefix for output files")

    # Process arguments
    args = parser.parse_args()
    #
    # Parse a config ini-style file
    #
    verbose = args.verbose

    here = os.path.realpath('.')

    config_file = args.config_file
    config = configparser.ConfigParser(defaults={'here': here})
    config.read(args.config_file)
    # Make output dir
    (tag, sampleOutdirForDelly) = mad.makeOutputDir(args, "DellyDir")
    # Create Logger if verbose
    loggeroutput = sampleOutdirForDelly + "/" + args.outprefix + "_iCallSV.log"
    '''
    logger.basicConfig(
        filename=loggeroutput,
        filemode='w',
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    '''
    # create logger with 'spam_application'
    logger = logging.getLogger("iCallSV")
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(loggeroutput)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    # Print if Verbose mode is on
    if(verbose):
        logger.info("iCallSV:Verbose mode on")
    if(verbose):
        logger.info('iCallSV:Reading configuration from %s', config_file)
    runDellyTag = True
    analysisFiles = []
    analysisType = ["DEL", "DUP", "INV", "TRA"]
    if(tag == False):
        for analysis in analysisType:
            at = analysis.lower()
            outputVcf = sampleOutdirForDelly + "/" + args.caseId + "_" + at + ".vcf"
            if(os.path.isfile(outputVcf)):
                tag = True
                runDellyTag = False
                analysisFiles.append(outputVcf)
            else:
                tag = False
    # Decide based on tag what to do
    if(tag):
        if(verbose):
            logger.info(
                'iCallSV:Output of delly for %s will be written in %s',
                args.caseId,
                sampleOutdirForDelly)
        # Run Delly and get the raw calls
        if(runDellyTag):
            (del_vcf, dup_vcf, inv_vcf, tra_vcf) = lrd.launch_delly_for_different_analysis_type(
                args, config, sampleOutdirForDelly)
        else:
            del_vcf, dup_vcf, inv_vcf, tra_vcf = analysisFiles
        # Run Filter Delly and get filtered calls
        (filter_del_vcf,
         filter_dup_vcf,
         filter_inv_vcf,
         filter_tra_vcf) = lfd.launch_filterdellycalls_for_different_analysis_type(args,
                                                                                   config,
                                                                                   sampleOutdirForDelly,
                                                                                   del_vcf,
                                                                                   dup_vcf,
                                                                                   inv_vcf,
                                                                                   tra_vcf)
        #Check and Combine all VCF to a single VCF file
        listOfFilteredVCFfiles = [filter_del_vcf, filter_dup_vcf, filter_inv_vcf, filter_tra_vcf]
        listOfFilteredVCFfiles_copy = list(listOfFilteredVCFfiles)
        for vcffile in listOfFilteredVCFfiles_copy:
            if( not os.path.isfile(vcffile)):
                if(verbose):
                    logger.warning("VCF file %s does not exists.", vcffile)
                listOfFilteredVCFfiles.remove(vcffile)
        
        combinedVCF = sampleOutdirForDelly + "/" + args.caseId + "_allSVFiltered.vcf"
        combinedVCF = cvcf.run(listOfFilteredVCFfiles, combinedVCF, verbose)
        # Check if VCF file is empty
        hasRecords = False
        with open(combinedVCF, 'r') as filecontent:
            if any(not line.startswith("#") for line in filecontent):
                hasRecords = True
            else:
                hasRecords = False
        # If there are VCF records do further analysis.
        if(hasRecords):
            combinedAnnPrefix = args.caseId + "_allSVFiltered"
            combinedTargetSeqView = args.caseId + "_allSVFiltered_tsvInput.txt"
            combinedTargetSeqViewCscore = args.caseId + "_allSVFiltered_cScore.txt"
            # convert vcf files to tab-delimited using vcf2tab
            combinedTAB = dvcf2tab.vcf2tab(combinedVCF, sampleOutdirForDelly, verbose)
            # Annotate using iAnnotateSV
            combinedAnnTAB = annSV.run(
                config.get("Python", "PYTHON"),
                config.get("iAnnotateSV", "ANNOSV"),
                config.get("iAnnotateSV", "GENOMEBUILD"),
                int(config.get("iAnnotateSV", "DISTANCE")),
                config.get("iAnnotateSV", "CANONICALTRANSCRIPTFILE"),
                config.get("iAnnotateSV", "UNIPROTFILE"),
                config.get("iAnnotateSV", "CosmicCensus"),
                config.get("iAnnotateSV", "CosmicFusionCounts"),
                config.get("iAnnotateSV", "RepeatRegionAnnotation"),
                config.get("iAnnotateSV", "DGvAnnotations"),
                combinedTAB, combinedAnnPrefix, sampleOutdirForDelly)
            # convert vcf to targetseqviewformat
            combinedTargetSeqView = dvcf2tsv.Convert2targetSeqView(
                args.caseId,
                os.path.basename(os.path.abspath(args.caseBam)),
                os.path.basename(os.path.abspath(args.caseBam)),
                combinedVCF,
                sampleOutdirForDelly,
                combinedTargetSeqView)
            # Get Confidence score using targetSeqView
            combinedTargetSeqViewCscore = rtsv.run(
                config.get(
                    "R", "RHOME"), config.get(
                    "TargetSeqView", "CalculateConfidenceScore"), 5, os.path.dirname(
                    os.path.abspath(
                        args.caseBam)), combinedTargetSeqView, config.get(
                        "TargetSeqView", "GENOMEBUILD"), int(
                            config.get(
                                "TargetSeqView", "ReadLength")), sampleOutdirForDelly, combinedTargetSeqViewCscore)
            # Merge Results from vcf, tab and targetseqview
            finalFile = mff.run(
                args.caseId,
                args.controlId,
                combinedVCF,
                combinedAnnTAB,
                combinedTargetSeqViewCscore,
                sampleOutdirForDelly,
                args.outprefix,
                args.verbose)
            # Filter Final Results
            finalFilterResults = fas.run(
                finalFile, sampleOutdirForDelly, args.outprefix, config.get(
                    "BlackListRegions", "BlackListGenes"), args.verbose, config.get(
                    "HotSpotRegions", "GenesToKeep"))

        else:
            if(verbose):
                logger.warn(
                    "All Records have been filtered in standard filtered step. Thus we will exit the program and not proceed.")
            outputFile = os.path.join(sampleOutdirForDelly, args.outprefix + "_final.txt")
            hp.make_empty_outputfile(outputFile)
            if(verbose):
                logger.info("Thank you for using iCallSV.")
            sys.exit(0)
    else:
        if(verbose):
            logger.fatal(
                "The output directory for the %s already exists. Please delete %s folder and rerun",
                args.caseId,
                sampleOutdirForDelly)
        sys.exit(1)


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("iCallSV:Elapsed time was %g seconds", totaltime)
