#!/usr/local/bin/python2.7
# encoding: utf-8
'''
iCallSV.iCallSV -- wrapper to run iCallSV package

iCallSV.iCallSV is a wrapper to the iCallSV package which facilitates calling structural variants from Next Generation Sequencing methods such as Illumina

It defines classes_and_methods

@author:     Ronak H Shah

@copyright:  2015-2016 Ronak H Shah. All rights reserved.

@license:    Apache License 2.0

@contact:    rons.shah@gmail.com
@deffield    updated: Updated
'''

import sys
import os
import time
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import ConfigParser as configparser
import logging
import make_analysis_dir as mad
import launch_Run_Delly as lrd
import dellyVcf2Tab as dvcf2tab
import combineVCF as cvcf
import Run_iAnnotateSV as annSV
import dellyVcf2targetSeqView as dvcf2tsv
import Run_targetSeqView as rtsv

__all__ = []
__version__ = 0.1
__date__ = '2015-03-30'
__updated__ = '2015-12-20'


def main(argv=None):  # IGNORE:C0111
    """Command line options."""

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
                        default=True, help="set verbosity level [default: %(default)s]")
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
        "-caseId",
        "--caseId",
        action="store",
        dest="caseId",
        required=True,
        metavar='caseID',
        help="Id of the case to be analyzed, this will be the sub-folder")
     parser.add_argument(
        "-controlId",
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
    # Create Logger if verbose
    loggeroutput = args.outdir + "/" + args.outprefix + "_iCallSV.log"
    logging.basicConfig(
        filename=loggeroutput,
        filemode='w',
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
    # Print if Verbose mode is on
    if(verbose):
        logging.info("iCallSV:Verbose mode on")
    here = os.path.realpath('.')

    config_file = args.config_file
    if(verbose):
        logging.info('iCallSV:Reading configuration from %s', config_file)
    config = configparser.ConfigParser(defaults={'here': here})
    config.read(args.config_file)
    (tag, sampleOutdirForDelly) = mad.makeOutputDir(args, "DellyDir")
    # Decide based on tag what to do
    if(tag):
        if(verbose):
            logging.info(
                'iCallSV:Output of delly for %s will be written in %s',
                args.caseId,
                sampleOutdirForDelly)
        # Run Delly and get the raw calls
        (del_vcf, dup_vcf, inv_vcf, tra_vcf) = lrd.launch_delly_for_different_analysis_type(
            args, config, sampleOutdirForDelly)
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
        # Combine all VCF to a single VCF file
        listOfFilteredVCFfiles = [filter_del_vcf, filter_dup_vcf, filter_inv_vcf, filter_tra_vcf]
        combinedVCF = sampleOutdirForDelly + "/" + args.caseId + "_allSVFiltered.vcf"
        combinedAnnVCF = args.caseId + "_allAnnotatedSVFiltered.tab"
        combinedTargetSeqView = args.caseId + "_allSVFiltered_tsvInput.txt"
        combinedTargetSeqViewCscore = args.caseId + "_allSVFiltered_cScore.txt"
        combinedVCF = cvcf.run(listOfFilteredVCFfiles, combinedVCF, verbose)
        # convert vcf files to tab-delimited using vcf2tab
        combinedTAB = dvcf2tab.vcf2tab(combinedVCF, sampleOutdirForDelly, verbose)
        # Annotate using iAnnotateSV
        combinedAnnTAB = annSV.run(
            config.get(
                "Python", "PYTHON"), config.get(
                "iAnnotateSV", "ANNOSV"), config.get(
                "iAnnotateSV", "GENOMEBUILD"), int(config.get(
                    "iAnnotateSV", "DISTANCE")), config.get(
                        "iAnnotateSV", "CANONICALTRANSCRIPTFILE"), combinedTab, combinedAnnVCF, sampleOutdirForDelly)
        # convert vcf to targetseqviewformat
        combinedTargetSeqView = dvcf2tsv.Convert2targetSeqView(
            args.caseId,
            args.caseBam,
            args.caseBam,
            combinedVCF,
            sampleOutdirForDelly,
            combinedTargetSeqView)
        # Get Confidence score using targetSeqView
        combinedTargetSeqViewCscore = rtsv.run(config.get("R", "RHOME"), 
                                               config.get("TargetSeqView", "CalculateConfidenceScore"), 
                                               5, args.caseBam, combinedTargetSeqView, 
                                               config.get("TargetSeqView", "GENOMEBUILD"), 
                                               int(config.get("TargetSeqView", "ReadLength")), 
                                               sampleOutdirForDelly, combinedTargetSeqViewCscore)
        # Merge Results from vcf, tab and targetseqview
    else:
        if(verbose):
            logging.fatal(
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
