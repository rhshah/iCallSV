'''
Created on Dec 23, 2015
Description: Merge VCF, iAnnotateSV tab and targetSeqView tab file into a single tab-delimited file

::Input::
sampleName: Name of the sample that has the structural abberations
sampleBamName: Name of the bam file.
sampleSplitBaName: Name of the split bam file (Use bam file if you dont have split bam file)
vcfFile: Input Delly VCF file for the conversion
outputDir: Directory to write the output file
outputFileName: Name of the output File
::Output::
outputFile: File with following header
"TumorId\tNormalId\tChr1\tPos1\tChr2\tPos2\tSV_Type\tGene1\tGene2\tTranscript1\tTranscript2\tSite1Description\tSite2Description\tFusion\tConfidence\tComments\tConnection_Type\tSV_LENGTH\tMAPQ\tPairEndReadSupport\tSplitReadSupport\tBrkptType\tConsensusSequence\tTumorVariantCount\tTumorSplitVariantCount\tTumorReadCount\tTumorGenotypeQScore\tNormalVariantCount\tNormalSplitVariantCount\tNormalReadCount\tNormalGenotypeQScore\n";

@author: Ronak H Shah
'''