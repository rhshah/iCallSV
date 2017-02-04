"""
helper
~~~~~~~~~~~~~~~

:Description: helper has many utilities for iCallSV
"""
'''
Created on February 2, 2017
Description: helper has many utilities for iCallSV
@author: Ronak H Shah

A) make_empty_output: Will make an empty output file with header
::Input::
outFile: string containing path to write the output
::Output::
outputFile: File with following header
"TumorId\tNormalId\tChr1\tPos1\tChr2\tPos2\tSV_Type\tGene1\tGene2\tTranscript1\tTranscript2\tSite1Description\tSite2Description\tFusion\tProbabilityScore\tConfidence\tComments\tConnection_Type\tSV_LENGTH\tMAPQ\tPairEndReadSupport\tSplitReadSupport\tBrkptType\tConsensusSequence\tTumorVariantCount\tTumorSplitVariantCount\tTumorReadCount\tTumorGenotypeQScore\tNormalVariantCount\tNormalSplitVariantCount\tNormalReadCount\tNormalGenotypeQScorerepName-repClass-repFamily:-site1\trepName-repClass-repFamily:-site2\tCC_Chr_Band\tCC_Tumour_Types(Somatic)\tCC_Cancer_Syndrome\tCC_Mutation_Type\tCC_Translocation_Partner\tDGv_Name-DGv_VarType-site1\tDGv_Name-DGv_VarType-site2\n";
'''

import os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp() #So that matplotlib doesnot complain stale file handle
try:
    import pandas as pd
except ImportError, e:
    print "helper: pandas is not installed, please install pandas as it is required to run the mapping."
    sys.exit(1)

def make_empty_outputfile(outFile):
	outDF = pd.DataFrame(
			columns=[
				"TumorId",
				"NormalId",
				"Chr1",
				"Pos1",
				"Chr2",
				"Pos2",
				"SV_Type",
				"Gene1",
				"Gene2",
				"Transcript1",
				"Transcript2",
				"Site1Description",
				"Site2Description",
				"Fusion",
				"ProbabilityScore",
				"Confidence",
				"Comments",
				"Connection_Type",
				"SV_LENGTH",
				"MAPQ",
				"PairEndReadSupport",
				"SplitReadSupport",
				"BrkptType",
				"ConsensusSequence",
				"TumorReferenceCount",
				"TumorSplitReferenceCount",
				"TumorVariantCount",
				"TumorSplitVariantCount",
				"TumorReadCount",
				"TumorGenotypeQScore",
				"NormalReferenceCount",
				"NormalSplitReferenceCount",
				"NormalVariantCount",
				"NormalSplitVariantCount",
				"NormalReadCount",
				"NormalGenotypeQScore",
				"Cosmic_Fusion_Counts",
				"repName-repClass-repFamily:-site1",
				"repName-repClass-repFamily:-site2",
				"CC_Chr_Band",
				"CC_Tumour_Types(Somatic)",
				"CC_Cancer_Syndrome",
				"CC_Mutation_Type",
				"CC_Translocation_Partner",
				"DGv_Name-DGv_VarType-site1",
				"DGv_Name-DGv_VarType-site2"])

	outDF.to_csv(outFile, sep='\t', index=False)