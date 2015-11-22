'''
Created on November 19, 2015
Description: This module will make directory structure for running analysis
@author: Ronak H Shah
'''
'''
::Inputs::
args: Arguments passsed to iCallSV
config: configuration file passed to iCallSV

'''

def launch_delly_for_different_analysis_type(args,config)
	#Run Delly for Deletion

        del_vcf = rd.run(
        delly=config.get("SVcaller","DELLY" ),
        analysisType="DEL",
        reference=config.get("ReferenceFasta","REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.patientId,
        mapq=config.get("ParametersToRunDelly","MAPQ" ),
        excludeRegions=config.get("ExcludeRegion","EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=True)
        #Run Delly for duplication
        dup_vcf = rd.run(
        delly=config.get("SVcaller","DELLY" ),
        analysisType="DUP",
        reference=config.get("ReferenceFasta","REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.patientId,
        mapq=config.get("ParametersToRunDelly","MAPQ" ),
        excludeRegions=config.get("ExcludeRegion","EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=True)
        #Run Delly for inversion
        inv_vcf = rd.run(
        delly=config.get("SVcaller","DELLY" ),
        analysisType="INV",
        reference=config.get("ReferenceFasta","REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.patientId,
        mapq=config.get("ParametersToRunDelly","MAPQ" ),
        excludeRegions=config.get("ExcludeRegion","EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=True)
        #Run Delly for Translocation
        bnd_vcf = rd.run(
        delly=config.get("SVcaller","DELLY" ),
        analysisType="TRA",
        reference=config.get("ReferenceFasta","REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.patientId,
        mapq=config.get("ParametersToRunDelly","MAPQ" ),
        excludeRegions=config.get("ExcludeRegion","EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=True)
        #Run Delly for insertion
        ins_vcf = rd.run(
        delly=config.get("SVcaller","DELLY" ),
        analysisType="INS",
        reference=config.get("ReferenceFasta","REFFASTA"),
        controlBam=args.controlBam,
        caseBam=args.caseBam,
        caseId=args.patientId,
        mapq=config.get("ParametersToRunDelly","MAPQ" ),
        excludeRegions=config.get("ExcludeRegion","EXREGIONS"),
        outputdir=sampleOutdirForDelly,
        verbose=verbose,
        debug=True)