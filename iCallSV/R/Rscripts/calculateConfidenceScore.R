#Created on 03/18/2015
#Description: this will run targetSeqView funtion to predict confidnece score of an SV
#@author:Ronak H Shah
args = commandArgs(TRUE)
calculateConfidenceScore<-function(nodes,bamFilePath,svFile,build,readLength,outdir,outsvFileName) {
	suppressMessages(library(targetSeqView))
	print("Correctly Loaded the Package.")
	registerDoMC(nodes)
	candidateSVs <- new("candidates")
	bamFilePath(candidateSVs) <- bamFilePath
	candidatesFileName(candidateSVs) <- svFile
	build(candidateSVs) <- build
	readLength(candidateSVs) <- readLength
	mmRate(candidateSVs) <- precomputedTargetCapture100bpMMRate()
	indelRate(candidateSVs) <- precomputedTargetCapture100bpIndelRate()
	print("Finished Making the candidates object. Will Now Score all Events.")
	candidateSVs <- fullScoreAndView(candidateSVs, verbose = TRUE, findSplitReads = TRUE)
	print ("Completed Scoring All Events.")
	svScore=candidateSVs@fullScore
	numberOfEvents = length(svScore)
	svData = read.delim(svFile,header=TRUE, sep="\t")
	svData$ProbabilityScore<-0
	print("Making Pdf for each event\n")
	for (i in 1:numberOfEvents ){
		sampleName = svData$SampleDesc[i]
		chr1 = svData$Chr1[i]
		start1 = svData$Start1[i]
		chr2 = svData$Chr2[i]
		start2 = svData$Start2[i]
		outFileName = paste(sampleName,chr1,start1,chr2,start2,"targetSeqViewPlot.pdf",sep="_")
		outFile = paste(outdir,outFileName,sep = "/")
		plotSV(candidateSVs, indices = i, pdfname=outFile)
		svData$ProbabilityScore[i]<-svScore[i]
	}
	outsvFile<-paste(outdir,outsvFileName,sep="/")
	write.table(svData, outsvFile, sep="\t",row.names=FALSE,quote=FALSE)
	print ("Finished Running targetSeqView on the given structural variants File.")
	
}
nodes = as.integer(args[1])
bamFilePath = args[2]
svFile = args[3]
build = args[4]
readLength = as.integer(args[5])
outdir = args[6]
outsvFileName = args[7]
calculateConfidenceScore(nodes,bamFilePath,svFile,build,readLength,outdir,outsvFileName)
