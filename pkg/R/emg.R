## Read EBI TAD

require(sads)

##############################################################################################

## Functions for working with EMP data

read.project.csv<-function(fileName,projectID,...) {
    summ=read.csv(fileName,stringsAsFactors=FALSE,...)
    attr(summ,"project.id")=projectID
    summ
}

getProjectSummary<-function(projectID) {
    url=paste("https://www.ebi.ac.uk/metagenomics/projects",projectID,"overview/doExport",sep="/")
    summ=read.project.csv(url,projectID)
    summ
}

projectSamples<-function(summ) {
    unique(sort(summ$Sample.ID))
}

projectRuns<-function(summ) {
    summ$Run.ID
}

runsBySample<-function(summ,sampleID) {
    summ$Run.ID[summ$Sample.ID==sampleID]
}

otu.url<-function(summ,runID) {
    runData=summ[summ$Run.ID==runID,]
    projectID=attr(summ,"project.id")
    url=paste("https://www.ebi.ac.uk/metagenomics//projects",projectID,"samples",runData["Sample.ID"],"runs",runID,"results/versions",sprintf("%.1f",runData["Release.version"]),"taxonomy/OTU-TSV",sep="/")
    url
}

read.otu.tsv<-function(fileName,...) {
    otu = read.delim(fileName,header=FALSE,skip=2,colClasses=c("character","numeric","character"),stringsAsFactors=FALSE,...)
    names(otu) = c("OTU","Count","Tax")
    otu[order(-otu$Count),]
}

getRunOtu<-function(summ,runID,verb=FALSE,plot.preston=FALSE) {
    url=otu.url(summ,runID)
    if (verb)
        message(runID)
    otu=read.otu.tsv(url)
    if (plot.preston)
        plot(octav(otu$Count),main=paste("Preston plot for",runID))
    otu
}

mergeOtu<-function(otu1,otu2) {
  stack=rbind(otu1[,1:2],otu2[,1:2])
  comb=tapply(stack$Count,stack$OTU,sum)
  otu=data.frame(OTU=as.vector(rownames(comb)),Count=as.vector(comb))
  otu[order(-otu$Count),]
}

getSampleOtu<-function(summ,sampleID,verb=TRUE,plot.preston=FALSE) {
    runs=runsBySample(summ,sampleID)
    runData=list()
    for (run in runs) {
        if (verb)
            message(paste(run,", ",sep=""),appendLF=FALSE)
        otu=getRunOtu(summ,run)
        if (plot.preston)
            plot(octav(otu$Count),main=paste("Preston plot for",run))
        runData[[run]]=otu
    }
    if (verb)
        message("END.")
    Reduce(mergeOtu,runData)
}


##############################################################################################



## eof

