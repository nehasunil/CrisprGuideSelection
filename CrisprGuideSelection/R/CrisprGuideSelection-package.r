#setFocusRegions(BEDfile focusTable) - sets format of focusRegions to have 4 columns, 
#with the fourth column being the name of the region
setFocusRegions <- function(focusTable)
{
  focusRegions=focusTable[,1:3]
  focusRegions[,4]=paste(as.character(focusRegions[,1]),":",focusRegions[,2],"-",focusRegions[,3],sep='')
  colnames(focusRegions)=c("chrom", "chromStart","chromEnd", "region")
  return(focusRegions)
}

#selectModule() - calculates and returns which module is the most cell-type specific. 
#
#The method selectModule() looks through all the Roadmap + Encode modules 
#and chooses the module with the highest percentage of its contents from the focus cell type. 
#If the highest percentage is less than test_minFractionFocusInModule, the module is deemed
#not cell-type specific and NA is returned
selectModule <- function(modDir)
{
  files <- list.files(path=modDir)
  
  clust=read.table(paste(modDir,files[1],sep=''))
  colnames(clust)[1:3]=c("chrom","chromStart","chromEnd")
  clust_regions=paste(as.character(clust$chrom),":",clust$chromStart,"-",clust$chromEnd,sep='')
  intersection=intersect(focusRegions$region,clust_regions)
  count=length(intersection)
  percentage=count/nrow(clust)
  dist=cbind(count,percentage)
  
  for(i in 2:length(files))
  {
    clust=read.table(paste(modDir,files[i],sep=''))
    colnames(clust)[1:3]=c("chrom", "chromStart","chromEnd")
    clust_regions=paste(as.character(clust$chrom),":",clust$chromStart,"-",clust$chromEnd,sep='')
    intersection=intersect(focusRegions$region,clust_regions)
    count=length(intersection)
    percentage=count/nrow(clust)
    dist=rbind(dist,c(count,percentage))
  }
  dist=as.data.frame(dist)
  if (max(dist$percentage)>test_minFractionFocusInModule)
  {
    maxInd=which.max(dist$percentage)
    module=read.table(paste(modDir,files[maxInd],sep=''))
    return(module)
  }
}

#calcModuleWeightingForSampling() - calculates and returns weighting of regions based on module gene ontologies
#parameter = data frame with gene ontology in first column and weighting for each GO in second column
#
#In order to give certain gene ontologies preference over others, 
#users can input weights for each region. The method calcWeightingForSampling() allows users 
#to assign a weight to each type of gene ontology and the method assigns a weight to each region 
#based on its gene ontology. For each gene ontology of interest, the method finds all modules with that gene ontology.
#For each module, the method finds all the intersections with focusRegions and sets the weighting vector 
#to the value assigned by the go parameter.
calcModuleWeightingForSampling <- function(go,go_table)
{
  #order go in ascending order according to weights
  colnames(go)=c("geneOntology","weight")
  go=as.data.frame(go[order(go$weight),])
  if(nrow(go)>0)
  {
    weighting = rep(0,times=nrow(focusRegions))
    for(i in 1:nrow(go))
    {
      modules=go_table[which(as.character(go_table[,2])==go[i,"geneOntology"]),1]
      for(k in 1:length(modules))
      {
        clust=read.table(paste(moduleDirectory,"cluster_",modules[k],".bed.gz",sep=''))
        colnames(clust)[1:3]=c("chrom", "chromStart","chromEnd")
        clust_regions=paste(as.character(clust$chrom),":",clust$chromStart,"-",clust$chromEnd,sep='')
        intersection=which(focusRegions$region %in% clust_regions)
        weighting[intersection]=go[i,"weight"]
      }
    }
    return(weighting)
  }
  else
  {
    return(rep(1,times=nrow(focusRegions)))
  }
}

#sampleRegions(vector weighted) - samples specified number of regions according to input values 
#
#The method sampleRegions() uses the selected module to calculate cell type specific regions. 
#The parameter for weights is a vector that corresponds to the assigned weights for each region 
#in focusRegions. The method first chooses the appropriate number of cell-type specific regions 
#according the user-inputed value fractionCellTypeSpecific from the module returned by selectModule(). 
#The code checks that there are sufficient cell-type specific regions to be chosen from. 
#If there are not enough, the number of cell-type specific and non cell-type specific regions are adjusted.
#Then, the non-specific regions are sampled from the regions in focusRegions that are not in the selected module. 
#The method returns a data frame with the first four columns of the specific and non-specific regions.
sampleRegions <- function(weighted)
{
  numCTSpecific=max_FractionModuleSpecific*numRegions
  numNonSpecific=numRegions-numCTSpecific
  
  focusRegionsWithWeighting=cbind(focusRegions,weighted)
  colnames(focusRegionsWithWeighting)=c("chrom","chromStart","chromEnd","region","weighting")
  
  ctSpecific=selectModule(moduleDirectory)
  ctSpecific=ctSpecific[,1:3]
  colnames(ctSpecific)[1:3]=c("chrom","chromStart","chromEnd")
  ctSpecific$region=paste(as.character(ctSpecific$chrom),":",ctSpecific$chromStart,"-",ctSpecific$chromEnd,sep='')
  
  if (nrow(ctSpecific) < numCTSpecific)
  {
    numCTspecific=nrow(ctSpecific)
    numNonSpecific = numRegions-numCTSpecific
    print("There are not enough cell type specific regions")
    print(paste(numCTspecific,"cell type specific regions will be sampled"))
    print(paste(numNonSpecific,"nonspecific regions will be sampled"))
  }
  
  specificSample=ctSpecific[sample(1:nrow(ctSpecific),numCTSpecific),]
  
  rest=focusRegionsWithWeighting
  intersections= which(focusRegionsWithWeighting$region %in% ctSpecific$region)
  rest=rest[-intersections,]
  
  if(sum(as.numeric(rest$weighting))==0)
  {
    print("There are no regions of the gene ontology of interest outside of the cell-type specific sample")
    return(specificSample[,c("chrom","chromStart","chromEnd","region")])
  }
  
  nonSpecificSample=rest[sample(1:nrow(rest),size=numNonSpecific,prob=rest[,ncol(rest)]),]
  sampledRegions=rbind(specificSample[,c("chrom","chromStart","chromEnd","region")],nonSpecificSample[,c("chrom","chromStart","chromEnd","region")])
  rownames(sampledRegions)=c(1:nrow(sampledRegions))
  colnames(sampledRegions)=c("chrom","chromStart","chromEnd","region")
  return (sampledRegions)
}



#findPossibleGuides(int guideLength, dataframe sampledRegions)
#returns all sequences with length guideLength from sampledRegions that start with CC or end with GG
findPossibleGuides <- function(sampledRegions,outputDir)
{
  require(BSgenome.Hsapiens.UCSC.hg19)
  genome <- BSgenome.Hsapiens.UCSC.hg19
  
  guides=c("chrom","chromStart","chromEnd","region")
  guidesNeg=c("chrom","chromStart","chromEnd","region")
  for (i in 1:nrow(sampledRegions))
  {
    for (guideLength in 20:23)
    {
      reg=subseq(genome[[as.character(sampledRegions[i,1])]],sampledRegions[i,2]-1,sampledRegions[i,3])
      #positive strand
      guideSeq=paste('G',paste(rep('N',times=guideLength-3),sep='',collapse=''),'GG',collapse='',sep='')
      guidesSeq=unique(as.character(matchPattern(guideSeq,reg,fixed=F)))
      
      match=matchPattern(guideSeq,reg,fixed=F)
      
      for(k in 1:length(match))
      {
        start = start(match)[k] + sampledRegions[i,2]
        end = end(match)[k] + sampledRegions[i,2] + 1
        guides=rbind(guides,c(as.character(sampledRegions[i,1]),start,end,as.character(sampledRegions[i,4])))
      }
      
      #negative strand
      guideSeqNeg=paste('CC',paste(rep('N',times=guideLength-3),sep='',collapse=''),'C',collapse='',sep='')
      guidesSeqNeg=unique(as.character(matchPattern(guideSeqNeg,reg,fixed=F)))
      guidesSeqNegRC=as.character(reverseComplement(DNAString(guidesSeqNeg[1])))
      for(j in 1:length(guidesSeqNeg))
      {
        guidesSeqNegRC=c(guidesSeqNegRC,as.character(reverseComplement(DNAString(guidesSeqNeg[j]))))
      }
      
      matchNeg=matchPattern(guideSeqNeg,reg,fixed=F)
      
      for(k in 1:length(matchNeg))
      {
        start = start(matchNeg)[k] + sampledRegions[i,2]
        end = end(matchNeg)[k] + sampledRegions[i,2] + 1
        guidesNeg=rbind(guidesNeg,c(as.character(sampledRegions[i,1]),start,end,as.character(sampledRegions[i,4])))
      }
    }
  }
  
  write.table(guidesSeq,paste(outputDir,"guideSeqPos.txt",sep=''),sep='\t')
  print("Guide sequences on positive strand saved to:")
  print(paste(outputDir,"guideSeqPos.txt",sep=''))
  write.table(guidesSeqNegRC,paste(outputDir,"guideSeqNeg.txt",sep=''),sep='\t')
  print("Guide sequences on negative strand saved to:")
  print(paste(outputDir,"guideSeqNeg.txt",sep=''))
  
  guides=guides[2:nrow(guides),]
  colnames(guides)=c("chrom","chromStart","chromEnd","sampledRegion")
  rownames(guides)=c(1:nrow(guides))
  write.table(guides,paste(outputDir,"guidePos.bed",sep=''),sep='\t')
  print("Guide BED file for positive strand saved to:")
  print(paste(outputDir,"guidePos.bed",sep=''))
  guidesNeg=guidesNeg[2:nrow(guidesNeg),]
  colnames(guidesNeg)=c("chrom","chromStart","chromEnd","sampledRegion")
  rownames(guidesNeg)=c(1:nrow(guidesNeg))
  write.table(guidesNeg,paste(outputDir,"guideNeg.bed",sep=''),sep='\t')
  print("Guide BED file for negative strand saved to:")
  print(paste(outputDir,"guideNeg.bed",sep=''))
}