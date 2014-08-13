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
      modules=go_table[which(as.character(go$weight)==go[i,"geneOntology"]),1]
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
  
  sampledRegions=rbind(specificSample[,c("chrom","chromStart","chromEnd","region")],nonSpecificSample[,c("chrom","chromStart","chromEnd","region")])
  rownames(sampledRegions)=c(1:nrow(sampledRegions))
  colnames(sampledRegions)=c("chrom","chromStart","chromEnd","region")
  return (sampledRegions)
}



#findPossibleGuides(int guideLength, dataframe sampledRegions)
#returns all sequences with length guideLength from sampledRegions that start with CC or end with GG
findPossibleGuides <- function(guideLength, sampledRegions)
{
  require(BSgenome.Hsapiens.UCSC.hg19)
  genome <- BSgenome.Hsapiens.UCSC.hg19
  
  guides=c("chr","start","end","region")
  for(i in 1:nrow(sampledRegions))
  {
    reg=subseq(genome[[as.character(sampledRegions[i,1])]],sampledRegions[i,2],sampledRegions[i,3])
    match <- matchPattern("GG", reg, max.mismatch=0)
    for(k in 1:length(match))
    {
      if(start(match)[k]>(guideLength-1))
      {
        start = end(match)[k] + sampledRegions[i,2]-(guideLength-1)
        end = end(match)[k] + sampledRegions[i,2]
        guides=rbind(guides,c(as.character(sampledRegions[i,1]),start,end,as.character(sampledRegions[i,4])))
      }
    }
  }
  guides=guides[2:nrow(guides),]
  return(guides)
}