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
  
  ctSpecific=selectModule()
  #for each module value, look for weighting value on focusRegions and append
  ctSpecific$Weighting=0
  for(i in 1:nrow(ctSpecific))
  {
    index=match(ctSpecific[i,4],focusRegionsWithWeighting[,4])
    if(!is.na(index))
    {
      ctSpecific[i,ncol(ctSpecific)]=as.numeric(as.character((focusRegionsWithWeighting[index,ncol(focusRegionsWithWeighting)])))
    }
  }
  
  if (nrow(ctSpecific) < numCTSpecific)
  {
    numCTspecific=nrow(ctSpecific)
    numNonSpecific = numRegions-numCTSpecific
    print("There are not enough cell type specific regions")
    print(paste(numCTspecific,"cell type specific regions will be sampled"))
    print(paste(numNonSpecific,"nonspecific regions will be sampled"))
  }
  specificSample=ctSpecific[sample(1:nrow(ctSpecific),numCTSpecific,prob=ctSpecific[,ncol(ctSpecific)]),]
  rest=focusRegionsWithWeighting
  for(i in 1:nrow(ctSpecific))
  {
    index=match(paste(ctSpecific[i,1],":",ctSpecific[i,2],"-",ctSpecific[i,3],sep=''),rest[,4])
    if(!is.na(index))
    {
      rest=rest[-index,]
    }
  }
  nonSpecificSample=rest[sample(1:nrow(rest),size=numNonSpecific,prob=rest[,ncol(rest)]),]
  return (rbind(specificSample[,1:4],nonSpecificSample[,1:4]))
}