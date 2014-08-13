focusRegions = read.table("/srv/gsfs0/projects/kundaje/commonRepository/epigenomeRoadmap/integrative/regulatoryRegions/ENCODERoadmap/reg2map/HoneyBadger2_release/DNase/p10/enh/BED_files_per_sample/regions_enh_E123.bed.gz")
moduleDirectory= "/srv/gsfs0/projects/kundaje/commonRepository/epigenomeRoadmap/integrative/regulatoryRegions/ENCODERoadmap/reg2map/HoneyBadger2_release/DNase/p10/enh/BED_files_per_cluster/"
go_table=read.table("/srv/gsfs0/projects/kundaje/users/summerStudents/2014/nsunil/project2/go_table.txt")
test_minFractionFocusInModule=0.9
max_FractionModuleSpecific=0.7
numRegions=9000

GOOfInterest=c("GO:0040007","GO:0051649","GO:0008037","GO:0030163","GO:0000082",
  "GO:0071702","GO:0008380","GO:0019439","GO:0019221","GO:0010608",
  "GO:0019058","GO:0022402")
#clusters: 12,4,10,15,17,23,31,35,38,42,47,54,62
weightsForGO=rep(50,3,1,6,12,4,37,1,1,2,14,12)
GOInput=as.data.frame(cbind(GOOfInterest,weightsForGO))

focusRegions=setFocusRegions(focusRegions)
sampledRegions=sampleRegions(calcModuleWeightingForSampling(GOInput,go_table))

#16,12
#15,54,35