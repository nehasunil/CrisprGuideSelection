focusRegions = read.table("/srv/gsfs0/projects/kundaje/commonRepository/epigenomeRoadmap/integrative/regulatoryRegions/ENCODERoadmap/reg2map/HoneyBadger2_release/DNase/p10/enh/BED_files_per_sample/regions_enh_E123.bed.gz")
moduleDirectory= "/srv/gsfs0/projects/kundaje/commonRepository/epigenomeRoadmap/integrative/regulatoryRegions/ENCODERoadmap/reg2map/HoneyBadger2_release/DNase/p10/enh/BED_files_per_cluster/"
go_table=read.table("/srv/gsfs0/projects/kundaje/users/summerStudents/2014/nsunil/project2/go_table.txt")
test_minFractionfocusRegionsInModule=0.9
max_FractionModuleSpecific=0.7
numRegions=9000

GOOfInterest=c(GO:0040007, and all those under it?)
weightsForGO=rep(1,times=length(GOOfInterest))
GOInput=cbind(GOOfInterest,weightsForGO)

focusRegions=setFocusRegions(focusRegions)
selectedRegions=selectRegions(calcModuleWeightingForSampling(GOInput))