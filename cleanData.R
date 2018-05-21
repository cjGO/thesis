## Data Script for processing the raw data in the mapping information and probe information datafiles

## Script subsets each cross data by chromosome with correlations



library(dplyr)
library(readr)
library(ggplot2)

options(stringsAsFactors=FALSE)

setwd('C:/Users/Colt/Desktop/shared/may21/')
#setwd('/media/shared/may21/')
mapping = read_csv('../mappedCoDom.csv')
coDom = read_csv('../CoDomEDIT.csv') #raw coDom.csv from GplusE folder with 1 fix for incorrect column names

'%!in%' <- function(x,y)!('%in%'(x,y))



chrList = makeChrList()
makeChrList = function(){
  chrIds = c()
  for (i in c('A','B','D')){
    for (x in 1:7){
      output = paste(x,i,'output',sep='')
      chrIds = c(chrIds,paste(x,i,sep=''))
      
    }
  }
  return(chrIds)
}

remainCrosses = c('E1_','E2_','E3_','E4_','E5_','E6_','E7_','E8_','E9_','E10_','E11_',
                  'K1_','K2_','K3_','K4_','K5_','K6_','K7_','K8_','K9_','K10_','K11_',
                  'L1_','L2_','L3_','L4_','L5_','L6_','L7_','L8_','L9_','L10_','L11_',
                  'R1_','R2_','R3_','R4_','R5_','R6_','R7_','R8_','R9_','R10_','R11_')

remainCrosses = c('E1_','E2_','E3_','E4_','E5_','E6_','E7_','E8_','E9_','E10_','E11_',
                  'K1_','K8_','K10_','K11_',
                  'L1_','L2_','L3_','L4_','L6_','L7_','L8_','L9_','L10_',
                  'R1_','R2_','R3_','R6_','R7_','R8_','R10_','R11_')
crossData2 = c()
for (i in remainCrosses){
  crossData2 = c(crossData2,paste0(i,'DATA'))
}
probeData = c()
for ( i in remainCrosses){
  probeData = c(probeData,paste0(i,'PROBES'))
}
#Data Frame with Information for MAPPED PROBES only
consensus = mapping[c(1,13,14)]
colnames(consensus) = c('ProbeID','Chromosome','Location')
#Removes unmapped probes in CoDom marker set
mappedCoDom = coDom[colnames(coDom) %in% consensus$ProbeID]
mappedCoDom = cbind(coDom[,1],mappedCoDom)
#makes vector with each element corresponding to each chromosome
chrList = makeChrList()

orderCon = list()
for (chr in seq(21)){
  hit = filter(consensus,consensus$Chromosome==chrList[chr]) %>% arrange(Location)
  orderCon[[chr]] = as.data.frame(hit)
}


crossOnly = filter(mappedCoDom, grepl('E1_',CrossID))
# all the columns in chromosome 1 in order
c1probes = orderCon[[1]]$ProbeID
hits = crossOnly[colnames(crossOnly) %in% c1probes]
hits = hits[c1probes]


#gets each cross-chromosome in a list in the correct order
for (cross in remainCrosses){
  print(cross)
  thisCross = list()
  secondList = list()
  crossOnly = filter(mappedCoDom, grepl(cross,CrossID))
  for (chr in seq(21)){
    c1probes = orderCon[[chr]]$ProbeID
    hits = crossOnly[colnames(crossOnly) %in% c1probes]
    hits = hits[c1probes]
    #print(hits[1:5,1:5])
    thisCross[[chr]] = hits
    secondList[[chr]] = orderCon[[chr]]
  }
  assign(paste0(cross,'DATA'),secondList)
  assign(paste0(cross,'PROBES'),thisCross)
}

addCorrData = function(data,probes){
  
  for (c in seq(21)){
    c1 = probes[[c]]
    numProbes = ncol(probes[[c]])
    corrCol = length(numProbes)
    for (i in 2:numProbes){
      corrCol[i] = abs(cor(c1[,i],c1[,(i-1)],use = 'pairwise.complete.obs')) 
    }
    data[[c]]$correlation = corrCol
  }
  return(data)
}
#E2_DATA = addCorrData(E2_DATA,E2_PROBES)

massCor = function(crossData,probeData){
  for (x in 1:length(crossData)){
    print(x)
    assign(crossData[x],addCorrData(get(crossData[x]),get(probeData[x])),envir = globalenv())
  }
}
massCor(crossData2,probeData)





pedigreeSimple = function(p1,p2,off){
  # makes a basic pedigree.txt file for 2 parents only
  numParents = nrow(p1) + nrow(p2)
  numOffspring = nrow(off)
  total = numParents + numOffspring
  print(total)
  ped = matrix(nrow=total,ncol=3)
  print(dim(ped))
  ped[,1] = 1:nrow(ped) # fills the first column which simply contains the numeric ID's of each individual
  ped[1:numParents,2:3] = c(0,0)
  ped[(numParents+1):nrow(ped),2]=1
  ped[(numParents+1):nrow(ped),3]=2
  print(dim(ped))
  write.table(ped,row.names = F,col.names = F,file='Pedigree.txt')
}






