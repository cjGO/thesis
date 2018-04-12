#for probe order based on maps

library(dplyr)
library(readr)

setwd('C:/Users/Colt/Documents/masters/RCode/')
mapping = read_csv('mappedCoDom.csv')
coDom = read_csv('CoDomEDIT.csv') #raw coDom.csv from GplusE folder with 1 fix for incorrect column names

#mapping[13] = 'consensus'

mead = function(df){
  #head() for dataframes/matrix
  df[1:10,1:5]
}

getMapped = function(df){
  # Input: the mapping information CSV file [mapping]
  # removes rows(probes) that don't have consensus locations
  df = df[!is.na(df$X14),] # remove columns where the consensus location is unknown
  return(df)
}


minConsensusInfo = function(df){
  # Input : getMapped output
  # Output : same information with just the relevant columns of : ProbeID | Consensus Chr | Consenus Location
  df = cbind(df$ProbeID,df$X13,df$X14)
  df = as.data.frame(df)
  colnames(df) = c('SNPid','Chromosome','Position')
  return(df)
}

getCoDom = function(coDom,mapped){
    # grabs the column names with probe id's in the mapped coDom markers
  markers = colnames(coDom[2:length(coDom)])
  return(markers)
}

mappedCoDom = function(coDom,mappedDF){
  #this removes the probes in the codominant marker set that do not have consensus locations
  CrossID = coDom$CrossID
  mappedProbes = mappedDF$ProbeID
  mappedCo = coDom[,colnames(coDom) %in% mappedProbes]
  mappedCo = cbind(CrossID,mappedCo)
  return(mappedCo)
}


# coDom = coDom
# crossid = coDom$CrossID
# mappedProbes = mapping$ProbeID
# 
# mappedCo = coDom[,colnames(coDom) %in% mappedProbes]
# mappedCo = cbind(crossid,mappedCo)
# return(mappedCo)


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

probeLocater = function(cpDF,chrID){
  #cpDF = called probes Data Frame ( has been filtered )
  #chrID = vector corresponding to each chromosome
  #returns a LIST with the data frames for each chromosome 1-7 [A], 8-14[B], 15-21[D]
  allchr = list()
  for (i in 1:21){
    hit = filter(cpDF, Chromosome == chrID[i]) %>% arrange(Position)
    print(class(hit))
    allchr[i] = as.data.frame(hit)
    #hit = as.character(hit$SNPid)
    #outputx = paste(chrID[i],'order',sep='')
  }
  return(allchr)
}





overlapMapDom = function(domProbes, orderMapped){
  for (i in 1:21){
    orderMapped[[i]] = intersect(orderMapped[[i]],domProbes)
  }
  return(orderMapped)
}


LEprobes = function(cross){
  # makes a vector containing probes where the corresponding SNP is fixed or missing in a cross
  # input: cross ( all rows pertaining to a single cross eg. E1_02,E1_03... with columns as probes)
  dropProbes = c()
  count = 0
  for (probe in colnames(cross[2:length(cross)])){
    count = count + 1
    print(count)
    if (cross[probe] %>% distinct() %>% rownames() %>% length() == 1){
      print(probe)
      dropProbes = c(dropProbes,probe)
    }
  }
  return(dropProbes)
}

LEprobes2 = function(cross){
  # makes a vector containing probes where the corresponding SNP is fixed or missing in a cross
  # input: cross ( all rows pertaining to a single cross eg. E1_02,E1_03... with columns as probes)
  dropProbes = c()
  count = 0
  for (i in length(cross)){
    if (length(rownames(distinct(cross[i])))){
      dropProbes = c(dropProbes,colnames(cross[i]))
    }
  }
  return(dropProbes)
}

simpleApply = function(x){
  length(unique(x))
 }


filterNoise = function(co){
  sd = sapply(co,simpleApply)
  tfn = sd>1
  cot = co[,tfn==TRUE]
  return(cot)
}

# x = filter(filteredCoDom, grepl('E1_',CrossID))
# sd = sapply(x,simpleApply)
# tfn = sd>1
# x2=x[,tfn==TRUE]
# 



LDprobesFilter = function(cross,droppedProbes){
  #for each population this will filter the probes in LD so the correlations can be computed
  LDprobes = cross[,colnames(cross[,1:length(cross)])%!in%droppedProbes]
  
  return(LDprobes)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

finallyLD = function(cleanDF,crossSubset,chrNumber){
  desiredProbes = cleanDF[[chrNumber]]
  desiredProbes = desiredProbes %in% colnames(cleanDF)
  p = length(desiredProbes)
  count = c()
  for (probe in (2:p)){
    xb = desiredProbes[probe-1]
    xf = desiredProbes[probe]
    xa = desiredProbes[probe+1]
    print(paste(xf,xb))
    bProbe = grep(xb,colnames(crossSubset))
    focalProbe = grep(xf,colnames(crossSubset))
    aProbe = grep(xa,colnames(crossSubset))
    c1 = cor(crossSubset[focalProbe],crossSubset[bProbe],use='complete.obs')
    c2 = cor(crossSubset[focalProbe],crossSubset[aProbe],use='complete.obs')
    c3 = c1 + c2 / 2
    count = c(count,c3)
  }
  return(count)
}

compatable = function(c1,c2){
  #checks if 2 probe data can have cor() calculated
  check = cbind(c1,c2)
  tst = any(complete.cases(check) == TRUE)
  return(tst)
}




finallyLD2 = function(cleanDF,crossSubset,chrNumber){
  desiredProbes = cleanDF[[chrNumber]]
  desiredProbes = desiredProbes %in% colnames(cleanDF)
  p = length(desiredProbes)
  count = c()
  for (probe in (2:p)){
    xb = desiredProbes[probe-1]
    xf = desiredProbes[probe]
    xa = desiredProbes[probe+1]
    print(paste(xf,xb))
    bProbe = grep(xb,colnames(crossSubset))
    focalProbe = grep(xf,colnames(crossSubset))
    aProbe = grep(xa,colnames(crossSubset))
  }
}
    #c1 = cor(crossSubset[focalProbe],crossSubset[bProbe],use='complete.obs')
    #c2 = cor(crossSubset[focalProbe],crossSubset[aProbe],use='complete.obs')
    #c3 = c1 + c2 / 2
    #count = c(count,c3)

mappedCoDom = function(coDom,mappedDF){
  #this removes the probes in the codominant marker set that do not have consensus locations
  CrossID = coDom$CrossID
  mappedProbes = mappedDF$ProbeID
  mappedCo = coDom[,colnames(coDom) %in% mappedProbes]
  mappedCo = cbind(CrossID,mappedCo)
  return(mappedCo)
}

overlapMapDom = function(crossSubset, orderMapped){
  for (i in 1:21){
    orderMapped[[i]] = intersect(dirty[[1]],colnames(f2[2:length(f2)]))
  }
  return(orderMapped)
}

coDom=coDom
mapping=mapping

getMap = getMapped(mapping) #filters the 35k probes w/o mapping info
minCon = minConsensusInfo(getMap)
getCoDomIDs = getCoDom(coDom)
mappedCoDomx = mappedCoDom(coDom,getMap) ## Rows = crosses / columns = mapped probes


filteredCoDom = mappedCoDom(coDom,mapping) # removes the probes w/o mapping info in the coDom.csv
filteredCoDomx = filterNoise(filteredCoDom) # filterNoise works on SUBSETS not the entire coDom.csv
#for single chromosome test
f = filter(filteredCoDom, grepl('E1_',CrossID))
dim(f)
f2 = filterNoise(f) #removes probes that are completely fixed or missing in each cross's subset
dim(f2)





domProbes = getCoDom(coDom)
#badProbes = LEprobes(f)

coDom = coDom
cList = makeChrList() 

df = mapping %>% getMapped() %>% minConsensusInfo()
dirty = probeLocater(df,cList) # list with ordered probes for each of the 21 chromosomes
clean = overlapMapDom(f2,dirty)
