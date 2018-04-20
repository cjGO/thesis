#for probe order based on maps

library(dplyr)
library(readr)
library(ggplot2)

options(stringsAsFactors=FALSE)

setwd('C:/Users/Colt/Documents/masters/RCode/')
mapping = read_csv('mappedCoDom.csv')
coDom = read_csv('CoDomEDIT.csv') #raw coDom.csv from GplusE folder with 1 fix for incorrect column names

keptX = c('dickens','diego','E')
keptX = c('croft','kielder','E')
keptX = c('evolution','kielder','E')
keptX = c('dickens','santiago','E')
keptX = c('dickens','revelation','E')
keptX = c('dickens','kielder','E')
keptX = c('dickens','leeds','E')
keptX = c('evolution','revelation','E')
keptX = c('croft','revelation','E')
keptX = c('kielder','revelation','E')
keptX = c('leeds','revelation','E')





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
    hit = filter(cpDF, Chromosome == chrID[i])
    hit2 = arrange(hit, as.numeric(as.character(hit$Position)))
    #print(class(hit))
    allchr[[i]] = as.data.frame(hit2)
    #hit = as.character(hit$SNPid)
    #outputx = paste(chrID[i],'order',sep='')
  }
  return(allchr)
}

probeLocater2 = function(cpDF,chrID){
  #cpDF = called probes Data Frame ( has been filtered )
  #chrID = vector corresponding to each chromosome
  #returns a LIST with the data frames for each chromosome 1-7 [A], 8-14[B], 15-21[D]
  allchr = list()
  for (i in 1:21){
    hit = filter(cpDF, Chromosome == chrID[i])
    hit2 = arrange(hit, as.numeric(as.character(hit$Position)))
    #print(class(hit))
    allchr[[i]] = as.data.frame(hit2)
    #hit = as.character(hit$SNPid)
    #outputx = paste(chrID[i],'order',sep='')
  }
  return(allchr)
}


overlapMapDom = function(probeMin, chrSubs){
  ## Pulls the relevant rows for each cross in each
  # domprobes = all the probeIDs in the coDom dataset
  # chrSubs = each Chromosome probe subset 
  orderMapped=list()
  for (i in 1:21){
    #print(i)
    #print(chrSubs[[i]])
    orderMapped[[i]] = chrSubs[[i]][chrSubs[[i]]$SNPid %in% probeMin,]
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

haldane = function(p1,p2){
  d = p2 - p1 # haldane
  calc = (1/2)*(1-exp(-2*(d)/100))
  #calc2 = (1-exp(-2*(d)))/2
  return(calc)
}

finallyLD2 = function(cleanDF,crossSubset,chrNumber){
  desiredProbes = cleanDF[[chrNumber]]$SNPid # probes ID's on given chromosome 1 - 21
  #print(length(desiredProbes))
  
  
  desiredLengths = cleanDF[[chrNumber]]$Position
  f3 = crossSubset[,colnames(crossSubset) %in% desiredProbes] # pulls columns w/ probes on chosen chromosome
  f4 = f3[as.vector(desiredProbes)] #puts the columns in the correct order based on mappinig data
  
  snpP = c()
  corP = c()
  lenP = c()
  halP = c()
  testR = c()
  p = length(desiredProbes)
  #print(p)
  for (i in 2:(p-1)){
    #print(i)
    xb = as.vector(unlist((f4[i-1])))#, digits=4 #probe behind
    xf = as.vector(unlist((f4[i]))) #probe focal
    xa = as.vector(unlist((f4[i+1]))) #probe ahead
    hit1 = cor(xb,xf,use="pairwise.complete")
    hit2 = cor(xa,xf,use="pairwise.complete")
    hit3 = (abs(hit1)+abs(hit2))/2
    
    hb = as.vector(cleanDF[[chrNumber]]$Position[i-1])
    hb = as.numeric(hb,digits=4)
    hf = as.vector(cleanDF[[chrNumber]]$Position[i])
    hf = as.numeric(hf,digits=4)
    #ha = as.integer(cleanDF[[chrNumber]]$Position[i+1])
    hal1 = haldane(hb,hf)
    #print(hb)
    #hal2 = haldane(hf,ha)
    #hal3 = max(hal1+hal2) - 1
    
    #testAdd = paste(hb,hf,sep='+')
    #testR = c(testR,as.integer(testAdd))
    #testR = c(testR,as.character(testAdd))
    
    #hal3 = (hal1+hal2)/2
    
    #print(paste(hit1,hit2,hit3))
    #print(as.numeric(cleanDF[[chrNumber]]$Position[i]))
  
    snpP = c(snpP,as.vector(cleanDF[[chrNumber]]$SNPid[i]))
    lenP = c(lenP,(as.vector(cleanDF[[chrNumber]]$Position[i])))
    corP = c(corP,hit3)
    halP = c(halP, (1-hal1))
  }
  outDF = data.frame(snpP,lenP,corP,halP)##testR)
  return(outDF)
  #return(outP[!is.na(outP)])
}


finallyLD3 = function(cleanDF,crossSubset,chrNumber){
  desiredProbes = cleanDF[[chrNumber]]$SNPid # probes ID's on given chromosome 1 - 21
  #print(length(desiredProbes))
  
  
  desiredLengths = cleanDF[[chrNumber]]$Position
  f3 = crossSubset[,colnames(crossSubset) %in% desiredProbes] # pulls columns w/ probes on chosen chromosome
  f4 = f3[as.vector(desiredProbes)] #puts the columns in the correct order based on mappinig data
  
  snpP = c()
  corP = c()
  lenP = c()
  halP = c()
  testR = c()
  p = length(desiredProbes)
  #print(p)
  for (i in 2:(p-1)){
    #print(i)
    xb = as.vector(unlist((f4[i-1])))#, digits=4 #probe behind
    xf = as.vector(unlist((f4[i]))) #probe focal
    xa = as.vector(unlist((f4[i+1]))) #probe ahead
    hit1 = cor(xb,xf,use="pairwise.complete")
    hit2 = cor(xa,xf,use="pairwise.complete")
    hit3 = (abs(hit1)+abs(hit2))/2
    
    hb = as.vector(cleanDF[[chrNumber]]$Position[i-1])
    hb = as.numeric(hb,digits=4)
    hf = as.vector(cleanDF[[chrNumber]]$Position[i])
    hf = as.numeric(hf,digits=4)
    #ha = as.integer(cleanDF[[chrNumber]]$Position[i+1])
    hal1 = haldane(hb,hf)
    #print(hb)
    #hal2 = haldane(hf,ha)
    #hal3 = max(hal1+hal2) - 1
    
    #testAdd = paste(hb,hf,sep='+')
    #testR = c(testR,as.integer(testAdd))
    #testR = c(testR,as.character(testAdd))
    
    #hal3 = (hal1+hal2)/2
    
    #print(paste(hit1,hit2,hit3))
    #print(as.numeric(cleanDF[[chrNumber]]$Position[i]))
    
    snpP = c(snpP,as.vector(cleanDF[[chrNumber]]$SNPid[i]))
    lenP = c(lenP,(as.vector(cleanDF[[chrNumber]]$Position[i])))
    corP = c(corP,hit3)
    halP = c(halP, (1-hal1))
  }
  outDF = data.frame(snpP,lenP,corP,halP)##testR)
  return(outDF)
  #return(outP[!is.na(outP)])
}
singleLD = function(cleanDF,crossSubset,chrNumber){
  desiredProbes = cleanDF[[chrNumber]]$SNPid # probes ID's on given chromosome 1 - 21
  #print(length(desiredProbes))
  
  
  desiredLengths = cleanDF[[chrNumber]]$Position
  f3 = crossSubset[,colnames(crossSubset) %in% desiredProbes] # pulls columns w/ probes on chosen chromosome
  f4 = f3[as.vector(desiredProbes)] #puts the columns in the correct order based on mappinig data
  
  diffL = c()
  snpP = c()
  corP = c()
  lenP = c()
  halP = c()
  testR = c()
  p = length(desiredProbes)
  #print(p)
  for (i in 2:(p-1)){
    #print(i)
    xb = as.vector(unlist((f4[i-1])))#, digits=4 #probe behind
    xf = as.vector(unlist((f4[i]))) #probe focal
    hit1 = abs(cor(xb,xf,use="pairwise.complete"))
    
    hb = as.vector(cleanDF[[chrNumber]]$Position[i-1])
    hb = as.numeric(hb,digits=4)
    hf = as.vector(cleanDF[[chrNumber]]$Position[i])
    hf = as.numeric(hf,digits=4)
    #ha = as.integer(cleanDF[[chrNumber]]$Position[i+1])
    hal1 = haldane(hb,hf)
    
    lengthDiff = hf - hb
    
    diffL = c(diffL,lengthDiff)
    snpP = c(snpP,as.vector(cleanDF[[chrNumber]]$SNPid[i]))
    lenP = c(lenP,(as.vector(cleanDF[[chrNumber]]$Position[i])))
    corP = c(corP,hit1)
    halP = c(halP, (1-hal1))
  }
  outDF = data.frame(snpP,lenP,corP,halP,diffL)##testR)
  return(outDF)
  #return(outP[!is.na(outP)])
}


#x = finallyLD3(dirty,f,1)
#y = finallyLD3(clean,f,1)
# x
# plot(x$halP)

mappedCoDom = function(coDom,mappedDF){
  #this removes the probes in the codominant marker set that do not have consensus locations
  CrossID = coDom$CrossID
  mappedProbes = mappedDF$ProbeID
  mappedCo = coDom[,colnames(coDom) %in% mappedProbes]
  mappedCo = cbind(CrossID,mappedCo)
  return(mappedCo)
}

breedtag = function(){
  xID = c()
  BreederIDs = c('E','K','R','L')
  CrossNumbers = 1:11
  for (bid in BreederIDs){
    for (cn in CrossNumbers){
      tag=(paste0(bid,cn,'_'))
      xID=c(xID,tag)
    }
  }
  return(xID)
}

overlapMapDom2 = function(probeMin, chrSubs){
  ## Pulls the relevant rows for each cross in each(
  # domprobes = all the probeIDs in the coDom dataset
  # chrSubs = each Chromosome probe subset [dirty]
  
  #need to make it so that if the probes are not in the dirty list that they are only
  #converted to 
  
  orderMapped=list()
  for (i in 1:21){
    #print(i)
    #print(chrSubs[[i]])
    orderMapped[[i]] = chrSubs[[i]]
    
    
    orderMapped[[i]] = chrSubs[[i]][chrSubs[[i]]$SNPid %in% probeMin,]
    
    
  }
  return(orderMapped)
}


cleandirty = function(x,y){
  ### gets relevant information in the correct index for comparison
  x = mutate(x,new_corr=NA)
  rownames(x) = x$snpP
  rownames(y) = y$snpP
  
  for (i in rownames(y)){
    realCorr = y[i,]$corP
    x[i,]$new_corr = realCorr
  }
  
  return(x)
}

coDom=coDom
mapping=mapping

getMap = getMapped(mapping) #filters the 35k probes w/o mapping info
minCon = minConsensusInfo(getMap) 
cList = makeChrList()
filteredCoDom = mappedCoDom(coDom,mapping) # removes the probes w/o mapping info in the coDom.csv


## ^^^^^^^ Same for all
## vvvvvvv cross dependent
f = filter(filteredCoDom, grepl('E10_',CrossID))
dim(f)
fprobes = colnames(f[1,])


f2 = filterNoise(f) #removes probes that are completely fixed or missing in each cross's subset
dim(f2)
f2probes = colnames(f2[1,])
dirty = probeLocater(minCon,cList)
clean = overlapMapDom(f2probes,dirty) # gives us the probe | chromosome | position ... of a given cross for every chromosome
x = finallyLD3(dirty,f2,1)



collectData = function(btag){
  
  getMap = getMapped(mapping) #filters the 35k probes w/o mapping info
  minCon = minConsensusInfo(getMap) 
  cList = makeChrList()
  filteredCoDom = mappedCoDom(coDom,mapping) # removes the probes w/o mapping info in the coDom.csv
  
  ## goes through each cross and collects the relevant information
  allDATA = list()
  count = 0
  chro = 1:21
  for( b in btag){
    thisCross = list()
    count = 0
    nam = paste0(b,'DATA')
    f = filter(filteredCoDom,grepl(b,CrossID))
    f2 = filterNoise(f)
    f2probes = colnames(f2[1,])
    dirty = probeLocater(minCon,cList)
    clean = overlapMapDom(f2probes,dirty)
    for (c in 1:21){
      count = count + 1
      x = singleLD(dirty,f,c)
      y = singleLD(clean,f,c)
      z = cleandirty(x,y)
      
      
      goodProbes5 = filter(z,new_corr > 0.5)$snpP
      goodProbes9 = filter(z,new_corr > 0.9)$snpP
      goodProbes0 = filter(z,new_corr > 0)$snpP
      
      passCorr5 = z$snpP %in% goodProbes5
      passCorr9 = z$snpP %in% goodProbes9
      passCorr0 = z$snpP %in% goodProbes0
      
      #print(keepP)
      #print(z$snpP)
      #print(class(z$snpP))
      #print(class(keepP))
      #@corrxT = z$snpP %in% goodProbes$keepP
      z = select(z,snpP,lenP,halP,new_corr,diffL)
      z = cbind(z,passCorr0,passCorr5,passCorr9)
      
      #print(x)
      thisCross[[count]] = z
      #assign(nyam,x,envir=globalenv())
    
      }
    print(nam)
    assign(nam,thisCross,envir=globalenv())
  }
  #print(allDATA[1])
  #return(allDATA)
}
btag= breedtag()
collectData('E1_')

"""
###
f2probes = colnames(f2[1,])
dirty = probeLocater(minCon,cList)
clean = overlapMapDom(f2probes,dirty)
####
plot(x)
hist(x)



par(mfrow=c(3,2))

## make graph showing haldane and correlations

ggplot()



# x2 = finallyLD3(dirty,f,1)
# y2 = finallyLD3(clean,f,1)
# z = cleandirty(x2,y2)

singleCrossNChr = function(Letter,Number){
  par(mfrow=c(3,4))
  for (i in 1:11){
    v = paste0(Letter,i)
    corr_adj = get(paste(v,'DATA',sep='_'))
    plot(x=corr_adj[[Number]]$diffL,corr_adj[[Number]]$new_corr)
  }
}
singleCrossNChr('K',5)

twoCrossOneC = function(L1,L2,L3,L4,Num){
  par(mfrow=c(3,4))
  for (i in 1:3){
    v1 = paste0(L1,i)
    v2 = paste0(L2,i)
    v3 =  paste0(L3,i)
    v4 = paste0(L4,i)
    adj1 = get(paste(v1,'DATA',sep='_'))
    adj2 = get(paste(v2,'DATA',sep='_'))
    adj3 = get(paste(v3,'DATA',sep='_'))
    adj4 = get(paste(v4,'DATA',sep='_'))
    
    plot(x=adj1[[Num]]$diffL,adj1[[Num]]$new_corr,main =paste('E',i,'-',cid[Num]))
    plot(x=adj2[[Num]]$diffL,adj2[[Num]]$new_corr,main=paste("L",i,'-',cid[Num]))
    plot(x=adj3[[Num]]$diffL,adj3[[Num]]$new_corr,main=paste('R',i,'-',cid[Num]))
    plot(x=adj4[[Num]]$diffL,adj4[[Num]]$new_corr,main=paste('K',i,'-',cid[Num]))
    #plot(x=adj1[[Num]]$lenP,adj1[[Num]]$new_corr,main ='E')
    #plot(x=adj2[[Num]]$lenP,adj2[[Num]]$new_corr,main="L")
    #plot(x=adj3[[Num]]$lenP,adj3[[Num]]$new_corr,main='R')
    #plot(x=adj4[[Num]]$lenP,adj4[[Num]]$new_corr,main='K')
  }
}

# [[Num]]$diffL
twoCrossOneC('E','L','K','R',17)
""
for (i in 1:21){
  twoCrossOneC('E','L','K','R',i)
  
}
"""