# filtering script


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

cid = makeChrList()

PROCESScorr = function(Cross, COR){
    ## pass a chr number and a cross ID to get the probes which are above a correlation threshold
  
  rollingProbes = c() ## keeps track of the probess passing all tests
  data = paste0(Cross,'_DATA')
  dataBack = get(data) # entire list for each cross
  print(data)
  print(class(dataBack))
  
  for (i in 1:5){
    #data = paste0(Cross,'_DATA')
    #print(data)
    #print(chrNum)
    #print(class(dataBack))
    #dataBack = get(data) # entire list for each cross
    #print(dataBack[[i]])
    goodprobes = filter(dataBack[[i]], new_corr>COR)
    #print(goodprobes)
    keepP = goodprobes$snpP # as vector
    corrT = dataBack[[i]]$snpP %in% keepP
    dataBack[[i]] = cbind(dataBack[[i]],corrT)
    print(dataBack[[i]])
    
    assign(data[[i]],dataBack[[i]],envir = globalenv())
    
    
    #assign(data,thisCross,envir=globalenv())
    
    
    #add a T/F column for rows(probes) w/ corr above 0.5
  }
}

PROCESScorr('E1',0.5)

Ech1 = filterProbes(1,'E',.2)
Lch1 = filterProbes(1,'L')
Kch1 = filterProbes(1,'K')
Rch1 = filterProbes(1,'R')
totCh1 = c(Ech1,Lch1,Kch1,Rch1)
length(unique(totCh1)) # 367 probes have at least 1 cross above 0.5
length(filter(minCon,Chromosome=='1A')$Chromosome)




findGood = function(minCon,CORR){

  # Makes a List which contains the probe ID's of the high corrrelation
    
  bid = c('E','L','K','R')
  cid = makeChrList()
  corrz = c(0.5,0.7,0.9)
  
  outTable = data.frame()
  goodLoc = list()
  holder = c()
  for (i in 1:21){
    holder = c()
    for (c in CORR){
      Ech = filterProbes(i,'E',c)
      Lch = filterProbes(i,'L',c)
      Kch = filterProbes(i,'K',c)
      Rch = filterProbes(i,'R',c)
      totCh = c(Ech,Lch,Kch,Rch)
      hitID = unique(totCh)
      #hits = length(unique(totCh))
      #all = length(filter(minCon,Chromosome== cid[i])$Chromosome)
      #percent = ( hits / all ) * 100
      #print(paste('cor: ',c,' ','chr',cid[i], ':',hits,'from out of', all, 'which is percent ->',percent,sep=' '))
      holder = c(holder,hitID)
      
    }
    goodLoc[[i]] = holder
    }
  return(goodLoc)
  }
  
x=  findGood(minCon,.5)


cid = makeChrList()


goodMap = function(goodPro,dirty){
  cid = cid # from makeChrList()
  for (i in 1:21){
    iDirty = dirty[[i]]
    iGood = goodPro[[i]]
    chromosomeLength = max(iDirty$Position)
    iPulled = iDirty[iDirty$SNPid %in% iGood,]
    
    mapped = iDirty$SNPid %in% iGood
    iClean=cbind(iDirty,mapped)
    index= rownames(iClean)
    xClean = cbind(iClean,index)
    
    
    plot(xClean$Position,xClean$index,
         col = as.factor(xClean$mapped),
         pch=20,
         cex=.5,
         xlab = 'Consensus Position',
         ylab = 'Probe Index',
         main = paste('Correlation Pass ',cid[i]))
    
    legend("bottomright",
           legend=levels(as.factor(xClean$mapped)),
           col=c('black','red'),
           pch=20)
    
  }
}



goodPro = findGood(minCon,.5)
goodMap(goodPro,dirty)




















































































####
#junk code
#vvvvvvvvv

idirty = dirty[[1]]
iGood = goodPro[[1]]
chromosomeLength = max(iDirty$Position)
#iPulled = iDirty[iDirty$SNPid %in% iGood,]
mapped = idirty$SNPid %in% iGood 
iClean=cbind(idirty,mapped)
index = rownames(iClean)
xClean = cbind(iClean,index)
qplot(xClean$Position,xClean$index)


fixx = xClean %>%
  mutate(Position = factor(Position)) %>%
  arrange(Position)
qplot(Position,index,col=mapping)


goodMap(goodPro,dirty)

pheno <- pheno %>%
  mutate(Genotype =  factor(Genotipe, levels = Gorder)) %>%
  arrange(Genotype) 

positionsFix = iClean %>%
  mutate(Position = factor(Position,levels= index)) %>%
  arrange(Position)


xClean = cbind(iClean,index)
qplot(index,Position,data=xClean)



Ech1 = stringentProbes(1,'E')
Lch1 = filterProbes(1,'L')
Kch1 = filterProbes(1,'K')
Rch1 = filterProbes(1,'R')
totCh1 = c(Ech1,Lch1,Kch1,Rch1)
length(unique(totCh1)) # 367 probes have at least 1 cross above 0.5
length(filter(minCon,Chromosome=='1A')$Chromosome)
