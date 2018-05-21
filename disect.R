xList = list()
xList[[1]] = c('E1_','Dickens_E1','Diego','') # Multiple Dickens, using Dickens_E1 for all .. for now
xList[[2]] = c('E2_','Croft_E10','Kielder_E10','')
xList[[3]] = c('E3_','Evolution_E8','Kielder_E8','')  
xList[[4]] = c('E4_','Dickens_E1','Santiago_E1','')  
xList[[5]] = c('E5_','Dickens_E2','Revelation_E2','')  
xList[[6]] = c('E6_','Dickens_E4','Kielder_E10','')
xList[[7]] = c('E7_','Dickens_E1','Leeds_E5','')  
xList[[8]] = c('E8_','Evolution_E8','Revelation_E11','')  
xList[[9]] = c('E9_','Croft_E10','Revelation_E9','')  
xList[[10]] = c('E10_','Kielder_E10','Revelation_E11','')  
xList[[11]] = c('E11_','Leeds_E6','Revelation_E11','')


'%!in%' <- function(x,y)!('%in%'(x,y))

kmeansBins = function(chromosome,k,crosses){
  ## Divides up the bins for a given CHROMOSOME NUMBER with location information based on K - clusters with the kmeans algorithm
  all = c()
  for (cross in crosses){
    print(cross)
    #print(cross)
    x = get(cross)
    x = x[[chromosome]][c(1,2,3,4)]
    x = filter(x,correlation>0.5)
    all = rbind(all,x)
  }
  all = all[1:3]
  all = unique(all)
  all = arrange(all,Location)
  cluster=kmeans(all$Location,k)$cluster
  all = cbind(all,cluster)
  #  plot(x=(all$Location),y=(1:length(all$Location)),cex=.8,pch = all$cluster,col=all$cluster,ylab='Index',xlab='centiMorgan', main=paste0('Chr',unique(all$Chromosome),' [clusters=',max(all$cluster),']'))
  
  return(all)
}

probePicker = function(kbins,Nprobes){
  # selects NPROBES from each cluste in kmeansBins outputs...
  # HOW? randomly for now
  # adds these probes to a vector to be used for masking in imputation
  # NOTE : if a cluster only has 1 probe in it and you pick 10 probes, you only get that 1 probe from that cluster
  numClusters = max(kbins$cluster)
  maskProbes = c()
  for (c in seq(numClusters)){
    clust = filter(kbins,cluster==c)
    if (Nprobes > dim(clust)[1]){
      maskProbes = c(maskProbes,clust$ProbeID)
    }else{
      maskProbes = c(maskProbes,sample(clust$ProbeID,Nprobes,replace=FALSE))  
    }
    
  }
  return(maskProbes)
}

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







c = 3
num = 3 # from xList



crossName = xList[[num]][1]
parent1 = xList[[num]][2]
parent2 = xList[[num]][3]
crossChoose = paste0(crossName,'DATA') # string with the name of the DATA file for chosen cross

p1row = filter(mappedCoDom,grepl(parent1,CrossID))
p2row = filter(mappedCoDom,grepl(parent2,CrossID))
offspring = filter(mappedCoDom,grepl(crossName,CrossID))

dataFile = get(paste0(crossName,'DATA'))

hitProbes = orderCon[[c]]$ProbeID # gets the probes on chosen chromosome in the correct order
p1rowX = p1row[colnames(p1row) %in% hitProbes] # filters the probes in parent 1 for those on CHOSEN CHROMOSOME
p1rowX = p1rowX[hitProbes] # puts the correct probes in order
p2rowX = p2row[colnames(p2row) %in% hitProbes]
p2rowX = p2rowX[hitProbes]
offspringX = offspring[colnames(offspring) %in% hitProbes]
offspringX = offspringX[hitProbes]

#processes the data for a given cross into the imputation format (unmasked form)
imputeData = rbind(p1rowX,p2rowX,offspringX)
imputeData = (imputeData + 1) # 0/1/2 format for parents(first 2 rows) and offspring(rest of rows)
indIDs = 1:nrow(imputeData)
imputeData = cbind(indIDs,imputeData)
imputeData[is.na(imputeData)] = 9
###################################################


# we'll test 8-10 markers per chromsome 


write.table(imputeData2,col.names=F,row.names=F,file='MaskedGenotypes.txt')
write.table(tail(imputeData,n=(nrow(imputeData)-2)),col.names=F,row.names=F,file='TrueGenotypes0.txt')

# random masked probes
allProbes = colnames(imputeData)
allBins = unique(length())
goodProbes =(filter(dataFile[[c]],correlation>0.5))
goodBins = length(unique(goodProbes$Location))

randomProbes = sample(allProbes,5)
naiveBins = kmeansBins(c,5,crossData2)
naiveProbes = probePicker(naiveBins,2)
customBins = kmeansBins(c,4,crossChoose)
customProbes = probePicker(customBins,2)

df = E3_DATA[[3]]
dfN = df[df$ProbeID %in% naiveProbes,] # no informative probes
dfC = df[df$ProbeID %in% customProbes,] # 2 informative probes
dfR = df[df$ProbeID %in% randomProbes,] # no informative probes

## this is jacked -> just suppose to be special 'histograms' that
## also show if the probes are informative or not
dfN2 = dfN
dfN2[is.na(dfN2)] = 0
dfN2$correlation = dfN2$correlation > 0.5
gn = ggplot(dfN2,aes(Location)) + geom_bar(aes(fill=correlation)) + xlim(c(0,250)) + scale_fill_discrete(limits = c('TRUE', 'FALSE')) + ggtitle('Naive')

dfC2 = dfC
dfC2[is.na(dfC2)] = 0
dfC2$correlation = dfC2$correlation > 0.5
gc = ggplot(dfC2,aes(Location)) + geom_bar(aes(fill=correlation)) + xlim(c(0,250)) + scale_fill_discrete(limits = c('TRUE', 'FALSE'))+ ggtitle('Custom')

dfR2 = dfR
dfR2[is.na(dfR2)] = 0
dfR2$correlation = dfR2$correlation > 0.5
gr = ggplot(dfR2,aes(Location)) + geom_bar(aes(fill=correlation)) + xlim(c(0,250)) + scale_fill_discrete(limits = c('TRUE', 'FALSE'))+ ggtitle('Random')

gn
gr
gc

filter(E3_DATA[[3]],ProbeID== 'AX-94684907')
filter(imputeData,colnames(imputeData) == 'AX-94684907')

mat2[,'saturn']

imputeData[,'AX-94619024']
imputeData[,'AX-94446053']
imputeData[,'AX-95630061']







