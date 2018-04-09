###
#
# R library for Master's of Theis
#
###

# Mapping population Analysis

LEprobes = function(cross){
  # makes a vector containing probes which the corresponding SNP is fixed or missing
  dropProbes = c()
  for (probe in colnames(crosses[2:length(crosses)])){
    if (crosses %>% select(probe) %>% distinct() %>% rownames() %>% length() == 1){
      dropProbes = c(dropProbes,probe)
    }
  }
  return(dropProbes)
}

LDprobesFilter = function(cross,droppedProbes){
  #for each population this will filter the probes in LD so the correlations can be computed
  LDprobes = cross[,colnames(cross[,1:length(cross)])%!in%droppedProbes]
  
  return(LDprobes)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

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
    hit = arrange(filter(df, Chromosome == chrID[i]), Position)
    print(class(hit))
    allchr[i] = as.data.frame(hit)
  }
  return(allchr)
}