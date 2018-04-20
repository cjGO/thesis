





rollingprobes = list()
  for (i in 1:21){
    chrProbes = c()
    print(i) 
    for ( b in c('E','L','K','R')) {
      for (x in 1:11){
        word = paste0(b,x,'_DATA')
        getData = get(word)
        probes = filter(getData[[i]], passCorr5 == TRUE)
        probes = probes$snpP
        chrProbes = c(chrProbes,probes)
      }  
    rollingprobes[[i]] = unique(chrProbes)
    }
}

# rolling probes has the chromosomse * probes to keep

coDom[,colnames(coDom) %in% rollingprobes[[1]]]
