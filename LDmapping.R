#for probe order based on maps

mapping = readxl::read_xlsx('Supplementary_file_3.xlsx')
#mapping[13] = 'consensus'

df = mapping[!is.na(mapping$Position),] # remove columns where the consensus location is unknown


ndf = cbind(df$`35K SNPId`,df$`Consensus_ Chromosome`,df$Position)
ndf = as.data.frame(ndf)
colnames(ndf) = c('SNPid','Chromosome','Position')


numbers = 1:7
letters = c('A','B','D')

chrIds = c()
for (i in c('A','B','D')){
  for (x in 1:7){
#    print(paste(i,x,sep=''))
    output = paste(x,i,'output',sep='')
    chrIds = c(chrIds,paste(x,i,sep=''))
#    assign(output,chrIds)
  }
}

#filter(ndf,Chromosome == '7A')

#orderedProbes = function(ndf,chrIds){
  ##orders probes based on consensus map
for (i in 1:21){
  hit = filter(ndf, Chromosome == chrIds[i])
  outputx = paste(chrIds[i],'order',sep='')
  assign(outputx,as.data.frame(hit))
  }
#}

orderedProbes(ndf,chrIds)


