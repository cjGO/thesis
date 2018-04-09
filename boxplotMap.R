#boxplot for mapping populations with data in the consensus population


mapping = readxl::read_xlsx('Supplementary_file_3.xlsx')
mapping[13] = 'consensus'

consensusMapped = filter(mapping, "Consensus_ Chromosome" >= 0)
df = consensusMapped


AxC = df$cM
SxR = df$cM__1
CSxP = df$cM__2
OxS = df$cM__3
AxP = df$cM__4
conM = df$Position

mapNames = c('AxC','SxR','CSxP','OxS','AxP')

xdiff = function(map,con){
  mapData = c()
  for (i in 1:length(map)){
    if (is.na(map[i]) == FALSE){
      diffX = con[i] - map[i]
      mapData=c(mapData,diffX)
    }
  }
  return(mapData)
}

AxC = xdiff(AxC,conM)
length(AxC)  = 9171

SxR = xdiff(SxR,conM)
length(SxR)  = 9171
#sum(is.na(S))

CSxP = xdiff(CSxP,conM)
length(CSxP)  = 9171
#sum(is.na(CS))

OxS = xdiff(OxS,conM)
length(OxS)  = 9171
#sum(is.na(O))

AxP = xdiff(AxP,conM)
length(AxP) = 9171
#sum(is.na(AP))

Ax
locations = data.frame(SxR,OxS,CSxP,AxP,AxC)

require(reshape2)
ggplot(data = melt(locations), aes(x=variable, y=value, alpha=.01)) + geom_boxplot(aes(fill=variable)) + coord_flip()
#meltLoc = melt(locations)
