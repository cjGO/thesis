


#############
#
# Author : Colton Gowan
#
# About : This file will contain the python functions created through
#         the duration of my Master's of Science Thesis for bioinformatics at the Roslin Institte
#

###########

### iInput:  a mapping csv where columns are ; PROBE-ID|PROBE-ID|MappingPop#1-Chromosome|MappingPop#1-Location|MappingPop#2-Chromosome|MappingPop#2-Location|...|...|General Consensus-Chromosome|General Consensus - Location

import sys
import re
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy



inFile = sys.argv[1]



def csvLocationsOnly(CSV):
    #Extract Rows with Consensus locations
    consensusProbes = {}
    
    with open(inFile,'r') as i:
        entries = i.readlines()
    for row in entries[1:]:
        rowSplit = row.split(',')
        probeID = rowSplit[0]
        if len(rowSplit[12]) >= 2: ## filter out probes missing mapping information
            rowSplit[13] = rowSplit[13].rstrip('\n')
            allData = [] # list to contain the mapping information from each pop.
            #for i in range(2,12,2):
               # allData.append(rowSplit[i] + '=' + rowSplit[i+1])
            allData = probeID , rowSplit[2]+'='+rowSplit[3],rowSplit[4]+'='+rowSplit[5],rowSplit[6]+'='+rowSplit[7],rowSplit[8]+'='+rowSplit[9],rowSplit[10]+'='+rowSplit[11],rowSplit[12]+'='+rowSplit[13]
            ## converts to List with [ProbeID, MappingPopChr+Loc,.]
            consensusProbes[probeID] = allData
    return consensusProbes

#'AX-94938785': ('AX-94938785=', '=', '=', '2D=282.28', '=', '2D=295.13')}




def splitPopulations(probes):
    #subsets the data provided for each of the 5 populations and consensus map, so each pop has its own dict with dict[probeID]=chromosome:location
    #inpt: 
    popList = {}
    for i in range(0,6):
        popList[i] = []
    for val in probes:
        mapInfo = (probes[val])
        dataAdd = []
        probeID = mapInfo[0]
        for i in range(0,6):
            if len(mapInfo[i]) > 3: # means population has provide data
                popList[i].append(probeID + '#' +mapInfo[i]) # make a place holder a list to fill up, then add it at end so probe iid can be included)
    return(popList)

def getPopData(popdata,N):
    #gets the locations and # of chromosomes probes have been mapped to each population
    wheatChr = []
    for i in range(1,8):
        for g in ['A','B','D']:
            wheatChr.append(str(i)+g)
    for i in range(0,6):
        for x in popdata[N]:
            ID,dat = x.split('+')
            Chr,Loc = dat.split('=')
            chrNum[i].append(Chr)
    return chrNum



def makePopDict():
    popList = {}
    for i in range(0,6):
        popList[i] = []
    return popList

conP = csvLocationsOnly(inFile)
popData = splitPopulations(conP)
mapInfo =  conP
#print(conP)
print(popData[1])

