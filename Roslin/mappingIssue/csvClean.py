


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
            allData = rowSplit[2]+'='+rowSplit[3],rowSplit[4]+'='+rowSplit[5],rowSplit[6]+'='+rowSplit[7],rowSplit[8]+'='+rowSplit[9],rowSplit[10]+'='+rowSplit[11],rowSplit[12]+'='+rowSplit[13]
            ## converts to List with [ProbeID, MappingPopChr+Loc,.]
            consensusProbes[probeID] = allData
    return consensusProbes

def splitPopulations(probes):
    #subsets the data provided for each of the 5 populations and consensus map
    pop0=[]
    pop1=[]
    pop2=[]
    pop3=[]
    pop4=[]
    pop5 = []
    popList = {}
    for i in range(0,6):
        popList[i] = []

    for val in probes:
        mapInfo = (probes[val])
        for i in range(0,6):
            if len(mapInfo[i]) > 2: # means population has provide data
                popList[i].append(mapInfo[i])
    print(popList[1])




conP = csvLocationsOnly(inFile)
splitPopulations(conP)




