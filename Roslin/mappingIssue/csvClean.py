


#############
#
# Author : Colton Gowan
#
# About : This file will contain the python functions created through
    #         the duration of my Master's of Science Thesis for bioinformatics at the Roslin Instittepyth
#
###########

### iInput:  a mapping csv where columns are ; PROBE-ID|PROBE-ID|MappingPop#1-Chromosome|MappingPop#1-Location|MappingPop#2-Chromosome|MappingPop#2-Location|...|...|General Consensus-Chromosome|General Consensus - Location

import sys
import re
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy


inFile = sys.argv[1]
domData = sys.argv[2]


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
    for i in range(0,7):
        popList[i] = []
    for val in probes:
        mapInfo = (probes[val])
        dataAdd = []
        probeID = mapInfo[0]
        for i in range(0,7):
            if len(mapInfo[i]) > 3: # means population has provide data
                popList[i].append(probeID + '#' +mapInfo[i]) # make a place holder a list to fill up, then add it at end so probe iid can be included)
    return(popList)

def getPopData(popdata,N):
    #gets the locations and # of chromosomes probes have been mapped to each population
    wheatChr = []
    chrNum = {}
    for i in range(1,8):
        for g in ['A','B','D']:
            wheatChr.append(str(i)+g)
    for i in range(0,6):
        for x in popdata[N]:
            ID,dat = x.split('#')
            chrNum[ID] = []
            Chr,Loc = dat.split('=')
            chrNum[ID].append(Chr+'-'+Loc)
    return chrNum

def makePopDict():
    popList = {}
    for i in range(0,6):
        popList[i] = []
    return popList


def mapPopConDiff(con,pop):
    #calculates the difference in location between the consensus and a given mapping population
    
    #first get all they keys or probes in the given population
    #second collect the locations of those probes in both the consensus and mapping pop
    #third calculate the difference
    #put the differences in a data frame which may be used to visualize in a histogram
    dataDiff = []
    presentProbes = pop.keys()
    for probe in presentProbes:
        popChr,popLoc = pop[probe][0].split('-')
        conChr,conLoc = con[probe][0].split('-')
        diff = float(conLoc) - float(popLoc)
        dataDiff.append(diff)
    return dataDiff

def avgDiff(consensus):
    # this will get the average difference between the consensus and provided map probes
    avgDiff = []
    for i in consensus:
        total = 0
        number = 0
        conChr,consensusLoc = consensus[i][6].split('=')
        for pop in consensus[i][1:6]:
            if len(pop) > 3:
               Chr,Loc = pop.split('=')
               number += 1
               total += float(Loc)
        avg = total / number
        avg = float(consensusLoc) - avg 
        avgDiff.append(avg)
    return avgDiff


def debugPopulation(popSubset,popLocations):
    assert(len(popSubset) == len(popLocations))

def collectIDs(markers):
    #extracts the probe ID's and places them in a lists
    idsOnly = []
    with open(markers,'r') as i:
        x = i.readlines()
    for i in x:
        idsOnly.append(i.rstrip())
    return idsOnly

def findOverlap(newData,conP):
    #subsets te probes + info between a file and the entire mapping info (e.g. we have a file with CoDominant marker ID's and we want to get the mapping info for just those ID's)
    overlapped = {}
    consensusList = conP.keys()
    for valid in newData:
        if valid in consensusList:
            overlapped[valid] = conP[valid]
    return overlapped




conP = csvLocationsOnly(inFile)
popData = splitPopulations(conP)


codom = collectIDs(domData)
x = findOverlap(codom,conP)


#def backToCSV(marker):
    #converts the data back into CSV file

with open(inFile,'r') as i:
    entries = i.readlines()
coDomIDs = x.keys()
for row in entries[1:]:
    rowSplit = row.split(',')
    probeID = rowSplit[0]
    if probeID in coDomIDs:
        print(row.rstrip())













"""
coDomMark = findOverlap(codom,conP)
coDomPop = splitPopulations(coDomMark)

pop1 = getPopData(popData,1)
pop2 = getPopData(popData,2)
pop3 = getPopData(popData,3)
pop4 = getPopData(popData,4)
pop5 = getPopData(popData,5)
popC = getPopData(popData,6)

pop1Locations = mapPopConDiff(popC,pop1)
pop2Locations = mapPopConDiff(popC,pop2)
pop3Locations = mapPopConDiff(popC,pop3)
pop4Locations = mapPopConDiff(popC,pop4)
pop5Locations = mapPopConDiff(popC,pop5)
avgLocations = avgDiff(conP)

fig, axes = plt.subplots(nrows=1, ncols=2)
ax0, ax1 = axes.flatten()

allData=[pop1Locations,pop2Locations,pop3Locations,pop4Locations,pop5Locations,avgLocations]
colors = ['green','red','tan','blue','magenta','black']
labels = ['AxC','SxR','OxS','AxP','CSxP','Consensus']
num_bins=20

ax0.hist(allData,20,histtype='step',color=colors,label=labels)
ax0.legend(prop={'size':10})
ax0.set_title('All probes with Consensus')

coDomMark = findOverlap(codom,conP)
coDomPop = splitPopulations(coDomMark)
popData=coDomPop
pop1 = getPopData(popData,1)
pop2 = getPopData(popData,2)
pop3 = getPopData(popData,3)
pop4 = getPopData(popData,4)
pop5 = getPopData(popData,5)
popC = getPopData(popData,6)

#pop1=AxC
#pop2=SxR
#pop3=OxS
#pop4=AxP
#pop5=CSxP

pop1Locations = mapPopConDiff(popC,pop1)
pop2Locations = mapPopConDiff(popC,pop2)
pop3Locations = mapPopConDiff(popC,pop3)
pop4Locations = mapPopConDiff(popC,pop4)
pop5Locations = mapPopConDiff(popC,pop5)
avgLocations = avgDiff(coDomMark)


allData=[pop1Locations,pop2Locations,pop3Locations,pop4Locations,pop5Locations,avgLocations]
colors = ['green','red','tan','blue','magenta','black']
labels = ['AxC','SxR','OxS','AxP','CSxP','Consensus']
num_bins=20

ax1.hist(allData,20,histtype='step',color=colors,label=labels)
ax1.legend(prop={'size':10})
ax1.set_title('CoDominant marker Overlap')

fig.tight_layout()
plt.show()
"""
