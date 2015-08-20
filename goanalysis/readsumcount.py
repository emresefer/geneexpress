import sys
import os
import math

def readMiRNAs():
    """
    """
    fname = "LCM-mir.txt"
    mirs = []
    with open(fname,"r") as infile:
        for line in infile:
            uselines = line.split("\r")
            for uind,useline in enumerate(uselines[1:]):
                if uind in range(0,5) or uind >= 604:
                   continue 
                splitted = uselines[uind+1].split("\t")
                gene = splitted[0]
                mirs.append(gene)
    return mirs

def readClusters():
    """
    """
    clust2mir = {}
    cid = None
    clustfname = "clusters.txt"
    with open(clustfname,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if not line.startswith("m"):
               cid = int(line)
               clust2mir[cid] = []
            else:
               clust2mir[cid].append(line)
    return clust2mir


datamirs = readMiRNAs()
print len(datamirs)
clust2mir = readClusters()

mir2gene = {}
count = 0
fname = "Summary_Counts.txt"
with open(fname,"r") as infile:
    for line in infile:
        line = line.rstrip()
        if count == 0:
           count += 1
           continue
        splitted = line.split("\t")
        gene,mirname = splitted[1],splitted[-3]
        if splitted[3] != "10090":
           continue
        if mirname not in datamirs:
           continue  
        mir2gene.setdefault(mirname,set())
        mir2gene[mirname].add(gene)
print len(mir2gene.keys())
avggene = 0.0
for geneim in mir2gene.keys():
    avggene += len(mir2gene[geneim])
avggene /= float(len(mir2gene.keys()))
print avggene    

genedir = "clustgenes"
if not os.path.exists(genedir):
   os.makedirs(genedir)
for clustid,curmirs in clust2mir.items():
    curgenes = []
    for curmir in curmirs:
        if not mir2gene.has_key(curmir):
           continue 
        curgenes.extend(mir2gene[curmir])
    print clustid,len(curgenes)
    genepath = "{0}/{1}.txt".format(genedir,clustid)
    if os.path.exists(genepath):
       continue 
    with open(genepath,"w") as outfile:
       for curgene in curgenes:
           outfile.write("{0}\n".format(curgene))      
          




     
        
