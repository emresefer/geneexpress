import sys
import os

def readBPOnto(fpath):
    """
    """
    ontos = []
    with open(fpath,"r") as infile:
        for line in infile:
            line = line.rstrip()
            ontos.append(line)
    return ontos        
            
ontofpath = "annos.txt"
ontos = readBPOnto(ontofpath)

resfolder = "sortedsubsetontology10"
#resfolder = "sortedsubsetontology5"
newresfolder = "bp_{0}".format(resfolder)
if not os.path.exists(newresfolder):
   os.makedirs(newresfolder) 
for fname in os.listdir(resfolder):
    if fname.find(".tsv")==-1:
       continue 
    fpath = "{0}/{1}".format(resfolder,fname)
    newfpath = "{0}/{1}".format(newresfolder,fname)
    with open(fpath,"r") as infile:
        with open(newfpath,"w") as outfile:
            for line in infile:
                line = line.rstrip()
                if line.startswith("#") or line.startswith("OVERREP") or line.startswith("N"):
                   outfile.write("{0}\n".format(line))
                else:
                   splitted = line.split("\t")
                   if splitted[5] in ontos:
                      outfile.write("{0}\n".format(line))
print "done"                       
exit(1)

fname = "gene_association.mgi"
outfname = "bp_gene_association.mgi"
with open(fname,"r") as infile:
    with open(outfname,"w") as outfile:
       for line in infile:
           line = line.rstrip()
           if line.startswith("!"):
              outfile.write("{0}\n".format(line))
              continue
           splitted = line.split("\t")
           assert splitted[8] in ["P","F","C"]
           if splitted[8] == "P":
              outfile.write("{0}\n".format(line))  
                   
