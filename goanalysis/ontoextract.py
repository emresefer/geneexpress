import os
import sys
import networkx as nx
import math

fname = "go.obo"
start = False
goid = None
funcs = []
with open(fname,"r") as infile:
    for line in infile:
        line = line.rstrip()
        if line.startswith("ontology:"):
           start = True
           continue
        if start:
           if line.startswith("id:"):
              goid = line.split("id: ")[1].replace(" ","")
           elif line.startswith("namespace:") and line.find("biological_process")!=-1:
              funcs.append(goid)

annopath = "annos.txt"
with open(annopath,"w") as outfile:
    for anno in funcs:
        outfile.write("{0}\n".format(anno))
