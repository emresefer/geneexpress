import sys
import os
import math
import scipy
from copy import deepcopy
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
import matplotlib
import random
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import cPickle as cpickle
import gzip
import scipy.cluster



def makeplotMainAll(xvals,yvaldictout,points,knots,gene,outspl,rempoints,remyvals,plotpath,allknot,allspline,usetimes,trixvals,triyvals):    
    """makes main plot
    Args:
       xvals,yvaldictout,points:
       knots,gene,outspl:
       rempoints,remyvals:
       plotpath:
       allknot,allspline:
    Returns:
    """
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 30 
    MARKERSIZE = 25
    DPI = 300

    plotpath = plotpath.replace(".","__")
    times = usetimes
    yvals = [item for algo in yvaldictout.keys() for item in yvaldictout[algo]] + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] + remyvals + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for ttime in times] 
    ymax,ymin = max(yvals), min(yvals)
    #plt.ylim(ymin-1.0,ymax+0.1)
    plt.ylim(ymin-0.05,ymax+0.01)
    plt.xlim(0.0,max(xvals)+0.5)
    locpos = 4
    xlabel = "Days"
    ylabel = ""
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    symmap = {}
    symmap["marker"] = {"spline":"p", "SplineFit":"*", "Knots":"p"}
    symmap["colors"] = {"spline":"r", "SplineFit":"g", "Knots":"k"}
    symmap["labels"] = {"spline":"Fitted Spline", "SplineFit":"Selected points {0}".format(len(points)), "Knots":"Knots {0}".format(len(knots))}

    assert sorted(knots) == knots
    knotvals = [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] 
    for algo in yvaldictout.keys():
        if algo in ["Knots"]:
           continue 
        plt.plot(points,yvaldictout[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    #plt.plot(knots,knotvals,marker=symmap["marker"]["Knots"],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"]["Knots"],label=symmap["labels"]["Knots"])
    plt.plot(usetimes,outspl(usetimes),symmap["colors"]["spline"],lw=3,label=symmap["labels"]["spline"])
    plt.plot(rempoints,remyvals,marker='o',markersize=MARKERSIZE,linestyle='None',color='b',label='Rem. points')

    allyvals = [scipy.interpolate.splev(knot,allspline._eval_args, der=0, ext=0) for tknot in allknot]
    #plt.plot(usetimes,allspline(usetimes),'y',lw=6,label='All')
    
    knotys = [scipy.interpolate.splev(tknot,allspline._eval_args, der=0, ext=0) for tknot in allknot]
    plt.plot(allknot,knotys,marker='+',markersize=MARKERSIZE,linestyle='None',color='b',label='Knots all')
       
    ax = plt.axes()        
    ax.xaxis.grid()
    #plt.plot(trixvals[0:2],triyvals[0:2],lw=3,color='y',label='Linear fit')
    #plt.plot(trixvals[1:],triyvals[1:],lw=3,color='y',label='Linear fit')    
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13) 
    plt.savefig(plotpath, dpi=DPI)
    
    

def makeplotMain(xvals,yvaldictout,points,knots,gene,outspl,rempoints,remyvals,plotpath,usetimes,trixvals,triyvals,legendbound=None):    
    """makes main plot
    Args:
       xvals,yvaldictout,points:
       knots,gene,outspl:
       rempoints,remyvals:
       plotpath:
    Returns:
    """
    plt.clf()
    plt.rc('font', family='serif', size=35)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 20 
    MARKERSIZE = 25
    DPI = 300

    plotpath = plotpath.replace(".","__")
    times = usetimes
    yvals = [item for algo in yvaldictout.keys() for item in yvaldictout[algo]] + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] + remyvals + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for ttime in times]
    if legendbound != None:
       ymax,ymin = legendbound[1]+0.1, legendbound[0]-0.1
       plt.ylim(ymin,ymax)
    else:       
       ymax,ymin = max(yvals), min(yvals)
       plt.ylim(ymin-0.05,ymax+0.01)
    plt.xlim(0.0,max(xvals)+0.5)
    locpos = 4
    xlabel = "Days"
    ylabel = ""
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    symmap = {}
    symmap["marker"] = {"spline":"p", "SplineFit":"*", "Knots":"p"}
    symmap["colors"] = {"spline":"r", "SplineFit":"g", "Knots":"k"}
    symmap["labels"] = {"spline":"Fitted Spline", "SplineFit":"Selected points {0}".format(len(points)), "Knots":"Knots {0}".format(len(knots))}

    assert sorted(knots) == knots
    knotvals = [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] 
    for algo in yvaldictout.keys():
        plt.plot(points,yvaldictout[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    plt.plot(knots,knotvals,marker=symmap["marker"]["Knots"],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"]["Knots"],label=symmap["labels"]["Knots"])
    plt.plot(usetimes,outspl(usetimes),symmap["colors"]["spline"],lw=3,label=symmap["labels"]["spline"])
    plt.plot(rempoints,remyvals,marker='o',markersize=MARKERSIZE,linestyle='None',color='b',label='Rem. points')
    #plt.plot([0.5,7.0],triyvals[0:2],lw=3,color='y',label='Linear fit')
    #plt.plot([7.0,28.0],triyvals[1:],lw=3,color='y',label='Linear fit')
     
    ax = plt.axes()        
    ax.xaxis.grid()  
    if legendbound == None: 
       plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13) 
    plt.savefig(plotpath, dpi=DPI)
    

def featureData(fname):
    """read feature data
    Args:
    Returns:
    """
    gene2pos,pos2gene = {}, {}
    count = 0
    with open(fname,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if count == 0:
               count += 1
               continue 
            splitted = line.split(",")
            genestr, posstr = [item.replace("\"","") for item in splitted[0],splitted[3]]
            gene = genestr.split(".")[1]
            chro,curpos = posstr.split(".")
            gene2pos.setdefault(gene,[])
            gene2pos[gene].append((chro,curpos))
            pos2gene[(chro,curpos)] = gene
    return gene2pos,pos2gene

def readSampleTarget(fname):
    """
    """
    name2time,time2name = {},{}
    count = 0
    with open(fname,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if count == 0:
               count += 1
               continue 
            namestr,timestr = line.split(",")[0:2]
            name2time[namestr] = timestr
            time2name[timestr] = namestr
    return name2time,time2name  

def readData(fname):
    """reads data
    """
    data = []
    ind2name = []
    count,nacount = 0,0
    ind2gene = []
    with open(fname,"r") as infile:
        for line in infile:
            line = line.rstrip()
            splitted = line.split(",")
            if count == 0:
               count += 1
               ind2gene = [item.replace("\"","") for item in splitted[2:]]
               continue
            ind2name.append(splitted[1].replace("\"",""))
            vals = [item.replace("\"","") for item in splitted[2:]]
            nacount += sum([1 for tval in vals if tval == "N/A"])
            data.append([float(item) if item!="N/A" else 0.0 for item in vals])
    print "data has {0} N/As".format(nacount)
    return data,ind2name,ind2gene

def chroAnalyze(ind2gene):
    """meth data in terms of chromosome analyze
    Args:
       ind2gene:
    Returns:   
    """
    chro2spot = {}
    for gene in ind2gene:
        chro = int(gene.split(".")[0].replace("chr",""))
        chro2spot.setdefault(chro,[])
        chro2spot[chro].append(gene)
    fname = "MGIBatchReport_20150828_200311.txt"
    gene2pos = {}
    count = 0
    with open(fname,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if count == 0:
               count += 1
               continue
            splitted = line.split("\t")
            gene,chro,start,end = splitted[0], splitted[-4], int(splitted[-2]), int(splitted[-1])
            gene2pos.setdefault(gene,[])
            gene2pos[gene].append((chro,start,end))
    gene2found = {}        
    for chro in chro2spot.keys():
        for tgene in chro2spot[chro]:
            tpos = int(tgene.split(".")[1])
            found = []
            for tgene2,vallist in gene2pos.items():
                for tchro,tstart,tend in vallist:
                    if str(chro) == str(tchro) and tpos >= tstart and tpos <= tend:
                       found.append(tgene2)
            gene2found[tgene] = found           
    return gene2found
    

def makeData(data,ind2gene,ind2name,name2time):
    """makes readable data out of all files
    Args:
       data:
       ind2gene,ind2name:
       name2time: 
    Returns:
       moddata:
    """
    moddata = []
    genecount,timecount = len(data[0]), len(data)
    for gind in xrange(genecount):
        curgene = ind2gene[gind]
        curdata = {}
        for tind in xrange(timecount):
            rowname = ind2name[tind]
            timestr = name2time[rowname]
            if not timestr.startswith("P"):
               continue 
            ttime = float(timestr.replace("P","").replace("a","").replace("b","").replace("c",""))
            curdata.setdefault(ttime,[])
            curdata[ttime].append(round(data[tind][gind]/100.0,5))
        medval = np.median(curdata[0.5])
        curdata = {ttime: [round(item-medval,5) for item in curdata[ttime]] for ttime in curdata.keys()}
        moddata.append(deepcopy(curdata))   
    return moddata

def transformData(moddata):
    """modify data in terms of log transform
    Args:
       moddata:
    Returns:
       logmoddata: 
    """
    logmoddata = []
    for rowdata in moddata:
        preval = math.log(np.median(rowdata[0.5])+0.0001,2.0)
        curlogval = []
        for ttime in sorted(rowdata.keys()):
            meanval = math.log(np.median(rowdata[ttime])+0.0001,2.0)
            curlogval.append(meanval-preval)
        logmoddata.append(list(curlogval))
    return logmoddata

def saveData(moddata,outfname):
    """save data to file
    Args:
       moddata:
       outfname:
    Returns:
    """
    outtimes = sorted(moddata[0].keys())
    with open(outfname,"w") as outfile:
        outfile.write("\t".join(["X"]+[str(item) for item in outtimes])+"\n")
        for rowdata in moddata:
            outvals = [np.median(rowdata[outtime]) for outtime in outtimes]
            outfile.write("\t".join(["X"]+[str(item) for item in outvals])+"\n")


def estQual(tusespl,x2val,points,rempoints):
    """estimate solution quality
    Args:
       usespl,x2val:
       points,rempoints:
    """
    sumval = 0.0
    for rpoint in rempoints:
        cury = scipy.interpolate.splev(rpoint,tusespl._eval_args, der=0, ext=0)
        curyval = sum([(cury-item)**2 for item in x2val[rpoint]])
        sumval += curyval/float(len(x2val[rpoint]))    
    #sumval = sum([((scipy.interpolate.splev(rpoint,usespl._eval_args, der=0, ext=0)-item)**2)/float(len(x2val[rpoint])) for rpoint in rempoints for item in x2val[rpoint]])
    #sumval = sum([(scipy.interpolate.splev(rpoint,usespl._eval_args, der=0, ext=0)-x2val[rpoint])**2 for rpoint in rempoints])
    return sumval                  

            
def runGreedy(times,usedata,yvallist,count,weights,initmode="change"):
    """runs greedy method
    Args:
       times,usedata:
       count,weights:
    Returns:   
    """
    def convertForm(curtimes,usedata):
        """converts form
        Args:
           curtimes:
        Returns:
           yvallist:   
        """        		
        tyvallist = []
        for cind,cdata in enumerate(usedata):
            calldata = []
            for tind,ttime in enumerate(curtimes):
                calldata.extend(usedata[cind][ttime])
            tyvallist.append(list(calldata))
            #assert len(calldata) == len(usetimes)
        return tyvallist

    if False:
     globtimes = set()
     for cind,cdata in enumerate(usedata):
        globtimes |= set(usedata[cind].keys())
     for cind,cdata in enumerate(usedata):
        assert len(globtimes) == len(usedata[cind].keys())
     usetimes = [ttime for tind,ttime in enumerate(times) for tind2 in xrange(len(usedata[0][ttime]))]
     yvallist = convertForm(times) #??
        
    fixedpoints = [0.5,26.0]
     
    reglambda = 0.0001
    if True:
       avgyvallist = []
       for cind,cdata in enumerate(usedata):
           cavgdata = [np.median(list(usedata[cind][ttime])) for tind,ttime in enumerate(times)]
           avgyvallist.append(list(cavgdata)) 
       avgchange = {}
       for tind in xrange(1,len(times)-1):
           ysumval = 0.0
           for tyval in avgyvallist:
               ysumval += abs(tyval[tind] - tyval[tind-1]) + abs(tyval[tind+1] - tyval[tind])
           avgchange[times[tind]] = ysumval
       for fixedpoint in fixedpoints:
           if avgchange.has_key(fixedpoint):
              del avgchange[fixedpoint]  
       points = [ttime[0] for ttime in sorted(avgchange.items(), key=lambda x: x[1], reverse=True)][0:count-2]
       
       points = fixedpoints + points
       points = sorted(points)
       #multipoints = [point for point in points for tind2 in xrange(len(usedata[0][point]))]
       rempoints = list(set(times) - set(points))

       mapval = {tval:tind for tind,tval in enumerate(times)}
       x2val = [{tval:avgyvallist[yind][tind] for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
       x2valblock = [{tval:list(usedata[yind][tval]) for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
       yvals = [[avgyvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
       #yvals = convertForm(points)
        
       minval = 10000000.0
       tind = 0
       y2knots = []
       outsplines = []    
       sumval = 0.0
       for yind,curyvals in enumerate(yvals):
           spl = scipy.interpolate.UnivariateSpline(points, curyvals, s=reglambda, k=3)
           if reglambda == 0.0:
              for tpon in points:
                  assert abs(scipy.interpolate.splev(tpon,spl._eval_args, der=0, ext=0)-x2val[yind][tpon]) < 0.00001
           outsplines.append(deepcopy(spl))
           y2knots.append(list(spl.get_knots()))
           tsumval = weights[yind]*estQual(spl,x2valblock[yind],points,rempoints)
           sumval += tsumval

       while True:
         print sumval,sumval/float(len(yvallist)*(len(times)-count))
         minsol = None
         for addpoint in rempoints:
            for delpoint in points:
                if delpoint in fixedpoints:
                   continue 
                newrempoints, newpoints = list(rempoints), list(points)
                newrempoints.remove(addpoint)
                newrempoints.append(delpoint)
                newpoints.remove(delpoint)
                newpoints.append(addpoint)
                newpoints = sorted(newpoints) 
                
                tempoutsplines = []
                globsumval = 0.0
                for yind,curyvals in enumerate(yvals):
                    inyvals = [x2val[yind][rpoint] for rpoint in newpoints]
                    inspl = scipy.interpolate.UnivariateSpline(newpoints, inyvals, s=reglambda, k=3)
                    
                    #spl2 = scipy.interpolate.InterpolatedUnivariateSpline(sortednewpoints, sortedinyvals,k=3)
                    #for tpon in newpoints:
                    #    assert abs(spl2.__call__(tpon)-x2val[yind][tpon]) < 0.00001
                    #    #assert abs(scipy.interpolate.splev(tpon,spl2._eval_args, der=0, ext=0)-x2val[yind][tpon]) < 0.00001
                    if reglambda == 0.0:
                       for tpon in newpoints:
                           assert abs(scipy.interpolate.splev(tpon,inspl._eval_args, der=0, ext=0)-x2val[yind][tpon])<0.00001
                    tempoutsplines.append(deepcopy(inspl))
                    tsumval = weights[yind]*estQual(inspl,x2valblock[yind],newpoints,newrempoints)
                    globsumval += tsumval
                if globsumval < sumval:
                   sumval = globsumval
                   minsol = (addpoint,delpoint)
                   outsplines = list(tempoutsplines)
                    
         if minsol == None:
           break
        
         rempoints.remove(minsol[0])
         rempoints.append(minsol[1])
         points.remove(minsol[1])
         points.append(minsol[0])
         points = sorted(points)
         yvals = [[avgyvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(avgyvallist))]
       assert points == sorted(points)
       print points
       #exit(1)

       rempoints = list(set(times) - set(points))
       yvals = [[avgyvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
       tsumval = 0.0
       y2knots = []
       outsplines = []
       for yind,curyvals in enumerate(yvals):
         inyvals = [x2val[yind][rpoint] for rpoint in points]
         tinspl = scipy.interpolate.UnivariateSpline(points, inyvals, s=reglambda, k=3)
         outsplines.append(deepcopy(tinspl))
         y2knots.append(list(tinspl.get_knots()))
         locsumval = weights[yind]*estQual(tinspl,x2valblock[yind],points,rempoints)
         tsumval += locsumval
       print tsumval
       print sumval/float(len(yvallist)*(len(times)-count))
       print tsumval/float(len(yvallist)*(len(times)-count))
       #exit(1)

       #assert abs(tsumval-sumval) < 0.0001     
       avgsum=0.0
       for yind in xrange(len(y2knots)):
         avgsum += len(y2knots[yind])
       print avgsum
       print avgsum/float(len(y2knots))
       print "infom"
       #exit(1)
       return sumval,sumval/float(len(yvallist)*(len(times)-count)),points,yvals,y2knots,outsplines


    rempoints = list(set(times) - set(fixedpoints))
    mapval = {tval:tind for tind,tval in enumerate(times)}
    x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
    x2valblock = [{tval:list(usedata[yind][tval]) for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
    yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in fixedpoints] for yind in xrange(len(yvallist))]
    reglambda = 75.0 #1.0 #50.0 #1.0
    #1.0 -> 805
    y2knots = []
    outsplines = []
    #[0.5, 2.5, 6.0, 11.0, 19.0, 28.0].
    sumval = 100000000000.0
    minsol = None

    #userempoints = [1.5, 1.0, 2.0, 3.0, 4.0, 5.0, 4.5, 9.0, 10.0, 11.0, 13.0, 15.0, 17.0, 18.0, 21.0, 25.0, 26.0, 8.5, 6.5, 6.0, 3.5, 2.5, 7.5, 5.5]
    #print len(userempoints)
    #exit(1)

    #area 0.391676970384
    #userempoints = [1.0, 8.5, 12.0, 23.0]
    userempoints = [1.0, 6.0, 13.5, 23.0]
    #userempoints = rempoints
    allvals = [selpoints for selpoints in itertools.combinations(userempoints, count-2)]

    reglambda = 50.0
    fixedpoints = [0.5, 1.5, 2.5, 3.0, 4.5, 6.0, 7.5, 8.5, 9.5, 10.0, 12.5, 18.0, 24.0, 28.0]
    allvals = [[1.0]]
    sol2val = {}
    random.shuffle(allvals)
    count = 0
    for selpoints in allvals:
        count += 1
        if count == 5000:
           break 
        #tflag = True
        #for mypoint in selpoints:
        #    if mypoint not in [1.5,15.0,17.0,19.0]:
            #if mypoint not in [17.0, 19.0, 11.0, 11.5, 6.0, 1.0,3.0]:
            #if mypoint not in [15.0, 19.0, 9.5, 9.0, 8.0, 3.0, 13.0, 7.5, 7.0, 10.5, 11.0, 6.0, 1.0,4.0,6.5]:
        #    if mypoint not in [15.0, 19.0, 9.5, 9.0, 8.0, 7.5, 7.0, 5.5, 4.5, 3.5, 3.0,2.0,1.5,10.5, 11.0, 6.0, 1.0,4.0,6.5,2.5,5.0,10.0]:
            #if mypoint not in [1.0, 3.0, 4.0, 6.0, 15.0, 6.5,11.0, 19.0,13.0,17.0]:
        #       tflag = False
        #       break
        #if not tflag:
        #   continue 
        if random.random() <= 1.0:
           print count
           print sumval,minsol,sumval/float(len(yvallist)*(len(times)-len(selpoints)-len(fixedpoints)))
           avglens = [len(list(tspl.get_knots())) for tspl in outsplines]
           print np.mean(avglens)
           print "multiple sols:"
           import operator
           sorted_data = sorted(sol2val.items(), key=operator.itemgetter(1))
           for keystr,optval in sorted_data[0:10]:
               print sorted([float(item) for item in keystr.split("_")]),optval,optval/float(len(yvals)*len(times))
        
        usepoints = sorted(list(selpoints) + fixedpoints)
        userempoints = list(set(rempoints) - set(selpoints))
        tempoutsplines = []
        globsumval = 0.0
        curknots = []

        #multipoints = [upoint for upoint in usepoints for tind2 in xrange(len(usedata[0][upoint]))]
        #multirempoints = [rpoint for rpoint in userempoints for tind2 in xrange(len(usedata[0][rpoint]))]
        #newyvals = convertForm(usepoints,usedata)
        #newremyvals = convertForm(userempoints,usedata)

        for yind,curyvals in enumerate(yvals):
            #cursumval,time2err,knots = runRFit(multipoints,curyvals,multirempoints,newremyvals[yind],usedata[yind])
            #modknots = changeRange(knots,times[0],times[-1])
            inyvals = [x2val[yind][rpoint] for rpoint in usepoints]
            inspl = scipy.interpolate.UnivariateSpline(usepoints, inyvals, s=reglambda, k=3)
            tempoutsplines.append(deepcopy(inspl))
            tsumval = weights[yind]*estQual(inspl,x2valblock[yind],usepoints,userempoints)
            globsumval += tsumval
            curknots.append(list(inspl.get_knots()))

            for tind,ttime in enumerate(usepoints):
                globsumval += np.std(list(usedata[yind][ttime]))**2
            
        print "inside: ",globsumval
        sol2val["_".join([str(titem) for titem in selpoints])] = globsumval
        if globsumval < sumval:
           sumval = globsumval
           minsol = list(selpoints)
           outsplines = list(tempoutsplines)
           y2knots = deepcopy(curknots)
    #print sol2val
    #print "done"
    #print y2knots
    
    if False:
     pts1 = [0.5,28.0]
     pts2 = [0.5,3.5,6.0,11.0, 16.0, 22.0, 28.0]
     pts3 = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0]
     pts4 = [0.5,7.0,28.0]

     def getArea(sentpoints,myyvals):
         """
         """
         locaream = 0.0
         for pind in xrange(len(sentpoints)-1):
            time1,time2 = sentpoints[pind:pind+2]
            y1,y2 = myyvals[pind:pind+2]
            if y1*y2 >= 0.0:
               locaream += abs(y2+y1)*(time2-time1)/2.0
            elif y1 <= 0.0 and y2 >= 0.0:
               ratio = (time2-time1)/float(y2-y1)
               breakp = time1+(ratio*abs(y1))
               locaream += (breakp-time1)*abs(y1)*0.5
               locaream += (time2-breakp)*abs(y2)*0.5
            elif y1 >= 0.0 and y2 <= 0.0:
               ratio = (time2-time1)/float(y1-y2)
               breakp = time1+(ratio*y1)
               locaream += (breakp-time1)*y1*0.5
               locaream += (time2-breakp)*abs(y2)*0.5
         return locaream
               
     sparea = 0.0
     linarea = 0.0
     y2lin,y2sp = {},{}
     y2lin1,y2lin2,y2lin3,y2lin4 = {},{},{},{}
     sampoints = np.arange(min(times),max(times)+0.01,0.01)
     points = sorted(list(minsol) + fixedpoints)
     for yind,curyvals in enumerate(yvals):
        inyvals = [x2val[yind][rpoint] for rpoint in points]
        inyvals1 = [x2val[yind][rpoint] for rpoint in pts1]
        inyvals2 = [x2val[yind][rpoint] for rpoint in pts2]
        inyvals3 = [x2val[yind][rpoint] for rpoint in pts3]
        inyvals4 = [x2val[yind][rpoint] for rpoint in pts4]
        loclinarea = getArea(points,inyvals)
        loclinarea1 = getArea(pts1,inyvals1)
        loclinarea2 = getArea(pts2,inyvals2)
        loclinarea3 = getArea(pts3,inyvals3)
        loclinarea4 = getArea(pts4,inyvals4)
        
        locsparea = 0.0    
        for samind in xrange(len(sampoints)-1):
            time1,time2 = sampoints[samind], sampoints[samind+1]
            tinspl = scipy.interpolate.UnivariateSpline(points, inyvals, s=reglambda, k=3)
            y1 = scipy.interpolate.splev(time1,tinspl._eval_args, der=0, ext=0)
            y2 = scipy.interpolate.splev(time2,tinspl._eval_args, der=0, ext=0)
            if y1*y2 < 0.0:
               print yind
            curaream = abs(y2+y1)*(time2-time1)/2.0
            if curaream >= 0.0:   
               locsparea += abs(y2+y1)*(time2-time1)/2.0
        y2sp[yind] = locsparea
        y2lin[yind] = loclinarea
        y2lin1[yind] = loclinarea1
        y2lin2[yind] = loclinarea2
        y2lin3[yind] = loclinarea3
        y2lin4[yind] = loclinarea4
        sparea += locsparea
        linarea += loclinarea
     difsp = {yind:abs(y2sp[yind]-y2lin[yind]) for yind in y2sp.keys()}
     difsp1 = {yind:abs(y2sp[yind]-y2lin1[yind]) for yind in y2sp.keys()}
     difsp2 = {yind:abs(y2sp[yind]-y2lin2[yind]) for yind in y2sp.keys()}
     difsp3 = {yind:abs(y2sp[yind]-y2lin3[yind]) for yind in y2sp.keys()}
     difsp4 = {yind:abs(y2sp[yind]-y2lin4[yind]) for yind in y2sp.keys()}
      
     import operator
     sorteddifs = sorted(difsp.items(), key=operator.itemgetter(1),reverse=True)
     sorteddifs1 = sorted(difsp1.items(), key=operator.itemgetter(1),reverse=True)
     sorteddifs2 = sorted(difsp2.items(), key=operator.itemgetter(1),reverse=True)
     sorteddifs3 = sorted(difsp3.items(), key=operator.itemgetter(1),reverse=True)
     sorteddifs4 = sorted(difsp4.items(), key=operator.itemgetter(1),reverse=True)

     #topfolder = "splineplots"
     #for fname in os.listdir(topfolder):
     #    if fname == ".DS_Store":
     #       continue 
     #    if len(fname.split("_")) == 3:
     #       continue
     #    newfname = "_".join(fname.split("_")[1:])
     #    code = "mv splineplots/{0} splineplots/{1}".format(fname, newfname)
     #    os.system(code)
             
     sortedfnames = []
     for myind,myval in sorteddifs4:
         gene = gene2ind[myind].lower()
         topfolder = "splineplots"
         flag = False
         for fname in os.listdir(topfolder):
             remgene = fname.split("_")[0].lower()
             if gene == remgene:
                flag = True
                foundfname = fname
                break
         if gene == "c/ebp":
            foundfname = "C,EBP_13_all.png"
            flag = True    
         if flag:   
            sortedfnames.append(foundfname)
         else:
            print gene
         #newfoundfname = "{0}_{1}".format(myind+1,foundfname)   
         #code = "mv splineplots/{0} splineplots/{1}".format(foundfname, newfoundfname)
         #os.system(code)
     print sortedfnames[0:4]
     sortedfnames.reverse()      
     base1,base2 = 1425712280, 1425712292       
     for bind,fname in enumerate(sortedfnames):
         fpath = "splineplots/{0}".format(fname)
         os.utime(fpath,(base1+bind,base2+bind))       
     print sortedfnames[0:4]
     sortedfnames.reverse() 
     for myind,fname in enumerate(sortedfnames):
         if fname.find("(")!=-1 and fname.find(")")!=-1:
            newfname = "{0}_{1}".format(myind+1,fname.split("(")[0]+fname.split(")")[1])
         else:
            newfname = "{0}_{1}".format(myind+1,fname) 
         code = "mv splineplots/{0} splineplots/{1}".format(fname, newfname)
         os.system(code)
     #print sorteddifs4
     print len(sorteddifs4)
     exit(1)
     
     print sparea
     print linarea
     if not os.path.exists("areadif"):
       os.makedirs("areadif")    
     with open("areadif/geneareadif.txt","w") as outfile:
        for sind,sval in sorteddifs:
            outfile.write("{0}\t{1}\n".format(gene2ind[sind],sval))
     with open("areadif/versus2linear.txt","w") as outfile:
        for sind,sval in sorteddifs1:
            outfile.write("{0}\t{1}\n".format(gene2ind[sind],sval))
     with open("areadif/versus7linear.txt","w") as outfile:
        for sind,sval in sorteddifs2:
            outfile.write("{0}\t{1}\n".format(gene2ind[sind],sval))
     with open("areadif/versus28linear.txt","w") as outfile:
        for sind,sval in sorteddifs3:
            outfile.write("{0}\t{1}\n".format(gene2ind[sind],sval))
     with open("areadif/select3.txt","w") as outfile:
        for sind,sval in sorteddifs4:
            outfile.write("{0}\t{1}\n".format(gene2ind[sind],sval))        
                            
    import operator
    sorted_data = sorted(sol2val.items(), key=operator.itemgetter(1))
    for keystr,optval in sorted_data[0:10]:
        print sorted([float(item) for item in keystr.split("_")]),optval,optval/float(len(yvals)*len(times))
    #exit(1)
     
    print "here"
    print minsol
    print globsumval    
    points = sorted(list(minsol) + fixedpoints)
    rempoints = list(set(times) - set(points))
    yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]  
    tsumval = 0.0
    for yind,curyvals in enumerate(yvals):
        inyvals = [x2val[yind][rpoint] for rpoint in points]
        tinspl = scipy.interpolate.UnivariateSpline(points, inyvals, s=reglambda, k=3)
        locsumval = weights[yind]*estQual(tinspl,x2valblock[yind],points,rempoints)
        tsumval += locsumval
    #print tsumval
    #print sumval    
    #assert abs(tsumval-sumval) < 0.0001
    print minsol
    print points
    print sumval
    print sumval/float(len(yvallist)*(len(times)-len(selpoints)-len(fixedpoints)))
    avglens = [len(list(tspl.get_knots())) for tspl in outsplines]
    print np.mean(avglens)
    #exit(1)
    return sumval,sumval/float(len(yvallist)*(len(times)-len(points))),points,yvals,y2knots,outsplines       




plotfolder = "splineplots"
if not os.path.exists(plotfolder):
   os.makedirs(plotfolder)
weightmode = "uni"
   
fname = "featureDataforCommonsites.csv"
gene2pos,pos2gene = featureData(fname)

print gene2pos.keys()
count = 0
for gene in gene2pos.keys():
    print gene,len(gene2pos[gene])
    count += len(gene2pos[gene])
print count    

fname = "pDataaug24.csv"            
name2time,time2name = readSampleTarget(fname)
times = sorted(list(set([float(item.replace("P","").replace("a","").replace("b","").replace("c","")) for item in time2name.keys() if item.startswith("P")])))
usetimes = list(times)

fname = "mergedData_CommonSites_Aug25.csv"
data,ind2name,ind2gene = readData(fname)

moddata = makeData(data,ind2gene,ind2name,name2time)

yvallist = deepcopy(moddata)
x2val = [{tval:np.median(yvallist[yind][tval]) for tind,tval in enumerate(usetimes)} for yind in xrange(len(yvallist))]
weights = [1.0]*len(yvallist)
sumval, avgsumval, points, yvals, y2knots, outsplines = runGreedy(times,moddata,yvallist,4,weights)

if True:
    rempoints = list(set(usetimes) - set(points))
     
    print "selected points are: "       
    print points

    avgyvallist = []
    for cind,cdata in enumerate(moddata):
        cavgdata = [np.median(list(moddata[cind][ttime])) for tind,ttime in enumerate(times)]
        avgyvallist.append(list(cavgdata))
           
    print "out: ",count,sumval,avgsumval            
    
    knotlens = {}
    knotarr = []
    for knots in y2knots:
        assert sorted(knots) == knots
        knotlens.setdefault(len(knots),0)
        knotlens[len(knots)] += 1
        knotarr.append(len(knots))
        
    #regs = [5.0,10.0,20.0,7.5,8.0,9.0,6.0,12.5,15.0,17.5,4.5,4.0,3.5,3.0,2.75,2.5,2.0,2.25,1.5,1.75,1.0,1.25,0.75,0.5,0.25]
    regs = [0.0001,0.01,0.1]
    mapval = {tval:tind for tind,tval in enumerate(usetimes)}
    allyvals = [[avgyvallist[yind][mapval[tpoint]] for tpoint in usetimes] for yind in xrange(len(yvallist))]

    avgdist = 0.0
    knotlens = {}
    for knots in y2knots:
        assert sorted(knots) == knots
        curavgdist = 0.0
        for ind1 in xrange(len(knots)-1):
            curavgdist += abs(knots[ind1+1]-knots[ind1])
        avgdist += curavgdist/float(len(knots)-1)
        knotlens.setdefault(len(knots),0)
        knotlens[len(knots)] += 1
    print "knot len distribution"
    print knotlens
    
    if True:
     print "plotting starts"

     collects = []
     for gind,youts in enumerate(yvals):
         outspl = outsplines[gind]
         tyvals = [item for item in youts] + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in y2knots[gind]] + [x2val[gind][rpoint] for rpoint in rempoints] #+ [scipy.interpolate.splev(ttime,outspl._eval_args, der=0, ext=0) for ttime in times]
         collects.extend(tyvals)
     print min(collects), max(collects)
     #exit(1)
         
     for gind,youts in enumerate(yvals):
        foundlambda, foundknots,foundspl = None, None, None
        mindifval = 1000.0
        minrealval = None
        for treglambda in regs:
            allspl = scipy.interpolate.UnivariateSpline(usetimes, allyvals[gind], s=treglambda, k=3)
            inferknots = allspl.get_knots()
            if len(inferknots) == knotarr[gind]:
               foundlambda = treglambda
               foundknots = list(inferknots)
               foundspl = deepcopy(allspl)
               break
            difval = abs(len(inferknots)-knotarr[gind])
            if difval < mindifval:
               mindifval = difval
               minrealval = len(inferknots)
        if foundlambda != None:
           print ind2gene[gind], len(foundknots),foundlambda
        else:
           print ind2gene[gind], knotarr[gind], mindifval, minrealval

        yvaldictout = {"SplineFit": list(youts)}
        remyvals = [x2val[gind][rpoint] for rpoint in rempoints]
        gene = ind2gene[gind]
        plotpath = "{0}/{1}/{2}_{3}".format(plotfolder,weightmode,gene.replace("/",","),len(points))
        print plotpath
        if os.path.exists(plotpath+".png"):
           print "not writing" 
           continue


        trixvals = []
        triyvals = []
        #trixvals = [0.5,7.0,28.0]
        #triyvals = [avgyvallist[gind][mapval[trixval]] for trixval in trixvals]

        if len(y2knots[gind]) > 4:
           continue 
        makeplotMain(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,plotpath,usetimes,trixvals,triyvals)
        #makeplotMain(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,plotpath,usetimes,(min(collects),max(collects)))
           
        allplotpath = "{0}/{1}_{2}_all".format(plotfolder,gene.replace("/",","),len(points))
        if os.path.exists(allplotpath+".png"):
           continue 
        if foundlambda != None:
           makeplotMainAll(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,allplotpath,foundknots,foundspl,usetimes,trixvals,triyvals)


exit(1)



print len(moddata)
print len(pos2gene)



logmoddata = transformData(moddata)
allvals = [item for rowdata in logmoddata for item in rowdata]
print max(allvals)
print np.mean(allvals)

allvals = [item for rowdata in moddata for ttime in rowdata.keys() for item in rowdata[ttime]]
print max(allvals)
print np.mean(allvals)

exit(1)




allvals = [item for rowdata in data for item in rowdata]
print max(allvals)
print min(allvals)
print np.mean(allvals)
exit(1)


#gene matching
allgenes = []
with open("mrnagenes.txt","r") as infile:
    for line in infile:
        line = line.rstrip()
        allgenes.append(line.lower())       
for tgene in gene2pos.keys():
    tgene = tgene.lower()
    found = []
    for tgene2 in allgenes:
        if tgene2.find(tgene)!=-1 or tgene.find(tgene2)!=-1:
           found.append(tgene2) 
    print tgene,found

#for tgene in allgenes:
#    if tgene.find("zfp")!=-1 or tgene.find("Kiaa")!=-1 or tgene.find("zf")!=-1:
#       print tgene
        
gene2found = chroAnalyze(ind2gene)
print gene2found
exit(1)
print pos2gene


