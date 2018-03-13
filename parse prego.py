# Copyright (c) 2018, Ron Zeira, Tel Aviv University
# All rights reserved.

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from test_graph_ilp import *

###############################
# Main function for reading the karyotypes of PREGO and analyzing it
def readIntervalGraphFile(fileName,interesting = range(1,24)):
    print fileName
    f = open(fileName,'r')
    intrevalEdges = {}
    refEdges = {}
    variantEdges = {}
    nodes = set()
    lines = f.readlines()
    f.close()
    for i in range(len(lines)):
        lines[i] = lines[i].rstrip()
    intSection = False
    refSection = False
    varSection = False
    for line in lines:
        if line.startswith("interval edges:"):
            intSection = True
            continue
        if line.startswith("ref edges:"):
            intSection = False
            refSection = True
            continue
        if line.startswith("variant edges:"):
            refSection = False
            varSection = True
            continue
        if line.startswith("source edges:"):
            varSection = False
            continue
        if line == "":
            continue
        l = line.split()
        first = l[0]
        if len(l)==7:
            second = l[3]
            val = l[6]
        if len(l)==3:
            second = l[1]
            val = l[2]
        if intSection:
            if int(val)>0:
                intrevalEdges[first,second] = int(val)
        if refSection:
            if int(val)>0:
                refEdges[first,second] = int(val)
        if varSection:
            if int(val)>0:
                variantEdges[first,second] = int(val)
        if int(val)>0:
            nodes.add(first)
            nodes.add(second)

    #print nodes
    #print intrevalEdges
    #print refEdges
    #print variantEdges
    G=nx.Graph()
    for u,v in intrevalEdges.keys():
        if intrevalEdges[u,v]>0:
            if G.has_edge(u,v):
                G.add_edge(u,v,weight=intrevalEdges[u,v]+G.get_edge_data(u,v)['weight'])
            else:
                G.add_edge(u,v,weight=intrevalEdges[u,v])
    for u,v in refEdges.keys():
        if refEdges[u,v]>0:
            if G.has_edge(u,v):
                G.add_edge(u,v,weight=refEdges[u,v]+G.get_edge_data(u,v)['weight'])
            else:
                G.add_edge(u,v,weight=refEdges[u,v])
    for u,v in variantEdges.keys():
        if variantEdges[u,v]>0:
            if G.has_edge(u,v):
                G.add_edge(u,v,weight=variantEdges[u,v]+G.get_edge_data(u,v)['weight'])
            else:
                G.add_edge(u,v,weight=variantEdges[u,v])

    # cc = nx.connected_component_subgraphs(G)
    # cc = []
    # #print cc
    # chromDist = {}
    # for c in cc:
    #     if len(c)>2:
    #         chroms = getChromosomes(c)
    #         if len(chroms)==1:
    #             #if chroms[0]!=11:
    #             #    continue
    #             print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    #             print "Chromosome:", chroms
    #             cn,adj = parseChromosome(chroms[0],intrevalEdges,refEdges,variantEdges)
    #             chromCont = chromosomeContainer()
    #             chromCont.setMatrices(cn,adj)
    #             print "Original:"
    #             print chromCont
    #             print "Compressing:"
    #             chromCont.compressChromosome()
    #             print chromCont
    #             print "Removing TD:"
    #             chromCont.compressTandemDuplications()
    #             print chromCont
    #             print "-----------------------------------------------------------"
    #             d = searchDistance(createDiploidChrome(chromCont.n),chromCont,10)
    #             chromDist[chroms[0]] = d
    #             #intervalGraphToMatrixOneChromosome(c)
    #             #plotIntervalGraph(c,intersectEdgeLists(intrevalEdges.keys(),c.edges()),intersectEdgeLists(refEdges.keys(),c.edges()),intersectEdgeLists(variantEdges.keys(),c.edges()))
    # #for chrom in range(1,23):
    # #    cn,adj = parseChromosome(chrom,intrevalEdges,refEdges,variantEdges)
    all = chromCC(variantEdges)
    chromDist = {}
    for chromList in all:
        if len(set(chromList).intersection(set(interesting)))==0:
            continue
        print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        chromCont = parseChromList(chromList,intrevalEdges,refEdges,variantEdges)
        print "Compressing:"
        chromCont.compressChromosome()
        print chromCont
        td=0
        print "Removing TD:"
        td = chromCont.compressTandemDuplications()
        print chromCont
        print "Compressing 2:"
        chromCont.compressChromosome()
        print chromCont
        print "-----------------------------------------------------------"


        chr = removeLongestDuplicatedChromosomes(chromCont)
        # chr_dip = chr.getDiploid(False,1)
        # print chr_dip
        # dd= searchDistance(chr_dip,chr,15,graphFileName="OV"+fileName[fileName.find("OV")+2]+"_Chr"+".".join([str(x) for x in chromList])+"BFB.graphml")
        # print dd
        # return
        if isinstance(chr,chromosomeContainer):
            i=0
            while (chromCont.cn[chromCont.telGenes]>=2).all():
                chromCont -= chr
                if (chromCont.cn[chromCont.telGenes]<2).any():
                    chromCont += chr
                    break
                i+=1
        print "removed ",i," chromosomes"
        print chromCont


        #print chromContOrg
        #return searchDistance(chromContOrg,chromCont,16,0,False,[],graphFileName="temp.graphml")

        dip = chromCont.getDiploid()
        lb = np.sum(np.logical_and(dip.adj==0,chromCont.adj>0))/4
        print "LB (BP/2) = ", lb
        impossible = impossiblePairs(dip,chromCont)
        graphFile = "OV"+fileName[fileName.find("OV")+2]+"_Chr"+".".join([str(x) for x in chromList])+".graphml"
        d = searchDistance(dip,chromCont,15,lb,False,impossible,graphFileName=graphFile)
        if d>=0:
            chromDist[tuple(chromList)] = d+td
        else:
            chromDist[tuple(chromList)] = d
    print "======================================================================"
    for chrom in chromDist:
        print "Chromosome ", chrom, ", distance ", chromDist[chrom]

def chromCC(variantEdges):
    cc = []
    d = {}
    for chrom in range(1,23):
        d[chrom] = set([chrom])
    for e1,e2 in variantEdges:
        c1 = getChrome(e1)
        c2 = getChrome(e2)
        s1 = d[c1]
        s2 = d[c2]
        union = s1.union(s2)
        for c3 in union:
            d[c3] = union
    all = set()
    for chrom in range(1,23):
        all.add(tuple(sorted(list(d[chrom]))))
        #print chrom, d[chrom]
    all = list(all)
    all = [list(t) for t in all]
    #print all
    return all

def parseChromList(chromLst,intrevalEdges,refEdges,variantEdges):
    print "Chromosomes:", chromLst
    chromCont = chromosomeContainer()
    offsets = []
    for chrom in chromLst:
        cn,adj = parseChromosome(chrom,intrevalEdges,refEdges,variantEdges)
        chromContNew = chromosomeContainer()
        chromContNew.setMatrices(cn,adj)
        off = chromCont.mergeKar(chromContNew)
        #print off
        offsets.append(off)
    if len(chromLst)>1:
        extChromList = []
        for chrom in chromLst:
            extChromList.append(getChromExt(intrevalEdges,chrom))
        otherEdges = crossChromosomeEdges(chromLst,extChromList,variantEdges)
        print "other edges",otherEdges
        newEdges = []
        for n1,c1,n2,c2,count in otherEdges:
            n1 += offsets[chromLst.index(c1)]
            n2 += offsets[chromLst.index(c2)]
            newEdges.append((n1,n2,count))
        print "new edges",newEdges
        chromCont.addAdjList(newEdges)
    print "Karyotype"
    print chromCont
    return chromCont

def crossChromosomeEdges(chromLst,extChromList,variantEdges):
    otherEdges = []
    for e1,e2 in variantEdges:
        if getChrome(e1)!=getChrome(e2) and getChrome(e1) in chromLst and getChrome(e2) in chromLst:
            c1 = chromLst.index(getChrome(e1))
            c2 = chromLst.index(getChrome(e2))
            i1 = extChromList[c1].index(e1)+2
            i2 = extChromList[c2].index(e2)+2
            otherEdges.append((i1,getChrome(e1),i2,getChrome(e2),variantEdges[e1,e2]))
    return otherEdges

def parseChromosome(chrom,intrevalEdges,refEdges,variantEdges):
    segments = getChromSegments(intrevalEdges,chrom)
    print "segments:", segments
    extremities = getChromExt(intrevalEdges,chrom)
    print "extremities:",extremities
    cn = np.zeros(len(segments)+2,dtype=int)
    adjacencies = np.zeros((len(extremities)+4,len(extremities)+4),dtype=int)
    for s in segments:
        a,b = getSegmentsExtFromList(s,extremities)
        cn[segments.index(s)+1] = intrevalEdges[a,b]
    cn[0] = cn[1]
    cn[-1] = cn[-2]
    #print cn
    for e1 in extremities:
        for e2 in extremities:
            #if (e1,e2) in refEdges:
            #    print "ref edge ",(e1,e2),refEdges[e1,e2]
            #if (e1,e2) in variantEdges:
            #    print "variant edge ",(e1,e2),variantEdges[e1,e2]
            adjacencies[extremities.index(e1)+2,extremities.index(e2)+2] = refEdges.get((e1,e2),refEdges.get((e2,e1),0))+variantEdges.get((e1,e2),variantEdges.get((e2,e1),0))
    adjacencies[1,2] = cn[0]
    adjacencies[2,1] = cn[0]
    adjacencies[-2,-3] = cn[-1]
    adjacencies[-3,-2] = cn[-1]
    ## inter chromosome connections
    otherEdges = []
    for e1,e2 in variantEdges:
        if (getChrome(e1)==chrom and getChrome(e2)!=chrom) or (getChrome(e1)!=chrom and getChrome(e2)==chrom):
            if  getChrome(e1)==chrom:
                p = (e1,e2,variantEdges[e1,e2])
            if getChrome(e2)==chrom:
                p = (e2,e1,variantEdges[e1,e2])
            otherEdges.append(p)
    #print "Other edges: ", otherEdges
    #print adjacencies
    return cn,adjacencies

def getChromSegments(intrevalEdges,chrom):
    s = set()
    for e1,e2 in intrevalEdges.keys():
        if getChrome(e1)==chrom:
            s.add(getSegment(e1))
    l = list(s)
    l.sort()
    return l

def getChromExt(intrevalEdges,chrom):
    s = set()
    for e1,e2 in intrevalEdges.keys():
        if getChrome(e1)==chrom:
            s.add(e1)
        if getChrome(e1)==chrom:
            s.add(e2)
    l = list(s)
    l.sort(key=parseNodeName)
    return l

def getSegmentsExtFromList(s,l):
    ret = []
    for e in l:
        if getSegment(e)==s:
            ret.append(e)
    ret.sort()
    return tuple(ret)

def intersectEdgeLists(l1,l2):
    l = []
    for u,v in l1:
        if (u,v) in l2:
            l.append((u,v))
        if (v,u) in l2:
            l.append((v,u))
    return l

def parseNodeName(s):
    l = s.split('_')
    chrome = int(l[1])
    start = int(l[2])
    segment = int(l[0][:-1])
    side = l[0][-1]
    return segment,side,chrome,start

def getChrome(s):
    return parseNodeName(s)[2]

def getSegment(s):
    return parseNodeName(s)[0]

def getChromosomes(g):
    s = set()
    for v in g:
        s.add(getChrome(v))
    return list(s)

def intervalGraphToMatrixOneChromosome(g):
    node2data = {}
    for v in g:
        node2data[v] = parseNodeName(v)
    nodes = node2data.keys()
    nodes.sort(key=lambda v: node2data[v][3])
    print nodes
    genes = []
    genes2pos = {}
    for v in nodes:
        segment,side,chrome,start = node2data[v]
        if not segment in genes:
            genes.append(segment)
            genes2pos[segment] = start
        else:
            genes2pos[segment] = (genes2pos[segment],start)
    print genes
    print genes2pos
    extremities = []
    extremities2genepos = {}
    for v in nodes:
        segment,side,chrome,start = node2data[v]
        extremities.append((segment,side))
        extremities2genepos[segment,side] = start
    cn = np.zeros(len(genes),dtype=int)
    adj = np.zeros((len(extremities),len(extremities)),dtype=int)
    for u,v in g.edges():
        segment1,side1 = node2data[u][0:2]
        segment2,side2 = node2data[v][0:2]
        if segment1==segment2 and side1!=side2:
            i = genes.index(segment1)
            cn[i] += g.get_edge_data(u,v)['weight']
        else:
            i = extremities.index((segment1,side1))
            j = extremities.index((segment2,side2))
            adj[i,j] += g.get_edge_data(u,v)['weight']
            adj[j,i] = adj[i,j]
    print cn
    print extremities
    print extremities2genepos
    print adj

def plotIntervalGraph(G,intrevalEdges,refEdges,variantEdges):
    plt.figure()
    pos=nx.spring_layout(G)
    #pos=nx.circular_layout(G)
    nx.draw_networkx_nodes(G,pos,node_size=700)
    #nx.draw_networkx_edges(G,pos,width=6)
    nx.draw_networkx_edges(G,pos,edgelist=intrevalEdges,width=6,edge_color='r')
    nx.draw_networkx_edges(G,pos,edgelist=refEdges,width=6)
    nx.draw_networkx_edges(G,pos,edgelist=variantEdges,width=6,edge_color='b')
    nx.draw_networkx_labels(G,pos,font_size=10,font_family='sans-serif')
    plt.axis('off')
    plt.show()


#for paper
#readIntervalGraphFile("PREGO/OV5_MIN10_ALL.coverage.final.merged.noXY.results",[2,17])
#readIntervalGraphFile("PREGO/OV1_MIN10_ALL.coverage.final.merged.noXY.results",[11,20])
#readIntervalGraphFile("PREGO/OV2_MIN10_ALL.coverage.final.merged.noXY.results",[12,16,18])

#test
readIntervalGraphFile("PREGO/OV5_MIN10_ALL.coverage.final.merged.noXY.results",range(1,24))