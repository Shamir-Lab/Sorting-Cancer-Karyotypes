# Copyright (c) 2018, Ron Zeira, Tel Aviv University
# All rights reserved.

from gurobipy import *
import numpy as np
import multiprocessing

####

def genesToExtremities(genes):
    l = []
    for i in genes:
        l.append(getGeneTail(i))
        l.append(getGeneHead(i))
    return l

###########
# representation as a function of n

def getNfromChrom(chrom):
    return max(np.abs(chrom))

def getNfromChromList(chromList):
    return max([getNfromChrom(chrom) for chrom in chromList])

def nToGenesOnly(n):
    return range(1,n+1)

def leftTelomereGene(firstGene=1):
    return firstGene-1

def rightTelomereGene(n):
    return n+1

def nToTelomereGenes(n):
    return [leftTelomereGene(),rightTelomereGene(n)]

def nToGenesAndTels(n):
    genes = nToGenesOnly(n)
    return [leftTelomereGene()]+genes+[rightTelomereGene(n)]

def nToNwithTels(n):
    return n+2

def leftTelomereExt():
    return getGeneRightExt(leftTelomereGene())

def rightTelomereExt(n):
    return getGeneLeftExt(rightTelomereGene(n))

def nToTelomereExtremities(n):
    return [leftTelomereExt(),rightTelomereExt(n)]

def nToExtremitiesOnly(n):
    return genesToExtremities(nToGenesOnly(n))

def nToExtremitiesAndTels(n):
    ext = nToExtremitiesOnly(n)
    return [leftTelomereExt()]+ext+[rightTelomereExt(n)]

def nToNumberExtWithTels(n):
    return 2*n+4


#####
# gene/extremity conversion

def getGeneTail(gene):
    return 2*gene

def getGeneHead(gene):
    return 2*gene+1

def getGeneLeftExt(gene):
    if gene>=0:
        return getGeneTail(abs(gene))
    else:
        return getGeneHead(abs(gene))

def getGeneRightExt(gene):
    if gene>=0:
        return getGeneHead(abs(gene))
    else:
        return getGeneTail(abs(gene))

def getOtherExt(ext):
    if ext%2:
        return ext-1
    else:
        return ext+1

def getGeneFromExt(ext):
    if ext%2:
        return (ext-1)/2
    else:
        return ext/2

####
#printing funcs

def getExtSideStr(ext):
    if ext%2:
        return "h"
    else:
        return "t"

def extToGeneSide(ext):
    return str(getGeneFromExt(ext))+getExtSideStr(ext)

def adjToTuple(adj):
    return extToGeneSide(adj[0]),extToGeneSide(adj[1])

def adjListToTupleList(l):
    return [adjToTuple(adj) for adj in l]

#####
#new conversion funcs

def chromToCn(chrom,n=None):
    if n is None:
        n = getNfromChrom(chrom)
    cn = np.zeros(nToNwithTels(n),dtype=int)
    for gene in np.abs(chrom):
        cn[gene] += 1
    cn[leftTelomereGene()] = 1
    cn[rightTelomereGene(n)] = 1
    return cn

def chromListToCn(chromList,n=None):
    if n is None:
        n = getNfromChromList(chromList)
    cn = np.zeros(nToNwithTels(n),dtype=int)
    for chrom in chromList:
        cn += chromToCn(chrom,n)
    return cn

def chromToAdjMat(chrom,n=None):
    if n is None:
        n = getNfromChrom(chrom)
    mat = np.zeros((nToNumberExtWithTels(n),nToNumberExtWithTels(n)),dtype=int)
    for i in range(len(chrom)-1):
        leftGene = chrom[i]
        rightGene = chrom[i+1]
        leftExt = getGeneRightExt(leftGene)
        rightExt = getGeneLeftExt(rightGene)
        if leftExt!=rightExt:
            mat[leftExt,rightExt] += 1
            mat[rightExt,leftExt] += 1
        else:
            mat[leftExt,rightExt] += 1
    mat[leftTelomereExt(),getGeneLeftExt(chrom[0])] = 1
    mat[getGeneLeftExt(chrom[0]),leftTelomereExt()] = 1
    mat[rightTelomereExt(n),getGeneRightExt(chrom[-1])] = 1
    mat[getGeneRightExt(chrom[-1]),rightTelomereExt(n)] = 1
    return mat

def chromListToAdjMat(chromList,n=None):
    if n is None:
        n = getNfromChromList(chromList)
    mat = np.zeros((nToNumberExtWithTels(n),nToNumberExtWithTels(n)),dtype=int)
    for chrom in chromList:
        mat += chromToAdjMat(chrom,n)
    return mat

def adjListToMat(adjList,m):
    mat = np.zeros((m,m),dtype=int)
    for tup in adjList:
        mat[tup[0],tup[1]] = tup[2]
        mat[tup[1],tup[0]] = tup[2]
    return mat

def increaseIntListBy(l,d):
    return list(np.array(l,dtype=int)+d)

################################################
## A class representing karyotypes

class chromosomeContainer:
    def __init__(self,karyotype=None,n=None):
        if karyotype:
            self.karyotype = karyotype[:]
            if n is None:
                self.n = getNfromChromList(karyotype)
            else:
                self.n = n
            self.genes = nToGenesOnly(self.n)
            self.extremities = nToExtremitiesOnly(self.n)

            self.telExt = nToTelomereExtremities(self.n)
            self.telGenes = nToTelomereGenes(self.n)
            self.extAndTels = nToExtremitiesAndTels(self.n)
            self.geneAndTels = nToGenesAndTels(self.n)

            self.cn = chromListToCn(karyotype,self.n)
            self.adj = chromListToAdjMat(karyotype,self.n)
        else:
            self.n = -2

    def setMatrices(self,cn,adj,telGenes=None):
        if not telGenes:
            telGenes = [0,len(cn)-1]
        self.n = len(cn) - len(telGenes)
        self.karyotype = None
        self.genes = list(set(range(len(cn))).difference(set(telGenes)))
        self.extremities = genesToExtremities(self.genes)

        self.telExt = sorted([getGeneHead(telGenes[i]) for i in range(len(telGenes)) if i%2==0]+[getGeneTail(telGenes[i]) for i in range(len(telGenes)) if i%2==1])
        self.telGenes = telGenes[:]
        self.geneAndTels = range(len(cn))
        self.extAndTels = sorted(self.telExt+self.extremities)

        self.cn = cn
        self.adj = adj

    def mergeKar(self,other):
        self.karyotype = None

        if self.n>0:
            geneOffset = nToNwithTels(self.n)
            extOffset = nToNumberExtWithTels(self.n)
            self.genes += increaseIntListBy(other.genes,geneOffset)
            self.extremities += increaseIntListBy(other.extremities,extOffset)
            self.telExt += increaseIntListBy(other.telExt,extOffset)
            self.telGenes += increaseIntListBy(other.telGenes,geneOffset)
            self.extAndTels += increaseIntListBy(other.extAndTels,extOffset)
            self.geneAndTels += increaseIntListBy(other.geneAndTels,geneOffset)
            cnTmp = np.zeros(len(self.geneAndTels),dtype=int)
            adjTmp = np.zeros((self.adj.shape[0]+other.adj.shape[0],self.adj.shape[1]+other.adj.shape[1]),dtype=int)
            cnTmp[:len(self.cn)] = self.cn
            cnTmp[len(self.cn):] = other.cn
            adjTmp[:self.adj.shape[0],:self.adj.shape[1]] = self.adj
            adjTmp[self.adj.shape[0]:,self.adj.shape[1]:] = other.adj
            self.cn = cnTmp
            self.adj = adjTmp
            self.n = other.n+geneOffset
            return extOffset
        else:
            self.genes = other.genes[:]
            self.extremities = other.extremities[:]
            self.telExt = other.telExt[:]
            self.telGenes = other.telGenes[:]
            self.extAndTels = other.extAndTels[:]
            self.geneAndTels = other.geneAndTels[:]
            self.cn = other.cn.copy()
            self.adj = other.adj.copy()
            self.n = other.n
            return 0

    def __len__(self):
        return self.n
    def __repr__(self):
        s = ""
        if self.karyotype:
            s += "Karyotype: " + str(self.karyotype) + "\n"
        s += "Genes: " + str(self.genes) + "\n"
        s += "Extremities: " + str(self.extremities) + "\n"
        s += "Telomere Extremities: "+ str(self.telExt) + "\n"
        s += "cn: " + str(self.cn) + "\n"
        #s += "adj: " + str(self.adj) + "\n"
        s += "adj: " + str(adjMatToAdjList(self.adj)) + "\n"
        return s
    def compressChromosome(self):
        newCn,newAdj, newTels = compressKaryotypeMat(self.cn,self.adj,self.telGenes)
        self.setMatrices(newCn,newAdj,newTels)
    def addAdjList(self,lst):
        for tup in lst:
            if tup[0]!=tup[1]:
                self.adj[tup[0],tup[1]]+=tup[2]
                self.adj[tup[1],tup[0]]+=tup[2]
            else:
                self.adj[tup[0],tup[1]]+=tup[2]
    def changePath(self,geneList,delta):
        self.cn[geneList] += delta
        for i in range(len(geneList)-1):
            self.adj[getGeneHead(geneList[i]),getGenetail(geneList[i]+1)] += delta
            self.adj[getGenetail(geneList[i]+1),getGeneHead(geneList[i])] += delta
    def compressTandemDuplications(self):
        newCn,newAdj,td = compressTandemDulications(self.cn,self.adj)
        self.setMatrices(newCn,newAdj,self.telGenes)
        return td
    def getDiploid(self,allChroms=True,ploidy=2):
        dip = chromosomeContainer()
        cnTmp = ploidy*np.ones(len(self.geneAndTels),dtype=int)
        adjTmp = np.zeros((self.adj.shape[0],self.adj.shape[1]),dtype=int)
        for gene in self.genes:
            previous = gene-1
            next = gene + 1
            adjTmp[getGeneHead(previous),getGeneTail(gene)] = ploidy
            adjTmp[getGeneTail(gene),getGeneHead(previous)] = ploidy
            adjTmp[getGeneHead(gene),getGeneTail(next)] = ploidy
            adjTmp[getGeneTail(next),getGeneHead(gene)] = ploidy
        if allChroms:
            dip.setMatrices(cnTmp,adjTmp,self.telGenes)
            return dip
        for i in range(0,len(self.telGenes),2):
            left,right = self.telGenes[i],self.telGenes[i+1]
            if (self.cn[left:right+1]==0).all():
                cnTmp[left:right+1]=0
                adjTmp[getGeneTail(left):getGeneHead(right)+1,getGeneTail(left):getGeneHead(right)+1] = 0
        dip.setMatrices(cnTmp,adjTmp,self.telGenes)
        return dip
    def __sub__(self, other):
        new = chromosomeContainer()
        new.setMatrices(self.cn-other.cn,self.adj-other.adj,self.telGenes)
        return new
    def __add__(self, other):
        new = chromosomeContainer()
        new.setMatrices(self.cn+other.cn,self.adj+other.adj,self.telGenes)
        return new
    def __mul__(self, other):
        new = chromosomeContainer()
        new.setMatrices(other*self.cn,other*self.adj,self.telGenes)
        return new
    def __rmul__(self, other):
        return self*other
    def __eq__(self, other):
        return (self.cn==other.cn).all() and (self.adj==other.adj).all()
    def __le__(self, other):
        return (self.cn<=other.cn).all() and (self.adj<=other.adj).all()
    def __ge__(self, other):
        return (self.cn>=other.cn).all() and (self.adj>=other.adj).all()
    def __lt__(self, other):
        return (self.cn<other.cn).all() and (self.adj<other.adj).all()
    def __gt__(self, other):
        return (self.cn>other.cn).all() and (self.adj>other.adj).all()
    def __contains__(self, item):
        return self>=item
    def has_distance_to(self,other):
        return ((self.cn>0)>=(other.cn>0)).all()
    def copy(self):
        new = chromosomeContainer()
        new.setMatrices(self.cn,self.adj,self.telGenes)
        return new
    def printGraphML(self,fileName,genePath=[],edgePath=[],start=True,end=True,vert=0,event=""):
        nodeColors = ["#FFC0CB","#FF1493","#9400D3","#7A67EE","#6495ED","#C6E2FF","#1E90FF","#00B2EE","#00E5EE","#00FA9A","#43CD80","#00FF00","#6B8E23","#C0FF3E","#CDCDC1","#FAFAD2","#FFFF00","#FFD700","#B8860B","#FFA500","#FF6103","#B22222","#FF3030","#FFC0FF"]
        if start:
            f = open(fileName,'w')
            f.write(getHeader())
        else:
            f = open(fileName,'a')

        #Iterate through and draw nodes
        pos = 0
        #lastChrm = 0
        chrm = 1
        numTels = 0
        for aNode in self.extAndTels:
            name = extToGeneSide(aNode)
            #if lastChrm != chrm:
            #    pos = 0
            kind = ""
            if aNode in self.extremities:
                kind = "internal"
            if aNode in self.telExt:
                kind = "telomeric"
                if numTels==0:
                    numTels +=1
                    chrm +=1
                    pos = 0
                else:
                    numTels = 0
            nodeGraphML = getNodeGraphML(chrm,name,kind,pos,getGeneFromExt(aNode) in genePath,vert)
            f.write(nodeGraphML)
            pos += 1
            #lastChrm = chrm

        if not end:
            nodeGraphML = getNodeGraphML(chrm+1,event,event,0,False,vert)
            f.write(nodeGraphML)

        #Iterate through interval edges
        for gene in self.genes:
            n1 = extToGeneSide(getGeneTail(gene))
            n2 = extToGeneSide(getGeneHead(gene))
            edgeGraphML = drawEdgeGraphML(n1,n2,self.cn[gene],"interval",gene in genePath,vert);
            f.write(edgeGraphML)

        #Iterate through adjacencies
        for e1 in self.extAndTels:
            for e2 in self.extAndTels:
                if e1 > e2:
                    continue
                if self.adj[e1,e2]==0:
                    continue
                n1 = extToGeneSide(e1)
                n2 = extToGeneSide(e2)
                if e1+1==e2 and getGeneFromExt(e1)!=getGeneFromExt(e2):
                    edgeGraphML = drawEdgeGraphML(n1,n2,self.adj[e1,e2],"reference",(e1,e2) in edgePath or (e2,e1) in edgePath,vert)
                else:
                    edgeGraphML = drawEdgeGraphML(n1,n2,self.adj[e1,e2],"variant",(e1,e2) in edgePath or (e2,e1) in edgePath,vert)
                f.write(edgeGraphML)

        if end:
            f.write(getFooter())

        f.close()

    def ilpResultsToCnAdj(self,cnDic,adjDic,k):
        cn = self.cn.copy()
        adj = self.adj.copy()
        for i in range(len(cn)):
            if hasattr(cnDic[k,i],'x'):
                cn[i] = cnDic[k,i].x
            else:
                cn[i] = cnDic[k,i]
        for i in range(len(adj)):
            for j in range(len(adj)):
                if (k,i,j) in adjDic:
                    if hasattr(adjDic[k,i,j],'x'):
                        adj[i,j] = adjDic[k,i,j].x
                        adj[j,i] = adjDic[k,i,j].x
                    else:
                        adj[i,j] = adjDic[k,i,j]
                        adj[j,i] = adjDic[k,i,j]
                elif (k,j,i) in adjDic:
                    if hasattr(adjDic[k,j,i],'x'):
                        adj[i,j] = adjDic[k,j,i].x
                        adj[j,i] = adjDic[k,j,i].x
                    else:
                        adj[i,j] = adjDic[k,j,i]
                        adj[j,i] = adjDic[k,j,i]
        new = chromosomeContainer()
        new.setMatrices(cn,adj,self.telGenes)
        return new


def adjMatToAdjList(mat):
    n = len(mat)
    adj = []
    for i in range(n):
        for j in range(i,n):
            if mat[i,j]>0:
                adj.append((i,j,mat[i,j]))
    return adj

def createDiploidChrome(n):
    return chromosomeContainer([range(1,n+1),range(1,n+1)],n)

def adjMatToLists(mat):
    n = len(mat)
    adj = [[] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if mat[i,j]>0:
                adj[i].append(j)
    return adj

###################################################

def getHeader():
    sb = []
    sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    sb.append("<!-- Test file -->\n")
    sb.append("<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns/graphml\"\n")
    sb.append("\txmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n\t")
    sb.append("xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/graphml\n")
    sb.append("\thttp://www.yworks.com/xml/schema/graphml/1.0/ygraphml.xsd\"\n")
    sb.append("\txmlns:y=\"http://www.yworks.com/xml/graphml\">\n")
    sb.append("\t<key id=\"d0\" for=\"node\" yfiles.type=\"nodegraphics\"/>\n")
    sb.append("\t<key id=\"d1\" for=\"node\" attr.name=\"BooleanValue\" ")
    sb.append("attr.type=\"boolean\"/>\n")
    sb.append("\t<key id=\"d2\" for=\"edge\" yfiles.type=\"edgegraphics\"/>\n")
    sb.append("\t<key id=\"d3\" for=\"edge\" attr.name=\"IntValue\" attr.type=")
    sb.append("\"int\"/>\n")
    sb.append("\t<key id=\"d4\" for=\"graph\" attr.name=\"StringValue\"")
    sb.append(" attr.type=\"string\"/>\n")
    sb.append("\t<graph id=\"Main\" edgedefault=\"directed\">\n")
    return "".join(sb)

def getNodeGraphML(chrm,name,kind,pos,dashed=False,vert=0):
    nodeColors =  ["#FFC0CB","#FF1493","#6495ED","#00FA9A","#FFC0FF"] #["#FFC0CB","#FF1493","#9400D3","#7A67EE","#6495ED","#C6E2FF","#1E90FF","#00B2EE","#00E5EE","#00FA9A","#43CD80","#00FF00","#6B8E23","#C0FF3E","#CDCDC1","#FAFAD2","#FFFF00","#FFD700","#B8860B","#FFA500","#FF6103","#B22222","#FF3030","#FFC0FF"]
    if kind == "internal":
        shape = "rectangle"
    elif kind == "telomeric":
        shape = "trapezoid"
    else:
        shape = "triangle2"
        #shape = "oval"
    sb = []
    nodeLabel = name
    nodeName = "nodeLoc" + nodeLabel +":"+ str(vert)
    sb.append("\t\t<data key=\"d4\">nodes</data>\n")
    sb.append("\t\t<node id=\""+ nodeName  + "\">\n")
    sb.append("\t\t\t<data key=\"d0\">\n")
    sb.append("\t\t\t\t<y:ShapeNode>\n")
    sb.append("\t\t\t\t\t<y:Geometry x=\"")
    if kind == "internal" or kind == "telomeric":
        #sb.append(str(pos*200));
        sb.append(str(pos*300))
    else:
        #sb.append(str(200))
        sb.append(str(300))
    sb.append(".0\" y=\"")
    #sb.append(str((vert+chrm)*100))
    sb.append(str((vert+chrm)*100))
    if kind == "internal" or kind == "telomeric":
        #sb.append(".0\" width=\"100.0\" height=\"20.0\"/>\n")
        sb.append(".0\" width=\"150.0\" height=\"40.0\"/>\n")
    else:
        #sb.append(".0\" width=\"300.0\" height=\"50.0\"/>\n")
        sb.append(".0\" width=\"400.0\" height=\"75.0\"/>\n")
    sb.append("\t\t\t\t\t<y:Fill color=\"")
    if kind == "internal" or kind == "telomeric":
        sb.append(nodeColors[chrm-1])
    else:
        #sb.append(nodeColors[17]);
        sb.append("#FFD700")
    sb.append("\" transparent=\"false\"/>\n")
    sb.append("\t\t\t\t\t<y:BorderStyle ")
    if not dashed:
        sb.append("type=\"line\"")
    else:
        sb.append("type=\"dashed\"")
    #sb.append(" width=\"3.0\"")
    sb.append(" width=\"5.0\"")
    sb.append(" color=\"#000000\"/>\n")
    #sb.append("\t\t\t\t\t<y:NodeLabel fontSize=\"20\" x=\"-20\" y=\"18\">")
    if kind == "internal" or kind == "telomeric":
        sb.append("\t\t\t\t\t<y:NodeLabel fontSize=\"35\" modelName=\"internal\" modelPosition=\"c\">")
    else:
        sb.append("\t\t\t\t\t<y:NodeLabel fontSize=\"35\" modelName=\"internal\" modelPosition=\"t\">")
    sb.append(nodeLabel)
    sb.append("</y:NodeLabel>\n")
    sb.append("\t\t\t\t\t<y:Shape type=\"")
    sb.append(shape)
    sb.append("\"/>\n")
    sb.append("\t\t\t\t</y:ShapeNode>\n")
    sb.append("\t\t\t</data>\n")
    sb.append("\t\t</node>\n")
    return "".join(sb)

def drawEdgeGraphML(n1,n2,count, kind,dashed=False,vert=0):
    sb = []
    if kind == "interval":
        color = "#FF0000"
    elif kind == "reference":
        color = "#000000"
    else:
        color = "#0000FF"

    ID = kind + ":" + n1 + ":" + n2 + ":" + str(vert)

    sb.append("\t\t<edge id=\"")
    sb.append(ID)
    sb.append("\" source=\"nodeLoc")
    sb.append(n1+ ":" + str(vert))
    sb.append("\" target=\"nodeLoc")
    sb.append(n2+ ":" + str(vert))
    sb.append("\">\n")
    sb.append("\t\t\t<data key=\"d2\">\n")

    if kind == "variant" and n1!=n2:
        sb.append("\t\t\t\t<y:ArcEdge>\n")
    else:
        sb.append("\t\t\t\t<y:PolyLineEdge>\n")
    sb.append("\t\t\t\t\t<y:LineStyle ")
    if not dashed:
        sb.append("type=\"line\"")
    else:
        sb.append("type=\"dashed\"")
    #sb.append(" width=\"2.0\" ")
    sb.append(" width=\"5.0\" ")
    sb.append("color=\"")
    sb.append(color)
    sb.append("\"/>\n")
    sb.append("\t\t\t\t\t<y:Path sx=\"0.0\" sy=\"0")
    sb.append(".0\" tx=\"0.0\" ty=\"0")
    sb.append(".0\"/>\n")
    sb.append("\t\t\t\t\t<y:BendStyle smoothed=\"true\"/>\n")
    sb.append("\t\t\t\t\t<y:Arrows source=\"none\" target=\"")
    sb.append("none")
    sb.append("\"/>\n")
    sb.append("\t\t\t\t\t<y:EdgeLabel ")
    sb.append("y=\"4\" ")
    sb.append("visible=\"true\"")
    #sb.append(" fontSize=\"20\"")
    sb.append(" fontSize=\"45\"")
    sb.append(" preferredPlacement=")
    sb.append("\"target\">")
    sb.append(str(count))
    sb.append("</y:EdgeLabel>\n")
    if kind=="variant" and n1!=n2:
        sb.append("\t\t\t\t</y:ArcEdge>\n")
    else:
        sb.append("\t\t\t\t</y:PolyLineEdge>\n")
    sb.append("\t\t\t</data>\n")
    sb.append("\t\t</edge>\n")
    return "".join(sb)

def getFooter():
    sb = []
    sb.append("\t</graph>\n")
    sb.append("</graphml>\n")
    return "".join(sb)

#################################################
def alternatingBFS(adj1,adj2,s):
    #print "Running BFS from ", s
    n = len(adj1)
    q = []
    visited = [False for i in range(n)]
    d = [-1 for i in range(n)]
    p = [None for i in range(n)]

    visited[s] = True
    d[s] = 0
    p[s] = None

    q.append(s)
    while len(q)>0:
        u = q.pop(0)
        if d[u]%2==0:
            adj = adj1
        else:
            adj = adj2
        for v in adj[u]:
            if not visited[v]:
                visited[v] = True
                d[v] = d[u] + 1
                p[v] = u
                q.append(v)
    possible = range(n)
    for u in range(n):
        if u==s or d[u]==-1:
            continue
        if d[u]%2==1 and u in adj2[s]:
            v = u
            #print "Cycle:"
            while v!=s:
                if v in possible:
                    possible.remove(v)
                #print v
                v=p[v]
            #print s
            if s in possible:
                possible.remove(s)
    possible = sorted(possible)
    return possible

def impossiblePairs(source,target):
    #n =source.n
    extremities = source.extAndTels
    sourceAdj = source.adj
    sourceCn = source.cn
    targetAdj = target.adj
    targetCn = target.cn
    #leftTelExt =  leftTelomereExt()
    #rightTelExt = rightTelomereExt(n)
    pairs = []
    for e1 in extremities:
        for e2 in extremities:
            if e2<e1:
                continue
            pairs.append((e1,e2))
    #pairs.remove((leftTelExt,rightTelExt))
    for t1 in source.telExt:
        for t2 in source.telExt:
            if t2<t1:
                continue
            pairs.remove((t1,t2))
    sourceAdjList = adjMatToLists(sourceAdj)
    targetAdjList = adjMatToLists(targetAdj)
    impossible = []
    for e1,e2 in pairs:
        if (sourceAdj[e1,e2] == 0 and targetAdj[e1,e2] ==0 and
                    sourceCn[getGeneFromExt(e1)] == sourceCn[getGeneFromExt(e2)] and
                    sourceCn[getGeneFromExt(e2)] == targetCn[getGeneFromExt(e1)] and
                    targetCn[getGeneFromExt(e1)] == targetCn[getGeneFromExt(e2)] and
                    #((sourceAdj[e1,:]>0)==(targetAdj[e1,:]>0)).all() and ((sourceAdj[e2,:]>0)==(targetAdj[e2,:]>0)).all() and
                    e2 in alternatingBFS(sourceAdjList,targetAdjList,e1) and e2 in alternatingBFS(targetAdjList,sourceAdjList,e1)):
            impossible.append((e1,e2))
    print "Zeroed",len(impossible),"pairs and that's", 100.0*len(impossible)/len(pairs), "%"
    return impossible

#np.sum(sourceAdj[e1,:])==np.sum(targetAdj[e1,:]) and np.sum(sourceAdj[e2,:])==np.sum(targetAdj[e2,:]) and \

def shortestSimplePath(s,t,genome):
    n = nToNumberExtWithTels(len(genome))
    adj = adjMatToLists(genome.adj)
    q = []
    visited = [False for i in range(n)]
    d = [-1 for i in range(n)]
    p = [None for i in range(n)]
    prevEdge = [None for i in range(n)]

    visited[s] = True
    d[s] = 0
    p[s] = None
    prevEdge[s] = None

    def updateDS(v,u):
        visited[v] = True
        d[v] = d[u] + 1
        p[v] = u
        q.append(v)

    q.append(s)
    while len(q)>0:
        u = q.pop(0)
        for v in adj[u]:
            if not visited[v] and prevEdge[u]!="adj":
                updateDS(v,u)
                prevEdge[v] = "adj"
        v = getOtherExt(u)
        if v>=n: continue
        if not visited[v] and len(adj[u])>0 and prevEdge[u]!="gene":
                updateDS(v,u)
                prevEdge[v] = "gene"

    path = []
    u = t
    while u!=s:
        path.append(u)
        path.append(prevEdge[u])
        u = p[u]
    path.append(u)
    path.reverse()
    #print path
    return path

def pathHasCnBound(path,genome,bound=1):
    for i in range(0,len(path)-2,2):
        u,kind,v = path[i:i+3]
        cn = 0
        if kind=="adj":
            cn = genome.adj[u,v]
        if kind=="gene":
            cn = genome.cn[getGeneFromExt(u)]
        if cn<bound:
            return False
    return True


def hasDuplicatedChomosomes(genome):
    path = shortestSimplePath(leftTelomereExt(),rightTelomereExt(genome.n),genome)
    return pathHasCnBound(path,genome,2)

def hasTheSameChromosome(g1,g2):
    path = shortestSimplePath(leftTelomereExt(),rightTelomereExt(g1.n),g1)
    return pathHasCnBound(path,g2)

def decomposeChromosomes(genome):
    n = genome.n
    genesOnly = genome.genes
    extOnly = genome.extremities
    extAndTels = genome.extAndTels
    adj = genome.adj
    cn = genome.cn
    leftTelExt =  leftTelomereExt()
    rightTelExt = rightTelomereExt(n)
    telExt = genome.telExt
    chroms = cn[0]

    print "before:",genome

    if chroms<2:
        return n,genome

    model = Model("decomposeChromosomes")

    fAdj = {}
    for k in range(chroms):
        for e1 in extAndTels:
            for e2 in extAndTels:
                fAdj[k,e1,e2] = model.addVar(obj=0, vtype=GRB.INTEGER, lb=0, ub=adj[e1,e2], name='fAdj['+str(k)+','+str(e1)+','+str(e2)+']')
    fGene = {}
    for k in range(chroms):
        for g in genesOnly:
            fGene[k,getGeneTail(g),getGeneHead(g)] = model.addVar(obj=0, vtype=GRB.INTEGER, lb=0, ub=cn[g], name='fGene['+str(k)+','+str(getGeneTail(g))+','+str(getGeneHead(g))+']')
            fGene[k,getGeneHead(g),getGeneTail(g)] = model.addVar(obj=0, vtype=GRB.INTEGER, lb=0, ub=cn[g], name='fGene['+str(k)+','+str(getGeneHead(g))+','+str(getGeneTail(g))+']')

    model.modelSense = GRB.MINIMIZE

    model.setParam('OutputFlag', False )

    model.update()

    #constraints
    for k in range(chroms):
        model.addConstr(quicksum(fAdj[k,leftTelExt,e] for e in extAndTels)==1)
        model.addConstr(quicksum(fAdj[k,e,rightTelExt] for e in extAndTels)==1)
        for e1 in extOnly:
            model.addConstr(quicksum(fAdj[k,e1,e2] for e2 in extAndTels)==fGene[k,getOtherExt(e1),e1])
            model.addConstr(quicksum(fAdj[k,e2,e1] for e2 in extAndTels)==fGene[k,e1,getOtherExt(e1)])

    for e1 in extAndTels:
        for e2 in extAndTels:
            if e1 != e2:
                model.addConstr(quicksum(fAdj[k,e1,e2]+fAdj[k,e2,e1] for k in range(chroms))==adj[e1,e2])
            else:
                model.addConstr(quicksum(fAdj[k,e1,e1] for k in range(chroms))==adj[e1,e1])
    for g in genesOnly:
        model.addConstr(quicksum(fGene[k,getGeneTail(g),getGeneHead(g)]+fGene[k,getGeneHead(g),getGeneTail(g)] for k in range(chroms))==cn[g])

    if chroms>1:
        model.setObjective(quicksum((fGene[0,getGeneTail(g),getGeneHead(g)]-fGene[1,getGeneTail(g),getGeneHead(g)])*(fGene[0,getGeneTail(g),getGeneHead(g)]-fGene[1,getGeneTail(g),getGeneHead(g)]) for g in genesOnly)+
                           quicksum((fAdj[0,e1,e2]-fAdj[1,e1,e2])*(fAdj[0,e1,e2]-fAdj[1,e1,e2]) for e1 in extAndTels for e2 in extAndTels),GRB.MINIMIZE)

    model.optimize()

    if model.status == GRB.status.INFEASIBLE:
        print "Infeasible!!!"
        #model.computeIIS()
        #model.write(r"D:\ronzeira\model.ilp")
        return -1

    chrom_list = []
    for k in range(chroms):
        #print "---------------------------------"
        #print "Chromosome copy ",k
        adjNew = adj.copy()
        cnNew = cn.copy()
        for e1 in extAndTels:
            for e2 in extAndTels:
                adjNew[e1,e2] = fAdj[k,e1,e2].x+fAdj[k,e2,e1].x
        for g in genesOnly:
            cnNew[g] = fGene[k,getGeneTail(g),getGeneHead(g)].x+fGene[k,getGeneHead(g),getGeneTail(g)].x
        cnNew[0] = 1
        cnNew[-1] = 1
        cc = chromosomeContainer()
        cc.setMatrices(cnNew,adjNew)
        chrom_list.append(cc)
        #print cc
    #print "---------------------------------"
    #print "Genome without first chromosome "
    new = genome-chrom_list[0]
    #print model.objVal
    print "after:",new
    if model.objVal==0:
        return model.objVal,new
    else:
        return model.objVal,genome

def removeLongestDuplicatedChromosomes(genome):
    genesOnly = genome.genes
    extOnly = genome.extremities
    extAndTels = genome.extAndTels
    adj = genome.adj
    cn = genome.cn
    telExt = genome.telExt
    telGenes = genome.telGenes
    chroms = sum(cn[telGenes])/2

    if chroms<=2:
        return -1

    print "before:",genome

    model = Model("removeLongestDuplicatedChromosomes")

    s = -1
    fAdj = {}
    for e1 in extAndTels:
        for e2 in extAndTels:
            fAdj[e1,e2] = model.addVar(obj=0, vtype=GRB.INTEGER, lb=0, ub=min(adj[e1,e2],1))
    for t in telExt:
        fAdj[s,t] = model.addVar(obj=1, vtype=GRB.INTEGER, lb=0, ub=min(cn[getGeneFromExt(t)],1))
        fAdj[t,s] = model.addVar(obj=1, vtype=GRB.INTEGER, lb=0, ub=min(cn[getGeneFromExt(t)],1))
    fGene = {}
    for g in genesOnly:
        fGene[getGeneTail(g),getGeneHead(g)] = model.addVar(obj=1, vtype=GRB.INTEGER, lb=0, ub=min(cn[g],1))
        fGene[getGeneHead(g),getGeneTail(g)] = model.addVar(obj=1, vtype=GRB.INTEGER, lb=0, ub=min(cn[g],1))

    model.modelSense = GRB.MAXIMIZE

    model.setParam('OutputFlag', False )

    model.update()

    #constraints
    model.addConstr(quicksum(fAdj[s,t] for t in telExt)==1)
    model.addConstr(quicksum(fAdj[t,s] for t in telExt)==1)
    for e1 in extOnly:
        model.addConstr(quicksum(fAdj[e1,e2] for e2 in extAndTels)==fGene[getOtherExt(e1),e1])
        model.addConstr(quicksum(fAdj[e2,e1] for e2 in extAndTels)==fGene[e1,getOtherExt(e1)])
    for t in telExt:
        model.addConstr(quicksum(fAdj[t,e] for e in extOnly)==fAdj[s,t])
        model.addConstr(quicksum(fAdj[e,t] for e in extOnly)==fAdj[t,s])

    for e1 in extAndTels:
        for e2 in extAndTels:
            if e1 != e2:
                model.addConstr(2*(fAdj[e1,e2]+fAdj[e2,e1])<=adj[e1,e2])
            else:
                model.addConstr(2*fAdj[e1,e1]<=adj[e1,e1])

    for g in genesOnly:
        model.addConstr(2*(fGene[getGeneTail(g),getGeneHead(g)]+fGene[getGeneHead(g),getGeneTail(g)])<=cn[g])

    for t in telExt:
        model.addConstr(2*(fAdj[s,t]+fAdj[t,s])<=cn[getGeneFromExt(t)])

    model.optimize()

    if model.status == GRB.status.INFEASIBLE:
        print "Infeasible!!!"
        #model.computeIIS()
        #model.write(r"D:\ronzeira\model.ilp")
        return -1

    print "objective " , model.objVal

    print "---------------------------------"
    print "Chromosome "
    adjNew = adj.copy()
    cnNew = cn.copy()
    for e1 in extAndTels:
        for e2 in extAndTels:
            if e1!=e2:
                adjNew[e1,e2] = (fAdj[e1,e2].x+fAdj[e2,e1].x)
            else:
                adjNew[e1,e2] = fAdj[e1,e2].x
    for g in genesOnly:
        cnNew[g] = (fGene[getGeneTail(g),getGeneHead(g)].x+fGene[getGeneHead(g),getGeneTail(g)].x)
    for t in telExt:
        if fAdj[s,t].x>0:
            t1 = t
        if fAdj[t,s].x>0:
            t2 = t
    if t1!=t2:
        cnNew[getGeneFromExt(t1)] = 1
        cnNew[getGeneFromExt(t2)] = 1
    else:
        cnNew[getGeneFromExt(t1)] = 2
    for t in telExt:
        if t!=t1 and t!=t2:
            cnNew[getGeneFromExt(t)] = 0
    cc = chromosomeContainer()
    cc.setMatrices(cnNew,adjNew,telGenes)
    print cc
    #print "---------------------------------"
    #print "Genome without one chromosome "
    #new = genome-cc
    #print "after:",new
    #return new
    return cc


#########################################################################

def twoTupleToThreeTuple(k,t):
    return k,t[0],t[1]

def orderPair(e1,e2):
    if e1<=e2:
        return e1,e2
    else:
        return e2,e1

def createPairList(l):
    pairs = []
    for i in l:
        for j in l:
            if j<i:
                continue
            pairs.append((i,j))
    return pairs

def kAndOrderedPair(k,e1,e2):
    t1,t2 = orderPair(e1,e2)
    return k,t1,t2
####################################################
## The main function to search for a scenario from source to taeger with distance d

def multiStepSortingIlp(source,target,d=1,check=False,impossible=[],graphFileName=""):
    #get data from inputs
    #n = source.n
    genesOnly = source.genes
    extOnly = source.extremities
    extAndTels = source.extAndTels
    numNodes = len(extAndTels)
    genesAndTels = source.geneAndTels
    sourceAdj = source.adj
    sourceCn = source.cn
    targetAdj = target.adj
    targetCn = target.cn
    allNodePairs = createPairList(extAndTels)
    extNodePairs = createPairList(extOnly)
    #leftTelExt =  leftTelomereExt()
    #rightTelExt = rightTelomereExt(n)
    telExt = source.telExt
    telNodePairs = createPairList(telExt)

    model = Model("MultiStep")

    #Vairables:
    #Amplifications
    dup_adj = {}
    for k in range(d):
        for e1,e2 in allNodePairs:
            dup_adj[k,e1,e2] = model.addVar(obj=0, vtype=GRB.BINARY, name='dup_adj['+str(k)+','+str(e1)+','+str(e2)+']')

    dup_gene = {}
    for k in range(d):
        for g in genesAndTels:
            dup_gene[k,g] = model.addVar(obj=0, vtype=GRB.BINARY, name='dup_gene['+str(k)+','+str(g)+']')

    duplication = {}
    for k in range(d):
        for e1,e2 in extNodePairs:
            duplication[k,e1,e2] = model.addVar(obj=1, vtype=GRB.BINARY, name='duplication['+str(k)+','+str(e1)+','+str(e2)+']')
        for t1,t2 in telNodePairs:
            duplication[k,t1,t2] = model.addVar(obj=1, vtype=GRB.BINARY, name='duplication['+str(k)+','+str(t1)+','+str(t2)+']')

    #deletion
    del_adj = {}
    for k in range(d):
        for e1,e2 in allNodePairs:
            del_adj[k,e1,e2] = model.addVar(obj=0, vtype=GRB.BINARY, name='del_adj['+str(k)+','+str(e1)+','+str(e2)+']')
        for t1,t2 in telNodePairs:
            del_adj[k,t1,t2] = 0

    del_gene = {}
    for k in range(d):
        for g in genesAndTels:
            del_gene[k,g] = model.addVar(obj=0, vtype=GRB.BINARY, name='del_gene['+str(k)+','+str(g)+']')

    deletion = {}
    for k in range(d):
        for e1,e2 in allNodePairs:
            deletion[k,e1,e2] = model.addVar(obj=1, vtype=GRB.BINARY, name='deletion['+str(k)+','+str(e1)+','+str(e2)+']')

    #DCJ
    cut = {}
    for k in range(d):
        for e1 in extAndTels:
            for e2 in extAndTels:
                cut[k,e1,e2] = model.addVar(obj=0.25, vtype=GRB.BINARY, name='cut['+str(k)+','+str(e1)+','+str(e2)+']')

    join = {}
    for k in range(d):
        for e1 in extAndTels:
            for e2 in extAndTels:
                join[k,e1,e2] = model.addVar(obj=0.25, vtype=GRB.BINARY, name='join['+str(k)+','+str(e1)+','+str(e2)+']')

    addTel = {}
    for k in range(d):
        for t1,t2 in telNodePairs:
            addTel[k,t1,t2] = model.addVar(obj=0, vtype=GRB.BINARY, name='addTel['+str(k)+','+str(e1)+','+str(e2)+']')

    #CN and adj
    cnVec = {}
    for k in range(1,d):
        for g in genesAndTels:
            cnVec[k,g] = model.addVar(obj=0, vtype=GRB.INTEGER, lb = 0 , name='cnVec['+str(k)+','+str(g)+']')
            #cnVec[k,g] = model.addVar(obj=0, lb = 0, name='cnVec['+str(k)+','+str(g)+']')

    for g in genesAndTels:
        cnVec[0,g] = int(sourceCn[g])
        cnVec[d,g] = int(targetCn[g])

    adjMat = {}
    for k in range(1,d):
        for e1,e2 in allNodePairs:
            adjMat[k,e1,e2] = model.addVar(obj=0, vtype=GRB.INTEGER, lb =0 , name='adjMat['+str(k)+','+str(e1)+','+str(e2)+']')
            #adjMat[k,e1,e2] = model.addVar(obj=0, lb =0 , name='adjMat['+str(k)+','+str(e1)+','+str(e2)+']')

    for e1,e2 in allNodePairs:
        adjMat[0,e1,e2] = int(sourceAdj[e1,e2])
        adjMat[d,e1,e2] = int(targetAdj[e1,e2])

    #test flow
    f = {}
    r = {}
    q = {}
    b = {}
    for k in range(d):
        for e1 in extAndTels:
            for e2 in extAndTels:
                f[k,e1,e2] = model.addVar(obj=0, vtype=GRB.INTEGER, lb=-numNodes, ub=numNodes, name='f['+str(k)+','+str(e1)+','+str(e2)+']')
                #f[k,e1,e2] = model.addVar(obj=0, vtype=GRB.CONTINUOUS, lb=-2*(n+1), ub=2*(n+1), name='f['+str(k)+','+str(e1)+','+str(e2)+']')
            r[k,e1] = model.addVar(obj=0, vtype=GRB.BINARY, name='r['+str(k)+','+str(e1)+']')
            q[k,e1] = model.addVar(obj=0, vtype=GRB.INTEGER, lb=0, name='q['+str(k)+','+str(e1)+']')
        b[k] = model.addVar(obj=0, vtype=GRB.INTEGER, lb=0, name='b['+str(k)+']')
    fd = {}
    rd = {}
    qd = {}
    bd = {}
    for k in range(d):
        for e1 in extAndTels:
            for e2 in extAndTels:
                fd[k,e1,e2] = model.addVar(obj=0, vtype=GRB.INTEGER, lb=-numNodes, ub=numNodes, name='fd['+str(k)+','+str(e1)+','+str(e2)+']')
            rd[k,e1] = model.addVar(obj=0, vtype=GRB.BINARY, name='rd['+str(k)+','+str(e1)+']')
            qd[k,e1] = model.addVar(obj=0, vtype=GRB.INTEGER, lb=0, name='qd['+str(k)+','+str(e1)+']')
        bd[k] = model.addVar(obj=0, vtype=GRB.INTEGER, lb=0, name='bd['+str(k)+']')

    if check:
        for k in range(1,d-1):
            for t1,t2 in telNodePairs:
                adjMat[k,t1,t2] = 0
                #join[k,t1,t2] = 0
                #join[k,t2,t1] = 0
                #cut[k,t1,t2] = 0
                #cut[k,t2,t1] = 0
                #del_adj[k,t1,t2] = 0
                #dup_adj[k,t1,t2] = 0
        for e1,e2 in allNodePairs:
            if (e1,e2) in impossible:
                for k in range(1,d-1):
                    adjMat[k,e1,e2] = 0
                    join[k,e1,e2] = 0
                    join[k,e2,e1] = 0
                    cut[k,e1,e2] = 0
                    cut[k,e2,e1] = 0
                    del_adj[k,e1,e2] = 0
                    deletion[k,e1,e2] = 0
                    dup_adj[k,e1,e2] = 0
                    duplication[k,e1,e2] = 0
                    #f[k,e1,e2] = 0
                    #fd[k,e1,e2] = 0
                for k in [0]:
                    join[k,e1,e2] = 0
                    join[k,e2,e1] = 0
                    cut[k,e1,e2] = 0
                    cut[k,e2,e1] = 0
                    del_adj[k,e1,e2] = 0
                    deletion[k,e1,e2] = 0
                    dup_adj[k,e1,e2] = 0
                    duplication[k,e1,e2] = 0
                    #f[k,e1,e2] = 0
                    #fd[k,e1,e2] = 0


        #print "Zeroed",count,"pairs and that's", 100*float(count)/len(pairs), "%"

    model.modelSense = GRB.MINIMIZE

    model.setParam('OutputFlag', False )

    model.setParam(GRB.Param.Threads, multiprocessing.cpu_count() - 2)
    model.setParam(GRB.Param.MIPGapAbs, 0.9)

    model.setParam(GRB.Param.BestObjStop, d)
    #model.setParam(GRB.Param.Cutoff, d)

    model.update()

    #constraints
    for k in range(d):
        #degree preservation
        for e in extAndTels:
            model.addConstr(cnVec[k,getGeneFromExt(e)]==quicksum(adjMat[kAndOrderedPair(k,e,e2)] for e2 in extAndTels)+adjMat[k,e,e])

        # #at most one duplication
        # model.addConstr(quicksum(duplication[twoTupleToThreeTuple(k,p)] for p in extNodePairs)+quicksum(duplication[twoTupleToThreeTuple(k,p)] for p in telNodePairs)<=1)
        #
        # #at most one deletion
        # model.addConstr(quicksum(deletion[twoTupleToThreeTuple(k,p)] for p in allNodePairs)<=1)
        #
        # #at most 2 double cut and join
        # model.addConstr(quicksum(cut[k,e1,e2] for e1 in extAndTels for e2 in extAndTels)<=2)
        #
        # #operations should be disjoint
        # for p in extNodePairs:
        #     model.addConstr(cut[twoTupleToThreeTuple(k,p)]+join[twoTupleToThreeTuple(k,p)]+duplication[twoTupleToThreeTuple(k,p)]+deletion[twoTupleToThreeTuple(k,p)]+del_adj[twoTupleToThreeTuple(k,p)]+dup_adj[twoTupleToThreeTuple(k,p)],'<=',1)
        # for p in telNodePairs:
        #     model.addConstr(addTel[twoTupleToThreeTuple(k,p)]+duplication[twoTupleToThreeTuple(k,p)]+deletion[twoTupleToThreeTuple(k,p)]<=1)
        # for e in extOnly:
        #     for t in telExt:
        #         model.addConstr(cut[kAndOrderedPair(k,e,t)]+join[kAndOrderedPair(k,e,t)]+deletion[kAndOrderedPair(k,e,t)]+del_adj[kAndOrderedPair(k,e,t)]+dup_adj[kAndOrderedPair(k,e,t)],'<=',1)


        #Only one operation in each iteration
        model.addConstr(quicksum(duplication[twoTupleToThreeTuple(k,p)] for p in telNodePairs)+\
            quicksum(duplication[twoTupleToThreeTuple(k,p)] for p in extNodePairs)+\
            quicksum(deletion[twoTupleToThreeTuple(k,p)] for p in allNodePairs)+ \
            0.5*quicksum(cut[k,e1,e2] for e1 in extAndTels for e2 in extAndTels)<=1)

        #add at most one telomere pair
        model.addConstr(quicksum(addTel[twoTupleToThreeTuple(k,p)] for p in telNodePairs)<=1)

        #DCJ
        #can't cut if there is no adjacency
        for e1 in extAndTels:
            for e2 in extAndTels:
                if e1 in telExt and e2 in telExt:
                    continue
                if e1!=e2:
                    model.addConstr(cut[k,e1,e2]+cut[k,e2,e1],'<=',adjMat[kAndOrderedPair(k,e1,e2)])
                else:
                    model.addConstr(cut[k,e1,e2],'<=',adjMat[kAndOrderedPair(k,e1,e2)])
        #the number of cuts=the number of joins for each adj
        for e1 in extAndTels:
            model.addConstr(quicksum(cut[k,e1,e2] for e2 in extAndTels)+quicksum(cut[k,e2,e1] for e2 in extAndTels)==quicksum(join[k,e1,e2] for e2 in extAndTels)+quicksum(join[k,e2,e1] for e2 in extAndTels))

        #add tels
        for p in telNodePairs:
            model.addConstr(cut[twoTupleToThreeTuple(k,p)] == addTel[twoTupleToThreeTuple(k,p)])
            if p[0]!=p[1]:
                cut[k,p[1],p[0]] = 0

        #DUPLICATION
        #can't cut if there is no gene
        for g in genesAndTels:
            model.addConstr(dup_gene[k,g],'<=',cnVec[k,g])
        #can't amplify if there is no adjacency
        for e1,e2 in allNodePairs:
            model.addConstr(dup_adj[k,e1,e2],'<=',adjMat[k,e1,e2])

        #for every node one in and one out
        for e in extOnly:
            model.addConstr(dup_gene[k,getGeneFromExt(e)] == quicksum(duplication[kAndOrderedPair(k,e,e2)] for e2 in extOnly) + quicksum(dup_adj[kAndOrderedPair(k,e,e2)] for e2 in extAndTels))
        for t in telExt:
            model.addConstr(dup_gene[k,getGeneFromExt(t)] == quicksum(duplication[kAndOrderedPair(k,t,t2)] for t2 in telExt))
            model.addConstr(quicksum(duplication[kAndOrderedPair(k,t,t2)] for t2 in telExt) == quicksum(dup_adj[kAndOrderedPair(k,t,e2)] for e2 in extAndTels))

        #test flow
        model.addConstr(b[k] == 2*quicksum(dup_gene[k,g] for g in genesOnly) + quicksum(dup_gene[k,getGeneFromExt(t)] for t in telExt))
        for e1 in extAndTels:
            for e2 in extAndTels:
                model.addConstr(f[k,e1,e2],'==',-f[k,e2,e1])
                if e1==e2:
                    f[k,e1,e2]=0
                    #model.addConstr(f[k,e1,e2]==0)
                elif getGeneFromExt(e1)==getGeneFromExt(e2):
                    model.addConstr(f[k,e1,e2],'<=',numNodes*dup_gene[k,getGeneFromExt(e1)])
                else:
                    model.addConstr(f[k,e1,e2],'<=',numNodes*dup_adj[kAndOrderedPair(k,e1,e2)])

            model.addConstr(quicksum(f[k,e1,e2] for e2 in extAndTels)==dup_gene[k,getGeneFromExt(e1)]-q[k,e1])
            if e1 in telExt:
                model.addConstr(r[k,e1]==quicksum(duplication[kAndOrderedPair(k,e1,t2)] for t2 in telExt if t2>e1))
            else:
                model.addConstr(r[k,e1]==quicksum(duplication[kAndOrderedPair(k,e1,e2)] for e2 in extOnly if e2>e1))
            model.addConstr(q[k,e1]<=b[k])
            model.addConstr(q[k,e1],'<=',numNodes*r[k,e1])
            model.addConstr(q[k,e1]>=b[k]-numNodes*(1-r[k,e1]))

        #DELETION
        #can't delete if there is no adjacency
        for e1,e2 in allNodePairs:
            model.addConstr(del_adj[k,e1,e2],'<=',adjMat[k,e1,e2])

        #can't cut if there is no gene
        for g in genesAndTels:
            model.addConstr(del_gene[k,g],'<=',cnVec[k,g])

        #for each ext, an adjacency is deleted iff the gene was deleted or it is the end of the deletion
        for e1 in extOnly:
            model.addConstr(quicksum(del_adj[kAndOrderedPair(k,e1,e2)] for e2 in extAndTels)==\
                            del_gene[k,getGeneFromExt(e1)]+\
                            quicksum(deletion[kAndOrderedPair(k,e1,e2)] for e2 in extAndTels))
        for t in telExt:
            model.addConstr(del_gene[k,getGeneFromExt(t)] == quicksum(deletion[kAndOrderedPair(k,t,t2)] for t2 in telExt))
            model.addConstr(quicksum(deletion[kAndOrderedPair(k,t,e2)] for e2 in extAndTels) == quicksum(del_adj[kAndOrderedPair(k,t,e2)] for e2 in extAndTels))

        #test flow
        model.addConstr(bd[k] == 2*quicksum(del_gene[k,g] for g in genesOnly) + quicksum(del_gene[k,getGeneFromExt(t)] for t in telExt))
        for e1 in extAndTels:
            for e2 in extAndTels:
                model.addConstr(fd[k,e1,e2],'==',-fd[k,e2,e1])
                if e1==e2:
                    fd[k,e1,e2]=0
                elif getGeneFromExt(e1)==getGeneFromExt(e2):
                    model.addConstr(fd[k,e1,e2],'<=',numNodes*del_gene[k,getGeneFromExt(e1)])
                else:
                    model.addConstr(fd[k,e1,e2],'<=',numNodes*del_adj[kAndOrderedPair(k,e1,e2)])
            model.addConstr(quicksum(fd[k,e1,e2] for e2 in extAndTels)==del_gene[k,getGeneFromExt(e1)]-qd[k,e1])
            model.addConstr(rd[k,e1]==quicksum(deletion[kAndOrderedPair(k,e1,e2)] for e2 in extAndTels if e2>e1))
            model.addConstr(qd[k,e1]<=bd[k])
            model.addConstr(qd[k,e1],'<=',numNodes*rd[k,e1])
            model.addConstr(qd[k,e1]>=bd[k]-numNodes*(1-rd[k,e1]))

        #MODIFY
        #modify genome
        for g in genesAndTels:
            if g in genesOnly:
                model.addConstr(cnVec[k+1,g] == cnVec[k,g] + dup_gene[k,g] - del_gene[k,g])
            else:
                if getGeneLeftExt(g) in telExt:
                    t = getGeneLeftExt(g)
                if getGeneRightExt(g) in telExt:
                    t = getGeneRightExt(g)
                model.addConstr(cnVec[k+1,g] == cnVec[k,g] + dup_gene[k,g] - del_gene[k,g] + quicksum(addTel[kAndOrderedPair(k,t,t2)] for t2 in telExt if t2!=t) + 2*addTel[k,t,t])

        #modify extremities
        for e1,e2 in extNodePairs:
            if e1!=e2:
                model.addConstr(adjMat[k+1,e1,e2] ,'==', adjMat[k,e1,e2] + dup_adj[k,e1,e2] + duplication[k,e1,e2] - del_adj[k,e1,e2] + deletion[k,e1,e2] - cut[k,e1,e2] - cut[k,e2,e1] + join[k,e1,e2] + join[k,e2,e1])
            else:
                model.addConstr(adjMat[k+1,e1,e2] ,'==', adjMat[k,e1,e2] + dup_adj[k,e1,e2] + duplication[k,e1,e2] - del_adj[k,e1,e2] + deletion[k,e1,e2] - cut[k,e1,e2] + join[k,e1,e2])

        for t in telExt:
            for e in extOnly:
                model.addConstr(adjMat[kAndOrderedPair(k+1,t,e)] ,'==', adjMat[kAndOrderedPair(k,t,e)] + deletion[kAndOrderedPair(k,t,e)] + dup_adj[kAndOrderedPair(k,t,e)] - del_adj[kAndOrderedPair(k,t,e)] - cut[kAndOrderedPair(k,t,e)] + join[kAndOrderedPair(k,t,e)])

    #check if a gene needs to be deleted or amplified at least once
    #adding CNTP constraints
    # for g in genesOnly:
    #     model.addConstr(quicksum(dup_gene[k,g] for k in range(d))-quicksum(del_gene[k,g] for k in range(d))==int(targetCn[g]-sourceCn[g]))
    #     if targetCn[g]>sourceCn[g]:
    #        model.addConstr(quicksum(dup_gene[k,g] for k in range(d))>=targetCn[g]-sourceCn[g])
    #     if sourceCn[g]>targetCn[g]:
    #        model.addConstr(quicksum(del_gene[k,g] for k in range(d))>=sourceCn[g]-targetCn[g])
        # for k in range(d):
        #    model.addConstr(cnVec[k,g]>=quicksum(adjMat[kAndOrderedPair(k,getGeneTail(g),e)] for e in extremities)+adjMat[k,getGeneTail(g),getGeneTail(g)])
        #    model.addConstr(cnVec[k,g]>=quicksum(adjMat[kAndOrderedPair(k,getGeneHead(g),e)] for e in extremities)+adjMat[k,getGeneHead(g),getGeneHead(g)])

    #solve
    model.optimize()

    if model.status == GRB.status.INFEASIBLE:
        print "Infeasible!!!"
        #model.computeIIS()
        #model.write(r"D:\ronzeira\model.ilp")
        return -1

    if model.status == GRB.status.CUTOFF:
        print "Cutoff!!!"
        #model.computeIIS()
        #model.write(r"D:\ronzeira\model.ilp")
        return -1

    if model.status == GRB.status.USER_OBJ_LIMIT:
        print "Objective limit!!!"

    offset = (len(source.telGenes)/2 + 1)

    #construct solution
    print "Total distance:", model.objVal
    if model.objVal==0:
        return model.objVal
    for k in range(d):
        print "---------------------------------"
        print "Step ",k
        duplicated_genes = []
        deleted_genes = []
        for g in genesAndTels:
            if dup_gene[k,g].x==1:
                duplicated_genes.append(g)
            if del_gene[k,g].x==1:
                deleted_genes.append(g)

        duplicated_adjacencies = []
        deleted_adjacencies = []
        cut_adjacencies = []
        join_adjacencies = []
        duplications = []
        deletions = []
        for p in allNodePairs:
            if hasattr(dup_adj[twoTupleToThreeTuple(k,p)],'x') and dup_adj[twoTupleToThreeTuple(k,p)].x==1:
                duplicated_adjacencies.append(p)
            if hasattr(del_adj[twoTupleToThreeTuple(k,p)],'x') and del_adj[twoTupleToThreeTuple(k,p)].x==1:
                deleted_adjacencies.append(p)
            if hasattr(cut[twoTupleToThreeTuple(k,p)],'x') and cut[twoTupleToThreeTuple(k,p)].x==1:
                cut_adjacencies.append(p)
            if hasattr(join[twoTupleToThreeTuple(k,p)],'x') and join[twoTupleToThreeTuple(k,p)].x==1:
                join_adjacencies.append(p)
            if hasattr(cut[k,p[1],p[0]],'x') and cut[k,p[1],p[0]].x==1 and p[0]!=p[1]:
                cut_adjacencies.append((p[1],p[0]))
            if hasattr(join[k,p[1],p[0]],'x') and join[k,p[1],p[0]].x==1 and p[0]!=p[1]:
                join_adjacencies.append((p[1],p[0]))
            if twoTupleToThreeTuple(k,p) in duplication and hasattr(duplication[twoTupleToThreeTuple(k,p)],'x') and duplication[twoTupleToThreeTuple(k,p)].x==1:
                duplications.append(p)
            if hasattr(deletion[twoTupleToThreeTuple(k,p)],'x') and deletion[twoTupleToThreeTuple(k,p)].x==1:
                deletions.append(p)

        if len(graphFileName)>0:
            event = ""
            if len(cut_adjacencies)>0:
                event += "DCJ "
            if len(duplicated_genes)>0:
                event += "Duplication "
            if len(deleted_genes)>0:
                event += "Deletion "
            kGenome = source.ilpResultsToCnAdj(cnVec,adjMat,k)
            kGenome.printGraphML(graphFileName,duplicated_genes+deleted_genes,duplicated_adjacencies+deleted_adjacencies+cut_adjacencies,k==0,False,offset*k,event)

        if duplications:
            print "Duplication anchor adjacencies:", adjListToTupleList(duplications)
            # print "b=",b[k].x
            # for e in extAndTels:
            #     if hasattr(r[k,e],'x') and r[k,e].x==1:
            #         print "r[",e,"]"
            #     if q[k,e].x>0:
            #         print "q[",e,"]=",q[k,e].x
            # for e1 in extAndTels:
            #     for e2 in extAndTels:
            #         if hasattr(f[k,e1,e2],'x') and f[k,e1,e2].x>0:
            #             print "f[",e1,",",e2,"]=",f[k,e1,e2].x
        if duplicated_genes:
            print "Duplicated genes:", duplicated_genes
        if duplicated_adjacencies:
            print "Duplicated adjacencies:", adjListToTupleList(duplicated_adjacencies)
        if deletions:
            print "Deletions anchor adjacencies:", adjListToTupleList(deletions)
        if deleted_genes:
            print "Deleted genes", deleted_genes
        if deleted_adjacencies:
            print "Deleted adjacencies:", adjListToTupleList(deleted_adjacencies)
        if cut_adjacencies:
            print "Cut adjacencies:", adjListToTupleList(cut_adjacencies)
        if join_adjacencies:
            print "Join adjacencies:", adjListToTupleList(join_adjacencies)

    print "---------------------------------"
    print "Total distance:", model.objVal

    if len(graphFileName)>0:
        target.printGraphML(graphFileName,[],[],False,True,offset*d)

    return model.objVal

###########################################################

#Simulation function
import random

def selectRandomSegment(kar,lim=None):
    n = len(kar)
    start = random.randint(0,n-1)
    end = random.randint(0,n-1)
    if end<start:
        start,end = end,start
    return start,end

def selectRandomSegmentWithLimit(kar,lim=None,ignoreDuplicatedGenes = False):
    n = len(kar)
    if lim is None:
        lim = n-1
    start = random.randint(0, n-1)
    length = random.randint(0, lim-1)
    while (start+length>n-1) or (hasDuplicatedGenes(kar[start:start+length+1]) and ignoreDuplicatedGenes):
        start = random.randint(0, n-1)
        length = random.randint(0, lim-1)
    return start,start+length

def reverseRandomSegment(kar,lim=None):
    start,end = selectRandomSegmentWithLimit(kar,lim,False)
    print "Reverse ",start, " to ", end
    return reverseSegment(kar,start,end)

def reverseSegment(kar,start,end):
    seg = []
    for i in range(start,end+1):
        seg.append(-kar[i])
    seg.reverse()
    print "Reversing ",kar[start:end+1]
    new  = kar[:start] + seg + kar[end+1:]
    return new

def duplicateRandomSegment(kar,lim=None):
    start,end = selectRandomSegmentWithLimit(kar,lim,True)
    print "Duplicate ",start, " to ", end
    return duplicateSegment(kar,start,end)

def duplicateSegment(kar,start,end):
    print "Duplicating ",kar[start:end+1]
    new = kar[:end+1] + kar[start:end+1] + kar[end+1:]
    return new

def deleteRandomSegment(kar,lim=None):
    start,end = selectRandomSegmentWithLimit(kar,lim,True)
    print "Delete ",start, " to ", end
    return deleteSegment(kar,start,end)

def deleteSegment(kar,start,end):
    print "Deleting ",kar[start:end+1]
    new = kar[:start] + kar[end+1:]
    return new

def traslocateRandom(c1,c2):
    s1 = random.randint(0, len(c1)-1)
    s2 = random.randint(0, len(c2)-1)
    print "Traslocation ", s1, "in first, ", s2, " in second"
    print "->",c1[:s1],c2[s2:], "and ",c2[:s2],c1[s1:]
    return c1[:s1]+c2[s2:],c2[:s2]+c1[s1:]

def hasDuplicatedGenes(l):
    d = {}
    for a in np.abs(l):
        if a in d:
            return True
        d[a] = 1
    return False

def deletedGenes(kar):
    n = getNfromChromList(kar)
    genes = range(n)
    karGenes = np.sort(np.unique(np.abs(kar)))
    deleted = sorted(list(set(genes) - set(karGenes)))
    return deleted

def revList(l):
    new = []
    for i in l:
        new.append(-i)
    new.reverse()
    return new

def findCompressedSegments(cn,adjMat,exclude=1):
    n = len(cn)-exclude
    #n = len(cn)
    missingGenes = (cn==0).nonzero()[0]
    segments = []
    i = exclude
    #i=0
    while i<n:
        seg = [i]
        j = i+1
        if i in missingGenes:
            while j in missingGenes and j<n:
                seg.append(j)
                j += 1
        else:
            while j<n and \
            adjMat[getGeneHead(j-1),getGeneTail(j)]>0 and \
            adjMat[getGeneHead(j-1),:].sum()==adjMat[getGeneHead(j-1),getGeneTail(j)] and \
            adjMat[:,getGeneTail(j)].sum()==adjMat[getGeneHead(j-1),getGeneTail(j)] and \
                    cn[j]==cn[j-1]:
                seg.append(j)
                j += 1
        i = j
        segments.append(seg)
    return segments

def segmentsToDict(segments):
    d = {}
    for i in range(1,len(segments)+1):
        seg = segments[i-1]
        d[seg[0]] = (i,seg)
        print "Segment ",i, " is ",seg[0]," to ",seg[-1]
        d[-seg[-1]] = (-i,revList(seg))
    return d

def compressKarWithDict(kar,d):
    i = 0
    newKar = []
    while i<len(kar):
        k,seg = d[kar[i]]
        for j in range(i,i+len(seg)):
            if seg[j-i]!=kar[j]:
                #print "Error: i=",i
                return None
        i = j+1
        newKar.append(k)
    return newKar

def compressChromList(chromList):
    n = getNfromChromList(chromList)
    adjMat = chromListToAdjMat(chromList,n)
    cn = chromListToCn(chromList,n)
    segments = findCompressedSegments(cn,adjMat)
    #print segments
    d = segmentsToDict(segments)
    nTag = len(segments)
    newChromList = []
    for kar in chromList:
        newChromList.append(compressKarWithDict(kar,d))
    #print newChromList
    return newChromList,[range(1,nTag+1),range(1,nTag+1)]

def translateSegments(segments,telGenes):
    newSeg = []
    newTel = []
    for i in range(len(segments)):
        seg = segments[i]
        if seg[0] in telGenes:
            newSeg.append(seg[0:1])
            seg = seg[1:]
            newTel.append(len(newSeg)-1)
        if len(seg) >0:
            if seg[-1] in telGenes:
                if len(seg[:-1])>0:
                    newSeg.append(seg[:-1])
                newSeg.append(seg[-1:])
                newTel.append(len(newSeg)-1)
            else:
                newSeg.append(seg)
    genesAndTels = range(len(newSeg))
    genesOnly = list(set(genesAndTels).difference(set(newTel)))
    #print "newSeg",newSeg
    #print "newTel",newTel
    #print "genesOnly",genesOnly
    return newSeg, genesOnly, newTel

def compressKaryotypeMat(cn,adjMat,telGenes=None):
    if not telGenes:
        telGenes = [0,len(cn)-1]
    segments = findCompressedSegments(cn,adjMat,0)
    newSeg, genesOnly, newTel = translateSegments(segments,telGenes)
    nTag = len(newSeg)
    if nTag==len(cn):
        return cn,adjMat,telGenes
    newCn = np.zeros(len(newSeg),dtype=int)
    for i in range(nTag):
        seg = newSeg[i]
        print "Segment ",i, " is ",seg[0]," to ",seg[-1]
        newCn[i] = cn[seg[0]]
    #print newCn
    newAdj = np.zeros((2*nTag,2*nTag),dtype=int)
    for i in range(nTag):
        for j in range(nTag):
            seg1 = newSeg[i]
            seg2 = newSeg[j]
            newAdj[getGeneHead(i),getGeneHead(j)] = adjMat[getGeneHead(seg1[-1]),getGeneHead(seg2[-1])]
            newAdj[getGeneHead(i),getGeneTail(j)] = adjMat[getGeneHead(seg1[-1]),getGeneTail(seg2[0])]
            newAdj[getGeneTail(i),getGeneHead(j)] = adjMat[getGeneTail(seg1[0]),getGeneHead(seg2[-1])]
            newAdj[getGeneTail(i),getGeneTail(j)] = adjMat[getGeneTail(seg1[0]),getGeneTail(seg2[0])]
    #print adjMatToAdjList(newAdj)
    return newCn,newAdj,newTel

def compressTandemDulications(cn,adjMat):
    td = 0
    genes = range(len(cn))
    for g in genes:
        if adjMat[getGeneHead(g),getGeneTail(g)]>0:
            print "Removing ", adjMat[getGeneHead(g),getGeneTail(g)], " tandem duplications from gene ",g
            td += adjMat[getGeneHead(g),getGeneTail(g)]
            cn[g] -= adjMat[getGeneHead(g),getGeneTail(g)]
            adjMat[getGeneHead(g),getGeneTail(g)] = 0
            adjMat[getGeneTail(g),getGeneHead(g)] = 0
    return cn,adjMat,td

def randomizeKaryotype(n=100,d=1):
    kar = [range(1,n+1),range(1,n+1)]
    for i in range(d):
        #print "Before",kar
        if len(kar)==0:
            break
        p = random.random()
        if p<0.2 or len(kar)==1:
            k=0
        else:
            k = random.randint(1,len(kar)-1)
        p = random.random()
        if p<0.05 and not hasDuplicatedGenes(kar[k]):
            print "Duplicating chromosome ", k
            kar.append(kar[k][:])
        elif p<0.075:
            print "Deleting chromosome ", k
            kar = kar[:k]+kar[k+1:]
        elif p < 0.15 and len(kar)>1:
            j=k
            while j==k:
                j = random.randint(0,len(kar)-1)
            kar[k],kar[j] = traslocateRandom(kar[k],kar[j])
        else:
            j = random.randint(1,3)
            if j==1:
                kar[k] = deleteRandomSegment(kar[k],int(0.1*n))
            elif j==2:
                kar[k] = duplicateRandomSegment(kar[k],int(0.3*n))
            elif j==3:
                kar[k] = reverseRandomSegment(kar[k],int(0.3*n))
        #print "After",kar
    if len(kar)==0:
        print "Reached an empty karyotype, restarting"
        return randomizeKaryotype(n,d)
    print kar
    kar,org = compressChromList(kar)
    print kar
    return kar,org

###################################################

#search function
import time

def searchDistance(source,target,bound,lb=1,check=False,impossible=[],offset=0,graphFileName=""):
    print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    start_search = time.time()
    for i in range(max(1,lb),bound+1):
        print "*****************************************"
        start_iter = time.time()
        d = multiStepSortingIlp(source,target,i,check,impossible,graphFileName)
        d += offset
        end_iter = time.time()
        print "Iteration ",i, " time :" ,end_iter-start_iter, " seconds"
        if d>=0:
            print "#####################################"
            #print "Found distance: ",d, " with bound ",i
            end_search = time.time()
            #print "Total search time: ", end_search-start_search, " seconds. Optimal distance: ",d, " , Simulated distance: ",bound , " , Source length: ", len(source), " , Target length: ", len(target)
            print "Total search time: ", end_search-start_search, " seconds. Optimal distance: ",d, " , Simulated distance: ",bound , " , Segments: ", len(source), " , Constraint: ", check
            return d
    print "Didn't find distance within bound ",bound
    end_search = time.time()
    print "Total search time: ", end_search-start_search, " seconds. Optimal distance: ","None", " , Simulated distance: ",bound, " , Segments: ", len(source), " , Constraint: ", check
    return -1

def simplifyGenome(org,kar):
    impossible = impossiblePairs(org,kar)
    if len(impossible)>0 or org==kar:
        return 0,impossible,kar
    new = kar
    d = 0
    while True:
        if org==new:
            break
        t,new = decomposeChromosomes(new)
        if t==0:
            d+=1
            continue
        else:
            break
    impossible = impossiblePairs(org,new)
    return d,impossible,new

##########################################
#test zone


#kar1 = chromosomeContainer([[1,2],[1,2]])
#kar1.addAdjList([(3,4,-2),(3,3,1),(4,4,1)])

#kar2 = chromosomeContainer([[1,2],[1,2]])
#kar3 = chromosomeContainer([[1,2],[1,2]])
#dip2 = createDiploidChrome(kar2.n)
#kar1.mergeKar(kar2)
#kar1.addAdjList([(3,4,-1),(11,12,-1),(3,11,1),(4,12,1)])
#kar1.addAdjList([(1,2,-1),(3,11,-1),(10,9,-1)])
#kar1.cn += np.array([-1,-1,0,0,-1,-1,0,0])
#kar1.addAdjList([(14,13,-1),(12,4,-1),(5,6,-1)])
#kar1.cn += np.array([0,0,-1,-1,0,0,-1,-1])
#kar1.mergeKar(kar2)
#kar1.mergeKar(kar3)
#kar1.printGraphML("test.graphml",[0,1,4,5],[(1,2),(3,11),(10,9)])
#print kar1
# kar1.addAdjList([(3,4,-1),(11,12,-1),(3,11,1),(4,12,1)])
# kar1.addAdjList([(19,20,-1),(11,12,-1),(19,11,1),(20,12,1)])
# dip1 = kar1.getDiploid()
# #dip1.mergeKar(dip2)
# print kar1
# print dip1
# multiStepSortingIlp(dip1,kar1,2)

#kar1 = chromosomeContainer([[1,-3,1,-3,4],[1,2,3,4]])
#dip1 = kar1.getDiploid()
#multiStepSortingIlp(dip1,kar1,2)

# kar1 = chromosomeContainer([[1,2,2,3],[1,2,3]])
# dip1 = kar1.getDiploid()
# multiStepSortingIlp(kar1,dip1,1)

# kar1 = chromosomeContainer([[1,3],[1,3]])
# dip1 = kar1.getDiploid()
# print kar1.has_distance_to(dip1)
# multiStepSortingIlp(dip1,kar1,2)
# multiStepSortingIlp(kar1,dip1,2)

#####################################################
#Simulations

def runSimulations(first_d=2,max_d=1,repetitions=5):
    for i in range(first_d,max_d+1):
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "Starting simulating distance d=",i
        for j in range(repetitions):
            print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
            print "Simulating instance ",j, " with distance ",i
            kar,org = randomizeKaryotype(100,i)
            print "The karyotype: ",kar, " Original: ", org
            kar = chromosomeContainer(kar,getNfromChromList(org))
            org = chromosomeContainer(org,getNfromChromList(org))
            lb = np.sum(np.logical_and(org.adj==0,kar.adj>0))/4
            print "LB (BP/2) = ", lb
            #print "Unconstrained @@@"
            #d = searchDistance(org,kar,i,lb)

            print "+++++++++++++++++++++++++++++++++"
            print "Constrained @@@"
            #impossible = impossiblePairs(org,kar)
            #d1,impossible,kar = simplifyGenome(org,kar)
            #print "Simplified duplications: ",d1, ", Impossible pairs: ",len(impossible)
            #d2 = searchDistance(org,kar,i,lb,True,impossible)
            d2 = searchDistance(org,kar,i,lb,False)
            #print "Total distance: ", d1+d2, " , Simulated distance: ",i
        #     if d2==-1:
        #         break
        #     if abs(d - d2)>0.01:
        #         print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        #         print "Found counter example!!!"
        #         break
        # if d2==-1:
        #     break
        # if abs(d - d2)>0.01:
        #     break

#runSimulations(first_d=1,max_d=8,repetitions=20)