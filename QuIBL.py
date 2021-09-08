#Quantifying Introgression via Branch Lengths
#Miriam Miyagi

from ete3 import Tree
import itertools as itt
import sys as sys
import numpy as np
from math import log
import ConfigParser
import csv
from joblib import Parallel, delayed
import multiprocessing

class tripletT:
	def __init__(self, tb1,tb2,tb3):
	#Input in the form [[outgroup1,branch lengths,...],..]
		self.taxaSet=[tb1[0],tb2[0],tb3[0]]
		self.branchMat=[tb1[1:],tb2[1:],tb3[1:]]
		self.canonOG=None
		self.models={tb1[0]:((0,0),(0,0),0),tb2[0]:((0,0),(0,0),0),tb3[0]:((0,0),(0,0),0)}#dict()
		self.BIC={tb1[0]:0,tb2[0]:0,tb3[0]:0}#dict()
		self.null={tb1[0]:(0,0),tb2[0]:(0,0),tb3[0]:(0,0)}#dict()

	def setCanonOutgroup(self,taxon):
	#Sets the `true' outgroup to `taxon'
		if taxon in self.taxaSet:
			self.canonOG=taxon
		else: print('Error: Not a valid outgroup for this triplet.')
	
	def branches(self,taxon):
	#Returns the set of branch lengths for trees in the topology with outgroup `taxon'
		IND=self.taxaSet.index(taxon)
		return self.branchMat[IND]

	def setModel(self,outgroup, cSet, pi, lmbd):
	#Updates the dictionary of models for `outgroup'
		if outgroup in self.taxaSet:
			self.models[outgroup]=(cSet,pi,lmbd)
		else: print('Error: Not a valid outgroup for this triplet.')
	def setBIC(self,outgroup, bic):
	#Updates the BIC value for `outgroup'
		if outgroup in self.taxaSet:
			self.BIC[outgroup]=bic
		else: print('Error: Not a valid outgroup for this triplet.')
	def setNull(self,outgroup,ntup):
	#Updates the one distribution model for `outgroup'
		if outgroup in self.taxaSet:
			self.null[outgroup]=ntup
		else: print('Error: Not a valid outgroup for this triplet.')


def readin_Newick(filepath):
#Reads in a file with newick trees on separate lines
	treeList=[]
	returnlist=[]
	f=open(filepath).read()
	composite=f.split(';')
	for item in composite:
		if item.strip():
			returnlist.append((item+';').strip('\n'))
	for tree in returnlist:
		treeList.append(Tree(tree))
	return treeList

def tripSetGen(treeSet,outgroup):
#Generates the set of triplets to check.
	leaves=treeSet[0].get_leaf_names()
	leaves.remove(outgroup)
	return list(itt.combinations(leaves,3))


def getTripBranches(treeList,canonOut):
#Goes through the trees in the input and calculates/sorts the internal branch lengths.
	setOfTriplets=[]
	triples=tripSetGen(treeList,canonOut)
	lenIt=sum(1 for _ in triples)
	#output=np.zeros((4,lenIt))
	output=[[None,[],[],[]] for i in range(lenIt)]
	for counter,tree in enumerate(treeList):
		tree.set_outgroup(canonOut)
		dist=0
		#tree.set_outgroup(canonOut)
		for index,triplet in enumerate(triples):
			output[index][0]=triplet
			tempTree=tree.copy('newick')
			tempTree.prune(triplet,preserve_branch_length=True)
			if len(tempTree.expand_polytomies())>1:
				continue
			common=tempTree.get_leaves()[0].get_common_ancestor(tempTree.get_leaves())
			for x in common.get_children():
				if not (x.is_leaf()):
					dist=x.get_distance(common)
				else:
					tempInd=triplet.index(x.get_leaf_names()[0])
			if dist>0:
				output[index][tempInd+1].append(dist)
	for tripp in output:
		setOfTriplets.append(tripletT([tripp[0][0]]+tripp[1],[tripp[0][1]]+tripp[2],[tripp[0][2]]+tripp[3]))
	return setOfTriplets
	
		
def conPDF(x, C, lmbd):
#Continuous version PDF
	if x>0 and x>=C*lmbd:
		y=0.5/lmbd*np.exp(-x/lmbd)*(1+np.exp(C))
	elif x>0 and C*lmbd>x:
		y=np.sinh(x/lmbd)/(lmbd*(-1+np.exp(C)))
	else:
		y=0
	return y

def qCalc(x,cArray,pArray,lmbd):
#Calculates new `q' tuple for x
	xq=[]
	for c in range(0, len(cArray)):
		xq.append(pArray[c]*conPDF(x,cArray[c],lmbd))
	xq=[l/sum(xq) for l in xq]#xq/sum(xq)
	return list(xq)

def modelLogLik(XSet, cArray, pArray, lmbd):
#Calculates the log-likelihood of a set of branches given a model
	tempSum=0
	for x in XSet:
		L=0
		for i in range(0,len(pArray)):
			L+=pArray[i]*conPDF(x,cArray[i],lmbd)
		tempSum+=log(L)
	return tempSum

def weightedLogLik(XQList, cArray,lmbd):
#Calculates the conditional likelihood on the Q values for each datapoint
	L=0
	for tup in XQList:
		temp=0
		for indi,pival in enumerate(tup[1]):
			temp+=pival*(conPDF(tup[0],cArray[indi],lmbd))
		L+=log(temp)
	return L

def cSearch(XSet, XQList, index, lmbd, cArray):
#Return the best C value from XSet
	rescaledXSet=[x/lmbd for x in XSet]
	tempLik=[]
	testCArray=cArray[:]
	for j,testVal in enumerate(rescaledXSet):
		testCArray[index]=testVal
		tempLik.append(weightedLogLik(XQList, testCArray, lmbd))
		testCArray[index]=cArray[index]
	op=rescaledXSet[np.argmax(tempLik)]
	return op

def pUpdate(XQList):
#Updates the mixture proportions conditional on the Q values
	newpArr=[]
	for i in range(0,len(XQList[1])):
		temp=[]
		for x in XQList:
			temp.append(x[1][i])
		newpArr.append(np.mean(temp))
	return newpArr

def XQUpdate(XSet, cArray, pArray, lmbd):
#Returns a list of the data and their associated Q values
	XQList=[]
	for x in XSet:
		outXQ=[x]
		outXQ.append(qCalc(x, cArray, pArray,lmbd))
		XQList.append(outXQ)
	return XQList

def calcDeriv(XQList, cArray, lmbd):
#Calculates the derivative of the likelihood w.r.t. lambda for use in the gradient ascent.
	deriv=0
	if len(cArray)==1:
		alpha=cArray[0]*lmbd
		for x in XQList:
			deriv+=(x[0]-lmbd)/pow(lmbd,2)
		return deriv	
	elif len(cArray)==2:
		alpha=0
		beta=cArray[1]*lmbd
		for x in XQList:
			val=x[0]
			if val<beta:
				p=x[1][1]#x[2]
				q=x[1][0]#x[1]
				deriv+=(val+(beta/(-1+np.exp(beta/lmbd)))-(np.exp(2*val/lmbd)*p*(2*val-beta)+(p+2*q)*beta)/((-1+np.exp(2*val/lmbd))*p+2*(-1+np.exp(beta/lmbd))*q)-lmbd)/pow(lmbd,2)
			if val>=beta:
				p=x[1][1]#x[2]
				q=x[1][0]#x[1]
				deriv+=(x[0]-beta+((p+2*q)*beta)/(p+np.exp(beta/lmbd)*p+2*q)-lmbd)/pow(lmbd,2)
		return deriv

def gradAscent(XSet, XQList, cArray, pArray, lmbd, stepScale, numSteps, threshold):
#Conducts gradient ascent to find the optimal lambda- returns this and the search scaling factor `stepScale.' Runs until either the change in likelihood is less than `threshold' or there have been `numSteps' iterations.
	if len(cArray)==1:
		return np.mean(XSet), stepScale
	elif len(cArray)==2:
		likArray=[0]
		likArray.append(weightedLogLik(XQList, cArray, lmbd))
		sc=0
		while abs(likArray[-1]-likArray[-2])>=threshold*abs(likArray[-1]) and sc<numSteps:
			newlmbd=lmbd+calcDeriv(XQList, cArray, lmbd)*stepScale
			if newlmbd<=0.0005:
				stepScale=stepScale/2
				sc+=1
			else:
				newcArray=[x*(lmbd/newlmbd) for x in cArray]
				newLik=weightedLogLik(XQList,newcArray,newlmbd)
				if newLik > likArray[-1]:
					lmbd=newlmbd
					cArray=newcArray
					likArray.append(newLik)
					sc+=1
				else:
					stepScale=stepScale/2
					sc+=1
		return lmbd, stepScale
	
	elif len(cArray)>=3:
		print('Error: C array is too long- only K=1 or K=2 currently supported.')
		return None



def exMax(tripletSet, K, threshold, numSteps, tempScale):
#Runs the expectation maximization scheme
	for triple in tripletSet:
		for outG in triple.taxaSet:
			stepScale=tempScale
			branchData=triple.branches(outG)
			if len(branchData)<2:
				continue
			cArray=list(np.zeros(K))
			pArray=[1./K]*K
			lmbd=np.mean(branchData)
			for x in range(1,K):
				cArray[x]=x
			steps=0
			while steps<numSteps:
				XQList=XQUpdate(branchData, cArray, pArray, lmbd)
				pArray=pUpdate(XQList)
				for indy in range(1,K):
					cArray[indy]=cSearch(branchData, XQList, indy, lmbd, cArray)
				storedlmbd=lmbd
				lmbd,stepScale=gradAscent(branchData, XQList, cArray, pArray, lmbd, stepScale, numSteps, threshold)
				cArray=[x*(storedlmbd/lmbd) for x in cArray]
				steps+=1
			triple.setModel(outG, list(cArray), list(pArray), lmbd)
			triple.setBIC(outG, np.log(len(branchData))*(2*K-1)-2*(modelLogLik(branchData, list(cArray), list(pArray), lmbd)))
			triple.setNull(outG,oneDistNull(branchData))
		print(str(triple.taxaSet)+' is complete.')
	return tripletSet

def PLexMax(triple, K, threshold, numSteps, tempScale):
#Runs the expectation maximization scheme
	for outG in triple.taxaSet:
		stepScale=tempScale
		branchData=triple.branches(outG)
		if len(branchData)<2:
			continue
		cArray=list(np.zeros(K))
		pArray=[1./K]*K
		lmbd=np.mean(branchData)
		for x in range(1,K):
			cArray[x]=x
		steps=0
		while steps<numSteps:
			XQList=XQUpdate(branchData, cArray, pArray, lmbd)
			pArray=pUpdate(XQList)
			for indy in range(1,K):
				cArray[indy]=cSearch(branchData, XQList, indy, lmbd, cArray)
			storedlmbd=lmbd
			lmbd,stepScale=gradAscent(branchData, XQList, cArray, pArray, lmbd, stepScale, numSteps, threshold)
			cArray=[x*(storedlmbd/lmbd) for x in cArray]
			steps+=1
		triple.setModel(outG, list(cArray), list(pArray), lmbd)
		triple.setBIC(outG, np.log(len(branchData))*(2*K-1)-2*(modelLogLik(branchData, list(cArray), list(pArray), lmbd)))
		triple.setNull(outG,oneDistNull(branchData))
	print(str(triple.taxaSet)+' is complete.')
	return triple


def oneDistNull(branchData):
	lmbd=np.mean(branchData)
	bic=np.log(len(branchData))-2*(modelLogLik(branchData, [0.0], [1.0], lmbd))
	return (lmbd,bic)
	

def outputFormatter(outputDict,inputDict):
	#num_cores=multiprocessing.cpu_count()
	K=int(inputDict['numdistributions'])
	lthresh=float(inputDict['likelihoodthresh'])
	numsteps=int(inputDict['numsteps'])
	gAScalar=float(inputDict['gradascentscalar'])
	canonOut=str(inputDict['totaloutgroup'])
	trees=getTripBranches(readin_Newick(inputDict['treefile']),canonOut)
	multi=bool(inputDict['multiproc'])
	corecap=int(inputDict['maxcores'])
	num_cores=min(multiprocessing.cpu_count(),corecap)
	print('Analysis Started')
	if multi:
		tripletSet=Parallel(n_jobs=num_cores)(delayed(PLexMax)(triple,K,lthresh,numsteps,gAScalar) for triple in trees)
	else:
		tripletSet=exMax(getTripBranches(readin_Newick(inputDict['treefile']),canonOut), int(inputDict['numdistributions']), float(inputDict['likelihoodthresh']), int(inputDict['numsteps']), float(inputDict['gradascentscalar']))
	with open(outputDict['outputpath'],'w') as csv_out:
		fieldnames=[]
		fieldnames=['triplet','outgroup','C1','C2','mixprop1', 'mixprop2','lambda2Dist', 'lambda1Dist', 'BIC2Dist', 'BIC1Dist','count']
		out_writer=csv.DictWriter(csv_out,quoting=csv.QUOTE_MINIMAL, fieldnames=fieldnames)
		out_writer.writeheader()
		for triple in tripletSet:
			for y in triple.taxaSet:
				#print triple.models.get(y)[0][0]
				out_writer.writerow({'triplet':triple.taxaSet[0]+'_'+triple.taxaSet[1]+'_'+triple.taxaSet[2], 'outgroup': y, 'C1': triple.models.get(y)[0][0], 'C2': triple.models.get(y)[0][1],'mixprop1': triple.models.get(y)[1][0],'mixprop2': triple.models.get(y)[1][1], 'lambda2Dist': triple.models.get(y)[2], 'lambda1Dist': triple.null.get(y)[0],'BIC2Dist': triple.BIC.get(y), 'BIC1Dist': triple.null.get(y)[1], 'count':len(triple.branches(y))})
	


def inputReader(filepath):
	#Reads the input file and passes settings to outputFormatter
	config=ConfigParser.ConfigParser()
	config.read(filepath)
	inputDict={}
	for flag in config.options('Input'):
		inputDict[flag]=config.get('Input',flag)

	if 'treefile' not in inputDict:
		print('Error: No input trees specified.')
		return 
	outputDict={}
	for flag in config.options('Output'):
		outputDict[flag]=config.get('Output',flag)
	outputFormatter(outputDict,inputDict)
		

#print exMax(getTripBranches(readin_Newick(sys.argv[1])), 2, 0.01, 5, 0.5)[0].models
inputReader(sys.argv[1])






