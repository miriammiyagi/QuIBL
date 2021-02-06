#Quantifying Introgression via Branch Lengths - Cython Version
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
import cyquibl_mod as cyq

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

def tripSetGen(treeSet):
#Generates the set of triplets to check.
	leaves=treeSet[0].get_leaf_names()
	return list(itt.combinations(leaves,3))


def getTripBranches(treeList,canonOut):
#Goes through the trees in the input and calculates/sorts the internal branch lengths.
	setOfTriplets=[]
	triples=tripSetGen(treeList)
	lenIt=sum(1 for _ in triples)
	output=np.zeros((4,lenIt))
	output=[[None,[],[],[]] for i in range(lenIt)]
	for counter,tree in enumerate(treeList):
		tree.set_outgroup(canonOut)
		if len(tree.expand_polytomies())>1:
			print('Tree '+str(counter)+' skipped due to polytomy.')
			continue
		dist=0
		for index,triplet in enumerate(triples):
			output[index][0]=triplet
			tempTree=tree.copy('newick')
			tempTree.prune(triplet,preserve_branch_length=True)
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
				XQList=cyq.XQUpdate(branchData, cArray, pArray, lmbd)
				pArray=cyq.pUpdate(XQList)
				for indy in range(1,K):
					cArray[indy]=cyq.cSearch(branchData, XQList, indy, lmbd, cArray)
				storedlmbd=lmbd
				lmbd,stepScale=cyq.gradAscent(branchData, XQList, cArray, pArray, lmbd, stepScale, numSteps, threshold)
				cArray=[x*(storedlmbd/lmbd) for x in cArray]
				steps+=1
			triple.setModel(outG, list(cArray), list(pArray), lmbd)
			triple.setBIC(outG, np.log(len(branchData))*(2*K-1)-2*(cyq.modelLogLik(branchData, list(cArray), list(pArray), lmbd)))
			triple.setNull(outG,oneDistNull(branchData))
	return tripletSet

def PLexMax(triple, K, threshold, numSteps, tempScale):
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
			XQList=cyq.XQUpdate(branchData, cArray, pArray, lmbd)
			pArray=cyq.pUpdate(XQList)
			for indy in range(1,K):
				cArray[indy]=cyq.cSearch(branchData, XQList, indy, lmbd, cArray)
			storedlmbd=lmbd
			lmbd,stepScale=cyq.gradAscent(branchData, XQList, cArray, pArray, lmbd, stepScale, numSteps, threshold)
			cArray=[x*(storedlmbd/lmbd) for x in cArray]
			steps+=1
		triple.setModel(outG, list(cArray), list(pArray), lmbd)
		triple.setBIC(outG, np.log(len(branchData))*(2*K-1)-2*(cyq.modelLogLik(branchData, list(cArray), list(pArray), lmbd)))
		triple.setNull(outG,oneDistNull(branchData))
	return triple


def oneDistNull(branchData):
	lmbd=np.mean(branchData)
	bic=np.log(len(branchData))-2*(cyq.modelLogLik(branchData, [0.0], [1.0], lmbd))
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

inputReader(sys.argv[1])
