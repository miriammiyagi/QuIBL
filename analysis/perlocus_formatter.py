#Formats output from QuIBL with a file with a list of loci to give locus-per-locus introgression estimates given a triplet.
#Input form is: python perlocus_formatter.py treefile_path quibl_out_csv_path triplet output_filepath canonical_outgroup
#Michael Miyagi
#6/6/2019

import csv
import sys
from ete3 import Tree
import numpy as np


def getModel(path,testTrip):
	out=[]
	with open(path) as modelfile:
		models=csv.DictReader(modelfile)
		for row in models:
			if row['triplet']==testTrip:
				out.append(row)
	return out

def conPDF(x, C, lmbd):
	if x>0 and x>=C*lmbd:
		y=0.5/lmbd*np.exp(-x/lmbd)*(1+np.exp(C))
	elif x>0 and C*lmbd>x:
		y=np.sinh(x/lmbd)/(lmbd*(-1+np.exp(C)))
	else:
		y=0
	return y

def readin_Newick(filepath):
	treeList=[]
	returnlist=[]
	f=open(filepath).read()
	composite=f.split(';')
	for item in composite:
		if item.strip():
			returnlist.append((item+';').strip('\n'))
	for tree in returnlist:
		treeList.append(Tree(tree))
	return treeList,returnlist

def getTripBranches(treeList,tripOfInt,model,canonOut):
	output=[]#[[tripOfInt[0]],[tripOfInt[1]],[tripOfInt[2]]]
	refMod=[]
	for outg in tripOfInt:
		for mod in model:
			if mod['outgroup']==outg:
				refMod.append(mod)
	#dist=0
	for tind,tree in enumerate(treeList):
		tree.set_outgroup(canonOut)
		#print tree.expand_polytomies()
		if len(tree.expand_polytomies())>1:
			print 'Tree '+str(tind)+' skipped due to polytomy.'
			continue
		dist=0
		tempTree=tree.copy('newick')
		tempTree.prune(tripOfInt,preserve_branch_length=True)
		common=tempTree.get_leaves()[0].get_common_ancestor(tempTree.get_leaves())
		for x in common.get_children():
			if not (x.is_leaf()):
				dist=x.get_distance(common)
			else:
				tempInd=tripOfInt.index(x.get_leaf_names()[0])
		if dist>0:
			pint=float(refMod[tempInd]['mixprop2'])*conPDF(dist, float(refMod[tempInd]['C2']), float(refMod[tempInd]['lambda2Dist']))
			pils=float(refMod[tempInd]['mixprop1'])*conPDF(dist, float(refMod[tempInd]['C1']), float(refMod[tempInd]['lambda2Dist']))
			iprob=pint/(pint+pils)
			output.append((tripOfInt[tempInd],tind,dist,iprob))
	return output

def writeOut(treesray,inpath,tripOfInt,outpath,canonOut):
	treelist,textrees=readin_Newick(treesray)
	model=getModel(inpath,tripOfInt)
	tripList=tripOfInt.split("_")
	output=getTripBranches(treelist,tripList,model,canonOut)
	with open(outpath,'w') as csv_out:
		fieldnames=['tree','outgroup','branchLength', 'intro_prob']
		out_writer=csv.DictWriter(csv_out,quoting=csv.QUOTE_MINIMAL, fieldnames=fieldnames)
		out_writer.writeheader()
		for locus in output:
			out_writer.writerow({'tree':textrees[locus[1]],'outgroup':locus[0],'branchLength':locus[2], 'intro_prob': locus[3]})
		
writeOut(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])


