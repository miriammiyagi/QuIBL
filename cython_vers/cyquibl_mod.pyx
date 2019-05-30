import numpy as np
from libc.math cimport exp,sinh,cosh,log,pow
cimport numpy as np

def conPDF(double x, double C, double lmbd):
	cdef double y
	if x>0 and x>=C*lmbd:
		y=0.5/lmbd*exp(-x/lmbd)*(1+exp(C))
	elif x>0 and C*lmbd>x:
		y=sinh(x/lmbd)/(lmbd*(-1+exp(C)))
	else:
		y=0
	return y

def qCalc(double x, list cArray, list pArray, double lmbd):
	cdef list xq=[]
	for c in range(0,len(cArray)):
		xq.append(pArray[c]*conPDF(x,cArray[c],lmbd))
	xq=[l/sum(xq) for l in xq]
	return xq

def modelLogLik(list XSet, list cArray, list pArray, double lmbd):
	cdef double tempSum=0
	for x in XSet:
		L=0
		for i in range(0,len(pArray)):
			L+=pArray[i]*conPDF(x,cArray[i],lmbd)
		tempSum+=log(L)
	return tempSum

def weightedLogLik(list XQList, list cArray, double lmbd):
	cdef double L=0
	for tup in XQList:
		temp=0
		for indi,pival in enumerate(tup[1]):
			temp+=pival*(conPDF(tup[0],cArray[indi],lmbd))
		L+=log(temp)
	return L
		
def cSearch(list XSet, list XQList, int index, double lmbd, list cArray):
	cdef list rescaledXSet
	cdef list tempLik=[]
	cdef list testCArray=cArray[:]
	rescaledXSet=[x/lmbd for x in XSet]
	for j,testVal in enumerate(rescaledXSet):
		testCArray[index]=testVal
		tempLik.append(weightedLogLik(XQList, testCArray, lmbd))
		testCArray[index]=cArray[index]
	op=rescaledXSet[np.argmax(tempLik)]
	return op

def pUpdate(list XQList):
	cdef list newpArr=[]
	for i in range(0,len(XQList[1])):
		temp=[]
		for x in XQList:
			temp.append(x[1][i])
		newpArr.append(np.mean(temp))
	return newpArr

def XQUpdate(list XSet, list cArray, list pArray, double lmbd):
	cdef list XQList=[]
	for x in XSet:
		outXQ=[x]
		outXQ.append(qCalc(x, cArray, pArray,lmbd))
		XQList.append(outXQ)
	return XQList

def calcDeriv(XQList, cArray, lmbd):
	cdef double val
	cdef double alpha
	cdef double beta
	cdef double deriv=0
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
				p=x[1][1]
				q=x[1][0]
				deriv+=(val+(beta/(-1+exp(beta/lmbd)))-(exp(2*val/lmbd)*p*(2*val-beta)+(p+2*q)*beta)/((-1+exp(2*val/lmbd))*p+2*(-1+exp(beta/lmbd))*q)-lmbd)/pow(lmbd,2)
			if val>=beta:
				p=x[1][1]
				q=x[1][0]
				deriv+=(x[0]-beta+((p+2*q)*beta)/(p+exp(beta/lmbd)*p+2*q)-lmbd)/pow(lmbd,2)
		return deriv

def gradAscent(list XSet,list XQList,list cArray,list pArray,double lmbd,double stepScale,int numSteps, double threshold):
	if len(cArray)==1:
		return np.mean(XSet), stepScale
	elif len(cArray)==2:
		likArray=[0]
		likArray.append(weightedLogLik(XQList, cArray, lmbd))
		sc=0
		while abs(likArray[-1]-likArray[-2])>=threshold*abs(likArray[-1]) and sc<numSteps:
			newlmbd=lmbd+calcDeriv(XQList, cArray, lmbd)*stepScale
			if newlmbd<=0:
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
		print 'Error: C array is too long- only K=1 or K=2 currently supported.'
		return None



