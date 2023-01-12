import os
import sys
def Read(fname,nrow):
	rez=open(fname,"r").read().split("\n")
	rez.remove("")
	return rez[nrow:]
def RemoveVoid(M):
	"""
	Print a matrix without empty rows/columns
	"""
	rIDX={}# Dict for old - new index conversion
	rIDX[0]={}
	rIDX[1]={}
	z=0
	for i in range(M.shape[0]):
		if sum(M[i])!=0.0:
			rIDX[0][z]=i
			z+=1
	z=0
	for i in range(M.shape[1]):
		if sum(M.transpose()[i])!=0.0:
			rIDX[1][z]=i
			z+=1
	nr_,nc_=len(rIDX[0].keys()),len(rIDX[1].keys())
	M_=np.zeros((nr_,nc_))
	for i in range(nr_):
		i_=rIDX[0][i]
		for j in range(nc_):
			j_=rIDX[1][j]
			M_[i,j]=M[i_,j_]
	return M_
