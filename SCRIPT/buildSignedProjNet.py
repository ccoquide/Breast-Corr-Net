import os
import sys
from UTILS import Read
import numpy as np
from numpy import array
from scipy import sparse
from collections import defaultdict as ddict
"""
This script constructs matrices corresponding to the column (RNA) and row (DNA) projections of the related bipartite network. All the corresponding networks, built from projections and the initial bipartite network, are link signed. It means that any I-J link has a signed weight +/- rho(I,J) where rho(A,B) is the pearson correlation coeffcient related to DNA I and RNA J. 
Input:
	1) Matrix element [col index] [row index] [Pearson corr. coeff.]
	2) Threshold sigma such that only matrix element equal or greater than sigma are kept and used for the projections
Output: DNA and RNA Projections with following format
[col index] [row index] [positive weight] [negative weight]
"""
#We want to build projectied mode-1 network from a bipartite (directed) network
#We consider here signed bipartite network. In such a context, if I-a and J_a are two positive links, then the path I-J will have the sign of the product. So if I,J share a node a, such that I-a is positive and J-a is negative, then the projection I-J will be negative. I use complex matrix in order to store information about sign. B_{ij}=a+ib where a is the positivity and b the negativity. B_{ij}=1 means a positive link is present and B_{ij}=i means we have a negative links. In this way, the adjacency matrix A=BB^{T} is such that abs(Re(A_{ij}))= number of positive links, and abs(Im'A_{ij}) is the number of negative links
def Sign(f):
	"""
	function taking as input a signed float, which is a correlation coeff. and give as output a complex number such that z=1+0.i encodes a positivie correlation coeff. and z=0+1.i a negative correlation coeff.
	"""
	if (abs(f)/f)>0:
		return 1+0J
	else:
		return 0+1J
def Sign_bis(f):
	if (abs(f)/f)>0:
		return abs(f)+0J
	else:
		return 0+abs(f)*1J
def SignedBipN_Builder(netfn, thrs, thrs_opt):
	"""
	This function builds the complex matrix representing the signe bipartite network
	Input:
		1) Filename of the matrix element file, whose element are either positive or negative Pearson correlations coeffcients
		2) Threshold sigma
	"""
	l=Read(netfn,0)
	tmp=l[0].split("\t")# First row the file contains ncol and nrow of the matrix
	ncol=int(tmp[0])
	nrow=int(tmp[1])
	shape=(nrow,ncol)
	l.remove(l[0])
	Nl=int(l[0])# Number of non-zero matrix element
	l.remove(l[0])
	J=[]# List of col index such that B[,j] is non-zero
	I=[]# List of row index such that B[i,] is non-zero
	V=[]# List of matrix non-zero element
	trans_proof_opt=False# Default value
	if thrs_opt == False:
		thrs_=float(thrs)
	if thrs_opt == True:
		if thrs == "brute-force":
			thrs_=1./(2**0.5)
		if thrs == "trans-proof":
			trans_proof_opt=True
	print("(N_1,N_2) = ",shape)
	for e in l:
		tmp=e.split("\t")
		if trans_proof_opt==False:
			if abs(float(tmp[2]))>thrs_:# Threshold filtering
				J.append(int(tmp[0])-1)
				I.append(int(tmp[1])-1)
				V.append(Sign(float(tmp[2])))# Sign function transform a positive/negative correlation coeff. into a complex number
		else:
			J.append(int(tmp[0])-1)
			I.append(int(tmp[1])-1)
			V.append(Sign_bis(float(tmp[2])))
	if trans_proof_opt==False:
		return sparse.coo_matrix((V,(I,J)),shape=shape), ncol, nrow
	else:
		return [J, I, V], ncol, nrow# return a sparse complex matrix object encoding the signed bipartite network
def M_prod(M_1,M_2):
	return M_1*M_2
def AltProd(J,I,V,JB,IB,VB,n):
	"""
		This function takes as in put two square matrices A, and B of size (n,m) and (m,n) respectively and in output the Matrix M issued from an alternative Matrix product such that Every partial element A[i,k]*B[k,j] contributing for M[i,j] element are filtered such that A[i,k]**2 + B[j,j]**2 > 1. This ensures the correlation transitivity.
	"""
	J_=[]# J indices for non-zero elements of the final matrix
	I_=[]# I indices for non-zero elements of the final matrix
	V_=[]# Non-zero element of the final matrix
	print(n)
	print(max(JB), max(IB))
	for k in range(len(J)):
		i,j,v=I[k],J[k],V[k]
		tmpv=[0+0J]*n
		for k_ in range(len(IB)):
			i_,j_,v_=IB[k_],JB[k_],VB[k_]
			if j==i_:
				vr,vi=abs(np.real(v)),abs(np.imag(v))
				v_r,v_i=abs(np.real(v_)),abs(np.imag(v_))
				if(((vr+vi)**2+(v_r+v_i)**2)>1.0):
					tmpv[j_]+=(v/max(vr,vi))*(v_/max(v_r,v_i))
		for j_ in range(len(tmpv)):
			if(tmpv[j_] != 0+0J):
				I_.append(i)
				J_.append(j_)
				V_.append(tmpv[j_])
		del(i_)
		del(j_)
		del(v_)
		del(tmpv)
	return sparse.coo_matrix((V_,(I_,J_)),shape=(n,n))

def Proj(M):
	"""
	Simple function computing either the row or column projection from a complex square matrix
	"""
	tmp=np.dot(M,M.transpose())# Compute the part of the projection where positive signed links are conserved
	tmp_=np.dot(M,(M.transpose().conjugate()))# Compute the part of the projection where negative signed links are conserved
	return np.real(tmp_)+np.imag(tmp)*1j
def AltProj(J,I,V,n):
	tmp=AltProd(J,I,V,I,J,V,n)# Compute the part of the projection where positive signed links are conserved
	V_=[]
	for v in V:
		V_.append(v.conjugate())
	tmp_=AltProd(J,I,V,I,J,V_,n)# Compute the part of the projection where negative signed links are conserved
	return np.real(tmp_)+np.imag(tmp)*1j	
def ProjectTo(space, M, ncol, nrow, thrs):
	"""
	This function takes as input, a string related to the space we are projecting onto
	it is either 'row' or 'col'
	"""
	if thrs!= "trans-proof":
		if space == 'row':
			return Proj(M)# Compute the projection and
		if space == 'col':
			return Proj(M.transpose())# Compute the projection, here the transpose is required
		if space not in ["row", "col"]:
			print("Error with the space projection")
			print("please choose a space option string in [\"row\", \"col\"]")
			exit()
	if thrs=="trans-proof":
		if space == 'row':
			return AltProj(M[0],M[1],M[2], nrow)# Compute the projection and
		if space == 'col':
			return AltProj(M[1], M[0], M[2], ncol)# Compute the projection, here the transpose is required
		if space not in ["row", "col"]:
			print("Error with the space projection")
			print("please choose a space option string in [\"row\", \"col\"]")
			exit()
def GetNz(M):
	"""
	Take the upper triangle sub matrix with non-zero element
	"""
	print("Extracting upper triangle sub matrix...")
	return sparse.find(sparse.triu(M,1))
	print("done")
def WriteProj(fname,M):
	"""
	This function takes as input the output filename and the complex matrix M describing one of the possible projections
	"""
	nN=0#Sum of weight of negative signed links
	nP=0#Sum of weight of positive signed links
	nM=0#Sum of mixed signed links weight. A mixed signed link is when pair i,j has both a neagtive and positive signed links
	ww=open(fname,"w")
	print("Extracting non-zero elements location")
	i,j,v=GetNz(M)# This function return the lists i, j and v, encoding for row indices, col indices describing upper-triangle non-zero matrix element, and the list of the corresponding complex values
	tmps=sorted(range(len(i)), key = lambda k: (j[k], i[k]))# This permits to print j,i indices sorted in increasing ordering
##########################
##Progression bar setting#
##########################
	if len(i)>=1000:
		nn=len(i)//100
	else:
		nn=1000
##########################
##Progression bar setting#
##########################
	n=M.get_shape()[0]# Dimension of the matrix
	ww.write(str(n)+"\n")
	ww.write(str(len(i))+"\n")#Number of non-zero matrix elements
	nLp=0
	nLn=0
	for link in v:
		if np.real(link)>0:
			nLp+=1# Counter for number of positive links
		if np.imag(link)>0:
			nLn+=1# Counter for number of negative links
	ww.write(str(nLp)+"\t"+str(nLn)+"\n")# Adding this information in the output matrix element file
	print("Start writing in", fname)
	NS=ddict(lambda: 0)#Dictionnary of nodes : [Node ID] -> integer in {1,2}. 1 is denotes a node involved in only negative/positive links, 2 denotes a node involved in at least one both signed links in the projection
	for x in range(len(tmps)):
		if x == 0 or x%nn == 0:
			print(x/nn, "%", end="\r")# Progression bar (percentage of matrix completion)		    
		z=tmps[x]
		tmp=v[z]# Matrix element
		ww.write(str(j[z]+1)+"\t"+str(i[z]+1)+"\t"+str(np.real(tmp))+"\t"+str(np.imag(tmp))+"\n")# writting the matrix element in the output file
		if np.imag(tmp)*np.real(tmp)!=0:
			nM+=np.imag(tmp)+np.real(tmp)# Counter for Mixed signed links weight
		if np.imag(tmp)!=0:# Counter for the number of negative signed links nodes are involved in
			NS[j[z]+1]+=1J# We consider both nodes since the network is undirected
			NS[i[z]+1]+=1J
		if np.real(tmp)!=0:# Counter for the number of positive signed links nodes are involved in
			NS[i[z]+1]+=1# We consider both nodes since the network is undirected
			NS[j[z]+1]+=1
		nP+=np.real(tmp)
		nN+=np.imag(tmp)
	cnt=0
	for n in NS:
		if ((np.real(NS[n])!=0) + (np.imag(NS[n])!=0) == 2):
			cnt+=1# Counter of nodes involved in mixed signed links
#	print(x/nn, "%")
	print("done")
	ww.close()
	print("Writing .info")
	print(nN, "total negative weight (nN)")
	print(nP, "total positive weight (nP)")
	print(nM, "mixed weight")
	print(cnt, "nodes involved in both negative and positive networks")
	print(len(i), "links in total")
	ww=open(fname.replace(".dat", ".info"),"w")
	ww.write("nN="+str(nN)+"\n"+"nP="+str(nP)+"\n"+"nM="+str(nM)+"\n"+str(cnt)+"nodes involved in both neg/pos networks\n")
	ww.write("nL="+str(len(i))+"\n")
	ww.close()
###########################################
###########################################
print("Constructing Signed Bipartite Matrix...")
net=sys.argv[1]
thrs=sys.argv[2]
thrs_opt=False# Default value, if true it means we don't use a decimal threshold but an option
#Option brute-force set a threshold a 1/sqrt(2), in this way we are sure about correlations transitivity
#Option trans-proof will take into account only pairs of correlations such that the sum of their square is greater than 1
if thrs[1]!=".":# Means we give a threshold option
	if thrs not in ["brute-force", "trans-proof"]:
		print("Please choose a decimal threshold in [0,1] or an option in [\"brute-force\",\"trans-proof\"]")
		exit()
	else:
		thrs_opt=True
if thrs_opt==False:
	print("thrs = ", thrs)
else:
	print("opt = ", thrs)
B,ncol,nrow=SignedBipN_Builder(net, thrs, thrs_opt)#	Building Signed Bipartite adjacency matrix from the matrix element file
if thrs=="trans-proof":
	print(B[0][:10])
	print(B[1][:10])
	print(B[2][:10])
print("done")
print("ROW-space projection (ADN) ...")
A=ProjectTo('row', B, ncol, nrow, thrs)# Prohection onto the row space
WriteProj(net.replace(".dat","-"+thrs+"-proj-ADN.dat"), A)# Writing the output file
print(nrow, "nodes in the ROW-space projection")
print("COL-space projection (ARN) ...")
A=ProjectTo('col', B, ncol, nrow, thrs)# Projection onto the col space
WriteProj(net.replace(".dat", "-"+thrs+"-proj-ARN.dat"), A)# Writing the output file
print(ncol, "nodes in the COL-space projection")
print("done")
