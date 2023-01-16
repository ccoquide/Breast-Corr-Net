import os
import sys
import numpy as np
from UTILS import Read
import matplotlib.pyplot as plt
'''
This python script permits to construct the bi-adjacency matrix showing the nestedness structure of the network nodes. It uses the Fitness and Complexity metrics presented in Tacchella et al. 2012. The analogy with this paper is the following. The Methylated DNA sites play the role of countries and the RNAs play the role of products. Columns and rows are reordered based on the Fitness and Complexity scores. Fitness is a measure associated to DNA nodes and Complexity to RNA nodes.

Inuts:
	-Bi-adjacency matrix file (str)
	-Option (str): "unit" or "weight". We can either consider links between nodes with a unitary weight, or with the correlation coefficients.
	-Option 2 (str): The Fitness and Complexity vectors noramlization that can be either "sum", the 1-norm, or "average", each component is divided by the averaged component.
	-Threshold (float): The threshold used to filter the bi-adjacency matrix elemnts
Outputs:
	-The network of positive and negative correlations, in two separated files.
	-The nodes Fitness and Complexity scores, in two separated files.
	-The two panels plot of the unreordered and reoredred bi-adjacency matrix, in two separated files. Top (bottom) panel is associated to positive (negative) correlations.

Type "python3 FitnessComplexity.py filename unit sum 0.2" to apply this script to the bi-adjacency matrix described in the file filename by considering binary entries (1: link exists, 0: link doesn't exist) keeping all correlations whose absolute value is greater than 0.2. The Fitness and Complexity vector are unit-sum normalized.
	
'''
def buildB(fname,opt,thrs,title_c,title_r):
	'''
	Construct the filtered bi-adjacency matrix from fname filtered using a threshold thrs
	Depending on the threshold, a node of the bipartite network can become an isolated node. We only keep non isolated nodes, therefore we need to reindex them.
	'''
	B=[]
	F=Read(fname,0)
	nc,nr=F[0].split("\t")
	nc=int(nc)
	nr=int(nr)
	COL={}# Dict of column nodes with new indexes. 
	ROW={}# Dict of row nodes with new indexes.
	t_COL={}# column nodes label
	t_ROW={}# row nodes label
	COL["+"]={}# Column nodes involved in positive correlations
	COL["-"]={}#Column nodes involved in negative correlations
	ROW["+"]={}# Row nodes involved in positive corr.
	ROW["-"]={}# Row nodes involved in negative corr.
	t_COL["+"]={}
	t_COL["-"]={}
	t_ROW["+"]={}
	t_ROW["-"]={}
	print(nc, "RNA")
	print(nr, "DNA")
	F.remove(F[0])
	nl=int(F[0])
	F.remove(F[0])
	for i in range(nl):
		tmp=F[i].split("\t")
		v=float(tmp[2])
		if abs(v)>thrs:
			if v > 0.0:
				COL["+"][tmp[0]]=1
				ROW["+"][tmp[1]]=1
			else:
				COL["-"][tmp[0]]=1
				ROW["-"][tmp[1]]=1
	for s in ["+","-"]:
		z=1
		for c in COL[s]:
			t_COL[s][z]=title_c[int(c)-1]# Label of the non isolated nodes
			COL[s][c]=z# Assigning the new index associated to non isolated nodes
			z+=1
	for s in ["+","-"]:
		z=1
		for r in ROW[s]:
			t_ROW[s][z]=title_r[int(r)-1]
			ROW[s][r]=z
			z+=1
	for i in range(nl):
		tmp=F[i].split("\t")
		v=float(tmp[2])
		if abs(v)>thrs:
			if opt == "weight":
				if v> 0.0:
					B.append(str(COL["+"][tmp[0]])+"\t"+str(ROW["+"][tmp[1]])+"\t"+str(v))
				else:
					B.append(str(COL["-"][tmp[0]])+"\t"+str(ROW["-"][tmp[1]])+"\t"+str(v))
			if opt == "unit":
				if v > 0.0:
					B.append(str(COL["+"][tmp[0]])+"\t"+str(ROW["+"][tmp[1]])+"\t"+"1.0")
				else:
					B.append(str(COL["-"][tmp[0]])+"\t"+str(ROW["-"][tmp[1]])+"\t"+"-1.0")
	NC=(len(COL["+"].keys()),len(COL["-"].keys()))# Tuple representing the number of non isolated column nodes involved in positive (negative) correlations
	NR=(len(ROW["+"].keys()),len(ROW["-"].keys()))# The same, but for the row nodes
	del(F)
	return B,NC,NR, t_COL, t_ROW
def conv(F,Q,F_,Q_,it):
	"""
	Convergence test for the power method permitting to find the Fitness and Complexity of the nodes
	Test if (v_n - v_n+1) is close to the null vector
	"""
	eps1=10**-14
	eps2=10**-14
	TST=[False,False]
	cnt1=0.0
	cnt2=0.0
	for i in range(F.shape[0]):
		cnt1+=(F[i]-F_[i])**2
	for i in range(Q.shape[0]):
		cnt2+=(Q[i]-Q_[i])**2
	cnt1=cnt1**0.5
	cnt2=cnt2**0.5
	if ((it == 0) | (it%20==0)):
		print(cnt1, cnt2)
	if cnt1<=eps1:
		TST[0]=True
	if cnt2<=eps2:
		TST[1]=True
	if ((TST[0]==True) & (TST[1]==True)):
		print("Converge at it = ", it)
		return True
	else:
		return False
def RecursiveFQ(B,NC,NR, norm):
	"""
	Power method for computing the Fitness/Complexity metric
	"""
	CONV={}
	SIGN={}
	CONV[True]=False# If it becomes true it means that the recursive process has converged for positive correlations
	CONV[False]=False# In case of negative correlations
	SIGN[True]="+"
	SIGN[False]="-"
	F={}
	F_={}#Tempory Fitness used to test the convergence
	Q={}
	Q_={}#Tempory Fitness used to test the convergence
	if NR[1]!=0:
		F[False]=np.zeros((NR[1],))# Initialization of Fitness metric for DNA for negative corr.
		F_[False]=np.ones((NR[1],))
		F_[False]=F_[False]/F_[False].sum()
	if NR[0]!=0:
		F[True]=np.zeros((NR[0],))# Initialization of Fitness metric for DNA for positive corr.
		F_[True]=np.ones((NR[0],))
		F_[True]=F_[True]/F_[True].sum()
	if NC[1]!=0:
		Q[False]=np.zeros((NC[1],))# Initialization of Complexity metric for RNA for negative corr.
		Q_[False]=np.ones((NC[1],))
		Q_[False]=Q_[False]/Q_[False].sum()
	if NC[0]!=0:
		Q[True]=np.zeros((NC[0],))# Initialization of Complexity metric for RNA for positive corr.
		Q_[True]=np.ones((NC[0],))
		Q_[True]=Q_[True]/Q_[True].sum()
	IT={}# Itreration dict
	IT[True]=0
	IT[False]=0
	PRINT={}
	PRINT[True]=True
	PRINT[False]=True
	it=0
	if norm == "average":
		while((CONV[True] == False) & (CONV[False] == False)):
			if ((it == 0) | (it%20==0)):
				print("it : ", it)
			for l in B:
				tmp=l.split("\t")
				R,D,v=int(tmp[0])-1,int(tmp[1])-1,float(tmp[2])
				F[v>0.0][D]+=Q_[v>0.0][R]*abs(v)
				Q[v>0.0][R]+=(abs(v)/F_[v>0.0][D])
			for q in Q:
				for i in range(Q[q].shape[0]):
					if Q[q][i]!=0:
						Q[q][i]=1/Q[q][i]
				if Q[q].sum()!=0.0:
					Q[q]=Q[q]/(1.0*(Q[q].sum()/Q[q].shape[0]))
			for f in F:
				if F[f].sum()!=0.0:
					F[f]=F[f]/(1.0*(F[f].sum()/F[f].shape[0]))
			if((it==0) | (it%20==0)):
				for s in F:
					print("sum for Fs = ", F[s].sum())
				for s in Q:
					print("sum for Qs = ", Q[s].sum())
			for c in CONV:
				if c in F.keys():
					if CONV[c] == False:
						IT[c]+=1
						CONV[c]=conv(F[c], Q[c], F_[c], Q_[c],it)
						for i in range(F[c].shape[0]):
							F_[c][i]=F[c][i]
						for i in range(Q[c].shape[0]):
							Q_[c][i]=Q[c][i]
					else:
						if PRINT[c]:
							print("Q and F  has converged in case of ", SIGN[c], " in ", IT[c], " iterations")
							PRINT[c]=False
			it+=1
	if norm == "sum":
		while((CONV[True] == False) & (CONV[False] == False)):
			if ((it == 0) | (it%20==0)):
				print("it : ", it)
			for l in B:
				tmp=l.split("\t")
				R,D,v=int(tmp[0])-1,int(tmp[1])-1,float(tmp[2])
				F[v>0.0][D]+=Q_[v>0.0][R]*abs(v)
				Q[v>0.0][R]+=(abs(v)/F_[v>0.0][D])
			for q in Q:
				for i in range(Q[q].shape[0]):
					if Q[q][i]!=0:
						Q[q][i]=1/Q[q][i]
				if Q[q].sum()!=0.0:
					Q[q]=Q[q]/(1.0*(Q[q].sum()))
			for f in F:
				if F[f].sum()!=0.0:
					F[f]=F[f]/(1.0*(F[f].sum()))
			if((it==0) | (it%20==0)):
				for s in F:
					print("sum for Fs = ", F[s].sum())
				for s in Q:
					print("sum for Qs = ", Q[s].sum())
			for c in CONV:
				if c in F.keys():
					if CONV[c] == False:
						IT[c]+=1
						CONV[c]=conv(F[c], Q[c], F_[c], Q_[c],it)
						for i in range(F[c].shape[0]):
							F_[c][i]=F[c][i]
						for i in range(Q[c].shape[0]):
							Q_[c][i]=Q[c][i]
					else:
						if PRINT[c]:
							print("Q and F  has converged in case of ", SIGN[c], " in ", IT[c], " iterations")
							PRINT[c]=False
			it+=1
	print("Q and F (norm = )", norm, "converged for both correlation signs in it = ", IT[True], " for ", SIGN[True], " and it = ", IT[False], " for ", SIGN[False])
	del(F)
	del(Q)
	#Set true zero values
#	for s in F_:
#		for i in range(F_[s].shape[0]):
#			if F_[s][i]<=10**(-14):
#				F_[s][i]=0.0
#	for s in Q_:
#		for i in range(Q_[s].shape[0]):
#			if Q_[s][i]<=10**(-14):
#				Q_[s][i]=0.0
	return F_, Q_
def PrintFQMatrix(B,Q,F,RQ,RF,NC,NR,opt,thrs, t_COL, t_ROW, norm):
	"""
	Create the two panels plot related to the columns and rows reordered bi-adjacency matrix
	input: 
		-The list of filtered links
		-Dictionary of non isolated column nodes associating their index with their Complexity rank
		-Dict. of non isolated row nodes associating their index with their Fitness rank
		-Number of column nodes (for positive / negative correlations)
		-Number of row nodes (for positive/ negative corr.)
		-opt : "unit" or "weight"
		-threshold used to filter the raw matrix
		-Label of non isolated column nodes
		-Label of non isolated row nodes
		-normalisation option: "average" or "sum"
	output:
		-List of links in two separated .net files
		-Pdf file representing the positive correlations matrix and negative correlations matrix
	"""
	ww=open("ADN_ARN_FQ-Matrix-"+opt+"-"+str(thrs)+"-positive.net","w")
	ww2=open("ADN_ARN_FQ-Matrix-"+opt+"-"+str(thrs)+"-negative.net","w")
	ww.write("RNA_ID\tDNA_ID\t"+opt+"\t"+"RNA_Name"+"\t"+"DNA_Name"+"\t"+"Complexity_Score"+"\t"+"Fitness_Score"+"\t"+"\n")
	ww2.write("RNA_ID\tDNA_ID\t"+opt+"\t"+"RNA_Name"+"\t"+"DNA_Name"+"\t"+"Complexity_Score"+"\t"+"Fitness_Score"+"\t"+"\n")
	Mp=np.zeros((NR[0],NC[0]))
	Mn=np.zeros((NR[1],NC[1]))
	for t in B:
		tmp=t.split("\t")
		j_,i_,v=int(tmp[0])-1,int(tmp[1])-1,float(tmp[2])
		if v>0.0:
			i=RF["+"][i_]
			j=RQ["+"][j_]
			Mp[i,j]=v
			ww.write(str(j+1)+"\t"+str(i+1)+"\t"+str(v)+"\t"+t_COL["+"][j_+1]+"\t"+t_ROW["+"][i_+1]+"\t"+str(Q[True][j_])+"\t"+str(F[True][i_])+"\n")
		else:
			i=RF["-"][i_]
			j=RQ["-"][j_]
			Mn[i,j]=abs(v)
			ww2.write(str(j+1)+"\t"+str(i+1)+"\t"+str(v)+"\t"+t_COL["-"][j_+1]+"\t"+t_ROW["-"][i_+1]+"\t"+str(Q[False][j_])+"\t"+str(F[False][i_])+"\n")
	ww.close()
	ww2.close()
	fig, ax = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(10,10))
	fig.tight_layout()
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=1, hspace=1)
	if Mp.shape[0]>0:
		ax1=ax[0].matshow(Mp.transpose(), aspect='auto')# The matrix is transposed for improving the figure reading
	if Mn.shape[0]>0:
		ax2=ax[1].matshow(Mn.transpose(), aspect='auto')
	ax[0].tick_params(axis="y", labelsize=10, pad=2)
	ax[0].tick_params(axis="x", labelsize=10, pad=2)
	ax[0].set_xlabel(r'RNA', fontsize=15)
	ax[0].set_ylabel(r'DNA', fontsize=15)
	ax[0].annotate(r'Corr. $+$', xy=(0.8,0.8), xycoords='axes fraction', color="black",  bbox=dict(facecolor='white', edgecolor='black'), fontsize="15")
	ax[0].annotate(r'Threshold = '+str(thrs)+", norm = "+norm, xy=(0.2,0.8), xycoords='axes fraction', color="black",  bbox=dict(facecolor='white', edgecolor='black'), fontsize="15")
	ax[1].tick_params(axis="y", labelsize=10, pad=2)
	ax[1].tick_params(axis="x", labelsize=10, pad=2)
	ax[1].set_xlabel(r'RNA', fontsize=15)
	ax[1].set_ylabel(r'DNA', fontsize=15)
	ax[1].annotate(r'Corr. $-$', xy=(0.8,0.8), xycoords='axes fraction', color="black",  bbox=dict(facecolor='white', edgecolor='black'), fontsize="15")
	plt.savefig("ADN-ARN_FQ-Matrix-"+norm+"-"+opt+"-"+str(thrs)+".pdf", bbox_inches = 'tight', pad_inches = 0)
def PrintMatrix(B,NC,NR,opt, thrs):
	"""
	Create the two panels plot related to the raw bi-adjacency matrix
	input:
	        -The list of filtered links
	        -Number of column nodes (for positive / negative correlations)
	        -Number of row nodes (for positive/ negative corr.)
	        -opt : "unit" or "weight"
	        -threshold used to filter the raw matrix
	output:
	        -Pdf file representing the positive correlations matrix and negative correlations matrix
	"""
	Mp=np.zeros((NR[0],NC[0]))
	Mn=np.zeros((NR[1],NC[1]))
	for t in B:
		tmp=t.split("\t")
		j_,i_,v=int(tmp[0])-1,int(tmp[1])-1,float(tmp[2])
		if v > 0.0:
			Mp[i_,j_]=v
		else:
			Mn[i_,j_]=abs(v)
	fig, ax = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(10,10))
	fig.tight_layout()
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=1, hspace=1)
	if Mp.shape[0]>0:
		ax1=ax[0].matshow(Mp.transpose(), aspect='auto')
	if Mn.shape[0]>0:
		ax2=ax[1].matshow(Mn.transpose(), aspect='auto')
	ax[0].tick_params(axis="y", labelsize=10, pad=2)
	ax[0].tick_params(axis="x", labelsize=10, pad=2)
	ax[0].set_xlabel(r'RNA', fontsize=15)
	ax[0].set_ylabel(r'DNA', fontsize=15)
	ax[0].annotate(r'Corr. +', xy=(0.8,0.8), xycoords='axes fraction', color="black",  bbox=dict(facecolor='white', edgecolor='black'), fontsize="15")
	ax[0].annotate(r'Threshold = '+str(thrs), xy=(0.2,0.8), xycoords='axes fraction', color="black",  bbox=dict(facecolor='white', edgecolor='black'), fontsize="15")
	ax[1].tick_params(axis="y", labelsize=10, pad=2)
	ax[1].tick_params(axis="x", labelsize=10, pad=2)
	ax[1].set_xlabel(r'RNA', fontsize=15)
	ax[1].set_ylabel(r'DNA', fontsize=15)
	ax[1].annotate(r'Corr. -', xy=(0.8,0.8), xycoords='axes fraction', color="black",  bbox=dict(facecolor='white', edgecolor='black'), fontsize="15")
	plt.savefig("ADN-ARN_Matrix-"+opt+"-"+str(thrs)+".pdf", bbox_inches = 'tight', pad_inches = 0)
#############
#           #
#############
fname=sys.argv[1]# filename of the row bi-adjacency matrix
#The first two lines of the file must contains the number of columns and rows respectively. Then the second line contains only the number of links. All other lines describe links. The two first elements are associated to the column and row node index respectively. Note that nodes indexes are in the range (1,N) in the file. The last element is the correlation coefficient associated to this link
opt=sys.argv[2]# Option 1 : "unit" or "weight"
norm=sys.argv[3]# Option2 : "sum" or "average" vector normalization
if opt not in ["unit", "weight"]:
	print("Please choose an option in [\"unit\", \"weight\"]")
	exit()
if norm not in ["average", "sum"]:
	print("Please choose a normalization option in [\"average\",\"sum\"]")
	exit()
thrs=float(sys.argv[4])#Threshold used for the raw matrix filtering
title_c=Read("ARN_titles.dat",0)
title_r=Read("ADN_titles.dat",0)
B,NC,NR,t_COL,t_ROW=buildB(fname,opt,thrs,title_c,title_r)# Extraction of the filtered list of links with associated information
F,Q=RecursiveFQ(B,NC,NR,norm)# Computing the Fitness and Complexity of DNA and RNA.
RF={}#Rank of nodes in terms of Fitness score
RQ={}#Rank of nodes in terms of Complexit score
tmpF={}
tmpQ={}
if True in F.keys():
	RF["+"]={}
	tmpF["+"]=sorted(range(F[True].shape[0]), key = lambda k: F[True][k], reverse=True)
if False in F.keys():
	RF["-"]={}
	tmpF["-"]=sorted(range(F[False].shape[0]), key = lambda k: F[False][k], reverse=True)
if True in Q.keys():
	RQ["+"]={}
	tmpQ["+"]=sorted(range(Q[True].shape[0]), key = lambda k: Q[True][k], reverse=True)
if False in Q.keys():
	RQ["-"]={}
	tmpQ["-"]=sorted(range(Q[False].shape[0]), key = lambda k: Q[False][k], reverse=True)
for s in tmpF:
	z=0
	for i in tmpF[s]:
		RF[s][i]=z
		z+=1
for s in tmpQ:
	z=0
	for i in tmpQ[s]:
		RQ[s][i]=z
		z+=1
ww=open("Fitness-DNA-"+opt+"-"+str(thrs)+".dat","w")# Output consisting in the list of row nodes with their associated Fitness score and rank
ww2=open("Complexity-RNA-"+opt+"-"+str(thrs)+".dat","w")# Output consisting in the list of column nodes with their associated Complexity socre and rank
ww.write("DNA\tFitness Score\tFitness Rank\tSign\n")
ww2.write("RNA\tComplexity Score\tComplexity Rank\tSign\n")
for k in range(len(tmpF["+"])):
	i_=tmpF["+"][k]
	ww.write(t_ROW["+"][i_+1] + "\t" + str(F[True][i_]) + "\t" + str(k+1)+"\t"+"+"+"\n")
for k in range(len(tmpF["-"])):
        i_=tmpF["-"][k]
        ww.write(t_ROW["-"][i_+1] + "\t" + str(F[False][i_]) + "\t" + str(k+1)+"\t"+"-"+"\n")
ww.close()
for k in range(len(tmpQ["+"])):
        i_=tmpQ["+"][k]
        ww2.write(t_COL["+"][i_+1] + "\t" + str(Q[True][i_]) + "\t" + str(k+1)+"\t"+"+"+"\n")
for k in range(len(tmpF["-"])):
        i_=tmpQ["-"][k]
        ww2.write(t_COL["-"][i_+1] + "\t" + str(Q[False][i_]) + "\t" + str(k+1)+"\t"+"-"+"\n")
ww2.close()
PrintMatrix(B,NC,NR,opt,thrs)# Plotting the raw matrix
PrintFQMatrix(B,Q,F,RQ,RF,NC,NR,opt,thrs, t_COL, t_ROW, norm)# Plotting the column and row reorganized matrix 
