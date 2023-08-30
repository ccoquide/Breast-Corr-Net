import os
import sys
from UTILS import Read
"""
Transform a matrix from R write function into a tsv list of matrix element [col ID] [row ID] [value]
Input:
	1) file with names of col indices as first row, then m rows of lenght n+1, where m is the num. of rows and m the num. of columns. First element of each row is the name of the corresponding row index
	2) Filename of the output matrix
"""
fname=sys.argv[1]# Filename of input
output=sys.argv[2]#Filename of the output
M=Read(fname,0)# List of column index names
print(M[0])
COL_NAMES=M[0].split("\t")
ROW_NAMES=[]
M.remove(M[0])
i=0
ww=open(output,"w")
OUT=[]
ncol=0
nrow=0
LIST=[]#Contains all ADN-ADN pairs such that their pearson correlation >0.5
for l in M:
	ROW_NAMES.append(l.split("\t")[0])
	tmp=l.split("\t")[1:]
	ncol=len(tmp)
	for j in range(len(tmp)):
#		if j+1!=i+1:
		OUT.append(str(j+1)+"\t"+str(i+1)+"\t"+tmp[j]+"\n")
		if abs(float(tmp[j]))>0.5:
			LIST.append((str(j+1), str(i+1), tmp[j]))
	i+=1
nrow=i
print("Matrix with ", nrow, " rows, and", ncol, " columns")
ww.write(str(ncol)+"\t"+str(nrow)+"\n")
nl=len(OUT)# Number of non-negative matrix element
print("we have ", str(nl), "non-zero matrix elements.")
ww.write(str(nl)+"\n")
for l in OUT:
	ww.write(l)
del(OUT)
ww.close()
ww1=open("ADN_NAMES.dat","w")
ww2=open("ARN_NAMES.dat","w")
for n in COL_NAMES:
	ww2.write(n+"\n")
for n in ROW_NAMES:
	ww1.write(n+"\n")
ww1.close()
ww2.close()
print(len(LIST))
for i in LIST:
	print(i)
