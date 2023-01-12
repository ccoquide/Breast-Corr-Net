import os
import sys
def Read(fname,nrow):
	rez=open(fname,"r").read().split("\n")
	rez.remove("")
	return rez[nrow:]
