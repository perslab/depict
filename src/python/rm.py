#! /usr/bin/python

import os

path = "/home/projects/depict/git/DEPICT/data/backgrounds"
with open("/tmp/ll2",'r') as infile:
	for line in infile.readlines():
		os.system("rm {}/{}/*clumped".format(path,line.strip()))	
		os.system("rm {}/{}/*nosex".format(path,line.strip()))	
		os.system("rm {}/{}/*log".format(path,line.strip()))	
		#os.system("gzip {}/{}/*txt".format(path,line.strip()))

