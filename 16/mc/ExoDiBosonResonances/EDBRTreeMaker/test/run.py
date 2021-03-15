#!/usr/bin/python
import time
import os

for ii in range(0,1):
  file=os.listdir("./")
  filetmp=os.listdir("./")
  wrong=[]
  failed=[]
  SUBMITFAILED=[]
  
  for f in filetmp:
  	if "crab_" not in f or ".sh" in f or "_R0-" not in f: file.remove(f);
  for f in file:
  	os.system("crab status "+f)
