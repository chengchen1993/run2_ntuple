#!/usr/bin/python
import time
import os

for ii in range(0,100000000):
  file=os.listdir("./")
  filetmp=os.listdir("./")
  wrong=[]
  failed=[]
  SUBMITFAILED=[]
  
  for f in filetmp:
  	if "crab_" not in f or ".sh" in f or "_R0-" in f: file.remove(f);
  for f in file:
  	os.system("crab status "+f +"> tmp.log")
  	lines=[]
  	filelog=open("tmp.log")
  	for line in filelog:
  		lines.append(line);
        if len(lines)>10:
        	if "100.0%" not in lines[10] or "running" in lines[10]:wrong.append(f);
        	if "failed" in lines[9] or "failed" in lines[10] or "RESUBMITFAILED" in lines[3]:failed.append(f)
        if len(lines)>3:
                if "SUBMITFAILED" in lines[3]:SUBMITFAILED.append(f)
  	os.system("rm tmp.log")
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  print "\n"
  for i  in failed:
          os.system("crab resubmit %s"%(i))
  #os.system("sleep 3 m")
  for i  in wrong:
          os.system("crab status %s"%(i))
  time.sleep( 600 )
