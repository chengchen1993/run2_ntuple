import os
file=os.listdir("./")
filetmp=os.listdir("./")
wrong=[]
failed=[]
SUBMITFAILED=[]

for f in filetmp:
	if "crab_" not in f or ".sh" in f or "_R0-" not in f: file.remove(f);
print file
for f in file:
	os.system("crab status "+f +"> tmp2.log")
        print f
	lines=[]
	filelog=open("tmp2.log")
	for line in filelog:
                print line
		lines.append(line);
#	if "100.0%" not in lines[10] or "running" in lines[10]:wrong.append(f);
#        if "100.0%" not in lines[len(lines)-1] or "running" in lines[len(lines)-1]:wrong.append(f);
#        rm=0
 #       if len(lines)>10:
  #             if "finished     		100.0%" in lines[9] or "finished     		100.0%" in lines[10]: 
   #                         rm=1
    #                        os.system("rm -rf "+f)
     #                       file.remove(f)
      #                      continue
        rm_and_submit=0

        if len(lines)<=3:
                  rm_and_submit=1
        if len(lines)>3 and len(lines)<=8:
                  if "SUBMITFAILED" in lines[3] or "SUBMITFAILED" in lines[4]: rm_and_submit=1
        if len(lines)>8:
                  if "Cannot retrieve the status" in lines[8] or "SUBMITFAILED" in lines[3] or "SUBMITFAILED" in lines[4] or "is not responding" in lines[9]: rm_and_submit=1

        resubmit=0
        if len(lines)>3 and len(lines)<=10:
                  if "RESUBMITFAILED" in lines[3]: resubmit=1
        if len(lines)>11:
                  if "failed" in lines[9] or "failed" in lines[10]  or "failed" in lines[11] or "RESUBMITFAILED" in lines[3]: resubmit=1

        if rm_and_submit==1:
                           SUBMITFAILED.append(f)
        else:
                           if resubmit==1: failed.append(f)
        os.system("rm tmp2.log")
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
for i  in SUBMITFAILED:
        os.system("rm -rf %s"%(i))
        os.system("crab submit %s"%('crab3_analysis'+i[5:12]+'_'+(i[12:])[:-4]+'.py'))
os.system("sleep 1 m")
for i  in failed:
        os.system("crab resubmit %s"%(i))
os.system("sleep 1 m")
for i  in file:
        os.system("crab status %s"%(i))
