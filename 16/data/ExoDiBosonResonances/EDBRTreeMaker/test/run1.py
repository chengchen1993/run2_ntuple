import os
lines=[]
file=open("test.txt")
for line in file:
                lines.append(line);
for line in lines:
                line=line.split(" ");
                if "/store/data/Run2016" in line[len(line)-1]:
                   print (line[len(line)-1])[0:108]
