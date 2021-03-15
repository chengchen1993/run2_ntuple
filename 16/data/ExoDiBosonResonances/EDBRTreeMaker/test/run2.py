import os
lines=[]
file=open("list.txt")
for line in file:
                lines.append(line);
i=0
for line in lines:
                i=i+1
                os.system("cp H_test.py "+str(i)+"H.py") 
                print "sed -i 's#ZZZZZZZZ#"+line[:-1]+"#g' "+str(i)+"H.py"
                print "sed -i 's/pickup/pickup"+str(i)+"H/g' "+str(i)+"H.py"  
