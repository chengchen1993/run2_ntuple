f=open("tt_pweight.txt","r")
num=[]
for line in f:
    if "weight id" in line:
         line= line.split('"')
         num.append(int(line[1]))
print(len(num))
print(num)
