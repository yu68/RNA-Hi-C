import sys

pair=open(sys.argv[1],"r")
output=open(sys.argv[2],"w")

length = 10
if len(sys.argv) > 3:
    length = int(sys.argv[3])

i=1
for line in pair:
  if i%2==0 and i%4!=0:
    line=line.strip()[length:]
  print >>output,line.strip()
  i+=1
