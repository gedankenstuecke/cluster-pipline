import sys

infile = open(sys.argv[1],"r")
data = infile.readlines()
outdata = []

for i in data:
    a = i.split("\t")
    outdata.append(a)

for i in outdata:
    print ">"+i[0]
    print i[5][:-1]
