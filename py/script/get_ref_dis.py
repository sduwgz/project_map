#!/usr/bin/python

def main():
    f = open("ref.cmap")
    for i in range(6):
        f.readline();
   
    pos = {}
    dis = {}
    line = f.readline().strip()
    data = line.split()
    while int(data[4]) == 1:
        pos[int(data[3])] = float(data[5])
        line = f.readline().strip()
        data = line.split()
    f.close()
    for idx in pos:
        if idx+1 not in pos:
            break
        dis[idx] = pos[idx+1] - pos[idx]
        print idx+1,idx,dis[idx]
    #for idx in dis:
     #   print "dis[%d]"%(idx),dis[idx]

main()        
