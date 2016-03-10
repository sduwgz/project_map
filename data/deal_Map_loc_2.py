#!/usr/bin/python

def read_First_Map():
    f = open("mole.map.miss")
    moleMap = {}
    line = f.readline().strip()
    data = line.split()      
    while line:
        assert len(data) == 2
        moleIdx = int(data[0])
        line = f.readline().strip()
        data = line.split()   
        moleMap[moleIdx] = []
        for i in range(len(data)):
            moleMap[moleIdx].append(int(data[i])) 
        line = f.readline().strip()
        data = line.split()   
    #for idx in moleMap:
     #   print str(idx) + " " + str(moleMap[idx])
    return moleMap

def read_mean():
    f = open("mean_68_iter4")
    moleMap = {}
    moleMapEnd = {}
    line = f.readline().strip()
    data = line.split()      
    while line:
        moleIdx = int(data[0])
        moleMap[moleIdx] = int(data[2])
        moleMapEnd[moleIdx] = int(data[3])
        line = f.readline().strip()
        data = line.split()   
    return moleMap,moleMapEnd

def main():
    moleMapmore = read_First_Map()
    moleMapone,moleMaponeEnd = read_mean()
    f = open("Map_loc_idx")
    print "begin cheak"
    line = f.readline().strip()
    data = line.split()      
    delmole = 0;
    delmoleshort = 0;
    delwrong = 0;
    while line:
        if len(data) < 8: 
            line = f.readline().strip()
            data = line.split() 
            delmole = delmole+1;
            continue
        if int(data[7]) - int(data[6]) < 10: 
            line = f.readline().strip()
            data = line.split() 
            delmoleshort = delmoleshort + 1;
            continue
        moleIdx = int(data[0])
        if data[3] == '-':
            moleIdx = -moleIdx
        moleStart = int(data[6])
        moleEnd = int(data[7])
        if moleIdx not in moleMapmore:
            line = f.readline().strip()
            data = line.split() 
            delwrong = delwrong+1;
            continue
        if moleStart not in moleMapmore[moleIdx]:
            line = f.readline().strip()
            data = line.split() 
            delwrong = delwrong+1;
            continue
        if moleStart != moleMapone[moleIdx] and moleEnd != moleMaponeEnd[moleIdx]:
            print "map diff"+str(moleIdx)           
        line = f.readline().strip()
        data = line.split() 
    print "del mole num:"+str(delmole) 
    print "del mole short num:"+str(delmoleshort) 
    print "del mole wrong num:"+str(delwrong) 




main()
