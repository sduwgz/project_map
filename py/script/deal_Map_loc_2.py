#!/usr/bin/python
import os

def read_First_Map():
    f = open("data/mole.map.miss")
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

def getFileList(filePath):
    filePath = str(filePath)
    if filePath == "":
        return []
    #filePath = filePath.replace("/","\\")
    #if filePath[-1] != "\\":
     #   filePath = filePath + "\\"
    a = os.listdir(filePath)
    b = [x for x in a if os.path.isfile(filePath+x)]
    return b


def read_mean(filename):
    f = open("new_test_map/"+filename)
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
    fileList = getFileList("./new_test_map/") 
    for filename in fileList:
        moleMapone,moleMaponeEnd = read_mean(filename)
        f = open("data/Map_loc_idx")
        print "begin cheak"
        line = f.readline().strip()
        data = line.split()      
        delmole = 0;
        delmoleshort = 0;
        delwrong = 0;
        out = open("new_diss/" + filename+"_diss","w")
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
            if moleStart <= moleMapone[moleIdx] and moleEnd >= moleMaponeEnd[moleIdx]:
                line = f.readline().strip()
                data = line.split() 
                continue
            if moleStart == moleMapone[moleIdx] or moleEnd == moleMaponeEnd[moleIdx]:
                line = f.readline().strip()
                data = line.split() 
                continue
            out.write("map diff"+str(moleIdx)+"\n")
            print "map diff"+str(moleIdx)           
            line = f.readline().strip()
            data = line.split() 
        out.write("del mole num:"+str(delmole) +"\n")
        out.write("del mole short num:"+str(delmoleshort)+"\n")
        out.write("del mole wrong num:"+str(delwrong))




main()
