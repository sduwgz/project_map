#!/usr/bin/env python
#-*- coding : utf-8 -*-

__author__ = "Li Yanbo"


#according to align start location and end location get start Enzyme id and end Enzyme id

def readRef():
    f = open("ref.cmap")
    for i in range(6):
        f.readline()
   
    enzyme = {}
  
    while True:
        line = f.readline().strip()
        if not line: break;
        data = line.split()
        enzymeIdx = int(data[3])
        enzymeLoc = int(float(data[5]))
        enzyme[enzymeIdx] = enzymeLoc
    f.close()
    with open("Enzyme_Idx_Loc","w") as out:
        for idx in enzyme: out.write(str(idx) + "\t" + str(enzyme[idx]) + "\n")
    return enzyme 

if __name__ == "__main__":
    enzyme = readRef()
    f = open("map_result","r") #align start location and end location
    f.readline()

    out = open("Map_Idx_In_Ref","w")
    while True:
        line = f.readline().strip()
        if not line: break;
        data = line.split()
        moleIndex = int(data[0])
        out.write(str(moleIndex) + "\t")
        moleStartMatchLoc = int(data[4])
        moleEndMatchLoc = int(data[5])
        for idx in enzyme:
            if enzyme[idx] == moleStartMatchLoc:
                out.write(str(idx) + "\t")
            if enzyme[idx] == moleEndMatchLoc:
                out.write(str(idx))
                break
        out.write("\n")     
    f.close()
    out.close()



