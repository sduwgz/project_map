#!/usr/bin/python

def readRightLocIdx():
    f = open("right_loc_idx");
    line = f.readline().strip()
    data = line.split()      
    moleStart = {}
    moleEnd = {}
    while line:
        moleIdx = int(data[0])
        score = float(data[1])
        if score > 20:
            moleStart[moleIdx] = int(data[6])
            moleEnd[moleIdx] = int(data[7])
        line = f.readline().strip()
        data = line.split()  
    f.close()
    return moleStart,moleEnd

def readMqrMap():
    moleStart,moleEnd = readRightLocIdx()
    out = open("right.map","w")
    f = open("MQR.map")
    for i in range(9):
        f.readline();
    line = f.readline().strip()
    data = line.split()      
    while line:
        moleIdx = int(data[1])
        if moleIdx in moleStart:
            diff = moleStart[moleIdx] - int(data[19])
            if moleEnd[moleIdx] - int(data[len(data)-1]) == diff+1 and len(data)>48:
                out.write(str(moleIdx)+"\t")
                for i in range((len(data)-18)/3):
                    #print data[i*3+18]
                    out.write(str(data[i*3+18])+"\t")
                    if data[i*3+19] != '-1' and data[i*3+20] != '-1':
                        #print int(data[i*3+19]) + diff
                        #print int(data[i*3+20]) + diff
                        out.write(str(int(data[i*3+19])+diff)+"\t"+str(int(data[i*3+20])+diff)+"\t")
                    else:
                        #print "-1\t-1"
                        out.write("-1\t-1\t")
                out.write("\n")    
        line = f.readline().strip()
        data = line.split()    
    f.close()


def main():
    readMqrMap()

main()        
