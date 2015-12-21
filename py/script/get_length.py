#!/usr/bin/python

def readInputBnx():
    f = open("input.bnx");
    for i in range(19):
        f.readline();
    moleLenExp = {}
    molecular = {}
    line = f.readline().strip()
    data = line.split()      
    while line:
        moleIdx = int(data[1])
        moleLenExp[moleIdx] = float(data[2])
        line = f.readline().strip()
        data = line.split()  
        assert int(data[0]) == 1
        molecular[moleIdx] = []
        for i in range(len(data)-2):
            d = float(data[i+1])
            molecular[moleIdx].append(d)
        f.readline()
        f.readline()
        line = f.readline().strip()
        data = line.split()    
    f.close()
    # for idx in moleLenExp:
       # print idx,moleLenExp[idx]
    return moleLenExp,molecular

def readMqrMap():
    moleLenExp,molecular = readInputBnx()
    f = open("MQR.map")
    for i in range(9):
        f.readline();
    moleLenAct = {}
    moleLenActM = {}
    moleLenExpM = {}
    line = f.readline().strip()
    data = line.split()      
    while line:
        moleIdx = int(data[1])
        moleLenAct[moleIdx] = float(data[8]) - float(data[7])
        moleLenActM[moleIdx] = float(data[10]) - float(data[9])
        startInMole = int(data[18])
        endInMole = int(data[-3])
        moleLenExpM[moleIdx] = abs(molecular[moleIdx][startInMole-1] - molecular[moleIdx][endInMole-1]) 
        line = f.readline().strip()
        data = line.split()    
    f.close()

    #for idx in moleLenAct:
       # print idx,moleLenAct[idx]
    return moleLenExp,moleLenAct,moleLenExpM,moleLenActM 

def main():
    out = open("mole_length","w")
    moleLenExp,moleLenAct,moleLenExpM,moleLenActM = readMqrMap()
    error = {}
    errorM = {}
    for idx in moleLenAct:
        assert (idx in moleLenExp) and (idx in moleLenAct) and (idx in moleLenExpM) and (idx in moleLenActM)
        error[idx] = moleLenAct[idx] - moleLenExp[idx]
        errorM[idx] = moleLenActM[idx] - moleLenExpM[idx]
        out.write(str(idx) + "\t" 
                + str(moleLenAct[idx]) + "\t" + str(moleLenExp[idx]) + "\t" + str(error[idx]) + "\t" 
                + str(moleLenActM[idx]) + "\t" + str(moleLenExpM[idx]) + "\t" + str(errorM[idx]) +"\n")

main()        
