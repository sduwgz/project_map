#!/usr/bin/python

def read_new_ref():
    f = open("new_re.fa")
    f.readline()
    line = f.readline()
    n = line.count("GCTCTTC")
    i = 0
    pos = -1
    case = 1
    while i < n:
        pos = line.index("GCTCTTC",pos+1)
        print pos,case
        i += 1
    n = line.count("GAAGAGC")
    i = 0
    pos = -1
    case = -1
    while i < n:
        pos = line.index("GAAGAGC",pos+1)
        print pos,case
        i += 1
   
    n = line.count("GCTCTTG")
    i = 0
    pos = -1
    case = 2
    while i < n:
        pos = line.index("GCTCTTG",pos+1)
        print pos,case
        i += 1
    n = line.count("CAAGAGC")
    i = 0
    pos = -1
    case = -2
    while i < n:
        pos = line.index("CAAGAGC",pos+1)
        print pos,case
        i += 1


def main():
    read_new_ref()
main()    
