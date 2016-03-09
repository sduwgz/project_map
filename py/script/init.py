import sys
import string
from itertools import ifilter, imap

def read_ref(f):
    pre_pos = 0
    ref_list = []
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        if line.startswith('#'):
            continue
        else:
            pos = int(float(line.split()[5]))
            ref_list.append(pos - pre_pos)
            pre_pos = pos
    return ref_list

def read_moles(f):
    mole_set = {}
    state = 0
    mole_id = 0
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        if line.startswith('#'):
            continue
        else:
            if state == 0:
                assert(line.startswith('0'))
                mole_id = int(line.split()[1])
                mole_set[mole_id] = []
                state = 1
                continue
            if state == 1:
                #assert(line.startswith('1'))
                if len(line.split()) < 5:
                    state = 2
                    continue
                pre_pos = int(float(line.split()[1]))
                for p in line.split()[2 : -1]:
                    pos = int(float(p))
                    mole_set[mole_id].append(pos - pre_pos)
                    pre_pos = pos
                    state = 2
                continue
            if state == 2:
                state = 3
                continue
            if state == 3:
                state = 0
                continue
                
    return mole_set

def read_map(f):
    return 1
'''
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        mole_set = {}
        ref = []
        mole_id = int(line.split()[1])
        direct = line.split()[6]
        pos = 18
        while pos < len(line):
            left = pos + 1
            right = pos + 2
            if left == -1 or right ==  -1:
                pos += 3
                continue
            start = line.split()[pos]
            pos += 3

'''     
