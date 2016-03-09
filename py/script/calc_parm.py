import sys
import string
from itertools import ifilter, imap
import numpy

def read_result(f):
    settled_mole = []
    for line in f.readlines():
        mole_id = int(line.split()[0])
        settled_mole.append(mole_id)
    return settled_mole

if __name__ == '__main__':
    with file(sys.argv[1]) as f:
        settled_mole = read_result(f)
    state = 0
    insert = []
    insert = []
    delta = []
    count = 0
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, sys.stdin)):
        if len(line.split()) == 7:
            mole_id = int(line.split()[0])
            if mole_id in settled_mole:
                state = 1
            else:
                state = 0
        else:
            count += 1
            if state == 1:
                insert.append(int(line.split()[1]) - int(line.split()[0]))
                insert.append(int((int(line.split()[3]) - int(line.split()[2]) + 0.0) / int(line.split()[4]) * 10000 + 0.5))
                delta.append(int(line.split()[4]) - int(line.split()[5]))
                if int(line.split()[3]) - int(line.split()[2]) > 0:
                    #print '%d %d'%(int(line.split()[3]) - int(line.split()[2]), int(line.split()[4]))
                    print '%f'%((int(line.split()[3]) - int(line.split()[2]) + 0.0) / int(line.split()[4]) * 10000)
                else:
                    print 0
            else:
                continue
    print count
    insert_array = numpy.array(insert)
    stat = [0] * 6
    for i in insert:
        stat[i] += 1
    print stat
    print insert_array.mean()
    print insert_array.var()
'''
    delta_array = numpy.array(delta)
    
    print delta_array.mean()
    print delta_array.var() ** 0.5
    stat = [0] * 7
    for i in insert:
        stat[i] += 1
    print stat
    
    insert_array = numpy.array(delete)

    stat = [0] * 20
    for i in insert:
        stat[i] += 1
    print stat
    print insert_array.mean()
    print insert_array.var()
'''        
